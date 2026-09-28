"""Codon-optimize proteins for synthesis, and gate each fragment on manufacturability.

Two halves:

1. **Optimise locally** with DnaChisel (MIT): maximise CAI, minimise repeated
   k-mers and hairpins, and forbid the patterns that break either synthesis or
   expression (homopolymers, polymerase pausing sites, and — for a bacterial
   host — internal Shine-Dalgarno / strong RBS / Chi site / cryptic starts).
2. **Gate on manufacturability**: score each candidate locally (cheap
   heuristics), then ask the vendor for the real complexity score and, if it is
   at or above the threshold, re-optimise with a higher repeat penalty.  The
   ``kmers_weight`` ramp is what actually buys synthesisability.

The recipe follows the Baker lab's SAPP/DMX release (``JB/domesticator.py`` +
``JB/idt.py`` of github.com/bwicky/SAPP_DMX, MIT), published with Qian, Milles,
Wicky, Ragotte et al., *Nat Commun* 2026.  Deviations from upstream and the
reasoning behind every setting are recorded in ``docs/design.md``.

This module is the algorithmic core and holds no I/O: see ``cli.py``.
"""

from __future__ import annotations

import argparse
import copy
import json
import platform
import random
import sys
import time
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


# Heavy import kept at module level on purpose: this tool is useless without it,
# and `import codon_optimize` is never done from the torch-only code paths.
import dnachisel as dc
from dnachisel import DnaOptimizationProblem, Location, NoSolutionError

# --- patterns that hurt synthesis or expression -----------------------------
# Shared by every species; each entry is (pattern, why).
GENERIC_AVOID = (
    ("AAAAA", "poly-A (terminator-like, synthesis)"),
    ("TTTTT", "poly-T (terminator-like, synthesis)"),
    ("CCCCCC", "poly-C (synthesis)"),
    ("GGGGGG", "poly-G (synthesis)"),
    ("ATCTGTT", "T7/T3 RNA polymerase pausing site"),
    ("GGRGGT", "G-quadruplex-like"),
)

# Any bacterial host: keep the mRNA from initiating anywhere but the 5' end.
# Meaningless in a eukaryote, which has no Shine-Dalgarno and initiates by
# cap-dependent scanning.
BACTERIAL_AVOID = (
    ("GGAGG", "internal Shine-Dalgarno"),
    ("TAAGGAG", "strong RBS"),
)

# E. coli only: the Chi site is recognised by RecBCD, which other genera lack.
ECOLI_AVOID = (
    ("GCTGGTGG", "Chi site (RecBCD)"),
)

# Which host class a codon table belongs to.  Unknown keys are treated as
# eukaryotic, i.e. only GENERIC_AVOID applies — the conservative choice, since
# applying a bacterial initiation rule to a eukaryote is simply wrong, whereas
# omitting one from an unlisted bacterium only loses an optimisation.
BACTERIAL_SPECIES = frozenset({"e_coli", "b_subtilis"})


def is_bacterial(species: str) -> bool:
    """True when internal-SD / RBS avoidance applies to this host."""
    return any(species.startswith(s) for s in BACTERIAL_SPECIES)


def is_ecoli(species: str) -> bool:
    """True when the E. coli-only patterns (Chi, cryptic starts) apply."""
    return species.startswith("e_coli")

# Cryptic start sites: a G/A-rich stretch 5-7 nt upstream of ATG/GTG/TTG.
# Dropped first when a sequence has no solution under the full constraint set.
ECOLI_ALT_START = (
    "RRRRRNNNNNDTG",
    "RRRRRNNNNNNDTG",
    "RRRRRNNNNNNNDTG",
)

# No restriction-site avoidance by default: the batch is not cloned by Golden
# Gate.  Pass --avoid with both strands of a site if that ever changes
# (BsaI would be "GGTCTC GAGACC").
DEFAULT_ENZYME_SITES: tuple[str, ...] = ()

# GC bands, as used for TWIST/IDT-style dsDNA synthesis.
GC_GLOBAL = (0.25, 0.65)
GC_WINDOW = (0.35, 0.65, 50)

# IDT's own published trouble thresholds for dsDNA fragments, used as the local
# pre-gate so we spend API calls only on plausible candidates.
IDT_GC_LIMITS = (0.25, 0.75)
IDT_MAX_AT_RUN = 10
IDT_MAX_GC_RUN = 6

STOP_CODONS = ("TAA", "TGA", "TAG")

# The window that matters, in nt relative to the A of ATG.  Not a guess: Kudla
# et al. (Science 2009) ran a moving-window analysis over 154 synonymous GFP
# variants and found -4..+37 — the ~30-nt ribosome binding site centred on the
# start codon — whose folding energy explained 44% of the variance in protein
# level (r = 0.66).  Notably the best window did *not* overlap the
# Shine-Dalgarno sequence, so this is a start-codon effect, not an SD one.
FOOTPRINT = (-4, 37)
# Codons of the ORF covered by that window (+37 nt ~ 12 codons).  Their AT
# content is a vector-independent proxy for the same thing: AT-rich pairs weakly,
# so it cannot form the structure the window is scored for.  This is the
# quantitative form of the cell-free rule of thumb "make the first ~10 codons
# AT-rich", and unlike dg_open_5p it needs no knowledge of the UTR.
FIVE_PRIME_CODONS = 12
# How much of the ORF to fold: the footprint plus room for competing structures.
# The 5' side is NOT a window — it is the real transcript start (see --utr5).
FOLD_ORF_NT = 60
# Codons that stage 1 decides and stage 2 must not touch: the whole folded window,
# not just the footprint.  Measured: re-choosing only the 24 nt *downstream* of the
# footprint moved a fixed head's dg_open_5p by 9-10 kcal/mol — as much as the head
# itself — so freezing 12 codons would let stage 2 undo stage 1's work.
HEAD_CODONS = FOLD_ORF_NT // 3

# IUPAC ambiguity codes, for checking the ambiguous patterns (GGRGGT, the
# cryptic-start motifs) against a candidate head in stage 1.
IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "[AG]", "Y": "[CT]",
         "S": "[GC]", "W": "[AT]", "K": "[GT]", "M": "[AC]", "B": "[CGT]",
         "D": "[AGT]", "H": "[ACT]", "V": "[ACG]", "N": "[ACGT]"}


class MinimizeNumKmers(dc.Specification):
    """Penalise repeated k-mers — the knob that buys synthesisability.

    Ported from ``JB/domesticator.py`` of the SAPP/DMX release (MIT).  The score
    is negative and proportional to the fraction of the sequence sitting inside
    non-unique k-mers, so raising ``boost`` trades CAI for fewer repeats.
    """

    best_possible_score = 0

    def __init__(self, k: int = 8, location=None, boost: float = 1.0) -> None:
        self.location = location
        self.k = k
        self.boost = boost

    def initialize_on_problem(self, problem, role=None):
        return self._copy_with_full_span_if_no_location(problem)

    def evaluate(self, problem):
        sequence = self.location.extract_sequence(problem.sequence)
        kmers = [sequence[i:i + self.k] for i in range(len(sequence) - self.k)]
        n_non_unique = sum(c for _, c in Counter(kmers).items() if c > 1)
        score = -(float(self.k) * n_non_unique) / len(sequence)
        return dc.SpecEvaluation(
            self, problem, score=score, locations=[self.location],
            message="Score: %.02f (%d non-unique %d-mers)"
                    % (score, n_non_unique, self.k),
        )

    def label_parameters(self):
        return [("k", str(self.k))]

    def short_label(self):
        return f"Avoid {self.k}mers {self.boost}"

    def __str__(self):
        return "MinimizeNum%dmers" % self.k


# ---------------------------------------------------------------------------
# optimisation
# ---------------------------------------------------------------------------
def build_constraints(
    orf_location: Location,
    full_location: Location,
    stop_location: Location | None,
    species: str,
    enzyme_sites: tuple[str, ...],
    with_alt_start: bool,
    max_gc_start: float | None = None,
) -> list:
    """Constraint list for one ORF.

    ``orf_location`` is the translated span, ``full_location`` additionally
    covers the stop codon: patterns and GC are enforced across the whole thing
    so nothing bad is created at the junction.

    NB upstream's ``constraints_easier`` is deep-copied *before* the GC
    constraints are appended, so its relaxed retry silently drops the GC bands
    (we saw 50-bp windows reach 72% GC that way).  Here the relaxed set differs
    from the strict one only by the alternative-start patterns.
    """
    cons = [dc.EnforceTranslation(location=orf_location)]
    if stop_location is not None:
        cons.append(dc.AvoidChanges(location=stop_location))

    for pattern, _why in GENERIC_AVOID:
        cons.append(dc.AvoidPattern(pattern, location=full_location))
    for site in enzyme_sites:
        cons.append(dc.AvoidPattern(site, location=full_location))

    if is_bacterial(species):
        for pattern, _why in BACTERIAL_AVOID:
            cons.append(dc.AvoidPattern(pattern, location=full_location))
    if is_ecoli(species):
        for pattern, _why in ECOLI_AVOID:
            cons.append(dc.AvoidPattern(pattern, location=full_location))
        if with_alt_start:
            for pattern in ECOLI_ALT_START:
                cons.append(dc.AvoidPattern(pattern, location=full_location))

    lo, hi = GC_GLOBAL
    cons.append(dc.EnforceGCContent(mini=lo, maxi=hi, location=full_location))
    wlo, whi, window = GC_WINDOW
    cons.append(dc.EnforceGCContent(mini=wlo, maxi=whi, window=window,
                                    location=full_location))

    # Hold the first FIVE_PRIME_CODONS codons AT-rich.  This is the *intervention*
    # form of the cell-free rule of thumb, and it does what merely ranking
    # candidates cannot: measured on a 270-aa bacterial ORF, capping head GC at 0.40 raised
    # AT content 0.42-0.47 -> 0.61 and dropped the best dg_open_5p from 11.1 to
    # 8.3 kcal/mol for ~0.011 of CAI.  Composition, so it needs no vector context.
    if max_gc_start is not None:
        head = min(FIVE_PRIME_CODONS * 3, orf_location.end)
        cons.append(dc.EnforceGCContent(maxi=max_gc_start, location=Location(0, head)))
    return cons


def reverse_translate(
    aa_sequence: str,
    species: str = "e_coli",
    kmers_weight: float = 10.0,
    enzyme_sites: tuple[str, ...] = DEFAULT_ENZYME_SITES,
    stop_codon: str | None = "TAA",
    max_tries: int = 20,
    max_gc_start: float | None = None,
    frozen_head: str | None = None,
) -> tuple[str, bool, bool]:
    """Optimise a DNA sequence coding for ``aa_sequence``.

    Returns ``(dna, alt_start_relaxed, at_rich_relaxed)``.  The ORF carries
    ``stop_codon`` when one is given; the stop itself is frozen, but patterns and
    GC are still enforced across the coding/stop junction.

    DnaChisel's search is stochastic, so a failed attempt is retried, and the
    attempts form a ladder that gives up constraints in order of least value:
    first the cryptic-start-site patterns (known to be unsatisfiable for some
    sequences), then the AT-rich head.  Whichever a solution needed is reported,
    so a construct that could not get its AT-rich start is visible rather than
    silently different from the rest of the batch.
    """
    bad = sorted(set(aa_sequence) - set("ACDEFGHIKLMNPQRSTVWY"))
    if bad:
        raise ValueError(f"non-standard residue(s) {bad} — cannot reverse translate")

    naive = dc.reverse_translate(aa_sequence)
    orf_len = len(aa_sequence) * 3
    if frozen_head:
        # Stage 2: keep stage 1's window verbatim, optimise only what follows.
        if dc.translate(frozen_head) != aa_sequence[:len(frozen_head) // 3]:
            raise ValueError("frozen_head does not translate to the ORF's start")
        naive = frozen_head + naive[len(frozen_head):]
    if stop_codon:
        naive = naive + stop_codon
    full_len = len(naive)

    orf_location = Location(0, orf_len)
    full_location = Location(0, full_len)
    stop_location = Location(orf_len, full_len) if stop_codon else None

    objectives = [
        MinimizeNumKmers(k=8, boost=kmers_weight, location=orf_location),
        dc.AvoidHairpins(boost=1.0, location=orf_location),
        dc.MaximizeCAI(species=species, boost=1.0, location=orf_location),
    ]

    # Relaxation ladder, most constrained first. Alt-start suppression is dropped
    # before the AT-rich head because the head buys more expression.
    rungs = [(True, max_gc_start)]
    if species == "e_coli":
        rungs.append((False, max_gc_start))
    if max_gc_start is not None:
        rungs.append((rungs[-1][0], None))
    per_rung = max(1, max_tries // len(rungs))

    solutions: list[tuple[str, bool, bool]] = []
    scores: list[float] = []
    last_error: NoSolutionError | None = None
    for attempt in range(max_tries):
        with_alt_start, gc_start = rungs[min(attempt // per_rung, len(rungs) - 1)]
        constraints = build_constraints(
            orf_location, full_location, stop_location, species, enzyme_sites,
            with_alt_start=with_alt_start, max_gc_start=gc_start,
        )
        if frozen_head:
            constraints.append(dc.AvoidChanges(location=Location(0, len(frozen_head))))
        problem = DnaOptimizationProblem(
            naive, constraints=constraints,
            objectives=copy.deepcopy(objectives), logger=None,
        )
        try:
            problem.resolve_constraints_by_random_mutations()
            problem.optimize()
            problem.resolve_constraints(final_check=True)
        except NoSolutionError as exc:
            last_error = exc
            continue
        solutions.append((problem.sequence, not with_alt_start,
                          max_gc_start is not None and gc_start is None))
        scores.append(problem.objectives_evaluations().scores_sum())
        break

    if not solutions:
        raise NoSolutionError(
            f"no solution after {max_tries} attempts", last_error.problem
            if last_error is not None else None,
        )
    # Upstream takes argmin here; objective scores are <= 0 and larger is
    # better, so take the best (max) when several attempts succeeded.
    best = int(np.argmax(scores))
    return solutions[best]


# ---------------------------------------------------------------------------
# local metrics (vendor-free pre-gate)
# ---------------------------------------------------------------------------
def gc_fraction(seq: str) -> float:
    return (seq.count("G") + seq.count("C")) / len(seq)


def window_gc_range(seq: str, window: int = 50) -> tuple[float, float]:
    if len(seq) <= window:
        g = gc_fraction(seq)
        return g, g
    values = [gc_fraction(seq[i:i + window]) for i in range(len(seq) - window + 1)]
    return min(values), max(values)


def longest_homopolymer(seq: str, bases: str = "ACGT") -> int:
    """Longest run of a *single* repeated base, restricted to ``bases``.

    IDT's published limits are homopolymeric ("10 or more As and Ts", "6 or more
    Gs and Cs"), i.e. AAAAAAAAAA — not a mixed A/T or G/C stretch, which is what
    the GC-window constraint already covers.
    """
    best = 0
    run = 0
    prev = ""
    for cur in seq:
        run = run + 1 if cur == prev else 1
        prev = cur
        if cur in bases:
            best = max(best, run)
    return best


def non_unique_kmer_fraction(seq: str, k: int = 8) -> float:
    kmers = [seq[i:i + k] for i in range(len(seq) - k + 1)]
    if not kmers:
        return 0.0
    return sum(c for c in Counter(kmers).values() if c > 1) / len(kmers)


def local_metrics(seq: str) -> dict:
    lo, hi = window_gc_range(seq)
    return {
        "length_bp": len(seq),
        "gc": round(gc_fraction(seq), 4),
        "gc_win50_min": round(lo, 4),
        "gc_win50_max": round(hi, 4),
        "homopolymer_max": longest_homopolymer(seq),
        "at_run_max": longest_homopolymer(seq, "AT"),
        "gc_run_max": longest_homopolymer(seq, "GC"),
        "nonunique_8mer_frac": round(non_unique_kmer_fraction(seq), 4),
    }


def five_prime_accessibility(
    orf: str,
    utr5: str = "",
    fold_orf_nt: int = FOLD_ORF_NT,
) -> dict:
    """How hard it is to keep the ribosome footprint single-stranded.

    In *E. coli* protein output is usually limited by translation **initiation**,
    not elongation: the 30S subunit can only load where the RBS and start codon
    are unpaired, so a hairpin over that window suppresses expression however
    good the downstream codons are.  Kudla et al. (*Science* 2009) made this
    concrete — 154 synonymous GFP variants spanned 250-fold in protein level,
    codon bias did not correlate, and mRNA folding stability near the ribosome
    binding site explained more than half the variance.

    ``dg_open_5p`` is the cost (kcal/mol, >= 0) of forcing the footprint to stay
    unpaired: MFE with that region constrained open, minus the unconstrained MFE.
    0 means the start is already accessible; large means it is sequestered.
    ``mfe_5p`` is the window's own MFE, for context.

    The scored window is Kudla's ``FOOTPRINT`` (-4..+37), which is almost entirely
    *inside* the ORF — only 4 nt of it are upstream of the ATG.  So the ORF's own
    first ~12 codons dominate the number, and ``at_frac_5p`` (their AT content)
    captures the same effect without needing any vector sequence at all.

    ``utr5``, when given, must be **the real 5' untranslated region of the
    transcript**: the sequence from the transcription start site up to (not
    including) the ATG, as it will exist in the assembled plasmid.  It is used
    whole, never truncated, because the molecule that folds in the cell begins at
    the TSS — shortening it deletes base-pairing partners and can make a
    sequestered start look open, or the reverse, unpredictably.  It refines the
    number rather than rescuing it: how much the UTR adds is itself contested
    (Kudla's best window avoids the Shine-Dalgarno entirely, while *Mol Cell*
    2018 argues SD occlusion by the coding region is the mechanism), and a T7-style
    system that decouples transcription from translation reportedly damps the
    synonymous effect on mRNA level.

    It is emphatically **not** the iVEC homology arm.  The arm is the DNA that
    overlaps the vector; it need not begin at the TSS, and since it is a copy of
    vector sequence, concatenating arm + ORF can count part of the UTR twice.
    The fragment we order is not the transcript: after recombination the arm
    merges into the vector, and what gets transcribed is vector-UTR + ORF.

    With no ``utr5`` the ORF is folded alone.  Since the decisive window is mostly
    intra-ORF that is a reasonable approximation, not a broken one — it just
    cannot see structure formed across the UTR junction.  ``utr5_bp`` records
    which case produced the value.
    """
    # Vector-independent, and computable without ViennaRNA: the AT content of the
    # codons inside the window. See FIVE_PRIME_CODONS.
    head = orf[:FIVE_PRIME_CODONS * 3]
    at_frac_5p = round(1.0 - gc_fraction(head), 4) if head else None

    try:
        import RNA  # ViennaRNA; optional
    except ImportError:
        return {"dg_open_5p": None, "mfe_5p": None, "utr5_bp": None,
                "at_frac_5p": at_frac_5p}

    utr = utr5.strip().upper()
    window = (utr + orf[:fold_orf_nt]).upper().replace("T", "U")
    start = len(utr)                                   # index of the A of AUG
    lo = max(0, start + FOOTPRINT[0])
    hi = min(len(window), start + FOOTPRINT[1])

    free = RNA.fold_compound(window)
    _, mfe_free = free.mfe()
    opened = RNA.fold_compound(window)
    for i in range(lo + 1, hi + 1):                    # ViennaRNA is 1-indexed
        opened.hc_add_up(i)
    _, mfe_open = opened.mfe()
    return {
        "dg_open_5p": round(mfe_open - mfe_free, 2),
        "mfe_5p": round(mfe_free, 2),
        "utr5_bp": len(utr),
        "at_frac_5p": at_frac_5p,
    }


def _forbidden_regexes(species: str, enzyme_sites: tuple[str, ...]) -> list:
    """Every sequence pattern the ORF must not contain, as compiled regexes."""
    import re

    patterns = [p for p, _ in GENERIC_AVOID] + list(enzyme_sites)
    if is_bacterial(species):
        patterns += [p for p, _ in BACTERIAL_AVOID]
    if is_ecoli(species):
        patterns += [p for p, _ in ECOLI_AVOID] + list(ECOLI_ALT_START)
    return [re.compile("".join(IUPAC[b] for b in p)) for p in patterns]


def optimize_head(
    aa_head: str,
    utr5: str = "",
    species: str = "e_coli",
    enzyme_sites: tuple[str, ...] = (),
    n_samples: int = 8000,
    rng: random.Random | None = None,
    n_best: int = 25,
) -> list[dict]:
    """Stage 1 — choose the folded window's codons for maximum start accessibility.

    Searches synonymous codon choices for the first ``HEAD_CODONS`` codons and
    returns the one with the lowest ``dg_open_5p``.  This optimises the quantity
    that matters *directly*, instead of constraining a proxy for it: composition
    turned out not to track the fold at all (along the measured Pareto front, an
    AT fraction of 0.44 appeared at both 6.5 and 10.3 kcal/mol, and 0.56 at both
    1.2 and 8.1), so a GC cap is a blunt instrument for this job.

    It is also far cheaper than sampling heads by re-optimising whole ORFs, which
    is what drawing ``--n-candidates`` full solutions amounts to: one head costs a
    ~2 ms fold instead of a ~1.4 s DnaChisel run, so thousands are affordable
    where tens were not.  On a 270-aa bacterial ORF this reached 0.0 kcal/mol where 30
    full-ORF draws had bottomed out at 3.8.

    Returns the best ``n_best`` heads, most accessible first, **not** a single
    winner: a head optimal in isolation can paint stage 2 into a corner. Measured
    on an 88-aa ORF, the dG-1.0 head left stage 2 with no solution at 10 of 10
    seeds while the dG-1.2 runner-up succeeded at 10 of 10 — a difference of 0.2
    kcal/mol against a total loss of feasibility. The caller takes the first head
    stage 2 can actually build on.

    Candidates must satisfy the same pattern and GC-window constraints stage 2
    enforces, since stage 2 cannot fix a frozen head.
    """
    import python_codon_tables as pct

    rng = rng or random.Random(0)
    table = pct.get_codons_table(species)
    synonyms = {aa: sorted(codons) for aa, codons in table.items() if aa != "*"}
    missing = sorted(set(aa_head) - set(synonyms))
    if missing:
        raise ValueError(f"no codons for residue(s) {missing}")

    forbidden = _forbidden_regexes(species, enzyme_sites)
    wlo, whi, window = GC_WINDOW
    scored: list[dict] = []
    seen: set[str] = set()
    evaluated = 0

    for _ in range(n_samples):
        head = "".join(rng.choice(synonyms[aa]) for aa in aa_head)
        if head in seen:
            continue
        seen.add(head)
        if any(r.search(head) for r in forbidden):
            continue
        if not GC_GLOBAL[0] <= gc_fraction(head) <= GC_GLOBAL[1]:
            continue
        lo, hi = window_gc_range(head, window)
        if lo < wlo or hi > whi:
            continue
        evaluated += 1
        acc = five_prime_accessibility(head, utr5, fold_orf_nt=len(head))
        scored.append({"head": head, **acc})

    if not scored:
        raise NoSolutionError(f"no valid head in {n_samples} samples", None)
    scored.sort(key=lambda c: c["dg_open_5p"])
    for candidate in scored[:n_best]:
        candidate["evaluated"] = evaluated
    return scored[:n_best]


def gate_span(orf: str, pad: str) -> str:
    """The part of a fragment re-optimising can actually change.

    Homology arms are excluded entirely.  They are fixed by the vector and they
    legitimately carry features we forbid inside an ORF — a vector-supplied RBS
    *is* a Shine-Dalgarno, a His tag is ~67% GC — so judging the pre-gate on
    them, or even on a junction margin, produces flags no objective can clear
    and burns the whole ramp on them.  Anything real at the junction shows up in
    IDT's score, which is computed on the true fragment, and the ramp can still
    respond by changing the ORF end.
    """
    return orf + pad


def local_flags(metrics: dict, total_length: int, min_len: int, max_len: int) -> list[str]:
    """Reasons a vendor is likely to push back, from published dsDNA limits.

    ``metrics`` describe the gate span; length is judged on the whole fragment.
    """
    flags = []
    if total_length < min_len:
        flags.append(f"short<{min_len}")
    if total_length > max_len:
        flags.append(f"long>{max_len}")
    if not IDT_GC_LIMITS[0] <= metrics["gc"] <= IDT_GC_LIMITS[1]:
        flags.append("gc_out_of_range")
    if metrics["gc_win50_max"] > GC_WINDOW[1] or metrics["gc_win50_min"] < GC_WINDOW[0]:
        flags.append("gc_window")
    if metrics["at_run_max"] >= IDT_MAX_AT_RUN:
        flags.append("at_run")
    if metrics["gc_run_max"] >= IDT_MAX_GC_RUN:
        flags.append("gc_run")
    return flags


def codon_adaptation_index(dna: str, species: str = "e_coli") -> float:
    """CAI of a coding sequence against the host's codon usage table.

    Single-codon amino acids (Met, Trp) and stops carry no information and are
    excluded, as in Sharp & Li's definition.
    """
    import python_codon_tables as pct

    table = pct.get_codons_table(species)
    weights: dict[str, float] = {}
    for aa, codons in table.items():
        if aa == "*" or len(codons) < 2:
            continue
        best = max(codons.values())
        if best <= 0:
            continue
        for codon, freq in codons.items():
            weights[codon] = max(freq / best, 1e-4)

    logs = [np.log(weights[dna[i:i + 3]])
            for i in range(0, len(dna) - len(dna) % 3, 3)
            if dna[i:i + 3] in weights]
    return float(np.exp(np.mean(logs))) if logs else float("nan")


# ---------------------------------------------------------------------------
# eBlock assembly
# ---------------------------------------------------------------------------
def make_pad(
    n_bases: int,
    preceding: str,
    enzyme_sites: tuple[str, ...],
    rng: random.Random,
    max_tries: int = 2000,
) -> str:
    """Neutral filler appended after the stop codon to reach the vendor minimum.

    Generated, not hand-picked, so it is reproducible from ``--seed``: balanced
    GC, no run longer than 3, no ``ATG`` (nothing to initiate on), and none of
    the forbidden patterns — checked across the junction with ``preceding``.
    """
    if n_bases <= 0:
        return ""
    forbidden = [p for p, _ in GENERIC_AVOID] + list(enzyme_sites) + ["ATG"]
    # The pad is non-coding filler. Keep every host's initiation motifs out of
    # it whatever the host: there is no cost, and a pad that can initiate is a
    # hazard in any expression system.
    forbidden += [p for p, _ in BACTERIAL_AVOID] + [p for p, _ in ECOLI_AVOID]
    tail = preceding[-12:]
    for _ in range(max_tries):
        pad = "".join(rng.choice("ACGT") for _ in range(n_bases))
        if not 0.40 <= gc_fraction(pad) <= 0.60:
            continue
        if longest_homopolymer(pad) > 3:
            continue
        # Start on a different base so the pad cannot extend a run the ORF ends
        # on (whose own length the ORF constraints already bound).
        if tail and pad[0] == tail[-1]:
            continue
        if any(site in tail + pad for site in forbidden):
            continue
        if non_unique_kmer_fraction(pad) > 0.0:
            continue
        return pad
    raise RuntimeError(f"could not generate a clean {n_bases} bp pad")


def assemble(orf: str, pad: str, adapter5: str, adapter3: str) -> str:
    """Full ordered fragment: homology arms stay terminal, pad sits inside them.

    The arms (iVEC) must match the vector exactly, so they are appended verbatim
    and never enter the optimisation.  The pad goes between the stop codon and
    the 3' arm — outside an arm it would leave a non-homologous flap, and on the
    5' side it would disturb the RBS-to-start spacing.
    """
    return adapter5 + orf + pad + adapter3


# ---------------------------------------------------------------------------
# per-gene driver
# ---------------------------------------------------------------------------
def optimize_gene(
    name: str,
    aa_sequence: str,
    args,
    client,
    rng: random.Random,
) -> dict:
    """Ramp the repeat penalty until the fragment clears the complexity gate.

    DnaChisel's search is stochastic, so each ramp step draws
    ``--n-candidates`` independent solutions and keeps the one whose start codon
    is *least* structured (lowest ``dg_open_5p``).  Codon choice barely moves
    expression through CAI, but it moves it a great deal through 5' mRNA
    structure, and this is the only place in the pipeline that can steer it.

    The pick is only as good as the folded molecule: without ``--utr5`` it ranks
    ORF 5' ends with no transcript context and is provisional.
    """
    weights = np.linspace(args.starting_kmers_weight, 100.0, args.n_steps)
    overhead = len(args.adapter5) + len(args.adapter3)
    best: dict | None = None

    # Stage 1, once per gene: pick the folded window, then freeze it. Doing this
    # here rather than inside the ramp is the point — the ramp only has to fix
    # synthesis complexity, and it can no longer perturb the start's accessibility.
    head = None
    head_rank = 0
    if args.head_samples > 0 and five_prime_accessibility("ATGAAA")["dg_open_5p"] is not None:
        n_head = min(HEAD_CODONS, len(aa_sequence))
        candidates = optimize_head(
            aa_sequence[:n_head], utr5=args.utr5, species=args.species,
            enzyme_sites=tuple(args.avoid), n_samples=args.head_samples, rng=rng,
        )
        # Take the most accessible head stage 2 can actually complete: the very
        # best can be infeasible downstream, at a dG cost of a fraction of a
        # kcal/mol to step past it.
        for rank, candidate in enumerate(candidates, start=1):
            np.random.seed(args.seed)
            try:
                reverse_translate(
                    aa_sequence, species=args.species, kmers_weight=float(weights[0]),
                    enzyme_sites=tuple(args.avoid), stop_codon=args.stop or None,
                    max_tries=args.max_attempts, frozen_head=candidate["head"],
                )
            except (NoSolutionError, ValueError):
                continue
            head, head_rank = candidate["head"], rank
            print(f"    stage 1: {candidate['evaluated']} valid heads -> "
                  f"dG_open {candidate['dg_open_5p']:.1f} kcal/mol, "
                  f"rank {rank} feasible ({n_head} codons frozen)")
            break
        else:
            print(f"    stage 1: none of {len(candidates)} heads left stage 2 a "
                  f"solution — falling back to single-stage optimisation")

    def rank(row: dict) -> tuple:
        """Lower is better: clean, then accessible 5', then complexity, repeats.

        ``dg_open_5p`` leads because it is the quantity Kudla measured; AT content
        of the same window breaks ties and carries the ranking on its own when
        ViennaRNA is absent (negated, since AT-rich is what we want).
        """
        score = row.get("idt_score")
        dg = row.get("dg_open_5p")
        at = row.get("at_frac_5p")
        return (
            len([f for f in row["local_flags"].split(";") if f]),
            float("inf") if score is None else score,
            float("inf") if dg is None else dg,
            -at if at is not None else 0.0,
            row["nonunique_8mer_frac"],
        )

    def remember(row: dict) -> None:
        nonlocal best
        if best is None or rank(row) < rank(best):
            best = row

    def build(step: int, weight: float, draw: int) -> dict | None:
        """One candidate: optimise, pad, assemble, measure."""
        # Distinct per (step, draw) so a rerun reproduces the same candidates.
        np.random.seed(args.seed + step * 1000 + draw)  # DnaChisel uses numpy
        try:
            orf, relaxed, at_relaxed = reverse_translate(
                aa_sequence, species=args.species, kmers_weight=float(weight),
                enzyme_sites=tuple(args.avoid), stop_codon=args.stop or None,
                max_tries=args.max_attempts,
                max_gc_start=(None if head is not None or args.max_gc_start >= 1.0
                              else args.max_gc_start),
                frozen_head=head,
            )
        except (NoSolutionError, ValueError) as exc:
            print(f"    step {step}.{draw}: no solution ({type(exc).__name__})")
            return None

        pad_len = max(0, args.min_length - (len(orf) + overhead))
        pad = make_pad(pad_len, orf, tuple(args.avoid), rng)
        fragment = assemble(orf, pad, args.adapter5, args.adapter3)

        # Reported metrics describe the whole ordered fragment (what the vendor
        # sees); flags are judged on the optimisable span only.
        metrics = local_metrics(fragment)
        flags = local_flags(local_metrics(gate_span(orf, pad)),
                            len(fragment), args.min_length, args.max_length)
        return {
            "name": name, "sequence": fragment, "orf_bp": len(orf),
            "pad_bp": len(pad), "cai": round(codon_adaptation_index(orf, args.species), 4),
            "kmers_weight": round(float(weight), 1), "steps_used": step,
            "alt_start_relaxed": relaxed, "at_rich_relaxed": at_relaxed,
            "head_codons": (len(head) // 3) if head else 0,
            "head_rank": head_rank,
            "local_flags": ";".join(flags),
            **five_prime_accessibility(orf, args.utr5, args.fold_orf_nt),
            **metrics,
        }

    for step, weight in enumerate(weights, start=1):
        drawn = [c for c in (build(step, weight, d) for d in range(args.n_candidates))
                 if c is not None]
        if not drawn:
            continue

        # Prefer a locally clean candidate; among equals, the most accessible 5'.
        clean = [c for c in drawn if not c["local_flags"]]
        ordered = sorted(clean or drawn, key=rank)
        if args.n_candidates > 1 and clean:
            dgs = [c["dg_open_5p"] for c in clean if c["dg_open_5p"] is not None]
            if dgs:
                print(f"    step {step}: {len(clean)}/{len(drawn)} clean, "
                      f"dG_open 5' {min(dgs):.1f}..{max(dgs):.1f} — "
                      f"taking {min(dgs):.1f} kcal/mol")
        candidate = ordered[0]
        flags = [f for f in candidate["local_flags"].split(";") if f]

        candidate["idt_score"] = None
        candidate["idt_pass"] = None

        # Cheap gate first: a locally flagged candidate is not worth an API call.
        if flags:
            print(f"    step {step}: local flags {flags} — re-optimising")
            remember(candidate)
            continue

        if client is None:
            return candidate

        # One API call per ramp step: only the top-ranked candidate is scored, and
        # a failure is answered by raising the repeat penalty rather than by
        # scoring this step's runners-up (which differ mainly in 5' structure).
        score = client.score(candidate["sequence"], name)
        candidate["idt_score"] = score
        candidate["idt_pass"] = bool(score < args.idt_score)
        if candidate["idt_pass"]:
            print(f"    step {step}: complexity {score:.1f} < {args.idt_score} — accepted")
            return candidate
        print(f"    step {step}: complexity {score:.1f} >= {args.idt_score} "
              f"— re-optimising")
        remember(candidate)

    if best is None:
        raise RuntimeError(f"{name}: no candidate produced in {args.n_steps} steps")
    # Ramp exhausted: keep the best seen, scored for real if it never was.
    if client is not None and best["idt_score"] is None:
        best["idt_score"] = client.score(best["sequence"], name)
        best["idt_pass"] = bool(best["idt_score"] < args.idt_score)
    print(f"    ramp exhausted — keeping best from step {best['steps_used']} "
          f"(complexity {best['idt_score']}, flags '{best['local_flags']}')")
    return best
