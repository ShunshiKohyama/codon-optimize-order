# Changelog

This tool designs sequences people pay to have synthesised, so a version is not
decoration: the same parameters can produce different DNA after an upgrade.
`<prefix>.meta.json` records the version of this package and of everything that
moves the output (DnaChisel, python_codon_tables, ViennaRNA, numpy). **Pin a tag
when a batch matters**, and keep the `meta.json` with the order:

```bash
pip install "git+https://github.com/ShunshiKohyama/codon-optimize-order@v0.1.1"
```

Versions follow [semantic versioning](https://semver.org/) with one addition:
**a change that alters designed sequences under unchanged parameters is a major
change**, even if no interface moves.

## v0.1.1 — 2026-10-13

**Fixed**

- Stage 1 (head selection) crashed with
  `NoSolutionError: no valid head in 8000 samples` for N termini whose every
  synonymous head matches an *E. coli* cryptic-start pattern. Seen with the mCherry
  family (`…GEEDNM…`: a Met after Glu-Glu-Asp-Asn). Stage 2 already dropped those
  patterns when they were unsatisfiable; stage 1 had no equivalent, and the
  exception was not caught.
  Heads that fail **only** the cryptic-start patterns are now held back and used
  only when no head passes the full set. The row is then marked
  `alt_start_relaxed`, and stage 2 starts below its strict rung
  (`reverse_translate(..., alt_start=False)`) instead of failing at it first. If even
  the relaxed filter leaves no head, the gene falls back to single-stage
  optimisation with a printed note instead of crashing. See `docs/design.md`.

**Added**

- `relax_alt_start` (`optimize_head`) and `alt_start` (`reverse_translate`);
  `relax_alt_start=False` restores the old all-or-nothing filter.
- `tests/test_stage1_relaxation.py` (`PYTHONPATH=src python -m unittest discover -s tests`).

**Why this is a patch, not a major change.** No designed sequence changes for an
input that has a head passing the full pattern set: the same samples are drawn and
ranked. Checked under the same seed against v0.1.0 — `examples/fragments.csv`
(2 rows, every output column identical) and a 14-protein panel (ubiquitin,
thioredoxin, the sfGFP example, ten random 120-mers, one mCherry-family protein):
13 of 13 inputs that ran before give byte-identical sequences, and the mCherry-family
protein, which crashed, now completes. Sequences for inputs that used to crash are
new by definition and depend on this version: keep the `meta.json`.

## v0.1.0 — 2026-09-28

First release. Extracted from a project pipeline and made host- and
vendor-generic; the algorithmic core is unchanged from the version it came from,
verified by regenerating four previously ordered fragments byte for byte under
the same seed.

**Added**

- One row per fragment as input, with per-row flanks, stop and 5'UTR, so a batch
  can mix constructs that need different flanks in a single run and a single
  order sheet.
- Host-specific constraints split three ways: generic patterns for every host,
  internal Shine–Dalgarno and strong RBS for any bacterium, Chi site and cryptic
  start codons for *E. coli* only.
- `<prefix>.meta.json` records this package's version and its output-determining
  dependencies alongside every parameter.
- READMEs in English and Japanese, with the pitfalls that cost real time, and
  `docs/design.md` with the rationale and the rejected alternatives.

**Fixed**

- Padding no longer lands in translated sequence. The pad sits between the ORF
  and the 3' flank, which is correct only when the ORF carries its own stop;
  with an empty `stop` the filler would be translated. The tool cannot know
  where inside a verbatim flank the stop falls, so it refuses and says how to
  resolve it.
- `--out-prefix` no longer overwrites `--input` when the two names coincide.
- The relaxed retry keeps the GC bands. Upstream deep-copies its constraint list
  before appending them, so its retry silently drops them — 50-bp windows
  reached 72 % GC that way.

**Known limits**

- The complexity gate is implemented for one vendor's API; the local heuristics
  are vendor-agnostic, the threshold is not.
- Eukaryotic hosts get the CAI table but none of the handling they actually need
  (Kozak context, CpG, splice sites), and the 5' stage should be switched off
  there. See "Initiation, not stability" in the README.
