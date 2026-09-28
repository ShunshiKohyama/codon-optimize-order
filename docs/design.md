# Design rationale

Why this tool makes the choices it does, and what was rejected on the way. Read
this before changing a parameter: most of the defaults are the *conclusion* of an
argument, and several are deliberately not the obvious value.

## The problem

Ordering synthetic genes has two failure modes that pull in opposite directions.
Optimise hard for codon usage and you produce repeat-rich, structured sequences
the vendor will not build. Optimise for manufacturability alone and you ship
genes that express badly. The manual loop — optimise on the vendor's site, read
the complexity colour, hand-fix the ones that fail — resolves this by trial and
error, unreproducibly, one gene at a time.

## Decision

**Optimise locally with DnaChisel; gate on the vendor's real complexity score;
ramp the repeat penalty until the gate clears.**

1. Reverse-translate with DnaChisel: objectives `MaximizeCAI`, `AvoidHairpins`,
   `MinimizeNumKmers(k=8)`; constraints on homopolymers, polymerase pausing
   sites, G-quadruplex-like motifs, and — for a bacterial host — internal
   Shine–Dalgarno, strong RBS, Chi site and cryptic start codons, plus global
   and windowed GC bands.
2. Score candidates with **local heuristics first**, so only a locally clean
   candidate costs an API call.
3. Ask the vendor for the real complexity score. Accept below the threshold.
4. On failure, ramp `kmers_weight` 10 → 100 over up to 10 steps and re-optimise.
   **Repeat suppression is the knob that actually buys synthesisability** — not
   GC, not CAI.

### Settled parameters

| parameter | value | why |
|---|---|---|
| complexity gate | **< 7** | matches upstream's default and the vendor's green band; the hard refusal limit is far higher, so this is a comfort margin, not a cliff |
| one fragment per protein | **not fused** | a bicistronic fragment carrying two genes was considered and rejected: it halves the fragment count but couples two genes' failure modes, and any redesign of one forces a resynthesis of both |
| stop codon | appended and **frozen** during optimisation | a standalone gene needs its own stop; `stop=""` exists for a vector-supplied C-terminal fusion |
| restriction sites | **none by default** | which sites matter is a property of the cloning method, not the gene. Avoiding a site you do not use costs CAI for nothing |
| cloning flanks | **none by default**, verbatim when given | they must match the vector exactly, so they never enter the optimisation |
| short-gene padding | filler after the stop, up to `--min-length` | vendors have a minimum fragment length that short ORFs fall under |
| seed | `0` | the search is stochastic; a fixed seed makes a batch reproduce byte for byte |

## The two-stage design

**Stage 1 chooses the first 20 codons for 5′ accessibility and freezes them;
stage 2 optimises the rest.**

The order matters. Running one combined optimisation lets the repeat penalty
undo the 5′ choice: re-choosing only the 24 nt *downstream* of the initiation
footprint moved a fixed head's `dg_open_5p` by 9–10 kcal/mol — as much as the
head itself. So the frozen region is the whole folded window (20 codons), not
just the footprint (12 codons).

### Why 5′ structure and not CAI

**CAI does not predict expression; initiation does.** Kudla et al.
(*Science* 2009) ran a moving-window analysis over 154 synonymous GFP variants
and found that the folding energy of **−4…+37 nt** around the start codon
explained **44 % of the variance** in protein level (*r* = 0.66). Codon
adaptation explained far less. Notably, the best window did *not* overlap the
Shine–Dalgarno sequence, so this is a start-codon accessibility effect rather
than an SD-occlusion one.

That reframes what codon choice is *for*. It is not a lever on translation
elongation worth optimising to three decimal places; it is the only lever the
pipeline has on initiation.

### AT content: an intervention, not a criterion

An AT-rich 5′ end pairs weakly and so cannot form the structure the window is
scored for. But **AT content works as an intervention and fails as a selection
criterion**: capping GC over the first 12 codons reliably lowers structure,
while ranking candidates by their AT content does not rank them by
`dg_open_5p`. So the cap is enforced (`--max-gc-start 0.40`) and the ranking is
done on the folding energy itself.

The cap is **not** loosened to 0.45. One gene that struggles at 0.40 is not
evidence the cap is too tight — it is evidence that gene needs a different
head, which stage 1 will find. Loosening a global default to rescue one
sequence trades every other gene's initiation for it.

## Rejected

**A `dg_open_5p` target or threshold.** There is no calibrated value to set.
Kudla's *r* = 0.66 is a correlation, not a dose–response, and the number moves
with the fold window, the UTR and the ViennaRNA parameter set. What is
actionable is the **rank within a batch** (`dg_open_pct`): "the most structured
start of these 48" is a thing you can act on. An absolute threshold would
manufacture a decision boundary that the evidence does not support.

**A CAI floor on the frozen head.** Spending head codons on accessibility costs
CAI, and the cost is larger for short proteins because the frozen 20 codons are
a bigger fraction of them. A floor was considered and rejected: it would buy a
number this document has already argued not to read, at the price of the one
mechanism that does predict expression. The trade is accepted as the point of
the design.

**CAI as a covariate for downstream modelling.** It is recorded as
*provenance* — what the optimiser achieved — not as a feature. Treating it as an
explanatory variable would smuggle back in the assumption the design rejects.

**Splitting over-long fragments.** Fragments above `--max-length` are flagged,
not split. Splitting is a cloning decision with assembly consequences; the tool
that chose your codons should not silently change your assembly strategy.

**A `dg_open_5p` floor worth chasing.** Below roughly 2 kcal/mol the differences
between candidate heads stop being meaningful. Reporting more precision there
invites over-reading noise.

## Known limits

**Host coverage is coarse.** `--species` swaps the CAI table, but the
host-dependent constraints follow a three-way classification: *E. coli* gets
everything, other listed bacteria get the initiation patterns without the
*E. coli*-only Chi site and cryptic-start motifs, and anything else gets only
the generic patterns. A eukaryotic host additionally needs Kozak context, CpG
and splice-site handling, none of which is implemented — and stage 1 should be
switched off there, since cap-dependent scanning is not the mechanism Kudla
measured.

**The 5′ numbers need the assembled construct.** Without `utr5` the fold sees
the ORF alone. The decisive window is mostly intra-ORF, so this is a fair
approximation, but it cannot see structure formed across the UTR junction. The
vector context is the part that matters and the part the tool cannot infer.

**Padding placement is not general.** The filler sits between the ORF and the 3′
flank, which is only correct when the ORF carries its own stop. With `stop`
empty the tool refuses to pad rather than guess where inside a verbatim flank
the stop falls.

**The gate is vendor-specific.** The complexity score comes from one vendor's
API. The local heuristics are vendor-agnostic; the threshold is not.

## Credits

Follows the Baker lab's SAPP/DMX release (`JB/domesticator.py` + `JB/idt.py`,
MIT), published with Qian, Milles, Wicky, Ragotte et al., *Nature
Communications* (2026). Deviations from upstream:

- the relaxed retry keeps the GC bands (upstream deep-copies its constraint list
  before appending them, so its retry silently drops them — 50-bp windows
  reached 72 % GC that way);
- local pre-screening before any API call;
- the two-stage 5′-structure design, which upstream does not have.
