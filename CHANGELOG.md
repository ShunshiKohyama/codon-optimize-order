# Changelog

This tool designs sequences people pay to have synthesised, so a version is not
decoration: the same parameters can produce different DNA after an upgrade.
`<prefix>.meta.json` records the version of this package and of everything that
moves the output (DnaChisel, python_codon_tables, ViennaRNA, numpy). **Pin a tag
when a batch matters**, and keep the `meta.json` with the order:

```bash
pip install "git+https://github.com/ShunshiKohyama/codon-optimize-order@v0.1.0"
```

Versions follow [semantic versioning](https://semver.org/) with one addition:
**a change that alters designed sequences under unchanged parameters is a major
change**, even if no interface moves.

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
