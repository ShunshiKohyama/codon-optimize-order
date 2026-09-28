# codon-order

Codon-optimize proteins for synthesis, and gate every fragment on whether the
vendor can actually build it.

Replaces the manual loop — optimise on the vendor's website, read the complexity
colour, redo the hard ones by hand — with one reproducible command. Give it a
table of proteins and their flanking sequences; get back ordered fragments, a
paste-ready order sheet, and a record of every parameter that produced them.

*[日本語版 → README.ja.md](README.ja.md)*

```bash
codon-order --input fragments.csv --out-prefix batch1
```

```
2 fragments; default host=e_coli; gate=complexity < 7.0
sfGFP_untagged  (238 aa, host=e_coli)
    stage 1: 2571 valid heads -> dG_open 3.1 kcal/mol, rank 1 feasible (20 codons frozen)
    sfGFP_untagged: 717 bp  CAI=0.9043  dG_open5'=3.1  complexity=4.0  pad=0  flags=''
...
wrote batch1.csv
wrote batch1.order.csv   <- paste into the vendor's bulk entry form
wrote batch1.meta.json
```

## Install

```bash
conda env create -f environment.yml && conda activate codon-order
# or, into an existing environment:
pip install -e .
```

Needs Python ≥ 3.11. ViennaRNA is a hard requirement for the 5′-structure stage;
without it that stage is skipped and the run still completes.

## Input

One row per fragment, so a batch can mix constructs that need different flanks —
a tagged protein and an untagged one, or both halves of a two-plasmid system — in
a single run and a single order sheet.

| column | required | meaning |
|---|---|---|
| `name` | ✔ | fragment name; becomes the vendor's `Name` column |
| `protein` | ✔ | amino-acid sequence to reverse-translate |
| `adapter5` | | prepended verbatim, never optimised |
| `adapter3` | | appended verbatim, never optimised |
| `stop` | | stop codon to append; **empty when a flank supplies it** |
| `utr5` | | real 5′UTR of the assembled construct, for folding only |
| `species` | | per-row codon table, if the batch is not one host |

Anything omitted falls back to the matching command-line flag, so a uniform batch
needs only `name` and `protein`. See `examples/fragments.csv`.

## What it actually does

**1. Optimise locally** with [DnaChisel](https://github.com/Edinburgh-Genome-Foundry/DnaChisel):
maximise CAI, minimise repeated k-mers and hairpins, and forbid the patterns that
break either synthesis or expression — homopolymers, polymerase pausing sites,
and for a bacterial host internal Shine–Dalgarno, strong RBS, Chi site and
cryptic start codons.

**2. Choose the start, then optimise the rest.** Stage 1 searches synonymous
codings of the first 20 codons for the *least structured* 5′ end and freezes
them; stage 2 optimises everything downstream with the head locked. Codon choice
barely moves expression through CAI. It moves it a great deal through 5′ mRNA
structure, and this is the only place in the pipeline that can steer it.

**3. Gate on manufacturability.** Score each candidate with cheap local
heuristics, then ask the vendor for the real complexity score. If the fragment is
at or above the threshold, raise the repeat penalty and re-optimise. Repeat
suppression is the knob that actually buys synthesisability.

The rationale for every setting, and what was rejected, is in
[`docs/design.md`](docs/design.md). Read it before changing a parameter.

## Reproducibility

DnaChisel's search is stochastic. `--seed` (default `0`) makes a batch reproduce
byte for byte, and `<prefix>.meta.json` records every parameter, the input path,
the versions and the gate outcome — enough to regenerate the exact fragments you
ordered months later.

## Complexity gate

Optional. Without credentials the optimisation still runs end to end and the
vendor columns are left empty; check the scores on the vendor's site before
ordering, or supply credentials:

```bash
cp .env.example .env     # then fill it in; .env is git-ignored
codon-order-complexity --check-auth
```

Currently implemented for IDT's SciTools Plus API (eBlock complexity screening).

## Pitfalls

These are the ones that cost real time. Each is a trap the tool cannot detect
for you.

### A flank that supplies the start codon

If the upstream partner of a fusion ends in `…ATG`, that ATG *is* the protein's
initiator methionine. Put it in `adapter5` and **start `protein` at residue 2**,
or you get two methionines. Nothing can optimise a start codon, so this costs
nothing — it just has to be stated.

Check your design by translating the assembled fragment before ordering. The
tool's own output is a single column; it will not notice a doubled residue.

### `stop` and padding interact

The pad that brings a short fragment up to `--min-length` sits between the ORF
and the 3′ flank, which is correct only when the ORF carries its own stop. With
`stop` empty — a flank supplies the stop, or the protein continues into a
C-terminal fusion — the same filler would be **translated**. The tool refuses to
pad in that case and tells you to lengthen the 3′ flank yourself. It cannot know
where inside a verbatim flank the stop falls.

### `--n-candidates 1` does not switch the 5′ stage off

It stops stage 2 from *re-picking* among solutions. Stage 1 still searches and
freezes the head. To disable the 5′ machinery entirely, use `--head-samples 0`.

### Turn stage 1 off when the fragment is not the translation start

An internal fragment of a fusion — everything downstream of an N-terminal tag —
has no initiation event of its own. Freezing its first 20 codons for structure
protects a region no ribosome loads on. The cost is small (a per cent or two of
CAI) but the benefit is zero, and the same goes for a eukaryotic host, where
initiation is cap-dependent scanning and the bacterial result does not transfer.
Use `--head-samples 0`.

### `dg_open_5p` is not calibrated

Treat the absolute value as meaningless and the **rank within your batch**
(`dg_open_pct`, 100 = most structured) as the thing to act on. The number moves
with the fold window, the UTR and the ViennaRNA parameter set. It also has a
floor: below roughly 2 kcal/mol the differences between candidates stop meaning
anything.

### Without `utr5` the 5′ numbers are provisional

The decisive window straddles the start codon, so folding the ORF alone is a fair
approximation, but it cannot see structure formed across the UTR junction. If you
know the assembled construct's transcription start, pass it.

### Vendor auth returns 500, not 401

IDT's token endpoint answers a *failed password grant* with
`500 {"Message":"An error has occurred."}` — the same response for a wrong
password and for a user that does not exist. A bad client gives
`400 {"error":"invalid_client"}` instead, so:

- `invalid_client` → your `IDT_CLIENT_ID` / `IDT_CLIENT_SECRET` are wrong
- `500` → the client is fine; the **account login** is being rejected

Do not debug this by retrying — the API tells you nothing more and repeated
failures risk locking the account. Verify the username and password by logging
into the website, and try the account username and the e-mail address (either may
be the one the API wants).

### Other hosts need more than `--species`

`--species` swaps the CAI table (`python_codon_tables` ships
`e_coli`, `b_subtilis`, `s_cerevisiae`, `h_sapiens`, `m_musculus`,
`d_melanogaster`, `c_elegans`, `g_gallus`). The host-dependent *constraints*
follow a coarser classification:

| host | what applies |
|---|---|
| `e_coli` | everything: generic + internal SD / strong RBS + Chi site + cryptic starts |
| other listed bacteria | generic + internal SD / strong RBS |
| anything else | generic patterns only |

For a eukaryote that is not enough on its own: you also want Kozak context,
CpG and splice-site handling, and you should pass `--head-samples 0`. Add the
species to `BACTERIAL_SPECIES` in `src/codon_order/optimize.py` if you are
working with a bacterium that is not listed.

## Credits

The recipe follows the Baker lab's SAPP/DMX release
([`JB/domesticator.py`](https://github.com/bwicky/SAPP_DMX) + `JB/idt.py`, MIT),
published with **Qian, Milles, Wicky, Ragotte et al., *Nature Communications*
(2026)**. The deviations from upstream are deliberate and are listed in
[`docs/design.md`](docs/design.md).

The 5′-structure criterion rests on **Kudla et al., *Science* 324:255 (2009)**,
whose moving-window analysis over 154 synonymous GFP variants found that the
folding energy of −4…+37 nt around the start codon explained 44 % of the
variance in protein level (*r* = 0.66) — and that the best window did *not*
overlap the Shine–Dalgarno sequence.

MIT licensed. See [LICENSE](LICENSE).
