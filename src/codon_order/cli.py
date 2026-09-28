"""Command line: read a table of proteins, write ordered fragments.

Input is one row per fragment, so a batch can mix constructs that need different
flanks — a tagged protein and an untagged one, or the two halves of a two-plasmid
system — in a single run and a single order sheet.  Columns:

==============  ========================================================
``name``        fragment name; becomes the vendor's Name column (required)
``protein``     amino-acid sequence to reverse-translate (required)
``adapter5``    prepended verbatim, never optimised (optional)
``adapter3``    appended verbatim, never optimised (optional)
``stop``        stop codon to append to the ORF; empty when a flank
                supplies it, or when the protein continues into a
                C-terminal fusion (optional, default ``--stop``)
``utr5``        real 5'UTR of the assembled construct, transcription
                start to the A of ATG, for folding only (optional)
``species``     per-row codon table, if the batch is not one host
                (optional, default ``--species``)
==============  ========================================================

Anything not given in a row falls back to the matching command-line flag, so a
uniform batch needs only ``name`` and ``protein``.

**If a flank supplies the start codon** — a fusion where the upstream partner
ends in the linker and the ATG — put the protein's own first residue in the
flank and start ``protein`` at residue 2.  Nothing can optimise a start codon,
so this is not a loss; it just has to be stated explicitly.
"""

from __future__ import annotations

import argparse
import json
import platform
import random
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from . import __version__
from . import optimize as core
from .idt import ComplexityClient, credentials_available, load_env


def dependency_versions() -> dict[str, str]:
    """Versions of everything that can move the output.

    DnaChisel drives the search, python_codon_tables supplies the CAI weights,
    ViennaRNA supplies the folding parameters and numpy seeds part of the search.
    A batch is only reproducible against the same set.
    """
    out: dict[str, str] = {}
    for mod, attr in (("dnachisel", "__version__"), ("python_codon_tables", "__version__"),
                      ("RNA", "__version__"), ("numpy", "__version__")):
        try:
            out[mod] = str(getattr(__import__(mod), attr, "unknown"))
        except ImportError:
            out[mod] = "not installed"
    return out

REQUIRED = ("name", "protein")
OPTIONAL = ("adapter5", "adapter3", "stop", "utr5", "species")


class RowArgs:
    """``args`` with this row's flanks patched in, for the per-gene driver."""

    def __init__(self, base, row):
        self.__dict__.update(vars(base))
        for field in OPTIONAL:
            value = row.get(field)
            if value is not None and str(value) != "nan":
                setattr(self, field, str(value))


def read_table(path: Path) -> pd.DataFrame:
    sep = "\t" if path.suffix.lower() in (".tsv", ".tab") else ","
    df = pd.read_csv(path, sep=sep, dtype=str).fillna("")
    df.columns = [c.strip().lower() for c in df.columns]
    missing = [c for c in REQUIRED if c not in df.columns]
    if missing:
        raise SystemExit(f"{path}: missing required column(s): {', '.join(missing)}")
    unknown = set(df.columns) - set(REQUIRED) - set(OPTIONAL)
    if unknown:
        print(f"note: ignoring unrecognised column(s): {', '.join(sorted(unknown))}")
    df["protein"] = df["protein"].str.strip().str.rstrip("*").str.upper()
    empty = df[(df["name"] == "") | (df["protein"] == "")]
    if len(empty):
        raise SystemExit(f"{path}: {len(empty)} row(s) have an empty name or protein")
    dupes = df["name"][df["name"].duplicated()].tolist()
    if dupes:
        raise SystemExit(f"{path}: duplicate name(s): {', '.join(dupes)}")
    return df


def check_padding(name: str, row_args, orf_len: int) -> None:
    """Refuse to pad a fragment whose stop codon is not at the end of the ORF.

    The pad sits between the ORF and the 3' flank, which is correct only when the
    ORF already carries its stop: the filler then lands in untranslated sequence.
    With ``stop`` empty the flank carries the stop — or the protein continues into
    a C-terminal fusion — and the same filler would be *translated*, appending
    junk residues or inserting them into the fusion.  The tool cannot know where
    inside a verbatim flank the stop falls, so it stops instead of guessing.
    """
    overhead = len(row_args.adapter5) + len(row_args.adapter3)
    if max(0, row_args.min_length - (orf_len + overhead)) and not row_args.stop:
        raise SystemExit(
            f"{name}: needs padding to reach --min-length {row_args.min_length}, but "
            f"'stop' is empty for this row, so the pad would land inside translated "
            f"sequence.\nFix it in one of these ways:\n"
            f"  - lengthen the 3' flank yourself, after its stop codon\n"
            f"  - set 'stop' for this row, if the ORF really should end there\n"
            f"  - lower --min-length if the vendor allows a shorter fragment"
        )


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True, type=Path,
                    help="CSV/TSV of fragments (see the columns above).")
    ap.add_argument("--out-prefix", default="codon_order",
                    help="Output files are <prefix>.csv, <prefix>.order.csv, "
                         "<prefix>.meta.json.")
    ap.add_argument("--species", default="e_coli",
                    help="Default codon table (a python_codon_tables key).")
    ap.add_argument("--idt-score", type=float, default=7.0,
                    help="Accept a fragment below this Total Complexity Score.")
    ap.add_argument("--skip-idt-query", action="store_true",
                    help="Optimise locally only; leave the vendor columns empty.")
    ap.add_argument("--n-steps", type=int, default=10)
    ap.add_argument("--starting-kmers-weight", type=float, default=10.0)
    ap.add_argument("--max-attempts", type=int, default=20)
    ap.add_argument("--n-candidates", type=int, default=6,
                    help="Solutions drawn per ramp step; the least-structured "
                         "start wins. 1 keeps stage 2 from re-picking (it does "
                         "NOT switch stage 1 off -- that is --head-samples 0).")
    ap.add_argument("--utr5", default="",
                    help="Default 5'UTR for folding: transcription start up to, "
                         "not including, the ATG of the ASSEMBLED construct. Not "
                         "a flank -- see 'adapter5'.")
    ap.add_argument("--head-samples", type=int, default=8000,
                    help=f"Stage 1: synonymous heads sampled for the first "
                         f"{core.HEAD_CODONS} codons. 0 disables stage 1 -- do that "
                         "for a eukaryotic host, and for any fragment that is not "
                         "itself the translation start.")
    ap.add_argument("--max-gc-start", type=float, default=0.40)
    ap.add_argument("--fold-orf-nt", type=int, default=core.FOLD_ORF_NT)
    ap.add_argument("--min-length", type=int, default=300,
                    help="Vendor minimum; shorter fragments are padded after the stop.")
    ap.add_argument("--max-length", type=int, default=1500,
                    help="Vendor maximum; longer fragments are flagged, not split.")
    ap.add_argument("--stop", default="TAA", choices=[*core.STOP_CODONS, ""],
                    help="Default stop codon ('' when a flank supplies it).")
    ap.add_argument("--avoid", nargs="*", default=list(core.DEFAULT_ENZYME_SITES),
                    help="Restriction sites to keep out, both strands. "
                         "For BsaI: --avoid GGTCTC GAGACC")
    ap.add_argument("--adapter5", default="", help="Default 5' flank, verbatim.")
    ap.add_argument("--adapter3", default="", help="Default 3' flank, verbatim.")
    ap.add_argument("--seed", type=int, default=0,
                    help="DnaChisel's search is stochastic; the same seed "
                         "reproduces a batch byte for byte.")
    ap.add_argument("--limit", type=int, default=0, help="Only the first N rows.")
    args = ap.parse_args(argv)

    out_csv = Path(f"{args.out_prefix}.csv").resolve()
    if out_csv == args.input.resolve():
        raise SystemExit(
            f"--out-prefix would write over the input file ({args.input}).\n"
            f"Pick a different prefix: the run writes <prefix>.csv, "
            f"<prefix>.order.csv and <prefix>.meta.json."
        )

    df = read_table(args.input)
    if args.limit:
        df = df.head(args.limit)

    client = None
    if not args.skip_idt_query:
        if credentials_available():
            client = ComplexityClient(load_env())
        else:
            print("note: no vendor credentials in .env — optimising locally only "
                  "(same as --skip-idt-query). The sequences are still produced; "
                  "only the vendor complexity columns stay empty. See the README "
                  "section 'Complexity gate'.")

    if args.head_samples > 0 and not core.is_bacterial(args.species):
        print(f"note: stage 1 picks the least-structured start, which rests on a "
              f"bacterial initiation result (Kudla 2009). '{args.species}' is not "
              f"a listed bacterium — consider --head-samples 0. See docs/design.md.")

    rng = random.Random(args.seed)
    gate = "off" if client is None else f"complexity < {args.idt_score}"
    print(f"{len(df)} fragments; default host={args.species}; gate={gate}\n")

    rows, started = [], time.time()
    for _, row in df.iterrows():
        name = row["name"]
        row_args = RowArgs(args, row)
        aa = row["protein"]
        print(f"{name}  ({len(aa)} aa, host={row_args.species})")
        check_padding(name, row_args, len(aa) * 3 + len(row_args.stop))
        record = core.optimize_gene(name, aa, row_args, client, rng)
        record.update({
            "aa_len": len(aa),
            "species": row_args.species,
            "adapter5_bp": len(row_args.adapter5),
            "adapter3_bp": len(row_args.adapter3),
            "stop": row_args.stop,
        })
        rows.append(record)
        print(f"    {name}: {record['length_bp']} bp  CAI={record['cai']}  "
              f"dG_open5'={record['dg_open_5p']}  "
              f"complexity={record.get('idt_score')}  "
              f"pad={record['pad_bp']}  flags='{record['local_flags']}'")

    elapsed = time.time() - started
    column_order = [
        "name", "species", "aa_len", "orf_bp", "pad_bp", "adapter5_bp",
        "adapter3_bp", "stop", "length_bp", "cai", "dg_open_5p", "dg_open_pct",
        "at_frac_5p", "mfe_5p", "utr5_bp", "gc", "gc_win50_min", "gc_win50_max",
        "homopolymer_max", "at_run_max", "gc_run_max", "nonunique_8mer_frac",
        "idt_score", "idt_pass", "local_flags", "kmers_weight", "steps_used",
        "alt_start_relaxed", "at_rich_relaxed", "head_codons", "head_rank",
        "sequence",
    ]
    out = pd.DataFrame(rows)
    out = out[[c for c in column_order if c in out.columns and c != "dg_open_pct"]]
    # An absolute dg_open_5p has no calibrated meaning; its rank within the batch
    # does. 100 = the most structured start of this batch, which is what you act on.
    out["dg_open_pct"] = (out["dg_open_5p"].rank(pct=True) * 100).round(0)
    out = out[[c for c in column_order if c in out.columns]]

    prefix = Path(args.out_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(f"{prefix}.csv", index=False)
    out[["name", "sequence"]].rename(
        columns={"name": "Name", "sequence": "Sequence"}
    ).to_csv(f"{prefix}.order.csv", index=False)
    meta = {
        "written": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "elapsed_s": round(elapsed, 1),
        "input": str(args.input),
        "fragments": len(out),
        # Parameters alone do not pin an output: the search, the codon tables and
        # the folding parameters all live in versioned code. Record them, or
        # "reproducible" is only true until something is upgraded.
        "codon_order_version": __version__,
        "dependencies": dependency_versions(),
        "python": platform.python_version(),
        "platform": platform.platform(),
        "parameters": {k: v for k, v in vars(args).items() if k != "input"},
        "gate": {"applied": client is not None,
                 "threshold": args.idt_score,
                 "passed": int(out["idt_pass"].sum()) if client is not None else None},
    }
    meta["parameters"]["out_prefix"] = str(args.out_prefix)
    Path(f"{prefix}.meta.json").write_text(json.dumps(meta, indent=2, default=str) + "\n")

    print(f"\n{len(out)} fragments in {elapsed:.0f}s; "
          f"CAI {out.cai.min():.2f}-{out.cai.max():.2f}; "
          f"padded {int((out.pad_bp > 0).sum())}; "
          f"locally flagged {int((out.local_flags != '').sum())}")
    long = out[out.length_bp > args.max_length]
    if len(long):
        print(f"WARNING: {len(long)} fragment(s) exceed --max-length "
              f"{args.max_length}: {', '.join(long.name)}")
    print(f"wrote {prefix}.csv")
    print(f"wrote {prefix}.order.csv   <- paste into the vendor's bulk entry form")
    print(f"wrote {prefix}.meta.json")


if __name__ == "__main__":
    sys.exit(main())
