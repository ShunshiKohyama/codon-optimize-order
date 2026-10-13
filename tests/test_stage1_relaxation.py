"""Stage 1 (head selection) must be able to give up the cryptic-start patterns, as stage 2 does.

Run:  PYTHONPATH=src python -m unittest discover -s tests -v
"""
import argparse
import random
import unittest

import dnachisel as dc
from dnachisel import NoSolutionError

from codon_order import optimize as core

# First 20 residues of the mCherry family (...GEEDNM...): every synonymous head matches a cryptic-start pattern.
MCHERRY_HEAD = "MVSKGEEDNMAIIKEFMRFK"
MCHERRY_90 = ("MVSKGEEDNMAIIKEFMRFKVHMEGSVNGHEFEIEGEGEGRPYEGTQTAKLKVTKGGPLPFAWDILSPQFMYGSKAYVKHPADIPDYLKLSFPEGFKWERVM")
# sfGFP, from examples/fragments.csv: some heads pass the full pattern set.
SFGFP_HEAD = "MSKGEELFTGVVPILVELDG"


def make_args(**over):
    base = dict(n_steps=2, starting_kmers_weight=10.0, adapter5="", adapter3="", head_samples=800, utr5="", species="e_coli",
                avoid=[], seed=0, stop="TAA", max_attempts=20, max_gc_start=0.40, fold_orf_nt=core.FOLD_ORF_NT,
                min_length=300, max_length=1500, n_candidates=1)
    base.update(over)
    return argparse.Namespace(**base)


class HeadRelaxation(unittest.TestCase):
    def test_unsatisfiable_head_is_returned_relaxed(self):
        got = core.optimize_head(MCHERRY_HEAD, n_samples=800, rng=random.Random(0))
        self.assertTrue(got)
        for c in got:
            self.assertTrue(c["alt_start_relaxed"])
            self.assertEqual(dc.translate(c["head"]), MCHERRY_HEAD)

    def test_old_all_or_nothing_behaviour_is_available(self):
        with self.assertRaises(NoSolutionError):
            core.optimize_head(MCHERRY_HEAD, n_samples=800, rng=random.Random(0), relax_alt_start=False)

    def test_relaxed_heads_still_obey_every_other_pattern(self):
        regexes = core._forbidden_regexes("e_coli", (), with_alt_start=False)
        for c in core.optimize_head(MCHERRY_HEAD, n_samples=800, rng=random.Random(0)):
            self.assertFalse(any(r.search(c["head"]) for r in regexes), c["head"])

    def test_nothing_changes_when_a_head_passes_the_full_set(self):
        a = core.optimize_head(SFGFP_HEAD, n_samples=800, rng=random.Random(0))
        b = core.optimize_head(SFGFP_HEAD, n_samples=800, rng=random.Random(0), relax_alt_start=False)
        self.assertTrue(a)
        self.assertEqual([c["head"] for c in a], [c["head"] for c in b])
        self.assertEqual([c["dg_open_5p"] for c in a], [c["dg_open_5p"] for c in b])
        self.assertFalse(any(c["alt_start_relaxed"] for c in a))

    def test_non_ecoli_hosts_are_not_relaxed(self):
        # the cryptic-start patterns are E. coli only, so there is nothing to relax elsewhere
        got = core.optimize_head(MCHERRY_HEAD, species="b_subtilis", n_samples=800, rng=random.Random(0))
        self.assertTrue(all(c["alt_start_relaxed"] is False for c in got))


class StageTwo(unittest.TestCase):
    def test_alt_start_false_skips_the_strict_rung(self):
        head = core.optimize_head(MCHERRY_HEAD, n_samples=800, rng=random.Random(0))[0]["head"]
        dna, alt_relaxed, _ = core.reverse_translate(MCHERRY_90, frozen_head=head, alt_start=False, max_tries=4)
        self.assertTrue(alt_relaxed)
        self.assertTrue(dna.startswith(head))
        self.assertEqual(dc.translate(dna[:-3]), MCHERRY_90)

    def test_default_ladder_is_unchanged(self):
        # with alt_start left at its default the first rung is still the strict one
        dna, alt_relaxed, _ = core.reverse_translate(MCHERRY_90[:30], max_tries=4)
        self.assertEqual(dc.translate(dna[:-3]), MCHERRY_90[:30])


class OptimizeGene(unittest.TestCase):
    def test_relaxed_head_is_used_and_reported(self):
        rec = core.optimize_gene("mcherry_like", MCHERRY_90, make_args(), None, random.Random(0))
        self.assertEqual(rec["head_codons"], core.HEAD_CODONS)
        self.assertTrue(rec["alt_start_relaxed"])
        self.assertEqual(dc.translate(rec["sequence"][:len(MCHERRY_90) * 3]), MCHERRY_90)

    def test_falls_back_to_single_stage_when_stage_one_finds_nothing(self):
        original = core.optimize_head

        def none_found(*a, **k):
            raise NoSolutionError("no valid head in 800 samples", None)

        core.optimize_head = none_found
        try:
            rec = core.optimize_gene("no_head", MCHERRY_90, make_args(), None, random.Random(0))
        finally:
            core.optimize_head = original
        self.assertEqual(rec["head_codons"], 0)
        self.assertEqual(dc.translate(rec["sequence"][:len(MCHERRY_90) * 3]), MCHERRY_90)


if __name__ == "__main__":
    unittest.main()
