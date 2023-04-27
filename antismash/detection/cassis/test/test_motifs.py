# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import os
from shutil import copy

from Bio.Seq import Seq

from antismash.common import path
from antismash.detection.cassis.promoters import Promoter, CombinedPromoter
from antismash.detection.cassis.motifs import generate_motifs, Motif, filter_meme_results, filter_fimo_results

from .test_cassis import CassisTestCore, convert_newline


def read_generated_expected_file(generated_file, expected_file):
    """Read generated and expected files and save to string variables"""
    generated_string = ""
    with open(generated_file, encoding="utf-8") as handle:
        generated_string = handle.read()
    generated_string = convert_newline(generated_string.rstrip())

    expected_string = ""
    with open(path.get_full_path(__file__, "data", expected_file), encoding="utf-8") as handle:
        expected_string = handle.read()
    expected_string = convert_newline(expected_string.rstrip())

    return [generated_string, expected_string]


class TestMotifs(CassisTestCore):
    def test_get_promoter_sets(self):
        meme_dir = os.path.join(self.options.output_dir, "meme")
        anchor_promoter = 5
        promoters = [Promoter("gene1", 1, 1, seq=Seq("acgtacgtacgtacgt")),
                     Promoter("gene2", 2, 2, seq=Seq("acgtacgtacgtacgt")),
                     CombinedPromoter("gene3", "gene4", 3, 4, seq=Seq("acgtacgtacgtacgt")),
                     Promoter("gene5", 5, 5, seq=Seq("acgtacgtacgtacgt")),
                     Promoter("gene6", 6, 6, seq=Seq("acgtacgtacgtacgt")),
                     # promoter with index=5 --> anchor promoter
                     Promoter("gene7", 7, 7, seq=Seq("acgtacgtacgtacgt")),
                     Promoter("gene8", 8, 8, seq=Seq("acgtacgtacgtacgt")),
                     Promoter("gene9", 9, 9, seq=Seq("acgtacgtacgtacgt"))]

        expected_motifs = [Motif(plus, minus) for plus in range(3) for minus in range(3-plus, 6)]
        self.assertEqual(generate_motifs(meme_dir, anchor_promoter, promoters), expected_motifs)

    def test_filter_meme_results(self):
        meme_dir = os.path.join(self.options.output_dir, "meme")
        anchor = "AFUA_6G09660"
        promoter_sets = [Motif(0, 3)]
        motif = Motif(0, 3, score=3.9e+003)
        motif.seqs = ["TTTCGACCCGTC",
                      "TTTCAAACCGTC",
                      "TTTTGATTCGTC",
                      "TTTTGACCGGTC",
                      "TTTTAGACGGTC",
                      "TTTTACCTCGTC",
                      "TCTCGATCCGTC",
                      "TTTCTATCCGTT",
                      "TTTTGGACCGCC",
                      "ATTTGGCCTGTC",
                      "TGTTGTCTCGTC",
                      "TTTGAGGCCGTC",
                      "TTGTATTCTGTC",
                      "TTTCTTCCTGTT"]
        expected_motifs = [motif]

        # this is a "real" MEME output file, I was too lazy to create my own fake XML file
        source = path.get_full_path(__file__, "data", "real_meme.xml")
        target = os.path.join(meme_dir, "+00_-03")
        if not os.path.exists(target):
            os.makedirs(target)
        copy(source, os.path.join(target, "meme.xml"))  # overwrite meme.xml if exists

        self.assertEqual(list(filter_meme_results(meme_dir, promoter_sets, anchor)), expected_motifs)
        binding_sites, expected_binding_sites = read_generated_expected_file(
            os.path.join(meme_dir, "+00_-03", "binding_sites.fasta"), "expected_binding_sites.fasta")
        self.assertEqual(binding_sites, expected_binding_sites)

    def test_filter_fimo_results(self):
        fimo_dir = os.path.join(self.options.output_dir, "fimo")
        motifs = [Motif(0, 3)]
        # gene2 will be the anchor promoter
        anchor_promoter = 1
        promoters = []
        for i in range(1, 16):
            promoters.append(Promoter(f"gene{i}", i * 10, i * 10 + 4))
        # need certain amount of promoters, otherwise the proportion of
        # promoters with a motif (motif frequency) will be too high --> error
        expected_motifs = [Motif(0, 3, hits={"gene1": 1, "gene2": 2})]

        # fake FIMO output file, corresponding to expected_motifs
        source = path.get_full_path(__file__, "data", "fake_short_fimo.txt")
        target = os.path.join(fimo_dir, "+00_-03")
        if not os.path.exists(target):
            os.makedirs(target)
        copy(source, os.path.join(target, "fimo.txt"))  # overwrite fimo.txt if exists

        found_motifs = filter_fimo_results(motifs, fimo_dir, promoters, anchor_promoter)
        assert found_motifs == expected_motifs
        bs_per_promoter, expected_bs_per_promoter = read_generated_expected_file(
            os.path.join(target, "bs_per_promoter.csv"), "expected_bs_per_promoter.csv")
        self.assertEqual(bs_per_promoter, expected_bs_per_promoter)
