# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""Test suite for the cassis cluster detection plugin"""

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import os
from shutil import copy
from tempfile import TemporaryDirectory
import unittest

from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation
from minimock import mock, restore

from antismash.common import path, secmet
from antismash.common.test import helpers
from antismash.config import build_config, destroy_config
from antismash.detection import cassis
from antismash.detection.cassis.runners import run_fimo


# helper methods
def convert_newline(string):
    """Convert all line endings to \n for OS independency"""
    return "\n".join(string.splitlines())


def read_generated_expected_file(generated_file, expected_file):
    """Read generated and expected files and save to string variables"""
    generated_string = ""
    with open(generated_file) as handle:
        generated_string = handle.read()
    generated_string = convert_newline(generated_string.rstrip())

    expected_string = ""
    with open(path.get_full_path(__file__, "data", expected_file)) as handle:
        expected_string = handle.read()
    expected_string = convert_newline(expected_string.rstrip())

    return [generated_string, expected_string]


def create_fake_record():
    """Set up a fake sequence record"""
    seq_record = helpers.DummyRecord(seq=Seq("acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta" * 196))
    seq_record.name = "test"
    locations = [FeatureLocation(100, 300, strand=1), FeatureLocation(101, 299, strand=-1),
                 FeatureLocation(250, 350, strand=1), FeatureLocation(500, 1000, strand=1),
                 FeatureLocation(1111, 1500, strand=-1), FeatureLocation(2000, 2200, strand=-1),
                 FeatureLocation(2999, 4000, strand=1), FeatureLocation(4321, 5678, strand=1),
                 FeatureLocation(6660, 9000, strand=-1)]
    for i in range(9):
        cds = helpers.DummyCDS(locus_tag="gene" + str(i+1))
        cds.location = locations[i]
        seq_record.add_cds_feature(cds)
        seq_record.add_gene(secmet.Gene(locations[i], locus_tag="gene" + str(i+1)))
        if i == 3 or i == 5:
            cds.sec_met = secmet.qualifiers.SecMetQualifier({"faked"}, [])
            cds.gene_functions.add(secmet.feature.GeneFunction.CORE, "testtool", "dummy")

    return seq_record


class TestCassisMethods(unittest.TestCase):
    def setUp(self):
        self.tempdir = TemporaryDirectory(prefix="as_cassis")
        self.options = build_config(["--cpus", "2", "--output-dir", self.tempdir.name],
                                    isolated=True, modules=[cassis])

    def tearDown(self):
        destroy_config()
        self.tempdir.cleanup()
        for subdir in ["meme", "fimo"]:
            assert not os.path.exists(path.get_full_path(__file__, subdir))
        assert not os.path.isfile(path.get_full_path(__file__, "data", "test_promoter_positions.csv"))
        assert not os.path.isfile(path.get_full_path(__file__, "data", "test_promoter_sequences.fasta"))

    def test_promoter_range(self):
        self.assertEqual(len(cassis.PROMOTER_RANGE), 250)

    def test_get_anchor_gene_names(self):
        anchor_genes = ["gene4", "gene6"]
        seq_record = create_fake_record()
        self.assertEqual(cassis.get_anchor_gene_names(seq_record), anchor_genes)

    def test_ignore_overlapping(self):
        expected_not_ignored = ["gene1", "gene4", "gene5", "gene6", "gene7", "gene8", "gene9"]
        expected_ignored = ["gene2", "gene3"]
        seq_record = create_fake_record()
        not_ignored, ignored = cassis.ignore_overlapping(seq_record.get_genes())

        self.assertEqual([x.locus_tag for x in ignored], expected_ignored)
        self.assertEqual([x.locus_tag for x in not_ignored], expected_not_ignored)

    def test_get_promoters(self):
        upstream_tss = 1000
        downstream_tss = 50
        seq_record = create_fake_record()
        genes, ignored_genes = cassis.ignore_overlapping(seq_record.get_genes())  # ignore ignored_genes

        # see cassis/promoterregions.png for details
        # [[start_prom1, end_prom1], [start_prom2, end_prom2], ...]
        expected_promoters = [
            [0, 150],
            [301, 550],
            [1450, 1999],
            [2150, 3049],
            [4001, 4371],
            [8950, 9603],
        ]

        promoters = cassis.get_promoters(seq_record, genes, upstream_tss, downstream_tss)
        self.assertEqual(list(map(lambda x: [x.start, x.end], promoters)), expected_promoters)
        cassis.write_promoters_to_file(self.options.output_dir, seq_record.name, promoters)
        # read expected files and save to string variable
        expected_sequences_file = ""
        with open(path.get_full_path(__file__, "data", "expected_promoter_sequences.fasta")) as handle:
            expected_sequences_file = handle.read()
        expected_sequences_file = convert_newline(expected_sequences_file.rstrip())

        expected_positions_file = ""
        with open(path.get_full_path(__file__, "data", "expected_promoter_positions.csv")) as handle:
            expected_positions_file = handle.read()
        expected_positions_file = convert_newline(expected_positions_file.rstrip())

        # read test files and save to string variable
        sequences_file = ""
        with open(os.path.join(self.options.output_dir, seq_record.name + "_promoter_sequences.fasta")) as handle:
            sequences_file = handle.read()
        sequences_file = convert_newline(sequences_file.rstrip())

        positions_file = ""
        with open(os.path.join(self.options.output_dir, seq_record.name + "_promoter_positions.csv")) as handle:
            positions_file = handle.read()
        positions_file = convert_newline(positions_file.rstrip())

        self.assertEqual(sequences_file, expected_sequences_file)
        self.assertEqual(positions_file, expected_positions_file)

    def test_get_anchor_promoter(self):
        anchor = "gene3"
        promoters = [cassis.Promoter("gene1", 1, 1),
                     cassis.Promoter("gene2", 2, 2),
                     cassis.CombinedPromoter("gene3", "gene4", 3, 4),
                     cassis.Promoter("gene5", 5, 5)]
        self.assertEqual(cassis.get_anchor_promoter_index(anchor, promoters), 2)

    def test_get_promoter_sets(self):
        meme_dir = os.path.join(self.options.output_dir, "meme")
        anchor_promoter = 5
        promoters = [cassis.Promoter("gene1", 1, 1, seq=Seq("acgtacgtacgtacgt")),
                     cassis.Promoter("gene2", 2, 2, seq=Seq("acgtacgtacgtacgt")),
                     cassis.CombinedPromoter("gene3", "gene4", 3, 4, seq=Seq("acgtacgtacgtacgt")),
                     cassis.Promoter("gene5", 5, 5, seq=Seq("acgtacgtacgtacgt")),
                     cassis.Promoter("gene6", 6, 6, seq=Seq("acgtacgtacgtacgt")),
                     # promoter with index=5 --> anchor promoter
                     cassis.Promoter("gene7", 7, 7, seq=Seq("acgtacgtacgtacgt")),
                     cassis.Promoter("gene8", 8, 8, seq=Seq("acgtacgtacgtacgt")),
                     cassis.Promoter("gene9", 9, 9, seq=Seq("acgtacgtacgtacgt"))]

        expected_promoter_sets = [cassis.Motif(plus, minus) for plus in range(3) for minus in range(3-plus, 6)]
        self.assertEqual(cassis.get_promoter_sets(meme_dir, anchor_promoter, promoters), expected_promoter_sets)

    def test_filter_meme_results(self):
        meme_dir = os.path.join(self.options.output_dir, "meme")
        anchor = "AFUA_6G09660"
        promoter_sets = [cassis.Motif(0, 3)]
        motif = cassis.Motif(0, 3, score=3.9e+003)
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

        self.assertEqual(list(cassis.filter_meme_results(meme_dir, promoter_sets, anchor)), expected_motifs)
        binding_sites, expected_binding_sites = read_generated_expected_file(
            os.path.join(meme_dir, "+00_-03", "binding_sites.fasta"), "expected_binding_sites.fasta")
        self.assertEqual(binding_sites, expected_binding_sites)

    def test_filter_fimo_results(self):
        fimo_dir = os.path.join(self.options.output_dir, "fimo")
        motifs = [cassis.Motif(0, 3)]
        # gene2 will be the anchor promoter
        anchor_promoter = 1
        promoters = []
        for i in range(1, 16):
            promoters.append(cassis.Promoter("gene%d" % i, i * 10, i * 10 + 4))
        # need certain amount of promoters, otherwise the proportion of
        # promoters with a motif (motif frequency) will be too high --> error
        expected_motifs = [cassis.Motif(0, 3, hits={"gene1": 1, "gene2": 2})]

        # fake FIMO output file, corresponding to expected_motifs
        source = path.get_full_path(__file__, "data", "fake_short_fimo.txt")
        target = os.path.join(fimo_dir, "+00_-03")
        if not os.path.exists(target):
            os.makedirs(target)
        copy(source, os.path.join(target, "fimo.txt"))  # overwrite fimo.txt if exists

        found_motifs = cassis.filter_fimo_results(motifs, fimo_dir, promoters, anchor_promoter)
        assert found_motifs == expected_motifs
        bs_per_promoter, expected_bs_per_promoter = read_generated_expected_file(
            os.path.join(target, "bs_per_promoter.csv"), "expected_bs_per_promoter.csv")
        self.assertEqual(bs_per_promoter, expected_bs_per_promoter)

    def test_get_islands(self):
        motifs = [cassis.Motif(0, 3, hits={"gene1": 1, "gene2": 2}),
                  cassis.Motif(0, 4, hits={"gene2": 3, "gene4": 2, "gene5": 1})]
        # gene2 will be the anchor promoter
        anchor_promoter = 1
        promoters = []
        for i in range(1, 7):
            promoters.append(cassis.Promoter("gene%d" % i, i * 10, i * 10 + 4))
        # resulting in 2 different islands (this example)
        # promoter (pos): 1 2 3 4 5 6
        # binding sites:  1 2 0 0 0 0
        # island:         |-|
        first_island = cassis.Island(promoters[0], promoters[1], motifs[0])
        # promoter (pos): 1 2 3 4 5 6
        # binding sites:  0 3 0 2 1 0
        # island:           |---|
        second_island = cassis.Island(promoters[1], promoters[4], motifs[1])
        expected_islands = [first_island, second_island]
        assert cassis.get_islands(anchor_promoter, motifs, promoters) == expected_islands

    def test_sort_by_abundance(self):
        islands = []

        # island 1: [gene1 -- gene2]
        motif = cassis.Motif(0, 3, score=3, hits={"gene1": 1, "gene2": 1})
        islands.append(cassis.Island(cassis.Promoter("gene1", 1, 1), cassis.Promoter("gene2", 2, 2), motif))
        # island 2: [gene2 -- gene5]
        motif = cassis.Motif(3, 0, score=2, hits={"gene2": 1, "gene3": 1, "gene4": 1, "gene5": 1})
        islands.append(cassis.Island(cassis.Promoter("gene2", 2, 2), cassis.Promoter("gene5", 5, 5), motif))
        # island 3: [gene1 -- gene5]
        motif = cassis.Motif(3, 3, score=1, hits={"gene1": 1, "gene2": 1, "gene3": 1, "gene4": 1, "gene5": 1})
        islands.append(cassis.Island(cassis.Promoter("gene1", 1, 1), cassis.Promoter("gene5", 5, 5), motif))

        # left border: 2x gene1, 1x gene2
        # right border: 2x gene5, 1x gene2

        expected_clusters = []
        # cluster 1: [gene1 -- gene5] --> abundance 2+2 (most abundant)
        start = cassis.ClusterMarker("gene1", cassis.Motif(3, 3, score=1))
        start.abundance = 2
        end = cassis.ClusterMarker("gene5", cassis.Motif(3, 3, score=1))
        end.abundance = 2
        expected_clusters.append(cassis.ClusterPrediction(start, end))
        # cluster 3: [gene2 -- gene5] --> abundance 1+2, score 2+1 (better/lower)
        start = cassis.ClusterMarker("gene2", cassis.Motif(3, 0, score=2))
        start.abundance = 1
        end = cassis.ClusterMarker("gene5", cassis.Motif(3, 3, score=1))
        end.abundance = 2
        expected_clusters.append(cassis.ClusterPrediction(start, end))
        # cluster 2: [gene1 -- gene2] --> abundance 2+1, score 1+3 (worse, higher)
        start = cassis.ClusterMarker("gene1", cassis.Motif(3, 3, score=1))
        start.abundance = 2
        end = cassis.ClusterMarker("gene2", cassis.Motif(0, 3, score=3))
        end.abundance = 1
        expected_clusters.append(cassis.ClusterPrediction(start, end))
        # cluster 4: [gene2 -- gene2] --> abundance 1+1
        start = cassis.ClusterMarker("gene2", cassis.Motif(3, 0, score=2))
        start.abundance = 1
        end = cassis.ClusterMarker("gene2", cassis.Motif(0, 3, score=3))
        end.abundance = 1
        expected_clusters.append(cassis.ClusterPrediction(start, end))
        # abundance: as high as possible
        # score: as low as possible

        self.assertEqual(cassis.sort_by_abundance(islands), expected_clusters)

    def test_check_cluster_predictions(self):
        seq_record = create_fake_record()
        promoters = [cassis.Promoter("gene1", 1, 5),
                     cassis.Promoter("gene2", 6, 10),
                     cassis.CombinedPromoter("gene3", "gene4", 11, 15)]
        ignored_genes = [  # see captured logging
            secmet.Gene(FeatureLocation(1, 5), locus_tag="gene5")
        ]
        clusters = [cassis.ClusterPrediction(cassis.ClusterMarker("gene1", cassis.Motif(3, 3, score=1)),
                                             cassis.ClusterMarker("gene4", cassis.Motif(3, 3, score=1)))]
        expected = [cassis.ClusterPrediction(cassis.ClusterMarker("gene1", cassis.Motif(3, 3, score=1)),
                                             cassis.ClusterMarker("gene4", cassis.Motif(3, 3, score=1)))]
        expected[0].start.promoter = "gene1"
        expected[0].end.promoter = "gene3+gene4"
        expected[0].genes = 4
        expected[0].promoters = 3

        assert cassis.check_cluster_predictions(clusters, seq_record, promoters, ignored_genes) == expected
        # -------------------- >> begin captured logging << --------------------
        # root: INFO: Best prediction (most abundant): 'gene1' -- 'gene4'
        # root: INFO: Upstream cluster border located at or next to sequence
        #               record border, prediction could have been truncated by record border
        # root: INFO: Downstream cluster border located at or next to sequence
        #               record border, prediction could have been truncated by record border
        # root: INFO: Gene 'gene2' is part of the predicted cluster, but it is
        #               overlapping with another gene and was ignored
        # root: INFO: Gene 'gene2' could have effected the cluster prediction
        # --------------------- >> end captured logging << ---------------------

    def test_cleanup_outdir(self):
        anchor_genes = ["gene1", "gene4"]
        cluster = cassis.ClusterPrediction(cassis.ClusterMarker("gene1", cassis.Motif(3, 3, score=1)),
                                           cassis.ClusterMarker("gene4", cassis.Motif(3, 3, score=1)))
        cluster.start.promoter = "gene1"
        cluster.end.promoter = "gene3+gene4"
        cluster.genes = 4
        cluster.promoters = 3
        cluster_predictions = {"gene1": [cluster]}

        # create some empty test dirs, which should be deleted during the test
        # prediction! --> keep!
        os.makedirs(os.path.join(self.options.output_dir, "meme", "gene1", "+03_-03"))
        # prediction! --> keep!
        os.makedirs(os.path.join(self.options.output_dir, "fimo", "gene1", "+03_-03"))
        # no prediction --> delete
        os.makedirs(os.path.join(self.options.output_dir, "meme", "gene1", "+04_-04"))
        # no prediction --> delete
        os.makedirs(os.path.join(self.options.output_dir, "fimo", "gene1", "+04_-04"))
        # no prediction --> delete
        os.makedirs(os.path.join(self.options.output_dir, "meme", "gene4", "+03_-03"))
        # no prediction --> delete
        os.makedirs(os.path.join(self.options.output_dir, "fimo", "gene4", "+03_-03"))
        # prediction for this gene, but not from this motif --> delete
        os.makedirs(os.path.join(self.options.output_dir, "meme", "gene4", "+04_-04"))
        # prediction for this gene, but not from this motif --> delete
        os.makedirs(os.path.join(self.options.output_dir, "fimo", "gene4", "+04_-04"))

        cassis.cleanup_outdir(anchor_genes, cluster_predictions, self.options)

        # assert kept directories
        self.assertTrue("gene1" in os.listdir(os.path.join(self.options.output_dir, "meme")))
        self.assertTrue("gene1" in os.listdir(os.path.join(self.options.output_dir, "fimo")))
        self.assertTrue("+03_-03" in os.listdir(os.path.join(self.options.output_dir, "meme", "gene1")))
        self.assertTrue("+03_-03" in os.listdir(os.path.join(self.options.output_dir, "fimo", "gene1")))

        # assert deleted directories
        self.assertTrue("gene4" not in os.listdir(os.path.join(self.options.output_dir, "meme")))
        self.assertTrue("gene4" not in os.listdir(os.path.join(self.options.output_dir, "fimo")))
        self.assertTrue("+04_-04" not in os.listdir(os.path.join(self.options.output_dir, "meme", "gene1")))
        self.assertTrue("+04_-04" not in os.listdir(os.path.join(self.options.output_dir, "fimo", "gene1")))


class TestMotifRepresentation(unittest.TestCase):
    def test_conversion(self):
        motif = cassis.Motif(3, 3)
        assert motif.pairing_string == "+03_-03"
        motif.plus = 4
        assert motif.pairing_string == "+04_-03"
        motif.minus = 2
        assert motif.pairing_string == "+04_-02"


class TestPromoters(unittest.TestCase):
    def test_promoter_id(self):
        assert cassis.Promoter("gene1", 1, 5).get_id() == "gene1"
        assert cassis.CombinedPromoter("gene1", "gene2", 1, 5).get_id() == "gene1+gene2"


class TestCassisStorageMethods(unittest.TestCase):
    def test_store_promoters(self):
        promoters = [cassis.Promoter("gene1", 10, 20, seq=Seq("cgtacgtacgt")),
                     cassis.Promoter("gene2", 30, 40, seq=Seq("cgtacgtacgt")),
                     cassis.CombinedPromoter("gene3", "gene4", 50, 60, seq=Seq("cgtacgtacgt"))]
        record_with_promoters = create_fake_record()
        cassis.store_promoters(promoters, record_with_promoters)  # add ("store") promoters to seq_record

        record_without_promoters = create_fake_record()  # just the same, without adding promoters

        # test if store_promoters changed any non-promoter feature (should not!)  # TODO

        # test promoter features
        expected_count = record_without_promoters.get_feature_count() + len(promoters)
        assert expected_count == record_with_promoters.get_feature_count()
        for i in range(len(promoters)):
            feature = record_with_promoters.get_generics()[i]
            assert feature.type == "promoter"
            assert feature.get_qualifier("seq") == ("cgtacgtacgt",)

        # especially test bidirectional promoter feature (third promoter, last feature)
        last_promoter = record_with_promoters.get_generics()[-1]
        assert last_promoter.get_qualifier("locus_tag") == ("gene3", "gene4")
        assert last_promoter.notes == ["bidirectional promoter"]

    def test_store_clusters(self):
        # this test is similar to test_store_promoters
        anchor = "gene3"

        start_marker = cassis.ClusterMarker("gene1", cassis.Motif(3, 3, score=1))
        start_marker.promoter = "gene1"
        start_marker.abundance = 2
        end_marker = cassis.ClusterMarker("gene4", cassis.Motif(3, 3, score=1))
        end_marker.promoter = "gene3+gene4"
        assert end_marker.abundance == 1
        first_cluster = cassis.ClusterPrediction(start_marker, end_marker)
        first_cluster.promoters = 3
        first_cluster.genes = 4

        start_marker = cassis.ClusterMarker("gene1", cassis.Motif(4, 4, score=1))
        start_marker.promoter = "gene1"
        assert start_marker.abundance == 1
        end_marker = cassis.ClusterMarker("gene5", cassis.Motif(4, 4, score=1))
        end_marker.promoter = "gene5"
        assert end_marker.abundance == 1
        second_cluster = cassis.ClusterPrediction(start_marker, end_marker)
        second_cluster.promoters = 3
        second_cluster.genes = 4

        clusters = [first_cluster, second_cluster]

        record_with_clusters = create_fake_record()
        record_without_clusters = create_fake_record()  # just the same, without adding clusters

        borders = cassis.create_cluster_borders(anchor, clusters, record_with_clusters)
        assert record_with_clusters.get_feature_count() == record_without_clusters.get_feature_count()

        for border in borders:
            record_with_clusters.add_cluster_border(border)

        # test if store_clusters changed any non-cluster feature (should not!)  # TODO

        # test cluster features
        assert record_without_clusters.get_feature_count() + len(clusters) == record_with_clusters.get_feature_count()
        for i, cluster in enumerate(clusters):
            cluster_border = record_with_clusters.get_cluster_borders()[i]
            self.assertEqual(cluster_border.type, "cluster_border")
            self.assertEqual(cluster_border.tool, "cassis")
            self.assertEqual(cluster_border.get_qualifier("anchor"), (anchor,))
            self.assertEqual(cluster_border.get_qualifier("genes"), (cluster.genes,))
            self.assertEqual(cluster_border.get_qualifier("promoters"), (cluster.promoters,))
            self.assertEqual(cluster_border.get_qualifier("gene_left"), (cluster.start.gene,))
            self.assertEqual(cluster_border.get_qualifier("gene_right"), (cluster.end.gene,))
            # don't test all feature qualifiers, only some


class TestCassisUtils(unittest.TestCase):
    def setUp(self):
        self.tempdir = TemporaryDirectory(prefix="as_cassis")
        self.options = build_config(["--cpus", "2", "--output-dir", self.tempdir.name],
                                    isolated=True, modules=[cassis])

    def tearDown(self):
        destroy_config()
        self.tempdir.cleanup()

    def test_run_fimo(self):
        seq_record = create_fake_record()
        meme_dir = os.path.join(self.options.output_dir, "meme")
        fimo_dir = os.path.join(self.options.output_dir, "fimo")
        fimo_subdirs = ["+00_-03", "+03_-00"]

        # fimo input (i)
        # --> sequences to search in
        copy(
            path.get_full_path(__file__, "data", "expected_promoter_sequences.fasta"),
            os.path.join(self.options.output_dir, seq_record.name + "_promoter_sequences.fasta")
        )
        # fimo input (ii)
        # --> meme output (meme.html; motifs to search for)
        # --> file derived from meme output (binding_sites.fasta;
        #      only additional information for users; still checked by fimo)
        for subdir in fimo_subdirs:
            if not os.path.exists(os.path.join(meme_dir, subdir)):
                os.makedirs(os.path.join(meme_dir, subdir))
            copy(path.get_full_path(__file__, "data", "fake_meme.html"),
                 os.path.join(meme_dir, subdir, "meme.html"))
            copy(path.get_full_path(__file__, "data", "fake_binding_sites.fasta"),
                 os.path.join(meme_dir, subdir, "binding_sites.fasta"))

        run_fimo(meme_dir, fimo_dir, seq_record, self.options, verbose=True)

        for subdir in fimo_subdirs:
            fimo_result, expected_fimo_result = read_generated_expected_file(
                                                    os.path.join(self.options.output_dir, "fimo", subdir, "fimo.txt"),
                                                    "fake_long_fimo.txt")

            self.assertEqual(fimo_result, expected_fimo_result)
