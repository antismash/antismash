# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""Test suite for the cassis cluster detection plugin"""

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring,consider-using-with

from argparse import Namespace
import os
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation

from antismash.common import path, secmet
from antismash.common.test import helpers
from antismash.config import build_config, destroy_config
from antismash.detection import cassis
from antismash.detection.cassis.cluster_prediction import ClusterMarker
from antismash.detection.cassis.motifs import Motif
from antismash.detection.cassis.promoters import Promoter, CombinedPromoter, get_promoters

cassis.VERBOSE_DEBUG = True


def convert_newline(string):
    """Convert all line endings to \n for OS independency"""
    return "\n".join(string.splitlines())


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
        if i in (3, 5):
            cds.gene_functions.add(secmet.qualifiers.GeneFunction.CORE, "testtool", "dummy", "product")

    return seq_record


class CassisTestCore(unittest.TestCase):
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


class TestCassisMethods(CassisTestCore):
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

        promoters = get_promoters(seq_record, genes, upstream_tss, downstream_tss)
        self.assertEqual(list(map(lambda x: [x.start, x.end], promoters)), expected_promoters)
        cassis.write_promoters_to_file(self.options.output_dir, seq_record.name, promoters)
        # read expected files and save to string variable
        expected_sequences_file = ""
        with open(path.get_full_path(__file__, "data", "expected_promoter_sequences.fasta"),
                  encoding="utf-8") as handle:
            expected_sequences_file = handle.read()
        expected_sequences_file = convert_newline(expected_sequences_file.rstrip())

        expected_positions_file = ""
        with open(path.get_full_path(__file__, "data", "expected_promoter_positions.csv"),
                  encoding="utf-8") as handle:
            expected_positions_file = handle.read()
        expected_positions_file = convert_newline(expected_positions_file.rstrip())

        # read test files and save to string variable
        sequences_file = ""
        with open(os.path.join(self.options.output_dir, seq_record.name + "_promoter_sequences.fasta"),
                  encoding="utf-8") as handle:
            sequences_file = handle.read()
        sequences_file = convert_newline(sequences_file.rstrip())

        positions_file = ""
        with open(os.path.join(self.options.output_dir, seq_record.name + "_promoter_positions.csv"),
                  encoding="utf-8") as handle:
            positions_file = handle.read()
        positions_file = convert_newline(positions_file.rstrip())

        self.assertEqual(sequences_file, expected_sequences_file)
        self.assertEqual(positions_file, expected_positions_file)

    def test_cleanup_outdir(self):
        anchor_genes = ["gene1", "gene4"]
        cluster = cassis.ClusterPrediction(ClusterMarker("gene1", Motif(3, 3, score=1)),
                                           ClusterMarker("gene4", Motif(3, 3, score=1)))
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
        motif = Motif(3, 3)
        assert motif.pairing_string == "+03_-03"
        motif.plus = 4
        assert motif.pairing_string == "+04_-03"
        motif.minus = 2
        assert motif.pairing_string == "+04_-02"


class TestPromoters(unittest.TestCase):
    def test_promoter_id(self):
        assert Promoter("gene1", 1, 5).get_id() == "gene1"
        assert CombinedPromoter("gene1", "gene2", 1, 5).get_id() == "gene1+gene2"


class TestCassisStorageMethods(unittest.TestCase):
    def test_store_promoters(self):
        promoters = [Promoter("gene1", 10, 20, seq=Seq("cgtacgtacgt")),
                     Promoter("gene2", 30, 40, seq=Seq("cgtacgtacgt")),
                     CombinedPromoter("gene3", "gene4", 50, 60, seq=Seq("cgtacgtacgt"))]
        record_with_promoters = create_fake_record()
        cassis.store_promoters(promoters, record_with_promoters)  # add ("store") promoters to seq_record

        record_without_promoters = create_fake_record()  # just the same, without adding promoters

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

    def test_store_subregions(self):
        # this test is similar to test_store_promoters
        anchor = "gene3"

        start_marker = ClusterMarker("gene1", Motif(3, 3, score=1))
        start_marker.promoter = "gene1"
        start_marker.abundance = 2
        end_marker = ClusterMarker("gene4", Motif(3, 3, score=1))
        end_marker.promoter = "gene3+gene4"
        assert end_marker.abundance == 1
        first_cluster = cassis.ClusterPrediction(start_marker, end_marker)
        first_cluster.promoters = 3
        first_cluster.genes = 4

        start_marker = ClusterMarker("gene1", Motif(4, 4, score=1))
        start_marker.promoter = "gene1"
        assert start_marker.abundance == 1
        end_marker = ClusterMarker("gene5", Motif(4, 4, score=1))
        end_marker.promoter = "gene5"
        assert end_marker.abundance == 1
        second_cluster = cassis.ClusterPrediction(start_marker, end_marker)
        second_cluster.promoters = 3
        second_cluster.genes = 4

        # order reversed because subregions are ordered by length when starts are the same
        region_predictions = [second_cluster, first_cluster]

        record_with_subregions = create_fake_record()
        record_without_subregions = create_fake_record()  # just the same, without adding subregions

        subregions = cassis.create_subregions(anchor, region_predictions, record_with_subregions)
        assert record_with_subregions.get_feature_count() == record_without_subregions.get_feature_count()

        for region in subregions:
            record_with_subregions.add_subregion(region)

        # test subregion features
        expected_count = record_without_subregions.get_feature_count() + len(subregions)
        assert record_with_subregions.get_feature_count() == expected_count
        for i, region in enumerate(region_predictions):
            subregion = record_with_subregions.get_subregions()[i]
            self.assertEqual(subregion.type, "subregion")
            self.assertEqual(subregion.tool, "cassis")
            self.assertEqual(subregion.label, anchor)
            self.assertEqual(subregion.get_qualifier("genes"), (region.genes,))
            self.assertEqual(subregion.get_qualifier("promoters"), (region.promoters,))
            self.assertEqual(subregion.get_qualifier("gene_left"), (region.start.gene,))
            self.assertEqual(subregion.get_qualifier("gene_right"), (region.end.gene,))
            # don't test all feature qualifiers, only some


class TestSubRegionFiltering(unittest.TestCase):
    def create_sub(self, start, end, anchor):
        return secmet.SubRegion(FeatureLocation(start, end), label=anchor, tool="cassis")

    def test_empty(self):
        assert cassis.filter_subregions([]) == []

    def test_unrelated_overlap(self):
        subs = [self.create_sub(300, 500, "A"),
                self.create_sub(400, 600, "B")]
        assert cassis.filter_subregions(subs) == subs

    def test_unrelated_contains(self):
        subs = [self.create_sub(300, 500, "A"),
                self.create_sub(400, 450, "B")]
        assert cassis.filter_subregions(subs) == subs

    def test_related_overlap(self):
        subs = [self.create_sub(300, 500, "A"),
                self.create_sub(400, 600, "A")]
        assert cassis.filter_subregions(subs) == subs

    def test_related_contains(self):
        subs = [self.create_sub(300, 600, "A"),
                self.create_sub(400, 500, "A")]
        assert cassis.filter_subregions(subs) == [subs[0]]

    def test_multiple(self):
        subs = [self.create_sub(300, 600, "A"),
                self.create_sub(400, 500, "A"),
                self.create_sub(400, 700, "A"),
                self.create_sub(400, 500, "B"),
                self.create_sub(300, 600, "B")]

        res = cassis.filter_subregions(subs)
        assert res == [subs[0], subs[4], subs[2]]


class TestResults(unittest.TestCase):
    def setUp(self):
        self.old_max_perc = cassis.MAX_PERCENTAGE
        self.old_max_gap = cassis.MAX_GAP_LENGTH

    def tearDown(self):
        cassis.MAX_PERCENTAGE = self.old_max_perc
        cassis.MAX_GAP_LENGTH = self.old_max_gap

    def test_base(self):
        results = cassis.CassisResults("test")
        assert results.record_id == "test"
        assert results.subregions == []
        assert results.promoters == []

    def test_regeneration(self):
        record = create_fake_record()
        results = cassis.CassisResults(record.id)
        # create a prediction, since it will generate a border with many extra qualifiers
        start_marker = ClusterMarker("gene1", Motif(3, 3, score=1))
        start_marker.promoter = "gene1"
        start_marker.abundance = 2
        end_marker = ClusterMarker("gene4", Motif(3, 3, score=1))
        end_marker.promoter = "gene3+gene4"
        assert end_marker.abundance == 1
        cluster = cassis.ClusterPrediction(start_marker, end_marker)
        results.subregions = cassis.create_subregions("gene1", [cluster], record)
        assert results.subregions

        results.promoters = [Promoter("gene1", 10, 20, seq=Seq("cgtacgtacgt")),
                             Promoter("gene2", 30, 40, seq=Seq("cgtacgtacgt")),
                             CombinedPromoter("gene3", "gene4", 50, 60, seq=Seq("cgtacgtacgt"))]

        round_trip = cassis.regenerate_previous_results(results.to_json(), record, None)
        assert isinstance(round_trip, cassis.CassisResults)
        assert len(results.subregions) == len(round_trip.subregions)
        for old, new in zip(results.subregions, round_trip.subregions):
            assert old.location == new.location
            assert old.to_biopython()[0].qualifiers == new.to_biopython()[0].qualifiers
        assert round_trip.promoters == results.promoters

    def test_changed_max_percentage(self):
        record = create_fake_record()
        json = cassis.CassisResults(record.id).to_json()
        assert isinstance(cassis.regenerate_previous_results(json, record, None),
                          cassis.CassisResults)
        cassis.MAX_PERCENTAGE += 5
        assert cassis.regenerate_previous_results(json, record, None) is None

    def test_changed_max_gap_length(self):
        record = create_fake_record()
        json = cassis.CassisResults(record.id).to_json()
        assert isinstance(cassis.regenerate_previous_results(json, record, None),
                          cassis.CassisResults)
        cassis.MAX_GAP_LENGTH += 1
        assert cassis.regenerate_previous_results(json, record, None) is None

    def test_not_same_record(self):
        record = create_fake_record()
        record.id = "A"
        other = create_fake_record()
        other.id = "B"
        json = cassis.CassisResults(record.id).to_json()
        assert isinstance(cassis.regenerate_previous_results(json, record, None),
                          cassis.CassisResults)
        assert cassis.regenerate_previous_results(json, other, None) is None

    def test_run_on_record_skips_work(self):
        record = create_fake_record()
        record.id = "real"
        results = cassis.CassisResults(record.id)

        with patch.object(cassis, "detect", return_value=cassis.CassisResults("fake")):
            assert cassis.run_on_record(record, results, None).record_id == "real"
            assert cassis.run_on_record(record, None, None).record_id == "fake"


class TestVersioning(unittest.TestCase):
    def setUp(self):
        self.config = Namespace()
        self.config.executables = Namespace()

    def check_with_version(self, fimo_version, meme_version):
        self.config.executables.meme = "meme"
        self.config.executables.fimo = "fimo"
        with patch.object(cassis.subprocessing, "run_meme_version", return_value=meme_version):
            with patch.object(cassis.subprocessing, "run_fimo_version", return_value=fimo_version):
                return cassis.check_prereqs(self.config)

    def test_missing(self):
        messages = cassis.check_prereqs(self.config)
        assert len(messages) == 2
        for message in messages:
            assert "Failed to locate executable" in message

    def test_correct_version(self):
        messages = self.check_with_version(fimo_version="4.11.2", meme_version="4.11.2")
        messages = self.check_with_version(fimo_version="5.5.2", meme_version="4.11.2")
        assert not messages

    def test_incorrect_version(self):
        fimo = "4.11.0"
        meme = "5.5.2"
        messages = self.check_with_version(fimo_version=fimo, meme_version=meme)
        assert messages == [
            "Incompatible MEME version, expected 4.11.2 but found 5.5.2",
            "Incompatible FIMO version, expected 4.11.2 or later but found 4.11.0",
        ]
        assert len(messages) == 2
        assert f"expected 4.11.2 but found {meme}" in messages[0]
        assert f"expected 4.11.2 or later but found {fimo}" in messages[1]
