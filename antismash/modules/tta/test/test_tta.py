# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

from argparse import Namespace
import unittest

from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation

from antismash.common.secmet.record import Record
from antismash.common.test.helpers import DummyCDS, DummyProtocluster, get_simple_options
from antismash.modules import tta


class TtaTest(unittest.TestCase):
    def setUp(self):
        # locations:            VVV         VVV
        record = Record(Seq("ATGTTATGAGGGTCATAACAT"))

        record.add_cds_feature(DummyCDS(0, 9, strand=1))
        record.add_cds_feature(DummyCDS(12, 21, strand=-1))

        cluster = DummyProtocluster(start=0, end=21)
        record.add_protocluster(cluster)
        record.create_candidate_clusters()
        record.create_regions()
        # if these aren't correct, the tests will fail
        assert len(cluster.cds_children) == 2
        assert len(record.get_regions()) == 1
        for cds in record.get_cds_features():
            assert cds.is_contained_by(cluster)
            assert cds.extract(record.seq) == "ATGTTATGA", str(cds.location)

        self.record = record

    def test_check_prereqs(self):
        """Test tta.check_prereqs()"""
        assert not tta.check_prereqs(None)

    def test_options(self):
        options = Namespace()
        for threshold in [0., 0.5, 1.0]:
            options.tta_threshold = threshold
            assert not tta.check_options(options)

        options.tta_threshold = 1.01
        assert tta.check_options(options) == ["Supplied threshold is out of range 0 to 1: 1.01"]
        options.tta_threshold = -0.01
        assert tta.check_options(options) == ["Supplied threshold is out of range 0 to 1: -0.01"]

    def test_detect(self):
        """Test tta.detect()"""
        self.assertEqual(len(self.record.get_cds_features()) + len(self.record.get_protoclusters()), 3)
        options = get_simple_options(tta, ["--tta-threshold", "0"])

        detected = tta.detect(self.record, options)
        self.assertEqual(len(detected), 2)
        # make sure features not added yet
        features = self.record.get_generics()
        self.assertEqual(len(features), 0)
        # add to record and make sure they're there
        detected.add_to_record(self.record)
        features = self.record.get_generics()
        self.assertEqual(len(features), 2)

        for feature in features:
            assert feature.type == "misc_feature"

        fw_tta, rv_tta = features[0], features[1]
        self.assertEqual(fw_tta.location.start, 3)
        self.assertEqual(fw_tta.location.end, 6)
        self.assertEqual(fw_tta.strand, 1)

        self.assertEqual(rv_tta.location.start, 15)
        self.assertEqual(rv_tta.location.end, 18)
        self.assertEqual(rv_tta.strand, -1)

    def test_feature_creation(self):
        fw_loc = FeatureLocation(210, 300, strand=1)
        fw_feature = SeqFeature(fw_loc, type='CDS')
        results = tta.tta.TTAResults('dummy', gc_content=1, threshold=0.65)
        ret = results.new_feature_from_other(fw_feature, 12)
        self.assertEqual(ret.strand, 1)
        self.assertEqual(ret.location.start, 222)
        self.assertEqual(ret.location.end, 225)

        rv_loc = FeatureLocation(210, 300, strand=-1)
        rv_feature = SeqFeature(rv_loc, type='CDS')
        ret = results.new_feature_from_other(rv_feature, 12)
        self.assertEqual(ret.strand, -1)
        self.assertEqual(ret.location.start, 285)
        self.assertEqual(ret.location.end, 288)
