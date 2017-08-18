import unittest

from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
from antismash.common import deprecated
from antismash.common.secmet.record import Record
from antismash.common.secmet.feature import CDSFeature, Cluster
from antismash.common.test.helpers import DummyCDS
from antismash.config import args
from antismash.modules import tta

class TtaTest(unittest.TestCase):
    def setUp(self):
        # locations:            VVV         VVV
        record = Record(Seq("ATGTTATGAGGGTCATAACAT", generic_dna))

        record.add_cds_feature(DummyCDS(0, 9, strand=1))
        record.add_cds_feature(DummyCDS(12, 21, strand=-1))

        cluster_loc = FeatureLocation(0, 21)
        cluster = Cluster(cluster_loc, 0, 0, [])
        record.add_cluster(cluster)
        # if these aren't correct, the tests will fail
        assert len(cluster.cds_children) == 2
        for cds in record.get_cds_features():
            assert cds.overlaps_with(cluster)
            assert cds.cluster == cluster, str(cds.location)
            assert cds.extract(record.seq) == "ATGTTATGA", str(cds.location)

        self.record = record


    def test_check_prereqs(self):
        """Test tta.check_prereqs()"""
        assert not tta.check_prereqs()


    def test_detect(self):
        """Test tta.detect()"""
        self.assertEqual(len(self.record.get_cds_features()) + len(self.record.get_clusters()), 3)
        options = args.simple_options(tta, ["--tta"])
        detected = tta.detect(self.record, options)
        self.assertEqual(len(detected), 2)
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
        results = tta.tta.TTAResults('dummy')
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
