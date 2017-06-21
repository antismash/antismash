import unittest

from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
from antismash.common import deprecated
from antismash.config import args
from antismash.modules import tta

class TtaTest(unittest.TestCase):
    def setUp(self):
        sequence = Seq("ATGTTATGAGGGTCATAACAT", generic_dna)
        record = SeqRecord(seq=sequence, id="test", name="test")

        fw_loc = FeatureLocation(0, 9, strand=1)
        fw_feature = SeqFeature(fw_loc, type="CDS")
        record.features.append(fw_feature)

        rv_loc = FeatureLocation(12, 21, strand=-1)
        rv_feature = SeqFeature(rv_loc, type="CDS")
        record.features.append(rv_feature)

        cluster_loc = FeatureLocation(0, 21)
        cluster_feature = SeqFeature(cluster_loc, type="cluster")
        record.features.append(cluster_feature)
        self.record = record


    def test_check_prereqs(self):
        """Test tta.check_prereqs()"""
        assert not tta.check_prereqs()


    def test_detect(self):
        """Test tta.detect()"""
        self.assertEqual(len(self.record.features), 3)
        options = args.simple_options(tta, ["--tta"])
        detected = tta.detect(self.record, options)
        self.assertEqual(len(detected), 2)
        self.assertEqual(len(self.record.features), 5)

        features = deprecated.get_all_features_of_type(self.record, 'misc_feature')
        self.assertEqual(len(features), 2)

        fw_tta, rv_tta = features[0], features[1]
        self.assertEqual(fw_tta.location.start, 3)
        self.assertEqual(fw_tta.location.end, 6)
        self.assertEqual(fw_tta.strand, 1)

        self.assertEqual(rv_tta.location.start, 15)
        self.assertEqual(rv_tta.location.end, 18)
        self.assertEqual(rv_tta.strand, -1)

    def test__create_tta_feature(self):
        """Test tta._create_tta_feature()"""
        fw_loc = FeatureLocation(210, 300, strand=1)
        fw_feature = SeqFeature(fw_loc, type='CDS')
        ret = tta.tta._create_tta_feature(fw_feature, 12)
        self.assertEqual(ret.strand, 1)
        self.assertEqual(ret.location.start, 222)
        self.assertEqual(ret.location.end, 225)

        rv_loc = FeatureLocation(210, 300, strand=-1)
        rv_feature = SeqFeature(rv_loc, type='CDS')
        ret = tta.tta._create_tta_feature(rv_feature, 12)
        self.assertEqual(ret.strand, -1)
        self.assertEqual(ret.location.start, 285)
        self.assertEqual(ret.location.end, 288)
