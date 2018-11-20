# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import unittest

from antismash.common import utils
from antismash.common.secmet import FeatureLocation
from antismash.common.test.helpers import DummyCDS, DummyRecord


class TestRobustProteinAnalysis(unittest.TestCase):
    def test_init(self):
        """Test RobustProteinAnalysis initialisation"""
        rpa = utils.RobustProteinAnalysis("MAGICHAT", ignore_invalid=True)
        self.assertIsInstance(rpa, utils.RobustProteinAnalysis)

        rpa = utils.RobustProteinAnalysis("MAGICHAT", ignore_invalid=False)
        self.assertIsInstance(rpa, utils.RobustProteinAnalysis)

        for bad_invalid in ["none", None, 3, []]:
            with self.assertRaises(TypeError):
                utils.RobustProteinAnalysis("MAGICHAT", ignore_invalid=bad_invalid)

    def test_uppercase(self):
        """Test RobustProteinAnalysis converts passed sequence to upper case"""
        rpa = utils.RobustProteinAnalysis("Magichat")
        assert rpa.original_sequence == "MAGICHAT"
        assert rpa.sequence == "MAGICHAT"

    def test_molecular_weight_ignore(self):
        """ Test RobustProteinAnalysis.molecular_weight() calculates
            correct weight when ignoring invalids"""
        rpa = utils.RobustProteinAnalysis("MAGICXHAT")
        self.assertEqual(802.9621, rpa.molecular_weight())  # default is True

        rpa = utils.RobustProteinAnalysis("MAGICXHAT", ignore_invalid=True)
        self.assertEqual(802.9621, rpa.molecular_weight())

    def test_molecular_weight_average(self):
        """ Test RobustProteinAnalysis.molecular_weight() calculates
            correct weight when not ignoring invalids
        """
        rpa = utils.RobustProteinAnalysis("MAGICXHAT", ignore_invalid=False)
        self.assertEqual(912.9621, rpa.molecular_weight())


class TestSignatureBuilding(unittest.TestCase):
    def test_extract_by_reference_positions(self):
        sig = utils.extract_by_reference_positions("ABC-DE-F", "A-BC-DEF", [0, 1, 3, 4])
        assert sig == "ACE-"
        sig = utils.extract_by_reference_positions("ABCDF", "ABCDE", [0, 1, 3, 4])
        assert sig == "ABDF"


class TestDistanceCalculations(unittest.TestCase):
    def setUp(self):
        self.record = DummyRecord(seq="A" * 100000)
        self.query = self.create_cds(50000, 50000, ["query_gene_prof"])
        self.record.add_cds_feature(self.query)

    def create_cds(self, start, end, profiles, strand=1):
        cds = DummyCDS(start, end, strand=strand, translation="A")
        for profile in profiles:
            cds.sec_met.add_domains([cds.sec_met.Domain(profile, 1e-5, 20.5, 2, "test")])
        return cds

    def test_empty_record(self):
        self.record._cds_features.clear()
        assert utils.distance_to_pfam(self.record, self.query, []) == -1

    def test_self_hit(self):
        assert utils.distance_to_pfam(self.record, self.query, ["query_gene_prof"]) == 0

    def test_simple_before(self):
        cds = self.create_cds(29000, 30000, profiles=["left20k"])
        self.record.add_cds_feature(cds)
        assert utils.distance_to_pfam(self.record, self.query, ["notleft20k"]) == -1
        assert utils.distance_to_pfam(self.record, self.query, ["left20k"]) == 20000

    def test_simple_after(self):
        cds = self.create_cds(60000, 63000, profiles=["right10k"])
        self.record.add_cds_feature(cds)
        assert utils.distance_to_pfam(self.record, self.query, ["notright10k"]) == -1
        assert utils.distance_to_pfam(self.record, self.query, ["right10k"]) == 10000

    def test_outside_before(self):
        cds = self.create_cds(5000, 9999, profiles=["outside"])
        self.record.add_cds_feature(cds)
        assert utils.distance_to_pfam(self.record, self.query, ["outside"]) == -1

    def test_outside_after(self):
        cds = self.create_cds(90001, 95000, profiles=["outside"])
        self.record.add_cds_feature(cds)
        assert utils.distance_to_pfam(self.record, self.query, ["outside"]) == -1

    def test_edge_overlap_after(self):
        cds = self.create_cds(90000, 91000, profiles=["r.edge"])
        self.record.add_cds_feature(cds)
        assert utils.distance_to_pfam(self.record, self.query, ["r.edge"]) == -1

        cds.location = FeatureLocation(89999, 91000, strand=1)
        assert utils.distance_to_pfam(self.record, self.query, ["r.edge"]) == 39999

        cds.location = FeatureLocation(89999, 91000, strand=-1)
        assert utils.distance_to_pfam(self.record, self.query, ["r.edge"]) == 39999

    def test_edge_overlap_before(self):
        cds = self.create_cds(9000, 10000, profiles=["l.edge"])
        self.record.add_cds_feature(cds)
        assert utils.distance_to_pfam(self.record, self.query, ["l.edge"]) == -1

        cds.location = FeatureLocation(9000, 10001, strand=1)
        assert utils.distance_to_pfam(self.record, self.query, ["l.edge"]) == 39999

        cds.location = FeatureLocation(9000, 10001, strand=-1)
        assert utils.distance_to_pfam(self.record, self.query, ["l.edge"]) == 39999

    def test_with_no_secmet(self):
        cds = self.create_cds(55000, 60000, profiles=[])
        cds.sec_met = None
        self.record.add_cds_feature(cds)
        assert utils.distance_to_pfam(self.record, self.query, ["test"]) == -1
