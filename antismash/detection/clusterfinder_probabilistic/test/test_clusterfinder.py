# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import json
import unittest

from Bio.Seq import Seq

from antismash.common.secmet.features import CDSFeature, PFAMDomain, FeatureLocation
from antismash.common.test.helpers import DummyRecord
from antismash.config import build_config, update_config, destroy_config
from antismash.detection import clusterfinder_probabilistic as clusterfinder


class ClusterFinderTest(unittest.TestCase):
    def setUp(self):
        self.config = build_config(["--cf-create-clusters",
                                    "--cf-mean-threshold", "0.6",
                                    "--cf-min-cds", "5",
                                    "--cf-min-pfams", "5"], modules=[clusterfinder],
                                   isolated=True)
        update_config({"enabled_cluster_types": []})

        self.record = DummyRecord(seq=Seq("A" * 2000))
        for start, end, probability, pfam_id in [(10, 20, 0.1, 'PF77777'),
                                                 (30, 40, 0.3, 'PF00106'),
                                                 (50, 60, 0.4, 'PF00107'),
                                                 (60, 70, 0.7, 'PF00109'),
                                                 (70, 80, 0.98, 'PF08484'),
                                                 (90, 100, 0.8, 'PF02401'),
                                                 (100, 110, 0.32, 'PF04369'),
                                                 (110, 120, 1.0, 'PF00128'),
                                                 (130, 140, 0.2, 'PF77776'),
                                                 (500, 505, None, 'PF77775'),
                                                 (1010, 1020, 0.1, 'PF77774'),
                                                 (1030, 1040, 0.3, 'PF00106'),
                                                 (1050, 1060, 0.4, 'PF00107'),
                                                 (1060, 1070, 0.7, 'PF00109'),
                                                 (1070, 1080, 0.98, 'PF08484'),
                                                 (1090, 1100, 0.8, 'PF02401'),
                                                 (1100, 1110, 0.32, 'PF04369'),
                                                 (1110, 1120, 1.0, 'PF00128')]:
            location = FeatureLocation(start, end, strand=1)
            self.record.add_cds_feature(CDSFeature(location, locus_tag=str(start), translation="A"))
            pfam = PFAMDomain(location, "dummy_description", protein_start=start + 1,
                              protein_end=end-1, identifier=pfam_id, tool="test")
            pfam.domain_id = "pfam_%d" % start
            pfam.probability = probability
            self.record.add_pfam_domain(pfam)

    def tearDown(self):
        destroy_config()

    def test_options(self):
        for val in [-1., -.01, -0.0001, 1.0001, 2.]:
            update_config({"cf_threshold": val})
            assert len(clusterfinder.check_options(self.config)) == 1

    def test_check_prereqs(self):
        self.assertEqual([], clusterfinder.check_prereqs())

    def test_find_nr_cds(self):
        left = (0, 5)
        newpos, num = clusterfinder.probabilistic.find_nr_cds(left, self.record)
        self.assertEqual(left, newpos)
        self.assertEqual(0, num)

        right = (150, 160)
        newpos, count = clusterfinder.probabilistic.find_nr_cds(right, self.record)
        assert newpos == right
        assert count == 0

        middle = (35, 115)
        newpos, count = clusterfinder.probabilistic.find_nr_cds(middle, self.record)
        assert newpos == (30, 120)
        assert count == 7

        small = (501, 504)
        newpos, count = clusterfinder.probabilistic.find_nr_cds(small, self.record)
        assert newpos == (500, 505)
        assert count == 1

    def test_find_probabilistic_clusters(self):
        ret = clusterfinder.find_probabilistic_clusters(self.record, self.config)
        assert len(ret) == 2
        assert ret[0].location.start == 30
        assert ret[0].location.end == 120
        assert ret[1].location.start == 1030
        assert ret[1].location.end == 1120

    def test_no_overlaps(self):
        results = clusterfinder.generate_results(self.record, self.config)
        areas = self.record.get_subregions()

        assert len(results.areas) == 2
        assert len(areas) == 2
        assert list(areas) == results.areas

        area = areas[0]
        assert area.location.start == 30
        assert area.location.end == 120
        self.assertAlmostEqual(0.6429, area.probability, places=4)

        area = areas[1]
        assert area.location.start == 1030
        assert area.location.end == 1120
        self.assertAlmostEqual(0.6429, area.probability, places=4)

    def test_full_loop(self):
        results = clusterfinder.run_on_record(self.record, None, self.config)
        old_areas = len(results.areas)

        json_string = json.dumps(results.to_json())
        regenned = clusterfinder.regenerate_previous_results(json.loads(json_string), self.record, self.config)
        assert len(regenned.areas) == old_areas
        for old, new in zip(results.areas, regenned.areas):
            assert old.location == new.location
            self.assertAlmostEqual(old.probability, new.probability, places=6)


class TestResults(unittest.TestCase):
    def setUp(self):
        self.config = build_config(["--cf-create-clusters"], modules=[clusterfinder], isolated=True)
        self.record = DummyRecord()

    def regen(self, results):
        return clusterfinder.regenerate_previous_results(results, self.record, self.config)

    def test_empty(self):
        assert self.regen(None) is None
        assert self.regen({}) is None

    def test_schema_mismatch(self):
        assert self.regen({"schema": 0}) is None

    def test_creation_change(self):
        results = clusterfinder.ClusterFinderResults(self.record.id, [], create=False)
        with self.assertRaisesRegex(ValueError, "detection results have changed"):
            self.regen(json.loads(json.dumps(results.to_json())))
