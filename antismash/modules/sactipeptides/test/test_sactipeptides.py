# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import unittest

from antismash.common.test.helpers import DummyCDS, DummyRecord
from antismash.modules.sactipeptides.specific_analysis import (
    SactiResults,
)


class TestCDSDuplication(unittest.TestCase):
    def setUp(self):
        self.record = DummyRecord(seq="A" * 100)
        self.cds = DummyCDS()

    def test_existing(self):
        self.record.add_cds_feature(self.cds)
        results = SactiResults(self.record.id)
        results.add_cds(self.cds)
        assert len(self.record.get_cds_features()) == 1
        results.add_to_record(self.record)
        assert len(self.record.get_cds_features()) == 1

    def test_new(self):
        results = SactiResults(self.record.id)
        results.add_cds(self.cds)
        assert len(self.record.get_cds_features()) == 0
        results.add_to_record(self.record)
        assert len(self.record.get_cds_features()) == 1

    def test_double_protocluster(self):
        results = SactiResults(self.record.id)
        assert len(results._new_cds_features) == 0
        results.add_cds(self.cds)
        assert len(results._new_cds_features) == 1
        results.add_cds(self.cds)
        assert len(results._new_cds_features) == 1
        results.add_cds(DummyCDS(locus_tag="different"))
        assert len(results._new_cds_features) == 2
