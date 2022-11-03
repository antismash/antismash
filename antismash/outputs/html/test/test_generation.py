# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import unittest

from antismash.main import get_all_modules
from antismash.common.secmet.test.helpers import DummySubRegion
from antismash.common.test.helpers import DummyRecord
from antismash.config import build_config, destroy_config, update_config
from antismash.outputs.html.generator import generate_webpage


class TestOutput(unittest.TestCase):
    def setUp(self):
        self.options = build_config([], modules=get_all_modules())

    def tearDown(self):
        destroy_config()

    def test_subregions_only(self):
        # construct a region with only subregions
        record = DummyRecord(seq="A" * 10)
        record.add_subregion(DummySubRegion())
        record.create_regions()

        # add information normally added by record processing
        record.record_index = 1
        update_config({"triggered_limit": False})

        # make sure HTML is generated without errors
        assert generate_webpage([record], [{}], self.options, get_all_modules())
