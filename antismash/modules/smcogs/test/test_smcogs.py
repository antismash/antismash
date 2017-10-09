# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import unittest

from antismash.common.secmet import GeneFunction
from antismash.modules.smcogs import classify


class TestSMCOGLoad(unittest.TestCase):
    def test_load(self):
        # this mostly just tests that the cog annotation file isn't corrupted
        annotations = classify.load_cog_annotations()
        assert len(annotations) == 301
        for key, function in annotations.items():
            assert isinstance(function, GeneFunction), "cog annotation %s has bad type" % key
