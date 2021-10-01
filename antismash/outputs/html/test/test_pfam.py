# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import unittest

from antismash.outputs.html.visualisers.generic_domains import (
    PFAM_BIOSYNTHETIC as BIOSYNTHETIC,
    PFAM_REGULATORY as REGULATORY,
    PFAM_TRANSPORT as TRANSPORT,
)

class TestClassifications(unittest.TestCase):
    def test_biosynth_regulatory_overlap(self):
        assert not BIOSYNTHETIC.intersection(REGULATORY)

    def test_biosynth_transport_overlap(self):
        assert not BIOSYNTHETIC.intersection(TRANSPORT)

    def test_regulatory_transport_overlap(self):
        assert not REGULATORY.intersection(TRANSPORT)
