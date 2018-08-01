# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import unittest

from antismash.common import path
from antismash.detection import hmm_detection


class TestClusterCSS(unittest.TestCase):
    def test_css_matches_rules(self):
        defined_clusters = set(name.lower() for name in hmm_detection.get_supported_cluster_types())
        available_classes = set()
        base_classes = {"hybrid"}  # a special case used at the javascript level
        less = path.get_full_path(__file__, "..", "css", "secmet.less")
        with open(less) as handle:
            for line in handle.readlines():
                if line.startswith('.'):
                    class_ = line[1:].split()[0]
                    available_classes.add(class_.lower())
                elif line.strip().startswith('.'):  # a reference/base class
                    class_ = line[1:].split()[0].split('(', 1)[0]
                    # now it's .thing, with a possible trailing ;
                    class_ = class_[1:].strip(";")
                    base_classes.add(class_.lower())
        missing_css = defined_clusters - available_classes
        assert not missing_css
        # allow for the extra base classes and hybrids
        extra_css = available_classes - defined_clusters - base_classes
        # and clusterfinders clustertypes
        extra_css -= {'cf_putative', 'cf_fatty_acid', 'cf_saccharide'}
        assert not extra_css
