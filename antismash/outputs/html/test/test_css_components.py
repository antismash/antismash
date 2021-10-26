# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import unittest

from antismash.common import path
from antismash.detection import hmm_detection


class TestClusterCSS(unittest.TestCase):
    def test_css_matches_rules(self):
        rules = hmm_detection._get_rules("loose")
        available_classes = set()
        base_classes = {
            "hybrid",  # a special case used at the javascript level
            "unknown",  # for regions containing only subregions
        }
        less = path.get_full_path(__file__, "..", "css", "secmet.scss")
        with open(less) as handle:
            for line in handle.readlines():
                if line.startswith('.'):
                    class_ = line[1:].split()[0]
                    available_classes.add(class_)
        missing_css = [f"{rule.name} (category: {rule.category})"
                       for rule in rules
                       if not available_classes.intersection({rule.name, rule.category})]
        assert not missing_css
        # allow for the extra base classes and hybrids
        products = {rule.name for rule in rules}
        categories = {rule.category for rule in rules}
        extra_css = available_classes.difference(products, categories, base_classes)
        assert not extra_css
