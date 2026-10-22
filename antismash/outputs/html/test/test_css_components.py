# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import re
import unittest
from unittest.mock import Mock

from antismash.common import path
from antismash.detection import hmm_detection
from antismash.outputs.html import js


class TestClusterCSS(unittest.TestCase):
    def test_regions_have_shared_class(self):
        region = Mock(product_categories={"PKS"})
        region.get_unique_protoclusters.return_value = [Mock(product="T1PKS")]
        assert js.get_region_css(region) == "PKS T1PKS secmet"

        region.product_categories = {"PKS", "NRPS"}
        assert js.get_region_css(region) == "hybrid secmet"

        region.product_categories = {"unknown"}
        region.get_unique_protoclusters.return_value = []
        assert js.get_region_css(region) == "unknown secmet"

    def test_css_matches_rules(self):
        rules = hmm_detection._get_rules("loose")
        available_classes = set()
        base_classes = {
            "hybrid",  # a special case used at the javascript level
            "unknown",  # for regions containing only subregions
        }
        css_path = path.get_full_path(__file__, "..", "css", "secmet.css")
        with open(css_path, encoding="utf-8") as handle:
            contents = handle.read()
        expected_properties = {
            "--secmet-fg", "--secmet-bg", "--secmet-highlight", "--secmet-hover-bg",
        }
        for selectors, declarations in re.findall(r"([^{}]+)\{([^{}]*)\}", contents):
            properties = set(re.findall(r"(--secmet-[a-z-]+)\s*:", declarations))
            if properties:
                assert properties == expected_properties, selectors.strip()
                available_classes.update(re.findall(r"\.([A-Za-z0-9_-]+)", selectors))

        assert re.search(
            r"\.secmet\s*\{\s*background-color:\s*var\(--secmet-bg\)", contents,
        )

        missing_css = [f"{rule.name} (category: {rule.category})"
                       for rule in rules
                       if not available_classes.intersection({rule.name, rule.category})]
        assert not missing_css
        # allow for the extra base classes and hybrids
        products = {rule.name for rule in rules}
        categories = {rule.category for rule in rules}
        extra_css = available_classes.difference(products, categories, base_classes)
        assert not extra_css
