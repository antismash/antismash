# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring,consider-using-with

from io import StringIO
import unittest

from antismash.common.hmm_rule_parser.categories import (
    _PARSED_FILES,
    _parse_categories,
    parse_category_file,
)
from antismash.common import json


class TestCategoryParsing(unittest.TestCase):
    def setUp(self):
        self.data = {
            "A": {"description": "some desc", "version": 1},
            "B": {"description": "some other desc", "version": 1},
        }
        _PARSED_FILES.clear()

    def tearDown(self):
        _PARSED_FILES.clear()

    def test_no_parents(self):
        cats = {cat.name: cat for cat in _parse_categories(self.data)}
        for name, value in self.data.items():
            assert cats[name].description == value["description"]
            assert cats[name].version == value["version"]
            assert cats[name].name == name
            assert cats[name].parent is None
        assert len(cats) == len(self.data)

    def test_missing_parents(self):
        parent = "Z"
        assert parent not in self.data  # should also be missing
        self.data["A"]["parent"] = parent
        with self.assertRaisesRegex(ValueError, "refers to unknown category"):
            _parse_categories(self.data)

    def test_with_parents(self):
        self.data["A"]["parent"] = "B"
        cats = {cat.name: cat for cat in _parse_categories(self.data)}
        for name, value in self.data.items():
            assert cats[name].description == value["description"]
            assert cats[name].version == value["version"]
        assert cats["A"].parent == cats["B"]
        assert cats["B"].parent is None

    def test_caching(self):
        def make_fake_handle():
            handle = StringIO(json.dumps(self.data))
            handle.name = "test"  # this exists on real files
            return handle

        assert not _PARSED_FILES
        parsed = parse_category_file(make_fake_handle())
        assert len(parsed) == 2
        assert len(_PARSED_FILES) == 1
        reparsed = parse_category_file(make_fake_handle())
        assert len(_PARSED_FILES) == 1
        assert reparsed is parsed
        assert len(parsed) == 2  # check it doesn't duplicate
