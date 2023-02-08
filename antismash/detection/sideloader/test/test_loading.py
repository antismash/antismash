# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring,too-many-public-methods

import json
import unittest

from antismash.common import errors, path
from antismash.detection.sideloader import general, loader

GOOD_FILE = path.get_full_path(__file__, "data", "good.json")
with open(GOOD_FILE, encoding="utf-8") as _handle:
    GOOD_DATA = _handle.read()


class TestValidation(unittest.TestCase):
    def test_full_good(self):
        results = loader.load_validated_json(GOOD_FILE, general._SCHEMA_FILE)
        assert results["tool"]["name"] == "test tool name"

    def test_bad_json(self):
        test_file = path.get_full_path(__file__, "data", "bad.json")
        with self.assertRaisesRegex(errors.AntismashInputError, "Expecting ',' delimiter"):
            loader.load_validated_json(test_file, general._SCHEMA_FILE)


class TestSchemaValidation(unittest.TestCase):
    def setUp(self):
        self.data = json.loads(GOOD_DATA)

    def is_valid(self, data=None, schema=None):
        if data is None:
            data = self.data
        if schema is None:
            schema = general._SCHEMA_FILE

        try:
            loader._ensure_valid(data, schema)
            return True
        except errors.AntismashInputError as err:
            assert "invalid sideload annotation" in str(err)
        return False

    def test_empty(self):
        assert not self.is_valid({})

    def test_bad_name(self):
        self.data["tool"].pop("name")
        assert not self.is_valid()

    def test_no_records(self):
        self.data["records"].clear()
        assert not self.is_valid()

    def test_detail_key(self):
        tool_config = self.data["tool"]["configuration"]
        for bad in [" key", "a", "_key", "ke,y", ""]:
            tool_config[bad] = "value"
            assert not self.is_valid(), bad
            tool_config.pop(bad)

        for good in ["good.key", "good_key", "good-key", "good-key-0"]:
            tool_config[good] = "value"
            assert self.is_valid(), good
            tool_config.pop(good)

    def test_detail_value(self):
        tool_config = self.data["tool"]["configuration"]
        for bad in ["_", " value", 0.7, None, 5, [], {}, "=b"]:
            tool_config["test"] = bad
            assert not self.is_valid(), bad

        for good in ["0.7", ".7", "-5", "-.5", "a=b,c=d", "symbols !@#$%^&*()", ["non-empty list"]]:
            tool_config["test"] = good
            assert self.is_valid(), good
