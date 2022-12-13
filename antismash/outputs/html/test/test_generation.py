# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import os
import unittest
from unittest.mock import Mock, patch

from helperlibs.wrappers.io import TemporaryDirectory

from antismash.main import get_all_modules
from antismash.common import html_renderer, path
from antismash.common.secmet.test.helpers import DummySubRegion
from antismash.common.test.helpers import DummyRecord
from antismash.config import build_config, destroy_config, update_config
from antismash.outputs.html import generator
from antismash.outputs.html.generator import (
    generate_webpage,
    build_antismash_js_url,
)


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
        with TemporaryDirectory() as temp_dir:
            update_config({"output_dir": temp_dir})
            assert generate_webpage([record], [{}], self.options, get_all_modules())


class TestJavascriptPresence(unittest.TestCase):
    def setUp(self):
        self.version = html_renderer.get_antismash_js_version()

    def build(self, return_values):
        options = Mock()
        options.database_dir = "dummy_data"
        with patch.object(path, "locate_file", side_effect=return_values) as patched:
            url = build_antismash_js_url(options)
            call_args = patched.call_args_list

        # ensure it wasn't called more times than expected
        assert len(call_args) == len(return_values)

        expected_args = [  # enough to cover all cases, only those used are checked
            path.get_full_path(generator.__file__, "js", "antismash.js"),
            os.path.join(options.database_dir, "as-js", self.version, "antismash.js"),
        ]
        for called, expected in zip(call_args, expected_args):
            assert os.path.realpath(called[0][0]) == os.path.realpath(expected)
        return url

    def test_fallback_to_web(self):
        url = self.build([None, None])
        assert url == html_renderer.get_antismash_js_url()

    def test_fallback_to_data(self):
        url = self.build([None, "data_dir"])
        assert url == "js/antismash.js"

    def test_fallback_to_in_tree(self):
        url = self.build(["local"])
        assert url == "js/antismash.js"
