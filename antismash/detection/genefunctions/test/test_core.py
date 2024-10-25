# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import unittest

import antismash
from antismash.common import json
from antismash.common.secmet.qualifiers import GeneFunction
from antismash.common.test.helpers import DummyRecord
from antismash.config import build_config, destroy_config
from antismash.detection import genefunctions
from antismash.detection.genefunctions.tools.core import (
    HMMFunctionResults,
    HMMHit,
    Tool,
)


# allow for existing test code to run without modifications by adding a quick helper
# that converts old-style to new-style, including keyword-only arguments
def build_hit(identifier, start, end, evalue, bitscore, query_id="dummy_query", **kwargs):
    return HMMHit(query_id=query_id, reference_id=identifier, query_start=start, query_end=end,
                  evalue=evalue, bitscore=bitscore, **kwargs)


class TestOptions(unittest.TestCase):
    def setup_options(self, args):
        return build_config(args, isolated=True, modules=antismash.get_all_modules())

    def tearDown(self):
        destroy_config()

    def test_prereqs(self):
        options = self.setup_options([])  # required to setup database_dir and others
        assert genefunctions.check_prereqs(options) == []

    def test_default(self):
        options = self.setup_options([])
        assert not options.minimal
        assert genefunctions.check_options(options) == []
        assert genefunctions.is_enabled(options)

    def test_disabled(self):
        options = self.setup_options(["--minimal"])
        assert options.minimal
        assert not options.genefunctions_enabled
        assert genefunctions.check_options(options) == []
        assert not genefunctions.is_enabled(options)

    def test_minimal_but_enabled(self):
        options = self.setup_options(["--minimal", "--enable-genefunctions"])
        assert options.minimal
        assert options.genefunctions_enabled
        assert genefunctions.check_options(options) == []
        assert genefunctions.is_enabled(options)


class TestFunctionResults(unittest.TestCase):
    def setUp(self):
        self.res_class = HMMFunctionResults
        hits = {"cds1": build_hit("desc1", 0, 100, 2.3e-126, 416),
                "cds2": build_hit("desc2", 5, 60, 3e-16, 20),
                }
        mapping = {"cds1": GeneFunction.TRANSPORT,
                   "cds2": GeneFunction.REGULATORY,
                   }
        self.record = DummyRecord()

        self.results = self.res_class(tool="toolname", best_hits=hits,
                                      function_mapping=mapping)

    def test_results_reconstruction(self):
        def check_results(results):
            assert results.tool == "toolname"
            assert isinstance(results.best_hits["cds1"], HMMHit)
            assert results.best_hits["cds1"].reference_id == 'desc1'
            assert results.best_hits["cds2"].bitscore == 20
            assert results.function_mapping["cds2"] == GeneFunction.REGULATORY

        hits = {"cds1": build_hit("desc1", 0, 100, 2.3e-126, 416),
                "cds2": build_hit("desc2", 5, 60, 3e-16, 20),
                }
        mapping = {"cds1": GeneFunction.TRANSPORT,
                   "cds2": GeneFunction.REGULATORY,
                   }
        results = self.res_class(tool="toolname", best_hits=hits,
                                 function_mapping=mapping)
        check_results(results)

        data = json.loads(json.dumps(results.to_json()))

        reconstructed = self.res_class.from_json(data)
        check_results(reconstructed)
        assert reconstructed.best_hits["cds1"].reference_id == hits["cds1"].reference_id


class TestToolTracking(unittest.TestCase):
    def setUp(self):
        self.results = genefunctions.AllFunctionResults(record_id="rec")
        self.tool = Tool(name="dummy", classify=lambda x: None)

    def test_unknown_tool(self):
        tool_results = HMMFunctionResults(tool=self.tool.name, best_hits={}, function_mapping={})
        with self.assertRaisesRegex(NotImplementedError, "no handling for"):
            self.results.add_tool_results(self.tool, tool_results)
