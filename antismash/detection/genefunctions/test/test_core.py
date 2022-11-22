# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import json
import unittest

import antismash
from antismash.common.hmmscan_refinement import HMMResult
from antismash.common.secmet.qualifiers import GeneFunction
from antismash.common.test.helpers import DummyRecord
from antismash.config import build_config, destroy_config
from antismash.detection import genefunctions


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


class SimpleResultsMixin:
    """ Requires:
         self.res_class  # the results class being tested
         self.results  # an instance of the results class
         self.record  # the record results were constructed with
    """

    def test_bad_record_id(self):
        data = json.loads(json.dumps(self.results.to_json()))
        self.record.id = "other"
        assert self.res_class.from_json(data, self.record) is None

        del data["record_id"]
        assert self.res_class.from_json(data, self.record) is None

    def test_bad_schema_version(self):
        data = json.loads(json.dumps(self.results.to_json()))
        data["schema_version"] += 1
        assert self.res_class.from_json(data, self.record) is None

        del data["schema_version"]
        assert self.res_class.from_json(data, self.record) is None


class TestFunctionResults(unittest.TestCase, SimpleResultsMixin):
    def setUp(self):
        self.res_class = genefunctions.core.FunctionResults
        hits = {"cds1": HMMResult("desc1", 0, 100, 2.3e-126, 416),
                "cds2": HMMResult("desc2", 5, 60, 3e-16, 20),
                }
        mapping = {"cds1": GeneFunction.TRANSPORT,
                   "cds2": GeneFunction.REGULATORY,
                   }
        self.record = DummyRecord()
        self.record.id = "rec_id"

        self.results = self.res_class(self.record.id, "toolname", best_hits=hits,
                                      function_mapping=mapping)

    def test_results_reconstruction(self):
        def check_results(results):
            assert results.record_id == "rec_id"
            assert results.tool == "toolname"
            assert isinstance(results.best_hits["cds1"], HMMResult)
            assert results.best_hits["cds1"].hit_id == 'desc1'
            assert results.best_hits["cds2"].bitscore == 20
            assert results.function_mapping["cds2"] == GeneFunction.REGULATORY

        hits = {"cds1": HMMResult("desc1", 0, 100, 2.3e-126, 416),
                "cds2": HMMResult("desc2", 5, 60, 3e-16, 20),
                }
        mapping = {"cds1": GeneFunction.TRANSPORT,
                   "cds2": GeneFunction.REGULATORY,
                   }
        results = self.res_class("rec_id", "toolname", best_hits=hits,
                                 function_mapping=mapping)
        check_results(results)

        data = json.loads(json.dumps(results.to_json()))

        record = DummyRecord()
        record.id = "rec_id"
        reconstructed = self.res_class.from_json(data, record)
        check_results(reconstructed)
        assert reconstructed.best_hits["cds1"].hit_id == hits["cds1"].hit_id


class TestAllFunctionResults(unittest.TestCase, SimpleResultsMixin):
    def setUp(self):
        self.res_class = genefunctions.AllFunctionResults
        self.record = DummyRecord()
        self.record.id = "rec_id"
        self.results = self.res_class(self.record.id)

    def test_invalid_type(self):
        with self.assertRaises(AssertionError):
            self.results.add_tool_results({"some": "results"})

    def test_duplicate_tools(self):
        single = genefunctions.core.FunctionResults(self.record.id, "toolname", {}, {})
        self.results.add_tool_results(single)

        with self.assertRaisesRegex(ValueError, "results already exist for tool"):
            self.results.add_tool_results(single)

    def test_results_reconstruction(self):
        def check_results(results, mapping):
            assert results.record_id == "rec_id"
            assert list(results._tools) == ["toolname"]
            assert results._tools["toolname"].tool == "toolname"
            assert results._tools["toolname"].function_mapping == mapping

        mapping = {"cds1": GeneFunction.TRANSPORT}
        tool_result = genefunctions.core.FunctionResults(self.record.id, "toolname", {}, mapping)
        self.results.add_tool_results(tool_result)
        check_results(self.results, mapping)

        data = json.loads(json.dumps(self.results.to_json()))
        assert data["tools"]
        reconstructed = self.res_class.from_json(data, self.record)
        check_results(reconstructed, mapping)

        reconstructed = genefunctions.regenerate_previous_results(data, self.record, None)
        check_results(reconstructed, mapping)
