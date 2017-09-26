# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

import os
import unittest

from helperlibs.bio import seqio
from helperlibs.wrappers.io import TemporaryDirectory

from antismash.common import secmet, path, subprocessing
import antismash.common.test.helpers as helpers
from antismash.config import args, get_config, update_config
from antismash.modules import smcogs

class Base(unittest.TestCase):
    def setUp(self):
        options = args.build_parser(modules=[smcogs]).parse_args(self.get_args())
        self.old_config = get_config().__dict__
        self.options = update_config(options)

        assert smcogs.check_prereqs() == []
        assert smcogs.check_options(self.options) == []
        assert smcogs.is_enabled(self.options)

        self.record = self.build_record(helpers.get_path_to_nisin_with_detection())

        def serial_run_func(function, arg_sets, _timeout=None):
            for arg_set in arg_sets:
                function(*arg_set)
        self.old_parallel = subprocessing.parallel_function
        subprocessing.parallel_function = serial_run_func

    def tearDown(self):
        subprocessing.parallel_function = self.old_parallel
        update_config({})

    def get_args(self):
        return ["--minimal", "--enable-smcogs"]

    def build_record(self, genbank):
        # construct a working record
        with open(genbank) as handle:
            seq_record = seqio.read(handle, "genbank")
        record = secmet.Record.from_biopython(seq_record)
        assert record.get_clusters()
        assert record.get_cluster(0).cds_children
        return record


class TestClassification(Base):
    def test_classifier(self):
        expected = open(path.get_full_path(__file__, "data", "nisin.txt")).readlines()
        with TemporaryDirectory(change=True):
            results = smcogs.run_on_record(self.record, None, self.options)
            contents = open("smcogs/smcogs.txt").readlines()
            assert contents == expected
            json = results.to_json()
            assert smcogs.SMCOGResults.from_json(json).to_json() == json

class TestTreeGeneration(Base):
    def get_args(self):
        return super().get_args() + ["--smcogs-trees"]

    def test_trees(self):
        with TemporaryDirectory(change=True):
            results = smcogs.run_on_record(self.record, None, self.options)
            assert len(results.tree_images) == 7
            for image in results.tree_images.values():
                assert os.path.exists(os.path.join(results.relative_tree_path, image))

            # test the results function properly
            json = results.to_json()
            assert smcogs.SMCOGResults.from_json(json).to_json() == json
            assert smcogs.check_previous_results(json, self.record, self.options).to_json() == json

            for cds in self.record.get_cluster(0).cds_children:
                hit = results.best_hits.get(cds.get_name())
                if hit:
                    assert not cds.notes
            results.add_to_record(self.record)
            for cds in self.record.get_cluster(0).cds_children:
                hit = results.best_hits.get(cds.get_name())
                if not hit:
                    continue
                assert cds.get_name() in results.tree_images
                assert len(cds.notes) == 2
