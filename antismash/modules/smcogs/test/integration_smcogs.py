# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import glob
import os
import unittest

from helperlibs.bio import seqio
from helperlibs.wrappers.io import TemporaryDirectory

import antismash
from antismash.common import secmet, path, subprocessing
import antismash.common.test.helpers as helpers
from antismash.config import get_config, update_config, destroy_config, build_config
from antismash.main import read_data, regenerate_results_for_record
from antismash.modules import smcogs


class Base(unittest.TestCase):
    def setUp(self):
        options = build_config(self.get_args(), isolated=True,
                                   modules=antismash.get_all_modules())
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
        destroy_config()
        update_config(self.old_config)

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
            assert smcogs.regenerate_previous_results(json, self.record, self.options).to_json() == json

            for cds in self.record.get_cluster(0).cds_children:
                hit = results.best_hits.get(cds.get_name())
                if hit:
                    assert not cds.notes
                    assert cds.gene_function in [secmet.feature.GeneFunction.OTHER,
                                                 secmet.feature.GeneFunction.CORE]
            results.add_to_record(self.record)
            for cds in self.record.get_cluster(0).cds_children:
                if cds.sec_met:
                    continue  # no sense checking, because we don't do anything with it
                hit = results.best_hits.get(cds.get_name())
                if not hit:
                    assert cds.gene_function == secmet.feature.GeneFunction.OTHER
                    continue
                assert cds.get_name() in results.tree_images
                assert len(cds.notes) == 1
                assert cds.gene_function != secmet.feature.GeneFunction.OTHER

    def test_trees_complete(self):
        with TemporaryDirectory() as output_dir:
            args = ["--minimal", "--smcogs-trees", "--output-dir", output_dir, helpers.get_path_to_nisin_genbank()]
            options = build_config(args, isolated=True, modules=antismash.get_all_modules())
            antismash.run_antismash(helpers.get_path_to_nisin_genbank(), options)


            with open(os.path.join(output_dir, "nisin.json")) as res_file:
                assert "antismash.modules.smcogs" in res_file.read()


            tree_files = list(glob.glob(os.path.join(output_dir, "smcogs", "*.png")))
            assert len(tree_files) == 7
            sample_tree = tree_files[0]

            # regen the results
            update_config({"reuse_results": os.path.join(output_dir, "nisin.json")})
            prior_results = read_data(None, options)
            record = prior_results.records[0]
            results = prior_results.results[0]
            regenned = regenerate_results_for_record(record, options, [smcogs], results)
            smcogs_results = regenned["antismash.modules.smcogs"]
            assert len(smcogs_results.tree_images) == 7
            assert os.path.exists(sample_tree)

            os.unlink(sample_tree)
            assert not os.path.exists(sample_tree)

            # attempt to regen the results, the deleted tree image will prevent it
            prior_results = read_data(None, options)
            record = prior_results.records[0]
            results = prior_results.results[0]
            regenned = regenerate_results_for_record(record, options, [smcogs], results)
            assert "antismash.modules.smcogs" not in regenned
