# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""
    Helper objects and functions for testing antismash
"""

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=protected-access,missing-docstring

# and some extra silencing because most of these classes are simple stubs
# pylint: disable=too-few-public-methods

import os
import stat

from helperlibs.wrappers.io import TemporaryDirectory

import antismash
from antismash.common import serialiser, path
from antismash.common.hmmscan_refinement import HMMResult
from antismash.common.module_results import ModuleResults
from antismash.common.secmet import Record
from antismash.common.secmet.test.helpers import (  # for import by others, pylint: disable=unused-import
    DummyAntismashDomain,
    DummyCandidateCluster,
    DummyCDS,
    DummyCDSMotif,
    DummyFeature,
    DummyPFAMDomain,
    DummyProtocluster,
    DummyRecord,
    DummyRegion,
)
from antismash.config import update_config
from antismash.config.args import build_parser
from antismash.main import get_all_modules, canonical_base_filename


def get_simple_options(module, args):
    modules = get_all_modules()
    if module is not None:
        modules = [module]
    return build_parser(from_config_file=False, modules=modules).parse_args(args)


class FakeHSP:
    def __init__(self, start, end, score):
        self.query_start = start
        self.query_end = end
        self.bitscore = score


class FakeHSPHit:
    "class for generating a HSP like datastructure"
    def __init__(self, query_id, hit_id, hit_start=0, hit_end=10, bitscore=10, evalue=1e-10, seeds=1):
        self.query_id = query_id
        self.hit_id = hit_id
        self.hit_start = hit_start
        self.hit_end = hit_end
        self.bitscore = bitscore
        self.evalue = evalue
        self.seeds = seeds

    def __repr__(self):
        return f"FakeHSP({vars(self)})"


class DummyHMMResult(HMMResult):
    def __init__(self, label="dummy", start=1, end=10, evalue=1, bitscore=1, internal_hits=None):
        super().__init__(label, start, end, evalue, bitscore, internal_hits=internal_hits)


class FakeHit:
    """ For generating hmmpfam2-like results """
    def __init__(self, start, end, score, desc):
        self.hsps = [FakeHSP(start, end, score)]
        self.description = desc

    def __iter__(self):
        return iter(self.hsps)


def get_path_to_nisin_genbank():
    file_path = __file__
    for _i in range(3):
        file_path = os.path.dirname(file_path)
    return os.path.join(file_path, 'test', 'integration', 'data', 'nisin.gbk')


def get_path_to_nisin_fasta():
    return path.get_full_path(__file__, "data", "nisin.fasta")


def get_path_to_nisin_with_detection():
    return get_path_to_nisin_genbank().replace('nisin', 'nisin_postdetection')


def get_path_to_balhymicin_genbank():
    file_path = __file__
    for _i in range(3):
        file_path = os.path.dirname(file_path)
    return os.path.join(file_path, 'test/integration/data/Y16952.gbk')


def _make_dir_deletable(base):
    # +rwx for all directories
    mode = stat.S_IRUSR | stat.S_IWUSR | stat.S_IXUSR
    os.chmod(base, mode)
    for root, dirs, _files in os.walk(base):
        for subdir in dirs:
            os.chmod(os.path.join(root, subdir), mode)


def run_and_regenerate_results_for_module(input_file, module, options,
                                          expected_record_count=1, callback=None):
    """ Runs antismash end to end over the given file with the given options
        and returns the given modules regenerated results

        if callback is supplied, it will be called with the output directory path
        as an argument before the output directory is cleared
    """
    with TemporaryDirectory(change=True) as tempdir:
        orig_output = options.output_dir
        update_config({"output_dir": tempdir})
        base_filename = canonical_base_filename(input_file, tempdir, options)
        json_filename = base_filename + ".json"
        assert not os.path.exists(json_filename)
        try:
            antismash.main.run_antismash(input_file, options)
        except Exception:
            update_config({"output_dir": orig_output})
            raise
        update_config({"output_dir": orig_output})
        results = serialiser.AntismashResults.from_file(json_filename)
        # remove things that were added by results, because otherwise the add isn't tested by detection
        # result regeneration
        for record in results.records:
            record.strip_antismash_annotations()
        if callback:
            callback(tempdir)
        # and while the genbank output still exists, grab that and check it's readable
        assert len(Record.from_genbank(base_filename + ".gbk")) == expected_record_count
        # for cases where sections of the source tree have no owner write permissions,
        # any directories that are copied from the source tree will also have no write permissions,
        # in which case the removal of the tempdir will fail unless those are changed
        _make_dir_deletable(tempdir)

    # not the responsibility of modules, but if it's wrong then everything is
    assert len(results.results) == expected_record_count
    assert len(results.records) == expected_record_count
    # ensure all detection stages add their relevant parts
    modules_to_regenerate = antismash.main.get_detection_modules()
    final = []
    for record, rec_results in zip(results.records, results.results):
        regenerate_results_for_record(record, options, modules_to_regenerate, module, rec_results)
        # post (other) detection has run, regenerate (since they may need regions etc)
        final.append(module.regenerate_previous_results(rec_results.get(module.__name__), record, options))
    for res in final:
        assert isinstance(res, ModuleResults)
    if expected_record_count == 1:
        return final[0]
    return final


def regenerate_results_for_record(record, options, modules, excluded_module, json_results):
    module_results = {}
    for module in modules:
        if module is excluded_module:
            continue
        json = json_results.get(module.__name__)
        if not json:
            continue
        mod_results = module.regenerate_previous_results(json, record, options)
        mod_results.add_to_record(record)
        if isinstance(mod_results, antismash.common.module_results.DetectionResults):
            for protocluster in mod_results.get_predicted_protoclusters():
                record.add_protocluster(protocluster)
            for region in mod_results.get_predicted_subregions():
                record.add_subregion(region)
        module_results[module.__name__] = mod_results
    record.create_candidate_clusters()
    record.create_regions()
    return module_results
