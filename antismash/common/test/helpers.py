# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""
    Helper objects and functions for testing antismash
"""

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

# and some extra silencing because most of these classes are simple stubs
# pylint: disable=too-few-public-methods

import os

from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation
from helperlibs.wrappers.io import TemporaryDirectory

import antismash
from antismash.common import serialiser, module_results, path
from antismash.common.secmet import Cluster, CDSFeature, Feature, Record
from antismash.config import update_config
from antismash.config.args import build_parser
from antismash.main import get_all_modules


class DummyFeature(Feature):
    def __init__(self, start, end, strand=1, feature_type="none"):
        super().__init__(FeatureLocation(start, end, strand), feature_type=feature_type)


class DummyCDS(CDSFeature):
    counter = 0

    def __init__(self, start=0, end=7, strand=1, locus_tag=None, translation=None):
        if not translation:
            translation = "A"*(abs(start-end))
        if not locus_tag:
            locus_tag = "dummy_locus_tag_%d" % DummyCDS.counter
            DummyCDS.counter += 1
        super().__init__(FeatureLocation(start, end, strand), translation=translation,
                         locus_tag=locus_tag)
        assert self.get_accession() == locus_tag, self.get_accession()


class DummyCluster(Cluster):
    def __init__(self, start, end, strand=1):
        cutoff = 10
        extent = 10
        product = ["dummy"]
        super().__init__(FeatureLocation(start, end, strand), cutoff, extent,
                         product)


def get_simple_options(module, args):
    modules = get_all_modules()
    if module is not None:
        modules = [module]
    return build_parser(from_config_file=False, modules=modules).parse_args(args)


class DummyRecord(Record):
    "class for generating a Record like data structure"
    def __init__(self, features=None, seq='FAKESEQ'):
        if isinstance(seq, str):
            seq = Seq(seq)
        super().__init__(seq)
        if features:
            for feature in features:
                self.add_feature(feature)
        self.record_index = 0


class FakeHSP(object):
    def __init__(self, start, end, score):
        self.query_start = start
        self.query_end = end
        self.bitscore = score


class FakeHSPHit:
    "class for generating a HSP like datastructure"
    def __init__(self, query_id, hit_id, hit_start, hit_end, bitscore, evalue):
        self.query_id = query_id
        self.hit_id = hit_id
        self.hit_start = hit_start
        self.hit_end = hit_end
        self.bitscore = bitscore
        self.evalue = evalue

    def __repr__(self):
        return "FakeHSP({})".format(str(vars(self)))


class FakeHit(object):
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


def run_and_regenerate_results_for_module(input_file, module, options,
                                          expected_record_count, callback=None):
    """ Runs antismash end to end over the given file with the given options
        and returns the given modules regenerated results

        if callback is supplied, it will be called with the output directory path
        as an argument before the output directory is cleared
    """
    with TemporaryDirectory(change=True) as tempdir:
        orig_output = options.output_dir
        update_config({"output_dir": tempdir})
        json_filename = os.path.join(options.output_dir, os.path.basename(input_file).rsplit('.', 1)[0] + ".json")
        assert not os.path.exists(json_filename)
        try:
            antismash.main.run_antismash(input_file, options)
        except:
            update_config({"output_dir": orig_output})
            raise
        update_config({"output_dir": orig_output})
        results = serialiser.AntismashResults.from_file(json_filename)
        if callback:
            callback(tempdir)
    # not the responsibility of modules, but if it's wrong then everything is
    assert len(results.results) == expected_record_count
    assert len(results.records) == expected_record_count
    if expected_record_count == 1:
        regenerated = antismash.main.regenerate_results_for_record(results.records[0],
                                     options, [module], results.results[0])
        print(list(regenerated))
        final = regenerated[module.__name__]
        assert isinstance(final, module_results.ModuleResults)
    else:
        regenerated = [antismash.main.regenerate_results_for_record(record, options,
                        [module], res) for record, res in zip(results.records, results.results)]
        final = [result[module.__name__] for result in regenerated]
        for res in final:
            assert isinstance(res, module_results.ModuleResults)
    return final
