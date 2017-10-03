# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""
    Helper objects for testing antismash
"""

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

# and some extra silencing because most of these classes are simple stubs
# pylint: disable=too-few-public-methods

import os

from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation

from antismash.common.secmet import Cluster, CDSFeature, Feature, Record
from antismash.config.args import build_parser
from antismash.main import get_all_modules

class DummyFeature(Feature):
    def __init__(self, start, end, strand=1, feature_type="none"):
        super().__init__(FeatureLocation(start, end, strand), feature_type=feature_type)

class DummyCDS(CDSFeature):
    counter = 0
    def __init__(self, start, end, strand=1, locus_tag=None, translation=None):
        trans = "A"*(abs(start-end))
        if not locus_tag:
            locus_tag = "dummy_locus_tag_%d" % DummyCDS.counter
            DummyCDS.counter += 1
        if translation:
            self._translation = translation
        super().__init__(FeatureLocation(start, end, strand), translation=trans,
                         locus_tag=locus_tag)

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
        super().__init__(Seq(seq))
        if features:
            for feature in features:
                self.add_feature(feature)
        self.record_index = 0

def get_path_to_nisin_genbank():
    path = __file__
    for _i in range(3):
        path = os.path.dirname(path)
    return os.path.join(path, 'test', 'integration', 'data', 'nisin.gbk')

def get_path_to_nisin_with_detection():
    return get_path_to_nisin_genbank().replace('nisin', 'nisin_postdetection')

def get_path_to_balhymicin_genbank():
    path = __file__
    for _i in range(3):
        path = os.path.dirname(path)
    return os.path.join(path, 'test/integration/data/Y16952.gbk')
