# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

import os

from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation
from antismash.common.secmet import Record, Cluster, CDSFeature, Feature

class DummyFeature(Feature):
    def __init__(self, start, end, strand=1):
        super().__init__(FeatureLocation(start, end, strand), feature_type="none")

class DummyCDS(CDSFeature):
    def __init__(self, start, end, strand=1):
        trans = "dummy_translation"
        locus_tag = "dummy_locus_tag"
        super().__init__(FeatureLocation(start, end, strand), translation=trans,
                         locus_tag=locus_tag)

class DummyCluster(Cluster):
    def __init__(self, start, end, strand=1):
        cutoff = 10
        extent = 10
        product = ["dummy"]
        super().__init__(FeatureLocation(start, end, strand), cutoff, extent,
                         product)

class FakeSeq(object):
    "class for generating a Seq like datastructure"
    def __init__(self, seq):
        self.seq = seq

    def translate(self, dummy):
        return self.seq

    def __str__(self):
        return self.seq

class FakeRecord(object):
    "class for generating a seq_record like data structure"
    def __init__(self, features=None, seq='FAKESEQ', real_seq=False):
        self.record_index = 0
        if features is None:
            features = []
        self.features = features
        if real_seq:
            self.seq = Seq(seq)
        else:
            self.seq = FakeSeq(seq)

    def __len__(self):
        """ returns the largest location of all features, so as to not break
            when new features are added to tests that extend past a hardcoded
            value

            if no features exist yet, returns the length of the sequence
        """
        if not self.features:
            return len(self.seq)
        return max(max(feature.location.end, feature.location.start) for feature in self.features)

    def get_cds_features(self):
        res = []
        for feature in self.features:
            if feature.type == "CDS":
                res.append(feature)
        return res

    def add_cluster(self, cluster):
        self.features.append(cluster)

    def get_clusters(self):
        res = []
        for feature in self.features:
            if feature.type == "cluster":
                res.append(feature)
        return res

class FakeFeature(object):
    "class for generating a SeqFeature like datastructure"
    def __init__(self, feature_type, location=None, qualifiers=None):
        self.type = feature_type
        self.qualifiers = { "translation" : ["trans"]}
        if qualifiers:
            self.qualifiers.update(qualifiers)
        self.location = location

    def __getattr__(self, attr):
        return self.__dict__.get(attr, [])

    def extract(self, seq):
        return seq

    def __repr__(self):
        return "FakeFeature(%r, %r, %r)" % (self.location, self.type,
                                            self.qualifiers)


def get_path_to_nisin_genbank():
    path = __file__
    for i in range(3):
        path = os.path.dirname(path)
    return os.path.join(path, 'test/integration/data/nisin.gbk')
