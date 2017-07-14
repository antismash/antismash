# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

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
    def __init__(self, features=None, seq='FAKESEQ'):
        if features is None:
            features = []
        self.features = features
        self.seq = FakeSeq(seq)

    def __len__(self):
        """ returns the largest location of all features, so as to not break
            when new features are added to tests that extend past a hardcoded
            value
        """
        return max(max(feature.location.end, feature.location.start) for feature in self.features)

class FakeFeature(object):
    "class for generating a SeqFeature like datastructure"
    def __init__(self, feature_type, location=None, qualifiers=None):
        self.type = feature_type
        self.qualifiers = {} if qualifiers is None else qualifiers
        self.location = location

    def extract(self, seq):
        return seq

    def __repr__(self):
        return "FakeFeature(%r, %r, %r)" % (self.location, self.type,
                                            self.qualifiers)
