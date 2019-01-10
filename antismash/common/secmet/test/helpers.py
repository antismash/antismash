# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Simple constructors for complicated features to simplify testing """

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

from ..features import CDSFeature, CDSMotif
from ..locations import FeatureLocation


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

    # bypass translation checks for these artifical CDSs
    @property
    def translation(self):
        return self._translation

    @translation.setter
    def translation(self, translation):
        self._translation = translation  # pylint: disable=attribute-defined-outside-init


class DummyCDSMotif(CDSMotif):
    counter = 0

    def __init__(self, start=0, end=6, strand=1, tool=None, domain_id=None):
        super().__init__(FeatureLocation(start, end, strand), tool)
        if not domain_id:
            domain_id = "dummy_domain%d_%d_%d" % (DummyCDSMotif.counter, start, end)
            DummyCDSMotif.counter += 1
        self.domain_id = domain_id
