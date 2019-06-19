# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Simple constructors for complicated features to simplify testing """

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

from ..features import (
    AntismashDomain,
    CDSFeature,
    CDSMotif,
    Protocluster,
    Feature,
    PFAMDomain,
    Region,
    SubRegion,
    CandidateCluster,
)
from ..features.candidate_cluster import CandidateClusterKind
from ..locations import FeatureLocation


class DummyAntismashDomain(AntismashDomain):
    counter = 0

    def __init__(self, start=0, end=3, strand=1, location=None, domain_id=None, tool="test_tool"):
        if location is None:
            location = FeatureLocation(start, end, strand=strand)
        if not domain_id:
            domain_id = "test_asDom_%d" % DummyAntismashDomain.counter
            DummyAntismashDomain.counter += 1
        super().__init__(location, tool=tool)
        self.domain_id = domain_id


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


class DummyFeature(Feature):
    def __init__(self, start=0, end=3, strand=1, feature_type="none"):
        super().__init__(FeatureLocation(start, end, strand), feature_type=feature_type)


class DummyProtocluster(Protocluster):
    def __init__(self, start=None, end=None, core_start=0, core_end=1,  # pylint: disable=too-many-arguments
                 core_location=None, tool="test", product="test-product",
                 cutoff=10, neighbourhood_range=10, high_priority_product=True):
        if core_location is None:
            core_location = FeatureLocation(core_start, core_end)
        if start is None:
            start = max(0, core_location.start - neighbourhood_range)
        if end is None:
            end = core_location.end + neighbourhood_range
        surrounds = FeatureLocation(start, end)
        super().__init__(core_location, surrounds, tool, product, cutoff,
                         neighbourhood_range, high_priority_product)


class DummyPFAMDomain(PFAMDomain):
    counter = 0

    def __init__(self, start=0, end=3, location=None, description="desc",  # pylint: disable=too-many-arguments
                 protein_start=0, protein_end=1, identifier=None, tool="test",
                 domain_id=None):
        if location is None:
            location = FeatureLocation(start, end)
        if identifier is None:
            identifier = "PF00001"
        super().__init__(location, description, protein_start, protein_end, identifier, tool)
        self.domain_id = domain_id or "dummy_pfam_%d" % DummyPFAMDomain.counter
        DummyPFAMDomain.counter += 1


class DummyRegion(Region):
    def __init__(self, candidate_clusters=None, subregions=None):
        if candidate_clusters is None:
            candidate_clusters = [DummyCandidateCluster()]
        if subregions is None:
            subregions = [DummySubRegion()]
        super().__init__(candidate_clusters, subregions)


class DummySubRegion(SubRegion):
    def __init__(self, start=0, end=10, location=None, tool="test",
                 label="test_label", probability=0.):
        if location is None:
            location = FeatureLocation(start, end)
        super().__init__(location, label=label, tool=tool, probability=probability)


class DummyCandidateCluster(CandidateCluster):
    def __init__(self, clusters=None, kind=None):
        if clusters is None:
            clusters = [DummyProtocluster()]
        if not kind:
            if len(clusters) == 1:
                kind = CandidateClusterKind.SINGLE
            else:
                kind = CandidateClusterKind.INTERLEAVED
        super().__init__(kind, clusters)
