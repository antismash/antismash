# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Simple constructors for complicated features to simplify testing """

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=protected-access,missing-docstring

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
from ..record import Record


class DummyAntismashDomain(AntismashDomain):
    counter = 0

    def __init__(self, start=0, end=3, strand=1, location=None, domain_id=None, tool="test_tool",
                 protein_start=0, protein_end=1, protein_location=None, locus_tag="dummyCDS",
                 domain=None):
        if location is None:
            location = FeatureLocation(start, end, strand=strand)
        if protein_location is None:
            protein_location = FeatureLocation(protein_start, protein_end)
        if not domain_id:
            domain_id = f"test_asDom_{DummyAntismashDomain.counter}"
            DummyAntismashDomain.counter += 1
        super().__init__(location, tool, protein_location, locus_tag, domain=domain)
        self.domain_id = domain_id


class DummyCDS(CDSFeature):
    counter = 0

    def __init__(self, start=0, end=7, strand=1, locus_tag=None, translation=None):
        if not translation:
            translation = "A"*(abs(start-end))
        if not locus_tag:
            locus_tag = f"dummy_locus_tag_{DummyCDS.counter}"
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

    def __init__(self, start=0, end=6, strand=1, tool="test", domain_id=None,
                 protein_start=0, protein_end=1, protein_location=None, locus_tag="dummyCDS"):
        if not protein_location:
            protein_location = FeatureLocation(protein_start, protein_end)
        super().__init__(FeatureLocation(start, end, strand), locus_tag, protein_location, tool)
        if not domain_id:
            domain_id = f"dummy_domain{DummyCDSMotif.counter}_{start}_{end}"
            DummyCDSMotif.counter += 1
        self.domain_id = domain_id


class DummyFeature(Feature):
    def __init__(self, start=0, end=3, strand=1, feature_type="none"):
        super().__init__(FeatureLocation(start, end, strand), feature_type=feature_type)


class DummyProtocluster(Protocluster):
    def __init__(self, start=None, end=None, core_start=0, core_end=1,  # pylint: disable=too-many-arguments
                 core_location=None, tool="test", product="test-product",
                 cutoff=10, neighbourhood_range=10, high_priority_product=True,
                 product_category="TEST-CATEGORY"):
        if core_location is None:
            core_location = FeatureLocation(core_start, core_end)
        if start is None:
            start = max(0, core_location.start - neighbourhood_range)
        if end is None:
            end = core_location.end + neighbourhood_range
        surrounds = FeatureLocation(start, end)
        super().__init__(core_location, surrounds, tool, product, cutoff,
                         neighbourhood_range, high_priority_product,
                         product_category=product_category)


class DummyPFAMDomain(PFAMDomain):
    counter = 0

    def __init__(self, start=0, end=3, location=None, description="desc",  # pylint: disable=too-many-arguments
                 protein_start=0, protein_end=1, identifier=None, tool="test", domain=None,
                 domain_id=None, protein_location=None, locus_tag="dummyCDS"):
        if location is None:
            location = FeatureLocation(start, end)
        if identifier is None:
            identifier = "PF00001"
        if protein_location is None:
            protein_location = FeatureLocation(protein_start, protein_end)
        super().__init__(location, description, protein_location, identifier, tool, locus_tag, domain=domain)
        self.domain_id = domain_id or f"dummy_pfam_{DummyPFAMDomain.counter}"
        DummyPFAMDomain.counter += 1


class DummyRecord(Record):
    "class for generating a Record like data structure"
    def __init__(self, features=None, seq='AGCTACGT', taxon='bacteria',
                 record_id=None):
        if features:
            max_feature_coordinate = max(feature.location.end for feature in features)
            seq = seq * max(1, max_feature_coordinate // len(seq))
        super().__init__(seq, transl_table=11 if taxon == 'bacteria' else 1)
        if features:
            for feature in features:
                self.add_feature(feature)
        self.record_index = 0
        if record_id is not None:
            self.id = record_id


class DummyRegion(Region):
    def __init__(self, candidate_clusters=None, subregions=None):
        if candidate_clusters is None:
            candidate_clusters = [DummyCandidateCluster()]
        if subregions is None:
            subregions = [DummySubRegion()]
        super().__init__(candidate_clusters, subregions)


class DummySubRegion(SubRegion):
    def __init__(self, start=0, end=10, location=None, tool="test",
                 label="test_label"):
        if location is None:
            location = FeatureLocation(start, end)
        super().__init__(location, label=label, tool=tool)


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
