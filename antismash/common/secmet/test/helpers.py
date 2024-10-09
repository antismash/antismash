# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Simple constructors for complicated features to simplify testing """

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=protected-access,missing-docstring

from Bio.Seq import Seq

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
from ..locations import CompoundLocation, FeatureLocation, Location
from ..record import Record, Seq


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

    def __init__(self, start=0, end=7, strand=1, locus_tag=None, translation=None, location=None):
        if not translation:
            translation = "A"*(abs(start-end))
        if not locus_tag:
            locus_tag = f"dummy_locus_tag_{DummyCDS.counter}"
            DummyCDS.counter += 1
        if location is None:
            location = FeatureLocation(start, end, strand)
        super().__init__(location, translation=translation,
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
                 product_category="TEST-CATEGORY", record_length=None):
        if core_location is None:
            if core_start > core_end:
                assert record_length
                core_location = CompoundLocation([FeatureLocation(core_start, record_length),
                                                  FeatureLocation(0, core_end)])
            else:
                core_location = FeatureLocation(core_start, core_end)
        core_start = core_location.parts[0].start
        core_end = core_location.parts[-1].end
        if start is None:
            start = core_start - neighbourhood_range
            if start < 0 and record_length is not None:
                start += record_length
            elif start < 0:
                start = 0
        if end is None:
            end = core_end + neighbourhood_range
            if record_length is not None and end > record_length:
                end %= record_length
        if start > end:
            assert record_length
            surrounds = CompoundLocation([FeatureLocation(start, record_length, 1),
                                          FeatureLocation(0, end, 1)])
        else:
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
    def __init__(self, features=None, seq=None, taxon='bacteria',
                 record_id=None, *, length=None, circular=False):
        if features and seq is None:
            length = length or max(f.location.end for f in features)
            seq = "AGCTACGT"
        self.length = length
        if seq is None:
            seq = "AGCTACGT"
        if isinstance(seq, str):
            seq = Seq(seq)
        super().__init__(seq, transl_table=11 if taxon == 'bacteria' else 1)
        if circular:
            self.make_circular()
        if features:
            for feature in features:
                self.add_feature(feature)
        self.record_index = 0
        if record_id is not None:
            self.id = record_id
        # pylint doesn't recognise the superclass has a sequence
        # so trick it into believing it exists
        self.seq = Seq(seq)

    @property
    def features(self):
        return self.get_all_features()

    def make_circular(self, circular=True):
        if not self._record.annotations:
            self._record.annotations = {}
        if circular:
            self._record.annotations["topology"] = "circular"
            assert self.is_circular()
        else:
            self._record.annotations["topology"] = "linear"
            assert not self.is_circular()

    def rotate(self, cut_point: int, padding: int = 0) -> None:
        return rotate(self, cut_point, padding=padding)

    def __len__(self):  # override the length so the sequence doesn't necessarily have to exist
        return self.length or super().__len__()


class DummyRegion(Region):
    def __init__(self, candidate_clusters=None, subregions=None, start=0, end=100):
        if not candidate_clusters and not subregions:
            if candidate_clusters is None:
                candidate_clusters = [DummyCandidateCluster(start=start, end=end)]
            if subregions is None:
                subregions = [DummySubRegion(start=start, end=end)]
        super().__init__(candidate_clusters, subregions)


class DummySubRegion(SubRegion):
    def __init__(self, start=0, end=10, location=None, tool="test",
                 label="test_label", record_length=None):
        if location is None:
            if record_length is not None and start > end:
                location = CompoundLocation([
                    FeatureLocation(start, record_length, 1),
                    FeatureLocation(0, end, 1),
                ])
            else:
                location = FeatureLocation(start, end, 1)
        super().__init__(location, label=label, tool=tool)


class DummyCandidateCluster(CandidateCluster):
    def __init__(self, clusters=None, kind=None, start=0, end=100, **kwargs):
        if clusters is None:
            clusters = [DummyProtocluster(start=start, end=end, record_length=kwargs.get("circular_wrap_point"))]
        if not kind:
            if len(clusters) == 1:
                kind = CandidateClusterKind.SINGLE
            else:
                kind = CandidateClusterKind.INTERLEAVED
        if "circular_wrap_point" not in kwargs and any(cluster.crosses_origin() for cluster in clusters):
            kwargs["circular_wrap_point"] = max(cluster.location.parts[0].end for cluster in clusters)
        super().__init__(kind, clusters, **kwargs)

    def get_candidate_cluster_number(self):
        # prevent failures when testing candidates in isolation from records
        try:
            return super().get_candidate_cluster_number()
        except ValueError:
            return -1


def rotate(record: Record, cut_point: int, padding: int = 0) -> None:
    if len(record) < cut_point:
        raise ValueError("Invalid cut location: outside record")
    seq = record.seq  # pylint struggles with recognising this exists, so # pylint: disable=no-member
    new_seq = seq[cut_point:] + Seq("A" * padding) + seq[:cut_point]
    assert len(new_seq) == len(seq) + padding

    for feature in record.all_features:
        if feature.type == "source":
            feature.location = FeatureLocation(
                feature.location.start,
                feature.location.end + padding,
                feature.location.strand
            )
            continue
        location = feature.location
        if cut_point in location:
            location = split_location_on_cut(location, cut_point)
        location = shift_location(location, cut_point, padding, len(record))
        feature.location = location
    record._cds_features = sorted(record._cds_features)
    record.seq = Seq(new_seq)
    record._record.annotations["topology"] = "circular"


def split_location_on_cut(location: Location, cut_point: int) -> CompoundLocation:
    new_parts = []
    for part in location.parts:
        if cut_point not in part:
            new_parts.append(part)
            continue
        new_parts.append(FeatureLocation(part.start, cut_point, part.strand))
        new_parts.append(FeatureLocation(cut_point, part.end, part.strand))
    return CompoundLocation(new_parts)


def test_split():
    assert split_location_on_cut(FeatureLocation(5, 12, -1), 7).parts == [
        FeatureLocation(5, 7, -1), FeatureLocation(7, 12, -1)
    ]
    compound = CompoundLocation([FeatureLocation(5, 12, 1), FeatureLocation(15, 21, 1)])
    assert split_location_on_cut(compound, 13) == compound

# ensure the splitting function works as expected, since other tests rely on it
# and since it won't run automatically due to file naming, manually test it here
test_split()

def shift_location(location: Location, cut_point: int, padding: int, rec_len: int) -> Location:
    new_parts = []
    for part in location.parts:
        if part.start < cut_point:
            new_parts.append(FeatureLocation(part.start + (rec_len - cut_point) + padding,
                                             part.end + (rec_len - cut_point) + padding,
                                             part.strand))
        else:
            new_parts.append(FeatureLocation(part.start - cut_point,
                                             part.end - cut_point,
                                             part.strand))
    if len(new_parts) == 1:
        return new_parts[0]
    return CompoundLocation(new_parts)


# ensure the shift function works as expected, since other tests rely on it
assert shift_location(FeatureLocation(7, 8, -1), 15, 10, 20) == FeatureLocation(22, 23, -1)
assert shift_location(FeatureLocation(16, 17, 1), 15, 10, 20) == FeatureLocation(1, 2, 1)
