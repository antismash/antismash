# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Implements a more readable and efficient version of biopython's SeqRecord.

    Most attributes are wrapped, with the record features being secmet features
    and accessible and already paritioned so accessing specific types of feature
    is more efficient.

    Does not, and will not, subclass SeqRecord to avoid mixing code types.
"""


import bisect
from collections import Counter, defaultdict, OrderedDict
import itertools
import logging
from typing import (
    Any,
    Dict,
    Iterable,
    List,
    Optional,
    Sequence,
    Tuple,
    Type,
    TypeVar,
    Union,
    cast,
)
from zlib import crc32

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord, UndefinedSequenceError

from .errors import SecmetInvalidInputError
from .features import (
    AntismashDomain,
    CandidateCluster,
    CDSFeature,
    CDSMotif,
    Domain,
    Feature,
    Gene,
    Module,
    PFAMDomain,
    Prepeptide,
    Protocluster,
    Region,
    Source,
    SubRegion,
)
from .features import CDSCollection
from .features.candidate_cluster import create_candidates_from_protoclusters

from .locations import (
    CompoundLocation,
    FeatureLocation,
    Location,
    connect_locations,
    get_distance_between_locations,
    location_bridges_origin,
    locations_overlap,
    ensure_valid_locations,
    remove_redundant_exons,
)

T = TypeVar("T", bound="Record")


ANTISMASH_SPECIFIC_TYPES: List[str] = [feature_class.FEATURE_TYPE for feature_class in [  # type: ignore
    AntismashDomain,
    CandidateCluster,
    Module,
    Protocluster,
    Region,
    SubRegion,
]]


class Record:
    """A record containing secondary metabolite clusters"""
    # slots not for space, but to stop use as a horrible global
    __slots__ = ["_record", "_seq", "skip", "_cds_features", "_cds_by_name",
                 "_cds_by_location",
                 "_protoclusters", "original_id", "_cds_motifs", "_pfam_domains",
                 "_antismash_domains", "_protocluster_numbering", "_nonspecific_features",
                 "record_index", "_genes", "_transl_table", "_domains_by_name",
                 "_pfams_by_cds_name", "_genes_by_name", "_modules",
                 "_candidate_clusters", "_candidate_clusters_numbering",
                 "_subregions", "_subregion_numbering",
                 "_regions", "_region_numbering", "_antismash_domains_by_tool",
                 "_antismash_domains_by_cds_name", "_gc_content",
                 "_cds_cache", "_cds_cache_dirty", "_sources",
                 ]

    def __init__(self, seq: Union[Seq, str] = "", *,
                 transl_table: int = 1, gc_content: float = -1.,
                 **kwargs: Any,
                 ) -> None:
        # prevent paths from being used as a sequence
        assert not set("./\\").issubset(set(seq)), "Invalid sequence provided"
        if isinstance(seq, str):
            seq = Seq(seq)
        assert isinstance(seq, Seq)
        self._record = SeqRecord(seq, **kwargs)
        self.record_index: Optional[int] = None
        self.original_id = None
        self.skip: Optional[str] = None  # TODO: move to yet another abstraction layer?

        self._genes: List[Gene] = []
        self._genes_by_name: Dict[str, List[Gene]] = defaultdict(list)

        self._cds_features: List[CDSFeature] = []
        self._cds_by_name: Dict[str, CDSFeature] = {}
        self._cds_by_location: Dict[str, CDSFeature] = {}
        self._cds_cache: tuple[CDSFeature, ...] = tuple()
        self._cds_cache_dirty: bool = False

        self._cds_motifs: List[CDSMotif] = []

        self._pfam_domains: List[PFAMDomain] = []
        self._pfams_by_cds_name: Dict[str, List[PFAMDomain]] = defaultdict(list)

        self._antismash_domains: List[AntismashDomain] = []
        self._antismash_domains_by_tool: Dict[str, List[AntismashDomain]] = defaultdict(list)
        self._antismash_domains_by_cds_name: Dict[str, List[AntismashDomain]] = defaultdict(list)

        self._modules: List[Module] = []

        # includes PFAMDomains and AntismashDomains
        self._domains_by_name: Dict[str, Domain] = {}  # for use as x[domain.get_name()] = domain

        self._nonspecific_features: List[Feature] = []

        self._protoclusters: List[Protocluster] = []
        self._protocluster_numbering: Dict[Protocluster, int] = {}

        self._candidate_clusters: List[CandidateCluster] = []
        self._candidate_clusters_numbering: Dict[CandidateCluster, int] = {}

        self._subregions: List[SubRegion] = []
        self._subregion_numbering: Dict[SubRegion, int] = {}

        self._regions: List[Region] = []
        self._region_numbering: Dict[Region, int] = {}

        self._transl_table = int(transl_table)
        self._gc_content: float = gc_content

        self._sources: List[Source] = []

    def __getattr__(self, attr: str) -> Any:
        # passthroughs to the original SeqRecord
        if attr in ["id", "seq", "description", "name", "annotations", "dbxrefs"]:
            return getattr(self._record, attr)
        if attr in Record.__slots__:
            return self.__getattribute__(attr)
        raise AttributeError(f"Record has no attribute {attr!r}")

    def __setattr__(self, attr: str, value: Any) -> None:
        # passthroughs to the original SeqRecord
        if attr in ["id", "seq", "description", "name"]:
            setattr(self._record, attr, value)
            return
        if attr in ["annotations"]:
            assert isinstance(value, dict)
            for key, val in value.items():
                self.add_annotation(key, val)
            return
        # something Record owns or shouldn't have
        try:
            super().__setattr__(attr, value)
        except AttributeError:
            raise AttributeError("Record does not support dynamically adding attributes")

    @property
    def transl_table(self) -> int:
        """ Returns the default translation table used throughout the Record """
        return self._transl_table

    def is_circular(self) -> bool:
        """ Returns True if the genome is circular """
        return self._record.annotations.get("topology", "").lower() == "circular"

    def add_annotation(self, key: str, value: List) -> None:
        """Adding annotations in Record"""
        if not isinstance(key, str) or not isinstance(value, (str, list)):
            raise ValueError('Key and Value are not in right format')
        self._record.annotations[key] = value

    def __len__(self) -> int:
        return len(self._record)

    def has_name(self, name: str) -> bool:
        """ Returns True if the given name matches any of the identifiers of the
            record
        """
        if name == self.id:
            return True
        if self.original_id is not None:
            return name == self.original_id
        return False

    # protocluster manipulation
    def add_protocluster(self, cluster: Protocluster) -> None:
        """ Add the given cluster to the record,
            causes cluster-CDS pairing to be recalculated """
        assert isinstance(cluster, Protocluster), type(cluster)
        assert cluster.location.start >= 0, cluster
        assert cluster.location.end <= len(self), f"{cluster} > {len(self)}"
        index = bisect.bisect_left(self._protoclusters, cluster)
        self._protoclusters.insert(index, cluster)
        cluster.parent_record = self
        # update numbering
        for i in range(index, len(self._protoclusters)):
            self._protocluster_numbering[self._protoclusters[i]] = i + 1  # 1-indexed
        for cds in self.get_cds_features_within_location(cluster.location):
            cluster.add_cds(cds)

    def get_protoclusters(self) -> Tuple[Protocluster, ...]:
        """A list of secondary metabolite clusters present in the record"""
        return tuple(self._protoclusters)

    def get_protocluster(self, index: int) -> Protocluster:
        """ Get the cluster with the given cluster number """
        return self._protoclusters[index - 1]  # change from 1-indexed to 0-indexed

    def get_protocluster_number(self, protocluster: Protocluster) -> int:
        """Returns protocluster number of a protocluster feature (1-indexed)
        """
        number = self._protocluster_numbering.get(protocluster)
        if number is None:
            raise ValueError("Protocluster not contained in record")
        return number

    def clear_protoclusters(self) -> None:
        """ Removes all Protocluster features.
            This also removes CandidateCluster features and resets all Regions.
        """
        self._protoclusters.clear()
        self.clear_candidate_clusters()

    # candidate cluster manipulation
    def add_candidate_cluster(self, cluster: CandidateCluster) -> None:
        """ Add the given CandidateCluster to the record """
        assert isinstance(cluster, CandidateCluster), type(cluster)
        assert cluster.location.start >= 0, cluster
        assert cluster.location.end <= len(self), f"{cluster} > {len(self)}"
        index = bisect.bisect_left(self._candidate_clusters, cluster)
        self._candidate_clusters.insert(index, cluster)
        cluster.parent_record = self
        # update numbering
        for i in range(index, len(self._candidate_clusters)):
            self._candidate_clusters_numbering[self._candidate_clusters[i]] = i + 1  # 1-indexed
        for cds in self.get_cds_features_within_location(cluster.location):
            cluster.add_cds(cds)

    def get_candidate_clusters(self) -> Tuple[CandidateCluster, ...]:
        """A list of secondary metabolite clusters present in the record"""
        return tuple(self._candidate_clusters)

    def get_candidate_cluster(self, index: int) -> CandidateCluster:
        """ Get the cluster with the given cluster number """
        return self._candidate_clusters[index - 1]  # change from 1-indexed to 0-indexed

    def get_candidate_cluster_number(self, cluster: CandidateCluster) -> int:
        """Returns candidate_clusters number of a CandidateCluster feature (1-indexed)
        """
        number = self._candidate_clusters_numbering.get(cluster)
        if number is None:
            raise ValueError("CandidateCluster not contained in record")
        return number

    def clear_candidate_clusters(self) -> None:
        """ Removes all CandidateCluster features and resets all Regions
        """
        for candidate_clusters in self._candidate_clusters:
            for cluster in candidate_clusters.protoclusters:
                cluster.parent = None
        self._candidate_clusters.clear()
        if self._regions:  # only recreate if they existed to start with
            self.clear_regions()
            self.create_regions()

    # subregion manipulation
    def add_subregion(self, subregion: SubRegion) -> None:
        """ Add the given SubRegion to the record
        """
        assert isinstance(subregion, SubRegion), type(subregion)
        assert subregion.location.start >= 0, subregion
        assert subregion.location.end <= len(self), f"{subregion} > {len(self)}"
        index = bisect.bisect_left(self._subregions, subregion)
        self._subregions.insert(index, subregion)
        subregion.parent_record = self
        # update numbering
        for i in range(index, len(self._subregions)):
            self._subregion_numbering[self._subregions[i]] = i + 1  # 1-indexed
        for cds in self.get_cds_features_within_location(subregion.location):
            subregion.add_cds(cds)

    def get_subregions(self) -> Tuple[SubRegion, ...]:
        """A list of subregions present in the record"""
        return tuple(self._subregions)

    def get_subregion(self, index: int) -> SubRegion:
        """ Get the subregion with the given cluster number """
        return self._subregions[index - 1]  # change from 1-indexed to 0-indexed

    def get_subregion_number(self, region: SubRegion) -> int:
        """Returns subregion number of a SubRegion feature (1-indexed)
        """
        number = self._subregion_numbering.get(region)
        if number is None:
            raise ValueError("SubRegion not contained in record")
        return number

    def clear_subregions(self) -> None:
        """ Removes all SubRegion features and resets all Regions depending on them
        """
        self._subregions.clear()
        if self._regions:  # only recreate if they existed to start with
            self.clear_regions()
            self.create_regions()

    # region manipulation
    def add_region(self, region: Region) -> None:
        """ Add the given Region to the record
        """
        assert isinstance(region, Region), type(region)
        assert region.location.start >= 0, region
        assert region.location.end <= len(self), f"{region} > {len(self)}"
        index = 0
        for i, existing_region in enumerate(self._regions):  # TODO: fix performance
            if region.overlaps_with(existing_region):
                logging.error("existing %s overlaps with new %s", existing_region, region)
                raise ValueError("regions cannot overlap")
            if region < existing_region:
                index = i  # before
                break
            index = i + 1  # after
        self._regions.insert(index, region)
        region.parent_record = self
        # update numbering
        for i in range(index, len(self._regions)):
            self._region_numbering[self._regions[i]] = i + 1  # 1-indexed
        # link any relevant CDS features
        for cds in self.get_cds_features_within_location(region.location):
            region.add_cds(cds)
            cds.region = region

    def get_regions(self) -> Tuple[Region, ...]:
        """ The Region features in the record representing regions of interest """
        return tuple(self._regions)

    def get_region(self, index: int) -> Region:
        """ Get the Region with the given cluster number """
        return self._regions[index - 1]  # change from 1-indexed to 0-indexed

    def get_region_number(self, region: Region) -> int:
        """Returns subregion number of a Region feature (1-indexed)
        """
        number = self._region_numbering.get(region)
        if number is None:
            raise ValueError("Region not contained in record")
        return number

    def clear_regions(self) -> None:
        "Remove all Region features"
        for region in self._regions:
            for cds in region.cds_children:
                cds.region = None
            for cluster in region.candidate_clusters:
                cluster.parent = None
            for subregion in region.subregions:
                subregion.parent = None
        self._regions.clear()
    # end of large manipulation sections

    def get_genes(self) -> Tuple[Gene, ...]:
        """ Returns all Gene features (different to CDSFeatures) in the record """
        return tuple(self._genes)

    def get_genes_by_name(self, name: str) -> List[Gene]:
        """ Return the genes with a given name """
        return list(self._genes_by_name[name])

    def get_cds_features(self) -> Tuple[CDSFeature, ...]:
        """A list of secondary metabolite clusters present in the record"""
        if self._cds_cache_dirty or not self._cds_features:
            self._cds_cache = tuple(self._cds_features)
            self._cds_cache_dirty = False
        return self._cds_cache

    def get_cds_name_mapping(self) -> Dict[str, CDSFeature]:
        """A dictionary mapping CDS name to CDS feature"""
        return dict(self._cds_by_name)

    def get_cds_by_name(self, name: str) -> CDSFeature:
        """ Return the CDS with the given name """
        return self._cds_by_name[name]

    def get_domain_by_name(self, name: str) -> Domain:
        """ Return the Domain with the given name """
        try:
            return self._domains_by_name[name]
        except KeyError:
            raise KeyError(f"record {self.id} contains no domain named {name}")

    def get_cds_motifs(self) -> Tuple[CDSMotif, ...]:
        """A list of secondary metabolite CDS_motifs present in the record"""
        return tuple(self._cds_motifs)

    def clear_cds_motifs(self) -> None:
        "Remove all CDSMotif features"
        for motif in self._cds_motifs:
            del self._domains_by_name[motif.get_name()]
        self._cds_motifs.clear()

    def get_pfam_domains(self) -> Tuple[PFAMDomain, ...]:
        """A list of secondary metabolite PFAM_domains present in the record"""
        return tuple(self._pfam_domains)

    def get_pfam_domains_in_cds(self, cds: Union[str, CDSFeature]) -> Tuple[PFAMDomain, ...]:
        """ Returns a list of PFAMDomains contained by a CDSFeature. Either the
            CDS name or the CDSFeature itself can be used to specify which CDS.
        """
        if isinstance(cds, CDSFeature):
            cds_name = cds.get_name()
        elif isinstance(cds, str):
            cds_name = cds
        else:
            raise TypeError(f"CDS must be a string or CDSFeature, not {type(cds)}")
        return tuple(self._pfams_by_cds_name[cds_name])

    def clear_pfam_domains(self) -> None:
        """ Remove all PFAMDomain features """
        self._pfams_by_cds_name.clear()
        # remove pfams only from the domains mapping
        for domain in self._pfam_domains:
            del self._domains_by_name[domain.get_name()]
        self._pfam_domains.clear()

    def get_antismash_domains(self) -> Tuple[AntismashDomain, ...]:
        """A list of secondary metabolite aSDomains present in the record"""
        return tuple(self._antismash_domains)

    def get_antismash_domains_by_tool(self, tool: str) -> Tuple[AntismashDomain, ...]:
        """A list of secondary metabolite aSDomains present in the record
           filtered by the given tool
        """
        return tuple(self._antismash_domains_by_tool.get(tool, []))

    def get_antismash_domains_in_cds(self, cds: Union[str, CDSFeature]) -> Tuple[AntismashDomain, ...]:
        """ A list of secondary metabolite aSDomains contained by a CDSFeature. Either the
            CDS name or the CDSFeature itself can be used to specify which CDS.
        """
        if isinstance(cds, CDSFeature):
            cds_name = cds.get_name()
        elif isinstance(cds, str):
            cds_name = cds
        else:
            raise TypeError(f"CDS must be a string or CDSFeature, not {type(cds)}")
        return tuple(self._antismash_domains_by_cds_name[cds_name])

    def clear_antismash_domains(self) -> None:
        "Remove all AntismashDomain features"
        # remove antismash domains only from the domains mapping
        for domain in self._antismash_domains:
            del self._domains_by_name[domain.get_name()]
        self._antismash_domains.clear()
        self._antismash_domains_by_tool.clear()
        self._antismash_domains_by_cds_name.clear()

    def get_generics(self) -> Tuple:
        """A list of secondary metabolite generics present in the record"""
        return tuple(self._nonspecific_features)

    def get_misc_feature_by_type(self, label: str) -> Tuple[Feature, ...]:
        """Returns a tuple of all generic features with a type matching label"""
        if label in ["protocluster", CandidateCluster.FEATURE_TYPE, "CDS", "CDSmotif", "subregion",
                     "region", "PFAM_domain", "aSDomain", "aSProdPred"]:
            raise ValueError(f"Use the appropriate get_* type instead for {label}")
        return tuple(i for i in self.get_generics() if i.type == label)

    @property
    def all_features(self) -> Iterable[Feature]:
        """ An iterator over all the features contained by the record """
        # in order of least references to other features
        return itertools.chain(
            self._sources,
            self._nonspecific_features,
            self._genes,
            self._cds_features,
            self._cds_motifs,
            self._antismash_domains,
            self._pfam_domains,
            self._modules,
            self._subregions,
            self._protoclusters,
            self._candidate_clusters,
            self._regions,
        )

    def get_cds_features_within_location(self, location: Location,
                                         with_overlapping: bool = False) -> List[CDSFeature]:
        """ Returns all CDS features within the given location

            Arguments:
                location: the location to use as a range
                with_overlapping: whether to include features which overlap the
                                  edges of the range

            Returns:
                a list of CDSFeatures, ordered by earliest position in feature location
        """
        results: list[CDSFeature] = []
        # shortcut if no CDS features exist
        if not self._cds_features:
            return results

        # compound locations need each chunk to be handled separately
        if len(location.parts) > 1:
            features: list[CDSFeature] = []
            # this is not particularly efficient, but it gets very complicated very quickly
            for part in location.parts:
                found = self.get_cds_features_within_location(part, with_overlapping=True)
                features.extend(f for f in found if f not in features)
            return [f for f in features if f.is_contained_by(location)]

        def find_start_in_list(location: Location, features: List[CDSFeature],
                               include_overlaps: bool) -> int:
            """ Find the earliest feature that starts before the location
                (and ends before, if include_overlaps is True)
            """
            dummy = Feature(location, feature_type='dummy')
            index = bisect.bisect_left(features, dummy)
            while index > 0 and features[index - 1].location.start == location.start:
                index -= 1
            if include_overlaps:
                while index >= 1 and features[index - 1].overlaps_with(dummy):
                    index -= 1
            return index

        if location.start < 0:
            assert isinstance(location, FeatureLocation)
            location = FeatureLocation(0, max(1, location.end))

        index = find_start_in_list(location, self._cds_features, with_overlapping)
        while index < len(self._cds_features):
            feature = self._cds_features[index]
            if feature.is_contained_by(location):
                results.append(feature)
            elif with_overlapping and feature.overlaps_with(location):
                results.append(feature)
            elif index + 1 < len(self._cds_features) and self._cds_features[index + 1].is_contained_by(feature):
                pass
            else:
                break
            index += 1
        return results

    def to_biopython(self) -> SeqRecord:
        """Returns a Bio.SeqRecord instance of the record"""
        bio_features: List[SeqFeature] = []
        for feature in sorted(self.all_features):
            bio_features.extend(feature.to_biopython())
        if "molecule_type" not in self._record.annotations:
            self._record.annotations["molecule_type"] = "DNA"
        return SeqRecord(self.seq, id=self._record.id, name=self._record.name,
                         description=self._record.description,
                         dbxrefs=self.dbxrefs, features=bio_features,
                         annotations=self._record.annotations,
                         letter_annotations=self._record.letter_annotations)

    def get_feature_count(self) -> int:
        """ Returns the total number of features contained in the record. """
        collections = [self._cds_features, self._protoclusters,
                       self._candidate_clusters, self._cds_motifs,
                       self._pfam_domains, self._antismash_domains,
                       self._nonspecific_features,
                       self._genes, self._regions, self._subregions,
                       self._sources]
        return sum(map(len, cast(List[List[Feature]], collections)))

    def add_gene(self, gene: Gene) -> None:
        """ Adds a Gene feature to the record """
        assert isinstance(gene, Gene), type(gene)
        ensure_valid_locations([gene], self.is_circular(), len(self))
        self._genes.append(gene)
        self._genes_by_name[gene.get_name()].append(gene)

    def add_cds_feature(self, cds_feature: CDSFeature) -> None:
        """ Add the given CDSFeature to the record,
            causes cluster-CDS pairing to be recalculated """
        assert isinstance(cds_feature, CDSFeature), type(cds_feature)
        # ensure it has a translation
        if not cds_feature.translation:
            raise ValueError(f"Missing translation for {cds_feature}")
        if len(cds_feature.translation) > 100000:
            raise ValueError(f"Translation too large: {cds_feature} {len(cds_feature.translation)}")
        location_key = str(cds_feature.location)
        if location_key in self._cds_by_location:
            raise SecmetInvalidInputError(
                f"Multiple CDS features have the same location: {cds_feature.location}")
        if cds_feature.get_name() in self._cds_by_name:
            error = SecmetInvalidInputError(
                f"multiple CDS features have the same name for mapping: {cds_feature.get_name()}"
            )
            # handle cases like splice variants
            if not cds_feature.locus_tag:
                raise error

            existing = self._cds_by_name[cds_feature.get_name()]
            if not (cds_feature.overlaps_with(existing) or
                    any(map(cds_feature.overlaps_with, self.get_genes_by_name(cds_feature.get_name())))):
                raise error

            new = f"{cds_feature.locus_tag}_{_location_checksum(cds_feature)}"
            cds_feature.locus_tag = new
            assert cds_feature.get_name() not in self._cds_by_name
        # only modify the record once all the checks are complete, otherwise
        # an exception can be caught leaving the state partially modified
        index = bisect.bisect_left(self._cds_features, cds_feature)
        self._cds_cache_dirty = True
        self._cds_features.insert(index, cds_feature)
        self._link_cds_to_parent(cds_feature)
        self._cds_by_location[location_key] = cds_feature
        self._cds_by_name[cds_feature.get_name()] = cds_feature

    def add_cds_motif(self, motif: Union[CDSMotif, Prepeptide]) -> None:
        """ Add the given CDSMotif to the record """
        assert isinstance(motif, (CDSMotif, Prepeptide)), f"{type(motif)}, {motif.type}"
        ensure_valid_locations([motif], self.is_circular(), len(self))
        self._cds_motifs.append(motif)
        assert motif.get_name(), f"motif {motif} has no identifiers"
        if motif.get_name() in self._domains_by_name:
            raise SecmetInvalidInputError(
                f"multiple Domain features have the same name for mapping: {motif.get_name()}"
            )
        self._domains_by_name[motif.get_name()] = motif
        if isinstance(motif, Prepeptide):
            assert motif.tool is not None

    def add_pfam_domain(self, pfam_domain: PFAMDomain) -> None:
        """ Add the given PFAMDomain to the record and links it in the parent CDS """
        assert isinstance(pfam_domain, PFAMDomain)
        ensure_valid_locations([pfam_domain], self.is_circular(), len(self))
        assert pfam_domain.get_name()
        self._pfam_domains.append(pfam_domain)
        if pfam_domain.get_name() in self._domains_by_name:
            raise SecmetInvalidInputError(
                f"multiple Domain features have the same name for mapping: {pfam_domain.get_name()}"
            )
        self._domains_by_name[pfam_domain.get_name()] = pfam_domain
        if pfam_domain.locus_tag:
            self._pfams_by_cds_name[pfam_domain.locus_tag].append(pfam_domain)

    def add_antismash_domain(self, antismash_domain: AntismashDomain) -> None:
        """ Add the given AntismashDomain to the record """
        assert isinstance(antismash_domain, AntismashDomain)
        ensure_valid_locations([antismash_domain], self.is_circular(), len(self))
        assert antismash_domain.get_name()
        assert antismash_domain.tool
        self._antismash_domains.append(antismash_domain)
        self._antismash_domains_by_tool[antismash_domain.tool].append(antismash_domain)
        if antismash_domain.get_name() in self._domains_by_name:
            raise SecmetInvalidInputError(
                f"multiple Domain features have the same name for mapping: {antismash_domain.get_name()}"
            )
        self._domains_by_name[antismash_domain.get_name()] = antismash_domain
        if antismash_domain.locus_tag:
            self._antismash_domains_by_cds_name[antismash_domain.locus_tag].append(antismash_domain)

    def add_module(self, module: Module) -> None:
        """ Add the given Module to the record """
        assert isinstance(module, Module)
        ensure_valid_locations([module], self.is_circular(), len(self))

        parents = set()
        for domain in module.domains:
            if domain.get_name() not in self._domains_by_name:
                raise ValueError(f"domain contained in module is not contained in record: {domain}")
            if not domain.locus_tag:
                raise ValueError(f"domain contained in module does not refer to a CDS: {domain}")
            parents.add(domain.locus_tag)

        for parent in parents:
            cds = self.get_cds_by_name(parent)
            if not cds:
                raise ValueError(f"domain contained in module refers to missing CDS: {parent}")
            cds.add_module(module)
        self._modules.append(module)

    def clear_modules(self) -> None:
        """ Removes all Modules from the record """
        self._modules.clear()

    def get_modules(self) -> Tuple[Module, ...]:
        """ Returns all Modules in the record """
        return tuple(self._modules)

    def add_source(self, source: Source) -> None:
        """ Add the given Source to the record """
        assert isinstance(source, Source)
        ensure_valid_locations([source], self.is_circular(), len(self))

        self._sources.append(source)

    def clear_sources(self) -> None:
        """ Removes all Sources from the record """
        self._sources.clear()

    def get_sources(self) -> Tuple[Source, ...]:
        """ Returns all Sources in the record """
        return tuple(self._sources)

    def has_multiple_sources(self) -> bool:
        """ Returns True if the record contains multiple Source features """
        return len(self._sources) > 1

    def add_feature(self, feature: Feature) -> None:
        """ Adds a Feature or any subclass to the relevant list """
        assert isinstance(feature, Feature), type(feature)
        ensure_valid_locations([feature], self.is_circular(), len(self))
        if isinstance(feature, Protocluster):
            self.add_protocluster(feature)
        elif isinstance(feature, CandidateCluster):
            self.add_candidate_cluster(feature)
        elif isinstance(feature, SubRegion):
            self.add_subregion(feature)
        elif isinstance(feature, Region):
            self.add_region(feature)
        elif isinstance(feature, CDSFeature):
            self.add_cds_feature(feature)
        elif isinstance(feature, Gene):
            self.add_gene(feature)
        elif isinstance(feature, (CDSMotif, Prepeptide)):
            self.add_cds_motif(feature)
        elif isinstance(feature, PFAMDomain):
            self.add_pfam_domain(feature)
        elif isinstance(feature, AntismashDomain):
            self.add_antismash_domain(feature)
        elif isinstance(feature, Module):
            self.add_module(feature)
        elif isinstance(feature, Source):
            self.add_source(feature)
        else:
            self._nonspecific_features.append(feature)

    def add_biopython_feature(self, feature: SeqFeature) -> None:
        """ Convert a biopython feature to SecMet feature, then add it to the
            record.
        """
        if feature.type == CDSFeature.FEATURE_TYPE:
            self.add_cds_feature(CDSFeature.from_biopython(feature, record=self))
        elif feature.type == Gene.FEATURE_TYPE:
            self.add_gene(Gene.from_biopython(feature, record=self))
        elif feature.type == Protocluster.FEATURE_TYPE:
            self.add_protocluster(Protocluster.from_biopython(feature, record=self))
        elif feature.type == "proto_core":
            # discard this, as info contained in it is in "protocluster" features
            pass
        elif feature.type == CDSMotif.FEATURE_TYPE:
            # skip component parts of prepeptides and regenerate from the core
            prepeptide = feature.qualifiers.get("prepeptide", [""])[0]
            if prepeptide:
                if prepeptide != "core":
                    return
                self.add_cds_motif(Prepeptide.from_biopython(feature, record=self))
                return
            motif = CDSMotif.from_biopython(feature, record=self)
            if not motif.domain_id and not motif.created_by_antismash:
                counter = 1
                template = f"non_aS_motif_{motif.location.start}_{motif.location.end}_{{}}"
                while template.format(counter) in self._domains_by_name:
                    counter += 1
                motif.domain_id = template.format(counter)
            self.add_cds_motif(motif)
        elif feature.type == PFAMDomain.FEATURE_TYPE:
            self.add_pfam_domain(PFAMDomain.from_biopython(feature, record=self))
        elif feature.type == AntismashDomain.FEATURE_TYPE:
            self.add_antismash_domain(AntismashDomain.from_biopython(feature, record=self))
        elif feature.type == CandidateCluster.FEATURE_TYPE:
            self.add_candidate_cluster(CandidateCluster.from_biopython(feature, record=self))
        elif feature.type == Region.FEATURE_TYPE:
            self.add_region(Region.from_biopython(feature, record=self))
        elif feature.type == SubRegion.FEATURE_TYPE:
            self.add_subregion(SubRegion.from_biopython(feature, record=self))
        elif feature.type == Module.FEATURE_TYPE:
            self.add_module(Module.from_biopython(feature, record=self))
        elif feature.type == Source.FEATURE_TYPE:
            self.add_source(Source.from_biopython(feature, record=self))
        else:
            self.add_feature(Feature.from_biopython(feature))

    @classmethod
    def from_biopython(cls: Type[T], seq_record: SeqRecord, taxon: str,
                       discard_antismash_features: bool = False, **kwargs: Any) -> T:
        """ Constructs a new Record instance from a biopython SeqRecord,
            also replaces biopython SeqFeatures with Feature subclasses

            If 'discard_antismash_features' is True, then feature types that are
            specific to antiSMASH that fail to convert will instead become generic
            'Feature'
        """
        postponed_features: Dict[str, Tuple[Type[Feature], List[SeqFeature]]] = OrderedDict()
        for kind in [CandidateCluster, Region, Module]:  # type: Type[Feature]
            postponed_features[kind.FEATURE_TYPE] = (kind, [])

        assert isinstance(seq_record, SeqRecord), type(seq_record)
        molecule_type = seq_record.annotations.get("molecule_type", "DNA")
        if not molecule_type.upper().endswith("DNA"):
            raise SecmetInvalidInputError(f"{molecule_type} records are not supported")
        if seq_record.seq and not Record.is_nucleotide_sequence(seq_record.seq):
            raise SecmetInvalidInputError("protein records are not supported")
        transl_table = 1  # standard
        if str(taxon) == "bacteria":
            transl_table = 11  # bacterial, archea, plant plastid code
        record = cls(seq=seq_record.seq, transl_table=transl_table, **kwargs)
        record._record = seq_record
        # because is_circular() can't be used reliably at this stage due to fasta files
        can_be_circular = taxon == "bacteria"
        try:
            ensure_valid_locations(seq_record.features, can_be_circular, len(seq_record.seq))
        except ValueError as err:
            raise SecmetInvalidInputError(f"{seq_record.id}: {err}") from err

        for feature in seq_record.features:
            if feature.location.ref or feature.location.ref_db:
                for ref in [feature.location.ref, feature.location.ref_db]:
                    if ref and ref != seq_record.id:
                        raise SecmetInvalidInputError(f"feature references another sequence: {feature.ref}")
                # to handle a biopython issue, set the references to None
                feature.ref = None
                feature.ref_db = None

            # prefilter some NCBI Pfam hits locations that are generated poorly
            if all([can_be_circular, feature.type == "misc_feature",
                    location_bridges_origin(feature.location, allow_reversing=False)]):
                feature.location = remove_redundant_exons(feature.location)

            if feature.type in postponed_features:
                # again, since CDSs need translations, skip if too small or not a CDS
                if feature.type != "CDS" or len(feature) >= 3:
                    postponed_features[feature.type][1].append(feature)
                continue
            # again, since CDSs need translations, skip if too small or not a CDS
            if feature.type != "CDS" or len(feature) >= 3:
                try:
                    record.add_biopython_feature(feature)
                except ValueError as err:
                    if feature.type not in ANTISMASH_SPECIFIC_TYPES or not discard_antismash_features:
                        raise SecmetInvalidInputError(str(err)) from err

        for feature_class, features in postponed_features.values():
            for feature in features:
                try:
                    record.add_feature(feature_class.from_biopython(feature, record=record))
                except ValueError as err:
                    # all postponed features are antismash-specific
                    assert feature.type in ANTISMASH_SPECIFIC_TYPES, feature.type
                    if not discard_antismash_features:
                        raise SecmetInvalidInputError(str(err)) from err
        return record

    @staticmethod
    def from_genbank(filepath: str, taxon: str = "bacteria") -> List["Record"]:
        """ Reads a genbank file and creates a Record instance for each record
            contained in the file.
        """
        records = []
        for bio in SeqIO.parse(filepath, "genbank"):
            records.append(Record.from_biopython(bio, taxon))
        return records

    def _link_cds_to_parent(self, cds: CDSFeature) -> None:
        """ connect the given CDS to any collection that contains it, if any """
        assert isinstance(cds, CDSFeature)

        # since regions cannot overlap, bisection is possible
        left = bisect.bisect_left(self._regions, cds)
        right = bisect.bisect_right(self._regions, cds, lo=left)
        for region in self._regions[left - 1:right + 1]:
            if cds.is_contained_by(region):
                region.add_cds(cds)
                cds.region = region

        # for other collections, since they may overlap heavily, exhaustive search required
        other_collections: Sequence[Sequence[CDSCollection]] = [
            self._protoclusters,
            self._candidate_clusters,
            self._subregions
        ]
        for collections in other_collections:
            for collection in collections:
                if cds.is_contained_by(collection):
                    collection.add_cds(cds)

    def get_aa_translation_from_location(self, location: FeatureLocation,
                                         transl_table: Union[str, int] = None) -> Seq:
        """ Obtain the translation for a feature based on its location """
        if location.end > len(self.seq):
            raise ValueError("location outside available sequence")
        if transl_table is None:
            transl_table = self._transl_table
        extracted = location.extract(self.seq).replace("-", "")
        if len(extracted) % 3 != 0:
            extracted = extracted[:-(len(extracted) % 3)]
        seq = extracted.translate(to_stop=True, table=transl_table)
        if not seq:
            # go past stop codons and hope for something to work with
            seq = extracted.translate(table=transl_table)

        # replace ambiguous proteins with an explicit unknown
        string_version = str(seq)
        for invalid in "*BJOUZ":
            string_version = string_version.replace(invalid, "X")

        if "-" in str(seq):
            seq = Seq(str(seq).replace("-", ""))

        return Seq(string_version)

    def get_cds_features_within_regions(self) -> list[CDSFeature]:
        """ Returns all CDS features in the record that are located within a
            region of interest
        """
        features: List[CDSFeature] = []
        for region in self._regions:
            features.extend(region.cds_children)
        return features

    def connect_locations(self, locations: list[Location], *, disable_wrapping: bool = False) -> Location:
        """ Combines the given locations into a single contiguous location, unless the record
            crosses the origin, in which case the resulting location may be split over the origin/

            Arguments:
                locations: the locations to combine
                disable_wrapping: if provided, explicitly disables forming a location over the origin

            Returns:
                the single location that covers all the input locations
        """
        wrap_point = len(self) if self.is_circular() and not disable_wrapping else None
        return connect_locations(locations, wrap_point=wrap_point)

    def create_candidate_clusters(self) -> int:
        """ Takes all Cluster instances and constructs CandidateClusters that cover
            each Cluster. Each combination of overlapping clusters will create
            a CandidateCluster, including a single cluster.

            Returns:
                the number of candidate_clusters created
        """
        if not self._protoclusters:
            return 0

        wrap_point = len(self) if self.is_circular() else None
        candidate_clusters = create_candidates_from_protoclusters(self._protoclusters, circular_wrap_point=wrap_point)

        for candidate_cluster in sorted(candidate_clusters):
            self.add_candidate_cluster(candidate_cluster)

        return len(candidate_clusters)

    def create_regions(self, candidate_clusters: List[CandidateCluster] = None,
                       subregions: List[SubRegion] = None) -> int:
        """ Creates Region features based on contained CandidateClusters and SubRegions
            and returns the number of regions created. Regions will not overlap.

            If supplied, parameters will override the Records own candidate_clusters
            and subregions.
        """
        if candidate_clusters is None:
            candidate_clusters = self._candidate_clusters
        if subregions is None:
            subregions = self._subregions

        if not candidate_clusters and not subregions:
            return 0

        areas: List[CDSCollection] = []
        areas.extend(candidate_clusters)
        areas.extend(subregions)
        areas.sort()

        wrap_point = len(self) if self.is_circular() else None

        # find all overlapping sets
        sections: list[tuple[Location, list[CDSCollection]]] = []

        location: Location = areas[0].location
        included_areas = [areas[0]]

        for area in areas[1:]:
            if not area.overlaps_with(location):
                sections.append((location, included_areas))
                location = area.location
                included_areas = [area]
            else:
                location = connect_locations([area.location, location], wrap_point=wrap_point)
                included_areas.append(area)

        # finalise the last, unterminated section
        sections.append((location, included_areas))

        # then handle any cases over cross-origin overlap of sections by merging first and last
        if len(sections) > 1:
            first_location, first_areas = sections[0]
            last_location, last_areas = sections[-1]
            if locations_overlap(first_location, last_location):
                sections.pop()
                location = connect_locations([first_location, last_location], wrap_point=wrap_point)
                for area in last_areas:
                    if area not in first_areas:
                        first_areas.append(area)
                sections[0] = (location, first_areas)

        # finally, create and add a region for each section
        regions_added = len(sections)
        for _, areas in sections:
            candidates = []
            subs = []
            for area in areas:
                if isinstance(area, CandidateCluster):
                    candidates.append(area)
                else:
                    assert isinstance(area, SubRegion), type(area)
                    subs.append(area)
            self.add_region(Region(candidates, subs))

        return regions_added

    def extend_location(self, location: Location, distance: int) -> Location:
        """ Constructs a new location which covers the given distance from the
            given location, capped at record limits unless the record is circular,
            in which case it wraps around.

            Arguments:
                location: the location to create an extended version of
                distance: the distance to extend the location, applies in both directions

            Returns:
                a new location, covering the requested area(s)
        """

        parts = list(location.parts)  # the original location must not be changed
        maximum = len(self)
        if location.strand == -1:
            parts.reverse()

        # catch any case where both sides extend past the edges *and then overlap*
        # in which case a simple bounding to the record edges is enough
        n0 = parts[0].start
        n1 = parts[-1].end
        ns = n0 - distance
        ne = n1 + distance
        if self.is_circular() and ns < 0 and ns + maximum <= ne:
            parts[0] = FeatureLocation(0, parts[0].end, location.strand)
            parts[-1] = FeatureLocation(parts[-1].start, maximum, location.strand)
            # merge any parts that now overlap with each other
            if ns < 0:
                upper = FeatureLocation(ns + maximum, maximum, location.strand)
                merged = False
                while parts and locations_overlap(parts[-1], upper):
                    merged = True
                    upper = FeatureLocation(min(parts[-1].start, upper.start), maximum, location.strand)
                    parts.pop()
                if merged:
                    parts = [upper] + parts
            if ne > maximum:
                lower = FeatureLocation(0, ne % maximum, location.strand)
                merged = False
                while parts and locations_overlap(parts[0], lower):
                    merged = True
                    lower = FeatureLocation(0, max(parts[0].end, lower.end), location.strand)
                    parts.pop(0)
                if merged:
                    parts.append(lower)

            # if there's only one part remaining, it should be reduced to a whole-record location
            if len(parts) == 1:
                assert parts == [FeatureLocation(0, maximum, location.strand)], parts
                return parts[0]

            if location.strand == -1:
                parts.reverse()
            if len(location.parts) > 1:
                return CompoundLocation(parts, operator=location.operator)
            return CompoundLocation(parts)

        # if the wrap around goes so far as to overlap other parts, merge them
        while len(parts) > 1 and locations_overlap(parts[0], parts[-1]):
            first = parts[0]
            second = parts[-1]
            parts[0] = FeatureLocation(min(first.start, second.start), max(first.end, second.end),
                                       first.strand if first.strand == second.strand else 0)
            parts.pop()

        start_part = parts[0]
        if start_part.start - distance < 0 and self.is_circular():
            parts[0] = FeatureLocation(0, start_part.end, parts[0].strand)
            parts.insert(0, FeatureLocation(min(maximum + (start_part.start - distance), maximum),
                                            maximum, start_part.strand))
        else:
            parts[0] = FeatureLocation(max(0, start_part.start - distance), start_part.end, start_part.strand)

        end_part = parts[-1]
        if end_part.end + distance > maximum and self.is_circular():
            parts[-1] = FeatureLocation(end_part.start, maximum, end_part.strand)
            parts.append(FeatureLocation(0, min(end_part.end + distance - maximum, maximum), end_part.strand))
        else:
            parts[-1] = FeatureLocation(end_part.start, min(end_part.end + distance, maximum), end_part.strand)

        # if the wrap around goes so far as to overlap other parts, merge them
        while len(parts) > 1 and locations_overlap(parts[0], parts[-1]):
            first = parts[0]
            second = parts[-1]
            parts[0] = FeatureLocation(min(first.start, second.start), max(first.end, second.end),
                                       first.strand if first.strand == second.strand else 0)
            parts.pop()

        if len(parts) == 1:
            return parts[0]
        if location.strand == -1:
            parts.reverse()
        # reuse the operator of the original, if it was compound to start with
        if len(location.parts) > 1:
            return CompoundLocation(parts, operator=location.operator)
        # otherwise use the default
        return CompoundLocation(parts)

    def get_distance_between_features(self, first: Feature, second: Feature) -> int:
        """ Returns the shortest distance between the two given features, crossing
            the origin if the record is circular.

            Overlapping features are considered to have zero distance.
        """
        return self.get_distance_between_locations(first.location, second.location)

    def get_distance_between_locations(self, first: Location, second: Location) -> int:
        """ Returns the shortest distance between the two given locations, crossing
            the origin if the record is circular.

            Overlapping locations are considered to have zero distance.
        """
        if self.is_circular():
            return get_distance_between_locations(first, second, wrap_point=len(self))
        return get_distance_between_locations(first, second)

    def get_nrps_pks_cds_features(self) -> List[CDSFeature]:
        """ Returns a list of all CDS features within Clusters that contain at least
            one NRPS/PKS domain.
        """
        return [feature for feature in self.get_cds_features_within_regions() if feature.nrps_pks.domains]

    def get_gc_content(self) -> float:
        """ Calculate the GC content of the record's sequence """
        if not self.seq:
            raise ValueError("Cannot calculate GC content of empty sequence")
        if self._gc_content < 0:  # not cached
            counter = Counter(str(self.seq))
            gc_count = counter['G'] + counter['C'] + counter['g'] + counter['c']
            self._gc_content = gc_count / len(self)
        return self._gc_content

    def strip_antismash_annotations(self) -> None:
        """ Removes all antismash features and annotations from the record """
        self.clear_protoclusters()
        self.clear_candidate_clusters()
        self.clear_subregions()
        self.clear_regions()
        self.clear_antismash_domains()
        self.clear_pfam_domains()
        self.clear_modules()

        # clean up antiSMASH-created CDSMotifs, but leave the rest
        motifs = self.get_cds_motifs()
        self.clear_cds_motifs()
        for motif in motifs:
            if not motif.created_by_antismash:
                self.add_cds_motif(motif)

        # clean up antiSMASH annotations in CDS features
        for feature in self.get_cds_features():
            feature.strip_antismash_annotations()

    @staticmethod
    def is_nucleotide_sequence(sequence: Union[Seq, str]) -> bool:
        """ Determines if a sequence is a nucleotide sequence based on content.

            Arguments:
                sequence: the sequence to check, either a string or Bio.Seq

            Returns:
                True if more than 80% of characters are nucleotide bases
        """
        try:
            other = str(sequence).lower()
        except UndefinedSequenceError:
            return False
        for char in "acgtn":
            other = other.replace(char, "")
        return len(other) < 0.2 * len(sequence)

    def to_genbank(self, filename: str) -> None:
        """ Writes the record to the given path in GenBank format

            Arguments:
                filename: the file path to write to
        """
        SeqIO.write([self.to_biopython()], filename, "genbank")


def _calculate_crc32(string: str) -> str:
    """ Calculates the crc32 checksum of an input string and returns the resulting checksum in hex.

        Arguments:
            string: The string to generate the checksum for

        Returns:
            A string containing the hexadecimal representation of the crc32 checksum
    """
    checksum = crc32(string.encode("utf-8"))
    return f"{checksum:x}"


def _location_checksum(feature: Feature) -> str:
    """ Calculate the checksum of the given feature's location and return the resulting value in a hex string.

        Arguments:
            feature: A Feature object to generate the location checksum for

        Returns:
            A string containing the feature's location checksum in hex
    """
    return _calculate_crc32(str(feature.location))
