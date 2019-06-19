# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Implements a more readable and efficient version of biopython's SeqRecord.

    Most attributes are wrapped, with the record features being secmet features
    and accessible and already paritioned so accessing specific types of feature
    is more efficient.

    Does not, and will not, subclass SeqRecord to avoid mixing code types.
"""


import bisect
from collections import Counter, defaultdict
import logging
from typing import Any, Dict, List, Tuple, Union, cast
from typing import Optional, Sequence, Set  # comment hints # pylint: disable=unused-import

from Bio import Alphabet, SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from Bio.SeqRecord import SeqRecord

from .errors import SecmetInvalidInputError
from .features import (
    AntismashDomain,
    CandidateCluster,
    CDSFeature,
    CDSMotif,
    Domain,
    Feature,
    Gene,
    PFAMDomain,
    Prepeptide,
    Protocluster,
    Region,
    SubRegion,
)
from .features import CDSCollection  # comment hints, pylint: disable=unused-import
from .features.candidate_cluster import create_candidates_from_protoclusters

from .locations import (
    location_bridges_origin,
    split_origin_bridging_location,
    combine_locations,
    ensure_valid_locations,
)


class Record:
    """A record containing secondary metabolite clusters"""
    # slots not for space, but to stop use as a horrible global
    __slots__ = ["_record", "_seq", "skip", "_cds_features", "_cds_by_name",
                 "_protoclusters", "original_id", "_cds_motifs", "_pfam_domains",
                 "_antismash_domains", "_protocluster_numbering", "_nonspecific_features",
                 "record_index", "_genes", "_transl_table", "_domains_by_name",
                 "_pfams_by_cds_name",
                 "_candidate_clusters", "_candidate_clusters_numbering",
                 "_subregions", "_subregion_numbering",
                 "_regions", "_region_numbering"]

    def __init__(self, seq: Union[Seq, str] = "", transl_table: int = 1, **kwargs: Any) -> None:
        # prevent paths from being used as a sequence
        assert not set("./\\").issubset(set(seq)), "Invalid sequence provided"
        if isinstance(seq, str):
            seq = Seq(seq)
        self._record = SeqRecord(seq, **kwargs)
        self.record_index = None  # type: Optional[int]
        self.original_id = None
        self.skip = None  # type: Optional[str] # TODO: move to yet another abstraction layer?
        self._genes = []  # type: List[Gene]

        self._cds_features = []  # type: List[CDSFeature]
        self._cds_by_name = {}  # type: Dict[str, CDSFeature]

        self._cds_motifs = []  # type: List[CDSMotif]

        self._pfam_domains = []  # type: List[PFAMDomain]
        self._pfams_by_cds_name = defaultdict(list)  # type: Dict[str, List[PFAMDomain]]

        self._antismash_domains = []  # type: List[AntismashDomain]

        # includes PFAMDomains and AntismashDomains
        self._domains_by_name = {}  # type: Dict[str, Domain]  # for use as x[domain.get_name()] = domain

        self._nonspecific_features = []  # type: List[Feature]

        self._protoclusters = []  # type: List[Protocluster]
        self._protocluster_numbering = {}  # type: Dict[Protocluster, int]

        self._candidate_clusters = []  # type: List[CandidateCluster]
        self._candidate_clusters_numbering = {}  # type: Dict[CandidateCluster, int]

        self._subregions = []  # type: List[SubRegion]
        self._subregion_numbering = {}  # type: Dict[SubRegion, int]

        self._regions = []  # type: List[Region]
        self._region_numbering = {}  # type: Dict[Region, int]

        self._transl_table = int(transl_table)

    def __getattr__(self, attr: str) -> Any:
        # passthroughs to the original SeqRecord
        if attr in ["id", "seq", "description", "name", "annotations", "dbxrefs"]:
            return getattr(self._record, attr)
        if attr in Record.__slots__:
            return self.__getattribute__(attr)
        raise AttributeError("Record has no attribute '%s'" % attr)

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

    # protocluster manipulation
    def add_protocluster(self, cluster: Protocluster) -> None:
        """ Add the given cluster to the record,
            causes cluster-CDS pairing to be recalculated """
        assert isinstance(cluster, Protocluster), type(cluster)
        assert cluster.location.start >= 0, cluster
        assert cluster.location.end <= len(self), "%s > %d" % (cluster, len(self))
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
        assert cluster.location.end <= len(self), "%s > %d" % (cluster, len(self))
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
        assert subregion.location.end <= len(self), "%s > %d" % (subregion, len(self))
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
        assert region.location.end <= len(self), "%s > %d" % (region, len(self))
        index = 0
        for i, existing_region in enumerate(self._regions):  # TODO: fix performance
            if region.overlaps_with(existing_region):
                logging.error("existing %s overlaps with new %s", existing_region, region)
                raise ValueError("regions cannot overlap")
            if region < existing_region:
                index = i  # before
                break
            else:
                index = i + 1  # after
        self._regions.insert(index, region)
        region.parent_record = self
        # update numbering
        for i in range(index, len(self._regions)):
            self._region_numbering[self._regions[i]] = i + 1  # 1-indexed
        # link any relevant CDS features
        self._link_region_to_cds_features(region)

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

    def get_cds_features(self) -> Tuple[CDSFeature, ...]:
        """A list of secondary metabolite clusters present in the record"""
        return tuple(self._cds_features)

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
            raise KeyError("record %s contains no domain named %s" % (self.id, name))

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
            raise TypeError("CDS must be a string or CDSFeature, not %s" % type(cds))
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

    def clear_antismash_domains(self) -> None:
        "Remove all AntismashDomain features"
        # remove antismash domains only from the domains mapping
        for domain in self._antismash_domains:
            del self._domains_by_name[domain.get_name()]
        self._antismash_domains.clear()

    def get_generics(self) -> Tuple:
        """A list of secondary metabolite generics present in the record"""
        return tuple(self._nonspecific_features)

    def get_misc_feature_by_type(self, label: str) -> Tuple[Feature, ...]:
        """Returns a tuple of all generic features with a type matching label"""
        if label in ["protocluster", CandidateCluster.FEATURE_TYPE, "CDS", "CDSmotif", "subregion",
                     "region", "PFAM_domain", "aSDomain", "aSProdPred"]:
            raise ValueError("Use the appropriate get_* type instead for %s" % label)
        return tuple(i for i in self.get_generics() if i.type == label)

    def get_all_features(self) -> List[Feature]:
        """ Returns all features
            note: This is slow, if only a specific type is required, use
                  the other get_*() functions
        """
        features = list(self.get_generics())
        features.extend(self._genes)
        features.extend(self.get_protoclusters())
        features.extend(self.get_candidate_clusters())
        features.extend(self.get_subregions())
        features.extend(self.get_regions())
        features.extend(self.get_cds_features())
        features.extend(self.get_cds_motifs())
        features.extend(self.get_antismash_domains())
        features.extend(self.get_pfam_domains())
        return features

    def get_cds_features_within_location(self, location: FeatureLocation,
                                         with_overlapping: bool = False) -> List[CDSFeature]:
        """ Returns all CDS features within the given location

            Arguments:
                location: the location to use as a range
                with_overlapping: whether to include features which overlap the
                                  edges of the range

            Returns:
                a list of CDSFeatures, ordered by earliest position in feature location
        """
        def find_start_in_list(location: FeatureLocation, features: List[CDSFeature],
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

        results = []  # type: List[CDSFeature]
        # shortcut if no CDS features exist
        if not self._cds_features:
            return results
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
        features = self.get_all_features()
        bio_features = []  # type: List[SeqFeature]
        for feature in sorted(features):
            bio_features.extend(feature.to_biopython())
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
                       self._genes, self._regions, self._subregions]
        return sum(map(len, cast(List[List[Feature]], collections)))

    def add_gene(self, gene: Gene) -> None:
        """ Adds a Gene feature to the record """
        assert isinstance(gene, Gene), type(gene)
        self._genes.append(gene)

    def add_cds_feature(self, cds_feature: CDSFeature) -> None:
        """ Add the given CDSFeature to the record,
            causes cluster-CDS pairing to be recalculated """
        assert isinstance(cds_feature, CDSFeature), type(cds_feature)
        # ensure it has a translation
        if not cds_feature.translation:
            raise ValueError("Missing translation info for %s" % cds_feature)
        index = bisect.bisect_left(self._cds_features, cds_feature)
        self._cds_features.insert(index, cds_feature)
        self._link_cds_to_parent(cds_feature)
        if cds_feature.get_name() in self._cds_by_name:
            error = SecmetInvalidInputError("multiple CDS features have the same name for mapping: %s" %
                                            cds_feature.get_name())
            # handle cases like splice variants
            if not cds_feature.locus_tag or not cds_feature.protein_id:
                raise error
            if not cds_feature.overlaps_with(self._cds_by_name[cds_feature.get_name()]):
                raise error
            new = "%s_%s" % (cds_feature.locus_tag, cds_feature.protein_id)
            if new in self._cds_by_name:
                raise error
            cds_feature.locus_tag = new
            assert cds_feature.get_name() not in self._cds_by_name
        self._cds_by_name[cds_feature.get_name()] = cds_feature

    def add_cds_motif(self, motif: Union[CDSMotif, Prepeptide]) -> None:
        """ Add the given CDSMotif to the record """
        assert isinstance(motif, (CDSMotif, Prepeptide)), "%s, %s" % (type(motif), motif.type)
        self._cds_motifs.append(motif)
        assert motif.get_name(), "motif %s has no identifiers" % motif
        if motif.get_name() in self._domains_by_name:
            raise SecmetInvalidInputError("multiple Domain features have the same name for mapping: %s" %
                                          motif.get_name())
        self._domains_by_name[motif.get_name()] = motif
        if isinstance(motif, Prepeptide):
            assert motif.tool is not None

    def add_pfam_domain(self, pfam_domain: PFAMDomain) -> None:
        """ Add the given PFAMDomain to the record and links it in the parent CDS """
        assert isinstance(pfam_domain, PFAMDomain)
        assert pfam_domain.get_name()
        self._pfam_domains.append(pfam_domain)
        if pfam_domain.get_name() in self._domains_by_name:
            raise SecmetInvalidInputError("multiple Domain features have the same name for mapping: %s" %
                                          pfam_domain.get_name())
        self._domains_by_name[pfam_domain.get_name()] = pfam_domain
        if pfam_domain.locus_tag:
            self._pfams_by_cds_name[pfam_domain.locus_tag].append(pfam_domain)

    def add_antismash_domain(self, antismash_domain: AntismashDomain) -> None:
        """ Add the given AntismashDomain to the record """
        assert isinstance(antismash_domain, AntismashDomain)
        assert antismash_domain.get_name()
        self._antismash_domains.append(antismash_domain)
        if antismash_domain.get_name() in self._domains_by_name:
            raise SecmetInvalidInputError("multiple Domain features have the same name for mapping: %s" %
                                          antismash_domain.get_name())
        self._domains_by_name[antismash_domain.get_name()] = antismash_domain

    def add_feature(self, feature: Feature) -> None:
        """ Adds a Feature or any subclass to the relevant list """
        assert isinstance(feature, Feature), type(feature)
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
        else:
            self._nonspecific_features.append(feature)

    def add_biopython_feature(self, feature: SeqFeature) -> None:
        """ Convert a biopython feature to SecMet feature, then add it to the
            record.
        """
        if feature.type == 'CDS':
            self.add_cds_feature(CDSFeature.from_biopython(feature, record=self))
        elif feature.type == 'gene':
            self.add_gene(Gene.from_biopython(feature))
        elif feature.type == "protocluster":
            self.add_protocluster(Protocluster.from_biopython(feature))
        elif feature.type == "proto_core":
            # discard this, as info contained in it is in "protocluster" features
            pass
        elif feature.type == 'CDS_motif':
            # skip component parts of prepeptides and regenerate from the core
            prepeptide = feature.qualifiers.get("prepeptide", [""])[0]
            if prepeptide:
                if prepeptide != "core":
                    return
                self.add_cds_motif(Prepeptide.from_biopython(feature))
                return
            motif = CDSMotif.from_biopython(feature)
            if not motif.domain_id and not motif.created_by_antismash:
                counter = 1
                template = "non_aS_motif_%d_%d_{}" % (motif.location.start, motif.location.end)
                while template.format(counter) in self._domains_by_name:
                    counter += 1
                motif.domain_id = template.format(counter)
            self.add_cds_motif(motif)
        elif feature.type == 'PFAM_domain':
            self.add_pfam_domain(PFAMDomain.from_biopython(feature))
        elif feature.type == 'aSDomain':
            self.add_antismash_domain(AntismashDomain.from_biopython(feature))
        elif feature.type == CandidateCluster.FEATURE_TYPE:
            raise ValueError("CandidateCluster features cannot be directly added from biopython")
        elif feature.type == 'region':
            raise ValueError("Region features cannot be directly added from biopython")
        elif feature.type == 'subregion':
            self.add_subregion(SubRegion.from_biopython(feature))
        else:
            self.add_feature(Feature.from_biopython(feature))

    @staticmethod
    def from_biopython(seq_record: SeqRecord, taxon: str) -> "Record":
        """ Constructs a new Record instance from a biopython SeqRecord,
            also replaces biopython SeqFeatures with Feature subclasses
        """
        postponed_features = {
            "region": [],
            CandidateCluster.FEATURE_TYPE: [],
        }  # type: Dict[str, SeqFeature]

        assert isinstance(seq_record, SeqRecord), type(seq_record)
        if seq_record.seq and isinstance(seq_record.seq, Seq):
            if isinstance(seq_record.seq.alphabet, Alphabet.ProteinAlphabet):
                raise SecmetInvalidInputError("protein records are not supported")
        transl_table = 1  # standard
        if str(taxon) == "bacteria":
            transl_table = 11  # bacterial, archea, plant plastid code
        record = Record(transl_table=transl_table)
        record._record = seq_record  # pylint: disable=protected-access
        # because is_circular() can't be used reliably at this stage due to fasta files
        can_be_circular = taxon == "bacteria"
        try:
            ensure_valid_locations(seq_record.features, can_be_circular, len(seq_record.seq))
        except ValueError as err:
            raise SecmetInvalidInputError(str(err))

        for feature in seq_record.features:
            if feature.ref or feature.ref_db:
                for ref in [feature.ref, feature.ref_db]:
                    if ref and ref != seq_record.id:
                        raise SecmetInvalidInputError("feature references another sequence: (%s)" % feature.ref)
                # to handle a biopython issue, set the references to None
                feature.ref = None
                feature.ref_db = None

            locations_adjusted = False
            name_modified = False

            if can_be_circular and location_bridges_origin(feature.location, allow_reversing=False):
                locations_adjusted = True
                original_location = feature.location
                try:
                    lower, upper = split_origin_bridging_location(feature.location)
                except ValueError as err:
                    raise SecmetInvalidInputError(str(err)) from err

                if feature.type in ['CDS', 'gene']:
                    name_modified = True
                    original_gene_name = feature.qualifiers.get("gene", [None])[0]
                    gene_name = original_gene_name
                    locus_tag = feature.qualifiers.get("locus_tag", [None])[0]
                    # if neither exist, set the gene name as it is less precise
                    # in meaning
                    if not gene_name and not locus_tag:
                        gene_name = "bridge"

                # nuke any translation, since it's now out of date
                feature.qualifiers.pop('translation', None)

                # add a separate feature for the upper section
                if len(upper) > 1:
                    feature.location = CompoundLocation(upper, original_location.operator)
                else:
                    feature.location = upper[0]
                if name_modified:
                    if gene_name:
                        feature.qualifiers["gene"] = [gene_name + "_UPPER"]
                    if locus_tag:
                        feature.qualifiers["locus_tag"] = [locus_tag + "_UPPER"]

                # since CDSs need translations, skip if too small or not a CDS
                if feature.type != "CDS" or len(feature) >= 3:
                    try:
                        record.add_biopython_feature(feature)
                    except ValueError as err:
                        raise SecmetInvalidInputError(str(err)) from err

                # adjust the current feature to only be the lower section
                if len(lower) > 1:
                    feature.location = CompoundLocation(lower, original_location.operator)
                else:
                    feature.location = lower[0]
                if name_modified:
                    if gene_name:
                        feature.qualifiers["gene"] = [gene_name + "_LOWER"]
                    if locus_tag:
                        feature.qualifiers["locus_tag"] = [locus_tag + "_LOWER"]

            if feature.type in postponed_features:
                # again, since CDSs need translations, skip if too small or not a CDS
                if feature.type != "CDS" or len(feature) >= 3:
                    postponed_features[feature.type].append(feature)
                continue
            # again, since CDSs need translations, skip if too small or not a CDS
            if feature.type != "CDS" or len(feature) >= 3:
                try:
                    record.add_biopython_feature(feature)
                except ValueError as err:
                    raise SecmetInvalidInputError(str(err)) from err

            # reset back to how the feature looked originally
            if locations_adjusted:
                feature.location = original_location
                if name_modified:
                    if not locus_tag:
                        feature.qualifiers.pop("locus_tag", "")
                    else:
                        feature.qualifiers["locus_tag"][0] = locus_tag
                    if not original_gene_name:
                        feature.qualifiers.pop("gene", "")
                    else:
                        feature.qualifiers["gene"][0] = original_gene_name
        try:
            for feature in postponed_features[CandidateCluster.FEATURE_TYPE]:
                record.add_feature(CandidateCluster.from_biopython(feature).convert_to_real_feature(record))
            for feature in postponed_features["region"]:
                record.add_feature(Region.from_biopython(feature).convert_to_real_feature(record))
        except ValueError as err:
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
        other_collections = [self._protoclusters, self._candidate_clusters,
                             self._subregions]  # type: Sequence[Sequence[CDSCollection]]
        for collections in other_collections:
            for collection in collections:
                if cds.is_contained_by(collection):
                    collection.add_cds(cds)

    def _link_region_to_cds_features(self, region: Region) -> None:
        """ connect the given cluster to every CDS feature within it's range """
        assert isinstance(region, Region)
        # quickly find the first cds with equal start
        index = bisect.bisect_left(self._cds_features, region)
        # move backwards until we find one that doesn't overlap
        while index >= 1 and self._cds_features[index - 1].is_contained_by(region):
            index -= 1
        # move forwards, adding to the cluster until a cds doesn't overlap
        while index < len(self._cds_features):
            cds = self._cds_features[index]
            if not cds.is_contained_by(region):
                break
            region.add_cds(cds)
            cds.region = region
            index += 1

    def get_aa_translation_from_location(self, location: FeatureLocation,
                                         transl_table: Union[str, int] = None) -> Seq:
        """ Obtain the translation for a feature based on its location """
        if location.end > len(self.seq):
            raise ValueError("location outside available sequence")
        if transl_table is None:
            transl_table = self._transl_table
        extracted = location.extract(self.seq).ungap('-')
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
        seq = Seq(string_version, Alphabet.generic_protein)

        if "-" in str(seq):
            seq = Seq(str(seq).replace("-", ""), Alphabet.generic_protein)
        return seq

    def get_cds_features_within_regions(self) -> List[CDSFeature]:  # pylint: disable=invalid-name
        """ Returns all CDS features in the record that are located within a
            region of interest
        """
        features = []  # type: List[CDSFeature]
        for region in self._regions:
            features.extend(region.cds_children)
        return features

    def create_candidate_clusters(self) -> int:
        """ Takes all Cluster instances and constructs CandidateClusters that cover
            each Cluster. Each combination of overlapping clusters will create
            a CandidateCluster, including a single cluster.

            Returns:
                the number of candidate_clusters created
        """
        if not self._protoclusters:
            return 0

        candidate_clusters = create_candidates_from_protoclusters(self._protoclusters)

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

        areas = []  # type: List[CDSCollection]
        areas.extend(candidate_clusters)
        areas.extend(subregions)
        areas.sort()

        region_location = FeatureLocation(max(0, areas[0].location.start),
                                          min(areas[0].location.end, len(self)))

        candidates = []
        subs = []
        if isinstance(areas[0], CandidateCluster):
            candidates.append(areas[0])
        else:
            assert isinstance(areas[0], SubRegion), type(areas[0])
            subs.append(areas[0])

        regions_added = 0
        for area in areas[1:]:
            if area.overlaps_with(region_location):
                region_location = combine_locations(area.location, region_location)
                if isinstance(area, CandidateCluster):
                    candidates.append(area)
                else:
                    assert isinstance(area, SubRegion), type(area)
                    subs.append(area)
                continue
            # no overlap means new region
            self.add_region(Region(candidates, subs))
            regions_added += 1
            region_location = area.location
            candidates = []
            subs = []
            if isinstance(area, CandidateCluster):
                candidates.append(area)
            else:
                assert isinstance(area, SubRegion), type(area)
                subs.append(area)

        # add the final region being built
        self.add_region(Region(candidates, subs))
        regions_added += 1

        return regions_added

    def get_nrps_pks_cds_features(self) -> List[CDSFeature]:
        """ Returns a list of all CDS features within Clusters that contain at least
            one NRPS/PKS domain.
        """
        return [feature for feature in self.get_cds_features_within_regions() if feature.nrps_pks.domains]

    def get_gc_content(self) -> float:
        """ Calculate the GC content of the record's sequence """
        if not self.seq:
            raise ValueError("Cannot calculate GC content of empty sequence")
        counter = Counter(str(self.seq))
        gc_count = counter['G'] + counter['C'] + counter['g'] + counter['c']
        return gc_count / len(self)

    def strip_antismash_annotations(self) -> None:
        """ Removes all antismash features and annotations from the record """
        self.clear_protoclusters()
        self.clear_candidate_clusters()
        self.clear_subregions()
        self.clear_regions()
        self.clear_antismash_domains()
        self.clear_pfam_domains()

        # clean up antiSMASH-created CDSMotifs, but leave the rest
        motifs = self.get_cds_motifs()
        self.clear_cds_motifs()
        for motif in motifs:
            if not motif.created_by_antismash:
                self.add_cds_motif(motif)

        # clean up antiSMASH annotations in CDS features
        for feature in self.get_cds_features():
            feature.strip_antismash_annotations()
