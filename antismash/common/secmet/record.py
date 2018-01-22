# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Implements a more readable and efficient version of biopython's SeqRecord.

    Most attributes are wrapped, with the record features being secmet features
    and accessible and already paritioned so accessing specific types of feature
    is more efficient.

    Does not, and will not, subclass SeqRecord to avoid mixing code types.
"""


import bisect
import logging
from typing import Dict, List, Tuple, Union

import Bio.Alphabet
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord

from .feature import Feature, CDSFeature, CDSMotif, AntismashDomain, Cluster, \
                     PFAMDomain, ClusterBorder, Prepeptide, Gene


class Record:
    """A record containing secondary metabolite clusters"""
    # slots not for space, but to stop use as a horrible global
    __slots__ = ["_record", "_seq", "skip", "_cds_features", "_cds_by_name",
                 "_clusters", "_cds_by_accession", "original_id",
                 "_cluster_borders", "_cds_motifs", "_pfam_domains", "_antismash_domains",
                 "_cluster_numbering", "_nonspecific_features", "record_index",
                 "_genes"]

    def __init__(self, seq="", **kwargs):
        self._record = SeqRecord(seq, **kwargs)
        self.record_index = None
        self.original_id = None
        self.skip = False  # TODO: move to yet another abstraction layer?
        self._genes = []
        self._cds_features = []
        self._cds_by_accession = {}
        self._cds_by_name = {}
        self._clusters = []
        self._cluster_borders = []
        self._cds_motifs = []
        self._pfam_domains = []
        self._antismash_domains = []
        self._cluster_numbering = {}
        self._nonspecific_features = []

    def __getattr__(self, attr):
        # passthroughs to the original SeqRecord
        if attr in ["id", "seq", "description", "name", "annotations", "dbxrefs"]:
            return getattr(self._record, attr)
        if attr in Record.__slots__:
            return self.__getattribute__(attr)
        raise AttributeError("Record has no attribute '%s'" % attr)

    def __setattr__(self, attr, value):
        # passthroughs to the original SeqRecord
        if attr in ["id", "seq", "description", "name"]:
            return setattr(self._record, attr, value)
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

    def add_annotation(self, key: str, value: List) -> None:
        """Adding annotations in Record"""
        if not isinstance(key, str) or not isinstance(value, (str, list)):
            raise ValueError('Key and Value are not in right format')
        self._record.annotations[key] = value

    def __len__(self):
        return len(self._record)

    def get_clusters(self) -> Tuple:
        """A list of secondary metabolite clusters present in the record"""
        return tuple(self._clusters)

    def get_cluster(self, index: int) -> Cluster:
        """ Get the cluster with the given cluster number """
        return self._clusters[index - 1]  # change from 1-indexed to 0-indexed

    def clear_clusters(self) -> None:
        "Remove all Cluster features and reset CDS linking"
        self._clusters.clear()
        for cds in self.get_cds_features():
            cds.cluster = None

    def get_cluster_borders(self) -> Tuple:
        "Return all ClusterBorder features"
        return tuple(self._cluster_borders)

    def clear_cluster_borders(self) -> None:
        "Remove all ClusterBorder features"
        self._cluster_borders.clear()
        for cluster in self._clusters:
            cluster.borders = []

    def get_genes(self) -> Tuple:
        """ Returns all Gene features (different to CDSFeatures) in the record """
        return tuple(self._genes)

    def get_cds_features(self) -> Tuple:
        """A list of secondary metabolite clusters present in the record"""
        return tuple(self._cds_features)

    def get_cds_accession_mapping(self) -> Dict[str, CDSFeature]:
        """A dictionary mapping CDS accession to CDS feature"""
        return dict(self._cds_by_accession)

    def get_cds_name_mapping(self) -> Dict[str, CDSFeature]:
        """A dictionary mapping CDS name to CDS feature"""
        return dict(self._cds_by_name)

    def get_cds_motifs(self) -> Tuple:
        """A list of secondary metabolite CDS_motifs present in the record"""
        return tuple(self._cds_motifs)

    def clear_cds_motifs(self) -> None:
        "Remove all CDSMotif features"
        self._cds_motifs.clear()

    def get_pfam_domains(self) -> Tuple:
        """A list of secondary metabolite PFAM_domains present in the record"""
        return tuple(self._pfam_domains)

    def get_antismash_domains(self) -> Tuple:
        """A list of secondary metabolite aSDomains present in the record"""
        return tuple(self._antismash_domains)

    def clear_antismash_domains(self) -> None:
        "Remove all AntismashDomain features"
        self._antismash_domains.clear()

    def get_generics(self) -> Tuple:
        """A list of secondary metabolite generics present in the record"""
        return tuple(self._nonspecific_features)

    def get_misc_feature_by_type(self, label: str) -> Tuple:
        """Returns a tuple of all generic features with a type matching label"""
        if label in ["cluster", "cluster_border", "CDS", "CDSmotif",
                     "PFAM_domain", "aSDomain", "aSProdPred"]:
            raise ValueError("Use the appropriate get_* type instead for %s" % label)
        return tuple(i for i in self.get_generics() if i.type == label)

    def get_all_features(self) -> List[Feature]:
        """ Returns all features
            note: This is slow, if only a specific type is required, use
                  the other get_*() functions
        """
        features = list(self.get_generics())
        features.extend(self._genes)
        features.extend(self.get_clusters())
        features.extend(self.get_cluster_borders())
        features.extend(self.get_cds_features())
        features.extend(self.get_cds_motifs())
        features.extend(self.get_antismash_domains())
        features.extend(self.get_pfam_domains())
        return features

    def get_cds_features_within_location(self, location, with_overlapping=False) -> List[Feature]:
        """ Returns all CDS features within the given location

            Arguments:
                location: the location to use as a range
                with_overlapping: whether to include features which overlap the
                                  edges of the range

            Returns:
                a list of CDSFeatures, ordered by earliest position in feature location
        """
        def find_start_in_list(location, features, include_overlaps: bool) -> int:
            """ Find the earliest feature that starts before the location
                (and ends before, if include_overlaps is True)
            """
            dummy = Feature(location, feature_type='dummy')
            index = bisect.bisect_left(features, dummy)
            while index > 0 and features[index -1].location.start == location.start:
                index -= 1
            if include_overlaps:
                while index >= 1 and features[index - 1].overlaps_with(dummy):
                    index -= 1
            return index

        results = []
        index = find_start_in_list(location, self._cds_features, with_overlapping)
        while index < len(self._cds_features):
            feature = self._cds_features[index]
            if feature.is_contained_by(location):
                results.append(feature)
            elif with_overlapping and feature.overlaps_with(location):
                results.append(feature)
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

    def get_cluster_number(self, cluster: Cluster) -> int:
        """Returns cluster number of a cluster feature (1-indexed)
            param cluster : A ClusterFeature instance
        """
        number = self._cluster_numbering.get(cluster)
        if number is None:
            raise ValueError("Cluster not contained in record")
        return number

    def get_feature_count(self) -> int:
        """ Returns the total number of features contained in the record. """
        return sum(map(len, [self._cds_features, self._clusters,
                             self._cluster_borders, self._cds_motifs,
                             self._pfam_domains, self._antismash_domains,
                             self._cluster_numbering, self._nonspecific_features,
                             self._genes]))

    def add_cluster(self, cluster: Cluster) -> None:
        """ Add the given cluster to the record,
            causes cluster-CDS pairing to be recalculated """
        assert isinstance(cluster, Cluster)
        assert cluster.location.start >= 0, cluster
        assert cluster.location.end <= len(self), "%s > %d" % (cluster, len(self))
        index = 0
        for i, existing_cluster in enumerate(self._clusters):  # TODO: fix performance
            if cluster.overlaps_with(existing_cluster):
                logging.error("existing %s overlaps with new %s", existing_cluster, cluster)
                raise ValueError("Clusters cannot overlap")
            if cluster < existing_cluster:
                index = i  # before
                break
            else:
                index = i + 1  # after
        self._clusters.insert(index, cluster)
        cluster.parent_record = self
        # update numbering
        for i in range(index, len(self._clusters)):
            self._cluster_numbering[self._clusters[i]] = i + 1  # 1-indexed
        # link any relevant CDS features
        self._link_cluster_to_cds_features(cluster)
        for cluster_border in self._cluster_borders:
            if cluster_border.is_contained_by(cluster):
                cluster.borders.append(cluster_border)

    def add_cluster_border(self, cluster_border: ClusterBorder) -> None:
        """ Add the given cluster_border to the record and to the cluster it
            overlaps with
        """
        assert isinstance(cluster_border, ClusterBorder)
        assert cluster_border.location.start >= 0
        assert cluster_border.location.end <= len(self)
        self._cluster_borders.append(cluster_border)
        # TODO fix performance
        for cluster in self._clusters:
            if cluster_border.is_contained_by(cluster):
                cluster.borders.append(cluster_border)
                break

    def add_gene(self, gene: Gene) -> None:
        """ Adds a Gene feature to the record """
        assert isinstance(gene, Gene), type(gene)
        self._genes.append(gene)

    def add_cds_feature(self, cds_feature: CDSFeature) -> None:
        """ Add the given CDSFeature to the record,
            causes cluster-CDS pairing to be recalculated """
        assert isinstance(cds_feature, CDSFeature), type(cds_feature)
        # provide a unique id over all records
        cds_feature.unique_id = "%s.%d-%d" % (self.id, cds_feature.location.start,
                                              cds_feature.location.end)
        # ensure it has a translation
        if not cds_feature.translation:
            if not self.seq:
                logging.error('No amino acid sequence in input entry for CDS %s, '
                              'and no nucleotide sequence provided to translate it from.',
                              cds_feature.unique_id)
                raise ValueError("Missing sequence info for CDS %s" % cds_feature.unique_id)
            cds_feature.translation = self.get_aa_translation_of_feature(cds_feature)
        index = bisect.bisect_left(self._cds_features, cds_feature)
        self._cds_features.insert(index, cds_feature)
        self._link_cds_to_parent(cds_feature)
        if cds_feature.get_accession() in self._cds_by_accession:
            raise ValueError("Multiple CDS features have the same accession for mapping: %s" %
                             cds_feature.get_accession())
        if cds_feature.get_name() in self._cds_by_name:
            raise ValueError("Multiple CDS features have the same name for mapping: %s" %
                             cds_feature.get_name())
        self._cds_by_accession[cds_feature.get_accession()] = cds_feature
        self._cds_by_name[cds_feature.get_name()] = cds_feature

    def remove_cds_feature(self, cds_feature: CDSFeature) -> None:
        """ Removes a CDS feature from a record and any associated cluster. """
        assert isinstance(cds_feature, CDSFeature)
        if cds_feature.cluster:
            del cds_feature.cluster.cds_children[cds_feature]
        del self._cds_by_accession[cds_feature.get_accession()]
        del self._cds_by_name[cds_feature.get_name()]
        self._cds_features.remove(cds_feature)

    def add_cds_motif(self, motif: Union[CDSMotif, Prepeptide]) -> None:
        """ Add the given cluster to the record """
        assert isinstance(motif, (CDSMotif, Prepeptide)), "%s, %s" % (type(motif), motif.type)
        self._cds_motifs.append(motif)

    def add_pfam_domain(self, pfam_domain: PFAMDomain) -> None:
        """ Add the given cluster to the record """
        assert isinstance(pfam_domain, PFAMDomain)
        self._pfam_domains.append(pfam_domain)

    def add_antismash_domain(self, antismash_domain: AntismashDomain) -> None:
        """ Add the given cluster to the record """
        assert isinstance(antismash_domain, AntismashDomain)
        self._antismash_domains.append(antismash_domain)

    def add_feature(self, feature: Feature) -> None:
        """ Adds a Feature or any subclass to the relevant list """
        assert isinstance(feature, Feature)
        if isinstance(feature, Cluster):
            self.add_cluster(feature)
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
            # TODO: check insertion order of clusters
            self.add_cds_feature(CDSFeature.from_biopython(feature))
        elif feature.type == 'gene':
            self.add_gene(Gene.from_biopython(feature))
        elif feature.type == 'cluster':
            self.add_cluster(Cluster.from_biopython(feature))  # TODO: fix performance
        elif feature.type == 'CDS_motif':
            self.add_cds_motif(CDSMotif.from_biopython(feature))
        elif feature.type == 'PFAM_domain':
            self.add_pfam_domain(PFAMDomain.from_biopython(feature))
        elif feature.type == 'aSDomain':
            self.add_antismash_domain(AntismashDomain.from_biopython(feature))
        elif feature.type == 'cluster_border':
            self.add_cluster_border(ClusterBorder.from_biopython(feature))
        else:
            self.add_feature(Feature.from_biopython(feature))

    @staticmethod
    def from_biopython(seq_record: SeqRecord) -> "Record":  # string because forward decl
        """ Constructs a new Record instance from a biopython SeqRecord,
            also replaces biopython SeqFeatures with Feature subclasses
        """
        assert isinstance(seq_record, SeqRecord)
        record = Record()
        record._record = seq_record
        for feature in seq_record.features:
            record.add_biopython_feature(feature)
        return record

    def _link_cds_to_parent(self, cds: CDSFeature) -> None:
        """ connect the given CDS to the cluster that contains it, if any """
        assert isinstance(cds, CDSFeature)
        left = bisect.bisect_left(self._clusters, cds)
        right = bisect.bisect_right(self._clusters, cds, lo=left)
        for cluster in self._clusters[left - 1:right + 1]:
            if cds.is_contained_by(cluster):
                cluster.add_cds(cds)
                cds.cluster = cluster

    def _link_cluster_to_cds_features(self, cluster: Cluster) -> None:
        """ connect the given cluster to every CDS feature within it's range """
        assert isinstance(cluster, Cluster)
        # quickly find the first cds with equal start
        index = bisect.bisect_left(self._cds_features, cluster)
        # move backwards until we find one that doesn't overlap
        while index >= 1 and self._cds_features[index - 1].is_contained_by(cluster):
            index -= 1
        # move forwards, adding to the cluster until a cds doesn't overlap
        while index < len(self._cds_features):
            cds = self._cds_features[index]
            if not cds.is_contained_by(cluster):
                break
            cluster.add_cds(cds)
            cds.cluster = cluster  # TODO: allow for multiple parent clusters?
            index += 1

    def get_aa_translation_of_feature(self, feature: Feature) -> Seq:
        """ Obtain content for translation qualifier for specific CDS feature in sequence record"""
        extracted = feature.extract(self.seq).ungap('-')
        if len(extracted) % 3 != 0:
            extracted = extracted[:-(len(extracted) % 3)]
        seq = extracted.translate(to_stop=True)
        if not seq:
            # go past stop codons and hope for something to work with
            seq = extracted.translate()
        if "*" in str(seq):
            seq = Seq(str(seq).replace("*", "X"), Bio.Alphabet.generic_protein)
        if "-" in str(seq):
            seq = Seq(str(seq).replace("-", ""), Bio.Alphabet.generic_protein)
        return seq

    def get_cds_features_within_clusters(self) -> List[CDSFeature]:  # pylint: disable=invalid-name
        """ Returns all CDS features in the record that are located within a
            cluster
        """
        features = []  # type: List[CDSFeature]
        for cluster in self._clusters:
            features.extend(cluster.cds_children)
        return features

    def write_cluster_specific_genbanks(self, output_dir: str = None) -> None:
        """ Write out a set genbank files, each containing a single cluster

            Arguments:
                output_dir: the directory to place the files, if None the
                            current working directory is used
        """
        bio_record = self.to_biopython()
        for cluster in self._clusters:
            cluster.write_to_genbank(directory=output_dir, record=bio_record)

    def create_clusters_from_borders(self) -> int:
        """ Takes all ClusterBorder instances and constructs Clusters that cover
            each ClusterBorder. If a cluster would overlap with another, the
            clusters are merged.

            Returns:
                the number of clusters created
        """
        if not self._cluster_borders:
            return 0
        borders = sorted(self._cluster_borders)
        cluster = Cluster(borders[0].location, borders[0].cutoff,
                          borders[0].extent, borders[0].products)
        if borders[0].rules:
            cluster.detection_rules = borders[0].rules

        clusters_added = 0

        for border in borders[1:]:
            if border.overlaps_with(cluster):
                start = max(cluster.location.start, border.location.start)
                end = max(cluster.location.end, border.location.end)
                cluster.location = FeatureLocation(start, end)
                cluster.products.extend(product for product in border.products if product not in cluster.products)
                cluster.detection_rules.extend(rule for rule in border.rules if rule not in cluster.detection_rules)
            else:
                cluster.contig_edge = cluster.location.start == 0 or cluster.location.end == len(self.seq)
                self.add_cluster(cluster)
                clusters_added += 1
                cluster = Cluster(border.location, border.cutoff, border.extent,
                                  border.products)
                if border.rules:
                    cluster.detection_rules = border.rules

        # add the final cluster being built if it wasn't added already
        cluster.contig_edge = cluster.location.start == 0 or cluster.location.end == len(self.seq)
        self.add_cluster(cluster)
        clusters_added += 1

        return clusters_added
