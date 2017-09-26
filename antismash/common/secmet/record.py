# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

import bisect
import logging

import Bio.Alphabet
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .feature import Feature, CDSFeature, CDSMotif, AntismashDomain, Cluster, \
                     PFAMDomain, ClusterBorder, Prepeptide

class Record:
    """A record containing secondary metabolite clusters"""
    # slots not for space, but to stop use as a horrible global
    __slots__ = ["_record", "_seq", "skip", "_cds_features", "_cds_mapping", "_clusters",
                 "_cluster_borders", "_cds_motifs", "_pfam_domains", "_antismash_domains",
                 "_cluster_numbering", "_nonspecific_features", "record_index"]
    def __init__(self, seq=None, **kwargs):
        self._record = SeqRecord(seq, **kwargs)
        self.record_index = None
        self.skip = False #TODO: move to yet another abstraction layer?
        self._cds_features = []
        self._cds_mapping = {} # maps CDS accession to CDS feature
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

    def add_annotation(self, key, value):
        """Adding annotations in Record"""
        if not isinstance(key, str) or not isinstance(value, (str, list)):
            raise ValueError('Key and Value are not in right format')
        self._record.annotations[key] = value

    def __len__(self):
        return len(self._record)

    def get_clusters(self):
        """A list of secondary metabolite clusters present in the record"""
        return tuple(self._clusters)

    def get_cluster(self, index):
        """ Get the cluster with the given cluster number """
        return self._clusters[index - 1] # change from 1-indexed to 0-indexed

    def clear_clusters(self):
        "Remove all Cluster features and reset CDS linking"
        self._clusters.clear()
        for cds in self.get_cds_features():
            cds.cluster = None

    def get_cluster_borders(self):
        "Return all ClusterBorder features"
        return tuple(self._cluster_borders)

    def clear_cluster_borders(self):
        "Remove all ClusterBorder features"
        self._cluster_borders.clear()
        for cluster in self._clusters:
            cluster.borders = []

    def get_cds_features(self):
        """A list of secondary metabolite clusters present in the record"""
        return tuple(self._cds_features)

    def get_cds_mapping(self):
        """A dictionary of CDS accession to CDS feature"""
        return dict(self._cds_mapping)

    def get_cds_motifs(self):
        """A list of secondary metabolite CDS_motifs present in the record"""
        return tuple(self._cds_motifs)

    def clear_cds_motifs(self):
        "Remove all CDSMotif features"
        self._cds_motifs.clear()

    def get_pfam_domains(self):
        """A list of secondary metabolite PFAM_domains present in the record"""
        return tuple(self._pfam_domains)

    def get_antismash_domains(self):
        """A list of secondary metabolite aSDomains present in the record"""
        return tuple(self._antismash_domains)

    def clear_antismash_domains(self):
        "Remove all AntismashDomain features"
        self._antismash_domains.clear()

    def get_generics(self):
        """A list of secondary metabolite generics present in the record"""
        return tuple(self._nonspecific_features)

    def get_misc_feature_by_type(self, label):
        """Returns a tuple of all generic features with a type matching label"""
        if label in ["cluster", "cluster_border", "CDS", "CDSmotif",
                     "PFAM_domain", "aSDomain", "aSProdPred"]:
            raise ValueError("Use the appropriate get_* type instead for %s" % label)
        return tuple(i for i in self.get_generics() if i.type == label)

    def get_all_features(self):
        """ Returns all features
            note: This is slow, if only a specific type is required, use
                  the other get_*() functions
        """
        features = list(self.get_generics())
        features.extend(self.get_clusters())
        features.extend(self.get_cluster_borders())
        features.extend(self.get_cds_features())
        features.extend(self.get_cds_motifs())
        features.extend(self.get_antismash_domains())
        features.extend(self.get_pfam_domains())
        return features

    def to_biopython(self):
        """Returns a Bio.SeqRecord instance of the record"""
        features = self.get_all_features()
        bio_features = []
        for feature in sorted(features):
            bio_features.extend(feature.to_biopython())
        return SeqRecord(self.seq, id=self._record.id, name=self._record.name,
                         description=self._record.description,
                         dbxrefs=self.dbxrefs, features=bio_features,
                         annotations=self._record.annotations,
                         letter_annotations=self._record.letter_annotations)

    def get_cluster_number(self, cluster):
        """Returns cluster number of a cluster feature (1-indexed)
            param cluster : A ClusterFeature instance
        """
        number = self._cluster_numbering.get(cluster)
        if number is None:
            raise ValueError("Cluster not contained in record")
        return number

    def get_feature_count(self):
        """ Returns the total number of features contained in the record. """
        return sum(map(len, [self._cds_features, self._clusters,
                        self._cluster_borders, self._cds_motifs,
                        self._pfam_domains, self._antismash_domains,
                        self._cluster_numbering, self._nonspecific_features]))

    def add_cluster(self, cluster):
        """ Add the given cluster to the record,
            causes cluster-CDS pairing to be recalculated """
        assert isinstance(cluster, Cluster)
        index = 0
        for i, existing_cluster in enumerate(self._clusters): # TODO: fix performance
            if cluster.overlaps_with(existing_cluster):
                raise ValueError("Clusters cannot overlap")
            if cluster < existing_cluster:
                index = i # before
                break
            else:
                index = i + 1 # after
        self._clusters.insert(index, cluster)
        cluster.parent_record = self
        # update numbering
        for i in range(index, len(self._clusters)):
            self._cluster_numbering[self._clusters[i]] = i + 1 # 1-indexed
        # link any relevant CDS features
        self._link_cluster_to_cds_features(cluster)
        for cluster_border in self._cluster_borders:
            if cluster.overlaps_with(cluster_border):
                if cluster_border.parent is not None:
                    raise ValueError("A cluster border is overlapping with two clusters")
                cluster.borders.append(cluster_border)
                cluster_border.parent = cluster
                break

    def add_cluster_border(self, cluster_border):
        assert isinstance(cluster_border, ClusterBorder)
        self._cluster_borders.append(cluster_border)
        # TODO fix performance
        for cluster in self._clusters:
            if cluster.overlaps_with(cluster_border):
                cluster.borders.append(cluster_border)
                cluster_border.parent = cluster
                break

    def add_cds_feature(self, cds_feature):
        """ Add the given cluster to the record,
            causes cluster-CDS pairing to be recalculated """
        assert isinstance(cds_feature, CDSFeature)

        # ensure it has a translation
        if not cds_feature.translation:
            cds_feature.translation = self.get_aa_translation_of_feature(cds_feature)
        index = bisect.bisect_left(self._cds_features, cds_feature)
        self._cds_features.insert(index, cds_feature)
        self._link_cds_to_parent(cds_feature)
        if cds_feature.get_accession() in self._cds_mapping:
            logging.critical("Multiple CDS features have the same accession for mapping")
        self._cds_mapping[cds_feature.get_accession()] = cds_feature
        cds_feature.unique_id = self.id + str(cds_feature.location)

    def add_cds_motif(self, motif):
        """ Add the given cluster to the record """
        assert isinstance(motif, (CDSMotif, Prepeptide)), "%s, %s" %(type(motif), motif.type)
        self._cds_motifs.append(motif)

    def add_pfam_domain(self, pfam_domain):
        """ Add the given cluster to the record """
        assert isinstance(pfam_domain, PFAMDomain)
        self._pfam_domains.append(pfam_domain)

    def add_antismash_domain(self, antismash_domain):
        """ Add the given cluster to the record """
        assert isinstance(antismash_domain, AntismashDomain)
        self._antismash_domains.append(antismash_domain)

    def add_feature(self, feature):
        """ Adds a Feature or any subclass to the relevant list """
        assert isinstance(feature, Feature)
        if isinstance(feature, Cluster):
            self.add_cluster(feature)
        elif isinstance(feature, CDSFeature):
            self.add_cds_feature(feature)
        elif isinstance(feature, (CDSMotif, Prepeptide)):
            self.add_cds_motif(feature)
        elif isinstance(feature, PFAMDomain):
            self.add_pfam_domain(feature)
        elif isinstance(feature, AntismashDomain):
            self.add_antismash_domain(feature)
        else:
            self._nonspecific_features.append(feature)

    def add_biopython_feature(self, feature):
        if feature.type == 'CDS':
            # TODO: check insertion order of clusters
            self.add_cds_feature(CDSFeature.from_biopython(feature))
        elif feature.type == 'cluster':
            self.add_cluster(Cluster.from_biopython(feature)) # TODO: fix performance
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
    def from_biopython(seq_record):
        """ Constructs a new Record instance from a biopython SeqRecord,
            also replaces biopython SeqFeatures with Feature subclasses
        """
        assert isinstance(seq_record, SeqRecord)
        record = Record()
        record._record = seq_record
        for feature in seq_record.features:
            record.add_biopython_feature(feature)
        return record

    def _link_cds_to_parent(self, cds):
        assert isinstance(cds, CDSFeature)
        left = bisect.bisect_left(self._clusters, cds)
        right = bisect.bisect_right(self._clusters, cds, lo=left)
        for cluster in self._clusters[left - 1:right + 1]:
            if cds.is_contained_by(cluster):
                cluster.add_cds(cds)
                cds.cluster = cluster

    def _link_cluster_to_cds_features(self, cluster):
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
            cds.cluster = cluster # TODO: allow for multiple parent clusters
            index += 1

    def get_aa_translation_of_feature(self, feature):
        """Obtain content for translation qualifier for specific CDS feature in sequence record"""
        seq = feature.extract(self.seq).ungap('-').translate(to_stop=True)
        if not seq:
            seq = feature.extract(self.seq).ungap('-').translate()
        if "*" in str(seq):
            seq = Seq(str(seq).replace("*", "X"), Bio.Alphabet.generic_protein)
        if "-" in str(seq):
            seq = Seq(str(seq).replace("-", ""), Bio.Alphabet.generic_protein)
        return seq

    def write_cluster_specific_genbanks(self, output_dir=None):
        bio_record = self.to_biopython()
        for cluster in self._clusters:
            cluster.write_to_genbank(directory=output_dir, record=bio_record)
