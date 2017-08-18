# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

import bisect
import logging

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .feature import Feature, CDSFeature, CDSMotif, AntismashDomain, Cluster, \
                     PFAMDomain, ClusterBorder

class _BisectHelper:
    def __init__(self, features):
        self.features = features

    def __len__(self):
        return len(self.features)

    def __getitem__(self, index):
        loc = self.features[index].location
        return (loc.start, loc.end)

class Record:
    """A record containing secondary metabolite clusters"""

    def __init__(self, seq=None, **kwargs):
        self._record = SeqRecord(seq, **kwargs)
        self.skip = False #TODO: move to yet another abstraction layer?
        self._cds_features = []
        self._clusters = []
        self._cluster_borders = []
        self._cds_motifs = []
        self._pfam_domains = []
        self._antismash_domains = []
        self._cluster_numbering = {}
        self._nonspecific_features = []

    def __getattr__(self, attr):
        # passthroughs to the original SeqRecord
        if attr in ["id", "seq", "description", "name", "annotations"]:
            return getattr(self._record, attr)
        if attr not in self.__dict__:
            raise AttributeError("Record has no attribute: %s" % attr)
        # something of ours
        return self.__dict__[attr]

    def __setattr__(self, attr, value):
        # passthroughs to the original SeqRecord
        if attr in ["id", "seq", "description", "name"]:
            setattr(self._record, attr, value)
        if attr in ["annotations"]:
            assert isinstance(value, dict)
            for key, val in value.items():
                self.add_annotation(key, val)
        else:
            super().__setattr__(attr, value)

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

    def get_cds_features(self):
        """A list of secondary metabolite clusters present in the record"""
        return tuple(self._cds_features)

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
                         dbxrefs=self._record.dbxrefs, features=bio_features,
                         annotations=self._record.annotations,
                         letter_annotations=self._record.letter_annotations)

    def get_cluster_number(self, cluster):
        """Returns cluster number of a cluster feature
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
            self._cluster_numbering[self._clusters[i]] = i
        # link any relevant CDS features
        self._link_cluster_to_cds_features(cluster)

    def add_cluster_border(self, cluster_border):
        assert isinstance(cluster_border, ClusterBorder)
        self._cluster_borders.append(cluster_border)

    def add_cds_feature(self, cds_feature):
        """ Add the given cluster to the record,
            causes cluster-CDS pairing to be recalculated """
        assert isinstance(cds_feature, CDSFeature)
        index = bisect.bisect_left(self._cds_features, cds_feature)
        self._cds_features.insert(index, cds_feature)
        self._link_cds_to_parent(cds_feature)

    def add_cds_motif(self, motif):
        """ Add the given cluster to the record """
        assert isinstance(motif, CDSMotif), "%s, %s" %(type(motif), motif.type)
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
        elif isinstance(feature, CDSMotif):
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
            if cds.overlaps_with(cluster):
                if cds not in cluster.cds_children: # TODO: fix performance
                    cluster.cds_children.append(cds)
                cds.cluster = cluster

    def _link_cluster_to_cds_features(self, cluster):
        assert isinstance(cluster, Cluster)
        # quickly find the first cds with equal start
        index = bisect.bisect_left(self._cds_features, cluster)
        # move backwards until we find one that doesn't overlap
        cluster_start = cluster.location.start
        while 1 <= index and self._cds_features[index - 1].location.end > cluster_start:
            index -= 1
        # move forwards, adding to the cluster until a cds doesn't overlap
        while index < len(self._cds_features):
            cds = self._cds_features[index]
            if not cds.overlaps_with(cluster):
                break
            cluster.cds_children.append(cds)
            cds.cluster = cluster # TODO: allow for multiple parent clusters
            index += 1
