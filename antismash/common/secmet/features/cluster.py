# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" A class for cluster features """

from collections import OrderedDict
import logging
import os
from typing import Dict, Iterable, List, Optional
from typing import Any  # comment hints # pylint: disable=unused-import
import warnings

from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord
from helperlibs.bio import seqio

from .cds_feature import CDSFeature
from .cluster_border import ClusterBorder  # comment hints # pylint: disable=unused-import
from .feature import Feature, FeatureLocation


class Cluster(Feature):
    """ A feature representing a cluster. Tracks which CDS features belong to it"""
    __slots__ = ["_extent", "_cutoff", "_products", "contig_edge",
                 "detection_rules", "smiles_structure",
                 "clusterblast", "knownclusterblast", "subclusterblast",
                 "parent_record", "cds_children", "borders"]

    def __init__(self, location: FeatureLocation, cutoff: int, extent: int, products: List[str]) -> None:
        super().__init__(location, feature_type="cluster",
                         created_by_antismash=True)

        self._extent = int(extent)
        self._cutoff = int(cutoff)
        self._products = []  # type: List[str]
        for product in products:
            self.add_product(product)

        self.contig_edge = None  # type: Optional[bool] # hmm_detection borderpredict
        self.detection_rules = []  # type: List[str]
        self.smiles_structure = None  # type: Optional[str] # SMILES string

        self.clusterblast = None  # type: Any  # TODO: actually a clusterblast result, remove
        self.knownclusterblast = None  # type: Any  # TODO: actually a clusterblast result, remove
        self.subclusterblast = None  # type: Any  # TODO: actually a clusterblast result, remove

        # for runtime management
        self.parent_record = None  # type: Any  # TODO: optional Record, but that's a circular dependency
        self.cds_children = OrderedDict()  # type: Dict[CDSFeature, None]
        self.borders = []  # type: List[ClusterBorder]

    @property
    def products(self) -> Iterable[str]:
        """ The products of a cluster """
        return tuple(self._products)

    def add_product(self, product: str) -> None:
        """ Add the given product to the cluster's list of products """
        assert product and isinstance(product, str), str(product)
        self._products.append(product)

    def get_cluster_number(self) -> int:
        """ Returns the cluster number which the parent record uses to refer to
            this cluster. """
        if not self.parent_record:
            raise ValueError("Cluster not contained in record")
        return self.parent_record.get_cluster_number(self)

    def trim_overlapping(self) -> None:
        """ Shrinks the cluster, where possible, to exclude any features which
            overlap with the edges of the cluster.
            Any feature fully contained before shrinking will still be fully
            contained.
        """
        if not self.parent_record:
            logging.warning("Trimming cluster which does not belong to a record")
            return
        features = self.parent_record.get_cds_features_within_location(self.location, with_overlapping=True)
        # don't trim if there's no features to trim by
        if not features:
            return

        # find the deepest feature that only overlaps at the beginning
        previous = None
        index = 0
        current = features[index]
        # track where to trim to
        start = self.location.start
        while current.overlaps_with(self) and not current.is_contained_by(self):
            start = max([start, current.location.start, current.location.end])
            previous = current
            index += 1
            if index >= len(features):
                current = None
                break
            current = features[index]

        # don't cause a contained feature to now overlap only
        if previous and current:
            start = min([start, current.location.start, current.location.end])

        # find the deepest feature that only overlaps at the end
        # but skip any indices already covered in the lead search
        lead_index = index
        previous = None
        index = len(features) - 1
        current = features[index]
        # track where to trim to
        end = self.location.end
        while index > lead_index and current.overlaps_with(self) and not current.is_contained_by(self):
            end = min([end, current.location.start, current.location.end])
            previous = current
            index -= 1
            if index < 0:
                current = None
                break
            current = features[index]

        # but don't cause a contained feature to now overlap only
        if previous and current:
            end = max([end, current.location.start, current.location.end])

        # finally, do the trim itself
        new_loc = FeatureLocation(start, end, self.location.strand)
        if self.location.start != start or self.location.end != end:
            logging.debug("Cluster %d trimming location from %s to %s",
                          self.get_cluster_number(), self.location, new_loc)
        # make sure the size is never increased
        assert self.location.start <= start < end <= self.location.end
        self.location = new_loc

        for cds in self.cds_children:
            assert cds.is_contained_by(self), "cluster trimming removed wholly contained CDS"

    def add_cds(self, cds: CDSFeature) -> None:
        """ Adds a CDSFeature to the cluster """
        assert isinstance(cds, CDSFeature)
        assert cds.is_contained_by(self), "cds %s outside cluster %s" % (cds, self)
        self.cds_children[cds] = None

    @property
    def cutoff(self) -> int:
        """ The maximal distance between genes when defining the cluster.
            The distance between core genes after definition will likely be
            smaller."""
        return self._cutoff

    @cutoff.setter
    def cutoff(self, cutoff: int) -> None:
        if cutoff is not None:
            cutoff = int(cutoff)
        self._cutoff = cutoff

    @property
    def extent(self) -> int:
        """ The distance the cluster extends from the first and last genes which
            from which the cluster was defined.
        """
        return self._extent

    @extent.setter
    def extent(self, extent: int) -> None:
        if extent is not None:
            self._extent = int(extent)

    def get_product_string(self) -> str:
        """ Returns the cluster's products as a single string """
        assert None not in self._products, self._products
        return "-".join(self._products)

    @property
    def probability(self) -> Optional[float]:
        """ The cluster probability, if relevant. """
        probabilities = {border.probability for border in self.borders}
        # if all agree on the probability
        if len(probabilities) == 1:
            return list(probabilities)[0]
        # if they disagree, return None again
        return None

    @staticmethod
    def from_biopython(bio_feature: SeqFeature, feature: "Cluster" = None,  # type: ignore
                       leftovers: Dict[str, List[str]] = None) -> "Cluster":
        if leftovers is None:
            leftovers = Feature.make_qualifiers_copy(bio_feature)
        # grab mandatory qualifiers and create the class
        cutoff = int(leftovers.pop("cutoff")[0])
        extent = int(leftovers.pop("extension")[0])
        products = leftovers.pop("product")[0].split("-")
        if not feature:
            feature = Cluster(bio_feature.location, cutoff, extent, products)
        # take the detection rules from "note"
        # first check it exists
        index = -1
        for i, note in enumerate(leftovers.get("note", [])):
            if note.startswith("Detection rule"):
                index = i
                break
        # then split it into the relevant pieces
        if index > -1:
            products = []
            rules = []
            text = leftovers["note"].pop(index)
            text = text.split(":", 1)[1]  # strip the leadin
            products_rules = text.split(";")  # separate rules
            for product_rule in products_rules:
                if not product_rule:
                    continue
                assert ": " in product_rule, product_rule
                product, rule = product_rule.split(": ")
                rule = rule[1:-1]  # strip ( )
                products.append(product.strip())
                rules.append(rule)
            feature.detection_rules = rules

        # grab optional qualifiers
        contig_edge = leftovers.pop("contig_edge", [None])[0]
        if not contig_edge:
            feature.contig_edge = None
        else:
            feature.contig_edge = contig_edge == "True"
        if "structure" in leftovers:
            feature.smiles_structure = leftovers.pop("structure")[0]
        # grab optional parent qualifiers
        updated = super(Cluster, feature).from_biopython(bio_feature, feature, leftovers)
        assert updated is feature, "feature changed: %s -> %s" % (feature, updated)
        assert isinstance(updated, Cluster)
        return updated

    def to_biopython(self, qualifiers: Optional[Dict[str, List[str]]] = None) -> List[SeqFeature]:
        mine = OrderedDict()  # type: Dict[str, List[str]]
        mine["cutoff"] = [str(self.cutoff)]
        mine["extension"] = [str(self.extent)]
        if self.contig_edge is not None:
            mine["contig_edge"] = [str(self.contig_edge)]
        if self.smiles_structure is not None:
            mine["structure"] = [self.smiles_structure]
        assert isinstance(self._products, list), type(self.products)
        mine["product"] = [self.get_product_string()]
        rule_text = ["Detection rule(s) for this cluster type:"]
        assert isinstance(self.detection_rules, list), type(self.detection_rules)
        for product, rule in zip(self.products, self.detection_rules):
            rule_text.append("%s: (%s);" % (product, rule))
        if qualifiers:
            mine.update(qualifiers)
        if "note" not in mine:
            mine["note"] = []
        mine["note"].append(" ".join(rule_text))
        return super().to_biopython(mine)

    def write_to_genbank(self, filename: str = None, directory: str = None, record: SeqRecord = None) -> None:
        """ Writes a genbank file containing only the information contained
            within the Cluster.
        """
        if not filename:
            filename = "%s.cluster%03d.gbk" % (self.parent_record.id, self.get_cluster_number())
        if directory:
            filename = os.path.join(directory, filename)

        if record is None:
            record = self.parent_record.to_biopython()
        assert isinstance(record, SeqRecord)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            cluster_record = record[self.location.start:self.location.end]

        cluster_record.annotations["date"] = record.annotations.get("date", '')
        cluster_record.annotations["source"] = record.annotations.get("source", '')
        cluster_record.annotations["organism"] = record.annotations.get("organism", '')
        cluster_record.annotations["taxonomy"] = record.annotations.get("taxonomy", [])
        cluster_record.annotations["data_file_division"] = record.annotations.get("data_file_division", 'UNK')
        # our cut-out clusters are always linear
        cluster_record.annotations["topology"] = "linear"

        seqio.write([cluster_record], filename, 'genbank')
