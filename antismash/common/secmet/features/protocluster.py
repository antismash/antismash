# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" A class for protocluster features """

from typing import Dict, List, Optional, Set

from Bio.SeqFeature import SeqFeature

from .cds_feature import CDSFeature
from .cdscollection import CDSCollection
from .feature import Feature, FeatureLocation
from ..locations import location_from_string
from ..qualifiers.t2pks import T2PKSQualifier
from ..qualifiers.gene_functions import GeneFunction


class Protocluster(CDSCollection):
    """ A feature which marks a specific region of a record as interesting.
        Protoclusters are only those determined by a rule-based method of detection
        and with a defined product.

        The core location covers the first and last CDS feature that caused the
        protocluster to be formed, while the surrounding location includes the context.
    """
    core_seqfeature_type = "proto_core"
    __slots__ = ["core_location", "detection_rule", "product", "tool", "cutoff",
                 "_definition_cdses", "neighbourhood_range", "t2pks"]

    def __init__(self, core_location: FeatureLocation, surrounding_location: FeatureLocation,
                 tool: str, product: str, cutoff: int, neighbourhood_range: int,
                 detection_rule: str) -> None:
        super().__init__(surrounding_location, feature_type="protocluster")
        # cluster-wide
        self.detection_rule = detection_rule
        if not product.replace("-", "").replace("_", "").isalnum() or product[0] in "-_" or product[-1] in "-_":
            raise ValueError("invalid protocluster product: %s" % product)
        self.product = product
        self.tool = tool

        # core specific
        self.core_location = core_location  # type: FeatureLocation
        self.cutoff = cutoff
        self._definition_cdses = set()  # type: Set[CDSFeature]

        # neighbourhood specific
        self.neighbourhood_range = neighbourhood_range

        # analysis annotations
        self.t2pks = None  # type: Optional[T2PKSQualifier]

    def __str__(self) -> str:
        return "Protocluster(%s, product=%s)" % (self.location, self.product)

    def get_protocluster_number(self) -> int:
        """ Returns the protoclusters's numeric ID, only guaranteed to be consistent
            when the same protoclusters are defined in the parent record
        """
        if not self._parent_record:
            raise ValueError("Protocluster not in a record")
        return self._parent_record.get_protocluster_number(self)

    @property
    def definition_cdses(self) -> Set[CDSFeature]:
        """ Returns the set of CDSFeatures responsible for the creation of this protocluster """
        return set(self._definition_cdses)

    def add_cds(self, cds: CDSFeature) -> None:
        super().add_cds(cds)
        if not cds.is_contained_by(self.core_location):
            return
        cores = cds.gene_functions.get_by_function(GeneFunction.CORE)
        if any(core.product == self.product for core in cores):
            self._definition_cdses.add(cds)

    def to_biopython(self, qualifiers: Optional[Dict[str, List[str]]] = None) -> List[SeqFeature]:
        common = {
            "neighbourhood": [str(self.neighbourhood_range)],
            "cutoff": [str(self.cutoff)],
            "product": [self.product],
            "aStool": [self.tool],
            "detection_rule": [self.detection_rule]
        }
        if self._parent_record:
            common["protocluster_number"] = [str(self.get_protocluster_number())]

        shared_qualifiers = dict(qualifiers) if qualifiers else {}
        shared_qualifiers.update(common)
        if self.t2pks:
            shared_qualifiers.update(self.t2pks.to_biopython_qualifiers())

        core_qualifiers = dict(self._qualifiers)
        core_qualifiers.update(shared_qualifiers)
        core_feature = SeqFeature(self.core_location, type=self.core_seqfeature_type)
        for key, val in sorted(core_qualifiers.items()):
            core_feature.qualifiers[key] = val
            # this doesn't go through to the Feature annotations, so add the tool explicitly
            core_feature.qualifiers["tool"] = ["antismash"]

        shared_qualifiers["core_location"] = [str(self.core_location)]
        neighbourhood_feature = super().to_biopython(shared_qualifiers)[0]
        return [neighbourhood_feature, core_feature]

    @staticmethod
    def from_biopython(bio_feature: SeqFeature, feature: "Protocluster" = None,  # type: ignore
                       leftovers: Dict[str, List[str]] = None) -> "Protocluster":
        assert bio_feature.type == "protocluster"
        if leftovers is None:
            leftovers = Feature.make_qualifiers_copy(bio_feature)
        # grab mandatory qualifiers and create the class
        neighbourhood_range = int(leftovers.pop("neighbourhood")[0])
        cutoff = int(leftovers.pop("cutoff")[0])
        product = leftovers.pop("product")[0]
        tool = leftovers.pop("aStool")[0]
        rule = leftovers.pop("detection_rule")[0]
        core_location = location_from_string(leftovers.pop("core_location")[0])
        if not feature:
            feature = Protocluster(core_location, bio_feature.location,
                              tool, product, cutoff, neighbourhood_range, rule)

        # remove run-specific info
        leftovers.pop("protocluster_number", "")

        # rebuild analysis annotations
        feature.t2pks = T2PKSQualifier.from_biopython_qualifiers(leftovers)

        # grab optional parent qualifiers
        updated = super(Protocluster, feature).from_biopython(bio_feature, feature, leftovers)
        assert updated is feature, "feature changed: %s -> %s" % (feature, updated)
        assert isinstance(updated, Protocluster)
        return updated
