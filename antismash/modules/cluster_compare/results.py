# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Contains the results classes for the ClusterCompare module """

import dataclasses
from typing import Any, Dict, List, Optional

from antismash.common.module_results import ModuleResults
from antismash.common.secmet import Record

from .data_structures import (
    Hit,
    ReferenceProtocluster,
    ReferenceRegion,
    ReferenceScorer,
    ScoresByRegion,
)


@dataclasses.dataclass
class ScoresByProtocluster:
    """ A base class for containing details about scoring of query and reference
        area pairings
    """
    kind = "unknown"
    details: Any

    def to_json(self) -> Dict[str, Any]:
        """ Constructs a JSON representation of the instance """
        raise NotImplementedError

    @classmethod
    def from_json(cls, data: Dict[str, Any], record: Record,
                  references: Dict[str, ReferenceRegion]) -> "ScoresByProtocluster":
        """ Reconstructs a ScoresByProtocluster instance from a JSON representation

            Requires the matching Record instance to rebuild links to CDS features,
            along with a mapping of ReferenceRegions to rebuilds links to references.
        """
        raise NotImplementedError


@dataclasses.dataclass
class ProtoToProtoScores(ScoresByProtocluster):
    """ Contains details of the scorings of reference protoclusters to query
        protoclusters
    """
    kind = "ProtoToProto"
    details: Dict[int, Dict[ReferenceRegion, Dict[ReferenceProtocluster, ReferenceScorer]]]

    def to_json(self) -> Dict[str, Any]:
        details: Dict[int, Dict[str, Dict[str, Dict[str, Any]]]] = {}
        for proto, ref_results in self.details.items():
            details[proto] = {}
            for ref_region, proto_results in ref_results.items():
                scores: Dict[str, Dict[str, Any]] = {}
                details[proto][ref_region.get_identifier()] = scores
                for ref_proto, scorer in proto_results.items():
                    scores[ref_proto.get_identifier()] = scorer.to_json()
        data = {
            "kind": self.kind,
            "details": details,
        }
        return data

    @classmethod
    def from_json(cls, data: Dict[str, Any], record: Record,
                  references: Dict[str, ReferenceRegion]) -> "ProtoToProtoScores":
        assert data["kind"] == cls.kind
        details: Dict[int, Dict[ReferenceRegion, Dict[ReferenceProtocluster, ReferenceScorer]]] = {}
        for proto, proto_results in data["details"].items():
            proto_num = int(proto)
            details[proto_num] = {}
            for ref_region, ref_proto_results in proto_results.items():
                converted: Dict[ReferenceProtocluster, ReferenceScorer] = {}
                details[proto_num][references[ref_region]] = converted
                for ref_proto_id, result in ref_proto_results.items():
                    scorer = ReferenceScorer.from_json(result, record, references)
                    for protocluster in references[ref_region].protoclusters:
                        if protocluster.get_identifier() == ref_proto_id:
                            converted[protocluster] = scorer
                            break
                    else:
                        raise ValueError(f"reference region {ref_proto_id} does not contain {ref_proto_id}")
        return cls(details)


@dataclasses.dataclass
class RegionToRegionScores(ScoresByProtocluster):
    """ Contains details of the scorings of reference regions to query regions
    """
    kind = "RegionToRegion"
    details: List[ReferenceScorer]

    def to_json(self) -> Dict[str, Any]:
        data = {
            "kind": self.kind,
            "details": [scorer.to_json() for scorer in self.details],
        }
        return data

    @classmethod
    def from_json(cls, data: Dict[str, Any], record: Record,
                  references: Dict[str, ReferenceRegion]) -> "RegionToRegionScores":
        return cls([ReferenceScorer.from_json(chunk, record, references) for chunk in data["details"]])


class ProtoToRegionScores(ScoresByProtocluster):
    """ Contains details of the scorings of reference protoclusters to query
        regions
    """
    kind = "ProtoToRegion"
    details: Dict[int, Dict[ReferenceRegion, ReferenceScorer]]

    def to_json(self) -> Dict[str, Any]:
        scores_by_proto = {}
        for proto_num, proto_results in self.details.items():
            for_proto: Dict[str, Dict[str, Dict[str, Any]]] = {}
            scores_by_proto[proto_num] = for_proto
            for ref_region, ref_proto_results in proto_results.items():
                for_proto[ref_region.get_identifier()] = ref_proto_results.to_json()

        data = {
            "kind": self.kind,
            "details": scores_by_proto,
        }
        return data

    @classmethod
    def from_json(cls, data: Dict[str, Any], record: Record,
                  references: Dict[str, ReferenceRegion]) -> "ProtoToRegionScores":
        details: Dict[int, Dict[ReferenceRegion, ReferenceScorer]] = {}
        for proto, proto_results in data["details"].items():
            proto_num = int(proto)
            details[proto_num] = {}
            for ref_region, ref_result in proto_results.items():
                scorer = ReferenceScorer.from_json(ref_result, record, references)
                details[proto_num][references[ref_region]] = scorer
        return cls(details)


@dataclasses.dataclass
class VariantResults:
    """ Contains the results for a particular search variant of a database and
        the query record.
    """
    variant_name: str
    scores_by_region: ScoresByRegion
    details: ScoresByProtocluster
    hits_by_region: Dict[ReferenceRegion, Dict[str, Hit]]

    def __post_init__(self) -> None:
        if not self.variant_name:
            raise ValueError("variant name cannot be empty")

    def to_json(self) -> Dict[str, Any]:
        """ Constructs a JSON representation of the VariantResults """
        scores = {}
        for ref, score in self.scores_by_region:
            scores[ref.get_identifier()] = score

        json_regions = {}
        for ref, _ in self.scores_by_region:
            as_json = ref.to_json()
            as_json["cds_mapping"] = ref.cds_mapping
            as_json["accession"] = ref.accession
            json_regions[ref.get_identifier()] = as_json

        json_hits = {}
        for ref, hits in self.hits_by_region.items():
            json_hits[ref.get_identifier()] = {name: hit.to_json() for name, hit in hits.items()}

        return {
            "name": self.variant_name,
            "scores_by_region": scores,
            "reference_regions": json_regions,
            "details": self.details.to_json(),
            "hits": json_hits,
        }

    @classmethod
    def from_json(cls, data: Dict[str, Any], record: Record) -> "VariantResults":
        """ Reconstructs a VariantResults instance from a JSON representation.

            Requires the matching Record instance to rebuild links to CDS features
        """
        ref_regions = {}
        for name, raw in data["reference_regions"].items():
            accession = raw["accession"]
            cds_mapping = raw["cds_mapping"]
            ref_regions[name] = ReferenceRegion.from_json(accession, raw, cds_mapping)
        scores_by_region: ScoresByRegion = []
        for name, score in data["scores_by_region"].items():
            scores_by_region.append((ref_regions[name], score))

        hits: Dict[ReferenceRegion, Dict[str, Hit]] = {}
        for name, raw_hits in data["hits"].items():
            hits[ref_regions[name]] = {hit_name: Hit.from_json(raw_hit, record)
                                       for hit_name, raw_hit in raw_hits.items()}

        details: ScoresByProtocluster
        kind = data["details"]["kind"]
        assert kind
        for option in [RegionToRegionScores, ProtoToProtoScores, ProtoToRegionScores]:
            assert issubclass(option, ScoresByProtocluster)
            if option.kind == kind:
                details = option.from_json(data["details"], record, ref_regions)
                break
        else:
            raise ValueError(f"unknown scoring type: '{kind}'")
        return cls(data["name"], scores_by_region, details, hits)


@dataclasses.dataclass
class DatabaseResults:
    """ Contains all results for a single reference database and the query record """
    name: str
    url: Optional[str]
    description: Optional[str]
    by_region: Dict[int, Dict[str, VariantResults]]

    def to_json(self) -> Dict[str, Any]:
        """ Constructs a JSON representation of the DatabaseResults """
        results: Dict[str, Any] = {
            "name": self.name,
            "by_region": {},
        }

        for region, variant_results in self.by_region.items():
            converted = {variant: result.to_json() for variant, result in variant_results.items()}
            results["by_region"][str(region)] = converted

        if self.url:
            results["url_format"] = self.url
        if self.description:
            results["description"] = self.description

        return results

    @classmethod
    def from_json(cls, data: Dict[str, Any], record: Record) -> "DatabaseResults":
        """ Reconstructs a DatabaseResults instance from a JSON representation.

            Requires the matching Record instance to rebuild links to CDS features
        """
        name = data["name"]
        url = data.get("url_format")
        description = data.get("description")
        results = {}
        for region, variant_results in data["by_region"].items():
            converted = {}
            for variant, result in variant_results.items():
                converted[variant] = VariantResults.from_json(result, record)
            results[int(region)] = converted
        return cls(name, url, description, results)


class ClusterCompareResults(ModuleResults):
    """ The results of cluster comparison """
    _schema_version = 1

    def __init__(self, record_id: str, by_database: Optional[Dict[str, DatabaseResults]]) -> None:
        super().__init__(record_id)
        # region number -> database name/description -> variant -> result
        self.by_database = by_database or {}

    def to_json(self) -> Dict[str, Any]:
        converted: Dict[str, Any] = {name: result.to_json() for name, result in self.by_database.items()}
        return {
            "record_id": self.record_id,
            "schema_version": self._schema_version,
            "db_results": converted,
        }

    @staticmethod
    def from_json(json: Dict[str, Any], record: Record) -> Optional["ClusterCompareResults"]:
        if ClusterCompareResults._schema_version != json["schema_version"]:
            raise ValueError("schema version mismatch in clustercompare results")
        if record.id != json["record_id"]:
            raise ValueError("record ID mismatch in clustercompare results")
        converted = {name: DatabaseResults.from_json(db, record) for name, db in json["db_results"].items()}
        return ClusterCompareResults(record.id, converted)

    def add_to_record(self, record: Record) -> None:
        pass
