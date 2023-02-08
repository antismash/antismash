# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" A collection of data structures for use in the rest of the module"""

import dataclasses
from enum import Enum, unique
import json
from typing import Any, Dict, List, Optional, Tuple

from antismash.common.secmet import CDSFeature, Record
from antismash.common.secmet.locations import location_from_string, FeatureLocation, Location


@unique
class Mode(Enum):
    """ The mode of an analysis, as searching for a reference area in a query area
        and a query in a reference area can have wildly different results.
    """
    QUERY_IN_REFERENCE = 0
    REFERENCE_IN_QUERY = 1
    BEST = 2

    def __str__(self) -> str:
        if self == Mode.REFERENCE_IN_QUERY:
            return "RiQ"
        if self == Mode.QUERY_IN_REFERENCE:
            return "QiR"
        return str(self.name).lower()

    @staticmethod
    def from_string(string: str) -> "Mode":
        """ Returns a Mode instance that has the matching string representation """
        for mode in Mode:
            if string == str(mode):
                return mode
        raise ValueError(f"unknown Mode string: '{string}'")


@dataclasses.dataclass(frozen=True)
class DBConfig:
    """ A simple container for passing around the configuration of a specific
        reference database.
    """
    name: str
    comparisons: List[str]
    path: str
    mode: Mode
    description: Optional[str] = None
    url: Optional[str] = None

    @classmethod
    def from_json(cls, data: Dict[str, Any], database_dir: str) -> "DBConfig":
        """ Reconstructs a DBConfig instance from a JSON representation.
            Requires the current config database directory to handle replacments.
        """
        data["mode"] = Mode.from_string(data["mode"])
        data["path"] = data["path"].replace("$datadir", database_dir)
        return cls(**data)

    @classmethod
    def from_file(cls, filename: str, database_dir: str) -> "DBConfig":
        """ Constructs a DBConfig instance from a file containing JSON.
            Requires the current config database directory to handle replacments.
        """
        with open(filename, encoding="utf-8") as handle:
            raw = json.load(handle)
        return cls.from_json(raw, database_dir)


SubComponents = Dict[Any, int]


@dataclasses.dataclass(frozen=True)
class Components:
    """ Contains the various components used throughout the component scoring
        section of the module
    """
    nrps: SubComponents
    pks: SubComponents
    secmet: SubComponents
    functions: SubComponents


@dataclasses.dataclass(frozen=True)
class Hit:
    """ A hit between a query CDS and a reference CDS """
    reference_record: str
    reference_id: str
    cds: CDSFeature
    percent_identity: float  # 0. to 100.
    blast_score: float
    percent_coverage: float  # 0. to 100.
    evalue: float

    def to_json(self) -> Dict[str, Any]:
        """ Returns a JSON representation of a Hit """
        return {
            "reference_record": self.reference_record,
            "reference_id": self.reference_id,
            "query_cds_name": self.cds.get_name(),
            "percent_identity": self.percent_identity,
            "blast_score": self.blast_score,
            "percent_coverage": self.percent_coverage,
            "evalue": self.evalue,
        }

    @classmethod
    def from_json(cls, data: Dict[str, Any], record: Record) -> "Hit":
        """ Reconstructs a Hit instance from a JSON representation

            Requires the matching Record instance in order to repopulate CDS links.
        """
        copy = dict(data)
        copy["cds"] = record.get_cds_by_name(copy.pop("query_cds_name"))
        return cls(**copy)


class ReferenceCDS:
    """ A CDS within a reference area """
    def __init__(self, name: str, function: str, components: Dict[str, List[Any]],
                 location: Location) -> None:
        self.name = name
        self.function = function
        self.components = components
        self.location = location

    def overlaps_with(self, area: "ReferenceArea") -> bool:
        """ Returns True if this CDS overlaps with the given ReferenceArea """
        return self.location.end > area.start and self.location.start < area.end

    @classmethod
    def from_json(cls, name: str, data: Dict[str, Any]) -> "ReferenceCDS":
        """ Reconstructs a ReferenceCDS from a JSON representation """
        return cls(name, data["function"], data["components"], location_from_string(data["location"]))

    def to_json(self) -> Dict[str, Any]:
        """ Constructs a JSON representation of a ReferenceCDS """
        return {
            "function": self.function,
            "components": self.components,
            "location": str(self.location),
        }

    def __str__(self) -> str:
        return f"ReferenceCDS({self.name}, {self.location})"

    def get_minimal_json(self) -> Dict[str, Any]:
        """ Returns the minimum information required as JSON for the purposes of
            drawing this CDS in HTML output
        """
        # to match IOrf in antismash-js for the purposes of drawing
        return {
            "locus_tag": self.name,
            "start": self.location.start,
            "end": self.location.end,
            "strand": self.location.strand,
            "function": self.function,
        }


class ReferenceArea:
    """ A generic area within a reference record, e.g. a region or protocluster """
    def __init__(self, accession: str, start: int, end: int, cds_mapping: Dict[str, str],
                 cdses: Dict[str, ReferenceCDS], products: List[str]) -> None:
        self.accession = accession
        self.start = start
        self.end = end
        self.cds_mapping = cds_mapping
        self.cdses = {name: cds for name, cds in cdses.items() if cds.overlaps_with(self)}
        self.products = products
        self._components: Optional[Components] = None

    def get_product_string(self) -> str:
        """ Returns a single string of the area's product(s) """
        return ", ".join(self.products)

    def get_identifier(self) -> str:
        """ Returns a user friendly identifier for the area within its parent
            record
        """
        return f"{self.accession}: {self.start}-{self.end}"

    def get_component_data(self) -> Optional[Components]:
        """ Returns the area's Component data """
        return self._components

    def set_component_data(self, data: Components) -> None:
        """ Sets the component data of the area to avoid (re)calculating when unnecessary """
        if not isinstance(data, Components):
            raise TypeError(f"component data expected to be 'Components', but was '{type(data)}'")
        self._components = data

    def to_json(self) -> Dict[str, Any]:
        """ Constructs a JSON representation of a ReferenceArea """
        return {
            "accession": self.accession,
            "cdses": {name: cds.to_json() for name, cds in self.cdses.items()},
            "cds_mapping": self.cds_mapping,
            "start": self.start,
            "end": self.end,
            "products": self.products,
        }


class ReferenceProtocluster(ReferenceArea):
    """ A single protocluster within a reference record """
    def __init__(self, accession: str, start: int, end: int, cds_mapping: Dict[str, str],
                 cdses: Dict[str, ReferenceCDS], cores: List[ReferenceCDS], product: str) -> None:
        super().__init__(accession, start, end, cds_mapping, cdses, [product])
        self.cores = cores
        self.product = product

    @classmethod
    def from_json(cls, accession: str, data: Dict[str, Any], cdses: Dict[str, ReferenceCDS],
                  cds_mapping: Dict[str, str]) -> "ReferenceProtocluster":
        """ Reconstructs a ReferenceProtocluster from a JSON representation """
        cores = [cdses[core] for core in data["core_cdses"]]
        location = location_from_string(data["location"])
        return cls(accession, location.start, location.end, cds_mapping, cdses, cores, data["product"])

    def to_json(self) -> Dict[str, Any]:
        return {
            "core_cdses": [cds.name for cds in self.cores],
            "product": self.product,
            "location": str(FeatureLocation(self.start, self.end)),
        }

    def __str__(self) -> str:
        return f"ReferenceProtocluster({self.start}, {self.end}, {self.product})"

    def __repr__(self) -> str:
        return str(self)


class ReferenceRegion(ReferenceArea):
    """ A single region within a reference record """
    def __init__(self, accession: str, start: int, end: int, protoclusters: List[ReferenceProtocluster],
                 cdses: Dict[str, ReferenceCDS], products: List[str], cds_mapping: Dict[str, str],
                 description: str, organism: str) -> None:
        super().__init__(accession, start, end, cds_mapping, cdses, products)
        self.protoclusters = protoclusters
        self.description = description
        self.organism = organism

    @classmethod
    def from_json(cls, accession: str, data: Dict[str, Any], cds_mapping: Dict[str, str]) -> "ReferenceRegion":
        """ Reconstructs a ReferenceRegion from a JSON representation """
        cdses = {name: ReferenceCDS.from_json(name, cds) for name, cds in data["cdses"].items()}

        return cls(
            accession,
            data["start"],
            data["end"],
            [ReferenceProtocluster.from_json(accession, proto, cdses, cds_mapping) for proto in data["protoclusters"]],
            cdses,
            data["products"],
            cds_mapping,
            data.get("description", ""),
            data.get("organism", ""),
        )

    def to_json(self) -> Dict[str, Any]:
        result = super().to_json()
        if self.organism:
            result["organism"] = self.organism
        if self.description:
            result["description"] = self.description
        result["protoclusters"] = [proto.to_json() for proto in self.protoclusters]
        return result

    def __repr__(self) -> str:
        return f"ReferenceRegion({self.accession}, {self.start}, {self.end}, {self.products})"


class ReferenceRecord:
    """ A complete reference record, for the purposes of tracking its various regions
        and vastly simplifying member CDS referencing in the database's matching
        FASTA file
    """
    def __init__(self, accession: str, regions: List[ReferenceRegion], cds_mapping: Dict[str, str]) -> None:
        self.accession = accession
        self.regions = regions
        self.cds_mapping = cds_mapping

    @classmethod
    def from_json(cls, accession: str, data: Dict[str, Any]) -> "ReferenceRecord":
        """ Reconstructs a ReferenceRecord from a JSON representation """
        regions = [ReferenceRegion.from_json(accession, region, data["cds_mapping"]) for region in data["regions"]]
        return ReferenceRecord(accession, regions, data["cds_mapping"])


# for caching database metadata
_LOADED_DATA: Dict[str, Dict[str, ReferenceRecord]] = {}


def load_data(filename: str) -> Dict[str, ReferenceRecord]:
    """ Loads a database from file, caching the results to avoid further
        loads when not necessary

        Arguments:
            filename: the path of the file to load

        Returns:
            a dictionary mapping reference accession to ReferenceRecord instance

    """
    if filename not in _LOADED_DATA:
        with open(filename, encoding="utf-8") as handle:
            raw = json.loads(handle.read())
        result = {}
        for accession, record in raw.items():
            result[accession] = ReferenceRecord.from_json(accession, record)
        _LOADED_DATA[filename] = result
    return _LOADED_DATA[filename]


class ReferenceScorer:
    """ Keeps track of scoring information and hits for a reference area """
    def __init__(self, best_hits: Dict[str, Hit], reference: ReferenceRegion,
                 identity: float, order: float, component: Optional[float],
                 mode: Mode, final_score: Optional[float] = None) -> None:
        self.accession = reference.accession
        self.hits_by_gene = best_hits
        for hit in best_hits.values():
            assert hit.reference_record == reference.accession, f"{hit.reference_record} != {reference.accession}"
        self.identity = identity
        self.order = order
        self.component = component
        self.reference = reference
        self._final_score: Optional[float] = final_score
        self.mode = mode

    def to_json(self) -> Dict[str, Any]:
        """ Constructs a JSON representation of a ReferenceScorer """
        return {
            "identity": self.identity,
            "order": self.order,
            "component": self.component,
            "ref_id": self.reference.get_identifier(),
            "final_score": self._final_score,
            "mode": str(self.mode),
            "hits": {name: hit.to_json() for name, hit in self.hits_by_gene.items()},
        }

    @classmethod
    def from_json(cls, data: Dict[str, Any], record: Record,
                  references: Dict[str, ReferenceRegion]) -> "ReferenceScorer":
        """ Regenerate a ReferenceScorer instance from a JSON representation """
        hits = {name: Hit.from_json(hit, record) for name, hit in data["hits"].items()}
        ref = references[data["ref_id"]]
        return ReferenceScorer(hits, ref, data["identity"], data["order"],
                               data["component"], Mode.from_string(data["mode"]),
                               final_score=data["final_score"])

    @property
    def final_score(self) -> float:
        """ Calculates the singular score of a result from its metric scores """
        if self._final_score is None:
            metrics = [self.identity]
            if len(self.hits_by_gene) > 1:
                metrics.append(self.order)
            if self.component is not None:
                metrics.append(self.component)
            score = 1.
            for metric in metrics:
                score *= metric
            self._final_score = score ** (1/len(metrics))
        return self._final_score

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        order = f"{self.order:.2f}" if len(self.hits_by_gene) > 1 else "N/A"
        components = f"{self.component:.2f}" if self.component is not None else "N/A"
        return f"ReferenceScorer({self.accession}: raw_id={self.identity:.2f}, order={order}, comp={components})"

    def table_string(self) -> str:
        """ Generates a string suitable for a tooltip of a particular result """
        order = f"{self.order:.2f}" if len(self.hits_by_gene) > 1 else "N/A"
        components = f"{self.component:.2f}" if self.component is not None else "N/A"
        return f"{self.final_score:.2f} (id:{self.identity:.2f}, order:{order}, components:{components})"

    def __lt__(self, other: Any) -> bool:
        if not isinstance(other, ReferenceScorer):
            raise TypeError(f"cannot compare ReferenceScorer to {type(other)}")
        return self.final_score < other.final_score


HitsByReference = Dict[ReferenceRegion, Dict[str, List[Hit]]]
ScoresByRegion = List[Tuple[ReferenceRegion, float]]
