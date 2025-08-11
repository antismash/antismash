# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Results classes for the sideloader module
"""

from abc import ABC as AbstractBaseClass, abstractmethod
import dataclasses
from typing import Any, Dict, List, Union

from antismash.common.secmet import Protocluster, Record, SubRegion
from antismash.common.secmet.features.protocluster import SideloadedProtocluster
from antismash.common.secmet.features.subregion import SideloadedSubRegion
from antismash.common.secmet.locations import (
    CompoundLocation,
    FeatureLocation,
    Location,
)
from antismash.common.module_results import DetectionResults


@dataclasses.dataclass(frozen=True)
class SideloadSimple:
    """ A simple frozen container for the sideload-simple argument values """
    accession: str
    start: int
    end: int


def _qualifier_mapping(original: Dict[str, Any]) -> Dict[str, List[str]]:
    """ Converts a dictionary mapping strings to any type to
        a dictionary mapping strings to strings
    """
    result: Dict[str, List[str]] = {}
    for key, val in original.items():
        if isinstance(val, str):
            val = [val]
        else:
            val = list(map(str, val))
        result[key] = val
    return result


@dataclasses.dataclass(frozen=True, eq=True, order=True, unsafe_hash=True)
class Tool:
    """ A collection of details about the tool responsible for annotations """
    name: str
    version: str
    description: str
    configuration: Dict[str, List[str]]

    def __post_init__(self) -> None:
        for i in self.name:
            if not i.isalpha() and i not in "_- ":
                raise ValueError(f"invalid character in tool name: '{i}'")

    def to_json(self) -> Dict[str, Any]:
        """ Convert a Tool instance to JSON """
        return dataclasses.asdict(self)

    @classmethod
    def from_json(cls, raw: Dict[str, Any]) -> "Tool":
        """ Rebuild a Tool instance from JSON """
        return cls(
            str(raw["name"]),
            str(raw["version"]),
            str(raw.get("description", "")),
            _qualifier_mapping(raw.get("configuration", {})),
        )


class _AreaAnnotation(AbstractBaseClass):
    @abstractmethod
    def build_location(self) -> Location:
        """ Constructs a location object from the sideloaded details """
        raise NotImplementedError


class SubRegionAnnotation(_AreaAnnotation):
    """ A class for containing arbitrary data about a sideloaded subregion annotation """
    kind = "subregion"

    def __init__(self, start: int, end: int, label: str, tool: Tool, details: dict[str, list[str]],
                 *, circular_origin: int = None,
                 ) -> None:
        if not circular_origin and end <= start:
            raise ValueError(f"{label}: area end ({end:,}) must be greater than area start ({start:,})")
        if circular_origin:
            if circular_origin < 0:
                raise ValueError("circular origin cannot be negative")
            if start > circular_origin:
                raise ValueError("start cannot be later than the circular origin")
        self.circular_origin = circular_origin
        self.start = start
        self.end = end
        self.label = label
        self.details = details
        assert isinstance(tool, Tool)
        self.tool = tool

    def __len__(self) -> int:
        if self.circular_origin and self.start > self.end:
            return self.circular_origin - self.start + self.end
        return self.end - self.start

    def build_location(self) -> Location:
        """ Constructs a location object from the sideloaded details """
        if self.circular_origin and self.start > self.end:
            return CompoundLocation([
                FeatureLocation(self.start, self.circular_origin, 1),
                FeatureLocation(0, self.end, 1),
            ])
        return FeatureLocation(self.start, self.end)

    def to_secmet(self) -> SideloadedSubRegion:
        """ Constructs a SideloadedSubRegion instance from the annotation """
        return SideloadedSubRegion(self.build_location(), self.tool.name, label=self.label,
                                   extra_qualifiers=self.details)

    def to_json(self) -> Dict[str, Any]:
        """ Converts the annotation to JSON (separate from the JSON schema)"""
        result = dict(vars(self))
        result["tool"] = self.tool.to_json()
        return result

    @classmethod
    def from_json(cls, raw: dict[str, Any], circular_origin: int = None) -> "SubRegionAnnotation":
        """ Reconstructs an annotation from JSON (separate from the JSON schema)"""
        return cls(
            int(raw["start"]),
            int(raw["end"]),
            str(raw["label"]),
            Tool.from_json(raw["tool"]),
            _qualifier_mapping(raw.get("details", {})),
            circular_origin=circular_origin,
        )

    @classmethod
    def from_schema_json(cls, raw: dict[str, Any], tool: Tool, *, circular_origin: int = None,
                         ) -> "SubRegionAnnotation":
        """ Builds an instance from JSON matching the protocluster schema, along with a tool """
        raw["tool"] = tool.to_json()
        return cls.from_json(raw, circular_origin=circular_origin)

    def __eq__(self, other: Any) -> bool:
        if not isinstance(other, SubRegionAnnotation):
            return False
        return vars(self) == vars(other)

    def __str__(self) -> str:
        return f"Subregion({self.tool.name}, {self.start}-{self.end}, {self.label})"


class ProtoclusterAnnotation(_AreaAnnotation):
    """ A class for containing arbitrary data about a sideloaded protocluster annotation """
    kind = "protocluster"

    def __init__(self, core_start: int, core_end: int, product: str, tool: Tool, details: Dict[str, List[str]],
                 neighbourhood_left: int = 0, neighbourhood_right: int = 0,
                 *, circular_origin: int = None,
                 ) -> None:
        self.circular_origin = circular_origin
        if not circular_origin:
            if core_end <= core_start:
                raise ValueError(f"{product}: area end ({core_end:,}) must be greater than area start ({core_start:,})")
            if neighbourhood_left < 0 or neighbourhood_right < 0:
                raise ValueError("neighbourhoods must be annotated as absolute distances, not relative")
            if core_start - neighbourhood_left < 0:
                raise ValueError("neighbourhoods cannot extend out of a record")
        elif circular_origin < 0:
            raise ValueError("circular origin cannot be negative")
        elif core_start > circular_origin:
            raise ValueError("start cannot be later than the circular origin")
        self.core_start = core_start
        self.core_end = core_end
        self.product = product
        assert isinstance(tool, Tool)
        self.tool = tool
        self.details = details
        self.neighbourhood_left = neighbourhood_left
        self.neighbourhood_right = neighbourhood_right

    @property
    def start(self) -> int:
        """ The start of the protocluster, including neighbourhood """
        if self.circular_origin:
            return (self.core_start + self.circular_origin - self.neighbourhood_left) % self.circular_origin
        return self.core_start - self.neighbourhood_left

    @property
    def end(self) -> int:
        """ The end of the protocluster, including neighbourhood """
        if self.circular_origin:
            return (self.core_end + self.circular_origin + self.neighbourhood_right) % self.circular_origin
        return self.core_end + self.neighbourhood_right

    @property
    def label(self) -> str:
        """ The product of the protocluster, for compatibility purposes """
        return self.product

    def __len__(self) -> int:
        if self.circular_origin and self.end < self.start:
            return self.circular_origin - self.start + self.end
        return self.end - self.start

    def to_json(self) -> Dict[str, Any]:
        """ Converts the annotation to JSON (separate from the JSON schema)"""
        res = dict(vars(self))
        res["tool"] = self.tool.to_json()
        return res

    def build_core_location(self) -> Location:
        """ Constructs a location object from the sideloaded details for the protocluster's core """
        if self.circular_origin and self.core_start > self.core_end:
            return CompoundLocation([
                FeatureLocation(self.core_start, self.circular_origin, 1),
                FeatureLocation(0, self.core_end, 1),
            ])
        return FeatureLocation(self.core_start, self.core_end)

    def build_location(self) -> Location:
        """ Constructs a location object from the sideloaded details for the
            protocluster's neighbourhood
        """
        if self.circular_origin and self.start > self.end:
            return CompoundLocation([
                FeatureLocation(self.start, self.circular_origin, 1),
                FeatureLocation(0, self.end, 1),
            ])
        return FeatureLocation(self.start, self.end)

    def to_secmet(self) -> SideloadedProtocluster:
        """ Constructs a SideloadedProtocluster instance from the annotation """
        return SideloadedProtocluster(
            self.build_core_location(),
            self.build_location(),
            self.tool.name,
            self.product,
            neighbourhood_range=max(self.neighbourhood_left, self.neighbourhood_right),
        )

    @classmethod
    def from_json(cls, raw: dict[str, Any], *, circular_origin: int = None) -> "ProtoclusterAnnotation":
        """ Reconstructs an annotation from JSON (separate from the JSON schema)"""
        return cls(
            int(raw["core_start"]),
            int(raw["core_end"]),
            str(raw["product"]),
            Tool.from_json(raw["tool"]),
            _qualifier_mapping(raw.get("details", {})),
            int(raw.get("neighbourhood_left", 0)),
            int(raw.get("neighbourhood_right", 0)),
            circular_origin=circular_origin,
        )

    @classmethod
    def from_schema_json(cls, raw: Dict[str, Any], tool: Tool,
                         *, circular_origin: int = None,
                         ) -> "ProtoclusterAnnotation":
        """ Builds an instance from JSON matching the protocluster schema, along with a tool """
        raw["tool"] = tool.to_json()
        return cls.from_json(raw, circular_origin=circular_origin)

    def __eq__(self, other: Any) -> bool:
        if not isinstance(other, ProtoclusterAnnotation):
            return False
        return vars(self) == vars(other)

    def __str__(self) -> str:
        return f"Protocluster({self.tool.name}, {self.start}-{self.end}, {self.product})"


class SideloadedResults(DetectionResults):
    """ A wrapper allowing for different kinds of sideloaded annotations """
    schema_version = 1

    def __init__(self, record_id: str, subregions: List[SubRegionAnnotation],
                 protoclusters: List[ProtoclusterAnnotation]) -> None:
        super().__init__(record_id)
        self.subregions = subregions
        self.protoclusters = protoclusters

    def to_json(self) -> Dict[str, Any]:
        return {
            "record_id": self.record_id,
            "schema_version": self.schema_version,
            "protoclusters": [proto.to_json() for proto in self.protoclusters],
            "subregions": [sub.to_json() for sub in self.subregions],
        }

    @staticmethod
    def from_json(json: Dict[str, Any], record: Record) -> "SideloadedResults":
        if json["schema_version"] != SideloadedResults.schema_version:
            raise ValueError("Detection results have changed. No results can be reused.")

        assert json["record_id"] == record.id
        origin = len(record) if record.is_circular() else None
        subs = [SubRegionAnnotation.from_json(sub, circular_origin=origin) for sub in json["subregions"]]
        protos = [ProtoclusterAnnotation.from_json(proto, circular_origin=origin) for proto in json["protoclusters"]]
        return SideloadedResults(json["record_id"], subs, protos)

    def get_areas(self) -> List[Union[ProtoclusterAnnotation, SubRegionAnnotation]]:
        """ Returns a list of all sideloaded annotations """
        areas: List[Union[ProtoclusterAnnotation, SubRegionAnnotation]] = []
        areas.extend(self.subregions)
        areas.extend(self.protoclusters)
        # sort consistently with secmet's region.get_sideloaded_areas()
        return sorted(areas, key=lambda x: (x.start, -len(x), x.tool.name))

    def get_predicted_subregions(self) -> List[SubRegion]:
        return [sub.to_secmet() for sub in self.subregions]

    def get_predicted_protoclusters(self) -> List[Protocluster]:
        return [proto.to_secmet() for proto in self.protoclusters]

    def add_to_record(self, record: Record) -> None:
        for subregion in self.subregions:
            record.add_subregion(subregion.to_secmet())
        for protocluster in self.protoclusters:
            record.add_protocluster(protocluster.to_secmet())
