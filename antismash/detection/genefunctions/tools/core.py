# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Core functions and classes for detection gene functions """

from abc import ABC as AbstractBaseClass, abstractmethod
from dataclasses import dataclass, field
import logging
from typing import Any, Callable, Dict, Generic, Iterable, Mapping, Optional, Self, TypeVar

from antismash.common import fasta, subprocessing, utils
from antismash.common.hmmscan_refinement import refine_hmmscan_results
from antismash.common.html_renderer import Markup
from antismash.common.json import JSONBase
from antismash.common.secmet import CDSFeature, Record
from antismash.common.secmet.qualifiers import ECGroup, GeneFunction
from antismash.config import ConfigType
from antismash.config.args import ModuleArgs

T = TypeVar("T", bound="Hit")


@dataclass(kw_only=True, slots=True)
class Hit(JSONBase):
    """ A simple container for a particular query and reference pair, with details
        of subfunctions, if any are known, along with optional descriptions.
    """
    query_id: str
    reference_id: str
    subfunctions: list[str] = field(default_factory=list)
    description: str = ""

    def get_full_description(self) -> str:
        """ Constructs a more detailed description of the hit """
        if self.description:
            return f"{self.reference_id}: {self.description}"
        return self.reference_id

    def get_html_fragment(self, metadata: dict[str, Any] = None) -> Markup:  # pylint: disable=unused-argument
        """ Constructs a small HTML-compatible fragment to describe the hit """
        return Markup(f"{self.query_id}: {self.get_full_description()}")

    @classmethod
    def from_json(cls, data: dict[str, Any]) -> Self:
        """ Reconstructs an instance from the given data """
        return cls(**data)


@dataclass(kw_only=True, slots=True)
class HMMHit(Hit):
    """ A variant of Hit, specific to HMMer profile hits, containing the basic
        information typically contained in HMM results
    """

    bitscore: float
    evalue: float
    query_start: Optional[int] = None
    query_end: Optional[int] = None

    @classmethod
    def from_json(cls, data: dict[str, Any]) -> Self:
        """ Reconstructs an instance from the given data """
        return cls(**data)


class _ResultBase(AbstractBaseClass, Generic[T]):
    @staticmethod
    @abstractmethod
    def regenerate_hits(hits_by_name: dict[str, dict[str, Any]]) -> dict[str, T]:
        """ Regenerates the hits contained by this generic, as instances of
            the classes contained types
        """
        raise NotImplementedError

    @abstractmethod
    def add_to_record(self, record: Record) -> None:
        """ Adds relevant contents of the results to the given record """
        raise NotImplementedError

    @abstractmethod
    def build_html_fragments(self) -> list[Markup]:
        """ Builds a list of HTML fragments, one for each contained hit """
        raise NotImplementedError

    @classmethod
    @abstractmethod
    def from_json(cls, data: Dict[str, Any]) -> Optional[Self]:
        """ Reconstructs an instance from the given data """
        raise NotImplementedError

    @abstractmethod
    def to_json(self) -> dict[str, Any]:
        """ Converts the instance into a JSON-compatible format """
        raise NotImplementedError


class FunctionResults(_ResultBase[T]):
    """ A container for functions detected by a particular tool """
    def __init__(self, *, tool: str, best_hits: dict[str, T],
                 function_mapping: Mapping[str, GeneFunction],
                 group_mapping: Mapping[str, list[ECGroup]] | None = None,
                 subfunction_mapping: Mapping[str, list[str]] | None = None,
                 ) -> None:
        self.tool = tool
        self.best_hits = best_hits
        self.function_mapping = function_mapping
        self.group_mapping: Mapping[str, list[ECGroup]] = group_mapping or {}
        self.subfunction_mapping: Mapping[str, list[str]] = subfunction_mapping or {}

    def add_to_record(self, record: Record) -> None:
        """ Annotate resistance genes in CDS features """
        logging.debug("annotating CDS features with %s info: %d CDSes",
                      self.tool, len(self.best_hits))
        for feature_name, result in self.best_hits.items():
            function: GeneFunction = self.function_mapping[feature_name]
            feature = record.get_cds_by_name(feature_name)
            feature.gene_functions.add(function, self.tool, result.get_full_description())

    def build_html_fragments(self) -> list[Markup]:
        return [hit.get_html_fragment() for hit in self.best_hits.values()]

    @classmethod
    def from_json(cls, data: Dict[str, Any]) -> Optional[Self]:
        mapping = {cds_name: GeneFunction.from_string(function)
                   for cds_name, function in data.pop("function_mapping").items()}
        hits = cls.regenerate_hits(data.pop("best_hits", {}))
        result = cls(tool=data.pop("tool"), best_hits=hits, function_mapping=mapping, **data)
        return result

    def to_json(self) -> Dict[str, Any]:
        return {
            "tool": self.tool,
            "best_hits": dict(self.best_hits.items()),
            "function_mapping": {name: str(function) for name, function in self.function_mapping.items()},
            "group_mapping": dict(self.group_mapping),
            "subfunction_mapping": dict(self.subfunction_mapping),
        }


class SimpleFunctionResults(FunctionResults[Hit]):
    """ A results container for tools with minimally detailed hits """
    @staticmethod
    def regenerate_hits(hits_by_name: dict[str, dict[str, Any]]) -> dict[str, Hit]:
        return {name: Hit.from_json(hit) for name, hit in hits_by_name.items()}


class HMMFunctionResults(FunctionResults[HMMHit]):
    """ A results container for tools based on HMMer hits """
    @staticmethod
    def regenerate_hits(hits_by_name: dict[str, dict[str, Any]]) -> dict[str, HMMHit]:
        return {name: HMMHit.from_json(hit) for name, hit in hits_by_name.items()}


def _dummy_add(_: ModuleArgs) -> None:
    """ A dummy function used as a default for Tool instances """


def _dummy_check(_: ConfigType) -> list[str]:
    """ A dummy function used as a default for Tool instances """
    return []


def _dummy_prepare(_: bool) -> list[str]:
    """ A dummy function used as a default for Tool instances """
    return []


@dataclass(kw_only=True)
class Tool:
    """ A container for a gene function detection tool's functionality

        name: the name of the tool
        classify: a function that, when called, provides the core information
                  about any gene functions the tool finds
        add_arguments: an optional function to add additional arguments to the general module arguments
        check_options: an optional function to check any runtime options specific to the tool
        prepare_data: an optional function to prepare any relevant data (e.g. pressing HMM profiles)
        check_prereqs: an optional function to check any prerequisites for runtime (e.g. required binaries)
    """
    name: str
    classify: Callable[[Iterable[CDSFeature], ConfigType], FunctionResults[Any]]
    add_arguments: Callable[[ModuleArgs], None] = _dummy_add
    check_options: Callable[[ConfigType], list[str]] = _dummy_check
    prepare_data: Callable[[bool], list[str]] = _dummy_prepare
    check_prereqs: Callable[[ConfigType], list[str]] = _dummy_check


def scan_profiles_for_functions(cds_features: Iterable[CDSFeature], database: str,
                                hmmscan_opts: Optional[list[str]] = None) -> dict[str, HMMHit]:
    """ Finds possible classifications for the provided genes.

        Arguments:
            cds_features: a list of CDSFeatures to classify
            database: the path to the database to check
            hmmscan_opts: a list of extra options to provide to hmmscan

        Returns:
            a dictionary mapping CDS name to a list of HMM hits
    """
    search_fasta = fasta.get_fasta_from_features(cds_features)
    results = subprocessing.run_hmmscan(database, search_fasta, hmmscan_opts)
    hmm_lengths = utils.get_hmm_lengths(database)
    hmm_results = refine_hmmscan_results(results, hmm_lengths)

    best_hits: Dict[str, HMMHit] = {}

    for cds in cds_features:
        cds_name = cds.get_name()
        hits = hmm_results.get(cds_name)
        if not hits:
            continue
        best = hits[0]
        best_hits[cds_name] = HMMHit(
            query_id=cds_name,
            reference_id=best.hit_id,
            bitscore=best.bitscore,
            evalue=best.evalue,
            query_start=best.query_start,
            query_end=best.query_end,
        )
    return best_hits
