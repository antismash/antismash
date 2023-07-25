# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" NRPS/PKS Module features """

from enum import Enum, unique
from typing import Any, Dict, List, Optional, Tuple, Type, TypeVar

from .antismash_domain import AntismashDomain
from .feature import Feature, FeatureLocation, SeqFeature

T = TypeVar("T", bound="Module")


@unique
class ModuleType(Enum):
    """ An Enum representing the type of a Module.
        Allows for more flexible conversion and more robust value constraints.
    """
    UNKNOWN = 0
    NRPS = 1
    PKS = 2

    def __str__(self) -> str:
        return str(self.name).lower()

    @staticmethod
    def from_string(label: str) -> "ModuleType":
        """ Converts a string to a ModuleType instance when possible.
            Raises an error if not possible.
        """
        for value in ModuleType:
            if str(value) == label:
                return value
        raise ValueError(f"unknown module type: {label}")


class Module(Feature):
    """ A feature containing one or more AntismashDomain features used to represent module-like
        structures such as polyketide/NRPS modules. While individual domains are contained by a CDS,
        modules may contain domains from more than one CDS, but this should only be used for
        situations like split-CDS trans-AT PKS modules.

        The location of Module covers the entire area of it's child domains.
    """
    __slots__ = ["_parent_cds_names", "_domains", "_substrate_monomer_pairs",
                 "_complete", "_is_starter", "_is_final", "_module_type", "_is_iterative",
                 ]
    types = ModuleType
    FEATURE_TYPE = "aSModule"

    def __init__(self, domains: List[AntismashDomain], module_type: ModuleType = ModuleType.UNKNOWN,
                 complete: bool = False, starter: bool = False, final: bool = False,
                 iterative: bool = False) -> None:
        if not domains:
            raise ValueError("at least one domain required in module")

        if len({domain.location.strand for domain in domains}) != 1:
            raise ValueError("domains within a module cannot be on different strands")

        # if the parent CDSes are on the reverse strand, the domains are in the wrong order
        # so ensure they're sorted here as they appear in the translation
        strand = domains[0].location.strand
        reverse = strand == -1
        domains = sorted(domains, key=lambda x: x.location.start, reverse=reverse)
        if reverse:
            location = FeatureLocation(domains[-1].location.start, domains[0].location.end, strand=strand)
        else:
            location = FeatureLocation(domains[0].location.start, domains[-1].location.end, strand=strand)
        super().__init__(location, self.FEATURE_TYPE, created_by_antismash=True)

        if not isinstance(module_type, ModuleType):
            raise TypeError(f"module_type must be a ModuleType instance, not {type(module_type)}")

        # parent CDS name order should match that of the domains provided
        self._parent_cds_names: List[str] = []
        for domain in domains:
            if domain.locus_tag not in self._parent_cds_names:
                self._parent_cds_names.append(domain.locus_tag)
        self._domains = domains
        self._substrate_monomer_pairs: List[Tuple[str, str]] = []
        self._complete = complete
        self._is_starter = starter
        self._is_final = final
        self._module_type = module_type
        self._is_iterative = iterative

    @property
    def domains(self) -> Tuple[AntismashDomain, ...]:
        """ Returns a sequence containing all the domains within the module """
        return tuple(self._domains)

    @property
    def module_type(self) -> ModuleType:
        """ Returns the type of the module, e.g. NRPS, PKS """
        return self._module_type

    @property
    def monomers(self) -> Tuple[Tuple[str, str], ...]:
        """ Returns a sequence of substrate and monomer pairs, with the monomer
            being the substrate as modified by this module. These values are
            added by use of the `Module.add_monomer()` method.
        """
        return tuple(self._substrate_monomer_pairs)

    @property
    def parent_cds_names(self) -> Tuple[str, ...]:
        """ Returns a sequence of parent CDS names as embedded in the domains
            provided at time of construction. The sequence is sorted alphabetically.
        """
        return tuple(self._parent_cds_names)

    @property
    def protein_location(self) -> FeatureLocation:
        """ Returns the protein location of the module within a CDS. If the module
            contains domains from multiple modules, a ValueError is raised
        """
        if self.is_multigene_module():
            raise ValueError("cannot generate protein location for multi-CDS module")
        return FeatureLocation(self._domains[0].protein_location.start, self._domains[-1].protein_location.end)

    def get_parent_protein_location(self, parent: str) -> FeatureLocation:
        """ Returns the location within the specified parent for multi-CDS modules """
        if parent not in self._parent_cds_names:
            raise ValueError(f"module {self} has no parent named {parent}")
        domains = [domain for domain in self._domains if domain.locus_tag == parent]
        return FeatureLocation(domains[0].protein_location.start, domains[-1].protein_location.end)

    def add_monomer(self, substrate: str, monomer: str) -> None:
        """ Adds a substrate and the monomer produced by this module with that
            substrate
        """
        if not substrate:
            raise ValueError("substrate is required")
        if not monomer:
            raise ValueError("monomer is required")
        self._substrate_monomer_pairs.append((substrate, monomer))

    def is_complete(self) -> bool:
        """ Returns True if the module is considered complete """
        return self._complete

    def is_starter_module(self) -> bool:
        """ Returns True if the module is an explicit starter module for multi-module
            constructs
        """
        return self._is_starter

    def is_final_module(self) -> bool:
        """ Returns True if the module is an explicit termination module for multi-module
            constructs
        """
        return self._is_final

    def is_iterative(self) -> bool:
        """ Returns True if the module is considered iterative """
        return self._is_iterative

    def is_multigene_module(self) -> bool:
        """ Returns True if the module contains domains from multiple CDS features """
        return len(self._parent_cds_names) > 1

    def get_substrate_monomer_pairs(self) -> Tuple[Tuple[str, str], ...]:
        """ Returns the substrate/monomer pairings as a tuple of tuples """
        return tuple(self._substrate_monomer_pairs)

    def to_biopython(self, qualifiers: Dict[str, Any] = None) -> List[SeqFeature]:
        new: Dict[str, Optional[List[str]]] = {
            "domains": [domain.get_name() for domain in self.domains],
            "locus_tags": sorted(self._parent_cds_names),
            "type": [str(self.module_type)],
        }

        if self.is_complete():
            new["complete"] = None
        else:
            new["incomplete"] = None

        if self.is_starter_module():
            new["starter_module"] = None
        if self.is_final_module():
            new["final_module"] = None
        if self.is_iterative():
            new["iterative"] = None

        if self._substrate_monomer_pairs:
            new["monomer_pairings"] = [f"{sub} -> {mon}" for sub, mon in self._substrate_monomer_pairs]

        new.update(qualifiers or {})

        return super().to_biopython(new)

    @classmethod
    def from_biopython(cls: Type[T], bio_feature: SeqFeature, feature: T = None,
                       leftovers: Dict[str, List[str]] = None, record: Any = None) -> T:
        if leftovers is None:
            leftovers = Feature.make_qualifiers_copy(bio_feature)
        leftovers.pop("locus_tags", None)  # remove, it's generated data once in a record
        try:
            domain_names = leftovers.pop("domains")
            module_type = ModuleType.from_string(leftovers.pop("type")[0])
        except KeyError as err:
            raise ValueError(f"module at {bio_feature.location} missing expected qualifier: {err}")
        except ValueError as err:
            raise ValueError(f"module at {bio_feature.location}: {err}")

        complete = "complete" in leftovers
        if not complete and "incomplete" not in leftovers:
            raise ValueError(f"module at {bio_feature.location} missing one of: complete, incomplete")
        leftovers.pop("complete", "")
        leftovers.pop("incomplete", "")

        if not record:
            raise ValueError("record instance required for regenerating Module instance from biopython")
        # biopython parsing doesn't properly handle string crossing two lines due to length
        # despite writing them in the first place
        # so ensure that any inserted space is removed
        domain_names = [domain.replace(" ", "") for domain in domain_names]
        try:
            domains = [record.get_domain_by_name(domain) for domain in domain_names]
        except KeyError as err:
            raise ValueError(f"record does not contain domain referenced by module at {bio_feature.location}: {err}")

        starter = "starter_module" in leftovers
        if starter:
            leftovers.pop("starter_module")
        final = "final_module" in leftovers
        if final:
            leftovers.pop("final_module")
        iterative = "iterative" in leftovers
        if iterative:
            leftovers.pop("iterative")

        module = cls(domains, module_type, complete, starter, final, iterative)

        raw_monomers = leftovers.pop("monomer_pairings", [])
        for raw in raw_monomers:
            try:
                substrate, monomer = raw.split(" -> ")
            except ValueError:
                raise ValueError(f"module at {bio_feature.location} has invalid monomer pairing: {raw!r}")
            substrate = substrate.strip()
            monomer = monomer.strip()
            if not substrate or not monomer:
                raise ValueError(f"module at {bio_feature.location} has invalid monomer pairing: {raw!r}")
            module.add_monomer(substrate, monomer)

        return module

    def __repr__(self) -> str:
        return f"Module({self.location}, {self.module_type}: {[dom.domain for dom in self.domains]}"
