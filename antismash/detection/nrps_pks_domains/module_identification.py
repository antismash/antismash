# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Functions to find module borders within a CDS and classes to represent them

Modules are expected to have certain components, and different configurations
have different requirements to be considered a complete module. Non-canonical
modules should never be considered complete.

Modules are generally expected to have the following domain layout ([] denotes optional):
[starter] loader [modification, ...] carrier_protein [finalisation]

TransAT PKS modules are special and can have the following layout:
starter(specifically Trans-AT-KS) [modification, ...] carrier_protein [modification(specifically KR)] [finalisation]

Stereoisomers (D-) will always be the first part if they are present, since it
seems to be most common in Norine's list of monomers.

Non-carbon methylations are represented as either OMe or NMe, Norine mostly uses
that form but occasionally has MeO.
"""

from typing import Any, Dict, Iterator, List
from typing import Optional  # comment hints, pylint: disable=unused-import

from antismash.common.hmmscan_refinement import HMMResult

ADENYLATIONS = {
    "AMP-binding",
    "A-OX",
}
ACYLTRANSFERASES = {
    "PKS_AT",
}
CONDENSATIONS = {
    "Cglyc",
    "Condensation_DCL",
    "Condensation_LCL",
    "Condensation_Starter",
    "Condensation_Dual",
    "Heterocyclization",
}
ENDS = {
    "Epimerization",
    "Thioesterase",
    "TD",
}
KETOSYNTHASES = {
    "PKS_KS",
}
MODIFIERS = {
    "PKS_DH",
    "PKS_DH2",
    "PKS_DHt",
    "PKS_KR",
    "PKS_ER",
    "cMT",
    "nMT",
    "oMT",
}
CARRIER_PROTEINS = {
    "ACP",
    "ACP_beta",
    "PCP",
    "PP-binding"
}
ALTERNATE_STARTERS = {
    "CAL_domain",
    "SAT",
}
IGNORED = {
    "ACPS",  # activates ACP domains, has no function without them
    "Aminotran_1_2", "Aminotran_3", "Aminotran_4", "Aminotran_5",  # no useful function yet
    "B",
    "ECH",
    "F",
    "FkbH",
    "GNAT",
    "Hal",
    "NAD_binding_4",
    "NRPS-COM_Cterm",  # external to modules, but important for polymers
    "NRPS-COM_Nterm",  # external to modules, but important for polymers
    "PKS_Docking_Cterm",  # external to modules, but important for polymers
    "PKS_Docking_Nterm",  # external to modules, but important for polymers
    "Polyketide_cyc", "Polyketide_cyc2",  # type-II PKS specific
    "PS",
    "TIGR01720",  # NRPS domain, between an Epimerase and the next Condensation
    "TIGR02353",  # NRPS terminal domain of unknown function
    "X",
}
SPECIAL = {
    "Trans-AT_docking",
}

CLASSIFICATIONS = {
    "A": ADENYLATIONS,
    "AT": ACYLTRANSFERASES,
    "C": CONDENSATIONS,
    "S": ALTERNATE_STARTERS,
    "E": ENDS,
    "KS": KETOSYNTHASES,
    "+":  MODIFIERS,
    "CP": CARRIER_PROTEINS,
    "!": SPECIAL,
    "ignore": IGNORED,
}


class IncompatibleComponentError(ValueError):
    """ For better communication when evaluating the addition of a component """
    pass


class Component:
    """ A component of a module, represents a single domain.
        A subtype can be optionally supplied to differentiate between
        types of domains (e.g. a Trans-AT variant of a KS domain)
    """
    def __init__(self, domain: HMMResult, subtype: str = "") -> None:
        self._domain = domain
        self.classification = classify(domain.hit_id)
        self.subtype = subtype
        assert self.classification

    @property
    def domain(self) -> HMMResult:
        """ The domain the component was built from """
        return self._domain

    @property
    def label(self) -> str:
        """ The name of the domain the component was built from, e.g. PKS_KS """
        return self._domain.hit_id

    def is_adenylation(self) -> bool:
        """ Returns True if the component can function as an adenylation domain """
        return self.label in ADENYLATIONS

    def is_acyltransferase(self) -> bool:
        """ Returns True if the component can function as an acyltransferase domain """
        return self.label in ACYLTRANSFERASES

    def is_condensation(self) -> bool:
        """ Returns True if the component can function as a condensation domain """
        return self.label in CONDENSATIONS

    def is_starter(self) -> bool:
        """ Returns True if the component can function as a starter domain """
        return any(self.label in collection for collection in (
            CONDENSATIONS,
            KETOSYNTHASES,
            ADENYLATIONS,
            ALTERNATE_STARTERS
        ))

    def is_loader(self) -> bool:
        """ Returns True if the component can function as a loader domain """
        return self.is_acyltransferase() or self.is_adenylation()

    def is_modification(self) -> bool:
        """ Returns True if the component can function as a modification domain """
        return self.label in MODIFIERS

    def is_carrier_protein(self) -> bool:
        """ Returns True if the component can function as a carrier protein (e.g. ACP)"""
        return self.label in CARRIER_PROTEINS

    def is_end(self) -> bool:
        """ Returns True if the component terminates the module """
        return self.label in ENDS

    def is_ignored(self) -> bool:
        """ Returns True if the component is ignored for the purposes of modules """
        return self.label in IGNORED

    def is_special(self) -> bool:
        """ Returns True if the component has a special function for specific modules """
        return self.label in SPECIAL

    def is_pks_specific(self) -> bool:
        """ Returns True if the component is specific to PKS modules """
        if self.label.startswith("PKS"):
            return True
        return self.label in ACYLTRANSFERASES or self.label in KETOSYNTHASES

    def is_nrps_specific(self) -> bool:
        """ Returns True if the component is specific to NRPS modules """
        return self.label in ADENYLATIONS or self.label in CONDENSATIONS

    def __str__(self) -> str:
        return self.classification + (("(%s)" % self.subtype) if self.subtype else "")

    def to_json(self) -> Dict[str, Any]:
        """ Generate a JSON representation of the component """
        result = {
            "domain": self._domain.to_json(),
        }  # type: Dict[str, Any]
        if self.subtype:
            result["subtype"] = self.subtype
        return result

    @classmethod
    def from_json(cls, data: Dict[str, Any]) -> "Component":
        """ Construct a component from a JSON representation """
        subtype = data.get("subtype", "")
        assert isinstance(subtype, str), subtype
        return cls(HMMResult.from_json(data["domain"]), subtype)


class Module:
    """ A class for constructing an NRPS/PKS module based on the various domains
        added.
    """
    def __init__(self) -> None:
        self._components = []  # type: List[Component]
        self._starter = None  # type: Optional[Component]
        self._loader = None  # type: Optional[Component]
        self._modifications = []  # type: List[Component]
        self._carrier_protein = None  # type: Optional[Component]
        self._end = None  # type: Optional[Component]
        self._others = []  # type: List[Component]

    def to_json(self) -> Dict[str, Any]:
        """ Generate a JSON representation of the module """
        return {"components": [comp.to_json() for comp in self._components]}

    @classmethod
    def from_json(cls, data: Dict[str, Any]) -> "Module":
        """ Construct a module from a JSON representation """
        module = cls()
        for value in data["components"]:
            module.add_component(Component.from_json(value))
        return module

    def is_pks(self) -> bool:
        """ Returns True if the module is a PKS module """
        return bool(self._starter and self._starter.is_pks_specific()
                    or self._loader and self._loader.is_pks_specific())

    def is_nrps(self) -> bool:
        """ Returns True if the module is a NRPS module """
        return bool(self._starter and self._starter.is_nrps_specific()
                    or self._loader and self._loader.is_nrps_specific())

    def is_trans_at(self) -> bool:
        """ Returns True if the module is Trans-AT variant of a PKS module """
        return bool(self._starter and self._starter.subtype == "Trans-AT-KS" and not self._loader)

    def is_iterative(self) -> bool:
        """ Returns True if the module is an iterative variant of a PKS module """
        return bool(self._starter and self._starter.subtype == "Iterative-KS")

    def ensure_suitable(self, component: Component) -> None:  # pylint: disable=too-many-branches
        """ Raises a ValueError if adding the given component would be an issue """
        if component.is_ignored() or component.is_special():
            return

        if self._end:
            raise IncompatibleComponentError("adding extra component after end")

        if component.is_starter() and not component.is_loader():
            if self._components:
                raise IncompatibleComponentError("starter after other components")

        elif component.is_loader():
            if self._loader:
                raise IncompatibleComponentError("duplicate loader")
            if self._starter:
                if self._starter.is_pks_specific() and component.is_nrps_specific():
                    raise IncompatibleComponentError("adding a NRPS loader to a PKS starter")
                if self._starter.is_nrps_specific() and component.is_pks_specific():
                    raise IncompatibleComponentError("adding a PKS loader to a NRPS starter")
            if self._end or self._carrier_protein or self._modifications:
                raise IncompatibleComponentError("bad ordering, loader after other non-starter components")

        elif component.is_modification():
            if self._end:
                raise IncompatibleComponentError("bad ordering, modification domain after end")
            if self._carrier_protein and not (self.is_trans_at() and component.domain.hit_id == "PKS_KR"):
                raise IncompatibleComponentError("bad ordering, modification domain after carrier protein")

        elif component.is_carrier_protein():
            if self._carrier_protein:
                raise IncompatibleComponentError("duplicate carrier protein")

        elif component.is_end():
            assert not self._end

        else:
            raise IncompatibleComponentError("unhandled %s" % component)

    def add_component(self, component: Component) -> None:
        """ Adds the given component to the module, raising an error if it
            would be invalid
        """
        if component.is_ignored():
            return
        try:
            self.ensure_suitable(component)
        except IncompatibleComponentError as err:
            raise IncompatibleComponentError("cannot add %s to %s: %s" % (component, self, str(err)))
        if component.is_starter() and not self._starter:
            assert not self._starter, self
            self._starter = component
            if component.is_loader():
                assert not self._loader, self
                self._loader = component
        elif component.is_loader():
            assert not self._loader, self
            self._loader = component
        elif component.is_modification():
            self._modifications.append(component)
        elif component.is_carrier_protein():
            assert not self._carrier_protein
            self._carrier_protein = component
        elif component.is_end():
            assert not self._end
            self._end = component
        else:
            self._others.append(component)

        self._components.append(component)

    def is_complete(self) -> bool:
        """ Returns True if the module is considered complete """
        if self._starter and self._loader and self._carrier_protein:
            return True
        return bool(not self._loader and self._starter and self._carrier_protein and self.is_trans_at())

    def is_terminated(self) -> bool:
        """ Returns True if the module has a finalising domain (e.g. an epimerase) """
        return self._end is not None

    def is_termination_module(self) -> bool:
        """ Returns True if the module has contains domain terminating the entire
            assembly line (e.g. a thioesterase)
        """
        return bool(self._end and self._end.label in ["Thioesterase", "TD"])

    def is_starter_module(self) -> bool:
        """ Returns True if the module has an explicit starter unit or has
            no condensation/KS domain
        """
        return bool(self._starter and self._starter.label in {"Condensation_Starter"}.union(ALTERNATE_STARTERS)
                    or self._starter and self._starter is self._loader)

    def is_empty(self) -> bool:
        """ Returns True if the module has no components """
        return not self._components

    def __str__(self) -> str:
        return "[%s]" % (",".join(str(component) for component in self._components))

    def __iter__(self) -> Iterator[Component]:
        for component in self._components:
            yield component

    @property
    def start(self) -> int:
        """ Returns the starting protein position of the first component in the module """
        assert self._components
        return self._components[0].domain.query_start

    @property
    def end(self) -> int:
        """ Returns the starting protein position of the first component in the module,
            excluding any product finalising domain (e.g. a thioesterase)
        """
        assert self._components
        if self._end and len(self._components) > 1 and self._end.label in ["TD", "Thioesterase"]:
            return self._components[-2].domain.query_end
        return self._components[-1].domain.query_end

    def get_monomer(self, base: str = "", fallback: bool = False) -> str:  # pylint: disable=too-many-branches
        """ Builds a monomer including modifications from the given base. """
        state = []  # type: List[str]
        if self.is_trans_at() and not base:
            base = "mal"

        if not base:
            return ""

        if self.is_pks():
            if any(mod.label == "PKS_KR" for mod in self._modifications):
                conversions = {"mal": "ohmal", "mmal": "ohmmal", "mxmal": "ohmxmal", "emal": "ohemal"}
                base = conversions.get(base, base)

            if any(mod.label in ["PKS_DH", "PKS_DH2", "PKS_DHt"] for mod in self._modifications):
                conversions = {"ohmal": "ccmal", "ohmmal": "ccmmal", "ohmxmal": "ccmxmal", "ohemal": "ccemal"}
                base = conversions.get(base, base)

            if any(mod.label == "PKS_ER" for mod in self._modifications):
                conversions = {"ccmal": "redmal", "ccmmal": "redmmal", "ccmxmal": "redmxmal", "ccemal": "redemal"}
                base = conversions.get(base, base)

        for mod in self._modifications:
            if mod.label == "nMT":
                state.append("NMe")
            elif mod.label == "cMT":
                state.append("Me")
            elif mod.label == "oMT":
                state.append("OMe")

        if base.endswith("mmal"):
            state.append("Me")
            base = base.replace("mmal", "mal", 1)

        if base in ["X", "pk"] and not fallback:
            base = "?"

        state.append(base)

        if self._end and self._end.label == "Epimerization":
            state.insert(0, "D")
        return "-".join(state)


def classify(profile_name: str) -> str:
    """ Classifies a profile name, raising an exception if classification not
        possible

        Arguments:
            profile_name: the name of the profile used to find a domain

        Returns:
            a classification
    """
    for key, val in CLASSIFICATIONS.items():
        if profile_name in val:
            return key
    raise ValueError("could not classify domain: %s" % profile_name)


def build_modules_for_cds(domains: List[HMMResult], ks_subtypes: List[HMMResult]) -> List[Module]:
    """ Constructs a list of modules for a CDS based on the domains provided

        Arguments:
            domains: a list of HMMResults, one for each domain found
            ks_subtypes: a list of HMMResults, one for each PKS_KS domain given in domains

        Returns:
            a list of modules
    """
    domains = sorted(domains, key=lambda x: x.query_start)
    modules = [Module()]
    subtypes = iter(ks_subtypes)
    sub = ""
    for domain in domains:
        component = Component(domain)
        assert component.classification, "missing classification for %s" % domain.hit_id
        if component.classification == "KS":
            sub = next(subtypes).hit_id
        else:
            sub = ""
        component = Component(domain, sub)
        # start a new module if we have an explicit starter
        if component.is_starter() and not component.is_loader() and not modules[-1].is_empty():
            modules.append(Module())
        try:
            modules[-1].add_component(component)
        except IncompatibleComponentError:
            modules.append(Module())
            modules[-1].add_component(component)

    if modules[-1].is_empty():
        modules.pop()
    for module in modules:
        assert not module.is_empty()
    return modules
