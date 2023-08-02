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

from dataclasses import dataclass
from typing import Any, Dict, Iterator, List, Optional, Sequence, Tuple

from antismash.common.hmmscan_refinement import HMMResult
from antismash.common.secmet import CDSFeature

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
    "cAT",  # not a typical acyltransferase, closer to a thioesterase
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
    "Beta_elim_lyase",
    "LPG_synthase_C",
    "TauD",
}
CARRIER_PROTEINS = {
    "ACP",
    "ACP_beta",
    "PCP",
    "PKS_PP",
    "PP-binding"
}
ALTERNATE_STARTERS = {
    "CAL_domain",
    "SAT",
}
NON_MODULE = {  # external to modules, but important to find
    "NRPS-COM_Cterm",
    "NRPS-COM_Nterm",
    "PKS_Docking_Cterm",
    "PKS_Docking_Nterm",
}
OTHER = {
    "ACPS",  # activates ACP domains, has no function without them
    "Aminotran_1_2", "Aminotran_3", "Aminotran_4", "Aminotran_5",  # no useful function yet
    "B",
    "ECH",
    "F",
    "FkbH",
    "GNAT",
    "Hal",
    "NAD_binding_4",
    "Polyketide_cyc", "Polyketide_cyc2",  # type-II PKS specific
    "PS",
    "PT",  # fungal nonreducing PKS product template domain
    "TIGR02353",  # NRPS terminal domain of unknown function
    "X",
}
SPECIAL = {
    "Trans-AT_docking",
    "TIGR01720",  # NRPS domain, between an Epimerase and the next Condensation
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
    ".": OTHER,
    "ignore": NON_MODULE,
}

DOUBLE_TRANSPORTER_CASES = {
    ("LPG_synthase_C", "Beta_elim_lyase"),  # e.g. AF484556.1 in LmnJ, doi:10.1038/s41467-021-25798-8
}


class IncompatibleComponentError(ValueError):
    """ For better communication when evaluating the addition of a component """
    pass  # pylint: disable=unnecessary-pass


class Component:
    """ A component of a module, represents a single domain. """
    def __init__(self, domain: HMMResult, cds_name: str) -> None:
        self._domain = domain
        self.classification = classify(domain.hit_id)
        assert cds_name
        self.locus = cds_name
        assert self.classification

    @property
    def domain(self) -> HMMResult:
        """ The domain the component was built from """
        return self._domain

    @property
    def label(self) -> str:
        """ The name of the domain the component was built from, e.g. PKS_KS """
        return self._domain.hit_id

    @property
    def subtype(self) -> Optional[str]:
        """ The subtype of the domain, if any """
        if len(self.domain.detailed_names) < 2:
            return None
        return self.domain.detailed_names[1]

    @property
    def subtypes(self) -> List[str]:
        """ All subtypes of the domain, if any """
        return self.domain.detailed_names[1:]

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
            ACYLTRANSFERASES,
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
        return self.label in NON_MODULE

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
        return self.classification + "".join(f"({sub})" for sub in self.subtypes)

    def to_json(self) -> Dict[str, Any]:
        """ Generate a JSON representation of the component """
        result: Dict[str, Any] = {
            "domain": self._domain.to_json(),
            "locus": self.locus,
        }
        return result

    @classmethod
    def from_json(cls, data: Dict[str, Any]) -> "Component":
        """ Construct a component from a JSON representation """
        return cls(HMMResult.from_json(data["domain"]), data["locus"])


class Module:
    """ A class for constructing an NRPS/PKS module based on the various domains
        added.
    """
    def __init__(self, first_in_cds: bool = False) -> None:
        self._components: List[Component] = []
        self._starter: Optional[Component] = None
        self._loader: Optional[Component] = None
        self._modifications: List[Component] = []
        self._carrier_protein: Optional[Component] = None
        self._end: Optional[Component] = None
        self._others: List[Component] = []
        self._first_in_cds = first_in_cds
        self._unambiguous_accept = 0  # handles lookahead acceptance

    def to_json(self) -> Dict[str, Any]:
        """ Generate a JSON representation of the module """
        return {
            "components": [comp.to_json() for comp in self._components],
            "first_in_cds": self._first_in_cds,
        }

    @classmethod
    def from_json(cls, data: Dict[str, Any]) -> "Module":
        """ Construct a module from a JSON representation """
        module = cls(data.get("first_in_cds", True))  # default to true for backwards-compatibility
        components = [Component.from_json(comp) for comp in data["components"]]
        for i, component in enumerate(components):
            module.add_component(component, components[i + 1:])
        return module

    def is_pks(self) -> bool:
        """ Returns True if the module is a PKS module """
        return bool(any((comp.is_pks_specific() for comp in self._components)))

    def is_nrps(self) -> bool:
        """ Returns True if the module is a NRPS module """
        return bool(self._starter and self._starter.is_nrps_specific()
                    or self._loader and self._loader.is_nrps_specific())

    def is_trans_at(self) -> bool:
        """ Returns True if the module is Trans-AT variant of a PKS module """
        # since there's some alternatives, start by checking the bare minimum
        if not (self.is_pks() and self._starter and not self._loader):
            return False
        # if the KS is specifically Trans-AT, that's good enough
        if self._starter.subtype == "Trans-AT-KS":
            return True
        # otherwise, since the KS subtype may not be accurate enough, look for an ATd
        return any(comp.domain.hit_id == "Trans-AT_docking" for comp in self._others)

    def is_iterative(self) -> bool:
        """ Returns True if the module is an iterative variant of a PKS module """
        return bool(self._starter and self._starter.subtype == "Iterative-KS")

    def ensure_suitable(self, component: Component, lookahead: Sequence[Component]) -> None:  # pylint: disable=too-many-branches
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
                upcoming = [comp.domain.hit_id for comp in lookahead]
                valid = False
                for case in DOUBLE_TRANSPORTER_CASES:
                    if list(case) == upcoming[:len(case)]:
                        valid = True
                        break
                if not valid:
                    raise IncompatibleComponentError("duplicate carrier protein")

        elif component.is_end():
            assert not self._end

    def add_component(self, component: Component, lookahead: Sequence[Component]) -> None:
        """ Adds the given component to the module, raising an error if it
            would be invalid
        """
        assert component.classification in CLASSIFICATIONS, f"invalid classification: {component.classification}"
        if component.is_ignored():
            return
        if self._unambiguous_accept > 0:
            self._unambiguous_accept -= 1
        else:
            try:
                self.ensure_suitable(component, lookahead)
            except IncompatibleComponentError as err:
                raise IncompatibleComponentError(f"cannot add {component} to {self}: {err}")
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
            if not self._carrier_protein:
                self._carrier_protein = component
            else:
                upcoming = [comp.domain.hit_id for comp in lookahead]
                longest = 0
                for case in DOUBLE_TRANSPORTER_CASES:
                    if list(case) == upcoming[:len(case)]:
                        longest = len(case)
                if longest > 0:
                    self._unambiguous_accept = longest
                    self._others.append(component)
        elif component.is_end():
            assert not self._end
            self._end = component
        else:
            self._others.append(component)

        self._components.append(component)

    def is_complete(self) -> bool:
        """ Returns True if the module is considered complete """
        # if missing a proper starter and it would be ok as the first module, ensure it's the first module
        if self._starter and self._starter is self._loader and not self._first_in_cds:
            return False
        # otherwise, if it has the three vital parts
        if self._starter and self._loader and self._carrier_protein:
            return True
        # lastly, if it's transAT and has a carrier protein, it's ok
        return bool(self.is_trans_at() and self._carrier_protein)

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
                    # cannot have loader-only module if it's not the first module in the CDS
                    or self._starter and self._starter is self._loader and self._first_in_cds)

    def is_empty(self) -> bool:
        """ Returns True if the module has no components """
        return not self._components

    def __str__(self) -> str:
        core = ",".join(str(component) for component in self._components)
        return f"[{core}]"

    def __iter__(self) -> Iterator[Component]:
        for component in self._components:
            yield component

    @property
    def components(self) -> Tuple[Component, ...]:
        """ Returns the components of the module """
        return tuple(self._components)

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
        state: List[str] = []
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
    raise ValueError(f"could not classify domain: {profile_name}")


def build_modules_for_cds(domains: List[HMMResult], cds_name: str) -> List[Module]:
    """ Constructs a list of modules for a CDS based on the domains provided

        Arguments:
            domains: a list of HMMResults, one for each domain found
            cds_name: the name of the CDS feature the domains were found in

        Returns:
            a list of modules
    """
    domains = sorted(domains, key=lambda x: x.query_start)
    modules = [Module(first_in_cds=True)]
    components = [Component(domain, cds_name) for domain in domains]
    for i, component in enumerate(components):
        assert component.classification, f"missing classification for {component.domain.hit_id}"
        component = Component(component.domain, cds_name)
        # start a new module if we have an explicit starter
        if component.is_starter() and not component.is_loader() and not modules[-1].is_empty():
            modules.append(Module())
        try:
            modules[-1].add_component(component, components[i+1:i+3])
        except IncompatibleComponentError:
            modules.append(Module())
            modules[-1].add_component(component, [])

    if modules[-1].is_empty():
        modules.pop()
    for module in modules:
        assert not module.is_empty()
    return modules


@dataclass
class CDSModuleInfo:
    """ Used for bundling relevant details together about a specific CDS and
        modules within it
    """
    cds: CDSFeature
    modules: List[Module]


def combine_modules(current: CDSModuleInfo, previous: CDSModuleInfo) -> Optional[Module]:
    """ Combines trailing/leading incomplete modules of sequential CDS features
        into single modules, provided the resulting module would be valid and the
        features are on the same strand.

        The module following the leading module of the next CDS can also be
        merged, providing it would normally form a single module with the merged
        module. E.g. a post-CP KR domain in trans-AT PKS modules.

        This function modifies the given objects lists of modules.

        Arguments:
            current: the CDSModuleInfo instance of the earlier CDS on the sequence
            previous: the CDSModuleInfo instance of the later CDS on the sequence

        Returns:
            the newly merged Module or None
    """
    # both strands must match
    if current.cds.location.strand != previous.cds.location.strand:
        return None
    # and both must have modules
    if not current.modules or not previous.modules:
        return None
    # ensure the two args are ordered as expected
    previous, current = sorted([current, previous], key=lambda info: info.cds,
                               reverse=current.cds.location.strand == -1)
    head = previous.modules[-1]
    tail = current.modules[0]
    # both modules must be incomplete
    if head.is_complete() or tail.is_complete():
        return None
    # and avoid creating hybrid modules
    if head.is_pks() and tail.is_nrps() or head.is_nrps() and tail.is_pks():
        return None
    module = Module()
    for i, component in enumerate(head):
        module.add_component(component, head.components[i + 1:])
    try:
        for i, component in enumerate(tail):
            module.add_component(component, tail.components[i + 1:])
    except IncompatibleComponentError:
        return None
    # if it's still incomplete after the merge, discard it
    if not module.is_complete():
        return None

    # finally, replace the existing leading fragment with the new full module
    # and remove the tail fragment
    previous.modules[-1] = module
    current.modules.pop(0)

    # since it's reached this far, check if the next module in the latter CDS
    # might merge in too (e.g. trailing KR domain for a trans-AT PKS module that
    # didn't originally get detected as such due to the split)
    if not current.modules:
        return module

    next_module = current.modules[0]
    # is it a KR following a split trans_AT module?
    # e.g. AM746336 (in both kirAII and kirAV) and AF484556.1 (in LnmJ)
    if module.is_trans_at() and len(next_module.components) == 1 \
            and next_module.components[0].domain.hit_id == "PKS_KR":
        module.add_component(next_module.components[0], [])
        current.modules.pop(0)

    return module
