# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Contains the relevant functions for calculating similarity scores of the
    various components of two areas
"""

from collections import defaultdict
from typing import (
    Sequence,
    Optional,
)

from antismash.common.secmet import CDSFeature
from antismash.common.secmet.features.module import ModuleType
from antismash.common.secmet.qualifiers.gene_functions import GeneFunction

from .data_structures import ReferenceArea, Components, Mode, SubComponents


def calculate_component_score(query_components: Components, reference: ReferenceArea,
                              mode: Mode = Mode.BEST) -> Optional[float]:
    """ Calculates a similarity score (0..1) of the query components and reference
        area.

        For accurate results, both query and reference should have been annotated
        by the same version of antiSMASH.

        Arguments:
            query_components: a Components instance for the query area
            reference: the ReferenceArea to compare the query to
            mode: which mode to use, as query in reference and reference in query
                  have very different results

        Returns:
            a float ranging from 0 to 1, inclusive, or None if there are no
            meaningful components to compare
    """
    ref = gather_reference_components(reference)
    # if there are no non-function reference components, don't give a score at all
    if not any([ref.nrps, ref.pks, ref.secmet]):
        if mode != Mode.QUERY_IN_REFERENCE:
            return None
    score = compare(ref, query_components, mode)
    return score


def compare(ref: Components, query: Components, mode: Mode) -> float:
    """ Compares two Component instances and generates the score for the two

        Arguments:
            ref: the reference Components
            query: the query Components
            mode: the mode to use

        Returns:
            a float ranging from 0 to 1, inclusive
    """
    # compare all the component subsections
    nrps = compare_combos(ref.nrps, query.nrps, mode)
    pks = compare_combos(ref.pks, query.pks, mode)
    secmet = compare_combos(ref.secmet, query.secmet, mode)
    functions = compare_combos(ref.functions, query.functions, mode)

    assert 0 <= secmet <= 1, secmet
    assert 0 <= functions <= 1, functions

    # for modules, set the expected maxmium values and weightings based on mode
    if mode == Mode.REFERENCE_IN_QUERY:
        max_modules = sum(ref.pks.values()) + sum(ref.nrps.values())
        nrps_weighting = sum(ref.nrps.values()) / (max_modules or 1)
    elif mode == Mode.QUERY_IN_REFERENCE:
        max_modules = sum(query.pks.values()) + sum(query.nrps.values())
        nrps_weighting = sum(query.nrps.values()) / (max_modules or 1)
    else:
        # best of both, so take the minimum of each (will tend to poor results)
        max_modules = min(sum(query.pks.values()) + sum(query.nrps.values()),
                          sum(ref.pks.values()) + sum(ref.nrps.values()))
        nrps_weighting = min(sum(query.nrps.values()), sum(ref.nrps.values())) / (max_modules or 1)

    # if there's no modules, ignore them for the calculation
    if not max_modules:
        return (secmet + functions) / 2

    if 0.001 <= nrps_weighting <= 0.999:
        modules = nrps * nrps_weighting + pks * (1-nrps_weighting)
    elif nrps_weighting >= 0.999:
        modules = nrps
    else:  # some modules but no NRPS modules means PKS modules
        modules = pks
    assert 0 <= modules <= 1, modules
    return sum([modules, secmet, functions]) / 3


def compare_combos(ref: SubComponents, query: SubComponents, mode: Mode) -> float:
    """ Compares two SubComponents and scores them. The two should be of the
        same general kind if any meaningful score is to be calculated.

        Each individual comparison within the SubComponent uses exact matching.

        Arguments:
            ref: the reference SubComponents
            query: the query SubComponents
            mode: the comparison mode

        Returns:
            a float ranging from 0 to 1, inclusive
    """
    # if there are no components in either input, skip the set operations
    if not ref:
        if mode == Mode.QUERY_IN_REFERENCE:
            return 0.
        return 1.
    if not query:
        if mode == Mode.REFERENCE_IN_QUERY:
            return 0.
        return 1.

    # gather the baseline of what kind of combinations are present in the components
    ref_combos = set(ref)
    query_combos = set(query)

    # set the maximum possible counts by mode
    if mode == Mode.REFERENCE_IN_QUERY:
        max_possible = sum(ref.values())
    elif Mode.QUERY_IN_REFERENCE:
        max_possible = sum(query.values())
    else:
        max_possible = min(sum(query.values()), sum(ref.values()))
    assert max_possible

    # find how many of the maximum possible are present
    found = 0
    for combo in ref_combos.intersection(query_combos):
        # use the minimum to avoid overscoring a particular combination
        # e.g. 4 in query and 1 in reference in reference_in_query mode
        # should only count as matching the 1 from the reference
        found += min(query[combo], ref[combo])

    score = found / max_possible
    assert 0. <= score <= 1., score
    return score


def gather_reference_components(reference: ReferenceArea) -> Components:
    """ Build a Components instance from a ReferenceArea.

        If the ReferenceArea has an existing result, it will be used for
        performance reasons. If not, the result of the gather will be stored
        in the ReferenceArea after calculation.

        Arguments:
            reference: the ReferenceArea to gather components from

        Returns:
            a Components instance
    """
    # this function could be simplified/skipped if the counts were stored in the database
    # but it's much easier to tweak the calculation logic in code rather than having
    # to regenerate the database every time

    # check for any existing result first
    existing = reference.get_component_data()
    if existing is not None:
        return existing

    # NRPS modules
    nrps: SubComponents = defaultdict(int)
    # PKS modules
    pks: SubComponents = defaultdict(int)
    # domains as per hmm_detection's findings
    secmet: SubComponents = defaultdict(int)
    # gene functions as per hmm_detection and genefunctions detection modules
    functions: SubComponents = defaultdict(int)

    # check each CDS for relevant components
    for cds in reference.cdses.values():
        # tracking a genefunction of 'other' is only an estimate of size, so skip those
        if cds.function != "other":
            functions[cds.function] += 1

        if cds.components["secmet"]:
            for domain in cds.components["secmet"]:
                secmet[domain] += 1

        for module in cds.components["modules"]:
            # while incomplete modules might be useful,
            # they also contribute a great deal of noise
            if not module["complete"]:
                continue

            if module["type"] == "pks":
                target = pks
            elif module["type"] == "nrps":
                target = nrps
            else:
                raise ValueError(f"unknown module type: {module['type']}")
            # modules are exact full matches, so no modifications needed
            target[tuple(module["domains"])] += 1

    results = Components(nrps, pks, secmet, functions)
    # save the results to avoid future calculations for the same reference
    reference.set_component_data(results)
    return results


def gather_query_components(area_features: Sequence[CDSFeature]) -> Components:
    """ Build a components instance from a set of CDSFeatures.

        Arguments:
            area_features: the CDSFeatures to gather components from

        Returns:
            a Components instance
    """
    # NRPS modules
    nrps: SubComponents = defaultdict(int)
    # PKS modules
    pks: SubComponents = defaultdict(int)
    # domains as per hmm_detection's findings
    secmet: SubComponents = defaultdict(int)
    # gene functions as per hmm_detection and genefunctions detection modules
    functions: SubComponents = defaultdict(int)

    for cds in area_features:
        # tracking a genefunction of 'other' is only an estimate of size, so skip those
        if cds.gene_function != GeneFunction.OTHER:
            functions[str(cds.gene_function)] += 1

        if cds.sec_met:
            for domain_name in cds.sec_met.domain_ids:
                secmet[domain_name] += 1

        for module in cds.modules:
            if not module.is_complete():
                continue
            if module.module_type == ModuleType.PKS:
                target = pks
            elif module.module_type == ModuleType.NRPS:
                target = nrps
            else:
                raise ValueError(f"unknown module type: {module.module_type}")
            domains = []
            # tweak the domain list of a module
            # NOTE: it's important to do these same steps in the database generation
            for domain in module.domains:
                if domain.domain is None:
                    continue
                domain_name = domain.domain
                # as PKS_KS subtypes aren't stored, don't keep the Condensation subtypes
                if domain_name.startswith("Condensation_"):
                    domain_name = domain_name.split("_", 1)[0]
                domains.append(domain_name)
            target[tuple(domains)] += 1

    return Components(nrps, pks, secmet, functions)
