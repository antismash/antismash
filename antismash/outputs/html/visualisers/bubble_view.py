# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""
An aggregation module for HTML output of NRPS/PKS modules, pulling in information
from:
 - detection.nrps_pks_domains: for the modules and domains themselves
 - modules.nrps_pks: for gene order (when not co-linear) modules per candidate cluster
"""

from collections import defaultdict
from dataclasses import dataclass
from typing import Any, Dict, List, Set

from antismash.common import path
from antismash.common.html_renderer import (
    HTMLSections,
    FileTemplate,
    Markup,
    docs_link,
)
from antismash.common.layers import OptionsLayer, RecordLayer, RegionLayer
from antismash.common.module_results import ModuleResults
from antismash.common.secmet import Module, Record, Region
from antismash.detection import nrps_pks_domains
from antismash.detection.nrps_pks_domains import (
    ModularDomain,
    NRPSPKSDomains as DetectionResults,
)
from antismash.detection.nrps_pks_domains.domain_drawing import (
    get_css_class_and_abbreviation,
)
from antismash.detection.nrps_pks_domains.domain_identification import HMMResult
from antismash.detection.nrps_pks_domains.module_identification import (
    CARRIER_PROTEINS,
    MODIFIERS,
    SPECIAL,
)
from antismash.modules import nrps_pks
from antismash.modules.nrps_pks.html_output import NrpspksLayer
from antismash.modules.nrps_pks.results import (
    CandidateClusterPrediction,
    NRPS_PKS_Results as AnalysisResults,
)


assert DetectionResults.schema_version == 4, "nrps_pks_domains results version mismatch, update required"
assert AnalysisResults.schema_version == 3, "nrps_pks results version mismatch, update required"


def has_enough_results(_record: Record, region: Region, results: Dict[str, ModuleResults]) -> bool:
    """ Checks if enough information is present to create at least one
        output visualisation HTML section

        Arguments:
            record: the parent record
            region: the region to check
            results: the results of all detection and analysis modules

        Returns:
            True if an HTML section would be created
    """
    if not nrps_pks.will_handle(region.products, region.product_categories) or nrps_pks.__name__ not in results:
        return False
    analysis_results = results[nrps_pks.__name__]
    assert isinstance(analysis_results, AnalysisResults)
    for candidate in analysis_results.region_predictions.get(region.get_region_number(), []):
        if candidate.polymer:
            return True
    return False


def generate_html(region_layer: RegionLayer, all_results: Dict[str, ModuleResults],
                  record_layer: RecordLayer, options_layer: OptionsLayer) -> HTMLSections:
    """ Generate the sidepanel HTML with results from the NRPS/PKS module """
    results = all_results[nrps_pks.__name__]
    assert isinstance(results, AnalysisResults)

    html = HTMLSections("nrps_pks")

    nrps_layer = NrpspksLayer(results, region_layer.region_feature, record_layer)

    tooltip = Markup(
        f"Shows module structures for each candidate cluster in {docs_link('NRPS', 'glossary/#nrps')} "
        f"and {docs_link('PKS', 'glossary/#t1pks')} regions. "
        "<br>"
        "Genes are shown in predicted order, and are only present when containing at least one complete module. "
        "<br>"
        f"A domain glossary is available {docs_link('here', 'modules/nrps_pks_domains/')}, "
        "and an explanation of the visualisation is available "
        f"{docs_link('here', 'modules/nrps_pks_modules/#bubblemodule-only-drawing')}. "
    )
    template = FileTemplate(path.get_full_path(__file__, "templates", "bubble_view.html"))
    section = template.render(record=record_layer, region=nrps_layer, docs_url=options_layer.urls.docs_baseurl,
                              tooltip_text=tooltip)
    html.add_detail_section("NRPS/PKS modules", section, "nrps-pks-bubbles")

    return html


def generate_javascript_data(record: Record, region: Region, results: Dict[str, ModuleResults]) -> Dict[str, Any]:
    """ Generates JSON data for the javascript to draw relevant results in HTML output

        Arguments:
            record: the relevant Record for the results
            region: the specific Region to generate data for
            results: the results that need data extracted

        Returns:
            a JSON-friendly dictionary with the relevant data
    """
    nrps_pks_results = results[nrps_pks.__name__]
    assert isinstance(nrps_pks_results, AnalysisResults)
    detection_results = results[nrps_pks_domains.__name__]
    assert isinstance(detection_results, DetectionResults)

    data: Dict[str, Any] = {}
    inactive_kr_domains = set()
    for domain_id, predictions in nrps_pks_results.domain_predictions.items():
        activity = predictions.get("kr_activity")
        if not activity:
            continue
        if not isinstance(activity, nrps_pks.data_structures.SimplePrediction):
            raise TypeError(f"{AnalysisResults} format is not as expected")
        if activity.prediction == "inactive":
            inactive_kr_domains.add(domain_id)
    hit_by_domain_by_cds_name = {}
    for cds, cds_result in detection_results.cds_results.items():
        flipped = {v: k for k, v in cds_result.domain_features.items()}
        hit_by_domain_by_cds_name[cds.get_name()] = flipped
    for result in nrps_pks_results.region_predictions[region.get_region_number()]:
        cand_json = _gen_js_data_for_candidate(record, result, inactive_kr_domains, hit_by_domain_by_cds_name)
        data[f"CC{result.candidate_cluster_number}"] = cand_json
    return data


def build_domain_json(profile_name: str, domain: ModularDomain, inactive: bool) -> Dict[str, Any]:
    """ Gather relevant data for drawing the domain

        Arguments:
            profile_name: the name of the profile the domain was hit with
            domain: the domain in question
            inactive: whether the domain is predicted to be inactive

        Returns:
            a JSON-friendly dictionary of the gathered data
    """
    # reuse the NRPS/PKS CSS for domain colouring and labelling (where it fits)
    css_class, name = get_css_class_and_abbreviation(profile_name)

    # though label the normally unlabelled carrier proteins
    if profile_name in CARRIER_PROTEINS:
        name = "CP"
    # and if the abbreviation is too long, don't label it at all
    if len(name) > 3:
        name = ""

    terminal_docking = ""
    if profile_name.endswith("_Nterm"):
        terminal_docking = "start"
    elif profile_name.endswith("_Cterm"):
        terminal_docking = "end"

    # if it's a 'special' domain, it'll draw differently
    special = terminal_docking or profile_name in SPECIAL

    # if it's a modifier, it'll be drawn at different heights
    modification = profile_name in MODIFIERS

    # then flesh out the description for tooltips
    description = profile_name
    types = domain.subtypes
    if types:
        if types[0] == profile_name:
            types = types[1:]
        if types:
            description += f"({')('.join(types)})"
    if inactive:
        description += " - inactive"

    return {
        "name": name,
        "description": description,
        "modifier": modification,
        "special": special,
        "cds": domain.locus_tag,
        "css": css_class,
        "inactive": inactive,
        "start": int(domain.protein_location.start),
        "terminalDocking": terminal_docking,
    }


@dataclass
class SimpleModule:
    """ A basic datastructure to simplify JSON conversion of most values """
    domains: List[Dict[str, Any]]
    complete: bool
    cds_position: int
    polymer: str = ""
    iterative: bool = False

    def to_json(self) -> Dict[str, Any]:
        """ Returns drawing-relevant JSON for the instance"""
        return {
            "domains": self.domains,
            "complete": self.complete,
            "iterative": self.iterative,
            "polymer": self.polymer,
        }


def _gen_js_data_for_candidate(record: Record, result: CandidateClusterPrediction,
                               inactive_domains: Set[str],
                               hit_by_domain_by_cds: Dict[str, Dict[ModularDomain, HMMResult]]
                               ) -> Dict[str, Any]:
    """ Generates module drawing JSON data for a specific candidate

        Arguments:
            record: the relevant Record for the results
            result: the result for a specific candidate cluster
            inactive_domains: a set of domain ids that are inactive
            hit_by_domain_by_cds: the results for NRPS/PKS domain detection,
                                  organised by CDS name and the ModularDomain

        Returns:
            a JSON-friendly dictionary with the relevant data
    """
    order = result.ordering
    if not order:
        return {}
    module_features: List[Module] = []
    for name in order:
        cds_modules = record.get_cds_by_name(name).modules
        if any(mod.is_complete() for mod in cds_modules):
            module_features.extend(cds_modules)

    cds_positions = {name: pos for pos, name in enumerate(order)}

    processed_pairs = set()
    modules: List[SimpleModule] = []
    domains_in_modules_by_cds = defaultdict(set)
    for module in module_features:
        if module.is_multigene_module():
            if module.parent_cds_names in processed_pairs:
                continue
            processed_pairs.add(module.parent_cds_names)

        domains = []
        for dom in module.domains:
            # some domains may be in cross-CDS modules where the parent CDS isn't in the candidate
            if dom.locus_tag not in cds_positions:
                continue
            assert isinstance(dom, ModularDomain)
            inactive = dom.domain_id in inactive_domains
            profile_name = hit_by_domain_by_cds[dom.locus_tag][dom].hit_id
            domains.append(build_domain_json(profile_name, dom, inactive))
            domains_in_modules_by_cds[dom.locus_tag].add(dom)
        polymer = module.get_substrate_monomer_pairs()[0][1] if module.is_complete() else ""

        # if not all of the parent CDS features of this module are in the candidate, it's not complete
        contained = list(filter(lambda parent: parent in cds_positions, module.parent_cds_names))
        complete = module.is_complete() and len(contained) == len(module.parent_cds_names)

        # if the leading CDS isn't in the candidate, use the name of the trailing CDS
        name = domains[0]["cds"]
        if name not in cds_positions:
            name = domains[-1]["cds"]

        modules.append(SimpleModule(domains, complete, cds_positions[name],
                                    polymer, module.is_iterative()))

    extras = []
    for name in order:
        if name not in domains_in_modules_by_cds:
            continue
        existing = domains_in_modules_by_cds[name]
        hits_by_domain = hit_by_domain_by_cds[name]
        for domain in hit_by_domain_by_cds[name]:
            if domain in existing:
                continue
            profile_name = hits_by_domain[domain].hit_id
            extras.append(SimpleModule([build_domain_json(profile_name, domain, False)], False, cds_positions[name]))

    modules = sorted(modules + extras, key=lambda x: (x.cds_position, x.domains[0]["start"]))

    return {
        "modules": [module.to_json() for module in modules],
    }
