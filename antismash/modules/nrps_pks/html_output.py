# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Generates HTML and JSON for the nrps_pks module """

import logging
import re
from typing import Any, Dict, Iterable, List, Optional, Set

from antismash.common import path
from antismash.common.html_renderer import HTMLSections, FileTemplate
from antismash.common.layers import RegionLayer, RecordLayer, OptionsLayer
from antismash.common.secmet import CDSFeature, Region, CandidateCluster

from .results import NRPS_PKS_Results, CandidateClusterPrediction, UNKNOWN


def will_handle(_products: List[str], categories: Set[str]) -> bool:
    """ Returns true if one or more relevant products or product categories are present """
    return bool(categories.intersection({"NRPS", "PKS"}))


def generate_html(region_layer: RegionLayer, results: NRPS_PKS_Results,
                  record_layer: RecordLayer, options_layer: OptionsLayer) -> HTMLSections:
    """ Generate the sidepanel HTML with results from the NRPS/PKS module """
    html = HTMLSections("nrps_pks")

    nrps_layer = NrpspksLayer(results, region_layer.region_feature, record_layer)

    features_with_domain_predictions: Dict[str, List[str]] = {}
    for domain_name, consensus in results.consensus.items():
        if not consensus:
            continue
        domain = record_layer.get_domain_by_name(domain_name)
        features_with_domain_predictions[domain.locus_tag] = []

    for feature_name, monomers in features_with_domain_predictions.items():
        for domain in record_layer.get_cds_by_name(feature_name).nrps_pks.domains:
            monomer = results.consensus.get(domain.feature_name)
            if monomer:
                monomers.append(monomer)

    prod_tt = ("Shows estimated product structure and polymer for each candidate cluster in the region. "
               "To show the product, click on the expander or the candidate cluster feature drawn in the overview. "
               )
    mon_tt = ("Shows the predicted substrates for each adenylation domain and acyltransferase within genes. "
              "Each gene prediction can be expanded to view detailed predictions of each domain. "
              "Each prediction can be expanded to view the predictions by tool "
              " (and, for some tools, further expanded for extra details). "
              )

    template_components = []
    # only show a polymer tab if there's polymers
    if nrps_layer.has_any_polymer():
        template_components.append(("products.html", "NRPS/PKS products", "nrps_pks_products", prod_tt))

    # always include monomers, if available
    if features_with_domain_predictions:
        template_components.append(("monomers.html", "NRPS/PKS substrates", "nrps_pks_monomers", mon_tt))

    for filename, name, class_name, tooltip in template_components:
        template = FileTemplate(path.get_full_path(__file__, "templates", filename))
        section = template.render(record=record_layer,
                                  region=nrps_layer,
                                  results=results,
                                  relevant_features=features_with_domain_predictions,
                                  options=options_layer,
                                  tooltip=tooltip)
        html.add_sidepanel_section(name, section, class_name)

    return html


def filter_norine_as(monomers: List[str], be_strict: bool = False) -> List[str]:
    """ Remove PKS and unknown substrate predictions
        use be_strict = True to also filter nrp/X
    """
    filtered_list = []
    bad_monomers = {'pk', UNKNOWN, 'hydrophilic', 'hydrophobic', 'mal', 'mmal'}
    if be_strict:
        bad_monomers = bad_monomers.union({'nrp', 'x'})
    for monomer in monomers:
        monomer = monomer.lower()
        if monomer in bad_monomers:
            continue
        assert '|' not in monomer, monomer
        if monomer != 'x' and not be_strict:
            monomer += '*'
        filtered_list.append(monomer)
    return filtered_list


def get_norine_url_for_specificities(specificities: List[List[str]],
                                     be_strict: bool = True) -> Optional[str]:
    """ generate NORINE URL string for direct querying from array of specificity predictions
        use be_strict=False to add * after each monomer
    """
    modules = []
    if not specificities:
        return None

    for domain_specificity_list in specificities:
        assert isinstance(domain_specificity_list, list), type(domain_specificity_list)
        if not domain_specificity_list:
            logging.error("retrieved empty domain list for assembling per protein NORINE link string")
            raise RuntimeError("empty domain list for NORINE generation")
        # always remove X for this particular URL/query
        chunks = [chunk for chunk in domain_specificity_list if chunk != 'X']
        chunks = filter_norine_as(chunks, be_strict=be_strict)
        if not chunks:
            continue
        query = "|".join(chunks)
        if len(chunks) > 1:
            query = "[" + query + "]"
        modules.append(query)

    if not modules:
        return None

    return "http://bioinfo.lifl.fr/norine/fingerPrintSearch.jsp?nrps1=" + ",".join(modules)


class CandidateClusterLayer:
    """ A helper for the HTML output for a candidate_cluster """
    def __init__(self, candidate_cluster: CandidateCluster, result: CandidateClusterPrediction) -> None:
        self.location = candidate_cluster.location
        self.number = candidate_cluster.get_candidate_cluster_number()
        self.transatpks = "transatpks" in candidate_cluster.products
        self.result = result
        self.products = "-".join(candidate_cluster.products)
        self.kind = str(candidate_cluster.kind).replace("_", " ")

    def __getattr__(self, attr: str) -> Any:
        if hasattr(self.result, attr):
            return getattr(self.result, attr)
        return super().__getattribute__(attr)

    def get_warning(self) -> str:
        """ A caveat for structure prediction accuracy """
        if not self.polymer:
            return ""
        core = ("Rough prediction of core scaffold based on assumed %s;"
                " tailoring reactions not taken into account")
        if self.domain_docking_used:
            detail = "PKS linker matching"
        else:
            detail = "PKS/NRPS colinearity"

        return core % detail

    def get_norine_url(self, be_strict: bool = True) -> str:
        """ Get a NORINE URL string for direct querying
            use be_strict=False to add * after each monomer"""

        monomers_per_protein_list = re.findall("\\(.*?\\)", self.polymer)
        i = 1
        nrpslist = []
        for monomers_per_protein in monomers_per_protein_list:
            monomers = monomers_per_protein[1:-1].split(" - ")

            monomers = filter_norine_as(monomers, be_strict=be_strict)
            if monomers:
                nrpslist.append("nrps" + str(i) + "=" + ",".join(monomers))
            i += 1
        urlstring = "http://bioinfo.lifl.fr/norine/fingerPrintSearch.jsp?"+"&".join(nrpslist)
        return urlstring


class NrpspksLayer(RegionLayer):
    """ A wrapper for RegionLayer that adds some specific sections for NRPS/PKS
        domains and structures.
    """
    def __init__(self, results: NRPS_PKS_Results, region_feature: Region, record: RecordLayer) -> None:
        self.url_strict: Dict[str, str] = {}  # gene name -> url
        self.url_relaxed: Dict[str, str] = {}  # gene name -> url
        self._build_urls(region_feature.cds_children)
        super().__init__(record, region_feature)
        assert isinstance(results, NRPS_PKS_Results), type(results)
        self.results = results

        region_number = region_feature.get_region_number()
        self.candidate_clusters: List[CandidateClusterLayer] = []
        for candidate_cluster_pred in results.region_predictions.get(region_number, []):
            candidate_cluster = record.get_candidate_cluster(candidate_cluster_pred.candidate_cluster_number)
            self.candidate_clusters.append(CandidateClusterLayer(candidate_cluster, candidate_cluster_pred))

    @property
    def sidepanel_features(self) -> List[str]:
        """ Returns a list of relevant CDSFeature names """
        return [feature.get_name() for feature in self.region_feature.cds_children if feature.nrps_pks]

    def _build_urls(self, cds_features: Iterable[CDSFeature]) -> None:
        for feature in cds_features:
            if not feature.nrps_pks:
                continue

            feature_name = feature.get_name()

            per_cds_predictions = []

            for domain in feature.nrps_pks.domains:
                if domain.name not in ["AMP-binding", "A-OX"]:
                    continue
                per_a_domain_predictions = set()
                for possibilities in domain.get_predictions().values():
                    for possibility in possibilities.split(","):
                        per_a_domain_predictions.add(possibility)
                per_cds_predictions.append(list(per_a_domain_predictions))

            if not per_cds_predictions:
                continue
            url = get_norine_url_for_specificities(per_cds_predictions)
            if url:
                self.url_strict[feature_name] = url
                url = get_norine_url_for_specificities(per_cds_predictions, be_strict=False)
                if url:
                    self.url_relaxed[feature_name] = url

    def has_any_polymer(self) -> bool:
        """ Does the region contain at least one candidate_cluster with a polymer set """
        return any(candidate.polymer for candidate in self.candidate_clusters)
