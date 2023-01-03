# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Generates HTML and JSON for the nrps_pks_domains module

    Generated content will include a prediction consensus in the drawn domains,
    but will not cover details or options. That is left to the predicting module
    to describe in the sidepanel.
"""

import itertools
import string
from typing import Any, Iterator, Optional, Union

from antismash.common import path
from antismash.common.html_renderer import (
    FileTemplate,
    HTMLSections,
    Markup,
    docs_link,
    replace_with,
)
from antismash.common.json import JSONDomain, JSONOrf, JSONModule
from antismash.common.layers import RegionLayer, RecordLayer
from antismash.common.module_results import ModuleResults
from antismash.common.secmet import CDSFeature, Module, Record, Region
from antismash.common.secmet.qualifiers import NRPSPKSQualifier
from antismash.config import ConfigType
from antismash.detection import nrps_pks_domains
from antismash.modules import nrps_pks

_UNLABLED_DOMAINS = set([
    "PCP",
    "ACP",
    "ACP_beta",
    "ACPS",
    "LPG_synthase_C",
    "NRPS-COM_Nterm",
    "NRPS-COM_Cterm",
    "PKS_Docking_Nterm",
    "PKS_Docking_Cterm",
    "PKS_PP",
    "PP-binding",
    "Trans-AT_docking",
    "Polyketide_cyc",
    "Polyketide_cyc2",
    "TIGR01720",
    "TIGR02353",
])
_ABBREVIATED_DOMAINS = {
    "Aminotran_1_2": "AmT",
    "Aminotran_3": "AmT",
    "Aminotran_4": "AmT",
    "Aminotran_5": "AmT",
    "AMP-binding": "A",
    "A-OX": "A-OX",
    "Beta_elim_lyase": "SH",
    "Cglyc": "C",
    "Condensation_DCL": "C",
    "Condensation_LCL": "C",
    "Condensation_Starter": "C",
    "Condensation_Dual": "C",
    "Heterocyclization": "C",
    "Epimerization": "E",
    "Thioesterase": "TE",
    "PKS_KS": "KS",
    "PKS_AT": "AT",
    "PKS_KR": "KR",
    "PKS_DH": "DH",
    "PKS_DH2": "DH",
    "PKS_DHt": "DHt",
    "PKS_ER": "ER",
    "Interface": "I",
    "IBH_Asp": "βH",
}
_CLASS_BY_ABBREVIATION = {
    "A": "adenylation",
    "A-OX": "adenylation",
    "AT": "acyltransferase",
    "C": "condensation",
    "E": "epimerase",
    "TD": "terminal",
    "TE": "terminal",
    "KS": "ketosynthase",
    "KR": "mod-kr",
    "DH": "mod-dh",
    "DHt": "mod-dh",
    "DH2": "mod-dh",
    "ER": "mod-er",
    "SH": "mod-sh",
}
_CLASS_BY_NAME = {
    "PCP": "transport",
    "ACP": "transport",
    "ACP_beta": "transport",
    "PKS_PP": "transport",
    "PP-binding": "transport",
    "NRPS-COM_Nterm": "docking",
    "NRPS-COM_Cterm": "docking",
    "PKS_Docking_Nterm": "docking",
    "PKS_Docking_Cterm": "docking",
    "Trans-AT_docking": "docking",
    "TD": "terminal",
}
_NAPDOS_REFERENCES = {
    "C": "all_C_190411_273.faa",
    "KS": "all_KS_191020_1877.faa",
}


def will_handle(_products: list[str], product_categories: set[str]) -> bool:
    """ Returns true if one or more relevant products or product categories are present """
    return bool(product_categories.intersection({"NRPS", "PKS"}))


def has_enough_results(_record: Record, region: Region, results: dict[str, ModuleResults]) -> bool:
    """ Checks if enough information is present to create at least one
        output visualisation HTML section

        Arguments:
            record: the parent record
            region: the region to check
            results: the results of all detection and analysis modules

        Returns:
            True if an HTML section would be created
    """
    # if no interesting modules have results for the whole record, get out early
    if not any(results.get(module.__name__) for module in [nrps_pks_domains, nrps_pks]):
        return False
    return bool(region.product_categories.intersection({"NRPS", "PKS"}))


def _get_domain_abbreviation(domain_name: str) -> str:
    """ Convert full domain name to abbreviation (if any) for HTML display) """
    if domain_name in _UNLABLED_DOMAINS:
        return ""
    return _ABBREVIATED_DOMAINS.get(domain_name, domain_name.split("_", 1)[0])


def _get_domain_class(abbreviation: str, domain_name: str) -> str:
    """ Convert full abbreviation (if any) or domain name to an HTML class for styling) """
    if abbreviation:
        res = _CLASS_BY_ABBREVIATION.get(abbreviation, "other")
    else:
        res = _CLASS_BY_NAME.get(domain_name, "other")
    return f"jsdomain-{res}"


def get_css_class_and_abbreviation(domain_name: str) -> tuple[str, str]:
    """ Convert a full domain name to a pair of CSS class and abbrevation """
    abbrevation = _get_domain_abbreviation(domain_name)
    css_class = _get_domain_class(abbrevation, domain_name)
    return css_class, abbrevation


def _build_paras_link(cds: CDSFeature, domain: NRPSPKSQualifier.Domain, signature: str) -> str:
    link = (
        "http://paras.bioinformatics.nl/api/submit_quick?"
        f"signature1={signature}"
        f"&amp;name1={cds.get_name()}"
        f"&amp;start1={domain.start}"
        f"&amp;end1={domain.end}"
    )
    return f'<a class="external-link" href="{ link }" target="_blank">Repredict substrate with PARAS</a>'


def _parse_domain(record: Record, domain: NRPSPKSQualifier.Domain,
                  feature: CDSFeature, signature: str = "",
                  ) -> JSONDomain:
    """ Convert a NRPS/PKS domain string to a dict useable by json.dumps

        Arguments:
            record: the Record containing the domain
            domain: the NRPSPKSQualifier.Domain in question
            feature: the CDSFeature that the domain belongs to

        Returns:
            a populated JSONDomain instance
    """
    predictions = list(domain.get_predictions().items())

    # Create url_link to NaPDoS for C and KS domains
    napdoslink = ""
    domainseq = str(feature.translation)[domain.start:domain.end]
    base = ("https://npdomainseeker.sdsc.edu/cgi-bin/process_request_napdos2.cgi?"
            "query_type=aa&amp;ref_seq_file={0}"
            "&amp;Sequence=%3E{0}_domain_from_antiSMASH%0D{1}")
    if domain.name == "PKS_KS":
        napdoslink = base.format(_NAPDOS_REFERENCES["KS"], replace_with("sequence"))
    elif "Condensation" in domain.name:
        napdoslink = base.format(_NAPDOS_REFERENCES["C"], replace_with("sequence"))
    blastlink = ("http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE=Proteins"
                 "&amp;PROGRAM=blastp&amp;BLAST_PROGRAMS=blastp"
                 f"&amp;QUERY={replace_with('sequence')}"
                 "&amp;LINK_LOC=protein&amp;PAGE_TYPE=BlastSearch")
    extra_links = []
    if signature and ("AMP-binding" in domain.name or "A-OX" in domain.name):
        extra_links.append(_build_paras_link(feature, domain, signature))

    dna_sequence = ""
    for as_domain in record.get_antismash_domains_in_cds(feature):
        if as_domain.get_name() == domain.feature_name:
            dna_sequence = as_domain.extract(record.seq)
            break
    assert dna_sequence
    css, abbreviation = get_css_class_and_abbreviation(domain.name)
    return JSONDomain.from_domain(
        domain, predictions, napdoslink, blastlink, domainseq, dna_sequence,
        abbreviation, css,
        extra_links=extra_links,
    )


def _build_module_js(module: Module, cds: CDSFeature, match_ids: dict[tuple[str, ...], str],
                     match_gen: Iterator[str]) -> JSONModule:
    """ Builds and returns a JSONModule isntance to match the given module """
    monomer = ""
    if module.monomers:
        monomer = module.monomers[0][1]
        if monomer.endswith("pk"):
            monomer = monomer[:-2] + "?"
        if monomer.endswith("X"):
            monomer = monomer[:-1] + "?"
    multi_cds: Optional[str] = None
    match_id: Optional[str] = None
    protein_location = module.get_parent_protein_location(cds.get_name())
    # determine which CDS is the head if the module crosses two CDS features
    if len(module.parent_cds_names) > 1:
        if module.parent_cds_names[0] == cds.get_name():
            multi_cds = "head"
        else:
            multi_cds = "tail"
        if module.parent_cds_names not in match_ids:
            match_ids[module.parent_cds_names] = next(match_gen)
        match_id = match_ids[module.parent_cds_names]
    js_module = JSONModule(protein_location.start, protein_location.end,
                           module.is_complete(), module.is_iterative(), monomer,
                           multi_cds=multi_cds, match_id=match_id)
    return js_module


def has_domain_details(region: Union[Region, RegionLayer]) -> bool:
    """ Returns True if there are domain details to be had for the given cluster """
    for cds in region.cds_children:
        if cds.nrps_pks:
            return True
    return False


def domains_have_predictions(region: Union[Region, RegionLayer]) -> bool:
    """ Returns True if any domain in the region has a prediction made for it """
    for feature in region.cds_children:
        for domain in feature.nrps_pks.domains:
            if "substrate consensus" in domain.get_predictions():
                return True
    return False


def generate_javascript_data(record: Record, region: Region,
                             results: dict[str, ModuleResults]) -> dict[str, Any]:
    """ Generate the javascript data required for interactive visualisation """
    orfs: list[JSONOrf] = []
    match_ids: dict[tuple[str, ...], str] = {}

    def match_id_generator() -> Iterator[str]:
        """ Generates match names as A, B, .. Z, AA, AB, .. AZ, AAA, ... """
        for size in itertools.count(start=1):
            for i in itertools.product(string.ascii_uppercase, repeat=size):
                yield "".join(i)

    match_gen = iter(match_id_generator())

    analysis_results = results.get(nrps_pks.__name__)
    if analysis_results:
        assert isinstance(analysis_results, nrps_pks.NRPS_PKS_Results)
        full_predictions = analysis_results.domain_predictions
    else:
        full_predictions = {}

    for feature in region.cds_children:
        if not feature.nrps_pks:
            continue
        js_orf = JSONOrf.from_cds(feature)
        for domain in feature.nrps_pks.domains:
            signatures = full_predictions.get(domain.feature_name, {}).get("nrpys")
            signature = ""
            if signatures:
                assert isinstance(signatures, nrps_pks.results.PredictorSVMResult)
                signature = signatures.aa34
            parsed = _parse_domain(record, domain, feature, signature)
            js_orf.add_domain(parsed)

        for module in feature.modules:
            js_module = _build_module_js(module, feature, match_ids, match_gen)
            js_orf.add_module(js_module)
        orfs.append(js_orf)

    return {
        "id": RegionLayer.build_anchor_id(region),
        "orfs": orfs,
    }


def generate_html(region_layer: RegionLayer, _results: dict[str, ModuleResults],
                  record_layer: RecordLayer, _options: ConfigType,
                  ) -> HTMLSections:
    """ Generate the details section of NRPS/PKS domains in the main HTML output """
    template = FileTemplate(path.get_full_path(__file__, "templates", "nrps_pks_domains.html"))
    html = HTMLSections("nrps_pks")
    if not has_domain_details(region_layer):
        return html

    # hide lids by default if none have predictions (e.g. in a minimal run)
    hide_lids = not domains_have_predictions(region_layer)

    tooltip = Markup(
        f"Shows {docs_link('NRPS', 'glossary/#nrps')}- and {docs_link('PKS', 'glossary/#t1pks')}-"
        "related domains for each feature that contains them. "
        "Click on each domain for more information about the domain's location, "
        "consensus monomer prediction, and other details.<br>"
        f"A domain glossary is available {docs_link('here', 'modules/nrps_pks_domains/')}, "
        f"and an explanation of the visualisation is available {docs_link('here', 'modules/nrps_pks_modules/')}."
    )

    section = template.render(has_domain_details=has_domain_details, region=region_layer,
                              record=record_layer,
                              tooltip_text=tooltip, hide_lids=hide_lids)
    html.add_detail_section("NRPS/PKS domains", section)
    return html
