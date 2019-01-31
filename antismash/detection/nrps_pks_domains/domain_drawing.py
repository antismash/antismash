# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Generates HTML and JSON for the nrps_pks_domains module

    Generated content will include a prediction consensus in the drawn domains,
    but will not cover details or options. That is left to the predicting module
    to describe in the sidepanel.
"""

from typing import Dict, List, Union

from antismash.common import path
from antismash.common.html_renderer import FileTemplate, HTMLSections
from antismash.common.json import JSONDomain, JSONOrf
from antismash.common.layers import RegionLayer, RecordLayer, OptionsLayer
from antismash.common.module_results import ModuleResults
from antismash.common.secmet import CDSFeature, Record, Region
from antismash.common.secmet.qualifiers import NRPSPKSQualifier


def will_handle(products: List[str]) -> bool:
    """ Returns true if one or more relevant products are present """
    return bool(set(products).intersection({"nrps", "t1pks", "t2pks", "transatpks",
                                            "nrpsfragment", "otherks"}))


def _parse_domain(record: Record, domain: NRPSPKSQualifier.Domain,
                  feature: CDSFeature) -> JSONDomain:
    """ Convert a NRPS/PKS domain string to a dict useable by json.dumps

        Arguments:
            record: the Record containing the domain
            domain: the NRPSPKSQualifier.Domain in question
            feature: the CDSFeature that the domain belongs to

        Returns:
            a populated JSONDomain instance
    """
    predictions = list(domain.predictions.items())

    # Create url_link to NaPDoS for C and KS domains
    napdoslink = ""
    domainseq = str(feature.translation)[domain.start:domain.end]
    base = ("http://napdos.ucsd.edu/cgi-bin/process_request.cgi?"
            "query_type=aa&amp;ref_seq_file=all_{0}_public_12062011.faa"
            "&amp;Sequence=%3E{0}_domain_from_antiSMASH%0D{1}")
    if domain.name == "PKS_KS":
        napdoslink = base.format("KS", domainseq)
    elif "Condensation" in domain.name:
        napdoslink = base.format("C", domainseq)
    blastlink = ("http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE=Proteins"
                 "&amp;PROGRAM=blastp&amp;BLAST_PROGRAMS=blastp"
                 "&amp;QUERY={}"
                 "&amp;LINK_LOC=protein&amp;PAGE_TYPE=BlastSearch").format(domainseq)

    dna_sequence = feature.extract(record.seq)
    return JSONDomain(domain, predictions, napdoslink, blastlink, domainseq, dna_sequence)


def generate_js_domains(region: Region, record: Record) -> Dict[str, Union[str, List[JSONOrf]]]:
    """ Creates a JSON-like structure for domains, used by javascript in
        drawing the domains
    """
    orfs = []  # type: List[JSONOrf]
    for feature in region.cds_children:
        if not feature.nrps_pks:
            continue
        js_orf = JSONOrf(feature)
        for domain in feature.nrps_pks.domains:
            js_orf.add_domain(_parse_domain(record, domain, feature))
        orfs.append(js_orf)

    return {'id': RegionLayer.build_anchor_id(region),
            'orfs': orfs}


def has_domain_details(region: Region) -> bool:
    """ Returns True if there are domain details to be had for the given cluster """
    for cds in region.cds_children:
        if cds.nrps_pks:
            return True
    return False


def generate_html(region_layer: RegionLayer, _results: ModuleResults,
                  _record_layer: RecordLayer, _options_layer: OptionsLayer
                  ) -> HTMLSections:
    """ Generate the details section of NRPS/PKS domains in the main HTML output """
    template = FileTemplate(path.get_full_path(__file__, 'templates', 'details.html'))
    section = template.render(has_domain_details=has_domain_details, region=region_layer)
    html = HTMLSections("nrps_pks")
    html.add_detail_section("NRPS/PKS domains", section)
    return html
