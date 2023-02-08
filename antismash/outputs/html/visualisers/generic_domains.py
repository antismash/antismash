# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" HTML and JSON generation for modules using HMMer to generate domains """

from collections import defaultdict
from dataclasses import dataclass
from typing import Any, Dict, List, Optional

from antismash.common import path
from antismash.common.module_results import ModuleResults

from antismash.common.hmmer import HmmerResults
from antismash.common.html_renderer import HTMLSections, FileTemplate
from antismash.common.layers import OptionsLayer, RecordLayer, RegionLayer
from antismash.common.secmet import PFAMDomain, Record, Region
from antismash.detection import cluster_hmmer, full_hmmer, tigrfam
from antismash.outputs.html import js

assert HmmerResults.schema_version == 2, "HmmerResults version mismatch, update required"

# a simplified and modified version of https://github.com/Alexamk/decRiPPter/tree/master/data/domains
with open(path.get_full_path(__file__, "data", "pfam_transport.txt"), encoding="utf-8") as _handle:
    PFAM_TRANSPORT = set(_handle.read().splitlines())
with open(path.get_full_path(__file__, "data", "pfam_biosynthetic.txt"), encoding="utf-8") as _handle:
    PFAM_BIOSYNTHETIC = set(_handle.read().splitlines())
with open(path.get_full_path(__file__, "data", "pfam_regulatory.txt"), encoding="utf-8") as _handle:
    PFAM_REGULATORY = set(_handle.read().splitlines())


@dataclass
class ToolInfo:
    """ A simple container for various information about each tool, used in
    rendering HTML.
    """
    title: str
    _tooltip: str = ""
    intro: str = ""
    url: str = ""

    @property
    def tooltip(self) -> str:
        """ Returns the tooltip for the tool, doing any string formatting required """
        if "{title}" in self._tooltip:
            return self._tooltip.format(title=self.title)
        return self._tooltip

    def __str__(self) -> str:
        return self.title


HELP = ("Shows {title} domains found in each gene within the region. "
        "Click on each domain for more information about the domain's "
        "accession, location, and description. "
        )

PFAM_TOOL = ToolInfo("Pfam", HELP + "Domains with a bold border have Gene Ontology information.",
                     url="https://www.ebi.ac.uk/interpro/entry/pfam/$ACCESSION")
TOOLS = {
    tigrfam: ToolInfo("TIGRFAM", HELP,
                      url="https://www.ncbi.nlm.nih.gov/genome/annotation_prok/evidence/$ACCESSION"),
}


def region_has_pfams(record: Record, region: Region) -> bool:
    """ Returns true if the region contains any CDS features with a Pfam annotation """
    return any(record.get_pfam_domains_in_cds(cds) for cds in region.cds_children)


def region_has_tool_results(record: Record, region: Region, results: Optional[HmmerResults]) -> bool:
    """ Returns true if the region contains any CDS features with an annotation
        from the given tool
    """
    if not results:
        return False
    assert results.tool
    for cds in region.cds_children:
        domains = record.get_antismash_domains_in_cds(cds)
        if any(results.tool.lower() == (domain.tool or "").lower() for domain in domains):
            return True
    return False


def generate_tool_data(record: Record, region: Region, results: Optional[HmmerResults]) -> List[Dict[str, Any]]:
    """ Generates the javascript data for a specific tool's results within the
        given region

        Arguments:
            record: the record the region belongs to
            region: the region to build data for
            results: the results of the tool

        Returns:
            a JSON-friendly data blob of a list of objects, each having the keys
                'id', 'seqLength', and 'domains', the last of which is a list
                of domain object JSON
    """

    if not results:
        return []

    cds_names = {cds.get_name() for cds in region.cds_children}

    domains_by_cds = defaultdict(list)
    for hit in results.hits:
        if hit.locus_tag in cds_names:
            domains_by_cds[hit.locus_tag].append({
                "start": hit.protein_start,
                "end": hit.protein_end,
                "name": hit.identifier,
                "description": hit.description,
                "accession": hit.identifier,
                "evalue": f"{hit.evalue:g}",
                "score": f"{hit.score:.1f}",
                "html_class": "generic-type-other",
            })

    cds_info = []
    for cds_name, domains in domains_by_cds.items():
        cds_info.append({
            "id": cds_name,
            "seqLength": len(record.get_cds_by_name(cds_name).translation),
            "domains": domains,
        })
    return cds_info


def has_enough_results(record: Record, region: Region, results: Dict[str, ModuleResults]) -> bool:
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
    if not any(results.get(module.__name__) for module in [cluster_hmmer, full_hmmer] + list(TOOLS)):
        return False

    if region_has_pfams(record, region):
        return True

    # otherwise, check if other module results are present in the region
    for module in TOOLS:
        tool_result = results.get(module.__name__)
        if not tool_result:
            continue
        assert isinstance(tool_result, HmmerResults)
        if region_has_tool_results(record, region, tool_result):
            return True
    return False


def generate_html(region_layer: RegionLayer, all_results: Dict[str, ModuleResults],
                  record_layer: RecordLayer, _options_layer: OptionsLayer) -> HTMLSections:
    """ Generate the HTML sections for all relevant module results
    """

    html = HTMLSections("hmmer-domains")

    template = FileTemplate(path.get_full_path(__file__, "templates", "generic_domains.html"))

    record = record_layer.seq_record
    region = region_layer.region_feature

    if region_has_pfams(record, region):
        section = template.render(region=region_layer, record=record_layer, info=PFAM_TOOL)
        html.add_detail_section(f"{PFAM_TOOL.title} domains", section, f"{PFAM_TOOL.title.lower()}-details")

    for tool, info in TOOLS.items():
        results = all_results.get(tool.__name__)
        if not results:
            continue
        assert isinstance(results, HmmerResults)
        if region_has_tool_results(record, region, results):
            section = template.render(region=region_layer, record=record_layer, info=info)
            html.add_detail_section(f"{info.title} domains", section, f"{info.title.lower()}-details")
    assert html.detail_sections, f"missing detail sections for Pfams and {list(TOOLS.keys())}"
    return html


def generate_javascript_data(record: Record, region: Region, results: Dict[str, ModuleResults]) -> List[Dict[str, Any]]:
    """ Converts plain HMMer domain annotations from detection and analysis modules
        into a JSON format for javascript to use

        Arguments:
            record: the record in question
            region: the region for which to generate domain data
            results: the results of all modules

        Returns:
            a JSON-friendly dictionary
    """
    tools = []

    if region_has_pfams(record, region):
        pfam_results = results.get(cluster_hmmer.__name__, results.get(full_hmmer.__name__))
        assert isinstance(pfam_results, HmmerResults)
        tools.append({
            "name": PFAM_TOOL.title.lower(),
            "data": generate_pfam_data(record, region),
            "url": PFAM_TOOL.url,
        })

    for tool, info in TOOLS.items():
        tool_results = results.get(tool.__name__)
        if not tool_results:
            continue
        assert isinstance(tool_results, HmmerResults)
        tools.append({
            "name": info.title.lower(),
            "data": generate_tool_data(record, region, tool_results),
            "url": info.url,
        })

    return tools


def generate_pfam_data(record: Record, region: Region) -> List[Dict[str, Any]]:
    """ Generates Pfam-specific data, including gene ontology information if present """
    def convert_pfam(pfam: PFAMDomain) -> Dict[str, Any]:
        """ Converts a single Pfam to it's JSON representation """
        if pfam.identifier in PFAM_TRANSPORT:
            colour = "transport"
        elif pfam.identifier in PFAM_REGULATORY:
            colour = "regulatory"
        elif pfam.identifier in PFAM_BIOSYNTHETIC:
            colour = "biosynthetic"
        else:
            colour = "other"
        base = {
            "start": pfam.protein_location.start,
            "end": pfam.protein_location.end,
            "go_terms": js.build_pfam2go_links(pfam.gene_ontologies),
            "html_class": f"generic-type-{colour}",
            "name": pfam.domain,
            "accession": pfam.identifier,  # the EBI URLs fail if including the version
            "description": pfam.description,
            "evalue": f"{pfam.evalue:g}",
            "score": f"{pfam.score:.1f}",
        }
        return base

    cdses = []
    for cds in region.cds_children:
        pfams = []
        for pfam in record.get_pfam_domains_in_cds(cds):
            pfams.append(convert_pfam(pfam))
        if pfams:
            cdses.append({
                "id": cds.get_name(),
                "seqLength": len(cds.translation),
                "domains": pfams,
            })
    return cdses
