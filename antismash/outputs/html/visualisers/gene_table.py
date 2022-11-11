# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" HTML and JSON generation for a gene/CDS table summary """

from typing import Any, Dict

from antismash.common import path
from antismash.common.html_renderer import (
    build_blastp_link,
    FileTemplate,
    HTMLSections,
    selected_cds_marker,
    WILDCARD_TEMPLATE,
)
from antismash.common.module_results import ModuleResults
from antismash.common.layers import OptionsLayer, RecordLayer, RegionLayer
from antismash.common.secmet import Record, Region


TOOLTIP = (
    "A brief tabular summary of genes/CDS features within the region.<br>"
    "Filtering the table will also search biosynthetic profiles and gene function data. "
    "If enabled, the overview will then zoom to show the area covered by the "
    "filtered selection.<br>"
    "Genes selected in the region drawing above will be marked in the table with "
    "an indicator to the left of the gene name."
)


def has_enough_results(_record: Record, region: Region, _results: Dict[str, ModuleResults]) -> bool:
    """ Checks if enough information is present to create at least one
        output visualisation HTML section

        Arguments:
            record: the parent record
            region: the region to check
            results: the results of all detection and analysis modules

        Returns:
            True if an HTML section would be created
    """
    return bool(region.cds_children)


def generate_html(region_layer: RegionLayer, all_results: Dict[str, ModuleResults],  # pylint: disable=unused-argument
                  record_layer: RecordLayer, _options_layer: OptionsLayer) -> HTMLSections:  # pylint: disable=unused-argument
    """ Generates the HTML sections for all relevant module results
    """
    html = HTMLSections("gene-table")

    template = FileTemplate(path.get_full_path(__file__, "templates", "gene_table.html"))

    section = template.render(region=region_layer, tooltip=TOOLTIP)
    html.add_detail_section("Gene overview", section, "gene-table-details")

    return html


def generate_javascript_data(record: Record, region: Region,  # pylint: disable=unused-argument
                             results: Dict[str, ModuleResults]) -> Dict[str, Any]:  # pylint: disable=unused-argument
    """ Generates the javascript data for all relevant module results
    """
    data: Dict[str, Any] = {
        "blast_template": build_blastp_link(WILDCARD_TEMPLATE.format("translation"), "BlastP"),
        "selected_template": selected_cds_marker(WILDCARD_TEMPLATE.format("locus_tag")),
        "orfs": {},
    }
    for cds in region.cds_children:
        detail = {
            "functions": [{
                "description": function.description.split("(Score", 1)[0],
                "function": str(function.function),
                "tool": function.tool if function.tool != "rule-based-clusters" else "biosynthetic profile",
            } for function in cds.gene_functions]
        }
        data["orfs"][cds.get_name()] = detail
    return data
