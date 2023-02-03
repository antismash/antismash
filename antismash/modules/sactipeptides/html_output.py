# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Handles HTML output for sactipeptides """

from collections import defaultdict
from typing import Dict, List, Set

from antismash.common import comparippson, path
from antismash.common.html_renderer import HTMLSections, FileTemplate, Markup
from antismash.common.secmet import Prepeptide
from antismash.common.layers import RegionLayer, RecordLayer, OptionsLayer

from .specific_analysis import SactiResults


def will_handle(products: List[str], _product_categories: Set[str]) -> bool:
    """ Returns true if one or more relevant products or product categories are present """
    return 'sactipeptide' in products


def generate_html(region_layer: RegionLayer, results: SactiResults,
                  record_layer: RecordLayer, options_layer: OptionsLayer) -> HTMLSections:
    """ Generates HTML for the module """
    html = HTMLSections("sactipeptides")

    motifs_by_locus: Dict[str, List[Prepeptide]] = defaultdict(list)
    for locus, motifs in results.motifs_by_locus.items():
        for motif in motifs:
            if motif.is_contained_by(region_layer.region_feature):
                motifs_by_locus[locus].append(motif)

    motifs_by_locus_by_core = {}
    for locus, motifs in motifs_by_locus.items():
        motifs_by_core = defaultdict(list)
        for motif in motifs:
            core = motif.core
            motifs_by_core[core].append(motif)
        motifs_by_locus_by_core[locus] = motifs_by_core

    detail_tooltip = Markup("<br>".join([
        (
            "Lists the possible core peptides for each biosynthetic enzyme. "
            "Each core peptide shows the leader and core peptide sequences. "
        ),
        comparippson.get_tooltip_text(),
    ]))
    template = FileTemplate(path.get_full_path(__file__, "templates", "details.html"))
    details = template.render(record=record_layer,
                              options=options_layer,
                              motifs_by_locus=motifs_by_locus_by_core,
                              comparippson_results=results.comparippson_results,
                              tooltip=detail_tooltip)
    html.add_detail_section("Sactipeptides", details)

    side_tooltip = ("Lists the possible core peptides in the region. "
                    "Each core peptide lists its RODEO score and predicted core sequence.")
    template = FileTemplate(path.get_full_path(__file__, "templates", "sidepanel.html"))
    sidepanel = template.render(record=record_layer,
                                options=options_layer,
                                motifs_by_locus=motifs_by_locus,
                                tooltip=side_tooltip)
    html.add_sidepanel_section("Sactipeptides", sidepanel)

    return html
