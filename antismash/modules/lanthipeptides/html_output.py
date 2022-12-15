# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Manages HTML construction for the Lanthipeptide module
"""

from collections import defaultdict
from typing import List, Set

from antismash.common import comparippson, path
from antismash.common.html_renderer import HTMLSections, FileTemplate, Markup
from antismash.common.layers import RegionLayer, RecordLayer, OptionsLayer

from .specific_analysis import LanthiResults


def will_handle(products: List[str], _product_categories: Set[str]) -> bool:
    """ Returns true if one or more relevant products or product categories are present """
    return any(map(lambda p: p.startswith("lanthipeptide"), products))


def generate_html(region_layer: RegionLayer, results: LanthiResults,
                  _record_layer: RecordLayer, _options_layer: OptionsLayer
                  ) -> HTMLSections:
    """ Generates HTML output for the module """
    html = HTMLSections("lanthipeptides")

    comparison = results.comparippson_results
    motifs_by_locus = results.get_motifs_for_region(region_layer.region_feature)

    motifs_by_locus_by_core = {}
    for locus, motifs in motifs_by_locus.items():
        motifs_by_core = defaultdict(list)
        for motif in motifs:
            core = motif.core
            motifs_by_core[core].append(motif)
        motifs_by_locus_by_core[locus] = motifs_by_core

    detail_tooltip = Markup("<br>".join([
        "Shows the possible core peptides for each biosynthetic enzyme.",
        comparippson.get_tooltip_text(),
    ]))
    template = FileTemplate(path.get_full_path(__file__, "templates", "details.html"))
    html.add_detail_section("Lanthipeptides", template.render(results=motifs_by_locus_by_core, tooltip=detail_tooltip,
                                                              comparippson=comparison))

    side_tooltip = ("Lists the possible core peptides in the region. "
                    "Each core peptide lists the number of lanthionine bridges, possible molecular weights, "
                    "and the scores for cleavage site prediction and RODEO.")
    template = FileTemplate(path.get_full_path(__file__, "templates", "sidepanel.html"))
    html.add_sidepanel_section("Lanthipeptides", template.render(results=motifs_by_locus, tooltip=side_tooltip))

    return html
