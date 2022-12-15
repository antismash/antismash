# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" HTML generation for the thiopeptides module """

from typing import Dict, List, Set

from antismash.common import comparippson, path
from antismash.common.html_renderer import HTMLSections, FileTemplate, Markup
from antismash.common.layers import RegionLayer, RecordLayer, OptionsLayer
from antismash.common.secmet import Prepeptide

from .specific_analysis import ThioResults


def will_handle(products: List[str], _product_categories: Set[str]) -> bool:
    """ Returns true if one or more relevant products or product categories are present """
    return 'thiopeptide' in products


def generate_html(region_layer: RegionLayer, results: ThioResults,
                  record_layer: RecordLayer, options_layer: OptionsLayer
                  ) -> HTMLSections:
    """ Generates HTML for the module """
    html = HTMLSections("thiopeptides")

    if not results:
        return html

    motifs = results.get_motifs_for_region(region_layer.region_feature)
    motif_groups: Dict[str, List[Prepeptide]] = {motif.core: [] for motif in motifs}
    for motif in motifs:
        motif_groups[motif.core].append(motif)

    detail_tooltip = Markup("<br>".join([
        (
            "Lists the possible core peptides for each biosynthetic enzyme. "
            "Predicted tail sequences are also shown, if present. "
        ),
        comparippson.get_tooltip_text(),
    ]))
    template = FileTemplate(path.get_full_path(__file__, "templates", "details.html"))
    details = template.render(record=record_layer,
                              motif_groups=motif_groups.values(),
                              comparippson_results=results.comparippson_results,
                              options=options_layer,
                              tooltip=detail_tooltip)
    html.add_detail_section("Thiopeptides", details)

    side_tooltip = ("Lists the possible core peptides in the region. "
                    "Each core peptide lists its possible molecular weights "
                    "and the scores for cleavage site prediction and RODEO. "
                    "If relevant, other features, such as macrocycle and amidation, will also be listed.")
    template = FileTemplate(path.get_full_path(__file__, "templates", "sidepanel.html"))
    sidepanel = template.render(record=record_layer,
                                motifs=motifs,
                                options=options_layer,
                                tooltip=side_tooltip)
    html.add_sidepanel_section("Thiopeptides", sidepanel)
    return html
