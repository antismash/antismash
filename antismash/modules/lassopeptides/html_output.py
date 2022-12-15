# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Manages HTML construction for the Lassopeptide module
"""

from typing import Dict, List, Set

from antismash.common import comparippson, path
from antismash.common.html_renderer import HTMLSections, FileTemplate, Markup
from antismash.common.layers import RecordLayer, RegionLayer, OptionsLayer

from .specific_analysis import LassoResults, Prepeptide


def will_handle(products: List[str], _product_categories: Set[str]) -> bool:
    """ Returns true if one or more relevant products or product categories are present """
    return 'lassopeptide' in products


def generate_html(region_layer: RegionLayer, results: LassoResults,
                  record_layer: RecordLayer, _options_layer: OptionsLayer) -> HTMLSections:
    """ Generates HTML for the module """
    html = HTMLSections("lassopeptides")

    loci = []
    for locus in results.motifs_by_locus:
        if record_layer.get_cds_by_name(locus).is_contained_by(region_layer.region_feature):
            loci.append(locus)

    by_locus: Dict[str, Dict[str, List[Prepeptide]]] = {}
    for locus in loci:
        motifs = results.motifs_by_locus[locus]
        motifs_by_core: Dict[str, List[Prepeptide]] = {motif.core: [] for motif in motifs}
        for motif in motifs:
            motifs_by_core[motif.core].append(motif)
        by_locus[locus] = motifs_by_core

    detail_tooltip = Markup("<br>".join([
        (
            "Lists the possible core peptides for each biosynthetic enzyme. "
            "Each core peptide shows the leader and core peptide sequences. "
        ),
        comparippson.get_tooltip_text(),
    ]))

    template = FileTemplate(path.get_full_path(__file__, "templates", "details.html"))
    section = template.render(motifs_by_locus=by_locus, tooltip=detail_tooltip,
                              comparippson_results=results.comparippson_results)
    html.add_detail_section("Lasso peptides", section)

    side_tooltip = ("Lists the possible core peptides in the region. "
                    "Each core peptide lists the number of disulfide bridges, possible molecular weights, "
                    "and the scores for cleavage site prediction and RODEO.")
    template = FileTemplate(path.get_full_path(__file__, "templates", "sidepanel.html"))
    html.add_sidepanel_section("Lasso peptides", template.render(motifs_by_locus=by_locus, tooltip=side_tooltip))

    return html
