# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Manages HTML construction for the Atropopeptide module
"""

from antismash.common import comparippson, path
from antismash.common.html_renderer import HTMLSections, FileTemplate, Markup
from antismash.common.layers import RecordLayer, RegionLayer, OptionsLayer

from .specific_analysis import AtropoResults, Prepeptide


def will_handle(products: list[str], _product_categories: set[str]) -> bool:
    """ Returns true if one or more relevant products or product categories are present """
    return 'atropopeptide' in products


def generate_html(region_layer: RegionLayer, results: AtropoResults,
                  record_layer: RecordLayer, _options_layer: OptionsLayer) -> HTMLSections:
    """ Generates HTML for the module """
    html = HTMLSections("atropopeptides")

    loci = []
    for locus in results.motifs_by_locus:
        if record_layer.get_cds_by_name(locus).is_contained_by(region_layer.region_feature):
            loci.append(locus)

    by_locus: dict[str, dict[str, list[Prepeptide]]] = {}
    for locus in loci:
        motifs = results.motifs_by_locus[locus]
        motifs_by_core: dict[str, list[Prepeptide]] = {motif.core: [] for motif in motifs}
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
    html.add_detail_section("Atropopeptides", section)

    return html

