# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Manages HTML construction for the Lanthipeptide module
"""

from typing import List

from antismash.common import path
from antismash.common.html_renderer import HTMLSections, FileTemplate
from antismash.common.layers import RegionLayer, RecordLayer, OptionsLayer

from .specific_analysis import LanthiResults


def will_handle(products: List[str]) -> bool:
    """ Returns True if the products provided are relevant to the module """
    return 'lanthipeptide' in products


def generate_html(region_layer: RegionLayer, results: LanthiResults,
                  _record_layer: RecordLayer, _options_layer: OptionsLayer
                  ) -> HTMLSections:
    """ Generates HTML output for the module """
    html = HTMLSections("lanthipeptides")

    detail_tooltip = ("Lists the possible core peptides for each biosynthetic enzyme, including the predicted class. "
                      "Each core peptide shows the leader and core peptide sequences, separated by a dash.")
    template = FileTemplate(path.get_full_path(__file__, "templates", "details.html"))
    motifs = results.get_motifs_for_region(region_layer.region_feature)
    html.add_detail_section("Lanthipeptides", template.render(results=motifs, tooltip=detail_tooltip))

    side_tooltip = ("Lists the possible core peptides in the region. "
                    "Each core peptide lists the number of lanthionine bridges, possible molecular weights, "
                    "and the scores for cleavage site prediction and RODEO.")
    template = FileTemplate(path.get_full_path(__file__, "templates", "sidepanel.html"))
    motifs = results.get_motifs_for_region(region_layer.region_feature)
    html.add_sidepanel_section("Lanthipeptides", template.render(results=motifs, tooltip=side_tooltip))

    return html
