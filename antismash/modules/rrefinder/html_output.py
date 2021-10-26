# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Manages HTML construction for the RREFinder module
"""

from typing import List, Set

from antismash.common import path
from antismash.common.html_renderer import HTMLSections, FileTemplate
from antismash.common.layers import RegionLayer, RecordLayer, OptionsLayer

from .rrefinder import RREFinderResults


def will_handle(_products: List[str], product_categories: Set[str]) -> bool:
    """ Returns true if one or more relevant products or product categories are present """
    return "RiPP" in product_categories


def generate_html(region_layer: RegionLayer, results: RREFinderResults,
                  _record_layer: RecordLayer, _options_layer: OptionsLayer
                  ) -> HTMLSections:
    """ Generates HTML output for the module """
    html = HTMLSections("rrefinder")

    side_tooltip = ("RREfinder results sidepanel.")
    template = FileTemplate(path.get_full_path(__file__, "templates", "sidepanel.html"))

    protoclusters = []
    for proto in region_layer.get_unique_protoclusters():
        if proto.get_protocluster_number() in results.hits_by_protocluster:
            protoclusters.append(proto)

    if protoclusters:
        section = template.render(results=results, protoclusters=protoclusters, tooltip=side_tooltip)
        html.add_sidepanel_section("RREFinder", section, class_name="RREfinder")

    return html
