# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Manages HTML construction for the RREFinder module
"""

from typing import List

from antismash.common import path
from antismash.common.html_renderer import HTMLSections, FileTemplate
from antismash.common.layers import RegionLayer, RecordLayer, OptionsLayer

from .rrefinder import is_ripp, RREFinderResults


def will_handle(products: List[str]) -> bool:
    """ Returns True if the products provided are relevant to the module """
    return any(is_ripp(product) for product in products)


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
