# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Handles HTML output for the Type II PKS module """

from typing import List, Set

from antismash.common import path
from antismash.common.html_renderer import FileTemplate, HTMLSections, docs_link
from antismash.common.layers import RegionLayer, RecordLayer, OptionsLayer

from .results import T2PKSResults


def will_handle(products: List[str], _product_categories: Set[str]) -> bool:
    """ Returns true if one or more relevant products or product categories are present """
    return "T2PKS" in products


def generate_html(region_layer: RegionLayer, results: T2PKSResults,
                  _record_layer: RecordLayer, _options_layer: OptionsLayer) -> HTMLSections:
    """ Generate the sidepanel HTML with results from the type II PKS module """
    html = HTMLSections("t2pks")

    predictions = []
    for cluster in region_layer.get_unique_protoclusters():
        if cluster.product == "T2PKS":
            predictions.append(results.cluster_predictions[cluster.get_protocluster_number()])

    template = FileTemplate(path.get_full_path(__file__, "templates", "sidepanel.html"))

    tooltip_content = (
        "Predictions of starter units, elongations, product classes, "
        "and potential molecular weights for type II PKS clusters."
        f"<br>More detailed information is available {docs_link('here', 'modules/t2pks')}."
    )
    html.add_sidepanel_section("Type II PKS", template.render(predictions=predictions,
                                                              tooltip_content=tooltip_content))

    return html
