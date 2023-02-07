# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Handles HTML generation for the clusterblast variants """

from typing import Any, List, Set

from antismash.common import path
from antismash.common.html_renderer import HTMLSections, FileTemplate, Markup
from antismash.common.layers import RegionLayer, RecordLayer, OptionsLayer

from .results import ClusterBlastResults


def will_handle(_products: List[str], _product_categories: Set[str]) -> bool:
    """ Relevant to every region, so return True for every product """
    return True


def generate_html(region_layer: RegionLayer, results: ClusterBlastResults,
                  record_layer: RecordLayer, options_layer: OptionsLayer
                  ) -> HTMLSections:
    """ Generates the HTML sections for all variants of clusterblast
    """

    html = HTMLSections("clusterblast")
    region = region_layer.region_feature

    base_tooltip = ("Shows %s that are similar to the current region. Genes marked with the "
                    "same colour are interrelated. White genes have no relationship.<br>"
                    "Click on reference genes to show details of similarities to "
                    "genes within the current region.")

    if options_layer.cb_general or region.clusterblast is not None:
        tooltip = base_tooltip % "regions from the antiSMASH database"
        tooltip += "<br>Click on an accession to open that entry in the antiSMASH database (if applicable)."
        div = generate_div(region_layer, record_layer, options_layer, "clusterblast", tooltip)
        html.add_detail_section("ClusterBlast", div, "clusterblast")
    if options_layer.cb_knownclusters or region.knownclusterblast is not None:
        assert results.knowncluster and results.knowncluster.data_version, "missing MIBiG data version"
        additional = {
            "version": results.knowncluster.data_version,
        }
        tooltip = base_tooltip % "clusters from the MiBIG database"
        tooltip += "<br>Click on an accession to open that entry in the MiBIG database."
        div = generate_div(region_layer, record_layer, options_layer, "knownclusterblast", tooltip,
                           additional)
        html.add_detail_section("KnownClusterBlast", div, "knownclusterblast")
    if options_layer.cb_subclusters or region.subclusterblast is not None:
        tooltip = base_tooltip % "sub-cluster units"
        div = generate_div(region_layer, record_layer, options_layer, "subclusterblast", tooltip)
        html.add_detail_section("SubClusterBlast", div, "subclusterblast")

    return html


def generate_div(region_layer: RegionLayer, record_layer: RecordLayer,
                 options_layer: OptionsLayer, search_type: str,
                 tooltip: str, additional: dict[str, Any] = None) -> Markup:
    """ Generates the specific HTML section of the body for a given variant of
        clusterblast
    """
    if additional is None:
        additional = {}
    template = FileTemplate(path.get_full_path(__file__, "templates", f"{search_type}.html"))
    return template.render(record=record_layer, region=region_layer, options=options_layer,
                           tooltip=tooltip, **additional)
