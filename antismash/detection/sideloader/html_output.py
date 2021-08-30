# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Handles HTML output for the Type II PKS module """

from collections import defaultdict
from typing import List, Set

from antismash.common import path
from antismash.common.module_results import ModuleResults
from antismash.common.html_renderer import FileTemplate, HTMLSections
from antismash.common.layers import RegionLayer, RecordLayer, OptionsLayer

from .data_structures import SideloadedResults


def will_handle(_products: List[str], _product_categories: Set[str]) -> bool:
    """ Relevant to every region, so return True for every product """
    return True


def generate_html(region_layer: RegionLayer, results: ModuleResults,
                  _record_layer: RecordLayer, _options_layer: OptionsLayer) -> HTMLSections:
    """ Generates the HTML sections for all sideloaded annotations """
    assert isinstance(results, SideloadedResults)
    template = FileTemplate(path.get_full_path(__file__, "templates", "general.html"))
    tooltip_content = ("This annotation was made externally by %r and not by antiSMASH")

    html = HTMLSections("sideloaded")
    tools_by_name = {}
    areas_by_tool_name = defaultdict(list)
    for area in results.get_areas():
        if not region_layer.location.start <= area.start <= region_layer.location.end:
            continue
        areas_by_tool_name[area.tool.name].append(area)
        tools_by_name[area.tool.name] = area.tool
    for tool_name, areas in areas_by_tool_name.items():
        # avoid HTML class names containing spaces
        html_name = tool_name.replace(" ", "-")
        html.add_detail_section(tool_name,
                                template.render(
                                    tool=tools_by_name[tool_name],
                                    areas=areas,
                                    tooltip_content=tooltip_content % tool_name,
                                ),
                                html_name)
    return html
