# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Manages HTML construction for the Lanthipeptide module
"""

from typing import List

from jinja2 import FileSystemLoader, Environment, StrictUndefined

from antismash.common import path
from antismash.common.layers import RegionLayer, RecordLayer, OptionsLayer

from .specific_analysis import LanthiResults


def will_handle(products: List[str]) -> bool:
    """ Returns True if the products provided are relevant to the module """
    return 'lanthipeptide' in products


def generate_details_div(region_layer: RegionLayer, results: LanthiResults,
                         _record_layer: RecordLayer, _options_layer: OptionsLayer) -> str:
    """ Generates a HTML div for the main page of results """
    if not results:
        return ""
    env = Environment(loader=FileSystemLoader(path.get_full_path(__file__, "templates")),
                      autoescape=True, undefined=StrictUndefined)
    template = env.get_template('details.html')
    motifs = results.get_motifs_for_region(region_layer.region_feature)
    details_div = template.render(results=motifs)
    return details_div


def generate_sidepanel(region_layer: RegionLayer, results: LanthiResults,
                       _record_layer: RecordLayer, _options_layer: OptionsLayer) -> str:
    """ Generates a div for the sidepanel results """
    env = Environment(loader=FileSystemLoader(path.get_full_path(__file__, "templates")),
                      autoescape=True, undefined=StrictUndefined)
    template = env.get_template('sidepanel.html')
    if not results:
        return ""
    motifs = results.get_motifs_for_region(region_layer.region_feature)
    sidepanel = template.render(results=motifs)
    return sidepanel
