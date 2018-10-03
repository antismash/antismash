# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Manages HTML construction for the Lanthipeptide module
"""

from typing import List

from jinja2 import FileSystemLoader, Environment, StrictUndefined

from antismash.common import path
from antismash.common.layers import RegionLayer, RecordLayer, OptionsLayer
from antismash.common.secmet import Prepeptide, Region

from .specific_analysis import LanthiResults


def will_handle(products: List[str]) -> bool:
    """ Returns True if the products provided are relevant to the module """
    return 'lanthipeptide' in products


class LanthipeptideLayer(RegionLayer):
    """ An extended RegionLayer for holding a list of LanthipeptideMotifs """
    def __init__(self, record: RecordLayer, region_feature: Region) -> None:
        RegionLayer.__init__(self, record, region_feature)
        self.motifs = []  # type: List[Prepeptide]
        for motif in self.record.seq_record.get_cds_motifs():
            if not isinstance(motif, Prepeptide):
                continue
            if not motif.is_contained_by(self.region_feature):
                continue
            if motif.peptide_class == "lanthipeptide":
                self.motifs.append(motif)


def generate_details_div(region_layer: RegionLayer, results: LanthiResults,
                         record_layer: RecordLayer, options_layer: OptionsLayer) -> str:
    """ Generates a HTML div for the main page of results """
    lanthi_layer = LanthipeptideLayer(record_layer, region_layer.region_feature)
    if not results:
        return ""
    env = Environment(loader=FileSystemLoader(path.get_full_path(__file__, "templates")),
                      autoescape=True, undefined=StrictUndefined)
    template = env.get_template('details.html')
    motifs = results.get_motifs_for_region(region_layer.region_feature)
    details_div = template.render(results=motifs)
    return details_div


def generate_sidepanel(region_layer: RegionLayer, results: LanthiResults,
                       record_layer: RecordLayer, options_layer: OptionsLayer) -> str:
    """ Generates a div for the sidepanel results """
    env = Environment(loader=FileSystemLoader(path.get_full_path(__file__, "templates")),
                      autoescape=True, undefined=StrictUndefined)
    template = env.get_template('sidepanel.html')
    region = LanthipeptideLayer(record_layer, region_layer.region_feature)
    if not results:
        return ""
    motifs = results.get_motifs_for_region(region_layer.region_feature)
    sidepanel = template.render(results=motifs)
    return sidepanel
