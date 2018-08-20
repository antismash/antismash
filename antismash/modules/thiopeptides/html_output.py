# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" HTML generation for the thiopeptides module """

from typing import List

from jinja2 import FileSystemLoader, Environment, StrictUndefined

from antismash.common import path
from antismash.common.layers import RegionLayer, RecordLayer, OptionsLayer
from antismash.common.secmet import Region

from .specific_analysis import ThioResults, ThiopeptideMotif


def will_handle(products: List[str]) -> bool:
    """ HTML generation only occurs if this function reutrns True """
    return 'thiopeptide' in products


class ThiopeptideLayer(RegionLayer):
    """ A wrapper of RegionLayer to allow for tracking the ThiopeptideMotifs """
    def __init__(self, record: RecordLayer, results: ThioResults, region_feature: Region) -> None:
        RegionLayer.__init__(self, record, region_feature)
        self.motifs = []  # type: List[ThiopeptideMotif]
        for motif in results.motifs:
            if motif.is_contained_by(self.region_feature) and isinstance(motif, ThiopeptideMotif):
                self.motifs.append(motif)


def generate_details_div(region_layer: RegionLayer, results: ThioResults,
                         record_layer: RecordLayer, options_layer: OptionsLayer) -> str:
    """ Generates the HTML details section from the ThioResults instance """
    env = Environment(loader=FileSystemLoader(path.get_full_path(__file__, "templates")),
                      autoescape=True, undefined=StrictUndefined)
    template = env.get_template('details.html')
    details_div = template.render(record=record_layer,
                                  cluster=ThiopeptideLayer(record_layer, results, region_layer.region_feature),
                                  options=options_layer)
    return details_div


def generate_sidepanel(region_layer: RegionLayer, results: ThioResults,
                       record_layer: RecordLayer, options_layer: OptionsLayer) -> str:
    """ Generates the HTML sidepanel section from the ThioResults instance """
    env = Environment(loader=FileSystemLoader(path.get_full_path(__file__, "templates")),
                      autoescape=True, undefined=StrictUndefined)
    template = env.get_template('sidepanel.html')
    cluster = ThiopeptideLayer(record_layer, results, region_layer.region_feature)
    record = record_layer
    sidepanel = template.render(record=record,
                                cluster=cluster,
                                options=options_layer)
    if cluster.motifs:
        return sidepanel
    return ""
