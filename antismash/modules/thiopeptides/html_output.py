# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" HTML generation for the thiopeptides module """

from typing import List

from jinja2 import FileSystemLoader, Environment, StrictUndefined

from antismash.common import path
from antismash.common.layers import ClusterLayer, RecordLayer, OptionsLayer
from antismash.common.secmet import Cluster

from .specific_analysis import ThioResults, ThiopeptideMotif


def will_handle(products: List[str]) -> bool:
    """ HTML generation only occurs if this function reutrns True """
    return 'thiopeptide' in products


class ThiopeptideLayer(ClusterLayer):
    """ A wrapper of ClusterLayer to allow for tracking the ThiopeptideMotifs """
    def __init__(self, record: RecordLayer, results: ThioResults, cluster_feature: Cluster) -> None:
        ClusterLayer.__init__(self, record, cluster_feature)
        self.motifs = []  # type: List[ThiopeptideMotif]
        for motif in results.motifs:
            if motif.is_contained_by(self.cluster_feature):
                self.motifs.append(motif)


def generate_details_div(cluster_layer: ClusterLayer, results: ThioResults,
                         record_layer: RecordLayer, options_layer: OptionsLayer) -> str:
    """ Generates the HTML details section from the ThioResults instance """
    env = Environment(loader=FileSystemLoader(path.get_full_path(__file__, "templates")),
                      autoescape=True, undefined=StrictUndefined)
    template = env.get_template('details.html')
    details_div = template.render(record=record_layer,
                                  cluster=ThiopeptideLayer(record_layer, results, cluster_layer.cluster_feature),
                                  options=options_layer)
    return details_div


def generate_sidepanel(cluster_layer: ClusterLayer, results: ThioResults,
                       record_layer: RecordLayer, options_layer: OptionsLayer) -> str:
    """ Generates the HTML sidepanel section from the ThioResults instance """
    env = Environment(loader=FileSystemLoader(path.get_full_path(__file__, "templates")),
                      autoescape=True, undefined=StrictUndefined)
    template = env.get_template('sidepanel.html')
    cluster = ThiopeptideLayer(record_layer, results, cluster_layer.cluster_feature)
    record = record_layer
    sidepanel = template.render(record=record,
                                cluster=cluster,
                                options=options_layer)
    if cluster.motifs:
        return sidepanel
    return ""
