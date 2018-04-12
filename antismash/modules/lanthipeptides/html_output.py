# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Manages HTML construction for the Lanthipeptide module
"""

from typing import List

from jinja2 import FileSystemLoader, Environment, StrictUndefined

from antismash.common import path
from antismash.common.layers import ClusterLayer, RecordLayer, OptionsLayer
from antismash.common.secmet import Prepeptide, Cluster

from .specific_analysis import LanthiResults


def will_handle(products: List[str]) -> bool:
    """ Returns True if the products provided are relevant to the module """
    return 'lanthipeptide' in products


class LanthipeptideLayer(ClusterLayer):
    """ An extended ClusterLayer for holding a list of LanthipeptideMotifs """
    def __init__(self, record: RecordLayer, cluster_feature: Cluster) -> None:
        ClusterLayer.__init__(self, record, cluster_feature)
        self.motifs = []  # type: List[Prepeptide]
        for motif in self.record.seq_record.get_cds_motifs():
            if not isinstance(motif, Prepeptide):
                continue
            if not motif.is_contained_by(self.cluster_feature):
                continue
            if motif.peptide_class == "lanthipeptide":
                self.motifs.append(motif)


def generate_details_div(cluster_layer: ClusterLayer, results: LanthiResults,
                         record_layer: RecordLayer, options_layer: OptionsLayer) -> str:
    """ Generates a HTML div for the main page of results """
    lanthi_layer = LanthipeptideLayer(record_layer, cluster_layer.cluster_feature)
    if not results:
        return ""
    env = Environment(loader=FileSystemLoader(path.get_full_path(__file__, "templates")),
                      autoescape=True, undefined=StrictUndefined)
    template = env.get_template('details.html')
    motifs_in_cluster = {}
    for locus in results.clusters.get(cluster_layer.get_cluster_number(), []):
        motifs_in_cluster[locus] = results.motifs_by_locus[locus]
    details_div = template.render(record=record_layer,
                                  cluster=lanthi_layer,
                                  options=options_layer,
                                  results=motifs_in_cluster)
    return details_div


def generate_sidepanel(cluster_layer: ClusterLayer, results: LanthiResults,
                       record_layer: RecordLayer, options_layer: OptionsLayer) -> str:
    """ Generates a div for the sidepanel results """
    env = Environment(loader=FileSystemLoader(path.get_full_path(__file__, "templates")),
                      autoescape=True, undefined=StrictUndefined)
    template = env.get_template('sidepanel.html')
    cluster = LanthipeptideLayer(record_layer, cluster_layer.cluster_feature)
    if not results:
        return ""
    motifs_in_cluster = {}
    for locus in results.clusters.get(cluster_layer.get_cluster_number(), []):
        motifs_in_cluster[locus] = results.motifs_by_locus[locus]
    sidepanel = template.render(record=record_layer,
                                cluster=cluster,
                                options=options_layer,
                                results=motifs_in_cluster)
    return sidepanel
