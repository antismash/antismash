# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

import logging

from jinja2 import FileSystemLoader, Environment, StrictUndefined

from antismash.common.layers import ClusterLayer
from antismash.common.secmet import Prepeptide

def will_handle(product):
    return product.find('lanthipeptide') > -1

class LanthipeptideLayer(ClusterLayer):
    def __init__(self, cluster, record, cluster_feature):
        ClusterLayer.__init__(self, cluster, record, cluster_feature)
        self.motifs = []
        for motif in self.record.seq_record.get_cds_motifs():
            if not isinstance(motif, Prepeptide):
                continue
            if motif.location.start < self.start or \
               motif.location.end > self.end: #TODO cleanup repeated cases
                continue
            if motif.peptide_type == "lanthipeptide":
                self.motifs.append(motif)

def generate_details_div(cluster_layer, record_layer, options_layer):
    env = Environment(
        loader=FileSystemLoader(['antismash/modules/lanthipeptides/templates']),
        autoescape=True, undefined=StrictUndefined)
    template = env.get_template('details.html')
    details_div = template.render(record=record_layer,
                           cluster=LanthipeptideLayer(cluster_layer.cluster, record_layer, cluster_layer.cluster_rec),
                           options=options_layer)
    return details_div


def generate_sidepanel(cluster_layer, record_layer, options_layer):
    env = Environment(
        loader=FileSystemLoader(['antismash/modules/lanthipeptides/templates']),
        autoescape=True, undefined=StrictUndefined)
    template = env.get_template('sidepanel.html')
    cluster = LanthipeptideLayer(cluster_layer.cluster, record_layer, cluster_layer.cluster_rec)
    if not cluster.motifs:
        return ""
    record = record_layer
    sidepanel = template.render(record=record,
                                cluster=cluster,
                                options=options_layer)
    return sidepanel
