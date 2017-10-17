# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

import logging

from jinja2 import FileSystemLoader, Environment, StrictUndefined

from antismash.common import path
from antismash.common.layers import ClusterLayer
from antismash.common.secmet import Prepeptide


def will_handle(products):
    return 'lanthipeptide' in products


class LanthipeptideLayer(ClusterLayer):
    def __init__(self, record, cluster_feature):
        ClusterLayer.__init__(self, record, cluster_feature)
        self.motifs = []
        for motif in self.record.seq_record.get_cds_motifs():
            if not isinstance(motif, Prepeptide):
                continue
            if not motif.is_contained_by(self.cluster_rec):
                continue
            if motif.peptide_type == "lanthipeptide":
                self.motifs.append(motif)


def generate_details_div(cluster_layer, results, record_layer, options_layer):
    lanthi_layer = LanthipeptideLayer(record_layer, cluster_layer.cluster_rec)
    if not (not options_layer.minimal or options_layer.lanthipeptides_enabled
            or lanthi_layer.motifs):
        return ""
    env = Environment(loader=FileSystemLoader(path.get_full_path(__file__, "templates")),
                      autoescape=True, undefined=StrictUndefined)
    template = env.get_template('details.html')
    details_div = template.render(record=record_layer,
                                  cluster=lanthi_layer,
                                  options=options_layer)
    return details_div


def generate_sidepanel(cluster_layer, results, record_layer, options_layer):
    env = Environment(
        loader=FileSystemLoader(path.get_full_path(__file__, "templates")),
        autoescape=True, undefined=StrictUndefined)
    template = env.get_template('sidepanel.html')
    cluster = LanthipeptideLayer(record_layer, cluster_layer.cluster_rec)
    if not cluster.motifs:
        return ""
    record = record_layer
    sidepanel = template.render(record=record,
                                cluster=cluster,
                                options=options_layer)
    return sidepanel
