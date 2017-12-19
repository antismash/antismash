# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

from jinja2 import FileSystemLoader, Environment, StrictUndefined

from antismash.common import path
from antismash.common.secmet import Prepeptide
from antismash.common.layers import ClusterLayer


def will_handle(products):
    return 'sactipeptide' in products


def generate_details_div(cluster_layer, results, record_layer, options_layer):
    env = Environment(loader=FileSystemLoader([path.get_full_path(__file__, 'templates')]),
                      autoescape=True, undefined=StrictUndefined)
    template = env.get_template('details.html')
    motifs_in_cluster = {}
    for locus in results.clusters.get(cluster_layer.get_cluster_number(), []):
        motifs_in_cluster[locus] = results.motifs_by_locus[locus]
    details_div = template.render(record=record_layer,
                                  cluster=SactipeptideLayer(record_layer, cluster_layer.cluster_rec),
                                  options=options_layer,
                                  results=motifs_in_cluster)
    return details_div


def generate_sidepanel(cluster_layer, results, record_layer, options_layer):
    env = Environment(loader=FileSystemLoader([path.get_full_path(__file__, 'templates')]),
                      autoescape=True, undefined=StrictUndefined)
    template = env.get_template('sidepanel.html')
    cluster = SactipeptideLayer(record_layer, cluster_layer.cluster_rec)
    motifs_in_cluster = {}
    for locus in results.clusters.get(cluster_layer.get_cluster_number(), []):
        motifs_in_cluster[locus] = results.motifs_by_locus[locus]
    sidepanel = template.render(record=record_layer,
                                cluster=cluster,
                                options=options_layer,
                                results=motifs_in_cluster)
    return sidepanel


class SactipeptideLayer(ClusterLayer):
    """ An extended ClusterLayer for holding a list of LanthipeptideMotifs """
    def __init__(self, record, cluster_feature):
        ClusterLayer.__init__(self, record, cluster_feature)
        self.motifs = []
        for motif in self.record.seq_record.get_cds_motifs():
            if not isinstance(motif, Prepeptide):
                continue
            if not motif.is_contained_by(self.cluster_rec):
                continue
            if motif.peptide_class == "sactipeptide":
                self.motifs.append(motif)
