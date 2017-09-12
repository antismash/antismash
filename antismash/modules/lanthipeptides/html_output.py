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
        self.core_peptides_predictions = generate_core_predictions(self.motifs)
        self.core_peptides_data = generate_core_data(self.motifs)

def generate_core_data(motifs):
    data = []
    for peptide in motifs:
        leader_seq = peptide.leader_seq
        core_seq = peptide.core_seq
        pred_class = peptide.peptide_class
        gene_id = peptide.get_name()
        data.append((leader_seq, core_seq, pred_class, gene_id))
    return data

def generate_core_predictions(motifs):
    predictions = []
    for core in motifs:
        gene_id = core.get_name()
        mass = "%.1f" % core.monoisotopic_mass
        mol_weight = "%.1f" % core.molecular_weight
        bridges = core.lan_bridges
        pred_class = core.peptide_class
        score = core.score
        rodeo_score = core.rodeo_score
        mods = core.get_modifications()
        alt_weights = core.alternative_weights
        predictions.append((gene_id, mass, mol_weight, bridges, pred_class, score, rodeo_score, alt_weights, mods))
    return predictions

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
    record = record_layer
    sidepanel = template.render(record=record,
                                cluster=cluster,
                                options=options_layer)
    if cluster.motifs:
        return sidepanel
    return ""
