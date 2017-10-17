# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

import logging
import os
import re

from jinja2 import FileSystemLoader, Environment, StrictUndefined

from antismash.common import path
from antismash.common.layers import ClusterLayer, RecordLayer, OptionsLayer

from .results import NRPS_PKS_Results


def will_handle(products):
    return set(products).intersection({"nrps", "t1pks", "t2pks", "transatpks", "other", "otherks"})


def generate_js_domains(cluster, seq_record, results, options):
    assert isinstance(results, NRPS_PKS_Results), type(results)
    cluster_feature = seq_record.get_cluster(cluster['idx'])
    cluster = NrpspksLayer(results, cluster_feature, RecordLayer(seq_record, results, OptionsLayer(options)))
    return cluster.get_domain_details()


def generate_details_div(cluster_layer, results, record_layer, options_layer):
    env = Environment(loader=FileSystemLoader([path.get_full_path(__file__, 'templates')]),
                      autoescape=True, undefined=StrictUndefined)
    template = env.get_template('details.html')
    cluster = NrpspksLayer(results, cluster_layer.cluster_rec, record_layer)
    details_div = template.render(record=record_layer,
                                  cluster=cluster,
                                  options=options_layer)
    return details_div


def generate_sidepanel(cluster_layer, results, record_layer, options_layer):
    env = Environment(loader=FileSystemLoader([path.get_full_path(__file__, 'templates')]),
                      autoescape=True, undefined=StrictUndefined)
    template = env.get_template('sidepanel.html')
    cluster = NrpspksLayer(results, cluster_layer.cluster_rec, record_layer)
    sidepanel = template.render(record=record_layer,
                                cluster=cluster,
                                options=options_layer)
    return sidepanel


class NrpspksLayer(ClusterLayer):
    def __init__(self, results, cluster_rec, record):
        self.url_strict = {}  # gene name -> url
        self.url_relaxed = {}  # gene name -> url
        super().__init__(record, cluster_rec)
        self.transatpks = False
        assert isinstance(results, NRPS_PKS_Results), list(results)
        self.results = results
        self.monomer, self.used_domain_docking = results.cluster_predictions.get(cluster_rec.get_cluster_number(), ("N/A", False))

    def get_domain_details(self):
        js_domains = []
        js_cluster_domains = {'id': "cluster-%s-details" % self.idx, 'orfs': []}
        features = self.cluster_rec.cds_children

        for feature in features:
            if not feature.nrps_pks:
                continue

            sequence = feature.translation

            js_orf = {
                'id': feature.get_name(),
                'sequence': sequence,
                'domains': [],
            }

            for domain in feature.nrps_pks.domains:
                js_domain = self.parse_domain(domain, feature)
                if js_domain:
                    js_orf['domains'].append(js_domain)

            if js_orf['domains']:
                js_cluster_domains['orfs'].append(js_orf)

        if js_cluster_domains['orfs']:
            js_domains.append(js_cluster_domains)

        return js_domains

    @property
    def warning(self):
        core = ("Rough prediction of core scaffold based on assumed %s;"
                " tailoring reactions not taken into account")
        if self.used_domain_docking:
            detail = "PKS linker matching"
        else:
            detail = "PKS/NRPS colinearity"

        return core % detail

    @property
    def sidepanel_predictions(self):
        sidepanel_predictions = {}
        features = self.cluster_rec.cds_children
        for feature in features:
            if not feature.nrps_pks:
                continue

            per_CDS_predictions = []
            gene_id = feature.get_name()

            for domain in feature.nrps_pks.domains:
                per_Adomain_predictions = []
                preds = self.parse_substrate_predictions(domain.predictions)
                if not preds:
                    continue

                if domain.name == "AMP-binding":
                    for key, val in preds:
                        if key in ["SNN score", "%ID to nearest neighbour"]:
                            continue
                        if key == "PrediCAT":
                            val = val.rsplit("-")[-1]
                        values = self.filter_norine_as(val.split(","))
                        if values:
                            per_Adomain_predictions.extend(val.split(","))
                    per_Adomains_predictions_unique = list(set(per_Adomain_predictions))
                    per_CDS_predictions.append(per_Adomains_predictions_unique)

                if gene_id not in sidepanel_predictions:
                    sidepanel_predictions[gene_id] = []
                sidepanel_predictions[gene_id].append(preds)

            if per_CDS_predictions:
                self.url_strict[gene_id] = self.get_norine_url_for_specArray(per_CDS_predictions)
                if self.url_strict[gene_id]:
                    self.url_relaxed[gene_id] = self.get_norine_url_for_specArray(per_CDS_predictions, be_strict=False)

        return sidepanel_predictions

    def is_nrps(self):
        return 'nrps' in self.cluster_rec.products

    def parse_substrate_predictions(self, nrps_predictions):
        "Parse the substrate predictions from the NRPS/PKS domain string"
        predictions = []

        for method, pred in sorted(nrps_predictions.items()):
            if method == "SNN score":
                pass
            elif method == "PID to NN":
                method = "%ID to nearest neighbour"
            elif pred == "no_call":
                pred = "N/A"
            predictions.append((method, pred))

        return predictions

    def get_structure_image_url(self):
        "Get the relative url to the structure image"
        expected = os.path.join("structures",
                                "genecluster%d.png" % self.cluster_rec.get_cluster_number())
        abs_path = os.path.join(self.record.options.output_dir, expected)
        if os.path.exists(abs_path):
            return expected
        return 'images/nostructure_icon.png'

    def get_monomer_prediction(self):
        "Get the monomer prediction of the cluster"
        return self.monomer

    def map_as_names_to_norine(self, as_name):
        """ Just a dictionary helper function to map antiSMASH amino acid nomenclature to NORINE"""

        as_replacement_dict = {'bht': 'bOH-Tyr',
                               'dhb': 'diOH-Bz',
                               'iva': 'Ival',
                               'pip': 'Hpr',
                               'sal': 'diOH-Bz',
                               'nrp': 'X',
                               'dpg': 'dhpg'}  # TODO: verify this one
        return as_replacement_dict.get(as_name, as_name)

    def filter_norine_as(self, as_list, be_strict=False):
        """ Remove PKS and unknown substrate predictions
            use be_strict = True to also filter nrp/X
        """
        filtered_list = []
        bad_monomers = {'pk', 'N/A', 'hydrophilic', 'hydrophobic', 'mal', 'mmal'}
        if be_strict:
            bad_monomers = bad_monomers.union({'nrp', 'X'})
        for monomer in as_list:
            for sub_monomer in monomer.split("|"):  # if no pipe, list of 1
                if sub_monomer in bad_monomers:
                    continue
                filtered_list.append(self.map_as_names_to_norine(sub_monomer))
        return filtered_list

    def get_norine_url_for_cluster(self, be_strict=True):
        """ Get a NORINE URL string for direct querying
            use be_strict=False to add * after each monomer"""

        monomer_string = self.get_monomer_prediction()
        monomers_per_protein_list = re.findall("\\(.*?\\)", monomer_string)
        i = 1
        nrpslist = []
        for monomers_per_protein in monomers_per_protein_list:
            monomers = monomers_per_protein[1:-1].split("-")

            if be_strict:
                monomers = [self.map_as_names_to_norine(element.lower()) for element in self.filter_norine_as(monomers, be_strict=True)]
            else:
                monomers = [self.map_as_names_to_norine(element.lower()) + "*" for element in self.filter_norine_as(monomers)]
            # Norine doesn't allow "x*" as a "relaxed" monomer, so we have to replace this with "x"
            monomers = list(map(lambda x: "x" if x == "x*" else x, monomers))

            if monomers:
                nrpslist.append("nrps" + str(i) + "=" + ",".join(monomers))
            i += 1
        urlstring = "http://bioinfo.lifl.fr/norine/fingerPrintSearch.jsp?"+"&".join(nrpslist)
        return urlstring

    def get_norine_url_for_specArray(self, specArray, be_strict=True):
        """ generate NORINE URL string for direct querying from array of specificity predictions
            use be_strict=False to add * after each monomer
        """
        modulelist = []
        if specArray:
            for domain_specificity_list in specArray:
                if not domain_specificity_list:
                    logging.error("retrieved empty domain list for assembling per protein NORINE link string")
                    raise RuntimeError("empty domain list for NORINE generation")
                if len(domain_specificity_list) == 1:
                    if be_strict:
                        modulelist.append(self.map_as_names_to_norine(domain_specificity_list[0]))
                    else:
                        modulelist.append(self.map_as_names_to_norine(domain_specificity_list[0])+"*")
                else:  # > 1
                    if be_strict:
                        modulelist.append("["+"|".join([self.map_as_names_to_norine(element) for element in self.filter_norine_as(domain_specificity_list, be_strict=True)])+"]")
                    else:
                        # we have to use be_strict to remove X from the list of predictions, as otherwise consenus: nrp will always match
                        monomers = [self.map_as_names_to_norine(element)+"*" for element in self.filter_norine_as(domain_specificity_list, be_strict=True)]
                        # Norine doesn't allow "x*" as a "relaxed" monomer, so we have to replace this with "x"
                        monomers = map(lambda x: "x" if x == "x*" else x, monomers)
                        modulelist.append("["+"|".join(monomers)+"]")
            urlstring = "http://bioinfo.lifl.fr/norine/fingerPrintSearch.jsp?nrps1="+",".join(modulelist)
            return urlstring

    def parse_domain(self, domain, feature):
        "Convert a NRPS/PKS domain string to a dict useable by json.dumps"
        predictions = self.parse_substrate_predictions(domain.predictions)

        # Create url_link to NaPDoS for C and KS domains
        napdoslink = ""
        domainseq = str(feature.get_aa_sequence())[domain.start:domain.end]
        dna_sequence = ''
        if self.record.options.input_type == 'nucl':
            dna_sequence = str(feature.extract(self.record.seq_record.seq))
        base = "http://napdos.ucsd.edu/cgi-bin/process_request.cgi?query_type=aa&amp;ref_seq_file=all_{0}_public_12062011.faa&amp;Sequence=%3E{0}_domain_from_antiSMASH%0D{1}"
        if domain.name == "PKS_KS":
            napdoslink = base.format("KS", domainseq)
        elif "Condensation" in domain.name:
            napdoslink = base.format("C", domainseq)
        blastlink = "http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE=Proteins&amp;PROGRAM=blastp&amp;BLAST_PROGRAMS=blastp&amp;QUERY=" + domainseq + "&amp;LINK_LOC=protein&amp;PAGE_TYPE=BlastSearch"

        js_domain = {'type': domain.name, 'start': domain.start,
                     'end': domain.end, 'predictions': predictions,
                     'napdoslink': napdoslink, 'blastlink': blastlink,
                     'sequence': domainseq, 'dna_sequence': dna_sequence}
        return js_domain
