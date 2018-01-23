# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Generates HTML and JSON for the nrps_pks module """


import logging
import os
import re
from typing import Dict, List, Optional, Tuple, Union

from jinja2 import FileSystemLoader, Environment, StrictUndefined

from antismash.common import path
from antismash.common.layers import ClusterLayer, RecordLayer, OptionsLayer
from antismash.common.secmet import CDSFeature
from antismash.common.secmet.qualifiers import NRPSPKSQualifier

from .results import NRPS_PKS_Results


class JSONBase(dict):
    """ A base class for JSON-serialisable objects """
    def __init__(self, keys):
        super().__init__()
        self._keys = keys

    def __getitem__(self, key):
        return getattr(self, key)

    def items(self):
        for key in self._keys:
            yield (key, getattr(self, key))

    def values(self):
        for key in self._keys:
            yield getattr(self, key)

    def __len__(self):
        return len(self._keys)


class JSONDomain(JSONBase):
    """ A JSON-serialisable object for simplifying domain datatypes throughout this file """
    def __init__(self, domain, predictions, napdos_link, blast_link, sequence, dna):
        super().__init__(['type', 'start', 'end', 'predictions', 'napdoslink',
                      'blastlink', 'sequence', 'dna_sequence'])
        self.type = str(domain.name)
        self.start = int(domain.start)
        self.end = int(domain.end)
        self.predictions = predictions
        self.napdoslink = str(napdos_link)
        self.blastlink = str(blast_link)
        self.sequence = str(sequence)
        self.dna_sequence = str(dna)


class JSONOrf(JSONBase):
    """ A JSON-serialisable object for simplifying ORF datatypes throughout this file """
    def __init__(self, feature):
        super().__init__(['id', 'sequence', 'domains'])
        self.sequence = feature.translation
        self.id = feature.get_name()
        self.domains = []

    def add_domain(self, domain: JSONDomain) -> None:
        """ Add a JSONDomain to the list of domains in this ORF """
        assert isinstance(domain, JSONDomain)
        self.domains.append(domain)


def will_handle(products: List[str]) -> bool:
    """ Returns true if one or more relevant products are present """
    return bool(set(products).intersection({"nrps", "t1pks", "t2pks", "transatpks",
                                            "nrpsfragment", "otherks"}))


def generate_js_domains(cluster, record, results, options
                        ) -> Optional[Dict[str, Union[str, List[JSONDomain]]]]:
    """ Generate """
    assert isinstance(results, NRPS_PKS_Results), type(results)
    cluster_feature = record.get_cluster(cluster['idx'])
    cluster = NrpspksLayer(results, cluster_feature, RecordLayer(record, results, OptionsLayer(options)))
    return cluster.get_domain_details()


def generate_details_div(cluster_layer, results, record_layer, options_layer) -> str:
    """ Generate the main HTML with results from the NRPS/PKS module """
    env = Environment(loader=FileSystemLoader(path.get_full_path(__file__, 'templates')),
                      autoescape=True, undefined=StrictUndefined)
    template = env.get_template('details.html')
    cluster = NrpspksLayer(results, cluster_layer.cluster_rec, record_layer)
    details_div = template.render(record=record_layer,
                                  cluster=cluster,
                                  options=options_layer)
    return details_div


def generate_sidepanel(cluster_layer, results, record_layer, options_layer) -> str:
    """ Generate the sidepanel HTML with results from the NRPS/PKS module """
    env = Environment(loader=FileSystemLoader(path.get_full_path(__file__, 'templates')),
                      autoescape=True, undefined=StrictUndefined)
    template = env.get_template('sidepanel.html')
    cluster = NrpspksLayer(results, cluster_layer.cluster_rec, record_layer)
    sidepanel = template.render(record=record_layer,
                                cluster=cluster,
                                options=options_layer)
    return sidepanel


def map_as_names_to_norine(as_name: str) -> str:
    """ Maps antiSMASH amino acid nomenclature to NORINE """

    as_replacement_dict = {'bht': 'bOH-Tyr',
                           'dhb': 'diOH-Bz',
                           'iva': 'Ival',
                           'pip': 'Hpr',
                           'sal': 'diOH-Bz',
                           'nrp': 'X',
                           'dpg': 'dhpg'}  # TODO: verify this one
    return as_replacement_dict.get(as_name, as_name)


def parse_substrate_predictions(nrps_predictions: Dict[str, str]) -> List[Tuple[str, str]]:
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


def filter_norine_as(as_list, be_strict=False) -> List[str]:
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
            filtered_list.append(map_as_names_to_norine(sub_monomer))
    return filtered_list


def get_norine_url_for_specificities(specificities, be_strict=True) -> Optional[str]:
    """ generate NORINE URL string for direct querying from array of specificity predictions
        use be_strict=False to add * after each monomer
    """
    modulelist = []
    if not specificities:
        return None

    for domain_specificity_list in specificities:
        if not domain_specificity_list:
            logging.error("retrieved empty domain list for assembling per protein NORINE link string")
            raise RuntimeError("empty domain list for NORINE generation")
        if len(domain_specificity_list) == 1:
            if be_strict:
                modulelist.append(map_as_names_to_norine(domain_specificity_list[0]))
            else:
                modulelist.append(map_as_names_to_norine(domain_specificity_list[0])+"*")
        else:  # > 1
            if be_strict:
                modules = []
                for element in filter_norine_as(domain_specificity_list, be_strict=True):
                    modules.append(map_as_names_to_norine(element))
                modulelist.append("[" + "|".join(modules) + "]")
            else:
                # we have to use be_strict to remove X from the list of predictions,
                # as otherwise consenus: nrp will always match
                monomers = [map_as_names_to_norine(element)+"*"
                            for element in filter_norine_as(domain_specificity_list, be_strict=True)]
                # Norine doesn't allow "x*" as a "relaxed" monomer, so we have to replace this with "x"
                monomers = list(map(lambda x: "x" if x == "x*" else x, monomers))
                modulelist.append("["+"|".join(monomers)+"]")
    return "http://bioinfo.lifl.fr/norine/fingerPrintSearch.jsp?nrps1=" + ",".join(modulelist)


class NrpspksLayer(ClusterLayer):
    """ A """
    def __init__(self, results, cluster_rec, record):
        self.url_strict = {}  # gene name -> url
        self.url_relaxed = {}  # gene name -> url
        super().__init__(record, cluster_rec)
        self.transatpks = False
        assert isinstance(results, NRPS_PKS_Results), list(results)
        self.results = results

        cluster_number = cluster_rec.get_cluster_number()
        default_prediction = ("N/A", False)
        self.monomer, self.used_domain_docking = results.cluster_predictions.get(cluster_number, default_prediction)

    def get_domain_details(self) -> Optional[Dict[str, Union[str, List[JSONOrf]]]]:
        """ Creates a JSON-like structure for domains, used by javascript in
            drawing the domains
        """
        orfs = []  # type: List[JSONOrf]
        for feature in self.cluster_rec.cds_children:
            if not feature.nrps_pks:
                continue
            js_orf = JSONOrf(feature)
            for domain in feature.nrps_pks.domains:
                js_orf.add_domain(self.parse_domain(domain, feature))
            orfs.append(js_orf)

        if orfs:
            return {'id': "cluster-%s-details" % self.get_cluster_number(),
                    'orfs': orfs}

        return None

    @property
    def warning(self) -> str:
        """ A caveat for structure prediction accuracy """
        # TODO: don't create warning if no structure image provided
        core = ("Rough prediction of core scaffold based on assumed %s;"
                " tailoring reactions not taken into account")
        if self.used_domain_docking:
            detail = "PKS linker matching"
        else:
            detail = "PKS/NRPS colinearity"

        return core % detail

    @property
    def sidepanel_predictions(self) -> Dict[str, List[List[Tuple[str, str]]]]:
        """ The sidepanel predictions for gene domains, includes URLs for
            related searches
        """
        sidepanel_predictions = {}  # type: Dict[str, List[List[Tuple[str, str]]]]
        features = self.cluster_rec.cds_children
        for feature in features:
            if not feature.nrps_pks:
                continue

            per_cds_predictions = []
            gene_id = feature.get_name()

            for domain in feature.nrps_pks.domains:
                per_a_domain_predictions = []
                preds = parse_substrate_predictions(domain.predictions)
                if not preds:
                    continue

                if domain.name == "AMP-binding":
                    for key, val in preds:
                        if key in ["SNN score", "%ID to nearest neighbour"]:
                            continue
                        if key == "PrediCAT":
                            val = val.rsplit("-")[-1]
                        values = filter_norine_as(val.split(","))
                        if values:
                            per_a_domain_predictions.extend(val.split(","))
                    unique_a_domain_predictions = list(set(per_a_domain_predictions))
                    per_cds_predictions.append(unique_a_domain_predictions)

                if gene_id not in sidepanel_predictions:
                    sidepanel_predictions[gene_id] = []
                sidepanel_predictions[gene_id].append(preds)

            if per_cds_predictions:
                self.url_strict[gene_id] = get_norine_url_for_specificities(per_cds_predictions)
                if self.url_strict[gene_id]:
                    self.url_relaxed[gene_id] = get_norine_url_for_specificities(per_cds_predictions,
                                                                 be_strict=False)

        return sidepanel_predictions

    def is_nrps(self) -> bool:
        """ is the cluster a NRPS or NRPS hybrid """
        return 'nrps' in self.cluster_rec.products

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
                monomers = [map_as_names_to_norine(element.lower())
                            for element in filter_norine_as(monomers, be_strict=True)]
            else:
                monomers = [map_as_names_to_norine(element.lower()) + "*"
                            for element in filter_norine_as(monomers)]
            # Norine doesn't allow "x*" as a "relaxed" monomer, so we have to replace this with "x"
            monomers = list(map(lambda x: "x" if x == "x*" else x, monomers))

            if monomers:
                nrpslist.append("nrps" + str(i) + "=" + ",".join(monomers))
            i += 1
        urlstring = "http://bioinfo.lifl.fr/norine/fingerPrintSearch.jsp?"+"&".join(nrpslist)
        return urlstring

    def parse_domain(self, domain: NRPSPKSQualifier.Domain, feature: CDSFeature
                     ) -> JSONDomain:
        "Convert a NRPS/PKS domain string to a dict useable by json.dumps"
        predictions = parse_substrate_predictions(domain.predictions)

        # Create url_link to NaPDoS for C and KS domains
        napdoslink = ""
        domainseq = str(feature.translation)[domain.start:domain.end]
        base = ("http://napdos.ucsd.edu/cgi-bin/process_request.cgi?"
                "query_type=aa&amp;ref_seq_file=all_{0}_public_12062011.faa"
                "&amp;Sequence=%3E{0}_domain_from_antiSMASH%0D{1}")
        if domain.name == "PKS_KS":
            napdoslink = base.format("KS", domainseq)
        elif "Condensation" in domain.name:
            napdoslink = base.format("C", domainseq)
        blastlink = ("http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE=Proteins"
                     "&amp;PROGRAM=blastp&amp;BLAST_PROGRAMS=blastp"
                     "&amp;QUERY={}"
                     "&amp;LINK_LOC=protein&amp;PAGE_TYPE=BlastSearch").format(domainseq)

        dna_sequence = feature.extract(self.record.seq_record.seq)
        return JSONDomain(domain, predictions, napdoslink, blastlink, domainseq, dna_sequence)
