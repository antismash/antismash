# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Generates HTML and JSON for the nrps_pks module """

import logging
import re
from typing import Any, Dict, Iterator, Iterable, List, Optional, Tuple, Union

from jinja2 import FileSystemLoader, Environment, StrictUndefined

from antismash.common import path
from antismash.common.layers import ClusterLayer, RecordLayer, OptionsLayer
from antismash.common.secmet import CDSFeature, Record, Cluster
from antismash.common.secmet.qualifiers import NRPSPKSQualifier
from antismash.config import ConfigType

from .results import NRPS_PKS_Results, UNKNOWN


class JSONBase(dict):
    """ A base class for JSON-serialisable objects """
    def __init__(self, keys: List[str]) -> None:
        super().__init__()
        self._keys = keys

    def __getitem__(self, key: str) -> Any:
        return getattr(self, key)

    def items(self) -> Iterator[Tuple[str, Any]]:  # type: ignore
        for key in self._keys:
            yield (key, getattr(self, key))

    def values(self) -> Iterator[Any]:  # type: ignore
        for key in self._keys:
            yield getattr(self, key)

    def __len__(self) -> int:
        return len(self._keys)


class JSONDomain(JSONBase):
    """ A JSON-serialisable object for simplifying domain datatypes throughout this file """
    def __init__(self, domain: NRPSPKSQualifier.Domain, predictions: List[Tuple[str, str]], napdos_link: str,
                 blast_link: str, sequence: str, dna: str) -> None:
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
    def __init__(self, feature: CDSFeature) -> None:
        super().__init__(['id', 'sequence', 'domains'])
        self.sequence = feature.translation
        self.id = feature.get_name()
        self.domains = []  # type: List[JSONDomain]

    def add_domain(self, domain: JSONDomain) -> None:
        """ Add a JSONDomain to the list of domains in this ORF """
        assert isinstance(domain, JSONDomain)
        self.domains.append(domain)


def will_handle(products: List[str]) -> bool:
    """ Returns true if one or more relevant products are present """
    return bool(set(products).intersection({"nrps", "t1pks", "t2pks", "transatpks",
                                            "nrpsfragment", "otherks"}))


def generate_js_domains(cluster_json: Dict[str, Any], record: Record, results: NRPS_PKS_Results,
                        options: ConfigType) -> Optional[Dict[str, Union[str, List[JSONOrf]]]]:
    """ Generate """
    assert isinstance(results, NRPS_PKS_Results), type(results)
    cluster_feature = record.get_cluster(cluster_json['idx'])
    cluster = NrpspksLayer(results, cluster_feature, RecordLayer(record, results, OptionsLayer(options)))
    return cluster.get_domain_details()


def generate_details_div(cluster_layer: ClusterLayer, results: NRPS_PKS_Results,
                         record_layer: RecordLayer, options_layer: OptionsLayer) -> str:
    """ Generate the main HTML with results from the NRPS/PKS module """
    env = Environment(loader=FileSystemLoader(path.get_full_path(__file__, 'templates')),
                      autoescape=True, undefined=StrictUndefined)
    template = env.get_template('details.html')
    cluster = NrpspksLayer(results, cluster_layer.cluster_feature, record_layer)
    details_div = template.render(record=record_layer,
                                  cluster=cluster,
                                  options=options_layer)
    return details_div


def generate_sidepanel(cluster_layer: ClusterLayer, results: NRPS_PKS_Results,
                       record_layer: RecordLayer, options_layer: OptionsLayer) -> str:
    """ Generate the sidepanel HTML with results from the NRPS/PKS module """
    env = Environment(loader=FileSystemLoader(path.get_full_path(__file__, 'templates')),
                      autoescape=True, undefined=StrictUndefined)
    template = env.get_template('sidepanel.html')
    nrps_layer = NrpspksLayer(results, cluster_layer.cluster_feature, record_layer)
    sidepanel = template.render(record=record_layer,
                                cluster=nrps_layer,
                                results=results,
                                options=options_layer)
    return sidepanel


def map_as_name_to_norine(as_name: str) -> str:
    """ Maps antiSMASH amino acid nomenclature to NORINE """

    as_replacement_dict = {'bht': 'bOH-Tyr',
                           'dhb': 'diOH-Bz',
                           'iva': 'Ival',
                           'pip': 'Hpr',
                           'sal': 'diOH-Bz',
                           'nrp': 'X',
                           'dpg': 'dhpg'}  # TODO: verify this one
    return as_replacement_dict.get(as_name, as_name)


def filter_norine_as(monomers: List[str], be_strict: bool = False) -> List[str]:
    """ Remove PKS and unknown substrate predictions
        use be_strict = True to also filter nrp/X
    """
    filtered_list = []
    bad_monomers = {'pk', UNKNOWN, 'hydrophilic', 'hydrophobic', 'mal', 'mmal'}
    if be_strict:
        bad_monomers = bad_monomers.union({'nrp', 'X'})
    for monomer in monomers:
        assert '|' not in monomer, monomer
        filtered_list.append(map_as_name_to_norine(monomer))
    return filtered_list


def get_norine_url_for_specificities(specificities: List[List[str]],
                                     be_strict: bool = True) -> Optional[str]:
    """ generate NORINE URL string for direct querying from array of specificity predictions
        use be_strict=False to add * after each monomer
    """
    modules = []
    if not specificities:
        return None

    for domain_specificity_list in specificities:
        assert isinstance(domain_specificity_list, list), type(domain_specificity_list)
        if not domain_specificity_list:
            logging.error("retrieved empty domain list for assembling per protein NORINE link string")
            raise RuntimeError("empty domain list for NORINE generation")
        wildcard = "*"
        separator = "|"
        if be_strict:
            wildcard = ""
        # always be strict to remove X and nrp
        chunks = [monomer for monomer in filter_norine_as(domain_specificity_list, be_strict=True)]
        query = (wildcard + separator).join(chunks) + wildcard
        if len(domain_specificity_list) > 1:
            query = "[" + query + "]"
        # lastly, remove the Norine-forbidden x*
        if be_strict:
            query.replace("x*", "x")
        modules.append(query)
    return "http://bioinfo.lifl.fr/norine/fingerPrintSearch.jsp?nrps1=" + ",".join(modules)


class NrpspksLayer(ClusterLayer):
    """ A wrapper for ClusterLayer that adds some specific sections for NRPS/PKS
        domains and structures.
    """
    # TODO: some of these need to shift to nrps_pks_domains output, the remainder
    #       can move into the NRPS_PKS_Results object instead.
    def __init__(self, results: NRPS_PKS_Results, cluster_feature: Cluster, record: RecordLayer) -> None:
        self.url_strict = {}  # type: Dict[str, str]  # gene name -> url
        self.url_relaxed = {}  # type: Dict[str, str]  # gene name -> url
        self._build_urls(cluster_feature.cds_children)
        super().__init__(record, cluster_feature)
        self.transatpks = False
        assert isinstance(results, NRPS_PKS_Results), type(results)
        self.results = results

        cluster_number = cluster_feature.get_cluster_number()
        default_prediction = ("N/A", False)
        self.monomer, self.used_domain_docking = results.cluster_predictions.get(cluster_number, default_prediction)

    def get_domain_details(self) -> Optional[Dict[str, Union[str, List[JSONOrf]]]]:
        """ Creates a JSON-like structure for domains, used by javascript in
            drawing the domains
        """
        orfs = []  # type: List[JSONOrf]
        for feature in self.cluster_feature.cds_children:
            if not feature.nrps_pks:
                continue
            js_orf = JSONOrf(feature)
            for domain in feature.nrps_pks.domains:
                js_orf.add_domain(self.parse_domain(domain, feature))
            orfs.append(js_orf)

        if orfs:
            return {'id': self.anchor_id,
                    'orfs': orfs}

        return None

    @property
    def warning(self) -> str:
        """ A caveat for structure prediction accuracy """
        if not self.monomer:
            return ""
        core = ("Rough prediction of core scaffold based on assumed %s;"
                " tailoring reactions not taken into account")
        if self.used_domain_docking:
            detail = "PKS linker matching"
        else:
            detail = "PKS/NRPS colinearity"

        return core % detail

    @property
    def sidepanel_features(self) -> List[str]:
        """ Returns a list of relevant CDSFeature names """
        return sorted(feature.get_name() for feature in self.cluster_feature.cds_children if feature.nrps_pks)

    def _build_urls(self, cds_features: Iterable[CDSFeature]) -> None:
        for feature in cds_features:
            if not feature.nrps_pks:
                continue

            feature_name = feature.get_name()

            per_cds_predictions = []

            for domain in feature.nrps_pks.domains:
                if domain.name not in ["AMP-binding", "A-OX"]:
                    continue
                per_a_domain_predictions = set()
                for possibilities in domain.predictions.values():
                    for possibility in filter_norine_as(possibilities.split(","), be_strict=False):
                        per_a_domain_predictions.add(map_as_name_to_norine(possibility))
                per_cds_predictions.append(list(per_a_domain_predictions))

            if not per_cds_predictions:
                continue
            self.url_strict[feature_name] = get_norine_url_for_specificities(per_cds_predictions)
            if self.url_strict[feature_name]:
                self.url_relaxed[feature_name] = get_norine_url_for_specificities(per_cds_predictions,
                                                                                  be_strict=False)

    def is_nrps(self) -> bool:
        """ is the cluster a NRPS or NRPS hybrid """
        return 'nrps' in self.cluster_feature.products

    def get_monomer_prediction(self) -> str:
        "Get the monomer prediction of the cluster"
        return self.monomer

    def get_norine_url_for_cluster(self, be_strict: bool = True) -> str:
        """ Get a NORINE URL string for direct querying
            use be_strict=False to add * after each monomer"""

        monomer_string = self.get_monomer_prediction()
        monomers_per_protein_list = re.findall("\\(.*?\\)", monomer_string)
        i = 1
        nrpslist = []
        for monomers_per_protein in monomers_per_protein_list:
            monomers = monomers_per_protein[1:-1].split("-")

            if be_strict:
                monomers = [map_as_name_to_norine(element.lower())
                            for element in filter_norine_as(monomers, be_strict=True)]
            else:
                monomers = [map_as_name_to_norine(element.lower()) + "*"
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
        predictions = list(domain.predictions.items())

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
