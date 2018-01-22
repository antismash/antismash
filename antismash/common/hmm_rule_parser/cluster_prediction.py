# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Detects specific domains and defines clusters based on domains detected
"""

from collections import defaultdict
import logging
from typing import Dict, List, Set, Tuple

from Bio.SearchIO._model.hsp import HSP

from antismash.common import fasta, path
from antismash.common.secmet.feature import ClusterBorder, SecMetQualifier, CDSFeature, \
                                            GeneFunction, FeatureLocation
from antismash.common.subprocessing import run_hmmsearch
from antismash.common.hmm_rule_parser import rule_parser
from antismash.common.signature import get_signature_profiles


class Domain:
    """ A simple container for the information needed to create a domain """
    def __init__(self, res: HSP, nseeds: str, tool: str):
        self.query_id = str(res.query_id)
        self.evalue = float(res.evalue)
        self.bitscore = float(res.bitscore)
        self.nseeds = str(nseeds)
        self.tool = tool

    def __repr__(self):
        return str(self)

    def __str__(self):
        ret = "{} E-value: {}, bitscore: {}, seeds: {}"
        return ret.format(self.query_id, self.evalue, self.bitscore, self.nseeds)


class CDSResults:
    def __init__(self, cds: CDSFeature, domains: List[Domain],
                 definition_domains: Dict[str, Set[str]]) -> None:
        """ Arguments:
                cds: the CDSFeature these results were based on
                domains: a list of Domains that were found in the CDSFeature
                definition_domains: a dictionary mapping cluster type to
                        domain names used to define that cluster, from this CDS
        """
        self.cds = cds
        self.domains = domains
        assert domains, "CDSResults not possible without some domains"
        assert isinstance(definition_domains, dict), type(definition_domains)
        # empty definition domains is ok
        self.definition_domains = definition_domains

    def annotate(self, products: Set[str]) -> None:
        if not self.cds.sec_met:
            self.cds.sec_met = SecMetQualifier(set(self.definition_domains), self.domains)
        else:
            self.cds.sec_met.add_products(set(products))
            self.cds.sec_met.add_domains(self.domains)
        all_matching = set()
        for cluster_type, matching_domains in self.definition_domains.items():
            all_matching.update(matching_domains)
            for domain in matching_domains:
                # TODO: use tool name instead of 'cluster_definition'
                self.cds.gene_functions.add(GeneFunction.CORE, "cluster_definition",
                                            "%s: %s" % (cluster_type, domain))

        # and add all detected domains as ADDITIONAL if not CORE
        for secmet_domain in self.cds.sec_met.domains:
            if secmet_domain.query_id in all_matching:
                continue
            self.cds.gene_functions.add(GeneFunction.ADDITIONAL, secmet_domain.tool,
                                        secmet_domain)


class DetectionResults:
    def __init__(self, cds_by_cluster: Dict[ClusterBorder, List[CDSResults]]) -> None:
        self.cds_by_cluster = cds_by_cluster

    @property
    def borders(self) -> List[ClusterBorder]:
        return list(self.cds_by_cluster)

    def annotate_cds_features(self) -> None:
        for border, cds_results in self.cds_by_cluster.items():
            for cds_result in cds_results:
                cds_result.annotate(border.products)


def find_clusters(record, cds_by_cluster_type, rules_by_name) -> List[ClusterBorder]:
    """ Detects gene clusters based on the identified core genes """
    clusters = []  # type: List[ClusterBorder]

    cds_feature_by_name = record.get_cds_name_mapping()

    for cluster_type, cds_names in cds_by_cluster_type.items():
        cds_features = sorted([cds_feature_by_name[cds] for cds in cds_names])
        rule = rules_by_name[cluster_type]
        cutoff = rule.cutoff
        extent = rule.extent
        start, end = sorted([cds_features[0].location.start, cds_features[0].location.end])
        cluster = ClusterBorder(FeatureLocation(start, end), tool="rule-based-clusters",
                                cutoff=cutoff, extent=extent, products=[cluster_type])
        assert cds_features[0].is_contained_by(cluster)
        assert cds_features[0] in record.get_cds_features_within_location(cluster.location)
        clusters.append(cluster)
        for cds in cds_features[1:]:
            feature_start, feature_end = sorted([cds.location.start, cds.location.end])
            dummy_location = FeatureLocation(cluster.location.start - cutoff, cluster.location.end + cutoff)
            if cds.is_contained_by(dummy_location):
                start = min(feature_start, start)
                end = max(feature_end, end)
                cluster.location = FeatureLocation(start, end)
            else:
                start = feature_start
                end = feature_end
                cluster = ClusterBorder(FeatureLocation(start, end), tool="rule-based-clusters",
                                        cutoff=cutoff, extent=extent, products=[cluster_type])
                clusters.append(cluster)

    for cluster in clusters:
        cluster.rules = [str(rules_by_name[product].conditions) for product in cluster.products]
        if cluster.location.start < 0:
            cluster.location = FeatureLocation(0, cluster.location.end)
            cluster.contig_edge = True
        if cluster.location.end > len(record):
            cluster.location = FeatureLocation(cluster.location.start, len(record))
            cluster.contig_edge = True

    logging.debug("%d rule-based cluster(s) found in record", len(clusters))
    return clusters


def hsp_overlap_size(first, second) -> int:
    """ Find the size of an overlapping region of two HSPs.

        Args:
            first: a HSP instance
            second: a HSP instance

        Returns:
            The size of the overlap in bases or zero if there is no overlap
    """
    assert first.hit_start < first.hit_end
    assert second.hit_start < second.hit_end
    segment_start = max(first.hit_start, second.hit_start)
    segment_end = min(first.hit_end, second.hit_end)
    return max(0, segment_end - segment_start)


def filter_results(results: List[HSP], results_by_id: Dict[str, List[HSP]], filter_file: str,
                   signature_names: Set[str]) -> Tuple[List[HSP], Dict[str, List[HSP]]]:
    """ Filter results by comparing scores of different models """
    # TODO: filterfile = path.get_full_path(__file__, "filterhmmdetails.txt")
    for line in open(filter_file, "r"):
        line = line.strip()
        equivalence_group = set(line.split(","))
        unknown = equivalence_group - signature_names
        if unknown:
            raise ValueError("Equivalence group contains unknown identifiers: %s" % (unknown))
        removed_ids = set()  # type: Set[int]
        for cds, cdsresults in results_by_id.items():
            # Check if multiple competing HMM hits are present
            hits = set(hit.query_id for hit in cdsresults)
            if len(hits & equivalence_group) < 2:
                continue
            # Identify overlapping hits
            overlapping_groups = []  # type: List[Set[HSP]]
            for hit in cdsresults:
                for otherhit in cdsresults:
                    if hit == otherhit or hsp_overlap_size(hit, otherhit) <= 20:
                        continue
                    new_group_needed = True
                    pairing = {hit, otherhit}
                    for group in overlapping_groups:
                        if pairing & group:
                            group.update(pairing)
                            new_group_needed = False
                    if new_group_needed:
                        overlapping_groups.append(pairing)

            # find the best in each group
            for group in overlapping_groups:
                # start with one of them
                best = list(group)[0]
                for hit in group:
                    if not best or hit.bitscore > best.bitscore:
                        best = hit
                # remove the rest
                for hit in group:
                    if hit != best:
                        if id(hit) not in removed_ids:
                            del results[results.index(hit)]
                            del results_by_id[cds][results_by_id[cds].index(hit)]
                            removed_ids.add(id(hit))
            assert results_by_id[cds]  # should always have one remaining
    return results, results_by_id


def filter_result_multiple(results: List[HSP], results_by_id: Dict[str, HSP]) -> Tuple[List[HSP], Dict[str, HSP]]:
    """ Filter multiple results of the same model within a gene """
    for cds, hits in results_by_id.items():
        query_scores = {}  # type: Dict[str, Tuple[int, HSP, float]]
        for i, hit in enumerate(hits):
            if query_scores.get(hit.query_id, (0, 0, -1))[2] < hit.bitscore:
                query_scores[hit.query_id] = (i, hit, hit.bitscore)
        best_hits = set(info[:2] for info in query_scores.values())
        results_by_id[cds] = [i[1] for i in sorted(best_hits)]
    results.clear()
    for cds in results_by_id:
        results.extend(results_by_id[cds])
    results.sort(key=lambda hit: hit.hit_start)
    return results, results_by_id


def create_rules(enabled_cluster_types: List[str], rule_file: str, signature_names: Set[str]
                 ) -> List[rule_parser.DetectionRule]:
    """ Creates DetectionRule instances from the default rules file

        Args:
            enabled_cluster_types: A list of type names.
                    Only types in this list will have a DetectionRule.
            rule_file: A path to a file containing cluster rules to use.

        Returns:
            A list of DetectionRules.
    """
    rules = []
    with open(rule_file, "r") as ruledata:
        parser = rule_parser.Parser("".join(ruledata.readlines()), signature_names)
    for rule in parser.rules:
        if rule.name in enabled_cluster_types:
            rules.append(rule)
    return rules


def apply_cluster_rules(record, results_by_id: Dict[str, HSP], feature_by_id: Dict[str, CDSFeature],
                        rules: List[rule_parser.DetectionRule]) -> Dict[str, Dict[str, Set[str]]]:
    """
        Run detection rules over each CDS and classify them if relevant.
        A CDS can satisfy multiple rules. If so, all rules satisfied
        will form part of the type string, separated by '-'.

        The 'other' type has a lower precedence than other rules and a hit with
        the 'other' rule will be ignored if another rule is also satisfied.

        Args:
            record: the record being checked
            results_by_id: A dict of CDS ID to a list of HSP results
            feature_by_id: A dict of CDS ID to CDS Feature
            rules: A list of DetectionRule instances

        Returns:
            A dictionary mapping CDS ID to
                a dictionary mapping cluster type string to
                    a set of domains used to determine the cluster.
    """
    if not results_by_id:  # TODO: update docstring
        return {}, {}

    cds_with_hits = sorted(results_by_id, key=lambda gene_id: feature_by_id[gene_id].location.start)

    def calculate_distance(first, second):
        """ Calculate the distance between two FeatureLocations """
        first_start, first_end = sorted([first.start, first.end])
        second_start, second_end = sorted([second.start, second.end])
        return min(abs(first_end - second_start), abs(second_end - first_start),
                   abs(first_start - second_start), abs(first_end - second_end))

    cds_domains_by_cluster_type = {}
    cluster_type_hits = defaultdict(set)
    for cds in cds_with_hits:
        feature = feature_by_id[cds]
        feature_start, feature_end = sorted([feature.location.start, feature.location.end])
        results = []  # type: List[str]
        rule_texts = []
        info_by_range = {}  # type: Dict[int, Tuple[Dict[str, CDSFeature], Dict[str, List[HSP]]]]
        domain_matches = set()  # type: Set[str]
        domains_by_cluster = {}  # type: Dict[str, Set[str]]
        for rule in rules:
            if rule.cutoff not in info_by_range:
                # TODO: improve performance
                nearby = record.get_cds_features_within_location(FeatureLocation(feature_start - rule.cutoff, feature_end + rule.cutoff), with_overlapping=True)
                nearby_features = {neighbour.get_name(): neighbour for neighbour in nearby}
                nearby_results = {neighbour: results_by_id[neighbour] for neighbour in nearby_features if neighbour in results_by_id}
                info_by_range[rule.cutoff] = (nearby_features, nearby_results)
            nearby_features, nearby_results = info_by_range[rule.cutoff]
            matching = rule.detect(cds, nearby_features, nearby_results)
            if matching.met and matching.matches:
                domains_by_cluster[rule.name] = matching.matches
                results.append(rule.name)
                rule_texts.append(rule.reconstruct_rule_text())
                domain_matches.update(matching.matches)
                cluster_type_hits[rule.name].add(cds)
        if domains_by_cluster:
            cds_domains_by_cluster_type[cds] = domains_by_cluster
    return cds_domains_by_cluster_type, cluster_type_hits


def detect_borders_and_signatures(record, signature_file: str, seeds_file: str,
                                  rules_file: str, filter_file: str, tool: str,
                                  options) -> None:
    """ Compares all CDS features in a record with HMM signatures and generates
        Cluster features based on those hits and the current cluster detection
        rules.

        Arguments:
            record: the record to analyse
            signature_file: a tab separated file; each row being a single HMM reference
                        with columns: label, description, minimum score cutoff, hmm path
            seeds_file: the file containing all HMM profiles
            rules_file: the file containing all the rules to use for cluster definition
            filter_file: a file containing equivalence sets of HMMs
            tool: the name of the tool providing the HMMs (e.g. clusterfinder, rule_based_clusters)
            options: antismash Config
    """
    enabled_cluster_types = options.enabled_cluster_types
    feature_by_id = record.get_cds_name_mapping()
    # if there's no CDS features, don't try to do anything
    if not feature_by_id:
        return None
    full_fasta = fasta.get_fasta_from_record(record)
    sig_by_name = {sig.name: sig for sig in get_signature_profiles(signature_file)}
    rules = create_rules(enabled_cluster_types, rules_file, set(sig_by_name))
    results = []
    results_by_id = {}  # type: Dict[str, HSP]

    runresults = run_hmmsearch(seeds_file, full_fasta, use_tempfile=True)
    for runresult in runresults:
        acc = runresult.accession.split('.')[0]
        # Store result if it is above cut-off
        for hsp in runresult.hsps:
            if hsp.query_id in sig_by_name:
                sig = sig_by_name[hsp.query_id]
            elif acc in sig_by_name:
                sig = sig_by_name[acc]
            else:
                raise ValueError('Failed to find signature for ID %s / ACC %s' % (
                                                    hsp.query_id, acc))
            if hsp.bitscore > sig.cutoff:
                results.append(hsp)
                if hsp.hit_id not in results_by_id:
                    results_by_id[hsp.hit_id] = [hsp]
                else:
                    results_by_id[hsp.hit_id].append(hsp)

    # Filter results by comparing scores of different models (for PKS systems)
    results, results_by_id = filter_results(results, results_by_id, filter_file, set(sig_by_name))

    # Filter multiple results of the same model in one gene
    results, results_by_id = filter_result_multiple(results, results_by_id)

    # Use rules to determine gene clusters
    cds_domains_by_cluster, cluster_type_hits = apply_cluster_rules(record, results_by_id, feature_by_id, rules)

    # Find number of sequences on which each pHMM is based
    num_seeds_per_hmm = get_sequence_counts(signature_file)

    # Save final results to record
    rules_by_name = {rule.name: rule for rule in rules}
    clusters = find_clusters(record, cluster_type_hits, rules_by_name)

    cds_results_by_cluster = {}
    for cluster in clusters:
        record.add_cluster_border(cluster)
        cds_results = []
        cluster_extent = FeatureLocation(cluster.location.start - cluster.extent,
                                         cluster.location.end + cluster.extent)
        for cds in record.get_cds_features_within_location(cluster_extent):
            domains = [Domain(res, num_seeds_per_hmm[res.query_id], tool) for res in results_by_id.get(cds.get_name(), [])]
            if domains:
                cds_results.append(CDSResults(cds, domains, cds_domains_by_cluster.get(cds.get_name(), {})))
        cds_results_by_cluster[cluster] = cds_results

    results = DetectionResults(cds_results_by_cluster)
    return results


def get_sequence_counts(details_file: str) -> Dict[str, str]:
    """ Gets the number of sequences/seeds used to generate each HMM signature

        Arguments:
            detail_file: a file containing all HMMs

        Returns:
            a dictionary mapping HMM name to the number of sequences used to
                generate it
    """
    result = {}
    for hmm in get_signature_profiles(details_file):
        for line in open(path.get_full_path(details_file, hmm.hmm_file), 'r'):
            if line.startswith('NSEQ '):
                result[hmm.name] = line[6:].strip()
                break
        # TODO: ideally this shouldn't ever happen, clean up inputs and change to error
        if hmm.name not in result:
            result[hmm.name] = "?"

    return result
