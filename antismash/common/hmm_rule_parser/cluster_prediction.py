# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Detects specific domains and defines clusters based on domains detected
"""

from collections import defaultdict
import logging
from typing import Any, Dict, List, Set, Tuple

from Bio.SearchIO._model.hsp import HSP

from antismash.common import fasta, path, serialiser
from antismash.common.secmet import Record, Cluster, CDSFeature, FeatureLocation
from antismash.common.secmet.locations import location_contains_other
from antismash.common.secmet.qualifiers import GeneFunction, SecMetQualifier
from antismash.common.subprocessing import run_hmmsearch
from antismash.common.hmm_rule_parser import rule_parser
from antismash.common.signature import get_signature_profiles


class CDSResults:
    """ Tracks the detection results for a single CDS """
    def __init__(self, cds: CDSFeature, domains: List[SecMetQualifier.Domain],
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

    def annotate(self, tool: str) -> None:
        """ Annotates a CDSFeature with the results gathered """
        all_matching = set()
        if not self.cds.sec_met:
            self.cds.sec_met = SecMetQualifier(self.domains)
        else:
            all_matching.update(set(self.cds.sec_met.domain_ids))
            self.cds.sec_met.add_domains(self.domains)
        for cluster_type, matching_domains in self.definition_domains.items():
            all_matching.update(matching_domains)
            for domain in matching_domains:
                self.cds.gene_functions.add(GeneFunction.CORE, tool, domain, cluster_type)

        # and add all detected domains as ADDITIONAL if not CORE
        for secmet_domain in self.cds.sec_met.domains:
            if secmet_domain.name in all_matching:
                continue
            self.cds.gene_functions.add(GeneFunction.ADDITIONAL, secmet_domain.tool,
                                        secmet_domain.name)

    def to_json(self) -> Dict[str, Any]:
        """ Constructs a JSON representation of a CDSResults instance"""
        json = {"cds_name": self.cds.get_name(),
                "domains": [domain.to_json() for domain in self.domains],
                "definition_domains": {key: list(val) for key, val in self.definition_domains.items()}
                }
        return json

    @staticmethod
    def from_json(json: Dict[str, Any], record: Record) -> "CDSResults":
        """ Constructs a CDSResults instance from a JSON representation """
        domains = []
        for json_domain in json["domains"]:
            domains.append(SecMetQualifier.Domain.from_json(json_domain))

        cds = record.get_cds_by_name(json["cds_name"])
        definition_domains = {key: set(val) for key, val in json["definition_domains"].items()}

        return CDSResults(cds, domains, definition_domains)


class RuleDetectionResults:
    """ A container for the all results of running the cluster prediction """
    def __init__(self, cds_by_cluster: Dict[Cluster, List[CDSResults]],
                 tool: str) -> None:
        self.cds_by_cluster = cds_by_cluster
        self.tool = str(tool)

    @property
    def clusters(self) -> List[Cluster]:
        """ A list of Clusters predicted """
        return list(self.cds_by_cluster)

    def annotate_cds_features(self) -> None:
        """ Annotate relevant CDS features with the HMM information detected """
        for cds_results in self.cds_by_cluster.values():
            for cds_result in cds_results:
                cds_result.annotate(self.tool)

    def to_json(self) -> Dict[str, Any]:
        """ Constructs a JSON representation from the RuleDetectionResults instance """
        cds_results_json = []  # type: List[Tuple[Dict[str, Any], List[Dict[str, Any]]]]
        json = {"tool": self.tool,
                "cds_by_cluster": cds_results_json}

        for cluster, cds_results in self.cds_by_cluster.items():
            json_cluster = serialiser.feature_to_json(cluster.to_biopython()[0])
            json_cds_results = [result.to_json() for result in cds_results]
            cds_results_json.append((json_cluster, json_cds_results))

        return json

    @staticmethod
    def from_json(json: Dict[str, Any], record: Record) -> "RuleDetectionResults":
        """ Constructs a RuleDetectionResults instance from a JSON representation """
        cds_by_cluster = {}
        for json_cluster, json_cds_results in json["cds_by_cluster"]:
            cluster = Cluster.from_biopython(serialiser.feature_from_json(json_cluster))
            cds_results = [CDSResults.from_json(result_json, record) for result_json in json_cds_results]
            cds_by_cluster[cluster] = cds_results

        return RuleDetectionResults(cds_by_cluster, json["tool"])


def remove_redundant_clusters(clusters: List[Cluster],
                              rules_by_name: Dict[str, rule_parser.DetectionRule]
                              ) -> List[Cluster]:
    """ Removes clusters which have superiors covering the same (or larger) region
    """
    clusters_by_rule = defaultdict(list)  # type: Dict[str, List[Cluster]]
    for cluster in clusters:
        clusters_by_rule[cluster.product].append(cluster)

    trimmed_clusters = []
    for cluster in clusters:
        rule_name = cluster.product
        is_redundant = False
        for superior in rules_by_name[rule_name].superiors:
            for other_cluster in clusters_by_rule.get(superior, []):
                if location_contains_other(other_cluster.core_location, cluster.core_location):
                    is_redundant = True
                    break
            if is_redundant:
                break
        if not is_redundant:
            trimmed_clusters.append(cluster)
    return trimmed_clusters


def find_clusters(record: Record, cds_by_cluster_type: Dict[str, Set[str]],
                  rules_by_name: Dict[str, rule_parser.DetectionRule]) -> List[Cluster]:
    """ Detects gene clusters based on the identified core genes """
    clusters = []  # type: List[Cluster]

    cds_feature_by_name = record.get_cds_name_mapping()

    for cluster_type, cds_names in cds_by_cluster_type.items():
        cds_features = sorted([cds_feature_by_name[cds] for cds in cds_names])
        rule = rules_by_name[cluster_type]
        cutoff = rule.cutoff
        core_location = cds_features[0].location
        for cds in cds_features[1:]:
            if cds.overlaps_with(FeatureLocation(core_location.start - cutoff,
                                                 core_location.end + cutoff)):
                core_location = FeatureLocation(min(cds.location.start, core_location.start),
                                                max(cds.location.end, core_location.end))
                assert core_location.start >= 0 and core_location.end <= len(record)
                continue
            # create the previous cluster and start a new location
            surrounds = FeatureLocation(max(0, core_location.start - rule.extent),
                                        min(core_location.end + rule.extent, len(record)))
            surrounding_cdses = record.get_cds_features_within_location(surrounds, with_overlapping=False)
            real_start = min(contained.location.start for contained in surrounding_cdses)
            real_end = max(contained.location.end for contained in surrounding_cdses)
            surrounds = FeatureLocation(real_start, real_end)
            clusters.append(Cluster(core_location, surrounding_location=surrounds,
                                    tool="rule-based-clusters", cutoff=cutoff,
                                    neighbourhood_range=rule.extent, product=cluster_type,
                                    detection_rule=str(rule.conditions)))
            core_location = cds.location

        # finalise the last cluster
        surrounds = FeatureLocation(max(0, core_location.start - rule.extent),
                                    min(core_location.end + rule.extent, len(record)))
        clusters.append(Cluster(core_location, surrounding_location=surrounds,
                                tool="rule-based-clusters", cutoff=cutoff,
                                neighbourhood_range=rule.extent, product=cluster_type,
                                detection_rule=str(rule.conditions)))

    # fit to record if outside
    for cluster in clusters:
        contained = FeatureLocation(max(0, cluster.location.start),
                                    min(cluster.location.end, len(record)))
        if contained != cluster.location:
            cluster.location = contained

    clusters = remove_redundant_clusters(clusters, rules_by_name)

    logging.debug("%d rule-based cluster(s) found in record", len(clusters))
    return clusters


def hsp_overlap_size(first: HSP, second: HSP) -> int:
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


def create_rules(rule_file: str, signature_names: Set[str]
                 ) -> List[rule_parser.DetectionRule]:
    """ Creates DetectionRule instances from the default rules file

        Args:
            rule_file: A path to a file containing cluster rules to use.

        Returns:
            A list of DetectionRules.
    """
    rules = []
    with open(rule_file, "r") as ruledata:
        parser = rule_parser.Parser("".join(ruledata.readlines()), signature_names)
    for rule in parser.rules:
        rules.append(rule)
    return rules


def apply_cluster_rules(record: Record, results_by_id: Dict[str, List[HSP]],
                        rules: List[rule_parser.DetectionRule]
                        ) -> Tuple[Dict[str, Dict[str, Set[str]]],
                                   Dict[str, Set[str]]]:
    """
        Run detection rules over each CDS and classify them if relevant.
        A CDS can satisfy multiple rules. If so, all rules satisfied
        will form part of the type string, separated by '-'.

        The 'other' type has a lower precedence than other rules and a hit with
        the 'other' rule will be ignored if another rule is also satisfied.

        Args:
            record: the record being checked
            results_by_id: A dict of CDS ID to a list of HSP results
            rules: A list of DetectionRule instances

        Returns:
            A tuple of
                a dictionary mapping CDS ID to
                    a dictionary mapping cluster type string to
                        a set of domains used to determine the cluster
                and a dictionary mapping rule name to
                    a set of CDS feature names that matched the rule
    """
    if not results_by_id:
        return {}, {}

    cds_with_hits = sorted(results_by_id, key=lambda gene_id: record.get_cds_by_name(gene_id).location.start)

    cds_domains_by_cluster_type = {}
    cluster_type_hits = defaultdict(set)  # type: Dict[str, Set[str]]
    for cds_name in cds_with_hits:
        feature = record.get_cds_by_name(cds_name)
        feature_start, feature_end = sorted([feature.location.start, feature.location.end])
        results = []  # type: List[str]
        rule_texts = []
        info_by_range = {}  # type: Dict[int, Tuple[Dict[str, CDSFeature], Dict[str, List[HSP]]]]
        domain_matches = set()  # type: Set[str]
        domains_by_cluster = {}  # type: Dict[str, Set[str]]
        for rule in rules:
            if rule.cutoff not in info_by_range:
                location = FeatureLocation(feature_start - rule.cutoff, feature_end + rule.cutoff)
                nearby = record.get_cds_features_within_location(location, with_overlapping=True)
                nearby_features = {neighbour.get_name(): neighbour for neighbour in nearby}
                nearby_results = {neighbour: results_by_id[neighbour]
                                  for neighbour in nearby_features if neighbour in results_by_id}
                info_by_range[rule.cutoff] = (nearby_features, nearby_results)
            nearby_features, nearby_results = info_by_range[rule.cutoff]
            matching = rule.detect(cds_name, nearby_features, nearby_results)
            if matching.met and matching.matches:
                domains_by_cluster[rule.name] = matching.matches
                results.append(rule.name)
                rule_texts.append(rule.reconstruct_rule_text())
                domain_matches.update(matching.matches)
                cluster_type_hits[rule.name].add(cds_name)
        if domains_by_cluster:
            cds_domains_by_cluster_type[cds_name] = domains_by_cluster
    return cds_domains_by_cluster_type, cluster_type_hits


def detect_clusters_and_signatures(record: Record, signature_file: str, seeds_file: str,
                                   rules_file: str, filter_file: str, tool: str) -> RuleDetectionResults:
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
    """
    full_fasta = fasta.get_fasta_from_record(record)
    # if there's no CDS features, don't try to do anything
    if not full_fasta:
        return RuleDetectionResults({}, tool)
    sig_by_name = {sig.name: sig for sig in get_signature_profiles(signature_file)}
    rules = create_rules(rules_file, set(sig_by_name))
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
    cds_domains_by_cluster, cluster_type_hits = apply_cluster_rules(record, results_by_id, rules)

    # Find number of sequences on which each pHMM is based
    num_seeds_per_hmm = get_sequence_counts(signature_file)

    # Save final results to record
    rules_by_name = {rule.name: rule for rule in rules}
    clusters = find_clusters(record, cluster_type_hits, rules_by_name)
    strip_inferior_domains(cds_domains_by_cluster, rules_by_name)

    cds_results_by_cluster = {}
    for cluster in clusters:
        cds_results = []
        for cds in record.get_cds_features_within_location(cluster.location):
            domains = []
            for hsp in results_by_id.get(cds.get_name(), []):
                domains.append(SecMetQualifier.Domain(hsp.query_id, hsp.evalue, hsp.bitscore,
                                                      num_seeds_per_hmm[hsp.query_id], tool))
            if domains:
                cds_results.append(CDSResults(cds, domains, cds_domains_by_cluster.get(cds.get_name(), {})))
        cds_results_by_cluster[cluster] = cds_results

    return RuleDetectionResults(cds_results_by_cluster, tool)


def strip_inferior_domains(cds_domains_by_cluster: Dict[str, Dict[str, Set[str]]],
                           rules_by_name: Dict[str, rule_parser.DetectionRule]) -> None:
    """ Remove any domain hits for each inferior rule within a CDS that the CDS also
        satisfies the rule's superior.

        Modifies cds_domains_by_cluster in place.
    """
    for domains_by_cluster in cds_domains_by_cluster.values():
        all_satisfied = set(domains_by_cluster)
        for product in all_satisfied:
            rule = rules_by_name[product]
            if set(rule.superiors).intersection(all_satisfied):
                domains_by_cluster.pop(product)


def get_sequence_counts(details_file: str) -> Dict[str, int]:
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
                result[hmm.name] = int(line[6:].strip())
                break
        if hmm.name not in result:
            raise ValueError("Unknown number of seeds for hmm file: %s" % details_file)

    return result
