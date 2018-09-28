# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Detects specific domains and defines clusters based on domains detected
"""

from collections import defaultdict
import logging
from typing import Any, Dict, List, Set, Tuple

from Bio.SearchIO._model.hsp import HSP
from Bio.SeqFeature import CompoundLocation

from antismash.common import fasta, path, serialiser
from antismash.common.secmet import Record, ClusterBorder, CDSFeature, Feature, FeatureLocation
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

    def annotate(self, product: str, tool: str) -> None:
        """ Annotates a CDSFeature with the results gathered """
        all_matching = set()
        if not self.cds.sec_met:
            self.cds.sec_met = SecMetQualifier(set(self.definition_domains), self.domains)
        else:
            all_matching.update(set(self.cds.sec_met.domain_ids))
            self.cds.sec_met.add_products({product})
            self.cds.sec_met.add_domains(self.domains)
        for cluster_type, matching_domains in self.definition_domains.items():
            all_matching.update(matching_domains)
            for domain in matching_domains:
                self.cds.gene_functions.add(GeneFunction.CORE, tool,
                                            "%s: %s" % (cluster_type, domain))

        # and add all detected domains as ADDITIONAL if not CORE
        for secmet_domain in self.cds.sec_met.domains:
            if secmet_domain.query_id in all_matching:
                continue
            self.cds.gene_functions.add(GeneFunction.ADDITIONAL, secmet_domain.tool,
                                        str(secmet_domain))

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
    def __init__(self, cds_by_cluster: Dict[ClusterBorder, List[CDSResults]],
                 tool: str) -> None:
        self.cds_by_cluster = cds_by_cluster
        self.tool = str(tool)

    @property
    def borders(self) -> List[ClusterBorder]:
        """ A list of ClusterBorders predicted """
        return list(self.cds_by_cluster)

    def annotate_cds_features(self) -> None:
        """ Annotate relevant CDS features with the HMM information detected """
        for border, cds_results in self.cds_by_cluster.items():
            for cds_result in cds_results:
                cds_result.annotate(border.product, self.tool)

    def to_json(self) -> Dict[str, Any]:
        """ Constructs a JSON representation from the RuleDetectionResults instance """
        cds_results_json = []  # type: List[Tuple[Dict[str, Any], List[Dict[str, Any]]]]
        json = {"tool": self.tool,
                "cds_by_cluster": cds_results_json}

        for border, cds_results in self.cds_by_cluster.items():
            json_border = serialiser.feature_to_json(border.to_biopython()[0])
            json_cds_results = [result.to_json() for result in cds_results]
            cds_results_json.append((json_border, json_cds_results))

        return json

    @staticmethod
    def from_json(json: Dict[str, Any], record: Record) -> "RuleDetectionResults":
        """ Constructs a RuleDetectionResults instance from a JSON representation """
        cds_by_cluster = {}
        for json_border, json_cds_results in json["cds_by_cluster"]:
            border = ClusterBorder.from_biopython(serialiser.feature_from_json(json_border))
            cds_results = [CDSResults.from_json(result_json, record) for result_json in json_cds_results]
            cds_by_cluster[border] = cds_results

        return RuleDetectionResults(cds_by_cluster, json["tool"])


def remove_redundant_borders(borders: List[ClusterBorder],
                             rules_by_name: Dict[str, rule_parser.DetectionRule]
                             ) -> List[ClusterBorder]:
    """ Removes clusters which have superiors covering the same (or larger) region
    """
    borders_by_rule = defaultdict(list)  # type: Dict[str, List[ClusterBorder]]
    for border in borders:
        borders_by_rule[border.product].append(border)

    trimmed_borders = []
    for border in borders:
        rule_name = border.product
        is_redundant = False
        for superior in rules_by_name[rule_name].superiors:
            for other_border in borders_by_rule.get(superior, []):
                if border.is_contained_by(other_border):
                    is_redundant = True
                    break
            if is_redundant:
                break
        if not is_redundant:
            trimmed_borders.append(border)
    return trimmed_borders


def find_clusters(record: Record, cds_by_cluster_type: Dict[str, Set[str]],
                  rules_by_name: Dict[str, rule_parser.DetectionRule]) -> List[ClusterBorder]:
    """ Detects gene borders based on the identified core genes """
    borders = []  # type: List[ClusterBorder]

    cds_feature_by_name = record.get_cds_name_mapping()

    for cluster_type, cds_names in cds_by_cluster_type.items():
        cds_features = sorted([cds_feature_by_name[cds] for cds in cds_names])
        rule = rules_by_name[cluster_type]
        cutoff = rule.cutoff
        extent = rule.extent
        for cds_feature_sublocation in sorted ([ cds_feature_sublocation for cds_feature in cds_features for cds_feature_sublocation in cds_feature.location.parts ], key=lambda cds_feature_sublocation: cds_feature_sublocation.start ):
            if len(borders) > 0:
                for dummy_sub_location in borders[-1].location.parts:
                    dummy_sub_location_feature = Feature(FeatureLocation(dummy_sub_location.start - cutoff, dummy_sub_location.end + cutoff), feature_type="dummy")
                    if dummy_sub_location_feature.overlaps_with(cds_feature_sublocation) and borders[-1].product == cluster_type:
                        # no need to have CompoundLocation here, the only exception would be ori-split borders, but they are handled elsewhere
                        borders[-1].location = FeatureLocation(min(dummy_sub_location.start, cds_feature_sublocation.start), max(dummy_sub_location.end, cds_feature_sublocation.end))
                    else:
                        borders.append(ClusterBorder(FeatureLocation(cds_feature_sublocation.start, cds_feature_sublocation.end), tool="rule-based-clusters",
                                                cutoff=cutoff, extent=extent, product=cluster_type))
            else:
                borders.append(ClusterBorder(FeatureLocation(cds_feature_sublocation.start, cds_feature_sublocation.end), tool="rule-based-clusters",
                                        cutoff=cutoff, extent=extent, product=cluster_type))

    for border in borders:
        border.rule = str(rules_by_name[border.product].conditions)
        if border.location.start < 0:
            border.location = FeatureLocation(0, border.location.end)
            border.contig_edge = True
        if border.location.end > len(record):
            border.location = FeatureLocation(border.location.start, len(record))
            border.contig_edge = True

    borders = remove_redundant_borders(borders, rules_by_name)

    # check if first and last borders of each metabolite were supposed to be together on a circular record
    if len(borders) > 1 and record.is_circular():
        for i in reversed(range(len(borders))): # changed according to https://stackoverflow.com/a/1207485
            for j in reversed(range(len(borders)-1)):
                if len(record) - borders[j].location.end + borders[i].location.start < max(borders[i].cutoff,borders[j].cutoff) and borders[j].product == borders[i].product:
                    borders[i].location = CompoundLocation([part for part in borders[j].location.parts] + [part for part in borders[i].location.parts])
                    # we can safely remove, because we will be never returning to [i] in none of the loops, inner or outer
                    borders.remove(borders[i])

    logging.debug("%d rule-based border(s) found in record", len(borders))
    return borders


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
                # TODO: improve performance
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


def detect_borders_and_signatures(record: Record, signature_file: str, seeds_file: str,
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
        return None
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
        record.add_cluster_border(cluster)
        cds_results = []
        cluster_extent = FeatureLocation(cluster.location.start - cluster.extent,
                                         cluster.location.end + cluster.extent)
        for cds in record.get_cds_features_within_location(cluster_extent):
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
