# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Detects specific domains and defines clusters based on domains detected
"""

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
        assert domains
        assert isinstance(definition_domains, dict), type(definition_domains)
        self.definition_domains = definition_domains
        assert self.definition_domains

    def annotate(self, products: Set[str]) -> None:
        if not self.domains or self.definition_domains:
            return
        if not self.cds.sec_met:
            self.cds.sec_met = SecMetQualifier(set(self.definition_domains), self.domains)
        else:
            self.cds.sec_met.add_products(products)
            self.cds.sec_met.add_domains(self.domains)
        all_matching = set()
        for cluster_type, matching_domains in self.definition_domains.items():
            all_matching.update(matching_domains)
            for domain in matching_domains:
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


def find_clusters(record, cds_domains_by_cluster, rules_by_name) -> List[ClusterBorder]:
    """ Detects gene clusters based on the identified core genes """
    clusters = []  # type: List[ClusterBorder]

    for feature in record.get_cds_features():
        within_cutoff = False
        feature_types = list(cds_domains_by_cluster.get(feature.get_name(), []))
        feature_start = min(feature.location.start, feature.location.end)
        feature_end = max(feature.location.start, feature.location.end)
        if not feature_types:
            continue
        feature_cutoff = -1
        feature_extension = -1
        for name in feature_types:
            rule = rules_by_name[name]
            feature_cutoff = max(feature_cutoff, rule.cutoff)
            feature_extension = max(feature_extension, rule.extent)
        cluster = None

        if clusters:
            cluster = clusters[-1]
            cluster_end = cluster.location.end
            # Check cutoff
            cutoff = max(cluster.cutoff, feature_cutoff)
            cutoff = max(cutoff, cluster.extent + feature_extension)
            within_cutoff = feature_start <= cluster_end + cutoff

        # start a new cluster if this is too far from the previous
        if not within_cutoff:
            if clusters:
                # Finalize the previous extended cluster
                cluster = clusters[-1]
                cluster.location = FeatureLocation(max(0, cluster.location.start - cluster.extent),
                                                   min(len(record), cluster.location.end + cluster.extent))
            # Create new cluster
            new_cluster = ClusterBorder(FeatureLocation(feature_start, feature_end), tool="rule-based clusters",
                                        cutoff=feature_cutoff, extent=feature_extension, products=feature_types)
            clusters.append(new_cluster)
            cluster = clusters[-1]

        # Update cluster
        start = min(cluster.location.start, feature_start)
        end = max(cluster.location.end, feature_end)
        # if the extents would overlap, but not cutoffs, then trim the later
        # cluster to start only where the previous ends
        # (occurs in CP006259, clusters 1 and 2)
        if len(clusters) > 1:
            previous = clusters[-2]
            if previous.location.end > start - cluster.extent:
                # pad with the current extent, since that will be removed when
                # the cluster is finalised
                start = previous.location.end + cluster.extent
        cluster.location = FeatureLocation(start, end)
        cluster.cutoff = max(cluster.cutoff, feature_cutoff)
        cluster.extent = max(cluster.extent, feature_extension)
        cluster.products = list(set(cluster.products) | set(feature_types))
        if len(cluster.products) > 1:
            cluster.products = sorted(list(filter(lambda prod: prod != "other", cluster.products)))

    if clusters:
        # Finalize the last extended cluster
        cluster = clusters[-1]
        extension = cluster.extent
        cluster.location = FeatureLocation(max(0, cluster.location.start - extension),
                                           min(len(record), cluster.location.end + extension))

    # TODO: Add a note to specify whether a cluster lies on the contig/scaffold edge or not

    for cluster in clusters:
        cluster.rules = [str(rules_by_name[product].conditions) for product in cluster.products]
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


def calculate_nearby_features(names: List[str], position: int, cutoff: int,
                              distances: List[int], features: Dict[str, CDSFeature],
                              results: Dict[str, List[HSP]]
                              ) -> Tuple[Dict[str, CDSFeature], Dict[str, List[HSP]]]:
    """ Find features within a specific distance

        Args:
            names: A list of CDS names.
            position: The index of the current CDS in the distances list.
            cutoff: The maximum distance to include.
            distances: An ordered list of distances from the first CDS in the list
            features: A dict of feature ID to Feature instance
            results: A dict of feature ID to list of HSP results

        Returns:
            A tuple of dictionaries. The dictionaries are strictly subsets of
            the features and results args.
    """
    nearby_features = {names[position]: features[names[position]]}

    center = distances[position]

    before = position - 1
    while before >= 0 and center - distances[before] <= cutoff:
        neighbour = names[before]
        nearby_features[neighbour] = features[neighbour]
        before -= 1

    after = position + 1
    while after < len(distances) and distances[after] - center <= cutoff:
        neighbour = names[after]
        nearby_features[neighbour] = features[neighbour]
        after += 1

    nearby_results = {name: results[name] for name in nearby_features}

    return nearby_features, nearby_results


def apply_cluster_rules(results_by_id: Dict[str, HSP], feature_by_id: Dict[str, CDSFeature],
                        rules) -> Dict[str, Dict[str, Set[str]]]:
    """
        Run detection rules over each CDS and classify them if relevant.
        A CDS can satisfy multiple rules. If so, all rules satisfied
        will form part of the type string, separated by '-'.

        The 'other' type has a lower precedence than other rules and a hit with
        the 'other' rule will be ignored if another rule is also satisfied.

        Args:
            results_by_id: A dict of CDS ID to a list of HSP results
            feature_by_id: A dict of CDS ID to CDS Feature
            rules: A list of DetectionRule instances

        Returns:
            A dictionary mapping CDS ID to
                a dictionary mapping cluster type string to
                    a set of domains used to determine the cluster.
    """
    if not results_by_id:
        return {}

    cds_with_hits = sorted(results_by_id, key=lambda gene_id: feature_by_id[gene_id].location.start)

    def calculate_distance(first, second):
        """ Calculate the distance between two FeatureLocations """
        first_start, first_end = sorted([first.start, first.end])
        second_start, second_end = sorted([second.start, second.end])
        return min(abs(first_end - second_start), abs(second_end - first_start),
                   abs(first_start - second_start), abs(first_end - second_end))

    first_cds = feature_by_id[cds_with_hits[0]]
    distances = [calculate_distance(first_cds.location, feature_by_id[cds].location) for cds in cds_with_hits]

    cds_domains_by_cluster_type = {}

    for i, cds in enumerate(cds_with_hits):
        results = []  # type: List[str]
        rule_texts = []
        info_by_range = {}  # type: Dict[int, Tuple[Dict[str, CDSFeature], Dict]]
        domain_matches = set()  # type: Set[str]
        domains_by_cluster = {}  # type: Dict[str, Set[str]]
        for rule in rules:
            if results and rule.name == "other":
                continue
            if rule.cutoff not in info_by_range:
                info_by_range[rule.cutoff] = calculate_nearby_features(cds_with_hits,
                        i, rule.cutoff, distances, feature_by_id, results_by_id)
            nearby_features, nearby_results = info_by_range[rule.cutoff]
            matching = rule.detect(cds, nearby_features, nearby_results)
            if matching.met and matching.matches:
                domains_by_cluster[rule.name] = matching.matches
                results.append(rule.name)
                rule_texts.append(rule.reconstruct_rule_text())
                domain_matches.update(matching.matches)
        if domains_by_cluster:
            cds_domains_by_cluster_type[cds] = domains_by_cluster
    return cds_domains_by_cluster_type


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
    cds_domains_by_cluster = apply_cluster_rules(results_by_id, feature_by_id, rules)

    # Find number of sequences on which each pHMM is based
    num_seeds_per_hmm = get_sequence_counts(signature_file)

    # Save final results to record
    rules_by_name = {rule.name: rule for rule in rules}
    clusters = find_clusters(record, cds_domains_by_cluster, rules_by_name)

    cds_results_by_cluster = {}
    for cluster in clusters:
        record.add_cluster_border(cluster)
        cds_results = []
        for cds in record.get_cds_features_within_location(cluster.location):
            domains = [Domain(res, num_seeds_per_hmm.get(res.query_id, "?"), tool) for res in results_by_id.get(cds.get_name(), [])]
            if domains and cds.get_name() in cds_domains_by_cluster:
                cds_results.append(CDSResults(cds, domains, cds_domains_by_cluster[cds.get_name()]))
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
