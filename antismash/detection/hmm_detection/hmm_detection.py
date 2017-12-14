# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

import logging
from typing import Dict, List, Set, Tuple

from Bio.SearchIO._model.hsp import HSP

from antismash.common import path, fasta
from antismash.common.secmet.feature import Cluster, SecMetQualifier, CDSFeature, GeneFunction
from antismash.common.deprecated import FeatureLocation
from antismash.common.subprocessing import run_hmmsearch
from antismash.detection.hmm_detection import rule_parser
from antismash.detection.hmm_detection.signatures import get_signature_profiles, get_signature_names


def find_clusters(record, cds_domains_by_cluster, rules_by_name) -> List[Cluster]:
    """ Detects gene clusters based on the identified core genes """
    clusters = []  # type: List[Cluster]

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
                # Finalize the last extended cluster
                cluster = clusters[-1]
                cluster.location = FeatureLocation(max(0, cluster.location.start - cluster.extent),
                                                   min(len(record), cluster.location.end + cluster.extent))
            # Create new cluster
            new_cluster = Cluster(FeatureLocation(feature_start, feature_end), feature_cutoff, feature_extension, feature_types)
            clusters.append(new_cluster)
            cluster = clusters[-1]

        # Update cluster
        cluster.location = FeatureLocation(min(cluster.location.start, feature_start), max(cluster.location.end, feature_end))
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

    # Add a note to specify whether a cluster lies on the contig/scaffold edge or not
    for cluster in clusters:
        edge = cluster.location.start == 0 or cluster.location.end == len(record)
        cluster.contig_edge = edge
        record.add_cluster(cluster)
        cluster.trim_overlapping()
    logging.info("%d cluster(s) found in record", len(clusters))
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


def filter_results(results: List[HSP], results_by_id: Dict[str, List[HSP]]) -> Tuple[List[HSP], Dict[str, List[HSP]]]:
    """ Filter results by comparing scores of different models """
    for line in open(path.get_full_path(__file__, "filterhmmdetails.txt"), "r"):
        line = line.strip()
        equivalence_group = set(line.split(","))
        unknown = equivalence_group - set(get_signature_names())
        if unknown:
            raise ValueError("Equivalence group contains unknown identifiers: %s" % (
                    unknown))
        removed_ids = set()
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


def create_rules(enabled_cluster_types: List[str]) -> List[rule_parser.DetectionRule]:
    """ Creates DetectionRule instances from the default rules file

        Args:
            enabled_cluster_types: A list of type names.
                    Only types in this list will have a DetectionRule.

        Returns:
            A list of DetectionRules.
    """
    rules = []
    # TODO: as4: We should move all user-customizable files into config subdirectory;
    # the rulefiles are redundant also in hmm_detection_dblookup
    with open(path.get_full_path(__file__, "cluster_rules.txt"), "r") as ruledata:
        parser = rule_parser.Parser("".join(ruledata.readlines()))
    # the 'other' rule is a special case, make sure it's last so we can skip it
    # if other rules hit first
    other = None
    for rule in parser.rules:
        if rule.name == "other":
            other = rule
            continue
        if rule.name in enabled_cluster_types:
            rules.append(rule)
    if other and "other" in enabled_cluster_types:
        rules.append(other)
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
        cds_domains_by_cluster_type[cds] = domains_by_cluster
    return cds_domains_by_cluster_type


def detect_signature_genes(record, options) -> None:
    """ Compares all CDS features in a record with HMM signatures and generates
        Cluster features based on those hits and the current cluster detection
        rules

        Arguments:
            record: the record to analyse
            options: antismash Config
    """
    enabled_cluster_types = options.enabled_cluster_types
    feature_by_id = record.get_cds_name_mapping()
    # if there's no genes, don't try to do anything
    if not feature_by_id:
        return None
    full_fasta = fasta.get_fasta_from_record(record)
    rules = create_rules(enabled_cluster_types)
    results = []
    sig_by_name = {}
    results_by_id = {}  # type: Dict[str, HSP]
    for sig in get_signature_profiles():
        sig_by_name[sig.name] = sig

    runresults = run_hmmsearch(path.get_full_path(__file__, 'data/bgc_seeds.hmm'),
                               full_fasta, use_tempfile=True)
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
    results, results_by_id = filter_results(results, results_by_id)

    # Filter multiple results of the same model in one gene
    results, results_by_id = filter_result_multiple(results, results_by_id)

    # Use rules to determine gene clusters
    cds_domains_by_cluster = apply_cluster_rules(results_by_id, feature_by_id, rules)

    # Find number of sequences on which each pHMM is based
    nseqdict = get_sequence_counts()

    # Save final results to record
    rules_by_name = {rule.name: rule for rule in rules}
    find_clusters(record, cds_domains_by_cluster, rules_by_name)

    for cds in results_by_id:
        feature = feature_by_id[cds]
        cluster_products = []
        if feature.cluster:
            cluster_products = feature.cluster.products
        _update_sec_met_entry(feature, results_by_id[cds],
                              cds_domains_by_cluster[cds], nseqdict,
                              cluster_products)

    # Add details of gene cluster detection to cluster features
    store_detection_details(rules_by_name, record)


def get_sequence_counts() -> Dict[str, str]:
    """ Gets the number of sequences/seeds used to generate each HMM signature

        Arguments:
            None

        Returns:
            a dictionary mapping HMM name to the number of sequences used to
                generate it
    """
    result = {}
    for hmm in get_signature_profiles():
        hmmfile = hmm.hmm_file
        for line in open(hmmfile, 'r'):
            if line.startswith('NSEQ '):
                result[hmm.name] = line[6:].strip()
                break
        # TODO: ideally this shouldn't ever happen, clean up inputs and change to error
        if hmm.name not in result:
            result[hmm.name] = "?"

    return result


def store_detection_details(rules, record) -> None:
    """ Add a note qualifier to every Cluster feature containing the detection
        type and the detection rule used to determine that type.

        Args:
            rules: A dictionary of the detection rules.
            record: The SeqRecord to modify the features of

        Returns:
            None. All changes are made in place.
    """
    for cluster in record.get_clusters():
        assert cluster.type == "cluster"
        cluster.detection_rules = [str(rules[product].conditions) for product in cluster.products]


class SecMetResult():
    """ A simple container for the information needed to create a domain """
    def __init__(self, res, nseeds):
        self.query_id = res.query_id
        self.evalue = res.evalue
        self.bitscore = res.bitscore
        self.nseeds = nseeds

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return "{} (E-value: {}, bitscore: {}, seeds: {})".format(
                self.query_id, self.evalue, self.bitscore, self.nseeds)


def _update_sec_met_entry(feature: CDSFeature, results: List[HSP],
                          definition_domains: Dict[str, Set[str]],
                          num_seeds: Dict[str, str],
                          cluster_products: List[str]) -> None:
    """ Add or updates the secondary metabolite information of a feature

        Arguments:
            feature: the feature to update
            results: a dictionary mapping gene id to HSP results
            definition_domains: a dictionary mapping cluster type to matching domains
            num_seeds: a dictionary mapping signature ids to number of seeds
                        generating the signature

        Returns:
            None
    """
    domains = [SecMetResult(res, num_seeds.get(res.query_id, "?")) for res in results]

    full_name = "-".join(list(definition_domains)) or "none"

    feature.sec_met = SecMetQualifier(full_name, domains)
    all_matching = set()
    for cluster_type, matching_domains in definition_domains.items():
        # don't add 'other' as CORE if another cluster applies
        if len(cluster_products) > 1 and cluster_type == "other":
            continue
        all_matching.update(matching_domains)
        for domain in matching_domains:
            feature.gene_functions.add(GeneFunction.CORE, "cluster_definition",
                                       "%s: %s" % (cluster_type, domain))

    # all all detected domains as ADDITIONAL if not CORE
    for secmet_domain in feature.sec_met.domains:
        if secmet_domain.query_id in all_matching:
            continue
        # skip if already added (e.g. reusing results)
        found = False
        for function in feature.gene_functions.get_by_tool("hmm_detection") or []:
            if function.description == str(secmet_domain):
                found = True
                break
        if not found:
            feature.gene_functions.add(GeneFunction.ADDITIONAL, "hmm_detection",
                                       secmet_domain)
