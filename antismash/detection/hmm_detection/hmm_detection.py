# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

import logging
from collections import defaultdict
from typing import Dict, Tuple, TypeVar

from antismash.common import path, deprecated
from antismash.common.secmet.feature import Cluster, SecMetQualifier, CDSFeature, GeneFunction
from antismash.common.deprecated import FeatureLocation
from antismash.common.subprocessing import run_hmmsearch
from antismash.detection.hmm_detection import rule_parser
from antismash.detection.hmm_detection.signatures import get_signature_profiles, get_signature_names

T = TypeVar("T")


def find_clusters(record, rules):
    """ Detects gene clusters based on the identified core genes """
    clusters = []

    for feature in record.get_cds_features():
        within_cutoff = False
        if not feature.sec_met:
            continue
        feature_type = feature.sec_met.clustertype
        feature_start = min(feature.location.start, feature.location.end)
        feature_end = max(feature.location.start, feature.location.end)
        if feature_type == "none":
            continue
        feature_cutoff = -1
        feature_extension = -1
        for name in feature_type.split("-"):
            rule = rules[name]
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
            new_cluster = Cluster(FeatureLocation(feature_start, feature_end), feature_cutoff, feature_extension, feature_type.split("-"))
            clusters.append(new_cluster)
            cluster = clusters[-1]

        # Update cluster
        cluster.location = FeatureLocation(min(cluster.location.start, feature_start), max(cluster.location.end, feature_end))
        cluster.cutoff = max(cluster.cutoff, feature_cutoff)
        cluster.extent = max(cluster.extent, feature_extension)
        cluster.products = list(set(cluster.products) | set(feature_type.split('-')))
        if len(cluster.products) > 1:
            cluster.products = list(filter(lambda prod: prod != "other", cluster.products))

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
    logging.info("%d cluster(s) found in record", len(clusters))


def hsp_overlap_size(first, second):
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


def filter_results(results, results_by_id):
    """ Filter results by comparing scores of different models """
    for line in open(path.get_full_path(__file__, "filterhmmdetails.txt"), "r"):
        line = line.strip()
        equivalence_group = set(line.split(","))
        unknown = equivalence_group - set(get_signature_names())
        if unknown:
            raise ValueError("Equivalence group contains unknown identifiers: %s" % (
                    unknown))
        for cds, cdsresults in results_by_id.items():
            # Check if multiple competing HMM hits are present
            hits = set(hit.query_id for hit in cdsresults)
            if len(hits & equivalence_group) < 2:
                continue
            # Identify overlapping hits
            overlapping_groups = []
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

            for group in overlapping_groups:
                # find the best
                best = None
                for hit in group:
                    if not best or hit.bitscore > best.bitscore:
                        best = hit
                # remove the rest
                for hit in group:
                    if hit != best:
                        del results[results.index(hit)]
                        del results_by_id[cds][results_by_id[cds].index(hit)]
            assert results_by_id[cds]  # should always have one remaining
    return results, results_by_id


def filter_result_multiple(results, results_by_id):
    """ Filter multiple results of the same model within a gene """
    for cds, hits in results_by_id.items():
        query_scores = {}
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


def filter_result_overlapping_genes(results, results_by_id, overlaps, feature_by_id):
    """ Filter results of overlapping genes
        (only gene with the best score can retain its result)
    """
    assert all(feature.type == "CDS" for feature in feature_by_id.values())
    filterhmm_list = []
    overlap_id_with_result = {}
    for line in open(path.get_full_path(__file__, "filterhmmdetails.txt"), "r").read().split("\n"):
        filterhmms = line.split(",")
        if filterhmms not in filterhmm_list:
            filterhmm_list.append(filterhmms)
    for cds in results_by_id:
        if overlaps[cds] not in overlap_id_with_result:
            overlap_id_with_result[overlaps[cds]] = [cds]
        elif cds not in overlap_id_with_result[overlaps[cds]]:
            overlap_id_with_result[overlaps[cds]].append(cds)
    for overlap_id in overlap_id_with_result:
        best_hit_scores = {}
        for cds in overlap_id_with_result[overlap_id]:
            for hit in results_by_id[cds]:
                feature = feature_by_id[hit.hit_id]
                dist = abs(feature.location.end - feature.location.start)
                if hit.query_id not in best_hit_scores or best_hit_scores[hit.query_id] < dist:
                    best_hit_scores[hit.query_id] = dist
        for cds in overlap_id_with_result[overlap_id]:
            to_delete = []
            for hit in results_by_id[cds]:
                feature = feature_by_id[hit.hit_id]
                dist = abs(feature.location.end - feature.location.start)
                if dist < best_hit_scores[hit.query_id]:
                    to_delete.append(hit)
                else:  # filter for filterhmmdetails.txt
                    for filterhmms in filterhmm_list:
                        if hit.query_id not in filterhmms:
                            continue
                        for similar_hit in filterhmms:
                            if similar_hit not in best_hit_scores:
                                continue
                            if dist < best_hit_scores[similar_hit]:
                                to_delete.append(hit)
                                break
            for hit in to_delete:
                del results[results.index(hit)]
                del results_by_id[cds][results_by_id[cds].index(hit)]
                if len(results_by_id[cds]) < 1:
                    del results_by_id[cds]
    return results, results_by_id


def create_rules(enabled_cluster_types):
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
    for rule in parser.rules:
        if rule.name in enabled_cluster_types:
            rules.append(rule)
    return rules


def calculate_nearby_features(names, position, cutoff, distances, features,
                              results: Dict[str, T]) -> Tuple[Dict[str, CDSFeature], Dict[str, T]]:
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


def apply_cluster_rules(results_by_id, feature_by_id, rules) -> Dict[str, str]:
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
            A dictionary of CDS ID to type string.
    """
    if not results_by_id:
        return defaultdict(lambda: "none")

    type_results = {}
    cds_with_hits = sorted(results_by_id, key=lambda gene_id: feature_by_id[gene_id].location.start)

    def calculate_distance(first, second):
        """ Calculate the distance between two FeatureLocations """
        first_start, first_end = sorted([first.start, first.end])
        second_start, second_end = sorted([second.start, second.end])
        return min(abs(first_end - second_start), abs(second_end - first_start),
                   abs(first_start - second_start), abs(first_end - second_end))

    first_cds = feature_by_id[cds_with_hits[0]]
    distances = [calculate_distance(first_cds.location, feature_by_id[cds].location) for cds in cds_with_hits]

    for i, cds in enumerate(cds_with_hits):
        results = []
        info_by_range = {}
        for rule in rules:
            if rule.cutoff not in info_by_range:
                info_by_range[rule.cutoff] = calculate_nearby_features(cds_with_hits,
                        i, rule.cutoff, distances, feature_by_id, results_by_id)
            nearby_features, nearby_results = info_by_range[rule.cutoff]
            if rule.detect(cds, nearby_features, nearby_results):
                results.append(rule.name)
        if len(results) > 1 and "other" in results:
            results.remove("other")
        if results:
            type_results[cds] = "-".join(results)
            feature = feature_by_id[cds]
            feature.gene_functions.add(GeneFunction.CORE, "cluster_definition", results_by_id[cds][0].query_id)
        else:
            type_results[cds] = "none"
    return type_results


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
    full_fasta = deprecated.get_multifasta(record)
    rules = create_rules(enabled_cluster_types)
    results = []
    sig_by_name = {}
    results_by_id = {}
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

    # Get overlap tables (for overlap filtering etc)
    overlaps = get_overlaps_table(record)

    # Filter results by comparing scores of different models (for PKS systems)
    results, results_by_id = filter_results(results, results_by_id)

    # Filter results of overlapping genes (only for plants)
    if options.taxon == 'plants':
        results, results_by_id = filter_result_overlapping_genes(results,
                                         results_by_id, overlaps, feature_by_id)

    # Filter multiple results of the same model in one gene
    results, results_by_id = filter_result_multiple(results, results_by_id)

    # Use rules to determine gene clusters
    typedict = apply_cluster_rules(results_by_id, feature_by_id, rules)

    # Find number of sequences on which each pHMM is based
    nseqdict = get_sequence_counts()

    # Save final results to record
    for cds in results_by_id:
        feature = feature_by_id[cds]
        _update_sec_met_entry(feature, results_by_id[cds], typedict[cds], nseqdict)

    rules = {rule.name: rule for rule in rules}
    find_clusters(record, rules)

    # Find additional NRPS/PKS genes in gene clusters
    add_additional_nrpspks_genes(typedict, results_by_id, record, nseqdict)
    # Add details of gene cluster detection to cluster features
    store_detection_details(rules, record)
    # If all-orfs option on, remove irrelevant short orfs
    if options.genefinding_tool == "all-orfs":
        remove_irrelevant_allorfs(record)


def get_sequence_counts() -> Dict[str, int]:
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
        # TODO: ideally this shouldn't ever happen
        if hmm.name not in result:
            result[hmm.name] = "?"

    return result


def remove_irrelevant_allorfs(record):
    """ Remove auto-orf features without unique sec_met qualifiers;
        remove glimmer ORFs overlapping with sec_met auto-orfs not caught by Glimmer
    """
    logging.critical("remove_irrelevant_allorfs() not yet tested")
    auto_orf_features = []
    other_features = []
    for feature in record.get_cds_features():
        if 'auto-all-orf' in feature.notes:
            auto_orf_features.append(feature)
        else:
            other_features.append(feature)
    to_delete = []
    for autofeature in auto_orf_features:
        if "sec_met" not in autofeature.qualifiers:
            to_delete.append(autofeature)
            continue
        glimmer_has_sec_met = False
        for otherfeature in other_features:
            if autofeature.overlaps_with(otherfeature) and otherfeature.sec_met:
                to_delete.append(autofeature)
                glimmer_has_sec_met = True
        if not glimmer_has_sec_met:
            for otherfeature in other_features:
                if autofeature.overlaps_with(otherfeature) and not otherfeature.secmet:
                    to_delete.append(otherfeature)
    logging.critical("attempting to remove features in remove_irrelevant_allorfs")
    featurenrs = []
    idx = 0
    for feature in record.features:
        if feature in to_delete:
            featurenrs.append(idx)
        idx += 1
    featurenrs.reverse()
    for featurenr in featurenrs:
        del record.features[featurenr]


def add_additional_nrpspks_genes(typedict, results_by_id, record, nseqdict):
    logging.critical("add_additional_nrpspks_genes() overriding existing sec_met")
    nrpspksdomains = ["PKS_KS", "PKS_AT", "ATd", "ene_KS", "mod_KS", "hyb_KS",
                      "itr_KS", "tra_KS", "Condensation", "AMP-binding", "A-OX"]
    clustercdsfeatures = record.get_cds_features_within_clusters()
    othercds_with_results = []
    for cds in clustercdsfeatures:
        gene_id = cds.get_name()
        if gene_id in results_by_id and typedict[gene_id] == "none":
            othercds_with_results.append((cds, gene_id))
    for cds, gene_id in othercds_with_results:
        cdsresults = [res.query_id for res in results_by_id[gene_id]]
        if len(set(nrpspksdomains) & set(cdsresults)) >= 1:
            _update_sec_met_entry(cds, results_by_id[gene_id], "other", nseqdict)


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


def _update_sec_met_entry(feature, results, clustertype, num_seeds) -> None:
    """ Add or updates the secondary metabolite information of a feature

        Arguments:
            feature: the feature to update
            results: a dictionary mapping gene id to HSP results
            clustertype: a string containing cluster types separated by -
            num_seeds: a dictionary mapping signature ids to number of seeds
                        generating the signature

        Returns:
            None
    """
    domains = [SecMetResult(res, num_seeds.get(res.query_id, "?")) for res in results]

    feature.sec_met = SecMetQualifier(clustertype, domains)
    if not feature.gene_functions.get_by_tool("cluster_definition"):
        feature.gene_functions.add(GeneFunction.ADDITIONAL, "hmm_detection", feature.sec_met.domains[0])


def get_overlaps_table(record) -> Dict[str, int]:
    """
        Identify overlapping genes by overlap group.
        Genes with no overlap will have their own group.

        Args:
            record: The record to search.

        Returns:
            A dictionary of gene id to int representing overlap group

    """
    overlaps = []
    overlap_by_id = {}
    features = record.get_cds_features()
    if len(features) < 1:
        return overlap_by_id
    features = sorted(features, key=lambda feature: feature.location.start)
    i = 0
    j = i + 1
    cds_queue = []
    while i <= j < len(features):
        cds = features[i]
        ncds = features[j]
        # TODO: ensure this works with opposite strands of genes
        if cds.location.end <= ncds.location.start + 1:
            overlaps.append([])
            cds_queue.append(cds)
            for cds in cds_queue:
                overlap_by_id[cds.get_name()] = len(overlaps) - 1
                overlaps[-1].append(cds)
            cds_queue = []
            i = j
        else:
            if cds.location.end < ncds.location.end:
                cds_queue.append(cds)
                i = j
            else:
                cds_queue.append(ncds)
        j += 1
    overlaps.append([])
    cds_queue.append(features[i])
    for cds in cds_queue:
        overlap_by_id[cds.get_name()] = len(overlaps) - 1
        overlaps[-1].append(cds)
    return overlap_by_id
