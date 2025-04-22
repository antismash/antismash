# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Detects specific domains and defines clusters based on domains detected
"""

import bisect
from collections import defaultdict
import dataclasses
import logging
import math
from typing import Any, Dict, FrozenSet, Iterable, Iterator, List, Optional, Set, Tuple, Union

from antismash.common import fasta, serialiser
from antismash.common.hmmscan_refinement import HSP
from antismash.common.secmet import CDSFeature, Feature, FeatureLocation, Protocluster, Record
from antismash.common.secmet.locations import (
    CompoundLocation,
    Location,
    location_bridges_origin,
    location_contains_other,
    locations_overlap,
    make_forwards,
)
from antismash.common.secmet.qualifiers import GeneFunction, SecMetQualifier
from antismash.common.subprocessing import run_hmmsearch
from antismash.common.hmm_rule_parser import rule_parser
from antismash.common.signature import get_signature_profiles, HmmSignature, Signature

from .structures import DynamicHit, DynamicProfile, HMMerHit, Multipliers, ProfileHit


GenericSets = Union[Iterable[set[str]], Iterable[FrozenSet[str]]]


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

    schema_version = 4

    def __init__(self, cds_by_cluster: Dict[Protocluster, List[CDSResults]],
                 tool: str, cdses_outside_clusters: List[CDSResults],
                 multipliers: Multipliers) -> None:
        self.cds_by_cluster = cds_by_cluster
        self.tool = str(tool)
        self.cdses_outside_clusters = cdses_outside_clusters
        self.multipliers = multipliers

    @property
    def protoclusters(self) -> List[Protocluster]:
        """ A list of Protoclusters predicted """
        return list(self.cds_by_cluster)

    def annotate_cds_features(self) -> None:
        """ Annotate relevant CDS features with the HMM information detected """
        for cds_results in self.cds_by_cluster.values():
            for cds_result in cds_results:
                cds_result.annotate(self.tool)
        for cds_result in self.cdses_outside_clusters:
            cds_result.annotate(self.tool)

    def to_json(self) -> Dict[str, Any]:
        """ Constructs a JSON representation from the RuleDetectionResults instance """
        cds_results_json: List[Tuple[Dict[str, Any], List[Dict[str, Any]]]] = []
        json = {
            "schema_version": self.schema_version,
            "tool": self.tool,
            "cds_by_protocluster": cds_results_json,
            "outside_protoclusters": [result.to_json() for result in self.cdses_outside_clusters],
            "multipliers": dataclasses.asdict(self.multipliers),
        }

        for cluster, cds_results in self.cds_by_cluster.items():
            json_cluster = serialiser.feature_to_json(cluster.to_biopython()[0])
            json_cds_results = [result.to_json() for result in cds_results]
            cds_results_json.append((json_cluster, json_cds_results))

        return json

    @staticmethod
    def from_json(json: Dict[str, Any], record: Record) -> Optional["RuleDetectionResults"]:
        """ Constructs a RuleDetectionResults instance from a JSON representation """
        if RuleDetectionResults.schema_version != json.get("schema_version", 1):
            return None

        cds_by_cluster = {}
        for json_cluster, json_cds_results in json["cds_by_protocluster"]:
            cluster = Protocluster.from_biopython(serialiser.feature_from_json(json_cluster))
            cds_results = [CDSResults.from_json(result_json, record) for result_json in json_cds_results]
            cds_by_cluster[cluster] = cds_results

        cdses_outside = [CDSResults.from_json(chunk, record) for chunk in json["outside_protoclusters"]]
        multipliers = Multipliers(**json["multipliers"])

        return RuleDetectionResults(cds_by_cluster, json["tool"], cdses_outside, multipliers)


def get_equivalence_groups_from_file(filename: str, signature_names: set[str]) -> list[set[str]]:
    """ Reads a file of equivalence groups.

        Arguments:
            filename: the path to the file containing the equivalence sets
            signature_names: the set of known signature names that the file will be
                             referring to

        Returns:
            a list of sets, one for each equivalence group.
    """
    groups = []
    with open(filename, "r", encoding="utf-8") as handle:
        for line in handle.readlines():
            line = line.strip()
            equivalence_group = set(line.split(","))
            unknown = equivalence_group - signature_names
            if unknown:
                raise ValueError(f"Equivalence group contains unknown identifiers: {unknown}")
            groups.append(equivalence_group)
    return groups


@dataclasses.dataclass
class Ruleset:
    """ A container for all information related to a particular set of rules """
    _rules: tuple[rule_parser.DetectionRule, ...]
    hmm_profiles: dict[str, HmmSignature]
    database_file: str
    valid_categories: set[str]
    tool: str
    multipliers: Multipliers = dataclasses.field(default_factory=Multipliers)
    dynamic_profiles: dict[str, DynamicProfile] = dataclasses.field(default_factory=dict)
    # the following should really be keyword only, but that needs 3.10+
    equivalence_groups: dataclasses.InitVar[GenericSets] = None
    filter_file: dataclasses.InitVar[str] = None

    def __post_init__(self, equivalence_groups: Iterable[Union[set[str], frozenset[str]]] = None,
                      filter_file: str = None) -> None:
        # add an additional field to allow fetching rules by name
        self._rules_by_name = {rule.name: rule for rule in self.rules}
        # if there were duplicate names, then report them
        if len(self._rules_by_name) != len(self.rules):
            counts: dict[str, int] = defaultdict(int)
            for rule in self.rules:
                counts[rule.name] += 1
            duplicated = {name for name, count in counts.items() if count > 1}
            raise ValueError(f"rule names must be unique: {duplicated}")

        # and another additional field that is a combination of both HMM and dynamic profiles
        if not self.hmm_profiles and not self.dynamic_profiles:
            raise ValueError("at least one HMM or dynamic profile is required")
        self.all_profiles: dict[str, Signature] = dict(self.hmm_profiles)
        self.all_profiles.update(self.dynamic_profiles)
        # again, find and report and duplicate names
        overlaps = set(self.hmm_profiles).intersection(set(self.dynamic_profiles))
        if overlaps:
            raise ValueError(f"HMM profiles and dynamic profiles overlap: {overlaps}")

        # build the equivalence groups from file if necessary
        if filter_file is not None:
            if equivalence_groups is not None:
                raise ValueError("Only one of 'filter_file' and 'equivalence_groups' can be provided")
            equivalence_groups = get_equivalence_groups_from_file(filter_file, set(self.all_profiles))
        if equivalence_groups is None:
            raise ValueError("One of 'filter_file' or 'equivalence_groups' must be provided")
        self._equivalence_groups = tuple(frozenset(group) for group in equivalence_groups)

        # update ranges by provided multipliers
        for rule in self._rules_by_name.values():
            rule.cutoff = int(rule.cutoff * self.multipliers.cutoff)
            rule.neighbourhood = int(rule.neighbourhood * self.multipliers.neighbourhood)

    @property
    def rules(self) -> tuple[rule_parser.DetectionRule, ...]:
        """ Returns the rules available within the ruleset """
        return self._rules

    def get_equivalence_groups(self) -> tuple[FrozenSet[str], ...]:
        """ Returns the profile equivalence groups use to filter out similar profiles
        """
        return self._equivalence_groups

    def get_rule_by_name(self, name: str) -> rule_parser.DetectionRule:
        """ Returns the rule with the given name.
        """
        rule = self._rules_by_name.get(name)
        if not rule:
            raise ValueError(f"Unknown rule: {name}")
        return rule

    def get_rule_names(self) -> set[str]:
        """ Returns a set of all rule names within the ruleset
        """
        return set(self._rules_by_name)

    def copy_with_replacements(self, **kwargs: Any) -> "Ruleset":
        """ Copies this instance, with those values supplied replacing the existing values.
        """
        # dataclasses.replace handles everything but InitVar members, which need to be explicitly provided
        # some original init-only members have to be grabbed from elsewhere in the instance,
        # so set what they would be here
        defaults = {
            "equivalence_groups": list(self._equivalence_groups),
        }
        if "rules" in kwargs:
            kwargs["_rules"] = kwargs.pop("rules")
        # then find all the relevant fields and use either the specific default above
        # or the default the field defines
        for name, field in self.__dataclass_fields__.items():  # pylint: disable=no-member # at least for 2.15
            if isinstance(field.type, dataclasses.InitVar):
                if name not in kwargs:
                    kwargs[name] = defaults.get(name, field.default)

        return dataclasses.replace(self, **kwargs)

    @classmethod
    def from_files(cls, signature_file: str, seeds: str, rule_files: list[str],
                   categories: set[str], filter_file: str, tool: str,
                   *,
                   dynamic_profiles: dict[str, DynamicProfile] = None,
                   multipliers: Multipliers = None,
                   ) -> "Ruleset":
        """ Generates an instance directly from files of signatures, equivalance groups,
            and rules

            Arguments:
                signature_file: the path to a file containing signature information
                seeds: the path to an HMM file containing all profiles mentioned in signatures
                rule_files: one or more files containing rules, ordered so that those with
                            no dependency on the others come first
                categories: a set of categories referenced by the rules
                filter_file: a path to a file containing equivalence groups
                tool: the name of the tool that will be using the ruleset

                *Optional arguments*
                dynamic_profiles: any dynamic profiles used, as a mapping of profile
                                  name to profile instance
                multipliers: rule distance multipliers
        """
        multipliers = multipliers or Multipliers()
        dynamic_profiles = dynamic_profiles or {}
        signatures = {sig.name: sig for sig in get_signature_profiles(signature_file)}
        all_signatures = set(dynamic_profiles).union(signatures)
        rules = create_rules(rule_files, all_signatures, categories, multipliers)

        return cls(tuple(rules), signatures, seeds, categories, tool,
                   dynamic_profiles=dynamic_profiles,
                   multipliers=multipliers,
                   filter_file=filter_file,
                   )


def remove_redundant_protoclusters(clusters: List[Protocluster],
                                   rules_by_name: Dict[str, rule_parser.DetectionRule],
                                   record: Record,
                                   ) -> List[Protocluster]:
    """ Removes clusters which have superiors covering the same (or larger) region
    """
    clusters_by_rule: Dict[str, List[Protocluster]] = defaultdict(list)
    for cluster in clusters:
        clusters_by_rule[cluster.product].append(cluster)

    def get_first_and_last(cluster: Protocluster) -> tuple[CDSFeature, CDSFeature]:
        cores = record.get_cds_features_within_location(cluster.core_location)
        return cores[0], cores[-1]

    trimmed_clusters = []
    for cluster in clusters:
        rule_name = cluster.product
        is_redundant = False
        first_core_cds, last_core_cds = get_first_and_last(cluster)
        for superior in rules_by_name[rule_name].superiors:
            for other_cluster in clusters_by_rule.get(superior, []):
                if other_cluster.core_location.contains(cluster.core_location):
                    is_redundant = True
                    continue
                other_first, other_last = get_first_and_last(other_cluster)
                if other_last < first_core_cds:
                    continue
                if last_core_cds < other_first:
                    continue
                is_redundant = True
                break
            if is_redundant:
                break
        if not is_redundant:
            trimmed_clusters.append(cluster)
    return trimmed_clusters


def apply_extenders(clusters: list[Protocluster],
                    rules_by_name: dict[str, rule_parser.DetectionRule],
                    record: Record, results_by_id: dict[str, list[ProfileHit]],
                    cds_domains_by_cluster: dict[str, dict[str, set[str]]],
                    ) -> list[Protocluster]:
    """ Extends the given protoclusters using the defining rule's extension
        information, if present and possible.

        Arguments:
            clusters: the protoclusters to try extending
            rules_by_name: a mapping of rule name to DetectionRule instance
            record: the record the clusters belong to
            results_by_id: a mapping of CDS name to profile hits for the full record
            cds_domains_by_cluster: a mapping of CDS ID to
                    a mapping of cluster type string to
                        a set of domains used to determine the cluster

        Returns:
            a new list of extended clusters
    """

    def mark_extendable(cdses: Iterable[CDSFeature], previous: CDSFeature,
                        rule: rule_parser.DetectionRule, core: FeatureLocation,
                        ) -> Iterator[CDSFeature]:
        for cds in cdses:
            # skip over any existing cores
            if cds.is_contained_by(core):
                continue
            # stop if last match is too far away from the current CDS
            if record.get_distance_between_locations(cds.location, previous.location) > rule.cutoff:
                break
            extendable = rule.can_extend_to(cds, results_by_id.get(cds.get_name(), []))
            if extendable:
                previous = cds  # update the previous match for distance checking
                cds_domains_by_cluster[cds.get_name()][rule.name].update(extendable.matches)
                yield cds  # let the caller do what it wants with the CDS

    def cycle(items: tuple[CDSFeature, ...], index: int, direction: int = 1) -> Iterator[CDSFeature]:
        if direction == -1:
            sections: list[Iterable[CDSFeature]] = [reversed(items[:index]), reversed(items[index+1:])]
        else:
            sections = [items[index:], items[:index]]

        if not record.is_circular():
            sections = [sections[0]]

        for section in sections:
            for item in section:
                yield item

    results = []
    cdses = record.get_cds_features()
    for cluster in clusters:
        core = cluster.core_location
        rule = rules_by_name[cluster.product]
        index = bisect.bisect_left(cdses, core)

        core_cdses = record.get_cds_features_within_location(cluster.core_location)

        # since the extender conditions might not be a subset of the main rule,
        # mark any contained and unmarked genes that satisfy the extender condition
        for cds in core_cdses:
            if cds in cluster.definition_cdses:
                continue
            extendable = rule.can_extend_to(cds, results_by_id.get(cds.get_name(), []))
            if extendable:
                cds_domains_by_cluster[cds.get_name()][rule.name].update(extendable.matches)

        # for each direction, the domain/result updates will be handled by the marker
        # but since the core location updates are direction dependent, they still must be handled

        for cds in mark_extendable(cycle(cdses, index, direction=-1), core_cdses[0], rule, core):
            core = record.connect_locations([cds.location, core])

        for cds in mark_extendable(cycle(cdses, index), core_cdses[-1], rule, core):
            core = record.connect_locations([cds.location, core])

        # the new core should be at least as large as the old core
        assert location_contains_other(core, cluster.core_location), f"{cluster.core_location} not in {core}"
        # update locations
        surrounds = _extend_area_location(core, rule.neighbourhood, record, force_cross_origin=True)
        results.append(Protocluster(core, surrounding_location=surrounds,
                                    tool="rule-based-clusters", cutoff=rule.cutoff,
                                    neighbourhood_range=rule.neighbourhood, product=rule.name,
                                    detection_rule=str(rule.conditions), product_category=rule.category))
    return results


def _extend_area_location(location: Location, distance: int, record: Record,
                          *, force_cross_origin: bool = False) -> Location:
    """ A restrictive wrapper of Record.extend_location(), specifically for areas
        and not CDS features. No area should have more than two parts when
        overlapping the origin in circular genomes and in non-circular genomes,
        there will only ever be one part.

        Arguments:
            location: the location to extend, at most two parts and in the forward strand
            distance: the distance to extend the location, both before and after
            record: the record the location belongs in
            force_cross_origin: whether to keep the result as a cross-origin location,
                                in cases where the initial location crosses the origin
                                and the resulting location would cover the full record

        Returns:
            a new location, covering the requested distance, with two parts iff it crosses the origin
    """
    # make sure that the incoming location is acceptable
    assert location.strand is not None or len(location.parts) == 1
    max_parts = 1
    if record.is_circular():
        max_parts = 2
        if len(location.parts) > 1 and not location_bridges_origin(location):
            raise ValueError("Areas with existing compound locations must bridge the origin")
        distance = min(distance, (len(record) - len(location)) // 2 + 1)

    if len(location.parts) > max_parts:
        raise ValueError("Area has too many sub-locations for record type")

    result = location
    # the extension will be much simpler if the strand is always fowards
    if location.strand != 1:
        result = make_forwards(location)

    # extend the location
    result = record.extend_location(result, distance)

    if location.crosses_origin() and len(result) == len(record) \
            and not result.crosses_origin() and force_cross_origin:
        # set the edges of the location as the midpoint of the initial location's end and start
        start = location.parts[0].start
        end = location.parts[-1].end
        padding = math.floor((start - end) / 2)
        mid = padding + end
        result = CompoundLocation([
            FeatureLocation(mid, len(record), result.strand),
            FeatureLocation(0, mid - 1, result.strand),
        ])

    # too many parts means something has gone wrong prior to this function
    if len(result.parts) > max_parts:
        raise ValueError("Area has too many sub-locations for record type")

    return result


def find_protoclusters(record: Record, cds_by_cluster_type: Dict[str, Set[str]],
                       rules_by_name: dict[str, rule_parser.DetectionRule],
                       results_by_id: dict[str, list[ProfileHit]],
                       cds_domains_by_cluster: dict[str, dict[str, set[str]]],
                       ) -> list[Protocluster]:
    """ Detects gene clusters based on the identified core genes

        Arguments:
            record: the record to find protoclusters within
            cds_by_cluster_type: a dictionary mapping rule name to
                                 a set of CDS names that matched that rule
            rules_by_name: a dictionary mapping rule names to DetectionRule instances
            results_by_id: a dictinoary mapping CDS name to a list of profile hits for that CDS
            cds_domains_by_cluster: a mapping of CDS ID to
                    a mapping of cluster type string to
                        a set of domains used to determine the cluster

        Returns:
            a list of Protocluster instances that were found
    """
    clusters: List[Protocluster] = []

    for cluster_type, cds_names in cds_by_cluster_type.items():
        cores = []
        # because sets aren't ordered, sort the features for consistency
        cds_features = sorted(record.get_cds_by_name(cds) for cds in cds_names)
        cross_origin = (feature for feature in cds_features if location_bridges_origin(feature.location))
        cds_features = sorted([feature for feature in cds_features if not location_bridges_origin(feature.location)])
        rule = rules_by_name[cluster_type]
        cutoff = rule.cutoff
        for cross in cross_origin:
            cores.append(Feature(record.connect_locations(cross.location.parts), "temp"))

        while cds_features:
            # since the features are sorted already, groups will automatically connect
            # except for possibly the first and last in a circular record
            cds = cds_features.pop(0)
            if not cores:
                cores.append(Feature(record.connect_locations([cds.location]), "temp"))
                continue
            previous = cores[-1]
            dummy = Feature(_extend_area_location(previous.location, cutoff, record), "temp")
            assert len(dummy.location) >= len(previous.location)
            if cds.overlaps_with(dummy):
                previous.location = record.connect_locations([previous.location, cds.location])
                continue
            cores.append(Feature(record.connect_locations([cds.location]), "temp"))
        assert cores
        # check for over-origin clusters within cutoff, which should only be first and last cores
        first = cores[0]
        last = cores[-1]
        if record.is_circular() and first is not last and first.location.start > last.location.start:
            if record.get_distance_between_features(first, last) < cutoff:
                last = cores.pop()
                first.location = record.connect_locations([last.location, first.location])
        # form a protocluster for each core
        for core in cores:
            core_location = core.location
            surrounds = _extend_area_location(core_location, rule.neighbourhood,
                                              record, force_cross_origin=True)
            clusters.append(Protocluster(core_location, surrounding_location=surrounds,
                                         tool="rule-based-clusters", cutoff=cutoff,
                                         neighbourhood_range=rule.neighbourhood, product=cluster_type,
                                         detection_rule=str(rule.conditions), product_category=rule.category))

    clusters = apply_extenders(clusters, rules_by_name, record, results_by_id, cds_domains_by_cluster)

    clusters = remove_redundant_protoclusters(clusters, rules_by_name, record)

    clusters = merge_over_origin(clusters, record)

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


def filter_results(results: list[HSP], results_by_id: Dict[str, list[HSP]],
                   equivalence_groups: GenericSets,
                   ) -> tuple[list[HSP], dict[str, list[HSP]]]:
    """ Filter results by comparing scores of different models """
    for equivalence_group in equivalence_groups:
        removed_ids: Set[int] = set()
        for cds, cdsresults in results_by_id.items():
            # Check if multiple competing HMM hits are present
            hits = set(hit.query_id for hit in cdsresults)
            if len(hits & equivalence_group) < 2:
                continue
            # Identify overlapping hits
            overlapping_groups: List[Set[HSP]] = []
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
        query_scores: Dict[str, Tuple[int, HSP, float]] = {}
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


def merge_over_origin(clusters: list[Protocluster], record: Record) -> list[Protocluster]:
    """ Merges all protoclusters of the same type that have overlapping core locations,
        or the cores are within cutoff distance of each other.

        Arguments:
            clusters: all the clusters to consider for merging
            record: the parent record for circular operations

        Returns:
            a new collection of non-overlapping protoclusters
    """
    # gather up protoclusters of the same type, as pairs of cluster and core location + cutoff
    by_product: dict[str, list[tuple[Protocluster, Location]]] = defaultdict(list)
    for cluster in clusters:
        by_product[cluster.product].append((cluster, record.extend_location(cluster.core_location, cluster.cutoff)))

    def merge_pair(first: Protocluster, second: Protocluster) -> Protocluster:
        core_location = record.connect_locations([first.core_location, second.core_location])
        surrounding_location = record.extend_location(core_location, first.neighbourhood_range)
        # in cases where the core doesn't span the full region, but the neighbourhood does
        # then the neighbourhood must also cross the origin, so split the difference between
        # the sides and meet (almost) in the middle
        if all([
            len(surrounding_location) == len(record),
            not surrounding_location.crosses_origin(),
            core_location.crosses_origin(),
        ]):
            halfway = (core_location.parts[0].start - core_location.parts[1].end) // 2
            surrounding_location = CompoundLocation([
                FeatureLocation(core_location.parts[0].start - halfway, len(record), 1),
                FeatureLocation(0, core_location.parts[1].end + halfway - 1, 1),
            ])
        return Protocluster(
            core_location=core_location,
            surrounding_location=surrounding_location,
            tool=first.tool,
            cutoff=first.cutoff,
            neighbourhood_range=first.neighbourhood_range,
            product=first.product,
            detection_rule=first.detection_rule,
            product_category=first.product_category,
        )

    results = []
    # merge all overlapping clusters
    for group in by_product.values():
        if len(group) < 2:
            results.extend(group)
            continue
        group.sort(key=lambda x: x[1].start)  # start of each cutoff-extended core
        new_clusters = [group[0]]
        prev_cluster, prev_location = group[0]
        for cluster, location in group[1:]:
            # compare just the core with the extended location, otherwise the cutoff is doubled
            if locations_overlap(cluster.core_location, prev_location):
                merged = merge_pair(prev_cluster, cluster)
                new_clusters[-1] = (merged, record.extend_location(merged.core_location, cluster.cutoff))
                prev_cluster, prev_location = new_clusters[-1]
            else:
                new_clusters.append((cluster, location))
                prev_cluster = cluster
                prev_location = location

        results.extend(new_clusters)

    return sorted(cluster for cluster, _ in results)


def create_rules(rule_files: list[str], signature_names: Set[str],
                 valid_categories: Set[str],
                 multipliers: Multipliers = None,
                 ) -> List[rule_parser.DetectionRule]:
    """ Creates DetectionRule instances from the default rules file

        Args:
            rule_files: a list of paths to files containing cluster rules
            signature_names: the set of all known profile/signature names
            valid_categories: the set of all valid rule categories
            multipliers: distance multipliers to apply to rules

        Returns:
            A list of DetectionRules.
    """
    aliases: dict[str, list[rule_parser.Token]] = {}
    rules: list[rule_parser.DetectionRule] = []
    multipliers = multipliers or Multipliers()
    for rule_file in rule_files:
        with open(rule_file, "r", encoding="utf-8") as ruledata:
            parser = rule_parser.Parser("".join(ruledata.readlines()), signature_names,
                                        valid_categories, rules, existing_aliases=aliases,
                                        multipliers=multipliers)
        aliases.update(parser.aliases)
        rules = parser.rules  # since they include all previous files, replace and don't extend
    return rules


def apply_cluster_rules(record: Record, results_by_id: Dict[str, List[ProfileHit]],
                        rules: Iterable[rule_parser.DetectionRule]
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

    cds_domains_by_cluster_type: Dict[str, Dict[str, Set[str]]] = defaultdict(lambda: defaultdict(set))
    cluster_type_hits: Dict[str, Set[str]] = defaultdict(set)
    for cds_name in cds_with_hits:
        feature = record.get_cds_by_name(cds_name)
        rule_texts = []
        info_by_range: Dict[int, Tuple[Dict[str, CDSFeature], Dict[str, List[HSP]]]] = {}
        for rule in rules:
            if rule.cutoff not in info_by_range:
                location = record.connect_locations([feature.location])
                assert len(location.parts) <= 2, location
                location = _extend_area_location(location, rule.cutoff, record)
                circular_origin = len(record) if len(location.parts) > 1 and record.is_circular() else 0
                nearby = record.get_cds_features_within_location(location, with_overlapping=True)
                nearby_features = {neighbour.get_name(): neighbour for neighbour in nearby}
                nearby_results = {neighbour: results_by_id[neighbour]
                                  for neighbour in nearby_features if neighbour in results_by_id}
                info_by_range[rule.cutoff] = (nearby_features, nearby_results)
            nearby_features, nearby_results = info_by_range[rule.cutoff]
            matching = rule.detect(cds_name, nearby_features, nearby_results, circular_origin=circular_origin)
            if matching.met and matching.matches:
                cds_domains_by_cluster_type[cds_name][rule.name].update(matching.matches)
                rule_texts.append(rule.reconstruct_rule_text())
                cluster_type_hits[rule.name].add(cds_name)
                for other_cds, other_matches in matching.ancillary_hits.items():
                    cluster_type_hits[rule.name].add(other_cds)
                    cds_domains_by_cluster_type[other_cds][rule.name].update(other_matches)
    return cds_domains_by_cluster_type, cluster_type_hits


def find_hmmer_hits(record: Record, sig_by_name: Dict[str, HmmSignature],
                    hmmer_db: str,
                    equivalence_groups: GenericSets,
                    ) -> dict[str, list[ProfileHit]]:
    """ Finds hits for HMMer profiles in the given record

        Arguments:
            record: the record to analyse
            sig_by_name: a dictionary mapping profile name to Signature instance
            hmmer_db: the path to the HMMer database to find hits with
            equivalence_groups: a list of equivalence sets of HMMs

        Returns:
            a dictionary mapping CDS name to a list of ProfileHit instances found
            in that CDS
    """
    results = []
    results_by_id: Dict[str, HSP] = {}
    runresults = run_hmmsearch(hmmer_db, fasta.get_fasta_from_record(record), use_tempfile=True)
    for runresult in runresults:
        acc = runresult.accession.split('.')[0]
        # Store result if it is above cut-off
        for hsp in runresult.hsps:
            if hsp.query_id in sig_by_name:
                sig = sig_by_name[hsp.query_id]
            elif acc in sig_by_name:
                sig = sig_by_name[acc]
            else:
                raise ValueError(f"Failed to find signature for ID {hsp.query_id} / ACC {acc}")
            if hsp.bitscore > sig.cutoff:
                results.append(hsp)
                if hsp.hit_id not in results_by_id:
                    results_by_id[hsp.hit_id] = [hsp]
                else:
                    results_by_id[hsp.hit_id].append(hsp)

    # Filter results by comparing scores of different models (for PKS systems)
    results, results_by_id = filter_results(results, results_by_id, equivalence_groups)

    # Filter multiple results of the same model in one gene
    results, results_by_id = filter_result_multiple(results, results_by_id)

    by_id: Dict[str, List[ProfileHit]] = defaultdict(list)
    for hsp in results:
        by_id[hsp.hit_id].append(HMMerHit.from_hsp(hsp, sig_by_name[hsp.query_id].seed_count))

    return by_id


def build_results(clusters: list[Protocluster], record: Record, tool: str,
                  results_by_id: dict[str, list[HSP]],
                  cds_domains_by_cluster: dict[str, dict[str, set[str]]],
                  annotate_existing_subregions: bool,
                  multipliers: Multipliers,
                  ) -> RuleDetectionResults:
    """ Builds a RuleDetectionResults instance from the provided details

        Arguments:
            clusters: a list of protoclusters found
            record: the record that was analysed
            tool: the name of the tool providing the HMMs (e.g. rule_based_clusters)
            cds_domains_by_cluster: a mapping of CDS ID to
                    a mapping of cluster type string to
                        a set of domains used to determine the cluster
            annotate_existing_subregions: if True, subregions already present in the record
                    will have domains annotated even if no protocluster is found
            multipliers: distance multipliers applied to the rules used

        Returns:
            the created RuleDetectionResults instance
    """

    def get_domains_for_cds(cds: CDSFeature) -> List[SecMetQualifier.Domain]:
        domains = []
        for hit in results_by_id.get(cds.get_name(), []):
            domains.append(SecMetQualifier.Domain(hit.query_id, hit.evalue, hit.bitscore,
                                                  hit.seeds, tool))
        return domains

    all_cds_results: dict[str, CDSResults] = {}

    cds_results_by_cluster = {}
    cdses_with_annotations = set()
    for cluster in clusters:
        cds_results = []
        for cds in record.get_cds_features_within_location(cluster.location):
            domains = get_domains_for_cds(cds)
            if domains:
                def_domains = cds_domains_by_cluster.get(cds.get_name(), {}).get(cluster.product, set())
                cds_result = all_cds_results.get(cds.get_name())
                if not cds_result:
                    cds_result = CDSResults(cds, domains, {})
                    all_cds_results[cds.get_name()] = cds_result
                cds_result.definition_domains[cluster.product] = def_domains

                cds_results.append(cds_result)
                cdses_with_annotations.add(cds)
        cds_results_by_cluster[cluster] = cds_results

    # add detected profile annotations for any existing subregions, if enabled
    cds_results_outside_clusters = []
    if annotate_existing_subregions:
        for subregion in record.get_subregions():
            for cds in subregion.cds_children:
                if cds in cdses_with_annotations:
                    continue
                domains = get_domains_for_cds(cds)
                if domains:
                    cds_results_outside_clusters.append(CDSResults(cds, domains, {}))
                    cdses_with_annotations.add(cds)

    return RuleDetectionResults(cds_results_by_cluster, tool, cds_results_outside_clusters,
                                multipliers)


def detect_protoclusters_and_signatures(record: Record, ruleset: Ruleset,
                                        annotate_existing_subregions: bool = True,
                                        ) -> RuleDetectionResults:
    """ Compares all CDS features in a record with HMM signatures and generates
        Protocluster features based on those hits and the current protocluster detection
        rules.

        Arguments:
            record: the record to analyse
            ruleset: the Ruleset instance defining what profiles exist
            annotate_existing_subregions: if True, subregions already present in the record
                    will have domains annotated even if no protocluster is found

        Returns:
            an instance of RuleDetectionResults
    """
    # if there's no CDS features, don't try to do anything
    if not record.get_cds_features():
        return RuleDetectionResults({}, ruleset.tool, [], ruleset.multipliers)

    # defaults in case of no HMMer profiles
    results_by_id: Dict[str, List[ProfileHit]] = {}

    # get the HMMer profile results, if relevant
    if ruleset.hmm_profiles:
        results_by_id.update(find_hmmer_hits(record, ruleset.hmm_profiles, ruleset.database_file,
                                             ruleset.get_equivalence_groups()))

    # gather dynamic hits and merge them with HMMer results
    dynamic_results = find_dynamic_hits(record, list(ruleset.dynamic_profiles.values()), results_by_id)
    for name, dynamic_hits in dynamic_results.items():
        if name not in results_by_id:
            results_by_id[name] = []
        results_by_id[name].extend(dynamic_hits)

    # Use rules to determine gene clusters
    cds_domains_by_cluster, cluster_type_hits = apply_cluster_rules(record, results_by_id, ruleset.rules)

    # annotate everything in detected protoclusters
    rules_by_name = {rule.name: rule for rule in ruleset.rules}
    clusters = find_protoclusters(record, cluster_type_hits, rules_by_name, results_by_id, cds_domains_by_cluster)
    strip_inferior_domains(cds_domains_by_cluster, rules_by_name)

    return build_results(clusters, record, ruleset.tool, results_by_id, cds_domains_by_cluster,
                         annotate_existing_subregions, ruleset.multipliers)


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


def find_dynamic_hits(record: Record, dynamic_profiles: list[DynamicProfile],
                      hmmer_hits: dict[str, list[ProfileHit]]) -> dict[str, list[DynamicHit]]:
    """ Finds hits for dynamic profiles

    Arguments:
        record: the Record to search
        dynamic_profiles: the dynamic profiles to find hits with
        hmmer_hits: existing HMMer profile hits for use in dynamic profiles

    Returns:
        a dictionary mapping CDS name to list of DynamicHit
    """
    results: Dict[str, List[DynamicHit]] = defaultdict(list)
    for profile in dynamic_profiles:
        for name, hits in profile.find_hits(record, hmmer_hits).items():
            results[name].extend(hits)
    return results
