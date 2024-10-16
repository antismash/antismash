# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""Atropopeptides detection module

"""

from collections import defaultdict
import logging
from typing import Any, Iterable, Optional

from antismash.common import (
    all_orfs,
    comparippson,
    module_results,
)
from antismash.common.secmet import Record, CDSFeature, Prepeptide
from antismash.common.secmet.errors import SecmetInvalidInputError
from antismash.common.secmet.locations import location_from_string
from antismash.common.secmet.qualifiers.gene_functions import GeneFunction
from antismash.config import get_config
from antismash.detection.hmm_detection.dynamic_profiles.atropopeptide_p450 import (
    MAX_DIST,
    MIN_LEN,
    MOTIF,
)


class AtropoResults(module_results.ModuleResults):
    schema_version = 1

    def __init__(self, record_id: str, comparippson_results: comparippson.MultiDBResults | None = None) -> None:
        super().__init__(record_id)
        # keep new CDS features
        self._new_cds_features: set[CDSFeature] = set()
        # keep new CDSMotifs by the gene they match to
        # e.g. self.motifs_by_locus[gene_locus] = [motif1, motif2..]
        self.motifs_by_locus: dict[str, list[Prepeptide]] = defaultdict(list)
        # keep clusters and which genes in them had precursor hits
        # e.g. self.clusters[cluster_number] = {gene1_locus, gene2_locus}
        self.clusters: dict[int, set[str]] = defaultdict(set)
        self.comparippson_results: Optional[comparippson.MultiDBResults] = comparippson_results

    def to_json(self) -> dict[str, Any]:
        cds_features = [(str(feature.location),
                         feature.get_name()) for feature in self._new_cds_features]
        motifs = {}
        for locus, locus_motifs in self.motifs_by_locus.items():
            motifs[locus] = [motif.to_json() for motif in locus_motifs]
        return {
            "record_id": self.record_id,
            "schema_version": AtropoResults.schema_version,
            "motifs": motifs,
            "new_cds_features": cds_features,
            "protoclusters": {key: list(val) for key, val in self.clusters.items()},
            "comparippson_results": self.comparippson_results.to_json() if self.comparippson_results else None,
        }

    @staticmethod
    def from_json(json: dict[str, Any], record: Record) -> Optional["AtropoResults"]:
        if json.get("schema_version") != AtropoResults.schema_version:
            logging.warning("Discarding Atropopeptide results, schema version mismatch")
            return None
        comp_results = comparippson.MultiDBResults.from_json(json["comparippson_results"])
        results = AtropoResults(json["record_id"], comparippson_results=comp_results)
        for locus, motifs in json["motifs"].items():
            for motif in motifs:
                results.motifs_by_locus[locus].append(Prepeptide.from_json(motif))
        results.clusters = {int(key): set(val) for key, val in json["protoclusters"].items()}
        for location, name in json["new_cds_features"]:
            loc = location_from_string(location)
            cds = all_orfs.create_feature_from_location(record, loc, label=name)
            results.add_cds(cds)
        return results

    def add_to_record(self, record: Record) -> None:
        existing = record.get_cds_name_mapping()
        for feature in self._new_cds_features:
            # since a precursor may be found by other RiPP modules
            if feature.get_name() not in existing:
                record.add_cds_feature(feature)

        for motifs in self.motifs_by_locus.values():
            for motif in motifs:
                record.add_cds_motif(motif)

    def add_cds(self, cds: CDSFeature) -> None:
        """ Add a newly found CDS feature that will be added to the record """
        # if already added by another protocluster, don't double up
        for existing in self._new_cds_features:
            if cds.get_name() == existing.get_name():
                return
        self._new_cds_features.add(cds)

def specific_analysis(record: Record) -> AtropoResults:
    """ Runs the full atropopeptide analysis over the given record

        Arguments:
            record: the Record instance to analyse

        Returns:
            A populated LassoResults object
    """
    results = AtropoResults(record.id)
    motif_count = 0
    existing_motifs: dict[str, set[str]] = defaultdict(set)
    found_cdses: set[str] = set()
    for cluster in record.get_protoclusters():
        if not cluster.product.startswith('atropopeptide'):
            continue

        precursor_candidates = list(cluster.cds_children)
        # Find candidate ORFs that are not yet annotated
        extra_orfs = all_orfs.find_all_orfs(record, cluster, min_length=MIN_LEN)
        for i, orf in enumerate(extra_orfs):
            new_orf = all_orfs.get_trimmed_orf(orf, record, min_length=MIN_LEN, label=orf.get_name())
            if new_orf:
                extra_orfs[i] = new_orf
        precursor_candidates.extend(extra_orfs)

        for candidate in precursor_candidates:
            if not candidate.sec_met or not set(candidate.sec_met.domain_ids).intersection({"atropopeptide_p450"}):
                continue

            neighbours = find_neighbours_in_range(candidate, precursor_candidates)
            for hit in neighbours:
                if hit.get_name() in found_cdses:
                    continue
                found = MOTIF.search(hit.translation)
                if found is None:
                    continue

                # leader ends 2 AA after the motif
                leader_end = found.end() + 2
                if len(hit.translation) <= leader_end:
                    logging.debug("Skipping too small atropopeptide precursor %s (%s)", hit.get_name(), hit.translation)
                    continue

                found_cdses.add(hit.get_name())

                # regardless of whether it's a duplication or not, track via protocluster
                results.clusters[cluster.get_protocluster_number()].add(candidate.get_name())

                # but if it is in multiple protoclusters, don't add it more than once
                if hit.get_name() in existing_motifs[candidate.get_name()]:
                    continue
                hit.gene_functions.add(GeneFunction.ADDITIONAL, "atropopeptides",
                             "predicted atropopeptide")


                existing_motifs[candidate.get_name()].add(hit.get_name())
                results.motifs_by_locus[candidate.get_name()].append(Prepeptide(
                    location=hit.location,
                    peptide_class="atropopeptide",
                    locus_tag=hit.get_name(),
                    leader=hit.translation[:leader_end],
                    core=hit.translation[leader_end:],
                    tool="antiSMASH",
                ))
                motif_count += 1
                # track new CDSFeatures if found with all_orfs
                if hit.region is None:
                    try:
                        results.add_cds(hit)
                    except SecmetInvalidInputError:
                        pass

    cores = {}
    for motifs in results.motifs_by_locus.values():
        for motif in motifs:
            cores[motif.get_name()] = motif.core
    results.comparippson_results = comparippson.compare_precursor_cores(cores, get_config())

    logging.debug("Atropopeptide module marked %d motifs", motif_count)
    return results


def find_neighbours_in_range(center: CDSFeature,
                             candidates: Iterable[CDSFeature]) -> list[CDSFeature]:
    """ Restrict a set of genes to those within precursor range of a central
        gene.

        Arguments:
            center: the gene to find the neighbours of
            candidates: the genes to filter by range

        Returns:
            a list of genes within range, with the same ordering as the input
    """
    neighbours = []
    for candidate in candidates:
        if candidate == center:
            continue
        if candidate < center:
            if center.location.start - candidate.location.start <= MAX_DIST:
                neighbours.append(candidate)
        else:
            if candidate.location.end - center.location.end <= MAX_DIST:
                neighbours.append(candidate)
    return neighbours



