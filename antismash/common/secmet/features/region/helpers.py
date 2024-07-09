# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Contains a feature covering a region, along with helper/utility functions for
    manipulating those features.
"""

from copy import deepcopy
from dataclasses import dataclass
from typing import Any, IO
import warnings

from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord
from helperlibs.bio import seqio

from ..abstract import AbstractRegion
from ..candidate_cluster import CandidateCluster
from ...locations import (
    FeatureLocation,
    Location,
    build_location_from_others,
    location_from_string,
)
from ..protocluster import Protocluster
from ..subregion import SubRegion


@dataclass(kw_only=True)
class RegionData:
    """ A circular-import safe container with a subset of full region information """
    start: int
    end: int
    candidate_clusters: tuple[CandidateCluster, ...]
    subregions: tuple[SubRegion, ...]


def _build_annotations(region: RegionData, original_annotations: dict[str, Any]) -> dict[str, Any]:
    """ Builds a new set of annotations from the given annotations and the relevant
        region data.

        Arguments:
            region: the details of the region being annotated
            original_annotations: the annotations of the full-genome record, as in a SeqRecord

        Returns:
            the annotations in a SeqRecord-compatible dictionary
    """
    # since modifications will be necessary, make and use a copy
    annotations = deepcopy(original_annotations)

    # while there ought to be existing antiSMASH structured comments, allow for them to be missing
    annotations.setdefault("structured_comment", {})
    annotations["structured_comment"].setdefault("antiSMASH-Data", {})
    antismash_comment = annotations["structured_comment"]["antiSMASH-Data"]

    # any additions here need to be added in the order they should appear

    antismash_comment["NOTE"] = "This is a single region extracted from a larger record!"
    antismash_comment["Orig. start"] = str(region.start)
    antismash_comment["Orig. end"] = str(region.end)

    return annotations


def _build_base_record(region: RegionData, record: SeqRecord) -> SeqRecord:
    """ Extracts a minimally adjusted new record for a region from a larger record.

        Arguments:
            region: the data of the region being extracted
            record: the parent record from which to extract a record for the region

        Returns:
            a new record covering only the given region
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        region_record = record[region.start:region.end]

    return region_record


def _adjust_motif(feature: SeqFeature, region: RegionData) -> None:
    """ Adjusts a prepeptide/CDS motif to refer to coordinates relative to the
        region start.

        Arguments:
            feature: the feature to adjust
            region: the data of the relevant region
    """
    for qual in ["leader_location", "tail_location"]:
        if qual not in feature.qualifiers:
            continue
        loc = location_from_string(feature.qualifiers[qual][0])
        parts = []
        for part in loc.parts:
            new_start = part.start - region.start
            new_end = part.end - region.start
            parts.append(FeatureLocation(new_start, new_end, part.strand))
        feature.qualifiers[qual] = [str(build_location_from_others(parts))]


def _adjust_protocluster(feature: SeqFeature, protocluster: Protocluster,
                         region: RegionData, new_number: int,
                         ) -> None:
    """ Adjusts a protocluster's references to be relative to the region start.

        Arguments:
            feature: the feature to adjust
            protocluster: the secmet protocluster feature the feature was built from
            region: the data of the relevant region
            new_number: the new protocluster number
    """
    # update core location qualifier first, if it's not the core feature
    if feature.type == Protocluster.FEATURE_TYPE:
        location = protocluster.core_location
        new_location = location.clone_with_offset(-region.start)
        feature.qualifiers["core_location"] = [str(new_location)]
    # then protocluster number
    feature.qualifiers["protocluster_number"] = [str(new_number)]


def _adjust_features(region: RegionData, region_record: SeqRecord) -> None:
    """ Adjusts any relevant features to be relative to the region start.

        Arguments:
            region: the data of the relevant region
            region_record: the record extracted for the given region
            record: the parent record of the region
    """
    protoclusters_by_original_number = {}
    for candidate in region.candidate_clusters:
        for protocluster in candidate.protoclusters:
            protoclusters_by_original_number[protocluster.get_protocluster_number()] = protocluster

    if region.candidate_clusters:
        first_candidate_cluster = min(cc.get_candidate_cluster_number() for cc in region.candidate_clusters)
        first_cluster = min(cluster.get_protocluster_number() for cluster in protoclusters_by_original_number.values())
    else:
        first_candidate_cluster = 0
        first_cluster = 0
    first_subregion = min(sub.get_subregion_number() for sub in region.subregions) if region.subregions else 0

    for feature in region_record.features:
        if feature.type == AbstractRegion.FEATURE_TYPE:
            candidates = feature.qualifiers.get("candidate_cluster_numbers")
            if not candidates:
                continue
            candidates = [str(int(num) - first_candidate_cluster + 1) for num in candidates]
            feature.qualifiers["candidate_cluster_numbers"] = candidates
        elif feature.type == CandidateCluster.FEATURE_TYPE:
            new = str(int(feature.qualifiers["candidate_cluster_number"][0]) - first_candidate_cluster + 1)
            feature.qualifiers["candidate_cluster_number"] = [new]
            new_clusters = [str(int(num) - first_cluster + 1) for num in feature.qualifiers["protoclusters"]]
            feature.qualifiers["protoclusters"] = new_clusters
        elif feature.type in [Protocluster.FEATURE_TYPE, "proto_core"]:
            original_number = int(feature.qualifiers["protocluster_number"][0])
            new_number = original_number - first_cluster + 1
            _adjust_protocluster(feature, protoclusters_by_original_number[original_number],
                                 region, new_number)
        elif feature.type == "subregion":
            new = str(int(feature.qualifiers["subregion_number"][0]) - first_subregion + 1)
            feature.qualifiers["subregion_number"] = [new]
        elif feature.type == "CDS_motif":
            _adjust_motif(feature, region)


def write_to_genbank(region: RegionData, record: SeqRecord, handle: IO) -> None:
    """ Writes a genbank file containing only the information contained
        within the region, linearising the new record if the region crossed the origin.

        Arguments:
            region: the data of the region
            record: the parent record of the region
            handle: the file handle in which the resulting genbank will be written
    """
    assert isinstance(record, SeqRecord), type(record)
    # some location modifications may be necessary in cross-origin regions,
    # which means the originals must be kept and changes reverted
    original_locations: dict[int, Location] = {id(feature): feature.location for feature in record.features}

    region_record = _build_base_record(region, record)
    _adjust_features(region, region_record)

    region_record.annotations = _build_annotations(region, record.annotations)

    seqio.write([region_record], handle, "genbank")

    # undo any location modifications
    for feature in record.features:
        feature.location = original_locations[id(feature)]
