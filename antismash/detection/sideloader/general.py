# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Generic sideloading """

import itertools
from typing import List, Optional

from antismash.common import path
from antismash.common.errors import AntismashInputError
from antismash.common.secmet import FeatureLocation, Record

from .data_structures import (
    ProtoclusterAnnotation,
    SideloadSimple,
    SideloadedResults,
    SubRegionAnnotation,
    Tool,
)
from .loader import load_validated_json

_SCHEMA_FILE = path.get_full_path(__file__, "schemas", "general", "schema.json")


def load_single_record_annotations(annotation_files: List[str], record: Record,
                                   manual: Optional[SideloadSimple],
                                   cds_markers: Optional[List[str]] = None,
                                   cds_marker_padding: int = 20000) -> SideloadedResults:
    """ Loads generic subregion/protocluster annotations from JSON files.

        Arguments:
            annotation_file: the paths to the JSON files containing annotations
            record_id: a record id to restrict annotation generation to
            manual: an optional manually specified area
            cds_markers: a list of CDS locus tags to manually generate subregions around
            cds_marker_padding: the size to include, in nucleotides, on each side of any CDS marker

        Returns:
            a GenericAnnotations instance containing all information from the
            annotation files
    """
    subregions: List[SubRegionAnnotation] = []
    protoclusters: List[ProtoclusterAnnotation] = []
    for annotations_file in annotation_files:
        raw = load_validated_json(annotations_file, _SCHEMA_FILE)
        tool = Tool.from_json(raw["tool"])
        for json_record in raw["records"]:
            name = json_record["name"]
            if not record.has_name(name):
                continue
            for area in json_record.get("subregions", []):
                subregions.append(SubRegionAnnotation.from_schema_json(area, tool))
            for area in json_record.get("protoclusters", []):
                protoclusters.append(ProtoclusterAnnotation.from_schema_json(area, tool))

    tool = Tool("manual", "N/A", "command line argument", {})
    if manual and record.has_name(manual.accession):
        subregion = SubRegionAnnotation(manual.start, min(manual.end, len(record.seq)), "", tool, {})
        subregions.append(subregion)

    if cds_markers:
        for name in cds_markers:
            try:
                cds = record.get_cds_by_name(name)
            except KeyError:
                continue
            start = max(cds.location.start - cds_marker_padding, 0)
            end = min(cds.location.end + cds_marker_padding, len(record.seq))
            subregion = SubRegionAnnotation(start, end, name, tool, {})
            subregions.append(subregion)

    for area in itertools.chain(protoclusters, subregions):
        location = FeatureLocation(area.start, area.end)
        if not record.get_cds_features_within_location(location):
            raise AntismashInputError(f"sideloaded area contains no complete CDS features in {record.id}: {area}")

    return SideloadedResults(record.id, subregions, protoclusters)
