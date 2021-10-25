# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Generic sideloading """

from typing import List, Optional

from antismash.common import path
from antismash.common.secmet import Record

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
                                   manual: Optional[SideloadSimple]) -> SideloadedResults:
    """ Loads generic subregion/protocluster annotations from JSON files.

        Arguments:
            annotation_file: the paths to the JSON files containing annotations
            record_id: a record id to restrict annotation generation to
            manual: an optional manually specified area

        Returns:
            a GenericAnnotations instance containing all information from the
            annotation files
    """
    subregions: List[SubRegionAnnotation] = []
    protoclusters:List[ProtoclusterAnnotation] = []
    for annotations_file in annotation_files:
        raw = load_validated_json(annotations_file, _SCHEMA_FILE)
        tool = Tool.from_json(raw["tool"])
        for json_record in raw["records"]:
            name = json_record["name"]
            if name != record.id:
                continue
            for area in json_record.get("subregions", []):
                subregions.append(SubRegionAnnotation.from_schema_json(area, tool))
            for area in json_record.get("protoclusters", []):
                protoclusters.append(ProtoclusterAnnotation.from_schema_json(area, tool))

    if manual and manual.accession == record.id:
        tool = Tool("manual", "N/A", "command line argument", {})
        subregion = SubRegionAnnotation(manual.start, manual.end, "", tool, {})
        subregions.append(subregion)

    return SideloadedResults(record.id, subregions, protoclusters)
