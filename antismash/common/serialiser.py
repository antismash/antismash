# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" A collection of functions and classes for serialising the result of
    running antismash analyses.
"""

import bz2
from collections import defaultdict, OrderedDict
import logging
import os
from typing import Any, Dict, IO, List, Union

from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, Reference
from Bio.SeqRecord import SeqRecord

from antismash.common import json
from antismash.common.module_results import ModuleResults
from antismash.common.secmet import Record
from antismash.common.secmet.locations import location_from_string

# Schema version changes:
# 1: initial state
# 2: added field: records.areas
# 3: added field to records.areas.protoclusters objects: category
# 4: added field to records: gc_content


class AntismashResults:
    """ A single repository of all results of an antismash run, including input
        filename, records and individual module results
    """
    SCHEMA_VERSION = 3

    # the key must accept all the schemas in values as valid
    # this will typically only be useful for backwards compatibility
    # typically different schema versions will have different detection rules,
    # so this won't really be relevant for reusing results between runs
    COMPATIBLE_SCHEMAS = defaultdict(set, {
        2: {1},  # records.areas isn't used
        3: {2, 1},  # a subset of records.areas changed, but still isn't used
        4: {3, 2, 1},  # gc_content will be missing, be can always be recalculated
    })

    def __init__(self, input_file: str, records: List[Record],
                 results: List[Dict[str, Union[ModuleResults, Dict[str, Any]]]],
                 version: str, timings: Dict[str, Dict[str, float]] = None,
                 taxon: str = "bacteria") -> None:
        self.input_file = input_file
        self.records = records
        self.results = results
        self.version = version
        self.timings_by_record = timings or {}  # {record_id : {module name: time}}
        self.taxon = taxon

    @staticmethod
    def from_file(handle: Union[str, IO]) -> "AntismashResults":
        """ Regenerates an instance of AntismashResults from JSON representation
            in a file
        """
        if isinstance(handle, str):
            _, ext = os.path.splitext(handle)
            if ext == ".bz2":
                handle = bz2.open(handle, "rt", encoding="utf-8")
            else:
                handle = open(handle, "r", encoding="utf-8")  # pylint: disable=consider-using-with
        try:
            data = json.loads(handle.read())
        except json.JSONDecodeError:
            raise ValueError(f"Cannot load results to reuse from {handle.name}, "
                             "is it an antiSMASH result JSON file?")
        schema = data.get("schema", 1)
        current = AntismashResults.SCHEMA_VERSION
        if schema != current and schema not in AntismashResults.COMPATIBLE_SCHEMAS[current]:
            raise ValueError(
                f"schema mismatch in previous results: expected {AntismashResults.SCHEMA_VERSION}"
                f", found {data.get('schema')}"
            )
        version = data["version"]
        input_file = data["input_file"]
        taxon = data.get("taxon", "bacteria")
        records = [record_from_json(rec, taxon) for rec in data["records"]]
        for record, rec_json in zip(records, data["records"]):
            if "original_id" in rec_json:
                record.original_id = rec_json["original_id"]
        results = [rec["modules"] for rec in data["records"]]
        return AntismashResults(input_file, records, results, version, taxon=taxon)

    def to_json(self) -> Dict[str, Any]:
        """ Constructs a JSON representation of the instance """
        res: Dict[str, Any] = OrderedDict()
        res["version"] = self.version
        res["input_file"] = self.input_file
        res["records"] = dump_records(self.results, self.records)
        res["timings"] = self.timings_by_record
        res["taxon"] = self.taxon
        res["schema"] = self.SCHEMA_VERSION
        return res

    def write_to_file(self, handle: Union[str, IO]) -> None:
        """ Writes a JSON representation of the instance to the given filename
            or handle
        """
        # do the JSON conversions before overwriting the previous file
        try:
            converted = json.dumps(self.to_json())
        except TypeError as error:
            message = f"Failed to convert JSON results: {error}"
            logging.error(message)
            raise TypeError(message) from error
        if isinstance(handle, str):
            handle = open(handle, "w", encoding="utf-8")  # pylint: disable=consider-using-with
        handle.write(converted)


def dump_records(results: List[Dict[str, Union[Dict[str, Any], ModuleResults]]],
                 secmet_records: List[Record], handle: Union[str, IO] = None) -> List[Dict[str, Any]]:
    """ Converts a list of records and a list of results to a JSON object.

        Arguments:
            results: a matching list of results to convert
            secmet_records: a matching list of Records for easier data access
            handle: a filename or file-like object to write the resulting JSON
                    object to

        Returns:
            the JSON object created
    """
    data = []
    assert isinstance(results, list)
    for i, secmet in enumerate(secmet_records):
        result = results[i]
        record = secmet.to_biopython()
        json_record = record_to_json(record)
        json_record["areas"] = gather_record_areas(secmet)
        if secmet.original_id:
            json_record["original_id"] = secmet.original_id
        json_record["gc_content"] = secmet.get_gc_content()
        modules: Dict[str, Dict] = OrderedDict()
        if result:
            logging.debug("Record %s has results for modules: %s", record.id,
                          ", ".join([mod.rsplit('.', 1)[-1] for mod, resultv in result.items() if resultv]))
        for module, m_results in result.items():
            logging.debug("Converting %s results to json", module)
            if m_results is None:
                logging.debug("%s results didn't exist", module)
                continue
            if isinstance(m_results, ModuleResults):
                modules[module] = m_results.to_json()
            else:
                raise TypeError(f"Module results for module {module} are of invalid type: {type(m_results)}")
        json_record["modules"] = modules
        data.append(json_record)

    if handle is None:
        return data

    # only wipe existing data if we have a valid file afterwards
    try:
        new_contents = json.dumps(data)
    except TypeError:
        logging.error("Error converting json data: %s", data)
        raise
    if isinstance(handle, str):
        handle = open(handle, "w", encoding="utf-8")  # pylint: disable=consider-using-with
    handle.write(new_contents)
    return data


def gather_record_areas(record: Record) -> List[Dict[str, Any]]:
    """ Constructs a JSON friendly representation of areas in Record

        Arguments:
            record: the Record instance to gather areas from

        Returns:
            a list of JSON appropriate dictionaries, one for each region
    """
    regions = []
    for region in record.get_regions():
        region_json = {
            "start": region.location.start,
            "end": region.location.end,
            "products": region.products,
            "protoclusters": {},
            "candidates": [],
            "subregions": [],
        }
        regions.append(region_json)

        protoclusters_by_obj = {proto: i for i, proto in enumerate(region.get_unique_protoclusters())}
        for proto, i in protoclusters_by_obj.items():
            region_json["protoclusters"][i] = {
                "category": proto.product_category,
                "start": proto.location.start,
                "end": proto.location.end,
                "core_start": proto.core_location.start,
                "core_end": proto.core_location.end,
                "product": proto.product,
                "tool": proto.tool,
            }
        for candidate in region.candidate_clusters:
            region_json["candidates"].append({
                "start": candidate.location.start,
                "end": candidate.location.end,
                "kind": str(candidate.kind),
                "protoclusters": [protoclusters_by_obj[proto] for proto in candidate.protoclusters],
            })
        for subregion in region.subregions:
            region_json["subregions"].append({
                "start": subregion.location.start,
                "end": subregion.location.end,
                "label": subregion.label,
                "tool": subregion.tool,
            })
    return regions


def record_to_json(record: SeqRecord) -> Dict[str, Any]:
    """ Constructs a JSON object representing a SeqRecord """
    def annotations_to_json(annotations: Dict) -> Dict[str, Any]:
        """ Converts the 'annotations' member of a SeqRecord """
        res = dict(annotations)
        res["references"] = []
        for reference in annotations.get("references", []):
            ref = dict(reference.__dict__)
            ref["location"] = [str(loc) for loc in ref["location"]]
            res["references"].append(ref)
        return res

    result: Dict[str, Any] = OrderedDict()
    result["id"] = record.id
    result["seq"] = sequence_to_json(record.seq)
    result["features"] = list(map(feature_to_json, record.features))
    result["name"] = record.name
    result["description"] = record.description
    result["dbxrefs"] = record.dbxrefs
    result["annotations"] = annotations_to_json(record.annotations)
    result["letter_annotations"] = record.letter_annotations
    return result


def record_from_json(data: Union[str, Dict], taxon: str) -> Record:
    """ Rebuilds a SeqRecord from JSON """
    if isinstance(data, str):
        data = json.loads(data)
    assert isinstance(data, dict)

    # build extra optional values for records that may or may not exist in the data
    # only being present when set
    kwargs = {}

    gc_content = data.get("gc_content")
    if gc_content is not None:
        kwargs["gc_content"] = gc_content

    def rebuild_references(annotations: Dict) -> Dict[str, List[Reference]]:
        """ Rebuilds the SeqRecord 'references' annotation from JSON """
        bases = annotations["references"]
        refs = []
        for ref in bases:
            new_reference = Reference()
            new_reference.__dict__ = ref
            new_reference.location = [location_from_string(loc) for loc in ref["location"]]
            refs.append(new_reference)
        annotations["references"] = refs
        return annotations
    seqb = SeqRecord(sequence_from_json(data["seq"]),
                     id=data["id"],
                     name=data["name"],
                     description=data["description"],
                     dbxrefs=data["dbxrefs"],
                     features=list(map(feature_from_json, data["features"])),
                     annotations=rebuild_references(data["annotations"]),
                     letter_annotations=data["letter_annotations"])
    return Record.from_biopython(seqb, taxon, **kwargs)


def sequence_to_json(sequence: Seq) -> Dict[str, str]:
    """ Constructs a JSON object that represents a Seq sequence """
    return {
        "data": str(sequence),
        "alphabet": "DNA",  # for compatibility, removable on schema version change
    }


def sequence_from_json(data: Union[str, Dict]) -> Seq:
    """ Reconstructs a Seq sequence from JSON """
    if isinstance(data, str):
        data = json.loads(data)
    assert isinstance(data, dict)
    return Seq(data["data"])


def feature_to_json(feature: SeqFeature) -> Dict[str, Any]:
    """ Creates a JSON representation of a SeqFeature """
    return {"location": str(feature.location),
            "type": feature.type,
            "qualifiers": feature.qualifiers}


def feature_from_json(data: Union[str, Dict]) -> SeqFeature:
    """ Converts a JSON representation of a feature into a SeqFeature """
    if isinstance(data, str):
        data = json.loads(data)
    assert isinstance(data, dict)
    return SeqFeature(location=location_from_string(data["location"]),
                      type=data["type"],
                      qualifiers=data["qualifiers"])
