# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" A collection of functions and classes for serialising the result of
    running antismash analyses.
"""

from collections import OrderedDict
import json
import logging
from typing import Any, Dict, IO, List, Union

import Bio.Alphabet
import Bio.Alphabet.IUPAC
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, Reference
from Bio.SeqRecord import SeqRecord

from antismash.common.module_results import ModuleResults
from antismash.common.secmet import Record
from antismash.common.secmet.locations import location_from_string


class AntismashResults:
    """ A single repository of all results of an antismash run, including input
        filename, records and individual module results
    """
    SCHEMA_VERSION = 1

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
            handle = open(handle, "r")
        data = json.loads(handle.read())
        if data.get("schema", 1) != AntismashResults.SCHEMA_VERSION:
            raise ValueError("schema mismatch in previous results: expected %s, found %s" % (
                                AntismashResults.SCHEMA_VERSION, data.get("schema")))
        version = data["version"]
        input_file = data["input_file"]
        taxon = data.get("taxon", "bacteria")
        records = [Record.from_biopython(record_from_json(rec), taxon) for rec in data["records"]]
        results = [rec["modules"] for rec in data["records"]]
        return AntismashResults(input_file, records, results, version, taxon=taxon)

    def to_json(self) -> Dict[str, Any]:
        """ Constructs a JSON representation of the instance """
        res = OrderedDict()  # type: Dict[str, Any]
        res["version"] = self.version
        res["input_file"] = self.input_file
        biopython = [rec.to_biopython() for rec in self.records]
        res["records"] = dump_records(biopython, self.results)
        res["timings"] = self.timings_by_record
        res["taxon"] = self.taxon
        res["schema"] = self.SCHEMA_VERSION
        return res

    def write_to_file(self, handle: Union[str, IO]) -> None:
        """ Writes a JSON representation of the instance to the given filename
            or handle
        """
        if isinstance(handle, str):
            handle = open(handle, "w")
        handle.write(json.dumps(self.to_json()))


def dump_records(records: List[SeqRecord], results: List[Dict[str, Union[Dict[str, Any], ModuleResults]]],
                 handle: Union[str, IO] = None) -> List[Dict[str, Any]]:
    """ Converts a list of records and a list of results to a JSON object.

        Arguments:
            records: a list of records to convert
            results: a matching list of results to convert
            handle: a filename or file-like object to write the resulting JSON
                    object to

        Returns:
            the JSON object created
    """
    data = []
    assert isinstance(results, list)
    for record, result in zip(records, results):
        json_record = record_to_json(record)
        modules = OrderedDict()  # type: Dict[str, Dict]
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
                raise TypeError("Module results for module %s are of invalid type: %s" % (module, type(m_results)))
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
        handle = open(handle, "w")
    handle.write(new_contents)
    return data


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

    result = OrderedDict()  # type: Dict[str, Any]
    result["id"] = record.id
    result["seq"] = sequence_to_json(record.seq)
    result["features"] = list(map(feature_to_json, record.features))
    result["name"] = record.name
    result["description"] = record.description
    result["dbxrefs"] = record.dbxrefs
    result["annotations"] = annotations_to_json(record.annotations)
    result["letter_annotations"] = record.letter_annotations
    return result


def record_from_json(data: Union[str, Dict]) -> SeqRecord:
    """ Rebuilds a SeqRecord from JSON """
    if isinstance(data, str):
        data = json.loads(data)
    assert isinstance(data, dict)

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

    return SeqRecord(sequence_from_json(data["seq"]),
                     id=data["id"],
                     name=data["name"],
                     description=data["description"],
                     dbxrefs=data["dbxrefs"],
                     features=list(map(feature_from_json, data["features"])),
                     annotations=rebuild_references(data["annotations"]),
                     letter_annotations=data["letter_annotations"])


def sequence_to_json(sequence: Seq) -> Dict[str, str]:
    """ Constructs a JSON object that represents a Seq sequence """
    return {"data": str(sequence),
            "alphabet": str(sequence.alphabet).rsplit('()')[0]}  # DNA() -> DNA


def sequence_from_json(data: Union[str, Dict]) -> Seq:
    """ Reconstructs a Seq sequence from JSON """
    if isinstance(data, str):
        data = json.loads(data)
    assert isinstance(data, dict)
    alphabet = data["alphabet"]
    if "IUPAC" in alphabet:
        alphabet_class = getattr(Bio.Alphabet.IUPAC, alphabet)
    else:
        alphabet_class = getattr(Bio.Alphabet, alphabet)
    return Seq(data["data"], alphabet=alphabet_class())


def feature_to_json(feature: SeqFeature) -> Dict[str, Any]:
    """ Creates a JSON representation of a SeqFeature """
    return {"location": str(feature.location),
            "type": feature.type,
            "id": feature.id,
            "qualifiers": feature.qualifiers}


def feature_from_json(data: Union[str, Dict]) -> SeqFeature:
    """ Converts a JSON representation of a feature into a SeqFeature """
    if isinstance(data, str):
        data = json.loads(data, object_pairs_hook=OrderedDict)
    assert isinstance(data, dict)
    return SeqFeature(location=location_from_string(data["location"]),
                      type=data["type"],
                      id=data["id"],
                      qualifiers=data["qualifiers"])
