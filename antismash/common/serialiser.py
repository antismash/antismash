# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

from collections import OrderedDict
import json
import logging

import Bio.Alphabet
import Bio.Alphabet.IUPAC
from Bio.Seq import Seq
from Bio.SeqFeature import ExactPosition, BeforePosition, AfterPosition, \
                           UnknownPosition, FeatureLocation, CompoundLocation, \
                           SeqFeature, Reference
from Bio.SeqRecord import SeqRecord

from antismash.common.module_results import ModuleResults

def write_records(records, results, handle):
    data = []
    for record, result in zip(records, results):
        json_record = record_to_json(record)
        modules = {}
        logging.debug("Record %s has results for modules: %s", record.id,
                      ", ".join([mod for mod, resultv in result.get("modules", {}).items() if resultv]))
        for module, m_results in result.get("modules", {}).items():
            logging.debug("Converting %s results to json", module)
            if not m_results:
                continue
            if isinstance(m_results, ModuleResults):
                modules[module] = m_results.to_json()
            elif isinstance(m_results, dict): # TODO :preferably no branching here
                logging.critical("module results was a dict, not ModuleResults")
                modules[module] = m_results
            else:
                raise TypeError("Module results for module %s are of invalid type: %s" % (module, type(m_results)))
        json_record["modules"] = modules
        data.append(json_record)
    # only wipe existing data if we have a valid file afterwards
    new_contents = json.dumps(data)
    if isinstance(handle, str):
        handle = open(handle, "w")
    handle.write(new_contents)

def read_records(handle):
    if isinstance(handle, str):
        handle = open(handle, "r")
    contents = handle.read()
    if not contents:
        raise ValueError("Results file contains no information")
    data = json.loads(contents, object_pairs_hook=OrderedDict)
    return list(map(record_from_json, data))

def record_to_json(record):
    def annotations_to_json(annotations):
        res = dict(annotations)
        res["references"] = []
        for reference in annotations["references"]:
            ref = dict(reference.__dict__)
            ref["location"] = [location_to_json(loc) for loc in ref["location"]]
            res["references"].append(ref)
        return res

    result = {}
    result["id"] = record.id
    result["seq"] = sequence_to_json(record.seq)
    result["features"] = list(map(feature_to_json, record.features))
    result["name"] = record.name
    result["description"] = record.description
    result["dbxrefs"] = record.dbxrefs
    result["annotations"] = annotations_to_json(record.annotations)
    result["letter_annotations"] = record.letter_annotations
    return result

def record_from_json(data):
    if isinstance(data, str):
        data = json.loads(data)

    def rebuild_references(annotations):
        bases = annotations["references"]
        refs = []
        for ref in bases:
            new_reference = Reference()
            new_reference.__dict__ = ref
            new_reference.location = [location_from_json(loc) for loc in ref["location"]]
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


def sequence_to_json(sequence):
    return {"data" : str(sequence),
            "alphabet" : str(sequence.alphabet).rsplit('()')[0]} # DNA() -> DNA

def sequence_from_json(data):
    if isinstance(data, str):
        data = json.loads(data)
    alphabet = data["alphabet"]
    if "IUPAC" in alphabet:
        alphabet_class = getattr(Bio.Alphabet.IUPAC, alphabet)
    else:
        alphabet_class = getattr(Bio.Alphabet, alphabet)
    return Seq(data["data"], alphabet=alphabet_class())


def feature_to_json(feature):
    return json.dumps({"location" : location_to_json(feature.location),
                       "type" : feature.type,
                       "id" : feature.id,
                       "qualifiers" : feature.qualifiers
                      })

def feature_from_json(data):
    if isinstance(data, str):
        data = json.loads(data, object_pairs_hook=OrderedDict)
    return SeqFeature(location=location_from_json(data["location"]),
                      type=data["type"],
                      id=data["id"],
                      qualifiers=data["qualifiers"])


def location_to_json(location):
    return str(location)

def location_from_json(data):
    """
        Converts from json representation, e.g. [<1:6](-), to FeatureLocation
        or CompoundLocation
    """
    def parse_position(string):
        if string[0] == '<':
            return BeforePosition(int(string[1:]))
        if string[0] == '>':
            return AfterPosition(int(string[1:]))
        if string == "UnknownPosition()":
            return UnknownPosition()
        return ExactPosition(int(string))

    def parse_single_location(string):
        start = parse_position(string[1:].split(':', 1)[0]) # [<1:6](-) -> <1
        end = parse_position(string.split(':', 1)[1].split(']', 1)[0]) # [<1:6](-) -> 6

        strand = string[-2] # [<1:6](-) -> -
        if strand == '-':
            strand = -1
        elif strand == '+':
            strand = 1
        elif strand == '?':
            strand = 0
        elif '(' not in string:
            strand = None
        else:
            raise ValueError("Cannot identify strand in location: %s", string)


        return FeatureLocation(start, end, strand=strand)

    assert isinstance(data, str)

    if '{' not in data:
        return parse_single_location(data)

    # otherwise it's a compound location
    # join{[1:6](+), [10:16](+)} -> ("join", "[1:6](+), [10:16](+)")
    operator, locations = data[:-1].split('{', 1)

    locations = [parse_single_location(part) for part in locations.split(', ')]
    return CompoundLocation(locations, operator=operator)
