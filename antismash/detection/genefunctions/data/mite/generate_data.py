#!/usr/bin/env python3
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Downloads a MITE release and extracts relevant data for use in antiSMASH """


import argparse
import hashlib
import json
import os
import sys
from typing import Any
from urllib import request
import zipfile

from helperlibs.wrappers.io import TemporaryDirectory

VERSION_INFO = {
    "1.3": {
        "schema": "1.4",
        "source_url": "https://zenodo.org/records/13963495/files/mite-standard/mite_data-1.3.zip",
        "md5sum": "8b0ba08f82b90e8ff5528ea5a932e035",
        "base_ref_url": "https://mite.bioinformatics.nl/repository/{accession}",
    },
}

FUNCTION_MAP: dict[str, int] = {
    "Acetylation": 2,
    "Acylation": 6,
    "Amination": 2,
    "Biaryl bond formation": 6,
    "Carboxylation": 2,
    "Cyclization": 1,
    "Deamination": 3,
    "Decarboxylation": 4,
    "Dehydration": 4,
    "Demethylation": 2,
    "Deoxygenation": 1,
    "Monooxygenation": 1,
    "Dioxygenation": 1,
    "Epimerization": 5,
    "Glycosylation": 2,
    "Halogenation": 2,
    "Heterocyclization": 1,
    "Hydrolysis": 3,
    "Hydroxylation": 1,
    "Macrocyclization": 6,
    "Macrolactam formation": 6,
    "Methylation": 2,
    "Oxidation": 1,
    "Phosphorylation": 2,
    "Prenylation": 2,
    "Reduction": 1,
    "Sulfation": 2,
}


def _get_checksum(filename: str, chunk_size: int = 2 ** 20) -> str:
    md5 = hashlib.md5()
    with open(filename, "rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            md5.update(chunk)

    return md5.hexdigest()


def _download_file(url: str, filename: str, chunk_size: int = 1024 ** 2) -> None:
    if os.path.exists(filename):
        return
    with request.urlopen(url) as req:
        with open(filename, "wb") as handle:
            while True:
                try:
                    chunk = req.read(chunk_size)
                    if not chunk:
                        break
                    handle.write(chunk)
                except IOError:
                    raise RuntimeError("ERROR: Download interrupted.")


def _combine_fastas(entire_zip: zipfile.ZipFile, items: list[zipfile.ZipInfo]) -> str:
    text = []
    for item in items:
        text.append(entire_zip.read(item.filename).decode())
    return "\n".join(text)


def _convert_json(original: dict[str, Any], schema: str) -> dict[str, Any]:
    if schema in ["1.4"]:
        versions = [tuple(map(int, release["version"].split("."))) for release in original["changelog"]["releases"]]
        functions = set()
        for reaction in original["reactions"]:
            functions.update(reaction["tailoring"])
        groups = set([FUNCTION_MAP[f] for f in functions if f != "Other"])
        return {
            "accession": original["accession"],
            "description": original["enzyme"]["description"],
            "functions": sorted(functions),
            "groups": sorted(groups),
            "version": ".".join(map(str, sorted(versions)[-1])),
        }
    raise ValueError(f"unhandled schema version: {schema}")


def _extract(unzipped: zipfile.ZipFile, output_base_dir: str,
             version: str, details: dict[str, str]) -> None:
    schema_version = details["schema"]
    json_data = os.path.join("mite_data", "data", "MITE")
    jsons = {}
    fasta_data = os.path.join("mite_data", "fasta")
    fastas = {}
    for item in unzipped.infolist():
        if item.file_size <= 0:
            continue
        if json_data in item.filename:
            raw = json.loads(unzipped.read(item))
            if raw["status"] == "active":
                jsons[raw["accession"]] = _convert_json(raw, schema_version)
        elif fasta_data in item.filename and item.filename.endswith(".fasta"):
            fastas[os.path.splitext(os.path.basename(item.filename))[0]] = item

    with open(os.path.join(output_base_dir, version, "metadata.json"), "w", encoding="UTF-8") as handle:
        json.dump({
            "name": "MITE",
            "version": version,
            "entries": jsons,
            "url": details["base_ref_url"]
        }, handle)

    with open(os.path.join(output_base_dir, version, "mite.fasta"), "w", encoding="UTF-8") as handle:
        for accession in jsons:
            for line in unzipped.read(fastas[accession]).decode().splitlines():
                if ">MITE" in line:
                    line = line.split()[0]
                handle.write(line)
                handle.write("\n")


def build(version: str, output_base_dir: str) -> None:
    """ Writes relevant data files for the given version into
        the given output directory, in a subdirectory of the version.
    """
    details = VERSION_INFO.get(version)
    if not details:
        raise ValueError(f"no details available for version {version}")
    os.makedirs(os.path.join(output_base_dir, version), exist_ok=True)

    target_filename = "mite_data.zip"
    with TemporaryDirectory(change=True):
        _download_file(details["source_url"], filename=target_filename)
        if details["md5sum"] != _get_checksum(target_filename):
            raise ValueError(f"Checksum mismatch with downloaded file: {target_filename}")
        with zipfile.ZipFile(target_filename) as unzipped:
            _extract(unzipped, output_base_dir, version, details)


def main(args: argparse.Namespace) -> None:  # pylint: disable=missing-function-docstring
    try:
        build(args.version, os.path.abspath(args.output_dir))
    except Exception as err:  # pylint: disable=broad-exception-caught
        print(f"Error processing data: {str(err)}", file=sys.stderr)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("version", help="The version of MITE to fetch")
    parser.add_argument(
        "-o", "--output-dir", default="mite_data",
        help="directory to place directory of versioned data",
    )
    main(parser.parse_args())
