#!/usr/bin/env python3
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Builds databases of various types for the CompaRiPPson module """

import argparse
from dataclasses import dataclass
import glob
import os
import sys
from typing import Any, Dict, List

from antismash.common import json


@dataclass
class Metadata:
    """ Required metadata for databases"""
    description_format: str
    fields: List[str]
    id_format: str
    name: str
    url: str = ""


MIBIG_METADATA = Metadata(
    description_format="@type@: @compounds@",
    fields=["accession", "compounds", "type", "locus"],
    id_format="@accession@ (@locus@)",
    name="MIBiG",
    url="https://mibig.secondarymetabolites.org/go/@accession@",
)

ASDB_METADATA = Metadata(
    description_format="@type@: @locus@",
    fields=["accession", "type", "locus", "start", "end"],
    id_format="@accession@",
    name="antiSMASH-DB",
    url="https://antismash-db.secondarymetabolites.org/area?record=@accession@&start=@start@&end=@end@",
)


@dataclass
class Entry:
    """ The basic entry details to flesh out the sequence information """
    accession: str
    locus: str
    type: str
    counter: int  # used only in generating a unique identifier
    sequence: str

    @property
    def fasta_name(self) -> str:
        """ A minimal, but unique, name for use in FASTA headers """
        return f"{self.counter}|{self.accession}"

    def to_json(self, extras: Dict[str, str] = None) -> Dict[str, str]:
        """ Converts the instance into a JSON friendly dictionary

            Arguments:
                extras: a dictionary of additional information to pack into,
                        or override, the base information

            Returns:
                a dictionary mapping field name to value
        """
        keys = ["accession", "locus", "type"]
        data = {key: getattr(self, key) for key in keys}
        if extras:
            data.update(extras)
        return data


@dataclass
class AntismashEntry(Entry):
    """ A antiSMASH-DB-specific entry to allow for start/end positions for
        linking to results
    """
    start: int
    end: int

    def to_json(self, extras: Dict[str, str] = None) -> Dict[str, str]:
        return super().to_json(extras={
            "start": str(self.start),
            "end": str(self.end),
        })


@dataclass
class MibigEntry(Entry):
    """ A MIBiG-specific entry to allow for compounds """
    compounds: str

    def to_json(self, extras: Dict[str, str] = None) -> Dict[str, str]:
        return super().to_json(extras={
            "compounds": self.compounds,
        })


def gather_entries(files: List[str]) -> List[Entry]:
    """ Gathers entries from antiSMASH output JSON files.

        Arguments:
            files: a list of filenames to gather from

        Returns:
            a list of entries
    """
    entries: List[Entry] = []

    def get_entries_from_file(json_path: str) -> None:
        """ Gather entries from a single file """
        with open(json_path, encoding="utf-8") as handle:
            data = json.load(handle)
        for record in data["records"]:
            accession = record["id"]
            # restore truncated names if possible
            if accession.endswith(".."):
                accession = record["annotations"]["accessions"][0]
            for module in RIPP_MODULES:
                results = record["modules"].get(module.__name__)
                if not results:
                    continue
                name = module.__name__.rsplit(".", 1)[-1]
                ripp_type = name.title()
                # some RiPP modules keep motifs in different layouts
                if name == "thiopeptides":
                    motifs = results["motifs"]
                elif name in ["lanthipeptides", "lassopeptides", "sactipeptides", "thiopeptides"]:
                    motifs = []
                    for _locus, motifs_for_locus in results["motifs"].items():
                        motifs.extend(motifs_for_locus)
                else:
                    raise ValueError(f"Unknown RiPP module: {name}, {module}")
                for motif in motifs:
                    assert isinstance(motif, dict), (module.__name__, motif)
                    location = location_from_string(motif["location"])
                    locus = motif["locus_tag"]
                    entries.append(AntismashEntry(accession, locus, ripp_type, len(entries) + 1,
                                                  motif["core"], location.start, location.end))

    for filename in files:
        get_entries_from_file(filename)
    return entries


def build_mibig_compound_string(compound_objects: List[Dict[str, Any]]) -> str:
    """ Builds a nice human-readable string from a list of compounds for use in
        a description of each MIBiG entry.

        Arguments:
            compound_objects: a list of compound objects from the MIBiG entry

        Returns:
            a nicely formatted string with all the compound names
    """
    compounds = [obj["compound"] for obj in compound_objects]
    assert compounds
    if len(compounds) == 1:
        compound_line = compounds[0]
    elif len(compounds) == 2:
        compound_line = " or ".join(compounds)
    else:
        compound_line = ", or".join([", ".join(compounds[:-1]), compounds[-1]])
    return compound_line


def gather_mibig(mibig_jsons: List[str]) -> List[Entry]:
    """ Gathers entries from a list of MIBiG JSON filenames.

        Arguments:
            mibig_json: a list of JSON filenames to use

        Returns:
            a list of entries
    """
    entries: List[Entry] = []

    for filename in mibig_jsons:
        with open(filename, encoding="utf-8") as handle:
            data = json.load(handle)
        if data["cluster"]["status"] != "active":
            continue
        ripp = data["cluster"].get("ripp", {})
        if not ripp:
            continue
        accession = data["cluster"]["mibig_accession"]
        compounds = build_mibig_compound_string(data["cluster"]["compounds"])
        ripp_type = ripp.get("subclass")
        precursors = ripp.get("precursor_genes", [])
        cores = {}
        for precursor in precursors:
            gene = precursor["gene_id"]
            for i, seq in enumerate(precursor["core_sequence"]):
                if seq in cores:
                    cores[seq] = "multiple genes"
                elif len(precursor["core_sequence"]) > 1:
                    cores[seq] = f"{gene} (core {i + 1})"
                else:
                    cores[seq] = gene
        for seq, gene in cores.items():
            entries.append(MibigEntry(accession, gene, ripp_type, len(entries) + 1, seq, compounds))

    return entries


def write_data(entries: List[Entry], metadata: Metadata, version: str,
               base_path: str,
               fasta_name: str = "cores.fa", metadata_name: str = "metadata.json") -> None:
    """ Writes the given entries and metadata out to a directory, creating it if necessary.
        Output path will be built from arugments as 'base_path/version/fasta_name' and
        'base_path/version/metadata_name'.

        Arguments:
            entries: the list of entries
            metadata: the core metadata dictionary defining fields, identifiers, descriptions,
                      and so on for use by antiSMASH
            version: the version number of the database, expected (but not enforce)
                     to be in a format like "4.1", this will be the name of the
                     subdirectory containing the output files
            base_path: the directory to contain the resulting database subdirectory in
            fasta_name: the filename to use for the FASTA file output
            metadata_name: the filename to use for the metadata file output

        Returns:
            None
    """
    assert base_path, "missing base output path"
    if os.path.sep in args.version:
        raise ValueError(f"Version contains path separators: {args.version}")

    real_dir = os.path.join(base_path, version)

    if not os.path.exists(real_dir):
        os.makedirs(real_dir)
    elif not os.path.isdir(real_dir):
        raise ValueError(f"Output path exists, but is not a directory: {real_dir}")

    with open(os.path.join(real_dir, fasta_name), "w", encoding="utf-8") as handle:
        for entry in entries:
            handle.write(f">{entry.fasta_name}\n{entry.sequence}\n")

    metadata_path = os.path.join(real_dir, metadata_name)
    with open(metadata_path, "w", encoding="utf-8") as handle:
        # not using json.dump here, because it'd be nice to have a line per entry for compactness,
        # while also having more than just a single line file

        # entire file opening
        handle.write("{\n")
        # write the core DB info as single line per value
        for key, value in vars(metadata).items():
            handle.write(f' "{key}": {json.dumps(value)},\n')
        handle.write(f' "version": "{version}",\n')
        # entries opening
        handle.write(' "entries": {\n')
        # all the individual entries, one per line
        entry_lines = []
        for entry in entries:
            entry_lines.append(f'  "{entry.counter}": {json.dumps(entry.to_json())}')
        handle.write(",\n".join(entry_lines))
        # entries closing
        handle.write('\n }\n')
        # entire file closing
        handle.write("}")


def main(output_dir: str, version: str, *, file_list: List[str] = None, data_dir: str = "",
         mibig: bool = False, asdb: bool = False, metadata: Metadata = None) -> None:
    """ The main entrypoint for building a database.

        One, but not both, of file_list and data_dir must be provided for the input files.

        Arguments:
            output_dir: the base path to use for writing the database output
            version: the version to annotate in metadata and to use as a subdirectory in output_dir
            file_list: a list of paths to files to use as inputs
            data_dir: a path to a directory containing inputs
            metadata: the metadata to use
            mibig: whether to parse MIBiG JSON instead of antiSMASH JSON,
                   also includes compounds in entries and defaults metadata to
                   a compatible setup
            asdb: whether to default to aS-DB compatible metadata

        Returns:
            None
    """
    if not output_dir:
        raise ValueError("Output path cannot be missing/empty")
    if mibig and asdb:
        raise ValueError("Both mibig and asdb set, only one can be valid")
    if mibig:
        if not metadata:
            metadata = MIBIG_METADATA
        entries = gather_mibig(file_list or glob.glob(os.path.join(data_dir, "*.json")))
    else:
        if not metadata:
            if asdb:
                metadata = ASDB_METADATA
            else:
                raise ValueError("Metadata missing or empty")
        entries = gather_entries(file_list or glob.glob(os.path.join(data_dir, "*", "*.json")))

    write_data(entries, metadata, version, output_dir)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--data", type=str, default="",
                        help=(
                            "Directory containing antiSMASH output directories"
                            " (or MIBiG JSON, if --mibig also given)"
                        ))
    parser.add_argument("--lof", default=None, type=argparse.FileType("r"),
                        help="A file containiner paths to JSON files, one per line")
    parser.add_argument("--version", type=str, required=True,
                        help="Specify version of the database")
    parser.add_argument("--mibig", default=False, action="store_true",
                        help="Parse MIBiG JSON instead of antiSMASH results")
    parser.add_argument("--asdb", default=False, action="store_true",
                        help="Use default antiSMASH-DB compatible metadata")
    parser.add_argument("--output-base-dir", type=str, required=True,
                        help="Base directory to add files to, a subdirectory with --version will be used")
    parser.add_argument("--metadata-core", default=None, type=argparse.FileType("r"))

    args = parser.parse_args()

    if not args.data and not args.lof:
        print("Neither '--data' or '--lof' was provided, one is required", file=sys.stderr)
        sys.exit(1)
    if args.data and args.lof:
        print("Both '--data' and '--lof' were provided, only one can be used", file=sys.stderr)
        sys.exit(1)

    if args.asdb:
        # having this import here avoids needing it at all for mibig mode
        # the "from" method of importing is required to work around a mypy issue
        from antismash.modules import (
            lanthipeptides,
            lassopeptides,
            sactipeptides,
            thiopeptides,
        )
        from antismash.common.secmet.locations import location_from_string
        RIPP_MODULES = [
            lanthipeptides,
            lassopeptides,
            sactipeptides,
            thiopeptides,
        ]

    try:
        main(args.output_base_dir, args.version, data_dir=args.data, mibig=args.mibig, asdb=args.asdb,
             file_list=args.lof.read().splitlines() if args.lof else None,
             metadata=Metadata(**json.load(args.metadata_core)) if args.metadata_core else None)
    except (ValueError, OSError) as err:
        print(err, file=sys.stderr)
        sys.exit(1)
