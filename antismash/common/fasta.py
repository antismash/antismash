# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" A collection of functions supporting the FASTA format
"""

from collections import OrderedDict
import logging
from typing import Dict, List

from antismash.common.secmet import CDSFeature, Record


def get_fasta_from_features(features: List[CDSFeature], numeric_names=False) -> str:
    """ Extract multi-protein FASTA from provided features """
    all_fastas = []
    if not numeric_names:
        for feature in features:
            all_fastas.append(">%s\n%s" % (feature.get_name(), feature.translation))
    else:
        for i, feature in enumerate(features):
            all_fastas.append(">%s\n%s" % (i, feature.translation))
    return "\n".join(all_fastas)


def get_fasta_from_record(record: Record) -> str:
    """ Extract multi-protein FASTA from all CDS features in sequence record """
    features = record.get_cds_features()
    all_fastas = []
    for feature in features:
        gene_id = feature.get_name()
        fasta_seq = feature.translation
        all_fastas.append(">%s\n%s" % (gene_id, fasta_seq))
    return "\n".join(all_fastas)


def write_fasta(names: List[str], seqs: List[str], filename: str) -> None:
    """ Write name/seq pairs to file in FASTA format """
    out_file = open(filename, "w")
    for name, seq in zip(names, seqs):
        out_file.write(">%s\n%s\n" % (name, seq))
    out_file.close()


def read_fasta(filename: str) -> Dict[str, str]:
    """ reads a fasta file into a dict: id -> sequence, returns the dict """
    ids = []
    sequence_info = []
    with open(filename, "r") as fasta:
        current_seq = []  # type: List[str]
        for line in fasta:
            line = line.strip()
            if not line:
                continue
            if line[0] == '>':
                ids.append(line[1:].replace(" ", "_"))
                if current_seq:
                    sequence_info.append("".join(current_seq))
                    current_seq.clear()
            else:
                if not ids:
                    raise ValueError("Sequence before identifier in fasta file")
                if not line.replace("-", "z").isalpha():
                    raise ValueError("Sequence contains non-alphabetic characters")
                current_seq.append(line)
    if current_seq:
        sequence_info.append("".join(current_seq))
    if len(ids) != len(sequence_info):
        raise ValueError("Fasta files contains different counts of sequences and ids")
    if not ids:
        logging.debug("Fasta file %s contains no sequences", filename)
        raise ValueError("Fasta file contains no sequences")
    return OrderedDict(zip(ids, sequence_info))
