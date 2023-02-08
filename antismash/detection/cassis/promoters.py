# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Promoter-related functions and classes for CASSIS """

import logging
import os
from typing import Any, Dict, List, Set

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from antismash.common.secmet import Record, Gene


class DuplicatePromoterError(Exception):
    '''Thrown when running into valid but duplicate promoter sequences during runtime'''
    pass  # pylint: disable=unnecessary-pass


class Promoter:
    """ Contains all the relevant info and helpers for promoters """
    def __init__(self, gene_name: str, start: int, end: int, seq: str = "") -> None:
        self.gene_name = str(gene_name)
        self.start = int(start)
        self.end = int(end)
        self.seq = seq

    def get_id(self) -> str:
        """ Returns the id of the promoter """
        return self.gene_name

    def get_gene_names(self) -> List[str]:
        """ Returns a list of gene names attached to this promoter """
        return [self.gene_name]

    def __len__(self) -> int:
        if self.seq is None:
            raise ValueError("Requesting length of a promoter sequence which hasn't been set")
        return len(self.seq)

    def __str__(self) -> str:
        return f"Promoter({self.get_id()!r}, {self.start}, {self.end})"

    def __repr__(self) -> str:
        return str(self)

    def __eq__(self, other: Any) -> bool:
        return (isinstance(other, Promoter)
                and self.gene_name == other.gene_name
                and self.start == other.start
                and self.end == other.end
                and self.seq == other.seq)

    def to_json(self) -> Dict[str, Any]:
        """ Serialises a Promoter to a dictionary for use in JSON """
        return {"gene": self.gene_name, "start": self.start, "end": self.end,
                "seq": str(self.seq), "type": "Promoter"}

    @staticmethod
    def from_json(json: Dict[str, Any]) -> "Promoter":
        """ Deserialises a Promoter from a dictionary """
        return Promoter(str(json["gene"]), int(json["start"]), int(json["end"]), seq=str(json["seq"]))


class CombinedPromoter(Promoter):
    """ A promoter class for cases where two genes are involved """
    def __init__(self, first_gene: str, second_gene: str, start: int, end: int, seq: str = "") -> None:
        super().__init__(first_gene, start, end, seq=seq)
        self.second_gene = str(second_gene)

    def get_id(self) -> str:
        return f"{self.gene_name}+{self.second_gene}"

    def get_gene_names(self) -> List[str]:
        return super().get_gene_names() + [self.second_gene]

    def __eq__(self, other: Any) -> bool:
        return (isinstance(other, CombinedPromoter)
                and super().__eq__(other)
                and self.second_gene == other.second_gene)

    def __str__(self) -> str:
        return super().__str__().replace("Promoter(", "CombinedPromoter(")

    def to_json(self) -> Dict[str, Any]:
        core = super().to_json()
        core["second_gene"] = self.second_gene
        core["type"] = "CombinedPromoter"
        return core

    @staticmethod
    def from_json(json: Dict[str, Any]) -> "CombinedPromoter":
        return CombinedPromoter(str(json["gene"]), str(json["second_gene"]),
                                int(json["start"]), int(json["end"]), seq=str(json["seq"]))


def get_promoters(record: Record, genes: List[Gene],
                  upstream_tss: int, downstream_tss: int) -> List[Promoter]:
    """Compute promoter sequences for each gene in the sequence record.

       For explanation of these numbers see file promoterregions.png
    """
    logging.debug("Computing promoter sequences")

    min_promoter_length = 6
    max_promoter_length = (upstream_tss + downstream_tss) * 2 + 1

    record_seq_length = len(record.seq)
    promoters = []
    promoter_ids: Set[str] = set()
    invalid = 0

    skip = False  # helper var for shared promoter of bidirectional genes
    for i, gene in enumerate(genes):
        gene_first = gene.location.start
        gene_last = gene.location.end
        assert gene_first < gene_last

        downstream_inside_gene = gene_last > gene_first + downstream_tss
        if gene.location.strand not in [1, -1]:
            raise ValueError(f"Gene {gene.get_name()} has unknown strand: {gene.location.strand}")

        is_special_case = (gene.location.strand == -1 and len(genes) > i + 1 and genes[i+1].location.strand == 1
                           and gene_last + upstream_tss >= genes[i+1].location.start - upstream_tss)

        if skip:  # two genes share the same promoter --> did computation with first gene, skip second gene
            skip = False
            continue

        if len(genes) == 1:
            if gene.location.strand == 1:
                if downstream_inside_gene:  # 1, 2
                    end = gene_first + downstream_tss
                else:  # 3, 7
                    end = gene_last
                start = max(0, gene_first - upstream_tss)
            elif gene.location.strand == -1:
                if downstream_inside_gene:  # 4, 5
                    start = gene_last - downstream_tss
                else:  # 6, 8
                    start = gene_first
                end = min(record_seq_length - 1, gene_last + upstream_tss)
            current_promoter = Promoter(gene.get_name(), start, end)

        # first gene of the record AND NOT special case #9
        elif i == 0 and not is_special_case:
            if gene.location.strand == 1:
                if downstream_inside_gene:  # 1, 2
                    end = gene_first + downstream_tss
                else:  # 3, 7
                    end = gene_last
                start = max(0, gene_first - upstream_tss)
            elif gene.location.strand == -1:
                other_first = genes[i+1].location.start
                other_in_upstream = other_first <= gene_last + upstream_tss
                if downstream_inside_gene:  # 4, 5
                    start = gene_last - downstream_tss
                else:  # 6, 8
                    start = gene_first
                if other_in_upstream:  # 5, 8
                    end = other_first - 1
                else:  # 4, 6
                    end = gene_last + upstream_tss
            current_promoter = Promoter(gene.get_name(), start, end)

        # last gene of record
        elif i == len(genes) - 1:
            if gene.location.strand == 1:
                other_last = genes[i-1].location.end
                other_inside_upstream = other_last < gene_first - upstream_tss
                if downstream_inside_gene:  # 1, 2
                    end = gene_first + downstream_tss
                else:  # 3, 7
                    end = gene_last
                if other_inside_upstream:  # 1, 3
                    start = gene_first - upstream_tss
                else:  # 2, 7
                    start = other_last + 1
            elif gene.location.strand == -1:
                if downstream_inside_gene:  # 4, 5
                    start = gene_last - downstream_tss
                else:  # 6, 8
                    start = gene_first
                end = min(record_seq_length - 1, gene_last + upstream_tss)
            current_promoter = Promoter(gene.get_name(), start, end)

        # special-case 9, bidirectional promoters
        elif is_special_case:
            other_first = genes[i+1].location.start
            other_last = genes[i+1].location.end
            other_downstream_in_other = other_last > other_first + downstream_tss
            if downstream_inside_gene:  # 9.. (1+4), (3+4)
                start = gene_last - downstream_tss
            else:  # 9.. (1+6), (3+6)
                start = gene_first
            if other_downstream_in_other:  # 9.. (1+4), (1+6)
                end = other_first + downstream_tss
            else:  # 9.. (3+4), (3+6)
                end = other_last
            current_promoter = CombinedPromoter(gene.get_name(), genes[i+1].get_name(), start, end)

            skip = True

        # "normal" cases
        else:
            if gene.location.strand == 1:
                other_last = genes[i-1].location.end
                other_inside_upstream = other_last < gene_first - upstream_tss
                if downstream_inside_gene:  # 1, 2
                    end = gene_first + downstream_tss
                else:  # 3, 7
                    end = gene_last
                if other_inside_upstream:  # 1, 3
                    start = gene_first - upstream_tss
                else:  # 2, 7
                    start = other_last + 1
            elif gene.location.strand == -1:
                other_first = genes[i+1].location.start
                other_in_upstream = other_first <= gene_last + upstream_tss
                if downstream_inside_gene:  # 5, 4
                    start = gene_last - downstream_tss
                else:  # 6, 8
                    start = gene_first
                if other_in_upstream:  # 5, 8
                    end = other_first - 1
                else:  # 4, 6
                    end = gene_last + upstream_tss
            current_promoter = Promoter(gene.get_name(), start, end)

        # ensure all promoters have locations within the record
        current_promoter.start = max(0, current_promoter.start)
        current_promoter.end = min(record_seq_length - 1, current_promoter.end)

        current_promoter.seq = record.seq[current_promoter.start:current_promoter.end + 1]

        if is_invalid_promoter_sequence(current_promoter, min_promoter_length, max_promoter_length):
            invalid += 1
            continue

        promoters.append(current_promoter)

        # check if promoter IDs are unique
        if current_promoter.get_id() in promoter_ids:
            logging.error("Promoter %r occurs at least twice. This may be caused by overlapping gene annotations",
                          current_promoter.get_id())
            raise DuplicatePromoterError
        promoter_ids.add(current_promoter.get_id())

    if invalid:
        logging.debug("Ignoring %d promoters due to invalid promoter sequences", invalid)

    logging.debug("Found %d promoter sequences for %d genes", len(promoters), len(genes))
    return promoters


def is_invalid_promoter_sequence(promoter: Promoter, min_length: int, max_length: int) -> bool:
    """ Checks if a promoter's sequence is valid.

        A valid promoter sequence must be within the given length range and
        contain at least one of each A, C, G, and T.
    """
    promoter_sequence = promoter.seq.upper()

    invalid_promoter_sequence = ""

    # check if promoter length is valid
    if not min_length <= len(promoter) <= max_length:
        invalid_promoter_sequence = "length"
        logging.warning("Promoter %r is invalid (length is %s)",
                        promoter.get_id(), len(promoter))

    # check if a, c, g and t occur at least once in the promoter sequence
    elif "A" not in promoter_sequence:
        invalid_promoter_sequence = "A"
    elif "C" not in promoter_sequence:
        invalid_promoter_sequence = "C"
    elif "G" not in promoter_sequence:
        invalid_promoter_sequence = "G"
    elif "T" not in promoter_sequence:
        invalid_promoter_sequence = "T"

    if invalid_promoter_sequence:
        if invalid_promoter_sequence != "length":
            # especially SiTaR doesn't like missing bases
            logging.warning("Promoter %r is invalid (sequence without %r)",
                            promoter.get_id(), invalid_promoter_sequence)

        # more details for debug logging
        logging.debug("Invalid promoter %r\n start %s\n end %s\n length %s\n",
                      promoter.get_id(), promoter.start,
                      promoter.end, len(promoter))

    return bool(invalid_promoter_sequence)


def write_promoters_to_file(output_dir: str, prefix: str, promoters: List[Promoter]) -> None:
    """ Write the promoters to the given file using the given prefix.

        The prefix helps in having separate files for multi-record inputs with
        a singular output directory.
    """
    # pylint: disable=consider-using-with
    # positions file
    pos_handle = open(os.path.join(output_dir, prefix + "_promoter_positions.csv"), "w", encoding="utf-8")
    pos_handle.write("\t".join(["#", "promoter", "start", "end", "length"]) + "\n")
    # sequences file
    seq_handle = open(os.path.join(output_dir, prefix + "_promoter_sequences.fasta"), "w", encoding="utf-8")
    # pylint: enable=consider-using-with

    for i, promoter in enumerate(promoters):
        # write promoter positions to file
        pos_handle.write("\t".join(map(str, [i + 1, promoter.get_id(),
                                             promoter.start + 1, promoter.end + 1,
                                             len(promoter)])) + "\n")

        # write promoter sequences to file
        SeqIO.write(SeqRecord(promoter.seq, id=promoter.get_id(),
                    description=f"length={len(promoter)}bp"),
                    seq_handle,
                    "fasta")
    pos_handle.close()
    seq_handle.close()


def get_anchor_promoter_index(anchor: str, promoters: List[Promoter]) -> int:
    """Find the name of the promoter which includes the anchor gene"""
    # the promoter ID is not (necessarily) equal to the anchor ID!
    for i, promoter in enumerate(promoters):
        if anchor in promoter.get_gene_names():
            return i

    raise ValueError(f"No promoter exists for the given anchor: {anchor}")
