# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Promoter-related functions and classes for CASSIS """

import logging
from typing import List, Union

from Bio.Seq import Seq

from antismash.common.secmet import Record


class DuplicatePromoterError(Exception):
    '''Thrown when running into valid but duplicate promoter sequences during runtime'''
    pass


class Promoter:
    """ Contains all the relevant info and helpers for promoters """
    def __init__(self, gene_name: str, start: int, end: int, seq: Union[Seq, str] = None) -> None:
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
        return "Promoter(%r, %d, %d)" % (self.get_id(), self.start, self.end)


class CombinedPromoter(Promoter):
    """ A promoter class for cases where two genes are involved """
    def __init__(self, first_gene: str, second_gene: str, start: int, end: int, seq: Union[Seq, str] = None) -> None:
        super().__init__(first_gene, start, end, seq=seq)
        self.second_gene = str(second_gene)

    def get_id(self) -> str:
        return "{}+{}".format(self.gene_name, self.second_gene)

    def get_gene_names(self) -> List[str]:
        return super().get_gene_names() + [self.second_gene]


def get_promoters(record: Record, genes, upstream_tss: int, downstream_tss: int) -> List[Promoter]:
    """Compute promoter sequences for each gene in the sequence record"""
    logging.debug("Computing promoter sequences")

    def first(location):
        return min(location.start, location.end)

    def last(location):
        return max(location.start, location.end)

    min_promoter_length = 6
    max_promoter_length = (upstream_tss + downstream_tss) * 2 + 1

    record_seq_length = len(record.seq)
    promoters = []
    invalid = 0

    skip = False  # helper var for shared promoter of bidirectional genes
    for i, gene in enumerate(genes):
        gene_first = first(gene.location)
        gene_last = last(gene.location)

        downstream_inside_gene = gene_last > gene_first + downstream_tss
        if gene.location.strand not in [1, -1]:
            raise ValueError("Gene %s has unknown strand: %s" % (gene.get_name(), gene.location.strand))

        if skip:  # two genes share the same promoter --> did computation with first gene, skip second gene
            skip = False
            # TODO: should this have a continue?

        elif len(genes) == 1:  # only one gene within record
            if gene.location.strand == 1:
                upstream_inside_record = gene_first - upstream_tss >= 0
                if downstream_inside_gene:
                    if upstream_inside_record:
                        # 1 (for explanation of these numbers see file promoterregions.png)
                        # fuzzy (>|<) gene locations will be transformed to exact promoter locations
                        # we could save the fuzzy locations for promoters, too, via a FeatureLocation object
                        # but we use/calculate with the exact promoter locations anyway, here and later on
                        promoters.append(Promoter(gene.get_name(), gene_first - upstream_tss, gene_first + downstream_tss))
                    else:
                        # 2
                        promoters.append(Promoter(gene.get_name(), 0, gene_first + downstream_tss))
                else:
                    if upstream_inside_record:
                        # 3
                        promoters.append(Promoter(gene.get_name(), gene_first - upstream_tss, gene_last))
                    else:
                        # 7
                        promoters.append(Promoter(gene.get_name(), 0, gene_last))

            elif gene.location.strand == -1:
                upstream_inside_record = gene_last + upstream_tss <= record_seq_length
                if downstream_inside_gene:
                    if upstream_inside_record:
                        # 4
                        promoters.append(Promoter(gene.get_name(), gene_last - downstream_tss, gene_last + upstream_tss))
                    else:
                        # 5
                        promoters.append(Promoter(gene.get_name(), gene_last - downstream_tss, record_seq_length))
                else:
                    if upstream_inside_record:
                        # 6
                        promoters.append(Promoter(gene.get_name(), gene_first, gene_last + upstream_tss))
                    else:
                        # 8
                        promoters.append(Promoter(gene.get_name(), gene_first, record_seq_length))

        # first gene of the record AND NOT special case #9
        elif (i == 0 and not (gene.location.strand == -1
                              and genes[i+1].location.strand == 1
                              and gene_last + upstream_tss >= first(genes[i+1].location) - upstream_tss)):
            if gene.location.strand == 1:
                within_record = gene_first - upstream_tss >= 0
                if downstream_inside_gene:
                    if within_record:
                        # 1
                        promoters.append(Promoter(gene.get_name(), gene_first - upstream_tss, gene_first + downstream_tss))
                    else:
                        # 2
                        promoters.append(Promoter(gene.get_name(), 0, gene_first + downstream_tss))
                else:
                    if within_record:
                        # 3
                        promoters.append(Promoter(gene.get_name(), gene_first - upstream_tss, gene_last))
                    else:
                        # 7
                        promoters.append(Promoter(gene.get_name(), 0, gene_last))

            elif gene.location.strand == -1:
                other_first = first(genes[i+1].location)
                other_in_upstream = other_first <= gene_last + upstream_tss
                if downstream_inside_gene:
                    if other_in_upstream:
                        # 5
                        promoters.append(Promoter(gene.get_name(), gene_last - downstream_tss, other_first - 1))
                    else:
                        # 4
                        promoters.append(Promoter(gene.get_name(), gene_last - downstream_tss, gene_last + upstream_tss))
                else:
                    if other_in_upstream:
                        # 8
                        promoters.append(Promoter(gene.get_name(), gene_first, other_first - 1))
                    else:
                        # 6
                        promoters.append(Promoter(gene.get_name(), gene_first, gene_last + upstream_tss))

        # last gene of record
        elif i == len(genes) - 1 and not skip:
            if gene.location.strand == 1:
                other_last = last(genes[i-1].location)
                other_inside_upstream = other_last < gene_first - upstream_tss
                if downstream_inside_gene:
                    if other_inside_upstream:
                        # 1
                        promoters.append(Promoter(gene.get_name(), gene_first - upstream_tss, gene_first + downstream_tss))
                    else:
                        # 2
                        promoters.append(Promoter(gene.get_name(), other_last + 1, gene_first + downstream_tss))
                else:
                    if other_inside_upstream:
                        # 3
                        promoters.append(Promoter(gene.get_name(), gene_first - upstream_tss, gene_last))
                    else:
                        # 7
                        promoters.append(Promoter(gene.get_name(), other_last + 1, gene_last))

            elif gene.location.strand == -1:
                upstream_inside_record = gene_last + upstream_tss <= record_seq_length
                if downstream_inside_gene:
                    if upstream_inside_record:
                        # 4
                        promoters.append(Promoter(gene.get_name(), gene_last - downstream_tss, gene_last + upstream_tss))
                    else:
                        # 5
                        promoters.append(Promoter(gene.get_name(), gene_last - downstream_tss, record_seq_length))
                else:
                    if upstream_inside_record:
                        # 6
                        promoters.append(Promoter(gene.get_name(), gene_first, gene_last + upstream_tss))
                    else:
                        # 8
                        promoters.append(Promoter(gene.get_name(), gene_first, record_seq_length))

        # special-case 9, bidirectional promoters
        elif (gene.location.strand == -1
                and genes[i+1].location.strand == 1
                and gene_last + upstream_tss >= first(genes[i+1].location) - upstream_tss):
            other_first = first(genes[i+1].location)
            other_last = last(genes[i+1].location)
            other_downstream_in_other = other_last > other_first + downstream_tss
            if downstream_inside_gene:
                if other_downstream_in_other:
                    # 9 (1+4)
                    promoters.append(CombinedPromoter(gene.get_name(), genes[i+1].get_name(), gene_last - downstream_tss, other_first + downstream_tss))
                else:
                    # 9 (3+4)
                    promoters.append(CombinedPromoter(gene.get_name(), genes[i+1].get_name(), gene_last - downstream_tss, other_last))
            else:
                if other_downstream_in_other:
                    # 9 (1+6)
                    promoters.append(CombinedPromoter(gene.get_name(), genes[i+1].get_name(), gene_first, other_first + downstream_tss))
                else:
                    # 9 (3+6)
                    promoters.append(CombinedPromoter(gene.get_name(), genes[i+1].get_name(), gene_first, other_last))

            skip = True

        # "normal" cases
        elif not skip:
            if gene.location.strand == 1:
                other_last = last(genes[i-1].location)
                other_inside_upstream = other_last < gene_first - upstream_tss
                if downstream_inside_gene:
                    if other_inside_upstream:
                        # 1
                        promoters.append(Promoter(gene.get_name(), gene_first - upstream_tss, gene_first + downstream_tss))
                    else:
                        # 2
                        promoters.append(Promoter(gene.get_name(), other_last + 1, gene_first + downstream_tss))
                else:
                    if other_inside_upstream:
                        # 3
                        promoters.append(Promoter(gene.get_name(), gene_first - upstream_tss, gene_last))
                    else:
                        # 7
                        promoters.append(Promoter(gene.get_name(), other_last + 1, gene_last))

            elif gene.location.strand == -1:
                other_first = first(genes[i+1].location)
                other_in_upstream = other_first <= gene_last + upstream_tss
                if downstream_inside_gene:
                    if other_in_upstream:
                        # 5
                        promoters.append(Promoter(gene.get_name(), gene_last - downstream_tss, other_first - 1))
                    else:
                        # 4
                        promoters.append(Promoter(gene.get_name(), gene_last - downstream_tss, gene_last + upstream_tss))
                else:
                    if other_in_upstream:
                        # 8
                        promoters.append(Promoter(gene.get_name(), gene_first, other_first - 1))
                    else:
                        # 6
                        promoters.append(Promoter(gene.get_name(), gene_first, gene_last + upstream_tss))

        # negative start position or stop position "beyond" record --> might happen in very small records
        if promoters[-1].start < 0:
            promoters[-1].start = 0
        if promoters[-1].end > record_seq_length - 1:
            promoters[-1].end = record_seq_length - 1

        # write promoter positions and sequences to file
        if not skip:
            promoter_sequence = record.seq[promoters[-1].start:promoters[-1].end + 1]
            promoter_length = len(promoter_sequence)

            invalid_promoter_sequence = ""

            # check if promoter length is valid
            if promoter_length < min_promoter_length or promoter_length > max_promoter_length:
                invalid_promoter_sequence = "length"

            # check if a, c, g and t occur at least once in the promoter sequence
            elif "A" not in promoter_sequence.upper():
                invalid_promoter_sequence = "A"
            elif "C" not in promoter_sequence.upper():
                invalid_promoter_sequence = "C"
            elif "G" not in promoter_sequence.upper():
                invalid_promoter_sequence = "G"
            elif "T" not in promoter_sequence.upper():
                invalid_promoter_sequence = "T"

            if invalid_promoter_sequence:
                invalid += 1

                if invalid_promoter_sequence == "length":
                    logging.warning("Promoter %r is invalid (length is %s)",
                                    promoters[-1].get_id(), promoter_length)
                else:
                    # especially SiTaR doesn't like such missings
                    logging.warning("Promoter %r is invalid (sequence without %r)",
                                    promoters[-1].get_id(), invalid_promoter_sequence)

                # more details for debug logging
                logging.debug("Invalid promoter %r\n start %s\n end %s\n length %s\n",
                              promoters[-1].get_id(), promoters[-1].start,
                              promoters[-1].end, promoter_length)

                promoters.pop()  # remove last (invalid!) promoter

            else:
                promoters[-1].seq = promoter_sequence

        # check if promoter IDs are unique
        if len(promoters) >= 2 and promoters[-1].get_id() == promoters[-2].get_id():
            logging.error("Promoter %r occurs at least twice. This may be caused by overlapping gene annotations",
                          promoters[-1].get_id())
            raise DuplicatePromoterError

    if invalid:
        logging.debug("Ignoring %d promoters due to invalid promoter sequences", invalid)

    logging.debug("Found %d promoter sequences for %d genes", len(promoters), len(genes))
    return promoters
