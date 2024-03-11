# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Gene finding using Prodigal
"""

import logging
from typing import Iterable

from antismash.common.secmet import CDSFeature, Record
from antismash.common.secmet.features.cds_feature import MAX_TRANSLATION_LENGTH
from antismash.common.secmet.locations import (
    CompoundLocation,
    FeatureLocation,
    Location,
    location_contains_other,
)
from antismash.common.subprocessing.prodigal import (
    run_prodigal as exec_prodigal,
    ProdigalGene,
)


def _build_location_from_prodigal(gene: ProdigalGene, max_length: int) -> Location:
    start = gene.start
    end = gene.end
    strand = gene.strand

    if start > max_length:
        raise ValueError("start coordinate exceeds max length")

    if start < end <= max_length:
        return FeatureLocation(start, end, strand=strand)

    end %= max_length
    return CompoundLocation([
        FeatureLocation(start, max_length, strand),
        FeatureLocation(0, end, strand),
    ][::strand])  # flip exon order for the reverse strand



def _filter_genes(genes_found: Iterable[ProdigalGene], record: Record) -> Iterable[Location]:
    locations = []

    for gene in genes_found:
        # discard any duplications resulting from circular wrapping
        if gene.start > len(record) and gene.end > len(record):
            continue

        length = gene.end - gene.start
        if length > (MAX_TRANSLATION_LENGTH * 3) * .9:
            logging.warning("Ignoring potential gene too long for dependencies: %dkb", length // 1000)
            continue

        locations.append(_build_location_from_prodigal(gene, len(record)))

    # prodigal will create a gene at the start of the contig if there's a valid stop codon after it
    # so remove it, if it's contained by a larger cross-origin gene found at the extension area
    if record.is_circular() and len(locations) > 1 and location_contains_other(locations[-1], locations[0]):
        locations.pop(0)

    return locations


def run_prodigal(record: Record) -> None:
    """ Run progidal to annotate prokaryotic sequences
    """
    seq = str(record.seq)
    if record.is_circular():
        # add a buffer in which to find genes that cross the origin
        extension = min(int(len(seq) * .9), 50_000)
        seq = f"{seq}{seq[:extension]}"

    locations = _filter_genes(exec_prodigal(seq), record)
    count = 0
    for location in locations:
        translation = record.get_aa_translation_from_location(location)
        feature = CDSFeature(location, locus_tag=f"ctg{record.record_index}_{count + 1}",
                             translation=translation, translation_table=record.transl_table)
        record.add_cds_feature(feature)
        count += 1

    logging.debug("prodigal found %d CDS features", count)
