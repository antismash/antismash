#!/usr/bin/env python3

""" Script to generate a general Clusterblast dataset from existing antiSMASH outputs

    Taken almost-verbatim from https://github.com/antismash/generate-databases, where it
    was regularly out of date.
"""


import argparse
from dataclasses import dataclass
import os
import sys
from typing import (
    IO,
    Optional,
    Self,
)

import antismash
from antismash.common import secmet
from antismash.detection import hmm_detection


DEFAULT_AS_OPTIONS = antismash.config.build_config(["--minimal"], modules=antismash.main.get_all_modules())
antismash.config.update_config({"all_enabled_modules": [hmm_detection]})


@dataclass
class AsdbLocus:   # pylint: disable=too-many-instance-attributes
    """ A single gene within a reference """
    start: int
    end: int
    strand: int
    accession: str
    annotation: str
    sequence: str
    identifier: str
    gene_id: Optional[str] = None
    locus_tag: Optional[str] = None
    protein_accession: Optional[str] = None
    draw_start: int = 0
    draw_end: int = 0

    def __post_init__(self) -> None:
        if not any([self.locus_tag, self.gene_id, self.protein_accession]):
            raise RuntimeError("No valid identifier for feature")
        if not self.sequence:
            raise RuntimeError("No valid sequence")
        # strip any ambiguous positions
        self.start = int(self.start)
        self.end = int(self.end)
        self.draw_start = self.draw_start or self.start
        self.draw_end = self.draw_end or self.end

    @classmethod
    def from_secmet(cls, feature: secmet.CDSFeature, accession: str,
                    draw_start: int = 0, draw_end: int = 0,
                    ) -> Self:
        """ Create an instance from the given CDSFeature

            Arguments:
                feature: the secmet feature to convert
                accession: the record's accession
                draw_start: the effective coordinates within the parent region,
                            for when the parent region crosses the origin
                draw_end: the effective coordinates within the parent region,
                          for when the parent region crosses the origin
        """
        annotation = feature.product or "(unknown)"

        optional_values: dict[str, Optional[str]] = {}

        if feature.locus_tag:
            locus_tag = feature.locus_tag
            if locus_tag.find('allorf') > -1:
                locus_tag = f"{accession}_{locus_tag}"
            optional_values['locus_tag'] = locus_tag

        optional_values['gene_id'] = feature.gene
        optional_values['protein_accession'] = feature.protein_id

        return cls(feature.start, feature.end, feature.strand, accession, annotation, feature.translation,
                   feature.get_name(), draw_start=draw_start, draw_end=draw_end,
                   **optional_values)

    @property
    def safe_locus_tag(self) -> str:
        """ The locus tag, or the legacy value when not present """
        return self.locus_tag or "no_locus_tag"

    @property
    def safe_accession(self) -> str:
        """ The protein accession, or the legacy value when not present """
        return self.protein_accession or self.identifier

    @property
    def safe_annotation(self) -> str:
        """" User-friendly description of the gene """
        ann = self.annotation.replace(" ", "_")
        if ann == "(unknown)":
            ann = ""
        return ann


@dataclass
class AsdbRegion:
    """ A region within an antiSMASH result """
    accession: str
    start: int
    end: int
    products: str
    description: str
    loci: list[AsdbLocus]

    @classmethod
    def from_secmet(cls, record: secmet.Record, region: secmet.Region,
                    ) -> Optional[Self]:
        """ Creates an instance from the given secmet Region """
        accession = record.annotations['accessions'][0]

        description = record.description

        if len(region.location) > 500000:
            print("Region too large:", record.id, region)
            return None

        loci: list[AsdbLocus] = []
        if region.crosses_origin():
            for cds in region.cds_children.pre_origin:
                loci.append(AsdbLocus.from_secmet(cds, accession,
                                                  draw_start=cds.start, draw_end=cds.end))
            for cds in region.cds_children.cross_origin:
                loci.append(AsdbLocus.from_secmet(cds, accession,
                                                  draw_start=cds.start, draw_end=cds.end + len(record)))
            for cds in region.cds_children.post_origin:
                loci.append(AsdbLocus.from_secmet(cds, accession,
                                                  draw_end=cds.end + len(record),
                                                  draw_start=cds.start + len(record)))
        else:
            loci.extend(AsdbLocus.from_secmet(cds, accession) for cds in region.cds_children)

        return cls(accession, int(region.start), int(region.end),
                   region.get_product_string(), description, loci)

    def write_proteins(self, handle: IO) -> None:
        """Write all loci to handle in FASTA format with long headers in legacy format."""
        for locus in self.loci:
            strand = "+" if locus.strand == 1 else "-"
            header = (
                # region identifier/location
                f">{locus.accession}|c{self.start}-{self.end}|"
                # gene location
                f"{locus.start + 1}-{locus.end}|{strand}|linearised{locus.draw_start + 1}-{locus.draw_end}|"
                # gene labels
                f"{locus.safe_locus_tag}|{locus.safe_annotation}|{locus.safe_accession}\n"
            )
            handle.write(header)
            handle.write(locus.sequence)
            handle.write('\n')

    def write_regions(self, handle: IO) -> None:
        """Write region info in tabular format."""
        elements = [
            self.accession,
            self.description,
            "c{s.start}-{s.end}".format(s=self),
            self.products,
            ";".join(map(lambda loc: loc.safe_locus_tag, self.loci)),
            ";".join(map(lambda loc: loc.safe_accession, self.loci)),
        ]
        handle.write("\t".join(elements) + "\n")


def run(files: list[str], output_dir: str) -> None:
    """ Generates a dataset from the given inputs.

        Arguments:
            files: a list of GBK antiSMASH output files to include
            output_dir: the directory in which to place the clusterblast-compatible database
    """

    regions: list[AsdbRegion] = []

    for filename in files:
        for record in secmet.Record.from_genbank(filename):
            for secmet_region in record.get_regions():
                if secmet_region.contig_edge:
                    continue
                converted = AsdbRegion.from_secmet(record, secmet_region)
                if converted:
                    regions.append(converted)

    regions.sort(key=lambda x: (x.accession, x.start))

    os.makedirs(output_dir, exist_ok=True)
    with open(os.path.join(output_dir, "clusters.txt"), "w", encoding="utf-8") as handle:
        for region in regions:
            region.write_regions(handle)
    with open(os.path.join(output_dir, "proteins.fasta"), "w", encoding="utf-8") as handle:
        for region in regions:
            region.write_proteins(handle)


def main() -> int:
    """ The main entry point, handling file finding when not explicitly given filenames """
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--dir", metavar="DIRECTORY",
                       help="Directory containing the antiSMASH DB output folders")
    group.add_argument("--lof", metavar="FILE",
                       help="List containing paths of antiSMASH GBK output files to parse")
    parser.add_argument("--output-dir", metavar="DIRECTORY", required=True,
                        help="Output directory")

    args = parser.parse_args()

    file_list: list[str] = []

    if args.dir is not None:
        for root, _, files in os.walk(args.dir):
            for filename in files:
                if filename.endswith(".gbk") and "region" not in filename:
                    full_name = os.path.join(root, filename)
                    file_list.append(full_name)
                    print("found", full_name)
    else:
        with open(args.lof, "r", encoding="utf-8") as handle:
            file_list = list(map(str.strip, handle.readlines()))
    if not file_list:
        print("No files listed or found", file=sys.stderr)
        return 1

    run(file_list, args.output_dir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
