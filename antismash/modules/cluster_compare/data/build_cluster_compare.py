#!/usr/bin/env python3
import argparse
import glob
import os
from typing import (
    Any,
    Dict,
    IO,
    List,
)

from antismash.common import json, secmet


class Counter:
    def __init__(self, start: int = 0) -> None:
        self.current = start

    def next(self) -> int:
        self.current += 1
        return self.current


def convert_module(module: secmet.Module) -> Dict[str, Any]:
    domains = []
    for domain in module.domains:
        assert domain.domain
        # don't include Condensation domain subtypes, since PKS subtypes aren't included
        if domain.domain.startswith("Condensation_"):
            domains.append("Condensation")
        else:
            domains.append(domain.domain)
    return {
        "domains": domains,
        "type": str(module.module_type),
        "complete": module.is_complete(),
    }


def convert_cds(cds: secmet.CDSFeature) -> Dict[str, Any]:
    function = secmet.GeneFunction.OTHER
    if cds.gene_function == secmet.GeneFunction.CORE:
        function = secmet.GeneFunction.CORE
    else:
        annotations = cds.gene_functions.get_by_tool("smcogs")
        if not annotations:
            annotations = cds.gene_functions.get_by_tool("rule-based-clusters")
        if annotations:
            function = annotations[0].function

    result = {
        "location": str(cds.location),
        "function": str(function),
        "components": {
            "secmet": [] if not cds.sec_met else cds.sec_met.domain_ids,
            "modules": [convert_module(module) for module in cds.modules],
        },
        "draw_start": cds.start,
        "draw_end": cds.end,
    }
    return result


def convert_protocluster(protocluster: secmet.Protocluster) -> Dict[str, Any]:
    result = {
        "product": protocluster.product,
        "core_cdses": [cds.get_name() for cds in protocluster.definition_cdses],
        "location": str(protocluster.location),
    }
    return result


def convert_region(region: secmet.Region, cds_mapping: Dict[int, str], cds_index: Counter, fasta: IO,
                   *, record_length: int = 0, accession_override: str = "") -> dict[str, Any]:
    # trim leading/trailing intergenic areas, accounting for possible cross-origin regions
    if region.crosses_origin():
        start = min(cds.start for cds in region.cds_children.pre_origin + region.cds_children.cross_origin)
        end = max(cds.end for cds in region.cds_children.post_origin + region.cds_children.cross_origin)
    else:
        start = min(cds.start for cds in region.cds_children)
        end = max(cds.end for cds in region.cds_children)

    def _convert_cds(feature: secmet.CDSFeature) -> dict[str, Any]:
        converted = convert_cds(feature)
        if region.crosses_origin():
            if feature.crosses_origin() or feature.end < feature.start:
                converted["draw_end"] = feature.end + record_length
            if feature in region.cds_children.post_origin:
                converted["draw_start"] = feature.start + record_length
                converted["draw_end"] = feature.end + record_length
        assert converted["draw_start"] < converted["draw_end"]
        return converted

    result = {
        "products": region.products,
        "protoclusters": [convert_protocluster(pc) for pc in region.get_unique_protoclusters()],
        "cdses": {cds.get_name(): _convert_cds(cds) for cds in region.cds_children},
        "start": start,
        "end": end,
    }
    if start > end:
        result["draw_end"] = end + record_length
    for cds in region.cds_children:
        index = cds_index.next()
        accession = accession_override or region.parent_record.id
        fasta.write(f">{accession}|{index}\n{cds.translation}\n")
        cds_mapping[index] = cds.get_name()
    return result


def convert_record(record: secmet.Record, fasta: IO, skip_contig_edge: bool = True, accession_override: str = "") -> Dict[str, Any]:
    result = {
        "regions": [],
        "cds_mapping": {},
    }  # type: Dict[str, Any]
    cds_index = Counter()
    for region in record.get_regions():
        if skip_contig_edge and region.contig_edge:
            continue
        result["regions"].append(convert_region(region, result["cds_mapping"], cds_index, fasta, record_length=len(record), accession_override=accession_override))
    return result


def convert_all(input_dir: str, output_dir: str) -> None:
    files = [name for name in glob.glob(os.path.join(input_dir, "*", "*.gbk")) if "region" not in name]
    assert files
    result = {}
    with open(os.path.join(output_dir, "proteins.fasta"), "w") as fasta:
        for filename in files:
            for record in secmet.Record.from_genbank(filename):
                record.id = record.annotations['structured_comment']['antiSMASH-Data'].get('Original ID', record.id)
                if record.id.startswith("c000"):
                    print(filename, "has bad record IDs")
                    continue
                result[record.id] = convert_record(record, fasta)
    with open(os.path.join(output_dir, "data.json"), "w") as handle:
        handle.write(json.dumps(result, indent=True))


def convert_all_mibig(input_dir: str, output_dir: str, accessions: List[str]) -> None:
    if not accessions:
        raise ValueError("no valid MIBiG accessions, a database could not be generated")
    result = {}
    with open(os.path.join(output_dir, "proteins.fasta"), "w") as fasta:
        for accession in accessions:
            # get mibig data
            json_path = os.path.join(input_dir, accession, "annotations.json")
            if not os.path.exists(json_path):
                continue
            with open(json_path) as handle:
                mibig_data = json.load(handle)
            if mibig_data["status"] != "active":
                continue
            genbank = os.path.join(input_dir, accession, accession.split('.')[0] + ".gbk")
            try:
                for record in secmet.Record.from_genbank(genbank):
                    record.id = record.annotations['structured_comment']['antiSMASH-Data'].get('Original ID', record.id)
                    converted = convert_record(record, fasta, skip_contig_edge=False, accession_override=accession)
                    classes = [c["class"] for c in mibig_data["biosynthesis"]["classes"]]
                    converted["regions"][0]["products"] = list(sorted(set(classes)))
                    if "other" in classes:
                        index = converted["regions"][0]["products"].index("other")
                        converted["regions"][0]["products"][index] += " (%s)" % mibig_data["biosynthesis"]["classes"][classes.index("other")]["subclass"]
                    converted["regions"][0]["organism"] = mibig_data["taxonomy"]["name"]
                    converted["regions"][0]["description"] = ", ".join([compound["name"] for compound in mibig_data["compounds"]])
                    result[accession] = converted
            except secmet.errors.SecmetInvalidInputError as err:
                print(accession, "failed:", err, list(mibig_data))
            except KeyError as err:
                print(accession, "is invalid:", err)
    with open(os.path.join(output_dir, "data.json"), "w") as handle:
        json.dump(result, handle, indent=True, sort_keys=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "input_path",
        type=str,
        help=("The directory containing antiSMASH result subdirectories"
              " (or MIBiG directories with MIBiG JSON and antiSMASH 'generated'"
              " subdirectories)")
    )
    parser.add_argument(
        "--output",
        default="data",
        metavar="DIR",
        type=str,
        help="The directory to place the resulting database files (default: %(default)s)"
    )
    parser.add_argument(
        "--mibig",
        default="",
        metavar="ACCESSIONS_FILE",
        type=str,
        help="Run in MIBiG mode including the accessions in the given file",
    )
    args = parser.parse_args()

    if not os.path.exists(args.output):
        os.makedirs(args.output)
    if args.mibig:
        with open(args.mibig) as _handle:
            convert_all_mibig(os.path.abspath(args.input_path), args.output,
                              _handle.read().splitlines())
    else:
        convert_all(os.path.abspath(args.input_path), args.output)
