#!/usr/bin/env python3
import argparse
import glob
import json
import os
from typing import (
    Any,
    Dict,
    IO,
    List,
)

from antismash.common import secmet


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
    }
    return result


def convert_protocluster(protocluster: secmet.Protocluster) -> Dict[str, Any]:
    result = {
        "product": protocluster.product,
        "core_cdses": [cds.get_name() for cds in protocluster.definition_cdses],
        "location": str(protocluster.location),
    }
    return result


def convert_region(region: secmet.Region, cds_mapping: Dict[int, str], cds_index: Counter, fasta: IO) -> Dict[str, Any]:
    result = {
        "products": region.products,
        "protoclusters": [convert_protocluster(pc) for pc in region.get_unique_protoclusters()],
        "cdses": {cds.get_name(): convert_cds(cds) for cds in region.cds_children},
        "start": min(cds.location.start for cds in region.cds_children),  # trim any intergenic areas
        "end": max(cds.location.end for cds in region.cds_children),
    }
    for cds in region.cds_children:
        index = cds_index.next()
        fasta.write(">%s|%d\n%s\n" % (region.parent_record.id, index, cds.translation))
        cds_mapping[index] = cds.get_name()
    return result


def convert_record(record: secmet.Record, fasta: IO, skip_contig_edge: bool = True) -> Dict[str, Any]:
    result = {
        "regions": [],
        "cds_mapping": {},
    }  # type: Dict[str, Any]
    cds_index = Counter()
    for region in record.get_regions():
        if skip_contig_edge and region.contig_edge:
            continue
        result["regions"].append(convert_region(region, result["cds_mapping"], cds_index, fasta))
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
        handle.write(json.dumps(result, indent=1))


def convert_all_mibig(input_dir: str, output_dir: str, accessions: List[str]) -> None:
    if not accessions:
        raise ValueError("no valid MIBiG accessions, a database could not be generated")
    result = {}
    with open(os.path.join(output_dir, "proteins.fasta"), "w") as fasta:
        for accession in accessions:
            # get mibig data
            json_path = os.path.join(input_dir, accession, accession + ".json")
            with open(json_path) as handle:
                mibig_data = json.load(handle)
            genbank = os.path.join(input_dir, accession, "generated", accession + ".gbk")
            try:
                for record in secmet.Record.from_genbank(genbank):
                    record.id = record.annotations['structured_comment']['antiSMASH-Data'].get('Original ID', record.id)
                    converted = convert_record(record, fasta, skip_contig_edge=False)
                    converted["regions"][0]["products"] = mibig_data["cluster"]["biosyn_class"]
                    if mibig_data["cluster"].get("other", {}).get("subclass"):
                        index = converted["regions"][0]["products"].index("Other")
                        converted["regions"][0]["products"][index] += " (%s)" % mibig_data["cluster"]["other"]["subclass"]
                    converted["regions"][0]["organism"] = mibig_data["cluster"]["organism_name"]
                    converted["regions"][0]["description"] = ", ".join([compound["compound"] for compound in mibig_data["cluster"]["compounds"]])
                    result[record.id] = converted
            except secmet.errors.SecmetInvalidInputError as err:
                print(accession, "failed:", err, list(mibig_data))
            except KeyError as err:
                print(accession, "is invalid:", err)
    with open(os.path.join(output_dir, "data.json"), "w") as handle:
        json.dump(result, handle, indent=1, sort_keys=True)


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
