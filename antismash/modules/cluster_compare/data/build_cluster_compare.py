#!/usr/bin/env python3
import glob
import json
import os
import sys
from typing import (
    Any,
    Dict,
    IO,
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
    result = {
        "location": str(cds.location),
        "function": str(cds.gene_function),
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


def convert_all_mibig(input_dir: str, output_dir: str) -> None:
    files = [name for name in glob.glob(os.path.join(input_dir, "*", "*.gbk")) if "region" not in name]
    assert files
    result = {}
    with open(os.path.join(output_dir, "proteins.fasta"), "w") as fasta:
        for filename in files:
            # get mibig data
            with open(os.path.join(filename.split(".")[0] + ".json")) as handle:
                mibig_data = json.load(handle)
            try:
                for record in secmet.Record.from_genbank(filename):
                    record.id = record.annotations['structured_comment']['antiSMASH-Data'].get('Original ID', record.id)
                    converted = convert_record(record, fasta, skip_contig_edge=False)
                    converted["regions"][0]["products"] = mibig_data["cluster"]["biosyn_class"]
                    if mibig_data["cluster"].get("other", {}).get("subclass"):
                        converted["regions"][0]["products"] += " (%s)" % mibig_data["cluster"]["other"]["subclass"]
                    converted["regions"][0]["organism"] = mibig_data["cluster"]["organism_name"]
                    converted["regions"][0]["description"] = ", ".join([compound["compound"] for compound in mibig_data["cluster"]["compounds"]])
                    result[record.id] = converted
            except secmet.errors.SecmetInvalidInputError as err:
                print(filename, "failed:", err, list(mibig_data))
            except KeyError as err:
                print(filename, "is invalid:", err)
    with open(os.path.join(output_dir, "data.json"), "w") as handle:
        handle.write(json.dumps(result, indent=1))


if __name__ == "__main__":
    mibig = "--mibig" in sys.argv
    if mibig:
        sys.argv.pop(sys.argv.index("--mibig"))
    if len(sys.argv) == 2:
        output = os.path.join(os.getcwd(), "data")
    elif len(sys.argv) == 3:
        output = sys.argv[2]
    else:
        print((f"Usage: {sys.argv[0]} directory_containing_antismash_output_dirs db_output_dir --mibig\n"
               "Outputs to a directory named 'data' if not specified"))
        sys.exit(1)
    if not os.path.exists(output):
        os.makedirs(output)
    if mibig:
        convert_all_mibig(os.path.abspath(sys.argv[1]), output)
    else:
        convert_all(os.path.abspath(sys.argv[1]), output)
