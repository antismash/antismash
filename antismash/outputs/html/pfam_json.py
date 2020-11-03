# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Pfam-related JSON conversion for javascript drawing """

from typing import Any, Dict, List

from antismash.common import path
from antismash.common.secmet import PFAMDomain, Record
from antismash.outputs.html import js


TRANSPORT = set(open(path.get_full_path(__file__, "data", "pfam_transport.txt")).read().splitlines())
BIOSYNTHETIC = set(open(path.get_full_path(__file__, "data", "pfam_biosynthetic.txt")).read().splitlines())
REGULATORY = set(open(path.get_full_path(__file__, "data", "pfam_regulatory.txt")).read().splitlines())


def gather_pfam_json(records: List[Record]) -> Dict[str, Dict[str, List[Dict[str, Any]]]]:
    """ Converts Pfam annotations from detection and analysis modules into a JSON
        format for javascript to use

        Arguments:
            records: the list of records containing Pfams to convert

        Returns:
            a dictionary mapping region anchor to
                a dictionary mapping "pfamOrfs" to
                    a list of dictionaries each representing a CDS with pfams
    """
    def convert_pfam(pfam: PFAMDomain) -> Dict[str, Any]:
        """ Converts a single Pfam to it's JSON representation """
        if pfam.identifier in TRANSPORT:
            colour = "transport"
        elif pfam.identifier in REGULATORY:
            colour = "regulatory"
        elif pfam.identifier in BIOSYNTHETIC:
            colour = "biosynthetic"
        else:
            colour = "other"
        # any changes here have to be matched in antismash-js
        return {
            "start": pfam.protein_location.start,
            "end": pfam.protein_location.end,
            "name": pfam.domain,
            "accession": pfam.full_identifier,
            "description": pfam.description,
            "translation": pfam.translation,
            "evalue": f"{pfam.evalue:g}",
            "score": f"{pfam.score:.1f}",
            "go_terms": js.build_pfam2go_links(pfam.gene_ontologies),
            "html_class": f"pfam-type-{colour}",
        }

    pfams_by_region = {}
    for record_number, record in enumerate(records):
        for region in record.get_regions():
            cdses = []
            for cds in region.cds_children:
                pfams = []
                for pfam in record.get_pfam_domains_in_cds(cds):
                    pfams.append(convert_pfam(pfam))
                if pfams:
                    cdses.append({
                        "id": cds.get_name(),
                        "seqLength": len(cds.translation),
                        "pfams": pfams,
                    })
            if cdses:
                pfams_by_region[f"r{record_number + 1}c{region.get_region_number()}"] = {"pfamOrfs": cdses}
    return pfams_by_region
