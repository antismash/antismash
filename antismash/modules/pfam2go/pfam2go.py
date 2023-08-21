# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""Maps Pfam IDs of Pfam domains found in the record to Gene Ontology terms, using the Pfam domains to Gene Ontology
term mapping supplied on geneontology.org. Current mapping used: version 02/24/2018"""

import logging
from collections import defaultdict
from typing import Any, Dict, List, Optional

from antismash.common import path
from antismash.common.module_results import ModuleResults
from antismash.common.secmet.features import PFAMDomain
from antismash.common.secmet.qualifiers import GOQualifier
from antismash.common.secmet.record import Record

DATA_FILE = path.get_full_path(__file__, 'data', 'pfam2go.txt')


class GeneOntology:  # pylint: disable=too-few-public-methods
    """A single Gene Ontology term; holds Gene Ontology ID and its human-readable description."""
    def __init__(self, _id: str, description: str) -> None:
        if not _id.startswith('GO:'):
            raise ValueError(f"Invalid Gene Ontology ID: {_id}")
        self.id = _id
        assert isinstance(description, str)
        self.description = description

    def __str__(self) -> str:
        return self.id


class GeneOntologies:
    """A collection of all Gene Ontology terms (as GeneOntology objects) for a Pfam ID."""
    def __init__(self, pfam: str, gos: List[GeneOntology]) -> None:
        self.pfam = str(pfam)
        assert self.pfam.startswith('PF')
        assert isinstance(gos, list) and gos
        self.go_entries = gos

    def __str__(self) -> str:
        return str([str(go_entry) for go_entry in self.go_entries])

    def as_dict(self) -> Dict[str, str]:
        """Return a dictionary of IDs to descriptions for all GeneOntology objects in the collection."""
        return {go_entry.id: go_entry.description for go_entry in self.go_entries}


class Pfam2GoResults(ModuleResults):
    """Holds results for Pfam to Gene Ontology module."""
    schema_version = 1

    def __init__(self, record_id: str, pfam_domains_with_gos: Dict[PFAMDomain, List[GeneOntologies]]) -> None:
        super().__init__(record_id)
        #  store mapping of PFAM domain ID and GO terms
        self.pfam_domains_with_gos = pfam_domains_with_gos

    def get_all_gos(self, pfam_domain: PFAMDomain) -> Dict[str, str]:
        """"Get all Gene Ontology IDs and descriptions for a PFAMDomain.
        Arguments:
            pfam_domain: PFAMDomain for which to get the ID/description pairings
        Returns:
            A dictionary mapping Gene Ontology IDs to descriptions.
        """
        all_gos = {}
        for ontologies in self.pfam_domains_with_gos[pfam_domain]:
            all_gos.update(ontologies.as_dict())
        return all_gos

    def add_to_record(self, record: Record) -> None:
        """Add GeneOntologies objects to the respective PFAMDomains.
        Arguments:
            record: Record to which to add GeneOntologies
        """
        if record.id != self.record_id:
            raise ValueError("Record to store in and record analysed don't match")
        for domain in self.pfam_domains_with_gos:
            domain.gene_ontologies = GOQualifier(self.get_all_gos(domain))

    def to_json(self) -> Dict[str, Any]:
        """ Construct a JSON representation of this instance """
        pfam_ontologies: Dict[str, Dict[str, str]] = {}
        json_out = {"pfams": pfam_ontologies,
                    "record_id": self.record_id,
                    "schema_version": Pfam2GoResults.schema_version}
        for all_ontologies in self.pfam_domains_with_gos.values():
            for ontologies in all_ontologies:
                pfam_ontologies[ontologies.pfam] = ontologies.as_dict()
        return json_out

    @staticmethod
    def from_json(json: Dict[str, Any], record: Record) -> Optional["Pfam2GoResults"]:
        """ Constructs a new Pfam2GoResults instance from a json format and the
            original record analysed.

            Arguments:
                json: JSON representation of Pfam2GoResults
                record: Record analysed

            Returns:
                A Pfam2GoResults instance constructed from the record and the JSON
        """
        if json["schema_version"] != Pfam2GoResults.schema_version:
            logging.warning("Schema version mismatch, discarding Pfam2GO results")
            return None
        all_pfam_ids_to_ontologies: Dict[PFAMDomain, List[GeneOntologies]] = defaultdict(list)
        for domain in record.get_pfam_domains():
            id_without_version = domain.identifier
            if id_without_version in json["pfams"]:
                all_ontology = [GeneOntology(go_id, go_description)
                                for go_id, go_description in json["pfams"][id_without_version].items()]
                all_pfam_ids_to_ontologies[domain].append(GeneOntologies(id_without_version, all_ontology))
        results = Pfam2GoResults(record.id, all_pfam_ids_to_ontologies)
        return results


def construct_mapping(mapfile: str) -> Dict[str, GeneOntologies]:
    """Read a file mapping Pfam IDs to Gene Ontology terms, then convert to a dictionary matching Pfam IDs to
    collections of all Gene Ontology terms for these IDs as GeneOntologies objects.
    The mapping file must be in the following format:
    Pfam:pfam_id pfam_symbol > GO:go_readable ; go_id

    Arguments:
        mapfile: the path of the file containing the Pfam ID to GO mappings

    Returns:
        A dictionary mapping a Pfam ID to a GeneOntologies object containing GeneOntology representations of all GO
        terms matched to this ID.
    """
    results = {}
    gene_ontology_per_pfam: Dict[str, List[GeneOntology]] = defaultdict(list)
    with open(path.get_full_path(__file__, mapfile), "r", encoding="utf-8") as pfam_map:
        for line in pfam_map:
            if line.startswith('!'):
                continue
            both_info = line.split(' > ')
            pfam_info = both_info[0]
            go_info = both_info[1].strip()
            pfam_split = pfam_info.split(' ')
            pfam_id = pfam_split[0].replace('Pfam:', '')
            go_split = go_info.split(' ; ')
            go_id = go_split[1]
            go_readable = go_split[0].replace('GO:', '')
            gene_ontology_per_pfam[pfam_id].append(GeneOntology(go_id, go_readable))
    for pfam, ontology_list in gene_ontology_per_pfam.items():
        results[pfam] = GeneOntologies(pfam, ontology_list)
    return results


def get_gos_for_pfams(record: Record) -> Dict[PFAMDomain, List[GeneOntologies]]:
    """ Find Gene Ontology terms for a record's Pfam domains.

    Arguments:
        record: Record instance to annotate with Gene Ontology information

    Returns:
        A dictionary mapping a specific PFAMDomain instance to a list of GeneOntologies within the PFAMDomain.
    """
    pfam_domains_with_gos: Dict[PFAMDomain, List[GeneOntologies]] = defaultdict(list)
    pfams = record.get_pfam_domains()
    full_gomap_as_ontologies = construct_mapping(DATA_FILE)
    if not pfams:
        logging.debug('No Pfam domains found in record, cannot create Pfam to Gene Ontology mapping')
    for pfam in pfams:
        gene_ontologies_for_pfam = full_gomap_as_ontologies.get(pfam.identifier)
        if gene_ontologies_for_pfam:
            pfam_domains_with_gos[pfam].append(gene_ontologies_for_pfam)
    return pfam_domains_with_gos
