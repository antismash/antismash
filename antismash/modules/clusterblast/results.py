# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Contains the results class for each variant of clusterblast """


from collections import OrderedDict
import logging
import os
from typing import Any, Dict, IO, List, Optional, Tuple

from antismash.common.module_results import ModuleResults
from antismash.common.layers import AbstractRelatedArea
from antismash.common.path import changed_directory
from antismash.common.secmet import Record, Region
from antismash.config import ConfigType, get_config

from .data_structures import Score, Query, Subject, ReferenceCluster, Protein, MibigEntry

_CLUSTER_LIMIT = 50


def get_result_limit() -> int:
    """ Returns the hard limit of cluster matches to report """
    return _CLUSTER_LIMIT


def get_display_limit() -> int:
    """ Returns the maximum number of matches to display in visualisations """
    return get_config().cb_nclusters


class KnownHitSummary(AbstractRelatedArea):
    """ Stores some information about a hit from known-CB in a handy to access way """
    def __init__(self, bgc_id: str, name: str, cluster_number: str, similarity: int,
                 cluster_type: str) -> None:
        super().__init__()
        self._bgc_id = bgc_id
        self._name = name
        self._cluster_number = cluster_number
        self._similarity = similarity
        self._cluster_type = cluster_type

    @property
    def identifier(self) -> str:
        return self._bgc_id

    @property
    def description(self) -> str:
        return self._name

    @property
    def product(self) -> str:
        return self._cluster_type

    @property
    def similarity_percentage(self) -> int:
        return self._similarity

    @property
    def url(self) -> str:
        # at time of writing, knownclusterblast hits are always MIBiG
        return f"https://mibig.secondarymetabolites.org/go/{self.identifier}"


class RegionResult:
    """ Stores results for a specific cluster in a record, for a particular
        flavour of clusterblast.
    """
    __slots__ = ["region", "ranking", "total_hits", "prefix", "displayed_reference_proteins"]

    def __init__(self, region: Region, ranking: List[Tuple[ReferenceCluster, Score]],
                 reference_proteins: Dict[str, Protein], prefix: str) -> None:
        """ Arguments:
                cluster: the cluster feature
                ranking: a list of tuples in the form (ReferenceCluster, Score)
                reference_proteins: used to generate details for SVG output,
                                    only relevant portions are stored
                prefix: an identifier for use in marking SVGs such that
                        javascript on the results page can differentiate between
                        types of clusterblast
        """
        record_prefix = (region.parent_record.original_id or region.parent_record.id).split(".", 1)[0]
        # remove self-hits
        if prefix != "subclusterblast":
            ranking = list(filter(lambda pair: pair[0].accession != record_prefix, ranking))
        self.region = region
        self.ranking = ranking[:get_result_limit()]  # [(ReferenceCluster, Score),...]
        self.total_hits = len(ranking)
        self.prefix = prefix
        self.displayed_reference_proteins = {}

        # for the SVG portion, limit the ranking to the display limit
        display_limit = get_display_limit()
        display_ranking = self.ranking[:display_limit]

        for ref_cluster, _ in display_ranking:
            for name in ref_cluster.tags:
                prot = reference_proteins[f"{ref_cluster.accession}_{name}"]
                self.displayed_reference_proteins[prot.full_name] = prot

        assert len(display_ranking) <= display_limit

    def get_best_match(self) -> Optional[KnownHitSummary]:
        """ Returns the single best match from knownclusterblast hits, if any """
        if not self.ranking:
            return None
        reference, score = self.ranking[0]
        return KnownHitSummary(reference.accession, reference.description,
                               reference.cluster_label,
                               score.similarity, reference.cluster_type)

    def jsonify(self) -> Dict[str, Any]:
        """ Convert the object into a simple dictionary for use in storing
            results.

            The function from_json() should reconstruct a new and equal
            RegionResult from the results of this function.

            Returns:
                a dict containing the object data in basic types
        """
        ranking = []
        for cluster, score in self.ranking:
            scoring = {key.lstrip("_"): getattr(score, key) for key in score.__slots__ if key != "scored_pairings"}
            json_cluster = {key: getattr(cluster, key) for key in cluster.__slots__}
            scoring["pairings"] = [(query.entry, query.index, vars(subject))
                                   for query, subject in score.scored_pairings]
            ranking.append((json_cluster, scoring))

        result = {"region_number": self.region.get_region_number(),
                  "total_hits": self.total_hits,
                  "ranking": ranking,
                  "prefix": self.prefix}

        return result

    @staticmethod
    def from_json(json: Dict[str, Any], record: Record,
                  reference_proteins: Dict[str, Protein]) -> "RegionResult":
        """ Convert a simple dictionary into a new RegionResult object.

            The function RegionResult.jsonify() should reconstruct the data
            provided here.

            Arguments:
                json: the dict of data to construct with
                record: the record used to create the data
                reference_proteins: a dict mapping protein name to Protein,
                                    used instead of duplicated storing of
                                    many Protiens

            Returns:
                a dict containing the object data in basic types
        """
        ranking = []
        for cluster, details in json["ranking"]:
            ref_cluster = ReferenceCluster(cluster["accession"], cluster["cluster_label"],
                                           cluster["proteins"], cluster["description"],
                                           cluster["cluster_type"], cluster["tags"])
            score = Score()
            pairings = details["pairings"]
            score.similarity = details.pop("similarity")
            for key, val in details.items():
                if key == "pairings":
                    continue
                setattr(score, key, val)
            for pairing in pairings:  # entry, index, dict of subject
                query = Query(pairing[0], pairing[1])
                subject = Subject.from_dict(pairing[2])
                score.scored_pairings.append((query, subject))
            ranking.append((ref_cluster, score))

        region = record.get_region(json["region_number"])
        result = RegionResult(region, ranking, reference_proteins, json["prefix"])
        result.total_hits = json["total_hits"]
        return result


class GeneralResults(ModuleResults):
    """ A variant-agnostic results class for clusterblast variants """
    schema_version = 5

    def __init__(self, record_id: str, search_type: str = "clusterblast",
                 data_version: str = None) -> None:
        assert record_id and isinstance(record_id, str)
        assert search_type and isinstance(search_type, str)
        super().__init__(record_id)
        self.region_results: List[RegionResult] = []
        self.search_type = search_type
        # keep here instead of duplicating in clusters
        # and only keep those that are relevant instead of 7 million
        self.proteins_of_interest: Dict[str, Protein] = OrderedDict()
        # hold mappings of cluster number -> protein name -> mibig entries
        self.mibig_entries: Dict[int, Dict[str, List[MibigEntry]]] = {}
        self.data_version: Optional[str] = data_version

    def add_region_result(self, result: RegionResult, reference_clusters: Dict[str, ReferenceCluster],
                          reference_proteins: Dict[str, Protein]) -> None:
        """ Add a result for a specific cluster. Reference proteins and clusters
            required in order to write the relevant information. Only those used
            are required.

            Arguments:
                result: the RegionResult to add
                reference_clusters: a dictionary mapping reference cluster name to ReferenceCluster
                reference_proteins: a dictionary mapping reference protein name to Protein

            Returns:
                None
        """
        assert isinstance(result, RegionResult)
        assert reference_clusters and reference_proteins  # {str:ReferenceCluster}, {str:Protein}
        self.region_results.append(result)
        # keep data on proteins relevant to this cluster
        for cluster, _ in result.ranking:
            cluster_label = cluster.get_name()
            protein_names = reference_clusters[cluster_label].tags
            for protein_name in protein_names:
                protein_name = f"{cluster.accession}_{protein_name}"
                self.proteins_of_interest[protein_name] = reference_proteins[protein_name]

    def write_to_file(self, record: Record, options: ConfigType) -> None:
        """ Write the results to a text file """
        for cluster in self.region_results:
            write_clusterblast_output(options, record, cluster,
                                      self.proteins_of_interest,
                                      searchtype=self.search_type)

    def to_json(self) -> Dict[str, Any]:
        if not self.region_results:
            return {}
        data = {"record_id": self.record_id,
                "schema_version": self.schema_version,
                "results": [res.jsonify() for res in self.region_results],
                "proteins": [protein.to_json() for protein in self.proteins_of_interest.values()],
                "search_type": self.search_type}
        if self.data_version:
            data["data_version"] = self.data_version
        if self.mibig_entries:
            entries = {}
            for cluster_number, proteins in self.mibig_entries.items():
                cluster = {}
                for protein, protein_entries in proteins.items():
                    cluster[protein] = [protein_entry.values for protein_entry in protein_entries]
                entries[cluster_number] = cluster
            data["mibig_entries"] = entries
        return data

    def add_to_record(self, _record: Record) -> None:
        pass

    @staticmethod
    def from_json(json: Dict[str, Any], record: Record) -> "GeneralResults":
        if json["schema_version"] != GeneralResults.schema_version:
            raise ValueError(f"Incompatible results schema version, expected {GeneralResults.schema_version}")
        assert record.id == json["record_id"]
        data_version = json.get("data_version")
        result = GeneralResults(json["record_id"], search_type=json["search_type"],
                                data_version=data_version)
        for prot in json["proteins"]:
            protein = Protein.from_json(prot)
            result.proteins_of_interest[protein.full_name] = protein
        for region_result in json["results"]:
            result.region_results.append(RegionResult.from_json(region_result,
                                         record, result.proteins_of_interest))
        if "mibig_entries" in json:
            entries: Dict[int, Dict[str, List[MibigEntry]]] = {}
            for cluster_number, proteins in json["mibig_entries"].items():
                entries[int(cluster_number)] = {}
                for protein_name, protein_entries in proteins.items():
                    entries[int(cluster_number)][protein_name] = [MibigEntry(*entry) for entry in protein_entries]
            result.mibig_entries = entries
        return result


class ClusterBlastResults(ModuleResults):
    """ An aggregate results container for all variants of clusterblast """
    schema_version = 2

    def __init__(self, record_id: str) -> None:
        super().__init__(record_id)
        self.general: Optional[GeneralResults] = None
        self.subcluster: Optional[GeneralResults] = None
        self.knowncluster: Optional[GeneralResults] = None

    def to_json(self) -> Dict[str, Any]:
        assert self.general or self.subcluster or self.knowncluster
        result = {"schema_version": self.schema_version,
                  "record_id": self.record_id}
        for name, subresult in [("general", self.general),
                                ("subcluster", self.subcluster),
                                ("knowncluster", self.knowncluster)]:
            if subresult and subresult.region_results:
                result[name] = subresult.to_json()
        return result

    @staticmethod
    def from_json(json: Dict[str, Any], record: Record) -> "ClusterBlastResults":
        if json["schema_version"] != ClusterBlastResults.schema_version:
            raise ValueError(f"Incompatible results schema version, expected {ClusterBlastResults.schema_version}")
        results = ClusterBlastResults(json["record_id"])
        for attr in ["general", "subcluster", "knowncluster"]:
            if attr in json:
                logging.debug("Regenerating clusterblast subresults: %s", attr)
                setattr(results, attr, GeneralResults.from_json(json[attr], record))
        assert results.general or results.subcluster or results.knowncluster
        return results

    def add_to_record(self, record: Record) -> None:
        for result in [self.general, self.subcluster, self.knowncluster]:
            if result is not None:
                result.add_to_record(record)

    def write_outputs(self, record: Record, options: ConfigType) -> None:
        for subresult in [self.general, self.knowncluster, self.subcluster]:
            if subresult:
                subresult.write_to_file(record, options)


def write_clusterblast_output(options: ConfigType, record: Record,
                              cluster_result: RegionResult, proteins: Dict[str, Protein],
                              searchtype: str = "clusterblast") -> None:
    """ Writes a text form of clusterblast results to file.

        Arguments:
            options: the antismash config
            record: the record that the results came from
            cluster_result: the RegionResult object to write information about
            proteins: a dict mapping protein name to Protein
            searchtype: the name of the module of which to write results for

        Returns:
            None
    """
    assert isinstance(proteins, dict)

    region_number = cluster_result.region.get_region_number()
    filename = f"{record.id}_c{region_number}.txt"

    with changed_directory(_get_output_dir(options, searchtype)):
        with open(filename, "w", encoding="utf-8") as out_file:
            _write_output(out_file, record, cluster_result, proteins)


def _write_output(out_file: IO, record: Record, cluster_result: RegionResult,
                  proteins: Dict[str, Protein]) -> None:
    ranking = cluster_result.ranking
    # Output for each hit: table of genes and locations of input cluster,
    # table of genes and locations of hit cluster, table of hits between the clusters
    out_file.write("ClusterBlast scores for " + record.id + "\n")
    out_file.write("\nTable of genes, locations, strands and annotations of query cluster:\n")
    for i, cds in enumerate(cluster_result.region.cds_children):
        if cds.strand == 1:
            strand = "+"
        else:
            strand = "-"
        out_file.write("\t".join([cds.get_name(), str(int(cds.location.start)),
                                  str(int(cds.location.end)), strand, cds.product]) + "\t\n")
    out_file.write("\n\nSignificant hits: \n")
    for i, cluster_and_score in enumerate(ranking):
        cluster = cluster_and_score[0]
        out_file.write(f"{i + 1}. {cluster.accession}\t{cluster.description}\n")

    out_file.write("\n\nDetails:")
    for i, cluster_and_score in enumerate(ranking):
        cluster, score = cluster_and_score
        nrhits = score.hits
        out_file.write("\n\n>>\n")
        out_file.write(f"{i + 1}. {cluster.accession}\n")
        out_file.write(f"Source: {cluster.description}\n")
        out_file.write(f"Type: {cluster.cluster_type}\n")
        out_file.write(f"Number of proteins with BLAST hits to this cluster: {nrhits}\n")
        out_file.write(f"Cumulative BLAST score: {score.blast_score}\n\n")
        out_file.write("Table of genes, locations, strands and annotations of subject cluster:\n")
        for protein_name in cluster.proteins:
            protein = proteins.get(protein_name)
            if protein:
                out_file.write(str(protein))
        out_file.write("\nTable of Blast hits (query gene, subject gene,"
                       " %identity, blast score, %coverage, e-value):\n")
        if score.scored_pairings:
            for query, subject in score.scored_pairings:
                out_file.write(f"{query.id}\t{subject.get_table_string()}\n")
        else:
            out_file.write("data not found\n")
        out_file.write("\n")


def _get_output_dir(options: ConfigType, searchtype: str) -> str:
    assert searchtype in ["clusterblast", "subclusterblast", "knownclusterblast"], searchtype
    output_dir = os.path.join(options.output_dir, searchtype)

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    if not os.path.isdir(output_dir):
        raise RuntimeError(f"{output_dir} exists as a file, but must be a directory")
    return output_dir
