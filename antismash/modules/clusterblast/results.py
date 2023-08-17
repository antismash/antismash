# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Contains the results class for each variant of clusterblast """


from collections import OrderedDict
import logging
import os
from typing import Any, Dict, List, Optional, Tuple

from antismash.common.module_results import ModuleResults
from antismash.common.path import changed_directory
from antismash.common.secmet import Record, Region
from antismash.config import ConfigType, get_config

from .data_structures import Score, Query, Subject, ReferenceCluster, Protein, MibigEntry
from .svg_builder import ClusterSVGBuilder

_CLUSTER_LIMIT = 50


def get_result_limit() -> int:
    """ Returns the hard limit of cluster matches to report """
    return _CLUSTER_LIMIT


class KnownHitSummary:
    """ Stores some information about a hit from known-CB in a handy to access way """
    def __init__(self, bgc_id: str, name: str, cluster_number: str, similarity: int,
                 cluster_type: str) -> None:
        self.bgc_id = str(bgc_id)
        self.name = str(name)
        self.cluster_number = cluster_number
        self.similarity = int(similarity)
        self.cluster_type = cluster_type


class RegionResult:
    """ Stores results for a specific cluster in a record, for a particular
        flavour of clusterblast.
    """
    __slots__ = ["region", "ranking", "total_hits", "svg_builder", "prefix"]

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
        self.region = region
        self.ranking = ranking[:get_result_limit()]  # [(ReferenceCluster, Score),...]
        self.total_hits = len(ranking)
        self.prefix = prefix
        # for the SVG portion, limit the ranking to the display limit
        display_limit = get_config().cb_nclusters
        # omitting any self-hits in the display
        display_ranking = self.ranking[:display_limit]
        if prefix != "subclusterblast":
            record_prefix = (region.parent_record.original_id or region.parent_record.id).split(".", 1)[0]
            display_ranking = list(filter(lambda pair: pair[0].accession != record_prefix, display_ranking))
            if len(display_ranking) < display_limit < len(self.ranking) - 1:
                display_ranking.append(self.ranking[display_limit])
        assert len(display_ranking) <= display_limit
        self.svg_builder = ClusterSVGBuilder(region, display_ranking, reference_proteins, prefix)

    def update_cluster_descriptions(self, search_type: str) -> None:
        """ Rebuilds the cluster's result descriptions.
            For knownclusterblast, this includes the accessions of the hits.
        """
        if search_type != "knownclusterblast":  # TODO clean this up
            setattr(self.region, search_type, self.svg_builder.get_cluster_descriptions())
            return
        hits = []
        for cluster in self.svg_builder.hits:
            hits.append(KnownHitSummary(cluster.accession, cluster.description,
                                        cluster.ref_cluster_number,
                                        cluster.similarity, cluster.cluster_type))
        self.region.knownclusterblast = hits

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
            scoring = {key: getattr(score, key) for key in score.__slots__ if key != "scored_pairings"}
            json_cluster = {key: getattr(cluster, key) for key in cluster.__slots__}
            scoring["pairings"] = [(query.entry, query.index, vars(subject))
                                   for query, subject in score.scored_pairings]
            # add the similarity score, since it isn't held in the other data structures
            unique_hit_count = len({subject.name for _, subject in score.scored_pairings})
            scoring["similarity"] = int(100 * unique_hit_count / len(cluster.proteins))
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
            # similarity isn't stored in memory, it's only in the JSON
            assert not hasattr(score, "similarity")  # to make sure the above comment stays accurate
            details.pop("similarity")
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

    def write_svg_files(self, svg_dir: str, prefix: str) -> List[str]:
        """ Write all generated SVG files, one overview SVG and one for each
            ReferenceCluster pairing. The overview SVG will have _all in
            the filename and the individual pairings will have a counter.

            Arguments:
                svg_dir: the directory to save the files in
                prefix: the file prefix to use (e.g. 'knownclusterblast')

            Returns:
                a list of filenames created
        """
        region_num = self.region.get_region_number()
        record_index = self.region.parent_record.record_index
        filename = f"{prefix}_r{record_index}c{region_num}_all.svg"
        with open(os.path.join(svg_dir, filename), "w", encoding="utf-8") as handle:
            handle.write(self.svg_builder.get_overview_contents(width=800, height=50 + 50 * len(self.svg_builder.hits)))

        files = []
        for i in range(len(self.svg_builder.hits)):
            filename = f"{prefix}_r{record_index}c{region_num}_{i + 1}.svg"  # 1-indexed
            with open(os.path.join(svg_dir, filename), "w", encoding="utf-8") as handle:
                handle.write(self.svg_builder.get_pairing_contents(i, width=800, height=230))
            files.append(filename)
        return files


class GeneralResults(ModuleResults):
    """ A variant-agnostic results class for clusterblast variants """
    schema_version = 3

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
                self.proteins_of_interest[protein_name] = reference_proteins[protein_name]

    def write_to_file(self, record: Record, options: ConfigType) -> None:
        """ Write the results to a text file """
        for cluster in self.region_results:
            write_clusterblast_output(options, record, cluster,
                                      self.proteins_of_interest,
                                      searchtype=self.search_type)

    def write_svg_files(self, svg_dir: str) -> None:
        """ Write the SVG files for each cluster with results """
        for cluster_result in self.region_results:
            cluster_result.write_svg_files(svg_dir, self.search_type)

    def to_json(self) -> Dict[str, Any]:
        if not self.region_results:
            return {}
        data = {"record_id": self.record_id,
                "schema_version": self.schema_version,
                "results": [res.jsonify() for res in self.region_results],
                "proteins": [{key: getattr(protein, key) for key in protein.__slots__}
                             for protein in self.proteins_of_interest.values()],
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
        for cluster_result in self.region_results:
            cluster_result.update_cluster_descriptions(self.search_type)

    @staticmethod
    def from_json(json: Dict[str, Any], record: Record) -> "GeneralResults":
        current = GeneralResults.schema_version
        # since 2 only added an optional data version, 2 and 1 are kind of compatible
        # so if there's a mismatch of version, but it's these two, let them through
        special_case = current == 2 and json["schema_version"] == 1
        if json["schema_version"] != current and not special_case:
            raise ValueError(f"Incompatible results schema version, expected {GeneralResults.schema_version}")
        assert record.id == json["record_id"]
        data_version = json.get("data_version")
        result = GeneralResults(json["record_id"], search_type=json["search_type"],
                                data_version=data_version)
        for prot in json["proteins"]:
            protein = Protein(prot["name"], prot["locus_tag"], prot["location"],
                              prot["strand"], prot["annotations"])
            result.proteins_of_interest[protein.locus_tag] = protein
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
        self.internal_homology_groups: Dict[int, List[List[str]]] = {}

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

    def write_svg_files(self, svg_dir: str) -> None:
        """ Write the SVG files for each clusterblast variant if it exists """
        for result in [self.general, self.subcluster, self.knowncluster]:
            if result:
                result.write_svg_files(svg_dir)


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
        _write_output(filename, record, cluster_result, proteins)


def _write_output(filename: str, record: Record, cluster_result: RegionResult,
                  proteins: Dict[str, Protein]) -> None:
    ranking = cluster_result.ranking
    # Output for each hit: table of genes and locations of input cluster,
    # table of genes and locations of hit cluster, table of hits between the clusters
    out_file = open(filename, "w", encoding="utf-8")  # pylint: disable=consider-using-with
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
    out_file.close()


def _get_output_dir(options: ConfigType, searchtype: str) -> str:
    assert searchtype in ["clusterblast", "subclusterblast", "knownclusterblast"], searchtype
    output_dir = os.path.join(options.output_dir, searchtype)

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    if not os.path.isdir(output_dir):
        raise RuntimeError(f"{output_dir} exists as a file, but must be a directory")
    return output_dir
