# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

from collections import OrderedDict
import logging
import os
from typing import Dict, List

from antismash.common.module_results import ModuleResults

from .data_structures import Score, Query, Subject, ReferenceCluster, Protein, MibigEntry
from .svg_builder import ClusterSVGBuilder

_CLUSTER_LIMIT = 50

def get_result_limit():
    return _CLUSTER_LIMIT

class ClusterResult:
    """ Stores results for a specific cluster in a record, for a particular
        flavour of clusterblast.
    """
    __slots__ = ["cluster", "ranking", "total_hits", "svg_builder", "prefix"]
    def __init__(self, cluster, ranking, reference_proteins, prefix):
        """ Arguments:
                cluster: the cluster feature
                ranking: a list of tuples in the form (ReferenceCluster, Score)
                reference_proteins: used to generate details for SVG output,
                                    only relevant portions are stored
                prefix: an identifier for use in marking SVGs such that
                        javascript on the results page can differentiate between
                        types of clusterblast
        """
        self.cluster = cluster # Cluster
        self.ranking = ranking[:get_result_limit()] # [(ReferenceCluster, Score),...]
        self.total_hits = len(ranking)
        self.svg_builder = ClusterSVGBuilder(cluster, ranking, reference_proteins, prefix)
        self.prefix = prefix

    def update_cluster_descriptions(self, search_type):
        if search_type != "knownclusterblast": #TODO clean this up
            setattr(self.cluster, search_type, self.svg_builder.get_cluster_descriptions())
            return
        self.cluster.knownclusterblast = list(zip(self.svg_builder.get_cluster_descriptions(),
                                 self.svg_builder.get_cluster_accessions()))

    def jsonify(self) -> Dict:
        """ Convert the object into a simple dictionary for use in storing
            results.

            The function from_json() should reconstruct a new and equal
            ClusterResult from the results of this function.

            Returns:
                a dict containing the object data in basic types
        """
        result = {}
        result["cluster_number"] = self.cluster.get_cluster_number()
        result["total_hits"] = self.total_hits
        ranking = []
        for cluster, score in self.ranking:
            scoring = {key : getattr(score, key) for key in score.__slots__ if key != "scored_pairings"}
            json_cluster = {key : getattr(cluster, key) for key in cluster.__slots__}
            scoring["pairings"] = [(query.entry, query.index, vars(subject)) for query, subject in score.scored_pairings]
            ranking.append((json_cluster, scoring))
        result["ranking"] = ranking
        result["prefix"] = self.prefix
        return result

    @staticmethod
    def from_json(json, record, reference_proteins):
        """ Convert a simple dictionary into a new ClusterResult object.

            The function ClusterResult.jsonify() should reconstruct the data
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
            for key, val in details.items():
                if key == "pairings":
                    continue
                setattr(score, key, val)
            for pairing in pairings: # entry, index, dict of subject
                query = Query(pairing[0], pairing[1])
                subject = Subject.from_dict(pairing[2])
                score.scored_pairings.append((query, subject))
            ranking.append((ref_cluster, score))

        cluster = record.get_cluster(json["cluster_number"])
        result = ClusterResult(cluster, ranking, reference_proteins, json["prefix"])
        result.total_hits = json["total_hits"]
        return result

    def write_svg_files(self, svg_dir, prefix) -> List[str]:
        """ Write all generated SVG files, one overview SVG and one for each
            ReferenceCluster pairing. The overview SVG will have _all in
            the filename and the individual pairings will have a counter.

            Arguments:
                svg_dir: the directory to save the files in
                prefix: the file prefix to use (e.g. 'knownclusterblast')

            Returns:
                a list of filenames created
        """
        cluster_num = self.cluster.get_cluster_number()

        filename = "%s%d_all.svg" % (prefix, cluster_num)
        with open(os.path.join(svg_dir, filename), "w") as handle:
            handle.write(self.svg_builder.get_overview_contents(width=800, height=50 + 50 * len(self.svg_builder.hits)))

        files = []
        for i in range(len(self.svg_builder.hits)):
            filename = "%s%d_%d.svg" % (prefix, cluster_num, i + 1) # 1-index
            with open(os.path.join(svg_dir, filename), "w") as handle:
                handle.write(self.svg_builder.get_pairing_contents(i, width=800, height=230))
            files.append(filename)
        return files

class GeneralResults(ModuleResults):
    schema_version = 1
    def __init__(self, record_id, search_type="clusterblast") -> None:
        assert record_id and isinstance(record_id, str)
        assert search_type and isinstance(search_type, str)
        super().__init__(record_id)
        self.cluster_results = [] # [ClusterResult, ...]
        self.search_type = search_type
        # keep here instead of duplicating in clusters
        # and only keep those that are relevant instead of 7 million
        self.proteins_of_interest = OrderedDict() # protein name -> Protein
        self.clusters_of_interest = OrderedDict() # cluster name -> ReferenceCluster
        self.mibig_entries = None

    def add_cluster_result(self, result, reference_clusters, reference_proteins):
        assert isinstance(result, ClusterResult)
        assert reference_clusters and reference_proteins # {str:ReferenceCluster}, {str:Protein}
        self.cluster_results.append(result)
        # keep data on proteins relevant to this cluster
        for cluster, _ in result.ranking:
            cluster_label = cluster.get_name()
            protein_names = reference_clusters[cluster_label].proteins
            self.clusters_of_interest[cluster_label] = cluster
            for protein_name in protein_names:
                self.proteins_of_interest[protein_name] = reference_proteins[protein_name]

    def write_to_file(self, record, options):
        for cluster in self.cluster_results:
            write_clusterblast_output(options, record, cluster,
                    self.proteins_of_interest, searchtype=self.search_type)

    def write_svg_files(self, svg_dir):
        for cluster_result in self.cluster_results:
            cluster_result.write_svg_files(svg_dir, self.search_type)

    def to_json(self):
        if not self.cluster_results:
            return None
        data = {"record_id" : self.record_id,
                "schema_version" : self.schema_version,
                "results" : [res.jsonify() for res in self.cluster_results],
                "proteins" : [{key : getattr(protein, key) for key in protein.__slots__} for protein in self.proteins_of_interest.values()],
                "search_type" : self.search_type}
        if self.mibig_entries is not None:
            entries = {}
            for cluster_number, proteins in self.mibig_entries.items():
                cluster = {}
                for protein, protein_entries in proteins.items():
                    cluster[protein] = [protein_entry.values for protein_entry in protein_entries]
                entries[cluster_number] = cluster
            data["mibig_entries"] = entries
        return data

    def add_to_record(self):
        for cluster_result in self.cluster_results:
            cluster_result.update_cluster_descriptions(self.search_type)

    @staticmethod
    def from_json(json, record):
        if json["schema_version"] != GeneralResults.schema_version:
            raise ValueError("Incompatible results schema version, expected %d" \
                    % GeneralResults.schema_version)
        assert record.id == json["record_id"]
        result = GeneralResults(json["record_id"], search_type=json["search_type"])
        for prot in json["proteins"]:
            protein = Protein(prot["name"], prot["locus_tag"], prot["location"],
                              prot["strand"], prot["annotations"])
            result.proteins_of_interest[protein.name] = protein
        for cluster_result in json["results"]:
            result.cluster_results.append(ClusterResult.from_json(cluster_result,
                                           record, result.proteins_of_interest))
        if "mibig_entries" in json:
            entries = {}
            for cluster_number, proteins in json["mibig_entries"].items():
                entries[int(cluster_number)] = {}
                for protein, protein_entries in proteins.items():
                    entries[int(cluster_number)][protein] = [MibigEntry(*entry) for entry in protein_entries]
            result.mibig_entries = entries
        return result

class ClusterBlastResults(ModuleResults):
    schema_version = 1
    def __init__(self, record_id):
        super().__init__(record_id)
        self.general = None
        self.subcluster = None
        self.knowncluster = None

    def to_json(self):
        assert self.general or self.subcluster or self.knowncluster
        result = {"schema_version" : self.schema_version,
                  "record_id" : self.record_id}
        for attr in ["general", "subcluster", "knowncluster"]:
            val = getattr(self, attr)
            if val and val.cluster_results:
                result[attr] = val.to_json()
        return result

    @staticmethod
    def from_json(json, record):
        if json["schema_version"] != ClusterBlastResults.schema_version:
            raise ValueError("Incompatible results schema version, expected %d" \
                    % ClusterBlastResults.schema_version)
        results = ClusterBlastResults(json["record_id"])
        for attr in ["general", "subcluster", "knowncluster"]:
            if attr in json:
                logging.debug("Regenerating clusterblast subresults: %s", attr)
                setattr(results, attr, GeneralResults.from_json(json[attr], record))
        assert results.general or results.subcluster or results.knowncluster
        return results

    def add_to_record(self, record):
        for result in [self.general, self.subcluster, self.knowncluster]:
            if result is not None:
                result.add_to_record()

    def write_svg_files(self, svg_dir):
        for result in [self.general, self.subcluster, self.knowncluster]:
            if result:
                result.write_svg_files(svg_dir)

def write_clusterblast_output(options, seq_record, cluster_result, proteins,
                              searchtype="clusterblast"):
    assert isinstance(proteins, dict)

    cluster_number = cluster_result.cluster.get_cluster_number()
    ranking = cluster_result.ranking

    #Output for each hit: table of genes and locations of input cluster, table of genes and locations of hit cluster, table of hits between the clusters
    currentdir = os.getcwd()
    os.chdir(_get_output_dir(options, searchtype))

    out_file = open("cluster" + str(cluster_number) + ".txt", "w")
    out_file.write("ClusterBlast scores for " + seq_record.id + "\n")
    out_file.write("\nTable of genes, locations, strands and annotations of query cluster:\n")
    for i, cds in enumerate(cluster_result.cluster.cds_children):
        if cds.strand == 1:
            strand = "+"
        else:
            strand = "-"
        out_file.write("\t".join([cds.get_accession(), str(int(cds.location.start)), str(int(cds.location.end)), strand, cds.product]) + "\t\n")
    out_file.write("\n\nSignificant hits: \n")
    for i, cluster_and_score in enumerate(ranking):
        cluster = cluster_and_score[0]
        out_file.write("{}. {}\t{}\n".format(i + 1, cluster.accession, cluster.description))

    out_file.write("\n\nDetails:")
    for i, cluster_and_score in enumerate(ranking):
        cluster, score = cluster_and_score
        nrhits = score.hits
        out_file.write("\n\n>>\n")
        out_file.write("{}. {}\n".format(i + 1, cluster.accession))
        out_file.write("Source: {}\n".format(cluster.description))
        out_file.write("Type: {}\n".format(cluster.cluster_type))
        out_file.write("Number of proteins with BLAST hits to this cluster: %d\n" % nrhits)
        out_file.write("Cumulative BLAST score: %d\n\n" % score.blast_score)
        out_file.write("Table of genes, locations, strands and annotations of subject cluster:\n")
        for protein_name in cluster.proteins:
            protein = proteins.get(protein_name)
            if protein:
                out_file.write(str(protein))
        out_file.write("\nTable of Blast hits (query gene, subject gene, %identity, blast score, %coverage, e-value):\n")
        if score.scored_pairings:
            for query, subject in score.scored_pairings:
                out_file.write("{}\t{}\n".format(query.id, subject.get_table_string()))
        else:
            out_file.write("data not found\n")
        out_file.write("\n")
    out_file.close()
    os.chdir(currentdir)

def _get_output_dir(options, searchtype):
    assert searchtype in ["clusterblast", "subclusterblast", "knownclusterblast"], searchtype
    output_dir = os.path.join(options.output_dir, searchtype)

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    if not os.path.isdir(output_dir):
        raise RuntimeError("%s exists as a file, but required as a directory" % output_dir)
    return output_dir
