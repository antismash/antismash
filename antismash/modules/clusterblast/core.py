# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

import logging
import os
from collections import OrderedDict

import antismash.common.deprecated as utils
from antismash.common.module_results import ModuleResults
import antismash.common.path as path
import antismash.common.subprocessing as subprocessing
from antismash.config.args import Config
from helperlibs.wrappers.io import TemporaryDirectory

_CLUSTER_LIMIT = 50

def get_result_limit():
    return _CLUSTER_LIMIT

class ClusterResult:
    __slots__ = ["cluster", "ranking", "total_hits"]
    def __init__(self, cluster, ranking):
        self.cluster = cluster # Cluster
        self.ranking = ranking[:get_result_limit()] # [(str, Score),...]
        self.total_hits = len(ranking)

    def jsonify(self):
        result = {}
        result["cluster_number"] = self.cluster.get_cluster_number()
        result["total_hits"] = self.total_hits
        ranking = []
        for cluster_name, score in self.ranking:
            scoring = {key : getattr(score, key) for key in score.__slots__ if key != "scored_pairings"}
            scoring["pairings"] = [(query.entry, query.index, vars(subject)) for query, subject in score.scored_pairings]
            ranking.append((cluster_name, scoring))
        result["ranking"] = ranking
        return result

    @staticmethod
    def from_json(json, record):
        ranking = []
        for cluster_name, details in json["ranking"]:
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
            ranking.append((cluster_name, score))

        cluster = record.get_cluster(json["cluster_number"])
        result = ClusterResult(cluster, ranking)
        result.total_hits = json["total_hits"]
        return result


class GeneralResults(ModuleResults):
    schema_version = 1
    def __init__(self, record_id, search_type="clusterblast"):
        super().__init__(record_id)
        self.cluster_results = []
        self.search_type = search_type
        # keep here instead of duplicating in clusters
        # and only keep those that are relevant instead of 7 million
        self.proteins_of_interest = {}

    def add_cluster_result(self, result, reference_clusters, reference_proteins):
        assert isinstance(result, ClusterResult)
        self.cluster_results.append(result)
        # keep data on proteins relevant to this cluster
        for cluster_label, _ in result.ranking:
            for protein_name in reference_clusters[cluster_label][0]:
                self.proteins_of_interest[protein_name] = reference_proteins.get(protein_name)

    def write_to_file(self, record, reference_clusters, options):
        for cluster in self.cluster_results:
            write_clusterblast_output(options, record, cluster,
                    reference_clusters, self.proteins_of_interest,
                    searchtype=self.search_type)

    def to_json(self):
        if not self.cluster_results:
            return None
        return {"record_id" : self.record_id,
                "schema_version" : self.schema_version,
                "results" : [res.jsonify() for res in self.cluster_results]}

    @staticmethod
    def from_json(json, record):
        if json["schema_version"] != GeneralResults.schema_version:
            raise ValueError("Incompatible results schema version, expected %d" \
                    % GeneralResults.schema_version)
        result = GeneralResults(json["record_id"])
        for cluster_result in json["results"]:
            result.cluster_results.append(ClusterResult.from_json(cluster_result, record))
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
                setattr(results, attr, GeneralResults.from_json(json[attr], record))
        assert results.general or results.subcluster or results.knowncluster
        return results

    def add_to_record(self, record):
        for result in [self.general, self.subcluster, self.knowncluster]:
            result.add_to_record(record)

def runblast(query, target):
    command = ["blastp", "-db", target, "-query", query, "-outfmt", "6",
               "-max_target_seqs", "10000", "-evalue", "1e-05",
               "-out", query.split(".")[0] + ".out"]
    res = subprocessing.execute(command)
    if res.return_code:
        raise RuntimeError("blastp run failed: %s..." % res.stderr[-200:])


def run_diamond(query, target, tempdir, options):
    command = [
        "diamond", "blastp",
        "--db", target,
        "--threads", str(options.cpus),
        "--query", query,
        "--compress", "0",
        "--max-target-seqs", "10000",
        "--evalue", "1e-05",
        "--out", "input.out",
        "--outfmt", "6",  # 6 is blast tabular format, just as in blastp
        "--tmpdir", tempdir
    ]
    return subprocessing.execute(command)


def make_blastdb(inputfile, dbname):
    command = ["makeblastdb", "-in", inputfile, "-out", dbname, "-dbtype", "prot"]
    subprocessing.execute(command)

def load_geneclusters(searchtype):
    #Load gene cluster database into memory
    options = Config()
    clusterblastdir = os.path.join(options.database_dir, 'clusterblast')
    subclusterblastdir = os.path.join(path.get_full_path(__file__, "data"), "sub")
    knownclusterblastdir = os.path.join(path.get_full_path(__file__, "data"), "known")

    if searchtype == "clusterblast":
        logging.info("ClusterBlast: Loading gene clusters database into memory...")
        geneclustersfile = os.path.join(clusterblastdir, "geneclusters.txt")
    elif searchtype == "subclusterblast":
        logging.info("SubClusterBlast: Loading gene clusters database into memory...")
        geneclustersfile = os.path.join(subclusterblastdir, "subclusters.txt")
    elif searchtype == "knownclusterblast":
        logging.info("KnownClusterBlast: Loading gene clusters database into memory...")
        geneclustersfile = os.path.join(knownclusterblastdir, "knownclusters.txt")
    geneclustersfile = open(geneclustersfile, "r")
    filetext = geneclustersfile.read()
    lines = [line for line in filetext.split("\n") if "\t" in line]
    clusters = {}
    for i in lines:
        tabs = i.split("\t")
        accession = tabs[0]
        clusterdescription = tabs[1]
        clusternr = tabs[2]
        clustertype = tabs[3]
        clustername = accession + "_" + clusternr
        clustertags = tabs[4].split(";")
        clusterprots = tabs[5].split(";")
        clusters[clustername] = [clusterprots, clusterdescription, clustertype, clustertags]
    return clusters

class Protein(object):
    """ Holds details of a protein. Members cannot be added dynamically.
    """
    # At time of writing this class will be instantiated ~7 million times per
    # clusterblast invocation.
    # With those numbers, the memory use of the class without slots: 2.3 Gb
    #                                                and with slots: 0.4 Gb
    __slots__ = ("name", "locus_tag", "location", "strand", "annotations")
    def __init__(self, name, locus_tag, location, strand, annotations):
        self.name = name
        self.locus_tag = locus_tag
        self.location = location
        self.strand = strand
        self.annotations = annotations

    def __str__(self):
        if len(self.location.split("-")) != 2:
            raise ValueError("Invalid location in Protein: %s"%self.location)
        tag = self.locus_tag
        if tag == "no_locus_tag":
            tag = self.name
        locations = self.location.replace("-", "\t")
        return "{}\t{}\t{}\t{}\t{}\n".format(tag, self.name, locations,
                                                 self.strand, self.annotations)

def load_geneclusterproteins(accessions, searchtype):
    options = Config()
    clusterblastdir = os.path.join(options.database_dir, 'clusterblast')
    subclusterblastdir = os.path.join(path.get_full_path(__file__, "data"), "sub")
    knownclusterblastdir = os.path.join(path.get_full_path(__file__, "data"), "known")
    #Load gene cluster database proteins info into memory
    if searchtype == "clusterblast":
        logging.info("ClusterBlast: Loading gene cluster database proteins into " \
        "memory...")
        gclusterprotsfile = os.path.join(clusterblastdir, "geneclusterprots.fasta")
    elif searchtype == "subclusterblast":
        logging.info("SubClusterBlast: Loading gene cluster database proteins into " \
        "memory...")
        gclusterprotsfile = os.path.join(subclusterblastdir, "subclusterprots.fasta")
    elif searchtype == "knownclusterblast":
        logging.info("KnownClusterBlast: Loading gene cluster database proteins into " \
        "memory...")
        gclusterprotsfile = os.path.join(knownclusterblastdir, "knownclusterprots.fasta")

    proteins = {}

    with open(gclusterprotsfile, 'r') as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line or line[0] != ">":
                continue
            tabs = line.split("|")
            locustag = tabs[4]
            if locustag in accessions:
                locustag = "h_" + locustag
            location = tabs[2]
            strand = tabs[3]
            annotations = tabs[5]
            name = tabs[6]
            proteins[name] = Protein(name, locustag, location, strand, annotations)
    return proteins

def load_clusterblast_database(seq_record, searchtype="clusterblast"):
    accessions = set()
    for cds in seq_record.get_cds_features():
        acc = cds.get_accession()
        accessions.add(acc)
    clusters = load_geneclusters(searchtype)
    proteins = load_geneclusterproteins(accessions, searchtype)
    return clusters, proteins

def create_blast_inputs(cluster):
    #Create input fasta files for BLAST search
    names = []
    seqs = []
    for cds in cluster.cds_children:
        if cds.strand == 1:
            strand = "+"
        else:
            strand = "-"
        fullname = "|".join(["input", "c" + str(cluster.get_cluster_number()), \
                             str(int(cds.location.start)) + "-" + \
                             str(int(cds.location.end)), \
                             strand, cds.get_accession(), cds.product])
        names.append(fullname)
        seqs.append(str(utils.get_aa_sequence(cds)))

    return names, seqs

def run_internal_blastsearch():
    #Run and parse BLAST search
    make_blastdb("internal_input.fasta", "internal_input.fasta")
    runblast("internal_input.fasta", "internal_input.fasta")
    blastoutput = open("internal_input.out", "r").read()
    return blastoutput

def uniqueblasthitfilter(blastlines):
    #Filter for best blast hits (of one query on each subject)
    query_subject_combinations = set()
    blastlines2 = []
    for tabs in blastlines:
        # if it doesn't even have both values, it's not a hit so skip it
        if len(tabs) < 2:
            continue
        query = tabs[0]
        subject = tabs[1]
        query_subject_combination = (query, subject)
        if query_subject_combination not in query_subject_combinations:
            query_subject_combinations.add(query_subject_combination)
            blastlines2.append(tabs)
    return blastlines2

class Subject:
    def __init__(self, name, genecluster, start, end, strand, annotation,
                 perc_ident, blastscore, perc_coverage, evalue, locus_tag):
        self.name = name
        self.genecluster = genecluster
        self.start = start
        self.end = end
        self.strand = strand
        self.annotation = annotation
        self.perc_ident = int(perc_ident)
        self.blastscore = int(blastscore)
        self.perc_coverage = float(perc_coverage)
        self.evalue = float(evalue)
        self.locus_tag = locus_tag

    def get_table_string(self):
        return "\t".join([str(i) for i in [self.name, self.perc_ident,
                                           self.blastscore, self.perc_coverage,
                                           self.evalue]])

    @staticmethod
    def from_dict(data):
        args = []
        for key in ["name", "genecluster", "start", "end", "strand",
                    "annotation", "perc_ident", "blastscore", "perc_coverage",
                    "evalue", "locus_tag"]:
            args.append(data[key])
# stop pylint throwing errors because it can't work out the splat operator
# pylint: disable=no-value-for-parameter
        return Subject(*args)
# pylint: enable=no-value-for-parameter


class Query:
    def __init__(self, entry, index):
        parts = entry.split("|")
        self.cluster_number = int(parts[1][1:]) # c1 -> 1
        self.id = parts[4]
        self.entry = entry
        self.subjects = OrderedDict()
        self.cluster_name_to_subjects = {}
        self.index = index

    def add_subject(self, subject):
        self.subjects[subject.name] = subject
        if subject.genecluster not in self.cluster_name_to_subjects:
            self.cluster_name_to_subjects[subject.genecluster] = []
        self.cluster_name_to_subjects[subject.genecluster].append(subject)

    def get_subjects_by_cluster(self, cluster_name):
        return self.cluster_name_to_subjects.get(cluster_name, [])

def parse_subject(tabs, seqlengths, geneclustergenes, seq_record):
    if len(tabs) < 12:
        logging.error("Malformed blast pairing: %s", "\t".join(tabs))
    query = tabs[0]
    subject_parts = tabs[1].split("|")
    subject = subject_parts[4]
    if subject == "no_locus_tag":
        subject = subject_parts[6]
    if subject in geneclustergenes:
        subject = "h_" + subject
    if len(subject_parts) > 6:
        locustag = subject_parts[6]
    else:
        locustag = ""
    genecluster = "{}_{}".format(subject_parts[0], subject_parts[1])
    start, end = subject_parts[2].split("-")[:2]
    strand = subject_parts[3]
    annotation = subject_parts[5]
    perc_ident = int(float(tabs[2]) + 0.5)
    evalue = str(tabs[10])
    blastscore = int(float(tabs[11])+0.5)
    query_key = query.split("|")[4]
    if query_key in seqlengths:
        perc_coverage = (float(tabs[3]) / seqlengths[query_key]) * 100
    else:
        feature_by_id = utils.get_feature_dict_protein_id(seq_record)
        seqlength = len(utils.get_aa_sequence(feature_by_id[query_key]))
        perc_coverage = (float(tabs[3]) / seqlength) * 100
    return Subject(subject, genecluster, start, end, strand, annotation,
                   perc_ident, blastscore, perc_coverage, evalue, locustag)

def parse_all_clusters(blasttext, minseqcoverage, minpercidentity, seq_record):
    """ Parses blast results, groups into results by cluster number

        blasttext: the output from diamond in blast format
        minseqcoverage: the exclusive lower bound of sequence coverage for a match
        minpercidentity: the exclusive lower bound of identity similarity for a match
        seq_record: used to get all gene ids in the cluster, and used as a
                backup to fetch sequence length if missing from seqlengths
    """
    seqlengths = get_cds_lengths(seq_record)
    geneclustergenes = [cds.get_accession() for cds in utils.get_cds_features_within_clusters(seq_record)]
    queries = OrderedDict()
    clusters = OrderedDict()
    blastlines = uniqueblasthitfilter([line.split("\t") for line in blasttext.rstrip().split("\n")])
    current_query = None
    queries_by_cluster_number = {}
    clusters_by_query_cluster_number = {}

    for tabs in blastlines:
        query = tabs[0]
        subject = parse_subject(tabs, seqlengths, geneclustergenes, seq_record)

        # only process the pairing if limits met
        if subject.perc_ident <= minpercidentity \
                or subject.perc_coverage <= minseqcoverage:
            continue

        new_query = query not in queries

        if new_query:
            current_query = Query(query, len(queries))
            cluster_number = current_query.cluster_number
            # is it a new cluster number? if so, reset collections
            if cluster_number not in queries_by_cluster_number:
                queries = OrderedDict()
                clusters = OrderedDict()
                # reset query index, since we started a new collection
                current_query.index = 0
                # link them
                queries_by_cluster_number[cluster_number] = queries
                clusters_by_query_cluster_number[cluster_number] = clusters
            # finally, add the query to the current tracker
            queries[query] = current_query

        if subject.genecluster not in clusters:
            clusters[subject.genecluster] = []
        clusters[subject.genecluster].append(current_query)

        # link the subject to the query
        current_query.add_subject(subject)

    return clusters_by_query_cluster_number, queries_by_cluster_number

def blastparse(blasttext, minseqcoverage, minpercidentity, seq_record):
    """ blasttext: the output from diamond in blast format
        minseqcoverage: the exclusive lower bound of sequence coverage for a match
        minpercidentity: the exclusive lower bound of identity similarity for a match
        seq_record: used to get all gene ids in the cluster, and used as a
                backup to fetch sequence length if missing from seqlengths
    """
    seqlengths = get_cds_lengths(seq_record)
    geneclustergenes = [cds.get_accession() for cds in utils.get_cds_features_within_clusters(seq_record)]
    queries = OrderedDict()
    clusters = OrderedDict()
    blastlines = uniqueblasthitfilter([line.split("\t") for line in blasttext.rstrip().split("\n")])
    current_query = None

    for tabs in blastlines:
        query = tabs[0]
        subject = parse_subject(tabs, seqlengths, geneclustergenes, seq_record)

        # only process the pairing if limits met
        if subject.perc_ident <= minpercidentity \
                or subject.perc_coverage <= minseqcoverage:
            continue

        new_query = query not in queries
        new_hit = subject.genecluster not in clusters

        if new_query:
            current_query = Query(query, len(queries))
            queries[query] = current_query

        if new_hit:
            clusters[subject.genecluster] = []
        clusters[subject.genecluster].append(current_query)

        # link the subject to the query
        current_query.add_subject(subject)

    return queries, clusters

def get_cds_lengths(seq_record):
    seqlengths = {}
    for cds in seq_record.get_cds_features():
        seqlength = len(str(utils.get_aa_sequence(cds)))
        seqlengths[cds.get_accession()] = seqlength
    return seqlengths

def find_internal_orthologous_groups(queries, clusternames):
    #find and store internal homologs
    groups = []
    for name in clusternames:
        if name not in queries:
            groups.append([name.split("|")[4]])
            continue
        query = queries[name]
        new_group = []
        for hit in query.subjects:
            if hit.startswith("h_"):
                new_group.append(hit[2:])
            else:
                new_group.append(hit)
        if query.id not in new_group:
            new_group.append(query.id)
        for i, other_group in enumerate(groups):
            for new_member in new_group:
                if new_member in other_group:
                    del groups[i]
                    for member in other_group:
                        if member not in new_group:
                            new_group.append(member)
                    break
        new_group.sort()
        groups.append(new_group)
    return groups

def internal_homology_blast(seq_record):
    """ Run BLAST on gene cluster proteins of each cluster on itself to find
        internal homologs
        store groups of homologs - including singles - in a dictionary
        as a list of lists accordingly
    """
    with TemporaryDirectory(change=True):
        logging.info("Finding internal homologs in each gene cluster..")
        internalhomologygroups = {}
        for genecluster in seq_record.get_clusters():
            cluster_number = genecluster.get_cluster_number()
            iqueryclusternames, iqueryclusterseqs = create_blast_inputs(genecluster)
            utils.writefasta(iqueryclusternames, iqueryclusterseqs, "internal_input.fasta")
            blastoutput = run_internal_blastsearch()
            queries, _ = blastparse(blastoutput, 25, 30, seq_record)
            groups = find_internal_orthologous_groups(queries, iqueryclusternames)
            internalhomologygroups[cluster_number] = groups
    return internalhomologygroups

def read_clusterblast_output(options):
    blastoutput = []
    for i in range(options.cpus):
        with open("input" + str(i) + ".out", "r") as handle:
            output = handle.read()
        blastoutput.append(output)
    return "".join(blastoutput)

def write_raw_clusterblastoutput(output_dir, blast_output, search_type="clusterblast"):
    if search_type == "clusterblast":
        filename = "clusterblastoutput.txt"
    elif search_type == "subclusterblast":
        filename = "subclusterblastoutput.txt"
    elif search_type == "knownclusterblast":
        filename = "knownclusterblastoutput.txt"
    with open(os.path.join(output_dir, filename), "w") as handle:
        handle.write(blast_output)

class Score:
    __slots__ = ("hits", "core_gene_hits", "blast_score", "synteny_score",
                 "core_bonus", "scored_pairings")
    def __init__(self):
        self.hits = 0
        self.core_gene_hits = 0
        self.blast_score = 0
        self.synteny_score = 0
        self.core_bonus = 0
        self.scored_pairings = []

    def sort_score(self):
        """ the pre-existing, unexplained sort weighting """
        if self.core_gene_hits:
            self.core_bonus = 3
        return (self.hits
                + self.core_bonus
                + self.core_gene_hits / 100.
                + self.blast_score / 10e7
                + self.synteny_score)

def parse_clusterblast_dict(queries, clusters, cluster_name, allcoregenes):
    result = Score()
    hitpositions = []
    hitposcorelist = []
    cluster_locii = clusters[cluster_name][0]
    for query in queries:
        querynrhits = 0
        for subject in query.get_subjects_by_cluster(cluster_name):
            assert cluster_name == subject.genecluster
            if subject.locus_tag not in cluster_locii:
                continue
            index_pair = [query.index, cluster_locii.index(subject.locus_tag)]
            if index_pair in hitpositions:
                continue
            querynrhits += 1
            result.blast_score += subject.blastscore
            result.scored_pairings.append([query, subject])
            hitpositions.append(index_pair)
        if querynrhits:
            result.hits += 1
            hit_value = 0
            if query.id in allcoregenes:
                result.core_gene_hits += 1
                hit_value = 1
            hitposcorelist.extend([hit_value]*querynrhits)
    return result, hitpositions, hitposcorelist

def find_clusterblast_hitsgroups(hitpositions):
    #Find groups of hits
    groups = {}
    for pair in hitpositions:
        if pair[0] not in groups:
            groups[pair[0]] = [pair[1]]
        else:
            groups[pair[0]].append(pair[1])
    return groups

def calculate_synteny_score(hitgroupsdict, hitpositions, hitposcorelist):
    scored_queries = set()
    scored_hits = set()
    # Calculate synteny score; give score only if more than one hits
    # (otherwise no synteny possible), and only once for every query gene and every hit gene
    synteny_score = 0
    for i, pos in enumerate(hitpositions[:-1]):
        i += 1 # we want a 1-index, not 0-index
        pair = hitpositions[i]
        # Check if a gene homologous to this gene has already been
        # scored for synteny in the previous entry
        if pos[1] in hitgroupsdict[pair[0]] or pos[0] in scored_queries or pos[1] in scored_hits:
            continue
        if (abs(pos[0] - pair[0]) < 2) and abs(pos[0] - pair[0]) == abs(pos[1] - pair[1]):
            synteny_score += 1
            if hitposcorelist[i - 1] == 1 or hitposcorelist[i] == 1:
                synteny_score += 1
            scored_queries.add(pos[0])
            scored_hits.add(pos[1])
    return synteny_score

def score_clusterblast_output(clusters, allcoregenes, cluster_names_to_queries):
    #Score BLAST output on all gene clusters
    #Rank gene cluster hits based on 1) number of protein hits covering >25% sequence length or at least 100aa alignment, with >30% identity and 2) cumulative blast score
    #Find number of protein hits and cumulative blast score for each gene cluster
    results = {}
    for cluster_name, queries in cluster_names_to_queries.items():
        result, hitpositions, hitposcorelist = parse_clusterblast_dict(queries, clusters, cluster_name, allcoregenes)
        if result.hits <= 1:
            continue
        hitgroupsdict = find_clusterblast_hitsgroups(hitpositions)
        result.synteny_score = calculate_synteny_score(hitgroupsdict, hitpositions, hitposcorelist)
        # ensure at least two different subjects were found
        initial = hitpositions[0][1]
        for _, subject in hitpositions[1:]:
            if subject != initial:
                results[cluster_name] = result
                break
    #Sort gene clusters by score
    return sorted(results.items(), reverse=True, key=lambda x: x[1].sort_score())

def _get_output_dir(options, searchtype):
    assert searchtype in ["clusterblast", "subclusterblast", "knownclusterblast"]
    output_dir = os.path.join(options.output_dir, searchtype)

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    if not os.path.isdir(output_dir):
        raise RuntimeError("%s exists as a file, but required as a directory" % output_dir)
    return output_dir


def write_clusterblast_output(options, seq_record, cluster_result, clusters,
                              proteins, searchtype="clusterblast"):
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
    top_hits = ranking[:100]
    for i, cluster_and_result in enumerate(top_hits):
        cluster = cluster_and_result[0]
        out_file.write("{}. {}\t{}\n".format(i + 1, cluster, clusters[cluster][1]))

    out_file.write("\n\nDetails:")
    for i, cluster_and_result in enumerate(top_hits):
        cluster, result = cluster_and_result
        nrhits = result.hits
        out_file.write("\n\n>>\n")
        out_file.write("{}. {}\n".format(i + 1, cluster))
        out_file.write("Source: {}\n".format(clusters[cluster][1]))
        out_file.write("Type: {}\n".format(clusters[cluster][2]))
        out_file.write("Number of proteins with BLAST hits to this cluster: %d\n" % nrhits)
        out_file.write("Cumulative BLAST score: %d\n\n" % result.blast_score)
        out_file.write("Table of genes, locations, strands and annotations of subject cluster:\n")
        clusterproteins = clusters[cluster][0]
        for protein_name in clusterproteins:
            protein = proteins.get(protein_name)
            if protein:
                out_file.write(str(protein))
        out_file.write("\nTable of Blast hits (query gene, subject gene, %identity, blast score, %coverage, e-value):\n")
        if result.scored_pairings:
            for query, subject in result.scored_pairings:
                out_file.write("{}\t{}\n".format(query.id, subject.get_table_string()))
        else:
            out_file.write("data not found\n")
        out_file.write("\n")
    out_file.close()
    os.chdir(currentdir)
