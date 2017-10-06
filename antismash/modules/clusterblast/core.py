# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

from collections import defaultdict, OrderedDict
import logging
import os
from typing import Dict, List, Set, Tuple # pylint: disable=unused-import

from helperlibs.wrappers.io import TemporaryDirectory

from antismash.common import path, subprocessing, deprecated as utils, secmet
from antismash.config import get_config

from .data_structures import Subject, Query, Protein, ReferenceCluster, Score

def get_core_gene_ids(record) -> Set[str]: #TODO: consider moving into secmet
    """ Fetches all gene accessions of genes with CORE gene function from all
        clusters in a record

        Arguments:
            record: the record to gather gene ids from

        Returns:
            a set containing all core gene names
    """
    cores = []
    for gene in record.get_cds_features_within_clusters():
        if gene.gene_function == secmet.GeneFunction.CORE:
            cores.append(gene.get_accession())
    return cores

def runblast(query, database) -> str:
    """ Runs blastp, comparing the given query with the given database

        An output file will be created, using the name of the query but with the
        extension changed to .out

        Arguments:
            query: the path of query sequence file
            target: the path of the database to compare to

        Returns:
            the name of the created output file
    """
    out_file = query.rsplit(".", 1)[0] + ".out"
    command = ["blastp", "-db", database, "-query", query, "-outfmt", "6",
               "-max_target_seqs", "10000", "-evalue", "1e-05",
               "-out", out_file]
    res = subprocessing.execute(command)
    if not res.successful():
        raise RuntimeError("blastp run failed: %s..." % res.stderr[-200:])
    return out_file

def run_diamond(query, database, tempdir, options) -> str:
    """ Runs diamond, comparing the given query to the given database

        Arguments:
            query: the path of query sequence file
            target: the path of the database to compare to
            tempdir: the path of a temporary directory for diamond to use
            options: antismash Config

        Returns:
            the name of the output file created
    """
    logging.debug("Running external command: diamond")
    command = [
        "diamond", "blastp",
        "--db", database,
        "--threads", str(options.cpus),
        "--query", query,
        "--compress", "0",
        "--max-target-seqs", "10000",
        "--evalue", "1e-05",
        "--out", "input.out",
        "--outfmt", "6",  # 6 is blast tabular format, just as in blastp
        "--tmpdir", tempdir
    ]
    result = subprocessing.execute(command)
    if not result.successful():
        raise RuntimeError("diamond failed to run: %s -> %s" % (command, result.stderr[-100:]))
    return "input.out"


def make_blastdb(inputfile, db_prefix) -> subprocessing.RunResult:
    """ Runs makeblastdb on the inputs to create a blast protein database

        makeblastdb will create 3 files with the given prefix and the extensions:
            .pin, .phr, .psq

        Arguments:
            inputfile: the input filename
            db_prefix: the prefix to use for the created database

        Returns:
            a subprocessing.RunResult instance
    """
    command = ["makeblastdb", "-in", inputfile, "-out", db_prefix, "-dbtype", "prot"]
    result = subprocessing.execute(command)
    if not result.successful():
        raise RuntimeError("makeblastdb failed to run: %s -> %s" % (command, result.stderr[-100:]))
    return result


def load_reference_clusters(searchtype) -> Dict[str, ReferenceCluster]:
    """ Load gene cluster database

        Arguments:
            searchtype: determines which database to use, allowable values:
                            clusterblast, subclusterblast, knownclusterblast

        Returns:
            a dictionary mapping reference cluster name to ReferenceCluster
            instance
    """
    options = get_config()

    if searchtype == "clusterblast":
        logging.info("ClusterBlast: Loading gene clusters database into memory...")
        data_dir = os.path.join(options.database_dir, 'clusterblast')
        reference_cluster_file = os.path.join(data_dir, "geneclusters.txt")
    elif searchtype == "subclusterblast":
        logging.info("SubClusterBlast: Loading gene clusters database into memory...")
        data_dir = path.get_full_path(__file__, "data", "sub")
        reference_cluster_file = os.path.join(data_dir, "subclusters.txt")
    elif searchtype == "knownclusterblast":
        logging.info("KnownClusterBlast: Loading gene clusters database into memory...")
        data_dir = path.get_full_path(__file__, "data", "known")
        reference_cluster_file = os.path.join(data_dir, "knownclusters.txt")
    with open(reference_cluster_file, "r") as handle:
        filetext = handle.read()
    lines = [line for line in filetext.splitlines() if "\t" in line]
    clusters = {}
    for i in lines:
        tabs = i.split("\t")
        accession = tabs[0]
        description = tabs[1]
        cluster_number = tabs[2]
        cluster_type = tabs[3]
        tags = tabs[4].split(";")
        proteins = tabs[5].split(";")
        if not proteins[-1]:
            proteins.pop(-1)
        cluster = ReferenceCluster(accession, cluster_number, proteins,
                                   description, cluster_type, tags)
        clusters[cluster.get_name()] = cluster
    return clusters

def load_reference_proteins(accessions, searchtype) -> Dict[str, Protein]:
    """ Load protein database

        Arguments:
            accessions: a set of all CDS names in the record to avoid collisions
            searchtype: determines which database to use, allowable values:
                            clusterblast, subclusterblast, knownclusterblast
        Returns:
            a dictionary mapping protein name to Protein instance
    """
    options = get_config()
    if searchtype == "clusterblast":
        logging.info("ClusterBlast: Loading gene cluster database proteins into memory...")
        data_dir = os.path.join(options.database_dir, 'clusterblast')
        protein_file = os.path.join(data_dir, "geneclusterprots.fasta")
    elif searchtype == "subclusterblast":
        logging.info("SubClusterBlast: Loading gene cluster database proteins into memory...")
        data_dir = path.get_full_path(__file__, "data", "sub")
        protein_file = os.path.join(data_dir, "subclusterprots.fasta")
    elif searchtype == "knownclusterblast":
        logging.info("KnownClusterBlast: Loading gene cluster database proteins into memory...")
        data_dir = path.get_full_path(__file__, "data", "known")
        protein_file = os.path.join(data_dir, "knownclusterprots.fasta")

    proteins = {}
    with open(protein_file, 'r') as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line or line[0] != ">":
                continue
            tabs = line.split("|")
            locustag = tabs[4]
            if locustag in accessions:
                locustag = "h_" + locustag # TODO: needs to be actually unique
            location = tabs[2]
            strand = tabs[3]
            annotations = tabs[5]
            name = tabs[6]
            proteins[name] = Protein(name, locustag, location, strand, annotations)
    return proteins

def load_clusterblast_database(record, searchtype="clusterblast") -> Tuple[Dict[str, ReferenceCluster], Dict[str, Protein]]:
    """ Load clusterblast database

        Arguments:
            record: the current record (to pull CDS names from to avoid collisions)
            searchtype: determines which database to use, allowable values:
                            clusterblast, subclusterblast, knownclusterblast
        Returns:
            a tuple of:
                a dictionary mapping cluster name to Cluster instance
                a dictionary mapping protein name to Protein instance
    """
    accessions = set()
    for cds in record.get_cds_features():
        acc = cds.get_accession()
        accessions.add(acc)
    clusters = load_reference_clusters(searchtype)
    proteins = load_reference_proteins(accessions, searchtype)
    return clusters, proteins

def create_blast_inputs(cluster) -> Tuple[List[str], List[str]]:
    """ Creates fasta file contents for the cluster's CDS features

        Arguments:
            cluster: the secmet.Cluster to pull data from

        Returns:
            a tuple of:
                a list of CDS names
                a matching list of CDS sequences
    """
    names = []
    seqs = []
    for cds in cluster.cds_children:
        if cds.strand == 1:
            strand = "+"
        else:
            strand = "-"
        fullname = "|".join(["input", "c%d" % cluster.get_cluster_number(),
                             "%d-%d" % (cds.location.start, cds.location.end),
                             strand, cds.get_accession(), cds.product])
        names.append(fullname)
        seqs.append(cds.get_aa_sequence())

    return names, seqs

def run_internal_blastsearch(query_filename) -> str:
    """ Constructs a blast database from the query and runs blastp on it

        Arguments:
            query_filename: the path of the query fasta file

        Returns:
            a string containing all blastp output
    """
    # TODO... why are query and database the same?
    make_blastdb(query_filename, "internal_input.fasta")
    runblast("internal_input.fasta", "internal_input.fasta")
    with open("internal_input.out", "r") as handle:
        blastoutput = handle.read()
    return blastoutput

def remove_duplicate_hits(blastlines) -> List[List[str]]:
    """ Filter for best blast hits (of one query on each subject)
    """
    query_subject_combinations = set() # type: Set[Tuple[Query, Subject]]
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

def parse_subject(tabs, seqlengths, accessions, record) -> Subject:
    """ Parses a blast-formatted subject line and converts to Subject instance

        Arguments:
            tabs: a list of line parts, the original line split on tabs
            seqlengths: a dictionary of CDS accession to CDS length for calculating
                        percentage of hit coverage
            accessions: a set of CDS locus tags to avoid collisions
            record: the Record to operate on

        Returns:
            a Subject instance
    """
    if len(tabs) < 12:
        logging.error("Malformed blast pairing: %s", "\t".join(tabs))
    query = tabs[0]
    subject_parts = tabs[1].split("|")
    subject = subject_parts[4]
    if subject == "no_locus_tag":
        subject = subject_parts[6]
    if subject in accessions:
        subject = "h_" + subject # TODO should be changed when the other h_ alteration is
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
    blastscore = int(float(tabs[11]) + 0.5)
    cds_accession = query.split("|")[4]
    if cds_accession in seqlengths:
        perc_coverage = (float(tabs[3]) / seqlengths[cds_accession]) * 100
    else:
        feature_by_id = record.get_cds_accession_mapping()
        seqlength = len(feature_by_id[cds_accession].get_aa_sequence())
        perc_coverage = (float(tabs[3]) / seqlength) * 100
    return Subject(subject, genecluster, start, end, strand, annotation,
                   perc_ident, blastscore, perc_coverage, evalue, locustag)

def parse_all_clusters(blasttext, record, min_seq_coverage, min_perc_identity
            ) -> Tuple[Dict[int, Dict[str, List[Query]]], Dict[int, Dict[str, Query]]]:
    """ Parses blast results, groups into results by cluster number

        Arguments:
            blasttext: the output from diamond in blast format
            record: used to get all gene ids in the cluster, and used as a
                    backup to fetch sequence length if missing from seqlengths
            min_seq_coverage: the exclusive lower bound of sequence coverage for a match
            min_perc_identity: the exclusive lower bound of identity similarity for a match

        Returns: a tuple of
                    a dictionary mapping record cluster number to a
                        dictionary of reference clusters names to a
                            list of query hits
                    a dictionary mapping record cluster number to a
                        dictionary of query name to Query instance
    """
    seqlengths = get_cds_lengths(record)
    geneclustergenes = [cds.get_accession() for cds in record.get_cds_features_within_clusters()]
    queries = OrderedDict() # type: Dict[str, Query]
    clusters = OrderedDict() # type: Dict[str, List[Query]]
    blastlines = remove_duplicate_hits([line.split("\t") for line in blasttext.rstrip().splitlines()])
    current_query = None
    queries_by_cluster_number = {} # type: Dict[int, Dict[str, Query]]
    clusters_by_query_cluster_number = {} # type: Dict[int, Dict[str, List[Query]]]

    for tabs in blastlines:
        query = tabs[0]
        subject = parse_subject(tabs, seqlengths, geneclustergenes, record)

        # only process the pairing if limits met
        if subject.perc_ident <= min_perc_identity \
                or subject.perc_coverage <= min_seq_coverage:
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

def blastparse(blasttext, record, min_seq_coverage=-1, min_perc_identity=-1) -> Tuple[Dict[str, Query], Dict[str, List[Query]]]:
    """ Parses blast output into a usable form, limiting to a single best hit
        for every query. Results can be further trimmed by minimum thresholds of
        both coverage and percent identity.

        Arguments:
            blasttext: the output from diamond in blast format
            record: used to get all gene ids in the cluster, and used as a
                    backup to fetch sequence length if missing from seqlengths
            min_seq_coverage: the exclusive lower bound of sequence coverage for a match
            min_perc_identity: the exclusive lower bound of identity similarity for a match

        Returns:
            a tuple of
                a dictionary mapping query id to Query instance
                a dictionary mapping cluster number to
                    a list of Query instances from that cluster
    """
    seqlengths = get_cds_lengths(record)
    accessions = [cds.get_accession() for cds in record.get_cds_features_within_clusters()]
    queries = OrderedDict() # type: Dict[str, Query]
    clusters = OrderedDict() # type: Dict[str, List[Query]]
    blastlines = remove_duplicate_hits([line.split("\t") for line in blasttext.rstrip().split("\n")])
    current_query = None

    for tabs in blastlines:
        query = tabs[0]
        subject = parse_subject(tabs, seqlengths, accessions, record)

        # only process the pairing if limits met
        if subject.perc_ident <= min_perc_identity \
                or subject.perc_coverage <= min_seq_coverage:
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

def get_cds_lengths(record) -> Dict[str, int]:
    """ Calculates the lengths of all CDS features in a Record.

        Arguments:
            record: the Record to gather CDS features from

        Returns:
            a dictionary mapping CDS accession to length of the CDS
    """
    lengths = {}
    for cds in record.get_cds_features():
        lengths[cds.get_accession()] = len(cds.get_aa_sequence())
    return lengths

def find_internal_orthologous_groups(queries, cluster_names) -> List[List[str]]:
    """ Finds internal orthologous groups from blast queries, each cluster
        being a distinct group

        Arguments:
            queries: the queries to build groups from
            cluster_names: the names of the clusters used to construct the queries

        Returns:
            a list of groups, each list being
                a list of query ids
    """ #TODO check description
    groups = []
    for name in cluster_names:
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
                    del groups[i] # TODO: stop changing during iteration
                    for member in other_group:
                        if member not in new_group:
                            new_group.append(member)
                    break
        new_group.sort()
        groups.append(new_group)
    return groups

def internal_homology_blast(record) -> Dict[int, List[List[str]]]:
    """ Run BLAST on gene cluster proteins of each cluster on itself to find
        internal homologs
        store groups of homologs - including singles - in a dictionary
        as a list of lists accordingly

        Arguments:
            record: the Record to generate groups from

        Returns:
            a dictionary mapping cluster_number to
                a list containing distinct groups represented by
                    lists of query ids
    """
    with TemporaryDirectory(change=True):
        logging.info("Finding internal homologs in each gene cluster...")
        internalhomologygroups = {}
        for cluster in record.get_clusters():
            cluster_number = cluster.get_cluster_number()
            iquerycluster_names, iqueryclusterseqs = create_blast_inputs(cluster)
            query_filename = "internal_input.fasta"
            utils.writefasta(iquerycluster_names, iqueryclusterseqs, query_filename)
            blastoutput = run_internal_blastsearch(query_filename)
            queries, _ = blastparse(blastoutput, record, min_seq_coverage=25,
                                    min_perc_identity=30)
            groups = find_internal_orthologous_groups(queries, iquerycluster_names)
            internalhomologygroups[cluster_number] = groups
    return internalhomologygroups

def write_raw_clusterblastoutput(output_dir, blast_output, prefix="clusterblast") -> str:
    """ Writes blast output to file

        NOTE: the output filename will not change for different records, if
              multiple records are in an input and the same output directory
              is used, the file written here will be overwritten

        Arguments:
            output_dir: the path of the directory to store results
            blast_output: the output of blast as a single string
            search_type: the prefix of the filename to create

        Returns:
            the name of the file written
    """
    filename = "%soutput.txt" % prefix
    with open(os.path.join(output_dir, filename), "w") as handle:
        handle.write(blast_output)
    return filename

def parse_clusterblast_dict(queries, clusters, cluster_name, allcoregenes) -> Tuple[Score, List[Tuple[int, int]], List[bool]]:
    """ Generates a score for a cluster, based on the queries and clusters,
        along with the pairings of subjects and queries used to determine that
        score.

        Arguments:
            queries: the queries to determine the score with
            clusters: the clusters to separate queries by
                      (a dict of name -> ReferenceCluster)
            cluster_name: the name of the cluster to score
            allcoregenes: a set of gene ids for bonus scoring

        Returns:
            a tuple of
                a Score instance
                a list of Query, Subject pairs ordered by best pairing
                a list containing an boolean for each pairing, being
                        True only if one of the queries in the pairings was
                        a core gene
    """
    result = Score()
    hitpositions = [] # type: List[Tuple[int, int]]
    hitposcorelist = []
    cluster_locii = clusters[cluster_name][0]
    for query in queries:
        querynrhits = 0
        for subject in query.get_subjects_by_cluster(cluster_name):
            assert cluster_name == subject.genecluster
            if subject.locus_tag not in cluster_locii:
                continue
            index_pair = (query.index, cluster_locii.index(subject.locus_tag))
            if index_pair in hitpositions:
                continue
            querynrhits += 1
            result.blast_score += subject.blastscore
            result.scored_pairings.append((query, subject))
            hitpositions.append(index_pair)
        if querynrhits:
            result.hits += 1
            hit_value = False
            if query.id in allcoregenes:
                result.core_gene_hits += 1
                hit_value = True
            hitposcorelist.extend([hit_value]*querynrhits)
    return result, hitpositions, hitposcorelist

def find_clusterblast_hitsgroups(hitpositions: List[Tuple[int, int]]) -> Dict[int, List[int]]:
    """ Finds groups of hits within a set of pairings

        Arguments:
            hitpositions: a list of Query-Subject index pairs to form groups from

        Returns:
            a dictionary mapping each unique Query to
                a list of Subjects that hit that Query
    """
    groups = defaultdict(list) # type: Dict
    for query, subject in hitpositions:
        groups[query].append(subject)
    return groups

def calculate_synteny_score(hitgroupsdict, hitpositions, core_genes_found) -> int:
    """ Calculate the synteny score of a ranking

        Scores will only be non-zero if more than one hit exists.
        Each query and each hit can only be considered once for scoring.

        Arguments:
            hitgroupsdict: a dictionary mapping each unique Query to
                               a list of Subjects that hit that Query
            hitpositions: a list of Query, Subject index pairs ordered by best pairing
            core_genes_found: a list of booleans for each pair in hitpositions,
                                representing whether a core gene was found in
                                that pairing

        Returns:
            a int representing the synteny score
    """
    scored_queries = set() # type: Set[int]
    scored_hits = set() # type: Set[int]
    synteny_score = 0
    for i, pos in enumerate(hitpositions[:-1]):
        query, subject = pos
        next_query, next_subject = hitpositions[i + 1]
        # Check if a gene homologous to this gene has already been
        # scored for synteny in the previous entry
        if subject in hitgroupsdict[next_query] or query in scored_queries or subject in scored_hits:
            continue
        if (abs(query - next_query) < 2) and abs(query - next_query) == abs(subject - next_subject):
            synteny_score += 1
            if core_genes_found[i] or core_genes_found[i + 1]:
                synteny_score += 1
            scored_queries.add(query)
            scored_hits.add(subject)
    return synteny_score

def score_clusterblast_output(clusters, allcoregenes, cluster_names_to_queries) -> List[Tuple[ReferenceCluster, Score]]:
    """ Generate scores for all cluster matches

        Scores are calculated as S = h + H + s + S + B, where
            h = number of queries with significant hits
            H = number of core genes with significant hits
            s = number of gene pairs with conserved synteny
            S = number of gene pairs with conserved synteny and a core gene involved
            B = core gene bonus (3 if core gene has hit in subject cluster, else 0)

        Rank gene cluster hits based on
            1) number of protein hits
            2) cumulative blast score

        Arguments:
            clusters: a dict mapping reference cluster name to ReferenceCluster
            allcoregenes: a set of all CDS accessions from the current record
            cluster_names_to_queries: a dictionary mapping record cluster number
                                        to a dictionary mapping query name to
                                        Query instance

        Returns:
            A list of ReferenceCluster-Score pairs, sorted in order of
            decreasing score.
    """
    results = {}
    for cluster_name, queries in cluster_names_to_queries.items():
        result, hitpositions, hitposcorelist = parse_clusterblast_dict(queries, clusters, cluster_name, allcoregenes)
        if result.hits <= 1:
            continue
        hitgroupsdict = find_clusterblast_hitsgroups(hitpositions)
        # combines both synteny scores
        result.synteny_score = calculate_synteny_score(hitgroupsdict, hitpositions, hitposcorelist)
        # ensure at least two different subjects were found
        initial = hitpositions[0][1]
        for _, subject in hitpositions[1:]:
            if subject != initial:
                results[clusters[cluster_name]] = result
                break
    #Sort gene clusters by score
    return sorted(results.items(), reverse=True, key=lambda x: x[1].sort_score())

def write_fastas_with_all_genes(clusters, filename, partitions=1) -> List[str]:
    """ Write fasta files containing all genes in all clusters in a
        blast friendly form.

        If partitioning the data into multiple files, the index of the partition
        will be included in the filename before the extension, e.g.
        input.fasta -> input0.fasta, input1.fasta, ...

        Arguments:
            clusters - a list of clusters to find genes in
            filename - the filename to use for the file
            partitions - the number of files to create (approx. equally sized)

        Returns:
            a list containing filenames of the written files
    """
    if not isinstance(partitions, int):
        raise TypeError("Partitions must be an int greater than 0")
    elif partitions < 1:
        raise ValueError("Partitions must be greater than 0")
    all_names, all_seqs = [], []
    for cluster in clusters:
        names, seqs = create_blast_inputs(cluster)
        all_names.extend(names)
        all_seqs.extend(seqs)
    if not (all_names and all_seqs):
        raise ValueError("Diamond search space contains no sequences")
    if partitions == 1:
        utils.writefasta(all_names, all_seqs, filename)
        return [filename]

    chunk_filename = "%d".join(os.path.splitext(filename))
    size = len(all_names) // partitions
    for i in range(partitions):
        chunk_names = all_names[i * size:(i + 1) * size]
        chunk_seqs = all_seqs[i * size:(i + 1) * size]
        if i == partitions - 1:
            chunk_names = all_names[i * size:]
            chunk_seqs = all_seqs[i * size:]
        utils.writefasta(chunk_names, chunk_seqs, chunk_filename % i)
    return [chunk_filename % i for i in range(partitions)]
