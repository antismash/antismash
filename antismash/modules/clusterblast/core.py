# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Shared functions across variants of clusterblast """

from collections import defaultdict, OrderedDict
import logging
import os
from tempfile import NamedTemporaryFile
from typing import Dict, Iterable, List, Set, Sequence, Tuple

from helperlibs.wrappers.io import TemporaryDirectory

from antismash.common import path, subprocessing, fasta, secmet
from antismash.common.subprocessing.diamond import \
    check_diamond_files as check_clusterblast_files  # used by other components # pylint: disable=unused-import

from antismash.common.subprocessing.blast import run_makeblastdb as make_blastdb
from antismash.config import get_config

from .data_structures import Subject, Query, Protein, ReferenceCluster, Score

_SHIPPED_DATA_DIR = path.get_full_path(__file__, "data")


def get_core_gene_ids(record: secmet.Record) -> Set[str]:  # TODO: consider moving into secmet
    """ Fetches all gene accessions of genes with CORE gene function from all
        clusters in a record

        Arguments:
            record: the record to gather gene ids from

        Returns:
            a set containing all core gene names
    """
    cores = set()
    for gene in record.get_cds_features_within_regions():
        if gene.gene_function == secmet.GeneFunction.CORE:
            cores.add(gene.get_name())
    return cores


def run_blast(query: str, database: str) -> str:
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
        raise RuntimeError(f"blastp run failed: ...{res.stderr[-200:]!r}")
    return out_file


def run_diamond_on_all_regions(regions: Sequence[secmet.Region], database: str) -> str:
    """ Runs diamond, comparing all features in the given regions to the given database

        Arguments:
            regions: the regions to use features from
            database: the path of the database to compare to

        Returns:
            diamond's output from stdout
    """
    logging.info("Comparing regions to reference database")
    extra_args = [
        "--compress", "0",
        "--max-target-seqs", "10000",
        "--evalue", "1e-05",
        "--outfmt", "6",  # 6 is blast tabular format, just as in blastp
    ]
    with NamedTemporaryFile() as temp_file:
        write_fastas_with_all_genes(regions, temp_file.name)
        stdout = subprocessing.run_diamond_search(temp_file.name, database, mode="blastp", opts=extra_args)
    return stdout


def load_reference_clusters(searchtype: str) -> Dict[str, ReferenceCluster]:
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
        logging.info("ClusterBlast: Loading gene cluster database into memory...")
        data_dir = os.path.join(options.database_dir, 'clusterblast')
    elif searchtype == "subclusterblast":
        logging.info("SubClusterBlast: Loading gene cluster database into memory...")
        data_dir = os.path.join(_SHIPPED_DATA_DIR, "sub")
    elif searchtype == "knownclusterblast":
        logging.info("KnownClusterBlast: Loading gene cluster database into memory...")
        kcb_root = os.path.join(options.database_dir, "knownclusterblast")
        version = path.find_latest_database_version(kcb_root)
        data_dir = os.path.join(kcb_root, version)

    reference_cluster_file = os.path.join(data_dir, "clusters.txt")
    with open(reference_cluster_file, "r", encoding="utf-8") as handle:
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


def load_reference_proteins(searchtype: str) -> Dict[str, Protein]:
    """ Load protein database

        Arguments:
            searchtype: determines which database to use, allowable values:
                            clusterblast, subclusterblast, knownclusterblast
        Returns:
            a dictionary mapping protein name to Protein instance
    """
    options = get_config()
    if searchtype == "clusterblast":
        logging.info("ClusterBlast: Loading gene cluster database proteins into memory...")
        data_dir = os.path.join(options.database_dir, 'clusterblast')
    elif searchtype == "subclusterblast":
        logging.info("SubClusterBlast: Loading gene cluster database proteins into memory...")
        data_dir = os.path.join(_SHIPPED_DATA_DIR, "sub")
    elif searchtype == "knownclusterblast":
        logging.info("KnownClusterBlast: Loading gene cluster database proteins into memory...")
        kcb_root = os.path.join(options.database_dir, "knownclusterblast")
        version = path.find_latest_database_version(kcb_root)
        data_dir = os.path.join(kcb_root, version)

    protein_file = os.path.join(data_dir, "proteins.fasta")
    proteins = {}
    with open(protein_file, "r", encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line or line[0] != ">":
                continue
            # some lines are malformed, so always split the name off the annotation
            # e.g. >x|y|1-2|-|z|Urea_carboxylase_{ECO:0000313|EMBL:CCF11062.1}|CRH36422
            tabs = line.split("|", 5)
            annotations, name = tabs[5].rsplit("|", 1)
            locustag = tabs[4]
            location = tabs[2]
            strand = tabs[3]
            proteins[locustag] = Protein(name, locustag, location, strand, annotations)
    return proteins


def load_clusterblast_database(searchtype: str = "clusterblast"
                               ) -> Tuple[Dict[str, ReferenceCluster], Dict[str, Protein]]:
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
    clusters = load_reference_clusters(searchtype)
    proteins = load_reference_proteins(searchtype)

    return clusters, proteins


def create_blast_inputs(region: secmet.Region) -> Tuple[List[str], List[str]]:
    """ Creates fasta file contents for the cluster's CDS features

        Arguments:
            region: the secmet.Region to pull data from

        Returns:
            a tuple of:
                a list of CDS names
                a matching list of CDS sequences
    """
    names = []
    seqs = []
    for cds in region.cds_children:
        if cds.strand == 1:
            strand = "+"
        else:
            strand = "-"
        fullname = "|".join(["input", f"c{region.get_region_number()}",
                             f"{cds.location.start:d}-{cds.location.end:d}",
                             strand, cds.get_name(), cds.product])
        names.append(fullname)
        seqs.append(cds.translation)

    return names, seqs


def run_internal_blastsearch(query_filename: str) -> str:
    """ Constructs a blast database from the query and runs blastp on it to
        find internal matches

        Arguments:
            query_filename: the path of the query fasta file

        Returns:
            a string containing all blastp output
    """
    make_blastdb(query_filename, "internal_input.fasta")
    run_blast("internal_input.fasta", "internal_input.fasta")
    with open("internal_input.out", "r", encoding="utf-8") as handle:
        blastoutput = handle.read()
    return blastoutput


def remove_duplicate_hits(blast_lines: List[List[str]]) -> List[List[str]]:
    """ Filter for best blast hits (of one query on each subject)

        Arguments:
            blast_lines: a list of lists, each inner list being a single line of
                         the blast input after splitting it up

        Returns:
            a subset of the input, keeping only the first hit for each pairing
    """
    query_subject_combinations: Set[Tuple[str, str]] = set()
    deduplicated = []
    for tabs in blast_lines:
        # if it doesn't even have both values, it's not a hit so skip it
        if len(tabs) < 2:
            continue
        query = tabs[0]
        subject = tabs[1]
        query_subject_combination = (query, subject)
        if query_subject_combination not in query_subject_combinations:
            query_subject_combinations.add(query_subject_combination)
            deduplicated.append(tabs)
    return deduplicated


def parse_subject(line_parts: List[str], seqlengths: Dict[str, int],
                  record: secmet.Record) -> Subject:
    """ Parses a blast-formatted subject line and converts to Subject instance

        Arguments:
            line_parts: a list of line parts, the original line split on line_parts
            seqlengths: a dictionary of CDS accession to CDS length for calculating
                        percentage of hit coverage
            names: a set of CDS names to avoid collisions
            record: the Record to operate on

        Returns:
            a Subject instance
    """
    if len(line_parts) < 12:
        logging.error("Malformed blast pairing: %s", "\t".join(line_parts))
    query = line_parts[0]
    subject_parts = line_parts[1].split("|")
    subject = subject_parts[4]
    if subject == "no_locus_tag":
        subject = subject_parts[6]
    if len(subject_parts) > 6:
        locustag = subject_parts[6]
    else:
        locustag = ""
    genecluster = f"{subject_parts[0]}_{subject_parts[1]}"
    start, end = subject_parts[2].split("-")[:2]
    strand = subject_parts[3]
    annotation = subject_parts[5]
    perc_ident = int(float(line_parts[2]) + 0.5)
    evalue = float(line_parts[10])
    blastscore = int(float(line_parts[11]) + 0.5)
    cds_name = query.split("|")[4]
    if cds_name in seqlengths:
        perc_coverage = (float(line_parts[3]) / seqlengths[cds_name]) * 100
    else:
        seqlength = len(record.get_cds_by_name(cds_name).translation)
        perc_coverage = (float(line_parts[3]) / seqlength) * 100
    return Subject(subject, genecluster, int(start), int(end), strand, annotation,
                   perc_ident, blastscore, perc_coverage, evalue, locustag)


def parse_all_clusters(blasttext: str, record: secmet.Record, min_seq_coverage: float, min_perc_identity: float
                       ) -> Tuple[Dict[int, Dict[str, List[Query]]],
                                  Dict[int, Dict[str, Query]]]:
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
    queries: Dict[str, Query] = OrderedDict()
    clusters: Dict[str, List[Query]] = OrderedDict()
    blastlines = remove_duplicate_hits([line.split("\t") for line in blasttext.rstrip().splitlines()])
    current_query = None
    queries_by_cluster_number: Dict[int, Dict[str, Query]] = {}
    clusters_by_query_cluster_number: Dict[int, Dict[str, List[Query]]] = {}

    for tabs in blastlines:
        query = tabs[0]
        subject = parse_subject(tabs, seqlengths, record)

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

        assert current_query is not None

        if subject.genecluster not in clusters:
            clusters[subject.genecluster] = []
        clusters[subject.genecluster].append(current_query)

        # link the subject to the query
        current_query.add_subject(subject)

    return clusters_by_query_cluster_number, queries_by_cluster_number


def blastparse(blasttext: str, record: secmet.Record, min_seq_coverage: float = -1.,
               min_perc_identity: float = -1.) -> Tuple[Dict[str, Query], Dict[str, List[Query]]]:
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
    queries: Dict[str, Query] = OrderedDict()
    clusters: Dict[str, List[Query]] = OrderedDict()
    blastlines = remove_duplicate_hits([line.split("\t") for line in blasttext.rstrip().split("\n")])
    current_query = None

    for tabs in blastlines:
        query = tabs[0]
        subject = parse_subject(tabs, seqlengths, record)

        # only process the pairing if limits met
        if subject.perc_ident <= min_perc_identity \
                or subject.perc_coverage <= min_seq_coverage:
            continue

        new_query = query not in queries
        new_hit = subject.genecluster not in clusters

        if new_query:
            current_query = Query(query, len(queries))
            queries[query] = current_query
        assert current_query is not None

        if new_hit:
            clusters[subject.genecluster] = []
        clusters[subject.genecluster].append(current_query)

        # link the subject to the query
        current_query.add_subject(subject)

    return queries, clusters


def get_cds_lengths(record: secmet.Record) -> Dict[str, int]:
    """ Calculates the lengths of each CDS feature in a Record.

        Arguments:
            record: the Record to gather CDS features from

        Returns:
            a dictionary mapping CDS accession to length of the CDS
    """
    lengths = {}
    for cds in record.get_cds_features():
        lengths[cds.get_name()] = len(cds.translation)
    return lengths


def find_internal_orthologous_groups(queries: Dict[str, Query], cluster_names: List[str]) -> List[List[str]]:
    """ Finds internal orthologous groups from blast queries, each cluster
        being a distinct group

        Arguments:
            queries: the queries to build groups from
            cluster_names: the names of the clusters used to construct the queries,
                           e.g. rec_id|c1|52644-53715|+|protein_id|description

        Returns:
            a list of groups, each list being
                a list of query ids in alphabetical order
    """

    groups = []
    for name in cluster_names:
        if name not in queries:
            groups.append(set([name.split("|")[4]]))
            continue
        query = queries[name]
        new_group = {query.id}.union(set(query.subjects))
        redundant_groups = []
        for i, other_group in enumerate(groups):
            if new_group.intersection(other_group):
                redundant_groups.append(i)
                new_group.update(other_group)

        for i in reversed(redundant_groups):
            del groups[i]

        groups.append(new_group)
    return [sorted(list(i)) for i in groups]


def internal_homology_blast(record: secmet.Record) -> Dict[int, List[List[str]]]:
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
        internalhomologygroups: Dict[int, List[List[str]]] = {}
        for region in record.get_regions():
            region_number = region.get_region_number()
            if not region.cds_children:
                internalhomologygroups[region_number] = []
                continue
            iquerycluster_names, iqueryclusterseqs = create_blast_inputs(region)
            query_filename = "internal_input.fasta"
            fasta.write_fasta(iquerycluster_names, iqueryclusterseqs, query_filename)
            blastoutput = run_internal_blastsearch(query_filename)
            queries, _ = blastparse(blastoutput, record, min_seq_coverage=25,
                                    min_perc_identity=30)
            groups = find_internal_orthologous_groups(queries, iquerycluster_names)
            internalhomologygroups[region_number] = groups
    return internalhomologygroups


def parse_clusterblast_dict(queries: List[Query], clusters: Dict[str, ReferenceCluster],
                            cluster_label: str, allcoregenes: Set[str]
                            ) -> Tuple[Score, List[Tuple[int, int]], List[bool]]:
    """ Generates a score for a cluster, based on the queries and clusters,
        along with the pairings of subjects and queries used to determine that
        score.

        Arguments:
            queries: the queries to determine the score with
            clusters: the clusters to separate queries by
                      (a dict of cluster label -> ReferenceCluster)
            cluster_label: the number of the cluster to score
            allcoregenes: a set of gene ids for bonus scoring

        Returns:
            a tuple of
                a Score instance
                a list of Query, Subject pairs ordered by best pairing
                a list containing a boolean for each pairing, being
                        True only if one of the queries in the pairings was
                        a core gene
    """
    result = Score()
    hitpositions: List[Tuple[int, int]] = []
    hitposcorelist = []
    cluster_locii = clusters[cluster_label].tags
    for query in queries:
        querynrhits = 0
        for subject in query.get_subjects_by_cluster(cluster_label):
            assert cluster_label == subject.genecluster
            if subject.name not in cluster_locii:
                continue
            index_pair = (query.index, cluster_locii.index(subject.name))
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
    groups: Dict[int, List[int]] = defaultdict(list)
    for query, subject in hitpositions:
        groups[query].append(subject)
    return groups


def calculate_synteny_score(hitgroups: Dict[int, List[int]],
                            hitpositions: List[Tuple[int, int]],
                            core_genes_found: List[bool]) -> int:
    """ Calculate the synteny score of a ranking

        Scores will only be non-zero if more than one hit exists.
        Each query and each hit can only be considered once for scoring.

        Arguments:
            hitgroups: a dictionary mapping each unique Query index to
                               a list of Subject indices that hit that Query
            hitpositions: a list of Query, Subject index pairs ordered by best pairing
            core_genes_found: a list of booleans for each pair in hitpositions,
                                representing whether a core gene was found in
                                that pairing

        Returns:
            a int representing the synteny score
    """
    scored_queries: Set[int] = set()
    scored_hits: Set[int] = set()
    synteny_score = 0
    for i, pos in enumerate(hitpositions[:-1]):
        query, subject = pos
        next_query, next_subject = hitpositions[i + 1]
        # Check if a gene homologous to this gene has already been
        # scored for synteny in the previous entry
        if subject in hitgroups[next_query] or query in scored_queries or subject in scored_hits:
            continue
        if (abs(query - next_query) < 2) and abs(query - next_query) == abs(subject - next_subject):
            synteny_score += 1
            if core_genes_found[i] or core_genes_found[i + 1]:
                synteny_score += 1
            scored_queries.add(query)
            scored_hits.add(subject)
    return synteny_score


def score_clusterblast_output(clusters: Dict[str, ReferenceCluster], allcoregenes: Set[str],
                              cluster_names_to_queries: Dict[str, List[Query]]
                              ) -> List[Tuple[ReferenceCluster, Score]]:
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
            cluster_names_to_queries: a dictionary mapping reference cluster name
                                        to a list of matching Query instances

        Returns:
            A list of ReferenceCluster-Score pairs, sorted in order of
            decreasing score.
    """
    results = {}
    for cluster_label, queries in cluster_names_to_queries.items():
        single_gene_reference = len(clusters[cluster_label].proteins) == 1
        result, hitpositions, hitposcorelist = parse_clusterblast_dict(queries, clusters, cluster_label, allcoregenes)
        if not result.hits:
            continue
        if not single_gene_reference and result.hits <= 1:
            continue
        hitgroups = find_clusterblast_hitsgroups(hitpositions)
        # combines both synteny scores
        result.synteny_score = calculate_synteny_score(hitgroups, hitpositions, hitposcorelist)
        initial = hitpositions[0][1]
        # if only one gene in reference, use it
        if single_gene_reference and len(hitpositions) == 1:
            results[clusters[cluster_label]] = result
        # otherwise ensure at least two different subjects were found
        for _, subject in hitpositions[1:]:
            if subject != initial:
                results[clusters[cluster_label]] = result
                break
    # Sort gene clusters by score
    return sorted(results.items(), reverse=True, key=lambda x: x[1].sort_score())


def write_fastas_with_all_genes(regions: Iterable[secmet.Region], filename: str,
                                partitions: int = 1) -> List[str]:
    """ Write fasta files containing all genes in all clusters in a
        blast friendly form.

        If partitioning the data into multiple files, the index of the partition
        will be included in the filename before the extension, e.g.
        input.fasta -> input0.fasta, input1.fasta, ...

        Arguments:
            regions: an iterable of clusters to find genes in
            filename: the filename to use for the file
            partitions: the number of files to create (approx. equally sized)

        Returns:
            a list containing filenames of the written files
    """
    if not isinstance(partitions, int):
        raise TypeError("Partitions must be an int greater than 0")
    if partitions < 1:
        raise ValueError("Partitions must be greater than 0")

    all_names, all_seqs = [], []
    for region in regions:
        names, seqs = create_blast_inputs(region)
        all_names.extend(names)
        all_seqs.extend(seqs)
    if not (all_names and all_seqs):
        raise ValueError("Diamond search space contains no sequences")
    if partitions == 1:
        fasta.write_fasta(all_names, all_seqs, filename)
        return [filename]

    chunk_filename = "%d".join(os.path.splitext(filename))
    size = len(all_names) // partitions
    for i in range(partitions):
        chunk_names = all_names[i * size:(i + 1) * size]
        chunk_seqs = all_seqs[i * size:(i + 1) * size]
        if i == partitions - 1:
            chunk_names = all_names[i * size:]
            chunk_seqs = all_seqs[i * size:]
        fasta.write_fasta(chunk_names, chunk_seqs, chunk_filename % i)
    return [chunk_filename % i for i in range(partitions)]
