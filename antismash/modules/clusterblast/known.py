# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

import logging
import os
from helperlibs.wrappers.io import TemporaryDirectory

import antismash.common.deprecated as utils
import antismash.common.path as path

from .core import parse_all_clusters, \
                  load_clusterblast_database, create_blast_inputs, run_diamond, \
                  write_raw_clusterblastoutput, score_clusterblast_output
from .results import ClusterResult, GeneralResults, write_clusterblast_output

# Tuple is ( binary_name, optional)
_required_binaries = [
    ('blastp', False),
    ('makeblastdb', False),
    ('diamond', False),
]

_required_files = [
    ('knownclusterprots.fasta', False),
    ('knownclusterprots.dmnd', False),
    ('knownclusters.txt', False)
]

def _get_datafile_path(filename):
    data_dir = path.get_full_path(__file__, 'data')
    return os.path.join(data_dir, 'known', filename)

def check_known_prereqs(options):
    "Check if all required applications are around"
    failure_messages = []
    for binary_name, optional in _required_binaries:
        if path.locate_executable(binary_name) is None and not optional:
            failure_messages.append("Failed to locate file: %r" % binary_name)

    for file_name, optional in _required_files:
        if path.locate_file(_get_datafile_path(file_name)) is None and not optional:
            failure_messages.append("Failed to locate file: %r" % file_name)

    return failure_messages

def run_knownclusterblast_on_record(seq_record, options):
    logging.info('Running known cluster search')
    clusters, proteins = load_clusterblast_database(seq_record, searchtype="knownclusterblast")
    return perform_knownclusterblast(options, seq_record, clusters, proteins)

def perform_knownclusterblast(options, seq_record, clusters, proteins):
    # Run BLAST on gene cluster proteins of each cluster and parse output
    logging.info("Running DIAMOND knowncluster searches..")
    results = GeneralResults(seq_record.id, search_type="knownclusterblast")

    all_names, all_seqs = [], []
    for cluster in seq_record.get_clusters():
        names, seqs = create_blast_inputs(cluster)
        all_names.extend(names)
        all_seqs.extend(seqs)
    if not (all_names and all_seqs):
        raise RuntimeError("Diamond search space contains no sequences")
    with TemporaryDirectory(change=True) as tempdir:
        utils.writefasta([qcname.replace(" ", "_") for qcname in all_names],
                         all_seqs, "input.fasta")
        run_diamond("input.fasta", _get_datafile_path('knownclusterprots'),
                    tempdir, options)
        with open("input.out", 'r') as handle:
            blastoutput = handle.read()
        write_raw_clusterblastoutput(options.output_dir, blastoutput,
                                     search_type="knownclusterblast")
    minseqcoverage = 40
    minpercidentity = 45
    clusters_by_number, _ = parse_all_clusters(blastoutput, minseqcoverage,
                                              minpercidentity, seq_record)

    allcoregenes = seq_record.get_cds_features()
    for genecluster in seq_record.get_clusters():
        clusternumber = genecluster.get_cluster_number()
        cluster_names_to_queries = clusters_by_number.get(clusternumber, {})
        ranking = score_clusterblast_output(clusters, allcoregenes, cluster_names_to_queries)
        # store results
        cluster_result = ClusterResult(genecluster, ranking, proteins)
        results.add_cluster_result(cluster_result, clusters, proteins)

        write_clusterblast_output(options, seq_record, cluster_result, proteins,
                                  searchtype="knownclusterblast")
    results.mibig_entries = mibig_protein_homology(blastoutput, seq_record, clusters, options)
    return results

class MibigEntry:
    def __init__(self, gene_id, gene_description, mibig_cluster,
                mibig_product, percent_id, blast_score, coverage, evalue):
        self.gene_id = gene_id
        self.gene_description = gene_description
        self.mibig_id = mibig_cluster.split("_c")[0]
        self.mibig_product = mibig_product
        self.percent_id = float(percent_id)
        self.blast_score = float(blast_score)
        self.coverage = float(coverage)
        self.evalue = float(evalue)

    @property
    def values(self):
        return [self.gene_id, self.gene_description, self.mibig_id,
                self.mibig_product, self.percent_id, self.blast_score,
                self.coverage, self.evalue]

    def __str__(self):
        return "%s\n" % "\t".join(str(val) for val in self.values)

def mibig_protein_homology(blastoutput, seq_record, clusters, options):
    """ Constructs a mapping of gene to MiBiG hits
        Returns a dict of dicts of lists, accessed by:
            mibig_entries[cluster_number][gene_accession] = list of MibigEntry
    """
    minseqcoverage = 20
    minpercidentity = 20
    _, queries_by_cluster = parse_all_clusters(blastoutput, minseqcoverage,
                                               minpercidentity, seq_record)
    mibig_entries = {}

    for cluster in seq_record.get_clusters():
        cluster_number = cluster.get_cluster_number()
        queries = queries_by_cluster.get(cluster_number, {})
        cluster_entries = {}
        # Since the BLAST query was only for proteins in the cluster just need to iterate through the keys
        for cluster_protein in queries.values():
            protein_name = cluster_protein.id
            protein_entries = []
            for subject in cluster_protein.subjects.values():
                entry = MibigEntry(subject.locus_tag, subject.annotation,
                                   subject.genecluster,
                                   clusters[subject.genecluster].cluster_type,
                                   subject.perc_ident, subject.blastscore,
                                   subject.perc_coverage, subject.evalue)
                protein_entries.append(entry)
            cluster_entries[protein_name] = protein_entries
        if cluster_entries:
            mibig_entries[cluster_number] = cluster_entries
    return mibig_entries
