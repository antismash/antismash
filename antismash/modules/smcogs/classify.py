# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

import os

from antismash.common import subprocessing, path, deprecated
from antismash.common.hmmscan_refinement import refine_hmmscan_results
from antismash.common.secmet import GeneFunction


def classify_genes(cds_features):
    smcogs_fasta = deprecated.get_specific_multifasta(cds_features)
    smcogs_opts = ["-E", "1E-6"]
    hmm_file = path.get_full_path(__file__, "data", "smcogs.hmm")
    smcogs_results = subprocessing.run_hmmscan(hmm_file, smcogs_fasta, smcogs_opts)
    hmm_lengths = deprecated.hmmlengths(hmm_file)
    return refine_hmmscan_results(smcogs_results, hmm_lengths)


def load_cog_annotations():
    "Load the smCOG type annotations from a file"
    mapping = {
        'B': GeneFunction.ADDITIONAL,  # 'biosynthetic-additional',
        'T': GeneFunction.TRANSPORT,  # 'transport',
        'R': GeneFunction.REGULATORY,  # 'regulatory',
    }
    annotations = {}
    for line in open(path.get_full_path(__file__, 'data', 'cog_annotations.txt'), 'r'):
        cog, _, key = line.strip().split('\t', 3)
        annotations[cog] = mapping.get(key, GeneFunction.OTHER)

    return annotations


def write_smcogs_file(hmm_results, cds_features, nrpspks_genes, options):
    nrpspks_names = set(feature.get_name() for feature in nrpspks_genes)
    # TODO don't overwrite with multiple records
    smcogfile = open(os.path.join(options.output_dir, "smcogs", "smcogs.txt"), "w")
    for feature in sorted(cds_features, key=lambda feat: feat.location.start):
        gene_id = feature.get_name()
        if gene_id in nrpspks_names:
            continue
        if gene_id in hmm_results:
            hits = hmm_results[gene_id]
            smcogfile.write(">> %s\n" % gene_id)
            smcogfile.write("name\tstart\tend\te-value\tscore\n")
            smcogfile.write("** smCOG hits **\n")
            for hit in hits:
                smcogfile.write("\t".join([hit.hit_id, str(hit.query_start),
                                           str(hit.query_end), str(hit.evalue),
                                           str(hit.bitscore)]) + "\n")
            smcogfile.write("\n\n")
    smcogfile.close()
