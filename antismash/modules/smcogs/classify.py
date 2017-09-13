# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

import os

from antismash.common import subprocessing, path, deprecated
from antismash.common.hmmscan_parser import refine_hmmscan_results
from antismash.common.secmet import GeneFunction

def classify_genes(cds_features):
    smcogs_fasta = deprecated.get_specific_multifasta(cds_features)
    smcogs_opts = ["-E", "1E-6"]
    hmm_file = path.get_full_path(__file__, "data/smcogs.hmm")
    smcogs_results = subprocessing.run_hmmscan(hmm_file, smcogs_fasta, smcogs_opts)
    hmm_lengths = deprecated.hmmlengths(hmm_file)
    return refine_hmmscan_results(smcogs_results, hmm_lengths)


def load_cog_annotations():
    "Load the smCOG type annotations from a file"
    mapping = {
        'B': GeneFunction.ADDITIONAL, #'biosynthetic-additional',
        'T': GeneFunction.TRANSPORT, #'transport',
        'R': GeneFunction.REGULATORY, #'regulatory',
    }
    annotations = {}
    for line in open(path.get_full_path(__file__, 'data/cog_annotations.txt'), 'r'):
        cog, _, key = line.strip().split('\t', 3)
        annotations[cog] = mapping.get(key, GeneFunction.OTHER)

    return annotations

def write_smcogs_file(hmm_results, cds_features, nrpspks_genes, options):
    nrpspks_names = set([feature.get_name() for feature in nrpspks_genes])
    smcogfile = open(os.path.join(options.output_dir, "smcogs", "smcogs.txt"), "w")
    for feature in cds_features:
        k = feature.get_name()
        if k not in nrpspks_names:
            if k in hmm_results:
                l = hmm_results[k]
                smcogfile.write(">> " + k + "\n")
                smcogfile.write("name\tstart\tend\te-value\tscore\n") # TODO: convert to results
                smcogfile.write("** smCOG hits **\n")
                for i in l:
                    smcogfile.write(str(i[0]) + "\t" + str(i[1]) + "\t" + str(i[2]) + "\t" + str(i[3]) + "\t" + str(i[4]) + "\n")
                smcogfile.write("\n\n")
    smcogfile.close()
