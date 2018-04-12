# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" The classification section of the smCOG module. Categorises gene function
    according to a curated set of HMM profiles.
"""

import os
from typing import Dict, Iterable, List

from antismash.common import subprocessing, path, fasta, utils
from antismash.common.hmmscan_refinement import refine_hmmscan_results, HMMResult
from antismash.common.secmet import GeneFunction, CDSFeature
from antismash.config import ConfigType


def classify_genes(cds_features: List[CDSFeature]) -> Dict[str, List[HMMResult]]:
    """ Finds possible classifications for the provided genes.

        Arguments:
            cds_features: a list of CDSFeatures to classify

        Returns:
            a dictionary mapping CDS name to a list of HMMResult instances of
                classifications
    """
    smcogs_fasta = fasta.get_fasta_from_features(cds_features)
    smcogs_opts = ["-E", "1E-6"]
    hmm_file = path.get_full_path(__file__, "data", "smcogs.hmm")
    smcogs_results = subprocessing.run_hmmscan(hmm_file, smcogs_fasta, smcogs_opts)
    hmm_lengths = utils.get_hmm_lengths(hmm_file)
    return refine_hmmscan_results(smcogs_results, hmm_lengths)


def load_cog_annotations() -> Dict[str, GeneFunction]:
    """ Load the smCOG type annotations from a file, returns a dictionary mapping
        smCOG id to the gene function of that smCOG.
    """
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


def write_smcogs_file(hmm_results: Dict[str, List[HMMResult]], cds_features: Iterable[CDSFeature],
                      nrpspks_genes: Iterable[CDSFeature], options: ConfigType) -> None:
    """ Writes a text file containing the smCOG results to the output directory
        defined in options
    """
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
