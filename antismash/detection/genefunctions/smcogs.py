# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" The classification section of the smCOG module. Categorises gene function
    according to a curated set of HMM profiles.
"""

from typing import Dict, List

from antismash.common import path
from antismash.common.hmmscan_refinement import HMMResult
from antismash.common.secmet import GeneFunction, CDSFeature
from antismash.config import ConfigType

from .core import FunctionResults, scan_for_functions


def classify(record_id: str, cds_features: List[CDSFeature],  # an API, so hide unused warning
             options: ConfigType) -> FunctionResults:  # pylint: disable=unused-argument
    """ Finds possible classifications for the provided CDS features.

        Arguments:
            cds_features: a list of CDSFeatures to classify

        Returns:
            a dictionary mapping CDS name to a list of HMMResult instances of
                classifications
    """

    hmm_file = path.get_full_path(__file__, "data", "smcogs.hmm")
    hits = scan_for_functions(cds_features, hmm_file, hmmscan_opts=["-E", "1E-16"])
    ids_to_function = build_function_mapping()
    cds_name_to_function = {}
    for cds_name, result in hits.items():
        smcog_id = result.hit_id.split(":", 1)[0]
        cds_name_to_function[cds_name] = ids_to_function[smcog_id]
        hits[cds_name] = HMMResult(result.hit_id.replace("_", " "), result.query_start,
                                   result.query_end, result.evalue, result.bitscore)
    return FunctionResults(record_id, "smcogs", hits, cds_name_to_function)


def build_function_mapping() -> Dict[str, GeneFunction]:
    """ Load the smCOG gene function mapping from a file.

        Returns:
             a dictionary mapping smCOG id to its matching gene function
    """
    mapping = {
        'B': GeneFunction.ADDITIONAL,  # 'biosynthetic-additional',
        'T': GeneFunction.TRANSPORT,  # 'transport',
        'R': GeneFunction.REGULATORY,  # 'regulatory',
    }
    annotations = {}
    with open(path.get_full_path(__file__, "data", "cog_annotations.txt"), "r", encoding="utf-8") as handle:
        lines = handle.readlines()
    for line in lines:
        cog, _desc, key = line.strip().split('\t', 3)
        annotations[cog] = mapping.get(key, GeneFunction.OTHER)

    return annotations
