# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Individual tools/methods for detection of gene functions """

import os
from typing import List

from antismash.common.secmet import CDSFeature, Record
from antismash.common.secmet.qualifiers import GeneFunction
from antismash.config import ConfigType

from .core import FunctionResults, scan_for_functions
from .smcogs import classify as smcogs_classification


def run_tools(record: Record, options: ConfigType) -> List[FunctionResults]:
    """ Runs all tools over CDS features within regions

        Arguments:
            record: the record to check
            options: antiSMASH config

        Returns:
            a list of FunctionResults objects, one for each tool
    """

    cds_features = record.get_cds_features_within_regions()
    functions = [
        find_resistance,
        smcogs_classification,
    ]
    return [func(record.id, cds_features, options) for func in functions]


def find_resistance(record_id: str, cds_features: List[CDSFeature], options: ConfigType) -> FunctionResults:
    """ Finds resistance genes with the Resfam core database
        as supplied by: http://www.dantaslab.org/resfams.

        Further information about the generation of the Resfam database can be found
        in:

        Gibson MK, Forsberg KJ, Dantas G.
        Improved annotation of antibiotic resistance functions reveals microbial
            resistomes cluster by ecology.
        The ISME Journal. 2014, doi:ISMEJ.2014.106

        Arguments:
            record_id: the name of the record
            cds_features: the CDS features to check for functions
            options: the antiSMASH config

        Returns:
            a FunctionResults instance containing the best hits for each gene, if any
    """
    hmm_file = os.path.join(options.database_dir, "resfam", "Resfams.hmm")
    hits = scan_for_functions(cds_features, hmm_file, hmmscan_opts=["--cut_ga"])
    mapping = {cds_name: GeneFunction.RESISTANCE for cds_name in hits}
    return FunctionResults(record_id, "resist", hits, mapping)
