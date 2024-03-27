# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

from typing import Optional

from antismash.common.secmet import Record
from antismash.detection.genefunctions.halogenases.halogenases import TailoringEnzymes
from antismash.detection.genefunctions.halogenases.flavin_dependent.substrate_analysis import fdh_specific_analysis

def specific_analysis(record: Record)\
    -> Optional[list[TailoringEnzymes]]:
    return fdh_specific_analysis(record)
