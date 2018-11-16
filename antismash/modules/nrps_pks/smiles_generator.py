# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Generation of structure predictions for product compounds in both SMILES and
    image format, where possible
"""

import logging
from typing import Dict

import antismash.common.path as path


def gen_smiles_from_pksnrps(compound_pred: str) -> str:
    """ Generates the SMILES string for a specific compound prediction """
    smiles = ""

    if not compound_pred:
        return smiles

    residues = compound_pred.replace("(", "").replace(")", "").replace(" + ", " ").replace(" - ", " ").split()
    # Counts the number of malonate and its derivatives in polyketides
    mal_count = 0
    for i in residues:
        if "mal" in i:
            mal_count += 1

    # Reflecting reduction states of ketide groups starting at beta carbon of type 1 polyketide
    if "pk" in residues and "mal" in residues[-1]:
        residues.pop(residues.index('pk')+1)
        residues.append('pks-end1')
    elif mal_count == len(residues):
        if residues[0] == "mal":
            residues[0] = "pks-start1"  # TODO why replace and not insert?
        if residues[-1] == "ccmal":
            residues.append('pks-end2')

    if len(residues) > 1:
        # Conventionally aaSMILES was used;
        # chirality expressed with "@@" causes indigo error
        aa_smiles = load_smiles()

        for monomer in residues:
            if monomer in aa_smiles:
                smiles += aa_smiles[monomer]
            elif '|' in monomer:
                logging.debug("Substituting 'nrp' for combined monomer %r", monomer)
                smiles += aa_smiles['nrp']
            else:
                logging.debug("No SMILES mapping for unknown monomer %r", monomer)
    return smiles


def load_smiles() -> Dict[str, str]:
    """Load smiles from a dictionary mapping residues to SMILES string"""
    aa_smiles = {}  # type: Dict[str, str]

    smiles_monomer = open(path.get_full_path(__file__, 'data', 'aaSMILES.txt'), 'r')

    for line in smiles_monomer.readlines():
        line = line.strip()
        if not line or line.startswith('#') or line == "END":
            continue
        smiles = line.split()
        assert len(smiles) == 2, "Invalid smiles line {!r}".format(line)
        assert smiles[0] not in aa_smiles, "%s contained twice in smiles data" % smiles[0]
        aa_smiles[smiles[0]] = smiles[1]

    smiles_monomer.close()
    return aa_smiles
