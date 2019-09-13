# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Generation of structure predictions for product compounds in both SMILES and
    image format, where possible
"""

import logging
from typing import Dict, Iterator, List, Tuple

import antismash.common.path as path

_STARTING_BONDS = {
    "C": 4,
    "c": 4,
    "O": 2,
    "N": 3,
}


class Atom:
    """ A construction that represents an atom within a SMILES context.

        Arguments:
            symbol: the symbol for the atom itself, e.g. "O", "Cl"
            bonds_to_left: an optional value indicating the number of bonds to the previous
                           atom in the chain
            identiifer: an identifier for other atoms to refer to this one
    """
    def __init__(self, symbol: str, bonds_to_left: int = 0, identifier: str = "") -> None:
        self.symbol = symbol
        self._available_bonds = _STARTING_BONDS.get(symbol, -1)
        self.bonds_to_left = bonds_to_left
        self.bonds_to_right = 0
        self.identifier = identifier
        self.branches = []  # type: List[List[Atom]]
        self.references_out = []  # type: List[str]
        self.references_in = []  # type: List[Atom]

    @property
    def available_bonds(self) -> int:
        """ Returns the number of bonds currently used by hydrogen that are available
            for methylations and so on
        """
        # treat aromatic carbon rings differently for now, since it's really hard to split
        if self.symbol == "c":
            return 0
        available = self._available_bonds
        # reduce by the bond count to the previous atom in line
        available -= self.bonds_to_left
        # then to the next atom in line
        available -= self.bonds_to_right
        # then reduce by all of the bonds up from all branches
        available -= sum(atoms[0].bonds_to_left for atoms in self.branches)
        # and finally, reduce for any circular links
        available -= len(self.references_out) + len(self.references_in)
        return available

    def to_smiles(self) -> str:
        """ Returns a SMILES string to represent this atom and its branches,
            including any references to atoms not contained within those branches
        """
        bond_symbol = {3: "#", 2: '='}.get(self.bonds_to_left, "")
        smiles = [bond_symbol, self.symbol]
        if self.identifier:
            smiles.append(self.identifier)
        smiles.extend(self.references_out)
        if not self.branches:
            return "".join(smiles)
        for branch in self.branches:
            smiles.append("(%s)" % "".join(atom.to_smiles() for atom in branch))
        return "".join(smiles)

    def flatten(self) -> List["Atom"]:
        """ Returns a list of atoms in branches as if traversed depth-first,
            always begins with the current atom
        """
        res = [self]
        for branch in self.branches:
            for atom in branch:
                res.extend(atom.flatten())
        return res


class Bonds:
    """ Represents an entire SMILES string with multiple atoms.

        Arguments:
            smiles: the smiles string to build from
    """
    def __init__(self, smiles: str) -> None:
        self._atoms_by_identifier = {}  # type: Dict[str, Atom]
        self._top_level_atoms = self._parse_smiles(smiles)  # type: List[Atom]

    def _parse_smiles(self, smiles: str) -> List[Atom]:
        """ Parse a SMILES string into a list of Atom instances """

        def chain(smiles: str) -> Tuple[List[Atom], str]:
            """ Recursively parse branches in SMILES and return found Atoms and leftover SMILES """
            atoms = []  # type: List[Atom]
            while smiles:
                symbol = smiles[0]
                smiles = smiles[1:]

                # handle control structures
                if symbol == ")":
                    return atoms, smiles

                if symbol == "(":
                    subchains, smiles = chain(smiles)
                    atoms[-1].branches.append(subchains)
                    continue

                # ring identifiers/references
                if symbol.isdigit():
                    if not atoms or isinstance(atoms[-1], list):
                        raise ValueError("invalid smiles: %s" % (symbol + smiles))
                    # if it's new, then it's a declaration
                    if symbol not in self._atoms_by_identifier:
                        atoms[-1].identifier = symbol
                        self._atoms_by_identifier[symbol] = atoms[-1]
                        continue
                    # otherwise it's a reference
                    atoms[-1].references_out.append(symbol)
                    self._atoms_by_identifier[symbol].references_in.append(atoms[-1])
                    continue

                # multicharacter symbols
                if symbol == "C" and smiles and smiles[0] == "l":
                    symbol = "Cl"
                    smiles = smiles[1:]

                # different bond indicators
                current_bond = 1
                if symbol in "#=":
                    if symbol == "#":
                        current_bond = 3
                    elif symbol == "=":
                        current_bond = 2
                    symbol = smiles[0]
                    smiles = smiles[1:]
                if atoms:
                    atoms[-1].bonds_to_right = current_bond

                # finally, construct and add the atom
                atom = Atom(symbol, bonds_to_left=current_bond)
                atoms.append(atom)
            return atoms, smiles

        atoms, smiles = chain(smiles)
        # change the first atom in the smiles to have no bonds to left instead of 1
        atoms[0].bonds_to_left = 0
        assert not smiles, smiles
        return atoms

    def to_smiles(self) -> str:
        """ Build a new SMILES string from the current state """
        return "".join(atom.to_smiles() for atom in self._top_level_atoms)

    def __iter__(self) -> Iterator[Atom]:
        for atom in self._top_level_atoms:
            for item in atom.flatten():
                yield item


def gen_smiles_from_pksnrps(compound_pred: str) -> str:
    """ Generates the SMILES string for a specific compound prediction """
    smiles = ""

    if not compound_pred:
        return smiles

    residues = compound_pred.replace("(", "").replace(")", "").replace(" + ", " ").replace(" - ", " ").split()
    # Counts the number of malonate and its derivatives in polyketides
    mal_count = 0
    for residue in residues:
        if "mal" in residue:
            mal_count += 1

    # Reflecting reduction states of ketide groups starting at beta carbon of type 1 polyketide
    if residues[0] == "pk" and "mal" in residues[-1]:
        residues.pop(residues.index('pk')+1)
        residues.append('pks-end1')
    elif mal_count == len(residues):
        if residues[0] == "mal":
            residues[0] = "pks-start1"  # TODO why replace and not insert?
        if residues[-1] == "ccmal":
            residues.append('pks-end2')

    aa_smiles = load_smiles()

    for i, monomer in enumerate(residues):
        lower_monomer = monomer.lower()
        partial_lower = lower_monomer.split("-")[-1]
        if lower_monomer in aa_smiles:
            smiles_chunk = aa_smiles[lower_monomer]
        elif partial_lower in aa_smiles:
            smiles_chunk = aa_smiles[partial_lower]
        elif '|' in monomer:
            logging.debug("Substituting 'X' for combined monomer %r", monomer)
            smiles_chunk = aa_smiles['x']
        else:
            logging.debug("No SMILES mapping for unknown monomer %r", monomer)
            continue
        # trim the trailing O for all but the last smiles chunk
        if i < len(residues) - 1 and smiles_chunk.endswith("C(=O)O"):
            smiles_chunk = smiles_chunk[:-1]
        smiles += smiles_chunk
    return smiles


def load_smiles() -> Dict[str, str]:
    """Load smiles from a dictionary mapping residues to SMILES string"""
    aa_smiles = {}  # type: Dict[str, str]

    smiles_monomer = open(path.get_full_path(__file__, 'data', 'aaSMILES.txt'), 'r')

    for line in smiles_monomer.readlines():
        line = line.split("#", 1)[0].strip()
        if not line:
            continue
        smiles = line.split()
        assert len(smiles) == 2, "Invalid smiles line {!r}".format(line)
        assert smiles[0].lower() not in aa_smiles, "%s contained twice in smiles data" % smiles[0]
        aa_smiles[smiles[0].lower()] = smiles[1]

    smiles_monomer.close()
    return aa_smiles
