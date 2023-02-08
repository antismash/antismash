# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Generation of structure predictions for product compounds in both SMILES and
    image format, where possible
"""

import logging
from typing import Dict, Iterator, List, Tuple

from antismash.common import path

_STARTING_BONDS = {
    "C": 4,
    "c": 4,
    "O": 2,
    "N": 3,
}


def _load_smiles() -> Dict[str, str]:
    """Load smiles from a dictionary mapping residues to SMILES string"""
    aa_smiles: Dict[str, str] = {}

    with open(path.get_full_path(__file__, "data", "aaSMILES.txt"), "r", encoding="utf-8") as handle:
        for line in handle:
            line = line.split("#", 1)[0].strip()
            if not line:
                continue
            smiles = line.split()
            assert len(smiles) == 2, f"Invalid smiles line {line!r}"
            assert smiles[0].lower() not in aa_smiles, f"{smiles[0]} contained twice in smiles data"
            aa_smiles[smiles[0].lower()] = smiles[1]

    return aa_smiles


_SMILES = _load_smiles()


def get_all_smiles() -> Dict[str, str]:
    """Get smiles from the cached dictionary"""
    return _SMILES


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
        self.branches: List[List[Atom]] = []
        self.references_out: List[str] = []
        self.references_in: List[Atom] = []

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
            smiles.append(f"({''.join(atom.to_smiles() for atom in branch)})")
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
        self._atoms_by_identifier: Dict[str, Atom] = {}
        self._top_level_atoms: List[Atom] = self._parse_smiles(smiles)

    def _parse_smiles(self, smiles: str) -> List[Atom]:
        """ Parse a SMILES string into a list of Atom instances """

        def chain(smiles: str) -> Tuple[List[Atom], str]:
            """ Recursively parse branches in SMILES and return found Atoms and leftover SMILES """
            atoms: List[Atom] = []
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
                        raise ValueError(f"invalid smiles: {symbol + smiles}")
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


def methylate(smiles: str, variant: str) -> str:
    """ Add a methyl group to the given smiles at the first available position.
        Methylations of the first and last atom will not happen, they are considered
        the backbone and the methylation should occur on the sidechain.

        If an appropriate attachment point cannot be found, the smiles will be
        returned unchanged.

        Arguments:
            smiles: the SMILES string to methylate
            variant: which methylation variant to use, allowable values: C, N, O

        Returns:
            the newly methylated smiles, or the original smiles if no methylation
            is possible
    """
    if variant not in "CNO":
        raise ValueError(f"expected methylation variant of C, N, or O, not {variant}")
    bonds = Bonds(smiles)
    atoms = list(bonds)
    start = 1
    if atoms[0].symbol == "N" and atoms[1].symbol == "C":
        start = 2  # skip the C too
    for atom in atoms[start:-1]:
        if atom.symbol == variant and atom.available_bonds > 0:
            atom.branches.append([Atom("C")])
            break
    return bonds.to_smiles()


def gen_smiles_from_pksnrps(components: List[Tuple[str, str, List[str]]]) -> str:
    """ Generates the SMILES string for a specific compound prediction

        Arguments:
            a list of tuples, each tuple containing
                predicated substrate
                predicated monomer built from the substrate
                a list of domain names

        Returns:
            a SMILES string
    """
    if not components:
        return ""

    chunks = []

    def get_smiles_chunk(substrate: str, monomer: str, domains: List[str]) -> str:
        """ Fetch (and modify if relevant) the SMILES for a particular monomer/substrate """
        # if smiles exist for the monomer, use them by preference
        if monomer.lower() in _SMILES:
            return _SMILES[monomer.lower()]

        # if the substrate doesn't exist in the smiles set, abort
        if substrate.lower() not in _SMILES:
            logging.debug("No SMILES mapping for unknown monomer %r", monomer)
            return ""

        # otherwise take the substrate and apply known modifications
        smiles_chunk = _SMILES[substrate.lower()]
        for domain in domains:
            for prefix in "cno":
                if domain == f"{prefix}MT":
                    smiles_chunk = methylate(smiles_chunk, prefix.upper())
        return smiles_chunk

    # trim the trailing O for all but the last smiles chunk
    for substrate, monomer, modifications in components[:-1]:
        smiles_chunk = get_smiles_chunk(substrate, monomer, modifications)

        if smiles_chunk.endswith("C(=O)O"):
            smiles_chunk = smiles_chunk[:-1]
        chunks.append(smiles_chunk)

    chunks.append(get_smiles_chunk(*components[-1]))

    # if the last chunk was a PKS, add the end
    if components[-1][1].endswith("mal"):
        chunks.append(_SMILES["pks-end2"])

    return "".join(chunks)
