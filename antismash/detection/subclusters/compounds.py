# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Compound metadata for the subcluster detection module """

from dataclasses import dataclass
from typing import Optional

from antismash.common import path


@dataclass(frozen=True)
class CompoundInfo:
    """ Chemical metadata associated with a subcluster rule

        Attributes:
            name: the name of the compound
            smiles: the SMILES string of the compound, if one is known
            classification: the chemical classes the compound belongs to
    """
    name: str
    smiles: Optional[str]
    classification: list[str]


_COMPOUND_CACHE: dict[str, CompoundInfo] = {}


def _ensure_compounds_loaded() -> None:
    """ Loads the compound details from disk into the cache, if not already loaded """
    if _COMPOUND_CACHE:
        return
    filename = path.get_full_path(__file__, "data", "compound_details.txt")
    _COMPOUND_CACHE.update(_read_compounds(filename))


def get_compound(rule_name: str) -> CompoundInfo:
    """ Returns the compound metadata for a single rule, loading the details
        from disk on the first call

        Every rule is expected to have a matching entry in the details file.

        Arguments:
            rule_name: the name of the detection rule

        Returns:
            the compound associated with the rule
    """
    _ensure_compounds_loaded()
    return _COMPOUND_CACHE[rule_name]


def _read_compounds(detail_file: str) -> dict[str, CompoundInfo]:
    """ Parses a compound details file into a mapping of rule name to compound

        The file is expected to be tab separated, each row being a single
        compound with the columns: rule name, compound name, SMILES string, and
        a semicolon separated classification. Lines starting with '#' are
        treated as comments and ignored.

        Arguments:
            detail_file: the path of the file to parse

        Returns:
            a dictionary mapping rule name to the relevant compound
    """
    compounds: dict[str, CompoundInfo] = {}
    with open(detail_file, encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                raise ValueError(f"Invalid line in {detail_file}: {line!r}")
            rule_name, name, smiles, classification_raw = parts[:4]
            compounds[rule_name] = CompoundInfo(
                name=name,
                smiles=smiles or None,
                classification=[c.strip() for c in classification_raw.split(";") if c.strip()],
            )
    return compounds
