# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" HMM signatures for the subcluster detection module """

from typing import Optional

from antismash.common import path
from antismash.common.signature import HmmSignature


# the description of every signature used by the module
DETAILS_FILE = path.get_full_path(__file__, "data", "hmmdetails.txt")
# the combined profile database built from those signatures by prepare_data()
AGGREGATE_HMM_FILE = path.get_full_path(__file__, "data", "subcluster_seeds.hmm")

_SIGNATURE_CACHE: dict[str, "SubclusterHmmSignature"] = {}


class SubclusterHmmSignature(HmmSignature):
    """ An HMM signature, extended with the accession of the source profile, if it has one """

    def __init__(self, name: str, description: str, cutoff: int,
                 hmm_path: str, seed_count: int = 0, *,
                 accession: Optional[str] = None) -> None:
        super().__init__(name, description, cutoff, hmm_path, seed_count)
        self.accession = accession


def _ensure_signatures_loaded() -> None:
    """ Loads the subcluster HMM signatures from disk into the cache, if not
        already loaded
    """
    if _SIGNATURE_CACHE:
        return
    signatures = _read_signatures(DETAILS_FILE)
    _SIGNATURE_CACHE.update((signature.name, signature) for signature in signatures)
    if len(_SIGNATURE_CACHE) != len(signatures):
        raise ValueError(f"Duplicate signature names in {DETAILS_FILE}")


def get_signatures() -> dict[str, SubclusterHmmSignature]:
    """ Returns all subcluster HMM signatures, loading them from disk on the
        first call

        Returns:
            a dictionary mapping signature name to the relevant signature
    """
    _ensure_signatures_loaded()
    return _SIGNATURE_CACHE


def _read_signatures(detail_file: str) -> list[SubclusterHmmSignature]:
    """ Generates subcluster HMM signatures from a file

        The file is expected to be tab separated, each row being a single HMM
        reference with the columns: name, description, minimum score cutoff,
        HMM path, and an optional accession. Paths in the file are assumed to be
        relative to the file itself. Lines starting with '#' are treated as comments and
        ignored.

        Arguments:
            detail_file: the path of the file to parse

        Returns:
            a list of SubclusterHmmSignatures
    """
    bad_lines: list[str] = []
    signatures: list[SubclusterHmmSignature] = []
    with open(detail_file, "r", encoding="utf-8") as data:
        for line in data.read().split("\n"):
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) not in (4, 5):
                bad_lines.append(line)
                continue
            name, desc, cutoff, filename = parts[:4]
            # an empty accession column is treated the same as a missing one
            accession = parts[4] if len(parts) == 5 and parts[4] else None
            signatures.append(SubclusterHmmSignature(
                name,
                desc,
                int(cutoff),
                path.get_full_path(detail_file, filename),
                accession=accession)
            )

    if bad_lines:
        raise ValueError("Invalid lines in HMM detail file (first 10):\n%s"
                         % "\n".join(bad_lines[:10]))

    return signatures
