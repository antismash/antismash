# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" JSON-friendly classes explicitly for use by the javascript drawing libraries
"""

from typing import Any, Iterator, List, Tuple

from antismash.common.secmet.features import CDSFeature
from antismash.common.secmet.qualifiers import NRPSPKSQualifier


class JSONBase(dict):
    """ A base class for JSON-serialisable objects """
    def __init__(self, keys: List[str]) -> None:
        super().__init__()
        self._keys = keys

    def __getitem__(self, key: str) -> Any:
        return getattr(self, key)

    def items(self) -> Iterator[Tuple[str, Any]]:  # type: ignore
        for key in self._keys:
            yield (key, getattr(self, key))

    def values(self) -> Iterator[Any]:  # type: ignore
        for key in self._keys:
            yield getattr(self, key)

    def __len__(self) -> int:
        return len(self._keys)


class JSONDomain(JSONBase):
    """ A JSON-serialisable object for simplifying domain datatypes throughout this file """
    def __init__(self, domain: NRPSPKSQualifier.Domain, predictions: List[Tuple[str, str]], napdos_link: str,
                 blast_link: str, sequence: str, dna: str) -> None:
        super().__init__(['type', 'start', 'end', 'predictions', 'napdoslink',
                          'blastlink', 'sequence', 'dna_sequence'])
        self.type = str(domain.name)
        self.start = int(domain.start)
        self.end = int(domain.end)
        self.predictions = predictions
        self.napdoslink = str(napdos_link)
        self.blastlink = str(blast_link)
        self.sequence = str(sequence)
        self.dna_sequence = str(dna)


class JSONOrf(JSONBase):
    """ A JSON-serialisable object for simplifying ORF datatypes throughout this file """
    def __init__(self, feature: CDSFeature) -> None:
        super().__init__(['id', 'sequence', 'domains'])
        self.sequence = feature.translation
        self.id = feature.get_name()
        self.domains = []  # type: List[JSONDomain]

    def add_domain(self, domain: JSONDomain) -> None:
        """ Add a JSONDomain to the list of domains in this ORF """
        assert isinstance(domain, JSONDomain)
        self.domains.append(domain)
