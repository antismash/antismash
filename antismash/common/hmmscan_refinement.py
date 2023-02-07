# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" A collection of functions for refining hmm domain fragments into best hits.
"""

from collections import defaultdict
from typing import Any, Dict, Iterable, List, Set, Tuple, Union

from Bio.SearchIO import QueryResult
from Bio.SearchIO._model.hsp import HSP


class HMMResult:
    """ A variant of HSP that allows for operations between multiple instances
        along with simplified operations e.g. len()
    """
    __slots__ = ["_hit_id", "_query_start", "_query_end", "_evalue", "_bitscore", "_internal_hits"]

    def __init__(self, hit_id: str, start: int, end: int, evalue: float,
                 bitscore: float, *, internal_hits: Iterable["HMMResult"] = None) -> None:
        self._hit_id = hit_id
        self._query_start = int(start)
        self._query_end = int(end)
        self._evalue = float(evalue)
        self._bitscore = float(bitscore)
        self._internal_hits: List[HMMResult] = []
        if internal_hits is not None:
            self.add_internal_hits(internal_hits)

    @property
    def hit_id(self) -> str:
        """ Returns the name of the matching profile """
        return self._hit_id

    @property
    def query_start(self) -> int:
        """ Returns the start position within the query's translation """
        return self._query_start

    @property
    def query_end(self) -> int:
        """ Returns the end position within the query's translation """
        return self._query_end

    @property
    def evalue(self) -> float:
        """ Returns the e-value of the hit """
        return self._evalue

    @property
    def bitscore(self) -> float:
        """ Returns the bitscore of the hit """
        return self._bitscore

    @property
    def internal_hits(self) -> Tuple["HMMResult", ...]:
        """ Any hits contained by this hit """
        return tuple(self._internal_hits)

    @property
    def detailed_names(self) -> List[str]:
        """ A list of hit ids, one for each depth of internal hit, stopping
            at the first depth for which there are no hits or multiple hits.
            Includes the current instance as the first element.
        """
        names = [self.hit_id]
        hits = self._internal_hits
        while len(hits) == 1:
            names.append(hits[0].hit_id)
            hits = hits[0]._internal_hits
        return names

    def add_internal_hits(self, hits: Iterable["HMMResult"]) -> None:
        """ Add hits within this hit """
        for hit in hits:
            if not hit.overlaps_with(self):
                raise ValueError(f"subdomain is not co-located with parent: {hit}, {self}")
        self._internal_hits.extend(hits)

    def __len__(self) -> int:
        return self.query_end - self.query_start

    def merge(self, other: "HMMResult") -> "HMMResult":
        """ Creates a new HMMResult instance from this instance and the
            provided instance.
        """
        assert self.hit_id == other.hit_id
        if self.query_start < other.query_start:
            start, end = self.query_start, other.query_end
        else:
            start, end = other.query_start, self.query_end
        return HMMResult(self.hit_id, start, end,
                         min(self.evalue, other.evalue),
                         max(self.bitscore, other.bitscore))

    def is_contained_by(self, other: "HMMResult") -> bool:
        """ Returns True if this instance is contained within the provided instance """
        if not isinstance(other, HMMResult):
            return False
        return other.query_start <= self.query_start < self.query_end <= other.query_end

    def overlaps_with(self, other: "HMMResult") -> bool:
        """ Returns True if this instance overlaps with the provided instance """
        if not isinstance(other, HMMResult):
            return False
        return other.query_end > self.query_start and self.query_end > other.query_start

    def to_json(self) -> Dict[str, Union[str, int, float]]:
        """ Converts the instance into a dictionary for use in json formats """
        data = {key.lstrip("_"): getattr(self, key) for key in self.__slots__}
        internal = data.pop("internal_hits")
        if internal:
            data["internal_hits"] = [hit.to_json() for hit in internal]
        return data

    @staticmethod
    def from_json(data: Dict[str, Any]) -> "HMMResult":
        """ Rebuilds a HMMResult instance from a JSON representation """
        internal_hits = [HMMResult.from_json(hit) for hit in data.get("internal_hits", [])]
        return HMMResult(str(data["hit_id"]), int(data["query_start"]), int(data["query_end"]),
                         float(data["evalue"]), float(data["bitscore"]), internal_hits=internal_hits)

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        subtypes = ""
        if self.detailed_names[1:]:
            subtypes = f", subtypes=[{', '.join(self.detailed_names[1:])}]"
        return (f"HMMResult({self.hit_id}, {self.query_start}, {self.query_end}, "
                f"evalue={self.evalue:g}, bitscore={self.bitscore}{subtypes})")

    def __eq__(self, other: Any) -> bool:
        if not isinstance(other, HMMResult):
            return False
        for attr in self.__slots__:
            if getattr(self, attr) != getattr(other, attr):
                return False
        return True

    def __hash__(self) -> int:
        return hash((self._hit_id, self._query_start, self._query_end, self._evalue, self._bitscore))


def _remove_overlapping(results: List[HMMResult], hmm_lengths: Dict[str, int]) -> List[HMMResult]:
    """ Strip domain list of overlapping domains,
        only keeping those with the highest scores

        Domains with an overlap of 20% or less aren't considered to be overlapping
    """
    non_overlapping = [results[0]]
    for result in results[1:]:
        previous = non_overlapping[-1]
        maxoverlap = 0.20 * max([hmm_lengths[result.hit_id],
                                 hmm_lengths[previous.hit_id]])
        if result.query_start < (previous.query_end - maxoverlap):
            # if the current result scores higher, replace the previous one
            if result.bitscore > previous.bitscore:
                non_overlapping[-1] = result
        else:
            non_overlapping.append(result)
    return non_overlapping


def _remove_incomplete(domains: List[HMMResult], hmm_lengths: Dict[str, int],
                       threshold: float = 0.5, fallback: float = 1./3.) -> List[HMMResult]:
    """ Removes all incomplete fragments for a domain type that are less than
        the threshold. If this would remove all hits for the domain type, then
        use the fallback as the threshold for the largest incomplete threshold.
        Will only strip all fragments if the largest fragment of this type is
        smaller than both threshold and fallback threshold.
    """
    assert fallback <= threshold
    complete = []
    for domain in domains:
        domainlength = hmm_lengths[domain.hit_id]
        if len(domain) > (threshold * domainlength):
            complete.append(domain)
    if complete:
        return complete

    # if none matched, just take the longest hit over the fallback size
    longest = 0.
    longest_index = 0
    for i, domain in enumerate(domains):
        domain_length = hmm_lengths[domain.hit_id]
        proportional_length = len(domain) / domain_length
        if proportional_length > longest:
            longest = proportional_length
            longest_index = i

    if longest > fallback:
        return [domains[longest_index]]

    # if still none, take a regulator if one exists
    for domain in domains:
        if "regulator" in domain.hit_id:
            return [domain]
    # ran out of fallbacks, return nothing
    return []


def _merge_domain_list(domains: List[HMMResult], hmm_lengths: Dict[str, int]) -> List[HMMResult]:
    """ Merges domains of the same kind if they would not be too long """
    categories: Dict[str, List[HMMResult]] = defaultdict(list)
    for domain in domains:
        categories[domain.hit_id].append(domain)
    remaining = []
    for category in categories.values():
        merged = category[0]
        max_span = 1.5 * hmm_lengths[merged.hit_id]  # TODO: use a more specific check
        for other in category[1:]:
            # only merge if the hit spans of the two domains are small enough
            if other.query_end - merged.query_start < max_span:
                merged = merged.merge(other)
            else:
                merged = other
        remaining.append(merged)
    return sorted(remaining, key=lambda result: result.query_start)


def _merge_immediate_neigbours(domains: List[HMMResult], hmm_lengths: Dict[str, int]) -> List[HMMResult]:
    result = [domains[0]]
    for domain in domains[1:]:
        if domain.hit_id != result[-1].hit_id:
            result.append(domain)
            continue
        # only merge if the hit spans of the two domains are small enough
        if domain.query_end - result[-1].query_start < 1.5 * hmm_lengths[domain.hit_id]:
            result[-1] = result[-1].merge(domain)
        else:
            result.append(domain)
    return result


def gather_by_query(results: List[HSP]) -> Dict[str, Set[HMMResult]]:
    """ Generates a mapping of query id to all HMMResults for that query

        Arguments:
            results: a list of HSP fragments as parsed by Bio's SearchIO

        Returns:
            a dictionary mapping query gene id to a set of HMMResults within that
            gene
    """
    results_by_id: Dict[str, Set[HMMResult]] = defaultdict(set)
    for result in results:
        for hsp in result.hsps:
            results_by_id[hsp.query_id].add(HMMResult(hsp.hit_id, hsp.query_start,
                                                      hsp.query_end, hsp.evalue,
                                                      hsp.bitscore))
    return results_by_id


def refine_hmmscan_results(hmmscan_results: List[QueryResult], hmm_lengths: Dict[str, int],
                           neighbour_mode: bool = False) -> Dict[str, List[HMMResult]]:
    """ Processes a list of QueryResult objects (from SearchIO.parse(..., 'hmmer3-text'))
            - merges domain fragments of the same ID
            - keeps only best hits from overlaps
            - removes incomplete domains

        Arguments:
            hmmscan_results: a list of QueryResult objects from Bio's SearchIO
            hmm_lengths: a dictionary mapping hmm id to length
            neighbour_mode: if on, does overlap removal before merge and merges
                            only when the next result has the same hit_id

        Returns:
            a mapping of gene name to list of HMMResults
    """
    results_by_id = gather_by_query(hmmscan_results)
    refined_results = {}
    for cds, results in results_by_id.items():
        refined = sorted(list(results), key=lambda result: result.query_start)
        if neighbour_mode:
            # Only keep best hits for overlapping domains
            refined = _remove_overlapping(refined, hmm_lengths)
            # Merge domain fragments which are really one domain
            refined = _merge_immediate_neigbours(refined, hmm_lengths)
        else:
            # Merge domain fragments which are really one domain
            refined = _merge_domain_list(refined, hmm_lengths)
            # Only keep best hits for overlapping domains
            refined = _remove_overlapping(refined, hmm_lengths)
        # Remove incomplete domains (covering less than 60% of total domain hmm length)
        refined = _remove_incomplete(refined, hmm_lengths)
        if refined:
            refined_results[cds] = refined

    return refined_results
