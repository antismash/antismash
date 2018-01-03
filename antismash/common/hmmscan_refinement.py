# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" A collection of functions for refining hmm domain fragments into best hits.
"""

from collections import defaultdict
from typing import Dict, Set, Union


def _remove_overlapping(results, hmm_lengths):
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


def _remove_incomplete(domains, hmm_lengths, threshold=0.5, fallback=1./3.):
    complete = []
    for i in domains:
        domainlength = hmm_lengths[i.hit_id]
        if len(i) > (threshold * domainlength):
            complete.append(i)
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


def _merge_domain_list(domainlist):
    categories = defaultdict(list)
    for domain in domainlist:
        categories[domain.hit_id].append(domain)
    remaining = []
    for category in categories.values():
        merged = category[0]
        for other in category[1:]:
            merged = merged.merge(other)
        remaining.append(merged)
    return sorted(remaining, key=lambda result: result.query_start)


def _merge_immediate_neigbours(domains):
    result = [domains[0]]
    for domain in domains[1:]:
        if domain.hit_id != result[-1].hit_id:
            result.append(domain)
            continue
        result[-1] = result[-1].merge(domain)
    return result


class HMMResult:
    """ A variant of HSP that allows for operations between multiple instances
        along with simplified operations e.g. len()
    """
    __slots__ = ["hit_id", "query_start", "query_end", "evalue", "bitscore"]

    def __init__(self, hit_id: str, start: int, end: int, evalue: float,
                 bitscore: float) -> None:
        self.hit_id = hit_id
        self.query_start = int(start)
        self.query_end = int(end)
        self.evalue = float(evalue)
        self.bitscore = float(bitscore)

    def __len__(self):
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

    def to_json(self) -> Dict[str, Union[str, int, float]]:
        """ Converts the instance into a dictionary for use in json formats """
        return {key: getattr(self, key) for key in self.__slots__}

    def __repr__(self):
        return str(self)

    def __str__(self):
        return "HMMResult(%s, %d, %d, evalue=%g, bitscore=%g)" % (self.hit_id,
                   self.query_start, self.query_end, self.evalue, self.bitscore)


def gather_by_query(results) -> Dict[str, Set[HMMResult]]:
    """ Generates a mapping of query id to all HMMResults for that query """
    results_by_id = defaultdict(set)  # type: Dict[str, Set[HMMResult]]
    for result in results:
        for hsp in result.hsps:
            results_by_id[hsp.query_id].add(HMMResult(hsp.hit_id, hsp.query_start,
                                                      hsp.query_end, hsp.evalue,
                                                      hsp.bitscore))
    return results_by_id


def refine_hmmscan_results(hmmscan_results, hmm_lengths, neighbour_mode=False):
    """ Processes a list of QueryResult objects (from SearchIO.parse(..., 'hmmer3-text'))
            - merges domain fragments of the same ID
            - keeps only best hits from overlaps
            - removes incomplete domains

        Arguments:
            hmmscan_results: a list of QueryResult objects
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
            refined = _merge_immediate_neigbours(refined)
        else:
            # Merge domain fragments which are really one domain
            refined = _merge_domain_list(refined)
            # Only keep best hits for overlapping domains
            refined = _remove_overlapping(refined, hmm_lengths)
        # Remove incomplete domains (covering less than 60% of total domain hmm length)
        refined = _remove_incomplete(refined, hmm_lengths)
        if refined:
            refined_results[cds] = refined

    return refined_results
