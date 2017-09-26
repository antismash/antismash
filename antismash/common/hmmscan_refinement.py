# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

from collections import defaultdict
import logging

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


class HMMResult:
    __slots__ = ["hit_id", "query_start", "query_end", "evalue", "bitscore"]
    def __init__(self, hit_id, start, end, evalue, bitscore):
        self.hit_id = hit_id
        self.query_start = int(start)
        self.query_end = int(end)
        self.evalue = float(evalue)
        self.bitscore = float(bitscore)

    def __len__(self):
        return self.query_end - self.query_start

    def merge(self, other):
        assert self.hit_id == other.hit_id
        if self.query_start < other.query_start:
            start, end = self.query_start, other.query_end
        else:
            start, end = other.query_start, self.query_end
        return HMMResult(self.hit_id, start, end,
                         min(self.evalue, other.evalue),
                         max(self.bitscore, other.bitscore))

    def to_json(self):
        return {key : getattr(self, key) for key in self.__slots__}

    def __str__(self):
        return "HMMResult(%s, %d, %d, evalue=%g, bitscore=%g)" % (self.hit_id,
                   self.query_start, self.query_end, self.evalue, self.bitscore)

def gather_by_query(results):
    results_by_id = defaultdict(set)
    for result in results:
        for hsp in result.hsps:
            results_by_id[hsp.query_id].add(HMMResult(hsp.hit_id, hsp.query_start,
                                                      hsp.query_end, hsp.evalue,
                                                      hsp.bitscore))
    return results_by_id


def refine_hmmscan_results(hmmscan_results, hmm_lengths):
    """ Processes a list of QueryResult objects (from SearchIO.parse(..., 'hmmer3-text'))
            - merges domain fragments of the same ID
            - keeps only best hits from overlaps
            - removes incomplete domains

        Returns a mapping of gene name to HMMResult
    """
    results_by_id = gather_by_query(hmmscan_results)
    refined_results = {}
    for cds, results in results_by_id.items():
        refined = sorted(list(results), key=lambda result: result.query_start)
        #Merge domain fragments which are really one domain
        refined = _merge_domain_list(refined)
        #Only keep best hits for overlapping domains
        refined = _remove_overlapping(refined, hmm_lengths)
        #Remove incomplete domains (covering less than 60% of total domain hmm length)
        refined = _remove_incomplete(refined, hmm_lengths)
        if refined:
            refined_results[cds] = refined

    return refined_results
