# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

from collections import defaultdict

def _strip_overlapping(results, hmm_lengths):
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

def _remove_incomplete(domainlist, hmm_lengths):
    complete = []
    for i in domainlist:
        domainlength = hmm_lengths[i.hit_id]
        if len(i) > (0.5 * domainlength):
            complete.append(i)
    return complete

def _merge_domain_list(domainlist, hmm_lengths):
    merged = [domainlist[0]]
    for current in domainlist[1:]:
        previous = merged[-1]
        domainlength = hmm_lengths[current.hit_id]
        if current.hit_id == previous.hit_id \
                and len(previous) + len(current) < 1.5 * domainlength:
            merged[-1] = previous.merge(current)
        else:
            merged.append(current)
    return merged

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
                         self.evalue * other.evalue,
                         self.bitscore + other.bitscore)

    def to_json(self):
        return {key : getattr(key) for key in self.__slots__}

def refine_hmmscan_results(hmmscan_results, hmm_lengths):
    """ Processes a list of QueryResult objects (from SearchIO.parse(..., 'hmmer3-text'))
            - keeps only best hits from overlaps
            - merges domain fragments of the same ID
            - removes incomplete domains

        Returns a mapping of gene name to HMMResult
    """
    results_by_id = defaultdict(set)
    for results in hmmscan_results:
        for hsp in results.hsps:
            results_by_id[hsp.query_id].add(HMMResult(hsp.hit_id, hsp.query_start,
                                                      hsp.query_end, hsp.evalue,
                                                      hsp.bitscore))
    refined_results = {}
    for cds, results in results_by_id.items():
        refined = sorted(list(results), key=lambda result: result.query_start)
        #Only keep best hits for overlapping domains
        refined = _strip_overlapping(refined, hmm_lengths)
        #Merge domain fragments which are really one domain
        refined = _merge_domain_list(refined, hmm_lengths)
        #Remove incomplete domains (covering less than 60% of total domain hmm length)
        refined = _remove_incomplete(refined, hmm_lengths)

        refined_results[cds] = refined

    return refined_results
