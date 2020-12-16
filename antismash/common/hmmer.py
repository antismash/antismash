# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Common functionality for finding and marking PFAMDomains within a record. """

from collections import defaultdict
from dataclasses import dataclass
import logging
import os
from typing import Any, Dict, Iterable, List, Optional, Set, Tuple

from antismash.common import fasta, module_results, path, pfamdb, subprocessing
from antismash.common.secmet import Record, CDSFeature
from antismash.common.secmet.features import FeatureLocation, PFAMDomain
from antismash.common.secmet.locations import location_from_string
from antismash.config import get_config


@dataclass(frozen=True)
class HmmerHit:  # pylint: disable=too-many-instance-attributes
    """ A class containing information about a HMMer hit

        Protein locations are expected to be python slice-style coordinates,
        with an inclusive start and an exclusive end.
    """
    location: str
    label: str
    locus_tag: str
    domain: str
    evalue: float
    score: float
    identifier: str
    description: str
    protein_start: int
    protein_end: int
    translation: str

    def __post_init__(self) -> None:
        if self.protein_start >= self.protein_end:
            raise ValueError("HMMer hit has inverted start and end: %r" % self)
        if len(self.translation) != len(self):
            raise ValueError("translation length does not match protein location size")

    def __len__(self) -> int:
        return self.protein_end - self.protein_start

    def to_json(self) -> Dict[str, Any]:
        """ Returns a JSON-ready representation of the instance """
        return dict(vars(self))

    @staticmethod
    def from_json(data: Dict[str, Any]) -> "HmmerHit":
        """ Reconstructs an instance from a JSON representation """
        return HmmerHit(**data)


class HmmerResults(module_results.ModuleResults):
    """ Results for hmmer-based detection """
    schema_version = 2

    def __init__(self, record_id: str, evalue: float, score: float,
                 database: str, tool: str, hits: List[HmmerHit]) -> None:
        super().__init__(record_id)
        self.hits = list(hits)
        self.evalue = float(evalue)
        self.score = float(score)
        self.database = str(database)
        self.tool = str(tool)

    def to_json(self) -> Dict[str, Any]:
        json = {"hits": [hit.to_json() for hit in self.hits],
                "record id": self.record_id,
                "schema": self.schema_version, "max evalue": self.evalue,
                "min score": self.score, "database": self.database,
                "tool": self.tool}
        return json

    @staticmethod
    def from_json(json: Dict[str, Any], record: Record) -> Optional["HmmerResults"]:
        """ Regenerate the results from JSON.
            If max_evalue or min_score aren't equal or narrower than those the
            results were generated with, the results will be discarded.
        """

        if record.id != json.get("record id"):
            logging.warning("Hmmer results are for different record, discarding previous results")
            return None

        if json.get("schema") != HmmerResults.schema_version:
            logging.warning("Hmmer results are for different result schema, discarding previous results")
            return None

        evalue = json.get("max evalue")
        score = json.get("min score")
        if evalue is None or score is None:
            raise ValueError("Invalid Hmmer result values")
        assert isinstance(score, float) and isinstance(evalue, float)

        hits = json.get("hits")
        if not isinstance(hits, list):
            raise TypeError("FullHmmer results contain unexpected types")
        hits = [HmmerHit(**hit) for hit in hits]

        return HmmerResults(record.id, evalue, score, json["database"], json["tool"], hits)

    def refilter(self, max_evalue: float, min_score: float) -> "HmmerResults":
        """ Trims the results to stricter thresholds for score and E-value """
        if max_evalue > self.evalue:
            raise ValueError("cannot refilter to a more lenient evalue: %s -> %s" % (self.evalue, max_evalue))
        if min_score < self.score:
            raise ValueError("cannot refilter to a more lenient score: %s -> %s" % (self.score, min_score))
        self.hits = [hit for hit in self.hits if hit.score >= min_score and hit.evalue <= max_evalue]
        self.evalue = max_evalue
        self.score = min_score
        return self

    def add_to_record(self, record: Record) -> None:
        """ Adds the hits as PFAMDomains to the given record """
        db_version = pfamdb.get_db_version_from_path(self.database)
        for i, hit in enumerate(self.hits):
            protein_location = FeatureLocation(hit.protein_start, hit.protein_end)
            pfam_feature = PFAMDomain(location_from_string(hit.location),
                                      description=hit.description, protein_location=protein_location,
                                      identifier=hit.identifier, tool=self.tool, locus_tag=hit.locus_tag)
            for key in ["label", "locus_tag", "domain", "evalue",
                        "score", "translation"]:
                setattr(pfam_feature, key, getattr(hit, key))
            pfam_feature.database = db_version
            pfam_feature.detection = "hmmscan"
            pfam_feature.domain_id = "{}_{}_{:04d}".format(self.tool, pfam_feature.locus_tag, i + 1)
            record.add_pfam_domain(pfam_feature)


def build_hits(record: Record, hmmscan_results: List, min_score: float,
               max_evalue: float, database: str) -> List[HmmerHit]:
    """ Builds HmmerHits from the given hmmscan results

        Arguments:
            record: the Record being scanned
            hmmscan_results: the results of Bio.SearchIO.parse
            min_score: a minimum allowable bitscore for hits (exclusive)
            max_evalue: a maximum allowable evalue for hits (exclusive)
            database: the name of the database used to find the hits

        Returns:
            a list of HmmerHit instances
    """
    hits = []
    feature_by_id = record.get_cds_name_mapping()

    for result in hmmscan_results:
        for hsp in result.hsps:
            if hsp.bitscore <= min_score or hsp.evalue >= max_evalue:
                continue

            feature = feature_by_id[hsp.query_id]
            location = feature.get_sub_location_from_protein_coordinates(hsp.query_start, hsp.query_end)

            hit = {"location": str(location),
                   "label": result.id, "locus_tag": feature.get_name(),
                   "domain": hsp.hit_id, "evalue": hsp.evalue, "score": hsp.bitscore,
                   "translation": feature.translation[hsp.query_start:hsp.query_end],
                   "identifier": pfamdb.get_pfam_id_from_name(hsp.hit_id, database),
                   "description": hsp.hit_description, "protein_start": hsp.query_start,
                   "protein_end": hsp.query_end,
                   }
            hits.append(HmmerHit(**hit))
    return hits


def remove_overlapping(hits: List[HmmerHit], cutoffs: Dict[str, float],
                       overlap_limit: int = 10) -> List[HmmerHit]:
    """ Filters a list of HmmerHits so that no overlapping hits remain.
        If two hits overlap by more than the specified limit, the best is kept
        based on the following, in order:
        - highest bitscore (normalised by profile cutoff)
        - longest hit
        - earliest start
        - identifier of the profile hit in ascending order

        Arguments:
            hits: the list of hits to filter
            cutoffs: a dictionary mapping profile identifier to cutoff as a float
            overlap_limit: the number of overlapping aminos required to be filtered

        Returns:
            a list of HmmerHits, sorted by protein start position
    """

    if not hits:
        assert 0
        return []

    try:
        normalised = {hit: cutoffs[hit.identifier] / hit.score for hit in hits}
    except KeyError as err:
        raise ValueError("cutoff mapping does not contain all hit identifiers") from err

    def ranking_stats(hit: HmmerHit) -> Tuple[float, float, int, str]:
        """ In order of importance:
                highest normalised score, longest hit, earliest start, identifier

            Score and hit are inverted in order for all components to sort the
            same way. An ascending sort will have the highest scoring hit first.

            Hits with the same values for the above are the same hit and can
            ranked equivalently with no loss of information.
        """
        return normalised[hit], 1/len(hit), hit.protein_start, hit.identifier

    # ensure ordering
    hits = sorted(hits, key=lambda hit: hit.protein_start)

    # build sets of overlapping hits, no two sets should overlap
    groups: List[Set[HmmerHit]] = []
    current = {hits[0]}
    max_current = hits[0].protein_end  # the end point of the set for simple checking
    for hit in hits:
        if max_current - overlap_limit < hit.protein_start:
            groups.append(current)
            current = {hit}
            max_current = hit.protein_end
            continue
        current.add(hit)
        max_current = max(max_current, hit.protein_end)
    # not forgetting to add the in-progress set
    groups.append(current)

    # the filter stage
    cleaned = []
    for group in groups:
        best_of: List[HmmerHit] = []
        # start with the highest score and work down
        for hit in sorted(group, key=ranking_stats):
            for other in best_of:
                start = other.protein_start + overlap_limit
                end = other.protein_end - overlap_limit
                # abort if its overlap with a higher ranked hit
                if hit.protein_start <= end and hit.protein_end >= start:
                    break
            else:  # no overlap with existing best
                best_of.append(hit)
        cleaned.extend(best_of)
        assert best_of
    assert cleaned
    return sorted(cleaned, key=lambda hit: hit.protein_start)


def run_hmmer(record: Record, features: Iterable[CDSFeature], max_evalue: float,
              min_score: float, database: str, tool: str, filter_overlapping: bool = True,
              use_cut_tc: bool = True) -> HmmerResults:
    """ Build hmmer results for the given features

        Arguments:
            record: the Record instance to run hmmer over
            features: the list of CDSFeatures to run over specifically
            max_evalue: a maximum evalue allowed for hits (exclusive)
            min_evalue: a minimum evalue allowed for hits (exclusive)
            database: the database to search for hits within
            tool: the name of the specific tool calling into this module
            use_cut_tc: whether to use threshold cutoff as specified in profiles
    """
    if not os.path.exists(database):
        raise ValueError("Given database does not exist: %s" % database)
    query_sequence = fasta.get_fasta_from_features(features)
    opts: List[str] = []
    if use_cut_tc:
        opts.append('--cut_tc')
    hmmscan_results = subprocessing.run_hmmscan(database, query_sequence, opts=opts)
    hits = build_hits(record, hmmscan_results, min_score, max_evalue, database)
    if filter_overlapping:
        results_by_cds = defaultdict(list)
        for hit in hits:
            results_by_cds[hit.locus_tag].append(hit)
        cutoffs = pfamdb.get_pfam_cutoffs(database)
        hits = []
        for locus_hits in results_by_cds.values():
            hits.extend(remove_overlapping(locus_hits, cutoffs))
    return HmmerResults(record.id, max_evalue, min_score, database, tool, hits)


def ensure_database_pressed(filepath: str, return_not_raise: bool = False) -> List[str]:
    """ Ensures that the given HMMer database exists and that the hmmpress
        generated files aren't out of date.

        Arguments:
            filepath: the path to the HMMer database
            return_not_raise: whether to catch errors and return their messages as strings

        Returns:
            any encountered error messages, will never be populated without return_not_raise == True
    """
    components = ["{}{}".format(filepath, ext) for ext in ['.h3f', '.h3i', '.h3m', '.h3p']]

    if path.is_outdated(components, filepath):
        logging.info("%s components missing or obsolete, re-pressing database", filepath)
        if "hmmpress" not in get_config().executables:
            msg = "Failed to hmmpress {!r}: cannot find executable for hmmpress".format(filepath)
            if not return_not_raise:
                raise RuntimeError(msg)
            return [msg]

        result = subprocessing.run_hmmpress(filepath)
        if not result.successful():
            msg = "Failed to hmmpress {!r}: {}".format(filepath, result.stderr)
            if not return_not_raise:
                raise RuntimeError(msg)
            return [msg]
    return []
