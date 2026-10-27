# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Construction of the view objects describing each detected subcluster

    Predictions are derived from the rule detection results, so they aren't
    stored with the results and are rebuilt when the results are reconstructed.
"""
from collections import defaultdict
from dataclasses import dataclass
from typing import Optional

from antismash.common.hmm_rule_parser.rule_parser import DetectionRule
from antismash.common.hmm_rule_parser.cluster_prediction import CDSResults, RuleDetectionResults
from antismash.common.secmet.locations import Location

from .compounds import CompoundInfo, get_compound
from .ruleset import get_ruleset
from .signatures import get_signatures


@dataclass(frozen=True)
class CDSDomainHit:
    """ A single signature match found within a specific CDS

        Attributes:
            domain_name: the name of the matching signature
            domain_description: a description of the signature, if it has one
            domain_accession: the accession of the signature, if it has one
            cds_name: the name of the CDS the match was found in
            evalue: the E-value of the match
            bitscore: the bitscore of the match
    """
    domain_name: str
    domain_description: Optional[str]
    domain_accession: Optional[str]
    cds_name: str
    evalue: float
    bitscore: float


class SubclusterPrediction:
    """ A single predicted subcluster, used as a view object for HTML rendering

        Attributes:
            rule: the detection rule whose conditions were met
            location: the location of the protocluster that produced this prediction,
                      i.e. the matching core along with any neighbourhood the rule defines
            compound: the compound associated with the detection rule
            domain_hits: a flat list of every domain hit, each paired with its CDS name
    """

    def __init__(
            self,
            *,
            rule: DetectionRule,
            location: Location,
            compound: CompoundInfo,
            domain_hits: list[CDSDomainHit],
    ) -> None:
        self.rule = rule
        self.location = location
        self.compound = compound
        self.domain_hits = domain_hits

    @property
    def rule_name(self) -> str:
        """ The name of the matching detection rule """
        return self.rule.name

    @property
    def conditions_str(self) -> str:
        """ The conditions of the matching rule, without the outermost brackets """
        text = str(self.rule.conditions)
        if text.startswith("(") and text.endswith(")"):
            return text[1:-1]
        return text

    @property
    def domain_hits_by_cds(self) -> dict[str, list[CDSDomainHit]]:
        """ The domain hits of the subcluster, grouped by the name of the CDS
            they were found in
        """
        hits_by_cds: dict[str, list[CDSDomainHit]] = defaultdict(list)
        for hit in self.domain_hits:
            hits_by_cds[hit.cds_name].append(hit)
        return dict(hits_by_cds)

    @property
    def domain_hits_by_domain(self) -> dict[str, list[CDSDomainHit]]:
        """ The domain hits of the subcluster, grouped by the name of the
            matching signature
        """
        hits_by_domain: dict[str, list[CDSDomainHit]] = defaultdict(list)
        for hit in self.domain_hits:
            hits_by_domain[hit.domain_name].append(hit)
        return dict(hits_by_domain)

    def extends_beyond(self, location: Location) -> tuple[bool, bool]:
        """ Returns whether the subcluster extends past the start and the end of
            the given location, e.g. that of a region

            Arguments:
                location: the location to compare to

            Returns:
                a tuple of whether the subcluster extends past the start and
                whether it extends past the end of the location
        """
        # the first and last parts are the sides of a location crossing the origin,
        # and the end of a location is exclusive, so the last base is end - 1
        return (self.location.parts[0].start not in location,
                self.location.parts[-1].end - 1 not in location)

    def __repr__(self) -> str:
        return (
            f"SubclusterPrediction(rule_name={self.rule_name!r}, "
            f"location={self.location.start}-{self.location.end}, "
            f"cds_count={len(self.domain_hits_by_cds)})"
        )


def _build_domain_hits(rule_name: str, cds_results: list[CDSResults]) -> list[CDSDomainHit]:
    """ Builds the domain hits that fired the given rule, flattened across the
        given CDS results

        Arguments:
            rule_name: the name of the rule the hits must have contributed to
            cds_results: the CDS results to gather the hits from

        Returns:
            a list of domain hits
    """
    signatures = get_signatures()
    hits: list[CDSDomainHit] = []
    for cds_result in cds_results:
        cds_name = cds_result.cds.get_name()
        fired = cds_result.definition_domains.get(rule_name, set())
        domains_by_name = {domain.name: domain for domain in cds_result.domains}
        for domain_name in sorted(fired):
            domain = domains_by_name[domain_name]
            signature = signatures[domain_name]
            hits.append(CDSDomainHit(
                domain_name=domain_name,
                domain_description=signature.description,
                domain_accession=signature.accession,
                cds_name=cds_name,
                evalue=domain.evalue,
                bitscore=domain.bitscore,
            ))
    return hits


def build_predictions(rule_results: RuleDetectionResults, strictness: str,
                      ) -> list[SubclusterPrediction]:
    """ Builds the view objects for each detected subcluster

        Arguments:
            rule_results: the detection results to build the predictions from
            strictness: the strictness level detection was run at, which determines
                        the ruleset the matching rules are taken from

        Returns:
            a list of predictions, one per detected subcluster
    """
    ruleset = get_ruleset(strictness)
    predictions: list[SubclusterPrediction] = []
    for protocluster, cds_results in rule_results.cds_by_cluster.items():
        rule = ruleset.get_rule_by_name(protocluster.product)
        predictions.append(SubclusterPrediction(
            rule=rule,
            location=protocluster.location,
            compound=get_compound(rule.name),
            domain_hits=_build_domain_hits(rule.name, cds_results),
        ))
    return predictions
