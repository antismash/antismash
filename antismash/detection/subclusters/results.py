# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Result handling for the subcluster detection module """

import logging
from typing import Any, Optional, Self

from antismash.common.hmm_rule_parser.cluster_prediction import RuleDetectionResults
from antismash.common.module_results import DetectionResults
from antismash.common.secmet import Record, Region, SubRegion

from .predictions import SubclusterPrediction, build_predictions
from .subregions import SubRegionMode, build_subregions


class SubclusterDetectionResults(DetectionResults):
    """ Results for the subcluster detection module """

    schema_version = 1  # increment when the JSON format changes

    def __init__(  # pylint: disable=too-many-arguments
            self,
            record_id: str,
            rule_results: RuleDetectionResults,
            rule_names: set[str],
            strictness: str,
            predictions: list[SubclusterPrediction],
            subregions: list[SubRegion],
            subregion_mode: SubRegionMode,
    ) -> None:
        super().__init__(record_id)
        self.rule_results = rule_results
        self.rule_names = rule_names
        self.strictness = strictness
        self.subregion_mode = subregion_mode
        self.predictions = predictions
        self.subregions = subregions

    def get_predictions_for_region(self, region: Region) -> list[SubclusterPrediction]:
        """ Returns all predictions overlapping the given region

            Arguments:
                region: the region to find predictions for

            Returns:
                a list of predictions, which may be empty
        """
        return [
            prediction for prediction in self.predictions
            if region.overlaps_with(prediction.location)
        ]

    def get_predicted_subregions(self) -> list[SubRegion]:
        """ Returns the subregion features for the detected subclusters

            These are built when detection runs, since they depend on the areas
            found by the other detection modules, and are added to the record
            during the region-formation step of the main pipeline.
        """
        return self.subregions

    def to_json(self) -> dict[str, Any]:
        """ Returns a JSON-ready representation of the instance """
        return {
            "schema_version": self.schema_version,
            "record_id": self.record_id,
            "strictness": self.strictness,
            "subregion_mode": self.subregion_mode,
            "rule_names": sorted(self.rule_names),
            "rule_results": self.rule_results.to_json(),
        }

    @classmethod
    def from_json(cls, json: dict[str, Any], record: Record) -> Optional[Self]:
        """ Regenerates the results from a JSON representation

            The results will be discarded if either the schema version of the
            results themselves or that of the contained rule detection results
            no longer matches.

            Arguments:
                json: the JSON representation to rebuild from
                record: the record the results were generated for

            Returns:
                the rebuilt instance, or None if the previous results are no
                longer valid
        """
        if json.get("schema_version") != cls.schema_version:
            logging.debug(
                "Discarding subcluster results: schema version %s != %s",
                json.get("schema_version"), cls.schema_version,
            )
            return None

        rule_results = RuleDetectionResults.from_json(json["rule_results"], record)
        if rule_results is None:
            logging.debug("Discarding subcluster results: rule detection result schema changed")
            return None

        # neither the predictions nor the subregions are stored, as they're
        # derived from the rule results and, for the subregions, the areas found
        # by the other detection modules, so both are rebuilt here
        subregion_mode = SubRegionMode(json["subregion_mode"])
        return cls(
            record_id=json["record_id"],
            rule_results=rule_results,
            rule_names=set(json["rule_names"]),
            strictness=json["strictness"],
            subregion_mode=subregion_mode,
            predictions=build_predictions(rule_results, json["strictness"]),
            subregions=build_subregions(record, rule_results.protoclusters,
                                        tool=rule_results.tool, mode=subregion_mode),
        )
