# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.


"""Adapt antiSMASH inputs and results for fungal boundary analysis."""

from __future__ import annotations

from copy import deepcopy
from dataclasses import asdict
import logging
from typing import Any

from antismash.common import path
from antismash.common.html_renderer import FileTemplate, HTMLSections
from antismash.common.layers import OptionsLayer, RecordLayer, RegionLayer
from antismash.common.module_results import DetectionResults
from antismash.common.secmet import FeatureLocation, Record, SubRegion
from antismash.config import ConfigType
from antismash.config.args import ModuleArgs
from antismash.detection import DetectionStage

from . import fundaries as fundaries

NAME = "fungal_boundary"
SHORT_DESCRIPTION = "Predict fungal cluster boundaries using Fundaries"
DETECTION_STAGE = DetectionStage.AREA_REFINEMENT


class FungalBoundaryResults(DetectionResults):
    """Expose a Fundaries report through antiSMASH's detection results API."""

    schema_version = 5

    def __init__(self, report: fundaries.FundariesReport) -> None:
        """Keep the predictor's report as the source of all boundary data."""
        super().__init__(report.contig)
        self.report = report

    def to_json(self) -> dict[str, Any]:
        """
            Save Fundaries report as a json, including boundary predictions
            and gene rows, to nested dictionaries without modifying the report.

            Arguments:
                none

            Returns:
                a JSON-serializable dictionary containing the schema version,
                record identifier, and complete Fundaries report
        """
        return {
            "schema_version": self.schema_version,
            "record_id": self.record_id,
            "report": asdict(self.report),
        }

    @staticmethod
    def from_json(json: dict[str, Any],
                  record: Record) -> FungalBoundaryResults | None:
        """ Restore fungal boundary results saved by to_json().
            so to_json() --> from_json()

            Reconstruct the Fundaries report, boundary predictions, and gene rows
            without modifying the saved data.

            Arguments:
                json: saved results with the schema version, record
                      identifier, and Fundaries report
                record: the record whose identifier is used to validate the results

            Returns:
                the reconstructed results, or None
        """
        if json.get("schema_version") != FungalBoundaryResults.schema_version:
            return None
        if json.get("record_id") != record.id:
            return None
        report_data = deepcopy(json["report"])
        for prediction in report_data["predictions"]:
            prediction["motif_type"] = fundaries.MotifType(prediction["motif_type"])
            prediction["gene_table"] = [fundaries.GeneRow(**gene)
                                        for gene in prediction["gene_table"]]
        report = fundaries.FundariesReport(
            contig=report_data["contig"],
            predictions=[fundaries.BoundaryPrediction(**prediction)
                         for prediction in report_data["predictions"]],
        )
        return FungalBoundaryResults(report)

    def get_predicted_subregions(self) -> list[SubRegion]:
        """Create boundary features only when antiSMASH requests them."""
        return [
            SubRegion(FeatureLocation(prediction.start, prediction.end),
                      tool=NAME, label=prediction.anchor_id)
            for prediction in self.report.predictions
        ]


def get_arguments() -> ModuleArgs:
    """Register a separate enable switch for this new wrapper."""
    args = ModuleArgs("Fungal boundary options", "fungalnew")
    args.add_analysis_toggle("fungalnew", dest="fungalnew",
                             action="store_true", default=False,
                             help="Predict fungal cluster boundaries "
                                  "using Fundaries.")
    return args


def is_enabled(options: ConfigType) -> bool:
    """Return whether this wrapper is enabled."""
    return options.fungalnew


def check_options(options: ConfigType) -> list[str]:
    """Restrict the enabled module to fungal analysis."""
    if is_enabled(options) and options.taxon != "fungi":
        return ["Fungal boundary prediction requires fungal sequences."]
    return []


def check_prereqs(_options: ConfigType) -> list[str]:
    """Fundaries uses antiSMASH's existing Python dependencies."""
    return []


def regenerate_previous_results(
        previous: dict[str, Any],
        record: Record,
        _options: ConfigType) -> FungalBoundaryResults | None:
    """Restore saved results through antiSMASH's standard hook."""
    return FungalBoundaryResults.from_json(previous, record)


def run_on_record(record: Record, results: FungalBoundaryResults | None,
                  _options: ConfigType) -> FungalBoundaryResults:
    """Pass the antiSMASH record directly to Fundaries."""
    if results is not None:
        return results

    logging.info("Predicting fungal boundaries with Fundaries for %s",
                 record.id)
    report = fundaries.main(record)
    return FungalBoundaryResults(report)


def generate_html(region_layer: RegionLayer, results: FungalBoundaryResults,
                  _record_layer: RecordLayer,
                  options_layer: OptionsLayer) -> HTMLSections:
    """Render the report's predictions for this region's anchors."""
    html = HTMLSections(NAME)
    if not options_layer.fungalnew:
        return html
    anchors = {sub.label for sub in region_layer.subregions if sub.tool == NAME}
    predictions = [prediction for prediction in results.report.predictions
                   if prediction.anchor_id in anchors]
    template = FileTemplate(path.get_full_path(__file__, "templates", "fungal_boundary.html"))
    html.add_detail_section(NAME, template.render(predictions=predictions),
                            "fungal_boundary-details")
    return html


def will_handle(_products: list[str], _product_categories: set[str]) -> bool:
    """Allow motif details for any region product type."""
    return True
