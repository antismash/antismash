# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Detect fungal BGC boundaries based on shared 5' UTR sequence motifs. """

from __future__ import annotations

from typing import Any

from antismash.common.html_renderer import HTMLSections
from antismash.common.layers import OptionsLayer, RecordLayer, RegionLayer
from antismash.common.secmet import Record
from antismash.config import ConfigType
from antismash.config.args import ModuleArgs
from antismash.detection import DetectionStage

from .prediction import BoundaryDetectionResults, predict_boundaries

NAME = "fungal_boundaries"
SHORT_DESCRIPTION = "Predict fungal cluster boundaries"
DETECTION_STAGE = DetectionStage.AREA_REFINEMENT


def get_arguments() -> ModuleArgs:
    """Build the module arguments."""
    args = ModuleArgs("BGC boundaries prediction", "boundaries")
    args.add_analysis_toggle("boundaries-fungal", dest="boundaries_fungal",
                             action="store_true", default=False,
                             help="Predict fungal cluster boundaries.")
    return args


def is_enabled(options: ConfigType) -> bool:
    """Return whether the module is enabled."""
    return options.boundaries_fungal


def check_options(options: ConfigType) -> list[str]:
    """Restrict the enabled module to fungal analysis."""
    if is_enabled(options) and options.taxon != "fungi":
        return ["Fungal boundary prediction requires fungal sequences."]
    return []


def check_prereqs(_options: ConfigType) -> list[str]:
    """Check that all prerequisites are satisfied."""
    return []


def regenerate_previous_results(
        previous: dict[str, Any],
        record: Record,
        _options: ConfigType) -> BoundaryDetectionResults | None:
    """Rebuild previous results"""
    return BoundaryDetectionResults.from_json(previous, record)


def run_on_record(record: Record, results: BoundaryDetectionResults | None,
                  _options: ConfigType) -> BoundaryDetectionResults:
    """Predict fungal BGC boundaries on the provided record."""
    if results is not None:
        return results

    return predict_boundaries(record)