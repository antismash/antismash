"""Adapt antiSMASH inputs and results for fungal boundary analysis."""

from __future__ import annotations

from copy import deepcopy
import logging
from typing import Any

from antismash.common import path
from antismash.common.html_renderer import FileTemplate, HTMLSections
from antismash.common.layers import OptionsLayer, RecordLayer, RegionLayer
from antismash.common.module_results import DetectionResults
from antismash.common.secmet import FeatureLocation, GeneFunction, Record, SubRegion
from antismash.common.serialiser import feature_from_json, feature_to_json
from antismash.config import ConfigType
from antismash.config.args import ModuleArgs
from antismash.detection import DetectionStage

from . import fundaries

NAME = "fungal_boundary"
SHORT_DESCRIPTION = "Predict fungal cluster boundaries using Fundaries"
DETECTION_STAGE = DetectionStage.AREA_REFINEMENT


class FungalBoundaryResults(DetectionResults):
    """Store boundary bars and motif evidence for one contig.

    ``motifs[anchor]`` contains motif, repeat_type, p_empirical, and genes.
    Each gene is a dictionary with locus and total_motif_hits, including zeros.
    """

    schema_version = 2

    def __init__(self, record_id: str) -> None:
        """Create empty prediction collections for one contig."""
        super().__init__(record_id)
        self.subregions: list[SubRegion] = []
        self.motifs: dict[str, dict[str, Any]] = {}

    def to_json(self) -> dict[str, Any]:
        """Serialize boundaries and motif evidence for antiSMASH."""
        return {
            "schema_version": self.schema_version,
            "record_id": self.record_id,
            "subregions": [feature_to_json(sub.to_biopython()[0])
                           for sub in self.subregions],
            "motifs": deepcopy(self.motifs),
        }

    @staticmethod
    def from_json(json: dict[str, Any],
                  record: Record) -> FungalBoundaryResults | None:
        """Restore saved results with the same schema and contig identifier."""
        if json.get("schema_version") != FungalBoundaryResults.schema_version:
            return None
        if json.get("record_id") != record.id:
            return None
        results = FungalBoundaryResults(record.id)
        results.subregions = [
            SubRegion.from_biopython(feature_from_json(sub))
            for sub in json["subregions"]
        ]
        results.motifs = deepcopy(json["motifs"])
        return results

    def get_subregion_groups(self) -> list[list[SubRegion]]:
        """Group identical boundaries and motifs while
        retaining each anchor."""
        groups: dict[tuple[str, ...], list[SubRegion]] = {}
        for subregion in sorted(self.subregions, key=lambda sub: (
                sub.location.start, sub.location.end, sub.label)):
            summary = self.motifs.get(subregion.label, {})
            motif = summary.get("motif", "")
            key = (str(subregion.location), subregion.tool,
                   summary.get("repeat_type", ""),
                   motif, "" if motif else subregion.label)
            groups.setdefault(key, []).append(subregion)
        return list(groups.values())

    def get_predicted_subregions(self) -> list[SubRegion]:
        """Supply one boundary bar per display group to antiSMASH."""
        return [SubRegion(group[0].location, tool=NAME,
                          label=", ".join(sub.label for sub in group))
                for group in self.get_subregion_groups()]


def prepare_input(record: Record, options: ConfigType) -> dict[str, Any]:
    """Collect one contig's inputs for inspection.

    Keep antiSMASH's sequence and feature objects in memory. Gene features and
    CDS features remain separate; mapping core CDSs to genes belongs to the
    prediction implementation. The options include antiSMASH's CPU budget.
    """
    cds_features = record.get_cds_features()
    anchors = [cds.get_name() for cds in cds_features
               if cds.gene_function == GeneFunction.CORE]
    return {
        "record_id": record.id,
        "sequence": record.seq,
        "genes": record.get_genes(),
        "cds_features": cds_features,
        "anchor_names": list(dict.fromkeys(anchors)),
        "options": options,
    }


def build_results(record_id: str,
                  predictions: list[dict[str, Any]]) -> FungalBoundaryResults:
    """Convert normalized predictor output into antiSMASH results.

    Each prediction supplies anchor, start, end, motif, repeat_type,
    p_empirical, and genes. Coordinates describe the expanded display interval
    and are zero-based and end-exclusive. Gene counts are uncapped totals.
    """
    results = FungalBoundaryResults(record_id)
    for prediction in predictions:
        anchor = prediction["anchor"]
        location = FeatureLocation(prediction["start"], prediction["end"])
        results.subregions.append(SubRegion(location, tool=NAME, label=anchor))
        results.motifs[anchor] = {
            "motif": prediction["motif"],
            "repeat_type": prediction["repeat_type"],
            "p_empirical": prediction.get("p_empirical"),
            "genes": deepcopy(prediction["genes"]),
        }
    return results


def results_from_report(
        record_id: str,
        report: fundaries.FundariesReport) -> FungalBoundaryResults:
    """Convert accepted Fundaries predictions, including their
    expanded gene tables."""
    predictions = []
    for prediction in report.results:
        if prediction["status"] != "accepted":
            continue
        genes = prediction["gene_table"]
        predictions.append({
            "anchor": prediction["core_gene"],
            "start": min(gene.start for gene in genes),
            "end": max(gene.end for gene in genes),
            "motif": prediction["motif"],
            "repeat_type": prediction["type"],
            "p_empirical": prediction["p_empirical"],
            "genes": [{"locus": gene.name, "total_motif_hits": gene.motif_hits}
                      for gene in genes],
        })
    return build_results(record_id, predictions)


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
    """Predict fungal boundaries for the record's core CDS anchors."""
    if results is not None:
        return results
    if not any(cds.gene_function == GeneFunction.CORE
               for cds in record.get_cds_features()):
        return build_results(record.id, [])
    logging.info("Predicting fungal boundaries with Fundaries for %s",
                 record.id)
    report = fundaries.main(record)
    logging.info("Fundaries analysed %d core genes "
                 "and accepted %d predictions",
                 report.analysed, report.accepted)
    return results_from_report(record.id, report)


def generate_html(region_layer: RegionLayer, results: FungalBoundaryResults,
                  _record_layer: RecordLayer,
                  options_layer: OptionsLayer) -> HTMLSections:
    """Render motif evidence using the existing fungal boundary template."""
    html = HTMLSections(NAME)
    if not options_layer.fungalnew:
        return html
    locations = [sub.location for sub in region_layer.subregions
                 if sub.tool == NAME]
    predictions = []
    for group in results.get_subregion_groups():
        if group[0].location not in locations:
            continue
        anchors = [(sub.label, results.motifs.get(sub.label)) for sub in group]
        predictions.append((anchors, results.motifs.get(group[0].label)))
    tooltip = (
        "Predictions with identical boundaries, motif types, and motifs share "
        "one bar and table. Empirical p-values are shown for each anchor. "
        "Gene tables include uncapped promoter counts and zero-count genes. "
        "N followed by a number denotes the number of spacer bases "
        "between repeat arms."
    )
    template = FileTemplate(path.get_full_path(__file__, "templates",
                                               "fungal_boundary.html"))
    html.add_detail_section(NAME, template.render(predictions=predictions,
                                                  tooltip=tooltip),
                            "fungal_boundary-details")
    return html


def will_handle(_products: list[str], _product_categories: set[str]) -> bool:
    """Allow motif details for any region product type."""
    return True
