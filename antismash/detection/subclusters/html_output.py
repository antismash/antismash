# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Manages HTML output for the subcluster detection module """

from typing import Any

from antismash.common import path
from antismash.common.html_renderer import FileTemplate, HTMLSections, Markup
from antismash.common.layers import RecordLayer, RegionLayer
from antismash.common.secmet import Record, Region
from antismash.config import ConfigType

from .results import SubclusterDetectionResults


def will_handle(_products: list[str], _categories: set[str]) -> bool:
    """ Relevant to every region, so return True for every product """
    return True


def generate_html(region_layer: RegionLayer, results: SubclusterDetectionResults,
                  _record_layer: RecordLayer, _options: ConfigType) -> HTMLSections:
    """ Generates the HTML detail sections for the subclusters within a region """
    predictions = results.get_predictions_for_region(region_layer.region_feature)

    tooltip = Markup("Subclusters are sets of genes responsible for producing"
                     " a specific chemical moiety.")

    template = FileTemplate(path.get_full_path(__file__, "templates", "details.html"))
    section = template.render(predictions=predictions, tooltip=tooltip,
                              anchor=region_layer.anchor_id)

    html = HTMLSections(name="subclusters")
    html.add_detail_section("Subclusters", section, class_name="subclusters")
    return html


def generate_javascript_data(_record: Record, region: Region,
                             results: SubclusterDetectionResults) -> list[Any]:
    """ Generates JSON data for the javascript to draw subclusters in HTML output

        Arguments:
            _record: the record the results belong to, unused as the region refers to it
            region: the specific region to generate data for
            results: the results of the module

        Returns:
            a list of JSON-friendly dictionaries, one for each subcluster in the region
    """
    region_anchor = RegionLayer.build_anchor_id(region)
    template = FileTemplate(path.get_full_path(__file__, "templates", "cds_tooltip.html"))

    predictions = results.get_predictions_for_region(region)

    data = []
    for i, prediction in enumerate(predictions, start=1):
        cds_results = []
        for cds_name, hits in prediction.domain_hits_by_cds.items():
            domains = []
            for hit in hits:
                domains.append({
                    "name": hit.domain_name,
                    "description": hit.domain_description,
                    "accession": hit.domain_accession,
                    "evalue": hit.evalue,
                    "bitscore": hit.bitscore,
                })
            cds_results.append({
                "cds": cds_name,
                "domains": domains,
                "tooltip": template.render(domains=domains, cds_name=cds_name),
            })

        entry: dict[str, Any] = {
            "identifier": f"subclusters-svg-{region_anchor}-sc{i}",
            "cds_results": cds_results,
        }
        entry["extends_left"], entry["extends_right"] = prediction.extends_beyond(region.location)
        data.append(entry)

    return data
