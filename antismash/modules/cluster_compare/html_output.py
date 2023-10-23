# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Handles HTML generation for the clusterblast variants """

from typing import Any, Dict, List, Optional, Set, Tuple

from antismash.common import path
from antismash.common.html_renderer import FileTemplate, HTMLSections, Markup
from antismash.common.layers import OptionsLayer, RegionLayer, RecordLayer
from antismash.common.secmet import locations, Record, Region

from .data_structures import ReferenceRegion
from .results import ClusterCompareResults, ScoresByProtocluster


DISPLAY_LIMIT = 10


def will_handle(_products: List[str], _product_categories: Set[str]) -> bool:
    """ Relevant to every region, so return True for every product """
    return True


def generate_html(region_layer: RegionLayer, results: ClusterCompareResults,
                  record_layer: RecordLayer, _options_layer: OptionsLayer
                  ) -> HTMLSections:
    """ Generates the HTML sections for all variants
    """

    html = HTMLSections("cluster-compare")
    base_tooltip = ("Shows areas that are similar to the current region to a reference database.<br>"
                    "Mouseover a score cell in the table to get a breakdown of how "
                    "the score was calculated.")

    for label, db_results in results.by_database.items():
        tooltip = base_tooltip
        if db_results.description:
            tooltip += f"{db_results.description}<br>"
        if db_results.url:
            tooltip += f"<br>Click on an accession to open that entry in the {db_results.name} database."
        variant_results = db_results.by_region.get(region_layer.get_region_number(), {})
        divs: List[Tuple[str, str, Markup]] = []

        # for the MIBiG DB specifically, always show region to region by default
        if label == "MIBiG":
            pairs = sorted(variant_results.items(), key=lambda name: (1 if name == "RegionToRegion" else 0, name))
        else:
            pairs = sorted(variant_results.items())

        for variant, result in pairs:
            scores = result.scores_by_region[:DISPLAY_LIMIT]
            if not scores:
                continue
            scores_by_proto = result.details.details
            tag = variant.replace(" ", "-")
            search_type = "row"
            kind = "Protocluster to Region"
            if "ProtoToProto" in variant:
                kind = "Protocluster to Protocluster"
                search_type = "matrix"
            elif "RegionToRegion" in variant:
                kind = "Region to Region"
                search_type = "single"
                assert isinstance(scores_by_proto, list)
                scores_by_proto = scores_by_proto[:DISPLAY_LIMIT]
            div = generate_div(tag, region_layer, record_layer, search_type,
                               tooltip, scores, scores_by_proto,
                               len(divs) == 0, label, db_results.url)
            divs.append((tag, kind, div))
        template = FileTemplate(path.get_full_path(__file__, "templates", "gathered.html"))
        markup = template.render(variants=divs, class_name=label, description="Similar gene clusters",
                                 tooltip=tooltip, anchor=region_layer.anchor_id)
        html.add_detail_section(f"{label} comparison", markup, label+"-cluster-compare")

    return html


def generate_div(tag: str, region_layer: RegionLayer, record_layer: RecordLayer,
                 template_name: str, tooltip: str, results: List[Tuple[ReferenceRegion, float]],
                 proto_results: ScoresByProtocluster, active: bool,
                 label: str, url: Optional[str]) -> Markup:
    """ Generates a single div within the details body for different kinds of
        comparisons against a particular database.

        Arguments:
            tag: the type of analysis
            region_layer: the relevant RegionLayer
            record_layer: the relevant RecordLayer
            template_name: the name of the template file to use
            tooltip: the tooltip to use for the section
            results: a ranked list of ReferenceRegion and scores
            proto_results: details usd in calculating the region ranking
            active: whether this particular div should be the default drawn
            label: the name of the database to use within the div
            url: a optional URL to use for linking a reference region externally

        Returns:
            a Markup instance of the generated div
    """
    template = FileTemplate(path.get_full_path(__file__, "templates", f"{template_name}.html"))
    return template.render(tag=tag, record=record_layer, region=region_layer,
                           tooltip=tooltip, results=results, proto_results=proto_results,
                           extra_class="comparison-container-active" if active else "",
                           class_name=label, url=url)


def generate_javascript_data(_record: Record, region: Region, results: ClusterCompareResults) -> Dict[str, Any]:
    """ Generates JSON data for the javascript to draw relevant results in HTML output

        Arguments:
            record: the relevant Record for the results
            region: the specific Region to generate data for
            results: the ClusterCompareResults that need data extracted

        Returns:
            a JSON-friendly dictionary with the relevant data
    """
    data: Dict[str, Any] = {
    }
    for label, db_results in results.by_database.items():
        data[label] = {}
        variant_results = db_results.by_region.get(region.get_region_number(), {})
        for variant, result in sorted(variant_results.items()):
            scores = sorted(result.scores_by_region, key=lambda x: x[1], reverse=True)[:DISPLAY_LIMIT]
            if not scores:
                continue

            variant_data: Dict[str, Dict[str, Any]] = {
                "reference_clusters": {}
            }
            data[label][variant] = variant_data

            for reference, _ in scores:
                ref_entry: Dict[str, Any] = {
                    "start": reference.start,
                    "end": reference.end,
                    "links": [],  # added to afterwards
                    "reverse": False,  # potentially changed later
                }
                genes = {}
                for cds in reference.cdses.values():
                    gene_json = cds.get_minimal_json()
                    gene_json["linked"] = {}
                    genes[cds.name] = gene_json
                variant_data["reference_clusters"][reference.get_identifier()] = ref_entry

                mismatching_strands = 0
                for ref_cds_id, hit in result.hits_by_region.get(reference, {}).items():
                    assert locations.locations_overlap(hit.cds.location, region.location)
                    query_cds = hit.cds
                    query_point = query_cds.location.start + (query_cds.location.end - query_cds.location.start) // 2
                    ref_cds = reference.cdses[ref_cds_id]
                    subject_point = ref_cds.location.start + (ref_cds.location.end - ref_cds.location.start) // 2
                    if query_cds.location.strand != ref_cds.location.strand:
                        mismatching_strands += 1
                    genes[ref_cds.name]["linked"][region.get_region_number()] = query_cds.get_name()
                    ref_entry["links"].append({
                        "query": query_cds.get_name(),
                        "subject": ref_cds.name,
                        "query_loc": query_point,
                        "subject_loc": subject_point,
                    })
                ref_entry["reverse"] = mismatching_strands > len(ref_entry["links"]) / 2
                ref_entry["genes"] = list(genes.values())
    return data
