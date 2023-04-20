# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Manages HTML construction for the TFBS finder module
"""

import bisect
from typing import Any, Sequence

from Bio.Seq import Seq

from antismash.common import path
from antismash.common.html_renderer import FileTemplate, HTMLSections, Markup, docs_link
from antismash.common.layers import OptionsLayer, RegionLayer, RecordLayer
from antismash.common.secmet import CDSFeature, Feature, FeatureLocation, Record, Region
from antismash.common.secmet.locations import location_contains_other

from .tfbs_finder import Confidence, TFBSFinderResults


AMBIGUITY_BASES = {
    "K": "TG",
    "Y": "TC",
    "S": "CG",
    "W": "AT",
    "R": "AG",
    "M": "CA",
}
PRE_SEQUENCE_SIZE = 3
POST_SEQUENCE_SIZE = 3


def will_handle(_products: list[str], _categories: set[str]) -> bool:
    """ Returns true if one or more relevant products or product categories are present """
    return True


def generate_html(region_layer: RegionLayer, results: TFBSFinderResults,
                  _record_layer: RecordLayer, _options_layer: OptionsLayer
                  ) -> HTMLSections:
    """ Generates HTML output for the module """
    html = HTMLSections("tfbs-finder")

    template = FileTemplate(path.get_full_path(__file__, "templates", "sidepanel.html"))
    body = FileTemplate(path.get_full_path(__file__, "templates", "details.html"))

    all_region_hits = results.get_hits_by_region(region_layer.get_region_number())

    if not all_region_hits:
        return html
    descriptions = {hit.name: hit.description for hit in all_region_hits}

    tooltip = "Shows descriptions for Transcription Factor Binding Site models"
    section = template.render(results=results, descriptions=descriptions, tooltip=tooltip)
    html.add_sidepanel_section("TFBS Finder", section, class_name="tfbs-finder")

    weak = results.get_hits_by_region(region_layer.get_region_number(), confidence=Confidence.WEAK)
    other = results.get_hits_by_region(region_layer.get_region_number(), confidence=Confidence.MEDIUM,
                                       allow_better=True)
    assert other or weak, all_region_hits

    tooltip = Markup(
        "Detailed information for Transcription Factor Binding Site hits, including "
        "surrounding genes where applicable. <br>"
        "The strand of the match is shown by an arrow underneath the match."
        "<br>"
        f"Detailed documentation is available {docs_link('here', 'modules/tfbs')}."
    )
    details = body.render(hits=other, weak_hits=weak, results=results,
                          tooltip=tooltip, anchor=region_layer.anchor_id)
    html.add_detail_section("TFBS Finder", details, class_name="tfbs-finder")

    return html


def get_sequence_matches(query: str, consensus: str) -> list:
    """ Compares hit sequence to the consensus sequence and report
        matching nucleotides, including ambiguity codes.

        Arguments:
            query: the sequence from the query record
            consensus: the consensus sequence from the PWM

        Returns:
            a list of booleans, one for each position in the input sequences
    """
    if len(query) != len(consensus):
        raise ValueError("query and consensus sequences must be the same length")
    matches = []
    for base, cons in zip(query, consensus):
        matches.append(base == cons or cons == "N" or base in AMBIGUITY_BASES.get(cons, ""))
    return matches


def add_neighbouring_genes(hit: dict[str, Any], genes: Sequence[CDSFeature]) -> dict[str, Any]:
    """ Adds neighbouring gene information to a JSON-ready representation of a hit

        Arguments:
            hit: the hit to add information to
            genes: the CDS features within the relevant region (sorted as from secmet)

        Returns:
            the given hit instance with modifications
    """
    start = hit["start"]
    end = hit["end"]
    dummy = Feature(FeatureLocation(start, end), feature_type="dummy")
    index = bisect.bisect_left(genes, dummy)
    if index == 0:
        hit["left"] = None
    else:
        left = genes[index - 1]
        assert start > left.location.start
        hit["left"] = {
            "name": left.get_name(),
            "location": left.location.end,
            "strand": left.location.strand
        }
    contained_by_left = index > 0 and location_contains_other(left.location, dummy.location)
    hit["contained_by_left"] = contained_by_left
    if contained_by_left:
        hit["right"] = {
            "name": left.get_name(),
            "location": left.location.end,
            "strand": left.location.strand
        }
        hit["left"]["location"] = left.location.start
    else:
        if index >= len(genes):
            hit["right"] = None
        else:
            right = genes[index]
            # some gene annotations can be smaller than the hit, in which case
            # keep it as the "mid" gene and look further for the next gene
            if right.is_contained_by(dummy):
                hit["right"] = None  # in case there is no right
                hit["mid"] = {
                    "name": right.get_name(),
                    "location": right.location.start,
                    "length": right.location.end - right.location.start,
                    "strand": right.location.strand,
                }
                # if the record is poorly annotated enough to have multiple
                # extremely short genes, only show one and move on to the next
                # meaningful one
                while right.is_contained_by(dummy):
                    index += 1
                    # stop if there's no more genes
                    if index >= len(genes):
                        return hit
                    right = genes[index]
            assert end < right.location.end
            hit["right"] = {
                "name": right.get_name(),
                "location": right.location.start,
                "strand": right.location.strand
            }
    if hit["contained_by_left"]:
        assert hit["right"] and hit["left"]
    return hit


def generate_javascript_data(record: Record, region: Region, results: TFBSFinderResults
                             ) -> list[dict[str, Any]]:
    """ Generates JSON data for javascript visualisation """

    hits_in_region = results.get_hits_by_region(region.get_region_number())
    if not hits_in_region:
        return []

    converted: list[dict[str, Any]] = []
    cds_features = region.cds_children
    for hit in hits_in_region:
        end = hit.start + len(hit)
        # don't extend context before the start of the record
        prefix_size = min(PRE_SEQUENCE_SIZE, hit.start)
        # neither should the context extend past the end of the record
        suffix_size = min(POST_SEQUENCE_SIZE, len(record.seq) - end)
        core = str(record.seq[hit.start:end])
        consensus = hit.consensus
        if hit.strand == -1:
            consensus = str(Seq(consensus).reverse_complement())

        data = {
            "name": hit.name,
            "start": hit.start,
            "end": end,
            "confidence": str(hit.confidence),
            "score": hit.score,
            "presequence": str(record.seq[hit.start - prefix_size: hit.start]),
            "sequence": core,
            "postsequence": str(record.seq[end:end + suffix_size]),
            "target": consensus,
            "strand": hit.strand,
            "matches": get_sequence_matches(core, consensus),
        }
        converted.append(add_neighbouring_genes(data, cds_features))
    return converted
