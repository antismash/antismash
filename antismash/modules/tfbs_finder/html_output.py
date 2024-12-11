# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Manages HTML construction for the TFBS finder module
"""

import bisect
from typing import Any, Optional, Sequence

from Bio.Seq import Seq

from antismash.common import path
from antismash.common.html_renderer import FileTemplate, HTMLSections, Markup, docs_link
from antismash.common.layers import OptionsLayer, RegionLayer, RecordLayer
from antismash.common.secmet import CDSFeature, Feature, Record, Region
from antismash.common.secmet.locations import (
    CompoundLocation,
    FeatureLocation,
    location_contains_other,
)

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

    weak = results.get_hits_by_region(region_layer.get_region_number(), confidence=Confidence.WEAK)
    other = results.get_hits_by_region(region_layer.get_region_number(), confidence=Confidence.MEDIUM,
                                       allow_better=True)
    assert other or weak, all_region_hits

    tooltip = "Shows descriptions for Transcription Factor Binding Site models"
    glossary_order = sorted(other, key=lambda x: x.name)
    descriptions = {hit.name: {"description": hit.description, "link": hit.link} for hit in glossary_order}
    section = template.render(results=results, descriptions=descriptions, tooltip=tooltip)
    html.add_sidepanel_section("TFBS Finder", section, class_name="tfbs-finder")

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


def _find_in_contiguous_area(start: int, end: int, genes: Sequence[CDSFeature],
                             ) -> tuple[Optional[CDSFeature], Optional[CDSFeature], Optional[CDSFeature]]:
    """ Finds the CDS features immediately next to, within, or containing the site with the
        given coordinates, if they exist.

        Arguments:
            start: the start coordinate of the hit
            end: the end coordinate of the hit
            region: the Region feature containing the site

        Returns:
            a tuple of
                the CDS to the left (which may fully or partially contain the site) or None
                the last CDS fully contained by the site, if any exist at all
                the CDS to the right (which may fully or partially contain the site) or None
    """
    dummy = Feature(FeatureLocation(start, end), feature_type="dummy")
    index = bisect.bisect_left(genes, dummy)
    left = None
    mid = None
    right = None
    if index != 0:
        left = genes[index - 1]
    while index < len(genes) and genes[index].is_contained_by(dummy):
        mid = genes[index]
        right = genes[index]
        index += 1
    if index < len(genes):
        right = genes[index]
    return left, mid, right


def find_neighbours(start: int, end: int, region: Region,
                    ) -> tuple[Optional[CDSFeature], Optional[CDSFeature], Optional[CDSFeature]]:
    """ Finds the CDS features immediately next to, within, or containing the site with the
        given coordinates, if they exist. If the region crosses the origin, the neighbouring
        features may include those over the origin from the site.

        Arguments:
            start: the start coordinate of the hit
            end: the end coordinate of the hit
            region: the Region feature containing the site

        Returns:
            a tuple of
                the CDS to the left (which may fully or partially contain the site) or None
                the last CDS fully contained by the site, if any exist at all
                the CDS to the right (which may fully or partially contain the site) or None
    """
    if not region.crosses_origin():
        return _find_in_contiguous_area(start, end, region.cds_children)
    pre_origin = region.cds_children.pre_origin
    cross_origin = region.cds_children.cross_origin
    post_origin = region.cds_children.post_origin

    left = None
    mid = None
    right = None
    # start with cross-origin sites, which are simpler
    if start > end:
        # a dummy cross-origin location, the inner dimensions don't matter here
        location = CompoundLocation([FeatureLocation(start, start + 1), FeatureLocation(0, end)])
        # starting with the cross_origin, since those may contain the site,
        # but they can't be outer limits
        if cross_origin:
            for i, cds in enumerate(cross_origin):
                if location_contains_other(cds.location, location):
                    left = cds
                    if i < len(cross_origin) - 1:
                        right = cross_origin[i]
                    elif post_origin:
                        right = post_origin[0]
                elif cds.is_contained_by(location):
                    mid = cds
            return left, mid, right
        # if between genes on either side, that's easier
        # then wrapping to pre-origin for left side
        if pre_origin:
            left = pre_origin[-1]
        # or wrapping to post-origin for right side
        if post_origin:
            right = post_origin[0]
        return left, mid, right

    # the site itself doesn't cross the origin, but the region might
    location = FeatureLocation(start, end)
    # so use the normal linear/non-crossing function for each half, as initial values
    left, mid, right = ((a or b) for a, b in zip(_find_in_contiguous_area(start, end, pre_origin),
                                                 _find_in_contiguous_area(start, end, post_origin)))

    # then look over the origin to fill in a missing value,
    # depending on which side of the origin the site is closest to
    max_coord = region.location.end
    if not left and pre_origin and start < max_coord // 2:
        left = pre_origin[-1]

    if not right and post_origin and start > max_coord // 2:
        right = post_origin[0]

    return left, mid, right


def add_neighbouring_genes(hit: dict[str, Any], left: Optional[CDSFeature], mid: Optional[CDSFeature],
                           right: Optional[CDSFeature],
                           ) -> dict[str, Any]:
    """ Adds neighbouring gene information to a JSON-ready representation of a hit

        Arguments:
            hit: the hit to add information to
            genes: the CDS features within the relevant region (sorted as from secmet)

        Returns:
            the given hit instance with modifications
    """
    start = hit["start"]
    end = hit["end"]
    hit["left"] = None
    hit["right"] = None

    if start > end:  # cross-origin
        location = CompoundLocation([FeatureLocation(start, start + 3), FeatureLocation(0, 3)])
    else:
        location = FeatureLocation(start, end)
    if left:
        hit["left"] = {
            "name": left.get_name(),
            "location": int(left.end),
            "strand": left.location.strand,
        }
        if location_contains_other(left.location, location):
            hit["contained_by_left"] = True

    if mid:
        length = mid.end - mid.start
        if mid.crosses_origin():
            length = mid.start + mid.end
        hit["mid"] = {
            "name": mid.get_name(),
            "location": int(mid.start),
            "length": length,
            "strand": mid.location.strand
        }

    if right:
        hit["right"] = {
            "name": right.get_name(),
            "location": int(right.start),
            "strand": right.location.strand
        }

    return hit


def generate_javascript_data(record: Record, region: Region, results: TFBSFinderResults
                             ) -> list[dict[str, Any]]:
    """ Generates JSON data for javascript visualisation """

    hits_in_region = results.get_hits_by_region(region.get_region_number())
    if not hits_in_region:
        return []

    converted: list[dict[str, Any]] = []
    for hit in hits_in_region:
        start = hit.start % len(record)
        end = start + len(hit)
        core = str(record.seq[start:end])
        if record.is_circular() and end > len(record):
            end %= len(record)
            core = str(record.seq[:end] + record.seq[start:])
        assert len(record) >= end
        # don't extend context before the start of the record
        prefix_size = min(PRE_SEQUENCE_SIZE, hit.start)
        # neither should the context extend past the end of the record
        suffix_size = min(POST_SEQUENCE_SIZE, len(record.seq) - end)
        consensus = hit.consensus
        if hit.strand == -1:
            consensus = str(Seq(consensus).reverse_complement())

        data = {
            "name": hit.name,
            "start": hit.start,
            "end": end,
            "confidence": str(hit.confidence),
            "species": str(hit.species),
            "link": str(hit.link),
            "score": hit.score,
            "presequence": str(record.seq[start - prefix_size: start]),
            "sequence": core,
            "postsequence": str(record.seq[end:end + suffix_size]),
            "target": consensus,
            "strand": hit.strand,
            "matches": get_sequence_matches(core, consensus),
        }
        left, mid, right = find_neighbours(hit.start, end, region)
        converted.append(add_neighbouring_genes(data, left, mid, right))
    return converted
