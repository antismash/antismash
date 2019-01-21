# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Handles HTML output for sactipeptides """

from collections import defaultdict
from typing import List
from typing import Dict  # comment hint, pylint: disable=unused-import

from antismash.common import path
from antismash.common.html_renderer import HTMLSections, FileTemplate
from antismash.common.secmet import Prepeptide, Region
from antismash.common.secmet import CDSMotif # comment hint, pylint: disable=unused-import
from antismash.common.layers import RegionLayer, RecordLayer, OptionsLayer

from .specific_analysis import SactiResults


def will_handle(products: List[str]) -> bool:
    """ Returns True if this module can handle the provided region products """
    return 'sactipeptide' in products


def generate_html(region_layer: RegionLayer, results: SactiResults,
                  record_layer: RecordLayer, options_layer: OptionsLayer) -> HTMLSections:
    """ Generates HTML for the module """
    html = HTMLSections("sactipeptides")

    motifs_in_region = defaultdict(list)  # type: Dict[str, List[CDSMotif]]
    for locus, motifs in results.motifs_by_locus.items():
        for motif in motifs:
            if motif.is_contained_by(region_layer.region_feature):
                motifs_in_region[locus].append(motif)

    sacti_layer = SactipeptideLayer(record_layer, region_layer.region_feature)

    template = FileTemplate(path.get_full_path(__file__, "templates", "details.html"))
    details = template.render(record=record_layer,
                              region=sacti_layer,
                              options=options_layer,
                              results=motifs_in_region)
    html.add_detail_section("Sactipeptides", details)

    template = FileTemplate(path.get_full_path(__file__, "templates", "sidepanel.html"))
    sidepanel = template.render(record=record_layer,
                                region=sacti_layer,
                                options=options_layer,
                                results=motifs_in_region)
    html.add_sidepanel_section("Sactipeptides", sidepanel)

    return html


class SactipeptideLayer(RegionLayer):
    """ An extended RegionLayer for holding a list of LanthipeptideMotifs """
    def __init__(self, record: RecordLayer, region_feature: Region) -> None:
        RegionLayer.__init__(self, record, region_feature)
        self.motifs = []  # type: List[Prepeptide]
        for motif in self.record.seq_record.get_cds_motifs():
            if not isinstance(motif, Prepeptide):
                continue
            if not motif.is_contained_by(self.region_feature):
                continue
            if motif.peptide_class == "sactipeptide":
                self.motifs.append(motif)
