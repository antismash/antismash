# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Handles HTML output for sactipeptides """

from typing import List

from jinja2 import FileSystemLoader, Environment, StrictUndefined

from antismash.common import path
from antismash.common.secmet import Prepeptide, Region
from antismash.common.layers import RegionLayer, RecordLayer, OptionsLayer

from .specific_analysis import SactiResults


def will_handle(products: List[str]) -> bool:
    """ Returns True if this module can handle the provided region products """
    return 'sactipeptide' in products


def generate_details_div(region_layer: RegionLayer, results: SactiResults,
                         record_layer: RecordLayer, options_layer: OptionsLayer) -> str:
    """ Generates the main page section of HTML with any results for the given
        region """
    env = Environment(loader=FileSystemLoader(path.get_full_path(__file__, 'templates')),
                      autoescape=True, undefined=StrictUndefined)
    template = env.get_template('details.html')
    motifs_in_region = {}
    for locus in results.regions.get(region_layer.get_region_number(), []):
        motifs_in_region[locus] = results.motifs_by_locus[locus]
    details_div = template.render(record=record_layer,
                                  region=SactipeptideLayer(record_layer, region_layer.region_feature),
                                  options=options_layer,
                                  results=motifs_in_region)
    return details_div


def generate_sidepanel(region_layer: RegionLayer, results: SactiResults,
                       record_layer: RecordLayer, options_layer: OptionsLayer) -> str:
    """ Generates the sidepanel section of HTML with any results for the given
        region """
    env = Environment(loader=FileSystemLoader(path.get_full_path(__file__, 'templates')),
                      autoescape=True, undefined=StrictUndefined)
    template = env.get_template('sidepanel.html')
    region = SactipeptideLayer(record_layer, region_layer.region_feature)
    motifs_in_region = {}
    for locus in results.regions.get(region_layer.get_region_number(), []):
        motifs_in_region[locus] = results.motifs_by_locus[locus]
    sidepanel = template.render(record=record_layer,
                                region=region,
                                options=options_layer,
                                results=motifs_in_region)
    return sidepanel


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
