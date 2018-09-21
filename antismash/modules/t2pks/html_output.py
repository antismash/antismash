# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Handles HTML output for the Type II PKS module """

from typing import List

from jinja2 import FileSystemLoader, Environment, StrictUndefined

from antismash.common import path
from antismash.common.layers import RegionLayer, RecordLayer, OptionsLayer

from .results import T2PKSResults


def will_handle(products: List[str]) -> bool:
    """ Returns true if one or more relevant products are present """
    return "t2pks" in products


def generate_sidepanel(region_layer: RegionLayer, results: T2PKSResults,
                       _record_layer: RecordLayer, _options_layer: OptionsLayer) -> str:
    """ Generate the sidepanel HTML with results from the type II PKS module """
    env = Environment(loader=FileSystemLoader(path.get_full_path(__file__, 'templates')),
                      autoescape=True, undefined=StrictUndefined)
    template = env.get_template('sidepanel.html')
    predictions = []
    for supercluster in region_layer.superclusters:
        for cluster in supercluster.clusters:
            if cluster.product == "t2pks":
                predictions.append(results.cluster_predictions[cluster.get_cluster_number()])
    return template.render(predictions=predictions)
