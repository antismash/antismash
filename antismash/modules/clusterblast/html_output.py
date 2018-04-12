# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Handles HTML generation for the clusterblast variants """

from typing import List

from jinja2 import FileSystemLoader, Environment, StrictUndefined

from antismash.common import path
from antismash.common.layers import ClusterLayer, RecordLayer, OptionsLayer

from .results import ClusterBlastResults


def will_handle(_products: List[str]) -> bool:
    """ Clusterblast is relevant to every cluster, so return True for every
        product """
    return True


def generate_details_div(cluster_layer: ClusterLayer, results: ClusterBlastResults,
                         record_layer: RecordLayer, options_layer: OptionsLayer) -> str:
    """ Generates the HTML sections of the body details for all variants
        of clusterblast
    """
    cluster = cluster_layer.cluster_feature
    divs = []
    if options_layer.cb_general or cluster.clusterblast is not None:
        divs.append(generate_div(cluster_layer, results, record_layer, options_layer, "clusterblast"))
    if options_layer.cb_knownclusters or cluster.knownclusterblast is not None:
        divs.append(generate_div(cluster_layer, results, record_layer, options_layer, "knownclusterblast"))
    if options_layer.cb_subclusters or cluster.subclusterblast is not None:
        divs.append(generate_div(cluster_layer, results, record_layer, options_layer, "subclusterblast"))
    return "\n".join(divs)


def generate_div(cluster_layer: ClusterLayer, _results: ClusterBlastResults,
                 record_layer: RecordLayer, options_layer: OptionsLayer, search_type: str) -> str:
    """ Generates the specific HTML section of the body for a given variant of
        clusterblast
    """
    template_path = path.get_full_path(__file__, "templates")
    env = Environment(loader=FileSystemLoader(template_path), autoescape=True,
                      undefined=StrictUndefined)
    template = env.get_template('%s.html' % search_type)
    details_div = template.render(record=record_layer,
                                  cluster=cluster_layer,
                                  options=options_layer)
    return details_div
