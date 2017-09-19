# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

from jinja2 import FileSystemLoader, Environment, StrictUndefined

def will_handle(product):
    return True

def generate_details_div(cluster_layer, record_layer, options_layer):
    cluster = cluster_layer.cluster_rec
    divs = []
    if options_layer.cb_general or cluster.clusterblast is not None:
        divs.append(generate_div(cluster_layer, record_layer, options_layer, "clusterblast"))
    if options_layer.cb_knownclusters or cluster.knownclusterblast is not None:
        divs.append(generate_div(cluster_layer, record_layer, options_layer, "knownclusterblast"))
    if options_layer.cb_subclusters or cluster.subclusterblast is not None:
        divs.append(generate_div(cluster_layer, record_layer, options_layer, "subclusterblast"))
    return "\n".join(divs)

def generate_div(cluster_layer, record_layer, options_layer, search_type):
    env = Environment(
        loader=FileSystemLoader(['antismash/modules/clusterblast/templates']),
        autoescape=True, undefined=StrictUndefined)
    template = env.get_template('%s.html' % search_type)
    details_div = template.render(record=record_layer,
                           cluster=cluster_layer,
                           options=options_layer)
    return details_div
