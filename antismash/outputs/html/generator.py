# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

import json
import string
import os

import jinja2
from jinja2 import FileSystemLoader, Environment

from antismash.common import path, deprecated as utils
from antismash.common.layers import RecordLayer, OptionsLayer
from antismash.outputs.html import js


def write_geneclusters_js(records, output_dir, extra_data):
    with open(os.path.join(output_dir, 'geneclusters.js'), 'w') as result:
        geneclusters = {}
        for record in records:
            for cluster in record['clusters']:
                idx = cluster['idx']
                geneclusters['cluster-%s' % idx] = cluster
        result.write('var geneclusters = %s;\n' % json.dumps(geneclusters, indent=4))

        js_domains = {}
        for domain in extra_data['js_domains']:
            idx = domain['id'].split('-')[1]
            js_domains['cluster-%s' % idx] = domain
        result.write('var details_data = %s;\n' % json.dumps(js_domains, indent=4))


def generate_webpage(seq_records, results, options):

    generate_searchgtr_htmls(seq_records, options)

    json_records = js.convert_records(seq_records, results, options)

    extra_data = dict(js_domains=[])

    for i, record in enumerate(seq_records):
        json_record = json_records[i]
        json_record['seq_id'] = "".join(char for char in json_record['seq_id'] if char in string.printable)
        for json_cluster in json_record['clusters']:  # json clusters
            from antismash import get_analysis_modules  # TODO break circular dependency
            handlers = find_plugins_for_cluster(get_analysis_modules(), json_cluster)
            for handler in handlers:
                # if there's no results for the module, don't let it try
                if handler.__name__ not in results[i]:
                    continue
                if "generate_js_domains" in dir(handler):
                    domains = handler.generate_js_domains(json_cluster, record, results[i][handler.__name__], options)
                    if domains:
                        extra_data['js_domains'].append(domains)

    write_geneclusters_js(json_records, options.output_dir, extra_data)

    # jinja
    with open(os.path.join(options.output_dir, 'index.html'), 'w') as result:
        env = Environment(autoescape=True, trim_blocks=True, lstrip_blocks=True,
                          undefined=jinja2.StrictUndefined,
                          loader=FileSystemLoader(path.get_full_path(__file__)))
        template = env.get_template('index.html')
        options_layered = OptionsLayer(options)
        records = [RecordLayer(record, result, options_layered) for record, result in zip(seq_records, results)]

        records_written = sum([len(record.seq_record.get_clusters()) for record in records])
        aux = template.render(records=records, options=options_layered,
                              utils=utils, extra_data=extra_data,
                              records_written=records_written,
                              config=options)
        result.write(aux)


def find_plugins_for_cluster(plugins, cluster):
    "Find a specific plugin responsible for a given gene cluster type"
    products = cluster['products']
    handlers = []
    for plugin in plugins:
        if not hasattr(plugin, 'will_handle'):
            continue
        if plugin.will_handle(products):
            handlers.append(plugin)
    return handlers


def load_searchgtr_search_form_template():
    """Create folder for SEARCHGTR HTML files, load search form template """
    template = open(path.get_full_path(__file__, "searchgtr_form.html"), "r")
    template = template.read()
    template = template.replace("\r", "\n")
    template_parts = template.split("FASTASEQUENCE")
    return template_parts


def generate_searchgtr_htmls(seq_records, options):
    """ Generate lists of COGs that are glycosyltransferases or transporters """
    gtrcoglist = ['SMCOG1045', 'SMCOG1062', 'SMCOG1102']
    searchgtrformtemplateparts = load_searchgtr_search_form_template()
    # TODO store somewhere sane
    js.searchgtr_links = {}
    for seq_record in seq_records:
        for feature in seq_record.get_cds_features():
            smcog_functions = feature.gene_functions.get_by_tool("smcogs")
            if not smcog_functions:
                continue
            smcog = smcog_functions[0].description.split(":")[0]
            if smcog not in gtrcoglist:
                continue
            html_dir = os.path.join(options.output_dir, "html")
            if not os.path.exists(html_dir):
                os.mkdir(html_dir)
            formfileloc = os.path.join(html_dir, feature.get_name() + "_searchgtr.html")
            link_loc = os.path.join("html", feature.get_name() + "_searchgtr.html")
            gene_id = feature.get_name()
            js.searchgtr_links[seq_record.id + "_" + gene_id] = link_loc
            with open(formfileloc, "w") as formfile:
                specificformtemplate = searchgtrformtemplateparts[0].replace("GlycTr", gene_id)
                formfile.write(specificformtemplate)
                formfile.write("%s\n%s" % (gene_id, feature.translation))
                formfile.write(searchgtrformtemplateparts[1])
