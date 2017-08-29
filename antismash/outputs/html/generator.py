# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

import json
import logging
import os

import jinja2
from jinja2 import FileSystemLoader, Environment

import antismash.common.deprecated as utils
import antismash.common.path as path
from antismash.common.layers import RecordLayer, OptionsLayer
from antismash.outputs.html import js

def write_geneclusters_js(records, output_dir, extra_data):
    with open(os.path.join(output_dir, 'geneclusters.js'), 'w') as h:
        geneclusters = {}
        for record in records:
            for cluster in record['clusters']:
                idx = cluster['idx']
                geneclusters['cluster-%s' % idx] = cluster
        h.write('var geneclusters = %s;\n' %
                json.dumps(geneclusters, indent=4))

        js_domains = {}
        for domain in extra_data['js_domains']:
            idx = domain['id'].split('-')[1]
            js_domains['cluster-%s' % idx] = domain
        h.write('var details_data = %s;\n' %
                json.dumps(js_domains, indent=4))


def generate_webpage(seq_records, results, options):

    generate_searchgtr_htmls(seq_records, options)

    records = js.convert_records(seq_records, results, options)

    extra_data = dict(js_domains=[], clusterblast_clusters=[],
                      subclusterblast_clusters=[], knownclusterblast_clusters=[])


    for i, record in enumerate(records):
        record['seq_id'] = utils.ascii_string(record['seq_id'])
        for cluster in record['clusters']:
            from antismash import gather_modules #TODO break circular dependency
            handlers = find_plugins_for_cluster(gather_modules(), cluster)
            for handler in handlers:
                if "generate_js_domains" in dir(handler):
                    handler.generate_js_domains(cluster, seq_records[i], options,
                                                          extra_data['js_domains'])


    write_geneclusters_js(records, options.output_dir, extra_data)

    # jinja
    with open(os.path.join(options.output_dir, 'index.html'), 'w') as t:
        env = Environment(autoescape=True, trim_blocks=True, lstrip_blocks=True,
                undefined=jinja2.StrictUndefined,
                loader=FileSystemLoader(['antismash/outputs/html']))
        template = env.get_template('index.html')
        options_layered = OptionsLayer(options)
        records = [RecordLayer(record, options_layered) for record in seq_records]

        records_written = sum([len(record.clusters) for record in records])
        aux = template.render(records=records, options=options_layered,
                              utils=utils, extra_data=extra_data,
                              records_written=records_written,
                              config=options)
        t.write(aux)

def find_plugins_for_cluster(plugins, cluster):
    "Find a specific plugin responsible for a given gene cluster type"
    product = cluster['type']
    handlers = []
    for plugin in plugins:
        if not hasattr(plugin, 'will_handle'):
            continue
        if plugin.will_handle(product):
            handlers.append(plugin)
    return handlers


def get_detection_rules(cluster_rec):
    rules = []
    for note in cluster_rec.qualifiers['note']:
        if note.startswith("Detection rule(s)"):
            rules.extend([rule.strip().replace('&', '&amp;') for rule in note[41:].split(';')])

    return rules

def load_searchgtr_search_form_template():
    #Create folder for SEARCHGTR HTML files, load search form template
    searchgtrformtemplate = open(path.get_full_path(__file__, "searchgtr_form.html"), "r")
    searchgtrformtemplate = searchgtrformtemplate.read()
    searchgtrformtemplate = searchgtrformtemplate.replace("\r", "\n")
    searchgtrformtemplateparts = searchgtrformtemplate.split("FASTASEQUENCE")
    return searchgtrformtemplateparts

def generate_searchgtr_htmls(seq_records, options):
    #Generate lists of COGs that are glycosyltransferases or transporters
    gtrcoglist = ['SMCOG1045', 'SMCOG1062', 'SMCOG1102']
    searchgtrformtemplateparts = load_searchgtr_search_form_template()
    #TODO store somewhere sane
    js.searchgtr_links = {}
    for seq_record in seq_records:
        smcogdict, _ = utils.get_smcog_annotations(seq_record)
        for feature in seq_record.get_cds_features():
            gene_id = feature.get_name()
            if gene_id not in smcogdict:
                continue
            smcog = smcogdict[gene_id]
            if smcog not in gtrcoglist:
                continue
            if not os.path.exists(options.full_outputfolder_path + os.sep + "html"):
                os.mkdir(options.full_outputfolder_path + os.sep + "html")
            formfileloc = options.full_outputfolder_path + os.sep + "html" + os.sep + feature.get_name() + "_searchgtr.html"
            link_loc = "html" + os.sep + feature.get_name() + "_searchgtr.html"
            js.searchgtr_links[seq_record.id + "_" + gene_id] = link_loc
            with open(formfileloc, "w") as formfile:
                specificformtemplate = searchgtrformtemplateparts[0].replace("GlycTr", gene_id)
                formfile.write(specificformtemplate)
                formfile.write("%s\n%s" % (gene_id, utils.get_aa_sequence(feature)))
                formfile.write(searchgtrformtemplateparts[1])
