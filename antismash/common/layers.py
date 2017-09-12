# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

import functools
import logging
import os

class RecordLayer:
    def __init__(self, seq_record, options):
        from antismash.outputs.html.js import convert_record, load_cog_annotations # TODO break this circular dependency
        self.seq_record = seq_record
        self.options = options
        self.record = convert_record(self.seq_record, load_cog_annotations(), self.options.options) #TODO stop this from being called again
        self.clusters = []
        for cluster in self.record['clusters']:
            self.clusters.append(ClusterLayer(cluster, self, self.seq_record.get_cluster(cluster['idx'])))

    @property
    def seq_name(self):
         return self.seq_record.name

    @property
    def number_clusters(self):
        " returns the number of clusters of the record "
        return len(self.record['clusters'])

    @property
    def seq_id(self):
        return self.record['seq_id']

    @property
    def orig_id(self):
        return self.record['orig_id']

    @property
    def name(self):
      pass

    @property
    def input_type(self):
        return self.options.input_type

    @property
    def has_details(self):
        return functools.reduce((lambda x, y: x or y),
                    map(lambda cluster: cluster.has_details(), self.clusters))

    def get_from_record(self):
        " returns the text to be displayed in the overview table > separator-text "
        if self.input_type == 'nucl':
            id_text = self.seq_id

            if self.orig_id:
                if len(self.orig_id) < 40:
                    id_text += " (original name was: %s)" % self.orig_id
                else:
                    id_text += " (original name was: %s...)" % self.orig_id[:60]
            text_from = ('The following clusters are from record %s:' % id_text)
        else:
            text_from = ('The following clusters were found in your protein input:')
        return text_from

    def no_result_note(self):
        " returns the 'no result note' to be displayed in overview page/table "
        if self.options.input_type == 'nucl':
            no_result_note = 'No secondary metabolite clusters were found in the input sequence(s)'
        else:
            no_result_note = 'No secondary metabolite biosynthesis proteins were found in the input sequence(s)'
        return no_result_note


class ClusterLayer:
    def __init__(self, cluster, record, cluster_rec):
        self.record = record
        self.cluster = cluster
        self.handlers = []
        self.cluster_rec = cluster_rec
        self.cluster_blast = []
        self.knowncluster_blast = []
        self.subcluster_blast = []
        if self.record.options.cb_knownclusters:
            self.knowncluster_blast_generator()
        if self.record.options.cb_subclusters:
            self.subcluster_blast_generator()
        if self.record.options.cb_general:
            self.cluster_blast_generator()

        self.find_plugins_for_cluster()
        self.has_details = self.determine_has_details()
        self.has_sidepanel = self.determine_has_sidepanel()
        self.probability = "BROKEN" #TODO
        self.has_domain_alignment = self.determine_has_domain_alignment()

    @property
    def type(self):
        return self.cluster_rec.get_product_string()

    @property
    def subtype(self):
        " for hybrid clusters, returns the list with the subtypes "
        return self.cluster_rec.products

    @property
    def hybrid(self):
        " returns if cluster is hybrid"
        return len(self.cluster_rec.products) > 1

    @property
    def idx(self):
        return self.cluster_rec.get_cluster_number()

    @property
    def start(self):
        return int(self.cluster_rec.location.start) + 1

    @property
    def end(self):
        return int(self.cluster_rec.location.end)

    @property
    def cluster_rec_qualifiers(self):
        return self.cluster_rec.to_biopython()[0].qualifiers

    @property
    def knowncluster(self):
        return self.cluster_rec.knownclusterblast

    @property
    def BGCid(self):
        return format(self.cluster['BGCid'].split('_')[0])


    @property
    def detection_rules(self):
        return self.cluster_rec.detection_rules

    def description_text(self):
        " returns the gene cluster description "
        if self.record.options.input_type == 'nucl':
            description_text = self.record.seq_name +' - Gene Cluster %s. Type = %s. Location: %s - %s nt. ' %(self.idx, self.type, self.start, self.end)
        else:
            description_text = self.record.seq_id + '- Gene Cluster %s. Type = %s. ' %(self.idx, self.type)
        if 'probability' in self.cluster:
            description_text += 'ClusterFinder probability: %s. ' % self.probability
        description_text += 'Click on genes for more information.'

        return description_text

    def cluster_blast_generator(self): # TODO: deduplicate
        if self.record.options.cb_general:
            top_hits = self.cluster_rec.clusterblast[:self.record.options.nclusters]
        for i, label in enumerate(top_hits):
            i += 1 # 1-indexed
            svg_file = os.path.join('svg', 'clusterblast%s_%s.svg' % (self.idx, i))
            opt_text = 'Cluster %s hit %s %s' % (self.idx, i, label)
            self.cluster_blast.append((opt_text, svg_file))


    def knowncluster_blast_generator(self):
        if self.record.options.cb_knownclusters:
            top_hits = self.cluster_rec.knownclusterblast[:self.record.options.nclusters]
        for i, label_pair in enumerate(top_hits):
            i += 1 # 1-indexed
            label, bgc_id = label_pair
            svg_file = os.path.join('svg', 'knownclusterblast%s_%s.svg' % (self.idx, i))
            opt_text = 'Cluster %s hit %s %s' % (self.idx, i, label)
            self.knowncluster_blast.append((opt_text, svg_file))


    def subcluster_blast_generator(self):
        if self.record.options.cb_subclusters:
            assert self.cluster_rec.subclusterblast is not None, self.cluster_rec.location
            top_hits = self.cluster_rec.subclusterblast[:self.record.options.nclusters]
        for i, label in enumerate(top_hits):
            i += 1 # since one-indexed
            svg_file = os.path.join('svg', 'subclusterblast%s_%s.svg' % (self.idx, i))
            opt_text = 'Cluster %s hit %s %s' % (self.idx, i, label)
            self.subcluster_blast.append((opt_text, svg_file))

    def find_plugins_for_cluster(self):
        "Find a specific plugin responsible for a given gene cluster type"
        product = self.type
        for plugin in self.record.options.plugins:
            if not hasattr(plugin, 'will_handle'):
                continue
            if plugin.will_handle(product):
                self.handlers.append(plugin)
        return self.handlers

    def determine_has_details(self):
        self.has_details = False
        for handler in self.handlers:
            if "generate_details_div" in dir(handler):
                self.has_details = True
        return self.has_details

    def determine_has_sidepanel(self):
        self.has_sidepanel = False
        for handler in self.handlers:
            if "generate_sidepanel" in dir(handler):
                self.has_sidepanel = True
        return self.has_sidepanel

    def determine_has_domain_alignment(self):
        self.has_domain_alignment = False
        for handler in self.handlers:
            if "generate_domain_alignment_div" in dir(handler):
                self.has_domain_alignment = True
        return self.has_domain_alignment


class OptionsLayer():
    def __init__(self, options):
        self.options = options

    @property
    def plugins(self):
        return self.options.all_enabled_modules

    @property
    def input_type(self):
        return self.options.input_type

    @property
    def taxon(self):
        return self.options.taxon

    @property
    def nclusters(self):
        return self.options.cb_nclusters

    @property
    def cb_general(self):
        return self.options.cb_general

    @property
    def cb_knownclusters(self):
        return self.options.cb_knownclusters

    @property
    def cb_subclusters(self):
        return self.options.cb_subclusters

    @property
    def output_dir(self):
        return self.options.output_dir

    def get_options(self):
        return self.options

    @property
    def has_docking(self):
        docking = False
        if "docking" in self.options:
            docking = True
        return docking

    @property
    def smcogs(self):
        return self.options.smcogs

    @property
    def tta(self):
        return self.options.tta

    @property
    def cassis(self):
        return self.options.cassis

    @property
    def borderpredict(self):
        return self.options.borderpredict

    @property
    def limit(self):
        return self.options.limit

    @property
    def triggered_limit(self):
        logging.critical("OptionsLayer.triggered_limit will always be False")
        return False #self.options.triggered_limit

    @property
    def transatpks_da(self):
        logging.critical("OptionsLayer.transatpks_da will always be False")
        return False # self.options.transatpks_da

    @property
    def logfile(self):
        return 'logfile' in self.options

    def download_logfile(self):
        logfile_path = os.path.abspath(self.options.logfile)
        if os.path.dirname(logfile_path).startswith(self.options.output_dir):
            rel_logfile_path = logfile_path[len(self.options.output_dir)+1:]
            return rel_logfile_path
        return None

    @property
    def base_url(self):
        if self.options.taxon == "fungi":
            return self.options.urls.fungi_baseurl
        return self.options.urls.bacteria_baseurl
