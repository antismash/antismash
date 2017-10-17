# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

import logging
import os


class RecordLayer:
    def __init__(self, seq_record, results, options):
        self.results = results
        self.seq_record = seq_record
        self.options = options
        self.clusters = []
        for cluster in seq_record.get_clusters():
            self.clusters.append(ClusterLayer(self, cluster))

    @property
    def seq_name(self):
        return self.seq_record.name

    @property
    def number_clusters(self):
        " returns the number of clusters of the record "
        return len(self.seq_record.get_clusters())

    @property
    def seq_id(self):
        return self.seq_record.id

    @property
    def orig_id(self):
        logging.critical("using dummy orig_id in RecordLayer")
        return "DUMMY ORIG_ID"

    def get_from_record(self):
        " returns the text to be displayed in the overview table > separator-text "
        if self.options.input_type == 'nucl':
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
    def __init__(self, record, cluster_rec):
        self.record = record
        self.handlers = []
        self.cluster_rec = cluster_rec
        self.cluster_blast = []
        self.knowncluster_blast = []
        self.subcluster_blast = []
        if self.cluster_rec.knownclusterblast:
            self.knowncluster_blast_generator()
        if self.cluster_rec.subclusterblast:
            self.subcluster_blast_generator()
        if self.cluster_rec.clusterblast:
            self.cluster_blast_generator()

        self.find_plugins_for_cluster()
        self.has_details = self.determine_has_details()
        self.has_sidepanel = self.determine_has_sidepanel()
        self.probability = "BROKEN"  # TODO when clusterfinder returns

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
    def knowncluster(self):
        return self.cluster_rec.knownclusterblast

    @property
    def BGCid(self):
        if not self.cluster_rec.knownclusterblast:
            return "-"
        return format(self.cluster_rec.knownclusterblast[0][1])

    @property
    def detection_rules(self):
        return ["%s: %s" % (product, rules) for product, rules in zip(self.cluster_rec.products, self.cluster_rec.detection_rules)]

    def description_text(self):
        " returns the gene cluster description "
        if self.record.options.input_type == 'nucl':
            description_text = self.record.seq_name \
                    + ' - Gene Cluster %s. Type = %s. Location: %s - %s nt. ' % (
                            self.idx, self.type, self.start, self.end)
        else:
            description_text = self.record.seq_id + '- Gene Cluster %s. Type = %s. ' % (self.idx, self.type)
        if self.probability != "BROKEN":  # TODO: real value check
            description_text += 'ClusterFinder probability: %s. ' % self.probability
        description_text += 'Click on genes for more information.'

        return description_text

    def cluster_blast_generator(self):  # TODO: deduplicate
        top_hits = self.cluster_rec.clusterblast[:self.record.options.cb_nclusters]
        for i, label in enumerate(top_hits):
            i += 1  # 1-indexed
            svg_file = os.path.join('svg', 'clusterblast%s_%s.svg' % (self.idx, i))
            self.cluster_blast.append((label, svg_file))

    def knowncluster_blast_generator(self):
        top_hits = self.cluster_rec.knownclusterblast[:self.record.options.cb_nclusters]
        for i, label_pair in enumerate(top_hits):
            i += 1  # 1-indexed
            label = label_pair[0]
            svg_file = os.path.join('svg', 'knownclusterblast%s_%s.svg' % (self.idx, i))
            self.knowncluster_blast.append((label, svg_file))

    def subcluster_blast_generator(self):
        assert self.cluster_rec.subclusterblast is not None, self.cluster_rec.location
        top_hits = self.cluster_rec.subclusterblast[:self.record.options.cb_nclusters]
        for i, label in enumerate(top_hits):
            i += 1  # since one-indexed
            svg_file = os.path.join('svg', 'subclusterblast%s_%s.svg' % (self.idx, i))
            self.subcluster_blast.append((label, svg_file))

    def find_plugins_for_cluster(self):
        "Find a specific plugin responsible for a given gene cluster type"
        products = self.subtype
        for plugin in self.record.options.plugins:
            if not hasattr(plugin, 'will_handle'):
                continue
            if plugin.will_handle(products):
                self.handlers.append(plugin)
        return self.handlers

    def determine_has_details(self):
        self.has_details = False
        for handler in self.handlers:
            if "generate_details_div" in dir(handler):
                self.has_details = True
                break
        return self.has_details

    def determine_has_sidepanel(self):
        self.has_sidepanel = False
        for handler in self.handlers:
            if "generate_sidepanel" in dir(handler):
                self.has_sidepanel = True
                break
        return self.has_sidepanel


class OptionsLayer:
    def __init__(self, options):
        self.options = options

    def __getattr__(self, attr):
        if attr in self.__dict__:
            return super().__getattribute__(attr)
        return getattr(self.options, attr)

    @property
    def plugins(self):
        from antismash.main import get_all_modules
        return get_all_modules()

    @property
    def smcogs(self):
        return not self.options.minimal or self.options.smcogs_enabled or self.options.smcogs_trees  # TODO work out a better way of doing this

    @property
    def triggered_limit(self):
        return self.options.triggered_limit

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
