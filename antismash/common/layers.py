# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" A collection of classes for simplification of HTML templating.

    Each layer type is a different construct with helper functions that don't
    really belong in the original objects.
"""

import os
from typing import Any, List, Optional, Tuple

from antismash.config import ConfigType, get_config
from antismash.common.module_results import ModuleResults
from antismash.common.secmet import Record, Region
from antismash.custom_typing import AntismashModule


class OptionsLayer:
    """ A layer for the global Config options. """
    def __init__(self, options: ConfigType) -> None:
        self.options = options

    def __getattr__(self, attr: str) -> Any:
        if attr in self.__dict__:
            return super().__getattribute__(attr)
        return getattr(self.options, attr)

    @property
    def plugins(self) -> List[AntismashModule]:
        """ A list of all modules """
        from antismash.main import get_all_modules
        return get_all_modules()

    def download_logfile(self) -> Optional[str]:
        """ Returns the path of the logfile, if it was created (otherwise None) """
        logfile_path = os.path.abspath(self.logfile)
        if os.path.dirname(logfile_path).startswith(self.output_dir):
            return logfile_path[len(self.output_dir) + 1:]
        return None

    @property
    def base_url(self) -> str:
        """ Returns the 'home' URL for fungismash/antismash """
        if self.options.taxon == "fungi":
            return self.options.urls.fungi_baseurl
        return self.options.urls.bacteria_baseurl


class RecordLayer:
    """ A layer for Record instances """
    def __init__(self, record: Record, results: Optional[ModuleResults], options: OptionsLayer) -> None:
        self.results = results
        self.seq_record = record
        self.options = options
        self.regions = []  # type: List[RegionLayer]
        for region in record.get_regions():
            self.regions.append(RegionLayer(self, region))

    def __getattr__(self, attr: str) -> Any:
        if attr in self.__dict__:
            return super().__getattribute__(attr)
        return getattr(self.seq_record, attr)

    def get_from_record(self) -> str:
        " returns the text to be displayed in the overview table > separator-text "

        current_id = self.seq_record.id

        orig_id = self.seq_record.original_id

        if orig_id:
            if len(orig_id) < 40:
                current_id += " (original name was: %s)" % orig_id
            else:
                current_id += " (original name was: %s...)" % orig_id[:60]

        source = self.annotations.get("source", "")
        if source:
            source = " (%s)" % source

        return 'The following regions are from record %s%s:' % (current_id, source)


class RegionLayer:
    """ A layer for Region instances, contains special members for result of
        the clusterblast and clusterfinder modules
    """
    def __init__(self, record: RecordLayer, region_feature: Region) -> None:
        assert isinstance(region_feature, Region), type(region_feature)
        assert region_feature.parent_record
        self.record = record  # type: RecordLayer
        self.anchor_id = self.build_anchor_id(region_feature)
        self.handlers = []  # type: List[AntismashModule]
        self.region_feature = region_feature  # type: Region

        self.cluster_blast = []  # type: List[Tuple[str, str]]
        self.knowncluster_blast = []  # type: List[Tuple[str, str]]
        self.subcluster_blast = []  # type: List[Tuple[str, str]]
        if self.region_feature.knownclusterblast:
            self.knowncluster_blast = self.knowncluster_blast_generator()
        if self.region_feature.subclusterblast:
            self.subcluster_blast = self.subcluster_blast_generator()
        if self.region_feature.clusterblast:
            self.cluster_blast = self.cluster_blast_generator()

        self.find_plugins_for_region()
        self.has_details = self.determine_has_details()
        self.has_sidepanel = self.determine_has_sidepanel()

    def __getattr__(self, attr: str) -> Any:
        if attr in self.__dict__:
            return super().__getattribute__(attr)
        return getattr(self.region_feature, attr)

    @property
    def best_knowncluster_type(self) -> str:
        """ The type of the best hit from knownclusterblast, if it was run """
        if not self.knownclusterblast:
            return ""
        result = self.knownclusterblast[0].cluster_type
        if result == "nrps":
            return "NRPS"
        return result

    @property
    def best_knowncluster_name(self) -> str:
        """ The name of the best hit from knownclusterblast, if it was run """
        if not self.knownclusterblast:
            return ""
        return self.knownclusterblast[0].name.replace("_", " ")

    @property
    def best_knowncluster_similarity(self) -> int:
        """ The percentage similarity of the best hit from knownclusterblast as
            an integer """
        if not self.knownclusterblast:
            return 0
        return self.knownclusterblast[0].similarity

    @property
    def bgc_id(self) -> str:
        """ The BGC id of the best hit from knownclusterblast, if it was run """
        if not self.region_feature.knownclusterblast:
            return ""
        return format(self.region_feature.knownclusterblast[0].bgc_id)

    @property
    def bgc_cluster_number(self) -> Optional[int]:
        """ The cluster number of the BGC id of the best knownclusterblast hit """
        if not self.region_feature.knownclusterblast:
            return None
        return self.region_feature.knownclusterblast[0].cluster_number

    @property
    def detection_rules(self) -> List[str]:
        """ The details of rules that caused the Region to be defined """
        return ["%s: %s" % (product, rules) for product, rules in
                zip(self.region_feature.products, self.region_feature.detection_rules)]

    def description_text(self) -> str:
        """ returns the Region description """
        description_text = 'Location: %s - %s nt. ' % (self.location.start + 1, self.location.end)
        if get_config().cf_create_clusters and self.probabilities:
            description_text += 'ClusterFinder probabilities: %s. ' % self.probabilities

        return description_text

    def cluster_blast_generator(self) -> List[Tuple[str, str]]:  # TODO: deduplicate
        """ Generates the details to use for clusterblast results """
        assert self.region_feature.clusterblast
        top_hits = self.region_feature.clusterblast[:self.record.options.cb_nclusters]
        results = []
        for i, label in enumerate(top_hits):
            i += 1  # 1-indexed
            svg_file = os.path.join('svg', 'clusterblast_r%dc%d_%s.svg' % (
                            self.record.record_index, self.get_region_number(), i))
            results.append((label, svg_file))
        return results

    def knowncluster_blast_generator(self) -> List[Tuple[str, str]]:
        """ Generates the details to use for knownclusterblast results """
        assert self.region_feature.knownclusterblast
        top_hits = self.region_feature.knownclusterblast[:self.record.options.cb_nclusters]
        results = []
        for i, summary in enumerate(top_hits):
            i += 1  # 1-indexed
            svg_file = os.path.join('svg', 'knownclusterblast_r%dc%d_%s.svg' % (
                            self.record.record_index, self.get_region_number(), i))
            results.append((summary.name, svg_file))
        return results

    def subcluster_blast_generator(self) -> List[Tuple[str, str]]:
        """ Generates the details to use for subclusterblast results """
        assert self.region_feature.subclusterblast
        assert self.region_feature.subclusterblast is not None, self.region_feature.location
        top_hits = self.region_feature.subclusterblast[:self.record.options.cb_nclusters]
        results = []
        for i, label in enumerate(top_hits):
            i += 1  # since one-indexed
            svg_file = os.path.join('svg', 'subclusterblast_r%dc%d_%s.svg' % (
                            self.record.record_index, self.get_region_number(), i))
            results.append((label, svg_file))
        return results

    def find_plugins_for_region(self) -> List[AntismashModule]:
        "Find a specific plugin responsible for a given Region type"
        for plugin in self.record.options.plugins:
            if not hasattr(plugin, 'will_handle'):
                continue
            if plugin.will_handle(self.products):
                self.handlers.append(plugin)
        return self.handlers

    def determine_has_details(self) -> bool:
        """ Sets has_details to be True if at least one plugin might create
            detail output for the Region
        """
        for handler in self.handlers:
            if hasattr(handler, "generate_details_div"):
                return True
        return False

    def determine_has_sidepanel(self) -> bool:
        """ Sets has_details to be True if at least one plugin might create
            output for the Region sidepanel
        """
        for handler in self.handlers:
            if hasattr(handler, "generate_sidepanel"):
                return True
        return False

    @staticmethod
    def build_anchor_id(region: Region) -> str:
        """ Builds a consistent HTML anchor identifier for a Region """
        return "r{}c{}".format(region.parent_record.record_index,
                               region.get_region_number())
