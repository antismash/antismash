# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" A collection of classes for simplification of HTML templating.

    Each layer type is a different construct with helper functions that don't
    really belong in the original objects.
"""

from abc import ABC as AbstractClass, abstractmethod
import os
from typing import Any, List, Optional

from antismash.config import ConfigType
from antismash.common.html_renderer import Markup
from antismash.common.module_results import ModuleResults
from antismash.common.secmet import Record, Region
from antismash.custom_typing import AntismashModule


class AbstractRelatedArea(AbstractClass):
    """ An interface for modules to use in order to avoid circular import issues """
    @property
    @abstractmethod
    def identifier(self) -> str:
        """ The unique identifier of the related area, e.g. an accession. """

    @property
    @abstractmethod
    def description(self) -> str:
        """ A description of the related area. """

    @property
    def product(self) -> str:
        """ The product of the related area. This will be an empty string if not relevant. """
        return ""

    @property
    def similarity_percentage(self) -> Optional[int]:
        """ The similarity of the related area to the current region, as a percentage in the range 0 - 100.
            If a percentage is not appropriate, the value will be None.
        """
        return None

    @property
    def url(self) -> str:
        """ The URL of the related area. This will be an empty string if not relevant. """
        return ""


class BlankRelatedArea(AbstractRelatedArea):
    """ A class to represent the lack of a related area """
    @property
    def identifier(self) -> str:
        return ""

    @property
    def description(self) -> str:
        return ""


class OptionsLayer:
    """ A layer for the global Config options. """
    def __init__(self, options: ConfigType, modules: List[AntismashModule]) -> None:
        self.options = options
        self.plugins = modules

    def __getattr__(self, attr: str) -> Any:
        if attr in self.__dict__:
            return super().__getattribute__(attr)
        return getattr(self.options, attr)

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

    @property
    def output_basename(self) -> str:
        """ Returns the base filename for the output files """
        return os.path.basename(self.options.output_basename)


class RecordLayer:
    """ A layer for Record instances """
    def __init__(self, record: Record, results: Optional[ModuleResults], options: OptionsLayer) -> None:
        self.results = results
        self.seq_record = record
        self.options = options
        self.regions: List[RegionLayer] = []
        for region in record.get_regions():
            self.regions.append(RegionLayer(self, region))

    def __getattr__(self, attr: str) -> Any:
        if attr in self.__dict__:
            return super().__getattribute__(attr)
        return getattr(self.seq_record, attr)

    def get_name(self) -> str:
        """ Returns the ID of a record, with any extra notation included """
        name = self.seq_record.id
        if self.seq_record.has_multiple_sources():
            sources = self.seq_record.get_sources()
            name += f" (combined with {len(sources) - 1} other{'s' if len(sources) > 2 else ''})"
        return name

    def get_from_record(self) -> Markup:
        """ Returns the text to be displayed in the HTML overview page """

        current_id = f"<strong>{self.get_name()}</strong>"

        orig_id = self.seq_record.original_id

        if orig_id:
            if len(orig_id) < 40:
                current_id += f" (original name was: {orig_id})"
            else:
                current_id += f" (original name was: {orig_id[:60]}...)"

        source = self.annotations.get("source", "")
        if source:
            source = f" ({source})"

        return Markup(f"{current_id}{source}")


class RegionLayer:
    """ A layer for Region instances, contains special members for result of
        the clusterblast module
    """
    def __init__(self, record: RecordLayer, region_feature: Region) -> None:
        assert isinstance(region_feature, Region), type(region_feature)
        assert region_feature.parent_record
        self.record: RecordLayer = record
        self.anchor_id = self.build_anchor_id(region_feature)
        self.handlers: List[AntismashModule] = []
        self.region_feature: Region = region_feature

        self.most_related_area: AbstractRelatedArea = BlankRelatedArea()

        self.find_plugins_for_region()
        self.has_details = self.determine_has_details()
        self.has_sidepanel = self.determine_has_sidepanel()

    def __getattr__(self, attr: str) -> Any:
        if attr in self.__dict__:
            return super().__getattribute__(attr)
        return getattr(self.region_feature, attr)

    @property
    def detection_rules(self) -> List[str]:
        """ The details of rules that caused the Region to be defined """
        return [f"{product}: {rules}" for product, rules in
                zip(self.region_feature.products, self.region_feature.detection_rules)]

    def description_text(self) -> str:
        """ returns the Region description """
        description_text = (
            f"Location: {self.location.start + 1:,d} - {self.location.end:,d} nt."
            f" (total: {len(self.location):,d} nt)"
        )
        return description_text

    def find_plugins_for_region(self) -> List[AntismashModule]:
        "Find a specific plugin responsible for a given Region type"
        for plugin in self.record.options.plugins:
            if not hasattr(plugin, 'will_handle'):
                continue
            if plugin.will_handle(self.products, self.product_categories):
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

    def has_subregion_by_tool(self, tool: str) -> bool:
        """ Returns True if any subregion in the region was created with the
            given tool name
        """
        return any(sub.tool == tool for sub in self.subregions)

    @staticmethod
    def build_anchor_id(region: Region) -> str:
        """ Builds a consistent HTML anchor identifier for a Region """
        return f"r{region.parent_record.record_index}c{region.get_region_number()}"
