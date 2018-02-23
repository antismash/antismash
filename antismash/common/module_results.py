# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Provides a base for all module results objects to enable reusing previous
    antismash run results.
"""

from typing import Dict, List

from .secmet import Record
from .secmet.feature import ClusterBorder


class ModuleResults:
    """
        For storage of individual module results. Should be subclassed by
        modules for a consistent API.
    """
    __slots__ = ["record_id"]

    def __init__(self, record_id: str) -> None:
        self.record_id = record_id

    def to_json(self) -> Dict:
        """
            Converts the contained results into a json structure of simple types

            Returns a dict that should be parsable to recreate
             the ModuleResults instance
        """
        raise NotImplementedError()

    @staticmethod
    def from_json(json: Dict, record: Record) -> "ModuleResults":
        """
            Converts a json structure back to a ModuleResults instance

            The ModuleResults instance returned should be able to regenerate the
            provided json by use of .to_json()
        """
        raise NotImplementedError()

    def add_to_record(self, record: Record) -> None:
        """
            Stores relevant information from the results in the given record
        """
        raise NotImplementedError()

    def generate_images(self, image_dir: str) -> None:
        """
            Generates svgs and saves them to the output image directory
        """
        pass


class DetectionResults(ModuleResults):
    """ Stores results for detection modules.

        add_to_record() no longer requires overriding
        get_predictions() should be overridden if the module predicts cluster borders
    """
    def add_to_record(self, record: Record) -> None:
        pass

    def get_predictions(self) -> List[ClusterBorder]:
        """ Returns a list of ClusterBorder features predicted by the module """
        return []
