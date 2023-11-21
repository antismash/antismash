# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Abstract classes/interfaces for various feature types, allowing circular imports to be avoided
"""

from abc import ABC as AbstractClass, abstractmethod

from Bio.SeqRecord import SeqRecord

from .feature import Feature


class AbstractRegion(Feature, AbstractClass):
    """ An interface for Region features, with additional methods that will be required.
        Child collections are represented as simple Features here instead of CDSCollections,
        to avoid circularity with CDSFeatures. Any concrete implementation of the class will
        need to swap in the more accurate cases.
    """

    @property
    @abstractmethod
    def detection_rules(self) -> list[str]:
        """ Returns a list of unique detection rules collected from all
            contained CandidateClusters
        """

    @abstractmethod
    def get_region_number(self) -> int:
        """ Returns the region's numeric ID, only guaranteed to be consistent
            when the same clusters and subregions are defined in the parent record
        """

    @abstractmethod
    def get_product_string(self) -> str:
        """ Returns a string of all unique products collected from all child collections
        """

    @property
    @abstractmethod
    def products(self) -> list[str]:
        """ Returns a list of unique products collected from all child containers
        """

    @property
    @abstractmethod
    def product_categories(self) -> set[str]:
        """ Returns a list of unique product categories collected from all child collections
        """

    @abstractmethod
    def write_to_genbank(self, filename: str = None, directory: str = None,
                         record: SeqRecord = None,
                         ) -> None:
        """ Writes a genbank file containing only the information contained
            within the Region.
        """
