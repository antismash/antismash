# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

class ModuleResults:
    """
        For storage of individual module results. Should be subclassed by
        modules for a consistent API.
    """
    def __init__(self, record_id):
        self.record_id = record_id

    def to_json(self):
        """
            Converts the contained results into a json structure of simple types

            Returns a dict that should be parsable to recreate
             the ModuleResults instance
        """
        raise NotImplementedError()

    @staticmethod
    def from_json(json):
        """
            Converts the json structure back to a ModuleResults instance

            The ModuleResults instance returned should be able to regenerate the
            provided json by use of .to_json()
        """
        raise NotImplementedError()

    def add_to_record(self, record):
        """
            Stores relevant information from the results in the given record
        """
        raise NotImplementedError()
