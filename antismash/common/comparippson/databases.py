# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Contains all the relevant database-related objects and functions for
    CompaRiPPson.

    Databases are supplied as a JSON file with database metadata and a FASTA file
    with simplified identifiers in the form "unique_id|accession or other
    additional context". Only the "unique_id" portion is used within the module,
    the remainder is for context and checking and should be small enough to not
    cause other issues with external dependencies.

    Database metadata includes required base information, then an "entries" object
    that contains all the relevant data for each entry in the FASTA file, keyed
    by that unique identifier.
"""

from dataclasses import dataclass, field
import glob
import json
import os
from typing import Any, Dict, List, Tuple

from antismash.common.html_renderer import Markup
from antismash.config import ConfigType, get_config

from .data_structures import Hit


@dataclass
class ComparippsonDatabase:
    """ Contains the vital information for a database, along with helper functions
        for generating output ina  the way the database specifies.

        The database metadata is cached once first accessed, as it is expected that
        the same database will be used by multiple RiPP modules. This cache can
        be explicitly cleared if necessary.
    """
    name: str
    version: str
    url: str
    id_format: str
    description_format: str
    fields: List[str]
    dir_name: str
    _reference_cache: Dict[str, Dict[str, str]] = field(default_factory=dict)

    def get_fasta_path(self, options: ConfigType) -> str:
        """ Returns the path to the FASTA file containing reference sequences """
        return os.path.join(options.database_dir, "comparippson", self.dir_name, self.version, "cores.fa")

    def get_metadata_path(self, options: ConfigType) -> str:
        """ Returns the path to the JSON file containing database metadata """
        return os.path.join(options.database_dir, "comparippson", self.dir_name, self.version, "metadata.json")

    @property
    def reference_cache(self) -> Dict[str, Dict[str, str]]:
        """ A cache of the metadata information for the database """
        if not self._reference_cache:
            with open(self.get_metadata_path(get_config()), encoding="utf-8") as handle:
                data = json.load(handle)
            self._reference_cache.update(data["entries"])
        return self._reference_cache

    def clear_cache(self) -> None:
        """ Clears the metadata cache """
        self._reference_cache.clear()

    def get_entry(self, identifier: str) -> Dict[str, str]:
        """ Returns the metadata matching the identifier for a particular sequence
            in the FASTA file for the database
        """
        return self.reference_cache[identifier]

    @classmethod
    def from_metadata_file(cls, metadata_path: str, dir_name: str) -> "ComparippsonDatabase":
        """ Builds an instance from a metadata file, along with the directory name
            for situations where the database name in the metadata is differently
            cased to the path on disk.
        """
        with open(metadata_path, encoding="utf-8") as handle:
            return cls.from_metadata(json.load(handle), dir_name)

    @classmethod
    def from_metadata(cls, data: Dict[str, Any], dir_name: str) -> "ComparippsonDatabase":
        """ Builds an instance from metadata, along with the directory name
            for situations where the database name in the metadata is differently
            cased to the path on disk.
        """
        fields = set(data["fields"])
        for key, entry in data["entries"].items():
            present = set(entry.keys())
            if present == fields:
                continue
            if present.issubset(fields):
                raise ValueError(f"CompaRiPPson DB entry has missing fields: entry {key}: {fields.difference(present)}")
            if fields.issubset(present):
                raise ValueError(f"CompaRiPPson DB entry has unknown fields: entry {key}: {present.difference(fields)}")
        data.pop("entries")
        return cls(**data, dir_name=dir_name)

    def to_json(self) -> Dict[str, str]:
        """ Converts the object, excluding caches, to a JSON-friendly dictionary """
        data = {key: val for key, val in vars(self).items() if key not in [
            "_reference_cache",
        ]}
        return data

    @classmethod
    def from_json(cls, data: Dict[str, Any]) -> "ComparippsonDatabase":
        """ Rebuilds an instance from a JSON-friendly dictionary """
        return cls(data["name"], data["version"], data["url"],
                   data["id_format"], data["description_format"], data["fields"],
                   data["dir_name"])

    def _build_text(self, fields: Dict[str, str], template: str) -> str:
        """ Fills in the given template with the fields provided """
        for field_name in self.fields:
            template = template.replace(f"@{field_name}@", fields[field_name])
        return template

    def build_identifier_for_hit(self, hit: Hit) -> Markup:
        """ Fills in the identifier template of the database with field values
            from the given hit. """
        return Markup.escape(self._build_text(hit.reference_fields, self.id_format))

    def build_description_for_hit(self, hit: Hit) -> Markup:
        """ Fills in the description template of the database with field values
            from the given hit. """
        return Markup.escape(self._build_text(hit.reference_fields, self.description_format))

    def build_url_for_hit(self, hit: Hit) -> Markup:
        """ Fills in the URL template of the database with field values from the
            given hit. """
        return Markup(self._build_text(hit.reference_fields, self.url))


_DATABASES: List[ComparippsonDatabase] = []


def find_latest_database_version(database_dir: str) -> str:
    """ Finds the most up-to-date database version in the given directory.
        Versions are expected to be in a XY.Z format, e.g. 3.0

        Arguments:
            database_dir: the path to the database directory

        Returns:
            the latest version number as a string, e.g. "3.0"
    """
    contents = glob.glob(os.path.join(database_dir, "*"))
    potentials: List[Tuple[float, str]] = []
    for name in contents:
        # only names in the form 2.0, 3.1, etc are valid
        try:
            version = os.path.basename(name)
            potentials.append((float(version), version))
        except ValueError:
            raise ValueError(f"Incompatible database version naming: {name}")
    if not potentials:
        raise ValueError(f"No matching database in location {database_dir}")
    latest = sorted(potentials)[-1]
    return latest[1]


def get_databases(options: ConfigType) -> List[ComparippsonDatabase]:
    """ Returns all databases found in the database directory. The results are
        cached, so modifiying the directory and calling the function again will
        not update results.

        Arguments:
            options: the antiSMASH config object from which the root of the database
                     directory is used

        Returns:
            a list of ComparippsonDatabase objects, one for each database found
    """
    if not _DATABASES:
        for subdir in glob.iglob(os.path.join(options.database_dir, "comparippson", "*")):
            version = find_latest_database_version(subdir)
            cores = os.path.join(subdir, version, "cores.fa")
            metadata = os.path.join(subdir, version, "metadata.json")
            if not os.path.exists(cores):
                raise RuntimeError(f"CompaRiPPson database directory present but missing data: {cores}")
            if not os.path.exists(metadata):
                raise RuntimeError(f"CompaRiPPson database directory present but missing data: {metadata}")
            db = ComparippsonDatabase.from_metadata_file(metadata, dir_name=os.path.basename(subdir))
            _DATABASES.append(db)
    for db in _DATABASES:
        assert isinstance(db, ComparippsonDatabase)
    return _DATABASES
