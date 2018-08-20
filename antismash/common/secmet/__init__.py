# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" A more accessible and defined set of data structures for interacting with
    a sequence and its annotations.

    Can convert to and from Biopython, but that use is intended only to simplify
    reading and writing to standard output formats (e.g. genbank).
"""

from .record import Record
from .features import (
    AntismashDomain,
    CDSFeature,
    CDSMotif,
    Cluster,
    Feature,
    FeatureLocation,
    Gene,
    PFAMDomain,
    Prepeptide,
    Region,
    SubRegion,
    SuperCluster,
)
from .qualifiers import GeneFunction
