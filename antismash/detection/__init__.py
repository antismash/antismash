# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""
Contains modules based on detecting base information to be potentially used by
multiple analysis or output modules.

Detection modules are run in three stages:
- detection over the whole genome (e.g. Pfam domains) and not contributing to areas
     (unless later stages use the information found, of course)
- area formation, creating either subregions or protoclusters
- area refinement, creating subregions based on added protoclusters
- detection of particular features within an area
already formed.
"""

import enum

DetectionStage = enum.Enum("DetectionStage", [
    "FULL_GENOME",
    "AREA_FORMATION",
    "AREA_REFINEMENT",
    "PER_AREA",
])
