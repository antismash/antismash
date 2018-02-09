# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" antiSMASH (antibiotics & Secondary Metabolite Analysis Shell):

    Genome mining and secondary metabolite biosynthetic gene cluster detection
    and analysis.

    The expected entry point as a library is run_antismash().
"""

from antismash.main import run_antismash, get_detection_modules, \
                           get_analysis_modules, get_output_modules, \
                           get_all_modules, __version__
