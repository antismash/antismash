# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

import unittest

from antismash.main import verify_options, get_detection_modules, get_analysis_modules
from antismash.config.args import build_parser

class TestAntismash(unittest.TestCase):
    def setUp(self):
        self.all_modules = get_detection_modules() + get_analysis_modules()
        self.default_options = build_parser(modules=self.all_modules).parse_args([])

    def test_default_options(self):
        # default options should work with all normal modules
        options = self.default_options
        assert verify_options(options, self.all_modules)

        # adding an incompatibility should not be ok
        options.tta = True
        options.input_type = 'prot'
        assert not verify_options(options, self.all_modules)

        # and just to be sure the above isn't just because tta
        options.input_type = 'nucl'
        assert verify_options(options, self.all_modules)
