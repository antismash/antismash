# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

import unittest

from antismash.main import verify_options, gather_modules
from antismash.config.args import build_parser

class TestAntismash(unittest.TestCase):
    def setUp(self):
        self.default_options = build_parser(modules=gather_modules(with_genefinding=True)).parse_args()

    def test_default_options(self):
        # default options should work with all normal modules
        options = self.default_options
        assert verify_options(options)

        # adding an incompatibility should not be ok
        options.tta = True
        options.input_type = 'prot'
        assert not verify_options(options)

        # and just to be sure the above isn't just because tta
        options.input_type = 'nucl'
        assert verify_options(options)
