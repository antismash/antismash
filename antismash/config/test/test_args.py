# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

import unittest
from argparse import Namespace

from antismash import gather_modules
import antismash.config.args as args

class TestConfig(unittest.TestCase):
    def setUp(self):
        self.core_parser = args.build_parser()
        self.default_parser = args.build_parser(modules=gather_modules())

    def tearDown(self):
        args.Config().__dict__.clear()

    def test_invalid_args(self):
        with self.assertRaises(SystemExit):
            self.core_parser.parse_args(["--an-invalid-switch"])

    def test_valid_args(self):
        # make sure args go through to the Namespace default object
        options = self.core_parser.parse_args(['--taxon', 'fungi'])
        assert options.taxon == 'fungi'
        # make sure they propagate to the Config singleton
        config = args.Config(options)
        assert config.taxon == 'fungi'

    def test_namespace_initialisation(self):
        # test intialisation from namespace
        namespace = Namespace()
        namespace.taxon = 'fungi'
        config = args.Config(namespace)
        assert config.taxon == 'fungi'
        # a new constructor should keep the value
        assert args.Config().taxon == 'fungi'

    def test_dict_initialisation(self):
        config = args.Config({'taxon': 'fungi'})
        assert config.taxon == 'fungi'
        # a new constructor should keep the value
        assert args.Config().taxon == 'fungi'

    def test_assignment_proofing(self):
        config = args.Config({'taxon': 'fungi'})
        assert config.taxon == 'fungi'
        # changing values in a Config object is invalid
        with self.assertRaises(RuntimeError):
            config.taxon = 'bacteria'
        # and verify it wasn't changed
        assert config.taxon == 'fungi'

    def test_get(self):
        config = args.Config({'a' : 1, 'b' : None})
        # check attribute and get are the same
        assert config.a == config.get('a')
        # check default values function
        assert config.get('b', 3) == None # since b exists
        assert config.get('c') is None # since c doesn't
        assert config.get('c', 3) == 3 # now with default as 3
