# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

from argparse import ArgumentError, Namespace
import os
import unittest

from helperlibs.wrappers.io import TemporaryDirectory

from antismash import get_all_modules
from antismash.config import get_config, update_config, destroy_config
from antismash.config import args


class TestConfig(unittest.TestCase):
    def setUp(self):
        self.core_parser = args.build_parser()
        modules = get_all_modules()
        self.default_parser = args.build_parser(modules=modules)

    def tearDown(self):
        destroy_config()

    def test_invalid_args(self):
        with self.assertRaises(SystemExit):
            self.core_parser.parse_args(["--an-invalid-switch"])

    def test_valid_args(self):
        # make sure args go through to the Namespace default object
        options = self.core_parser.parse_args(['--taxon', 'fungi'])
        assert options.taxon == 'fungi'
        # make sure they propagate to the Config singleton
        config = update_config(options)
        assert config.taxon == 'fungi'

    def test_namespace_initialisation(self):
        # test intialisation from namespace
        namespace = Namespace()
        namespace.taxon = 'fungi'
        config = update_config(namespace)
        assert config.taxon == 'fungi'
        # a new constructor should keep the value
        assert get_config().taxon == 'fungi'

    def test_dict_initialisation(self):
        config = update_config({'taxon': 'fungi'})
        assert config.taxon == 'fungi'
        # a new constructor should keep the value
        assert get_config().taxon == 'fungi'

    def test_assignment_proofing(self):
        config = update_config({'taxon': 'fungi'})
        assert config.taxon == 'fungi'
        # changing values in a Config object is invalid
        with self.assertRaises(RuntimeError):
            config.taxon = 'bacteria'
        # and verify it wasn't changed
        assert config.taxon == 'fungi'

    def test_get(self):
        config = update_config({'a': 1, 'b': None})
        # check attribute and get are the same
        assert config.a == config.get('a')
        # check default values function
        assert config.get('b', 3) is None  # since b exists
        assert config.get('c') is None  # since c doesn't
        assert config.get('c', 3) == 3  # now with default as 3

    def test_config_files(self):
        # TODO change some values in the file generated and check they're kept
        with TemporaryDirectory(change=True):
            args.build_parser()
            parser = args.build_parser(modules=get_all_modules())
            default_options = parser.parse_args([])
            parser.write_to_config_file("default_options.cfg")

            parser = args.build_parser(modules=get_all_modules(), from_config_file=True)
            from_file = parser.parse_args(["@default_options.cfg"])
            assert vars(default_options) == vars(from_file)

    def test_paths(self):
        options = self.core_parser.parse_args(["--reuse-results", "local"])
        assert os.sep in options.reuse_results


class TestModuleArgs(unittest.TestCase):
    def test_bad_values(self):
        with self.assertRaisesRegex(ValueError, "Argument prefixes must be alphanumeric"):
            mod_args = args.ModuleArgs('test options', 'prefix-has-dash')
        with self.assertRaisesRegex(ValueError, "Argument prefixes cannot start with numbers"):
            mod_args = args.ModuleArgs('test options', '2prefix')
        with self.assertRaisesRegex(TypeError, "Argument prefix must be a string"):
            mod_args = args.ModuleArgs('test options', 7)
        with self.assertRaisesRegex(ValueError, "Argument prefixes must be at least 2 chars"):
            args.ModuleArgs('test options', '')
        with self.assertRaisesRegex(ValueError, "Argument group must have a title"):
            args.ModuleArgs('', 'prefix')
        mod_args = args.ModuleArgs('test args', 'test')
        with self.assertRaisesRegex(ValueError, "Options must have a name"):
            mod_args.add_option('',
                                dest='some_place',
                                type=int,
                                default=0,
                                help="help")

        with self.assertRaisesRegex(ValueError, "Destination for option cannot contain hyphens"):
            mod_args.add_option('name',
                                dest='some-place',
                                type=int,
                                default=0,
                                help="help")
        with self.assertRaisesRegex(ValueError, "Arguments must have a default"):
            mod_args.add_option('name',
                                dest='some-place',
                                type=int,
                                help="help")
        with self.assertRaisesRegex(ValueError, "Arguments must have a type"):
            mod_args.add_option('name',
                                dest='some-place',
                                default="",
                                help="help")
        with self.assertRaisesRegex(ValueError, "Arguments must have a help"):
            mod_args.add_option('name',
                                dest='some-place',
                                default="",
                                type=str)
        with self.assertRaisesRegex(ValueError, "Arguments must have a destination"):
            mod_args.add_option('name',
                                help='help',
                                default="",
                                type=str)

    def test_bad_option_names(self):
        mod_args = args.ModuleArgs('test args', 'test')
        mod_args.add_option('--test', default="", type=str, help="no", dest="test")
        with self.assertRaisesRegex(ArgumentError, "argument --test: conflicting option string: --test"):
            mod_args.add_option('--test', default="", type=str, help="no", dest="test")

    def test_good_options(self):
        mod_args = args.ModuleArgs('test args', 'test')
        mod_args.add_option('test', default="", type=str, help="no", dest="test")
        mod_args.add_option('other', default="x", type=str, help="no", dest="test_other")
        parser = args.AntismashParser(parents=[mod_args])
        options = parser.parse_args(["--test", "thing"])
        assert options.test_other == "x"
        assert options.test == "thing"
        options = parser.parse_args(["--test", "1", "--test-other", "y"])
        assert options.test == "1"
        assert options.test_other == "y"

        mod_args = args.ModuleArgs('test', 't2pks')


# pytest only here for simple output capturing
def test_help(capsys):
    args.build_parser(modules=get_all_modules()).print_help()
    out, err = capsys.readouterr()
    assert "--minimal" not in out

    args.build_parser(modules=get_all_modules()).print_help(show_all=True)
    out_all, err_all = capsys.readouterr()
    assert err == err_all and not err
    # make sure show_all does something
    assert out != out_all
    assert "--minimal" in out_all
