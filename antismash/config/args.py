# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

import argparse
import multiprocessing
import os
import threading
from collections import defaultdict

class AntiSmashParser(argparse.ArgumentParser):
    """Custom argument parser for antiSMASH
    """
    _show_all = False
    _displayGroup = {}

    def __init__(self, *args, **kwargs):
        """Initialisation method for the parser class"""
        self.parents = kwargs.get("parents", [])
        if self.parents:
            kwargs["parents"] = [parent.parser for parent in self.parents]
        kwargs["add_help"] = False
        super().__init__(*args, **kwargs)

    def add_argument_group(self, *args, **kwargs):
        basic = kwargs.get("basic", True) # True unless deliberately set False
        if not args:
            group = kwargs["title"]
        else:
            group = args[0]
        if group not in self._displayGroup:
            self._displayGroup[group] = []
        if "basic" in kwargs:
            del kwargs["basic"]
        if basic:
            self._displayGroup[group].extend(["basic"])
        if "param" in kwargs:
            self._displayGroup[group].extend(kwargs["param"])
            del kwargs["param"]
        return super().add_argument_group(*args, **kwargs)

    def print_help(self, file=None, show_all=False):
        self._show_all = show_all
        super(AntiSmashParser, self).print_help(file)

    def write_to_config_file(self, filename):
        def construct_arg_text(arg, dests):
            if arg.dest in ["help", "help_showall"] or arg.dest in dests:
                return []
            # skip postional args
            if not arg.option_strings:
                return []
            lines = []
            dests.add(arg.dest)
            flag = sorted(arg.option_strings, key=len, reverse=True)[0].lstrip('-')
            if "%(default)" in arg.help:
                help_text = arg.help % {"default": arg.default}
            else:
                help_text = arg.help
            lines.append("## {}".format(help_text))
            if arg.choices:
                lines.append("## Possible choices: {}".format(",".join(arg.choices)))
            if arg.default is not None:
                if not arg.default:
                    lines.append("#{}".format(flag))
                else:
                    lines.append("{} {}".format(flag, arg.default))
            else:
                lines.append("{}".format(flag))
            return ["\n".join(lines), "\n"]

        outfile = open(filename, 'w')
        dests = set()
        titles = defaultdict(lambda: defaultdict(lambda: []))
        for parent in self.parents:
            if getattr(parent, 'parser'):
                for arg in parent.parser._actions:
                    titles[parent.title][parent.prefix].extend(construct_arg_text(arg, dests))

        for arg in self._actions:
            titles["Core options"]["DUMMY"].extend(construct_arg_text(arg, dests))

        for title, prefixes in titles.items():
            banner = "#"*10 + "\n"
            outfile.writelines([banner,
                                "#\t{}\n".format(title),
                               banner, "\n"])
            for prefix, lines in prefixes.items():
                if prefix != "DUMMY":
                    outfile.write("[{}]\n".format(prefix))
                outfile.write("\n".join(lines))
            outfile.write("\n\n\n")
        outfile.close()

    def format_help(self):
        """Custom help format"""
        help_text = """
########### antiSMASH ver. {version} #############

{usage}

{args}
--------
Options
--------
{opts}
""".format(version="5 alpha", usage=self.format_usage(), args=self._get_args_text(), opts=self._get_opts_text())
        return help_text

    def format_usage(self):
        if self._show_all:
            formatter = self._get_formatter()
            formatter.add_usage(self.usage, self._actions, self._mutually_exclusive_groups)
            return formatter.format_help()
        return "usage: {prog} [-h] [options ..] sequence".format(prog=self.prog) + "\n"

    def _get_args_text(self):
        # fetch arg lists using formatter
        formatter = self._get_formatter()
        for action_group in self._action_groups:
            if action_group.title == "positional arguments":
                formatter.start_section("arguments")
                formatter.add_arguments(action_group._group_actions)
                formatter.end_section()
                break
        return formatter.format_help()

    def _get_opts_text(self):
        # fetch opt lists using formatter
        formatter = self._get_formatter()
        for action_group in self._action_groups:
            if action_group.title in ["optional arguments"]:
                formatter.add_arguments(action_group._group_actions)
        for action_group in self._action_groups:
            if action_group.title not in ["optional arguments", "positional arguments"]:
                show_opt = self._show_all
                if not show_opt:
                    if "basic" in self._displayGroup[action_group.title]:
                        show_opt = True
# TODO: keep lines or not?
#                    elif len(list(set(sys.argv) & set(self._displayGroup[action_group.title]))) > 0:
#                        show_opt = True
                if show_opt:
                    formatter.start_section(action_group.title)
                    if action_group.description is None:
                        action_group.description = ''
                    formatter.add_text(action_group.description)
                    formatter.add_arguments(action_group._group_actions)
                    formatter.end_section()
        return formatter.format_help()

    def convert_arg_line_to_args(self, line):
        """ overrides original to properly parse config files """
        line = line.strip()
        if not line:
            return []
        if line[0] in ["#", '[']:
            return []
        args = line.split()
        args[0] = "--" + args[0]
        return args

class FullPathAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, os.path.abspath(values))

class ModuleArgs:
    def __init__(self, title, prefix, override_safeties=False, always_on=True, # TODO: remove always_on when hmm_detection becomes core
                     enabled_by_default=False):
        self.title = title
        self.parser = AntiSmashParser(add_help=False)
# TODO: keep lines or not?
#        if kwargs and set(kwargs) not in {"override_safeties", "enabled_by_default"}:
#            raise ValueError("Unknown keyword arguments: %s" % list(kwargs))
        self.override = override_safeties #kwargs.get("override_safeties")
        self.enabled_by_default = enabled_by_default #kwargs.get("enabled_by_default")
        self.always_enabled = always_on
        self.group = self.parser.add_argument_group(title="Additional analysis", basic=self.override) # TODO
        self.options = self.parser.add_argument_group(title=title, basic=False)
        if not isinstance(prefix, str):
            raise TypeError("Argument prefix must be a string")
        self.prefix = prefix
        if len(self.prefix) < 2 and not self.override:
            raise ValueError("Argument prefixes must be at least 2 chars")
        self.skip_type_check = self.override
        self.single_arg = False
        self.args = []

    def add_option(self, name, *args, **kwargs):
        self._add_argument(self.options, name, *args, **kwargs)

    def add_analysis_toggle(self, name, *args, **kwargs):
        self._add_argument(self.group, name, *args, **kwargs)

    def _add_argument(self, group, name, *args, **kwargs):
        if self.single_arg:
            raise ValueError("Cannot add more arguments after an argument that is just the prefix name")
        if "default" not in kwargs:
            raise ValueError("Arguments must have a default, (default=)")
        if "type" not in kwargs:
            if "action" in kwargs: # most actions make types optional
                self.skip_type_check = True
                argtype = None
            else:
                raise ValueError("Arguments must have a type, (type=)")
        else:
            argtype = kwargs["type"]
        if "help" not in kwargs:
            raise ValueError("Arguments must have a description, (help=)")
        if "dest" not in kwargs:
            raise ValueError("Arguments must have a destination, (dest=)")
        default = kwargs["default"]
        dest = kwargs["dest"]
        self.verify_required(default, argtype, kwargs["help"])

        name, dest = self.process_names(name, dest)
        kwargs["dest"] = dest
        self.args.append(group.add_argument(name, *args, **kwargs))

    def verify_required(self, default, declared_type, description):
        if default is None:
            raise ValueError("Argument defaults may not be None")
        if not description:
            raise ValueError("Arguments must have a help string")

        if self.skip_type_check:
            return

        if declared_type is None:
            raise ValueError("Arguments must have a type declared")
        if not isinstance(default, declared_type):
            raise TypeError("Argument default doesn't match chosen type")

    def process_names(self, name, dest):
        if isinstance(name, list):
            raise TypeError("Module arguments may not have short-forms")
        try:
            name = str(name)
            dest = str(dest)
        except:
            raise TypeError("Argument name and dest must be strings")
        if self.override:
            return name, dest
        # ensure all destinations and flags start with the prefix
        if name.lstrip("-") == self.prefix:
            # not having anything else is ok if it's the only arg
            if self.args:
                raise ValueError("Arg name must not be just prefix with multiple args")
            self.single_arg = True
            if not dest:
                dest = self.prefix
        elif not name.lstrip("-").startswith(self.prefix + "-"):
            name = "--{}-{}".format(self.prefix, name.lstrip("-"))

        if not dest:
            dest = name.lstrip("--").replace("-", "_")
        elif dest == self.prefix:
            if not self.single_arg:
                print(name, dest, self.prefix)
                raise ValueError("Destination must include more information than the prefix")
        elif not dest.startswith(self.prefix + "_"):
            dest = "{}_{}".format(self.prefix, dest)
        if "-" in dest:
            raise ValueError("Destination for option cannot contain hyphens")

        return name, dest


def build_parser(from_config_file=False, modules=None):
    parents = [basic_options(), output_options(), advanced_options(),
               debug_options()]
    minimal = specific_debugging(modules)
    if minimal:
        parents.extend(minimal)
    if modules is not None:
        parents.extend(module.get_arguments() for module in modules)

    if from_config_file:
        parser = AntiSmashParser(parents=parents, fromfile_prefix_chars="@")
    else:
        parser = AntiSmashParser(parents=parents)

    #positional arguments
    parser.add_argument('sequences',
                        metavar='sequence',
                        nargs="*",
                        help="GenBank/EMBL/FASTA file(s) containing DNA.")

    #optional non-grouped arguments
    parser.add_argument('-h', '--help',
                        dest='help',
                        action='store_true',
                        default=False,
                        help="Show this help text.")
    parser.add_argument('--help-showall',
                        dest='help_showall',
                        action='store_true',
                        default=False,
                        help="Show full lists of arguments on this help text.")
    parser.add_argument('-c', '--cpus',
                        dest='cpus',
                        type=int,
                        default=multiprocessing.cpu_count(),
                        help="How many CPUs to use in parallel. (default: %(default)s)")

    ## grouped arguments

    return parser

def basic_options():
#    parser = AntiSmashParser(add_help=False)
    group = ModuleArgs("Basic analysis options", '', override_safeties=True)
#    group = parser.add_argument_group('Basic analysis options', '', basic=True)

    group.add_option('--taxon',
                       dest='taxon',
                       default='bacteria',
                       choices=['bacteria', 'fungi'],
                       type=str,
                       help="Taxonomic classification of input sequence. (default: %(default)s)")

    group.add_option('--input-type',
                       dest='input_type',
                       default='nucl',
                       choices=['nucl', 'prot'],
                       type=str,
                       help="Determine input type: amino acid sequence(s) or nucleotide sequence(s). (default: %(default)s)")
    return group

def output_options():
    group = ModuleArgs("Output options", 'output')
    group.add_option('--output-dir',
                       dest='output_dir',
                       default="",
                       type=str,
                       action=FullPathAction,
                       help="Directory to write results to.")
    return group

def advanced_options():
#    parser = AntiSmashParser(add_help=False)
    group = ModuleArgs("Advanced options", '', override_safeties=True)
    group.add_option('--reuse-results',
                       dest='reuse_results',
                       type=str,
                       action=FullPathAction,
                       default="",
                       metavar="PATH",
                       help="Use the previous results from the specified json datafile")
    group.add_option('--limit',
                       dest="limit",
                       type=int,
                       default=-1,
                       help="Only process the first <limit> records (default: %(default)s). -1 to disable")
    group.add_option('--minlength',
                       dest="minlength",
                       type=int,
                       default=1000,
                       help="Only process sequences larger than <minlength> (default: 1000).")
    group.add_option('--start',
                       dest='start',
                       type=int,
                       default=-1,
                       help="Start analysis at nucleotide specified.")
    group.add_option('--end',
                       dest='end',
                       type=int,
                       default=-1,
                       help="End analysis at nucleotide specified")
    group.add_option('--databases',
                       dest='database_dir',
                       default=os.path.join(os.path.dirname(os.path.dirname(__file__)), 'databases'),
                       type=str,
                       help="Root directory of the databases.")
    return group

def debug_options():
    group = ModuleArgs("Debugging & Logging options", '', override_safeties=True)
#    group = parser.add_argument_group("Debugging & Logging options", '', basic=True)
    group.add_option('-v', '--verbose',
                       dest='verbose',
                       action='store_true',
                       default=False,
                       help="Print verbose status information to stderr.")
    group.add_option('-d', '--debug',
                       dest='debug',
                       action='store_true',
                       default=False,
                       help="Print debugging information to stderr.")
    group.add_option('--logfile',
                       dest='logfile',
                       default=argparse.SUPPRESS,
                       type=str,
                       help="Also write logging output to a file.")
    group.add_option('--statusfile',
                       dest='statusfile',
                       default=argparse.SUPPRESS,
                       type=str,
                       help="Write the current status to a file.")
    group.add_option('--list-plugins',
                       dest='list_plugins',
                       action='store_true',
                       default=False,
                       help="List all available sec. met. detection modules.")
    group.add_option('--check-prereqs',
                       dest='check_prereqs_only',
                       action='store_true',
                       default=False,
                       help="Just check if all prerequisites are met.")
    group.add_option('--limit-to-record',
                       dest='limit_to_record',
                       default="",
                       metavar="record_id",
                       type=str,
                       help="Limit analysis to the record with ID record_id")
    group.add_option('--skip-cleanup',
                       dest='skip_cleanup',
                       action='store_true',
                       default=False,
                       help="Don't clean up temporary result files")
    group.add_option('-V', '--version',
                       dest='version',
                       action='store_true',
                       default=False,
                       help="Display the version number and exit.")
    return group

def specific_debugging(modules):
    if not modules:
        return None
    # only relevant for modules that are disabled by --minimal
    relevant_modules = []
    for module in modules:
        args = module.get_arguments()
        if args.always_enabled or not args.enabled_by_default:
            continue
        relevant_modules.append(module)
    if not relevant_modules:
        return None

    group = ModuleArgs('Debugging options for cluster-specific analyses', '', override_safeties=True)
    group.add_option('--minimal',
                       dest='minimal',
                       action='store_true',
                       default=False,
                       help="Only run core detection modules, no analysis modules unless explicitly enabled")
    errors = []
    for module in relevant_modules:
        try:
            group.add_option('--enable-%s' % (module.NAME),
                               dest='enabled_specific_plugins',
                               action='append_const',
                               const=module.NAME,
                               default=False,
                               help="Enable %s (default: enabled, unless --minimal is specified)" % module.SHORT_DESCRIPTION)
        except AttributeError as err:
            errors.append(str(err).replace("'module' object", module.__name__))
    if errors:
        raise AttributeError("\n\t".join([''] + errors))
    return group

def simple_options(module, args):
    modules = []
    if module is not None:
        modules = [module]
    return build_parser(from_config_file=False, modules=modules).parse_args(args)

class Config():
    __singleton = None
    __lock = threading.Lock()
    class _Config():
        def __init__(self, indict):
            if indict:
                self.__dict__.update(indict)

        def get(self, key, default=None):
            return self.__dict__.get(key, default)

        def __getattr__(self, attr):
            return self.__dict__[attr]

        def __setattr__(self, attr, value):
            raise RuntimeError("Config options can't be set directly")

        def __iter__(self):
            for i in self.__dict__.items():
                yield i

        def __repr__(self):
            return str(self)

        def __str__(self):
            return str(dict(self))

        def __len__(self):
            return len(self.__dict__)

    def __new__(cls, namespace=None):
        if namespace is None:
            values = {}
        elif isinstance(namespace, dict):
            values = namespace
        else:
            values = namespace.__dict__
        Config.__lock.acquire()
        if Config.__singleton is None:
            Config.__singleton = Config._Config(values)
        else:
            Config.__singleton.__dict__.update(values)
        Config.__lock.release()
        return Config.__singleton
