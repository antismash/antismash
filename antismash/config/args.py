# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

import argparse
from collections import defaultdict
import multiprocessing
import os


class AntismashParser(argparse.ArgumentParser):
    """ Custom argument parser for antiSMASH to override default help output

        Also helps with module level arguments

        Important changes
            parents argument: A list of ModuleArgs instead of other parsers
            add_help argument: Always overridden to False
    """
    def __init__(self, *args, **kwargs):
        # while keeping parents as ModuleArg instances is important
        # ArgumentParser itself expects parsers, so modify kwargs accordingly
        self.parents = kwargs.get("parents", [])
        if self.parents:
            kwargs["parents"] = [parent.parser for parent in self.parents]
        # ensure default -h/--help is disabled, since help is more complex here
        kwargs["add_help"] = False
        self._show_all = False
        # include in basic help any parent requesting it
        self._basic_help_groups = set(parent.title for parent in self.parents if parent.basic)
        # and since the additional analysis group doesn't propagate while empty
        self._basic_help_groups.add("Additional analysis")
        super().__init__(*args, **kwargs)

    def add_argument_group(self, title, description=None, basic=False, **kwargs):
        """ Overrides original to enable tracking of which help category a group
            is in.

            Arguments:
                title: the group label for showing help
                description: text with further information about the group
                basic: whether or not the group should be shown with just --help
                **kwargs: catches any argparse specific args that we don't use

            Returns:
                the group added
        """
        if not title:
            raise ValueError("Argument groups must have a group label/title")
        # if one group is added with a label and shows basic help, show for all
        # with the same label
        if basic:
            self._basic_help_groups.add(title)
        return super().add_argument_group(title, description, **kwargs)

    def print_help(self, file=None, show_all=False) -> None:  # arg added, so pylint: disable=arguments-differ
        """ Overrides parent print_help() to be able to pass through whether all
            help should be shown or not.

            Arguments:
                file: the file to write the help to, defaults to stdout
                show_all: True if all help should be shown, otherwise only basic

            Returns:
                None
        """
        self._show_all = show_all
        super().print_help(file)

    def write_to_config_file(self, filename, values=None) -> None:
        """ Write the default options to file in a form that can be parsed again

            Arguments:
                filename: the filename to write to
                values: a Namespace object with values to use instead of defaults
        """
        def construct_arg_text(arg, processed_destinations, values):
            """ Only construct if not already processed and not a help or input
                arg
            """
            if arg.dest in ["help", "help_showall", "write_config_file",
                            "reuse_results", "check_prereqs_only",
                            "list_plugins", "version", "genefinding_gff3"] \
                    or arg.dest in processed_destinations:
                return []
            # skip postional args
            if not arg.option_strings:
                return []
            lines = []
            processed_destinations.add(arg.dest)
            flag = sorted(arg.option_strings, key=len, reverse=True)[0].lstrip('-')
            # fill in the current default if relevant
            if "%(default)" in arg.help:
                help_text = arg.help % {"default": arg.default}
            else:
                help_text = arg.help
            lines.append("## {}".format(help_text))
            # add the set of possible choices, if relevant
            if arg.choices:
                lines.append("## Possible choices: {}".format(",".join(arg.choices)))
            # start with default value
            value = arg.default
            # and if it was supplied, use it
            if arg.dest in values:
                value = values.get(arg.dest)
            # convert options with list values to strings
            was_list = isinstance(value, list)
            if value and isinstance(value, list):
                value = ",".join(value)
            if isinstance(arg, argparse._StoreTrueAction):
                # don't write the value if it only changes to True,
                # otherwise the True will be considered a positional arg
                # and causes problems with the sequence
                state = "#"
                if value:
                    state = ""
                lines.append("{}{}".format(state, flag))
            elif not value:
                default = arg.default
                if not default:
                    default = arg.metavar
                lines.append("#{} {}".format(flag, default))
            else:
                default = ""
                if value == arg.default or was_list and value.split(",") == arg.default:
                    default = "#"
                lines.append("{}{} {}".format(default, flag, value))
            return ["\n".join(lines)]

        if not values:
            values = argparse.Namespace()

        outfile = open(filename, 'w')
        dests = set()  # set of processed destinations
        titles = defaultdict(lambda: defaultdict(list))
        for parent in sorted(self.parents, key=lambda group: group.title):
            if getattr(parent, 'parser'):
                for arg in parent.parser.get_actions():
                    titles[parent.title][parent.prefix].extend(construct_arg_text(arg, dests, values))

        for arg in self._actions:
            titles["Core options"]["core"].extend(construct_arg_text(arg, dests, values))

        for title, prefixes in titles.items():
            banner = "#"*10 + "\n"
            outfile.writelines([banner, "#\t{}\n".format(title), banner, "\n"])
            for prefix, lines in prefixes.items():
#                if prefix:
#                    outfile.write("[{}]\n".format(prefix))
                outfile.write("\n".join(lines))
                outfile.write("\n")
            outfile.write("\n\n")
        outfile.close()

    def format_help(self) -> str:
        """Custom help formatter"""
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

    def format_usage(self) -> str:
        """ Custom usage generator """
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
                    if action_group.title in self._basic_help_groups:
                        show_opt = True
                if show_opt:
                    formatter.start_section(action_group.title)
                    if action_group.description is None:
                        action_group.description = ''
                    formatter.add_text(action_group.description)
                    formatter.add_arguments(action_group._group_actions)
                    formatter.end_section()
        return formatter.format_help()

    def convert_arg_line_to_args(self, arg_line):
        """ overrides original to properly parse config files """
        line = arg_line.strip()
        if not line:
            return []
        # skip comments and section labels
        if line[0] in ["#", '[']:
            return []
        args = line.split()
        # prepend -- so the parser recognises it properly
        args[0] = "--" + args[0]
        return args

    def get_actions(self):
        """ a getter for _actions operated on by ArgumentParser, which may
            change behaviour """
        return tuple(self._actions)


class FullPathAction(argparse.Action):  # pylint: disable=too-few-public-methods
    """ An argparse.Action to ensure provided paths are absolute. """
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, os.path.abspath(values))


class ModuleArgs:
    """ The vehicle for adding module specific arguments in sane groupings.
        Each module should have a unique prefix for their arguments, for clarity
        and safety. The prefix only has to be supplied when constructing a
        ModuleArgs instance, it will be automatically applied to each argument
        added.

        To add arguments, use the following:
            for an analysis option (e.g. --pref-general):
                ModuleArgs.add_analysis_toggle('general', ...)
            for a module config option (e.g. --pref-max-size):
                ModuleArgs.add_option('max-size', ...)

        Arguments:
            title: The label shown on the help screen
            prefix: A prefix to use for all options used by the module
            override_safeties: If True, overrides many safety and consistency
                    checks, only exists to help placeholders and will be removed
            enabled_by_default: whether the module's analysis will be on by
                    default or whether it has to be specifically enabled,
                    if True, it will create a --enable-* arg to be used if
                    --minimal is provided.
            basic_help: Whether help for module options (not analysis toggles)
                    should be shown in the basic --help output.
    """
    def __init__(self, title, prefix, override_safeties=False,
                 enabled_by_default=False, basic_help=False):
        if not title:
            raise ValueError("Argument group must have a title")
        self.title = title
        self.parser = AntismashParser()
        self.override = override_safeties
        self.enabled_by_default = enabled_by_default
        # options for the module
        self.options = self.parser.add_argument_group(title=title, basic=basic_help)
        # the main analysis toggle(s)
        self.group = self.parser.add_argument_group(title="Additional analysis",
                                                    basic=True)
        if not isinstance(prefix, str):
            raise TypeError("Argument prefix must be a string")
        self.prefix = prefix
        if len(self.prefix) < 2 and not self.override:
            raise ValueError("Argument prefixes must be at least 2 chars")
        if self.prefix and not self.prefix.isalpha():
            raise ValueError("Argument prefixes must only be alphabetic")
        self.skip_type_check = self.override
        self.single_arg = False
        self.args = []
        self.basic = basic_help

    def add_option(self, name, *args, **kwargs) -> None:
        if not name:
            raise ValueError("Options must have a name")
        self._add_argument(self.options, name, *args, **kwargs)

    def add_analysis_toggle(self, name, *args, **kwargs) -> None:
        """ Add a simple on-off option to appear in the "Additional analysis"
            section. Every module that isn't running by default must have one
            of these arguments.
        """
        self._add_argument(self.group, name, *args, **kwargs)

    def _add_argument(self, group, name, *args, **kwargs) -> None:
        if self.single_arg:
            raise ValueError("Cannot add more arguments after an argument that is just the prefix name")
        # prevent the option name being considered destination by argparse
        if not name.startswith("-"):
            name = "-%s%s" % ("-" if len(name) > 1 else "", name)
        # most actions make types optional, so handle that
        if "type" not in kwargs and "action" in kwargs:
            self.skip_type_check = True
        self.verify_required(kwargs)

        name, kwargs["dest"] = self.process_names(name, kwargs["dest"])
        self.args.append(group.add_argument(name, *args, **kwargs))

    def verify_required(self, kwargs) -> None:
        """ Checks that all the important sections have arguments provided.

            If self.skip_type_check is truthy, extra checks on default types
            will be ignored.

            Raises relevant exceptions for any problems.

            Arguments:
                kwargs: the keyword args to check

            Returns:
                None
        """
        if kwargs.get("default") is None:
            raise ValueError("Arguments must have a default, (default=)")
        if not kwargs.get("help"):
            raise ValueError("Arguments must have a help string, (help=)")
        if not kwargs.get("dest"):
            raise ValueError("Arguments must have a destination, (dest=)")

        if self.skip_type_check:
            return

        declared_type = kwargs.get("type")
        if declared_type is None:
            raise ValueError("Arguments must have a type declared, (type=)")
        if not isinstance(kwargs["default"], declared_type):
            raise TypeError("Argument default doesn't match chosen type")

    def process_names(self, name, dest):
        if isinstance(name, list):
            raise TypeError("Module arguments may not have short-forms")
        if not isinstance(name, str) or not isinstance(dest, str):
            raise TypeError("Argument name and dest must be strings")

        # skip the remaining safety checks if set up
        # required for some core options only
        if self.override:
            return name, dest

        # ensure all destinations and flags start with the prefix
        if name.lstrip("-") == self.prefix:
            # not having anything else is ok if it's the only arg
            if self.args:
                raise ValueError("Arg name must not be just prefix if supporting multiple args")
            self.single_arg = True
            if not dest:
                dest = self.prefix
        elif not name.lstrip("-").startswith(self.prefix + "-"):
            name = "--{}-{}".format(self.prefix, name.lstrip("-"))

        if not dest:
            dest = name.lstrip("--").replace("-", "_")
        elif dest == self.prefix:
            if not self.single_arg:
                raise ValueError("Destination must include more information than the prefix")
        elif not dest.startswith(self.prefix + "_"):
            dest = "{}_{}".format(self.prefix, dest)
        if "-" in dest:
            raise ValueError("Destination for option cannot contain hyphens")

        return name, dest


def build_parser(from_config_file=False, modules=None) -> AntismashParser:
    """ Constructs an AntismashParser with all the default options antismash
        requires for proper operation, along with any added by the provided
        modules.

        Arguments:
            from_config_file: whether to allow loading from file with args like
                                @file
            modules: a list of modules/objects implementing get_arguments() for
                        construction of the arguments themselves

        Returns:
            an AntismashParser instance with options for all provided modules
    """
    parents = [basic_options(), output_options(), advanced_options(),
               debug_options()]
    minimal = specific_debugging(modules)
    if minimal:
        parents.append(minimal)
    if modules is not None:
        parents.extend(module.get_arguments() for module in modules)

    if from_config_file:
        parser = AntismashParser(parents=parents, fromfile_prefix_chars="@")
    else:
        parser = AntismashParser(parents=parents)

    # positional arguments
    parser.add_argument('sequences',
                        metavar='SEQUENCE',
                        nargs="*",
                        help="GenBank/EMBL/FASTA file(s) containing DNA.")

    # optional non-grouped arguments
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
    return parser


def basic_options():
    group = ModuleArgs("Basic analysis options", '', override_safeties=True, basic_help=True)

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
    group = ModuleArgs("Output options", 'output', basic_help=True)
    group.add_option('--output-dir',
                     dest='output_dir',
                     default="",
                     type=str,
                     action=FullPathAction,
                     help="Directory to write results to.")
    return group


def advanced_options():
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
                     help="Only process sequences larger than <minlength> (default: %(default)d).")
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
                     metavar="PATH",
                     type=str,
                     help="Root directory of the databases (default: %(default)s).")
    group.add_option('--write-config-file',
                     dest='write_config_file',
                     default="",
                     metavar="PATH",
                     action=FullPathAction,
                     type=str,
                     help="Write a config file to the supplied path")
    return group


def debug_options():
    group = ModuleArgs("Debugging & Logging options", '', override_safeties=True)
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
                     default="",
                     metavar="PATH",
                     type=str,
                     help="Also write logging output to a file.")
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
                     metavar="RECORD_ID",
                     type=str,
                     help="Limit analysis to the record with ID record_id")
    group.add_option('-V', '--version',
                     dest='version',
                     action='store_true',
                     default=False,
                     help="Display the version number and exit.")
    group.add_option('--profiling',
                     dest='profile',
                     action='store_true',
                     default=False,
                     help="Generate a profiling report, disables multiprocess python.")
    return group


def specific_debugging(modules):
    if not modules:
        return None
    # only relevant for modules that are disabled by --minimal
    relevant_modules = []
    for module in modules:
        args = module.get_arguments()
        if not args.enabled_by_default:
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
                             dest='%s_enabled' % (module.NAME),
                             action='store_true',
                             default=False,
                             help="Enable %s (default: enabled, unless --minimal is specified)" % module.SHORT_DESCRIPTION)
        except AttributeError as err:
            errors.append(str(err).replace("'module' object", module.__name__))
    if errors:
        raise AttributeError("\n\t".join([''] + errors))
    return group
