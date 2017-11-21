# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""

Rules must be provided to Parser as a single strings, one string for each rule.
All whitespace is treated only as a separator for symbols, so spaces or tabs are
fine. Having extra whitespace is also ok.

Line comments are possible with #, e.g. 'rule #comment'. Commented sections last
until a newline character is found (if any).

The minimum condition is always at the cluster level, if the count required is
great than the number of options, multiple clusters will have to provide one
of the options. A negated minimum is the equivalent of requiring less than the
given number of options.
Duplicate ids in the list will cause an error.

The minscore condition requires the presence of a signature with a bitscore
greater than or equal to the provided value, e.g.:
    minscore(a, 150)
A negated minscore condition will be true for a CDS if the signature is present
and has a score less than the value or the signature is is not present at all.

And and or both operate at both a CDS and cluster level. To limit to within a
single CDS use the cds options, e.g.:
    cds(a and b)
    cds(a and b) and c
    cds(a and (b or c))

cds(a and b) and cds(b and c) would be satisfied by two unique CDSs or a
single CDS with a, b and c

The logical operator precedence is 'not' > 'and' > 'or'. If you need to,
adding groups can collect lower precedence operators first, e.g.:
    a and (b or c).

The symbol 'cluster' cannot be a valid identifier to prevent confusion with the
old rule system where 'cluster' was a keyword. Using it will raise an error.

Rules cannot contain entirely negated conditions, as these cannot anchor clusters.
At least one positive requirement must exist.
Examples of negated rules:
    not a
    not a and not c
    not (a or c)
Examples of rules with positive requirements:
    a
    a and not c
    (a or c)
    (a or not c) # an 'or not' combination is effectively ignored for anchoring


Grammar hints:
    (e) -> all of
    [e] -> optional
    {e} -> repeated (+ >=1, * >= 0)
    label:e -> reference label for e

The grammar itself:
    LIST_OPEN = '['
    LIST_CLOSE = ']'
    GROUP_OPEN = '('
    GROUP_CLOSE = ')'
    INT = [-]{[1-9]}*[0-9]
    ID = [a-zA-Z]{[a-zA-Z0-9_-]}*
    MINIMUM_LABEL = 'minimum'
    CDS_LABEL = 'cds'
    SCORE_LABEL = 'minscore'

    BINARY_OP = ('and' | 'or')
    UNARY_OP = 'not'
    OPS = BINARY_OP | UNARY_OP

    RULE_MARKER = "RULE"
    COMMENT_MARKER = "COMMENT"
    CUTOFF_MARKER = "CUTOFF"
    EXTENT_MARKER = "EXTENT"
    CONDITIONS_MARKER = "CONDITIONS"


    RULE = RULE_MARKER classification:ID [COMMENT_MARKER:COMMENTS]
            CUTOFF_MARKER cutoff:INT EXTENT_MARKER extension:INT
            CONDITIONS_MARKER conditions:CONDITIONS

    CONDITIONS = CONDITION {BINARY_OP CONDITIONS}*;
    CONDITION =  [UNARY_OP] ( ID | CDS | MINIMUM | CONDITION_GROUP );
    CDS_CONDITION = [UNARY_OP] ID {BINARY_OP CDS_CONDITION}*;
    CONDITION_GROUP = GROUP_OPEN CONDITIONS GROUP_CLOSE;
    MINIMUM = MINIMUM_LABEL GROUP_OPEN
                count:INT COMMA
                LIST
              GROUP_CLOSE;
    SCORE = SCORE_LABEL GROUP_OPEN
                ID COMMA
                score:INT
              GROUP_CLOSE;
    LIST = LIST_OPEN contents:COMMA_SEPARATED_IDS LIST_CLOSE;
    COMMA_SEPARATED_IDS = ID {COMMA ID}*;
    CDS = GROUP_OPEN CDS_CONDITION GROUP_CLOSE;

cds(a and b) and not a and not cds(b or c))

Condition examples:
    a
    a and b
    (a and b)
    a or b
    not a and b
    not (a and b)
    a and not b
    not (a and not b)
    a or (b and c)
    (a or b) and cds(x and y)
    cds(x and (b or c))
    minimum(3, [a,b,c,d]) and (a or b) # 3 of a,b,c or d, with one being a or b
    minimum(3, [a,b,c,d]) and a and b # 3 of a,b,c or d, with one being a and one being b

Complete examples:
    RULE t1pks
        CUTOFF 20
        EXTENT 20
        CONDITIONS cds(PKS_AT and (PKS_KS or ene_KS
                                   or mod_KS or hyb_KS
                                   or itr_KS or tra_KS))

"""

import string
from enum import IntEnum
from operator import xor  # so type hints can be bool and not int
from typing import List, Set, Union

from antismash.detection.hmm_detection import signatures


class RuleSyntaxError(SyntaxError):
    """ Specifically for errors resulting from bad syntax in the rules being
        parsed.
    """
    pass


def is_legal_identifier(identifier) -> bool:
    """ Returns true if the identifier matches the form:
        [a-zA-Z]{[a-zA-Z0-9_-]}*
    """
    if not identifier[0].isalpha():
        return False
    for char in identifier:
        if not (char.isalpha() or char.isdigit() or char in ['_', '-']):
            return False
    # for now, keep cluster reserved to avoid confusion with previous system
    if identifier == "cluster":
        return False
    # to avoid confusion with minscore
    if identifier == "score":
        return False
    # keywords
    return True


def find_condition_identifiers(tokens) -> Set[str]:
    """ Finds and returns all identifiers within a condition section"""
    identifiers = set()
    in_conditions = False
    for token in list(tokens):
        if token.type == TokenTypes.CONDITIONS:
            in_conditions = True
        elif token.type.is_a_rule_keyword():
            in_conditions = False
        elif in_conditions and token.type == TokenTypes.IDENTIFIER:
            identifiers.add(token.identifier)
    return identifiers


class TokenTypes(IntEnum):
    """ Creates distinct values for each token. Each terminal has a value,
        along with variables like identifiers and numerical values
    """
    GROUP_OPEN = 1
    GROUP_CLOSE = 2
    LIST_OPEN = 3
    LIST_CLOSE = 4
    IDENTIFIER = 6
    MINIMUM = 7
    CDS = 8
    AND = 9
    OR = 10  # pylint thinks it's too short, so pylint: disable=invalid-name
    NOT = 11
    INT = 12
    COMMA = 13
    SCORE = 14
    RULE = 15
    COMMENT = 17
    CUTOFF = 18
    EXTENT = 19
    CONDITIONS = 20

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return str(self.name).lower()  # the str() is for pylint's sake

    def get_hit_string(self):
        return str(self)

    @classmethod
    def classify(cls, text):
        """ Returns the type of a token or raises an error if no classification
            possible
        """
        classification = Tokeniser.mapping.get(text)
        if classification is None:
            if text.isdigit() or text[0] == '-' and text[1:].isdigit():
                classification = cls.INT
            elif is_legal_identifier(text):
                classification = cls.IDENTIFIER
            else:
                raise RuleSyntaxError("Unclassifiable token: %s" % text)
        return classification

    def is_a_rule_keyword(self):
        """ Returns True if the token is a rule structure keyword such as
            RULE, COMMENT, CONDITIONS, etc
        """
        return 15 <= self.value <= 21


class Tokeniser:
    """ Converts a text block into a list of tokens
    """
    mapping = {"(": TokenTypes.GROUP_OPEN, ")": TokenTypes.GROUP_CLOSE,
               "[": TokenTypes.LIST_OPEN, "]": TokenTypes.LIST_CLOSE,
               "and": TokenTypes.AND, "or": TokenTypes.OR,
               "not": TokenTypes.NOT, ",": TokenTypes.COMMA,
               "minimum": TokenTypes.MINIMUM, "cds": TokenTypes.CDS,
               "minscore": TokenTypes.SCORE, "RULE": TokenTypes.RULE,
               "CONDITIONS": TokenTypes.CONDITIONS,
               "COMMENT": TokenTypes.COMMENT, "CUTOFF": TokenTypes.CUTOFF,
               "EXTENT": TokenTypes.EXTENT}

    def __init__(self, text):
        self.text = text
        self.tokens = []
        self.tokenise()
        self.current_symbol = []

    def tokenise(self):
        """ Does the work of separating tokens """
        self.current_symbol = []
        global_position = 0
        position = 0
        line = 1
        while global_position < len(self.text):
            char = self.text[global_position]
            # whitespace of any kind separates symbols of interest
            if char in string.whitespace:
                self._finalise(line, position)
                if char == "\n":
                    position = 0
                    line += 1
                else:
                    position += 1
                global_position += 1
                continue
            # these are their own tokens always, and always single chars
            if char in Tokeniser.mapping:
                self._finalise(line, position)
                self.tokens.append(Token(char, line, position))
            # part of a multi-char symbol
            elif char.isalnum() or char in ['-', '_']:
                self.current_symbol.append(char)
            # a comment, don't parse anything in the line
            elif char == "#":
                self._finalise(line, position)
                try:
                    # skip to after a newline if it exists
                    global_position = self.text.index('\n', global_position)
                    position = 0
                except ValueError:
                    # no newline after #, so we're finished
                    return
            else:
                raise RuleSyntaxError("Unexpected character in rule: %s\n%s\n%s^"
                                      % (char, self.text, " "*position))
            position += 1
            global_position += 1
        self._finalise(line, len(self.text))

    def _finalise(self, line, position):
        """ convert the current collection of chars into a Token """
        if not self.current_symbol:
            return
        self.tokens.append(Token("".join(self.current_symbol), line, position))
        self.current_symbol.clear()


class Token:
    """ Keeps the token details, the text, where it is in the total text block,
        and what type """
    def __init__(self, token_text: str, line_number: int, position: int) -> None:
        self.token_text = token_text
        self.type = TokenTypes.classify(token_text)
        self.line_number = line_number
        self.position = int(position)
        if len(token_text) > 1:
            self.position -= len(token_text)

    def __getattr__(self, key):
        if key == 'value':
            if self.type != TokenTypes.INT:
                raise AttributeError("Token is not numeric")
            return int(self.token_text)
        elif key == 'identifier':
            if self.type != TokenTypes.IDENTIFIER:
                raise AttributeError("Token has no identifier")
            return self.token_text
        return self.__dict__[key]

    def __repr__(self):
        if self.type == TokenTypes.IDENTIFIER:
            return "'{}'".format(self.token_text)
        if self.type == TokenTypes.INT:
            return self.token_text
        return str(self.type)


class Details:
    """ Keeps together a collection of useful contextual information when
        parsing
    """
    def __init__(self, cds, feats, results, cutoff):
        self.cds = cds  # str, name of cds that is being classified
        self.features_by_id = feats  # { id : feature }
        self.results_by_id = results  # { id : HSP list }
        self.possibilities = set(res.query_id for res in results.get(cds, []))
        self.cutoff = cutoff  # int

    def in_range(self, cds, other) -> bool:
        """ returns True if the two Locations are within cutoff distance

            this may be redundant if inputs are already limited, but here
            for safety
        """
        cds_start, cds_end = sorted([cds.start, cds.end])
        other_start, other_end = sorted([other.start, other.end])
        distance = min(abs(cds_end - other_start), abs(other_end - cds_start),
                       abs(cds_end - other_end), abs(other_start - cds_start))
        return distance < self.cutoff

    def just_cds(self, cds_of_interest) -> "Details":
        """ creates a new Details object with cds_of_interest as the focus

            the largest impact is Details.possibilities is updated
        """
        return Details(cds_of_interest, self.features_by_id, self.results_by_id,
                       self.cutoff)


class ConditionMet:
    def __init__(self, met: bool, matches: Union[Set, "ConditionMet"] = None) -> None:
        assert isinstance(met, bool)
        self.met = met
        self.matches = set()  # type: Set[str]
        if isinstance(matches, ConditionMet):
            self.matches = matches.matches
        elif matches:
            self.matches = matches
        assert isinstance(self.matches, set), type(matches)

    def __bool__(self):
        return self.met


class Conditions:
    """ The base condition case. Can, and probably will, contain sub conditions
        (e.g. (a and b or c) style groups).
    """
    def __init__(self, negated: bool, sub_conditions=None) -> None:
        self.negated = negated
        self.hits = 0
        assert self.negated in [False, True]
        if sub_conditions is None:
            sub_conditions = []
        self.sub_conditions = sub_conditions
        # just make sure that it's empty or that we have binary ops (a OR b..)
        assert not sub_conditions or len(sub_conditions) % 2 == 1
        if self.sub_conditions:
            assert all(isinstance(sub, Conditions) for sub in self.operands)
            assert all(isinstance(sub, TokenTypes) for sub in self.operators)
            for operator in self.operators:
                assert operator in [TokenTypes.AND, TokenTypes.OR]
            unique_operands = set()  # type: Set[str]
            for operand in map(str, self.operands):
                if operand in unique_operands:
                    raise ValueError("Rule contains repeated condition: %s\nfrom rule %s"
                                     % (operand, self))
                unique_operands.add(operand)

    @property
    def operands(self):
        return self.sub_conditions[::2]

    @property
    def operators(self):
        return self.sub_conditions[1::2]

    def are_subconditions_satisfied(self, details: Details, local_only=False) -> ConditionMet:
        """
            local_only limits the search to the single CDS in details
        """
        if len(self.sub_conditions) == 1:
            sub = self.sub_conditions[0].get_satisfied(details, local_only)
            return ConditionMet(sub.met, sub)

        # since ANDs are bound together, all subconditions we have here are ORs
        assert all(operator == TokenTypes.OR for operator in self.operators)
        # which means a simple any() will cover us
        sub_results = [sub.get_satisfied(details, local_only) for sub in self.operands]
        matching = set()  # type: Set[str]
        met = False
        for sub_result in sub_results:
            matching |= sub_result.matches
            met |= sub_result.met
        return ConditionMet(met, matching)

    def get_satisfied(self, details: Details, local_only=False) -> ConditionMet:
        satisfied = self.is_satisfied(details, local_only)
        if satisfied and satisfied.matches:
            self.hits += 1
        return satisfied

    def is_satisfied(self, details: Details, local_only=False) -> ConditionMet:
        """ Returns True if this condition is satisfied.
            Should be overridden in subclasses to suit their specific case.
        """
        subs = self.are_subconditions_satisfied(details, local_only)
        return ConditionMet(xor(self.negated, subs.met), subs)

    def get_hit_string(self) -> str:
        prefix = "not " if self.negated else ""
        if len(self.sub_conditions) == 1 \
                and not isinstance(self.sub_conditions[0], AndCondition):
            return "{}*({}{})".format(self.hits, prefix, self.sub_conditions[0].get_hit_string())
        return "{}*({}{})".format(self.hits, prefix, " ".join(sub.get_hit_string() for sub in self.sub_conditions))

    def contains_positive_condition(self) -> bool:
        # don't bother checking subs if this itself is negated
        if self.negated:
            return False
        # if no subconditions no need to check
        if not self.sub_conditions:
            return True
        return any(sub.contains_positive_condition() for sub in self.operands)

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        prefix = "not " if self.negated else ""
        if len(self.sub_conditions) == 1 \
                and not isinstance(self.sub_conditions[0], AndCondition):
            return "{}{}".format(prefix, self.sub_conditions[0])
        return "{}({})".format(prefix, " ".join(map(str, self.sub_conditions)))


class AndCondition(Conditions):
    """ Represents any number of chained AND operators.

        AndConditions can't be negated, as that is based on their operands, e.g.
        'a not and b' is bad sytax and 'not (a and b)' is negating the parent
        conditions group 'not (...)'.
        """

    def __init__(self, subconditions):
        super().__init__(False, subconditions)

    def is_satisfied(self, details: Details, local_only=False) -> ConditionMet:
        results = [sub.get_satisfied(details, local_only) for sub in self.operands]
        matched = set()  # type: Set[str]
        met = True
        for result in results:
            matched |= result.matches
            met = met and result.met
        return ConditionMet(met, matched)

    def get_hit_string(self) -> str:
        return "{}*({})".format(self.hits, " and ".join(op.get_hit_string() for op in self.operands))

    def __str__(self):
        return " and ".join(map(str, self.operands))


class MinimumCondition(Conditions):
    """ Represents the minimum() condition type"""
    def __init__(self, negated, count, options):
        self.count = count
        self.options = set(options)
        if len(self.options) != len(options):
            raise ValueError("Minimum conditions cannot have repeated options")
        if count < 1:
            raise ValueError("Minimum conditions must have a required count > 0")
        super().__init__(negated)

    def is_satisfied(self, details, local_only=False) -> ConditionMet:
        """ local_only is ignored here, since a MinimumCondition can't be inside
            a CDSCondition
        """
        hits = self.options.intersection(set(details.possibilities))
        hit_count = len(hits)
        if hit_count >= self.count:
            return ConditionMet(not self.negated, hits)

        current_cds = details.features_by_id[details.cds]

        # check to see if the remaining hits are in nearby CDSs
        for other_id, other_feature in details.features_by_id.items():
            if other_id == details.cds:
                continue
            if not details.in_range(current_cds.location, other_feature.location):
                continue
            other_options = [r.query_id for r in details.results_by_id.get(other_id, [])]
            hit_count += len(self.options.intersection(set(other_options)))
            if hit_count >= self.count:
                return ConditionMet(not self.negated, hits)
        return ConditionMet(self.negated, hits)

    def get_hit_string(self) -> str:
        return "{}*{}".format(self.hits, str(self))

    def __str__(self):
        return "{}minimum({}, [{}])".format("not " if self.negated else "",
                self.count, ", ".join(sorted(list(self.options))))


class CDSCondition(Conditions):
    """ Represents the cds() condition type """
    def is_satisfied(self, details, local_only=False) -> ConditionMet:
        # all child conditions have to be within a single CDS
        # so we force local_only to True
        satisfied_internally = super().are_subconditions_satisfied(details.just_cds(details.cds),
                local_only=True)
        # start with the current cds (and end if local_only or satisifed)
        if local_only or satisfied_internally:
            return ConditionMet(xor(self.negated, satisfied_internally.met), satisfied_internally)

        results = satisfied_internally.met
        # negative matches must also ensure all neighbours are negative
        start_feature = details.features_by_id[details.cds]
        for cds, feature in details.features_by_id.items():
            if details.cds == cds:
                continue
            if not details.in_range(start_feature.location, feature.location):
                continue
            results = results or super().are_subconditions_satisfied(details.just_cds(cds), local_only=True).met

        return ConditionMet(xor(results, self.negated))

    def get_hit_string(self) -> str:
        prefix = "not " if self.negated else ""
        return "{}*({}cds({}))".format(self.hits, prefix, " ".join(sub.get_hit_string() for sub in self.sub_conditions))

    def __str__(self):
        prefix = "not " if self.negated else ""
        return "{}cds({})".format(prefix, " ".join(map(str, self.sub_conditions)))


class SingleCondition(Conditions):
    """ Represents a single domain type that must be present """
    def __init__(self, negated, name):
        self.name = name
        super().__init__(negated)

    def is_satisfied(self, details: Details, local_only=False) -> ConditionMet:
        found_in_cds = self.name in details.possibilities
        # do we only care about this CDS? then use the smaller set
        if local_only:
            return ConditionMet(xor(self.negated, found_in_cds), set([self.name]))

        # we found all we were looking for, or we aren't allowed to look further
        if found_in_cds:
            cond = ConditionMet(xor(self.negated, found_in_cds), set([self.name]))
            return cond

        cds_feature = details.features_by_id[details.cds]
        # look at neighbours in range
        for other, other_hits in details.results_by_id.items():
            if other == details.cds:
                continue
            other_location = details.features_by_id[other].location
            if not details.in_range(cds_feature.location, other_location):
                continue
            other_possibilities = [res.query_id for res in other_hits]
            if self.name in other_possibilities:
                # a positive match, so we can exit early
                if not self.negated:
                    return ConditionMet(True)

        # if negated and we failed to find anything, that's a good thing
        return ConditionMet(self.negated)

    def get_hit_string(self) -> str:
        if self.negated:
            return "{}*(not {})".format(self.hits, self.name)
        return "{}*{}".format(self.hits, self.name)

    def __str__(self):
        return "{}{}".format("not " if self.negated else "", self.name)


class ScoreCondition(Conditions):
    """ Represents the minscore() condition """
    def __init__(self, negated: bool, name: str, score: int) -> None:
        self.name = name
        self.score = score
        super().__init__(negated)

    def is_satisfied(self, details: Details, local_only=False) -> ConditionMet:
        """ local_only is ignored since a ScoreCondition can't be inside a
            CDSCondition
        """
        # do we only care about this CDS? then use the smaller set
        found_in_cds = False
        match = set()
        if self.name in details.possibilities:
            for result in details.results_by_id[details.cds]:
                if result.query_id == self.name and result.bitscore >= self.score:
                    found_in_cds = True
                    match.add(self.name)
                    break

        if found_in_cds:
            return ConditionMet(not self.negated, match)

        cds_feature = details.features_by_id[details.cds]
        # look at neighbours in range
        for other, other_hits in details.results_by_id.items():
            other_location = details.features_by_id[other].location
            if not details.in_range(cds_feature.location, other_location):
                continue
            other_possibilities = [res.query_id for res in other_hits]
            if self.name in other_possibilities:
                for result in other_hits:
                    # a positive match, so we can exit early
                    if result.query_id == self.name and result.bitscore >= self.score:
                        return ConditionMet(not self.negated)

        # if negated and we failed to find anything, that's a good thing
        return ConditionMet(self.negated)

    def get_hit_string(self) -> str:
        if self.negated:
            return "{}*({})".format(self.hits, str(self))
        return "{}*{}".format(self.hits, str(self))

    def __str__(self):
        return "{}minscore({}, {})".format("not " if self.negated else "", self.name, self.score)


class DetectionRule:
    """ Contains all information about a rule, i.e.
            name: the label given for the rule
            cutoff: the cutoff used to construct clusters (in bases)
            extent: the extent to use to include nearby CDS features (in bases)
            conditions: the conditions that potential clusters have to satisfy
            comments: any comments provided in the rule
        """
    def __init__(self, name, cutoff, extent, conditions, comments=""):
        self.name = name
        self.cutoff = cutoff
        self.extent = extent
        self.conditions = conditions
        if not conditions.contains_positive_condition():
            raise ValueError("A rule's conditions must contain at least one positive requirement")
        self.comments = comments
        self.hits = 0

    def contains_positive_condition(self) -> bool:
        return self.conditions.contains_positive_condition()

    def detect(self, cds, feature_by_id, results_by_id) -> ConditionMet:
        """ Returns True if a cluster can be formed around this CDS
            using this rule
        """
        details = Details(cds, feature_by_id, results_by_id, self.cutoff)
        results = self.conditions.get_satisfied(details)
        if results and results.matches:
            self.hits += 1
        return results

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        condition_text = str(self.conditions)
        # strip off outer parens if they exist
        if condition_text[0] == "(" and condition_text[-1] == ')':
            condition_text = condition_text[1:-1]
        return "{}\t{}\t{}\t{}".format(self.name, self.cutoff // 1000,
                                       self.extent // 1000, condition_text)

    def reconstruct_rule_text(self) -> str:
        """ Generate a string that can be tokenised and parsed to recreate this
            rule
        """
        condition_text = str(self.conditions)
        # strip off outer parens if they exist
        if condition_text[0] == "(" and condition_text[-1] == ')':
            condition_text = condition_text[1:-1]
        comments = ""
        if self.comments:
            comments = "COMMENTS" + self.comments + " "
        return "RULE {} {}CUTOFF {} EXTENT {} CONDITIONS {}".format(self.name,
                    comments, self.cutoff // 1000, self.extent // 1000, condition_text)

    def get_hit_string(self) -> str:
        return self.conditions.get_hit_string()[3:-1]


# a typedef, even positions will be a Conditions instance, odd will be TokenTypes
ConditionList = List[Union[Conditions, TokenTypes]]  # pylint: disable=invalid-name


class Parser:
    """ Responsible for parsing an entire block of text. Rules parsed from the
        text are stored in the .rules member.
    """
    def __init__(self, text):
        self.lines = text.splitlines()
        self.rules = []
        self.current_line = 1
        self.current_token = None
        identifiers = set()
        tokens = Tokeniser(text.expandtabs()).tokens
        # gather all signature identifiers from condition blocks
        identifiers = find_condition_identifiers(tokens)
        # start the iterator up for parsing
        self.tokens = iter(tokens)
        try:
            self.current_token = next(self.tokens)
        except StopIteration:
            raise ValueError("No rules to parse")
        rule_names = set()
        while self.current_token and self.current_token.type == TokenTypes.RULE:
            rule = self._parse_rule()
            if rule.name in rule_names:
                raise ValueError("Multiple rules specified for the same rule name")
            rule_names.add(rule.name)
            self.rules.append(rule)
        if self.current_token:
            raise RuleSyntaxError("Expected RULE but found %s\n%s\n%s%s" % (
                                    self.current_token.type,
                                    self.lines[self.current_line - 1],
                                    " " * self.current_token.position, "^"))
        # verify gathered signature idenfiers exist as signatures
        unknown = identifiers - signatures.get_signature_names()
        if unknown:
            raise ValueError("Rules contained identifers without signatures: %s" % ", ".join(sorted(list(unknown))))

    def _consume(self, expected: TokenTypes) -> Token:
        self.current_line = self.current_token.line_number
        if self.current_token is None:
            raise RuleSyntaxError("Unexpected end of rule, expected %s" % expected)
        if self.current_token.type != expected:
            raise RuleSyntaxError("Expected %s but found %s\n%s\n%s%s" % (
                    expected, self.current_token.type,
                    self.lines[self.current_line - 1],
                    " "*self.current_token.position, "^"))
        consumed = self.current_token
        try:
            self.current_token = next(self.tokens)
            self.current_line = self.current_token.line_number
        except StopIteration:
            self.current_token = None
        return consumed

    def _consume_int(self) -> int:
        return self._consume(TokenTypes.INT).value

    def _consume_identifier(self) -> str:
        return self._consume(TokenTypes.IDENTIFIER).identifier

    def _parse_rule(self) -> DetectionRule:
        """ RULE = RULE_MARKER classification:ID [COMMENT_MARKER:COMMENTS]
                    CUTOFF_MARKER cutoff:INT EXTENT_MARKER extension:INT
                    CONDITIONS_MARKER conditions:CONDITIONS
        """
        self._consume(TokenTypes.RULE)
        rule_name = self._consume_identifier()
        comments = ""
        if self.current_token.type == TokenTypes.COMMENT:
            comments = self._parse_comments()
        self._consume(TokenTypes.CUTOFF)
        cutoff = self._consume_int() * 1000  # convert from kilobases
        self._consume(TokenTypes.EXTENT)
        extent = self._consume_int() * 1000
        self._consume(TokenTypes.CONDITIONS)
        conditions = Conditions(False, self._parse_conditions())
        if self.current_token is not None and self.current_token.type != TokenTypes.RULE:
            raise RuleSyntaxError("Unexpected symbol %s\n%s\n%s%s" % (
                    self.current_token.type,
                    self.lines[self.current_line - 1],
                    " "*self.current_token.position, "^"))
        return DetectionRule(rule_name, cutoff, extent, conditions,
                             comments=comments)

    def _parse_comments(self) -> str:
        """ COMMENTS = COMMENTS_MARKER comments
        """
        comment_tokens = []
        self._consume(TokenTypes.COMMENT)
        while not self.current_token.type == TokenTypes.CUTOFF:
            comment_tokens.append(self.current_token)
            self.current_token = next(self.tokens)
        return " ".join([token.token_text for token in comment_tokens])

    def _is_not(self) -> bool:
        """ [UNARY_OP]
        """
        negated = self.current_token.type == TokenTypes.NOT
        if negated:
            self._consume(TokenTypes.NOT)
        return negated

    def _parse_ands(self, lvalue, allow_cds) -> AndCondition:
        """ CONDITION and CONDITION { and CONDITION}
            ^ lvalue being passed in
        """
        and_conditions = [lvalue]
        and_conditions.append(self._consume(TokenTypes.AND).type)
        and_conditions.append(self._parse_single_condition(allow_cds))
        while self.current_token and self.current_token.type == TokenTypes.AND:
            and_conditions.append(self._consume(TokenTypes.AND).type)
            next_condition = self._parse_single_condition(allow_cds)
            and_conditions.append(next_condition)
        return AndCondition(and_conditions)

    def _parse_conditions(self, allow_cds=True, is_group=False) -> ConditionList:
        """    CONDITIONS = CONDITION {BINARY_OP CONDITIONS}*;
        """
        conditions = []  # type: List[Union[Conditions, TokenTypes]]
        lvalue = self._parse_single_condition(allow_cds)
        append_lvalue = True  # capture the lvalue if it's the only thing
        while self.current_token and self.current_token.type in [TokenTypes.AND,
                                                                 TokenTypes.OR]:
            if self.current_token.type == TokenTypes.AND:
                conditions.append(self._parse_ands(lvalue, allow_cds))
                append_lvalue = False
            else:
                if append_lvalue:
                    conditions.append(lvalue)
                self._consume(TokenTypes.OR)
                conditions.append(TokenTypes.OR)
                lvalue = self._parse_single_condition(allow_cds)
                append_lvalue = True
        if append_lvalue:
            conditions.append(lvalue)

        if is_group:
            if self.current_token is None:
                raise RuleSyntaxError("Unexpected end of rule, expected )\n%s" % (
                        self.lines[self.current_line - 1]))
            if self.current_token.type != TokenTypes.GROUP_CLOSE:
                raise RuleSyntaxError("Expected the end of a group, found %s\n%s\n%s^" % (
                        self.current_token,
                        self.lines[self.current_line - 1],
                        " " * self.current_token.position))
        elif self.current_token and self.current_token.type == TokenTypes.RULE:
            # this rule has ended since another is beginning
            return conditions
        elif self.current_token is not None:
            raise RuleSyntaxError("Unexpected symbol, found %s\n%s\n%s^" % (
                    self.current_token.token_text,
                    self.lines[self.current_line - 1],
                    " " * self.current_token.position))
        return conditions

    def _parse_single_condition(self, allow_cds) -> Conditions:
        """
            CONDITION = [UNARY_OP] ( ID | CONDITION_GROUP | MINIMUM | CDS );
            or we're in a CDS (i.e. allow_cds == False)
            CDS_CONDITION = [UNARY_OP] ID {BINARY_OP CDS_CONDITION}*;

        """
        negated = self._is_not()
        if self.current_token is None:
            raise RuleSyntaxError("Rules cannot end in not")
        if self.current_token.type == TokenTypes.GROUP_OPEN:
            return Conditions(negated, self._parse_group(allow_cds))
        elif allow_cds and self.current_token.type == TokenTypes.MINIMUM:
            return self._parse_minimum(negated=negated)
        elif allow_cds and self.current_token.type == TokenTypes.CDS:
            return CDSCondition(negated, self._parse_cds())
        elif self.current_token.type == TokenTypes.SCORE:
            return self._parse_score(negated=negated)
        return SingleCondition(negated, self._consume_identifier())

    def _parse_score(self, negated=False) -> ScoreCondition:
        """
            SCORE = minscore GROUP_OPEN ID COMMA INT GROUP_CLOSE
            e.g. minscore(trsC, 150)
        """
        self._consume(TokenTypes.SCORE)
        self._consume(TokenTypes.GROUP_OPEN)
        ident = self._consume_identifier()
        assert isinstance(ident, str)
        self._consume(TokenTypes.COMMA)
        score = self._consume_int()
        assert isinstance(score, int)
        self._consume(TokenTypes.GROUP_CLOSE)
        return ScoreCondition(negated, ident, score)

    def _parse_cds(self) -> ConditionList:
        """
            CDS = GROUP_OPEN CDS_CONDITION GROUP_CLOSE;
        """
        cds_token = self.current_token
        self._consume(TokenTypes.CDS)
        self._consume(TokenTypes.GROUP_OPEN)
        conditions = self._parse_conditions(allow_cds=False, is_group=True)
        if not conditions:
            raise RuleSyntaxError("cds conditions must have contents:\n%s\n%s^" % (
                    self.lines[self.current_line - 1],
                    " " * cds_token.position))
        self._consume(TokenTypes.GROUP_CLOSE)
        return conditions

    def _parse_group(self, allow_cds) -> ConditionList:
        """
            CONDITION_GROUP = GROUP_OPEN CONDITIONS GROUP_CLOSE;
        """
        self._consume(TokenTypes.GROUP_OPEN)
        conditions = self._parse_conditions(allow_cds, is_group=True)
        self._consume(TokenTypes.GROUP_CLOSE)
        return conditions

    def _parse_minimum(self, negated=False) -> MinimumCondition:
        """
            MINIMUM = MINIMUM_LABEL GROUP_OPEN
                  count:INT COMMA
                  LIST COMMA
                  GROUP_CLOSE;
        """
        initial_token = self.current_token
        self._consume(TokenTypes.MINIMUM)
        self._consume(TokenTypes.GROUP_OPEN)
        count = self._consume_int()
        self._consume(TokenTypes.COMMA)
        options = self._parse_list()
        self._consume(TokenTypes.GROUP_CLOSE)
        if count < 0:
            raise ValueError("minimum count must be greater than zero: \n%s\n%s^" % (
                             self.lines[self.current_line - 1],
                             " "*initial_token.position))
        return MinimumCondition(negated, count, options)

    def _parse_list(self) -> List[str]:
        """
            LIST = LIST_OPEN contents:COMMA_SEPARATED_IDS LIST_CLOSE;
        """
        self._consume(TokenTypes.LIST_OPEN)
        contents = [self._consume_identifier()]
        while self.current_token.type == TokenTypes.COMMA:
            self._consume(TokenTypes.COMMA)
            contents.append(self._consume_identifier())
        self._consume(TokenTypes.LIST_CLOSE)
        return contents
