# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" The base Pairing class and all default promoter ranges for CASSIS """


class Pairing:
    """ A common component throughout CASSIS is a range given by number of
        promoters upstream and number of promoters downstream. This base class
        is used to reduce duplication and simplify string generation for labels.
    """
    def __init__(self, plus: int, minus: int) -> None:
        self.plus = int(plus)
        self.minus = int(minus)

    @property
    def pairing_string(self) -> str:
        """ A string representing the range of this pairing. """
        return f"+{self.plus:02d}_-{self.minus:02d}"

    def __str__(self) -> str:
        return self.pairing_string


# all possible promoter sets for motif detection
# plus --> include <plus> promoters downstream the anchor gene's promoter
# minus --> include <minus> promoters upstream the anchor gene's promoter
PROMOTER_RANGE = [Pairing(plus, minus) for plus in range(16) for minus in range(16) if plus + minus >= 3]
