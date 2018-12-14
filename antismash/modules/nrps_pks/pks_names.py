# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Helpers for converting PKS names into their short or long forms, e.g.
        mal -> Malonyl-CoA
"""

LONG_TO_SHORT = {'Malonyl-CoA': 'mal', 'Methylmalonyl-CoA': 'mmal', 'Methoxymalonyl-CoA': 'mxmal',
                 'Ethylmalonyl-CoA': 'emal', 'Isobutyryl-CoA': 'isobut', '2-Methylbutyryl-CoA': '2metbut',
                 'trans-1,2-CPDA': 'trans-1,2-CPDA', 'Acetyl-CoA': 'Acetyl-CoA', 'Benzoyl-CoA': 'benz',
                 'Propionyl-CoA': 'prop', '3-Methylbutyryl-CoA': '3metbut',
                 'CE-Malonyl-CoA': 'cemal', '2-Rhyd-Malonyl-CoA': '2Rhydmal', 'CHC-CoA': 'CHC-CoA',
                 'inactive': 'inactive'}
SHORT_TO_LONG = {val: key for key, val in LONG_TO_SHORT.items()}


def get_long_form(name: str, default: str = None) -> str:
    """ Convert a short form PKS monomer name into a long form name

        Arguments:
            name: the short form to convert
            default: the value to use if there is no mapping for the short form
                     (defaults to None, which reuses the short form)

        Returns:
            the long name if a mapping exists, otherwise the default value
    """
    return SHORT_TO_LONG.get(name, default or name)


def get_short_form(name: str, default: str = None) -> str:
    """ Convert a long form PKS monomer name into a short form name

        Arguments:
            name: the long form to convert
            default: the value to use if there is no mapping for the long form
                     (defaults to None, which reuses the long form)

        Returns:
            the short name if a mapping exists, otherwise the default value
    """
    return LONG_TO_SHORT.get(name, default or name)
