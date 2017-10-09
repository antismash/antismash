# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.


class Signature(object):
    """Secondary metabolite signature"""
    def __init__(self, name, _type, description, cutoff, path):
        self.name = name
        self.type = _type
        self.description = description
        self.cutoff = cutoff
        self.path = path
