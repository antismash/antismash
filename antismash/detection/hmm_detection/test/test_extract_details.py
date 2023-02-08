# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring,consider-using-with

from io import StringIO
import unittest

import antismash.detection.hmm_detection.data.extract_details as core


NO_TC_TEMPLATE = """HMMER3/f [3.1b2 | February 2015]
NAME  {name}
ACC   EX01234.5
DESC  {description}
LENG  120
NSEQ  240
//
"""

TC_TEMPLATE = """HMMER3/f [3.1b2 | February 2015]
NAME  {name}
ACC   EX01234.5
DESC  {description}
LENG  120
NSEQ  240
TC    {cutoff} 12.34;
//
"""


def create_profile(name: str, description: str, cutoff: float = -1.0) -> StringIO:
    """Create a fake HMM profile handle"""
    if cutoff >= 0:
        data = TC_TEMPLATE.format(name=name, description=description, cutoff=cutoff)
    else:
        data = NO_TC_TEMPLATE.format(name=name, description=description)

    handle = StringIO(data)
    handle.name = "path/to/fake.hmm"

    return handle


class ExtractDetailsTest(unittest.TestCase):
    def test_all_defaults(self):
        profile = create_profile("DEFAULTS", "All default values", 23.45)
        ret = core.run(profile)
        expected = "DEFAULTS\tAll default values\t23\tfake.hmm"
        assert ret == expected

    def test_no_cutoff(self):
        profile = create_profile("NO_CUTOFF", "No cutoff in this one")
        with self.assertRaisesRegex(RuntimeError, "No TC line found and no cutoff specified"):
            core.run(profile)

        profile.seek(0)
        ret = core.run(profile, cutoff=42)
        expected = "NO_CUTOFF\tNo cutoff in this one\t42\tfake.hmm"
        assert ret == expected

    def test_custom_name(self):
        profile = create_profile("IGNORED_NAME", "The original name is ignored", 23.45)
        ret = core.run(profile, name="CUSTOM_NAME")
        expected = "CUSTOM_NAME\tThe original name is ignored\t23\tfake.hmm"
        assert ret == expected

    def test_custom_descripton(self):
        profile = create_profile("IGNORED_DESC", "You won't read this", 23.45)
        ret = core.run(profile, description="Custom description here")
        expected = "IGNORED_DESC\tCustom description here\t23\tfake.hmm"
        assert ret == expected

    def test_all_custom(self):
        profile = create_profile("IGNORED", "Also ignored", 23.45)
        ret = core.run(profile, "CUSTOM_NAME", "Custom description", 17)
        expected = "CUSTOM_NAME\tCustom description\t17\tfake.hmm"
        assert ret == expected

    def test_tigr(self):
        profile = create_profile("TIGR01234", "Something: In TIGR formatting", 23.45)
        ret = core.run(profile)
        expected = "TIGR01234\tIn TIGR formatting\t23\tfake.hmm"
        assert ret == expected
