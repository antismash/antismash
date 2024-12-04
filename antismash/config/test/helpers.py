# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Helper functions for tests that need to make config-related decisions """

import pytest

from antismash.config.executables import find_executable_path


skip_without_fimo = pytest.mark.skipif(
    not find_executable_path("fimo"),
    reason="need fimo to run"
)

skip_without_fimo_and_meme = pytest.mark.skipif(
    not find_executable_path("fimo") and not find_executable_path("meme"),
    reason="need fimo and meme to run",
)
