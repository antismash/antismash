# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" A collection of functions and classes for running external commands.
"""

from .base import (
    execute,
    parallel_execute,
    parallel_function,
    RunResult,
)

from .blast import (
    run_blastp,
    run_blastp_version,
    run_makeblastdb_version,
)

from .diamond import (
    run_diamond,
    run_diamond_makedb,
    run_diamond_search,
    run_diamond_version,
)

from .hmmpfam import (
    run_hmmpfam2,
    run_hmmpfam2_version,
)

from .hmmpress import (
    run_hmmpress,
    run_hmmpress_version,
)

from .hmmscan import (
    run_hmmscan,
    run_hmmscan_version,
)

from .hmmsearch import (
    run_hmmsearch,
    run_hmmsearch_version,
)

from .java import (
    run_java_version,
)


from .memesuite import (
    run_fimo_simple,
    run_fimo_version,
    run_meme_version,
)

from .muscle import (
    run_muscle_single,
    run_muscle_version,
)

from .prodigal import (
    run_prodigal_version,
)
