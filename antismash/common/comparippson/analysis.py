# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Contains the analysis portions of CompaRiPPson """

from collections import defaultdict
from dataclasses import dataclass, field
from tempfile import NamedTemporaryFile
from typing import Any, Dict, List, IO, Optional, Tuple

from antismash.common import path
from antismash.common.html_renderer import (
    FileTemplate,
    Markup,
    RIPP_CLASSES,
    spanned_sequence,
)
from antismash.common.subprocessing import execute
from antismash.config import ConfigType

from .data_structures import Hit, JsonConvertible
from .databases import ComparippsonDatabase as ComparippsonDB, get_databases


@dataclass(frozen=True)
class DBResults(JsonConvertible):
    """ Results for a specific database, with hits as a dictionary of query name
        to list of Hit objects
    """
    database: ComparippsonDB
    hits: Dict[str, List[Hit]]
    aliases: Dict[str, str]

    def to_json(self) -> Dict[str, Any]:
        return {
            "database": self.database.to_json(),
            "hits": {key: [hit.to_json() for hit in hits] for key, hits in self.hits.items()},
            "aliases": dict(self.aliases),
        }

    @classmethod
    def from_json(cls, data: Dict[str, Any]) -> "DBResults":
        hits = {key: [Hit.from_json(hit) for hit in hits] for key, hits in data["hits"].items()}
        return cls(ComparippsonDB.from_json(data["database"]), hits, data["aliases"])


@dataclass(frozen=True)
class ResultsByQuery:
    """ Results for a single query from multiple databases, with hits as a
        dictionary of database name to list of Hit objects
    """
    name: str
    hits: Dict[str, List[Hit]]


@dataclass
class MultiDBResults(JsonConvertible):
    """ Results for multiple queries from multiple databases
    """
    db_results: List[DBResults]
    aliases: Dict[str, str]
    by_query: Optional[Dict[str, ResultsByQuery]] = None
    db_by_name: Dict[str, ComparippsonDB] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.by_query is not None:
            return
        self.by_query = {}
        query_hits: Dict[str, Dict[str, List[Hit]]] = defaultdict(dict)
        for db_result in self.db_results:
            db = db_result.database
            self.db_by_name[db.name] = db
            for query, hits in db_result.hits.items():
                assert db.name not in query_hits[query]
                query_hits[query][db.name] = hits
        for query, hits_by_db in query_hits.items():
            self.by_query[query] = ResultsByQuery(query, hits_by_db)

    def to_json(self) -> Dict[str, Any]:
        return {
            "db_results": [results.to_json() for results in self.db_results],
            "aliases": self.aliases,
        }

    @classmethod
    def from_json(cls, data: Dict[str, Any]) -> "MultiDBResults":
        db_results = [DBResults.from_json(hits) for hits in data["db_results"]]
        return cls(db_results, data["aliases"])

    def get_html_for_query(self, name: str, colour_subset: str = None) -> Markup:
        """ Generates the relevant HTML for the given query

            Arguments:
                name: the name of the query as used to initially generate the results
                colour_subset: an optional subset of RIPP_CLASSES keys, in cases
                               where the defaults are not wanted

            Returns:
                a Markup instance with the constructed HTML
        """
        assert self.by_query is not None
        query = self.aliases.get(name, name)
        display_name = name
        if self.aliases.get(query):
            display_name = f"{name} and others with a shared core"

        if query not in self.by_query:
            return Markup("")

        if colour_subset:
            for char in colour_subset:
                if char not in RIPP_CLASSES:
                    expected = "".join(RIPP_CLASSES)
                    raise ValueError(f"Unknown character(s) '{char}' in subset of '{expected}'")

        def colour_sequence(core: str) -> Markup:
            if colour_subset is None:
                return spanned_sequence(core, class_mapping=RIPP_CLASSES)
            subset = {char: RIPP_CLASSES[char] for char in colour_subset}
            return spanned_sequence(core, class_mapping=subset)

        template = FileTemplate(path.get_full_path(__file__, "templates", "generic.html"))
        chunks = []
        for db_name, hits in self.by_query[query].hits.items():
            db = self.db_by_name[db_name]
            assert isinstance(db, ComparippsonDB), type(db)
            hits = sorted(hits, key=lambda x: (-x.match_count, x.reference.sequence,
                                               db.build_identifier_for_hit(x)))
            groups = []
            same = []
            current = hits[0]
            for hit in hits[1:]:
                if hit.reference.sequence == current.reference.sequence:
                    same.append(hit)
                    continue
                if same:
                    tail = f"and {len(same)} other"
                    if len(same) > 1:
                        tail += "s"
                else:
                    tail = ""
                groups.append((current, tail, same))
                current = hit
                same = []
            groups.append((current, f"and {len(same)} other{'s' if len(same) > 1 else ''}" if same else "", same))
            chunks.append(template.render(coloured_ripp_sequence=colour_sequence,
                                          name=display_name, db=db, groups=groups))
        return Markup("".join(chunks))


def _parse_simple_hits(handle: IO, threshold: float, database: ComparippsonDB) -> Dict[str, List[Hit]]:
    """ Parses Hits using the matching format """
    hits = defaultdict(list)
    for line in handle:
        parts = line.strip().split()
        ref_id = parts[Hit.REFERENCE_ID_INDEX].split("|", 1)[0]
        entry = database.get_entry(ref_id)
        hit = Hit.from_simple_blast_line(parts, entry)
        if hit.similarity > threshold:
            hits[hit.query_name].append(hit)
    return hits


def run_simple_blastp(database: ComparippsonDB, query_sequences: Dict[str, str],
                      threshold: float, options: ConfigType) -> Dict[str, List[Hit]]:
    """ Runs blastp over multiple query sequences with the given database.

        Arguments:
            database: the database to compare to
            query_sequences: query sequences as a dictionary of name to sequence
            threshold: minimum similarity as a fraction (i.e. 0.0 to 1.0)
            options: antiSMASH config object

        Returns:
            a dictionary mapping query name to a list of Hit objects
    """
    command = [
        options.executables.blastp,
        "-num_threads", str(options.cpus),
        "-db", database.get_fasta_path(options),
        "-outfmt", Hit.BLAST_FORMAT,
    ]

    fasta = "\n".join([f">{key}\n{val}" for key, val in query_sequences.items()])

    with NamedTemporaryFile() as temp:
        result = execute(command, stdin=fasta, stdout=temp)
        if not result.successful():
            error_line = result.stderr.splitlines()[-1]
            raise RuntimeError(f"blastp returned {result.return_code}: {error_line}")
        with open(temp.name, encoding="utf-8") as handle:
            hits = _parse_simple_hits(handle, threshold, database)
    return hits


def filter_duplicates(cores: Dict[str, str]) -> Tuple[Dict[str, str], Dict[str, str]]:
    """ Constructs a new dictionary containing only unique cores, along with a
        mapping of names from removed core to name of core that was kept.

        Arguments:
            cores: a dictionary mapping name to core sequence

        Returns:
            a tuple of
                a dictionary mapping name to core sequence
                a dictionary mapping removed name to name in first return value
    """
    core_to_name: Dict[str, str] = {}
    copies = {}
    for name, core in cores.items():
        if core in core_to_name:
            copies[name] = core_to_name[core]
        else:
            core_to_name[core] = name
    return {name: core for core, name in core_to_name.items()}, copies


def compare_precursor_cores(cores: Dict[str, str], options: ConfigType, threshold: float = 0.05,
                            databases: List[ComparippsonDB] = None) -> MultiDBResults:
    """ Compares the given cores to databases of cores.

        Arguments:
            cores: a dictionary of mapping name to amino sequence
            options: the antiSMASH config object
            threshold: the minimum similarity threshold for hits (in the range 0.0 - 1.0)
            databases: the list of database objects to use, defaults to all databases
                       present in the data directory

        Returns:
            a MultiDBResults object with all hits
    """
    results = []
    deduplicated_cores, aliases = filter_duplicates(cores)
    for db in databases or get_databases(options):
        hits = run_simple_blastp(db, deduplicated_cores, threshold, options)
        results.append(DBResults(db, hits, aliases))
    return MultiDBResults(results, aliases)
