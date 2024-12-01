# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Gene function classification based on the MITE database """

from collections import defaultdict
from dataclasses import dataclass
import pathlib
import tempfile
from typing import Any, Iterable, Optional, Self

from antismash.common import fasta, json, path
from antismash.common.secmet import GeneFunction, CDSFeature, ECGroup
from antismash.common.subprocessing import diamond
from antismash.config import ConfigType, get_config
from antismash.config.args import ModuleArgs

from .core import FunctionResults, Hit, Markup, Tool

TOOL_NAME = "MITE"
MINIMUM_IDENTITY = 60.


@dataclass(frozen=True, slots=True)
class MiteEntry:
    """ A container for a particular MITE entry's metadata  """
    accession: str
    description: str
    function: GeneFunction
    groups: tuple[ECGroup, ...]
    subfunctions: tuple[str, ...]
    version: str

    @classmethod
    def from_json(cls, data: dict[str, Any]) -> Self:
        """ Reconstructs an instance from the given data """
        return cls(
            accession=data["accession"],
            description=data["description"],
            function=GeneFunction.ADDITIONAL,
            groups=tuple(ECGroup(g) for g in data["groups"]),
            subfunctions=tuple(data["functions"]),
            version=data["version"],
        )


@dataclass(frozen=True, kw_only=True)
class Dataset:
    """ A container for information about a specific version of MITE's metadata """
    version: str
    database_path: str
    entries: dict[str, MiteEntry]
    url: str

    def get_entry(self, accession: str) -> MiteEntry:
        """ Returns the Mite metadata for the entry with the given accession """
        if accession not in self.entries:
            raise ValueError(f"no such entry: {accession}")
        return self.entries[accession]

    @classmethod
    def from_files(cls, base_dir: pathlib.Path) -> Self:
        """ Builds an instance from the given database file(s) """
        metadata_file = base_dir / "metadata.json"
        with metadata_file.open() as handle:
            raw = json.load(handle)
        entries = {str(accession): MiteEntry.from_json(entry) for accession, entry in raw["entries"].items()}
        return cls(
            database_path=str(base_dir / "mite.fasta"),
            entries=entries,
            url=raw["url"],
            version=raw["version"],
        )


_DATASETS: dict[str, Dataset] = {}


@dataclass(kw_only=True, slots=True)
class MiteHit(Hit):
    """ Contains information about a hit against a particular MITE entry """
    identity: float
    evalue: float
    bitscore: float

    def get_full_description(self) -> str:
        return f"{self.reference_id}: {self.description} ({int(self.identity)}% identity)"

    def get_html_fragment(self, metadata: dict[str, Any] | None = None, *, hide_id: bool = False) -> Markup:
        components = []
        if not hide_id:
            components.append(f"{self.query_id}: ")
        url = metadata.get("url", "").format(accession=self.reference_id) if metadata else ""
        if url:
            components.append(f'<a href="{url}">')
        components.append(self.reference_id)
        if url:
            components.append("</a>")
        components.append(f" (identity: {self.identity:0.2f}%) {self.description}")
        return Markup("".join(components))


class Results(FunctionResults[MiteHit]):
    """ Results of MITE gene function classification """
    def __init__(self, *, version: str, url: str, **kwargs: Any) -> None:
        super().__init__(**kwargs)
        self.version = version
        self.url = url

    def build_html_fragments(self) -> list[Markup]:
        metadata = {
            "url": self.url,
        }
        return [hit.get_html_fragment(metadata) for hit in self.best_hits.values()]

    @staticmethod
    def regenerate_hits(hits_by_name: dict[str, dict[str, Any]]) -> dict[str, MiteHit]:
        return {name: MiteHit.from_json(hit) for name, hit in hits_by_name.items()}

    def to_json(self) -> dict[str, Any]:
        data = super().to_json()
        data["version"] = self.version
        data["url"] = self.url
        return data

    @classmethod
    def from_json(cls, data: dict[str, Any]) -> Self:
        hits = cls.regenerate_hits(data.pop("best_hits", {}))
        mapping = {cds_name: GeneFunction.from_string(function)
                   for cds_name, function in data.pop("function_mapping").items()}
        return cls(version=data.pop("version"), url=data.pop("url"), best_hits=hits, function_mapping=mapping,
                   **data)


def add_arguments(args: ModuleArgs) -> None:
    """ Adds additional command-line arguments to the given argument group """
    args.add_option(
        "mite-version",
        dest="mite_version",
        type=str,
        default="latest",
        help="MITE database version number to use (e.g. 1.3) (default: %(default)s)."
    )


def check_prereqs(options: ConfigType) -> list[str]:
    """ Checks that all prerequisites are satisfied """
    problems = []
    for binary_name in ["diamond"]:
        if binary_name not in options.executables:
            problems.append(f"Failed to locate file: {binary_name!r}")

    try:
        available = False
        if get_dataset(options.genefunctions_mite_version):
            available = True
    except (FileNotFoundError, RuntimeError):
        pass

    if not available:
        problems.append(f"requested MITE version ({options.genefunctions_mite_version}) not present")

    # don't try to prepare data if there's already missing binaries or data
    if not problems:
        problems.extend(prepare_data(logging_only=True))

    return problems


def _get_blast_hits(cds_features: Iterable[CDSFeature], dataset: Dataset) -> list[diamond.Hit]:
    with tempfile.NamedTemporaryFile() as query_file:
        query_file.write(fasta.get_fasta_from_features(cds_features).encode())
        query_file.flush()
        raw = diamond.run_diamond_search(query_file.name, dataset.database_path, opts=[
            "--id", str(MINIMUM_IDENTITY),
        ])
    return diamond.parse_tabular_output(raw)


def classify(cds_features: Iterable[CDSFeature], options: ConfigType,
             ) -> Results:
    """ Finds possible classifications for the provided CDS features.

        Arguments:
            cds_features: the CDSFeatures to classify
            options: the config object

        Returns:
            a results instance with best hits and classification for each CDS
    """
    version = options.genefunctions_mite_version
    dataset = get_dataset(version)

    all_hits = _get_blast_hits(cds_features, dataset)

    # gather the hits
    blast_hits_by_cds_name = defaultdict(list)
    for hit in all_hits:
        subfunctions = dataset.entries[hit.reference_id].subfunctions
        blast_hits_by_cds_name[hit.query_id].append(MiteHit(
            query_id=hit.query_id,
            reference_id=hit.reference_id,
            description=", ".join(subfunctions),
            subfunctions=list(subfunctions),
            identity=hit.identity or 0.,
            bitscore=hit.bitscore or 0.,
            evalue=hit.evalue or 0.,
        ))
    hits_by_cds_name: dict[str, Hit] = {}
    function_mapping = {}
    group_mapping = {}
    subfunction_mapping = {}
    for cds_name, blast_hits in blast_hits_by_cds_name.items():
        if not blast_hits:
            continue
        ordered = sorted(blast_hits, reverse=True, key=lambda x: x.identity)
        best = ordered[0]
        entry = dataset.entries[best.reference_id]
        hits_by_cds_name[cds_name] = best
        function_mapping[cds_name] = entry.function
        group_mapping[cds_name] = entry.groups
        subfunction_mapping[cds_name] = entry.subfunctions

    return Results(tool=TOOL_NAME,
                   best_hits=hits_by_cds_name,
                   function_mapping=function_mapping,
                   group_mapping=group_mapping,
                   subfunction_mapping=subfunction_mapping,
                   version=dataset.version,
                   url=dataset.url)


def get_dataset(version: str = "latest") -> Dataset:
    """ Returns the metadata of the specified version. If the version is not defined,
        the version used will be the latest available in the database directory.

        Data is loaded from disk only on the first request.
    """
    if version in _DATASETS:
        return _DATASETS[version]

    data_dir = pathlib.Path(get_config().database_dir) / "mite"

    if version == "latest":
        requested = pathlib.Path(path.find_latest_database_version(str(data_dir),
                                                                   required_file_pattern="mite.fasta"))
        _DATASETS[str(requested)] = Dataset.from_files(data_dir / requested)
        _DATASETS[str(version)] = _DATASETS[str(requested)]
        version = str(requested)

    if version not in _DATASETS:
        _DATASETS[version] = Dataset.from_files(data_dir / version)

    return _DATASETS[version]


def prepare_data(logging_only: bool = False) -> list[str]:
    """ Prepares the MITE dataset

        Returns:
            a list of error messages
    """
    base_data_dir = get_config().database_dir
    if "mounted_at_runtime" in base_data_dir:
        return []
    mite = pathlib.Path(base_data_dir) / "mite"
    # check if the dataset is even available

    versions = [version for version in mite.glob("*") if version.is_dir()]
    if not versions:
        return ["No MITE versions available"]

    errors = []
    for version in versions:
        metadata = version / "metadata.json"
        sequences = version / "mite.fasta"
        db = version / "mite.db.dmnd"
        # ensure present
        if not metadata.exists():
            errors.append(f"{TOOL_NAME} directory for {version} present but missing metadata file")
        if not sequences.exists():
            errors.append(f"{TOOL_NAME} directory for {version} present but missing sequence file")
        errors.extend(diamond.check_diamond_files(str(metadata), str(sequences), str(db), logging_only=True))
    if errors and not logging_only:
        raise ValueError(f"Could not prepare data for {TOOL_NAME}: {', '.join(errors)}")
    return errors


def regenerate_results(data: Optional[dict[str, Any]]) -> Optional[Results]:
    """ Regenerates the tool's results, if possible """
    if not data:
        return None
    return Results.from_json(data)


TOOL = Tool(
    name=TOOL_NAME,
    add_arguments=add_arguments,
    check_prereqs=check_prereqs,
    classify=classify,
    prepare_data=prepare_data,
)
