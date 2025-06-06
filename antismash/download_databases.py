# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""Script to download antiSMASH data files not shipped with the code."""

import argparse
import gzip
import hashlib
import lzma
import os
import pathlib
import sys
import tarfile
from typing import Any, Type
from urllib import error as urlerror
from urllib import request

import antismash
from antismash.common.hmmer import ensure_database_pressed
from antismash.common.html_renderer import (
    get_antismash_js_version,
    get_antismash_js_url,
)
from antismash.common import json

PFAM_LATEST_VERSION = "35.0"
PFAM_LATEST_URL = f"https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam{PFAM_LATEST_VERSION}/Pfam-A.hmm.gz"
PFAM_LATEST_ARCHIVE_CHECKSUM = "48ec2d1123c84046b00279eae1fb3d5be1b578e6221453f329d16954c89d0d35"
PFAM_LATEST_CHECKSUM = "8d3e2ffa785f91ee0e24a3994d2dcfff6f382e3cf663784a47688e7d95297fee"

CLUSTERBLAST_URL = "https://dl.secondarymetabolites.org/releases/clusterblast/clusterblast_4.0.tar.xz"
CLUSTERBLAST_ARCHIVE_CHECKSUM = "dc90a171751472d88088e12f52446c49fe0efff437d92110072151f523eb3ee3"
CLUSTERBLAST_FASTA_CHECKSUM = "55ea9e62dcc0c083b42f8209a48889f710e888db72d5433f4429d4a6ec1b946b"

KCB_URL = "https://dl.secondarymetabolites.org/releases/knownclusterblast/kcb_4.0.tar.xz"
KCB_ARCHIVE_CHECKSUM = "80f860cf947f6b0024a0414968484ba36de321991f70a11ca87e0ec034ae3233"
KCB_FASTA_CHECKSUM = "65b9a9bc8700d8fef54b80a2e57d88e7ca8e05f1b57a848bb6a5e8f47a3070b9"

CLUSTERCOMPARE_DBS = {
    "mibig": {
        "url": "https://dl.secondarymetabolites.org/releases/clustercompare/cc_mibig_4.0.tar.xz",
        "archive": "df009291f71d313d0e98260d1e9baae42db4efe1972d4fef502fd05355062ff5",
        "fasta": "6cad44254d8ad7e2e735fcd15eb2e40b526150db62d8dff68001e533c5aa22d2",
    },
}

COMPARIPPSON_DBS = {
    "asdb": {
        "archive_hash": "87796c22f5d13946222f94240e1d5896a42a7039907abce2ad858e9f83c8ac73",
        "fasta_hash": "af2a94ec525d75ff1dc31a9d3e74ec1653a05047546ac0961eac1d350e407fbe",
        "metadata_hash": "0235bfde62b6060d1cd790cba478228a268fea8e8fb36d582a20534cf785cd88",
        "url": "https://dl.secondarymetabolites.org/releases/comparippson/asdb_4.0.tar.xz",
        "version": "4.0",
    },
    "mibig": {
        "archive_hash": "cadeb37aaea894bb68833d672a853c0ae515644b0faf9107a2b45becae9ff2e5",
        "fasta_hash": "f5ece2ed3880e40b849aa999ed76445d796a71b97681021343b2a83ccf0dfde0",
        "metadata_hash": "4c591dae7be95622870530a2d911698350fceb951eba4d13121c5dc70923e365",
        "url": "https://dl.secondarymetabolites.org/releases/comparippson/mibig_4.0.tar.xz",
        "version": "4.0",
    },
}

MITE = {
    "url": "https://dl.secondarymetabolites.org/releases/mite/mite_1.3.tar.xz",
    "archive": "2a3e456700a720a078a00c6f3d0c9d25fb51274fa92f5e57e747cd4b3c2aa76d",
    "version": "1.3",
    "sequences": "c75144d41974b371fefc9875ee11dfc227839614cc1ce370ba4838c961354e6a",
    "metadata": "4a83153df89ae754533ad7701c7407fc5d4b50fd09fae3654e52cc14e43aa202",
}

RESFAM_URL = "https://dl.secondarymetabolites.org/releases/resfams/Resfams.hmm.gz"
RESFAM_ARCHIVE_CHECKSUM = "ae1e4be5a4a4ef5f669d36daddb94020df0ec41c448a21c642bdeca43abed342"
RESFAM_CHECKSUM = "780a0afdaa8c7a85d1b31367da66363c5b86ab8bc4a7288a8598f137cf55fbc0"

TIGRFAM_URL = "https://dl.secondarymetabolites.org/releases/tigrfam/TIGRFam.hmm.gz"
TIGRFAM_ARCHIVE_CHECKSUM = "0747a5efa85d594cd598e67a04611328ab865e6a9401eb143fce6e93edf22bd0"
TIGRFAM_CHECKSUM = "b905785603e1015bbe78d0235d0a6a131345bae35aaab7afbd8f7b59f7798842"

STACHELHAUS_URL = "https://dl.secondarymetabolites.org/releases/stachelhaus/1.1/signatures.tsv.xz"
STACHELHAUS_ARCHIVE_CHECKSUM = "90ea6d56503be456f1919769821fec884aff95e7a7290b32a5d8761d3ecd9f98"
STACHELHAUS_CHECKSUM = "eda79ad1e2c12791e70cb819721c109697295ecab5e8ba253ae658e1afb58cc5"

NRPS_SVM_URL = "https://dl.secondarymetabolites.org/releases/nrps_svm/2.0/models.tar.xz"
NRPS_SVM_ARCHIVE_CHECKSUM = "501bea0cb5048405c6bbf3e6cb05975e927117c4a080f02f452f86c0eb08c132"
NRPS_SVM_CHECKSUMS_COLLECTION = "9e53bfd3ab4076d52fcb651831561a2bb8325eff0a7a074f8287fe0ba4e78f54"

TRANSATOR_URL = "https://dl.secondarymetabolites.org/releases/transATor/transATor_2023.02.23.tar.xz"
TRANSATOR_ARCHIVE_CHECKSUM = "1e6cc5ab5aa0c70e26da47097520a49be754967f330c2857895b211e7601751c"
TRANSATOR_CHECKSUM = "6744065e91d3a1a43664420b38031ce02c4a762b3dffeed2162a1640e9acfceb"

LOCAL_FILE_PATH = os.path.abspath(os.path.dirname(__file__))

CHUNK = 128 * 1024


class DownloadError(RuntimeError):
    """Exception to throw when downloads fail."""

    pass  # pylint: disable=unnecessary-pass


def get_remote_filesize(url: str) -> int:
    """Get the file size of the remote file."""
    try:
        with request.urlopen(request.Request(url, method="HEAD")) as usock:
            dbfilesize = usock.info().get("Content-Length", "0")
    except urlerror.URLError:
        dbfilesize = "0"

    dbfilesize = int(dbfilesize)  # db file size in bytes
    return dbfilesize


def get_free_space(folder: str) -> int:
    """Return folder/drive free space (in bytes)."""
    return os.statvfs(folder).f_bfree * os.statvfs(folder).f_frsize


def check_diskspace(file_url: str) -> None:
    """Check if sufficient disk space is available."""
    dbfilesize = get_remote_filesize(file_url)
    free_space = get_free_space(".")
    if free_space < dbfilesize:
        raise DownloadError(
            "ERROR: Insufficient disk space available (required: {dbfilesize}, free: {free_space})."
        )


def download_file(url: str, filename: str) -> str:
    """Download a file."""
    try:
        req = request.urlopen(url)  # pylint: disable=consider-using-with
    except urlerror.URLError:
        raise DownloadError("ERROR: File not found on server.\nPlease check your internet connection.")

    # use 1 because we want to divide by the expected size, can't use 0
    expected_size = int(req.info().get("Content-Length", "1"))

    basename = os.path.basename(filename)
    dirname = os.path.dirname(filename)
    if not os.path.isdir(dirname):
        os.makedirs(dirname)

    overall = 0
    with open(filename, "wb") as handle:
        while True:
            try:
                chunk = req.read(CHUNK)
                if not chunk:
                    print("")
                    break
                overall += len(chunk)
                print(
                    f"\rDownloading {basename}: {(overall / expected_size) * 100:5.2f}% downloaded.",
                    end="",
                )
                handle.write(chunk)
            except IOError:
                raise DownloadError("ERROR: Download interrupted.")
    return filename


def checksum(filename: str, chunksize: int = 2 ** 20) -> str:
    """Get the SHA256 checksum of a file."""
    sha = hashlib.sha256()
    with open(filename, "rb") as handle:
        for chunk in iter(lambda: handle.read(chunksize), b""):
            sha.update(chunk)

    return sha.hexdigest()


def unzip_file(filename: str, decompressor: Any, error_type: Type[Exception]) -> str:
    """Decompress a compressed file."""
    newfilename, _ = os.path.splitext(filename)
    basename = os.path.basename(filename)
    try:
        zipfile = decompressor.open(filename, "rb")
        chunksize = 128 * 1024
        with open(newfilename, "wb") as handle:
            while True:
                try:
                    chunk = zipfile.read(chunksize)
                    if not chunk:
                        break
                    handle.write(chunk)
                except IOError:
                    raise DownloadError("ERROR: Unzipping interrupted.")
    except error_type:
        raise RuntimeError(
            f"Error extracting {basename}. Please try to extract it manually."
        )
    print(f"Extraction of {basename} finished successfully.")
    return newfilename


def untar_file(filename: str) -> None:
    """Extract a TAR/GZ file."""
    basename = filename.rpartition(os.sep)[2]
    try:
        with tarfile.open(filename) as tar:
            tar.extractall(path=filename.rpartition(os.sep)[0])
    except tarfile.ReadError:
        print(
            f"ERROR: Error extracting {basename}. Please try to extract it manually."
        )
        return
    print(f"Extraction of {basename} finished successfully.")


def delete_file(filename: str) -> None:
    """Delete a file."""
    try:
        os.remove(filename)
    except OSError:
        pass


def present_and_checksum_matches(filename: str, sha256sum: str) -> bool:
    """Check if a file is present and the checksum matches."""
    if os.path.exists(filename):
        print(f"Creating checksum of {os.path.basename(filename)}")
        csum = checksum(filename)
        if csum == sha256sum:
            return True
    return False


def download_if_not_present(url: str, filename: str, sha256sum: str) -> None:
    """Download a file if it's not present or checksum doesn't match."""
    # If we are missing the archive file, go and download
    if not present_and_checksum_matches(filename, sha256sum):
        download_file(url, filename)

    print(f"Creating checksum of {os.path.basename(filename)}")
    csum = checksum(filename)
    if csum != sha256sum:
        raise DownloadError(
            f"Error downloading {filename}, sha256sum mismatch. Expected {sha256sum}, got {csum}."
        )


def download_antismash_js(db_dir: str) -> None:
    """ Downloads the latest relevant version of the antiSMASH javascript """
    version = get_antismash_js_version()
    url = get_antismash_js_url()
    download_file(url, os.path.join(db_dir, "as-js", version, "antismash.js"))


def download_mite(db_dir: str) -> None:
    """ Downloads the latest converted MITE dataset """
    archive_filename = os.path.join(db_dir, "mite", "mite.tar.xz")
    version = MITE["version"]
    subdir = os.path.join(db_dir, "mite", version)
    if all([
        present_and_checksum_matches(os.path.join(subdir, "metadata.json"), MITE["metadata"]),
        present_and_checksum_matches(os.path.join(subdir, "mite.fasta"), MITE["sequences"]),
    ]):
        print(f"MITE {MITE['version']} already present")
        return

    print("Downloading MITE", MITE["version"])
    check_diskspace(MITE["url"])
    download_if_not_present(MITE["url"], archive_filename, MITE["archive"])
    filename = unzip_file(archive_filename, lzma, lzma.LZMAError)
    untar_file(filename)
    for name in ["metadata.json", "mite.fasta"]:
        assert os.path.exists(os.path.join(subdir, name))
    delete_file(filename)
    delete_file(filename + ".xz")


def download_pfam(db_dir: str, url: str, version: str, archive_checksum: str, db_checksum: str) -> None:
    """Download and compile the PFAM database."""
    archive_filename = os.path.join(db_dir, "pfam", version, "Pfam-A.hmm.gz")
    db_filename = os.path.splitext(archive_filename)[0]

    if present_and_checksum_matches(db_filename, db_checksum):
        print("PFAM file present and ok for version", version)
        return

    print("Downloading PFAM version", version)
    check_diskspace(url)
    download_if_not_present(url, archive_filename, archive_checksum)
    filename = unzip_file(archive_filename, gzip, gzip.zlib.error)  # type: ignore
    ensure_database_pressed(filename)
    delete_file(filename + ".gz")


def download_resfam(db_dir: str) -> None:
    """Download and sanitise the Resfam database."""
    archive_filename = os.path.join(db_dir, "resfam", "Resfams.hmm.gz")
    filename = os.path.splitext(archive_filename)[0]

    if present_and_checksum_matches(filename, RESFAM_CHECKSUM):
        print("Resfam file present and checked")
        return

    print("Downloading Resfam database")
    check_diskspace(RESFAM_URL)
    download_if_not_present(RESFAM_URL, archive_filename, RESFAM_ARCHIVE_CHECKSUM)
    filename = unzip_file(archive_filename, gzip, gzip.zlib.error)  # type: ignore
    delete_file(filename + ".gz")
    ensure_database_pressed(filename)


def download_tigrfam(db_dir: str) -> None:
    """Download the TIGRFam database."""
    archive_filename = os.path.join(db_dir, "tigrfam", "TIGRFam.hmm.gz")
    filename = os.path.splitext(archive_filename)[0]

    if present_and_checksum_matches(filename, TIGRFAM_CHECKSUM):
        print("TIGRFam database present and checked")
    else:
        print("Downloading TIGRFam database")
        check_diskspace(TIGRFAM_URL)
        download_if_not_present(TIGRFAM_URL, archive_filename, TIGRFAM_ARCHIVE_CHECKSUM)
        filename = unzip_file(archive_filename, gzip, gzip.zlib.error)  # type: ignore
        delete_file(archive_filename)

    ensure_database_pressed(filename)


def download_clusterblast(db_dir: str) -> None:
    """Download the clusterblast database."""
    archive_filename = os.path.join(db_dir, CLUSTERBLAST_URL.rpartition("/")[2])
    fasta_filename = os.path.join(db_dir, "clusterblast", "proteins.fasta")

    if present_and_checksum_matches(fasta_filename, CLUSTERBLAST_FASTA_CHECKSUM):
        print("ClusterBlast fasta file present and checked")
        return

    print("Downloading ClusterBlast database.")
    check_diskspace(CLUSTERBLAST_URL)
    download_if_not_present(CLUSTERBLAST_URL, archive_filename, CLUSTERBLAST_ARCHIVE_CHECKSUM)
    filename = unzip_file(archive_filename, lzma, lzma.LZMAError)
    untar_file(filename)
    delete_file(filename)
    delete_file(filename + ".xz")


def download_clustercompare(db_dir: str, name: str, url: str, archive: str, fasta: str) -> None:
    """Download a ClusterCompare database."""
    version = url.rsplit("_", 1)[1].split(".tar", 1)[0]  # e.g. cc_mibig_3.1.tar.xz -> 3.1
    archive_filename = os.path.join(db_dir, "clustercompare", name, url.rpartition("/")[2])
    fasta_filename = os.path.join(db_dir, "clustercompare", name, version, "proteins.fasta")

    if present_and_checksum_matches(fasta_filename, fasta):
        print(f"ClusterCompare {name} FASTA file present and checked")
        return

    print(f"Downloading ClusterCompare {name} database.")
    check_diskspace(url)
    download_if_not_present(url, archive_filename, archive)
    filename = unzip_file(archive_filename, lzma, lzma.LZMAError)
    untar_file(filename)
    assert os.path.exists(fasta_filename)
    delete_file(filename)
    delete_file(filename + ".xz")


def download_comparippson_db(db_dir: str, name: str, url: str, version: str,
                             archive_hash: str, fasta_hash: str, metadata_hash: str) -> None:
    """Download a CompaRiPPson database."""
    archive_filename = os.path.join(db_dir, "comparippson", name, version, url.rpartition("/")[2])
    fasta_filename = os.path.join(db_dir, "comparippson", name, version, "cores.fa")
    metadata_filename = os.path.join(db_dir, "comparippson", name, version, "metadata.json")

    if all([
        present_and_checksum_matches(fasta_filename, fasta_hash),
        present_and_checksum_matches(metadata_filename, metadata_hash),
    ]):
        print(f"CompaRiPPson {name} {version} files present and checked")
        return

    print(f"Downloading CompaRiPPson {name} {version} database.")
    check_diskspace(url)
    download_if_not_present(url, archive_filename, archive_hash)
    filename = unzip_file(archive_filename, lzma, lzma.LZMAError)
    untar_file(filename)
    assert os.path.exists(fasta_filename)
    assert os.path.exists(metadata_filename)
    delete_file(filename)
    delete_file(filename + ".xz")


def download_knownclusterblast(db_dir: str) -> None:
    """Download the knownclusterblast database."""
    version = KCB_URL.rsplit("_", 1)[1].split(".tar", 1)[0]  # e.g. kcb_3.1.tar.xz -> 3.1
    archive_filename = os.path.join(db_dir, "knownclusterblast", KCB_URL.rpartition("/")[2])
    fasta_filename = os.path.join(db_dir, "knownclusterblast", version, "proteins.fasta")

    if present_and_checksum_matches(fasta_filename, KCB_FASTA_CHECKSUM):
        print("KnownClusterBlast fasta file present and checked")
        return

    print("Downloading KnownClusterBlast database.")
    check_diskspace(KCB_URL)
    download_if_not_present(KCB_URL, archive_filename, KCB_ARCHIVE_CHECKSUM)
    filename = unzip_file(archive_filename, lzma, lzma.LZMAError)
    untar_file(filename)
    delete_file(filename)
    delete_file(filename + ".xz")


def download_stachelhaus(db_dir: str) -> None:
    """Download the Stachelhaus signatures file."""
    version, name_part = STACHELHAUS_URL.split("/")[-2:]
    dirname = os.path.join(db_dir, "nrps_pks", "stachelhaus", version)
    archive_filename = os.path.join(dirname, name_part)
    signatures_filename, _ = os.path.splitext(archive_filename)

    if present_and_checksum_matches(signatures_filename, STACHELHAUS_CHECKSUM):
        print("Stachelhaus signature table file present and checked")
        return

    print("Downloading Stachelhaus signature table.")
    check_diskspace(STACHELHAUS_URL)
    download_if_not_present(STACHELHAUS_URL, archive_filename, STACHELHAUS_ARCHIVE_CHECKSUM)
    filename = unzip_file(archive_filename, lzma, lzma.LZMAError)
    delete_file(filename + ".xz")


def download_nprs_svm(db_dir: str) -> bool:
    """Download NRPS SVM models"""
    version, name_part = NRPS_SVM_URL.split("/")[-2:]
    dirname = pathlib.Path(db_dir) / "nrps_pks" / "svm" / version
    archive_filename = dirname / name_part
    checksum_file = dirname / "checksums.json"

    if present_and_checksum_matches(str(checksum_file), NRPS_SVM_CHECKSUMS_COLLECTION):
        with checksum_file.open("r", encoding="utf-8") as handle:
            checksums: dict[str, str] = json.load(handle)
        extras_found = False
        for name in dirname.rglob("*.mdl"):
            rel_name = name.relative_to(dirname)
            if str(rel_name) not in checksums:
                print("Please remove superfluous file", rel_name)
                extras_found = True
                continue
        if extras_found:
            print("Found superfluous files, please remove and try again.")
            return False

        needs_download = False
        for fname, csum in checksums.items():
            file = dirname / fname
            if not file.is_file():
                print("Missing file", fname)
                needs_download = True
                continue

            if csum != checksum(str(file)):
                print("Checksum mismatch for", fname)
                needs_download = True
                continue
    else:
        needs_download = True

    if not needs_download:
        print("NRPS SVM models present and checked.")
        return True

    print("Downloading NRPS SVM models.")
    check_diskspace(NRPS_SVM_URL)
    download_if_not_present(NRPS_SVM_URL, str(archive_filename), NRPS_SVM_ARCHIVE_CHECKSUM)
    tar = unzip_file(str(archive_filename), lzma, lzma.LZMAError)
    untar_file(tar)
    delete_file(tar)
    assert checksum_file.exists()
    delete_file(str(archive_filename))

    return True


def download_transator(db_dir: str) -> None:
    """Download the Stachelhaus signatures file."""
    version = TRANSATOR_URL.rsplit("_", 1)[1].split(".tar", 1)[0]  # e.g. transator_2022.11.16.tar.xz -> 2022.11.16
    dirname = os.path.join(db_dir, "nrps_pks", "transATor", version)
    archive_path = os.path.join(dirname, TRANSATOR_URL.rpartition("/")[-1])
    profile_path = os.path.join(dirname, "transATor.hmm")

    if present_and_checksum_matches(profile_path, TRANSATOR_CHECKSUM):
        print("TransATor profiles present and checked")
        return

    print("Downloading TransATor profiles.")
    check_diskspace(TRANSATOR_URL)
    download_if_not_present(TRANSATOR_URL, archive_path, TRANSATOR_ARCHIVE_CHECKSUM)
    tar = unzip_file(archive_path, lzma, lzma.LZMAError)
    untar_file(tar)
    delete_file(tar)
    assert os.path.exists(profile_path)
    delete_file(archive_path)


def download(args: argparse.Namespace) -> bool:
    """Download all the large external databases needed, along with non-python
       components.
    """
    download_antismash_js(args.database_dir)

    # grab the latest pfam
    download_pfam(
        args.database_dir,
        PFAM_LATEST_URL,
        PFAM_LATEST_VERSION,
        PFAM_LATEST_ARCHIVE_CHECKSUM,
        PFAM_LATEST_CHECKSUM,
    )

    download_mite(args.database_dir)
    download_resfam(args.database_dir)

    download_tigrfam(args.database_dir)

    download_clusterblast(args.database_dir)
    download_knownclusterblast(args.database_dir)

    for name, details in CLUSTERCOMPARE_DBS.items():
        download_clustercompare(args.database_dir, name, **details)
    for name, details in COMPARIPPSON_DBS.items():
        download_comparippson_db(args.database_dir, name, **details)

    download_stachelhaus(args.database_dir)
    if not download_nprs_svm(args.database_dir):
        return False
    download_transator(args.database_dir)
    return True


def _main() -> None:
    """ Downloads, decompresses, and compiles large databases. Also ensures
        antiSMASH's module data is prepared.
    """
    # Small dance to grab the antiSMASH config for the database dir.
    # All the modules are required to parse the config file,
    # and any executable paths defined should be kept.
    all_modules = antismash.get_detection_modules() + antismash.get_analysis_modules()
    config = antismash.config.build_config(args=[], parser=None, isolated=False, modules=all_modules)

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--database-dir",
        default=config.database_dir,
        metavar="DIR",
        help="Base directory for the antiSMASH databases (default: %(default)s).",
    )

    args = parser.parse_args()
    antismash.config.update_config({"database_dir": args.database_dir})
    if not download(args):
        print("Errors occurred while downloading, aborted.")
        sys.exit(1)
    try:
        print("Pre-building all databases...")
        antismash.main.prepare_module_data()
        print("done.")
    except Exception as err:  # pylint: disable=broad-except
        print("Error encountered while preparing module data:", str(err), file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    _main()
