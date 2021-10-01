# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""Script to download and install Pfam and ClusterBlast databases."""

import argparse
import gzip
import hashlib
import lzma
import os
import sys
import tarfile
from typing import Any, Type
from urllib import error as urlerror
from urllib import request

import antismash
from antismash.common.hmmer import ensure_database_pressed
from antismash.common.subprocessing import execute

PFAM_LATEST_VERSION = "34.0"
PFAM_LATEST_URL = "ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam{}/Pfam-A.hmm.gz".format(
    PFAM_LATEST_VERSION
)
PFAM_LATEST_ARCHIVE_CHECKSUM = "b18a98bb1f92b9afb7533fda43baac98d9c69e5b97e68592a43fc00d671ea3f6"
PFAM_LATEST_CHECKSUM = "bfc669e58fdafad0950a52bb9020d6e23c87aa174f4f5bc4ebdac8bb9c30b634"

CLUSTERBLAST_URL = "https://dl.secondarymetabolites.org/releases/clusterblast/clusterblast_20201113.tar.xz"
CLUSTERBLAST_ARCHIVE_CHECKSUM = "bf45276f034615b0827627d16426fec2aca0464391bc41c76ccfae81049ba95a"
CLUSTERBLAST_FASTA_CHECKSUM = "ea43624407ae399cd6e78bf2e6868e9ff0b9581f333645063212d0ef239c7f5b"

CLUSTERCOMPARE_DBS = {
    "mibig": {
        "url": "https://dl.secondarymetabolites.org/releases/clustercompare/cc_mibig_20201210.tar.xz",
        "archive": "25bc9ce7d09f8325d93c0345367def76147e0c6c4c8bd5adbd2ecc79a1827170",
        "fasta": "3fe75f4a95823297ea14817da328000f160f15e5e78335c629812a6cc0c64dcf",
    },
}

RESFAM_URL = "http://dantaslab.wustl.edu/resfams/Resfams.hmm.gz"
RESFAM_ARCHIVE_CHECKSUM = "82e9325283b999b1fb1351502b2d12194561c573d9daef3e623e905c1af66fd6"
RESFAM_LINES = 132214
RESFAM_SIZE = 20123005

TIGRFAM_URL = "https://dl.secondarymetabolites.org/releases/tigrfam/TIGRFam.hmm.gz"
TIGRFAM_ARCHIVE_CHECKSUM = "0747a5efa85d594cd598e67a04611328ab865e6a9401eb143fce6e93edf22bd0"
TIGRFAM_CHECKSUM = "b905785603e1015bbe78d0235d0a6a131345bae35aaab7afbd8f7b59f7798842"

LOCAL_FILE_PATH = os.path.abspath(os.path.dirname(__file__))

CHUNK = 128 * 1024


class DownloadError(RuntimeError):
    """Exception to throw when downloads fail."""

    pass  # pylint: disable=unnecessary-pass


def get_remote_filesize(url: str) -> int:
    """Get the file size of the remote file."""
    try:
        usock = request.urlopen(request.Request(url, method="HEAD"))
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
            "ERROR: Insufficient disk space available (required: %d, free: %d)." % (dbfilesize, free_space)
        )


def download_file(url: str, filename: str) -> str:
    """Download a file."""
    try:
        req = request.urlopen(url)
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
                    "\rDownloading {}: {:5.2f}% downloaded.".format(
                        basename, (overall / expected_size) * 100
                    ),
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
            "Error extracting %s. Please try to extract it manually." % (os.path.basename(filename))
        )
    print("Extraction of %s finished successfully." % (os.path.basename(filename)))
    return newfilename


def untar_file(filename: str) -> None:
    """Extract a TAR/GZ file."""
    try:
        tar = tarfile.open(filename)
        tar.extractall(path=filename.rpartition(os.sep)[0])
        tar.close()
    except tarfile.ReadError:
        print(
            "ERROR: Error extracting %s. Please try to extract it manually."
            % (filename.rpartition(os.sep)[2])
        )
        return
    print("Extraction of %s finished successfully." % (filename.rpartition(os.sep)[2]))


def delete_file(filename: str) -> None:
    """Delete a file."""
    try:
        os.remove(filename)
    except OSError:
        pass


def present_and_checksum_matches(filename: str, sha256sum: str) -> bool:
    """Check if a file is present and the checksum matches."""
    if os.path.exists(filename):
        print("Creating checksum of %s" % os.path.basename(filename))
        csum = checksum(filename)
        if csum == sha256sum:
            return True
    return False


def present_and_line_count_matches(filename: str, lines: int) -> bool:
    """Check if a file is present and the number of lines matches."""
    if not os.path.exists(filename):
        return False

    _count = -1
    for _count, _ in enumerate(open(filename, "r")):
        pass
    _count += 1

    return _count == lines


def present_and_size_matches(filename: str, size: int) -> bool:
    """Check if a file is present and the file size matches."""
    if not os.path.exists(filename):
        return False

    read_size = os.stat(filename).st_size
    return size == read_size


def download_if_not_present(url: str, filename: str, sha256sum: str) -> None:
    """Download a file if it's not present or checksum doesn't match."""
    # If we are missing the archive file, go and download
    if not present_and_checksum_matches(filename, sha256sum):
        download_file(url, filename)

    print("Creating checksum of %s" % os.path.basename(filename))
    csum = checksum(filename)
    if csum != sha256sum:
        raise DownloadError(
            "Error downloading %s, sha256sum mismatch. Expected %s, got %s." % (filename, sha256sum, csum)
        )


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
    url = RESFAM_URL

    # checksum of existing not matched because it has a convert timestamp in it
    # So check size and line count as an approximation
    if present_and_size_matches(filename, RESFAM_SIZE) and \
       present_and_line_count_matches(filename, RESFAM_LINES):
        print("Resfams database present and checked")
        return

    print("Downloading Resfam database")
    check_diskspace(url)
    download_if_not_present(url, archive_filename, RESFAM_ARCHIVE_CHECKSUM)
    filename = unzip_file(archive_filename, gzip, gzip.zlib.error)  # type: ignore
    delete_file(filename + ".gz")
    # remove tabs
    converted = execute(["hmmconvert", filename])
    print("Ensuring all cutoffs are present")
    # add TC to those entries missing them
    # calculated as 10% less than the minimum scoring hit in their own group
    missing_cutoffs = {
        "RF0174": int(374 * 0.9),
        "RF0172": int(85 * 0.9),
        "RF0173": int(295 * 0.9),
        "RF0168": int(691 * 0.9),
    }
    with open(filename, "w") as handle:
        lines = list(converted.stdout)
        i = 0
        while i < len(lines):
            # find an accession
            while i < len(lines) and not lines[i].startswith("ACC"):
                handle.write(lines[i])
                i += 1
            # end of file with no new accession
            if i >= len(lines):
                break
            # write the accession line itself
            handle.write(lines[i])

            # add the cutoffs if missing
            acc = lines[i].split()[1]
            if acc not in missing_cutoffs:
                continue
            value = missing_cutoffs[acc]
            # an accession of interest, so add cutoffs in the same place as others
            while not lines[i].startswith("CKSUM"):
                handle.write(lines[i])
                i += 1
            # write the CKSUM line
            handle.write(lines[i])
            # and finally add the cutoffs
            for cutoff in ["GA", "TC", "NC"]:
                handle.write("%s    %d.00 %d.00\n" % (cutoff, value, value))
            i += 1

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
    archive_filename = os.path.join(db_dir, "clustercompare", url.rpartition("/")[2])
    fasta_filename = os.path.join(db_dir, "clustercompare", name, "proteins.fasta")

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


def download(args: argparse.Namespace) -> None:
    """Download all the large external databases needed."""
    # grab the latest pfam
    download_pfam(
        args.database_dir,
        PFAM_LATEST_URL,
        PFAM_LATEST_VERSION,
        PFAM_LATEST_ARCHIVE_CHECKSUM,
        PFAM_LATEST_CHECKSUM,
    )

    download_resfam(args.database_dir)

    download_tigrfam(args.database_dir)

    download_clusterblast(args.database_dir)
    for name, details in CLUSTERCOMPARE_DBS.items():
        download_clustercompare(args.database_dir, name, **details)


def _main() -> None:
    """ Downloads, decompresses, and compiles large databases. Also ensures
        antiSMASH's module data is prepared.
    """
    # Small dance to grab the antiSMASH config for the database dir.
    # We don't actually want to keep anything else, but we need to load all the
    # modules to make sure we can parse the file.
    all_modules = antismash.get_detection_modules() + antismash.get_analysis_modules()
    config = antismash.config.build_config(args=[], parser=None, isolated=True, modules=all_modules)

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--database-dir",
        default=config.database_dir,
        metavar="DIR",
        help="Base directory for the antiSMASH databases (default: %(default)s).",
    )

    args = parser.parse_args()
    antismash.config.update_config({"database_dir": args.database_dir})
    download(args)
    try:
        print("Pre-building all databases...")
        antismash.main.prepare_module_data()
        print("done.")
    except Exception as err:  # pylint: disable=broad-except
        print("Error encountered while preparing module data:", str(err), file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    _main()
