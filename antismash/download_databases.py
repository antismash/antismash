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

PFAM27_URL = "ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam27.0/Pfam-A.hmm.gz"
PFAM27_ARCHIVE_CHECKSUM = "b29bc2c54db8090531df0361a781b8d7397f60ebedc0c36a16e7d45e999cc329"
PFAM27_CHECKSUM = "ea35d7e4029b9d34eb8422ae69a230a2a5d25a52a8e207425791e8fdeb38aac8"

PFAM_LATEST_VERSION = "31.0"
PFAM_LATEST_URL = "ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam{}/Pfam-A.hmm.gz".format(
    PFAM_LATEST_VERSION
)
PFAM_LATEST_ARCHIVE_CHECKSUM = "33f2dfab4d695fc3d2337119debbfcd8c801aaf8e2c312bd738c105a84007973"
PFAM_LATEST_CHECKSUM = "40b742a1bebca10f6b86c1364624dbdc34f79890954777eb95f8a48428de39d3"

CLUSTERBLAST_URL = "https://dl.secondarymetabolites.org/releases/4.0.0/clusterblast_20170105_v8_31.tar.xz"
CLUSTERBLAST_ARCHIVE_CHECKSUM = "c9c1f2c07ce97ea453345034727555928e7b0b8907f700805bb6a912865bb315"
CLUSTERBLAST_DMND_CHECKSUM = "388df3e711b3049ad851bfc8bd45ec292a3808907f048e6a7e5f4a25b90699f8"

RESFAM_URL = "http://dantaslab.wustl.edu/resfams/Resfams.hmm.gz"
RESFAM_ARCHIVE_CHECKSUM = "82e9325283b999b1fb1351502b2d12194561c573d9daef3e623e905c1af66fd6"
RESFAM_LINES = 132214
RESFAM_SIZE = 20123005

LOCAL_FILE_PATH = os.path.abspath(os.path.dirname(__file__))

CHUNK = 128 * 1024


class DownloadError(RuntimeError):
    """Exception to throw when downloads fail."""

    pass


def get_remote_filesize(url: str) -> int:
    """Get the file size of the remote file."""
    try:
        usock = request.urlopen(request.Request(url, method="HEAD"))
        dbfilesize = usock.info().get("Content-Length", "0")  # type: ignore
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


def download_clusterblast(db_dir: str) -> None:
    """Download the clusterblast database."""
    archive_filename = os.path.join(db_dir, CLUSTERBLAST_URL.rpartition("/")[2])
    dmnd_filename = os.path.join(db_dir, "clusterblast", "geneclusterprots.dmnd")

    if present_and_checksum_matches(dmnd_filename, CLUSTERBLAST_DMND_CHECKSUM):
        print("ClusterBlast database present and checked")
        return

    print("Downloading ClusterBlast database.")
    check_diskspace(CLUSTERBLAST_URL)
    download_if_not_present(CLUSTERBLAST_URL, archive_filename, CLUSTERBLAST_ARCHIVE_CHECKSUM)
    filename = unzip_file(archive_filename, lzma, lzma.LZMAError)
    untar_file(filename)
    delete_file(filename)
    delete_file(filename + ".xz")


def download(args: argparse.Namespace) -> None:
    """Download all the large external databases needed."""
    # ClusterFinder is stuck to PFAM 27.0, so always grab that
    download_pfam(args.database_dir, PFAM27_URL, "27.0", PFAM27_ARCHIVE_CHECKSUM, PFAM27_CHECKSUM)

    # And also grab the latest
    download_pfam(
        args.database_dir,
        PFAM_LATEST_URL,
        PFAM_LATEST_VERSION,
        PFAM_LATEST_ARCHIVE_CHECKSUM,
        PFAM_LATEST_CHECKSUM,
    )

    download_resfam(args.database_dir)

    download_clusterblast(args.database_dir)


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
    download(args)
    try:
        antismash.main.prepare_module_data()
    except Exception as err:  # pylint: disable=broad-except
        print("Error encountered while preparing module data:", str(err), file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    _main()
