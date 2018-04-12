# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""Script to download and install Pfam and ClusterBlast databases."""

from urllib import request, error as urlerror
import tarfile
import gzip
import hashlib
import lzma
import os
import platform
import ctypes

from antismash.common.subprocessing import execute


PFAM27_URL = "ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam27.0/Pfam-A.hmm.gz"
PFAM27_ARCHIVE_CHECKSUM = "b29bc2c54db8090531df0361a781b8d7397f60ebedc0c36a16e7d45e999cc329"
PFAM27_CHECKSUM = "ea35d7e4029b9d34eb8422ae69a230a2a5d25a52a8e207425791e8fdeb38aac8"

PFAM_LATEST_VERSION = "31.0"
PFAM_LATEST_URL = "ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/" \
                  "Pfam{}/Pfam-A.hmm.gz".format(PFAM_LATEST_VERSION)
PFAM_LATEST_ARCHIVE_CHECKSUM = "33f2dfab4d695fc3d2337119debbfcd8c801aaf8e2c312bd738c105a84007973"
PFAM_LATEST_CHECKSUM = "40b742a1bebca10f6b86c1364624dbdc34f79890954777eb95f8a48428de39d3"

CLUSTERBLAST_URL = "https://dl.secondarymetabolites.org/releases/4.0.0/clusterblast_20170105_v8_31.tar.xz"
CLUSTERBLAST_ARCHIVE_CHECKSUM = "c9c1f2c07ce97ea453345034727555928e7b0b8907f700805bb6a912865bb315"
CLUSTERBLAST_DMND_CHECKSUM = "388df3e711b3049ad851bfc8bd45ec292a3808907f048e6a7e5f4a25b90699f8"

LOCAL_FILE_PATH = os.path.abspath(os.path.dirname(__file__))

CHUNK = 128 * 1024


class DownloadError(RuntimeError):
    """Exception to throw when downloads fail."""

    pass


def get_remote_filesize(url):
    """Get the file size of the remote file."""
    try:
        usock = request.urlopen(request.Request(url, method='HEAD'))
        dbfilesize = usock.info().get('Content-Length', '0')
    except urlerror.URLError:
        dbfilesize = '0'

    dbfilesize = int(dbfilesize)  # db file size in bytes
    return dbfilesize


def get_free_space(folder):
    """Return folder/drive free space (in bytes)."""
    if platform.system() == 'Windows':
        free_bytes = ctypes.c_ulonglong(0)
        ctypes.windll.kernel32.GetDiskFreeSpaceExW(ctypes.c_wchar_p(folder), None, None, ctypes.pointer(free_bytes))
        return free_bytes.value
    return os.statvfs(folder).f_bfree * os.statvfs(folder).f_frsize


def check_diskspace(file_url):
    """Check if sufficient disk space is available."""
    dbfilesize = get_remote_filesize(file_url)
    free_space = get_free_space(".")
    if free_space < dbfilesize:
        raise DownloadError('ERROR: Insufficient disk space available (required: %d, free: %d).' %
                            (dbfilesize, free_space))


def download_file(url, filename):
    """Download a file."""
    try:
        req = request.urlopen(url)
    except urlerror.URLError:
        raise DownloadError('ERROR: File not found on server.\nPlease check your internet connection.')

    # use 1 because we want to divide by the expected size, can't use 0
    expected_size = int(req.info().get('Content-Length', '1'))

    basename = os.path.basename(filename)
    dirname = os.path.dirname(filename)
    if not os.path.isdir(dirname):
        os.makedirs(dirname)

    overall = 0
    with open(filename, 'wb') as fp:
        while True:
            try:
                chunk = req.read(CHUNK)
                if not chunk:
                    print('')
                    break
                overall += len(chunk)
                print("\rDownloading {}: {:5.2f}% downloaded.".format(
                      basename, (overall / expected_size) * 100), end='')
                fp.write(chunk)
            except IOError:
                raise DownloadError('ERROR: Download interrupted.')
    return filename


def checksum(filename, chunksize=2 ** 20):
    """Get the SHA256 checksum of a file."""
    sha = hashlib.sha256()
    with open(filename, 'rb') as fh:
        for chunk in iter(lambda: fh.read(chunksize), b''):
            sha.update(chunk)

    return sha.hexdigest()


def unzip_file(filename, decompressor, error_type):
    """Decompress a compressed file."""
    newfilename, _ = os.path.splitext(filename)
    try:
        zipfile = decompressor.open(filename, 'rb')
        chunksize = 128 * 1024
        with open(newfilename, 'wb') as fp:
            while True:
                try:
                    chunk = zipfile.read(chunksize)
                    if not chunk:
                        break
                    fp.write(chunk)
                except IOError:
                    raise DownloadError('ERROR: Unzipping interrupted.')
    except error_type:
        print("ERROR: Error extracting %s. Please try to extract it manually." %
              (os.path.basename(filename)))
        return
    print("Extraction of %s finished successfully." % (os.path.basename(filename)))
    return newfilename


def untar_file(filename):
    """Extract a TAR/GZ file."""
    try:
        tar = tarfile.open(filename)
        tar.extractall(path=filename.rpartition(os.sep)[0])
        tar.close()
    except tarfile.ReadError:
        print("ERROR: Error extracting %s. Please try to extract it manually." %
              (filename.rpartition(os.sep)[2]))
        return
    print("Extraction of %s finished successfully." % (filename.rpartition(os.sep)[2]))


# TODO: use common function?
def compile_pfam(filename):
    """Compile a HMMer database with hmmpress."""
    command = ['hmmpress', '-f', filename]
    execute(command)


def delete_file(filename):
    """Delete a file."""
    try:
        os.remove(filename)
    except OSError:
        pass


def present_and_checksum_matches(filename, sha256sum):
    """Check if a file is present and the checksum matches."""
    if os.path.exists(filename):
        print("Creating checksum of %s" % os.path.basename(filename))
        csum = checksum(filename)
        if csum == sha256sum:
            return True
    return False


def download_if_not_present(url, filename, sha256sum):
    """Download a file if it's not present or checksum doesn't match."""
    # If we are missing the archive file, go and download
    if not present_and_checksum_matches(filename, sha256sum):
        download_file(url, filename)

    print("Creating checksum of %s" % os.path.basename(filename))
    csum = checksum(filename)
    if csum != sha256sum:
        raise DownloadError("Error downloading %s, sha256sum mismatch. Expected %s, got %s." %
                            (filename, sha256sum, csum))


def download_pfam(url, version, archive_checksum, db_checksum):
    """Download and compile the PFAM database."""
    archive_filename = os.path.join(LOCAL_FILE_PATH, "databases", "pfam", version, "Pfam-A.hmm.gz")
    db_filename = os.path.splitext(archive_filename)[0]

    if present_and_checksum_matches(db_filename, db_checksum):
        print("PFAM file present and ok for version", version)
        return

    print("Downloading PFAM version", version)
    check_diskspace(url)
    download_if_not_present(url, archive_filename, archive_checksum)
    filename = unzip_file(archive_filename, gzip, gzip.zlib.error)
    compile_pfam(filename)
    delete_file(filename + ".gz")


def download_clusterblast():
    """Download the clusterblast database."""
    database_dir = os.path.join(LOCAL_FILE_PATH, "databases",)
    archive_filename = os.path.join(database_dir, CLUSTERBLAST_URL.rpartition('/')[2])
    dmnd_filename = os.path.join(database_dir, "clusterblast", "geneclusterprots.dmnd")

    if present_and_checksum_matches(dmnd_filename, CLUSTERBLAST_DMND_CHECKSUM):
        print("ClusterBlast databse present and checked")
        return

    print("Downloading ClusterBlast database.")
    check_diskspace(CLUSTERBLAST_URL)
    download_if_not_present(CLUSTERBLAST_URL, archive_filename, CLUSTERBLAST_ARCHIVE_CHECKSUM)
    filename = unzip_file(archive_filename, lzma, lzma.LZMAError)
    untar_file(filename)
    delete_file(filename)
    delete_file(filename + ".xz")


def main():
    """Download and compile all the large external databases needed."""
    # ClusterFinder is stuck to PFAM 27.0, so always grab that
    download_pfam(PFAM27_URL, "27.0", PFAM27_ARCHIVE_CHECKSUM, PFAM27_CHECKSUM)

    # And also grab the latest
    download_pfam(PFAM_LATEST_URL, PFAM_LATEST_VERSION, PFAM_LATEST_ARCHIVE_CHECKSUM, PFAM_LATEST_CHECKSUM)

    download_clusterblast()

    # hmmpress the NRPS/PKS specific databases
    compile_pfam(os.path.join(LOCAL_FILE_PATH, "detection", "nrps_pks_domains", "data", "abmotifs.hmm"))
    compile_pfam(os.path.join(LOCAL_FILE_PATH, "detection", "nrps_pks_domains", "data", "dockingdomains.hmm"))
    compile_pfam(os.path.join(LOCAL_FILE_PATH, "detection", "nrps_pks_domains", "data", "ksdomains.hmm"))
    compile_pfam(os.path.join(LOCAL_FILE_PATH, "detection", "nrps_pks_domains", "data", "nrpspksdomains.hmm"))
    # TODO: re-add a compile call for SANDPUMA once that is in

    # hmmpress the smcog specific database
    compile_pfam(os.path.join(LOCAL_FILE_PATH, "modules", "smcogs", "data", "smcogs.hmm"))

    # TODO: Press the bgc_seeds once there is an antismash.common function for it


if __name__ == "__main__":
    main()
