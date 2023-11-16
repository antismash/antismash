"""Setuptools magic to install antiSMASH."""
import glob
import os
import subprocess
import sys

from setuptools import setup, find_packages
from setuptools.command.test import test as TestCommand


def read(fname):
    """Read a file from the current directory."""
    return open(os.path.join(os.path.dirname(__file__), fname), encoding="utf-8").read()


long_description = read('README.md')

install_requires = [
    'brawn == 1.0.1',
    'numpy == 1.26.2',
    'biopython == 1.78',
    'helperlibs == 0.2.1',
    'jinja2 == 3.1.2',
    'joblib == 1.3.2',
    'jsonschema==4.11.0',
    'markupsafe == 2.1.3',
    'nrpys == 0.1.1',
    'pysvg-py3 == 0.2.2.post3',
    'bcbio-gff == 0.7.0',
    'libsass == 0.22.0',
    'matplotlib == 3.8.1',
    'scipy == 1.11.3',
    'scikit-learn == 1.3.2',
    'MOODS-python == 1.9.4.1'
]

tests_require = [
    'pytest >= 7.2.0, < 8',
    'coverage',
    'pylint == 3.0.2',
    'mypy == 0.982',  # for consistent type checking
]


def read_version():
    """Read the version from the appropriate place in the library."""
    with open(os.path.join('antismash', 'main.py'), 'r', encoding="utf-8") as handle:
        for line in handle:
            if line.startswith('__version__'):
                return line.split('=')[-1].strip().strip('"')
    raise ValueError("unable to find version")


def find_data_files():
    """Setuptools package_data globbing is stupid, so make this work ourselves."""
    data_files = []
    for pathname in glob.glob("antismash/**/*", recursive=True):
        if pathname.endswith('.pyc'):
            continue
        if pathname.endswith('.py'):
            continue
        if '__pycache__' in pathname:
            continue
        if pathname[:-1].endswith('.hmm.h3'):
            continue
        if pathname.endswith('bgc_seeds.hmm'):
            continue
        pathname = glob.escape(pathname)
        pathname = pathname[10:]
        data_files.append(pathname)
    if "HARDCODE_ANTISMASH_GIT_VERSION" in os.environ:
        version_file = os.path.join('antismash', 'git_hash')
        with open(version_file, 'wt', encoding="utf-8") as handle:
            try:
                git_version = subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD'],
                                                      universal_newlines=True).strip()
                changes = subprocess.check_output(['git', 'status', '--porcelain'],
                                                  universal_newlines=True).splitlines()
                if len(changes) != 0:
                    git_version += "(changed)"
                handle.write(git_version)
            except (OSError, subprocess.CalledProcessError):
                pass
        data_files.append(version_file)
    return data_files


class PyTest(TestCommand):
    """Allow running tests via python setup.py test."""

    def finalize_options(self):
        """Test command magic."""
        TestCommand.finalize_options(self)
        self.test_args = []
        self.test_suite = True  # pylint: disable=attribute-defined-outside-init

    def run_tests(self):
        """Run tests."""
        import pytest  # pylint: disable=import-outside-toplevel
        errcode = pytest.main(self.test_args)
        sys.exit(errcode)


setup(
    name="antismash",
    python_requires='>=3.9',
    version=read_version(),
    packages=find_packages(exclude=["*.tests", "*.tests.*", "tests.*", "tests"]),
    package_data={
        'antismash': find_data_files(),
    },
    author='antiSMASH development team',
    author_email='antismash@secondarymetabolites.org',
    description='The antibiotics and Secondary Metabolites Analysis Shell.',
    long_description=long_description,
    long_description_content_type='text/markdown',
    install_requires=install_requires,
    tests_require=tests_require,
    entry_points={
        'console_scripts': [
            'download-antismash-databases=antismash.download_databases:_main',
            'antismash=antismash.__main__:entrypoint',
        ],
    },
    cmdclass={'test': PyTest},
    url='https://github.com/antismash/antismash',
    license='GNU Affero General Public License v3 or later (AGPLv3+)',
    classifiers=[
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: GNU Affero General Public License v3 or later (AGPLv3+)',
        'Operating System :: OS Independent',
    ],
    extras_require={
        'testing': tests_require,
    },
)
