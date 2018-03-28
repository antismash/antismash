"""Setuptools magic to install antiSMASH."""
import os
from setuptools import setup
from setuptools.command.test import test as TestCommand
import sys


def read(fname):
    """Read a file from the current directory."""
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


if os.path.exists('README.rst'):
    long_description = read('README.rst')
else:
    long_description = read('README.md')

install_requires = [
    'argparse',
    'cssselect',
    'numpy',
    'biopython >= 1.70',
    'helperlibs',
    'pysvg-py3',
    'pyExcelerator',
    'bcbio-gff',
    'networkx',
    'pandas',
    'matplotlib',
    'scipy',
    'scikit-learn',
]

tests_require = [
    'pytest',
    'minimock',
    'coverage',
    'pylint',
]


def read_version():
    """Read the version fromt he appropriate place in the library."""
    for line in open(os.path.join('antismash', 'main.py'), 'r'):
        if line.startswith('__version__'):
            return line.split('=')[-1].strip().strip('"')


class PyTest(TestCommand):
    """Allow running tests via python setup.py test."""

    def finalize_options(self):
        """Test command magic."""
        TestCommand.finalize_options(self)
        self.test_args = []
        self.test_suite = True

    def run_tests(self):
        """Run tests."""
        import pytest
        errcode = pytest.main(self.test_args)
        sys.exit(errcode)


setup(
    name="antismash",
    version=read_version(),
    packages=['antismash'],
    author='antiSMASH development team',
    author_email='antismash@secondarymetabolites.org',
    description='The antibiotics and Secondary Metabolites Analysis Shell.',
    long_description=long_description,
    install_requires=install_requires,
    tests_require=tests_require,
    entry_points={
        'console_scripts': [
            'download-antismash-databases=antismash.download_databases:main',
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
