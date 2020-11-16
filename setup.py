from setuptools import Extension
from setuptools import setup, find_packages

with open("README.md", "r") as readme:
    LONG_DESCRIPTION = readme.read()

CLASSIFIERS = [
    "Environment :: Console",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
    "Natural Language :: English",
    "Operating System :: POSIX :: Linux",
    "Operating System :: MacOS :: MacOS X",
    "Programming Language :: C++",
    "Programming Language :: Python :: 3.6",
    "Programming Language :: Python :: 3.7",
    "Programming Language :: Python :: 3.8",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
    "Development Status :: 5 - Production/Stable"
]

SETUP_METADATA = \
    {
        "name": "treesapp",
        "version": "0.9.5",
        "description": "TreeSAPP is a functional and taxonomic annotation tool for genomes and metagenomes.",
        "long_description": LONG_DESCRIPTION,
        "long_description_content_type": "text/markdown",
        "author": "Connor Morgan-Lang",
        "author_email": "c.morganlang@gmail.com",
        "url": "https://github.com/hallamlab/TreeSAPP",
        "license": "GPL-3.0",
        "packages": find_packages(),
        "include_package_data": True,
        "entry_points": {'console_scripts': ['treesapp = treesapp.__main__:main']},
        "classifiers": CLASSIFIERS,
        "ext_modules": [Extension("_tree_parser",
                                  sources=["treesapp/extensions/tree_parsermodule.cpp"],
                                  language="c++"),
                        Extension("_fasta_reader",
                                  sources=["treesapp/extensions/fasta_reader.cpp"],
                                  depends=["treesapp/include/fasta_reader.hpp"],
                                  language="c++",
                                  include_dirs=["./treesapp/include"])
                        ],
        "install_requires": ["samsum", "pyfastx",
                             "pygtrie", "six",
                             "biopython", "ete3",
                             "numpy", "scipy", "scikit-learn", "joblib",
                             "seaborn", "matplotlib", "tqdm", "packaging"],
        "setup_requires": [
            "setuptools>=50.3.1"
        ],
        "extras_require": {
            'test': ['pytest', 'pytest-cov'],
        }
    }

setup(**SETUP_METADATA)
