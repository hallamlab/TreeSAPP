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
    "Programming Language :: Python :: 3.5",
    "Programming Language :: Python :: 3.6",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
]

SETUP_METADATA = \
    {
        "name": "treesapp",
        "version": "0.2.6",
        "description": "TreeSAPP is a functional and taxonomic annotation tool",
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
                                  sources=["include/tree_parsermodule.cpp"],
                                  language="c++"),
                        Extension("_fasta_reader",
                                  sources=["include/fasta_reader.cpp"],
                                  depends=["include/fasta_reader.hpp"],
                                  language="c++",
                                  include_dirs=["./include/"])
                        ],
        "install_requires": ["pygtrie>=2.3", "ete3", "numpy", "biopython>=1.68", "scipy", "six"]
    }

setup(**SETUP_METADATA)
