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
                     "version": "0.2.2",
                     "description": "TreeSAPP is a functional and taxonomic annotation tool",
                     "long_description": LONG_DESCRIPTION,
                     "long_description_content_type": "text/markdown",
                     "author": "Connor Morgan-Lang",
                     "author_email": "c.morganlang@gmail.com",
                     "url": "https://github.com/hallamlab/TreeSAPP",
                     "license": "GPL-3.0",
                     "packages": find_packages(),
                     "classifiers": CLASSIFIERS,
                     "ext_modules": [Extension("_tree_parser",
                                               sources=["sub_binaries/TreeSAPP_extensions/tree_parsermodule.cpp"],
                                               language="c++",
                                               include_dirs=["./sub_binaries/TreeSAPP_extensions"]),
                                     Extension("_fasta_reader",
                                               sources=["sub_binaries/TreeSAPP_extensions/fasta_reader.cpp"],
                                               language="c++",
                                               include_dirs=["./sub_binaries/TreeSAPP_extensions"])
                                     ],
                     "install_requires": ["pygtrie>=2.3", "ete3", "numpy", "biopython>=1.68", "scipy"]
               }

setup(**SETUP_METADATA)
