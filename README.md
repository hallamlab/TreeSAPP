# TreeSAPP: Tree-based Sensitive and Accurate Phylogenetic Profiler

[![Build Status](https://travis-ci.org/hallamlab/TreeSAPP.svg?branch=master)](https://travis-ci.org/hallamlab/TreeSAPP)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/1ae0c5b20743430497e1397f6e91839f)](https://www.codacy.com/manual/TreeSAPP/TreeSAPP?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=hallamlab/TreeSAPP&amp;utm_campaign=Badge_Grade)
[![PyPI version](https://badge.fury.io/py/treesapp.svg)](https://badge.fury.io/py/treesapp)

[![Anaconda-Server Badge](https://anaconda.org/bioconda/treesapp/badges/installer/conda.svg)](https://conda.anaconda.org/bioconda)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/treesapp/badges/platforms.svg)](https://anaconda.org/bioconda/treesapp)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/treesapp/badges/version.svg)](https://anaconda.org/bioconda/treesapp)

Connor Morgan-Lang, Ryan McLaughlin, Grace Zhang, Kevin Chan, Zachary Armstrong, and Steven J. Hallam

## Overview:

TreeSAPP is a python package for phylogenetically annotating genomes and metagenomes.
 Here is a diagram of the workflow:

![alt text](https://github.com/hallamlab/TreeSAPP/blob/master/dev_utils/pipeline_figure_horizontal.png)

## Installation:

TreeSAPP supports Python versions 3.5, 3.6, 3.7 and 3.8.

### Conda
TreeSAPP and most of its dependencies can be installed in its own environment using conda.

```bash
conda create -n treesapp_cenv -c bioconda -c conda-forge treesapp
conda activate treesapp_cenv
```
If you plan on building your own reference packages you will also require [USEARCH](https://www.drive5.com/usearch/).

### Singularity
If you're working in an HPC environment and don't have conda installed, we also have a 
[singularity](https://sylabs.io/guides/3.5/user-guide/) container available:
```bash
singularity pull library://cmorganl/default/treesapp
singularity exec treesapp.sif
```

### PyPI
The most recent version of TreeSAPP is hosted on the Python Package Index (PyPI) and can be installed using `pip install treesapp`.
 Alternatively you can install the latest development version of TreeSAPP locally with `git clone`.
 In either case we recommend installing within a virtual environment using the python package [`virtualenv`](https://virtualenv.pypa.io/en/latest/).
```
cd ~/bin
virtualenv ~/bin/treesapp_venv
source ~/bin/treesapp_venv/bin/activate
git clone https://github.com/hallamlab/TreeSAPP.git
cd TreeSAPP/
python setup.py sdist
pip install dist/treesapp*.tar.gz
```

If you opted to install TreeSAPP either using `pip` or by cloning the development version from GitHub you will need to
[install dependencies](https://github.com/hallamlab/TreeSAPP/blob/master/docs/dep_install.md) that you do not already
have installed (i.e. they will need to be in you're environment's path).

## Running TreeSAPP

To list all the sub-commands run `treesapp`.

To test the `assign` workflow, run:
```
treesapp assign -i ~/bin/TreeSAPP/test_data/marker_test_suite.faa -m prot --trim_align -o assign_test -t M0701,M0702,M0705
```

To assign sequences in your genome of interest:
```
treesapp assign -i Any.fasta -o ~/path/to/output/directory/
```
As in the previous command, we recommend using the `--trim_align` flag and increasing the number of threads to use with `-n`.


## Tutorials

If we do not yet have a reference package for a gene you are interested in,
please try [building a new reference package](https://github.com/hallamlab/TreeSAPP/blob/master/docs/Marker_package_creation_tutorial.md).
Of course, if you run into any problems or would like to collaborate on building many reference packages
don't hesitate to email us or create a new issue with an 'enhancement' label.

To determine whether the sequences used to build your new reference package are what you think they are,
 and whether it might unexpectedly annotate homologous sequences,
 see the [purity tutorial](https://github.com/hallamlab/TreeSAPP/blob/master/docs/purity_tutorial.md).

If you are working with a particularly complex reference package, from an orthologous group for example, or have extra
 phylogenetic information you'd like to include in your classifications,
 try [annotating extra features](https://github.com/hallamlab/TreeSAPP/blob/master/docs/layering_tutorial.md) with `treesapp layer`.


## Yet to come:
- [Interpreting `treesapp assign` results]()
- [Evaluating classification accuracy]()
- [Taxonomically decorating trees for iTOL]()
- [Terraform](https://github.com/hallamlab/TreeSAPP/tree/master/terraform)

