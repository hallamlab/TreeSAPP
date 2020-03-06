# TreeSAPP: Tree-based Sensitive and Accurate Phylogenetic Profiler

[![Build Status](https://travis-ci.org/hallamlab/TreeSAPP.svg?branch=master)](https://travis-ci.org/hallamlab/TreeSAPP)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/1ae0c5b20743430497e1397f6e91839f)](https://www.codacy.com/manual/TreeSAPP/TreeSAPP?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=hallamlab/TreeSAPP&amp;utm_campaign=Badge_Grade)

Connor Morgan-Lang, Ryan McLaughlin, Grace Zhang, Kevin Chan, Zachary Armstrong, and Steven J. Hallam

## Overview:

TreeSAPP is a python package for phylogenetically annotating genomes and metagenomes.
 Here is diagram of the workflow:

![alt text](https://github.com/hallamlab/TreeSAPP/blob/master/dev_utils/pipeline_figure_horizontal.png)

## Download and installation:

TreeSAPP can be installed using conda:

```bash
conda install -c conda-forge -c bioconda treeesapp
```
If you're working in an HPC environment and don't have conda installed, we also have a 
[singularity](https://sylabs.io/guides/3.5/user-guide/) container available:

```bash
singularity pull library://cmorganl/treesapp
singularity exec treesapp.sif
```

Finally, if you want to install the latest version of TreeSAPP locally,
 you can use `git clone` to pull down the latest version.
We recommend using a virtual environment using the python package [`virtualenv`](https://virtualenv.pypa.io/en/latest/) 
while installing TreeSAPP and all dependencies.

```
cd ~/bin
virtualenv ~/bin/treesapp_venv
source ~/bin/treesapp_venv/bin/activate
pip install treesapp
make rpkm
```

However, the pipeline will not run without several dependencies.

### Downloading dependencies:

If you do not already have the dependencies for TreeSAPP installed on your computer,
 we've listed how to easily download and install each one below. Good luck!

#### RAxML
A simple `git clone` of their [GitHub page](https://github.com/stamatak/standard-RAxML) should work
for Linux and Mac operating systems.
From here, consult the README file in the standard-RAxML directory for installation instructions using make. 
We highly recommend *only using release 8.2.12* as older versions were found to not estimate pendant distances of placements as accurately.
However, the executable MUST be named `raxmlHPC` or it will not be found by TreeSAPP!

#### HMMER
TreeSAPP uses HMMER for identifying marker gene sequences in proteins and genomes.
The latest version (v3.3) is available at http://hmmer.org/.
Download it from there and follow their installation guide under DOCUMENTATION.

#### Prodigal
Prodigal (version 2.6.3) can be downloaded from the [GitHub page](https://github.com/hyattpd/Prodigal).
Follow the [installation guide](https://github.com/hyattpd/Prodigal/wiki/installation) on their GitHub wiki to install.
There is an upcoming version 3 so these links may become outdated soon!

#### MAFFT
 MAFFT multiple alignment software is only required for creating and updating reference packages,
 it is not a part of the main workflow. Therefore, feel free to skip installing MAFFT unless you plan on
 doing either one of those tasks. If not, here is the [MAFFT webpage](https://mafft.cbrc.jp/alignment/software/).
 Download and installation instructions are available from there.

#### USEARCH
 The current version of TreeSAPP uses USEARCH for multiple clustering stages,
 such as when building reference packages. Fortunately, no application should require
 huge amounts of RAM so we can use the free, 32-bit version available at the
 [Robert Edgar's drive5 website](https://drive5.com/usearch/download.html).

#### OD-Seq
 [OD-Seq](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-015-0702-1) 
 is used for detecting mis-annotated or "outliers" in multiple sequence alignments when building new reference packages.
  It can be installed into TreeSAPP as a part of `make all` or in isolation with `make odseq`.
  Source files can also be downloaded from the University College Dublin's website using this
   [link](http://www.bioinf.ucd.ie/download/od-seq.tar.gz).

#### Finishing up
I hope that wasn't too painful. If you think you have installed everything, try running `treesapp info`!
It will check for the required executables up front and you will be
quickly notified if some are missing, or at least TreeSAPP is unable to find them.
In the case you do not have sudo permissions to move these executables to a globally-available directory (e.g. /usr/local/bin/),
you can copy them to `treesapp_venv/lib/python*/site-packages/treesapp-*.egg/treesapp/sub_binaries/` and TreeSAPP will be able to find and use them.

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
though, as in the previous assign command, we recommend using the `--trim_align` flag,
 and increasing the number of threads and processors to use with `-n`.


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

The easiest way to get started with TreeSAPP is by using [Terraform](https://github.com/hallamlab/TreeSAPP/tree/master/terraform)
 to provision a Google Cloud Platform instance with TreeSAPP and all its dependencies.
 __This is outdated and scripts supporting this are under repair__

Yet to come:
- [Evaluating classification accuracy]()
- [Taxonomically decorating trees for iTOL]()

