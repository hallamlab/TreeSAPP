# TreeSAPP: Tree-based Sensitive and Accurate Phylogenetic Profiling

Connor Morgan-Lang, Kishori M. Konwar, Zachary Armstrong, Young C. Song, and Steven J. Hallam

## Overview:

A python pipeline for identifying and mapping marker gene sequences onto reference phylogenetic trees.
 Here is a workflow diagram:

![alt text]()

## Download and installation:

Otherwise, carry on. The first command will clone the TreeSAPP repository
from the GitHub page while the second and third compile C++ and Python extensions
required by the pipeline.

```
git clone git@github.com:hallamlab/TreeSAPP.git
make
make install
```

### Downloading dependencies:

If you do not already have the dependencies for TreeSAPP installed on your computer,
 we've listed how to easily download and install each one below. Good luck!

#### RAxML
A simple `git clone` of their [GitHub page](https://github.com/stamatak/standard-RAxML) should work
for Linux and Mac operating systems. From here, consult the README file in the standard-RAxML directory for
installation instructions using make.
We have tested several versions and found no problems from V.7.1 to the most recent release as of 
December 1st, 2015. However, the executable MUST be named `raxmlHPC` or it will not be found by TreeSAPP!
If you find an incompatibility please notify us through the Issues feed!

#### HMMER
TreeSAPP uses HMMER for identifying marker gene sequences in proteins and genomes.
The latest version (HMMER3.1b2) is available at http://hmmer.org/.
Download it from there and follow their installation guide under DOCUMENTATION.

#### Prodigal
Prodigal (version 2.6.3) can be downloaded from the [GitHub page](https://github.com/hyattpd/Prodigal).
Follow the [installation guide](https://github.com/hyattpd/Prodigal/wiki/installation) on their GitHub wiki to install.
There is an upcoming version 3 so these links may become outdated soon!

#### trimal
The software and installation instructions are available on the
[trimAl GitHub](https://github.com/scapella/trimal).

#### MAFFT
 MAFFT multiple alignment software is only required for creating and updating reference packages,
 it is not a part of the main workflow. Therefore, feel free to skip installing MAFFT unless you plan on
 doing either one of those tasks. If not, here is the [MAFFT webpage](https://mafft.cbrc.jp/alignment/software/).
 Download and installation instructions are available from there.

#### Finishing up
I hope that wasn't too painful. If you think you have installed everything, try running TreeSAPP!
It will check for the required executables up front and you will be
quickly notified if some are missing, or at least TreeSAPP is unable to find them.

### Running TreeSAPP

To list all the options with brief help statements `./treesapp.py -h`.

To perform a basic run with only the required arguments:
```
./treesapp.py -i input.fasta -o ~/path/to/output/directory/
```
Executables are automatically detected in both the $PATH and in the
sub_binaries/mac or sub_binaries/ubuntu, depending on your OS. However, if your executables
are together elsewhere, TreeSAPP can be directed to them with `--executables`.


### Other tutorials of interest

For an easy introduction to TreeSAPP,
please refer to the [Google Cloud Platform tutorial](https://github.com/hallamlab/TreeSAPP/blob/fastsearch/GCP_tutorial.md).

Yet to come:

- [Building reference packages]()
- [Evaluating classification accuracy]()
- [Annotating extra features]()

