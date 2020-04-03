# TreeSAPP dependencies

If you do not already have the dependencies for TreeSAPP installed on your computer,
 we've listed how to download and install each one below.
 
TreeSAPP can find these if they are already installed somewhere in your environment's PATH (e.g. /usr/local/bin/).
 Your path variable can be printed by typing `echo $PATH` on Linux and Mac systems.
 
Alternatively, if you want to keep these versions separate from those in your PATH you can also store these in your TreeSAPP installation's
 `sub_binaries/` directory. If you cloned TreeSAPP into your `~/bin/` directory,
 then copy all dependency executables into `~/bin/TreeSAPP/sub_binaries/`.
 Things a little more complicated if TreeSAPP was installed via `pip`.
 First, identify the location of your `treesapp` executable with `which treesapp`.
 An example location is `~/treesapp_venv/bin/treesapp`.
 Relative to here, `sub_binaries/` is in `~/treesapp_venv/lib/python3.7/site-packages/treesapp/sub_binaries/`.
 Your path will be different if you are using a Python 3 version other than 3.7.

## Required

### RAxML
A simple `git clone` of their [GitHub page](https://github.com/stamatak/standard-RAxML) should work
for Linux and Mac operating systems.
From here, consult the README file in the standard-RAxML directory for installation instructions using make. 
We highly recommend *only using release 8.2.12* as older versions were found to not estimate pendant distances of placements as accurately.
However, the executable MUST be named `raxmlHPC` or it will not be found by TreeSAPP!

**Recommended version**: 8.2.12

**Citation**: Stamatakis, A. (2006). RAxML-VI-HPC: Maximum likelihood-based phylogenetic analyses with thousands of taxa and mixed models.
 Bioinformatics, 22(21), 2688–2690. https://doi.org/10.1093/bioinformatics/btl446

### HMMER
TreeSAPP uses HMMER for identifying marker gene sequences in proteins and genomes
 and performing a profile alignment to include query sequences in a multiple sequence alignment prior to phylogenetic placement.
The latest version is available at http://hmmer.org/.
Download it from there and follow their installation guide under DOCUMENTATION.

**Recommended version**: 3.3

**Citation**: Eddy, S. R. (1998). Profile hidden Markov models.
 Bioinformatics (Oxford, England), 14(9), 755–763. https://doi.org/btb114

### Prodigal
Prodigal is used for ORF prediction and can be downloaded from the [GitHub page](https://github.com/hyattpd/Prodigal).
Follow the [installation guide](https://github.com/hyattpd/Prodigal/wiki/installation) on their GitHub wiki to install.

**Recommended version**: 2.6.3

**Citation**: Hyatt, D., Chen, G.-L., Locascio, P. F., Land, M. L., Larimer, F. W., & Hauser, L. J. (2010). Prodigal: prokaryotic gene recognition and translation initiation site identification.
 BMC Bioinformatics, 11, 119. https://doi.org/10.1186/1471-2105-11-119

### BMGE
BMGE is used for selecting phylogenetically informative regions from multiple sequence alignments.
 This is optionally used prior to building reference phylogenies and phylogenetic placement due to significant
 reductions in compute time. It is distributed with TreeSAPP (in `treesapp/sub_binaries/`) 
 but can also be installed using [conda](https://anaconda.org/phylofisher/bmge).
 Old download links are no longer functional.

**Recommended version**: 1.12

**Citation**: Criscuolo, A., & Gribaldo, S. (2010). BMGE (Block Mapping and Gathering with Entropy): A new software for selection of phylogenetic informative regions from multiple sequence alignments.
 BMC Evolutionary Biology, 10(1). https://doi.org/10.1186/1471-2148-10-210

## Optional

If you would like to build or update reference packages you will also need to install `FastTree`, `MAFFT`, `USEARCH` and `OD-Seq`.

You will need to install `BWA` if you have FASTQ files that you would like to derive relative abundance values of classified sequences from.

### FastTree

FastTree can be used to build reference trees instead of RAxML by invoking "fast-mode" with the flag '--fast' in `treesapp create`.
 In practice, we haven't observed a drastic decrease in classification performance between RAxML and FastTree so its completely okay to use it in our opinion.
 FastTree can be installed using [conda](https://anaconda.org/bioconda/fasttree) or 
by following the installation instructions at http://microbesonline.org/fasttree/#Install.

**Recommended version**: 2.1.10 Double precision

**Citation**: Price, M. N., Dehal, P. S., & Arkin, A. P. (2010). FastTree 2 – Approximately Maximum-Likelihood Trees for Large Alignments.
 PLoS ONE, 5(3), e9490. https://doi.org/10.1371/journal.pone.0009490

### MAFFT
 MAFFT multiple alignment software is only required for creating and updating reference packages
 (`treesapp create` and `treesapp update`, respectively); it is not a part of the classification workflow.
 Therefore, feel free to skip installing MAFFT unless you plan on doing either one of those tasks.
 If not, here is the [MAFFT webpage](https://mafft.cbrc.jp/alignment/software/).
 Download and installation instructions are available from there.

**Recommended version**: 7.407

**Citation**: Katoh, K., & Standley, D. M. (2013). MAFFT multiple sequence alignment software version 7: Improvements in performance and usability.
 Molecular Biology and Evolution, 30(4), 772–780. https://doi.org/10.1093/molbev/mst010

### USEARCH
 TreeSAPP uses USEARCH for clustering sequences while building reference packages.
 Fortunately, no application should require huge amounts of RAM so we can use the free, 32-bit version available at
 [Robert Edgar's USEARCH website](https://drive5.com/usearch/download.html).
 You will receive a compiled binary file for your distribution.

**Recommended version**: 11.0.667

**Citation**: Edgar, R. C. (2010). Search and clustering orders of magnitude faster than BLAST.
 Bioinformatics, 26(19), 2460–2461. https://doi.org/10.1093/bioinformatics/btq461

### OD-Seq
 OD-Seq is used for detecting mis-annotated or "outliers" in multiple sequence alignments when building new reference packages.
 Source files can be downloaded from the University College Dublin's website using this
 [link](http://www.bioinf.ucd.ie/download/od-seq.tar.gz).
 It can be compiled by either `make all` or in isolation with `make odseq`.

**Recommended version**: 1.0

**Citation**: Jehl, P., Sievers, F., & Higgins, D. G. (2015). OD-seq: Outlier detection in multiple sequence alignments.
 BMC Bioinformatics, 16(1), 1–11. https://doi.org/10.1186/s12859-015-0702-1

### BWA

BWA MEM is used for mapping short reads to classified DNA open reading frames (ORFs) if ORF prediction was performed
 and the '--rpkm' flag was used in `treesapp assign` or `treesapp abundance` was called.
BWA can be installed using [conda](https://anaconda.org/bioconda/bwa) or by following the instructions on the
 [GitHub page](https://github.com/lh3/bwa).

**Recommended version**: 0.7.17

**Citation**: Li, H. (2013). Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM.
 ArXiv Preprint ArXiv, 00(00), 3. https://doi.org/arXiv:1303.3997

## Finishing up
I hope that wasn't too painful. If you think you have installed everything, try running `treesapp info`!
It will check for the required executables up front and you will be quickly notified if some are missing or TreeSAPP is unable to find them.
