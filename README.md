# TreeSAPP: Tree-based Sensitive and Accurate Phylogenetic Profiler

![tests](https://github.com/hallamlab/TreeSAPP/workflows/tests/badge.svg)
[![PyPI version](https://badge.fury.io/py/treesapp.svg)](https://badge.fury.io/py/treesapp)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/treesapp/badges/version.svg)](https://anaconda.org/bioconda/treesapp)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/treesapp/badges/platforms.svg)](https://anaconda.org/bioconda/treesapp)
[![Docker Repository on Quay](https://quay.io/repository/hallamlab/treesapp/status "Docker Repository on Quay")](https://quay.io/repository/hallamlab/treesapp)

[![Codacy Badge](https://api.codacy.com/project/badge/Grade/b1937000c13040e8bba62f46e954796e)](https://www.codacy.com/gh/hallamlab/TreeSAPP?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=hallamlab/TreeSAPP&amp;utm_campaign=Badge_Grade)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/treesapp/README.html)
[![Python version](https://img.shields.io/pypi/pyversions/treesapp.svg)](https://img.shields.io/pypi/pyversions/)
[![codecov](https://codecov.io/gh/hallamlab/TreeSAPP/branch/master/graph/badge.svg)](https://codecov.io/gh/hallamlab/TreeSAPP)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/treesapp/badges/downloads.svg)](https://anaconda.org/bioconda/treesapp)

## Overview

TreeSAPP is a python package for functional and taxonomic annotation of proteins
 from genomes and metagenomes using phylogenetic placement.

## Quick start

We recommend installing TreeSAPP into its own conda environment with the following command:

```bash
conda create -n treesapp_cenv -c bioconda -c conda-forge treesapp
conda activate treesapp_cenv
```

To list all the sub-commands run `treesapp`.

To test the `assign` workflow, run:
```bash
treesapp assign -i TreeSAPP/tests/test_data/marker_test_suite.faa -m prot --trim_align -o assign_test -t McrA,DsrAB
```

To classify sequences in your genome of interest:
```bash
treesapp assign -i my.fasta -o ~/path/to/output/directory/
```

TreeSAPP comes installed with 33 reference packages involved in a variety of biogeochemical and cellular processes.
We also have many more reference packages available on our [RefPkgs repository](https://github.com/hallamlab/RefPkgs)
and you can view the complete list [here](https://github.com/hallamlab/RefPkgs/wiki/refpkgs).

## Tutorials

All of our tutorials are available on the [GitHub wiki](https://github.com/hallamlab/TreeSAPP/wiki) page.
Here are some specific tutorial examples:

If we do not yet have a reference package for a gene you are interested in,
please try [building a new reference package](https://github.com/hallamlab/TreeSAPP/wiki/Building-reference-packages-with-TreeSAPP).
Of course, if you run into any problems or would like to collaborate on building many reference packages
don't hesitate to email us or create a new issue with an 'enhancement' label.

To determine whether the sequences used to build your new reference package are what you think they are,
 and whether it might unexpectedly annotate homologous sequences,
 see the [purity tutorial](https://github.com/hallamlab/TreeSAPP/wiki/Testing-the-functional-purity-of-reference-packages).

If you are working with a particularly complex reference package, from an orthologous group for example, or have extra
 phylogenetic information you'd like to include in your classifications,
 try [annotating extra features](https://github.com/hallamlab/TreeSAPP/wiki/Layering-annotations-onto-classifications) with `treesapp layer`.

## Citation

If you found TreeSAPP useful in your work, please cite the following paper:

Morgan-Lang, C., McLaughlin, R., Armstrong, Z., Zhang, G., Chan, K., & Hallam, S. J. (2020). 
[TreeSAPP: The Tree-based Sensitive and Accurate Phylogenetic Profiler](https://doi.org/10.1093/bioinformatics/btaa588). 
Bioinformatics, 1–8.

This was brought to you by the team:

- Connor Morgan-Lang ([cmorganl](https://github.com/cmorganl), maintainer)
- Ryan McLaughlin ([McGlock](https://github.com/McGlock))
- Grace Zhang ([grace72](https://github.com/gracez72))
- Kevin Chan ([kevinxchan](https://github.com/kevinxchan))
- Zachary Armstrong
- Steven J. Hallam

### References

If you're feeling extra citation-happy, please consider citing the following works as well:

- Eddy, S. R. (1998). Profile hidden Markov models. Bioinformatics (Oxford, England), 14(9), 755–763.
- Criscuolo, A., & Gribaldo, S. (2010). BMGE (Block Mapping and Gathering with Entropy): A new software for selection of phylogenetic informative regions from multiple sequence alignments. BMC Evolutionary Biology, 10(1).
- Kozlov, A. M., Darriba, D., Flouri, T., Morel, B., & Stamatakis, A. (2019). RAxML-NG: a fast, scalable and user-friendly tool for maximum likelihood phylogenetic inference. Bioinformatics, 35(21), 4453–4455.
- Barbera, P., Kozlov, A. M., Czech, L., Morel, B., & Stamatakis, A. (2018). EPA-ng: Massively Parallel Evolutionary Placement of Genetic Sequences. Systematic Biology, 0(0), 291658.
