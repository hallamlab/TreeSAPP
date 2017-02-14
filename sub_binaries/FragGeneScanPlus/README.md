#FragGeneScan-Plus for scalable high-throughput short-read open reading frame prediction

Dongjae Kim, Aria S. Hahn, Shang-Ju Wu, Niels W. Hanson, Kishori M. Konwar, and Steven J. Hallam 

A fundamental step in the analysis of environmental sequence information is the prediction of potential genes or open reading frames (ORFs) encoding the metabolic potential of individual cells and entire microbial communities. FragGeneScan, a software designed to predict intact and incomplete ORFs on short sequencing reads combines codon usage bias, sequencing error models and start/stop codon patterns in a hidden Markov model to find the most likely path of hidden states from a given input sequence, provides a promising route for gene recovery in environmental datasets with incomplete assemblies. However, the current implementation of FragGeneScan does not scale efficiently with increasing input data size. This limits application of FragGeneScan to contemporary environmental datasets that can exceed 100s of Gb. Here, we present FragGeneScan-Plus, an improved implementation of the FragGeneScan gene prediction model that leverages algorithmic thread synchronization and efficient in-memory data management to utilize multiple CPU cores without blocking I/O operations. FragGeneScan-Plus can process data approximately 5-times faster than FragGeneScan using a single core and approximately 50-times faster using eight hyper-threaded cores when benchmarked against simulated and real world environmental datasets of varying complexity.

More information can be found on the [Wiki](https://github.com/hallamlab/FragGeneScanPlus/wiki).

Performance analysis and code can be found in [fgsp_results/](fgsp_results/)

If using FragGeneScanPlus for academic work, please cite:

Dongjae Kim, Aria S. Hahn, Shang-Ju Wu, Niels W. Hanson, Kishori M. Konwar, Steven J. Hallam. *FragGeneScan+: high-throughput short-read gene prediction*, Proceedings of the 2015 IEEE Conference on Computational Intelligence in Bioinformatics and Computational Biology (CIBCB 2015), Niagara Falls, Canada, August 12-15, 2015. [doi:10.1109/cibcb.2015.7300341](http://dx.doi.org/10.1109/cibcb.2015.7300341).
