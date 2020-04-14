# Overview

TreeSAPP's gene-centric classification workflow assigns both a gene name and taxonomic label to proteins. It is capable of classifying either genomic or proteomic sequences in this manner, and they can be derived from genomes (isolates, single-cell amplified or metagenome-assembled) or metagenomes. The TreeSAPP subcommand used for classifying sequences is `treesapp assign`. Nucleotide reference packages are currently not supported.

***

# Workflow description

The standard workflow begins with conceptually translating open reading frames (ORFs) using Prodigal **[1]**.
The protein ORFs are then queried against a set of hidden Markov models (HMMs) using `hmmsearch` from the HMMER suite **[2]**. Proteins that are deemed homologous are then divided into groups based on the regions they aligned to in the HMM profile and these groups are used in a profile alignment with `hmmalign` to generate a multiple sequence alignment with reference and query sequences. RAxML's evolutionary placement algorithm (EPA) is used to insert the query sequences into the reference phylogeny **[3]**. Each placement, called a PQuery, is assigned a likelihood weight ratio, and distal and pendant length distances **[4]**. These are used for further filtering out non-homologous sequences.

PQueries are assigned a taxonomic label based on the lowest common ancestor (LCA) of their descendents. Furthermore, a taxonomic rank is recommended by a linear model that models the correlation between taxonomic rank and phylogenetic distance. This is useful in situations where a PQuery is placed on an edge with few children such that the LCA is highly resolved but the phylogenetic distance between the query and ancestral reference state is large.

Finally, classified sequences are written to a classification table and PQueries are concatenated into a single JPlace file for each reference package's classifications. The latter can be used for visualization with iTOL **[5]**.

Here is a diagram of the workflow:
![alt text](https://github.com/hallamlab/TreeSAPP/blob/master/docs/assign_workflow.png)

***

# Usage

`treesapp assign` only requires a single FASTA-formatted file containing either genomic or protein sequences.

### Basic usage

By default TreeSAPP assumes the inputs are genomic sequences and will run Prodigal.

```bash
treesapp assign -i my.fasta -o treesapp_output/
```

If you have protein sequences you can skip this step by including `-m prot`.

```
treesapp assign -i proteins.faa -o treesapp_output/ -m prot
```

### Recommended arguments

In addition to the basics, we tend to also extract phylogenetic informative regions from multiple sequence alignments using BMGE with the flag `--trim_align` **[6]**. This speeds up the phylogenetic placement stage without significantly harming classification performance.

Another favourite of ours is jacking up the number of threads available to Prodigal (FASTA is split into chunks and Prodigal is ran in parallel), HMMER, and RAxML with `--num_procs`.

For a list of all available arguments use `treesapp assign -h`.

### Targeted classifications

`treesapp assign` identifies homologs using all reference packages by default, although in many circumstances this is not necessary or desired. The argument for selecting the subset of reference packages to query is `-t`. This requires one or more comma-separated reference package codes, **not** the short-form gene or protein name. These codes can be found by running `treesapp info -v`. The bottom chunk of standard output lists the gene/protein names mapped to reference package codes, descriptions and other information. Here is an example:

```
S0001 -> DsrAB, prot, functional, Dissimilatory sulfite reductase alpha/beta subunits [EC:1.8.99.5], 16_May_2019
G0101 -> GH101, prot, functional, , 25_Jul_2019
G0163 -> GH163, prot, functional, , 25_Jul_2019
G0020 -> GH20, prot, functional, Glycoside hydrolase family 20, 26_Feb_2019
C0094 -> GH94, prot, functional, Glycoside hydrolase family 94, 16_Dec_2018
H0001 -> HydA, prot, functional, H-cluster of [FeFe]-hydrogenase (PF02906 & PF02256), 29_Jun_2019
M0701 -> McrA, prot, functional, Methyl coenzyme M reductase alpha subunit, 29_Apr_2019
M0702 -> McrB, prot, functional, Methyl coenzyme M reductase beta subunit, 05_Dec_2018
M0705 -> McrG, prot, functional, Methyl coenzyme M reductase gamma subunit, 05_Dec_2018
```

To analyze just the McrABG subunits (last three rows), you would include `-t M0701,M0702,M0705` in your `treesapp assign` command.

### Abundance estimates

`treesapp assign` is also capable of assigning FPKM abundance values to classified sequences when the `--fpkm` flag is invoked. At least one FASTQ file is also required (using the `-r` argument) and a second if the reads are paired end and not interleaved (path provided to `-2`). The `--pairing` option is used to specify whether the library is paired-end ('pe') or single-end ('se').

BWA is used to index the classified nucleotide ORF sequences then map reads to this index **[7]**. [`samsum`](https://github.com/hallamlab/samsum), a Python package dependency, parses the SAM file written and calculates FPKM values of each ORF.

This routine is only available if a genome was provided and Prodigal was ran - we cannot map DNA reads to protein sequences.

***

# Outputs

### Classification table

The table containing classifications, that is the output that most people would be interested in, is `final_outputs/marker_contig_map.tsv` under the output directory. Each row contains a single classified sequence with the gene/protein name and taxonomic label it was classified as. Chances are, you are just interested in the first six columns.

The column names and descriptions follow:

* `Sample`: Name of the sample the ORF originated from. This is taken from the name of the FASTA file provided to `-i`/`--fasta_input` except when the `--fpkm` flag is used, in which case the name is taken from the FASTQ file(s).
* `Query`: Name of the query protein sequence or conceptually translated ORF sequence.
* `Marker`: Name of the reference package (i.e. gene or protein name) the query sequence was classified as. If you are unfamiliar with this name, `treesapp info -v` can give you a brief description.
* `Length`: Length of the query sequence that was aligned to the reference profile HMM and subsequently used for phylogenetic placement. This could be, and usually is, a subsequence of the full-length ORF or protein sequence.
* `Taxonomy`: Taxonomic label assigned to the query sequence based on LCA alone. Should only be used for development and testing purposes.
* `Confident_Taxonomy`: Final taxonomic label assigned to the query sequence.
* `Abundance`: The abundance of the query sequence in the dataset. If `--fpkm` was not used it is set to 1, otherwise it is a fragments per kilobase per million mappable reads (FPKM) value.
* `iNode`: The internal node number representing where on the reference phylogeny the query sequence was placed.
* `LWR`: Likelihood weight ratio assigned to the query by RAxML-EPA.
* `EvoDist`: Sum of all phylogenetic distances calculated by RAxML-EPA as well as the average node-to-tip distances.
* `Distances`: A comma-separated list of the distances that summed to equal `EvoDist`

### iTOL inputs

Within the root output directory, there is a directory called `iTOL_output/`. Within there is a directory for each of the reference package names that query sequences were assigned to and within *those* are files you can use to visualize your placements in [iTOL](https://itol.embl.de/).

If you have an iTOL account already you should login first and once you are in your workspace click-and-drag the JPlace file from one of the iTOL output refpkg directories into a workspace. If you do not have an iTOL account (and don't want to get one) you can navigate to the "Upload" button at the top of the page and upload the JPlace file through that.

Once you are viewing the phylogeny you should toggle on the "Phylogenetic Placements" dataset. Then, to finish uploading the annotation files, click-and-drag the remaining files in that directory (e.g. McrA_labels.txt, McrA_colour_strip.txt, etc.) onto the webpage *et voila*, you have a b-e-a-utiful phylogeny showing you where sequences in your dataset were placed.

***

# References

**[1]** Hyatt, D., Chen, G.-L., Locascio, P. F., Land, M. L., Larimer, F. W., & Hauser, L. J. (2010). Prodigal: prokaryotic gene recognition and translation initiation site identification. BMC Bioinformatics, 11, 119. https://doi.org/10.1186/1471-2105-11-119

**[2]** Eddy, S. R. (1998). Profile hidden Markov models. Bioinformatics (Oxford, England), 14(9), 755–763. https://doi.org/btb114

**[3]** Berger, S. A., Krompass, D., & Stamatakis, A. (2011). Performance, accuracy, and web server for evolutionary placement of short sequence reads under maximum likelihood. Systematic Biology, 60(3), 291–302. https://doi.org/10.1093/sysbio/syr010

**[4]** Matsen, F. A., Hoffman, N. G., Gallagher, A., & Stamatakis, A. (2012). A format for phylogenetic placements. PLoS ONE, 7(2), 1–4. https://doi.org/10.1371/journal.pone.0031009

**[5]** Letunic, I., & Bork, P. (2019). Interactive Tree Of Life (iTOL) v4: recent updates and new developments. Nucleic Acids Research, 47(W1), W256–W259. https://doi.org/10.1093/nar/gkz239

**[6]** Criscuolo, A., & Gribaldo, S. (2010). BMGE (Block Mapping and Gathering with Entropy): A new software for selection of phylogenetic informative regions from multiple sequence alignments. BMC Evolutionary Biology, 10(1). https://doi.org/10.1186/1471-2148-10-210

**[7]** Li, H. (2013). Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. ArXiv Preprint ArXiv, 00(00), 3. https://doi.org/arXiv:1303.3997