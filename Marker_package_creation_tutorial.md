# Building new reference packages for TreeSAPP

This tutorial is meant for users who would like to analyze a marker gene
 with TreeSAPP that currently doesn't have an associated reference package.
 This tutorial is also useful to people who find a reference package is
 completely useless to them and would rather start over than update.

## Ingredients

The script we will need to use for this tutorial is `create_treesapp_refpkg.py`.
 It is packaged within the TreeSAPP GitHub repository and has been used extensively to
build all the reference packages available. The most basic usage of this pipeline
is visualized below:

![alt text](https://github.com/hallamlab/TreeSAPP/blob/master/dev_utils/Create_TreeSAPP_RefPkg_pipeline.png)

As you can see `create_treesapp_refpkg.py` depends on several other open source software.
Many of these are also dependencies of `treesapp.py` (e.g. HMMER, RAxML)
but some are not. These are: MAFFT, TrimAl, USEARCH and potentially
FastTree if you're interested in `fast` mode as well as the Python package BioPython (>1.67).

Last but not least: the input sequences. To promote transparency and reproducibility,
TreeSAPP assumes that the sequences (either nucleotide or protein) are accessioned
on the NCBI and therefore, available through Entrez through BioPython's
Entrez API. This allows TreeSAPP to automatically determine the
taxonomic lineage assigned to each sequence without this information being provided by other means
 (manual or in a mapping file... basically leading to more work for you!).
 However, seeing as folks will invariably want to work with the latest and greatest sequences,
 and don't like being restricted there is a way to semi-automatically include lineage information
 in the FASTA input for as-of-yet unaccessions sequences.

## Step 1. Downloading the sequences

### What databases can I retrieve sequences from?

Throughout the development and publication of TreeSAPP I have frequented
[FunGene, the functional gene pipeline and repository](http://fungene.cme.msu.edu/).
This database provides curated sequences for many popular genes,
beyond just those of biogeochemical relevance. FunGene is also great
since they have hidden Markov models (HMMs) for every gene available,
providing a sense of what the full length gene is and allowing users to
filter the sequences based on the percentage of an HMM covered.
Since we're building a reference tree I suggest downloading all full
 length and near-full length genes, so setting that HMM coverage parameter to between 60 and 90%.

Apart from FunGene, the NCBI itself is another great option though you
must be extremely wary of the sequences it returns. Totally unrelated
sequences are often returned so while aggressive length and taxonomic filters are your friend,
manual review is required. I'd reserve this for "experts only".

### Alternative databases

FunGene and the NCBI are my favourite options, though this leaves out
the vast majority of biological sequence databases for no reason other than simplicity.
 IMG, KEGG, UniProt, Pfam, and other similar databases are all great but,
  in my experience I have been unable to find an API for these to retrieve lineage information,
 and therefore this needs to be provided in the header.

## Step 2. Creating the reference package

A "reference package" consists of 5 files:
A multiple sequence alignment file in FASTA format, a HMM, a table mapping tree node numbers to taxa,
 a tree in NEWICK format, and a NEWICK tree with bootstrap values.

A basic command could look something like this:
```
$ ./create_treesapp_refpkg.py \
 --fasta_input rpoB_proteins.faa --output_dir rpoB_TreeSAPP_create \
 --code_name rpoB --identity 90 \
 --num_threads 8 --verbose
```

## Step 3. Integration

Congrats, you've created your TreeSAPP reference package! The final step
is to run those commands that were printed as `create_treesapp_refpkg.py` finished.
They should look something like this:
```
To integrate these data for use in TreeSAPP, the following steps must be performed:
1. Include properly formatted 'denominator' codes in data/tree_data/cog_list.tsv and data/tree_data/ref_build_parameters.tsv
2. $ cp rpoB_create_85/TreeSAPP_files_rpoB/tax_ids_rpoB.txt data/tree_data/
3. $ cp rpoB_create_85/TreeSAPP_files_rpoB/rpoB_tree.txt data/tree_data/
4. $ cp rpoB_create_85/TreeSAPP_files_rpoB/rpoB.hmm data/hmm_data/
```

## Creating multiple versions

A significant time sink is waiting for TreeSAPP to download the lineages
 associated with each sequence accession. This time can be saved if you're
 building multiple reference packages from the same FASTA file, say by
 clustering at different sequence similarities.
 After creating the first reference package just make the output directories
 for the other versions and copy the `accessions_id_lineage_map.tsv` to
 those directories. Here is an example:

```
$ ./create_treesapp_refpkg.py \
 -i rpoB_proteins.faa -o rpoB_create_90/ \
 -T 4 -c rpoB -p 90 --cluster --verbose --taxa_lca
$ mkdir rpoB_create_80/
$ mkdir rpoB_create_70/
$ cp rpoB_create_90/accessions_id_lineage_map.tsv rpoB_create_80/
$ cp rpoB_create_90/accessions_id_lineage_map.tsv rpoB_create_70/
```

Following this, run `create_treesapp_refpkg.py` for the other versions.
