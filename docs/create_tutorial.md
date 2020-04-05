# Building new reference packages for TreeSAPP

This tutorial is meant for users who would like to analyze a gene family lacking an associated TreeSAPP reference package (refpkg).
 This tutorial is also useful to people who find a refpkg is
 completely useless to them and would rather start over than use `treesapp update`.

## Ingredients

The subcommand we will need to use for this tutorial is `treesapp create`.
This subcommand's workflow is visualized below:

![alt text](https://github.com/hallamlab/TreeSAPP/blob/master/docs/create_workflow.png)

As you can see `treesapp create` depends on several other open source software.
Most of these are also dependencies of `treesapp assign` (HMMER, RAxML, BMGE) but some are not (FastTree, MAFFT, USEARCH, OD-Seq).

To promote transparency and reproducibility, TreeSAPP assumes that the sequences (either nucleotide or protein)
 are accessioned on the NCBI and therefore available through BioPython's Entrez API.
This allows TreeSAPP to automatically determine the taxonomic lineage assigned to each sequence without this information
 being provided. Hopefully this means less work for you!
However, seeing as folks will invariably want to work with the latest and greatest sequences there is a way to 
 include lineage information in the FASTA input for as-of-yet unaccessions sequences (below).

## Step 1. Downloading the sequences

### What databases can I retrieve sequences from?

Throughout the development of TreeSAPP I have frequented
 [FunGene](http://fungene.cme.msu.edu/), the functional gene pipeline and repository.
This database provides curated sequences (of varying reliability) for many popular function and taxonomic anchor genes.
FunGene is also great since they have hidden Markov models (HMMs) for every gene available,
 providing a sense of what the full length gene is and allowing users to
 filter the sequences based on the percentage of an HMM covered.
Since we're building a reference tree I suggest downloading all full
 length and near-full length genes, so setting that HMM coverage parameter to between 60 and 90%.

Apart from FunGene, the NCBI itself is another great option though you must be extremely wary of the sequences it returns.
 Non-homologous sequences are often returned by queries so while aggressive length and taxonomic filters are your friend,
manual review is required.

### Alternative databases

FunGene and the NCBI are my favourite options, though this leaves out
the vast majority of biological sequence databases for no reason other than simplicity.
IMG, KEGG, UniProt, Pfam, and other similar databases are all great but they lack an API to retrieve lineage information,
 and therefore this needs to be provided in the header.
 
 If you want to include sequences that haven't been uploaded and accessioned in Entrez,
  you can do so by modifying the sequence's FASTA header to follow the format:
 >SeqID lineage=cellular organisms; Domain; Phylum; Class [Organism_name]
 
 where `SeqID` should be replaced with a temporary, unique accession or ID (e.g. AMH87091),
  "Domain; Phylum; Class" need to be replaced with the appropriate values for the organism this sequence was derived from, 
  and `Organism_name` should be replaced with the appropriate organism name, such as 'Candidatus Moimicrobe mangepain'.
  No quotes required, though.

## Step 2. Creating the reference package

The most basic refpkg consists of 4 files:
a multiple sequence alignment file in FASTA format, a HMM, a tree in NEWICK format, and a table mapping identifiers
 in the phylogenetic tree to taxa. In addition to these files, the same NEWICK tree with bootstrap files can also be generated.

A basic command could look something like this:
```
$ treesapp create \
 --fasta_input rpoB_proteins.faa --output_dir rpoB_TreeSAPP_create \
 --refpkg_name rpoB --cluster --identity 0.90 \
 --num_procs 8 --molecule prot
```

Other options that some may find useful are `--profile`, `--guarantee`, `--accession2taxid`, `--fast`, and `--headless`.

`--profile` can be used if you have a HMM that you want to use to extract domains from the input sequences.
For example, maybe input sequences contain multiple different domains or multiple copies of the same domain.
In some cases this is desired while in others it could complicate phylogenetic inference and the profile HMM,
 and may lead to false positive classifications.
It is advised to use an HMM that models the domain you're interested in and to build a refpkg from these sequences only.
This is also handy if you don't fully trust the sequences you're using to build a refpkg.

`--fast` invokes "fast mode" where FastTree is used to build the phylogenetic tree instead of RAxML.
 The main advantage here is speed. Currently a tree with bootstrap values is not generated in fast-mode.
We have not observed a significant difference in classification performance between RAxML and FastTree phylogenies,
 so if you are building a large refpkg I suggest using fast mode.

`--headless` is useful when building refpkgs on a server and don't want to be prompted for any reason.
Without `--headless`, when the `--cluster` flag is used the sequences of each cluster are presented to the user for them to manually select
the sequence with the best annotation, sequence length, etc.
 There are cases when this isn't feasible and instead the longest sequence (cluster centroid) is used.

`--guarantee` can be used to ensure that sequences make it into the final set of reference sequences,
 regardless of whether they would be normally removed by filters and clustering.
It's argument is a FASTA file with the sequences that you want to make it through.
These do not need to be in the FASTA provided to `--fasta_input`.

`--accession2taxid` is a path to an NCBI accession2taxid file for more rapid accession-to-lineage mapping.
These files can be downloaded from the NCBI's ftp site [here](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/).
We use prot.accession2taxid.
Note: its only worthwhile if you have many hundreds or thousands of sequences in your input that need to be assigned a lineage.
Otherwise it is probably faster to just wait for TreeSAPP to send queries to the Entrez database.

## Step 3. Integration

The final step is to run those final commands that were printed as `treesapp create` completed.
They should look something like this:
```
To integrate these data for use in TreeSAPP, the following steps must be performed:
1. Include properly formatted 'denominator' codes in data/tree_data/ref_build_parameters.tsv
2. $ cp rpoB_create_85/TreeSAPP_files_rpoB/tax_ids_rpoB.txt TreeSAPP/treesapp/data/tree_data/
3. $ cp rpoB_create_85/TreeSAPP_files_rpoB/rpoB_tree.txt TreeSAPP/treesapp/data/tree_data/
4. $ cp rpoB_create_85/TreeSAPP_files_rpoB/rpoB.hmm TreeSAPP/treesapp/data/hmm_data/
```

These need to be executed outside of `treesapp create` since it looks for all refpkg files in subdirectories
of a single location. These are not copied automatically for two reasons: 
1. Prevent over-writing existing refpkgs that users may want to retain (sentimental value, etc.)
2. Give the user an opportunity to inspect and test the refpkg before (blindly) using it for analyses

This brings us, dear reader, to evaluating the quality of refpkgs with [`treesapp purity`](https://github.com/hallamlab/TreeSAPP/blob/master/docs/purity_tutorial.md).

An upcoming tutorial will show you how to estimate the taxonomic classification accuracy with your new refpkg using `treesapp evaluate`. 
  
## Creating multiple versions

A significant time sink is waiting for TreeSAPP to download the lineages associated with each sequence accession.
This time can be saved if you're building multiple refpkgs from the same FASTA file,
 say by  clustering at different sequence similarities.
After creating the first refpkg just make the output directories for the other versions and copy the
 `accessions_id_lineage_map.tsv` to those directories. Here is an example:

```
$ treesapp create \
 -i rpoB_proteins.faa -o rpoB_create_90/ \
 -n 4 -c rpoB -p 0.90 --cluster --taxa_lca
$ mkdir rpoB_create_80/
$ mkdir rpoB_create_70/
$ cp rpoB_create_90/intermediates/accessions_id_lineage_map.tsv rpoB_create_80/
$ cp rpoB_create_90/intermediates/accessions_id_lineage_map.tsv rpoB_create_70/
```

Following this, run `treesapp create` for the other versions.
