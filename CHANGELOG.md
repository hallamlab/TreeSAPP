## [0.11.4] - 2022-05-22
### Added
- Centroid inference for pOTUs based on the midpoint, or balance point, of all cluster members.
- A table summarizing the intra-cluster evolutionary distances ('phylotu_cluster_stats.tsv').
- 

### Fixed
- Estimation of local-alignment distances for a taxonomic rank - now only considers sequences from monophyletic taxa.
- `--min_seq_length` always overrides the minimum profile HMM proportion threshold in `treesapp update`.

### Changed
- The minimum sequence length of 30 (AA) has been removed be default, but can still be used as before with `--min_seq_length`.
Results likely will not change as more stringent filtering thresholds were already applied in downstream steps.
- `--pc` was removed from `treesapp create`'s argument list

## [0.11.3] - 2021-12-15
### Added
- Option to use pairwise local-alignment clustering (with MMSeqs2) in `treesapp phylotu`

### Fixed
- Bad error statement when estimating _alpha_ in `treesapp phylotu`
- Can append unannotated features to "Unknown" label if already present in taxa_map.
- Problem updating a reference package with sequences from UniProt (>sp|... header).

### Changed
- Checkpointing is improved in `treesapp assign`.
  It is able to pick up outputs at any stage and decide what needs to be ran for each reference package.
  Reference package targets can be modified between reruns.

## [0.11.2] - 2021-06-12
### Added
- Option '--unknown_colour' for `treesapp colour` where a colour for "Unknown" features or taxa are included in the iTOL files.
- New options for pre-clustering the classified sequences using either Barberra et al.'s placement-space method or 
  pairwise alignment to speed up pOTU inference. Controlled with the "-p/--pre_mode" argument
- Dynamic evolutionary distance threshold for query sequences based on branch lengths descendent from placement position 
- RecA, RadA and RpoB reference packages being distributed as part of the core set
- The new '--query_coverage' command-line parameter is available in `treesapp assign` and drastically improves precision
  and recall in conjunction with '--hmm_coverage'. Both are set to 80% by default.
- '--delete' flag added to `treesapp phylotu` to optionally remove all intermediate files and directories.
  Useful for de novo methods when multiple phylogenies are inferred.
- Silent mode in `treesapp assign` can be activated by the '--silent' flag.
  No logging to console but log file is still populated.

### Fixed
- `treesapp package edit` assigns a leaf node only to the most resolved feature annotation
- Estimating `treesapp phylotu`'s alpha threshold is improved
- Setting distal and pendant lengths during aelw summary allows placements to be correctly filtered
- Final rank of a query sequence's assigned taxonomic lineage is not adjusted with aELW placement summary
- Detecting input sequence type for `treesapp evaluate`

### Changed
- Non-taxonomic features are coloured in alphabetical order (according to the palette used) in `treesapp colour`
- iTOL colour-strip files dataset labels are now the feature name
- Users are warned if multiple feature annotations are assigned to a leaf node during `treesapp package edit`
- `treesapp phylotu`'s _de novo_ pOTU workflow adds the most related reference sequences when inferring each phylogeny 
  to handle truncated query sequences

## [0.11.0] - 2021-04-27
TreeSAPP version 0.11.0 changes how users store and interact with reference package feature annotations.
These feature annotations are clade-specific labels that indicate some extra-taxonomic features that are characteristic of sequences in the reference package.

For example, in the particulate methane monooxygenase and ammonia monooxygenase subunit A reference package, XmoA, 
the feature annotations indicate which paralog is represented by a clade (PmoA, AmoA, EmoA, etc.)
As another example, the methyl coenzyme M reductase subunit A (McrA) reference package contains feature annotations for
each pathway of methanogenesis that is used by the different clades.

We recommend updating to this version, and updating reference packages you have created.

### Added
- A new attribute called 'feature_annotations' has been introduced to reference packages.
  It can store what was previously saved to iTOL-compatible annotation files by `treesapp colour`.
- `treesapp package edit` accepts a taxonomy-phenotype mapping file to populate the feature_annotations attribute.
  See [Wiki](https://github.com/hallamlab/TreeSAPP/wiki/Reference-package-operations) for details.
- `treesapp update` with automatically propagate feature annotations from the original reference package by mapping
  the reference sequences through their unique descriptions (organism name and accession).
-  `treesapp package view tree` will print a Newick tree with each leaf node's accession and description.
- `treesapp abundance` creates a simple_bar.txt file for each sample analyzed.
- Ability to automatically detect the sequence type based on the input provided.
- PQuery classification data is stored in each reference package in the 'training_df' attribute as a pandas.DataFrame.
- Improved query sequence filtering by phylogenetic placement information in `treesapp update`
- Now able to update a reference package's 'lineage_ids' attribute with `treesapp package edit`
- `treesapp create` is able to accept multiple fasta files through --fastx_input and concatenate them into the one 
  file used to build the reference package.

### Fixed
- Segmentation fault from Prodigal is no longer possible as `treesapp assign` verifies input presence earlier.
- `treesapp purity` bug where the reference package path was not correctly passed to `treesapp assign` if in the same directory
- Calculation of tree coverage in `treesapp purity`

### Changed
- Renamed the classification table made by `treesapp assign` (and used by subcommands like `layer`) 'classifications.tsv'.
- The reference package attribute 'refpkg_code' is automatically set and
  does not need to be changed as it is guaranteed to be unique.
- The reference package disband path has been changed to just the reference package code.
- `treesapp colour` accesses and uses the 'feature_annotations' to write iTOL-compatible annotation files
  (i.e. colour_strip.txt and colours_styles.txt). It no longer accepts taxonomy-phenotype tables.
- `treesapp layer` uses the 'feature_annotations' attribute in reference packages to annotate classified sequences.
- The versioned sequence accessions (or first split for unformatted sequence headers) are used in the
  ReferencePackage lineage_ids attribute. This ensures unique sequence IDs and helps with iterative updates.

## [0.10.4] - 2021-03-25
### Fixed
- Checkpoint determination in `treesapp abundance`
- '--report append' and 'report update' was not working properly in `treesapp abundance`.
  Fixed by deduplicating PQueries prior to appending. 

### Changed
- Checks whether all FASTQ file paths exist earlier in `treesapp abundance`

## [0.10.3] - 2021-03-23
### Fixed
- Replaced duplicate SAM file paths for unique ones when multiple fastqs are provided to `treesapp abundance`
- Prevent `treesapp abundance` from overwriting `treesapp assign` outputs when '--overwrite' is used

## [0.10.2] - 2021-03-08
### Added
- New flag '--deduplicate' for `treesapp create` to remove redundant sequences within 99.9% similarity before querying Entrez.
- New option for `treesapp assign` called '--hmm_coverage' that allows users to control minimum percentage of a
  profile HMM that a query sequence's alignment must cover

### Fixed
- Mapping some EggNOG identifiers with seqs2lineage mapping process

### Changed
- Reference package training is now optional with `treesapp update`. A fasta file is no longer required.

## [0.10.1] - 2021-02-15
### Added
- Control of relative abundance metric for `treesapp abundance` to use with `treesapp assign`

### Fixed
- Using all predicted ORFs when calculating abundance values for proper TPM values
- Fixed TypeError when `treesapp abundance` is used with single-end or interleaved FASTQ files

### Changed
- Using miniconda for installing dependencies in GitHub Action 'tests'

## [0.10.0] - 2021-02-09
### Added
- (#33) Checkpointing to all major (i.e. time-consuming) subcommands: `abundance`, `assign`, `create`, `purity`, `train`
- Added `stages` argument to `treesapp abundance` to control checkpoint
- Able to control the maximum number of examples used and SVC kernel for training through `treesapp create`
- (#66) Documentation for `treesapp evaluate`
- A new _de novo_ clustering mode for `treesapp phylotu` will infer a new tree of just query sequences
- A new `treesapp phylotu` output mapping classified sequences to their cluster
- Able to append results from `treesapp abundance` to classification table for multiple read datasets

### Fixed
- (#71) Ability to rerun and append results from `treesapp evaluate` has been restored
- File paths with square brackets and parentheses no longer trip up HMMER or RAxML-NG in various modules.
  Some single-quotes were needed.
- Properly truncate sequences in `treesapp evaluate` with '--length' argument
- Can add bipartition support to JPlace trees for visualization. Bipartitions were not formatted properly.
- File paths of reference package trees with support values

### Changed
- `treesapp evaluate` uses a tqdm progress bar instead of printing updates to stdout when classifying
- Removed the prodigal header tags when ORFs are predicted by TreeSAPP. Could lead to significant RAM reduction.
- Flag to activate relative abundance calculation in `treesapp assign` changed from '--rpkm' to '--rel_abund'

## [0.9.8] - 2021-01-15
### Added
- More unit tests for file parsers and functions for `treesapp assign`
- More rigorous integrative test for `mcc_calculator.py` with higher coverage.

### Fixed
- A KeyError in `treesapp update` caused when ORFs on the same sequence classified to the same refpkg.
  Their dictionary keys were overwritten in simulate_entrez_records().
- Mapping sequence names to lineages provided in a seqs2lineage file.
  ORF names containing parentheses were not being escaped, messing up regular expression matching.
- (#64) Able to overwrite reference packages in the same directory
- `mcc_calculator.py` was not calculating TP, TN and FP correctly.

### Changed
- Only `max_lwr` mode in `treesapp assign` will use the linear model for taxonomic rank recommendation
- Warning in `treesapp update` if all query sequences are shorter than the minimum sequence length threshold
- Log files don't contain the sequence names removed from FASTA objects at various steps,
  potentially significantly reducing the size of these files.
- `mcc_calculator.py` summarizes the number of queries with missing taxonomic lineages in a single warning.
- Log writes to stderr while `treesapp package view` writes reference package data to sys.stdout.

## [0.9.7] - 2021-01-06
### Added
- New option in `treesapp package` called 'rename' to change the name of a reference package
- MMSeqs2 and MAFFT versions included in `treesapp info` output
- GitHub Actions for CI
- Unit tests

### Fixed
- (#63) An annoying warning that a fasta with classified nucleotide sequences couldn't be created for protein sequence input.
- (#62) Corrected FPKM values.
- (#58) List available reference packages (link to [RefPkgs](https://github.com/hallamlab/RefPkgs) in README).

### Changed
- Updated required version of `samsum` to 0.1.4

## [0.9.0] - 2020-09-15

This release is for version 0.9.0 and covers many new features. It is strongly recommended that everyone updates to this version.

There is __no__ backwards compatability of between this and older versions.
All reference packages built using older versions will need to be remade from scratch.

### Added
- Uses EPA-NG and RAxML-NG for phylogenetic placement and inference, respectively. Using EPA-NG drops the runtime drastically.
- (#48) Reference packages are stored as a single, pickled file.
- The subcommand treesapp package has been specifically designed to allow users to still interact with this binary file.
- (#49) Reduced RAM usage by not loading query FASTA, rather loading just headers and loading only the query sequences that matched a reference package's profile HMM.
- Using the RefPkgs repository for version control of all non-core reference packages. Available for everyone to contribute to!
- Profile HMMs are dereplicated at the genus rank for more sensitive profiles.
- Users are able to update reference packages with GTDB-tk lineages by using the --seqs2lineage argument.
- The subcommand treesapp train fully supports checkpointing. Checkpointing in other subcommands is still to come (but well on its way).
- We finally have some sort of a test suite.

### Fixed
- (#51) Error while parsing some alignments from hmmsearch

### Changed
- Format of the classification table (final_outputs/marker_contig_map.tsv) has been changed. Only a single column for the recommended taxonomy and the hmmsearch-derived E-value is reported.
- Reference trees are automatically rooted using ETE3's 'set_outgroup' function with the farthest node. Polytomies are also automatically resolved when the tree is built using FastTree.
- In treesapp create (and therefore treesapp update) reference sequence outlier detection and removal using OD-Seq has been made optional and can be requested using the '--outdet_align' flag.