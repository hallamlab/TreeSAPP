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