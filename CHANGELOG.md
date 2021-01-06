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
