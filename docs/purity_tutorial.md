# Testing the functional purity of reference packages

This documentation covers the subcommand `treesapp purity`, used for determining the orthologous groups, genes, acivities etc.
 present in a reference package (RefPkg) by assigning sequences from a specially-formatted FASTA file.

This is an important step when annotating genomes and metagenomes using a gene-centric approach where the breadth of
reference sequences just isn't as comprehensive as public biological sequence repositories such as RefSeq.
When using a single RefPkg, the major hazard is annotating marginally homologous sequences
 because the more similar sequences are absent from the reference set.
This problem may be mitigated by using multiple RefPkgs that cover this sequence space or using a RefPkg that includes all homologs.
An example of the latter approach is the RefPkg DsrAB where we include both the oxidative and reductive forms of DsrA and DsrB
 (paralogs themselves) as well as a clade each of AsrC, nitrite and sulphite reductases and other COG2221-related.

Moreover, `treesapp purity` is a further validation step in cases where the repositories
 used to gather your candidate reference sequences were questionably curated.
This method can report homologs that you may not have been aware of and highlights instances where
 the HMM used is not as specific as we'd like*. 

## Ingredients

As mentioned, a specially-formatted FASTA file is required to perform the most basic purity-analysis.
In this FASTA, the functional information for each group must be the prefix of every header
 and the unique sequence identifiers can follow.
Here is an example snippet where sequences belong to the families TIGR00001, TIGR00014 and TIGR04570:
```
>TIGR00001_SP|P94976|RL35_MYCTU/2-64
PKAKTHSGASKRFRRTGT.GKIVRQKANRRHLLEHKPSTRT.RRLDGRTVVAANDTKRVT.SLLNG
>TIGR00014_GP|12544000|emb|CAC26382.1/3-114
VTIFHNPRCSTSRNTLAYLRDKDIEPEIVQYLKDTPTASELKELFNTLGIPV.HDGIRTREAEYTELGLS.PETPETELIDAIVAHPRLLQRPIVVTAKGARIARPKIDVIDSI
>TIGR04570_GB|CAE77370|MYCMY/1-87
MNEFNLAKDKTMISKIFKKIPWFYHLIFFLIGLVVGLLFQFLRVKTFAFPYFFILFFAVLLTYCVLFIIISPMIKQNWFIKRVKNEK
```

Users are able to use `TreeSAPP/dev_utils/TIGRFAM_seed_named.faa` that includes all the TIGRFAM seed sequences.
This file was created by concatenating all [TIGRFAM protein sequence files](ftp://ftp.jcvi.org/pub/data/TIGRFAMs/)
 together while prepending headers with their respective TIGRFAM identifier, separated with an underscore.
Additionally a tabular file with 3 columns, like that of `TreeSAPP/dev_utils/TIGRFAM_info.tsv`,
 can be used to provide the detailed functional information.
This was created using an ugly bash one-liner from the TIGRFAM INFO files.
 
We are not tied to either of these formats so if you have strong feelings to improve these,
 let us know on the issues page!
We are thinking of ways to make it easier for folks to create their own FASTAs of characterized sequences.
 
Crucially, though, these sequences _must be trustworthy_.
This analysis is only as good as the inputs, so if the reference sequences' activities haven't been characterized there
is no point in using them.
This is why we are only using the seed sequences where these sequence annotations are definitely correct. 

## Usage

Use `treesapp purity -h` to print the usage and command-line arguments available to you. 

In this example we will evaluate the purity of the RefPkg for the alpha subunit of methyl-coenzyme M reductase (McrA).
Since we have provided the characterized sequences and information table in `TreeSAPP/dev_utils`
 it is very easy to run after installing TreeSAPP.
The following command is ran from the TreeSAPP directory after cloning from GitHub and installing. 

```bash
treesapp purity -i dev_utils/TIGRFAM_seed_named.faa -x dev_utils/TIGRFAM_info.tsv \
-r McrA -o McrA_purity -p treesapp/data/ -m prot -n 2
```

`-p treesapp/data/` specifies the path to the RefPkg, indicating that this RefPkg was created and
 individual HMM, alignment, tree and lineage information files were copied to their respective directories
  (details can be found in the [treesapp create tutorial](https://github.com/hallamlab/TreeSAPP/blob/master/docs/Marker_package_creation_tutorial.md)).

## Output

The output should be:
```
Summarizing assignments for reference package McrA
Ortholog	Hits	Leaves	Tree-coverage	Description
--------------------------------------------------------------------------------
TIGR03256	4	4	1.8	methyl-coenzyme M reductase, alpha subunit
```

The _hits_ column is the number of characterized query sequences that were mapped to the tree and classified by `treesapp assign`.
The _leaves_ column is the number of reference sequences that were transitively assigned as that OG
 (by being a descendent of the node in the RefPkg's tree where an input sequence was placed).
This assignment is similar to how TreeSAPP assigns taxonomy to query sequences.

Also, all reference sequences in the RefPkg that were classified are listed in the log file (e.g. McrA_purity/TreeSAPP_purity_log.txt ),
 under the ortholog group (OG; first column in the above table) they were assigned. For example:

```
TIGR03256:
        Methanosarcina barkeri | CAA68357
        Methanothermus fervidus DSM 2088 | ADP77583
        Methanothermobacter thermautotrophicus | CAA30639
        Methanopyrus kandleri AV19 | AAM01870
``` 

These results indicate that the McrA RefPkg is pure. Also, since _hits_ is equal to _leaves_ and the _Tree-coverage_
 is low these characterized query sequences are very similar to the RefPkg references.

## Notes

\* This can be improved where reference sequences contain domains that are not tightly coupled to
 the activity (or activities) desired from sequences in the RefPkg.
Just trim that stuff off and rerun `treesapp create`! Unfortunately much of this is manual.
Feel free to contact us for help by email!