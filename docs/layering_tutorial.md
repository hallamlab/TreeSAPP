# Layering annotations onto classification outputs

This tutorial is to show the usefulness, and sometimes necessity, of `treesapp layer` to improve the classifications
 of TreeSAPP beyond functional and taxonomic annotations. 
With this module, users are able to use annotation files (colours_style.txt) from iTOL to add anything from 
 taxonomic guild to metabolic information, or really any aspect of a gene that is phylogenetically conserved,
 to the marker_contig_map.tsv file from `treesapp assign`.

## Ingredients

To use `treesapp layer`, you will need outputs created by `treesapp assign` as well as 
 annotation files created by annotating reference trees in iTOL.
If you have TreeSAPP cloned from GitHub somewhere, you can probably use a command like this to generate outputs we can layer.

```bash
treesapp assign -i ~/bin/TreeSAPP/test_data/marker_test_suite.faa \
-o marker_test/ \
-m prot -n 2 -t M0701,S0001 --trim_align
```

In this example, we will be using annotation files that come with TreeSAPP, 
 in the 'TreeSAPP/treesapp/data/iTOL_data/' directory.

## Usage

We have classified McrA (M0701) and DsrAB (S0001) sequences found in the file 'marker_test_suite.faa' so
 the next step is to add the metabolic annotations for McrA and the gene annotations for DsrAB.
The DsrAB reference package was created from seed sequences in the [FunGene](http://fungene.cme.msu.edu/) database for
 DsrA and DsrB, together in a single tree.
This works because they are paralogs that have evolved to form the holoenzyme
 necessary for dissimilatory sulfite reduction.

```bash
treesapp layer --treesapp_output marker_test/
```

Because we are using annotation files that are provided in the default `--annot_dir` directory,
 we do not need to provide any other arguments.
TreeSAPP looks in the directory for any files that do not end in the '_colours_style.txt' or '_colour_strip.txt' suffix
 and attempts to use these to layer the annotations.
Still, if you were analyzing a different reference package and had an annotation file outside of the default `--annot_dir`
 you could either change the path of `--annot_dir` or combine the annotation files in `--annot_dir` and yours with `-c/--colours_style`.
 
```bash
treesapp layer --treesapp_output marker_test/ -c DUF2112_Guilds.txt
```

TreeSAPP would then use all the annotation files in the default directory as well as 'DUF2112_Guilds.txt'.

## Creating your own annotation files

To create the files like [McrA_Metabolism.txt](https://github.com/hallamlab/TreeSAPP/blob/master/treesapp/data/iTOL_data/McrA_Metabolism.txt),
 you will need the reference tree from a reference package built using `treesapp create`.
This process is interactive, thanks to [iTOL](https://itol.embl.de/)(!),
 but can still be pretty time consuming especially 
 if there are misannotated sequences or gene sequences genes you didn't expect to find in your reference package.
For example, the DsrAB reference package contains not only proteins of DsrA and DsrB, but also AsrC (TIGR02912), 
 nitrite and sulphite reductases and many other sequences that were not expected by me (Connor) when first building the reference package.
Because they were unexpected, I also had to first figure out what they were!

Here is the most efficient process I've come up with.
It is unfortunately necessary to colour nodes in iTOL with the unlabelled tree as
 this is the only way that the annotation file can also be used in iTOL and not just by `treesapp layer`.
Begin by uploading the tree twice: once with no labels (just the numeric identifiers for each leaf),
 and again with the labels so you know who is who.
With the labels guiding you, on the unlabelled tree click on the lowest common ancestor nodes of clades that encompass
 a feature you are attempting to annotate.
From the drop-down menu, select 'Color' -> 'New colour range' and choose a colour and name for the feature.
If the feature has already been coloured elsewhere in the tree you are able to use 'Color' -> 'Add to existing range'.

Once you have finished colouring all the clades in the tree select 'Export' in the 'Controls' panel.
Change the format to 'Colors and styles annotation' and click 'export'. Now the last and most important step:
 name the file with the reference package name followed by an underscore, then the name of the feature and use a '.txt' extension.
Voila! Your own annotation file. 
