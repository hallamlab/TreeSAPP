# MLTreeMap: A python implementation accompanied by algorithmic improvements

Kishori M. Konwar, Young C. Song, Connor Morgan-Lang, and Steven J. Hallam

## Overview:

A unified python software for generating phylogenetic trees from genomic sequences. 

## Download and installation:

```
git clone git@github.com:hallamlab/MLTreeMap.git
```
The exectutables for the required softwares are included in the MLTreeMap/sub_binaries/mac directory.
The Linux executables are in MLTreeMap/sub_binaries/ubuntu. If these do not work out-of-the-box, instructions
on installing the dependencies specific to your machine are included below.
### Downloading dependencies for Linux:

#### RAxML:
A simple `git clone` of their [github page](https://github.com/stamatak/standard-RAxML) should work.
We have tested several versions and found no problems from V.7.1 to the most recent release as of 
December 1st, 2015. However, the executable MUST be named `raxmlHPC` or it will not be found by MLTreeMap!
If you have find an incompatability please notify us through the Issues feed!

#### Gblocks:
For Linux/x86:
```
cd path/to/MLTreeMap/sub_binaries/
wget http://molevol.cmima.csic.es/castresana/Gblocks/Gblocks_Linux_0.91b.tar.Z
tar xzf Gblocks_Linux_0.91b.tar.Z
rm Gblocks
ln -s Gblocks_0.91b/Gblocks ./
```

#### Genewise:
You can either install the package with apt-get:
```
sudo apt-get install wise
```
Or try to install wise from source
```
wget http://www.ebi.ac.uk/~birney/wise2/wise2.4.1.tar.gz
tar xzf wise2.4.1.tar.gz
rm wise2.4.1.tar.gz
cd wise2.4.1/src/
make all
```
If you have problems involving `getline` being previously declared follow these instructions:
http://www.langebio.cinvestav.mx/bioinformatica/jacob/?p=709

#### HMMER
hmmalign is the only HMMER module required by MLTreeMap, but HMMER3 is incompatible with this
version of MLTreeMap. HMMER 2.4 works and can be downloaded from
http://hmmer.janelia.org/download.html.

### Running MLTreeMap

To list all the options with brief help statements `./mltreemap.py -h`.

To perform a basic run with only the required arguments:
```
./mltreemap.py -i input.fasta -o ~/path/to/output/directory/
```
As you can see, executables are automatically detected in both the $PATH and in the
sub_binaries/mac or sub_binaries/ubuntu, depending on your OS. However, if your executables
are elsewhere, MLTreeMap can be directed to them through the `--executables` argument.
