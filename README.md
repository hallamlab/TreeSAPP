# MLTreeMap: A python implementation accompanied by algorithmic improvements

Kishori M. Konwar, Young C. Song, Connor Morgan-Lang, and Steven J. Hallam

## Overview:

A unified python software for generating phylogenetic trees from genomic sequences. 

## Download and installation:

```
git clone git@github.com:hallamlab/MLTreeMap.git
```
This is the only command that is required for machines running MacOSX. The exectutables for the
required softwares are included in the MLTreeMap/sub_binaries directory. However, for Linux machines
these executables will not work. They must be installed separately and ideally in the sub_binaries
directory (with the included executables removed).
### Downloading dependencies for Linux:
#### Gblocks:
For Linux/x86:
```
cd path/to/MLTreeMap/sub_binaries/
wget http://molevol.cmima.csic.es/castresana/Gblocks/Gblocks_Linux_0.91b.tar.Z
tar xzf Gblocks_Linux_0.91b.tar.Z
rm Gblocks
ln -s Gblocks_0.91b/Gblocks ./
```

#### LAST:
```
wget http://last.cbrc.jp/last-588.zip
unzip last-588.zip
cd last-588/
make
cd ../
ln -s last-588/src/lastal ./
ln -s last-588/src/lastdb ./
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
If you have problems with involving `getline` being previously declared follow these instructions:
http://www.langebio.cinvestav.mx/bioinformatica/jacob/?p=709

#### HMMER
hmmalign is the only HMMER module required by MLTreeMap, but HMMER3 is incompatible with this
version of MLTreeMap. HMMER 2.4 works and can be downloaded from
http://hmmer.janelia.org/download.html.


