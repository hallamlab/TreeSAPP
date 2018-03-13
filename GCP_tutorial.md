# Running TreeSAPP on Google Cloud Platform

This tutorial is meant to allow people to run TreeSAPP in the cloud with as few installation steps as possible
and without the worry of interfering with current software installations.

__Note__: this page is incomplete but the commands listed below will work.
I will be adding additional details on connecting to data buckets and selecting an appropriate
instance for TreeSAPP analyses.


## Installation:

Begin by creating a new instance on [Google Cloud Platform](https://cloud.google.com/).
If you do not already have an account you will need to create an account. Forunately, Google does offer
student accounts with free compute time to familiarize yourself with the platform!

My cloud instance is running Ubuntu 16.04 with 24 CPUs and 22Gb of RAM.

Once you have this instance started, log on to the node via `ssh`.

From here, it is just a matter of downloading and installing the dependencies of TreeSAPP
before compiling TreeSAPP itself.

1. `sudo apt-get update; sudo apt-get install python3-pip`
2. `sudo pip3 install numpy biopython`
3. `git clone https://github.com/hallamlab/TreeSAPP.git; cd TreeSAPP; make; make install`
4. `tar xzf sub_binaries/wise2.2.0.tar.gz; cd wise2.2.0/src/; make all`
5. `mv bin/genewise ~/TreeSAPP/sub_binaries/; cd ~/TreeSAPP/; rm -r wise2.2.0/`
6. `tar xzf sub_binaries/Gblocks_Linux_0.91b.tar.Z; cp Gblocks_0.91b/Gblocks sub_binaries/; rm -r Gblocks_*`
7. `sudo apt-get install infernal`
8. `sudo apt-get install ncbi-blast+`
9. `sudo apt-get install default-jre-headless`
10. `wget https://www.drive5.com/muscle/downloads3.8.31/muscle3.8.31_i86linux64.tar.gz; tar xzf muscle3.8.31_i86linux64.tar.gz; mv muscle3.8.31_i86linux64 sub_binaries/muscle; rm muscle*gz`
11. `tar xzf sub_binaries/hmmer-2.4i.tar.gz; cd hmmer-2.4i/; ./configure; make; make install; cp src/hmmbuild ../sub_binaries/; cd -`

If something did not complete successfully please create an issue on the GitHUb page!
As I mentioned, this was successful using Ubuntu 16.04 but the basics should work for most other operating systems supported by GCP.
Of course, the `apt` package manager will need to be swapped if using a different OS (e.g. `yum` for RedHat) and the package names will
likely change.

## Test the installation

At this point, everything should be installed with no errors (warnings are *probably* okay).

To ensure TreeSAPP will work with your data I have several test fasta files, any one of which should validate the installation.

Here is an example:

`./TreeSAPP/treesapp.py -i TreeSAPP/test_data/marker_test_suite.fna -o marker_test -T 8 --verbose`

This will analyze a FASTA file containing nucleotide sequences of many different marker genes.
The `-o` argument specifies the name of the output directory, `-T` indicates 8 threads should be used when possible,
and we would like a verbose runtime log with `--verbose`.

For more information on the parameters available, use `~/TreeSAPP/treesapp.p -h`.

## Analyzing your data

This section will describe how to mount a "bucket" for your data to your compute instance in the case you selected
a small hard disk to save money (which I suggest).
