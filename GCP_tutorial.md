# Running TreeSAPP on Google Cloud Platform

This tutorial is meant to allow people to run TreeSAPP in the cloud with as few installation steps as possible
and without the worry of interfering with current software installations.

__Note__: this page is a work in progress.
I will be adding additional details on connecting to data buckets and more!

## Account setup:

Begin by creating a new instance on [Google Cloud Platform](https://cloud.google.com/).
If you do not already have an account you will need to create an account (quick guide [here]()).

To create a new instance, navigate to your GCP "Dashboard" within the "Console".
It should look something like this:

![dash_pic](/home/connor/Pictures/MICB425_treesapp/dashboard.png "Caption")

Under the "Resources" tile, select "Compute Engine".
At the top of this new page, there will be a dashboard with a "CREATE INSTANCE" tab. Here is where it gets interesting!

Pick a name and zone (I chose 'us-west1-b'). Now is the first of two important options: Machine type.
TreeSAPP has relatively low RAM requirements but does benefit greatly from parallelization.
I suggest either one of the 8 or 16 high CPU options (code named n1-highcpu-8 or n1-highcpu-16, respectively)
from the dropdown menu for simplicity, though feel free to create a custom type to suite your needs.
You do not want to deploy a container to this instance. The second important option is the "Boot disk":
I recommend selecting "Ubuntu 16.04 LTS" since TreeSAPP is known to work with this OS.
Keeping the disk type "Standard persistent disk", increase the Size to 40 Gb.
Finally, change the access scope from "Allow default access" to "Allow full access to all Cloud APIs".

You may have noticed an estimated cost on the right side of the page.
This is an estimated cost for running the instance continuously for a whole month.
Fortunately, Google does offer $300 USD in free compute time to familiarize yourself with the platform!
Click the "Create" button at the bottom of the page to spin up this instance.

## Connecting to the instance

This is relatively straightforward once Google's GCP utility suite, gsutil, is installed.
Here is the [link](https://cloud.google.com/storage/docs/gsutil_install) to the quick installation guide.

Once this is installed, the command to connect to your server is:
`gcloud compute ssh my_instance` where my_instance is the name you chose while creating your instance.

## Buckets of fun

GCP does not appreciate storing data on their compute servers. For this reason
they have implemented data storage on specialized storage servers.
We are able to create a partition on their storage servers that meet our needs,
and these partitions are called buckets.

You can create a bucket in the GCP console by navigating to Products & services (top left menu icon),
then Storage. This will bring you to a page where you can create a bucket.

Once you have a bucket, these can be easily mounted to a compute instance using `gcsfuse`.
`gcsfuse` does not come pre-installed,
though the [process is trivial](https://github.com/GoogleCloudPlatform/gcsfuse/blob/master/docs/installing.md).

## Installation:

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


## Finishing up

When you are finished processing all your samples, make sure you STOP your instance
so you don't continue incurring charges from google to your credit card
(or at least it stops chewing up your free compute time!). To do this, navigate back to your "Compute Engine" console,
selectyour instance and hit the "STOP" button at the top of your page.
 It may take up to a minute to successfully shut it down and I suggest making sure it does.
That's it!