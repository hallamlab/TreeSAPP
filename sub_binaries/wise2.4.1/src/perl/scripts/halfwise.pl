#!/usr/local/bin/perl

#
# These are the things you might want to change. Scroll down for docs on the program
# or preferable go pod2text halfwise.pl
#

# so that the testing works at sanger centre. You could change this obviously 
# on your site, though it may be easier to install bioperl correctly

use lib '/nfs/disk92/PerlSource/Bioperl/Releases/bioperl-0.04.2/';


BEGIN {
  eval {
    use Bio::Tools::Blast;
  };
  if ( $@ ) {
    print STDERR ("\nThe bioperl distribution has not be installed.\n Read the installation instructions at the top of this script going pod2text halfwise.pl\nThe bioperl distribution can be found from http://bio.perl.org/\nThe halfwise script has been tested with 0.04.x releases\n");
    exit(1);
  }
}


$defopts  = "-init wing -cace -cut 25 -aln 200 -quiet";
$half_db  = "/somewhere/halfwise.db";

# this is for testing @sanger
#$half_db  = "/nfs/disk21/birney/prog/wise2/perl/scripts/halfwise.db";

$hmmfetch = "hmmfetch";
$genewise = "genewisedb";
$hmmdb    = "Pfam";


#
# The documentation is written as pod. 
# to view it nicely use pod2text halfwise.pl
# or to make html use pod2html
#


=head1 NAME - halfwise

halfwise is a script designed to compare a large DNA sequence to Pfam
sensibly. It does the following steps:

=over

=item blast search

Compares the DNA sequence to a representive database of protein sequences
from Pfam

=item HMM retrieval

Retrieves the Pfam HMMs to compare to the DNA sequence using a liberal
threshold

=item GeneWise comparison

Retrieves those HMMs and uses genewisedb to compare the DNA sequence to the
Pfam HMMs.

=back

The blast search drastically reduces the number of HMMs needed to be considered
while not losing much sensitive due to a relatively comprehensive protein database
and a liberal cutoff

=head1 USAGE

Once everything is installed properely, using halfwise should be quite simple. Go

 halfwise.pl

to get basic help and 

 halfwise.pl -help 

for more detailed options.

The basic mode of running is just to take a DNA sequence file and write to stdout, like

 halfwise.pl cosmid.fasta > cosmid.hlf

For more description of the options you can use with the genewisedb program (and this
includes a large number of output options) go 

  halfwise.pl -help

and also read the postscript documentation (wise2.ps) that came in this distribution

=head1 INSTALLATION

For halfwise to work you need a number of different software packages installed. This
means that halfwise is a very robust and working script (and didnt take long to write!)
but, you are going to have to hop around installing things first. Sorry!

Everything installs very easily! All the components required for halfwise can be found at
ftp://ftp.sanger.ac.uk/pub/birney/wise2/halfwise, so you can make a one stop shop!

You need

=over

=item Wise2.1.16 or higher level in the 2.1 release series

You should have it now, as you are reading this!

=item HMMer2.1.1 or higher

This is for the HMM retrieval system. Available from http://hmmer.wustl.edu or
the above ftp site (the HMMer2 directory in the Wise2 package does not have the
correct facilities).

=item bioperl central distribution, 0.04.2 or higher in the 0.04 release series

This is available from http://bio.perl.org or the above ftp site

=item The blast executables

available from http://ncbi.nlm.nih.gov or http://blast.wustl.edu/ (I did not make
these available from the ftp site as you probably have them already installed!)

=item Pfam and halfwise databases

available from the above ftp site

=back

For the software packages, follow the instructions for how to install them. Briefly

=over 

=item Wise2.1.16

type 'make all' in the wise2.1.16 directory and then copy the contents of the bin
directory to /usr/local/bin or somewhere on your path

=item HMMer2.1.1

type './configure' followed by 'make' in the hmmer-2.1.1 directory and then 'make install'
(assumming you want to install in /usr/local/bin)

=item bioperl-0.04.2

It is a standard CPAN type module. Go 'perl Makefile.PL', 'make', 'make install' in the
bioperl-0.04.2 directory. you do not need to install the Compiled extension, but if you 
want to use other parts of bioperl, you may need them, so there is no reason not to ;)

=back

For the databases, pick up the Pfam database and the halfwise.db. Put them somewhere where
you usually put databases. setdb or pressdb the halfwise.db ready for blast (depending
on what flavour of blast you are using - your database manager should know this). Also
the easiest thing to do is to setenv HMMERDB to the directory where the Pfam database is.
You then need to index the Pfam database going 'hmmindex Pfam'.

finally you need to adjust the variables at the top the halfwise script so that everything
works. 

This all seems complicated but its not so bad. Here is a check list.

  o pick up and install wise2.1.16
    - make all in wise2.1.16
    - cp bin/* /usr/local/bin
    - setenv WISECONFIGDIR to xxx/wise2.1.16/wisecfg

  o pick up and install hmmer2.1.1
    - ./configure
    - make
    - make install

  o pick up and install bioperl
    - perl Makefile.PL
    - make
    - make install

  o check you have blast installed

  o pick up halfwise.db and Pfam
    - setdb halfwise.db
    - hmmindex Pfam
    - setenv HMMERDB /somewhere/data/

  o change the location of the databases etc
    at the top of this script

If you have any problems, just email me <birney@sanger.ac.uk>

=cut
  
 

#
# Make sure halfwise is there...
#

$blast_db_dir = $ENV{'BLASTDB'};


if( ! -e $half_db && ! (defined $blast_db_dir && -e "$blast_db_dir/$half_db"  )) {
  die("The halfwise database [$half_db] has not been correctly indicated.\nPlease read the installation instructions at the top of this script\n");
}

$filename = shift; # dna sequence we hope!

if( !defined $filename ) {
    print "halfwise <dna-sequence-fasta> genewisedb-options\n";
    print "searches Pfam using initial blast followed by genewise\n";
    print "    The default genewise options are $defopts\n\n";
    print "Genewisedb options are listed by going halfwise -help\n";
    exit(1);
}

# if help - print the help out from genewisedb...

if( $filename =~ "-help" ) {
    print "halfwise <dna-sequence-fasta> genewisedb-options\n";
    print "    The default genewise options are $defopts\nGenewisedb options\n\n";
    open(GDB,"$genewise -help |");
    while(<GDB>) {
	print;
    }
    close(GDB) || die "Could not close pipe to genewisedb executable [$genewise]. Installation problem?";
    exit(1);
}

$verbose = 0;
	
if( @ARGV > 0 ) {
    $opts = join(' ',@ARGV);
} else {
    $opts = $defopts;
}

	   
$verbose && print STDOUT "Doing Blast...\n";
# make the blast output
# BLASTMAT=/nfs/disk100/pubseq/blastdb/

#$ENV{'BLASTMAT'} = "/nfs/disk100/pubseq/blastdb/";


if( system("blastx $half_db $filename -V=10000 -B=10000 > $filename.blx.$$") != 0 ) {
  die "There blastx search did not complete correctly\n";
}

$verbose && print STDOUT "Done\n";

eval {
    $blast = Bio::Tools::Blast->new( 
				    -file => "$filename.blx.$$",
				    -parse => 1,
				    -signif => '1e-3',
				    -strict => 1,
				    );
};
if( $@ ) {
    print "#No hits found in the blast search\n";
    exit;
}


# ok ... now loop over the HSPs and print them out!


foreach $hit ( $blast->hits ) {
    $name = $hit->name;
    $verbose && print STDOUT "Reading $name\n";
    if( !($name =~ /^([^\/]+)\//) ) {
	warn("Bad hit name $name");
    }
    $acc = $1;
    
    $hash{$acc} = 1;
}

# build minidb of the HMMs

open(TEMPDB,">$filename.dbhmm.$$");
$count = 0;
foreach $acc ( keys %hash ) {
    $count++;
    $verbose && print STDOUT "Loading $acc\n";
    open(GETZ,"$hmmfetch $hmmdb $acc |");
    while(<GETZ>) {
	print TEMPDB $_;
    }
    close(GETZ);
}
close(TEMPDB);
$verbose && print STDOUT "Collected $count hmms\n";

# run genewisedb

open(GDB,"$genewise -pfam $filename.dbhmm.$$ -dnas $filename $opts |") || die "Could not open genewise - ERROR! $!";

while(<GDB>) {
    print;
}
close(GDB) || die "Could not open genewise pipe (well close it anyway) $! $?";

unlink("$filename.dbhmm.$$");
unlink("$filename.blx.$$");

exit(0); # for james...





