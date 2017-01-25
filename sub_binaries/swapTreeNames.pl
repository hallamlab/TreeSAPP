#!/usr/bin/perl -w
#swapTreeNames.pl by Young C. Song (Written on July 22nd, 2009)

use strict;
use Getopt::Long;

my $result;
my $TREE;
my $LIST;
my $OUTPUT;

if(defined $ARGV[0]) {
    $result=GetOptions(
   "tree=s"=>\$TREE,
   "list=s"=>\$LIST,
   "output=s"=>\$OUTPUT
   );
}

if(!$result) {
   die("Usage: $0 [-t tree_file] [-l list_file] [-o output_file]
   
   This program reads TREE file and LIST file, and changes the name of the sequences in your
   TREE file.  The LIST file should be consisting of two columns, each separated by tab character.
   The first column should contain the names of sequences as shown in your TREE file.  The second
   column should contain the names of sequences, which will replace those in the first column.
   
     E.g.  First Col	Second Col
           DQ521862	    GOM clone SMI1-GC205-mcr50
           AB176928	    Kuroshima clone Kuro-mcrA-4.02
           AJ937707	    PMMV clone PMMV-mcrA377
           CP000099	    Methanosarcina barkeri str. Fusaro
           
   NOTE: Please avoid using round brackets in your names in the second column (square brackets
   are fine)\n"
      );
}#end if

my %nameHash;
my @names = ();

open(LIST, "<$LIST") or die("$LIST: $!");

while(<LIST>) {

   chomp;
   
   my $listLine = $_;
   
   if($listLine =~ /(.+)\t(.+)/){
      @names = split("\t", $listLine);
      $nameHash{$names[0]} = $names[1];
   }#end if
}#end while
   
close LIST;

open(TREEFILE, "<$TREE") or die("$TREE: $!");;
open(OUTFILE, ">$OUTPUT") or die("$OUTPUT: $!");;

while(<TREEFILE>) {

   chomp;
   my $treeLine = $_;
   
   foreach my $firstName(sort(keys %nameHash)) {
      
      if($treeLine =~ m/$firstName/) {
         $treeLine =~ s/$firstName/$nameHash{$firstName}/;   
      }#end if
   }#end foreach $arbName
   print OUTFILE $treeLine."\n"; 
}#end while

close OUTFILE;
close TREEFILE;
