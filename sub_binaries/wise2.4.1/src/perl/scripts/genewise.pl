#!/usr/local/bin/perl

use Wise2;

my $pro_file  = shift; # first argument from @ARGV
my $dna_file  = shift; # second argument @ARGV

if( !defined $dna_file ) {
    die "genewise.pl <protein-seq-fasta> <dna-seq>\n";
}

#
# read in inputs. Read in first as generic 'Sequence' objects
# and then converted to specific 'Protein' or 'Genomic' type
# objects
# 

open(PRO,$pro_file) || die "Could not open $pro_file!";
$seq = &Wise2::Sequence::read_fasta_Sequence(\*PRO);
$pro = &Wise2::Protein::Protein_from_Sequence($seq);

if( $pro == 0 ) {
        # can't interpolate function calls <sigh>
        die sprintf("Could not make protein from sequence %s!",$seq->name());
}

open(DNA,$dna_file) || die "Could not open $pro_file!";
$seq = &Wise2::Sequence::read_fasta_Sequence(\*DNA);
$gen = &Wise2::Genomic_from_Sequence($seq);
	
if( $gen == 0 ) {
        # can't interpolate function calls <sigh>
        die sprintf("Could not genomic from sequence %s!",$seq->name());
}

#
# Read in data structures needed for 
# genewise type algorthim
#
# These will be automatically read from WISECONFIGDIR if necessary.
#
#

$gf   = &Wise2::read_GeneFrequency21_file("human.gf");
$ct   = &Wise2::CodonTable::read_CodonTable_file("codon.table");
$rmd  = &Wise2::RandomModelDNA_std();
$alg  = &Wise2::gwrap_alg_type_from_string("623");
$gp   = &Wise2::GeneParameter21_wrap($gf,0.00001,0.000001,$rmd,0,1,1,$ct,0.999,0.99,0.8,0.99,0.4);

# sets memory amount for main memory

&Wise2::change_max_BaseMatrix_kbytes(10000); # 10 Megabytes.


#
# these are for the protein part of the comparison
#

$cm   = &Wise2::CompMat::read_Blast_file_CompMat("blosum62.bla");
$rm   = &Wise2::default_RandomModel();

# do it!

$pg = new Wise2::PotentialGene; # empty potential gene

$alb = &Wise2::AlnBlock_from_protein_genewise_wrap($pro,$gen,$cm,-12,-2,$gp,$rmd,$rmd,$alg,0,1,$rm,0.9,$pg);

&Wise2::protein2genomic_ascii_display($alb,$pro,$gen,$ct,15,50,STDOUT);




