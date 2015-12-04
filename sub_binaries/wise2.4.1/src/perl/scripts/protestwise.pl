#!/usr/local/bin/perl

#
# protestwise.pl <protein-seq-fasta> <dna-seq-fasta>\n
# produces on STDOUT a new protein sequence which is the 
# DNA sequence 'fixed' by the comparison to the protein sequence.

# in particular frameshift errors get mapped to X

# written by James cuff (james@ebi.ac.uk)
# Hacked by Ewan (birney@sanger.ac.uk). Talk to ewan
# first about the script.


use Wise2;

my $pro_file  = shift; # first argument from @ARGV
my $dna_file  = shift; # second argument @ARGV

if( !defined $dna_file ) {
    die "ProtESTwise.pl <protein-seq-fasta> <dna-seq>\nProduces output of the DNA sequence\n'fixed' by the comparison on stdout";
}

# read in inputs. Read in first as generic 'Sequence' objects
# and then converted to specific 'Protein' or 'cdna' type
# objects
 
open(PRO,$pro_file) || die "Could not open $pro_file!";
$seq = &Wise2::Sequence::read_fasta_Sequence(\*PRO);
$pro = &Wise2::Protein::Protein_from_Sequence($seq);

if( $pro == 0 ) {
    # can't interpolate function calls <sigh>
    die sprintf("Could not make protein from sequence %s!",$seq->name());
}

open(DNA,$dna_file) || die "Could not open $pro_file!";
$seq = &Wise2::Sequence::read_fasta_Sequence(\*DNA);
$cdna = &Wise2::cDNA::cDNA_from_Sequence($seq);

if( $cdna == 0 ) {
    # can't interpolate function calls <sigh>
    die sprintf("Could not genomic from sequence %s!",$seq->name());
}

# Read in data structures needed for 
# estwise type algorthim
#
# These will be automatically read from WISECONFIGDIR if necessary.

# this is the indel rate
$cp  = &Wise2::flat_cDNAParser(0.001);

# codon table
$ct  = &Wise2::CodonTable::read_CodonTable_file("codon.table");

#this means we are not using any codon bias
$cm  = &Wise2::flat_CodonMapper($ct);

# this is the substitution error
$cm->sprinkle_errors_over_CodonMapper(0.001);

# random model needed if we are not using syn
$rmd = &Wise2::RandomModelDNA_std();

# means estwise3 algorithm. Not obvious!
$alg = 0;

# sets memory amount for main memory
&Wise2::change_max_BaseMatrix_kbytes(100000); # 10 Megabytes.

# these are for the protein part of the comparison
$comp   = &Wise2::CompMat::read_Blast_file_CompMat("blosum62.bla");
$rm     = &Wise2::default_RandomModel();

# do it!
$alb = &Wise2::AlnBlock_from_Protein_estwise_wrap($pro,$cdna,$cp,$cm,$ct,$comp,-12,-2,0,$rmd,$alg,$rm,1,0.9);
$proseq = "";


#
# This is where we get clever!
#
# The for loops across the alignments. The protein sequence is in $alc->alu(0). The
# DNA sequence is in $alc->alu(1). We are interested in codons in the DNA sequence
# and turns those into amino acids. Sequence insertions or deletions become X's
#

for($alc=$alb->start();$alc->at_end() != 1;$alc = $alc->next()) {
    


    if( $alc->alu(1)->text_label() =~ /^INSERT$/ ) {
	next; # skip protein inserts relative to the DNA sequence
	# NB different from SEQUENCE_INSERTION.
    }

    if( $alc->alu(1)->text_label() =~ /CODON/ ) {	
	# get out sequence from $start to $end
	# $start and $end are in bio coordinates
	$start = $alc->alu(1)->start+1; 
	$end   = $alc->alu(1)->end+1;                
	$dnatemp = "";

	for($x=$start;$x < $end;$x++){
	    $tmp = &Wise2::cDNA::cDNA_seqchar($cdna,$x);
	    $dnatemp=$dnatemp.$tmp;
	}	

	
	$temp = $ct->aminoacid_from_seq($dnatemp);

	# if codon has an N, then set the residue to unk X,
	# we could be clever about this and work out what 
	# it is likely to be, but hell...

	$temp =~ s/x/X/;

	$proseq .= $temp;
    } else {
	# deletion or insertion of a base
	$proseq .= 'X';
    }

    
}

# make the new protein sequence and 
# dump it to stdout

$namecdna = $cdna->baseseq()->name();
$new = &Wise2::new_Sequence_from_strings($namecdna,$proseq);
$new->write_fasta(STDOUT);
