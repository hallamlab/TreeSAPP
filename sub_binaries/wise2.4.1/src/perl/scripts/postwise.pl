#!/usr/local/bin/perl


#
# Skeleton script to emulate postwise
#   written by Ewan - expecting major hacking from other people ;)
#
#

use Wise2;

#
# By intialising objects to 0, the Wise2 libraries have 
# an outside chance of catching unintialised variables.
# also delibrately 'non-blesses' them, making them uncovertable
# for some functions...
#

my $pgl = 0;  # potential gene list
my $gen = 0;  # genomic object
my $genr= 0; # genomic region object

my $gf= 0;   # gene frequencies
my $comp=0; # comparison matrix

my $gpara=0; # holder for all gene parameters
my $rmd=0;   # dna random model
my $rm=0;    # protein random model

my $pro_file  = shift;
my $dna_file  = shift;

if( !defined $dna_file ) {
    die "postwise.pl <protein-seq-fasta> <dna-seq>\n";
}

#
# read in inputs
#

$gen = &Wise2::Genomic::read_fasta_file_Genomic($dna_file);

if( $gen == 0 ) {
    die "Could not read $dna_file as a Genomic Sequence";
}

# if you want to get at the sequence, here it is...

# $seq = $gen->baseseq();
# $seq->write_fasta(STDOUT);

#
# Allocate an empty PotentialGeneList
#

$pgl = new Wise2::PotentialGeneList;

open(PRO,$pro_file) || die "Could not open $pro_file!";

while (1) {
    $seq = &Wise2::Sequence::read_fasta_Sequence(\*PRO);
    if( $seq == 0 ) {
	last; # end of file
    }

    $pro = &Wise2::Protein::Protein_from_Sequence($seq);

    if( $pro == 0 ) {
	# can't interpolate function calls <sigh>
	die sprintf("Could not make protein from sequence %s!",seq->name());
    }
    $pg  = new Wise2::PotentialGene;
    
    $pg->set_homolog($pro); # sets this protein as a homolog to be used
    

    #
    # This should be set - obviously - to where you think the gene should start/end
    # if end < start then it is on the opposite strand.
    #

    $pg->set_guess_start(1);
    
    # good OOP... but - descending 3 func calls to get out
    # a simple attribute. *don't* do this in major loops
    
    $pg->set_guess_end ($gen->baseseq()->end()); 

    #
    # add this to potential gene list
    #
    
    $pgl->add_pg($pg);

}

#
# these are for the gene model
#

#
# We should have a short cut for this in Wise2.pm probably.
#

$gf   = &Wise2::GeneFrequency21::read_GeneFrequency21_file("human.gf");
$ct   = &Wise2::CodonTable::read_CodonTable_file("codon.table");
$rmd  = &Wise2::RandomModelDNA::RandomModelDNA_std();
$alg  = &Wise2::gwrap_alg_type_from_string("623");
$gp   = &Wise2::GeneParameter21_wrap($gf,0.00001,0.000001,$rmd,0,1,1,$ct,0.999,0.99,0.8,0.99,0.4);

#destroy memory held by Gene frequency object

# sets memory amount for main memory

&Wise2::change_max_BaseMatrix_kbytes(10000); # 10 Megabytes.

$gf = 0;

#
# these are for the protein part of the comparison
#

$cm   = &Wise2::CompMat::read_Blast_file_CompMat("blosum62.bla");
$rm   = &Wise2::RandomModel::default_RandomModel();

#
# Make a new genomic region from the Genomic object
# GenomicRegion is like a 'GenomeSequence' (vaguely) in ACeDB - it is
# our outer-layer context of a region of genome, hopefully containing some
# sequence ;)
#

$gr   = &Wise2::GenomicRegion::new_GenomicRegion($gen);

# print STDERR "Going to resolve! $gr $g2 $gen $pgl $alg $comp $gp $rmd\n";

#
# resolve those genes!
#
# This is the major call -> into the gwrap package...
#

$count = &Wise2::resolve_PotentialGenes_on_GenomicRegion($gr,$pgl,$alg,$alg,0,0,$cm,-12,-2,$gp,$rmd,$rmd," ",0,20);

#
# shows as ascii text
#

$gr->show_pretty_GenomicRegion(STDOUT);

#
# shows as ACeDB objects
#

$gr->show_ace_GenomicRegion($gen->baseseq->name,STDOUT);

#
# Gets out the alignments. Shows them.
#

foreach $pg ( $pgl->each_pg() ) {
    
    $alb = $pg->alb();
    $pro = $pg->homolog();


    if( $alb == 0 ) {
	warn sprintf("At potential gene %s - no alignment",$pg->homolog()->baseseq()->name());
	next; 
    }


    $gtemp = $gr->genomic();
    $gtemp = $gtemp->truncate_Genomic($pg->guess_start,$pg->guess_end);

    &Wise2::protein2genomic_ascii_display($alb,$pro,$gtemp,$ct,15,50,STDOUT);

    print "\n\n//\n";
}

#
# Finished!
# 


