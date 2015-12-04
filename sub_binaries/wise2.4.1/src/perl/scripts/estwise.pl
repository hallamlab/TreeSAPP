#!/usr/local/bin/perl

=head1 NAME

estwise.pl

=head1 SYNOPSIS

estwise.pl <protein-seq-fasta> <dna-db-fasta>

=head1 DESCRIPTION

compares a protein sequence to a DNA database, allowing for
frame shift errors

=head1 OPTIONS

This perl script is just an example. For more options you
should read the perl script and customise it for your own
applications

=cut

use Wise2;
use strict;
# so we can use STDOUT as a bare word.
no strict 'subs';

my $pro_file = shift;
my $dna_file = shift;

if( !defined $dna_file ) {
    system("perldoc $0");
    exit(1);
}

# read in the protein file

open(PRO,$pro_file) || die "Could not open $pro_file!";
my $seq = &Wise2::Sequence::read_fasta_Sequence(PRO);
my $pro  = &Wise2::Protein::Protein_from_Sequence($seq);

# Convert the protein to a profile-HMM (called a threestatemodel in wise2 speak)

my $comp = &Wise2::CompMat::read_Blast_file_CompMat("blosum62.bla");
my $rm   = &Wise2::default_RandomModel();
my $tsm  = &Wise2::ThreeStateModel::ThreeStateModel_from_half_bit_Sequence($pro,$comp,$rm,-12,-2);

# make a profile-HMM database from the single profile-HMM

my $tdb  = &Wise2::ThreeStateDB::new_single_ThreeStateDB($tsm,$rm);

# destroy tsm memory...

$tsm = 0;

# load up dna file for as a simple sequence db

my $sdb = &Wise2::single_fasta_SequenceDB($dna_file);
if( $sdb == 0 ) {
    die "Could not read in $dna_file!";
}

# convert to a typed (cDNA) db

my $cdb = &Wise2::new_cDNADB($sdb);

# load parameters to the estwise search

# this is the indel rate
my $cp  = &Wise2::flat_cDNAParser(0.001);

# codon table
my $ct  = &Wise2::CodonTable::read_CodonTable_file("codon.table");

#this means we are not using any codon bias

my $cm  = &Wise2::flat_CodonMapper($ct);

# this is the substitution error

$cm->sprinkle_errors_over_CodonMapper(0.001);

# random model needed if we are not using syn

my $rmd = &Wise2::RandomModelDNA_std();

# means estwise3 algorithm

my $alg = &Wise2::alg_estwrap_from_string("333");

my $dbsi = &Wise2::new_serial_DBSearchImpl();

#print STDERR "Launching ($tdb,$cdb,$cp,$cm,$rmd,1,$alg,0) \n";

my $hs = &Wise2::Hscore_from_TSM_estwise($tdb,$cdb,$cp,$cm,$rmd,1,$alg,0,0.9,-1,1,$dbsi);

# simple basic show
#$hs->show(STDOUT);

$hs->sort_Hscore_by_score();


#
# This is how to loop over hscore object. Notice 
# that we don't make a DataScore object until we really have
# to - this is because the memory managment in hscore is
# customised - meaning retrieval of a DataScore means an entire
# allocation and copy of contents (usually it is just up'ing a reference count)
#
# Do not EVER loop over hs->datascore(...) putting everything into 
# an array (for a nice foreach call). It will kill your machine for
# large databases. 
#

my $len = $hs->length();

my($i,$ds,$alb,$cdna);
for($i=0;$i<$len;$i++) {

    if( &Wise2::Score2Bits($hs->score($i)) < 25 ) {
	last;
    }

    $ds = $hs->datascore($i);
    print STDERR sprintf("$i %12s %12s [%s] %4.2f\n",$ds->query()->name(),$ds->target()->name,$ds->target->is_reversed == 1 ? "-" : "+",&Wise2::Score2Bits($ds->score()));

}


#
# Get out the sequence of the first ones up to 5
#
# Align them and report back
#

for($i=0;$i<$len && $i < 5;$i++) {

    if( &Wise2::Score2Bits($hs->score($i)) < 25 ) {
	last;
    }

    $ds = $hs->datascore($i);
    $cdna = $cdb->get_entry($ds->target);
    if( ! $cdna ) {
	warn("Can't get entry!");
	next;
    }

    $alb= &Wise2::AlnBlock_from_Protein_estwise_wrap($pro,$cdna,$cp,$cm,$ct,$comp,-12,-2,0,$rmd,0,$rm,0.9,1);
    &Wise2::protcdna_ascii_display($alb,$pro->baseseq->seq,$pro->baseseq->name,$pro->baseseq->offset,$cdna,$ct,15,50,0,STDOUT);

}

#explicitly reap memory - helps bug tracking!

$tsm = 0;
$hs = 0;
$ds = 0;
#$tdb =0;

#$seq =0;
$pro = 0;
$comp = 0;
$rm = 0;
$cp = 0;
$ct = 0;


print STDERR "Finished!\n";

