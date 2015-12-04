#!/usr/local/bin/perl -w


=head1 NAME 

blastwise

=head1 SYNOPSIS

blastwise swiss myseq.fa

Blasts your fasta dna sequence against swissprot, reads in the 
results and then uses the GeneWise algorithm to make
matches in the correct place

=head1 OPTIONS

To get other output options, please edit this script to your own
needs. This script is really just to seed the process for your own
site: I cant second guess all the things you want to do with the
output, so you give it a go.

This assummes that the database you are searching *is* swissprot and
you have SRS indexed it into the database called swissprot.  If you
are running it on different database with a different indexer you will
need to change the executable used to retrieve the sequence

=cut

use strict;
no strict 'subs';

my $db  = shift;
my $seqfile = shift;

if( !defined $seqfile ) {
    system("pod2text $0");
    exit(1);
}

my $BLASTEXE = 'blastx';     # might be different on your system of course!
my $TEMPFILE = "temp.bw"; # feel free to change this!
my $genefile = "human.gf";

# for debugging @sanger centre

use lib '/nfs/disk92/PerlSource/Bioperl/Releases/bioperl-0.05/';

BEGIN {
    eval {
	require Wise2;
    };
    if ( $@ ) {
	print STDERR ("\nThe C-compiled engine for Wise2 has not been installed.\n Please read the installation instructions for Wise2 for building it\n\n - bascially make perl at the top level\n\n");
	exit(1);
    }
    eval {
	require Bio::Tools::Blast;
    };
    if ( $@ ) {
	print STDERR ("\nThe bioperl Blast object is not installed. Please install it from http://bio.perl.org/\n");
	exit(1);
    }
}


my ($blast);

# make genewise parameter objects
# before we try to make blast

my $gf   = &Wise2::read_GeneFrequency21_file($genefile) || die "Could not read $genefile for parameters" ;
my $ct   = &Wise2::CodonTable::read_CodonTable_file("codon.table") || die "Could not read codon table";
my $rmd  = &Wise2::RandomModelDNA_std() || die "Could not make standard random model (weird!)";
my $alg  = &Wise2::gwrap_alg_type_from_string("623") || die "Could not make algorithm";
my $gp   = &Wise2::GeneParameter21_wrap($gf,0.00001,0.000001,$rmd,0,1,1,$ct,0.999,0.99,0.8,0.99,0.4) || die "Could not make gene parameters";

# sets memory amount for main memory

&Wise2::change_max_BaseMatrix_kbytes(100000); # 10 Megabytes.


#
# these are for the protein part of the comparison
#

my $cm   = &Wise2::CompMat::read_Blast_file_CompMat("blosum62.bla");
my $rm   = &Wise2::default_RandomModel();

#
# read in the DNA sequence before anything goes wrong!
#

open(DNA,$seqfile) || die "Could not open $seqfile!";
my $seq = &Wise2::Sequence::read_fasta_Sequence(\*DNA);
my $gen = &Wise2::Genomic_from_Sequence($seq);
close(DNA);
	
if( $gen == 0 ) {
        # can't interpolate function calls <sigh>
        die sprintf("Could not genomic from sequence %s!",$seq->name());
}


# we try to run blast
if( ! -e $TEMPFILE ) {
    if ( system("$BLASTEXE $db $seqfile > $TEMPFILE") != 0 ) {
	die("Could not run blast [$BLASTEXE $db $seqfile > $TEMPFILE] $! $?");
    }
} else {
    print STDERR "Skipping making blast";
}

# read it into Blast object. We have to catch an exception if there
# is nothing on old versions of bioperl. 

eval {
    $blast = Bio::Tools::Blast->new( 
				    -file => "$TEMPFILE",
				    -parse => 1,
				    -signif => '1e-3',
				    -strict => 1,
				    );
};
if( $@ ) {
    print "#No hits found in the blast search\n";
    exit;
}

# get rid of TEMPFILE now.

# comment this out for debugging.

unlink($TEMPFILE);


#
# Ok. Ready to rock and roll. Loop over blast hits. Take the hit
# out - ask whether its highest HSP is already in a region which has a gene. 
#
# If not, build genewise object.
#

#
# This is the data structure which we will store all the genes in as we make 
# them
#

my $mother = &Wise2::GenomicRegion::new_GenomicRegion($gen);


my ($make_aln,$hit,@hsps,$hsp,$pgl,$pg,$name,@genes,$pro);


#
# This is the main loop. Each time around it selects a single protein
# to add into the 
#
#
#

$make_aln = 0;
foreach $hit ( $blast->hits ) {
    $name = $hit->name();
    @hsps = sort { $b->score <=> $a->score } $hit->hsps();
    $hsp = shift @hsps;

    if( $make_aln > 10 ) {
	last;
    }

    $make_aln++;

    if( &covered_by_genes($hsp,@genes) ) {
	next;
    }


    # ok - run genewise
    
    # we need the protein...
    # Here is the pipe you should change if you are not using getz
    
    open(PRO,"getz -sf fasta -d '[swissprot-id:$name]' |") || die "Could not fork getz pipe";
    $seq = &Wise2::Sequence::read_fasta_Sequence(\*PRO);
    close(PRO) || die "Could not open pipe to getz with [getz -sf fasta -d '[swissprot-id:$name]' |].\n Perhaps you want a different indexing executable?";
    
    if( $seq ) {
	$pro = &Wise2::Protein::Protein_from_Sequence($seq);
    } else {
	warn("Could not read $name");
	next;
    }


    $pgl = new Wise2::PotentialGeneList; # empty list.

    $pg = new Wise2::PotentialGene; # empty potential gene
    $pgl->add_pg($pg);
    $pg->set_homolog($pro);

    #
    # These should be changed to be the start/end point 
    # in the DNA sequence.
    #

    @hsps = sort { $a->start('query') <=> $b->start('query') } $hit->hsps();
    $hsp = shift @hsps;
    my $dstart = $hsp->start('query');

    if( ($dstart - 500) > 0 ) {
	$dstart -= 500;
    } else {
	$dstart = 0;
    }

    @hsps = sort { $b->end('query') <=> $a->end('query') } $hit->hsps();
    $hsp = shift @hsps;
    my $dend = $hsp->end('query');

    if( $dend + 500 < $gen->baseseq->len ) {
	$dend += 500;
    } else {
	$dend = $gen->baseseq->len;
    }

    if( $hsp->strand('query') eq 'Plus' ) {
	$pg->set_guess_start($dstart);
	$pg->set_guess_end($dend);
    } else {
	$pg->set_guess_end($dstart);
	$pg->set_guess_start($dend);
    }

    &Wise2::resolve_PotentialGenes_on_GenomicRegion($mother,$pgl,$alg,$alg,20,20,$cm,-12,-2,$gp,$rmd,$rmd," ",0,25);

    # transfer genes into hash.
    @genes = ();
    foreach my $gene ( $mother->each_gene() ) {
	push(@genes,$gene);
    }

    # get out alignment

    my $alb = $pgl->pg(0)->alb;

    # alignment coordinates are relative to the truncated genomic DNA 
    # sequence, not the absolute. So we need to truncate the genomic
    # DNA before we put it into the system

    # print STDERR "Truncating ", $pgl->pg(0)->guess_start, ":", $pgl->pg(0)->guess_end, "\n";


    # mystic -1's to map from bio coordinates to C coordinates
    my $gent;
    if( $pgl->pg(0)->guess_start < $pgl->pg(0)->guess_end ) {
	$gent = $gen->truncate_Genomic($pgl->pg(0)->guess_start-1,$pgl->pg(0)->guess_end);
    } else {
	$gent = $gen->truncate_Genomic($pgl->pg(0)->guess_start,$pgl->pg(0)->guess_end-1);
    }

    &Wise2::protein2genomic_ascii_display($alb,$pro,$gent,$ct,15,50,STDOUT);

}

#
# show all the genes?
#

$mother->show_ace_GenomicRegion($gen->baseseq->name(),STDOUT);




#### sub routines ####

#
# Determines is this hsp is actually already covered.
#

sub covered_by_genes {
    my $hsp = shift;
    my (@genes) = @_;
    my ($start,$end,$gene,$strand,$mid);

    ($start,$end) = $hsp->range('query');
    ($strand)     = $hsp->strand('query');

    $mid = $start + ($end - $start)/2;

    foreach $gene ( @genes ) {
	if( $gene->start > $gene->end ) {
	    if( $strand eq 'Plus' ) {
		next;
	    }
	} else {
	    if( $strand eq 'Minus' ) {
		next;
	    }
	}

	# now look for internal overlap

	if( $strand eq 'Plus' ) {
	    if( $mid > $gene->start && $mid < $gene->end ) {
		# is covered
		return 1;
	    } 
	} else {
	    if( $mid > $gene->end && $mid < $gene->start ) {
		# is covered
		return 1;
	    } 
	}
    }

    # not covered by any gene - get out!
    
    return 0;
}




