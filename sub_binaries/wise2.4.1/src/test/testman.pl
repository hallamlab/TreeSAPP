#!/usr/local/bin/perl

if( $opt eq 'cf' ) {
    $conf = shift;
    $ENV{'WISECONFIGDIR'} = $conf;
}


#
# Reads in test files with test name, do and diff directives for each unit
#

@testarr = &read_all_testset(\*STDIN);
open(RES,">test.results");

$error = &do_all_testset(\@testarr,\*RES,\*STDERR);

if( $error == 1 ) {
    print STDERR "Tests failed\nPLEASE: email ewan birney (birney\@sanger.ac.uk) with the test.results file\n\n";
} else {
    print STDERR "Tests passed\n";
}




sub do_all_testset {
    my $testarr = shift;
    my $results_file = shift;
    my $talk_file = shift;
    my @array;
    my $error = 0;

    @array = @{$testarr};

    foreach $test (@array) {
	if( &do_testset($test,$results_file,$talk_file) == 1 ) {
	    $error = 1;
	}
    }

    return $error;
}

sub do_testset {
    my $test = shift;
    my $result_file = shift;
    my $talk_file = shift;
    my ($name,$do,$diff,$desc);
    my	$error =0;

    $name  = $test->{'name'};
    $do  = $test->{'do'};
    $do = "../bin/$do";
    $diff  = $test->{'diff'};
    $desc  = $test->{'desc'};

    if( $talk_file ) {
	print $talk_file ">>> Doing $name $desc\n>>> Calling [$do]\n";
    }

    if( ($exit = system($do)) != 0 ) {
	print $talk_file ">>> Bad exit status - $exit (system - $!)\n";
	$error =1;
    } else {

	open(DIFF,"diff $diff|");
	while(<DIFF>) {
	    chop;
	    if( /^[><].*wise.*\d/ ) {
		next;
	    }
	    if( /^[><].*\$[nN]ame\$/ ) {
		next;
	    }
	    if( /unreleased/ ) {
		next;
	    }

	    if( /^[<>]/ ) {
		print $result_file "Error in $name: $_\n";
		$error =1;
	    }
	}
    }
	
    if( $talk_file ) {
	if( $error ) {
	    print $talk_file ">>> $name ... Failed\n";
	} else {
	    print $talk_file ">>> $name ... Passed\n";
	}
    }

    return $error;
    
}
    
    

sub read_all_testset {
    my $file = shift;
    my @array;
    my $h;

    while( ($h = &read_testset($file)) ) {
	push(@array,$h);
    }

    return @array;
}

sub read_testset {
    my $file = shift;
    my $hash = {};
    my ($name,$do,$diff,$desc);

    
    while(<$file>) {
	chop;
	/^#/ && next;
	/^\/\// && last;
	/^name\s+(\S+)/ && do { $name = $1; next; };
	/^desc\s+(.*)$/ && do { $desc = $1; next; };
	/^do\s+(.*)$/ && do { $do = $1; next; };
	/^diff\s+(\S+\s+\S+)/ && do { $diff = $1; next; };
	/^\s+$/ && next;
	warn("Could not understand $_");
    }

    if( defined $name && defined $do && defined $diff) {
	$hash->{'name'} = $name;
	$hash->{'do'} = $do;
	$hash->{'diff'} = $diff;
	$hash->{'desc'} = $desc;

	return $hash;
    } else {
	return undef;
    }
}


	
	

