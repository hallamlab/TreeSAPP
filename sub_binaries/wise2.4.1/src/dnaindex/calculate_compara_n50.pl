
my %len = (
	   '1' =>  210989981,
	   '10' =>  127014353,
	   '10_random' =>  593526,
	   '11' =>  129303236,
	   '11_random' =>  0,
	   '12' =>  125990086,
	   '12_random' =>  294697,
	   '13' =>  93269031,
	   '13_random' =>  262883,
	   '14' =>  86368090,
	   '14_random' =>  418315,
	   '15' =>  74363214,
	   '15_random' =>  118395,
	   '16' =>  73084034,
	   '16_random' =>  1338419,
	   '17' =>  73859773,
	   '17_random' =>  219907,
	   '18' =>  72035922,
	   '18_random' =>  316662,
	   '19' =>  53798088,
	   '19_random' =>  0,
	   '1_random' =>  2549662,
	   '2' =>  221356626,
	   '20' =>  61555802,
	   '21' =>  33823978,
	   '22' =>  33785749,
	   '2_random' =>  735116,
	   '3' =>  185824536,
	   '3_random' =>  582736,
	   '4' =>  169058312,
	   '4_random' =>  1250162,
	   '5' =>  166043302,
	   '5_random' =>  1959775,
	   '6' =>  158660857,
	   '7' =>  147591730,
	   '7_random' =>  605209,
	   '8' =>  124912316,
	   '8_random' =>  0,
	   '9' =>  109839313,
	   '9_random' =>  643757,
	   'NA_random' =>  7926661,
	   'UL_random' =>  11960799,
	   'X' =>  124872836,
	   'X_random' =>  5188694,
	   'Y' =>  21805613,
);



use strict;

my $total = 0;
my @lengths;
my $q;

while( <> ) {
  /chr/ || next;
  #1024 chr7	0	1024	chr7	171263	170239	-
  my ($len,$a,$startx,$endx,$b,$starty,$endy) = split;

  if( $total < $endx ) {
    $total = $endx;
  }

  $q = $a;
  push(@lengths,$len);
}

$q =~ s/query_chr//;

$total = $len{$q};

#$total = 3_000_000_000;

@lengths = sort { $b <=> $a } @lengths;
my $top = $lengths[0];

my $run = 0;
foreach my $len ( @lengths ) {
  if( $run+$len > $total / 2 ) {
    print "$q N50 $len LONGEST $top LENGTH $total\n";
    last;
  }
  $run += $len;
  #print STDERR "got $run with $len\n";
}
  
