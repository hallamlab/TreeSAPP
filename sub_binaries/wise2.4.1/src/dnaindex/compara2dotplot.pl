use strict;

my $shift_x = 20;
my $shift_y = 20;

my $length_x = 500;
my $length_y = 500;

my $max_x = 200_000_000;
my $max_y = 200_000_000;

my $scale = 10_000_000;

my $slice = 1;
my $slice_x_start = 16_900_000;
my $slice_y_start = 17_400_000;
my $slice_x_end   = 17_600_000;



if( $slice == 1 ) {
  $max_x = ($slice_x_end - $slice_x_start);
  $max_y = $max_x;
}

print <<EOF;
%!PS-Adobe-2.0
% Created by compara2dotplot. Ewan Birney <birney\@ebi.ac.uk>
/Helvetica findfont 12 scalefont setfont
0.25 setlinewidth
/Helvetica findfont 8 scalefont setfont
EOF

#if( $rotate ) {
#   print "0 700 translate\n";
#}

# grid lines.

print "$shift_x $shift_y moveto\n";
print $shift_x + $length_x," $shift_y lineto stroke\n";

print "$shift_x $shift_y moveto\n";
print "$shift_x ",$shift_y+$length_y," lineto stroke\n";


my $i =0;

if( $slice == 0 ) {
  while( $i < $max_x ) {
    my $xpos = $shift_x + ($i*$length_x / $max_x);
    print $xpos," $shift_y moveto\n";
    print $xpos," ",$shift_y-2," lineto stroke\n";
    print $xpos," ",$shift_y-10," moveto\n";
    my $number = $i / 10_000_000;
    print "($number) show\n";
    $i += $scale;
  }
  
  $i =0;
  while( $i < $max_y ) {
    my $ypos = $shift_y + ($i*$length_y / $max_y);
    print "$shift_x $ypos moveto\n";
    print $shift_x-2," $ypos lineto stroke\n";
    print $shift_x-13," $ypos moveto\n";
    my $number = $i / 10_000_000;
    print "($number) show\n";
    $i += $scale;
  }
}


while( <> ) {
  /chr/ || next;
  #1024 chr7	0	1024	chr7	171263	170239	-
  my ($len,$a,$startx,$endx,$b,$starty,$endy) = split;
  
  if( $len < 500 ) {
    next;
  }

  if( $slice ) {
    if( $startx > $slice_x_end || $endx < $slice_x_start ) {
      next;
    }

    $startx -= $slice_x_start;
    $endx -= $slice_x_start;
  
    $starty -= $slice_y_start;
    $endy -=  $slice_y_start;
  }

  my $coord_x_start = $shift_x + ($startx * $length_x / $max_x);
  my $coord_x_end   = $shift_x + ($endx * $length_x / $max_x);


  my $coord_y_start = $shift_x + ($starty * $length_y / $max_y);
  my $coord_y_end   = $shift_y + ($endy * $length_y / $max_y);

  
  print "$coord_x_start $coord_y_start moveto\n";
  print "$coord_x_end   $coord_y_end   lineto\n";
  print "stroke\n";

}

print "showpage\n";
