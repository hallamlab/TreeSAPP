

use Bio::SeqIO;


my $seqin = Bio::SeqIO->new( -format => 'fasta');
my $seqout = Bio::SeqIO->new( -file => '>-' );


my $seq_size = shift;

if( ! defined $seq_size ) {
  print STDERR "no read size set, going to use 500\n";
  $seq_size = 500;
}

my $coverage = shift;
if( !defined $coverage ) {
  $coverage = 5;
}


my $seq = $seqin->next_seq();

my $error = 0;

my $rev = 1;

my $insert = 1;

my $del = 1;

my @alphabet = ('A' , 'T' , 'G' , 'C' );

my $reads = $seq->length * 5 / $seq_size;

foreach my $i ( 1 .. $reads ) {
  $start = int(rand($seq->length-$seq_size));
  my $subseq = $seq->trunc($start+1,$start+$seq_size);
  $subseq->id(($seq->id.".".$start.".".($start+$seq_size)));


  if( $error != 0 ) {
      my $seqstr = $subseq->seq();
      foreach my $j ( 1 ..$seq_size ) {
	  if( $random == 1 && int(rand(1000)) < 2 ) {
	      my $new  = $alphabet[int(rand(4))];
	      print STDERR "[$i] Randomising a position $j... ",substr($seqstr,$j,1)," to $new on ".$subseq->id,"\n";
	      substr($seqstr,$j,1) = $new;
	  }

      }
      $subseq->seq($seqstr);
  }

  if( $rev != 0 ) {
      if( int(rand(10)) < 5 ) {
	  $seq = $seq->revcom();
      }
  }

  $seqout->write_seq($subseq);
}
