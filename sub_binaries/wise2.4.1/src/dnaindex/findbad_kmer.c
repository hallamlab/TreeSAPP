
#include "../models/version.h"
#include "kmer_direct.h"
#include "kmer_hash.h"
#include "kmer_glib_index.h"
#include "singleseqspace.h"

#include "kmer_count.h"

char * program_name = "findbad_kmer";


void show_version(FILE * ofp)
{
  fprintf(ofp,"%s\nVersion: %s\nReleased: %s\nCompiled: %s\n",program_name,VERSION_NUMBER,RELEASE_DAY,COMPILE_DATE);
  fprintf(ofp,"\nThis program is freely distributed under a Gnu Public License\n");
  fprintf(ofp,"The source code is copyright (c) EMBL and others\n");
  fprintf(ofp,"There is no warranty, implied or otherwise on the performance of this program\n");
  fprintf(ofp,"For more information read the GNULICENSE file in the distribution\n\n");
  fprintf(ofp,"Credits: Ewan Birney <birney@ebi.ac.uk>\n");
  exit(63);   
}

void show_help(FILE * ofp)
{
  fprintf(ofp,"%s fasta-file\n",program_name);
  fprintf(ofp,"Kmer options\n");
  fprintf(ofp,"  -min_count 4            minimum count for reporting\n");
  fprintf(ofp,"  -kmer [12]              kmer size for kmer-graph\n");
  fprintf(ofp,"  -[no]use_hash [default no] use 28bit hash for large kmer sizes\n");
  fprintf(ofp,"  -[no]use_glib [default no] use glib based pure hash, up to 15mers\n");

}





int main(int argc,char ** argv)
{
  KmerCountAllocator * kca;
  KmerCount * count;

  char buffer[200];

  KmerIndexInterface * kii;
  FILE * ifp;
  Sequence * seq;
  Sequence * tseq;
  char c;

  int kmer_size = 12;
  boolean use_28bithash = 0;
  boolean use_glib = 0;
  boolean show_stdout = TRUE;

  int min_count = 4;
 
  int rev;

  SinglePosSpace * sps;

  
  kmer_t next_number;
  long int i;
  int j;
  char * seq_str;

  strip_out_integer_argument(&argc,argv,"kmer",&kmer_size);

  strip_out_boolean_def_argument(&argc,argv,"use_hash",&use_28bithash);
  strip_out_boolean_def_argument(&argc,argv,"use_glib",&use_glib);

  strip_out_integer_argument(&argc,argv,"min_count",&min_count);

  if( argc != 2 ) {
    show_help(stdout);
    exit(12);
  }

  kca = new_KmerCountAllocator();

  if( use_28bithash == TRUE )  {
    kii = new_KmerHashIndex(kmer_size);
  } else if( use_glib == TRUE ) {
    kii = new_interface_KmerGlibIndex(kmer_size);
  } else {
    kii = new_KmerDirectIndex(kmer_size);
  }

  sps = new_SinglePosSpace(1,2000);


  ifp = openfile(argv[1],"r");
  if( ifp == NULL ) {
    fatal("Could not open file %s",argv[1]);
  }

  while( (seq = read_fasta_Sequence(ifp)) != NULL ) {
    fprintf(stderr,"Handling sequence %s\n",seq->name);

    for(rev =0;rev < 2;rev++) {
      if( rev == 1 ) {
	tseq = seq;
	seq = reverse_complement_Sequence(seq);
	free_Sequence(tseq);
      }


      for(i=0;i<seq->len - kii->kmer_size;i++) {
	for(j=0;j< kii->kmer_size+1;j++) {
	  c = toupper(seq->seq[(int)(i+j)]);
	  if( c != 'A' && c != 'T' && c != 'G' && c != 'C' ) {
	    break;
	  }
	}
	if( j < kii->kmer_size+1 ) {
	  continue;
	}
	seq_str = seq->seq;
	
	/*  fprintf(stderr,"Handling sequence %s.. %d\n",seq->name,i);*/
	
	next_number = forward_dna_number_from_string(seq_str+i, kii->kmer_size);
	
	count = (KmerCount*) (*kii->retrieve_by_kmer)(kii->handle,next_number);
	if( count != NULL ) {
	  count->count++;
	} else {
	  count = new_KmerCount_KmerCountAllocator(kca);
	  count->count = 1;
	  (*kii->insert_by_kmer)(kii->handle,next_number,count);
	}
	
      }
    }
    free_Sequence(seq);

  }

  if( show_stdout ) {
    next_number = -1;

    while( (next_number = (*kii->next_filled_kmer)(kii->handle,next_number)) != -1 ) {
      /*      fprintf(stderr,"Traversed with %ld\n",next_number); */
      count = (KmerCount*) (*kii->retrieve_by_kmer)(kii->handle,next_number);
      if( count->count >= min_count ) {
	reverse_map_dna_number(next_number,kii->kmer_size,buffer);
	buffer[kii->kmer_size] = '\0';
	fprintf(stdout,"%s %d\n",buffer,count->count);
      }
    }
  }

  return 0;

}
