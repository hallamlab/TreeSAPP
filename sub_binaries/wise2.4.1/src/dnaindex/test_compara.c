#include "comparapath.h"
#include "largeseqreader.h"
#include "kmer_direct.h"
#include "kmer_hash.h"
#include <sys/time.h>
#include "hitlist.h"

void print_timing(
struct timeval *tv1,
struct timeval *tv2,
const char *str
){
  double diff = 1000000.0*(tv2->tv_sec - tv1->tv_sec) + (tv2->tv_usec - tv1->tv_usec);
  diff /= 1000000.0;

  fprintf(stderr, "TIMING: %8.3f %s\n", diff, str);
}


int main(int argc,char ** argv)
{
  struct timeval tv1, tv2, tv3, tv4, tv5, tv6, tv7, tv8, tv9, tv10;

  ComparaIndex * ci;
  Sequence * seq1;
  Sequence * seq2;
  int i;

  HitList * hl;
  HitListOutputImpl * hloi;
  LinearHSPmanager * lm;

  ComparaLinkStart * q;
  ComparaLinkStart * t;

  KmerIndexInterface * kii;

  HSPset * out;

  long int spline;

  kii = new_KmerDirectIndex(12);
  
  /*  kii = new_KmerHashIndex(22);*/
  

  hloi = new_HitListOutputImpl_from_argv(&argc,argv);

  ci = new_ComparaIndex(kii);

  /*  fprintf(stderr,"Ratio to genome is %.2f\n",(double)3000000000 / (double) ci->index_len);*/

  gettimeofday(&tv1, NULL);
  seq1 = read_large_dna_Sequence_file(argv[1],10000000);

  gettimeofday(&tv2, NULL);
  seq2 = read_large_dna_Sequence_file(argv[2],10000000);

  printf("add_Sequence_ComparaIndex(ci,seq1,0,1000,0)\n");
  gettimeofday(&tv3, NULL);
  q = add_Sequence_ComparaIndex(ci,seq1,0,5000000,0);

  printf("add_Sequence_ComparaIndex(ci,seq2,0,5000000,0)\n");
  gettimeofday(&tv4, NULL);
  t = add_Sequence_ComparaIndex(ci,seq2,1,5000000,0);

  printf("insert_revcom_Splines(ci)\n");

  gettimeofday(&tv5, NULL);
  spline = insert_revcom_Splines(ci,q);

  fprintf(stdout,"Inserted %ld splines\n",spline);

  printf("HSPset_from_ComparaIndex(ci,q)\n");
  gettimeofday(&tv6, NULL);

  out = HSPset_from_ComparaIndex(ci,q);

  lm = LinearHSPmanager_alloc_std();
  add_LinearHSPmanager(lm,out);

  printf("HitList_from_LinearHSPmanager(lm)\n");
  gettimeofday(&tv7, NULL);
  /*  hl = HitList_from_LinearHSPmanager(lm);*/

  printf("show_stats_ComparaIndex(ci,q,stdout)\n");
  gettimeofday(&tv8, NULL);
  show_stats_ComparaIndex(ci,q,stdout);

  printf("show_HSPset(out,stdout)\n");
  gettimeofday(&tv9, NULL);
  show_HSPset(out,stdout); 

  gettimeofday(&tv10, NULL);

  fflush(stdout);

  /*  free_KmerHashIndex(kii->handle);*/
  /*  show_HitList_HitListOutputImpl(hloi,hl,stdout);*/

  print_timing(&tv1, &tv2, "read seq1");
  print_timing(&tv2, &tv3, "read seq2");
  print_timing(&tv3, &tv4, "add seq1");
  print_timing(&tv4, &tv5, "add seq2");
  print_timing(&tv5, &tv6, "insert_revcom_Splines");
  print_timing(&tv6, &tv7, "HSPset_from_ComparaIndex");
  print_timing(&tv7, &tv8, "HitList_from_LinearHSPmanager");
  print_timing(&tv8, &tv9, "show_stats_ComparaIndex");
  print_timing(&tv9, &tv10, "show_HSPset");
  print_timing(&tv1, &tv10, "TOTAL");
  fflush(stderr);

}
