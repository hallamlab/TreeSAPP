#include "pairbase.h"
#include "pairbaseseq.h"
#include "complexevalset.h"


char * seq1 = "ATTGGGGGTGACCTGTTAGT";
char * seq2 = "AGGGTGGGGGACTTGATAGT";


int main(int argc,char ** argv)
{

  PairBaseSeq * seq;
  ComplexSequence * cs;
  ComplexSequenceEval * splice5;
  ComplexSequenceEval * splice3;
  SeqAlign * sa;

  pairbase_type pb;


  pb = MAKE_PAIRBASE(BASE_A,BASE_T);

  if( anchor_base_from_pairbase(pb) != BASE_A ) {
    printf("not ok 1\n");
  } else {
    printf("ok 1\n");
  }

  if( informant_base_from_pairbase(pb) != BASE_T ) { 
    printf("not ok 2\n");
  } else {
    printf("ok 2\n");
  }


  seq = new_PairBaseSeq_strings(seq1,seq2);

  if( seq == NULL ) {
    printf("not ok 3\n");
  } else {
    printf("ok 3\n");
  }

  sa = SeqAlign_from_PairBaseSeq(seq);

  if( sa == NULL ) {
    printf("not ok 4\n");
  } else {
    printf("ok 4\n");
  }


  if( strcmp(sa->seq[0]->seq,seq1) != 0 ) {
    printf("not ok 5\n");
    printf("Got %s\n",sa->seq[0]->seq);
  } else {
    printf("ok 5\n");
  }

  if( strcmp(sa->seq[1]->seq,seq2) != 0 ) {
    printf("not ok 6\n");
    printf("Got %s\n",sa->seq[1]->seq);
  } else {
    printf("ok 6\n");
  }

  splice5 = stupid_5SS();
  splice3 = stupid_3SS();

  cs = ComplexSequence_from_PairBaseSeq(seq,splice5,splice3);

  printf("ok 7\n");

  show_ComplexSequence(cs,stdout);
}




