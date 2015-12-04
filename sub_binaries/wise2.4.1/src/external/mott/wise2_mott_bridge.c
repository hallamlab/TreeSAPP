
#include "wise2_mott_bridge.h"
#include "mott_api.h"

#include "compmat.h"

static boolean has_init = FALSE;


static void init_mott_evalue_standard_sw(void) 
{
  CompMat * mat;
  int ** matrix;
  int i,j;

  matrix = calloc(100,sizeof(int*));
  for(i='A';i<'Z'+1;i++) {
    matrix[i] = calloc(100,sizeof(int));
  }
  
  mat = read_Blast_file_CompMat("blosum62.bla");

  for(i='A';i<='Z';i++) {
    for(j='A';j<='Z';j++) {
      matrix[i][j] = mat->comp[i-'A'][j-'A'];
    }
  }

  if( (i=InitPvaluesMott(matrix,12,2,1.0e-5)) != 1 ) {
    fatal("Unable to initiate Mott scores. [%d] Have to bail!",i);
  }

  free_CompMat(mat);

  has_init = TRUE;
}

static double calc_evalue_Mott(void * data,Sequence * a,Sequence * b,int raw_score,int database_size)
{
  double prob;
  int ok;

  if( has_init == FALSE ) {
    init_mott_evalue_standard_sw();
  }

  prob = SW_PValueMott(raw_score,a->seq,b->seq,a->len,b->len,&ok);

  return prob * database_size;

}

static double calc_bits_Mott(void * data,int lengtha,int lengthb,int raw_score)
{
  return (raw_score/2.0);
}

static char * attribution_Mott(void * data)
{
  return "Mott R (2000) Accurate Formula for P-values of gapped local sequence\n and profile alignments  J. Mol Biol. 300:649-659";
}

static void free_data_noop(void * data)
{
  return;
}

SearchStatInterface * new_Mott_SearchStatInterface(void)
{
  SearchStatInterface * out;

  out = SearchStatInterface_alloc();

  out->calc_evalue = calc_evalue_Mott;
  out->calc_bits   = calc_bits_Mott;
  out->attribution = attribution_Mott;
  out->free_data   = free_data_noop;

  return out;
}


