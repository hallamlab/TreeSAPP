#include "orthoset.h"


int main(int argc,char ** argv)
{
  OrthoGraph * g;
  OrthoGraph * m;

  OrthoGraphLinkPara * p;

  g = OrthoGraph_alloc_std();

  read_pairwise_OrthoSeqPos(g,stdin);

  dump_simple_OrthoGraph(g,stdout);


  m = merged_OrthoGraph_from_pairwise(g);

  fprintf(stdout,"//\n");

  dump_simple_OrthoGraph(m,stdout);

  fprintf(stdout,"//\n");

  p = OrthoGraphLinkPara_alloc();

  add_OrthoLinks_to_OrthoGraph(m,p);

  dump_simple_OrthoGraph(m,stdout);
  
}

