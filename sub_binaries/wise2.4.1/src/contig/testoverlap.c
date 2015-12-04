#include "overlap.h"





int main ( int argc, char ** argv)
{
  ContigAssembly * ca;
  FILE * ifp;

  ifp = fopen(argv[1],"r");


  ca = ContigAssembly_read(ifp);

  ContigAssembly_validate(ca);

  OverlapCollection_sort_by_score(ca->overlap);

  ContigAssembly_dump(ca,stdout);

}
