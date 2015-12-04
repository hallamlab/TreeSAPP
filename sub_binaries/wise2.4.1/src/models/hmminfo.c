#define WISE2_CROSS_HMMER2
#include "wise2xhmmer2.h"

#include "threestatemodel.h"



main(int argc,char ** argv)
{
  int i;
  ThreeStateModel * tsm;
  int show_all = 1;

  tsm = HMMer2_read_ThreeStateModel(argv[1]);

  for(i=0;i<tsm->len;i++) {
    fprintf(stdout,"%d\t%.4f\t%c\n",i,information_from_ThreeStateUnit(tsm->unit[i],tsm->rm),tsm->unit[i]->display_char);
  }


  if( show_all == 1 ) {
    for(i=0;i<26;i++) {
     fprintf(stdout,"%d\t%c\t%f\n",i,(char)('A'+i),tsm->rm->aminoacid[i]);
    }
  }


}

