#ifdef _cplusplus
extern "C" {
#endif
#include "telegraph.h"


/* Function:  write_Telegraph_grammar(*gm,*ofp)
 *
 * Descrip:    Write out a telegraph pseudo-grammar
 *
 *
 * Arg:         *gm [UNKN ] Undocumented argument [GenericMatrix]
 * Arg:        *ofp [UNKN ] Undocumented argument [FILE]
 *
 */
# line 17 "telegraph.dy"
void write_Telegraph_grammar(GenericMatrix *gm,FILE *ofp)
{
  int i;
  int j;
  int k;

  fprintf(ofp,"/* pseudo-grammar generated from dynamite */\n");
  fprintf(ofp,"\n");
  fprintf(ofp,"grammar %s {\n\n",gm->name);

  for(i=0;i<gm->len;i++) {
    fprintf(ofp,"  state %s;\n",gm->state[i]->name);
  }

  fprintf(ofp,"\n\n");

  
  for(i=0;i<gm->len;i++) {
    for(j=0;j<gm->state[i]->len;j++) {
      fprintf(ofp,"  %12s -> ",gm->state[i]->source[j]->state_source);
      for(k=0;k < gm->state[i]->source[j]->offi;k++) {
	fprintf(ofp,"Q");
      } 
      fprintf(ofp," ");
      for(k=0;k < gm->state[i]->source[j]->offj;k++) {
	fprintf(ofp,"T");
      } 
      fprintf(ofp," %12s { %s + %s};\n",gm->state[i]->name,gm->state[i]->source[j]->source_expr,gm->state[i]->source_expr == NULL ? "0" : gm->state[i]->source_expr);
    }
  }

  fprintf(ofp,"\n\n}\n\n");


}


# line 51 "telegraph.c"

#ifdef _cplusplus
}
#endif
