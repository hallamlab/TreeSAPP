#include "compugen.h"



/*
 * This function is passed in as a pointer to a function to
 * strcat_ExprTree_Scoped found in type.c
 *
 * It terminates the parse in the node which statisifies the
 * k position (the range looping system)
 *
 */
boolean place_k_loop(ExprTree * node,char * buffer,void * data)
{
  OneModelTrans * trans;
  trans = (OneModelTrans * )data;

  if( trans->wseq_node_to_replace == node ) {
    if( trans->map_func == NULL ) {
      strcat(buffer,"k");
    } else {
      strcat(buffer,trans->map_func);
      strcat(buffer,"(k)");
    }
    return TRUE;
  }

  return FALSE;
}

  
boolean write_cugen_model(OneModel *om,DYNFILE * dfp)
{
  FuncInfo * fi;
  char   tprf[CUGEN_MAX_PRF_NUM_LEN+1], w_base[CUGEN_MAX_PRF_NUM_LEN+1];
  char * tseq, * wseq;
  int    i;

  /*** prepare function information ***/

  
  fi = FuncInfo_named_from_varstr(FI_CALLABLE,"write_OneModel_model_%s",om->name);
  add_line_to_Ftext(fi->ft,"This function writes the onemodel datastructure",om->name);
  add_break_to_Ftext(fi->ft);
  add_break_to_Ftext(fi->ft);
  start_function_FuncInfo(fi,dfp,"boolean write_OneModel_model_%s(FILE * ofp)",om->name);

  add_break(dfp);
  
  expr(dfp,"fprintf(ofp,\"DY_MODEL : %s;\\n\");",om->name);
  expr(dfp,"fprintf(ofp,\"data_types : %s;\\n\");",om->datatype_str);

  if ( om->statenames[0] == NULL ) {
    warn ( "There are no states in the model" );
    return FALSE;
  }
  expr(dfp,"fprintf(ofp,\"state_names = { %s\");",om->statenames[0]);
  for ( i = 1; (i < CUGEN_MAX_STATE_NAMES) && (om->statenames[i] != NULL); i++ )
    expr(dfp,"fprintf(ofp,\", %s\");",om->statenames[i]);
  expr(dfp,"fprintf(ofp,\" };\\n\");");

  if ( om->semistatenames[0] != NULL ) {
    expr(dfp,"fprintf(ofp,\"semistate_names = { %s\");",om->semistatenames[0]);
    for ( i = 1; (i < CUGEN_MAX_STATE_NAMES) && (om->semistatenames[i] != NULL); i++ )
      expr(dfp,"fprintf(ofp,\", %s\");",om->semistatenames[i]);
    expr(dfp,"fprintf(ofp,\" };\\n\");");
  }

  if ( om->midstatenames[0] != NULL ) {
    expr(dfp,"fprintf(ofp,\"midstate_names = { %s\");",om->midstatenames[0]);
    for ( i = 1; (i < CUGEN_MAX_STATE_NAMES) && (om->midstatenames[i] != NULL); i++ )
      expr(dfp,"fprintf(ofp,\", %s\");",om->midstatenames[i]);
    expr(dfp,"fprintf(ofp,\" };\\n\");");
  }

  expr(dfp,"fprintf(ofp,\"nat_query = genprof;\\n\");");

  expr(dfp,"fprintf(ofp,\"line_size = %d;\\n\");", om->prof_line_len);

  if ( om->seqfuncs[0] != NULL ) {
    expr(dfp,"fprintf(ofp,\"sequence : %s\");",om->seqfuncs[0]);
    for ( i = 1; (i < CUGEN_MAX_FUNCS) && (om->seqfuncs[i] != NULL); i++ )
      expr(dfp,"fprintf(ofp,\", %s\");",om->seqfuncs[i]);
    expr(dfp,"fprintf(ofp,\";\\n\");");
  }

  expr(dfp,"fprintf(ofp,\"dy_transition from       to         dx dy t_prf dprf t_seq      dseq w_base w_seq      dwx dwy xlabel     ylabel\\n\");");

  for ( i = 0; (i < CUGEN_MAX_ONEMODEL_TRANS) && (om->trans[i] != NULL); i++ ) {

    if ( (om->trans[i]->tprf != NULL) || (om->trans[i]->tpos_tprf == om->prof_line_len - 1) )
      sprintf ( tprf, "%d", om->trans[i]->tpos_tprf );
    else
      strcpy ( tprf, "-" );

    if ( om->trans[i]->tseq == NULL )
      tseq = "-";
    else
      tseq = om->trans[i]->tseq;

    if ( om->trans[i]->wprf == NULL )
      strcpy ( w_base, "-" );
    else
      sprintf ( w_base, "%d", om->trans[i]->tpos_wprf );

     if ( om->trans[i]->wseq == NULL )
      wseq = "-";
    else
      wseq = om->trans[i]->wseq;
   

    expr(dfp,"fprintf(ofp,\"%-13s %-10s %-10s %2d %2d %-5s %4d %-10s %4d %-6s %-10s %3d %3d %-10s %-10s\\n\");",
	 om->trans[i]->name, om->trans[i]->from_state, om->trans[i]->to_state, om->trans[i]->dx,
	 om->trans[i]->dy, tprf, om->trans[i]->dprf, tseq, om->trans[i]->dseq, w_base, wseq,
	 om->trans[i]->dwx, om->trans[i]->dwy, om->trans[i]->xlabel, om->trans[i]->ylabel);
  }

  expr(dfp,"fprintf(ofp,\"end_transition\\n\");"); 
 
  expr(dfp,"fprintf(ofp,\"endmodel\\n\");");  

  close_function(dfp);
  add_break(dfp);

  return TRUE;
}

  
boolean write_cugen_profile(OneModel * om,MethodTypeSet * mts,DYNFILE * dfp)
{
  FuncInfo * fi;

  char buffer[MAXLINE]; /* place to put strcat function */
  int i;
  int range;

  /*** prepare function information ***/

  
  fi = FuncInfo_named_from_varstr(FI_CALLABLE,"write_OneModel_profile_%s",om->name);
  add_line_to_Ftext(fi->ft,"This function writes the onemodel datastructure",om->name);
  add_break_to_Ftext(fi->ft);
  add_break_to_Ftext(fi->ft);


  /*  add_args_GenericMatrix_FuncInfo(fi,om->arg_str); */

  
  start_function_FuncInfo(fi,dfp,"boolean write_OneModel_profile_%s(%s * mat,FILE * ofp)",om->name,om->name);
  /*  expr(dfp,"%s * mat",om->name); */
  expr(dfp,"int genericprofile[%d];", om->prof_line_len);
  expr(dfp,"int i,j,k;");

  /** started function, now see what basematrix says about memory **/
  add_break(dfp);

  /*  expr(dfp,"mat = allocate_%s_only(%s);",om->name,om->chain_str); */

  expr(dfp,"fprintf(ofp,\"#DYNAMITE\\n\");");
  expr(dfp,"fprintf(ofp,\"%d\\n\");", om->prof_line_len);

  expr(dfp,"for(i=0;i<mat->leni;i++)");
 
  startbrace_tag(dfp,"Over all positions in query");
 
  /* Filling in the 4*states_num first places of the profile */
  expr(dfp,"genericprofile [%d] = %d;", OM_TMPSTART, 0);
  expr(dfp,"genericprofile [%d] = %d;", OM_TMPSTART+1, 0);
  expr(dfp,"genericprofile [%d] = %d;", OM_TMPSTART+2, CUGEN_MININF);
  expr(dfp,"genericprofile [%d] = %d;", OM_TMPSTART+3, CUGEN_MININF);
  
  expr(dfp,"genericprofile [%d] = %d;", OM_TMPEND*4, CUGEN_MININF);
  expr(dfp,"genericprofile [%d] = %d;", OM_TMPEND*4+1, CUGEN_MININF);
  expr(dfp,"genericprofile [%d] = %d;", OM_TMPEND*4+2, 0);
  expr(dfp,"genericprofile [%d] = %d;", OM_TMPEND*4+3, 0);
  
  expr(dfp,"for(j=%d;j<%d;j++)", 4 * OM_TMP_STATES_NUM, 4 * (om->states_num + om->semistates_num));
  hang_expr(dfp,"genericprofile [j] = %d;", CUGEN_MININF);

  for(i=0;om->trans[i] != NULL;i++) {
    auto OneModelTrans * trans = om->trans[i];

    add_block_comment(dfp,"Doing state %s to %s",trans->from_state,trans->to_state);
    /* we already have the component parts! Hurray! */


    if( trans->tprf != NULL ) {
      buffer[0] ='\0';
      if(strcat_ExprTree_Scoped(trans->tprf,buffer,om->sc,mts,NULL,NULL,NULL) == 1 ) {
	warn("Parse error in building compugen generic profile structure");
	return FALSE;
      }
      
      /*	if( range != 1 ) {
		warn("You have a tprf line with a range != 1. Ooops! [%d]",range);
		}
      */
      
      
      
      /* simply allocate it to this tpos */
      
      expr(dfp,"genericprofile [%d] = %s;",trans->tpos_tprf,buffer);
    }
    buffer[0]='\0';
    range = -1;

    if( trans->wprf != NULL ) {


      if( strcat_ExprTree_Scoped(trans->wprf,buffer,om->sc,mts,NULL,place_k_loop,(void *)trans) == 1 ) {
	warn("Parse error in building compugen generic profile structure");
	return FALSE;
      }
      
      /* loop over the range. It has already been entered as k */
      
      expr(dfp,"for(k=0;k<%d;k++)",trans->range);
      hang_expr(dfp,"genericprofile [%d+k] = %s;",trans->tpos_wprf,buffer);
    }
  }

  expr(dfp,"genericprofile [%d] = %d;", om->prof_line_len - 1, CUGEN_MININF);

  expr(dfp,"fprintf(ofp,\"a\");");
  expr(dfp,"for(j=0;j<%d;j++)", om->prof_line_len);
  startbrace(dfp);
  expr(dfp,"fprintf(ofp,\" %%.2f\", (float)(genericprofile[j]));");
  expr(dfp,"if( j%%20 == 0 && j != 0)");
  hang_expr(dfp,"fprintf(ofp,\"\\n\");");
  closebrace(dfp); /* loop over the print out */
  expr(dfp,"fprintf(ofp,\"\\n\");");
  closebrace(dfp); /* the loop over i */
  close_function(dfp);
  add_break(dfp);

  return TRUE;
}

					    


boolean write_cugen_funcs(OneModel * om,MethodTypeSet * mts,DYNFILE * dfp)
{
  boolean ret = TRUE;
  Scope * sc;

  sc = std_Dynamite_Scope();


  if( write_cugen_model(om,dfp) == FALSE ) {
    ret = FALSE;
  }

  if( write_cugen_profile(om,mts,dfp) == FALSE ) {
    ret = FALSE;
  }

  return ret;
}
