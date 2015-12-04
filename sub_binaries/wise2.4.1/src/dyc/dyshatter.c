#ifdef _cplusplus
extern "C" {
#endif
#include "dyshatter.h"


/* Function:  write_shatter_functions(dfp,gm,dpi)
 *
 * Descrip:    Writes shatter functions
 *
 *
 * Arg:        dfp [UNKN ] Undocumented argument [DYNFILE *]
 * Arg:         gm [UNKN ] Undocumented argument [GenericMatrix *]
 * Arg:        dpi [UNKN ] Undocumented argument [DPImplementation *]
 *
 */
# line 19 "dyshatter.dy"
void write_shatter_functions(DYNFILE * dfp,GenericMatrix * gm,DPImplementation * dpi)
{

  macro(dfp,"#define %s_SHATTER_SPECIAL(matrix,i,j,STATE) matrix->shatter->special[STATE][j]",gm->name);
  macro(dfp,"#define %s_SHATTER_MATRIX(matrix,i,j,STATE)  fetch_cell_value_ShatterMatrix(mat->shatter,i,j,STATE)",gm->name);

  add_break(dfp);
	
  write_shatter_read_func(dfp,gm);
  
  write_shatter_access_funcs(dfp,gm);

  matrix_calculate_shatter_func(dfp,gm);

}



# line 37 "dyshatter.dy"
void write_shatter_access_funcs(DYNFILE * dfp,GenericMatrix * gm)
{
  FuncInfo * fi;

  fi = FuncInfo_named_from_varstr(FI_INTERNAL,"%s_shatter_access_main",gm->name);
  start_function_FuncInfo(fi,dfp,"int %s_shatter_access_main(%s * mat,int i,int j,int state)",gm->name,gm->name);
  expr(dfp,"return %s_SHATTER_MATRIX(mat,i,j,state);",gm->name);
  close_function(dfp);
  add_break(dfp);

  fi = FuncInfo_named_from_varstr(FI_INTERNAL,"%s_shatter_access_special",gm->name);
  start_function_FuncInfo(fi,dfp,"int %s_shatter_access_special(%s * mat,int i,int j,int state)",gm->name,gm->name);
  expr(dfp,"return %s_SHATTER_SPECIAL(mat,i,j,state);",gm->name);
  close_function(dfp);
  add_break(dfp);
}
  

# line 55 "dyshatter.dy"
void write_shatter_read_func(DYNFILE * dfp,GenericMatrix * gm)
{

  FuncInfo * fi;

  fi = FuncInfo_named_from_varstr(FI_CALLABLE,"PackAln_read_Shatter_%s",gm->name);
  add_line_to_Ftext(fi->ft,"Reads off PackAln from shatter matrix structure",gm->name);

  start_function_FuncInfo(fi,dfp,"PackAln * PackAln_read_Shatter_%s(%s * mat)",gm->name,gm->name);
  expr(dfp,"%s_access_func_holder holder",gm->name);
  add_break(dfp);
  expr(dfp,"holder.access_main    = %s_shatter_access_main;",gm->name);
  expr(dfp,"holder.access_special = %s_shatter_access_special;",gm->name);

  expr(dfp,"assert(mat)");
  expr(dfp,"assert(mat->shatter)");
  expr(dfp,"return PackAln_read_generic_%s(mat,holder);",gm->name);
 
  close_function(dfp);
  add_break(dfp);
}



/* Function:  matrix_calculate_shatter_func(dfp,gm)
 *
 * Descrip:    for shatter
 *
 *
 * Arg:        dfp [UNKN ] Undocumented argument [DYNFILE *]
 * Arg:         gm [UNKN ] Undocumented argument [GenericMatrix *]
 *
 */
# line 82 "dyshatter.dy"
void matrix_calculate_shatter_func(DYNFILE * dfp,GenericMatrix * gm)
{
  FuncInfo * fi;
  ArgInfo * ai;
  CellSignatureSet * csig;
  int i;
  int j;

  fi = FuncInfo_named_from_varstr(FI_CALLABLE,"calculate_shatter_%s",gm->name);
  add_line_to_Ftext(fi->ft,"This function calculates the %s matrix when in shatter mode",gm->name);
  
  ai =  ArgInfo_in_FuncInfo_from_varstr(fi,"mat");

  start_function_FuncInfo(fi,dfp,"boolean calculate_shatter_%s(%s * mat,DPEnvelope * dpenv)",gm->name,gm->name);
  expr(dfp,"int i;");
  expr(dfp,"int j;");
  expr(dfp,"int k;");
  expr(dfp,"int should_calc");
  expr(dfp,"int leni;");
  expr(dfp,"int lenj;");
  expr(dfp,"int tot;");
  expr(dfp,"int num;");
  expr(dfp,"int starti;");
  expr(dfp,"int startj;");
  expr(dfp,"int endi;");
  expr(dfp,"int endj;");

  add_break(dfp);

  csig = CellSignatureSet_from_GenericMatrix(gm);
  
  expr(dfp,"int * SIG_0_0;");
  for(i=0;i<csig->len;i++) {
    expr(dfp,"int * SIG_%d_%d;",csig->sig[i]->offi,csig->sig[i]->offj);
  }

  add_break(dfp);



  expr(dfp,"leni = mat->leni;");
  expr(dfp,"lenj = mat->lenj;");
  add_break(dfp);

  expr(dfp,"mat->shatter = new_ShatterMatrix(dpenv,%d,lenj,%d);",gm->len,gm->spec_len);


  expr(dfp,"prepare_DPEnvelope(dpenv)");

  expr(dfp,"starti = dpenv->starti");
  expr(dfp,"if( starti < 0 )");
  hang_expr(dfp,"starti = 0");

  expr(dfp,"startj = dpenv->startj");
  expr(dfp,"if( startj < 0 )");
  hang_expr(dfp,"startj = 0");

  expr(dfp,"endi = dpenv->endi");
  expr(dfp,"if( endi > mat->leni )");
  hang_expr(dfp,"endi = mat->leni");
  
  expr(dfp,"endj = dpenv->endj");
  expr(dfp,"if( endj > mat->lenj )");
  hang_expr(dfp,"endj = mat->lenj");
  

  expr(dfp,"tot = (endi-starti) * (endj-startj);");
  expr(dfp,"num = 0;");
  

  /*** see if there any specials to specials to do ***/
  
  add_break(dfp);
  expr(dfp,"start_reporting(\"%s Matrix calculation: \");",gm->name);


  expr(dfp,"for(j=startj;j<endj;j++)");
  startbrace(dfp);
  expr(dfp,"auto int score");
  expr(dfp,"auto int temp");
  expr(dfp,"for(i=starti;i<endi;i++)");
  startbrace(dfp);

  add_block_comment(dfp,"Check if is in envelope - code identical to is_in_DPEnvelope, but aggressively inlined here for speed");

  expr(dfp,"should_calc = 0;");
  expr(dfp,"for(k=0;k<dpenv->len;k++)");
  startbrace(dfp);
  expr(dfp,"auto DPUnit * u;");
  expr(dfp,"u = dpenv->dpu[k];");
  expr(dfp,"switch(u->type)");
  startbrace(dfp);
  expr(dfp,"case DPENV_RECT :");
  startcase(dfp);
  expr(dfp,"if( i >= u->starti && j >= u->startj && i < (u->starti+u->height) && j < (u->startj+u->length))");
  hang_expr(dfp,"should_calc = 1");
  expr(dfp,"break;");
  closecase(dfp);
  expr(dfp,"case DPENV_DIAG :");
  startcase(dfp);
  expr(dfp,"if(  abs( (i-j) - (u->starti-u->startj)) <= u->height && i+j >= u->starti+u->startj && i+j+u->length >= u->starti+u->startj)");
  hang_expr(dfp,"should_calc = 1");
  expr(dfp,"break;");
  closecase(dfp);
  closebrace(dfp);
  expr(dfp,"if( should_calc == 1 )");
  hang_expr(dfp,"break;");
  closebrace(dfp);

  expr(dfp,"if( should_calc == 0)");
  hang_expr(dfp,"continue;");
  

  add_break(dfp);
  

  expr(dfp,"SIG_0_0 = fetch_cell_from_ShatterMatrix(mat->shatter,i,j);");
  for(i=0;i<csig->len;i++) {
    expr(dfp,"SIG_%d_%d = fetch_cell_from_ShatterMatrix(mat->shatter,i-%d,j-%d);",csig->sig[i]->offi,csig->sig[i]->offj,csig->sig[i]->offi,csig->sig[i]->offj);
  }

  add_break(dfp);

  write_score_block_shatter(dfp,gm,"SHATTER_","mat","SHATTER_SPECIAL",TRUE);


  closebrace(dfp);
  /**** if there are any specials do them here ****/

  write_special_block(dfp,gm,"SHATTER_","SHATTER_SPECIAL",NULL);
  
  closebrace(dfp);
  
  /*** stop reporting ***/
  
  expr(dfp,"stop_reporting()");
    
  expr(dfp,"return TRUE");
  
  close_function(dfp);
  
  add_break(dfp);


}


# line 229 "dyshatter.dy"
void write_score_block_shatter(DYNFILE * dfp,GenericMatrix * gm,char * matrixtag,char * pointertag,char * specialtag,boolean use_special)
{
  register int i;
  register int j;
  register int k;
  


  for(i=0;i<gm->len;i++) {
    auto CellState * state;
    state = gm->state[i];
    
    add_break(dfp);

    if( state->footprint_start > 1 || state->footprint_end < 0)   {
      add_block_comment(dfp,"State %s has a footprint of %d - %d",state->name,state->footprint_start,state->footprint_end);
      expr(dfp,"if( SEQENDWITHIN(%d) != TRUE || SEQSTARTWITHIN(%d) != TRUE )",state->footprint_end,state->footprint_start);
      startbrace_tag(dfp,"Footprint exists");
    }



    add_block_comment(dfp,"For state %s",state->name);
    add_block_comment(dfp,"setting first movement to score",state->name);
    


    if( state->source[0]->position != SOURCE_POS_ALL) {
      add_block_comment(dfp,"Has restricted position");
      expr(dfp,"if( %s )",source_allowed_statement(state->source[0]->position,state->source[0]->offi,state->source[0]->offj));
      startbrace(dfp);
    }

    /*********************************************************************/
    /* this line looks like                                              */
    /*   score = ProteinMatrix_EXPL_MATRIX(mat,i-1,j-1,MATCH) + xxxxx    */
    /*********************************************************************/
    
    expr(dfp,"score = SIG_%d_%d[%s] + %s",state->source[0]->offi,state->source[0]->offj,
	state->source[0]->state_source,state->source[0]->calc_expr);


    if( state->source[0]->isspecial == TRUE ) {
      fatal("Cannot have a special to matrix transition as the first transition");
    }

    if( state->source[0]->position != SOURCE_POS_ALL) {
      closebrace(dfp);
    }

    
    if( dfp->code_debug_level > 5) {
      expr(dfp,"if( score > IMPOSSIBLY_HIGH_SCORE )");
      hang_expr(dfp,"log_full_error(WARNING,5,\"[%%4d][%%4d] State %s source %s Impossibly high score [%%d]\",i,j,score);",state->name,state->source[0]->state_source);
    }

    if( dfp->code_debug_level > 100 ) {
      expr(dfp,"fprintf(stderr,\"MATRIX: [%%4d][%%4d] State %s source %s got score %%d\\n\",i,j,score);",state->name,state->source[0]->state_source);
    }

    /**** ok this is to stop underflow, but is v.v.v. hacky ****/
    /*** removing underflow hack 
    expr(dfp,"if(score < (-10000000) )");
    hang_expr(dfp,"score = (-10000000)");
    ****/
    
    /****************************************/
    /* now we do if then on score and temp  */
    /****************************************/
    
    for(j=1;j<state->len;j++)	{

      if( use_special == FALSE && state->source[j]->isspecial == TRUE ) 
	continue; /** don't use the special! **/

      if( state->source[j]->position != SOURCE_POS_ALL) {
	add_block_comment(dfp,"Has restricted position");
	expr(dfp,"if( %s )",source_allowed_statement(state->source[j]->position,state->source[j]->offi,state->source[j]->offj));
	startbrace(dfp);
      }



      add_block_comment(dfp,"From state %s to state %s",state->source[j]->state_source,
			state->name);	
      if( state->source[j]->isspecial == TRUE )
	expr(dfp,"temp = %s_%s(%s,i-%d,j-%d,%s) + %s",gm->name,specialtag,pointertag,
	     state->source[j]->offi,state->source[j]->offj,state->source[j]->state_source,
	     state->source[j]->calc_expr);

      else	expr(dfp,"temp = SIG_%d_%d[%s] + %s",
		     state->source[j]->offi,state->source[j]->offj,
		     state->source[j]->state_source,
		     state->source[j]->calc_expr);
     


      if( dfp->code_debug_level > 5) {
	expr(dfp,"if( temp > IMPOSSIBLY_HIGH_SCORE )");
      hang_expr(dfp,"log_full_error(WARNING,5,\"[%%4d][%%4d] State %s source %s Impossibly high score [%%d]\",i,j,temp);",state->name,state->source[j]->state_source);
      }

      if( dfp->code_debug_level > 100 ) {
	expr(dfp,"fprintf(stderr,\"MATRIX: [%%4d][%%4d] State %s source %s got score %%d\\n\",i,j,temp);",state->name,state->source[0]->state_source);
      }

 
      /**** ok this is to stop underflow, but is v.v.v. hacky ****/
      /**** removing underflow hack
      expr(dfp,"if(score < (-10000000) )");
      hang_expr(dfp,"score = (-10000000)");
      ****/

			
      /** if we have a specified calcfunc - use it here **/
      if(gm->calcfunc != NULL ) {
	expr(dfp,"score = %s(score,temp);",gm->calcfunc);
      } else{

	expr(dfp,"if( temp  > score )");
	startbrace(dfp);
	expr(dfp,"score = temp;");
	/** ok for shadow matrix should put things in here */
	closebrace(dfp);

      }
      if( state->source[j]->position != SOURCE_POS_ALL) {
	closebrace(dfp);
      }

    }
    
    /************************/
    /* finished blocks      */
    /* put in global calc   */
    /************************/
    add_break(dfp);
    add_block_comment(dfp,"Ok - finished max calculation for %s",state->name);
    add_block_comment(dfp,"Add any movement independant score and put away"); 
    
    if( state->calc_expr != NULL)
      expr(dfp," score += %s",state->calc_expr); 

    /***************************/			
    /* put away score          */
    /***************************/
    
    expr(dfp," SIG_0_0[%s] = score;",state->name);


    if( use_special == FALSE ) {
      add_block_comment(dfp,"Finished calculating state %s",state->name);
      continue;
    }

    
    
    /************************/
    /* for each special     */
    /* thats has this as    */
    /* source we have to    */
    /* update               */
    /************************/
    
    for(j=0;j<gm->spec_len;j++) {
      auto CellState * specstate;
      specstate = gm->special[j];
      
      
      for(k=0;k<specstate->len;k++) {
	if( strcmp(specstate->source[k]->state_source,state->name) == 0) {
	  /********************************/
	  /* is a special source!         */
	  /********************************/
	  add_break(dfp);
	  add_block_comment(dfp,"state %s is a source for special %s",state->name,specstate->name);


	  if( specstate->source[k]->position != SOURCE_POS_ALL) {
	    add_block_comment(dfp,"Has restricted position");
	    expr(dfp,"if( %s )",source_allowed_statement(specstate->source[k]->position,specstate->source[k]->offi,specstate->source[k]->offj));
	    startbrace(dfp);
	  }


	  expr(dfp,"temp = score + (%s) + (%s) ",specstate->source[k]->calc_expr,specstate->calc_expr == NULL ? "0" : specstate->calc_expr );


	  if(gm->calcfunc != NULL ) {
	    expr(dfp,"%s_%s(%s,i,j,%s) = %s(%s_%s(%s,i,j,%s),temp);",
		 gm->name,specialtag,pointertag,specstate->name,
		 gm->calcfunc,
		 gm->name,specialtag,pointertag,specstate->name);
	    
	  } else{

	    expr(dfp,"if( temp > %s_%s(%s,i,j,%s) ) ",gm->name,specialtag,pointertag,specstate->name);
	    startbrace(dfp);
	    expr(dfp,"%s_%s(%s,i,j,%s) = temp",gm->name,specialtag,pointertag,specstate->name);
	    
	    closebrace(dfp);
	    add_break(dfp);

	  }
	  if( specstate->source[k]->position != SOURCE_POS_ALL) {
	    closebrace(dfp);
	  }
	  
	} /* end of if this special state was a source for previous guy */
	

      } /* end for each source of the special state */

      if( dfp->code_debug_level > 4) {
	expr(dfp,"if( %s_%s(mat,0,j,%s) > IMPOSSIBLY_HIGH_SCORE )",gm->name,specialtag,specstate->name);
	hang_expr(dfp,"log_full_error(WARNING,5,\"[%%4d][%%4d] Special state %s Impossibly high score [%%d] found\",i,j,%s_%s(mat,0,j,%s));",specstate->name,gm->name,specialtag,specstate->name);
      }


    }  /* end for each special state */
    

    
    if( state->footprint_start < 0 || state->footprint_end > 1 )
      {
	closebrace(dfp);
      }
    
    
    add_break(dfp);
    add_block_comment(dfp,"Finished calculating state %s",state->name);
    
  } /* end of for each state */
	
}
  
# line 469 "dyshatter.c"

#ifdef _cplusplus
}
#endif
