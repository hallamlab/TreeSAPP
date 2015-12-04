#ifdef _cplusplus
extern "C" {
#endif
#include "dyna2.h"



 /*********************************/
 /* Access and checking functions */
 /*********************************/

# line 147 "dyna2.dy"
CellSignatureSet * CellSignatureSet_from_GenericMatrix(GenericMatrix * gm)
{
  int i;
  int j;
  int k;
  CellSignatureSet * out;
  CellSignature * sig;

  out = CellSignatureSet_alloc_std();

  for(i=0;i<gm->len;i++) {
    for(j=0;j<gm->state[i]->len;j++) {
      
      for(k=0;k<out->len;k++) {
	if( gm->state[i]->source[j]->offi == out->sig[k]->offi &&
	    gm->state[i]->source[j]->offj == out->sig[k]->offj) {
	  break;
	}
      }
      if( k >= out->len ) {
	/* new */
	sig = CellSignature_alloc();
	sig->offi = gm->state[i]->source[j]->offi;
	sig->offj = gm->state[i]->source[j]->offj;
	add_CellSignatureSet(out,sig);
      }
    }
  }

  return out;
}

# line 179 "dyna2.dy"
boolean  can_do_threads(GenericMatrix * gm)
{
  if( gm->qtype != NULL ) {
    if( gm->qtype->is_thread_safe == FALSE ) 
      return FALSE;
  }

  if( gm->ttype != NULL ) {
    if( gm->ttype->is_thread_safe == FALSE ) 
      return FALSE;
  }

  return TRUE;
}
    

# line 195 "dyna2.dy"
CellState * CellState_from_str(GenericMatrix * gm,char * str)
{
  register int i;
  
  for(i=0;i<gm->len;i++)
    if( strcmp(gm->state[i]->name,str) == 0)
      return gm->state[i];
  
  
  for(i=0;i<gm->spec_len;i++)
    if( strcmp(gm->special[i]->name,str) == 0)
      return gm->special[i];
  
  
  return NULL;
}


/* Function:  prepare_matrix(gm,mts,dycw,failing_errors)
 *
 * Descrip:    main function to check GenericMatrix onced parsed
 *
 *             checks
 *               state defaults
 *               state/source cross references
 *               labels
 *               calc epxressions
 *               types and type migration
 *               calc parsing
 *
 *
 * Arg:                    gm [RW   ] GenericMatrix to be checked [GenericMatrix *]
 * Arg:                   mts [READ ] Type and Method Scope [MethodTypeSet *]
 * Arg:                  dycw [UNKN ] Undocumented argument [DycWarning *]
 * Arg:        failing_errors [READ ] Calc line parser on which errors fail [ParseError]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 228 "dyna2.dy"
boolean prepare_matrix(GenericMatrix * gm,MethodTypeSet * mts,DycWarning * dycw,ParseError failing_errors)
{
  boolean ret = TRUE;
  ParseError pe;

  push_errormsg_stack("In preparing matrix %s",gm->name);

  (void) calc_footprint(gm);

  /* this must happen first */
  if( cross_reference_state_and_source(gm) == FALSE ) {
    warn("Unable to cross reference state and source");
    ret = FALSE;
  }

  /* now this to provide score setting */

  if( check_start_end(gm) == FALSE ) {
    warn("Start/End points are faulty");
    ret = FALSE;
  }

  /* set global score */

  gm->defscore_all_states = stringalloc("NEGI");

  if( perculate_state_defaults(gm) == FALSE )  {
    warn("Unable to perculate state defaults");
    ret = FALSE;
  }


  if( assign_source_no(gm) == FALSE ) {
    warn("Weird! unable to assign source numbers!");
    ret = FALSE;
  }
 
  if( prepare_labels(gm) == FALSE ) {
    warn("Unable to prepare labels");
    ret = FALSE;
  }

  if( handle_names(gm) == FALSE ) {
    warn("Unable to handle all the resource/query/target objects");
    ret = FALSE;
  }

  if( check_source_positions(gm) == FALSE ) {
    warn("Unable to resolve positions");
    ret = FALSE;
  }

  if( calc_window(gm) == FALSE ) {
    warn("Unable to calculate window");
    ret = FALSE;
  }


  if( check_cell_refs(gm) == FALSE ) {
    warn("Unable to resolve all the cell offset refs into correct offsets");
    ret = FALSE;
  }


  if( make_StructHolder_for_GenericMatrix(gm,mts) == FALSE ) {
    warn("Unable to build structure holder for GenericMatrix");
    ret = FALSE;
  }

 

  if( ret == FALSE ) {
    warn("Failing simple cross-checks, aborting before calc-line parsing");
    pop_errormsg_stack();
    return FALSE;
  }

  pe = parse_calc_line_GenericMatrix(gm,mts,dycw);
  if( pe & failing_errors ) {
    warn("Failed to parse calc lines");
    fprintf(stderr,"\nThe following parser errors were considered fatal:\n");
    complain_ParseError_to_file(pe & failing_errors,stderr);

    if( (pe & (~failing_errors)) != 0 ) { 
      fprintf(stderr,"\nThe following parser errors are only warnings:\n");
      complain_ParseError_to_file(pe & (~failing_errors),stderr);
    }

    fprintf(stderr,"\n");

    ret = FALSE;
  }

  gm->mts = hard_link_MethodTypeSet(mts);

  /***
    Very hacky to conform to old gm code: 

    type migration moves logical types to real. This is ok for now,
    but it does mean that in the function writing code we have
    lost the logical type information.

    So... probably should have the mts in the middle of function
    code.

    ***/

  if( ret == TRUE )
    handle_type_migration(gm,mts);

  pop_errormsg_stack();

  return ret;
}

/* Function:  assign_source_no(gm)
 *
 * Descrip:    Adds a unique transition number for CellSource
 *
 *
 * Arg:        gm [UNKN ] Undocumented argument [GenericMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 346 "dyna2.dy"
boolean assign_source_no(GenericMatrix * gm)
{
  int i,j;
  int no = 0;

  for(i=0;i<gm->len;i++)
    for(j=0;j<gm->state[i]->len;j++)
      gm->state[i]->source[j]->trans_no = no++;

  for(i=0;i<gm->spec_len;i++)
    for(j=0;j<gm->special[i]->len;j++)
      gm->special[i]->source[j]->trans_no = no++;
  
  return TRUE;
}

/* Function:  add_GenericMatrix_Scope(sc,gm)
 *
 * Descrip:    Adds the current generic matrix variables to the scope with 
 *             "mat->" as scope resolver
 *
 *
 * Arg:        sc [UNKN ] Undocumented argument [Scope *]
 * Arg:        gm [UNKN ] Undocumented argument [GenericMatrix *]
 *
 */
# line 367 "dyna2.dy"
void add_GenericMatrix_Scope(Scope * sc,GenericMatrix * gm)
{
  int i;
  ScopeUnit * su;


  su = ScopeUnit_from_nat(NULL,"i","","int");
  su->no_accept = stringalloc(gm->target->name);
  add_Scope(sc,su);
  su = ScopeUnit_from_nat(NULL,"j","","int");
  su->no_accept = stringalloc(gm->query->name);
  add_Scope(sc,su);
  su = ScopeUnit_from_nat(NULL,"mat","","*Matrix*");
  add_Scope(sc,su);
  
  su = ScopeUnit_from_nat(NULL,gm->query->name,"mat->",gm->query->element_type);
  su->no_accept = stringalloc("j");

  add_Scope(sc,su);
  su = ScopeUnit_from_nat(NULL,gm->target->name,"mat->",gm->target->element_type);
  su->no_accept = stringalloc("i");

  add_Scope(sc,su);

  for(i=0;i<gm->res_len;i++) {
    su = ScopeUnit_from_nat(NULL,gm->resource[i]->name,"mat->",gm->resource[i]->element_type);
    add_Scope(sc,su);
  }
  
  for(i=0;i<gm->ev_len;i++) {
    su = ScopeUnit_from_nat(NULL,gm->ev[i]->name,"",gm->ev[i]->type);
    add_Scope(sc,su);
  }


}

/* Function:  parse_calc_line_GenericMatrix(gm,mts,dycw)
 *
 * Descrip:    Function which actually does the loop over calc lines.
 *             Importantly calls /allocd_calc_line from type.dy for
 *             the main parsing
 *
 *
 * Arg:          gm [UNKN ] Undocumented argument [GenericMatrix *]
 * Arg:         mts [UNKN ] Undocumented argument [MethodTypeSet *]
 * Arg:        dycw [UNKN ] Undocumented argument [DycWarning *]
 *
 * Return [UNKN ]  Undocumented return value [ParseError]
 *
 */
# line 410 "dyna2.dy"
ParseError parse_calc_line_GenericMatrix(GenericMatrix * gm,MethodTypeSet * mts,DycWarning * dycw)
{
  ParseError pe  = 0;
  char * temp;
  int i;
  int j;
  Scope * sc;
  ExprTree * expr;

  sc = std_Dynamite_Scope();
  
  add_GenericMatrix_Scope(sc,gm);

  gm->sc = sc;

  for(i=0;i<gm->len;i++) {
    if( gm->state[i]->calc_expr != NULL ) {
      push_errormsg_stack("In parsing calc line for state [%s] (source ind.)",gm->state[i]->name);
      temp = allocd_calc_line(gm->state[i]->calc_expr,sc,mts,dycw,&pe,&expr);
      if( temp == NULL ) {
	temp = stringalloc(gm->state[i]->calc_expr);
      }
      gm->state[i]->source_expr = gm->state[i]->calc_expr;
      gm->state[i]->calc_expr = temp;
      gm->state[i]->etr      = expr;

      pop_errormsg_stack();

    }

    for(j=0;j<gm->state[i]->len;j++) {

      push_errormsg_stack("In parsing calc line for state [%s] source [%s]",gm->state[i]->name,gm->state[i]->source[j]->state_source);

      temp = allocd_calc_line(gm->state[i]->source[j]->calc_expr,sc,mts,dycw,&pe,&expr);
      if( temp == NULL ) {
	temp = stringalloc(gm->state[i]->source[j]->calc_expr);
      }
      gm->state[i]->source[j]->source_expr = gm->state[i]->source[j]->calc_expr;
      gm->state[i]->source[j]->calc_expr = temp;
      gm->state[i]->source[j]->etr = expr;


      pop_errormsg_stack();
    }
  }

  for(i=0;i<gm->spec_len;i++) {
    if( gm->special[i]->calc_expr != NULL ) {
      push_errormsg_stack("In parsing calc line for state [%s] (source ind.)",gm->special[i]->name);
      temp = allocd_calc_line(gm->special[i]->calc_expr,sc,mts,dycw,&pe,&expr);
      if( temp == NULL ) {
	temp = stringalloc(gm->special[i]->calc_expr);
      }
      gm->special[i]->source_expr = gm->special[i]->calc_expr;
      gm->special[i]->calc_expr = temp;
      gm->special[i]->etr      = expr;

      pop_errormsg_stack();

    }

    for(j=0;j<gm->special[i]->len;j++) {

      push_errormsg_stack("In parsing calc line for state [%s] source [%s]",gm->special[i]->name,gm->special[i]->source[j]->state_source);

      temp = allocd_calc_line(gm->special[i]->source[j]->calc_expr,sc,mts,dycw,&pe,&expr);
      if( temp == NULL ) {
	temp = stringalloc(gm->special[i]->source[j]->calc_expr);
      }
      gm->special[i]->source[j]->source_expr = gm->special[i]->source[j]->calc_expr;
      gm->special[i]->source[j]->calc_expr = temp;
      gm->special[i]->source[j]->etr      = expr;


      pop_errormsg_stack();
    }
  }

  return pe;
}


/* Function:  calc_window(gm)
 *
 * Descrip:    calculates window_i and window_j
 *
 *
 * Arg:        gm [UNKN ] Undocumented argument [GenericMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 497 "dyna2.dy"
boolean calc_window(GenericMatrix * gm)
{
  register int i;
  register int j;
  int ilen = 0;
  int jlen = 0;


  for(i=0;i<gm->len;i++)
    {
      auto CellState * state;
      state = gm->state[i];
      for(j=0;j<state->len;j++)
	{
	  if( state->source[j]->offi > ilen )
	    ilen = state->source[j]->offi;
	  if( state->source[j]->offj > jlen )
	    jlen = state->source[j]->offj;
	}
    }


  gm->window_i = ilen;
  gm->window_j = jlen;

  return TRUE;
}

/* Function:  check_start_end(gm)
 *
 * Descrip:    checks we have a start + end (and only 1 each!)
 *             and sets start's defscore to 0
 *
 *
 * Arg:        gm [UNKN ] Undocumented argument [GenericMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 529 "dyna2.dy"
boolean check_start_end(GenericMatrix * gm)
{
  int i;
  int j;
  boolean ret = TRUE;

  for(i=0;i<gm->spec_len;i++) {
    if( gm->special[i]->is_end == TRUE ) {
      if( gm->special[i]->is_start == TRUE ) {
	warn("Trying to make state %s both the start and the end!",gm->special[i]->name);
	ret = FALSE;
      }
      break;
    }
  } 

  if( i == gm->spec_len ) {
    for(i=0;i<gm->spec_len;i++) {
      if( strcmp(gm->special[i]->name,"END") == 0 ) {
	warn("You have not got a !end special, but you do have a END special state. Presuming that you wanted to make that the END ;)");
	gm->special[i]->is_end = TRUE;
	break;
      }
    }
    if( i == gm->spec_len) {
      warn("You have no end special state. Impossible matrix");
      ret = FALSE;
    }
  } else {
    /*** check there are no more end states ***/

    for(j=i+1;j<gm->spec_len;j++) {
      if( gm->special[j]->is_end == TRUE ) {
	warn("Special state %s is also end (as well as %s)",gm->special[j]->name,gm->special[i]->name);
      }
    }
  } /*** end of else ***/

  /*** done end, do start ***/


  for(i=0;i<gm->spec_len;i++) {
    if( gm->special[i]->is_start == TRUE ) {
      if( gm->special[i]->is_end == TRUE ) {
	warn("Trying to make state %s both the start and the end!",gm->special[i]->name);
	ret = FALSE;
      }
      break;
    }
  } 

  if( i == gm->spec_len ) {
    for(i=0;i<gm->spec_len;i++) {
      if( strcmp(gm->special[i]->name,"START") == 0 ) {
	warn("You have not got a !start special, but you do have a START special state. Presuming that you wanted to make that the START ;)");
	gm->special[i]->is_start = TRUE;
	break;
      }
    }
    if( i == gm->spec_len) {
      warn("You have no start special state. Impossible matrix");
      ret = FALSE;
    }
  } else {
    /*** check there are no more start states ***/

    for(j=i+1;j<gm->spec_len;j++) {
      if( gm->special[j]->is_start == TRUE ) {
	warn("Special state %s is also start (as well as %s)",gm->special[j]->name,gm->special[i]->name);
      }
    }
    /*** set score for start = 0 **/

    gm->special[i]->def_score = stringalloc("0");
  } /*** end of else ***/

  return ret;
}


/* Function:  check_source_positions(gm)
 *
 * Descrip:    checks the top/bottom/left/right source
 *             positions
 *
 *
 * Arg:        gm [UNKN ] Undocumented argument [GenericMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 613 "dyna2.dy"
boolean check_source_positions(GenericMatrix * gm)
{
  int i;
  int j;
  boolean ret = TRUE;

  for(i=0;i<gm->len;i++){
    auto CellState * s;
    s = gm->state[i];
    for(j=0;j<s->len;j++) {
      if( j == 0 && s->source[0]->isspecial == TRUE ) {
	warn("For [%s] source [%s] your first transition is from a special. This is not allowed currently to make code generation easier. There must be a non-special transition to place here. could you rearrange your dy file?",s->name,s->source[0]->state_source);
	ret = FALSE;
      }

      if( s->source[j]->position < 16 ) {
	ret = FALSE;
	warn("For [%s] source [%s] you have an impossible position vector (top/left/bottom/right 1/2/4/8) %d",s->name,s->source[j]->state_source,s->source[j]->position);
      }
      else {
	if( s->source[j]->position == SOURCE_POS_ALL ) 
	  continue; /* fine */
	if( s->source[j]->isspecial == FALSE ) {
	  warn("For [%s] source [%s] you have indicated an edge to calculate on, but it is a non special source",s->name,s->source[j]->state_source);
	  ret = FALSE;
	}
	else if( s->source[j]->position != SOURCE_POS_TOPLEFT && s->source[j]->position != SOURCE_POS_TOPLEFT ) {
	  /*** allow? ***/
	  /*
	  warn("Really sorry, you have a lone top bottem left or right source which is not yet implemented. talk to Ewan for a work around");
	  ret = FALSE;
	  */
	}
      }
    }
  }

  for(i=0;i<gm->spec_len;i++){
    auto CellState * s;
    s = gm->special[i];
    for(j=0;j<s->len;j++) {
      if( s->source[j]->position < 16 ) {
	ret = FALSE;
	warn("For [%s] source [%s] you have an impossible position vector (top/left/bottom/right 1/2/4/8) %d",s->name,s->source[j]->state_source,s->source[j]->position);
      }
      else {
	if( s->source[j]->position == SOURCE_POS_ALL ) 
	  continue; /* fine */
	else if( (s->source[j]->position == SOURCE_POS_TOPLEFT || s->source[j]->position == SOURCE_POS_TOP || s->source[j]->position == SOURCE_POS_BOTTOM || s->source[j]->position == SOURCE_POS_BOTTOMRIGHT) && s->source[j]->isspecial == TRUE) {
	  warn("State [%s] source [%s] has top/bottom tags - but it is a special to special. Surely you only mean left (no top/bottom concept to specials)",s->name,s->source[j]->state_source);

	  ret = FALSE;
	}
      }
    }
  }

  return ret;
}


/* Function:  calc_footprint(gm)
 *
 * Descrip:    calculates footprint
 *
 *             now useless
 *
 *
 * Arg:        gm [UNKN ] Undocumented argument [GenericMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 680 "dyna2.dy"
boolean calc_footprint(GenericMatrix * gm)
{
  register int i;
  register int foot = 1;

  for(i=0;i<gm->len;i++)
    if( gm->state[i]->footprint_end > foot )
      foot = gm->state[i]->footprint_end;

  gm->footprint = foot;

  return TRUE;
}

/* Function:  cross_reference_state_and_source(gm)
 *
 * Descrip:    makes sure each source has a state
 *
 *
 * Arg:        gm [UNKN ] Undocumented argument [GenericMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 698 "dyna2.dy"
boolean cross_reference_state_and_source(GenericMatrix * gm)
{
  register int i;
  register int j;

  for(i=0;i<gm->len;i++) {
    gm->state[i]->state_number = i;
  }
  j = gm->len;
  
  for(i=0;i<gm->spec_len;i++) {
    gm->special[i]->state_number = i+j;
  }

  for(i=0;i<gm->len;i++) {
    auto CellState * state;
    auto CellState * temp;
    state = gm->state[i];
    for(j=0;j<state->len;j++) {
      if( (temp=CellState_from_str(gm,state->source[j]->state_source)) == NULL){
	warn("In matrix %s - State %s asks for source %s but there is no State of that name",gm->name,state->name,state->source[j]->state_source);
	return FALSE;
      }
      if( temp->is_special_i == TRUE || temp->is_special_j == TRUE)
	state->source[j]->isspecial = TRUE;
      state->source[j]->from_state_no = temp->state_number;
    }
  }

  for(i=0;i<gm->spec_len;i++) {
    auto CellState * state;
    auto CellState * temp;
    state = gm->special[i];
    for(j=0;j<state->len;j++) {

      if( (temp=CellState_from_str(gm,state->source[j]->state_source)) == NULL)
	{
	  warn("In matrix %s - State %s asks for source %s but there is no State of that name",gm->name,state->name,state->source[j]->state_source);
	  return FALSE;
	}
      if( temp->is_special_i == TRUE || temp->is_special_j == TRUE) {
	state->source[j]->isspecial = TRUE;
	state->specialtospecial = TRUE;
	gm->specialtospecial = TRUE;
      }
    }
  }

  return TRUE;
}


# line 750 "dyna2.dy"
boolean check_cell_refs(GenericMatrix * gm)
{
  boolean ret = TRUE;
  register int i;
  register int j;

  for(i=0;i<gm->len;i++) {
    for(j=0;j<gm->state[i]->len;j++) {
      if( gm->state[i]->source[j]->offi == 0 && gm->state[i]->source[j]->offj == 0 ) {
	warn("In state %s (source %s), both offi and offj are zero: dynamite cannot currently handle cell internal references",gm->state[i]->name,gm->state[i]->source[j]->state_source);
	ret = FALSE;
      }
      if( gm->state[i]->source[j]->offi < 0 || gm->state[i]->source[j]->offj < 0 ) {
	warn("In state %s, offi,offj [%d][%d] has some negative indices: offi and offj are always positive, probably you should just strip off the negative signs",gm->state[i]->name,gm->state[i]->source[j]->offi,gm->state[i]->source[j]->offj);
	ret = FALSE;
      }
    }
  }

  for(i=0;i<gm->spec_len;i++) {
    for(j=0;j<gm->special[i]->len;j++) {
      if( gm->special[i]->source[j]->isspecial == TRUE ) {
	if( gm->special[i]->source[j]->offj == 0 ) {
	  warn("In special state %s, source %s (also special) got a offj of zero, an impossible reference",gm->special[i]->name,gm->special[i]->source[j]->state_source);
	  ret = FALSE;
	}
      } else {
	if( gm->special[i]->source[j]->offj != 0 ) {
	  warn("In special state %s, source %s (a main matrix cell) got an offset of non zero (%d). This cannot be modelled at the moment",gm->special[i]->name,gm->special[i]->source[j]->state_source, gm->special[i]->source[j]->offj);
	  ret = FALSE;
	}
      }

      if( gm->special[i]->source[j]->offj  < 0 ) {
	warn("In special state %s, source %s, offj is %d (Negative!), and you can't have negative offsets",
gm->special[i]->name,gm->special[i]->source[j]->state_source);
	ret = FALSE;
      }
    }
  }

  return ret;
}


# line 795 "dyna2.dy"
boolean prepare_labels(GenericMatrix * gm)
{
  boolean all_labelled = TRUE;
  boolean ret = TRUE;
  register int i;
  register int j;

  for(i=0;i<gm->len;i++) {
    auto CellState * state;
    state = gm->state[i];
    for(j=0;j<state->len;j++) {
      if( state->source[j]->query_label == NULL ) {
	warn("For state %s, source %s, missing a query_label",state->name,state->source[j]->state_source);
	ret = FALSE;
      }
      if( state->source[j]->target_label == NULL ) {
	warn("For state %s, source %s, missing a target_label",state->name,state->source[j]->state_source);
	ret = FALSE;
      }
    }
  }

  for(i=0;i<gm->spec_len;i++) {
    auto CellState * state;
    state = gm->special[i];
    for(j=0;j<state->len;j++) {
      if( state->source[j]->query_label == NULL ) {
	warn("For state %s, source %s, missing a query_label",state->name,state->source[j]->state_source);
	ret = FALSE;
      }
      if( state->source[j]->target_label == NULL ) {
	warn("For state %s, source %s, missing a target_label",state->name,state->source[j]->state_source);
	ret = FALSE;
      }
    }
  }

  gm->canlabel = all_labelled;

  return ret;
}

# line 837 "dyna2.dy"
boolean handle_names(GenericMatrix * gm)
{
  boolean ret = TRUE;
  register int i;


  if( gm->query->name == NULL ) {
    warn("Your query [type %s]  had no name: calling it 'query'",CKS(gm->query->element_type));
    gm->query->name = stringalloc("query"); 
  }

  if( gm->target->name == NULL ) {
    warn("Your target [type %s] had no name: calling it 'target'",CKS(gm->target->element_type));
    gm->target->name = stringalloc("target");
  }

  for(i=0;i<gm->res_len;i++) {
    if( gm->resource[i]->name == NULL ) {
      warn("Resource number %d had no name... cannot process",i);
      ret = FALSE;
    }
    if( gm->resource[i]->element_type == NULL ) {
      warn("Resource number %d had no type... cannot process",i);
      ret = FALSE;
    }

  }

  return ret;
}


# line 869 "dyna2.dy"
boolean handle_type_migration(GenericMatrix * gm,MethodTypeSet * mts)
{
  register int i;
  Type * t;

  if( gm->query_len == NULL ) {
    gm->query_len = length_string_from_GenericMatrix_type(gm->query->element_type);
  }

  if( gm->target_len == NULL ) {
    gm->target_len = length_string_from_GenericMatrix_type(gm->target->element_type);
  }

  if( (t = Type_from_name(mts,gm->query->element_type)) != NULL) {
    ckfree(gm->query->element_type);
    gm->query->element_type = stringalloc(t->real);
    gm->qtype = t;
  }

  if( (t = Type_from_name(mts,gm->target->element_type)) != NULL) {
    ckfree(gm->target->element_type);
    gm->target->element_type = stringalloc(t->real);
    gm->ttype = t;
  }

  for(i=0;i<gm->res_len;i++) {
    if( (t = Type_from_name(mts,gm->resource[i]->element_type)) != NULL) {
      ckfree(gm->resource[i]->element_type);
      gm->resource[i]->element_type = stringalloc(t->real);
    }
  }

  return TRUE;
}


# line 905 "dyna2.dy"
boolean perculate_state_defaults(GenericMatrix * gm)
{
  register int i;
  register int j;
  boolean ret = TRUE;

  for(i=0;i<gm->len;i++) {
    auto CellState * state;
    state = gm->state[i];

    if( state->def_score == NULL) {
      if( gm->defscore_all_states == NULL ) {
	warn("State %s has no specific default score and you have not specified a global default score",state->name);
	ret = FALSE;
      }
      state->def_score = stringalloc( gm->defscore_all_states);
    }
    for(j=0;j<state->len;j++) {
      if( state->source[j]->offi == -1 )
	state->source[j]->offi = state->offi;
      if( state->source[j]->offj == -1 )
	state->source[j]->offj = state->offj;
      
      /*** now do labels ***/
      if( state->source[j]->query_label == NULL && state->query_label != NULL )
	state->source[j]->query_label = stringalloc(state->query_label);
      
      if( state->source[j]->target_label == NULL && state->target_label != NULL )
	state->source[j]->target_label = stringalloc(state->target_label);

     
    }
    
  }
  
  for(i=0;i<gm->spec_len;i++) {
    auto CellState * state;
    state = gm->special[i];
    
    gm->special[i]->is_special_i = TRUE;
    
    if( state->def_score == NULL) {
      if( gm->defscore_all_states == NULL ) {
	warn("State %s has no specific default score and you have not specified a global default score",state->name);
	ret = FALSE;
      }
      state->def_score = stringalloc( gm->defscore_all_states);
    }
    for(j=0;j<state->len;j++) {
      if( state->source[j]->offi == -1 ) 
	state->source[j]->offi = state->offi;
      if( state->source[j]->offj == -1 )
	state->source[j]->offj = state->offj;
      
      /*** now do labels ***/
      if( state->source[j]->query_label == NULL && state->query_label != NULL )
	state->source[j]->query_label = stringalloc(state->query_label);
      
      if( state->source[j]->target_label == NULL && state->target_label != NULL )
	state->source[j]->target_label = stringalloc(state->target_label);
      
      
    }

  }
  

  return ret;
}


# line 976 "dyna2.dy"
boolean make_StructHolder_for_GenericMatrix(GenericMatrix * gm,MethodTypeSet * mts)
{
  StructHolder * out;
  StructElement * temp;
  register int i;

  out = StructHolder_alloc_std();


  if( out == NULL ) {
    warn("Unable to build generic matrix %s's structure",gm->name);
    return FALSE;
  }

  out->name = stringalloc(gm->name);

  temp = basic_add_StructElement(out,"basematrix","BaseMatrix *");
  temp->def=stringalloc("NULL");

  if( temp == NULL ) {
    warn("Serious problem - unable to add one of the structelements in make_StructHolder_for_GenericMatrix");
    return FALSE;
  }

  temp = basic_add_StructElement(out,"shatter","ShatterMatrix *");
  temp->def=stringalloc("NULL");

  temp = basic_add_StructElement(out,"leni","int");
  temp->def=stringalloc("0");

  if( temp == NULL ) {
    warn("Serious problem - unable to add one of the structelements in make_StructHolder_for_GenericMatrix");
    return FALSE;
  }

  temp = basic_add_StructElement(out,"lenj","int");
  temp->def=stringalloc("0");

  if( temp == NULL ) {
    warn("Serious problem - unable to add one of the structelements in make_StructHolder_for_GenericMatrix");
    return FALSE;
  }



  /*** ok... this is now typed! ***/

  temp = StructElement_from_MethodTypeSet(gm->query->name,gm->query->element_type,mts);
  add_StructHolder(out,temp);
  temp = StructElement_from_MethodTypeSet(gm->target->name,gm->target->element_type,mts);
  add_StructHolder(out,temp);
  for(i=0;i<gm->res_len;i++) {
    temp = StructElement_from_MethodTypeSet(gm->resource[i]->name,gm->resource[i]->element_type,mts);
    temp->isfunc   = gm->resource[i]->isfunc;
    add_StructHolder(out,temp);
  }

  gm->sh = out;

  return TRUE;
}


# line 1039 "dyna2.dy"
StructElement * StructElement_for_GenericMatrix_type(char * name,char * element)
{
  StructElement * out;


  if( name == NULL || element == NULL ) {
    warn("In trying to make a struct element for GenericMatrix type, got some NULL fields... this is not good");
    return NULL;
  }

  out = (StructElement *) ckalloc (sizeof(StructElement));

  if( strcmp(element,"PROTEIN") == 0 || strcmp(element,"GENOMIC") == 0 || strcmp(element,"CDNA") == 0 ) {
    out->element_type = stringalloc("ComplexSequence *");
  }

  else {
    out->element_type = stringalloc(element);
  }

  out->name = stringalloc(name);
  out->islinked = TRUE;

  return out;
}

# line 1065 "dyna2.dy"
char * length_string_from_GenericMatrix_type(char * element)
{
  if( strcmp(element,"PROTEIN") == 0 )
    return stringalloc("seq->len");
  if( strcmp(element,"GENOMIC") == 0 || strcmp(element,"CDNA") == 0 || strcmp(element,"DNA") == 0)
    return stringalloc("seq->len");
  else {
    warn("Cannot automatically find length string for element type [%s] (perhaps it is not a logical type). You"
" should have a field:len in this line if it is not a logical type. Assuming [len] anyway",element);
    return stringalloc("len");
  }

}

# line 1079 "dyna2.dy"
boolean can_interpret_type(char * type)
{
  if( strcmp(type,"PROTEIN") == 0 || strcmp(type,"GENOMIC") == 0 || strcmp(type,"CDNA") == 0)
    return TRUE;
  return FALSE;
}
 
# line 1086 "dyna2.dy"
char * interpret_type(char * type)
{
  if( strcmp(type,"PROTEIN") == 0 || strcmp(type,"GENOMIC") == 0 || strcmp(type,"CDNA") == 0)
    return stringalloc("ComplexSequence *");

  return NULL;
}


# line 1095 "dyna2.dy"
CellState * start_CellState_from_GenericMatrix(GenericMatrix * gm)
{
  register int i;

  for(i=0;i<gm->spec_len;i++) {
    if(gm->special[i]->is_start == TRUE ) 
      return gm->special[i];
  }

  return NULL;
}

# line 1107 "dyna2.dy"
CellState * end_CellState_from_GenericMatrix(GenericMatrix * gm)
{
  register int i;

  for(i=0;i<gm->spec_len;i++) {
    if(gm->special[i]->is_end == TRUE ) 
      return gm->special[i];
  }

  return NULL;
}

 /***************************/
 /* noddy display functions */
 /***************************/

# line 1123 "dyna2.dy"
void show_GenericMatrix(GenericMatrix * gm,char padchar,FILE * ofp)
{
  register int i;

  fprintf(ofp,"Matrix: %s\n",gm->name);
  
  fprintf(ofp,"Query Type [%s] Name [%s]\n",gm->query->element_type,gm->query->name);
  
  fprintf(ofp,"Target Type [%s] Name [%s]\n",gm->target->element_type,gm->target->name);
  
  for(i=0;i<gm->res_len;i++)
    fprintf(ofp,"Resource Type[%s] Name [%s]\n",gm->resource[i]->element_type,gm->resource[i]->name);
  
  for(i=0;i<gm->len;i++)
    show_CellState(gm->state[i],padchar,1,ofp);
}


# line 1141 "dyna2.dy"
void show_CellState(CellState * cell,char padchar,int num,FILE * ofp)
{
  register int i;
  fprintf(ofp,"\nFor state %s\n",cell->name);
  for(i=0;i<num;i++)
    fputc(padchar,ofp);
  fprintf(ofp,"def_score: %s\n",cell->def_score == NULL ? "ANY [NULL STRING]" : cell->def_score);

  for(i=0;i<num;i++)
    fputc(padchar,ofp);
  fprintf(ofp,"calculation for each score: %s\n",CKS(cell->calc_expr));

  for(i=0;i<num;i++)
    fputc(padchar,ofp);
  fprintf(ofp,"offi: %d\n",cell->offi);

  for(i=0;i<num;i++)
    fputc(padchar,ofp);
  fprintf(ofp,"offj: %d\n",cell->offj);


  for(i=0;i<cell->len;i++)
    show_CellSource(cell->source[i],padchar,num+1,ofp);
}

# line 1166 "dyna2.dy"
void show_CellSource(CellSource * cell,char padchar,int num,FILE * ofp)
{
  register int i;
  for(i=0;i<num;i++)
    fputc(padchar,ofp);
  fprintf(ofp,"state_source: %s\n",cell->state_source == NULL ? "ANY [NULL STRING]" : cell->state_source);
  for(i=0;i<num;i++)
    fputc(padchar,ofp);
  fprintf(ofp,"offi: %d\n",cell->offi);

  for(i=0;i<num;i++)
    fputc(padchar,ofp);
  fprintf(ofp,"offj: %d\n",cell->offj);


}


/*********************/
/* Parsing functions */
/*********************/


# line 1189 "dyna2.dy"
GenericMatrix * read_GenericMatrix(FILE * ifp)
{
  char buffer[MAXLINE];
  GenericMatrix * out = NULL;

  while( get_watched_line(buffer,MAXLINE,ifp) != NULL )
    {
      if( strstartcmp(buffer,"matrix") == 0)
	{
	  out = read_GenericMatrix_line(buffer,ifp);
	  break;
	}
    }

  if( out == NULL )
    warn("Unable to read GenericMatrix");

  while( get_watched_line(buffer,MAXLINE,ifp) != NULL )
    {
      if( strstartcmp(buffer,"%{") == 0)
	break;
    }

  return out;
}


# line 1216 "dyna2.dy"
GenericMatrix * read_GenericMatrix_line(char * line,FILE * ifp)
{
  GenericMatrix * out;
  ExternVariable * ev;
  CellState * temp;
  char buffer[MAXLINE];
  char ** base;
  char ** splitstr;
  char * runner;
  StructElement * res_temp;
  boolean isspecial = FALSE;
  CollapsableLabel * cal;
  



  base=splitstr=breakstring(line,spacestr);

  if( strcmp(*splitstr,"matrix") != 0) {
    log_full_error(WARNING,0,"In parsing the line starting %s it had no matrix tag!",line);
    return NULL;
  }
  splitstr++;

  if( *splitstr == NULL ) {
    log_full_error(WARNING,0,"In parsing the first matrix line there was no matrix name: must have matrix <name>",line);
    return NULL;
  }


  out = GenericMatrix_alloc_std();

  if( out == NULL )
    return NULL;

  out->name = stringalloc(*splitstr);

  push_errormsg_stack("In reading matrix definition for %s",out->name);

  ckfree(base);


  while( get_watched_line(buffer,MAXLINE,ifp) != NULL) {
    chop_newline(buffer);

    if( strwhitestartcmp(buffer,"#",spacestr) == 0 )
      continue;

    if( only_whitespace(buffer,spacestr) == TRUE)
      continue;

    if( strstartcmp(buffer,"endmatrix") == 0)
      break;			
    if( strwhitestartcmp(buffer,"endmatrix",spacestr) == 0 ) {
      warn("endmatrix tag not flush to the start of the line. Ok, but not in the specification");
      break;
    }
    
    if( strstartcmp(buffer,"end") == 0 ) {
      warn("got an 'end' tag [%s] but expecting a 'endmatrix'. Considering this a parse failure");
      goto error;
    }
    
    
    /*** ok, proper parsing now ****/
    
    if( strwhitestartcmp(buffer,"collapse",spacestr) == 0 ) {
      cal = read_CollapsableLabel_line(buffer);
      if( cal == NULL ) {
	warn("Cannot read Collapsable label line. Ignoring collapsed label");
	continue;
      }
      add_cal_GenericMatrix(out,cal);
      continue;
    }

    else if( strwhitestartcmp(buffer,"extern",spacestr) == 0 ){
      ev = read_ExternVariable_line(buffer);
      if( ev == NULL ) {
	warn("unable to read an Extern line");
      }
      else add_ev_GenericMatrix(out,ev);
      continue;
    } else if( strwhitestartcmp(buffer,"state",spacestr) == 0) {
      
      
      if( strstr(buffer,"!special") != NULL )
	isspecial = TRUE;
      else if ( strstr(buffer,"SPECIAL") != NULL ) {
	warn("In state [%s], got a SPECIAL tag. This has been replaced with !special for consistency. Please change ;)",buffer);
	isspecial = TRUE;
      }
      else	isspecial = FALSE;
      
      /*********************************************/
      /* this function actually reads in the state */
      /* block                                     */
      /*********************************************/
      
      temp = read_CellState_line(buffer,ifp);
      
      /* read in state */
      
      if( temp == NULL ) {
	/* warning already issued, just chain */
	/*	warn("unable to read line for CellState in GenericMatrix %s - going to fail",out->name); */
	goto error;
      }
      if( isspecial == FALSE) {
	if( add_GenericMatrix(out,temp) == FALSE ) {
	  warn("Able to read - but unable to add - line for CellState in GeericMatrix %s - going to return now",out->name);
	  goto error;
	}
      } else	{
	if( add_spec_GenericMatrix(out,temp) == FALSE ) {
	  warn("Able to read - but unable to add - line for CellState in GeericMatrix %s - going to return now",out->name);
	  goto error;
	}
      }
      continue;		/* should not chop up line! */
    }


    /**** OK not a state line, hence a "processed" line ****/
    
    /* split up line, look at first word, decide what to do */
    /* - probably loop through the rest of the words to     */
    /* to read out - these are the for(splitstr++; etc loops*/ 
    
    
    base=splitstr=breakstring(buffer,spacestr);
    
    
    /* NB, base free'd at the end of the if/else switch */
    

    if( strwhitestartcmp(*splitstr,"query",spacestr) == 0) {
      if( out->query_name != NULL ) {
	log_full_error(WARNING,0,"This is the second time to specify a query - only one allowed: ignoring [%s]",buffer);
	continue;
      }

      /*** allocate's memory etc ready */
      out->query = StructElement_alloc();
      
      for(splitstr++;*splitstr;splitstr++) {
	  if( (runner=string_from_quoted_equality(*splitstr)) == NULL ) {
	    warn("You have specified a modifier [%s] to query but it has either no '=' sign or no quoted argument. The '=' character should be flush to both the tag and the quoted (using \")  argument",CKS(*splitstr));
	    continue;
	  }
	  if( strstartcmp(*splitstr,"name") == 0 ) {
	    out->query->name = runner;
	  }
	  else if ( strstartcmp(*splitstr,"field:name") ==  0) {
	    out->query_name = runner;
	  }
	  else if ( strstartcmp(*splitstr,"field:len") == 0) {
	    out->query_len = runner;
	  }
	  else if ( strstartcmp(*splitstr,"type") == 0) {
	    out->query->element_type = runner;
	  }
	  else	{
	    warn("Got modify %s=%s fine for tag query - but don't know what to do with it!",*splitstr,runner);
			  ckfree(runner);
	  }
	    } /* end of query modifers */
    } /* end of query if */
    else if( strwhitestartcmp(*splitstr,"target",spacestr) == 0)
      {
	if( out->target_name != NULL ) {
	  log_full_error(WARNING,0,"This is the second time to specify a target - only one allowed Ignoring [%s]",buffer);
	}
	/*** allocate's memory etc ready */
	
	
	out->target = StructElement_alloc();
	
	for(splitstr++;*splitstr;splitstr++) {
	    if( (runner=string_from_quoted_equality(*splitstr)) == NULL ) {
	      warn("You have specified a modifier [%s] to target but it has either no '=' sign or no quoted argument. The '=' character should be flush to both the tag and the quoted (using \")  argument",CKS(*splitstr));
	      continue;
	    }
	    if( strstartcmp(*splitstr,"name") == 0 ) {
		out->target->name = runner;
	    }
	    else if ( strstartcmp(*splitstr,"field:name") ==  0) {
		out->target_name = runner;
	      }
	    else if ( strstartcmp(*splitstr,"field:len") == 0)
	      {
		out->target_len = runner;
	      }
	    else if ( strstartcmp(*splitstr,"type") == 0)
	      {
		out->target->element_type = runner;
	      }
	    else	{
	      warn("Got modify %s=%s fine for tag target - but don't know what to do with it!",*splitstr,runner);
	      ckfree(runner);
	    }
	  } /* end of target modifers */
      } /* end of target if */
    else if( strstartcmp(*splitstr,"resource") == 0)
      {
	/*** allocate's memory etc ready */
	res_temp = StructElement_alloc();
			
	add_res_GenericMatrix(out,res_temp);

	for(splitstr++;*splitstr;splitstr++)
	  {
	    if( (runner=string_from_quoted_equality(*splitstr)) == NULL )
	      {
		warn("You have specified a modifier [%s] to resource but it has either no '=' sign or no quoted argument. The '=' character should be flush to both the tag and the quoted (using \")  argument",CKS(*splitstr));
		continue;
	      }
	    if( strstartcmp(*splitstr,"name") == 0 )
	      {
		res_temp->name = runner;
	      }
	    else if ( strstartcmp(*splitstr,"type") == 0)
	      {
		res_temp->element_type = runner;
	      }
	    else	{
	      warn("Got modifier %s=%s fine for tag resource - but don't know what to do with it!",*splitstr,runner);
	      ckfree(runner);
	    }
	  } /* end of resource modifers */
      } /* end of resource if */
    else if( strstartcmp(*splitstr,"globaldefaultscore") == 0)
      {
	warn("No need for globaldefaultscore lines anymore");
      }
    else if ( strstartcmp(*splitstr,"calcfunc") == 0) {
      if( *++splitstr == NULL ) {
	warn("Got a calcfunc tag with no function!");
      } else {
	out->calcfunc = stringalloc(*splitstr);
      }
    }
    else {
      warn("Could not understand line in matrix parse");
    }


		
    ckfree(base);
    }

  pop_errormsg_stack();
	
  return out;

  error :
  pop_errormsg_stack();

  out = free_GenericMatrix(out);
  return NULL;
}

/* Function:  read_ExternVariable_line(line)
 *
 * Descrip:    reads line like extern name="xxx" type="xxx"
 *
 *
 * Arg:        line [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [ExternVariable *]
 *
 */
# line 1481 "dyna2.dy"
ExternVariable * read_ExternVariable_line(char * line)
{
  ExternVariable * out;
  char ** base;
  char ** brk;
  char * nameq = NULL;
  char * tq = NULL;

  if( strwhitestartcmp(line,"extern",spacestr) != 0 ) {
    warn("Tried to pass read_ExternVariable_line without an extern tag. Nope!");
    return NULL;
  }

  base = brk = breakstring(line,spacestr);

  for(;*brk != NULL;brk++) {
    if( strcmp(*brk,"extern") == 0 )
      continue;
    if( strstartcmp(*brk,"name") == 0 ) {
      if( (nameq = string_from_quoted_equality(*brk)) == NULL) {
	warn("In reading extern line, got a name tag, but no argument. The tag should have no whitespace between aname and equals");
	continue;
      }
    } else if ( strstartcmp(*brk,"type") == 0 ) {
      if( (tq = string_from_quoted_equality(*brk)) == NULL) {
	warn("In reading extern line, got a type tag, but no argument. The tag should have no whitespace between the type and =");
	continue;
      }
    } else {
      warn("Did not understand tag [%s] in extern line",*brk);
    }
  }
  
  ckfree(base);


  out = ExternVariable_alloc();

  out->name = nameq;
  out->type = tq;

  return out;
}
  

# line 1526 "dyna2.dy"
CollapsableLabel * read_CollapsableLabel_line(char * line)
{
  CollapsableLabel * out;
  char * runner;
  char * run2;

  if( strwhitestartcmp(line,"collapse",spacestr) != 0 ) {
    warn("Tried to pass read_CollapsableLabel_line a no collapse line. Problem!");
    return NULL;
  }

  runner = strtok(line,spacestr);
  runner = strtok(NULL,spacestr);
  run2   = strtok(NULL,spacestr);
  if( runner == NULL || run2 == NULL) {
    warn("Collapsable line has no collapsable label. Ooops");
    return NULL;
  }

  out = CollapsableLabel_alloc();
  if( out == NULL)
    return NULL;

  out->query = stringalloc(runner);
  out->target = stringalloc(run2);
  
  return out;

}


# line 1557 "dyna2.dy"
CellState  * read_CellState_line(char * line,FILE * ifp)
{
  CellState * out;
  CellSource * temp;
  char buffer[MAXLINE];
  char ** base;
  char ** splitstr;
  char * runner;
	

  /*** allocate and die if no memory ***/
  /*** warnings should be issued in  ***/
  /*** alloc                         ***/


  /*** check state and name and Loop through the current line given ****/

  base=splitstr=breakstring(line,spacestr);


  if( strcmp(*splitstr,"state") != 0) {
    log_full_error(WARNING,0,"In parsing the line starting %s it had no state tag!",line);
    return NULL;
  }

  splitstr++;


  if( *splitstr == NULL ) {
    log_full_error(WARNING,0,"In parsing the line source [%s] there was no state name",line);
    return NULL;
  }


  out = CellState_alloc_std();

  if( out == NULL )
    return NULL;

  out->name = stringalloc(*splitstr);

  push_errormsg_stack("In Reading state %s",out->name);


  /*** looping through first line ***/

  for(splitstr++;*splitstr;splitstr++) {
    if( strstartcmp(*splitstr,"offi") == 0) {
      runner=string_from_quoted_equality(*splitstr);
      if( runner == NULL ) {
	warn("Unable to read offi in state %s. The tag should look like offi=\"<number>\" with no whitespace",out->name);
	continue;
      }
      out->offi = atoi(runner);
      ckfree(runner);
    }
    else if ( strstartcmp(*splitstr,"offj") == 0) {
      runner=string_from_quoted_equality(*splitstr);
      if( runner == NULL ) {
	warn("Unable to read offj in state %s. The tag should look like offj=\"<number>\" with no whitespace",out->name);
	continue;
      }
      out->offj = atoi(runner);
      ckfree(runner);
    }
    else if ( strstartcmp(*splitstr,"defscore") == 0)
      {
	warn("No need for defscore lines anymore");
      }			
    else if ( strstartcmp(*splitstr,"calc") == 0)
      {
	runner=string_from_quoted_equality(*splitstr);
	if( runner == NULL ) {
	    warn("Unable to read source independent calc line in state %s. Remember that there should be no white space between calc and the \"string\", ie calc=\"gap\" ",out->name);
	    continue;
	  }
	out->calc_expr=runner;
      }
    else if ( strstartcmp(*splitstr,"!special") == 0) {
      out->is_special_i=TRUE;
    }
    else if ( strstartcmp(*splitstr,"SPECIAL") == 0) {
      out->is_special_i=TRUE;
    }
    
    else if ( strstartcmp(*splitstr,"!end") == 0) {
      out->is_end = TRUE;
    }
    else if ( strstartcmp(*splitstr,"!start") == 0) {
      out->is_start = TRUE;
    }
    else {
      warn("Parse error in state %s - cannot make sense of %s",out->name,*splitstr);
    }
  }
  
  ckfree(base);


  while( get_watched_line(buffer,MAXLINE,ifp) != NULL) {
    chop_newline(buffer);

    if( strwhitestartcmp(buffer,"#",spacestr) == 0 )
      continue;
    
    if( strwhitestartcmp(buffer,"endstate",spacestr) == 0 )
      break;

    if( only_whitespace(buffer,spacestr) == TRUE)
      continue;

    if( strwhitestartcmp(buffer,"end",spacestr) == 0 ) {
      warn("Got an end line [%s] but expecting endstate. Will fail.",buffer);
      goto error;
    }
    if( strwhitestartcmp(buffer,"state",spacestr) == 0 ) {
      warn("Got the line [%s], a state start line inside a state. Expect you forgot an endstate.",buffer);
      goto error;
    }


    if( strwhitestartcmp(buffer,"query_label",spacestr) == 0 ) {
      base = splitstr = breakstring(buffer,spacestr);
      splitstr++;
      if( *splitstr == NULL ) {
	warn("Picked up query_label tag but no query label in state %s",out->name);
      }
      else out->query_label = stringalloc(*splitstr);
      ckfree(base);
      continue;
    }


    if( strwhitestartcmp(buffer,"target_label",spacestr) == 0 ) {
      base = splitstr = breakstring(buffer,spacestr);
      splitstr++;
      if( *splitstr == NULL ) {
	warn("Picked up target_label tag but no target label in state %s",out->name);
      }
      else out->target_label = stringalloc(*splitstr);
      ckfree(base);
      continue;
    }

    if( strwhitestartcmp(buffer,"calc",spacestr) == 0 ) {
      base = splitstr = breakstring(buffer,spacestr);
      if( out->calc_expr != NULL ) {
	warn("Already picked up a calc line [%s]. Replacing with %s\n",out->calc_expr,runner);
	ckfree(out->calc_expr);
      }
      
      out->calc_expr = string_from_quoted_equality(*splitstr);
      ckfree(base);
      continue;
    }

    else if( strstr(buffer,"source") == NULL) {
      warn("Parse error in reading state %s - cannot interpret [%s]",out->name,buffer);
      continue;
    }

    temp = read_CellSource_line(buffer,ifp);
    
    /*		fprintf(stderr,"Have read source line!\n"); */

    if( temp == NULL ) {
      /* warning already issued, just chain back up */
      /* warn("unable to read line for CellSource in CellState %s - going to fail parser",out->name); */
      goto error;
    }

    if( add_CellState(out,temp) == FALSE ) {
      warn("Able to read - but unable to add - line for CellSource in CellState %s - going to return now",out->name);
      goto error;
    }
    
  }
	
  pop_errormsg_stack();

  return out;

  error :
    
  pop_errormsg_stack();
  out = free_CellState(out);
  return NULL;

}

# line 1747 "dyna2.dy"
int source_bit2pos(int bit)
{
  if( bit == 0 ) 
    return SOURCE_POS_ALL;

  if( (bit & SOURCE_TOP_BIT) ) {
    if( (bit & SOURCE_LEFT_BIT) ) {
      return SOURCE_POS_TOPLEFT;
    }
    if( (bit & ~SOURCE_TOP_BIT) != 0 ) {
	return bit;
    }
    else 
      return SOURCE_POS_TOP;
  
  }

  if( (bit & SOURCE_LEFT_BIT) ) {
    if( (bit & ~SOURCE_LEFT_BIT) != 0 )
      return bit;
    else return SOURCE_POS_LEFT;
  }


  if( (bit & SOURCE_BOTTOM_BIT) ) {
    if( (bit & SOURCE_RIGHT_BIT) ) {
      return SOURCE_POS_BOTTOMRIGHT;
    }
    if( (bit & ~SOURCE_BOTTOM_BIT) != 0 ) {
	return bit;
    }
    else 
      return SOURCE_POS_BOTTOM;
  
  }

  if( (bit & SOURCE_RIGHT_BIT) ) {
    if( (bit & ~SOURCE_RIGHT_BIT) != 0 )
      return bit;
    else return SOURCE_POS_RIGHT;
  }

  warn("Got a HIDEOUS error in source_bit2pos");
  return SOURCE_POS_ALL;
}


# line 1794 "dyna2.dy"
CellSource * read_CellSource_line(char * line,FILE * ifp)
{
  CellSource * out;
  char buffer[MAXLINE];
  char * runner;	
  char * temp;
  char ** base;
  char ** splitstr;
  int posbits = 0; /** use with SOURCE_TOP_BIT etc **/

	
  /* parse first line */

  base=splitstr=breakstring(line,spacestr);

  if( strcmp(*splitstr,"source") != 0)  {
    warn("In parsing the line starting %s it had no source tag!",line);
    return NULL;
  }
  splitstr++;
  
  if( *splitstr == NULL )  {
    warn("In parsing the line source [%s] there was no source tag",line);
    return NULL;
  }
  
  out = CellSource_alloc();

  if( out == NULL)
    return NULL;

  out->state_source = stringalloc(*splitstr);


  push_errormsg_stack("In reading source %s",out->state_source);

  for(splitstr++;*splitstr;splitstr++)
    {
      if( strstartcmp(*splitstr,"offi") == 0)
	{
	  runner=string_from_quoted_equality(*splitstr);
	  if( runner == NULL )
	    {
	      warn("Unable to read offi in Source %s",line);
	      continue;
	    }
	  out->offi = atoi(runner);
	  ckfree(runner);
	}
      else if ( strstartcmp(*splitstr,"offj") == 0)
	{
	  runner=string_from_quoted_equality(*splitstr);
	  if( runner == NULL )
	    {
	      warn("Unable to read offi in Source %s",line);
	      continue;
	    }
	  out->offj = atoi(runner);
	  ckfree(runner);
	}
      else if ( strcmp(*splitstr,"!top") == 0 ) {
	posbits = (posbits | SOURCE_TOP_BIT);
      }
      else if ( strcmp(*splitstr,"!left") == 0 ) {
	posbits = (posbits | SOURCE_LEFT_BIT);
      }
      else if ( strcmp(*splitstr,"!right") == 0 ) {
	posbits = (posbits | SOURCE_RIGHT_BIT);
      }
      else if ( strcmp(*splitstr,"!bottom") == 0 ) {
	posbits = (posbits | SOURCE_BOTTOM_BIT);
      }
      
      else	{
		  warn("Parse error in source %s - cannot make sense of %s",line,*splitstr);
		}
    }

  out->position = source_bit2pos(posbits);
  ckfree(base);

  while( get_watched_line(buffer,MAXLINE,ifp) != NULL) {
    chop_newline(buffer);


    if( strwhitestartcmp(buffer,"#",spacestr) == 0 )
      continue;
    if( only_whitespace(buffer,spacestr) == TRUE)
      continue;
    
    if( strwhitestartcmp(buffer,"endsource",spacestr) == 0 )
      break;
    if( strwhitestartcmp(buffer,"end",spacestr) == 0 ) {
      warn("you have a end line [%s] but expecting an endsource line. Will fail",buffer);
      goto error;
    }

    if( strwhitestartcmp(buffer,"source",spacestr) == 0 ) {
      warn("In reading a source, got a source tag [%s]. Expect you forgot an endsource. Will fail",buffer);
      goto error;
    }
    
    
    if( strwhitestartcmp(buffer,"query_label",spacestr) == 0 ) {
      base = splitstr = breakstring(buffer,spacestr);
      splitstr++;
      if( *splitstr == NULL ) {
	warn("Picked up query_label tag but no query label in source %s",out->state_source);
      }
      else out->query_label = stringalloc(*splitstr);
      ckfree(base);
      continue;
    }


    else if( strwhitestartcmp(buffer,"target_label",spacestr) == 0 ) {
      base = splitstr = breakstring(buffer,spacestr);
      splitstr++;
      if( *splitstr == NULL ) {
	warn("Picked up target_label tag but no target label in source %s",out->state_source);
      }
      else out->target_label = stringalloc(*splitstr);
      ckfree(base);
      continue;
    }
    
    else if( strwhitestartcmp(buffer,"calc",spacestr) == 0 ) {
      temp = read_calc_line(buffer);

      if( temp == NULL ) {
	warn("unable to read calc line in Source %s - going to return now",out->state_source);
	return NULL;
      }
      out->calc_expr = stringalloc(temp);
    }
    
    else {
      warn("Cannot understand the line [%s] on CellSource",buffer);
    }
    
  }
  
  pop_errormsg_stack();
  
  return out;

  error :

  pop_errormsg_stack();

    out = free_CellSource(out);
  return NULL;
}

# line 1948 "dyna2.dy"
char * read_calc_line(char * buffer)
{
  char * runner;
  /* get to = */
  for(;*buffer && *buffer != '=';buffer++)
    ;
  /* get to " */
  for(;*buffer && *buffer != '"';buffer++)
    ;
  if( *buffer == '\0' ) {
    warn("In reading calc line string, got to the end of the buffer without gettint to a \"");
    return NULL;
  }

  runner = buffer+1;
  for(buffer++;*buffer && *buffer != '"';buffer++)
    ;
  *buffer='\0';
  return runner;
}

# line 1913 "dyna2.c"
/* Function:  hard_link_CellSignature(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [CellSignature *]
 *
 * Return [UNKN ]  Undocumented return value [CellSignature *]
 *
 */
CellSignature * hard_link_CellSignature(CellSignature * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a CellSignature object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  CellSignature_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [CellSignature *]
 *
 */
CellSignature * CellSignature_alloc(void) 
{
    CellSignature * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(CellSignature *) ckalloc (sizeof(CellSignature))) == NULL)  {  
      warn("CellSignature_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->offi = 0;   
    out->offj = 0;   


    return out;  
}    


/* Function:  free_CellSignature(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [CellSignature *]
 *
 * Return [UNKN ]  Undocumented return value [CellSignature *]
 *
 */
CellSignature * free_CellSignature(CellSignature * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a CellSignature obj. Should be trappable"); 
      return NULL;   
      }  


#ifdef PTHREAD   
    assert(pthread_mutex_lock(&(obj->dynamite_mutex)) == 0); 
#endif   
    if( obj->dynamite_hard_link > 1)     {  
      return_early = 1;  
      obj->dynamite_hard_link--; 
      }  
#ifdef PTHREAD   
    assert(pthread_mutex_unlock(&(obj->dynamite_mutex)) == 0);   
#endif   
    if( return_early == 1)   
      return NULL;   


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_CellSignatureSet(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_CellSignatureSet
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [CellSignature **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_CellSignatureSet(CellSignature ** list,int i,int j)  
{
    CellSignature * temp;    
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_CellSignatureSet(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_CellSignatureSet which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [CellSignature **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_CellSignatureSet(CellSignature ** list,int left,int right,int (*comp)(CellSignature * ,CellSignature * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_CellSignatureSet(list,left,(left+right)/2); 
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_CellSignatureSet (list,++last,i);   
      }  
    swap_CellSignatureSet (list,left,last);  
    qsort_CellSignatureSet(list,left,last-1,comp);   
    qsort_CellSignatureSet(list,last+1,right,comp);  
}    


/* Function:  sort_CellSignatureSet(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_CellSignatureSet
 *
 *
 * Arg:         obj [UNKN ] Object containing list [CellSignatureSet *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_CellSignatureSet(CellSignatureSet * obj,int (*comp)(CellSignature *, CellSignature *)) 
{
    qsort_CellSignatureSet(obj->sig,0,obj->len-1,comp);  
    return;  
}    


/* Function:  expand_CellSignatureSet(obj,len)
 *
 * Descrip:    Really an internal function for add_CellSignatureSet
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [CellSignatureSet *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_CellSignatureSet(CellSignatureSet * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_CellSignatureSet called with no need");   
      return TRUE;   
      }  


    if( (obj->sig = (CellSignature ** ) ckrealloc (obj->sig,sizeof(CellSignature *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_CellSignatureSet, returning FALSE"); 
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_CellSignatureSet(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [CellSignatureSet *]
 * Arg:        add [OWNER] Object to add to the list [CellSignature *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_CellSignatureSet(CellSignatureSet * obj,CellSignature * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_CellSignatureSet(obj,obj->len + CellSignatureSetLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->sig[obj->len++]=add;    
    return TRUE; 
}    


/* Function:  flush_CellSignatureSet(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [CellSignatureSet *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_CellSignatureSet(CellSignatureSet * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->sig[i] != NULL)   {  
        free_CellSignature(obj->sig[i]); 
        obj->sig[i] = NULL;  
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  CellSignatureSet_alloc_std(void)
 *
 * Descrip:    Equivalent to CellSignatureSet_alloc_len(CellSignatureSetLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [CellSignatureSet *]
 *
 */
CellSignatureSet * CellSignatureSet_alloc_std(void) 
{
    return CellSignatureSet_alloc_len(CellSignatureSetLISTLENGTH);   
}    


/* Function:  CellSignatureSet_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [CellSignatureSet *]
 *
 */
CellSignatureSet * CellSignatureSet_alloc_len(int len) 
{
    CellSignatureSet * out; /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = CellSignatureSet_alloc()) == NULL) 
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->sig = (CellSignature ** ) ckcalloc (len,sizeof(CellSignature *))) == NULL)  {  
      warn("Warning, ckcalloc failed in CellSignatureSet_alloc_len");    
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_CellSignatureSet(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [CellSignatureSet *]
 *
 * Return [UNKN ]  Undocumented return value [CellSignatureSet *]
 *
 */
CellSignatureSet * hard_link_CellSignatureSet(CellSignatureSet * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a CellSignatureSet object: passed a NULL object");    
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  CellSignatureSet_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [CellSignatureSet *]
 *
 */
CellSignatureSet * CellSignatureSet_alloc(void) 
{
    CellSignatureSet * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(CellSignatureSet *) ckalloc (sizeof(CellSignatureSet))) == NULL)    {  
      warn("CellSignatureSet_alloc failed ");    
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->sig = NULL; 
    out->len = out->maxlen = 0;  


    return out;  
}    


/* Function:  free_CellSignatureSet(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [CellSignatureSet *]
 *
 * Return [UNKN ]  Undocumented return value [CellSignatureSet *]
 *
 */
CellSignatureSet * free_CellSignatureSet(CellSignatureSet * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a CellSignatureSet obj. Should be trappable");  
      return NULL;   
      }  


#ifdef PTHREAD   
    assert(pthread_mutex_lock(&(obj->dynamite_mutex)) == 0); 
#endif   
    if( obj->dynamite_hard_link > 1)     {  
      return_early = 1;  
      obj->dynamite_hard_link--; 
      }  
#ifdef PTHREAD   
    assert(pthread_mutex_unlock(&(obj->dynamite_mutex)) == 0);   
#endif   
    if( return_early == 1)   
      return NULL;   
    if( obj->sig != NULL)    {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->sig[i] != NULL) 
          free_CellSignature(obj->sig[i]);   
        }  
      ckfree(obj->sig);  
      }  


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_CellSource(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [CellSource *]
 *
 * Return [UNKN ]  Undocumented return value [CellSource *]
 *
 */
CellSource * hard_link_CellSource(CellSource * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a CellSource object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  CellSource_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [CellSource *]
 *
 */
CellSource * CellSource_alloc(void) 
{
    CellSource * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(CellSource *) ckalloc (sizeof(CellSource))) == NULL)    {  
      warn("CellSource_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->state_source = NULL;    
    out->offi = -1;  
    out->offj = -1;  
    out->calc_expr = NULL;   
    out->source_expr = NULL; 
    out->etr = NULL; 
    out->isspecial = FALSE;  
    out->query_label = NULL; 
    out->target_label = NULL;    
    out->position = SOURCE_POS_ALL;  
    out->trans_no = 0;   
    out->from_state_no = 0;  


    return out;  
}    


/* Function:  free_CellSource(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [CellSource *]
 *
 * Return [UNKN ]  Undocumented return value [CellSource *]
 *
 */
CellSource * free_CellSource(CellSource * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a CellSource obj. Should be trappable");    
      return NULL;   
      }  


#ifdef PTHREAD   
    assert(pthread_mutex_lock(&(obj->dynamite_mutex)) == 0); 
#endif   
    if( obj->dynamite_hard_link > 1)     {  
      return_early = 1;  
      obj->dynamite_hard_link--; 
      }  
#ifdef PTHREAD   
    assert(pthread_mutex_unlock(&(obj->dynamite_mutex)) == 0);   
#endif   
    if( return_early == 1)   
      return NULL;   
    if( obj->state_source != NULL)   
      ckfree(obj->state_source);     
    if( obj->calc_expr != NULL)  
      ckfree(obj->calc_expr);    
    if( obj->source_expr != NULL)    
      ckfree(obj->source_expr);  
    if( obj->etr != NULL)    
      free_ExprTree(obj->etr);   
    if( obj->query_label != NULL)    
      ckfree(obj->query_label);  
    if( obj->target_label != NULL)   
      ckfree(obj->target_label);     


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_CellState(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_CellState
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [CellSource **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_CellState(CellSource ** list,int i,int j)  
{
    CellSource * temp;   
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_CellState(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_CellState which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [CellSource **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_CellState(CellSource ** list,int left,int right,int (*comp)(CellSource * ,CellSource * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_CellState(list,left,(left+right)/2);    
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_CellState (list,++last,i);  
      }  
    swap_CellState (list,left,last); 
    qsort_CellState(list,left,last-1,comp);  
    qsort_CellState(list,last+1,right,comp); 
}    


/* Function:  sort_CellState(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_CellState
 *
 *
 * Arg:         obj [UNKN ] Object containing list [CellState *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_CellState(CellState * obj,int (*comp)(CellSource *, CellSource *)) 
{
    qsort_CellState(obj->source,0,obj->len-1,comp);  
    return;  
}    


/* Function:  expand_CellState(obj,len)
 *
 * Descrip:    Really an internal function for add_CellState
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [CellState *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_CellState(CellState * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_CellState called with no need");  
      return TRUE;   
      }  


    if( (obj->source = (CellSource ** ) ckrealloc (obj->source,sizeof(CellSource *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_CellState, returning FALSE");    
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_CellState(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [CellState *]
 * Arg:        add [OWNER] Object to add to the list [CellSource *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_CellState(CellState * obj,CellSource * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_CellState(obj,obj->len + CellStateLISTLENGTH) == FALSE) 
        return FALSE;    
      }  


    obj->source[obj->len++]=add; 
    return TRUE; 
}    


/* Function:  flush_CellState(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [CellState *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_CellState(CellState * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->source[i] != NULL)    {  
        free_CellSource(obj->source[i]); 
        obj->source[i] = NULL;   
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  CellState_alloc_std(void)
 *
 * Descrip:    Equivalent to CellState_alloc_len(CellStateLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [CellState *]
 *
 */
CellState * CellState_alloc_std(void) 
{
    return CellState_alloc_len(CellStateLISTLENGTH); 
}    


/* Function:  CellState_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [CellState *]
 *
 */
CellState * CellState_alloc_len(int len) 
{
    CellState * out;/* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = CellState_alloc()) == NULL)    
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->source = (CellSource ** ) ckcalloc (len,sizeof(CellSource *))) == NULL) {  
      warn("Warning, ckcalloc failed in CellState_alloc_len");   
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_CellState(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [CellState *]
 *
 * Return [UNKN ]  Undocumented return value [CellState *]
 *
 */
CellState * hard_link_CellState(CellState * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a CellState object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  CellState_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [CellState *]
 *
 */
CellState * CellState_alloc(void) 
{
    CellState * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(CellState *) ckalloc (sizeof(CellState))) == NULL)  {  
      warn("CellState_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->name = NULL;    
    out->def_score = NULL;   
    out->calc_expr = NULL;   
    out->source_expr = NULL; 
    out->offi = 0;   
    out->offj = 0;   
    out->is_special_i = FALSE;   
    out->is_special_j = FALSE;   
    out->is_end = FALSE; 
    out->is_start = FALSE;   
    out->specialtospecial = FALSE;   
    out->source = NULL;  
    out->len = out->maxlen = 0;  
    out->query_char = NULL;  
    out->target_char = NULL; 
    out->footprint_start = 0;    
    out->footprint_end = 1;  
    out->query_label = NULL; 
    out->target_label = NULL;    
    out->position = SOURCE_POS_ALL;  
    out->etr = NULL; 
    out->state_number = -1;  


    return out;  
}    


/* Function:  free_CellState(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [CellState *]
 *
 * Return [UNKN ]  Undocumented return value [CellState *]
 *
 */
CellState * free_CellState(CellState * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a CellState obj. Should be trappable"); 
      return NULL;   
      }  


#ifdef PTHREAD   
    assert(pthread_mutex_lock(&(obj->dynamite_mutex)) == 0); 
#endif   
    if( obj->dynamite_hard_link > 1)     {  
      return_early = 1;  
      obj->dynamite_hard_link--; 
      }  
#ifdef PTHREAD   
    assert(pthread_mutex_unlock(&(obj->dynamite_mutex)) == 0);   
#endif   
    if( return_early == 1)   
      return NULL;   
    if( obj->name != NULL)   
      ckfree(obj->name);     
    if( obj->def_score != NULL)  
      ckfree(obj->def_score);    
    if( obj->calc_expr != NULL)  
      ckfree(obj->calc_expr);    
    if( obj->source_expr != NULL)    
      ckfree(obj->source_expr);  
    if( obj->source != NULL) {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->source[i] != NULL)  
          free_CellSource(obj->source[i]);   
        }  
      ckfree(obj->source);   
      }  
    if( obj->query_char != NULL) 
      ckfree(obj->query_char);   
    if( obj->target_char != NULL)    
      ckfree(obj->target_char);  
    if( obj->query_label != NULL)    
      ckfree(obj->query_label);  
    if( obj->target_label != NULL)   
      ckfree(obj->target_label);     
    if( obj->etr != NULL)    
      free_ExprTree(obj->etr);   


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_CollapsableLabel(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [CollapsableLabel *]
 *
 * Return [UNKN ]  Undocumented return value [CollapsableLabel *]
 *
 */
CollapsableLabel * hard_link_CollapsableLabel(CollapsableLabel * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a CollapsableLabel object: passed a NULL object");    
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  CollapsableLabel_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [CollapsableLabel *]
 *
 */
CollapsableLabel * CollapsableLabel_alloc(void) 
{
    CollapsableLabel * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(CollapsableLabel *) ckalloc (sizeof(CollapsableLabel))) == NULL)    {  
      warn("CollapsableLabel_alloc failed ");    
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->query = NULL;   
    out->target = NULL;  


    return out;  
}    


/* Function:  free_CollapsableLabel(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [CollapsableLabel *]
 *
 * Return [UNKN ]  Undocumented return value [CollapsableLabel *]
 *
 */
CollapsableLabel * free_CollapsableLabel(CollapsableLabel * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a CollapsableLabel obj. Should be trappable");  
      return NULL;   
      }  


#ifdef PTHREAD   
    assert(pthread_mutex_lock(&(obj->dynamite_mutex)) == 0); 
#endif   
    if( obj->dynamite_hard_link > 1)     {  
      return_early = 1;  
      obj->dynamite_hard_link--; 
      }  
#ifdef PTHREAD   
    assert(pthread_mutex_unlock(&(obj->dynamite_mutex)) == 0);   
#endif   
    if( return_early == 1)   
      return NULL;   
    if( obj->query != NULL)  
      ckfree(obj->query);    
    if( obj->target != NULL) 
      ckfree(obj->target);   


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_ExternVariable(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [ExternVariable *]
 *
 * Return [UNKN ]  Undocumented return value [ExternVariable *]
 *
 */
ExternVariable * hard_link_ExternVariable(ExternVariable * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a ExternVariable object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  ExternVariable_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ExternVariable *]
 *
 */
ExternVariable * ExternVariable_alloc(void) 
{
    ExternVariable * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(ExternVariable *) ckalloc (sizeof(ExternVariable))) == NULL)    {  
      warn("ExternVariable_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->name = NULL;    
    out->type = NULL;    


    return out;  
}    


/* Function:  free_ExternVariable(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [ExternVariable *]
 *
 * Return [UNKN ]  Undocumented return value [ExternVariable *]
 *
 */
ExternVariable * free_ExternVariable(ExternVariable * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a ExternVariable obj. Should be trappable");    
      return NULL;   
      }  


#ifdef PTHREAD   
    assert(pthread_mutex_lock(&(obj->dynamite_mutex)) == 0); 
#endif   
    if( obj->dynamite_hard_link > 1)     {  
      return_early = 1;  
      obj->dynamite_hard_link--; 
      }  
#ifdef PTHREAD   
    assert(pthread_mutex_unlock(&(obj->dynamite_mutex)) == 0);   
#endif   
    if( return_early == 1)   
      return NULL;   
    if( obj->name != NULL)   
      ckfree(obj->name);     
    if( obj->type != NULL)   
      ckfree(obj->type);     


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_GenericMatrix(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_GenericMatrix
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [CellState **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_GenericMatrix(CellState ** list,int i,int j)  
{
    CellState * temp;    
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_GenericMatrix(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_GenericMatrix which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [CellState **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_GenericMatrix(CellState ** list,int left,int right,int (*comp)(CellState * ,CellState * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_GenericMatrix(list,left,(left+right)/2);    
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_GenericMatrix (list,++last,i);  
      }  
    swap_GenericMatrix (list,left,last); 
    qsort_GenericMatrix(list,left,last-1,comp);  
    qsort_GenericMatrix(list,last+1,right,comp); 
}    


/* Function:  sort_GenericMatrix(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_GenericMatrix
 *
 *
 * Arg:         obj [UNKN ] Object containing list [GenericMatrix *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_GenericMatrix(GenericMatrix * obj,int (*comp)(CellState *, CellState *)) 
{
    qsort_GenericMatrix(obj->state,0,obj->len-1,comp);   
    return;  
}    


/* Function:  expand_GenericMatrix(obj,len)
 *
 * Descrip:    Really an internal function for add_GenericMatrix
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [GenericMatrix *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_GenericMatrix(GenericMatrix * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_GenericMatrix called with no need");  
      return TRUE;   
      }  


    if( (obj->state = (CellState ** ) ckrealloc (obj->state,sizeof(CellState *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_GenericMatrix, returning FALSE");    
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_GenericMatrix(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [GenericMatrix *]
 * Arg:        add [OWNER] Object to add to the list [CellState *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_GenericMatrix(GenericMatrix * obj,CellState * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_GenericMatrix(obj,obj->len + GenericMatrixLISTLENGTH) == FALSE) 
        return FALSE;    
      }  


    obj->state[obj->len++]=add;  
    return TRUE; 
}    


/* Function:  flush_GenericMatrix(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [GenericMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_GenericMatrix(GenericMatrix * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->state[i] != NULL) {  
        free_CellState(obj->state[i]);   
        obj->state[i] = NULL;    
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  swap_spec_GenericMatrix(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_spec_GenericMatrix
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [CellState **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_spec_GenericMatrix(CellState ** list,int i,int j)  
{
    CellState * temp;    
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_spec_GenericMatrix(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_spec_GenericMatrix which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [CellState **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_spec_GenericMatrix(CellState ** list,int left,int right,int (*comp)(CellState * ,CellState * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_spec_GenericMatrix(list,left,(left+right)/2);   
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_spec_GenericMatrix (list,++last,i); 
      }  
    swap_spec_GenericMatrix (list,left,last);    
    qsort_spec_GenericMatrix(list,left,last-1,comp); 
    qsort_spec_GenericMatrix(list,last+1,right,comp);    
}    


/* Function:  sort_spec_GenericMatrix(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_spec_GenericMatrix
 *
 *
 * Arg:         obj [UNKN ] Object containing list [GenericMatrix *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_spec_GenericMatrix(GenericMatrix * obj,int (*comp)(CellState *, CellState *)) 
{
    qsort_spec_GenericMatrix(obj->special,0,obj->spec_len-1,comp);   
    return;  
}    


/* Function:  expand_spec_GenericMatrix(obj,len)
 *
 * Descrip:    Really an internal function for add_spec_GenericMatrix
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [GenericMatrix *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_spec_GenericMatrix(GenericMatrix * obj,int len) 
{


    if( obj->spec_maxlen > obj->spec_len )   {  
      warn("expand_GenericMatrixspec_ called with no need"); 
      return TRUE;   
      }  


    if( (obj->special = (CellState ** ) ckrealloc (obj->special,sizeof(CellState *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_GenericMatrix, returning FALSE");    
      return FALSE;  
      }  
    obj->spec_maxlen = len;  
    return TRUE; 
}    


/* Function:  add_spec_GenericMatrix(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [GenericMatrix *]
 * Arg:        add [OWNER] Object to add to the list [CellState *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_spec_GenericMatrix(GenericMatrix * obj,CellState * add) 
{
    if( obj->spec_len >= obj->spec_maxlen)   {  
      if( expand_spec_GenericMatrix(obj,obj->spec_len + GenericMatrixLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->special[obj->spec_len++]=add;   
    return TRUE; 
}    


/* Function:  flush_spec_GenericMatrix(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [GenericMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_spec_GenericMatrix(GenericMatrix * obj) 
{
    int i;   


    for(i=0;i<obj->spec_len;i++) { /*for i over list length*/ 
      if( obj->special[i] != NULL)   {  
        free_CellState(obj->special[i]); 
        obj->special[i] = NULL;  
        }  
      } /* end of for i over list length */ 


    obj->spec_len = 0;   
    return i;    
}    


/* Function:  swap_res_GenericMatrix(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_res_GenericMatrix
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [StructElement  **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_res_GenericMatrix(StructElement  ** list,int i,int j)  
{
    StructElement  * temp;   
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_res_GenericMatrix(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_res_GenericMatrix which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [StructElement  **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_res_GenericMatrix(StructElement  ** list,int left,int right,int (*comp)(StructElement  * ,StructElement  * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_res_GenericMatrix(list,left,(left+right)/2);    
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_res_GenericMatrix (list,++last,i);  
      }  
    swap_res_GenericMatrix (list,left,last); 
    qsort_res_GenericMatrix(list,left,last-1,comp);  
    qsort_res_GenericMatrix(list,last+1,right,comp); 
}    


/* Function:  sort_res_GenericMatrix(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_res_GenericMatrix
 *
 *
 * Arg:         obj [UNKN ] Object containing list [GenericMatrix *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_res_GenericMatrix(GenericMatrix * obj,int (*comp)(StructElement  *, StructElement  *)) 
{
    qsort_res_GenericMatrix(obj->resource,0,obj->res_len-1,comp);    
    return;  
}    


/* Function:  expand_res_GenericMatrix(obj,len)
 *
 * Descrip:    Really an internal function for add_res_GenericMatrix
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [GenericMatrix *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_res_GenericMatrix(GenericMatrix * obj,int len) 
{


    if( obj->res_maxlen > obj->res_len )     {  
      warn("expand_GenericMatrixres_ called with no need");  
      return TRUE;   
      }  


    if( (obj->resource = (StructElement  ** ) ckrealloc (obj->resource,sizeof(StructElement  *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_GenericMatrix, returning FALSE");    
      return FALSE;  
      }  
    obj->res_maxlen = len;   
    return TRUE; 
}    


/* Function:  add_res_GenericMatrix(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [GenericMatrix *]
 * Arg:        add [OWNER] Object to add to the list [StructElement  *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_res_GenericMatrix(GenericMatrix * obj,StructElement  * add) 
{
    if( obj->res_len >= obj->res_maxlen) {  
      if( expand_res_GenericMatrix(obj,obj->res_len + GenericMatrixLISTLENGTH) == FALSE) 
        return FALSE;    
      }  


    obj->resource[obj->res_len++]=add;   
    return TRUE; 
}    


/* Function:  flush_res_GenericMatrix(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [GenericMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_res_GenericMatrix(GenericMatrix * obj) 
{
    int i;   


    for(i=0;i<obj->res_len;i++)  { /*for i over list length*/ 
      if( obj->resource[i] != NULL)  {  
        free_StructElement(obj->resource[i]);    
        obj->resource[i] = NULL; 
        }  
      } /* end of for i over list length */ 


    obj->res_len = 0;    
    return i;    
}    


/* Function:  swap_cal_GenericMatrix(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_cal_GenericMatrix
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [CollapsableLabel **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_cal_GenericMatrix(CollapsableLabel ** list,int i,int j)  
{
    CollapsableLabel * temp; 
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_cal_GenericMatrix(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_cal_GenericMatrix which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [CollapsableLabel **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_cal_GenericMatrix(CollapsableLabel ** list,int left,int right,int (*comp)(CollapsableLabel * ,CollapsableLabel * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_cal_GenericMatrix(list,left,(left+right)/2);    
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_cal_GenericMatrix (list,++last,i);  
      }  
    swap_cal_GenericMatrix (list,left,last); 
    qsort_cal_GenericMatrix(list,left,last-1,comp);  
    qsort_cal_GenericMatrix(list,last+1,right,comp); 
}    


/* Function:  sort_cal_GenericMatrix(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_cal_GenericMatrix
 *
 *
 * Arg:         obj [UNKN ] Object containing list [GenericMatrix *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_cal_GenericMatrix(GenericMatrix * obj,int (*comp)(CollapsableLabel *, CollapsableLabel *)) 
{
    qsort_cal_GenericMatrix(obj->cal,0,obj->cal_len-1,comp); 
    return;  
}    


/* Function:  expand_cal_GenericMatrix(obj,len)
 *
 * Descrip:    Really an internal function for add_cal_GenericMatrix
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [GenericMatrix *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_cal_GenericMatrix(GenericMatrix * obj,int len) 
{


    if( obj->cal_maxlen > obj->cal_len )     {  
      warn("expand_GenericMatrixcal_ called with no need");  
      return TRUE;   
      }  


    if( (obj->cal = (CollapsableLabel ** ) ckrealloc (obj->cal,sizeof(CollapsableLabel *)*len)) == NULL)     {  
      warn("ckrealloc failed for expand_GenericMatrix, returning FALSE");    
      return FALSE;  
      }  
    obj->cal_maxlen = len;   
    return TRUE; 
}    


/* Function:  add_cal_GenericMatrix(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [GenericMatrix *]
 * Arg:        add [OWNER] Object to add to the list [CollapsableLabel *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_cal_GenericMatrix(GenericMatrix * obj,CollapsableLabel * add) 
{
    if( obj->cal_len >= obj->cal_maxlen) {  
      if( expand_cal_GenericMatrix(obj,obj->cal_len + GenericMatrixLISTLENGTH) == FALSE) 
        return FALSE;    
      }  


    obj->cal[obj->cal_len++]=add;    
    return TRUE; 
}    


/* Function:  flush_cal_GenericMatrix(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [GenericMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_cal_GenericMatrix(GenericMatrix * obj) 
{
    int i;   


    for(i=0;i<obj->cal_len;i++)  { /*for i over list length*/ 
      if( obj->cal[i] != NULL)   {  
        free_CollapsableLabel(obj->cal[i]);  
        obj->cal[i] = NULL;  
        }  
      } /* end of for i over list length */ 


    obj->cal_len = 0;    
    return i;    
}    


/* Function:  swap_ev_GenericMatrix(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_ev_GenericMatrix
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [ExternVariable   **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_ev_GenericMatrix(ExternVariable   ** list,int i,int j)  
{
    ExternVariable   * temp; 
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_ev_GenericMatrix(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_ev_GenericMatrix which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [ExternVariable   **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_ev_GenericMatrix(ExternVariable   ** list,int left,int right,int (*comp)(ExternVariable   * ,ExternVariable   * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_ev_GenericMatrix(list,left,(left+right)/2); 
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_ev_GenericMatrix (list,++last,i);   
      }  
    swap_ev_GenericMatrix (list,left,last);  
    qsort_ev_GenericMatrix(list,left,last-1,comp);   
    qsort_ev_GenericMatrix(list,last+1,right,comp);  
}    


/* Function:  sort_ev_GenericMatrix(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_ev_GenericMatrix
 *
 *
 * Arg:         obj [UNKN ] Object containing list [GenericMatrix *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_ev_GenericMatrix(GenericMatrix * obj,int (*comp)(ExternVariable   *, ExternVariable   *)) 
{
    qsort_ev_GenericMatrix(obj->ev,0,obj->ev_len-1,comp);    
    return;  
}    


/* Function:  expand_ev_GenericMatrix(obj,len)
 *
 * Descrip:    Really an internal function for add_ev_GenericMatrix
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [GenericMatrix *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_ev_GenericMatrix(GenericMatrix * obj,int len) 
{


    if( obj->ev_maxlen > obj->ev_len )   {  
      warn("expand_GenericMatrixev_ called with no need");   
      return TRUE;   
      }  


    if( (obj->ev = (ExternVariable   ** ) ckrealloc (obj->ev,sizeof(ExternVariable   *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_GenericMatrix, returning FALSE");    
      return FALSE;  
      }  
    obj->ev_maxlen = len;    
    return TRUE; 
}    


/* Function:  add_ev_GenericMatrix(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [GenericMatrix *]
 * Arg:        add [OWNER] Object to add to the list [ExternVariable   *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_ev_GenericMatrix(GenericMatrix * obj,ExternVariable   * add) 
{
    if( obj->ev_len >= obj->ev_maxlen)   {  
      if( expand_ev_GenericMatrix(obj,obj->ev_len + GenericMatrixLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->ev[obj->ev_len++]=add;  
    return TRUE; 
}    


/* Function:  flush_ev_GenericMatrix(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [GenericMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_ev_GenericMatrix(GenericMatrix * obj) 
{
    int i;   


    for(i=0;i<obj->ev_len;i++)   { /*for i over list length*/ 
      if( obj->ev[i] != NULL)    {  
        free_ExternVariable(obj->ev[i]); 
        obj->ev[i] = NULL;   
        }  
      } /* end of for i over list length */ 


    obj->ev_len = 0; 
    return i;    
}    


/* Function:  GenericMatrix_alloc_std(void)
 *
 * Descrip:    Equivalent to GenericMatrix_alloc_len(GenericMatrixLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GenericMatrix *]
 *
 */
GenericMatrix * GenericMatrix_alloc_std(void) 
{
    return GenericMatrix_alloc_len(GenericMatrixLISTLENGTH); 
}    


/* Function:  GenericMatrix_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [GenericMatrix *]
 *
 */
GenericMatrix * GenericMatrix_alloc_len(int len) 
{
    GenericMatrix * out;/* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = GenericMatrix_alloc()) == NULL)    
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->state = (CellState ** ) ckcalloc (len,sizeof(CellState *))) == NULL)    {  
      warn("Warning, ckcalloc failed in GenericMatrix_alloc_len");   
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    if((out->special = (CellState ** ) ckcalloc (len,sizeof(CellState *))) == NULL)  {  
      warn("Warning, ckcalloc failed in GenericMatrix_alloc_len");   
      return NULL;   
      }  
    out->spec_len = 0;   
    out->spec_maxlen = len;  


    if((out->resource = (StructElement  ** ) ckcalloc (len,sizeof(StructElement  *))) == NULL)   {  
      warn("Warning, ckcalloc failed in GenericMatrix_alloc_len");   
      return NULL;   
      }  
    out->res_len = 0;    
    out->res_maxlen = len;   


    if((out->cal = (CollapsableLabel ** ) ckcalloc (len,sizeof(CollapsableLabel *))) == NULL)    {  
      warn("Warning, ckcalloc failed in GenericMatrix_alloc_len");   
      return NULL;   
      }  
    out->cal_len = 0;    
    out->cal_maxlen = len;   


    if((out->ev = (ExternVariable   ** ) ckcalloc (len,sizeof(ExternVariable   *))) == NULL) {  
      warn("Warning, ckcalloc failed in GenericMatrix_alloc_len");   
      return NULL;   
      }  
    out->ev_len = 0; 
    out->ev_maxlen = len;    


    return out;  
}    


/* Function:  hard_link_GenericMatrix(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GenericMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [GenericMatrix *]
 *
 */
GenericMatrix * hard_link_GenericMatrix(GenericMatrix * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a GenericMatrix object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  GenericMatrix_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GenericMatrix *]
 *
 */
GenericMatrix * GenericMatrix_alloc(void) 
{
    GenericMatrix * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(GenericMatrix *) ckalloc (sizeof(GenericMatrix))) == NULL)  {  
      warn("GenericMatrix_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->name = NULL;    
    out->type = 0;   
    out->state = NULL;   
    out->len = out->maxlen = 0;  
    out->special = NULL; 
    out->spec_len = out->spec_maxlen = 0;    
    out->query = NULL;   
    out->query_name = NULL;  
    out->query_len = NULL;   
    out->target = NULL;  
    out->target_name = NULL; 
    out->target_len = NULL;  
    out->resource = NULL;    
    out->res_len = out->res_maxlen = 0;  
    out->cal = NULL; 
    out->cal_len = out->cal_maxlen = 0;  
    out->ev = NULL;  
    out->ev_len = out->ev_maxlen = 0;    
    out->defscore_all_states = NULL; 
    out->window_i = 0;   
    out->window_j = 0;   
    out->footprint = 1;  
    out->cansearch = FALSE;  
    out->canlabel = FALSE;   
    out->specialtospecial = FALSE;   
    out->sh = NULL;  
    out->calcfunc = FALSE;   
    out->sc = FALSE; 
    out->mts = NULL; 


    return out;  
}    


/* Function:  free_GenericMatrix(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GenericMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [GenericMatrix *]
 *
 */
GenericMatrix * free_GenericMatrix(GenericMatrix * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a GenericMatrix obj. Should be trappable"); 
      return NULL;   
      }  


#ifdef PTHREAD   
    assert(pthread_mutex_lock(&(obj->dynamite_mutex)) == 0); 
#endif   
    if( obj->dynamite_hard_link > 1)     {  
      return_early = 1;  
      obj->dynamite_hard_link--; 
      }  
#ifdef PTHREAD   
    assert(pthread_mutex_unlock(&(obj->dynamite_mutex)) == 0);   
#endif   
    if( return_early == 1)   
      return NULL;   
    if( obj->name != NULL)   
      ckfree(obj->name);     
    if( obj->state != NULL)  {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->state[i] != NULL)   
          free_CellState(obj->state[i]); 
        }  
      ckfree(obj->state);    
      }  
    if( obj->special != NULL)    {  
      for(i=0;i<obj->spec_len;i++)   {  
        if( obj->special[i] != NULL) 
          free_CellState(obj->special[i]);   
        }  
      ckfree(obj->special);  
      }  
    if( obj->query != NULL)  
      free_StructElement(obj->query);    
    if( obj->query_name != NULL) 
      ckfree(obj->query_name);   
    if( obj->query_len != NULL)  
      ckfree(obj->query_len);    
    /* obj->qtype is linked in */ 
    if( obj->target != NULL) 
      free_StructElement(obj->target);   
    if( obj->target_name != NULL)    
      ckfree(obj->target_name);  
    if( obj->target_len != NULL) 
      ckfree(obj->target_len);   
    /* obj->ttype is linked in */ 
    if( obj->resource != NULL)   {  
      for(i=0;i<obj->res_len;i++)    {  
        if( obj->resource[i] != NULL)    
          free_StructElement(obj->resource[i]);  
        }  
      ckfree(obj->resource); 
      }  
    if( obj->cal != NULL)    {  
      for(i=0;i<obj->cal_len;i++)    {  
        if( obj->cal[i] != NULL) 
          free_CollapsableLabel(obj->cal[i]);    
        }  
      ckfree(obj->cal);  
      }  
    if( obj->ev != NULL) {  
      for(i=0;i<obj->ev_len;i++) {  
        if( obj->ev[i] != NULL)  
          free_ExternVariable(obj->ev[i]);   
        }  
      ckfree(obj->ev);   
      }  
    if( obj->defscore_all_states != NULL)    
      ckfree(obj->defscore_all_states);  
    if( obj->sh != NULL) 
      free_StructHolder(obj->sh);    
    if( obj->calcfunc != NULL)   
      ckfree(obj->calcfunc);     
    if( obj->sc != NULL) 
      free_Scope(obj->sc);   
    if( obj->mts != NULL)    
      free_MethodTypeSet(obj->mts);  


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif
