#include "compugen.h"
#include "dynafunc.h" /* only for get_argstr and get_chainstr */

#ifndef MAX
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#endif

static void            update_state_names(GenericMatrix * gm, OneModel * om,
					  char ** start_state_name, char ** end_state_name);

static char          * copy_str(char * in_str);

static compugen_rc_t   update_transitions(GenericMatrix * gm, OneModel * om,
					  char * start_state_name, char * end_state_name);

static OneModelTrans * update_transition_fields(CellSource * from_state, CellState * to_state,
						int transition_num, char * start_state_name,
						char * end_state_name, int to_semistate);

static OneModelTrans * update_transition_from_semistate(CellState * state, int transition_num);

static compugen_rc_t   parse_calc_line(ExprTree * expr_tree, OneModel * om, int transition_num,
				       seq_func_table_t seq_func_table,
				       char * query_name, char * target_name);

static compugen_rc_t   update_calc_om_data(ExprTree * tree, OneModel * om, int transition_num,
					   seq_func_table_t seq_func_table,
					   char * query_name, char * target_name);

static compugen_rc_t   update_tseq(ExprTree * tree, OneModel * om, int transition_num,
				   seq_func_table_t seq_func_table);

static compugen_rc_t   update_tprf(ExprTree * tree, OneModel * om, int transition_num);

static compugen_rc_t   update_tprf(ExprTree * tree, OneModel * om, int transition_num);

static compugen_rc_t   update_w(ExprTree * tree, OneModel * om, int transition_num,
				seq_func_table_t seq_func_table, char * target_name);

static compugen_rc_t   update_wseq(ExprTree * tree, OneModel * om, int transition_num,
				   seq_func_table_t seq_func_table);

static ExprTree      * find_wseq(ExprTree * tree, char * target_name);

static int             get_seq_diff(ExprTree * tree);

static int             is_word(ExprTree * tree, char * str);

static int             is_diff(ExprTree * tree);

static char          * get_seq_func_name(ExprTree * tree);

static int             get_seq_func_index(char * dy_func_name, seq_func_table_t seq_func_table);

static int             find_trees_num(ExprTree * tree);

static int             divide_expr_tree(ExprTree * tree, ExprTree * trees_array[], int trees_num);

static calc_type_t     find_calc_type(ExprTree * tree, char * query_name, char * target_name);

static void            find_prf_seq(ExprTree * tree, int * is_prf, int * is_seq,
				    char * query_name, char * target_name);

static compugen_rc_t   add_seqfunc(char * new_seqfunc, char * seq_funcs[]);

static OneModelTrans * OneModelTrans_alloc(void);

static OneModel      * OneModel_alloc(void);

static void            make_seq_func_table ( seq_func_table_t * table );

static char          * find_datatypes(char * target_type);

static char          * create_semistate_name ( char * state_name );


static int tpos;


/***************************************************************************
 *
 * FUNCTION    : OneModel_from_GenericMatrix
 *
 * DESCRIPTION : Filling One Model structure according to generic matrix.
 *               
 * INPUT       : GenericMatrix * gm
 *
 * OUTPUT      : None.
 *
 * RETURN      : OneModel * - The One Model structure.
 *
 ***************************************************************************/

OneModel * OneModel_from_GenericMatrix(GenericMatrix * gm)
{
  OneModel * out;
  char     * start_state_name = NULL;
  char     * end_state_name = NULL;

  out = OneModel_alloc();

  /* Updating model name */
  out->name = copy_str ( gm->name );

  /* Updating data typs */
  if ( (out->datatype_str = find_datatypes(gm->ttype->logical)) == NULL )
    return NULL;

  /* Updating states and midstates names */
  update_state_names ( gm, out, &start_state_name, &end_state_name );

  /* Updating states number */
  out->states_num = gm->len + OM_TMP_STATES_NUM;

  tpos = 4 * (out->states_num + out->semistates_num);

  out->arg_str = get_argstr_GenericMatrix(gm);
  out->chain_str = get_chainstr_GenericMatrix(gm);
  out->sc = gm->sc;

  /* Updating transition table */
  if ( update_transitions ( gm, out, start_state_name, end_state_name ) != COMPUGEN_OK ) 
    return NULL;

  out->prof_line_len = tpos;

  if ( start_state_name != NULL )
    free ( start_state_name );
  if ( end_state_name != NULL )
    free ( end_state_name );

  return out;
}

/***************************************************************************
 *
 * FUNCTION    : update_state_names
 *
 * DESCRIPTION : Updating states semistates and midstates names in One Model
 *               structure.
 *               
 * INPUT       : GenericMatrix * gm
 *
 * OUTPUT      : OneModel * om
 *               char ** start_state_name
 *               char ** end_state_name
 *
 * RETURN      : None.
 *
 ***************************************************************************/

static void update_state_names(GenericMatrix * gm, OneModel * om,
			       char ** start_state_name, char ** end_state_name)
{
  int   midstates_num = 0;
  int   i;
 
  /* Updating TMPSTART and TMPEND state names */
  om->statenames[OM_TMPSTART] = copy_str ( "TMPSTART" );
  om->statenames[OM_TMPEND] = copy_str ( "TMPEND" );
 
  /* Updating regular state names */
  for ( i = 0; i < gm->len; i++ ) {
    om->statenames[OM_TMP_STATES_NUM+i] = copy_str ( gm->state[i]->name );
    /* Updating semistate name if exists */
    if ( gm->state[i]->etr != NULL ) {
      om->semistatenames[om->semistates_num] = create_semistate_name ( gm->state[i]->name );
      (om->semistates_num)++;
    }
  }

  /* Updating midstate names */
  for ( i = 0; i < gm->spec_len; i++ ) {
    if ( gm->special[i]->is_start )
      *start_state_name = copy_str ( gm->special[i]->name );
    else if ( gm->special[i]->is_end )
      *end_state_name = copy_str ( gm->special[i]->name );
    else {
      om->midstatenames[midstates_num] = copy_str ( gm->special[i]->name );
      midstates_num++;
      /* Updating semistate name if exists */
      if ( gm->special[i]->etr != NULL ) {
	om->semistatenames[om->semistates_num] = create_semistate_name ( gm->special[i]->name );
	(om->semistates_num)++;
      }
    }
  }
}

/***************************************************************************
 *
 * FUNCTION    : copy_str
 *
 * DESCRIPTION : Allocating a string with the length of the input string, 
 *               and copying the input string to the new string.
 *               
 * INPUT       : char * in_str
 *
 * OUTPUT      : None.
 *
 * RETURN      : char * - The new string.
 *
 ***************************************************************************/

static char * copy_str(char * in_str)
{
  char * out_str;

  if ( (out_str = (char *)malloc(strlen(in_str) + 1)) == NULL ) {
    perror ( "copy_str : malloc" );
    exit ( 1 );
  }
  strcpy ( out_str, in_str );

  return out_str;
}

/***************************************************************************
 *
 * FUNCTION    : update_transitions
 *
 * DESCRIPTION : Updating One Model transition structure for all the
 *               transitions.
 *               
 * INPUT       : GenericMatrix * gm
 *               char * start_state_name
 *               char * end_state_name
 *
 * OUTPUT      : OneModel * om
 *
 * RETURN      : compugen_rc_t
 *
 ***************************************************************************/

static compugen_rc_t update_transitions(GenericMatrix * gm, OneModel * om,
					char * start_state_name, char * end_state_name)
{
  int                i, j;
  seq_func_table_t   seq_func_table;
  compugen_rc_t      rc;
  int                transition_num = OM_TMP_STATES_NUM;
  char               temp_name[20];
  int                to_semistate;

  make_seq_func_table ( &seq_func_table );

  /* Updating transition START -> TMPSTART and TMPEND -> END */

  om->trans[OM_TMPSTART] = OneModelTrans_alloc();
  om->trans[OM_TMPSTART]->from_state = copy_str ( "START" );
  om->trans[OM_TMPSTART]->to_state = copy_str ( "TMPSTART" );
  om->trans[OM_TMPSTART]->name = copy_str ( "trans0" );
  om->trans[OM_TMPSTART]->xlabel = copy_str ( "-" );
  om->trans[OM_TMPSTART]->ylabel = copy_str ( "-" );

  om->trans[OM_TMPEND] = OneModelTrans_alloc();
  om->trans[OM_TMPEND]->from_state = copy_str ( "TMPEND" );
  om->trans[OM_TMPEND]->to_state = copy_str ( "END" );
  om->trans[OM_TMPEND]->name = copy_str ( "trans1" );
  om->trans[OM_TMPEND]->xlabel = copy_str ( "-" );
  om->trans[OM_TMPEND]->ylabel = copy_str ( "-" );

  /* Updating all transitions */

  /* Transitions to states */

  for ( i = 0; i < gm->len; i++ ) {

    /* Checking if there is a semistate - if it is adding a transition semistate -> state */
    if ( gm->state[i]->etr != NULL ) {
      om->trans[transition_num] = update_transition_from_semistate(gm->state[i], transition_num);
      if ( (rc = parse_calc_line(gm->state[i]->etr, om, transition_num, seq_func_table,
				 gm->query->name, gm->target->name)) != COMPUGEN_OK )
	return rc;
      transition_num++;
      to_semistate = TRUE;
    }
    else
      to_semistate = FALSE;

    /* Adding all transitions to state */
    for ( j = 0; j < gm->state[i]->len; j++ ) {
      push_errormsg_stack("Reconciling Compugen functions for %s to %s",
			  gm->state[i]->source[j]->state_source,gm->state[i]->name);
      om->trans[transition_num] = update_transition_fields ( gm->state[i]->source[j],
							     gm->state[i], transition_num,
							     start_state_name, end_state_name,
							     to_semistate );
      if ( (rc = parse_calc_line(gm->state[i]->source[j]->etr, om, transition_num,
				 seq_func_table, gm->query->name, gm->target->name))
	   != COMPUGEN_OK )
	return rc;
      pop_errormsg_stack();
      transition_num++;
    }
  }

  /* Transitions to midstates */

  for ( i = 0; i < gm->spec_len; i++ ) {

    /* Checking if there is a semistate - if it is adding a transition semistate -> state */
    if ( gm->special[i]->etr != NULL ) {
      om->trans[transition_num] = update_transition_from_semistate(gm->special[i], transition_num);
      if ( (rc = parse_calc_line(gm->special[i]->etr, om, transition_num, seq_func_table,
				 gm->query->name, gm->target->name)) != COMPUGEN_OK )
	return rc;
      transition_num++;
      to_semistate = TRUE;
    }
    else
      to_semistate = FALSE;

    /* Adding all transitions to state */
    for ( j = 0; j < gm->special[i]->len; j++ ) {
      push_errormsg_stack("Reconciling Compugen functions for %s to %s",
			  gm->special[i]->source[j]->state_source,gm->special[i]->name);
      om->trans[transition_num] = update_transition_fields ( gm->special[i]->source[j],
							     gm->special[i], transition_num,
							     start_state_name, end_state_name,
							     to_semistate );
      if ( (rc = parse_calc_line(gm->special[i]->source[j]->etr, om, transition_num,
				 seq_func_table, gm->query->name, gm->target->name))
	   != COMPUGEN_OK )
	return rc;
      pop_errormsg_stack();
      transition_num++;
    }
  }

  /* A transition TMPSTART -> TMPSTART - needed for One Model */
  om->trans[transition_num] = OneModelTrans_alloc();
  om->trans[transition_num]->from_state = copy_str ( "TMPSTART" );
  om->trans[transition_num]->to_state = copy_str ( "TMPSTART" );
  sprintf ( temp_name, "trans%d", transition_num );
  om->trans[transition_num]->name = copy_str ( temp_name );
  om->trans[transition_num]->xlabel = copy_str ( "-" );
  om->trans[transition_num]->ylabel = copy_str ( "-" );
  om->trans[transition_num]->dx = om->trans[transition_num]->dy = 1;
  om->trans[transition_num]->tpos_tprf = tpos;
  tpos++;

  return COMPUGEN_OK;
}

/***************************************************************************
 *
 * FUNCTION    : update_transition_fields
 *
 * DESCRIPTION : Allocating a One Model transition structure and filling in
 *               all the fields that don't depend on the calc expression.
 *               
 * INPUT       : CellSource * from_state
 *               CellState * to_state
 *               int transition_num
 *               char * start_state_name
 *               char * end_state_name
 *               int to_semistate
 *
 * OUTPUT      : None.
 *
 * RETURN      : OneModelTrans *
 *
 ***************************************************************************/

static OneModelTrans * update_transition_fields(CellSource * from_state, CellState * to_state,
						int transition_num, char * start_state_name,
						char * end_state_name, int to_semistate)
{
  char            temp_name[20];
  OneModelTrans * out;

  out = OneModelTrans_alloc();

  if ( (start_state_name != NULL) && (strcmp(from_state->state_source, start_state_name) == 0) )
    out->from_state = copy_str ( "TMPSTART" );
  else
    out->from_state = copy_str ( from_state->state_source );

  if ( to_semistate )
    out->to_state = create_semistate_name ( to_state->name );
  else if ( (end_state_name != NULL) && (strcmp(to_state->name, end_state_name) == 0) )
    out->to_state = copy_str ( "TMPEND" );
  else
    out->to_state = copy_str ( to_state->name );

  sprintf ( temp_name, "trans%d", transition_num );
  out->name = copy_str ( temp_name );

  out->dx = (from_state->offj < 0) ? MAX(to_state->offj, 0) : from_state->offj;
  out->dy = (from_state->offi < 0) ? MAX(to_state->offi, 0) : from_state->offi;

  out->xlabel = copy_str ( "sequence" );
  out->ylabel = copy_str ( "sequence" );

  return out;
}

/***************************************************************************
 *
 * FUNCTION    : update_transition_from_semistate
 *
 * DESCRIPTION : Allocating a One Model transition structure for a
 *               transition from a semistate to a state, and filling in all
 *               the fields that don't depend on the calc expression.
 *               
 * INPUT       : CellState * state
 *               int transition_num
 *
 * OUTPUT      : None.
 *
 * RETURN      : OneModelTrans *
 *
 ***************************************************************************/

static OneModelTrans * update_transition_from_semistate(CellState * state, int transition_num)
{
  char            temp_name[20];
  OneModelTrans * out;

  out = OneModelTrans_alloc();

  out->from_state = create_semistate_name ( state->name );
  out->to_state = copy_str ( state->name );

  sprintf ( temp_name, "trans%d", transition_num );
  out->name = copy_str ( temp_name );

  out->dx = 0;
  out->dy = 0;

  out->xlabel = copy_str ( "-" );
  out->ylabel = copy_str ( "-" );

  return out;
}

/***************************************************************************
 *
 * FUNCTION    : parse_calc_line
 *
 * DESCRIPTION : 
 *               
 * INPUT       : ExprTree * expr_tree
 *               int transition_num
 *               seq_func_table_t seq_func_table
 *               char * query_name
 *               char * target_name
 *
 * OUTPUT      : OneModel * om
 *
 * RETURN      : compugen_rc_t
 *
 ***************************************************************************/

static compugen_rc_t parse_calc_line(ExprTree * expr_tree, OneModel * om, int transition_num,
				     seq_func_table_t seq_func_table,
				     char * query_name, char * target_name)
{
  int             i;
  ExprTree     ** trees;
  int             trees_num;
  compugen_rc_t   rc;

  /* If root is not a statement - error */
  if ( (expr_tree->type != ETR_STATEMENT) || (expr_tree->nochild != 1) ) {
    fprintf ( stderr, "parse_calc_line : Expression tree root has illegal type\n" );
    exit ( 1 );
  }

  /* Finding trees number and allocating trees */

  trees_num = find_trees_num(expr_tree->child[0]);

  if ( trees_num == 0 )
    return COMPUGEN_OK;

  if ( trees_num > MAX_CALC_TYPES ) { /* Until we handle this */
    warn ( "Too many trees in calc expression" );
    return COMPUGEN_TOO_MANY_CALC_FUNCS;
  }
  
  if ( (trees = (ExprTree **)calloc(trees_num, sizeof(ExprTree *))) == NULL ) {
    perror ( "parse_calc_line : calloc" );
    exit ( 1 );
  }

  /* Dividing expression tree into sub-trees */
  divide_expr_tree ( expr_tree->child[0], trees, 0 );

  /* Updating OM data */
  for ( i = 0; i < trees_num; i++ ) {
    if ( (rc = update_calc_om_data(trees[i], om, transition_num, seq_func_table,
				   query_name, target_name)) != COMPUGEN_OK )
      return rc;
  }

  free ( trees );

  return COMPUGEN_OK;
}

/***************************************************************************
 *
 * FUNCTION    : update_calc_om_data
 *
 * DESCRIPTION : 
 *               
 * INPUT       : ExprTree * tree
 *               OneModel * om
 *               int transition_num
 *               seq_func_table_t seq_func_table
 *
 * OUTPUT      : None.
 *
 * RETURN      : compugen_rc_t - OK / error code.
 *
 ***************************************************************************/

static compugen_rc_t update_calc_om_data(ExprTree * tree, OneModel * om, int transition_num,
					 seq_func_table_t seq_func_table,
					 char * query_name, char * target_name)
{
  switch ( find_calc_type(tree, query_name, target_name) ) {
  case CALC_SEQ:
    return update_tseq(tree, om, transition_num, seq_func_table);
    break;
  case CALC_PRF:
    return update_tprf(tree, om, transition_num);
    break;
  case CALC_W:
    return update_w(tree, om, transition_num, seq_func_table, target_name);
    break;
  }
}

/***************************************************************************
 *
 * FUNCTION    : update_tseq
 *
 * DESCRIPTION : 
 *               
 * INPUT       : ExprTree * tree
 *               OneModel * om
 *               int transition_num
 *               seq_func_table_t seq_func_table
 *
 * OUTPUT      : None.
 *
 * RETURN      : compugen_rc_t - OK / error code.
 *
 ***************************************************************************/

static compugen_rc_t update_tseq(ExprTree * tree, OneModel * om, int transition_num,
				 seq_func_table_t seq_func_table)
{
  int    index;
  char * dy_tseq_name;
  
  if ( om->trans[transition_num]->tseq != NULL ) {
    warn ( "More than one tseq" );
    return COMPUGEN_TOO_MANY_CALC_FUNCS;
  }

  if ( (dy_tseq_name = get_seq_func_name(tree)) == NULL ) {
    warn ( "tseq name was not found in tree" );
    return COMPUGEN_ILLEGAL_SEQ_EXPR;
  }

  if ( (index = get_seq_func_index(dy_tseq_name, seq_func_table)) == -1 ) {
    warn ( "tseq name was not found in tseq table %s", dy_tseq_name );
    return COMPUGEN_ILLEGAL_SEQ_EXPR;
  }

  om->trans[transition_num]->tseq = copy_str ( seq_func_table.func[index].om_func_name );

  om->trans[transition_num]->dseq = get_seq_diff ( tree );

  add_seqfunc ( om->trans[transition_num]->tseq, om->seqfuncs );

  return COMPUGEN_OK;
}

/***************************************************************************
 *
 * FUNCTION    : update_tprf
 *
 * DESCRIPTION : 
 *               
 * INPUT       : ExprTree * tree
 *               OneModel * om
 *               int transition_num
 *
 * OUTPUT      : None.
 *
 * RETURN      : compugen_rc_t - OK / error code.
 *
 ***************************************************************************/

static compugen_rc_t update_tprf(ExprTree * tree, OneModel * om, int transition_num)
{
  if ( om->trans[transition_num]->tprf != NULL ) {
    warn ( "more than one tprf" );
    return COMPUGEN_TOO_MANY_CALC_FUNCS;
  }

  om->trans[transition_num]->tprf = tree;
  om->trans[transition_num]->tpos_tprf = tpos;
  om->trans[transition_num]->range = 1;
  tpos++;

  return COMPUGEN_OK;
}

/***************************************************************************
 *
 * FUNCTION    : update_tprf
 *
 * DESCRIPTION : 
 *               
 * INPUT       : ExprTree * tree
 *               OneModel * om
 *               int transition_num
 *               seq_func_table_t seq_func_table
 *
 * OUTPUT      : None.
 *
 * RETURN      : compugen_rc_t - OK / error code.
 *
 ***************************************************************************/

static compugen_rc_t update_w(ExprTree * tree, OneModel * om, int transition_num,
			      seq_func_table_t seq_func_table, char * target_name)
{
  compugen_rc_t   rc;
  ExprTree      * wseq_tree;

  if ( om->trans[transition_num]->wprf != NULL ) {
    warn ( "more than one wprf" );
    return COMPUGEN_TOO_MANY_CALC_FUNCS;
  }

  om->trans[transition_num]->wprf = tree;
  om->trans[transition_num]->tpos_wprf = tpos;

  if ( (wseq_tree = find_wseq(tree, target_name)) == NULL ) {
    warn ( "No wseq was found" );
    return COMPUGEN_NO_WSEQ;
  }

  if ( (rc = update_wseq(wseq_tree, om, transition_num, seq_func_table)) != COMPUGEN_OK )
    return rc;

  tpos += om->trans[transition_num]->range;

  return COMPUGEN_OK;
}

/***************************************************************************
 *
 * FUNCTION    : update_wseq
 *
 * DESCRIPTION : 
 *               
 * INPUT       : ExprTree * tree
 *               OneModel * om
 *               int transition_num
 *               seq_func_table_t seq_func_table
 *
 * OUTPUT      : None.
 *
 * RETURN      : compugen_rc_t - OK / error code.
 *
 ***************************************************************************/

static compugen_rc_t update_wseq(ExprTree * tree, OneModel * om, int transition_num,
				 seq_func_table_t seq_func_table)
{
  int    index;
  char * dy_wseq_name;
  
  if ( om->trans[transition_num]->wseq != NULL ) {
    warn ( "More than one wseq" );
    return COMPUGEN_TOO_MANY_CALC_FUNCS;
  }

  if ( (dy_wseq_name = get_seq_func_name(tree)) == NULL ) {
    warn ( "wseq name was not found in tree");
    return COMPUGEN_ILLEGAL_SEQ_EXPR;
  }

  if ( (index = get_seq_func_index(dy_wseq_name, seq_func_table)) == -1 ) {
    warn ( "wseq name was not found in table %s", dy_wseq_name );
    return COMPUGEN_ILLEGAL_SEQ_EXPR;
  }

  om->trans[transition_num]->wseq = copy_str ( seq_func_table.func[index].om_func_name );
  om->trans[transition_num]->wseq_node_to_replace = tree;
  om->trans[transition_num]->dwx = get_seq_diff ( tree );
  om->trans[transition_num]->range = seq_func_table.func[index].range;
  if( seq_func_table.func[index].map_func != NULL ) 
    om->trans[transition_num]->map_func = copy_str(seq_func_table.func[index].map_func);

  add_seqfunc ( om->trans[transition_num]->wseq, om->seqfuncs );

  return COMPUGEN_OK;
}

/***************************************************************************
 *
 * FUNCTION    : find_wseq
 *
 * DESCRIPTION : 
 *               
 * INPUT       : ExprTree * tree
 *
 * OUTPUT      : None.
 *
 * RETURN      : ExprTree * - wseq expression tree.
 *
 ***************************************************************************/

static ExprTree * find_wseq ( ExprTree * tree, char * target_name )
{
  ExprTree * wseq_tree;
  int        i;

  if ( (tree->type == ETR_METHOD) &&
       (tree->nochild == 2) &&
       (tree->child[1]->type == ETR_COMMALIST) ) {
    for ( i = 0; i < tree->child[1]->nochild; i++ ) {
      if ( is_word(tree->child[1]->child[i], target_name) )
	return tree;
    }
  }

  for ( i = 0; i < tree->nochild; i++ ) {
    if ( (wseq_tree = find_wseq(tree->child[i], target_name)) != NULL )
      return wseq_tree;
  }

  return NULL;
}

/***************************************************************************
 *
 * FUNCTION    : get_seq_diff
 *
 * DESCRIPTION : 
 *               
 * INPUT       : ExprTree * tree
 *
 * OUTPUT      : None.
 *
 * RETURN      : int - Sequence index difference.
 *
 ***************************************************************************/

static int get_seq_diff ( ExprTree * tree )
{
  int   i;

  if ( (tree->type != ETR_METHOD) || (tree->nochild < 2) ||
       (tree->child[1]->type != ETR_COMMALIST) ) {
    warn ( "Illegal tseq/wseq expression" );
    return 0;
  }

  for ( i = 0; i < tree->child[1]->nochild; i++ ) {
    if ( is_word(tree->child[1]->child[i], "j") )
      return 0;
    else if ( is_diff(tree->child[1]->child[i]) &&
	      is_word(tree->child[1]->child[i]->child[0], "j") &&
	      (tree->child[1]->child[i]->child[2]->type == ETR_NUMBER) )
      return (atoi(tree->child[1]->child[i]->child[2]->word));      
  }

  return 0;
}

/***************************************************************************
 *
 * FUNCTION    : is_word
 *
 * DESCRIPTION : 
 *               
 * INPUT       : ExprTree * tree
 *               char * str
 *
 * OUTPUT      : None.
 *
 * RETURN      : int - TRUE / FALSE.
 *
 ***************************************************************************/

static int is_word ( ExprTree * tree, char * str )
{
  return ( (tree->type == ETR_TAG) &&
	   (tree->nochild == 1) &&
	   (tree->child[0]->type == ETR_NAME) &&
	   (strcmp(tree->child[0]->word, str) == 0) );
}

/***************************************************************************
 *
 * FUNCTION    : is_diff
 *
 * DESCRIPTION : 
 *               
 * INPUT       : ExprTree * tree
 *
 * OUTPUT      : None.
 *
 * RETURN      : int - TRUE / FALSE.
 *
 ***************************************************************************/

static int is_diff ( ExprTree * tree )
{
  return ( (tree->type == ETR_EXPRESSION) &&
	   (tree->nochild == 3) &&
	   (tree->child[1]->type == ETR_OPERATOR) &&
	   (strcmp(tree->child[1]->word, "-") == 0) );
}

/***************************************************************************
 *
 * FUNCTION    : get_seq_func_name
 *
 * DESCRIPTION : 
 *               
 * INPUT       : ExprTree * tree
 *
 * OUTPUT      : None.
 *
 * RETURN      : char * - Sequence function name.
 *
 ***************************************************************************/

static char * get_seq_func_name ( ExprTree * tree )
{
  if ( (tree->type != ETR_METHOD) || (tree->nochild < 1) || (tree->child[0]->type != ETR_TAG) ||
       (tree->child[0]->nochild != 1) || (tree->child[0]->child[0]->type != ETR_NAME) )
    return NULL;

  return tree->child[0]->child[0]->word;
}

/***************************************************************************
 *
 * FUNCTION    : get_seq_func_index
 *
 * DESCRIPTION : 
 *               
 * INPUT       : char * dy_func_name
 *               seq_func_table_t seq_func_table
 *
 * OUTPUT      : None.
 *
 * RETURN      : int - The index of the sequence function in the table.
 *
 ***************************************************************************/

static int get_seq_func_index ( char * dy_func_name, seq_func_table_t seq_func_table )
{
  int   i;

  for ( i = 0; i < seq_func_table.funcs_num; i++ ) {
    if ( strcmp(seq_func_table.func[i].dy_func_name, dy_func_name) == 0 )
      return i;
  }
  return -1;
}


/***************************************************************************
 *
 * FUNCTION    : find_trees_num
 *
 * DESCRIPTION : Finding the number of sub-trees separated by '+'.
 *               
 * INPUT       : ExprTree * tree
 *
 * OUTPUT      : None.
 *
 * RETURN      : int - The number of trees.
 *
 ***************************************************************************/

static int find_trees_num ( ExprTree * tree )
{

  if ( (tree == NULL) || ((tree->type == ETR_NUMBER) && (strcmp(tree->word, "0") == 0)) )
    return 0;

  if ( is_sum_tree(tree) )
    return ( find_trees_num(tree->child[0]) + find_trees_num(tree->child[2]) );

  return 1;
}

/***************************************************************************
 *
 * FUNCTION    : divide_expr_tree
 *
 * DESCRIPTION : Dividing expression tree to sub-trees of the parts of the
 *               expression separated by '+'.
 *               
 * INPUT       : ExprTree * tree
 *               int trees_num
 *
 * OUTPUT      : ExprTree * trees_array[]
 *
 * RETURN      : int - New number of trees.
 *
 ***************************************************************************/

static int divide_expr_tree ( ExprTree * tree, ExprTree * trees_array[], int trees_num )
{
  if ( (tree == NULL) || ((tree->type == ETR_NUMBER) && (strcmp(tree->word, "0") == 0)) )
    return trees_num;
    
  if ( is_sum_tree(tree) ) {
    trees_num = divide_expr_tree ( tree->child[0], trees_array, trees_num );
    trees_num = divide_expr_tree ( tree->child[2], trees_array, trees_num );
  }
  else {
    trees_array[trees_num] = tree;
    trees_num++;
  }

  return trees_num;
}

/***************************************************************************
 *
 * FUNCTION    : find_calc_type
 *
 * DESCRIPTION : Checking if the expression depends on profile, sequence, or
 *               both.
 *               
 * INPUT       : ExprTree * tree
 *
 * OUTPUT      : None.
 *
 * RETURN      : calc_type_t - CALC_PRF / CALC_SEQ / CALC_W.
 *
 ***************************************************************************/

static calc_type_t find_calc_type(ExprTree * tree, char * query_name, char * target_name)
{
  int   is_prf = FALSE;
  int   is_seq = FALSE;

  find_prf_seq ( tree, &is_prf, &is_seq, query_name, target_name );

  if ( is_seq ) {
    if ( is_prf )
      return CALC_W;
    else
      return CALC_SEQ;
  }
  return CALC_PRF;
}

/***************************************************************************
 *
 * FUNCTION    : find_prf_seq
 *
 * DESCRIPTION : Looking for query_name, target_name, "i" and "j" in the
 *               expression tree.
 *               
 * INPUT       : ExprTree * tree
 *
 * OUTPUT      : int * is_prf
 *               int * is_seq
 *
 * RETURN      : None.
 *
 ***************************************************************************/

static void find_prf_seq(ExprTree * tree, int * is_prf, int * is_seq,
			 char * query_name, char * target_name)
{
  int i;

  if ( tree == NULL )
    return;

  if ( is_word(tree, query_name) || is_word(tree, "i") )
    *is_prf = TRUE;
  else if ( is_word(tree, target_name) || is_word(tree, "j") )
      *is_seq = TRUE;

  for ( i = 0; i < tree->nochild; i++ )
    find_prf_seq ( tree->child[i], is_prf, is_seq, query_name, target_name );
}

/***************************************************************************
 *
 * FUNCTION    : add_seqfunc
 *
 * DESCRIPTION : 
 *               
 * INPUT       : char * new_seqfunc
 *               char * seq_funcs[]
 *
 * OUTPUT      : None.
 *
 * RETURN      : compugen_rc_t
 *
 ***************************************************************************/

static compugen_rc_t add_seqfunc ( char * new_seqfunc, char * seq_funcs[] ) 
{
  int   i;

  for ( i = 0; i < CUGEN_MAX_FUNCS; i++ ) {
    if ( seq_funcs[i] == NULL ) {
      seq_funcs[i] = copy_str ( new_seqfunc );
      return COMPUGEN_OK;
    }
    else if ( strcmp(seq_funcs[i], new_seqfunc) == 0 )
      return COMPUGEN_OK;
  }

  warn ( "More than %d seqfuncs", CUGEN_MAX_FUNCS );
  return COMPUGEN_TOO_MANY_CALC_FUNCS;
}


/***************************************************************************
 *
 * FUNCTION    : 
 *
 * DESCRIPTION : 
 *               
 * INPUT       : 
 *
 * OUTPUT      : 
 *
 * RETURN      : 
 *
 ***************************************************************************/

static OneModelTrans * OneModelTrans_alloc(void)
{
  OneModelTrans * out;

  out = (OneModelTrans *) malloc(sizeof(OneModelTrans));

  out->tprf = out->wprf = out->wseq_node_to_replace = NULL;

  out->tseq = out->wseq = out->from_state = out->to_state = out->name = out->map_func = NULL;

  out->tpos_tprf = out->tpos_wprf = out->range = out->dseq = out->dwx = out->dwy = out->dprf = 0;

  return out;
}


/***************************************************************************
 *
 * FUNCTION    : OneModel_alloc
 *
 * DESCRIPTION : 
 *               
 * INPUT       : None.
 *
 * OUTPUT      : None.
 *
 * RETURN      : OneModel *
 *
 ***************************************************************************/

static OneModel * OneModel_alloc(void)
{
  int i;
  OneModel * out;

  out= (OneModel * ) malloc(sizeof(OneModel));
  
  for(i=0;i<CUGEN_MAX_STATE_NAMES;i++) {
    out->statenames[i] = out->semistatenames[i] = out->midstatenames[i] = NULL;
  }

  for(i=0;i<CUGEN_MAX_FUNCS;i++) {
    out->seqfuncs[i] = NULL;
  }

  for(i=0;i<CUGEN_MAX_ONEMODEL_TRANS;i++) {
    out->trans[i] = NULL;
  }

  out->name = NULL;

  out->datatype_str = NULL;

  out->states_num = out->semistates_num = 0;

  return out;
}

/***************************************************************************
 *
 * FUNCTION    : make_seq_func_table
 *
 * DESCRIPTION : 
 *               
 * INPUT       : None.
 *
 * OUTPUT      : seq_func_table_t * table
 *
 * RETURN      : None.
 *
 ***************************************************************************/

static void make_seq_func_table ( seq_func_table_t * table )
{
  table->funcs_num = 3;

  table->func = (seq_func_t *)calloc(3, sizeof(seq_func_t));

  table->func[0].dy_func_name = copy_str ( "AMINOACID" );
  table->func[0].om_func_name = copy_str ( "ID_SEQ" );
  table->func[0].range = 32;
  table->func[0].map_func = copy_str ( "map_onemodel_aa" );

  table->func[1].dy_func_name = copy_str ( "CDNA_BASE" );
  table->func[1].om_func_name = copy_str ( "ID_SEQ" );
  table->func[1].range = 32;
  table->func[0].map_func = copy_str ( "map_onemodel_base" );

  table->func[2].dy_func_name = copy_str ( "CDNA_CODON" );
  table->func[2].om_func_name = copy_str ( "TRIPLETS_SEQ" );
  table->func[2].range = 125;
  table->func[0].map_func = NULL;
}

/***************************************************************************
 *
 * FUNCTION    : find_datatypes
 *
 * DESCRIPTION : 
 *               
 * INPUT       : char * target_type
 *
 * OUTPUT      : None.
 *
 * RETURN      : char * - One model data types.
 *
 ***************************************************************************/

static char * find_datatypes(char * target_type)
{
  if ( target_type == NULL )
    return NULL;

  if ( strcmp(target_type, "PROTEIN") == 0 )
    return "PP";

  if ( (strcmp(target_type, "DNA") == 0) ||
       (strcmp(target_type, "CDNA") == 0) ||
       (strcmp(target_type, "GENOMIC") == 0) )
    return "PN";
}

/***************************************************************************
 *
 * FUNCTION    : create_semistate_name
 *
 * DESCRIPTION : 
 *               
 * INPUT       : char * state_name
 *
 * OUTPUT      : None.
 *
 * RETURN      : char * - The semistate name.
 *
 ***************************************************************************/

static char * create_semistate_name ( char * state_name )
{
  char * semistate_name;

  if ( (semistate_name = (char *)malloc(strlen(state_name) + 5)) == NULL ) {
    perror ( "create_semistate_name : malloc" );
    exit ( 1 );
  }

  strcpy ( semistate_name, "SEMI" );
  strcat ( semistate_name, state_name );

  return semistate_name;
}
