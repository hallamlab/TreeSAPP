#ifdef _cplusplus
extern "C" {
#endif
#include "dy2om.h"


    /***
      Globals used to communicate with yacc parser

      These are defined in type.dy. Here we will define
      them as externs. 

      ***/

 extern char * calc_lex_string;
 extern int stringpos;
 extern ExprTree * root;


    /*
      Conversion functions from internal dynamite datastruct
      to OneModel structure
      */

# line 93 "dy2om.dy"
OneModel * OneModel_from_GenericMatrix(GenericMatrix * gm,MethodTypeSet * mts)
{
    OneModel * out;
    OmTrans * trans;
    Scope * sc;
    int i,j;

    out = OneModel_alloc_std();

    sc = std_Dynamite_Scope();

    add_GenericMatrix_Scope(sc,gm);

    /*** add states ***/

    for(i=0;i<gm->len;i++) {
	for(j=0;j<gm->state[i]->len;j++) {
	    trans = OmTrans_from_CellSource(0,gm->state[i]->source[j],mts,sc);
	    add_OneModel(out,trans); /* adds trans to out->trans list */
	}
    }

    /*** should now do specials ***/


    free_Scope(sc);

    return out;

}






# line 129 "dy2om.dy"
OmTrans * OmTrans_from_CellSource(int base_point,CellSource * s,MethodTypeSet * mts,Scope * sc)
{
    OmTrans * out;
    ExprTree * et;
    CugenYaccReturn ret;

    out = OmTrans_alloc();

    out->dx = s->offi;
    out->dy = s->offj;
    
    /** this is the tricky part **/

    /** copied from type.dy allocd_calc_line **/

    /*** set globals that we share with lex/yacc ***/

    calc_lex_string = s->source_expr;
    stringpos= 0;
    root = NULL;

    /*** parse it ***/
    
    yyparse();

    /*** check root, if NULL... otta here ***/

    if( root == NULL ) {
	fatal("Nasty error, we cannot parse this calc line in making OneModel transition");
    }
    
    find_toplevel_name(root);

    /** debugging stuff for now **/

    printf("Doing calc line [%s]\n",s->source_expr);

    print_ExprTree(root);

    /** ok, just enter the function correctly **/

    ret = descend_ExprTree_Cugen(root,mts,TRUE,sc,out);
    if( ret.is_doable == 0  ) {
      printf("Can't do this calc");
    }
    else { 
      printf("Can do this calc");
    }



    printf("End\n\n");



    return out;
}




# line 190 "dy2om.dy"
OmState * OmState_from_CellState(CellState * cs)
{
    OmState * out;

    out = OmState_alloc();

    /* warn if NULL...oops! */

    out->name = stringalloc(cs->name);
    if( cs->is_special_i == TRUE )
	out->type = OmState_TYPE_MIDSTATE;
    else out->type = OmState_TYPE_STATE;

    return out;
}

/* Descent of yacc grammar to get
   Compugen model 
   */


/* Function:  descend_ExprTree_Cugen(et,mts,is_attachable,sc,ot)
 *
 * Descrip:    Recursive parser of the yacc grammar for
 *             compugen.
 *
 *             If it is impossible for compugen to do, then
 *             is_doable == 0 and the descent should abort asap.
 *
 *             If not, at the lowest level of being able to map
 *             to a compugen number, it is mapped. 
 *
 *             targetname should be entire scope.
 *
 *
 * Arg:                   et [UNKN ] Undocumented argument [ExprTree *]
 * Arg:                  mts [UNKN ] Undocumented argument [MethodTypeSet *]
 * Arg:        is_attachable [UNKN ] Undocumented argument [boolean]
 * Arg:                   sc [UNKN ] Undocumented argument [Scope *]
 * Arg:                   ot [UNKN ] Undocumented argument [OmTrans *]
 *
 * Return [UNKN ]  Undocumented return value [CugenYaccReturn]
 *
 */
# line 223 "dy2om.dy"
CugenYaccReturn descend_ExprTree_Cugen(ExprTree * et,MethodTypeSet * mts,boolean is_attachable,Scope * sc,OmTrans * ot)
{
  int i;
  CugenYaccReturn out;
  CugenYaccReturn child;
  ScopeUnit * su;
  OmUnit * unit;


  set_CugenYaccReturn(&out);

  switch(et->type) {
  case ETR_STATEMENT : 
    /** descend into grammar **/

    for(i=0;i<et->nochild;i++) {
      child = descend_ExprTree_Cugen(et->child[i],mts,is_attachable,sc,ot);
      if( child.is_doable == 0 ) {
	out.is_doable = 0; 
	return out; /** get out of here now! **/
      }
    }

    /** nothing else, if we are here, then we are ok! **/
    out.is_doable = 1;

    return out;

  case ETR_EXPRESSION : 
    /** expression is *at the moment* only 3 things wide **/

    if( et->nochild != 3) {
      warn("Annoying. Sadly we have written the compugen parser to only have 3 piece expressions. You have %d piece ones. So... have to fail",et->nochild);
      out.is_doable = 0;
      return out;
    }

    /** we can only handle pluses **/

    if( et->child[1]->type != ETR_OPERATOR || et->child[1]->word[0] != '+' ) {
      warn("Sadly the compugen parser can only handle simple plus expressions. you have a middle expression of type [%d] with character [%s]",et->child[1]->type,CKS(et->child[1]->word));
      out.is_doable = 0;
      return out;
    }

    /** ok, now descend into 0 and 2 child **/

    child = descend_ExprTree_Cugen(et->child[0],mts,is_attachable,sc,ot);
    if( child.is_doable == 0 ) {
      out.is_doable = 0; 
      return out; /** get out of here now! **/
    }
    
      
	  
    
    child = descend_ExprTree_Cugen(et->child[2],mts,is_attachable,sc,ot);
    if( child.is_doable == 0 ) {
      out.is_doable = 0; 
      return out; /** get out of here now! **/
    }

    /** fine! **/
    
    out.is_doable = 1;
    return out;

    /** end of ETR_EXPRESSION **/
  case ETR_NUMBER :
    
    /** too easy. This is a non-index, query associated thing **/

    out.is_doable = 1;
    out.is_index  = 0;
    out.is_query  = 1;
    
    if( is_attachable == FALSE ) {
      /** we are in an array/method call, just chain out for them to know **/

      return out;
    }

    /** ok, lets attach it now! **/
    
    if( ot->tprf_unit != NULL ) {
      warn("This is annoying problem. We have a very simple grammar involved in compugen things, and can't cope with 2 tprfs. We should!!");
      out.is_doable = 0;
      return out;
    }

    unit = OmUnit_alloc();

    /*** should put away C loop now ***/
    
    unit->base = (-1); /** unassigned for the moment **/
    unit->length = 1;
    unit->is_tprf = 1;
    out.unit = unit;
    ot->tprf_unit = unit;

    return out;

    /*** end of NUMBER ***/

  case ETR_NAME :

    
    /** If it is top level, need to scope to query/target/resource/extern **/

    if( (et->attrib & IS_TOPLEVEL) == IS_TOPLEVEL) {
      
      su = ScopeUnit_from_Scope(sc,et->word);
      if( su == NULL ) {
	warn("Currently we can't cope with implicit externs in compugen funcs");
	out.is_doable = 0;
	return out;
      }
      
      if( su->scope_type == SCOPE_TARGET ) {
	out.is_doable   = 1;
	out.is_variable = 1;
	out.is_query    = 0;
	return out;
      }
      
      /*** otherwise, if this attachable, then - lets attach to tprf **/
      out.is_doable = 1;	
      out.is_index  = 0;
      out.is_variable = 0;
      out.is_query  = 1;
    
      if( is_attachable == FALSE ) {
	/** we are in an array/method call, just chain out for them to know **/
	out.is_variable = 1;
	return out;
      }


      /** we should look at the type here... 
	is this really something we can "add"**/

      /** ok, lets attach it now! **/
    
      if( ot->tprf_unit != NULL ) {
	warn("This is annoying problem. We have a very simple grammar involved in compugen things, and can't cope with 2 tprfs. We should!!");
	out.is_doable = 0;
	return out;
      }

      unit = OmUnit_alloc();

      /*** should put away C loop now ***/
      
      unit->base = (-1); /** unassigned for the moment **/
      unit->length = 1;
      unit->is_tprf = 1;
      out.unit = unit;
      ot->tprf_unit = unit;
      
      return out;
      
      /*** end of if TOP_LEVEL  ***/
    }

    /** is not top level... do we care? **/

    if( is_attachable == TRUE ) {
      warn("I feel bad about this: we have a non-top level name which you claim is attachable. Should not be!");
      
    }
    return out;

    /** end of ETR_NAME **/

  case ETR_TAG :

    /** right... need to loop into tag... **/

    /** for the moment, handle only 1 length tags. Should be most **/

    if( et->nochild > 1 ) {
      warn("Can't currently cope with tags with more than one child. Yikes!");
      out.is_doable = 0;
      return out;
    }

    child = descend_ExprTree_Cugen(et->child[0],mts,is_attachable,sc,ot);
    return child; /** get out of here now! **/

  case ETR_ARRAY :

    /** ok, lets see the first part of the array **/

    child = descend_ExprTree_Cugen(et->child[0],mts,FALSE,sc,ot);
    if( child.is_doable == 0 )
      return child;

    if( child.is_query == 0 ) {
      /** is a target type. No! **/

      warn("Certainly can't array into target type for compugen port");
      out.is_doable = 0;
      return out;
    }

    /** should be a query type, not indexed. Lets hope so! **/

    /** ok, now the expression into the array. Lets see what it is **/

    child = descend_ExprTree_Cugen(et->child[1],mts,FALSE,sc,ot);

    printf("Child of array into is %d\n",child.is_doable);
    if( child.is_doable == 0 )
      return child;


    /** if child is query ... then we have a tprf ***/

    if( child.is_query == 1 ) {

      if( ot->tprf_unit != NULL ) {
	warn("This is annoying problem. We have a very simple grammar involved in compugen things, and can't cope with 2 tprfs. We should!!");
	out.is_doable = 0;
	return out;
      }

      unit = OmUnit_alloc();

      /*** should put away C loop now ***/
      
      unit->base = (-1); /** unassigned for the moment **/
      unit->length = 1;
      unit->is_tprf = 1;
      out.unit = unit;
      ot->tprf_unit = unit;

      out.is_doable = 1;

      return out;
    } else {
      warn("Not sure what to do with non query typed array expressions");
      return out;
    }
      
    warn("Should not have got here!!!!");
    return out;

    

  default :
    warn("Fell into an uncatched ETR_ type [%d]. Does not look good for you!",et->type);
    out.is_doable = 0;
    return out;
  } /* end of switch */

}


   
      

# line 484 "dy2om.dy"
void set_CugenYaccReturn(CugenYaccReturn * ret)
{
  ret->is_doable = 0;
  ret->is_index  = 0;
  ret->is_query  = 0;
  ret->unit      = NULL;
}





    /*
      Debugging/test functions, not going to be
      used in final port 
      
      */



/* Function:  write_OneModel(om,ofp)
 *
 * Descrip:    writes the model definition
 *
 *
 * Arg:         om [UNKN ] Undocumented argument [OneModel *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 507 "dy2om.dy"
void write_OneModel(OneModel * om,FILE * ofp)
{
    int i;

    fprintf(ofp,"FILE_TYPE=MODEL;\n");

    /*** boring stuff ***/

    fprintf(ofp,"transition_name from to  dx dy t_prf dprf t_seq dseq w_base w_seq dwx dwy\n");

    for(i=0;i<om->len;i++) {
	write_OmTrans(i,om->trans[i],ofp);
    }

}



/* Function:  write_OmTrans(number,ot,ofp)
 *
 * Descrip:    writes one transition line
 *
 *
 * Arg:        number [UNKN ] Undocumented argument [int]
 * Arg:            ot [UNKN ] Undocumented argument [OmTrans *]
 * Arg:           ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 528 "dy2om.dy"
void write_OmTrans(int number,OmTrans * ot,FILE * ofp)
{
    fprintf(ofp,"Trans%d  %s %s %d %d %d %d %d %d %d %d %d %d\n",number,
	    ot->from->name,
	    ot->to->name,
	    ot->dx,
	    ot->dy,
	    ot->t_prf,
	    ot->dprf,
	    ot->t_seq,
	    ot->dseq,
	    ot->w_base,	   
	    ot->w_seq,
	    ot->dwx,
	    ot->dwy );
}

# line 501 "dy2om.c"
/* Function:  hard_link_OmState(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [OmState *]
 *
 * Return [UNKN ]  Undocumented return value [OmState *]
 *
 */
OmState * hard_link_OmState(OmState * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a OmState object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  OmState_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [OmState *]
 *
 */
OmState * OmState_alloc(void) 
{
    OmState * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(OmState *) ckalloc (sizeof(OmState))) == NULL)  {  
      warn("OmState_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->name = NULL;    
    out->type = OmState_TYPE_UNKNOWN;    


    return out;  
}    


/* Function:  free_OmState(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [OmState *]
 *
 * Return [UNKN ]  Undocumented return value [OmState *]
 *
 */
OmState * free_OmState(OmState * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a OmState obj. Should be trappable");   
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


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_OmUnit(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [OmUnit *]
 *
 * Return [UNKN ]  Undocumented return value [OmUnit *]
 *
 */
OmUnit * hard_link_OmUnit(OmUnit * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a OmUnit object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  OmUnit_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [OmUnit *]
 *
 */
OmUnit * OmUnit_alloc(void) 
{
    OmUnit * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(OmUnit *) ckalloc (sizeof(OmUnit))) == NULL)    {  
      warn("OmUnit_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->c_string = NULL;    
    out->base = 0;   
    out->length = 0; 
    out->is_tprf = 0;    


    return out;  
}    


/* Function:  free_OmUnit(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [OmUnit *]
 *
 * Return [UNKN ]  Undocumented return value [OmUnit *]
 *
 */
OmUnit * free_OmUnit(OmUnit * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a OmUnit obj. Should be trappable");    
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
    if( obj->c_string != NULL)   
      ckfree(obj->c_string);     


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_OmTrans(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [OmTrans *]
 *
 * Return [UNKN ]  Undocumented return value [OmTrans *]
 *
 */
OmTrans * hard_link_OmTrans(OmTrans * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a OmTrans object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  OmTrans_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [OmTrans *]
 *
 */
OmTrans * OmTrans_alloc(void) 
{
    OmTrans * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(OmTrans *) ckalloc (sizeof(OmTrans))) == NULL)  {  
      warn("OmTrans_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->dx = 0; 
    out->dy = 0; 
    out->t_prf = -1; 
    out->dprf = 0;   
    out->t_seq = -1; 
    out->dseq = 0;   
    out->w_base = -1;    
    out->w_seq = 0;  
    out->dwx = 0;    
    out->dwy = 0;    
    out->tprf_unit = NULL;   
    out->wbase_unit = NULL;  
    out->current_base = 0;   


    return out;  
}    


/* Function:  free_OmTrans(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [OmTrans *]
 *
 * Return [UNKN ]  Undocumented return value [OmTrans *]
 *
 */
OmTrans * free_OmTrans(OmTrans * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a OmTrans obj. Should be trappable");   
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
    /* obj->from is linked in */ 
    /* obj->to is linked in */ 
    if( obj->tprf_unit != NULL)  
      free_OmUnit(obj->tprf_unit);   
    if( obj->wbase_unit != NULL) 
      free_OmUnit(obj->wbase_unit);  


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_OmTransFunc(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [OmTransFunc *]
 *
 * Return [UNKN ]  Undocumented return value [OmTransFunc *]
 *
 */
OmTransFunc * hard_link_OmTransFunc(OmTransFunc * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a OmTransFunc object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  OmTransFunc_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [OmTransFunc *]
 *
 */
OmTransFunc * OmTransFunc_alloc(void) 
{
    OmTransFunc * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(OmTransFunc *) ckalloc (sizeof(OmTransFunc))) == NULL)  {  
      warn("OmTransFunc_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->name = NULL;    
    out->me = NULL;  


    return out;  
}    


/* Function:  free_OmTransFunc(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [OmTransFunc *]
 *
 * Return [UNKN ]  Undocumented return value [OmTransFunc *]
 *
 */
OmTransFunc * free_OmTransFunc(OmTransFunc * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a OmTransFunc obj. Should be trappable");   
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
    if( obj->me != NULL) 
      free_Method(obj->me);  


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_OneModel(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_OneModel
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [OmTrans **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_OneModel(OmTrans ** list,int i,int j)  
{
    OmTrans * temp;  
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_OneModel(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_OneModel which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [OmTrans **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_OneModel(OmTrans ** list,int left,int right,int (*comp)(OmTrans * ,OmTrans * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_OneModel(list,left,(left+right)/2); 
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_OneModel (list,++last,i);   
      }  
    swap_OneModel (list,left,last);  
    qsort_OneModel(list,left,last-1,comp);   
    qsort_OneModel(list,last+1,right,comp);  
}    


/* Function:  sort_OneModel(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_OneModel
 *
 *
 * Arg:         obj [UNKN ] Object containing list [OneModel *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_OneModel(OneModel * obj,int (*comp)(OmTrans *, OmTrans *)) 
{
    qsort_OneModel(obj->trans,0,obj->len-1,comp);    
    return;  
}    


/* Function:  expand_OneModel(obj,len)
 *
 * Descrip:    Really an internal function for add_OneModel
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [OneModel *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_OneModel(OneModel * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_OneModel called with no need");   
      return TRUE;   
      }  


    if( (obj->trans = (OmTrans ** ) ckrealloc (obj->trans,sizeof(OmTrans *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_OneModel, returning FALSE"); 
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_OneModel(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [OneModel *]
 * Arg:        add [OWNER] Object to add to the list [OmTrans *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_OneModel(OneModel * obj,OmTrans * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_OneModel(obj,obj->len + OneModelLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->trans[obj->len++]=add;  
    return TRUE; 
}    


/* Function:  flush_OneModel(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [OneModel *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_OneModel(OneModel * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->trans[i] != NULL) {  
        free_OmTrans(obj->trans[i]); 
        obj->trans[i] = NULL;    
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  swap_st_OneModel(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_st_OneModel
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [OmState **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_st_OneModel(OmState ** list,int i,int j)  
{
    OmState * temp;  
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_st_OneModel(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_st_OneModel which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [OmState **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_st_OneModel(OmState ** list,int left,int right,int (*comp)(OmState * ,OmState * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_st_OneModel(list,left,(left+right)/2);  
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_st_OneModel (list,++last,i);    
      }  
    swap_st_OneModel (list,left,last);   
    qsort_st_OneModel(list,left,last-1,comp);    
    qsort_st_OneModel(list,last+1,right,comp);   
}    


/* Function:  sort_st_OneModel(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_st_OneModel
 *
 *
 * Arg:         obj [UNKN ] Object containing list [OneModel *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_st_OneModel(OneModel * obj,int (*comp)(OmState *, OmState *)) 
{
    qsort_st_OneModel(obj->state,0,obj->st_len-1,comp);  
    return;  
}    


/* Function:  expand_st_OneModel(obj,len)
 *
 * Descrip:    Really an internal function for add_st_OneModel
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [OneModel *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_st_OneModel(OneModel * obj,int len) 
{


    if( obj->st_maxlen > obj->st_len )   {  
      warn("expand_OneModelst_ called with no need");    
      return TRUE;   
      }  


    if( (obj->state = (OmState ** ) ckrealloc (obj->state,sizeof(OmState *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_OneModel, returning FALSE"); 
      return FALSE;  
      }  
    obj->st_maxlen = len;    
    return TRUE; 
}    


/* Function:  add_st_OneModel(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [OneModel *]
 * Arg:        add [OWNER] Object to add to the list [OmState *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_st_OneModel(OneModel * obj,OmState * add) 
{
    if( obj->st_len >= obj->st_maxlen)   {  
      if( expand_st_OneModel(obj,obj->st_len + OneModelLISTLENGTH) == FALSE) 
        return FALSE;    
      }  


    obj->state[obj->st_len++]=add;   
    return TRUE; 
}    


/* Function:  flush_st_OneModel(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [OneModel *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_st_OneModel(OneModel * obj) 
{
    int i;   


    for(i=0;i<obj->st_len;i++)   { /*for i over list length*/ 
      if( obj->state[i] != NULL) {  
        free_OmState(obj->state[i]); 
        obj->state[i] = NULL;    
        }  
      } /* end of for i over list length */ 


    obj->st_len = 0; 
    return i;    
}    


/* Function:  swap_tf_OneModel(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_tf_OneModel
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [OmTransFunc **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_tf_OneModel(OmTransFunc ** list,int i,int j)  
{
    OmTransFunc * temp;  
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_tf_OneModel(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_tf_OneModel which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [OmTransFunc **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_tf_OneModel(OmTransFunc ** list,int left,int right,int (*comp)(OmTransFunc * ,OmTransFunc * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_tf_OneModel(list,left,(left+right)/2);  
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_tf_OneModel (list,++last,i);    
      }  
    swap_tf_OneModel (list,left,last);   
    qsort_tf_OneModel(list,left,last-1,comp);    
    qsort_tf_OneModel(list,last+1,right,comp);   
}    


/* Function:  sort_tf_OneModel(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_tf_OneModel
 *
 *
 * Arg:         obj [UNKN ] Object containing list [OneModel *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_tf_OneModel(OneModel * obj,int (*comp)(OmTransFunc *, OmTransFunc *)) 
{
    qsort_tf_OneModel(obj->tfunc,0,obj->tf_len-1,comp);  
    return;  
}    


/* Function:  expand_tf_OneModel(obj,len)
 *
 * Descrip:    Really an internal function for add_tf_OneModel
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [OneModel *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_tf_OneModel(OneModel * obj,int len) 
{


    if( obj->tf_maxlen > obj->tf_len )   {  
      warn("expand_OneModeltf_ called with no need");    
      return TRUE;   
      }  


    if( (obj->tfunc = (OmTransFunc ** ) ckrealloc (obj->tfunc,sizeof(OmTransFunc *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_OneModel, returning FALSE"); 
      return FALSE;  
      }  
    obj->tf_maxlen = len;    
    return TRUE; 
}    


/* Function:  add_tf_OneModel(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [OneModel *]
 * Arg:        add [OWNER] Object to add to the list [OmTransFunc *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_tf_OneModel(OneModel * obj,OmTransFunc * add) 
{
    if( obj->tf_len >= obj->tf_maxlen)   {  
      if( expand_tf_OneModel(obj,obj->tf_len + OneModelLISTLENGTH) == FALSE) 
        return FALSE;    
      }  


    obj->tfunc[obj->tf_len++]=add;   
    return TRUE; 
}    


/* Function:  flush_tf_OneModel(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [OneModel *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_tf_OneModel(OneModel * obj) 
{
    int i;   


    for(i=0;i<obj->tf_len;i++)   { /*for i over list length*/ 
      if( obj->tfunc[i] != NULL) {  
        free_OmTransFunc(obj->tfunc[i]); 
        obj->tfunc[i] = NULL;    
        }  
      } /* end of for i over list length */ 


    obj->tf_len = 0; 
    return i;    
}    


/* Function:  OneModel_alloc_std(void)
 *
 * Descrip:    Equivalent to OneModel_alloc_len(OneModelLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [OneModel *]
 *
 */
OneModel * OneModel_alloc_std(void) 
{
    return OneModel_alloc_len(OneModelLISTLENGTH);   
}    


/* Function:  OneModel_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [OneModel *]
 *
 */
OneModel * OneModel_alloc_len(int len) 
{
    OneModel * out; /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = OneModel_alloc()) == NULL) 
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->trans = (OmTrans ** ) ckcalloc (len,sizeof(OmTrans *))) == NULL)    {  
      warn("Warning, ckcalloc failed in OneModel_alloc_len");    
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    if((out->state = (OmState ** ) ckcalloc (len,sizeof(OmState *))) == NULL)    {  
      warn("Warning, ckcalloc failed in OneModel_alloc_len");    
      return NULL;   
      }  
    out->st_len = 0; 
    out->st_maxlen = len;    


    if((out->tfunc = (OmTransFunc ** ) ckcalloc (len,sizeof(OmTransFunc *))) == NULL)    {  
      warn("Warning, ckcalloc failed in OneModel_alloc_len");    
      return NULL;   
      }  
    out->tf_len = 0; 
    out->tf_maxlen = len;    


    return out;  
}    


/* Function:  hard_link_OneModel(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [OneModel *]
 *
 * Return [UNKN ]  Undocumented return value [OneModel *]
 *
 */
OneModel * hard_link_OneModel(OneModel * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a OneModel object: passed a NULL object");    
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  OneModel_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [OneModel *]
 *
 */
OneModel * OneModel_alloc(void) 
{
    OneModel * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(OneModel *) ckalloc (sizeof(OneModel))) == NULL)    {  
      warn("OneModel_alloc failed ");    
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->name = NULL;    
    out->state_num = 0;  
    out->semistate_num = 0;  
    out->midstate_num = 0;   
    out->profile_line_size = 0;  
    out->trans = NULL;   
    out->len = out->maxlen = 0;  
    out->state = NULL;   
    out->st_len = out->st_maxlen = 0;    
    out->tfunc = NULL;   
    out->tf_len = out->tf_maxlen = 0;    


    return out;  
}    


/* Function:  free_OneModel(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [OneModel *]
 *
 * Return [UNKN ]  Undocumented return value [OneModel *]
 *
 */
OneModel * free_OneModel(OneModel * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a OneModel obj. Should be trappable");  
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
    if( obj->trans != NULL)  {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->trans[i] != NULL)   
          free_OmTrans(obj->trans[i]);   
        }  
      ckfree(obj->trans);    
      }  
    if( obj->state != NULL)  {  
      for(i=0;i<obj->st_len;i++) {  
        if( obj->state[i] != NULL)   
          free_OmState(obj->state[i]);   
        }  
      ckfree(obj->state);    
      }  
    if( obj->tfunc != NULL)  {  
      for(i=0;i<obj->tf_len;i++) {  
        if( obj->tfunc[i] != NULL)   
          free_OmTransFunc(obj->tfunc[i]);   
        }  
      ckfree(obj->tfunc);    
      }  


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif
