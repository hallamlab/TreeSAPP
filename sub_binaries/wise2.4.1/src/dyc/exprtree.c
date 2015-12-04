#ifdef _cplusplus
extern "C" {
#endif
#include "exprtree.h"

# line 99 "exprtree.dy"
void parentfy_ExprTree(ExprTree * et)
{
  int i;

  for(i=0;i<et->nochild;i++) {
    et->child[i]->parent = et;
    et->child[i]->position_in_parent = i;
    parentfy_ExprTree(et->child[i]);
  }

 
}

# line 112 "exprtree.dy"
void declared_ExprTree(ExprTree * et)
{
  fprintf(stderr,"Declaring...");

  et->attrib = et->attrib | IS_DECLARED;

}
  
  
    
# line 122 "exprtree.dy"
void find_toplevel_name(ExprTree * et)
{
  /** et should be an expression or a tag **/

  int i;

  for(i=0;i<et->nochild;i++) {
    switch(et->child[i]->type) {
    case ETR_NAME :
      if( et->position_in_parent == 0 || !(et->parent->type == ETR_REFERENCE || et->parent->type == ETR_STRUCTREF) ) {
	/*      	printf("Found as a name -> [%d][%d]",et->parent->position_in_parent,et->parent->type);
	print_ExprTree(et->child[i]);
	printf("\n");
	*/
	et->child[i]->attrib |= IS_TOPLEVEL;
      } else {
       
	/*	printf("Ignoring name %d,%d",et->parent->position_in_parent,et->parent->type);
	print_ExprTree(et->child[i]);
	printf("\n");
	*/

	et->child[i]->attrib &= ~IS_TOPLEVEL;

      }
      break;
    default :
/*
      printf("Going into");
      print_ExprTree(et->child[i]);
      printf("\n");
*/
      find_toplevel_name(et->child[i]);
      break;
    }
  }
}


# line 161 "exprtree.dy"
void strcat_ExprTree(ExprTree * ExprTree,char * buffer)
{
  int i;

  switch(ExprTree->type) {
  case ETR_NUMBER : strcat(buffer,ExprTree->word); return;
  case ETR_OPERATOR : strcat(buffer,ExprTree->word); return;
  case ETR_EXPRESSION : strcat(buffer,"(");
    for(i=0;i< ExprTree->nochild;i++)
      strcat_ExprTree(ExprTree->child[i],buffer);
    strcat(buffer,")");
    return;
  case ETR_STATEMENT :
    for(i=0;i< ExprTree->nochild;i++)
      strcat_ExprTree(ExprTree->child[i],buffer);
    return;
  case ETR_NAME : strcat(buffer,ExprTree->word);
    return;
  case ETR_ARRAY : 
    strcat_ExprTree(ExprTree->child[0],buffer);
    strcat(buffer,"[");
    strcat_ExprTree(ExprTree->child[1],buffer);
    strcat(buffer,"]");
    return;
  case ETR_TAG : 
    for(i=0;i< ExprTree->nochild;i++)
      strcat_ExprTree(ExprTree->child[i],buffer);
    return;
  case ETR_STRUCTREF :
    strcat_ExprTree(ExprTree->child[0],buffer);
    strcat(buffer,ExprTree->child[1]->word);
    strcat_ExprTree(ExprTree->child[2],buffer);
    return;
  case ETR_REFERENCE :
    strcat(buffer,ExprTree->child[0]->word);
    strcat_ExprTree(ExprTree->child[1],buffer);
    return;
  case ETR_METHOD :
    fprintf(stderr,"So buffer is now [%s]\n",buffer);
    strcat_ExprTree(ExprTree->child[0],buffer);
    strcat(buffer,"(");
    fprintf(stderr,"So buffer is now [%s]\n",buffer);
    strcat_ExprTree(ExprTree->child[1],buffer);

    /*    for(i=0;i<ExprTree->child[1]->nochild;i++) 
      strcat_ExprTree(ExprTree->child[1]->child[i],buffer);
      */

    strcat(buffer,")");
    return;
  case ETR_COMMALIST :
    for(i=0;i<ExprTree->nochild;i++) {
      strcat_ExprTree(ExprTree->child[i],buffer);
      if( i != ExprTree->nochild-1)
	strcat(buffer,",");
    }
    return;

  default :
      warn("In trying to make Expr string, got unobtainable type!");
  }
}
  
# line 224 "exprtree.dy"
void print_ExprTree(ExprTree * ExprTree)
{
  int i;

  switch(ExprTree->type) {
  case ETR_NUMBER : printf("#{%s}",ExprTree->word); return;
  case ETR_OPERATOR : printf("Op{%s}",ExprTree->word); return;
  case ETR_EXPRESSION : printf("Expression (");
    for(i=0;i< ExprTree->nochild;i++)
      print_ExprTree(ExprTree->child[i]);
    printf(")");
    return;
  case ETR_STATEMENT : printf("Top level:");
    for(i=0;i< ExprTree->nochild;i++)
      print_ExprTree(ExprTree->child[i]);
    printf("\n");
    return;
  case ETR_NAME : printf("N{%s}",ExprTree->word);
    for(i=0;i< ExprTree->nochild;i++)
      print_ExprTree(ExprTree->child[i]);
    return;
  case ETR_ARRAY : printf("{ArrayInto{");
    print_ExprTree(ExprTree->child[0]);
    printf("} of [");
    print_ExprTree(ExprTree->child[1]);
    printf("]}");
    return;
  case ETR_TAG : 
    printf("Tag{");
    for(i=0;i< ExprTree->nochild;i++)
      print_ExprTree(ExprTree->child[i]);
    printf("}");
    return;
  case ETR_STRUCTREF :
    printf("struct{%s}(",ExprTree->child[1]->word);
    print_ExprTree(ExprTree->child[0]);
    printf(":");
    print_ExprTree(ExprTree->child[2]);
    printf(")");
    return;
  case ETR_REFERENCE :
    printf("Ref{");  
    puts(ExprTree->child[0]->word);
    print_ExprTree(ExprTree->child[1]);
    printf("}\n");
    return;
  case ETR_METHOD :
    if( ExprTree->attrib & IS_DECLARED ) 
      printf("DEC:");

    printf("Method{");
    for(i=0;i<ExprTree->child[0]->nochild;i++) 
      print_ExprTree(ExprTree->child[0]->child[i]);
    if( ExprTree->nochild == 1 )
      printf("}(void)");
    else {
      printf("}(");
      for(i=0;i<ExprTree->child[1]->nochild;i++) 
	print_ExprTree(ExprTree->child[1]->child[i]);
      printf(")");
    }
    return;
  case ETR_COMMALIST :
    printf("List{");
    for(i=0;i<ExprTree->nochild;i++) 
      print_ExprTree(ExprTree->child[i]);
    printf("}");
    return;

  default :
      printf("Unprintable!");
  }
}


/*** declarations ***/

# line 301 "exprtree.dy"
ExprTree * new_ExprTree_decl_method(ExprTree * name,ExprTree * list)
{
  ExprTree * out;
  
  out = new_ExprTree();

  out->type = ETR_DECL_METHOD;
  out->child[0]=name;
  out->child[1]=list;
  out->nochild = 2;

  return out;
}

# line 315 "exprtree.dy"
ExprTree * new_ExprTree_decl_variable(ExprTree * type,ExprTree * name)
{
  ExprTree * out;
  
  out = new_ExprTree();

  out->type = ETR_DECL_VARIABLE;
  out->child[0]=type;
  out->child[1]=name;
  out->nochild = 2;

  return out;
}

# line 329 "exprtree.dy"
ExprTree * add_to_decl_list_ExprTree(ExprTree * list,ExprTree * add)
{
  if( list->type != ETR_DECL_LIST ) {
    warn("Attempting to add to a non commalist %d",list->type);
    return list;
  }

  add_ExprTree(list,add);

  return list;
}

# line 341 "exprtree.dy"
ExprTree * new_ExprTree_decl_list(ExprTree * start)
{
  ExprTree * out;
  
  out = new_ExprTree();

  out->type = ETR_DECL_LIST;
  add_ExprTree(out,start);

  return out;
}

# line 353 "exprtree.dy"
ExprTree * new_ExprTree_struct_ref(ExprTree * left, ExprTree * ref,ExprTree * right)
{
  ExprTree * out;
  
  out = new_ExprTree();

  out->type = ETR_STRUCTREF;
  out->child[0]=left;
  out->child[1]=ref;
  out->child[2]=right;
  out->nochild = 3;

  return out;
}

# line 368 "exprtree.dy"
ExprTree * new_ExprTree_ref(char op,ExprTree * right)
{
  ExprTree * out;
  
  out = new_ExprTree();

  out->type = ETR_REFERENCE;
  out->child[0]=new_ExprTree_token(op);
  out->child[1]=right;
  out->nochild = 2;

  return out;
}

# line 382 "exprtree.dy"
ExprTree * add_to_commalist_ExprTree(ExprTree * list,ExprTree * add)
{
  if( list->type != ETR_COMMALIST ) {
    warn("Attempting to add to a non commalist %d",list->type);
    return list;
  }

  add_ExprTree(list,add);

  return list;
}

# line 394 "exprtree.dy"
ExprTree * new_ExprTree_commalist(ExprTree * start)
{
  ExprTree * out;
  
  out = new_ExprTree();

  out->type = ETR_COMMALIST;
  add_ExprTree(out,start);

  return out;
}

# line 406 "exprtree.dy"
ExprTree * new_ExprTree_method(ExprTree * one,ExprTree * other)
{
  ExprTree * out;

  out = ExprTree_alloc();
  out->type = ETR_METHOD;

  /* printf("This one %d\n",one->nochild); */

  add_ExprTree(out,one);

  if(other == NULL) {
    return out;
  }

  if(other->type != ETR_COMMALIST ) {
    warn("Attempting to have a non commalist argument for a method");
  }

  add_ExprTree(out,other);

  return out;
}

# line 430 "exprtree.dy"
ExprTree * new_ExprTree_tag_from_name(ExprTree * name)
{

  ExprTree * out;

  out= new_ExprTree();
  out->type = ETR_TAG;

  add_ExprTree(out,name);


  return out;
}

# line 444 "exprtree.dy"
ExprTree * new_ExprTree_array(ExprTree * tag,ExprTree * expr) 
{
  ExprTree * out;

  out = new_ExprTree();
  out->type = ETR_ARRAY;
  out->child[0] = tag;
  out->child[1] = expr;
  out->nochild = 2;
  return out;
}

# line 456 "exprtree.dy"
ExprTree * new_ExprTree_binary_expr(ExprTree * left,char op,ExprTree * rgt)
{
  ExprTree * out;

  
  out = new_ExprTree();

  out->type = ETR_EXPRESSION;
  out->child[0] = left;
  out->child[1] = new_ExprTree_token(op);
  out->child[2] = rgt;
  out->nochild = 3;

  return out;
}

# line 472 "exprtree.dy"
ExprTree * new_ExprTree_token(char t)
{
  ExprTree * out;
  char buf[2];

  buf[0] = t;
  buf[1] = '\0';

  out= new_ExprTree();
  out->type= ETR_OPERATOR;
  out->token = t;
  out->word = stringalloc(buf);

  return out;
}

# line 488 "exprtree.dy"
boolean add_ExprTree(ExprTree * one,ExprTree * child)
{
  if( one->nochild+1 > EXPRTREE_MAXCHILD ) {
    warn("Overflow in ExprTree at %d children",EXPRTREE_MAXCHILD);
    return FALSE;
  }

  one->child[one->nochild++] = child;
  return TRUE;
}

# line 499 "exprtree.dy"
ExprTree * new_ExprTree(void)
{
  int i;
  ExprTree * out;
  
  out = ExprTree_alloc();
  for(i=0;i<128;i++)
    out->child[i]= NULL;

 
  return out;
}



# line 440 "exprtree.c"
/* Function:  hard_link_ExprTree(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [ExprTree *]
 *
 * Return [UNKN ]  Undocumented return value [ExprTree *]
 *
 */
ExprTree * hard_link_ExprTree(ExprTree * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a ExprTree object: passed a NULL object");    
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  ExprTree_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ExprTree *]
 *
 */
ExprTree * ExprTree_alloc(void) 
{
    ExprTree * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(ExprTree *) ckalloc (sizeof(ExprTree))) == NULL)    {  
      warn("ExprTree_alloc failed ");    
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    /* child[EXPRTREE_MAXCHILD] is an array: no default possible */ 
    out->nochild = 0;    
    out->token = 'u';    
    out->word = NULL;    
    out->type = 0;   
    out->attrib = 0; 
    out->parent = NULL;  
    out->position_in_parent = 0; 


    return out;  
}    


/* Function:  free_ExprTree(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [ExprTree *]
 *
 * Return [UNKN ]  Undocumented return value [ExprTree *]
 *
 */
ExprTree * free_ExprTree(ExprTree * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a ExprTree obj. Should be trappable");  
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
    /* Unable to make free function for obj->child[EXPRTREE_MAXCHILD] */ 
    if( obj->word != NULL)   
      ckfree(obj->word);     
    /* Unable to make free function for obj->parent */ 


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif
