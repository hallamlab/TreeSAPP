#ifdef _cplusplus
extern "C" {
#endif
#include "standardout.h"




# line 24 "standardout.dy"
void show_help_StandardOutputOptions(FILE * ofp)
{
  fprintf(ofp,"Standard alignment outputs\n");
  fprintf(ofp,"   -alb          show align block format\n");
  fprintf(ofp,"   -pal          show raw alignment\n");
  fprintf(ofp,"   -calb         show cumlative align block\n");
  fprintf(ofp,"   -cpal         show cumlative raw alignemnt\n");

}


# line 35 "standardout.dy"
StandardOutputOptions * new_StandardOutputOptions_from_argv(int * argc,char ** argv)
{
  StandardOutputOptions * out;

  out = StandardOutputOptions_alloc();

  out->show_alb = strip_out_boolean_argument(argc,argv,"alb");
  out->show_pal = strip_out_boolean_argument(argc,argv,"pal");
  out->show_cumlative_alb = strip_out_boolean_argument(argc,argv,"calb");
  out->show_cumlative_pal = strip_out_boolean_argument(argc,argv,"cpal");

  return out;
}


# line 50 "standardout.dy"
void show_StandardOutputOptions(StandardOutputOptions * out,AlnBlock * alb,PackAln * pal,char * divide_str,FILE * ofp)
{
  assert(out);
  assert(alb);
  assert(pal);
  assert(ofp);
  assert(divide_str);

  if( out->show_alb == TRUE ) {
    mapped_ascii_AlnBlock(alb,Score2Bits,0,ofp);
    fprintf(ofp,"%s\n",divide_str);
  }

  if( out->show_cumlative_alb == TRUE ) {
    mapped_ascii_AlnBlock(alb,Score2Bits,1,ofp);
    fprintf(ofp,"%s\n",divide_str);
  }


  if( out->show_cumlative_pal == TRUE ) {
    show_bits_and_cumlative_PackAln(pal,ofp);
    fprintf(ofp,"%s\n",divide_str);
  }

  if( out->show_pal == TRUE ) {
    show_simple_PackAln(pal,ofp);
    fprintf(ofp,"%s\n",divide_str);
  }

  return;
}
# line 68 "standardout.c"
/* Function:  hard_link_StandardOutputOptions(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [StandardOutputOptions *]
 *
 * Return [UNKN ]  Undocumented return value [StandardOutputOptions *]
 *
 */
StandardOutputOptions * hard_link_StandardOutputOptions(StandardOutputOptions * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a StandardOutputOptions object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  StandardOutputOptions_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [StandardOutputOptions *]
 *
 */
StandardOutputOptions * StandardOutputOptions_alloc(void) 
{
    StandardOutputOptions * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(StandardOutputOptions *) ckalloc (sizeof(StandardOutputOptions))) == NULL)  {  
      warn("StandardOutputOptions_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->show_alb = FALSE;   
    out->show_pal = FALSE;   
    out->show_cumlative_alb = FALSE; 
    out->show_cumlative_pal = FALSE; 


    return out;  
}    


/* Function:  free_StandardOutputOptions(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [StandardOutputOptions *]
 *
 * Return [UNKN ]  Undocumented return value [StandardOutputOptions *]
 *
 */
StandardOutputOptions * free_StandardOutputOptions(StandardOutputOptions * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a StandardOutputOptions obj. Should be trappable"); 
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



#ifdef _cplusplus
}
#endif
