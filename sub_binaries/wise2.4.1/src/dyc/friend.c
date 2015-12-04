#ifdef _cplusplus
extern "C" {
#endif
#include "friend.h"


# line 15 "friend.dy"
Friend * read_Friend_line(char * line,FILE * ifp)
{
  Friend * out;
  char * runner;

  if( strstartcmp(line,"friend") != 0 ) {
    warn("Trying to read [%s] as a friend line? I think not!");
  }

  runner = strtok(line,spacestr);
  runner = strtok(NULL,spacestr);

  if( runner == NULL ) {
    warn("Could not find any friend for a friend line. Sad!");
  }

  out = Friend_alloc();

  out->name = stringalloc(runner);

  return out;
}


# line 39 "friend.dy"
void write_Friend_header(DYNFILE * dfp,Friend * fr)
{
  if( dfp->package_name != NULL ) {
    fprintf(dfp->head,"#ifndef DYNAMITE_DEFINED_%s\n",fr->name);
    fprintf(dfp->head,"typedef struct %s%s %s%s;\n",dfp->package_name,fr->name,dfp->package_name,fr->name);
    fprintf(dfp->head,"#define %s %s%s\n",fr->name,dfp->package_name,fr->name);
    fprintf(dfp->head,"#define DYNAMITE_DEFINED_%s\n#endif\n\n",fr->name);
  } else 
  fprintf(dfp->head,"#ifndef DYNAMITE_DEFINED_%s\ntypedef struct %s %s;\n#define DYNAMITE_DEFINED_%s\n#endif\n\n",fr->name,fr->name,fr->name,fr->name);
}
# line 42 "friend.c"
/* Function:  hard_link_Friend(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [Friend *]
 *
 * Return [UNKN ]  Undocumented return value [Friend *]
 *
 */
Friend * hard_link_Friend(Friend * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a Friend object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  Friend_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Friend *]
 *
 */
Friend * Friend_alloc(void) 
{
    Friend * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(Friend *) ckalloc (sizeof(Friend))) == NULL)    {  
      warn("Friend_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->name = NULL;    


    return out;  
}    


/* Function:  free_Friend(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [Friend *]
 *
 * Return [UNKN ]  Undocumented return value [Friend *]
 *
 */
Friend * free_Friend(Friend * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a Friend obj. Should be trappable");    
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



#ifdef _cplusplus
}
#endif
