#ifdef _cplusplus
extern "C" {
#endif
#include "randomdb.h"

/* Function:  Sequence_from_RandomProteinDB(rndp)
 *
 * Descrip:    Makes a new random sequence from a
 *             RandomProteinDB
 *
 *
 * Arg:        rndp [UNKN ] Undocumented argument [RandomProteinDB *]
 *
 * Return [UNKN ]  Undocumented return value [Sequence *]
 *
 */
# line 44 "randomdb.dy"
Sequence * Sequence_from_RandomProteinDB(RandomProteinDB * rndp)
{
  int i;
  Sequence * out;
  char * name = "RandomSequence";
  char buffer[512];
  char * seq;

  if( rndp->use_flat_length == FALSE ) {
    warn("Ooops have not implemented gaussian lengths yet. Yikes!");
    return NULL;
  }


  seq = ckcalloc(rndp->length+1,sizeof(char));
  for(i=0;i<rndp->length+1;i++) {
    seq[i] = draw_random_aa_RandomModel(rndp->emission);
  }

  sprintf(buffer,"%s%d",name,rndp->num++);
  out = new_Sequence_from_strings(name,seq);
  ckfree(seq);
  return out;
}

/* Function:  new_flat_RandomProteinDB(rm,length)
 *
 * Descrip:    Makes a new flat RandomProteinDB
 *
 *
 * Arg:            rm [SOFT ] RandomModel [RandomModel *]
 * Arg:        length [UNKN ] length of protein to produce [int]
 *
 * Return [UNKN ]  Undocumented return value [RandomProteinDB *]
 *
 */
# line 75 "randomdb.dy"
RandomProteinDB * new_flat_RandomProteinDB(RandomModel * rm,int length)
{
  RandomProteinDB * out;
  
  out = RandomProteinDB_alloc();
  out->use_flat_length = TRUE;
  out->emission = hard_link_RandomModel(rm);
  out->length = length;

  return out;
}


/* Function:  Sequence_from_RandomDNADB(rndp)
 *
 * Descrip:    Makes a new random sequence from a
 *             RandomDNADB
 *
 *
 * Arg:        rndp [UNKN ] Undocumented argument [RandomDNADB *]
 *
 * Return [UNKN ]  Undocumented return value [Sequence *]
 *
 */
# line 92 "randomdb.dy"
Sequence * Sequence_from_RandomDNADB(RandomDNADB * rndp)
{
  int i;
  Sequence * out;
  char * name = "RandomSequence";
  char buffer[512];
  char * seq;

  if( rndp->use_flat_length == FALSE ) {
    warn("Ooops have not implemented gaussian lengths yet. Yikes!");
    return NULL;
  }


  seq = ckcalloc(rndp->length+1,sizeof(char));
  for(i=0;i<rndp->length+1;i++) {
    seq[i] = draw_random_base_RandomModelDNA(rndp->emission);
  }

  sprintf(buffer,"%s%d",name,rndp->num++);
  out = new_Sequence_from_strings(buffer,seq);
  ckfree(seq);
  return out;
}

/* Function:  new_flat_RandomDNADB(rm,length)
 *
 * Descrip:    Makes a new flat RandomDNADB
 *
 *
 * Arg:            rm [SOFT ] RandomModel emission [RandomModelDNA *]
 * Arg:        length [UNKN ] length of DNA to produce [int]
 *
 * Return [UNKN ]  Undocumented return value [RandomDNADB *]
 *
 */
# line 123 "randomdb.dy"
RandomDNADB * new_flat_RandomDNADB(RandomModelDNA * rm,int length)
{
  RandomDNADB * out;
  
  out = RandomDNADB_alloc();
  out->use_flat_length = TRUE;
  out->emission = hard_link_RandomModelDNA(rm);
  out->length = length;
  out->num = 1;

  return out;
}

# line 121 "randomdb.c"
/* Function:  hard_link_RandomProteinDB(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [RandomProteinDB *]
 *
 * Return [UNKN ]  Undocumented return value [RandomProteinDB *]
 *
 */
RandomProteinDB * hard_link_RandomProteinDB(RandomProteinDB * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a RandomProteinDB object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  RandomProteinDB_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [RandomProteinDB *]
 *
 */
RandomProteinDB * RandomProteinDB_alloc(void) 
{
    RandomProteinDB * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(RandomProteinDB *) ckalloc (sizeof(RandomProteinDB))) == NULL)  {  
      warn("RandomProteinDB_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->use_flat_length = TRUE; 
    out->length = 0; 
    out->length_dist = NULL; 
    out->emission = NULL;    
    out->num = 0;    


    return out;  
}    


/* Function:  free_RandomProteinDB(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [RandomProteinDB *]
 *
 * Return [UNKN ]  Undocumented return value [RandomProteinDB *]
 *
 */
RandomProteinDB * free_RandomProteinDB(RandomProteinDB * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a RandomProteinDB obj. Should be trappable");   
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
    if( obj->length_dist != NULL)    
      free_Histogram(obj->length_dist);  
    if( obj->emission != NULL)   
      free_RandomModel(obj->emission);   


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_RandomDNADB(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [RandomDNADB *]
 *
 * Return [UNKN ]  Undocumented return value [RandomDNADB *]
 *
 */
RandomDNADB * hard_link_RandomDNADB(RandomDNADB * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a RandomDNADB object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  RandomDNADB_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [RandomDNADB *]
 *
 */
RandomDNADB * RandomDNADB_alloc(void) 
{
    RandomDNADB * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(RandomDNADB *) ckalloc (sizeof(RandomDNADB))) == NULL)  {  
      warn("RandomDNADB_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->use_flat_length = TRUE; 
    out->length = 0; 
    out->length_dist = NULL; 
    out->emission = NULL;    
    out->num = 0;    


    return out;  
}    


/* Function:  free_RandomDNADB(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [RandomDNADB *]
 *
 * Return [UNKN ]  Undocumented return value [RandomDNADB *]
 *
 */
RandomDNADB * free_RandomDNADB(RandomDNADB * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a RandomDNADB obj. Should be trappable");   
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
    if( obj->length_dist != NULL)    
      free_Histogram(obj->length_dist);  
    if( obj->emission != NULL)   
      free_RandomModelDNA(obj->emission);    


    ckfree(obj); 
    return NULL; 
}    


/* Function:  replace_use_flat_length_RandomProteinDB(obj,use_flat_length)
 *
 * Descrip:    Replace member variable use_flat_length
 *             For use principly by API functions
 *
 *
 * Arg:                    obj [UNKN ] Object holding the variable [RandomProteinDB *]
 * Arg:        use_flat_length [OWNER] New value of the variable [boolean]
 *
 * Return [SOFT ]  member variable use_flat_length [boolean]
 *
 */
boolean replace_use_flat_length_RandomProteinDB(RandomProteinDB * obj,boolean use_flat_length) 
{
    if( obj == NULL)     {  
      warn("In replacement function use_flat_length for object RandomProteinDB, got a NULL object"); 
      return FALSE;  
      }  
    obj->use_flat_length = use_flat_length;  
    return TRUE; 
}    


/* Function:  access_use_flat_length_RandomProteinDB(obj)
 *
 * Descrip:    Access member variable use_flat_length
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [RandomProteinDB *]
 *
 * Return [SOFT ]  member variable use_flat_length [boolean]
 *
 */
boolean access_use_flat_length_RandomProteinDB(RandomProteinDB * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function use_flat_length for object RandomProteinDB, got a NULL object");    
      return FALSE;  
      }  
    return obj->use_flat_length;     
}    


/* Function:  replace_length_RandomProteinDB(obj,length)
 *
 * Descrip:    Replace member variable length
 *             For use principly by API functions
 *
 *
 * Arg:           obj [UNKN ] Object holding the variable [RandomProteinDB *]
 * Arg:        length [OWNER] New value of the variable [int]
 *
 * Return [SOFT ]  member variable length [boolean]
 *
 */
boolean replace_length_RandomProteinDB(RandomProteinDB * obj,int length) 
{
    if( obj == NULL)     {  
      warn("In replacement function length for object RandomProteinDB, got a NULL object");  
      return FALSE;  
      }  
    obj->length = length;    
    return TRUE; 
}    


/* Function:  access_length_RandomProteinDB(obj)
 *
 * Descrip:    Access member variable length
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [RandomProteinDB *]
 *
 * Return [SOFT ]  member variable length [int]
 *
 */
int access_length_RandomProteinDB(RandomProteinDB * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function length for object RandomProteinDB, got a NULL object"); 
      return 0;  
      }  
    return obj->length;  
}    


/* Function:  replace_length_dist_RandomProteinDB(obj,length_dist)
 *
 * Descrip:    Replace member variable length_dist
 *             For use principly by API functions
 *
 *
 * Arg:                obj [UNKN ] Object holding the variable [RandomProteinDB *]
 * Arg:        length_dist [OWNER] New value of the variable [Histogram *]
 *
 * Return [SOFT ]  member variable length_dist [boolean]
 *
 */
boolean replace_length_dist_RandomProteinDB(RandomProteinDB * obj,Histogram * length_dist) 
{
    if( obj == NULL)     {  
      warn("In replacement function length_dist for object RandomProteinDB, got a NULL object"); 
      return FALSE;  
      }  
    obj->length_dist = length_dist;  
    return TRUE; 
}    


/* Function:  access_length_dist_RandomProteinDB(obj)
 *
 * Descrip:    Access member variable length_dist
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [RandomProteinDB *]
 *
 * Return [SOFT ]  member variable length_dist [Histogram *]
 *
 */
Histogram * access_length_dist_RandomProteinDB(RandomProteinDB * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function length_dist for object RandomProteinDB, got a NULL object");    
      return NULL;   
      }  
    return obj->length_dist;     
}    


/* Function:  replace_emission_RandomProteinDB(obj,emission)
 *
 * Descrip:    Replace member variable emission
 *             For use principly by API functions
 *
 *
 * Arg:             obj [UNKN ] Object holding the variable [RandomProteinDB *]
 * Arg:        emission [OWNER] New value of the variable [RandomModel *]
 *
 * Return [SOFT ]  member variable emission [boolean]
 *
 */
boolean replace_emission_RandomProteinDB(RandomProteinDB * obj,RandomModel * emission) 
{
    if( obj == NULL)     {  
      warn("In replacement function emission for object RandomProteinDB, got a NULL object");    
      return FALSE;  
      }  
    obj->emission = emission;    
    return TRUE; 
}    


/* Function:  access_emission_RandomProteinDB(obj)
 *
 * Descrip:    Access member variable emission
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [RandomProteinDB *]
 *
 * Return [SOFT ]  member variable emission [RandomModel *]
 *
 */
RandomModel * access_emission_RandomProteinDB(RandomProteinDB * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function emission for object RandomProteinDB, got a NULL object");   
      return NULL;   
      }  
    return obj->emission;    
}    


/* Function:  replace_num_RandomProteinDB(obj,num)
 *
 * Descrip:    Replace member variable num
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [RandomProteinDB *]
 * Arg:        num [OWNER] New value of the variable [int]
 *
 * Return [SOFT ]  member variable num [boolean]
 *
 */
boolean replace_num_RandomProteinDB(RandomProteinDB * obj,int num) 
{
    if( obj == NULL)     {  
      warn("In replacement function num for object RandomProteinDB, got a NULL object"); 
      return FALSE;  
      }  
    obj->num = num;  
    return TRUE; 
}    


/* Function:  access_num_RandomProteinDB(obj)
 *
 * Descrip:    Access member variable num
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [RandomProteinDB *]
 *
 * Return [SOFT ]  member variable num [int]
 *
 */
int access_num_RandomProteinDB(RandomProteinDB * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function num for object RandomProteinDB, got a NULL object");    
      return 0;  
      }  
    return obj->num;     
}    


/* Function:  replace_use_flat_length_RandomDNADB(obj,use_flat_length)
 *
 * Descrip:    Replace member variable use_flat_length
 *             For use principly by API functions
 *
 *
 * Arg:                    obj [UNKN ] Object holding the variable [RandomDNADB *]
 * Arg:        use_flat_length [OWNER] New value of the variable [boolean]
 *
 * Return [SOFT ]  member variable use_flat_length [boolean]
 *
 */
boolean replace_use_flat_length_RandomDNADB(RandomDNADB * obj,boolean use_flat_length) 
{
    if( obj == NULL)     {  
      warn("In replacement function use_flat_length for object RandomDNADB, got a NULL object"); 
      return FALSE;  
      }  
    obj->use_flat_length = use_flat_length;  
    return TRUE; 
}    


/* Function:  access_use_flat_length_RandomDNADB(obj)
 *
 * Descrip:    Access member variable use_flat_length
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [RandomDNADB *]
 *
 * Return [SOFT ]  member variable use_flat_length [boolean]
 *
 */
boolean access_use_flat_length_RandomDNADB(RandomDNADB * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function use_flat_length for object RandomDNADB, got a NULL object");    
      return FALSE;  
      }  
    return obj->use_flat_length;     
}    


/* Function:  replace_length_RandomDNADB(obj,length)
 *
 * Descrip:    Replace member variable length
 *             For use principly by API functions
 *
 *
 * Arg:           obj [UNKN ] Object holding the variable [RandomDNADB *]
 * Arg:        length [OWNER] New value of the variable [int]
 *
 * Return [SOFT ]  member variable length [boolean]
 *
 */
boolean replace_length_RandomDNADB(RandomDNADB * obj,int length) 
{
    if( obj == NULL)     {  
      warn("In replacement function length for object RandomDNADB, got a NULL object");  
      return FALSE;  
      }  
    obj->length = length;    
    return TRUE; 
}    


/* Function:  access_length_RandomDNADB(obj)
 *
 * Descrip:    Access member variable length
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [RandomDNADB *]
 *
 * Return [SOFT ]  member variable length [int]
 *
 */
int access_length_RandomDNADB(RandomDNADB * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function length for object RandomDNADB, got a NULL object"); 
      return 0;  
      }  
    return obj->length;  
}    


/* Function:  replace_length_dist_RandomDNADB(obj,length_dist)
 *
 * Descrip:    Replace member variable length_dist
 *             For use principly by API functions
 *
 *
 * Arg:                obj [UNKN ] Object holding the variable [RandomDNADB *]
 * Arg:        length_dist [OWNER] New value of the variable [Histogram *]
 *
 * Return [SOFT ]  member variable length_dist [boolean]
 *
 */
boolean replace_length_dist_RandomDNADB(RandomDNADB * obj,Histogram * length_dist) 
{
    if( obj == NULL)     {  
      warn("In replacement function length_dist for object RandomDNADB, got a NULL object"); 
      return FALSE;  
      }  
    obj->length_dist = length_dist;  
    return TRUE; 
}    


/* Function:  access_length_dist_RandomDNADB(obj)
 *
 * Descrip:    Access member variable length_dist
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [RandomDNADB *]
 *
 * Return [SOFT ]  member variable length_dist [Histogram *]
 *
 */
Histogram * access_length_dist_RandomDNADB(RandomDNADB * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function length_dist for object RandomDNADB, got a NULL object");    
      return NULL;   
      }  
    return obj->length_dist;     
}    


/* Function:  replace_emission_RandomDNADB(obj,emission)
 *
 * Descrip:    Replace member variable emission
 *             For use principly by API functions
 *
 *
 * Arg:             obj [UNKN ] Object holding the variable [RandomDNADB *]
 * Arg:        emission [OWNER] New value of the variable [RandomModelDNA *]
 *
 * Return [SOFT ]  member variable emission [boolean]
 *
 */
boolean replace_emission_RandomDNADB(RandomDNADB * obj,RandomModelDNA * emission) 
{
    if( obj == NULL)     {  
      warn("In replacement function emission for object RandomDNADB, got a NULL object");    
      return FALSE;  
      }  
    obj->emission = emission;    
    return TRUE; 
}    


/* Function:  access_emission_RandomDNADB(obj)
 *
 * Descrip:    Access member variable emission
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [RandomDNADB *]
 *
 * Return [SOFT ]  member variable emission [RandomModelDNA *]
 *
 */
RandomModelDNA * access_emission_RandomDNADB(RandomDNADB * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function emission for object RandomDNADB, got a NULL object");   
      return NULL;   
      }  
    return obj->emission;    
}    


/* Function:  replace_num_RandomDNADB(obj,num)
 *
 * Descrip:    Replace member variable num
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [RandomDNADB *]
 * Arg:        num [OWNER] New value of the variable [int]
 *
 * Return [SOFT ]  member variable num [boolean]
 *
 */
boolean replace_num_RandomDNADB(RandomDNADB * obj,int num) 
{
    if( obj == NULL)     {  
      warn("In replacement function num for object RandomDNADB, got a NULL object"); 
      return FALSE;  
      }  
    obj->num = num;  
    return TRUE; 
}    


/* Function:  access_num_RandomDNADB(obj)
 *
 * Descrip:    Access member variable num
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [RandomDNADB *]
 *
 * Return [SOFT ]  member variable num [int]
 *
 */
int access_num_RandomDNADB(RandomDNADB * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function num for object RandomDNADB, got a NULL object");    
      return 0;  
      }  
    return obj->num;     
}    



#ifdef _cplusplus
}
#endif
