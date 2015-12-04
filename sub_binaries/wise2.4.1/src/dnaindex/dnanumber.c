#ifdef _cplusplus
extern "C" {
#endif
#include "dnanumber.h"

# line 23 "dnanumber.dy"
DnaNumberSequence * new_DnaNumberSequence(Sequence * seq,int nmer_size)
{
  int i;
  DnaNumberSequence * out;

  assert(seq);
  out = DnaNumberSequence_alloc();
  out->seq = malloc(sizeof(DnaNumber)*seq->len);
  
  for(i=0;i<seq->len;i++) 
    out->seq[i] = dna_number_from_string(seq->seq+i,nmer_size);
  
  out->orig = hard_link_Sequence(seq);
  out->len = seq->len-nmer_size-1;

  return out;
}

# line 41 "dnanumber.dy"
char first_char_from_dnanumber(int dnanumber,int nmer_size,int flipped)
{
  int base = 1;
  int basepair;
  int power;
  int i;

  if( flipped == 1 ) {
    /* first number is therefore lowest bit */
    basepair = dnanumber % 4;
    /*    fprintf(stderr,"Going to use %d as number\n",basepair);*/
    basepair = complement_base(basepair);

    return char_from_base(basepair);
  } else {
    for(i=0;i<nmer_size-1;i++) 
      base *= 4;
    
    basepair = (int) (dnanumber / base);
    
    return char_from_base(basepair);
  }

}

# line 66 "dnanumber.dy"
DnaNumber dna_number_from_string(char * str,int nmer_size)
{
  int i;
  int base = 1;
  DnaNumber out;
  int forward;
  int backward;

  out.flipped = 2;
  out.number  = 0;

  for(i=0;i<nmer_size-1;i++) 
    base *= 4;

  for(i=0;i<nmer_size;i++) {
    forward = base_from_char(str[i]);
    backward = complement_base(base_from_char(str[nmer_size-1-i]));
    
    if( forward == BASE_N || backward == BASE_N ) {
      return out; 
    }

    if( forward > backward ) {
      out.flipped = 0;
      break;
    } 
    if( backward > forward ) {
      out.flipped = 1;
      break;
    }
    
  }

  assert(out.flipped != 2);

  if( out.flipped == 0 ) {
    for(i=0;i<nmer_size;i++) {
      out.number += base * base_from_char(str[i]);
      base = base / 4;
    }
  } else {
    for(i=0;i<nmer_size;i++) {
      /*      fprintf(stderr,"For position %d, [%d], using %d [%c]as complemented base\n",i,base,complement_base(base_from_char(str[nmer_size-1-i])),str[nmer_size-1-i]); */

      out.number += base * complement_base(base_from_char(str[nmer_size-1-i]));
      base = base / 4;
    }
  }

  return out;

}

/* Function:  free_DnaNumber(dn)
 *
 * Descrip:    dummy free
 *
 *
 * Arg:        dn [UNKN ] Undocumented argument [DnaNumber *]
 *
 * Return [UNKN ]  Undocumented return value [DnaNumber *]
 *
 */
# line 122 "dnanumber.dy"
DnaNumber * free_DnaNumber(DnaNumber * dn)
{
  free(dn);
}


# line 119 "dnanumber.c"
/* Function:  hard_link_DnaNumberSequence(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [DnaNumberSequence *]
 *
 * Return [UNKN ]  Undocumented return value [DnaNumberSequence *]
 *
 */
DnaNumberSequence * hard_link_DnaNumberSequence(DnaNumberSequence * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a DnaNumberSequence object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  DnaNumberSequence_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DnaNumberSequence *]
 *
 */
DnaNumberSequence * DnaNumberSequence_alloc(void) 
{
    DnaNumberSequence * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(DnaNumberSequence *) ckalloc (sizeof(DnaNumberSequence))) == NULL)  {  
      warn("DnaNumberSequence_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->seq = NULL; 
    out->len = 0;    
    out->orig = NULL;    


    return out;  
}    


/* Function:  free_DnaNumberSequence(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [DnaNumberSequence *]
 *
 * Return [UNKN ]  Undocumented return value [DnaNumberSequence *]
 *
 */
DnaNumberSequence * free_DnaNumberSequence(DnaNumberSequence * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a DnaNumberSequence obj. Should be trappable"); 
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
    if( obj->seq != NULL)    
      free_DnaNumber(obj->seq);  
    if( obj->orig != NULL)   
      free_Sequence(obj->orig);  


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif
