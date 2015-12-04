#ifdef _cplusplus
extern "C" {
#endif
#include "shotgun.h"




/* Function:  generate_shotgun_reads(shot,input,ofp)
 *
 * Descrip:    Generates a file of Shotgun reads from a particular
 *             sequence randomly
 *
 *
 * Arg:         shot [UNKN ] Undocumented argument [ShotgunPara *]
 * Arg:        input [UNKN ] Undocumented argument [Sequence *]
 * Arg:          ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 27 "shotgun.dy"
void generate_shotgun_reads(ShotgunPara * shot,Sequence * input,FILE * ofp)
{
  int i;
  int j;
  int pos;
  
  Sequence * rev;

  Sequence * out;
  char * seqstr;
  char buffer[MAXLINE];


  init_random();

  seqstr = calloc(shot->read_length+1,sizeof(char));
  rev = reverse_complement_Sequence(input);
 
  for(i=0;i<shot->number;i++) {
    pos = random_integer(input->len-shot->insert_size);
    fprintf(stderr,"position at %d\n",pos);

    for(j=0;j<shot->read_length;j++) {
      seqstr[j] = input->seq[pos+j];
    }
    
    seqstr[j] = '\0';
    sprintf(buffer,"Read.%d.f",i);
    out = new_Sequence_from_strings(buffer,seqstr);
    write_fasta_Sequence(out,ofp);

    if( shot->forward_only == 1 ) {
      continue;
    }

    for(j=0;j<shot->read_length;j++) {
      seqstr[j] = rev->seq[rev->len-(pos+shot->insert_size+j)];
    }
    
    seqstr[j] = '\0';
    sprintf(buffer,"Read.%d.r",i);
    out = new_Sequence_from_strings(buffer,seqstr);
    write_fasta_Sequence(out,ofp);
  }

  free_Sequence(rev);
  free(seqstr);

}








# line 75 "shotgun.c"
/* Function:  hard_link_ShotgunPara(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [ShotgunPara *]
 *
 * Return [UNKN ]  Undocumented return value [ShotgunPara *]
 *
 */
ShotgunPara * hard_link_ShotgunPara(ShotgunPara * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a ShotgunPara object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  ShotgunPara_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ShotgunPara *]
 *
 */
ShotgunPara * ShotgunPara_alloc(void) 
{
    ShotgunPara * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(ShotgunPara *) ckalloc (sizeof(ShotgunPara))) == NULL)  {  
      warn("ShotgunPara_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->read_length = 0;    
    out->insert_size = 0;    
    out->number = 0; 
    out->forward_only = 0;   


    return out;  
}    


/* Function:  free_ShotgunPara(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [ShotgunPara *]
 *
 * Return [UNKN ]  Undocumented return value [ShotgunPara *]
 *
 */
ShotgunPara * free_ShotgunPara(ShotgunPara * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a ShotgunPara obj. Should be trappable");   
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
