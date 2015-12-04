#ifdef _cplusplus
extern "C" {
#endif
#include "cdna.h"



/* Function:  truncate_cDNA(cdna,start,stop)
 *
 * Descrip:    Truncates a cDNA sequence. Basically uses
 *             the /magic_trunc_Sequence function (of course!)
 *
 *             It does not alter cdna, rather it returns a new
 *             sequence with that truncation
 *
 *
 * Arg:         cdna [READ ] cDNA that is truncated [cDNA *]
 * Arg:        start [UNKN ] Undocumented argument [int]
 * Arg:         stop [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [cDNA *]
 *
 */
# line 45 "cdna.dy"
cDNA * truncate_cDNA(cDNA * cdna,int start,int stop)
{
  return cDNA_from_Sequence(magic_trunc_Sequence(cdna->baseseq,start,stop));
}

/* Function:  read_fasta_file_cDNA(filename)
 *
 * Descrip:    Reads a fasta file assumming that it is cDNA. 
 *             Will complain if it is not, and return NULL.
 *
 *
 * Arg:        filename [UNKN ] filename to be opened and read [char *]
 *
 * Return [UNKN ]  Undocumented return value [cDNA *]
 *
 */
# line 56 "cdna.dy"
cDNA * read_fasta_file_cDNA(char * filename)
{
  Sequence * seq;

  seq = read_fasta_file_Sequence(filename);
  if( seq == NULL ) {
    return NULL;
  }

  return cDNA_from_Sequence(seq);
}


/* Function:  read_fasta_cDNA(ifp)
 *
 * Descrip:    Reads a fasta file assumming that it is cDNA. 
 *             Will complain if it is not, and return NULL.
 *
 *
 * Arg:        ifp [UNKN ] file point to be read from [FILE *]
 *
 * Return [UNKN ]  Undocumented return value [cDNA *]
 *
 */
# line 75 "cdna.dy"
cDNA * read_fasta_cDNA(FILE * ifp)
{
  Sequence * seq;

  seq = read_fasta_Sequence(ifp);
  if( seq == NULL ) {
    return NULL;
  }

  return cDNA_from_Sequence(seq);
}

/* Function:  read_efetch_cDNA(estr)
 *
 * Descrip:    Reads a efetch specified query
 *             Uses, of course /read_efetch_Sequence
 *
 *
 * Arg:        estr [READ ] efetch string which is read [char *]
 *
 * Return [UNKN ]  Undocumented return value [cDNA *]
 *
 */
# line 93 "cdna.dy"
cDNA * read_efetch_cDNA(char * estr)
{
  return cDNA_from_Sequence(read_efetch_Sequence(estr));
}

/* Function:  read_SRS_cDNA(srsquery)
 *
 * Descrip:    Reads a SRS sequence using srs4 syntax.
 *             Uses, of course, /read_SRS_Sequence
 *
 *
 *
 * Arg:        srsquery [READ ] string query representing SRS name [char *]
 *
 * Return [UNKN ]  Undocumented return value [cDNA *]
 *
 */
# line 105 "cdna.dy"
cDNA * read_SRS_cDNA(char * srsquery)
{
  return cDNA_from_Sequence(read_SRS_Sequence(srsquery));
}


/* Function:  cDNA_name(cdna)
 *
 * Descrip:    Returns the name of the cDNA
 *
 *
 * Arg:        cdna [UNKN ] Undocumented argument [cDNA *]
 *
 * Return [UNKN ]  Undocumented return value [char *]
 *
 */
# line 115 "cdna.dy"
char * cDNA_name(cDNA * cdna)
{
  return cdna->baseseq->name;
}

/* Function:  cDNA_length(cdna)
 *
 * Descrip:    Returns the length of the cDNA
 *
 *
 * Arg:        cdna [UNKN ] Undocumented argument [cDNA *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
# line 124 "cdna.dy"
int cDNA_length(cDNA * cdna)
{
  return cdna->baseseq->len;
}

/* Function:  cDNA_seqchar(cdna,pos)
 *
 * Descrip:    Returns sequence character at this position.
 *
 *
 * Arg:        cdna [UNKN ] cDNA [cDNA *]
 * Arg:         pos [UNKN ] position in cDNA to get char [int]
 *
 * Return [UNKN ]  Undocumented return value [char]
 *
 */
# line 135 "cdna.dy"
char cDNA_seqchar(cDNA * cdna,int pos)
{
  return cdna->baseseq->seq[pos];
}


/* Function:  cDNA_from_Sequence(seq)
 *
 * Descrip:    makes a new cDNA from a Sequence. It 
 *             owns the Sequence memory, ie will attempt a /free_Sequence
 *             on the structure when /free_cDNA is called
 *
 *             If you want to give this cDNA this Sequence and
 *             forget about it, then just hand it this sequence and set
 *             seq to NULL (no need to free it). If you intend to use 
 *             the sequence object elsewhere outside of the cDNA datastructure
 *             then use cDNA_from_Sequence(/hard_link_Sequence(seq))
 *
 *
 *
 * Arg:        seq [OWNER] Sequence to make cDNA from [Sequence *]
 *
 * Return [UNKN ]  Undocumented return value [cDNA *]
 *
 */
# line 155 "cdna.dy"
cDNA * cDNA_from_Sequence(Sequence * seq)
{
  cDNA * out;
  int conv;

  if( seq == NULL ) {
    warn("Trying to make a cdna sequence from a NULL baseseq.");
    return NULL;
  }

  if( is_dna_Sequence(seq) == FALSE ) {
    warn("Trying to make a cDNA sequence from a non cDNA base sequence [%s].",seq->name);
    return NULL;
  }

  uppercase_Sequence(seq);

  force_to_dna_Sequence(seq,1.0,&conv);
 
  if( conv != 0 ) {
    log_full_error(INFO,0,"In making %s a cdna sequence, converted %d bases (%2.1f%%) to N's from non ATGCN",seq->name,conv,(double)conv*100/(double)seq->len);
  }

  out = cDNA_alloc();

  out->baseseq = seq;

  return out;
}




# line 199 "cdna.c"
/* Function:  hard_link_cDNA(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [cDNA *]
 *
 * Return [UNKN ]  Undocumented return value [cDNA *]
 *
 */
cDNA * hard_link_cDNA(cDNA * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a cDNA object: passed a NULL object");    
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  cDNA_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [cDNA *]
 *
 */
cDNA * cDNA_alloc(void) 
{
    cDNA * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(cDNA *) ckalloc (sizeof(cDNA))) == NULL)    {  
      warn("cDNA_alloc failed ");    
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->baseseq = NULL; 


    return out;  
}    


/* Function:  free_cDNA(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [cDNA *]
 *
 * Return [UNKN ]  Undocumented return value [cDNA *]
 *
 */
cDNA * free_cDNA(cDNA * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a cDNA obj. Should be trappable");  
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
    if( obj->baseseq != NULL)    
      free_Sequence(obj->baseseq);   


    ckfree(obj); 
    return NULL; 
}    


/* Function:  replace_baseseq_cDNA(obj,baseseq)
 *
 * Descrip:    Replace member variable baseseq
 *             For use principly by API functions
 *
 *
 * Arg:            obj [UNKN ] Object holding the variable [cDNA *]
 * Arg:        baseseq [OWNER] New value of the variable [Sequence *]
 *
 * Return [SOFT ]  member variable baseseq [boolean]
 *
 */
boolean replace_baseseq_cDNA(cDNA * obj,Sequence * baseseq) 
{
    if( obj == NULL)     {  
      warn("In replacement function baseseq for object cDNA, got a NULL object");    
      return FALSE;  
      }  
    obj->baseseq = baseseq;  
    return TRUE; 
}    


/* Function:  access_baseseq_cDNA(obj)
 *
 * Descrip:    Access member variable baseseq
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [cDNA *]
 *
 * Return [SOFT ]  member variable baseseq [Sequence *]
 *
 */
Sequence * access_baseseq_cDNA(cDNA * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function baseseq for object cDNA, got a NULL object");   
      return NULL;   
      }  
    return obj->baseseq;     
}    



#ifdef _cplusplus
}
#endif
