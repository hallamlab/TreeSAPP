#ifdef _cplusplus
extern "C" {
#endif
#include "dna.h"



/* Function:  truncate_DNA(dna,start,stop)
 *
 * Descrip:    Truncates a DNA sequence. Basically uses
 *             the /trunc_Sequence function (of course!)
 *
 *             It does not alter dna, rather it returns a new
 *             sequence with that truncation
 *
 *
 * Arg:          dna [READ ] DNA that is truncated [DNA *]
 * Arg:        start [UNKN ] Undocumented argument [int]
 * Arg:         stop [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [DNA *]
 *
 */
# line 30 "dna.dy"
DNA * truncate_DNA(DNA * dna,int start,int stop)
{
  return DNA_from_Sequence(trunc_Sequence(dna->baseseq,start,stop));
}

/* Function:  read_fasta_file_DNA (filename)
 *
 * Descrip:    Reads a fasta file assumming that it is DNA. 
 *             Will complain if it is not, and return NULL.
 *
 *
 * Arg:        filename [UNKN ] filename to be opened and read [char *]
 *
 * Return [UNKN ]  Undocumented return value [DNA *]
 *
 */
# line 41 "dna.dy"
DNA * read_fasta_file_DNA (char * filename)
{
  Sequence * seq;

  seq = read_fasta_file_Sequence(filename);
  if( seq == NULL ) {
    return NULL;
  }

  return DNA_from_Sequence(seq);
}


/* Function:  read_fasta_DNA (ifp)
 *
 * Descrip:    Reads a fasta file assumming that it is DNA. 
 *             Will complain if it is not, and return NULL.
 *
 *
 * Arg:        ifp [UNKN ] file point to be read from [FILE *]
 *
 * Return [UNKN ]  Undocumented return value [DNA *]
 *
 */
# line 60 "dna.dy"
DNA * read_fasta_DNA (FILE * ifp)
{
  Sequence * seq;

  seq = read_fasta_Sequence(ifp);
  if( seq == NULL ) {
    return NULL;
  }

  return DNA_from_Sequence(seq);
}

/* Function:  read_efetch_DNA(estr)
 *
 * Descrip:    Reads a efetch specified query
 *             Uses, of course /read_efetch_Sequence
 *
 *
 * Arg:        estr [READ ] efetch string which is read [char *]
 *
 * Return [UNKN ]  Undocumented return value [DNA *]
 *
 */
# line 78 "dna.dy"
DNA * read_efetch_DNA(char * estr)
{
  return DNA_from_Sequence(read_efetch_Sequence(estr));
}

/* Function:  read_SRS_DNA(srsquery)
 *
 * Descrip:    Reads a SRS sequence using srs4 syntax.
 *             Uses, of course, /read_SRS_Sequence
 *
 *
 *
 * Arg:        srsquery [READ ] string query representing SRS name [char *]
 *
 * Return [UNKN ]  Undocumented return value [DNA *]
 *
 */
# line 90 "dna.dy"
DNA * read_SRS_DNA(char * srsquery)
{
  return DNA_from_Sequence(read_SRS_Sequence(srsquery));
}


/* Function:  DNA_name (dna)
 *
 * Descrip:    Returns the name of the DNA
 *
 *
 * Arg:        dna [UNKN ] Undocumented argument [DNA *]
 *
 * Return [UNKN ]  Undocumented return value [char *]
 *
 */
# line 100 "dna.dy"
char * DNA_name (DNA * dna)
{
  return dna->baseseq->name;
}

/* Function:  DNA_length (dna)
 *
 * Descrip:    Returns the length of the DNA
 *
 *
 * Arg:        dna [UNKN ] Undocumented argument [DNA *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
# line 109 "dna.dy"
int DNA_length (DNA * dna)
{
  return dna->baseseq->len;
}

/* Function:  DNA_seqchar(dna,pos)
 *
 * Descrip:    Returns sequence character at this position.
 *
 *
 * Arg:        dna [UNKN ] DNA [DNA *]
 * Arg:        pos [UNKN ] position in DNA to get char [int]
 *
 * Return [UNKN ]  Undocumented return value [char]
 *
 */
# line 120 "dna.dy"
char DNA_seqchar(DNA * dna,int pos)
{
  return dna->baseseq->seq[pos];
}


/* Function:  DNA_from_Sequence(seq)
 *
 * Descrip:    makes a new DNA from a Sequence. It 
 *             owns the Sequence memory, ie will attempt a /free_Sequence
 *             on the structure when /free_DNA is called
 *
 *             If you want to give this DNA this Sequence and
 *             forget about it, then just hand it this sequence and set
 *             seq to NULL (no need to free it). If you intend to use 
 *             the sequence object elsewhere outside of the DNA datastructure
 *             then use DNA_from_Sequence(/hard_link_Sequence(seq))
 *
 *
 *
 * Arg:        seq [UNKN ] Sequence to make DNA from [Sequence *]
 *
 * Return [UNKN ]  Undocumented return value [DNA *]
 *
 */
# line 140 "dna.dy"
DNA * DNA_from_Sequence(Sequence * seq)
{
  DNA * out;


  if( is_dna_Sequence(seq) == FALSE ) {
    warn("Trying to make a DNA sequence from a non DNA base sequence [%s].",seq->name);
    return NULL;
  }

  out = DNA_alloc();

  out->baseseq = seq;

  return out;
}




# line 186 "dna.c"
/* Function:  hard_link_DNA(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [DNA *]
 *
 * Return [UNKN ]  Undocumented return value [DNA *]
 *
 */
DNA * hard_link_DNA(DNA * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a DNA object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  DNA_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DNA *]
 *
 */
DNA * DNA_alloc(void) 
{
    DNA * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(DNA *) ckalloc (sizeof(DNA))) == NULL)  {  
      warn("DNA_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->baseseq = NULL; 


    return out;  
}    


/* Function:  free_DNA(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [DNA *]
 *
 * Return [UNKN ]  Undocumented return value [DNA *]
 *
 */
DNA * free_DNA(DNA * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a DNA obj. Should be trappable");   
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



#ifdef _cplusplus
}
#endif
