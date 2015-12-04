#ifdef _cplusplus
extern "C" {
#endif
#include "protein.h"



/* Function:  truncate_Protein(pro,start,stop)
 *
 * Descrip:    Truncates a protein sequence. Basically uses
 *             the /trunc_Sequence function (of course!)
 *
 *             It does not alter pro, rather it returns a new
 *             sequence with that truncation
 *
 *
 * Arg:          pro [READ ] Protein that is truncated [Protein *]
 * Arg:        start [UNKN ] Undocumented argument [int]
 * Arg:         stop [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [Protein *]
 *
 */
# line 45 "protein.dy"
Protein * truncate_Protein(Protein * pro,int start,int stop)
{
  return Protein_from_Sequence(trunc_Sequence(pro->baseseq,start,stop));
}

/* Function:  read_fasta_file_Protein (filename)
 *
 * Descrip:    Reads a fasta file assumming that it is protein. 
 *             Will complain if it is not, and return NULL.
 *
 *
 * Arg:        filename [UNKN ] filename to be opened and read [char *]
 *
 * Return [UNKN ]  Undocumented return value [Protein *]
 *
 */
# line 56 "protein.dy"
Protein * read_fasta_file_Protein (char * filename)
{
  Protein * out;
  Sequence * seq;

  seq = read_fasta_file_Sequence(filename);
  if( seq == NULL ) {
    return NULL;
  }
  
  out =  Protein_from_Sequence(seq);

  return out;
}


/* Function:  read_fasta_Protein (ifp)
 *
 * Descrip:    Reads a fasta file assumming that it is protein. 
 *             Will complain if it is not, and return NULL.
 *
 *
 * Arg:        ifp [UNKN ] file point to be read from [FILE *]
 *
 * Return [UNKN ]  Undocumented return value [Protein *]
 *
 */
# line 78 "protein.dy"
Protein * read_fasta_Protein (FILE * ifp)
{
  Sequence * seq;

  seq = read_fasta_Sequence(ifp);
  if( seq == NULL ) {
    return NULL;
  }

  return Protein_from_Sequence(seq);
}

/* Function:  read_efetch_Protein(estr)
 *
 * Descrip:    Reads a efetch specified query
 *             Uses, of course /read_efetch_Sequence
 *
 *
 * Arg:        estr [READ ] efetch string which is read [char *]
 *
 * Return [UNKN ]  Undocumented return value [Protein *]
 *
 */
# line 96 "protein.dy"
Protein * read_efetch_Protein(char * estr)
{
  return Protein_from_Sequence(read_efetch_Sequence(estr));
}

/* Function:  read_SRS_Protein(srsquery)
 *
 * Descrip:    Reads a SRS sequence using srs4 syntax.
 *             Uses, of course, /read_SRS_Sequence
 *
 *
 *
 * Arg:        srsquery [READ ] string query representing SRS name [char *]
 *
 * Return [UNKN ]  Undocumented return value [Protein *]
 *
 */
# line 108 "protein.dy"
Protein * read_SRS_Protein(char * srsquery)
{
  return Protein_from_Sequence(read_SRS_Sequence(srsquery));
}


/* Function:  Protein_name (pr)
 *
 * Descrip:    Returns the name of the protein
 *
 *
 * Arg:        pr [UNKN ] Undocumented argument [Protein *]
 *
 * Return [UNKN ]  Undocumented return value [char *]
 *
 */
# line 118 "protein.dy"
char * Protein_name (Protein * pr)
{
  return pr->baseseq->name;
}

/* Function:  Protein_length (pr)
 *
 * Descrip:    Returns the length of the protein
 *
 *
 * Arg:        pr [UNKN ] Undocumented argument [Protein *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
# line 127 "protein.dy"
int Protein_length (Protein * pr)
{
  return pr->baseseq->len;
}

/* Function:  Protein_seqchar(pr,pos)
 *
 * Descrip:    Returns sequence character at this position.
 *
 *
 * Arg:         pr [UNKN ] Protein [Protein *]
 * Arg:        pos [UNKN ] position in protein to get char [int]
 *
 * Return [UNKN ]  Undocumented return value [char]
 *
 */
# line 138 "protein.dy"
char Protein_seqchar(Protein * pr,int pos)
{
  return pr->baseseq->seq[pos];
}


/* Function:  Protein_from_Sequence(seq)
 *
 * Descrip:    makes a new protein from a Sequence. It 
 *             owns the Sequence memory, ie will attempt a /free_Sequence
 *             on the structure when /free_Protein is called
 *
 *             If you want to give this protein this Sequence and
 *             forget about it, then just hand it this sequence and set
 *             seq to NULL (no need to free it). If you intend to use 
 *             the sequecne object elsewhere outside of the Protein datastructure
 *             then use Protein_from_Sequence(/hard_link_Sequence(seq))
 *
 *
 *
 * Arg:        seq [OWNER] Sequence to make protein from [Sequence *]
 *
 * Return [UNKN ]  Undocumented return value [Protein *]
 *
 */
# line 158 "protein.dy"
Protein * Protein_from_Sequence(Sequence * seq)
{
  Protein * out;


  if( seq == NULL ) {
    warn("Attempting to make a protein from a NULL sequence");
    return NULL;
  }

  if( is_protein_Sequence(seq) == FALSE ) {
    warn("Trying to make a protein sequence from a non protein base sequence [%s].",seq->name);
    return NULL;
  }

  out = Protein_alloc();
  
  out->baseseq = seq;

  return out;
}




# line 194 "protein.c"
/* Function:  hard_link_Protein(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [Protein *]
 *
 * Return [UNKN ]  Undocumented return value [Protein *]
 *
 */
Protein * hard_link_Protein(Protein * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a Protein object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  Protein_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Protein *]
 *
 */
Protein * Protein_alloc(void) 
{
    Protein * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(Protein *) ckalloc (sizeof(Protein))) == NULL)  {  
      warn("Protein_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->baseseq = NULL; 


    return out;  
}    


/* Function:  free_Protein(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [Protein *]
 *
 * Return [UNKN ]  Undocumented return value [Protein *]
 *
 */
Protein * free_Protein(Protein * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a Protein obj. Should be trappable");   
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


/* Function:  replace_baseseq_Protein(obj,baseseq)
 *
 * Descrip:    Replace member variable baseseq
 *             For use principly by API functions
 *
 *
 * Arg:            obj [UNKN ] Object holding the variable [Protein *]
 * Arg:        baseseq [OWNER] New value of the variable [Sequence *]
 *
 * Return [SOFT ]  member variable baseseq [boolean]
 *
 */
boolean replace_baseseq_Protein(Protein * obj,Sequence * baseseq) 
{
    if( obj == NULL)     {  
      warn("In replacement function baseseq for object Protein, got a NULL object"); 
      return FALSE;  
      }  
    obj->baseseq = baseseq;  
    return TRUE; 
}    


/* Function:  access_baseseq_Protein(obj)
 *
 * Descrip:    Access member variable baseseq
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [Protein *]
 *
 * Return [SOFT ]  member variable baseseq [Sequence *]
 *
 */
Sequence * access_baseseq_Protein(Protein * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function baseseq for object Protein, got a NULL object");    
      return NULL;   
      }  
    return obj->baseseq;     
}    



#ifdef _cplusplus
}
#endif
