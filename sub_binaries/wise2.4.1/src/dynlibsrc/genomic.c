#ifdef _cplusplus
extern "C" {
#endif
#include "genomic.h"



/* Function:  truncate_Genomic(gen,start,stop)
 *
 * Descrip:    Truncates a Genomic sequence. Basically uses
 *             the /magic_trunc_Sequence function (of course!)
 *
 *             It does not alter gen, rather it returns a new
 *             sequence with that truncation
 *
 *             Handles repeat information correctly.
 *
 *
 * Arg:          gen [READ ] Genomic that is truncated [Genomic *]
 * Arg:        start [UNKN ] Undocumented argument [int]
 * Arg:         stop [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [Genomic *]
 *
 */
# line 60 "genomic.dy"
Genomic * truncate_Genomic(Genomic * gen,int start,int stop)
{
  Genomic * out;
  GenomicRepeat * tempr;
  int tst,tend;
  int i;


  out = Genomic_from_Sequence(magic_trunc_Sequence(gen->baseseq,start,stop));
  assert(out);

  for(i=0;i<gen->len;i++) {
    if( start < stop ) {
      /* then we have to figure out if this repeat overlaps, and if so, put it in */
      if( gen->repeat[i]->start > stop || gen->repeat[i]->end < start ) {
	continue;
      }

      /* overlaps in some way */
      tst = gen->repeat[i]->start - start;
      tend = gen->repeat[i]->end - start;

      if( tst < 0 ) {
	tst =0;
      }
      if( tend > (stop-start) ) {
	tend = (stop-start);
      }
    } else {
      /* then we have to figure out if this repeat overlaps, and if so, put it in */
      if( gen->repeat[i]->start > start || gen->repeat[i]->end < stop ) {
	continue;
      }

      /* overlaps in some way */
      tend =  start - gen->repeat[i]->start;
      tst  = start - gen->repeat[i]->end;

      if( tst < 0 ) {
	tst =0;
      }
      if( tend > (start-stop) ) {
	tend = (start-stop);
      }
    }

    tempr = GenomicRepeat_alloc();
    tempr->start = tst;
    tempr->end = tend;
    add_Genomic(out,tempr);
    
  }

  return out;
    
}


/* Function:  reverse_complement_Genomic(gen)
 *
 * Descrip:    Reverse Complements s Genomic sequence. 
 *
 *             Handles repeat information correctly
 *
 *
 * Arg:        gen [READ ] Genomic that is revomcp [Genomic *]
 *
 * Return [UNKN ]  Undocumented return value [Genomic *]
 *
 */
# line 125 "genomic.dy"
Genomic * reverse_complement_Genomic(Genomic * gen)
{
  return truncate_Genomic(gen,gen->baseseq->len,0);
}


/* Function:  read_fasta_file_Genomic(filename,length_of_N)
 *
 * Descrip:    Reads a fasta file assumming that it is Genomic. 
 *             Will complain if it is not, and return NULL.
 *
 *
 * Arg:           filename [UNKN ] filename to be opened and read [char *]
 * Arg:        length_of_N [UNKN ] length of N to be considered repeat. -1 means none [int]
 *
 * Return [UNKN ]  Undocumented return value [Genomic *]
 *
 */
# line 138 "genomic.dy"
Genomic * read_fasta_file_Genomic(char * filename,int length_of_N)
{
  Sequence * seq;

  seq = read_fasta_file_Sequence(filename);
  if( seq == NULL ) {
    warn("Unable to read sequence from %s [%d], so no genomic",filename,seq);
    return NULL;
  }


  if( length_of_N < 0 ) 
    return Genomic_from_Sequence(seq);
  else 
    return Genomic_from_Sequence_Nheuristic(seq,length_of_N);
 
}


/* Function:  read_fasta_Genomic(ifp,length_of_N)
 *
 * Descrip:    Reads a fasta file assumming that it is Genomic. 
 *             Will complain if it is not, and return NULL.
 *
 *
 * Arg:                ifp [UNKN ] file point to be read from [FILE *]
 * Arg:        length_of_N [UNKN ] length of N to be considered repeat. -1 means none [int]
 *
 * Return [UNKN ]  Undocumented return value [Genomic *]
 *
 */
# line 164 "genomic.dy"
Genomic * read_fasta_Genomic(FILE * ifp,int length_of_N)
{
  Sequence * seq;

  seq = read_fasta_Sequence(ifp);
  if( seq == NULL ) {
    return NULL;
  }

  if( length_of_N < 0 ) 
    return Genomic_from_Sequence(seq);
  else 
    return Genomic_from_Sequence_Nheuristic(seq,length_of_N);

}


/* Function:  Genomic_name(gen)
 *
 * Descrip:    Returns the name of the Genomic
 *
 *
 * Arg:        gen [UNKN ] Undocumented argument [Genomic *]
 *
 * Return [UNKN ]  Undocumented return value [char *]
 *
 */
# line 185 "genomic.dy"
char * Genomic_name(Genomic * gen)
{
  return gen->baseseq->name;
}

/* Function:  Genomic_length(gen)
 *
 * Descrip:    Returns the length of the Genomic
 *
 *
 * Arg:        gen [UNKN ] Undocumented argument [Genomic *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
# line 194 "genomic.dy"
int Genomic_length(Genomic * gen)
{
  return gen->baseseq->len;
}

/* Function:  Genomic_seqchar(gen,pos)
 *
 * Descrip:    Returns sequence character at this position.
 *
 *
 * Arg:        gen [UNKN ] Genomic [Genomic *]
 * Arg:        pos [UNKN ] position in Genomic to get char [int]
 *
 * Return [UNKN ]  Undocumented return value [char]
 *
 */
# line 205 "genomic.dy"
char Genomic_seqchar(Genomic * gen,int pos)
{
  return gen->baseseq->seq[pos];
}


/* Function:  Genomic_from_Sequence_Nheuristic(seq,length_of_N)
 *
 * Descrip:    makes a new genomic from a Sequence, but
 *             assummes that all the N runs greater than
 *             a certain level are actually repeats.
 *
 *
 * Arg:                seq [UNKN ] Undocumented argument [Sequence *]
 * Arg:        length_of_N [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [Genomic *]
 *
 */
# line 216 "genomic.dy"
Genomic * Genomic_from_Sequence_Nheuristic(Sequence * seq,int length_of_N)
{
  Genomic * out;
  GenomicRepeat * temp;
  char * run, *run2;
  
  out = Genomic_from_Sequence(seq);
  if( out == NULL ) {
    warn("Could not make a new Genomic Sequence from Sequence in Nheuristic");
    return NULL;
  }

  for(run=strchr(seq->seq,'N');run != NULL;) {
    for(run2 = run; *run2 && *run2 == 'N';run2++)
      ;
    if( run2 - run > length_of_N) {
      temp = GenomicRepeat_alloc();
      add_Genomic(out,temp);
      temp->start = run - seq->seq;
      temp->end   = run2 - seq->seq;
    } 
    run = strchr(run2,'N');
  }

  return out;
}
  

/* Function:  Genomic_from_Sequence(seq)
 *
 * Descrip:    makes a new genomic from a Sequence. It 
 *             owns the Sequence memory, ie will attempt a /free_Sequence
 *             on the structure when /free_Genomic is called
 *
 *             If you want to give this genomic this Sequence and
 *             forget about it, then just hand it this sequence and set
 *             seq to NULL (no need to free it). If you intend to use 
 *             the sequence object elsewhere outside of the Genomic datastructure
 *             then use Genomic_from_Sequence(/hard_link_Sequence(seq))
 *
 *             This is part of a strict typing system, and therefore
 *             is going to convert all non ATGCNs to Ns. You will lose
 *             information here.
 *
 *
 * Arg:        seq [OWNER] Sequence to make genomic from [Sequence *]
 *
 * Return [UNKN ]  Undocumented return value [Genomic *]
 *
 */
# line 261 "genomic.dy"
Genomic * Genomic_from_Sequence(Sequence * seq)
{
  Genomic * out;
  int conv;

  if( seq == NULL ) {
    warn("Cannot make a genomic sequence from a NULL sequence");
    return NULL;
  }


  if( is_dna_Sequence(seq) == FALSE ) {
    warn("Trying to make a genomic sequence from a non genomic base sequence [%s] Type is %d.",seq->name,seq->type);
    return NULL;
  }

  uppercase_Sequence(seq);

  force_to_dna_Sequence(seq,1.0,&conv);
 
  if( conv != 0 ) {
    log_full_error(INFO,0,"In making %s a genomic sequence, converted %d bases (%2.1f%%) to N's from non ATGCN",seq->name,conv,(double)conv*100/(double)seq->len);
  }

  out = Genomic_alloc_std();

  out->baseseq = seq;

  return out;
}


/* Function:  show_Genomic(gen,ofp)
 *
 * Descrip:    For debugging
 *
 *
 * Arg:        gen [UNKN ] Undocumented argument [Genomic *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 297 "genomic.dy"
void show_Genomic(Genomic * gen,FILE * ofp)
{
  int i;
  assert(gen);
  write_fasta_Sequence(gen->baseseq,ofp);
  fprintf(ofp,"%d repeats\n",gen->len);
  for(i=0;i<gen->len;i++) {
    if( gen->baseseq->offset < gen->baseseq->end ) 
      fprintf(ofp,"Repeat from %d - %d  Biologically: %d - %d \n",gen->repeat[i]->start,gen->repeat[i]->end,gen->baseseq->offset+gen->repeat[i]->start,gen->baseseq->offset -1 + gen->repeat[i]->end);
    else 
      fprintf(ofp,"Repeat from %d - %d  Biologically: %d - %d \n",gen->repeat[i]->start,gen->repeat[i]->end,gen->baseseq->offset - gen->repeat[i]->start,gen->baseseq->offset - gen->repeat[i]->end +1);
  }
  
  
}

# line 318 "genomic.c"
/* Function:  hard_link_GenomicRepeat(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GenomicRepeat *]
 *
 * Return [UNKN ]  Undocumented return value [GenomicRepeat *]
 *
 */
GenomicRepeat * hard_link_GenomicRepeat(GenomicRepeat * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a GenomicRepeat object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  GenomicRepeat_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GenomicRepeat *]
 *
 */
GenomicRepeat * GenomicRepeat_alloc(void) 
{
    GenomicRepeat * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(GenomicRepeat *) ckalloc (sizeof(GenomicRepeat))) == NULL)  {  
      warn("GenomicRepeat_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->start = 0;  
    out->end = 0;    
    out->type = NULL;    


    return out;  
}    


/* Function:  free_GenomicRepeat(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GenomicRepeat *]
 *
 * Return [UNKN ]  Undocumented return value [GenomicRepeat *]
 *
 */
GenomicRepeat * free_GenomicRepeat(GenomicRepeat * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a GenomicRepeat obj. Should be trappable"); 
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
    if( obj->type != NULL)   
      ckfree(obj->type);     


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_Genomic(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_Genomic
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [GenomicRepeat **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_Genomic(GenomicRepeat ** list,int i,int j)  
{
    GenomicRepeat * temp;    
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_Genomic(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_Genomic which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [GenomicRepeat **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_Genomic(GenomicRepeat ** list,int left,int right,int (*comp)(GenomicRepeat * ,GenomicRepeat * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_Genomic(list,left,(left+right)/2);  
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_Genomic (list,++last,i);    
      }  
    swap_Genomic (list,left,last);   
    qsort_Genomic(list,left,last-1,comp);    
    qsort_Genomic(list,last+1,right,comp);   
}    


/* Function:  sort_Genomic(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_Genomic
 *
 *
 * Arg:         obj [UNKN ] Object containing list [Genomic *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_Genomic(Genomic * obj,int (*comp)(GenomicRepeat *, GenomicRepeat *)) 
{
    qsort_Genomic(obj->repeat,0,obj->len-1,comp);    
    return;  
}    


/* Function:  expand_Genomic(obj,len)
 *
 * Descrip:    Really an internal function for add_Genomic
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [Genomic *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_Genomic(Genomic * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_Genomic called with no need");    
      return TRUE;   
      }  


    if( (obj->repeat = (GenomicRepeat ** ) ckrealloc (obj->repeat,sizeof(GenomicRepeat *)*len)) == NULL)     {  
      warn("ckrealloc failed for expand_Genomic, returning FALSE");  
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_Genomic(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [Genomic *]
 * Arg:        add [OWNER] Object to add to the list [GenomicRepeat *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_Genomic(Genomic * obj,GenomicRepeat * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_Genomic(obj,obj->len + GenomicLISTLENGTH) == FALSE) 
        return FALSE;    
      }  


    obj->repeat[obj->len++]=add; 
    return TRUE; 
}    


/* Function:  flush_Genomic(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [Genomic *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_Genomic(Genomic * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->repeat[i] != NULL)    {  
        free_GenomicRepeat(obj->repeat[i]);  
        obj->repeat[i] = NULL;   
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  Genomic_alloc_std(void)
 *
 * Descrip:    Equivalent to Genomic_alloc_len(GenomicLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Genomic *]
 *
 */
Genomic * Genomic_alloc_std(void) 
{
    return Genomic_alloc_len(GenomicLISTLENGTH); 
}    


/* Function:  Genomic_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [Genomic *]
 *
 */
Genomic * Genomic_alloc_len(int len) 
{
    Genomic * out;  /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = Genomic_alloc()) == NULL)  
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->repeat = (GenomicRepeat ** ) ckcalloc (len,sizeof(GenomicRepeat *))) == NULL)   {  
      warn("Warning, ckcalloc failed in Genomic_alloc_len"); 
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_Genomic(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [Genomic *]
 *
 * Return [UNKN ]  Undocumented return value [Genomic *]
 *
 */
Genomic * hard_link_Genomic(Genomic * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a Genomic object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  Genomic_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Genomic *]
 *
 */
Genomic * Genomic_alloc(void) 
{
    Genomic * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(Genomic *) ckalloc (sizeof(Genomic))) == NULL)  {  
      warn("Genomic_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->baseseq = NULL; 
    out->repeat = NULL;  
    out->len = out->maxlen = 0;  


    return out;  
}    


/* Function:  free_Genomic(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [Genomic *]
 *
 * Return [UNKN ]  Undocumented return value [Genomic *]
 *
 */
Genomic * free_Genomic(Genomic * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a Genomic obj. Should be trappable");   
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
    if( obj->repeat != NULL) {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->repeat[i] != NULL)  
          free_GenomicRepeat(obj->repeat[i]);    
        }  
      ckfree(obj->repeat);   
      }  


    ckfree(obj); 
    return NULL; 
}    


/* Function:  replace_baseseq_Genomic(obj,baseseq)
 *
 * Descrip:    Replace member variable baseseq
 *             For use principly by API functions
 *
 *
 * Arg:            obj [UNKN ] Object holding the variable [Genomic *]
 * Arg:        baseseq [OWNER] New value of the variable [Sequence *]
 *
 * Return [SOFT ]  member variable baseseq [boolean]
 *
 */
boolean replace_baseseq_Genomic(Genomic * obj,Sequence * baseseq) 
{
    if( obj == NULL)     {  
      warn("In replacement function baseseq for object Genomic, got a NULL object"); 
      return FALSE;  
      }  
    obj->baseseq = baseseq;  
    return TRUE; 
}    


/* Function:  access_baseseq_Genomic(obj)
 *
 * Descrip:    Access member variable baseseq
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [Genomic *]
 *
 * Return [SOFT ]  member variable baseseq [Sequence *]
 *
 */
Sequence * access_baseseq_Genomic(Genomic * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function baseseq for object Genomic, got a NULL object");    
      return NULL;   
      }  
    return obj->baseseq;     
}    


/* Function:  access_repeat_Genomic(obj,i)
 *
 * Descrip:    Access members stored in the repeat list
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the list [Genomic *]
 * Arg:          i [UNKN ] Position in the list [int]
 *
 * Return [SOFT ]  Element of the list [GenomicRepeat *]
 *
 */
GenomicRepeat * access_repeat_Genomic(Genomic * obj,int i) 
{
    if( obj == NULL)     {  
      warn("In accessor function repeat for object Genomic, got a NULL object"); 
      return NULL;   
      }  
    if( obj->len <= i )  {  
      warn("In accessor function repeat for object Genomic, index %%d is greater than list length %%d",i,obj->len);  
      return NULL;   
      }  
    return obj->repeat[i];   
}    


/* Function:  length_repeat_Genomic(obj)
 *
 * Descrip:    discover the length of the list
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the list [Genomic *]
 *
 * Return [UNKN ]  length of the list [int]
 *
 */
int length_repeat_Genomic(Genomic * obj) 
{
    if( obj == NULL)     {  
      warn("In length function repeat for object Genomic, got a NULL object");   
      return -1;     
      }  
    return obj->len;     
}    


/* Function:  replace_start_GenomicRepeat(obj,start)
 *
 * Descrip:    Replace member variable start
 *             For use principly by API functions
 *
 *
 * Arg:          obj [UNKN ] Object holding the variable [GenomicRepeat *]
 * Arg:        start [OWNER] New value of the variable [int]
 *
 * Return [SOFT ]  member variable start [boolean]
 *
 */
boolean replace_start_GenomicRepeat(GenomicRepeat * obj,int start) 
{
    if( obj == NULL)     {  
      warn("In replacement function start for object GenomicRepeat, got a NULL object"); 
      return FALSE;  
      }  
    obj->start = start;  
    return TRUE; 
}    


/* Function:  access_start_GenomicRepeat(obj)
 *
 * Descrip:    Access member variable start
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [GenomicRepeat *]
 *
 * Return [SOFT ]  member variable start [int]
 *
 */
int access_start_GenomicRepeat(GenomicRepeat * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function start for object GenomicRepeat, got a NULL object");    
      return 0;  
      }  
    return obj->start;   
}    


/* Function:  replace_end_GenomicRepeat(obj,end)
 *
 * Descrip:    Replace member variable end
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [GenomicRepeat *]
 * Arg:        end [OWNER] New value of the variable [int]
 *
 * Return [SOFT ]  member variable end [boolean]
 *
 */
boolean replace_end_GenomicRepeat(GenomicRepeat * obj,int end) 
{
    if( obj == NULL)     {  
      warn("In replacement function end for object GenomicRepeat, got a NULL object");   
      return FALSE;  
      }  
    obj->end = end;  
    return TRUE; 
}    


/* Function:  access_end_GenomicRepeat(obj)
 *
 * Descrip:    Access member variable end
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [GenomicRepeat *]
 *
 * Return [SOFT ]  member variable end [int]
 *
 */
int access_end_GenomicRepeat(GenomicRepeat * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function end for object GenomicRepeat, got a NULL object");  
      return 0;  
      }  
    return obj->end;     
}    


/* Function:  replace_type_GenomicRepeat(obj,type)
 *
 * Descrip:    Replace member variable type
 *             For use principly by API functions
 *
 *
 * Arg:         obj [UNKN ] Object holding the variable [GenomicRepeat *]
 * Arg:        type [OWNER] New value of the variable [char *]
 *
 * Return [SOFT ]  member variable type [boolean]
 *
 */
boolean replace_type_GenomicRepeat(GenomicRepeat * obj,char * type) 
{
    if( obj == NULL)     {  
      warn("In replacement function type for object GenomicRepeat, got a NULL object");  
      return FALSE;  
      }  
    obj->type = type;    
    return TRUE; 
}    


/* Function:  access_type_GenomicRepeat(obj)
 *
 * Descrip:    Access member variable type
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [GenomicRepeat *]
 *
 * Return [SOFT ]  member variable type [char *]
 *
 */
char * access_type_GenomicRepeat(GenomicRepeat * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function type for object GenomicRepeat, got a NULL object"); 
      return NULL;   
      }  
    return obj->type;    
}    



#ifdef _cplusplus
}
#endif
