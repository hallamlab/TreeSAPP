#ifdef _cplusplus
extern "C" {
#endif
#include "complexconsensi.h"

/* Function:  show_ComplexConsensi(cc,*ofp)
 *
 * Descrip:    shows complexconsensi in vaguely human form
 *
 *
 * Arg:          cc [UNKN ] Undocumented argument [ComplexConsensi *]
 * Arg:        *ofp [UNKN ] Undocumented argument [FILE]
 *
 */
# line 31 "complexconsensi.dy"
void show_ComplexConsensi(ComplexConsensi * cc,FILE *ofp)
{
  register int i;

  for(i=0;i<cc->len;i++)
    show_ComplexConsensusWord(cc->ccw[i],ofp);

}

/* Function:  show_ComplexConsensusWord(ccw,ofp)
 *
 * Descrip:    shows an individual ccword
 *
 *
 * Arg:        ccw [UNKN ] Undocumented argument [ComplexConsensusWord *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 45 "complexconsensi.dy"
void show_ComplexConsensusWord(ComplexConsensusWord * ccw,FILE * ofp)
{
  fprintf(ofp,"%s %4.4f %d\n",ccw->pattern,ccw->p,ccw->score);
}


/* Function:  word_from_ComplexConsensi(word,cc)
 *
 * Descrip:    Way of getting probabilities out of a consensus.
 *             If it is in the Consensus, then gets prob, otherwise 0.
 *
 *
 * Arg:        word [UNKN ] Undocumented argument [char *]
 * Arg:          cc [UNKN ] Undocumented argument [ComplexConsensi *]
 *
 * Return [UNKN ]  Undocumented return value [Probability]
 *
 */
# line 56 "complexconsensi.dy"
Probability word_from_ComplexConsensi(char * word,ComplexConsensi * cc)
{
  ComplexConsensusWord * ccw;

  ccw = best_ComplexConsensusWord(word,cc);

  if( ccw == NULL )
    return 0.0;

  return ccw->p;
}

/* Function:  score_from_ComplexConsensi(word,cc)
 *
 * Descrip:    Way of getting scores out of a consensus.
 *             If it is in the Consensus, then gets score, otherwise NEGI.
 *
 *
 * Arg:        word [UNKN ] Undocumented argument [char *]
 * Arg:          cc [UNKN ] Undocumented argument [ComplexConsensi *]
 *
 * Return [UNKN ]  Undocumented return value [Score]
 *
 */
# line 73 "complexconsensi.dy"
Score score_from_ComplexConsensi(char * word,ComplexConsensi * cc)
{
  ComplexConsensusWord * ccw;

  ccw = best_ComplexConsensusWord(word,cc);

  if( ccw == NULL )
    return NEGI;

  return ccw->score;
}

/* Function:  best_ComplexConsensusWord(word,cc)
 *
 * Descrip:    Finds the best (highest) match to this word
 *
 *
 * Arg:        word [UNKN ] Undocumented argument [char *]
 * Arg:          cc [UNKN ] Undocumented argument [ComplexConsensi *]
 *
 * Return [UNKN ]  Undocumented return value [ComplexConsensusWord *]
 *
 */
# line 89 "complexconsensi.dy"
ComplexConsensusWord * best_ComplexConsensusWord(char * word,ComplexConsensi * cc)
{
  register int i;

  for(i=0;i<cc->len;i++) {
    if( word_matches_ComplexConsensusWord(word,cc->ccw[i]) == TRUE)
      return cc->ccw[i];
  }

  return NULL;
}
    

/* Function:  word_matches_ComplexConsensusWord(word,ccw)
 *
 * Descrip:    Core of the matching system. Checks that word matches 
 *             ComplexConsensusWord. This says that '-' matches anything.
 *             Issues a warning if hits a '\0' in word
 *
 *
 * Arg:        word [UNKN ] Undocumented argument [char *]
 * Arg:         ccw [UNKN ] Undocumented argument [ComplexConsensusWord *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 108 "complexconsensi.dy"
boolean word_matches_ComplexConsensusWord(char * word,ComplexConsensusWord * ccw)
{
  return strcmp_with_dashes(word,ccw->pattern);
} 

/* Function:  strcmp_with_dashes(word,pattern)
 *
 * Descrip:    The matching 'function' used by
 *             /word_matches_ComplexConsensusWord
 *
 *
 * Arg:           word [UNKN ] Undocumented argument [char *]
 * Arg:        pattern [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 118 "complexconsensi.dy"
boolean strcmp_with_dashes(char * word,char * pattern)
{
  char * runner;
  char * base = word;

  for(runner = pattern;*runner && *word;word++,runner++) {
    if( *runner == '-') {
      continue;
    }

    if( *runner != *word ) {
      break;
    }

  }

  if( *word == '\0' && *runner != '\0') {
    warn("Tried to match a word [%s] which is shorter than the pattern [%s]",base,pattern);
    return FALSE;
  }

  if( *runner == '\0') {
    /* end of pattern */
    return TRUE;
  }

  return FALSE;
}




 /************************/
 /* I/O                  */
 /************************/


/* Function:  read_ComplexConsensi_file(filename)
 *
 * Descrip:    Reads a file containing the ComplexConsensi.
 *             Not every useful, as most times these consensi
 *             are in one file, with other things
 *
 *
 * Arg:        filename [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [ComplexConsensi *]
 *
 */
# line 161 "complexconsensi.dy"
ComplexConsensi * read_ComplexConsensi_file(char * filename)
{
  FILE * ifp;
  ComplexConsensi * out;

  ifp = openfile(filename,"r");
  if( ifp == NULL ) {
    warn("Could not open file %s for ComplexConsensi",filename);
    return NULL;
  }

  out = read_ComplexConsensi(ifp);

  fclose(ifp);

  return out;
}


/* Function:  read_ComplexConsensi(ifp)
 *
 * Descrip:    Reads on ComplexConsensi from the FILE ifp.
 *
 *
 * Arg:        ifp [UNKN ] input filestream [FILE *]
 *
 * Return [UNKN ]  Undocumented return value [ComplexConsensi *]
 *
 */
# line 185 "complexconsensi.dy"
ComplexConsensi * read_ComplexConsensi(FILE * ifp)
{
  ComplexConsensi * out;
  ComplexConsensusWord * temp;
  char buffer[MAXLINE];


  out = ComplexConsensi_alloc_std();
  if( out == NULL )
    return NULL;


  while( fgets(buffer,MAXLINE,ifp) != NULL ) {
    if( buffer[0] == '#' ) 
      continue;
    temp = read_ComplexConsensusWord_line(buffer);
    if( temp != NULL ) 
      add_ComplexConsensi(out,temp);
  }

  return out;
}
    

/* Function:  read_ComplexConsensusWord_line(line)
 *
 * Descrip:    Reads a single ccword from a line
 *
 *
 * Arg:        line [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [ComplexConsensusWord *]
 *
 */
# line 214 "complexconsensi.dy"
ComplexConsensusWord * read_ComplexConsensusWord_line(char * line)
{
  ComplexConsensusWord * out;
  char * pattern;
  char * score;


  pattern = strtok(line,spacestr);
  score   = strtok(NULL,spacestr);

  if( pattern == NULL || score == NULL ) {
    warn("Could not read a ComplexConsenusWord... ooops");
    return NULL;
  }


  out = ComplexConsensusWord_alloc();

  if( out == NULL)
    return NULL;

  out->pattern = stringalloc(pattern);
  out->score = (double) atof(score);

  return out;
}


# line 270 "complexconsensi.c"
/* Function:  hard_link_ComplexConsensusWord(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [ComplexConsensusWord *]
 *
 * Return [UNKN ]  Undocumented return value [ComplexConsensusWord *]
 *
 */
ComplexConsensusWord * hard_link_ComplexConsensusWord(ComplexConsensusWord * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a ComplexConsensusWord object: passed a NULL object");    
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  ComplexConsensusWord_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ComplexConsensusWord *]
 *
 */
ComplexConsensusWord * ComplexConsensusWord_alloc(void) 
{
    ComplexConsensusWord * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(ComplexConsensusWord *) ckalloc (sizeof(ComplexConsensusWord))) == NULL)    {  
      warn("ComplexConsensusWord_alloc failed ");    
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->pattern = NULL; 
    out->score = 0;  
    out->p = 0.0;    


    return out;  
}    


/* Function:  free_ComplexConsensusWord(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [ComplexConsensusWord *]
 *
 * Return [UNKN ]  Undocumented return value [ComplexConsensusWord *]
 *
 */
ComplexConsensusWord * free_ComplexConsensusWord(ComplexConsensusWord * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a ComplexConsensusWord obj. Should be trappable");  
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
    if( obj->pattern != NULL)    
      ckfree(obj->pattern);  


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_ComplexConsensi(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_ComplexConsensi
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [ComplexConsensusWord **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_ComplexConsensi(ComplexConsensusWord ** list,int i,int j)  
{
    ComplexConsensusWord * temp; 
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_ComplexConsensi(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_ComplexConsensi which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [ComplexConsensusWord **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_ComplexConsensi(ComplexConsensusWord ** list,int left,int right,int (*comp)(ComplexConsensusWord * ,ComplexConsensusWord * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_ComplexConsensi(list,left,(left+right)/2);  
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_ComplexConsensi (list,++last,i);    
      }  
    swap_ComplexConsensi (list,left,last);   
    qsort_ComplexConsensi(list,left,last-1,comp);    
    qsort_ComplexConsensi(list,last+1,right,comp);   
}    


/* Function:  sort_ComplexConsensi(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_ComplexConsensi
 *
 *
 * Arg:         obj [UNKN ] Object containing list [ComplexConsensi *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_ComplexConsensi(ComplexConsensi * obj,int (*comp)(ComplexConsensusWord *, ComplexConsensusWord *)) 
{
    qsort_ComplexConsensi(obj->ccw,0,obj->len-1,comp);   
    return;  
}    


/* Function:  expand_ComplexConsensi(obj,len)
 *
 * Descrip:    Really an internal function for add_ComplexConsensi
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [ComplexConsensi *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_ComplexConsensi(ComplexConsensi * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_ComplexConsensi called with no need");    
      return TRUE;   
      }  


    if( (obj->ccw = (ComplexConsensusWord ** ) ckrealloc (obj->ccw,sizeof(ComplexConsensusWord *)*len)) == NULL)     {  
      warn("ckrealloc failed for expand_ComplexConsensi, returning FALSE");  
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_ComplexConsensi(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [ComplexConsensi *]
 * Arg:        add [OWNER] Object to add to the list [ComplexConsensusWord *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_ComplexConsensi(ComplexConsensi * obj,ComplexConsensusWord * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_ComplexConsensi(obj,obj->len + ComplexConsensiLISTLENGTH) == FALSE) 
        return FALSE;    
      }  


    obj->ccw[obj->len++]=add;    
    return TRUE; 
}    


/* Function:  flush_ComplexConsensi(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [ComplexConsensi *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_ComplexConsensi(ComplexConsensi * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->ccw[i] != NULL)   {  
        free_ComplexConsensusWord(obj->ccw[i]);  
        obj->ccw[i] = NULL;  
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  ComplexConsensi_alloc_std(void)
 *
 * Descrip:    Equivalent to ComplexConsensi_alloc_len(ComplexConsensiLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ComplexConsensi *]
 *
 */
ComplexConsensi * ComplexConsensi_alloc_std(void) 
{
    return ComplexConsensi_alloc_len(ComplexConsensiLISTLENGTH); 
}    


/* Function:  ComplexConsensi_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [ComplexConsensi *]
 *
 */
ComplexConsensi * ComplexConsensi_alloc_len(int len) 
{
    ComplexConsensi * out;  /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = ComplexConsensi_alloc()) == NULL)  
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->ccw = (ComplexConsensusWord ** ) ckcalloc (len,sizeof(ComplexConsensusWord *))) == NULL)    {  
      warn("Warning, ckcalloc failed in ComplexConsensi_alloc_len"); 
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_ComplexConsensi(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [ComplexConsensi *]
 *
 * Return [UNKN ]  Undocumented return value [ComplexConsensi *]
 *
 */
ComplexConsensi * hard_link_ComplexConsensi(ComplexConsensi * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a ComplexConsensi object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  ComplexConsensi_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ComplexConsensi *]
 *
 */
ComplexConsensi * ComplexConsensi_alloc(void) 
{
    ComplexConsensi * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(ComplexConsensi *) ckalloc (sizeof(ComplexConsensi))) == NULL)  {  
      warn("ComplexConsensi_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->ccw = NULL; 
    out->len = out->maxlen = 0;  
    out->name = NULL;    


    return out;  
}    


/* Function:  free_ComplexConsensi(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [ComplexConsensi *]
 *
 * Return [UNKN ]  Undocumented return value [ComplexConsensi *]
 *
 */
ComplexConsensi * free_ComplexConsensi(ComplexConsensi * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a ComplexConsensi obj. Should be trappable");   
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
    if( obj->ccw != NULL)    {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->ccw[i] != NULL) 
          free_ComplexConsensusWord(obj->ccw[i]);    
        }  
      ckfree(obj->ccw);  
      }  
    if( obj->name != NULL)   
      ckfree(obj->name);     


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif
