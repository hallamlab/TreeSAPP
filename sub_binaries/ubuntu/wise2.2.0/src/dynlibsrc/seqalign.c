#ifdef _cplusplus
extern "C" {
#endif
#include "seqalign.h"


/* Function:  ColumnCount_from_SeqAlign(sa,col)
 *
 * Descrip:    Gives you a column count
 *
 *             You are supposed to have enough sense to
 *             cache this on the caller if you need it again
 *
 *
 * Arg:         sa [UNKN ] Undocumented argument [SeqAlign *]
 * Arg:        col [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [ColumnCount *]
 *
 */
# line 48 "seqalign.dy"
ColumnCount * ColumnCount_from_SeqAlign(SeqAlign * sa,int col)
{
  ColumnCount * out;
  double count[27] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
		       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
		       0.0, 0.0, 0.0, 0.0, 0.0, 0.0 , 0.0};
  int i;
  
  if( sa == NULL || sa->seq[0]->len < col ) {
    warn("Cannot make column count at %d - out of range/seqalign null",col);
    return NULL;
  }

  for(i=0;i<sa->len;i++) {

    if( sa->seq[i]->len < col ) {
      warn("For seqalign sequence %d (%s), although first seq ok, this doesn't extend to %d. Error.",i,sa->seq[i]->name,col);
      return NULL;
    }
    if ( isalpha(sa->seq[i]->seq[col]) ) {
      count[toupper(sa->seq[i]->seq[col]) - 'A']++;
    } else {
      count[26]++;
    }

  }

  out = ColumnCount_alloc();
  for(i=0;i<27;i++) 
    out->count[i] = count[i];

  return out;
}

/* Function:  write_selex_SeqAlign(sa,name_len,block_len,ofp)
 *
 * Descrip:    Writes selex format
 *
 *
 * Arg:               sa [UNKN ] Undocumented argument [const SeqAlign *]
 * Arg:         name_len [UNKN ] Undocumented argument [int]
 * Arg:        block_len [UNKN ] Undocumented argument [int]
 * Arg:              ofp [UNKN ] Undocumented argument [FILE *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 85 "seqalign.dy"
boolean write_selex_SeqAlign(const SeqAlign * sa,int name_len,int block_len,FILE * ofp)
{
  int col;
  int i;
  int j;
  int k;


  for(col=0;col<sa->seq[0]->len;) {
    for(j=0;j<sa->len;j++) {
      for(k=0;k<name_len && sa->seq[j]->name[k] != '\0';k++) {
	fputc(sa->seq[j]->name[k],ofp);
      }
      if( sa->seq[j]->name[k] == '\0' && k >= name_len ) {
	warn("In printing selex alignment, name %s is truncated at %d",sa->seq[j]->name,k);
      }

      /* put in the rest of k */
      for(;k<name_len;k++)
	fputc(' ',ofp);

      /* add one space */
      fputc(' ',ofp);

      /* do the block */

      for(i=col;(i - col) < block_len && i < sa->seq[j]->len ; i++) {
	fputc(sa->seq[j]->seq[i],ofp);
      }
      fputc('\n',ofp);
    } /* end over all sequences for this block */

    col += block_len;
    fprintf(ofp,"\n\n");
  } 

  return TRUE;
}

/* Function:  read_selex_SeqAlign(ifp)
 *
 * Descrip:    Reads in selex (Stockholm) alignment
 *
 *             At the moment ignores all #= stuff on the
 *             sequence
 *
 *             Read HMMER documentation for a definition of
 *             the format
 *
 *
 * Arg:        ifp [UNKN ] Undocumented argument [FILE *]
 *
 * Return [UNKN ]  Undocumented return value [SeqAlign *]
 *
 */
# line 133 "seqalign.dy"
SeqAlign * read_selex_SeqAlign(FILE * ifp)
{
  char buffer[4096];
  SeqAlign * out;
  Sequence * temp;
  int i;
  boolean isfirst;
  boolean isfseq;
  boolean isend;
  char * name;
  char * seq;
  char * strip;
  int pos;

  out = SeqAlign_alloc_std();

  while( fgets(buffer,4096,ifp) != NULL ) {
    fprintf(stderr,"reading %s in selex alignment\n",buffer);
    if( buffer[0] == '#' ) 
      continue; /* skips over # */
    if( buffer[0] == '/' && buffer[1] == '/' )
      break;
    if( !isspace(buffer[0]) ) {
      /** in a block **/

      i = 0; /* this is the first sequence */
      if( out->len == 0 ) {
	isfirst = TRUE;
	isfseq  = TRUE;
      } else {
	isfirst = FALSE;
	isfseq  = FALSE;
      }

      for(;;) {
	/* strip trailing white space */

	for(strip = buffer + strlen(buffer) - 1;
	    strip > buffer && isspace(*strip);strip--)
	  ;

        strip++;
	*strip = '\0';
	
	for(name = seq = buffer;!isspace(*seq) && *seq; seq++)
	  ;

	/** if '\0' - worry **/
	
	if( !*seq ) {
	  warn("For name %s in selex alignment, no sequence!",name);
	  return out;
	}
	*seq = '\0';
	if( isfirst == TRUE ) {
	  temp = empty_Sequence_from_dynamic_memory(stringalloc(name));
	  add_SeqAlign(out,temp);
	} else {
	  if ( i >= out->len ) {
	    warn("For sequence %s, got over previous blocks size",name);
	    return NULL;
	  }

	  temp = out->seq[i++];
	}

	/* if is fseq, then find the next character */

	if( isfseq == TRUE ) {
	  for(seq++;*seq && isspace(*seq);seq++)
	    ;
	  if( *seq == '\0' ) {
	    warn("In first sequence %s, got to end of line with no sequence!",name);
	    return NULL;
	  }

	  isfseq = FALSE;
	  pos = seq - buffer;
	} 
	/* else pos is as normal */
	
	/* map all non alnum chars to '-' */
	
	for(seq = buffer+pos;*seq;seq++) {
	  if( !isalnum(*seq) ) {
	    *seq = '-';
	  }
	}

	/* add to sequence */
	add_string_to_Sequence(temp,buffer+pos);
	/* get next line */

	while( (seq = fgets(buffer,4096,ifp)) != NULL ) {
	  if( buffer[0] != '#' ) 
	    break;

	  if( buffer[0] == '/' && buffer[1] == '/' )
	    break;

	}
	if( seq == NULL ) {
	  isend = TRUE;
	  break;
	}

	/* if it is blank, break out of block */
	if( !isalnum(buffer[0]) ) {
	  break;
	}

	if( buffer[0] == '/' && buffer[1] == '/' )
	  break;

      }

      if( isend == TRUE )
	break;
    } /* end of was a non space line */

    /* else - ignore it !*/

  }

  return out;
}

      



# line 250 "seqalign.c"
/* Function:  swap_SeqAlign(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_SeqAlign
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [Sequence **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_SeqAlign(Sequence ** list,int i,int j)  
{
    Sequence * temp; 
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_SeqAlign(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_SeqAlign which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [Sequence **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_SeqAlign(Sequence ** list,int left,int right,int (*comp)(Sequence * ,Sequence * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_SeqAlign(list,left,(left+right)/2); 
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_SeqAlign (list,++last,i);   
      }  
    swap_SeqAlign (list,left,last);  
    qsort_SeqAlign(list,left,last-1,comp);   
    qsort_SeqAlign(list,last+1,right,comp);  
}    


/* Function:  sort_SeqAlign(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_SeqAlign
 *
 *
 * Arg:         obj [UNKN ] Object containing list [SeqAlign *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_SeqAlign(SeqAlign * obj,int (*comp)(Sequence *, Sequence *)) 
{
    qsort_SeqAlign(obj->seq,0,obj->len-1,comp);  
    return;  
}    


/* Function:  expand_SeqAlign(obj,len)
 *
 * Descrip:    Really an internal function for add_SeqAlign
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [SeqAlign *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_SeqAlign(SeqAlign * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_SeqAlign called with no need");   
      return TRUE;   
      }  


    if( (obj->seq = (Sequence ** ) ckrealloc (obj->seq,sizeof(Sequence *)*len)) == NULL)     {  
      warn("ckrealloc failed for expand_SeqAlign, returning FALSE"); 
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_SeqAlign(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [SeqAlign *]
 * Arg:        add [OWNER] Object to add to the list [Sequence *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_SeqAlign(SeqAlign * obj,Sequence * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_SeqAlign(obj,obj->len + SeqAlignLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->seq[obj->len++]=add;    
    return TRUE; 
}    


/* Function:  flush_SeqAlign(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [SeqAlign *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_SeqAlign(SeqAlign * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->seq[i] != NULL)   {  
        free_Sequence(obj->seq[i]);  
        obj->seq[i] = NULL;  
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  SeqAlign_alloc_std(void)
 *
 * Descrip:    Equivalent to SeqAlign_alloc_len(SeqAlignLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SeqAlign *]
 *
 */
SeqAlign * SeqAlign_alloc_std(void) 
{
    return SeqAlign_alloc_len(SeqAlignLISTLENGTH);   
}    


/* Function:  SeqAlign_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [SeqAlign *]
 *
 */
SeqAlign * SeqAlign_alloc_len(int len) 
{
    SeqAlign * out; /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = SeqAlign_alloc()) == NULL) 
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->seq = (Sequence ** ) ckcalloc (len,sizeof(Sequence *))) == NULL)    {  
      warn("Warning, ckcalloc failed in SeqAlign_alloc_len");    
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_SeqAlign(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [SeqAlign *]
 *
 * Return [UNKN ]  Undocumented return value [SeqAlign *]
 *
 */
SeqAlign * hard_link_SeqAlign(SeqAlign * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a SeqAlign object: passed a NULL object");    
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  SeqAlign_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SeqAlign *]
 *
 */
SeqAlign * SeqAlign_alloc(void) 
{
    SeqAlign * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(SeqAlign *) ckalloc (sizeof(SeqAlign))) == NULL)    {  
      warn("SeqAlign_alloc failed ");    
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
    out->name = NULL;    
    out->seq = NULL; 
    out->len = out->maxlen = 0;  


    return out;  
}    


/* Function:  free_SeqAlign(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [SeqAlign *]
 *
 * Return [UNKN ]  Undocumented return value [SeqAlign *]
 *
 */
SeqAlign * free_SeqAlign(SeqAlign * obj) 
{
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a SeqAlign obj. Should be trappable");  
      return NULL;   
      }  


    if( obj->dynamite_hard_link > 1)     {  
      obj->dynamite_hard_link--; 
      return NULL;   
      }  
    if( obj->name != NULL)   
      ckfree(obj->name);     
    if( obj->seq != NULL)    {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->seq[i] != NULL) 
          free_Sequence(obj->seq[i]);    
        }  
      ckfree(obj->seq);  
      }  


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_ColumnCount(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [ColumnCount *]
 *
 * Return [UNKN ]  Undocumented return value [ColumnCount *]
 *
 */
ColumnCount * hard_link_ColumnCount(ColumnCount * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a ColumnCount object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  ColumnCount_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ColumnCount *]
 *
 */
ColumnCount * ColumnCount_alloc(void) 
{
    ColumnCount * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(ColumnCount *) ckalloc (sizeof(ColumnCount))) == NULL)  {  
      warn("ColumnCount_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
    /* count[27] is an array: no default possible */ 


    return out;  
}    


/* Function:  free_ColumnCount(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [ColumnCount *]
 *
 * Return [UNKN ]  Undocumented return value [ColumnCount *]
 *
 */
ColumnCount * free_ColumnCount(ColumnCount * obj) 
{


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a ColumnCount obj. Should be trappable");   
      return NULL;   
      }  


    if( obj->dynamite_hard_link > 1)     {  
      obj->dynamite_hard_link--; 
      return NULL;   
      }  


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif
