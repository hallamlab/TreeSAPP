#ifdef _cplusplus
extern "C" {
#endif
#include "packaln.h"

/* Function:  read_simple_PackAln_file(file)
 *
 * Descrip:    Reads in a PackAln from a file in show_simple_PackAln
 *
 *
 * Arg:        file [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
# line 57 "packaln.dy"
PackAln * read_simple_PackAln_file(char * file)
{
  FILE * ifp;
  PackAln * out;

  if( (ifp=openfile(file,"r")) == NULL ) {
    warn("Could not open %s for PackAln file reading",file);
  }

  out = read_simple_PackAln(ifp);

  fclose(ifp);

  return out;
  
}

/* Function:  read_simple_PackAln(ifp)
 *
 * Descrip:    Reads in PackAln from file format in show_simple_PackAln
 *
 *
 * Arg:        ifp [UNKN ] Undocumented argument [FILE *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
# line 77 "packaln.dy"
PackAln * read_simple_PackAln(FILE * ifp)
{
  PackAln * out;
  PackAlnUnit * unit;
  char buffer[MAXLINE];

  assert(ifp);

  out = PackAln_alloc_std();

  /* score line */
  fgets(buffer,MAXLINE,ifp);

  while( fgets(buffer,MAXLINE,ifp) != NULL ) {
    if( strstartcmp(buffer,"//") == 0 ) {
      break;
    }

    unit = PackAlnUnit_alloc();
    
    sscanf(buffer,"Position i:[%d] j:[%d] State:[%d] Score: %d",&unit->i,&unit->j,&unit->state,&unit->score);

    add_PackAln(out,unit);
  }

  return out;
}


/* Function:  show_simple_PackAlnUnit(pau,ofp)
 *
 * Descrip:    shows packalnunit very simply ;)
 *
 *
 * Arg:        pau [UNKN ] Undocumented argument [PackAlnUnit *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 110 "packaln.dy"
void show_simple_PackAlnUnit(PackAlnUnit * pau,FILE * ofp)
{
  fprintf(ofp,"Position i:[%d] j:[%d] State:[%d] Score: %d\n",pau->i,pau->j,pau->state,pau->score);
}

/* Function:  show_text_PackAlnUnit(state_to_char,pau,ofp)
 *
 * Descrip:    shows packalnunit with the text mapping
 *
 *
 * Arg:        state_to_char [UNKN ] Undocumented argument [NullString]
 * Arg:                  pau [UNKN ] Undocumented argument [PackAlnUnit *]
 * Arg:                  ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 119 "packaln.dy"
void show_text_PackAlnUnit(PackAlnUnit * pau,char * (*state_to_char)(int),FILE * ofp)
{
  fprintf(ofp,"Position i:[%4d] j:[%4d] State:[%2d] Score:[%3d] [%s]\n",pau->i,pau->j,pau->state,pau->score,(*state_to_char)(pau->state));
}

/* Function:  show_bits_and_cumlative_PackAln(pal,ofp)
 *
 * Descrip:    Shows packaln as: 
 *
 *             i,j,state,score,bits,cumlative-score,cumlative-bits
 *
 *             cumlative score and cumlative bits are useful sometimes
 *
 *
 * Arg:        pal [UNKN ] Undocumented argument [PackAln *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 131 "packaln.dy"
void show_bits_and_cumlative_PackAln(PackAln * pal,FILE * ofp)
{
  int i;
  int cs = 0;
  double cb = 0.0;


  fprintf(ofp,"Score %d\n",pal->score);
  for(i=0;i<pal->len;i++) {
    auto PackAlnUnit * pau;
    pau = pal->pau[i];

    cs += pal->pau[i]->score;
    cb += Score2Bits(pal->pau[i]->score);
    fprintf(ofp,"i [%4d] j [%4d] state [%2d] score [%4d] bits [%2.2f] Score-CF [%6d] Bits-CF[%4.2f]\n",pau->i,pau->j,pau->state,pau->score,Score2Bits(pau->score),cs,cb);
  }

} 

 
/* Function:  show_simple_PackAln(pal,ofp)
 *
 * Descrip:    shows packaln with a pretty verbose debugging 
 *             format
 *
 *
 * Arg:        pal [UNKN ] Undocumented argument [PackAln *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 155 "packaln.dy"
void show_simple_PackAln(PackAln * pal,FILE * ofp)
{
  register int i;

  fprintf(ofp,"Score %d\n",pal->score);
  for(i=0;i<pal->len;i++)
    show_simple_PackAlnUnit(pal->pau[i],ofp);
}

/* Function:  show_text_PackAln(state_to_char,pal,ofp)
 *
 * Descrip:    shows packaln with a pretty verbose debugging 
 *             format, but with a conversion function from state number to
 *             a string
 *
 *
 * Arg:        state_to_char [UNKN ] Undocumented argument [NullString]
 * Arg:                  pal [UNKN ] Undocumented argument [PackAln *]
 * Arg:                  ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 169 "packaln.dy"
void show_text_PackAln(PackAln * pal,char * (*state_to_char)(int),FILE * ofp)
{
  register int i;

  fprintf(ofp,"Score %d\n",pal->score);
  for(i=0;i<pal->len;i++)
    show_text_PackAlnUnit(pal->pau[i],state_to_char,ofp);
}

/* Function:  invert_PackAln(pal)
 *
 * Descrip:    inverts the packaln so that the last unit is the first
 *             etc. Because most alignments are read backwards this
 *             is useful
 *
 *
 * Arg:        pal [UNKN ] PackAln to be inverted  [PackAln *]
 *
 */
# line 185 "packaln.dy"
void invert_PackAln(PackAln * pal) 
{
  PackAlnUnit ** temp;
  register int i;

  /*** there are better ways to do this! ***/

  temp = (PackAlnUnit **) ckcalloc(pal->len,sizeof(PackAlnUnit *));
  
  for(i=0;i<pal->len;i++) 
    temp[i] = pal->pau[pal->len-1-i];
  for(i=0;i<pal->len;i++)
    pal->pau[i] = temp[i];

  ckfree(temp);
}



# line 194 "packaln.c"
/* Function:  hard_link_PackAlnUnit(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [PackAlnUnit *]
 *
 * Return [UNKN ]  Undocumented return value [PackAlnUnit *]
 *
 */
PackAlnUnit * hard_link_PackAlnUnit(PackAlnUnit * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a PackAlnUnit object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  PackAlnUnit_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [PackAlnUnit *]
 *
 */
PackAlnUnit * PackAlnUnit_alloc(void) 
{
    PackAlnUnit * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(PackAlnUnit *) ckalloc (sizeof(PackAlnUnit))) == NULL)  {  
      warn("PackAlnUnit_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->i = 0;  
    out->j = 0;  
    out->state = 0;  
    out->score = 0;  


    return out;  
}    


/* Function:  free_PackAlnUnit(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [PackAlnUnit *]
 *
 * Return [UNKN ]  Undocumented return value [PackAlnUnit *]
 *
 */
PackAlnUnit * free_PackAlnUnit(PackAlnUnit * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a PackAlnUnit obj. Should be trappable");   
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


/* Function:  swap_PackAln(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_PackAln
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [PackAlnUnit **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_PackAln(PackAlnUnit ** list,int i,int j)  
{
    PackAlnUnit * temp;  
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_PackAln(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_PackAln which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [PackAlnUnit **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_PackAln(PackAlnUnit ** list,int left,int right,int (*comp)(PackAlnUnit * ,PackAlnUnit * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_PackAln(list,left,(left+right)/2);  
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_PackAln (list,++last,i);    
      }  
    swap_PackAln (list,left,last);   
    qsort_PackAln(list,left,last-1,comp);    
    qsort_PackAln(list,last+1,right,comp);   
}    


/* Function:  sort_PackAln(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_PackAln
 *
 *
 * Arg:         obj [UNKN ] Object containing list [PackAln *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_PackAln(PackAln * obj,int (*comp)(PackAlnUnit *, PackAlnUnit *)) 
{
    qsort_PackAln(obj->pau,0,obj->len-1,comp);   
    return;  
}    


/* Function:  expand_PackAln(obj,len)
 *
 * Descrip:    Really an internal function for add_PackAln
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [PackAln *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_PackAln(PackAln * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_PackAln called with no need");    
      return TRUE;   
      }  


    if( (obj->pau = (PackAlnUnit ** ) ckrealloc (obj->pau,sizeof(PackAlnUnit *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_PackAln, returning FALSE");  
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_PackAln(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [PackAln *]
 * Arg:        add [OWNER] Object to add to the list [PackAlnUnit *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_PackAln(PackAln * obj,PackAlnUnit * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_PackAln(obj,obj->len + PackAlnLISTLENGTH) == FALSE) 
        return FALSE;    
      }  


    obj->pau[obj->len++]=add;    
    return TRUE; 
}    


/* Function:  flush_PackAln(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [PackAln *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_PackAln(PackAln * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->pau[i] != NULL)   {  
        free_PackAlnUnit(obj->pau[i]);   
        obj->pau[i] = NULL;  
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  PackAln_alloc_std(void)
 *
 * Descrip:    Equivalent to PackAln_alloc_len(PackAlnLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * PackAln_alloc_std(void) 
{
    return PackAln_alloc_len(PackAlnLISTLENGTH); 
}    


/* Function:  PackAln_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * PackAln_alloc_len(int len) 
{
    PackAln * out;  /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = PackAln_alloc()) == NULL)  
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->pau = (PackAlnUnit ** ) ckcalloc (len,sizeof(PackAlnUnit *))) == NULL)  {  
      warn("Warning, ckcalloc failed in PackAln_alloc_len"); 
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_PackAln(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [PackAln *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * hard_link_PackAln(PackAln * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a PackAln object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  PackAln_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * PackAln_alloc(void) 
{
    PackAln * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(PackAln *) ckalloc (sizeof(PackAln))) == NULL)  {  
      warn("PackAln_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->pau = NULL; 
    out->len = out->maxlen = 0;  
    out->score = 0;  


    return out;  
}    


/* Function:  free_PackAln(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [PackAln *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * free_PackAln(PackAln * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a PackAln obj. Should be trappable");   
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
    if( obj->pau != NULL)    {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->pau[i] != NULL) 
          free_PackAlnUnit(obj->pau[i]); 
        }  
      ckfree(obj->pau);  
      }  


    ckfree(obj); 
    return NULL; 
}    


/* Function:  access_pau_PackAln(obj,i)
 *
 * Descrip:    Access members stored in the pau list
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the list [PackAln *]
 * Arg:          i [UNKN ] Position in the list [int]
 *
 * Return [SOFT ]  Element of the list [PackAlnUnit *]
 *
 */
PackAlnUnit * access_pau_PackAln(PackAln * obj,int i) 
{
    if( obj == NULL)     {  
      warn("In accessor function pau for object PackAln, got a NULL object");    
      return NULL;   
      }  
    if( obj->len <= i )  {  
      warn("In accessor function pau for object PackAln, index %%d is greater than list length %%d",i,obj->len); 
      return NULL;   
      }  
    return obj->pau[i];  
}    


/* Function:  length_pau_PackAln(obj)
 *
 * Descrip:    discover the length of the list
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the list [PackAln *]
 *
 * Return [UNKN ]  length of the list [int]
 *
 */
int length_pau_PackAln(PackAln * obj) 
{
    if( obj == NULL)     {  
      warn("In length function pau for object PackAln, got a NULL object");  
      return -1;     
      }  
    return obj->len;     
}    


/* Function:  replace_score_PackAln(obj,score)
 *
 * Descrip:    Replace member variable score
 *             For use principly by API functions
 *
 *
 * Arg:          obj [UNKN ] Object holding the variable [PackAln *]
 * Arg:        score [OWNER] New value of the variable [int]
 *
 * Return [SOFT ]  member variable score [boolean]
 *
 */
boolean replace_score_PackAln(PackAln * obj,int score) 
{
    if( obj == NULL)     {  
      warn("In replacement function score for object PackAln, got a NULL object");   
      return FALSE;  
      }  
    obj->score = score;  
    return TRUE; 
}    


/* Function:  access_score_PackAln(obj)
 *
 * Descrip:    Access member variable score
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [PackAln *]
 *
 * Return [SOFT ]  member variable score [int]
 *
 */
int access_score_PackAln(PackAln * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function score for object PackAln, got a NULL object");  
      return 0;  
      }  
    return obj->score;   
}    


/* Function:  replace_i_PackAlnUnit(obj,i)
 *
 * Descrip:    Replace member variable i
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [PackAlnUnit *]
 * Arg:          i [OWNER] New value of the variable [int]
 *
 * Return [SOFT ]  member variable i [boolean]
 *
 */
boolean replace_i_PackAlnUnit(PackAlnUnit * obj,int i) 
{
    if( obj == NULL)     {  
      warn("In replacement function i for object PackAlnUnit, got a NULL object");   
      return FALSE;  
      }  
    obj->i = i;  
    return TRUE; 
}    


/* Function:  access_i_PackAlnUnit(obj)
 *
 * Descrip:    Access member variable i
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [PackAlnUnit *]
 *
 * Return [SOFT ]  member variable i [int]
 *
 */
int access_i_PackAlnUnit(PackAlnUnit * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function i for object PackAlnUnit, got a NULL object");  
      return 0;  
      }  
    return obj->i;   
}    


/* Function:  replace_j_PackAlnUnit(obj,j)
 *
 * Descrip:    Replace member variable j
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [PackAlnUnit *]
 * Arg:          j [OWNER] New value of the variable [int]
 *
 * Return [SOFT ]  member variable j [boolean]
 *
 */
boolean replace_j_PackAlnUnit(PackAlnUnit * obj,int j) 
{
    if( obj == NULL)     {  
      warn("In replacement function j for object PackAlnUnit, got a NULL object");   
      return FALSE;  
      }  
    obj->j = j;  
    return TRUE; 
}    


/* Function:  access_j_PackAlnUnit(obj)
 *
 * Descrip:    Access member variable j
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [PackAlnUnit *]
 *
 * Return [SOFT ]  member variable j [int]
 *
 */
int access_j_PackAlnUnit(PackAlnUnit * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function j for object PackAlnUnit, got a NULL object");  
      return 0;  
      }  
    return obj->j;   
}    


/* Function:  replace_state_PackAlnUnit(obj,state)
 *
 * Descrip:    Replace member variable state
 *             For use principly by API functions
 *
 *
 * Arg:          obj [UNKN ] Object holding the variable [PackAlnUnit *]
 * Arg:        state [OWNER] New value of the variable [int]
 *
 * Return [SOFT ]  member variable state [boolean]
 *
 */
boolean replace_state_PackAlnUnit(PackAlnUnit * obj,int state) 
{
    if( obj == NULL)     {  
      warn("In replacement function state for object PackAlnUnit, got a NULL object");   
      return FALSE;  
      }  
    obj->state = state;  
    return TRUE; 
}    


/* Function:  access_state_PackAlnUnit(obj)
 *
 * Descrip:    Access member variable state
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [PackAlnUnit *]
 *
 * Return [SOFT ]  member variable state [int]
 *
 */
int access_state_PackAlnUnit(PackAlnUnit * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function state for object PackAlnUnit, got a NULL object");  
      return 0;  
      }  
    return obj->state;   
}    


/* Function:  replace_score_PackAlnUnit(obj,score)
 *
 * Descrip:    Replace member variable score
 *             For use principly by API functions
 *
 *
 * Arg:          obj [UNKN ] Object holding the variable [PackAlnUnit *]
 * Arg:        score [OWNER] New value of the variable [int]
 *
 * Return [SOFT ]  member variable score [boolean]
 *
 */
boolean replace_score_PackAlnUnit(PackAlnUnit * obj,int score) 
{
    if( obj == NULL)     {  
      warn("In replacement function score for object PackAlnUnit, got a NULL object");   
      return FALSE;  
      }  
    obj->score = score;  
    return TRUE; 
}    


/* Function:  access_score_PackAlnUnit(obj)
 *
 * Descrip:    Access member variable score
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [PackAlnUnit *]
 *
 * Return [SOFT ]  member variable score [int]
 *
 */
int access_score_PackAlnUnit(PackAlnUnit * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function score for object PackAlnUnit, got a NULL object");  
      return 0;  
      }  
    return obj->score;   
}    



#ifdef _cplusplus
}
#endif
