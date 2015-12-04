#ifdef _cplusplus
extern "C" {
#endif
#include "ftext.h"

/* Function:  add_line_to_Ftext(ft,str,)
 *
 * Descrip:    adds a vsprintf'd line (below MAXLINE length!) to
 *             the Ftext.
 *
 *
 * Arg:         ft [UNKN ] Undocumented argument [Ftext *]
 * Arg:        str [UNKN ] Undocumented argument [char *]
 * Arg:            [UNKN ] Undocumented argument [.]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 26 "ftext.dy"
boolean add_line_to_Ftext(Ftext * ft,char * str,...)
{
  char buffer[MAXLINE];
  va_list ap;


  va_start(ap,str);
  vsprintf(buffer,str,ap);
  va_end(ap);

  return add_Fblock(ft->fb[ft->len-1],stringalloc(buffer));
}

/* Function:  add_break_to_Ftext(ft)
 *
 * Descrip:    puts in a break into the Ftext
 *
 *
 * Arg:        ft [UNKN ] Undocumented argument [Ftext *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 42 "ftext.dy"
boolean add_break_to_Ftext(Ftext * ft)
{
  Fblock * temp;

  temp = Fblock_alloc_std();
  return add_Ftext(ft,temp);
}

/* Function:  single_Ftext_from_str(str)
 *
 * Descrip:    Makes a complete Ftext from just this string:
 *
 *
 * Arg:        str [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [Ftext *]
 *
 */
# line 53 "ftext.dy"
Ftext * single_Ftext_from_str(char * str)
{
  Fblock * temp;
  Ftext * out;

  out = Ftext_alloc_len(1);

  temp = Fblock_alloc_len(1);
  
  add_Fblock(temp,stringalloc(str));
  add_Ftext(out,temp);

  return out;
}
  
/* Function:  show_eddystyle_Ftext(ft,header,depth,ofp,blank_text)
 *
 * Descrip:    shows Ftext as a C style comment with
 *               *
 *               * 
 *               *
 *             indenting: 
 *
 *             will return number of lines printed.
 *
 *
 * Arg:                ft [READ ] Ftext to be shown [Ftext *]
 * Arg:            header [READ ] Header for the first line, eg, "description:" [char *]
 * Arg:             depth [READ ] depth of from * to text [int]
 * Arg:               ofp [UNKN ] output file [FILE *]
 * Arg:        blank_text [READ ] if non NULL, what to put if ft is empty, Can be NULL. [char *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
# line 83 "ftext.dy"
int show_eddystyle_Ftext(Ftext * ft,char * header,int depth,FILE * ofp,char * blank_text)
{
  int len = 0;
  int i;
  int j;
  int st = 1;
  char space[]="                                                                    ";

  if( ft == NULL || ft->len ==  0 || ft->fb[0]->len == 0) {
    fprintf(ofp," * %*s %s\n",depth-9,header,blank_text == NULL ? "No text" : blank_text);
    return 1;
  }


  fprintf(ofp," * %*s %s\n",-(depth-4),header,ft->fb[0]->line[0]);
  len++;

  for(i=0;i<ft->len;i++) {
    for(j=0;j<ft->fb[i]->len;j++) {
      if( st == 1 ) {
	st =0;
	continue;
      }
      fprintf(ofp," *%.*s%s\n",depth-2,space,ft->fb[i]->line[j]);
      len++;
    }
    fprintf(ofp," *\n");
    len++;
  }

  return len;
}
      

/* Function:  latex_Ftext(ft,ofp)
 *
 * Descrip:    Provides a latex dump of some text.
 *
 *             Lines that start flush to the left
 *             are given as paragraphs
 *
 *             Line that are indented are made verbatim
 *
 *
 * Arg:         ft [UNKN ] Undocumented argument [Ftext *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 125 "ftext.dy"
void latex_Ftext(Ftext * ft,FILE * ofp)
{
  int i,j;
  boolean inverbatim = FALSE;

  for(i=0;i<ft->len;i++) {
    for(j=0;j<ft->fb[i]->len;j++) {
      if( isspace(ft->fb[i]->line[j][0]) ) {
	if( inverbatim == FALSE ) {
	  fprintf(ofp,"\\begin{verbatim}\n");
	  inverbatim = TRUE;
	} 
      } else if ( inverbatim == TRUE ) {
	fprintf(ofp,"\\end{verbatim}\n");
	inverbatim = FALSE;
      } 
      
      fprintf(ofp,"%s\n",ft->fb[i]->line[j]);
    }
    fprintf(ofp,"\n\n");
  }

  /* if we are in verbatim - get out! */
  if ( inverbatim == TRUE ) 
    fprintf(ofp,"\\end{verbatim}\n");
}


      
/* Function:  dump_Ftext(ft,*ofp)
 *
 * Descrip:    stupid function which gives flat dump of ftext
 *
 *
 * Arg:          ft [UNKN ] Undocumented argument [Ftext *]
 * Arg:        *ofp [UNKN ] Undocumented argument [FILE]
 *
 */
# line 157 "ftext.dy"
void dump_Ftext(Ftext * ft,FILE *ofp)
{
  int i;

  for(i=0;i<ft->len;i++) {
    dump_Fblock_str("",ft->fb[i],ofp);
    fprintf(ofp,"\n");
  }

}

/* Function:  dump_Ftext_pre(pre,ft,*ofp)
 *
 * Descrip:    stupid function which gives flat dump of ftext,
 *             with a start of pre 
 *
 *
 * Arg:         pre [UNKN ] Undocumented argument [char *]
 * Arg:          ft [UNKN ] Undocumented argument [Ftext *]
 * Arg:        *ofp [UNKN ] Undocumented argument [FILE]
 *
 */
# line 172 "ftext.dy"
void dump_Ftext_pre(char * pre,Ftext * ft,FILE *ofp)
{
  int i;

  for(i=0;i<ft->len;i++) {
    dump_Fblock_str(pre,ft->fb[i],ofp);
    fprintf(ofp,"\n");
  }

}

/* Function:  dump_Fblock_str(pre,fb,ofp)
 *
 * Descrip:    sub of /dump_Ftext
 *
 *
 * Arg:        pre [UNKN ] Undocumented argument [char *]
 * Arg:         fb [UNKN ] Undocumented argument [Fblock *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 187 "ftext.dy"
void dump_Fblock_str(char * pre,Fblock * fb,FILE * ofp)
{
  int i;

  for(i=0;i<fb->len;i++) {
    fprintf(ofp,"%s%s\n",pre,fb->line[i]);
  }
}

/* Function:  read_Ftext(buffer,fgets_func,maxlen,ifp,endpoint)
 *
 * Descrip:    reads in lines until it hits endpoint.
 *             A bit of an internal: is going to use buffer and maxlen
 *             as its buffer (this is so it is not fixed to one length of
 *             buffer). You probably have abuffer in your calling function.
 *
 *
 * Arg:            buffer [WRITE] pointer to a char * buffer of maxlen that can be written in [char *]
 * Arg:        fgets_func [UNKN ] Undocumented argument [NullString]
 * Arg:            maxlen [READ ] maximum size of buffer to read [int]
 * Arg:               ifp [UNKN ] input file [FILE *]
 * Arg:          endpoint [READ ] a string which is used (using /strstartcmp) as a tag for the end of text [char *]
 *
 * Return [UNKN ]  Undocumented return value [Ftext *]
 *
 */
# line 207 "ftext.dy"
Ftext * read_Ftext(char * buffer,int maxlen,FILE * ifp,char * endpoint,char * (*fgets_func)(char *,int,FILE *))
{
  Ftext * out;
  Fblock * temp;

  out = Ftext_alloc_std();
  if( fgets_func == NULL  )
    fgets_func = fgets;

  while((temp= read_Fblock(buffer,maxlen,ifp,endpoint,fgets_func)) != NULL ) {
    add_Ftext(out,temp);
    
    if( strstartcmp(buffer,endpoint) == 0 )
      return out;
  }

  warn("Got a NULL Fblock in reading a Ftext. Not good news!");
  return out;
}


/* Function:  read_Fblock(buffer,fgets_func,maxlen,ifp,endpoint)
 *
 * Descrip:    Really an internal for read_Ftext
 *
 *
 * Arg:            buffer [UNKN ] Undocumented argument [char *]
 * Arg:        fgets_func [UNKN ] Undocumented argument [NullString]
 * Arg:            maxlen [UNKN ] Undocumented argument [int]
 * Arg:               ifp [UNKN ] Undocumented argument [FILE *]
 * Arg:          endpoint [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [Fblock *]
 *
 */
# line 232 "ftext.dy"
Fblock * read_Fblock(char * buffer,int maxlen,FILE * ifp,char * endpoint,char * (*fgets_func)(char *,int,FILE *))
{
  Fblock * out;

  out = Fblock_alloc_std();


  while((*fgets_func)(buffer,maxlen,ifp) != NULL ) {
    if( strstartcmp(buffer,endpoint) == 0 )
      return out;
    if( strstartcmp(buffer,"\n") == 0 ) 
      return out;
    buffer[strlen(buffer)-1] ='\0'; /*** strips off '\n' ***/

    /*    fprintf(stderr,"Adding %s",buffer); */

    add_Fblock(out,stringalloc(buffer));

  }

  warn("Got to then end of the file in reading a Fblock. Not a good sign!");
  return out;
}    




# line 301 "ftext.c"
/* Function:  swap_Fblock(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_Fblock
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [char **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_Fblock(char ** list,int i,int j)  
{
    char * temp; 
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_Fblock(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_Fblock which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [char **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_Fblock(char ** list,int left,int right,int (*comp)(char * ,char * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_Fblock(list,left,(left+right)/2);   
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_Fblock (list,++last,i); 
      }  
    swap_Fblock (list,left,last);    
    qsort_Fblock(list,left,last-1,comp); 
    qsort_Fblock(list,last+1,right,comp);    
}    


/* Function:  sort_Fblock(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_Fblock
 *
 *
 * Arg:         obj [UNKN ] Object containing list [Fblock *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_Fblock(Fblock * obj,int (*comp)(char *, char *)) 
{
    qsort_Fblock(obj->line,0,obj->len-1,comp);   
    return;  
}    


/* Function:  expand_Fblock(obj,len)
 *
 * Descrip:    Really an internal function for add_Fblock
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [Fblock *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_Fblock(Fblock * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_Fblock called with no need"); 
      return TRUE;   
      }  


    if( (obj->line = (char ** ) ckrealloc (obj->line,sizeof(char *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_Fblock, returning FALSE");   
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_Fblock(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [Fblock *]
 * Arg:        add [OWNER] Object to add to the list [char *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_Fblock(Fblock * obj,char * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_Fblock(obj,obj->len + FblockLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->line[obj->len++]=add;   
    return TRUE; 
}    


/* Function:  flush_Fblock(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [Fblock *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_Fblock(Fblock * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->line[i] != NULL)  {  
        ckfree(obj->line[i]);    
        obj->line[i] = NULL; 
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  Fblock_alloc_std(void)
 *
 * Descrip:    Equivalent to Fblock_alloc_len(FblockLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Fblock *]
 *
 */
Fblock * Fblock_alloc_std(void) 
{
    return Fblock_alloc_len(FblockLISTLENGTH);   
}    


/* Function:  Fblock_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [Fblock *]
 *
 */
Fblock * Fblock_alloc_len(int len) 
{
    Fblock * out;   /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = Fblock_alloc()) == NULL)   
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->line = (char ** ) ckcalloc (len,sizeof(char *))) == NULL)   {  
      warn("Warning, ckcalloc failed in Fblock_alloc_len");  
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_Fblock(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [Fblock *]
 *
 * Return [UNKN ]  Undocumented return value [Fblock *]
 *
 */
Fblock * hard_link_Fblock(Fblock * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a Fblock object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  Fblock_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Fblock *]
 *
 */
Fblock * Fblock_alloc(void) 
{
    Fblock * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(Fblock *) ckalloc (sizeof(Fblock))) == NULL)    {  
      warn("Fblock_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->line = NULL;    
    out->len = out->maxlen = 0;  


    return out;  
}    


/* Function:  free_Fblock(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [Fblock *]
 *
 * Return [UNKN ]  Undocumented return value [Fblock *]
 *
 */
Fblock * free_Fblock(Fblock * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a Fblock obj. Should be trappable");    
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
    if( obj->line != NULL)   {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->line[i] != NULL)    
          ckfree(obj->line[i]);  
        }  
      ckfree(obj->line); 
      }  


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_Ftext(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_Ftext
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [Fblock **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_Ftext(Fblock ** list,int i,int j)  
{
    Fblock * temp;   
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_Ftext(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_Ftext which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [Fblock **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_Ftext(Fblock ** list,int left,int right,int (*comp)(Fblock * ,Fblock * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_Ftext(list,left,(left+right)/2);    
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_Ftext (list,++last,i);  
      }  
    swap_Ftext (list,left,last); 
    qsort_Ftext(list,left,last-1,comp);  
    qsort_Ftext(list,last+1,right,comp); 
}    


/* Function:  sort_Ftext(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_Ftext
 *
 *
 * Arg:         obj [UNKN ] Object containing list [Ftext *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_Ftext(Ftext * obj,int (*comp)(Fblock *, Fblock *)) 
{
    qsort_Ftext(obj->fb,0,obj->len-1,comp);  
    return;  
}    


/* Function:  expand_Ftext(obj,len)
 *
 * Descrip:    Really an internal function for add_Ftext
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [Ftext *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_Ftext(Ftext * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_Ftext called with no need");  
      return TRUE;   
      }  


    if( (obj->fb = (Fblock ** ) ckrealloc (obj->fb,sizeof(Fblock *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_Ftext, returning FALSE");    
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_Ftext(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [Ftext *]
 * Arg:        add [OWNER] Object to add to the list [Fblock *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_Ftext(Ftext * obj,Fblock * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_Ftext(obj,obj->len + FtextLISTLENGTH) == FALSE) 
        return FALSE;    
      }  


    obj->fb[obj->len++]=add; 
    return TRUE; 
}    


/* Function:  flush_Ftext(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [Ftext *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_Ftext(Ftext * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->fb[i] != NULL)    {  
        free_Fblock(obj->fb[i]); 
        obj->fb[i] = NULL;   
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  Ftext_alloc_std(void)
 *
 * Descrip:    Equivalent to Ftext_alloc_len(FtextLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Ftext *]
 *
 */
Ftext * Ftext_alloc_std(void) 
{
    return Ftext_alloc_len(FtextLISTLENGTH); 
}    


/* Function:  Ftext_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [Ftext *]
 *
 */
Ftext * Ftext_alloc_len(int len) 
{
    Ftext * out;/* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = Ftext_alloc()) == NULL)    
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->fb = (Fblock ** ) ckcalloc (len,sizeof(Fblock *))) == NULL) {  
      warn("Warning, ckcalloc failed in Ftext_alloc_len");   
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_Ftext(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [Ftext *]
 *
 * Return [UNKN ]  Undocumented return value [Ftext *]
 *
 */
Ftext * hard_link_Ftext(Ftext * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a Ftext object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  Ftext_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Ftext *]
 *
 */
Ftext * Ftext_alloc(void) 
{
    Ftext * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(Ftext *) ckalloc (sizeof(Ftext))) == NULL)  {  
      warn("Ftext_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->fb = NULL;  
    out->len = out->maxlen = 0;  


    return out;  
}    


/* Function:  free_Ftext(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [Ftext *]
 *
 * Return [UNKN ]  Undocumented return value [Ftext *]
 *
 */
Ftext * free_Ftext(Ftext * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a Ftext obj. Should be trappable"); 
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
    if( obj->fb != NULL) {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->fb[i] != NULL)  
          free_Fblock(obj->fb[i]);   
        }  
      ckfree(obj->fb);   
      }  


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif
