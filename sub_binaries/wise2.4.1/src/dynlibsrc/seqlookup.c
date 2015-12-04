#ifdef _cplusplus
extern "C" {
#endif
#include "seqlookup.h"


/* Function:  free_SeqLookupInterface(sli)
 *
 * Descrip:    Frees SeqLookupInterface - overrides dynamite default
 *
 *
 * Arg:        sli [UNKN ] Undocumented argument [SeqLookupInterface *]
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupInterface *]
 *
 */
# line 112 "seqlookup.dy"
SeqLookupInterface * free_SeqLookupInterface(SeqLookupInterface * sli)
{
  int i;

  if( sli == NULL ) {
    return NULL;
  }

  for(i=0;i<sli->len;i++) {
    free_Sequence(sli->seq_store[i]);
  }

  (*sli->free_data)(sli->data);
  ckfree(sli);

  return NULL;

}


/* Function:  load_SequenceDB_SeqLookupLoadPara(p,db,sli)
 *
 * Descrip:    Loads a SequenceDB into a hash on the basis of the SeqLookupLoadPara
 *
 *             returns the number of sequences loaded
 *
 *
 * Arg:          p [UNKN ] Undocumented argument [SeqLookupLoadPara *]
 * Arg:         db [UNKN ] Undocumented argument [SequenceDB *]
 * Arg:        sli [UNKN ] Undocumented argument [SeqLookupInterface *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
# line 137 "seqlookup.dy"
int load_SequenceDB_SeqLookupLoadPara(SeqLookupLoadPara * p,SequenceDB * db,SeqLookupInterface * sli)
{
  int ret;
  int c;
  int i;
  Sequence * seq;
  int base[5];
  int start_base;
  int char_count = 0;


  for(i=0,start_base=1;i<5;i++) {
    base[i] = start_base;
    start_base = start_base * 26;
  }


  for(c=0,seq = init_SequenceDB(db,&ret); seq != NULL;seq = get_next_SequenceDB(db) ) {
    /* don't have to hardlink and then also free - just store here */
    c++;

    if( p->start_seq_load >= 0 && c-1 < p->start_seq_load ) {
      continue;
    }
    if( p->end_seq_load >= 0 && c-1 >= p->end_seq_load ) {
      break;
    }


    add_SeqLookupInterface(sli,seq);


    if( p->report_stagger >= 1 && c % p->report_stagger == 0 ) {
      info("Loaded %d sequences (%d characters, %.2f)... at %s\n",c,char_count,(double)char_count/(double)p->report_stagger,seq->name);
      char_count = 0;
    }

    char_count += seq->len;

    (*sli->add_seq)(sli->data,seq,p);

   
    if( p->truncate != 0 && c == p->truncate ) {
      info("Asked to truncate load after %d sequences\n",p->truncate);
      break;
    }
  }
    

  return c;

}




/* Function:  show_help_SeqLookupLoadPara(ofp)
 *
 * Descrip:    Shows help associated with Sequence loading
 *
 *
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 196 "seqlookup.dy"
void show_help_SeqLookupLoadPara(FILE * ofp)
{
  fprintf(ofp,"Sequence Index Loading parameters\n");
  fprintf(ofp,"   -seqloadtile             tiling skip of load (default 1)\n");
  fprintf(ofp,"   -seqloadreport <number>  report (using info) at what stagger rating (default none)\n");
  fprintf(ofp,"   -seqloadtrunc  <number>  truncate load after this number of sequeneces (useful for debugging)\n");
  fprintf(ofp,"   -seqloadstart  <number>  start position in database for seq load\n");
  fprintf(ofp,"   -seqloadend    <number>  end position in database for seq load\n");
  fprintf(ofp,"   -[no]seqloadlow          mark low complexity words for use with low complexity numbing\n");
 
  return;
}

/* Function:  new_SeqLookupLoadPara_from_argv(argc,argv)
 *
 * Descrip:    Builds new SeqLookup load from a command line
 *
 *
 * Arg:        argc [UNKN ] Undocumented argument [int *]
 * Arg:        argv [UNKN ] Undocumented argument [char **]
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupLoadPara *]
 *
 */
# line 212 "seqlookup.dy"
SeqLookupLoadPara * new_SeqLookupLoadPara_from_argv(int * argc,char ** argv)
{
  SeqLookupLoadPara * out;

  out = SeqLookupLoadPara_alloc();

  out->mark_low_complexity = 1;
  out->start_seq_load = -1;
  out->end_seq_load   = -1;

  strip_out_integer_argument(argc,argv,"seqloadtile",&out->tile_freq);

  strip_out_integer_argument(argc,argv,"seqloadreport",&out->report_stagger);

  strip_out_integer_argument(argc,argv,"seqloadtrunc",&out->truncate);

  strip_out_integer_argument(argc,argv,"seqloadstart",&out->start_seq_load);
  
  strip_out_integer_argument(argc,argv,"seqloadend",&out->end_seq_load);

  strip_out_boolean_def_argument(argc,argv,"seqloadlow",&out->mark_low_complexity);

  return out;
}



/* Function:  seq_number_dna_15mer_noN(seq)
 *
 * Descrip:    Function for DNA sequence to number on 15mers,
 *             Ns get mapped to -1
 *
 *
 * Arg:        seq [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
# line 243 "seqlookup.dy"
int seq_number_dna_15mer_noN(char * seq)
{
  int i;
  int ret = 0;
  int base = 1;
  int no = 0;

  for(i=0;i<15;i++) {
    no = base_from_char(seq[i]);
    if( no == BASE_N ) {
      return -1;
    }
    
    ret += base * no;
    base = base * 4;
  }

  return ret;
    

}

/* Function:  seq_number_dna_7mer_noN(seq)
 *
 * Descrip:    Function for DNA sequence to number on 15mers,
 *             Ns get mapped to -1
 *
 *
 * Arg:        seq [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
# line 269 "seqlookup.dy"
int seq_number_dna_7mer_noN(char * seq)
{
  int i;
  int ret = 0;
  int base = 1;
  int no = 0;

  for(i=0;i<7;i++) {
    no = base_from_char(seq[i]);
    if( no == BASE_N ) {
      return -1;
    }
    
    ret += base * no;
    base = base * 4;
  }

  return ret;
    
}


/* Function:  seq_number_aa_5mer(seq)
 *
 * Descrip:    Function for the amino acid to number on 5mers
 *
 *
 * Arg:        seq [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
# line 294 "seqlookup.dy"
int seq_number_aa_5mer(char * seq)
{
  int i;
  int ret = 0;
  int base = 1;
  int no = 0;

  for(i=0;i<5;i++) {
    no = toupper(seq[i])-'A';
    if( no > 26 || no < 0 ) {
      no = 'X'-'A';
    }
    ret += base * no;
    base = base * 26;
  }

  return ret;
}

/* Function:  flags_from_5aa_sequence(seq)
 *
 * Descrip:    returns simple lowcomplexity flag or 
 *             not for this sequence
 *
 *
 * Arg:        seq [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [char]
 *
 */
# line 317 "seqlookup.dy"
char flags_from_5aa_sequence(char * seq)
{
  int count =  0;
  
  if( seq[0] != seq[1] ) {
    count++;
  } 
  if( seq[1] != seq[2] ) {
    count++;
  } 
  if( seq[2] != seq[3] ) {
    count++;
  } 
  if( seq[3] != seq[4] ) {
    count++;
  } 
  if( seq[4] != seq[0] ) {
    count++;
  }

  if( count < 3 ) {
    return SEQLOOKUP_LOWCOMPLEXITY;
  } else {
    return 0;
  }

}

/* Function:  free_SeqLookupClientInterface(sli)
 *
 * Descrip:    Frees SeqLookupClientInterface - overrides dynamite default
 *
 *
 * Arg:        sli [UNKN ] Undocumented argument [SeqLookupClientInterface *]
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupClientInterface *]
 *
 */
# line 349 "seqlookup.dy"
SeqLookupClientInterface * free_SeqLookupClientInterface(SeqLookupClientInterface * sli)
{
  if( sli == NULL ) {
    return NULL;
  }

  (*sli->free_data)(sli->data);
  ckfree(sli);

  return NULL;

}

/* Function:  free_SeqLookupResultInterface(sli)
 *
 * Descrip:    Frees SeqLookupResultInterface - overrides dynamite default
 *
 *
 * Arg:        sli [UNKN ] Undocumented argument [SeqLookupResultInterface *]
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupResultInterface *]
 *
 */
# line 366 "seqlookup.dy"
SeqLookupResultInterface * free_SeqLookupResultInterface(SeqLookupResultInterface * sli)
{
  if( sli == NULL ) {
    return NULL;
  }

/*  fprintf(stderr,"Freeing results interface\n"); */
  (*sli->free_data)(sli->data);


  return NULL;

}





# line 340 "seqlookup.c"
/* Function:  hard_link_SeqLookupResultInterface(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [SeqLookupResultInterface *]
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupResultInterface *]
 *
 */
SeqLookupResultInterface * hard_link_SeqLookupResultInterface(SeqLookupResultInterface * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a SeqLookupResultInterface object: passed a NULL object");    
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  SeqLookupResultInterface_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupResultInterface *]
 *
 */
SeqLookupResultInterface * SeqLookupResultInterface_alloc(void) 
{
    SeqLookupResultInterface * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(SeqLookupResultInterface *) ckalloc (sizeof(SeqLookupResultInterface))) == NULL)    {  
      warn("SeqLookupResultInterface_alloc failed ");    
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->next = NULL;    
    out->is_more = NULL; 
    out->free_data = NULL;   


    return out;  
}    


/* Function:  hard_link_SeqLookupClientInterface(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [SeqLookupClientInterface *]
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupClientInterface *]
 *
 */
SeqLookupClientInterface * hard_link_SeqLookupClientInterface(SeqLookupClientInterface * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a SeqLookupClientInterface object: passed a NULL object");    
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  SeqLookupClientInterface_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupClientInterface *]
 *
 */
SeqLookupClientInterface * SeqLookupClientInterface_alloc(void) 
{
    SeqLookupClientInterface * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(SeqLookupClientInterface *) ckalloc (sizeof(SeqLookupClientInterface))) == NULL)    {  
      warn("SeqLookupClientInterface_alloc failed ");    
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->lookup = NULL;  
    out->is_populated = NULL;    
    out->free_data = NULL;   


    return out;  
}    


/* Function:  hard_link_SeqLookupLoadPara(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [SeqLookupLoadPara *]
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupLoadPara *]
 *
 */
SeqLookupLoadPara * hard_link_SeqLookupLoadPara(SeqLookupLoadPara * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a SeqLookupLoadPara object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  SeqLookupLoadPara_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupLoadPara *]
 *
 */
SeqLookupLoadPara * SeqLookupLoadPara_alloc(void) 
{
    SeqLookupLoadPara * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(SeqLookupLoadPara *) ckalloc (sizeof(SeqLookupLoadPara))) == NULL)  {  
      warn("SeqLookupLoadPara_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->type = SeqLookupLoad_Protein;   
    out->tile_freq = 1;  
    out->report_stagger = 0; 
    out->truncate = 0;   
    out->mark_low_complexity = 1;    
    out->start_seq_load = -1;    
    out->end_seq_load = -1;  


    return out;  
}    


/* Function:  free_SeqLookupLoadPara(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [SeqLookupLoadPara *]
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupLoadPara *]
 *
 */
SeqLookupLoadPara * free_SeqLookupLoadPara(SeqLookupLoadPara * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a SeqLookupLoadPara obj. Should be trappable"); 
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


/* Function:  swap_SeqLookupInterface(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_SeqLookupInterface
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [Sequence **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_SeqLookupInterface(Sequence ** list,int i,int j)  
{
    Sequence * temp; 
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_SeqLookupInterface(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_SeqLookupInterface which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [Sequence **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_SeqLookupInterface(Sequence ** list,int left,int right,int (*comp)(Sequence * ,Sequence * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_SeqLookupInterface(list,left,(left+right)/2);   
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_SeqLookupInterface (list,++last,i); 
      }  
    swap_SeqLookupInterface (list,left,last);    
    qsort_SeqLookupInterface(list,left,last-1,comp); 
    qsort_SeqLookupInterface(list,last+1,right,comp);    
}    


/* Function:  sort_SeqLookupInterface(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_SeqLookupInterface
 *
 *
 * Arg:         obj [UNKN ] Object containing list [SeqLookupInterface *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_SeqLookupInterface(SeqLookupInterface * obj,int (*comp)(Sequence *, Sequence *)) 
{
    qsort_SeqLookupInterface(obj->seq_store,0,obj->len-1,comp);  
    return;  
}    


/* Function:  expand_SeqLookupInterface(obj,len)
 *
 * Descrip:    Really an internal function for add_SeqLookupInterface
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [SeqLookupInterface *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_SeqLookupInterface(SeqLookupInterface * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_SeqLookupInterface called with no need"); 
      return TRUE;   
      }  


    if( (obj->seq_store = (Sequence ** ) ckrealloc (obj->seq_store,sizeof(Sequence *)*len)) == NULL)     {  
      warn("ckrealloc failed for expand_SeqLookupInterface, returning FALSE");   
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_SeqLookupInterface(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [SeqLookupInterface *]
 * Arg:        add [OWNER] Object to add to the list [Sequence *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_SeqLookupInterface(SeqLookupInterface * obj,Sequence * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_SeqLookupInterface(obj,obj->len + SeqLookupInterfaceLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->seq_store[obj->len++]=add;  
    return TRUE; 
}    


/* Function:  flush_SeqLookupInterface(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [SeqLookupInterface *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_SeqLookupInterface(SeqLookupInterface * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->seq_store[i] != NULL) {  
        free_Sequence(obj->seq_store[i]);    
        obj->seq_store[i] = NULL;    
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  SeqLookupInterface_alloc_std(void)
 *
 * Descrip:    Equivalent to SeqLookupInterface_alloc_len(SeqLookupInterfaceLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupInterface *]
 *
 */
SeqLookupInterface * SeqLookupInterface_alloc_std(void) 
{
    return SeqLookupInterface_alloc_len(SeqLookupInterfaceLISTLENGTH);   
}    


/* Function:  SeqLookupInterface_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupInterface *]
 *
 */
SeqLookupInterface * SeqLookupInterface_alloc_len(int len) 
{
    SeqLookupInterface * out;   /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = SeqLookupInterface_alloc()) == NULL)   
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->seq_store = (Sequence ** ) ckcalloc (len,sizeof(Sequence *))) == NULL)  {  
      warn("Warning, ckcalloc failed in SeqLookupInterface_alloc_len");  
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_SeqLookupInterface(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [SeqLookupInterface *]
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupInterface *]
 *
 */
SeqLookupInterface * hard_link_SeqLookupInterface(SeqLookupInterface * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a SeqLookupInterface object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  SeqLookupInterface_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupInterface *]
 *
 */
SeqLookupInterface * SeqLookupInterface_alloc(void) 
{
    SeqLookupInterface * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(SeqLookupInterface *) ckalloc (sizeof(SeqLookupInterface))) == NULL)    {  
      warn("SeqLookupInterface_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->get_client = NULL;  
    out->add_seq = NULL; 
    out->lookup_array_head = NULL;   
    out->add_direct_number = NULL;   
    out->free_data = NULL;   
    out->seq_store = NULL;   
    out->len = out->maxlen = 0;  


    return out;  
}    



#ifdef _cplusplus
}
#endif
