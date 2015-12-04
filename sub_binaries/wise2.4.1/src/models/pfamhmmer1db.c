#ifdef _cplusplus
extern "C" {
#endif
#include "pfamhmmer1db.h"


/* Function:  ThreeStateModel_from_name_PfamHmmer1DB(phd,name)
 *
 * Descrip:    reads a named model - akin to indexing
 *
 *
 * Arg:         phd [UNKN ] Undocumented argument [PfamHmmer1DB *]
 * Arg:        name [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateModel *]
 *
 */
# line 53 "pfamhmmer1db.dy"
ThreeStateModel * ThreeStateModel_from_name_PfamHmmer1DB(PfamHmmer1DB * phd,char * name)
{
  char buffer[512];
  ThreeStateModel * tsm;

  sprintf(buffer,"%s/%s.hmm",phd->dirname,name);

  /*  tsm = Wise2_read_ThreeStateModel_from_hmmer1_file(buffer); */
  tsm = HMMer2_read_ThreeStateModel(buffer);

  if( tsm == NULL ) {
    warn("Could not open Hmmer1 style hmm from Pfam db on file [%s]",buffer);
    return NULL;
  }

  if( tsm->name != NULL ) {
    ckfree(tsm->name);
  }
  tsm->name = stringalloc(name);

  /* ignore random stuff for the moment */

  return tsm;
}

/* Function:  read_next_TSM_PfamHmmer1DB(phd,return_status)
 *
 * Descrip:    reads the next threestatemodel
 *             for PfamHmmer1DB, placing the correct
 *             status into return_status( DB_RETURN_OK, etc).
 *
 *
 * Arg:                  phd [UNKN ] Undocumented argument [PfamHmmer1DB *]
 * Arg:        return_status [UNKN ] Undocumented argument [int *]
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateModel *]
 *
 */
# line 83 "pfamhmmer1db.dy"
ThreeStateModel * read_next_TSM_PfamHmmer1DB(PfamHmmer1DB * phd,int * return_status)
{
  ThreeStateModel * out;

  if( phd->cur >= phd->len ) {
    *return_status = DB_RETURN_END;
    return NULL;
  }
  out = read_TSM_from_PfamHmmer1Entry(phd->en[phd->cur++],phd->dirname);

  if( out == NULL ) {
    *return_status = DB_RETURN_ERROR;
    return NULL;
  } else {
    *return_status = DB_RETURN_OK;
    out->rm = hard_link_RandomModel(phd->def);
    return out;
  }
  return out;
}

/* Function:  read_TSM_from_PfamHmmer1Entry(en,dir)
 *
 * Descrip:    reads an individual HMMer model from the entry
 *             specification
 *
 *
 * Arg:         en [UNKN ] Undocumented argument [PfamHmmer1Entry *]
 * Arg:        dir [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateModel *]
 *
 */
# line 108 "pfamhmmer1db.dy"
ThreeStateModel * read_TSM_from_PfamHmmer1Entry(PfamHmmer1Entry * en,char * dir)
{
  char buffer[512];
  ThreeStateModel * tsm;

  sprintf(buffer,"%s/%s.hmm",dir,en->entryname);

  /*  tsm = Wise2_read_ThreeStateModel_from_hmmer1_file(buffer); */
  tsm = HMMer2_read_ThreeStateModel(buffer);

  if( tsm == NULL ) {
    warn("Could not open Hmmer1 style hmm from Pfam db on file [%s]",buffer);
    return NULL;
  }

  if( tsm->name != NULL ) {
    ckfree(tsm->name);
  }
  tsm->name = stringalloc(en->entryname);
  display_char_in_ThreeStateModel(tsm);

  /* ignore random stuff for the moment */

  if( en->is_hmmls == FALSE ) {
    force_weighted_local_model(tsm,1.0,1.0,1.0);
  } else {
    force_weighted_local_model(tsm,1.0,0.5,0.5);
  }

  return tsm;
}

/* Function:  PfamHmmer1DB_from_dirname(dirname)
 *
 * Descrip:    Makes a new PfamHmmer1DB from the dir name.
 *
 *             The directory should have a file HMMs which has
 *             entries like
 *
 *             rrm hmmls 12
 *
 *             with -r indicating a specialised random model.
 *
 *
 * Arg:        dirname [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [PfamHmmer1DB *]
 *
 */
# line 150 "pfamhmmer1db.dy"
PfamHmmer1DB * PfamHmmer1DB_from_dirname(char * dirname)
{
  char buffer[512];
  PfamHmmer1DB * out;
  FILE * ifp;
  char * runner;
  PfamHmmer1Entry * en;

  if( dirname == NULL ) {
    warn("passed through a NULL dirname into PfamHmmer1DB!");
    return NULL;
  }
  
  sprintf(buffer,"%s/HMMs",dirname);
  if( (ifp= openfile(buffer,"r")) == NULL ) {
    warn("Could not open %s as PfamHmmer1DB file list",buffer);
    return NULL;
  }

  out = PfamHmmer1DB_alloc_std();
  out->def = default_RandomModel();
  out->dirname = stringalloc(dirname);

  while( fgets(buffer,512,ifp) != NULL ) {
    if( (runner=strtok(buffer,spacestr)) == NULL) {
      continue; /* silently */
    }
    en  = PfamHmmer1Entry_alloc();
    en->entryname = stringalloc(runner);
    if( (runner=strtok(NULL,spacestr)) == NULL) {
      warn("Got a bad HMM.s line for a Pfam db. Skipping");
      free_PfamHmmer1Entry(en);
      continue;
    }
    if( strstr(runner,"hmmls") != NULL ) {
      en->is_hmmls = TRUE;
    } else {
      en->is_hmmls = FALSE;
    }


    if( strstr(runner,"-r") != NULL ) {
      en->is_random = TRUE;
    } else {
      en->is_random = FALSE;
    }

    if( (runner=strtok(NULL,spacestr)) == NULL) {
      warn("Got a bad HMM.s line for a Pfam db. Skipping");
      free_PfamHmmer1Entry(en);
      continue;
    }

    if( is_double_string(runner,&en->bits_cutoff) == FALSE ) {
      warn("%s does not look like a bits cutoff to me. Calling it 25",en->entryname,runner);
      en->bits_cutoff = 25;
    }
    add_PfamHmmer1DB(out,en);
  }

  fclose(ifp);
  return out;
}


# line 198 "pfamhmmer1db.c"
/* Function:  hard_link_PfamHmmer1Entry(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [PfamHmmer1Entry *]
 *
 * Return [UNKN ]  Undocumented return value [PfamHmmer1Entry *]
 *
 */
PfamHmmer1Entry * hard_link_PfamHmmer1Entry(PfamHmmer1Entry * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a PfamHmmer1Entry object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  PfamHmmer1Entry_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [PfamHmmer1Entry *]
 *
 */
PfamHmmer1Entry * PfamHmmer1Entry_alloc(void) 
{
    PfamHmmer1Entry * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(PfamHmmer1Entry *) ckalloc (sizeof(PfamHmmer1Entry))) == NULL)  {  
      warn("PfamHmmer1Entry_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->entryname = NULL;   
    out->is_random = FALSE;  
    out->is_hmmls = FALSE;   
    out->bits_cutoff = 0;    


    return out;  
}    


/* Function:  free_PfamHmmer1Entry(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [PfamHmmer1Entry *]
 *
 * Return [UNKN ]  Undocumented return value [PfamHmmer1Entry *]
 *
 */
PfamHmmer1Entry * free_PfamHmmer1Entry(PfamHmmer1Entry * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a PfamHmmer1Entry obj. Should be trappable");   
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
    if( obj->entryname != NULL)  
      ckfree(obj->entryname);    


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_PfamHmmer1DB(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_PfamHmmer1DB
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [PfamHmmer1Entry **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_PfamHmmer1DB(PfamHmmer1Entry ** list,int i,int j)  
{
    PfamHmmer1Entry * temp;  
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_PfamHmmer1DB(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_PfamHmmer1DB which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [PfamHmmer1Entry **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_PfamHmmer1DB(PfamHmmer1Entry ** list,int left,int right,int (*comp)(PfamHmmer1Entry * ,PfamHmmer1Entry * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_PfamHmmer1DB(list,left,(left+right)/2); 
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_PfamHmmer1DB (list,++last,i);   
      }  
    swap_PfamHmmer1DB (list,left,last);  
    qsort_PfamHmmer1DB(list,left,last-1,comp);   
    qsort_PfamHmmer1DB(list,last+1,right,comp);  
}    


/* Function:  sort_PfamHmmer1DB(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_PfamHmmer1DB
 *
 *
 * Arg:         obj [UNKN ] Object containing list [PfamHmmer1DB *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_PfamHmmer1DB(PfamHmmer1DB * obj,int (*comp)(PfamHmmer1Entry *, PfamHmmer1Entry *)) 
{
    qsort_PfamHmmer1DB(obj->en,0,obj->len-1,comp);   
    return;  
}    


/* Function:  expand_PfamHmmer1DB(obj,len)
 *
 * Descrip:    Really an internal function for add_PfamHmmer1DB
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [PfamHmmer1DB *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_PfamHmmer1DB(PfamHmmer1DB * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_PfamHmmer1DB called with no need");   
      return TRUE;   
      }  


    if( (obj->en = (PfamHmmer1Entry ** ) ckrealloc (obj->en,sizeof(PfamHmmer1Entry *)*len)) == NULL)     {  
      warn("ckrealloc failed for expand_PfamHmmer1DB, returning FALSE"); 
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_PfamHmmer1DB(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [PfamHmmer1DB *]
 * Arg:        add [OWNER] Object to add to the list [PfamHmmer1Entry *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_PfamHmmer1DB(PfamHmmer1DB * obj,PfamHmmer1Entry * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_PfamHmmer1DB(obj,obj->len + PfamHmmer1DBLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->en[obj->len++]=add; 
    return TRUE; 
}    


/* Function:  flush_PfamHmmer1DB(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [PfamHmmer1DB *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_PfamHmmer1DB(PfamHmmer1DB * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->en[i] != NULL)    {  
        free_PfamHmmer1Entry(obj->en[i]);    
        obj->en[i] = NULL;   
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  PfamHmmer1DB_alloc_std(void)
 *
 * Descrip:    Equivalent to PfamHmmer1DB_alloc_len(PfamHmmer1DBLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [PfamHmmer1DB *]
 *
 */
PfamHmmer1DB * PfamHmmer1DB_alloc_std(void) 
{
    return PfamHmmer1DB_alloc_len(PfamHmmer1DBLISTLENGTH);   
}    


/* Function:  PfamHmmer1DB_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [PfamHmmer1DB *]
 *
 */
PfamHmmer1DB * PfamHmmer1DB_alloc_len(int len) 
{
    PfamHmmer1DB * out; /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = PfamHmmer1DB_alloc()) == NULL) 
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->en = (PfamHmmer1Entry ** ) ckcalloc (len,sizeof(PfamHmmer1Entry *))) == NULL)   {  
      warn("Warning, ckcalloc failed in PfamHmmer1DB_alloc_len");    
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_PfamHmmer1DB(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [PfamHmmer1DB *]
 *
 * Return [UNKN ]  Undocumented return value [PfamHmmer1DB *]
 *
 */
PfamHmmer1DB * hard_link_PfamHmmer1DB(PfamHmmer1DB * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a PfamHmmer1DB object: passed a NULL object");    
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  PfamHmmer1DB_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [PfamHmmer1DB *]
 *
 */
PfamHmmer1DB * PfamHmmer1DB_alloc(void) 
{
    PfamHmmer1DB * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(PfamHmmer1DB *) ckalloc (sizeof(PfamHmmer1DB))) == NULL)    {  
      warn("PfamHmmer1DB_alloc failed ");    
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->en = NULL;  
    out->len = out->maxlen = 0;  
    out->dirname = NULL; 
    out->cur = 0;    
    out->def = NULL; 


    return out;  
}    


/* Function:  free_PfamHmmer1DB(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [PfamHmmer1DB *]
 *
 * Return [UNKN ]  Undocumented return value [PfamHmmer1DB *]
 *
 */
PfamHmmer1DB * free_PfamHmmer1DB(PfamHmmer1DB * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a PfamHmmer1DB obj. Should be trappable");  
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
    if( obj->en != NULL) {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->en[i] != NULL)  
          free_PfamHmmer1Entry(obj->en[i]);  
        }  
      ckfree(obj->en);   
      }  
    if( obj->dirname != NULL)    
      ckfree(obj->dirname);  
    if( obj->def != NULL)    
      free_RandomModel(obj->def);    


    ckfree(obj); 
    return NULL; 
}    


/* Function:  access_en_PfamHmmer1DB(obj,i)
 *
 * Descrip:    Access members stored in the en list
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the list [PfamHmmer1DB *]
 * Arg:          i [UNKN ] Position in the list [int]
 *
 * Return [SOFT ]  Element of the list [PfamHmmer1Entry *]
 *
 */
PfamHmmer1Entry * access_en_PfamHmmer1DB(PfamHmmer1DB * obj,int i) 
{
    if( obj == NULL)     {  
      warn("In accessor function en for object PfamHmmer1DB, got a NULL object");    
      return NULL;   
      }  
    if( obj->len <= i )  {  
      warn("In accessor function en for object PfamHmmer1DB, index %%d is greater than list length %%d",i,obj->len); 
      return NULL;   
      }  
    return obj->en[i];   
}    


/* Function:  length_en_PfamHmmer1DB(obj)
 *
 * Descrip:    discover the length of the list
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the list [PfamHmmer1DB *]
 *
 * Return [UNKN ]  length of the list [int]
 *
 */
int length_en_PfamHmmer1DB(PfamHmmer1DB * obj) 
{
    if( obj == NULL)     {  
      warn("In length function en for object PfamHmmer1DB, got a NULL object");  
      return -1;     
      }  
    return obj->len;     
}    


/* Function:  replace_dirname_PfamHmmer1DB(obj,dirname)
 *
 * Descrip:    Replace member variable dirname
 *             For use principly by API functions
 *
 *
 * Arg:            obj [UNKN ] Object holding the variable [PfamHmmer1DB *]
 * Arg:        dirname [OWNER] New value of the variable [char *]
 *
 * Return [SOFT ]  member variable dirname [boolean]
 *
 */
boolean replace_dirname_PfamHmmer1DB(PfamHmmer1DB * obj,char * dirname) 
{
    if( obj == NULL)     {  
      warn("In replacement function dirname for object PfamHmmer1DB, got a NULL object");    
      return FALSE;  
      }  
    obj->dirname = dirname;  
    return TRUE; 
}    


/* Function:  access_dirname_PfamHmmer1DB(obj)
 *
 * Descrip:    Access member variable dirname
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [PfamHmmer1DB *]
 *
 * Return [SOFT ]  member variable dirname [char *]
 *
 */
char * access_dirname_PfamHmmer1DB(PfamHmmer1DB * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function dirname for object PfamHmmer1DB, got a NULL object");   
      return NULL;   
      }  
    return obj->dirname;     
}    


/* Function:  replace_cur_PfamHmmer1DB(obj,cur)
 *
 * Descrip:    Replace member variable cur
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [PfamHmmer1DB *]
 * Arg:        cur [OWNER] New value of the variable [int]
 *
 * Return [SOFT ]  member variable cur [boolean]
 *
 */
boolean replace_cur_PfamHmmer1DB(PfamHmmer1DB * obj,int cur) 
{
    if( obj == NULL)     {  
      warn("In replacement function cur for object PfamHmmer1DB, got a NULL object");    
      return FALSE;  
      }  
    obj->cur = cur;  
    return TRUE; 
}    


/* Function:  access_cur_PfamHmmer1DB(obj)
 *
 * Descrip:    Access member variable cur
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [PfamHmmer1DB *]
 *
 * Return [SOFT ]  member variable cur [int]
 *
 */
int access_cur_PfamHmmer1DB(PfamHmmer1DB * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function cur for object PfamHmmer1DB, got a NULL object");   
      return 0;  
      }  
    return obj->cur;     
}    


/* Function:  replace_def_PfamHmmer1DB(obj,def)
 *
 * Descrip:    Replace member variable def
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [PfamHmmer1DB *]
 * Arg:        def [OWNER] New value of the variable [RandomModel *]
 *
 * Return [SOFT ]  member variable def [boolean]
 *
 */
boolean replace_def_PfamHmmer1DB(PfamHmmer1DB * obj,RandomModel * def) 
{
    if( obj == NULL)     {  
      warn("In replacement function def for object PfamHmmer1DB, got a NULL object");    
      return FALSE;  
      }  
    obj->def = def;  
    return TRUE; 
}    


/* Function:  access_def_PfamHmmer1DB(obj)
 *
 * Descrip:    Access member variable def
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [PfamHmmer1DB *]
 *
 * Return [SOFT ]  member variable def [RandomModel *]
 *
 */
RandomModel * access_def_PfamHmmer1DB(PfamHmmer1DB * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function def for object PfamHmmer1DB, got a NULL object");   
      return NULL;   
      }  
    return obj->def;     
}    


/* Function:  replace_entryname_PfamHmmer1Entry(obj,entryname)
 *
 * Descrip:    Replace member variable entryname
 *             For use principly by API functions
 *
 *
 * Arg:              obj [UNKN ] Object holding the variable [PfamHmmer1Entry *]
 * Arg:        entryname [OWNER] New value of the variable [char *]
 *
 * Return [SOFT ]  member variable entryname [boolean]
 *
 */
boolean replace_entryname_PfamHmmer1Entry(PfamHmmer1Entry * obj,char * entryname) 
{
    if( obj == NULL)     {  
      warn("In replacement function entryname for object PfamHmmer1Entry, got a NULL object");   
      return FALSE;  
      }  
    obj->entryname = entryname;  
    return TRUE; 
}    


/* Function:  access_entryname_PfamHmmer1Entry(obj)
 *
 * Descrip:    Access member variable entryname
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [PfamHmmer1Entry *]
 *
 * Return [SOFT ]  member variable entryname [char *]
 *
 */
char * access_entryname_PfamHmmer1Entry(PfamHmmer1Entry * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function entryname for object PfamHmmer1Entry, got a NULL object");  
      return NULL;   
      }  
    return obj->entryname;   
}    


/* Function:  replace_is_random_PfamHmmer1Entry(obj,is_random)
 *
 * Descrip:    Replace member variable is_random
 *             For use principly by API functions
 *
 *
 * Arg:              obj [UNKN ] Object holding the variable [PfamHmmer1Entry *]
 * Arg:        is_random [OWNER] New value of the variable [boolean]
 *
 * Return [SOFT ]  member variable is_random [boolean]
 *
 */
boolean replace_is_random_PfamHmmer1Entry(PfamHmmer1Entry * obj,boolean is_random) 
{
    if( obj == NULL)     {  
      warn("In replacement function is_random for object PfamHmmer1Entry, got a NULL object");   
      return FALSE;  
      }  
    obj->is_random = is_random;  
    return TRUE; 
}    


/* Function:  access_is_random_PfamHmmer1Entry(obj)
 *
 * Descrip:    Access member variable is_random
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [PfamHmmer1Entry *]
 *
 * Return [SOFT ]  member variable is_random [boolean]
 *
 */
boolean access_is_random_PfamHmmer1Entry(PfamHmmer1Entry * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function is_random for object PfamHmmer1Entry, got a NULL object");  
      return FALSE;  
      }  
    return obj->is_random;   
}    


/* Function:  replace_is_hmmls_PfamHmmer1Entry(obj,is_hmmls)
 *
 * Descrip:    Replace member variable is_hmmls
 *             For use principly by API functions
 *
 *
 * Arg:             obj [UNKN ] Object holding the variable [PfamHmmer1Entry *]
 * Arg:        is_hmmls [OWNER] New value of the variable [boolean]
 *
 * Return [SOFT ]  member variable is_hmmls [boolean]
 *
 */
boolean replace_is_hmmls_PfamHmmer1Entry(PfamHmmer1Entry * obj,boolean is_hmmls) 
{
    if( obj == NULL)     {  
      warn("In replacement function is_hmmls for object PfamHmmer1Entry, got a NULL object");    
      return FALSE;  
      }  
    obj->is_hmmls = is_hmmls;    
    return TRUE; 
}    


/* Function:  access_is_hmmls_PfamHmmer1Entry(obj)
 *
 * Descrip:    Access member variable is_hmmls
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [PfamHmmer1Entry *]
 *
 * Return [SOFT ]  member variable is_hmmls [boolean]
 *
 */
boolean access_is_hmmls_PfamHmmer1Entry(PfamHmmer1Entry * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function is_hmmls for object PfamHmmer1Entry, got a NULL object");   
      return FALSE;  
      }  
    return obj->is_hmmls;    
}    


/* Function:  replace_bits_cutoff_PfamHmmer1Entry(obj,bits_cutoff)
 *
 * Descrip:    Replace member variable bits_cutoff
 *             For use principly by API functions
 *
 *
 * Arg:                obj [UNKN ] Object holding the variable [PfamHmmer1Entry *]
 * Arg:        bits_cutoff [OWNER] New value of the variable [double]
 *
 * Return [SOFT ]  member variable bits_cutoff [boolean]
 *
 */
boolean replace_bits_cutoff_PfamHmmer1Entry(PfamHmmer1Entry * obj,double bits_cutoff) 
{
    if( obj == NULL)     {  
      warn("In replacement function bits_cutoff for object PfamHmmer1Entry, got a NULL object"); 
      return FALSE;  
      }  
    obj->bits_cutoff = bits_cutoff;  
    return TRUE; 
}    


/* Function:  access_bits_cutoff_PfamHmmer1Entry(obj)
 *
 * Descrip:    Access member variable bits_cutoff
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [PfamHmmer1Entry *]
 *
 * Return [SOFT ]  member variable bits_cutoff [double]
 *
 */
double access_bits_cutoff_PfamHmmer1Entry(PfamHmmer1Entry * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function bits_cutoff for object PfamHmmer1Entry, got a NULL object");    
      return 0;  
      }  
    return obj->bits_cutoff;     
}    



#ifdef _cplusplus
}
#endif
