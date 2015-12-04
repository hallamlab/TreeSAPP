#ifdef _cplusplus
extern "C" {
#endif
#include "sequencedb.h"



/* Function:  get_Sequence_from_SequenceDB(sdb,de)
 *
 * Descrip:    Quite a mindless function which retrieves sequences
 *             via indexes
 *
 *             Going to spend too much time in fopen if this is used
 *             too much
 *
 *
 * Arg:        sdb [UNKN ] Undocumented argument [SequenceDB *]
 * Arg:         de [UNKN ] Undocumented argument [DataEntry *]
 *
 * Return [UNKN ]  Undocumented return value [Sequence *]
 *
 */
# line 76 "sequencedb.dy"
Sequence * get_Sequence_from_SequenceDB(SequenceDB * sdb,DataEntry * de)
{
  FILE * ifp;
  Sequence * ret;

  if( de == NULL ) {
    warn("Cannot get sequence database entry with a null dataentry!");
    return NULL;
  }

  if( sdb == NULL ) {
    warn("Cannot get sequence database entry with a null sequence db!");
    return NULL;
  }

  if( de->filename == NULL ) {
    warn("Cannot get sequence database entry with no attached filename");
    return NULL;
  }


  /* actually, all our info is in dataentry */

  ifp = openfile(de->filename,"r");
  if( ifp == NULL ) {
    warn("Bad error - could not open database file %s for reading indexed sequence",de->filename);
    return NULL;
  }

  fseek(ifp,de->byte_position,SEEK_SET);

  switch(de->data[1]) {
  case SEQ_DB_FASTA : 
    ret = read_fasta_Sequence(ifp);
    break;
  default :
    warn("Unknown SequenceDB type [%d]",de->data[1]);
    ret = NULL;
  }

  fclose(ifp);

  return ret;
}



/* Function:  add_SequenceDB_info_DataEntry(sdb,de)
 *
 * Descrip:    A function which places data into dataentry so we can
 *             be guarenteed to retrieve it sometime.
 *
 *             It uses 0 and 1 points in the Data array.
 *
 *
 * Arg:        sdb [UNKN ] Undocumented argument [SequenceDB *]
 * Arg:         de [UNKN ] Undocumented argument [DataEntry *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 129 "sequencedb.dy"
boolean add_SequenceDB_info_DataEntry(SequenceDB * sdb,DataEntry * de)
{
  if( sdb == NULL || de == NULL ) {
    warn("Null objects being passed into add_SequenceDB_info_DataEntry. Can't be good!");
    return FALSE;
  }

  de->filename = sdb->fs[sdb->current_source]->filename; /* if there... */
  de->byte_position  = sdb->byte_position; /* of this sequence */
  de->data[1]  = sdb->fs[sdb->current_source]->format;

  return TRUE;
}
								      

/* Function:  close_SequenceDB(last,sdb)
 *
 * Descrip:    top level function that closes the SequenceDB
 *             after the last sequence is read.
 *
 *
 * Arg:        last [WRITE] Sequence object to be freed  [Sequence *]
 * Arg:         sdb [READ ] database to be closed [SequenceDB *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 151 "sequencedb.dy"
boolean close_SequenceDB(Sequence * last,SequenceDB * sdb)
{
  if( last != NULL )
    free_Sequence(last);

  if( sdb->sequence_no == 1 && sdb->has_warned_single == 0) {
    info("Your sequence database has only sequence in it. It is quite likely there was a more efficient way to run this");
    sdb->has_warned_single = 1;
  }
  /*** nothing else to do? ***/

  sdb->current_source = (-1);
  return TRUE;
}


/* Function:  init_SequenceDB(sdb,return_status)
 *
 * Descrip:    top level function that starts a database read on
 *             SequenceDB
 *
 *
 *
 * Arg:                  sdb [READ ] sequence database [SequenceDB *]
 * Arg:        return_status [WRITE] returns the database status as found in database.h [int *]
 *
 * Return [UNKN ]  Undocumented return value [Sequence *]
 *
 */
# line 175 "sequencedb.dy"
Sequence * init_SequenceDB(SequenceDB * sdb,int * return_status)
{
  Sequence * temp = NULL;
  int count;

  sdb->current_source = 0;
  sdb->sequence_no =0;
  load_next_fs_SequenceDB(sdb);

  if( sdb->seq_start != -1 && sdb->seq_end != -1 ) {
    for(count=0;count <= sdb->seq_start;count++) {
      temp = reload_SequenceDB(temp,sdb,return_status);
    }
    return temp;
  }

  return reload_SequenceDB(NULL,sdb,return_status);
}

/* Function:  reload_SequenceDB(last,sdb,return_status)
 *
 * Descrip:    top level function that reloads a sequence database
 *
 *
 *
 * Arg:                 last [WRITE] previous sequence to be used: will simply be freed at the moment [Sequence *]
 * Arg:                  sdb [UNKN ] sequence database [SequenceDB *]
 * Arg:        return_status [WRITE] returns the database status as found in database.h [int *]
 *
 * Return [UNKN ]  Undocumented return value [Sequence *]
 *
 */
# line 202 "sequencedb.dy"
Sequence * reload_SequenceDB(Sequence * last,SequenceDB * sdb,int * return_status)
{
  Sequence * out;
  int count = 0;

  /*
   * free last Sequence: if we did something clever with
   * memory, this is where we should do it 
   */

  if( last != NULL ) 
    free_Sequence(last);


  /* if there is a seq_end, then see whether this is the end */

  if( sdb->seq_end != -1 && sdb->sequence_no == sdb->seq_end ) {
    /* end */
    sdb->current_source = -1;
    *return_status = DB_RETURN_END;
    return NULL;
  }

  /** see if we can read a Sequence now **/


  if( (out = get_next_SequenceDB(sdb)) != NULL ) {
    *return_status = DB_RETURN_OK;
    sdb->sequence_no++;
    return out;
  }

 
  if( SequenceDB_at_end(sdb) == TRUE ) {
    if( close_last_fs_SequenceDB(sdb) == FALSE ) {
      warn("On file source [%d] [%s] could not close",sdb->current_source,sdb->fs[sdb->current_source]->filename);
      *return_status = DB_RETURN_ERROR;
      return NULL;
    }
    *return_status = DB_RETURN_END;
    return NULL;
  }

  /** ok, see if we can swap FileSources then **/

  for(;;) {
    if( close_last_fs_SequenceDB(sdb) == FALSE ) {
      warn("On file source [%d] [%s] could not close",sdb->current_source,sdb->fs[sdb->current_source]->filename);
      *return_status = DB_RETURN_ERROR;
      return NULL;
    }
    
    if( load_next_fs_SequenceDB(sdb) == FALSE ) {
      warn("On file source [%d] [%s] could not open the file",sdb->current_source+1,sdb->fs[sdb->current_source+1]->filename);
      *return_status = DB_RETURN_ERROR;
      return NULL;
    }
    
    
    if( (out = get_next_SequenceDB(sdb)) != NULL ) {
      *return_status = DB_RETURN_OK;
      return out;
    }
    count++;
    warn("Ok, don't like this, just loaded the next Filesource, and got no sequence. Nope!");

    if( SequenceDB_at_end(sdb) == TRUE ) {
      *return_status = DB_RETURN_END;
      return NULL;
    }

    if( count > 10 ) {
      /*** break out of infinite loop ***/

      warn("Too many failed reloads of databases, going to fail");
      *return_status = DB_RETURN_ERROR;
      return NULL;
    }


  } /** back for for(;;) **/

}
    
/* Function:  get_next_SequenceDB(sdb)
 *
 * Descrip:    Main switch around formats    
 *
 *
 * Arg:        sdb [UNKN ] Undocumented argument [SequenceDB *]
 *
 * Return [UNKN ]  Undocumented return value [Sequence *]
 *
 */
# line 290 "sequencedb.dy"
Sequence * get_next_SequenceDB(SequenceDB * sdb)
{
  
  /* remember the byte position now */

  sdb->byte_position = ftell(sdb->current_file);

  switch (sdb->fs[sdb->current_source]->format) {
  case SEQ_DB_FASTA : 
    return read_fasta_Sequence(sdb->current_file);
  default :
    warn("Unknown SequenceDB type [%d]",sdb->fs[sdb->current_source]->format);
    return NULL;
  }
}


/* Function:  SequenceDB_at_end(sdb)
 *
 * Descrip:    Tells you if the SequenceDB is actually ended
 *             in terms of no more FileSources to eat through
 *
 *
 * Arg:        sdb [UNKN ] Undocumented argument [SequenceDB *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 312 "sequencedb.dy"
boolean SequenceDB_at_end(SequenceDB * sdb)
{
  if( sdb->current_source == -1 ) {
    warn("Bad bug: asking when it has finished when you have not init'd seqdb %s",sdb->name);
    return TRUE;
  }

  if( sdb->current_source+1 < sdb->len ) {
    return FALSE;
  }

  return TRUE;
}


/* Function:  load_next_fs_SequenceDB(sdb)
 *
 * Descrip:    Opens or attaches next FileSource stream
 *
 *             Does not close anything - use /close_last_fs_SequenceDB
 *
 *
 * Arg:        sdb [UNKN ] Undocumented argument [SequenceDB *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 333 "sequencedb.dy"
boolean load_next_fs_SequenceDB(SequenceDB * sdb)
{
  FileSource * fs;

  if( sdb->current_source == -1 ) {
    warn("Bad bug: trying to close last source when you have not init'd seqdb %s",sdb->name);
    return FALSE;
  }

  if( sdb->current_source >= sdb->len ) {
    warn("Bad bug. Someone is trying to load the next fs file when there are none (has not tested with SequenceDB_at_end...). So. I will fail, but database is actually at the end");
    return FALSE;
  }
  

  fs = sdb->fs[sdb->current_source];

  if( fs->filename != NULL ) {
    if( (sdb->current_file = openfile(fs->filename,"r")) == NULL ) {
      warn("Could not open file [%s] for database [%s]",fs->filename,sdb->name);
      return FALSE;
    }
  } else {
    sdb->current_file = fs->input;
  }

  return TRUE;
}

/* Function:  close_last_fs_SequenceDB(sdb)
 *
 * Descrip:    closes the last FileSource: checks if it was a straight stream
 *             (in which case does not close)
 *
 *
 * Arg:        sdb [UNKN ] Undocumented argument [SequenceDB *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 367 "sequencedb.dy"
boolean close_last_fs_SequenceDB(SequenceDB * sdb)
{
  FileSource * fs;

  if( sdb->current_source == -1 ) {
    warn("Bad bug: trying to close last source when you have not init'd seqdb %s",sdb->name);
    return FALSE;
  }

  fs = sdb->fs[sdb->current_source];

  if( fs->filename != NULL ) {
    fclose(sdb->current_file);
  } else if( fs->input != NULL ) {
    warn("Can't handle closes on streams yet. Not sure what to do!");
  }
    
  

  sdb->current_source++;

  return TRUE;
}


  /*** I/O ****/

/* Function:  SequenceDB_from_FILE_and_format(input,format)
 *
 * Descrip:    makes a SequencDB from a straight file stream.
 *
 *             This means SequenceDB will *not* close it when
 *             the SequenceDB is closed.
 *
 *
 * Arg:         input [READ ] filestream [FILE *]
 * Arg:        format [UNKN ] format as defined by /word_to_format [int]
 *
 * Return [UNKN ]  Undocumented return value [SequenceDB *]
 *
 */
# line 403 "sequencedb.dy"
SequenceDB * SequenceDB_from_FILE_and_format(FILE * input,int format)
{
  SequenceDB * out;
  FileSource * fs;
  
  out = SequenceDB_alloc_len(1);

  fs = FileSource_from_FILE_and_format(input,format);

  add_SequenceDB(out,fs);

  return out;

}

/* Function:  FileSource_from_FILE_and_format(input,format)
 *
 * Descrip:    Makes a file source from a straigth stream
 *
 *
 * Arg:         input [UNKN ] Undocumented argument [FILE *]
 * Arg:        format [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [FileSource *]
 *
 */
# line 422 "sequencedb.dy"
FileSource * FileSource_from_FILE_and_format(FILE * input,int format)
{
  FileSource * fs;

  fs = FileSource_alloc();

  fs->input = input;
  fs->format = format;

  return fs;
}

/* Function:  single_fasta_SequenceDB(filename)
 *
 * Descrip:    pre-packed single fasta file db
 *
 *
 *
 * Arg:        filename [UNKN ] name of fastadb [char *]
 *
 * Return [UNKN ]  Undocumented return value [SequenceDB *]
 *
 */
# line 440 "sequencedb.dy"
SequenceDB * single_fasta_SequenceDB(char * filename)
{
  SequenceDB * out;
  FileSource * fs;
  
  if( touchfile(filename) == FALSE) {
    warn("Cannot make SequenceDB from an unopenable fileanme [%s]",filename);
    return NULL;
  }

  fs = FileSource_alloc();
  fs->filename = stringalloc(filename);
  fs->format   = SEQ_DB_FASTA;

  out = SequenceDB_alloc_len(1);
  out->seq_start = -1;
  out->seq_end   = -1;
  add_SequenceDB(out,fs);

  return out;
}
  
  

/* Function:  read_SequenceDB_line(line,ifp)
 *
 * Descrip:    Reads a SequenceDB definition from
 *
 *             seqdb <name>
 *             <filename> <format> <type>
 *             ...
 *             endseqdb
 *
 *
 *
 * Arg:        line [UNKN ] starting line (seqdb line) [char *]
 * Arg:         ifp [UNKN ] file input [FILE *]
 *
 * Return [UNKN ]  Undocumented return value [SequenceDB *]
 *
 */
# line 476 "sequencedb.dy"
SequenceDB * read_SequenceDB_line(char * line,FILE * ifp)
{
  SequenceDB * out = NULL;
  FileSource * fs;
  char buffer[MAXLINE];
  char * runner;


  if( strstartcmp(line,"seqdb") != 0 ) {
    warn("Attempting to read a sequence line without a seqdb start");
    return NULL;
  }

  runner = strtok(line,spacestr);
  runner = strtok(line,spacestr);

  if( runner == NULL ) {
    out->name = stringalloc("UnNamedDatabase");
  }
  else out->name = stringalloc(runner);


  out = SequenceDB_alloc_std();

  while( fgets(buffer,MAXLINE,ifp) != NULL  ){
    if( strstartcmp(buffer,"#") == 0 ) 
      continue;
    if( strstartcmp(buffer,"end") == 0 )
      break;
    fs = FileSource_from_line(buffer);
    if( fs != NULL )
      add_SequenceDB(out,fs);
  }

  return out;
}

/* Function:  word_to_format(word)
 *
 * Descrip:    converts char * to format for SequenceDB FileSources
 *
 *
 * Arg:        word [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
# line 516 "sequencedb.dy"
int word_to_format(char * word)
{
  if( strcmp(word,"fasta") == 0 ) {
    return SEQ_DB_FASTA;
  }

  return SEQ_DB_UNKNOWN;
}

    
/* Function:  FileSource_from_line(line)
 *
 * Descrip:    Reads line
 *             filename format type
 *
 *             where format is determined by /word_to_format
 *             and type is protein/dna
 *
 *
 * Arg:        line [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [FileSource *]
 *
 */
# line 534 "sequencedb.dy"
FileSource * FileSource_from_line(char * line)
{
  FileSource * out;
  char * runner;
  char * run2;
  char * run3;

  runner = strtok(line,spacestr);
  run2 = strtok(line,spacestr);
  run3 = strtok(line,spacestr);

  if( runner == NULL || run2 == NULL || run3 == NULL ) {
    warn("You have not provided a database source line");
    return NULL;
  }

 

  out = FileSource_alloc();

  out->filename = stringalloc(runner);

  if( (out->format = word_to_format(run2)) == SEQ_DB_UNKNOWN) {
    warn("For filename %s, the format [%s] is unknown to me",runner,run2);
  }

  return out;
}




# line 579 "sequencedb.c"
/* Function:  hard_link_FileSource(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [FileSource *]
 *
 * Return [UNKN ]  Undocumented return value [FileSource *]
 *
 */
FileSource * hard_link_FileSource(FileSource * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a FileSource object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  FileSource_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [FileSource *]
 *
 */
FileSource * FileSource_alloc(void) 
{
    FileSource * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(FileSource *) ckalloc (sizeof(FileSource))) == NULL)    {  
      warn("FileSource_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->filename = NULL;    
    out->format = 0; 
    out->type = 0;   


    return out;  
}    


/* Function:  free_FileSource(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [FileSource *]
 *
 * Return [UNKN ]  Undocumented return value [FileSource *]
 *
 */
FileSource * free_FileSource(FileSource * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a FileSource obj. Should be trappable");    
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
    if( obj->filename != NULL)   
      ckfree(obj->filename);     
    /* obj->input is linked in */ 


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_SequenceDB(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_SequenceDB
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [FileSource **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_SequenceDB(FileSource ** list,int i,int j)  
{
    FileSource * temp;   
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_SequenceDB(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_SequenceDB which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [FileSource **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_SequenceDB(FileSource ** list,int left,int right,int (*comp)(FileSource * ,FileSource * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_SequenceDB(list,left,(left+right)/2);   
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_SequenceDB (list,++last,i); 
      }  
    swap_SequenceDB (list,left,last);    
    qsort_SequenceDB(list,left,last-1,comp); 
    qsort_SequenceDB(list,last+1,right,comp);    
}    


/* Function:  sort_SequenceDB(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_SequenceDB
 *
 *
 * Arg:         obj [UNKN ] Object containing list [SequenceDB *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_SequenceDB(SequenceDB * obj,int (*comp)(FileSource *, FileSource *)) 
{
    qsort_SequenceDB(obj->fs,0,obj->len-1,comp); 
    return;  
}    


/* Function:  expand_SequenceDB(obj,len)
 *
 * Descrip:    Really an internal function for add_SequenceDB
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [SequenceDB *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_SequenceDB(SequenceDB * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_SequenceDB called with no need"); 
      return TRUE;   
      }  


    if( (obj->fs = (FileSource ** ) ckrealloc (obj->fs,sizeof(FileSource *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_SequenceDB, returning FALSE");   
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_SequenceDB(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [SequenceDB *]
 * Arg:        add [OWNER] Object to add to the list [FileSource *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_SequenceDB(SequenceDB * obj,FileSource * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_SequenceDB(obj,obj->len + SequenceDBLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->fs[obj->len++]=add; 
    return TRUE; 
}    


/* Function:  flush_SequenceDB(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [SequenceDB *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_SequenceDB(SequenceDB * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->fs[i] != NULL)    {  
        free_FileSource(obj->fs[i]); 
        obj->fs[i] = NULL;   
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  SequenceDB_alloc_std(void)
 *
 * Descrip:    Equivalent to SequenceDB_alloc_len(SequenceDBLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SequenceDB *]
 *
 */
SequenceDB * SequenceDB_alloc_std(void) 
{
    return SequenceDB_alloc_len(SequenceDBLISTLENGTH);   
}    


/* Function:  SequenceDB_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [SequenceDB *]
 *
 */
SequenceDB * SequenceDB_alloc_len(int len) 
{
    SequenceDB * out;   /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = SequenceDB_alloc()) == NULL)   
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->fs = (FileSource ** ) ckcalloc (len,sizeof(FileSource *))) == NULL) {  
      warn("Warning, ckcalloc failed in SequenceDB_alloc_len");  
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_SequenceDB(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [SequenceDB *]
 *
 * Return [UNKN ]  Undocumented return value [SequenceDB *]
 *
 */
SequenceDB * hard_link_SequenceDB(SequenceDB * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a SequenceDB object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  SequenceDB_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SequenceDB *]
 *
 */
SequenceDB * SequenceDB_alloc(void) 
{
    SequenceDB * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(SequenceDB *) ckalloc (sizeof(SequenceDB))) == NULL)    {  
      warn("SequenceDB_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->name = NULL;    
    out->fs = NULL;  
    out->len = out->maxlen = 0;  
    out->current_source = -1;    
    out->sequence_no = 0;    
    out->byte_position = 0;  
    out->has_warned_single = 0;  
    out->seq_start = -1; 
    out->seq_end = -1;   


    return out;  
}    


/* Function:  free_SequenceDB(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [SequenceDB *]
 *
 * Return [UNKN ]  Undocumented return value [SequenceDB *]
 *
 */
SequenceDB * free_SequenceDB(SequenceDB * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a SequenceDB obj. Should be trappable");    
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
    if( obj->name != NULL)   
      ckfree(obj->name);     
    if( obj->fs != NULL) {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->fs[i] != NULL)  
          free_FileSource(obj->fs[i]);   
        }  
      ckfree(obj->fs);   
      }  
    /* obj->current_file is linked in */ 


    ckfree(obj); 
    return NULL; 
}    


/* Function:  replace_name_SequenceDB(obj,name)
 *
 * Descrip:    Replace member variable name
 *             For use principly by API functions
 *
 *
 * Arg:         obj [UNKN ] Object holding the variable [SequenceDB *]
 * Arg:        name [OWNER] New value of the variable [char *]
 *
 * Return [SOFT ]  member variable name [boolean]
 *
 */
boolean replace_name_SequenceDB(SequenceDB * obj,char * name) 
{
    if( obj == NULL)     {  
      warn("In replacement function name for object SequenceDB, got a NULL object"); 
      return FALSE;  
      }  
    obj->name = name;    
    return TRUE; 
}    


/* Function:  access_name_SequenceDB(obj)
 *
 * Descrip:    Access member variable name
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [SequenceDB *]
 *
 * Return [SOFT ]  member variable name [char *]
 *
 */
char * access_name_SequenceDB(SequenceDB * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function name for object SequenceDB, got a NULL object");    
      return NULL;   
      }  
    return obj->name;    
}    


/* Function:  access_fs_SequenceDB(obj,i)
 *
 * Descrip:    Access members stored in the fs list
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the list [SequenceDB *]
 * Arg:          i [UNKN ] Position in the list [int]
 *
 * Return [SOFT ]  Element of the list [FileSource *]
 *
 */
FileSource * access_fs_SequenceDB(SequenceDB * obj,int i) 
{
    if( obj == NULL)     {  
      warn("In accessor function fs for object SequenceDB, got a NULL object");  
      return NULL;   
      }  
    if( obj->len <= i )  {  
      warn("In accessor function fs for object SequenceDB, index %%d is greater than list length %%d",i,obj->len);   
      return NULL;   
      }  
    return obj->fs[i];   
}    


/* Function:  length_fs_SequenceDB(obj)
 *
 * Descrip:    discover the length of the list
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the list [SequenceDB *]
 *
 * Return [UNKN ]  length of the list [int]
 *
 */
int length_fs_SequenceDB(SequenceDB * obj) 
{
    if( obj == NULL)     {  
      warn("In length function fs for object SequenceDB, got a NULL object");    
      return -1;     
      }  
    return obj->len;     
}    


/* Function:  replace_current_source_SequenceDB(obj,current_source)
 *
 * Descrip:    Replace member variable current_source
 *             For use principly by API functions
 *
 *
 * Arg:                   obj [UNKN ] Object holding the variable [SequenceDB *]
 * Arg:        current_source [OWNER] New value of the variable [int]
 *
 * Return [SOFT ]  member variable current_source [boolean]
 *
 */
boolean replace_current_source_SequenceDB(SequenceDB * obj,int current_source) 
{
    if( obj == NULL)     {  
      warn("In replacement function current_source for object SequenceDB, got a NULL object");   
      return FALSE;  
      }  
    obj->current_source = current_source;    
    return TRUE; 
}    


/* Function:  access_current_source_SequenceDB(obj)
 *
 * Descrip:    Access member variable current_source
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [SequenceDB *]
 *
 * Return [SOFT ]  member variable current_source [int]
 *
 */
int access_current_source_SequenceDB(SequenceDB * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function current_source for object SequenceDB, got a NULL object");  
      return 0;  
      }  
    return obj->current_source;  
}    


/* Function:  replace_current_file_SequenceDB(obj,current_file)
 *
 * Descrip:    Replace member variable current_file
 *             For use principly by API functions
 *
 *
 * Arg:                 obj [UNKN ] Object holding the variable [SequenceDB *]
 * Arg:        current_file [OWNER] New value of the variable [FILE *]
 *
 * Return [SOFT ]  member variable current_file [boolean]
 *
 */
boolean replace_current_file_SequenceDB(SequenceDB * obj,FILE * current_file) 
{
    if( obj == NULL)     {  
      warn("In replacement function current_file for object SequenceDB, got a NULL object"); 
      return FALSE;  
      }  
    obj->current_file = current_file;    
    return TRUE; 
}    


/* Function:  access_current_file_SequenceDB(obj)
 *
 * Descrip:    Access member variable current_file
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [SequenceDB *]
 *
 * Return [SOFT ]  member variable current_file [FILE *]
 *
 */
FILE * access_current_file_SequenceDB(SequenceDB * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function current_file for object SequenceDB, got a NULL object");    
      return NULL;   
      }  
    return obj->current_file;    
}    


/* Function:  replace_sequence_no_SequenceDB(obj,sequence_no)
 *
 * Descrip:    Replace member variable sequence_no
 *             For use principly by API functions
 *
 *
 * Arg:                obj [UNKN ] Object holding the variable [SequenceDB *]
 * Arg:        sequence_no [OWNER] New value of the variable [int]
 *
 * Return [SOFT ]  member variable sequence_no [boolean]
 *
 */
boolean replace_sequence_no_SequenceDB(SequenceDB * obj,int sequence_no) 
{
    if( obj == NULL)     {  
      warn("In replacement function sequence_no for object SequenceDB, got a NULL object");  
      return FALSE;  
      }  
    obj->sequence_no = sequence_no;  
    return TRUE; 
}    


/* Function:  access_sequence_no_SequenceDB(obj)
 *
 * Descrip:    Access member variable sequence_no
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [SequenceDB *]
 *
 * Return [SOFT ]  member variable sequence_no [int]
 *
 */
int access_sequence_no_SequenceDB(SequenceDB * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function sequence_no for object SequenceDB, got a NULL object"); 
      return 0;  
      }  
    return obj->sequence_no;     
}    


/* Function:  replace_byte_position_SequenceDB(obj,byte_position)
 *
 * Descrip:    Replace member variable byte_position
 *             For use principly by API functions
 *
 *
 * Arg:                  obj [UNKN ] Object holding the variable [SequenceDB *]
 * Arg:        byte_position [OWNER] New value of the variable [int]
 *
 * Return [SOFT ]  member variable byte_position [boolean]
 *
 */
boolean replace_byte_position_SequenceDB(SequenceDB * obj,int byte_position) 
{
    if( obj == NULL)     {  
      warn("In replacement function byte_position for object SequenceDB, got a NULL object");    
      return FALSE;  
      }  
    obj->byte_position = byte_position;  
    return TRUE; 
}    


/* Function:  access_byte_position_SequenceDB(obj)
 *
 * Descrip:    Access member variable byte_position
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [SequenceDB *]
 *
 * Return [SOFT ]  member variable byte_position [int]
 *
 */
int access_byte_position_SequenceDB(SequenceDB * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function byte_position for object SequenceDB, got a NULL object");   
      return 0;  
      }  
    return obj->byte_position;   
}    


/* Function:  replace_filename_FileSource(obj,filename)
 *
 * Descrip:    Replace member variable filename
 *             For use principly by API functions
 *
 *
 * Arg:             obj [UNKN ] Object holding the variable [FileSource *]
 * Arg:        filename [OWNER] New value of the variable [char *]
 *
 * Return [SOFT ]  member variable filename [boolean]
 *
 */
boolean replace_filename_FileSource(FileSource * obj,char * filename) 
{
    if( obj == NULL)     {  
      warn("In replacement function filename for object FileSource, got a NULL object"); 
      return FALSE;  
      }  
    obj->filename = filename;    
    return TRUE; 
}    


/* Function:  access_filename_FileSource(obj)
 *
 * Descrip:    Access member variable filename
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [FileSource *]
 *
 * Return [SOFT ]  member variable filename [char *]
 *
 */
char * access_filename_FileSource(FileSource * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function filename for object FileSource, got a NULL object");    
      return NULL;   
      }  
    return obj->filename;    
}    


/* Function:  replace_input_FileSource(obj,input)
 *
 * Descrip:    Replace member variable input
 *             For use principly by API functions
 *
 *
 * Arg:          obj [UNKN ] Object holding the variable [FileSource *]
 * Arg:        input [OWNER] New value of the variable [FILE *]
 *
 * Return [SOFT ]  member variable input [boolean]
 *
 */
boolean replace_input_FileSource(FileSource * obj,FILE * input) 
{
    if( obj == NULL)     {  
      warn("In replacement function input for object FileSource, got a NULL object");    
      return FALSE;  
      }  
    obj->input = input;  
    return TRUE; 
}    


/* Function:  access_input_FileSource(obj)
 *
 * Descrip:    Access member variable input
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [FileSource *]
 *
 * Return [SOFT ]  member variable input [FILE *]
 *
 */
FILE * access_input_FileSource(FileSource * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function input for object FileSource, got a NULL object");   
      return NULL;   
      }  
    return obj->input;   
}    


/* Function:  replace_format_FileSource(obj,format)
 *
 * Descrip:    Replace member variable format
 *             For use principly by API functions
 *
 *
 * Arg:           obj [UNKN ] Object holding the variable [FileSource *]
 * Arg:        format [OWNER] New value of the variable [int]
 *
 * Return [SOFT ]  member variable format [boolean]
 *
 */
boolean replace_format_FileSource(FileSource * obj,int format) 
{
    if( obj == NULL)     {  
      warn("In replacement function format for object FileSource, got a NULL object");   
      return FALSE;  
      }  
    obj->format = format;    
    return TRUE; 
}    


/* Function:  access_format_FileSource(obj)
 *
 * Descrip:    Access member variable format
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [FileSource *]
 *
 * Return [SOFT ]  member variable format [int]
 *
 */
int access_format_FileSource(FileSource * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function format for object FileSource, got a NULL object");  
      return 0;  
      }  
    return obj->format;  
}    


/* Function:  replace_type_FileSource(obj,type)
 *
 * Descrip:    Replace member variable type
 *             For use principly by API functions
 *
 *
 * Arg:         obj [UNKN ] Object holding the variable [FileSource *]
 * Arg:        type [OWNER] New value of the variable [int]
 *
 * Return [SOFT ]  member variable type [boolean]
 *
 */
boolean replace_type_FileSource(FileSource * obj,int type) 
{
    if( obj == NULL)     {  
      warn("In replacement function type for object FileSource, got a NULL object"); 
      return FALSE;  
      }  
    obj->type = type;    
    return TRUE; 
}    


/* Function:  access_type_FileSource(obj)
 *
 * Descrip:    Access member variable type
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [FileSource *]
 *
 * Return [SOFT ]  member variable type [int]
 *
 */
int access_type_FileSource(FileSource * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function type for object FileSource, got a NULL object");    
      return 0;  
      }  
    return obj->type;    
}    



#ifdef _cplusplus
}
#endif
