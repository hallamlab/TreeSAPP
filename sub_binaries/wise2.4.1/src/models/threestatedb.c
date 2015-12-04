#ifdef _cplusplus
extern "C" {
#endif
#include "threestatedb.h"

/* Function:  open_for_indexing_ThreeStateDB(mdb)
 *
 * Descrip:    opens database ready for calls to 
 *             /indexed_ThreeStateModel_ThreeStateDB
 *
 *
 *
 * Arg:        mdb [UNKN ] database to open for index calls [ThreeStateDB *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 92 "threestatedb.dy"
boolean open_for_indexing_ThreeStateDB(ThreeStateDB * mdb)
{
  switch(mdb->dbtype) {
  case TSMDB_SINGLE :
    return TRUE;
  case TSMDB_HMMER1PFAM :
    return TRUE; /* name is the index! */
  case TSMDB_PROTEIN :
    return TRUE; /* sequence db's don't need opening for indexing */
  case TSMDB_GENERIC :
    return ((*mdb->open_index_generic)(mdb));
  default :
    warn("Unknown threestatemodel db type.");
    return FALSE;
  }

}

/* Function:  close_for_indexing_ThreeStateDB(mdb)
 *
 * Descrip:    closes an indexable database
 *
 *
 * Arg:        mdb [UNKN ] Undocumented argument [ThreeStateDB *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 113 "threestatedb.dy"
boolean close_for_indexing_ThreeStateDB(ThreeStateDB * mdb)
{
  switch(mdb->dbtype) {
  case TSMDB_GENERIC :
    return ((*mdb->close_index_generic)(mdb));
  default :
    return TRUE;
  }
}

/* Function:  set_search_type_ThreeStateDB(tdb,type)
 *
 * Descrip:    Set the search type of this threestatedb...
 *
 *
 * Arg:         tdb [UNKN ] Undocumented argument [ThreeStateDB *]
 * Arg:        type [UNKN ] to set can be any of the modes found in threestatemodel [char *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 128 "threestatedb.dy"
boolean set_search_type_ThreeStateDB(ThreeStateDB * tdb,char * type)
{
  int ret;

  ret = threestatemodel_mode_from_string(type);

  if( ret == TSM_unknown ) {
    /* warning already issued */
    return FALSE;
  }

  tdb->type = ret;
  return TRUE;

}

  
/* Function:  indexed_ThreeStateModel_ThreeStateDB(mdb,en)
 *
 * Descrip:    Retrieves a model from a database which has been opened
 *             for indexing by /open_for_indexing_ThreeStateDB
 *
 *             The index information comes from the dataentry which should 
 *             have been from a search of the ThreeStateDB.
 *
 *
 * Arg:        mdb [UNKN ] database where this is indexed [ThreeStateDB *]
 * Arg:         en [UNKN ] dataentry to pull the model from [DataEntry *]
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateModel *]
 *
 */
# line 156 "threestatedb.dy"
ThreeStateModel * indexed_ThreeStateModel_ThreeStateDB(ThreeStateDB * mdb,DataEntry * en)
{
  Sequence * seq;
  Protein * pro;
  ThreeStateModel * tsm;


  switch(mdb->dbtype) {
  case TSMDB_SINGLE :
    return hard_link_ThreeStateModel(mdb->single);
  case TSMDB_HMMER1PFAM :
    tsm = ThreeStateModel_from_name_PfamHmmer1DB(mdb->phdb,en->name);
    set_startend_policy_ThreeStateModel(tsm,mdb->type,30,0.2);
    return tsm;

  case TSMDB_PROTEIN :
    seq = get_Sequence_from_SequenceDB(mdb->sdb,en);
    if( seq == NULL ) {
      warn("could not retrieve %s as a sequence from database",en->name);
      return NULL;
    }

    pro = Protein_from_Sequence(seq);

    if( pro == NULL ) {
      warn("Could not convert sequence to a protein. Exiting!");
      return NULL;
    }

    /* convert protein to threestatemodel */

    tsm = ThreeStateModel_from_half_bit_Sequence(pro,mdb->comp,mdb->rm,mdb->gap,mdb->ext);

    if( tsm == NULL ) {
      warn("Could not convert protein to threestatemode. Exiting!");
      free_Protein(pro);
      return NULL;
    }

    free_Protein(pro);
    /* DB status already set by seqdb */
    return tsm;
  case TSMDB_GENERIC :
    tsm = ((*mdb->index_generic)(mdb,en));
    if( tsm == NULL ) {
      return NULL;
    }
    /*   fprintf(stdout,"Setting %d as policy\n",mdb->type); */
    set_startend_policy_ThreeStateModel(tsm,mdb->type,30,0.2);

    return tsm;
    
  default : 
    warn("Unknown threestatedb type");
    return NULL;
  }

  warn("Should never get here - in threestatedb reload!");

  return NULL;

}

/* Function:  new_proteindb_ThreeStateDB(sdb,comp,gap,ext)
 *
 * Descrip:    makes a new ThreeStateDB from a
 *             sequencedb (better be protein!)
 *
 *
 *
 * Arg:         sdb [READ ] sequence database to use [SequenceDB *]
 * Arg:        comp [READ ] comparison matrix to use [CompMat *]
 * Arg:         gap [READ ] gap open penalty [int]
 * Arg:         ext [READ ] gap extensions penalty [int]
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateDB *]
 *
 */
# line 229 "threestatedb.dy"
ThreeStateDB * new_proteindb_ThreeStateDB(SequenceDB * sdb,CompMat * comp,int gap,int ext)
{
  ThreeStateDB * out;

  out = ThreeStateDB_alloc();
  out->sdb  =  hard_link_SequenceDB(sdb);
  out->comp =  hard_link_CompMat(comp);
  out->gap  = gap;
  out->ext  = ext;
  out->rm   = default_RandomModel();
  out->dbtype = TSMDB_PROTEIN;
  return out;

}
  
/* Function:  new_single_ThreeStateDB(tsm,rm)
 *
 * Descrip:    Making a new ThreeStateDB from a single
 *             model
 *
 *
 *
 * Arg:        tsm [READ ] a single ThreeStateModel [ThreeStateModel *]
 * Arg:         rm [READ ] random model to be used in comparisons.. [RandomModel *]
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateDB *]
 *
 */
# line 252 "threestatedb.dy"
ThreeStateDB * new_single_ThreeStateDB(ThreeStateModel * tsm,RandomModel * rm)
{
  ThreeStateDB * out;

  out = ThreeStateDB_alloc();
  out->single = hard_link_ThreeStateModel(tsm);
  if( tsm->rm == NULL ) 
    out->rm = hard_link_RandomModel(rm);
  else
    out->rm = hard_link_RandomModel(tsm->rm);
  out->dbtype = TSMDB_SINGLE;
  return out;

}

/* Function:  new_PfamHmmer1DB_ThreeStateDB(dirname)
 *
 * Descrip:    Makes a new PfamHmmer1DB from a filename
 *             indicating the directory
 *
 *
 * Arg:        dirname [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateDB *]
 *
 */
# line 271 "threestatedb.dy"
ThreeStateDB * new_PfamHmmer1DB_ThreeStateDB(char * dirname)
{
  ThreeStateDB * out;

  out = ThreeStateDB_alloc();

  out->phdb = PfamHmmer1DB_from_dirname(dirname);
  out->dbtype = TSMDB_HMMER1PFAM;
  
  return out;
}

  
/* Function:  dataentry_add_ThreeStateDB(de,tss,mdb)
 *
 * Descrip:    This function adds the internal entry information 
 *             (eg indexing point) into the dataentry
 *
 *
 * Arg:         de [UNKN ] Undocumented argument [DataEntry *]
 * Arg:        tss [UNKN ] Undocumented argument [ThreeStateScore *]
 * Arg:        mdb [UNKN ] Undocumented argument [ThreeStateDB *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 288 "threestatedb.dy"
boolean dataentry_add_ThreeStateDB(DataEntry * de,ThreeStateScore * tss,ThreeStateDB * mdb)
{

  switch(mdb->dbtype) {
  case TSMDB_SINGLE :
    de->name = stringalloc(mdb->single->name);
    return TRUE;
  case TSMDB_HMMER1PFAM :
    if( tss == NULL ) {
    } else {
      de->name = stringalloc(tss->name);
    }
    return TRUE;
  case TSMDB_PROTEIN :
    add_SequenceDB_info_DataEntry(mdb->sdb,de);
    return TRUE;
  case TSMDB_GENERIC :
    if( (*mdb->dataentry_add)(mdb,de) == FALSE ) {
      warn("Could not add dataentry info to the entry %s",tss->name);
      return FALSE;
    } else {
      return TRUE;
    }
    
  default : 
    warn("Unknown threestatedb type");
    return FALSE;
  }


  return TRUE;
}


/* Function:  open_ThreeStateDB(mdb)
 *
 * Descrip:    Open function for ThreeStateDB.
 *             An internal for this file but also
 *             used by, for example, GeneWiseDB that
 *             wants to get at the underlying models, 
 *             not the log-odds.
 *
 *
 * Arg:        mdb [UNKN ] Undocumented argument [ThreeStateDB *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 329 "threestatedb.dy"
boolean open_ThreeStateDB(ThreeStateDB * mdb)
{
  int ret;
  int return_status;
  ThreeStateModel * temp;
  int count;

  mdb->current_no = 0;

  switch( mdb->dbtype ) {

  case TSMDB_SINGLE :
    return TRUE; /* should be fine! */
  case TSMDB_HMMER1PFAM :
    if( mdb->phdb == NULL ) {
      warn("No hmmer1 db to open for threestatedb!");
      return FALSE;
    }
    mdb->phdb->cur = 0;
    break;
  case TSMDB_PROTEIN :
    if( mdb->sdb == NULL ) {
      warn("Attempting to open a protein tsm with no sequence db!");
      return FALSE;
    } 
    mdb->seq_cache = init_SequenceDB(mdb->sdb,&ret);
    if( ret == DB_RETURN_ERROR ) {
      return FALSE;
    }
    if( ret == DB_RETURN_END ) {
      warn("Due to some bad coding, can't cope with single protein databases in tsmdbs. oooops!");
    }
    break;
  case TSMDB_GENERIC :
    ((*mdb->open_generic)(mdb));
    break;
  default :
    warn("Got an unrecognisable tsm db type in opening tsm %d",mdb->dbtype);
    return FALSE;
  }
    
  if( mdb->hmm_model_start != -1 && mdb->hmm_model_end != -1 ) {
    for(count=0;count<mdb->hmm_model_start;count++) {
      temp = read_TSM_ThreeStateDB(mdb,&return_status);
      free_ThreeStateModel(temp);
    }
  }


  return TRUE;
}


/* Function:  read_TSM_ThreeStateDB(mdb,return_status)
 *
 * Descrip:    Reads a threestatemodel out from the 
 *             database. People will probably want the
 *             ThreeStateScore *not* the model, but some
 *             systems will want the model.
 *
 *
 * Arg:                  mdb [UNKN ] Undocumented argument [ThreeStateDB *]
 * Arg:        return_status [UNKN ] Undocumented argument [int *]
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateModel *]
 *
 */
# line 388 "threestatedb.dy"
ThreeStateModel * read_TSM_ThreeStateDB(ThreeStateDB * mdb,int * return_status)
{
  ThreeStateModel * tsm;
  Protein * pro;
  Sequence * seq;

  if( mdb->hmm_model_end != -1 && mdb->current_no == mdb->hmm_model_end ) {
    *return_status = DB_RETURN_END;
    return NULL;
  }

  mdb->current_no++;

  switch( mdb->dbtype ) {

  case TSMDB_SINGLE :
    *return_status = DB_RETURN_END;
    if( mdb->single->rm == NULL ) {
      warn("Threestate model without an internal random model!");  
      mdb->single->rm = hard_link_RandomModel(mdb->rm);
    }

    return hard_link_ThreeStateModel(mdb->single);
  case TSMDB_HMMER1PFAM :
    tsm= read_next_TSM_PfamHmmer1DB(mdb->phdb,return_status);
    set_startend_policy_ThreeStateModel(tsm,mdb->type,30,0.2);
    return tsm;

  case TSMDB_PROTEIN :
    if( mdb->seq_cache != NULL ) {
      /* just after an open. Should actually use this sequence, and flush the cache */
      pro = Protein_from_Sequence(hard_link_Sequence(mdb->seq_cache));
      mdb->seq_cache = free_Sequence(mdb->seq_cache);
      *return_status = DB_RETURN_OK;
    } else {

      /* reload a sequence from a database */
      seq = reload_SequenceDB(NULL,mdb->sdb,return_status);

      /* exit now if error */
      if( *return_status == DB_RETURN_ERROR ) {
	return NULL; /* might have leaked memory. Ugh! */
      } 

      /* if we get NULL... for the moment, silent flag end */

      if( seq == NULL ) {
	*return_status = DB_RETURN_END;
	return NULL;
      }

      pro = Protein_from_Sequence(seq);
    }
    if( pro == NULL ) {
      warn("Could not convert sequence to a protein. Exiting!");
      *return_status = DB_RETURN_ERROR;
      return NULL;
    }

    /* convert protein to threestatemodel */

    tsm = ThreeStateModel_from_half_bit_Sequence(pro,mdb->comp,mdb->rm,mdb->gap,mdb->ext);

    if( tsm == NULL ) {
      warn("Could not convert protein to threestatemode. Exiting!");
      free_Protein(pro);
      *return_status = DB_RETURN_ERROR;
      return NULL;
    }

    /* DB status already set by seqdb */
    return tsm;
  case TSMDB_GENERIC :
    tsm =  ((*mdb->reload_generic)(mdb,return_status));
    if( tsm == NULL ) {
      return NULL; /* means end of database */
    }
    set_startend_policy_ThreeStateModel(tsm,mdb->type,30,0.2);
    return tsm;

  default :
    warn("Got an unrecognisable tsm db type in read-load");
    return NULL;
  }


}
 
  
/* Function:  init_ThreeStateDB(mdb,return_status)
 *
 * Descrip:    Init function for ThreeStateDB
 *
 *             Is going to open file, read first model, complain if
 *             NULL, and convert to a score system.
 *
 *
 * Arg:                  mdb [RW   ] Model database [ThreeStateDB *]
 * Arg:        return_status [WRITE] return from database.h system [int *]
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateScore *]
 *
 */
# line 486 "threestatedb.dy"
ThreeStateScore * init_ThreeStateDB(ThreeStateDB * mdb,int * return_status)
{
  ThreeStateModel * tsm;
  ThreeStateScore * tss;


  if( open_ThreeStateDB(mdb) == FALSE) {
    warn("Could not open ThreeStateDB, hence could not init it!");
    *return_status = DB_RETURN_ERROR;
    return NULL;
  }

  tsm = read_TSM_ThreeStateDB(mdb,return_status);
  if( *return_status == DB_RETURN_ERROR) 
    return NULL;

  set_startend_policy_ThreeStateModel(tsm,mdb->type,30,0.2);

  fold_RandomModel_into_ThreeStateModel(tsm,mdb->rm);

  tss = ThreeStateScore_from_ThreeStateModel(tsm);

  free_ThreeStateModel(tsm);

  *return_status = DB_RETURN_OK;

  return tss;
}

/* Function:  reload_ThreeStateDB(prev,tss,mdb,return_status)
 *
 * Descrip:    reloads the ThreeStateDB.
 *
 *             Frees the previous score system (could recycle memory).
 *             Reads database. calls END if gets NULL from read_HMF_ThreeStateModel
 *
 *
 * Arg:                 prev [UNKN ] Undocumented argument [ThreeStateScore *]
 * Arg:                  tss [UNKN ] the previous score system [NullString]
 * Arg:                  mdb [UNKN ] model database system [ThreeStateDB *]
 * Arg:        return_status [WRITE] return from database.h system [int *]
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateScore *]
 *
 */
# line 525 "threestatedb.dy"
ThreeStateScore * reload_ThreeStateDB(ThreeStateScore * prev,ThreeStateDB * mdb,int * return_status)
{
  ThreeStateModel * tsm;
  ThreeStateScore * tss;

  free_ThreeStateScore(prev);

  if( mdb->dbtype == TSMDB_SINGLE ) {
    *return_status = DB_RETURN_END;
    return NULL;
  }


  tsm = read_TSM_ThreeStateDB(mdb,return_status);
  if( *return_status != DB_RETURN_OK) 
    return NULL;

  set_startend_policy_ThreeStateModel(tsm,mdb->type,30,0.2);

  fold_RandomModel_into_ThreeStateModel(tsm,mdb->rm);

  tss = ThreeStateScore_from_ThreeStateModel(tsm);

  free_ThreeStateModel(tsm);


  *return_status = DB_RETURN_OK;

  return tss;
}
  

/* Function:  close_ThreeStateDB(prev,mdb)
 *
 * Descrip:    closes ThreeStateDB
 *
 *             At the moment, only needs to free previous
 *             and close the file
 *
 *
 * Arg:        prev [UNKN ] the last ThreeStateScore to be freed [ThreeStateScore *]
 * Arg:         mdb [UNKN ] Model database [ThreeStateDB *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 566 "threestatedb.dy"
boolean close_ThreeStateDB(ThreeStateScore * prev,ThreeStateDB * mdb)
{

  if( prev != NULL )
    free_ThreeStateScore(prev);
  if( mdb == NULL ) {
    warn("Trying to close a NULL threestatedb - considering this an error!");
    return FALSE;
  }
  switch(mdb->dbtype) {
  case TSMDB_SINGLE :
    return TRUE;
  case TSMDB_HMMER1PFAM :
    fclose(mdb->current_file);
    mdb->current_file = NULL;
    return TRUE; /* name is the index! */
  case TSMDB_PROTEIN :
    close_SequenceDB(NULL,mdb->sdb);
    return TRUE; /* sequence db's don't need opening for indexing */
  case TSMDB_GENERIC :
    return ((*mdb->close_index_generic)(mdb));
  default :
    warn("Unknown threestatemodel db type.");
    return FALSE;
  }

  warn("Should never get here!");
  return FALSE;

}

# line 583 "threestatedb.c"
/* Function:  hard_link_ThreeStateDB(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [ThreeStateDB *]
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateDB *]
 *
 */
ThreeStateDB * hard_link_ThreeStateDB(ThreeStateDB * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a ThreeStateDB object: passed a NULL object");    
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  ThreeStateDB_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateDB *]
 *
 */
ThreeStateDB * ThreeStateDB_alloc(void) 
{
    ThreeStateDB * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(ThreeStateDB *) ckalloc (sizeof(ThreeStateDB))) == NULL)    {  
      warn("ThreeStateDB_alloc failed ");    
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->dbtype = TSMDB_UNKNOWN; 
    out->filename = NULL;    
    out->type = TSM_default; 
    out->rm = NULL;  
    out->byte_position = 0;  
    out->single = NULL;  
    out->phdb = NULL;    
    out->sdb = NULL; 
    out->comp = NULL;    
    out->gap = 0;    
    out->ext = 0;    
    out->seq_cache = NULL;   
    out->reload_generic = NULL;  
    out->open_generic = NULL;    
    out->close_generic = NULL;   
    out->dataentry_add = NULL;   
    out->open_index_generic = NULL;  
    out->index_generic = NULL;   
    out->close_index_generic = NULL; 
    out->hmm_model_start = -1;   
    out->hmm_model_end = -1; 
    out->current_no = 0; 


    return out;  
}    


/* Function:  free_ThreeStateDB(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [ThreeStateDB *]
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateDB *]
 *
 */
ThreeStateDB * free_ThreeStateDB(ThreeStateDB * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a ThreeStateDB obj. Should be trappable");  
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
    /* obj->current_file is linked in */ 
    if( obj->rm != NULL) 
      free_RandomModel(obj->rm);     
    if( obj->single != NULL) 
      free_ThreeStateModel(obj->single);     
    if( obj->phdb != NULL)   
      free_PfamHmmer1DB(obj->phdb);  
    if( obj->sdb != NULL)    
      free_SequenceDB(obj->sdb);     
    if( obj->comp != NULL)   
      free_CompMat(obj->comp);   
    if( obj->seq_cache != NULL)  
      free_Sequence(obj->seq_cache);     
    /* obj->reload_generic is a function pointer */ 
    /* obj->open_generic is a function pointer */ 
    /* obj->close_generic is a function pointer */ 
    /* obj->dataentry_add is a function pointer */ 
    /* obj->open_index_generic is a function pointer */ 
    /* obj->index_generic is a function pointer */ 
    /* obj->close_index_generic is a function pointer */ 
    /* obj->data is linked in */ 


    ckfree(obj); 
    return NULL; 
}    


/* Function:  replace_dbtype_ThreeStateDB(obj,dbtype)
 *
 * Descrip:    Replace member variable dbtype
 *             For use principly by API functions
 *
 *
 * Arg:           obj [UNKN ] Object holding the variable [ThreeStateDB *]
 * Arg:        dbtype [OWNER] New value of the variable [int]
 *
 * Return [SOFT ]  member variable dbtype [boolean]
 *
 */
boolean replace_dbtype_ThreeStateDB(ThreeStateDB * obj,int dbtype) 
{
    if( obj == NULL)     {  
      warn("In replacement function dbtype for object ThreeStateDB, got a NULL object"); 
      return FALSE;  
      }  
    obj->dbtype = dbtype;    
    return TRUE; 
}    


/* Function:  access_dbtype_ThreeStateDB(obj)
 *
 * Descrip:    Access member variable dbtype
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [ThreeStateDB *]
 *
 * Return [SOFT ]  member variable dbtype [int]
 *
 */
int access_dbtype_ThreeStateDB(ThreeStateDB * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function dbtype for object ThreeStateDB, got a NULL object");    
      return 0;  
      }  
    return obj->dbtype;  
}    


/* Function:  replace_filename_ThreeStateDB(obj,filename)
 *
 * Descrip:    Replace member variable filename
 *             For use principly by API functions
 *
 *
 * Arg:             obj [UNKN ] Object holding the variable [ThreeStateDB *]
 * Arg:        filename [OWNER] New value of the variable [char *]
 *
 * Return [SOFT ]  member variable filename [boolean]
 *
 */
boolean replace_filename_ThreeStateDB(ThreeStateDB * obj,char * filename) 
{
    if( obj == NULL)     {  
      warn("In replacement function filename for object ThreeStateDB, got a NULL object");   
      return FALSE;  
      }  
    obj->filename = filename;    
    return TRUE; 
}    


/* Function:  access_filename_ThreeStateDB(obj)
 *
 * Descrip:    Access member variable filename
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [ThreeStateDB *]
 *
 * Return [SOFT ]  member variable filename [char *]
 *
 */
char * access_filename_ThreeStateDB(ThreeStateDB * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function filename for object ThreeStateDB, got a NULL object");  
      return NULL;   
      }  
    return obj->filename;    
}    



#ifdef _cplusplus
}
#endif
