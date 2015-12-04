#ifdef _cplusplus
extern "C" {
#endif
#include "cdnadb.h"


/* Function:  show_Hscore_cDNADB(hs,ofp)
 *
 * Descrip:    shows the Hscore by the cDNADB information
 *
 *
 *
 * Arg:         hs [UNKN ] High Score structure [Hscore *]
 * Arg:        ofp [UNKN ] output file [FILE *]
 *
 */
# line 75 "cdnadb.dy"
void show_Hscore_cDNADB(Hscore * hs,FILE * ofp)
{
  int i;

  for(i=0;i<hs->len;i++)
    fprintf(ofp,"Query [%20s] Target [%20s] %d\n",hs->ds[i]->query->name,hs->ds[i]->target->name,hs->ds[i]->score);
}

/* Function:  get_cDNA_from_cDNADB(cdnadb,de)
 *
 * Descrip:    Gets cDNA sequence out from
 *             the cDNADB using the information stored in
 *             dataentry
 *
 *
 * Arg:        cdnadb [READ ] cDNA database [cDNADB *]
 * Arg:            de [READ ] DataEntry information  [DataEntry *]
 *
 * Return [UNKN ]  Undocumented return value [cDNA *]
 *
 */
# line 92 "cdnadb.dy"
cDNA * get_cDNA_from_cDNADB(cDNADB * cdnadb,DataEntry * de)
{
  Sequence * seq;
  Sequence * temp;

  if( cdnadb == NULL ) {
    warn("Cannot get entry from a null database");
    return NULL;
  }

  if( de == NULL ) {
    warn("Cannot get entry with a null dataentry");
    return NULL;
  }


  if( cdnadb->is_single_seq == TRUE ) {
    if( de->is_reversed == TRUE ) 
      return cDNA_from_Sequence(hard_link_Sequence(cdnadb->rev->seq));
    else 
      return cDNA_from_Sequence(hard_link_Sequence(cdnadb->forw->seq));
  }

  /* we need to get out the Sequence from seqdb */

  seq = get_Sequence_from_SequenceDB(cdnadb->sdb,de);
  if( seq == NULL ) {
    warn("Cannot get entry for %s from cDNA db",de->name);
    return NULL;
  }

  if( seq->type != SEQUENCE_DNA) {
    warn("Sequence from %s data entry doesn't look like DNA. Forcing it to",de->name);
  }

  force_to_dna_Sequence(seq,1.0,NULL);

  if( de->is_reversed == TRUE ) {
    temp = reverse_complement_Sequence(seq);
    free_Sequence(seq);
    seq = temp;
  }

  return cDNA_from_Sequence(seq);
}



/* Function:  dataentry_add_cDNADB(de,cs,cdnadb)
 *
 * Descrip:    adds information to dataentry from cDNADB
 *
 *             will eventually add file offset and format information,
 *             but this is handled by the SequenceDB mainly.
 *
 *
 * Arg:            de [UNKN ] Undocumented argument [DataEntry *]
 * Arg:            cs [UNKN ] Undocumented argument [ComplexSequence *]
 * Arg:        cdnadb [UNKN ] Undocumented argument [cDNADB *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 146 "cdnadb.dy"
boolean dataentry_add_cDNADB(DataEntry * de,ComplexSequence * cs,cDNADB * cdnadb)
{
  if( cs == NULL || cs->seq == NULL ) {
    warn("Adding a dataentry with a NULL complex sequence or null internal sequence. Nope!");
    return FALSE;
  }

  if( cdnadb->is_single_seq == FALSE) 
    add_SequenceDB_info_DataEntry(cdnadb->sdb,de);
  de->name = stringalloc(cs->seq->name);
  de->is_reversed = is_reversed_Sequence(cs->seq);
  return TRUE;
}


/* Function:  init_cDNADB(cdnadb,return_status)
 *
 * Descrip:    top level function which opens the cDNA database
 *
 *
 * Arg:               cdnadb [UNKN ] protein database [cDNADB *]
 * Arg:        return_status [WRITE] the status of the open from database.h [int *]
 *
 * Return [UNKN ]  Undocumented return value [ComplexSequence *]
 *
 */
# line 167 "cdnadb.dy"
ComplexSequence * init_cDNADB(cDNADB * cdnadb,int * return_status)
{
  ComplexSequence * cs;
  Sequence * seq;

  if( cdnadb->is_single_seq == TRUE) {
    *return_status = DB_RETURN_OK;
    cdnadb->done_forward = TRUE;
    return hard_link_ComplexSequence(cdnadb->forw);
    
  }

  /* is a seq db */

  seq = init_SequenceDB(cdnadb->sdb,return_status);

  if( seq == NULL || *return_status == DB_RETURN_ERROR || *return_status == DB_RETURN_END ) {
    return NULL; /** error already reported **/
  }

  if( force_to_dna_Sequence(seq,cdnadb->error_tol,NULL) == FALSE ) {
    warn("first sequence below error level, have to fail at the moment. Ooops...");
    free_Sequence(seq);
    *return_status = DB_RETURN_ERROR;
    return NULL;
  }

  cdnadb->current = seq;
  cdnadb->done_forward = TRUE;
  cs = new_ComplexSequence(seq,cdnadb->cses);
  if( cs == NULL ) {
    warn("Cannot make initial ComplexSequence. Unable to error catch this. Failing!");
    *return_status = DB_RETURN_ERROR;
    return NULL;
  }



  return cs;
}

/* Function:  reload_cDNADB(last,cdnadb,return_status)
 *
 * Descrip:    function which reloads the database
 *
 *
 * Arg:                 last [UNKN ] previous complex sequence, will be freed [ComplexSequence *]
 * Arg:               cdnadb [UNKN ] Undocumented argument [cDNADB *]
 * Arg:        return_status [WRITE] return_status of the load [int *]
 *
 * Return [OWNER]  a new ComplexSequence object [ComplexSequence *]
 *
 */
# line 215 "cdnadb.dy"
ComplexSequence * reload_cDNADB(ComplexSequence * last,cDNADB * cdnadb,int * return_status)
{
  ComplexSequence * cs;
  Sequence * seq,*temp;
  

  /** free Complex Sequence **/

  if ( last != NULL ) {
    free_ComplexSequence(last);
  }

  if( cdnadb->forward_only == TRUE) {
     temp = reload_SequenceDB(NULL,cdnadb->sdb,return_status);
     if ( *return_status  != DB_RETURN_OK ) {
         return NULL;
     }
    cs = new_ComplexSequence(temp,cdnadb->cses);
    return cs;
  }

  if( cdnadb->is_single_seq == TRUE ) {
    if( cdnadb->done_forward == TRUE ) {
      *return_status = DB_RETURN_OK;
      cdnadb->done_forward = FALSE;
      return hard_link_ComplexSequence(cdnadb->rev);
    } else {
      *return_status = DB_RETURN_END;
      return NULL;
    }
  }

  
  /** standard database **/


  if( cdnadb->done_forward == TRUE ) {
    if( cdnadb->current == NULL ) {
      warn("A bad internal cDNA db error - unable to find current sequence in db reload");
      *return_status = DB_RETURN_ERROR;
      return NULL;
    }

    temp = reverse_complement_Sequence(cdnadb->current);


    if( temp == NULL ) {
      warn("A bad internal cDNA db error - unable to reverse complements current");
      *return_status = DB_RETURN_ERROR;
      return NULL;
    }

    cs = new_ComplexSequence(temp,cdnadb->cses);

    if( cs == NULL ) {
      warn("A bad internal cDNA db error - unable to make complex sequence in db reload");
      *return_status = DB_RETURN_ERROR;
      return NULL;
    }

    free_Sequence(temp);
    cdnadb->done_forward = FALSE;
    return cs;
  }


  /* otherwise we have to get a new sequence */

  seq = reload_SequenceDB(NULL,cdnadb->sdb,return_status);

  if( seq == NULL || *return_status == DB_RETURN_ERROR || *return_status == DB_RETURN_END ) {
    return NULL; /** error already reported **/
  }

  uppercase_Sequence(seq);

  if( force_to_dna_Sequence(seq,cdnadb->error_tol,NULL) == FALSE ) {
    if( cdnadb->error_handling == CDNADB_READ_THROUGH ) {
      warn("Unable to map %s sequence to a cDNA sequence, but ignoring that for the moment...",seq->name);
      free_Sequence(seq);
      return reload_cDNADB(NULL,cdnadb,return_status);
    } else {
      warn("Unable to map %s sequence to a cDNA sequence. Failing",seq->name);
      *return_status = DB_RETURN_ERROR;
      return NULL;
    }
  }


  cs = new_ComplexSequence(seq,cdnadb->cses);
  if( cs == NULL ) {
    if( cdnadb->error_handling == CDNADB_READ_THROUGH ) {
      warn("Unable to map %s sequence to a cDNA sequence, but ignoring that for the moment...",seq->name);
      free_Sequence(seq);
      return reload_cDNADB(NULL,cdnadb,return_status);
    } else {
      warn("Unable to map %s sequence to a cDNA sequence. Failing",seq->name);
      *return_status = DB_RETURN_ERROR;
      return NULL;
    }
  }

  cdnadb->current = free_Sequence(cdnadb->current);
  cdnadb->current = seq;
  cdnadb->done_forward= TRUE;

  return cs;
}
  


/* Function:  close_cDNADB(cs,cdnadb)
 *
 * Descrip:    top level function which closes the cDNA database
 *
 *
 * Arg:            cs [UNKN ] last complex sequence  [ComplexSequence *]
 * Arg:        cdnadb [UNKN ] protein database [cDNADB *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 332 "cdnadb.dy"
boolean close_cDNADB(ComplexSequence * cs,cDNADB * cdnadb) 
{
  if( cdnadb->is_single_seq == TRUE ) {
    return TRUE;
  }

  if( cs != NULL)
    free_ComplexSequence(cs);

  if( cdnadb->current != NULL)
    cdnadb->current = free_Sequence(cdnadb->current);

  return close_SequenceDB(NULL,cdnadb->sdb);
}

/* Function:  new_cDNADB_from_single_seq(seq)
 *
 * Descrip:    To make a new cDNA database
 *             from a single cDNA Sequence with a eval system
 *
 *
 * Arg:        seq [UNKN ] sequence which as placed into cDNADB structure. [cDNA *]
 *
 * Return [UNKN ]  Undocumented return value [cDNADB *]
 *
 */
# line 353 "cdnadb.dy"
cDNADB * new_cDNADB_from_single_seq(cDNA * seq)
{
  ComplexSequence * cs,*cs_rev;
  Sequence * temp;
  ComplexSequenceEvalSet * cses;
  
  cses = default_cDNA_ComplexSequenceEvalSet();

  cs = new_ComplexSequence(seq->baseseq,cses);
  temp = reverse_complement_Sequence(seq->baseseq);
  cs_rev = new_ComplexSequence(temp,cses);
  free_Sequence(temp);
  free_ComplexSequenceEvalSet(cses);

  return new_cDNADB_from_forrev_cseq(cs,cs_rev);
}
  

/* Function:  new_cDNADB_from_forrev_cseq(cs,cs_rev)
 *
 * Descrip:    To make a new cDNA database
 *             from a single ComplexSequence
 *
 *
 * Arg:            cs [OWNER] complex sequence which is held. [ComplexSequence *]
 * Arg:        cs_rev [OWNER] complex sequence which is held. [ComplexSequence *]
 *
 * Return [UNKN ]  Undocumented return value [cDNADB *]
 *
 */
# line 379 "cdnadb.dy"
cDNADB * new_cDNADB_from_forrev_cseq(ComplexSequence * cs,ComplexSequence * cs_rev)
{
  cDNADB * out;


  out = cDNADB_alloc();
  out->is_single_seq = TRUE;
  out->forw = cs;
  out->rev = cs_rev;
  
  return out;
}


  
/* Function:  new_cDNADB(seqdb)
 *
 * Descrip:    To make a new cDNA database
 *
 *
 * Arg:        seqdb [READ ] sequence database [SequenceDB *]
 *
 * Return [UNKN ]  Undocumented return value [cDNADB *]
 *
 */
# line 399 "cdnadb.dy"
cDNADB * new_cDNADB(SequenceDB * seqdb)
{
  cDNADB * out;
  ComplexSequenceEvalSet * cses;
  

  if( seqdb == NULL ) {
    warn("No sequencedb - can't make a cDNADB!");
    return NULL;
  }

  /** should check sequence database **/
  cses = default_cDNA_ComplexSequenceEvalSet();

  if( cses->type != SEQUENCE_CDNA ) {
    warn("You can't make a cDNA database with a non SEQUENCE_cDNA cses type [%d]",cses->type);
    return NULL;
  }


  out = cDNADB_alloc();

  out->is_single_seq = FALSE;
  out->sdb  = hard_link_SequenceDB(seqdb);
  out->cses = hard_link_ComplexSequenceEvalSet(cses);
  free_ComplexSequenceEvalSet(cses);

  return out;
}

 
# line 404 "cdnadb.c"
/* Function:  hard_link_cDNADB(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [cDNADB *]
 *
 * Return [UNKN ]  Undocumented return value [cDNADB *]
 *
 */
cDNADB * hard_link_cDNADB(cDNADB * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a cDNADB object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  cDNADB_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [cDNADB *]
 *
 */
cDNADB * cDNADB_alloc(void) 
{
    cDNADB * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(cDNADB *) ckalloc (sizeof(cDNADB))) == NULL)    {  
      warn("cDNADB_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->is_single_seq = FALSE;  
    out->done_forward = FALSE;   
    out->forward_only = FALSE;   
    out->forw = NULL;    
    out->rev = NULL; 
    out->sdb = NULL; 
    out->current = NULL; 
    out->cses = NULL;    
    out->error_handling = CDNADB_READ_THROUGH;   
    out->error_tol = 0.7;    


    return out;  
}    


/* Function:  free_cDNADB(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [cDNADB *]
 *
 * Return [UNKN ]  Undocumented return value [cDNADB *]
 *
 */
cDNADB * free_cDNADB(cDNADB * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a cDNADB obj. Should be trappable");    
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
    if( obj->forw != NULL)   
      free_ComplexSequence(obj->forw);   
    if( obj->rev != NULL)    
      free_ComplexSequence(obj->rev);    
    if( obj->sdb != NULL)    
      free_SequenceDB(obj->sdb);     
    if( obj->current != NULL)    
      free_Sequence(obj->current);   
    if( obj->cses != NULL)   
      free_ComplexSequenceEvalSet(obj->cses);    


    ckfree(obj); 
    return NULL; 
}    


/* Function:  replace_is_single_seq_cDNADB(obj,is_single_seq)
 *
 * Descrip:    Replace member variable is_single_seq
 *             For use principly by API functions
 *
 *
 * Arg:                  obj [UNKN ] Object holding the variable [cDNADB *]
 * Arg:        is_single_seq [OWNER] New value of the variable [boolean]
 *
 * Return [SOFT ]  member variable is_single_seq [boolean]
 *
 */
boolean replace_is_single_seq_cDNADB(cDNADB * obj,boolean is_single_seq) 
{
    if( obj == NULL)     {  
      warn("In replacement function is_single_seq for object cDNADB, got a NULL object");    
      return FALSE;  
      }  
    obj->is_single_seq = is_single_seq;  
    return TRUE; 
}    


/* Function:  access_is_single_seq_cDNADB(obj)
 *
 * Descrip:    Access member variable is_single_seq
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [cDNADB *]
 *
 * Return [SOFT ]  member variable is_single_seq [boolean]
 *
 */
boolean access_is_single_seq_cDNADB(cDNADB * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function is_single_seq for object cDNADB, got a NULL object");   
      return FALSE;  
      }  
    return obj->is_single_seq;   
}    


/* Function:  replace_done_forward_cDNADB(obj,done_forward)
 *
 * Descrip:    Replace member variable done_forward
 *             For use principly by API functions
 *
 *
 * Arg:                 obj [UNKN ] Object holding the variable [cDNADB *]
 * Arg:        done_forward [OWNER] New value of the variable [boolean]
 *
 * Return [SOFT ]  member variable done_forward [boolean]
 *
 */
boolean replace_done_forward_cDNADB(cDNADB * obj,boolean done_forward) 
{
    if( obj == NULL)     {  
      warn("In replacement function done_forward for object cDNADB, got a NULL object"); 
      return FALSE;  
      }  
    obj->done_forward = done_forward;    
    return TRUE; 
}    


/* Function:  access_done_forward_cDNADB(obj)
 *
 * Descrip:    Access member variable done_forward
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [cDNADB *]
 *
 * Return [SOFT ]  member variable done_forward [boolean]
 *
 */
boolean access_done_forward_cDNADB(cDNADB * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function done_forward for object cDNADB, got a NULL object");    
      return FALSE;  
      }  
    return obj->done_forward;    
}    


/* Function:  replace_forward_only_cDNADB(obj,forward_only)
 *
 * Descrip:    Replace member variable forward_only
 *             For use principly by API functions
 *
 *
 * Arg:                 obj [UNKN ] Object holding the variable [cDNADB *]
 * Arg:        forward_only [OWNER] New value of the variable [boolean]
 *
 * Return [SOFT ]  member variable forward_only [boolean]
 *
 */
boolean replace_forward_only_cDNADB(cDNADB * obj,boolean forward_only) 
{
    if( obj == NULL)     {  
      warn("In replacement function forward_only for object cDNADB, got a NULL object"); 
      return FALSE;  
      }  
    obj->forward_only = forward_only;    
    return TRUE; 
}    


/* Function:  access_forward_only_cDNADB(obj)
 *
 * Descrip:    Access member variable forward_only
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [cDNADB *]
 *
 * Return [SOFT ]  member variable forward_only [boolean]
 *
 */
boolean access_forward_only_cDNADB(cDNADB * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function forward_only for object cDNADB, got a NULL object");    
      return FALSE;  
      }  
    return obj->forward_only;    
}    


/* Function:  replace_forw_cDNADB(obj,forw)
 *
 * Descrip:    Replace member variable forw
 *             For use principly by API functions
 *
 *
 * Arg:         obj [UNKN ] Object holding the variable [cDNADB *]
 * Arg:        forw [OWNER] New value of the variable [ComplexSequence *]
 *
 * Return [SOFT ]  member variable forw [boolean]
 *
 */
boolean replace_forw_cDNADB(cDNADB * obj,ComplexSequence * forw) 
{
    if( obj == NULL)     {  
      warn("In replacement function forw for object cDNADB, got a NULL object"); 
      return FALSE;  
      }  
    obj->forw = forw;    
    return TRUE; 
}    


/* Function:  access_forw_cDNADB(obj)
 *
 * Descrip:    Access member variable forw
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [cDNADB *]
 *
 * Return [SOFT ]  member variable forw [ComplexSequence *]
 *
 */
ComplexSequence * access_forw_cDNADB(cDNADB * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function forw for object cDNADB, got a NULL object");    
      return NULL;   
      }  
    return obj->forw;    
}    


/* Function:  replace_rev_cDNADB(obj,rev)
 *
 * Descrip:    Replace member variable rev
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [cDNADB *]
 * Arg:        rev [OWNER] New value of the variable [ComplexSequence *]
 *
 * Return [SOFT ]  member variable rev [boolean]
 *
 */
boolean replace_rev_cDNADB(cDNADB * obj,ComplexSequence * rev) 
{
    if( obj == NULL)     {  
      warn("In replacement function rev for object cDNADB, got a NULL object");  
      return FALSE;  
      }  
    obj->rev = rev;  
    return TRUE; 
}    


/* Function:  access_rev_cDNADB(obj)
 *
 * Descrip:    Access member variable rev
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [cDNADB *]
 *
 * Return [SOFT ]  member variable rev [ComplexSequence *]
 *
 */
ComplexSequence * access_rev_cDNADB(cDNADB * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function rev for object cDNADB, got a NULL object"); 
      return NULL;   
      }  
    return obj->rev;     
}    


/* Function:  replace_sdb_cDNADB(obj,sdb)
 *
 * Descrip:    Replace member variable sdb
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [cDNADB *]
 * Arg:        sdb [OWNER] New value of the variable [SequenceDB *]
 *
 * Return [SOFT ]  member variable sdb [boolean]
 *
 */
boolean replace_sdb_cDNADB(cDNADB * obj,SequenceDB * sdb) 
{
    if( obj == NULL)     {  
      warn("In replacement function sdb for object cDNADB, got a NULL object");  
      return FALSE;  
      }  
    obj->sdb = sdb;  
    return TRUE; 
}    


/* Function:  access_sdb_cDNADB(obj)
 *
 * Descrip:    Access member variable sdb
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [cDNADB *]
 *
 * Return [SOFT ]  member variable sdb [SequenceDB *]
 *
 */
SequenceDB * access_sdb_cDNADB(cDNADB * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function sdb for object cDNADB, got a NULL object"); 
      return NULL;   
      }  
    return obj->sdb;     
}    


/* Function:  replace_current_cDNADB(obj,current)
 *
 * Descrip:    Replace member variable current
 *             For use principly by API functions
 *
 *
 * Arg:            obj [UNKN ] Object holding the variable [cDNADB *]
 * Arg:        current [OWNER] New value of the variable [Sequence *]
 *
 * Return [SOFT ]  member variable current [boolean]
 *
 */
boolean replace_current_cDNADB(cDNADB * obj,Sequence * current) 
{
    if( obj == NULL)     {  
      warn("In replacement function current for object cDNADB, got a NULL object");  
      return FALSE;  
      }  
    obj->current = current;  
    return TRUE; 
}    


/* Function:  access_current_cDNADB(obj)
 *
 * Descrip:    Access member variable current
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [cDNADB *]
 *
 * Return [SOFT ]  member variable current [Sequence *]
 *
 */
Sequence * access_current_cDNADB(cDNADB * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function current for object cDNADB, got a NULL object"); 
      return NULL;   
      }  
    return obj->current;     
}    


/* Function:  replace_cses_cDNADB(obj,cses)
 *
 * Descrip:    Replace member variable cses
 *             For use principly by API functions
 *
 *
 * Arg:         obj [UNKN ] Object holding the variable [cDNADB *]
 * Arg:        cses [OWNER] New value of the variable [ComplexSequenceEvalSet *]
 *
 * Return [SOFT ]  member variable cses [boolean]
 *
 */
boolean replace_cses_cDNADB(cDNADB * obj,ComplexSequenceEvalSet * cses) 
{
    if( obj == NULL)     {  
      warn("In replacement function cses for object cDNADB, got a NULL object"); 
      return FALSE;  
      }  
    obj->cses = cses;    
    return TRUE; 
}    


/* Function:  access_cses_cDNADB(obj)
 *
 * Descrip:    Access member variable cses
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [cDNADB *]
 *
 * Return [SOFT ]  member variable cses [ComplexSequenceEvalSet *]
 *
 */
ComplexSequenceEvalSet * access_cses_cDNADB(cDNADB * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function cses for object cDNADB, got a NULL object");    
      return NULL;   
      }  
    return obj->cses;    
}    



#ifdef _cplusplus
}
#endif
