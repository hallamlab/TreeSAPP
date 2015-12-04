#ifdef _cplusplus
extern "C" {
#endif
#include "proteindb.h"


/* Function:  show_Hscore_ProteinDB(hs,ofp)
 *
 * Descrip:    shows the Hscore by the ProteinDB information
 *
 *
 *
 * Arg:         hs [UNKN ] High Score structure [Hscore *]
 * Arg:        ofp [UNKN ] output file [FILE *]
 *
 */
# line 61 "proteindb.dy"
void show_Hscore_ProteinDB(Hscore * hs,FILE * ofp)
{
  int i;

  for(i=0;i<hs->len;i++)
    fprintf(ofp,"Query [%20s] Target [%20s] %d\n",hs->ds[i]->query->name,hs->ds[i]->target->name,hs->ds[i]->score);
}

/* Function:  get_Protein_from_ProteinDB(prodb,de)
 *
 * Descrip:    Gets Protein sequence out from
 *             the proteindb using the information stored in
 *             dataentry
 *
 *
 * Arg:        prodb [READ ] ProteinDB database [ProteinDB *]
 * Arg:           de [READ ] DataEntry information  [DataEntry *]
 *
 * Return [UNKN ]  Undocumented return value [Protein *]
 *
 */
# line 78 "proteindb.dy"
Protein * get_Protein_from_ProteinDB(ProteinDB * prodb,DataEntry * de)
{
  Sequence * seq;

  if( prodb->is_single_seq == TRUE ) {
    return Protein_from_Sequence(hard_link_Sequence(prodb->single->seq));
  }

  /* we need to get out the Sequence from seqdb */

  seq = get_Sequence_from_SequenceDB(prodb->sdb,de);
  if( seq == NULL ) {
    warn("Cannot get entry for %s from Protein db",de->name);
    return NULL;
  }

  seq->type = SEQUENCE_PROTEIN; /* force to protein */

  return Protein_from_Sequence(seq);
}


/* Function:  dataentry_add_ProteinDB(de,cs,prodb)
 *
 * Descrip:    adds information to dataentry from ProteinDB
 *
 *             This information is the necessary information for
 *             the proteindb to find this sequence later
 *
 *
 * Arg:           de [UNKN ] Undocumented argument [DataEntry *]
 * Arg:           cs [UNKN ] Undocumented argument [ComplexSequence *]
 * Arg:        prodb [UNKN ] Undocumented argument [ProteinDB *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 106 "proteindb.dy"
boolean dataentry_add_ProteinDB(DataEntry * de,ComplexSequence * cs,ProteinDB * prodb)
{
  if( cs == NULL || cs->seq == NULL ){
     warn("Adding a dataentry with a NULL complex sequence or null internal sequence. Nope!");
     return FALSE;
  }

  if( prodb->is_single_seq == FALSE)
     add_SequenceDB_info_DataEntry(prodb->sdb,de);

  de->name = stringalloc(cs->seq->name);

  return TRUE;
}


/* Function:  init_ProteinDB(prodb,return_status)
 *
 * Descrip:    top level function which opens the protein database
 *
 *
 * Arg:                prodb [UNKN ] protein database [ProteinDB *]
 * Arg:        return_status [WRITE] the status of the open from database.h [int *]
 *
 * Return [UNKN ]  Undocumented return value [ComplexSequence *]
 *
 */
# line 128 "proteindb.dy"
ComplexSequence * init_ProteinDB(ProteinDB * prodb,int * return_status)
{
  ComplexSequence * cs;
  Sequence * seq;

  if( prodb->is_single_seq == TRUE) {
    *return_status = DB_RETURN_OK;
    return prodb->single;
  }

  seq = init_SequenceDB(prodb->sdb,return_status);

  if( seq == NULL || *return_status == DB_RETURN_ERROR || *return_status == DB_RETURN_END ) {
    return NULL; /** error already reported **/
  }

  if( prodb->test_dna == FALSE ) {
    seq->type = SEQUENCE_PROTEIN;
  } else {
    if( seq->type != SEQUENCE_PROTEIN ) {
      warn("For sequence %s, looks like a DNA sequence. Failing");
      *return_status = DB_RETURN_ERROR;
    }
  }

  cs = new_ComplexSequence(seq,prodb->cses);

  free_Sequence(seq);

  return cs;
}

/* Function:  reload_ProteinDB(last,prodb,return_status)
 *
 * Descrip:    function which reloads the database
 *
 *
 * Arg:                 last [UNKN ] previous complex sequence, will be freed [ComplexSequence *]
 * Arg:                prodb [UNKN ] Undocumented argument [ProteinDB *]
 * Arg:        return_status [WRITE] return_status of the load [int *]
 *
 * Return [UNKN ]  Undocumented return value [ComplexSequence *]
 *
 */
# line 166 "proteindb.dy"
ComplexSequence * reload_ProteinDB(ComplexSequence * last,ProteinDB * prodb,int * return_status)
{
  ComplexSequence * cs;
  Sequence * seq;
    
  if( prodb->is_single_seq == TRUE ) {
    *return_status = DB_RETURN_END;
    return NULL;
  }

  /** free Complex Sequence **/
  if( last != NULL ) {
    free_ComplexSequence(last);
  }

  seq = reload_SequenceDB(NULL,prodb->sdb,return_status);

  if( seq == NULL || *return_status == DB_RETURN_ERROR || *return_status == DB_RETURN_END ) {
    return NULL; /** error already reported **/
  }
  if( prodb->test_dna == FALSE ) {
    seq->type = SEQUENCE_PROTEIN;
  } else {
    if( seq->type != SEQUENCE_PROTEIN ) {
      warn("For sequence %s, looks like a DNA sequence. Failing");
      *return_status = DB_RETURN_ERROR;
    }
  }

  cs = new_ComplexSequence(seq,prodb->cses);

  free_Sequence(seq);

  return cs;
}
  


/* Function:  close_ProteinDB(cs,prodb)
 *
 * Descrip:    top level function which closes the protein database
 *
 *
 * Arg:           cs [UNKN ] last complex sequence  [ComplexSequence *]
 * Arg:        prodb [UNKN ] protein database [ProteinDB *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 210 "proteindb.dy"
boolean close_ProteinDB(ComplexSequence * cs,ProteinDB * prodb) 
{
  if( prodb->is_single_seq == TRUE ) {
    return TRUE;
  }

  if( cs != NULL) {
    free_ComplexSequence(cs);
  }

  return close_SequenceDB(NULL,prodb->sdb);
}

/* Function:  new_ProteinDB_from_single_seq(seq)
 *
 * Descrip:    To make a new protein database
 *             from a single Sequence with default amino acid mapping
 *
 *
 * Arg:        seq [UNKN ] sequence which as placed into ProteinDB structure. [Sequence *]
 *
 * Return [UNKN ]  Undocumented return value [ProteinDB *]
 *
 */
# line 229 "proteindb.dy"
ProteinDB * new_ProteinDB_from_single_seq(Sequence * seq)
{
  ComplexSequenceEvalSet * cses;
  ComplexSequence * cs;

  cses = default_aminoacid_ComplexSequenceEvalSet();

  cs = new_ComplexSequence(seq,cses);

  free_ComplexSequenceEvalSet(cses);

  return new_ProteinDB_from_single_cseq(cs);
}
  

/* Function:  new_ProteinDB_from_single_cseq(cs)
 *
 * Descrip:    To make a new protein database
 *             from a single ComplexSequence
 *
 *
 * Arg:        cs [UNKN ] complex sequence which is held. [ComplexSequence *]
 *
 * Return [UNKN ]  Undocumented return value [ProteinDB *]
 *
 */
# line 250 "proteindb.dy"
ProteinDB * new_ProteinDB_from_single_cseq(ComplexSequence * cs)
{
  ProteinDB * out;


  out = ProteinDB_alloc();
  out->is_single_seq = TRUE;
  out->single = cs;
  
  return out;
}

/* Function:  single_fasta_ProteinDB(filename)
 *
 * Descrip:    pre-packed single fasta protein database
 *
 *
 *
 * Arg:        filename [UNKN ] name of fasta file [char *]
 *
 * Return [UNKN ]  Undocumented return value [ProteinDB *]
 *
 */
# line 268 "proteindb.dy"
ProteinDB * single_fasta_ProteinDB(char * filename)
{
  return new_ProteinDB(single_fasta_SequenceDB(filename),default_aminoacid_ComplexSequenceEvalSet());
}

  
/* Function:  new_ProteinDB(seqdb,cses)
 *
 * Descrip:    To make a new protein database
 *
 *
 * Arg:        seqdb [UNKN ] sequence database [SequenceDB *]
 * Arg:         cses [UNKN ] protein evaluation set [ComplexSequenceEvalSet *]
 *
 * Return [UNKN ]  Undocumented return value [ProteinDB *]
 *
 */
# line 280 "proteindb.dy"
ProteinDB * new_ProteinDB(SequenceDB * seqdb,ComplexSequenceEvalSet * cses)
{
  ProteinDB * out;

  /** should check sequence database **/

  if( seqdb == NULL ) {
     warn("Cannot make ProteinDB from NULL SequenceDB object");
     return NULL;
  }

  if( cses->type != SEQUENCE_PROTEIN ) {
    warn("You can't make a protein database with a non SEQUENCE_PROTEIN cses type [%d]",cses->type);
    return NULL;
  }


  out = ProteinDB_alloc();

  out->sdb  = seqdb;
  out->cses = cses;

  return out;
}

 
# line 300 "proteindb.c"
/* Function:  hard_link_ProteinDB(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [ProteinDB *]
 *
 * Return [UNKN ]  Undocumented return value [ProteinDB *]
 *
 */
ProteinDB * hard_link_ProteinDB(ProteinDB * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a ProteinDB object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  ProteinDB_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ProteinDB *]
 *
 */
ProteinDB * ProteinDB_alloc(void) 
{
    ProteinDB * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(ProteinDB *) ckalloc (sizeof(ProteinDB))) == NULL)  {  
      warn("ProteinDB_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->is_single_seq = FALSE;  
    out->is_random_db = FALSE;   
    out->single = NULL;  
    out->sdb = NULL; 
    out->cses = NULL;    
    out->rnd = NULL; 
    out->test_dna = FALSE;   


    return out;  
}    


/* Function:  free_ProteinDB(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [ProteinDB *]
 *
 * Return [UNKN ]  Undocumented return value [ProteinDB *]
 *
 */
ProteinDB * free_ProteinDB(ProteinDB * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a ProteinDB obj. Should be trappable"); 
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
    if( obj->single != NULL) 
      free_ComplexSequence(obj->single);     
    if( obj->sdb != NULL)    
      free_SequenceDB(obj->sdb);     
    if( obj->cses != NULL)   
      free_ComplexSequenceEvalSet(obj->cses);    
    if( obj->rnd != NULL)    
      free_RandomProteinDB(obj->rnd);    


    ckfree(obj); 
    return NULL; 
}    


/* Function:  replace_is_single_seq_ProteinDB(obj,is_single_seq)
 *
 * Descrip:    Replace member variable is_single_seq
 *             For use principly by API functions
 *
 *
 * Arg:                  obj [UNKN ] Object holding the variable [ProteinDB *]
 * Arg:        is_single_seq [OWNER] New value of the variable [boolean]
 *
 * Return [SOFT ]  member variable is_single_seq [boolean]
 *
 */
boolean replace_is_single_seq_ProteinDB(ProteinDB * obj,boolean is_single_seq) 
{
    if( obj == NULL)     {  
      warn("In replacement function is_single_seq for object ProteinDB, got a NULL object"); 
      return FALSE;  
      }  
    obj->is_single_seq = is_single_seq;  
    return TRUE; 
}    


/* Function:  access_is_single_seq_ProteinDB(obj)
 *
 * Descrip:    Access member variable is_single_seq
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [ProteinDB *]
 *
 * Return [SOFT ]  member variable is_single_seq [boolean]
 *
 */
boolean access_is_single_seq_ProteinDB(ProteinDB * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function is_single_seq for object ProteinDB, got a NULL object");    
      return FALSE;  
      }  
    return obj->is_single_seq;   
}    


/* Function:  replace_is_random_db_ProteinDB(obj,is_random_db)
 *
 * Descrip:    Replace member variable is_random_db
 *             For use principly by API functions
 *
 *
 * Arg:                 obj [UNKN ] Object holding the variable [ProteinDB *]
 * Arg:        is_random_db [OWNER] New value of the variable [boolean]
 *
 * Return [SOFT ]  member variable is_random_db [boolean]
 *
 */
boolean replace_is_random_db_ProteinDB(ProteinDB * obj,boolean is_random_db) 
{
    if( obj == NULL)     {  
      warn("In replacement function is_random_db for object ProteinDB, got a NULL object");  
      return FALSE;  
      }  
    obj->is_random_db = is_random_db;    
    return TRUE; 
}    


/* Function:  access_is_random_db_ProteinDB(obj)
 *
 * Descrip:    Access member variable is_random_db
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [ProteinDB *]
 *
 * Return [SOFT ]  member variable is_random_db [boolean]
 *
 */
boolean access_is_random_db_ProteinDB(ProteinDB * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function is_random_db for object ProteinDB, got a NULL object"); 
      return FALSE;  
      }  
    return obj->is_random_db;    
}    


/* Function:  replace_single_ProteinDB(obj,single)
 *
 * Descrip:    Replace member variable single
 *             For use principly by API functions
 *
 *
 * Arg:           obj [UNKN ] Object holding the variable [ProteinDB *]
 * Arg:        single [OWNER] New value of the variable [ComplexSequence *]
 *
 * Return [SOFT ]  member variable single [boolean]
 *
 */
boolean replace_single_ProteinDB(ProteinDB * obj,ComplexSequence * single) 
{
    if( obj == NULL)     {  
      warn("In replacement function single for object ProteinDB, got a NULL object");    
      return FALSE;  
      }  
    obj->single = single;    
    return TRUE; 
}    


/* Function:  access_single_ProteinDB(obj)
 *
 * Descrip:    Access member variable single
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [ProteinDB *]
 *
 * Return [SOFT ]  member variable single [ComplexSequence *]
 *
 */
ComplexSequence * access_single_ProteinDB(ProteinDB * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function single for object ProteinDB, got a NULL object");   
      return NULL;   
      }  
    return obj->single;  
}    


/* Function:  replace_sdb_ProteinDB(obj,sdb)
 *
 * Descrip:    Replace member variable sdb
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [ProteinDB *]
 * Arg:        sdb [OWNER] New value of the variable [SequenceDB *]
 *
 * Return [SOFT ]  member variable sdb [boolean]
 *
 */
boolean replace_sdb_ProteinDB(ProteinDB * obj,SequenceDB * sdb) 
{
    if( obj == NULL)     {  
      warn("In replacement function sdb for object ProteinDB, got a NULL object");   
      return FALSE;  
      }  
    obj->sdb = sdb;  
    return TRUE; 
}    


/* Function:  access_sdb_ProteinDB(obj)
 *
 * Descrip:    Access member variable sdb
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [ProteinDB *]
 *
 * Return [SOFT ]  member variable sdb [SequenceDB *]
 *
 */
SequenceDB * access_sdb_ProteinDB(ProteinDB * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function sdb for object ProteinDB, got a NULL object");  
      return NULL;   
      }  
    return obj->sdb;     
}    


/* Function:  replace_cses_ProteinDB(obj,cses)
 *
 * Descrip:    Replace member variable cses
 *             For use principly by API functions
 *
 *
 * Arg:         obj [UNKN ] Object holding the variable [ProteinDB *]
 * Arg:        cses [OWNER] New value of the variable [ComplexSequenceEvalSet *]
 *
 * Return [SOFT ]  member variable cses [boolean]
 *
 */
boolean replace_cses_ProteinDB(ProteinDB * obj,ComplexSequenceEvalSet * cses) 
{
    if( obj == NULL)     {  
      warn("In replacement function cses for object ProteinDB, got a NULL object");  
      return FALSE;  
      }  
    obj->cses = cses;    
    return TRUE; 
}    


/* Function:  access_cses_ProteinDB(obj)
 *
 * Descrip:    Access member variable cses
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [ProteinDB *]
 *
 * Return [SOFT ]  member variable cses [ComplexSequenceEvalSet *]
 *
 */
ComplexSequenceEvalSet * access_cses_ProteinDB(ProteinDB * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function cses for object ProteinDB, got a NULL object"); 
      return NULL;   
      }  
    return obj->cses;    
}    


/* Function:  replace_rnd_ProteinDB(obj,rnd)
 *
 * Descrip:    Replace member variable rnd
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [ProteinDB *]
 * Arg:        rnd [OWNER] New value of the variable [RandomProteinDB *]
 *
 * Return [SOFT ]  member variable rnd [boolean]
 *
 */
boolean replace_rnd_ProteinDB(ProteinDB * obj,RandomProteinDB * rnd) 
{
    if( obj == NULL)     {  
      warn("In replacement function rnd for object ProteinDB, got a NULL object");   
      return FALSE;  
      }  
    obj->rnd = rnd;  
    return TRUE; 
}    


/* Function:  access_rnd_ProteinDB(obj)
 *
 * Descrip:    Access member variable rnd
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [ProteinDB *]
 *
 * Return [SOFT ]  member variable rnd [RandomProteinDB *]
 *
 */
RandomProteinDB * access_rnd_ProteinDB(ProteinDB * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function rnd for object ProteinDB, got a NULL object");  
      return NULL;   
      }  
    return obj->rnd;     
}    


/* Function:  replace_test_dna_ProteinDB(obj,test_dna)
 *
 * Descrip:    Replace member variable test_dna
 *             For use principly by API functions
 *
 *
 * Arg:             obj [UNKN ] Object holding the variable [ProteinDB *]
 * Arg:        test_dna [OWNER] New value of the variable [boolean]
 *
 * Return [SOFT ]  member variable test_dna [boolean]
 *
 */
boolean replace_test_dna_ProteinDB(ProteinDB * obj,boolean test_dna) 
{
    if( obj == NULL)     {  
      warn("In replacement function test_dna for object ProteinDB, got a NULL object");  
      return FALSE;  
      }  
    obj->test_dna = test_dna;    
    return TRUE; 
}    


/* Function:  access_test_dna_ProteinDB(obj)
 *
 * Descrip:    Access member variable test_dna
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [ProteinDB *]
 *
 * Return [SOFT ]  member variable test_dna [boolean]
 *
 */
boolean access_test_dna_ProteinDB(ProteinDB * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function test_dna for object ProteinDB, got a NULL object"); 
      return FALSE;  
      }  
    return obj->test_dna;    
}    



#ifdef _cplusplus
}
#endif
