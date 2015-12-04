#ifdef _cplusplus
extern "C" {
#endif
#include "wise2xhmmer2.h"

#ifndef HMMER2_EXTERN_DEFINED
#include "globals.h"
#endif

/* Function:  bootstrap_HMMer2(void)
 *
 * Descrip:    You need to call this function before you
 *             call any other to get the HMMer2 system working ;)
 *
 *
 *
 */
# line 45 "wise2xhmmer2.dy"
void bootstrap_HMMer2(void)
{
  SetAlphabet(hmmAMINO);
}

/* Function:  read_SeqAlign_HMMER(seqfile)
 *
 * Descrip:    This function read a single SeqAlign from
 *             the filename using the HMMER libraries
 *
 *
 * Arg:        seqfile [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [SeqAlign *]
 *
 */
# line 54 "wise2xhmmer2.dy"
SeqAlign * read_SeqAlign_HMMER(char * seqfile)
{
  AINFO ainfo;
  char ** data;
  int format;
  SeqAlign * out;
  int i;
  Sequence * temp;

  if (! SeqfileFormat(seqfile, &format, NULL) ) {
    switch (squid_errno) {
    case SQERR_NOFILE: 
      warn("Alignment file %s could not be opened for reading", seqfile);
      return NULL;
    case SQERR_FORMAT: 
    default:           
      warn("Failed to determine format of alignment file %s", seqfile);
    }
  }
  
  /* attempt to read it in now */

  if (! ReadAlignment(seqfile, format, &data, &ainfo)) {
    warn("Failed to read aligned sequence file %s", seqfile);
    return NULL;
  }

  /* make SeqAlign */

  out = SeqAlign_alloc_len(ainfo.nseq);
  
  for (i = 0; i < ainfo.nseq; i++) {
    temp = new_Sequence_from_strings(ainfo.sqinfo[i].name,data[i]);
    add_SeqAlign(out,temp);
  }

  return out;
}
    

  


/* Function:  HMMer2_read_ThreeStateModel(filename)
 *
 * Descrip:    This function reads a single HMM from
 *             filename for a ThreeStateModel
 *
 *
 * Arg:        filename [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateModel *]
 *
 */
# line 101 "wise2xhmmer2.dy"
ThreeStateModel * HMMer2_read_ThreeStateModel(char * filename)
{
  struct plan7_s *hmm;
  ThreeStateModel * out;
  HMMFILE * hmmfp;

  if( filename == NULL ) {
    warn("You have passed in a NULL filename into HMMer2_read_ThreeStateModel - unlikely to work!");
    return NULL;
  }

  hmmfp = HMMFileOpen(filename,"HMMERDB");

  if( hmmfp == NULL ) {
    warn("Could not open %s for HMM reading",filename);
    return NULL;
  }

  if( HMMFileRead(hmmfp, &hmm) == 0 ) {
    return NULL;
  }

  if( hmm == NULL ) {
    warn("Unable to read HMM file... ");
    return NULL;
  }
  
  out = ThreeStateModel_from_HMMer2(hmm);
  FreePlan7(hmm);

  return out;
}

/* Function:  HMMer2_ThreeStateDB(filename)
 *
 * Descrip:    This function makes a new ThreeStateDB
 *             from a filename for the HMMer2 system.
 *
 *
 * Arg:        filename [SOFT ] filename of the HMMer db [char *]
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateDB *]
 *
 */
# line 140 "wise2xhmmer2.dy"
ThreeStateDB * HMMer2_ThreeStateDB(char * filename)
{
  ThreeStateDB * out;

  out = ThreeStateDB_alloc();
  out->dbtype    = TSMDB_GENERIC;
  out->filename = stringalloc(filename);
  out->reload_generic = reload_ThreeStateDB_HMMer2;
  out->open_generic   = open_ThreeStateDB_HMMer2;
  out->close_generic  = close_ThreeStateDB_HMMer2;
  out->dataentry_add  = dataentry_ThreeStateDB_HMMer2;
  out->open_index_generic = open_for_indexing_ThreeStateDB_HMMer2;
  out->close_index_generic = close_for_indexing_ThreeStateDB_HMMer2;
  out->index_generic = index_ThreeStateDB_HMMer2;
    
  return out;
}

/* Function:  open_for_indexing_ThreeStateDB_HMMer2(tdb)
 *
 * Descrip:    This function is the generic wrapper for the
 *             threestatedb open for indexing 
 *
 *
 * Arg:        tdb [UNKN ] Undocumented argument [ThreeStateDB *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 163 "wise2xhmmer2.dy"
boolean open_for_indexing_ThreeStateDB_HMMer2(ThreeStateDB * tdb)
{
  return TRUE; /* nowt to do ! */
}

/* Function:  close_for_indexing_ThreeStateDB_HMMer2(tdb)
 *
 * Descrip:    This function is the generic wrapper for the
 *             threestatedb open for indexing 
 *
 *
 * Arg:        tdb [UNKN ] Undocumented argument [ThreeStateDB *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 173 "wise2xhmmer2.dy"
boolean close_for_indexing_ThreeStateDB_HMMer2(ThreeStateDB * tdb)
{
  return TRUE; /* nowt to do ! */
}


#ifdef HMMER_2_1_2 /* new style indexing */

 ThreeStateModel * index_ThreeStateDB_HMMer2(ThreeStateDB * tdb,DataEntry *de)
{
  struct plan7_s *hmm;
  ThreeStateModel * out;
  HMMFILE * hmmfp;

  hmmfp = HMMFileOpen(tdb->filename,"HMMERDB");
  if( hmmfp == NULL ) {
    warn("Could not open %s as HMMER file\n",tdb->filename);
    return NULL;
  }

  if( hmmfp->gsi == NULL ) {
    warn("Could not retrieve sequences from HMMER database as it is not indexed!");
    return NULL;
  }

  if( ! HMMFilePositionByName(hmmfp,de->name) ) {
    warn("Unable to find HMM %s in %s",de->name,tdb->filename);
    return NULL;
  }

  if( HMMFileRead(hmmfp, &hmm) == 0 ) {
    return NULL;
  }

  if( hmm == NULL ) {
    warn("Unable to read HMM file... ");
    return NULL;
  }
  
  out = ThreeStateModel_from_HMMer2(hmm);
  FreePlan7(hmm);
  HMMFileClose(hmmfp);

  return out;
}
#else /* use old style indexing */

/* Function:  index_ThreeStateDB_HMMer2(tdb,*de)
 *
 * Descrip:    This function is the generic wrapper for the
 *             threestatedb get indexed
 *
 *
 * Arg:        tdb [UNKN ] Undocumented argument [ThreeStateDB *]
 * Arg:        *de [UNKN ] Undocumented argument [DataEntry]
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateModel*]
 *
 */
# line 225 "wise2xhmmer2.dy"
ThreeStateModel*  index_ThreeStateDB_HMMer2(ThreeStateDB * tdb,DataEntry *de)
{
  struct plan7_s *hmm;
  ThreeStateModel * out;
  HMMFILE * hmmfp;


  hmmfp = HMMFileOpenFseek(tdb->filename,"HMMERDB",de->byte_position);

  if( hmmfp == NULL ) {
    warn("Could not open %s to byte position %d",tdb->filename,de->data[0]);
    return NULL;
  }

  if( HMMFileRead(hmmfp, &hmm) == 0 ) {
    return NULL;
  }

  if( hmm == NULL ) {
    warn("Unable to read HMM file... ");
    return NULL;
  }
  
  out = ThreeStateModel_from_HMMer2(hmm);
  FreePlan7(hmm);
  HMMFileClose(hmmfp);

  return out;
}
#endif /* HMMER_2_1_2 */

  
/* Function:  dataentry_ThreeStateDB_HMMer2(tdb,en)
 *
 * Descrip:    This function is the generic wrapper
 *             for the threestatedb info for hmmer2
 *
 *
 * Arg:        tdb [UNKN ] Undocumented argument [ThreeStateDB *]
 * Arg:         en [UNKN ] Undocumented argument [DataEntry *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 262 "wise2xhmmer2.dy"
boolean dataentry_ThreeStateDB_HMMer2(ThreeStateDB * tdb,DataEntry * en)
{
  HMMFILE * hmmfp;

  hmmfp = (HMMFILE *) (tdb->data);

#ifdef HMMER_2_1_2
  en->data[0] = 0;
  /* we rely on the fact that someone has already set the name */

#else
  en->data[0] = HMMFtell(hmmfp); 
  en->byte_position = tdb->byte_position; 
#endif

  return TRUE;
}


/* Function:  close_ThreeStateDB_HMMer2(tdb)
 *
 * Descrip:    This function is the generic wrapper
 *             for the threestatedb close for hmmer2
 *
 *
 * Arg:        tdb [UNKN ] Undocumented argument [ThreeStateDB *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 286 "wise2xhmmer2.dy"
boolean close_ThreeStateDB_HMMer2(ThreeStateDB * tdb)
{
  HMMFILE * hmmfp;

  hmmfp = (HMMFILE *) (tdb->data);

  HMMFileClose(hmmfp);
  return TRUE;
}


/* Function:  reload_ThreeStateDB_HMMer2(tdb,return_status)
 *
 * Descrip:    This function is the generic wrapper
 *             for the threestatedb relaod function for hmmer2
 *
 *
 * Arg:                  tdb [UNKN ] Undocumented argument [ThreeStateDB *]
 * Arg:        return_status [UNKN ] Undocumented argument [int *]
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateModel *]
 *
 */
# line 302 "wise2xhmmer2.dy"
ThreeStateModel * reload_ThreeStateDB_HMMer2(ThreeStateDB * tdb,int * return_status)
{
  struct plan7_s *hmm;
  ThreeStateModel * out;
  HMMFILE * hmmfp;

  hmmfp = (HMMFILE *) (tdb->data);

#ifndef HMMER_2_1_2
  tdb->byte_position = HMMFtell(hmmfp); 
#endif

  if( HMMFileRead(hmmfp, &hmm) == 0 ) {
    /* end of file */
    *return_status = DB_RETURN_END;
    return NULL;
  }

  if( hmm == NULL ) {
    warn("Unable to read HMM file... ");
    *return_status = DB_RETURN_ERROR;
  }
  
  *return_status = DB_RETURN_OK;

  out = ThreeStateModel_from_HMMer2(hmm);
  FreePlan7(hmm);

  return out;
}

  

/* Function:  open_ThreeStateDB_HMMer2(tdb)
 *
 * Descrip:    This function is the generic wrapper
 *             for the threestatedb open function
 *
 *
 * Arg:        tdb [UNKN ] Undocumented argument [ThreeStateDB *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 340 "wise2xhmmer2.dy"
boolean open_ThreeStateDB_HMMer2(ThreeStateDB * tdb) 
{
  if( tdb->data != NULL ) {
    warn("Attempting to open an already open tdb. Ugh!");
    return FALSE;
  }

  if( (tdb->data = (void *)HMMFileOpen(tdb->filename,"HMMERDB")) == NULL ) {
    warn("Could not open %s as a filename for hmmer2 directory",tdb->filename);
    return FALSE;
  }
  tdb->byte_position = 0;

  return TRUE;
}
  



/* Function:  ThreeStateModel_from_HMMer2(p7)
 *
 * Descrip:    This function converts a HMMer2 HMM into
 *             a Wise2 HMM (threestatemodel). 
 *
 *
 * Arg:        p7 [UNKN ] Undocumented argument [struct plan7_s *]
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateModel *]
 *
 */
# line 363 "wise2xhmmer2.dy"
ThreeStateModel * ThreeStateModel_from_HMMer2(struct plan7_s * p7)
{
  ThreeStateModel * out;
  ThreeStateUnit  * unit;
  int i,j;
  char * runner;

  if( p7 == NULL ) {
    warn("Can't make a ThreeStateModel from a NULL pointer... ");
    return NULL;
  }

  out = ThreeStateModel_alloc_len(p7->M);

  if( out == NULL ) 
    return NULL; /* already complained*/

  out->name = stringalloc(p7->name);


#ifdef HMMER_2_1_2
  if( p7->acc != NULL ) {
    out->accession = stringalloc(p7->acc);
  } else {
    out->accession = stringalloc(p7->name);
  } 
#else
  /* munge around accession number. Ugh! */

  if( p7->desc != NULL ) {
    runner = strtok(p7->desc,spacestr);
    if( runner != NULL && strstartcmp(runner,"PF") == 0 ) {
      out->accession = stringalloc(runner);
    } else {
      out->accession = stringalloc(p7->name);
    }
  }
#endif

  for(i=0;i<p7->M;i++) {

    unit = blank_ThreeStateUnit();
    add_ThreeStateModel(out,unit);


    if( i != p7->M-1) {
      unit->transition[TSM_MATCH2MATCH]= p7->t[i+1][TMM];
      unit->transition[TSM_MATCH2DELETE]= p7->t[i+1][TMD];
      unit->transition[TSM_MATCH2INSERT]= p7->t[i+1][TMI];
      unit->transition[TSM_INSERT2MATCH]= p7->t[i+1][TIM];
      unit->transition[TSM_INSERT2INSERT]= p7->t[i+1][TII];
      unit->transition[TSM_DELETE2MATCH]= p7->t[i+1][TDM];
      unit->transition[TSM_DELETE2DELETE]= p7->t[i+1][TDD];
      unit->transition[TSM_START2MATCH] = p7->begin[i+1];
      unit->transition[TSM_MATCH2END] =   p7->end[i+1];
    }


    if( i == p7->M-1) {
      /* all probability onto end */

      unit->transition[TSM_MATCH2END] =  1.0;
    }

    /* we have to map sean's probabilities over to ours */

    for(j=0;j<Alphabet_size;j++) {
      unit->match_emission[Alphabet[j]-'A'] = p7->mat[i+1][j];
    }

    if( i != p7->M-1) {
      for(j=0;j<Alphabet_size;j++) {
	unit->insert_emission[Alphabet[j]-'A'] = p7->ins[i+1][j];
      }
    } else {
      /* no insert position! - could leave as 0's? */
      /* but this does upset some sensible error tracking code */
      /* set it to 1/26 probs ;) */
      for(j=0;j<Alphabet_size;j++) {
	unit->insert_emission[Alphabet[j]-'A'] = p7->ins[1][j];
      }

    }
  }

  out->rm = RandomModel_alloc();

  for(i=0;i<26;i++)
    out->rm->aminoacid[i] = 1.0;

  for(j=0;j<Alphabet_size;j++) {
    out->rm->aminoacid[Alphabet[j]-'A'] = p7->null[j];
  }

  display_char_in_ThreeStateModel(out);

  return out;
}
  

  




# line 499 "wise2xhmmer2.c"

#ifdef _cplusplus
}
#endif
