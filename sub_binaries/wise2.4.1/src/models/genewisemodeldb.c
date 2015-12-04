#ifdef _cplusplus
extern "C" {
#endif
#include "genewisemodeldb.h"



/* Function:  init_GeneWiseDB(gdb,return_status)
 *
 * Descrip:    inits a genewise database. Remember this is used
 *             for both genomic and cdna searches
 *
 *
 *
 * Arg:                  gdb [UNKN ] Undocumented argument [GeneWiseDB *]
 * Arg:        return_status [UNKN ] Undocumented argument [int *]
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseScore *]
 *
 */
# line 47 "genewisemodeldb.dy"
GeneWiseScore * init_GeneWiseDB(GeneWiseDB * gdb,int * return_status)
{
  ThreeStateModel * tsm;
  GeneWise * gw;
  GeneWiseScore * gws;
 

  *return_status = DB_RETURN_ERROR;

  if( gdb->is_single == TRUE ) {
    *return_status = DB_RETURN_OK;
    return gdb->gws;
  }

  if( open_ThreeStateDB(gdb->tdb) == FALSE ) {
    warn("Could not open three state db, so no genewisemodel possible!");
    return NULL;
  }

  tsm = read_TSM_ThreeStateDB(gdb->tdb,return_status);

  if( *return_status == DB_RETURN_ERROR) {
    warn("Cannot read a ThreeStateModelDB for the GeneWiseDB. Problem!");
    return NULL;
  }

  gw  =  GeneWise_from_ThreeStateModel_GDB(tsm,gdb);
  gws =  GeneWiseScore_from_GeneWise(gw);

  free_ThreeStateModel(tsm);
  free_GeneWise(gw);


  return gws;
}

/* Function:  reload_GeneWiseDB(prev,gdb,return_status)
 *
 * Descrip:    Reloads a genewise database
 *
 *
 *
 * Arg:                 prev [UNKN ] Undocumented argument [GeneWiseScore *]
 * Arg:                  gdb [UNKN ] Undocumented argument [GeneWiseDB *]
 * Arg:        return_status [UNKN ] Undocumented argument [int *]
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseScore *]
 *
 */
# line 87 "genewisemodeldb.dy"
GeneWiseScore * reload_GeneWiseDB(GeneWiseScore * prev,GeneWiseDB * gdb,int * return_status)
{
  ThreeStateModel * tsm;
  GeneWise * gw;
  GeneWiseScore * gws;
 


  if( gdb->is_single == TRUE ) {
    *return_status = DB_RETURN_END;
    return NULL;
  }

  if( prev != NULL)
    free_GeneWiseScore(prev);


  *return_status = DB_RETURN_ERROR;

  tsm = read_TSM_ThreeStateDB(gdb->tdb,return_status);

  if( *return_status != DB_RETURN_OK) {
    if( *return_status == DB_RETURN_END) {
      return NULL;
    }

    warn("Cannot read a ThreeStateModelDB for the GeneWiseDB. Problem!");
    return NULL;
  }

  gw  =  GeneWise_from_ThreeStateModel_GDB(tsm,gdb);
  gws =  GeneWiseScore_from_GeneWise(gw);



  free_ThreeStateModel(tsm);
  free_GeneWise(gw);

  *return_status = DB_RETURN_OK;

  return gws;
}

/* Function:  close_GeneWiseDB(gws,gdb)
 *
 * Descrip:    closes a GeneWiseDB
 *
 *
 * Arg:        gws [UNKN ] Undocumented argument [GeneWiseScore *]
 * Arg:        gdb [UNKN ] Undocumented argument [GeneWiseDB *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 133 "genewisemodeldb.dy"
boolean close_GeneWiseDB(GeneWiseScore * gws,GeneWiseDB * gdb)
{
  if( gdb->is_single == TRUE ) {
    return TRUE;
  }

  if( gws != NULL)
    free_GeneWiseScore(gws);

  return close_ThreeStateDB(NULL,gdb->tdb);
}

/* Function:  dataentry_add_GeneWiseDB(de,gws,gdb)
 *
 * Descrip:    adds dataentry stuff to a query. Relies completely
 *             on threestatedb
 *
 *
 * Arg:         de [UNKN ] Undocumented argument [DataEntry *]
 * Arg:        gws [UNKN ] Undocumented argument [GeneWiseScore *]
 * Arg:        gdb [UNKN ] Undocumented argument [GeneWiseDB *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 149 "genewisemodeldb.dy"
boolean dataentry_add_GeneWiseDB(DataEntry * de,GeneWiseScore * gws,GeneWiseDB * gdb)
{
  if( gdb->is_single == TRUE) {
    if ( gdb->gws->name == NULL ) {
      warn("No name for a single GeneWiseDB unit.");
      de->name = stringalloc("NoName");
    } else {
      de->name = stringalloc(gdb->gws->name);
    }

    return TRUE;
  }

  
  /*de->name = stringalloc(gws->name);*/

  /* otherwise, pass it on to tdb */
  dataentry_add_ThreeStateDB(de,NULL,gdb->tdb);
  if( de->name == NULL ) {
    if( gws->name != NULL ) {
      de->name = stringalloc(gws->name);
    }
  }
  return TRUE;
}

/* Function:  GeneWise_from_ThreeStateModel_GDB(tsm,gdb)
 *
 * Descrip:    Makes a genewise models from the threestatemodel
 *             and parameters held in the GeneWiseDB
 *
 *
 * Arg:        tsm [UNKN ] Undocumented argument [ThreeStateModel *]
 * Arg:        gdb [UNKN ] Undocumented argument [GeneWiseDB *]
 *
 * Return [UNKN ]  Undocumented return value [GeneWise *]
 *
 */
# line 180 "genewisemodeldb.dy"
GeneWise * GeneWise_from_ThreeStateModel_GDB(ThreeStateModel * tsm,GeneWiseDB * gdb)
{
  int i;
  GeneWise * out;


  if( gdb->gpara != NULL ) {
    out=GeneWise_from_ThreeStateModel(tsm,gdb->gpara->gp,gdb->gpara->cm,gdb->allN,gdb->gpara->gwcm);
    if( out == NULL ) {
      return NULL;
    }

    if( gdb->gpara->modelled_splice == FALSE)
      flatten_balance_scores_GeneWise(out);
  } else {
    out = GeneWise_alloc_len(tsm->len);
    out->name = stringalloc(tsm->name);

    for(i=0;i<tsm->len;i++) {
      add_GeneWise(out,GeneWiseSegment_from_ThreeStateUnit(tsm->unit[i],gdb->cds_factor,gdb->cm,NULL,gdb->allN));
    }
  }


  if( gdb->is_syn == TRUE ) {
    if( tsm->rm == NULL ) {
      warn("Bad error. ThreeStateModel does not have a random model!");
      return free_GeneWise(out);
    }
    GeneWise_fold_in_synchronised_RandomModel(out,tsm->rm,gdb->cm,gdb->cm->ct,0.5);
  } else {
    GeneWise_fold_in_RandomModelDNA(out,gdb->rmd);
  }
  
  if( gdb->flat_insert == TRUE ) {
    check_flat_insert(out,1,0,gdb->cm->ct);
  }

  return out;
}

/* Function:  new_GeneWiseDB_cdna(syn,tdb,cp,cm,rmd,use_syn,flat_insert,allN)
 *
 * Descrip:    makes a new GeneWiseDB from its component parts,
 *             assumming a cDNA db.
 *
 *             All the objects are hard-linked internally, so you can, if you
 *             wish, free them once passing them into this function
 *
 *
 * Arg:                syn [UNKN ] if ture, use a synchronous coding model vs internally stored tdb rm's [NullString]
 * Arg:                tdb [UNKN ] three state model db to use [ThreeStateDB *]
 * Arg:                 cp [UNKN ] codon parser function to remove from match state [cDNAParser *]
 * Arg:                 cm [UNKN ] Undocumented argument [CodonMapper *]
 * Arg:                rmd [UNKN ] random model (dna) [RandomModelDNA *]
 * Arg:            use_syn [UNKN ] Undocumented argument [boolean]
 * Arg:        flat_insert [UNKN ] Undocumented argument [boolean]
 * Arg:               allN [UNKN ] Undocumented argument [Probability]
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseDB *]
 *
 */
# line 233 "genewisemodeldb.dy"
GeneWiseDB * new_GeneWiseDB_cdna(ThreeStateDB * tdb,cDNAParser * cp,CodonMapper * cm,RandomModelDNA * rmd,boolean use_syn,boolean flat_insert,Probability allN) 
{
  GeneWiseDB * out;

  out = GeneWiseDB_alloc();

  out->tdb = hard_link_ThreeStateDB(tdb);
  out->is_single = FALSE;
  
  out->rmd = hard_link_RandomModelDNA(rmd);
  out->cm = hard_link_CodonMapper(cm);
  out->cds_factor = (1.0 - removed_probability_from_cds_cdna(cp));
  out->is_syn = use_syn;
  out->allN = allN;
  out->flat_insert = flat_insert;
  return out;
}

/* Function:  new_GeneWiseDB(tdb,gp,rmd,use_syn,allN)
 *
 * Descrip:    makes a new GeneWiseDB from its component parts.
 *             All the objects are hard-linked internally, so you can, if you
 *             wish, free them once passing them into this function
 *
 *
 * Arg:            tdb [UNKN ] Undocumented argument [ThreeStateDB *]
 * Arg:             gp [UNKN ] Undocumented argument [GeneParameter21 *]
 * Arg:            rmd [UNKN ] Undocumented argument [RandomModelDNA *]
 * Arg:        use_syn [UNKN ] Undocumented argument [boolean]
 * Arg:           allN [UNKN ] Undocumented argument [Probability]
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseDB *]
 *
 */
# line 256 "genewisemodeldb.dy"
GeneWiseDB * new_GeneWiseDB(ThreeStateDB * tdb,GeneParameter21 * gp,RandomModelDNA * rmd,boolean use_syn,Probability allN)
{
  GeneWiseDB * out;

  out = GeneWiseDB_alloc();

  out->tdb = hard_link_ThreeStateDB(tdb);
  out->is_single = FALSE;
  
  out->rmd = hard_link_RandomModelDNA(rmd);
  out->gpara = hard_link_GeneParameter21(gp);
  out->cm = hard_link_CodonMapper(gp->cm);
  out->is_syn = use_syn;
  out->allN   = allN;
  return out;
}

/* Function:  new_single_GeneWiseDB(gws)
 *
 * Descrip:    makes a new GeneWiseDB from a single GeneWiseScore.
 *             It hard links it, so you should free it afterwards
 *             in its own scope.
 *
 *
 * Arg:        gws [UNKN ] Undocumented argument [GeneWiseScore *]
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseDB *]
 *
 */
# line 278 "genewisemodeldb.dy"
GeneWiseDB * new_single_GeneWiseDB(GeneWiseScore * gws)
{
  GeneWiseDB * out;

  out = GeneWiseDB_alloc();
  out->gws = hard_link_GeneWiseScore(gws);
  out->is_single = TRUE;

  return out;
}


/* Function:  init_GwLite_GeneWiseDB(gdb,return_status)
 *
 * Descrip:    inits a genewise database for gwlite models
 *
 *
 *
 *
 * Arg:                  gdb [UNKN ] Undocumented argument [GeneWiseDB *]
 * Arg:        return_status [UNKN ] Undocumented argument [int *]
 *
 * Return [UNKN ]  Undocumented return value [GwLiteScore *]
 *
 */
# line 295 "genewisemodeldb.dy"
GwLiteScore * init_GwLite_GeneWiseDB(GeneWiseDB * gdb,int * return_status)
{
  ThreeStateModel * tsm;
  GeneWise * gw;
  GwLite * gl;
  GwLiteScore * gws;
 

  *return_status = DB_RETURN_ERROR;

  if( gdb->is_single == TRUE ) {
    *return_status = DB_RETURN_OK;
    /* at the moment, die horribly */
    fatal("Apologies. We don't handle single genewise model databases for gwlite yet");
    return gdb->gwls;
  }

  if( open_ThreeStateDB(gdb->tdb) == FALSE ) {
    warn("Could not open three state db, so no genewisemodel possible!");
    return NULL;
  }

  tsm = read_TSM_ThreeStateDB(gdb->tdb,return_status);

  if( *return_status == DB_RETURN_ERROR) {
    warn("Cannot read a ThreeStateModelDB for the GeneWiseDB. Problem!");
    return NULL;
  }

  gw  =  GeneWise_from_ThreeStateModel_GDB(tsm,gdb);
  gl  =  GwLite_from_GeneWise(gw);
  gws =  GwLiteScore_from_GwLite(gl);

  free_ThreeStateModel(tsm);
  free_GeneWise(gw);
  free_GwLite(gl);

  return gws;
}

/* Function:  reload_GwLite_GeneWiseDB(prev,gdb,return_status)
 *
 * Descrip:    Reloads a genewise database for a GwLite database
 *
 *
 *
 * Arg:                 prev [UNKN ] Undocumented argument [GwLiteScore *]
 * Arg:                  gdb [UNKN ] Undocumented argument [GeneWiseDB *]
 * Arg:        return_status [UNKN ] Undocumented argument [int *]
 *
 * Return [UNKN ]  Undocumented return value [GwLiteScore *]
 *
 */
# line 339 "genewisemodeldb.dy"
GwLiteScore * reload_GwLite_GeneWiseDB(GwLiteScore * prev,GeneWiseDB * gdb,int * return_status)
{
  ThreeStateModel * tsm;
  GeneWise * gw;
  GwLite * gl;
  GwLiteScore * gws;
 


  if( gdb->is_single == TRUE ) {
    *return_status = DB_RETURN_END;
    return NULL;
  }

  if( prev != NULL)
    free_GwLiteScore(prev);


  *return_status = DB_RETURN_ERROR;

  tsm = read_TSM_ThreeStateDB(gdb->tdb,return_status);

  if( *return_status != DB_RETURN_OK) {
    if( *return_status == DB_RETURN_END) {
      return NULL;
    }

    warn("Cannot read a ThreeStateModelDB for the GeneWiseDB. Problem!");
    return NULL;
  }

  gw  =  GeneWise_from_ThreeStateModel_GDB(tsm,gdb);
  gl  =  GwLite_from_GeneWise(gw);
  gws =  GwLiteScore_from_GwLite(gl);

  free_ThreeStateModel(tsm);
  free_GeneWise(gw);
  free_GwLite(gl);

  *return_status = DB_RETURN_OK;

  return gws;
}

/* Function:  close_GwLite_GeneWiseDB(gws,gdb)
 *
 * Descrip:    closes a GeneWiseDB
 *
 *
 * Arg:        gws [UNKN ] Undocumented argument [GwLiteScore *]
 * Arg:        gdb [UNKN ] Undocumented argument [GeneWiseDB *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 386 "genewisemodeldb.dy"
boolean close_GwLite_GeneWiseDB(GwLiteScore * gws,GeneWiseDB * gdb)
{
  if( gdb->is_single == TRUE ) {
    return TRUE;
  }

  if( gws != NULL)
    free_GwLiteScore(gws);

  return close_ThreeStateDB(NULL,gdb->tdb);
}

/* Function:  dataentry_add_GwLite_GeneWiseDB(de,gws,gdb)
 *
 * Descrip:    adds dataentry stuff to a query. Relies completely
 *             on threestatedb
 *
 *
 * Arg:         de [UNKN ] Undocumented argument [DataEntry *]
 * Arg:        gws [UNKN ] Undocumented argument [GwLiteScore *]
 * Arg:        gdb [UNKN ] Undocumented argument [GeneWiseDB *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 402 "genewisemodeldb.dy"
boolean dataentry_add_GwLite_GeneWiseDB(DataEntry * de,GwLiteScore * gws,GeneWiseDB * gdb)
{
  if( gdb->is_single == TRUE) {
    if ( gdb->gws->name == NULL ) {
      warn("No name for a single GeneWiseDB unit.");
      de->name = stringalloc("NoName");
    } else {
      de->name = stringalloc(gdb->gws->name);
    }

    return TRUE;
  }

  
  /*de->name = stringalloc(gws->name);*/

  /* otherwise, pass it on to tdb */
  dataentry_add_ThreeStateDB(de,NULL,gdb->tdb);
  if( de->name == NULL ) {
    if( gws->name != NULL ) {
      de->name = stringalloc(gws->name);
    }
  }
  return TRUE;
}

# line 483 "genewisemodeldb.c"
/* Function:  hard_link_GeneWiseDB(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GeneWiseDB *]
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseDB *]
 *
 */
GeneWiseDB * hard_link_GeneWiseDB(GeneWiseDB * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a GeneWiseDB object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  GeneWiseDB_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseDB *]
 *
 */
GeneWiseDB * GeneWiseDB_alloc(void) 
{
    GeneWiseDB * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(GeneWiseDB *) ckalloc (sizeof(GeneWiseDB))) == NULL)    {  
      warn("GeneWiseDB_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->tdb = NULL; 
    out->cds_factor = 0; 
    out->cm = NULL;  
    out->rmd = NULL; 
    out->gws = NULL; 
    out->gpara = NULL;   
    out->is_single = FALSE;  
    out->is_syn = TRUE;  
    out->allN = 1.0; 
    out->flat_insert = FALSE;    
    out->gwls = NULL;    


    return out;  
}    


/* Function:  free_GeneWiseDB(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GeneWiseDB *]
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseDB *]
 *
 */
GeneWiseDB * free_GeneWiseDB(GeneWiseDB * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a GeneWiseDB obj. Should be trappable");    
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
    if( obj->tdb != NULL)    
      free_ThreeStateDB(obj->tdb);   
    if( obj->cm != NULL) 
      free_CodonMapper(obj->cm);     
    if( obj->rmd != NULL)    
      free_RandomModelDNA(obj->rmd);     
    if( obj->gws != NULL)    
      free_GeneWiseScore(obj->gws);  
    if( obj->gpara != NULL)  
      free_GeneParameter21(obj->gpara);  
    if( obj->gwls != NULL)   
      free_GwLiteScore(obj->gwls);   


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif
