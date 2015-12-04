#ifdef _cplusplus
extern "C" {
#endif
#include "gwquickdb.h"


/* Function:  init_GeneWiseQuickDB(gdb,return_status)
 *
 * Descrip:    inits a genewisequick database. 
 *
 *
 * Arg:                  gdb [UNKN ] Undocumented argument [GeneWiseQuickDB *]
 * Arg:        return_status [UNKN ] Undocumented argument [int *]
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseScoreFlat *]
 *
 */
# line 37 "gwquickdb.dy"
GeneWiseScoreFlat * init_GeneWiseQuickDB(GeneWiseQuickDB * gdb,int * return_status)
{
  GeneWiseScore * gws;
  
  gws = init_GeneWiseDB(gdb->gwdb,return_status);
  if( gws != NULL ) {
    return GeneWiseScoreFlat_from_GeneWiseScore(gws);
  }

}


/* Function:  reload_GeneWiseQuickDB(prev,gdb,return_status)
 *
 * Descrip:    Reloads a genewisequick database
 *
 *
 *
 * Arg:                 prev [UNKN ] Undocumented argument [GeneWiseScoreFlat *]
 * Arg:                  gdb [UNKN ] Undocumented argument [GeneWiseQuickDB *]
 * Arg:        return_status [UNKN ] Undocumented argument [int *]
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseScoreFlat *]
 *
 */
# line 53 "gwquickdb.dy"
GeneWiseScoreFlat * reload_GeneWiseQuickDB(GeneWiseScoreFlat * prev,GeneWiseQuickDB * gdb,int * return_status)
{
  GeneWiseScore * gws;

  if( prev != NULL )
    free_GeneWiseScoreFlat(prev);

  gws = reload_GeneWiseDB(NULL,gdb->gwdb,return_status);
  if( gws != NULL ) {
    return GeneWiseScoreFlat_from_GeneWiseScore(gws);
  }

}


/* Function:  close_GeneWiseQuickDB(gws,gdb)
 *
 * Descrip:    closes a GeneWiseDB
 *
 *
 * Arg:        gws [UNKN ] Undocumented argument [GeneWiseScoreFlat *]
 * Arg:        gdb [UNKN ] Undocumented argument [GeneWiseQuickDB *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 71 "gwquickdb.dy"
boolean close_GeneWiseQuickDB(GeneWiseScoreFlat * gws,GeneWiseQuickDB * gdb)
{
  if( gws != NULL )
    free_GeneWiseScoreFlat(gws);

  return close_GeneWiseDB(NULL,gdb->gwdb);
}

/* Function:  dataentry_add_GeneWiseQuickDB(de,gws,gdb)
 *
 * Descrip:    adds dataentry stuff to a query.
 *
 *
 * Arg:         de [UNKN ] Undocumented argument [DataEntry *]
 * Arg:        gws [UNKN ] Undocumented argument [GeneWiseScoreFlat *]
 * Arg:        gdb [UNKN ] Undocumented argument [GeneWiseQuickDB *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 82 "gwquickdb.dy"
boolean dataentry_add_GeneWiseQuickDB(DataEntry * de,GeneWiseScoreFlat * gws,GeneWiseQuickDB * gdb)
{
  return dataentry_add_GeneWiseDB(de,NULL,gdb->gwdb);
}

/* Function:  GeneWiseQuickDB_from_GeneWiseDB(gwdb)
 *
 * Descrip:    Makes a new genewisequickdb from a genewisemodeldb
 *
 *
 * Arg:        gwdb [READ ] genewisedb - hard links as it enters [GeneWiseDB *]
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseQuickDB *]
 *
 */
# line 92 "gwquickdb.dy"
GeneWiseQuickDB * GeneWiseQuickDB_from_GeneWiseDB(GeneWiseDB * gwdb)
{
  GeneWiseQuickDB * out;

  out = GeneWiseQuickDB_alloc();
  out->gwdb = hard_link_GeneWiseDB(gwdb);

  return out;
}

# line 108 "gwquickdb.c"
/* Function:  hard_link_GeneWiseQuickDB(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GeneWiseQuickDB *]
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseQuickDB *]
 *
 */
GeneWiseQuickDB * hard_link_GeneWiseQuickDB(GeneWiseQuickDB * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a GeneWiseQuickDB object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  GeneWiseQuickDB_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseQuickDB *]
 *
 */
GeneWiseQuickDB * GeneWiseQuickDB_alloc(void) 
{
    GeneWiseQuickDB * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(GeneWiseQuickDB *) ckalloc (sizeof(GeneWiseQuickDB))) == NULL)  {  
      warn("GeneWiseQuickDB_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->gwdb = NULL;    


    return out;  
}    


/* Function:  free_GeneWiseQuickDB(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GeneWiseQuickDB *]
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseQuickDB *]
 *
 */
GeneWiseQuickDB * free_GeneWiseQuickDB(GeneWiseQuickDB * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a GeneWiseQuickDB obj. Should be trappable");   
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
    if( obj->gwdb != NULL)   
      free_GeneWiseDB(obj->gwdb);    


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif
