#ifndef DYNAMITErandomdbHEADERFILE
#define DYNAMITErandomdbHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "wisebase.h"
#include "histogram.h"
#include "randommodel.h"



struct Wise2_RandomProteinDB {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    boolean use_flat_length;     
    int length;  
    Histogram * length_dist;     
    RandomModel * emission;  
    int num;     
    } ;  
/* RandomProteinDB defined */ 
#ifndef DYNAMITE_DEFINED_RandomProteinDB
typedef struct Wise2_RandomProteinDB Wise2_RandomProteinDB;
#define RandomProteinDB Wise2_RandomProteinDB
#define DYNAMITE_DEFINED_RandomProteinDB
#endif


struct Wise2_RandomDNADB {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    boolean use_flat_length;     
    int length;  
    Histogram * length_dist;     
    RandomModelDNA * emission;   
    int num;     
    } ;  
/* RandomDNADB defined */ 
#ifndef DYNAMITE_DEFINED_RandomDNADB
typedef struct Wise2_RandomDNADB Wise2_RandomDNADB;
#define RandomDNADB Wise2_RandomDNADB
#define DYNAMITE_DEFINED_RandomDNADB
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  Sequence_from_RandomProteinDB(rndp)
 *
 * Descrip:    Makes a new random sequence from a
 *             RandomProteinDB
 *
 *
 * Arg:        rndp [UNKN ] Undocumented argument [RandomProteinDB *]
 *
 * Return [UNKN ]  Undocumented return value [Sequence *]
 *
 */
Sequence * Wise2_Sequence_from_RandomProteinDB(RandomProteinDB * rndp);
#define Sequence_from_RandomProteinDB Wise2_Sequence_from_RandomProteinDB


/* Function:  new_flat_RandomProteinDB(rm,length)
 *
 * Descrip:    Makes a new flat RandomProteinDB
 *
 *
 * Arg:            rm [SOFT ] RandomModel [RandomModel *]
 * Arg:        length [UNKN ] length of protein to produce [int]
 *
 * Return [UNKN ]  Undocumented return value [RandomProteinDB *]
 *
 */
RandomProteinDB * Wise2_new_flat_RandomProteinDB(RandomModel * rm,int length);
#define new_flat_RandomProteinDB Wise2_new_flat_RandomProteinDB


/* Function:  Sequence_from_RandomDNADB(rndp)
 *
 * Descrip:    Makes a new random sequence from a
 *             RandomDNADB
 *
 *
 * Arg:        rndp [UNKN ] Undocumented argument [RandomDNADB *]
 *
 * Return [UNKN ]  Undocumented return value [Sequence *]
 *
 */
Sequence * Wise2_Sequence_from_RandomDNADB(RandomDNADB * rndp);
#define Sequence_from_RandomDNADB Wise2_Sequence_from_RandomDNADB


/* Function:  new_flat_RandomDNADB(rm,length)
 *
 * Descrip:    Makes a new flat RandomDNADB
 *
 *
 * Arg:            rm [SOFT ] RandomModel emission [RandomModelDNA *]
 * Arg:        length [UNKN ] length of DNA to produce [int]
 *
 * Return [UNKN ]  Undocumented return value [RandomDNADB *]
 *
 */
RandomDNADB * Wise2_new_flat_RandomDNADB(RandomModelDNA * rm,int length);
#define new_flat_RandomDNADB Wise2_new_flat_RandomDNADB


/* Function:  hard_link_RandomProteinDB(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [RandomProteinDB *]
 *
 * Return [UNKN ]  Undocumented return value [RandomProteinDB *]
 *
 */
RandomProteinDB * Wise2_hard_link_RandomProteinDB(RandomProteinDB * obj);
#define hard_link_RandomProteinDB Wise2_hard_link_RandomProteinDB


/* Function:  RandomProteinDB_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [RandomProteinDB *]
 *
 */
RandomProteinDB * Wise2_RandomProteinDB_alloc(void);
#define RandomProteinDB_alloc Wise2_RandomProteinDB_alloc


/* Function:  free_RandomProteinDB(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [RandomProteinDB *]
 *
 * Return [UNKN ]  Undocumented return value [RandomProteinDB *]
 *
 */
RandomProteinDB * Wise2_free_RandomProteinDB(RandomProteinDB * obj);
#define free_RandomProteinDB Wise2_free_RandomProteinDB


/* Function:  hard_link_RandomDNADB(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [RandomDNADB *]
 *
 * Return [UNKN ]  Undocumented return value [RandomDNADB *]
 *
 */
RandomDNADB * Wise2_hard_link_RandomDNADB(RandomDNADB * obj);
#define hard_link_RandomDNADB Wise2_hard_link_RandomDNADB


/* Function:  RandomDNADB_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [RandomDNADB *]
 *
 */
RandomDNADB * Wise2_RandomDNADB_alloc(void);
#define RandomDNADB_alloc Wise2_RandomDNADB_alloc


/* Function:  free_RandomDNADB(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [RandomDNADB *]
 *
 * Return [UNKN ]  Undocumented return value [RandomDNADB *]
 *
 */
RandomDNADB * Wise2_free_RandomDNADB(RandomDNADB * obj);
#define free_RandomDNADB Wise2_free_RandomDNADB


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
boolean Wise2_replace_length_dist_RandomProteinDB(RandomProteinDB * obj,Histogram * length_dist);
#define replace_length_dist_RandomProteinDB Wise2_replace_length_dist_RandomProteinDB
boolean Wise2_replace_use_flat_length_RandomDNADB(RandomDNADB * obj,boolean use_flat_length);
#define replace_use_flat_length_RandomDNADB Wise2_replace_use_flat_length_RandomDNADB
boolean Wise2_access_use_flat_length_RandomProteinDB(RandomProteinDB * obj);
#define access_use_flat_length_RandomProteinDB Wise2_access_use_flat_length_RandomProteinDB
boolean Wise2_access_use_flat_length_RandomDNADB(RandomDNADB * obj);
#define access_use_flat_length_RandomDNADB Wise2_access_use_flat_length_RandomDNADB
int Wise2_access_length_RandomProteinDB(RandomProteinDB * obj);
#define access_length_RandomProteinDB Wise2_access_length_RandomProteinDB
boolean Wise2_replace_length_RandomDNADB(RandomDNADB * obj,int length);
#define replace_length_RandomDNADB Wise2_replace_length_RandomDNADB
Histogram * Wise2_access_length_dist_RandomProteinDB(RandomProteinDB * obj);
#define access_length_dist_RandomProteinDB Wise2_access_length_dist_RandomProteinDB
int Wise2_access_length_RandomDNADB(RandomDNADB * obj);
#define access_length_RandomDNADB Wise2_access_length_RandomDNADB
RandomModel * Wise2_access_emission_RandomProteinDB(RandomProteinDB * obj);
#define access_emission_RandomProteinDB Wise2_access_emission_RandomProteinDB
boolean Wise2_replace_length_dist_RandomDNADB(RandomDNADB * obj,Histogram * length_dist);
#define replace_length_dist_RandomDNADB Wise2_replace_length_dist_RandomDNADB
int Wise2_access_num_RandomProteinDB(RandomProteinDB * obj);
#define access_num_RandomProteinDB Wise2_access_num_RandomProteinDB
Histogram * Wise2_access_length_dist_RandomDNADB(RandomDNADB * obj);
#define access_length_dist_RandomDNADB Wise2_access_length_dist_RandomDNADB
boolean Wise2_replace_use_flat_length_RandomProteinDB(RandomProteinDB * obj,boolean use_flat_length);
#define replace_use_flat_length_RandomProteinDB Wise2_replace_use_flat_length_RandomProteinDB
boolean Wise2_replace_emission_RandomDNADB(RandomDNADB * obj,RandomModelDNA * emission);
#define replace_emission_RandomDNADB Wise2_replace_emission_RandomDNADB
boolean Wise2_replace_num_RandomProteinDB(RandomProteinDB * obj,int num);
#define replace_num_RandomProteinDB Wise2_replace_num_RandomProteinDB
RandomModelDNA * Wise2_access_emission_RandomDNADB(RandomDNADB * obj);
#define access_emission_RandomDNADB Wise2_access_emission_RandomDNADB
boolean Wise2_replace_emission_RandomProteinDB(RandomProteinDB * obj,RandomModel * emission);
#define replace_emission_RandomProteinDB Wise2_replace_emission_RandomProteinDB
boolean Wise2_replace_num_RandomDNADB(RandomDNADB * obj,int num);
#define replace_num_RandomDNADB Wise2_replace_num_RandomDNADB
boolean Wise2_replace_length_RandomProteinDB(RandomProteinDB * obj,int length);
#define replace_length_RandomProteinDB Wise2_replace_length_RandomProteinDB
int Wise2_access_num_RandomDNADB(RandomDNADB * obj);
#define access_num_RandomDNADB Wise2_access_num_RandomDNADB

#ifdef _cplusplus
}
#endif

#endif
