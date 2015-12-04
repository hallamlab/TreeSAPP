#ifndef DYNAMITEgeneoutputHEADERFILE
#define DYNAMITEgeneoutputHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "dyna.h"
#include "geneutil.h"



struct Wise2_GeneOutputData {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    AlnBlock * alb;  
    PackAln  * pal;  
    GenomicRegion * gr;  
    Genomic * gen;   
    CodonTable * ct;     
    } ;  
/* GeneOutputData defined */ 
#ifndef DYNAMITE_DEFINED_GeneOutputData
typedef struct Wise2_GeneOutputData Wise2_GeneOutputData;
#define GeneOutputData Wise2_GeneOutputData
#define DYNAMITE_DEFINED_GeneOutputData
#endif


struct Wise2_GeneOutputPara {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    boolean show_genes;  
    boolean show_gff;    
    boolean show_trans;  
    boolean show_cdna;   
    boolean show_geneutr;    
    boolean show_alb;    
    boolean show_pal;    
    boolean show_debug;  
    char * divide_string;    
    } ;  
/* GeneOutputPara defined */ 
#ifndef DYNAMITE_DEFINED_GeneOutputPara
typedef struct Wise2_GeneOutputPara Wise2_GeneOutputPara;
#define GeneOutputPara Wise2_GeneOutputPara
#define DYNAMITE_DEFINED_GeneOutputPara
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  hard_link_GeneOutputData(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GeneOutputData *]
 *
 * Return [UNKN ]  Undocumented return value [GeneOutputData *]
 *
 */
GeneOutputData * Wise2_hard_link_GeneOutputData(GeneOutputData * obj);
#define hard_link_GeneOutputData Wise2_hard_link_GeneOutputData


/* Function:  GeneOutputData_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GeneOutputData *]
 *
 */
GeneOutputData * Wise2_GeneOutputData_alloc(void);
#define GeneOutputData_alloc Wise2_GeneOutputData_alloc


/* Function:  free_GeneOutputData(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GeneOutputData *]
 *
 * Return [UNKN ]  Undocumented return value [GeneOutputData *]
 *
 */
GeneOutputData * Wise2_free_GeneOutputData(GeneOutputData * obj);
#define free_GeneOutputData Wise2_free_GeneOutputData


/* Function:  hard_link_GeneOutputPara(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GeneOutputPara *]
 *
 * Return [UNKN ]  Undocumented return value [GeneOutputPara *]
 *
 */
GeneOutputPara * Wise2_hard_link_GeneOutputPara(GeneOutputPara * obj);
#define hard_link_GeneOutputPara Wise2_hard_link_GeneOutputPara


/* Function:  GeneOutputPara_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GeneOutputPara *]
 *
 */
GeneOutputPara * Wise2_GeneOutputPara_alloc(void);
#define GeneOutputPara_alloc Wise2_GeneOutputPara_alloc


/* Function:  free_GeneOutputPara(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GeneOutputPara *]
 *
 * Return [UNKN ]  Undocumented return value [GeneOutputPara *]
 *
 */
GeneOutputPara * Wise2_free_GeneOutputPara(GeneOutputPara * obj);
#define free_GeneOutputPara Wise2_free_GeneOutputPara


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
void Wise2_show_GeneOutput(GeneOutputData * data,GeneOutputPara * para,FILE * ofp);
#define show_GeneOutput Wise2_show_GeneOutput
double Wise2_id_map_func(int i);
#define id_map_func Wise2_id_map_func
GeneOutputPara * Wise2_new_GeneOutputPara_from_argv(int * argc,char ** argv);
#define new_GeneOutputPara_from_argv Wise2_new_GeneOutputPara_from_argv
void Wise2_show_help_GeneOutputPara(FILE * ofp);
#define show_help_GeneOutputPara Wise2_show_help_GeneOutputPara


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/

#ifdef _cplusplus
}
#endif

#endif
