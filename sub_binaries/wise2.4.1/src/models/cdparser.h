#ifndef DYNAMITEcdparserHEADERFILE
#define DYNAMITEcdparserHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "wisebase.h"
#include "probability.h"


enum cDNAParserTrans {
  PCD_INSERT_2_BASE = 0,
  PCD_INSERT_1_BASE,
  PCD_DELETE_2_BASE,
  PCD_DELETE_1_BASE,
  PCD_PARSER_TRANS_LEN };

/* Object cDNAParser
 *
 * Descrip: This object holds the (very few) extra
 *        transition information needed for the
 *        estwise algorithm. It is sort of like
 *        the 'gene model' part of sequencing
 *        error (but very very simple)
 *
 *
 *
 */
struct Wise2_cDNAParser {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    Probability trans[PCD_PARSER_TRANS_LEN];     
    } ;  
/* cDNAParser defined */ 
#ifndef DYNAMITE_DEFINED_cDNAParser
typedef struct Wise2_cDNAParser Wise2_cDNAParser;
#define cDNAParser Wise2_cDNAParser
#define DYNAMITE_DEFINED_cDNAParser
#endif


/* Object cDNAParserScore
 *
 * Descrip: This object is the score counter
 *        point to cDNAParser (which is
 *        in probabilities). 
 *
 *
 */
struct Wise2_cDNAParserScore {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    Score trans[PCD_PARSER_TRANS_LEN];   
    } ;  
/* cDNAParserScore defined */ 
#ifndef DYNAMITE_DEFINED_cDNAParserScore
typedef struct Wise2_cDNAParserScore Wise2_cDNAParserScore;
#define cDNAParserScore Wise2_cDNAParserScore
#define DYNAMITE_DEFINED_cDNAParserScore
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  removed_probability_from_cds_cdna(cdp)
 *
 * Descrip:    Makes a convienient sum over all the transition
 *             probabilities
 *
 *
 * Arg:        cdp [UNKN ] Undocumented argument [cDNAParser *]
 *
 * Return [UNKN ]  Undocumented return value [Probability]
 *
 */
Probability Wise2_removed_probability_from_cds_cdna(cDNAParser * cdp);
#define removed_probability_from_cds_cdna Wise2_removed_probability_from_cds_cdna


/* Function:  cDNAParserScore_from_cDNAParser(cdp)
 *
 * Descrip:    Makes a new Score object from its probability
 *             counterpart
 *
 *
 * Arg:        cdp [UNKN ] Undocumented argument [cDNAParser *]
 *
 * Return [UNKN ]  Undocumented return value [cDNAParserScore *]
 *
 */
cDNAParserScore * Wise2_cDNAParserScore_from_cDNAParser(cDNAParser * cdp);
#define cDNAParserScore_from_cDNAParser Wise2_cDNAParserScore_from_cDNAParser


/* Function:  flat_cDNAParser(p)
 *
 * Descrip:    Makes a flat (ie, indels of 1 or 2 == p)
 *             cDNA parser. This means that insertions
 *             and deletions of both 1 or 2 bases are
 *             all parameterised at the same probability
 *
 *
 *
 * Arg:        p [READ ] probability of an indel [Probability]
 *
 * Return [UNKN ]  Undocumented return value [cDNAParser *]
 *
 */
cDNAParser * Wise2_flat_cDNAParser(Probability p);
#define flat_cDNAParser Wise2_flat_cDNAParser


/* Function:  hard_link_cDNAParser(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [cDNAParser *]
 *
 * Return [UNKN ]  Undocumented return value [cDNAParser *]
 *
 */
cDNAParser * Wise2_hard_link_cDNAParser(cDNAParser * obj);
#define hard_link_cDNAParser Wise2_hard_link_cDNAParser


/* Function:  cDNAParser_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [cDNAParser *]
 *
 */
cDNAParser * Wise2_cDNAParser_alloc(void);
#define cDNAParser_alloc Wise2_cDNAParser_alloc


/* Function:  free_cDNAParser(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [cDNAParser *]
 *
 * Return [UNKN ]  Undocumented return value [cDNAParser *]
 *
 */
cDNAParser * Wise2_free_cDNAParser(cDNAParser * obj);
#define free_cDNAParser Wise2_free_cDNAParser


/* Function:  hard_link_cDNAParserScore(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [cDNAParserScore *]
 *
 * Return [UNKN ]  Undocumented return value [cDNAParserScore *]
 *
 */
cDNAParserScore * Wise2_hard_link_cDNAParserScore(cDNAParserScore * obj);
#define hard_link_cDNAParserScore Wise2_hard_link_cDNAParserScore


/* Function:  cDNAParserScore_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [cDNAParserScore *]
 *
 */
cDNAParserScore * Wise2_cDNAParserScore_alloc(void);
#define cDNAParserScore_alloc Wise2_cDNAParserScore_alloc


/* Function:  free_cDNAParserScore(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [cDNAParserScore *]
 *
 * Return [UNKN ]  Undocumented return value [cDNAParserScore *]
 *
 */
cDNAParserScore * Wise2_free_cDNAParserScore(cDNAParserScore * obj);
#define free_cDNAParserScore Wise2_free_cDNAParserScore


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/

#ifdef _cplusplus
}
#endif

#endif
