#ifndef DYNAMITErandommodelHEADERFILE
#define DYNAMITErandommodelHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif

#include "wisebase.h"
#include "probability.h"
#include "codon.h"
#include "sequence.h"


struct Wise2_RandomModel {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    Probability aminoacid[26];   
    char * name;     
    } ;  
/* RandomModel defined */ 
#ifndef DYNAMITE_DEFINED_RandomModel
typedef struct Wise2_RandomModel Wise2_RandomModel;
#define RandomModel Wise2_RandomModel
#define DYNAMITE_DEFINED_RandomModel
#endif


struct Wise2_RandomModelScoreaa {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    Score aminoacid[26];     
    char * name;     
    } ;  
/* RandomModelScoreaa defined */ 
#ifndef DYNAMITE_DEFINED_RandomModelScoreaa
typedef struct Wise2_RandomModelScoreaa Wise2_RandomModelScoreaa;
#define RandomModelScoreaa Wise2_RandomModelScoreaa
#define DYNAMITE_DEFINED_RandomModelScoreaa
#endif


struct Wise2_RandomCodonScore {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    Score codon[126];    
    char * name;     
    } ;  
/* RandomCodonScore defined */ 
#ifndef DYNAMITE_DEFINED_RandomCodonScore
typedef struct Wise2_RandomCodonScore Wise2_RandomCodonScore;
#define RandomCodonScore Wise2_RandomCodonScore
#define DYNAMITE_DEFINED_RandomCodonScore
#endif


struct Wise2_RandomCodon {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    Probability codon[126];  
    char * name;     
    } ;  
/* RandomCodon defined */ 
#ifndef DYNAMITE_DEFINED_RandomCodon
typedef struct Wise2_RandomCodon Wise2_RandomCodon;
#define RandomCodon Wise2_RandomCodon
#define DYNAMITE_DEFINED_RandomCodon
#endif


struct Wise2_RandomModelDNA {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    Probability base[5];     
    char * name;     
    } ;  
/* RandomModelDNA defined */ 
#ifndef DYNAMITE_DEFINED_RandomModelDNA
typedef struct Wise2_RandomModelDNA Wise2_RandomModelDNA;
#define RandomModelDNA Wise2_RandomModelDNA
#define DYNAMITE_DEFINED_RandomModelDNA
#endif


struct Wise2_RandomModelDNAScore {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    Score base[5];   
    char * name;     
    } ;  
/* RandomModelDNAScore defined */ 
#ifndef DYNAMITE_DEFINED_RandomModelDNAScore
typedef struct Wise2_RandomModelDNAScore Wise2_RandomModelDNAScore;
#define RandomModelDNAScore Wise2_RandomModelDNAScore
#define DYNAMITE_DEFINED_RandomModelDNAScore
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  draw_random_aa_RandomModel(rm)
 *
 * Descrip:    Draws an amino acid from the random distribution
 *
 *
 * Arg:        rm [UNKN ] Undocumented argument [RandomModel *]
 *
 * Return [UNKN ]  Undocumented return value [char]
 *
 */
char Wise2_draw_random_aa_RandomModel(RandomModel * rm);
#define draw_random_aa_RandomModel Wise2_draw_random_aa_RandomModel


/* Function:  draw_random_base_RandomModelDNA(rm)
 *
 * Descrip:    Draws a base from the random distribution
 *
 *
 * Arg:        rm [UNKN ] Undocumented argument [RandomModelDNA *]
 *
 * Return [UNKN ]  Undocumented return value [char]
 *
 */
char Wise2_draw_random_base_RandomModelDNA(RandomModelDNA * rm);
#define draw_random_base_RandomModelDNA Wise2_draw_random_base_RandomModelDNA


/* Function:  RandomCodonScore_from_RandomCodon(rc)
 *
 * Descrip:    Makes a score RandomCodon (log space)
 *             from a probability based random codon
 *
 *
 * Arg:        rc [UNKN ] Undocumented argument [RandomCodon *]
 *
 * Return [UNKN ]  Undocumented return value [RandomCodonScore *]
 *
 */
RandomCodonScore * Wise2_RandomCodonScore_from_RandomCodon(RandomCodon * rc);
#define RandomCodonScore_from_RandomCodon Wise2_RandomCodonScore_from_RandomCodon


/* Function:  flatten_RandomCodon(rc)
 *
 * Descrip:    Sets all probabilities to 1.0 - ie,
 *             odds them to themselves.
 *
 *             This is equivalent to saying that the randomcodon
 *             is being odd-ratioed to itself
 *
 *             Also equivalent of saying all the scores (in log
 *             space) will be 0
 *
 *
 * Arg:        rc [UNKN ] Undocumented argument [RandomCodon *]
 *
 */
void Wise2_flatten_RandomCodon(RandomCodon * rc);
#define flatten_RandomCodon Wise2_flatten_RandomCodon


/* Function:  fold_in_RandomModelDNA_into_RandomCodon(rc,rmd)
 *
 * Descrip:    Makes the randomcodon numbers become the odds ratio
 *             between their probabilitys and flat dna random model
 *             (0th order markov)
 *
 *
 * Arg:         rc [UNKN ] Undocumented argument [RandomCodon *]
 * Arg:        rmd [UNKN ] Undocumented argument [RandomModelDNA *]
 *
 */
void Wise2_fold_in_RandomModelDNA_into_RandomCodon(RandomCodon * rc,RandomModelDNA * rmd);
#define fold_in_RandomModelDNA_into_RandomCodon Wise2_fold_in_RandomModelDNA_into_RandomCodon


/* Function:  show_RandomCodonScore(rcs,ofp)
 *
 * Descrip:    shows RandomCodonScore
 *
 *             for debugging
 *
 *
 * Arg:        rcs [UNKN ] Undocumented argument [RandomCodonScore *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
void Wise2_show_RandomCodonScore(RandomCodonScore * rcs,FILE * ofp);
#define show_RandomCodonScore Wise2_show_RandomCodonScore


/* Function:  show_RandomModelDNAScore(rds,ofp)
 *
 * Descrip:    shows RandomModelsDNAScore
 *
 *             for debugging
 *
 *
 * Arg:        rds [UNKN ] Undocumented argument [RandomModelDNAScore *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
void Wise2_show_RandomModelDNAScore(RandomModelDNAScore * rds,FILE * ofp);
#define show_RandomModelDNAScore Wise2_show_RandomModelDNAScore


/* Function:  folded_RandomModelDNAScore_from_2RMD(dis,rnd)
 *
 * Descrip:    gives a odds ratio between two random models
 *
 *
 * Arg:        dis [UNKN ] Undocumented argument [RandomModelDNA *]
 * Arg:        rnd [UNKN ] Undocumented argument [RandomModelDNA *]
 *
 * Return [UNKN ]  Undocumented return value [RandomModelDNAScore *]
 *
 */
RandomModelDNAScore * Wise2_folded_RandomModelDNAScore_from_2RMD(RandomModelDNA * dis,RandomModelDNA * rnd);
#define folded_RandomModelDNAScore_from_2RMD Wise2_folded_RandomModelDNAScore_from_2RMD


/* Function:  RandomCodon_from_raw_CodonFrequency(codon[64],*ct)
 *
 * Descrip:    From raw counts (no adjustment to amino acids) of codons
 *             gives you a RandomCodon model
 *
 *             No prior is used (? perhaps should have a flat prior)
 *
 *             N's are handled correctly
 *
 *
 * Arg:        codon[64] [UNKN ] Undocumented argument [double]
 * Arg:              *ct [UNKN ] Undocumented argument [CodonTable]
 *
 * Return [UNKN ]  Undocumented return value [RandomCodon  *]
 *
 */
RandomCodon  * Wise2_RandomCodon_from_raw_CodonFrequency(double codon[64],CodonTable *ct);
#define RandomCodon_from_raw_CodonFrequency Wise2_RandomCodon_from_raw_CodonFrequency


/* Function:  flat_RandomCodon(ct)
 *
 * Descrip:    Makes a flat RandomCodon from CodonTable
 *
 *
 * Arg:        ct [UNKN ] Undocumented argument [CodonTable *]
 *
 * Return [UNKN ]  Undocumented return value [RandomCodon *]
 *
 */
RandomCodon * Wise2_flat_RandomCodon(CodonTable * ct);
#define flat_RandomCodon Wise2_flat_RandomCodon


/* Function:  RandomModelDNAScore_from_RandomModelDNA(rmd)
 *
 * Descrip:    Gives you a log space RandomModelDNAScore
 *             from a probability space one
 *
 *
 * Arg:        rmd [UNKN ] Undocumented argument [RandomModelDNA *]
 *
 * Return [UNKN ]  Undocumented return value [RandomModelDNAScore *]
 *
 */
RandomModelDNAScore * Wise2_RandomModelDNAScore_from_RandomModelDNA(RandomModelDNA * rmd);
#define RandomModelDNAScore_from_RandomModelDNA Wise2_RandomModelDNAScore_from_RandomModelDNA


/* Function:  RandomModelDNA_std(void)
 *
 * Descrip:    Returns a structure with 0.25 set in each place
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [RandomModelDNA *]
 *
 */
RandomModelDNA * Wise2_RandomModelDNA_std(void);
#define RandomModelDNA_std Wise2_RandomModelDNA_std


/* Function:  RandomModelDNA_std_human(void)
 *
 * Descrip:    Set human random model (slightly G/C)
 *
 *             Not sure where I got the numbers now. Ooops
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [RandomModelDNA *]
 *
 */
RandomModelDNA * Wise2_RandomModelDNA_std_human(void);
#define RandomModelDNA_std_human Wise2_RandomModelDNA_std_human


/* Function:  Score_Sequence_is_random(s,rms)
 *
 * Descrip:    Gives the score of a Sequence vs a random model
 *
 *
 * Arg:          s [UNKN ] Undocumented argument [Sequence *]
 * Arg:        rms [UNKN ] Undocumented argument [RandomModelScoreaa *]
 *
 * Return [UNKN ]  Undocumented return value [Score]
 *
 */
Score Wise2_Score_Sequence_is_random(Sequence * s,RandomModelScoreaa * rms);
#define Score_Sequence_is_random Wise2_Score_Sequence_is_random


/* Function:  Prob_Sequence_is_random(s,rm)
 *
 * Descrip:    Gives the probability of a Sequence vs a random model
 *
 *
 * Arg:         s [UNKN ] Undocumented argument [Sequence *]
 * Arg:        rm [UNKN ] Undocumented argument [RandomModel *]
 *
 * Return [UNKN ]  Undocumented return value [Probability]
 *
 */
Probability Wise2_Prob_Sequence_is_random(Sequence * s,RandomModel * rm);
#define Prob_Sequence_is_random Wise2_Prob_Sequence_is_random


/* Function:  RandomModelScoreaa_from_RandomModel(rm)
 *
 * Descrip:    Gives a score based RandomModel from a probability based one
 *
 *
 * Arg:        rm [UNKN ] Undocumented argument [RandomModel *]
 *
 * Return [UNKN ]  Undocumented return value [RandomModelScoreaa *]
 *
 */
RandomModelScoreaa * Wise2_RandomModelScoreaa_from_RandomModel(RandomModel * rm);
#define RandomModelScoreaa_from_RandomModel Wise2_RandomModelScoreaa_from_RandomModel


/* Function:  default_RandomModel(void)
 *
 * Descrip:    Gives a default random model numbers from
 *             swissprot34- via the HMMEr1 package
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [RandomModel *]
 *
 */
RandomModel * Wise2_default_RandomModel(void);
#define default_RandomModel Wise2_default_RandomModel


/* Function:  read_RandomModel(ifp)
 *
 * Descrip:    Reads a simplistic RandomModel file of
 *
 *             C 0.0123
 *
 *             etc type of format
 *
 *
 * Arg:        ifp [UNKN ] Undocumented argument [FILE *]
 *
 * Return [UNKN ]  Undocumented return value [RandomModel *]
 *
 */
RandomModel * Wise2_read_RandomModel(FILE * ifp);
#define read_RandomModel Wise2_read_RandomModel


/* Function:  hard_link_RandomModel(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [RandomModel *]
 *
 * Return [UNKN ]  Undocumented return value [RandomModel *]
 *
 */
RandomModel * Wise2_hard_link_RandomModel(RandomModel * obj);
#define hard_link_RandomModel Wise2_hard_link_RandomModel


/* Function:  RandomModel_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [RandomModel *]
 *
 */
RandomModel * Wise2_RandomModel_alloc(void);
#define RandomModel_alloc Wise2_RandomModel_alloc


/* Function:  free_RandomModel(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [RandomModel *]
 *
 * Return [UNKN ]  Undocumented return value [RandomModel *]
 *
 */
RandomModel * Wise2_free_RandomModel(RandomModel * obj);
#define free_RandomModel Wise2_free_RandomModel


/* Function:  hard_link_RandomModelScoreaa(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [RandomModelScoreaa *]
 *
 * Return [UNKN ]  Undocumented return value [RandomModelScoreaa *]
 *
 */
RandomModelScoreaa * Wise2_hard_link_RandomModelScoreaa(RandomModelScoreaa * obj);
#define hard_link_RandomModelScoreaa Wise2_hard_link_RandomModelScoreaa


/* Function:  RandomModelScoreaa_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [RandomModelScoreaa *]
 *
 */
RandomModelScoreaa * Wise2_RandomModelScoreaa_alloc(void);
#define RandomModelScoreaa_alloc Wise2_RandomModelScoreaa_alloc


/* Function:  free_RandomModelScoreaa(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [RandomModelScoreaa *]
 *
 * Return [UNKN ]  Undocumented return value [RandomModelScoreaa *]
 *
 */
RandomModelScoreaa * Wise2_free_RandomModelScoreaa(RandomModelScoreaa * obj);
#define free_RandomModelScoreaa Wise2_free_RandomModelScoreaa


/* Function:  hard_link_RandomCodonScore(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [RandomCodonScore *]
 *
 * Return [UNKN ]  Undocumented return value [RandomCodonScore *]
 *
 */
RandomCodonScore * Wise2_hard_link_RandomCodonScore(RandomCodonScore * obj);
#define hard_link_RandomCodonScore Wise2_hard_link_RandomCodonScore


/* Function:  RandomCodonScore_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [RandomCodonScore *]
 *
 */
RandomCodonScore * Wise2_RandomCodonScore_alloc(void);
#define RandomCodonScore_alloc Wise2_RandomCodonScore_alloc


/* Function:  free_RandomCodonScore(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [RandomCodonScore *]
 *
 * Return [UNKN ]  Undocumented return value [RandomCodonScore *]
 *
 */
RandomCodonScore * Wise2_free_RandomCodonScore(RandomCodonScore * obj);
#define free_RandomCodonScore Wise2_free_RandomCodonScore


/* Function:  hard_link_RandomCodon(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [RandomCodon *]
 *
 * Return [UNKN ]  Undocumented return value [RandomCodon *]
 *
 */
RandomCodon * Wise2_hard_link_RandomCodon(RandomCodon * obj);
#define hard_link_RandomCodon Wise2_hard_link_RandomCodon


/* Function:  RandomCodon_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [RandomCodon *]
 *
 */
RandomCodon * Wise2_RandomCodon_alloc(void);
#define RandomCodon_alloc Wise2_RandomCodon_alloc


/* Function:  free_RandomCodon(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [RandomCodon *]
 *
 * Return [UNKN ]  Undocumented return value [RandomCodon *]
 *
 */
RandomCodon * Wise2_free_RandomCodon(RandomCodon * obj);
#define free_RandomCodon Wise2_free_RandomCodon


/* Function:  hard_link_RandomModelDNA(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [RandomModelDNA *]
 *
 * Return [UNKN ]  Undocumented return value [RandomModelDNA *]
 *
 */
RandomModelDNA * Wise2_hard_link_RandomModelDNA(RandomModelDNA * obj);
#define hard_link_RandomModelDNA Wise2_hard_link_RandomModelDNA


/* Function:  RandomModelDNA_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [RandomModelDNA *]
 *
 */
RandomModelDNA * Wise2_RandomModelDNA_alloc(void);
#define RandomModelDNA_alloc Wise2_RandomModelDNA_alloc


/* Function:  free_RandomModelDNA(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [RandomModelDNA *]
 *
 * Return [UNKN ]  Undocumented return value [RandomModelDNA *]
 *
 */
RandomModelDNA * Wise2_free_RandomModelDNA(RandomModelDNA * obj);
#define free_RandomModelDNA Wise2_free_RandomModelDNA


/* Function:  hard_link_RandomModelDNAScore(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [RandomModelDNAScore *]
 *
 * Return [UNKN ]  Undocumented return value [RandomModelDNAScore *]
 *
 */
RandomModelDNAScore * Wise2_hard_link_RandomModelDNAScore(RandomModelDNAScore * obj);
#define hard_link_RandomModelDNAScore Wise2_hard_link_RandomModelDNAScore


/* Function:  RandomModelDNAScore_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [RandomModelDNAScore *]
 *
 */
RandomModelDNAScore * Wise2_RandomModelDNAScore_alloc(void);
#define RandomModelDNAScore_alloc Wise2_RandomModelDNAScore_alloc


/* Function:  free_RandomModelDNAScore(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [RandomModelDNAScore *]
 *
 * Return [UNKN ]  Undocumented return value [RandomModelDNAScore *]
 *
 */
RandomModelDNAScore * Wise2_free_RandomModelDNAScore(RandomModelDNAScore * obj);
#define free_RandomModelDNAScore Wise2_free_RandomModelDNAScore


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
char * Wise2_access_name_RandomModelDNA(RandomModelDNA * obj);
#define access_name_RandomModelDNA Wise2_access_name_RandomModelDNA
boolean Wise2_replace_name_RandomModel(RandomModel * obj,char * name);
#define replace_name_RandomModel Wise2_replace_name_RandomModel
boolean Wise2_replace_name_RandomModelDNA(RandomModelDNA * obj,char * name);
#define replace_name_RandomModelDNA Wise2_replace_name_RandomModelDNA
char * Wise2_access_name_RandomModel(RandomModel * obj);
#define access_name_RandomModel Wise2_access_name_RandomModel

#ifdef _cplusplus
}
#endif

#endif
