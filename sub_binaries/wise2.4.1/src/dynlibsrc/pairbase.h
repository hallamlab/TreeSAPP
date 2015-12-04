#ifndef DYNAMITEpairbaseHEADERFILE
#define DYNAMITEpairbaseHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "codon.h"
#include "probability.h"
#include "codonmatrix.h"
#include <stdio.h>

typedef char pairbase_type;
typedef int  pairbase_codon_type;

#define BASE_GAP  5
#define BASE_OPEN 6

#define IS_NOT_BASE(a) (a == BASE_GAP ? 1 : a == BASE_OPEN ? 1 : 0)
#define MAKE_PAIRBASE(anchor,informant) (anchor*7+informant)

#define PAIRBASE_LENGTH (7*7)

#define PAIRBASE_CODON_LENGTH (PAIRBASE_LENGTH*PAIRBASE_LENGTH*PAIRBASE_LENGTH)

struct Wise2_PairBaseModel {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    Probability base[PAIRBASE_LENGTH];   
    } ;  
/* PairBaseModel defined */ 
#ifndef DYNAMITE_DEFINED_PairBaseModel
typedef struct Wise2_PairBaseModel Wise2_PairBaseModel;
#define PairBaseModel Wise2_PairBaseModel
#define DYNAMITE_DEFINED_PairBaseModel
#endif


struct Wise2_PairBaseCodonModel {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    Probability codon[PAIRBASE_CODON_LENGTH];    
    } ;  
/* PairBaseCodonModel defined */ 
#ifndef DYNAMITE_DEFINED_PairBaseCodonModel
typedef struct Wise2_PairBaseCodonModel Wise2_PairBaseCodonModel;
#define PairBaseCodonModel Wise2_PairBaseCodonModel
#define DYNAMITE_DEFINED_PairBaseCodonModel
#endif


struct Wise2_PairBaseModelScore {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    Score base[PAIRBASE_LENGTH];     
    } ;  
/* PairBaseModelScore defined */ 
#ifndef DYNAMITE_DEFINED_PairBaseModelScore
typedef struct Wise2_PairBaseModelScore Wise2_PairBaseModelScore;
#define PairBaseModelScore Wise2_PairBaseModelScore
#define DYNAMITE_DEFINED_PairBaseModelScore
#endif


struct Wise2_PairBaseCodonModelScore {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    Score codon[PAIRBASE_CODON_LENGTH];  
    } ;  
/* PairBaseCodonModelScore defined */ 
#ifndef DYNAMITE_DEFINED_PairBaseCodonModelScore
typedef struct Wise2_PairBaseCodonModelScore Wise2_PairBaseCodonModelScore;
#define PairBaseCodonModelScore Wise2_PairBaseCodonModelScore
#define DYNAMITE_DEFINED_PairBaseCodonModelScore
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  diagonal_tweak_PairBaseCodonModel(m,ratio_on,ratio_off_positive,ratio_off_negative)
 *
 * Descrip:    Tweaks a PairBaseCodonModel with on and off diagonal ratios
 *
 *
 * Arg:                         m [UNKN ] Undocumented argument [PairBaseCodonModel *]
 * Arg:                  ratio_on [UNKN ] Undocumented argument [double]
 * Arg:        ratio_off_positive [UNKN ] Undocumented argument [double]
 * Arg:        ratio_off_negative [UNKN ] Undocumented argument [double]
 *
 */
void Wise2_diagonal_tweak_PairBaseCodonModel(PairBaseCodonModel * m,double ratio_on,double ratio_off_positive,double ratio_off_negative);
#define diagonal_tweak_PairBaseCodonModel Wise2_diagonal_tweak_PairBaseCodonModel


/* Function:  flatten_diagonal_PairBaseCodonModel(m,ct)
 *
 * Descrip:    flattens out the diagonal signal
 *
 *
 * Arg:         m [UNKN ] Undocumented argument [PairBaseCodonModel *]
 * Arg:        ct [UNKN ] Undocumented argument [CodonTable *]
 *
 */
void Wise2_flatten_diagonal_PairBaseCodonModel(PairBaseCodonModel * m,CodonTable * ct);
#define flatten_diagonal_PairBaseCodonModel Wise2_flatten_diagonal_PairBaseCodonModel


/* Function:  flatten_diagonal_PairBaseCodonModelScore(m,ct)
 *
 * Descrip:    flattens out the diagonal signal - for scores!
 *
 *
 * Arg:         m [UNKN ] Undocumented argument [PairBaseCodonModelScore *]
 * Arg:        ct [UNKN ] Undocumented argument [CodonTable *]
 *
 */
void Wise2_flatten_diagonal_PairBaseCodonModelScore(PairBaseCodonModelScore * m,CodonTable * ct);
#define flatten_diagonal_PairBaseCodonModelScore Wise2_flatten_diagonal_PairBaseCodonModelScore


/* Function:  zero_PairBaseModelScore(void)
 *
 * Descrip:    a 0 pairbasemodel score 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [PairBaseModelScore *]
 *
 */
PairBaseModelScore * Wise2_zero_PairBaseModelScore(void);
#define zero_PairBaseModelScore Wise2_zero_PairBaseModelScore


/* Function:  very_simple_PairBaseCodonModel(id,rnd,nonm,gap,ct)
 *
 * Descrip:    Makes a PairBaseCodonModel from just a one parameter! Wow!
 *
 *
 * Arg:          id [UNKN ] Undocumented argument [Probability]
 * Arg:         rnd [UNKN ] Undocumented argument [Probability]
 * Arg:        nonm [UNKN ] Undocumented argument [Probability]
 * Arg:         gap [UNKN ] Undocumented argument [Probability]
 * Arg:          ct [UNKN ] Undocumented argument [CodonTable *]
 *
 * Return [UNKN ]  Undocumented return value [PairBaseCodonModel *]
 *
 */
PairBaseCodonModel * Wise2_very_simple_PairBaseCodonModel(Probability id,Probability rnd,Probability nonm,Probability gap,CodonTable * ct);
#define very_simple_PairBaseCodonModel Wise2_very_simple_PairBaseCodonModel


/* Function:  make_flat_PairBaseCodonModel(cp,nonm,gap,ct)
 *
 * Descrip:    Makes a PairBaseCodonModel from a protein matrix - assumming a flat 
 *             mapping to CodonMatrix
 *
 *
 * Arg:          cp [UNKN ] Undocumented argument [CompProb *]
 * Arg:        nonm [UNKN ] Undocumented argument [Probability]
 * Arg:         gap [UNKN ] Undocumented argument [Probability]
 * Arg:          ct [UNKN ] Undocumented argument [CodonTable *]
 *
 * Return [UNKN ]  Undocumented return value [PairBaseCodonModel *]
 *
 */
PairBaseCodonModel * Wise2_make_flat_PairBaseCodonModel(CompProb * cp,Probability nonm,Probability gap,CodonTable * ct);
#define make_flat_PairBaseCodonModel Wise2_make_flat_PairBaseCodonModel


/* Function:  make_start_PairBaseCodonModelScore(ct)
 *
 * Descrip:    Makes a PairBaseCodonModel score for start codon
 *
 *
 * Arg:        ct [UNKN ] Undocumented argument [CodonTable *]
 *
 * Return [UNKN ]  Undocumented return value [PairBaseCodonModelScore *]
 *
 */
PairBaseCodonModelScore * Wise2_make_start_PairBaseCodonModelScore(CodonTable * ct);
#define make_start_PairBaseCodonModelScore Wise2_make_start_PairBaseCodonModelScore


/* Function:  make_stop_PairBaseCodonModelScore(ct)
 *
 * Descrip:    Makes a PairBaseCodonModel score for start codon
 *
 *
 * Arg:        ct [UNKN ] Undocumented argument [CodonTable *]
 *
 * Return [UNKN ]  Undocumented return value [PairBaseCodonModelScore *]
 *
 */
PairBaseCodonModelScore * Wise2_make_stop_PairBaseCodonModelScore(CodonTable * ct);
#define make_stop_PairBaseCodonModelScore Wise2_make_stop_PairBaseCodonModelScore


/* Function:  make_conserved_PairBaseCodonModel(cons,non_cons,aa_var,ct)
 *
 * Descrip:    Makes a PairBaseCodonModel for a particular character in the CodonTable
 *
 *
 * Arg:            cons [UNKN ] Undocumented argument [Probability]
 * Arg:        non_cons [UNKN ] Undocumented argument [Probability]
 * Arg:          aa_var [UNKN ] Undocumented argument [char]
 * Arg:              ct [UNKN ] Undocumented argument [CodonTable *]
 *
 * Return [UNKN ]  Undocumented return value [PairBaseCodonModel *]
 *
 */
PairBaseCodonModel * Wise2_make_conserved_PairBaseCodonModel(Probability cons,Probability non_cons,char aa_var,CodonTable * ct);
#define make_conserved_PairBaseCodonModel Wise2_make_conserved_PairBaseCodonModel


/* Function:  make_PairBaseCodonModel(codon_matrix,nonm,gap,ct)
 *
 * Descrip:    Makes a PairBaseCodonModel from a CodonMatrix and parameters
 *
 *
 * Arg:        codon_matrix [UNKN ] Undocumented argument [CodonMatrix *]
 * Arg:                nonm [UNKN ] Undocumented argument [Probability]
 * Arg:                 gap [UNKN ] Undocumented argument [Probability]
 * Arg:                  ct [UNKN ] Undocumented argument [CodonTable *]
 *
 * Return [UNKN ]  Undocumented return value [PairBaseCodonModel *]
 *
 */
PairBaseCodonModel * Wise2_make_PairBaseCodonModel(CodonMatrix * codon_matrix,Probability nonm,Probability gap,CodonTable * ct);
#define make_PairBaseCodonModel Wise2_make_PairBaseCodonModel


/* Function:  simple_PairBaseModel(iden,other,gap)
 *
 * Descrip:    Makes a pair base model from simple leading diagonal
 *
 *
 * Arg:         iden [UNKN ] Undocumented argument [Probability]
 * Arg:        other [UNKN ] Undocumented argument [Probability]
 * Arg:          gap [UNKN ] Undocumented argument [Probability]
 *
 * Return [UNKN ]  Undocumented return value [PairBaseModel *]
 *
 */
PairBaseModel * Wise2_simple_PairBaseModel(Probability iden,Probability other,Probability gap);
#define simple_PairBaseModel Wise2_simple_PairBaseModel


/* Function:  new_PairBaseCodonModelScore(pbcm)
 *
 * Descrip:    Makes a codon score from a codon model
 *
 *
 * Arg:        pbcm [UNKN ] Undocumented argument [PairBaseCodonModel *]
 *
 * Return [UNKN ]  Undocumented return value [PairBaseCodonModelScore *]
 *
 */
PairBaseCodonModelScore * Wise2_new_PairBaseCodonModelScore(PairBaseCodonModel * pbcm);
#define new_PairBaseCodonModelScore Wise2_new_PairBaseCodonModelScore


/* Function:  new_PairBaseModelScore(pbm)
 *
 * Descrip:    Makes a base score from a base model
 *
 *
 * Arg:        pbm [UNKN ] Undocumented argument [PairBaseModel *]
 *
 * Return [UNKN ]  Undocumented return value [PairBaseModelScore *]
 *
 */
PairBaseModelScore * Wise2_new_PairBaseModelScore(PairBaseModel * pbm);
#define new_PairBaseModelScore Wise2_new_PairBaseModelScore


/* Function:  show_PairBaseModelScore(sc,ofp)
 *
 * Descrip:    Debugging
 *
 *
 * Arg:         sc [UNKN ] Undocumented argument [PairBaseModelScore *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
void Wise2_show_PairBaseModelScore(PairBaseModelScore * sc,FILE * ofp);
#define show_PairBaseModelScore Wise2_show_PairBaseModelScore


/* Function:  show_PairBaseCodonModelScore(sc,ct,ofp)
 *
 * Descrip:    Debugging
 *
 *
 * Arg:         sc [UNKN ] Undocumented argument [PairBaseCodonModelScore *]
 * Arg:         ct [UNKN ] Undocumented argument [CodonTable *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
void Wise2_show_PairBaseCodonModelScore(PairBaseCodonModelScore * sc,CodonTable * ct,FILE * ofp);
#define show_PairBaseCodonModelScore Wise2_show_PairBaseCodonModelScore


/* Function:  reverse_pairbase_codon(codon)
 *
 * Descrip:    Inverts pairbase codon
 *
 *
 * Arg:        codon [UNKN ] Undocumented argument [pairbase_codon_type]
 *
 * Return [UNKN ]  Undocumented return value [pairbase_codon_type]
 *
 */
pairbase_codon_type Wise2_reverse_pairbase_codon(pairbase_codon_type codon);
#define reverse_pairbase_codon Wise2_reverse_pairbase_codon


/* Function:  complement_pairbase(b)
 *
 * Descrip:    complements a pairbase 
 *
 *
 * Arg:        b [UNKN ] Undocumented argument [pairbase_type]
 *
 * Return [UNKN ]  Undocumented return value [pairbase_type]
 *
 */
pairbase_type Wise2_complement_pairbase(pairbase_type b);
#define complement_pairbase Wise2_complement_pairbase


/* Function:  pairbase_codon_from_seq(seq)
 *
 * Descrip:    Makes a pairbase_codon from a pairbase_sequence
 *
 *
 * Arg:        seq [UNKN ] Undocumented argument [pairbase_type *]
 *
 * Return [UNKN ]  Undocumented return value [pairbase_codon_type]
 *
 */
pairbase_codon_type Wise2_pairbase_codon_from_seq(pairbase_type * seq);
#define pairbase_codon_from_seq Wise2_pairbase_codon_from_seq


/* Function:  decompose_pairbase_codon(t,a,b,c)
 *
 * Descrip:    Decomposes a pairbase codon
 *
 *
 * Arg:        t [UNKN ] Undocumented argument [pairbase_codon_type]
 * Arg:        a [UNKN ] Undocumented argument [pairbase_type *]
 * Arg:        b [UNKN ] Undocumented argument [pairbase_type *]
 * Arg:        c [UNKN ] Undocumented argument [pairbase_type *]
 *
 */
void Wise2_decompose_pairbase_codon(pairbase_codon_type t,pairbase_type * a,pairbase_type * b,pairbase_type * c);
#define decompose_pairbase_codon Wise2_decompose_pairbase_codon


/* Function:  anchor_base_from_pairbase(pairbase)
 *
 * Descrip:    Finds the anchor base from a pair base
 *
 *
 * Arg:        pairbase [UNKN ] Undocumented argument [pairbase_type]
 *
 * Return [UNKN ]  Undocumented return value [base]
 *
 */
base Wise2_anchor_base_from_pairbase(pairbase_type pairbase);
#define anchor_base_from_pairbase Wise2_anchor_base_from_pairbase


/* Function:  informant_base_from_pairbase(pairbase)
 *
 * Descrip:    Finds the informant base from a pair base
 *
 *
 * Arg:        pairbase [UNKN ] Undocumented argument [pairbase_type]
 *
 * Return [UNKN ]  Undocumented return value [base]
 *
 */
base Wise2_informant_base_from_pairbase(pairbase_type pairbase);
#define informant_base_from_pairbase Wise2_informant_base_from_pairbase


/* Function:  char_for_base(base)
 *
 * Descrip:    gives back the character for the base
 *
 *
 * Arg:        base [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [char]
 *
 */
char Wise2_char_for_base(int base);
#define char_for_base Wise2_char_for_base


/* Function:  hard_link_PairBaseModel(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [PairBaseModel *]
 *
 * Return [UNKN ]  Undocumented return value [PairBaseModel *]
 *
 */
PairBaseModel * Wise2_hard_link_PairBaseModel(PairBaseModel * obj);
#define hard_link_PairBaseModel Wise2_hard_link_PairBaseModel


/* Function:  PairBaseModel_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [PairBaseModel *]
 *
 */
PairBaseModel * Wise2_PairBaseModel_alloc(void);
#define PairBaseModel_alloc Wise2_PairBaseModel_alloc


/* Function:  free_PairBaseModel(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [PairBaseModel *]
 *
 * Return [UNKN ]  Undocumented return value [PairBaseModel *]
 *
 */
PairBaseModel * Wise2_free_PairBaseModel(PairBaseModel * obj);
#define free_PairBaseModel Wise2_free_PairBaseModel


/* Function:  hard_link_PairBaseCodonModel(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [PairBaseCodonModel *]
 *
 * Return [UNKN ]  Undocumented return value [PairBaseCodonModel *]
 *
 */
PairBaseCodonModel * Wise2_hard_link_PairBaseCodonModel(PairBaseCodonModel * obj);
#define hard_link_PairBaseCodonModel Wise2_hard_link_PairBaseCodonModel


/* Function:  PairBaseCodonModel_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [PairBaseCodonModel *]
 *
 */
PairBaseCodonModel * Wise2_PairBaseCodonModel_alloc(void);
#define PairBaseCodonModel_alloc Wise2_PairBaseCodonModel_alloc


/* Function:  free_PairBaseCodonModel(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [PairBaseCodonModel *]
 *
 * Return [UNKN ]  Undocumented return value [PairBaseCodonModel *]
 *
 */
PairBaseCodonModel * Wise2_free_PairBaseCodonModel(PairBaseCodonModel * obj);
#define free_PairBaseCodonModel Wise2_free_PairBaseCodonModel


/* Function:  hard_link_PairBaseModelScore(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [PairBaseModelScore *]
 *
 * Return [UNKN ]  Undocumented return value [PairBaseModelScore *]
 *
 */
PairBaseModelScore * Wise2_hard_link_PairBaseModelScore(PairBaseModelScore * obj);
#define hard_link_PairBaseModelScore Wise2_hard_link_PairBaseModelScore


/* Function:  PairBaseModelScore_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [PairBaseModelScore *]
 *
 */
PairBaseModelScore * Wise2_PairBaseModelScore_alloc(void);
#define PairBaseModelScore_alloc Wise2_PairBaseModelScore_alloc


/* Function:  free_PairBaseModelScore(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [PairBaseModelScore *]
 *
 * Return [UNKN ]  Undocumented return value [PairBaseModelScore *]
 *
 */
PairBaseModelScore * Wise2_free_PairBaseModelScore(PairBaseModelScore * obj);
#define free_PairBaseModelScore Wise2_free_PairBaseModelScore


/* Function:  hard_link_PairBaseCodonModelScore(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [PairBaseCodonModelScore *]
 *
 * Return [UNKN ]  Undocumented return value [PairBaseCodonModelScore *]
 *
 */
PairBaseCodonModelScore * Wise2_hard_link_PairBaseCodonModelScore(PairBaseCodonModelScore * obj);
#define hard_link_PairBaseCodonModelScore Wise2_hard_link_PairBaseCodonModelScore


/* Function:  PairBaseCodonModelScore_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [PairBaseCodonModelScore *]
 *
 */
PairBaseCodonModelScore * Wise2_PairBaseCodonModelScore_alloc(void);
#define PairBaseCodonModelScore_alloc Wise2_PairBaseCodonModelScore_alloc


/* Function:  free_PairBaseCodonModelScore(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [PairBaseCodonModelScore *]
 *
 * Return [UNKN ]  Undocumented return value [PairBaseCodonModelScore *]
 *
 */
PairBaseCodonModelScore * Wise2_free_PairBaseCodonModelScore(PairBaseCodonModelScore * obj);
#define free_PairBaseCodonModelScore Wise2_free_PairBaseCodonModelScore


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
