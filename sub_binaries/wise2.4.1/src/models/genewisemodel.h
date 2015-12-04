#ifndef DYNAMITEgenewisemodelHEADERFILE
#define DYNAMITEgenewisemodelHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "geneparser21.h"
#include "geneparameter.h"
#include "threestatemodel.h"
#include "codonmapper.h"
#include "cdparser.h"
#include "genefrequency.h"
#include "geneutil.h"

#define GeneWiseScoreLISTLENGTH 128
#define GeneWiseLISTLENGTH 128
#define MAX_PROTEIN_GENEWISE 4096

enum GeneWiseTransition {
  GW_MATCH2MATCH,
  GW_MATCH2INSERT,
  GW_MATCH2DELETE,
  GW_MATCH2END,
  GW_INSERT2MATCH,
  GW_INSERT2INSERT,
  GW_INSERT2DELETE,
  GW_INSERT2END,
  GW_DELETE2MATCH,
  GW_DELETE2INSERT,
  GW_DELETE2DELETE,
  GW_DELETE2END,
  GW_START2MATCH,
  GW_START2INSERT,
  GW_START2DELETE,
  GW_MATCH_BALANCE_5SS,
  GW_INSERT_BALANCE_5SS,
  GW_MATCH_BALANCE_3SS,
  GW_INSERT_BALANCE_3SS,
  GW_TRANSITION_LEN
};

#define GW_EMISSION_LEN 126




/* Object GeneWiseSegment
 *
 * Descrip: This is a particular HMM node, with
 *        match and insert emissions in the codon space
 *        and the transitions
 *
 *        intron/frameshifting transitions are stored
 *        in a different datastructure, as they are
 *        not position dependent
 *
 *
 */
struct Wise2_GeneWiseSegment {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    Probability match[GW_EMISSION_LEN];  
    Probability insert[GW_EMISSION_LEN];     
    Probability transition[GW_TRANSITION_LEN];   
    } ;  
/* GeneWiseSegment defined */ 
#ifndef DYNAMITE_DEFINED_GeneWiseSegment
typedef struct Wise2_GeneWiseSegment Wise2_GeneWiseSegment;
#define GeneWiseSegment Wise2_GeneWiseSegment
#define DYNAMITE_DEFINED_GeneWiseSegment
#endif


/* Object GeneWise
 *
 * Descrip: This is an expand HMM for codon
 *        matching, suitable for genewise and
 *        estwise type algorithms. It is simple
 *        a list of nodes
 *
 *
 */
struct Wise2_GeneWise {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    GeneWiseSegment ** seg;  
    int len;/* len for above seg  */ 
    int maxlen; /* maxlen for above seg */ 
    char * name;     
    } ;  
/* GeneWise defined */ 
#ifndef DYNAMITE_DEFINED_GeneWise
typedef struct Wise2_GeneWise Wise2_GeneWise;
#define GeneWise Wise2_GeneWise
#define DYNAMITE_DEFINED_GeneWise
#endif


/* Object GeneWiseScoreSegment
 *
 * Descrip: This is the log space equivalent
 *        of GeneWiseSegment
 *
 *
 */
struct Wise2_GeneWiseScoreSegment {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    Score match[GW_EMISSION_LEN];    
    Score insert[GW_EMISSION_LEN];   
    Score transition[GW_TRANSITION_LEN];     
    } ;  
/* GeneWiseScoreSegment defined */ 
#ifndef DYNAMITE_DEFINED_GeneWiseScoreSegment
typedef struct Wise2_GeneWiseScoreSegment Wise2_GeneWiseScoreSegment;
#define GeneWiseScoreSegment Wise2_GeneWiseScoreSegment
#define DYNAMITE_DEFINED_GeneWiseScoreSegment
#endif


/* Object GeneWiseScore
 *
 * Descrip: This is the log space equivalent
 *        of the GeneWise 
 *
 *
 */
struct Wise2_GeneWiseScore {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    GeneWiseScoreSegment ** seg;     
    int len;/* len for above seg  */ 
    int maxlen; /* maxlen for above seg */ 
    char * name;     
    } ;  
/* GeneWiseScore defined */ 
#ifndef DYNAMITE_DEFINED_GeneWiseScore
typedef struct Wise2_GeneWiseScore Wise2_GeneWiseScore;
#define GeneWiseScore Wise2_GeneWiseScore
#define DYNAMITE_DEFINED_GeneWiseScore
#endif


/* Object GeneWiseScoreFlat
 *
 * Descrip: This is a specialised datastructure
 *        which is equivalent to the GeneWiseScore
 *        object, but layed out more efficiently
 *        for memory lookup. The actual code is
 *        usually 10% faster. If you have a really
 *        large model however it might barf!
 *
 *
 */
struct Wise2_GeneWiseScoreFlat {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    GeneWiseScoreSegment * seg;  
    int len;     
    } ;  
/* GeneWiseScoreFlat defined */ 
#ifndef DYNAMITE_DEFINED_GeneWiseScoreFlat
typedef struct Wise2_GeneWiseScoreFlat Wise2_GeneWiseScoreFlat;
#define GeneWiseScoreFlat Wise2_GeneWiseScoreFlat
#define DYNAMITE_DEFINED_GeneWiseScoreFlat
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  pack_GeneWiseScore(gws)
 *
 * Descrip:    Packing up the GeneWise model into a byte structure
 *
 *
 * Arg:        gws [UNKN ] Undocumented argument [GeneWiseScore *]
 *
 * Return [UNKN ]  Undocumented return value [char *]
 *
 */
char * Wise2_pack_GeneWiseScore(GeneWiseScore * gws);
#define pack_GeneWiseScore Wise2_pack_GeneWiseScore


/* Function:  GeneWiseScoreFlat_from_GeneWiseScore(gws)
 *
 * Descrip:    This produces a flattened GeneWiseSegment structure
 *             for use in quick implementations (memory lookup
 *             is much better due to everything being a single
 *             piece of memory).
 *
 *
 * Arg:        gws [UNKN ] Undocumented argument [GeneWiseScore *]
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseScoreFlat *]
 *
 */
GeneWiseScoreFlat * Wise2_GeneWiseScoreFlat_from_GeneWiseScore(GeneWiseScore * gws);
#define GeneWiseScoreFlat_from_GeneWiseScore Wise2_GeneWiseScoreFlat_from_GeneWiseScore


/* Function:  free_GeneWiseScoreFlat(obj)
 *
 * Descrip:    Frees the GeneWiseScoreFlat datastructure
 *
 *             overrides the usual deconstructor
 *
 *
 * Arg:        obj [UNKN ] Undocumented argument [GeneWiseScoreFlat *]
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseScoreFlat *]
 *
 */
GeneWiseScoreFlat * Wise2_free_GeneWiseScoreFlat(GeneWiseScoreFlat * obj);
#define free_GeneWiseScoreFlat Wise2_free_GeneWiseScoreFlat


/* Function:  map_phase0_codons_AlnBlock_GeneWise(alb,gws,cseq)
 *
 * Descrip:    This function does something very
 *             sinister.
 *
 *             It maps the phase 0 introns which 
 *             have three base pairs added on the
 *             end. 
 *
 *             It actually changes the AlnBlock structure
 *
 *
 * Arg:         alb [UNKN ] Undocumented argument [AlnBlock *]
 * Arg:         gws [UNKN ] Undocumented argument [GeneWiseScore *]
 * Arg:        cseq [UNKN ] Undocumented argument [ComplexSequence *]
 *
 */
void Wise2_map_phase0_codons_AlnBlock_GeneWise(AlnBlock * alb,GeneWiseScore * gws,ComplexSequence * cseq);
#define map_phase0_codons_AlnBlock_GeneWise Wise2_map_phase0_codons_AlnBlock_GeneWise


/* Function:  flatten_balance_scores_GeneWise(gw)
 *
 * Descrip:    This function is make all the balance scores 0 (hence prob-ratio to 1).
 *
 *             Used when you are using naive models
 *
 *
 * Arg:        gw [UNKN ] genewise model to flatten [GeneWise *]
 *
 */
void Wise2_flatten_balance_scores_GeneWise(GeneWise * gw);
#define flatten_balance_scores_GeneWise Wise2_flatten_balance_scores_GeneWise


/* Function:  GeneWise_from_ThreeStateModel_cdna(tsm,cp,cm,allN)
 *
 * Descrip:    This function makes a 
 *             GeneWise model for the estwise
 *             type algorithms
 *
 *
 * Arg:         tsm [UNKN ] Undocumented argument [ThreeStateModel *]
 * Arg:          cp [UNKN ] Undocumented argument [cDNAParser *]
 * Arg:          cm [UNKN ] Undocumented argument [CodonMapper *]
 * Arg:        allN [UNKN ] Undocumented argument [Probability]
 *
 * Return [UNKN ]  Undocumented return value [GeneWise *]
 *
 */
GeneWise * Wise2_GeneWise_from_ThreeStateModel_cdna(ThreeStateModel * tsm,cDNAParser * cp,CodonMapper * cm,Probability allN);
#define GeneWise_from_ThreeStateModel_cdna Wise2_GeneWise_from_ThreeStateModel_cdna


/* Function:  GeneWise_from_ThreeStateModel_setfactor(tsm,factor,cm,allN)
 *
 * Descrip:    This function makes a 
 *             GeneWise model for the estwise
 *             type algorithms
 *
 *
 * Arg:           tsm [UNKN ] Undocumented argument [ThreeStateModel *]
 * Arg:        factor [UNKN ] Undocumented argument [Probability]
 * Arg:            cm [UNKN ] Undocumented argument [CodonMapper *]
 * Arg:          allN [UNKN ] Undocumented argument [Probability]
 *
 * Return [UNKN ]  Undocumented return value [GeneWise *]
 *
 */
GeneWise * Wise2_GeneWise_from_ThreeStateModel_setfactor(ThreeStateModel * tsm,Probability factor,CodonMapper * cm,Probability allN);
#define GeneWise_from_ThreeStateModel_setfactor Wise2_GeneWise_from_ThreeStateModel_setfactor


/* Function:  GeneWise_from_ThreeStateModel(tsm,gp,cm,allN,gwcm)
 *
 * Descrip:    This makes a genewise model from a 
 *             threestatemodel for the genewise type
 *             algorithms.
 *
 *             Notice you have to provide the gene parameters
 *             being used
 *
 *             Stop is now not used
 *
 *
 * Arg:         tsm [UNKN ] Undocumented argument [ThreeStateModel *]
 * Arg:          gp [UNKN ] Undocumented argument [GeneParser21 *]
 * Arg:          cm [UNKN ] Undocumented argument [CodonMapper *]
 * Arg:        allN [UNKN ] Undocumented argument [Probability]
 * Arg:        gwcm [UNKN ] Undocumented argument [GeneWiseCodonModel *]
 *
 * Return [UNKN ]  Undocumented return value [GeneWise *]
 *
 */
GeneWise * Wise2_GeneWise_from_ThreeStateModel(ThreeStateModel * tsm,GeneParser21 * gp,CodonMapper * cm,Probability allN,GeneWiseCodonModel * gwcm);
#define GeneWise_from_ThreeStateModel Wise2_GeneWise_from_ThreeStateModel


/* Function:  GeneWise_fold_in_synchronised_RandomModel(gw,rm,cm,*ct,stop_codon_background)
 *
 * Descrip:    This function places 'log-odd' scores of the
 *             genewise model assumming that the random model
 *             is a protein model with the codon mapper system
 *             added in, *and* that the path of the random model
 *             is synchronous with the query model.
 *
 *             It fudges stop codons with the stop score given
 *             as a probability.
 *
 *             In other words, this should give bits scores as
 *             if it was a protein, even though it is DNA
 *
 *
 * Arg:                           gw [UNKN ] Undocumented argument [GeneWise *]
 * Arg:                           rm [UNKN ] Undocumented argument [RandomModel *]
 * Arg:                           cm [UNKN ] Undocumented argument [CodonMapper *]
 * Arg:                          *ct [UNKN ] Undocumented argument [CodonTable]
 * Arg:        stop_codon_background [UNKN ] Undocumented argument [Probability]
 *
 */
void Wise2_GeneWise_fold_in_synchronised_RandomModel(GeneWise * gw,RandomModel * rm,CodonMapper * cm,CodonTable *ct,Probability stop_codon_background);
#define GeneWise_fold_in_synchronised_RandomModel Wise2_GeneWise_fold_in_synchronised_RandomModel


/* Function:  check_flat_insert(gw,should_force,should_warn,ct)
 *
 * Descrip:    This function checks that the insert model is bang on
 *             zero, forcing it to zero
 *
 *             Potentially it warns for non zeros as well
 *
 *
 * Arg:                  gw [UNKN ] Undocumented argument [GeneWise *]
 * Arg:        should_force [UNKN ] Undocumented argument [boolean]
 * Arg:         should_warn [UNKN ] Undocumented argument [boolean]
 * Arg:                  ct [UNKN ] Undocumented argument [CodonTable *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_check_flat_insert(GeneWise * gw,boolean should_force,boolean should_warn,CodonTable * ct);
#define check_flat_insert Wise2_check_flat_insert


/* Function:  GeneWise_fold_in_RandomModelDNA(gw,rmd)
 *
 * Descrip:    This function folds in a simple random model
 *             (single base position) into a genewise model
 *
 *
 * Arg:         gw [UNKN ] Undocumented argument [GeneWise *]
 * Arg:        rmd [UNKN ] Undocumented argument [RandomModelDNA *]
 *
 */
void Wise2_GeneWise_fold_in_RandomModelDNA(GeneWise * gw,RandomModelDNA * rmd);
#define GeneWise_fold_in_RandomModelDNA Wise2_GeneWise_fold_in_RandomModelDNA


/* Function:  Protein_from_GeneWise_AlnColumn(dna,is_random_AlnColumn,col,position_in_aln,last_column,ct)
 *
 * Descrip:    Produces a protein object from a genewise/estwise
 *             style label set, setting the last retrieved column
 *
 *
 * Arg:                        dna [UNKN ] Undocumented argument [Sequence *]
 * Arg:        is_random_AlnColumn [UNKN ] Undocumented argument [NullString]
 * Arg:                        col [UNKN ] Undocumented argument [AlnColumn *]
 * Arg:            position_in_aln [UNKN ] Undocumented argument [int]
 * Arg:                last_column [UNKN ] Undocumented argument [AlnColumn **]
 * Arg:                         ct [UNKN ] Undocumented argument [CodonTable *]
 *
 * Return [UNKN ]  Undocumented return value [Protein *]
 *
 */
Protein * Wise2_Protein_from_GeneWise_AlnColumn(Sequence * dna,AlnColumn * col,int position_in_aln,AlnColumn ** last_column,CodonTable * ct,boolean (*is_random_AlnColumn)(const AlnColumn *));
#define Protein_from_GeneWise_AlnColumn Wise2_Protein_from_GeneWise_AlnColumn


/* Function:  GeneWiseScore_from_GeneWise(gw)
 *
 * Descrip:    Makes a Score (log based) object from
 *             a probability based object
 *
 *
 * Arg:        gw [UNKN ] Undocumented argument [GeneWise *]
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseScore *]
 *
 */
GeneWiseScore * Wise2_GeneWiseScore_from_GeneWise(GeneWise * gw);
#define GeneWiseScore_from_GeneWise Wise2_GeneWiseScore_from_GeneWise


/* Function:  hard_link_GeneWiseSegment(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GeneWiseSegment *]
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseSegment *]
 *
 */
GeneWiseSegment * Wise2_hard_link_GeneWiseSegment(GeneWiseSegment * obj);
#define hard_link_GeneWiseSegment Wise2_hard_link_GeneWiseSegment


/* Function:  GeneWiseSegment_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseSegment *]
 *
 */
GeneWiseSegment * Wise2_GeneWiseSegment_alloc(void);
#define GeneWiseSegment_alloc Wise2_GeneWiseSegment_alloc


/* Function:  free_GeneWiseSegment(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GeneWiseSegment *]
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseSegment *]
 *
 */
GeneWiseSegment * Wise2_free_GeneWiseSegment(GeneWiseSegment * obj);
#define free_GeneWiseSegment Wise2_free_GeneWiseSegment


/* Function:  add_GeneWise(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [GeneWise *]
 * Arg:        add [OWNER] Object to add to the list [GeneWiseSegment *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_GeneWise(GeneWise * obj,GeneWiseSegment * add);
#define add_GeneWise Wise2_add_GeneWise


/* Function:  flush_GeneWise(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [GeneWise *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_GeneWise(GeneWise * obj);
#define flush_GeneWise Wise2_flush_GeneWise


/* Function:  GeneWise_alloc_std(void)
 *
 * Descrip:    Equivalent to GeneWise_alloc_len(GeneWiseLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GeneWise *]
 *
 */
GeneWise * Wise2_GeneWise_alloc_std(void);
#define GeneWise_alloc_std Wise2_GeneWise_alloc_std


/* Function:  GeneWise_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [GeneWise *]
 *
 */
GeneWise * Wise2_GeneWise_alloc_len(int len);
#define GeneWise_alloc_len Wise2_GeneWise_alloc_len


/* Function:  hard_link_GeneWise(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GeneWise *]
 *
 * Return [UNKN ]  Undocumented return value [GeneWise *]
 *
 */
GeneWise * Wise2_hard_link_GeneWise(GeneWise * obj);
#define hard_link_GeneWise Wise2_hard_link_GeneWise


/* Function:  GeneWise_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GeneWise *]
 *
 */
GeneWise * Wise2_GeneWise_alloc(void);
#define GeneWise_alloc Wise2_GeneWise_alloc


/* Function:  free_GeneWise(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GeneWise *]
 *
 * Return [UNKN ]  Undocumented return value [GeneWise *]
 *
 */
GeneWise * Wise2_free_GeneWise(GeneWise * obj);
#define free_GeneWise Wise2_free_GeneWise


/* Function:  hard_link_GeneWiseScoreSegment(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GeneWiseScoreSegment *]
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseScoreSegment *]
 *
 */
GeneWiseScoreSegment * Wise2_hard_link_GeneWiseScoreSegment(GeneWiseScoreSegment * obj);
#define hard_link_GeneWiseScoreSegment Wise2_hard_link_GeneWiseScoreSegment


/* Function:  GeneWiseScoreSegment_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseScoreSegment *]
 *
 */
GeneWiseScoreSegment * Wise2_GeneWiseScoreSegment_alloc(void);
#define GeneWiseScoreSegment_alloc Wise2_GeneWiseScoreSegment_alloc


/* Function:  free_GeneWiseScoreSegment(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GeneWiseScoreSegment *]
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseScoreSegment *]
 *
 */
GeneWiseScoreSegment * Wise2_free_GeneWiseScoreSegment(GeneWiseScoreSegment * obj);
#define free_GeneWiseScoreSegment Wise2_free_GeneWiseScoreSegment


/* Function:  add_GeneWiseScore(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [GeneWiseScore *]
 * Arg:        add [OWNER] Object to add to the list [GeneWiseScoreSegment *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_GeneWiseScore(GeneWiseScore * obj,GeneWiseScoreSegment * add);
#define add_GeneWiseScore Wise2_add_GeneWiseScore


/* Function:  flush_GeneWiseScore(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [GeneWiseScore *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_GeneWiseScore(GeneWiseScore * obj);
#define flush_GeneWiseScore Wise2_flush_GeneWiseScore


/* Function:  GeneWiseScore_alloc_std(void)
 *
 * Descrip:    Equivalent to GeneWiseScore_alloc_len(GeneWiseScoreLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseScore *]
 *
 */
GeneWiseScore * Wise2_GeneWiseScore_alloc_std(void);
#define GeneWiseScore_alloc_std Wise2_GeneWiseScore_alloc_std


/* Function:  GeneWiseScore_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseScore *]
 *
 */
GeneWiseScore * Wise2_GeneWiseScore_alloc_len(int len);
#define GeneWiseScore_alloc_len Wise2_GeneWiseScore_alloc_len


/* Function:  hard_link_GeneWiseScore(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GeneWiseScore *]
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseScore *]
 *
 */
GeneWiseScore * Wise2_hard_link_GeneWiseScore(GeneWiseScore * obj);
#define hard_link_GeneWiseScore Wise2_hard_link_GeneWiseScore


/* Function:  GeneWiseScore_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseScore *]
 *
 */
GeneWiseScore * Wise2_GeneWiseScore_alloc(void);
#define GeneWiseScore_alloc Wise2_GeneWiseScore_alloc


/* Function:  free_GeneWiseScore(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GeneWiseScore *]
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseScore *]
 *
 */
GeneWiseScore * Wise2_free_GeneWiseScore(GeneWiseScore * obj);
#define free_GeneWiseScore Wise2_free_GeneWiseScore


/* Function:  hard_link_GeneWiseScoreFlat(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GeneWiseScoreFlat *]
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseScoreFlat *]
 *
 */
GeneWiseScoreFlat * Wise2_hard_link_GeneWiseScoreFlat(GeneWiseScoreFlat * obj);
#define hard_link_GeneWiseScoreFlat Wise2_hard_link_GeneWiseScoreFlat


/* Function:  GeneWiseScoreFlat_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseScoreFlat *]
 *
 */
GeneWiseScoreFlat * Wise2_GeneWiseScoreFlat_alloc(void);
#define GeneWiseScoreFlat_alloc Wise2_GeneWiseScoreFlat_alloc


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
Probability Wise2_probability_of_this_codon(int codon,RandomModelDNA * rmd);
#define probability_of_this_codon Wise2_probability_of_this_codon
void Wise2_show_GeneWise(GeneWise * gw,FILE * ofp);
#define show_GeneWise Wise2_show_GeneWise
void Wise2_show_GeneWiseSegment(GeneWiseSegment * seg,FILE * ofp);
#define show_GeneWiseSegment Wise2_show_GeneWiseSegment
GeneWiseSegment * Wise2_GeneWiseSegment_from_ThreeStateUnit(ThreeStateUnit * tsu,Probability factor,CodonMapper * cm,GeneWiseCodonModel * gwcm,Probability allN);
#define GeneWiseSegment_from_ThreeStateUnit Wise2_GeneWiseSegment_from_ThreeStateUnit
Probability Wise2_Probability_of_codon(int codon,CodonTable * ct,Probability * aminoacid_26_array,Probability stop);
#define Probability_of_codon Wise2_Probability_of_codon
GeneWiseScoreSegment * Wise2_GeneWiseScoreSegment_from_GeneWiseSegment(GeneWiseSegment * prev,GeneWiseSegment * seg);
#define GeneWiseScoreSegment_from_GeneWiseSegment Wise2_GeneWiseScoreSegment_from_GeneWiseSegment
void Wise2_swap_GeneWise(GeneWiseSegment ** list,int i,int j) ;
#define swap_GeneWise Wise2_swap_GeneWise
void Wise2_qsort_GeneWise(GeneWiseSegment ** list,int left,int right,int (*comp)(GeneWiseSegment * ,GeneWiseSegment * ));
#define qsort_GeneWise Wise2_qsort_GeneWise
void Wise2_sort_GeneWise(GeneWise * obj,int (*comp)(GeneWiseSegment *, GeneWiseSegment *));
#define sort_GeneWise Wise2_sort_GeneWise
boolean Wise2_expand_GeneWise(GeneWise * obj,int len);
#define expand_GeneWise Wise2_expand_GeneWise
void Wise2_swap_GeneWiseScore(GeneWiseScoreSegment ** list,int i,int j) ;
#define swap_GeneWiseScore Wise2_swap_GeneWiseScore
void Wise2_qsort_GeneWiseScore(GeneWiseScoreSegment ** list,int left,int right,int (*comp)(GeneWiseScoreSegment * ,GeneWiseScoreSegment * ));
#define qsort_GeneWiseScore Wise2_qsort_GeneWiseScore
void Wise2_sort_GeneWiseScore(GeneWiseScore * obj,int (*comp)(GeneWiseScoreSegment *, GeneWiseScoreSegment *));
#define sort_GeneWiseScore Wise2_sort_GeneWiseScore
boolean Wise2_expand_GeneWiseScore(GeneWiseScore * obj,int len);
#define expand_GeneWiseScore Wise2_expand_GeneWiseScore

#ifdef _cplusplus
}
#endif

#endif
