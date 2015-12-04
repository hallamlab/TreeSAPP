#ifndef DYNAMITEthreestatemodelHEADERFILE
#define DYNAMITEthreestatemodelHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif

#include "dyna.h"
#include "randommodel.h"

#define ThreeStateModelLISTLENGTH 128
#define ThreeStateScoreLISTLENGTH 128

enum tsm_trans_type {
  TSM_MATCH2MATCH  = 0,
  TSM_MATCH2INSERT,  
  TSM_MATCH2DELETE,
  TSM_INSERT2MATCH,
  TSM_INSERT2INSERT,  
  TSM_INSERT2DELETE,
  TSM_DELETE2MATCH,
  TSM_DELETE2INSERT,  
  TSM_DELETE2DELETE,  
  TSM_START2MATCH,    
  TSM_START2INSERT,   
  TSM_START2DELETE,   
  TSM_MATCH2END,      
  TSM_INSERT2END,     
  TSM_DELETE2END };

#define TRANSITION_LEN 15

#define BLANK_SCORE  (-10000000)

extern char * std_alphabet;
#define ALPHABET_SIZE 26

typedef enum TSM_StartEndMode {
  TSM_unknown = 26,
  TSM_default,
  TSM_local,
  TSM_global,
  TSM_wing,
  TSM_endbiased
} TSM_StartEndMode;


/* Object ThreeStateUnit
 *
 * Descrip: This object is the probability version
 *        of hte common unit to profile HMMs, ie
 *        the match,insert,delete triple
 *
 *
 */
struct Wise2_ThreeStateUnit {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    Probability match_emission[ALPHABET_SIZE];   
    Probability insert_emission[ALPHABET_SIZE];  
    Probability transition[TRANSITION_LEN];  
    char display_char;   
    } ;  
/* ThreeStateUnit defined */ 
#ifndef DYNAMITE_DEFINED_ThreeStateUnit
typedef struct Wise2_ThreeStateUnit Wise2_ThreeStateUnit;
#define ThreeStateUnit Wise2_ThreeStateUnit
#define DYNAMITE_DEFINED_ThreeStateUnit
#endif


/* Object ThreeStateModel
 *
 * Descrip: This is profile-HMM object, similar to the
 *        SAM and HMMer plan9 architecture. 
 *
 *
 */
struct Wise2_ThreeStateModel {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    char * name;    /*  name of the model */ 
    ThreeStateUnit ** unit; /*  the actuall three state probs and emissions */ 
    int len;/* len for above unit  */ 
    int maxlen; /* maxlen for above unit */ 
    char * alphabet;    /*  alphabet used  */ 
    char * accession;   /*  accession number */ 
    double threshold;   /*  bits threshold (if sensible) */ 
    RandomModel * rm;   /*  Random model for the model: maybe NULL! */ 
    } ;  
/* ThreeStateModel defined */ 
#ifndef DYNAMITE_DEFINED_ThreeStateModel
typedef struct Wise2_ThreeStateModel Wise2_ThreeStateModel;
#define ThreeStateModel Wise2_ThreeStateModel
#define DYNAMITE_DEFINED_ThreeStateModel
#endif


/* Object ThreeStateScoreUnit
 *
 * Descrip: This is the Score (ie , log-odded)
 *        version of the ThreeStateUnit
 *
 *
 */
struct Wise2_ThreeStateScoreUnit {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    Score match[ALPHABET_SIZE];  
    Score insert[ALPHABET_SIZE];     
    Score trans[TRANSITION_LEN];     
    } ;  
/* ThreeStateScoreUnit defined */ 
#ifndef DYNAMITE_DEFINED_ThreeStateScoreUnit
typedef struct Wise2_ThreeStateScoreUnit Wise2_ThreeStateScoreUnit;
#define ThreeStateScoreUnit Wise2_ThreeStateScoreUnit
#define DYNAMITE_DEFINED_ThreeStateScoreUnit
#endif


/* Object ThreeStateScore
 *
 * Descrip: This is hte score version of the
 *        ThreeStateModel, ie this has its
 *        emissions and transitions log-odded
 *
 *
 */
struct Wise2_ThreeStateScore {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    ThreeStateScoreUnit ** unit;     
    int len;/* len for above unit  */ 
    int maxlen; /* maxlen for above unit */ 
    char * name;     
    char * accession;    
    } ;  
/* ThreeStateScore defined */ 
#ifndef DYNAMITE_DEFINED_ThreeStateScore
typedef struct Wise2_ThreeStateScore Wise2_ThreeStateScore;
#define ThreeStateScore Wise2_ThreeStateScore
#define DYNAMITE_DEFINED_ThreeStateScore
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  information_from_ThreeStateUnit(tsu,rm)
 *
 * Descrip:     Gets the information content (K-L divergence) vs a background for a position
 *
 *
 * Arg:        tsu [UNKN ] Undocumented argument [ThreeStateUnit *]
 * Arg:         rm [UNKN ] Undocumented argument [RandomModel *]
 *
 * Return [UNKN ]  Undocumented return value [double]
 *
 */
double Wise2_information_from_ThreeStateUnit(ThreeStateUnit * tsu,RandomModel * rm);
#define information_from_ThreeStateUnit Wise2_information_from_ThreeStateUnit


/* Function:  threestatemodel_mode_from_string(mode)
 *
 * Descrip:    gets out the mode from a string
 *
 *
 * Arg:        mode [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_threestatemodel_mode_from_string(char * mode);
#define threestatemodel_mode_from_string Wise2_threestatemodel_mode_from_string


/* Function:  set_startend_policy_ThreeStateModel(tsm,mode,wing_length,internal_bias)
 *
 * Descrip:    Sets the start/end policy on the basis of the mode
 *
 *
 * Arg:                  tsm [UNKN ] Undocumented argument [ThreeStateModel *]
 * Arg:                 mode [UNKN ] Undocumented argument [TSM_StartEndMode]
 * Arg:          wing_length [UNKN ] Undocumented argument [int]
 * Arg:        internal_bias [UNKN ] Undocumented argument [Probability]
 *
 */
void Wise2_set_startend_policy_ThreeStateModel(ThreeStateModel * tsm,TSM_StartEndMode mode,int wing_length,Probability internal_bias);
#define set_startend_policy_ThreeStateModel Wise2_set_startend_policy_ThreeStateModel


/* Function:  force_endbias_model(tsm,startend,internal)
 *
 * Descrip:    Makes start/end on probability and
 *             internal another
 *
 *             not probabilistically correct
 *
 *
 * Arg:             tsm [UNKN ] Undocumented argument [ThreeStateModel *]
 * Arg:        startend [UNKN ] Undocumented argument [double]
 * Arg:        internal [UNKN ] Undocumented argument [double]
 *
 */
void Wise2_force_endbias_model(ThreeStateModel * tsm,double startend,double internal);
#define force_endbias_model Wise2_force_endbias_model


/* Function:  force_global_model(tsm,prob_into_model)
 *
 * Descrip:    Makes start at position 0 and end at position end,
 *             no other positions being valid
 *
 *
 *
 * Arg:                    tsm [UNKN ] ThreeStateModel to be 'forced' [ThreeStateModel *]
 * Arg:        prob_into_model [UNKN ] Probability to start the model: for true global will be 1.0 [double]
 *
 */
void Wise2_force_global_model(ThreeStateModel * tsm,double prob_into_model) ;
#define force_global_model Wise2_force_global_model


/* Function:  force_wing_local_model(tsm,terminus_prob,wing_length)
 *
 * Descrip:    Places the balanace of the start probability
 *             at the start and then the rest spread evenly
 *             descending over the wing stretch of sequences
 *
 *             Same with the end probability
 *
 *
 * Arg:                  tsm [UNKN ] Undocumented argument [ThreeStateModel *]
 * Arg:        terminus_prob [UNKN ] the amount of probability to put on the real start and end [double]
 * Arg:          wing_length [UNKN ] the rest of the probability spread over this distance into the wing [int]
 *
 */
void Wise2_force_wing_local_model(ThreeStateModel * tsm,double terminus_prob,int wing_length);
#define force_wing_local_model Wise2_force_wing_local_model


/* Function:  force_weighted_local_model(tsm,prob_into_model,ratio_start,ratio_end)
 *
 * Descrip:    places the ratio of probability to start/end,
 *             and then distributes the rest over the start/end
 *
 *
 *
 * Arg:                    tsm [UNKN ] ThreeStateModel to be 'forced' [ThreeStateModel *]
 * Arg:        prob_into_model [UNKN ] Probability to start the model: for true global will be 1.0 [double]
 * Arg:            ratio_start [UNKN ] ratio of prob to unit 0 to the rest (1.0 means all goes to start) [double]
 * Arg:              ratio_end [UNKN ] ratio of prob to unit (last) to the rest (1.0 means all goes to the end) [double]
 *
 */
void Wise2_force_weighted_local_model(ThreeStateModel * tsm,double prob_into_model,double ratio_start,double ratio_end) ;
#define force_weighted_local_model Wise2_force_weighted_local_model


/* Function:  force_hmmfs_ThreeStateModel(tsm)
 *
 * Descrip:    places the probability at start end to precisely match
 *             hmmfs code.
 *
 *
 *
 * Arg:        tsm [UNKN ] ThreeStateModel to be 'forced' [ThreeStateModel *]
 *
 */
void Wise2_force_hmmfs_ThreeStateModel(ThreeStateModel * tsm);
#define force_hmmfs_ThreeStateModel Wise2_force_hmmfs_ThreeStateModel


/* Function:  ThreeStateScore_from_ThreeStateModel(tsm)
 *
 * Descrip:    Converts the three probability form to the score form.
 *
 *             There is real complications to this due to the fact that the prob
 *             model is a "push" model, that is MATCH2MATCH[i] means the match from i
 *             to i+1, whereas the DP routines rely on a "pull" model, ie
 *             MATCH2MATCH[i] is i-1 to i.
 *
 *             This routines does the conversion from push to pull, essentially offsetting
 *             the Match2 and Delete2 movements.
 *
 *
 * Arg:        tsm [READ ] ThreeStateModel probability form [ThreeStateModel *]
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateScore     *]
 *
 */
ThreeStateScore     * Wise2_ThreeStateScore_from_ThreeStateModel(ThreeStateModel * tsm);
#define ThreeStateScore_from_ThreeStateModel Wise2_ThreeStateScore_from_ThreeStateModel


/* Function:  show_ThreeStateModel(tsm,ofp)
 *
 * Descrip:    shows pretty ugly debugging type format.
 *
 *             Pretty useless
 *
 *
 * Arg:        tsm [UNKN ] Undocumented argument [ThreeStateModel *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
void Wise2_show_ThreeStateModel(ThreeStateModel * tsm,FILE * ofp);
#define show_ThreeStateModel Wise2_show_ThreeStateModel


/* Function:  pseudo_Protein_from_ThreeStateModel(tsm)
 *
 * Descrip:    Makes a protein sequence out of the display characters.
 *             Not very useful!
 *
 *
 *
 * Arg:        tsm [UNKN ] Undocumented argument [ThreeStateModel *]
 *
 * Return [UNKN ]  Undocumented return value [Protein *]
 *
 */
Protein * Wise2_pseudo_Protein_from_ThreeStateModel(ThreeStateModel * tsm);
#define pseudo_Protein_from_ThreeStateModel Wise2_pseudo_Protein_from_ThreeStateModel


/* Function:  ThreeStateModel_from_half_bit_Sequence(pro,mat,rm,gap,ext)
 *
 * Descrip:    Makes a local three-state-model from a sequence.  this is scary
 *             hackery, assumming that the matrix is half-bits and normalising in a
 *             *very* wrong way to get "probabilities" out.
 *
 *             Works though
 *
 *
 * Arg:        pro [READ ] protein sequence [Protein *]
 * Arg:        mat [READ ] comparison matrix to use [CompMat *]
 * Arg:         rm [READ ] random model which you assumme the matrix was built with [RandomModel *]
 * Arg:        gap [READ ] gap open penalty [int]
 * Arg:        ext [READ ] gap ext penalty [int]
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateModel *]
 *
 */
ThreeStateModel * Wise2_ThreeStateModel_from_half_bit_Sequence(Protein * pro,CompMat * mat,RandomModel * rm,int gap,int ext);
#define ThreeStateModel_from_half_bit_Sequence Wise2_ThreeStateModel_from_half_bit_Sequence


/* Function:  global_ThreeStateModel_from_half_bit_Sequence(pro,mat,rm,gap,ext)
 *
 * Descrip:    Makes a global three-state-model from a sequence.
 *             Like the local version, this is scary hackery, assumming
 *             that the matrix is half-bits and normalising in a *very*
 *             wrong way to get "probabilities" out.
 *
 *             Works though
 *
 *
 * Arg:        pro [READ ] protein sequence [Protein *]
 * Arg:        mat [READ ] comparison matrix to use [CompMat *]
 * Arg:         rm [READ ] random model which you assumme the matrix was built with [RandomModel *]
 * Arg:        gap [READ ] gap open penalty [int]
 * Arg:        ext [READ ] gap ext penalty [int]
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateModel *]
 *
 */
ThreeStateModel * Wise2_global_ThreeStateModel_from_half_bit_Sequence(Protein * pro,CompMat * mat,RandomModel * rm,int gap,int ext);
#define global_ThreeStateModel_from_half_bit_Sequence Wise2_global_ThreeStateModel_from_half_bit_Sequence


/* Function:  display_char_in_ThreeStateModel(tsm)
 *
 * Descrip:    Makes the display chars in the tsm come from
 *             the most likely amino acid in the match position
 *
 *
 * Arg:        tsm [UNKN ] Undocumented argument [ThreeStateModel *]
 *
 */
void Wise2_display_char_in_ThreeStateModel(ThreeStateModel * tsm);
#define display_char_in_ThreeStateModel Wise2_display_char_in_ThreeStateModel


/* Function:  fold_RandomModel_into_ThreeStateModel(tsm,rm)
 *
 * Descrip:    divides the emission and insertion scores
 *             by the random model for use in log-odd
 *             implementations
 *
 *
 * Arg:        tsm [UNKN ] Undocumented argument [ThreeStateModel *]
 * Arg:         rm [UNKN ] Undocumented argument [RandomModel *]
 *
 */
void Wise2_fold_RandomModel_into_ThreeStateModel(ThreeStateModel * tsm,RandomModel * rm);
#define fold_RandomModel_into_ThreeStateModel Wise2_fold_RandomModel_into_ThreeStateModel


/* Function:  add_sensible_start_end_global_for_HMMer(tsm)
 *
 * Descrip:    added start 0.5, 0.25,0.25 for hmmls type
 *
 *             Deprecated
 *
 *
 * Arg:        tsm [UNKN ] Undocumented argument [ThreeStateModel *]
 *
 */
void Wise2_add_sensible_start_end_global_for_HMMer(ThreeStateModel * tsm);
#define add_sensible_start_end_global_for_HMMer Wise2_add_sensible_start_end_global_for_HMMer


/* Function:  add_sensible_start_end_looping_for_HMMer(tsm,loop)
 *
 * Descrip:    adds start/end for loop probabilities coming
 *             in
 *
 *             Deprecated
 *
 *
 * Arg:         tsm [UNKN ] Undocumented argument [ThreeStateModel *]
 * Arg:        loop [UNKN ] Undocumented argument [Probability]
 *
 */
void Wise2_add_sensible_start_end_looping_for_HMMer(ThreeStateModel * tsm,Probability loop);
#define add_sensible_start_end_looping_for_HMMer Wise2_add_sensible_start_end_looping_for_HMMer


/* Function:  write_HMF_ThreeStateModel(tsm,ofp)
 *
 * Descrip:    writes the HMF format
 *
 *
 * Arg:        tsm [UNKN ] Undocumented argument [ThreeStateModel *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
void Wise2_write_HMF_ThreeStateModel(ThreeStateModel * tsm,FILE * ofp);
#define write_HMF_ThreeStateModel Wise2_write_HMF_ThreeStateModel


/* Function:  read_HMF_ThreeStateModel(ifp)
 *
 * Descrip:    reads the HMF format. Leaves the file ready
 *             to read the next format
 *
 *
 * Arg:        ifp [UNKN ] Undocumented argument [FILE *]
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateModel *]
 *
 */
ThreeStateModel * Wise2_read_HMF_ThreeStateModel(FILE * ifp);
#define read_HMF_ThreeStateModel Wise2_read_HMF_ThreeStateModel


/* Function:  write_HMMer_1_7_ascii_ThreeStateModel(tsm,ofp)
 *
 * Descrip:    writes a HMMer version 1.7 (also ok with 1.8) file
 *
 *
 * Arg:        tsm [UNKN ] Undocumented argument [ThreeStateModel *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
void Wise2_write_HMMer_1_7_ascii_ThreeStateModel(ThreeStateModel * tsm,FILE * ofp);
#define write_HMMer_1_7_ascii_ThreeStateModel Wise2_write_HMMer_1_7_ascii_ThreeStateModel


/* Function:  read_HMMer_1_7_ascii_file(filename)
 *
 * Descrip:    reads a HMMer ascii version 1.7 (1.8) file from filename.
 *
 *
 *
 * Arg:        filename [UNKN ] the name fo the hmmer file [char *]
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateModel *]
 *
 */
ThreeStateModel * Wise2_read_HMMer_1_7_ascii_file(char * filename);
#define read_HMMer_1_7_ascii_file Wise2_read_HMMer_1_7_ascii_file


/* Function:  read_HMMer_1_7_ascii(ifp)
 *
 * Descrip:    Basic function to read HMMer version 1.7(1.8) files.
 *
 *
 * Arg:        ifp [UNKN ] Undocumented argument [FILE *]
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateModel *]
 *
 */
ThreeStateModel * Wise2_read_HMMer_1_7_ascii(FILE * ifp) ;
#define read_HMMer_1_7_ascii Wise2_read_HMMer_1_7_ascii


/* Function:  hard_link_ThreeStateUnit(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [ThreeStateUnit *]
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateUnit *]
 *
 */
ThreeStateUnit * Wise2_hard_link_ThreeStateUnit(ThreeStateUnit * obj);
#define hard_link_ThreeStateUnit Wise2_hard_link_ThreeStateUnit


/* Function:  ThreeStateUnit_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateUnit *]
 *
 */
ThreeStateUnit * Wise2_ThreeStateUnit_alloc(void);
#define ThreeStateUnit_alloc Wise2_ThreeStateUnit_alloc


/* Function:  free_ThreeStateUnit(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [ThreeStateUnit *]
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateUnit *]
 *
 */
ThreeStateUnit * Wise2_free_ThreeStateUnit(ThreeStateUnit * obj);
#define free_ThreeStateUnit Wise2_free_ThreeStateUnit


/* Function:  add_ThreeStateModel(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [ThreeStateModel *]
 * Arg:        add [OWNER] Object to add to the list [ThreeStateUnit *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_ThreeStateModel(ThreeStateModel * obj,ThreeStateUnit * add);
#define add_ThreeStateModel Wise2_add_ThreeStateModel


/* Function:  flush_ThreeStateModel(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [ThreeStateModel *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_ThreeStateModel(ThreeStateModel * obj);
#define flush_ThreeStateModel Wise2_flush_ThreeStateModel


/* Function:  ThreeStateModel_alloc_std(void)
 *
 * Descrip:    Equivalent to ThreeStateModel_alloc_len(ThreeStateModelLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateModel *]
 *
 */
ThreeStateModel * Wise2_ThreeStateModel_alloc_std(void);
#define ThreeStateModel_alloc_std Wise2_ThreeStateModel_alloc_std


/* Function:  ThreeStateModel_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateModel *]
 *
 */
ThreeStateModel * Wise2_ThreeStateModel_alloc_len(int len);
#define ThreeStateModel_alloc_len Wise2_ThreeStateModel_alloc_len


/* Function:  hard_link_ThreeStateModel(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [ThreeStateModel *]
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateModel *]
 *
 */
ThreeStateModel * Wise2_hard_link_ThreeStateModel(ThreeStateModel * obj);
#define hard_link_ThreeStateModel Wise2_hard_link_ThreeStateModel


/* Function:  ThreeStateModel_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateModel *]
 *
 */
ThreeStateModel * Wise2_ThreeStateModel_alloc(void);
#define ThreeStateModel_alloc Wise2_ThreeStateModel_alloc


/* Function:  free_ThreeStateModel(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [ThreeStateModel *]
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateModel *]
 *
 */
ThreeStateModel * Wise2_free_ThreeStateModel(ThreeStateModel * obj);
#define free_ThreeStateModel Wise2_free_ThreeStateModel


/* Function:  hard_link_ThreeStateScoreUnit(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [ThreeStateScoreUnit *]
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateScoreUnit *]
 *
 */
ThreeStateScoreUnit * Wise2_hard_link_ThreeStateScoreUnit(ThreeStateScoreUnit * obj);
#define hard_link_ThreeStateScoreUnit Wise2_hard_link_ThreeStateScoreUnit


/* Function:  ThreeStateScoreUnit_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateScoreUnit *]
 *
 */
ThreeStateScoreUnit * Wise2_ThreeStateScoreUnit_alloc(void);
#define ThreeStateScoreUnit_alloc Wise2_ThreeStateScoreUnit_alloc


/* Function:  free_ThreeStateScoreUnit(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [ThreeStateScoreUnit *]
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateScoreUnit *]
 *
 */
ThreeStateScoreUnit * Wise2_free_ThreeStateScoreUnit(ThreeStateScoreUnit * obj);
#define free_ThreeStateScoreUnit Wise2_free_ThreeStateScoreUnit


/* Function:  add_ThreeStateScore(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [ThreeStateScore *]
 * Arg:        add [OWNER] Object to add to the list [ThreeStateScoreUnit *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_ThreeStateScore(ThreeStateScore * obj,ThreeStateScoreUnit * add);
#define add_ThreeStateScore Wise2_add_ThreeStateScore


/* Function:  flush_ThreeStateScore(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [ThreeStateScore *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_ThreeStateScore(ThreeStateScore * obj);
#define flush_ThreeStateScore Wise2_flush_ThreeStateScore


/* Function:  ThreeStateScore_alloc_std(void)
 *
 * Descrip:    Equivalent to ThreeStateScore_alloc_len(ThreeStateScoreLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateScore *]
 *
 */
ThreeStateScore * Wise2_ThreeStateScore_alloc_std(void);
#define ThreeStateScore_alloc_std Wise2_ThreeStateScore_alloc_std


/* Function:  ThreeStateScore_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateScore *]
 *
 */
ThreeStateScore * Wise2_ThreeStateScore_alloc_len(int len);
#define ThreeStateScore_alloc_len Wise2_ThreeStateScore_alloc_len


/* Function:  hard_link_ThreeStateScore(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [ThreeStateScore *]
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateScore *]
 *
 */
ThreeStateScore * Wise2_hard_link_ThreeStateScore(ThreeStateScore * obj);
#define hard_link_ThreeStateScore Wise2_hard_link_ThreeStateScore


/* Function:  ThreeStateScore_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateScore *]
 *
 */
ThreeStateScore * Wise2_ThreeStateScore_alloc(void);
#define ThreeStateScore_alloc Wise2_ThreeStateScore_alloc


/* Function:  free_ThreeStateScore(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [ThreeStateScore *]
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateScore *]
 *
 */
ThreeStateScore * Wise2_free_ThreeStateScore(ThreeStateScore * obj);
#define free_ThreeStateScore Wise2_free_ThreeStateScore


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
boolean Wise2_replace_name_ThreeStateModel(ThreeStateModel * obj,char * name);
#define replace_name_ThreeStateModel Wise2_replace_name_ThreeStateModel
char * Wise2_access_name_ThreeStateModel(ThreeStateModel * obj);
#define access_name_ThreeStateModel Wise2_access_name_ThreeStateModel
int Wise2_length_unit_ThreeStateModel(ThreeStateModel * obj);
#define length_unit_ThreeStateModel Wise2_length_unit_ThreeStateModel
boolean Wise2_replace_alphabet_ThreeStateModel(ThreeStateModel * obj,char * alphabet);
#define replace_alphabet_ThreeStateModel Wise2_replace_alphabet_ThreeStateModel
char Wise2_access_display_char_ThreeStateUnit(ThreeStateUnit * obj);
#define access_display_char_ThreeStateUnit Wise2_access_display_char_ThreeStateUnit
char * Wise2_access_alphabet_ThreeStateModel(ThreeStateModel * obj);
#define access_alphabet_ThreeStateModel Wise2_access_alphabet_ThreeStateModel
boolean Wise2_replace_display_char_ThreeStateUnit(ThreeStateUnit * obj,char display_char);
#define replace_display_char_ThreeStateUnit Wise2_replace_display_char_ThreeStateUnit
boolean Wise2_replace_accession_ThreeStateModel(ThreeStateModel * obj,char * accession);
#define replace_accession_ThreeStateModel Wise2_replace_accession_ThreeStateModel
RandomModel * Wise2_access_rm_ThreeStateModel(ThreeStateModel * obj);
#define access_rm_ThreeStateModel Wise2_access_rm_ThreeStateModel
char * Wise2_access_accession_ThreeStateModel(ThreeStateModel * obj);
#define access_accession_ThreeStateModel Wise2_access_accession_ThreeStateModel
ThreeStateUnit * Wise2_access_unit_ThreeStateModel(ThreeStateModel * obj,int i);
#define access_unit_ThreeStateModel Wise2_access_unit_ThreeStateModel
boolean Wise2_replace_threshold_ThreeStateModel(ThreeStateModel * obj,double threshold);
#define replace_threshold_ThreeStateModel Wise2_replace_threshold_ThreeStateModel
boolean Wise2_replace_rm_ThreeStateModel(ThreeStateModel * obj,RandomModel * rm);
#define replace_rm_ThreeStateModel Wise2_replace_rm_ThreeStateModel
double Wise2_access_threshold_ThreeStateModel(ThreeStateModel * obj);
#define access_threshold_ThreeStateModel Wise2_access_threshold_ThreeStateModel
ThreeStateScoreUnit * Wise2_ThreeStateScoreUnit_from_ThreeStateUnit(ThreeStateUnit * prev,ThreeStateUnit * tsu);
#define ThreeStateScoreUnit_from_ThreeStateUnit Wise2_ThreeStateScoreUnit_from_ThreeStateUnit
void Wise2_show_ThreeStateUnit(ThreeStateUnit * tsu,FILE * ofp);
#define show_ThreeStateUnit Wise2_show_ThreeStateUnit
ThreeStateUnit * Wise2_ThreeStateUnit_from_half_bit_aminoacid(char aa,CompMat * mat,RandomModel * rm,int gap,int ext);
#define ThreeStateUnit_from_half_bit_aminoacid Wise2_ThreeStateUnit_from_half_bit_aminoacid
char Wise2_display_char_for_ThreeStateUnit(ThreeStateUnit * tsu) ;
#define display_char_for_ThreeStateUnit Wise2_display_char_for_ThreeStateUnit
void Wise2_write_HMF_ThreeStateUnit(ThreeStateUnit * tsu,FILE * ofp);
#define write_HMF_ThreeStateUnit Wise2_write_HMF_ThreeStateUnit
ThreeStateUnit * Wise2_read_HMF_ThreeStateUnit(char * line,FILE * ifp);
#define read_HMF_ThreeStateUnit Wise2_read_HMF_ThreeStateUnit
void Wise2_write_HMMer_1_7_ascii_ThreeStateUnit(ThreeStateUnit * tsu,int no,FILE * ofp);
#define write_HMMer_1_7_ascii_ThreeStateUnit Wise2_write_HMMer_1_7_ascii_ThreeStateUnit
ThreeStateModel * Wise2_convert_push_model_to_pull_model(ThreeStateModel * tsm);
#define convert_push_model_to_pull_model Wise2_convert_push_model_to_pull_model
ThreeStateUnit  * Wise2_convert_push_unit_to_pull_unit(ThreeStateUnit * prev,ThreeStateUnit * this);
#define convert_push_unit_to_pull_unit Wise2_convert_push_unit_to_pull_unit
ThreeStateUnit * Wise2_read_HMMer_1_7_ThreeStateUnit(char * alphabet,FILE * ifp) ;
#define read_HMMer_1_7_ThreeStateUnit Wise2_read_HMMer_1_7_ThreeStateUnit
ThreeStateUnit * Wise2_blank_ThreeStateUnit(void);
#define blank_ThreeStateUnit Wise2_blank_ThreeStateUnit
void Wise2_swap_ThreeStateModel(ThreeStateUnit ** list,int i,int j) ;
#define swap_ThreeStateModel Wise2_swap_ThreeStateModel
void Wise2_qsort_ThreeStateModel(ThreeStateUnit ** list,int left,int right,int (*comp)(ThreeStateUnit * ,ThreeStateUnit * ));
#define qsort_ThreeStateModel Wise2_qsort_ThreeStateModel
void Wise2_sort_ThreeStateModel(ThreeStateModel * obj,int (*comp)(ThreeStateUnit *, ThreeStateUnit *));
#define sort_ThreeStateModel Wise2_sort_ThreeStateModel
boolean Wise2_expand_ThreeStateModel(ThreeStateModel * obj,int len);
#define expand_ThreeStateModel Wise2_expand_ThreeStateModel
void Wise2_swap_ThreeStateScore(ThreeStateScoreUnit ** list,int i,int j) ;
#define swap_ThreeStateScore Wise2_swap_ThreeStateScore
void Wise2_qsort_ThreeStateScore(ThreeStateScoreUnit ** list,int left,int right,int (*comp)(ThreeStateScoreUnit * ,ThreeStateScoreUnit * ));
#define qsort_ThreeStateScore Wise2_qsort_ThreeStateScore
void Wise2_sort_ThreeStateScore(ThreeStateScore * obj,int (*comp)(ThreeStateScoreUnit *, ThreeStateScoreUnit *));
#define sort_ThreeStateScore Wise2_sort_ThreeStateScore
boolean Wise2_expand_ThreeStateScore(ThreeStateScore * obj,int len);
#define expand_ThreeStateScore Wise2_expand_ThreeStateScore

#ifdef _cplusplus
}
#endif

#endif
