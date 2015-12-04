#ifndef DYNAMITEgenefrequencyHEADERFILE
#define DYNAMITEgenefrequencyHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "probability.h"
#include "codon.h"
#include "wisebase.h"
#include "dyna.h"
#include "randommodel.h"

enum {
  GF21_CENTRAL_STAY,
  GF21_PY_STAY,
  GF21_SPACER_STAY,
  GF21_NO_SPACER,
  GF21_INTRON_CORR_TERM,
  GENEFREQUENCY21_TRANSITION_LEN
};

#define GeneConsensusLISTLENGTH 128

enum {
  GeneConsensusType_5SS = 131,
  GeneConsensusType_3SS,
  GeneConsensusType_CDS,
  GeneConsensusType_Intron_Corr_Term,
  GeneConsensusType_Intron_emission,
  GeneConsensusType_Pyrimidine_emission,
  GeneConsensusType_Spacer_emission,
  GeneConsensusType_Central_stay,
  GeneConsensusType_Pyrimidine_stay,
  GeneConsensusType_Spacer_stay,
  GeneConsensusType_No_spacer,
  GeneConsensusType_Error
};


struct Wise2_GeneSingleCons {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    char * string;   
    double number;   
    } ;  
/* GeneSingleCons defined */ 
#ifndef DYNAMITE_DEFINED_GeneSingleCons
typedef struct Wise2_GeneSingleCons Wise2_GeneSingleCons;
#define GeneSingleCons Wise2_GeneSingleCons
#define DYNAMITE_DEFINED_GeneSingleCons
#endif


struct Wise2_GeneConsensus {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int center;  
    GeneSingleCons ** gsc;   
    int len;/* len for above gsc  */ 
    int maxlen; /* maxlen for above gsc */ 
    } ;  
/* GeneConsensus defined */ 
#ifndef DYNAMITE_DEFINED_GeneConsensus
typedef struct Wise2_GeneConsensus Wise2_GeneConsensus;
#define GeneConsensus Wise2_GeneConsensus
#define DYNAMITE_DEFINED_GeneConsensus
#endif


struct Wise2_GeneFrequency21 {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    GeneConsensus * ss5;     
    GeneConsensus * ss3;     
    double codon[64];    
    double central[4];   
    double py[4];    
    double spacer[4];    
    double transition[GENEFREQUENCY21_TRANSITION_LEN];   
    double cds_triplet[64]; /*  phase 0 */ 
    } ;  
/* GeneFrequency21 defined */ 
#ifndef DYNAMITE_DEFINED_GeneFrequency21
typedef struct Wise2_GeneFrequency21 Wise2_GeneFrequency21;
#define GeneFrequency21 Wise2_GeneFrequency21
#define DYNAMITE_DEFINED_GeneFrequency21
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  RandomCodon_from_cds_triplet(gf)
 *
 * Descrip:    makes a randomcodon probability emission
 *             from the counts in genefrequency
 *
 *
 * Arg:        gf [UNKN ] Undocumented argument [GeneFrequency21 *]
 *
 * Return [UNKN ]  Undocumented return value [RandomCodon *]
 *
 */
RandomCodon * Wise2_RandomCodon_from_cds_triplet(GeneFrequency21 * gf);
#define RandomCodon_from_cds_triplet Wise2_RandomCodon_from_cds_triplet


/* Function:  RandomModelDNA_from_central_GeneFrequency21(gf)
 *
 * Descrip:    Makes a random model from the central gene 
 *             model of an intron. Ideal for tieing intron
 *             state distribution to the randommodel
 *
 *
 * Arg:        gf [UNKN ] Undocumented argument [GeneFrequency21 *]
 *
 * Return [UNKN ]  Undocumented return value [RandomModelDNA *]
 *
 */
RandomModelDNA * Wise2_RandomModelDNA_from_central_GeneFrequency21(GeneFrequency21 * gf);
#define RandomModelDNA_from_central_GeneFrequency21 Wise2_RandomModelDNA_from_central_GeneFrequency21


/* Function:  ComplexConsensi_5SS_from_GeneFrequency(gf)
 *
 * Descrip:    makes 5'SS ComplexConsensi from GeneFrequency21 structure using
 *
 *               CCC|XXXXXXX score = no(5'SS with CCC|XXXXXXX) / no(CCC in cds).
 *
 *
 * Arg:        gf [UNKN ] Undocumented argument [GeneFrequency21 *]
 *
 * Return [UNKN ]  Undocumented return value [ComplexConsensi *]
 *
 */
ComplexConsensi * Wise2_ComplexConsensi_5SS_from_GeneFrequency(GeneFrequency21 * gf);
#define ComplexConsensi_5SS_from_GeneFrequency Wise2_ComplexConsensi_5SS_from_GeneFrequency


/* Function:  ComplexConsensi_3SS_from_GeneFrequency(gf)
 *
 * Descrip:    makes 3'SS ComplexConsensi from GeneFrequency21 structure using
 *
 *               ZZZ|CCC score = no(3'SS with ZZZ|CCC) / no(CCC in cds).
 *
 *
 * Arg:        gf [UNKN ] Undocumented argument [GeneFrequency21 *]
 *
 * Return [UNKN ]  Undocumented return value [ComplexConsensi *]
 *
 */
ComplexConsensi * Wise2_ComplexConsensi_3SS_from_GeneFrequency(GeneFrequency21 * gf);
#define ComplexConsensi_3SS_from_GeneFrequency Wise2_ComplexConsensi_3SS_from_GeneFrequency


/* Function:  CodonFrequency_from_GeneFrequency21(gf,ct)
 *
 * Descrip:    Builds a codon frequency table from raw counts
 *             in the counts file
 *
 *
 * Arg:        gf [UNKN ] Undocumented argument [GeneFrequency21 *]
 * Arg:        ct [UNKN ] Undocumented argument [CodonTable *]
 *
 * Return [UNKN ]  Undocumented return value [CodonFrequency *]
 *
 */
CodonFrequency * Wise2_CodonFrequency_from_GeneFrequency21(GeneFrequency21 * gf,CodonTable * ct);
#define CodonFrequency_from_GeneFrequency21 Wise2_CodonFrequency_from_GeneFrequency21


/* Function:  read_GeneFrequency21_file(filename)
 *
 * Descrip:    Opens the file with /openfile
 *
 *             Reads in a GeneFrequency (Mor-Ewan style)
 *
 *
 *
 * Arg:        filename [UNKN ] will open from WISECONFIGDIR etc via openfile [char *]
 *
 * Return [UNKN ]  a newly allocated structure [GeneFrequency21 *]
 *
 */
GeneFrequency21 * Wise2_read_GeneFrequency21_file(char * filename);
#define read_GeneFrequency21_file Wise2_read_GeneFrequency21_file


/* Function:  read_GeneFrequency21(ifp)
 *
 * Descrip:    Reads in a GeneFrequency (Mor-Ewan style)
 *             file from ifp
 *
 *
 * Arg:        ifp [UNKN ] file pointer [FILE *]
 *
 * Return [UNKN ]  a newly allocated structure [GeneFrequency21 *]
 *
 */
GeneFrequency21 * Wise2_read_GeneFrequency21(FILE * ifp);
#define read_GeneFrequency21 Wise2_read_GeneFrequency21


/* Function:  double_from_line(buffer)
 *
 * Descrip:    helper string function
 *
 *
 * Arg:        buffer [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [double]
 *
 */
double Wise2_double_from_line(char * buffer);
#define double_from_line Wise2_double_from_line


/* Function:  check_type_GeneFrequency(*line,ifp,center,phase)
 *
 * Descrip:    Pretty sneaky function 
 *
 *             give 
 *
 *               line starting with "type xxx"
 *               ifp  file pointer
 *               centre a &int for returning the centre value of the consensus if any.
 *               phase  a &int for returning the phase value of the consensus if any.
 *
 *             you get *back* the line with the line "begin consensus" or "number" in it.
 *
 *               
 *             It returns a GeneConsensusType
 *
 *             with GeneConsensusType_Error on error
 *
 *
 *
 * Arg:         *line [UNKN ] Undocumented argument [char]
 * Arg:           ifp [UNKN ] Undocumented argument [FILE *]
 * Arg:        center [UNKN ] Undocumented argument [int *]
 * Arg:         phase [UNKN ] Undocumented argument [int *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_check_type_GeneFrequency(char *line,FILE * ifp,int * center,int * phase)  ;
#define check_type_GeneFrequency Wise2_check_type_GeneFrequency


/* Function:  read_base_GeneConsensus(base_array,line,ifp)
 *
 * Descrip:    assummes base_array is 4 positions long
 *               
 *             line should have begin consensus on it and be of MAXLINE length as it will be used as the buffer.
 *               
 *             This does **not** check that you have filled up all 4 positions.
 *
 *
 * Arg:        base_array [UNKN ] Undocumented argument [double *]
 * Arg:              line [UNKN ] Undocumented argument [char*]
 * Arg:               ifp [UNKN ] Undocumented argument [FILE *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_read_base_GeneConsensus(double * base_array,char* line,FILE * ifp);
#define read_base_GeneConsensus Wise2_read_base_GeneConsensus


/* Function:  read_codon_GeneConsensus(codon_array,line,ifp)
 *
 * Descrip:    assummes codon_array is 64 positions long
 *               
 *             line should have begin consensus on it and be of MAXLINE length as it will be used as the buffer.
 *
 *             This does **not** check that you have filled up all 64 positions.
 *
 *
 * Arg:        codon_array [UNKN ] Undocumented argument [double *]
 * Arg:               line [UNKN ] Undocumented argument [char*]
 * Arg:                ifp [UNKN ] Undocumented argument [FILE *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_read_codon_GeneConsensus(double * codon_array,char* line,FILE * ifp);
#define read_codon_GeneConsensus Wise2_read_codon_GeneConsensus


/* Function:  read_line_GeneConsensus(line,ifp)
 *
 * Descrip:    Reads a single GeneConsensus from a file
 *
 *
 * Arg:        line [UNKN ] Undocumented argument [char *]
 * Arg:         ifp [UNKN ] Undocumented argument [FILE *]
 *
 * Return [UNKN ]  Undocumented return value [GeneConsensus  *]
 *
 */
GeneConsensus  * Wise2_read_line_GeneConsensus(char * line,FILE * ifp);
#define read_line_GeneConsensus Wise2_read_line_GeneConsensus


/* Function:  hard_link_GeneSingleCons(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GeneSingleCons *]
 *
 * Return [UNKN ]  Undocumented return value [GeneSingleCons *]
 *
 */
GeneSingleCons * Wise2_hard_link_GeneSingleCons(GeneSingleCons * obj);
#define hard_link_GeneSingleCons Wise2_hard_link_GeneSingleCons


/* Function:  GeneSingleCons_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GeneSingleCons *]
 *
 */
GeneSingleCons * Wise2_GeneSingleCons_alloc(void);
#define GeneSingleCons_alloc Wise2_GeneSingleCons_alloc


/* Function:  free_GeneSingleCons(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GeneSingleCons *]
 *
 * Return [UNKN ]  Undocumented return value [GeneSingleCons *]
 *
 */
GeneSingleCons * Wise2_free_GeneSingleCons(GeneSingleCons * obj);
#define free_GeneSingleCons Wise2_free_GeneSingleCons


/* Function:  add_GeneConsensus(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [GeneConsensus *]
 * Arg:        add [OWNER] Object to add to the list [GeneSingleCons *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_GeneConsensus(GeneConsensus * obj,GeneSingleCons * add);
#define add_GeneConsensus Wise2_add_GeneConsensus


/* Function:  flush_GeneConsensus(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [GeneConsensus *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_GeneConsensus(GeneConsensus * obj);
#define flush_GeneConsensus Wise2_flush_GeneConsensus


/* Function:  GeneConsensus_alloc_std(void)
 *
 * Descrip:    Equivalent to GeneConsensus_alloc_len(GeneConsensusLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GeneConsensus *]
 *
 */
GeneConsensus * Wise2_GeneConsensus_alloc_std(void);
#define GeneConsensus_alloc_std Wise2_GeneConsensus_alloc_std


/* Function:  GeneConsensus_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [GeneConsensus *]
 *
 */
GeneConsensus * Wise2_GeneConsensus_alloc_len(int len);
#define GeneConsensus_alloc_len Wise2_GeneConsensus_alloc_len


/* Function:  hard_link_GeneConsensus(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GeneConsensus *]
 *
 * Return [UNKN ]  Undocumented return value [GeneConsensus *]
 *
 */
GeneConsensus * Wise2_hard_link_GeneConsensus(GeneConsensus * obj);
#define hard_link_GeneConsensus Wise2_hard_link_GeneConsensus


/* Function:  GeneConsensus_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GeneConsensus *]
 *
 */
GeneConsensus * Wise2_GeneConsensus_alloc(void);
#define GeneConsensus_alloc Wise2_GeneConsensus_alloc


/* Function:  free_GeneConsensus(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GeneConsensus *]
 *
 * Return [UNKN ]  Undocumented return value [GeneConsensus *]
 *
 */
GeneConsensus * Wise2_free_GeneConsensus(GeneConsensus * obj);
#define free_GeneConsensus Wise2_free_GeneConsensus


/* Function:  hard_link_GeneFrequency21(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GeneFrequency21 *]
 *
 * Return [UNKN ]  Undocumented return value [GeneFrequency21 *]
 *
 */
GeneFrequency21 * Wise2_hard_link_GeneFrequency21(GeneFrequency21 * obj);
#define hard_link_GeneFrequency21 Wise2_hard_link_GeneFrequency21


/* Function:  GeneFrequency21_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GeneFrequency21 *]
 *
 */
GeneFrequency21 * Wise2_GeneFrequency21_alloc(void);
#define GeneFrequency21_alloc Wise2_GeneFrequency21_alloc


/* Function:  free_GeneFrequency21(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GeneFrequency21 *]
 *
 * Return [UNKN ]  Undocumented return value [GeneFrequency21 *]
 *
 */
GeneFrequency21 * Wise2_free_GeneFrequency21(GeneFrequency21 * obj);
#define free_GeneFrequency21 Wise2_free_GeneFrequency21


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
GeneSingleCons * Wise2_read_line_GeneSingleCons(char * line);
#define read_line_GeneSingleCons Wise2_read_line_GeneSingleCons
GeneFrequency21  * Wise2_untouched_GeneFrequency21(void);
#define untouched_GeneFrequency21 Wise2_untouched_GeneFrequency21


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
GeneConsensus * Wise2_access_ss5_GeneFrequency21(GeneFrequency21 * obj);
#define access_ss5_GeneFrequency21 Wise2_access_ss5_GeneFrequency21
boolean Wise2_replace_ss3_GeneFrequency21(GeneFrequency21 * obj,GeneConsensus * ss3);
#define replace_ss3_GeneFrequency21 Wise2_replace_ss3_GeneFrequency21
GeneConsensus * Wise2_access_ss3_GeneFrequency21(GeneFrequency21 * obj);
#define access_ss3_GeneFrequency21 Wise2_access_ss3_GeneFrequency21
boolean Wise2_replace_number_GeneSingleCons(GeneSingleCons * obj,double number);
#define replace_number_GeneSingleCons Wise2_replace_number_GeneSingleCons
char * Wise2_access_string_GeneSingleCons(GeneSingleCons * obj);
#define access_string_GeneSingleCons Wise2_access_string_GeneSingleCons
boolean Wise2_replace_center_GeneConsensus(GeneConsensus * obj,int center);
#define replace_center_GeneConsensus Wise2_replace_center_GeneConsensus
double Wise2_access_number_GeneSingleCons(GeneSingleCons * obj);
#define access_number_GeneSingleCons Wise2_access_number_GeneSingleCons
int Wise2_access_center_GeneConsensus(GeneConsensus * obj);
#define access_center_GeneConsensus Wise2_access_center_GeneConsensus
boolean Wise2_replace_ss5_GeneFrequency21(GeneFrequency21 * obj,GeneConsensus * ss5);
#define replace_ss5_GeneFrequency21 Wise2_replace_ss5_GeneFrequency21
GeneSingleCons * Wise2_access_gsc_GeneConsensus(GeneConsensus * obj,int i);
#define access_gsc_GeneConsensus Wise2_access_gsc_GeneConsensus
boolean Wise2_replace_string_GeneSingleCons(GeneSingleCons * obj,char * string);
#define replace_string_GeneSingleCons Wise2_replace_string_GeneSingleCons
int Wise2_length_gsc_GeneConsensus(GeneConsensus * obj);
#define length_gsc_GeneConsensus Wise2_length_gsc_GeneConsensus
double Wise2_nocds_from_ambiguous_codon(char * codon,double * codon_freq_array);
#define nocds_from_ambiguous_codon Wise2_nocds_from_ambiguous_codon
ComplexConsensusWord * Wise2_ComplexConsensusWord_from_string_and_prob(char * string,Probability p);
#define ComplexConsensusWord_from_string_and_prob Wise2_ComplexConsensusWord_from_string_and_prob
void Wise2_show_flat_GeneFrequency21(GeneFrequency21 * gf21,FILE * ofp);
#define show_flat_GeneFrequency21 Wise2_show_flat_GeneFrequency21
void Wise2_show_codon_emission(double * codon,FILE * ofp);
#define show_codon_emission Wise2_show_codon_emission
void Wise2_show_single_codon_emission(double no,int base4codon,FILE * ofp);
#define show_single_codon_emission Wise2_show_single_codon_emission
void Wise2_show_base_emission(double * base,FILE * ofp);
#define show_base_emission Wise2_show_base_emission
void Wise2_show_GeneConsensus(GeneConsensus * gc,FILE * ofp);
#define show_GeneConsensus Wise2_show_GeneConsensus
void Wise2_show_GeneSingleCons(GeneSingleCons * gsc,FILE * ofp);
#define show_GeneSingleCons Wise2_show_GeneSingleCons
boolean Wise2_skip_consensus(FILE * ifp);
#define skip_consensus Wise2_skip_consensus
int Wise2_string_to_GeneConsensusType(char * string);
#define string_to_GeneConsensusType Wise2_string_to_GeneConsensusType
void Wise2_swap_GeneConsensus(GeneSingleCons ** list,int i,int j) ;
#define swap_GeneConsensus Wise2_swap_GeneConsensus
void Wise2_qsort_GeneConsensus(GeneSingleCons ** list,int left,int right,int (*comp)(GeneSingleCons * ,GeneSingleCons * ));
#define qsort_GeneConsensus Wise2_qsort_GeneConsensus
void Wise2_sort_GeneConsensus(GeneConsensus * obj,int (*comp)(GeneSingleCons *, GeneSingleCons *));
#define sort_GeneConsensus Wise2_sort_GeneConsensus
boolean Wise2_expand_GeneConsensus(GeneConsensus * obj,int len);
#define expand_GeneConsensus Wise2_expand_GeneConsensus

#ifdef _cplusplus
}
#endif

#endif
