#ifndef DYNAMITEestwrapHEADERFILE
#define DYNAMITEestwrapHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "estwise3.h"
#include "estloop3.h"
#include "estfrag3.h"
#include "estslim3.h"
#include "estquick3.h"
#include "estslimloop.h"

typedef enum est_alg_type {
  ESTWISE_3,
  ESTLOOP_3,
  ESTSLIM_3,
  ESTQUICK_3,
  ESTFRAG_3,
  ESTSLIM_L
} est_alg_type;



    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  string_from_alg_estwrap(alg_type)
 *
 * Descrip:    Returns the string form of the algorithm
 *
 *
 * Arg:        alg_type [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [char *]
 *
 */
char * Wise2_string_from_alg_estwrap(int alg_type);
#define string_from_alg_estwrap Wise2_string_from_alg_estwrap


/* Function:  alg_estwrap_from_string(str)
 *
 * Descrip:    This function returns the algorithm type
 *             for an est search from the string
 *
 *
 * Arg:        str [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_alg_estwrap_from_string(char * str);
#define alg_estwrap_from_string Wise2_alg_estwrap_from_string


/* Function:  AlnBlock_from_Protein_estwise_wrap(pro,cdna,cp,cm,ct,comp,gap,ext,is_global,rmd,alg,rm,use_syn,allN,dpri,palpoi)
 *
 * Descrip:    This function is the guts for the est single alignment
 *             mode.
 *
 *             It uses /AlnBlock_from_TSM_estwise_wrap for the
 *             heavy part of the call
 *
 *
 * Arg:              pro [READ ] protein to be used in the comparison [Protein *]
 * Arg:             cdna [READ ] cdna to be compared to [cDNA *]
 * Arg:               cp [READ ] cdna parser indicating insertion deletion probabilities [cDNAParser *]
 * Arg:               cm [READ ] codon mapper indicating substitution errors etc [CodonMapper *]
 * Arg:               ct [READ ] codon table for the codon->amino acid mappings [CodonTable *]
 * Arg:             comp [READ ] comparison matrix to use [CompMat *]
 * Arg:              gap [UNKN ] gap penalty [int]
 * Arg:              ext [UNKN ] extension penalty [int]
 * Arg:        is_global [UNKN ] if true, global start-end in protein is used [boolean]
 * Arg:              rmd [UNKN ] random model of dna to use [RandomModelDNA *]
 * Arg:              alg [UNKN ] est algorithm type to use [int]
 * Arg:               rm [UNKN ] random protein model for use with compmat [RandomModel *]
 * Arg:          use_syn [UNKN ] if true, uses a synchronous coding model [boolean]
 * Arg:             allN [UNKN ] Undocumented argument [Probability]
 * Arg:             dpri [UNKN ] Undocumented argument [DPRunImpl *]
 * Arg:           palpoi [WRITE] the raw packed alignment output if wanted [PackAln **]
 *
 * Return [UNKN ]  Undocumented return value [AlnBlock *]
 *
 */
AlnBlock * Wise2_AlnBlock_from_Protein_estwise_wrap(Protein * pro,cDNA * cdna,cDNAParser * cp,CodonMapper * cm,CodonTable * ct,CompMat * comp,int gap,int ext,boolean is_global,RandomModelDNA * rmd,int alg,RandomModel * rm,boolean use_syn,Probability allN,DPRunImpl * dpri,PackAln ** palpoi);
#define AlnBlock_from_Protein_estwise_wrap Wise2_AlnBlock_from_Protein_estwise_wrap


/* Function:  AlnBlock_from_TSM_estwise_wrap(tsm,cdna,cp,cm,ct,rmd,alg,use_syn,force_flat_insert,allN,dpri,palpoi)
 *
 * Descrip:    This function is the basic wrap for Protein models
 *             vs cDNA sequences.
 *
 *
 * Arg:                      tsm [READ ] threestatemodel to be compared to the dna [ThreeStateModel *]
 * Arg:                     cdna [READ ] cdna to be compared to [cDNA *]
 * Arg:                       cp [READ ] cdna parser indicating insertion deletion probabilities [cDNAParser *]
 * Arg:                       cm [READ ] codon mapper indicating substitution errors etc [CodonMapper *]
 * Arg:                       ct [READ ] codon table for the codon->amino acid mappings [CodonTable *]
 * Arg:                      rmd [UNKN ] random model of dna to use [RandomModelDNA *]
 * Arg:                      alg [UNKN ] est algorithm type to use [int]
 * Arg:                  use_syn [UNKN ] if true, uses a synchronous coding model [boolean]
 * Arg:        force_flat_insert [UNKN ] Undocumented argument [boolean]
 * Arg:                     allN [UNKN ] Undocumented argument [Probability]
 * Arg:                     dpri [UNKN ] Undocumented argument [DPRunImpl *]
 * Arg:                   palpoi [WRITE] the raw packed alignment output if wanted [PackAln **]
 *
 * Return [UNKN ]  Undocumented return value [AlnBlock *]
 *
 */
AlnBlock * Wise2_AlnBlock_from_TSM_estwise_wrap(ThreeStateModel * tsm,cDNA * cdna,cDNAParser * cp,CodonMapper * cm,CodonTable * ct,RandomModelDNA * rmd,int alg,boolean use_syn,boolean force_flat_insert,Probability allN,DPRunImpl * dpri,PackAln ** palpoi);
#define AlnBlock_from_TSM_estwise_wrap Wise2_AlnBlock_from_TSM_estwise_wrap


/* Function:  Hscore_from_TSM_estwise(tdb,cdb,cp,cm,rmd,use_syn,alg,bits_cutoff,allN,flat_insert,report_level,die_on_error,dbsi)
 *
 * Descrip:    Runs a database search for the estwise set
 *             of algorithms
 *
 *
 * Arg:                 tdb [READ ] a three state model database [ThreeStateDB *]
 * Arg:                 cdb [READ ] a dna sequence database [cDNADB *]
 * Arg:                  cp [READ ] the codon parser for this comparison [cDNAParser *]
 * Arg:                  cm [READ ] the codon mapper for this comparison [CodonMapper *]
 * Arg:                 rmd [READ ] random model used for the dna sequence comparison [RandomModelDNA *]
 * Arg:             use_syn [UNKN ] whether a synchronous coding model should be used or not [boolean]
 * Arg:                 alg [UNKN ] algorithm to use [int]
 * Arg:         bits_cutoff [UNKN ] Undocumented argument [double]
 * Arg:                allN [UNKN ] Undocumented argument [Probability]
 * Arg:         flat_insert [UNKN ] Undocumented argument [boolean]
 * Arg:        report_level [UNKN ] Undocumented argument [int]
 * Arg:        die_on_error [UNKN ] if true, dies if there is an error [boolean]
 * Arg:                dbsi [UNKN ] Undocumented argument [DBSearchImpl *]
 *
 * Return [OWNER]  a newly allocated Hscore structure of the search [Hscore *]
 *
 */
Hscore * Wise2_Hscore_from_TSM_estwise(ThreeStateDB * tdb,cDNADB * cdb,cDNAParser * cp,CodonMapper * cm,RandomModelDNA * rmd,boolean use_syn,int alg,double bits_cutoff,Probability allN,boolean flat_insert,int report_level,boolean die_on_error,DBSearchImpl * dbsi);
#define Hscore_from_TSM_estwise Wise2_Hscore_from_TSM_estwise


/* Function:  write_mul_estwise_AlnBlock(alb,ct,ofp)
 *
 * Descrip:    writes an mul format protein multiple
 *             alignment from an AlnBlock with the first
 *             sequence an HMM/protein and ignored, and
 *             the second and subsequent sequences cdna sequences
 *             which are then translated into proteins
 *
 *             This relies considerably on the alb being made
 *             correctly, and if it is not, then god help you.
 *
 *             the estwisedb programs makes the alb correctly
 *
 *
 * Arg:        alb [UNKN ] Undocumented argument [AlnBlock *]
 * Arg:         ct [UNKN ] Undocumented argument [CodonTable *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
void Wise2_write_mul_estwise_AlnBlock(AlnBlock * alb,CodonTable * ct,FILE * ofp);
#define write_mul_estwise_AlnBlock Wise2_write_mul_estwise_AlnBlock


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
