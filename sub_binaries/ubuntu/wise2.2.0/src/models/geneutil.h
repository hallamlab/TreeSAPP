#ifndef DYNAMITEgeneutilHEADERFILE
#define DYNAMITEgeneutilHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "dyna.h"




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  new_GenomicRegion_from_GeneWise(gen,pseudo,alb)
 *
 * Descrip:    Makes a new GenomicRegion with the genes
 *             predicted from this AlnBlock
 *
 *             Really a wrapper around the add_Genes_to_GenomicRegion_GeneWise
 *             and other functionality
 *
 *
 * Arg:           gen [UNKN ] genomic sequence to use [Genomic *]
 * Arg:        pseudo [UNKN ] If true, predicts frameshifted genes as pseudogenes [boolean]
 * Arg:           alb [UNKN ] genewise alignment to predict genes from [AlnBlock *]
 *
 * Return [UNKN ]  a newly allocated structure [GenomicRegion *]
 *
 */
GenomicRegion * Wise2_new_GenomicRegion_from_GeneWise(Genomic * gen,boolean pseudo,AlnBlock * alb);
#define new_GenomicRegion_from_GeneWise Wise2_new_GenomicRegion_from_GeneWise


/* Function:  add_Genes_to_GenomicRegion_GeneWise(gr,org_start,org_end,alb,root,pseudo,make_name)
 *
 * Descrip:    Potential an Alnblock may have more
 *             than one gene due to looping models
 *
 *             This adds all the genes to gr
 *
 *
 *
 * Arg:               gr [UNKN ] genomic region to add genes to [GenomicRegion *]
 * Arg:        org_start [UNKN ] start point of the dna to which the alb was made from [int]
 * Arg:          org_end [UNKN ] end point of the dna to which the alb was made from [int]
 * Arg:              alb [UNKN ] logical label alignment [AlnBlock *]
 * Arg:             root [UNKN ] the second argument to make_name [char *]
 * Arg:           pseudo [UNKN ] If true, frameshifted genes are predicted as pseudo genes [boolean]
 * Arg:        make_name [FUNCP] a pointer to a function to actually make the name of the gene [char * (*make_name]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_add_Genes_to_GenomicRegion_GeneWise(GenomicRegion * gr,int org_start,int org_end,AlnBlock * alb,char * root,boolean pseudo,char * (*make_name)(Wise2_Genomic *,char *,int ,Wise2_Gene * ));
#define add_Genes_to_GenomicRegion_GeneWise Wise2_add_Genes_to_GenomicRegion_GeneWise


/* Function:  Gene_from_AlnColumn_GeneWise(alc,org_start,org_end,predict_pseudo_for_frameshift,end)
 *
 * Descrip:    A wrap for making a gene structure from
 *             an AlnBlock derived from one of the genewise
 *             methods
 *
 *
 * Arg:                                  alc [UNKN ] Alignment column in an AlnBlock produced by genewise [AlnColumn *]
 * Arg:                            org_start [UNKN ] offset that the genewise alignment was to the coordinate system [int]
 * Arg:                              org_end [UNKN ] emd that the genewise alignment was to the coordinate system [int]
 * Arg:        predict_pseudo_for_frameshift [UNKN ] Undocumented argument [boolean]
 * Arg:                                  end [WRITE] pointer to a AlnColumn * to say when it has finished with this gene [AlnColumn **]
 *
 * Return [UNKN ]  Undocumented return value [Gene *]
 *
 */
Gene * Wise2_Gene_from_AlnColumn_GeneWise(AlnColumn * alc,int org_start,int org_end,boolean predict_pseudo_for_frameshift,AlnColumn ** end);
#define Gene_from_AlnColumn_GeneWise Wise2_Gene_from_AlnColumn_GeneWise


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
boolean Wise2_is_frameshift_AlnColumn_genewise(const AlnColumn * alc);
#define is_frameshift_AlnColumn_genewise Wise2_is_frameshift_AlnColumn_genewise
boolean Wise2_is_random_AlnColumn_genewise(const AlnColumn * alc);
#define is_random_AlnColumn_genewise Wise2_is_random_AlnColumn_genewise
boolean Wise2_is_intron_AlnColumn_genewise(const AlnColumn * alc);
#define is_intron_AlnColumn_genewise Wise2_is_intron_AlnColumn_genewise

#ifdef _cplusplus
}
#endif

#endif
