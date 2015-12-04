#ifndef DYNAMITEgwrapHEADERFILE
#define DYNAMITEgwrapHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "dyna.h"
#include "genewise6.h"
#include "geneloop6.h"
#include "genestretch6.h"
#include "genewise4.h"
#include "gwlite.h"
#include "geneutil.h"

#include "genewise21.h"
#include "geneloop21.h" 

/**
#include "genelink21.h"
**/

#include "geneparameter.h"
#include "genefrequency.h"
#include "genestats.h"

#include "seqhit.h"
#include "estwise3.h"
#include "estloop3.h"

#include "genewisehsp.h"


enum GWWRAP_ALG_TYPE {
  GWWRAP_2193 = 12,
  GWWRAP_2193L,
  GWWRAP_2193I,
  GWWRAP_623,
  GWWRAP_623L,
  GWWRAP_421,
  GWWRAP_6LITE,
  GWWRAP_333,
  GWWRAP_333L,
  GWWRAP_623S,
  GWWRAP_623P
};




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  GeneParameter21_wrap(gf,subs_error,indel_error,rmd,use_modelled_codon,use_modelled_splice,tie_intron_prob,ct,rnd_loop,cds_loop,rnd_to_model,link_loop,link_to_model)
 *
 * Descrip:    A general wrap over the production of parameters for the
 *             GeneWise programs. The geneparameter21 holds all the parameters,
 *             and can be approximated for the 6:23 and 4:21 algorithms
 *
 *             This function is the best way to make a GeneParameter21 object
 *             as all the different options for how to make it or modify its
 *             contents are laid out as arguments to this function
 *
 *
 * Arg:                         gf [READ ] Gene Frequency data structure, holding counts for splice sites etc [GeneFrequency21 *]
 * Arg:                 subs_error [UNKN ] substitution error on the dna sequence [double]
 * Arg:                indel_error [UNKN ] rough estimate of the insertion/deletion per base error rate [double]
 * Arg:                        rmd [UNKN ] the random model of the DNA that is used [RandomModelDNA *]
 * Arg:         use_modelled_codon [UNKN ] if TRUE, model codon frequency [boolean]
 * Arg:        use_modelled_splice [UNKN ] if TRUE, make splice models from gf parameters [boolean]
 * Arg:            tie_intron_prob [UNKN ] Undocumented argument [boolean]
 * Arg:                         ct [UNKN ] codon table which is used for codon->aa mapping [CodonTable *]
 * Arg:                   rnd_loop [UNKN ] Undocumented argument [Probability]
 * Arg:                   cds_loop [UNKN ] Undocumented argument [Probability]
 * Arg:               rnd_to_model [UNKN ] Undocumented argument [Probability]
 * Arg:                  link_loop [UNKN ] Undocumented argument [Probability]
 * Arg:              link_to_model [UNKN ] Undocumented argument [Probability]
 *
 * Return [UNKN ]  A newly allocated structure [GeneParameter21 *]
 *
 */
GeneParameter21 * Wise2_GeneParameter21_wrap(GeneFrequency21 * gf,double subs_error,double indel_error,RandomModelDNA * rmd,boolean use_modelled_codon,boolean use_modelled_splice,boolean tie_intron_prob,CodonTable * ct,Probability rnd_loop,Probability cds_loop,Probability rnd_to_model,Probability link_loop,Probability link_to_model);
#define GeneParameter21_wrap Wise2_GeneParameter21_wrap


/* Function:  AlnBlock_from_protein_genewise_wrap(protein,pg,is_global,dna,comp,gap,ext,gpara,rmd,intergenic,alg,use_syn,rm,allN,startendmode,dpri,pal,gwp)
 *
 * Descrip:    A function which aligns a Protein sequecne to a Genomic sequence
 *             under the Comparison matrix comp and the gene paras in gpara.
 *
 *             This is the best function for accessing GeneWise functionality
 *             for a protein to dna comparison, allowing for introns.
 *
 *             To make the protein object, you will first read in a generic
 *             sequence object using something like read_fasta_Sequence and
 *             then convert it to a protein object using new_Protein_from_Sequence
 *
 *             To make the genomic object, you will first read in a generic
 *             sequence object using something like read_fasta_Sequence and
 *             then convert it to a genomic object using new_Genomic_from_Sequence
 *
 *             To make a CompMat object you will use read_Blast_file_CompMat
 *             from the compmat module. It is likely, if the Wise2 enviroment
 *             has been set up correctly that read_Blast_file_CompMat("blosum62.bla")
 *             will be fine. You should at the moment only use halfbit matrices
 *             (blosum62 is one such matrix)
 *
 *             To make the necessary random modules use the default construtors
 *             in the randommodel module
 *
 *             To make the gene parameter object use the GeneParameter21_wrap
 *             function found in this module. It will need GeneFrequencies
 *             read in using the read_GeneFrequency21_file function in
 *             the genefrequency module.  Again if Wise2 has been set up
 *             correctly, read_GeneFrequency21_file("human.gf") should work
 *
 *             To again a valid algorithm type use gwrap_alg_type_from_string
 *             found in this module. gwrap_alg_type_from_string("623") would
 *             be a good choice
 *
 *
 *             This function basically makes a threestatemodel (standard HMM) from
 *             the protein and the comparison matrix with the *scary* assumption that
 *             the comparison matrix is in half bit form. It then calls 
 *              /AlnBlock_from_TSM_genewise_wrap to do the nasty stuff. 
 *
 *
 * Arg:             protein [UNKN ] protein sequence used in the comparison [Protein *]
 * Arg:                  pg [READ ] Potential gene - could be NULL - if rough exon positions are known  [NullString]
 * Arg:           is_global [UNKN ] has now become flag for local/global/end-biased switch [NullString]
 * Arg:                 dna [UNKN ] genomic DNA sequence used  [Genomic *]
 * Arg:                comp [UNKN ] protein comparison matrix *in half bits* [CompMat *]
 * Arg:                 gap [UNKN ] gap penalty (negative) [int]
 * Arg:                 ext [UNKN ] extension penalty (negative) [int]
 * Arg:               gpara [UNKN ] Gene parameters. [GeneParameter21 *]
 * Arg:                 rmd [UNKN ] models to be compared to [RandomModelDNA *]
 * Arg:          intergenic [UNKN ] model of random dna between genes [RandomModelDNA *]
 * Arg:                 alg [UNKN ] algorithm type [int]
 * Arg:             use_syn [UNKN ] Undocumented argument [boolean]
 * Arg:                  rm [UNKN ] Undocumented argument [RandomModel *]
 * Arg:                allN [UNKN ] Undocumented argument [Probability]
 * Arg:        startendmode [UNKN ] Undocumented argument [TSM_StartEndMode]
 * Arg:                dpri [UNKN ] Undocumented argument [DPRunImpl *]
 * Arg:                 pal [WRITE] Raw alginment to be saved if non-NULL [PackAln **]
 * Arg:                 gwp [UNKN ] Undocumented argument [GeneWiseRunPara *]
 *
 * Return [UNKN ]  Undocumented return value [AlnBlock *]
 *
 */
AlnBlock * Wise2_AlnBlock_from_protein_genewise_wrap(Protein * protein,Genomic * dna,CompMat * comp,int gap,int ext,GeneParameter21 * gpara,RandomModelDNA * rmd,RandomModelDNA * intergenic,int alg,boolean use_syn,RandomModel * rm,Probability allN,TSM_StartEndMode startendmode,DPRunImpl * dpri,PackAln ** pal,GeneWiseRunPara * gwp);
#define AlnBlock_from_protein_genewise_wrap Wise2_AlnBlock_from_protein_genewise_wrap


/* Function:  AlnBlock_from_TSM_genewise_wrap(tsm,gen,gpara,rmd,intergenic,use_syn,alg,allN,flat_insert,dpri,palpoi,dpenv)
 *
 * Descrip:    A function which aligns a protein HMM (as found
 *             in my threestatemodel structure) to a genomic DNA 
 *             sequence. 
 *
 *             At the moment you are unlikely to be reading in the
 *             HMM structure yourself, so this is not something
 *             you will be doing.
 *
 *             The core algorithms for each method are found in
 *             genewise21/geneloop21 etc files. 
 *
 *
 *
 * Arg:                tsm [UNKN ] protein TSM to be used in the comparison [ThreeStateModel *]
 * Arg:                gen [UNKN ] genomic DNA sequence used  [Genomic *]
 * Arg:              gpara [UNKN ] Gene parameters. [GeneParameter21 *]
 * Arg:                rmd [UNKN ] models to be compared to [RandomModelDNA *]
 * Arg:         intergenic [UNKN ] model of random dna between genes [RandomModelDNA *]
 * Arg:            use_syn [UNKN ] use a synchronous null model [boolean]
 * Arg:                alg [UNKN ] algorithm type [int]
 * Arg:               allN [UNKN ] Undocumented argument [Probability]
 * Arg:        flat_insert [UNKN ] Undocumented argument [boolean]
 * Arg:               dpri [UNKN ] Undocumented argument [DPRunImpl *]
 * Arg:             palpoi [WRITE] Raw alginment to be saved if non-NULL [PackAln **]
 * Arg:              dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [AlnBlock *]
 *
 */
AlnBlock * Wise2_AlnBlock_from_TSM_genewise_wrap(ThreeStateModel * tsm,Genomic * gen,GeneParameter21 * gpara,RandomModelDNA * rmd,RandomModelDNA * intergenic,boolean use_syn,int alg,Probability allN,boolean flat_insert,DPRunImpl * dpri,PackAln ** palpoi,DPEnvelope * dpenv);
#define AlnBlock_from_TSM_genewise_wrap Wise2_AlnBlock_from_TSM_genewise_wrap


/* Function:  Hscore_from_TSM_genewise(tdb,gdb,gpara,rmd,intergenic,use_syn,alg,bits_cutoff,allN,report_level,die_on_error,flat_insert,dbsi)
 *
 * Descrip:    Runs a database search of the genewise algorithm. 
 *
 *             This makes a high score object which you can then use 
 *             to retrieve enteries as well as print out the top score (!)
 *
 *
 * Arg:                 tdb [READ ] a database of profileHMMs  [ThreeStateDB *]
 * Arg:                 gdb [READ ] a database of genomic sequence [GenomicDB *]
 * Arg:               gpara [READ ] geneparameters [GeneParameter21 *]
 * Arg:                 rmd [READ ] random model to be compared with in non syn mode [RandomModelDNA *]
 * Arg:          intergenic [READ ] random model of intergenic DNA (usually the same as rmd) [RandomModelDNA *]
 * Arg:             use_syn [UNKN ] use synchronous random model [boolean]
 * Arg:                 alg [UNKN ] algorithm type [int]
 * Arg:         bits_cutoff [UNKN ] cutoff in bits of the scores to store [double]
 * Arg:                allN [UNKN ] Undocumented argument [Probability]
 * Arg:        report_level [UNKN ] stagger rate of reporting progress on stderr  -1 means never [int]
 * Arg:        die_on_error [UNKN ] if true, exits on error (not used at the moment) [boolean]
 * Arg:         flat_insert [UNKN ] Undocumented argument [boolean]
 * Arg:                dbsi [UNKN ] Undocumented argument [DBSearchImpl *]
 *
 * Return [UNKN ]  a new Hscore object of the entire db search [Hscore *]
 *
 */
Hscore * Wise2_Hscore_from_TSM_genewise(ThreeStateDB * tdb,GenomicDB * gdb,GeneParameter21 * gpara,RandomModelDNA * rmd,RandomModelDNA * intergenic,boolean use_syn,int alg,double bits_cutoff,Probability allN,int report_level,boolean die_on_error,boolean flat_insert,DBSearchImpl * dbsi);
#define Hscore_from_TSM_genewise Wise2_Hscore_from_TSM_genewise


/* Function:  gwrap_alg_type_from_string(str)
 *
 * Descrip:    Gives you the integer interpretation from
 *             the string, which is one of
 *             2193 2193L, 623, 623L, 421, 2193LINK
 *
 *             This integer can then be passed into routines
 *             like AlnBlock_from_protein_genewise_wrap
 *
 *
 *
 * Arg:        str [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_gwrap_alg_type_from_string(char * str);
#define gwrap_alg_type_from_string Wise2_gwrap_alg_type_from_string


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
cDNAParserScore * Wise2_cDNAParserScore_from_GeneParser21Score(GeneParser21Score * gps);
#define cDNAParserScore_from_GeneParser21Score Wise2_cDNAParserScore_from_GeneParser21Score

#ifdef _cplusplus
}
#endif

#endif
