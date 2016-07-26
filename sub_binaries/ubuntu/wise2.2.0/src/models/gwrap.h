#ifndef DYNAMITEgwrapHEADERFILE
#define DYNAMITEgwrapHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "dyna.h"
#include "genewise6.h"
#include "geneloop6.h"
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

#include "seqhit.h"
#include "estwise3.h"
#include "estloop3.h"


#define PotentialGeneListLISTLENGTH 64
#define PotentialGeneLISTLENGTH 64
#define PotentialTranscriptLISTLENGTH 64

enum GWWRAP_ALG_TYPE {
  GWWRAP_2193 = 12,
  GWWRAP_2193L,
  GWWRAP_2193I,
  GWWRAP_623,
  GWWRAP_623L,
  GWWRAP_421,
  GWWRAP_6LITE,
  GWWRAP_333,
  GWWRAP_333L
};


/* Object PotentialExon
 *
 * Descrip: Data structure to represent exon information
 *        to be passed into dpenvelopes etc
 *
 *
 */
struct Wise2_PotentialExon {  
    int dynamite_hard_link;  
    int tstart;  
    int tend;    
    int qstart;  
    int qend;    
    } ;  
/* PotentialExon defined */ 
#ifndef DYNAMITE_DEFINED_PotentialExon
typedef struct Wise2_PotentialExon Wise2_PotentialExon;
#define PotentialExon Wise2_PotentialExon
#define DYNAMITE_DEFINED_PotentialExon
#endif


struct Wise2_PotentialTranscript {  
    int dynamite_hard_link;  
    PotentialExon ** pex;    
    int len;/* len for above pex  */ 
    int maxlen; /* maxlen for above pex */ 
    } ;  
/* PotentialTranscript defined */ 
#ifndef DYNAMITE_DEFINED_PotentialTranscript
typedef struct Wise2_PotentialTranscript Wise2_PotentialTranscript;
#define PotentialTranscript Wise2_PotentialTranscript
#define DYNAMITE_DEFINED_PotentialTranscript
#endif


/* Object PotentialGene
 *
 * Descrip: This data structure hopefully stores
 *        the necessary information for finding a
 *        gene via a gene wise type algorithm
 *
 *
 */
struct Wise2_PotentialGene {  
    int dynamite_hard_link;  
    int guess_start;     
    int guess_end;   
    PotentialTranscript ** pet;  
    int len;/* len for above pet  */ 
    int maxlen; /* maxlen for above pet */ 
    boolean is_global;   
    char * name;     
    Protein * homolog;  /*  one of these two will be used */ 
    ThreeStateModel * tsm;   
    AlnBlock * alb; /*  if we want to save the alignments. */ 
    double bitscore;     
    int slop_query;  
    int slop_target;     
    int query_length;    
    } ;  
/* PotentialGene defined */ 
#ifndef DYNAMITE_DEFINED_PotentialGene
typedef struct Wise2_PotentialGene Wise2_PotentialGene;
#define PotentialGene Wise2_PotentialGene
#define DYNAMITE_DEFINED_PotentialGene
#endif


/* Object PotentialGeneList
 *
 * Descrip: a list of potential genes.
 *
 *        Made either - 
 *             externally or
 *
 *             from a MSP crunch datastructure
 *
 *
 */
struct Wise2_PotentialGeneList {  
    int dynamite_hard_link;  
    PotentialGene ** pg;     
    int len;/* len for above pg  */ 
    int maxlen; /* maxlen for above pg */ 
    } ;  
/* PotentialGeneList defined */ 
#ifndef DYNAMITE_DEFINED_PotentialGeneList
typedef struct Wise2_PotentialGeneList Wise2_PotentialGeneList;
#define PotentialGeneList Wise2_PotentialGeneList
#define DYNAMITE_DEFINED_PotentialGeneList
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  read_PotentialGene_file(filename)
 *
 * Descrip:    reads in potential gene format having open the
 *             file
 *
 *
 * Arg:        filename [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [PotentialGene *]
 *
 */
PotentialGene * Wise2_read_PotentialGene_file(char * filename);
#define read_PotentialGene_file Wise2_read_PotentialGene_file


/* Function:  read_PotentialGene(*line,ifp)
 *
 * Descrip:    Reads in a potential gene format
 *
 *
 * Arg:        *line [UNKN ] Undocumented argument [char]
 * Arg:          ifp [UNKN ] Undocumented argument [FILE *]
 *
 * Return [UNKN ]  Undocumented return value [PotentialGene *]
 *
 */
PotentialGene * Wise2_read_PotentialGene(char *line,FILE * ifp);
#define read_PotentialGene Wise2_read_PotentialGene


/* Function:  read_PotentialTranscript(line,ifp)
 *
 * Descrip:    read in format ->
 *                pgene
 *                ptrans
 *                pexon start end qstart qend
 *                endptrans
 *                endpgene
 *
 *
 * Arg:        line [UNKN ] Undocumented argument [char *]
 * Arg:         ifp [UNKN ] Undocumented argument [FILE *]
 *
 * Return [UNKN ]  Undocumented return value [PotentialTranscript *]
 *
 */
PotentialTranscript * Wise2_read_PotentialTranscript(char * line,FILE * ifp);
#define read_PotentialTranscript Wise2_read_PotentialTranscript


/* Function:  DPEnvelope_from_PotentialGene(pg)
 *
 * Descrip:    Makes a DPEnv structure from a potential gene
 *             with Exons in it
 *
 *             If there are no potential transcripts, returns NULL,
 *             which is what you probably want
 *
 *
 * Arg:        pg [UNKN ] Undocumented argument [PotentialGene *]
 *
 * Return [UNKN ]  Undocumented return value [DPEnvelope *]
 *
 */
DPEnvelope * Wise2_DPEnvelope_from_PotentialGene(PotentialGene * pg);
#define DPEnvelope_from_PotentialGene Wise2_DPEnvelope_from_PotentialGene


/* Function:  add_PotentialTranscript_to_DPEnvelope(dpenv,pet,pg)
 *
 * Descrip:    Adds the potential exons in the Potential transcript to
 *             the dpenvelope.
 *
 *
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:          pet [UNKN ] Undocumented argument [PotentialTranscript *]
 * Arg:           pg [UNKN ] Undocumented argument [PotentialGene *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_PotentialTranscript_to_DPEnvelope(DPEnvelope * dpenv,PotentialTranscript * pet,PotentialGene * pg);
#define add_PotentialTranscript_to_DPEnvelope Wise2_add_PotentialTranscript_to_DPEnvelope


/* Function:  PotentialGeneList_from_DnaSequenceHitList(dsl,window,wing_length,min_score)
 *
 * Descrip:    This makes a PotentialGeneList from a
 *             DnaSequenceHitList, which is a module which,
 *             for example, abstracts the MSP crunch output.
 *
 *             The three parameters are:
 *                 window - what window size to consider a potential gene in
 *                 wing_length - length of wing sequences to add onto the start/end points of a hit
 *                 min_score - minimum score to trigger a potential gene.
 *
 *             The potential genes are selected as follows:
 *
 *                 foreach window
 *                        Take the best scoring segment.
 *                        if( > min_score) 
 *                             if( there_is_a_segment which start/end + wing_length overlaps + the same name)
 *                                  extend that segment
 *                             else
 *                                  make a new potential gene
 *
 *
 *
 * Arg:                dsl [UNKN ] Undocumented argument [DnaSequenceHitList *]
 * Arg:             window [UNKN ] Undocumented argument [int]
 * Arg:        wing_length [UNKN ] Undocumented argument [int]
 * Arg:          min_score [UNKN ] Undocumented argument [double]
 *
 * Return [UNKN ]  Undocumented return value [PotentialGeneList *]
 *
 */
PotentialGeneList * Wise2_PotentialGeneList_from_DnaSequenceHitList(DnaSequenceHitList * dsl,int window,int wing_length,double min_score);
#define PotentialGeneList_from_DnaSequenceHitList Wise2_PotentialGeneList_from_DnaSequenceHitList


/* Function:  read_PotentialGeneList_pgfasta_file(filename)
 *
 * Descrip:    reads a file with lines like
 *
 *              >ROA1_HUMAN:1234:3445
 *              <sequence>
 *
 *             As potential gene with a guess start point at 1234 and end
 *             at 3445
 *
 *
 * Arg:        filename [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [PotentialGeneList *]
 *
 */
PotentialGeneList * Wise2_read_PotentialGeneList_pgfasta_file(char * filename);
#define read_PotentialGeneList_pgfasta_file Wise2_read_PotentialGeneList_pgfasta_file


/* Function:  read_PotentialGeneList_pgfasta(ifp)
 *
 * Descrip:    reads a file with lines like
 *
 *             >ROA1_HUMAN:1234:3445
 *             <sequence>
 *
 *             As potential gene with a guess start point at 1234 and end
 *             at 3445
 *
 *
 * Arg:        ifp [UNKN ] Undocumented argument [FILE *]
 *
 * Return [UNKN ]  Undocumented return value [PotentialGeneList *]
 *
 */
PotentialGeneList * Wise2_read_PotentialGeneList_pgfasta(FILE * ifp);
#define read_PotentialGeneList_pgfasta Wise2_read_PotentialGeneList_pgfasta


/* Function:  resolve_PotentialGenes_on_GenomicRegion(gr,pgl,alg_protein,alg_hmm,prot_thr,hmm_thr,comp,gap,ext,gpara,rmd,inter,fetch_from_pipe,should_free,make_name,bit_cut_off,dpri)
 *
 * Descrip:    Takes the potential gene list, (made for example from
 *             MSP crunch file through DnaSequenceHitList object),
 *             and calculates genes on it.
 *
 *             This is the core functionality to postwise.
 *
 *             There are two basic modes to this routine that probably
 *             should be split up. In the first mode, the Potential Gene List
 *             has actual proteins/or HMMs (tsm - threestatemodels) in 
 *             memory in the structure. This then loops through the list
 *             calling genewise functions (going through the wrap functions
 *             in this module).
 *
 *             The other mode has no proteins in memory but a way of fetching
 *             proteins from a pipe using a string passed in. The string looks
 *             like a sprintf string where the name of the protein is the
 *             only thing to be substituted. For example
 *                 
 *                efetch -f %s
 *
 *             or
 *
 *                getz -d '[swissprot-id:%s]'
 *
 *             would be sensible examples. Remember that the command should produce
 *             things on stdout, and that you should (obviously) make sure that the
 *             indexer uses the same name as the name in, for example, the msp crunch
 *             file.
 *
 *             If should_free is true then it frees the protein sequence and any
 *             alignment after the calculation. Otherwise it stores both of these
 *             in the potentialgene structure.
 *
 *             For the business end of the algorithm, this function uses the
 *             /AlnBlock_from_protein_genewise_wrap and /AlnBlock_from_TSM_genewise_wrap
 *             functions in this module. 
 *                 
 *
 *
 * Arg:                     gr [UNKN ] genomic region to make genes on [GenomicRegion *]
 * Arg:                    pgl [UNKN ] potential gene list of homologs and start/end points [PotentialGeneList *]
 * Arg:            alg_protein [UNKN ] algorithm to use with protein comparisons [int]
 * Arg:                alg_hmm [UNKN ] algorithm to use with HMM comparisons [int]
 * Arg:               prot_thr [UNKN ] Threshold under which genes are not predicted [double]
 * Arg:                hmm_thr [UNKN ] Undocumented argument [double]
 * Arg:                   comp [UNKN ] Comparison matrix for protein sequences [CompMat *]
 * Arg:                    gap [UNKN ] Gap cost for protein sequences [int]
 * Arg:                    ext [UNKN ] Extension cost for protein sequences [int]
 * Arg:                  gpara [UNKN ] Gene parameters [GeneParameter21 *]
 * Arg:                    rmd [UNKN ] Random Model to compare against [RandomModelDNA *]
 * Arg:                  inter [UNKN ] Random Model for intergenic regions [RandomModelDNA *]
 * Arg:        fetch_from_pipe [UNKN ] For potential genes with just a name, a pipe to get sequences from [char *]
 * Arg:            should_free [UNKN ] Free the protein sequences after resolving it [boolean]
 * Arg:              make_name [FUNCP] a pointer to a function which makes the name of the gene [char *(*make_name]
 * Arg:            bit_cut_off [UNKN ] Undocumented argument [double]
 * Arg:                   dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_resolve_PotentialGenes_on_GenomicRegion(GenomicRegion * gr,PotentialGeneList * pgl,int alg_protein,int alg_hmm,double prot_thr,double hmm_thr,CompMat * comp,int gap,int ext,GeneParameter21 * gpara,RandomModelDNA * rmd,RandomModelDNA * inter,char * fetch_from_pipe,boolean should_free,char *(*make_name)(Wise2_Genomic * gen,char *,int,Wise2_Gene *),double bit_cut_off,DPRunImpl * dpri);
#define resolve_PotentialGenes_on_GenomicRegion Wise2_resolve_PotentialGenes_on_GenomicRegion


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


/* Function:  AlnBlock_from_protein_genewise_wrap(protein,dna,comp,gap,ext,gpara,rmd,intergenic,alg,is_global,use_syn,rm,allN,pg,dpri,pal)
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
 * Arg:           protein [UNKN ] protein sequence used in the comparison [Protein *]
 * Arg:               dna [UNKN ] genomic DNA sequence used  [Genomic *]
 * Arg:              comp [UNKN ] protein comparison matrix *in half bits* [CompMat *]
 * Arg:               gap [UNKN ] gap penalty (negative) [int]
 * Arg:               ext [UNKN ] extension penalty (negative) [int]
 * Arg:             gpara [UNKN ] Gene parameters. [GeneParameter21 *]
 * Arg:               rmd [UNKN ] models to be compared to [RandomModelDNA *]
 * Arg:        intergenic [UNKN ] model of random dna between genes [RandomModelDNA *]
 * Arg:               alg [UNKN ] algorithm type [int]
 * Arg:         is_global [UNKN ] has now become flag for local/global/end-biased switch [int]
 * Arg:           use_syn [UNKN ] Undocumented argument [boolean]
 * Arg:                rm [UNKN ] Undocumented argument [RandomModel *]
 * Arg:              allN [UNKN ] Undocumented argument [Probability]
 * Arg:                pg [READ ] Potential gene - could be NULL - if rough exon positions are known  [PotentialGene *]
 * Arg:              dpri [UNKN ] Undocumented argument [DPRunImpl *]
 * Arg:               pal [WRITE] Raw alginment to be saved if non-NULL [PackAln **]
 *
 * Return [UNKN ]  Undocumented return value [AlnBlock *]
 *
 */
AlnBlock * Wise2_AlnBlock_from_protein_genewise_wrap(Protein * protein,Genomic * dna,CompMat * comp,int gap,int ext,GeneParameter21 * gpara,RandomModelDNA * rmd,RandomModelDNA * intergenic,int alg,int is_global,boolean use_syn,RandomModel * rm,Probability allN,PotentialGene * pg,DPRunImpl * dpri,PackAln ** pal);
#define AlnBlock_from_protein_genewise_wrap Wise2_AlnBlock_from_protein_genewise_wrap


/* Function:  AlnBlock_from_TSM_genewise_wrap(tsm,gen,gpara,rmd,intergenic,use_syn,alg,allN,flat_insert,pg,dpri,palpoi)
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
 * Arg:                 pg [READ ] Potential gene - could be NULL - if rough exon positions are known [PotentialGene *]
 * Arg:               dpri [UNKN ] Undocumented argument [DPRunImpl *]
 * Arg:             palpoi [WRITE] Raw alginment to be saved if non-NULL [PackAln **]
 *
 * Return [UNKN ]  Undocumented return value [AlnBlock *]
 *
 */
AlnBlock * Wise2_AlnBlock_from_TSM_genewise_wrap(ThreeStateModel * tsm,Genomic * gen,GeneParameter21 * gpara,RandomModelDNA * rmd,RandomModelDNA * intergenic,boolean use_syn,int alg,Probability allN,boolean flat_insert,PotentialGene * pg,DPRunImpl * dpri,PackAln ** palpoi);
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


/* Function:  hard_link_PotentialExon(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [PotentialExon *]
 *
 * Return [UNKN ]  Undocumented return value [PotentialExon *]
 *
 */
PotentialExon * Wise2_hard_link_PotentialExon(PotentialExon * obj);
#define hard_link_PotentialExon Wise2_hard_link_PotentialExon


/* Function:  PotentialExon_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [PotentialExon *]
 *
 */
PotentialExon * Wise2_PotentialExon_alloc(void);
#define PotentialExon_alloc Wise2_PotentialExon_alloc


/* Function:  free_PotentialExon(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [PotentialExon *]
 *
 * Return [UNKN ]  Undocumented return value [PotentialExon *]
 *
 */
PotentialExon * Wise2_free_PotentialExon(PotentialExon * obj);
#define free_PotentialExon Wise2_free_PotentialExon


/* Function:  add_PotentialTranscript(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [PotentialTranscript *]
 * Arg:        add [OWNER] Object to add to the list [PotentialExon *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_PotentialTranscript(PotentialTranscript * obj,PotentialExon * add);
#define add_PotentialTranscript Wise2_add_PotentialTranscript


/* Function:  flush_PotentialTranscript(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [PotentialTranscript *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_PotentialTranscript(PotentialTranscript * obj);
#define flush_PotentialTranscript Wise2_flush_PotentialTranscript


/* Function:  PotentialTranscript_alloc_std(void)
 *
 * Descrip:    Equivalent to PotentialTranscript_alloc_len(PotentialTranscriptLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [PotentialTranscript *]
 *
 */
PotentialTranscript * Wise2_PotentialTranscript_alloc_std(void);
#define PotentialTranscript_alloc_std Wise2_PotentialTranscript_alloc_std


/* Function:  PotentialTranscript_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [PotentialTranscript *]
 *
 */
PotentialTranscript * Wise2_PotentialTranscript_alloc_len(int len);
#define PotentialTranscript_alloc_len Wise2_PotentialTranscript_alloc_len


/* Function:  hard_link_PotentialTranscript(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [PotentialTranscript *]
 *
 * Return [UNKN ]  Undocumented return value [PotentialTranscript *]
 *
 */
PotentialTranscript * Wise2_hard_link_PotentialTranscript(PotentialTranscript * obj);
#define hard_link_PotentialTranscript Wise2_hard_link_PotentialTranscript


/* Function:  PotentialTranscript_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [PotentialTranscript *]
 *
 */
PotentialTranscript * Wise2_PotentialTranscript_alloc(void);
#define PotentialTranscript_alloc Wise2_PotentialTranscript_alloc


/* Function:  free_PotentialTranscript(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [PotentialTranscript *]
 *
 * Return [UNKN ]  Undocumented return value [PotentialTranscript *]
 *
 */
PotentialTranscript * Wise2_free_PotentialTranscript(PotentialTranscript * obj);
#define free_PotentialTranscript Wise2_free_PotentialTranscript


/* Function:  add_PotentialGene(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [PotentialGene *]
 * Arg:        add [OWNER] Object to add to the list [PotentialTranscript *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_PotentialGene(PotentialGene * obj,PotentialTranscript * add);
#define add_PotentialGene Wise2_add_PotentialGene


/* Function:  flush_PotentialGene(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [PotentialGene *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_PotentialGene(PotentialGene * obj);
#define flush_PotentialGene Wise2_flush_PotentialGene


/* Function:  PotentialGene_alloc_std(void)
 *
 * Descrip:    Equivalent to PotentialGene_alloc_len(PotentialGeneLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [PotentialGene *]
 *
 */
PotentialGene * Wise2_PotentialGene_alloc_std(void);
#define PotentialGene_alloc_std Wise2_PotentialGene_alloc_std


/* Function:  PotentialGene_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [PotentialGene *]
 *
 */
PotentialGene * Wise2_PotentialGene_alloc_len(int len);
#define PotentialGene_alloc_len Wise2_PotentialGene_alloc_len


/* Function:  hard_link_PotentialGene(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [PotentialGene *]
 *
 * Return [UNKN ]  Undocumented return value [PotentialGene *]
 *
 */
PotentialGene * Wise2_hard_link_PotentialGene(PotentialGene * obj);
#define hard_link_PotentialGene Wise2_hard_link_PotentialGene


/* Function:  PotentialGene_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [PotentialGene *]
 *
 */
PotentialGene * Wise2_PotentialGene_alloc(void);
#define PotentialGene_alloc Wise2_PotentialGene_alloc


/* Function:  free_PotentialGene(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [PotentialGene *]
 *
 * Return [UNKN ]  Undocumented return value [PotentialGene *]
 *
 */
PotentialGene * Wise2_free_PotentialGene(PotentialGene * obj);
#define free_PotentialGene Wise2_free_PotentialGene


/* Function:  add_PotentialGeneList(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [PotentialGeneList *]
 * Arg:        add [OWNER] Object to add to the list [PotentialGene *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_PotentialGeneList(PotentialGeneList * obj,PotentialGene * add);
#define add_PotentialGeneList Wise2_add_PotentialGeneList


/* Function:  flush_PotentialGeneList(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [PotentialGeneList *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_PotentialGeneList(PotentialGeneList * obj);
#define flush_PotentialGeneList Wise2_flush_PotentialGeneList


/* Function:  PotentialGeneList_alloc_std(void)
 *
 * Descrip:    Equivalent to PotentialGeneList_alloc_len(PotentialGeneListLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [PotentialGeneList *]
 *
 */
PotentialGeneList * Wise2_PotentialGeneList_alloc_std(void);
#define PotentialGeneList_alloc_std Wise2_PotentialGeneList_alloc_std


/* Function:  PotentialGeneList_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [PotentialGeneList *]
 *
 */
PotentialGeneList * Wise2_PotentialGeneList_alloc_len(int len);
#define PotentialGeneList_alloc_len Wise2_PotentialGeneList_alloc_len


/* Function:  hard_link_PotentialGeneList(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [PotentialGeneList *]
 *
 * Return [UNKN ]  Undocumented return value [PotentialGeneList *]
 *
 */
PotentialGeneList * Wise2_hard_link_PotentialGeneList(PotentialGeneList * obj);
#define hard_link_PotentialGeneList Wise2_hard_link_PotentialGeneList


/* Function:  PotentialGeneList_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [PotentialGeneList *]
 *
 */
PotentialGeneList * Wise2_PotentialGeneList_alloc(void);
#define PotentialGeneList_alloc Wise2_PotentialGeneList_alloc


/* Function:  free_PotentialGeneList(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [PotentialGeneList *]
 *
 * Return [UNKN ]  Undocumented return value [PotentialGeneList *]
 *
 */
PotentialGeneList * Wise2_free_PotentialGeneList(PotentialGeneList * obj);
#define free_PotentialGeneList Wise2_free_PotentialGeneList


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
int Wise2_access_guess_end_PotentialGene(PotentialGene * obj);
#define access_guess_end_PotentialGene Wise2_access_guess_end_PotentialGene
boolean Wise2_replace_name_PotentialGene(PotentialGene * obj,char * name);
#define replace_name_PotentialGene Wise2_replace_name_PotentialGene
boolean Wise2_replace_homolog_PotentialGene(PotentialGene * obj,Protein * homolog);
#define replace_homolog_PotentialGene Wise2_replace_homolog_PotentialGene
boolean Wise2_replace_qend_PotentialExon(PotentialExon * obj,int qend);
#define replace_qend_PotentialExon Wise2_replace_qend_PotentialExon
Protein * Wise2_access_homolog_PotentialGene(PotentialGene * obj);
#define access_homolog_PotentialGene Wise2_access_homolog_PotentialGene
boolean Wise2_replace_alb_PotentialGene(PotentialGene * obj,AlnBlock * alb);
#define replace_alb_PotentialGene Wise2_replace_alb_PotentialGene
int Wise2_access_qstart_PotentialExon(PotentialExon * obj);
#define access_qstart_PotentialExon Wise2_access_qstart_PotentialExon
AlnBlock * Wise2_access_alb_PotentialGene(PotentialGene * obj);
#define access_alb_PotentialGene Wise2_access_alb_PotentialGene
PotentialGene * Wise2_access_pg_PotentialGeneList(PotentialGeneList * obj,int i);
#define access_pg_PotentialGeneList Wise2_access_pg_PotentialGeneList
boolean Wise2_replace_bitscore_PotentialGene(PotentialGene * obj,double bitscore);
#define replace_bitscore_PotentialGene Wise2_replace_bitscore_PotentialGene
boolean Wise2_replace_guess_start_PotentialGene(PotentialGene * obj,int guess_start);
#define replace_guess_start_PotentialGene Wise2_replace_guess_start_PotentialGene
double Wise2_access_bitscore_PotentialGene(PotentialGene * obj);
#define access_bitscore_PotentialGene Wise2_access_bitscore_PotentialGene
boolean Wise2_replace_guess_end_PotentialGene(PotentialGene * obj,int guess_end);
#define replace_guess_end_PotentialGene Wise2_replace_guess_end_PotentialGene
boolean Wise2_replace_slop_query_PotentialGene(PotentialGene * obj,int slop_query);
#define replace_slop_query_PotentialGene Wise2_replace_slop_query_PotentialGene
int Wise2_length_pet_PotentialGene(PotentialGene * obj);
#define length_pet_PotentialGene Wise2_length_pet_PotentialGene
int Wise2_access_slop_query_PotentialGene(PotentialGene * obj);
#define access_slop_query_PotentialGene Wise2_access_slop_query_PotentialGene
boolean Wise2_access_is_global_PotentialGene(PotentialGene * obj);
#define access_is_global_PotentialGene Wise2_access_is_global_PotentialGene
boolean Wise2_replace_slop_target_PotentialGene(PotentialGene * obj,int slop_target);
#define replace_slop_target_PotentialGene Wise2_replace_slop_target_PotentialGene
int Wise2_access_qend_PotentialExon(PotentialExon * obj);
#define access_qend_PotentialExon Wise2_access_qend_PotentialExon
boolean Wise2_replace_qstart_PotentialExon(PotentialExon * obj,int qstart);
#define replace_qstart_PotentialExon Wise2_replace_qstart_PotentialExon
int Wise2_access_slop_target_PotentialGene(PotentialGene * obj);
#define access_slop_target_PotentialGene Wise2_access_slop_target_PotentialGene
ThreeStateModel * Wise2_access_tsm_PotentialGene(PotentialGene * obj);
#define access_tsm_PotentialGene Wise2_access_tsm_PotentialGene
boolean Wise2_replace_query_length_PotentialGene(PotentialGene * obj,int query_length);
#define replace_query_length_PotentialGene Wise2_replace_query_length_PotentialGene
int Wise2_length_pg_PotentialGeneList(PotentialGeneList * obj);
#define length_pg_PotentialGeneList Wise2_length_pg_PotentialGeneList
int Wise2_access_query_length_PotentialGene(PotentialGene * obj);
#define access_query_length_PotentialGene Wise2_access_query_length_PotentialGene
PotentialTranscript * Wise2_access_pet_PotentialGene(PotentialGene * obj,int i);
#define access_pet_PotentialGene Wise2_access_pet_PotentialGene
PotentialExon * Wise2_access_pex_PotentialTranscript(PotentialTranscript * obj,int i);
#define access_pex_PotentialTranscript Wise2_access_pex_PotentialTranscript
char * Wise2_access_name_PotentialGene(PotentialGene * obj);
#define access_name_PotentialGene Wise2_access_name_PotentialGene
int Wise2_length_pex_PotentialTranscript(PotentialTranscript * obj);
#define length_pex_PotentialTranscript Wise2_length_pex_PotentialTranscript
void Wise2_show_PotentialGeneList(PotentialGeneList *pgl,FILE * ofp);
#define show_PotentialGeneList Wise2_show_PotentialGeneList
boolean Wise2_replace_tstart_PotentialExon(PotentialExon * obj,int tstart);
#define replace_tstart_PotentialExon Wise2_replace_tstart_PotentialExon
boolean Wise2_replace_is_global_PotentialGene(PotentialGene * obj,boolean is_global);
#define replace_is_global_PotentialGene Wise2_replace_is_global_PotentialGene
int Wise2_access_tstart_PotentialExon(PotentialExon * obj);
#define access_tstart_PotentialExon Wise2_access_tstart_PotentialExon
int Wise2_access_guess_start_PotentialGene(PotentialGene * obj);
#define access_guess_start_PotentialGene Wise2_access_guess_start_PotentialGene
boolean Wise2_replace_tend_PotentialExon(PotentialExon * obj,int tend);
#define replace_tend_PotentialExon Wise2_replace_tend_PotentialExon
boolean Wise2_replace_tsm_PotentialGene(PotentialGene * obj,ThreeStateModel * tsm);
#define replace_tsm_PotentialGene Wise2_replace_tsm_PotentialGene
int Wise2_access_tend_PotentialExon(PotentialExon * obj);
#define access_tend_PotentialExon Wise2_access_tend_PotentialExon
void Wise2_sort_PotentialExons_by_start(PotentialTranscript * pet);
#define sort_PotentialExons_by_start Wise2_sort_PotentialExons_by_start
int Wise2_compare_PotentialExons_start(PotentialExon * one,PotentialExon * two);
#define compare_PotentialExons_start Wise2_compare_PotentialExons_start
void Wise2_invsort_PotentialExons_by_start(PotentialTranscript * pet);
#define invsort_PotentialExons_by_start Wise2_invsort_PotentialExons_by_start
int Wise2_invcompare_PotentialExons_start(PotentialExon * one,PotentialExon * two);
#define invcompare_PotentialExons_start Wise2_invcompare_PotentialExons_start
cDNAParserScore * Wise2_cDNAParserScore_from_GeneParser21Score(GeneParser21Score * gps);
#define cDNAParserScore_from_GeneParser21Score Wise2_cDNAParserScore_from_GeneParser21Score
void Wise2_swap_PotentialTranscript(PotentialExon ** list,int i,int j) ;
#define swap_PotentialTranscript Wise2_swap_PotentialTranscript
void Wise2_qsort_PotentialTranscript(PotentialExon ** list,int left,int right,int (*comp)(PotentialExon * ,PotentialExon * ));
#define qsort_PotentialTranscript Wise2_qsort_PotentialTranscript
void Wise2_sort_PotentialTranscript(PotentialTranscript * obj,int (*comp)(PotentialExon *, PotentialExon *));
#define sort_PotentialTranscript Wise2_sort_PotentialTranscript
boolean Wise2_expand_PotentialTranscript(PotentialTranscript * obj,int len);
#define expand_PotentialTranscript Wise2_expand_PotentialTranscript
void Wise2_swap_PotentialGene(PotentialTranscript ** list,int i,int j) ;
#define swap_PotentialGene Wise2_swap_PotentialGene
void Wise2_qsort_PotentialGene(PotentialTranscript ** list,int left,int right,int (*comp)(PotentialTranscript * ,PotentialTranscript * ));
#define qsort_PotentialGene Wise2_qsort_PotentialGene
void Wise2_sort_PotentialGene(PotentialGene * obj,int (*comp)(PotentialTranscript *, PotentialTranscript *));
#define sort_PotentialGene Wise2_sort_PotentialGene
boolean Wise2_expand_PotentialGene(PotentialGene * obj,int len);
#define expand_PotentialGene Wise2_expand_PotentialGene
void Wise2_swap_PotentialGeneList(PotentialGene ** list,int i,int j) ;
#define swap_PotentialGeneList Wise2_swap_PotentialGeneList
void Wise2_qsort_PotentialGeneList(PotentialGene ** list,int left,int right,int (*comp)(PotentialGene * ,PotentialGene * ));
#define qsort_PotentialGeneList Wise2_qsort_PotentialGeneList
void Wise2_sort_PotentialGeneList(PotentialGeneList * obj,int (*comp)(PotentialGene *, PotentialGene *));
#define sort_PotentialGeneList Wise2_sort_PotentialGeneList
boolean Wise2_expand_PotentialGeneList(PotentialGeneList * obj,int len);
#define expand_PotentialGeneList Wise2_expand_PotentialGeneList

#ifdef _cplusplus
}
#endif

#endif
