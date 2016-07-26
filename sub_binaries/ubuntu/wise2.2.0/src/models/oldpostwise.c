#include "dyna.h"
#include "version.h" 
#include "genedisplay.h"

char * program_name = "postwise";

/*
 * program specific includes
 */
#include "gwrap.h"


/*
 * program specific variables
 */

char * pgl_file         = NULL; 
PotentialGeneList * pgl = NULL;

char * exblx_file       = NULL;
DnaSequenceHitList * dsl= NULL;

boolean use_msp_crunch  = TRUE;

char * fetch_string     = "efetch swir:%s";

char * dna_file         = NULL;
Genomic * gen           = NULL;
GenomicRegion * gr      = NULL;

char * gene_file        = NULL;
GeneFrequency21 * gf    = NULL;

char * matrix_file      = NULL;
CompMat * mat           = NULL;

char * gap_str          = "12";
int gap                 = 12;

char * ext_str          = "2";
int ext                 = 2;

char * codon_file       = NULL;
CodonTable * ct         = NULL;

char * output_file      = "-";
FILE * ofp              = NULL;

char * window_str       = "500";
int window              = 500;

char * wing_str         = "4000";
int wing                = 4000;

char * min_score_str    = "100.0";
double min_score        = 100.0;

char  * cut_off_str     = "0.0";
double cut_off          = 0.0;

char * palg_str         = "623";
int palg                = GWWRAP_623;

RandomModelDNA * rmd    = NULL;

boolean show_alignments  = FALSE;
boolean show_translation = FALSE;
boolean show_cdna        = FALSE;
boolean show_ace         = FALSE;
boolean show_gff         = FALSE;
boolean show_gene_str    = FALSE;
boolean show_pretty_gene = FALSE;
boolean show_raw_gene    = FALSE;

char * divide_str         = "//";
Probability rnd_loop      = 0.999;
Probability cds_loop      = 0.99;
Probability rnd_to_model  = (1 - 0.999) / 3;
Probability link_loop     = 0.98;
Probability link_to_model = (1- 0.98) / 3;

Probability subs_error    = 0.00001;
Probability indel_error   = 0.00001;

boolean model_codon       = FALSE;
boolean model_splice      = TRUE;
boolean tie_intron        = TRUE;

char * kbyte_str         = NULL;
int kbyte                = 10000; /* will be reset in build_defaults */

char * make_ace_gene_name(Genomic * gen,char * protein_name,int number_in_gene,Gene * gene)
{
  char buffer[128];
  char * alpha = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";
  int firstno,secondno;

  firstno = (int)floor(number_in_gene/26);
  number_in_gene -= firstno*26;

  sprintf(buffer,"%s.gw.%s.%c%c",gen->baseseq->name,protein_name,alpha[firstno],alpha[number_in_gene] );
  return stringalloc(buffer);
  
}

  
boolean show_output(void)
{
  int i;
  Protein * p;
  cDNA * c;
  Genomic * gentemp;
  Genomic * g;


  if( show_alignments == TRUE ) {
    for(i=0;i<pgl->len;i++) {
      if( pgl->pg[i]->homolog == NULL || pgl->pg[i]->alb == NULL ) {
	fprintf(ofp,"// Sorry - no prediction made, probably due to an error\n");
	continue;
      }

      if( pgl->pg[i]->bitscore < cut_off ) {
	fprintf(ofp,"// Bits score %f below cut off %f\n", pgl->pg[i]->bitscore,cut_off);
	continue;
      }

      gentemp = truncate_Genomic(gr->genomic,pgl->pg[i]->guess_start,pgl->pg[i]->guess_end);
      if( gentemp == NULL ) {
	warn("Could not make truncation in output display. Bad....");
	continue;
      }
      protgene_ascii_display(pgl->pg[i]->alb,pgl->pg[i]->homolog->baseseq->seq,pgl->pg[i]->homolog->baseseq->name,pgl->pg[i]->homolog->baseseq->offset,gentemp,ct,15,60,TRUE,ofp);
      free_Genomic(gentemp);
    }
    fprintf(ofp,"%s\n",divide_str);
  }

  if( show_ace == TRUE ) {
    show_ace_GenomicRegion(gr,gen->baseseq->name,ofp);
    fprintf(ofp,"%s\n",divide_str);
  }

  if( show_gff == TRUE ) {
    show_GFF_GenomicRegion(gr,gen->baseseq->name,"GeneWise",ofp);
    fprintf(ofp,"%s\n",divide_str);
  }

  if( show_translation == TRUE ) {
    for(i=0;i<gr->len;i++) {
      p = get_Protein_from_Translation(gr->gene[i]->transcript[0]->translation[0],ct);
      write_fasta_Sequence(p->baseseq,ofp);
    }
    fprintf(ofp,"%s\n",divide_str);
  }

  if( show_cdna == TRUE ) {
    for(i=0;i<gr->len;i++) {
      c = get_cDNA_from_Transcript(gr->gene[i]->transcript[0]);
      write_fasta_Sequence(c->baseseq,ofp);
    }
    fprintf(ofp,"%s\n",divide_str);
  }


  if( show_pretty_gene == TRUE ) {
    show_pretty_GenomicRegion(gr,ofp);
    fprintf(ofp,"%s\n",divide_str);
  }

  if( show_raw_gene == TRUE ) {
    show_GenomicRegion(gr,ofp);
    fprintf(ofp,"%s\n",divide_str);
  }
  
}


boolean build_gene_prediction(void)
{
  int count;
  GeneParameter21 * gp;

  /*  gp = GeneParameter21_wrap(gf,0.00001,0.00001,rmd,TRUE,TRUE,ct,0.1,0.1,0.1,0.1); */
  change_max_BaseMatrix_kbytes(kbyte);

  log_full_error(INFO,0,"I have %d Potential Genes to resolve",pgl->len);

  gp = GeneParameter21_wrap(gf,subs_error,indel_error,rmd,model_codon,model_splice,tie_intron,ct,rnd_loop,cds_loop,rnd_to_model,link_loop,link_to_model);

  count = resolve_PotentialGenes_on_GenomicRegion(gr,pgl,palg,GWWRAP_623L,10,10,mat,-gap,-ext,gp,rmd,rmd,fetch_string,show_alignments == TRUE ? FALSE : TRUE,make_ace_gene_name,cut_off);


  free_GeneParameter21(gp);

  return TRUE;
}

boolean build_objects(void)
{
  boolean ret = TRUE;
  FILE * ifp;


  if( is_integer_string(wing_str,&wing) == FALSE ) {
    warn("Wing amount [%s] is not an integer\n",wing_str);
    ret = FALSE;
  }

  if( is_integer_string(window_str,&window) == FALSE ) {
    warn("Window level [%s] is not an integer\n",window_str);
    ret = FALSE;
  }

  if( is_double_string(min_score_str,&min_score) == FALSE ) {
    warn("Minimum score [%s] is not a number\n",min_score_str);
    ret = FALSE;
  }

  if( is_double_string(cut_off_str,&cut_off) == FALSE ) {
    warn("Cut off string [%s] is not a number\n",cut_off_str);
    ret = FALSE;
  }
  
  if( use_msp_crunch == TRUE ) {
    ifp = openfile(exblx_file,"r");
    if( ifp == NULL ) {
      warn("Could not open %s as a exblx file",exblx_file);
      ret = FALSE;
    } else {
      dsl = read_MSPcrunch_DnaSequenceHitList(ifp);
      if( dsl == NULL ) {
	warn("Could not read MSP crunch data in %s",exblx_file);
	ret = FALSE;
      } else {
	/*fprintf(stderr,"Passing in %d and %d\n",window,wing); */
	/*show_DnaSequenceHitList(dsl,stdout);*/ 
	pgl = PotentialGeneList_from_DnaSequenceHitList(dsl,window,wing,min_score);
	if( pgl == NULL ) {
	  warn("Could not convert MSP crunch data in a potential gene list");
	  ret = FALSE;
	}
	/*show_PotentialGeneList(pgl,stdout);*/
      }
      
    } 
  } else {
  
    if( (pgl = read_PotentialGeneList_pgfasta_file(pgl_file)) == NULL ) {
    warn("Could not read Potential Genes in %s",pgl_file);
    ret = FALSE;
    }
  }

  if( pgl->len == 0 ) {
    warn("You have 0 Potential Genes. Considering this a bad thing!");
    ret = FALSE;
  }

  if( (gen = read_fasta_file_Genomic(dna_file)) == NULL) {
    ret = FALSE;
    warn("Could not read dna sequence file in %s",gf);
  } else {
    gr = new_GenomicRegion(gen);
  }

  if( (gf = read_GeneFrequency21_file(gene_file)) == NULL) {
    ret = FALSE;
    warn("Could not read a GeneFrequency file in %s",gene_file);
  }

  if( kbyte_str != NULL ) {
    if( is_integer_string(kbyte_str,&kbyte) == FALSE ) {
      warn("Could not get maximum kbyte number %s",kbyte_str);
      ret = FALSE;
    }
  }

  if( (palg = gwrap_alg_type_from_string(palg_str)) == -1 ) {
    warn("Cannot use %s as a valid algorithm type",palg_str);
    ret = FALSE;
  }

  if( (mat = read_Blast_file_CompMat(matrix_file)) == NULL) {
    ret = FALSE;
    warn("Could not read Comparison matrix file in %s",matrix_file);
  }

  if( (ct = read_CodonTable_file(codon_file)) == NULL) {
    ret = FALSE;
    warn("Could not read codon table file in %s",codon_file);
  }

  if( is_integer_string(gap_str,&gap) == FALSE ) {
    warn("Gap [%s] is not an integer\n",gap_str);
    ret = FALSE;
  }

  if( is_integer_string(ext_str,&ext) == FALSE ) {
    warn("Gap [%s] is not an integer\n",ext_str);
    ret = FALSE;
  }


  if( (ofp = openfile(output_file,"W")) ==  NULL) {
    warn("Could not open %s as an output file",output_file);
    ret = FALSE;
  }

  rmd = RandomModelDNA_std();
  return ret;
}
    
void build_defaults(void)
{
  gene_file = "human.gf";
  matrix_file = "blosum62.bla";
  codon_file = "codon.table";
  output_file = "-";

  kbyte = get_max_BaseMatrix_kbytes();
}

void show_short_help(void)
{
  fprintf(stdout,"%s (%s)\n",program_name,VERSION_NUMBER);
  fprintf(stdout,"This program is freely distributed under a GPL. See -version for more info\n");
  fprintf(stdout,"Copyright (c) GRL limited: portions of the code are from separate copyright\n\n");
  fprintf(stdout,"postwise <dna-file-fasta> <MSP-crunch -x output> '-' means stdin\n");
  fprintf(stdout," Options. In any order.\n");
  fprintf(stdout," Region selection [-window,-wing,-fetch,-min,-cut]\n");
  fprintf(stdout," Protein  [-gap,-ext,-matrix]\n");
  fprintf(stdout," Model    [-codon,-gene]\n");
  fprintf(stdout," Alg      [-palg,-kbyte]\n");
  fprintf(stdout," Output   [-pretty,-ace,-gff,-genes,-trans,-divide,-cdna]\n");
  fprintf(stdout," Standard [-help,-version,-silent,-quiet,-errorlog]\n");
  fprintf(stdout,"\nFor more help go %s -help.\n",program_name);
  fprintf(stdout,"\nSee WWW help at http://www.sanger.ac.uk/Software/Wise2/postwise.shtml\n");
  exit(63);   
}

void show_help(char * help_arg)
{
  fprintf(stdout,"%s (%s)\n",program_name,VERSION_NUMBER);
  fprintf(stdout,"postwise <dna-file-fasta> <msp-crunch -d output>\nOutput on stdout\n");
  /* program specific help */
  fprintf(stdout,"  -gap    [%3d]  gap penalty\n",gap);
  fprintf(stdout,"  -ext    [%3d]  extension penalty\n",gap);
  fprintf(stdout,"  -matrix [%s]  Comparison matrix\n",matrix_file);
  fprintf(stdout,"  -codon  [%s]  Codon file\n",codon_file);
  fprintf(stdout,"  -gene   [%s]  Gene parameter file\n",gene_file);
  fprintf(stdout,"  -fasta        Input file is fasta 'munged' format for proteins\n");
  fprintf(stdout,"Options for converting msp crunch output\n");
  fprintf(stdout,"  -fetch  [%s]  Pipe to open to fetch protein sequences\n",fetch_string);
  fprintf(stdout,"  -window [%s]  Window size of dna to consider\n",window_str);
  fprintf(stdout,"  -wing   [%s]  Size of 'wing' sequences on each end\n",wing_str);
  fprintf(stdout,"  -min    [%s]  Minimum blast score to trigger genewise\n",min_score_str);
  fprintf(stdout,"  -cut    [%s]  Minimum genewise score to predict on\n",cut_off_str);
  fprintf(stdout,"Output options\n");
  fprintf(stdout,"  -pretty       show pretty alignment format\n");
  fprintf(stdout,"  -ace          show ace subsequence model\n");
  fprintf(stdout,"  -gff          Gene feature format output\n");
  fprintf(stdout,"  -gener        show raw gene structure output\n");
  fprintf(stdout,"  -trans        show translation\n");
  fprintf(stdout,"  -cDNA         show cDNA\n");
  fprintf(stdout,"  -divide [%s]  divide string for outputs\n",divide_str);
  fprintf(stdout,"Algorithm options\n");
  fprintf(stdout,"  -palg   [%s]   Algorithm for protein alignments\n",palg_str);
  fprintf(stdout,"  -kbyte  [%5d]  Max number of kilobytes used in main calculation\n",kbyte);  fprintf(stdout,"Standard options\n");
  fprintf(stdout,"  -help      help\n  -version   show version and compile info\n  -silent    No messages on stderr\n  -quiet    No report on stderr\n  -erroroffstd No warning messages to stderr\n  -errorlog [file] Log warning messages to file\n");
  fprintf(stdout,"\nSee WWW help at http://www.sanger.ac.uk/Software/Wise2/postwise.shtml\n");
  exit(63);   
}

void show_version(void)
{
  fprintf(stdout,"%s\nVersion: %s\nReleased: %s\nCompiled: %s\n",program_name,VERSION_NUMBER,RELEASE_DAY,COMPILE_DATE);
  fprintf(stdout,"\nThis program is freely distributed under a Gnu Public License\n");
  fprintf(stdout,"The source code is copyright (c) GRL 1998 and others\n");
  fprintf(stdout,"There is no warranty, implied or otherwise on the performance of this program\n");
  fprintf(stdout,"For more information read the GNULICENSE file in the distribution\n\n");
  fprintf(stdout,"Credits: Ewan Birney <birney@sanger.ac.uk> wrote the core code.\n");
  fprintf(stdout,"         Portions of this code is from HMMer1, written by Sean Eddy\n");
  exit(63);   
}


int main(int argc,char ** argv) 
{
  int i;
  char * help;
  char * errlog;
  FILE * efp;
  char * temp;

  build_defaults();
  
  if( strip_out_boolean_argument(&argc,argv,"help") == TRUE ) {
    show_help(NULL);
  }

  if( (strip_out_boolean_argument(&argc,argv,"version")) == TRUE ) {
    show_version();
  }

  if( (strip_out_boolean_argument(&argc,argv,"silent")) == TRUE ) {
    erroroff(REPORT);
    erroroff(INFO);
    erroroff(WARNING);
  }

  if( (strip_out_boolean_argument(&argc,argv,"quiet")) == TRUE ) {
    erroroff(REPORT);
    erroroff(INFO);
  }

  if( (strip_out_boolean_argument(&argc,argv,"erroroffstd")) == TRUE ) {
    errorstderroff(WARNING);
  }


  if( (errlog=strip_out_assigned_argument(&argc,argv,"errlog")) != NULL ) {
    if( add_log_filename(errlog) == FALSE ) {
      fatal("Could not use %s as a error log file\n",errlog);
    } else {
      warn("Logging errors to %s as well as stderr",errlog);
      errorlogon(WARNING);
    }
  }

  if( (temp = strip_out_assigned_argument(&argc,argv,"gap")) != NULL )
    gap_str = temp;

  if( (temp = strip_out_assigned_argument(&argc,argv,"ext")) != NULL )
    ext_str = temp;

  if( (temp = strip_out_assigned_argument(&argc,argv,"matrix")) != NULL )
    matrix_file = temp;

  if( (temp = strip_out_assigned_argument(&argc,argv,"codon")) != NULL )
    codon_file = temp;

  if( (temp = strip_out_assigned_argument(&argc,argv,"gene")) != NULL )
    gene_file = temp;

  if( (temp = strip_out_assigned_argument(&argc,argv,"fetch")) != NULL )
    fetch_string = temp;

  if( (temp = strip_out_assigned_argument(&argc,argv,"window")) != NULL )
    window_str = temp;

  if( (temp = strip_out_assigned_argument(&argc,argv,"wing")) != NULL )
    wing_str = temp;

  if( (temp = strip_out_assigned_argument(&argc,argv,"min")) != NULL )
    min_score_str = temp;

  if( (temp = strip_out_assigned_argument(&argc,argv,"cut")) != NULL )
    cut_off_str = temp;

  if( (temp = strip_out_assigned_argument(&argc,argv,"palg")) != NULL )
    palg_str = temp;

  if( (temp = strip_out_assigned_argument(&argc,argv,"divide")) != NULL )
    divide_str = temp;

  if( strip_out_boolean_argument(&argc,argv,"ace") == TRUE ) 
    show_ace = TRUE;

  if( strip_out_boolean_argument(&argc,argv,"pretty") == TRUE ) 
    show_alignments = TRUE;

  if( strip_out_boolean_argument(&argc,argv,"gener") == TRUE ) 
    show_raw_gene = TRUE;

  if( strip_out_boolean_argument(&argc,argv,"gff") == TRUE ) 
    show_gff = TRUE;

  if( strip_out_boolean_argument(&argc,argv,"trans") == TRUE ) 
    show_translation = TRUE;

  if( strip_out_boolean_argument(&argc,argv,"cdna") == TRUE ) 
    show_cdna = TRUE;

  if( strip_out_boolean_argument(&argc,argv,"genes") == TRUE ) 
    show_pretty_gene = TRUE;

  if( strip_out_boolean_argument(&argc,argv,"fasta") == TRUE ) 
    use_msp_crunch = FALSE;


  if( (temp = strip_out_assigned_argument(&argc,argv,"kbyte")) != NULL )
    kbyte_str = temp;

  if( show_gff == FALSE && show_ace == FALSE && show_translation == FALSE && show_pretty_gene == FALSE ) {
    show_alignments  = TRUE;
    show_translation = TRUE;
    show_pretty_gene = TRUE;
  }

  strip_out_remaining_options_with_warning(&argc,argv);
  

  if( argc != 3 ) {
    warn("Wrong number of arguments (expect 2)! Arg line looked like (after option processing)");
    for(i=1;i<argc;i++) {
      fprintf(stderr,"%s\n",argv[i]);
    }
    show_short_help();

  }

  dna_file = argv[1];
  exblx_file = argv[2];
  pgl_file = argv[2]; /* if we are not using msp crunch files! */
   

  if( build_objects() == FALSE ) 
    fatal("Cannot build objects!");

  if( build_gene_prediction() == FALSE ) 
    fatal("Cannot build gene prediction. Quite possibly a bug");

  show_output();

}






