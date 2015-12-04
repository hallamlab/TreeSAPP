#include "alignwisedp.h"


#include "version.h"
#include "genestats.h"
#include "geneutil.h"
#include "standardout.h"


char * program_name = "alignwise";



void show_Gene_debug(AlnBlock * alb,Sequence * gen,CodonTable * ct,FILE * ofp)
{
  AlnColumn * alc;
  int is_reversed;
  char seq[5];

  for(alc = alb->start;alc != NULL;alc = alc->next ) {
    if( strstr(alc->alu[1]->text_label,"REV") != NULL ) {
      is_reversed = 1;
    } else {
      is_reversed = 0;
    }

    if( strstr(alc->alu[1]->text_label,"CODON") != NULL ) {
      if( is_reversed == 0 ) {
	fprintf(ofp,"%-2.2f  %-12s  %c%c%c  %c  [%d]\n",Score2Bits(alc->alu[0]->score[0]),
		alc->alu[1]->text_label,
		gen->seq[alc->alu[1]->start+1],
		gen->seq[alc->alu[1]->start+2],
		gen->seq[alc->alu[1]->start+3],
		aminoacid_from_seq(ct,gen->seq+alc->alu[1]->start+1),alc->alu[1]->start+1
		);
      } else {
	seq[0] = char_complement_base(gen->seq[alc->alu[1]->start+3]);
	seq[1] = char_complement_base(gen->seq[alc->alu[1]->start+2]);
	seq[2] = char_complement_base(gen->seq[alc->alu[1]->start+1]);
	seq[3] = '\0';
	fprintf(ofp,"%-2.2f  %-12s  %c%c%c  %c [%c%c%c]  [%d]\n",Score2Bits(alc->alu[0]->score[0]),
		alc->alu[1]->text_label,
		gen->seq[alc->alu[1]->start+1],
		gen->seq[alc->alu[1]->start+2],
		gen->seq[alc->alu[1]->start+3],
		aminoacid_from_seq(ct,seq),
		seq[0],
		seq[1],
		seq[2], alc->alu[1]->start+1
		);
      }
    } else {
      fprintf(ofp,"%-2.2f %-12s %d-%d\n",Score2Bits(alc->alu[0]->score[0]),
	      alc->alu[1]->text_label,alc->alu[1]->start,alc->alu[1]->end);
    }
  }


}

void show_version(FILE * ofp)
{
  fprintf(ofp,"%s\nVersion: %s\nReleased: %s\nCompiled: %s\n",program_name,VERSION_NUMBER,RELEASE_DAY,COMPILE_DATE);
  fprintf(ofp,"\nThis program is freely distributed under a Gnu Public License\n");
  fprintf(ofp,"The source code is copyright (c) EMBL and others\n");
  fprintf(ofp,"There is no warranty, implied or otherwise on the performance of this program\n");
  fprintf(ofp,"For more information read the GNULICENSE file in the distribution\n\n");
  fprintf(ofp,"Credits: Ewan Birney <birney@ebi.ac.uk>\n");
  exit(63);   
}

void show_help(FILE * ofp)
{
  fprintf(ofp,"%s multiple-alignment-as-fasta\n",program_name);

  show_help_AlignGeneModelParam(ofp);

  show_help_ShowGenomicRegionOptions(ofp);

  fprintf(ofp,"Debugging parameter space\n");
  fprintf(ofp,"   -show_aligngene - show odd-ratios for likelihoods\n");
  fprintf(ofp,"   -show_alignmodel [file] file of gene structure to show\n");
  fprintf(ofp,"   -singleexon [180] - minimum coding length to show single exon\n");
  fprintf(ofp,"   -multiexon  [90] - minimum coding length to show multi exon\n");


  show_help_DPRunImpl(ofp);

  show_help_StandardOutputOptions(ofp);

  show_standard_options(ofp);
}

int main(int argc,char ** argv)
{
  int i;

  DPRunImpl * dpri = NULL;

  FILE * ifp;
  SeqAlign   * al;


  AlignGeneModelParam * agmp;
  AlignGeneModel * agm;
  AlignGeneModelScore * alignscore;
  AlignGeneModelFrame * frame;
  

  PackAln * pal;
  AlnBlock * alb;
  AlnBlock * collapsed;

  Genomic * genomic = NULL;
  GenomicRegion * gr = NULL;
  GenomicRegion * temp_gr = NULL;

  StandardOutputOptions * std_opt;
  ShowGenomicRegionOptions * sgro;
  
  char * dump_packaln = NULL;
  char * read_packaln = NULL;
  FILE * packifp = NULL;

  boolean show_trans    = 1;
  boolean show_collapsed = 0;
  boolean show_gene_raw = 0;
  boolean show_gene_model = FALSE;
  int genedebug = 0;

  GenomicRegion * inputgr = NULL;
  char * input_gr_file = NULL;

  Score intronscore;
  Probability intronopen = 0.0001;

  Score genescore;
  Probability geneopen_factor = 0.00000005;

  int multiexon = 90;
  int singleexon = 180;
 

  dpri = new_DPRunImpl_from_argv(&argc,argv);
  if( dpri == NULL ) {
    fatal("Unable to build DPRun implementation. Bad arguments");
  }


  agmp = new_AlignGeneModelParam_from_argv(&argc,argv);

  std_opt = new_StandardOutputOptions_from_argv(&argc,argv);
  sgro = new_ShowGenomicRegionOptions_from_argv(&argc,argv);
  
  strip_out_float_argument(&argc,argv,"intron",&intronopen);
  strip_out_float_argument(&argc,argv,"genefactor",&geneopen_factor);
 
  strip_out_integer_argument(&argc,argv,"singleexon",&singleexon);
  strip_out_integer_argument(&argc,argv,"multiexon",&multiexon);

  strip_out_boolean_def_argument(&argc,argv,"show_aligngene",&show_gene_model);
  input_gr_file = strip_out_assigned_argument(&argc,argv,"show_alignmodel");

  strip_out_boolean_def_argument(&argc,argv,"genedebug",&genedebug);
  strip_out_boolean_def_argument(&argc,argv,"albcoll",&show_collapsed);

  dump_packaln = strip_out_assigned_argument(&argc,argv,"dump");
  read_packaln = strip_out_assigned_argument(&argc,argv,"recover");


  strip_out_standard_options(&argc,argv,show_help,show_version);
  if( argc != 2 ) {
    show_help(stdout);
    exit(12);
  }


  intronscore = Probability2Score(intronopen);

  genescore = Probability2Score(intronopen*geneopen_factor);

  info("Genescore is %2.4f and intron score is %2.4f\n",Score2Bits(genescore),Score2Bits(intronscore));

  if((ifp = openfile(argv[1],"r")) == NULL ) {
    fatal("Could not open file %s",argv[1]);
  }

  al = read_fasta_SeqAlign(ifp);
  
  info("Read in alignment, anchor %s, length %d species %d\n",al->seq[0]->name,al->seq[0]->len,al->len);

  assert(al);
  assert(al->len > 1);

  agm = create_AlignGeneModel(al,agmp);

  if( input_gr_file != NULL ) {
    if((ifp = openfile(input_gr_file,"r")) == NULL ) {
      fatal("Could not open file %s",input_gr_file);
    }
    inputgr = read_genes_GenomicRegion(ifp);
    if( inputgr->len <= 0 ) {
      warn("Alignment input model has no genes!");
    }
  }

  if( show_gene_model == TRUE ) {
    show_AlignGeneModel(agm,al,agmp->ct,inputgr,stdout,agmp);
    exit(0);
  }

    

  alignscore = AlignGeneModelScore_from_AlignGeneModel(agm);

  frame = AlignGeneModelFrame_alloc();
  frame->len = 1;

  if( read_packaln != NULL ) {
    packifp = openfile(read_packaln,"r");
    if( packifp == NULL ) {
      fatal("File %s is unopenable - ignoring dump command",dump_packaln);
    } else {
      pal = read_simple_PackAln(packifp);
    }
  } else {
    pal = PackAln_bestmemory_AlignWise(frame,alignscore,intronscore,genescore,NULL,dpri);
  }
  
  if( dump_packaln != NULL ) {
    packifp = openfile(dump_packaln,"w");
    if( packifp == NULL ) {
      warn("File %s is unopenable - ignoring dump command",dump_packaln);
    } else {
      show_simple_PackAln(pal,packifp);
    }
  }

  alb = convert_PackAln_to_AlnBlock_AlignWise(pal);

  genomic = Genomic_from_Sequence(al->seq[0]);
  gr = new_GenomicRegion(genomic);

  add_Genes_to_GenomicRegion_new(gr,alb);
  
  temp_gr = new_GenomicRegion_discard_short(gr,multiexon,singleexon);
  /*  free_GenomicRegion(gr);*/
  gr = temp_gr;



  show_StandardOutputOptions(std_opt,alb,pal,"//",stdout);

  if( genedebug ) {
    show_Gene_debug(alb,al->seq[0],agmp->ct,stdout); 
  }

  if( show_collapsed ) {
    collapsed = collapsed_AlnBlock(alb,1);
    mapped_ascii_AlnBlock(collapsed,Score2Bits,0,stdout);
    fprintf(stdout,"//\n");
  }

  show_GenomicRegionOptions(sgro,gr,agmp->ct,"//",stdout);

  return 0;
}
