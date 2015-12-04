#include "version.h"
#include "genestats.h"
#include "geneutil.h"
#include "standardout.h"
#include "pseudowise7.h"
#include "genedisplay.h"

#include "genewise6.h"

char * program_name = "pseudowise";


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
  fprintf(ofp,"%s query-protein-as-fasta query-cdna-as-fasta target-genomic-region\n",program_name);
  fprintf(ofp,"  -pretty   show pretty output\n");
  show_help_GeneModelParam(ofp);

  show_help_StandardOutputOptions(ofp);

  show_help_DPRunImpl(ofp);

  show_standard_options(ofp);
}

typedef struct {
  double gene_model;
  double pseudo_model;
  int synonymous;
  int nonsynonymous;
  int unlikely;
  int frameshift;
  int identical;
  int stop;
  int codon;
  int intron;
} PseudoGeneAssessment;

PseudoGeneAssessment pseudo_gene_assessment(AlnBlock * alb,Sequence * pep,Sequence * cdna,Sequence *gen,CompMat * protein,CodonTable * ct,double subs,double indel)
{

  int aa;
  int gen_pos;
  int pep_pos;

  int i,j;
  AlnColumn * alc;


  PseudoGeneAssessment out;

  out.gene_model    = 1.0;
  out.pseudo_model  = 1.0;
  out.synonymous    = 0;
  out.nonsynonymous = 0;
  out.identical     = 0;
  out.unlikely      = 0;
  out.frameshift    = 0;
  out.stop          = 0;
  out.codon         = 0;
  out.intron        = 0;

  for(alc=alb->start;alc!=NULL;alc = alc->next ) {
    gen_pos = alc->alu[1]->start+1;
    pep_pos = alc->alu[0]->start+1;

    if( strcmp(alc->alu[1]->text_label,"CODON") == 0 && strcmp(alc->alu[0]->text_label,"MATCH_STATE") == 0 ) {      
      out.codon++;
      if( is_stop_codon(codon_from_seq(gen->seq+gen_pos),ct) ) {
	out.stop++;
	continue;
      }

      aa = aminoacid_from_seq(ct,gen->seq+gen_pos);
     
      /* see if there has been a synomous codon change here */
      if( aa == pep->seq[pep_pos] ) {
	if( cdna->seq[pep_pos*3] != gen->seq[gen_pos] ||
	    cdna->seq[pep_pos*3+1] != gen->seq[gen_pos+1] ||
	    cdna->seq[pep_pos*3+2] != gen->seq[gen_pos+2] ) {
	  out.synonymous++;
	} else {
	  /* identical codon */
	  out.identical++;
	}
      } else {
	out.nonsynonymous++;
	if( protein->comp[aa][pep->seq[pep_pos]] < 0 ) {
	  out.unlikely++;
	}
      }
      
    } else if ( strstr(alc->alu[1]->text_label,"SEQUENCE_DELETION") != NULL ) {
      out.frameshift++;
    } else if ( strstr(alc->alu[1]->text_label,"SEQUENCE_INSERTION") != NULL ) {
      out.frameshift++;
    } else if( strstr(alc->alu[1]->text_label,"5SS") != NULL ) {
      out.intron++;
    }
  }

  return out;
}


int check_pep_cdna_synchrony(Sequence * cdna,Sequence * pep,CodonTable * ct)
{
  int i;

  assert(cdna);
  assert(pep);

  if( ! (is_dna_SequenceType(cdna->type)) ) {
    fatal("cDNA sequnece is not DNA!");
  }
  
  for(i=0;i<pep->len;i++) {
    if( pep->seq[i] != aminoacid_from_seq(ct,cdna->seq+(i*3)) ) {
      fatal("cDNA sequence does not match up with peptide sequence");
    }
  }

  return 0;
}


int main(int argc,char ** argv)
{
  DPRunImpl * dpri     = NULL;
  GeneModelParam * gmp = NULL;
  GeneModel * gm       = NULL;
  CodonTable * ct      = NULL;

  Sequence * prot      = NULL;
  Protein  * protein   = NULL;
  Sequence * cdna      = NULL;
  Sequence * genomic   = NULL;
  Genomic  * gen       = NULL;

  ComplexSequenceEvalSet * cses = NULL;
  ComplexSequence * cseq = NULL;

  ThreeStateModel * tsm;
  GeneWise * gw;

  GeneWiseScore * gws;
  GeneParser4Score * gps;
  GeneParser4 * gp;

  GenomicRegion * gr;

  CodonMapper  * cm = NULL;

  StandardOutputOptions * std_opt;

  boolean show_pretty = FALSE;

  PackAln * pal;
  AlnBlock * alb;
  
  PseudoGeneAssessment pa;

  boolean use_genewise= TRUE;

  char * matrix = "blosum62.bla";
  CompMat * mat;
  double indel_rate = 0.01;
  double subs_rate  = 0.01;
  int offset;
  RandomModel * rm;
  RandomModelDNA * rmd;


  dpri = new_DPRunImpl_from_argv(&argc,argv);
  if( dpri == NULL ) {
    fatal("Unable to build DPRun implementation. Bad arguments");
  }

  gmp = new_GeneModelParam_from_argv(&argc,argv);

  strip_out_boolean_def_argument(&argc,argv,"pretty",&show_pretty);

  std_opt = new_StandardOutputOptions_from_argv(&argc,argv);
  strip_out_boolean_def_argument(&argc,argv,"gw",&use_genewise);

  strip_out_standard_options(&argc,argv,show_help,show_version);
  if( argc != 4 ) {
    show_help(stdout);
    exit(12);
  }


  if((gm=GeneModel_from_GeneModelParam(gmp)) == NULL ) {
    fatal("Could not build gene model");
  }

  ct= read_CodonTable_file("codon.table");


  prot    = read_fasta_file_Sequence(argv[1]);
  cdna    = read_fasta_file_Sequence(argv[2]);
  genomic = read_fasta_file_Sequence(argv[3]);

  if( prot == NULL || cdna == NULL || genomic == NULL ) {
    fatal("need to have 3 sequences- protein guide, cdna and genomic");
  }


  offset = check_pep_cdna_synchrony(cdna,prot,ct);

  rm = default_RandomModel();
  rmd = RandomModelDNA_std();

  mat = read_Blast_file_CompMat(matrix);

  cm = flat_CodonMapper(ct);

  sprinkle_errors_over_CodonMapper(cm,subs_rate);

  protein = Protein_from_Sequence(prot);
  
  tsm = ThreeStateModel_from_half_bit_Sequence(protein,mat,rm,-11,-1);

  gw = GeneWise_from_ThreeStateModel(tsm,NULL,cm,0.0001,NULL);

  cses = new_ComplexSequenceEvalSet_from_GeneModel(gm);
  cseq = new_ComplexSequence(genomic,cses);

  gp  = std_GeneParser4(indel_rate,indel_rate*0.1);
 

  gps = GeneParser4Score_from_GeneParser4(gp);

  GeneWise_fold_in_RandomModelDNA(gw,rmd);

  gws = GeneWiseScore_from_GeneWise(gw);

  if( use_genewise ) {
    pal = PackAln_bestmemory_GeneWise6(gws,cseq,gps,NULL,dpri);
    alb = convert_PackAln_to_AlnBlock_GeneWise6(pal);
  } else {
    pal = PackAln_bestmemory_PseudoWise7(gws,cseq,gps,Probability2Score(indel_rate),NULL,dpri);
    alb = convert_PackAln_to_AlnBlock_PseudoWise7(pal);
  }
  

  show_StandardOutputOptions(std_opt,alb,pal,"//",stdout);

  pa = pseudo_gene_assessment(alb,
			      prot,
			      cdna,
			      genomic,
			      mat,ct,subs_rate,indel_rate);

  printf("Synonymous     : %d\n",pa.synonymous);
  printf("Nonsynonymous  : %d\n",pa.nonsynonymous);
  printf("Ka/Ks          : %.2f\n",(double)pa.nonsynonymous/(double)pa.synonymous);
  printf("Unlikely       : %d\n",pa.unlikely);
  printf("Identical      : %d\n",pa.identical);
  printf("Stop           : %d\n",pa.stop);
  printf("Total codons   : %d\n",pa.codon);
  printf("Frameshift     : %d\n",pa.frameshift);
  printf("Intron         : %d\n",pa.intron);

  printf("//\n");

  gen = Genomic_from_Sequence(genomic);
  gr  = new_GenomicRegion_from_GeneWise(gen,FALSE,alb);

  show_pretty_GenomicRegion(gr,0,stdout);

  printf("//\n");
  if( show_pretty) {
    protgene_ascii_display(alb,prot->seq,prot->name,prot->offset,gen,ct,15,50,FALSE, stdout);
  }


  return 0;

}
