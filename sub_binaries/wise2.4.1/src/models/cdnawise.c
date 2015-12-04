#include "dyna.h"
#include "genestats.h"
#include "geneparser4.h"
#include "cdnawise10.h"
#include "version.h" 
#include "geneutil.h"
#include "geneoutput.h"

char * program_name = "cdnawise";

char * codon_file = "codon.table";

void show_version(FILE * ofp)
{
  fprintf(ofp,"%s\nVersion: %s\nReleased: %s\nCompiled: %s\n",program_name,VERSION_NUMBER,RELEASE_DAY,COMPILE_DATE);


  fprintf(ofp,"\nThis program is freely distributed under a Gnu Public License\n");
  fprintf(ofp,"The source code is copyright (c) EMBL (2001) and others\n");
  fprintf(ofp,"There is no warranty, implied or otherwise on the performance of this program\n");
  fprintf(ofp,"For more information read the GNULICENSE file in the distribution\n\n");
  fprintf(ofp,"Credits: Ewan Birney <birney@sanger.ac.uk> wrote the core code.\n");
  exit(63);   
}


void show_help(FILE * ofp)
{

  fprintf(ofp,"%s (%s)\n",program_name,VERSION_NUMBER);
  fprintf(ofp,"cdnawise <cdna-file> <genomic-dna-file> in fasta format\n");
  fprintf(ofp,"cDna sequence options\n");
  fprintf(ofp,"  -s               start position in dna\n");
  fprintf(ofp,"  -t               end position in dna\n");
  fprintf(ofp,"Genomic Dna sequence options\n");
  fprintf(ofp,"  -u               start position in dna\n");
  fprintf(ofp,"  -v               end position in dna\n");
  fprintf(ofp,"Match algorithm options\n");
  fprintf(ofp,"  -g               gap open penalty\n");
  fprintf(ofp,"  -e               gap ext penalty\n");

  show_help_GeneModelParam(ofp);
  show_help_GeneOutputPara(ofp);

  show_help_DPRunImpl(ofp);
  show_standard_options(ofp);
}



int main(int argc,char ** argv)
{
  Sequence * cdna;
  Sequence * gen;
  Sequence * active_gen;
  Sequence * active_cdna;

  int i;
  int dstart = -1;
  int dend   = -1;

  int cstart = -1;
  int cend   = -1;

  int gapopen = 12;
  int gapext  = 2;

  CodonTable * ct = NULL;
  CodonMatrixScore * cm = NULL;
  RandomCodon * rndcodon = NULL;
  RandomCodonScore * rndcodonscore = NULL;
  DnaMatrix * dm   = NULL;

  DPRunImpl * dpri = NULL;
 
  GeneModel * gm;
  GeneModelParam * gmp;
  GeneStats * gs;
  GeneParser21 * gp21;
  GeneParser21Score * gp21s;
  GeneParser4Score * gp;

  GeneOutputData data;
  GeneOutputPara * geneout;


  ComplexSequenceEvalSet * cdna_cses;
  ComplexSequenceEvalSet * gen_cses;

  ComplexSequence * cs_cdna;
  ComplexSequence * cs_gen;
  
  Genomic * gent;
  GenomicRegion * gr;

  CompMat  * cmat;
  CompProb * cprob;
  char * matfile = "blosum62.bla";
  Protein * trans;

  PackAln * pal;
  AlnBlock * alb;

  FILE * ofp = stdout;

  dpri = new_DPRunImpl_from_argv(&argc,argv);
  gmp  = new_GeneModelParam_from_argv(&argc,argv);

  geneout = new_GeneOutputPara_from_argv(&argc,argv);

  strip_out_integer_argument(&argc,argv,"u",&dstart);
  strip_out_integer_argument(&argc,argv,"v",&dend);

  strip_out_integer_argument(&argc,argv,"s",&cstart);
  strip_out_integer_argument(&argc,argv,"t",&cend);

  strip_out_integer_argument(&argc,argv,"g",&gapopen);
  strip_out_integer_argument(&argc,argv,"e",&gapext);


  strip_out_standard_options(&argc,argv,show_help,show_version);


  ct = read_CodonTable_file(codon_file);

  cmat = read_Blast_file_CompMat(matfile);
  cprob = CompProb_from_halfbit(cmat);
  cm = naive_CodonMatrixScore_from_prob(ct,cprob);
  
  gm = GeneModel_from_GeneModelParam(gmp);

  cdna = read_fasta_file_Sequence(argv[1]);
  gen = read_fasta_file_Sequence(argv[2]);

  if( dstart != -1 || dend != -1 ) {
    if( dstart == -1 ) {
      dstart = 1;
    }
    if( dend == -1 ) {
      dend = gen->len;
    }
    active_gen = magic_trunc_Sequence(gen,dstart,dend);
  } else {
    active_gen = hard_link_Sequence(gen);
  }

  if( cstart != -1 || cend != -1 ) {
    if( cstart == -1 ) {
      cstart = 1;
    }
    if( cend == -1 ) {
      cend = gen->len;
    }
    active_cdna = magic_trunc_Sequence(gen,cstart,cend);
  } else {
    active_cdna = hard_link_Sequence(gen);
  }

  

  rndcodon = RandomCodon_from_raw_CodonFrequency(gm->codon,ct);
  fold_in_RandomModelDNA_into_RandomCodon(rndcodon,gm->rnd);

  rndcodonscore = RandomCodonScore_from_RandomCodon(rndcodon);

  assert(active_cdna);
  assert(active_gen);

  cdna_cses = default_cDNA_ComplexSequenceEvalSet();
  gen_cses  = new_ComplexSequenceEvalSet_from_GeneModel(gm);

  cs_cdna = new_ComplexSequence(active_cdna,cdna_cses);
  cs_gen  = new_ComplexSequence(active_gen,gen_cses);

  gp21 = std_GeneParser21();
  GeneParser21_fold_in_RandomModelDNA(gp21,gm->rnd);
  gp21s = GeneParser21Score_from_GeneParser21(gp21);
  gp = GeneParser4Score_from_GeneParser21Score(gp21s);
 
  dm = identity_DnaMatrix(Probability2Score(halfbit2Probability(1)),Probability2Score(halfbit2Probability(-1)));

  assert(cs_cdna);
  assert(cs_gen);
  assert(gp);
  assert(rndcodonscore);
  assert(dm);
  assert(dpri);
  
  /*  show_CodonMatrixScore(cm,ct,ofp);*/
  /*
  pal = PackAln_bestmemory_CdnaWise10(cs_cdna,cs_gen,gp,cm,rndcodonscore,dm,
				      Probability2Score(halfbit2Probability(-gapopen)),
				      Probability2Score(halfbit2Probability(-gapext)),
				      Probability2Score(halfbit2Probability(-5)),
				      +50,
				      NULL,
				      dpri);
  */
  pal = PackAln_bestmemory_CdnaWise10(cs_cdna,cs_gen,gp,cm,rndcodonscore,dm,
				      0,
				      0,
				      Probability2Score(halfbit2Probability(-5)),
				      +50,
				      NULL,
				      dpri);


  alb = convert_PackAln_to_AlnBlock_CdnaWise10(pal);

  gent = Genomic_from_Sequence(gen);
  assert(gent);

  gr = new_GenomicRegion(gent);
  assert(gr);


  add_Genes_to_GenomicRegion_GeneWise(gr,active_gen->offset,active_gen->end,alb,cdna->name,0,NULL);
				      

  data.pal = pal;
  data.alb = alb;
  data.gr  = gr;
  data.gen = gent;
  data.ct  = ct;

  show_GeneOutput(&data,geneout,stdout);

  return 0;
}

