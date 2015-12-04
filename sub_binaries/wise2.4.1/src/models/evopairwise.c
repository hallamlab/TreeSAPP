#define WISE2_CROSS_HMMER2
#include "wise2xhmmer2.h"

#include "threestatedp.h"
#include "threestateloop.h"
#include "seqaligndisplay.h"
#include "version.h"

#include "dnamatrix.h"


char * program_name = "evopairwise";


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
  fprintf(ofp,"%s align_cds_pair hmm_set\n",program_name);

  show_help_DPRunImpl(ofp);

  show_standard_options(ofp);

}

Score logl_pseudogene(char * ref,char * diff,DnaProbMatrix * dm)
{
  int i;
  Score s = 0;

  for(i=0;i<3;i++) {
    s += Probability2Score(dm->prob[base_from_char(ref[i])][base_from_char(diff[i])]);
  }

  return s;
}

Score logl_negative_selection(char * ref,char * diff,ThreeStateUnit * unit,CodonTable * ct,DnaProbMatrix * dm)
{
  int i;
  Score s = 0;
  char ref_aa;
  char diff_aa;

  /* we have to assess this position having changed */
  for(i=0;i<3;i++) {
    s += Probability2Score(dm->prob[base_from_char(ref[i])][base_from_char(diff[i])]);
  }

  /* if the position has not changed, then we know it could not have been selected */


  ref_aa  = aminoacid_from_seq(ct,ref);
  diff_aa = aminoacid_from_seq(ct,diff);

  if( ref_aa == diff_aa ) {
    return s;
  }

  /* else we add the difference in probability between the two amino acids */
  /*
  fprintf(stdout,"%c vs %c has %d plays %d for total of %d\n",ref_aa,diff_aa,
	  Probability2Score(unit->match_emission[ref_aa-'A']),
	  Probability2Score(unit->match_emission[diff_aa-'A']),
	  Probability2Score(unit->match_emission[diff_aa-'A'])  - Probability2Score(unit->match_emission[ref_aa-'A'])
	  );
  */


  s += Probability2Score(unit->match_emission[diff_aa-'A'])  - Probability2Score(unit->match_emission[ref_aa-'A']);

  return s;
}


Score logl_positive_selection(char * ref,char * diff,RandomModel * rm,CodonTable *ct,DnaProbMatrix * dm)
{
  int i;
  Score s = 0;
  char ref_aa;
  char diff_aa;

  /* we have to assess this position having changed */
  for(i=0;i<3;i++) {
    s += Probability2Score(dm->prob[base_from_char(ref[i])][base_from_char(diff[i])]);
  }

  /* if the position has not changed, then we know it would have been selected */


  ref_aa  = aminoacid_from_seq(ct,ref);
  diff_aa = aminoacid_from_seq(ct,diff);

  if( ref_aa == diff_aa ) {
    return s;
  }

  /* else we add the probability of seeing this amino acid*/

  s += Probability2Score(rm->aminoacid[diff_aa-'A']);

  return s;  
}

void show_verbose_evo(AlnBlock * alb,ThreeStateModel * tsm,Sequence * ref,Sequence * diff,CodonTable * ct,FILE * ofp)
{
  AlnColumn * alc;
  Protein * hmmp;

  Sequence * ref_trans;
  Sequence * diff_trans;

  DnaProbMatrix * negative_dm;
  DnaProbMatrix * pseudo_dm;
  
  int i;
  int count = 0;
  double est_mutation = 0.0;

  int dna_offset;

  Score total_pseudo = 0;
  Score total_neg = 0;
  Score pseudo = 0;
  Score neg = 0;

  int count_ref_positive = 0;
  int count_ref_negative = 0; 

  int count_ref_negative_0_5   = 0;
  int count_ref_negative_5_10  = 0;
  int count_ref_negative_10_15 = 0;

  int syn_sites = 0;
  int nonsyn_sites = 0;

  int syn_changes = 0;
  int nonsyn_changes = 0;

  int diff_score;

  char diff_aa;
  char ref_aa;

  int score_ratio = 0;
  Score score_neg_5  = Probability2Score(Bits2Probability(-5.0));
  Score score_neg_10 = Probability2Score(Bits2Probability(-10.0));


  int k;

  for(i=0;i<ref->len;i+=3) {

    /* if this has changed, then it is definitely non syn */
    if( aminoacid_from_seq(ct,ref->seq+i) != aminoacid_from_seq(ct,diff->seq+i)) {
      for(k=0;k<3;k++) {
	if( ref->seq[i+k] != diff->seq[i+k] ) {
	  nonsyn_changes++;
	}
      }
    } else {
      /* could still be syn change */
      for(k=0;k<3;k++) {
	if( ref->seq[i+k] != diff->seq[i+k] ) {
	  syn_changes++;
	}
      }
    }

    /* calculate the sites. There is always 2 non syn sites */

    nonsyn_sites += 2;

    if( four_fold_sites_CodonTable(ct,ref->seq+i) > 0 ) {
      syn_sites++;
    } else {
      nonsyn_sites += 1;
    } 
  }

  for(i=0;i<ref->len;i++) {
    if( ref->seq[i] != diff->seq[i] ) {
      count++;
    }
  }


  est_mutation = (double)count / (double)ref->len;


  pseudo_dm = DnaProbMatrix_from_match(1.0 - est_mutation,NMaskType_BANNED);
  negative_dm = DnaProbMatrix_from_match(1.0 - (est_mutation*2),NMaskType_BANNED);


  ref_trans = translate_Sequence(ref,ct);
  diff_trans = translate_Sequence(diff,ct);
  
  hmmp = pseudo_Protein_from_ThreeStateModel(tsm);

  for(alc=alb->start;alc != NULL;alc = alc->next) {
    /*    fprintf(stdout,"In position %s\n",alc->alu[0]->text_label); */
    if( strcmp(alc->alu[0]->text_label,"SEQUENCE") == 0 &&
	strcmp(alc->alu[1]->text_label,"SEQUENCE") == 0 ) {
      dna_offset = alc->alu[1]->end*3;

      pseudo = 	      logl_pseudogene(ref->seq+dna_offset,diff->seq+dna_offset,pseudo_dm);
      neg = 	      logl_negative_selection(ref->seq+dna_offset,diff->seq+dna_offset,tsm->unit[alc->alu[0]->end],ct,
					      pseudo_dm);

      /*
      fprintf(ofp,"Position %d [%c], vs %d [%c,%c] Scores Negative %d, Pseudo %d\n",
	      alc->alu[0]->end,hmmp->baseseq->seq[alc->alu[0]->end],
	      alc->alu[1]->end,ref_trans->seq[alc->alu[1]->end],diff_trans->seq[alc->alu[1]->end],
	      neg,
	      pseudo
	      );
      */

      ref_aa = ref_trans->seq[alc->alu[1]->end];
      diff_aa = diff_trans->seq[alc->alu[1]->end]; 
      if( ref_aa != diff_aa  ) {
	score_ratio += Probability2Score(tsm->unit[alc->alu[0]->end]->match_emission[ref_aa-'A']) - Probability2Score(tsm->unit[alc->alu[0]->end]->match_emission[diff_aa-'A']);

	diff_score = Probability2Score(tsm->unit[alc->alu[0]->end]->match_emission[ref_aa-'A']) - Probability2Score(tsm->unit[alc->alu[0]->end]->match_emission[diff_aa-'A']);
 
	if( diff_score < 0) {
	  count_ref_negative++;
	  if( diff_score > score_neg_5 ) {
	    count_ref_negative_0_5++;
	  } else if ( diff_score > score_neg_10 ) {
	    count_ref_negative_5_10++;
	  } else {
	    count_ref_negative_10_15++;
	  }
	} else {
	  count_ref_positive++;
	}

      }

      total_pseudo += pseudo;
      total_neg += neg;
    }
  }

  fprintf(ofp,"%s\t%s\t%.2f\t%d\t%d\t%d\t%d\t%d\n",ref->name,hmmp->baseseq->name,Score2Bits(score_ratio),
	  count_ref_positive,count_ref_negative,
	  count_ref_negative_0_5,
	  count_ref_negative_5_10,
	  count_ref_negative_10_15);


  /*
  fprintf(ofp,"%s,%s Total Pseudo %d vs Negative %d, Ratio %.4f  Positive %d Negative %d Score %.2f Syn %d Changes %d NonSyn %d Changes %d\n",ref->name,hmmp->baseseq->name,total_pseudo,total_neg,Score2Bits(total_neg-total_pseudo),count_ref_positive,count_ref_negative,Score2Bits(score_ratio),syn_sites,syn_changes,nonsyn_sites,nonsyn_changes);
  */

  free_Protein(hmmp);
	      
}

int main(int argc,char ** argv)
{
  int i;
  SequenceSet * in;
  Sequence * trans;
  ThreeStateDB * tsd;

  DPRunImpl * dpri;
  CodonTable * ct;

  int return_status;
  ThreeStateModel * tsm;
  ThreeStateScore * tss;
  Protein * hmmp;
  ComplexSequence * cs;
  ComplexSequenceEvalSet * cses;

  PackAln * pal;
  AlnBlock * alb;

  int show_align = 0;
  int show_alb   = 0;
  int show_verbose = 1;
  int show_trans = 0;

  ct = read_CodonTable_file("codon.table");

  cses = default_aminoacid_ComplexSequenceEvalSet();
  
  dpri = new_DPRunImpl_from_argv(&argc,argv);

  strip_out_boolean_def_argument(&argc,argv,"pretty",&show_align);

  strip_out_boolean_def_argument(&argc,argv,"alb",&show_alb);

  strip_out_boolean_def_argument(&argc,argv,"trans",&show_trans);

  if( argc != 3 ) {
    show_help(stdout);
    exit(63);
  }


  in = read_fasta_SequenceSet_file(argv[1]);

  tsd = HMMer2_ThreeStateDB(argv[2]);

  assert(in);
  assert(tsd);
  assert(in->len == 2);

  trans = translate_Sequence(in->set[0],ct);

  if( show_trans ) {
    write_fasta_Sequence(trans,stdout);
  }

  cs = new_ComplexSequence(trans,cses);

  open_ThreeStateDB(tsd);

  while( (tsm = read_TSM_ThreeStateDB(tsd,&return_status)) != NULL ) {
    fold_RandomModel_into_ThreeStateModel(tsm,tsm->rm);
    set_startend_policy_ThreeStateModel(tsm,TSM_local,10,1.0);

    tss = ThreeStateScore_from_ThreeStateModel(tsm);
    hmmp = pseudo_Protein_from_ThreeStateModel(tsm);

    pal = PackAln_bestmemory_ThreeStateLoop(tss,cs,NULL,dpri);
    alb = convert_PackAln_to_AlnBlock_ThreeStateLoop(pal);

    if( show_alb ) {
      show_flat_AlnBlock(alb,stdout);
    }

    if( show_align ) {
      write_pretty_seq_align(alb,hmmp->baseseq,trans,15,50,stdout);      
    }
    if( show_verbose ) {
      show_verbose_evo(alb,tsm,in->set[0],in->set[1],ct,stdout);
    }
      
  }
  

}
