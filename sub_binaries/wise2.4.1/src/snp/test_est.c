#include "genotype.h"
#include "locus_model_estimators.h"

int main(int argc,char ** argv)
{
  int i,j;
  FILE * ifp;
  BiGenotypeSet * bgs;
  LocusGenotypeCountSet * lgsc;
  LocusGenotypeCount * new_lgc;

  int randomise = 1;
  double sm;

  ModelProbabilitySet * mps_normal;
  ModelProbabilitySet * mps_first;
  ModelProbabilitySet * mps_second;
  ModelProbabilitySet * mps_het;
  
  double central;

  bgs = new_BiGenotypeSet();

  ifp = openfile(argv[1],"r");

  read_hapmap_genotype_file(bgs,"CEPH",ifp);

  fprintf(stderr," read CEPHs\n");

  ifp = openfile(argv[2],"r");

  read_hapmap_genotype_file(bgs,"YRI",ifp);

  fprintf(stderr," read YRIs\n");

  ifp = openfile(argv[3],"r");

  read_hapmap_genotype_file(bgs,"JAP",ifp);

  fprintf(stderr," read JAPs\n");
  
  ifp = openfile(argv[4],"r");

  read_hapmap_genotype_file(bgs,"HAN",ifp);

  fprintf(stderr," read HANs\n");
  

  fprintf(stderr,"Read in genotypes\n");

  lgsc = LocusGenotypeCountSet_from_BiGenotypeSet(bgs);

  fprintf(stderr,"converted to counts\n");

  /*  free_BiGenotypeSet(bgs); */

  for(i=0;i<lgsc->len;i++) {

    fprintf(stderr,"handling locus %d\n",i);

    if( seen_each_genotype_in_all_populations_LocusGenotypeCount(lgsc->lgc[i]) == FALSE ) {
      fprintf(stdout,"Locus %s 0 genotypes in some population\n",lgsc->lgc[i]->bl->locus_id);
      continue;
    }
 
    if( (sm = smallest_minor_allele_LocusGenotypeCount(lgsc->lgc[i])) < 0.10 ) {
      fprintf(stdout,"Locus %s has too small minor allele (%f)\n",lgsc->lgc[i]->bl->locus_id,
	      sm);
      continue;
    }


    if( randomise ) {
      new_lgc = resampled_LocusGenotypeCount(lgsc->lgc[i]);
    } else {
      new_lgc = hard_link_LocusGenotypeCount(lgsc->lgc[i]);
    }

    mps_normal = estimate_model_ModelProbabilitySet(new_lgc,0.2,10,LocusModel_Normal_Multi,0.4,20);
    mps_first  = estimate_model_ModelProbabilitySet(new_lgc,0.2,10,LocusModel_SingleHomozygous_First_Selection,0.4,20);
    mps_second = estimate_model_ModelProbabilitySet(new_lgc,0.2,10,LocusModel_SingleHomozygous_Second_Selection,0.4,20);
    mps_het   =  estimate_model_ModelProbabilitySet(new_lgc,0.2,10,LocusModel_Hetrozygous_Selection,0.4,20);
    
    fprintf(stdout,"Locus %s\tN: %f %.3f %.2f\t%.3f\tFirst (%.2f %.2f) %.3f\tHet: (%.2f %.2f) %.3f\tSecond (%.2f %.2f) %.3f\n",
	    new_lgc->bl->locus_id,Score2Bits(mps_normal->best->likelihood_score),sm,
	    chisquared_stat_LocusGenotypeCount(new_lgc),
	    mps_normal->best->estimate_first_freq[0],
	    mps_first->best->estimate_selection_hemi,
	    mps_first->best->estimate_first_freq[0],
	    Score2Bits(mps_first->best->likelihood_score - mps_normal->best->likelihood_score),
	    mps_het->best->estimate_selection_hemi,
	    mps_het->best->estimate_first_freq[0],
	    Score2Bits(mps_het->best->likelihood_score - mps_normal->best->likelihood_score),
	    mps_second->best->estimate_selection_hemi,
	    mps_first->best->estimate_first_freq[0],
	    Score2Bits(mps_second->best->likelihood_score - mps_normal->best->likelihood_score));

    fflush(stdout);

    free_ModelProbabilitySet(mps_normal);
    free_ModelProbabilitySet(mps_first);
    free_ModelProbabilitySet(mps_second);
    free_ModelProbabilitySet(mps_het);
	

    free_LocusGenotypeCount(new_lgc);

    /*
   fprintf(stdout,"Locus %s\t%f\t",new_lgc->bl->locus_id,Score2Bits(mps->best->likelihood_score));
    for(j=0;j<new_lgc->len;j++) {
      central = central_first_frequency_PopulationGenotypeCount(new_lgc->pgc[j]);
      fprintf(stdout,"Central: %f\tEstimate: %f\t",central,mps->best->estimate_first_freq[j]);
    }
    */

  }
 
}
