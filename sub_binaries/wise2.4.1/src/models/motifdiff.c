#include "transfactor.h"
#include "transregion.h"
#include "version.h"
#include "commandline.h"

char * program_name = "motifdiff";

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
  fprintf(ofp,"%s motif-library sequence\n",program_name);
  fprintf(ofp," -lr  motif library is in laurence's format\n");
  fprintf(ofp," -ben motif library is in ben's IUPAC format\n");
  fprintf(ofp," -[no]show_match  - show all matches that pass criteria, default no\n");
  fprintf(ofp," -[no]show_region - show only dense cluster regions, default yes\n");
  fprintf(ofp," -circular [no]   - for randomisation, number of positions to permute\n");

  show_help_TransFactorBuildPara(ofp);

  show_help_TransFactorMatchPara(ofp);

  show_standard_options(ofp);
}



int main(int argc,char ** argv)
{
  TransFactorSet * tfs;
  TransFactorSet * newtfs; /* only used as a temp when permuting */
  TransFactorBuildPara * tfb;
  TransFactorMatchSetCompara * comp;
  TransFactorMatchPara * matchp;
  TransFactorComparaPara * comparap;

  TransFactorMatchSet * one;
  TransFactorMatchSet * two;

  boolean use_laurence = FALSE;
  boolean use_ben = FALSE;
  boolean show_matchset = FALSE;
  boolean show_region  = FALSE;
  boolean show_compara = TRUE;
  boolean show_aln     = FALSE;

  int end_on_seq = 0;
  int seq_count = 0;

  int rotate_number = 0;

  int z;
  int is_gapped;

  int one_pos;
  int two_pos;
  int one_j;
  int two_j;
  int found_match;

  int count;
  int id;
  double perc_id;

  FILE * ifp;
  SeqAlign * sa;

  int i;

  tfb    = new_TransFactorBuildPara_from_argv(&argc,argv);
 
  matchp = new_TransFactorMatchPara_from_argv(&argc,argv);
  matchp->min_relative = 7.0;
  matchp->relative_prob_bits = 0.5;
  matchp->relative_prob = 0.5;
  matchp->type = TFM_RELATIVE;

  strip_out_boolean_def_argument(&argc,argv,"show_match",&show_matchset);
  strip_out_boolean_def_argument(&argc,argv,"show_region",&show_region);
  strip_out_integer_argument(&argc,argv,"circular",&rotate_number);


  if( strip_out_boolean_argument(&argc,argv,"lr") == TRUE ) {
    use_laurence = TRUE;
  }
  if( strip_out_boolean_argument(&argc,argv,"ben") == TRUE ) {
    use_ben = TRUE;
  }


  strip_out_standard_options(&argc,argv,show_help,show_version);


  if( argc != 3 ) {
    show_help(stdout);
    exit(12);
  }

  if( use_laurence == TRUE ) {
    tfs = read_laurence_TransFactorSet_file(argv[1]);
  } else if( use_ben == TRUE ) {
    tfs = read_ben_IUPAC_TransFactorSet_file(argv[1]);
  } else {
    tfs = read_TransFactorSet_file(argv[1]);
  }

  if( tfs->len == 0 ) {
    fatal("No transcription factors in set!");
  }

  build_TransFactorSet(tfs,tfb);

  if( rotate_number > 0 ) {
    fprintf(stdout,"PERMUTED RESULTS - %d rotation\n",rotate_number);
    newtfs = circular_permuted_TransFactorSet(tfs,rotate_number);
    free_TransFactorSet(tfs);
    tfs = newtfs;
  }

  ifp = openfile(argv[2],"r");
  if( ifp == NULL ) {
    fatal("Cannot open %s as input file",argv[2]);
  }

  fprintf(stderr,"Getting into main loop with %s\n",argv[2]);

  while( (sa = read_fasta_SeqAlign(ifp)) ) {
    
    for(i=0;i<sa->len;i++) {
      sa->seq[i]->type = SEQUENCE_DNA;
    }

    if( sa == NULL ) {
      fatal("Could not read fasta alignment in %s",argv[2]);
    }
    
    if( sa->len != 2 ) {
      fatal("Has to be a pairwise alignment");
    }

    fprintf(stderr,"About to calculate...\n");

    one = calculate_TransFactorMatchSet(sa->seq[0],tfs,matchp);

    fprintf(stderr,"calculated transfactor match (one)\n");

    two = calculate_TransFactorMatchSet(sa->seq[1],tfs,matchp);

    fprintf(stderr,"calculated transfactor match (two)\n");

    sort_by_start_TransFactorMatchSet(one);
    sort_by_start_TransFactorMatchSet(two);


    fprintf(stderr,"Sorted\n");
    one_pos = 0;
    two_pos = 0;

    

    
    
    for(i=0;i<sa->seq[0]->len;i++) {
      auto int seen_one = 0;
      auto int seen_two = 0;
      for(one_pos; one_pos < one->len && one->match[one_pos]->start < i;one_pos++) {
	;
      }
      for(two_pos; two_pos < two->len && two->match[two_pos]->start < i;two_pos++) {
	;
      }
      


      if( one_pos < one->len && one->match[one_pos]->start == i ) {
	seen_one = 1;
      } 
      if( two_pos < two->len && two->match[two_pos]->start == i ) {
	seen_two = 1;
      } 
  
      
      if( show_aln && (seen_one == 1 && seen_two == 1)) {
	fprintf(stdout,"Position %d %d,%d\n",i,seen_one,seen_two);
	fprintf(stdout,"    %.12s\n",one->target->seq+i);
	fprintf(stdout,"    %.12s\n",two->target->seq+i);
      }
      

      if( seen_one == 0 && seen_two == 0 ) {
	continue; /* nothing to do */
      }


      is_gapped = FALSE;
      z = i -3;
      if( z < 0 ) { 
	z = 0; 
      }

      for(;z<i+10 && z<sa->seq[0]->len;z++) {
	if( is_gapped_SeqAlign(sa,z) == TRUE ) {
	  is_gapped = TRUE;
	  break;
	}
      }

      if( seen_one == seen_two ) {
	id = 0;
	count = 0;
	for(z=one->match[one_pos]->start;z<one->match[one_pos]->end;z++) {
	  if( toupper(sa->seq[0]->seq[z]) == toupper(sa->seq[1]->seq[z]) ) {
	    id++;
	  } 
	  count++;
	}
	perc_id = (double)id/(double)count;
      }


      /*      fprintf(stdout," %d %d at %d,%d with %d vs %d\n",seen_one,seen_two,one_pos,two_pos,one->match[one_pos]->start,two->match[two_pos]->start);*/

      if( seen_one != seen_two ) {
	/* must have a mismatch */

	/* we no longer want to output mismatches as we are looking at relatives */
	/*
	if( one_pos >= one->len || two_pos >= two->len ) {
	  printf("%s Mismatch at %d AT_END\n",(is_gapped == TRUE ? "GAPPED" : "SOLID"), i);
	} else {
	  printf("%s Mismatch at %d [%d,%d] %s %d %s %d\n",
		 (is_gapped == TRUE ? "GAPPED" : "SOLID"),
		 i,seen_one,seen_two,
		 one->match[one_pos]->factor->name,one->match[one_pos]->start,
		 two->match[two_pos]->factor->name,two->match[two_pos]->start );
	}
	*/

      } else {
	/* if 1,1 still need to check. Will do this later */
	if( seen_one == seen_two && seen_one == 1 ) {
	  if( one->match[one_pos]->factor != two->match[two_pos]->factor ) {

	    /* more than one factor can match at the same position. Need to loop*/
	    found_match = 0;
	    for(one_j=one_pos;one_j < one->len && one->match[one_pos]->start == one->match[one_j]->start;one_j++) {
	      for(two_j=two_pos;two_j < two->len && two->match[two_pos]->start == two->match[two_j]->start;two_j++ ) {
		if( one->match[one_j]->factor == two->match[two_j]->factor && one->match[one_j]->strand == two->match[two_j]->strand ) {
		  found_match = 1;

		  printf("%s %.2f SAME %d %s %d %.2f %.2f %s %d %.2f %.2f\n",
			 (is_gapped == TRUE ? "GAPPED" : "SOLID"), perc_id,		   
			 i,
			 one->match[one_j]->factor->name,one->match[one_j]->start,
			 one->match[one_j]->bit_score,(one->match[one_j]->bit_score - Probability2Bits(one->match[one_j]->factor->max_prob)),
			 two->match[two_j]->factor->name,two->match[two_j]->start,
			 two->match[two_j]->bit_score,(two->match[two_j]->bit_score - Probability2Bits(two->match[two_j]->factor->max_prob)));
		  break;
			 
		} 
	      }
	    }
	    if( found_match == 0 ) {
	      /*      printf("mismatch-covered %d [%d,%d] %s %d %s %d\n",i,seen_one,seen_two,
		     one->match[one_pos]->factor->name,one->match[one_pos]->start,
		     two->match[two_pos]->factor->name,two->match[two_pos]->start );
	      */
	    }
	  } else {
	    if( one->match[one_pos]->strand == two->match[two_pos]->strand ) {
	      printf("%s %.2f SAME %d %s %d %.2f %.2f %s %d %.2f %.2f\n",
		     (is_gapped == TRUE ? "GAPPED" : "SOLID"), perc_id,		   
		     i,
		     one->match[one_pos]->factor->name,one->match[one_pos]->start,
		     one->match[one_pos]->bit_score,(one->match[one_pos]->bit_score - Probability2Bits(one->match[one_pos]->factor->max_prob)),
		     two->match[two_pos]->factor->name,two->match[two_pos]->start,
		     two->match[two_pos]->bit_score,(two->match[two_pos]->bit_score - Probability2Bits(two->match[two_pos]->factor->max_prob))
		     
		     );
	    }
	  }
	}
      }
      
      
    }
    
  }
  
}
