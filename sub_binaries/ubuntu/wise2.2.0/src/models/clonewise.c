
#include "clonewisedp.h"
#include "localclonewisedp.h"
#include "version.h"


/*
 * Yup - I know these are globals. It just ends up being
 * too much of a pain in the arse to coordinate the help printing
 * defaults with the parameters by passing around a void * structure
 * which you cast back just to prevent no globals. You want to do this -
 * feel free to send me a patch to this and also wisecommandline.[ch] and
 * port all the other wise main functions to use it!
 */

char * program_name = "genomewise";
int query_gap_start  = 2;
int query_gap_extend = 1;
int target_gap_start  = 2;
int target_gap_extend = 1;
int query_switch_cost = 10;
int match_score      = 10;
int mismatch_score   = -1;

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
  fprintf(ofp,"%s trusted-clone-path weak-clone-path\n",program_name);
  fprintf(ofp,"Algorithm\n");
  /*  fprintf(ofp,"   -spread          smell area to use in comparing clones\n"); */
  fprintf(ofp,"   -match           [%d] clone match score\n",match_score);
  fprintf(ofp,"   -mismatch        [%d] clone mismatch score\n",mismatch_score);
  fprintf(ofp,"   -wgap            [%d] weak    gap start cost\n",query_gap_start);
  fprintf(ofp,"   -wext            [%d] weak    gap extension\n",query_gap_extend);
  fprintf(ofp,"   -wswitch         [%d] weak    switch cost\n",query_switch_cost);
  fprintf(ofp,"   -tgap            [%d] trusted gap start cost\n",target_gap_start);
  fprintf(ofp,"   -text            [%d] trusted gap extension\n",target_gap_extend);

  fprintf(ofp,"Output\n");
  fprintf(ofp,"   -[no]path        show path format (default no)\n");
  fprintf(ofp,"   -[no]zip         show zip format (default yes)\n");
  fprintf(ofp,"   -[no]alb         show aln block format (default no)\n");
  fprintf(ofp,"   -[no]pal         show packaln format (default no)\n");

  show_standard_options(ofp);
}

double id(int i)
{
  return (double)i;
}


void debug_zip(AlnBlock * alb,MappedCloneSet * query,MappedCloneSet * target,FILE * ofp)
{
  AlnColumn * alc;
  MappedCloneSet * temp;
  int i;

  /* descend each block, printing left and right info */

  for(alc=alb->start;alc != NULL;alc = alc->next) {
    fprintf(ofp,"%4d %4d:%4d %12s:%12s -- ",alc->alu[0]->score[0],alc->alu[0]->end,alc->alu[1]->end,alc->alu[0]->text_label,alc->alu[1]->text_label);

    if( strstr(alc->alu[0]->text_label,"SKIP") != NULL || strstr(alc->alu[0]->text_label,"MATCH") != NULL ) {
      fprintf(ofp,"[");
      temp = subsection_MappedCloneSet(query,alc->alu[0]->end,alc->alu[0]->end,0);
      if( temp->len == 0 ) {
        fprintf(ofp,"]");
      }
      for(i=0;i<temp->len;i++) {
	fprintf(ofp,"%s%c",temp->clone[i]->clone_name,i+1 == temp->len ? ']' : ',');
      }
      free_MappedCloneSet(temp);
    } else {
      fprintf(ofp,"[]");
    }

    if( strstr(alc->alu[1]->text_label,"SKIP") != NULL || strstr(alc->alu[1]->text_label,"MATCH") != NULL ) {
      fprintf(ofp,"[");
      temp = subsection_MappedCloneSet(target,alc->alu[1]->end,alc->alu[1]->end,0);
      if( temp->len == 0 ) {
        fprintf(ofp,"]");
      }
      for(i=0;i<temp->len;i++) {
	fprintf(ofp,"%s%c",temp->clone[i]->clone_name,i+1 == temp->len ? ']' : ',');
      }
      free_MappedCloneSet(temp);
    } else {
      fprintf(ofp,"[]");
    }
    fprintf(ofp,"\n");
  }

}

void extended_path(AlnBlock * alb,MappedCloneSet * query,MappedCloneSet * target,FILE * ofp)
{
  AlnColumn * alc;
  AlnColumn * runner;

  int coord;
  int coord_runner;
  int startpos;
  int i;
  int j;
  int is_new;
  

  MappedCloneSet * temp;
  MappedCloneSet * seen_query;
  MappedCloneSet * seen_target;

  seen_query  = MappedCloneSet_alloc_std();
  seen_target = MappedCloneSet_alloc_std();


  coord = 0;
  for(alc=alb->start;alc != NULL;alc = alc->next) {

    /* deal with weak first */
    
    temp = subsection_MappedCloneSet(query,alc->alu[0]->end,alc->alu[0]->end,0);
    for(i=0;i<temp->len;i++) {
      is_new = 1;
      for(j=0;j<seen_query->len;j++) {
	if( temp->clone[i] == seen_query->clone[j] ) {
	  is_new = 0;
	  break;
	}
      }
      if( is_new == 0 ) {
	/* no need to print it out */
	continue;
      }

      add_MappedCloneSet(seen_query,hard_link_MappedClone(temp->clone[i]));

      /* time to print it out */
      startpos = coord;
      coord_runner = coord;

      
      for(runner = alc;runner != NULL && runner->alu[0]->end <= temp->clone[i]->end;runner=runner->next)
	coord_runner++;

      if( runner == NULL ) {
	continue;
      }
      
      fprintf(ofp,"QUERY\t%d\t%d\t%s\t%s\t%s\t%d\t%d\n",startpos,coord_runner,temp->clone[i]->clone_name,temp->clone[i]->accession,temp->clone[i]->contig,alc->alu[1]->end,runner->alu[1]->end);
    }
    free_MappedCloneSet(temp);
    /* Ooops. copy-and-paste. Bad Ewan! */

    temp = subsection_MappedCloneSet(target,alc->alu[1]->end,alc->alu[1]->end,0);
    for(i=0;i<temp->len;i++) {
      is_new = 1;
      for(j=0;j<seen_target->len;j++) {
	if( temp->clone[i] == seen_target->clone[j] ) {
	  is_new = 0;
	  break;
	}
      }
      if( is_new == 0 ) {
	/* no need to print it out */
	continue;
      }

      add_MappedCloneSet(seen_target,hard_link_MappedClone(temp->clone[i]));

      /* time to print it out */
      startpos = coord;
      coord_runner = coord;
      for(runner = alc;runner != NULL && runner->alu[1]->end <= temp->clone[i]->end;runner=runner->next)
	coord_runner++;
      /* if the end is not there - skip this one? */

      if( runner == NULL ) {
	continue;
      }

      fprintf(ofp,"TARGET\t%d\t%d\t%s\t%s\t%s\t%d\t%d\n",startpos,coord_runner,temp->clone[i]->clone_name,temp->clone[i]->accession,temp->clone[i]->contig,alc->alu[1]->end,runner->alu[1]->end);
    }
    free_MappedCloneSet(temp);
    
    coord++;
  }
}



int main (int argc,char ** argv)
{
  MappedCloneSet * trusted;
  MappedCloneSet * weak;
  MappedCloneMatch * match;


  FILE * in;
  int kbyte = 10000;
  PackAln  * pal;
  AlnBlock * alb;

  int spread = 30;

  boolean show_alb = 0;
  boolean show_pal = 0;
  boolean show_zip = 1;
  boolean show_path = 0;

  char * alg_string = "local";
  char * temp;
  char * divide_string = "//";

  strip_out_boolean_def_argument(&argc,argv,"alb",&show_alb);
  strip_out_boolean_def_argument(&argc,argv,"pal",&show_pal);
  strip_out_boolean_def_argument(&argc,argv,"zip",&show_zip);
  strip_out_boolean_def_argument(&argc,argv,"path",&show_path);

  strip_out_integer_argument(&argc,argv,"wgap",&query_gap_start);
  strip_out_integer_argument(&argc,argv,"wext",&query_gap_extend);
  strip_out_integer_argument(&argc,argv,"wswitch",&query_switch_cost);
  strip_out_integer_argument(&argc,argv,"tgap",&target_gap_start);
  strip_out_integer_argument(&argc,argv,"text",&target_gap_extend);
  strip_out_integer_argument(&argc,argv,"match",&match_score);
  strip_out_integer_argument(&argc,argv,"mismatch",&mismatch_score);
  temp =strip_out_assigned_argument(&argc,argv,"alg");
  if( temp != NULL ) {
    alg_string = temp;
  }

  strip_out_integer_argument(&argc,argv,"spread",&spread);
  strip_out_integer_argument(&argc,argv,"kbyte",&kbyte);

  strip_out_standard_options(&argc,argv,show_help,show_version);
  if( argc != 3 ) {
    show_help(stdout); 
    exit(12);
  }

  in = fopen(argv[1],"r");
  if( in == NULL ) {
    fatal("Unable to open %s",argv[1]);
  }
  trusted = read_MappedCloneSet(in);

  /*  fprintf(stderr,"first start %d\n",trusted->clone[0]->start);*/

  in = fopen(argv[2],"r");
  if( in == NULL ) {
    fatal("Unable to open %s",argv[2]);
  }

  change_max_BaseMatrix_kbytes(kbyte);

  weak = read_MappedCloneSet(in);

  synchronise_MappedCloneSets(trusted,weak);

  /*  fprintf(stderr,"score for 2,2 is %d\n",MappedCloneSet_match(weak,trusted,2,2,0,10,-5)); */
  
  match = new_MappedCloneMatch(weak,trusted,match_score,mismatch_score);

  fprintf(stderr,"Match matrix calculated\n");

  if( strcmp(alg_string,"global") == 0 ) {
    pal = PackAln_bestmemory_CloneWise(weak,trusted,match,-query_gap_start,-query_gap_extend,-target_gap_start,-target_gap_extend,spread,-query_switch_cost,NULL);
    alb = convert_PackAln_to_AlnBlock_CloneWise(pal);
  } else if ( strcmp(alg_string,"local") == 0 ) {
    pal = PackAln_bestmemory_LocalCloneWise(weak,trusted,match,-query_gap_start,-query_gap_extend,-target_gap_start,-target_gap_extend,spread,-query_switch_cost,NULL);
    alb = convert_PackAln_to_AlnBlock_LocalCloneWise(pal);
  } else {
    /* keep gcc happy */
    pal = NULL;
    alb = NULL;
    fatal("Not a proper algorithm string %s",alg_string);
  }

  if( show_path ) {
    extended_path(alb,weak,trusted,stdout);
    fprintf(stdout,"%s\n",divide_string);
  }

  if( show_zip ) {
    debug_zip(alb,weak,trusted,stdout);
    fprintf(stdout,"%s\n",divide_string);
  }

  if( show_alb ) {
    mapped_ascii_AlnBlock(alb,id,1,stdout);
    fprintf(stdout,"%s\n",divide_string);
  }

  if( show_pal ) {
    show_simple_PackAln(pal,stdout);
    fprintf(stdout,"%s\n",divide_string);
  }

  return 0;
}




