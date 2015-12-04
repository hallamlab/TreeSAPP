
#include "dnaalign.h"
#include "dnamatcher.h"
#include "version.h"



char * program_name = "amplimer_resolver";

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
  fprintf(ofp,"%s target-fasta-file amplimer-fasta-file\n",program_name);

  show_help_DnaMatchPara(ofp);

  show_help_HitListOutputImpl(ofp);

  show_standard_options(ofp);
}


Sequence * swapped_Sequence(Sequence * reference,HitList * hl)
{
  Sequence * out;
  int i;
  AlnColumn * alc;
  char buffer[512];
  char * swapped;
  int coord;
  int swap;
 
  sprintf(buffer,"swapped_%s",reference->name);

  out = new_Sequence_from_strings(buffer,reference->seq);

  for(i=0;i<reference->len;i++) {
    out->seq[i] = tolower(reference->seq[i]);
  }

  /* calloc sets things to 0 be defn */
  swapped = calloc(reference->len,sizeof(char));
  
  for(i=0;i<hl->len;i++) {
    fprintf(stderr,"Using... %s\n",hl->pair[i]->target->name);

    for(alc =hl->pair[i]->aln[0]->alb->start ;alc != NULL;alc = alc->next) {
      /*      fprintf(stderr," Looking at %d,%d\n",alc->alu[0]->end,alc->alu[1]->end);*/

      if( strcmp(alc->alu[0]->text_label,"SEQUENCE") == 0 &&
	  strcmp(alc->alu[1]->text_label,"SEQUENCE") == 0 ) {
	/*	fprintf(stderr,"MATCH Looking at %d,%d\n",alc->alu[0]->end,alc->alu[1]->end); */
	coord = alc->alu[0]->end;
	swap = hl->pair[i]->target->seq[alc->alu[1]->end];

	if( swapped[coord] == 1 ) {
	  warn("Reswapping position %d to %c from %c",coord,swap,out->seq[coord]);
	} else if( toupper(swap) != toupper(out->seq[coord]) ) {
	  info("Swapping position %d from %c to %c",coord,out->seq[coord],swap);
	}

    
	out->seq[coord] = toupper(swap);
	swapped[coord] = 1;
      }
    }
  }


  ckfree(swapped);
  return out;
}

Sequence * translate_swapped(Sequence * swapped) 
{
  CodonTable * ct;
  int i,j;
  Sequence * out;
  
  out = Sequence_alloc();
  out->name = stringalloc(swapped->name);
  out->seq = calloc(1+swapped->len/3,sizeof(char));

  ct = read_CodonTable_file("codon.table");

  for(i=0,j=0;i<swapped->len;i+=3,j++) {
    out->seq[j] = aminoacid_from_seq(ct,swapped->seq+i);
    if( isupper(swapped->seq[i]) && isupper(swapped->seq[i+1]) &&
	isupper(swapped->seq[i+2]) ) {
      out->seq[j] = toupper(out->seq[j]);
    } else{
      out->seq[j] = tolower(out->seq[j]);
    }
  }

  out->seq[j] = '\0';

  return out;

}


int main(int argc,char ** argv)
{
  DnaMatchPara * para;
  HitListOutputImpl * hitoutput;
  HitList * hitlist;

  Sequence * reference;
  Sequence * swap;
  Sequence * trans;
  SequenceSet * amplimers;
  int show_hitlist = 0;
  int show_swapped = 1;

  hitoutput = new_HitListOutputImpl_from_argv(&argc,argv);
  para = new_DnaMatchPara_from_argv(&argc,argv);

  strip_out_boolean_def_argument(&argc,argv,"hitlist",&show_hitlist);
  strip_out_boolean_def_argument(&argc,argv,"swapped",&show_hitlist);

  strip_out_standard_options(&argc,argv,show_help,show_version);
  if( argc != 3 ) {
    show_help(stdout);
    exit(12);
  }

  reference = read_fasta_file_Sequence(argv[1]);
  amplimers = read_fasta_SequenceSet_file(argv[2]);

  hitlist = HitList_from_Sequence_SequenceSet_DNA(reference,amplimers,para);
  
  if( show_hitlist ) {
    show_HitList_HitListOutputImpl(hitoutput,hitlist,stdout);
  }

  swap  = swapped_Sequence(reference,hitlist);
  trans = translate_swapped(swap);


  if( show_swapped ) {
    write_fasta_Sequence(swap,stdout);
    write_fasta_Sequence(trans,stdout);
  }

}
