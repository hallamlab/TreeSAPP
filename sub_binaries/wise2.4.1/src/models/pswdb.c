/* 
   pswdb.c
   
   database smith/waterman searching. 
   will also include abc model.
   pre-alpha

   Modelled on estwisedb.c

   Richard Copley. 6/12/1998
*/

#include "dyna.h"
#include "abc.h"
#include "pba.h"
#include "sw_wrap.h"
#include "seqaligndisplay.h"
#include "version.h"
#include <string.h>
#include <math.h>

char * program_name = "pswdb";

CompMat * mat           = NULL;
DBSearchImpl * dbsi     = NULL;
Hscore * hs             = NULL;
DPRunImpl * dpri        = NULL;

ProteinDB * tdb = NULL;
ProteinDB * qdb = NULL;

int gap = 12;
int ext = 2;
int a = 120;
int b = 10; 
int c = 3;
int report_stagger = 100;
int aln_number = 500;
int max_desc = 500;

char * querydb             = NULL;
char * targetdb            = NULL;
char * matrix_file      = "BLOSUM62.bla";
char * output_file 	= "-";

char * cutoff_str    = NULL;
double search_cutoff = 40;
double prob_cutoff = 10;

boolean show_histogram = TRUE;
boolean show_ids = FALSE;
boolean show_pretty = TRUE;
boolean use_abc = FALSE;
boolean use_pba = FALSE;
boolean use_query_db = FALSE;
boolean fit_histogram = TRUE;

/* makes an alignment anchored to the query sequence */
boolean make_anchored_alignment = FALSE; 

FILE * ofp = NULL;

void show_help(FILE * ofp)
{
   fprintf(ofp,"\npswdb <options> <query_db> <target_db>\nSeqs in fasta format\n"
          "\t-g gap penalty (default 12)\n"
          "\t-e ext penatly (default 2)\n"
          "\t-m comp matrix (default BLOSUM62.bla)\n"
          "\t-abc use the abc model\n"
          "\t-a   a penalty for above (default 120)\n"
          "\t-b   b penalty for above (default 10)\n"
          "\t-c   c penalty for above (default 3)\n"
          "\t-pba use the pba model\n"
          "\t-max_desc Number of one line descriptions (default = 500)\n"
          "\t-max_aln Number of alignments to show (default = 50)\n"
	  "\t-cut Search cutoff (score) (default = 40)\n"
          "\t-ids Output seq ids with alignments\n"
	  "\t-nohis do not fit histogram\n"
          );

   /* show dbsearch impl help */
   show_help_DBSearchImpl(ofp);

   /* show run time */
   show_help_DPRunImpl(ofp);

   /* standard options */
  fprintf(stdout,"Standard options\n");
  fprintf(stdout,"  -help      help\n  -version   show version and compile info\n  -silent    No messages on stderr\n  -quiet    No report on stderr\n  -erroroffstd No warning messages to stderr\n  -errorlog [file] Log warning messages to file\n");
  fprintf(stdout,"\nSee WWW help at http://www.somewhere.edu/\n");
  
}


void show_version(void)
{
  fprintf(stdout,"%s\nVersion: %s\nReleased: %s\nCompiled: %s\n",program_name,VERSION_NUMBER,RELEASE_DAY,COMPILE_DATE);
  fprintf(stdout,"\nThis program is freely distributed under a Gnu Public License\n");
  fprintf(stdout,"The source code is copyright (c) EMBL 1998 and others\n");
  fprintf(stdout,"There is no warranty, implied or otherwise on the performance of this program\n");
  fprintf(stdout,"For more information read the GNULICENSE file in the distribution\n\n");
  fprintf(stdout,"Credits: Richard Copley (richard.copley@embl-heidelberg.de) wrote the core code.\n");
  fprintf(stdout,"         Portions of this code was written by Ewan Birney\n");
  fprintf(stdout,"         Portions of this code was from HMMer2, written by Sean Eddy\n");
  fprintf(stdout,"\nFor more information on the abc algorithm, see: \n");
  fprintf(stdout,"       S. F. Altschul - Proteins 32:88-96 (1998)\n\n");
  exit(63);   
}

boolean show_output(void)
{

   double pcid;

   boolean ret = TRUE;

   int i;

   int length,nid;
   int qmatch;

   char aa1,aa2;

   Protein * qseq;
   Protein * tseq;

   AlnBlock * alb;
   AlnColumn * alc;

   /* sort by bit score */

   sort_Hscore_by_score(hs);

/*
   if( make_anchored_alignment == TRUE ) {
     anc = single_unit_AlnBlock(qseq->len,"SEQUENCE");
     set = SequenceSet_alloc_std();
   }
*/



   if( fit_histogram == FALSE || hs->his == NULL || hs->his->total < 1000 ) {

      info("Cannot fit histogram to a db smaller than 1,000");
      fprintf(ofp,"[Warning: Can't fit histogram to a db smaller than 1,000]\n\n");
      show_histogram = FALSE;

   } else {

      if( use_abc ){
         fit_Hscore_to_EVD(hs,3000);
      } else  {
         fit_Hscore_to_EVD(hs,300);
      }

   }

   fprintf(ofp,"\n#High Score list\n");
   fprintf(ofp,"-------------------------------------------------------------\n");

   for(i=0;i<hs->len;i++) {
      if( (hs->ds[i]->score < search_cutoff) || (i >= max_desc) || (hs->ds[i]->evalue > prob_cutoff)){
         aln_number = i;
         break;
      }
      fprintf(ofp,"Probe %-24s Dbase %-24s %c %4d %.2g\n",hs->ds[i]->query->name,hs->ds[i]->target->name,hs->ds[i]->target->is_reversed == TRUE ? '-' : '+',hs->ds[i]->score,hs->ds[i]->evalue);

   }


   fprintf(ofp,"\n\n#Alignments\n");
   fprintf(ofp,"-------------------------------------------------------------\n\n");

   for(i=0;i<hs->len;i++) {
      if( hs->ds[i]->score < search_cutoff ) {
         break;
      }


      if( i >= aln_number ) {
         break;
      }
       
      fprintf(ofp,"\n\n>%s vs %s [%d]\n",hs->ds[i]->query->name,hs->ds[i]->target->name,i+1 );

      qseq = get_Protein_from_ProteinDB(qdb,hs->ds[i]->query);
      tseq = get_Protein_from_ProteinDB(tdb,hs->ds[i]->target);

      if( use_abc ) {
         alb = Align_Proteins_ABC(qseq,tseq,mat,-a,-b,-c,NULL,dpri);
      } else { 
         alb = Align_Proteins_SmithWaterman(qseq,tseq,mat,-gap,-ext,NULL,dpri);
      }

      length = 0;
      nid = 0;
      qmatch = 0;
      for( alc = alb->start; alc != NULL; alc = alc->next) {
         length++;
         aa1 = Protein_seqchar(qseq,1+alc->alu[0]->start);
         aa2 = Protein_seqchar(tseq,1+alc->alu[1]->start);
         if(!strcmp(alc->alu[0]->text_label,"SEQUENCE")){
            qmatch++;
         }
         if((!strcmp(alc->alu[0]->text_label,"SEQUENCE")) &&
            (!strcmp(alc->alu[1]->text_label,"SEQUENCE")) &&
            (aa1 == aa2)){
            nid++;
         }
      }
      length--;
      pcid = 100.00 * ((double)nid/(double)length);


/*
      if( make_anchored_alignment == TRUE ) {
	alb->als[1]->data = (void *) tseq;
	add_SequenceSet(seqset,tseq);
	add_to_anchored_AlnBlock(anc,alb);
      }
*/

	
      if( show_pretty == TRUE ) {
         fprintf(ofp,"\nP: %.2g Sc: %3d aln_len: %d p_len: %d t_len: %d match: %4.1f%% ID: %4.1f%%\n\n",hs->ds[i]->evalue,hs->ds[i]->score,length,Protein_length(qseq),Protein_length(tseq),100.00*((double)qmatch/(double)(Protein_length(qseq))),pcid);
         if( show_ids == FALSE ) {
	   ckfree(qseq->baseseq->name);
	   qseq->baseseq->name = stringalloc("probe:");

	   ckfree(tseq->baseseq->name);
	   tseq->baseseq->name = stringalloc("target:");

	   /*** icky code. ***/
	   /*
            strcpy(qseq->baseseq->name,"probe:");
            strcpy(tseq->baseseq->name,"target:");
	   */
         }
         write_pretty_Protein_align(alb,qseq,tseq,15,50,stdout);
      }

      fprintf(ofp,"\n-------------------------------------------------------------\n\n");

      free_Protein(qseq);
      free_Protein(tseq);
      free_AlnBlock(alb);

   }


   if( show_histogram == TRUE ) {
      PrintASCIIHistogram(hs->his,ofp);
   }

   return ret;

}

boolean free_objects(void)
{

   if( tdb != NULL ) {
      tdb = free_ProteinDB(tdb);
   }
   if( qdb != NULL ) {
      qdb = free_ProteinDB(qdb);
   }
   if( hs != NULL ) {
      hs = free_Hscore(hs);
   }

   return TRUE;
}

boolean build_objects(void)
{

   boolean ret = TRUE;
   Sequence * qseq;

   if( (mat = read_Blast_file_CompMat(matrix_file)) == NULL) {
      warn("Could not read Comparison matrix file in %s",matrix_file);
      ret = FALSE;
   }


   if((use_abc == TRUE) && (strcmp(matrix_file,"abc.bla"))) {
      factor_CompMat(mat,10);
   }

   if( use_query_db ) { 
      if((qdb = single_fasta_ProteinDB(querydb)) == NULL ) {
         ret = FALSE;
         warn("Could not read Protein sequence in %s",querydb);
      }
   } else {
      if((qseq = read_fasta_file_Sequence(querydb)) == NULL ) {
         ret = FALSE;
         warn("Could not read single Protein sequence in %s",querydb);
      }
      qdb = new_ProteinDB_from_single_seq(qseq);

      /** free actual sequence now. It is held onto by the database */
      free_Sequence(qseq);
   }

   if((tdb = single_fasta_ProteinDB(targetdb)) == NULL ) {
      ret = FALSE;
      warn("Could not read Protein sequence in %s",targetdb);
   }

   if( (ofp = openfile(output_file,"W")) ==  NULL) {
      warn("Could not open %s as an output file",output_file);
      ret = FALSE;
   }

   return ret;

}

boolean search_db(void)
{

   info("Starting search....");

   if( use_abc ) {
      hs = Hscore_from_ProteinABC(qdb,tdb,mat,-a,-b,-c,search_cutoff,report_stagger,FALSE,dbsi);
   } else if( use_pba){
/*
      hs = Hscore_from_ProteinBA(qdb,tdb,mat,bentry,bexit,bfor_trans,b_self_trans,b3exit,search_cutoff,report_stgger,dbsi);
*/
   } else  {
      hs = Hscore_from_ProteinSW(qdb,tdb,mat,-gap,-ext,search_cutoff,report_stagger,FALSE,dbsi);
   }

   if( hs == NULL ) {
      return FALSE;
   }

   return TRUE;
}

boolean show_header(FILE * ofp)

{
  fprintf(ofp,"-------------------------------------------------------------\n");
  fprintf(ofp,"Wise2 - Protein vs. Protein\n");
  fprintf(ofp,"Program: %s version: %s released: %s\n",program_name,VERSION_NUMBER,RELEASE_DAY);
  fprintf(ofp,"This program is freely distributed under a Gnu Public License.\n");
  fprintf(ofp,"   See -version for more info on copyright\n");
  fprintf(ofp,"Bugs and credits to: Richard Copley <copley@embl-heidelberg.de>\n"); 
  fprintf(ofp,"                     Ewan Birney <birney@sanger.ac.uk>\n");
  fprintf(ofp,"-------------------------------------------------------------\n\n");
  if( use_abc ){
     fprintf(ofp,"Algorithm type: abc model\n");
     fprintf(ofp,"a:  %d\n",a);
     fprintf(ofp,"b:  %d\n",b);
     fprintf(ofp,"c:  %d\n",c);
  } else { 
     fprintf(ofp,"Algorithm type: Smith/Waterman\n");
     fprintf(ofp,"Gap open:       %d\n",gap);
     fprintf(ofp,"Gap extension:  %d\n",ext);
  }
  fprintf(ofp,"Query info from:     %s\n",querydb);
  fprintf(ofp,"Database info from:  %s\n",targetdb);
  fprintf(ofp,"Comp Matrix:         %s\n",matrix_file);

  return TRUE;
}

int main(int argc,char *argv[])
{
  char * errlog;

   if( strip_out_boolean_argument(&argc,argv,"h") == TRUE || strip_out_boolean_argument(&argc,argv,"-help") == TRUE) {
      show_help(stdout); 
      exit(1);
   }

   (void) strip_out_integer_argument(&argc,argv,"g",&gap);
   (void) strip_out_integer_argument(&argc,argv,"e",&ext);
   (void) strip_out_integer_argument(&argc,argv,"a",&a);
   (void) strip_out_integer_argument(&argc,argv,"b",&b);
   (void) strip_out_integer_argument(&argc,argv,"c",&c);
   (void) strip_out_integer_argument(&argc,argv,"max_desc",&max_desc);
   (void) strip_out_integer_argument(&argc,argv,"max_aln",&aln_number);
   (void) strip_out_float_argument(&argc,argv,"cut",&search_cutoff);
   use_abc = strip_out_boolean_argument(&argc,argv,"abc"); 
   use_pba = strip_out_boolean_argument(&argc,argv,"pba");
   show_ids = strip_out_boolean_argument(&argc,argv,"ids");
   use_query_db = strip_out_boolean_argument(&argc,argv,"db_vs_db");
   
   if( strip_out_boolean_argument(&argc,argv,"nohis") == TRUE ) {
     fit_histogram = FALSE;
   }

   matrix_file = strip_out_assigned_argument(&argc,argv,"m");
   if( matrix_file == NULL)
      matrix_file = "BLOSUM62.bla";

   /* database implementation stuff */
   dbsi = new_DBSearchImpl_from_argv(&argc,argv);

   /* run time */

   dpri = new_DPRunImpl_from_argv(&argc,argv);

   
   /* standard options */

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
   

   /*** processed all arguments. Should be left with sequences in argv[1] and argv[2] */

   if( argc != 3 ) {
     warn("Must have two arguments for sequence 1 and sequence 2 %d",argc);
     show_help(stdout);
     exit(1);
   }

   if(use_abc && use_pba){
     warn("Choose only one of pba or abc");
     show_help(stdout);
     exit(1);
   }

   if( use_pba ){
      warn("PBA algorithm not yet working");
      exit(1);
   }

   querydb = argv[1];
   targetdb = argv[2];

   if( build_objects() == FALSE)
      fatal("Could not build objects!");

   show_header(stdout); 

   if( search_db() == FALSE)
      warn("Could not search database");

   show_output(); 

   free_objects(); 


   return 0;
}
