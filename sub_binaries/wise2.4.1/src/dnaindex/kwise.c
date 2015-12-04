#include "../models/version.h"
#include "assembly_stream_cons.h"
#include "kmer_direct.h"
#include "kmer_hash.h"
#include "kmer_glib_index.h"

#include "kmer_assembly_untangler.h"
#include "kmer_assembly_error.h"
#include "kmer_assembly_contig.h"
#include "assembly.h"

char * program_name = "kwise";


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
  fprintf(ofp,"%s <arguments>\n",program_name);
  fprintf(ofp,"  assembles sequences into contigs\n");
  fprintf(ofp,"Kmer options\n");
  fprintf(ofp,"  -kmer [12]              kmer size for kmer-graph\n");
  fprintf(ofp,"  -[no]use_hash [default no] use 28bit hash for large kmer sizes\n");
  fprintf(ofp,"  -[no]use_glib [default no] use glib based pure hash, up to 15mers\n");
  fprintf(ofp,"Read entry control\n");
  fprintf(ofp,"  -read_no [number]       maximum number of reads to enter (useful for debugging)\n");
  show_help_AssemblyStreamConstructor(ofp);

  fprintf(ofp,"Statistics and Logging options\n");
  fprintf(ofp,"  -stats <filename>  write kmer index stats to file\n");

  show_standard_options(ofp);
}


int main(int argc,char ** argv)
{
  KmerAssemblyContigPara * p;

  KmerIndexInterface * kii;
  KmerAssemblyIndex * kai;

  int kmer_size = 12;
  boolean use_28bithash = 0;
  boolean use_glib = 0;
  int read_no = 0;
  long int r = 0;

  boolean add_rev = 1;
  boolean remove_error = 1;

  char * stats_file;
  FILE * stats = NULL;

  SinglePosSpace * sps;

  AssemblyStreamConstructor * acs;
  AssemblySequenceStream * ass;
  AssemblySequence * aseq;
  AssemblySequence * mirror;
 
  Assembly * assembly;

  char * fasta_read = NULL;
  FILE * fasta_read_fp = NULL;

  /* no options, and have to have at least one option */
  if( argc == 1 ) {
    show_help(stdout);
    exit(12);
  }


  acs = new_AssemblyStreamConstructor_from_argv(&argc,argv);

  strip_out_integer_argument(&argc,argv,"kmer",&kmer_size);

  strip_out_boolean_def_argument(&argc,argv,"use_hash",&use_28bithash);
  strip_out_boolean_def_argument(&argc,argv,"use_glib",&use_glib);

  strip_out_boolean_def_argument(&argc,argv,"both",&add_rev);
  strip_out_boolean_def_argument(&argc,argv,"error",&remove_error);

  stats_file = strip_out_assigned_argument(&argc,argv,"stats");

  fasta_read = strip_out_assigned_argument(&argc,argv,"fasta_read");

  strip_out_integer_argument(&argc,argv,"read_no",&read_no);

  strip_out_standard_options(&argc,argv,show_help,show_version);


  if( argc != 1 ) {
    show_help(stdout);
    exit(12);
  }

  if( fasta_read != NULL ) {
    fasta_read_fp = openfile(fasta_read,"w");
    if( fasta_read_fp == NULL ) {
      fatal("Asked for a read stream in %s, but unable to open as read file...",fasta_read);
    }
  }

  p = KmerAssemblyContigPara_alloc();
  p->minimum_len = 0;
  p->minimum_depth = 1;

  ass = new_AssemblySequenceStream_from_AssemblyStreamConstructor(acs);

  if( stats_file != NULL ) {
    stats = openfile(stats_file,"w");
    if( stats == NULL ) {
      fatal("Could not open stats file [%s]",stats_file);
    }
  }


  if( use_28bithash == TRUE )  {
    kii = new_KmerHashIndex(kmer_size);
  } else if( use_glib == TRUE ) {
    kii = new_interface_KmerGlibIndex(kmer_size);
  } else {
    kii = new_KmerDirectIndex(kmer_size);
  }

  sps = new_SinglePosSpace(1,2000);

  kai = new_KmerAssemblyIndex(kii,sps);

  while( (aseq = (*ass->next_AssemblySequence)(ass->handle)) ) {
    fprintf(stderr,"Loading in %s\n",aseq->seq->name);
    if( fasta_read_fp != NULL ) {
      write_fasta_Sequence(aseq->seq,fasta_read_fp);
    }

    add_AssemblySequence_KmerAssemblyIndex(kai,aseq,1000000);
    if( add_rev == 1 ) {
      mirror = mirrored_AssemblySequence(aseq);
      add_AssemblySequence_KmerAssemblyIndex(kai,mirror,1000000); 
    }

    r++;
    if( read_no > 0 && r > read_no) {
      break;
    }
  }

  
  if( remove_error == 1 ) {
    resolve_forward_errors_KmerAssembly(kai,1,1,250); 
    remove_errors_KmerAssemblyIndex(kai,1);
  }

  untangle_KmerAssembly(kai);

  mark_tangles_KmerAssembly(kai);  

  assembly = Assembly_from_KmerAssemblyIndex(kai,p);

  if( stats != NULL) 
    show_extensive_stats_KmerAssemblyIndex(kai,stats);


  dump_contigs_as_fasta_Assembly(assembly,stdout);

  return 0;
}

