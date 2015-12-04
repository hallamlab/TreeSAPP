
/* has to be before the others due to nasty namespace clashes */
#define WISE2_CROSS_HMMER2
#include "wise2xhmmer2.h"


#include "dyna.h"
#include "version.h" 

char * program_name = "estwise";

/*
 * program specific includes
 */

#include "estwrap.h"
#include "genedisplay.h"
#include "matchsum.h"


/*
 * program specific variables
 */

char * dna_seq_file  = NULL;
cDNA * cdna          = NULL;

char * tstart_str    = NULL;
int tstart           = -1;

char * tend_str      = NULL;
int tend             = -1;
boolean reverse      = FALSE;
boolean target_abs         = FALSE;

boolean use_tsm      = FALSE;


char * protein_file  = NULL;
Protein * pro        = NULL;

char * hmm_file       = NULL;
ThreeStateModel * tsm = NULL;
char * hmm_name       = NULL;

char * qstart_str    = NULL;
int qstart           = -1;

char * qend_str      = NULL;
int qend             = -1;


char * matrix_file      = NULL;
CompMat * mat           = NULL;

char * gap_str          = "12";
int gap                 = 12;

char * ext_str          = "2";
int ext                 = 2;

char * codon_file       = NULL;
CodonTable * ct         = NULL;

char * output_file      = "-";
FILE * ofp              = NULL;

RandomModelDNA * rmd    = NULL;

char * subs_string      = "0.01";
double subs_error       = 0.01;

char * indel_string      = "0.01";
double indel_error       = 0.01;


char * null_string       = "syn";
boolean use_syn          = FALSE;
boolean force_flat_insert= FALSE;
boolean should_set_hmm   = FALSE;


char * allN_string      = "1.0";
Probability allN        = 1.0;

char * startend_string   = "default";
int startend             = TSM_default;

int alg                  = ESTWISE_3;
char * alg_str           = NULL;


boolean show_PackAln     = FALSE;
boolean show_AlnBlock    = FALSE;
boolean show_pretty      = FALSE;
boolean show_pep         = FALSE;
boolean show_para        = FALSE;
boolean show_match_sum   = FALSE;

boolean do_both = TRUE;

char * main_block_str      = "50";
int main_block           = 50;

char * divide_str        = "//";

AlnBlock * alb;
PackAln  * pal;

cDNAParser * cp;
RandomModel * rm;
CodonMapper * cm;

MatchSummarySet * mss;

DPRunImpl * dpri = NULL;


void show_parameters(void)
{
  fprintf(ofp,"%s %s (%s release)\n",program_name,VERSION_NUMBER,RELEASE_DAY);
  fprintf(ofp,"This program is freely distributed under a GPL. See source directory\n");
  fprintf(ofp,"Copyright (c) GRL limited: portions of the code are from separate copyrights\n\n");
  if( use_tsm == FALSE) {
    fprintf(ofp,"Query protein:    %s\n",pro->baseseq->name);
    fprintf(ofp,"Comp Matrix:      %s\n",matrix_file);
    fprintf(ofp,"Gap open:         %d\n",gap);
    fprintf(ofp,"Gap extension:    %d\n",ext);
  }
  else 
    fprintf(ofp,"Query model:      %s\n",tsm->name);
  fprintf(ofp,"Start/End         %s\n",startend_string);
  fprintf(ofp,"Target Sequence   %s\n",cdna->baseseq->name);
  fprintf(ofp,"Strand:           %s\n",do_both == TRUE ? "both" : reverse == TRUE ? "reverse" : "forward");
  fprintf(ofp,"Codon Table:      %s\n",codon_file);
  fprintf(ofp,"Subs error:       %2.2g\n",subs_error);
  fprintf(ofp,"Indel error:      %2.2g\n",indel_error);
  fprintf(ofp,"Algorithm         %s\n",alg_str);
}

boolean show_output(void)
{
  Protein * ps;
  AlnColumn * alt;
  Protein * trans;

  if( show_pretty == TRUE ) {
    fprintf(ofp,"\n%s output\nScore %4.2f bits over entire alignment\n",program_name,Score2Bits(pal->score));
    if( use_syn == TRUE) {
      fprintf(ofp,"Score reported as bits over a synchronised coding sequence model\n\n");
    } else {
      fprintf(ofp,"Score reported as bits over a flat random model\n\n");
    }

    if( use_tsm == FALSE ) 
      protcdna_ascii_display(alb,pro->baseseq->seq,pro->baseseq->name,pro->baseseq->offset,cdna,ct,15,main_block,alg == ESTLOOP_3 ? TRUE : FALSE,ofp);
    else {
      ps = pseudo_Protein_from_ThreeStateModel(tsm);
      protcdna_ascii_display(alb,ps->baseseq->seq,ps->baseseq->name,ps->baseseq->offset,cdna,ct,15,main_block,alg == ESTLOOP_3 ? TRUE : FALSE,ofp);
      free_Protein(ps);
    }
    fprintf(ofp,"%s\n",divide_str);
  }

  if( show_match_sum == TRUE ) {
    show_MatchSummary_estwise_header(ofp);
    show_MatchSummarySet_estwise(mss,ofp);
    fprintf(ofp,"%s\n",divide_str);
  }


  if( show_pep == TRUE ) {
    alt = alb->start;
    for(;alt != NULL;) {
      trans = Protein_from_GeneWise_AlnColumn(cdna->baseseq,alt,1,&alt,ct,is_random_AlnColumn_genewise);
      if ( trans == NULL ) 
	break;
      write_fasta_Sequence(trans->baseseq,ofp);
      free_Protein(trans);
    }
    fprintf(ofp,"%s\n",divide_str);
  }

  if( show_AlnBlock == TRUE ) {
    mapped_ascii_AlnBlock(alb,Score2Bits,0,ofp);
    fprintf(ofp,"%s\n",divide_str);
  }

  if( show_PackAln == TRUE ) {
    show_simple_PackAln(pal,ofp);
    fprintf(ofp,"%s\n",divide_str);
  }
  
  return TRUE;
}


boolean free_temporary_objects(void)
{
  free_AlnBlock(alb);
  free_PackAln(pal);
  free_MatchSummarySet(mss);

  return TRUE;
}

boolean free_io_objects(void)
{
  if( use_tsm == TRUE) {
    free_ThreeStateModel(tsm);
  } else {
    free_Protein(pro);
  }

  free_CodonTable(ct);
  free_RandomModelDNA(rmd);
  free_cDNA(cdna);

  return TRUE;
}

void reverse_target(void)
{
  cDNA * cd_temp;

  cd_temp = cDNA_from_Sequence(reverse_complement_Sequence(cdna->baseseq));
  free_temporary_objects();
  free_cDNA(cdna);
  
  cdna = cd_temp;
}

boolean build_alignment(void)
{
  


  if( use_tsm == FALSE ) {
    alb =  AlnBlock_from_Protein_estwise_wrap(pro,cdna,cp,cm,ct,mat,-gap,-ext,startend == TSM_global ? 1 : 0,rmd,alg,rm,use_syn,allN,dpri,&pal);
  } else {
    alb =  AlnBlock_from_TSM_estwise_wrap(tsm,cdna,cp,cm,ct,rmd,alg,use_syn,force_flat_insert,allN,dpri,&pal);
  }

  if( alb == NULL )
    return FALSE;

  if( use_tsm ) {
    mss = MatchSummarySet_from_AlnBlock_estwise(alb,tsm->name,1,cdna->baseseq);
  } else {
    mss = MatchSummarySet_from_AlnBlock_estwise(alb,pro->baseseq->name,pro->baseseq->offset,cdna->baseseq);
  }

  return TRUE;
}


boolean build_objects(void)
{
  boolean ret = TRUE;
  Protein * pro_temp;
  cDNA    * cdna_temp;

  startend = threestatemodel_mode_from_string(startend_string);
  if( startend == TSM_unknown ) {
    warn("String %s was unable to converted into a start/end policy\n",startend_string);
    ret = FALSE;
  }


  if( tstart_str != NULL ) {
    if( is_integer_string(tstart_str,&tstart) == FALSE || tstart < 0) {
      warn("Could not make %s out as target start",tstart);
      ret = FALSE;
    }
  }

  if( is_integer_string(gap_str,&gap) == FALSE ) {
      warn("Could not make %s out as gap penalty (must be integer at the moment)",gap_str);
      ret = FALSE;
  }
  

  if( is_integer_string(ext_str,&ext) == FALSE ) {
    warn("Could not make %s out as gap penalty (must be integer at the moment)",ext_str);
    ret = FALSE;
  }
  

  if( tend_str != NULL ) {
    if( is_integer_string(tend_str,&tend) == FALSE || tend < 0) {
      warn("Could not make %s out as target end",tend);
      ret = FALSE;
    }
  }

  if( (cdna = read_fasta_file_cDNA(dna_seq_file)) == NULL ) {
    ret = FALSE;
    warn("Could not read cdnaomic sequence in %s",dna_seq_file);
  } else {
    if( tstart != -1 || tend != -1 ) {
      if( tstart == -1 )
	tstart = 0;
      if( tend == -1 ) 
	tend = cdna->baseseq->len;

      cdna_temp = truncate_cDNA(cdna,tstart,tend);
      if( cdna_temp == NULL ){
	ret = FALSE;
      } else {
	free_cDNA(cdna);
	cdna = cdna_temp;
      }
    }
  }

  if( reverse == TRUE ) {
    if( tstart > tend ) {
      warn("You have already reversed the DNA by using %d - %d truncation. Re-reversing",tstart,tend);
    }
    
    cdna_temp = cDNA_from_Sequence(reverse_complement_Sequence(cdna->baseseq));
    free_cDNA(cdna);
    cdna = cdna_temp;
  }

  if( target_abs == TRUE ) {
    cdna->baseseq->offset = 1;
    cdna->baseseq->end  = strlen(cdna->baseseq->seq);
  }
      

  if( qstart_str != NULL ) {
    if( is_integer_string(qstart_str,&qstart) == FALSE || qstart < 0) {
      warn("Could not make %s out as query start",qstart);
      ret = FALSE;
    }
  }

  if( qend_str != NULL ) {
    if( is_integer_string(qend_str,&qend) == FALSE || qend < 0) {
      warn("Could not make %s out as query end",qend);
      ret = FALSE;
    }
  }



  if( use_tsm == FALSE ) {

    if( startend != TSM_default && startend != TSM_global && startend != TSM_local ) {
      warn("Proteins can only have local/global startend policies set, not %s",startend_string);
      ret = FALSE;
    }

    if( (pro = read_fasta_file_Protein(protein_file)) == NULL ) {
      ret = FALSE;
      warn("Could not read Protein sequence in %s",protein_file);
    } else {
      if( qstart != -1 || qend != -1 ) {
	if( qstart == -1 )
	  qstart = 0;
	if( qend == -1 ) 
	  qend = pro->baseseq->len;

	pro_temp = truncate_Protein(pro,qstart,qend);
	if( pro_temp == NULL ){
	  ret = FALSE;
	} else {
	  free_Protein(pro);
	  pro = pro_temp;
	}
      }
    }
  } else {
    /** using a HMM **/

    tsm = HMMer2_read_ThreeStateModel(hmm_file);


    if( tsm == NULL ) {
      warn("Could not read hmm from %s\n",hmm_file);
      ret = FALSE;
    }  else {

      display_char_in_ThreeStateModel(tsm);
      if( hmm_name != NULL ) {
	ckfree(tsm->name);
	tsm->name = stringalloc(hmm_name);
      }
      
      if( tsm == NULL ) {
	warn("Could not read %s as a hmm",hmm_file);
      }
      
      /** have to set start/end **/
      /** already calculated    **/

      set_startend_policy_ThreeStateModel(tsm,startend,15,0.2);

    }

  } /* end of if HMM */


  if( strcmp(null_string,"syn") == 0 ) {
    use_syn = TRUE;
  } else if ( strcmp(null_string,"flat") == 0 ) {
    use_syn = FALSE;
  } else {
    warn("Cannot interpret [%s] as a null model string\n",null_string);
    ret = FALSE;
  }

  if( main_block_str != NULL ) {
    if( is_integer_string(main_block_str,&main_block) == FALSE ) {
      warn("Could not get maximum main_block number %s",main_block_str);
      ret = FALSE;
    }
  }
   
   
  if( alg_str != NULL ) {
    alg = alg_estwrap_from_string(alg_str);
  } else {
    if( use_tsm == TRUE ) {
      alg_str = "312L";
    } else {
      alg_str = "312";
    }
    alg = alg_estwrap_from_string(alg_str);
  }



  if( is_double_string(allN_string,&allN) == FALSE ) {
    warn("Could not convert %s to a double",allN_string);
    ret = FALSE;
  }

  if( is_double_string(subs_string,&subs_error) == FALSE ) {
    warn("Could not convert %s to a double",subs_error);
    ret = FALSE;
  }

  if( is_double_string(indel_string,&indel_error) == FALSE ) {
    warn("Could not convert %s to a double",indel_error);
    ret = FALSE;
  }

  if( (mat = read_Blast_file_CompMat(matrix_file)) == NULL) {
    ret = FALSE;
    warn("Could not read Comparison matrix file in %s",matrix_file);
  }

  if( (ct = read_CodonTable_file(codon_file)) == NULL) {
    ret = FALSE;
    warn("Could not read codon table file in %s",codon_file);
  }

  
  if( (ofp = openfile(output_file,"W")) ==  NULL) {
    warn("Could not open %s as an output file",output_file);
    ret = FALSE;
  }

  rmd = RandomModelDNA_std();

  rm = default_RandomModel();
  if( rm == NULL ) {
    warn("No random model. Ugh!");
    ret = FALSE;
  }

  cp = flat_cDNAParser(indel_error);
  cm = flat_CodonMapper(ct);
  sprinkle_errors_over_CodonMapper(cm,subs_error);

  return ret;

}

void build_defaults(void)
{
  codon_file = "codon.table";
  matrix_file = "BLOSUM62.bla";
  


}


void show_short_help(void)
{
  fprintf(stdout,"%s (%s)\n",program_name,VERSION_NUMBER);
  fprintf(stdout,"This program is freely distributed under a GPL. See -version for more info\n");
  fprintf(stdout,"Copyright (c) GRL limited: portions of the code are from separate copyright\n\n");
  fprintf(stdout,"estwise <protein-file> <dna-file> in fasta format\n");
  fprintf(stdout," Options. In any order, '-' as filename (for any input) means stdin\n");
  fprintf(stdout," Dna      [-u,-v,-trev,-tfor,-both,-tabs]\n");
  fprintf(stdout," Protein  [-s,-t,-g,-e,-m]\n");
  fprintf(stdout," HMM      [-hmmer,-hname]\n");
  fprintf(stdout," Model    [-codon,-subs,-indel,-null]\n Alg      [-kbyte,-alg]\n");
  fprintf(stdout," Output   [-pretty,-para,-sum,-alb,-pal,-block,-divide]\n");
  fprintf(stdout," Standard [-help,-version,-silent,-quiet,-errorlog]\n");
  fprintf(stdout,"\nFor more help go %s -help.\n",program_name);
  fprintf(stdout,"\nSee WWW help at http://www.sanger.ac.uk/Software/Wise2/\n");
  exit(63);   
}

void show_help(FILE * ofp)
{
  fprintf(ofp,"%s (%s)\n",program_name,VERSION_NUMBER);
  fprintf(ofp,"estwise <protein-file> <dna-file> in fasta format\n");
  /* program specific help */
  fprintf(ofp,"Dna sequence options\n");
  fprintf(ofp,"  -u               start position in dna\n");
  fprintf(ofp,"  -v               end position in dna\n");
  fprintf(ofp,"  -trev            reverse complement dna\n");
  fprintf(ofp,"  -tfor            use forward strands only\n");
  fprintf(ofp,"  -both [default]  do both strands\n");
  fprintf(ofp,"  -tabs            Positions reported as absolute to DNA\n");
  fprintf(ofp,"Protein comparison options\n");
  fprintf(ofp,"  -s               start position in protein\n");
  fprintf(ofp,"  -t               end   position in protein\n");
  fprintf(ofp,"  -[g]ap    [%3d]  gap penalty\n",gap);
  fprintf(ofp,"  -[e]xt    [%3d]  extension penalty\n",ext);
  fprintf(ofp,"  -[m]atrix [%s]   Comparison matrix\n",matrix_file);
  fprintf(ofp,"HMM options\n");
  fprintf(ofp,"  -hmmer           Protein file is HMMer file (version 2 compatible)\n");
  fprintf(ofp,"  -hname           Name of HMM rather than using the filename\n");
  fprintf(ofp,"Model options\n");
  fprintf(ofp,"  -init   [%s]    [default/global/local/wing] start end policy of HMM/protein match\n",startend_string);
  fprintf(ofp,"  -subs   [%2.2g] Substitution error rate\n",subs_error);
  fprintf(ofp,"  -indel  [%2.2g] Insertion/deletion error rate\n",indel_error);
  fprintf(ofp,"  -null   [syn/flat]   Random Model as synchronous or flat [default syn]\n");
  fprintf(ofp,"  -alln   [%s]   Probability of matching a NNN codon\n",allN_string);
  fprintf(ofp,"Algorithm options\n");
  fprintf(ofp,"  -alg    [333,333L,312,312L,333F]  Algorithm used\n");
  fprintf(ofp,"Output options [default -pretty -para]\n");
  fprintf(ofp,"  -pretty          show pretty ascii output\n");
  fprintf(ofp,"  -para            show parameters\n");
  fprintf(ofp,"  -sum             show summary information\n");
  fprintf(ofp,"  -pep             show protein translation, splicing frameshifts\n");
  fprintf(ofp,"  -alb             show logical AlnBlock alignment\n");
  fprintf(ofp,"  -pal             show raw matrix alignment\n");
  fprintf(ofp,"  -block  [%s]     Length of main block in pretty output\n",main_block_str);
  fprintf(ofp,"  -divide [%s]     divide string for multiple outputs\n",divide_str);
  show_help_DPRunImpl(ofp);
  show_standard_options(ofp);
  fprintf(ofp,"\nSee WWW help at http://www.sanger.ac.uk/Software/Wise2/\n");
  exit(63);   
}

void show_version(FILE * ofp)
{
  fprintf(ofp,"%s\nVersion: %s\nReleased: %s\nCompiled: %s\n",program_name,VERSION_NUMBER,RELEASE_DAY,COMPILE_DATE);
  fprintf(ofp,"\nThis program is freely distributed under a Gnu Public License\n");
  fprintf(ofp,"The source code is copyright (c) GRL 1998 and others\n");
  fprintf(ofp,"There is no warranty, implied or otherwise on the performance of this program\n");
  fprintf(ofp,"For more information read the GNULICENSE file in the distribution\n\n");
  fprintf(ofp,"Credits: Ewan Birney <birney@sanger.ac.uk> wrote the core code.\n");
  fprintf(ofp,"         Portions of this code is from HMMer2, written by Sean Eddy\n");
  exit(63);   
}


int main(int argc,char ** argv) 
{
  int i;
  char * temp;

  build_defaults();

  bootstrap_HMMer2();
  
  strip_out_standard_options(&argc,argv,show_help,show_version);

  dpri = new_DPRunImpl_from_argv(&argc,argv);

  if( (temp = strip_out_assigned_argument(&argc,argv,"gap")) != NULL )
    gap_str = temp;

  if( (temp = strip_out_assigned_argument(&argc,argv,"g")) != NULL )
    gap_str = temp;

  if( (temp = strip_out_assigned_argument(&argc,argv,"ext")) != NULL )
    ext_str = temp;

  if( (temp = strip_out_assigned_argument(&argc,argv,"init")) != NULL ) {
    startend_string = temp;
  }


  if( (temp = strip_out_assigned_argument(&argc,argv,"e")) != NULL )
    ext_str = temp;

  if( (temp = strip_out_assigned_argument(&argc,argv,"matrix")) != NULL )
    matrix_file = temp;

  if( (temp = strip_out_assigned_argument(&argc,argv,"m")) != NULL )
    matrix_file = temp;

  if( (temp = strip_out_assigned_argument(&argc,argv,"s")) != NULL )
    qstart_str = temp;

  if( (temp = strip_out_assigned_argument(&argc,argv,"t")) != NULL )
    qend_str = temp;

  if( (temp = strip_out_assigned_argument(&argc,argv,"u")) != NULL )
    tstart_str = temp;

  if( (temp = strip_out_assigned_argument(&argc,argv,"v")) != NULL )
    tend_str = temp;

  if( (temp = strip_out_assigned_argument(&argc,argv,"codon")) != NULL )
    codon_file = temp;

  if( (temp = strip_out_assigned_argument(&argc,argv,"null")) != NULL )
    null_string = temp;


  if( (temp = strip_out_assigned_argument(&argc,argv,"alg")) != NULL )
    alg_str = temp;


  if( (strip_out_boolean_argument(&argc,argv,"hmmer")) == TRUE )
    use_tsm = TRUE;


  if( (temp = strip_out_assigned_argument(&argc,argv,"hname")) != NULL )
    hmm_name = temp;

  if( (strip_out_boolean_argument(&argc,argv,"trev")) == TRUE ) {
    do_both = FALSE;
    reverse = TRUE;
  }

  if( (strip_out_boolean_argument(&argc,argv,"tfor")) == TRUE ) {
    if( reverse == TRUE ) {
      warn("You have specified both trev and tfor. Interpreting as both strands");
      do_both = TRUE;
    } else {
      do_both = FALSE;
      reverse = FALSE;
    }
  }


  if( (temp = strip_out_assigned_argument(&argc,argv,"subs")) != NULL )
    subs_string = temp;

  if( (temp = strip_out_assigned_argument(&argc,argv,"indel")) != NULL )
    indel_string = temp;

  if( (strip_out_boolean_argument(&argc,argv,"both")) == TRUE )
    do_both = TRUE;

  if( (strip_out_boolean_argument(&argc,argv,"tabs")) == TRUE )
    target_abs = TRUE;
   
  if( (strip_out_boolean_argument(&argc,argv,"pretty")) != FALSE )
    show_pretty = TRUE;

  if( (strip_out_boolean_argument(&argc,argv,"para")) != FALSE )
    show_para = TRUE;

  if( (strip_out_boolean_argument(&argc,argv,"sum")) != FALSE )
    show_match_sum = TRUE;

  if( (strip_out_boolean_argument(&argc,argv,"pep")) != FALSE )
    show_pep = TRUE;

  if( (strip_out_boolean_argument(&argc,argv,"alb")) != FALSE )
    show_AlnBlock = TRUE;

  if( (temp = strip_out_assigned_argument(&argc,argv,"alln")) != NULL )
    allN_string = temp;


  if( (strip_out_boolean_argument(&argc,argv,"pal")) != FALSE )
    show_PackAln = TRUE;


  if( (temp = strip_out_assigned_argument(&argc,argv,"divide")) != NULL )
    divide_str = temp;

  if( (temp = strip_out_assigned_argument(&argc,argv,"block")) != NULL )
    main_block_str = temp;

  strip_out_remaining_options_with_warning(&argc,argv);
  
  if( argc !=  3 ) {
    warn("Wrong number of arguments (expect 2)! Arg line looked like (after option processing)");
    for(i=1;i<argc;i++) {
      fprintf(stderr,"   %s\n",argv[i]);
    }
    show_short_help();
  }

  if(  show_match_sum == FALSE && show_pretty == FALSE && show_AlnBlock == FALSE && show_PackAln == FALSE && show_pep == FALSE ) {
    show_pretty = TRUE;
    show_para = TRUE;
  }
 
  dna_seq_file = argv[2];
  if( use_tsm == FALSE) 
    protein_file = argv[1];
  else 
    hmm_file  = argv[1];

  if( build_objects() == FALSE) 
    fatal("Could not build objects!");

  if( show_para == TRUE)
    show_parameters();

  if( build_alignment() == FALSE)
    fatal("Could not build alignment!");

  if( show_output() == FALSE)
    fatal("Could not show alignment. Sorry!");

  if( do_both == TRUE) {
    reverse_target();

    if( build_alignment() == FALSE)
      fatal("Could not build alignment!");

    if( show_output() == FALSE)
      fatal("Could not show alignment. Sorry!");
  }

  return 0;
}




