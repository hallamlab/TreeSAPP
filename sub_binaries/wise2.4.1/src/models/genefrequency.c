#ifdef _cplusplus
extern "C" {
#endif
#include "genefrequency.h"

/* Function:  RandomCodon_from_cds_triplet(gf)
 *
 * Descrip:    makes a randomcodon probability emission
 *             from the counts in genefrequency
 *
 *
 * Arg:        gf [UNKN ] Undocumented argument [GeneFrequency21 *]
 *
 * Return [UNKN ]  Undocumented return value [RandomCodon *]
 *
 */
# line 83 "genefrequency.dy"
RandomCodon * RandomCodon_from_cds_triplet(GeneFrequency21 * gf)
{
  int i;
  int a,b,c;
  double total;
  double lit;

  RandomCodon * out;

  out = RandomCodon_alloc();

  for(i=0,total=0;i<64;i++)
    total += gf->cds_triplet[i];
  
  for(i=0;i<125;i++) {
    if( has_random_bases(i) ) {
      lit = 0.0;
      for(a=0;a<4;a++)
	for(b=0;b<4;b++)
	  for(c=0;c<4;c++)
	    lit += gf->cds_triplet[base4_codon_from_codon(permute_possible_random_bases(i,a,b,c))];
      out->codon[i] = lit/(64*total);
    } else {
      out->codon[i] = gf->cds_triplet[base4_codon_from_codon(i)]/(total);
    }
  }

  return out;
}


/* Function:  RandomModelDNA_from_central_GeneFrequency21(gf)
 *
 * Descrip:    Makes a random model from the central gene 
 *             model of an intron. Ideal for tieing intron
 *             state distribution to the randommodel
 *
 *
 * Arg:        gf [UNKN ] Undocumented argument [GeneFrequency21 *]
 *
 * Return [UNKN ]  Undocumented return value [RandomModelDNA *]
 *
 */
# line 119 "genefrequency.dy"
RandomModelDNA * RandomModelDNA_from_central_GeneFrequency21(GeneFrequency21 * gf)
{
  RandomModelDNA * out;
  double total;

  total = sum_Probability_array(gf->central,4);


  out = RandomModelDNA_alloc();

  out->base[BASE_A] = gf->central[BASE_A]/total;
  out->base[BASE_T] = gf->central[BASE_T]/total;
  out->base[BASE_G] = gf->central[BASE_G]/total;
  out->base[BASE_C] = gf->central[BASE_C]/total;
  out->base[BASE_N] = 1.0;

  return out;
}

/* Function:  ComplexConsensi_5SS_from_GeneFrequency(gf)
 *
 * Descrip:    makes 5'SS ComplexConsensi from GeneFrequency21 structure using
 *
 *               CCC|XXXXXXX score = no(5'SS with CCC|XXXXXXX) / no(CCC in cds).
 *
 *
 * Arg:        gf [UNKN ] Undocumented argument [GeneFrequency21 *]
 *
 * Return [UNKN ]  Undocumented return value [ComplexConsensi *]
 *
 */
# line 143 "genefrequency.dy"
ComplexConsensi * ComplexConsensi_5SS_from_GeneFrequency(GeneFrequency21 * gf)
{
  register int i;
  double nocds;
  ComplexConsensi * out;

  out = ComplexConsensi_alloc_len(gf->ss5->len);



  for(i=0;i<gf->ss5->len;i++) {
    nocds = nocds_from_ambiguous_codon(gf->ss5->gsc[i]->string,gf->cds_triplet);
    if( nocds < 20 ) {
      warn("In making 5'SS consensi, got %g cds for codon %s ... not happy about this",nocds,gf->ss5->gsc[i]->string);
    }
    add_ComplexConsensi(out,ComplexConsensusWord_from_string_and_prob(gf->ss5->gsc[i]->string,gf->ss5->gsc[i]->number / nocds));
  }

  return out;
  
}

/* Function:  ComplexConsensi_3SS_from_GeneFrequency(gf)
 *
 * Descrip:    makes 3'SS ComplexConsensi from GeneFrequency21 structure using
 *
 *               ZZZ|CCC score = no(3'SS with ZZZ|CCC) / no(CCC in cds).
 *
 *
 * Arg:        gf [UNKN ] Undocumented argument [GeneFrequency21 *]
 *
 * Return [UNKN ]  Undocumented return value [ComplexConsensi *]
 *
 */
# line 170 "genefrequency.dy"
ComplexConsensi * ComplexConsensi_3SS_from_GeneFrequency(GeneFrequency21 * gf)
{
  register int i;
  double nocds;
  ComplexConsensi * out;

  out = ComplexConsensi_alloc_len(gf->ss3->len);

  for(i=0;i<gf->ss3->len;i++) {
    nocds = nocds_from_ambiguous_codon(gf->ss3->gsc[i]->string+3,gf->cds_triplet);
    if( nocds < 20 ) {
      warn("In making 3'SS consensi, got %g cds for codon %s ... not happy about this!",nocds,gf->ss3->gsc[i]->string);
    }
    add_ComplexConsensi(out,ComplexConsensusWord_from_string_and_prob(gf->ss3->gsc[i]->string,gf->ss3->gsc[i]->number / nocds));
  }

  return out;
  
}

/* Function:  nocds_from_ambiguous_codon(codon,codon_freq_array)
 *
 * Descrip:    helper function for above guys
 *
 *
 * Arg:                   codon [UNKN ] Undocumented argument [char *]
 * Arg:        codon_freq_array [UNKN ] Undocumented argument [double *]
 *
 * Return [UNKN ]  Undocumented return value [double]
 *
 */
# line 194 "genefrequency.dy"
double nocds_from_ambiguous_codon(char * codon,double * codon_freq_array)
{
  int factor = 1;
  int one;
  int two;
  int three;
  int i,j,k;
  double ret = 0.0;

  one = base_from_char(*codon == '-' ? 'N' : *codon);
  two = base_from_char(*(codon+1) == '-' ? 'N' : *(codon+1));
  three = base_from_char(*(codon+2) == '-' ? 'N' : *(codon+2));


  if(one == BASE_N)
    factor *= 4;
  if(two == BASE_N)
    factor *= 4;
  if(three == BASE_N)
    factor *= 4;
  

  for(i=0;i<4;i++)
    for(j=0;j<4;j++) 
      for(k=0;k<4;k++) 
	if( (one == i || one == BASE_N) && (two == j || two == BASE_N) && (three == k || three == BASE_N)) { 
	  ret += codon_freq_array[i*16+j*4+k];
	}
  

  ret = ret / factor;

  if( ret < 0.0000000000000001 ) {
    warn("For codon  %c%c%c we have a frequency of %g",*codon,*(codon+1),*(codon+2),ret);
    ret = 0.0000000000000001;
  }

  return ret;
}

/* Function:  ComplexConsensusWord_from_string_and_prob(string,p)
 *
 * Descrip:    convienent constructor
 *
 *
 * Arg:        string [UNKN ] Undocumented argument [char *]
 * Arg:             p [UNKN ] Undocumented argument [Probability]
 *
 * Return [UNKN ]  Undocumented return value [ComplexConsensusWord *]
 *
 */
# line 238 "genefrequency.dy"
ComplexConsensusWord * ComplexConsensusWord_from_string_and_prob(char * string,Probability p)
{
  ComplexConsensusWord * out;

  out = ComplexConsensusWord_alloc();

  out->pattern = stringalloc(string);
  out->p = p;
  out->score = Probability2Score(p);

  return out;
}


/* Function:  CodonFrequency_from_GeneFrequency21(gf,ct)
 *
 * Descrip:    Builds a codon frequency table from raw counts
 *             in the counts file
 *
 *
 * Arg:        gf [UNKN ] Undocumented argument [GeneFrequency21 *]
 * Arg:        ct [UNKN ] Undocumented argument [CodonTable *]
 *
 * Return [UNKN ]  Undocumented return value [CodonFrequency *]
 *
 */
# line 256 "genefrequency.dy"
CodonFrequency * CodonFrequency_from_GeneFrequency21(GeneFrequency21 * gf,CodonTable * ct)
{
  return CodonFrequence_from_raw_counts(gf->codon,ct);
}


/* Function:  show_flat_GeneFrequency21(gf21,ofp)
 *
 * Descrip:    For debugging
 *
 *
 * Arg:        gf21 [UNKN ] Undocumented argument [GeneFrequency21 *]
 * Arg:         ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 266 "genefrequency.dy"
void show_flat_GeneFrequency21(GeneFrequency21 * gf21,FILE * ofp)
{
  if( gf21->ss5 != NULL ) {
    fprintf(ofp,"5'SS\n");
    show_GeneConsensus(gf21->ss5,ofp);
  }
  if( gf21->ss3 != NULL ) {
    fprintf(ofp,"3'SS\n");
    show_GeneConsensus(gf21->ss3,ofp);
  }

  fprintf(ofp,"Codon frequency\n");
  show_codon_emission(gf21->codon,ofp);

  fprintf(ofp,"Central emission\n");
  show_base_emission(gf21->central,ofp);

  fprintf(ofp,"Pyrimidine emission\n");
  show_base_emission(gf21->py,ofp);

  fprintf(ofp,"Spacer emission\n");
  show_base_emission(gf21->spacer,ofp);

  fprintf(ofp,"Transitions\n");
  show_Probability_array(gf21->transition,GENEFREQUENCY21_TRANSITION_LEN,ofp);

  fprintf(ofp,"\n\n");
}

/* Function:  show_codon_emission(codon,ofp)
 *
 * Descrip:    For debugging
 *
 *
 * Arg:        codon [UNKN ] Undocumented argument [double *]
 * Arg:          ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 299 "genefrequency.dy"
void show_codon_emission(double * codon,FILE * ofp)
{
  register int i;

  fprintf(ofp,"begin consensus\n");

  for(i=0;i<64;i++) 
    show_single_codon_emission(codon[i],i,ofp);

  fprintf(ofp,"end consensus\n");
}

/* Function:  show_single_codon_emission(no,base4codon,ofp)
 *
 * Descrip:    For debugging
 *
 *
 * Arg:                no [UNKN ] Undocumented argument [double]
 * Arg:        base4codon [UNKN ] Undocumented argument [int]
 * Arg:               ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 315 "genefrequency.dy"
void show_single_codon_emission(double no,int base4codon,FILE * ofp)
{
  codon c;
  base one;
  base two;
  base three;

  c = codon_from_base4_codon(base4codon);

  all_bases_from_codon(c,&one,&two,&three);

  fprintf(ofp,"%c%c%c %.2f\n",char_from_base(one),char_from_base(two),char_from_base(three),no);
}

/* Function:  show_base_emission(base,ofp)
 *
 * Descrip:    For debugging
 *
 *
 * Arg:        base [UNKN ] Undocumented argument [double *]
 * Arg:         ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 333 "genefrequency.dy"
void show_base_emission(double * base,FILE * ofp)
{
  register int i;

  fprintf(ofp,"begin consensus\n");

  for(i=0;i<4;i++)
    fprintf(ofp,"%c %.2f\n",char_from_base(i),base[i]);

  fprintf(ofp,"end consensus\n");
}


/* Function:  show_GeneConsensus(gc,ofp)
 *
 * Descrip:    For debugging
 *
 *
 * Arg:         gc [UNKN ] Undocumented argument [GeneConsensus *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 350 "genefrequency.dy"
void show_GeneConsensus(GeneConsensus * gc,FILE * ofp)
{
  register int i;

  fprintf(ofp,"begin consensus\n");

  for(i=0;i<gc->len;i++)
    show_GeneSingleCons(gc->gsc[i],ofp);

  fprintf(ofp,"end consensus\n");
}


/* Function:  show_GeneSingleCons(gsc,ofp)
 *
 * Descrip:    For debugging
 *
 *
 * Arg:        gsc [UNKN ] Undocumented argument [GeneSingleCons *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 367 "genefrequency.dy"
void show_GeneSingleCons(GeneSingleCons * gsc,FILE * ofp)
{
  fprintf(ofp,"%s %f\n",gsc->string,gsc->number);
}




 /*** reading in ***/

/* Function:  read_GeneFrequency21_file(filename)
 *
 * Descrip:    Opens the file with /openfile
 *
 *             Reads in a GeneFrequency (Mor-Ewan style)
 *
 *
 *
 * Arg:        filename [UNKN ] will open from WISECONFIGDIR etc via openfile [char *]
 *
 * Return [UNKN ]  a newly allocated structure [GeneFrequency21 *]
 *
 */
# line 386 "genefrequency.dy"
GeneFrequency21 * read_GeneFrequency21_file(char * filename)
{
  GeneFrequency21 * out;
  FILE * ifp;

  ifp = openfile(filename,"r");

  if( ifp == NULL ) {
    warn("Could not open %s as a genefrequency file",filename);
    return NULL;
  }

  out = read_GeneFrequency21(ifp);

  fclose(ifp);

  return out;
}

/* Function:  read_GeneFrequency21(ifp)
 *
 * Descrip:    Reads in a GeneFrequency (Mor-Ewan style)
 *             file from ifp
 *
 *
 * Arg:        ifp [UNKN ] file pointer [FILE *]
 *
 * Return [UNKN ]  a newly allocated structure [GeneFrequency21 *]
 *
 */
# line 412 "genefrequency.dy"
GeneFrequency21 * read_GeneFrequency21(FILE * ifp)
{
  GeneFrequency21 * out;
  GeneConsensus   * temp;
  char buffer[MAXLINE];
  int phase;
  int center;
  int type;
  boolean err = FALSE;

  out = GeneFrequency21_alloc();

  while( fgets(buffer,MAXLINE,ifp) != NULL ) {
    if( buffer[0] == '#' )
      continue;


    if( strwhitestartcmp(buffer,"type",spacestr) == 0 ) {

      phase = 3; /** if no phase, assumme it is for all phases **/

      type = check_type_GeneFrequency(buffer,ifp,&center,&phase);
      
      switch(type) {
      case GeneConsensusType_5SS :
	if( phase == 3) {
	  temp = read_line_GeneConsensus(buffer,ifp);
	  temp->center = center;
	  out->ss5 = temp;
	}
	else {
	  if( skip_consensus(ifp) == FALSE ) {
	    warn("Unable to skip phase'd 5'SS information ... problem!");
	    break;
	  }
	}
	break;
      case GeneConsensusType_3SS :
	if( phase == 3) {
	  temp = read_line_GeneConsensus(buffer,ifp);
	  temp->center = center;
	  out->ss3 = temp;
	}
	else {
	  if( skip_consensus(ifp) == FALSE ) {
	    warn("Unable to skip phase'd 5'SS information ... problem!");
	    err = TRUE;
	  }
	}
	break;
      case GeneConsensusType_CDS :
	if( phase == 0) {
	  if( read_codon_GeneConsensus(out->codon,buffer,ifp) == FALSE ) {
	    warn("Unable to read codon information in GeneFrequency21... problem!");
	    break;
	  }
	}
	else if( phase == 3 ) {
	  /*** we need this! ***/
	  if( read_codon_GeneConsensus(out->cds_triplet,buffer,ifp) == FALSE ) {
	    warn("Unable to read codon information in GeneFrequency21... problem!");
	    break;
	  }
	}
	else { /** in a different phase **/
	  if( skip_consensus(ifp) == FALSE ) {
	    warn("Unable to skip phase'd CDS information ... problem!");
	    err = TRUE;
	  }
	}
	break;
      case GeneConsensusType_Intron_emission :
	if( phase == 3 ) {
	  if( read_base_GeneConsensus(out->central,buffer,ifp) == FALSE ) {
	    warn("Unable to read Intron emissions in genefrequency21 ... problem!");
	    err = TRUE;
	  }
	}
	else {
	  if( skip_consensus(ifp) == FALSE ) {
	    warn("Unable to skip phase'd CDS information ... problem!");
	    err = TRUE;
	  }
	}
	break;
      case GeneConsensusType_Pyrimidine_emission :
	if( phase == 3 ) {
	  if( read_base_GeneConsensus(out->py,buffer,ifp) == FALSE ) {
	    warn("Unable to read pyrimidine emissions in genefrequency21 ... problem!");
	    err = TRUE;
	  }
	}
	else {
	  if( skip_consensus(ifp) == FALSE ) {
	    warn("Unable to skip phase'd pyrimidine information ... problem!");
	    err = TRUE;
	  }
	}
	break;
      case GeneConsensusType_Spacer_emission :
	if( phase == 3 ) {
	  if( read_base_GeneConsensus(out->spacer,buffer,ifp) == FALSE ) {
	    warn("Unable to read spacer emissions in genefrequency21 ... problem!");
	    err = TRUE;
	  }
	}
	else {
	  if( skip_consensus(ifp) == FALSE ) {
	    warn("Unable to skip phase'd spacer information ... problem!");
	    err = TRUE;
	  }
	}
	break;
      case GeneConsensusType_Central_stay :
	out->transition[GF21_CENTRAL_STAY] = double_from_line(buffer);
	break;
      case GeneConsensusType_Pyrimidine_stay :
	out->transition[GF21_PY_STAY] = double_from_line(buffer);
	break;
      case GeneConsensusType_Spacer_stay :
	out->transition[GF21_SPACER_STAY] = double_from_line(buffer);
	break;
      case GeneConsensusType_No_spacer :
	out->transition[GF21_NO_SPACER] = double_from_line(buffer);
	break;
      case GeneConsensusType_Intron_Corr_Term :
	switch(phase) {
	case 0 :
/*	  out->transition[GF21_INTRON_CORR_TERM_0] = double_from_line(buffer); */
	  break;
	case 1 :
/*	  out->transition[GF21_INTRON_CORR_TERM_1] = double_from_line(buffer); */
	  break;
	case 2 :
/*	  out->transition[GF21_INTRON_CORR_TERM_2] = double_from_line(buffer); */
	  break;
	case 3 :
	  out->transition[GF21_INTRON_CORR_TERM] = double_from_line(buffer); 
	  break;
	default :
	  warn("Well... I have got some bad news for you. We found a phase of %d in Intron correction term. ",phase);
	  break;
	}
	break;
      default :
	warn("Got an unidenitifable type in GeneFrequency21 parse. Skippping");
	if( skip_consensus(ifp) == FALSE ) {
	  warn("Unable to skip phase'd 5'SS information ... problem!");
	  err = TRUE;
	}
	
      }

      if( err == TRUE ) {
	warn("You have had an unrecoverable error in GeneFrequency21 parsing");
	break;
      }

    }
    else {
      striptoprint(buffer);
      warn("Could not understand line [%s] in GeneFrequency21 parse",buffer);
    }
  }

  return out;
}
	
	
/* Function:  double_from_line(buffer)
 *
 * Descrip:    helper string function
 *
 *
 * Arg:        buffer [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [double]
 *
 */
# line 584 "genefrequency.dy"
double double_from_line(char * buffer)
{
  char * runner;
  char * end;
  double ret;

  runner = strtok(buffer,spacestr);

  if( runner == NULL ) {
    warn("Unable to read a number in double_from_line");
    return -1.0;
  }

  ret = strtod(runner,&end);

  if( end == runner || isalnum((int)*end) ) {
    warn("Bad conversion of string [%s] to double [%f] occured",runner,ret);
  }

  return ret;
}
  
/* Function:  skip_consensus(ifp)
 *
 * Descrip:    helper function for
 *             teh file parsing
 *
 *
 * Arg:        ifp [UNKN ] Undocumented argument [FILE *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 611 "genefrequency.dy"
boolean skip_consensus(FILE * ifp)
{
  char buffer[MAXLINE];


  while(fgets(buffer,MAXLINE,ifp) != NULL ) 
    if( strwhitestartcmp(buffer,"end",spacestr) == 0)
      break;

  if( feof(ifp) || ferror(ifp) )
    return FALSE;
  return TRUE;
}


/* Function:  check_type_GeneFrequency(*line,ifp,center,phase)
 *
 * Descrip:    Pretty sneaky function 
 *
 *             give 
 *
 *               line starting with "type xxx"
 *               ifp  file pointer
 *               centre a &int for returning the centre value of the consensus if any.
 *               phase  a &int for returning the phase value of the consensus if any.
 *
 *             you get *back* the line with the line "begin consensus" or "number" in it.
 *
 *               
 *             It returns a GeneConsensusType
 *
 *             with GeneConsensusType_Error on error
 *
 *
 *
 * Arg:         *line [UNKN ] Undocumented argument [char]
 * Arg:           ifp [UNKN ] Undocumented argument [FILE *]
 * Arg:        center [UNKN ] Undocumented argument [int *]
 * Arg:         phase [UNKN ] Undocumented argument [int *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
# line 644 "genefrequency.dy"
int check_type_GeneFrequency(char *line,FILE * ifp,int * center,int * phase)  
{
  int ret = GeneConsensusType_Error;
  char * runner;


  if( strwhitestartcmp(line,"type",spacestr) != 0 ) {
    warn("Attempting to check phase of consensus with no type line...");
    return GeneConsensusType_Error;
  }

  runner = strtok(line,spacestr);
  runner = strtok(NULL,spacestr);

  if( runner == NULL ) {
    warn("GeneFrequency type with no type. Can't read type, must set to error, but problem in later parsing");
    ret = GeneConsensusType_Error;
  }

  else {
    ret = string_to_GeneConsensusType(runner);
  }

  while( fgets(line,MAXLINE,ifp) != NULL ) {
    if( line[0] == '#' )
      continue;

    else if( strwhitestartcmp(line,"phase",spacestr) == 0 ) {
      runner = strtok(line,spacestr);
      runner = strtok(NULL,spacestr);

      if( runner == NULL ) {
	warn("Got phase line with no phase. Sad....");
	continue;
      }
      

      if( phase != NULL ) {
	if( strcmp(runner,"all") ==0 || strcmp(runner,"All") == 0)
	  *phase = 3;
	else *phase = atoi(runner);
      }
    }

    else if( strwhitestartcmp(line,"center",spacestr) == 0 || strwhitestartcmp(line,"centre",spacestr) == 0) {
      runner = strtok(line,spacestr);
      runner = strtok(NULL,spacestr);

      if( runner == NULL ) {
	warn("Got center line with no phase. Sad....");
	continue;
      }

      if( center != NULL ) {
	*center = atoi(runner);
      }
    }
    else {
      break;
    }
  }



  return ret;
}

/* Function:  string_to_GeneConsensusType(string)
 *
 * Descrip:    flips string 5SS to enum type
 *
 *
 * Arg:        string [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
# line 715 "genefrequency.dy"
int string_to_GeneConsensusType(char * string)
{


  if( strcmp(string,"5SS") == 0 )
    return GeneConsensusType_5SS;
  else if( strcmp(string,"3SS") == 0 )
    return GeneConsensusType_3SS;
  else if( strcmp(string,"CDS") == 0 )
    return GeneConsensusType_CDS;
  else if( strcmp(string,"Intron_Corr_Term") == 0 )
    return  GeneConsensusType_Intron_Corr_Term;
  else if( strcmp(string,"Intron_emission") == 0 )
    return  GeneConsensusType_Intron_emission;
  else if( strcmp(string,"Pyrimidine_emission") == 0 )
    return  GeneConsensusType_Pyrimidine_emission;
  else if( strcmp(string,"Spacer_emission") == 0 )
    return GeneConsensusType_Spacer_emission;
  else if( strcmp(string,"Central_Intron_Stay_Prob") == 0 )
    return GeneConsensusType_Central_stay;
  else if( strcmp(string,"Pyrimidine_Stay_Prob") == 0 )
    return GeneConsensusType_Pyrimidine_stay;
  else if( strcmp(string,"Spacer_Stay_Prob") == 0 )
    return GeneConsensusType_Spacer_stay;
  else if( strcmp(string,"No_Spacer_Prob") == 0 )
    return GeneConsensusType_No_spacer;
  else {
    warn("Could convert string [%s] into a gene frequency type",string);
    return GeneConsensusType_Error;
  }
}
  
/* Function:  read_base_GeneConsensus(base_array,line,ifp)
 *
 * Descrip:    assummes base_array is 4 positions long
 *               
 *             line should have begin consensus on it and be of MAXLINE length as it will be used as the buffer.
 *               
 *             This does **not** check that you have filled up all 4 positions.
 *
 *
 * Arg:        base_array [UNKN ] Undocumented argument [double *]
 * Arg:              line [UNKN ] Undocumented argument [char*]
 * Arg:               ifp [UNKN ] Undocumented argument [FILE *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 754 "genefrequency.dy"
boolean read_base_GeneConsensus(double * base_array,char* line,FILE * ifp)
{
  boolean ret = TRUE;
  int b;
  char * base;
  char * number;


  if( strwhitestartcmp(line,"begin",spacestr) != 0 || strstr(line,"consensus") == NULL ) {
    warn("In reading base GeneConsensus line, got no 'begin consensus' tag [%s]",line);
    return FALSE;
  }


  while( fgets(line,MAXLINE,ifp) != NULL ) {
    if( line[0] == '#' )
      continue;

    if( strwhitestartcmp(line,"end",spacestr) == 0 )
      break;

    base = strtok(line,spacestr);
    number = strtok(NULL,spacestr);

    if( base == NULL ) {
      warn("Found an uncommented line in base consensus with no leading base word");
      continue;
    }

    if( number == NULL ) {
      warn("For base %s, no number found",base);
      ret = FALSE;
      continue;
    }

    if( strlen(base) > 1 || (b=base_from_char(*base)) == BASE_N ) {
      warn("Could not interpret %s as an actual DNA base in read_base_GeneConsensus");
      ret = FALSE;
      continue;
    }

    base_array[b]= atof(number);

  }

  return ret;
}



/* Function:  read_codon_GeneConsensus(codon_array,line,ifp)
 *
 * Descrip:    assummes codon_array is 64 positions long
 *               
 *             line should have begin consensus on it and be of MAXLINE length as it will be used as the buffer.
 *
 *             This does **not** check that you have filled up all 64 positions.
 *
 *
 * Arg:        codon_array [UNKN ] Undocumented argument [double *]
 * Arg:               line [UNKN ] Undocumented argument [char*]
 * Arg:                ifp [UNKN ] Undocumented argument [FILE *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 811 "genefrequency.dy"
boolean read_codon_GeneConsensus(double * codon_array,char* line,FILE * ifp)
{
  boolean ret = TRUE;
  char * codon;
  char * number;


  if( strwhitestartcmp(line,"begin",spacestr) != 0 || strstr(line,"consensus") == NULL ) {
    warn("In reading codon GeneConsensus line, got no 'begin consensus' tag [%s]",line);
    return FALSE;
  }


  while( fgets(line,MAXLINE,ifp) != NULL ) {
    if( line[0] == '#' )
      continue;

    if( strwhitestartcmp(line,"end",spacestr) == 0 )
      break;

    codon = strtok(line,spacestr);
    number = strtok(NULL,spacestr);

    if( codon == NULL ) {
      warn("Found an uncommented line in codon consensus with no leading codon word");
      continue;
    }

    if( number == NULL ) {
      warn("For codon %s, no number found",codon);
      ret = FALSE;
      continue;
    }

    if( strchr(codon,'N') != NULL ) 
      continue;

    if( is_non_ambiguous_codon_seq(codon) == FALSE ) {
      warn("Codon %s is not really a codon... problem!");
      ret = FALSE;
      continue;
    }



    codon_array[base4_codon_from_seq(codon)]= atof(number);

  }

  return ret;
}
    
/* Function:  read_line_GeneConsensus(line,ifp)
 *
 * Descrip:    Reads a single GeneConsensus from a file
 *
 *
 * Arg:        line [UNKN ] Undocumented argument [char *]
 * Arg:         ifp [UNKN ] Undocumented argument [FILE *]
 *
 * Return [UNKN ]  Undocumented return value [GeneConsensus  *]
 *
 */
# line 866 "genefrequency.dy"
GeneConsensus  * read_line_GeneConsensus(char * line,FILE * ifp)
{
  GeneConsensus * out;
  GeneSingleCons * temp;
  char buffer[MAXLINE];
  char * runner;



  if( strwhitestartcmp(line,"begin",spacestr) != 0 ) {
    warn("Attempting to read a GeneConsensus structure with a line not starting with 'begin' [%s]",line);
    return NULL;
  }

  runner = strtok(line,spacestr);
  runner = strtok(NULL,spacestr);

  if( runner == NULL || strcmp(runner,"consensus") != 0 ) {
    warn("Attempting to read a GeneConsensus structure without a 'begin consensus' tag [%s]",line);
    return NULL;
  }

  out = GeneConsensus_alloc_std();

  while( fgets(buffer,MAXLINE,ifp) != NULL ) {
    if( buffer[0] == '#' )
      continue;
    if( strwhitestartcmp(buffer,"end",spacestr) == 0 )
      break;
    
   
    temp = read_line_GeneSingleCons(buffer);

    if( temp == NULL ) {
      warn("Unable to process GeneSingleCons line... dropping out...");
      break;
    }

    add_GeneConsensus(out,temp);
  }


  return out;
}
    


# line 913 "genefrequency.dy"
GeneSingleCons * read_line_GeneSingleCons(char * line)
{
  GeneSingleCons * out;
  char * runner;
  char * run2;


  runner = strtok(line,spacestr);
  run2   = strtok(NULL,spacestr);

  if( runner == NULL || run2 == NULL ) {
    warn("In read_line_GeneSingleCons was not give two different words in line [%s]",line);
    return NULL;
  }

  out = GeneSingleCons_alloc();

  out->string = stringalloc(runner);

  out->number = strtod(run2,&runner);

  if( runner == run2 || *runner != '\0' ) {
    warn("In read_line_GeneSingleCons, for string [%s], unable to convert the number [%s]",out->string,run2);
  }

  return out;
}


# line 942 "genefrequency.dy"
GeneFrequency21  * untouched_GeneFrequency21(void)
{
  register int i;
  GeneFrequency21 * out;

  out = GeneFrequency21_alloc();

  for(i=0;i<64;i++)
    out->codon[i] = (-1);

  for(i=0;i<4;i++)
    out->central[i] = out->py[i] = out->spacer[i] = (-1);
  
  for(i=0;i<GENEFREQUENCY21_TRANSITION_LEN;i++) 
    out->transition[i] = (-1.0);

  return out;
}



# line 1013 "genefrequency.c"
/* Function:  hard_link_GeneSingleCons(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GeneSingleCons *]
 *
 * Return [UNKN ]  Undocumented return value [GeneSingleCons *]
 *
 */
GeneSingleCons * hard_link_GeneSingleCons(GeneSingleCons * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a GeneSingleCons object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  GeneSingleCons_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GeneSingleCons *]
 *
 */
GeneSingleCons * GeneSingleCons_alloc(void) 
{
    GeneSingleCons * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(GeneSingleCons *) ckalloc (sizeof(GeneSingleCons))) == NULL)    {  
      warn("GeneSingleCons_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->string = NULL;  
    out->number = 0; 


    return out;  
}    


/* Function:  free_GeneSingleCons(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GeneSingleCons *]
 *
 * Return [UNKN ]  Undocumented return value [GeneSingleCons *]
 *
 */
GeneSingleCons * free_GeneSingleCons(GeneSingleCons * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a GeneSingleCons obj. Should be trappable");    
      return NULL;   
      }  


#ifdef PTHREAD   
    assert(pthread_mutex_lock(&(obj->dynamite_mutex)) == 0); 
#endif   
    if( obj->dynamite_hard_link > 1)     {  
      return_early = 1;  
      obj->dynamite_hard_link--; 
      }  
#ifdef PTHREAD   
    assert(pthread_mutex_unlock(&(obj->dynamite_mutex)) == 0);   
#endif   
    if( return_early == 1)   
      return NULL;   
    if( obj->string != NULL) 
      ckfree(obj->string);   


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_GeneConsensus(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_GeneConsensus
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [GeneSingleCons **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_GeneConsensus(GeneSingleCons ** list,int i,int j)  
{
    GeneSingleCons * temp;   
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_GeneConsensus(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_GeneConsensus which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [GeneSingleCons **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_GeneConsensus(GeneSingleCons ** list,int left,int right,int (*comp)(GeneSingleCons * ,GeneSingleCons * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_GeneConsensus(list,left,(left+right)/2);    
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_GeneConsensus (list,++last,i);  
      }  
    swap_GeneConsensus (list,left,last); 
    qsort_GeneConsensus(list,left,last-1,comp);  
    qsort_GeneConsensus(list,last+1,right,comp); 
}    


/* Function:  sort_GeneConsensus(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_GeneConsensus
 *
 *
 * Arg:         obj [UNKN ] Object containing list [GeneConsensus *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_GeneConsensus(GeneConsensus * obj,int (*comp)(GeneSingleCons *, GeneSingleCons *)) 
{
    qsort_GeneConsensus(obj->gsc,0,obj->len-1,comp); 
    return;  
}    


/* Function:  expand_GeneConsensus(obj,len)
 *
 * Descrip:    Really an internal function for add_GeneConsensus
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [GeneConsensus *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_GeneConsensus(GeneConsensus * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_GeneConsensus called with no need");  
      return TRUE;   
      }  


    if( (obj->gsc = (GeneSingleCons ** ) ckrealloc (obj->gsc,sizeof(GeneSingleCons *)*len)) == NULL)     {  
      warn("ckrealloc failed for expand_GeneConsensus, returning FALSE");    
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_GeneConsensus(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [GeneConsensus *]
 * Arg:        add [OWNER] Object to add to the list [GeneSingleCons *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_GeneConsensus(GeneConsensus * obj,GeneSingleCons * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_GeneConsensus(obj,obj->len + GeneConsensusLISTLENGTH) == FALSE) 
        return FALSE;    
      }  


    obj->gsc[obj->len++]=add;    
    return TRUE; 
}    


/* Function:  flush_GeneConsensus(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [GeneConsensus *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_GeneConsensus(GeneConsensus * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->gsc[i] != NULL)   {  
        free_GeneSingleCons(obj->gsc[i]);    
        obj->gsc[i] = NULL;  
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  GeneConsensus_alloc_std(void)
 *
 * Descrip:    Equivalent to GeneConsensus_alloc_len(GeneConsensusLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GeneConsensus *]
 *
 */
GeneConsensus * GeneConsensus_alloc_std(void) 
{
    return GeneConsensus_alloc_len(GeneConsensusLISTLENGTH); 
}    


/* Function:  GeneConsensus_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [GeneConsensus *]
 *
 */
GeneConsensus * GeneConsensus_alloc_len(int len) 
{
    GeneConsensus * out;/* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = GeneConsensus_alloc()) == NULL)    
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->gsc = (GeneSingleCons ** ) ckcalloc (len,sizeof(GeneSingleCons *))) == NULL)    {  
      warn("Warning, ckcalloc failed in GeneConsensus_alloc_len");   
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_GeneConsensus(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GeneConsensus *]
 *
 * Return [UNKN ]  Undocumented return value [GeneConsensus *]
 *
 */
GeneConsensus * hard_link_GeneConsensus(GeneConsensus * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a GeneConsensus object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  GeneConsensus_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GeneConsensus *]
 *
 */
GeneConsensus * GeneConsensus_alloc(void) 
{
    GeneConsensus * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(GeneConsensus *) ckalloc (sizeof(GeneConsensus))) == NULL)  {  
      warn("GeneConsensus_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->center = 0; 
    out->gsc = NULL; 
    out->len = out->maxlen = 0;  


    return out;  
}    


/* Function:  free_GeneConsensus(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GeneConsensus *]
 *
 * Return [UNKN ]  Undocumented return value [GeneConsensus *]
 *
 */
GeneConsensus * free_GeneConsensus(GeneConsensus * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a GeneConsensus obj. Should be trappable"); 
      return NULL;   
      }  


#ifdef PTHREAD   
    assert(pthread_mutex_lock(&(obj->dynamite_mutex)) == 0); 
#endif   
    if( obj->dynamite_hard_link > 1)     {  
      return_early = 1;  
      obj->dynamite_hard_link--; 
      }  
#ifdef PTHREAD   
    assert(pthread_mutex_unlock(&(obj->dynamite_mutex)) == 0);   
#endif   
    if( return_early == 1)   
      return NULL;   
    if( obj->gsc != NULL)    {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->gsc[i] != NULL) 
          free_GeneSingleCons(obj->gsc[i]);  
        }  
      ckfree(obj->gsc);  
      }  


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_GeneFrequency21(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GeneFrequency21 *]
 *
 * Return [UNKN ]  Undocumented return value [GeneFrequency21 *]
 *
 */
GeneFrequency21 * hard_link_GeneFrequency21(GeneFrequency21 * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a GeneFrequency21 object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  GeneFrequency21_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GeneFrequency21 *]
 *
 */
GeneFrequency21 * GeneFrequency21_alloc(void) 
{
    GeneFrequency21 * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(GeneFrequency21 *) ckalloc (sizeof(GeneFrequency21))) == NULL)  {  
      warn("GeneFrequency21_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->ss5 = NULL; 
    out->ss3 = NULL; 
    /* codon[64] is an array: no default possible */ 
    /* central[4] is an array: no default possible */ 
    /* py[4] is an array: no default possible */ 
    /* spacer[4] is an array: no default possible */ 
    /* transition[GENEFREQUENCY21_TRANSITION_LEN] is an array: no default possible */ 
    /* cds_triplet[64] is an array: no default possible */ 


    return out;  
}    


/* Function:  free_GeneFrequency21(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GeneFrequency21 *]
 *
 * Return [UNKN ]  Undocumented return value [GeneFrequency21 *]
 *
 */
GeneFrequency21 * free_GeneFrequency21(GeneFrequency21 * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a GeneFrequency21 obj. Should be trappable");   
      return NULL;   
      }  


#ifdef PTHREAD   
    assert(pthread_mutex_lock(&(obj->dynamite_mutex)) == 0); 
#endif   
    if( obj->dynamite_hard_link > 1)     {  
      return_early = 1;  
      obj->dynamite_hard_link--; 
      }  
#ifdef PTHREAD   
    assert(pthread_mutex_unlock(&(obj->dynamite_mutex)) == 0);   
#endif   
    if( return_early == 1)   
      return NULL;   
    if( obj->ss5 != NULL)    
      free_GeneConsensus(obj->ss5);  
    if( obj->ss3 != NULL)    
      free_GeneConsensus(obj->ss3);  


    ckfree(obj); 
    return NULL; 
}    


/* Function:  replace_ss5_GeneFrequency21(obj,ss5)
 *
 * Descrip:    Replace member variable ss5
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [GeneFrequency21 *]
 * Arg:        ss5 [OWNER] New value of the variable [GeneConsensus *]
 *
 * Return [SOFT ]  member variable ss5 [boolean]
 *
 */
boolean replace_ss5_GeneFrequency21(GeneFrequency21 * obj,GeneConsensus * ss5) 
{
    if( obj == NULL)     {  
      warn("In replacement function ss5 for object GeneFrequency21, got a NULL object"); 
      return FALSE;  
      }  
    obj->ss5 = ss5;  
    return TRUE; 
}    


/* Function:  access_ss5_GeneFrequency21(obj)
 *
 * Descrip:    Access member variable ss5
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [GeneFrequency21 *]
 *
 * Return [SOFT ]  member variable ss5 [GeneConsensus *]
 *
 */
GeneConsensus * access_ss5_GeneFrequency21(GeneFrequency21 * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function ss5 for object GeneFrequency21, got a NULL object");    
      return NULL;   
      }  
    return obj->ss5;     
}    


/* Function:  replace_ss3_GeneFrequency21(obj,ss3)
 *
 * Descrip:    Replace member variable ss3
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [GeneFrequency21 *]
 * Arg:        ss3 [OWNER] New value of the variable [GeneConsensus *]
 *
 * Return [SOFT ]  member variable ss3 [boolean]
 *
 */
boolean replace_ss3_GeneFrequency21(GeneFrequency21 * obj,GeneConsensus * ss3) 
{
    if( obj == NULL)     {  
      warn("In replacement function ss3 for object GeneFrequency21, got a NULL object"); 
      return FALSE;  
      }  
    obj->ss3 = ss3;  
    return TRUE; 
}    


/* Function:  access_ss3_GeneFrequency21(obj)
 *
 * Descrip:    Access member variable ss3
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [GeneFrequency21 *]
 *
 * Return [SOFT ]  member variable ss3 [GeneConsensus *]
 *
 */
GeneConsensus * access_ss3_GeneFrequency21(GeneFrequency21 * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function ss3 for object GeneFrequency21, got a NULL object");    
      return NULL;   
      }  
    return obj->ss3;     
}    


/* Function:  replace_center_GeneConsensus(obj,center)
 *
 * Descrip:    Replace member variable center
 *             For use principly by API functions
 *
 *
 * Arg:           obj [UNKN ] Object holding the variable [GeneConsensus *]
 * Arg:        center [OWNER] New value of the variable [int]
 *
 * Return [SOFT ]  member variable center [boolean]
 *
 */
boolean replace_center_GeneConsensus(GeneConsensus * obj,int center) 
{
    if( obj == NULL)     {  
      warn("In replacement function center for object GeneConsensus, got a NULL object");    
      return FALSE;  
      }  
    obj->center = center;    
    return TRUE; 
}    


/* Function:  access_center_GeneConsensus(obj)
 *
 * Descrip:    Access member variable center
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [GeneConsensus *]
 *
 * Return [SOFT ]  member variable center [int]
 *
 */
int access_center_GeneConsensus(GeneConsensus * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function center for object GeneConsensus, got a NULL object");   
      return 0;  
      }  
    return obj->center;  
}    


/* Function:  access_gsc_GeneConsensus(obj,i)
 *
 * Descrip:    Access members stored in the gsc list
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the list [GeneConsensus *]
 * Arg:          i [UNKN ] Position in the list [int]
 *
 * Return [SOFT ]  Element of the list [GeneSingleCons *]
 *
 */
GeneSingleCons * access_gsc_GeneConsensus(GeneConsensus * obj,int i) 
{
    if( obj == NULL)     {  
      warn("In accessor function gsc for object GeneConsensus, got a NULL object");  
      return NULL;   
      }  
    if( obj->len <= i )  {  
      warn("In accessor function gsc for object GeneConsensus, index %%d is greater than list length %%d",i,obj->len);   
      return NULL;   
      }  
    return obj->gsc[i];  
}    


/* Function:  length_gsc_GeneConsensus(obj)
 *
 * Descrip:    discover the length of the list
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the list [GeneConsensus *]
 *
 * Return [UNKN ]  length of the list [int]
 *
 */
int length_gsc_GeneConsensus(GeneConsensus * obj) 
{
    if( obj == NULL)     {  
      warn("In length function gsc for object GeneConsensus, got a NULL object");    
      return -1;     
      }  
    return obj->len;     
}    


/* Function:  replace_string_GeneSingleCons(obj,string)
 *
 * Descrip:    Replace member variable string
 *             For use principly by API functions
 *
 *
 * Arg:           obj [UNKN ] Object holding the variable [GeneSingleCons *]
 * Arg:        string [OWNER] New value of the variable [char *]
 *
 * Return [SOFT ]  member variable string [boolean]
 *
 */
boolean replace_string_GeneSingleCons(GeneSingleCons * obj,char * string) 
{
    if( obj == NULL)     {  
      warn("In replacement function string for object GeneSingleCons, got a NULL object");   
      return FALSE;  
      }  
    obj->string = string;    
    return TRUE; 
}    


/* Function:  access_string_GeneSingleCons(obj)
 *
 * Descrip:    Access member variable string
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [GeneSingleCons *]
 *
 * Return [SOFT ]  member variable string [char *]
 *
 */
char * access_string_GeneSingleCons(GeneSingleCons * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function string for object GeneSingleCons, got a NULL object");  
      return NULL;   
      }  
    return obj->string;  
}    


/* Function:  replace_number_GeneSingleCons(obj,number)
 *
 * Descrip:    Replace member variable number
 *             For use principly by API functions
 *
 *
 * Arg:           obj [UNKN ] Object holding the variable [GeneSingleCons *]
 * Arg:        number [OWNER] New value of the variable [double]
 *
 * Return [SOFT ]  member variable number [boolean]
 *
 */
boolean replace_number_GeneSingleCons(GeneSingleCons * obj,double number) 
{
    if( obj == NULL)     {  
      warn("In replacement function number for object GeneSingleCons, got a NULL object");   
      return FALSE;  
      }  
    obj->number = number;    
    return TRUE; 
}    


/* Function:  access_number_GeneSingleCons(obj)
 *
 * Descrip:    Access member variable number
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [GeneSingleCons *]
 *
 * Return [SOFT ]  member variable number [double]
 *
 */
double access_number_GeneSingleCons(GeneSingleCons * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function number for object GeneSingleCons, got a NULL object");  
      return 0;  
      }  
    return obj->number;  
}    



#ifdef _cplusplus
}
#endif
