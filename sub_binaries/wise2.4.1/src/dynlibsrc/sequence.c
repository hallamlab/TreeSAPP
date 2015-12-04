#ifdef _cplusplus
extern "C" {
#endif
#include "wisebase.h"
#include "sequence.h"

/* Function:  new_Sequence_from_strings(name,seq)
 *
 * Descrip:    Makes a new sequence from strings given. 
 *             Separate memory will be allocated for them
 *             and them copied into it.
 *
 *             They can be NULL, in which case 
 *             o  a dummy name SequenceName will be assigned
 *             o  No sequence placed and length of zero.
 *
 *             Though this is dangerous later on. 
 *
 *             The sequence type is calculated automatically using
 *             /best_guess_type. If you want a DNA sequence but are
 *             unsure of the content of, for example, IUPAC codes,
 *             please use /force_to_dna_Sequence before using the
 *             sequence. Most of the rest of dynamite relies on a
 *             five letter A,T,G,C,N alphabet, but this function
 *             will allow any sequence type to be stored, so please
 *             check if you want to save yourself alot of grief.
 *
 *             In perl and other interfaces, this is a much safer
 *             constructor than the raw "new" type
 *
 *
 * Arg:        name [READ ] name of sequence, memory is allocated for it. [char *]
 * Arg:         seq [READ ] char * of sequence, memory is allocated for it. [char *]
 *
 * Return [UNKN ]  Undocumented return value [Sequence *]
 *
 */
# line 124 "sequence.dy"
Sequence * new_Sequence_from_strings(char * name,char * seq)
{
  Sequence * out;
  
  out = Sequence_alloc();
  
  if( name == NULL) 
    name = "SequenceName";

  out->name = stringalloc(name);

  if( seq == NULL ) {
    out->seq = NULL;
    out->len = 0;
    return out;
  }

  out->seq = stringalloc(seq);
  out->len = strlen(out->seq);
  out->offset = 1;
  out->end = out->len;

  out->type = best_guess_type(out);

  return out;
}


/* Function:  looks_like_accession(name)
 *
 * Descrip:    Returns true if name looks like [A-Za-z]+[0-9]+
 *             This should be an accession number 
 *
 *
 * Arg:        name [READ ] name to be tested [char *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 158 "sequence.dy"
boolean looks_like_accession(char * name)
{
  char * run;
  
  for(run=name;*run && isalpha((int)*run);run++)
    ;
  if( *run == '\0')
    return FALSE;
  for(;*run && isalnum((int)*run) && !isalpha((int)*run);run++)
    ;
  if( *run == '\0')
    return TRUE;
  return FALSE;
}


/* Function:  make_len_type_Sequence(seq)
 *
 * Descrip:    makes seq->len and seq->end match the seq->seq
 *             length number. 
 *
 *             It also checks the type of the sequence with
 *             /best_guess_type
 *
 *
 * Arg:        seq [RW   ] Sequence object [Sequence *]
 *
 */
# line 184 "sequence.dy"
void make_len_type_Sequence(Sequence * seq)
{
  seq->len = strlen(seq->seq);

  seq->end = seq->len + seq->offset -1;
  seq->type = best_guess_type(seq);

}

/* Function:  best_guess_type(seq)
 *
 * Descrip:    Guesses DNA or protein, by adding all
 *             the A,T,G,C up and if len < 300 && > 95% or 
 *             len > 300 && > 75% then considers
 *             it to be DNA. NB - Ns now counted.
 *
 *
 * Arg:        seq [READ ] Sequence to be guessed [Sequence *]
 *
 * Return [OWNER]  SEQUENCE_DNA or SEQUENCE_PROTEIN [int]
 *
 */
# line 202 "sequence.dy"
int   best_guess_type(Sequence * seq)
{
  register int i;
  int ch[26];
  int no;
  int unplaced = 0;
  int total;

  for(i=0;i<26;i++) 
    ch[i] = 0;
  
  for(i=0;i<seq->len;i++) {
    /*fprintf(stderr,"character is %c and total is %d\n",seq->seq[i],unplaced);*/
    if( (no=(int)( toupper(seq->seq[i])-'A')) > 26 || no < 0 )
      unplaced++;
    else ch[no]++;
  }
  
  total = (ch['A' - 'A']+ch['T' - 'A'] + ch['G' - 'A'] + ch['C' - 'A'] + ch['N' - 'A']);

  if( seq->len < 300 ) {
    if( (double)total/(double)seq->len > 0.95 )
      return SEQUENCE_DNA;
    else return SEQUENCE_PROTEIN;
  } else {
    if( (double)total/(double)seq->len > 0.75 )
      return SEQUENCE_DNA;
    else return SEQUENCE_PROTEIN;
  }

}

/* Function:  Sequence_type_to_string(type)
 *
 * Descrip:    Converts sequence type (SEQUENCE_*) to a string
 *
 *
 * Arg:        type [UNKN ] type eg SEQUENCE_PROTEIN [int]
 *
 * Return [UNKN ]  Undocumented return value [char *]
 *
 */
# line 239 "sequence.dy"
char * Sequence_type_to_string(int type)
{
  switch (type ) {
  case SEQUENCE_UNKNOWN : return "Unknown type";
  case SEQUENCE_PROTEIN : return "Protein";
  case SEQUENCE_DNA : return "Dna";
  case SEQUENCE_CDNA : return "cDNA";
  case SEQUENCE_GENOMIC : return "Genomic";
  case SEQUENCE_EST : return "Est";
  case SEQUENCE_RNA  : return "RNA";
  default : return "Undefined!";
  }
}

/* Function:  uppercase_Sequence(seq)
 *
 * Descrip:    makes all the sequence uppercase
 *
 *
 * Arg:        seq [RW   ] Sequence to be uppercased [Sequence *]
 *
 */
# line 259 "sequence.dy"
void uppercase_Sequence(Sequence * seq)
{
  int i;

  for(i=0;i<seq->len;i++)
    seq->seq[i] = toupper((int)seq->seq[i]);
}

/* Function:  force_to_dna_Sequence(seq,fraction,number_of_conver)
 *
 * Descrip:    This 
 *              a) sees how many non ATGCN characters there are in Seq
 *              b) If the level is below fraction
 *                 a) flips non ATGC chars to N
 *                 b) writes number of conversions to number_of_conver
 *                 c) returns TRUE
 *              c) else returns FALSE
 *
 *             fraction of 0.0 means completely intolerant of errors
 *             fraction of 1.0 means completely tolerant of errors
 *
 *
 *
 * Arg:                     seq [RW   ] sequence object read and converted  [Sequence *]
 * Arg:                fraction [READ ] number 0..1 for tolerance of conversion [double]
 * Arg:        number_of_conver [WRITE] number of conversions actually made [int *]
 *
 * Return [READ ]  TRUE for conversion to DNA, FALSE if not [boolean]
 *
 */
# line 286 "sequence.dy"
boolean force_to_dna_Sequence(Sequence * seq,double fraction,int * number_of_conver)
{
  int count =0;
  int i;

  if( seq == NULL ) {
    warn("Attempting to force a sequence with no Sequence object!\n");
    return FALSE;
  } 
  if( seq->len <= 0 ) {
    warn("Trying to make a sequence with a length of %d. Bad news!",seq->len);
    return FALSE;
  }

  
  for(i=0;i<seq->len;i++) {
    /* if it is lower case, uppercase it! */
    seq->seq[i] = (char)toupper((int)seq->seq[i]);
    
    if( !is_valid_base_char(seq->seq[i]) ) {
      count++;
    }
  }


  if( ((double)count/(double)seq->len) < fraction ) {
    seq->type = SEQUENCE_DNA;
    if( count != 0 ) {
      for(i=0;i<seq->len;i++) {
	if( !is_valid_base_char(seq->seq[i]) ) {
	  seq->seq[i] = 'N';
	}
      }
    }
    if( number_of_conver != NULL ) {
      *number_of_conver = count;
    }
    return TRUE;
  }
  else {
    if( number_of_conver != NULL ) {
      *number_of_conver = count;
    }
    return FALSE;
  }
}
    
/* Function:  is_reversed_Sequence(seq)
 *
 * Descrip:    Currently the sequence object stores 
 *             reversed sequences as start > end.
 *
 *             This tests that and returns true if it is
 *
 *
 * Arg:        seq [READ ] sequence to test [Sequence *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 342 "sequence.dy"
boolean is_reversed_Sequence(Sequence * seq)
{
  if( seq->offset > seq->end)
    return TRUE;
  return FALSE;
}

/* Function:  translate_Sequence(dna,ct)
 *
 * Descrip:    This translates a DNA sequence to a protein.
 *             It assummes that it starts at first residue
 *             (use trunc_Sequence to chop a sequence up).
 *
 *
 * Arg:        dna [READ ] DNA sequence to be translated [Sequence *]
 * Arg:         ct [READ ] Codon table to do codon->aa mapping [CodonTable *]
 *
 * Return [OWNER]  new protein sequence [Sequence *]
 *
 */
# line 359 "sequence.dy"
Sequence * translate_Sequence(Sequence * dna,CodonTable * ct)
{
  Sequence * out;
  int i;
  int j;
  int len;
  char * seq;
  char * name;
  char buffer[512];

  if( is_dna_Sequence(dna) == FALSE) {
    warn("Trying to make a translation from a non DNA sequence... type is [%s]",Sequence_type_to_string(dna->type));
    return NULL;
  }

  len = dna->len/3 + 1;
  
  seq = ckcalloc(len,sizeof(char));

  sprintf(buffer,"%s.tr",dna->name == NULL ? "NoNameDNASeq" : dna->name);

  name = stringalloc(buffer);

  out  = Sequence_from_dynamic_memory(name,seq);

  for(i=0,j=0;i<dna->len-3;i+=3,j++) {
    out->seq[j] = aminoacid_from_seq(ct,dna->seq+i);
  }
  out->seq[j] = '\0';

  out->type  = SEQUENCE_PROTEIN;
  out->len   = strlen(out->seq);

  return out;
}

/* Function:  reverse_complement_Sequence(seq)
 *
 * Descrip:    This both complements and reverses a sequence,
 *             - a common wish!
 *
 *             The start/end are correct with respect to the start/end
 *             of the sequence (ie start = end, end = start).
 *
 *
 * Arg:        seq [READ ] Sequence to that is used to reverse (makes a new Sequence) [Sequence *]
 *
 * Return [OWNER]  new Sequence which is reversed [Sequence *]
 *
 */
# line 406 "sequence.dy"
Sequence * reverse_complement_Sequence(Sequence * seq)
{
  Sequence * out;
  int i;
  int j;


  if( is_dna_Sequence(seq) != TRUE ) {
    warn("Cannot reverse complement non-DNA sequence... type is %s",Sequence_type_to_string(seq->type));
    return NULL;
  }

  out = Sequence_from_static_memory(seq->name,seq->seq);
  
  for(j=0,i=seq->len-1;i >= 0;i--,j++) {
    if( seq->seq[i] == '.' || seq->seq[i] == '-' ) {
      out->seq[j] = seq->seq[i];
    } else {
      out->seq[j] = char_complement_base(seq->seq[i]);
    }
    /*fprintf(stderr,"In position %d placed %c from %c\n",j,out->seq[j],seq->seq[i]);*/
  }

  out->len = strlen(seq->seq);

  out->offset = seq->end;
  out->end   = seq->offset;
  out->type  = seq->type;
  

  return out;
}
/* Function:  magic_trunc_Sequence(seq,start,end)
 *
 * Descrip:    Clever function for dna sequences.
 *
 *             When start < end, truncates normally
 *
 *             when start > end, truncates end,start and then
 *             reverse complements.
 *
 *             ie. If you have a coordinate system where reverse 
 *             sequences are labelled in reverse start/end way,
 *             then this routine produces the correct sequence.
 *
 *
 * Arg:          seq [READ ] sequence that is the source to be truncated [Sequence *]
 * Arg:        start [READ ] start point [int]
 * Arg:          end [READ ] end point [int]
 *
 * Return [OWNER]  new Sequence which is truncated/reversed [Sequence *]
 *
 */
# line 456 "sequence.dy"
Sequence * magic_trunc_Sequence(Sequence * seq,int start,int end)
{
  Sequence * temp;
  Sequence * out;

  if( is_dna_Sequence(seq) == FALSE) {
    warn("Cannot magic truncate on a non DNA sequence... type is %s",Sequence_type_to_string(seq->type));
    return NULL;
  }

  if( start < 0 || end < 0 ) {
    warn("Attempting a magic truncation on indices which are less than zero [%d:%d]. Clearly impossible",start,end);
    return NULL;
  }

  if( start < end ) {
    return trunc_Sequence(seq,start,end);
  }
  else {
    temp = trunc_Sequence(seq,end,start);
    if( temp == NULL ) {
      warn("Unable to truncate sequence");
      return NULL;
    }


    out = reverse_complement_Sequence(temp);

    free_Sequence(temp);

    return out;
  }

}


/* Function:  trunc_Sequence(seq,start,end)
 *
 * Descrip:    truncates a sequence. It produces a new memory structure
 *             which is filled from sequence start to end.
 *
 *             Please notice
 *               
 *               Truncation is in C coordinates. That is
 *             the first residue is 0 and end is the number of the
 *             residue after the cut-point. In otherwords to 
 *             2 - 3 would be a single residue truncation. So - if
 *             you want to work in more usual, 'inclusive' molecular
 *             biology numbers, which start at 1, then you need to say
 *
 *               trunc_Sequence(seq,start-1,end);
 *
 *             (NB, should be (end - 1 + 1) = end for the last coordinate).
 *
 *               Truncation occurs against the *absolute* coordinate
 *             system of the Sequence, not the offset/end pair inside.
 *             So, this is a very bad error
 *              
 *               ** wrong code, and also leaks memory **
 *
 *               tru = trunc_Sequence(trunc_Sequence(seq,50,80),55,75); 
 *
 *             This the most portable way of doing this
 *
 *               temp = trunc_Sequence(seq,50,80);
 *
 *               tru  = trunc_Sequence(temp,55-temp->offset,75-temp->offset);
 *
 *               free_Sequence(temp);
 *
 *
 *
 * Arg:          seq [READ ] object holding the sequence to be truncated [Sequence *]
 * Arg:        start [READ ] start point of truncation [int]
 * Arg:          end [READ ] end point of truncation [int]
 *
 * Return [OWNER]  newly allocated sequence structure [Sequence *]
 *
 */
# line 532 "sequence.dy"
Sequence * trunc_Sequence(Sequence * seq,int start,int end)
{
  char * name;
  char * seqb;
  Sequence * out;

  if( start < 0 || end < 0 ) {
    warn("Attempting a truncation on indices which are less than zero [%d:%d]. Clearly impossible",start,end);
    return NULL;
  }
  
  if( end <= start ) {
    warn("Trying to truncate Sequence from %d - %d",start,end);
    return NULL;
  }
  
  if( end > seq->len ) {
    warn("Trying to truncate Sequecne %s from %d - %d when length is %d",
	 seq->name,start,end,seq->len);
    return NULL;
    }
  
  name = stringalloc(seq->name);
  
  seqb = (char *) ckcalloc(end-start+1,sizeof(char));
  
  memcpy(seqb,seq->seq+start,(end-start));
  seqb[end-start] = '\0';
  
  out = Sequence_from_dynamic_memory(name,seqb);

  out->len = strlen(out->seq);
  out->type = seq->type;
  out->offset = seq->offset+start;
  out->end  = seq->offset + end-1;

  return out;
}

/* Function:  read_SRS_db_Sequence(datastring,srsstring)
 *
 * Descrip:    A function for you to easily specify the sequence name
 *             and the database separately. Just concatonates the two
 *             strings with : betwqeen them. Therefore you should use
 *             "swisprot-id" for example as your datastring.
 *
 *             calls /read_SRS_Sequence
 *
 *
 * Arg:        datastring [READ ] string representing the database (swissprot-id) [char *]
 * Arg:         srsstring [READ ] string for the name (eg, ROA1_HUMAN) [char *]
 *
 * Return [UNKN ]  Undocumented return value [Sequence *]
 *
 */
# line 582 "sequence.dy"
Sequence * read_SRS_db_Sequence(char * datastring,char * srsstring)
{
  char buffer[256];
  
  sprintf(buffer,"%s:%s",datastring,srsstring);
  
  return read_SRS_Sequence(buffer);
}

/* Function:  read_SRS_Sequence(srsstring)
 *
 * Descrip:    reads SRS specified sequence. calls popoen
 *             with getz -f using srs4 syntax. Will only read
 *             the first sequence if there is more than one in the 
 *             SRS spec, and does not warn you about additional 
 *             sequences
 *
 *
 * Arg:        srsstring [READ ] srs spec'd string swissprot-id:ROA1_HUMAN [char *]
 *
 * Return [UNKN ]  Undocumented return value [Sequence *]
 *
 */
# line 600 "sequence.dy"
Sequence * read_SRS_Sequence(char * srsstring)
{
  FILE * pipe;
  char buffer[MAXLINE];
  Sequence * out;
  
  
  sprintf(buffer,"getz -d '[%s]' ",srsstring);
  
  pipe = popen(buffer,"r");
  
  if ( pipe == NULL ) {
    warn("Could not open %s as an SRS database string - probably no getz",srsstring);
    return NULL;
  }
  
  out = read_fasta_Sequence(pipe);
  
  pclose(pipe);
  
  return out;
}

/* Function:  read_efetch_Sequence(efetchstring)
 *
 * Descrip:    reads efetch specificed sequence. calls popen to
 *             efetch. A hack around accession numbers so that if the 
 *             thing looks like WP:acc number, calls it with -a...
 *             otherwise assummes you have both database and name in the
 *             efetchstring
 *
 *
 * Arg:        efetchstring [READ ] efetch valid string [char *]
 *
 * Return [UNKN ]  Undocumented return value [Sequence *]
 *
 */
# line 632 "sequence.dy"
Sequence * read_efetch_Sequence(char * efetchstring)
{
  FILE * pf;
  Sequence * out;
  char buffer[MAXLINE];

  if( strstartcmp(efetchstring,"WP:") != 0 && looks_like_accession(efetchstring+3) == TRUE) {
    sprintf(buffer,"efetch -f -a %s",efetchstring);
  }
  else {
    sprintf(buffer,"efetch -f %s",efetchstring);
  }

  pf = popen(buffer,"r");

  if( pf == NULL ) {
    warn("Could not open efetch pipe with [%s]",efetchstring);
    return NULL;
  }

  out = read_fasta_Sequence(pf);

  pclose(pf);

  return out;
}

/* Function:  read_fasta_file_Sequence(filename)
 *
 * Descrip:    Just a call
 *               a) open filename
 *               b) read sequence with /read_fasta_Sequence
 *               c) close file.
 *
 *
 * Arg:        filename [READ ] filename to open  [char *]
 *
 * Return [UNKN ]  Undocumented return value [Sequence *]
 *
 */
# line 667 "sequence.dy"
Sequence * read_fasta_file_Sequence(char * filename)
{
  Sequence * out;
  FILE * ifp;
  
  ifp = openfile(filename,"r");
  
  if( ifp == NULL ) {
    warn("Cannot open %s for read_fasta_file",filename);
    return NULL;
  }


  out = read_fasta_Sequence(ifp);
  
  fclose(ifp);
  
  return out;
}

/* Function:  read_Sequence_EMBL_seq(buffer,maxlen,ifp)
 *
 * Descrip:    reads the sequence part of an EMBL file.
 *
 *             This function can either take a file which 
 *             starts
 *
 *
 *
 * Arg:        buffer [RW   ] buffer containing the first line. [char *]
 * Arg:        maxlen [READ ] length of buffer [int]
 * Arg:           ifp [READ ] input file to read from [FILE *]
 *
 * Return [UNKN ]  Undocumented return value [Sequence *]
 *
 */
# line 698 "sequence.dy"
Sequence * read_Sequence_EMBL_seq(char * buffer,int maxlen,FILE * ifp)
{
  Sequence * out;
  char seqbuffer[SEQUENCEBLOCK];
  int i = 0;
  signed char c;

  if( !isalpha((int)*buffer) ) {
    warn("I don't like this - got a buffer of [%s] in reading an EMBL sequence section",buffer);
  }

  do {
    if( strstartcmp(buffer,"SQ") == 0 ) {
      break;
    }
  } while ( fgets(buffer,maxlen,ifp) != NULL);

  out = empty_Sequence_from_dynamic_memory(stringalloc("EMBLseq"));

  while( (c=fgetc(ifp)) != EOF ) {
    if( c == '/' && (c=fgetc(ifp)) == '/') 
      break; /*** ugly perhaps. what about single / lines? ***/

    if( isalpha(c) )
      seqbuffer[i++] = c;
    if( i > SEQUENCEBLOCK-2) {
      seqbuffer[i] = '\0';
      if( add_string_to_Sequence(out,seqbuffer) == FALSE )
	{
	  warn("Could not read full sequence of %s - returning\n",out->name);
	  return out;
	}
      i = 0;
    }
  }
  
  /* ok have to now put away final buffer read! */
  
  seqbuffer[i] = '\0';
  
  add_string_to_Sequence(out,seqbuffer);

  /** add back > if need be **/
  
  if( feof(ifp) || c != '/' ) {
    warn("In parsing an EMBL file got an poor ending of a sequence region");
  } else {
    while( (c=fgetc(ifp)) != '\n' && c != EOF )
      ;
  }
  
  make_len_type_Sequence(out);

  return out;
}

/* Function:  read_fasta_Sequence(ifp)
 *
 * Descrip:    reads a fasta file assumming pretty 
 *             standard sanity for the file layout
 *
 *
 * Arg:        ifp [UNKN ] Undocumented argument [FILE *]
 *
 * Return [UNKN ]  Undocumented return value [Sequence *]
 *
 */
# line 758 "sequence.dy"
Sequence * read_fasta_Sequence(FILE * ifp)
{
  char buffer[MAXLINE];
  char c; /* we need to do a nasty fgetc/ungetc to detect whether we have ended */
  char seqbuffer[SEQUENCEBLOCK];
  int i = 0;
  int j;
  int tax_id = -1;
  double weight = 1.0;
  boolean seen_weight = 0;
  char * temp;
  char * desc;
  char * name;
  Sequence * out;


  

  /* first line is '>' */

  if( fgets(buffer,MAXLINE,ifp) == NULL ) {
    return NULL;
  }

  if( buffer[strlen(buffer)-1] != '\n' ) {
    /* excessive description line */
    while( (c=fgetc(ifp)) != '\n' && c != EOF ) {
      ;
    }
    if( c == EOF ) {
      warn("In sequence file, overead line into end of file. Exiting");
      return NULL;
    }
  }


  if( (temp = strstr(buffer,"tax_id=")) != NULL ) {
    tax_id = atoi(temp+strlen("tax_id="));
  }

  if( (temp = strstr(buffer,"weight=")) != NULL ) {
    weight = atof(temp+strlen("weight="));
  }

  if( buffer[0] == '/' && buffer[1] == '/' ) {
    /* silent exit as if EOF */
    return NULL;
  }
  
  if( buffer[0] != '>' ) {
    warn("First character read not >, assumming is not fasta");
    return NULL;
  }

  /* some annoying formats have >\s+\S+ ... Grrrr! */

  for(i=1,name=buffer+1;*name && isspace(*name);name++,i++)
    ;
    
  /* now delimit this point */
  for(;!isspace(buffer[i]);i++)
    ;

  buffer[i] = '\0';

  /* now find description line if here */

  for(i++;isspace(buffer[i]) && buffer[i] != '\0';i++) {
    ;
  }
  if( buffer[i] != '\0' ) {
    desc = buffer+i;
    for(i++;buffer[i] != '\0' && buffer[i] != '\n';i++) {
      ;
    }
    buffer[i] = '\0';
  } else {
    desc = NULL;
  }


  out = empty_Sequence_from_dynamic_memory(stringalloc(name));
  if( desc != NULL ){
    out->desc = stringalloc(desc);
  }

  out->tax_id = tax_id;
  out->weight = weight;

  i =0;
  while( 1 ) {
    i =0;
    /* get/ungetc to see if we should read this line */
    c = fgetc(ifp);
    if( c == EOF || c == '>' || c == '/') {
      break;
    }
    seqbuffer[i++] = c;
    /* get the rest of the line */
    if( fgets(buffer,MAXLINE,ifp) == NULL ) {
      warn("Strangely truncated line in fasta file");
      break;
    }

    /*fprintf(stderr,"Getting into while loop.. with [%s]\n",buffer); */

    for(j=0;buffer[j] != '\0';j++) {
      if( isalpha(buffer[j]) || buffer[j] == '-' || buffer[j] == '.' || buffer[j] == '~' )
	seqbuffer[i++] = buffer[j];
    }
    seqbuffer[i] = '\0';

    if( add_string_to_Sequence(out,seqbuffer) == FALSE ) {
      warn("Could not read full sequence of %s - returning\n",out->name);
      return out;
    }
  }

  /* ok have to now put away final buffer read! */
  
  seqbuffer[i] = '\0';
  
  add_string_to_Sequence(out,seqbuffer);

  /* end of file or '>' or '/' */
  if( c == '>' || c == '/' )
    ungetc(c,ifp);

  make_len_type_Sequence(out);

  return out;
      
}



  
/* Function:  read_old_fasta_Sequence(ifp)
 *
 * Descrip:    reads the fasta file: format is
 *
 *             >name
 *             sequence
 *
 *             allocates a structure and puts in the
 *             sequence. Calls /make_len_type_Sequence to
 *             check type and length.
 *
 *             It leaves the '>' on the next fasta sequence
 *             for multiple sequence reading
 *
 *
 * Arg:        ifp [READ ] input file to read from [FILE *]
 *
 * Return [OWNER]  new Sequence structure  [Sequence *]
 *
 */
# line 911 "sequence.dy"
Sequence * read_old_fasta_Sequence(FILE * ifp)
{
  Sequence * out;
  char seqbuffer[SEQUENCEBLOCK];
  int i = 0;
  signed char c = EOF;
  
  if( feof(ifp) ) {
    /* fail silently */
    return NULL;
  }

  while( (c=fgetc(ifp)) != EOF &&  isspace(c)) 
    ;

  if( feof(ifp) )
    return NULL;

  if( c != '>' ) {
    warn("First letter read is not '>' - assumming it is not a fasta stream");
    return NULL;
  }

  
  if( c == EOF || feof(ifp) )
    return NULL;
  
  /*** ok = got to > ****/
  /*** read in name  ****/
  
  
  while( !isspace(c=fgetc(ifp)) && c != EOF)
    seqbuffer[i++]=c;
  
  if( c == EOF)
    return NULL;
  
  seqbuffer[i]='\0';
  
  /*** name now in sequence buffer - make sequence ***/
  
  out = empty_Sequence_from_dynamic_memory(stringalloc(seqbuffer));
  
  if( out == NULL )
    return NULL;


  /*** ok, suck in the rest of this line if necessary (ie something else on the 1st line) ***/

  while( c != EOF && c != '\n' )
    c=fgetc(ifp);

  
  /*** now read in sequence ***/
  
  for(i=0; !feof(ifp) && (c=fgetc(ifp)) != '>' && c != EOF;)
    {
      if( isalpha(c) || c == '-' || c == '.' || c == '~' )
	seqbuffer[i++] = c;
      if( i > SEQUENCEBLOCK-2)
	{
	  seqbuffer[i] = '\0';
	  if( add_string_to_Sequence(out,seqbuffer) == FALSE )
	    {
	      warn("Could not read full sequence of %s - returning\n",out->name);
	      return out;
	    }
	  i = 0;
	}
    }


  
  
  /* ok have to now put away final buffer read! */
  
  seqbuffer[i] = '\0';
  
  add_string_to_Sequence(out,seqbuffer);

  /** add back > if need be **/
  
  if( c == '>' )
    ungetc(c,ifp);
  
  make_len_type_Sequence(out);

  return out;
}

/* Function:  show_Sequence_residue_list(seq,start,end,ofp)
 *
 * Descrip:    shows a region of a sequence as
 *                124  A
 *                125  T
 *
 *             etc from start to end. The numbers
 *             are in C coordinates (ie, 0 is the first
 *             letter).
 *
 *             useful for debugging
 *
 *
 * Arg:          seq [READ ] Sequence to show [Sequence *]
 * Arg:        start [READ ] start of list [int]
 * Arg:          end [READ ] end of list [int]
 * Arg:          ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 1017 "sequence.dy"
void show_Sequence_residue_list(Sequence * seq,int start,int end,FILE * ofp)
{
  int i;
  
  for(i=start;i<end;i++) {
    fprintf(ofp,"%d %c\n",i,seq->seq[i]);
  }
}

/* Function:  add_string_to_Sequence(seq,more)
 *
 * Descrip:    New add_string_to_Sequence which should be much better at
 *             handling memory update
 *
 *
 * Arg:         seq [UNKN ] Undocumented argument [Sequence *]
 * Arg:        more [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 1030 "sequence.dy"
boolean add_string_to_Sequence(Sequence * seq,char * more)
{
  int strl;

  assert(seq);
  assert(more);
  assert(seq->maxlen > 0);

  strl = strlen(more);

  if( (seq->len + strl) > (seq->maxlen -2) ) {
    if( seq->maxlen < SEQUENCE_REALLOC_LINEAR ) {
      seq->seq = realloc( seq->seq,seq->maxlen*2);
      seq->maxlen = seq->maxlen *2;
    } else {
      seq->seq = realloc( seq->seq,seq->maxlen + SEQUENCE_REALLOC_LINEAR );
      seq->maxlen = seq->maxlen + SEQUENCE_REALLOC_LINEAR;
    }
  }
  assert(seq->seq);

  strcat(seq->seq,more);

  seq->len = strlen(seq->seq);

  return TRUE;
}


/* Function:  add_string_to_Sequence_old(seq,more)
 *
 * Descrip:    Dodgy function. This is meant to add the more
 *             sequence to seq (into ->seq). Not sure how stable
 *             this is. In theory it reallocates memory on
 *             the basis of ->maxlen.
 *
 *
 * Arg:         seq [UNKN ] Sequence to add sequence to [Sequence *]
 * Arg:        more [UNKN ] pointer to sequence to add [char *]
 *
 * Return [UNKN ]  TRUE if successful, FALSE if not [boolean]
 *
 */
# line 1070 "sequence.dy"
boolean add_string_to_Sequence_old(Sequence * seq,char * more)
{
  register int len;
  register int blocklen;
  void * temp;
  
  
  len = strlen(more)+1;
  
  if( len < seq->maxlen - seq->len )
    {
      /*** ok can add to this block! ****/
      strcat(seq->seq,more);
      seq->len = strlen(seq->seq);
      return TRUE;
    }
  
  
  /*** nope - need to realloc ****/
  
  len -= (seq->maxlen - seq->len); /* amount that needs to be realloc'd */
  blocklen = 1 + (int)(len / SEQUENCEBLOCK); /* number of blocks */
  blocklen *= SEQUENCEBLOCK;                 /* make that into bytes */
  blocklen += seq->maxlen;    /* final size of string */
  
  temp = ckrealloc ( seq->seq,blocklen);             
  
  
  if( temp == NULL )
    {
      warn("Sequence block error for sequence %s on blocklen %d\n",CKS(seq->name),blocklen);
      return FALSE;
    }
  
  seq->seq = (char *) temp; /* realloc moves the memory for us as well */
  seq->maxlen = blocklen;
  
  /*** copy in string ****/
  
  strcat(seq->seq,more);
  
  seq->len = strlen(seq->seq);
  
  return TRUE;
}

/* Function:  empty_Sequence_from_dynamic_memory(name)
 *
 * Descrip:    Only allocates sequence structure and name
 *
 *
 * Arg:        name [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [Sequence *]
 *
 */
# line 1119 "sequence.dy"
Sequence * empty_Sequence_from_dynamic_memory(char * name)
{
  Sequence * out;
  
  out = Sequence_alloc();
  
  if( out == NULL )
    return NULL;
  
  if( name == NULL )
    {
      warn("Attempting to make an empty sequence with no name: assigning dummy name");
      name = stringalloc("DummyName");
    }
  
  out->name = name;
  out->seq = (char *) ckcalloc (SEQUENCE_STARTSIZE,sizeof(char));
  out->maxlen = SEQUENCE_STARTSIZE;
  out->len = 0;
  
  return out;
}

/* Function:  Sequence_alloc_len(len)
 *
 * Descrip:    allocates sequence structure with enough
 *             length in char for len sequence.
 *
 *
 * Arg:        len [READ ] length of blank sequene space [int]
 *
 * Return [UNKN ]  Undocumented return value [Sequence *]
 *
 */
# line 1148 "sequence.dy"
Sequence * Sequence_alloc_len(int len)
{
  Sequence * out;

  out = Sequence_alloc();
  if( out == NULL)
    return NULL;

  out->seq = (char *) ckcalloc (len,sizeof(char));
  out->maxlen = out->len = len;

  return out;
}

/* Function:  Sequence_from_static_memory (name,seq)
 *
 * Descrip:    Allocates the sequence structure and memory for
 *             name and seq, copies them in.
 *
 *
 * Arg:        name [UNKN ] Undocumented argument [char *]
 * Arg:         seq [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [Sequence *]
 *
 */
# line 1167 "sequence.dy"
Sequence * Sequence_from_static_memory (char * name,char * seq)
{
  return Sequence_from_dynamic_memory(stringalloc(name),stringalloc(seq));
}

/* Function:  Sequence_from_dynamic_memory(name,seq)
 *
 * Descrip:    Allocates the sequence structure and simple attaches
 *             name and seq to the correct places
 *
 *
 * Arg:        name [UNKN ] name of sequence [char *]
 * Arg:         seq [UNKN ] a char * to correct sequence [char *]
 *
 * Return [UNKN ]  Undocumented return value [Sequence *]
 *
 */
# line 1180 "sequence.dy"
Sequence * Sequence_from_dynamic_memory(char * name,char * seq)
{
  Sequence * out;
  
  
  if( seq == NULL) {
      warn("Cannot make a sequence with no sequence!");
      return NULL;
    }
  
  if( name == NULL ) {
      warn("You are attempting to make a sequence with no name - assigning dummy name");
      name = stringalloc ("DummyName");
    }
  
  out = Sequence_alloc();
  
  if( out == NULL)
    return out;
  
  out->name = name;
  out->seq = seq;
  
  out->maxlen = out->len = strlen(seq);
  
  return out;
}

/* Function:  write_fasta_Sequence(seq,ofp)
 *
 * Descrip:    writes a fasta file of the form
 *             >name
 *             Sequence
 *
 *
 * Arg:        seq [READ ] sequence to be written [Sequence *]
 * Arg:        ofp [UNKN ] file to write to [FILE *]
 *
 */
# line 1217 "sequence.dy"
void write_fasta_Sequence(Sequence * seq,FILE * ofp)
{
  assert(seq);
  fprintf(ofp,">%s\n",seq->name == NULL ? "NoName_Null_string" : seq->name );
  show_line(seq->seq,60,ofp);
}

/* Function:  read_fasta_SequenceSet(ifp)
 *
 * Descrip:    reads a fasta file as a sequence set
 *
 *
 * Arg:        ifp [UNKN ] Undocumented argument [FILE *]
 *
 * Return [UNKN ]  Undocumented return value [SequenceSet *]
 *
 */
# line 1227 "sequence.dy"
SequenceSet * read_fasta_SequenceSet(FILE * ifp)
{
  SequenceSet * out;
  Sequence * in;

  out = SequenceSet_alloc_std();

  while( (in = read_fasta_Sequence(ifp)) != NULL ) {
    add_SequenceSet(out,in);
  }

  return out;
}

/* Function:  read_fasta_SequenceSet_file(filename)
 *
 * Descrip:    opens file and reads in sequence set
 *
 *
 * Arg:        filename [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [SequenceSet *]
 *
 */
# line 1244 "sequence.dy"
SequenceSet * read_fasta_SequenceSet_file(char * filename)
{
  SequenceSet * out;
  FILE * ifp;

  ifp = openfile(filename,"r");
  if( ifp == NULL ) {
    warn("Could not open %s as a file",filename);
    return NULL;
  }

  out = read_fasta_SequenceSet(ifp);

  fclose(ifp);

  return out;

}


# line 1274 "sequence.c"
/* Function:  hard_link_Sequence(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [Sequence *]
 *
 * Return [UNKN ]  Undocumented return value [Sequence *]
 *
 */
Sequence * hard_link_Sequence(Sequence * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a Sequence object: passed a NULL object");    
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  Sequence_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Sequence *]
 *
 */
Sequence * Sequence_alloc(void) 
{
    Sequence * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(Sequence *) ckalloc (sizeof(Sequence))) == NULL)    {  
      warn("Sequence_alloc failed ");    
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->name = NULL;    
    out->seq = NULL; 
    out->len = 0;    
    out->maxlen = 0; 
    out->offset = 1; 
    out->end = (-1); 
    out->type = SEQUENCE_UNKNOWN;    
    out->tax_id = 0; 
    out->weight = 1.0;   
    out->desc = NULL;    


    return out;  
}    


/* Function:  free_Sequence(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [Sequence *]
 *
 * Return [UNKN ]  Undocumented return value [Sequence *]
 *
 */
Sequence * free_Sequence(Sequence * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a Sequence obj. Should be trappable");  
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
    if( obj->name != NULL)   
      ckfree(obj->name);     
    if( obj->seq != NULL)    
      ckfree(obj->seq);  
    if( obj->desc != NULL)   
      ckfree(obj->desc);     


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_SequenceSet(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_SequenceSet
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [Sequence **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_SequenceSet(Sequence ** list,int i,int j)  
{
    Sequence * temp; 
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_SequenceSet(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_SequenceSet which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [Sequence **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_SequenceSet(Sequence ** list,int left,int right,int (*comp)(Sequence * ,Sequence * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_SequenceSet(list,left,(left+right)/2);  
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_SequenceSet (list,++last,i);    
      }  
    swap_SequenceSet (list,left,last);   
    qsort_SequenceSet(list,left,last-1,comp);    
    qsort_SequenceSet(list,last+1,right,comp);   
}    


/* Function:  sort_SequenceSet(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_SequenceSet
 *
 *
 * Arg:         obj [UNKN ] Object containing list [SequenceSet *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_SequenceSet(SequenceSet * obj,int (*comp)(Sequence *, Sequence *)) 
{
    qsort_SequenceSet(obj->set,0,obj->len-1,comp);   
    return;  
}    


/* Function:  expand_SequenceSet(obj,len)
 *
 * Descrip:    Really an internal function for add_SequenceSet
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [SequenceSet *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_SequenceSet(SequenceSet * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_SequenceSet called with no need");    
      return TRUE;   
      }  


    if( (obj->set = (Sequence ** ) ckrealloc (obj->set,sizeof(Sequence *)*len)) == NULL)     {  
      warn("ckrealloc failed for expand_SequenceSet, returning FALSE");  
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_SequenceSet(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [SequenceSet *]
 * Arg:        add [OWNER] Object to add to the list [Sequence *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_SequenceSet(SequenceSet * obj,Sequence * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_SequenceSet(obj,obj->len + SequenceSetLISTLENGTH) == FALSE) 
        return FALSE;    
      }  


    obj->set[obj->len++]=add;    
    return TRUE; 
}    


/* Function:  flush_SequenceSet(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [SequenceSet *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_SequenceSet(SequenceSet * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->set[i] != NULL)   {  
        free_Sequence(obj->set[i]);  
        obj->set[i] = NULL;  
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  SequenceSet_alloc_std(void)
 *
 * Descrip:    Equivalent to SequenceSet_alloc_len(SequenceSetLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SequenceSet *]
 *
 */
SequenceSet * SequenceSet_alloc_std(void) 
{
    return SequenceSet_alloc_len(SequenceSetLISTLENGTH); 
}    


/* Function:  SequenceSet_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [SequenceSet *]
 *
 */
SequenceSet * SequenceSet_alloc_len(int len) 
{
    SequenceSet * out;  /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = SequenceSet_alloc()) == NULL)  
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->set = (Sequence ** ) ckcalloc (len,sizeof(Sequence *))) == NULL)    {  
      warn("Warning, ckcalloc failed in SequenceSet_alloc_len"); 
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_SequenceSet(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [SequenceSet *]
 *
 * Return [UNKN ]  Undocumented return value [SequenceSet *]
 *
 */
SequenceSet * hard_link_SequenceSet(SequenceSet * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a SequenceSet object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  SequenceSet_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SequenceSet *]
 *
 */
SequenceSet * SequenceSet_alloc(void) 
{
    SequenceSet * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(SequenceSet *) ckalloc (sizeof(SequenceSet))) == NULL)  {  
      warn("SequenceSet_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->set = NULL; 
    out->len = out->maxlen = 0;  


    return out;  
}    


/* Function:  free_SequenceSet(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [SequenceSet *]
 *
 * Return [UNKN ]  Undocumented return value [SequenceSet *]
 *
 */
SequenceSet * free_SequenceSet(SequenceSet * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a SequenceSet obj. Should be trappable");   
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
    if( obj->set != NULL)    {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->set[i] != NULL) 
          free_Sequence(obj->set[i]);    
        }  
      ckfree(obj->set);  
      }  


    ckfree(obj); 
    return NULL; 
}    


/* Function:  replace_name_Sequence(obj,name)
 *
 * Descrip:    Replace member variable name
 *             For use principly by API functions
 *
 *
 * Arg:         obj [UNKN ] Object holding the variable [Sequence *]
 * Arg:        name [OWNER] New value of the variable [char *]
 *
 * Return [SOFT ]  member variable name [boolean]
 *
 */
boolean replace_name_Sequence(Sequence * obj,char * name) 
{
    if( obj == NULL)     {  
      warn("In replacement function name for object Sequence, got a NULL object");   
      return FALSE;  
      }  
    obj->name = name;    
    return TRUE; 
}    


/* Function:  access_name_Sequence(obj)
 *
 * Descrip:    Access member variable name
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [Sequence *]
 *
 * Return [SOFT ]  member variable name [char *]
 *
 */
char * access_name_Sequence(Sequence * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function name for object Sequence, got a NULL object");  
      return NULL;   
      }  
    return obj->name;    
}    


/* Function:  replace_seq_Sequence(obj,seq)
 *
 * Descrip:    Replace member variable seq
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [Sequence *]
 * Arg:        seq [OWNER] New value of the variable [char *]
 *
 * Return [SOFT ]  member variable seq [boolean]
 *
 */
boolean replace_seq_Sequence(Sequence * obj,char * seq) 
{
    if( obj == NULL)     {  
      warn("In replacement function seq for object Sequence, got a NULL object");    
      return FALSE;  
      }  
    obj->seq = seq;  
    return TRUE; 
}    


/* Function:  access_seq_Sequence(obj)
 *
 * Descrip:    Access member variable seq
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [Sequence *]
 *
 * Return [SOFT ]  member variable seq [char *]
 *
 */
char * access_seq_Sequence(Sequence * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function seq for object Sequence, got a NULL object");   
      return NULL;   
      }  
    return obj->seq;     
}    


/* Function:  replace_len_Sequence(obj,len)
 *
 * Descrip:    Replace member variable len
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [Sequence *]
 * Arg:        len [OWNER] New value of the variable [int]
 *
 * Return [SOFT ]  member variable len [boolean]
 *
 */
boolean replace_len_Sequence(Sequence * obj,int len) 
{
    if( obj == NULL)     {  
      warn("In replacement function len for object Sequence, got a NULL object");    
      return FALSE;  
      }  
    obj->len = len;  
    return TRUE; 
}    


/* Function:  access_len_Sequence(obj)
 *
 * Descrip:    Access member variable len
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [Sequence *]
 *
 * Return [SOFT ]  member variable len [int]
 *
 */
int access_len_Sequence(Sequence * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function len for object Sequence, got a NULL object");   
      return 0;  
      }  
    return obj->len;     
}    


/* Function:  replace_maxlen_Sequence(obj,maxlen)
 *
 * Descrip:    Replace member variable maxlen
 *             For use principly by API functions
 *
 *
 * Arg:           obj [UNKN ] Object holding the variable [Sequence *]
 * Arg:        maxlen [OWNER] New value of the variable [int]
 *
 * Return [SOFT ]  member variable maxlen [boolean]
 *
 */
boolean replace_maxlen_Sequence(Sequence * obj,int maxlen) 
{
    if( obj == NULL)     {  
      warn("In replacement function maxlen for object Sequence, got a NULL object"); 
      return FALSE;  
      }  
    obj->maxlen = maxlen;    
    return TRUE; 
}    


/* Function:  access_maxlen_Sequence(obj)
 *
 * Descrip:    Access member variable maxlen
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [Sequence *]
 *
 * Return [SOFT ]  member variable maxlen [int]
 *
 */
int access_maxlen_Sequence(Sequence * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function maxlen for object Sequence, got a NULL object");    
      return 0;  
      }  
    return obj->maxlen;  
}    


/* Function:  replace_offset_Sequence(obj,offset)
 *
 * Descrip:    Replace member variable offset
 *             For use principly by API functions
 *
 *
 * Arg:           obj [UNKN ] Object holding the variable [Sequence *]
 * Arg:        offset [OWNER] New value of the variable [int]
 *
 * Return [SOFT ]  member variable offset [boolean]
 *
 */
boolean replace_offset_Sequence(Sequence * obj,int offset) 
{
    if( obj == NULL)     {  
      warn("In replacement function offset for object Sequence, got a NULL object"); 
      return FALSE;  
      }  
    obj->offset = offset;    
    return TRUE; 
}    


/* Function:  access_offset_Sequence(obj)
 *
 * Descrip:    Access member variable offset
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [Sequence *]
 *
 * Return [SOFT ]  member variable offset [int]
 *
 */
int access_offset_Sequence(Sequence * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function offset for object Sequence, got a NULL object");    
      return 0;  
      }  
    return obj->offset;  
}    


/* Function:  replace_end_Sequence(obj,end)
 *
 * Descrip:    Replace member variable end
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [Sequence *]
 * Arg:        end [OWNER] New value of the variable [int]
 *
 * Return [SOFT ]  member variable end [boolean]
 *
 */
boolean replace_end_Sequence(Sequence * obj,int end) 
{
    if( obj == NULL)     {  
      warn("In replacement function end for object Sequence, got a NULL object");    
      return FALSE;  
      }  
    obj->end = end;  
    return TRUE; 
}    


/* Function:  access_end_Sequence(obj)
 *
 * Descrip:    Access member variable end
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [Sequence *]
 *
 * Return [SOFT ]  member variable end [int]
 *
 */
int access_end_Sequence(Sequence * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function end for object Sequence, got a NULL object");   
      return 0;  
      }  
    return obj->end;     
}    


/* Function:  replace_type_Sequence(obj,type)
 *
 * Descrip:    Replace member variable type
 *             For use principly by API functions
 *
 *
 * Arg:         obj [UNKN ] Object holding the variable [Sequence *]
 * Arg:        type [OWNER] New value of the variable [char]
 *
 * Return [SOFT ]  member variable type [boolean]
 *
 */
boolean replace_type_Sequence(Sequence * obj,char type) 
{
    if( obj == NULL)     {  
      warn("In replacement function type for object Sequence, got a NULL object");   
      return FALSE;  
      }  
    obj->type = type;    
    return TRUE; 
}    


/* Function:  access_type_Sequence(obj)
 *
 * Descrip:    Access member variable type
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [Sequence *]
 *
 * Return [SOFT ]  member variable type [char]
 *
 */
char access_type_Sequence(Sequence * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function type for object Sequence, got a NULL object");  
      return 'u';    
      }  
    return obj->type;    
}    


/* Function:  replace_tax_id_Sequence(obj,tax_id)
 *
 * Descrip:    Replace member variable tax_id
 *             For use principly by API functions
 *
 *
 * Arg:           obj [UNKN ] Object holding the variable [Sequence *]
 * Arg:        tax_id [OWNER] New value of the variable [int]
 *
 * Return [SOFT ]  member variable tax_id [boolean]
 *
 */
boolean replace_tax_id_Sequence(Sequence * obj,int tax_id) 
{
    if( obj == NULL)     {  
      warn("In replacement function tax_id for object Sequence, got a NULL object"); 
      return FALSE;  
      }  
    obj->tax_id = tax_id;    
    return TRUE; 
}    


/* Function:  access_tax_id_Sequence(obj)
 *
 * Descrip:    Access member variable tax_id
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [Sequence *]
 *
 * Return [SOFT ]  member variable tax_id [int]
 *
 */
int access_tax_id_Sequence(Sequence * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function tax_id for object Sequence, got a NULL object");    
      return 0;  
      }  
    return obj->tax_id;  
}    


/* Function:  replace_weight_Sequence(obj,weight)
 *
 * Descrip:    Replace member variable weight
 *             For use principly by API functions
 *
 *
 * Arg:           obj [UNKN ] Object holding the variable [Sequence *]
 * Arg:        weight [OWNER] New value of the variable [double]
 *
 * Return [SOFT ]  member variable weight [boolean]
 *
 */
boolean replace_weight_Sequence(Sequence * obj,double weight) 
{
    if( obj == NULL)     {  
      warn("In replacement function weight for object Sequence, got a NULL object"); 
      return FALSE;  
      }  
    obj->weight = weight;    
    return TRUE; 
}    


/* Function:  access_weight_Sequence(obj)
 *
 * Descrip:    Access member variable weight
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [Sequence *]
 *
 * Return [SOFT ]  member variable weight [double]
 *
 */
double access_weight_Sequence(Sequence * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function weight for object Sequence, got a NULL object");    
      return 0;  
      }  
    return obj->weight;  
}    


/* Function:  replace_desc_Sequence(obj,desc)
 *
 * Descrip:    Replace member variable desc
 *             For use principly by API functions
 *
 *
 * Arg:         obj [UNKN ] Object holding the variable [Sequence *]
 * Arg:        desc [OWNER] New value of the variable [char *]
 *
 * Return [SOFT ]  member variable desc [boolean]
 *
 */
boolean replace_desc_Sequence(Sequence * obj,char * desc) 
{
    if( obj == NULL)     {  
      warn("In replacement function desc for object Sequence, got a NULL object");   
      return FALSE;  
      }  
    obj->desc = desc;    
    return TRUE; 
}    


/* Function:  access_desc_Sequence(obj)
 *
 * Descrip:    Access member variable desc
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [Sequence *]
 *
 * Return [SOFT ]  member variable desc [char *]
 *
 */
char * access_desc_Sequence(Sequence * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function desc for object Sequence, got a NULL object");  
      return NULL;   
      }  
    return obj->desc;    
}    


/* Function:  access_set_SequenceSet(obj,i)
 *
 * Descrip:    Access members stored in the set list
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the list [SequenceSet *]
 * Arg:          i [UNKN ] Position in the list [int]
 *
 * Return [SOFT ]  Element of the list [Sequence *]
 *
 */
Sequence * access_set_SequenceSet(SequenceSet * obj,int i) 
{
    if( obj == NULL)     {  
      warn("In accessor function set for object SequenceSet, got a NULL object");    
      return NULL;   
      }  
    if( obj->len <= i )  {  
      warn("In accessor function set for object SequenceSet, index %%d is greater than list length %%d",i,obj->len); 
      return NULL;   
      }  
    return obj->set[i];  
}    


/* Function:  length_set_SequenceSet(obj)
 *
 * Descrip:    discover the length of the list
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the list [SequenceSet *]
 *
 * Return [UNKN ]  length of the list [int]
 *
 */
int length_set_SequenceSet(SequenceSet * obj) 
{
    if( obj == NULL)     {  
      warn("In length function set for object SequenceSet, got a NULL object");  
      return -1;     
      }  
    return obj->len;     
}    



#ifdef _cplusplus
}
#endif
