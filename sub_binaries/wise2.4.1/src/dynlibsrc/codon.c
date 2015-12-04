#ifdef _cplusplus
extern "C" {
#endif
#include "codon.h"


/* Function:  read_CodonTable_file(file)
 *
 * Descrip:    Opens filename, reads it as if a Ewan style
 *             codon table and closes.
 *
 *
 * Arg:        file [READ ] filename to open [char *]
 *
 * Return [OWNER]  A codon-table, NULL if error [CodonTable *]
 *
 */
# line 80 "codon.dy"
CodonTable * read_CodonTable_file(char * file)
{
  FILE * ifp;
  CodonTable * out;

  ifp = openfile(file,"r");

  if( ifp == NULL) {
    warn("Could not open file %s as codon table file",file);
    return NULL;
  }

  out = read_CodonTable(ifp);

  fclose(ifp);

  return out;
}
  
/* Function:  read_CodonTable(ifp)
 *
 * Descrip:    reads a codon table from a filestream in Ewan
 *             format.
 *
 *             As Ewan format is really bad and has no start/stop
 *             this will effectively read to the end of the file.
 *             Ooops.
 *
 *
 * Arg:        ifp [READ ] file input [FILE *]
 *
 * Return [UNKN ]  Undocumented return value [CodonTable *]
 *
 */
# line 109 "codon.dy"
CodonTable * read_CodonTable(FILE * ifp)
{
  char buffer[MAXLINE];
  CodonTable * out;
  codon c;
  char * runner;
  char * run2;

  out = CodonTable_alloc();
  memset(out->codon_str,'x',125);

  while( fgets(buffer,MAXLINE,ifp) != NULL) {
    if( buffer[0] == '#' || buffer[0] == '!')
      continue;
    runner = strtok(buffer,spacestr);
    run2 = strtok(NULL,spacestr);
    
    if( runner == NULL || run2 == NULL ){
      warn("Unable to read a line in codon table");
    }

    c = codon_from_seq(runner);

    out->codon_str[c] = *run2;
  }

  return out;
}
  

/* Function:  alloc_aminoacid_from_seq(ct,seq)
 *
 * Descrip:    Not very useful function: allocates a single amino
 *             acid (ie, buffer length one) so it can be freed later.
 *
 *
 * Arg:         ct [READ ] codon table [CodonTable *]
 * Arg:        seq [READ ] pointer to DNA Sequence chars [char *]
 *
 * Return [UNKN ]  Undocumented return value [char *]
 *
 */
# line 147 "codon.dy"
char * alloc_aminoacid_from_seq(CodonTable * ct,char * seq)
{
  char buf[2];

  buf[1] = '\0';

  buf[0] = aminoacid_from_codon(ct,codon_from_seq(seq));

  return stringalloc(buf);
}

/* Function:  aminoacid_from_seq(ct,seq)
 *
 * Descrip:    Returns the amino acid for this position in the DNA sequence
 *             Takes the pointer +1 and +2 points.
 *
 *             No error checks implemented. Probably a mistake ;)
 *
 *
 * Arg:         ct [READ ] codon table [CodonTable *]
 * Arg:        seq [READ ] pointer to DNA chars [char *]
 *
 * Return [UNKN ]  an amino acid char (A-Z) [aa]
 *
 */
# line 168 "codon.dy"
aa aminoacid_from_seq(CodonTable * ct,char * seq)
{
  return aminoacid_from_codon(ct,codon_from_seq(seq));
}

/* Function:  aminoacid_from_codon(ct,c)
 *
 * Descrip:    returns amino acid for this codon number (NB codon numbers 0-125)
 *
 *
 * Arg:        ct [READ ] codon table [CodonTable *]
 * Arg:         c [READ ] codon number [codon]
 *
 * Return [READ ]  aminoacid that is this codon (X for ambiguous, * for stop) [aa]
 *
 */
# line 180 "codon.dy"
aa aminoacid_from_codon(CodonTable * ct,codon c)
{
  return ct->codon_str[c];
}

/* Function:  aminoacid_no_from_codon(ct,c)
 *
 * Descrip:    a sister function to aminoacid_from_codon:
 *             returns amino acid number (0-26) for this codon number (0-125)
 *
 *
 * Arg:        ct [READ ] codon table [CodonTable *]
 * Arg:         c [READ ] codon number [codon]
 *
 * Return [READ ]  aminoacid number [0-26] for this codon [int]
 *
 */
# line 193 "codon.dy"
int aminoacid_no_from_codon(CodonTable * ct,codon c)
{
  return (ct->codon_str[c] - 'A');
}

/* Function:  is_stop_codon(c,ct)
 *
 * Descrip:    tells you whether this codon number is really a stop
 *             in this translation table
 *
 *
 * Arg:         c [READ ] codon number [codon]
 * Arg:        ct [READ ] codon table [CodonTable *]
 *
 * Return [UNKN ]  TRUE if is stop, FALSE otherwise [boolean]
 *
 */
# line 206 "codon.dy"
boolean is_stop_codon(codon c,CodonTable * ct)
{
  aa a;

  a = aminoacid_from_codon(ct,c);

  if( a == 'X' || a == '*') {
    return TRUE;
  }

  return FALSE;
}

/* Function:  is_non_ambiguous_codon_seq(seq)
 *
 * Descrip:    Tells you if this codon is a real codon
 *
 *
 * Arg:        seq [READ ] pointer to DNA sequence [char *]
 *
 * Return [UNKN ]  TRUE if real codon, FALSE if contains N's [boolean]
 *
 */
# line 225 "codon.dy"
boolean is_non_ambiguous_codon_seq(char * seq)
{
  if( *seq == '\0' || *(seq+1) == '\0' || *(seq+2) == '\0') {
    warn("Attempting to find a codon number is something less than 3 bases long!");
    return FALSE;
  }

  if( base_from_char(*(seq++)) == BASE_N)
    return FALSE;
  if( base_from_char(*(seq++)) == BASE_N)
    return FALSE;
  if( base_from_char(*(seq)) == BASE_N)
    return FALSE;

  return TRUE;
}

/* Function:  is_valid_aminoacid(ct,c)
 *
 * Descrip:    Tells you if this letter (c) is recognised as a valid amino acid
 *             in this codon table
 *
 *
 * Arg:        ct [READ ] Codon Table [CodonTable *]
 * Arg:         c [UNKN ] aminoacid [char]
 *
 * Return [UNKN ]  TRUE if valid, FALSE if not. [boolean]
 *
 */
# line 250 "codon.dy"
boolean is_valid_aminoacid(CodonTable * ct,char c)
{
  if( strchr(ct->codon_str,c) != NULL )
    return TRUE;
  else return FALSE;
}

/* Function:  is_valid_base_char(c)
 *
 * Descrip:    Tells you if the letter is A,T,C,G,N (NB, N is ok).
 *
 *
 * Arg:        c [READ ] base [char]
 *
 * Return [UNKN ]  TRUE if (ATGCN) FALSE otherwise [boolean]
 *
 */
# line 263 "codon.dy"
boolean is_valid_base_char(char c)
{
  if( c == 'A' || c == 'T' || c == 'G' || c == 'C' || c == 'N')
    return TRUE;
  return FALSE;
}

/* Function:  codon_from_base4_codon(c)
 *
 * Descrip:    maps a 0-63 codon to a 0-123 codon. Suprisingly useful.
 *
 *
 * Arg:        c [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [codon]
 *
 */
# line 273 "codon.dy"
codon codon_from_base4_codon(int c)
{
  base one;
  base two;
  base three;

  one = c / 16;
  c -= one*16;

  two = c/4;
  c -= 4*two;
  three = c;

  return 25*one+5*two+three;
}

/* Function:  base4_codon_from_codon(c)
 *
 * Descrip:    maps a 0-125 codon to a 0-63 codon.
 *
 *             If ambiguous then returns 64 having issued a warning.
 *
 *
 * Arg:        c [READ ] codon 0-125 [codon]
 *
 * Return [UNKN ]  base 4 codon (0-63) [int]
 *
 */
# line 297 "codon.dy"
int base4_codon_from_codon(codon c)
{
  base one;
  base two;
  base three;

  all_bases_from_codon(c,&one,&two,&three);

  if( one == BASE_N || two == BASE_N || three == BASE_N)  {
    /* GSS 07:07:2000 DISABLED WARNING */
    warn("Attempting to convert an ambiguous codon to base 64"
         " returning 64");
    return 64;
  }

  return one*16 + two*4 + three;
}
  
/* Function:  has_random_bases(c)
 *
 * Descrip:    Tests to see if this codon number has any N's in it
 *
 *
 * Arg:        c [READ ] codon number 0-124 [codon]
 *
 * Return [UNKN ]  TRUE if has N's , FALSE otherwise [boolean]
 *
 */
# line 321 "codon.dy"
boolean has_random_bases(codon c)
{
  base o;
  base w;
  base t;

  o = base_from_codon(c,1);
  w = base_from_codon(c,2);
  t = base_from_codon(c,3);

  if( o == BASE_N || w == BASE_N || t == BASE_N)
    return TRUE;

  return FALSE;
}

/* Function:  permute_possible_random_bases(c,one,two,three)
 *
 * Descrip:    Bizarely useful function for calculating ambiguity scores.
 *
 *             This takes the codon c, and for each possible base, 
 *             if it is N, replaces it with one, two or three.
 *
 *             If the base is not N, it remains the same
 *
 *
 * Arg:            c [READ ] codon number [codon]
 * Arg:          one [READ ] base to replace first position if N [base]
 * Arg:          two [READ ] base to replace second position if N [base]
 * Arg:        three [READ ] base to replace third position if N [base]
 *
 * Return [UNKN ]  codon number  [codon]
 *
 */
# line 351 "codon.dy"
codon permute_possible_random_bases(codon c,base one,base two,base three)
{
  base o;
  base w;
  base t;

  o = base_from_codon(c,1);
  w = base_from_codon(c,2);
  t = base_from_codon(c,3);

  if( o == BASE_N )
    o = one;
  if( w == BASE_N )
    w = two;
  if( t == BASE_N )
    t = three;

  return one*25+two*5+three;
}


/* Function:  all_bases_from_codon(c,one,two,three)
 *
 * Descrip:    Really an internal function, by useful enough to
 *             encourage outside use.
 *
 *             Takes codon c and breaks it into 3 base-numbers
 *
 *
 * Arg:            c [UNKN ] Undocumented argument [codon]
 * Arg:          one [UNKN ] Undocumented argument [base *]
 * Arg:          two [UNKN ] Undocumented argument [base *]
 * Arg:        three [UNKN ] Undocumented argument [base *]
 *
 */
# line 378 "codon.dy"
void all_bases_from_codon(codon c,base * one,base * two,base * three)
{
  base o;
  base t;
  base r;

  o = c/25;
  c -= o*25;

  t = c/5;
  c -= t*5;

  r = c;

  if( one != NULL )
    *one = o;
  if( two != NULL )
    *two = t;
  if( three != NULL )
    *three = r;

}

/* Function:  reverse_codon(c)
 *
 * Descrip:    Reverses codon. Takes a forward codon number and
 *             builds the inverted codon number
 *
 *
 * Arg:        c [UNKN ] Undocumented argument [codon]
 *
 * Return [UNKN ]  Undocumented return value [codon]
 *
 */
# line 405 "codon.dy"
codon reverse_codon(codon c)
{
  base one;
  base two;
  base three;

  all_bases_from_codon(c,&one,&two,&three);

  one = complement_base(one);
  two = complement_base(two);
  three = complement_base(three);
  
  return three*125 + two * 25 + one;
}


/* Function:  base_from_codon(c,pos)
 *
 * Descrip:    Probably not the best function to use for this, but 
 *             useful. Takes a codon and with pos being 1,2,3 gives
 *             you the firt,second of third base
 *
 *
 * Arg:          c [UNKN ] Undocumented argument [codon]
 * Arg:        pos [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [base]
 *
 */
# line 426 "codon.dy"
base  base_from_codon(codon c,int pos)
{
  base one;
  base two;
  base three;

  one = c/25;
  c -= one*25;

  two = c/5;
  c -= two*5;

  three = c;

  switch(pos) {
  case 1 : return one;
  case 2: return two;
  case 3: return three;
  default :
    warn("This is not good news: got a position asked which is not 1,2,3 [%d](BTW - first base is 1)",pos);
    return BASE_N;
  }

}
 
/* Function:  codon_from_seq(seq)
 *
 * Descrip:    takes an ASCII coded pointer to a 3 base pair
 *             sequence (it could be the part of a sequence: it only
 *             assummes that the seq points with 3 chars at pos 0,1,2 
 *             in C coordinates from seq. No NULL is required). It 
 *             ives back the codon as made from standard mapping, ie,
 *             25*base_1+5*base_2 + base3 being a number from 0-124 inc.
 *
 *
 * Arg:        seq [UNKN ] pointer to sequence of at least 3 chrs long. [char *]
 *
 * Return [UNKN ]  Undocumented return value [codon]
 *
 */
# line 461 "codon.dy"
codon codon_from_seq(char * seq)
{
  base one;
  base two;
  base three;

  if( !is_base(*seq) || !is_base(*(seq+1)) || !is_base (*(seq+2)) ) {
    warn("Attempting to translate bad base %c%c%c",seq[0],seq[1],seq[2]);
    one = BASE_N;
    two = BASE_N;
    three = BASE_N;
  } else {
    one = base_from_char(*seq);
    two = base_from_char(*(seq+1));
    three = base_from_char(*(seq+2));
  }

  return one*25+two*5+three;
}

/* Function:  base4_codon_from_seq(seq)
 *
 * Descrip:    Sometimes it is more useful to work in base64, ie, 
 *             non N. this functions does the same thing as 
 *             /codon_from_seq but produces a seq being
 *             16*base1 + 4 *base2 + base3
 *
 *
 * Arg:        seq [UNKN ] pointer to sequence of at least 3 chrs long [char *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
# line 489 "codon.dy"
int base4_codon_from_seq(char * seq)
{
  base one;
  base two;
  base three;

  one = base_from_char(*seq);
  two = base_from_char(*(seq+1));
  three = base_from_char(*(seq+2));

  if( one == BASE_N || two == BASE_N || three == BASE_N)
    return 64;

  else return one*16+two*4+three;
}

/* Function:  char_from_base(b)
 *
 * Descrip:    maps a base number (-04 inc) to A,T,G,C,N
 *
 *
 * Arg:        b [UNKN ] Undocumented argument [base]
 *
 * Return [UNKN ]  Undocumented return value [char]
 *
 */
# line 508 "codon.dy"
char char_from_base(base b)
{
  switch(b) {
  case BASE_A : return 'A';
  case BASE_G : return 'G';
  case BASE_C : return 'C';
  case BASE_T : return 'T';
  case BASE_N : return 'N';
  default : return '?';
  }

}
  
/* Function:  base_from_char(c)
 *
 * Descrip:    maps a char (atcgn) to number, 
 *             case insensitive, returns BASE_N
 *             if not atcgn
 *
 *
 * Arg:        c [UNKN ] Undocumented argument [char]
 *
 * Return [UNKN ]  Undocumented return value [base]
 *
 */
# line 526 "codon.dy"
base base_from_char(char c)
{
  c = (char)toupper((int)c);

  switch(c) {
  case 'A' : return BASE_A;
  case 'T' : return BASE_T;
  case 'G' : return BASE_G;
  case 'C' : return BASE_C;
  case 'N' : return BASE_N;
  default :  return BASE_N;
  }
}

/* Function:  char_complement_base(c)
 *
 * Descrip:    the char equivalent of /complement_base.
 *             this gives the complement in char of a base
 *             in char. Does not check for non ATGCN
 *
 *
 * Arg:        c [UNKN ] Undocumented argument [char]
 *
 * Return [UNKN ]  Undocumented return value [char]
 *
 */
# line 545 "codon.dy"
char char_complement_base(char c)
{
  if( c == '-' || c == '.' || c == '~' ) {
    return c;
  }

  return char_from_base(complement_base(base_from_char(c)));
}

/* Function:  complement_base(b)
 *
 * Descrip:    gives back the complement as a number
 *             ofthe base (given as a number)
 *
 *
 * Arg:        b [UNKN ] Undocumented argument [base]
 *
 * Return [UNKN ]  Undocumented return value [base]
 *
 */
# line 558 "codon.dy"
base complement_base(base b)
{
  switch(b) {
  case BASE_A : return BASE_T;
  case BASE_G : return BASE_C;
  case BASE_C : return BASE_G;
  case BASE_T : return BASE_A;
  case BASE_N : return BASE_N;
  default : return BASE_N;
  }

}



/* Function:  four_fold_sites_CodonTable(*ct,seq)
 *
 * Descrip:    returns the number of four fold degenerate
 *             sites in this codon
 *
 *
 * Arg:        *ct [UNKN ] Undocumented argument [CodonTable]
 * Arg:        seq [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
# line 577 "codon.dy"
int four_fold_sites_CodonTable(CodonTable *ct,char * seq)
{
  int site_count;
  int i;
  int j;
  char codon[3];
  char aa;
  

  assert(ct);

  codon[0] = seq[0];
  codon[1] = seq[1];
  codon[2] = seq[2];

  aa = aminoacid_from_seq(ct,codon);
  
  site_count = 0;

  for(i=0;i<3;i++) {
    for( j=0;j<4;j++) {
      codon[i] = char_from_base(j);
      if( aa != aminoacid_from_seq(ct,codon) ) {
	break;
      }
    }
    if( j >= 4 ) {
      site_count++;
    }
    codon[i] = seq[i];
  }

  return site_count;

}


# line 652 "codon.c"
/* Function:  hard_link_CodonTable(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [CodonTable *]
 *
 * Return [UNKN ]  Undocumented return value [CodonTable *]
 *
 */
CodonTable * hard_link_CodonTable(CodonTable * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a CodonTable object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  CodonTable_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [CodonTable *]
 *
 */
CodonTable * CodonTable_alloc(void) 
{
    CodonTable * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(CodonTable *) ckalloc (sizeof(CodonTable))) == NULL)    {  
      warn("CodonTable_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    /* codon_str[125] is an array: no default possible */ 
    out->name = NULL;    


    return out;  
}    


/* Function:  free_CodonTable(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [CodonTable *]
 *
 * Return [UNKN ]  Undocumented return value [CodonTable *]
 *
 */
CodonTable * free_CodonTable(CodonTable * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a CodonTable obj. Should be trappable");    
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


    ckfree(obj); 
    return NULL; 
}    


/* Function:  replace_name_CodonTable(obj,name)
 *
 * Descrip:    Replace member variable name
 *             For use principly by API functions
 *
 *
 * Arg:         obj [UNKN ] Object holding the variable [CodonTable *]
 * Arg:        name [OWNER] New value of the variable [char *]
 *
 * Return [SOFT ]  member variable name [boolean]
 *
 */
boolean replace_name_CodonTable(CodonTable * obj,char * name) 
{
    if( obj == NULL)     {  
      warn("In replacement function name for object CodonTable, got a NULL object"); 
      return FALSE;  
      }  
    obj->name = name;    
    return TRUE; 
}    


/* Function:  access_name_CodonTable(obj)
 *
 * Descrip:    Access member variable name
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [CodonTable *]
 *
 * Return [SOFT ]  member variable name [char *]
 *
 */
char * access_name_CodonTable(CodonTable * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function name for object CodonTable, got a NULL object");    
      return NULL;   
      }  
    return obj->name;    
}    



#ifdef _cplusplus
}
#endif
