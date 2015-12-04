#ifdef _cplusplus
extern "C" {
#endif
#include "compmat.h"


/* Function:  simple_CompProb(match,rnd)
 *
 * Descrip:    Makes a simple CompProb matrix 
 *
 *
 * Arg:        match [UNKN ] Undocumented argument [Probability]
 * Arg:          rnd [UNKN ] Undocumented argument [Probability]
 *
 * Return [UNKN ]  Undocumented return value [CompProb *]
 *
 */
# line 50 "compmat.dy"
CompProb * simple_CompProb(Probability match,Probability rnd)
{
  CompProb * p;
  Probability rem;
  int i,j;

  p = CompProb_alloc();
  rem = (1.0-match)/20;

  for(i=0;i<26;i++) {
    for(j=i;j<26;j++) {
      if( i == j ) {
	p->comp[i][j] = match/rnd;
      } else {
	p->comp[i][j] = rem/rnd;
	p->comp[j][i] = rem/rnd;
      }
    }
  }

  return p;
}

/* Function:  fold_column_RandomModel_CompProb(cp,rm)
 *
 * Descrip:    Folds a random model in over the columns
 *
 *
 * Arg:        cp [UNKN ] Undocumented argument [CompProb *]
 * Arg:        rm [UNKN ] Undocumented argument [RandomModel *]
 *
 */
# line 76 "compmat.dy"
void fold_column_RandomModel_CompProb(CompProb * cp,RandomModel * rm)
{
  int i,j;

  for(i=0;i<26;i++) {
    for(j=0;j<26;j++) {
      if( rm->aminoacid[j] < 0.00000000000001 ) {
        continue;
      }
      cp->comp[i][j] = cp->comp[i][j] / rm->aminoacid[j];
    }
  }
}


/* Function:  simple_aa_CompProb(match,set,rnd)
 *
 * Descrip:    Makes a simple CompProb with simple aa rules
 *
 *
 * Arg:        match [UNKN ] Undocumented argument [Probability]
 * Arg:          set [UNKN ] Undocumented argument [Probability]
 * Arg:          rnd [UNKN ] Undocumented argument [Probability]
 *
 * Return [UNKN ]  Undocumented return value [CompProb *]
 *
 */
# line 94 "compmat.dy"
CompProb * simple_aa_CompProb(Probability match,Probability set,Probability rnd)
{
  CompProb * p;
  Probability rem;
  int i,j,k;
  char * aa_set[] = { "FYW","LIV","RK","DE" };

  p = CompProb_alloc();
  rem = (1.0-match)/20;

  for(i=0;i<26;i++) {
    for(j=i;j<26;j++) {
      if( i == j ) {
	p->comp[i][j] = match/rnd;
      } else {
	for( k = 0;k<4;k++) {
	  if( strchr(aa_set[k],'A'+i) != NULL && strchr(aa_set[k],'A'+j) != NULL ) {
	    p->comp[i][j] = set/rnd;
	    p->comp[j][i] = set/rnd;
	    fprintf(stderr,"Setting %c,%c to %.2f\n",'A'+i,'A'+j,set/rnd);
	  } else {
	    p->comp[i][j] = rem/rnd;
	    p->comp[j][i] = rem/rnd;
	  }
	}
      }
    }
  }

  return p;
}

/* Function:  CompMat_from_CompProb(cp)
 *
 * Descrip:    Maps a CompProb to a CompMat going through
 *             Probability2Score
 *
 *
 * Arg:        cp [UNKN ] Undocumented argument [CompProb *]
 *
 * Return [UNKN ]  Undocumented return value [CompMat *]
 *
 */
# line 130 "compmat.dy"
CompMat * CompMat_from_CompProb(CompProb * cp)
{
  int i,j;
  CompMat * cm;

  cm = CompMat_alloc();

  for(i=0;i<26;i++) {
    for(j=0;j<26;j++) {
      cm->comp[i][j] = Probability2Score(cp->comp[i][j]);
    }
  }


  return cm;
}

/* Function:  CompProb_from_halfbit(cm)
 *
 * Descrip:    Maps a halfbit matrix to a prob matrix by rebasing
 *             etc. 
 *
 *             *Really* not sensible!
 *
 *
 * Arg:        cm [UNKN ] Undocumented argument [CompMat *]
 *
 * Return [UNKN ]  Undocumented return value [CompProb *]
 *
 */
# line 153 "compmat.dy"
CompProb * CompProb_from_halfbit(CompMat * cm)
{
  CompProb * out;
  int i,j;

  out = CompProb_alloc();
  for(i=0;i<26;i++) 
    for(j=0;j<26;j++)
      out->comp[i][j] = halfbit2Probability(cm->comp[i][j]);

  return out;
}

/* Function:  CompMat_from_halfbit(cm)
 *
 * Descrip:    flips a halfbit based matrix (eg, blosum62) into a score
 *             based matrix just by rebasing the log etc. 
 *
 *             Not a sensible function ...
 *
 *
 *
 * Arg:        cm [UNKN ] Undocumented argument [CompMat *]
 *
 * Return [UNKN ]  Undocumented return value [CompMat *]
 *
 */
# line 173 "compmat.dy"
CompMat * CompMat_from_halfbit(CompMat * cm)
{
  CompMat * out;
  int i,j;

  out = CompMat_alloc();
  for(i=0;i<26;i++) 
    for(j=0;j<26;j++)
      out->comp[i][j] = Probability2Score(halfbit2Probability(cm->comp[i][j]));

  return out;
}


/* Function:  factor_CompMat(cm,factor)
 *
 * Descrip:    multiples all the scores by the amount
 *
 *
 * Arg:            cm [UNKN ] compmat object [CompMat *]
 * Arg:        factor [UNKN ] amount to multiple by [int]
 *
 */
# line 193 "compmat.dy"
void factor_CompMat(CompMat * cm,int factor)
{
  int i,j;

  for(i=0;i<26;i++)
    for(j=0;j<26;j++)
      cm->comp[i][j] *= factor;

}


/* Function:  fail_safe_CompMat_access(cm,aa1,aa2)
 *
 * Descrip:    gives the fail form of the macro CompMat_AAMATCH which 
 *             checks that aa1 and a2 are sensible and that cm is not NULL.
 *
 *
 * Arg:         cm [UNKN ] compmat object [CompMat *]
 * Arg:        aa1 [UNKN ] first amino acid [int]
 * Arg:        aa2 [UNKN ] second amino acid [int]
 *
 * Return [UNKN ]  Undocumented return value [Score]
 *
 */
# line 212 "compmat.dy"
Score fail_safe_CompMat_access(CompMat * cm,int aa1,int aa2)
{
  if( cm == NULL) {
    warn("Attempting to access NULL CompMat.");
    return 0;
  }

  if( aa1 < 0 || aa1 >= 26 || aa2 < 0 || aa2 > 26) {
    warn("Attempting to access CompMat with aminoacid numbers %d:%d (they must be bound between 0:26, returning 0 here",aa1,aa2);
    return 0;
  }

  else return cm->comp[aa1][aa2];
}

/* Function:  write_Blast_CompMat(cm,ofp)
 *
 * Descrip:    writes a protien CompMat with a standard
 *             alphabet.
 *
 *
 * Arg:         cm [UNKN ] CompMat object [CompMat *]
 * Arg:        ofp [UNKN ] file to output [FILE *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 234 "compmat.dy"
boolean write_Blast_CompMat(CompMat * cm,FILE * ofp)
{
  return write_Blast_CompMat_alphabet(cm,"ARNDCQEGHILKMFPSTWYVBZX",ofp);
}

/* Function:  bad_CompMat_alphabet(al)
 *
 * Descrip:    checks that this string is ok for BLAST alphabet mappings.
 *
 *
 * Arg:        al [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 243 "compmat.dy"
boolean bad_CompMat_alphabet(char * al)
{
  char * runner;
  for(runner=al;*runner;runner++)
    if( !isalpha((int)*runner) && toupper((int)*runner) != *runner )  {
      warn("Attempting to use [%s] as a CompMat alphabet: needs all upper case, no spaced letters",al);
      return TRUE;
    }

  return FALSE;
}
    
/* Function:  write_Blast_CompMat_alphabet(cm,alphabet,ofp)
 *
 * Descrip:    actualy writes out the Blast CFormat. The alphabet is 
 *             what order you want the amino acids. If you want the 
 *             standard format use /write_Blast_CompMat
 *
 *
 * Arg:              cm [UNKN ] comp mat object [CompMat *]
 * Arg:        alphabet [UNKN ] string for alphabet to be used [char *]
 * Arg:             ofp [UNKN ] fileoutput [FILE *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 264 "compmat.dy"
boolean write_Blast_CompMat_alphabet(CompMat * cm,char * alphabet,FILE * ofp)
{
  char * runner;
  char * run2;
  int minnumbers[26];
  int len;
  int i;
  int min;

  if( bad_CompMat_alphabet(alphabet) == TRUE )
    return FALSE; /* warning already issued */

  fprintf(ofp,"# File made by *Wise CompMat library module\n");
  fprintf(ofp,"# Comparison matrix in BLAST format\n");
  fprintf(ofp,"# Usually matrices are given in half-bit units\n");
  fprintf(ofp,"# First line is alphabet, the * column is the lowest score\n");
  fprintf(ofp,"#      File Created [%s]\n",now_string());
  fprintf(ofp,"#      Matrix name  [%s]\n",cm->name == NULL ? "No Name" : cm->name);


  /*** print out alphabet wit correct spacing ***/
  fprintf(ofp," %c",alphabet[0]);
  for(runner=alphabet+1;*runner;runner++)
    fprintf(ofp,"  %c",*runner);
  fprintf(ofp,"  *\n");

  /*** print out each row: remember the minimum number for printing later **/

  for(runner=alphabet,len=0;*runner;runner++) {
    min = cm->comp[*runner-'A'][0];
    fprintf(ofp,"%- d",cm->comp[*runner-'A'][0]);
    for(run2=alphabet+1;*run2;run2++) {
      fprintf(ofp," %- d",cm->comp[*runner-'A'][*run2-'A']);
      if( cm->comp[*runner-'A'][*run2-'A'] < min )
	min = cm->comp[*runner-'A'][*run2-'A'];
    }

    minnumbers[len++] = min;
    fprintf(ofp," % d\n",min);
  }

  /*** final row... *  ***/

  fprintf(ofp,"% d",minnumbers[0]);
  for(i=1;i<len;i++)
    fprintf(ofp," % d",minnumbers[i]);

  fprintf(ofp,"  1\n");

  /*** finished! ***/

  return TRUE;
}

/* Function:  read_Blast_file_CompMat(filename)
 *
 * Descrip:    Opens file, reads matrix, closes file.
 *             calls /read_Blast_CompMat for the actual format
 *             reading. Uses /openfile to open the file,
 *             so will open from config files.
 *
 *
 * Arg:        filename [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [CompMat *]
 *
 */
# line 324 "compmat.dy"
CompMat * read_Blast_file_CompMat(char * filename)
{
  CompMat * out;
  FILE * ifp;

  ifp = openfile(filename,"r");

  if( ifp == NULL ) {
    warn("Could not open %s as a filename for read Blast matrix",filename);
    return NULL;
  }

  out = read_Blast_CompMat(ifp);
  if( out != NULL ) {
    out->name = stringalloc(filename);
  }

  fclose(ifp);

  return out;
}

/* Function:  read_Blast_CompMat(ifp)
 *
 * Descrip:    reads a BLAST format matrix and
 *             allocates a new ComMat structure.
 *
 *
 * Arg:        ifp [UNKN ] Undocumented argument [FILE *]
 *
 * Return [UNKN ]  Undocumented return value [CompMat *]
 *
 */
# line 350 "compmat.dy"
CompMat * read_Blast_CompMat(FILE * ifp)
{
  char buffer[MAXLINE];
  int alphabet[MAXLINE];
  char * runner;
  int len;
  int linenum;
  int row;
  CompMat * out;



  /*** 
    Skip over # lines...

    read first line: is alphabet ie
    A R T G .... *
    ***/

  while( fgets(buffer,MAXLINE,ifp) != NULL)
    if( buffer[0] != '#')
      break;

  /** loop over line, getting letters: warn if longer than one, or not a letter **/

  for(len=0,runner=strtok(buffer,spacestr);runner != NULL;runner=strtok(NULL,spacestr)) {
    if( *runner == '*' )
      break; /* end column */

    if( !isalpha((int)*runner) || strlen(runner) > 1 ) {
      warn("In read Blast comp mat, probably an error: trying to interpret [%s] as an amino acid",runner);
      return NULL;
    }
    
    alphabet[len++] = toupper((int)*runner) -'A';
  }


  out = blank_CompMat();
  linenum = 0;

  /** get len lines, each line, get len numbers and put them away **/

  while( fgets(buffer,MAXLINE,ifp) != NULL ) {

    if( linenum >= len )
      break;
    
    for(runner=strtok(buffer,spacestr),row = 0;runner != NULL && row < len;runner=strtok(NULL,spacestr),row++) {
      if( is_integer_string(runner,&out->comp[alphabet[linenum]][alphabet[row]]) == FALSE ) {
	warn("In read Blast comp mat, probably an error: trying to interpret [%s] as a number ... continuing",runner);
      }
    }
    linenum++;
  }

  return out;

}

/* Function:  read_Blast_file_CompProb(filename)
 *
 * Descrip:    Reads a BLAST format comp prob from file
 *
 *
 * Arg:        filename [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [CompProb *]
 *
 */
# line 413 "compmat.dy"
CompProb * read_Blast_file_CompProb(char * filename)
{
  FILE * ifp;
  CompProb * out;

  ifp = openfile(filename,"r");

  if( ifp == NULL ) {
    warn("Could not open %s for compprob",filename);
    return NULL;
  }

  out = read_Blast_CompProb(ifp);

  fclose(ifp);

  return out;
}


/* Function:  read_Blast_CompProb(ifp)
 *
 * Descrip:    reads a BLAST format matrix and
 *             allocates a new CompProb structure.
 *
 *
 * Arg:        ifp [READ ] file input [FILE *]
 *
 * Return [UNKN ]  newly allocated CompProb [CompProb *]
 *
 */
# line 440 "compmat.dy"
CompProb * read_Blast_CompProb(FILE * ifp)
{
  char buffer[MAXLINE];
  int alphabet[MAXLINE];
  char * runner;
  int len;
  int linenum;
  int row;
  CompProb * out;



  /*** 
    Skip over # lines...

    read first line: is alphabet ie
    A R T G .... *
    ***/

  while( fgets(buffer,MAXLINE,ifp) != NULL)
    if( buffer[0] != '#')
      break;

  /** loop over line, getting letters: warn if longer than one, or not a letter **/

  for(len=0,runner=strtok(buffer,spacestr);runner != NULL;runner=strtok(NULL,spacestr)) {
    if( *runner == '*' )
      break; /* end column */

    if( !isalpha((int)*runner) || strlen(runner) > 1 ) {
      warn("In read Blast comp mat, probably an error: trying to interpret [%s] as an amino acid",runner);
      return NULL;
    }
    
    alphabet[len++] = toupper((int)*runner) -'A';
  }


  out = blank_CompProb();
  linenum = 0;

  /** get len lines, each line, get len numbers and put them away **/

  while( fgets(buffer,MAXLINE,ifp) != NULL ) {
    if( strstartcmp(buffer,"//") == 0 )
      break;
    if( linenum >= len )
      break;
    
    for(runner=strtok(buffer,spacestr),row = 0;runner != NULL && row < len;runner=strtok(NULL,spacestr),row++) {
      if( is_double_string(runner,&out->comp[alphabet[linenum]][alphabet[row]]) == FALSE ) {
	warn("In read Blast comp prob, probably an error: trying to interpret [%s] as a number ... continuing",runner);
      }
    }
    linenum++;
  }

  return out;
}

/* Function:  blank_CompMat(void)
 *
 * Descrip:    makes a 0,0 matrix
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [CompMat *]
 *
 */
# line 503 "compmat.dy"
CompMat * blank_CompMat(void)
{
  register int i;
  register int j;
  CompMat * out;


  out = CompMat_alloc();
  if( out == NULL)
    return NULL;

  for(i=0;i<26;i++)
    for(j=0;j<26;j++)
      out->comp[i][j] = 0;

  return out;
}

/* Function:  blank_CompProb(void)
 *
 * Descrip:    makes a 1.0 prob matrix
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [CompProb *]
 *
 */
# line 524 "compmat.dy"
CompProb * blank_CompProb(void)
{
  register int i;
  register int j;
  CompProb * out;


  out = CompProb_alloc();
  if( out == NULL)
    return NULL;

  for(i=0;i<26;i++)
    for(j=0;j<26;j++)
      out->comp[i][j] = 1.0;

  return out;
}






# line 594 "compmat.c"
/* Function:  hard_link_CompProb(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [CompProb *]
 *
 * Return [UNKN ]  Undocumented return value [CompProb *]
 *
 */
CompProb * hard_link_CompProb(CompProb * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a CompProb object: passed a NULL object");    
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  CompProb_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [CompProb *]
 *
 */
CompProb * CompProb_alloc(void) 
{
    CompProb * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(CompProb *) ckalloc (sizeof(CompProb))) == NULL)    {  
      warn("CompProb_alloc failed ");    
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    /* comp[26][26] is an array: no default possible */ 
    out->name = NULL;    


    return out;  
}    


/* Function:  free_CompProb(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [CompProb *]
 *
 * Return [UNKN ]  Undocumented return value [CompProb *]
 *
 */
CompProb * free_CompProb(CompProb * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a CompProb obj. Should be trappable");  
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


/* Function:  hard_link_CompMat(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [CompMat *]
 *
 * Return [UNKN ]  Undocumented return value [CompMat *]
 *
 */
CompMat * hard_link_CompMat(CompMat * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a CompMat object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  CompMat_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [CompMat *]
 *
 */
CompMat * CompMat_alloc(void) 
{
    CompMat * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(CompMat *) ckalloc (sizeof(CompMat))) == NULL)  {  
      warn("CompMat_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    /* comp[26][26] is an array: no default possible */ 
    out->name = NULL;    


    return out;  
}    


/* Function:  free_CompMat(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [CompMat *]
 *
 * Return [UNKN ]  Undocumented return value [CompMat *]
 *
 */
CompMat * free_CompMat(CompMat * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a CompMat obj. Should be trappable");   
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


/* Function:  replace_name_CompMat(obj,name)
 *
 * Descrip:    Replace member variable name
 *             For use principly by API functions
 *
 *
 * Arg:         obj [UNKN ] Object holding the variable [CompMat *]
 * Arg:        name [OWNER] New value of the variable [char *]
 *
 * Return [SOFT ]  member variable name [boolean]
 *
 */
boolean replace_name_CompMat(CompMat * obj,char * name) 
{
    if( obj == NULL)     {  
      warn("In replacement function name for object CompMat, got a NULL object");    
      return FALSE;  
      }  
    obj->name = name;    
    return TRUE; 
}    


/* Function:  access_name_CompMat(obj)
 *
 * Descrip:    Access member variable name
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [CompMat *]
 *
 * Return [SOFT ]  member variable name [char *]
 *
 */
char * access_name_CompMat(CompMat * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function name for object CompMat, got a NULL object");   
      return NULL;   
      }  
    return obj->name;    
}    



#ifdef _cplusplus
}
#endif
