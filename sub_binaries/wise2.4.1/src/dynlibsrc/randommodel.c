#ifdef _cplusplus
extern "C" {
#endif
#include "randommodel.h"

/* Function:  draw_random_aa_RandomModel(rm)
 *
 * Descrip:    Draws an amino acid from the random distribution
 *
 *
 * Arg:        rm [UNKN ] Undocumented argument [RandomModel *]
 *
 * Return [UNKN ]  Undocumented return value [char]
 *
 */
# line 56 "randommodel.dy"
char draw_random_aa_RandomModel(RandomModel * rm)
{
  double draw,tot;
  int i;

  draw = random_0_to_1();
  for(tot=rm->aminoacid[0],i=0;draw > tot && i<26;tot += rm->aminoacid[++i]) 
    ;
  if( i >= 26 ) {
    warn("Weird - got a draw %f which outside of random model total %f\n",draw,tot);
    return '?';
  }
  return 'A'+i;
}

/* Function:  draw_random_base_RandomModelDNA(rm)
 *
 * Descrip:    Draws a base from the random distribution
 *
 *
 * Arg:        rm [UNKN ] Undocumented argument [RandomModelDNA *]
 *
 * Return [UNKN ]  Undocumented return value [char]
 *
 */
# line 74 "randommodel.dy"
char draw_random_base_RandomModelDNA(RandomModelDNA * rm)
{
  double draw,tot;
  int i;

  draw = random_0_to_1();
  for(tot=rm->base[0],i=0;draw > tot && i<4;tot += rm->base[++i]) 
    ;
  if( i >= 26 ) {
    warn("Weird - got a draw %f which outside of random model total %f\n",draw,tot);
    return '?';
  }
  return char_from_base(i);
}


/* Function:  RandomCodonScore_from_RandomCodon(rc)
 *
 * Descrip:    Makes a score RandomCodon (log space)
 *             from a probability based random codon
 *
 *
 * Arg:        rc [UNKN ] Undocumented argument [RandomCodon *]
 *
 * Return [UNKN ]  Undocumented return value [RandomCodonScore *]
 *
 */
# line 94 "randommodel.dy"
RandomCodonScore * RandomCodonScore_from_RandomCodon(RandomCodon * rc)
{
  RandomCodonScore  * out;

  out = RandomCodonScore_alloc();

  Probability2Score_move(rc->codon,out->codon,126);

  if( rc-> name != NULL)
    out->name = stringalloc(rc->name);

  return out;
}


/* Function:  flatten_RandomCodon(rc)
 *
 * Descrip:    Sets all probabilities to 1.0 - ie,
 *             odds them to themselves.
 *
 *             This is equivalent to saying that the randomcodon
 *             is being odd-ratioed to itself
 *
 *             Also equivalent of saying all the scores (in log
 *             space) will be 0
 *
 *
 * Arg:        rc [UNKN ] Undocumented argument [RandomCodon *]
 *
 */
# line 119 "randommodel.dy"
void flatten_RandomCodon(RandomCodon * rc)
{
  int i;

  for(i=0;i<125;i++)
    rc->codon[i] = 1.0;

}

/* Function:  fold_in_RandomModelDNA_into_RandomCodon(rc,rmd)
 *
 * Descrip:    Makes the randomcodon numbers become the odds ratio
 *             between their probabilitys and flat dna random model
 *             (0th order markov)
 *
 *
 * Arg:         rc [UNKN ] Undocumented argument [RandomCodon *]
 * Arg:        rmd [UNKN ] Undocumented argument [RandomModelDNA *]
 *
 */
# line 133 "randommodel.dy"
void fold_in_RandomModelDNA_into_RandomCodon(RandomCodon * rc,RandomModelDNA * rmd)
{
  register int one;
  register int two;
  register int three;

  if( rc == NULL || rmd == NULL ) {
    warn("Passed in NULL objects to fold_in_RandomModelDNA_into_RandomCodon");
  }

  for(one=0;one < 5;one++)
    for(two =0;two<5;two ++)
      for(three=0;three<5;three++)
	rc->codon[(one*25)+(two*5)+(three)] /= (rmd->base[one]*rmd->base[two]*rmd->base[three]);


}

/* Function:  show_RandomCodonScore(rcs,ofp)
 *
 * Descrip:    shows RandomCodonScore
 *
 *             for debugging
 *
 *
 * Arg:        rcs [UNKN ] Undocumented argument [RandomCodonScore *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 156 "randommodel.dy"
void show_RandomCodonScore(RandomCodonScore * rcs,FILE * ofp)
{
  register int i;

  for(i=0;i<125;i++) {
    fprintf(ofp,"Score %3d is %d\n",i,rcs->codon[i]);

  }
}

/* Function:  show_RandomModelDNAScore(rds,ofp)
 *
 * Descrip:    shows RandomModelsDNAScore
 *
 *             for debugging
 *
 *
 * Arg:        rds [UNKN ] Undocumented argument [RandomModelDNAScore *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 171 "randommodel.dy"
void show_RandomModelDNAScore(RandomModelDNAScore * rds,FILE * ofp)
{
  register int i;

  for(i=0;i<5;i++) {
    fprintf(ofp,"Base %d[%c], Score %d [prob %.2f]\n",i,char_from_base(i),rds->base[i],Score2Probability(rds->base[i]));
  }
}

/* Function:  folded_RandomModelDNAScore_from_2RMD(dis,rnd)
 *
 * Descrip:    gives a odds ratio between two random models
 *
 *
 * Arg:        dis [UNKN ] Undocumented argument [RandomModelDNA *]
 * Arg:        rnd [UNKN ] Undocumented argument [RandomModelDNA *]
 *
 * Return [UNKN ]  Undocumented return value [RandomModelDNAScore *]
 *
 */
# line 183 "randommodel.dy"
RandomModelDNAScore * folded_RandomModelDNAScore_from_2RMD(RandomModelDNA * dis,RandomModelDNA * rnd)
{
  int i;
  RandomModelDNAScore * out;

  out = RandomModelDNAScore_alloc();

  for(i=0;i<5;i++)
    out->base[i]= Probability2Score(dis->base[i]/rnd->base[i]);

  return out;
}



/* Function:  RandomCodon_from_raw_CodonFrequency(codon[64],*ct)
 *
 * Descrip:    From raw counts (no adjustment to amino acids) of codons
 *             gives you a RandomCodon model
 *
 *             No prior is used (? perhaps should have a flat prior)
 *
 *             N's are handled correctly
 *
 *
 * Arg:        codon[64] [UNKN ] Undocumented argument [double]
 * Arg:              *ct [UNKN ] Undocumented argument [CodonTable]
 *
 * Return [UNKN ]  Undocumented return value [RandomCodon  *]
 *
 */
# line 206 "randommodel.dy"
RandomCodon  * RandomCodon_from_raw_CodonFrequency(double codon[64],CodonTable *ct)
{

  RandomCodon * out;
  register int i;
  double total = 0.0;
  base one;
  base two;
  base three;
  int o,t,r;



  /** codon frequencies here *do* include protein amino acid freq...
    ie, they are raw frequencies, not adjusted for a codon table **/


  /** the only thing is that we have to figure out how to
    deal with N'd codons, which will just be summed over... **/


  out= RandomCodon_alloc();


  for(i=0;i<64;i++) {
    total += codon[i];
  }

  for(i=0;i<125;i++) {

    if( has_random_bases(i) == FALSE ) {

      out->codon[i] = codon[base4_codon_from_codon(i)]/total;
    }
    
    else {
      all_bases_from_codon(i,&one,&two,&three);
      
      if( one == BASE_N && two != BASE_N && three != BASE_N ) {
	for(o=0;o<4;o++)
	  out->codon[i] += (codon[o*16+two*4+three]/total);
      }
      else if( one == BASE_N && two == BASE_N && three != BASE_N) {
	for(o=0;o<4;o++)
	  for(t=0;t<4;t++)
	    out->codon[i] += (codon[o*16+t*4+three]/total);
      }
      else if( one == BASE_N && two == BASE_N && three == BASE_N) {
	for(o=0;o<4;o++)
	  for(t=0;t<4;t++)
	    for(r=0;r<4;r++)
	      out->codon[i] += (codon[o*16+t*4+r]/total);
      }
      else if( one != BASE_N && two == BASE_N && three != BASE_N) {
	for(t=0;t<4;t++)
	  out->codon[i] += (codon[one*16+t*4+three]/total);
      }
      else if( one != BASE_N && two == BASE_N && three == BASE_N) {
	for(t=0;t<4;t++)
	  for(r=0;r<4;r++)
	    out->codon[i] += (codon[one*16+t*4+r]/total);
      }
      else if( one != BASE_N && two != BASE_N && three == BASE_N) {
	for(r=0;r<4;r++)
	  out->codon[i] += (codon[one*16+two*4+r]/total);
      }
    }
  }

  out->codon[125] = 0.0;


  return out;
} 


/* Function:  flat_RandomCodon(ct)
 *
 * Descrip:    Makes a flat RandomCodon from CodonTable
 *
 *
 * Arg:        ct [UNKN ] Undocumented argument [CodonTable *]
 *
 * Return [UNKN ]  Undocumented return value [RandomCodon *]
 *
 */
# line 285 "randommodel.dy"
RandomCodon * flat_RandomCodon(CodonTable * ct)
{
  int i;

  RandomCodon * rc;

  rc = RandomCodon_alloc();

  for(i=0;i<125;i++) {
    if( is_stop_codon(i,ct) ) {
      rc->codon[i] = 0.0;
    } else {
      rc->codon[i] = 1.0 / 61.0;
    }
  }

  return rc;
}
  

/* Function:  RandomModelDNAScore_from_RandomModelDNA(rmd)
 *
 * Descrip:    Gives you a log space RandomModelDNAScore
 *             from a probability space one
 *
 *
 * Arg:        rmd [UNKN ] Undocumented argument [RandomModelDNA *]
 *
 * Return [UNKN ]  Undocumented return value [RandomModelDNAScore *]
 *
 */
# line 309 "randommodel.dy"
RandomModelDNAScore * RandomModelDNAScore_from_RandomModelDNA(RandomModelDNA * rmd)
{
  RandomModelDNAScore * out;
  register int i;

  out = RandomModelDNAScore_alloc();
  if( out == NULL )
    return NULL;

  for(i=0;i<5;i++) {
    out->base[i] = Probability2Score(rmd->base[i]);
  }


  return out;
}

/* Function:  RandomModelDNA_std(void)
 *
 * Descrip:    Returns a structure with 0.25 set in each place
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [RandomModelDNA *]
 *
 */
# line 329 "randommodel.dy"
RandomModelDNA * RandomModelDNA_std(void)
{
  register int i;
  RandomModelDNA * out;

  out = RandomModelDNA_alloc();
  if( out == NULL )
    return NULL;

  for(i=0;i<4;i++)
    out->base[i] = (1.0) / 4.0;

  out->base[4] = 1.0;

  return out;
}

/* Function:  RandomModelDNA_std_human(void)
 *
 * Descrip:    Set human random model (slightly G/C)
 *
 *             Not sure where I got the numbers now. Ooops
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [RandomModelDNA *]
 *
 */
# line 351 "randommodel.dy"
RandomModelDNA * RandomModelDNA_std_human(void)
{
  RandomModelDNA * out;

  out = RandomModelDNA_alloc();
  if( out == NULL )
    return NULL;

  out->base[BASE_A]= 0.245;
  out->base[BASE_G]= 0.251;
  out->base[BASE_C]= 0.253;
  out->base[BASE_T]= 0.248;

  out->base[BASE_N]= 1.0;

  return out;
}

/* Function:  Score_Sequence_is_random(s,rms)
 *
 * Descrip:    Gives the score of a Sequence vs a random model
 *
 *
 * Arg:          s [UNKN ] Undocumented argument [Sequence *]
 * Arg:        rms [UNKN ] Undocumented argument [RandomModelScoreaa *]
 *
 * Return [UNKN ]  Undocumented return value [Score]
 *
 */
# line 372 "randommodel.dy"
Score     Score_Sequence_is_random(Sequence * s,RandomModelScoreaa * rms)
{
  register int i;
  Score sc = 0;

  for(i=0;i<s->len;i++)
    sc += rms->aminoacid[s->seq[i]-'A'];

  return sc;
}

/* Function:  Prob_Sequence_is_random(s,rm)
 *
 * Descrip:    Gives the probability of a Sequence vs a random model
 *
 *
 * Arg:         s [UNKN ] Undocumented argument [Sequence *]
 * Arg:        rm [UNKN ] Undocumented argument [RandomModel *]
 *
 * Return [UNKN ]  Undocumented return value [Probability]
 *
 */
# line 386 "randommodel.dy"
Probability Prob_Sequence_is_random(Sequence * s,RandomModel * rm)
{
  register int i;
  Probability p = 1.0;
  

  for(i=0;i<s->len;i++) {
    p *= rm->aminoacid[s->seq[i]-'A'];
  }

  return p;
}
  

/* Function:  RandomModelScoreaa_from_RandomModel(rm)
 *
 * Descrip:    Gives a score based RandomModel from a probability based one
 *
 *
 * Arg:        rm [UNKN ] Undocumented argument [RandomModel *]
 *
 * Return [UNKN ]  Undocumented return value [RandomModelScoreaa *]
 *
 */
# line 403 "randommodel.dy"
RandomModelScoreaa * RandomModelScoreaa_from_RandomModel(RandomModel * rm)
{
  register int i;
  RandomModelScoreaa * out;

  out = RandomModelScoreaa_alloc();
  if( out == NULL )
    return NULL;

  for(i=0;i<26;i++) 
    out->aminoacid[i] = Probability2Score(rm->aminoacid[i]);

  return out;
}

/* Function:  default_RandomModel(void)
 *
 * Descrip:    Gives a default random model numbers from
 *             swissprot34- via the HMMEr1 package
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [RandomModel *]
 *
 */
# line 422 "randommodel.dy"
RandomModel * default_RandomModel(void)
{
  RandomModel * out;
  int i;

  out = RandomModel_alloc();

  for(i=0;i<26;i++)
    out->aminoacid[i] = 0.00000001;
  
  out->aminoacid['A' -'A'] = 0.08713;
  out->aminoacid['C' -'A'] = 0.03347;
  out->aminoacid['D' -'A'] = 0.04687;
  out->aminoacid['E' -'A'] = 0.04953;
  out->aminoacid['F' -'A'] = 0.03977;
  out->aminoacid['G' -'A'] = 0.08861;
  out->aminoacid['H' -'A'] = 0.03362;
  out->aminoacid['I' -'A'] = 0.03689;
  out->aminoacid['K' -'A'] = 0.08048;
  out->aminoacid['L' -'A'] = 0.08536;
  out->aminoacid['M' -'A'] = 0.01475;
  out->aminoacid['N' -'A'] = 0.04043;
  out->aminoacid['P' -'A'] = 0.05068;
  out->aminoacid['Q' -'A'] = 0.03826;
  out->aminoacid['R' -'A'] = 0.04090;
  out->aminoacid['S' -'A'] = 0.06958;
  out->aminoacid['T' -'A'] = 0.05854;
  out->aminoacid['V' -'A'] = 0.06472;
  out->aminoacid['W' -'A'] = 0.01049;
  out->aminoacid['Y' -'A'] = 0.02992;

  return out;
}

/* Function:  read_RandomModel(ifp)
 *
 * Descrip:    Reads a simplistic RandomModel file of
 *
 *             C 0.0123
 *
 *             etc type of format
 *
 *
 * Arg:        ifp [UNKN ] Undocumented argument [FILE *]
 *
 * Return [UNKN ]  Undocumented return value [RandomModel *]
 *
 */
# line 463 "randommodel.dy"
RandomModel * read_RandomModel(FILE * ifp)
{
  char buffer[MAXLINE];
  char c;
  float f;
  RandomModel * out;
  register int i;

  out = RandomModel_alloc();

  if( out == NULL ) 
    return NULL;
  
  for(i=0;i<26;i++)
    out->aminoacid[i] = 0.000001;

  while( fgets(buffer,MAXLINE,ifp) != NULL ) {
    if( buffer[0] == '!' )
      continue;
    sscanf(buffer,"%c %f",&c,&f);
    c = toupper((int)c);

    if( c-'A' < 0 || c-'A' > 26 ) {
      warn("Have picked up an awfully dodgy character [%c] in reading random model",c);
    }

    out->aminoacid[c-'A'] = f;
  }

  return out;
}

# line 551 "randommodel.c"
/* Function:  hard_link_RandomModel(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [RandomModel *]
 *
 * Return [UNKN ]  Undocumented return value [RandomModel *]
 *
 */
RandomModel * hard_link_RandomModel(RandomModel * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a RandomModel object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  RandomModel_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [RandomModel *]
 *
 */
RandomModel * RandomModel_alloc(void) 
{
    RandomModel * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(RandomModel *) ckalloc (sizeof(RandomModel))) == NULL)  {  
      warn("RandomModel_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    /* aminoacid[26] is an array: no default possible */ 
    out->name = NULL;    


    return out;  
}    


/* Function:  free_RandomModel(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [RandomModel *]
 *
 * Return [UNKN ]  Undocumented return value [RandomModel *]
 *
 */
RandomModel * free_RandomModel(RandomModel * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a RandomModel obj. Should be trappable");   
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


/* Function:  hard_link_RandomModelScoreaa(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [RandomModelScoreaa *]
 *
 * Return [UNKN ]  Undocumented return value [RandomModelScoreaa *]
 *
 */
RandomModelScoreaa * hard_link_RandomModelScoreaa(RandomModelScoreaa * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a RandomModelScoreaa object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  RandomModelScoreaa_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [RandomModelScoreaa *]
 *
 */
RandomModelScoreaa * RandomModelScoreaa_alloc(void) 
{
    RandomModelScoreaa * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(RandomModelScoreaa *) ckalloc (sizeof(RandomModelScoreaa))) == NULL)    {  
      warn("RandomModelScoreaa_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    /* aminoacid[26] is an array: no default possible */ 
    out->name = NULL;    


    return out;  
}    


/* Function:  free_RandomModelScoreaa(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [RandomModelScoreaa *]
 *
 * Return [UNKN ]  Undocumented return value [RandomModelScoreaa *]
 *
 */
RandomModelScoreaa * free_RandomModelScoreaa(RandomModelScoreaa * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a RandomModelScoreaa obj. Should be trappable");    
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


/* Function:  hard_link_RandomCodonScore(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [RandomCodonScore *]
 *
 * Return [UNKN ]  Undocumented return value [RandomCodonScore *]
 *
 */
RandomCodonScore * hard_link_RandomCodonScore(RandomCodonScore * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a RandomCodonScore object: passed a NULL object");    
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  RandomCodonScore_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [RandomCodonScore *]
 *
 */
RandomCodonScore * RandomCodonScore_alloc(void) 
{
    RandomCodonScore * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(RandomCodonScore *) ckalloc (sizeof(RandomCodonScore))) == NULL)    {  
      warn("RandomCodonScore_alloc failed ");    
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    /* codon[126] is an array: no default possible */ 
    out->name = NULL;    


    return out;  
}    


/* Function:  free_RandomCodonScore(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [RandomCodonScore *]
 *
 * Return [UNKN ]  Undocumented return value [RandomCodonScore *]
 *
 */
RandomCodonScore * free_RandomCodonScore(RandomCodonScore * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a RandomCodonScore obj. Should be trappable");  
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


/* Function:  hard_link_RandomCodon(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [RandomCodon *]
 *
 * Return [UNKN ]  Undocumented return value [RandomCodon *]
 *
 */
RandomCodon * hard_link_RandomCodon(RandomCodon * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a RandomCodon object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  RandomCodon_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [RandomCodon *]
 *
 */
RandomCodon * RandomCodon_alloc(void) 
{
    RandomCodon * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(RandomCodon *) ckalloc (sizeof(RandomCodon))) == NULL)  {  
      warn("RandomCodon_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    /* codon[126] is an array: no default possible */ 
    out->name = NULL;    


    return out;  
}    


/* Function:  free_RandomCodon(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [RandomCodon *]
 *
 * Return [UNKN ]  Undocumented return value [RandomCodon *]
 *
 */
RandomCodon * free_RandomCodon(RandomCodon * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a RandomCodon obj. Should be trappable");   
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


/* Function:  hard_link_RandomModelDNA(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [RandomModelDNA *]
 *
 * Return [UNKN ]  Undocumented return value [RandomModelDNA *]
 *
 */
RandomModelDNA * hard_link_RandomModelDNA(RandomModelDNA * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a RandomModelDNA object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  RandomModelDNA_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [RandomModelDNA *]
 *
 */
RandomModelDNA * RandomModelDNA_alloc(void) 
{
    RandomModelDNA * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(RandomModelDNA *) ckalloc (sizeof(RandomModelDNA))) == NULL)    {  
      warn("RandomModelDNA_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    /* base[5] is an array: no default possible */ 
    out->name = NULL;    


    return out;  
}    


/* Function:  free_RandomModelDNA(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [RandomModelDNA *]
 *
 * Return [UNKN ]  Undocumented return value [RandomModelDNA *]
 *
 */
RandomModelDNA * free_RandomModelDNA(RandomModelDNA * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a RandomModelDNA obj. Should be trappable");    
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


/* Function:  hard_link_RandomModelDNAScore(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [RandomModelDNAScore *]
 *
 * Return [UNKN ]  Undocumented return value [RandomModelDNAScore *]
 *
 */
RandomModelDNAScore * hard_link_RandomModelDNAScore(RandomModelDNAScore * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a RandomModelDNAScore object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  RandomModelDNAScore_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [RandomModelDNAScore *]
 *
 */
RandomModelDNAScore * RandomModelDNAScore_alloc(void) 
{
    RandomModelDNAScore * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(RandomModelDNAScore *) ckalloc (sizeof(RandomModelDNAScore))) == NULL)  {  
      warn("RandomModelDNAScore_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    /* base[5] is an array: no default possible */ 
    out->name = NULL;    


    return out;  
}    


/* Function:  free_RandomModelDNAScore(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [RandomModelDNAScore *]
 *
 * Return [UNKN ]  Undocumented return value [RandomModelDNAScore *]
 *
 */
RandomModelDNAScore * free_RandomModelDNAScore(RandomModelDNAScore * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a RandomModelDNAScore obj. Should be trappable");   
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


/* Function:  replace_name_RandomModelDNA(obj,name)
 *
 * Descrip:    Replace member variable name
 *             For use principly by API functions
 *
 *
 * Arg:         obj [UNKN ] Object holding the variable [RandomModelDNA *]
 * Arg:        name [OWNER] New value of the variable [char *]
 *
 * Return [SOFT ]  member variable name [boolean]
 *
 */
boolean replace_name_RandomModelDNA(RandomModelDNA * obj,char * name) 
{
    if( obj == NULL)     {  
      warn("In replacement function name for object RandomModelDNA, got a NULL object"); 
      return FALSE;  
      }  
    obj->name = name;    
    return TRUE; 
}    


/* Function:  access_name_RandomModelDNA(obj)
 *
 * Descrip:    Access member variable name
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [RandomModelDNA *]
 *
 * Return [SOFT ]  member variable name [char *]
 *
 */
char * access_name_RandomModelDNA(RandomModelDNA * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function name for object RandomModelDNA, got a NULL object");    
      return NULL;   
      }  
    return obj->name;    
}    


/* Function:  replace_name_RandomModel(obj,name)
 *
 * Descrip:    Replace member variable name
 *             For use principly by API functions
 *
 *
 * Arg:         obj [UNKN ] Object holding the variable [RandomModel *]
 * Arg:        name [OWNER] New value of the variable [char *]
 *
 * Return [SOFT ]  member variable name [boolean]
 *
 */
boolean replace_name_RandomModel(RandomModel * obj,char * name) 
{
    if( obj == NULL)     {  
      warn("In replacement function name for object RandomModel, got a NULL object");    
      return FALSE;  
      }  
    obj->name = name;    
    return TRUE; 
}    


/* Function:  access_name_RandomModel(obj)
 *
 * Descrip:    Access member variable name
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [RandomModel *]
 *
 * Return [SOFT ]  member variable name [char *]
 *
 */
char * access_name_RandomModel(RandomModel * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function name for object RandomModel, got a NULL object");   
      return NULL;   
      }  
    return obj->name;    
}    



#ifdef _cplusplus
}
#endif
