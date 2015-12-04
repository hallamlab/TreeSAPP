#ifdef _cplusplus
extern "C" {
#endif
#include "codonmapper.h"


/* Function:  CodonFrequence_from_raw_counts(codon,ct)
 *
 * Descrip:    Builds a codon frequency from raw counts as just an array
 *
 *
 * Arg:        codon [UNKN ] Undocumented argument [double *]
 * Arg:           ct [UNKN ] Undocumented argument [CodonTable *]
 *
 * Return [UNKN ]  Undocumented return value [CodonFrequency *]
 *
 */
# line 56 "codonmapper.dy"
CodonFrequency * CodonFrequence_from_raw_counts(double * codon,CodonTable * ct)
{
  double total[26];
  CodonFrequency * cf;
  register int i;
  int c;
  register int j;

  for(i=0;i<26;i++) {
    total[i] = 0.0;
    for(j=0;j<64;j++) {
      c = codon_from_base4_codon(j);
      if( ct->codon_str[c] == ('A' + i) ) {
	total[i] += codon[j];
      }
    }
  }

  cf = CodonFrequency_alloc();

  for(i=0;i<64;i++) {
    c = codon_from_base4_codon(i);

    if( is_stop_codon(c,ct)  )
      continue;

    
    if( codon[i] < 0.0000000001)
      cf->freq[i] = 0.0;
    else {
      if( total[ct->codon_str[c] -'A'] < 0.00000001 ) {
	warn("For codon %d, amino acid %c, we have no frequency",i,ct->codon_str[i]);
      }
      else cf->freq[i] = codon[i] / total[ct->codon_str[c]-'A'];
    }
  }

  return cf;
}



/* Function:  show_CodonMapper(cm,ofp)
 *
 * Descrip:    Shows codon mapper in vaguely human form
 *
 *
 * Arg:         cm [UNKN ] Undocumented argument [CodonMapper *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 102 "codonmapper.dy"
void show_CodonMapper(CodonMapper * cm,FILE * ofp)
{
  register int i;
  register int j;

  for(i=0;i<125;i++) { 
    fprintf(ofp,"[%3d][%c] %.2f",i,aminoacid_from_codon(cm->ct,i),cm->codon_map[i][0]);
    for(j=1;j<26;j++) 
      fprintf(ofp,",%.2f",cm->codon_map[i][j]);
    fprintf(ofp,"\n");
  }

}

/* Function:  show_CodonFrequency(cf,ct,ofp)
 *
 * Descrip:    Shows codon frequency in vaguely human form
 *
 *
 * Arg:         cf [UNKN ] Undocumented argument [CodonFrequency *]
 * Arg:         ct [UNKN ] Undocumented argument [CodonTable *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 120 "codonmapper.dy"
void show_CodonFrequency(CodonFrequency * cf,CodonTable * ct,FILE * ofp)
{
  int i;

  for(i=0;i<64;i++)
    fprintf(ofp,"[%3d][%c] %.2f\n",i,aminoacid_from_codon(ct,codon_from_base4_codon(i)),cf->freq[i]);

}


/* Function:  flat_CodonMapper(ct)
 *
 * Descrip:    Makes a CodonMapper with no codon bias
 *             or error possiblities from codon table
 *
 *
 *
 * Arg:        ct [UNKN ] Codon Table giving codon->aa info [CodonTable *]
 *
 * Return [UNKN ]  Undocumented return value [CodonMapper *]
 *
 */
# line 137 "codonmapper.dy"
CodonMapper * flat_CodonMapper(CodonTable * ct)
{
  CodonFrequency * cf;
  CodonMapper * out;

  cf = flat_CodonFrequency(ct);

  out = new_CodonMapper(ct,cf);

  free_CodonFrequency(cf);

  return out;
}


/* Function:  flat_CodonFrequency(ct)
 *
 * Descrip:    Makes a no-biased codon Frequency.
 *             Probabaly most used in /flat_CodonMapper
 *
 *
 *
 * Arg:        ct [UNKN ] Undocumented argument [CodonTable *]
 *
 * Return [UNKN ]  Undocumented return value [CodonFrequency *]
 *
 */
# line 158 "codonmapper.dy"
CodonFrequency * flat_CodonFrequency(CodonTable * ct)
{
  int number[26];
  CodonFrequency * out;
  register int i;


  out = CodonFrequency_alloc();

  for(i=0;i<26;i++)
    number[i] = 0;

  for(i=0;i<64;i++) 
    out->freq[i]= 0.0;

  for(i=0;i<125;i++)
    if( has_random_bases(i) == FALSE && is_stop_codon(i,ct) == FALSE) 
      number[aminoacid_no_from_codon(ct,i)]++;

  for(i=0;i<64;i++) {
    if( is_stop_codon(codon_from_base4_codon(i),ct) == FALSE && number[aminoacid_no_from_codon(ct,codon_from_base4_codon(i))] != 0) 
      out->freq[i] = 1.0 / number[aminoacid_no_from_codon(ct,codon_from_base4_codon(i))];
  }

  return out;
}
  
/***

***/

/* Function:  construct_amino_number_array(number,ct)
 *
 * Descrip:    Assummes number is an int * of length 26
 *
 *             Files up each position with the number of codons representing that aa
 *
 *
 *
 * Arg:        number [UNKN ] Undocumented argument [int *]
 * Arg:            ct [UNKN ] Undocumented argument [CodonTable *]
 *
 */
# line 197 "codonmapper.dy"
void construct_amino_number_array(int * number,CodonTable * ct)
{
  register int i;
  register int j;

  for(i=0;i<26;i++) {
    number[i] = 0;
    for(j=0;j<64;j++) {
      if( ct->codon_str[j] == 'A' + i)
	number[i]++;
    }
  }

}
  

/* Function:  map_codon_array_CodonMapper(codon_array,protein_array,stop,cm)
 *
 * Descrip:    Now defunct.
 *
 *
 * Arg:          codon_array [UNKN ] Undocumented argument [double *]
 * Arg:        protein_array [UNKN ] Undocumented argument [double *]
 * Arg:                 stop [UNKN ] Undocumented argument [double]
 * Arg:                   cm [UNKN ] Undocumented argument [CodonMapper *]
 *
 */
# line 218 "codonmapper.dy"
void map_codon_array_CodonMapper(double * codon_array,double * protein_array,double stop,CodonMapper * cm)
{
  register int i;

  for(i=0;i<125;i++) {
    if( is_stop_codon(i,cm->ct)== TRUE ) 
      codon_array[i] = stop;
    else codon_array[i] = map_codon_CodonMapper(i,protein_array,cm);
  }

}

/* Function:  true_map_codon_array_CodonMapper(codon_array,protein_array,cm)
 *
 * Descrip:    Takes an array of probabilities from 0-26 in protein array
 *             and writes into codon_array the adjusted probability from the
 *             codon mapper. Ie, maps a protein emission line to a codon emission
 *             line. This is the main use of CodonMapper.
 *
 *
 *
 * Arg:          codon_array [WRITE] array (0-124) for the codon probabilities to be placed [double *]
 * Arg:        protein_array [READ ] array (0-25) for the protein probabilities to be read from [const double *]
 * Arg:                   cm [UNKN ] Codon Mapper that provides the protein->codon mapping [CodonMapper *]
 *
 */
# line 241 "codonmapper.dy"
void true_map_codon_array_CodonMapper(double * codon_array,const double * protein_array,CodonMapper * cm)
{
  int i;

  for(i=0;i<125;i++) {
    codon_array[i] = map_codon_CodonMapper(i,protein_array,cm);
  }
}

/* Function:  map_codon_CodonMapper(codon,protein_array,cm)
 *
 * Descrip:    Does the mapping for a single codon from a protein_array
 *
 *
 * Arg:                codon [UNKN ] Undocumented argument [int]
 * Arg:        protein_array [UNKN ] Undocumented argument [const double *]
 * Arg:                   cm [UNKN ] Undocumented argument [CodonMapper *]
 *
 * Return [UNKN ]  Undocumented return value [double]
 *
 */
# line 255 "codonmapper.dy"
double map_codon_CodonMapper(int codon,const double * protein_array,CodonMapper * cm)
{
  register int i;
  double out = 0.0;

 
  if( cm->codon_map[codon][0] < 0.0 ) {
    warn("Attempting to map a codon with below zero prob in map_codon_CodonMapper. This is bad news....");
    return 0.0;
  }
  
  for(i=0;i<26;i++) {
    out += protein_array[i] * cm->codon_map[codon][i];
  }

  return out;
}

/* Function:  sprinkle_errors_over_CodonMapper(cm,error)
 *
 * Descrip:    Takes a codon mapper and assummes that the majority of errors
 *             are due to a single base change in the codon at probability error.
 *             Therefore, for each codon it adds error * prob(codon) * 0.25 to each 
 *             other codon one base away, taking away therefore the result.
 *
 *
 *
 * Arg:           cm [READ ] CodonMapper to be sprinkled [CodonMapper *]
 * Arg:        error [UNKN ] substitution error rate [double]
 *
 */
# line 283 "codonmapper.dy"
void sprinkle_errors_over_CodonMapper(CodonMapper * cm,double error)
{
  int i;
  int j;
  int k;
  base one;
  base two;
  base three;
  int new_codon;
  double scratch[125][26];


  /*
   * put all the codons into scratch, but 
   * subtracting 3*error which is the possibility
   * of an error.
   *
   * The self->self errors (eg, G to G) will be
   * added back later
   */

  for(i=0;i<125;i++) 
    for(j=0;j<26;j++) 
      scratch[i][j] = cm->codon_map[i][j] * (1-(3*error));


  /*
   * Now for each codon, find the single base change,
   * and add the probability for each amino acid onto it
   * factored by 0.25
   *
   */

  for(i=0;i<125;i++) {
    all_bases_from_codon(i,&one,&two,&three);

    if( one != BASE_N) {
      for(j=0;j<4;j++) {
	new_codon = j*25 + two * 5 + three;

	for(k=0;k<26;k++) 
	  scratch[i][k] += (cm->codon_map[new_codon][k] * error * 0.25);
      }
    }
    if( two != BASE_N) {
      for(j=0;j<4;j++) {
	
	new_codon = one*25 + j * 5 + three;

	for(k=0;k<26;k++) 
	  scratch[i][k] += (cm->codon_map[new_codon][k] * error * 0.25);
      }
    }
    if( three != BASE_N) {
      for(j=0;j<4;j++) {
	
	new_codon = one*25 + two * 5 + j;

	for(k=0;k<26;k++) 
	  scratch[i][k] += (cm->codon_map[new_codon][k] * error * 0.25);
      }
    }  
    
  }


  /*
   * Now map back to original memory
   */

  for(i=0;i<125;i++) 
    for(j=0;j<26;j++) 
      cm->codon_map[i][j] = scratch[i][j];



}
	


/* Function:  new_CodonMapper(ct,cf)
 *
 * Descrip:    The only way you should make a CodonMapper!
 *
 *             Makes a codon mapper from CodonTable and frequency
 *
 *
 * Arg:        ct [UNKN ] Undocumented argument [CodonTable *]
 * Arg:        cf [UNKN ] Undocumented argument [CodonFrequency *]
 *
 * Return [UNKN ]  Undocumented return value [CodonMapper *]
 *
 */
# line 369 "codonmapper.dy"
CodonMapper * new_CodonMapper(CodonTable * ct,CodonFrequency * cf)
{
  register int i;
  register int j;
  int k;
  base one;
  base two;
  base three;
  int base4;
  int oi,ti,ri;
  double total_freq;
  CodonMapper * out;
  

  out = CodonMapper_alloc();



  out->ct = hard_link_CodonTable(ct);

  for(i=0;i<125;i++) {
    for(j=0;j<26;j++) 
      out->codon_map[i][j] =0.0;


    if( has_random_bases(i) == FALSE ) {
      if( is_stop_codon(i,ct) == TRUE ) {
	for(k=0;k<26;k++)
	out->codon_map[i][k] = (0.0);
      }
      else {
	out->codon_map[i][aminoacid_no_from_codon(ct,i)] = cf->freq[base4_codon_from_codon(i)];
      }
    }
    else { /*** is a random base ***/


      /***
	sneaky stuff. What we want to do is loop over all possible
	codons, adding up their frequencies for the amino acids 
	they represent. This is done by looping over all possible
	bases for each position and then letting through ones 
	which either have an N at this position or is the actual base.

	***/



      all_bases_from_codon(i,&one,&two,&three);

      total_freq = 0.0;

      for(oi=0;oi<4;oi++)
	for(ti=0;ti<4;ti++) 
	  for(ri=0;ri<4;ri++) {
	    if( (one == BASE_N || one == oi) && (two == BASE_N || two == ti) && (three == BASE_N || three == ri) ) {

	      base4 = codon_from_base4_codon(oi*16+ti*4+ri);
	      if( !is_stop_codon(base4,ct) ) {
		out->codon_map[i][aminoacid_no_from_codon(ct,base4)] += cf->freq[base4_codon_from_codon(base4)];
	      }
	    } /* end of if one == BASE_N ||  one == oi */
	  }  /* end of for oi,ti,ri */

      
    } /* end of else */
  }

  return out;
}


# line 443 "codonmapper.c"
/* Function:  hard_link_CodonMapper(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [CodonMapper *]
 *
 * Return [UNKN ]  Undocumented return value [CodonMapper *]
 *
 */
CodonMapper * hard_link_CodonMapper(CodonMapper * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a CodonMapper object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  CodonMapper_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [CodonMapper *]
 *
 */
CodonMapper * CodonMapper_alloc(void) 
{
    CodonMapper * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(CodonMapper *) ckalloc (sizeof(CodonMapper))) == NULL)  {  
      warn("CodonMapper_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->ct = NULL;  
    /* codon_map[125][26] is an array: no default possible */ 


    return out;  
}    


/* Function:  free_CodonMapper(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [CodonMapper *]
 *
 * Return [UNKN ]  Undocumented return value [CodonMapper *]
 *
 */
CodonMapper * free_CodonMapper(CodonMapper * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a CodonMapper obj. Should be trappable");   
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
    if( obj->ct != NULL) 
      free_CodonTable(obj->ct);  


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_CodonFrequency(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [CodonFrequency *]
 *
 * Return [UNKN ]  Undocumented return value [CodonFrequency *]
 *
 */
CodonFrequency * hard_link_CodonFrequency(CodonFrequency * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a CodonFrequency object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  CodonFrequency_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [CodonFrequency *]
 *
 */
CodonFrequency * CodonFrequency_alloc(void) 
{
    CodonFrequency * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(CodonFrequency *) ckalloc (sizeof(CodonFrequency))) == NULL)    {  
      warn("CodonFrequency_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    /* freq[64] is an array: no default possible */ 


    return out;  
}    


/* Function:  free_CodonFrequency(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [CodonFrequency *]
 *
 * Return [UNKN ]  Undocumented return value [CodonFrequency *]
 *
 */
CodonFrequency * free_CodonFrequency(CodonFrequency * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a CodonFrequency obj. Should be trappable");    
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


    ckfree(obj); 
    return NULL; 
}    


/* Function:  replace_ct_CodonMapper(obj,ct)
 *
 * Descrip:    Replace member variable ct
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [CodonMapper *]
 * Arg:         ct [OWNER] New value of the variable [CodonTable *]
 *
 * Return [SOFT ]  member variable ct [boolean]
 *
 */
boolean replace_ct_CodonMapper(CodonMapper * obj,CodonTable * ct) 
{
    if( obj == NULL)     {  
      warn("In replacement function ct for object CodonMapper, got a NULL object");  
      return FALSE;  
      }  
    obj->ct = ct;    
    return TRUE; 
}    


/* Function:  access_ct_CodonMapper(obj)
 *
 * Descrip:    Access member variable ct
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [CodonMapper *]
 *
 * Return [SOFT ]  member variable ct [CodonTable *]
 *
 */
CodonTable * access_ct_CodonMapper(CodonMapper * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function ct for object CodonMapper, got a NULL object"); 
      return NULL;   
      }  
    return obj->ct;  
}    



#ifdef _cplusplus
}
#endif
