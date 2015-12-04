#ifdef _cplusplus
extern "C" {
#endif
#include "probability.h"


/* Function:  Probability_from_average_state_occupancy(length)
 *
 * Descrip:    for single state (exponetial decays) takes an average length
 *             and converts that to a probability that will produce that
 *             length (on average) for the state. NB... this *assumes* that
 *             you want a single state exp decay.
 *
 *
 * Arg:        length [UNKN ] average length of state [double]
 *
 * Return [UNKN ]  Undocumented return value [Probability]
 *
 */
# line 44 "probability.dy"
Probability Probability_from_average_state_occupancy(double length)
{
  return 1 - (1.0 / length);
}

/* Function:  state_occupancy_from_Probability(p)
 *
 * Descrip:    If you have a single state then this will tak
 *             the probability for the state->state transition and
 *             give you back the average length in the state
 *
 *
 * Arg:        p [UNKN ] probability of staying in the state [double]
 *
 * Return [UNKN ]  Undocumented return value [double]
 *
 */
# line 56 "probability.dy"
double      state_occupancy_from_Probability(double p)
{
  return 1 / (1-p);
}

/* Function:  Probability_logsum(one,two)
 *
 * Descrip:    gives back a score of the sum in
 *             probability space of the two scores.
 *
 *             This is the function verison of this
 *             code, which is not efficient *at all*
 *
 *
 * Arg:        one [UNKN ] Undocumented argument [Score]
 * Arg:        two [UNKN ] Undocumented argument [Score]
 *
 * Return [UNKN ]  Undocumented return value [Score]
 *
 */
# line 68 "probability.dy"
Score Probability_logsum(Score one,Score two)
{
  return Probability2Score(Score2Probability(one) + Score2Probability(two));
}

/* Function:  Probability2Bits(p)
 *
 * Descrip:    Gives back a Bits score for a probability
 *
 *
 * Arg:        p [UNKN ] Undocumented argument [Probability]
 *
 * Return [UNKN ]  Undocumented return value [Bits]
 *
 */
# line 76 "probability.dy"
Bits Probability2Bits(Probability p)
{
  return Score2Bits(Probability2Score(p));
}

/* Function:  show_Score_array(s,len,ofp)
 *
 * Descrip:    shows a score array as score, score ,score ...'
 *
 *
 * Arg:          s [UNKN ] Score array [Score *]
 * Arg:        len [UNKN ] length of array [int]
 * Arg:        ofp [UNKN ] output filestream [FILE *]
 *
 */
# line 88 "probability.dy"
void show_Score_array(Score * s,int len,FILE * ofp)
{
  register int i;

  fprintf(ofp,"\"%d",s[0]);
  for(i=1;i<len;i++) {
    fprintf(ofp,",%d",s[i]);
  }

  fprintf(ofp,"\"");

}
  
/* Function:  show_Probability_array_exp(p,len,ofp)
 *
 * Descrip:    shows a proability array in scientific notation.
 *
 *
 *
 * Arg:          p [UNKN ] probability array [Probability *]
 * Arg:        len [UNKN ] length of proability array [int]
 * Arg:        ofp [UNKN ] output filestream [FILE *]
 *
 */
# line 109 "probability.dy"
void show_Probability_array_exp(Probability * p,int len,FILE * ofp)
{
  int i;

  fprintf(ofp,"%4.4g",p[0]);
  for(i=1;i<len;i++) {
    fprintf(ofp,",%4.4g",p[i]);
  }

}

/* Function:  read_Probability_array(p,len,start_of_array)
 *
 * Descrip:    reads in a probability array of comma separated numbers.
 *             It calls /is_double_string to test whether the numbers are
 *             probabilities. It tries ito read in len numbers: if it runs out of
 *             commad separated guys it returns FALSE
 *
 *
 * Arg:                     p [UNKN ] Undocumented argument [Probability *]
 * Arg:                   len [UNKN ] Undocumented argument [int]
 * Arg:        start_of_array [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 126 "probability.dy"
boolean read_Probability_array(Probability * p,int len,char * start_of_array)
{
  char* runner;
  int no =0;
  
  for(runner=strtok(start_of_array,", ");runner != NULL;runner=strtok(NULL,", "),no++) {
    if( no >= len ) {
      return FALSE;
    }
    if( is_double_string(runner,&p[no]) == FALSE ) {
      return FALSE;
    }
  }
  return TRUE;
}

  
/* Function:  show_Probability_array(p,len,ofp)
 *
 * Descrip:    shows a proability array in %f notation.
 *
 *
 *
 * Arg:          p [UNKN ] probability array [Probability *]
 * Arg:        len [UNKN ] length of proability array [int]
 * Arg:        ofp [UNKN ] output filestream [FILE *]
 *
 */
# line 151 "probability.dy"
void show_Probability_array(Probability * p,int len,FILE * ofp)
{
  register int i;

  fprintf(ofp,"\"%f",p[0]);
  for(i=1;i<len;i++) {
    fprintf(ofp,",%f",p[i]);
  }

  fprintf(ofp,"\"");

}

/* Function:  sum_Probability_array(p,len)
 *
 * Descrip:    adds up the probability array given
 *
 *
 * Arg:          p [UNKN ] probability array  [Probability *]
 * Arg:        len [UNKN ] length of array [int]
 *
 * Return [UNKN ]  Undocumented return value [Probability]
 *
 */
# line 170 "probability.dy"
Probability sum_Probability_array(Probability * p,int len)
{
  register int i;
  Probability ret = 0.0;

  for(i=0;i<len;i++)
    ret += p[i];

  return ret;
}

/* Function:  set_Probability_array(set,p,len)
 *
 * Descrip:    Sets the probability array to p
 *
 *
 * Arg:        set [UNKN ] probability array to set [Probability *]
 * Arg:          p [UNKN ] probability to set it to  [Probability]
 * Arg:        len [UNKN ] length of probability array [int]
 *
 * Return [UNKN ]  Undocumented return value [Probability *]
 *
 */
# line 188 "probability.dy"
Probability * set_Probability_array(Probability * set,Probability p,int len)
{
  register int i;

  for(i=0;i<len;i++) 
    set[i] = p;

  return set;
}

/* Function:  Probability2Score_move(from,to,len)
 *
 * Descrip:    moves the probability array from to the (same length)
 *             score array to going through Probability2Score function
 *
 *
 * Arg:        from [UNKN ] probability array to get the numbers [Probability *]
 * Arg:          to [UNKN ] Score array to put the numbers [Score *]
 * Arg:         len [UNKN ] length of arrays [int]
 *
 * Return [UNKN ]  Undocumented return value [Score *]
 *
 */
# line 206 "probability.dy"
Score * Probability2Score_move(Probability * from,Score * to,int len)
{
  register int i;

  for(i=0;i<len;i++)
    to[i] = Probability2Score(from[i]);

  return to;
}

/* Function:  Probability_move(from,to,len)
 *
 * Descrip:    moves from to to 
 *
 *
 * Arg:        from [UNKN ] probability array with the numbers [const Probability *]
 * Arg:          to [UNKN ] probability array to be written into [Probability *]
 * Arg:         len [UNKN ] length [int]
 *
 * Return [UNKN ]  Undocumented return value [Probability *]
 *
 */
# line 223 "probability.dy"
Probability * Probability_move(const Probability * from,Probability * to,int len)
{
  register int i;

  for(i=0;i<len;i++)
    to[i] = from[i];

  return to;
}

/* Function:  Score_move(from,to,len)
 *
 * Descrip:    moves from to to 
 *
 *
 * Arg:        from [UNKN ] Score array with the numbers [const Score *]
 * Arg:          to [UNKN ] Score array to be written into [Score *]
 * Arg:         len [UNKN ] length [int]
 *
 * Return [UNKN ]  Undocumented return value [Score *]
 *
 */
# line 240 "probability.dy"
Score * Score_move(const Score * from,Score * to,int len)
{
  register int i;

  for(i=0;i<len;i++)
    to[i] = from[i];

  return to;
}


/* Function:  renormalise_Probability_array(array,len)
 *
 * Descrip:    Reasonably stupid function. Sums up probability array
 *             and then simply uses a linear renormalisation to get to the
 *             array adding to 1.0
 *
 *             returns the difference between the original sum and 1,0
 *
 *
 * Arg:        array [UNKN ] array to renormalise [Probability *]
 * Arg:          len [UNKN ] length of array [int]
 *
 * Return [UNKN ]  Undocumented return value [double]
 *
 */
# line 261 "probability.dy"
double renormalise_Probability_array(Probability * array,int len)
{
  register int i;
  double total;

  total = sum_Probability_array(array,len);

  for(i=0;i<len;i++)
    array[i] /= total;

  return (1.0 - total);
}


/* Function:  Probability_array_divide(to,top,bottem,len)
 *
 * Descrip:    divides one prob array by another pairwise
 *
 *
 * Arg:            to [UNKN ] probability array to be written into [Probability *]
 * Arg:           top [UNKN ] probability array to be divided [const Probability *]
 * Arg:        bottem [UNKN ] probability array which divides [const Probability *]
 * Arg:           len [UNKN ] length [int]
 *
 * Return [UNKN ]  Undocumented return value [Probability *]
 *
 */
# line 283 "probability.dy"
Probability * Probability_array_divide(Probability * to,const Probability * top,const Probability * bottem,int len)
{
  register int i;

  for(i=0;i<len;i++) 
    to[i] = top[i] / bottem[i];

  return to;
}


/* Function:  Probability_array_multiply(to,top,bottem,len)
 *
 * Descrip:    multiplies one prob array by another pairwise
 *
 *
 * Arg:            to [UNKN ] probability array to be written into [Probability *]
 * Arg:           top [UNKN ] probability array to be mulitpled [const Probability *]
 * Arg:        bottem [UNKN ] probability array to be mulitpled [const Probability *]
 * Arg:           len [UNKN ] length [int]
 *
 * Return [UNKN ]  Undocumented return value [Probability *]
 *
 */
# line 302 "probability.dy"
Probability * Probability_array_multiply(Probability * to,const Probability * top,const Probability * bottem,int len)
{
  register int i;

  for(i=0;i<len;i++) 
    to[i] = top[i] * bottem[i];

  return to;
}

/* Function:  Probability_array_add(to,top,bottem,len)
 *
 * Descrip:    sums one prob array with another pairwise
 *
 *
 * Arg:            to [UNKN ] probability array to be written into [Probability *]
 * Arg:           top [UNKN ] probability array to be summed [const Probability *]
 * Arg:        bottem [UNKN ] probability array to be summed [const Probability *]
 * Arg:           len [UNKN ] length [int]
 *
 * Return [UNKN ]  Undocumented return value [Probability *]
 *
 */
# line 320 "probability.dy"
Probability * Probability_array_add(Probability * to,const Probability * top,const Probability * bottem,int len)
{
  register int i;

  for(i=0;i<len;i++) 
    to[i] = top[i] + bottem[i];

  return to;
}


/* Function:  Probability_array_subtract(to,top,bottem,len)
 *
 * Descrip:    subtracts one prob array by another pairwise
 *
 *
 * Arg:            to [UNKN ] probability array to be written into [Probability *]
 * Arg:           top [UNKN ] probability array to be subtracted [const Probability *]
 * Arg:        bottem [UNKN ] probability array that subtracts [const Probability *]
 * Arg:           len [UNKN ] length [int]
 *
 * Return [UNKN ]  Undocumented return value [Probability *]
 *
 */
# line 339 "probability.dy"
Probability * Probability_array_subtract(Probability * to,const Probability * top,const Probability * bottem,int len)
{
  register int i;

  for(i=0;i<len;i++) 
    to[i] = top[i] - bottem[i];

  return to;
}


/* Function:  Score_array_add(to,top,bottem,len)
 *
 * Descrip:    sums one score array with another pairwise
 *
 *
 * Arg:            to [UNKN ] score array to be written into [Score *]
 * Arg:           top [UNKN ] score  array to be summed [Score *]
 * Arg:        bottem [UNKN ] score array to be summed [Score *]
 * Arg:           len [UNKN ] length [int]
 *
 * Return [UNKN ]  Undocumented return value [Score *]
 *
 */
# line 358 "probability.dy"
Score * Score_array_add(Score * to,Score * top,Score * bottem,int len)
{
  register int i;

  for(i=0;i<len;i++) 
    to[i] = top[i] + bottem[i];

  return to;
}


/* Function:  Score_array_subtract(to,top,bottem,len)
 *
 * Descrip:    subtracts one score array by another pairwise
 *
 *
 * Arg:            to [UNKN ] score array to be written into [Score *]
 * Arg:           top [UNKN ] score array to be subtracted [const Score *]
 * Arg:        bottem [UNKN ] score array that subtracts [const Score *]
 * Arg:           len [UNKN ] length [int]
 *
 * Return [UNKN ]  Undocumented return value [Score *]
 *
 */
# line 377 "probability.dy"
Score * Score_array_subtract(Score * to,const Score * top,const Score * bottem,int len)
{
  register int i;

  for(i=0;i<len;i++) 
    to[i] = top[i] - bottem[i];

  return to;
}

/* Function:  Score_Probability_sum(one,two)
 *
 * Descrip:    Badly implemented sum in probability
 *             space
 *
 *
 * Arg:        one [UNKN ] Undocumented argument [Score]
 * Arg:        two [UNKN ] Undocumented argument [Score]
 *
 * Return [UNKN ]  Undocumented return value [Score]
 *
 */
# line 391 "probability.dy"
Score Score_Probability_sum(Score one,Score two)
{
  return Probability2Score(Score2Probability(one) + Score2Probability(two));
}

/* Function:  Probability2Score(p)
 *
 * Descrip:    maps probabilities to scores. Deals
 *             sensibly with 0's.
 *
 *
 * Arg:        p [UNKN ] Undocumented argument [Probability]
 *
 * Return [UNKN ]  Undocumented return value [Score]
 *
 */
# line 400 "probability.dy"
Score Probability2Score(Probability p)
{
  if( p < PROBABILITY_MINIMUM )
    return NEGI;

  else return (Score) (INTEGER_FACTOR * (log(p)));
}

/* Function:  halfbit2Probability(half_bit)
 *
 * Descrip:    maps halfbits (log2(prob*2) to
 *             probabilities
 *
 *
 * Arg:        half_bit [UNKN ] Undocumented argument [double]
 *
 * Return [UNKN ]  Undocumented return value [Probability]
 *
 */
# line 412 "probability.dy"
Probability halfbit2Probability(double half_bit)
{

  if( half_bit <= NEGI + PROBABILITY_MINIMUM ) 
    return 0.0; /** just to guarentee conversion **/
  /*** ok half bit = 2 * log(P) / log(2.0) ***/

  return exp(half_bit * log(2.0) /2);
}

/* Function:  Bits2Probability(bits)
 *
 * Descrip:    maps halfbits (log2(prob*2) to
 *             probabilities
 *
 *
 * Arg:        bits [UNKN ] Undocumented argument [double]
 *
 * Return [UNKN ]  Undocumented return value [Probability]
 *
 */
# line 426 "probability.dy"
Probability Bits2Probability(double bits)
{

  if( bits <= NEGI ) 
    return 0.0; /** just to guarentee conversion **/
  /*** ok bits = log(P) / log(2.0) ***/

  return exp(bits * log(2.0));
}

/* Function:  Probability2halfbit(p)
 *
 * Descrip:    maps probabilities to halfbits.
 *             Deals with 0's sensibly
 *
 *
 * Arg:        p [UNKN ] Undocumented argument [Probability]
 *
 * Return [UNKN ]  Undocumented return value [double]
 *
 */
# line 440 "probability.dy"
double Probability2halfbit(Probability p)
{
  if( p < PROBABILITY_MINIMUM )
    return NEGI;
  return 2*log(p)/log(2.0);
}

  
/* Function:  Score2Probability(s)
 *
 * Descrip:    maps scores to probabilities
 *
 *
 * Arg:        s [UNKN ] Undocumented argument [Score]
 *
 * Return [UNKN ]  Undocumented return value [Probability]
 *
 */
# line 451 "probability.dy"
Probability Score2Probability(Score s)
{
  if( s <= NEGI )
    return 0.0;
  else return  exp((((double) s )/INTEGER_FACTOR ));
}

/* Function:  Score2Bits(s)
 *
 * Descrip:    maps scores to bits
 *
 *
 * Arg:        s [UNKN ] Undocumented argument [Score]
 *
 * Return [UNKN ]  Undocumented return value [Bits]
 *
 */
# line 461 "probability.dy"
Bits Score2Bits(Score s)
{
  return (Bits)((double)s / ((log((double)2.0))*INTEGER_FACTOR) );
}

# line 562 "probability.c"

#ifdef _cplusplus
}
#endif
