#ifndef DYNAMITEprobabilityHEADERFILE
#define DYNAMITEprobabilityHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif

#include "wisebase.h"

#define LOG_NEGATIVE_INFINITY (-100000000)
#define NEGI LOG_NEGATIVE_INFINITY
#define PROBABILITY_MINIMUM (0.0000000000000000001)

#define INTEGER_FACTOR 500

typedef double Probability;
typedef int    Score;
typedef double Bits;

#define MINIMUM_ERROR (0.00000001)



    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



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
Probability Wise2_Probability_from_average_state_occupancy(double length);
#define Probability_from_average_state_occupancy Wise2_Probability_from_average_state_occupancy


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
double Wise2_state_occupancy_from_Probability(double p);
#define state_occupancy_from_Probability Wise2_state_occupancy_from_Probability


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
Score Wise2_Probability_logsum(Score one,Score two);
#define Probability_logsum Wise2_Probability_logsum


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
Bits Wise2_Probability2Bits(Probability p);
#define Probability2Bits Wise2_Probability2Bits


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
void Wise2_show_Score_array(Score * s,int len,FILE * ofp);
#define show_Score_array Wise2_show_Score_array


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
void Wise2_show_Probability_array_exp(Probability * p,int len,FILE * ofp);
#define show_Probability_array_exp Wise2_show_Probability_array_exp


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
boolean Wise2_read_Probability_array(Probability * p,int len,char * start_of_array);
#define read_Probability_array Wise2_read_Probability_array


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
void Wise2_show_Probability_array(Probability * p,int len,FILE * ofp);
#define show_Probability_array Wise2_show_Probability_array


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
Probability Wise2_sum_Probability_array(Probability * p,int len);
#define sum_Probability_array Wise2_sum_Probability_array


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
Probability * Wise2_set_Probability_array(Probability * set,Probability p,int len);
#define set_Probability_array Wise2_set_Probability_array


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
Score * Wise2_Probability2Score_move(Probability * from,Score * to,int len);
#define Probability2Score_move Wise2_Probability2Score_move


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
Probability * Wise2_Probability_move(const Probability * from,Probability * to,int len);
#define Probability_move Wise2_Probability_move


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
Score * Wise2_Score_move(const Score * from,Score * to,int len);
#define Score_move Wise2_Score_move


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
double Wise2_renormalise_Probability_array(Probability * array,int len);
#define renormalise_Probability_array Wise2_renormalise_Probability_array


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
Probability * Wise2_Probability_array_divide(Probability * to,const Probability * top,const Probability * bottem,int len);
#define Probability_array_divide Wise2_Probability_array_divide


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
Probability * Wise2_Probability_array_multiply(Probability * to,const Probability * top,const Probability * bottem,int len);
#define Probability_array_multiply Wise2_Probability_array_multiply


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
Probability * Wise2_Probability_array_add(Probability * to,const Probability * top,const Probability * bottem,int len);
#define Probability_array_add Wise2_Probability_array_add


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
Probability * Wise2_Probability_array_subtract(Probability * to,const Probability * top,const Probability * bottem,int len);
#define Probability_array_subtract Wise2_Probability_array_subtract


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
Score * Wise2_Score_array_add(Score * to,Score * top,Score * bottem,int len);
#define Score_array_add Wise2_Score_array_add


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
Score * Wise2_Score_array_subtract(Score * to,const Score * top,const Score * bottem,int len);
#define Score_array_subtract Wise2_Score_array_subtract


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
Score Wise2_Score_Probability_sum(Score one,Score two);
#define Score_Probability_sum Wise2_Score_Probability_sum


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
Score Wise2_Probability2Score(Probability p);
#define Probability2Score Wise2_Probability2Score


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
Probability Wise2_halfbit2Probability(double half_bit);
#define halfbit2Probability Wise2_halfbit2Probability


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
Probability Wise2_Bits2Probability(double bits);
#define Bits2Probability Wise2_Bits2Probability


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
double Wise2_Probability2halfbit(Probability p);
#define Probability2halfbit Wise2_Probability2halfbit


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
Probability Wise2_Score2Probability(Score s);
#define Score2Probability Wise2_Score2Probability


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
Bits Wise2_Score2Bits(Score s);
#define Score2Bits Wise2_Score2Bits


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/

#ifdef _cplusplus
}
#endif

#endif
