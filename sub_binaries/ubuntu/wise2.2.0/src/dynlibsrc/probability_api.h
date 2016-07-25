

/* Helper functions in the module
 *
 * Wise2_Probability_from_average_state_occupancy
 * Wise2_state_occupancy_from_Probability
 * Wise2_Probability2Score
 * Wise2_Score2Probability
 * Wise2_Score2Bits
 * Wise2_halfbit2Probability
 *



/* These functions are not associated with an object */
/* Function:  Wise2_Probability_from_average_state_occupancy(length)
 *
 * Descrip:    for single state (exponetial decays) takes an average length
 *             and converts that to a probability that will produce that
 *             length (on average) for the state. NB... this *assumes* that
 *             you want a single state exp decay.
 *
 *
 * Arg:        length       average length of state [double]
 *
 * Returns Undocumented return value [Probability]
 *
 */
Probability Wise2_Probability_from_average_state_occupancy( double length);

/* Function:  Wise2_state_occupancy_from_Probability(p)
 *
 * Descrip:    If you have a single state then this will tak
 *             the probability for the state->state transition and
 *             give you back the average length in the state
 *
 *
 * Arg:        p            probability of staying in the state [double]
 *
 * Returns Undocumented return value [double]
 *
 */
double Wise2_state_occupancy_from_Probability( double p);

/* Function:  Wise2_Probability2Score(p)
 *
 * Descrip:    maps probabilities to scores. Deals
 *             sensibly with 0's.
 *
 *
 * Arg:        p            Undocumented argument [Probability]
 *
 * Returns Undocumented return value [Score]
 *
 */
Score Wise2_Probability2Score( Probability p);

/* Function:  Wise2_Score2Probability(s)
 *
 * Descrip:    maps scores to probabilities
 *
 *
 * Arg:        s            Undocumented argument [Score]
 *
 * Returns Undocumented return value [Probability]
 *
 */
Probability Wise2_Score2Probability( Score s);

/* Function:  Wise2_Score2Bits(s)
 *
 * Descrip:    maps scores to bits
 *
 *
 * Arg:        s            Undocumented argument [Score]
 *
 * Returns Undocumented return value [Bits]
 *
 */
Bits Wise2_Score2Bits( Score s);

/* Function:  Wise2_halfbit2Probability(half_bit)
 *
 * Descrip:    maps halfbits (log2(prob*2) to
 *             probabilities
 *
 *
 * Arg:        half_bit     Undocumented argument [double]
 *
 * Returns Undocumented return value [Probability]
 *
 */
Probability Wise2_halfbit2Probability( double half_bit);

