

MODULE = Wise2 PACKAGE = Wise2

Probability
Probability_from_average_state_occupancy(length)
	double length
	CODE:
	RETVAL = Wise2_Probability_from_average_state_occupancy(length);
	OUTPUT:
	RETVAL



double
state_occupancy_from_Probability(p)
	double p
	CODE:
	RETVAL = Wise2_state_occupancy_from_Probability(p);
	OUTPUT:
	RETVAL



Score
Probability2Score(p)
	Probability p
	CODE:
	RETVAL = Wise2_Probability2Score(p);
	OUTPUT:
	RETVAL



Probability
Score2Probability(s)
	Score s
	CODE:
	RETVAL = Wise2_Score2Probability(s);
	OUTPUT:
	RETVAL



Bits
Score2Bits(s)
	Score s
	CODE:
	RETVAL = Wise2_Score2Bits(s);
	OUTPUT:
	RETVAL



Probability
halfbit2Probability(half_bit)
	double half_bit
	CODE:
	RETVAL = Wise2_halfbit2Probability(half_bit);
	OUTPUT:
	RETVAL



