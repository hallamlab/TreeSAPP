#ifdef _cplusplus
extern "C" {
#endif
#include "wiserandom.h"

static boolean isinit = FALSE;

/* Function:  init_random(void)
 *
 * Descrip:    initates the random generator to
 *             time byte...
 *
 *
 *
 */
# line 20 "wiserandom.dy"
void init_random(void)
{
  srandom((long) time(NULL));
  isinit = TRUE;
}

/* Function:  random_integer(l)
 *
 * Descrip:    returns an integer between 0 and l
 *             though I don't think we will get 0
 *             very often. Hmmm.
 *
 *
 * Arg:        l [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
# line 31 "wiserandom.dy"
int random_integer(int l)
{
  double rand;
  
  rand = random_0_to_1();
  
  return (int) (l*rand);
}

/* Function:  random_0_to_1(void)
 *
 * Descrip:    returns a random number between
 *             0 and 1
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [double]
 *
 */
# line 44 "wiserandom.dy"
double random_0_to_1(void)
{
  double ret;
  long rnd;

  if( isinit == FALSE)
    init_random();
  isinit = TRUE;

  rnd = random();

  ret = (double)rnd / (double)LONG_MAX;

/*  fprintf(stderr," got... %D %f\n",rnd,ret);*/

  return ret;
}
	
# line 67 "wiserandom.c"

#ifdef _cplusplus
}
#endif
