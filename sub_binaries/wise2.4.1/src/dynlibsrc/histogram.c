#ifdef _cplusplus
extern "C" {
#endif
#include "histogram.h"


/* Function:  Evalue_from_Histogram(his,score)
 *
 * Descrip: No Description
 *
 * Arg:          his [UNKN ] Histogram object [Histogram *]
 * Arg:        score [UNKN ] score you want the evalue for [double]
 *
 * Return [UNKN ]  Undocumented return value [double]
 *
 */
# line 85 "histogram.dy"
double Evalue_from_Histogram(Histogram * his,double score)
{
   return ExtremeValueE(score,his->param[EVD_MU],his->param[EVD_LAMBDA],his->total);
}

/* Function:  new_Histogram(min,max,lumpsize)
 *
 * Descrip: No Description
 *
 * Arg:             min [UNKN ] minimum score (integer) [int]
 * Arg:             max [UNKN ] maximum score (integer) [int]
 * Arg:        lumpsize [UNKN ] when reallocating histogram, the reallocation amount [int]
 *
 * Return [UNKN ]  Undocumented return value [Histogram *]
 *
 */
# line 113 "histogram.dy"
Histogram * new_Histogram(int min, int max, int lumpsize)
{
  Histogram *h;
  int            newsize;
  int            i;

  newsize = max - min + 1;

  h = Histogram_alloc(); /* ewan changed */
  if( h == NULL ) {
    return NULL;
  }

  h->min       = min;
  h->max       = max;
  h->total     = 0;
  h->lowscore  = INT_MAX;
  h->highscore = INT_MIN;
  h->lumpsize  = lumpsize;
  h->histogram = (int *) ckcalloc (newsize,sizeof(int)); /* ewan changed */
  for (i = 0; i < newsize; i++) h->histogram[i] = 0;

  h->expect    = NULL;
  h->fit_type  = HISTFIT_NONE;

  h->param[EVD_WONKA] = 1.0;	/* just in case, make sure this initializes */

  return h;
}


/* Function:  UnfitHistogram(h)
 *
 * Descrip: No Description
 *
 * Arg:        h [UNKN ] Undocumented argument [Histogram *]
 *
 */
# line 156 "histogram.dy"
void UnfitHistogram(Histogram * h)
{
  if (h->expect != NULL) 
     ckfree(h->expect);
  h->expect   = NULL;
  h->fit_type = HISTFIT_NONE;
}



/* Function:  AddToHistogram(h,sc)
 *
 * Descrip: No Description
 *
 * Arg:         h [UNKN ] Undocumented argument [Histogram *]
 * Arg:        sc [UNKN ] Undocumented argument [float]
 *
 */
# line 182 "histogram.dy"
void AddToHistogram(Histogram * h, float sc)
{
  int score;
  int moveby;
  int prevsize;
  int newsize;
  int i;

  /* Adding to a histogram conflicts with existing fit:
   * prohibit this.
   */
  if (h->fit_type != HISTFIT_NONE)
    fatal("AddToHistogram(): Can't add to a fitted histogram\n");
  

  /* histogram bins are defined as:  score >= bin value, < bin+1 
   * -1.9 -> -2    -0.4 -> -1    1.9 -> 1
   * -2.1 -> -3     0.4 -> 0     2.1 -> 2
   */
  score = (int) floor(sc);

  /* Check to see if we must reallocate the histogram.
   */
  if (score < h->min)
    {
      prevsize = h->max - h->min + 1;
      moveby   = (h->min - score) + h->lumpsize;
      newsize  = prevsize + moveby;
      h->min  -= moveby;

      h->histogram = (int *) ckrealloc(h->histogram, sizeof(int) * newsize);
      if( h->histogram == NULL ) {
	  fatal("Unable to extend histogram... have to crash... sorry!");
	  }          
      memmove(h->histogram+moveby, h->histogram, sizeof(int) * prevsize);
      for (i = 0; i < moveby; i++)
	h->histogram[i] = 0;
    }
  else if (score > h->max)
    {
      prevsize = h->max - h->min + 1;
      h->max   = h->lumpsize + score;
      newsize  = h->max - h->min + 1;

      h->histogram = (int *) ckrealloc(h->histogram, sizeof(int) * newsize);
      if( h->histogram == NULL ) {
          fatal("Cannot realloc histogram... going to die... sorry!");
	  }
      for (i = prevsize; i < newsize; i++)
	h->histogram[i] = 0;
    }

  /* Bump the correct bin.
   * The bin number is score - h->min
   */
  h->histogram[score - h->min]++;
  h->total++;
  if (score < h->lowscore) h->lowscore   = score;
  if (score > h->highscore) h->highscore = score;

#if DEBUG
  fprintf(stderr, "AddToHistogram(): added %.1f; rounded to %d; in bin %d (%d-%d)\n",
	  sc, score, score-h->min, h->min, h->max);
#endif
  return;
}


/* Function:  PrintASCIIHistogram(h,fp)
 *
 * Descrip: No Description
 *
 * Arg:         h [UNKN ] histogram to print [Histogram *]
 * Arg:        fp [UNKN ] open file to print to (stdout works) [FILE *]
 *
 */
# line 268 "histogram.dy"
void PrintASCIIHistogram(Histogram * h,FILE * fp)
{
  int units;
  int maxbar;
  int num;
  int i, idx;
  char buffer[81];		/* output line buffer */
  int  pos;			/* position in output line buffer */
  int  lowbound, lowcount;	/* cutoffs on the low side  */
  int  highbound, highcount;	/* cutoffs on the high side */
  int  emptybins = 3;

  /* Find out how we'll scale the histogram.
   * We have 59 characters to play with on a
   * standard 80-column terminal display:
   * leading "%5d %6d %6d|" occupies 20 chars.
   * Save the peak position, we'll use it later.
   */
  maxbar = 0;
  for (i = h->lowscore - h->min; i <= h->highscore - h->min; i++)
    if (h->histogram[i] > maxbar) 
      {
	maxbar   = h->histogram[i];     /* max height    */
	lowbound = i + h->min;     	/* peak position */
      }

  /* Truncate histogram display on both sides, ad hoc fashion.
   * Start from the peak; then move out until we see <emptybins> empty bins,
   * and stop.
   */
  highbound = lowbound;		/* start at peak position */
  for (num = 0; lowbound > h->lowscore; lowbound--)
    {
      i = lowbound - h->min;
      if (h->histogram[i] > 0) { num = 0;               continue; } /* reset */
      if (++num == emptybins)  { lowbound += emptybins; break;    } /* stop  */
    }
  for (num = 0; highbound < h->highscore; highbound++)
    {
      i = highbound - h->min;
      if (h->histogram[i] > 0) { num = 0;                continue; } /* reset */
      if (++num == emptybins)  { highbound -= emptybins; break;    } /* stop  */
    }
				/* collect counts outside of bounds */
  for (lowcount = 0, i = h->lowscore - h->min; i <= lowbound - h->min; i++)
    lowcount += h->histogram[i];
  for (highcount = 0, i = h->highscore - h->min; i >= highbound - h->min; i--)
    highcount += h->histogram[i];

				/* maxbar might need raised now; then set our units  */
  if (lowcount  > maxbar) maxbar = lowcount;
  if (highcount > maxbar) maxbar = highcount;
  units = ((maxbar-1)/ 59) + 1;


  /* Print the histogram
   */
  fprintf(fp, "%5s %6s %6s  (one = represents %d sequences)\n", 
	  "score", "obs", "exp", units);
  fprintf(fp, "%5s %6s %6s\n", "-----", "---", "---");
  buffer[80] = '\0';
  buffer[79] = '\n';
  for (i = h->lowscore; i <= h->highscore; i++)
    {
      memset(buffer, ' ', 79 * sizeof(char));
      idx = i - h->min;

      /* Deal with special cases at edges
       */
      if      (i < lowbound)  continue;
      else if (i > highbound) continue;
      else if (i == lowbound && i != h->lowscore) 
	{
	  sprintf(buffer, "<%4d %6d %6s|", i+1, lowcount, "-");
	  if (lowcount > 0) {
	    num = 1+(lowcount-1) / units;
	    if (num > 60) fatal("oops - more than 60 somethings in printing... ");
	    for (pos = 20; num > 0; num--)  buffer[pos++] = '=';
	  }
	  fputs(buffer, fp);
	  continue;
	}
      else if (i == highbound && i != h->highscore)
	{
	  sprintf(buffer, ">%4d %6d %6s|", i, highcount, "-");
	  if (highcount > 0) {
	    num = 1+(highcount-1) / units;
	    for (pos = 20; num > 0; num--)  buffer[pos++] = '=';
	  }
	  fputs(buffer, fp);
	  continue;
	}

      /* Deal with most cases
       */
      if (h->fit_type != HISTFIT_NONE) 
	sprintf(buffer, "%5d %6d %6d|", 
		i, h->histogram[idx], (int) h->expect[idx]);
      else
	sprintf(buffer, "%5d %6d %6s|", i, h->histogram[idx], "-");
      buffer[20] = ' ';		/* sprintf writes a null char */

      /* Mark the histogram bar for observed hits
       */ 
      if (h->histogram[idx] > 0) {
	num = 1 + (h->histogram[idx]-1) / units;
	for (pos = 20; num > 0; num--)  buffer[pos++] = '=';
      }
	  
      /* Mark the theoretically expected value
       */
      if (h->fit_type != HISTFIT_NONE && (int) h->expect[idx] > 0)
	{
				/* "corrected" line */
#ifdef SRE_REMOVED
	  if (h->fit_type == HISTFIT_EVD)
	    {
	      pos = 20 + (int)(h->param[EVD_WONKA] * h->expect[idx] - 1) / units;
	      if (pos >= 78) pos = 78; /* be careful of buffer bounds */
	      buffer[pos] = 'o';
	    }
#endif
				/* true (uncorrected) line */
	  pos = 20 + (int)(h->expect[idx]-1) / units;
	  if (pos >= 78) pos = 78; /* be careful of buffer bounds */
	  buffer[pos] = '*';
	}

      /* Print the line
       */
      fputs(buffer, fp);
    }

  /* Print details about the statistics
   */
  switch (h->fit_type) {
  case HISTFIT_NONE:
    fprintf(fp, "\n\n%% No statistical fit available\n");
    break;
    
  case HISTFIT_EVD:
    fprintf(fp, "\n\n%% Statistical details of theoretical EVD fit:\n");
    fprintf(fp, "              mu = %10.4f\n", h->param[EVD_MU]);
    fprintf(fp, "          lambda = %10.4f\n", h->param[EVD_LAMBDA]);
#ifdef SRE_REMOVED
    fprintf(fp, "    fraction fit = %10.4f\n", h->param[EVD_WONKA]);
#endif
    fprintf(fp, "chi-sq statistic = %10.4f\n", h->chisq);
    fprintf(fp, "  P(chi-square)  = %10.4g\n", h->chip);
    break;

  case HISTFIT_GAUSSIAN:
    fprintf(fp, "\n\n%% Statistical details of theoretical Gaussian fit:\n");
    fprintf(fp, "            mean = %10.4f\n", h->param[GAUSS_MEAN]);
    fprintf(fp, "              sd = %10.4f\n", h->param[GAUSS_SD]);
    fprintf(fp, "chi-sq statistic = %10.4f\n", h->chisq);
    fprintf(fp, "  P(chi-square)  = %10.4g\n", h->chip);
    break;
  }    
  return;
}
  


/* Function:  EVDBasicFit(h)
 *
 * Descrip: No Description
 *
 * Arg:        h [UNKN ] histogram to fit [Histogram *]
 *
 */
# line 456 "histogram.dy"
void EVDBasicFit(Histogram * h)
{
  float *d;            /* distribution P(S < x)          */
  float *x;            /* x-axis of P(S<x) for Linefit() */
  int    hsize;
  int    sum;
  int    sc, idx;		/* loop indices for score or score-h->min   */
  float  slope, intercept;	/* m,b fit from Linefit()                   */
  float  corr;			/* correlation coeff of line fit, not used  */
  float  lambda, mu;		/* slope, intercept converted to EVD params */

  /* Allocations for x, y axes
   * distribution d runs from min..max with indices 0..max-min
   *     i.e. score - min = index into d, x, histogram, and expect
   */
  hsize = h->highscore - h->lowscore + 1;
  d         = (float *) ckalloc(sizeof(float) * hsize);
  x         = (float *) ckalloc(sizeof(float) * hsize);
  for (idx = 0; idx < hsize; idx++)
    d[idx] = x[idx] = 0.;

  /* Calculate P(S < x) distribution from histogram.
   * note off-by-one of sc, because histogram bin contains scores between
   * x and x+1. 
   */ 
  sum = 0;
  for (sc = h->lowscore; sc <= h->highscore; sc++)
    {
      sum += h->histogram[sc - h->min];
      d[sc - h->lowscore] = (float) sum / (float) h->total;
      x[sc - h->lowscore] = (float) (sc + 1);
    }

  /* Do a linear regression fit to the log[-log(P(S<x))] "line".
   * we have log[-log(1-P(S>x))]  = -lambda * x + lambda * mu
   * so lambda = -m  and mu = b/lambda
   */
				/* convert y axis to log[-log(P(S<x))]  */
  for (sc = h->lowscore; sc < h->highscore; sc++)
    d[sc - h->lowscore] = log(-1. * log(d[sc - h->lowscore]));

				/* do the linear regression */
  Linefit(x, d, hsize-1, &intercept, &slope, &corr);
				/* calc mu, lambda */
  lambda = -1. * slope;
  mu     = intercept / lambda;

  /* Set the EVD parameters in the histogram;
   * pass 2 for additional lost degrees of freedom because we fit mu, lambda.
   */
  ExtremeValueSetHistogram(h, mu, lambda, h->lowscore, h->highscore, 1.0, 2);

  free(x);
  free(d);
  return;
}


/* Function:  ExtremeValueFitHistogram(h,censor,high_hint)
 *
 * Descrip: No Description
 *
 * Arg:                h [UNKN ] histogram to fit [Histogram *]
 * Arg:           censor [UNKN ] TRUE to censor data left of the peak [int]
 * Arg:        high_hint [UNKN ] score cutoff; above this are real hits that arent fit [float]
 *
 * Return [UNKN ]  if fit is judged to be valid else 0 if fit is invalid (too few seqs.) [int]
 *
 */
# line 540 "histogram.dy"
int ExtremeValueFitHistogram(Histogram * h, int censor, float high_hint) 
{
  float *x;                     /* array of EVD samples to fit */
  int   *y;                     /* histogram counts            */ 
  int    n;			/* number of observed samples  */
  int    z;			/* number of censored samples  */
  int    hsize;			/* size of histogram           */
  float  lambda, mu;		/* new estimates of lambda, mu */
  int    sc;		        /* loop index for score        */
  int    lowbound;		/* lower bound of fitted region*/
  int    highbound;		/* upper bound of fitted region*/
  int    new_highbound;
  int    iteration;

  /* Determine lower bound on fitted region;
   * if we're censoring the data, choose the peak of the histogram.
   * if we're not, then we take the whole histogram.
   */
  lowbound = h->lowscore;
  if (censor) 
    {
      int max = -1;
      for (sc = h->lowscore; sc <= h->highscore; sc++)
	if (h->histogram[sc - h->min] > max) 
	  {
	    max      = h->histogram[sc - h->min];
	    lowbound = sc;
	  }
    }

  /* Determine initial upper bound on fitted region.
   */
  highbound = MIN(high_hint, h->highscore);

  /* Now, iteratively converge on our lambda, mu:
   */
  for (iteration = 0; iteration < 100; iteration++)
    {
      /* Construct x, y vectors.
       */
      x = NULL;
      y = NULL;
      hsize = highbound - lowbound + 1;
      if (hsize < 5) {
	warn("On iteration %d, got %d bins, which is not fitable",iteration,hsize);
	goto FITFAILED; /* require at least 5 bins or we don't fit */
      }


      x = ckalloc(sizeof(float) * hsize);
      y = ckalloc(sizeof(int)   * hsize);
      if( x == NULL || y == NULL ) {
          warn("Out of temporary memory for evd fitting... exiting with error, though I'd worry about this");
	  return 0;
	  }

      n = 0;
      for (sc = lowbound; sc <= highbound; sc++)
	{
	  x[sc-lowbound] = (float) sc + 0.5; /* crude, but tests OK */
	  y[sc-lowbound] = h->histogram[sc - h->min];
	  n             += h->histogram[sc - h->min];
	}

      if (n < 100) {
	warn("On iteration %d, got only %d points, which is not fitable",iteration,n);
	goto FITFAILED;  /* require fitting to at least 100 points */
      }


      /* If we're censoring, estimate z, the number of censored guys
       * left of the bound. Our initial estimate is crudely that we're
       * missing e^-1 of the total distribution (which would be exact
       * if we censored exactly at mu; but we censored at the observed peak).
       * Subsequent estimates are more exact based on our current estimate of mu.
       */
      if (censor)
	{
	  if (iteration == 0)
	    z = MIN(h->total-n, (int) (0.58198 * (float) n));
	  else
	    {
	      double psx;
	      psx = EVDDistribution((float) lowbound, mu, lambda);
	      z = MIN(h->total-n, (int) ((double) n * psx / (1. - psx)));
	    }
	}

      /* Do an ML fit
       */
      if (censor) {
	if (! EVDCensoredFit(x, y, hsize, z, (float) lowbound, &mu, &lambda)) {
	  warn("On iteration %d, unable to make maxlikehood evd fit with censor",iteration);
	  goto FITFAILED;
	}
      } else {
	if (! EVDMaxLikelyFit(x, y, hsize, &mu, &lambda)) {
	  warn("On iteration %d, unable to make maxlikehood evd fit without censor",iteration);
	  goto FITFAILED;
	}
      }


      /* Find the Eval = 1 point as a new highbound;
       * the total number of samples estimated to "belong" to the EVD is n+z  
       */
      new_highbound = (int)
	(mu - (log (-1. * log((double) (n+z-1) / (double)(n+z))) / lambda));

      free(x);
      free(y);
      if (new_highbound >= highbound) break; 
      highbound = new_highbound;
    }

  /* Set the histogram parameters;
   * - the wonka factor is n+z / h->total : e.g. that's the fraction of the
   *   hits that we expect to match the EVD, others are generally lower
   * - we fit from lowbound to highbound; thus we lose 2 degrees of freedom
   *   for fitting mu, lambda, but we get 1 back because we're unnormalized
   *   in this interval, hence we pass 2-1 = 1 as ndegrees.
   *   
   *   Mon Jan 19 06:18:14 1998: wonka = 1.0, temporarily disabled.
   */    
  ExtremeValueSetHistogram(h, mu, lambda, lowbound, highbound, 1.0, 1); 
  return 1;

FITFAILED:
  UnfitHistogram(h);
  if (x != NULL) free(x);
  if (y != NULL) free(y);
  return 0;
}

    
/* Function:  ExtremeValueSetHistogram(h,mu,lambda,lowbound,highbound,wonka,ndegrees)
 *
 * Descrip: No Description
 *
 * Arg:                h [UNKN ] the histogram to set [Histogram *]
 * Arg:               mu [UNKN ] mu location parameter                 [float]
 * Arg:           lambda [UNKN ] lambda scale parameter [float]
 * Arg:         lowbound [UNKN ] low bound of the histogram that was fit [float]
 * Arg:        highbound [UNKN ] high bound of histogram that was fit [float]
 * Arg:            wonka [UNKN ] fudge factor; fraction of hits estimated to be "EVD-like" [float]
 * Arg:         ndegrees [UNKN ] extra degrees of freedom to subtract in chi2 test: [int]
 *
 */
# line 702 "histogram.dy"
void ExtremeValueSetHistogram(Histogram * h, float mu, float lambda, float lowbound, float highbound, float wonka, int ndegrees)
{
  int   sc;
  int   hsize, idx;
  int   nbins;
  float delta;

  UnfitHistogram(h);
  h->fit_type          = HISTFIT_EVD;
  h->param[EVD_LAMBDA] = lambda;
  h->param[EVD_MU]     = mu;
  h->param[EVD_WONKA]  = wonka;

  hsize     = h->max - h->min + 1;
  h->expect = (float *) ckalloc(sizeof(float) * hsize);
  if( h->expect == NULL ) {
     fatal("Cannot make memory for expect thing... ");
     }
  for (idx = 0; idx < hsize; idx++)
    h->expect[idx] = 0.;

  /* Calculate the expected values for the histogram.
   */
  for (sc = h->min; sc <= h->max; sc++)
    h->expect[sc - h->min] =
      ExtremeValueE((float)(sc), h->param[EVD_MU], h->param[EVD_LAMBDA], 
		    h->total) -
      ExtremeValueE((float)(sc+1), h->param[EVD_MU], h->param[EVD_LAMBDA],
		    h->total);
  
  /* Calculate the goodness-of-fit (within whole region)
   */
  h->chisq = 0.;
  nbins    = 0;
  for (sc = lowbound; sc <= highbound; sc++)
    if (h->expect[sc-h->min] >= 5. && h->histogram[sc-h->min] >= 5)
      {
	delta = (float) h->histogram[sc-h->min] - h->expect[sc-h->min];
	h->chisq += delta * delta / h->expect[sc-h->min];
	nbins++;
      }

  /* Since we fit the whole histogram, there is at least 
   * one constraint on chi-square: the normalization to h->total.
   */
  if (nbins > 1 + ndegrees)
    h->chip = (float) IncompleteGamma((double)(nbins-1-ndegrees)/2., 
				      (double) h->chisq/2.);
  else
    h->chip = 0.;		
}



/* Function:  GaussianFitHistogram(h,high_hint)
 *
 * Descrip: No Description
 *
 * Arg:                h [UNKN ] histogram to fit [Histogram *]
 * Arg:        high_hint [UNKN ] score cutoff; above this are `real' hits that aren't fit [float]
 *
 * Return [UNKN ]  if fit is judged to be valid else 0 if fit is invalid (too few seqs.)            [int]
 *
 */
# line 777 "histogram.dy"
int GaussianFitHistogram(Histogram * h, float high_hint)
{
  float sum;
  float sqsum;
  float delta;
  int   sc;
  int   nbins;
  int   hsize, idx;
  
  /* Clear any previous fitting from the histogram.
   */
  UnfitHistogram(h);

  /* Determine if we have enough hits to fit the histogram;
   * arbitrarily require 1000.
   */
  if (h->total < 1000) { h->fit_type = HISTFIT_NONE; return 0; }

  /* Simplest algorithm for mean and sd;
   * no outlier detection yet (not even using high_hint)
   * 
   * Magic 0.5 correction is because our histogram is for
   * scores between x and x+1; we estimate the expectation
   * (roughly) as x + 0.5. 
   */
  sum = sqsum = 0.;
  for (sc = h->lowscore; sc <= h->highscore; sc++)
    {
      delta  = (float) sc + 0.5;
      sum   += (float) h->histogram[sc-h->min] * delta;
      sqsum += (float) h->histogram[sc-h->min] * delta * delta;
    }
  h->fit_type          = HISTFIT_GAUSSIAN;
  h->param[GAUSS_MEAN] = sum / (float) h->total;
  h->param[GAUSS_SD]   = sqrt((sqsum - (sum*sum/(float)h->total)) / 
			      (float)(h->total-1));
  
  /* Calculate the expected values for the histogram.
   * Note that the magic 0.5 correction appears again.
   * Calculating difference between distribution functions for Gaussian 
   * would be correct but hard.
   */
  hsize     = h->max - h->min + 1;
  h->expect = (float *) ckalloc(sizeof(float) * hsize);
  if( h->expect == NULL ) {
      fatal("Unable to allocate expect space in histogram... sorry!");
  }
  for (idx = 0; idx < hsize; idx++)
    h->expect[idx] = 0.;

  for (sc = h->min; sc <= h->max; sc++)
    {
      delta = (float) sc + 0.5 - h->param[GAUSS_MEAN];
      h->expect[sc - h->min] =
	(float) h->total * ((1. / (h->param[GAUSS_SD] * sqrt(2.*3.14159))) * 
        (exp(-1.* delta*delta / (2. * h->param[GAUSS_SD] * h->param[GAUSS_SD]))));
    }

  /* Calculate the goodness-of-fit (within region that was fitted)
   */
  h->chisq = 0.;
  nbins    = 0;
  for (sc = h->lowscore; sc <= h->highscore; sc++)
    if (h->expect[sc-h->min] >= 5. && h->histogram[sc-h->min] >= 5)
      {
	delta = (float) h->histogram[sc-h->min] - h->expect[sc-h->min];
	h->chisq += delta * delta / h->expect[sc-h->min];
	nbins++;
      }
	/* -1 d.f. for normalization; -2 d.f. for two free parameters */
  if (nbins > 3)
    h->chip = (float) IncompleteGamma((double)(nbins-3)/2., 
				      (double) h->chisq/2.);
  else
    h->chip = 0.;		

  return 1;
}


/* Function:  GaussianSetHistogram(h,mean,sd)
 *
 * Descrip: No Description
 *
 * Arg:           h [UNKN ] Undocumented argument [Histogram *]
 * Arg:        mean [UNKN ] Undocumented argument [float]
 * Arg:          sd [UNKN ] Undocumented argument [float]
 *
 */
# line 871 "histogram.dy"
void GaussianSetHistogram(Histogram * h, float mean, float sd)
{
  int   sc;
  int   hsize, idx;
  int   nbins;
  float delta;

  UnfitHistogram(h);
  h->fit_type          = HISTFIT_GAUSSIAN;
  h->param[GAUSS_MEAN] = mean;
  h->param[GAUSS_SD]   = sd;

  /* Calculate the expected values for the histogram.
   */
  hsize     = h->max - h->min + 1;
  h->expect = (float *) ckalloc(sizeof(float) * hsize);
  if( h->expect == NULL ) {
      fatal("Unable to allocate expect size in expected histogram...");
  }

  for (idx = 0; idx < hsize; idx++)
    h->expect[idx] = 0.;

  /* Note: ideally we'd use the Gaussian distribution function
   * to find the histogram occupancy in the window sc..sc+1. 
   * However, the distribution function is hard to calculate.
   * Instead, estimate the histogram by taking the density at sc+0.5.
   */
  for (sc = h->min; sc <= h->max; sc++)
    { 
      delta = ((float)sc + 0.5) - h->param[GAUSS_MEAN];
      h->expect[sc - h->min] =
	(float) h->total * ((1. / (h->param[GAUSS_SD] * sqrt(2.*3.14159))) * 
	    (exp(-1.*delta*delta / (2. * h->param[GAUSS_SD] * h->param[GAUSS_SD]))));
    }

  /* Calculate the goodness-of-fit (within whole region)
   */
  h->chisq = 0.;
  nbins    = 0;
  for (sc = h->lowscore; sc <= h->highscore; sc++)
    if (h->expect[sc-h->min] >= 5. && h->histogram[sc-h->min] >= 5)
      {
	delta = (float) h->histogram[sc-h->min] - h->expect[sc-h->min];
	h->chisq += delta * delta / h->expect[sc-h->min];
	nbins++;
      }
	/* -1 d.f. for normalization */
  if (nbins > 1)
    h->chip = (float) IncompleteGamma((double)(nbins-1)/2., 
				      (double) h->chisq/2.);
  else
    h->chip = 0.;		
}



/* Function:  EVDDensity(x,mu,lambda)
 *
 * Descrip: No Description
 *
 * Arg:             x [UNKN ] Undocumented argument [float]
 * Arg:            mu [UNKN ] Undocumented argument [float]
 * Arg:        lambda [UNKN ] Undocumented argument [float]
 *
 * Return [UNKN ]  Undocumented return value [double]
 *
 */
# line 942 "histogram.dy"
double EVDDensity(float x, float mu, float lambda)
{
  return (lambda * exp(-1. * lambda * (x - mu) 
		       - exp(-1. * lambda * (x - mu))));
}

/* Function:  EVDDistribution(x,mu,lambda)
 *
 * Descrip: No Description
 *
 * Arg:             x [UNKN ] Undocumented argument [float]
 * Arg:            mu [UNKN ] Undocumented argument [float]
 * Arg:        lambda [UNKN ] Undocumented argument [float]
 *
 * Return [UNKN ]  Undocumented return value [double]
 *
 */
# line 962 "histogram.dy"
double EVDDistribution(float x, float mu, float lambda)
{
  return (exp(-1. * exp(-1. * lambda * (x - mu))));
}

/* Function:  ExtremeValueP(x,mu,lambda)
 *
 * Descrip: No Description
 *
 * Arg:             x [UNKN ] score [float]
 * Arg:            mu [UNKN ] characteristic value of extreme value distribution [float]
 * Arg:        lambda [UNKN ] decay constant of extreme value distribution [float]
 *
 * Return [UNKN ]  P(S>x) [double]
 *
 */
# line 988 "histogram.dy"
double ExtremeValueP(float x, float mu, float lambda)
{
  double y;
			/* a roundoff issue arises; use 1 - e^-x --> x for small x */
  y = exp(-1. * lambda * (x - mu));
  if (y < 1e-7) return y;
  else          return (1.0 - exp(-1. * y));
}


/* Function:  ExtremeValueP2(x,mu,lambda,N)
 *
 * Descrip: No Description
 *
 * Arg:             x [UNKN ] score [float]
 * Arg:            mu [UNKN ] characteristic value of extreme value distribution [float]
 * Arg:        lambda [UNKN ] decay constant of extreme value distribution [float]
 * Arg:             N [UNKN ] number of trials (number of sequences) [int]
 *
 * Return [UNKN ]  P(S>x) for database of size N [double]
 *
 */
# line 1018 "histogram.dy"
double ExtremeValueP2(float x, float mu, float lambda, int N)
{
  double y;
  y = N * ExtremeValueP(x,mu,lambda);
  if (y < 1e-7) return y;
  else          return (1.0 - exp(-1. * y));
}

/* Function:  ExtremeValueE(x,mu,lambda,N)
 *
 * Descrip: No Description
 *
 * Arg:             x [UNKN ] score [float]
 * Arg:            mu [UNKN ] characteristic value of extreme value distribution [float]
 * Arg:        lambda [UNKN ] decay constant of extreme value distribution [float]
 * Arg:             N [UNKN ] number of trials (number of sequences) [int]
 *
 * Return [UNKN ]  E(S>x) for database of size N [double]
 *
 */
# line 1046 "histogram.dy"
double ExtremeValueE(float x, float mu, float lambda, int N)
{
  return (double)N * ExtremeValueP(x,mu,lambda);
}


/* Function:  EVDrandom(mu,lambda)
 *
 * Descrip: No Description
 *
 * Arg:            mu [UNKN ] Undocumented argument [float]
 * Arg:        lambda [UNKN ] Undocumented argument [float]
 *
 * Return [UNKN ]  Undocumented return value [float]
 *
 */
# line 1068 "histogram.dy"
float EVDrandom(float mu, float lambda)
{
  float p = 0.0;

  /* Very unlikely, but possible,
   * that sre_random() would give us exactly 0 or 1 
   */
  while (p == 0. || p == 1.) p = sre_random(); 
  return mu - log(-1. * log(p)) / lambda;
} 
 

/* Function:  Lawless416(ret_f,ret_df,y,*x,x,*y,n,lambda,*ret_f,*ret_df)
 *
 * Descrip: No Description
 *
 * Arg:          ret_f [WRITE] RETURN: 4.1.6 evaluated at lambda [NullString]
 * Arg:         ret_df [WRITE] RETURN: first derivative of 4.1.6 evaluated at lambda [NullString]
 * Arg:              y [UNKN ] NULL (or y-axis of a histogram) [NullString]
 * Arg:             *x [UNKN ] Undocumented argument [float]
 * Arg:              x [UNKN ] array of sample values (or x-axis of a histogram) [NullString]
 * Arg:             *y [UNKN ] Undocumented argument [int]
 * Arg:              n [UNKN ] number of samples (or number of histogram bins) [int]
 * Arg:         lambda [UNKN ] a lambda to test [float]
 * Arg:         *ret_f [UNKN ] Undocumented argument [float]
 * Arg:        *ret_df [UNKN ] Undocumented argument [float]
 *
 */
# line 1107 "histogram.dy"
void Lawless416(float *x, int *y, int n, float lambda, float *ret_f, float *ret_df)
{

  double esum;			/* \sum e^(-lambda xi)      */
  double xesum;			/* \sum xi e^(-lambda xi)   */
  double xxesum;		/* \sum xi^2 e^(-lambda xi) */
  double xsum;			/* \sum xi                  */
  double mult;			/* histogram count multiplier */
  double total;			/* total samples            */
  int i;


  esum = xesum = xsum  = xxesum = total = 0.;
  for (i = 0; i < n; i++)
    {
      mult = (y == NULL) ? 1. : (double) y[i];
      xsum   += mult * x[i];
      xesum  += mult * x[i] * exp(-1. * lambda * x[i]);
      xxesum += mult * x[i] * x[i] * exp(-1. * lambda * x[i]);
      esum   += mult * exp(-1. * lambda * x[i]);
      total  += mult;
    }
  *ret_f  = 1./lambda - xsum / total + xesum / esum;
  *ret_df = ((xesum / esum) * (xesum / esum))
    - (xxesum / esum)
    - (1. / (lambda * lambda));

  return;
}


/* Function:  Lawless422(*x,ret_df,y,ret_f,x,*y,n,z,c,lambda,*ret_f,*ret_df)
 *
 * Descrip: No Description
 *
 * Arg:             *x [UNKN ] Undocumented argument [float]
 * Arg:         ret_df [WRITE] RETURN: first derivative of 4.2.2 evaluated at lambda [NullString]
 * Arg:              y [UNKN ] NULL (or y axis of a histogram) [NullString]
 * Arg:          ret_f [WRITE] RETURN: 4.2.2 evaluated at lambda [NullString]
 * Arg:              x [UNKN ] array of sample values (or x axis of a histogram) [NullString]
 * Arg:             *y [UNKN ] Undocumented argument [int]
 * Arg:              n [UNKN ] number of observed samples (or number of histogram bins) [int]
 * Arg:              z [UNKN ] number of censored samples  [int]
 * Arg:              c [UNKN ] censoring value; all observed x_i >= c          [float]
 * Arg:         lambda [UNKN ] a lambda to test [float]
 * Arg:         *ret_f [UNKN ] Undocumented argument [float]
 * Arg:        *ret_df [UNKN ] Undocumented argument [float]
 *
 */
# line 1170 "histogram.dy"
void Lawless422(float *x, int *y, int n, int z, float c,float lambda, float *ret_f, float *ret_df)
{
  double esum;			/* \sum e^(-lambda xi)      + z term    */
  double xesum;			/* \sum xi e^(-lambda xi)   + z term    */
  double xxesum;		/* \sum xi^2 e^(-lambda xi) + z term    */
  double xsum;			/* \sum xi                  (no z term) */
  double mult;			/* histogram count multiplier */
  double total;			/* total samples            */
  int i;

  esum = xesum = xsum  = xxesum = total = 0.;
  for (i = 0; i < n; i++)
    {
      mult = (y == NULL) ? 1. : (double) y[i];
      xsum   += mult * x[i];
      esum   += mult *               exp(-1. * lambda * x[i]);
      xesum  += mult * x[i] *        exp(-1. * lambda * x[i]);
      xxesum += mult * x[i] * x[i] * exp(-1. * lambda * x[i]);
      total  += mult;
    }

  /* Add z terms for censored data
   */
  esum   += (double) z *         exp(-1. * lambda * c);
  xesum  += (double) z * c *     exp(-1. * lambda * c);
  xxesum += (double) z * c * c * exp(-1. * lambda * c);

  *ret_f  = 1./lambda - xsum / total + xesum / esum;
  *ret_df = ((xesum / esum) * (xesum / esum))
    - (xxesum / esum)
    - (1. / (lambda * lambda));

  return;
}



/* Function:  EVDMaxLikelyFit(ret_lambda,x,*x,c,ret_mu,*c,n,*ret_mu,*ret_lambda)
 *
 * Descrip: No Description
 *
 * Arg:         ret_lambda [WRITE] RETURN: ML estimate of lambda [NullString]
 * Arg:                  x [UNKN ] list of EVD distributed samples or x axis of histogram [NullString]
 * Arg:                 *x [UNKN ] Undocumented argument [float]
 * Arg:                  c [UNKN ] NULL, or y-axis of histogram [NullString]
 * Arg:             ret_mu [WRITE] RETURN: ML estimate of mu [NullString]
 * Arg:                 *c [UNKN ] Undocumented argument [int]
 * Arg:                  n [UNKN ] number of samples, or number of histogram bins  [int]
 * Arg:            *ret_mu [UNKN ] Undocumented argument [float]
 * Arg:        *ret_lambda [UNKN ] Undocumented argument [float]
 *
 * Return [UNKN ]  on success; 0 on any failure [int]
 *
 */
# line 1238 "histogram.dy"
int EVDMaxLikelyFit(float *x, int *c, int n, float *ret_mu, float *ret_lambda)
{
  float  lambda, mu;
  float  fx;			/* f(x)  */
  float  dfx;			/* f'(x) */
  double esum;                  /* \sum e^(-lambda xi) */ 
  double mult;
  double total;
  float  tol = 1e-5;
  int    i;

  /* 1. Find an initial guess at lambda: linear regression here?
   */
  lambda = 0.2;

  /* 2. Use Newton/Raphson to solve Lawless 4.1.6 and find ML lambda
   */
  for (i = 0; i < 100; i++)
    {
      Lawless416(x, c, n, lambda, &fx, &dfx);
      if (fabs(fx) < tol) break;             /* success */
      lambda = lambda - fx / dfx;	     /* Newton/Raphson is simple */
      if (lambda <= 0.) lambda = 0.001;      /* but be a little careful  */
    }

  /* 2.5: If we did 100 iterations but didn't converge, Newton/Raphson failed.
   *      Resort to a bisection search. Worse convergence speed
   *      but guaranteed to converge (unlike Newton/Raphson).
   *      We assume (!?) that fx is a monotonically decreasing function of x;
   *      i.e. fx > 0 if we are left of the root, fx < 0 if we
   *      are right of the root.
   */ 
  if (i == 100)
    {
      float left, right, mid;
      info("EVD maxlik fit - Newton/Raphson failed; switchover to bisection");

				/* First we need to bracket the root */
      lambda = right = left = 0.2;
      Lawless416(x, c, n, lambda, &fx, &dfx);
      if (fx < 0.) 
	{			/* fix right; search left. */
	  do {
	    left -= 0.1;
	    if (left < 0.) { 
	      info("failed to bracket root"); 
	      return 0; 
	    }
	    Lawless416(x, c, n, left, &fx, &dfx);
	  } while (fx < 0.);
	}
      else
	{			/* fix left; search right. */
	  do {
	    right += 0.1;
	    Lawless416(x, c, n, right, &fx, &dfx);
	    if (right > 100.) {
	      info("failed to bracket root"); 
	      return 0; 
	    }
	  } while (fx > 0.);
	}
			/* now we bisection search in left/right interval */
      for (i = 0; i < 100; i++)
	{
	  mid = (left + right) / 2.; 
	  Lawless416(x, c, n, mid, &fx, &dfx);
	  if (fabs(fx) < tol) break;             /* success */
	  if (fx > 0.)	left = mid;
	  else          right = mid;
	}
      if (i == 100) { 
	warn("even the bisection search failed"); 
	return 0; 
      }
      lambda = mid;
    }

  /* 3. Substitute into Lawless 4.1.5 to find mu
   */
  esum = 0.;
  for (i = 0; i < n; i++)
    {
      mult   = (c == NULL) ? 1. : (double) c[i];
      esum  += mult * exp(-1 * lambda * x[i]);
      total += mult;
    }
  mu = -1. * log(esum / total) / lambda;

  *ret_lambda = lambda;
  *ret_mu     = mu;   
  return 1;
}


/* Function:  EVDCensoredFit(ret_mu,x,ret_lambda,y,*x,*y,n,z,c,*ret_mu,*ret_lambda)
 *
 * Descrip: No Description
 *
 * Arg:             ret_mu [WRITE] RETURN: ML estimate of mu [NullString]
 * Arg:                  x [UNKN ] list of EVD distributed samples or x axis of histogram [NullString]
 * Arg:         ret_lambda [WRITE] RETURN: ML estimate of lambda [NullString]
 * Arg:                  y [UNKN ] NULL, or y axis of histogram [NullString]
 * Arg:                 *x [UNKN ] Undocumented argument [float]
 * Arg:                 *y [UNKN ] Undocumented argument [int]
 * Arg:                  n [UNKN ] number of observed samples,or number of histogram bins  [int]
 * Arg:                  z [UNKN ] number of censored samples [int]
 * Arg:                  c [UNKN ] censoring value (all x_i >= c) [float]
 * Arg:            *ret_mu [UNKN ] Undocumented argument [float]
 * Arg:        *ret_lambda [UNKN ] Undocumented argument [float]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
# line 1366 "histogram.dy"
int EVDCensoredFit(float *x, int *y, int n, int z, float c,float *ret_mu, float *ret_lambda)
{
  float  lambda, mu;
  float  fx;			/* f(x)  */
  float  dfx;			/* f'(x) */
  double esum;                  /* \sum e^(-lambda xi) */ 
  double mult;
  double total;
  float  tol = 1e-5;
  int    i;

  /* 1. Find an initial guess at lambda: linear regression here?
   */
  lambda = 0.2;

  /* 2. Use Newton/Raphson to solve Lawless 4.2.2 and find ML lambda
   */
  for (i = 0; i < 100; i++)
    {
      Lawless422(x, y, n, z, c, lambda, &fx, &dfx);
      if (fabs(fx) < tol) break;             /* success */
      lambda = lambda - fx / dfx;	     /* Newton/Raphson is simple */
      if (lambda <= 0.) lambda = 0.001;      /* but be a little careful  */
    }

 /* 2.5: If we did 100 iterations but didn't converge, Newton/Raphson failed.
   *      Resort to a bisection search. Worse convergence speed
   *      but guaranteed to converge (unlike Newton/Raphson).
   *      We assume (!?) that fx is a monotonically decreasing function of x;
   *      i.e. fx > 0 if we are left of the root, fx < 0 if we
   *      are right of the root.
   */ 
  if (i == 100)
    {
      float left, right, mid;
				/* First we need to bracket the root */
      info("EVDCensor fit: Newton/Raphson failed; switched to bisection");
      lambda = right = left = 0.2;
      Lawless422(x, y, n, z, c, lambda, &fx, &dfx);
      if (fx < 0.) 
	{			/* fix right; search left. */
	  do {
	    left -= 0.03;
	    if (left < 0.) { 
	      info("failed to bracket root"); 
	      return 0;
	    }
	    Lawless422(x, y, n, z, c, left, &fx, &dfx);
	  } while (fx < 0.);
	}
      else
	{			/* fix left; search right. */
	  do {
	    right += 0.1;
	    Lawless422(x, y, n, z, c, left, &fx, &dfx);
	    if (right > 100.) {
	      info("failed to bracket root");
	      return 0;
	    }
	  } while (fx > 0.);
	}
			/* now we bisection search in left/right interval */
      for (i = 0; i < 100; i++)
	{
	  mid = (left + right) / 2.; 
	  Lawless422(x, y, n, z, c, left, &fx, &dfx);
	  if (fabs(fx) < tol) break;             /* success */
	  if (fx > 0.)	left = mid;
	  else          right = mid;
	}
      if (i == 100) {
	info("even the bisection search failed");
	return 0;
      }
      lambda = mid;
    }

  /* 3. Substitute into Lawless 4.2.3 to find mu
   */
  esum =  total = 0.;
  for (i = 0; i < n; i++)
    {
      mult   = (y == NULL) ? 1. : (double) y[i];
      esum  += mult * exp(-1. * lambda * x[i]);
      total += mult;
    }
  esum += (double) z * exp(-1. * lambda * c);    /* term from censored data */
  mu = -1. * log(esum / total) / lambda;        

  *ret_lambda = lambda;
  *ret_mu     = mu;   
  return 1;
}

/* Function:  Linefit(ret_r,x,*x,ret_b,y,ret_a,*y,N,*ret_a,*ret_b,*ret_r)
 *
 * Descrip: No Description
 *
 * Arg:         ret_r [WRITE] RETURN: correlation coefficient    [NullString]
 * Arg:             x [UNKN ] x values of data [NullString]
 * Arg:            *x [UNKN ] Undocumented argument [float]
 * Arg:         ret_b [WRITE] RETURN: slope [NullString]
 * Arg:             y [UNKN ] y values of data                [NullString]
 * Arg:         ret_a [WRITE] RETURN: intercept [NullString]
 * Arg:            *y [UNKN ] Undocumented argument [float]
 * Arg:             N [UNKN ] number of data points [int]
 * Arg:        *ret_a [UNKN ] Undocumented argument [float]
 * Arg:        *ret_b [UNKN ] Undocumented argument [float]
 * Arg:        *ret_r [UNKN ] Undocumented argument [float]
 *
 * Return [UNKN ]  on success, 0 on failure. [int]
 *
 */
# line 1484 "histogram.dy"
int Linefit(float *x, float *y, int N, float *ret_a, float *ret_b, float *ret_r) 
{				
  float xavg, yavg;
  float sxx, syy, sxy;
  int   i;
  
  /* Calculate averages, xavg and yavg
   */
  xavg = yavg = 0.0;
  for (i = 0; i < N; i++)
    {
      xavg += x[i];
      yavg += y[i];
    }
  xavg /= (float) N;
  yavg /= (float) N;

  sxx = syy = sxy = 0.0;
  for (i = 0; i < N; i++)
    {
      sxx    += (x[i] - xavg) * (x[i] - xavg);
      syy    += (y[i] - yavg) * (y[i] - xavg);
      sxy    += (x[i] - xavg) * (y[i] - yavg);
    }
  *ret_b = sxy / sxx;
  *ret_a = yavg - xavg*(*ret_b);
  *ret_r = sxy / (sqrt(sxx) * sqrt(syy));
  return 1;
}


#define RANGE 268435456		/* 2^28        */
#define DIV   16384		/* sqrt(RANGE) */
#define MULT  72530821		/* my/Cathy's birthdays, x21, x even (Knuth)*/


static int sre_reseed = 0;	/* TRUE to reinit sre_random() */
static int sre_randseed = 666;	/* default seed for sre_random()   */

/* Function:  sre_random(void)
 *
 * Descrip: No Description
 *
 *
 * Return [UNKN ]  Undocumented return value [float]
 *
 */
# line 1553 "histogram.dy"
float sre_random(void)
{
  static long  rnd;
  static int   firsttime = 1;
  long         high1, low1;
  long         high2, low2; 

  if (sre_reseed || firsttime) 
    {
      sre_reseed = firsttime = 0;
      if (sre_randseed <= 0) sre_randseed = 666; /* seeds of zero break me */
      high1 = sre_randseed / DIV;  low1  = sre_randseed % DIV;
      high2 = MULT / DIV;          low2  = MULT % DIV;
      rnd = (((high2*low1 + high1*low2) % DIV)*DIV + low1*low2) % RANGE;
    }
   high1 = rnd / DIV;  low1  = rnd % DIV;
   high2 = MULT / DIV; low2  = MULT % DIV;
   rnd = (((high2*low1 + high1*low2) % DIV)*DIV + low1*low2) % RANGE;

   return ((float) rnd / (float) RANGE);  
}
#undef RANGE
#undef DIV
#undef MULT

/* Function:  IncompleteGamma(a,x)
 *
 * Descrip: No Description
 *
 * Arg:        a [UNKN ] Undocumented argument [double]
 * Arg:        x [UNKN ] Undocumented argument [double]
 *
 * Return [UNKN ]  Undocumented return value [double]
 *
 */
# line 1609 "histogram.dy"
double IncompleteGamma(double a, double x)
{
  int iter;			/* iteration counter */

  if (a <= 0.) fatal("IncompleteGamma(): a must be > 0");
  if (x <  0.) fatal("IncompleteGamma(): x must be >= 0");

  /* For x > a + 1 the following gives rapid convergence;
   * calculate 1 - P(a,x) = \frac{\Gamma(a,x)}{\Gamma(a)}:
   *     use a continued fraction development for \Gamma(a,x).
   */
  if (x > a+1) 
    {
      double oldp;		/* previous value of p    */
      double nu0, nu1;		/* numerators for continued fraction calc   */
      double de0, de1;		/* denominators for continued fraction calc */

      nu0 = 0.;			/* A_0 = 0       */
      de0 = 1.;			/* B_0 = 1       */
      nu1 = 1.;			/* A_1 = 1       */
      de1 = x;			/* B_1 = x       */

      oldp = nu1;
      for (iter = 1; iter < 100; iter++)
	{
	  /* Continued fraction development:
	   * set A_j = b_j A_j-1 + a_j A_j-2
	   *     B_j = b_j B_j-1 + a_j B_j-2
           * We start with A_2, B_2.
	   */
				/* j = even: a_j = iter-a, b_j = 1 */
				/* A,B_j-2 are in nu0, de0; A,B_j-1 are in nu1,de1 */
	  nu0 = nu1 + ((double)iter - a) * nu0;
	  de0 = de1 + ((double)iter - a) * de0;

				/* j = odd: a_j = iter, b_j = x */
				/* A,B_j-2 are in nu1, de1; A,B_j-1 in nu0,de0 */
	  nu1 = x * nu0 + (double) iter * nu1;
	  de1 = x * de0 + (double) iter * de1;

				/* rescale */
	  if (de1) 
	    { 
	      nu0 /= de1; 
	      de0 /= de1;
	      nu1 /= de1;
	      de1 =  1.;
	    }
				/* check for convergence */
	  if (fabs((nu1-oldp)/nu1) < 1.e-7)
	    return nu1 * exp(a * log(x) - x - Gammln(a));

	  oldp = nu1;
	}
      fatal("IncompleteGamma(): failed to converge using continued fraction approx");
    }
  else /* x <= a+1 */
    {
      double p;			/* current sum               */
      double val;		/* current value used in sum */

      /* For x <= a+1 we use a convergent series instead:
       *   P(a,x) = \frac{\gamma(a,x)}{\Gamma(a)},
       * where
       *   \gamma(a,x) = e^{-x}x^a \sum_{n=0}{\infty} \frac{\Gamma{a}}{\Gamma{a+1+n}} x^n
       * which looks appalling but the sum is in fact rearrangeable to
       * a simple series without the \Gamma functions:
       *   = \frac{1}{a} + \frac{x}{a(a+1)} + \frac{x^2}{a(a+1)(a+2)} ...
       * and it's obvious that this should converge nicely for x <= a+1.
       */
      
      p = val = 1. / a;
      for (iter = 1; iter < 10000; iter++)
	{
	  val *= x / (a+(double)iter);
	  p   += val;
	  
	  if (fabs(val/p) < 1.e-7)
	    return 1. - p * exp(a * log(x) - x - Gammln(a));
	}
      fatal("IncompleteGamma(): failed to converge using series approx");
    }
  /*NOTREACHED*/
  return 0.;
}


/* Function:  Gammln(x)
 *
 * Descrip: No Description
 *
 * Arg:        x [UNKN ] Undocumented argument [float]
 *
 * Return [UNKN ]  Undocumented return value [float]
 *
 */
# line 1715 "histogram.dy"
float Gammln(float x)
{
  int i;
  double xx, tx;
  double tmp, value;
  static double cof[11] = {
    4.694580336184385e+04,
    -1.560605207784446e+05,
    2.065049568014106e+05,
    -1.388934775095388e+05,
    5.031796415085709e+04,
    -9.601592329182778e+03,
    8.785855930895250e+02,
    -3.155153906098611e+01,
    2.908143421162229e-01,
    -2.319827630494973e-04,
    1.251639670050933e-10
  };
  
  /* Protect against x=0. We see this in Dirichlet code,
   * for terms alpha = 0. This is a severe hack but it is effective
   * and safe. (due to GJM)
   */ 
  if (x <= 0.0) return 999999.; 

  xx       = x - 1.0;
  tx = tmp = xx + 11.0;
  value    = 1.0;
  for (i = 10; i >= 0; i--)	/* sum least significant terms first */
    {
      value += cof[i] / tmp;
      tmp   -= 1.0;
    }
  value  = log(value);
  tx    += 0.5;
  value += 0.918938533 + (xx+0.5)*log(tx) - tx;
  return (float) value;
}


# line 1416 "histogram.c"
/* Function:  hard_link_Histogram(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [Histogram *]
 *
 * Return [UNKN ]  Undocumented return value [Histogram *]
 *
 */
Histogram * hard_link_Histogram(Histogram * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a Histogram object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  Histogram_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Histogram *]
 *
 */
Histogram * Histogram_alloc(void) 
{
    Histogram * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(Histogram *) ckalloc (sizeof(Histogram))) == NULL)  {  
      warn("Histogram_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->histogram = NULL;   
    out->min = 0;    
    out->max = 0;    
    out->highscore = 0;  
    out->lowscore = 0;   
    out->lumpsize = 0;   
    out->total = 0;  
    out->expect = NULL;  
    out->fit_type = 0;   
    /* param[3] is an array: no default possible */ 
    out->chisq = 0;  
    out->chip = 0;   


    return out;  
}    


/* Function:  free_Histogram(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [Histogram *]
 *
 * Return [UNKN ]  Undocumented return value [Histogram *]
 *
 */
Histogram * free_Histogram(Histogram * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a Histogram obj. Should be trappable"); 
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
    if( obj->histogram != NULL)  
      ckfree(obj->histogram);    
    if( obj->expect != NULL) 
      ckfree(obj->expect);   


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif
