#ifndef DYNAMITEhistogramHEADERFILE
#define DYNAMITEhistogramHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "wisebase.h"

#define HISTFIT_NONE     0	/* no fit done yet               */
#define HISTFIT_EVD      1	/* fit type = extreme value dist */
#define HISTFIT_GAUSSIAN 2	/* fit type = Gaussian           */
#define EVD_MU		 0	/* EVD fit parameter mu          */
#define EVD_LAMBDA       1	/* EVD fit parameter lambda      */
#define EVD_WONKA        2      /* EVD fit fudge factor          */
#define GAUSS_MEAN       0	/* Gaussian parameter mean       */
#define GAUSS_SD         1	/* Gaussian parameter std. dev.  */

#ifndef MIN
#define MIN(a,b) ((a)<(b) ? (a) : (b))
#endif

/* Object Histogram
 *
 * Descrip: This Object came from Sean Eddy excellent histogram package.
 *        He donated it free of all restrictions to allow it to be used
 *        in the Wise2 package without complicated licensing terms.
 *        He is a *very* nice man.
 *
 *        It was made into a dynamite object so that
 *           a) External ports to scripting languages would be trivial
 *           b) cooperation with future versions of histogram.c would be possible.
 *
 *        Here is the rest of the documentation from sean.
 *
 *        Keep a score histogram. 
 *
 *        The main implementation issue here is that the range of
 *        scores is unknown, and will go negative. histogram is
 *        a 0..max-min array that represents the range min..max.
 *        A given score is indexed in histogram array as score-min.
 *        The AddToHistogram function deals with dynamically 
 *        resizing the histogram array when necessary.
 *
 *
 */
struct Wise2_Histogram {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int*   histogram;   /*  counts of hits                      */ 
    int min;    /*  elem 0 of histogram == min          */ 
    int max;    /*  last elem of histogram == max       */ 
    int highscore;  /*  highest active elem has this score  */ 
    int lowscore;   /*  lowest active elem has this score   */ 
    int lumpsize;   /*  when resizing, overalloc by this    */ 
    int total;  /*  total # of hits counted             */ 
    float* expect;  /*  expected counts of hits             */ 
    int fit_type;   /*  flag indicating distribution type   */ 
    float param[3]; /*  parameters used for fits            */ 
    float chisq;    /*  chi-squared val for goodness of fit */ 
    float chip; /*  P value for chisquared              */ 
    } ;  
/* Histogram defined */ 
#ifndef DYNAMITE_DEFINED_Histogram
typedef struct Wise2_Histogram Wise2_Histogram;
#define Histogram Wise2_Histogram
#define DYNAMITE_DEFINED_Histogram
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



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
double Wise2_Evalue_from_Histogram(Histogram * his,double score);
#define Evalue_from_Histogram Wise2_Evalue_from_Histogram


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
Histogram * Wise2_new_Histogram(int min, int max, int lumpsize);
#define new_Histogram Wise2_new_Histogram


/* Function:  UnfitHistogram(h)
 *
 * Descrip: No Description
 *
 * Arg:        h [UNKN ] Undocumented argument [Histogram *]
 *
 */
void Wise2_UnfitHistogram(Histogram * h);
#define UnfitHistogram Wise2_UnfitHistogram


/* Function:  AddToHistogram(h,sc)
 *
 * Descrip: No Description
 *
 * Arg:         h [UNKN ] Undocumented argument [Histogram *]
 * Arg:        sc [UNKN ] Undocumented argument [float]
 *
 */
void Wise2_AddToHistogram(Histogram * h, float sc);
#define AddToHistogram Wise2_AddToHistogram


/* Function:  PrintASCIIHistogram(h,fp)
 *
 * Descrip: No Description
 *
 * Arg:         h [UNKN ] histogram to print [Histogram *]
 * Arg:        fp [UNKN ] open file to print to (stdout works) [FILE *]
 *
 */
void Wise2_PrintASCIIHistogram(Histogram * h,FILE * fp);
#define PrintASCIIHistogram Wise2_PrintASCIIHistogram


/* Function:  EVDBasicFit(h)
 *
 * Descrip: No Description
 *
 * Arg:        h [UNKN ] histogram to fit [Histogram *]
 *
 */
void Wise2_EVDBasicFit(Histogram * h);
#define EVDBasicFit Wise2_EVDBasicFit


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
int Wise2_ExtremeValueFitHistogram(Histogram * h, int censor, float high_hint) ;
#define ExtremeValueFitHistogram Wise2_ExtremeValueFitHistogram


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
void Wise2_ExtremeValueSetHistogram(Histogram * h, float mu, float lambda, float lowbound, float highbound, float wonka, int ndegrees);
#define ExtremeValueSetHistogram Wise2_ExtremeValueSetHistogram


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
int Wise2_GaussianFitHistogram(Histogram * h, float high_hint);
#define GaussianFitHistogram Wise2_GaussianFitHistogram


/* Function:  GaussianSetHistogram(h,mean,sd)
 *
 * Descrip: No Description
 *
 * Arg:           h [UNKN ] Undocumented argument [Histogram *]
 * Arg:        mean [UNKN ] Undocumented argument [float]
 * Arg:          sd [UNKN ] Undocumented argument [float]
 *
 */
void Wise2_GaussianSetHistogram(Histogram * h, float mean, float sd);
#define GaussianSetHistogram Wise2_GaussianSetHistogram


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
double Wise2_EVDDensity(float x, float mu, float lambda);
#define EVDDensity Wise2_EVDDensity


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
double Wise2_EVDDistribution(float x, float mu, float lambda);
#define EVDDistribution Wise2_EVDDistribution


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
double Wise2_ExtremeValueP(float x, float mu, float lambda);
#define ExtremeValueP Wise2_ExtremeValueP


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
double Wise2_ExtremeValueP2(float x, float mu, float lambda, int N);
#define ExtremeValueP2 Wise2_ExtremeValueP2


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
double Wise2_ExtremeValueE(float x, float mu, float lambda, int N);
#define ExtremeValueE Wise2_ExtremeValueE


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
float Wise2_EVDrandom(float mu, float lambda);
#define EVDrandom Wise2_EVDrandom


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
void Wise2_Lawless416(float *x, int *y, int n, float lambda, float *ret_f, float *ret_df);
#define Lawless416 Wise2_Lawless416


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
void Wise2_Lawless422(float *x, int *y, int n, int z, float c,float lambda, float *ret_f, float *ret_df);
#define Lawless422 Wise2_Lawless422


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
int Wise2_EVDMaxLikelyFit(float *x, int *c, int n, float *ret_mu, float *ret_lambda);
#define EVDMaxLikelyFit Wise2_EVDMaxLikelyFit


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
int Wise2_EVDCensoredFit(float *x, int *y, int n, int z, float c,float *ret_mu, float *ret_lambda);
#define EVDCensoredFit Wise2_EVDCensoredFit


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
int Wise2_Linefit(float *x, float *y, int N, float *ret_a, float *ret_b, float *ret_r) ;
#define Linefit Wise2_Linefit


/* Function:  sre_random(void)
 *
 * Descrip: No Description
 *
 *
 * Return [UNKN ]  Undocumented return value [float]
 *
 */
float Wise2_sre_random(void);
#define sre_random Wise2_sre_random


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
double Wise2_IncompleteGamma(double a, double x);
#define IncompleteGamma Wise2_IncompleteGamma


/* Function:  Gammln(x)
 *
 * Descrip: No Description
 *
 * Arg:        x [UNKN ] Undocumented argument [float]
 *
 * Return [UNKN ]  Undocumented return value [float]
 *
 */
float Wise2_Gammln(float x);
#define Gammln Wise2_Gammln


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
Histogram * Wise2_hard_link_Histogram(Histogram * obj);
#define hard_link_Histogram Wise2_hard_link_Histogram


/* Function:  Histogram_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Histogram *]
 *
 */
Histogram * Wise2_Histogram_alloc(void);
#define Histogram_alloc Wise2_Histogram_alloc


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
Histogram * Wise2_free_Histogram(Histogram * obj);
#define free_Histogram Wise2_free_Histogram


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
