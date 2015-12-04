/* ARIADNE V.1.0

 GAPSTAT - functions for calculating the statistics of gapped
 alignment scores 


 Copyright 

 Richard Mott 2000
 
 Wellcome Trust Centre For Human Genetics
 Roosevelt Drive
 Oxford OX3 7AD
 

*/

#include<math.h>
#include<limits.h>
#include<stdio.h>
#include<stdlib.h>
#include"gapstat.h"


/************************************************************************/
/* 1. */
/* Karlin-Altschul statistics for HSPs */

/************************************************************************/

int 
KarlinAltschulStatistics( int **matrix, double *freq1, double *freq2, double *lambda, double *Kminus, double *Kplus, double *H, double *r, double *s ) {

/* compute all necessary quantities for ungapped statistics 

   inputs are the matrix and symbol frequencies freq1, freq2

   outputs are 
   lambda            the exponential rate 
   Kminus, Kplus     bounds on K
   H                 the entropy (for edge corrections)
   r, s              quantities needed for GEM gapped statistics

*/

  double *h;
  int hmin, hmax;
  double R1, R2, R3;
  double mean;

  h = get_h( matrix, freq1, freq2, &hmin, &hmax, &mean );

  if ( (*lambda = solve_for_lambda( h, hmin, hmax )) > 0 ) {
    
    *H = entropy_H(h,hmin,hmax,*lambda);
    
    iglehart( h, hmin, hmax, *lambda, -1, 1, &R1, &R2, &R3, Kplus, Kminus );
    
    *Kplus /= R3;
    *Kminus /= R3;
    
    *s = R1;
    *r = R2;

    /*  printf("lambda=%8.6f K-=%8.6f K+=%8.6f H=%8.4f r=%8.4f s=%8.4f\n", *lambda,  *Kplus, *Kminus, *H, *r, *s); */
    return 1;
  }
  else
    return 0; /* error condition */
    }

double
lambda_func( double lambda, double *h, int hmin, int hmax) {

  int c;
  double d, l;
  double f = 0.0;

  l = exp(lambda);
  d = exp(lambda*hmin);

  for(c=hmin;c<=hmax;c++) {
    f += d*h[c];
    d *= l;
  }
  return f-1.0;
}

/* Solve the equation

1 = sum_{s} exp(lambda s) h(s) 

to get the exponential rate constant for HSPs */


double
solve_for_lambda( double *h, int hmin, int hmax  ) {

  double left=0.001, right=1.0, midpt;
  double fleft, fright;
  double fmidpt; /* lambda, f;*/

  fleft = lambda_func( left, h, hmin, hmax );
  fright = lambda_func( right, h, hmin, hmax );

  if ( fleft*fright > 0.0 ) {
    fprintf(stderr, "ERROR could not bracket root, %g %g   %g %g\n", left, fleft, right, fright);
    return -1.0;
  }

  while( right-left > 1.0e-6 ) {
    midpt = (left+right)/2.0;
    fmidpt = lambda_func( midpt, h, hmin, hmax );
    if ( fmidpt < 0.0 ) {
      left = midpt;
      fleft = fmidpt;
    }
    else {
      right = midpt;
      fright = fmidpt;
    }
  }

  return midpt;
}


double entropy_H( double *h, int hmin, int hmax, double lambda ) {

/* compute the entropy of a score */

  int i; /*, j;*/
  double H=0.0;

  double l = exp(lambda);
  double d = exp(lambda*hmin);
  
  for(i=hmin;i<=hmax;i++) {
      H += h[i]*d*i;
      d *= l;
  }

  return H;
}

double HSP_length_correction( double H, double K, double lambda, int len1, int len2 )
{

/* Altschul-Gish length correction formula */

  double L_HSP;

  L_HSP = log(K*len1*len2)/(H*lambda);
/*  printf("E_HSP %8.4f L_HSP %8.4f K_HSP %8.4f\n", H, L_HSP, log(K*len1*len2) ); */
  return L_HSP;
}


/* gapped alignment statistics  */

int gem_statistics( double gap_start, double gap_extend, double lambda, double s, double Kplus, double Kminus, double *alpha, double *theta1, double *theta2, double *K1, double *K2 ) {


  /* alpha and theta */

  /* if *alpha > 0 then don't compute it */
  if ( *alpha <= 0 ) { 
    if ( lambda*gap_extend > 0.0 ) 
      *alpha = 2*s*exp(-lambda*(gap_start+gap_extend))/(1-exp(-lambda*gap_extend));
    else 
      *alpha = 1.0e10;
  }

  if ( *alpha < 0.44 ) { /* logarithmic domain */
    *theta2 = gapped_theta( *alpha, &upper_objective );
    *theta1 = gapped_theta( *alpha, &lower_objective );

/*    printf("theta1 = %.4f theta2 = %.4f\n", *theta1, *theta2 );
    printf("lambda1 = %.4f lambda2 = %.4f\n", lambda* *theta1, lambda* *theta2 );
*/    
    /* upper bound */

    *K2 = Kplus * *theta2 * *theta2;

    /* lower bound */

    if ( *theta1 < 1.0 ) 
      *K1 = Kminus* *theta1 * *theta1 * (1-exp(- *alpha))/(1- *theta1);
    else 
      *K1 = Kminus;

    return 1;
  }
  else { /* linear domain */
    *theta2 = *theta1  = *K1 = *K2 = 1.0e-8;
    return 0;
  }
}


/************************************************************************/
/* 2. General Routines */
/************************************************************************/

/* compute the distribution h(x), the probablity that two residues
   have score x, given the current score matrix and residue frequency
   distributions */

double *get_h( int **matrix, double *freq1, double *freq2, int *hmin, int *hmax, double *mean ) {

  int m, n;
  double *distribution;

  *hmax = -10000;
  *hmin = +10000;
  *mean = 0.0;

  for(n=0;n<256;n++)
    for(m=0;m<256;m++) {
      if ( matrix[n][m] > *hmax )
        *hmax = matrix[n][m];
      if ( matrix[n][m] < *hmin )
        *hmin = matrix[n][m];
    }

  distribution = (double*)calloc((*hmax)- (*hmin)+1,sizeof(double))- *hmin;

  /*  printf("min %d max %d\n", *hmin, *hmax ); */

  for(n=0;n<256;n++)
    for(m=0;m<256;m++)
      distribution[matrix[m][n]] += freq1[n]*freq2[m];

  for(n=*hmin;n<=*hmax;n++)
    *mean += distribution[n]*n;

  /*  {
    double s=0.0;
    printf( "mean = %e\n", *mean );
    for(n=*hmin;n<=*hmax;n++) {
      printf("%d dist %e %e\n", n, distribution[n], s=s+distribution[n]);
    }
  }
  */
  return distribution;
}

/* Solve the Associated Dam Equation */

double associated_dam_eqn( double *distribution, int min, int max, double lambda, double 
**transient, int *cn) {

  int its=100;
  double *f = (double*)calloc(max-min+1,sizeof(double))-min;
  double **h = (double**)calloc(its,sizeof(double*));
  double *converged;
  double x, r=0.0;
  int t, i, Max=0, Max1, k, hi;
  double d, l;

  h[1] = (double*)calloc(max+1,sizeof(double));
  h[1][max] = 0.0;
  for(i=max-1;i>=0;i--)
    h[1][i] = h[1][i+1] + distribution[i+1];

  r=0.0;
  d = exp(lambda*min);
  l = exp(lambda);
  for(i=min;i<=max;i++) {
    f[i] = distribution[i]*d;
    d *= l;
    r += f[i];
  }

  d = 1.0;
  for(i=0;i<=max;i++){
    h[1][i] *= d;
    d *= l;
  }

  for(k=2;k<its;k++) {
    Max=k*max;
    Max1 = (k-1)*max;
    h[k] = (double*)calloc(k*(max+1),sizeof(double));
    for(t=0;t<=Max;t++) {
      hi =  t > max ? max : t ;
      h[k][t] =  t > max ? 0.0 : h[1][t] ;
      for(i=hi;i>=min && (t-i) <= Max1;i--) {
        h[k][t] += f[i]*h[k-1][t-i];
      }
    }
  }

  r = 0.0;
  k--;
  *cn=0;
  for(t=0;t<=Max;t++) {
    if ( h[k][t] > 0.0 ) {
      x = h[k-1][t]/h[k][t];
      if ( fabs(1.0-x) < 0.0001 ) {
        (*cn)++;
/*      printf("%5d %8.6f %8.6f converged\n", t, h[k][t], x ); */
      }
      else
        break;
    }
  }
  converged = (double*)calloc((*cn)+1,sizeof(double));
  for(t=0;t<=*cn;t++)
    converged[t] = h[k][t];

  *transient = converged;

  while(t>0 && h[k][t] < h[k][t-1]) t--;
  r = h[k][t];
  /*  printf("t=%d r=%g\n", t, r); */
  return r;
}


/* compute objective function for upper_bound */

double upper_objective( double alpha, double beta ) {

  double f = beta*log(alpha) + lgamma(1.0-beta);
  return f;
}

/* compute objective function for upper_bound */

double lower_objective( double alpha, double beta ) {

  double f = log(2.0) + beta*log(alpha/2.0) -log(1.0-beta) -log(2.0-beta);
  return f;
}

/* solve for theta(alpha) */

double gapped_theta( double alpha, double (*objective)( double, double) ) {

  /*double beta;*/
  double left, right, midpt;
  double fleft, fright, fmidpt;

  left = 0.001;
  right = 0.999;

  fleft = (*objective)( alpha, left );
  fright = (*objective)( alpha, right );

  if ( fleft*fright >= 0.0 ) {
    if ( alpha < 1.0e-3 ) 
      return 1.0; 
    else {
      fprintf(stderr, "ERROR could not bracket root %e  %e   %e %e alpha %e\n", left, fleft, right, fright, alpha);
      return 0.0;
    }
  }

  while( right-left > 1.0e-6 ) {

    midpt = (left+right)/2.0;
    fmidpt = (*objective)( alpha, midpt );

    if ( fleft < fright ) {
      if ( fmidpt < 0.0 ) {
        left = midpt;
        fleft = fmidpt;
      }
      else {
        right = midpt;
        fright = fmidpt;
      }
    }
    else {
      if ( fmidpt > 0.0 ) {
        left = midpt;
        fleft = fmidpt;
      }
      else {
        right = midpt;
        fright = fmidpt;
      }
    }
    /*    printf("left %e right %e\n", left, right );  */
  }
  return midpt;
}



void iglehart( double *h, int hmin, int hmax, double lambda, double xmu, int delta, double *R1, double *R2, double *R3, double *Kplus, double *Kminus ) {
  
/*
   Use Igleharts formulae to compute the constants R1, R2, R3 for a
   random walk with negative drift

   If S(n) = sum_{1}^{n} X(i), then

   A = sum_{n} S(n) has asymptotic distribution

   P( A > t ) ~ R1 exp( -lambda t )

   [R1 == s in GEM statistics ]

   If N is the time of first entry into negative values, then

   B = max{n=1}^{N} S(n)

   and

   P( B > t ) ~ R2 exp( -lambda t )

   [R2 == s in GEM stats]

   Finally R3 = E(N), expected time before hitting negative values

   In what follows:

   Xmu = E( X exp(lambda X) )

   ESN = 1 - E( exp( lambda S(N) ) )
   EN = E(N)

   
   Then R1 = [ 1 - E( exp( lambda S(N) ) ) ] / [ L * Xmu * E(N) ]
           = ESN / ( lambda * Xmu * EN )

   and R2 = R1 * [ 1 - E( exp( lambda S(N) ) ) ]
          = R1 * ESN
      
   R3 = E(N)

   also:

   L = lambda for cts distributions, or

   L = (exp(lambda*delta)-1)/delta for arithmetic distributions, where delta is the span

   EN = exp ( sum_{k>0} Pr( S(k)>=0 ) / k )
   
   ESN = exp ( -sum_{k>0} E( exp( lambda S(k)) ; S(k)<0 ) / k ) 
*/

  double Xmu;
  double ESN;
  double EN;
  double SN;
  double last_EN=1.0e10;
  double last_ESN=1.0e10;

  double *s, *t;

  int n, k, i;
  int min, max;
  int min1, max1;
  double mean, expect, total;
  double d, eld, l=exp(lambda);
  double *s_ptr=NULL;
  double *t_ptr;

/* Xmu */

  if ( xmu > 0.0 ) {
    Xmu = xmu;
  }
  else {
    Xmu = 0.0;
    mean = 0.0;
    d = exp(lambda*hmin);
    for(k=hmin;k<=hmax;k++) {
      Xmu += h[k]*k*d;
      mean += h[k]*d;
      d *= l;
    }

  }

/* EN and ESN */

  t_ptr = (double*)calloc(hmax-hmin+1,sizeof(double));
  t = t_ptr-hmin;

  for(k=hmin;k<=hmax;k++)
    t[k] = h[k];

  SN = 0.0;
  for(k=0;k<=hmax;k++)
    SN += t[k];
  EN = SN;

  ESN = 0.0;
  
  SN = 0.0;
  d = exp(lambda*hmin);
  for(k=hmin;k<0;k++) {
    SN += t[k]*d;
    d *= l;
  }
  
  ESN = SN;

  mean = 0.0;
  total=0.0;
  for(k=hmin;k<=hmax;k++) {
    mean += t[k]*k;
    total += t[k];
  }


  for(n=2;n<100;n++) { 

    min = n*hmin;
    max = n*hmax;

    s_ptr = (double*)calloc(max-min+1,sizeof(double));
    s = s_ptr - min;

    for(k=min;k<=max;k++) {
      min1 = hmin;
      if ( k-min1 > (n-1)*hmax ) min1 = k-(n-1)*hmax;
      max1 = hmax;
      if ( k-max1 < (n-1)*hmin ) max1 = k-(n-1)*hmin;
      
      for(i=min1;i<=max1;i++) {
	s[k] += h[i]*t[k-i];
      }
      
    }

    SN = 0.0;
    for(k=0;k<=max;k++)
      SN += s[k];
    SN /= n;
    EN += SN;

    SN = 0.0;
    d = exp(lambda*min);
    for(k=min;k<0;k++){ 
      SN += s[k]*d;
      d *= l;
    }
    
    SN /= n;

    ESN += SN;

    expect = 0.0;
    total =0.0;
    for(k=min;k<=max;k++) {
      expect += k*s[k];
      total += s[k];
    }


    if ( fabs(last_EN-EN) < 1.0e-4 && fabs(last_ESN-ESN) < 1.0e-4 ) {
      free(t_ptr);
      break;
    }

    last_EN = EN;
    last_ESN = ESN;

    free(t_ptr);
    t_ptr = s_ptr;
    t = s;
  }

  free(s_ptr);

  EN = exp(EN);
  ESN = exp(-ESN);

  eld = exp(lambda*delta);
  *Kplus = ESN / ( (1-1.0/eld) * Xmu * EN / delta );  /*  K+ */
  *Kminus = *R1 = ESN / ( (eld-1.0) * Xmu * EN / delta );  /* K- */

  *R2 = *R1 *ESN;

  *Kplus *= ESN;
  *Kminus *= ESN;

  *R3 = EN;

  /*  printf( "ESN %8.6f EN %8.6f Xmu %8.6f R1 %8.6f R2 %8.6f LD %8.6f\n",  
	 ESN, EN, Xmu, *R1, *R2, (exp(lambda*delta)-1)/delta ); */
}

 
/* Empirical GEM statistics */

double EmpiricalGEM( double lambda0, double K0, double H, double alpha, int len1, int len2, double *theta, double *logkappa, double score ) {

  double alpha2 = alpha*alpha;
  double Kappa, Lambda, pval, x, f;
  double L = len1*len2;
  f = log(L)*(1.0/(double)len1 + 1.0/(double)len2);

  *theta = 1.013 -2.61*alpha + f*( -0.76 + 9.34*alpha +1.12/H);
  *logkappa = 0.26 -18.92*alpha + f*(-1.76 + 32.69*alpha +  192.52*alpha2 +  3.24/H );

  Kappa = K0*exp(*logkappa);
  Lambda = lambda0*(*theta);
  
  x = Kappa*L*exp(-Lambda*score);
  
  if ( x < 1.0e-6 )
    pval = x;
  else
    pval = 1.0-exp(-x);
  
  /*  printf("alpha %.5f lambda0 %.5f K0 %.5f H %.5f f %.5f theta %.5f logkappa %.5f pval %e %5d %5d x %e\n", alpha, lambda0, K0, H, f, *theta, *logkappa, pval, len1, len2, x );  */
    
  return pval;
}


double AutoIncreaseGapPenalty( double s, double lambda0, double alpha_max, int *A, int *B ) {

  double alpha = 2*s*exp(-lambda0*(*A + *B))/(1-exp(-lambda0* *B));
  int its=0;

  while( its++ < 100 ) {
    (*A)++;
    alpha = 2*s*exp(-lambda0*(*A + *B))/(1-exp(-lambda0* *B));
    if ( alpha < alpha_max ) break;
    (*B)++;
    alpha = 2*s*exp(-lambda0*(*A + *B))/(1-exp(-lambda0* *B));
    if ( alpha < alpha_max ) break;
  }
  printf( "WARNING: gap penalty increased to %d+%dk\n", *A, *B );
  return alpha;
}

