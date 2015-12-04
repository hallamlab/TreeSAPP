#include<stdio.h>
#include<limits.h>
#include<math.h>
#include"gapstat.h"

double default_lambda0, default_K0, default_Kplus, default_H, default_r, default_s, default_alpha;
int count=0;
int gap_init=0;
int gap_extend=0;
int **matrix=NULL;
double pval_threshold=0;
int initialised=0;

int KarlinAltschulStatistics2( int **matrix, double *freq1, double *freq2, double *lambda, double *Kminus, double *Kplus, double *H, double *r, double *s );
double * PseudoResidueFrequencies2( char *seq, int len, int pseudo_len );
double *RobinsonResidueFrequencies2 (void);
double *GetHistogram( double *freq1, double *freq2, int **matrix, int *hmin, int *hmax, double *mean );


int InitPvaluesMott( int **Matrix, int Gap_init, int Gap_extend, double Pval_threshold ) {

  /* sets up the P-value calculation. Should be called once at the
     start of the program. The substitution matrix and gap penalties
     are passed in at this point and cached. 

     Note the a gap of length k costs

     Gap_init + k*Gap_extend;

     Matrix is a **int for which matrix[i][j] is the substitution
     score for symbols i, j. Only uppercase letters 'A' through 'Z'
     matter.
  
     Pval_threshold should be a fairly small value such as 1.0e-5,
     such that pvalues greater than this are only evaluated relative
     to the default protein sequence composition. This speeds things
     up a lot.
     
     return value is an int:

     0 => FAILURE, because mean substitution score is positive: You
     must change the substitutio matrix and re-evaluate InitPvaluesMott()

     1 => SUCCESS, everything is OK to commence the search.

     2 => WARNING: the computed P-values are likely to be inaccurate
     because the scoring scheme is lax. You are recommended to
     increase the gap penalty and re-evaluate InitPvaluesMott().

  */

  double *freq0 = RobinsonResidueFrequencies2();
  matrix = Matrix;
  gap_init = Gap_init;
  gap_extend = Gap_extend;
  pval_threshold = Pval_threshold;

  if ( (initialised = KarlinAltschulStatistics2( matrix, freq0, freq0, &default_lambda0, &default_K0, &default_Kplus, &default_H, &default_r, &default_s )) ) {
    default_alpha = 2*default_s*exp(-default_lambda0*(gap_init+gap_extend))/(1-exp(-default_lambda0*gap_extend));
    if ( default_alpha > 0.25 ) 
      initialised = 2;
  }

  free(freq0);
  return initialised;
}

double SW_PValueMott( double SWScore, char *seqA, char *seqB, int lengthA, int lengthB, int *ok ) {
 
  /* returns the PValue for the SWScore, for the two sequences.
     
  ok == 1 means the P-value is reliable
  ok == 0 means a problem was encountered - most probably the mean symbol match was positive for the sequences, so no stats can be computed
  ok == 2 means the Pvalue was computed but is probably unreliable becuase the scoring scheme is too lax for this comparison
  ok == -1 means the initialisation was not performed, and so no pvalue was computed.
  */
  
  double pval1, pval2;
  double theta, logkappa;
  double lambda0, K0, Kplus, H, r, s, alpha;
  
  *ok = -1;
  pval2 = 1.0;
  if ( initialised ) {
    pval1 = pval2 = EmpiricalGEM( default_lambda0, default_K0, default_H, default_alpha, lengthA, lengthB, &theta, &logkappa, SWScore );
    count++;
    *ok = 1;
    /*  printf ( "pval1 %e default_alpha %e threshold %e\n", pval1, default_alpha, pval_threshold ); */
    if ( pval1 < pval_threshold ) {
      double *freqA = PseudoResidueFrequencies2( seqA, lengthA, 100 );
      double *freqB = PseudoResidueFrequencies2( seqB, lengthB, 100 );
      if ( KarlinAltschulStatistics2( matrix, freqA, freqB, &lambda0, &K0, &Kplus, &H, &r, &s ) ) {
	alpha = 2*s*exp(-lambda0*(gap_init+gap_extend))/(1-exp(-lambda0*gap_extend));
	if ( alpha > 0.25 ) 
	  *ok = 2;
	else 
	  *ok = 1;
	pval2 = EmpiricalGEM( lambda0, K0, H, alpha, lengthA, lengthB, &theta, &logkappa, SWScore );
      }
      else {
	*ok = 0;
	pval2 = 1.0;
      }
      free(freqA);
      free(freqB);
    }
  }
  return pval2;
}
  
int KarlinAltschulStatistics2( int **matrix, double *freq1, double *freq2, double *lambda, double *Kminus, double *Kplus, double *H, double *r, double *s ) {

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

  h = GetHistogram( freq1, freq2, matrix, &hmin, &hmax, &mean );
  
  if ( mean < 0.0 && ( *lambda = solve_for_lambda( h, hmin, hmax ))  > 0 ) {
    
    *H = entropy_H(h,hmin,hmax,*lambda);
    
    iglehart( h, hmin, hmax, *lambda, -1, 1, &R1, &R2, &R3, Kplus, Kminus );
    
    *Kplus /= R3;
    *Kminus /= R3;
    
    *s = R1;
    *r = R2;
    
    free(h+hmin);

    /*    fprintf(stderr, "lambda=%8.6f K-=%8.6f K+=%8.6f H=%8.4f r=%8.4f s=%8.4f\n", *lambda,  *Kplus, *Kminus, *H, *r, *s);  */
    return 1;
  }
  else {
    int n;
    /*    fprintf( stderr, "mean : %e h distribution:\n", mean);
        for(n=hmin;n<=hmax;n++){
      fprintf(stderr, "%5d %e\n", n, h[n]); 
    } */
    free(h+hmin);
    return 0;
  }

}


double *GetHistogram( double *freq1, double *freq2, int **matrix, int *hmin, int *hmax, double *mean ) {
  
  int i, j, n;
  double *h;
  *hmin = INT_MAX;
  *hmax = INT_MIN;
  *mean = 0.0;

  for(i='A';i<='Z';i++) {
    for(j='A';j<='Z';j++) {
      if ( (n=matrix[i][j]) > *hmax ) 
	*hmax = n;
      if ( (n=matrix[i][j]) < *hmin ) 
	*hmin = n;
    }
  }
  
  h = (double*)calloc(*hmax-*hmin+1,sizeof(double))-*hmin;
  
  for(i='A';i<='Z';i++) 
    for(j='A';j<='Z';j++) {
      h[matrix[i][j]]+=freq1[i]*freq2[j];
    }


  for(j=*hmin;j<=*hmax;j++)
    *mean += h[j]*j;
  
  return h;
}


double * PseudoResidueFrequencies2( char *seq, int len, int pseudo_len ) {

  int n, c;
  double *freq = RobinsonResidueFrequencies2();
  double t=0.0;


  if ( len > 0 ) {

    for(c='A';c<='Z';c++) 
      freq[c] = pseudo_len*freq[c];

    for(n=0;n<len;n++) {
      freq[toupper(seq[n])]++;
    }

    t = 0;
    for(c='A';c<='Z';c++) 
      t += freq[c];

    if ( t > 0 )
      for(c='A';c<='Z';c++) {
	freq[c] /= t;
      }
  }
  
  return freq;
}



  double *
RobinsonResidueFrequencies2 (void) {

/* these are are amino-acid frequencies used by blast, from Robinson & Robinson */

  double *freq = (double*)calloc(256, sizeof(double));
  int c;
  double t=0.0;

  freq['A'] = 35155;
  freq['R'] = 23105;
  freq['N'] = 20212;
  freq['D'] = 24161;
  freq['C'] = 8669;
  freq['Q'] = 19208;
  freq['E'] = 28354;
  freq['G'] = 33229;
  freq['H'] = 9906;
  freq['I'] = 23161;
  freq['L'] = 40625;
  freq['K'] = 25872;
  freq['M'] = 10101;
  freq['F'] = 17367;
  freq['P'] = 23435;
  freq['S'] = 32070;
  freq['T'] = 26311;
  freq['W'] = 5990;
  freq['Y'] = 14488;
  freq['V'] = 29012;

  freq['*'] = 0.00+0;
  freq['X'] = 0.00+0;
  freq['x'] = 0.00+0;
  freq['.'] = 0.0;
  
  t = 0.0;
  for(c='A';c<='Z';c++) 
    t += freq[c];
  
  for(c='A';c<='Z';c++){
    freq[toupper(c)] = freq[c] /= t;
  }
  
  
  return freq;
}
