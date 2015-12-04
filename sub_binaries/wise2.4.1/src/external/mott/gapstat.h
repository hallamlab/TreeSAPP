void init_lambda( double lambda );
double e2lk( int k);


int 
KarlinAltschulStatistics( int **matrix, double *freq1, double *freq2, double *lambda, double *Kminus, double *Kplus, double *H, double *r, double *s ) ;

double
lambda_func( double lambda, double *h, int hmin, int hmax);

double
solve_for_lambda( double *h, int hmin, int hmax  );

double entropy_H( double *h, int hmin, int hmax, double lambda ) ;

double HSP_length_correction( double H, double K, double lambda, int len1, int len2);

void iglehart( double *h, int hmin, int hmax, double lambda, double xmu, int delta, double *R1, double *R2, double *R3, double *Kplus, double *Kminus );

double *get_h( int **matrix, double *freq1, double *freq2, int *hmin, int *hmax, double *mean );

double associated_dam_eqn( double *distribution, int min, int max, double lambda, double 
**transient, int *cn);

int gem_statistics( double gap_start, double gap_extend, double lambda, double s, double Kplus, double Kminus, double *alpha, double *theta1, double *theta2, double *K1, double *K2 );

double upper_objective( double alpha, double beta );
double lower_objective( double alpha, double beta );
double gapped_theta( double alpha, double (*objective)( double, double) );

double GEM_length_correction( double L_HSP, double alpha, int len1, int len2 );


double regression_coeffs( double lambda, double Kmin, double H, double alpha, int len1, int len2, double *theta, double *kappa, double score );

double empirical_GEM( double lambda0, double K0, double H, double alpha, int len1, int len2, double *theta, double *logkappa, double score );

double EmpiricalGEM( double lambda0, double K0, double H, double alpha, int len1, int len2, double *theta, double *logkappa, double score );

double AutoIncreaseGapPenalty( double s, double lambda0, double alpha_max, int *A, int *B );
