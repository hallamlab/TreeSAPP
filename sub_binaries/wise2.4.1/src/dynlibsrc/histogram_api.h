

/* Functions that create, manipulate or act on Histogram
 *
 * Wise2_UnfitHistogram
 * Wise2_AddToHistogram
 * Wise2_PrintASCIIHistogram
 * Wise2_EVDBasicFit
 * Wise2_ExtremeValueFitHistogram
 * Wise2_ExtremeValueSetHistogram
 * Wise2_GaussianFitHistogram
 * Wise2_GaussianSetHistogram
 * Wise2_Evalue_from_Histogram
 * Wise2_hard_link_Histogram
 * Wise2_Histogram_alloc
 * Wise2_free_Histogram [destructor]
 *
 */



/* Helper functions in the module
 *
 * Wise2_new_Histogram
 *

/* API for object Histogram */
/* Function:  Wise2_UnfitHistogram(h)
 *
 * Descrip: No Description
 *
 * Arg:        h            Undocumented argument [Wise2_Histogram *]
 *
 * Returns Undocumented return value [void]
 *
 */
void Wise2_UnfitHistogram( Wise2_Histogram * h);

/* Function:  Wise2_AddToHistogram(h,sc)
 *
 * Descrip: No Description
 *
 * Arg:        h            Undocumented argument [Wise2_Histogram *]
 * Arg:        sc           Undocumented argument [float]
 *
 * Returns Undocumented return value [void]
 *
 */
void Wise2_AddToHistogram( Wise2_Histogram * h,float sc);

/* Function:  Wise2_PrintASCIIHistogram(h,fp)
 *
 * Descrip: No Description
 *
 * Arg:        h            histogram to print [Wise2_Histogram *]
 * Arg:        fp           open file to print to (stdout works) [FILE *]
 *
 * Returns Undocumented return value [void]
 *
 */
void Wise2_PrintASCIIHistogram( Wise2_Histogram * h,FILE * fp);

/* Function:  Wise2_EVDBasicFit(h)
 *
 * Descrip: No Description
 *
 * Arg:        h            histogram to fit [Wise2_Histogram *]
 *
 * Returns Undocumented return value [void]
 *
 */
void Wise2_EVDBasicFit( Wise2_Histogram * h);

/* Function:  Wise2_ExtremeValueFitHistogram(h,censor,high_hint)
 *
 * Descrip: No Description
 *
 * Arg:        h            histogram to fit [Wise2_Histogram *]
 * Arg:        censor       TRUE to censor data left of the peak [int]
 * Arg:        high_hint    score cutoff; above this are real hits that arent fit [float]
 *
 * Returns if fit is judged to be valid else 0 if fit is invalid (too few seqs.) [int]
 *
 */
int Wise2_ExtremeValueFitHistogram( Wise2_Histogram * h,int censor,float high_hint);

/* Function:  Wise2_ExtremeValueSetHistogram(h,mu,lambda,lowbound,highbound,wonka,ndegrees)
 *
 * Descrip: No Description
 *
 * Arg:        h            the histogram to set [Wise2_Histogram *]
 * Arg:        mu           mu location parameter                 [float]
 * Arg:        lambda       lambda scale parameter [float]
 * Arg:        lowbound     low bound of the histogram that was fit [float]
 * Arg:        highbound    high bound of histogram that was fit [float]
 * Arg:        wonka        fudge factor; fraction of hits estimated to be "EVD-like" [float]
 * Arg:        ndegrees     extra degrees of freedom to subtract in chi2 test: [int]
 *
 * Returns Undocumented return value [void]
 *
 */
void Wise2_ExtremeValueSetHistogram( Wise2_Histogram * h,float mu,float lambda,float lowbound,float highbound,float wonka,int ndegrees);

/* Function:  Wise2_GaussianFitHistogram(h,high_hint)
 *
 * Descrip: No Description
 *
 * Arg:        h            histogram to fit [Wise2_Histogram *]
 * Arg:        high_hint    score cutoff; above this are `real' hits that aren't fit [float]
 *
 * Returns if fit is judged to be valid else 0 if fit is invalid (too few seqs.)            [int]
 *
 */
int Wise2_GaussianFitHistogram( Wise2_Histogram * h,float high_hint);

/* Function:  Wise2_GaussianSetHistogram(h,mean,sd)
 *
 * Descrip: No Description
 *
 * Arg:        h            Undocumented argument [Wise2_Histogram *]
 * Arg:        mean         Undocumented argument [float]
 * Arg:        sd           Undocumented argument [float]
 *
 * Returns Undocumented return value [void]
 *
 */
void Wise2_GaussianSetHistogram( Wise2_Histogram * h,float mean,float sd);

/* Function:  Wise2_Evalue_from_Histogram(his,score)
 *
 * Descrip: No Description
 *
 * Arg:        his          Histogram object [Wise2_Histogram *]
 * Arg:        score        score you want the evalue for [double]
 *
 * Returns Undocumented return value [double]
 *
 */
double Wise2_Evalue_from_Histogram( Wise2_Histogram * his,double score);

/* Function:  Wise2_hard_link_Histogram(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj          Object to be hard linked [Wise2_Histogram *]
 *
 * Returns Undocumented return value [Wise2_Histogram *]
 *
 */
Wise2_Histogram * Wise2_hard_link_Histogram( Wise2_Histogram * obj);

/* Function:  Wise2_Histogram_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Returns Undocumented return value [Wise2_Histogram *]
 *
 */
Wise2_Histogram * Wise2_Histogram_alloc();

/* This is the destructor function, ie, call this to free object*/
/* Function:  Wise2_free_Histogram(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj          Object that is free'd [Wise2_Histogram *]
 *
 * Returns Undocumented return value [Wise2_Histogram *]
 *
 */
Wise2_Histogram * Wise2_free_Histogram( Wise2_Histogram * obj);



/* These functions are not associated with an object */
/* Function:  Wise2_new_Histogram(min,max,lumpsize)
 *
 * Descrip: No Description
 *
 * Arg:        min          minimum score (integer) [int]
 * Arg:        max          maximum score (integer) [int]
 * Arg:        lumpsize     when reallocating histogram, the reallocation amount [int]
 *
 * Returns Undocumented return value [Wise2_Histogram *]
 *
 */
Wise2_Histogram * Wise2_new_Histogram( int min,int max,int lumpsize);

