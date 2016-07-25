

MODULE = Wise2 PACKAGE = Wise2::Histogram

void
UnfitHistogram(h)
	Wise2_Histogram * h
	CODE:
	Wise2_UnfitHistogram(h);



void
add(h,sc)
	Wise2_Histogram * h
	float sc
	CODE:
	Wise2_AddToHistogram(h,sc);



void
show(h,fp)
	Wise2_Histogram * h
	FILE * fp
	CODE:
	Wise2_PrintASCIIHistogram(h,fp);



void
EVDBasicFit(h)
	Wise2_Histogram * h
	CODE:
	Wise2_EVDBasicFit(h);



int
fit_EVD(h,censor,high_hint)
	Wise2_Histogram * h
	int censor
	float high_hint
	CODE:
	RETVAL = Wise2_ExtremeValueFitHistogram(h,censor,high_hint);
	OUTPUT:
	RETVAL



void
set_EVD(h,mu,lambda,lowbound,highbound,wonka,ndegrees)
	Wise2_Histogram * h
	float mu
	float lambda
	float lowbound
	float highbound
	float wonka
	int ndegrees
	CODE:
	Wise2_ExtremeValueSetHistogram(h,mu,lambda,lowbound,highbound,wonka,ndegrees);



int
fit_Gaussian(h,high_hint)
	Wise2_Histogram * h
	float high_hint
	CODE:
	RETVAL = Wise2_GaussianFitHistogram(h,high_hint);
	OUTPUT:
	RETVAL



void
set_Gaussian(h,mean,sd)
	Wise2_Histogram * h
	float mean
	float sd
	CODE:
	Wise2_GaussianSetHistogram(h,mean,sd);



double
evalue(his,score)
	Wise2_Histogram * his
	double score
	CODE:
	RETVAL = Wise2_Evalue_from_Histogram(his,score);
	OUTPUT:
	RETVAL



Wise2_Histogram *
hard_link_Histogram(obj)
	Wise2_Histogram * obj
	CODE:
	RETVAL = Wise2_hard_link_Histogram(obj);
	OUTPUT:
	RETVAL



Wise2_Histogram *
alloc()
	CODE:
	RETVAL = Wise2_Histogram_alloc();
	OUTPUT:
	RETVAL




Wise2_Histogram *
new(class)
	char * class
	PPCODE:
	Wise2_Histogram * out;
	out = Wise2_Histogram_alloc();
	ST(0) = sv_newmortal();
	sv_setref_pv(ST(0),class,(void*)out);
	XSRETURN(1);

void
DESTROY(obj)
	Wise2_Histogram * obj
	CODE:
	Wise2_free_Histogram(obj);



MODULE = Wise2 PACKAGE = Wise2

Wise2_Histogram *
new_Histogram(min,max,lumpsize)
	int min
	int max
	int lumpsize
	CODE:
	RETVAL = Wise2_new_Histogram(min,max,lumpsize);
	OUTPUT:
	RETVAL



