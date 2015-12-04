

MODULE = Wise2 PACKAGE = Wise2::DnaMatrix

Wise2_DnaMatrix *
hard_link_DnaMatrix(obj)
	Wise2_DnaMatrix * obj
	CODE:
	RETVAL = Wise2_hard_link_DnaMatrix(obj);
	OUTPUT:
	RETVAL



Wise2_DnaMatrix *
alloc()
	CODE:
	RETVAL = Wise2_DnaMatrix_alloc();
	OUTPUT:
	RETVAL




Wise2_DnaMatrix *
new(class)
	char * class
	PPCODE:
	Wise2_DnaMatrix * out;
	out = Wise2_DnaMatrix_alloc();
	ST(0) = sv_newmortal();
	sv_setref_pv(ST(0),class,(void*)out);
	XSRETURN(1);

void
DESTROY(obj)
	Wise2_DnaMatrix * obj
	CODE:
	Wise2_free_DnaMatrix(obj);



MODULE = Wise2 PACKAGE = Wise2::DnaProbMatrix

void
flat_null_DnaProbMatrix(dpm)
	Wise2_DnaProbMatrix * dpm
	CODE:
	Wise2_flat_null_DnaProbMatrix(dpm);



Wise2_DnaProbMatrix *
hard_link_DnaProbMatrix(obj)
	Wise2_DnaProbMatrix * obj
	CODE:
	RETVAL = Wise2_hard_link_DnaProbMatrix(obj);
	OUTPUT:
	RETVAL



Wise2_DnaProbMatrix *
alloc()
	CODE:
	RETVAL = Wise2_DnaProbMatrix_alloc();
	OUTPUT:
	RETVAL




Wise2_DnaProbMatrix *
new(class)
	char * class
	PPCODE:
	Wise2_DnaProbMatrix * out;
	out = Wise2_DnaProbMatrix_alloc();
	ST(0) = sv_newmortal();
	sv_setref_pv(ST(0),class,(void*)out);
	XSRETURN(1);

void
DESTROY(obj)
	Wise2_DnaProbMatrix * obj
	CODE:
	Wise2_free_DnaProbMatrix(obj);



MODULE = Wise2 PACKAGE = Wise2

Wise2_DnaMatrix *
identity_DnaMatrix(id_score,mismatch)
	Score id_score
	Score mismatch
	CODE:
	RETVAL = Wise2_identity_DnaMatrix(id_score,mismatch);
	OUTPUT:
	RETVAL



Wise2_DnaProbMatrix *
DnaProbMatrix_from_match(match,nmask_type)
	Probability match
	int nmask_type
	CODE:
	RETVAL = Wise2_DnaProbMatrix_from_match(match,nmask_type);
	OUTPUT:
	RETVAL



Wise2_DnaMatrix *
DnaMatrix_from_DnaProbMatrix(dpm)
	Wise2_DnaProbMatrix * dpm
	CODE:
	RETVAL = Wise2_DnaMatrix_from_DnaProbMatrix(dpm);
	OUTPUT:
	RETVAL



