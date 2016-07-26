

MODULE = Wise2 PACKAGE = Wise2::CodonMapper

void
sprinkle_errors_over_CodonMapper(cm,error)
	Wise2_CodonMapper * cm
	double error
	CODE:
	Wise2_sprinkle_errors_over_CodonMapper(cm,error);



Wise2_CodonMapper *
hard_link_CodonMapper(obj)
	Wise2_CodonMapper * obj
	CODE:
	RETVAL = Wise2_hard_link_CodonMapper(obj);
	OUTPUT:
	RETVAL



Wise2_CodonMapper *
alloc()
	CODE:
	RETVAL = Wise2_CodonMapper_alloc();
	OUTPUT:
	RETVAL



boolean
set_ct(obj,ct)
	Wise2_CodonMapper * obj
	Wise2_CodonTable * ct
	CODE:
	RETVAL = Wise2_replace_ct_CodonMapper(obj,Wise2_hard_link_CodonTable(ct));
	OUTPUT:
	RETVAL



Wise2_CodonTable *
ct(obj)
	Wise2_CodonMapper * obj
	INIT:
Wise2_CodonTable * temp;
	CODE:
	temp = Wise2_hard_link_CodonTable(Wise2_access_ct_CodonMapper(obj));
	RETVAL = temp;
	OUTPUT:
	RETVAL




Wise2_CodonMapper *
new(class)
	char * class
	PPCODE:
	Wise2_CodonMapper * out;
	out = Wise2_CodonMapper_alloc();
	ST(0) = sv_newmortal();
	sv_setref_pv(ST(0),class,(void*)out);
	XSRETURN(1);

void
DESTROY(obj)
	Wise2_CodonMapper * obj
	CODE:
	Wise2_free_CodonMapper(obj);



MODULE = Wise2 PACKAGE = Wise2

Wise2_CodonMapper *
flat_CodonMapper(ct)
	Wise2_CodonTable * ct
	CODE:
	RETVAL = Wise2_flat_CodonMapper(ct);
	OUTPUT:
	RETVAL



