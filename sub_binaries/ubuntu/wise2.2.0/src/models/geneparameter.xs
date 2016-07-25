

MODULE = Wise2 PACKAGE = Wise2::GeneParameter21

Wise2_GeneParameter21 *
hard_link_GeneParameter21(obj)
	Wise2_GeneParameter21 * obj
	CODE:
	RETVAL = Wise2_hard_link_GeneParameter21(obj);
	OUTPUT:
	RETVAL



Wise2_GeneParameter21 *
GeneParameter21_alloc_std()
	CODE:
	RETVAL = Wise2_GeneParameter21_alloc_std();
	OUTPUT:
	RETVAL




Wise2_GeneParameter21 *
new(class)
	char * class
	PPCODE:
	Wise2_GeneParameter21 * out;
	out = Wise2_GeneParameter21_alloc_std();
	ST(0) = sv_newmortal();
	sv_setref_pv(ST(0),class,(void*)out);
	XSRETURN(1);

void
DESTROY(obj)
	Wise2_GeneParameter21 * obj
	CODE:
	Wise2_free_GeneParameter21(obj);



MODULE = Wise2 PACKAGE = Wise2

