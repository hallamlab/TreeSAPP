

MODULE = Wise2 PACKAGE = Wise2::cDNAParser

Wise2_cDNAParser *
hard_link_cDNAParser(obj)
	Wise2_cDNAParser * obj
	CODE:
	RETVAL = Wise2_hard_link_cDNAParser(obj);
	OUTPUT:
	RETVAL



Wise2_cDNAParser *
alloc()
	CODE:
	RETVAL = Wise2_cDNAParser_alloc();
	OUTPUT:
	RETVAL




Wise2_cDNAParser *
new(class)
	char * class
	PPCODE:
	Wise2_cDNAParser * out;
	out = Wise2_cDNAParser_alloc();
	ST(0) = sv_newmortal();
	sv_setref_pv(ST(0),class,(void*)out);
	XSRETURN(1);

void
DESTROY(obj)
	Wise2_cDNAParser * obj
	CODE:
	Wise2_free_cDNAParser(obj);



MODULE = Wise2 PACKAGE = Wise2

Wise2_cDNAParser *
flat_cDNAParser(p)
	Probability p
	CODE:
	RETVAL = Wise2_flat_cDNAParser(p);
	OUTPUT:
	RETVAL



