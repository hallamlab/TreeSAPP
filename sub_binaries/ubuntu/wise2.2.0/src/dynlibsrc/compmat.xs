

MODULE = Wise2 PACKAGE = Wise2::CompMat

Score
fail_safe_CompMat_access(cm,aa1,aa2)
	Wise2_CompMat * cm
	int aa1
	int aa2
	CODE:
	RETVAL = Wise2_fail_safe_CompMat_access(cm,aa1,aa2);
	OUTPUT:
	RETVAL



boolean
write_Blast_CompMat(cm,ofp)
	Wise2_CompMat * cm
	FILE * ofp
	CODE:
	RETVAL = Wise2_write_Blast_CompMat(cm,ofp);
	OUTPUT:
	RETVAL



Wise2_CompMat *
read_Blast_file_CompMat(filename)
	char * filename
	CODE:
	RETVAL = Wise2_read_Blast_file_CompMat(filename);
	OUTPUT:
	RETVAL



Wise2_CompMat *
read_Blast_CompMat(ifp)
	FILE * ifp
	CODE:
	RETVAL = Wise2_read_Blast_CompMat(ifp);
	OUTPUT:
	RETVAL



Wise2_CompMat *
hard_link_CompMat(obj)
	Wise2_CompMat * obj
	CODE:
	RETVAL = Wise2_hard_link_CompMat(obj);
	OUTPUT:
	RETVAL



Wise2_CompMat *
alloc()
	CODE:
	RETVAL = Wise2_CompMat_alloc();
	OUTPUT:
	RETVAL



boolean
set_name(obj,name)
	Wise2_CompMat * obj
	char * name
	CODE:
	RETVAL = Wise2_replace_name_CompMat(obj,Wise2_stringalloc(name));
	OUTPUT:
	RETVAL



char *
name(obj)
	Wise2_CompMat * obj
	INIT:
	char * temp;
	CODE:
	temp = Wise2_stringalloc(Wise2_access_name_CompMat(obj));
	RETVAL = temp;
	OUTPUT:
	RETVAL




Wise2_CompMat *
new(class)
	char * class
	PPCODE:
	Wise2_CompMat * out;
	out = Wise2_CompMat_alloc();
	ST(0) = sv_newmortal();
	sv_setref_pv(ST(0),class,(void*)out);
	XSRETURN(1);

void
DESTROY(obj)
	Wise2_CompMat * obj
	CODE:
	Wise2_free_CompMat(obj);



MODULE = Wise2 PACKAGE = Wise2

