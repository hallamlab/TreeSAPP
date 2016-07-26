

MODULE = Wise2 PACKAGE = Wise2::RandomModelDNA

Wise2_RandomModelDNA *
hard_link_RandomModelDNA(obj)
	Wise2_RandomModelDNA * obj
	CODE:
	RETVAL = Wise2_hard_link_RandomModelDNA(obj);
	OUTPUT:
	RETVAL



Wise2_RandomModelDNA *
alloc()
	CODE:
	RETVAL = Wise2_RandomModelDNA_alloc();
	OUTPUT:
	RETVAL



boolean
set_name(obj,name)
	Wise2_RandomModelDNA * obj
	char * name
	CODE:
	RETVAL = Wise2_replace_name_RandomModelDNA(obj,Wise2_stringalloc(name));
	OUTPUT:
	RETVAL



char *
name(obj)
	Wise2_RandomModelDNA * obj
	INIT:
	char * temp;
	CODE:
	temp = Wise2_stringalloc(Wise2_access_name_RandomModelDNA(obj));
	RETVAL = temp;
	OUTPUT:
	RETVAL




Wise2_RandomModelDNA *
new(class)
	char * class
	PPCODE:
	Wise2_RandomModelDNA * out;
	out = Wise2_RandomModelDNA_alloc();
	ST(0) = sv_newmortal();
	sv_setref_pv(ST(0),class,(void*)out);
	XSRETURN(1);

void
DESTROY(obj)
	Wise2_RandomModelDNA * obj
	CODE:
	Wise2_free_RandomModelDNA(obj);



MODULE = Wise2 PACKAGE = Wise2::RandomModel

Wise2_RandomModel *
hard_link_RandomModel(obj)
	Wise2_RandomModel * obj
	CODE:
	RETVAL = Wise2_hard_link_RandomModel(obj);
	OUTPUT:
	RETVAL



Wise2_RandomModel *
alloc()
	CODE:
	RETVAL = Wise2_RandomModel_alloc();
	OUTPUT:
	RETVAL



boolean
set_name(obj,name)
	Wise2_RandomModel * obj
	char * name
	CODE:
	RETVAL = Wise2_replace_name_RandomModel(obj,Wise2_stringalloc(name));
	OUTPUT:
	RETVAL



char *
name(obj)
	Wise2_RandomModel * obj
	INIT:
	char * temp;
	CODE:
	temp = Wise2_stringalloc(Wise2_access_name_RandomModel(obj));
	RETVAL = temp;
	OUTPUT:
	RETVAL




Wise2_RandomModel *
new(class)
	char * class
	PPCODE:
	Wise2_RandomModel * out;
	out = Wise2_RandomModel_alloc();
	ST(0) = sv_newmortal();
	sv_setref_pv(ST(0),class,(void*)out);
	XSRETURN(1);

void
DESTROY(obj)
	Wise2_RandomModel * obj
	CODE:
	Wise2_free_RandomModel(obj);



MODULE = Wise2 PACKAGE = Wise2

Wise2_RandomModelDNA *
RandomModelDNA_std()
	CODE:
	RETVAL = Wise2_RandomModelDNA_std();
	OUTPUT:
	RETVAL



Wise2_RandomModel *
default_RandomModel()
	CODE:
	RETVAL = Wise2_default_RandomModel();
	OUTPUT:
	RETVAL



