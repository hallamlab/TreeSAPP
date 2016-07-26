

MODULE = Wise2 PACKAGE = Wise2::ThreeStateDB

Wise2_ThreeStateModel *
indexed_model(mdb,en)
	Wise2_ThreeStateDB * mdb
	Wise2_DataEntry * en
	CODE:
	RETVAL = Wise2_indexed_ThreeStateModel_ThreeStateDB(mdb,en);
	OUTPUT:
	RETVAL



Wise2_ThreeStateDB *
new_proteindb_ThreeStateDB(sdb,comp,gap,ext)
	Wise2_SequenceDB * sdb
	Wise2_CompMat * comp
	int gap
	int ext
	CODE:
	RETVAL = Wise2_new_proteindb_ThreeStateDB(sdb,comp,gap,ext);
	OUTPUT:
	RETVAL



Wise2_ThreeStateDB *
new_PfamHmmer1DB_ThreeStateDB(dirname)
	char * dirname
	CODE:
	RETVAL = Wise2_new_PfamHmmer1DB_ThreeStateDB(dirname);
	OUTPUT:
	RETVAL



Wise2_ThreeStateDB *
new_single_ThreeStateDB(tsm,rm)
	Wise2_ThreeStateModel * tsm
	Wise2_RandomModel * rm
	CODE:
	RETVAL = Wise2_new_single_ThreeStateDB(tsm,rm);
	OUTPUT:
	RETVAL



Wise2_ThreeStateDB *
hard_link_ThreeStateDB(obj)
	Wise2_ThreeStateDB * obj
	CODE:
	RETVAL = Wise2_hard_link_ThreeStateDB(obj);
	OUTPUT:
	RETVAL



Wise2_ThreeStateDB *
alloc()
	CODE:
	RETVAL = Wise2_ThreeStateDB_alloc();
	OUTPUT:
	RETVAL



boolean
set_dbtype(obj,dbtype)
	Wise2_ThreeStateDB * obj
	int dbtype
	CODE:
	RETVAL = Wise2_replace_dbtype_ThreeStateDB(obj,dbtype);
	OUTPUT:
	RETVAL



int
dbtype(obj)
	Wise2_ThreeStateDB * obj
	CODE:
	RETVAL = Wise2_access_dbtype_ThreeStateDB(obj);
	OUTPUT:
	RETVAL



boolean
set_filename(obj,filename)
	Wise2_ThreeStateDB * obj
	char * filename
	CODE:
	RETVAL = Wise2_replace_filename_ThreeStateDB(obj,Wise2_stringalloc(filename));
	OUTPUT:
	RETVAL



char *
filename(obj)
	Wise2_ThreeStateDB * obj
	INIT:
	char * temp;
	CODE:
	temp = Wise2_stringalloc(Wise2_access_filename_ThreeStateDB(obj));
	RETVAL = temp;
	OUTPUT:
	RETVAL




Wise2_ThreeStateDB *
new(class)
	char * class
	PPCODE:
	Wise2_ThreeStateDB * out;
	out = Wise2_ThreeStateDB_alloc();
	ST(0) = sv_newmortal();
	sv_setref_pv(ST(0),class,(void*)out);
	XSRETURN(1);

void
DESTROY(obj)
	Wise2_ThreeStateDB * obj
	CODE:
	Wise2_free_ThreeStateDB(obj);



MODULE = Wise2 PACKAGE = Wise2

