

MODULE = Wise2 PACKAGE = Wise2::RandomProteinDB

Wise2_RandomProteinDB *
hard_link_RandomProteinDB(obj)
	Wise2_RandomProteinDB * obj
	CODE:
	RETVAL = Wise2_hard_link_RandomProteinDB(obj);
	OUTPUT:
	RETVAL



Wise2_RandomProteinDB *
alloc()
	CODE:
	RETVAL = Wise2_RandomProteinDB_alloc();
	OUTPUT:
	RETVAL



boolean
set_use_flat_length(obj,use_flat_length)
	Wise2_RandomProteinDB * obj
	boolean use_flat_length
	CODE:
	RETVAL = Wise2_replace_use_flat_length_RandomProteinDB(obj,use_flat_length);
	OUTPUT:
	RETVAL



boolean
use_flat_length(obj)
	Wise2_RandomProteinDB * obj
	CODE:
	RETVAL = Wise2_access_use_flat_length_RandomProteinDB(obj);
	OUTPUT:
	RETVAL



boolean
set_length(obj,length)
	Wise2_RandomProteinDB * obj
	int length
	CODE:
	RETVAL = Wise2_replace_length_RandomProteinDB(obj,length);
	OUTPUT:
	RETVAL



int
length(obj)
	Wise2_RandomProteinDB * obj
	CODE:
	RETVAL = Wise2_access_length_RandomProteinDB(obj);
	OUTPUT:
	RETVAL



boolean
set_length_dist(obj,length_dist)
	Wise2_RandomProteinDB * obj
	Wise2_Histogram * length_dist
	CODE:
	RETVAL = Wise2_replace_length_dist_RandomProteinDB(obj,Wise2_hard_link_Histogram(length_dist));
	OUTPUT:
	RETVAL



Wise2_Histogram *
length_dist(obj)
	Wise2_RandomProteinDB * obj
	INIT:
Wise2_Histogram * temp;
	CODE:
	temp = Wise2_hard_link_Histogram(Wise2_access_length_dist_RandomProteinDB(obj));
	RETVAL = temp;
	OUTPUT:
	RETVAL



boolean
set_emission(obj,emission)
	Wise2_RandomProteinDB * obj
	Wise2_RandomModel * emission
	CODE:
	RETVAL = Wise2_replace_emission_RandomProteinDB(obj,Wise2_hard_link_RandomModel(emission));
	OUTPUT:
	RETVAL



Wise2_RandomModel *
emission(obj)
	Wise2_RandomProteinDB * obj
	INIT:
Wise2_RandomModel * temp;
	CODE:
	temp = Wise2_hard_link_RandomModel(Wise2_access_emission_RandomProteinDB(obj));
	RETVAL = temp;
	OUTPUT:
	RETVAL



boolean
set_num(obj,num)
	Wise2_RandomProteinDB * obj
	int num
	CODE:
	RETVAL = Wise2_replace_num_RandomProteinDB(obj,num);
	OUTPUT:
	RETVAL



int
num(obj)
	Wise2_RandomProteinDB * obj
	CODE:
	RETVAL = Wise2_access_num_RandomProteinDB(obj);
	OUTPUT:
	RETVAL




Wise2_RandomProteinDB *
new(class)
	char * class
	PPCODE:
	Wise2_RandomProteinDB * out;
	out = Wise2_RandomProteinDB_alloc();
	ST(0) = sv_newmortal();
	sv_setref_pv(ST(0),class,(void*)out);
	XSRETURN(1);

void
DESTROY(obj)
	Wise2_RandomProteinDB * obj
	CODE:
	Wise2_free_RandomProteinDB(obj);



MODULE = Wise2 PACKAGE = Wise2::RandomDNADB

Wise2_RandomDNADB *
hard_link_RandomDNADB(obj)
	Wise2_RandomDNADB * obj
	CODE:
	RETVAL = Wise2_hard_link_RandomDNADB(obj);
	OUTPUT:
	RETVAL



Wise2_RandomDNADB *
alloc()
	CODE:
	RETVAL = Wise2_RandomDNADB_alloc();
	OUTPUT:
	RETVAL



boolean
set_use_flat_length(obj,use_flat_length)
	Wise2_RandomDNADB * obj
	boolean use_flat_length
	CODE:
	RETVAL = Wise2_replace_use_flat_length_RandomDNADB(obj,use_flat_length);
	OUTPUT:
	RETVAL



boolean
use_flat_length(obj)
	Wise2_RandomDNADB * obj
	CODE:
	RETVAL = Wise2_access_use_flat_length_RandomDNADB(obj);
	OUTPUT:
	RETVAL



boolean
set_length(obj,length)
	Wise2_RandomDNADB * obj
	int length
	CODE:
	RETVAL = Wise2_replace_length_RandomDNADB(obj,length);
	OUTPUT:
	RETVAL



int
length(obj)
	Wise2_RandomDNADB * obj
	CODE:
	RETVAL = Wise2_access_length_RandomDNADB(obj);
	OUTPUT:
	RETVAL



boolean
set_length_dist(obj,length_dist)
	Wise2_RandomDNADB * obj
	Wise2_Histogram * length_dist
	CODE:
	RETVAL = Wise2_replace_length_dist_RandomDNADB(obj,Wise2_hard_link_Histogram(length_dist));
	OUTPUT:
	RETVAL



Wise2_Histogram *
length_dist(obj)
	Wise2_RandomDNADB * obj
	INIT:
Wise2_Histogram * temp;
	CODE:
	temp = Wise2_hard_link_Histogram(Wise2_access_length_dist_RandomDNADB(obj));
	RETVAL = temp;
	OUTPUT:
	RETVAL



boolean
set_emission(obj,emission)
	Wise2_RandomDNADB * obj
	Wise2_RandomModelDNA * emission
	CODE:
	RETVAL = Wise2_replace_emission_RandomDNADB(obj,Wise2_hard_link_RandomModelDNA(emission));
	OUTPUT:
	RETVAL



Wise2_RandomModelDNA *
emission(obj)
	Wise2_RandomDNADB * obj
	INIT:
Wise2_RandomModelDNA * temp;
	CODE:
	temp = Wise2_hard_link_RandomModelDNA(Wise2_access_emission_RandomDNADB(obj));
	RETVAL = temp;
	OUTPUT:
	RETVAL



boolean
set_num(obj,num)
	Wise2_RandomDNADB * obj
	int num
	CODE:
	RETVAL = Wise2_replace_num_RandomDNADB(obj,num);
	OUTPUT:
	RETVAL



int
num(obj)
	Wise2_RandomDNADB * obj
	CODE:
	RETVAL = Wise2_access_num_RandomDNADB(obj);
	OUTPUT:
	RETVAL




Wise2_RandomDNADB *
new(class)
	char * class
	PPCODE:
	Wise2_RandomDNADB * out;
	out = Wise2_RandomDNADB_alloc();
	ST(0) = sv_newmortal();
	sv_setref_pv(ST(0),class,(void*)out);
	XSRETURN(1);

void
DESTROY(obj)
	Wise2_RandomDNADB * obj
	CODE:
	Wise2_free_RandomDNADB(obj);



MODULE = Wise2 PACKAGE = Wise2

