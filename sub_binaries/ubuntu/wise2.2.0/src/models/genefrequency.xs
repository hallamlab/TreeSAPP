

MODULE = Wise2 PACKAGE = Wise2::GeneFrequency21

Wise2_GeneFrequency21 *
hard_link_GeneFrequency21(obj)
	Wise2_GeneFrequency21 * obj
	CODE:
	RETVAL = Wise2_hard_link_GeneFrequency21(obj);
	OUTPUT:
	RETVAL



Wise2_GeneFrequency21 *
alloc()
	CODE:
	RETVAL = Wise2_GeneFrequency21_alloc();
	OUTPUT:
	RETVAL



boolean
set_ss5(obj,ss5)
	Wise2_GeneFrequency21 * obj
	Wise2_GeneConsensus * ss5
	CODE:
	RETVAL = Wise2_replace_ss5_GeneFrequency21(obj,Wise2_hard_link_GeneConsensus(ss5));
	OUTPUT:
	RETVAL



Wise2_GeneConsensus *
ss5(obj)
	Wise2_GeneFrequency21 * obj
	INIT:
Wise2_GeneConsensus * temp;
	CODE:
	temp = Wise2_hard_link_GeneConsensus(Wise2_access_ss5_GeneFrequency21(obj));
	RETVAL = temp;
	OUTPUT:
	RETVAL



boolean
set_ss3(obj,ss3)
	Wise2_GeneFrequency21 * obj
	Wise2_GeneConsensus * ss3
	CODE:
	RETVAL = Wise2_replace_ss3_GeneFrequency21(obj,Wise2_hard_link_GeneConsensus(ss3));
	OUTPUT:
	RETVAL



Wise2_GeneConsensus *
ss3(obj)
	Wise2_GeneFrequency21 * obj
	INIT:
Wise2_GeneConsensus * temp;
	CODE:
	temp = Wise2_hard_link_GeneConsensus(Wise2_access_ss3_GeneFrequency21(obj));
	RETVAL = temp;
	OUTPUT:
	RETVAL




Wise2_GeneFrequency21 *
new(class)
	char * class
	PPCODE:
	Wise2_GeneFrequency21 * out;
	out = Wise2_GeneFrequency21_alloc();
	ST(0) = sv_newmortal();
	sv_setref_pv(ST(0),class,(void*)out);
	XSRETURN(1);

void
DESTROY(obj)
	Wise2_GeneFrequency21 * obj
	CODE:
	Wise2_free_GeneFrequency21(obj);



MODULE = Wise2 PACKAGE = Wise2::GeneConsensus

Wise2_GeneConsensus *
hard_link_GeneConsensus(obj)
	Wise2_GeneConsensus * obj
	CODE:
	RETVAL = Wise2_hard_link_GeneConsensus(obj);
	OUTPUT:
	RETVAL



Wise2_GeneConsensus *
GeneConsensus_alloc_std()
	CODE:
	RETVAL = Wise2_GeneConsensus_alloc_std();
	OUTPUT:
	RETVAL



boolean
set_center(obj,center)
	Wise2_GeneConsensus * obj
	int center
	CODE:
	RETVAL = Wise2_replace_center_GeneConsensus(obj,center);
	OUTPUT:
	RETVAL



int
center(obj)
	Wise2_GeneConsensus * obj
	CODE:
	RETVAL = Wise2_access_center_GeneConsensus(obj);
	OUTPUT:
	RETVAL



Wise2_GeneSingleCons *
gsc(obj,i)
	Wise2_GeneConsensus * obj
	int i
	INIT:
Wise2_GeneSingleCons * temp;
	CODE:
	temp = Wise2_hard_link_GeneSingleCons(Wise2_access_gsc_GeneConsensus(obj,i));
	RETVAL = temp;
	OUTPUT:
	RETVAL



int
length_gsc(obj)
	Wise2_GeneConsensus * obj
	CODE:
	RETVAL = Wise2_length_gsc_GeneConsensus(obj);
	OUTPUT:
	RETVAL



int
flush_gsc(obj)
	Wise2_GeneConsensus * obj
	CODE:
	RETVAL = Wise2_flush_GeneConsensus(obj);
	OUTPUT:
	RETVAL



boolean
add_gsc(obj,add)
	Wise2_GeneConsensus * obj
	Wise2_GeneSingleCons * add
	CODE:
	RETVAL = Wise2_add_GeneConsensus(obj,Wise2_hard_link_GeneSingleCons(add));
	OUTPUT:
	RETVAL




Wise2_GeneConsensus *
new(class)
	char * class
	PPCODE:
	Wise2_GeneConsensus * out;
	out = Wise2_GeneConsensus_alloc_std();
	ST(0) = sv_newmortal();
	sv_setref_pv(ST(0),class,(void*)out);
	XSRETURN(1);

void
DESTROY(obj)
	Wise2_GeneConsensus * obj
	CODE:
	Wise2_free_GeneConsensus(obj);

void
each_gsc(obj)
	Wise2_GeneConsensus * obj
	PPCODE:
	int i=0;
	int len;
	SV* temp;
	len = Wise2_length_gsc_GeneConsensus(obj);
	for(i=0;i<len;i++){
	  temp = sv_newmortal();
	  sv_setref_pv(temp, "Wise2::GeneSingleCons", (void*) (Wise2_hard_link_GeneSingleCons(Wise2_access_gsc_GeneConsensus(obj,i))));
	  XPUSHs(temp);
	  }
	XSRETURN(len);



MODULE = Wise2 PACKAGE = Wise2::GeneSingleCons

Wise2_GeneSingleCons *
hard_link_GeneSingleCons(obj)
	Wise2_GeneSingleCons * obj
	CODE:
	RETVAL = Wise2_hard_link_GeneSingleCons(obj);
	OUTPUT:
	RETVAL



Wise2_GeneSingleCons *
alloc()
	CODE:
	RETVAL = Wise2_GeneSingleCons_alloc();
	OUTPUT:
	RETVAL



boolean
set_string(obj,string)
	Wise2_GeneSingleCons * obj
	char * string
	CODE:
	RETVAL = Wise2_replace_string_GeneSingleCons(obj,Wise2_stringalloc(string));
	OUTPUT:
	RETVAL



char *
string(obj)
	Wise2_GeneSingleCons * obj
	INIT:
	char * temp;
	CODE:
	temp = Wise2_stringalloc(Wise2_access_string_GeneSingleCons(obj));
	RETVAL = temp;
	OUTPUT:
	RETVAL



boolean
set_number(obj,number)
	Wise2_GeneSingleCons * obj
	double number
	CODE:
	RETVAL = Wise2_replace_number_GeneSingleCons(obj,number);
	OUTPUT:
	RETVAL



double
number(obj)
	Wise2_GeneSingleCons * obj
	CODE:
	RETVAL = Wise2_access_number_GeneSingleCons(obj);
	OUTPUT:
	RETVAL




Wise2_GeneSingleCons *
new(class)
	char * class
	PPCODE:
	Wise2_GeneSingleCons * out;
	out = Wise2_GeneSingleCons_alloc();
	ST(0) = sv_newmortal();
	sv_setref_pv(ST(0),class,(void*)out);
	XSRETURN(1);

void
DESTROY(obj)
	Wise2_GeneSingleCons * obj
	CODE:
	Wise2_free_GeneSingleCons(obj);



MODULE = Wise2 PACKAGE = Wise2

Wise2_GeneFrequency21 *
read_GeneFrequency21_file(filename)
	char * filename
	CODE:
	RETVAL = Wise2_read_GeneFrequency21_file(filename);
	OUTPUT:
	RETVAL



Wise2_GeneFrequency21 *
read_GeneFrequency21(ifp)
	FILE * ifp
	CODE:
	RETVAL = Wise2_read_GeneFrequency21(ifp);
	OUTPUT:
	RETVAL



