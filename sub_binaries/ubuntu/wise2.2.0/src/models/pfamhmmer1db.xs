

MODULE = Wise2 PACKAGE = Wise2::PfamHmmer1DB

Wise2_PfamHmmer1DB *
hard_link_PfamHmmer1DB(obj)
	Wise2_PfamHmmer1DB * obj
	CODE:
	RETVAL = Wise2_hard_link_PfamHmmer1DB(obj);
	OUTPUT:
	RETVAL



Wise2_PfamHmmer1DB *
PfamHmmer1DB_alloc_std()
	CODE:
	RETVAL = Wise2_PfamHmmer1DB_alloc_std();
	OUTPUT:
	RETVAL



Wise2_PfamHmmer1Entry *
en(obj,i)
	Wise2_PfamHmmer1DB * obj
	int i
	INIT:
Wise2_PfamHmmer1Entry * temp;
	CODE:
	temp = Wise2_hard_link_PfamHmmer1Entry(Wise2_access_en_PfamHmmer1DB(obj,i));
	RETVAL = temp;
	OUTPUT:
	RETVAL



int
length_en(obj)
	Wise2_PfamHmmer1DB * obj
	CODE:
	RETVAL = Wise2_length_en_PfamHmmer1DB(obj);
	OUTPUT:
	RETVAL



int
flush_en(obj)
	Wise2_PfamHmmer1DB * obj
	CODE:
	RETVAL = Wise2_flush_PfamHmmer1DB(obj);
	OUTPUT:
	RETVAL



boolean
add_en(obj,add)
	Wise2_PfamHmmer1DB * obj
	Wise2_PfamHmmer1Entry * add
	CODE:
	RETVAL = Wise2_add_PfamHmmer1DB(obj,Wise2_hard_link_PfamHmmer1Entry(add));
	OUTPUT:
	RETVAL



boolean
set_dirname(obj,dirname)
	Wise2_PfamHmmer1DB * obj
	char * dirname
	CODE:
	RETVAL = Wise2_replace_dirname_PfamHmmer1DB(obj,Wise2_stringalloc(dirname));
	OUTPUT:
	RETVAL



char *
dirname(obj)
	Wise2_PfamHmmer1DB * obj
	INIT:
	char * temp;
	CODE:
	temp = Wise2_stringalloc(Wise2_access_dirname_PfamHmmer1DB(obj));
	RETVAL = temp;
	OUTPUT:
	RETVAL



boolean
set_cur(obj,cur)
	Wise2_PfamHmmer1DB * obj
	int cur
	CODE:
	RETVAL = Wise2_replace_cur_PfamHmmer1DB(obj,cur);
	OUTPUT:
	RETVAL



int
cur(obj)
	Wise2_PfamHmmer1DB * obj
	CODE:
	RETVAL = Wise2_access_cur_PfamHmmer1DB(obj);
	OUTPUT:
	RETVAL



boolean
set_def(obj,def)
	Wise2_PfamHmmer1DB * obj
	Wise2_RandomModel * def
	CODE:
	RETVAL = Wise2_replace_def_PfamHmmer1DB(obj,Wise2_hard_link_RandomModel(def));
	OUTPUT:
	RETVAL



Wise2_RandomModel *
def(obj)
	Wise2_PfamHmmer1DB * obj
	INIT:
Wise2_RandomModel * temp;
	CODE:
	temp = Wise2_hard_link_RandomModel(Wise2_access_def_PfamHmmer1DB(obj));
	RETVAL = temp;
	OUTPUT:
	RETVAL




Wise2_PfamHmmer1DB *
new(class)
	char * class
	PPCODE:
	Wise2_PfamHmmer1DB * out;
	out = Wise2_PfamHmmer1DB_alloc_std();
	ST(0) = sv_newmortal();
	sv_setref_pv(ST(0),class,(void*)out);
	XSRETURN(1);

void
DESTROY(obj)
	Wise2_PfamHmmer1DB * obj
	CODE:
	Wise2_free_PfamHmmer1DB(obj);

void
each_en(obj)
	Wise2_PfamHmmer1DB * obj
	PPCODE:
	int i=0;
	int len;
	SV* temp;
	len = Wise2_length_en_PfamHmmer1DB(obj);
	for(i=0;i<len;i++){
	  temp = sv_newmortal();
	  sv_setref_pv(temp, "Wise2::PfamHmmer1Entry", (void*) (Wise2_hard_link_PfamHmmer1Entry(Wise2_access_en_PfamHmmer1DB(obj,i))));
	  XPUSHs(temp);
	  }
	XSRETURN(len);



MODULE = Wise2 PACKAGE = Wise2::PfamHmmer1Entry

Wise2_PfamHmmer1Entry *
hard_link_PfamHmmer1Entry(obj)
	Wise2_PfamHmmer1Entry * obj
	CODE:
	RETVAL = Wise2_hard_link_PfamHmmer1Entry(obj);
	OUTPUT:
	RETVAL



Wise2_PfamHmmer1Entry *
alloc()
	CODE:
	RETVAL = Wise2_PfamHmmer1Entry_alloc();
	OUTPUT:
	RETVAL



boolean
set_entryname(obj,entryname)
	Wise2_PfamHmmer1Entry * obj
	char * entryname
	CODE:
	RETVAL = Wise2_replace_entryname_PfamHmmer1Entry(obj,Wise2_stringalloc(entryname));
	OUTPUT:
	RETVAL



char *
entryname(obj)
	Wise2_PfamHmmer1Entry * obj
	INIT:
	char * temp;
	CODE:
	temp = Wise2_stringalloc(Wise2_access_entryname_PfamHmmer1Entry(obj));
	RETVAL = temp;
	OUTPUT:
	RETVAL



boolean
set_is_random(obj,is_random)
	Wise2_PfamHmmer1Entry * obj
	boolean is_random
	CODE:
	RETVAL = Wise2_replace_is_random_PfamHmmer1Entry(obj,is_random);
	OUTPUT:
	RETVAL



boolean
is_random(obj)
	Wise2_PfamHmmer1Entry * obj
	CODE:
	RETVAL = Wise2_access_is_random_PfamHmmer1Entry(obj);
	OUTPUT:
	RETVAL



boolean
set_is_hmmls(obj,is_hmmls)
	Wise2_PfamHmmer1Entry * obj
	boolean is_hmmls
	CODE:
	RETVAL = Wise2_replace_is_hmmls_PfamHmmer1Entry(obj,is_hmmls);
	OUTPUT:
	RETVAL



boolean
is_hmmls(obj)
	Wise2_PfamHmmer1Entry * obj
	CODE:
	RETVAL = Wise2_access_is_hmmls_PfamHmmer1Entry(obj);
	OUTPUT:
	RETVAL



boolean
set_bits_cutoff(obj,bits_cutoff)
	Wise2_PfamHmmer1Entry * obj
	double bits_cutoff
	CODE:
	RETVAL = Wise2_replace_bits_cutoff_PfamHmmer1Entry(obj,bits_cutoff);
	OUTPUT:
	RETVAL



double
bits_cutoff(obj)
	Wise2_PfamHmmer1Entry * obj
	CODE:
	RETVAL = Wise2_access_bits_cutoff_PfamHmmer1Entry(obj);
	OUTPUT:
	RETVAL




Wise2_PfamHmmer1Entry *
new(class)
	char * class
	PPCODE:
	Wise2_PfamHmmer1Entry * out;
	out = Wise2_PfamHmmer1Entry_alloc();
	ST(0) = sv_newmortal();
	sv_setref_pv(ST(0),class,(void*)out);
	XSRETURN(1);

void
DESTROY(obj)
	Wise2_PfamHmmer1Entry * obj
	CODE:
	Wise2_free_PfamHmmer1Entry(obj);



MODULE = Wise2 PACKAGE = Wise2

