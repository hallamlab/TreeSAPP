

MODULE = Wise2 PACKAGE = Wise2::PackAln

void
show_simple_PackAln(pal,ofp)
	Wise2_PackAln * pal
	FILE * ofp
	CODE:
	Wise2_show_simple_PackAln(pal,ofp);



void
show_bits_and_cumlative_PackAln(pal,ofp)
	Wise2_PackAln * pal
	FILE * ofp
	CODE:
	Wise2_show_bits_and_cumlative_PackAln(pal,ofp);



Wise2_PackAln *
hard_link_PackAln(obj)
	Wise2_PackAln * obj
	CODE:
	RETVAL = Wise2_hard_link_PackAln(obj);
	OUTPUT:
	RETVAL



Wise2_PackAln *
PackAln_alloc_std()
	CODE:
	RETVAL = Wise2_PackAln_alloc_std();
	OUTPUT:
	RETVAL



Wise2_PackAlnUnit *
pau(obj,i)
	Wise2_PackAln * obj
	int i
	INIT:
Wise2_PackAlnUnit * temp;
	CODE:
	temp = Wise2_hard_link_PackAlnUnit(Wise2_access_pau_PackAln(obj,i));
	RETVAL = temp;
	OUTPUT:
	RETVAL



int
length_pau(obj)
	Wise2_PackAln * obj
	CODE:
	RETVAL = Wise2_length_pau_PackAln(obj);
	OUTPUT:
	RETVAL



int
flush_pau(obj)
	Wise2_PackAln * obj
	CODE:
	RETVAL = Wise2_flush_PackAln(obj);
	OUTPUT:
	RETVAL



boolean
add_pau(obj,add)
	Wise2_PackAln * obj
	Wise2_PackAlnUnit * add
	CODE:
	RETVAL = Wise2_add_PackAln(obj,Wise2_hard_link_PackAlnUnit(add));
	OUTPUT:
	RETVAL



boolean
set_score(obj,score)
	Wise2_PackAln * obj
	int score
	CODE:
	RETVAL = Wise2_replace_score_PackAln(obj,score);
	OUTPUT:
	RETVAL



int
score(obj)
	Wise2_PackAln * obj
	CODE:
	RETVAL = Wise2_access_score_PackAln(obj);
	OUTPUT:
	RETVAL




Wise2_PackAln *
new(class)
	char * class
	PPCODE:
	Wise2_PackAln * out;
	out = Wise2_PackAln_alloc_std();
	ST(0) = sv_newmortal();
	sv_setref_pv(ST(0),class,(void*)out);
	XSRETURN(1);

void
DESTROY(obj)
	Wise2_PackAln * obj
	CODE:
	Wise2_free_PackAln(obj);

void
each_pau(obj)
	Wise2_PackAln * obj
	PPCODE:
	int i=0;
	int len;
	SV* temp;
	len = Wise2_length_pau_PackAln(obj);
	for(i=0;i<len;i++){
	  temp = sv_newmortal();
	  sv_setref_pv(temp, "Wise2::PackAlnUnit", (void*) (Wise2_hard_link_PackAlnUnit(Wise2_access_pau_PackAln(obj,i))));
	  XPUSHs(temp);
	  }
	XSRETURN(len);



MODULE = Wise2 PACKAGE = Wise2::PackAlnUnit

Wise2_PackAlnUnit *
hard_link_PackAlnUnit(obj)
	Wise2_PackAlnUnit * obj
	CODE:
	RETVAL = Wise2_hard_link_PackAlnUnit(obj);
	OUTPUT:
	RETVAL



Wise2_PackAlnUnit *
alloc()
	CODE:
	RETVAL = Wise2_PackAlnUnit_alloc();
	OUTPUT:
	RETVAL



boolean
set_i(obj,i)
	Wise2_PackAlnUnit * obj
	int i
	CODE:
	RETVAL = Wise2_replace_i_PackAlnUnit(obj,i);
	OUTPUT:
	RETVAL



int
i(obj)
	Wise2_PackAlnUnit * obj
	CODE:
	RETVAL = Wise2_access_i_PackAlnUnit(obj);
	OUTPUT:
	RETVAL



boolean
set_j(obj,j)
	Wise2_PackAlnUnit * obj
	int j
	CODE:
	RETVAL = Wise2_replace_j_PackAlnUnit(obj,j);
	OUTPUT:
	RETVAL



int
j(obj)
	Wise2_PackAlnUnit * obj
	CODE:
	RETVAL = Wise2_access_j_PackAlnUnit(obj);
	OUTPUT:
	RETVAL



boolean
set_state(obj,state)
	Wise2_PackAlnUnit * obj
	int state
	CODE:
	RETVAL = Wise2_replace_state_PackAlnUnit(obj,state);
	OUTPUT:
	RETVAL



int
state(obj)
	Wise2_PackAlnUnit * obj
	CODE:
	RETVAL = Wise2_access_state_PackAlnUnit(obj);
	OUTPUT:
	RETVAL



boolean
set_score(obj,score)
	Wise2_PackAlnUnit * obj
	int score
	CODE:
	RETVAL = Wise2_replace_score_PackAlnUnit(obj,score);
	OUTPUT:
	RETVAL



int
score(obj)
	Wise2_PackAlnUnit * obj
	CODE:
	RETVAL = Wise2_access_score_PackAlnUnit(obj);
	OUTPUT:
	RETVAL




Wise2_PackAlnUnit *
new(class)
	char * class
	PPCODE:
	Wise2_PackAlnUnit * out;
	out = Wise2_PackAlnUnit_alloc();
	ST(0) = sv_newmortal();
	sv_setref_pv(ST(0),class,(void*)out);
	XSRETURN(1);

void
DESTROY(obj)
	Wise2_PackAlnUnit * obj
	CODE:
	Wise2_free_PackAlnUnit(obj);



MODULE = Wise2 PACKAGE = Wise2

