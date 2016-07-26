

MODULE = Wise2 PACKAGE = Wise2::AlnRange

Wise2_AlnRange *
hard_link_AlnRange(obj)
	Wise2_AlnRange * obj
	CODE:
	RETVAL = Wise2_hard_link_AlnRange(obj);
	OUTPUT:
	RETVAL



Wise2_AlnRange *
alloc()
	CODE:
	RETVAL = Wise2_AlnRange_alloc();
	OUTPUT:
	RETVAL



boolean
set_starti(obj,starti)
	Wise2_AlnRange * obj
	int starti
	CODE:
	RETVAL = Wise2_replace_starti_AlnRange(obj,starti);
	OUTPUT:
	RETVAL



int
starti(obj)
	Wise2_AlnRange * obj
	CODE:
	RETVAL = Wise2_access_starti_AlnRange(obj);
	OUTPUT:
	RETVAL



boolean
set_startj(obj,startj)
	Wise2_AlnRange * obj
	int startj
	CODE:
	RETVAL = Wise2_replace_startj_AlnRange(obj,startj);
	OUTPUT:
	RETVAL



int
startj(obj)
	Wise2_AlnRange * obj
	CODE:
	RETVAL = Wise2_access_startj_AlnRange(obj);
	OUTPUT:
	RETVAL



boolean
set_startstate(obj,startstate)
	Wise2_AlnRange * obj
	int startstate
	CODE:
	RETVAL = Wise2_replace_startstate_AlnRange(obj,startstate);
	OUTPUT:
	RETVAL



int
startstate(obj)
	Wise2_AlnRange * obj
	CODE:
	RETVAL = Wise2_access_startstate_AlnRange(obj);
	OUTPUT:
	RETVAL



boolean
set_stopi(obj,stopi)
	Wise2_AlnRange * obj
	int stopi
	CODE:
	RETVAL = Wise2_replace_stopi_AlnRange(obj,stopi);
	OUTPUT:
	RETVAL



int
stopi(obj)
	Wise2_AlnRange * obj
	CODE:
	RETVAL = Wise2_access_stopi_AlnRange(obj);
	OUTPUT:
	RETVAL



boolean
set_stopj(obj,stopj)
	Wise2_AlnRange * obj
	int stopj
	CODE:
	RETVAL = Wise2_replace_stopj_AlnRange(obj,stopj);
	OUTPUT:
	RETVAL



int
stopj(obj)
	Wise2_AlnRange * obj
	CODE:
	RETVAL = Wise2_access_stopj_AlnRange(obj);
	OUTPUT:
	RETVAL



boolean
set_stopstate(obj,stopstate)
	Wise2_AlnRange * obj
	int stopstate
	CODE:
	RETVAL = Wise2_replace_stopstate_AlnRange(obj,stopstate);
	OUTPUT:
	RETVAL



int
stopstate(obj)
	Wise2_AlnRange * obj
	CODE:
	RETVAL = Wise2_access_stopstate_AlnRange(obj);
	OUTPUT:
	RETVAL



boolean
set_startscore(obj,startscore)
	Wise2_AlnRange * obj
	int startscore
	CODE:
	RETVAL = Wise2_replace_startscore_AlnRange(obj,startscore);
	OUTPUT:
	RETVAL



int
startscore(obj)
	Wise2_AlnRange * obj
	CODE:
	RETVAL = Wise2_access_startscore_AlnRange(obj);
	OUTPUT:
	RETVAL



boolean
set_stopscore(obj,stopscore)
	Wise2_AlnRange * obj
	int stopscore
	CODE:
	RETVAL = Wise2_replace_stopscore_AlnRange(obj,stopscore);
	OUTPUT:
	RETVAL



int
stopscore(obj)
	Wise2_AlnRange * obj
	CODE:
	RETVAL = Wise2_access_stopscore_AlnRange(obj);
	OUTPUT:
	RETVAL




Wise2_AlnRange *
new(class)
	char * class
	PPCODE:
	Wise2_AlnRange * out;
	out = Wise2_AlnRange_alloc();
	ST(0) = sv_newmortal();
	sv_setref_pv(ST(0),class,(void*)out);
	XSRETURN(1);

void
DESTROY(obj)
	Wise2_AlnRange * obj
	CODE:
	Wise2_free_AlnRange(obj);



MODULE = Wise2 PACKAGE = Wise2::AlnRangeSet

Wise2_AlnRangeSet *
hard_link_AlnRangeSet(obj)
	Wise2_AlnRangeSet * obj
	CODE:
	RETVAL = Wise2_hard_link_AlnRangeSet(obj);
	OUTPUT:
	RETVAL



Wise2_AlnRangeSet *
AlnRangeSet_alloc_std()
	CODE:
	RETVAL = Wise2_AlnRangeSet_alloc_std();
	OUTPUT:
	RETVAL



boolean
set_score(obj,score)
	Wise2_AlnRangeSet * obj
	int score
	CODE:
	RETVAL = Wise2_replace_score_AlnRangeSet(obj,score);
	OUTPUT:
	RETVAL



int
score(obj)
	Wise2_AlnRangeSet * obj
	CODE:
	RETVAL = Wise2_access_score_AlnRangeSet(obj);
	OUTPUT:
	RETVAL



Wise2_AlnRange *
alr(obj,i)
	Wise2_AlnRangeSet * obj
	int i
	INIT:
Wise2_AlnRange * temp;
	CODE:
	temp = Wise2_hard_link_AlnRange(Wise2_access_alr_AlnRangeSet(obj,i));
	RETVAL = temp;
	OUTPUT:
	RETVAL



int
length_alr(obj)
	Wise2_AlnRangeSet * obj
	CODE:
	RETVAL = Wise2_length_alr_AlnRangeSet(obj);
	OUTPUT:
	RETVAL



int
flush_alr(obj)
	Wise2_AlnRangeSet * obj
	CODE:
	RETVAL = Wise2_flush_AlnRangeSet(obj);
	OUTPUT:
	RETVAL



boolean
add_alr(obj,add)
	Wise2_AlnRangeSet * obj
	Wise2_AlnRange * add
	CODE:
	RETVAL = Wise2_add_AlnRangeSet(obj,Wise2_hard_link_AlnRange(add));
	OUTPUT:
	RETVAL




Wise2_AlnRangeSet *
new(class)
	char * class
	PPCODE:
	Wise2_AlnRangeSet * out;
	out = Wise2_AlnRangeSet_alloc_std();
	ST(0) = sv_newmortal();
	sv_setref_pv(ST(0),class,(void*)out);
	XSRETURN(1);

void
DESTROY(obj)
	Wise2_AlnRangeSet * obj
	CODE:
	Wise2_free_AlnRangeSet(obj);

void
each_alr(obj)
	Wise2_AlnRangeSet * obj
	PPCODE:
	int i=0;
	int len;
	SV* temp;
	len = Wise2_length_alr_AlnRangeSet(obj);
	for(i=0;i<len;i++){
	  temp = sv_newmortal();
	  sv_setref_pv(temp, "Wise2::AlnRange", (void*) (Wise2_hard_link_AlnRange(Wise2_access_alr_AlnRangeSet(obj,i))));
	  XPUSHs(temp);
	  }
	XSRETURN(len);



MODULE = Wise2 PACKAGE = Wise2

