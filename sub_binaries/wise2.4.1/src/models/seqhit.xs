

MODULE = Wise2 PACKAGE = Wise2::DnaSequenceHitList

void
show_DnaSequenceHitList(dsl,ofp)
	Wise2_DnaSequenceHitList * dsl
	FILE * ofp
	CODE:
	Wise2_show_DnaSequenceHitList(dsl,ofp);



Wise2_DnaSequenceHitList *
read_MSPcrunch_DnaSequenceHitList(ifp)
	FILE * ifp
	CODE:
	RETVAL = Wise2_read_MSPcrunch_DnaSequenceHitList(ifp);
	OUTPUT:
	RETVAL



Wise2_DnaSequenceHitList *
hard_link_DnaSequenceHitList(obj)
	Wise2_DnaSequenceHitList * obj
	CODE:
	RETVAL = Wise2_hard_link_DnaSequenceHitList(obj);
	OUTPUT:
	RETVAL



Wise2_DnaSequenceHitList *
alloc()
	CODE:
	RETVAL = Wise2_DnaSequenceHitList_alloc();
	OUTPUT:
	RETVAL



boolean
set_forward(obj,forward)
	Wise2_DnaSequenceHitList * obj
	Wise2_SegmentHitList * forward
	CODE:
	RETVAL = Wise2_replace_forward_DnaSequenceHitList(obj,Wise2_hard_link_SegmentHitList(forward));
	OUTPUT:
	RETVAL



Wise2_SegmentHitList *
forward(obj)
	Wise2_DnaSequenceHitList * obj
	INIT:
Wise2_SegmentHitList * temp;
	CODE:
	temp = Wise2_hard_link_SegmentHitList(Wise2_access_forward_DnaSequenceHitList(obj));
	RETVAL = temp;
	OUTPUT:
	RETVAL



boolean
set_backward(obj,backward)
	Wise2_DnaSequenceHitList * obj
	Wise2_SegmentHitList * backward
	CODE:
	RETVAL = Wise2_replace_backward_DnaSequenceHitList(obj,Wise2_hard_link_SegmentHitList(backward));
	OUTPUT:
	RETVAL



Wise2_SegmentHitList *
backward(obj)
	Wise2_DnaSequenceHitList * obj
	INIT:
Wise2_SegmentHitList * temp;
	CODE:
	temp = Wise2_hard_link_SegmentHitList(Wise2_access_backward_DnaSequenceHitList(obj));
	RETVAL = temp;
	OUTPUT:
	RETVAL




Wise2_DnaSequenceHitList *
new(class)
	char * class
	PPCODE:
	Wise2_DnaSequenceHitList * out;
	out = Wise2_DnaSequenceHitList_alloc();
	ST(0) = sv_newmortal();
	sv_setref_pv(ST(0),class,(void*)out);
	XSRETURN(1);

void
DESTROY(obj)
	Wise2_DnaSequenceHitList * obj
	CODE:
	Wise2_free_DnaSequenceHitList(obj);



MODULE = Wise2 PACKAGE = Wise2::SegmentHitList

Wise2_SegmentHitList *
hard_link_SegmentHitList(obj)
	Wise2_SegmentHitList * obj
	CODE:
	RETVAL = Wise2_hard_link_SegmentHitList(obj);
	OUTPUT:
	RETVAL



Wise2_SegmentHitList *
SegmentHitList_alloc_std()
	CODE:
	RETVAL = Wise2_SegmentHitList_alloc_std();
	OUTPUT:
	RETVAL



Wise2_SegmentHit *
seghit(obj,i)
	Wise2_SegmentHitList * obj
	int i
	INIT:
Wise2_SegmentHit * temp;
	CODE:
	temp = Wise2_hard_link_SegmentHit(Wise2_access_seghit_SegmentHitList(obj,i));
	RETVAL = temp;
	OUTPUT:
	RETVAL



int
length_seghit(obj)
	Wise2_SegmentHitList * obj
	CODE:
	RETVAL = Wise2_length_seghit_SegmentHitList(obj);
	OUTPUT:
	RETVAL



int
flush_seghit(obj)
	Wise2_SegmentHitList * obj
	CODE:
	RETVAL = Wise2_flush_SegmentHitList(obj);
	OUTPUT:
	RETVAL



boolean
add_seghit(obj,add)
	Wise2_SegmentHitList * obj
	Wise2_SegmentHit * add
	CODE:
	RETVAL = Wise2_add_SegmentHitList(obj,Wise2_hard_link_SegmentHit(add));
	OUTPUT:
	RETVAL




Wise2_SegmentHitList *
new(class)
	char * class
	PPCODE:
	Wise2_SegmentHitList * out;
	out = Wise2_SegmentHitList_alloc_std();
	ST(0) = sv_newmortal();
	sv_setref_pv(ST(0),class,(void*)out);
	XSRETURN(1);

void
DESTROY(obj)
	Wise2_SegmentHitList * obj
	CODE:
	Wise2_free_SegmentHitList(obj);

void
each_seghit(obj)
	Wise2_SegmentHitList * obj
	PPCODE:
	int i=0;
	int len;
	SV* temp;
	len = Wise2_length_seghit_SegmentHitList(obj);
	for(i=0;i<len;i++){
	  temp = sv_newmortal();
	  sv_setref_pv(temp, "Wise2::SegmentHit", (void*) (Wise2_hard_link_SegmentHit(Wise2_access_seghit_SegmentHitList(obj,i))));
	  XPUSHs(temp);
	  }
	XSRETURN(len);



MODULE = Wise2 PACKAGE = Wise2::SegmentHit

Wise2_SegmentHit *
hard_link_SegmentHit(obj)
	Wise2_SegmentHit * obj
	CODE:
	RETVAL = Wise2_hard_link_SegmentHit(obj);
	OUTPUT:
	RETVAL



Wise2_SegmentHit *
alloc()
	CODE:
	RETVAL = Wise2_SegmentHit_alloc();
	OUTPUT:
	RETVAL



boolean
set_name(obj,name)
	Wise2_SegmentHit * obj
	char * name
	CODE:
	RETVAL = Wise2_replace_name_SegmentHit(obj,Wise2_stringalloc(name));
	OUTPUT:
	RETVAL



char *
name(obj)
	Wise2_SegmentHit * obj
	INIT:
	char * temp;
	CODE:
	temp = Wise2_stringalloc(Wise2_access_name_SegmentHit(obj));
	RETVAL = temp;
	OUTPUT:
	RETVAL



boolean
set_qstart(obj,qstart)
	Wise2_SegmentHit * obj
	int qstart
	CODE:
	RETVAL = Wise2_replace_qstart_SegmentHit(obj,qstart);
	OUTPUT:
	RETVAL



int
qstart(obj)
	Wise2_SegmentHit * obj
	CODE:
	RETVAL = Wise2_access_qstart_SegmentHit(obj);
	OUTPUT:
	RETVAL



boolean
set_qend(obj,qend)
	Wise2_SegmentHit * obj
	int qend
	CODE:
	RETVAL = Wise2_replace_qend_SegmentHit(obj,qend);
	OUTPUT:
	RETVAL



int
qend(obj)
	Wise2_SegmentHit * obj
	CODE:
	RETVAL = Wise2_access_qend_SegmentHit(obj);
	OUTPUT:
	RETVAL



boolean
set_tstart(obj,tstart)
	Wise2_SegmentHit * obj
	int tstart
	CODE:
	RETVAL = Wise2_replace_tstart_SegmentHit(obj,tstart);
	OUTPUT:
	RETVAL



int
tstart(obj)
	Wise2_SegmentHit * obj
	CODE:
	RETVAL = Wise2_access_tstart_SegmentHit(obj);
	OUTPUT:
	RETVAL



boolean
set_tend(obj,tend)
	Wise2_SegmentHit * obj
	int tend
	CODE:
	RETVAL = Wise2_replace_tend_SegmentHit(obj,tend);
	OUTPUT:
	RETVAL



int
tend(obj)
	Wise2_SegmentHit * obj
	CODE:
	RETVAL = Wise2_access_tend_SegmentHit(obj);
	OUTPUT:
	RETVAL



boolean
set_score(obj,score)
	Wise2_SegmentHit * obj
	double score
	CODE:
	RETVAL = Wise2_replace_score_SegmentHit(obj,score);
	OUTPUT:
	RETVAL



double
score(obj)
	Wise2_SegmentHit * obj
	CODE:
	RETVAL = Wise2_access_score_SegmentHit(obj);
	OUTPUT:
	RETVAL



boolean
set_next_hit(obj,next_hit)
	Wise2_SegmentHit * obj
	Wise2_SegmentHit * next_hit
	CODE:
	RETVAL = Wise2_replace_next_hit_SegmentHit(obj,Wise2_hard_link_SegmentHit(next_hit));
	OUTPUT:
	RETVAL



Wise2_SegmentHit *
next_hit(obj)
	Wise2_SegmentHit * obj
	INIT:
Wise2_SegmentHit * temp;
	CODE:
	temp = Wise2_hard_link_SegmentHit(Wise2_access_next_hit_SegmentHit(obj));
	RETVAL = temp;
	OUTPUT:
	RETVAL




Wise2_SegmentHit *
new(class)
	char * class
	PPCODE:
	Wise2_SegmentHit * out;
	out = Wise2_SegmentHit_alloc();
	ST(0) = sv_newmortal();
	sv_setref_pv(ST(0),class,(void*)out);
	XSRETURN(1);

void
DESTROY(obj)
	Wise2_SegmentHit * obj
	CODE:
	Wise2_free_SegmentHit(obj);



MODULE = Wise2 PACKAGE = Wise2

