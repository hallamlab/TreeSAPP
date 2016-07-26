

MODULE = Wise2 PACKAGE = Wise2::MatchSummarySet

Wise2_MatchSummarySet *
MatchSummarySet_from_AlnBlock_estwise(alb,qname,offset,target)
	Wise2_AlnBlock * alb
	char * qname
	int offset
	Wise2_Sequence * target
	CODE:
	RETVAL = Wise2_MatchSummarySet_from_AlnBlock_estwise(alb,qname,offset,target);
	OUTPUT:
	RETVAL



Wise2_MatchSummarySet *
MatchSummarySet_from_AlnBlock_genewise(alb,qname,protoff,target)
	Wise2_AlnBlock * alb
	char * qname
	int protoff
	Wise2_Sequence * target
	CODE:
	RETVAL = Wise2_MatchSummarySet_from_AlnBlock_genewise(alb,qname,protoff,target);
	OUTPUT:
	RETVAL



Wise2_MatchSummarySet *
hard_link_MatchSummarySet(obj)
	Wise2_MatchSummarySet * obj
	CODE:
	RETVAL = Wise2_hard_link_MatchSummarySet(obj);
	OUTPUT:
	RETVAL



Wise2_MatchSummarySet *
MatchSummarySet_alloc_std()
	CODE:
	RETVAL = Wise2_MatchSummarySet_alloc_std();
	OUTPUT:
	RETVAL



Wise2_MatchSummary *
ms(obj,i)
	Wise2_MatchSummarySet * obj
	int i
	INIT:
Wise2_MatchSummary * temp;
	CODE:
	temp = Wise2_hard_link_MatchSummary(Wise2_access_ms_MatchSummarySet(obj,i));
	RETVAL = temp;
	OUTPUT:
	RETVAL



int
length_ms(obj)
	Wise2_MatchSummarySet * obj
	CODE:
	RETVAL = Wise2_length_ms_MatchSummarySet(obj);
	OUTPUT:
	RETVAL



int
flush_ms(obj)
	Wise2_MatchSummarySet * obj
	CODE:
	RETVAL = Wise2_flush_MatchSummarySet(obj);
	OUTPUT:
	RETVAL



boolean
add_ms(obj,add)
	Wise2_MatchSummarySet * obj
	Wise2_MatchSummary * add
	CODE:
	RETVAL = Wise2_add_MatchSummarySet(obj,Wise2_hard_link_MatchSummary(add));
	OUTPUT:
	RETVAL




Wise2_MatchSummarySet *
new(class)
	char * class
	PPCODE:
	Wise2_MatchSummarySet * out;
	out = Wise2_MatchSummarySet_alloc_std();
	ST(0) = sv_newmortal();
	sv_setref_pv(ST(0),class,(void*)out);
	XSRETURN(1);

void
DESTROY(obj)
	Wise2_MatchSummarySet * obj
	CODE:
	Wise2_free_MatchSummarySet(obj);

void
each_ms(obj)
	Wise2_MatchSummarySet * obj
	PPCODE:
	int i=0;
	int len;
	SV* temp;
	len = Wise2_length_ms_MatchSummarySet(obj);
	for(i=0;i<len;i++){
	  temp = sv_newmortal();
	  sv_setref_pv(temp, "Wise2::MatchSummary", (void*) (Wise2_hard_link_MatchSummary(Wise2_access_ms_MatchSummarySet(obj,i))));
	  XPUSHs(temp);
	  }
	XSRETURN(len);



MODULE = Wise2 PACKAGE = Wise2::MatchSummary

Wise2_MatchSummary *
hard_link_MatchSummary(obj)
	Wise2_MatchSummary * obj
	CODE:
	RETVAL = Wise2_hard_link_MatchSummary(obj);
	OUTPUT:
	RETVAL



Wise2_MatchSummary *
alloc()
	CODE:
	RETVAL = Wise2_MatchSummary_alloc();
	OUTPUT:
	RETVAL



boolean
set_bits(obj,bits)
	Wise2_MatchSummary * obj
	double bits
	CODE:
	RETVAL = Wise2_replace_bits_MatchSummary(obj,bits);
	OUTPUT:
	RETVAL



double
bits(obj)
	Wise2_MatchSummary * obj
	CODE:
	RETVAL = Wise2_access_bits_MatchSummary(obj);
	OUTPUT:
	RETVAL



boolean
set_qname(obj,qname)
	Wise2_MatchSummary * obj
	char * qname
	CODE:
	RETVAL = Wise2_replace_qname_MatchSummary(obj,Wise2_stringalloc(qname));
	OUTPUT:
	RETVAL



char *
qname(obj)
	Wise2_MatchSummary * obj
	INIT:
	char * temp;
	CODE:
	temp = Wise2_stringalloc(Wise2_access_qname_MatchSummary(obj));
	RETVAL = temp;
	OUTPUT:
	RETVAL



boolean
set_tname(obj,tname)
	Wise2_MatchSummary * obj
	char * tname
	CODE:
	RETVAL = Wise2_replace_tname_MatchSummary(obj,Wise2_stringalloc(tname));
	OUTPUT:
	RETVAL



char *
tname(obj)
	Wise2_MatchSummary * obj
	INIT:
	char * temp;
	CODE:
	temp = Wise2_stringalloc(Wise2_access_tname_MatchSummary(obj));
	RETVAL = temp;
	OUTPUT:
	RETVAL



boolean
set_qstart(obj,qstart)
	Wise2_MatchSummary * obj
	int qstart
	CODE:
	RETVAL = Wise2_replace_qstart_MatchSummary(obj,qstart);
	OUTPUT:
	RETVAL



int
qstart(obj)
	Wise2_MatchSummary * obj
	CODE:
	RETVAL = Wise2_access_qstart_MatchSummary(obj);
	OUTPUT:
	RETVAL



boolean
set_qend(obj,qend)
	Wise2_MatchSummary * obj
	int qend
	CODE:
	RETVAL = Wise2_replace_qend_MatchSummary(obj,qend);
	OUTPUT:
	RETVAL



int
qend(obj)
	Wise2_MatchSummary * obj
	CODE:
	RETVAL = Wise2_access_qend_MatchSummary(obj);
	OUTPUT:
	RETVAL



boolean
set_tstart(obj,tstart)
	Wise2_MatchSummary * obj
	int tstart
	CODE:
	RETVAL = Wise2_replace_tstart_MatchSummary(obj,tstart);
	OUTPUT:
	RETVAL



int
tstart(obj)
	Wise2_MatchSummary * obj
	CODE:
	RETVAL = Wise2_access_tstart_MatchSummary(obj);
	OUTPUT:
	RETVAL



boolean
set_tend(obj,tend)
	Wise2_MatchSummary * obj
	int tend
	CODE:
	RETVAL = Wise2_replace_tend_MatchSummary(obj,tend);
	OUTPUT:
	RETVAL



int
tend(obj)
	Wise2_MatchSummary * obj
	CODE:
	RETVAL = Wise2_access_tend_MatchSummary(obj);
	OUTPUT:
	RETVAL



boolean
set_qintron(obj,qintron)
	Wise2_MatchSummary * obj
	int qintron
	CODE:
	RETVAL = Wise2_replace_qintron_MatchSummary(obj,qintron);
	OUTPUT:
	RETVAL



int
qintron(obj)
	Wise2_MatchSummary * obj
	CODE:
	RETVAL = Wise2_access_qintron_MatchSummary(obj);
	OUTPUT:
	RETVAL



boolean
set_qframeshift(obj,qframeshift)
	Wise2_MatchSummary * obj
	int qframeshift
	CODE:
	RETVAL = Wise2_replace_qframeshift_MatchSummary(obj,qframeshift);
	OUTPUT:
	RETVAL



int
qframeshift(obj)
	Wise2_MatchSummary * obj
	CODE:
	RETVAL = Wise2_access_qframeshift_MatchSummary(obj);
	OUTPUT:
	RETVAL



boolean
set_tintron(obj,tintron)
	Wise2_MatchSummary * obj
	int tintron
	CODE:
	RETVAL = Wise2_replace_tintron_MatchSummary(obj,tintron);
	OUTPUT:
	RETVAL



int
tintron(obj)
	Wise2_MatchSummary * obj
	CODE:
	RETVAL = Wise2_access_tintron_MatchSummary(obj);
	OUTPUT:
	RETVAL



boolean
set_tframeshift(obj,tframeshift)
	Wise2_MatchSummary * obj
	int tframeshift
	CODE:
	RETVAL = Wise2_replace_tframeshift_MatchSummary(obj,tframeshift);
	OUTPUT:
	RETVAL



int
tframeshift(obj)
	Wise2_MatchSummary * obj
	CODE:
	RETVAL = Wise2_access_tframeshift_MatchSummary(obj);
	OUTPUT:
	RETVAL




Wise2_MatchSummary *
new(class)
	char * class
	PPCODE:
	Wise2_MatchSummary * out;
	out = Wise2_MatchSummary_alloc();
	ST(0) = sv_newmortal();
	sv_setref_pv(ST(0),class,(void*)out);
	XSRETURN(1);

void
DESTROY(obj)
	Wise2_MatchSummary * obj
	CODE:
	Wise2_free_MatchSummary(obj);



MODULE = Wise2 PACKAGE = Wise2

