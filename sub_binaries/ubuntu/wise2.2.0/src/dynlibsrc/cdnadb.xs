

MODULE = Wise2 PACKAGE = Wise2::cDNADB

Wise2_cDNA *
get_entry(cdnadb,de)
	Wise2_cDNADB * cdnadb
	Wise2_DataEntry * de
	CODE:
	RETVAL = Wise2_get_cDNA_from_cDNADB(cdnadb,de);
	OUTPUT:
	RETVAL



Wise2_cDNADB *
hard_link_cDNADB(obj)
	Wise2_cDNADB * obj
	CODE:
	RETVAL = Wise2_hard_link_cDNADB(obj);
	OUTPUT:
	RETVAL



Wise2_cDNADB *
alloc()
	CODE:
	RETVAL = Wise2_cDNADB_alloc();
	OUTPUT:
	RETVAL



boolean
set_is_single_seq(obj,is_single_seq)
	Wise2_cDNADB * obj
	boolean is_single_seq
	CODE:
	RETVAL = Wise2_replace_is_single_seq_cDNADB(obj,is_single_seq);
	OUTPUT:
	RETVAL



boolean
is_single_seq(obj)
	Wise2_cDNADB * obj
	CODE:
	RETVAL = Wise2_access_is_single_seq_cDNADB(obj);
	OUTPUT:
	RETVAL



boolean
set_done_forward(obj,done_forward)
	Wise2_cDNADB * obj
	boolean done_forward
	CODE:
	RETVAL = Wise2_replace_done_forward_cDNADB(obj,done_forward);
	OUTPUT:
	RETVAL



boolean
done_forward(obj)
	Wise2_cDNADB * obj
	CODE:
	RETVAL = Wise2_access_done_forward_cDNADB(obj);
	OUTPUT:
	RETVAL



boolean
set_forward_only(obj,forward_only)
	Wise2_cDNADB * obj
	boolean forward_only
	CODE:
	RETVAL = Wise2_replace_forward_only_cDNADB(obj,forward_only);
	OUTPUT:
	RETVAL



boolean
forward_only(obj)
	Wise2_cDNADB * obj
	CODE:
	RETVAL = Wise2_access_forward_only_cDNADB(obj);
	OUTPUT:
	RETVAL



boolean
set_forw(obj,forw)
	Wise2_cDNADB * obj
	Wise2_ComplexSequence * forw
	CODE:
	RETVAL = Wise2_replace_forw_cDNADB(obj,Wise2_hard_link_ComplexSequence(forw));
	OUTPUT:
	RETVAL



Wise2_ComplexSequence *
forw(obj)
	Wise2_cDNADB * obj
	INIT:
Wise2_ComplexSequence * temp;
	CODE:
	temp = Wise2_hard_link_ComplexSequence(Wise2_access_forw_cDNADB(obj));
	RETVAL = temp;
	OUTPUT:
	RETVAL



boolean
set_rev(obj,rev)
	Wise2_cDNADB * obj
	Wise2_ComplexSequence * rev
	CODE:
	RETVAL = Wise2_replace_rev_cDNADB(obj,Wise2_hard_link_ComplexSequence(rev));
	OUTPUT:
	RETVAL



Wise2_ComplexSequence *
rev(obj)
	Wise2_cDNADB * obj
	INIT:
Wise2_ComplexSequence * temp;
	CODE:
	temp = Wise2_hard_link_ComplexSequence(Wise2_access_rev_cDNADB(obj));
	RETVAL = temp;
	OUTPUT:
	RETVAL



boolean
set_sdb(obj,sdb)
	Wise2_cDNADB * obj
	Wise2_SequenceDB * sdb
	CODE:
	RETVAL = Wise2_replace_sdb_cDNADB(obj,Wise2_hard_link_SequenceDB(sdb));
	OUTPUT:
	RETVAL



Wise2_SequenceDB *
sdb(obj)
	Wise2_cDNADB * obj
	INIT:
Wise2_SequenceDB * temp;
	CODE:
	temp = Wise2_hard_link_SequenceDB(Wise2_access_sdb_cDNADB(obj));
	RETVAL = temp;
	OUTPUT:
	RETVAL



boolean
set_current(obj,current)
	Wise2_cDNADB * obj
	Wise2_Sequence * current
	CODE:
	RETVAL = Wise2_replace_current_cDNADB(obj,Wise2_hard_link_Sequence(current));
	OUTPUT:
	RETVAL



Wise2_Sequence *
current(obj)
	Wise2_cDNADB * obj
	INIT:
Wise2_Sequence * temp;
	CODE:
	temp = Wise2_hard_link_Sequence(Wise2_access_current_cDNADB(obj));
	RETVAL = temp;
	OUTPUT:
	RETVAL



boolean
set_cses(obj,cses)
	Wise2_cDNADB * obj
	Wise2_ComplexSequenceEvalSet * cses
	CODE:
	RETVAL = Wise2_replace_cses_cDNADB(obj,Wise2_hard_link_ComplexSequenceEvalSet(cses));
	OUTPUT:
	RETVAL



Wise2_ComplexSequenceEvalSet *
cses(obj)
	Wise2_cDNADB * obj
	INIT:
Wise2_ComplexSequenceEvalSet * temp;
	CODE:
	temp = Wise2_hard_link_ComplexSequenceEvalSet(Wise2_access_cses_cDNADB(obj));
	RETVAL = temp;
	OUTPUT:
	RETVAL




Wise2_cDNADB *
new(class)
	char * class
	PPCODE:
	Wise2_cDNADB * out;
	out = Wise2_cDNADB_alloc();
	ST(0) = sv_newmortal();
	sv_setref_pv(ST(0),class,(void*)out);
	XSRETURN(1);

void
DESTROY(obj)
	Wise2_cDNADB * obj
	CODE:
	Wise2_free_cDNADB(obj);



MODULE = Wise2 PACKAGE = Wise2

Wise2_cDNADB *
new_cDNADB_from_single_seq(seq)
	Wise2_cDNA * seq
	CODE:
	RETVAL = Wise2_new_cDNADB_from_single_seq(seq);
	OUTPUT:
	RETVAL



Wise2_cDNADB *
new_cDNADB(seqdb)
	Wise2_SequenceDB * seqdb
	CODE:
	RETVAL = Wise2_new_cDNADB(seqdb);
	OUTPUT:
	RETVAL



