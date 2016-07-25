

MODULE = Wise2 PACKAGE = Wise2::GenomicDB

Wise2_Genomic *
get_Genomic_from_GenomicDB(gendb,de)
	Wise2_GenomicDB * gendb
	Wise2_DataEntry * de
	CODE:
	RETVAL = Wise2_get_Genomic_from_GenomicDB(gendb,de);
	OUTPUT:
	RETVAL



Wise2_GenomicDB *
hard_link_GenomicDB(obj)
	Wise2_GenomicDB * obj
	CODE:
	RETVAL = Wise2_hard_link_GenomicDB(obj);
	OUTPUT:
	RETVAL



Wise2_GenomicDB *
alloc()
	CODE:
	RETVAL = Wise2_GenomicDB_alloc();
	OUTPUT:
	RETVAL



boolean
set_is_single_seq(obj,is_single_seq)
	Wise2_GenomicDB * obj
	boolean is_single_seq
	CODE:
	RETVAL = Wise2_replace_is_single_seq_GenomicDB(obj,is_single_seq);
	OUTPUT:
	RETVAL



boolean
is_single_seq(obj)
	Wise2_GenomicDB * obj
	CODE:
	RETVAL = Wise2_access_is_single_seq_GenomicDB(obj);
	OUTPUT:
	RETVAL



boolean
set_done_forward(obj,done_forward)
	Wise2_GenomicDB * obj
	boolean done_forward
	CODE:
	RETVAL = Wise2_replace_done_forward_GenomicDB(obj,done_forward);
	OUTPUT:
	RETVAL



boolean
done_forward(obj)
	Wise2_GenomicDB * obj
	CODE:
	RETVAL = Wise2_access_done_forward_GenomicDB(obj);
	OUTPUT:
	RETVAL



boolean
set_forw(obj,forw)
	Wise2_GenomicDB * obj
	Wise2_ComplexSequence * forw
	CODE:
	RETVAL = Wise2_replace_forw_GenomicDB(obj,Wise2_hard_link_ComplexSequence(forw));
	OUTPUT:
	RETVAL



Wise2_ComplexSequence *
forw(obj)
	Wise2_GenomicDB * obj
	INIT:
Wise2_ComplexSequence * temp;
	CODE:
	temp = Wise2_hard_link_ComplexSequence(Wise2_access_forw_GenomicDB(obj));
	RETVAL = temp;
	OUTPUT:
	RETVAL



boolean
set_rev(obj,rev)
	Wise2_GenomicDB * obj
	Wise2_ComplexSequence * rev
	CODE:
	RETVAL = Wise2_replace_rev_GenomicDB(obj,Wise2_hard_link_ComplexSequence(rev));
	OUTPUT:
	RETVAL



Wise2_ComplexSequence *
rev(obj)
	Wise2_GenomicDB * obj
	INIT:
Wise2_ComplexSequence * temp;
	CODE:
	temp = Wise2_hard_link_ComplexSequence(Wise2_access_rev_GenomicDB(obj));
	RETVAL = temp;
	OUTPUT:
	RETVAL



boolean
set_sdb(obj,sdb)
	Wise2_GenomicDB * obj
	Wise2_SequenceDB * sdb
	CODE:
	RETVAL = Wise2_replace_sdb_GenomicDB(obj,Wise2_hard_link_SequenceDB(sdb));
	OUTPUT:
	RETVAL



Wise2_SequenceDB *
sdb(obj)
	Wise2_GenomicDB * obj
	INIT:
Wise2_SequenceDB * temp;
	CODE:
	temp = Wise2_hard_link_SequenceDB(Wise2_access_sdb_GenomicDB(obj));
	RETVAL = temp;
	OUTPUT:
	RETVAL



boolean
set_current(obj,current)
	Wise2_GenomicDB * obj
	Wise2_Genomic * current
	CODE:
	RETVAL = Wise2_replace_current_GenomicDB(obj,Wise2_hard_link_Genomic(current));
	OUTPUT:
	RETVAL



Wise2_Genomic *
current(obj)
	Wise2_GenomicDB * obj
	INIT:
Wise2_Genomic * temp;
	CODE:
	temp = Wise2_hard_link_Genomic(Wise2_access_current_GenomicDB(obj));
	RETVAL = temp;
	OUTPUT:
	RETVAL



boolean
set_cses(obj,cses)
	Wise2_GenomicDB * obj
	Wise2_ComplexSequenceEvalSet * cses
	CODE:
	RETVAL = Wise2_replace_cses_GenomicDB(obj,Wise2_hard_link_ComplexSequenceEvalSet(cses));
	OUTPUT:
	RETVAL



Wise2_ComplexSequenceEvalSet *
cses(obj)
	Wise2_GenomicDB * obj
	INIT:
Wise2_ComplexSequenceEvalSet * temp;
	CODE:
	temp = Wise2_hard_link_ComplexSequenceEvalSet(Wise2_access_cses_GenomicDB(obj));
	RETVAL = temp;
	OUTPUT:
	RETVAL




Wise2_GenomicDB *
new(class)
	char * class
	PPCODE:
	Wise2_GenomicDB * out;
	out = Wise2_GenomicDB_alloc();
	ST(0) = sv_newmortal();
	sv_setref_pv(ST(0),class,(void*)out);
	XSRETURN(1);

void
DESTROY(obj)
	Wise2_GenomicDB * obj
	CODE:
	Wise2_free_GenomicDB(obj);



MODULE = Wise2 PACKAGE = Wise2

Wise2_GenomicDB *
new_GenomicDB_from_single_seq(gen,cses,score_in_repeat_coding)
	Wise2_Genomic * gen
	Wise2_ComplexSequenceEvalSet * cses
	int score_in_repeat_coding
	CODE:
	RETVAL = Wise2_new_GenomicDB_from_single_seq(gen,cses,score_in_repeat_coding);
	OUTPUT:
	RETVAL



Wise2_GenomicDB *
new_GenomicDB(seqdb,cses,length_of_N,repeat_in_cds_score)
	Wise2_SequenceDB * seqdb
	Wise2_ComplexSequenceEvalSet * cses
	int length_of_N
	int repeat_in_cds_score
	CODE:
	RETVAL = Wise2_new_GenomicDB(seqdb,cses,length_of_N,repeat_in_cds_score);
	OUTPUT:
	RETVAL



