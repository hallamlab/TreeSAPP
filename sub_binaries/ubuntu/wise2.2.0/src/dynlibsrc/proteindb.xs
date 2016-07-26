

MODULE = Wise2 PACKAGE = Wise2::ProteinDB

Wise2_ProteinDB *
hard_link_ProteinDB(obj)
	Wise2_ProteinDB * obj
	CODE:
	RETVAL = Wise2_hard_link_ProteinDB(obj);
	OUTPUT:
	RETVAL



Wise2_ProteinDB *
alloc()
	CODE:
	RETVAL = Wise2_ProteinDB_alloc();
	OUTPUT:
	RETVAL



boolean
set_is_single_seq(obj,is_single_seq)
	Wise2_ProteinDB * obj
	boolean is_single_seq
	CODE:
	RETVAL = Wise2_replace_is_single_seq_ProteinDB(obj,is_single_seq);
	OUTPUT:
	RETVAL



boolean
is_single_seq(obj)
	Wise2_ProteinDB * obj
	CODE:
	RETVAL = Wise2_access_is_single_seq_ProteinDB(obj);
	OUTPUT:
	RETVAL



boolean
set_is_random_db(obj,is_random_db)
	Wise2_ProteinDB * obj
	boolean is_random_db
	CODE:
	RETVAL = Wise2_replace_is_random_db_ProteinDB(obj,is_random_db);
	OUTPUT:
	RETVAL



boolean
is_random_db(obj)
	Wise2_ProteinDB * obj
	CODE:
	RETVAL = Wise2_access_is_random_db_ProteinDB(obj);
	OUTPUT:
	RETVAL



boolean
set_single(obj,single)
	Wise2_ProteinDB * obj
	Wise2_ComplexSequence * single
	CODE:
	RETVAL = Wise2_replace_single_ProteinDB(obj,Wise2_hard_link_ComplexSequence(single));
	OUTPUT:
	RETVAL



Wise2_ComplexSequence *
single(obj)
	Wise2_ProteinDB * obj
	INIT:
Wise2_ComplexSequence * temp;
	CODE:
	temp = Wise2_hard_link_ComplexSequence(Wise2_access_single_ProteinDB(obj));
	RETVAL = temp;
	OUTPUT:
	RETVAL



boolean
set_sdb(obj,sdb)
	Wise2_ProteinDB * obj
	Wise2_SequenceDB * sdb
	CODE:
	RETVAL = Wise2_replace_sdb_ProteinDB(obj,Wise2_hard_link_SequenceDB(sdb));
	OUTPUT:
	RETVAL



Wise2_SequenceDB *
sdb(obj)
	Wise2_ProteinDB * obj
	INIT:
Wise2_SequenceDB * temp;
	CODE:
	temp = Wise2_hard_link_SequenceDB(Wise2_access_sdb_ProteinDB(obj));
	RETVAL = temp;
	OUTPUT:
	RETVAL



boolean
set_cses(obj,cses)
	Wise2_ProteinDB * obj
	Wise2_ComplexSequenceEvalSet * cses
	CODE:
	RETVAL = Wise2_replace_cses_ProteinDB(obj,Wise2_hard_link_ComplexSequenceEvalSet(cses));
	OUTPUT:
	RETVAL



Wise2_ComplexSequenceEvalSet *
cses(obj)
	Wise2_ProteinDB * obj
	INIT:
Wise2_ComplexSequenceEvalSet * temp;
	CODE:
	temp = Wise2_hard_link_ComplexSequenceEvalSet(Wise2_access_cses_ProteinDB(obj));
	RETVAL = temp;
	OUTPUT:
	RETVAL



boolean
set_rnd(obj,rnd)
	Wise2_ProteinDB * obj
	Wise2_RandomProteinDB * rnd
	CODE:
	RETVAL = Wise2_replace_rnd_ProteinDB(obj,Wise2_hard_link_RandomProteinDB(rnd));
	OUTPUT:
	RETVAL



Wise2_RandomProteinDB *
rnd(obj)
	Wise2_ProteinDB * obj
	INIT:
Wise2_RandomProteinDB * temp;
	CODE:
	temp = Wise2_hard_link_RandomProteinDB(Wise2_access_rnd_ProteinDB(obj));
	RETVAL = temp;
	OUTPUT:
	RETVAL



boolean
set_test_dna(obj,test_dna)
	Wise2_ProteinDB * obj
	boolean test_dna
	CODE:
	RETVAL = Wise2_replace_test_dna_ProteinDB(obj,test_dna);
	OUTPUT:
	RETVAL



boolean
test_dna(obj)
	Wise2_ProteinDB * obj
	CODE:
	RETVAL = Wise2_access_test_dna_ProteinDB(obj);
	OUTPUT:
	RETVAL




Wise2_ProteinDB *
new(class)
	char * class
	PPCODE:
	Wise2_ProteinDB * out;
	out = Wise2_ProteinDB_alloc();
	ST(0) = sv_newmortal();
	sv_setref_pv(ST(0),class,(void*)out);
	XSRETURN(1);

void
DESTROY(obj)
	Wise2_ProteinDB * obj
	CODE:
	Wise2_free_ProteinDB(obj);



MODULE = Wise2 PACKAGE = Wise2

Wise2_ProteinDB *
new_ProteinDB_from_single_seq(seq)
	Wise2_Sequence * seq
	CODE:
	RETVAL = Wise2_new_ProteinDB_from_single_seq(seq);
	OUTPUT:
	RETVAL



Wise2_ProteinDB *
single_fasta_ProteinDB(filename)
	char * filename
	CODE:
	RETVAL = Wise2_single_fasta_ProteinDB(filename);
	OUTPUT:
	RETVAL



Wise2_ProteinDB *
new_ProteinDB(seqdb,cses)
	Wise2_SequenceDB * seqdb
	Wise2_ComplexSequenceEvalSet * cses
	CODE:
	RETVAL = Wise2_new_ProteinDB(seqdb,cses);
	OUTPUT:
	RETVAL



