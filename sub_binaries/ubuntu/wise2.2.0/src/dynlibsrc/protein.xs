

MODULE = Wise2 PACKAGE = Wise2::Protein

Wise2_Protein *
Protein_from_Sequence(seq)
	Wise2_Sequence * seq
	CODE:
	RETVAL = Wise2_Protein_from_Sequence(Wise2_hard_link_Sequence(seq));
	OUTPUT:
	RETVAL



Wise2_Protein *
hard_link_Protein(obj)
	Wise2_Protein * obj
	CODE:
	RETVAL = Wise2_hard_link_Protein(obj);
	OUTPUT:
	RETVAL



Wise2_Protein *
alloc()
	CODE:
	RETVAL = Wise2_Protein_alloc();
	OUTPUT:
	RETVAL



boolean
set_baseseq(obj,baseseq)
	Wise2_Protein * obj
	Wise2_Sequence * baseseq
	CODE:
	RETVAL = Wise2_replace_baseseq_Protein(obj,Wise2_hard_link_Sequence(baseseq));
	OUTPUT:
	RETVAL



Wise2_Sequence *
baseseq(obj)
	Wise2_Protein * obj
	INIT:
Wise2_Sequence * temp;
	CODE:
	temp = Wise2_hard_link_Sequence(Wise2_access_baseseq_Protein(obj));
	RETVAL = temp;
	OUTPUT:
	RETVAL




Wise2_Protein *
new(class)
	char * class
	PPCODE:
	Wise2_Protein * out;
	out = Wise2_Protein_alloc();
	ST(0) = sv_newmortal();
	sv_setref_pv(ST(0),class,(void*)out);
	XSRETURN(1);

void
DESTROY(obj)
	Wise2_Protein * obj
	CODE:
	Wise2_free_Protein(obj);



MODULE = Wise2 PACKAGE = Wise2

