

MODULE = Wise2 PACKAGE = Wise2::cDNA

Wise2_cDNA *
truncate_cDNA(cdna,start,stop)
	Wise2_cDNA * cdna
	int start
	int stop
	CODE:
	RETVAL = Wise2_truncate_cDNA(cdna,start,stop);
	OUTPUT:
	RETVAL



Wise2_cDNA *
read_fasta_file_cDNA(filename)
	char * filename
	CODE:
	RETVAL = Wise2_read_fasta_file_cDNA(filename);
	OUTPUT:
	RETVAL



char *
cDNA_name(cdna)
	Wise2_cDNA * cdna
	CODE:
	RETVAL = Wise2_cDNA_name(cdna);
	OUTPUT:
	RETVAL



int
cDNA_length(cdna)
	Wise2_cDNA * cdna
	CODE:
	RETVAL = Wise2_cDNA_length(cdna);
	OUTPUT:
	RETVAL



char
cDNA_seqchar(cdna,pos)
	Wise2_cDNA * cdna
	int pos
	CODE:
	RETVAL = Wise2_cDNA_seqchar(cdna,pos);
	OUTPUT:
	RETVAL



Wise2_cDNA *
cDNA_from_Sequence(seq)
	Wise2_Sequence * seq
	CODE:
	RETVAL = Wise2_cDNA_from_Sequence(Wise2_hard_link_Sequence(seq));
	OUTPUT:
	RETVAL



Wise2_cDNA *
hard_link_cDNA(obj)
	Wise2_cDNA * obj
	CODE:
	RETVAL = Wise2_hard_link_cDNA(obj);
	OUTPUT:
	RETVAL



Wise2_cDNA *
alloc()
	CODE:
	RETVAL = Wise2_cDNA_alloc();
	OUTPUT:
	RETVAL



boolean
set_baseseq(obj,baseseq)
	Wise2_cDNA * obj
	Wise2_Sequence * baseseq
	CODE:
	RETVAL = Wise2_replace_baseseq_cDNA(obj,Wise2_hard_link_Sequence(baseseq));
	OUTPUT:
	RETVAL



Wise2_Sequence *
baseseq(obj)
	Wise2_cDNA * obj
	INIT:
Wise2_Sequence * temp;
	CODE:
	temp = Wise2_hard_link_Sequence(Wise2_access_baseseq_cDNA(obj));
	RETVAL = temp;
	OUTPUT:
	RETVAL




Wise2_cDNA *
new(class)
	char * class
	PPCODE:
	Wise2_cDNA * out;
	out = Wise2_cDNA_alloc();
	ST(0) = sv_newmortal();
	sv_setref_pv(ST(0),class,(void*)out);
	XSRETURN(1);

void
DESTROY(obj)
	Wise2_cDNA * obj
	CODE:
	Wise2_free_cDNA(obj);



MODULE = Wise2 PACKAGE = Wise2

