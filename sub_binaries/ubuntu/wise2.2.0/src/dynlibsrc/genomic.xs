

MODULE = Wise2 PACKAGE = Wise2::Genomic

Wise2_Genomic *
truncate_Genomic(gen,start,stop)
	Wise2_Genomic * gen
	int start
	int stop
	CODE:
	RETVAL = Wise2_truncate_Genomic(gen,start,stop);
	OUTPUT:
	RETVAL



char *
Genomic_name(gen)
	Wise2_Genomic * gen
	CODE:
	RETVAL = Wise2_Genomic_name(gen);
	OUTPUT:
	RETVAL



int
Genomic_length(gen)
	Wise2_Genomic * gen
	CODE:
	RETVAL = Wise2_Genomic_length(gen);
	OUTPUT:
	RETVAL



char
Genomic_seqchar(gen,pos)
	Wise2_Genomic * gen
	int pos
	CODE:
	RETVAL = Wise2_Genomic_seqchar(gen,pos);
	OUTPUT:
	RETVAL



void
show(gen,ofp)
	Wise2_Genomic * gen
	FILE * ofp
	CODE:
	Wise2_show_Genomic(gen,ofp);



Wise2_Genomic *
hard_link_Genomic(obj)
	Wise2_Genomic * obj
	CODE:
	RETVAL = Wise2_hard_link_Genomic(obj);
	OUTPUT:
	RETVAL



Wise2_Genomic *
Genomic_alloc_std()
	CODE:
	RETVAL = Wise2_Genomic_alloc_std();
	OUTPUT:
	RETVAL



boolean
set_baseseq(obj,baseseq)
	Wise2_Genomic * obj
	Wise2_Sequence * baseseq
	CODE:
	RETVAL = Wise2_replace_baseseq_Genomic(obj,Wise2_hard_link_Sequence(baseseq));
	OUTPUT:
	RETVAL



Wise2_Sequence *
baseseq(obj)
	Wise2_Genomic * obj
	INIT:
Wise2_Sequence * temp;
	CODE:
	temp = Wise2_hard_link_Sequence(Wise2_access_baseseq_Genomic(obj));
	RETVAL = temp;
	OUTPUT:
	RETVAL



Wise2_GenomicRepeat *
repeat(obj,i)
	Wise2_Genomic * obj
	int i
	INIT:
Wise2_GenomicRepeat * temp;
	CODE:
	temp = Wise2_hard_link_GenomicRepeat(Wise2_access_repeat_Genomic(obj,i));
	RETVAL = temp;
	OUTPUT:
	RETVAL



int
length_repeat(obj)
	Wise2_Genomic * obj
	CODE:
	RETVAL = Wise2_length_repeat_Genomic(obj);
	OUTPUT:
	RETVAL



int
flush_repeat(obj)
	Wise2_Genomic * obj
	CODE:
	RETVAL = Wise2_flush_Genomic(obj);
	OUTPUT:
	RETVAL



boolean
add_repeat(obj,add)
	Wise2_Genomic * obj
	Wise2_GenomicRepeat * add
	CODE:
	RETVAL = Wise2_add_Genomic(obj,Wise2_hard_link_GenomicRepeat(add));
	OUTPUT:
	RETVAL




Wise2_Genomic *
new(class)
	char * class
	PPCODE:
	Wise2_Genomic * out;
	out = Wise2_Genomic_alloc_std();
	ST(0) = sv_newmortal();
	sv_setref_pv(ST(0),class,(void*)out);
	XSRETURN(1);

void
DESTROY(obj)
	Wise2_Genomic * obj
	CODE:
	Wise2_free_Genomic(obj);

void
each_repeat(obj)
	Wise2_Genomic * obj
	PPCODE:
	int i=0;
	int len;
	SV* temp;
	len = Wise2_length_repeat_Genomic(obj);
	for(i=0;i<len;i++){
	  temp = sv_newmortal();
	  sv_setref_pv(temp, "Wise2::GenomicRepeat", (void*) (Wise2_hard_link_GenomicRepeat(Wise2_access_repeat_Genomic(obj,i))));
	  XPUSHs(temp);
	  }
	XSRETURN(len);



MODULE = Wise2 PACKAGE = Wise2::GenomicRepeat

Wise2_GenomicRepeat *
hard_link_GenomicRepeat(obj)
	Wise2_GenomicRepeat * obj
	CODE:
	RETVAL = Wise2_hard_link_GenomicRepeat(obj);
	OUTPUT:
	RETVAL



Wise2_GenomicRepeat *
alloc()
	CODE:
	RETVAL = Wise2_GenomicRepeat_alloc();
	OUTPUT:
	RETVAL



boolean
set_start(obj,start)
	Wise2_GenomicRepeat * obj
	int start
	CODE:
	RETVAL = Wise2_replace_start_GenomicRepeat(obj,start);
	OUTPUT:
	RETVAL



int
start(obj)
	Wise2_GenomicRepeat * obj
	CODE:
	RETVAL = Wise2_access_start_GenomicRepeat(obj);
	OUTPUT:
	RETVAL



boolean
set_end(obj,end)
	Wise2_GenomicRepeat * obj
	int end
	CODE:
	RETVAL = Wise2_replace_end_GenomicRepeat(obj,end);
	OUTPUT:
	RETVAL



int
end(obj)
	Wise2_GenomicRepeat * obj
	CODE:
	RETVAL = Wise2_access_end_GenomicRepeat(obj);
	OUTPUT:
	RETVAL



boolean
set_type(obj,type)
	Wise2_GenomicRepeat * obj
	char * type
	CODE:
	RETVAL = Wise2_replace_type_GenomicRepeat(obj,Wise2_stringalloc(type));
	OUTPUT:
	RETVAL



char *
type(obj)
	Wise2_GenomicRepeat * obj
	INIT:
	char * temp;
	CODE:
	temp = Wise2_stringalloc(Wise2_access_type_GenomicRepeat(obj));
	RETVAL = temp;
	OUTPUT:
	RETVAL




Wise2_GenomicRepeat *
new(class)
	char * class
	PPCODE:
	Wise2_GenomicRepeat * out;
	out = Wise2_GenomicRepeat_alloc();
	ST(0) = sv_newmortal();
	sv_setref_pv(ST(0),class,(void*)out);
	XSRETURN(1);

void
DESTROY(obj)
	Wise2_GenomicRepeat * obj
	CODE:
	Wise2_free_GenomicRepeat(obj);



MODULE = Wise2 PACKAGE = Wise2

Wise2_Genomic *
read_fasta_file_Genomic(filename,length_of_N)
	char * filename
	int length_of_N
	CODE:
	RETVAL = Wise2_read_fasta_file_Genomic(filename,length_of_N);
	OUTPUT:
	RETVAL



Wise2_Genomic *
read_fasta_Genomic(ifp,length_of_N)
	FILE * ifp
	int length_of_N
	CODE:
	RETVAL = Wise2_read_fasta_Genomic(ifp,length_of_N);
	OUTPUT:
	RETVAL



Wise2_Genomic *
Genomic_from_Sequence_Nheuristic(seq,length_of_N)
	Wise2_Sequence * seq
	int length_of_N
	CODE:
	RETVAL = Wise2_Genomic_from_Sequence_Nheuristic(seq,length_of_N);
	OUTPUT:
	RETVAL



Wise2_Genomic *
Genomic_from_Sequence(seq)
	Wise2_Sequence * seq
	CODE:
	RETVAL = Wise2_Genomic_from_Sequence(Wise2_hard_link_Sequence(seq));
	OUTPUT:
	RETVAL



