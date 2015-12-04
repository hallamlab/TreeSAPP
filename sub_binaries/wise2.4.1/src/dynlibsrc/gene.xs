

MODULE = Wise2 PACKAGE = Wise2::Gene

Wise2_Genomic *
get_Genomic_from_Gene(gene)
	Wise2_Gene * gene
	INIT:
Wise2_Genomic * temp;
	CODE:
	temp = Wise2_hard_link_Genomic(Wise2_get_Genomic_from_Gene(gene));
	RETVAL = temp;
	OUTPUT:
	RETVAL



void
show_pretty_Gene(ge,show_supporting,ofp)
	Wise2_Gene * ge
	boolean show_supporting
	FILE * ofp
	CODE:
	Wise2_show_pretty_Gene(ge,show_supporting,ofp);



Wise2_Gene *
hard_link_Gene(obj)
	Wise2_Gene * obj
	CODE:
	RETVAL = Wise2_hard_link_Gene(obj);
	OUTPUT:
	RETVAL



Wise2_Gene *
Gene_alloc_std()
	CODE:
	RETVAL = Wise2_Gene_alloc_std();
	OUTPUT:
	RETVAL



boolean
set_start(obj,start)
	Wise2_Gene * obj
	int start
	CODE:
	RETVAL = Wise2_replace_start_Gene(obj,start);
	OUTPUT:
	RETVAL



int
start(obj)
	Wise2_Gene * obj
	CODE:
	RETVAL = Wise2_access_start_Gene(obj);
	OUTPUT:
	RETVAL



boolean
set_end(obj,end)
	Wise2_Gene * obj
	int end
	CODE:
	RETVAL = Wise2_replace_end_Gene(obj,end);
	OUTPUT:
	RETVAL



int
end(obj)
	Wise2_Gene * obj
	CODE:
	RETVAL = Wise2_access_end_Gene(obj);
	OUTPUT:
	RETVAL



boolean
set_parent(obj,parent)
	Wise2_Gene * obj
	Wise2_GenomicRegion * parent
	CODE:
	RETVAL = Wise2_replace_parent_Gene(obj,Wise2_hard_link_GenomicRegion(parent));
	OUTPUT:
	RETVAL



Wise2_GenomicRegion *
parent(obj)
	Wise2_Gene * obj
	INIT:
Wise2_GenomicRegion * temp;
	CODE:
	temp = Wise2_hard_link_GenomicRegion(Wise2_access_parent_Gene(obj));
	RETVAL = temp;
	OUTPUT:
	RETVAL



boolean
set_genomic(obj,genomic)
	Wise2_Gene * obj
	Wise2_Genomic * genomic
	CODE:
	RETVAL = Wise2_replace_genomic_Gene(obj,Wise2_hard_link_Genomic(genomic));
	OUTPUT:
	RETVAL



Wise2_Genomic *
genomic(obj)
	Wise2_Gene * obj
	INIT:
Wise2_Genomic * temp;
	CODE:
	temp = Wise2_hard_link_Genomic(Wise2_access_genomic_Gene(obj));
	RETVAL = temp;
	OUTPUT:
	RETVAL



Wise2_Transcript *
transcript(obj,i)
	Wise2_Gene * obj
	int i
	INIT:
Wise2_Transcript * temp;
	CODE:
	temp = Wise2_hard_link_Transcript(Wise2_access_transcript_Gene(obj,i));
	RETVAL = temp;
	OUTPUT:
	RETVAL



int
length_transcript(obj)
	Wise2_Gene * obj
	CODE:
	RETVAL = Wise2_length_transcript_Gene(obj);
	OUTPUT:
	RETVAL



int
flush_transcript(obj)
	Wise2_Gene * obj
	CODE:
	RETVAL = Wise2_flush_Gene(obj);
	OUTPUT:
	RETVAL



boolean
add_transcript(obj,add)
	Wise2_Gene * obj
	Wise2_Transcript * add
	CODE:
	RETVAL = Wise2_add_Gene(obj,Wise2_hard_link_Transcript(add));
	OUTPUT:
	RETVAL



boolean
set_name(obj,name)
	Wise2_Gene * obj
	char * name
	CODE:
	RETVAL = Wise2_replace_name_Gene(obj,Wise2_stringalloc(name));
	OUTPUT:
	RETVAL



char *
name(obj)
	Wise2_Gene * obj
	INIT:
	char * temp;
	CODE:
	temp = Wise2_stringalloc(Wise2_access_name_Gene(obj));
	RETVAL = temp;
	OUTPUT:
	RETVAL



boolean
set_bits(obj,bits)
	Wise2_Gene * obj
	double bits
	CODE:
	RETVAL = Wise2_replace_bits_Gene(obj,bits);
	OUTPUT:
	RETVAL



double
bits(obj)
	Wise2_Gene * obj
	CODE:
	RETVAL = Wise2_access_bits_Gene(obj);
	OUTPUT:
	RETVAL



boolean
set_seqname(obj,seqname)
	Wise2_Gene * obj
	char * seqname
	CODE:
	RETVAL = Wise2_replace_seqname_Gene(obj,Wise2_stringalloc(seqname));
	OUTPUT:
	RETVAL



char *
seqname(obj)
	Wise2_Gene * obj
	INIT:
	char * temp;
	CODE:
	temp = Wise2_stringalloc(Wise2_access_seqname_Gene(obj));
	RETVAL = temp;
	OUTPUT:
	RETVAL



boolean
set_ispseudo(obj,ispseudo)
	Wise2_Gene * obj
	boolean ispseudo
	CODE:
	RETVAL = Wise2_replace_ispseudo_Gene(obj,ispseudo);
	OUTPUT:
	RETVAL



boolean
ispseudo(obj)
	Wise2_Gene * obj
	CODE:
	RETVAL = Wise2_access_ispseudo_Gene(obj);
	OUTPUT:
	RETVAL




Wise2_Gene *
new(class)
	char * class
	PPCODE:
	Wise2_Gene * out;
	out = Wise2_Gene_alloc_std();
	ST(0) = sv_newmortal();
	sv_setref_pv(ST(0),class,(void*)out);
	XSRETURN(1);

void
DESTROY(obj)
	Wise2_Gene * obj
	CODE:
	Wise2_free_Gene(obj);

void
each_transcript(obj)
	Wise2_Gene * obj
	PPCODE:
	int i=0;
	int len;
	SV* temp;
	len = Wise2_length_transcript_Gene(obj);
	for(i=0;i<len;i++){
	  temp = sv_newmortal();
	  sv_setref_pv(temp, "Wise2::Transcript", (void*) (Wise2_hard_link_Transcript(Wise2_access_transcript_Gene(obj,i))));
	  XPUSHs(temp);
	  }
	XSRETURN(len);



MODULE = Wise2 PACKAGE = Wise2

