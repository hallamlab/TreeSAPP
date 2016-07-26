

MODULE = Wise2 PACKAGE = Wise2::GenomicRegion

Wise2_GenomicRegion *
new_GenomicRegion(gen)
	Wise2_Genomic * gen
	CODE:
	RETVAL = Wise2_new_GenomicRegion(gen);
	OUTPUT:
	RETVAL



Wise2_GenomicRegion *
read_EMBL_GenomicRegion_file(filename)
	char * filename
	CODE:
	RETVAL = Wise2_read_EMBL_GenomicRegion_file(filename);
	OUTPUT:
	RETVAL



Wise2_GenomicRegion *
read_EMBL_GenomicRegion(ifp)
	FILE * ifp
	CODE:
	RETVAL = Wise2_read_EMBL_GenomicRegion(ifp);
	OUTPUT:
	RETVAL



boolean
add_Gene_to_GenomicRegion(gr,gene)
	Wise2_GenomicRegion * gr
	Wise2_Gene * gene
	CODE:
	RETVAL = Wise2_add_Gene_to_GenomicRegion(gr,gene);
	OUTPUT:
	RETVAL



void
show_ace_GenomicRegion(gr,seq_name,ofp)
	Wise2_GenomicRegion * gr
	char * seq_name
	FILE * ofp
	CODE:
	Wise2_show_ace_GenomicRegion(gr,seq_name,ofp);



void
show_pretty_GenomicRegion(gr,show_supporting,ofp)
	Wise2_GenomicRegion * gr
	boolean show_supporting
	FILE * ofp
	CODE:
	Wise2_show_pretty_GenomicRegion(gr,show_supporting,ofp);



void
write_Diana_FT_GenomicRegion(gr,ofp)
	Wise2_GenomicRegion * gr
	FILE * ofp
	CODE:
	Wise2_write_Diana_FT_GenomicRegion(gr,ofp);



void
write_Embl_FT_GenomicRegion(gr,ofp)
	Wise2_GenomicRegion * gr
	FILE * ofp
	CODE:
	Wise2_write_Embl_FT_GenomicRegion(gr,ofp);



Wise2_GenomicRegion *
hard_link_GenomicRegion(obj)
	Wise2_GenomicRegion * obj
	CODE:
	RETVAL = Wise2_hard_link_GenomicRegion(obj);
	OUTPUT:
	RETVAL



Wise2_GenomicRegion *
GenomicRegion_alloc_std()
	CODE:
	RETVAL = Wise2_GenomicRegion_alloc_std();
	OUTPUT:
	RETVAL



Wise2_Gene *
gene(obj,i)
	Wise2_GenomicRegion * obj
	int i
	INIT:
Wise2_Gene * temp;
	CODE:
	temp = Wise2_hard_link_Gene(Wise2_access_gene_GenomicRegion(obj,i));
	RETVAL = temp;
	OUTPUT:
	RETVAL



int
length_gene(obj)
	Wise2_GenomicRegion * obj
	CODE:
	RETVAL = Wise2_length_gene_GenomicRegion(obj);
	OUTPUT:
	RETVAL



int
flush_gene(obj)
	Wise2_GenomicRegion * obj
	CODE:
	RETVAL = Wise2_flush_GenomicRegion(obj);
	OUTPUT:
	RETVAL



boolean
add_gene(obj,add)
	Wise2_GenomicRegion * obj
	Wise2_Gene * add
	CODE:
	RETVAL = Wise2_add_GenomicRegion(obj,Wise2_hard_link_Gene(add));
	OUTPUT:
	RETVAL



boolean
set_genomic(obj,genomic)
	Wise2_GenomicRegion * obj
	Wise2_Genomic * genomic
	CODE:
	RETVAL = Wise2_replace_genomic_GenomicRegion(obj,Wise2_hard_link_Genomic(genomic));
	OUTPUT:
	RETVAL



Wise2_Genomic *
genomic(obj)
	Wise2_GenomicRegion * obj
	INIT:
Wise2_Genomic * temp;
	CODE:
	temp = Wise2_hard_link_Genomic(Wise2_access_genomic_GenomicRegion(obj));
	RETVAL = temp;
	OUTPUT:
	RETVAL




Wise2_GenomicRegion *
new(class)
	char * class
	PPCODE:
	Wise2_GenomicRegion * out;
	out = Wise2_GenomicRegion_alloc_std();
	ST(0) = sv_newmortal();
	sv_setref_pv(ST(0),class,(void*)out);
	XSRETURN(1);

void
DESTROY(obj)
	Wise2_GenomicRegion * obj
	CODE:
	Wise2_free_GenomicRegion(obj);

void
each_gene(obj)
	Wise2_GenomicRegion * obj
	PPCODE:
	int i=0;
	int len;
	SV* temp;
	len = Wise2_length_gene_GenomicRegion(obj);
	for(i=0;i<len;i++){
	  temp = sv_newmortal();
	  sv_setref_pv(temp, "Wise2::Gene", (void*) (Wise2_hard_link_Gene(Wise2_access_gene_GenomicRegion(obj,i))));
	  XPUSHs(temp);
	  }
	XSRETURN(len);



MODULE = Wise2 PACKAGE = Wise2

