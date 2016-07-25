

MODULE = Wise2 PACKAGE = Wise2::Exon

Wise2_Exon *
hard_link_Exon(obj)
	Wise2_Exon * obj
	CODE:
	RETVAL = Wise2_hard_link_Exon(obj);
	OUTPUT:
	RETVAL



Wise2_Exon *
Exon_alloc_std()
	CODE:
	RETVAL = Wise2_Exon_alloc_std();
	OUTPUT:
	RETVAL



boolean
set_start(obj,start)
	Wise2_Exon * obj
	int start
	CODE:
	RETVAL = Wise2_replace_start_Exon(obj,start);
	OUTPUT:
	RETVAL



int
start(obj)
	Wise2_Exon * obj
	CODE:
	RETVAL = Wise2_access_start_Exon(obj);
	OUTPUT:
	RETVAL



boolean
set_end(obj,end)
	Wise2_Exon * obj
	int end
	CODE:
	RETVAL = Wise2_replace_end_Exon(obj,end);
	OUTPUT:
	RETVAL



int
end(obj)
	Wise2_Exon * obj
	CODE:
	RETVAL = Wise2_access_end_Exon(obj);
	OUTPUT:
	RETVAL



boolean
set_used(obj,used)
	Wise2_Exon * obj
	boolean used
	CODE:
	RETVAL = Wise2_replace_used_Exon(obj,used);
	OUTPUT:
	RETVAL



boolean
used(obj)
	Wise2_Exon * obj
	CODE:
	RETVAL = Wise2_access_used_Exon(obj);
	OUTPUT:
	RETVAL



boolean
set_score(obj,score)
	Wise2_Exon * obj
	double score
	CODE:
	RETVAL = Wise2_replace_score_Exon(obj,score);
	OUTPUT:
	RETVAL



double
score(obj)
	Wise2_Exon * obj
	CODE:
	RETVAL = Wise2_access_score_Exon(obj);
	OUTPUT:
	RETVAL



Wise2_SupportingFeature *
sf(obj,i)
	Wise2_Exon * obj
	int i
	INIT:
Wise2_SupportingFeature * temp;
	CODE:
	temp = Wise2_hard_link_SupportingFeature(Wise2_access_sf_Exon(obj,i));
	RETVAL = temp;
	OUTPUT:
	RETVAL



int
length_sf(obj)
	Wise2_Exon * obj
	CODE:
	RETVAL = Wise2_length_sf_Exon(obj);
	OUTPUT:
	RETVAL



int
flush_sf(obj)
	Wise2_Exon * obj
	CODE:
	RETVAL = Wise2_flush_Exon(obj);
	OUTPUT:
	RETVAL



boolean
add_sf(obj,add)
	Wise2_Exon * obj
	Wise2_SupportingFeature * add
	CODE:
	RETVAL = Wise2_add_Exon(obj,Wise2_hard_link_SupportingFeature(add));
	OUTPUT:
	RETVAL



boolean
set_phase(obj,phase)
	Wise2_Exon * obj
	int phase
	CODE:
	RETVAL = Wise2_replace_phase_Exon(obj,phase);
	OUTPUT:
	RETVAL



int
phase(obj)
	Wise2_Exon * obj
	CODE:
	RETVAL = Wise2_access_phase_Exon(obj);
	OUTPUT:
	RETVAL




Wise2_Exon *
new(class)
	char * class
	PPCODE:
	Wise2_Exon * out;
	out = Wise2_Exon_alloc_std();
	ST(0) = sv_newmortal();
	sv_setref_pv(ST(0),class,(void*)out);
	XSRETURN(1);

void
DESTROY(obj)
	Wise2_Exon * obj
	CODE:
	Wise2_free_Exon(obj);

void
each_sf(obj)
	Wise2_Exon * obj
	PPCODE:
	int i=0;
	int len;
	SV* temp;
	len = Wise2_length_sf_Exon(obj);
	for(i=0;i<len;i++){
	  temp = sv_newmortal();
	  sv_setref_pv(temp, "Wise2::SupportingFeature", (void*) (Wise2_hard_link_SupportingFeature(Wise2_access_sf_Exon(obj,i))));
	  XPUSHs(temp);
	  }
	XSRETURN(len);



MODULE = Wise2 PACKAGE = Wise2::Transcript

Wise2_cDNA *
get_cDNA_from_Transcript(trs)
	Wise2_Transcript * trs
	INIT:
Wise2_cDNA * temp;
	CODE:
	temp = Wise2_hard_link_cDNA(Wise2_get_cDNA_from_Transcript(trs));
	RETVAL = temp;
	OUTPUT:
	RETVAL



Wise2_Transcript *
hard_link_Transcript(obj)
	Wise2_Transcript * obj
	CODE:
	RETVAL = Wise2_hard_link_Transcript(obj);
	OUTPUT:
	RETVAL



Wise2_Transcript *
Transcript_alloc_std()
	CODE:
	RETVAL = Wise2_Transcript_alloc_std();
	OUTPUT:
	RETVAL



Wise2_Exon *
exon(obj,i)
	Wise2_Transcript * obj
	int i
	INIT:
Wise2_Exon * temp;
	CODE:
	temp = Wise2_hard_link_Exon(Wise2_access_exon_Transcript(obj,i));
	RETVAL = temp;
	OUTPUT:
	RETVAL



int
length_exon(obj)
	Wise2_Transcript * obj
	CODE:
	RETVAL = Wise2_length_exon_Transcript(obj);
	OUTPUT:
	RETVAL



int
flush_exon(obj)
	Wise2_Transcript * obj
	CODE:
	RETVAL = Wise2_flush_ex_Transcript(obj);
	OUTPUT:
	RETVAL



boolean
add_exon(obj,add)
	Wise2_Transcript * obj
	Wise2_Exon * add
	CODE:
	RETVAL = Wise2_add_ex_Transcript(obj,Wise2_hard_link_Exon(add));
	OUTPUT:
	RETVAL



boolean
set_parent(obj,parent)
	Wise2_Transcript * obj
	Wise2_Gene * parent
	CODE:
	RETVAL = Wise2_replace_parent_Transcript(obj,Wise2_hard_link_Gene(parent));
	OUTPUT:
	RETVAL



Wise2_Gene *
parent(obj)
	Wise2_Transcript * obj
	INIT:
Wise2_Gene * temp;
	CODE:
	temp = Wise2_hard_link_Gene(Wise2_access_parent_Transcript(obj));
	RETVAL = temp;
	OUTPUT:
	RETVAL



Wise2_Translation *
translation(obj,i)
	Wise2_Transcript * obj
	int i
	INIT:
Wise2_Translation * temp;
	CODE:
	temp = Wise2_hard_link_Translation(Wise2_access_translation_Transcript(obj,i));
	RETVAL = temp;
	OUTPUT:
	RETVAL



int
length_translation(obj)
	Wise2_Transcript * obj
	CODE:
	RETVAL = Wise2_length_translation_Transcript(obj);
	OUTPUT:
	RETVAL



int
flush_translation(obj)
	Wise2_Transcript * obj
	CODE:
	RETVAL = Wise2_flush_Transcript(obj);
	OUTPUT:
	RETVAL



boolean
add_translation(obj,add)
	Wise2_Transcript * obj
	Wise2_Translation * add
	CODE:
	RETVAL = Wise2_add_Transcript(obj,Wise2_hard_link_Translation(add));
	OUTPUT:
	RETVAL



boolean
set_cDNA(obj,cDNA)
	Wise2_Transcript * obj
	Wise2_cDNA * cDNA
	CODE:
	RETVAL = Wise2_replace_cDNA_Transcript(obj,Wise2_hard_link_cDNA(cDNA));
	OUTPUT:
	RETVAL



Wise2_cDNA *
cDNA(obj)
	Wise2_Transcript * obj
	INIT:
Wise2_cDNA * temp;
	CODE:
	temp = Wise2_hard_link_cDNA(Wise2_access_cDNA_Transcript(obj));
	RETVAL = temp;
	OUTPUT:
	RETVAL




Wise2_Transcript *
new(class)
	char * class
	PPCODE:
	Wise2_Transcript * out;
	out = Wise2_Transcript_alloc_std();
	ST(0) = sv_newmortal();
	sv_setref_pv(ST(0),class,(void*)out);
	XSRETURN(1);

void
DESTROY(obj)
	Wise2_Transcript * obj
	CODE:
	Wise2_free_Transcript(obj);

void
each_exon(obj)
	Wise2_Transcript * obj
	PPCODE:
	int i=0;
	int len;
	SV* temp;
	len = Wise2_length_exon_Transcript(obj);
	for(i=0;i<len;i++){
	  temp = sv_newmortal();
	  sv_setref_pv(temp, "Wise2::Exon", (void*) (Wise2_hard_link_Exon(Wise2_access_exon_Transcript(obj,i))));
	  XPUSHs(temp);
	  }
	XSRETURN(len);

void
each_translation(obj)
	Wise2_Transcript * obj
	PPCODE:
	int i=0;
	int len;
	SV* temp;
	len = Wise2_length_translation_Transcript(obj);
	for(i=0;i<len;i++){
	  temp = sv_newmortal();
	  sv_setref_pv(temp, "Wise2::Translation", (void*) (Wise2_hard_link_Translation(Wise2_access_translation_Transcript(obj,i))));
	  XPUSHs(temp);
	  }
	XSRETURN(len);



MODULE = Wise2 PACKAGE = Wise2

