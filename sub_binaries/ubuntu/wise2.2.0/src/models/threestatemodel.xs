

MODULE = Wise2 PACKAGE = Wise2::ThreeStateModel

void
force_global_model(tsm,prob_into_model)
	Wise2_ThreeStateModel * tsm
	double prob_into_model
	CODE:
	Wise2_force_global_model(tsm,prob_into_model);



void
force_weighted_local_model(tsm,prob_into_model,ratio_start,ratio_end)
	Wise2_ThreeStateModel * tsm
	double prob_into_model
	double ratio_start
	double ratio_end
	CODE:
	Wise2_force_weighted_local_model(tsm,prob_into_model,ratio_start,ratio_end);



Wise2_ThreeStateModel *
ThreeStateModel_from_half_bit_Sequence(pro,mat,rm,gap,ext)
	Wise2_Protein * pro
	Wise2_CompMat * mat
	Wise2_RandomModel * rm
	int gap
	int ext
	CODE:
	RETVAL = Wise2_ThreeStateModel_from_half_bit_Sequence(pro,mat,rm,gap,ext);
	OUTPUT:
	RETVAL



void
write_HMMer_1_7_ascii_ThreeStateModel(tsm,ofp)
	Wise2_ThreeStateModel * tsm
	FILE * ofp
	CODE:
	Wise2_write_HMMer_1_7_ascii_ThreeStateModel(tsm,ofp);



Wise2_ThreeStateModel *
hard_link_ThreeStateModel(obj)
	Wise2_ThreeStateModel * obj
	CODE:
	RETVAL = Wise2_hard_link_ThreeStateModel(obj);
	OUTPUT:
	RETVAL



Wise2_ThreeStateModel *
ThreeStateModel_alloc_std()
	CODE:
	RETVAL = Wise2_ThreeStateModel_alloc_std();
	OUTPUT:
	RETVAL



boolean
set_name(obj,name)
	Wise2_ThreeStateModel * obj
	char * name
	CODE:
	RETVAL = Wise2_replace_name_ThreeStateModel(obj,Wise2_stringalloc(name));
	OUTPUT:
	RETVAL



char *
name(obj)
	Wise2_ThreeStateModel * obj
	INIT:
	char * temp;
	CODE:
	temp = Wise2_stringalloc(Wise2_access_name_ThreeStateModel(obj));
	RETVAL = temp;
	OUTPUT:
	RETVAL



Wise2_ThreeStateUnit *
unit(obj,i)
	Wise2_ThreeStateModel * obj
	int i
	INIT:
Wise2_ThreeStateUnit * temp;
	CODE:
	temp = Wise2_hard_link_ThreeStateUnit(Wise2_access_unit_ThreeStateModel(obj,i));
	RETVAL = temp;
	OUTPUT:
	RETVAL



int
length_unit(obj)
	Wise2_ThreeStateModel * obj
	CODE:
	RETVAL = Wise2_length_unit_ThreeStateModel(obj);
	OUTPUT:
	RETVAL



int
flush_unit(obj)
	Wise2_ThreeStateModel * obj
	CODE:
	RETVAL = Wise2_flush_ThreeStateModel(obj);
	OUTPUT:
	RETVAL



boolean
add_unit(obj,add)
	Wise2_ThreeStateModel * obj
	Wise2_ThreeStateUnit * add
	CODE:
	RETVAL = Wise2_add_ThreeStateModel(obj,Wise2_hard_link_ThreeStateUnit(add));
	OUTPUT:
	RETVAL



boolean
set_alphabet(obj,alphabet)
	Wise2_ThreeStateModel * obj
	char * alphabet
	CODE:
	RETVAL = Wise2_replace_alphabet_ThreeStateModel(obj,Wise2_stringalloc(alphabet));
	OUTPUT:
	RETVAL



char *
alphabet(obj)
	Wise2_ThreeStateModel * obj
	INIT:
	char * temp;
	CODE:
	temp = Wise2_stringalloc(Wise2_access_alphabet_ThreeStateModel(obj));
	RETVAL = temp;
	OUTPUT:
	RETVAL



boolean
set_accession(obj,accession)
	Wise2_ThreeStateModel * obj
	char * accession
	CODE:
	RETVAL = Wise2_replace_accession_ThreeStateModel(obj,Wise2_stringalloc(accession));
	OUTPUT:
	RETVAL



char *
accession(obj)
	Wise2_ThreeStateModel * obj
	INIT:
	char * temp;
	CODE:
	temp = Wise2_stringalloc(Wise2_access_accession_ThreeStateModel(obj));
	RETVAL = temp;
	OUTPUT:
	RETVAL



boolean
set_threshold(obj,threshold)
	Wise2_ThreeStateModel * obj
	double threshold
	CODE:
	RETVAL = Wise2_replace_threshold_ThreeStateModel(obj,threshold);
	OUTPUT:
	RETVAL



double
threshold(obj)
	Wise2_ThreeStateModel * obj
	CODE:
	RETVAL = Wise2_access_threshold_ThreeStateModel(obj);
	OUTPUT:
	RETVAL



boolean
set_rm(obj,rm)
	Wise2_ThreeStateModel * obj
	Wise2_RandomModel * rm
	CODE:
	RETVAL = Wise2_replace_rm_ThreeStateModel(obj,Wise2_hard_link_RandomModel(rm));
	OUTPUT:
	RETVAL



Wise2_RandomModel *
rm(obj)
	Wise2_ThreeStateModel * obj
	INIT:
Wise2_RandomModel * temp;
	CODE:
	temp = Wise2_hard_link_RandomModel(Wise2_access_rm_ThreeStateModel(obj));
	RETVAL = temp;
	OUTPUT:
	RETVAL




Wise2_ThreeStateModel *
new(class)
	char * class
	PPCODE:
	Wise2_ThreeStateModel * out;
	out = Wise2_ThreeStateModel_alloc_std();
	ST(0) = sv_newmortal();
	sv_setref_pv(ST(0),class,(void*)out);
	XSRETURN(1);

void
DESTROY(obj)
	Wise2_ThreeStateModel * obj
	CODE:
	Wise2_free_ThreeStateModel(obj);

void
each_unit(obj)
	Wise2_ThreeStateModel * obj
	PPCODE:
	int i=0;
	int len;
	SV* temp;
	len = Wise2_length_unit_ThreeStateModel(obj);
	for(i=0;i<len;i++){
	  temp = sv_newmortal();
	  sv_setref_pv(temp, "Wise2::ThreeStateUnit", (void*) (Wise2_hard_link_ThreeStateUnit(Wise2_access_unit_ThreeStateModel(obj,i))));
	  XPUSHs(temp);
	  }
	XSRETURN(len);



MODULE = Wise2 PACKAGE = Wise2::ThreeStateUnit

Wise2_ThreeStateUnit *
hard_link_ThreeStateUnit(obj)
	Wise2_ThreeStateUnit * obj
	CODE:
	RETVAL = Wise2_hard_link_ThreeStateUnit(obj);
	OUTPUT:
	RETVAL



Wise2_ThreeStateUnit *
alloc()
	CODE:
	RETVAL = Wise2_ThreeStateUnit_alloc();
	OUTPUT:
	RETVAL



boolean
set_display_char(obj,display_char)
	Wise2_ThreeStateUnit * obj
	char display_char
	CODE:
	RETVAL = Wise2_replace_display_char_ThreeStateUnit(obj,display_char);
	OUTPUT:
	RETVAL



char
display_char(obj)
	Wise2_ThreeStateUnit * obj
	CODE:
	RETVAL = Wise2_access_display_char_ThreeStateUnit(obj);
	OUTPUT:
	RETVAL




Wise2_ThreeStateUnit *
new(class)
	char * class
	PPCODE:
	Wise2_ThreeStateUnit * out;
	out = Wise2_ThreeStateUnit_alloc();
	ST(0) = sv_newmortal();
	sv_setref_pv(ST(0),class,(void*)out);
	XSRETURN(1);

void
DESTROY(obj)
	Wise2_ThreeStateUnit * obj
	CODE:
	Wise2_free_ThreeStateUnit(obj);



MODULE = Wise2 PACKAGE = Wise2

Wise2_ThreeStateModel *
read_HMMer_1_7_ascii_file(filename)
	char * filename
	CODE:
	RETVAL = Wise2_read_HMMer_1_7_ascii_file(filename);
	OUTPUT:
	RETVAL



Wise2_ThreeStateModel *
read_HMMer_1_7_ascii(ifp)
	FILE * ifp
	CODE:
	RETVAL = Wise2_read_HMMer_1_7_ascii(ifp);
	OUTPUT:
	RETVAL



