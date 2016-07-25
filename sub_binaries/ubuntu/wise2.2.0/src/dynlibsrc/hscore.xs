

MODULE = Wise2 PACKAGE = Wise2::Hscore

int
minimum_score_Hscore(hs)
	Wise2_Hscore * hs
	CODE:
	RETVAL = Wise2_minimum_score_Hscore(hs);
	OUTPUT:
	RETVAL



int
maximum_score_Hscore(hs)
	Wise2_Hscore * hs
	CODE:
	RETVAL = Wise2_maximum_score_Hscore(hs);
	OUTPUT:
	RETVAL



void
sort_Hscore_by_score(hs)
	Wise2_Hscore * hs
	CODE:
	Wise2_sort_Hscore_by_score(hs);



int
length(obj)
	Wise2_Hscore * obj
	CODE:
	RETVAL = Wise2_length_datascore_Hscore(obj);
	OUTPUT:
	RETVAL



Wise2_DataScore *
datascore(hs,i)
	Wise2_Hscore * hs
	int i
	CODE:
	RETVAL = Wise2_get_datascore_Hscore(hs,i);
	OUTPUT:
	RETVAL



int
score(hs,i)
	Wise2_Hscore * hs
	int i
	CODE:
	RETVAL = Wise2_get_score_Hscore(hs,i);
	OUTPUT:
	RETVAL



double
evalue(hs,i)
	Wise2_Hscore * hs
	int i
	CODE:
	RETVAL = Wise2_get_evalue_Hscore(hs,i);
	OUTPUT:
	RETVAL



void
show(hs,ofp)
	Wise2_Hscore * hs
	FILE * ofp
	CODE:
	Wise2_basic_show_Hscore(hs,ofp);



Wise2_Hscore *
hard_link_Hscore(obj)
	Wise2_Hscore * obj
	CODE:
	RETVAL = Wise2_hard_link_Hscore(obj);
	OUTPUT:
	RETVAL



Wise2_Hscore *
Hscore_alloc_std()
	CODE:
	RETVAL = Wise2_Hscore_alloc_std();
	OUTPUT:
	RETVAL




Wise2_Hscore *
new(class)
	char * class
	PPCODE:
	Wise2_Hscore * out;
	out = Wise2_Hscore_alloc_std();
	ST(0) = sv_newmortal();
	sv_setref_pv(ST(0),class,(void*)out);
	XSRETURN(1);

void
DESTROY(obj)
	Wise2_Hscore * obj
	CODE:
	Wise2_free_Hscore(obj);



MODULE = Wise2 PACKAGE = Wise2::DataScore

Wise2_DataScore *
hard_link_DataScore(obj)
	Wise2_DataScore * obj
	CODE:
	RETVAL = Wise2_hard_link_DataScore(obj);
	OUTPUT:
	RETVAL



Wise2_DataScore *
alloc()
	CODE:
	RETVAL = Wise2_DataScore_alloc();
	OUTPUT:
	RETVAL



boolean
set_query(obj,query)
	Wise2_DataScore * obj
	Wise2_DataEntry * query
	CODE:
	RETVAL = Wise2_replace_query_DataScore(obj,Wise2_hard_link_DataEntry(query));
	OUTPUT:
	RETVAL



Wise2_DataEntry *
query(obj)
	Wise2_DataScore * obj
	INIT:
Wise2_DataEntry * temp;
	CODE:
	temp = Wise2_hard_link_DataEntry(Wise2_access_query_DataScore(obj));
	RETVAL = temp;
	OUTPUT:
	RETVAL



boolean
set_target(obj,target)
	Wise2_DataScore * obj
	Wise2_DataEntry * target
	CODE:
	RETVAL = Wise2_replace_target_DataScore(obj,Wise2_hard_link_DataEntry(target));
	OUTPUT:
	RETVAL



Wise2_DataEntry *
target(obj)
	Wise2_DataScore * obj
	INIT:
Wise2_DataEntry * temp;
	CODE:
	temp = Wise2_hard_link_DataEntry(Wise2_access_target_DataScore(obj));
	RETVAL = temp;
	OUTPUT:
	RETVAL



boolean
set_score(obj,score)
	Wise2_DataScore * obj
	int score
	CODE:
	RETVAL = Wise2_replace_score_DataScore(obj,score);
	OUTPUT:
	RETVAL



int
score(obj)
	Wise2_DataScore * obj
	CODE:
	RETVAL = Wise2_access_score_DataScore(obj);
	OUTPUT:
	RETVAL



boolean
set_evalue(obj,evalue)
	Wise2_DataScore * obj
	double evalue
	CODE:
	RETVAL = Wise2_replace_evalue_DataScore(obj,evalue);
	OUTPUT:
	RETVAL



double
evalue(obj)
	Wise2_DataScore * obj
	CODE:
	RETVAL = Wise2_access_evalue_DataScore(obj);
	OUTPUT:
	RETVAL




Wise2_DataScore *
new(class)
	char * class
	PPCODE:
	Wise2_DataScore * out;
	out = Wise2_DataScore_alloc();
	ST(0) = sv_newmortal();
	sv_setref_pv(ST(0),class,(void*)out);
	XSRETURN(1);

void
DESTROY(obj)
	Wise2_DataScore * obj
	CODE:
	Wise2_free_DataScore(obj);



MODULE = Wise2 PACKAGE = Wise2::DataEntry

Wise2_DataEntry *
hard_link_DataEntry(obj)
	Wise2_DataEntry * obj
	CODE:
	RETVAL = Wise2_hard_link_DataEntry(obj);
	OUTPUT:
	RETVAL



Wise2_DataEntry *
alloc()
	CODE:
	RETVAL = Wise2_DataEntry_alloc();
	OUTPUT:
	RETVAL



boolean
set_name(obj,name)
	Wise2_DataEntry * obj
	char * name
	CODE:
	RETVAL = Wise2_replace_name_DataEntry(obj,Wise2_stringalloc(name));
	OUTPUT:
	RETVAL



char *
name(obj)
	Wise2_DataEntry * obj
	INIT:
	char * temp;
	CODE:
	temp = Wise2_stringalloc(Wise2_access_name_DataEntry(obj));
	RETVAL = temp;
	OUTPUT:
	RETVAL



boolean
set_is_reversed(obj,is_reversed)
	Wise2_DataEntry * obj
	boolean is_reversed
	CODE:
	RETVAL = Wise2_replace_is_reversed_DataEntry(obj,is_reversed);
	OUTPUT:
	RETVAL



boolean
is_reversed(obj)
	Wise2_DataEntry * obj
	CODE:
	RETVAL = Wise2_access_is_reversed_DataEntry(obj);
	OUTPUT:
	RETVAL




Wise2_DataEntry *
new(class)
	char * class
	PPCODE:
	Wise2_DataEntry * out;
	out = Wise2_DataEntry_alloc();
	ST(0) = sv_newmortal();
	sv_setref_pv(ST(0),class,(void*)out);
	XSRETURN(1);

void
DESTROY(obj)
	Wise2_DataEntry * obj
	CODE:
	Wise2_free_DataEntry(obj);



MODULE = Wise2 PACKAGE = Wise2

Wise2_Hscore *
std_score_Hscore(cut_off,report_stagger)
	int cut_off
	int report_stagger
	CODE:
	RETVAL = Wise2_std_score_Hscore(cut_off,report_stagger);
	OUTPUT:
	RETVAL



