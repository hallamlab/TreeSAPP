

MODULE = Wise2 PACKAGE = Wise2::AlnBlock

void
dump_ascii_AlnBlock(alb,ofp)
	Wise2_AlnBlock * alb
	FILE * ofp
	CODE:
	Wise2_dump_ascii_AlnBlock(alb,ofp);



Wise2_AlnBlock *
hard_link_AlnBlock(obj)
	Wise2_AlnBlock * obj
	CODE:
	RETVAL = Wise2_hard_link_AlnBlock(obj);
	OUTPUT:
	RETVAL



Wise2_AlnBlock *
AlnBlock_alloc_std()
	CODE:
	RETVAL = Wise2_AlnBlock_alloc_std();
	OUTPUT:
	RETVAL



boolean
set_start(obj,start)
	Wise2_AlnBlock * obj
	Wise2_AlnColumn * start
	CODE:
	RETVAL = Wise2_replace_start_AlnBlock(obj,Wise2_hard_link_AlnColumn(start));
	OUTPUT:
	RETVAL



Wise2_AlnColumn *
start(obj)
	Wise2_AlnBlock * obj
	INIT:
Wise2_AlnColumn * temp;
	CODE:
	temp = Wise2_hard_link_AlnColumn(Wise2_access_start_AlnBlock(obj));
	RETVAL = temp;
	OUTPUT:
	RETVAL



Wise2_AlnSequence *
seq(obj,i)
	Wise2_AlnBlock * obj
	int i
	INIT:
Wise2_AlnSequence * temp;
	CODE:
	temp = Wise2_hard_link_AlnSequence(Wise2_access_seq_AlnBlock(obj,i));
	RETVAL = temp;
	OUTPUT:
	RETVAL



int
length_seq(obj)
	Wise2_AlnBlock * obj
	CODE:
	RETVAL = Wise2_length_seq_AlnBlock(obj);
	OUTPUT:
	RETVAL



int
flush_seq(obj)
	Wise2_AlnBlock * obj
	CODE:
	RETVAL = Wise2_flush_AlnBlock(obj);
	OUTPUT:
	RETVAL



boolean
add_seq(obj,add)
	Wise2_AlnBlock * obj
	Wise2_AlnSequence * add
	CODE:
	RETVAL = Wise2_add_AlnBlock(obj,Wise2_hard_link_AlnSequence(add));
	OUTPUT:
	RETVAL



boolean
set_length(obj,length)
	Wise2_AlnBlock * obj
	int length
	CODE:
	RETVAL = Wise2_replace_length_AlnBlock(obj,length);
	OUTPUT:
	RETVAL



int
length(obj)
	Wise2_AlnBlock * obj
	CODE:
	RETVAL = Wise2_access_length_AlnBlock(obj);
	OUTPUT:
	RETVAL



boolean
set_score(obj,score)
	Wise2_AlnBlock * obj
	int score
	CODE:
	RETVAL = Wise2_replace_score_AlnBlock(obj,score);
	OUTPUT:
	RETVAL



int
score(obj)
	Wise2_AlnBlock * obj
	CODE:
	RETVAL = Wise2_access_score_AlnBlock(obj);
	OUTPUT:
	RETVAL




Wise2_AlnBlock *
new(class)
	char * class
	PPCODE:
	Wise2_AlnBlock * out;
	out = Wise2_AlnBlock_alloc_std();
	ST(0) = sv_newmortal();
	sv_setref_pv(ST(0),class,(void*)out);
	XSRETURN(1);

void
DESTROY(obj)
	Wise2_AlnBlock * obj
	CODE:
	Wise2_free_AlnBlock(obj);

void
each_seq(obj)
	Wise2_AlnBlock * obj
	PPCODE:
	int i=0;
	int len;
	SV* temp;
	len = Wise2_length_seq_AlnBlock(obj);
	for(i=0;i<len;i++){
	  temp = sv_newmortal();
	  sv_setref_pv(temp, "Wise2::AlnSequence", (void*) (Wise2_hard_link_AlnSequence(Wise2_access_seq_AlnBlock(obj,i))));
	  XPUSHs(temp);
	  }
	XSRETURN(len);



MODULE = Wise2 PACKAGE = Wise2::AlnColumn

boolean
at_end(alc)
	Wise2_AlnColumn * alc
	CODE:
	RETVAL = Wise2_at_end_AlnColumn(alc);
	OUTPUT:
	RETVAL



Wise2_AlnColumn *
hard_link_AlnColumn(obj)
	Wise2_AlnColumn * obj
	CODE:
	RETVAL = Wise2_hard_link_AlnColumn(obj);
	OUTPUT:
	RETVAL



Wise2_AlnColumn *
AlnColumn_alloc_std()
	CODE:
	RETVAL = Wise2_AlnColumn_alloc_std();
	OUTPUT:
	RETVAL



Wise2_AlnUnit *
alu(obj,i)
	Wise2_AlnColumn * obj
	int i
	INIT:
Wise2_AlnUnit * temp;
	CODE:
	temp = Wise2_hard_link_AlnUnit(Wise2_access_alu_AlnColumn(obj,i));
	RETVAL = temp;
	OUTPUT:
	RETVAL



int
length_alu(obj)
	Wise2_AlnColumn * obj
	CODE:
	RETVAL = Wise2_length_alu_AlnColumn(obj);
	OUTPUT:
	RETVAL



int
flush_alu(obj)
	Wise2_AlnColumn * obj
	CODE:
	RETVAL = Wise2_flush_AlnColumn(obj);
	OUTPUT:
	RETVAL



boolean
add_alu(obj,add)
	Wise2_AlnColumn * obj
	Wise2_AlnUnit * add
	CODE:
	RETVAL = Wise2_add_AlnColumn(obj,Wise2_hard_link_AlnUnit(add));
	OUTPUT:
	RETVAL



boolean
set_next(obj,next)
	Wise2_AlnColumn * obj
	Wise2_AlnColumn * next
	CODE:
	RETVAL = Wise2_replace_next_AlnColumn(obj,Wise2_hard_link_AlnColumn(next));
	OUTPUT:
	RETVAL



Wise2_AlnColumn *
next(obj)
	Wise2_AlnColumn * obj
	INIT:
Wise2_AlnColumn * temp;
	CODE:
	temp = Wise2_hard_link_AlnColumn(Wise2_access_next_AlnColumn(obj));
	RETVAL = temp;
	OUTPUT:
	RETVAL




Wise2_AlnColumn *
new(class)
	char * class
	PPCODE:
	Wise2_AlnColumn * out;
	out = Wise2_AlnColumn_alloc_std();
	ST(0) = sv_newmortal();
	sv_setref_pv(ST(0),class,(void*)out);
	XSRETURN(1);

void
DESTROY(obj)
	Wise2_AlnColumn * obj
	CODE:
	Wise2_free_AlnColumn(obj);

void
each_alu(obj)
	Wise2_AlnColumn * obj
	PPCODE:
	int i=0;
	int len;
	SV* temp;
	len = Wise2_length_alu_AlnColumn(obj);
	for(i=0;i<len;i++){
	  temp = sv_newmortal();
	  sv_setref_pv(temp, "Wise2::AlnUnit", (void*) (Wise2_hard_link_AlnUnit(Wise2_access_alu_AlnColumn(obj,i))));
	  XPUSHs(temp);
	  }
	XSRETURN(len);



MODULE = Wise2 PACKAGE = Wise2::AlnUnit

int
bio_start(alu)
	Wise2_AlnUnit * alu
	CODE:
	RETVAL = Wise2_bio_start_AlnUnit(alu);
	OUTPUT:
	RETVAL



int
bio_end(alu)
	Wise2_AlnUnit * alu
	CODE:
	RETVAL = Wise2_bio_end_AlnUnit(alu);
	OUTPUT:
	RETVAL



Wise2_AlnUnit *
hard_link_AlnUnit(obj)
	Wise2_AlnUnit * obj
	CODE:
	RETVAL = Wise2_hard_link_AlnUnit(obj);
	OUTPUT:
	RETVAL



Wise2_AlnUnit *
alloc()
	CODE:
	RETVAL = Wise2_AlnUnit_alloc();
	OUTPUT:
	RETVAL



boolean
set_start(obj,start)
	Wise2_AlnUnit * obj
	int start
	CODE:
	RETVAL = Wise2_replace_start_AlnUnit(obj,start);
	OUTPUT:
	RETVAL



int
start(obj)
	Wise2_AlnUnit * obj
	CODE:
	RETVAL = Wise2_access_start_AlnUnit(obj);
	OUTPUT:
	RETVAL



boolean
set_end(obj,end)
	Wise2_AlnUnit * obj
	int end
	CODE:
	RETVAL = Wise2_replace_end_AlnUnit(obj,end);
	OUTPUT:
	RETVAL



int
end(obj)
	Wise2_AlnUnit * obj
	CODE:
	RETVAL = Wise2_access_end_AlnUnit(obj);
	OUTPUT:
	RETVAL



boolean
set_label(obj,label)
	Wise2_AlnUnit * obj
	int label
	CODE:
	RETVAL = Wise2_replace_label_AlnUnit(obj,label);
	OUTPUT:
	RETVAL



int
label(obj)
	Wise2_AlnUnit * obj
	CODE:
	RETVAL = Wise2_access_label_AlnUnit(obj);
	OUTPUT:
	RETVAL



boolean
set_text_label(obj,text_label)
	Wise2_AlnUnit * obj
	char * text_label
	CODE:
	RETVAL = Wise2_replace_text_label_AlnUnit(obj,Wise2_stringalloc(text_label));
	OUTPUT:
	RETVAL



char *
text_label(obj)
	Wise2_AlnUnit * obj
	INIT:
	char * temp;
	CODE:
	temp = Wise2_stringalloc(Wise2_access_text_label_AlnUnit(obj));
	RETVAL = temp;
	OUTPUT:
	RETVAL



boolean
set_next(obj,next)
	Wise2_AlnUnit * obj
	Wise2_AlnUnit * next
	CODE:
	RETVAL = Wise2_replace_next_AlnUnit(obj,Wise2_hard_link_AlnUnit(next));
	OUTPUT:
	RETVAL



Wise2_AlnUnit *
next(obj)
	Wise2_AlnUnit * obj
	INIT:
Wise2_AlnUnit * temp;
	CODE:
	temp = Wise2_hard_link_AlnUnit(Wise2_access_next_AlnUnit(obj));
	RETVAL = temp;
	OUTPUT:
	RETVAL



boolean
set_in_column(obj,in_column)
	Wise2_AlnUnit * obj
	boolean in_column
	CODE:
	RETVAL = Wise2_replace_in_column_AlnUnit(obj,in_column);
	OUTPUT:
	RETVAL



boolean
in_column(obj)
	Wise2_AlnUnit * obj
	CODE:
	RETVAL = Wise2_access_in_column_AlnUnit(obj);
	OUTPUT:
	RETVAL



boolean
set_seq(obj,seq)
	Wise2_AlnUnit * obj
	Wise2_AlnSequence * seq
	CODE:
	RETVAL = Wise2_replace_seq_AlnUnit(obj,Wise2_hard_link_AlnSequence(seq));
	OUTPUT:
	RETVAL



Wise2_AlnSequence *
seq(obj)
	Wise2_AlnUnit * obj
	INIT:
Wise2_AlnSequence * temp;
	CODE:
	temp = Wise2_hard_link_AlnSequence(Wise2_access_seq_AlnUnit(obj));
	RETVAL = temp;
	OUTPUT:
	RETVAL




Wise2_AlnUnit *
new(class)
	char * class
	PPCODE:
	Wise2_AlnUnit * out;
	out = Wise2_AlnUnit_alloc();
	ST(0) = sv_newmortal();
	sv_setref_pv(ST(0),class,(void*)out);
	XSRETURN(1);

void
DESTROY(obj)
	Wise2_AlnUnit * obj
	CODE:
	Wise2_free_AlnUnit(obj);



MODULE = Wise2 PACKAGE = Wise2::AlnSequence

Wise2_AlnSequence *
hard_link_AlnSequence(obj)
	Wise2_AlnSequence * obj
	CODE:
	RETVAL = Wise2_hard_link_AlnSequence(obj);
	OUTPUT:
	RETVAL



Wise2_AlnSequence *
alloc()
	CODE:
	RETVAL = Wise2_AlnSequence_alloc();
	OUTPUT:
	RETVAL



boolean
set_start(obj,start)
	Wise2_AlnSequence * obj
	Wise2_AlnUnit * start
	CODE:
	RETVAL = Wise2_replace_start_AlnSequence(obj,Wise2_hard_link_AlnUnit(start));
	OUTPUT:
	RETVAL



Wise2_AlnUnit *
start(obj)
	Wise2_AlnSequence * obj
	INIT:
Wise2_AlnUnit * temp;
	CODE:
	temp = Wise2_hard_link_AlnUnit(Wise2_access_start_AlnSequence(obj));
	RETVAL = temp;
	OUTPUT:
	RETVAL



boolean
set_data_type(obj,data_type)
	Wise2_AlnSequence * obj
	int data_type
	CODE:
	RETVAL = Wise2_replace_data_type_AlnSequence(obj,data_type);
	OUTPUT:
	RETVAL



int
data_type(obj)
	Wise2_AlnSequence * obj
	CODE:
	RETVAL = Wise2_access_data_type_AlnSequence(obj);
	OUTPUT:
	RETVAL



boolean
set_data(obj,data)
	Wise2_AlnSequence * obj
	void * data
	CODE:
	RETVAL = Wise2_replace_data_AlnSequence(obj,data);
	OUTPUT:
	RETVAL



void *
data(obj)
	Wise2_AlnSequence * obj
	CODE:
	RETVAL = Wise2_access_data_AlnSequence(obj);
	OUTPUT:
	RETVAL



boolean
set_bio_start(obj,bio_start)
	Wise2_AlnSequence * obj
	int bio_start
	CODE:
	RETVAL = Wise2_replace_bio_start_AlnSequence(obj,bio_start);
	OUTPUT:
	RETVAL



int
bio_start(obj)
	Wise2_AlnSequence * obj
	CODE:
	RETVAL = Wise2_access_bio_start_AlnSequence(obj);
	OUTPUT:
	RETVAL



boolean
set_bio_end(obj,bio_end)
	Wise2_AlnSequence * obj
	int bio_end
	CODE:
	RETVAL = Wise2_replace_bio_end_AlnSequence(obj,bio_end);
	OUTPUT:
	RETVAL



int
bio_end(obj)
	Wise2_AlnSequence * obj
	CODE:
	RETVAL = Wise2_access_bio_end_AlnSequence(obj);
	OUTPUT:
	RETVAL




Wise2_AlnSequence *
new(class)
	char * class
	PPCODE:
	Wise2_AlnSequence * out;
	out = Wise2_AlnSequence_alloc();
	ST(0) = sv_newmortal();
	sv_setref_pv(ST(0),class,(void*)out);
	XSRETURN(1);

void
DESTROY(obj)
	Wise2_AlnSequence * obj
	CODE:
	Wise2_free_AlnSequence(obj);



MODULE = Wise2 PACKAGE = Wise2

