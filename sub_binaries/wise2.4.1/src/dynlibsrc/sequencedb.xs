

MODULE = Wise2 PACKAGE = Wise2::SequenceDB

boolean
close_SequenceDB(last,sdb)
	Wise2_Sequence * last
	Wise2_SequenceDB * sdb
	CODE:
	RETVAL = Wise2_close_SequenceDB(last,sdb);
	OUTPUT:
	RETVAL



Wise2_SequenceDB *
hard_link_SequenceDB(obj)
	Wise2_SequenceDB * obj
	CODE:
	RETVAL = Wise2_hard_link_SequenceDB(obj);
	OUTPUT:
	RETVAL



Wise2_SequenceDB *
SequenceDB_alloc_std()
	CODE:
	RETVAL = Wise2_SequenceDB_alloc_std();
	OUTPUT:
	RETVAL



boolean
set_name(obj,name)
	Wise2_SequenceDB * obj
	char * name
	CODE:
	RETVAL = Wise2_replace_name_SequenceDB(obj,Wise2_stringalloc(name));
	OUTPUT:
	RETVAL



char *
name(obj)
	Wise2_SequenceDB * obj
	INIT:
	char * temp;
	CODE:
	temp = Wise2_stringalloc(Wise2_access_name_SequenceDB(obj));
	RETVAL = temp;
	OUTPUT:
	RETVAL



Wise2_FileSource *
fs(obj,i)
	Wise2_SequenceDB * obj
	int i
	INIT:
Wise2_FileSource * temp;
	CODE:
	temp = Wise2_hard_link_FileSource(Wise2_access_fs_SequenceDB(obj,i));
	RETVAL = temp;
	OUTPUT:
	RETVAL



int
length_fs(obj)
	Wise2_SequenceDB * obj
	CODE:
	RETVAL = Wise2_length_fs_SequenceDB(obj);
	OUTPUT:
	RETVAL



int
flush_fs(obj)
	Wise2_SequenceDB * obj
	CODE:
	RETVAL = Wise2_flush_SequenceDB(obj);
	OUTPUT:
	RETVAL



boolean
add_fs(obj,add)
	Wise2_SequenceDB * obj
	Wise2_FileSource * add
	CODE:
	RETVAL = Wise2_add_SequenceDB(obj,Wise2_hard_link_FileSource(add));
	OUTPUT:
	RETVAL



boolean
set_current_source(obj,current_source)
	Wise2_SequenceDB * obj
	int current_source
	CODE:
	RETVAL = Wise2_replace_current_source_SequenceDB(obj,current_source);
	OUTPUT:
	RETVAL



int
current_source(obj)
	Wise2_SequenceDB * obj
	CODE:
	RETVAL = Wise2_access_current_source_SequenceDB(obj);
	OUTPUT:
	RETVAL



boolean
set_current_file(obj,current_file)
	Wise2_SequenceDB * obj
	FILE * current_file
	CODE:
	RETVAL = Wise2_replace_current_file_SequenceDB(obj,current_file);
	OUTPUT:
	RETVAL



FILE *
current_file(obj)
	Wise2_SequenceDB * obj
	CODE:
	RETVAL = Wise2_access_current_file_SequenceDB(obj);
	OUTPUT:
	RETVAL



boolean
set_sequence_no(obj,sequence_no)
	Wise2_SequenceDB * obj
	int sequence_no
	CODE:
	RETVAL = Wise2_replace_sequence_no_SequenceDB(obj,sequence_no);
	OUTPUT:
	RETVAL



int
sequence_no(obj)
	Wise2_SequenceDB * obj
	CODE:
	RETVAL = Wise2_access_sequence_no_SequenceDB(obj);
	OUTPUT:
	RETVAL



boolean
set_byte_position(obj,byte_position)
	Wise2_SequenceDB * obj
	int byte_position
	CODE:
	RETVAL = Wise2_replace_byte_position_SequenceDB(obj,byte_position);
	OUTPUT:
	RETVAL



int
byte_position(obj)
	Wise2_SequenceDB * obj
	CODE:
	RETVAL = Wise2_access_byte_position_SequenceDB(obj);
	OUTPUT:
	RETVAL




Wise2_SequenceDB *
new(class)
	char * class
	PPCODE:
	Wise2_SequenceDB * out;
	out = Wise2_SequenceDB_alloc_std();
	ST(0) = sv_newmortal();
	sv_setref_pv(ST(0),class,(void*)out);
	XSRETURN(1);

void
DESTROY(obj)
	Wise2_SequenceDB * obj
	CODE:
	Wise2_free_SequenceDB(obj);

void
each_fs(obj)
	Wise2_SequenceDB * obj
	PPCODE:
	int i=0;
	int len;
	SV* temp;
	len = Wise2_length_fs_SequenceDB(obj);
	for(i=0;i<len;i++){
	  temp = sv_newmortal();
	  sv_setref_pv(temp, "Wise2::FileSource", (void*) (Wise2_hard_link_FileSource(Wise2_access_fs_SequenceDB(obj,i))));
	  XPUSHs(temp);
	  }
	XSRETURN(len);



MODULE = Wise2 PACKAGE = Wise2::FileSource

Wise2_FileSource *
hard_link_FileSource(obj)
	Wise2_FileSource * obj
	CODE:
	RETVAL = Wise2_hard_link_FileSource(obj);
	OUTPUT:
	RETVAL



Wise2_FileSource *
alloc()
	CODE:
	RETVAL = Wise2_FileSource_alloc();
	OUTPUT:
	RETVAL



boolean
set_filename(obj,filename)
	Wise2_FileSource * obj
	char * filename
	CODE:
	RETVAL = Wise2_replace_filename_FileSource(obj,Wise2_stringalloc(filename));
	OUTPUT:
	RETVAL



char *
filename(obj)
	Wise2_FileSource * obj
	INIT:
	char * temp;
	CODE:
	temp = Wise2_stringalloc(Wise2_access_filename_FileSource(obj));
	RETVAL = temp;
	OUTPUT:
	RETVAL



boolean
set_input(obj,input)
	Wise2_FileSource * obj
	FILE * input
	CODE:
	RETVAL = Wise2_replace_input_FileSource(obj,input);
	OUTPUT:
	RETVAL



FILE *
input(obj)
	Wise2_FileSource * obj
	CODE:
	RETVAL = Wise2_access_input_FileSource(obj);
	OUTPUT:
	RETVAL



boolean
set_format(obj,format)
	Wise2_FileSource * obj
	int format
	CODE:
	RETVAL = Wise2_replace_format_FileSource(obj,format);
	OUTPUT:
	RETVAL



int
format(obj)
	Wise2_FileSource * obj
	CODE:
	RETVAL = Wise2_access_format_FileSource(obj);
	OUTPUT:
	RETVAL



boolean
set_type(obj,type)
	Wise2_FileSource * obj
	int type
	CODE:
	RETVAL = Wise2_replace_type_FileSource(obj,type);
	OUTPUT:
	RETVAL



int
type(obj)
	Wise2_FileSource * obj
	CODE:
	RETVAL = Wise2_access_type_FileSource(obj);
	OUTPUT:
	RETVAL




Wise2_FileSource *
new(class)
	char * class
	PPCODE:
	Wise2_FileSource * out;
	out = Wise2_FileSource_alloc();
	ST(0) = sv_newmortal();
	sv_setref_pv(ST(0),class,(void*)out);
	XSRETURN(1);

void
DESTROY(obj)
	Wise2_FileSource * obj
	CODE:
	Wise2_free_FileSource(obj);



MODULE = Wise2 PACKAGE = Wise2

Wise2_SequenceDB *
single_fasta_SequenceDB(filename)
	char * filename
	CODE:
	RETVAL = Wise2_single_fasta_SequenceDB(filename);
	OUTPUT:
	RETVAL



