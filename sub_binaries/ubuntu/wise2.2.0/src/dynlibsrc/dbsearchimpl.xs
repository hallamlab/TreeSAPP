

MODULE = Wise2 PACKAGE = Wise2::DBSearchImpl

char *
string(dbsi)
	Wise2_DBSearchImpl * dbsi
	INIT:
	char * temp;
	CODE:
	temp = Wise2_stringalloc(Wise2_impl_string_DBSearchImpl(dbsi));
	RETVAL = temp;
	OUTPUT:
	RETVAL



Wise2_DBSearchImpl *
hard_link_DBSearchImpl(obj)
	Wise2_DBSearchImpl * obj
	CODE:
	RETVAL = Wise2_hard_link_DBSearchImpl(obj);
	OUTPUT:
	RETVAL



Wise2_DBSearchImpl *
alloc()
	CODE:
	RETVAL = Wise2_DBSearchImpl_alloc();
	OUTPUT:
	RETVAL



boolean
set_type(obj,type)
	Wise2_DBSearchImpl * obj
	int type
	CODE:
	RETVAL = Wise2_replace_type_DBSearchImpl(obj,type);
	OUTPUT:
	RETVAL



int
type(obj)
	Wise2_DBSearchImpl * obj
	CODE:
	RETVAL = Wise2_access_type_DBSearchImpl(obj);
	OUTPUT:
	RETVAL



boolean
set_trace_level(obj,trace_level)
	Wise2_DBSearchImpl * obj
	int trace_level
	CODE:
	RETVAL = Wise2_replace_trace_level_DBSearchImpl(obj,trace_level);
	OUTPUT:
	RETVAL



int
trace_level(obj)
	Wise2_DBSearchImpl * obj
	CODE:
	RETVAL = Wise2_access_trace_level_DBSearchImpl(obj);
	OUTPUT:
	RETVAL



boolean
set_trace_file(obj,trace_file)
	Wise2_DBSearchImpl * obj
	FILE * trace_file
	CODE:
	RETVAL = Wise2_replace_trace_file_DBSearchImpl(obj,trace_file);
	OUTPUT:
	RETVAL



FILE *
trace_file(obj)
	Wise2_DBSearchImpl * obj
	CODE:
	RETVAL = Wise2_access_trace_file_DBSearchImpl(obj);
	OUTPUT:
	RETVAL



boolean
set_suggest_thread_no(obj,suggest_thread_no)
	Wise2_DBSearchImpl * obj
	int suggest_thread_no
	CODE:
	RETVAL = Wise2_replace_suggest_thread_no_DBSearchImpl(obj,suggest_thread_no);
	OUTPUT:
	RETVAL



int
suggest_thread_no(obj)
	Wise2_DBSearchImpl * obj
	CODE:
	RETVAL = Wise2_access_suggest_thread_no_DBSearchImpl(obj);
	OUTPUT:
	RETVAL



boolean
set_search_routine(obj,search_routine)
	Wise2_DBSearchImpl * obj
	int search_routine
	CODE:
	RETVAL = Wise2_replace_search_routine_DBSearchImpl(obj,search_routine);
	OUTPUT:
	RETVAL



int
search_routine(obj)
	Wise2_DBSearchImpl * obj
	CODE:
	RETVAL = Wise2_access_search_routine_DBSearchImpl(obj);
	OUTPUT:
	RETVAL




Wise2_DBSearchImpl *
new(class)
	char * class
	PPCODE:
	Wise2_DBSearchImpl * out;
	out = Wise2_DBSearchImpl_alloc();
	ST(0) = sv_newmortal();
	sv_setref_pv(ST(0),class,(void*)out);
	XSRETURN(1);

void
DESTROY(obj)
	Wise2_DBSearchImpl * obj
	CODE:
	Wise2_free_DBSearchImpl(obj);



MODULE = Wise2 PACKAGE = Wise2

Wise2_DBSearchImpl *
new_pthread_DBSearchImpl()
	CODE:
	RETVAL = Wise2_new_pthread_DBSearchImpl();
	OUTPUT:
	RETVAL



Wise2_DBSearchImpl *
new_serial_DBSearchImpl()
	CODE:
	RETVAL = Wise2_new_serial_DBSearchImpl();
	OUTPUT:
	RETVAL



