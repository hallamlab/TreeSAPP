

MODULE = Wise2 PACKAGE = Wise2::ComplexSequence

Wise2_ComplexSequence *
hard_link_ComplexSequence(obj)
	Wise2_ComplexSequence * obj
	CODE:
	RETVAL = Wise2_hard_link_ComplexSequence(obj);
	OUTPUT:
	RETVAL



Wise2_ComplexSequence *
alloc()
	CODE:
	RETVAL = Wise2_ComplexSequence_alloc();
	OUTPUT:
	RETVAL



boolean
set_type(obj,type)
	Wise2_ComplexSequence * obj
	int type
	CODE:
	RETVAL = Wise2_replace_type_ComplexSequence(obj,type);
	OUTPUT:
	RETVAL



int
type(obj)
	Wise2_ComplexSequence * obj
	CODE:
	RETVAL = Wise2_access_type_ComplexSequence(obj);
	OUTPUT:
	RETVAL



boolean
set_seq(obj,seq)
	Wise2_ComplexSequence * obj
	Wise2_Sequence * seq
	CODE:
	RETVAL = Wise2_replace_seq_ComplexSequence(obj,Wise2_hard_link_Sequence(seq));
	OUTPUT:
	RETVAL



Wise2_Sequence *
seq(obj)
	Wise2_ComplexSequence * obj
	INIT:
Wise2_Sequence * temp;
	CODE:
	temp = Wise2_hard_link_Sequence(Wise2_access_seq_ComplexSequence(obj));
	RETVAL = temp;
	OUTPUT:
	RETVAL




Wise2_ComplexSequence *
new(class)
	char * class
	PPCODE:
	Wise2_ComplexSequence * out;
	out = Wise2_ComplexSequence_alloc();
	ST(0) = sv_newmortal();
	sv_setref_pv(ST(0),class,(void*)out);
	XSRETURN(1);

void
DESTROY(obj)
	Wise2_ComplexSequence * obj
	CODE:
	Wise2_free_ComplexSequence(obj);



MODULE = Wise2 PACKAGE = Wise2::ComplexSequenceEvalSet

Wise2_ComplexSequenceEvalSet *
hard_link_ComplexSequenceEvalSet(obj)
	Wise2_ComplexSequenceEvalSet * obj
	CODE:
	RETVAL = Wise2_hard_link_ComplexSequenceEvalSet(obj);
	OUTPUT:
	RETVAL



Wise2_ComplexSequenceEvalSet *
ComplexSequenceEvalSet_alloc_std()
	CODE:
	RETVAL = Wise2_ComplexSequenceEvalSet_alloc_std();
	OUTPUT:
	RETVAL



boolean
set_type(obj,type)
	Wise2_ComplexSequenceEvalSet * obj
	int type
	CODE:
	RETVAL = Wise2_replace_type_ComplexSequenceEvalSet(obj,type);
	OUTPUT:
	RETVAL



int
type(obj)
	Wise2_ComplexSequenceEvalSet * obj
	CODE:
	RETVAL = Wise2_access_type_ComplexSequenceEvalSet(obj);
	OUTPUT:
	RETVAL



boolean
set_has_been_prepared(obj,has_been_prepared)
	Wise2_ComplexSequenceEvalSet * obj
	boolean has_been_prepared
	CODE:
	RETVAL = Wise2_replace_has_been_prepared_ComplexSequenceEvalSet(obj,has_been_prepared);
	OUTPUT:
	RETVAL



boolean
has_been_prepared(obj)
	Wise2_ComplexSequenceEvalSet * obj
	CODE:
	RETVAL = Wise2_access_has_been_prepared_ComplexSequenceEvalSet(obj);
	OUTPUT:
	RETVAL



boolean
set_left_window(obj,left_window)
	Wise2_ComplexSequenceEvalSet * obj
	int left_window
	CODE:
	RETVAL = Wise2_replace_left_window_ComplexSequenceEvalSet(obj,left_window);
	OUTPUT:
	RETVAL



int
left_window(obj)
	Wise2_ComplexSequenceEvalSet * obj
	CODE:
	RETVAL = Wise2_access_left_window_ComplexSequenceEvalSet(obj);
	OUTPUT:
	RETVAL



boolean
set_right_window(obj,right_window)
	Wise2_ComplexSequenceEvalSet * obj
	int right_window
	CODE:
	RETVAL = Wise2_replace_right_window_ComplexSequenceEvalSet(obj,right_window);
	OUTPUT:
	RETVAL



int
right_window(obj)
	Wise2_ComplexSequenceEvalSet * obj
	CODE:
	RETVAL = Wise2_access_right_window_ComplexSequenceEvalSet(obj);
	OUTPUT:
	RETVAL



boolean
set_left_lookback(obj,left_lookback)
	Wise2_ComplexSequenceEvalSet * obj
	int left_lookback
	CODE:
	RETVAL = Wise2_replace_left_lookback_ComplexSequenceEvalSet(obj,left_lookback);
	OUTPUT:
	RETVAL



int
left_lookback(obj)
	Wise2_ComplexSequenceEvalSet * obj
	CODE:
	RETVAL = Wise2_access_left_lookback_ComplexSequenceEvalSet(obj);
	OUTPUT:
	RETVAL




Wise2_ComplexSequenceEvalSet *
new(class)
	char * class
	PPCODE:
	Wise2_ComplexSequenceEvalSet * out;
	out = Wise2_ComplexSequenceEvalSet_alloc_std();
	ST(0) = sv_newmortal();
	sv_setref_pv(ST(0),class,(void*)out);
	XSRETURN(1);

void
DESTROY(obj)
	Wise2_ComplexSequenceEvalSet * obj
	CODE:
	Wise2_free_ComplexSequenceEvalSet(obj);



MODULE = Wise2 PACKAGE = Wise2

