

MODULE = Wise2 PACKAGE = Wise2

void
change_max_BaseMatrix_kbytes(new_kilo_number)
	int new_kilo_number
	CODE:
	Wise2_change_max_BaseMatrix_kbytes(new_kilo_number);



int
get_max_BaseMatrix_kbytes()
	CODE:
	RETVAL = Wise2_get_max_BaseMatrix_kbytes();
	OUTPUT:
	RETVAL



