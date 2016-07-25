

MODULE = Wise2 PACKAGE = Wise2::DnaStartEnd

Wise2_DnaStartEnd *
hard_link_DnaStartEnd(obj)
	Wise2_DnaStartEnd * obj
	CODE:
	RETVAL = Wise2_hard_link_DnaStartEnd(obj);
	OUTPUT:
	RETVAL



Wise2_DnaStartEnd *
alloc()
	CODE:
	RETVAL = Wise2_DnaStartEnd_alloc();
	OUTPUT:
	RETVAL




Wise2_DnaStartEnd *
new(class)
	char * class
	PPCODE:
	Wise2_DnaStartEnd * out;
	out = Wise2_DnaStartEnd_alloc();
	ST(0) = sv_newmortal();
	sv_setref_pv(ST(0),class,(void*)out);
	XSRETURN(1);

void
DESTROY(obj)
	Wise2_DnaStartEnd * obj
	CODE:
	Wise2_free_DnaStartEnd(obj);



MODULE = Wise2 PACKAGE = Wise2

Wise2_AlnBlock *
make_align_dnaalign(one,two,mat,se,qgap,qext,tgap,text,dpri)
	Wise2_Sequence * one
	Wise2_Sequence * two
	Wise2_DnaMatrix * mat
	Wise2_DnaStartEnd * se
	int qgap
	int qext
	int tgap
	int text
	Wise2_DPRunImpl * dpri
	CODE:
	RETVAL = Wise2_make_align_dnaalign(one,two,mat,se,qgap,qext,tgap,text,dpri);
	OUTPUT:
	RETVAL



Wise2_DnaStartEnd *
DnaStartEnd_from_policy(policy)
	char * policy
	CODE:
	RETVAL = Wise2_DnaStartEnd_from_policy(policy);
	OUTPUT:
	RETVAL



