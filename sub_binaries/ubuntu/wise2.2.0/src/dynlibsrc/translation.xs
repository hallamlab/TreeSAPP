

MODULE = Wise2 PACKAGE = Wise2::Translation

Wise2_Protein *
get_Protein_from_Translation(ts,ct)
	Wise2_Translation * ts
	Wise2_CodonTable * ct
	INIT:
Wise2_Protein * temp;
	CODE:
	temp = Wise2_hard_link_Protein(Wise2_get_Protein_from_Translation(ts,ct));
	RETVAL = temp;
	OUTPUT:
	RETVAL



Wise2_Translation *
hard_link_Translation(obj)
	Wise2_Translation * obj
	CODE:
	RETVAL = Wise2_hard_link_Translation(obj);
	OUTPUT:
	RETVAL



Wise2_Translation *
alloc()
	CODE:
	RETVAL = Wise2_Translation_alloc();
	OUTPUT:
	RETVAL



boolean
set_start(obj,start)
	Wise2_Translation * obj
	int start
	CODE:
	RETVAL = Wise2_replace_start_Translation(obj,start);
	OUTPUT:
	RETVAL



int
start(obj)
	Wise2_Translation * obj
	CODE:
	RETVAL = Wise2_access_start_Translation(obj);
	OUTPUT:
	RETVAL



boolean
set_end(obj,end)
	Wise2_Translation * obj
	int end
	CODE:
	RETVAL = Wise2_replace_end_Translation(obj,end);
	OUTPUT:
	RETVAL



int
end(obj)
	Wise2_Translation * obj
	CODE:
	RETVAL = Wise2_access_end_Translation(obj);
	OUTPUT:
	RETVAL



boolean
set_parent(obj,parent)
	Wise2_Translation * obj
	Wise2_Transcript * parent
	CODE:
	RETVAL = Wise2_replace_parent_Translation(obj,Wise2_hard_link_Transcript(parent));
	OUTPUT:
	RETVAL



Wise2_Transcript *
parent(obj)
	Wise2_Translation * obj
	INIT:
Wise2_Transcript * temp;
	CODE:
	temp = Wise2_hard_link_Transcript(Wise2_access_parent_Translation(obj));
	RETVAL = temp;
	OUTPUT:
	RETVAL



boolean
set_protein(obj,protein)
	Wise2_Translation * obj
	Wise2_Protein * protein
	CODE:
	RETVAL = Wise2_replace_protein_Translation(obj,Wise2_hard_link_Protein(protein));
	OUTPUT:
	RETVAL



Wise2_Protein *
protein(obj)
	Wise2_Translation * obj
	INIT:
Wise2_Protein * temp;
	CODE:
	temp = Wise2_hard_link_Protein(Wise2_access_protein_Translation(obj));
	RETVAL = temp;
	OUTPUT:
	RETVAL




Wise2_Translation *
new(class)
	char * class
	PPCODE:
	Wise2_Translation * out;
	out = Wise2_Translation_alloc();
	ST(0) = sv_newmortal();
	sv_setref_pv(ST(0),class,(void*)out);
	XSRETURN(1);

void
DESTROY(obj)
	Wise2_Translation * obj
	CODE:
	Wise2_free_Translation(obj);



MODULE = Wise2 PACKAGE = Wise2

