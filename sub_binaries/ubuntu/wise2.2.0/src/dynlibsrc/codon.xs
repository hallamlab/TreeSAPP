

MODULE = Wise2 PACKAGE = Wise2::CodonTable

Wise2_CodonTable *
read_CodonTable_file(file)
	char * file
	CODE:
	RETVAL = Wise2_read_CodonTable_file(file);
	OUTPUT:
	RETVAL



Wise2_CodonTable *
read_CodonTable(ifp)
	FILE * ifp
	CODE:
	RETVAL = Wise2_read_CodonTable(ifp);
	OUTPUT:
	RETVAL



aa
aminoacid_from_seq(ct,seq)
	Wise2_CodonTable * ct
	char * seq
	CODE:
	RETVAL = Wise2_aminoacid_from_seq(ct,seq);
	OUTPUT:
	RETVAL



aa
aminoacid_from_codon(ct,c)
	Wise2_CodonTable * ct
	codon c
	CODE:
	RETVAL = Wise2_aminoacid_from_codon(ct,c);
	OUTPUT:
	RETVAL



boolean
is_stop_codon(c,ct)
	codon c
	Wise2_CodonTable * ct
	CODE:
	RETVAL = Wise2_is_stop_codon(c,ct);
	OUTPUT:
	RETVAL



boolean
is_valid_aminoacid(ct,c)
	Wise2_CodonTable * ct
	char c
	CODE:
	RETVAL = Wise2_is_valid_aminoacid(ct,c);
	OUTPUT:
	RETVAL



Wise2_CodonTable *
hard_link_CodonTable(obj)
	Wise2_CodonTable * obj
	CODE:
	RETVAL = Wise2_hard_link_CodonTable(obj);
	OUTPUT:
	RETVAL



Wise2_CodonTable *
alloc()
	CODE:
	RETVAL = Wise2_CodonTable_alloc();
	OUTPUT:
	RETVAL



boolean
set_name(obj,name)
	Wise2_CodonTable * obj
	char * name
	CODE:
	RETVAL = Wise2_replace_name_CodonTable(obj,Wise2_stringalloc(name));
	OUTPUT:
	RETVAL



char *
name(obj)
	Wise2_CodonTable * obj
	INIT:
	char * temp;
	CODE:
	temp = Wise2_stringalloc(Wise2_access_name_CodonTable(obj));
	RETVAL = temp;
	OUTPUT:
	RETVAL




Wise2_CodonTable *
new(class)
	char * class
	PPCODE:
	Wise2_CodonTable * out;
	out = Wise2_CodonTable_alloc();
	ST(0) = sv_newmortal();
	sv_setref_pv(ST(0),class,(void*)out);
	XSRETURN(1);

void
DESTROY(obj)
	Wise2_CodonTable * obj
	CODE:
	Wise2_free_CodonTable(obj);



MODULE = Wise2 PACKAGE = Wise2

boolean
is_non_ambiguous_codon_seq(seq)
	char * seq
	CODE:
	RETVAL = Wise2_is_non_ambiguous_codon_seq(seq);
	OUTPUT:
	RETVAL



codon
codon_from_base4_codon(c)
	int c
	CODE:
	RETVAL = Wise2_codon_from_base4_codon(c);
	OUTPUT:
	RETVAL



int
base4_codon_from_codon(c)
	codon c
	CODE:
	RETVAL = Wise2_base4_codon_from_codon(c);
	OUTPUT:
	RETVAL



boolean
has_random_bases(c)
	codon c
	CODE:
	RETVAL = Wise2_has_random_bases(c);
	OUTPUT:
	RETVAL



codon
permute_possible_random_bases(c,one,two,three)
	codon c
	base one
	base two
	base three
	CODE:
	RETVAL = Wise2_permute_possible_random_bases(c,one,two,three);
	OUTPUT:
	RETVAL



base
base_from_codon(c,pos)
	codon c
	int pos
	CODE:
	RETVAL = Wise2_base_from_codon(c,pos);
	OUTPUT:
	RETVAL



codon
codon_from_seq(seq)
	char * seq
	CODE:
	RETVAL = Wise2_codon_from_seq(seq);
	OUTPUT:
	RETVAL



int
base4_codon_from_seq(seq)
	char * seq
	CODE:
	RETVAL = Wise2_base4_codon_from_seq(seq);
	OUTPUT:
	RETVAL



char
char_from_base(b)
	base b
	CODE:
	RETVAL = Wise2_char_from_base(b);
	OUTPUT:
	RETVAL



base
base_from_char(c)
	char c
	CODE:
	RETVAL = Wise2_base_from_char(c);
	OUTPUT:
	RETVAL



char
char_complement_base(c)
	char c
	CODE:
	RETVAL = Wise2_char_complement_base(c);
	OUTPUT:
	RETVAL



base
complement_base(b)
	base b
	CODE:
	RETVAL = Wise2_complement_base(b);
	OUTPUT:
	RETVAL



