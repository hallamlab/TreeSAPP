

MODULE = Wise2 PACKAGE = Wise2::Sequence

void
uppercase(seq)
	Wise2_Sequence * seq
	CODE:
	Wise2_uppercase_Sequence(seq);



boolean
force_to_dna(seq,fraction)
	Wise2_Sequence * seq
	double fraction
	CODE:
	RETVAL = Wise2_force_to_dna_Sequence(seq,fraction,NULL);
	OUTPUT:
	RETVAL



boolean
is_reversed(seq)
	Wise2_Sequence * seq
	CODE:
	RETVAL = Wise2_is_reversed_Sequence(seq);
	OUTPUT:
	RETVAL



Wise2_Sequence *
translate(dna,ct)
	Wise2_Sequence * dna
	Wise2_CodonTable * ct
	CODE:
	RETVAL = Wise2_translate_Sequence(dna,ct);
	OUTPUT:
	RETVAL



Wise2_Sequence *
revcomp(seq)
	Wise2_Sequence * seq
	CODE:
	RETVAL = Wise2_reverse_complement_Sequence(seq);
	OUTPUT:
	RETVAL



Wise2_Sequence *
magic_trunc(seq,start,end)
	Wise2_Sequence * seq
	int start
	int end
	CODE:
	RETVAL = Wise2_magic_trunc_Sequence(seq,start,end);
	OUTPUT:
	RETVAL



Wise2_Sequence *
trunc(seq,start,end)
	Wise2_Sequence * seq
	int start
	int end
	CODE:
	RETVAL = Wise2_trunc_Sequence(seq,start,end);
	OUTPUT:
	RETVAL



Wise2_Sequence *
read_fasta_file_Sequence(filename)
	char * filename
	CODE:
	RETVAL = Wise2_read_fasta_file_Sequence(filename);
	OUTPUT:
	RETVAL



Wise2_Sequence *
read_Sequence_EMBL_seq(buffer,maxlen,ifp)
	char * buffer
	int maxlen
	FILE * ifp
	CODE:
	RETVAL = Wise2_read_Sequence_EMBL_seq(buffer,maxlen,ifp);
	OUTPUT:
	RETVAL



Wise2_Sequence *
read_fasta_Sequence(ifp)
	FILE * ifp
	CODE:
	RETVAL = Wise2_read_fasta_Sequence(ifp);
	OUTPUT:
	RETVAL



void
show_debug(seq,start,end,ofp)
	Wise2_Sequence * seq
	int start
	int end
	FILE * ofp
	CODE:
	Wise2_show_Sequence_residue_list(seq,start,end,ofp);



void
write_fasta(seq,ofp)
	Wise2_Sequence * seq
	FILE * ofp
	CODE:
	Wise2_write_fasta_Sequence(seq,ofp);



void
validate(seq)
	Wise2_Sequence * seq
	CODE:
	Wise2_make_len_type_Sequence(seq);



Wise2_Sequence *
hard_link_Sequence(obj)
	Wise2_Sequence * obj
	CODE:
	RETVAL = Wise2_hard_link_Sequence(obj);
	OUTPUT:
	RETVAL



Wise2_Sequence *
alloc()
	CODE:
	RETVAL = Wise2_Sequence_alloc();
	OUTPUT:
	RETVAL



boolean
set_name(obj,name)
	Wise2_Sequence * obj
	char * name
	CODE:
	RETVAL = Wise2_replace_name_Sequence(obj,Wise2_stringalloc(name));
	OUTPUT:
	RETVAL



char *
name(obj)
	Wise2_Sequence * obj
	INIT:
	char * temp;
	CODE:
	temp = Wise2_stringalloc(Wise2_access_name_Sequence(obj));
	RETVAL = temp;
	OUTPUT:
	RETVAL



boolean
set_seq(obj,seq)
	Wise2_Sequence * obj
	char * seq
	CODE:
	RETVAL = Wise2_replace_seq_Sequence(obj,Wise2_stringalloc(seq));
	OUTPUT:
	RETVAL



char *
seq(obj)
	Wise2_Sequence * obj
	INIT:
	char * temp;
	CODE:
	temp = Wise2_stringalloc(Wise2_access_seq_Sequence(obj));
	RETVAL = temp;
	OUTPUT:
	RETVAL



boolean
set_len(obj,len)
	Wise2_Sequence * obj
	int len
	CODE:
	RETVAL = Wise2_replace_len_Sequence(obj,len);
	OUTPUT:
	RETVAL



int
len(obj)
	Wise2_Sequence * obj
	CODE:
	RETVAL = Wise2_access_len_Sequence(obj);
	OUTPUT:
	RETVAL



boolean
set_maxlen(obj,maxlen)
	Wise2_Sequence * obj
	int maxlen
	CODE:
	RETVAL = Wise2_replace_maxlen_Sequence(obj,maxlen);
	OUTPUT:
	RETVAL



int
maxlen(obj)
	Wise2_Sequence * obj
	CODE:
	RETVAL = Wise2_access_maxlen_Sequence(obj);
	OUTPUT:
	RETVAL



boolean
set_offset(obj,offset)
	Wise2_Sequence * obj
	int offset
	CODE:
	RETVAL = Wise2_replace_offset_Sequence(obj,offset);
	OUTPUT:
	RETVAL



int
offset(obj)
	Wise2_Sequence * obj
	CODE:
	RETVAL = Wise2_access_offset_Sequence(obj);
	OUTPUT:
	RETVAL



boolean
set_end(obj,end)
	Wise2_Sequence * obj
	int end
	CODE:
	RETVAL = Wise2_replace_end_Sequence(obj,end);
	OUTPUT:
	RETVAL



int
end(obj)
	Wise2_Sequence * obj
	CODE:
	RETVAL = Wise2_access_end_Sequence(obj);
	OUTPUT:
	RETVAL



boolean
set_type(obj,type)
	Wise2_Sequence * obj
	char type
	CODE:
	RETVAL = Wise2_replace_type_Sequence(obj,type);
	OUTPUT:
	RETVAL



char
type(obj)
	Wise2_Sequence * obj
	CODE:
	RETVAL = Wise2_access_type_Sequence(obj);
	OUTPUT:
	RETVAL



boolean
set_tax_id(obj,tax_id)
	Wise2_Sequence * obj
	int tax_id
	CODE:
	RETVAL = Wise2_replace_tax_id_Sequence(obj,tax_id);
	OUTPUT:
	RETVAL



int
tax_id(obj)
	Wise2_Sequence * obj
	CODE:
	RETVAL = Wise2_access_tax_id_Sequence(obj);
	OUTPUT:
	RETVAL



boolean
set_weight(obj,weight)
	Wise2_Sequence * obj
	double weight
	CODE:
	RETVAL = Wise2_replace_weight_Sequence(obj,weight);
	OUTPUT:
	RETVAL



double
weight(obj)
	Wise2_Sequence * obj
	CODE:
	RETVAL = Wise2_access_weight_Sequence(obj);
	OUTPUT:
	RETVAL



boolean
set_desc(obj,desc)
	Wise2_Sequence * obj
	char * desc
	CODE:
	RETVAL = Wise2_replace_desc_Sequence(obj,Wise2_stringalloc(desc));
	OUTPUT:
	RETVAL



char *
desc(obj)
	Wise2_Sequence * obj
	INIT:
	char * temp;
	CODE:
	temp = Wise2_stringalloc(Wise2_access_desc_Sequence(obj));
	RETVAL = temp;
	OUTPUT:
	RETVAL




Wise2_Sequence *
new(class)
	char * class
	PPCODE:
	Wise2_Sequence * out;
	out = Wise2_Sequence_alloc();
	ST(0) = sv_newmortal();
	sv_setref_pv(ST(0),class,(void*)out);
	XSRETURN(1);

void
DESTROY(obj)
	Wise2_Sequence * obj
	CODE:
	Wise2_free_Sequence(obj);



MODULE = Wise2 PACKAGE = Wise2::SequenceSet

Wise2_SequenceSet *
hard_link_SequenceSet(obj)
	Wise2_SequenceSet * obj
	CODE:
	RETVAL = Wise2_hard_link_SequenceSet(obj);
	OUTPUT:
	RETVAL



Wise2_SequenceSet *
SequenceSet_alloc_std()
	CODE:
	RETVAL = Wise2_SequenceSet_alloc_std();
	OUTPUT:
	RETVAL



Wise2_Sequence *
set(obj,i)
	Wise2_SequenceSet * obj
	int i
	INIT:
Wise2_Sequence * temp;
	CODE:
	temp = Wise2_hard_link_Sequence(Wise2_access_set_SequenceSet(obj,i));
	RETVAL = temp;
	OUTPUT:
	RETVAL



int
length_set(obj)
	Wise2_SequenceSet * obj
	CODE:
	RETVAL = Wise2_length_set_SequenceSet(obj);
	OUTPUT:
	RETVAL



int
flush_set(obj)
	Wise2_SequenceSet * obj
	CODE:
	RETVAL = Wise2_flush_SequenceSet(obj);
	OUTPUT:
	RETVAL



boolean
add_set(obj,add)
	Wise2_SequenceSet * obj
	Wise2_Sequence * add
	CODE:
	RETVAL = Wise2_add_SequenceSet(obj,Wise2_hard_link_Sequence(add));
	OUTPUT:
	RETVAL




Wise2_SequenceSet *
new(class)
	char * class
	PPCODE:
	Wise2_SequenceSet * out;
	out = Wise2_SequenceSet_alloc_std();
	ST(0) = sv_newmortal();
	sv_setref_pv(ST(0),class,(void*)out);
	XSRETURN(1);

void
DESTROY(obj)
	Wise2_SequenceSet * obj
	CODE:
	Wise2_free_SequenceSet(obj);

void
each_set(obj)
	Wise2_SequenceSet * obj
	PPCODE:
	int i=0;
	int len;
	SV* temp;
	len = Wise2_length_set_SequenceSet(obj);
	for(i=0;i<len;i++){
	  temp = sv_newmortal();
	  sv_setref_pv(temp, "Wise2::Sequence", (void*) (Wise2_hard_link_Sequence(Wise2_access_set_SequenceSet(obj,i))));
	  XPUSHs(temp);
	  }
	XSRETURN(len);



MODULE = Wise2 PACKAGE = Wise2

char *
Sequence_type_to_string(type)
	int type
	CODE:
	RETVAL = Wise2_Sequence_type_to_string(type);
	OUTPUT:
	RETVAL



Wise2_Sequence *
new_Sequence_from_strings(name,seq)
	char * name
	char * seq
	CODE:
	RETVAL = Wise2_new_Sequence_from_strings(name,seq);
	OUTPUT:
	RETVAL



