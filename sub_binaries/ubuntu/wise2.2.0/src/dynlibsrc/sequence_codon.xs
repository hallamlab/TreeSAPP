

MODULE = Wise2 PACKAGE = Wise2

Wise2_Sequence *
reverse_complement_Sequence(seq)
	Wise2_Sequence * seq
	CODE:
	RETVAL = Wise2_reverse_complement_Sequence(seq);
	OUTPUT:
	RETVAL



Wise2_Sequence *
magic_trunc_Sequence(seq,start,end)
	Wise2_Sequence * seq
	int start
	int end
	CODE:
	RETVAL = Wise2_magic_trunc_Sequence(seq,start,end);
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



