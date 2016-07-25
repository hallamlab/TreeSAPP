

MODULE = Wise2 PACKAGE = Wise2

boolean
write_pretty_str_align(alb,qname,query,tname,target,name,main,ofp)
	Wise2_AlnBlock * alb
	char * qname
	char * query
	char * tname
	char * target
	int name
	int main
	FILE * ofp
	CODE:
	RETVAL = Wise2_write_pretty_str_align(alb,qname,query,tname,target,name,main,ofp);
	OUTPUT:
	RETVAL



boolean
write_pretty_seq_align(alb,q,t,name,main,ofp)
	Wise2_AlnBlock * alb
	Wise2_Sequence * q
	Wise2_Sequence * t
	int name
	int main
	FILE * ofp
	CODE:
	RETVAL = Wise2_write_pretty_seq_align(alb,q,t,name,main,ofp);
	OUTPUT:
	RETVAL



boolean
write_pretty_Protein_align(alb,q,t,name,main,ofp)
	Wise2_AlnBlock * alb
	Wise2_Protein * q
	Wise2_Protein * t
	int name
	int main
	FILE * ofp
	CODE:
	RETVAL = Wise2_write_pretty_Protein_align(alb,q,t,name,main,ofp);
	OUTPUT:
	RETVAL



