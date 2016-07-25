

MODULE = Wise2 PACKAGE = Wise2

boolean
protein2genomic_ascii_display(alb,p,gen,ct,name,main,ofp)
	Wise2_AlnBlock * alb
	Wise2_Protein * p
	Wise2_Genomic * gen
	Wise2_CodonTable * ct
	int name
	int main
	FILE * ofp
	CODE:
	RETVAL = Wise2_protein2genomic_ascii_display(alb,p,gen,ct,name,main,ofp);
	OUTPUT:
	RETVAL



boolean
protcdna_ascii_display(alb,protsequence,protname,protoff,cdna,ct,name,main,mult,ofp)
	Wise2_AlnBlock * alb
	char * protsequence
	char * protname
	int protoff
	Wise2_cDNA * cdna
	Wise2_CodonTable * ct
	int name
	int main
	boolean mult
	FILE * ofp
	CODE:
	RETVAL = Wise2_protcdna_ascii_display(alb,protsequence,protname,protoff,cdna,ct,name,main,mult,ofp);
	OUTPUT:
	RETVAL



