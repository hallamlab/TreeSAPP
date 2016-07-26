

MODULE = Wise2 PACKAGE = Wise2

Wise2_Hscore *
Hscore_from_TSM_estwise(tdb,cdb,cp,cm,rmd,use_syn,alg,bits_cutoff,allN,flat_insert,report_level,die_on_error,dbsi)
	Wise2_ThreeStateDB * tdb
	Wise2_cDNADB * cdb
	Wise2_cDNAParser * cp
	Wise2_CodonMapper * cm
	Wise2_RandomModelDNA * rmd
	boolean use_syn
	int alg
	double bits_cutoff
	Probability allN
	boolean flat_insert
	int report_level
	boolean die_on_error
	Wise2_DBSearchImpl * dbsi
	CODE:
	RETVAL = Wise2_Hscore_from_TSM_estwise(tdb,cdb,cp,cm,rmd,use_syn,alg,bits_cutoff,allN,flat_insert,report_level,die_on_error,dbsi);
	OUTPUT:
	RETVAL



Wise2_AlnBlock *
AlnBlock_from_Protein_estwise_wrap(pro,cdna,cp,cm,ct,comp,gap,ext,is_global,rmd,alg,rm,use_syn,allN,dpri)
	Wise2_Protein * pro
	Wise2_cDNA * cdna
	Wise2_cDNAParser * cp
	Wise2_CodonMapper * cm
	Wise2_CodonTable * ct
	Wise2_CompMat * comp
	int gap
	int ext
	boolean is_global
	Wise2_RandomModelDNA * rmd
	int alg
	Wise2_RandomModel * rm
	boolean use_syn
	Probability allN
	Wise2_DPRunImpl * dpri
	CODE:
	RETVAL = Wise2_AlnBlock_from_Protein_estwise_wrap(pro,cdna,cp,cm,ct,comp,gap,ext,is_global,rmd,alg,rm,use_syn,allN,dpri,NULL);
	OUTPUT:
	RETVAL



Wise2_AlnBlock *
AlnBlock_from_TSM_estwise_wrap(tsm,cdna,cp,cm,ct,rmd,alg,use_syn,force_flat_insert,allN,dpri)
	Wise2_ThreeStateModel * tsm
	Wise2_cDNA * cdna
	Wise2_cDNAParser * cp
	Wise2_CodonMapper * cm
	Wise2_CodonTable * ct
	Wise2_RandomModelDNA * rmd
	int alg
	boolean use_syn
	boolean force_flat_insert
	Probability allN
	Wise2_DPRunImpl * dpri
	CODE:
	RETVAL = Wise2_AlnBlock_from_TSM_estwise_wrap(tsm,cdna,cp,cm,ct,rmd,alg,use_syn,force_flat_insert,allN,dpri,NULL);
	OUTPUT:
	RETVAL



int
alg_estwrap_from_string(str)
	char * str
	CODE:
	RETVAL = Wise2_alg_estwrap_from_string(str);
	OUTPUT:
	RETVAL



