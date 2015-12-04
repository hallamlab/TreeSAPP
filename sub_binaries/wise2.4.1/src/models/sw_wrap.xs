

MODULE = Wise2 PACKAGE = Wise2

Wise2_AlnBlock *
Align_strings_ProteinSmithWaterman(one,two,comp,gap,ext,dpenv,dpri)
	char * one
	char * two
	Wise2_CompMat * comp
	int gap
	int ext
	Wise2_DPEnvelope * dpenv
	Wise2_DPRunImpl * dpri
	CODE:
	RETVAL = Wise2_Align_strings_ProteinSmithWaterman(one,two,comp,gap,ext,dpenv,dpri);
	OUTPUT:
	RETVAL



Wise2_AlnBlock *
Align_Sequences_ProteinSmithWaterman(one,two,comp,gap,ext,dpenv,dpri)
	Wise2_Sequence * one
	Wise2_Sequence * two
	Wise2_CompMat * comp
	int gap
	int ext
	Wise2_DPEnvelope * dpenv
	Wise2_DPRunImpl * dpri
	CODE:
	RETVAL = Wise2_Align_Sequences_ProteinSmithWaterman(one,two,comp,gap,ext,dpenv,dpri);
	OUTPUT:
	RETVAL



Wise2_AlnBlock *
Align_Proteins_SmithWaterman(one,two,comp,gap,ext,dpenv,dpri)
	Wise2_Protein * one
	Wise2_Protein * two
	Wise2_CompMat * comp
	int gap
	int ext
	Wise2_DPEnvelope * dpenv
	Wise2_DPRunImpl * dpri
	CODE:
	RETVAL = Wise2_Align_Proteins_SmithWaterman(one,two,comp,gap,ext,dpenv,dpri);
	OUTPUT:
	RETVAL



Wise2_AlnBlock *
Align_Proteins_ABC(one,two,comp,a,b,c,dpenv,dpri)
	Wise2_Protein * one
	Wise2_Protein * two
	Wise2_CompMat * comp
	int a
	int b
	int c
	Wise2_DPEnvelope * dpenv
	Wise2_DPRunImpl * dpri
	CODE:
	RETVAL = Wise2_Align_Proteins_ABC(one,two,comp,a,b,c,dpenv,dpri);
	OUTPUT:
	RETVAL



Wise2_AlnBlock *
Align_Sequences_ProteinABC(one,two,comp,a,b,c,dpenv,dpri)
	Wise2_Sequence * one
	Wise2_Sequence * two
	Wise2_CompMat * comp
	int a
	int b
	int c
	Wise2_DPEnvelope * dpenv
	Wise2_DPRunImpl * dpri
	CODE:
	RETVAL = Wise2_Align_Sequences_ProteinABC(one,two,comp,a,b,c,dpenv,dpri);
	OUTPUT:
	RETVAL



Wise2_Hscore *
Hscore_from_ProteinSW(querydb,targetdb,comp,gap,ext,bits_cutoff,report_level,die_on_error,dbsi)
	Wise2_ProteinDB* querydb
	Wise2_ProteinDB* targetdb
	Wise2_CompMat* comp
	int gap
	int ext
	double bits_cutoff
	int report_level
	boolean die_on_error
	Wise2_DBSearchImpl* dbsi
	CODE:
	RETVAL = Wise2_Hscore_from_ProteinSW(querydb,targetdb,comp,gap,ext,bits_cutoff,report_level,die_on_error,dbsi);
	OUTPUT:
	RETVAL



Wise2_Hscore *
Hscore_from_ProteinABC(querydb,targetdb,comp,a,b,c,bits_cutoff,report_level,die_on_error,dbsi)
	Wise2_ProteinDB* querydb
	Wise2_ProteinDB* targetdb
	Wise2_CompMat* comp
	int a
	int b
	int c
	double bits_cutoff
	int report_level
	boolean die_on_error
	Wise2_DBSearchImpl* dbsi
	CODE:
	RETVAL = Wise2_Hscore_from_ProteinABC(querydb,targetdb,comp,a,b,c,bits_cutoff,report_level,die_on_error,dbsi);
	OUTPUT:
	RETVAL



Wise2_Hscore *
Hscore_from_ProteinBA(querydb,targetdb,comp,bentry,bexit,bfor_trans,b_self_trans,b3exit,bits_cutoff,report_level,dbsi)
	Wise2_ProteinDB* querydb
	Wise2_ProteinDB* targetdb
	Wise2_CompMat* comp
	Score bentry
	Score bexit
	Score bfor_trans
	Score b_self_trans
	Score b3exit
	double bits_cutoff
	int report_level
	Wise2_DBSearchImpl* dbsi
	CODE:
	RETVAL = Wise2_Hscore_from_ProteinBA(querydb,targetdb,comp,bentry,bexit,bfor_trans,b_self_trans,b3exit,bits_cutoff,report_level,dbsi);
	OUTPUT:
	RETVAL



