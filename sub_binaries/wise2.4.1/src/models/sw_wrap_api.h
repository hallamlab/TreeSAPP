

/* Helper functions in the module
 *
 * Wise2_Align_strings_ProteinSmithWaterman
 * Wise2_Align_Sequences_ProteinSmithWaterman
 * Wise2_Align_Proteins_SmithWaterman
 * Wise2_Align_Proteins_ABC
 * Wise2_Align_Sequences_ProteinABC
 * Wise2_Hscore_from_ProteinSW
 * Wise2_Hscore_from_ProteinABC
 * Wise2_Hscore_from_ProteinBA
 *



/* These functions are not associated with an object */
/* Function:  Wise2_Align_strings_ProteinSmithWaterman(one,two,comp,gap,ext,dpenv,dpri)
 *
 * Descrip:    This is the most *stupidly* abstracted view of two sequences
 *             getting aligned, being two strings.
 *
 *             It would be much better if you used Sequence objects or Protein
 *             objects to carry the proteins.
 *
 *
 * Arg:        one          string of the first sequence [char *]
 * Arg:        two          string of the second sequence [char *]
 * Arg:        comp         Comparison Matrix [Wise2_CompMat *]
 * Arg:        gap          gap penalty [int]
 * Arg:        ext          extension penalty [int]
 * Arg:        dpenv        Undocumented argument [Wise2_DPEnvelope *]
 * Arg:        dpri         Undocumented argument [Wise2_DPRunImpl *]
 *
 * Returns Undocumented return value [Wise2_AlnBlock *]
 *
 */
Wise2_AlnBlock * Wise2_Align_strings_ProteinSmithWaterman( char * one,char * two,Wise2_CompMat * comp,int gap,int ext,Wise2_DPEnvelope * dpenv,Wise2_DPRunImpl * dpri);

/* Function:  Wise2_Align_Sequences_ProteinSmithWaterman(one,two,comp,gap,ext,dpenv,dpri)
 *
 * Descrip:    This function is a mid-level abstraction of
 *             comparing two sequences, which could be
 *             generic types (eg DNA!). This is tested
 *             for and warnings are given but the alignment
 *             is still calculated. To prevent this test
 *             warning either make sure the Sequence types
 *             are set to PROTEIN or, better still, use the
 *             high level abstraction Align_Proteins_SmithWaterman
 *
 *             Otherwise this performs a standard smith waterman
 *             protein alignment...
 *
 *             To display the alignment use  write_pretty_seq_align
 *
 *
 * Arg:        one          First sequence to compare [Wise2_Sequence *]
 * Arg:        two          Second sequecne to compare [Wise2_Sequence *]
 * Arg:        comp         Comparison matrix to use [Wise2_CompMat *]
 * Arg:        gap          gap penalty. Must be negative or 0 [int]
 * Arg:        ext          ext penalty. Must be negative or 0 [int]
 * Arg:        dpenv        Undocumented argument [Wise2_DPEnvelope *]
 * Arg:        dpri         Undocumented argument [Wise2_DPRunImpl *]
 *
 * Returns new AlnBlock structure representing the alignment [Wise2_AlnBlock *]
 *
 */
Wise2_AlnBlock * Wise2_Align_Sequences_ProteinSmithWaterman( Wise2_Sequence * one,Wise2_Sequence * two,Wise2_CompMat * comp,int gap,int ext,Wise2_DPEnvelope * dpenv,Wise2_DPRunImpl * dpri);

/* Function:  Wise2_Align_Proteins_SmithWaterman(one,two,comp,gap,ext,dpenv,dpri)
 *
 * Descrip:    This is the most correct way of aligning two Proteins,
 *             using Protein objects, which can be assummed to be
 *             proteins with no objections
 *
 *             To display the alignment use write_pretty_Protein_align
 *
 *
 *
 * Arg:        one          Protein to align [Wise2_Protein *]
 * Arg:        two          Protein to align [Wise2_Protein *]
 * Arg:        comp         Comparison Matrix [Wise2_CompMat *]
 * Arg:        gap          gap penalty [int]
 * Arg:        ext          extension penalty [int]
 * Arg:        dpenv        Undocumented argument [Wise2_DPEnvelope *]
 * Arg:        dpri         Undocumented argument [Wise2_DPRunImpl *]
 *
 * Returns Undocumented return value [Wise2_AlnBlock *]
 *
 */
Wise2_AlnBlock * Wise2_Align_Proteins_SmithWaterman( Wise2_Protein * one,Wise2_Protein * two,Wise2_CompMat * comp,int gap,int ext,Wise2_DPEnvelope * dpenv,Wise2_DPRunImpl * dpri);

/* Function:  Wise2_Align_Proteins_ABC(one,two,comp,a,b,c,dpenv,dpri)
 *
 * Descrip:    Analogous to Align_Proteins_SmithWaterman for ABC model
 *
 *
 * Arg:        one          protein to align [Wise2_Protein *]
 * Arg:        two          protein to align [Wise2_Protein *]
 * Arg:        comp         comparison matrix [Wise2_CompMat *]
 * Arg:        a            generalized affine gap cost a [int]
 * Arg:        b            generalized affine gap cost b [int]
 * Arg:        c            generalized affine gap cost c [int]
 * Arg:        dpenv        Undocumented argument [Wise2_DPEnvelope *]
 * Arg:        dpri         Undocumented argument [Wise2_DPRunImpl *]
 *
 * Returns Undocumented return value [Wise2_AlnBlock *]
 *
 */
Wise2_AlnBlock * Wise2_Align_Proteins_ABC( Wise2_Protein * one,Wise2_Protein * two,Wise2_CompMat * comp,int a,int b,int c,Wise2_DPEnvelope * dpenv,Wise2_DPRunImpl * dpri);

/* Function:  Wise2_Align_Sequences_ProteinABC(one,two,comp,a,b,c,dpenv,dpri)
 *
 * Descrip:    Align_Sequences_ProteinABC
 *             this function is analogous to Align_Sequences_ProteinSmithWaterman
 *             but using the abc model
 *
 *
 * Arg:        one          Sequence to align [Wise2_Sequence *]
 * Arg:        two          Sequence to align [Wise2_Sequence *]
 * Arg:        comp         Comparison Matrix [Wise2_CompMat *]
 * Arg:        a            genearlized affine gap cost  [int]
 * Arg:        b            genearlized affine gap cost  [int]
 * Arg:        c            genearlized affine gap cost  [int]
 * Arg:        dpenv        Undocumented argument [Wise2_DPEnvelope *]
 * Arg:        dpri         Undocumented argument [Wise2_DPRunImpl *]
 *
 * Returns Undocumented return value [Wise2_AlnBlock *]
 *
 */
Wise2_AlnBlock * Wise2_Align_Sequences_ProteinABC( Wise2_Sequence * one,Wise2_Sequence * two,Wise2_CompMat * comp,int a,int b,int c,Wise2_DPEnvelope * dpenv,Wise2_DPRunImpl * dpri);

/* Function:  Wise2_Hscore_from_ProteinSW(querydb,targetdb,comp,gap,ext,bits_cutoff,report_level,die_on_error,dbsi)
 *
 * Descrip:    Runs a database psw search 
 *
 *
 * Arg:        querydb      query database  [Wise2_ProteinDB*]
 * Arg:        targetdb     target database [Wise2_ProteinDB*]
 * Arg:        comp         comparison matrix [Wise2_CompMat*]
 * Arg:        gap          gap penalty [int]
 * Arg:        ext          extension penalty [int]
 * Arg:        bits_cutoff   [double]
 * Arg:        report_level  [int]
 * Arg:        die_on_error  [boolean]
 * Arg:        dbsi          [Wise2_DBSearchImpl*]
 *
 * Returns Undocumented return value [Wise2_Hscore *]
 *
 */
Wise2_Hscore * Wise2_Hscore_from_ProteinSW( Wise2_ProteinDB* querydb,Wise2_ProteinDB* targetdb,Wise2_CompMat* comp,int gap,int ext,double bits_cutoff,int report_level,boolean die_on_error,Wise2_DBSearchImpl* dbsi);

/* Function:  Wise2_Hscore_from_ProteinABC(querydb,targetdb,comp,a,b,c,bits_cutoff,report_level,die_on_error,dbsi)
 *
 * Descrip:    Runs a database abc search 
 *
 *
 * Arg:        querydb      query database  [Wise2_ProteinDB*]
 * Arg:        targetdb     target database [Wise2_ProteinDB*]
 * Arg:        comp         comparison matrix [Wise2_CompMat*]
 * Arg:        a            generalized affine gap cost a [int]
 * Arg:        b            generalized affine gap cost b [int]
 * Arg:        c            generalized affine gap cost c [int]
 * Arg:        bits_cutoff   [double]
 * Arg:        report_level  [int]
 * Arg:        die_on_error  [boolean]
 * Arg:        dbsi          [Wise2_DBSearchImpl*]
 *
 * Returns Undocumented return value [Wise2_Hscore *]
 *
 */
Wise2_Hscore * Wise2_Hscore_from_ProteinABC( Wise2_ProteinDB* querydb,Wise2_ProteinDB* targetdb,Wise2_CompMat* comp,int a,int b,int c,double bits_cutoff,int report_level,boolean die_on_error,Wise2_DBSearchImpl* dbsi);

/* Function:  Wise2_Hscore_from_ProteinBA(querydb,targetdb,comp,bentry,bexit,bfor_trans,b_self_trans,b3exit,bits_cutoff,report_level,dbsi)
 *
 * Descrip:    Runs a database pba search
 *
 *
 * Arg:        querydb      query database [Wise2_ProteinDB*]
 * Arg:        targetdb     target database [Wise2_ProteinDB*]
 * Arg:        comp         comparison matrix [Wise2_CompMat*]
 * Arg:        bentry        [Score]
 * Arg:        bexit         [Score]
 * Arg:        bfor_trans    [Score]
 * Arg:        b_self_trans  [Score]
 * Arg:        b3exit        [Score]
 * Arg:        bits_cutoff   [double]
 * Arg:        report_level  [int]
 * Arg:        dbsi          [Wise2_DBSearchImpl*]
 *
 * Returns Undocumented return value [Wise2_Hscore *]
 *
 */
Wise2_Hscore * Wise2_Hscore_from_ProteinBA( Wise2_ProteinDB* querydb,Wise2_ProteinDB* targetdb,Wise2_CompMat* comp,Score bentry,Score bexit,Score bfor_trans,Score b_self_trans,Score b3exit,double bits_cutoff,int report_level,Wise2_DBSearchImpl* dbsi);

