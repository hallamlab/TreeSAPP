

/* Helper functions in the module
 *
 * Wise2_Hscore_from_TSM_estwise
 * Wise2_AlnBlock_from_Protein_estwise_wrap
 * Wise2_AlnBlock_from_TSM_estwise_wrap
 * Wise2_alg_estwrap_from_string
 *



/* These functions are not associated with an object */
/* Function:  Wise2_Hscore_from_TSM_estwise(tdb,cdb,cp,cm,rmd,use_syn,alg,bits_cutoff,allN,flat_insert,report_level,die_on_error,dbsi)
 *
 * Descrip:    Runs a database search for the estwise set
 *             of algorithms
 *
 *
 * Arg:        tdb          a three state model database [Wise2_ThreeStateDB *]
 * Arg:        cdb          a dna sequence database [Wise2_cDNADB *]
 * Arg:        cp           the codon parser for this comparison [Wise2_cDNAParser *]
 * Arg:        cm           the codon mapper for this comparison [Wise2_CodonMapper *]
 * Arg:        rmd          random model used for the dna sequence comparison [Wise2_RandomModelDNA *]
 * Arg:        use_syn      whether a synchronous coding model should be used or not [boolean]
 * Arg:        alg          algorithm to use [int]
 * Arg:        bits_cutoff  Undocumented argument [double]
 * Arg:        allN         Undocumented argument [Probability]
 * Arg:        flat_insert  Undocumented argument [boolean]
 * Arg:        report_level Undocumented argument [int]
 * Arg:        die_on_error if true, dies if there is an error [boolean]
 * Arg:        dbsi         Undocumented argument [Wise2_DBSearchImpl *]
 *
 * Returns a newly allocated Hscore structure of the search [Wise2_Hscore *]
 *
 */
Wise2_Hscore * Wise2_Hscore_from_TSM_estwise( Wise2_ThreeStateDB * tdb,Wise2_cDNADB * cdb,Wise2_cDNAParser * cp,Wise2_CodonMapper * cm,Wise2_RandomModelDNA * rmd,boolean use_syn,int alg,double bits_cutoff,Probability allN,boolean flat_insert,int report_level,boolean die_on_error,Wise2_DBSearchImpl * dbsi);

/* Function:  Wise2_AlnBlock_from_Protein_estwise_wrap(pro,cdna,cp,cm,ct,comp,gap,ext,is_global,rmd,alg,rm,use_syn,allN,dpri,palpoi)
 *
 * Descrip:    This function is the guts for the est single alignment
 *             mode.
 *
 *             It uses /AlnBlock_from_TSM_estwise_wrap for the
 *             heavy part of the call
 *
 *
 * Arg:        pro          protein to be used in the comparison [Wise2_Protein *]
 * Arg:        cdna         cdna to be compared to [Wise2_cDNA *]
 * Arg:        cp           cdna parser indicating insertion deletion probabilities [Wise2_cDNAParser *]
 * Arg:        cm           codon mapper indicating substitution errors etc [Wise2_CodonMapper *]
 * Arg:        ct           codon table for the codon->amino acid mappings [Wise2_CodonTable *]
 * Arg:        comp         comparison matrix to use [Wise2_CompMat *]
 * Arg:        gap          gap penalty [int]
 * Arg:        ext          extension penalty [int]
 * Arg:        is_global    if true, global start-end in protein is used [boolean]
 * Arg:        rmd          random model of dna to use [Wise2_RandomModelDNA *]
 * Arg:        alg          est algorithm type to use [int]
 * Arg:        rm           random protein model for use with compmat [Wise2_RandomModel *]
 * Arg:        use_syn      if true, uses a synchronous coding model [boolean]
 * Arg:        allN         Undocumented argument [Probability]
 * Arg:        dpri         Undocumented argument [Wise2_DPRunImpl *]
 * Arg:        palpoi       the raw packed alignment output if wanted [Wise2_PackAln **]
 *
 * Returns Undocumented return value [Wise2_AlnBlock *]
 *
 */
Wise2_AlnBlock * Wise2_AlnBlock_from_Protein_estwise_wrap( Wise2_Protein * pro,Wise2_cDNA * cdna,Wise2_cDNAParser * cp,Wise2_CodonMapper * cm,Wise2_CodonTable * ct,Wise2_CompMat * comp,int gap,int ext,boolean is_global,Wise2_RandomModelDNA * rmd,int alg,Wise2_RandomModel * rm,boolean use_syn,Probability allN,Wise2_DPRunImpl * dpri,Wise2_PackAln ** palpoi);

/* Function:  Wise2_AlnBlock_from_TSM_estwise_wrap(tsm,cdna,cp,cm,ct,rmd,alg,use_syn,force_flat_insert,allN,dpri,palpoi)
 *
 * Descrip:    This function is the basic wrap for Protein models
 *             vs cDNA sequences.
 *
 *
 * Arg:        tsm          threestatemodel to be compared to the dna [Wise2_ThreeStateModel *]
 * Arg:        cdna         cdna to be compared to [Wise2_cDNA *]
 * Arg:        cp           cdna parser indicating insertion deletion probabilities [Wise2_cDNAParser *]
 * Arg:        cm           codon mapper indicating substitution errors etc [Wise2_CodonMapper *]
 * Arg:        ct           codon table for the codon->amino acid mappings [Wise2_CodonTable *]
 * Arg:        rmd          random model of dna to use [Wise2_RandomModelDNA *]
 * Arg:        alg          est algorithm type to use [int]
 * Arg:        use_syn      if true, uses a synchronous coding model [boolean]
 * Arg:        force_flat_insert Undocumented argument [boolean]
 * Arg:        allN         Undocumented argument [Probability]
 * Arg:        dpri         Undocumented argument [Wise2_DPRunImpl *]
 * Arg:        palpoi       the raw packed alignment output if wanted [Wise2_PackAln **]
 *
 * Returns Undocumented return value [Wise2_AlnBlock *]
 *
 */
Wise2_AlnBlock * Wise2_AlnBlock_from_TSM_estwise_wrap( Wise2_ThreeStateModel * tsm,Wise2_cDNA * cdna,Wise2_cDNAParser * cp,Wise2_CodonMapper * cm,Wise2_CodonTable * ct,Wise2_RandomModelDNA * rmd,int alg,boolean use_syn,boolean force_flat_insert,Probability allN,Wise2_DPRunImpl * dpri,Wise2_PackAln ** palpoi);

/* Function:  Wise2_alg_estwrap_from_string(str)
 *
 * Descrip:    This function returns the algorithm type
 *             for an est search from the string
 *
 *
 * Arg:        str          Undocumented argument [char *]
 *
 * Returns Undocumented return value [int]
 *
 */
int Wise2_alg_estwrap_from_string( char * str);

