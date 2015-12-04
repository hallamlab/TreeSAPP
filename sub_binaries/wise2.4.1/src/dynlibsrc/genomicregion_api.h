

/* Functions that create, manipulate or act on GenomicRegion
 *
 * Wise2_new_GenomicRegion
 * Wise2_read_EMBL_GenomicRegion_file
 * Wise2_read_EMBL_GenomicRegion
 * Wise2_add_Gene_to_GenomicRegion
 * Wise2_show_ace_GenomicRegion
 * Wise2_show_pretty_GenomicRegion
 * Wise2_write_Diana_FT_GenomicRegion
 * Wise2_write_Embl_FT_GenomicRegion
 * Wise2_hard_link_GenomicRegion
 * Wise2_GenomicRegion_alloc_std
 * Wise2_access_gene_GenomicRegion
 * Wise2_length_gene_GenomicRegion
 * Wise2_flush_GenomicRegion
 * Wise2_add_GenomicRegion
 * Wise2_replace_genomic_GenomicRegion
 * Wise2_access_genomic_GenomicRegion
 * Wise2_free_GenomicRegion [destructor]
 *
 */

/* API for object GenomicRegion */
/* Function:  Wise2_new_GenomicRegion(gen)
 *
 * Descrip:    makes a genomicregion from a genomic sequence
 *
 *
 * Arg:        gen          Undocumented argument [Wise2_Genomic *]
 *
 * Returns Undocumented return value [Wise2_GenomicRegion *]
 *
 */
Wise2_GenomicRegion * Wise2_new_GenomicRegion( Wise2_Genomic * gen);

/* Function:  Wise2_read_EMBL_GenomicRegion_file(filename)
 *
 * Descrip:    Reads in both EMBL sequence and features 
 *
 *
 * Arg:        filename     Undocumented argument [char *]
 *
 * Returns Undocumented return value [Wise2_GenomicRegion *]
 *
 */
Wise2_GenomicRegion * Wise2_read_EMBL_GenomicRegion_file( char * filename);

/* Function:  Wise2_read_EMBL_GenomicRegion(ifp)
 *
 * Descrip:    Reads in both EMBL sequence and features 
 *
 *
 * Arg:        ifp          Undocumented argument [FILE *]
 *
 * Returns Undocumented return value [Wise2_GenomicRegion *]
 *
 */
Wise2_GenomicRegion * Wise2_read_EMBL_GenomicRegion( FILE * ifp);

/* Function:  Wise2_add_Gene_to_GenomicRegion(gr,gene)
 *
 * Descrip:    adds a Gene to this GenomicRegion, making
 *             sure that it parent/son relationship is ok
 *
 *
 * Arg:        gr           GenomicRegion to be added to [Wise2_GenomicRegion *]
 * Arg:        gene         Gene to be added [Wise2_Gene *]
 *
 * Returns Undocumented return value [boolean]
 *
 */
boolean Wise2_add_Gene_to_GenomicRegion( Wise2_GenomicRegion * gr,Wise2_Gene * gene);

/* Function:  Wise2_show_ace_GenomicRegion(gr,seq_name,ofp)
 *
 * Descrip:    shows ACeDB subsequence source.
 *
 *             Assummes
 *               a only one transcript per gene
 *               b only cds exons are used
 *
 *
 * Arg:        gr           Undocumented argument [Wise2_GenomicRegion *]
 * Arg:        seq_name     Undocumented argument [char *]
 * Arg:        ofp          Undocumented argument [FILE *]
 *
 * Returns Undocumented return value [void]
 *
 */
void Wise2_show_ace_GenomicRegion( Wise2_GenomicRegion * gr,char * seq_name,FILE * ofp);

/* Function:  Wise2_show_pretty_GenomicRegion(gr,show_supporting,ofp)
 *
 * Descrip: No Description
 *
 * Arg:        gr           Undocumented argument [Wise2_GenomicRegion *]
 * Arg:        show_supporting Undocumented argument [boolean]
 * Arg:        ofp          Undocumented argument [FILE *]
 *
 * Returns Undocumented return value [void]
 *
 */
void Wise2_show_pretty_GenomicRegion( Wise2_GenomicRegion * gr,boolean show_supporting,FILE * ofp);

/* Function:  Wise2_write_Diana_FT_GenomicRegion(gr,ofp)
 *
 * Descrip:    Writes Embl feature table for diana use. Does assumme that
 *             there is only one transcript per gene and only
 *             cds exons are used
 *
 *             Output like
 *
 *                FT   misc_feature       join(100..200)
 *
 *
 * Arg:        gr           Undocumented argument [Wise2_GenomicRegion *]
 * Arg:        ofp          Undocumented argument [FILE *]
 *
 * Returns Undocumented return value [void]
 *
 */
void Wise2_write_Diana_FT_GenomicRegion( Wise2_GenomicRegion * gr,FILE * ofp);

/* Function:  Wise2_write_Embl_FT_GenomicRegion(gr,ofp)
 *
 * Descrip:    Writes Embl feature table. Does assumme that
 *             there is only one transcript per gene and only
 *             cds exons are used
 *
 *             Output like
 *
 *                FT   CDS          join(100..200)
 *
 *
 * Arg:        gr           Undocumented argument [Wise2_GenomicRegion *]
 * Arg:        ofp          Undocumented argument [FILE *]
 *
 * Returns Undocumented return value [void]
 *
 */
void Wise2_write_Embl_FT_GenomicRegion( Wise2_GenomicRegion * gr,FILE * ofp);

/* Function:  Wise2_hard_link_GenomicRegion(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj          Object to be hard linked [Wise2_GenomicRegion *]
 *
 * Returns Undocumented return value [Wise2_GenomicRegion *]
 *
 */
Wise2_GenomicRegion * Wise2_hard_link_GenomicRegion( Wise2_GenomicRegion * obj);

/* Function:  Wise2_GenomicRegion_alloc_std(void)
 *
 * Descrip:    Equivalent to GenomicRegion_alloc_len(GenomicRegionLISTLENGTH)
 *
 *
 *
 * Returns Undocumented return value [Wise2_GenomicRegion *]
 *
 */
Wise2_GenomicRegion * Wise2_GenomicRegion_alloc_std();

/* Function:  Wise2_access_gene_GenomicRegion(obj,i)
 *
 * Descrip:    Access members stored in the gene list
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the list [Wise2_GenomicRegion *]
 * Arg:        i            Position in the list [int]
 *
 * Returns Element of the list [Wise2_Gene *]
 *
 */
Wise2_Gene * Wise2_access_gene_GenomicRegion( Wise2_GenomicRegion * obj,int i);

/* Function:  Wise2_length_gene_GenomicRegion(obj)
 *
 * Descrip:    discover the length of the list
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the list [Wise2_GenomicRegion *]
 *
 * Returns length of the list [int]
 *
 */
int Wise2_length_gene_GenomicRegion( Wise2_GenomicRegion * obj);

/* Function:  Wise2_flush_GenomicRegion(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj          Object which contains the list  [Wise2_GenomicRegion *]
 *
 * Returns Undocumented return value [int]
 *
 */
int Wise2_flush_GenomicRegion( Wise2_GenomicRegion * obj);

/* Function:  Wise2_add_GenomicRegion(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj          Object which contains the list [Wise2_GenomicRegion *]
 * Arg:        add          Object to add to the list [Wise2_Gene *]
 *
 * Returns Undocumented return value [boolean]
 *
 */
boolean Wise2_add_GenomicRegion( Wise2_GenomicRegion * obj,Wise2_Gene * add);

/* Function:  Wise2_replace_genomic_GenomicRegion(obj,genomic)
 *
 * Descrip:    Replace member variable genomic
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_GenomicRegion *]
 * Arg:        genomic      New value of the variable [Wise2_Genomic *]
 *
 * Returns member variable genomic [boolean]
 *
 */
boolean Wise2_replace_genomic_GenomicRegion( Wise2_GenomicRegion * obj,Wise2_Genomic * genomic);

/* Function:  Wise2_access_genomic_GenomicRegion(obj)
 *
 * Descrip:    Access member variable genomic
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_GenomicRegion *]
 *
 * Returns member variable genomic [Wise2_Genomic *]
 *
 */
Wise2_Genomic * Wise2_access_genomic_GenomicRegion( Wise2_GenomicRegion * obj);

/* This is the destructor function, ie, call this to free object*/
/* Function:  Wise2_free_GenomicRegion(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj          Object that is free'd [Wise2_GenomicRegion *]
 *
 * Returns Undocumented return value [Wise2_GenomicRegion *]
 *
 */
Wise2_GenomicRegion * Wise2_free_GenomicRegion( Wise2_GenomicRegion * obj);

