

/* Functions that create, manipulate or act on ThreeStateModel
 *
 * Wise2_force_global_model
 * Wise2_force_weighted_local_model
 * Wise2_ThreeStateModel_from_half_bit_Sequence
 * Wise2_write_HMMer_1_7_ascii_ThreeStateModel
 * Wise2_hard_link_ThreeStateModel
 * Wise2_ThreeStateModel_alloc_std
 * Wise2_replace_name_ThreeStateModel
 * Wise2_access_name_ThreeStateModel
 * Wise2_access_unit_ThreeStateModel
 * Wise2_length_unit_ThreeStateModel
 * Wise2_flush_ThreeStateModel
 * Wise2_add_ThreeStateModel
 * Wise2_replace_alphabet_ThreeStateModel
 * Wise2_access_alphabet_ThreeStateModel
 * Wise2_replace_accession_ThreeStateModel
 * Wise2_access_accession_ThreeStateModel
 * Wise2_replace_threshold_ThreeStateModel
 * Wise2_access_threshold_ThreeStateModel
 * Wise2_replace_rm_ThreeStateModel
 * Wise2_access_rm_ThreeStateModel
 * Wise2_free_ThreeStateModel [destructor]
 *
 */



/* Functions that create, manipulate or act on ThreeStateUnit
 *
 * Wise2_hard_link_ThreeStateUnit
 * Wise2_ThreeStateUnit_alloc
 * Wise2_replace_display_char_ThreeStateUnit
 * Wise2_access_display_char_ThreeStateUnit
 * Wise2_free_ThreeStateUnit [destructor]
 *
 */



/* Helper functions in the module
 *
 * Wise2_read_HMMer_1_7_ascii_file
 * Wise2_read_HMMer_1_7_ascii
 *

/* API for object ThreeStateModel */
/* Function:  Wise2_force_global_model(tsm,prob_into_model)
 *
 * Descrip:    Makes start at position 0 and end at position end,
 *             no other positions being valid
 *
 *
 *
 * Arg:        tsm          ThreeStateModel to be 'forced' [Wise2_ThreeStateModel *]
 * Arg:        prob_into_model Probability to start the model: for true global will be 1.0 [double]
 *
 * Returns Undocumented return value [void]
 *
 */
void Wise2_force_global_model( Wise2_ThreeStateModel * tsm,double prob_into_model);

/* Function:  Wise2_force_weighted_local_model(tsm,prob_into_model,ratio_start,ratio_end)
 *
 * Descrip:    places the ratio of probability to start/end,
 *             and then distributes the rest over the start/end
 *
 *
 *
 * Arg:        tsm          ThreeStateModel to be 'forced' [Wise2_ThreeStateModel *]
 * Arg:        prob_into_model Probability to start the model: for true global will be 1.0 [double]
 * Arg:        ratio_start  ratio of prob to unit 0 to the rest (1.0 means all goes to start) [double]
 * Arg:        ratio_end    ratio of prob to unit (last) to the rest (1.0 means all goes to the end) [double]
 *
 * Returns Undocumented return value [void]
 *
 */
void Wise2_force_weighted_local_model( Wise2_ThreeStateModel * tsm,double prob_into_model,double ratio_start,double ratio_end);

/* Function:  Wise2_ThreeStateModel_from_half_bit_Sequence(pro,mat,rm,gap,ext)
 *
 * Descrip:    Makes a local three-state-model from a sequence.  this is scary
 *             hackery, assumming that the matrix is half-bits and normalising in a
 *             *very* wrong way to get "probabilities" out.
 *
 *             Works though
 *
 *
 * Arg:        pro          protein sequence [Wise2_Protein *]
 * Arg:        mat          comparison matrix to use [Wise2_CompMat *]
 * Arg:        rm           random model which you assumme the matrix was built with [Wise2_RandomModel *]
 * Arg:        gap          gap open penalty [int]
 * Arg:        ext          gap ext penalty [int]
 *
 * Returns Undocumented return value [Wise2_ThreeStateModel *]
 *
 */
Wise2_ThreeStateModel * Wise2_ThreeStateModel_from_half_bit_Sequence( Wise2_Protein * pro,Wise2_CompMat * mat,Wise2_RandomModel * rm,int gap,int ext);

/* Function:  Wise2_write_HMMer_1_7_ascii_ThreeStateModel(tsm,ofp)
 *
 * Descrip:    writes a HMMer version 1.7 (also ok with 1.8) file
 *
 *
 * Arg:        tsm          Undocumented argument [Wise2_ThreeStateModel *]
 * Arg:        ofp          Undocumented argument [FILE *]
 *
 * Returns Undocumented return value [void]
 *
 */
void Wise2_write_HMMer_1_7_ascii_ThreeStateModel( Wise2_ThreeStateModel * tsm,FILE * ofp);

/* Function:  Wise2_hard_link_ThreeStateModel(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj          Object to be hard linked [Wise2_ThreeStateModel *]
 *
 * Returns Undocumented return value [Wise2_ThreeStateModel *]
 *
 */
Wise2_ThreeStateModel * Wise2_hard_link_ThreeStateModel( Wise2_ThreeStateModel * obj);

/* Function:  Wise2_ThreeStateModel_alloc_std(void)
 *
 * Descrip:    Equivalent to ThreeStateModel_alloc_len(ThreeStateModelLISTLENGTH)
 *
 *
 *
 * Returns Undocumented return value [Wise2_ThreeStateModel *]
 *
 */
Wise2_ThreeStateModel * Wise2_ThreeStateModel_alloc_std();

/* Function:  Wise2_replace_name_ThreeStateModel(obj,name)
 *
 * Descrip:    Replace member variable name
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_ThreeStateModel *]
 * Arg:        name         New value of the variable [char *]
 *
 * Returns member variable name [boolean]
 *
 */
boolean Wise2_replace_name_ThreeStateModel( Wise2_ThreeStateModel * obj,char * name);

/* Function:  Wise2_access_name_ThreeStateModel(obj)
 *
 * Descrip:    Access member variable name
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_ThreeStateModel *]
 *
 * Returns member variable name [char *]
 *
 */
char * Wise2_access_name_ThreeStateModel( Wise2_ThreeStateModel * obj);

/* Function:  Wise2_access_unit_ThreeStateModel(obj,i)
 *
 * Descrip:    Access members stored in the unit list
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the list [Wise2_ThreeStateModel *]
 * Arg:        i            Position in the list [int]
 *
 * Returns Element of the list [Wise2_ThreeStateUnit *]
 *
 */
Wise2_ThreeStateUnit * Wise2_access_unit_ThreeStateModel( Wise2_ThreeStateModel * obj,int i);

/* Function:  Wise2_length_unit_ThreeStateModel(obj)
 *
 * Descrip:    discover the length of the list
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the list [Wise2_ThreeStateModel *]
 *
 * Returns length of the list [int]
 *
 */
int Wise2_length_unit_ThreeStateModel( Wise2_ThreeStateModel * obj);

/* Function:  Wise2_flush_ThreeStateModel(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj          Object which contains the list  [Wise2_ThreeStateModel *]
 *
 * Returns Undocumented return value [int]
 *
 */
int Wise2_flush_ThreeStateModel( Wise2_ThreeStateModel * obj);

/* Function:  Wise2_add_ThreeStateModel(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj          Object which contains the list [Wise2_ThreeStateModel *]
 * Arg:        add          Object to add to the list [Wise2_ThreeStateUnit *]
 *
 * Returns Undocumented return value [boolean]
 *
 */
boolean Wise2_add_ThreeStateModel( Wise2_ThreeStateModel * obj,Wise2_ThreeStateUnit * add);

/* Function:  Wise2_replace_alphabet_ThreeStateModel(obj,alphabet)
 *
 * Descrip:    Replace member variable alphabet
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_ThreeStateModel *]
 * Arg:        alphabet     New value of the variable [char *]
 *
 * Returns member variable alphabet [boolean]
 *
 */
boolean Wise2_replace_alphabet_ThreeStateModel( Wise2_ThreeStateModel * obj,char * alphabet);

/* Function:  Wise2_access_alphabet_ThreeStateModel(obj)
 *
 * Descrip:    Access member variable alphabet
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_ThreeStateModel *]
 *
 * Returns member variable alphabet [char *]
 *
 */
char * Wise2_access_alphabet_ThreeStateModel( Wise2_ThreeStateModel * obj);

/* Function:  Wise2_replace_accession_ThreeStateModel(obj,accession)
 *
 * Descrip:    Replace member variable accession
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_ThreeStateModel *]
 * Arg:        accession    New value of the variable [char *]
 *
 * Returns member variable accession [boolean]
 *
 */
boolean Wise2_replace_accession_ThreeStateModel( Wise2_ThreeStateModel * obj,char * accession);

/* Function:  Wise2_access_accession_ThreeStateModel(obj)
 *
 * Descrip:    Access member variable accession
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_ThreeStateModel *]
 *
 * Returns member variable accession [char *]
 *
 */
char * Wise2_access_accession_ThreeStateModel( Wise2_ThreeStateModel * obj);

/* Function:  Wise2_replace_threshold_ThreeStateModel(obj,threshold)
 *
 * Descrip:    Replace member variable threshold
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_ThreeStateModel *]
 * Arg:        threshold    New value of the variable [double]
 *
 * Returns member variable threshold [boolean]
 *
 */
boolean Wise2_replace_threshold_ThreeStateModel( Wise2_ThreeStateModel * obj,double threshold);

/* Function:  Wise2_access_threshold_ThreeStateModel(obj)
 *
 * Descrip:    Access member variable threshold
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_ThreeStateModel *]
 *
 * Returns member variable threshold [double]
 *
 */
double Wise2_access_threshold_ThreeStateModel( Wise2_ThreeStateModel * obj);

/* Function:  Wise2_replace_rm_ThreeStateModel(obj,rm)
 *
 * Descrip:    Replace member variable rm
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_ThreeStateModel *]
 * Arg:        rm           New value of the variable [Wise2_RandomModel *]
 *
 * Returns member variable rm [boolean]
 *
 */
boolean Wise2_replace_rm_ThreeStateModel( Wise2_ThreeStateModel * obj,Wise2_RandomModel * rm);

/* Function:  Wise2_access_rm_ThreeStateModel(obj)
 *
 * Descrip:    Access member variable rm
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_ThreeStateModel *]
 *
 * Returns member variable rm [Wise2_RandomModel *]
 *
 */
Wise2_RandomModel * Wise2_access_rm_ThreeStateModel( Wise2_ThreeStateModel * obj);

/* This is the destructor function, ie, call this to free object*/
/* Function:  Wise2_free_ThreeStateModel(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj          Object that is free'd [Wise2_ThreeStateModel *]
 *
 * Returns Undocumented return value [Wise2_ThreeStateModel *]
 *
 */
Wise2_ThreeStateModel * Wise2_free_ThreeStateModel( Wise2_ThreeStateModel * obj);

/* API for object ThreeStateUnit */
/* Function:  Wise2_hard_link_ThreeStateUnit(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj          Object to be hard linked [Wise2_ThreeStateUnit *]
 *
 * Returns Undocumented return value [Wise2_ThreeStateUnit *]
 *
 */
Wise2_ThreeStateUnit * Wise2_hard_link_ThreeStateUnit( Wise2_ThreeStateUnit * obj);

/* Function:  Wise2_ThreeStateUnit_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Returns Undocumented return value [Wise2_ThreeStateUnit *]
 *
 */
Wise2_ThreeStateUnit * Wise2_ThreeStateUnit_alloc();

/* Function:  Wise2_replace_display_char_ThreeStateUnit(obj,display_char)
 *
 * Descrip:    Replace member variable display_char
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_ThreeStateUnit *]
 * Arg:        display_char New value of the variable [char]
 *
 * Returns member variable display_char [boolean]
 *
 */
boolean Wise2_replace_display_char_ThreeStateUnit( Wise2_ThreeStateUnit * obj,char display_char);

/* Function:  Wise2_access_display_char_ThreeStateUnit(obj)
 *
 * Descrip:    Access member variable display_char
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_ThreeStateUnit *]
 *
 * Returns member variable display_char [char]
 *
 */
char Wise2_access_display_char_ThreeStateUnit( Wise2_ThreeStateUnit * obj);

/* This is the destructor function, ie, call this to free object*/
/* Function:  Wise2_free_ThreeStateUnit(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj          Object that is free'd [Wise2_ThreeStateUnit *]
 *
 * Returns Undocumented return value [Wise2_ThreeStateUnit *]
 *
 */
Wise2_ThreeStateUnit * Wise2_free_ThreeStateUnit( Wise2_ThreeStateUnit * obj);



/* These functions are not associated with an object */
/* Function:  Wise2_read_HMMer_1_7_ascii_file(filename)
 *
 * Descrip:    reads a HMMer ascii version 1.7 (1.8) file from filename.
 *
 *
 *
 * Arg:        filename     the name fo the hmmer file [char *]
 *
 * Returns Undocumented return value [Wise2_ThreeStateModel *]
 *
 */
Wise2_ThreeStateModel * Wise2_read_HMMer_1_7_ascii_file( char * filename);

/* Function:  Wise2_read_HMMer_1_7_ascii(ifp)
 *
 * Descrip:    Basic function to read HMMer version 1.7(1.8) files.
 *
 *
 * Arg:        ifp          Undocumented argument [FILE *]
 *
 * Returns Undocumented return value [Wise2_ThreeStateModel *]
 *
 */
Wise2_ThreeStateModel * Wise2_read_HMMer_1_7_ascii( FILE * ifp);

