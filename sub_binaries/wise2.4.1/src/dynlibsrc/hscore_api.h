

/* Functions that create, manipulate or act on Hscore
 *
 * Wise2_minimum_score_Hscore
 * Wise2_maximum_score_Hscore
 * Wise2_sort_Hscore_by_score
 * Wise2_length_datascore_Hscore
 * Wise2_get_datascore_Hscore
 * Wise2_get_score_Hscore
 * Wise2_get_evalue_Hscore
 * Wise2_basic_show_Hscore
 * Wise2_hard_link_Hscore
 * Wise2_Hscore_alloc_std
 * Wise2_free_Hscore [destructor]
 *
 */



/* Functions that create, manipulate or act on DataScore
 *
 * Wise2_hard_link_DataScore
 * Wise2_DataScore_alloc
 * Wise2_replace_query_DataScore
 * Wise2_access_query_DataScore
 * Wise2_replace_target_DataScore
 * Wise2_access_target_DataScore
 * Wise2_replace_score_DataScore
 * Wise2_access_score_DataScore
 * Wise2_replace_evalue_DataScore
 * Wise2_access_evalue_DataScore
 * Wise2_free_DataScore [destructor]
 *
 */



/* Functions that create, manipulate or act on DataEntry
 *
 * Wise2_hard_link_DataEntry
 * Wise2_DataEntry_alloc
 * Wise2_replace_name_DataEntry
 * Wise2_access_name_DataEntry
 * Wise2_replace_is_reversed_DataEntry
 * Wise2_access_is_reversed_DataEntry
 * Wise2_free_DataEntry [destructor]
 *
 */



/* Helper functions in the module
 *
 * Wise2_std_score_Hscore
 *

/* API for object Hscore */
/* Function:  Wise2_minimum_score_Hscore(hs)
 *
 * Descrip:    gets the minimum score from Hscore
 *
 *
 * Arg:        hs           Undocumented argument [Wise2_Hscore *]
 *
 * Returns Undocumented return value [int]
 *
 */
int Wise2_minimum_score_Hscore( Wise2_Hscore * hs);

/* Function:  Wise2_maximum_score_Hscore(hs)
 *
 * Descrip:    gets the maximum score from Hscore
 *
 *
 * Arg:        hs           Undocumented argument [Wise2_Hscore *]
 *
 * Returns Undocumented return value [int]
 *
 */
int Wise2_maximum_score_Hscore( Wise2_Hscore * hs);

/* Function:  Wise2_sort_Hscore_by_score(hs)
 *
 * Descrip:    As it says, sorts the high score by its score
 *
 *
 * Arg:        hs           Hscore to be sorted [Wise2_Hscore *]
 *
 * Returns Undocumented return value [void]
 *
 */
void Wise2_sort_Hscore_by_score( Wise2_Hscore * hs);

/* Function:  Wise2_length_datascore_Hscore(obj)
 *
 * Descrip:    Returns the number of datascores in the hscore
 *             structure
 *
 *
 * Arg:        obj          Hscore object [Wise2_Hscore *]
 *
 * Returns Undocumented return value [int]
 *
 */
int Wise2_length_datascore_Hscore( Wise2_Hscore * obj);

/* Function:  Wise2_get_datascore_Hscore(hs,i)
 *
 * Descrip:    Returns the specific datascore held at this
 *             position.
 *
 *             This requires a considerable amount of memory
 *             duplication, so please dont process all your
 *             results by looping through this.
 *
 *
 * Arg:        hs           Hscore object [Wise2_Hscore *]
 * Arg:        i            position to be read [int]
 *
 * Returns New datascore object [Wise2_DataScore *]
 *
 */
Wise2_DataScore * Wise2_get_datascore_Hscore( Wise2_Hscore * hs,int i);

/* Function:  Wise2_get_score_Hscore(hs,i)
 *
 * Descrip: No Description
 *
 * Arg:        hs           Hscore object [Wise2_Hscore *]
 * Arg:        i            position to be read [int]
 *
 * Returns score  [int]
 *
 */
int Wise2_get_score_Hscore( Wise2_Hscore * hs,int i);

/* Function:  Wise2_get_evalue_Hscore(hs,i)
 *
 * Descrip:    Returns the evalue of the specific datascore held at this position.
 *
 *
 *
 * Arg:        hs           Hscore object [Wise2_Hscore *]
 * Arg:        i            position to be read [int]
 *
 * Returns evalue  [double]
 *
 */
double Wise2_get_evalue_Hscore( Wise2_Hscore * hs,int i);

/* Function:  Wise2_basic_show_Hscore(hs,ofp)
 *
 * Descrip:    The most baby-talk showing of Hscore
 *
 *
 * Arg:        hs           Undocumented argument [Wise2_Hscore *]
 * Arg:        ofp          Undocumented argument [FILE *]
 *
 * Returns Undocumented return value [void]
 *
 */
void Wise2_basic_show_Hscore( Wise2_Hscore * hs,FILE * ofp);

/* Function:  Wise2_hard_link_Hscore(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj          Object to be hard linked [Wise2_Hscore *]
 *
 * Returns Undocumented return value [Wise2_Hscore *]
 *
 */
Wise2_Hscore * Wise2_hard_link_Hscore( Wise2_Hscore * obj);

/* Function:  Wise2_Hscore_alloc_std(void)
 *
 * Descrip:    Equivalent to Hscore_alloc_len(HscoreLISTLENGTH)
 *
 *
 *
 * Returns Undocumented return value [Wise2_Hscore *]
 *
 */
Wise2_Hscore * Wise2_Hscore_alloc_std();

/* This is the destructor function, ie, call this to free object*/
/* Function:  Wise2_free_Hscore(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj          Object that is free'd [Wise2_Hscore *]
 *
 * Returns Undocumented return value [Wise2_Hscore *]
 *
 */
Wise2_Hscore * Wise2_free_Hscore( Wise2_Hscore * obj);

/* API for object DataScore */
/* Function:  Wise2_hard_link_DataScore(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj          Object to be hard linked [Wise2_DataScore *]
 *
 * Returns Undocumented return value [Wise2_DataScore *]
 *
 */
Wise2_DataScore * Wise2_hard_link_DataScore( Wise2_DataScore * obj);

/* Function:  Wise2_DataScore_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Returns Undocumented return value [Wise2_DataScore *]
 *
 */
Wise2_DataScore * Wise2_DataScore_alloc();

/* Function:  Wise2_replace_query_DataScore(obj,query)
 *
 * Descrip:    Replace member variable query
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_DataScore *]
 * Arg:        query        New value of the variable [Wise2_DataEntry *]
 *
 * Returns member variable query [boolean]
 *
 */
boolean Wise2_replace_query_DataScore( Wise2_DataScore * obj,Wise2_DataEntry * query);

/* Function:  Wise2_access_query_DataScore(obj)
 *
 * Descrip:    Access member variable query
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_DataScore *]
 *
 * Returns member variable query [Wise2_DataEntry *]
 *
 */
Wise2_DataEntry * Wise2_access_query_DataScore( Wise2_DataScore * obj);

/* Function:  Wise2_replace_target_DataScore(obj,target)
 *
 * Descrip:    Replace member variable target
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_DataScore *]
 * Arg:        target       New value of the variable [Wise2_DataEntry *]
 *
 * Returns member variable target [boolean]
 *
 */
boolean Wise2_replace_target_DataScore( Wise2_DataScore * obj,Wise2_DataEntry * target);

/* Function:  Wise2_access_target_DataScore(obj)
 *
 * Descrip:    Access member variable target
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_DataScore *]
 *
 * Returns member variable target [Wise2_DataEntry *]
 *
 */
Wise2_DataEntry * Wise2_access_target_DataScore( Wise2_DataScore * obj);

/* Function:  Wise2_replace_score_DataScore(obj,score)
 *
 * Descrip:    Replace member variable score
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_DataScore *]
 * Arg:        score        New value of the variable [int]
 *
 * Returns member variable score [boolean]
 *
 */
boolean Wise2_replace_score_DataScore( Wise2_DataScore * obj,int score);

/* Function:  Wise2_access_score_DataScore(obj)
 *
 * Descrip:    Access member variable score
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_DataScore *]
 *
 * Returns member variable score [int]
 *
 */
int Wise2_access_score_DataScore( Wise2_DataScore * obj);

/* Function:  Wise2_replace_evalue_DataScore(obj,evalue)
 *
 * Descrip:    Replace member variable evalue
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_DataScore *]
 * Arg:        evalue       New value of the variable [double]
 *
 * Returns member variable evalue [boolean]
 *
 */
boolean Wise2_replace_evalue_DataScore( Wise2_DataScore * obj,double evalue);

/* Function:  Wise2_access_evalue_DataScore(obj)
 *
 * Descrip:    Access member variable evalue
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_DataScore *]
 *
 * Returns member variable evalue [double]
 *
 */
double Wise2_access_evalue_DataScore( Wise2_DataScore * obj);

/* This is the destructor function, ie, call this to free object*/
/* Function:  Wise2_free_DataScore(obj)
 *
 * Descrip:    Correctly handles destruction of a datascore
 *
 *
 * Arg:        obj          Undocumented argument [Wise2_DataScore *]
 *
 * Returns Undocumented return value [Wise2_DataScore *]
 *
 */
Wise2_DataScore * Wise2_free_DataScore( Wise2_DataScore * obj);

/* API for object DataEntry */
/* Function:  Wise2_hard_link_DataEntry(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj          Object to be hard linked [Wise2_DataEntry *]
 *
 * Returns Undocumented return value [Wise2_DataEntry *]
 *
 */
Wise2_DataEntry * Wise2_hard_link_DataEntry( Wise2_DataEntry * obj);

/* Function:  Wise2_DataEntry_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Returns Undocumented return value [Wise2_DataEntry *]
 *
 */
Wise2_DataEntry * Wise2_DataEntry_alloc();

/* Function:  Wise2_replace_name_DataEntry(obj,name)
 *
 * Descrip:    Replace member variable name
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_DataEntry *]
 * Arg:        name         New value of the variable [char *]
 *
 * Returns member variable name [boolean]
 *
 */
boolean Wise2_replace_name_DataEntry( Wise2_DataEntry * obj,char * name);

/* Function:  Wise2_access_name_DataEntry(obj)
 *
 * Descrip:    Access member variable name
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_DataEntry *]
 *
 * Returns member variable name [char *]
 *
 */
char * Wise2_access_name_DataEntry( Wise2_DataEntry * obj);

/* Function:  Wise2_replace_is_reversed_DataEntry(obj,is_reversed)
 *
 * Descrip:    Replace member variable is_reversed
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_DataEntry *]
 * Arg:        is_reversed  New value of the variable [boolean]
 *
 * Returns member variable is_reversed [boolean]
 *
 */
boolean Wise2_replace_is_reversed_DataEntry( Wise2_DataEntry * obj,boolean is_reversed);

/* Function:  Wise2_access_is_reversed_DataEntry(obj)
 *
 * Descrip:    Access member variable is_reversed
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_DataEntry *]
 *
 * Returns member variable is_reversed [boolean]
 *
 */
boolean Wise2_access_is_reversed_DataEntry( Wise2_DataEntry * obj);

/* This is the destructor function, ie, call this to free object*/
/* Function:  Wise2_free_DataEntry(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj          Object that is free'd [Wise2_DataEntry *]
 *
 * Returns Undocumented return value [Wise2_DataEntry *]
 *
 */
Wise2_DataEntry * Wise2_free_DataEntry( Wise2_DataEntry * obj);



/* These functions are not associated with an object */
/* Function:  Wise2_std_score_Hscore(cut_off,report_stagger)
 *
 * Descrip:    This gives you a standard Hscore
 *             module with a cutoff in score
 *
 *
 * Arg:        cut_off      Undocumented argument [int]
 * Arg:        report_stagger Undocumented argument [int]
 *
 * Returns Undocumented return value [Wise2_Hscore *]
 *
 */
Wise2_Hscore * Wise2_std_score_Hscore( int cut_off,int report_stagger);

