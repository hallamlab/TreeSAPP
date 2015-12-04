

/* Functions that create, manipulate or act on ComplexSequence
 *
 * Wise2_hard_link_ComplexSequence
 * Wise2_ComplexSequence_alloc
 * Wise2_replace_type_ComplexSequence
 * Wise2_access_type_ComplexSequence
 * Wise2_replace_seq_ComplexSequence
 * Wise2_access_seq_ComplexSequence
 * Wise2_free_ComplexSequence [destructor]
 *
 */



/* Functions that create, manipulate or act on ComplexSequenceEvalSet
 *
 * Wise2_hard_link_ComplexSequenceEvalSet
 * Wise2_ComplexSequenceEvalSet_alloc_std
 * Wise2_replace_type_ComplexSequenceEvalSet
 * Wise2_access_type_ComplexSequenceEvalSet
 * Wise2_replace_has_been_prepared_ComplexSequenceEvalSet
 * Wise2_access_has_been_prepared_ComplexSequenceEvalSet
 * Wise2_replace_left_window_ComplexSequenceEvalSet
 * Wise2_access_left_window_ComplexSequenceEvalSet
 * Wise2_replace_right_window_ComplexSequenceEvalSet
 * Wise2_access_right_window_ComplexSequenceEvalSet
 * Wise2_replace_left_lookback_ComplexSequenceEvalSet
 * Wise2_access_left_lookback_ComplexSequenceEvalSet
 * Wise2_free_ComplexSequenceEvalSet [destructor]
 *
 */

/* API for object ComplexSequence */
/* Function:  Wise2_hard_link_ComplexSequence(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj          Object to be hard linked [Wise2_ComplexSequence *]
 *
 * Returns Undocumented return value [Wise2_ComplexSequence *]
 *
 */
Wise2_ComplexSequence * Wise2_hard_link_ComplexSequence( Wise2_ComplexSequence * obj);

/* Function:  Wise2_ComplexSequence_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Returns Undocumented return value [Wise2_ComplexSequence *]
 *
 */
Wise2_ComplexSequence * Wise2_ComplexSequence_alloc();

/* Function:  Wise2_replace_type_ComplexSequence(obj,type)
 *
 * Descrip:    Replace member variable type
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_ComplexSequence *]
 * Arg:        type         New value of the variable [int]
 *
 * Returns member variable type [boolean]
 *
 */
boolean Wise2_replace_type_ComplexSequence( Wise2_ComplexSequence * obj,int type);

/* Function:  Wise2_access_type_ComplexSequence(obj)
 *
 * Descrip:    Access member variable type
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_ComplexSequence *]
 *
 * Returns member variable type [int]
 *
 */
int Wise2_access_type_ComplexSequence( Wise2_ComplexSequence * obj);

/* Function:  Wise2_replace_seq_ComplexSequence(obj,seq)
 *
 * Descrip:    Replace member variable seq
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_ComplexSequence *]
 * Arg:        seq          New value of the variable [Wise2_Sequence *]
 *
 * Returns member variable seq [boolean]
 *
 */
boolean Wise2_replace_seq_ComplexSequence( Wise2_ComplexSequence * obj,Wise2_Sequence * seq);

/* Function:  Wise2_access_seq_ComplexSequence(obj)
 *
 * Descrip:    Access member variable seq
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_ComplexSequence *]
 *
 * Returns member variable seq [Wise2_Sequence *]
 *
 */
Wise2_Sequence * Wise2_access_seq_ComplexSequence( Wise2_ComplexSequence * obj);

/* This is the destructor function, ie, call this to free object*/
/* Function:  Wise2_free_ComplexSequence(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj          Object that is free'd [Wise2_ComplexSequence *]
 *
 * Returns Undocumented return value [Wise2_ComplexSequence *]
 *
 */
Wise2_ComplexSequence * Wise2_free_ComplexSequence( Wise2_ComplexSequence * obj);

/* API for object ComplexSequenceEvalSet */
/* Function:  Wise2_hard_link_ComplexSequenceEvalSet(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj          Object to be hard linked [Wise2_ComplexSequenceEvalSet *]
 *
 * Returns Undocumented return value [Wise2_ComplexSequenceEvalSet *]
 *
 */
Wise2_ComplexSequenceEvalSet * Wise2_hard_link_ComplexSequenceEvalSet( Wise2_ComplexSequenceEvalSet * obj);

/* Function:  Wise2_ComplexSequenceEvalSet_alloc_std(void)
 *
 * Descrip:    Equivalent to ComplexSequenceEvalSet_alloc_len(ComplexSequenceEvalSetLISTLENGTH)
 *
 *
 *
 * Returns Undocumented return value [Wise2_ComplexSequenceEvalSet *]
 *
 */
Wise2_ComplexSequenceEvalSet * Wise2_ComplexSequenceEvalSet_alloc_std();

/* Function:  Wise2_replace_type_ComplexSequenceEvalSet(obj,type)
 *
 * Descrip:    Replace member variable type
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_ComplexSequenceEvalSet *]
 * Arg:        type         New value of the variable [int]
 *
 * Returns member variable type [boolean]
 *
 */
boolean Wise2_replace_type_ComplexSequenceEvalSet( Wise2_ComplexSequenceEvalSet * obj,int type);

/* Function:  Wise2_access_type_ComplexSequenceEvalSet(obj)
 *
 * Descrip:    Access member variable type
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_ComplexSequenceEvalSet *]
 *
 * Returns member variable type [int]
 *
 */
int Wise2_access_type_ComplexSequenceEvalSet( Wise2_ComplexSequenceEvalSet * obj);

/* Function:  Wise2_replace_has_been_prepared_ComplexSequenceEvalSet(obj,has_been_prepared)
 *
 * Descrip:    Replace member variable has_been_prepared
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_ComplexSequenceEvalSet *]
 * Arg:        has_been_prepared New value of the variable [boolean]
 *
 * Returns member variable has_been_prepared [boolean]
 *
 */
boolean Wise2_replace_has_been_prepared_ComplexSequenceEvalSet( Wise2_ComplexSequenceEvalSet * obj,boolean has_been_prepared);

/* Function:  Wise2_access_has_been_prepared_ComplexSequenceEvalSet(obj)
 *
 * Descrip:    Access member variable has_been_prepared
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_ComplexSequenceEvalSet *]
 *
 * Returns member variable has_been_prepared [boolean]
 *
 */
boolean Wise2_access_has_been_prepared_ComplexSequenceEvalSet( Wise2_ComplexSequenceEvalSet * obj);

/* Function:  Wise2_replace_left_window_ComplexSequenceEvalSet(obj,left_window)
 *
 * Descrip:    Replace member variable left_window
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_ComplexSequenceEvalSet *]
 * Arg:        left_window  New value of the variable [int]
 *
 * Returns member variable left_window [boolean]
 *
 */
boolean Wise2_replace_left_window_ComplexSequenceEvalSet( Wise2_ComplexSequenceEvalSet * obj,int left_window);

/* Function:  Wise2_access_left_window_ComplexSequenceEvalSet(obj)
 *
 * Descrip:    Access member variable left_window
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_ComplexSequenceEvalSet *]
 *
 * Returns member variable left_window [int]
 *
 */
int Wise2_access_left_window_ComplexSequenceEvalSet( Wise2_ComplexSequenceEvalSet * obj);

/* Function:  Wise2_replace_right_window_ComplexSequenceEvalSet(obj,right_window)
 *
 * Descrip:    Replace member variable right_window
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_ComplexSequenceEvalSet *]
 * Arg:        right_window New value of the variable [int]
 *
 * Returns member variable right_window [boolean]
 *
 */
boolean Wise2_replace_right_window_ComplexSequenceEvalSet( Wise2_ComplexSequenceEvalSet * obj,int right_window);

/* Function:  Wise2_access_right_window_ComplexSequenceEvalSet(obj)
 *
 * Descrip:    Access member variable right_window
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_ComplexSequenceEvalSet *]
 *
 * Returns member variable right_window [int]
 *
 */
int Wise2_access_right_window_ComplexSequenceEvalSet( Wise2_ComplexSequenceEvalSet * obj);

/* Function:  Wise2_replace_left_lookback_ComplexSequenceEvalSet(obj,left_lookback)
 *
 * Descrip:    Replace member variable left_lookback
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_ComplexSequenceEvalSet *]
 * Arg:        left_lookback New value of the variable [int]
 *
 * Returns member variable left_lookback [boolean]
 *
 */
boolean Wise2_replace_left_lookback_ComplexSequenceEvalSet( Wise2_ComplexSequenceEvalSet * obj,int left_lookback);

/* Function:  Wise2_access_left_lookback_ComplexSequenceEvalSet(obj)
 *
 * Descrip:    Access member variable left_lookback
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_ComplexSequenceEvalSet *]
 *
 * Returns member variable left_lookback [int]
 *
 */
int Wise2_access_left_lookback_ComplexSequenceEvalSet( Wise2_ComplexSequenceEvalSet * obj);

/* This is the destructor function, ie, call this to free object*/
/* Function:  Wise2_free_ComplexSequenceEvalSet(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj          Object that is free'd [Wise2_ComplexSequenceEvalSet *]
 *
 * Returns Undocumented return value [Wise2_ComplexSequenceEvalSet *]
 *
 */
Wise2_ComplexSequenceEvalSet * Wise2_free_ComplexSequenceEvalSet( Wise2_ComplexSequenceEvalSet * obj);

