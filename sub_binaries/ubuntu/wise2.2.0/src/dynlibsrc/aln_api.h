

/* Functions that create, manipulate or act on AlnBlock
 *
 * Wise2_dump_ascii_AlnBlock
 * Wise2_hard_link_AlnBlock
 * Wise2_AlnBlock_alloc_std
 * Wise2_replace_start_AlnBlock
 * Wise2_access_start_AlnBlock
 * Wise2_access_seq_AlnBlock
 * Wise2_length_seq_AlnBlock
 * Wise2_flush_AlnBlock
 * Wise2_add_AlnBlock
 * Wise2_replace_length_AlnBlock
 * Wise2_access_length_AlnBlock
 * Wise2_replace_score_AlnBlock
 * Wise2_access_score_AlnBlock
 * Wise2_free_AlnBlock [destructor]
 *
 */



/* Functions that create, manipulate or act on AlnColumn
 *
 * Wise2_at_end_AlnColumn
 * Wise2_hard_link_AlnColumn
 * Wise2_AlnColumn_alloc_std
 * Wise2_access_alu_AlnColumn
 * Wise2_length_alu_AlnColumn
 * Wise2_flush_AlnColumn
 * Wise2_add_AlnColumn
 * Wise2_replace_next_AlnColumn
 * Wise2_access_next_AlnColumn
 * Wise2_free_AlnColumn [destructor]
 *
 */



/* Functions that create, manipulate or act on AlnUnit
 *
 * Wise2_bio_start_AlnUnit
 * Wise2_bio_end_AlnUnit
 * Wise2_hard_link_AlnUnit
 * Wise2_AlnUnit_alloc
 * Wise2_replace_start_AlnUnit
 * Wise2_access_start_AlnUnit
 * Wise2_replace_end_AlnUnit
 * Wise2_access_end_AlnUnit
 * Wise2_replace_label_AlnUnit
 * Wise2_access_label_AlnUnit
 * Wise2_replace_text_label_AlnUnit
 * Wise2_access_text_label_AlnUnit
 * Wise2_replace_next_AlnUnit
 * Wise2_access_next_AlnUnit
 * Wise2_replace_in_column_AlnUnit
 * Wise2_access_in_column_AlnUnit
 * Wise2_replace_seq_AlnUnit
 * Wise2_access_seq_AlnUnit
 * Wise2_free_AlnUnit [destructor]
 *
 */



/* Functions that create, manipulate or act on AlnSequence
 *
 * Wise2_hard_link_AlnSequence
 * Wise2_AlnSequence_alloc
 * Wise2_replace_start_AlnSequence
 * Wise2_access_start_AlnSequence
 * Wise2_replace_data_type_AlnSequence
 * Wise2_access_data_type_AlnSequence
 * Wise2_replace_data_AlnSequence
 * Wise2_access_data_AlnSequence
 * Wise2_replace_bio_start_AlnSequence
 * Wise2_access_bio_start_AlnSequence
 * Wise2_replace_bio_end_AlnSequence
 * Wise2_access_bio_end_AlnSequence
 * Wise2_free_AlnSequence [destructor]
 *
 */

/* API for object AlnBlock */
/* Function:  Wise2_dump_ascii_AlnBlock(alb,ofp)
 *
 * Descrip:    Dumps the alignment in rereadable ascii form.
 *
 *             Not really for human consumption
 *
 *
 * Arg:        alb          AlnBlock to dump [Wise2_AlnBlock *]
 * Arg:        ofp          File stream to dump to [FILE *]
 *
 * Returns Undocumented return value [void]
 *
 */
void Wise2_dump_ascii_AlnBlock( Wise2_AlnBlock * alb,FILE * ofp);

/* Function:  Wise2_hard_link_AlnBlock(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj          Object to be hard linked [Wise2_AlnBlock *]
 *
 * Returns Undocumented return value [Wise2_AlnBlock *]
 *
 */
Wise2_AlnBlock * Wise2_hard_link_AlnBlock( Wise2_AlnBlock * obj);

/* Function:  Wise2_AlnBlock_alloc_std(void)
 *
 * Descrip:    Equivalent to AlnBlock_alloc_len(AlnBlockLISTLENGTH)
 *
 *
 *
 * Returns Undocumented return value [Wise2_AlnBlock *]
 *
 */
Wise2_AlnBlock * Wise2_AlnBlock_alloc_std();

/* Function:  Wise2_replace_start_AlnBlock(obj,start)
 *
 * Descrip:    Replace member variable start
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_AlnBlock *]
 * Arg:        start        New value of the variable [Wise2_AlnColumn *]
 *
 * Returns member variable start [boolean]
 *
 */
boolean Wise2_replace_start_AlnBlock( Wise2_AlnBlock * obj,Wise2_AlnColumn * start);

/* Function:  Wise2_access_start_AlnBlock(obj)
 *
 * Descrip:    Access member variable start
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_AlnBlock *]
 *
 * Returns member variable start [Wise2_AlnColumn *]
 *
 */
Wise2_AlnColumn * Wise2_access_start_AlnBlock( Wise2_AlnBlock * obj);

/* Function:  Wise2_access_seq_AlnBlock(obj,i)
 *
 * Descrip:    Access members stored in the seq list
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the list [Wise2_AlnBlock *]
 * Arg:        i            Position in the list [int]
 *
 * Returns Element of the list [Wise2_AlnSequence *]
 *
 */
Wise2_AlnSequence * Wise2_access_seq_AlnBlock( Wise2_AlnBlock * obj,int i);

/* Function:  Wise2_length_seq_AlnBlock(obj)
 *
 * Descrip:    discover the length of the list
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the list [Wise2_AlnBlock *]
 *
 * Returns length of the list [int]
 *
 */
int Wise2_length_seq_AlnBlock( Wise2_AlnBlock * obj);

/* Function:  Wise2_flush_AlnBlock(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj          Object which contains the list  [Wise2_AlnBlock *]
 *
 * Returns Undocumented return value [int]
 *
 */
int Wise2_flush_AlnBlock( Wise2_AlnBlock * obj);

/* Function:  Wise2_add_AlnBlock(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj          Object which contains the list [Wise2_AlnBlock *]
 * Arg:        add          Object to add to the list [Wise2_AlnSequence *]
 *
 * Returns Undocumented return value [boolean]
 *
 */
boolean Wise2_add_AlnBlock( Wise2_AlnBlock * obj,Wise2_AlnSequence * add);

/* Function:  Wise2_replace_length_AlnBlock(obj,length)
 *
 * Descrip:    Replace member variable length
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_AlnBlock *]
 * Arg:        length       New value of the variable [int]
 *
 * Returns member variable length [boolean]
 *
 */
boolean Wise2_replace_length_AlnBlock( Wise2_AlnBlock * obj,int length);

/* Function:  Wise2_access_length_AlnBlock(obj)
 *
 * Descrip:    Access member variable length
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_AlnBlock *]
 *
 * Returns member variable length [int]
 *
 */
int Wise2_access_length_AlnBlock( Wise2_AlnBlock * obj);

/* Function:  Wise2_replace_score_AlnBlock(obj,score)
 *
 * Descrip:    Replace member variable score
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_AlnBlock *]
 * Arg:        score        New value of the variable [int]
 *
 * Returns member variable score [boolean]
 *
 */
boolean Wise2_replace_score_AlnBlock( Wise2_AlnBlock * obj,int score);

/* Function:  Wise2_access_score_AlnBlock(obj)
 *
 * Descrip:    Access member variable score
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_AlnBlock *]
 *
 * Returns member variable score [int]
 *
 */
int Wise2_access_score_AlnBlock( Wise2_AlnBlock * obj);

/* This is the destructor function, ie, call this to free object*/
/* Function:  Wise2_free_AlnBlock(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj          Object that is free'd [Wise2_AlnBlock *]
 *
 * Returns Undocumented return value [Wise2_AlnBlock *]
 *
 */
Wise2_AlnBlock * Wise2_free_AlnBlock( Wise2_AlnBlock * obj);

/* API for object AlnColumn */
/* Function:  Wise2_at_end_AlnColumn(alc)
 *
 * Descrip:    This tells you whether the AlnColumn is at the
 *             end without passing NULL's around
 *
 *
 *
 * Arg:        alc          AlnColumn [Wise2_AlnColumn *]
 *
 * Returns Undocumented return value [boolean]
 *
 */
boolean Wise2_at_end_AlnColumn( Wise2_AlnColumn * alc);

/* Function:  Wise2_hard_link_AlnColumn(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj          Object to be hard linked [Wise2_AlnColumn *]
 *
 * Returns Undocumented return value [Wise2_AlnColumn *]
 *
 */
Wise2_AlnColumn * Wise2_hard_link_AlnColumn( Wise2_AlnColumn * obj);

/* Function:  Wise2_AlnColumn_alloc_std(void)
 *
 * Descrip:    Equivalent to AlnColumn_alloc_len(AlnColumnLISTLENGTH)
 *
 *
 *
 * Returns Undocumented return value [Wise2_AlnColumn *]
 *
 */
Wise2_AlnColumn * Wise2_AlnColumn_alloc_std();

/* Function:  Wise2_access_alu_AlnColumn(obj,i)
 *
 * Descrip:    Access members stored in the alu list
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the list [Wise2_AlnColumn *]
 * Arg:        i            Position in the list [int]
 *
 * Returns Element of the list [Wise2_AlnUnit *]
 *
 */
Wise2_AlnUnit * Wise2_access_alu_AlnColumn( Wise2_AlnColumn * obj,int i);

/* Function:  Wise2_length_alu_AlnColumn(obj)
 *
 * Descrip:    discover the length of the list
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the list [Wise2_AlnColumn *]
 *
 * Returns length of the list [int]
 *
 */
int Wise2_length_alu_AlnColumn( Wise2_AlnColumn * obj);

/* Function:  Wise2_flush_AlnColumn(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj          Object which contains the list  [Wise2_AlnColumn *]
 *
 * Returns Undocumented return value [int]
 *
 */
int Wise2_flush_AlnColumn( Wise2_AlnColumn * obj);

/* Function:  Wise2_add_AlnColumn(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj          Object which contains the list [Wise2_AlnColumn *]
 * Arg:        add          Object to add to the list [Wise2_AlnUnit *]
 *
 * Returns Undocumented return value [boolean]
 *
 */
boolean Wise2_add_AlnColumn( Wise2_AlnColumn * obj,Wise2_AlnUnit * add);

/* Function:  Wise2_replace_next_AlnColumn(obj,next)
 *
 * Descrip:    Replace member variable next
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_AlnColumn *]
 * Arg:        next         New value of the variable [Wise2_AlnColumn *]
 *
 * Returns member variable next [boolean]
 *
 */
boolean Wise2_replace_next_AlnColumn( Wise2_AlnColumn * obj,Wise2_AlnColumn * next);

/* Function:  Wise2_access_next_AlnColumn(obj)
 *
 * Descrip:    Access member variable next
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_AlnColumn *]
 *
 * Returns member variable next [Wise2_AlnColumn *]
 *
 */
Wise2_AlnColumn * Wise2_access_next_AlnColumn( Wise2_AlnColumn * obj);

/* This is the destructor function, ie, call this to free object*/
/* Function:  Wise2_free_AlnColumn(obj)
 *
 * Descrip:    Specilased deconstructor needed because
 *             of linked list nature of the data structure
 *
 *
 * Arg:        obj          Undocumented argument [Wise2_AlnColumn *]
 *
 * Returns Undocumented return value [Wise2_AlnColumn *]
 *
 */
Wise2_AlnColumn * Wise2_free_AlnColumn( Wise2_AlnColumn * obj);

/* API for object AlnUnit */
/* Function:  Wise2_bio_start_AlnUnit(alu)
 *
 * Descrip:    Tells the bio-coordinate of the
 *             start point of this alnunit
 *
 *
 * Arg:        alu          Undocumented argument [Wise2_AlnUnit *]
 *
 * Returns Undocumented return value [int]
 *
 */
int Wise2_bio_start_AlnUnit( Wise2_AlnUnit * alu);

/* Function:  Wise2_bio_end_AlnUnit(alu)
 *
 * Descrip:    Tells the bio-coordinate of the
 *             end point of this alnunit
 *
 *
 * Arg:        alu          Undocumented argument [Wise2_AlnUnit *]
 *
 * Returns Undocumented return value [int]
 *
 */
int Wise2_bio_end_AlnUnit( Wise2_AlnUnit * alu);

/* Function:  Wise2_hard_link_AlnUnit(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj          Object to be hard linked [Wise2_AlnUnit *]
 *
 * Returns Undocumented return value [Wise2_AlnUnit *]
 *
 */
Wise2_AlnUnit * Wise2_hard_link_AlnUnit( Wise2_AlnUnit * obj);

/* Function:  Wise2_AlnUnit_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Returns Undocumented return value [Wise2_AlnUnit *]
 *
 */
Wise2_AlnUnit * Wise2_AlnUnit_alloc();

/* Function:  Wise2_replace_start_AlnUnit(obj,start)
 *
 * Descrip:    Replace member variable start
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_AlnUnit *]
 * Arg:        start        New value of the variable [int]
 *
 * Returns member variable start [boolean]
 *
 */
boolean Wise2_replace_start_AlnUnit( Wise2_AlnUnit * obj,int start);

/* Function:  Wise2_access_start_AlnUnit(obj)
 *
 * Descrip:    Access member variable start
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_AlnUnit *]
 *
 * Returns member variable start [int]
 *
 */
int Wise2_access_start_AlnUnit( Wise2_AlnUnit * obj);

/* Function:  Wise2_replace_end_AlnUnit(obj,end)
 *
 * Descrip:    Replace member variable end
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_AlnUnit *]
 * Arg:        end          New value of the variable [int]
 *
 * Returns member variable end [boolean]
 *
 */
boolean Wise2_replace_end_AlnUnit( Wise2_AlnUnit * obj,int end);

/* Function:  Wise2_access_end_AlnUnit(obj)
 *
 * Descrip:    Access member variable end
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_AlnUnit *]
 *
 * Returns member variable end [int]
 *
 */
int Wise2_access_end_AlnUnit( Wise2_AlnUnit * obj);

/* Function:  Wise2_replace_label_AlnUnit(obj,label)
 *
 * Descrip:    Replace member variable label
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_AlnUnit *]
 * Arg:        label        New value of the variable [int]
 *
 * Returns member variable label [boolean]
 *
 */
boolean Wise2_replace_label_AlnUnit( Wise2_AlnUnit * obj,int label);

/* Function:  Wise2_access_label_AlnUnit(obj)
 *
 * Descrip:    Access member variable label
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_AlnUnit *]
 *
 * Returns member variable label [int]
 *
 */
int Wise2_access_label_AlnUnit( Wise2_AlnUnit * obj);

/* Function:  Wise2_replace_text_label_AlnUnit(obj,text_label)
 *
 * Descrip:    Replace member variable text_label
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_AlnUnit *]
 * Arg:        text_label   New value of the variable [char *]
 *
 * Returns member variable text_label [boolean]
 *
 */
boolean Wise2_replace_text_label_AlnUnit( Wise2_AlnUnit * obj,char * text_label);

/* Function:  Wise2_access_text_label_AlnUnit(obj)
 *
 * Descrip:    Access member variable text_label
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_AlnUnit *]
 *
 * Returns member variable text_label [char *]
 *
 */
char * Wise2_access_text_label_AlnUnit( Wise2_AlnUnit * obj);

/* Function:  Wise2_replace_next_AlnUnit(obj,next)
 *
 * Descrip:    Replace member variable next
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_AlnUnit *]
 * Arg:        next         New value of the variable [Wise2_AlnUnit *]
 *
 * Returns member variable next [boolean]
 *
 */
boolean Wise2_replace_next_AlnUnit( Wise2_AlnUnit * obj,Wise2_AlnUnit * next);

/* Function:  Wise2_access_next_AlnUnit(obj)
 *
 * Descrip:    Access member variable next
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_AlnUnit *]
 *
 * Returns member variable next [Wise2_AlnUnit *]
 *
 */
Wise2_AlnUnit * Wise2_access_next_AlnUnit( Wise2_AlnUnit * obj);

/* Function:  Wise2_replace_in_column_AlnUnit(obj,in_column)
 *
 * Descrip:    Replace member variable in_column
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_AlnUnit *]
 * Arg:        in_column    New value of the variable [boolean]
 *
 * Returns member variable in_column [boolean]
 *
 */
boolean Wise2_replace_in_column_AlnUnit( Wise2_AlnUnit * obj,boolean in_column);

/* Function:  Wise2_access_in_column_AlnUnit(obj)
 *
 * Descrip:    Access member variable in_column
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_AlnUnit *]
 *
 * Returns member variable in_column [boolean]
 *
 */
boolean Wise2_access_in_column_AlnUnit( Wise2_AlnUnit * obj);

/* Function:  Wise2_replace_seq_AlnUnit(obj,seq)
 *
 * Descrip:    Replace member variable seq
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_AlnUnit *]
 * Arg:        seq          New value of the variable [Wise2_AlnSequence *]
 *
 * Returns member variable seq [boolean]
 *
 */
boolean Wise2_replace_seq_AlnUnit( Wise2_AlnUnit * obj,Wise2_AlnSequence * seq);

/* Function:  Wise2_access_seq_AlnUnit(obj)
 *
 * Descrip:    Access member variable seq
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_AlnUnit *]
 *
 * Returns member variable seq [Wise2_AlnSequence *]
 *
 */
Wise2_AlnSequence * Wise2_access_seq_AlnUnit( Wise2_AlnUnit * obj);

/* This is the destructor function, ie, call this to free object*/
/* Function:  Wise2_free_AlnUnit(obj)
 *
 * Descrip:    Specilased deconstructor needed because
 *             of linked list nature of the data structure
 *
 *
 * Arg:        obj          Undocumented argument [Wise2_AlnUnit *]
 *
 * Returns Undocumented return value [Wise2_AlnUnit *]
 *
 */
Wise2_AlnUnit * Wise2_free_AlnUnit( Wise2_AlnUnit * obj);

/* API for object AlnSequence */
/* Function:  Wise2_hard_link_AlnSequence(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj          Object to be hard linked [Wise2_AlnSequence *]
 *
 * Returns Undocumented return value [Wise2_AlnSequence *]
 *
 */
Wise2_AlnSequence * Wise2_hard_link_AlnSequence( Wise2_AlnSequence * obj);

/* Function:  Wise2_AlnSequence_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Returns Undocumented return value [Wise2_AlnSequence *]
 *
 */
Wise2_AlnSequence * Wise2_AlnSequence_alloc();

/* Function:  Wise2_replace_start_AlnSequence(obj,start)
 *
 * Descrip:    Replace member variable start
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_AlnSequence *]
 * Arg:        start        New value of the variable [Wise2_AlnUnit *]
 *
 * Returns member variable start [boolean]
 *
 */
boolean Wise2_replace_start_AlnSequence( Wise2_AlnSequence * obj,Wise2_AlnUnit * start);

/* Function:  Wise2_access_start_AlnSequence(obj)
 *
 * Descrip:    Access member variable start
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_AlnSequence *]
 *
 * Returns member variable start [Wise2_AlnUnit *]
 *
 */
Wise2_AlnUnit * Wise2_access_start_AlnSequence( Wise2_AlnSequence * obj);

/* Function:  Wise2_replace_data_type_AlnSequence(obj,data_type)
 *
 * Descrip:    Replace member variable data_type
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_AlnSequence *]
 * Arg:        data_type    New value of the variable [int]
 *
 * Returns member variable data_type [boolean]
 *
 */
boolean Wise2_replace_data_type_AlnSequence( Wise2_AlnSequence * obj,int data_type);

/* Function:  Wise2_access_data_type_AlnSequence(obj)
 *
 * Descrip:    Access member variable data_type
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_AlnSequence *]
 *
 * Returns member variable data_type [int]
 *
 */
int Wise2_access_data_type_AlnSequence( Wise2_AlnSequence * obj);

/* Function:  Wise2_replace_data_AlnSequence(obj,data)
 *
 * Descrip:    Replace member variable data
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_AlnSequence *]
 * Arg:        data         New value of the variable [void *]
 *
 * Returns member variable data [boolean]
 *
 */
boolean Wise2_replace_data_AlnSequence( Wise2_AlnSequence * obj,void * data);

/* Function:  Wise2_access_data_AlnSequence(obj)
 *
 * Descrip:    Access member variable data
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_AlnSequence *]
 *
 * Returns member variable data [void *]
 *
 */
void * Wise2_access_data_AlnSequence( Wise2_AlnSequence * obj);

/* Function:  Wise2_replace_bio_start_AlnSequence(obj,bio_start)
 *
 * Descrip:    Replace member variable bio_start
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_AlnSequence *]
 * Arg:        bio_start    New value of the variable [int]
 *
 * Returns member variable bio_start [boolean]
 *
 */
boolean Wise2_replace_bio_start_AlnSequence( Wise2_AlnSequence * obj,int bio_start);

/* Function:  Wise2_access_bio_start_AlnSequence(obj)
 *
 * Descrip:    Access member variable bio_start
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_AlnSequence *]
 *
 * Returns member variable bio_start [int]
 *
 */
int Wise2_access_bio_start_AlnSequence( Wise2_AlnSequence * obj);

/* Function:  Wise2_replace_bio_end_AlnSequence(obj,bio_end)
 *
 * Descrip:    Replace member variable bio_end
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_AlnSequence *]
 * Arg:        bio_end      New value of the variable [int]
 *
 * Returns member variable bio_end [boolean]
 *
 */
boolean Wise2_replace_bio_end_AlnSequence( Wise2_AlnSequence * obj,int bio_end);

/* Function:  Wise2_access_bio_end_AlnSequence(obj)
 *
 * Descrip:    Access member variable bio_end
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_AlnSequence *]
 *
 * Returns member variable bio_end [int]
 *
 */
int Wise2_access_bio_end_AlnSequence( Wise2_AlnSequence * obj);

/* This is the destructor function, ie, call this to free object*/
/* Function:  Wise2_free_AlnSequence(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj          Object that is free'd [Wise2_AlnSequence *]
 *
 * Returns Undocumented return value [Wise2_AlnSequence *]
 *
 */
Wise2_AlnSequence * Wise2_free_AlnSequence( Wise2_AlnSequence * obj);

