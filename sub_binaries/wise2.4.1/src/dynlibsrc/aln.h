#ifndef DYNAMITEalnHEADERFILE
#define DYNAMITEalnHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "wisebase.h"


#define AlnColumnLISTLENGTH 32
#define AlnBlockLISTLENGTH   32
#define AlnBlockListLISTLENGTH   32

#define AlnUnitSCORENUMBER 8

#ifndef DYNAMITE_DEFINED_AlnColumn
typedef struct Wise2_AlnColumn Wise2_AlnColumn;
#define AlnColumn Wise2_AlnColumn
#define DYNAMITE_DEFINED_AlnColumn
#endif

#ifndef DYNAMITE_DEFINED_AlnUnit
typedef struct Wise2_AlnUnit Wise2_AlnUnit;
#define AlnUnit Wise2_AlnUnit
#define DYNAMITE_DEFINED_AlnUnit
#endif

#ifndef DYNAMITE_DEFINED_AlnSequence
typedef struct Wise2_AlnSequence Wise2_AlnSequence;
#define AlnSequence Wise2_AlnSequence
#define DYNAMITE_DEFINED_AlnSequence
#endif

/* Object AlnUnit
 *
 * Descrip: This is the basic unit of the label alignment.
 *        It describes a single mark-up over one sequence:
 *        being a start, an end and a text_label.
 *
 *
 */
struct Wise2_AlnUnit {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int start;  /*  start position in the sequence */ 
    int end;    /*  end position in the sequence */ 
    int label;  /*  not used */ 
    char * text_label;  /*  text label of this position */ 
    AlnUnit * next; /*  next AlnUnit in this sequence */ 
    int score[AlnUnitSCORENUMBER];  /*  a series of scores for this position. */ 
    boolean in_column;  /*  not used  */ 
    AlnSequence * seq;   
    } ;  
/* AlnUnit defined */ 
#ifndef DYNAMITE_DEFINED_AlnUnit
typedef struct Wise2_AlnUnit Wise2_AlnUnit;
#define AlnUnit Wise2_AlnUnit
#define DYNAMITE_DEFINED_AlnUnit
#endif


/* Object AlnColumn
 *
 * Descrip: This is a coupling of AlnUnits from different sequences.
 *        Each AlnUnit is meant to represent *the equivalent* piece
 *        of biological information in some sense (ie, they are
 *        alignmed), even though quite possibly they are very 
 *        different types of information
 *
 *
 */
struct Wise2_AlnColumn {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    AlnUnit ** alu; /*  list of the AlnUnits in this column */ 
    int len;/* len for above alu  */ 
    int maxlen; /* maxlen for above alu */ 
    AlnColumn * next;   /*  the next AlnColumn in this block */ 
    } ;  
/* AlnColumn defined */ 
#ifndef DYNAMITE_DEFINED_AlnColumn
typedef struct Wise2_AlnColumn Wise2_AlnColumn;
#define AlnColumn Wise2_AlnColumn
#define DYNAMITE_DEFINED_AlnColumn
#endif


/* Object AlnSequence
 *
 * Descrip: Each Sequence in an AlnBlock is represented by one of these, and
 *        in many ways this is an orthogonal way of looking at the alignment
 *        than the AlnColumns. If you look at the alignment from one 
 *        AlnSequence you will just see the individual labels on this 
 *        sequence
 *
 *
 */
struct Wise2_AlnSequence {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    AlnUnit * start;    /*  the first AlnUnit of this sequence */ 
    int data_type;  /*  not used */ 
    void * data;    /*  not used - don't use! */ 
    int bio_start;  /*  start of this sequence in its 'bio' coordinates */ 
    int bio_end;    /*  end of this sequence in its 'bio' coordinates */ 
    } ;  
/* AlnSequence defined */ 
#ifndef DYNAMITE_DEFINED_AlnSequence
typedef struct Wise2_AlnSequence Wise2_AlnSequence;
#define AlnSequence Wise2_AlnSequence
#define DYNAMITE_DEFINED_AlnSequence
#endif


/* Object AlnBlock
 *
 * Descrip: AlnBlock is the main representation of alignments from Dynamite. Each
 *        AlnBlock represents any number of 'sequences', of any type, which share
 *        things in common. The alignment is represented by a series of /AlnColumns 
 *        (linked list) in which each AlnColumn has a series of AlnUnits, each 
 *        Unit being a start/end/text_label triple. Alternatively, one can see
 *        each sequence in isolation, and not ask what it is aligned to, but rather
 *        what labels it has on it. 
 *
 *
 */
struct Wise2_AlnBlock {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    AlnColumn * start;  /*  the first AlnColumn in the alignment */ 
    AlnSequence ** seq; /*  a list of AlnSequences in the alignment */ 
    int len;/* len for above seq  */ 
    int maxlen; /* maxlen for above seq */ 
    int length; /*  not used  */ 
    int score;  /*  not used */ 
    } ;  
/* AlnBlock defined */ 
#ifndef DYNAMITE_DEFINED_AlnBlock
typedef struct Wise2_AlnBlock Wise2_AlnBlock;
#define AlnBlock Wise2_AlnBlock
#define DYNAMITE_DEFINED_AlnBlock
#endif


/* Object AlnBlockList
 *
 * Descrip: No Description
 *
 */
struct Wise2_AlnBlockList {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    AlnBlock ** alb;     
    int len;/* len for above alb  */ 
    int maxlen; /* maxlen for above alb */ 
    } ;  
/* AlnBlockList defined */ 
#ifndef DYNAMITE_DEFINED_AlnBlockList
typedef struct Wise2_AlnBlockList Wise2_AlnBlockList;
#define AlnBlockList Wise2_AlnBlockList
#define DYNAMITE_DEFINED_AlnBlockList
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  collapsed_AlnBlock(alb,comp_row)
 *
 * Descrip:    This function builds a new "collapsed" AlnBlock
 *             on the similarity of one AlnSeq row 
 *
 *
 * Arg:             alb [UNKN ] Undocumented argument [AlnBlock *]
 * Arg:        comp_row [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [AlnBlock *]
 *
 */
AlnBlock * Wise2_collapsed_AlnBlock(AlnBlock * alb,int comp_row);
#define collapsed_AlnBlock Wise2_collapsed_AlnBlock


/* Function:  add_to_anchored_AlnBlock(growing,add)
 *
 * Descrip:    This function assummes that the first AlnSequence
 *             in the current and the adding AlnBlock are the same,
 *             and that you want to add the second AlnSequence of the
 *             second alnblock to the first, anchored on the first
 *             sequence.
 *
 *             In other words, this builds up an anchored alignment, on
 *             the first sequence
 *
 *             This AlnBlock, like many others, consumes the sequence
 *             in the second alnblock, unsetting things so that it is
 *             valid. This makes for some pretty hairy coding.
 *
 *
 * Arg:        growing [UNKN ] Undocumented argument [AlnBlock *]
 * Arg:            add [UNKN ] Undocumented argument [AlnBlock *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_to_anchored_AlnBlock(AlnBlock * growing,AlnBlock * add);
#define add_to_anchored_AlnBlock Wise2_add_to_anchored_AlnBlock


/* Function:  single_unit_AlnBlock(len,label)
 *
 * Descrip:    This function makes a new AlnBlock of length len
 *             and depth one, with each Block having the given
 *             label and start = 0 end = start +1;
 *
 *             It starts with a alu going from -1 to 0
 *
 *
 * Arg:          len [UNKN ] Undocumented argument [int]
 * Arg:        label [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [AlnBlock *]
 *
 */
AlnBlock * Wise2_single_unit_AlnBlock(int len,char * label);
#define single_unit_AlnBlock Wise2_single_unit_AlnBlock


/* Function:  split_AlnBlock(alb,is_spacer_column)
 *
 * Descrip:    This function splits an AlnBlock into
 *             separate alignments (stored in the resulting AlnBlockList).
 *
 *             The alb is split wherever there is a column that returns
 *             true to the is_spacer_column function, discarding these
 *             columns
 *
 *             This function completely destroys the AlnBlock object that
 *             is passed in, but to make sure that API functions dont
 *             get confused, the alb that is passed in is simply stripped
 *             of its AlnColumn information (so it is still a valid alb,
 *             if empty). But - beware - all the alncolumns might or might
 *             not be there, so dont pass in albs and hold on to anything
 *             inside them!
 *
 *
 * Arg:                     alb [UNKN ] Undocumented argument [AlnBlock *]
 * Arg:        is_spacer_column [UNKN ] Undocumented argument [NullString]
 *
 * Return [UNKN ]  Undocumented return value [AlnBlockList *]
 *
 */
AlnBlockList * Wise2_split_AlnBlock(AlnBlock * alb,boolean (*is_spacer_column)(const AlnColumn * alc));
#define split_AlnBlock Wise2_split_AlnBlock


/* Function:  score_line_from_AlnBlock(alb,seqno)
 *
 * Descrip:    gets the score out for a particular alb sequence
 *             line
 *
 *
 * Arg:          alb [UNKN ] Undocumented argument [AlnBlock *]
 * Arg:        seqno [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_score_line_from_AlnBlock(AlnBlock * alb,int seqno);
#define score_line_from_AlnBlock Wise2_score_line_from_AlnBlock


/* Function:  at_end_AlnColumn(alc)
 *
 * Descrip:    This tells you whether the AlnColumn is at the
 *             end without passing NULL's around
 *
 *
 *
 * Arg:        alc [READ ] AlnColumn [AlnColumn *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_at_end_AlnColumn(AlnColumn * alc);
#define at_end_AlnColumn Wise2_at_end_AlnColumn


/* Function:  bio_start_AlnUnit(alu)
 *
 * Descrip:    Tells the bio-coordinate of the
 *             start point of this alnunit
 *
 *
 * Arg:        alu [UNKN ] Undocumented argument [AlnUnit *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_bio_start_AlnUnit(AlnUnit * alu);
#define bio_start_AlnUnit Wise2_bio_start_AlnUnit


/* Function:  bio_end_AlnUnit(alu)
 *
 * Descrip:    Tells the bio-coordinate of the
 *             end point of this alnunit
 *
 *
 * Arg:        alu [UNKN ] Undocumented argument [AlnUnit *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_bio_end_AlnUnit(AlnUnit * alu);
#define bio_end_AlnUnit Wise2_bio_end_AlnUnit


/* Function:  swallow_AlnColumn_multiple(master,eaten,comp_func)
 *
 * Descrip:    This function will 'swallow' any number of AlnColumns as long
 *             as the comparison function of the labels match (the basic
 *             comp function would be something like strcmp(a,b) == 0 ? TRUE : FALSE)
 *             The columns are 'swallowed' into master and come from eaten. (these
 *             columns could be in the same linked list, though it only makes sense
 *             if the master is before the eaten).
 *
 *             It returns the first column that it could not swallow.
 *
 *             you use this to collapse regions of the label alignment.
 *
 *
 * Arg:           master [UNKN ] column which will eat other columns [AlnColumn *]
 * Arg:            eaten [UNKN ] column which will be consumed [AlnColumn *]
 * Arg:        comp_func [FUNCP] comparison function for label set [boolean (*comp_func]
 *
 * Return [UNKN ]  Undocumented return value [AlnColumn *]
 *
 */
AlnColumn * Wise2_swallow_AlnColumn_multiple(AlnColumn * master,AlnColumn * eaten,boolean (*comp_func)(char *,char *));
#define swallow_AlnColumn_multiple Wise2_swallow_AlnColumn_multiple


/* Function:  swallow_AlnColumn_number(master,eaten,num,comp_func)
 *
 * Descrip:    Basicaly the same as /swallow_AlnColumn_mulitple but there is a maximum number
 *             of columns it will swallow
 *
 *
 * Arg:           master [UNKN ] column which will eat other columns [AlnColumn *]
 * Arg:            eaten [UNKN ] column which will be consumed [AlnColumn *]
 * Arg:              num [UNKN ] max number of columns to eat [int]
 * Arg:        comp_func [FUNCP] comparison function for label set [boolean (*comp_func]
 *
 * Return [UNKN ]  number of columns eaten [int]
 *
 */
int Wise2_swallow_AlnColumn_number(AlnColumn * master,AlnColumn * eaten,int num,boolean (*comp_func)(char *,char *));
#define swallow_AlnColumn_number Wise2_swallow_AlnColumn_number


/* Function:  swallow_AlnColumn(master,eaten,comp_func)
 *
 * Descrip:    This is the function that actually does the 'swallowing'. It will
 *             try to swallow eaten into master. If comp_func does not give us an
 *             ok (actually using /can_swallow_AlnColumn it returns FALSE. Otherwise
 *             it moves on the end of AlnColumn in master to eaten and adds the 
 *             score of eaten to master.
 *
 *
 * Arg:           master [UNKN ] column which will eat  [AlnColumn *]
 * Arg:            eaten [UNKN ] column which will dissappear into master if eatable [AlnColumn *]
 * Arg:        comp_func [FUNCP] comparison function for labels [boolean (*comp_func]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_swallow_AlnColumn(AlnColumn * master,AlnColumn * eaten,boolean (*comp_func)(char *,char *));
#define swallow_AlnColumn Wise2_swallow_AlnColumn


/* Function:  replace_and_free_AlnColumn_with_one(start,end,insert)
 *
 * Descrip:    Linked list manipulation function
 *
 *             Puts insert between start and end, and free's from start->next
 *             onwards. *Beware* if start is linked to end before calling this
 *             function thsi wil free end and everything chained to it. Think
 *             before you call this!
 *
 *
 * Arg:         start [UNKN ] Undocumented argument [AlnColumn *]
 * Arg:           end [UNKN ] Undocumented argument [AlnColumn *]
 * Arg:        insert [UNKN ] Undocumented argument [AlnColumn *]
 *
 */
void Wise2_replace_and_free_AlnColumn_with_one(AlnColumn * start,AlnColumn * end,AlnColumn * insert);
#define replace_and_free_AlnColumn_with_one Wise2_replace_and_free_AlnColumn_with_one


/* Function:  replace_AlnColumn_with_one(start,end,insert)
 *
 * Descrip:    Linked list manipulation function
 *
 *             places insert between start and end. If start/end are not
 *             continuous then it will loop out the start/end region
 *
 *
 * Arg:         start [UNKN ] Undocumented argument [AlnColumn *]
 * Arg:           end [UNKN ] Undocumented argument [AlnColumn *]
 * Arg:        insert [UNKN ] Undocumented argument [AlnColumn *]
 *
 */
void Wise2_replace_AlnColumn_with_one(AlnColumn * start,AlnColumn * end,AlnColumn * insert);
#define replace_AlnColumn_with_one Wise2_replace_AlnColumn_with_one


/* Function:  insert_AlnColumn(start,insert)
 *
 * Descrip:    Linked list manipulation function
 *
 *             places insert just after start: links insert
 *             up to what start was linked to
 *
 *
 * Arg:         start [UNKN ] Undocumented argument [AlnColumn *]
 * Arg:        insert [UNKN ] Undocumented argument [AlnColumn *]
 *
 */
void Wise2_insert_AlnColumn(AlnColumn * start,AlnColumn * insert);
#define insert_AlnColumn Wise2_insert_AlnColumn


/* Function:  go_back_n_AlnColumn(alb,start,n)
 *
 * Descrip:    Linked list movement function
 *
 *             A nasty function to reverse up a singly linked list by going to
 *             the start and coming back until you are in the current position. yuk.
 *
 *
 * Arg:          alb [UNKN ] Undocumented argument [AlnBlock *]
 * Arg:        start [UNKN ] Undocumented argument [AlnColumn *]
 * Arg:            n [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [AlnColumn *]
 *
 */
AlnColumn * Wise2_go_back_n_AlnColumn(AlnBlock * alb,AlnColumn * start,int n);
#define go_back_n_AlnColumn Wise2_go_back_n_AlnColumn


/* Function:  dump_ascii_AlnBlock(alb,ofp)
 *
 * Descrip:    Dumps the alignment in rereadable ascii form.
 *
 *             Not really for human consumption
 *
 *
 * Arg:        alb [UNKN ] AlnBlock to dump [AlnBlock *]
 * Arg:        ofp [UNKN ] File stream to dump to [FILE *]
 *
 */
void Wise2_dump_ascii_AlnBlock(AlnBlock * alb,FILE * ofp);
#define dump_ascii_AlnBlock Wise2_dump_ascii_AlnBlock


/* Function:  mapped_ascii_AlnBlockList(alb,score_to_double,*al,ofp)
 *
 * Descrip:    Shows a list of AlnBlocks with an arbitary mapping
 *             of the score to some other system
 *
 *
 * Arg:                    alb [UNKN ] AlnBlock to dump [NullString]
 * Arg:        score_to_double [UNKN ] Undocumented argument [NullString]
 * Arg:                    *al [UNKN ] Undocumented argument [AlnBlockList]
 * Arg:                    ofp [UNKN ] File stream to dump to [FILE *]
 *
 */
void Wise2_mapped_ascii_AlnBlockList(AlnBlockList *al,double (*score_to_double)(int),FILE * ofp);
#define mapped_ascii_AlnBlockList Wise2_mapped_ascii_AlnBlockList


/* Function:  mapped_ascii_AlnBlock(alb,score_to_double,do_cumlative,ofp)
 *
 * Descrip:    Shows AlnBlock with an arbitary mapping of
 *             the score to some other system. 
 *
 *
 * Arg:                    alb [UNKN ] AlnBlock to dump [AlnBlock *]
 * Arg:        score_to_double [UNKN ] Undocumented argument [NullString]
 * Arg:           do_cumlative [UNKN ] Undocumented argument [int]
 * Arg:                    ofp [UNKN ] File stream to dump to [FILE *]
 *
 */
void Wise2_mapped_ascii_AlnBlock(AlnBlock * alb,double (*score_to_double)(int),int do_cumlative,FILE * ofp);
#define mapped_ascii_AlnBlock Wise2_mapped_ascii_AlnBlock


/* Function:  read_ascii_dump_AlnBlock(ifp)
 *
 * Descrip:    Reads an ascii dumped alignment
 *
 *
 * Arg:        ifp [UNKN ] File stream to read from [FILE *]
 *
 * Return [UNKN ]  Undocumented return value [AlnBlock *]
 *
 */
AlnBlock * Wise2_read_ascii_dump_AlnBlock(FILE * ifp);
#define read_ascii_dump_AlnBlock Wise2_read_ascii_dump_AlnBlock


/* Function:  show_flat_AlnBlock(alb,ofp)
 *
 * Descrip:    Shows the AlnBlock in vaguely human
 *             readable form
 *
 *
 * Arg:        alb [UNKN ] AlnBlock to show [AlnBlock *]
 * Arg:        ofp [UNKN ] output [FILE *]
 *
 */
void Wise2_show_flat_AlnBlock(AlnBlock * alb,FILE * ofp);
#define show_flat_AlnBlock Wise2_show_flat_AlnBlock


/* Function:  get_second_end_AlnColumn(alb)
 *
 * Descrip:    Not sure if this is used!
 *
 *
 * Arg:        alb [UNKN ] Undocumented argument [AlnBlock *]
 *
 * Return [UNKN ]  Undocumented return value [AlnColumn *]
 *
 */
AlnColumn * Wise2_get_second_end_AlnColumn(AlnBlock * alb);
#define get_second_end_AlnColumn Wise2_get_second_end_AlnColumn


/* Function:  get_end_AlnColumn(alb)
 *
 * Descrip:    To get to the last AlnColumn. If this was
 *             a doubly linked list, life would be much easier
 *
 *
 * Arg:        alb [UNKN ] Undocumented argument [AlnBlock *]
 *
 * Return [UNKN ]  Undocumented return value [AlnColumn *]
 *
 */
AlnColumn * Wise2_get_end_AlnColumn(AlnBlock * alb);
#define get_end_AlnColumn Wise2_get_end_AlnColumn


/* Function:  link_AlnUnits_AlnBlock(alb)
 *
 * Descrip:    Links up all AlnUnits to their parent
 *             sequences
 *
 *
 * Arg:        alb [UNKN ] Undocumented argument [AlnBlock *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_link_AlnUnits_AlnBlock(AlnBlock * alb);
#define link_AlnUnits_AlnBlock Wise2_link_AlnUnits_AlnBlock


/* Function:  new_pairwise_AlnColumn(void)
 *
 * Descrip:    Function as a constructor for the special
 *             case of a pairwise alignment. Makes an
 *             AlnColumn and puts in two AlnUnits all ready
 *             to be linked in.
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [AlnColumn *]
 *
 */
AlnColumn * Wise2_new_pairwise_AlnColumn(void);
#define new_pairwise_AlnColumn Wise2_new_pairwise_AlnColumn


/* Function:  free_AlnColumn(obj)
 *
 * Descrip:    Specilased deconstructor needed because
 *             of linked list nature of the data structure
 *
 *
 * Arg:        obj [UNKN ] Undocumented argument [AlnColumn *]
 *
 * Return [UNKN ]  Undocumented return value [AlnColumn *]
 *
 */
AlnColumn * Wise2_free_AlnColumn(AlnColumn * obj);
#define free_AlnColumn Wise2_free_AlnColumn


/* Function:  free_AlnUnit(obj)
 *
 * Descrip:    Specilased deconstructor needed because
 *             of linked list nature of the data structure
 *
 *
 * Arg:        obj [UNKN ] Undocumented argument [AlnUnit *]
 *
 * Return [UNKN ]  Undocumented return value [AlnUnit *]
 *
 */
AlnUnit * Wise2_free_AlnUnit(AlnUnit * obj);
#define free_AlnUnit Wise2_free_AlnUnit


/* Function:  hard_link_AlnUnit(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [AlnUnit *]
 *
 * Return [UNKN ]  Undocumented return value [AlnUnit *]
 *
 */
AlnUnit * Wise2_hard_link_AlnUnit(AlnUnit * obj);
#define hard_link_AlnUnit Wise2_hard_link_AlnUnit


/* Function:  AlnUnit_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [AlnUnit *]
 *
 */
AlnUnit * Wise2_AlnUnit_alloc(void);
#define AlnUnit_alloc Wise2_AlnUnit_alloc


/* Function:  add_AlnColumn(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [AlnColumn *]
 * Arg:        add [OWNER] Object to add to the list [AlnUnit *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_AlnColumn(AlnColumn * obj,AlnUnit * add);
#define add_AlnColumn Wise2_add_AlnColumn


/* Function:  flush_AlnColumn(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [AlnColumn *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_AlnColumn(AlnColumn * obj);
#define flush_AlnColumn Wise2_flush_AlnColumn


/* Function:  AlnColumn_alloc_std(void)
 *
 * Descrip:    Equivalent to AlnColumn_alloc_len(AlnColumnLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [AlnColumn *]
 *
 */
AlnColumn * Wise2_AlnColumn_alloc_std(void);
#define AlnColumn_alloc_std Wise2_AlnColumn_alloc_std


/* Function:  AlnColumn_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [AlnColumn *]
 *
 */
AlnColumn * Wise2_AlnColumn_alloc_len(int len);
#define AlnColumn_alloc_len Wise2_AlnColumn_alloc_len


/* Function:  hard_link_AlnColumn(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [AlnColumn *]
 *
 * Return [UNKN ]  Undocumented return value [AlnColumn *]
 *
 */
AlnColumn * Wise2_hard_link_AlnColumn(AlnColumn * obj);
#define hard_link_AlnColumn Wise2_hard_link_AlnColumn


/* Function:  AlnColumn_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [AlnColumn *]
 *
 */
AlnColumn * Wise2_AlnColumn_alloc(void);
#define AlnColumn_alloc Wise2_AlnColumn_alloc


/* Function:  hard_link_AlnSequence(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [AlnSequence *]
 *
 * Return [UNKN ]  Undocumented return value [AlnSequence *]
 *
 */
AlnSequence * Wise2_hard_link_AlnSequence(AlnSequence * obj);
#define hard_link_AlnSequence Wise2_hard_link_AlnSequence


/* Function:  AlnSequence_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [AlnSequence *]
 *
 */
AlnSequence * Wise2_AlnSequence_alloc(void);
#define AlnSequence_alloc Wise2_AlnSequence_alloc


/* Function:  free_AlnSequence(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [AlnSequence *]
 *
 * Return [UNKN ]  Undocumented return value [AlnSequence *]
 *
 */
AlnSequence * Wise2_free_AlnSequence(AlnSequence * obj);
#define free_AlnSequence Wise2_free_AlnSequence


/* Function:  add_AlnBlock(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [AlnBlock *]
 * Arg:        add [OWNER] Object to add to the list [AlnSequence *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_AlnBlock(AlnBlock * obj,AlnSequence * add);
#define add_AlnBlock Wise2_add_AlnBlock


/* Function:  flush_AlnBlock(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [AlnBlock *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_AlnBlock(AlnBlock * obj);
#define flush_AlnBlock Wise2_flush_AlnBlock


/* Function:  AlnBlock_alloc_std(void)
 *
 * Descrip:    Equivalent to AlnBlock_alloc_len(AlnBlockLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [AlnBlock *]
 *
 */
AlnBlock * Wise2_AlnBlock_alloc_std(void);
#define AlnBlock_alloc_std Wise2_AlnBlock_alloc_std


/* Function:  AlnBlock_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [AlnBlock *]
 *
 */
AlnBlock * Wise2_AlnBlock_alloc_len(int len);
#define AlnBlock_alloc_len Wise2_AlnBlock_alloc_len


/* Function:  hard_link_AlnBlock(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [AlnBlock *]
 *
 * Return [UNKN ]  Undocumented return value [AlnBlock *]
 *
 */
AlnBlock * Wise2_hard_link_AlnBlock(AlnBlock * obj);
#define hard_link_AlnBlock Wise2_hard_link_AlnBlock


/* Function:  AlnBlock_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [AlnBlock *]
 *
 */
AlnBlock * Wise2_AlnBlock_alloc(void);
#define AlnBlock_alloc Wise2_AlnBlock_alloc


/* Function:  free_AlnBlock(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [AlnBlock *]
 *
 * Return [UNKN ]  Undocumented return value [AlnBlock *]
 *
 */
AlnBlock * Wise2_free_AlnBlock(AlnBlock * obj);
#define free_AlnBlock Wise2_free_AlnBlock


/* Function:  add_AlnBlockList(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [AlnBlockList *]
 * Arg:        add [OWNER] Object to add to the list [AlnBlock *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_AlnBlockList(AlnBlockList * obj,AlnBlock * add);
#define add_AlnBlockList Wise2_add_AlnBlockList


/* Function:  flush_AlnBlockList(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [AlnBlockList *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_AlnBlockList(AlnBlockList * obj);
#define flush_AlnBlockList Wise2_flush_AlnBlockList


/* Function:  AlnBlockList_alloc_std(void)
 *
 * Descrip:    Equivalent to AlnBlockList_alloc_len(AlnBlockListLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [AlnBlockList *]
 *
 */
AlnBlockList * Wise2_AlnBlockList_alloc_std(void);
#define AlnBlockList_alloc_std Wise2_AlnBlockList_alloc_std


/* Function:  AlnBlockList_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [AlnBlockList *]
 *
 */
AlnBlockList * Wise2_AlnBlockList_alloc_len(int len);
#define AlnBlockList_alloc_len Wise2_AlnBlockList_alloc_len


/* Function:  hard_link_AlnBlockList(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [AlnBlockList *]
 *
 * Return [UNKN ]  Undocumented return value [AlnBlockList *]
 *
 */
AlnBlockList * Wise2_hard_link_AlnBlockList(AlnBlockList * obj);
#define hard_link_AlnBlockList Wise2_hard_link_AlnBlockList


/* Function:  AlnBlockList_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [AlnBlockList *]
 *
 */
AlnBlockList * Wise2_AlnBlockList_alloc(void);
#define AlnBlockList_alloc Wise2_AlnBlockList_alloc


/* Function:  free_AlnBlockList(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [AlnBlockList *]
 *
 * Return [UNKN ]  Undocumented return value [AlnBlockList *]
 *
 */
AlnBlockList * Wise2_free_AlnBlockList(AlnBlockList * obj);
#define free_AlnBlockList Wise2_free_AlnBlockList


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
AlnColumn * Wise2_access_start_AlnBlock(AlnBlock * obj);
#define access_start_AlnBlock Wise2_access_start_AlnBlock
AlnUnit * Wise2_access_alu_AlnColumn(AlnColumn * obj,int i);
#define access_alu_AlnColumn Wise2_access_alu_AlnColumn
int Wise2_length_alu_AlnColumn(AlnColumn * obj);
#define length_alu_AlnColumn Wise2_length_alu_AlnColumn
int Wise2_access_label_AlnUnit(AlnUnit * obj);
#define access_label_AlnUnit Wise2_access_label_AlnUnit
boolean Wise2_replace_text_label_AlnUnit(AlnUnit * obj,char * text_label);
#define replace_text_label_AlnUnit Wise2_replace_text_label_AlnUnit
boolean Wise2_replace_start_AlnBlock(AlnBlock * obj,AlnColumn * start);
#define replace_start_AlnBlock Wise2_replace_start_AlnBlock
boolean Wise2_replace_bio_end_AlnSequence(AlnSequence * obj,int bio_end);
#define replace_bio_end_AlnSequence Wise2_replace_bio_end_AlnSequence
char * Wise2_access_text_label_AlnUnit(AlnUnit * obj);
#define access_text_label_AlnUnit Wise2_access_text_label_AlnUnit
int Wise2_length_seq_AlnBlock(AlnBlock * obj);
#define length_seq_AlnBlock Wise2_length_seq_AlnBlock
boolean Wise2_replace_next_AlnUnit(AlnUnit * obj,AlnUnit * next);
#define replace_next_AlnUnit Wise2_replace_next_AlnUnit
int Wise2_access_length_AlnBlock(AlnBlock * obj);
#define access_length_AlnBlock Wise2_access_length_AlnBlock
AlnUnit * Wise2_access_next_AlnUnit(AlnUnit * obj);
#define access_next_AlnUnit Wise2_access_next_AlnUnit
int Wise2_access_score_AlnBlock(AlnBlock * obj);
#define access_score_AlnBlock Wise2_access_score_AlnBlock
boolean Wise2_replace_in_column_AlnUnit(AlnUnit * obj,boolean in_column);
#define replace_in_column_AlnUnit Wise2_replace_in_column_AlnUnit
AlnColumn * Wise2_access_next_AlnColumn(AlnColumn * obj);
#define access_next_AlnColumn Wise2_access_next_AlnColumn
boolean Wise2_access_in_column_AlnUnit(AlnUnit * obj);
#define access_in_column_AlnUnit Wise2_access_in_column_AlnUnit
int Wise2_access_start_AlnUnit(AlnUnit * obj);
#define access_start_AlnUnit Wise2_access_start_AlnUnit
boolean Wise2_replace_seq_AlnUnit(AlnUnit * obj,AlnSequence * seq);
#define replace_seq_AlnUnit Wise2_replace_seq_AlnUnit
int Wise2_access_end_AlnUnit(AlnUnit * obj);
#define access_end_AlnUnit Wise2_access_end_AlnUnit
AlnSequence * Wise2_access_seq_AlnUnit(AlnUnit * obj);
#define access_seq_AlnUnit Wise2_access_seq_AlnUnit
int Wise2_access_bio_end_AlnSequence(AlnSequence * obj);
#define access_bio_end_AlnSequence Wise2_access_bio_end_AlnSequence
boolean Wise2_replace_start_AlnSequence(AlnSequence * obj,AlnUnit * start);
#define replace_start_AlnSequence Wise2_replace_start_AlnSequence
boolean Wise2_replace_length_AlnBlock(AlnBlock * obj,int length);
#define replace_length_AlnBlock Wise2_replace_length_AlnBlock
AlnUnit * Wise2_access_start_AlnSequence(AlnSequence * obj);
#define access_start_AlnSequence Wise2_access_start_AlnSequence
boolean Wise2_replace_next_AlnColumn(AlnColumn * obj,AlnColumn * next);
#define replace_next_AlnColumn Wise2_replace_next_AlnColumn
boolean Wise2_replace_data_type_AlnSequence(AlnSequence * obj,int data_type);
#define replace_data_type_AlnSequence Wise2_replace_data_type_AlnSequence
boolean Wise2_replace_end_AlnUnit(AlnUnit * obj,int end);
#define replace_end_AlnUnit Wise2_replace_end_AlnUnit
int Wise2_access_data_type_AlnSequence(AlnSequence * obj);
#define access_data_type_AlnSequence Wise2_access_data_type_AlnSequence
AlnSequence * Wise2_access_seq_AlnBlock(AlnBlock * obj,int i);
#define access_seq_AlnBlock Wise2_access_seq_AlnBlock
boolean Wise2_replace_data_AlnSequence(AlnSequence * obj,void * data);
#define replace_data_AlnSequence Wise2_replace_data_AlnSequence
boolean Wise2_replace_start_AlnUnit(AlnUnit * obj,int start);
#define replace_start_AlnUnit Wise2_replace_start_AlnUnit
void * Wise2_access_data_AlnSequence(AlnSequence * obj);
#define access_data_AlnSequence Wise2_access_data_AlnSequence
boolean Wise2_replace_score_AlnBlock(AlnBlock * obj,int score);
#define replace_score_AlnBlock Wise2_replace_score_AlnBlock
boolean Wise2_replace_bio_start_AlnSequence(AlnSequence * obj,int bio_start);
#define replace_bio_start_AlnSequence Wise2_replace_bio_start_AlnSequence
boolean Wise2_replace_label_AlnUnit(AlnUnit * obj,int label);
#define replace_label_AlnUnit Wise2_replace_label_AlnUnit
int Wise2_access_bio_start_AlnSequence(AlnSequence * obj);
#define access_bio_start_AlnSequence Wise2_access_bio_start_AlnSequence
boolean Wise2_can_swallow_AlnColumn(AlnColumn * master,AlnColumn * eaten,boolean (*comp_func)(char *,char *));
#define can_swallow_AlnColumn Wise2_can_swallow_AlnColumn
boolean Wise2_identical_labels_in_AlnColumn(AlnColumn * one,AlnColumn * two,boolean (*comp_func)(char *,char *));
#define identical_labels_in_AlnColumn Wise2_identical_labels_in_AlnColumn
boolean Wise2_identical_labels_in_AlnUnits(AlnUnit * one,AlnUnit * two,boolean (*comp_func)(char *,char *));
#define identical_labels_in_AlnUnits Wise2_identical_labels_in_AlnUnits
AlnColumn * Wise2_read_dumped_ascii_AlnColumn_line(char * line);
#define read_dumped_ascii_AlnColumn_line Wise2_read_dumped_ascii_AlnColumn_line
void Wise2_show_flat_AlnColumn(AlnColumn * alc,FILE * ofp);
#define show_flat_AlnColumn Wise2_show_flat_AlnColumn
void Wise2_show_flat_AlnUnit(AlnUnit * alu,FILE * ofp);
#define show_flat_AlnUnit Wise2_show_flat_AlnUnit
AlnUnit * Wise2_read_flat_AlnUnit_line(char * line,int * ret_pos);
#define read_flat_AlnUnit_line Wise2_read_flat_AlnUnit_line
void Wise2_swap_AlnColumn(AlnUnit ** list,int i,int j) ;
#define swap_AlnColumn Wise2_swap_AlnColumn
void Wise2_qsort_AlnColumn(AlnUnit ** list,int left,int right,int (*comp)(AlnUnit * ,AlnUnit * ));
#define qsort_AlnColumn Wise2_qsort_AlnColumn
void Wise2_sort_AlnColumn(AlnColumn * obj,int (*comp)(AlnUnit *, AlnUnit *));
#define sort_AlnColumn Wise2_sort_AlnColumn
boolean Wise2_expand_AlnColumn(AlnColumn * obj,int len);
#define expand_AlnColumn Wise2_expand_AlnColumn
void Wise2_swap_AlnBlock(AlnSequence ** list,int i,int j) ;
#define swap_AlnBlock Wise2_swap_AlnBlock
void Wise2_qsort_AlnBlock(AlnSequence ** list,int left,int right,int (*comp)(AlnSequence * ,AlnSequence * ));
#define qsort_AlnBlock Wise2_qsort_AlnBlock
void Wise2_sort_AlnBlock(AlnBlock * obj,int (*comp)(AlnSequence *, AlnSequence *));
#define sort_AlnBlock Wise2_sort_AlnBlock
boolean Wise2_expand_AlnBlock(AlnBlock * obj,int len);
#define expand_AlnBlock Wise2_expand_AlnBlock
void Wise2_swap_AlnBlockList(AlnBlock ** list,int i,int j) ;
#define swap_AlnBlockList Wise2_swap_AlnBlockList
void Wise2_qsort_AlnBlockList(AlnBlock ** list,int left,int right,int (*comp)(AlnBlock * ,AlnBlock * ));
#define qsort_AlnBlockList Wise2_qsort_AlnBlockList
void Wise2_sort_AlnBlockList(AlnBlockList * obj,int (*comp)(AlnBlock *, AlnBlock *));
#define sort_AlnBlockList Wise2_sort_AlnBlockList
boolean Wise2_expand_AlnBlockList(AlnBlockList * obj,int len);
#define expand_AlnBlockList Wise2_expand_AlnBlockList

#ifdef _cplusplus
}
#endif

#endif
