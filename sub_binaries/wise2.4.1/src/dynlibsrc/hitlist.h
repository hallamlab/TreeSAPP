#ifndef DYNAMITEhitlistHEADERFILE
#define DYNAMITEhitlistHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "sequence.h"
#include "hsp.h"
#include "aln.h"
#include "compmat.h"
#include "btcanvas.h"
#include "asciibtcanvas.h"
#include "searchstatinterface.h"


#define HitListLISTLENGTH 256
#define HitPairLISTLENGTH 16

typedef enum HitListOutputFormat {
  HitListOutputFormatPseudoBlast = 34,
  HitListOutputFormatXML,
  HitListOutputFormatTab,
  HitListAlnCumlative,
  HitListOutputFormatUnknown
} HitListOutputFormat;

struct Wise2_HitAln {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int raw_score;   
    double bit_score;    
    double evalue;   
    AlnBlock * alb;  
    } ;  
/* HitAln defined */ 
#ifndef DYNAMITE_DEFINED_HitAln
typedef struct Wise2_HitAln Wise2_HitAln;
#define HitAln Wise2_HitAln
#define DYNAMITE_DEFINED_HitAln
#endif


struct Wise2_HitPair {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    Sequence * query;    
    Sequence * target;   
    int raw_score;   
    double bit_score;    
    double evalue;   
    HitAln ** aln;   
    int len;/* len for above aln  */ 
    int maxlen; /* maxlen for above aln */ 
    boolean target_reversed;     
    } ;  
/* HitPair defined */ 
#ifndef DYNAMITE_DEFINED_HitPair
typedef struct Wise2_HitPair Wise2_HitPair;
#define HitPair Wise2_HitPair
#define DYNAMITE_DEFINED_HitPair
#endif


struct Wise2_HitList {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    HitPair ** pair;     
    int len;/* len for above pair  */ 
    int maxlen; /* maxlen for above pair */ 
    CompMat * mat;   
    boolean (*write_btc_func)(AlnBlock *,Sequence *,Sequence *,btCanvas * btc);  
    char * stat_attrib;  
    } ;  
/* HitList defined */ 
#ifndef DYNAMITE_DEFINED_HitList
typedef struct Wise2_HitList Wise2_HitList;
#define HitList Wise2_HitList
#define DYNAMITE_DEFINED_HitList
#endif


struct Wise2_HitListOutputImpl {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    HitListOutputFormat type;    
    } ;  
/* HitListOutputImpl defined */ 
#ifndef DYNAMITE_DEFINED_HitListOutputImpl
typedef struct Wise2_HitListOutputImpl Wise2_HitListOutputImpl;
#define HitListOutputImpl Wise2_HitListOutputImpl
#define DYNAMITE_DEFINED_HitListOutputImpl
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  sort_HitList_by_score(hl)
 *
 * Descrip:    Sorts by score
 *
 *
 * Arg:        hl [UNKN ] Undocumented argument [HitList *]
 *
 */
void Wise2_sort_HitList_by_score(HitList * hl);
#define sort_HitList_by_score Wise2_sort_HitList_by_score


/* Function:  compare_HitPair_score(one,two)
 *
 * Descrip:    internal function to sort by score
 *
 *
 * Arg:        one [UNKN ] Undocumented argument [HitPair *]
 * Arg:        two [UNKN ] Undocumented argument [HitPair *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_compare_HitPair_score(HitPair * one,HitPair * two) ;
#define compare_HitPair_score Wise2_compare_HitPair_score


/* Function:  apply_SearchStat_to_HitList(hspm,ssi,database_size)
 *
 * Descrip:    Applies statistics across a hitlist
 *
 *
 * Arg:                 hspm [UNKN ] Undocumented argument [HitList *]
 * Arg:                  ssi [UNKN ] Undocumented argument [SearchStatInterface *]
 * Arg:        database_size [UNKN ] Undocumented argument [int]
 *
 */
void Wise2_apply_SearchStat_to_HitList(HitList * hspm,SearchStatInterface * ssi,int database_size);
#define apply_SearchStat_to_HitList Wise2_apply_SearchStat_to_HitList


/* Function:  HitList_from_LinearHSPmanager(lm)
 *
 * Descrip:    Converts a LinearHSPmanager into a HitList
 *
 *
 * Arg:        lm [UNKN ] Undocumented argument [LinearHSPmanager *]
 *
 * Return [UNKN ]  Undocumented return value [HitList *]
 *
 */
HitList * Wise2_HitList_from_LinearHSPmanager(LinearHSPmanager * lm);
#define HitList_from_LinearHSPmanager Wise2_HitList_from_LinearHSPmanager


/* Function:  HitPair_from_HSPset(set,mat)
 *
 * Descrip:    Builds a Hitpair from an HSP, not doing 
 *             alignment
 *
 *
 * Arg:        set [UNKN ] Undocumented argument [HSPset *]
 * Arg:        mat [UNKN ] Undocumented argument [CompMat *]
 *
 * Return [UNKN ]  Undocumented return value [HitPair *]
 *
 */
HitPair * Wise2_HitPair_from_HSPset(HSPset * set,CompMat * mat);
#define HitPair_from_HSPset Wise2_HitPair_from_HSPset


/* Function:  ungapped_AlnBlock_from_HSP(hsp,q,t,mat)
 *
 * Descrip:    Builds an expanded AlnBlock with one AlnColumn 
 *             per residue for an ungapped HSP
 *
 *
 * Arg:        hsp [UNKN ] Undocumented argument [HSP *]
 * Arg:          q [UNKN ] Undocumented argument [Sequence *]
 * Arg:          t [UNKN ] Undocumented argument [Sequence *]
 * Arg:        mat [UNKN ] Undocumented argument [CompMat *]
 *
 * Return [UNKN ]  Undocumented return value [AlnBlock *]
 *
 */
AlnBlock * Wise2_ungapped_AlnBlock_from_HSP(HSP * hsp,Sequence * q,Sequence * t,CompMat * mat);
#define ungapped_AlnBlock_from_HSP Wise2_ungapped_AlnBlock_from_HSP


/* Function:  new_HitListOutputImpl_from_argv(argc,argv)
 *
 * Descrip:    Builds a new HitListOutputFormat from commandline
 *
 *
 * Arg:        argc [UNKN ] Undocumented argument [int *]
 * Arg:        argv [UNKN ] Undocumented argument [char **]
 *
 * Return [UNKN ]  Undocumented return value [HitListOutputImpl *]
 *
 */
HitListOutputImpl * Wise2_new_HitListOutputImpl_from_argv(int * argc,char ** argv);
#define new_HitListOutputImpl_from_argv Wise2_new_HitListOutputImpl_from_argv


/* Function:  show_help_HitListOutputImpl(ofp)
 *
 * Descrip:    Shows help for HitList output
 *
 *
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
void Wise2_show_help_HitListOutputImpl(FILE * ofp);
#define show_help_HitListOutputImpl Wise2_show_help_HitListOutputImpl


/* Function:  show_HitList_HitListOutputImpl(hloi,hl,ofp)
 *
 * Descrip:    Shows a hitlist wrt to output impl
 *
 *
 * Arg:        hloi [UNKN ] Undocumented argument [HitListOutputImpl *]
 * Arg:          hl [UNKN ] Undocumented argument [HitList *]
 * Arg:         ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
void Wise2_show_HitList_HitListOutputImpl(HitListOutputImpl * hloi,HitList * hl,FILE * ofp);
#define show_HitList_HitListOutputImpl Wise2_show_HitList_HitListOutputImpl


/* Function:  write_alb_HitList(hl,ofp)
 *
 * Descrip:    Writes Alb output
 *
 *
 * Arg:         hl [UNKN ] Undocumented argument [HitList *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
void Wise2_write_alb_HitList(HitList * hl,FILE * ofp);
#define write_alb_HitList Wise2_write_alb_HitList


/* Function:  write_XML_HitList(hl,ofp)
 *
 * Descrip:    Writes XML output
 *
 *
 * Arg:         hl [UNKN ] Undocumented argument [HitList *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
void Wise2_write_XML_HitList(HitList * hl,FILE * ofp);
#define write_XML_HitList Wise2_write_XML_HitList


/* Function:  write_tab_HitList(hl,ofp)
 *
 * Descrip:    Writes tab delimited tab like output
 *
 *
 * Arg:         hl [UNKN ] Undocumented argument [HitList *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
void Wise2_write_tab_HitList(HitList * hl,FILE * ofp);
#define write_tab_HitList Wise2_write_tab_HitList


/* Function:  write_pseudoblast_HitList(hl,ofp)
 *
 * Descrip:    Writes pseudoblast output
 *
 *
 * Arg:         hl [UNKN ] Undocumented argument [HitList *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
void Wise2_write_pseudoblast_HitList(HitList * hl,FILE * ofp);
#define write_pseudoblast_HitList Wise2_write_pseudoblast_HitList


/* Function:  write_pretty_Seq_blast_align_btc(alb,one,two,btc)
 *
 * Descrip:    Chains up to char* level alignment writer
 *
 *
 * Arg:        alb [UNKN ] Undocumented argument [AlnBlock *]
 * Arg:        one [UNKN ] Undocumented argument [Sequence *]
 * Arg:        two [UNKN ] Undocumented argument [Sequence *]
 * Arg:        btc [UNKN ] Undocumented argument [btCanvas *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_write_pretty_Seq_blast_align_btc(AlnBlock * alb,Sequence * one,Sequence * two,btCanvas * btc);
#define write_pretty_Seq_blast_align_btc Wise2_write_pretty_Seq_blast_align_btc


/* Function:  write_pretty_str_blast_align_btc(alb,qname,query,tname,target,btc)
 *
 * Descrip:    This function writes precisely
 *             what you expect for a a simple alignment.
 *
 *             We can reuse this routine all over the place because 
 *             we dont use any hard coded structure for the
 *             query or the target sequence letters. ... but crap
 *             type checking it has to be said!
 *
 *             Also we use a generic btCanvas that could have
 *             any implementation underneath (eg, ASCII, postscript etc).
 *
 *
 * Arg:           alb [UNKN ] Undocumented argument [AlnBlock *]
 * Arg:         qname [UNKN ] Undocumented argument [char *]
 * Arg:         query [UNKN ] Undocumented argument [char *]
 * Arg:         tname [UNKN ] Undocumented argument [char *]
 * Arg:        target [UNKN ] Undocumented argument [char *]
 * Arg:           btc [UNKN ] Undocumented argument [btCanvas *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_write_pretty_str_blast_align_btc(AlnBlock * alb,char * qname,char * query,char * tname,char * target,btCanvas * btc);
#define write_pretty_str_blast_align_btc Wise2_write_pretty_str_blast_align_btc


/* Function:  hard_link_HitAln(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [HitAln *]
 *
 * Return [UNKN ]  Undocumented return value [HitAln *]
 *
 */
HitAln * Wise2_hard_link_HitAln(HitAln * obj);
#define hard_link_HitAln Wise2_hard_link_HitAln


/* Function:  HitAln_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [HitAln *]
 *
 */
HitAln * Wise2_HitAln_alloc(void);
#define HitAln_alloc Wise2_HitAln_alloc


/* Function:  free_HitAln(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [HitAln *]
 *
 * Return [UNKN ]  Undocumented return value [HitAln *]
 *
 */
HitAln * Wise2_free_HitAln(HitAln * obj);
#define free_HitAln Wise2_free_HitAln


/* Function:  add_HitPair(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [HitPair *]
 * Arg:        add [OWNER] Object to add to the list [HitAln *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_HitPair(HitPair * obj,HitAln * add);
#define add_HitPair Wise2_add_HitPair


/* Function:  flush_HitPair(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [HitPair *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_HitPair(HitPair * obj);
#define flush_HitPair Wise2_flush_HitPair


/* Function:  HitPair_alloc_std(void)
 *
 * Descrip:    Equivalent to HitPair_alloc_len(HitPairLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [HitPair *]
 *
 */
HitPair * Wise2_HitPair_alloc_std(void);
#define HitPair_alloc_std Wise2_HitPair_alloc_std


/* Function:  HitPair_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [HitPair *]
 *
 */
HitPair * Wise2_HitPair_alloc_len(int len);
#define HitPair_alloc_len Wise2_HitPair_alloc_len


/* Function:  hard_link_HitPair(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [HitPair *]
 *
 * Return [UNKN ]  Undocumented return value [HitPair *]
 *
 */
HitPair * Wise2_hard_link_HitPair(HitPair * obj);
#define hard_link_HitPair Wise2_hard_link_HitPair


/* Function:  HitPair_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [HitPair *]
 *
 */
HitPair * Wise2_HitPair_alloc(void);
#define HitPair_alloc Wise2_HitPair_alloc


/* Function:  free_HitPair(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [HitPair *]
 *
 * Return [UNKN ]  Undocumented return value [HitPair *]
 *
 */
HitPair * Wise2_free_HitPair(HitPair * obj);
#define free_HitPair Wise2_free_HitPair


/* Function:  add_HitList(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [HitList *]
 * Arg:        add [OWNER] Object to add to the list [HitPair *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_HitList(HitList * obj,HitPair * add);
#define add_HitList Wise2_add_HitList


/* Function:  flush_HitList(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [HitList *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_HitList(HitList * obj);
#define flush_HitList Wise2_flush_HitList


/* Function:  HitList_alloc_std(void)
 *
 * Descrip:    Equivalent to HitList_alloc_len(HitListLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [HitList *]
 *
 */
HitList * Wise2_HitList_alloc_std(void);
#define HitList_alloc_std Wise2_HitList_alloc_std


/* Function:  HitList_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [HitList *]
 *
 */
HitList * Wise2_HitList_alloc_len(int len);
#define HitList_alloc_len Wise2_HitList_alloc_len


/* Function:  hard_link_HitList(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [HitList *]
 *
 * Return [UNKN ]  Undocumented return value [HitList *]
 *
 */
HitList * Wise2_hard_link_HitList(HitList * obj);
#define hard_link_HitList Wise2_hard_link_HitList


/* Function:  HitList_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [HitList *]
 *
 */
HitList * Wise2_HitList_alloc(void);
#define HitList_alloc Wise2_HitList_alloc


/* Function:  free_HitList(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [HitList *]
 *
 * Return [UNKN ]  Undocumented return value [HitList *]
 *
 */
HitList * Wise2_free_HitList(HitList * obj);
#define free_HitList Wise2_free_HitList


/* Function:  hard_link_HitListOutputImpl(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [HitListOutputImpl *]
 *
 * Return [UNKN ]  Undocumented return value [HitListOutputImpl *]
 *
 */
HitListOutputImpl * Wise2_hard_link_HitListOutputImpl(HitListOutputImpl * obj);
#define hard_link_HitListOutputImpl Wise2_hard_link_HitListOutputImpl


/* Function:  HitListOutputImpl_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [HitListOutputImpl *]
 *
 */
HitListOutputImpl * Wise2_HitListOutputImpl_alloc(void);
#define HitListOutputImpl_alloc Wise2_HitListOutputImpl_alloc


/* Function:  free_HitListOutputImpl(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [HitListOutputImpl *]
 *
 * Return [UNKN ]  Undocumented return value [HitListOutputImpl *]
 *
 */
HitListOutputImpl * Wise2_free_HitListOutputImpl(HitListOutputImpl * obj);
#define free_HitListOutputImpl Wise2_free_HitListOutputImpl


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
void Wise2_swap_HitPair(HitAln ** list,int i,int j) ;
#define swap_HitPair Wise2_swap_HitPair
void Wise2_qsort_HitPair(HitAln ** list,int left,int right,int (*comp)(HitAln * ,HitAln * ));
#define qsort_HitPair Wise2_qsort_HitPair
void Wise2_sort_HitPair(HitPair * obj,int (*comp)(HitAln *, HitAln *));
#define sort_HitPair Wise2_sort_HitPair
boolean Wise2_expand_HitPair(HitPair * obj,int len);
#define expand_HitPair Wise2_expand_HitPair
void Wise2_swap_HitList(HitPair ** list,int i,int j) ;
#define swap_HitList Wise2_swap_HitList
void Wise2_qsort_HitList(HitPair ** list,int left,int right,int (*comp)(HitPair * ,HitPair * ));
#define qsort_HitList Wise2_qsort_HitList
void Wise2_sort_HitList(HitList * obj,int (*comp)(HitPair *, HitPair *));
#define sort_HitList Wise2_sort_HitList
boolean Wise2_expand_HitList(HitList * obj,int len);
#define expand_HitList Wise2_expand_HitList

#ifdef _cplusplus
}
#endif

#endif
