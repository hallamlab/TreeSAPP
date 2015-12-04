#ifndef DYNAMITEasciibtcanvasHEADERFILE
#define DYNAMITEasciibtcanvasHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "btcanvas.h"


struct Wise2_Ascii_btc_Data {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    FILE * ofp; /*  file to write to */ 
    int current_x;  /*  position in main line  */ 
    int paint_x;    /*  painting cursor on the line */ 
    int res_left;   /*  amount reserved on left */ 
    int main;   /*  main block amount */ 
    int res_right;  /*  amount reserved on right */ 
    char ** scratch;    /*  scratch pad lines  */ 
    int depth_scratch;  /*  depth of scratch pad for memory */ 
    boolean in_use;  
    btPasteArea * bpa;  /*  this is what we recycle.. */ 
    } ;  
/* Ascii_btc_Data defined */ 
#ifndef DYNAMITE_DEFINED_Ascii_btc_Data
typedef struct Wise2_Ascii_btc_Data Wise2_Ascii_btc_Data;
#define Ascii_btc_Data Wise2_Ascii_btc_Data
#define DYNAMITE_DEFINED_Ascii_btc_Data
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  ascii_btCanvas_usage(ofp)
 *
 * Descrip:      prints out usage for
 *
 *               ascii_btCanvas_from_commandline
 *
 *
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
void Wise2_ascii_btCanvas_usage(FILE * ofp);
#define ascii_btCanvas_usage Wise2_ascii_btCanvas_usage


/* Function:  ascii_btCanvas_from_commandline(argc,argv,default_left,default_main,default_right,ofp,height)
 *
 * Descrip:    Makes a ascii btCanvas from the command line,
 *             swallowing up options for
 *                -acleft
 *                -acright
 *                -acmain
 *
 *             If there are no options, builds with 15,50,5
 *
 *
 * Arg:                 argc [UNKN ] Undocumented argument [int *]
 * Arg:                 argv [UNKN ] Undocumented argument [char **]
 * Arg:         default_left [UNKN ] Undocumented argument [int]
 * Arg:         default_main [UNKN ] Undocumented argument [int]
 * Arg:        default_right [UNKN ] Undocumented argument [int]
 * Arg:                  ofp [UNKN ] Undocumented argument [FILE *]
 * Arg:               height [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [btCanvas *]
 *
 */
btCanvas * Wise2_ascii_btCanvas_from_commandline(int * argc,char ** argv,int default_left,int default_main,int default_right,FILE * ofp,int height);
#define ascii_btCanvas_from_commandline Wise2_ascii_btCanvas_from_commandline


/* Function:  new_Ascii_btCanvas(ofp,left,main,right,height)
 *
 * Descrip:    The only function specifically for Ascii bt Canvases.
 *
 *             Use this to make a new btCanvas. Then use functions like
 *             /get_paste_area_btCanvas to actually use it. Everything
 *             else is handled by data structures and pointer-to-functions
 *             which are hidden to you (and be thankful for that!)
 *
 *             The standard /free_btCanvas will free the hidden data structures
 *             as well
 *
 *
 * Arg:           ofp [UNKN ] FILE stream to write the ascii to [FILE *]
 * Arg:          left [UNKN ] amount of text to reserve on the left [int]
 * Arg:          main [UNKN ] Undocumented argument [int]
 * Arg:         right [UNKN ] amount of text to reserve on the right [int]
 * Arg:        height [UNKN ] height of block [int]
 *
 * Return [UNKN ]  Undocumented return value [btCanvas *]
 *
 */
btCanvas * Wise2_new_Ascii_btCanvas(FILE * ofp,int left,int main,int right,int height);
#define new_Ascii_btCanvas Wise2_new_Ascii_btCanvas


/* Function:  hard_link_Ascii_btc_Data(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [Ascii_btc_Data *]
 *
 * Return [UNKN ]  Undocumented return value [Ascii_btc_Data *]
 *
 */
Ascii_btc_Data * Wise2_hard_link_Ascii_btc_Data(Ascii_btc_Data * obj);
#define hard_link_Ascii_btc_Data Wise2_hard_link_Ascii_btc_Data


/* Function:  Ascii_btc_Data_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Ascii_btc_Data *]
 *
 */
Ascii_btc_Data * Wise2_Ascii_btc_Data_alloc(void);
#define Ascii_btc_Data_alloc Wise2_Ascii_btc_Data_alloc


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
btPasteArea * Wise2_next_Ascii_btpa(btCanvas * btc,int length);
#define next_Ascii_btpa Wise2_next_Ascii_btpa
boolean Wise2_paste_char_bt_Ascii(btPasteArea * bpa,int x,int y,char c,int format);
#define paste_char_bt_Ascii Wise2_paste_char_bt_Ascii
boolean Wise2_can_get_bt_Ascii(btCanvas * btc,int length);
#define can_get_bt_Ascii Wise2_can_get_bt_Ascii
btPasteArea * Wise2_get_ascii_left_btPasteArea(btCanvas * btc);
#define get_ascii_left_btPasteArea Wise2_get_ascii_left_btPasteArea
btPasteArea * Wise2_get_ascii_right_btPasteArea(btCanvas * btc);
#define get_ascii_right_btPasteArea Wise2_get_ascii_right_btPasteArea
boolean Wise2_ascii_next_line_btPasteArea(btCanvas * btc);
#define ascii_next_line_btPasteArea Wise2_ascii_next_line_btPasteArea
btPasteArea * Wise2_free_Ascii_btpa(btPasteArea * obj);
#define free_Ascii_btpa Wise2_free_Ascii_btpa
Ascii_btc_Data * Wise2_new_Ascii_btc_Data(FILE * ofp,int left,int main,int right,int height);
#define new_Ascii_btc_Data Wise2_new_Ascii_btc_Data
btCanvas * Wise2_free_Ascii_btc(btCanvas * btc);
#define free_Ascii_btc Wise2_free_Ascii_btc
Ascii_btc_Data * Wise2_free_Ascii_btc_Data(Ascii_btc_Data * obj);
#define free_Ascii_btc_Data Wise2_free_Ascii_btc_Data

#ifdef _cplusplus
}
#endif

#endif
