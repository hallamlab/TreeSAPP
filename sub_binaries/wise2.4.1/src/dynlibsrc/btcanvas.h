#ifndef DYNAMITEbtcanvasHEADERFILE
#define DYNAMITEbtcanvasHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "wisebase.h"


typedef enum btCanvasDirection {
  BC_UP,
  BC_DOWN,
  BC_RIGHT,
  BC_LEFT 
} btCanvasDirection ;

/*
 * Icky... icky.. icky
 *
 * Would you believe, some C compilers will not allow typdef's
 * inside function to pointer definitions in structures...
 *
 * Hence - hard coded package names!!!! Yuk!!!
 *
 */
#ifndef DYNAMITE_DEFINED_btPastArea
typedef struct Wise2_btPastArea Wise2_btPastArea;
#define btPastArea Wise2_btPastArea
#define DYNAMITE_DEFINED_btPastArea
#endif

#ifndef DYNAMITE_DEFINED_btCanvas
typedef struct Wise2_btCanvas Wise2_btCanvas;
#define btCanvas Wise2_btCanvas
#define DYNAMITE_DEFINED_btCanvas
#endif

struct Wise2_btPasteArea {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int type;    
    int height;  
    int length;  
    boolean (*paste_char)(struct Wise2_btPasteArea *,int,int,char,int);  
    void * canvas_data;  
    struct Wise2_btPasteArea * (*decons)(struct Wise2_btPasteArea * bt); 
    } ;  
/* btPasteArea defined */ 
#ifndef DYNAMITE_DEFINED_btPasteArea
typedef struct Wise2_btPasteArea Wise2_btPasteArea;
#define btPasteArea Wise2_btPasteArea
#define DYNAMITE_DEFINED_btPasteArea
#endif


/* Object btCanvas
 *
 * Descrip: No Description
 *
 */
struct Wise2_btCanvas {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int type;    
    int height;  
    int res_right;   
    int res_left;    
    boolean (*can_get_paste_area)(struct Wise2_btCanvas *,int);  
    boolean (*advance_line)(struct Wise2_btCanvas *);    
    struct Wise2_btPasteArea * (*get_paste_area)(struct Wise2_btCanvas *,int);   
    struct Wise2_btPasteArea * (*get_reserved_right)(struct Wise2_btCanvas *);   
    struct Wise2_btPasteArea * (*get_reserved_left)(struct Wise2_btCanvas *);    
    void * canvas_data;  
    struct Wise2_btCanvas * (*decons)(struct Wise2_btCanvas * btc);  
    } ;  
/* btCanvas defined */ 
#ifndef DYNAMITE_DEFINED_btCanvas
typedef struct Wise2_btCanvas Wise2_btCanvas;
#define btCanvas Wise2_btCanvas
#define DYNAMITE_DEFINED_btCanvas
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  can_get_paste_area_btCanvas(btc,length)
 *
 * Descrip:    This returns whether there is room for length
 *             of chars in this line. If not, advance a line
 *
 *
 * Arg:           btc [UNKN ] Undocumented argument [btCanvas *]
 * Arg:        length [UNKN ] length of block wanted [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_can_get_paste_area_btCanvas(btCanvas * btc,int length);
#define can_get_paste_area_btCanvas Wise2_can_get_paste_area_btCanvas


/* Function:  advance_line_btCanvas(btc)
 *
 * Descrip:    Advances the canvas to the next line
 *
 *
 * Arg:        btc [UNKN ] Undocumented argument [btCanvas *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_advance_line_btCanvas(btCanvas * btc);
#define advance_line_btCanvas Wise2_advance_line_btCanvas


/* Function:  get_reserved_left_btCanvas(btc)
 *
 * Descrip:    This gets a paste-able area at the left hand
 *             side of block
 *
 *
 * Arg:        btc [UNKN ] Undocumented argument [btCanvas *]
 *
 * Return [UNKN ]  Undocumented return value [btPasteArea *]
 *
 */
btPasteArea * Wise2_get_reserved_left_btCanvas(btCanvas * btc);
#define get_reserved_left_btCanvas Wise2_get_reserved_left_btCanvas


/* Function:  get_reserved_right_btCanvas(btc)
 *
 * Descrip:    This gets a paste-able area at the right hand
 *             side of block
 *
 *
 * Arg:        btc [UNKN ] Undocumented argument [btCanvas *]
 *
 * Return [UNKN ]  Undocumented return value [btPasteArea *]
 *
 */
btPasteArea * Wise2_get_reserved_right_btCanvas(btCanvas * btc);
#define get_reserved_right_btCanvas Wise2_get_reserved_right_btCanvas


/* Function:  get_paste_area_btCanvas(btc,length)
 *
 * Descrip:    This gets a paste-able area of a certain length from
 *             the current cursor. NULL on error
 *
 *
 * Arg:           btc [UNKN ] Undocumented argument [btCanvas *]
 * Arg:        length [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [btPasteArea *]
 *
 */
btPasteArea * Wise2_get_paste_area_btCanvas(btCanvas * btc,int length);
#define get_paste_area_btCanvas Wise2_get_paste_area_btCanvas


/* Function:  paste_substr_btPasteArea(bta,x,y,str,len,map_func,dir,format)
 *
 * Descrip:    This will paste the substr of length len (obviously this could
 *             not be a \0 terminated string) starting from str[0] into btCanvas
 *             If map_func is non NULL, it converts each of the characters 
 *             using the map_func first. This is surprisingly useful.
 *
 *             map_func can be NULL, in which case the characters are just
 *             used as seen
 *
 *             Length has to be less than 64 chars
 *
 *
 * Arg:             bta [UNKN ] PasteArea to be used [btPasteArea *]
 * Arg:               x [UNKN ] x position of start of string [int]
 * Arg:               y [UNKN ] y position of start of string [int]
 * Arg:             str [READ ] string to get the substr from to be pasted [const char *]
 * Arg:             len [READ ] length of the substr [int]
 * Arg:        map_func [FUNCP] Function to specifiy the mapping of chars from str to display [char (*map_func]
 * Arg:             dir [UNKN ] direction (BC_UP, BC_DOWN, BC_RIGHT, BC_LEFT) [btCanvasDirection]
 * Arg:          format [UNKN ] format (if recognised by the cnavas) [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_paste_substr_btPasteArea(btPasteArea * bta,int x,int y,const char * str,int len,char (*map_func)(char),btCanvasDirection dir,int format);
#define paste_substr_btPasteArea Wise2_paste_substr_btPasteArea


/* Function:  paste_string_btPasteArea(bta,x,y,str,dir,format)
 *
 * Descrip:    This will paste the '\0' terminated string into x,y in the direction 
 *             specified (up, down, left or right)
 *
 *
 * Arg:           bta [UNKN ] PasteArea to be used [btPasteArea *]
 * Arg:             x [UNKN ] x position of start of string [int]
 * Arg:             y [UNKN ] y position of start of string [int]
 * Arg:           str [READ ] string to be pasted [const char *]
 * Arg:           dir [UNKN ] direction (BC_UP, BC_DOWN, BC_RIGHT, BC_LEFT) [btCanvasDirection]
 * Arg:        format [UNKN ] format (if recognised by the cnavas) [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_paste_string_btPasteArea(btPasteArea * bta,int x,int y,const char * str,btCanvasDirection dir,int format);
#define paste_string_btPasteArea Wise2_paste_string_btPasteArea


/* Function:  paste_char_btPasteArea(bta,x,y,c,format)
 *
 * Descrip:    This will paste one character at x,y into the paste area.
 *             If the canvas understands the format style, it will take
 *             notice of it.
 *
 *
 * Arg:           bta [UNKN ] Undocumented argument [btPasteArea *]
 * Arg:             x [UNKN ] Undocumented argument [int]
 * Arg:             y [UNKN ] Undocumented argument [int]
 * Arg:             c [UNKN ] Undocumented argument [char]
 * Arg:        format [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_paste_char_btPasteArea(btPasteArea * bta,int x,int y,char c,int format) ;
#define paste_char_btPasteArea Wise2_paste_char_btPasteArea


/* Function:  free_btPasteArea(obj)
 *
 * Descrip:    Specialised deconstructor. Ensures the 
 *             data structures are freed
 *
 *
 * Arg:        obj [WRITE] btPasteArea obj to be destroyed [btPasteArea *]
 *
 * Return [UNKN ]  Undocumented return value [btPasteArea *]
 *
 */
btPasteArea * Wise2_free_btPasteArea(btPasteArea * obj);
#define free_btPasteArea Wise2_free_btPasteArea


/* Function:  free_btCanvas(obj)
 *
 * Descrip:    Specialised deconstructor. Ensures the 
 *             data structures are freed
 *
 *
 * Arg:        obj [WRITE] btCanvas obj to be destroyed [btCanvas *]
 *
 * Return [UNKN ]  Undocumented return value [btCanvas *]
 *
 */
btCanvas * Wise2_free_btCanvas(btCanvas * obj);
#define free_btCanvas Wise2_free_btCanvas


/* Function:  hard_link_btPasteArea(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [btPasteArea *]
 *
 * Return [UNKN ]  Undocumented return value [btPasteArea *]
 *
 */
btPasteArea * Wise2_hard_link_btPasteArea(btPasteArea * obj);
#define hard_link_btPasteArea Wise2_hard_link_btPasteArea


/* Function:  btPasteArea_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [btPasteArea *]
 *
 */
btPasteArea * Wise2_btPasteArea_alloc(void);
#define btPasteArea_alloc Wise2_btPasteArea_alloc


/* Function:  hard_link_btCanvas(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [btCanvas *]
 *
 * Return [UNKN ]  Undocumented return value [btCanvas *]
 *
 */
btCanvas * Wise2_hard_link_btCanvas(btCanvas * obj);
#define hard_link_btCanvas Wise2_hard_link_btCanvas


/* Function:  btCanvas_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [btCanvas *]
 *
 */
btCanvas * Wise2_btCanvas_alloc(void);
#define btCanvas_alloc Wise2_btCanvas_alloc


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/

#ifdef _cplusplus
}
#endif

#endif
