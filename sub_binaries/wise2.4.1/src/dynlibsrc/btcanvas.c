#ifdef _cplusplus
extern "C" {
#endif
#include "btcanvas.h"

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
# line 70 "btcanvas.dy"
boolean can_get_paste_area_btCanvas(btCanvas * btc,int length)
{

  return ((*btc->can_get_paste_area)(btc,length));
}


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
# line 81 "btcanvas.dy"
boolean advance_line_btCanvas(btCanvas * btc)
{

  return ((*btc->advance_line)(btc));
}

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
# line 92 "btcanvas.dy"
btPasteArea * get_reserved_left_btCanvas(btCanvas * btc)
{
  return ((*btc->get_reserved_left)(btc));
}

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
# line 102 "btcanvas.dy"
btPasteArea * get_reserved_right_btCanvas(btCanvas * btc)
{
  return ((*btc->get_reserved_right)(btc));
}


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
# line 113 "btcanvas.dy"
btPasteArea * get_paste_area_btCanvas(btCanvas * btc,int length)
{
  return ((*btc->get_paste_area)(btc,length));
}
  
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
# line 138 "btcanvas.dy"
boolean paste_substr_btPasteArea(btPasteArea * bta,int x,int y,const char * str,int len,char (*map_func)(char),btCanvasDirection dir,int format)
{
  char buf[64];
  const char *run;
  int i;

  for(run=str,i=0;i<len;i++,run++) {
    if( map_func == NULL ) {
      buf[i] = *run;
    } else {
      buf[i] = (*map_func)(*run);
    }
  }

  buf[i] = '\0';

  return paste_string_btPasteArea(bta,x,y,buf,dir,format);
}


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
# line 169 "btcanvas.dy"
boolean paste_string_btPasteArea(btPasteArea * bta,int x,int y,const char * str,btCanvasDirection dir,int format)
{
  const char * run;

  for(run=str;*run;run++) {
    if( paste_char_btPasteArea(bta,x,y,*run,format) == FALSE ) {
      warn("Unable to paste the word %s into btcanvas",str);
      return FALSE;
    }
    switch (dir) {
    case BC_UP : y--; break;
    case BC_DOWN : y++; break;
    case BC_RIGHT : x++; break;
    case BC_LEFT : x--; break;
    default :
      warn("You have not put in a valid direction into the paste_string function. Sod off... ");
      return FALSE;
    }
  }

  return TRUE;
}
    
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
# line 198 "btcanvas.dy"
boolean paste_char_btPasteArea(btPasteArea * bta,int x,int y,char c,int format) 
{
  if( x > bta->length || y > bta->height ) {
    warn("Trying to paste a character into an unpasteable position [%d,%d] into [%d,%d]",x,y,bta->length,bta->height);
    return FALSE;
  }

  return ((*bta->paste_char))(bta,x,y,c,format);

}

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
# line 216 "btcanvas.dy"
btPasteArea * free_btPasteArea(btPasteArea * obj)
{
  return ((*obj->decons))(obj);
}

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
# line 228 "btcanvas.dy"
btCanvas * free_btCanvas(btCanvas * obj)
{
  return ((*obj->decons))(obj);
}

# line 223 "btcanvas.c"
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
btPasteArea * hard_link_btPasteArea(btPasteArea * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a btPasteArea object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  btPasteArea_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [btPasteArea *]
 *
 */
btPasteArea * btPasteArea_alloc(void) 
{
    btPasteArea * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(btPasteArea *) ckalloc (sizeof(btPasteArea))) == NULL)  {  
      warn("btPasteArea_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->type = 0;   
    out->height = 0; 
    out->length = 0; 
    out->paste_char = NULL;  
    out->canvas_data = NULL; 
    out->decons = NULL;  


    return out;  
}    


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
btCanvas * hard_link_btCanvas(btCanvas * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a btCanvas object: passed a NULL object");    
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  btCanvas_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [btCanvas *]
 *
 */
btCanvas * btCanvas_alloc(void) 
{
    btCanvas * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(btCanvas *) ckalloc (sizeof(btCanvas))) == NULL)    {  
      warn("btCanvas_alloc failed ");    
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->type = 0;   
    out->height = 0; 
    out->res_right = 0;  
    out->res_left = 0;   
    out->can_get_paste_area = NULL;  
    out->advance_line = NULL;    
    out->get_paste_area = NULL;  
    out->get_reserved_right = NULL;  
    out->get_reserved_left = NULL;   
    out->canvas_data = NULL; 
    out->decons = NULL;  


    return out;  
}    



#ifdef _cplusplus
}
#endif
