#ifdef _cplusplus
extern "C" {
#endif
#include "asciibtcanvas.h"

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
# line 34 "asciibtcanvas.dy"
void ascii_btCanvas_usage(FILE * ofp)
{
  fprintf(ofp,"Options for ascii canvas\n");
  fprintf(ofp,"  -acleft  [integer]    size of left hand panel\n");
  fprintf(ofp,"  -acmain  [integer]    size of main panel\n");
  fprintf(ofp,"  -acright [integer]    size of right hand panel\n");
}

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
# line 51 "asciibtcanvas.dy"
btCanvas * ascii_btCanvas_from_commandline(int * argc,char ** argv,int default_left,int default_main,int default_right,FILE * ofp,int height)
{
  int left = default_left;
  int main_length = default_main;
  int right = default_right;

  strip_out_integer_argument(argc,argv,"acleft",&left);
  strip_out_integer_argument(argc,argv,"acright",&right);
  strip_out_integer_argument(argc,argv,"acmain",&main_length);

  return new_Ascii_btCanvas(ofp,left,main_length,right,height);
}


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
# line 81 "asciibtcanvas.dy"
btCanvas * new_Ascii_btCanvas(FILE * ofp,int left,int main,int right,int height)
{
  btCanvas * out;
  Ascii_btc_Data * d;

  out = btCanvas_alloc();

  d = new_Ascii_btc_Data(ofp,left,main,right,height);

  out->canvas_data = (void *) d;

  out->decons = free_Ascii_btc;

  /* now need to put in correct pointers to functions for
   * canvas implementation
   */
  out->can_get_paste_area = can_get_bt_Ascii;
  out->advance_line = ascii_next_line_btPasteArea;
  out->get_paste_area = next_Ascii_btpa;
  out->get_reserved_right = get_ascii_right_btPasteArea;
  out->get_reserved_left = get_ascii_left_btPasteArea;


  return out;

}



  
/* Function:  next_Ascii_btpa(btc,length)
 *
 * Descrip:    gets the next btPasteArea. Here we con
 *             people by simply passing out the 'bpa'
 *             held in canvas and never reallocating it
 *             (see /free_Ascii_btpa which is already attached to
 *             the bpa).
 *
 *
 * Arg:           btc [UNKN ] Undocumented argument [btCanvas *]
 * Arg:        length [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [btPasteArea *]
 *
 */
# line 119 "asciibtcanvas.dy"
btPasteArea * next_Ascii_btpa(btCanvas * btc,int length)
{
  Ascii_btc_Data * abd = NULL;
  
  abd = (Ascii_btc_Data*) btc->canvas_data;

  if( abd->in_use == TRUE ) { 
    warn("You are already using a btPasteArea on this canvas. Only one at a time! Probably you have not freed the btPasteArea before hand");
    return NULL;
  }
  
  if( abd->current_x + length > abd->main + abd->res_left ) {
    warn("Asking for more block than I can give you. You have not tested with can_get_paste_area. Bad boy!");
    return NULL;
  }

  abd->in_use = TRUE;
  abd->paint_x  = abd->current_x;
  abd->current_x += length;
  abd->bpa->length = length;
  return abd->bpa;
}


/* Function:  paste_char_bt_Ascii(bpa,x,y,c,format)
 *
 * Descrip:    The paste function. Going to get at all the info
 *             (obviously) through data
 *
 *
 * Arg:           bpa [UNKN ] Undocumented argument [btPasteArea *]
 * Arg:             x [UNKN ] Undocumented argument [int]
 * Arg:             y [UNKN ] Undocumented argument [int]
 * Arg:             c [UNKN ] Undocumented argument [char]
 * Arg:        format [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 148 "asciibtcanvas.dy"
boolean paste_char_bt_Ascii(btPasteArea * bpa,int x,int y,char c,int format)
{
  Ascii_btc_Data * abd = NULL;

  abd = (Ascii_btc_Data*) bpa->canvas_data;
  
/*  fprintf(stderr,"Printing at %d %d %c\n",x,y,c); */
  abd->scratch[y][abd->paint_x + x] = c;
  return TRUE;
}

/* Function:  can_get_bt_Ascii(btc,length)
 *
 * Descrip:    Can get function
 *
 *
 * Arg:           btc [UNKN ] Undocumented argument [btCanvas *]
 * Arg:        length [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 163 "asciibtcanvas.dy"
boolean can_get_bt_Ascii(btCanvas * btc,int length)
{
  Ascii_btc_Data * abd = NULL;

  abd = (Ascii_btc_Data*) btc->canvas_data;

/*  fprintf(stderr,"Current is %d to %d\n",abd->current_x,length); */

  if( abd->current_x + length >= abd->main + abd->res_left )
    return FALSE;


/*  fprintf(stderr,"Returning TRUE\n"); */

  return TRUE;
}

/* Function:  get_ascii_left_btPasteArea(btc)
 *
 * Descrip:    The get left area function.
 *             Again, con people into thinking that we are
 *             passing a 'live' bpa
 *
 *
 * Arg:        btc [UNKN ] Undocumented argument [btCanvas *]
 *
 * Return [UNKN ]  Undocumented return value [btPasteArea *]
 *
 */
# line 186 "asciibtcanvas.dy"
btPasteArea * get_ascii_left_btPasteArea(btCanvas * btc)
{
  Ascii_btc_Data * abd = NULL;
  int length;
 
  abd = (Ascii_btc_Data*) btc->canvas_data;

  length = abd->res_left;

  if( abd->in_use == TRUE ) { 
    warn("You are already using a btPasteArea on this canvas. Only one at a time! Probably you have not freed the btPasteArea before hand");
    return NULL;
  }
  
  abd->in_use = TRUE;
  abd->paint_x  = 0;
  abd->bpa->length = length;
  return abd->bpa;
}

/* Function:  get_ascii_right_btPasteArea(btc)
 *
 * Descrip:    The get right area function.
 *             Again, con people into thinking that we are
 *             passing a 'live' bpa
 *
 *
 * Arg:        btc [UNKN ] Undocumented argument [btCanvas *]
 *
 * Return [UNKN ]  Undocumented return value [btPasteArea *]
 *
 */
# line 212 "asciibtcanvas.dy"
btPasteArea * get_ascii_right_btPasteArea(btCanvas * btc)
{
  Ascii_btc_Data * abd = NULL;
  int length;
 
  abd = (Ascii_btc_Data*) btc->canvas_data;

  length = abd->res_left + abd->main;

  if( abd->in_use == TRUE ) { 
    warn("You are already using a btPasteArea on this canvas. Only one at a time! Probably you have not freed the btPasteArea before hand");
    return NULL;
  }
  
  abd->in_use = TRUE;
  abd->paint_x  = length;
  abd->bpa->length = length;
  return abd->bpa;
}

/* Function:  ascii_next_line_btPasteArea(btc)
 *
 * Descrip:    Advancement function.
 *
 *
 * Arg:        btc [UNKN ] Undocumented argument [btCanvas *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 236 "asciibtcanvas.dy"
boolean ascii_next_line_btPasteArea(btCanvas * btc)
{
  Ascii_btc_Data * abd = NULL;
  int i,len;
 
/*  fprintf(stderr,"Going to advance\n"); */
  abd = (Ascii_btc_Data*) btc->canvas_data;

  if( abd->in_use == TRUE ) { 
    warn("You are already using a btPasteArea on this canvas, and now you are asking to advance a line. Ouch");
    return FALSE;
  }
  len = abd->res_left + abd->main + abd->res_right;

  for(i=0;i<abd->depth_scratch;i++) {
/*    fprintf(stderr,"About to print %s\n",abd->scratch[i]);*/
    fputs(abd->scratch[i],abd->ofp);
  }
  fputs("\n\n",abd->ofp);

  for(i=0;i<abd->depth_scratch;i++) {
    memset(abd->scratch[i],' ',len);
  }
  
  abd->current_x = abd->res_left;

  return TRUE;
  
}


/* Function:  free_Ascii_btpa(obj)
 *
 * Descrip:    Deconstructor for a btPasteArea we are
 *             going to make. Will be attached on construction
 *
 *
 * Arg:        obj [UNKN ] Undocumented argument [btPasteArea *]
 *
 * Return [UNKN ]  Undocumented return value [btPasteArea *]
 *
 */
# line 272 "asciibtcanvas.dy"
btPasteArea * free_Ascii_btpa(btPasteArea * obj)
{
  Ascii_btc_Data * abd;

  /* all we have to do is set to FALSE the Ascii_btc_Data pointed to
     by canvas_data
     */

  abd = (Ascii_btc_Data*) obj->canvas_data;

  abd->in_use = FALSE;

  return NULL;
}
  


/* Function:  new_Ascii_btc_Data(ofp,left,main,right,height)
 *
 * Descrip:    makes new ascii data. NB. Notice allocation of
 *             'dummy' btPasteArea and of the scratch pad.
 *
 *
 * Arg:           ofp [UNKN ] Undocumented argument [FILE *]
 * Arg:          left [UNKN ] Undocumented argument [int]
 * Arg:          main [UNKN ] Undocumented argument [int]
 * Arg:         right [UNKN ] Undocumented argument [int]
 * Arg:        height [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [Ascii_btc_Data *]
 *
 */
# line 294 "asciibtcanvas.dy"
Ascii_btc_Data * new_Ascii_btc_Data(FILE * ofp,int left,int main,int right,int height)
{
  Ascii_btc_Data * out;
  int tot;
  int i;

  tot = left + main + right +2;

  out = Ascii_btc_Data_alloc();

  out->ofp = ofp;
  out->res_left = left;
  out->main = main;
  out->res_right = right;
  out->current_x= left;
  out->paint_x = left;

  
  out->scratch = (char **) ckcalloc(height,sizeof(char *));
  for(i=0;i<height;i++) {
    out->scratch[i] = (char *) ckcalloc(tot,sizeof(char));
    memset(out->scratch[i],' ',tot-2);
    out->scratch[i][tot-1] = '\0';
    out->scratch[i][tot-2] = '\n';
  }
  
  out->depth_scratch = height;
  out->in_use = FALSE;
  out->bpa = btPasteArea_alloc();
  out->bpa->height = height;
  out->bpa->canvas_data = (void *) out;
  out->bpa->decons = free_Ascii_btpa;
  out->bpa->paste_char = paste_char_bt_Ascii;

  
  
  return out;
}

  
/* Function:  free_Ascii_btc(btc)
 *
 * Descrip:    Deconstructor for the btcanvas we are going
 *             to make. This function will be attached to
 *             the btcanvas on construction
 *
 *
 * Arg:        btc [UNKN ] Undocumented argument [btCanvas *]
 *
 * Return [UNKN ]  Undocumented return value [btCanvas *]
 *
 */
# line 340 "asciibtcanvas.dy"
btCanvas * free_Ascii_btc(btCanvas * btc)
{
  btc->canvas_data = (void *) free_Ascii_btc_Data((Ascii_btc_Data *)(btc->canvas_data));
  ckfree(btc);
  return NULL;
}


/* Function:  free_Ascii_btc_Data(obj)
 *
 * Descrip:    Specialist deconstructor really for scratch pad
 *
 *
 * Arg:        obj [WRITE] Ascii_btc_Data to be zapped [Ascii_btc_Data *]
 *
 * Return [UNKN ]  Undocumented return value [Ascii_btc_Data *]
 *
 */
# line 355 "asciibtcanvas.dy"
Ascii_btc_Data * free_Ascii_btc_Data(Ascii_btc_Data * obj)
{
  int i;

  for(i=0;i<obj->depth_scratch;i++) {
    ckfree(obj->scratch[i]);
  }

  ckfree(obj->scratch);
  ckfree(obj->bpa); /* very subtle. bpa actually is a complete dummy obj
		     * calling free_btPasteArea would be **really** bad
		     */
  ckfree(obj);

  return NULL;
}




# line 429 "asciibtcanvas.c"
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
Ascii_btc_Data * hard_link_Ascii_btc_Data(Ascii_btc_Data * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a Ascii_btc_Data object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  Ascii_btc_Data_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Ascii_btc_Data *]
 *
 */
Ascii_btc_Data * Ascii_btc_Data_alloc(void) 
{
    Ascii_btc_Data * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(Ascii_btc_Data *) ckalloc (sizeof(Ascii_btc_Data))) == NULL)    {  
      warn("Ascii_btc_Data_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->ofp = NULL; 
    out->current_x = 0;  
    out->paint_x = 0;    
    out->res_left = 0;   
    out->main = 0;   
    out->res_right = 0;  
    out->scratch = NULL; 
    out->depth_scratch = 0;  
    out->in_use = FALSE; 
    out->bpa = NULL; 


    return out;  
}    



#ifdef _cplusplus
}
#endif
