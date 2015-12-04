#ifdef _cplusplus
extern "C" {
#endif
#include "fivestatemodel.h"


# line 95 "fivestatemodel.dy"
void dump_FiveStateModel(FiveStateModel * five,FILE * ofp)
{
  int i;

  assert(five);
  assert(ofp);

  for(i=0;i<five->len;i++) {
    fprintf(ofp,"Position %d %f %f\n",i,five->unit[i]->transition[FSM_OUTBOUND2INBOUND],
	    five->unit[i]->transition[FSM_OUTBOUND2END]);
  }

}

# line 109 "fivestatemodel.dy"
void dump_FiveStateScore(FiveStateScore * five,FILE * ofp)
{
  int i;

  assert(five);
  assert(ofp);

  for(i=0;i<five->len;i++) {
    fprintf(ofp,"Position %d %c %d %d %d %d\n",i,five->unit[i]->display_char,five->unit[i]->transition[FSM_OUTBOUND2INBOUND],
	    five->unit[i]->transition[FSM_OUTBOUND2END],five->unit[i]->transition[FSM_INBOUND2MATCH],five->unit[i]->transition[FSM_MATCH2OUTBOUND]);
  }

}

# line 123 "fivestatemodel.dy"
FiveStateFrameSet * read_FiveStateFrameSet_file(char * context,char * block_str)
{
  FiveStateFrameSet * out;
  char filename[MAXLINE];
  FILE * ifp;

  assert(context);
  assert(block_str);
  sprintf(filename,"%s/%s",context,block_str);

  ifp = openfile(filename,"r");
  if( ifp == NULL ) {
    warn("Could not open %s as filename for block structure file",filename);
    return NULL;
  }
  
  out = read_FiveStateFrameSet(context,ifp);

  fclose(ifp);
  return out;

}


# line 147 "fivestatemodel.dy"
FiveStateFrameSet * read_FiveStateFrameSet(char * context,FILE * ifp)
{
  char buffer[MAXLINE];
  FiveStateFrameSet * out;
  FiveStateFrame * temp;
  ThreeStateModel * tsm;
  char * file;
  char filename[MAXLINE];

  out = FiveStateFrameSet_alloc_std();
  temp = FiveStateFrame_alloc_std();
  
  add_FiveStateFrameSet(out,temp);
  
  while( fgets(buffer,MAXLINE,ifp) != NULL ) {
    if( strstartcmp(buffer,"//") == 0 ) {
	
      temp = FiveStateFrame_alloc_std();
      add_FiveStateFrameSet(out,temp);
    } else {
      file = strtok(buffer,spacestr);
      sprintf(filename,"%s/%s",context,file);

      tsm = Wise2_HMMer2_read_ThreeStateModel(filename);
      if( tsm == NULL ) 
	fatal("Could not read %s as a HMMer2 file - die horribly!",tsm);
      add_FiveStateFrame(temp,tsm);
    }
      
  }
  return out;
}

# line 180 "fivestatemodel.dy"
void blank_FiveStateUnit(FiveStateUnit * u)
{
  int i;

  for(i=0;i<ALPHABET_SIZE;i++) {
    u->match[i] = 1.0;
    u->insert[i] = 1.0;
  } 

  for(i=0;i<FSM_TRANSITION_LENGTH;i++) {
    u->transition[i] = 0.0;
  } 

}


# line 196 "fivestatemodel.dy"
FiveStateModel * FiveStateModel_from_FiveStateFrameSet(FiveStateFrameSet * frame)
{
  FiveStateModel * out;
  FiveStateUnit *  u;
  int len;
  int max;
  int i,j,k,l;
  int f;
  FiveStateFrameIndicator * fsfi;
  ThreeStateModel * tsm;
  ThreeStateUnit * donor;

  len = 0;
  for(i=0;i<frame->len;i++) {
    for(j=0;j<frame->set[i]->len;j++) {
      len += frame->set[i]->stack[j]->len;
    }
  }
  
  assert(len > 0);

  out = FiveStateModel_alloc_len(len);

  /* the first state is a dummy state to get
     probability to get to INBOUND */
  u = FiveStateUnit_alloc();
  blank_FiveStateUnit(u);
  u->display_char = '*';
  u->transition[FSM_START2INBOUND] = 1.0;
  u->transition[FSM_INBOUND2MATCH] = 1.0/frame->set[0]->len;
  u->transition[FSM_INBOUND2INBOUND] = 1.0 - u->transition[FSM_INBOUND2MATCH]; 

  fsfi = FiveStateFrameIndicator_alloc();
  fsfi->start = 1;
  fsfi->tsm = hard_link_ThreeStateModel(frame->set[0]->stack[0]);
  add_tsm_FiveStateModel(out,fsfi);

  add_FiveStateModel(out,u);


  f = 0;
  for(i=0;i<frame->len;i++) {
    for(j=0;j<frame->set[i]->len;j++) {
      for(k=0;k<frame->set[i]->stack[j]->len;k++){
	u = FiveStateUnit_alloc();
	blank_FiveStateUnit(u);

	if( fsfi->tsm == NULL ) {
	  fsfi->tsm = hard_link_ThreeStateModel(frame->set[i]->stack[j]);
	}

	donor = frame->set[i]->stack[j]->unit[k];

	for(l=0;l<ALPHABET_SIZE;l++) {
	  u->match[l] = donor->match_emission[l];
	  u->insert[l] = donor->insert_emission[l];
	}

	u->display_char = donor->display_char;
	add_FiveStateModel(out,u);

	u->transition[FSM_START2MATCH]   = 0.0;
	u->transition[FSM_START2INSERT]  = 0.0;
	u->transition[FSM_START2DELETE]  = 0.0;

	u->transition[FSM_MATCH2MATCH]   = donor->transition[TSM_MATCH2MATCH];
	u->transition[FSM_MATCH2INSERT]  = donor->transition[TSM_MATCH2INSERT];
	u->transition[FSM_MATCH2DELETE]  = donor->transition[TSM_MATCH2DELETE];
	u->transition[FSM_MATCH2END]     = 0.0;

	u->transition[FSM_DELETE2MATCH]  = donor->transition[TSM_DELETE2MATCH];
	u->transition[FSM_DELETE2INSERT] = donor->transition[TSM_DELETE2INSERT];
	u->transition[FSM_DELETE2DELETE] = donor->transition[TSM_DELETE2DELETE];
	u->transition[FSM_DELETE2END]    = 0.0;

	u->transition[FSM_INSERT2MATCH]  = donor->transition[TSM_INSERT2MATCH];
	u->transition[FSM_INSERT2INSERT] = donor->transition[TSM_INSERT2INSERT];
	u->transition[FSM_INSERT2DELETE] = donor->transition[TSM_INSERT2DELETE];
	u->transition[FSM_INSERT2END]    = 0.0;

	/* by default, set inbound and outbound to transition on at 1.0 */
	u->transition[FSM_OUTBOUND2OUTBOUND] = 1.0;
	u->transition[FSM_INBOUND2INBOUND]   = 1.0;


	/* all the interesting stuff happens at the end of HMM */
	if( k == frame->set[i]->stack[j]->len-1 ) {
	  /* end of this stack. transition out to outbound */
	  u->transition[FSM_MATCH2OUTBOUND] = 1.0;

	  /* set fsfi end point here */
	  fsfi->end = out->len;

	  
	  /* if this is also end of the stack, handle the switch from out to in */
	  if( j == frame->set[i]->len-1 ) {
	    u->transition[FSM_INBOUND2INBOUND] = 0.0;

	    if( i == frame->len-1 ) {
	      u->transition[FSM_MATCH2END] = 1.0;
	    }


	    /* because Dynamite has to have an offset in the model, need one
	       extra units to flip to inbound */

	    u = FiveStateUnit_alloc();
	    blank_FiveStateUnit(u);
	    u->display_char = '*';

	    u->transition[FSM_OUTBOUND2INBOUND] = 1.0;



	    add_FiveStateModel(out,u);


	    u = FiveStateUnit_alloc();
	    blank_FiveStateUnit(u);
	    u->display_char = '*';


	    /* if the end of the model */

	    if( i == frame->len-1 ) {
	      u->transition[FSM_OUTBOUND2END] = 1.0;
	      u->transition[FSM_INBOUND2END] = 1.0;
	    } else {
	      u->transition[FSM_INBOUND2MATCH] = 1.0/frame->set[i+1]->len;
	      u->transition[FSM_INBOUND2INBOUND] = 1.0 - u->transition[FSM_INBOUND2MATCH]; 
	    }
	    
	    add_FiveStateModel(out,u);
	    
	  } else {
	    /* the stack continues - adjust inbound to next point */
	    u->transition[FSM_INBOUND2MATCH]   = 1.0 / frame->set[i]->len;
	    u->transition[FSM_INBOUND2INBOUND] = 1.0 - u->transition[FSM_INBOUND2MATCH]; 
	  } 

	  /* allocate and add next fsfi */

	  fsfi = FiveStateFrameIndicator_alloc();
	  fsfi->start = out->len;
	  add_tsm_FiveStateModel(out,fsfi);

	}


      } /* end of the component hmm */
    } /* end of the stack */
  } /* end of the frame */
  
  return out;
}
  

# line 353 "fivestatemodel.dy"
void fold_RandomModel_into_FiveStateModel(FiveStateModel * fsm,RandomModel * rm)
{
  register int i;
  register int j;

  assert(fsm);
  assert(rm);
  for(i=0;i<fsm->len;i++) {
    auto FiveStateUnit * tsu;
    tsu = fsm->unit[i];

    for(j=0;j<26;j++) {
      if( rm->aminoacid[j] < 0.00000001 ) {
	warn("While trying to fold in random model, amino acid %d [%c] was below zero, ignoring",j,'A'+j);
	continue;
      }
      
      tsu->match[j]  = tsu->match[j] / rm->aminoacid[j];
      tsu->insert[j] = tsu->insert[j] / rm->aminoacid[j];
    }
  }
}



# line 378 "fivestatemodel.dy"
FiveStateModel * FiveStateModel_from_flat_ThreeStateModel(ThreeStateModel * tsm)
{
  FiveStateModel * out;
  FiveStateUnit * u;
  ThreeStateUnit * donor;
  int i,j;

  out = FiveStateModel_alloc_len(tsm->len);

  for(i=0;i<tsm->len;i++) {
    donor = tsm->unit[i];
    u = FiveStateUnit_alloc();
    for(j=0;j<ALPHABET_SIZE;j++) {
      u->match[j] = donor->match_emission[j];
      u->insert[j] = donor->insert_emission[j];
    }
    u->transition[FSM_START2MATCH]   = donor->transition[TSM_START2MATCH];
    u->transition[FSM_START2INSERT]  = donor->transition[TSM_START2INSERT];
    u->transition[FSM_START2DELETE]  = donor->transition[TSM_START2DELETE];

    u->transition[FSM_MATCH2MATCH]   = donor->transition[TSM_MATCH2MATCH];
    u->transition[FSM_MATCH2INSERT]  = donor->transition[TSM_MATCH2INSERT];
    u->transition[FSM_MATCH2DELETE]  = donor->transition[TSM_MATCH2DELETE];
    u->transition[FSM_MATCH2END]     = donor->transition[TSM_MATCH2END];

    u->transition[FSM_DELETE2MATCH]  = donor->transition[TSM_DELETE2MATCH];
    u->transition[FSM_DELETE2INSERT] = donor->transition[TSM_DELETE2INSERT];
    u->transition[FSM_DELETE2DELETE] = donor->transition[TSM_DELETE2DELETE];
    u->transition[FSM_DELETE2END]    = donor->transition[TSM_DELETE2END];

    u->transition[FSM_INSERT2MATCH]  = donor->transition[TSM_INSERT2MATCH];
    u->transition[FSM_INSERT2INSERT] = donor->transition[TSM_INSERT2INSERT];
    u->transition[FSM_INSERT2DELETE] = donor->transition[TSM_INSERT2DELETE];
    u->transition[FSM_INSERT2END]    = donor->transition[TSM_INSERT2END];

    u->transition[FSM_MATCH2INBOUND] = 0.0;
    u->transition[FSM_DELETE2INBOUND] = 0.0;
    u->transition[FSM_INSERT2INBOUND] = 0.0;

    u->transition[FSM_MATCH2INBOUND] = 0.0;
    u->transition[FSM_DELETE2INBOUND] = 0.0;
    u->transition[FSM_INSERT2INBOUND] = 0.0;

    u->display_char = donor->display_char;
    add_FiveStateModel(out,u);
  }
   
  return out;

}

/* Function:  pseudo_Protein_from_FiveStateModel(tsm)
 *
 * Descrip:    Makes a protein sequence out of the display characters.
 *             Not very useful!
 *
 *
 *
 * Arg:        tsm [UNKN ] Undocumented argument [FiveStateModel *]
 *
 * Return [UNKN ]  Undocumented return value [Protein *]
 *
 */
# line 434 "fivestatemodel.dy"
Protein * pseudo_Protein_from_FiveStateModel(FiveStateModel * tsm)
{
  int i;

  Sequence * seq;

  seq = Sequence_alloc();
  seq->name = stringalloc(tsm->name);
  seq->seq = ckcalloc(tsm->len+1,sizeof(char));

  for(i=0;i<tsm->len;i++) {
    seq->seq[i] = tsm->unit[i]->display_char;
  }
  seq->seq[i]='\0';
  make_len_type_Sequence(seq);
  seq->type = SEQUENCE_PROTEIN;

  return Protein_from_Sequence(seq);
}

# line 454 "fivestatemodel.dy"
FiveStateScore * FiveStateScore_from_FiveStateModel(FiveStateModel * fsm)
{
  FiveStateScore * out;
  int i;

  out = FiveStateScore_alloc_len(fsm->len);
  
  add_FiveStateScore(out,FiveStateScoreUnit_from_FiveStateUnit(NULL,fsm->unit[0]));

  for(i=1;i<fsm->len;i++) {
    add_FiveStateScore(out,FiveStateScoreUnit_from_FiveStateUnit(fsm->unit[i-1],fsm->unit[i]));
  }

  return out;
}

# line 470 "fivestatemodel.dy"
FiveStateScoreUnit * FiveStateScoreUnit_from_FiveStateUnit(FiveStateUnit * prev,FiveStateUnit * u)
{
  FiveStateScoreUnit * out;

  out = FiveStateScoreUnit_alloc();

  Probability2Score_move(u->match,out->match,ALPHABET_SIZE);
  Probability2Score_move(u->insert,out->insert,ALPHABET_SIZE);

  out->display_char = u->display_char;

  if( prev != NULL ) {
 
    out->transition[FSM_MATCH2MATCH]     = Probability2Score(prev->transition[FSM_MATCH2MATCH]);
    out->transition[FSM_INSERT2MATCH]    = Probability2Score(prev->transition[FSM_INSERT2MATCH]);
    out->transition[FSM_DELETE2MATCH]    = Probability2Score(prev->transition[FSM_DELETE2MATCH]);
    out->transition[FSM_INBOUND2MATCH]   = Probability2Score(prev->transition[FSM_INBOUND2MATCH]);
    out->transition[FSM_OUTBOUND2MATCH]  = Probability2Score(prev->transition[FSM_OUTBOUND2MATCH]);

    out->transition[FSM_MATCH2DELETE]    = Probability2Score(prev->transition[FSM_MATCH2DELETE]);
    out->transition[FSM_INSERT2DELETE]   = Probability2Score(prev->transition[FSM_INSERT2DELETE]);
    out->transition[FSM_DELETE2DELETE]   = Probability2Score(prev->transition[FSM_DELETE2DELETE]);
    out->transition[FSM_INBOUND2DELETE]  = Probability2Score(prev->transition[FSM_INBOUND2DELETE]);
    out->transition[FSM_OUTBOUND2DELETE] = Probability2Score(prev->transition[FSM_OUTBOUND2DELETE]);

    out->transition[FSM_MATCH2INBOUND]    = Probability2Score(prev->transition[FSM_MATCH2INBOUND]);
    out->transition[FSM_INSERT2INBOUND]   = Probability2Score(prev->transition[FSM_INSERT2INBOUND]);
    out->transition[FSM_DELETE2INBOUND]   = Probability2Score(prev->transition[FSM_DELETE2INBOUND]);
    out->transition[FSM_INBOUND2INBOUND]  = Probability2Score(prev->transition[FSM_INBOUND2INBOUND]);
    out->transition[FSM_OUTBOUND2INBOUND]  = Probability2Score(prev->transition[FSM_OUTBOUND2INBOUND]);

    out->transition[FSM_MATCH2OUTBOUND]    = Probability2Score(prev->transition[FSM_MATCH2OUTBOUND]);
    out->transition[FSM_INSERT2OUTBOUND]   = Probability2Score(prev->transition[FSM_INSERT2OUTBOUND]);
    out->transition[FSM_DELETE2OUTBOUND]   = Probability2Score(prev->transition[FSM_DELETE2OUTBOUND]);
    out->transition[FSM_OUTBOUND2OUTBOUND]  = Probability2Score(prev->transition[FSM_OUTBOUND2OUTBOUND]);

    out->transition[FSM_INBOUND2END] = Probability2Score(prev->transition[FSM_INBOUND2END]);
    out->transition[FSM_OUTBOUND2END] = Probability2Score(prev->transition[FSM_OUTBOUND2END]);
    out->transition[FSM_OUTBOUND2INBOUND] = Probability2Score(prev->transition[FSM_OUTBOUND2INBOUND]);


  } else {

    out->transition[FSM_MATCH2MATCH]     = NEGI;
    out->transition[FSM_INSERT2MATCH]    = NEGI;
    out->transition[FSM_DELETE2MATCH]    = NEGI;
    out->transition[FSM_INBOUND2MATCH]   = NEGI;
    out->transition[FSM_OUTBOUND2MATCH]  = NEGI;

    out->transition[FSM_MATCH2DELETE]    = NEGI;
    out->transition[FSM_INSERT2DELETE]   = NEGI;
    out->transition[FSM_DELETE2DELETE]   = NEGI;
    out->transition[FSM_INBOUND2DELETE]  = NEGI;
    out->transition[FSM_OUTBOUND2DELETE] = NEGI;

    out->transition[FSM_MATCH2INBOUND]    = NEGI;
    out->transition[FSM_INSERT2INBOUND]   = NEGI;
    out->transition[FSM_DELETE2INBOUND]   = NEGI;
    out->transition[FSM_INBOUND2INBOUND]  = NEGI;


    out->transition[FSM_MATCH2OUTBOUND]    = NEGI;
    out->transition[FSM_INSERT2OUTBOUND]   = NEGI;
    out->transition[FSM_DELETE2OUTBOUND]   = NEGI;
    out->transition[FSM_OUTBOUND2OUTBOUND] = NEGI;

    out->transition[FSM_INBOUND2END] = NEGI;
    out->transition[FSM_OUTBOUND2END] = NEGI;
    out->transition[FSM_OUTBOUND2INBOUND] = NEGI;

  }

  out->transition[FSM_MATCH2INSERT] = Probability2Score(u->transition[FSM_MATCH2INSERT]);
  out->transition[FSM_INSERT2INSERT] = Probability2Score(u->transition[FSM_INSERT2INSERT]);
  out->transition[FSM_DELETE2INSERT] = Probability2Score(u->transition[FSM_DELETE2INSERT]);
  out->transition[FSM_INBOUND2INSERT] = Probability2Score(u->transition[FSM_INBOUND2INSERT]);
  out->transition[FSM_OUTBOUND2INSERT] = Probability2Score(u->transition[FSM_OUTBOUND2INSERT]);
  
  out->transition[FSM_START2MATCH] = Probability2Score(u->transition[FSM_START2MATCH]);
  out->transition[FSM_START2INSERT] = Probability2Score(u->transition[FSM_START2INSERT]);
  out->transition[FSM_START2DELETE] = Probability2Score(u->transition[FSM_START2DELETE]);
  out->transition[FSM_START2INBOUND] = Probability2Score(u->transition[FSM_START2INBOUND]);
  
  out->transition[FSM_MATCH2END] = Probability2Score(u->transition[FSM_MATCH2END]);
  out->transition[FSM_INSERT2END] = Probability2Score(u->transition[FSM_INSERT2END]);
  out->transition[FSM_DELETE2END] = Probability2Score(u->transition[FSM_DELETE2END]);

  return out;
}





# line 491 "fivestatemodel.c"
/* Function:  swap_FiveStateFrame(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_FiveStateFrame
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [ThreeStateModel **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_FiveStateFrame(ThreeStateModel ** list,int i,int j)  
{
    ThreeStateModel * temp;  
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_FiveStateFrame(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_FiveStateFrame which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [ThreeStateModel **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_FiveStateFrame(ThreeStateModel ** list,int left,int right,int (*comp)(ThreeStateModel * ,ThreeStateModel * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_FiveStateFrame(list,left,(left+right)/2);   
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_FiveStateFrame (list,++last,i); 
      }  
    swap_FiveStateFrame (list,left,last);    
    qsort_FiveStateFrame(list,left,last-1,comp); 
    qsort_FiveStateFrame(list,last+1,right,comp);    
}    


/* Function:  sort_FiveStateFrame(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_FiveStateFrame
 *
 *
 * Arg:         obj [UNKN ] Object containing list [FiveStateFrame *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_FiveStateFrame(FiveStateFrame * obj,int (*comp)(ThreeStateModel *, ThreeStateModel *)) 
{
    qsort_FiveStateFrame(obj->stack,0,obj->len-1,comp);  
    return;  
}    


/* Function:  expand_FiveStateFrame(obj,len)
 *
 * Descrip:    Really an internal function for add_FiveStateFrame
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [FiveStateFrame *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_FiveStateFrame(FiveStateFrame * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_FiveStateFrame called with no need"); 
      return TRUE;   
      }  


    if( (obj->stack = (ThreeStateModel ** ) ckrealloc (obj->stack,sizeof(ThreeStateModel *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_FiveStateFrame, returning FALSE");   
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_FiveStateFrame(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [FiveStateFrame *]
 * Arg:        add [OWNER] Object to add to the list [ThreeStateModel *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_FiveStateFrame(FiveStateFrame * obj,ThreeStateModel * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_FiveStateFrame(obj,obj->len + FiveStateFrameLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->stack[obj->len++]=add;  
    return TRUE; 
}    


/* Function:  flush_FiveStateFrame(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [FiveStateFrame *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_FiveStateFrame(FiveStateFrame * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->stack[i] != NULL) {  
        free_ThreeStateModel(obj->stack[i]); 
        obj->stack[i] = NULL;    
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  FiveStateFrame_alloc_std(void)
 *
 * Descrip:    Equivalent to FiveStateFrame_alloc_len(FiveStateFrameLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [FiveStateFrame *]
 *
 */
FiveStateFrame * FiveStateFrame_alloc_std(void) 
{
    return FiveStateFrame_alloc_len(FiveStateFrameLISTLENGTH);   
}    


/* Function:  FiveStateFrame_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [FiveStateFrame *]
 *
 */
FiveStateFrame * FiveStateFrame_alloc_len(int len) 
{
    FiveStateFrame * out;   /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = FiveStateFrame_alloc()) == NULL)   
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->stack = (ThreeStateModel ** ) ckcalloc (len,sizeof(ThreeStateModel *))) == NULL)    {  
      warn("Warning, ckcalloc failed in FiveStateFrame_alloc_len");  
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_FiveStateFrame(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [FiveStateFrame *]
 *
 * Return [UNKN ]  Undocumented return value [FiveStateFrame *]
 *
 */
FiveStateFrame * hard_link_FiveStateFrame(FiveStateFrame * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a FiveStateFrame object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  FiveStateFrame_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [FiveStateFrame *]
 *
 */
FiveStateFrame * FiveStateFrame_alloc(void) 
{
    FiveStateFrame * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(FiveStateFrame *) ckalloc (sizeof(FiveStateFrame))) == NULL)    {  
      warn("FiveStateFrame_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->stack = NULL;   
    out->len = out->maxlen = 0;  


    return out;  
}    


/* Function:  free_FiveStateFrame(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [FiveStateFrame *]
 *
 * Return [UNKN ]  Undocumented return value [FiveStateFrame *]
 *
 */
FiveStateFrame * free_FiveStateFrame(FiveStateFrame * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a FiveStateFrame obj. Should be trappable");    
      return NULL;   
      }  


#ifdef PTHREAD   
    assert(pthread_mutex_lock(&(obj->dynamite_mutex)) == 0); 
#endif   
    if( obj->dynamite_hard_link > 1)     {  
      return_early = 1;  
      obj->dynamite_hard_link--; 
      }  
#ifdef PTHREAD   
    assert(pthread_mutex_unlock(&(obj->dynamite_mutex)) == 0);   
#endif   
    if( return_early == 1)   
      return NULL;   
    if( obj->stack != NULL)  {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->stack[i] != NULL)   
          free_ThreeStateModel(obj->stack[i]);   
        }  
      ckfree(obj->stack);    
      }  


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_FiveStateFrameSet(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_FiveStateFrameSet
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [FiveStateFrame **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_FiveStateFrameSet(FiveStateFrame ** list,int i,int j)  
{
    FiveStateFrame * temp;   
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_FiveStateFrameSet(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_FiveStateFrameSet which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [FiveStateFrame **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_FiveStateFrameSet(FiveStateFrame ** list,int left,int right,int (*comp)(FiveStateFrame * ,FiveStateFrame * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_FiveStateFrameSet(list,left,(left+right)/2);    
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_FiveStateFrameSet (list,++last,i);  
      }  
    swap_FiveStateFrameSet (list,left,last); 
    qsort_FiveStateFrameSet(list,left,last-1,comp);  
    qsort_FiveStateFrameSet(list,last+1,right,comp); 
}    


/* Function:  sort_FiveStateFrameSet(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_FiveStateFrameSet
 *
 *
 * Arg:         obj [UNKN ] Object containing list [FiveStateFrameSet *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_FiveStateFrameSet(FiveStateFrameSet * obj,int (*comp)(FiveStateFrame *, FiveStateFrame *)) 
{
    qsort_FiveStateFrameSet(obj->set,0,obj->len-1,comp); 
    return;  
}    


/* Function:  expand_FiveStateFrameSet(obj,len)
 *
 * Descrip:    Really an internal function for add_FiveStateFrameSet
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [FiveStateFrameSet *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_FiveStateFrameSet(FiveStateFrameSet * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_FiveStateFrameSet called with no need");  
      return TRUE;   
      }  


    if( (obj->set = (FiveStateFrame ** ) ckrealloc (obj->set,sizeof(FiveStateFrame *)*len)) == NULL)     {  
      warn("ckrealloc failed for expand_FiveStateFrameSet, returning FALSE");    
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_FiveStateFrameSet(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [FiveStateFrameSet *]
 * Arg:        add [OWNER] Object to add to the list [FiveStateFrame *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_FiveStateFrameSet(FiveStateFrameSet * obj,FiveStateFrame * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_FiveStateFrameSet(obj,obj->len + FiveStateFrameSetLISTLENGTH) == FALSE) 
        return FALSE;    
      }  


    obj->set[obj->len++]=add;    
    return TRUE; 
}    


/* Function:  flush_FiveStateFrameSet(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [FiveStateFrameSet *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_FiveStateFrameSet(FiveStateFrameSet * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->set[i] != NULL)   {  
        free_FiveStateFrame(obj->set[i]);    
        obj->set[i] = NULL;  
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  FiveStateFrameSet_alloc_std(void)
 *
 * Descrip:    Equivalent to FiveStateFrameSet_alloc_len(FiveStateFrameSetLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [FiveStateFrameSet *]
 *
 */
FiveStateFrameSet * FiveStateFrameSet_alloc_std(void) 
{
    return FiveStateFrameSet_alloc_len(FiveStateFrameSetLISTLENGTH); 
}    


/* Function:  FiveStateFrameSet_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [FiveStateFrameSet *]
 *
 */
FiveStateFrameSet * FiveStateFrameSet_alloc_len(int len) 
{
    FiveStateFrameSet * out;/* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = FiveStateFrameSet_alloc()) == NULL)    
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->set = (FiveStateFrame ** ) ckcalloc (len,sizeof(FiveStateFrame *))) == NULL)    {  
      warn("Warning, ckcalloc failed in FiveStateFrameSet_alloc_len");   
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_FiveStateFrameSet(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [FiveStateFrameSet *]
 *
 * Return [UNKN ]  Undocumented return value [FiveStateFrameSet *]
 *
 */
FiveStateFrameSet * hard_link_FiveStateFrameSet(FiveStateFrameSet * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a FiveStateFrameSet object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  FiveStateFrameSet_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [FiveStateFrameSet *]
 *
 */
FiveStateFrameSet * FiveStateFrameSet_alloc(void) 
{
    FiveStateFrameSet * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(FiveStateFrameSet *) ckalloc (sizeof(FiveStateFrameSet))) == NULL)  {  
      warn("FiveStateFrameSet_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->set = NULL; 
    out->len = out->maxlen = 0;  


    return out;  
}    


/* Function:  free_FiveStateFrameSet(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [FiveStateFrameSet *]
 *
 * Return [UNKN ]  Undocumented return value [FiveStateFrameSet *]
 *
 */
FiveStateFrameSet * free_FiveStateFrameSet(FiveStateFrameSet * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a FiveStateFrameSet obj. Should be trappable"); 
      return NULL;   
      }  


#ifdef PTHREAD   
    assert(pthread_mutex_lock(&(obj->dynamite_mutex)) == 0); 
#endif   
    if( obj->dynamite_hard_link > 1)     {  
      return_early = 1;  
      obj->dynamite_hard_link--; 
      }  
#ifdef PTHREAD   
    assert(pthread_mutex_unlock(&(obj->dynamite_mutex)) == 0);   
#endif   
    if( return_early == 1)   
      return NULL;   
    if( obj->set != NULL)    {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->set[i] != NULL) 
          free_FiveStateFrame(obj->set[i]);  
        }  
      ckfree(obj->set);  
      }  


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_FiveStateUnit(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [FiveStateUnit *]
 *
 * Return [UNKN ]  Undocumented return value [FiveStateUnit *]
 *
 */
FiveStateUnit * hard_link_FiveStateUnit(FiveStateUnit * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a FiveStateUnit object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  FiveStateUnit_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [FiveStateUnit *]
 *
 */
FiveStateUnit * FiveStateUnit_alloc(void) 
{
    FiveStateUnit * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(FiveStateUnit *) ckalloc (sizeof(FiveStateUnit))) == NULL)  {  
      warn("FiveStateUnit_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    /* match[ALPHABET_SIZE] is an array: no default possible */ 
    /* insert[ALPHABET_SIZE] is an array: no default possible */ 
    /* transition[FSM_TRANSITION_LENGTH] is an array: no default possible */ 
    out->display_char = 'u'; 


    return out;  
}    


/* Function:  free_FiveStateUnit(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [FiveStateUnit *]
 *
 * Return [UNKN ]  Undocumented return value [FiveStateUnit *]
 *
 */
FiveStateUnit * free_FiveStateUnit(FiveStateUnit * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a FiveStateUnit obj. Should be trappable"); 
      return NULL;   
      }  


#ifdef PTHREAD   
    assert(pthread_mutex_lock(&(obj->dynamite_mutex)) == 0); 
#endif   
    if( obj->dynamite_hard_link > 1)     {  
      return_early = 1;  
      obj->dynamite_hard_link--; 
      }  
#ifdef PTHREAD   
    assert(pthread_mutex_unlock(&(obj->dynamite_mutex)) == 0);   
#endif   
    if( return_early == 1)   
      return NULL;   


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_FiveStateFrameIndicator(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [FiveStateFrameIndicator *]
 *
 * Return [UNKN ]  Undocumented return value [FiveStateFrameIndicator *]
 *
 */
FiveStateFrameIndicator * hard_link_FiveStateFrameIndicator(FiveStateFrameIndicator * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a FiveStateFrameIndicator object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  FiveStateFrameIndicator_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [FiveStateFrameIndicator *]
 *
 */
FiveStateFrameIndicator * FiveStateFrameIndicator_alloc(void) 
{
    FiveStateFrameIndicator * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(FiveStateFrameIndicator *) ckalloc (sizeof(FiveStateFrameIndicator))) == NULL)  {  
      warn("FiveStateFrameIndicator_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->tsm = NULL; 
    out->start = 0;  
    out->end = 0;    


    return out;  
}    


/* Function:  free_FiveStateFrameIndicator(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [FiveStateFrameIndicator *]
 *
 * Return [UNKN ]  Undocumented return value [FiveStateFrameIndicator *]
 *
 */
FiveStateFrameIndicator * free_FiveStateFrameIndicator(FiveStateFrameIndicator * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a FiveStateFrameIndicator obj. Should be trappable");   
      return NULL;   
      }  


#ifdef PTHREAD   
    assert(pthread_mutex_lock(&(obj->dynamite_mutex)) == 0); 
#endif   
    if( obj->dynamite_hard_link > 1)     {  
      return_early = 1;  
      obj->dynamite_hard_link--; 
      }  
#ifdef PTHREAD   
    assert(pthread_mutex_unlock(&(obj->dynamite_mutex)) == 0);   
#endif   
    if( return_early == 1)   
      return NULL;   
    if( obj->tsm != NULL)    
      free_ThreeStateModel(obj->tsm);    


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_FiveStateModel(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_FiveStateModel
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [FiveStateUnit **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_FiveStateModel(FiveStateUnit ** list,int i,int j)  
{
    FiveStateUnit * temp;    
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_FiveStateModel(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_FiveStateModel which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [FiveStateUnit **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_FiveStateModel(FiveStateUnit ** list,int left,int right,int (*comp)(FiveStateUnit * ,FiveStateUnit * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_FiveStateModel(list,left,(left+right)/2);   
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_FiveStateModel (list,++last,i); 
      }  
    swap_FiveStateModel (list,left,last);    
    qsort_FiveStateModel(list,left,last-1,comp); 
    qsort_FiveStateModel(list,last+1,right,comp);    
}    


/* Function:  sort_FiveStateModel(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_FiveStateModel
 *
 *
 * Arg:         obj [UNKN ] Object containing list [FiveStateModel *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_FiveStateModel(FiveStateModel * obj,int (*comp)(FiveStateUnit *, FiveStateUnit *)) 
{
    qsort_FiveStateModel(obj->unit,0,obj->len-1,comp);   
    return;  
}    


/* Function:  expand_FiveStateModel(obj,len)
 *
 * Descrip:    Really an internal function for add_FiveStateModel
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [FiveStateModel *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_FiveStateModel(FiveStateModel * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_FiveStateModel called with no need"); 
      return TRUE;   
      }  


    if( (obj->unit = (FiveStateUnit ** ) ckrealloc (obj->unit,sizeof(FiveStateUnit *)*len)) == NULL)     {  
      warn("ckrealloc failed for expand_FiveStateModel, returning FALSE");   
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_FiveStateModel(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [FiveStateModel *]
 * Arg:        add [OWNER] Object to add to the list [FiveStateUnit *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_FiveStateModel(FiveStateModel * obj,FiveStateUnit * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_FiveStateModel(obj,obj->len + FiveStateModelLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->unit[obj->len++]=add;   
    return TRUE; 
}    


/* Function:  flush_FiveStateModel(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [FiveStateModel *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_FiveStateModel(FiveStateModel * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->unit[i] != NULL)  {  
        free_FiveStateUnit(obj->unit[i]);    
        obj->unit[i] = NULL; 
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  swap_tsm_FiveStateModel(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_tsm_FiveStateModel
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [FiveStateFrameIndicator **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_tsm_FiveStateModel(FiveStateFrameIndicator ** list,int i,int j)  
{
    FiveStateFrameIndicator * temp;  
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_tsm_FiveStateModel(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_tsm_FiveStateModel which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [FiveStateFrameIndicator **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_tsm_FiveStateModel(FiveStateFrameIndicator ** list,int left,int right,int (*comp)(FiveStateFrameIndicator * ,FiveStateFrameIndicator * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_tsm_FiveStateModel(list,left,(left+right)/2);   
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_tsm_FiveStateModel (list,++last,i); 
      }  
    swap_tsm_FiveStateModel (list,left,last);    
    qsort_tsm_FiveStateModel(list,left,last-1,comp); 
    qsort_tsm_FiveStateModel(list,last+1,right,comp);    
}    


/* Function:  sort_tsm_FiveStateModel(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_tsm_FiveStateModel
 *
 *
 * Arg:         obj [UNKN ] Object containing list [FiveStateModel *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_tsm_FiveStateModel(FiveStateModel * obj,int (*comp)(FiveStateFrameIndicator *, FiveStateFrameIndicator *)) 
{
    qsort_tsm_FiveStateModel(obj->frame,0,obj->tsm_len-1,comp);  
    return;  
}    


/* Function:  expand_tsm_FiveStateModel(obj,len)
 *
 * Descrip:    Really an internal function for add_tsm_FiveStateModel
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [FiveStateModel *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_tsm_FiveStateModel(FiveStateModel * obj,int len) 
{


    if( obj->tsm_maxlen > obj->tsm_len )     {  
      warn("expand_FiveStateModeltsm_ called with no need"); 
      return TRUE;   
      }  


    if( (obj->frame = (FiveStateFrameIndicator ** ) ckrealloc (obj->frame,sizeof(FiveStateFrameIndicator *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_FiveStateModel, returning FALSE");   
      return FALSE;  
      }  
    obj->tsm_maxlen = len;   
    return TRUE; 
}    


/* Function:  add_tsm_FiveStateModel(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [FiveStateModel *]
 * Arg:        add [OWNER] Object to add to the list [FiveStateFrameIndicator *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_tsm_FiveStateModel(FiveStateModel * obj,FiveStateFrameIndicator * add) 
{
    if( obj->tsm_len >= obj->tsm_maxlen) {  
      if( expand_tsm_FiveStateModel(obj,obj->tsm_len + FiveStateModelLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->frame[obj->tsm_len++]=add;  
    return TRUE; 
}    


/* Function:  flush_tsm_FiveStateModel(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [FiveStateModel *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_tsm_FiveStateModel(FiveStateModel * obj) 
{
    int i;   


    for(i=0;i<obj->tsm_len;i++)  { /*for i over list length*/ 
      if( obj->frame[i] != NULL) {  
        free_FiveStateFrameIndicator(obj->frame[i]); 
        obj->frame[i] = NULL;    
        }  
      } /* end of for i over list length */ 


    obj->tsm_len = 0;    
    return i;    
}    


/* Function:  FiveStateModel_alloc_std(void)
 *
 * Descrip:    Equivalent to FiveStateModel_alloc_len(FiveStateModelLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [FiveStateModel *]
 *
 */
FiveStateModel * FiveStateModel_alloc_std(void) 
{
    return FiveStateModel_alloc_len(FiveStateModelLISTLENGTH);   
}    


/* Function:  FiveStateModel_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [FiveStateModel *]
 *
 */
FiveStateModel * FiveStateModel_alloc_len(int len) 
{
    FiveStateModel * out;   /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = FiveStateModel_alloc()) == NULL)   
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->unit = (FiveStateUnit ** ) ckcalloc (len,sizeof(FiveStateUnit *))) == NULL) {  
      warn("Warning, ckcalloc failed in FiveStateModel_alloc_len");  
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    if((out->frame = (FiveStateFrameIndicator ** ) ckcalloc (len,sizeof(FiveStateFrameIndicator *))) == NULL)    {  
      warn("Warning, ckcalloc failed in FiveStateModel_alloc_len");  
      return NULL;   
      }  
    out->tsm_len = 0;    
    out->tsm_maxlen = len;   


    return out;  
}    


/* Function:  hard_link_FiveStateModel(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [FiveStateModel *]
 *
 * Return [UNKN ]  Undocumented return value [FiveStateModel *]
 *
 */
FiveStateModel * hard_link_FiveStateModel(FiveStateModel * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a FiveStateModel object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  FiveStateModel_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [FiveStateModel *]
 *
 */
FiveStateModel * FiveStateModel_alloc(void) 
{
    FiveStateModel * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(FiveStateModel *) ckalloc (sizeof(FiveStateModel))) == NULL)    {  
      warn("FiveStateModel_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->name = NULL;    
    out->unit = NULL;    
    out->len = out->maxlen = 0;  
    out->frame = NULL;   
    out->tsm_len = out->tsm_maxlen = 0;  


    return out;  
}    


/* Function:  free_FiveStateModel(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [FiveStateModel *]
 *
 * Return [UNKN ]  Undocumented return value [FiveStateModel *]
 *
 */
FiveStateModel * free_FiveStateModel(FiveStateModel * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a FiveStateModel obj. Should be trappable");    
      return NULL;   
      }  


#ifdef PTHREAD   
    assert(pthread_mutex_lock(&(obj->dynamite_mutex)) == 0); 
#endif   
    if( obj->dynamite_hard_link > 1)     {  
      return_early = 1;  
      obj->dynamite_hard_link--; 
      }  
#ifdef PTHREAD   
    assert(pthread_mutex_unlock(&(obj->dynamite_mutex)) == 0);   
#endif   
    if( return_early == 1)   
      return NULL;   
    if( obj->name != NULL)   
      ckfree(obj->name);     
    if( obj->unit != NULL)   {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->unit[i] != NULL)    
          free_FiveStateUnit(obj->unit[i]);  
        }  
      ckfree(obj->unit); 
      }  
    if( obj->frame != NULL)  {  
      for(i=0;i<obj->tsm_len;i++)    {  
        if( obj->frame[i] != NULL)   
          free_FiveStateFrameIndicator(obj->frame[i]);   
        }  
      ckfree(obj->frame);    
      }  


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_FiveStateScoreUnit(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [FiveStateScoreUnit *]
 *
 * Return [UNKN ]  Undocumented return value [FiveStateScoreUnit *]
 *
 */
FiveStateScoreUnit * hard_link_FiveStateScoreUnit(FiveStateScoreUnit * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a FiveStateScoreUnit object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  FiveStateScoreUnit_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [FiveStateScoreUnit *]
 *
 */
FiveStateScoreUnit * FiveStateScoreUnit_alloc(void) 
{
    FiveStateScoreUnit * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(FiveStateScoreUnit *) ckalloc (sizeof(FiveStateScoreUnit))) == NULL)    {  
      warn("FiveStateScoreUnit_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    /* match[ALPHABET_SIZE] is an array: no default possible */ 
    /* insert[ALPHABET_SIZE] is an array: no default possible */ 
    /* transition[FSM_TRANSITION_LENGTH] is an array: no default possible */ 
    out->display_char = 'u'; 


    return out;  
}    


/* Function:  free_FiveStateScoreUnit(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [FiveStateScoreUnit *]
 *
 * Return [UNKN ]  Undocumented return value [FiveStateScoreUnit *]
 *
 */
FiveStateScoreUnit * free_FiveStateScoreUnit(FiveStateScoreUnit * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a FiveStateScoreUnit obj. Should be trappable");    
      return NULL;   
      }  


#ifdef PTHREAD   
    assert(pthread_mutex_lock(&(obj->dynamite_mutex)) == 0); 
#endif   
    if( obj->dynamite_hard_link > 1)     {  
      return_early = 1;  
      obj->dynamite_hard_link--; 
      }  
#ifdef PTHREAD   
    assert(pthread_mutex_unlock(&(obj->dynamite_mutex)) == 0);   
#endif   
    if( return_early == 1)   
      return NULL;   


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_FiveStateScore(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_FiveStateScore
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [FiveStateScoreUnit **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_FiveStateScore(FiveStateScoreUnit ** list,int i,int j)  
{
    FiveStateScoreUnit * temp;   
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_FiveStateScore(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_FiveStateScore which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [FiveStateScoreUnit **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_FiveStateScore(FiveStateScoreUnit ** list,int left,int right,int (*comp)(FiveStateScoreUnit * ,FiveStateScoreUnit * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_FiveStateScore(list,left,(left+right)/2);   
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_FiveStateScore (list,++last,i); 
      }  
    swap_FiveStateScore (list,left,last);    
    qsort_FiveStateScore(list,left,last-1,comp); 
    qsort_FiveStateScore(list,last+1,right,comp);    
}    


/* Function:  sort_FiveStateScore(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_FiveStateScore
 *
 *
 * Arg:         obj [UNKN ] Object containing list [FiveStateScore *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_FiveStateScore(FiveStateScore * obj,int (*comp)(FiveStateScoreUnit *, FiveStateScoreUnit *)) 
{
    qsort_FiveStateScore(obj->unit,0,obj->len-1,comp);   
    return;  
}    


/* Function:  expand_FiveStateScore(obj,len)
 *
 * Descrip:    Really an internal function for add_FiveStateScore
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [FiveStateScore *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_FiveStateScore(FiveStateScore * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_FiveStateScore called with no need"); 
      return TRUE;   
      }  


    if( (obj->unit = (FiveStateScoreUnit ** ) ckrealloc (obj->unit,sizeof(FiveStateScoreUnit *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_FiveStateScore, returning FALSE");   
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_FiveStateScore(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [FiveStateScore *]
 * Arg:        add [OWNER] Object to add to the list [FiveStateScoreUnit *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_FiveStateScore(FiveStateScore * obj,FiveStateScoreUnit * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_FiveStateScore(obj,obj->len + FiveStateScoreLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->unit[obj->len++]=add;   
    return TRUE; 
}    


/* Function:  flush_FiveStateScore(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [FiveStateScore *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_FiveStateScore(FiveStateScore * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->unit[i] != NULL)  {  
        free_FiveStateScoreUnit(obj->unit[i]);   
        obj->unit[i] = NULL; 
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  FiveStateScore_alloc_std(void)
 *
 * Descrip:    Equivalent to FiveStateScore_alloc_len(FiveStateScoreLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [FiveStateScore *]
 *
 */
FiveStateScore * FiveStateScore_alloc_std(void) 
{
    return FiveStateScore_alloc_len(FiveStateScoreLISTLENGTH);   
}    


/* Function:  FiveStateScore_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [FiveStateScore *]
 *
 */
FiveStateScore * FiveStateScore_alloc_len(int len) 
{
    FiveStateScore * out;   /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = FiveStateScore_alloc()) == NULL)   
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->unit = (FiveStateScoreUnit ** ) ckcalloc (len,sizeof(FiveStateScoreUnit *))) == NULL)   {  
      warn("Warning, ckcalloc failed in FiveStateScore_alloc_len");  
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_FiveStateScore(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [FiveStateScore *]
 *
 * Return [UNKN ]  Undocumented return value [FiveStateScore *]
 *
 */
FiveStateScore * hard_link_FiveStateScore(FiveStateScore * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a FiveStateScore object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  FiveStateScore_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [FiveStateScore *]
 *
 */
FiveStateScore * FiveStateScore_alloc(void) 
{
    FiveStateScore * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(FiveStateScore *) ckalloc (sizeof(FiveStateScore))) == NULL)    {  
      warn("FiveStateScore_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->name = NULL;    
    out->unit = NULL;    
    out->len = out->maxlen = 0;  


    return out;  
}    


/* Function:  free_FiveStateScore(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [FiveStateScore *]
 *
 * Return [UNKN ]  Undocumented return value [FiveStateScore *]
 *
 */
FiveStateScore * free_FiveStateScore(FiveStateScore * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a FiveStateScore obj. Should be trappable");    
      return NULL;   
      }  


#ifdef PTHREAD   
    assert(pthread_mutex_lock(&(obj->dynamite_mutex)) == 0); 
#endif   
    if( obj->dynamite_hard_link > 1)     {  
      return_early = 1;  
      obj->dynamite_hard_link--; 
      }  
#ifdef PTHREAD   
    assert(pthread_mutex_unlock(&(obj->dynamite_mutex)) == 0);   
#endif   
    if( return_early == 1)   
      return NULL;   
    if( obj->name != NULL)   
      ckfree(obj->name);     
    if( obj->unit != NULL)   {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->unit[i] != NULL)    
          free_FiveStateScoreUnit(obj->unit[i]); 
        }  
      ckfree(obj->unit); 
      }  


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif
