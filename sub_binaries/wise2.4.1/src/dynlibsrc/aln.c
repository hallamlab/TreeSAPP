#ifdef _cplusplus
extern "C" {
#endif
#include "aln.h"


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
# line 111 "aln.dy"
AlnBlock * collapsed_AlnBlock(AlnBlock * alb,int comp_row)
{
  AlnBlock * out;
  AlnColumn * temp = NULL;
  AlnColumn * prev = NULL;
  AlnColumn * ori;
  int i;

  out = AlnBlock_alloc();
  out->start = NULL;

  assert(alb);
  assert(comp_row >= 0);
  assert(out);

  for(ori=alb->start;ori != NULL;ori = ori->next) {
    if( temp != NULL && strcmp(temp->alu[comp_row]->text_label,ori->alu[comp_row]->text_label) == 0 ) { 
      /* same column; eat it up */
      for(i=0;i<ori->len;i++) {
	temp->alu[i]->end        = ori->alu[i]->end;
	temp->alu[i]->score[0]   += ori->alu[i]->score[0];
      }
      
    } else {
      /* new column */
      prev = temp;
      temp = AlnColumn_alloc_len(ori->len);
      for(i=0;i<ori->len;i++) {
	temp->alu[i]= AlnUnit_alloc();
	temp->alu[i]->text_label = ori->alu[i]->text_label;
	temp->alu[i]->start      = ori->alu[i]->start;
	temp->alu[i]->end        = ori->alu[i]->end;
	temp->alu[i]->score[0]   = ori->alu[i]->score[0];
      }
      temp->len = ori->len;
      if( prev == NULL ) {
	out->start = temp;
      } else {
	prev->next = temp;
      }
    } 
      

  }

  return out;
}

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
# line 173 "aln.dy"
boolean add_to_anchored_AlnBlock(AlnBlock * growing,AlnBlock * add)
{
  AlnColumn * anc;
  AlnColumn * new;
  AlnColumn * insert;
  AlnUnit * temp;
  AlnUnit ** prev;
  int i;

  assert(growing);
  assert(add);
  anc = growing->start;
  new = add->start;

  add_AlnBlock(growing,add->seq[1]);
  /* this should prevent freeing the sequence we are about to use */
  add->len = 1;
  prev = &(add->seq[1]->start);
  /* wait until new is inside anc */
  for(;new->alu[0]->start != anc->alu[0]->start;anc = anc->next) {
    /*fprintf(stderr,"Moving to %d-%d\n",anc->alu[0]->start,anc->alu[0]->end);*/
    if( anc->alu[0]->start > new->alu[0]->start ) {
      warn("Somehow new has overrun anc, new %d-%d and anc %d-%d",new->alu[0]->start,new->alu[0]->end,anc->alu[0]->start,anc->alu[0]->end);
      return FALSE;
    }
    /* add a blank unit onto the anc */
    temp = AlnUnit_alloc();
    temp->score[0] = 0;
    temp->start = temp->end = 0;
    temp->text_label = "INSERT";
    add_AlnColumn(anc,temp);
    *prev = temp;
    prev = &(temp->next);
  }

  /* ok - now we are at the start. We loop over for each 'matched' state */

  *prev = new->alu[1];
  for(;new->next != NULL;) {
    
    /* if the anchor state is actually an insert here, then
     *  this is an insert not in the added sequence (if it were,
     * we would have added it last time through the loop)
     *
     *  So - skip insert points
     */

    if( anc->alu[0]->start == anc->alu[0]->end ) {
      temp = AlnUnit_alloc();
      temp->score[0] = 0;
      temp->start = temp->end = new->alu[1]->end;
      temp->text_label = "INSERT";
      *prev = temp;
      prev  = &(temp->next);
      temp->next = new->alu[1];
      add_AlnColumn(anc,temp);
      anc = anc->next;
      continue;
    }




    /* sanity check that we are at the same position on the two sequences */
    if( anc->alu[0]->start != new->alu[0]->start ||
	anc->alu[0]->end != new->alu[0]->end ) {
      warn("SANITY CHECK failed! new and anchor out of sync! anchor at %d-%d and new %d-%d",anc->alu[0]->start,anc->alu[0]->end,new->alu[0]->start,new->alu[0]->end);
      return FALSE;
    }

    add_AlnColumn(anc,new->alu[1]);
    /* see whether the next state is actually an insert */
    prev = &(new->alu[1]->next);
    for(new = new->next;new != NULL && new->alu[0]->start == new->alu[0]->end;new = new->next) {


      /*fprintf(stderr,"Adding in insert states\n");*/
      /* bail out if we have reached the end of anc */
      if( anc->next == NULL ) {
	warn("HELP - run out of add positions!");
	return FALSE;
      }
      
      /* does anc already have an insert state here? */
      if( anc->next->alu[0]->start == anc->next->alu[0]->end ) {
	/* yes! - just add to the end */
	anc = anc->next;
	add_AlnColumn(anc,new->alu[1]);
      } else {
	/* we have to drop in a new alncolumn */

	insert = AlnColumn_alloc_std();
	insert->next = anc->next;
	anc->next = insert;
	/*fprintf(stderr,"There are %d states to add to\n",anc->len);*/

	for(i=0;i<anc->len-1;i++) {
	  temp = AlnUnit_alloc();
	  temp->score[0] = 0;
	  temp->start = temp->end = anc->alu[i]->end;
	  temp->text_label = "INSERT";
	  temp->next        = anc->alu[i]->next;
	  anc->alu[i]->next = temp;
	  
	  add_AlnColumn(insert,temp);
	}
	anc = anc->next; /* could be anc = insert */
	/* add the new seq */
	add_AlnColumn(insert,new->alu[1]);

	prev = &(new->alu[1]->next);
	
      } /* end of else add in an insert */
      
    } /* end of looping over possible new inserts */

    /*fprintf(stderr,"moving on one in anchor (new already moved!) %d-%d\n",anc->alu[0]->start,anc->alu[0]->end);*/

    prev = &(new->alu[1]->next);
    anc = anc->next;

  } /* end of new alignment */

  /*fprintf(stderr,"Left alignment at %d-%d\n",anc->alu[0]->start,anc->alu[0]->end);*/

  /* now fill into the end of the anchored alignment */

  for(;anc != NULL;anc = anc->next ) {
    /* add a blank unit onto the anc */
    temp = AlnUnit_alloc();
    temp->score[0] = 0;
    temp->start = temp->end = 0;
    temp->text_label = "INSERT";
    add_AlnColumn(anc,temp);
    *prev = temp;
    prev = &(temp->next);
  }

  *prev = NULL;

  return TRUE;
	
}

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
# line 324 "aln.dy"
AlnBlock * single_unit_AlnBlock(int len,char * label)
{
  AlnBlock * out;
  AlnSequence * als;
  AlnColumn * alc;
  AlnUnit * alu;
  AlnColumn ** prev;
  AlnUnit ** aprev;
  int i;

  out = AlnBlock_alloc_std();
  als = AlnSequence_alloc();
  add_AlnBlock(out,als);

  prev = &(out->start);
  aprev = &(als->start);
  for(i=-1;i<len;i++) {
    alc = AlnColumn_alloc_std();
    *prev = alc;
    alu = AlnUnit_alloc();
    alu->start = i;
    alu->end = i+1;
    alu->text_label = label;
    add_AlnColumn(alc,alu);
    *aprev = alu;

    prev = &(alc->next);
    aprev = &(alu->next);

  }

  *prev = NULL;
  *aprev = NULL;
  
  return out;
}
    
      
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
# line 378 "aln.dy"
AlnBlockList * split_AlnBlock(AlnBlock * alb,boolean (*is_spacer_column)(const AlnColumn * alc))
{
  AlnBlockList * out;
  AlnBlock * temp;
  AlnColumn * alc;
  AlnColumn * prev;
  AlnSequence * als;
  int i;

  if( is_spacer_column == NULL ) {
    warn("You must pass in function for the spacer columns to split_AlnBlock");
  }

  if( alb == NULL ) {
    warn("Unable to split a NULL alnblock!");
  }


  out = AlnBlockList_alloc_std();
  assert(out);



  for(alc=alb->start;alc != NULL;) {
    for(;alc != NULL && (*is_spacer_column)(alc) == TRUE;alc = alc->next) {
      if( alc == NULL)
	break;

      /* we need to detach this alncolumn from the entire alignment 
       * as it will stop being valid. 
       */
      
      /* this is cheating, but it is still relevant to catch problems */
      if( alc->dynamite_hard_link > 1 ) {
	warn("Detaching a alncolumn in split_AlnBlock that is being held somewhere else. This is very dangerous");
       
	/* ok - give it a fighting chance of not crashing the client */

	for(i=0;i<alc->len;i++) {
	  /* alnunits don't point to next */
	  alc->alu[i]->next = NULL;
	}

	alc->next = NULL;
      }

      /* now - free it. It *may* be being held onto by
       * someone else. Pity him...
       */
      free_AlnColumn(alc);
    }

    /* alc is now the first alncolumn of a new alb */
    if( alc == NULL) {
      if( out->len == 0 ) {
	warn("Unable to find any alignments in AlnBlock splitting!");
      }
      break; /* will return */
    }
  
    temp = AlnBlock_alloc_std();
    add_AlnBlockList(out,temp);

    /* we need to make an exact copy of the AlnSequences */

    for(i=0;i<alb->len;i++) {
      assert(alc->alu[i]);
      als = AlnSequence_alloc();
      add_AlnBlock(temp,als);
      als->data_type = alb->seq[i]->data_type;
      als->data      = alb->seq[i]->data;
      als->bio_start = alb->seq[i]->bio_start;
      als->bio_end   = alb->seq[i]->bio_end;
      
      /* unit is the unit in this column */

      als->start = alc->alu[i];
    }

    /* attach alc as the first column to temp */

    temp->start = alc;

    /* keep walking along the list now until we hit
     * a random column. All the columns and units are
     * already correctly laid out in memory
     */

    for(prev = alc;alc != NULL && (*is_spacer_column)(alc) == FALSE;) {
      prev = alc;
      alc = alc->next;
    }

    /* prev is the last column. detach it */

    prev->next = NULL;
    for(i=0;i<prev->len;i++)
      prev->alu[i]->next = NULL;
    
    /* return to eating empty columns */
  }

  /* alb is now completely defunct. We need to make it
   *  impotent with regard to the old stuff it used to have
   *
   * this basically makes alb complete empty
   */

  alb->len =0;
  alb->start = NULL;

  return out;
}

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
# line 496 "aln.dy"
int score_line_from_AlnBlock(AlnBlock * alb,int seqno)
{
  AlnUnit * alu;
  int score = 0;

  assert(alb);
  if( seqno > alb->len ) {
    warn("Asking for a score number of a sequence more than alb. giving back 0");
    return 0;
  }

  for(alu = alb->seq[seqno]->start;alu != NULL;alu = alu->next) {
    score += alu->score[0];
  }

  return score;
}
  
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
# line 522 "aln.dy"
boolean at_end_AlnColumn(AlnColumn * alc)
{
  if( alc->next == NULL )
    return 1;
  return 0;
}


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
# line 535 "aln.dy"
int bio_start_AlnUnit(AlnUnit * alu)
{
  if( alu->seq == NULL ) {
    warn("A mis-connected AlnUnit (no sequence!) - Returning a start as if start == 1");
    return alu->start + 1 + 1;
  }

  return alu->start + alu->seq->bio_start + 1;
}

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
# line 550 "aln.dy"
int bio_end_AlnUnit(AlnUnit * alu)
{
  if( alu->seq == NULL ) {
    warn("A mis-connected AlnUnit (no sequence!) - Returning a start as if start == 1");
    return alu->end + 1;
  }

  return alu->end + alu->seq->bio_start;
}



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
# line 578 "aln.dy"
AlnColumn * swallow_AlnColumn_multiple(AlnColumn * master,AlnColumn * eaten,boolean (*comp_func)(char *,char *))
{
  register int i;


  for(i=0;eaten != NULL; eaten = eaten->next) {
    if( swallow_AlnColumn(master,eaten,comp_func) == FALSE ) 
      break;
  }

  return eaten;
}  


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
# line 602 "aln.dy"
int swallow_AlnColumn_number(AlnColumn * master,AlnColumn * eaten,int num,boolean (*comp_func)(char *,char *))
{
  register int i;

  for(i=0;i < num && eaten != NULL; eaten = eaten->next) {
    if( swallow_AlnColumn(master,eaten,comp_func) == FALSE ) 
      break;
  }

  return i;
}


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
# line 626 "aln.dy"
boolean swallow_AlnColumn(AlnColumn * master,AlnColumn * eaten,boolean (*comp_func)(char *,char *))
{
  register int i;

  if( can_swallow_AlnColumn(master,eaten,comp_func) == FALSE ) {
    return FALSE;
  }

  for(i=0;i<master->len;i++) {
    master->alu[i]->end = eaten->alu[i]->end;

    /*** ok, for the moment, only add the 1st score, but eventually ***/
    /*** we should add all of them...                               ***/

    master->alu[i]->score[0] += eaten->alu[i]->score[0];
  }

  

  return TRUE;
}

/* Function:  can_swallow_AlnColumn(comp_func,master,eaten)
 *
 * Descrip:    checks to see if two columns are mergable from comp_func.
 *               First uses /identical_labels_in_AlnColumn to see if labels can be merged
 *               Then checks that starts in master are greater than starts in eaten
 *
 *
 * Arg:        comp_func [UNKN ] Undocumented argument [NullString]
 * Arg:           master [UNKN ] Undocumented argument [AlnColumn *]
 * Arg:            eaten [UNKN ] Undocumented argument [AlnColumn *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 654 "aln.dy"
boolean can_swallow_AlnColumn(AlnColumn * master,AlnColumn * eaten,boolean (*comp_func)(char *,char *))
{
  register int i;

  if( identical_labels_in_AlnColumn(master,eaten,comp_func) == FALSE ) 
     return FALSE;

  for(i=0;i<master->len;i++) 
    if( master->alu[i]->start >= eaten->alu[i]->start) {
      warn("In trying to compare to AlnColumns, have some 'eatable' starts greater than master starts %d %d in row %d",master->alu[i]->start,eaten->alu[i]->start,i);
      return FALSE;
    }

  return TRUE;
}
    

/* Function:  identical_labels_in_AlnColumn(comp_func,one,two)
 *
 * Descrip:    checks to see if two AlnColumns has mergable labels by
 *             comp_func. calls /identical_labels_in_AlnUnits for
 *             the actual comparison.
 *
 *
 * Arg:        comp_func [UNKN ] Undocumented argument [NullString]
 * Arg:              one [UNKN ] Undocumented argument [AlnColumn *]
 * Arg:              two [UNKN ] Undocumented argument [AlnColumn *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 677 "aln.dy"
boolean identical_labels_in_AlnColumn(AlnColumn * one,AlnColumn * two,boolean (*comp_func)(char *,char *))
{
  register int i;

  if( one->len != two->len ) {
    warn("Attempting to see if two AlnColumns with *different numbers of units* %d,%d are identical...serious problem",one->len,two->len);
    return FALSE;
  }

  for(i=0;i<one->len;i++) {
    if( identical_labels_in_AlnUnits(one->alu[i],two->alu[i],comp_func) == FALSE )
      return FALSE;
  }

  return TRUE;
}

/* Function:  identical_labels_in_AlnUnits(comp_func,one,two)
 *
 * Descrip:    actually calls the comp_func for the label compairson
 *
 *
 * Arg:        comp_func [UNKN ] Undocumented argument [NullString]
 * Arg:              one [UNKN ] Undocumented argument [AlnUnit *]
 * Arg:              two [UNKN ] Undocumented argument [AlnUnit *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 698 "aln.dy"
boolean identical_labels_in_AlnUnits(AlnUnit * one,AlnUnit * two,boolean (*comp_func)(char *,char *))
{
  if( (*comp_func)(one->text_label,two->text_label) == TRUE )
    return TRUE;
  else return FALSE;
}
    

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
# line 714 "aln.dy"
void replace_and_free_AlnColumn_with_one(AlnColumn * start,AlnColumn * end,AlnColumn * insert)
{
  replace_AlnColumn_with_one(start,end,insert);

  free_AlnColumn(start->next);
}



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
# line 729 "aln.dy"
void replace_AlnColumn_with_one(AlnColumn * start,AlnColumn * end,AlnColumn * insert)
{
  start->next = insert;
  insert->next = end->next;
  end->next = NULL;

}


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
# line 744 "aln.dy"
void insert_AlnColumn(AlnColumn * start,AlnColumn * insert)
{
  insert->next = start->next;
  start->next = insert;
}


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
# line 757 "aln.dy"
AlnColumn * go_back_n_AlnColumn(AlnBlock * alb,AlnColumn * start,int n)
{
  /*** really quite hacky ****/

  AlnColumn * cursor;
  AlnColumn * back;
  register int i;

  for(i=0,cursor = alb->start;i < n && cursor != NULL && cursor != start;i++,cursor = cursor->next)
    ;

  if( i < n ) {
    return NULL; /** should I post an error? **/
  }

  for(back = alb->start;cursor != NULL && cursor != start;cursor = cursor->next, back = back->next) 
    ;

  if( cursor == NULL ) {
    warn("could not find you AlnColumn in AlnBlock at all... so could not get xxx positions back");
    return NULL;
  }

  return back;
}


 /*******************************/
 /* show functions for aln block*/
 /*                             */
 /* these are really basic and  */
 /* should use alndisplay for   */
 /* ascii display               */
 /*******************************/

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
# line 800 "aln.dy"
void dump_ascii_AlnBlock(AlnBlock * alb,FILE * ofp)
{
  int i;
  AlnColumn * alc;

  for(alc=alb->start;alc != NULL;alc = alc->next) {
    fprintf(ofp,"[%d:%d \"%s\" %d]",alc->alu[0]->start,alc->alu[0]->end,alc->alu[0]->text_label,alc->alu[0]->score[0]);
 
    for(i=1;i<alc->len;i++) {
      fprintf(ofp,",[%d:%d \"%s\" %d]",alc->alu[i]->start,alc->alu[i]->end,alc->alu[i]->text_label,alc->alu[i]->score[0]);
    }

    fprintf(ofp,"\n");
  }

  fprintf(ofp,"\\\n");

}

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
# line 826 "aln.dy"
void mapped_ascii_AlnBlockList(AlnBlockList *al,double (*score_to_double)(int),FILE * ofp)
{
  int i;
  assert(al);
  assert(ofp);
  assert(score_to_double);
  for(i=0;i<al->len;i++) {
    fprintf(ofp,"\n\nAlignment %d\n",i+1);
    mapped_ascii_AlnBlock(al->alb[i],score_to_double,0,ofp);
  }
}

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
# line 845 "aln.dy"
void mapped_ascii_AlnBlock(AlnBlock * alb,double (*score_to_double)(int),int do_cumlative,FILE * ofp)
{
  int i;
  AlnColumn * alc;
  int cuml = 0;

  if( alb == NULL || score_to_double == NULL || ofp == NULL ) {
    warn("Passing in null objects to mapped_ascii_AlnBlock - unable to show!");
    return;
  }


  for(alc=alb->start;alc != NULL;alc = alc->next) {
    cuml = cuml + alc->alu[0]->score[0];
    fprintf(ofp,"%3.2f ",(double)(*score_to_double)(alc->alu[0]->score[0]));
    fprintf(ofp,"[%d:%d \"%s\" %d]",alc->alu[0]->start,alc->alu[0]->end,alc->alu[0]->text_label,alc->alu[0]->score[0]);
 
    for(i=1;i<alc->len;i++) {
      fprintf(ofp,",[%d:%d \"%s\" %d]",alc->alu[i]->start,alc->alu[i]->end,alc->alu[i]->text_label,alc->alu[i]->score[0]);
    }

    if( do_cumlative == 1 ) {
      fprintf(ofp," {%3.2f}  ",(double)(*score_to_double)(cuml),cuml);
    }

    fprintf(ofp,"\n");
  }


}

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
# line 881 "aln.dy"
AlnBlock * read_ascii_dump_AlnBlock(FILE * ifp)
{
  char buffer[MAXLINE];
  AlnBlock * out;
  AlnColumn ** attach;
  AlnColumn * new;

  out = AlnBlock_alloc_std();

  attach = &out->start;

  while( fgets(buffer,MAXLINE,ifp) != NULL ) {
    if( strstartcmp(buffer,"//") == 0 )
      break;
    
    new = read_dumped_ascii_AlnColumn_line(buffer);
    if( new == NULL ){
      warn("Unable to read entire AlnBlock. Returning no alignment");
      free_AlnBlock(out);
      return NULL;
    }

    *attach = new;
    attach = &new->next;
  }
  return out;

}

/* Function:  read_dumped_ascii_AlnColumn_line(line)
 *
 * Descrip:    Reads one line of an ascii dumped alignment
 *
 *
 * Arg:        line [UNKN ] line to be read [char *]
 *
 * Return [UNKN ]  Undocumented return value [AlnColumn *]
 *
 */
# line 916 "aln.dy"
AlnColumn * read_dumped_ascii_AlnColumn_line(char * line)
{
  char ** base, **bstr;
  AlnUnit * alu;
  AlnColumn * out;

  int start;
  int end;
  int score;
  char buffer[128]; /* implies max label 128 characters. Worrying? */

  out = AlnColumn_alloc_len(2); /** most times we are reading in pairwise alignments **/

  /** split on comma, then sscanf each one **/

  base = bstr = breakstring_protect(line,",","\"");

  for(;*bstr != NULL;bstr++) {
    sscanf(*bstr,"[%d,%d \"%s\" %d",&start,&end,buffer,&score);
    alu = AlnUnit_alloc();

    alu->start = start;
    alu->end   = end;
    alu->score[0] = score;
    alu->text_label = stringalloc(buffer);
    
    add_AlnColumn(out,alu);
  }

  return out;
}
    

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
# line 957 "aln.dy"
void show_flat_AlnBlock(AlnBlock * alb,FILE * ofp)
{
  AlnColumn * alc;
  register int i;

  for(i=0,alc = alb->start;alc != NULL;alc = alc->next,i++) { 
    fprintf(ofp,"Column %d:\n",i);
    show_flat_AlnColumn(alc,ofp);
  }
}

/* Function:  show_flat_AlnColumn(alc,ofp)
 *
 * Descrip:    sub for show_flat_AlnBlock
 *
 *
 * Arg:        alc [UNKN ] Undocumented argument [AlnColumn *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 972 "aln.dy"
void show_flat_AlnColumn(AlnColumn * alc,FILE * ofp)
{
 register int i;
 for(i=0;i<alc->len;i++) {
   fprintf(ofp,"Unit %2d- ",i);
   show_flat_AlnUnit(alc->alu[i],ofp);
 }
 fprintf(ofp,"\n");

}

/* Function:  show_flat_AlnUnit(alu,ofp)
 *
 * Descrip:    sub for show_flat_AlnUnit
 *
 *
 * Arg:        alu [UNKN ] Undocumented argument [AlnUnit *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 987 "aln.dy"
void show_flat_AlnUnit(AlnUnit * alu,FILE * ofp)
{
  fprintf(ofp,"[%4d-%4d] [%s]\n",alu->start,alu->end,alu->text_label == NULL ? "No Label" : alu->text_label);
}


/* Function:  read_flat_AlnUnit_line(line,ret_pos)
 *
 * Descrip:    Not used currently. To read in the flat output
 *             format 
 *
 *
 * Arg:           line [UNKN ] Undocumented argument [char *]
 * Arg:        ret_pos [UNKN ] Undocumented argument [int *]
 *
 * Return [UNKN ]  Undocumented return value [AlnUnit *]
 *
 */
# line 998 "aln.dy"
AlnUnit * read_flat_AlnUnit_line(char * line,int * ret_pos)
{
  AlnUnit * out;
  char buffer[MAXLINE];
  int start;
  int end;
  int pos;
  int num;

  sscanf(line,"Unit %d- Start: [%d], End: [%d] Label: number [%d] text [%s]",&pos,&start,&end,&num,buffer);

  out = AlnUnit_alloc();
  if( out == NULL )
    return NULL;

  out->start = start;
  out->end = end;
  out->label = num;
  out->text_label = stringalloc(buffer);

  if( ret_pos != NULL )
    *ret_pos = pos;

  return out;
}

 /***********************************/
 /* movement functions around Aln's */
 /***********************************/

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
# line 1031 "aln.dy"
AlnColumn * get_second_end_AlnColumn(AlnBlock * alb)
{
  AlnColumn * end = NULL;
  AlnColumn * prev = NULL;

  for(end = alb->start;end->next != NULL;prev = end,end = end->next)
    ;

  return prev;
}
  

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
# line 1047 "aln.dy"
AlnColumn * get_end_AlnColumn(AlnBlock * alb)
{
  AlnColumn * end;

  for(end = alb->start;end->next != NULL;end = end->next)
    ;

  return end;
}

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
# line 1061 "aln.dy"
boolean link_AlnUnits_AlnBlock(AlnBlock * alb)
{
  AlnSequence * aseq;
  AlnUnit * au;
  int i;

  for(i=0;i<alb->len;i++) {
    aseq = alb->seq[i];
    for(au = aseq->start;au != NULL;au++) {
      au->seq = aseq;
    }
  }

  return TRUE;
}
  


 

 /************************************/
 /* constructors/deconstructors to   */
 /* deal with linked list aspect of  */
 /* data                             */
 /************************************/

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
# line 1093 "aln.dy"
AlnColumn * new_pairwise_AlnColumn(void)
{
  AlnColumn * out;

  out = AlnColumn_alloc_len(2);
  add_AlnColumn(out,AlnUnit_alloc());
  add_AlnColumn(out,AlnUnit_alloc());

  return out;
}

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
# line 1109 "aln.dy"
AlnColumn * free_AlnColumn(AlnColumn * obj)
{

  if( obj == NULL ) {
    warn("passed a NULL object into free_AlnColumn!");
    return NULL;
  }

  if( obj->dynamite_hard_link > 1 ) {
    obj->dynamite_hard_link--;
    return NULL;
  }

  if( obj->next != NULL )
    free_AlnColumn(obj->next);
  if( obj->alu != NULL )
    ckfree(obj->alu); /*** units free'd from linked list ***/

  ckfree(obj);

  return NULL;
}

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
# line 1137 "aln.dy"
AlnUnit * free_AlnUnit(AlnUnit * obj)
{
  if( obj == NULL )
    return NULL;

  if( obj->dynamite_hard_link > 1 ) {
    obj->dynamite_hard_link--;
    return NULL;
  }

  if( obj->next != NULL )
    free_AlnUnit(obj->next);

  ckfree(obj);
  return NULL;
}


# line 1224 "aln.c"
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
AlnUnit * hard_link_AlnUnit(AlnUnit * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a AlnUnit object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  AlnUnit_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [AlnUnit *]
 *
 */
AlnUnit * AlnUnit_alloc(void) 
{
    AlnUnit * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(AlnUnit *) ckalloc (sizeof(AlnUnit))) == NULL)  {  
      warn("AlnUnit_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->start = 0;  
    out->end = 0;    
    out->label = 0;  
    /* score[AlnUnitSCORENUMBER] is an array: no default possible */ 
    out->in_column = TRUE;   


    return out;  
}    


/* Function:  swap_AlnColumn(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_AlnColumn
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [AlnUnit **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_AlnColumn(AlnUnit ** list,int i,int j)  
{
    AlnUnit * temp;  
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_AlnColumn(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_AlnColumn which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [AlnUnit **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_AlnColumn(AlnUnit ** list,int left,int right,int (*comp)(AlnUnit * ,AlnUnit * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_AlnColumn(list,left,(left+right)/2);    
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_AlnColumn (list,++last,i);  
      }  
    swap_AlnColumn (list,left,last); 
    qsort_AlnColumn(list,left,last-1,comp);  
    qsort_AlnColumn(list,last+1,right,comp); 
}    


/* Function:  sort_AlnColumn(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_AlnColumn
 *
 *
 * Arg:         obj [UNKN ] Object containing list [AlnColumn *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_AlnColumn(AlnColumn * obj,int (*comp)(AlnUnit *, AlnUnit *)) 
{
    qsort_AlnColumn(obj->alu,0,obj->len-1,comp); 
    return;  
}    


/* Function:  expand_AlnColumn(obj,len)
 *
 * Descrip:    Really an internal function for add_AlnColumn
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [AlnColumn *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_AlnColumn(AlnColumn * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_AlnColumn called with no need");  
      return TRUE;   
      }  


    if( (obj->alu = (AlnUnit ** ) ckrealloc (obj->alu,sizeof(AlnUnit *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_AlnColumn, returning FALSE");    
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


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
/* will expand function if necessary */ 
boolean add_AlnColumn(AlnColumn * obj,AlnUnit * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_AlnColumn(obj,obj->len + AlnColumnLISTLENGTH) == FALSE) 
        return FALSE;    
      }  


    obj->alu[obj->len++]=add;    
    return TRUE; 
}    


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
int flush_AlnColumn(AlnColumn * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->alu[i] != NULL)   {  
        free_AlnUnit(obj->alu[i]);   
        obj->alu[i] = NULL;  
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  AlnColumn_alloc_std(void)
 *
 * Descrip:    Equivalent to AlnColumn_alloc_len(AlnColumnLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [AlnColumn *]
 *
 */
AlnColumn * AlnColumn_alloc_std(void) 
{
    return AlnColumn_alloc_len(AlnColumnLISTLENGTH); 
}    


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
AlnColumn * AlnColumn_alloc_len(int len) 
{
    AlnColumn * out;/* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = AlnColumn_alloc()) == NULL)    
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->alu = (AlnUnit ** ) ckcalloc (len,sizeof(AlnUnit *))) == NULL)  {  
      warn("Warning, ckcalloc failed in AlnColumn_alloc_len");   
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


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
AlnColumn * hard_link_AlnColumn(AlnColumn * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a AlnColumn object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  AlnColumn_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [AlnColumn *]
 *
 */
AlnColumn * AlnColumn_alloc(void) 
{
    AlnColumn * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(AlnColumn *) ckalloc (sizeof(AlnColumn))) == NULL)  {  
      warn("AlnColumn_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->alu = NULL; 
    out->len = out->maxlen = 0;  


    return out;  
}    


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
AlnSequence * hard_link_AlnSequence(AlnSequence * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a AlnSequence object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  AlnSequence_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [AlnSequence *]
 *
 */
AlnSequence * AlnSequence_alloc(void) 
{
    AlnSequence * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(AlnSequence *) ckalloc (sizeof(AlnSequence))) == NULL)  {  
      warn("AlnSequence_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->start = NULL;   
    out->data_type = 0;  
    out->bio_start = 1;  
    out->bio_end = -1;   


    return out;  
}    


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
AlnSequence * free_AlnSequence(AlnSequence * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a AlnSequence obj. Should be trappable");   
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
    if( obj->start != NULL)  
      free_AlnUnit(obj->start);  
    /* obj->data is linked in */ 


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_AlnBlock(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_AlnBlock
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [AlnSequence **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_AlnBlock(AlnSequence ** list,int i,int j)  
{
    AlnSequence * temp;  
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_AlnBlock(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_AlnBlock which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [AlnSequence **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_AlnBlock(AlnSequence ** list,int left,int right,int (*comp)(AlnSequence * ,AlnSequence * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_AlnBlock(list,left,(left+right)/2); 
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_AlnBlock (list,++last,i);   
      }  
    swap_AlnBlock (list,left,last);  
    qsort_AlnBlock(list,left,last-1,comp);   
    qsort_AlnBlock(list,last+1,right,comp);  
}    


/* Function:  sort_AlnBlock(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_AlnBlock
 *
 *
 * Arg:         obj [UNKN ] Object containing list [AlnBlock *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_AlnBlock(AlnBlock * obj,int (*comp)(AlnSequence *, AlnSequence *)) 
{
    qsort_AlnBlock(obj->seq,0,obj->len-1,comp);  
    return;  
}    


/* Function:  expand_AlnBlock(obj,len)
 *
 * Descrip:    Really an internal function for add_AlnBlock
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [AlnBlock *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_AlnBlock(AlnBlock * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_AlnBlock called with no need");   
      return TRUE;   
      }  


    if( (obj->seq = (AlnSequence ** ) ckrealloc (obj->seq,sizeof(AlnSequence *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_AlnBlock, returning FALSE"); 
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


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
/* will expand function if necessary */ 
boolean add_AlnBlock(AlnBlock * obj,AlnSequence * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_AlnBlock(obj,obj->len + AlnBlockLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->seq[obj->len++]=add;    
    return TRUE; 
}    


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
int flush_AlnBlock(AlnBlock * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->seq[i] != NULL)   {  
        free_AlnSequence(obj->seq[i]);   
        obj->seq[i] = NULL;  
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  AlnBlock_alloc_std(void)
 *
 * Descrip:    Equivalent to AlnBlock_alloc_len(AlnBlockLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [AlnBlock *]
 *
 */
AlnBlock * AlnBlock_alloc_std(void) 
{
    return AlnBlock_alloc_len(AlnBlockLISTLENGTH);   
}    


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
AlnBlock * AlnBlock_alloc_len(int len) 
{
    AlnBlock * out; /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = AlnBlock_alloc()) == NULL) 
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->seq = (AlnSequence ** ) ckcalloc (len,sizeof(AlnSequence *))) == NULL)  {  
      warn("Warning, ckcalloc failed in AlnBlock_alloc_len");    
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


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
AlnBlock * hard_link_AlnBlock(AlnBlock * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a AlnBlock object: passed a NULL object");    
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  AlnBlock_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [AlnBlock *]
 *
 */
AlnBlock * AlnBlock_alloc(void) 
{
    AlnBlock * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(AlnBlock *) ckalloc (sizeof(AlnBlock))) == NULL)    {  
      warn("AlnBlock_alloc failed ");    
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->start = NULL;   
    out->seq = NULL; 
    out->len = out->maxlen = 0;  
    out->length = 0; 
    out->score = 0;  


    return out;  
}    


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
AlnBlock * free_AlnBlock(AlnBlock * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a AlnBlock obj. Should be trappable");  
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
    if( obj->start != NULL)  
      free_AlnColumn(obj->start);    
    if( obj->seq != NULL)    {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->seq[i] != NULL) 
          free_AlnSequence(obj->seq[i]); 
        }  
      ckfree(obj->seq);  
      }  


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_AlnBlockList(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_AlnBlockList
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [AlnBlock **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_AlnBlockList(AlnBlock ** list,int i,int j)  
{
    AlnBlock * temp; 
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_AlnBlockList(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_AlnBlockList which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [AlnBlock **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_AlnBlockList(AlnBlock ** list,int left,int right,int (*comp)(AlnBlock * ,AlnBlock * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_AlnBlockList(list,left,(left+right)/2); 
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_AlnBlockList (list,++last,i);   
      }  
    swap_AlnBlockList (list,left,last);  
    qsort_AlnBlockList(list,left,last-1,comp);   
    qsort_AlnBlockList(list,last+1,right,comp);  
}    


/* Function:  sort_AlnBlockList(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_AlnBlockList
 *
 *
 * Arg:         obj [UNKN ] Object containing list [AlnBlockList *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_AlnBlockList(AlnBlockList * obj,int (*comp)(AlnBlock *, AlnBlock *)) 
{
    qsort_AlnBlockList(obj->alb,0,obj->len-1,comp);  
    return;  
}    


/* Function:  expand_AlnBlockList(obj,len)
 *
 * Descrip:    Really an internal function for add_AlnBlockList
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [AlnBlockList *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_AlnBlockList(AlnBlockList * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_AlnBlockList called with no need");   
      return TRUE;   
      }  


    if( (obj->alb = (AlnBlock ** ) ckrealloc (obj->alb,sizeof(AlnBlock *)*len)) == NULL)     {  
      warn("ckrealloc failed for expand_AlnBlockList, returning FALSE"); 
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


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
/* will expand function if necessary */ 
boolean add_AlnBlockList(AlnBlockList * obj,AlnBlock * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_AlnBlockList(obj,obj->len + AlnBlockListLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->alb[obj->len++]=add;    
    return TRUE; 
}    


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
int flush_AlnBlockList(AlnBlockList * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->alb[i] != NULL)   {  
        free_AlnBlock(obj->alb[i]);  
        obj->alb[i] = NULL;  
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  AlnBlockList_alloc_std(void)
 *
 * Descrip:    Equivalent to AlnBlockList_alloc_len(AlnBlockListLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [AlnBlockList *]
 *
 */
AlnBlockList * AlnBlockList_alloc_std(void) 
{
    return AlnBlockList_alloc_len(AlnBlockListLISTLENGTH);   
}    


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
AlnBlockList * AlnBlockList_alloc_len(int len) 
{
    AlnBlockList * out; /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = AlnBlockList_alloc()) == NULL) 
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->alb = (AlnBlock ** ) ckcalloc (len,sizeof(AlnBlock *))) == NULL)    {  
      warn("Warning, ckcalloc failed in AlnBlockList_alloc_len");    
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


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
AlnBlockList * hard_link_AlnBlockList(AlnBlockList * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a AlnBlockList object: passed a NULL object");    
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  AlnBlockList_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [AlnBlockList *]
 *
 */
AlnBlockList * AlnBlockList_alloc(void) 
{
    AlnBlockList * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(AlnBlockList *) ckalloc (sizeof(AlnBlockList))) == NULL)    {  
      warn("AlnBlockList_alloc failed ");    
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->alb = NULL; 
    out->len = out->maxlen = 0;  


    return out;  
}    


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
AlnBlockList * free_AlnBlockList(AlnBlockList * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a AlnBlockList obj. Should be trappable");  
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
    if( obj->alb != NULL)    {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->alb[i] != NULL) 
          free_AlnBlock(obj->alb[i]);    
        }  
      ckfree(obj->alb);  
      }  


    ckfree(obj); 
    return NULL; 
}    


/* Function:  replace_start_AlnBlock(obj,start)
 *
 * Descrip:    Replace member variable start
 *             For use principly by API functions
 *
 *
 * Arg:          obj [UNKN ] Object holding the variable [AlnBlock *]
 * Arg:        start [OWNER] New value of the variable [AlnColumn *]
 *
 * Return [SOFT ]  member variable start [boolean]
 *
 */
boolean replace_start_AlnBlock(AlnBlock * obj,AlnColumn * start) 
{
    if( obj == NULL)     {  
      warn("In replacement function start for object AlnBlock, got a NULL object");  
      return FALSE;  
      }  
    obj->start = start;  
    return TRUE; 
}    


/* Function:  access_start_AlnBlock(obj)
 *
 * Descrip:    Access member variable start
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [AlnBlock *]
 *
 * Return [SOFT ]  member variable start [AlnColumn *]
 *
 */
AlnColumn * access_start_AlnBlock(AlnBlock * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function start for object AlnBlock, got a NULL object"); 
      return NULL;   
      }  
    return obj->start;   
}    


/* Function:  access_seq_AlnBlock(obj,i)
 *
 * Descrip:    Access members stored in the seq list
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the list [AlnBlock *]
 * Arg:          i [UNKN ] Position in the list [int]
 *
 * Return [SOFT ]  Element of the list [AlnSequence *]
 *
 */
AlnSequence * access_seq_AlnBlock(AlnBlock * obj,int i) 
{
    if( obj == NULL)     {  
      warn("In accessor function seq for object AlnBlock, got a NULL object");   
      return NULL;   
      }  
    if( obj->len <= i )  {  
      warn("In accessor function seq for object AlnBlock, index %%d is greater than list length %%d",i,obj->len);    
      return NULL;   
      }  
    return obj->seq[i];  
}    


/* Function:  length_seq_AlnBlock(obj)
 *
 * Descrip:    discover the length of the list
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the list [AlnBlock *]
 *
 * Return [UNKN ]  length of the list [int]
 *
 */
int length_seq_AlnBlock(AlnBlock * obj) 
{
    if( obj == NULL)     {  
      warn("In length function seq for object AlnBlock, got a NULL object"); 
      return -1;     
      }  
    return obj->len;     
}    


/* Function:  replace_length_AlnBlock(obj,length)
 *
 * Descrip:    Replace member variable length
 *             For use principly by API functions
 *
 *
 * Arg:           obj [UNKN ] Object holding the variable [AlnBlock *]
 * Arg:        length [OWNER] New value of the variable [int]
 *
 * Return [SOFT ]  member variable length [boolean]
 *
 */
boolean replace_length_AlnBlock(AlnBlock * obj,int length) 
{
    if( obj == NULL)     {  
      warn("In replacement function length for object AlnBlock, got a NULL object"); 
      return FALSE;  
      }  
    obj->length = length;    
    return TRUE; 
}    


/* Function:  access_length_AlnBlock(obj)
 *
 * Descrip:    Access member variable length
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [AlnBlock *]
 *
 * Return [SOFT ]  member variable length [int]
 *
 */
int access_length_AlnBlock(AlnBlock * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function length for object AlnBlock, got a NULL object");    
      return 0;  
      }  
    return obj->length;  
}    


/* Function:  replace_score_AlnBlock(obj,score)
 *
 * Descrip:    Replace member variable score
 *             For use principly by API functions
 *
 *
 * Arg:          obj [UNKN ] Object holding the variable [AlnBlock *]
 * Arg:        score [OWNER] New value of the variable [int]
 *
 * Return [SOFT ]  member variable score [boolean]
 *
 */
boolean replace_score_AlnBlock(AlnBlock * obj,int score) 
{
    if( obj == NULL)     {  
      warn("In replacement function score for object AlnBlock, got a NULL object");  
      return FALSE;  
      }  
    obj->score = score;  
    return TRUE; 
}    


/* Function:  access_score_AlnBlock(obj)
 *
 * Descrip:    Access member variable score
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [AlnBlock *]
 *
 * Return [SOFT ]  member variable score [int]
 *
 */
int access_score_AlnBlock(AlnBlock * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function score for object AlnBlock, got a NULL object"); 
      return 0;  
      }  
    return obj->score;   
}    


/* Function:  access_alu_AlnColumn(obj,i)
 *
 * Descrip:    Access members stored in the alu list
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the list [AlnColumn *]
 * Arg:          i [UNKN ] Position in the list [int]
 *
 * Return [SOFT ]  Element of the list [AlnUnit *]
 *
 */
AlnUnit * access_alu_AlnColumn(AlnColumn * obj,int i) 
{
    if( obj == NULL)     {  
      warn("In accessor function alu for object AlnColumn, got a NULL object");  
      return NULL;   
      }  
    if( obj->len <= i )  {  
      warn("In accessor function alu for object AlnColumn, index %%d is greater than list length %%d",i,obj->len);   
      return NULL;   
      }  
    return obj->alu[i];  
}    


/* Function:  length_alu_AlnColumn(obj)
 *
 * Descrip:    discover the length of the list
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the list [AlnColumn *]
 *
 * Return [UNKN ]  length of the list [int]
 *
 */
int length_alu_AlnColumn(AlnColumn * obj) 
{
    if( obj == NULL)     {  
      warn("In length function alu for object AlnColumn, got a NULL object");    
      return -1;     
      }  
    return obj->len;     
}    


/* Function:  replace_next_AlnColumn(obj,next)
 *
 * Descrip:    Replace member variable next
 *             For use principly by API functions
 *
 *
 * Arg:         obj [UNKN ] Object holding the variable [AlnColumn *]
 * Arg:        next [OWNER] New value of the variable [AlnColumn *]
 *
 * Return [SOFT ]  member variable next [boolean]
 *
 */
boolean replace_next_AlnColumn(AlnColumn * obj,AlnColumn * next) 
{
    if( obj == NULL)     {  
      warn("In replacement function next for object AlnColumn, got a NULL object");  
      return FALSE;  
      }  
    obj->next = next;    
    return TRUE; 
}    


/* Function:  access_next_AlnColumn(obj)
 *
 * Descrip:    Access member variable next
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [AlnColumn *]
 *
 * Return [SOFT ]  member variable next [AlnColumn *]
 *
 */
AlnColumn * access_next_AlnColumn(AlnColumn * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function next for object AlnColumn, got a NULL object"); 
      return NULL;   
      }  
    return obj->next;    
}    


/* Function:  replace_start_AlnUnit(obj,start)
 *
 * Descrip:    Replace member variable start
 *             For use principly by API functions
 *
 *
 * Arg:          obj [UNKN ] Object holding the variable [AlnUnit *]
 * Arg:        start [OWNER] New value of the variable [int]
 *
 * Return [SOFT ]  member variable start [boolean]
 *
 */
boolean replace_start_AlnUnit(AlnUnit * obj,int start) 
{
    if( obj == NULL)     {  
      warn("In replacement function start for object AlnUnit, got a NULL object");   
      return FALSE;  
      }  
    obj->start = start;  
    return TRUE; 
}    


/* Function:  access_start_AlnUnit(obj)
 *
 * Descrip:    Access member variable start
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [AlnUnit *]
 *
 * Return [SOFT ]  member variable start [int]
 *
 */
int access_start_AlnUnit(AlnUnit * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function start for object AlnUnit, got a NULL object");  
      return 0;  
      }  
    return obj->start;   
}    


/* Function:  replace_end_AlnUnit(obj,end)
 *
 * Descrip:    Replace member variable end
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [AlnUnit *]
 * Arg:        end [OWNER] New value of the variable [int]
 *
 * Return [SOFT ]  member variable end [boolean]
 *
 */
boolean replace_end_AlnUnit(AlnUnit * obj,int end) 
{
    if( obj == NULL)     {  
      warn("In replacement function end for object AlnUnit, got a NULL object"); 
      return FALSE;  
      }  
    obj->end = end;  
    return TRUE; 
}    


/* Function:  access_end_AlnUnit(obj)
 *
 * Descrip:    Access member variable end
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [AlnUnit *]
 *
 * Return [SOFT ]  member variable end [int]
 *
 */
int access_end_AlnUnit(AlnUnit * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function end for object AlnUnit, got a NULL object");    
      return 0;  
      }  
    return obj->end;     
}    


/* Function:  replace_label_AlnUnit(obj,label)
 *
 * Descrip:    Replace member variable label
 *             For use principly by API functions
 *
 *
 * Arg:          obj [UNKN ] Object holding the variable [AlnUnit *]
 * Arg:        label [OWNER] New value of the variable [int]
 *
 * Return [SOFT ]  member variable label [boolean]
 *
 */
boolean replace_label_AlnUnit(AlnUnit * obj,int label) 
{
    if( obj == NULL)     {  
      warn("In replacement function label for object AlnUnit, got a NULL object");   
      return FALSE;  
      }  
    obj->label = label;  
    return TRUE; 
}    


/* Function:  access_label_AlnUnit(obj)
 *
 * Descrip:    Access member variable label
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [AlnUnit *]
 *
 * Return [SOFT ]  member variable label [int]
 *
 */
int access_label_AlnUnit(AlnUnit * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function label for object AlnUnit, got a NULL object");  
      return 0;  
      }  
    return obj->label;   
}    


/* Function:  replace_text_label_AlnUnit(obj,text_label)
 *
 * Descrip:    Replace member variable text_label
 *             For use principly by API functions
 *
 *
 * Arg:               obj [UNKN ] Object holding the variable [AlnUnit *]
 * Arg:        text_label [OWNER] New value of the variable [char *]
 *
 * Return [SOFT ]  member variable text_label [boolean]
 *
 */
boolean replace_text_label_AlnUnit(AlnUnit * obj,char * text_label) 
{
    if( obj == NULL)     {  
      warn("In replacement function text_label for object AlnUnit, got a NULL object");  
      return FALSE;  
      }  
    obj->text_label = text_label;    
    return TRUE; 
}    


/* Function:  access_text_label_AlnUnit(obj)
 *
 * Descrip:    Access member variable text_label
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [AlnUnit *]
 *
 * Return [SOFT ]  member variable text_label [char *]
 *
 */
char * access_text_label_AlnUnit(AlnUnit * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function text_label for object AlnUnit, got a NULL object"); 
      return NULL;   
      }  
    return obj->text_label;  
}    


/* Function:  replace_next_AlnUnit(obj,next)
 *
 * Descrip:    Replace member variable next
 *             For use principly by API functions
 *
 *
 * Arg:         obj [UNKN ] Object holding the variable [AlnUnit *]
 * Arg:        next [OWNER] New value of the variable [AlnUnit *]
 *
 * Return [SOFT ]  member variable next [boolean]
 *
 */
boolean replace_next_AlnUnit(AlnUnit * obj,AlnUnit * next) 
{
    if( obj == NULL)     {  
      warn("In replacement function next for object AlnUnit, got a NULL object");    
      return FALSE;  
      }  
    obj->next = next;    
    return TRUE; 
}    


/* Function:  access_next_AlnUnit(obj)
 *
 * Descrip:    Access member variable next
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [AlnUnit *]
 *
 * Return [SOFT ]  member variable next [AlnUnit *]
 *
 */
AlnUnit * access_next_AlnUnit(AlnUnit * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function next for object AlnUnit, got a NULL object");   
      return NULL;   
      }  
    return obj->next;    
}    


/* Function:  replace_in_column_AlnUnit(obj,in_column)
 *
 * Descrip:    Replace member variable in_column
 *             For use principly by API functions
 *
 *
 * Arg:              obj [UNKN ] Object holding the variable [AlnUnit *]
 * Arg:        in_column [OWNER] New value of the variable [boolean]
 *
 * Return [SOFT ]  member variable in_column [boolean]
 *
 */
boolean replace_in_column_AlnUnit(AlnUnit * obj,boolean in_column) 
{
    if( obj == NULL)     {  
      warn("In replacement function in_column for object AlnUnit, got a NULL object");   
      return FALSE;  
      }  
    obj->in_column = in_column;  
    return TRUE; 
}    


/* Function:  access_in_column_AlnUnit(obj)
 *
 * Descrip:    Access member variable in_column
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [AlnUnit *]
 *
 * Return [SOFT ]  member variable in_column [boolean]
 *
 */
boolean access_in_column_AlnUnit(AlnUnit * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function in_column for object AlnUnit, got a NULL object");  
      return FALSE;  
      }  
    return obj->in_column;   
}    


/* Function:  replace_seq_AlnUnit(obj,seq)
 *
 * Descrip:    Replace member variable seq
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [AlnUnit *]
 * Arg:        seq [OWNER] New value of the variable [AlnSequence *]
 *
 * Return [SOFT ]  member variable seq [boolean]
 *
 */
boolean replace_seq_AlnUnit(AlnUnit * obj,AlnSequence * seq) 
{
    if( obj == NULL)     {  
      warn("In replacement function seq for object AlnUnit, got a NULL object"); 
      return FALSE;  
      }  
    obj->seq = seq;  
    return TRUE; 
}    


/* Function:  access_seq_AlnUnit(obj)
 *
 * Descrip:    Access member variable seq
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [AlnUnit *]
 *
 * Return [SOFT ]  member variable seq [AlnSequence *]
 *
 */
AlnSequence * access_seq_AlnUnit(AlnUnit * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function seq for object AlnUnit, got a NULL object");    
      return NULL;   
      }  
    return obj->seq;     
}    


/* Function:  replace_start_AlnSequence(obj,start)
 *
 * Descrip:    Replace member variable start
 *             For use principly by API functions
 *
 *
 * Arg:          obj [UNKN ] Object holding the variable [AlnSequence *]
 * Arg:        start [OWNER] New value of the variable [AlnUnit *]
 *
 * Return [SOFT ]  member variable start [boolean]
 *
 */
boolean replace_start_AlnSequence(AlnSequence * obj,AlnUnit * start) 
{
    if( obj == NULL)     {  
      warn("In replacement function start for object AlnSequence, got a NULL object");   
      return FALSE;  
      }  
    obj->start = start;  
    return TRUE; 
}    


/* Function:  access_start_AlnSequence(obj)
 *
 * Descrip:    Access member variable start
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [AlnSequence *]
 *
 * Return [SOFT ]  member variable start [AlnUnit *]
 *
 */
AlnUnit * access_start_AlnSequence(AlnSequence * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function start for object AlnSequence, got a NULL object");  
      return NULL;   
      }  
    return obj->start;   
}    


/* Function:  replace_data_type_AlnSequence(obj,data_type)
 *
 * Descrip:    Replace member variable data_type
 *             For use principly by API functions
 *
 *
 * Arg:              obj [UNKN ] Object holding the variable [AlnSequence *]
 * Arg:        data_type [OWNER] New value of the variable [int]
 *
 * Return [SOFT ]  member variable data_type [boolean]
 *
 */
boolean replace_data_type_AlnSequence(AlnSequence * obj,int data_type) 
{
    if( obj == NULL)     {  
      warn("In replacement function data_type for object AlnSequence, got a NULL object");   
      return FALSE;  
      }  
    obj->data_type = data_type;  
    return TRUE; 
}    


/* Function:  access_data_type_AlnSequence(obj)
 *
 * Descrip:    Access member variable data_type
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [AlnSequence *]
 *
 * Return [SOFT ]  member variable data_type [int]
 *
 */
int access_data_type_AlnSequence(AlnSequence * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function data_type for object AlnSequence, got a NULL object");  
      return 0;  
      }  
    return obj->data_type;   
}    


/* Function:  replace_data_AlnSequence(obj,data)
 *
 * Descrip:    Replace member variable data
 *             For use principly by API functions
 *
 *
 * Arg:         obj [UNKN ] Object holding the variable [AlnSequence *]
 * Arg:        data [OWNER] New value of the variable [void *]
 *
 * Return [SOFT ]  member variable data [boolean]
 *
 */
boolean replace_data_AlnSequence(AlnSequence * obj,void * data) 
{
    if( obj == NULL)     {  
      warn("In replacement function data for object AlnSequence, got a NULL object");    
      return FALSE;  
      }  
    obj->data = data;    
    return TRUE; 
}    


/* Function:  access_data_AlnSequence(obj)
 *
 * Descrip:    Access member variable data
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [AlnSequence *]
 *
 * Return [SOFT ]  member variable data [void *]
 *
 */
void * access_data_AlnSequence(AlnSequence * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function data for object AlnSequence, got a NULL object");   
      return NULL;   
      }  
    return obj->data;    
}    


/* Function:  replace_bio_start_AlnSequence(obj,bio_start)
 *
 * Descrip:    Replace member variable bio_start
 *             For use principly by API functions
 *
 *
 * Arg:              obj [UNKN ] Object holding the variable [AlnSequence *]
 * Arg:        bio_start [OWNER] New value of the variable [int]
 *
 * Return [SOFT ]  member variable bio_start [boolean]
 *
 */
boolean replace_bio_start_AlnSequence(AlnSequence * obj,int bio_start) 
{
    if( obj == NULL)     {  
      warn("In replacement function bio_start for object AlnSequence, got a NULL object");   
      return FALSE;  
      }  
    obj->bio_start = bio_start;  
    return TRUE; 
}    


/* Function:  access_bio_start_AlnSequence(obj)
 *
 * Descrip:    Access member variable bio_start
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [AlnSequence *]
 *
 * Return [SOFT ]  member variable bio_start [int]
 *
 */
int access_bio_start_AlnSequence(AlnSequence * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function bio_start for object AlnSequence, got a NULL object");  
      return 0;  
      }  
    return obj->bio_start;   
}    


/* Function:  replace_bio_end_AlnSequence(obj,bio_end)
 *
 * Descrip:    Replace member variable bio_end
 *             For use principly by API functions
 *
 *
 * Arg:            obj [UNKN ] Object holding the variable [AlnSequence *]
 * Arg:        bio_end [OWNER] New value of the variable [int]
 *
 * Return [SOFT ]  member variable bio_end [boolean]
 *
 */
boolean replace_bio_end_AlnSequence(AlnSequence * obj,int bio_end) 
{
    if( obj == NULL)     {  
      warn("In replacement function bio_end for object AlnSequence, got a NULL object"); 
      return FALSE;  
      }  
    obj->bio_end = bio_end;  
    return TRUE; 
}    


/* Function:  access_bio_end_AlnSequence(obj)
 *
 * Descrip:    Access member variable bio_end
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [AlnSequence *]
 *
 * Return [SOFT ]  member variable bio_end [int]
 *
 */
int access_bio_end_AlnSequence(AlnSequence * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function bio_end for object AlnSequence, got a NULL object");    
      return 0;  
      }  
    return obj->bio_end;     
}    



#ifdef _cplusplus
}
#endif
