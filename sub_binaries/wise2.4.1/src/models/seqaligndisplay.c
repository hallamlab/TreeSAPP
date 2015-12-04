#ifdef _cplusplus
extern "C" {
#endif
#include "seqaligndisplay.h"

/* Function:  write_pretty_Protein_align(alb,q,t,name,main,ofp)
 *
 * Descrip:    This gives an interface into the
 *             alignment display using Protein
 *             objects
 *
 *
 * Arg:         alb [UNKN ] alignment structure [AlnBlock *]
 * Arg:           q [UNKN ] first sequence [Protein *]
 * Arg:           t [UNKN ] second sequence  [Protein *]
 * Arg:        name [UNKN ] length of the name block [int]
 * Arg:        main [UNKN ] length of the main block [int]
 * Arg:         ofp [UNKN ] output file [FILE *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 30 "seqaligndisplay.dy"
boolean write_pretty_Protein_align(AlnBlock * alb,Protein * q,Protein * t,int name,int main,FILE * ofp)
{
  if( alb == NULL || q == NULL || t == NULL ) {
    warn("NULL objects being passed into write_pretty_Protein_align");
    return FALSE;
  }
  return write_pretty_seq_align(alb,q->baseseq,t->baseseq,name,main,ofp);
}



/* Function:  write_pretty_seq_align(alb,q,t,name,main,ofp)
 *
 * Descrip:    This gives an interface into the alignment
 *             display using sequences and files. A more
 *             generic function is write_pretty_str_align
 *
 *
 * Arg:         alb [UNKN ] alignment structure [AlnBlock *]
 * Arg:           q [UNKN ] first sequence [Sequence *]
 * Arg:           t [UNKN ] second sequence  [Sequence *]
 * Arg:        name [UNKN ] length of the name block [int]
 * Arg:        main [UNKN ] length of the main block [int]
 * Arg:         ofp [UNKN ] output file [FILE *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 53 "seqaligndisplay.dy"
boolean write_pretty_seq_align(AlnBlock * alb,Sequence * q,Sequence * t,int name,int main,FILE * ofp)
{
  char qname[64];
  char tname[64];
  btCanvas * btc;

  if( alb == NULL || q == NULL || t == NULL ) {
    warn("NULL objects being passed into write_pretty_seq_align");
    return FALSE;
  }

  if( name > 64 ) {
    warn("Sorry - hard coded limited, can't have names longer than 64");
    return FALSE;
  }

  
  if(  strlen(q->name) > name ) {
    warn("Name %s is longer than allowed name block (%d). Truncating\n",q->name,name);
    strncpy(qname,q->name,name);
    qname[name] = '\0';
  } else {
    strcpy(qname,q->name);
  }

  if(  strlen(t->name) > name ) {
    warn("Name %s is longer than allowed name block (%d). Truncating\n",t->name,name);
    strncpy(tname,t->name,name);
    tname[name] = '\0';
  } else {
    strcpy(tname,t->name);
  }
  
  btc = new_Ascii_btCanvas(ofp,name+6,main,0,3); /*+6 in case we want to put in numbers */

  write_pretty_str_align_btc(alb,qname,q->seq,tname,t->seq,btc);

  /** destroy btc canvas **/

  free_btCanvas(btc);

  return TRUE;
}


/* Function:  write_pretty_str_align(alb,qname,query,tname,target,name,main,ofp)
 *
 * Descrip:    This gives an interface into the alignment
 *             display using strings and files.
 *
 *
 * Arg:           alb [UNKN ] alignment structure [AlnBlock *]
 * Arg:         qname [UNKN ] name of first sequence [char *]
 * Arg:         query [UNKN ] first sequence [char *]
 * Arg:         tname [UNKN ] name of second sequence [char *]
 * Arg:        target [UNKN ] second sequence [char *]
 * Arg:          name [UNKN ] length of the name block [int]
 * Arg:          main [UNKN ] length of the main block [int]
 * Arg:           ofp [UNKN ] output file [FILE *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 111 "seqaligndisplay.dy"
boolean write_pretty_str_align(AlnBlock * alb,char * qname,char * query,char * tname,char * target,int name,int main,FILE * ofp)
{
  boolean out;
  btCanvas * btc;
  
  btc = new_Ascii_btCanvas(ofp,name+6,main,0,3); /*+6 in case we want to put in numbers */

  out = write_pretty_str_align_btc(alb,qname,query,tname,target,btc);

  /** destroy btc canvas **/

  free_btCanvas(btc);

  return out;
}

/* Function:  write_pretty_str_align_btc(alb,qname,query,tname,target,btc)
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
# line 139 "seqaligndisplay.dy"
boolean write_pretty_str_align_btc(AlnBlock * alb,char * qname,char * query,char * tname,char * target,btCanvas * btc)
{
  AlnColumn * alc;
  AlnUnit * q;
  AlnUnit * t;
  char buffer[14];

  btPasteArea * btp;

  for(alc=alb->start;alc != NULL;) {

    /** put names in **/

    btp = get_reserved_left_btCanvas(btc);
    paste_string_btPasteArea(btp,0,0,qname,BC_RIGHT,0);
    paste_string_btPasteArea(btp,0,2,tname,BC_RIGHT,0);
    
    sprintf(buffer,"%d",alc->alu[0]->start+1+1);

    paste_string_btPasteArea(btp,12,0,buffer,BC_RIGHT,0);

    sprintf(buffer,"%d",alc->alu[1]->start+1+1);

    paste_string_btPasteArea(btp,12,2,buffer,BC_RIGHT,0);
    
    free_btPasteArea(btp);
    /** now loop over this block **/

    for(;alc != NULL &&  can_get_paste_area_btCanvas(btc,1) == TRUE;alc=alc->next) {
      
      q = alc->alu[0];
      t = alc->alu[1];

      /*
       * at the end, break
       */
      if( strcmp(q->text_label,"END") == 0 ) {
	alc = NULL;
	break;
      }

      if( strcmp(t->text_label,"LOOP") == 0 ) {
	advance_line_btCanvas(btc);
	for(;alc != NULL && strcmp(alc->alu[1]->text_label,"LOOP") == 0;alc = alc->next) {
	  ;
	}
	break;
      }

      /*
       * Get the paste area, length 1, depth will be 3
       */

      btp = get_paste_area_btCanvas(btc,1);

      /*
       * Write in the query sequence
       *
       */

      if( strcmp(q->text_label,"SEQUENCE") == 0 || strstr(q->text_label,"BOUND") != NULL ) {
	paste_char_btPasteArea(btp,0,0,((int)query[q->start+1]),0);
      } else if( strcmp(q->text_label,"UNMATCHED_SEQUENCE") == 0 ) {
	paste_char_btPasteArea(btp,0,0,tolower((int)query[q->start+1]),0);
      } else {
	/** is insert- we could check **/
	if( strcmp(q->text_label,"INSERT") != 0 ) {
	  warn("Got an uninterpretable label, %s",q->text_label);
	  paste_char_btPasteArea(btp,0,0,'?',0);
	} else {
	  paste_char_btPasteArea(btp,0,0,'-',0);
	}
      }

      /*
       * Write in the target sequence
       *
       */

      if( strcmp(t->text_label,"SEQUENCE") == 0 ) {
	paste_char_btPasteArea(btp,0,2,toupper((int)target[t->start+1]),0);
      } else if( strcmp(t->text_label,"UNMATCHED_SEQUENCE") == 0 ) {
	paste_char_btPasteArea(btp,0,2,tolower((int)target[t->start+1]),0);
      } else {
	/** is insert- we could check **/
	if( strcmp(t->text_label,"INSERT") != 0 ) {
	  warn("Got an uninterpretable label, %s",t->text_label);
	  paste_char_btPasteArea(btp,0,2,'?',0);
	} else {
	  paste_char_btPasteArea(btp,0,2,'-',0);
	}
      }

      /*
       * Match line
       */



      if( strcmp(q->text_label,"SEQUENCE") == 0 && strcmp(t->text_label,"SEQUENCE") == 0 ) {
	if( q->score[0] > 0 ) {
	  if( query[q->start+1] == target[t->start+1] ) {
	    paste_char_btPasteArea(btp,0,1,target[t->start+1],0);
	  } else {	   
	    paste_char_btPasteArea(btp,0,1,'+',0);
	  }
	}
      } else 
	paste_char_btPasteArea(btp,0,1,' ',0);
      
      free_btPasteArea(btp);

    } /* end of for this block */

    advance_line_btCanvas(btc);
  } /* end of for the alignment */

  return TRUE; /* we never returned false. Ooops! */
}



# line 271 "seqaligndisplay.c"

#ifdef _cplusplus
}
#endif
