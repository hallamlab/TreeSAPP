#ifdef _cplusplus
extern "C" {
#endif
#include "genedisplay.h"

/** #define DEBUG **/


/* Function:  protein2genomic_ascii_display(alb,p,gen,ct,name,main,ofp)
 *
 * Descrip:    shows the alignment in alb between protsequence and protname
 *             with genomic into ofp with pretty formatting
 *
 *
 * Arg:         alb [UNKN ] logical alignment [AlnBlock *]
 * Arg:           p [UNKN ] protein sequence [Protein *]
 * Arg:         gen [UNKN ] genomic dna to do the comparison [Genomic *]
 * Arg:          ct [UNKN ] codon table for translation [CodonTable *]
 * Arg:        name [UNKN ] length of name block [int]
 * Arg:        main [UNKN ] length of main block [int]
 * Arg:         ofp [UNKN ] output file [FILE *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 34 "genedisplay.dy"
boolean protein2genomic_ascii_display(AlnBlock * alb,Protein * p,Genomic * gen,CodonTable * ct,int name,int main,FILE * ofp)
{
  return protgene_ascii_display(alb,p->baseseq->seq,p->baseseq->name,p->baseseq->offset,gen,ct,name,main,FALSE,ofp);
}

/* Function:  protgene_ascii_display(alb,protsequence,protname,protoff,gen,ct,name,main,mult,ofp)
 *
 * Descrip:    shows the alignment in alb between protsequence and protname
 *             with genomic into ofp with pretty formatting
 *
 *
 * Arg:                 alb [UNKN ] logical alignment [AlnBlock *]
 * Arg:        protsequence [UNKN ] protein sequence - either real or an artifical consensus [char *]
 * Arg:            protname [UNKN ] name of the protein [char *]
 * Arg:             protoff [UNKN ] offset of the alb from the protein [int]
 * Arg:                 gen [UNKN ] genomic dna to do the comparison [Genomic *]
 * Arg:                  ct [UNKN ] codon table for translation [CodonTable *]
 * Arg:                name [UNKN ] length of name block [int]
 * Arg:                main [UNKN ] length of main block [int]
 * Arg:                mult [UNKN ] is multi-match [boolean]
 * Arg:                 ofp [UNKN ] output file [FILE *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 54 "genedisplay.dy"
boolean protgene_ascii_display(AlnBlock * alb,char * protsequence,char * protname,int protoff,Genomic * gen,CodonTable * ct,int name,int main,boolean mult,FILE * ofp)
{
  boolean ret;
  btCanvas * btc;

  btc = new_Ascii_btCanvas(ofp,name+6,main,0,6);
  ret = protdna_btc_display(alb,protsequence,protname,protoff,gen->baseseq,ct,name,main,btc,match_central_line_std,mult);
  free_btCanvas(btc);

  return ret;
}

/* Function:  protcdna_ascii_display(alb,protsequence,protname,protoff,cdna,ct,name,main,mult,ofp)
 *
 * Descrip:    shows the alignment in alb between protsequence and protname
 *             with cdna into ofp with pretty formatting
 *
 *
 * Arg:                 alb [UNKN ] logical alignment [AlnBlock *]
 * Arg:        protsequence [UNKN ] protein sequence - either real or an artifical consensus [char *]
 * Arg:            protname [UNKN ] name of the protein [char *]
 * Arg:             protoff [UNKN ] offset of the alb from the protein [int]
 * Arg:                cdna [UNKN ] cdna of the match [cDNA *]
 * Arg:                  ct [UNKN ] codon table for translation [CodonTable *]
 * Arg:                name [UNKN ] length of name block [int]
 * Arg:                main [UNKN ] length of main block [int]
 * Arg:                mult [UNKN ] is multi-match [boolean]
 * Arg:                 ofp [UNKN ] output file [FILE *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 81 "genedisplay.dy"
boolean protcdna_ascii_display(AlnBlock * alb,char * protsequence,char * protname,int protoff,cDNA * cdna,CodonTable * ct,int name,int main,boolean mult,FILE * ofp)
{
  boolean ret;
  btCanvas * btc;

  btc = new_Ascii_btCanvas(ofp,name+6,main,0,6);
  ret = protdna_btc_display(alb,protsequence,protname,protoff,cdna->baseseq,ct,name,main,btc,match_central_line_std,mult);
  free_btCanvas(btc);
  return ret;
}

# line 92 "genedisplay.dy"
char match_central_line_std(char hmm,int score,char seq)
{
  if( hmm == seq ) {
    return hmm;
  }

  if( score <= 0 ) {
    return ' ';
  }

  return '+';
}
    
# line 105 "genedisplay.dy"
boolean protdna_btc_display(AlnBlock * alb,char * protsequence,char * protname_in,int protoff,Sequence * dna,CodonTable * ct,int name,int main,btCanvas * btc,char (*match_central_line)(char,int,char),boolean multalign)
{
  AlnColumn * alc;
  AlnColumn * alc_temp,*alc_endscore;
  int a_phase, d_phase;
  int intron_number = 1;
  int aln_score;
  int aln_num = 1;
  btPasteArea * btp;
  char tempbuf[2];
  char protname[60];
  char dnaname[60];
  char dnatemp[4];
  char protc;
  char transc;
  boolean is_reversed = FALSE;
  boolean issplit;

  if( strlen(protname_in) > name ) {
    info("Name %s is longer than allowed name block (%d). Truncating\n",protname_in,name);
    strncpy(protname,protname_in,name);
    protname[name] = '\0';
  } else {
    strcpy(protname,protname_in);
  }

  if( strlen(dna->name) > name ) {
    info("Name %s is longer than allowed name block (%d). Truncating\n",dna->name,name);
    strncpy(dnaname,dna->name,name);
    dnaname[name] = '\0';
  } else {
    strcpy(dnaname,dna->name);
  }

  if( dna->offset > dna->end ) {
    is_reversed = TRUE;
  }
  
  for(alc=alb->start;alc != NULL;) {

    if ( strcmp(alc->alu[1]->text_label,"END") == 0 ) 
      break; /* end of alignment */


    for(;alc != NULL && is_random_AlnColumn_genewise(alc) == TRUE;alc = alc->next)
      ;

    if( alc == NULL)
      break; /* end of alignment */


      
    if( multalign == TRUE ) {
      /* get the end score */
      for(aln_score = 0,alc_endscore=alc;alc_endscore->next != NULL;alc_endscore = alc_endscore->next) {
	if( is_random_AlnColumn_genewise(alc_endscore) == TRUE)
	  break;
	aln_score += alc_endscore->alu[0]->score[0];
      }
      /*aln_score += alc_endscore->alu[0]->score[0];*/
      write_alignment_separator(btc,aln_num++,aln_score);
    }
      
      
    while( alc != NULL ) {


      write_name_start_stuff(btc,protname,protoff,dnaname,dna,name,alc);

      for(; alc != NULL;alc=alc->next ) {
	
	if( is_random_AlnColumn_genewise(alc) == TRUE ) 
	  break;

	if( strcmp(alc->alu[1]->text_label,"INSERT") == 0 ) {
	  if( can_get_paste_area_btCanvas(btc,1) == FALSE) 
	    break;  /* back to upper for, to place names and starts */
	  btp = get_paste_area_btCanvas(btc,1);
	  
	  paste_char_btPasteArea(btp,0,0,protsequence[alc->alu[0]->start+1],0);
	  paste_char_btPasteArea(btp,0,2,'-',0);
	  free_btPasteArea(btp);
	} else if ( strcmp(alc->alu[1]->text_label,"SEQUENCE_INSERTION") == 0 ||
		    strcmp(alc->alu[1]->text_label,"SEQUENCE_DELETION") == 0 ) {
	  if( can_get_paste_area_btCanvas(btc,1) == FALSE) 
	    break;  /* back to upper for, to place names and starts */
	  btp = get_paste_area_btCanvas(btc,1);
	  
	  if( strcmp(alc->alu[0]->text_label,"INSERT_STATE")== 0 ) {
	    paste_char_btPasteArea(btp,0,0,'-',0);
	  }
	  else {
	    paste_char_btPasteArea(btp,0,0,protsequence[alc->alu[0]->end],0);
	  }

	  sprintf(tempbuf,"%d",alc->alu[1]->end - alc->alu[1]->start);
	  paste_char_btPasteArea(btp,0,3,tempbuf[0],0);
	  paste_char_btPasteArea(btp,0,2,'!',0);

	  free_btPasteArea(btp);
	  
	} else if (strcmp(alc->alu[1]->text_label,"END") == 0 && strcmp(alc->alu[0]->text_label,"END") == 0) {
	  break; /* end of alignment */
	} else if ( strcmp(alc->alu[1]->text_label,"RANDOM_SEQUENCE") == 0 ) {
	  break;
	} else if( strcmp(alc->alu[1]->text_label,"CODON") == 0 ) {
	  
	  if( can_get_paste_area_btCanvas(btc,1) == FALSE) 
	    break;  /* back to upper for, to place names and starts */
	  
	  btp = get_paste_area_btCanvas(btc,1);
	  
	  if( strcmp(alc->alu[0]->text_label,"INSERT_STATE")== 0 ) {
	    write_codon_match(btp,'-',' ',alc->alu[1]->start+1,aminoacid_from_seq(ct,dna->seq+alc->alu[1]->start+1),dna->seq+alc->alu[1]->start+1);
	  } else if( strcmp(alc->alu[0]->text_label,"BEFORE_MATCH") == 0 ||
	             strcmp(alc->alu[0]->text_label,"AFTER_MATCH") == 0 ) {
	    write_codon_match(btp,'~',' ',alc->alu[1]->start+1,aminoacid_from_seq(ct,dna->seq+alc->alu[1]->start+1),dna->seq+alc->alu[1]->start+1);
	  } else {
	    write_codon_match(btp,protsequence[alc->alu[0]->end],(*match_central_line)(protsequence[alc->alu[0]->end],alc->alu[0]->score[0],aminoacid_from_seq(ct,dna->seq+alc->alu[1]->start+1)),alc->alu[1]->start+1,aminoacid_from_seq(ct,dna->seq+alc->alu[1]->start+1),dna->seq+alc->alu[1]->start+1);
	  }
	  
	  free_btPasteArea(btp);
	  
	  continue;
	} else if ( strstartcmp(alc->alu[1]->text_label,"5SS") == 0 )  {
	  
	  
	  /*
	   * intron stuff. Figure out the start and end, 
	   * then place the 5'SS Central and End.
	   * 
	   * If we can't fit in the intron, loop over 
	   * in this region before returning to higher loop. 
	   *
	   */
	  
	  if( strcmp(alc->alu[1]->text_label,"5SS_PHASE_0") == 0 ) {
	    d_phase = 0;
	  } else if ( strcmp(alc->alu[1]->text_label,"5SS_PHASE_1") == 0 ) {
	    d_phase = 1;
	  } else if ( strcmp(alc->alu[1]->text_label,"5SS_PHASE_2") == 0 ) {
	    d_phase = 2;
	  } else {
	    warn("No no no. You have a non 0,1,2 phase intron (god knows how!). Not displaying it %s",alc->alu[1]->text_label);
	    advance_line_btCanvas(btc);
	    return FALSE;
	  }
	  
	  alc_temp = alc->next;
	  
	  if( strcmp(alc_temp->alu[1]->text_label,"CENTRAL_INTRON") != 0 ) {
	    warn("Bad news. I have found a 5SS in your alignment, but it is not followed by a central intron node. Don't like it!");
	    advance_line_btCanvas(btc);
	    return FALSE;
	  }
	  
	  for(alc_temp = alc_temp->next ;alc_temp != NULL && strstartcmp(alc_temp->alu[1]->text_label,"3SS") != 0;alc_temp = alc_temp->next) 
	    ;
	  
	  if( alc_temp == NULL ) {
	    warn("Got to the end of the alignment in the middle of an intron from %s. Weird!",alc->alu[1]->text_label);
	    advance_line_btCanvas(btc);
	    return FALSE;
	  }

	  if( strcmp(alc_temp->alu[1]->text_label,"3SS_PHASE_0") == 0 ) {
	    a_phase = 0;
	  } else if ( strcmp(alc_temp->alu[1]->text_label,"3SS_PHASE_1") == 0 ) {
	    a_phase = 1;
	  } else if ( strcmp(alc_temp->alu[1]->text_label,"3SS_PHASE_2") == 0 ) {
	    a_phase = 2;
	  } else {
	    warn("No no no. You have a non 0,1,2 phase intron (god knows how!). Not displaying it %s",alc_temp->alu[1]->text_label);
	    advance_line_btCanvas(btc);
	    return FALSE;
	  }

	  /*
	   * At this point we have alc on 5SS alc_temp on 3SS.
	   *
	   * Check to see if we can place 5SS and Central intron piece
	   * on the line, if not advance.
	   *
	   */

	  if( can_get_paste_area_btCanvas(btc,d_phase+7+17) == FALSE) {
	    advance_line_btCanvas(btc);
	    
	    write_name_start_stuff(btc,protname,protoff,dnaname,dna,name,alc);
	  }
	  
	  /*** ok, if we can't get it now then we are fucked ***/
	  
	  if( can_get_paste_area_btCanvas(btc,d_phase+7+17) == FALSE) {
	    warn("You have specified a length of your main canvas too small. I need at least 23 characters long.");
	    advance_line_btCanvas(btc);
	    return FALSE;
	  }

	  btp = get_paste_area_btCanvas(btc,d_phase+7);
	  
	  /* ? split phase */
	  if( a_phase == 0 || (a_phase != d_phase ) ) {
	    protc = ' ';
	    transc = ' ';
	    dnatemp[0]= '\0';
	    issplit = FALSE; 
	  } else {

	    if( strcmp(alc_temp->alu[0]->text_label,"INSERT_STATE")== 0 ) {
	      protc = '-';
	    } else {
	      protc = protsequence[alc->alu[0]->start+1];
	    }

	    dnatemp[0] = tolower((int)dna->seq[alc->alu[1]->start+1]);
	    if( d_phase == 2) {
	      dnatemp[1] = tolower((int)dna->seq[alc->alu[1]->start+2]);
	    } else {
	      dnatemp[1] = tolower((int)dna->seq[alc_temp->alu[1]->end-1]);
	    }
	    dnatemp[2] = tolower((int)dna->seq[alc_temp->alu[1]->end]);
	    dnatemp[3] = '\0';

	    transc = aminoacid_from_seq(ct,dnatemp);
	    issplit = TRUE; 
	  }

	  write_5intron_match(btp,d_phase,7,dna->seq+alc->alu[1]->start+1);
	  free_btPasteArea(btp);
	  
	  btp = get_paste_area_btCanvas(btc,17);

	  if( is_reversed == FALSE ) 
	    write_intron_desc(btp,alc->alu[1]->start+1+d_phase+dna->offset,alc_temp->alu[1]->start+3+dna->offset,intron_number++,issplit,protc,transc,dnatemp);
	  else
	    write_intron_desc(btp,dna->offset - (alc->alu[1]->start+d_phase+1),dna->offset - (alc_temp->alu[1]->start+3),intron_number++,issplit,protc,transc,dnatemp);

	  free_btPasteArea(btp);


	  /* 
	   * written the start of the intron, now to deal with the
	   * acceptor. We need to loop here, because we might go over the
	   * line length. 
	   */
	  
	  alc = alc->next->next;  /*** move alc forward two columns ***/

	  while( alc != alc_temp ) {
	    for(; alc != alc_temp;alc = alc->next) { /** alc_temp is 3SS **/
	      if( strcmp(alc->alu[1]->text_label,"PYRIMIDINE_TRACT") == 0 ) {
		if( can_get_paste_area_btCanvas(btc,1) == FALSE ) 
		  break;
		btp = get_paste_area_btCanvas(btc,1);
		paste_char_btPasteArea(btp,0,3,dna->seq[alc->alu[1]->start+1],0);
		paste_char_btPasteArea(btp,0,4,'+',0);
		free_btPasteArea(btp);
	      } else if( strcmp(alc->alu[1]->text_label,"SPACER") == 0 ) {
		if( can_get_paste_area_btCanvas(btc,1) == FALSE ) 
		  break;
		btp = get_paste_area_btCanvas(btc,1);
		paste_char_btPasteArea(btp,0,3,dna->seq[alc->alu[1]->start+1],0);
		free_btPasteArea(btp);
	      } else {
		warn("Sorry, don't know how to print %s. Skipping...",alc->alu[1]->text_label);
	      }
	    }
	 
	    /** end for for loop **/

	    if ( alc == alc_temp ) {
	      break;
	    }
	    
	    /*** run out of space ***/
	    
	    advance_line_btCanvas(btc);

	    write_name_start_stuff(btc,protname,protoff,dnaname,dna,name,alc);
	  
	  } /** end of while still in central->3SS **/
	  
	  /*
	   * Now do 3SS 
	   *
	   */
	  
	  if( can_get_paste_area_btCanvas(btc,a_phase == 0 ? 3 : 3- a_phase + 3) == FALSE ) {
	    advance_line_btCanvas(btc);
	    write_name_start_stuff(btc,protname,protoff,dnaname,dna,name,alc);
	  }

	  if( a_phase != 0 ) {
	    btp = get_paste_area_btCanvas(btc,3 - a_phase + 3);
	    
	    write_3intron_match(btp,a_phase,3,dna->seq + alc->alu[1]->start+1);
	    
	    free_btPasteArea(btp);
	  } else {
	    btp = get_paste_area_btCanvas(btc,3);
	    write_3intron_match(btp,a_phase,3,dna->seq + alc->alu[1]->start+1);
	    free_btPasteArea(btp);
	  }
	  
	  /*
	   * Finished with intron !!!
	   */
	} else {
	  warn("Sorry, could not print the alignment %s:%s column",alc->alu[0]->text_label,alc->alu[1]->text_label);
	}

      } /*** in this loop ***/


      advance_line_btCanvas(btc);

      if( alc == NULL)
	break;

      if ( is_random_AlnColumn_genewise(alc) == TRUE) 
	break;

    } /* end of while over alignments */
  } /* end of foreach alignment */


  /*** end of print ! **/

  return TRUE;
} 

# line 437 "genedisplay.dy"
boolean write_alignment_separator(btCanvas * btc,int aln,int score)
{
  char buffer[64];
  btPasteArea * btp;

  sprintf(buffer,"Alignment %d Score %4.2f (Bits)",aln,Score2Bits(score));

  btp = get_paste_area_btCanvas(btc,strlen(buffer));

  paste_string_btPasteArea(btp,0,5,buffer,BC_RIGHT,0);

  free_btPasteArea(btp);

  advance_line_btCanvas(btc);

  return TRUE;
}

# line 455 "genedisplay.dy"
boolean write_name_start_stuff(btCanvas * btc,char * protname,int protoff,char * dnaname,Sequence * dna,int name_len,AlnColumn * alc)
{
  char buffer[64];
  btPasteArea * btp;

  btp = get_reserved_left_btCanvas(btc);
	  
  paste_string_btPasteArea(btp,0,0,protname,BC_RIGHT,0);    
  paste_string_btPasteArea(btp,0,3,dnaname,BC_RIGHT,0);
  
  sprintf(buffer,"%d",alc->alu[0]->start+1+protoff);
  paste_string_btPasteArea(btp,name_len+5-strlen(buffer),0,buffer,BC_RIGHT,0);
  
  if( dna->offset < dna->end ) 
    sprintf(buffer,"%d",alc->alu[1]->start+1+dna->offset);
  else 
    sprintf(buffer,"-%d",dna->offset - (alc->alu[1]->start+1)); 

  paste_string_btPasteArea(btp,name_len+5-strlen(buffer),3,buffer,BC_RIGHT,0);
  
  free_btPasteArea(btp);

  return TRUE;
}


# line 481 "genedisplay.dy"
boolean write_intron_desc(btPasteArea * btp,int start,int stop,int in_number,boolean is_split,char prot,char trans,char * dna)
{
  char buffer[32];

  if( is_split ) {
    sprintf(buffer,"  %c:%c[%s]  ",prot,trans,dna);
    paste_string_btPasteArea(btp,0,2,buffer,BC_RIGHT,0);
  }

  if( in_number < 1000 ) {
    sprintf(buffer,"  Intron %-3d  ",in_number);
  } else {
    sprintf(buffer,"  Intron ???  ",in_number);
  }

  paste_string_btPasteArea(btp,0,3,buffer,BC_RIGHT,0);

  if( start < 10000000 && stop < 10000000 ) 
    sprintf(buffer,"[%-7d:%7d]",start,stop);
  else 
    sprintf(buffer,"[???????:???????]");

  paste_string_btPasteArea(btp,0,4,buffer,BC_RIGHT,0);

  return TRUE;
}


# line 509 "genedisplay.dy"
boolean write_3intron_match(btPasteArea * btp,int phase,int length,char * seq)
{
  char buf[2];
  int i;
  int prl;

  prl = 3 - phase + length;
  if( phase == 0 ) {
    prl = 3;
  }


  sprintf(buf,"%d",phase);

  for(i=0;i<prl;i++) {
    paste_char_btPasteArea(btp,i,3,toupper((int)seq[i]),0);
  }


  for(i=0;i<length-2;i++) {
    paste_char_btPasteArea(btp,i,4,'-',0);
  }

  
  paste_char_btPasteArea(btp,i++,4,buf[0],0);
  paste_char_btPasteArea(btp,i,4,'>',0);


  return TRUE;
}
  

# line 541 "genedisplay.dy"
boolean write_5intron_match(btPasteArea * btp,int phase,int length,char * seq)
{
  char buf[2];
  int i;

  sprintf(buf,"%d",phase);

  for(i=0;i<phase+length;i++) {
    paste_char_btPasteArea(btp,i,3,toupper((int)seq[i]),0);
  }

  paste_char_btPasteArea(btp,phase,4,'<',0);
  paste_char_btPasteArea(btp,phase+1,4,buf[0],0);
  
  for(i=2;i<length;i++) {
    paste_char_btPasteArea(btp,i+phase,4,'-',0);
  }

  return TRUE;
}

# line 562 "genedisplay.dy"
boolean write_codon_match(btPasteArea * btp,char match_letter,char midline,int c_start,char aa,char * seq) 
{
  paste_char_btPasteArea(btp,0,0,match_letter,0);
  paste_char_btPasteArea(btp,0,1,midline,0);
  paste_char_btPasteArea(btp,0,2,aa,0);
  if( strchr("ATGCatgc",*seq) != NULL ) 
    paste_char_btPasteArea(btp,0,3,tolower((int)*seq),0);
  else if( !isalpha((int)*seq) ) {
    warn("Attempting to write a non alpha chacater as part of a dna sequence [%d]",(int)(*seq));
    paste_char_btPasteArea(btp,0,3,'?',0);
  } else 
    paste_char_btPasteArea(btp,0,3,toupper((int)*seq),0);

  seq++;

  if( strchr("ATGCatgc",*seq) != NULL ) 
    paste_char_btPasteArea(btp,0,4,tolower((int)*seq),0);
  else if( !isalpha((int)*seq) ) {
    warn("Attempting to write a non alpha chacater as part of a dna sequence [%d]",(int)(*seq));
    paste_char_btPasteArea(btp,0,3,'?',0);
  }
  else 
    paste_char_btPasteArea(btp,0,4,toupper((int)*seq),0);

  seq++;

  if( strchr("ATGCatgc",*seq) != NULL ) 
    paste_char_btPasteArea(btp,0,5,tolower((int)*seq),0);
  else if( !isalpha((int)*seq) ) {
    warn("Attempting to write a non alpha chacater as part of a dna sequence [%d]",(int)(*seq));
    paste_char_btPasteArea(btp,0,3,'?',0);
  }
  else 
    paste_char_btPasteArea(btp,0,5,toupper((int)*seq),0);

  return TRUE;
}





# line 610 "genedisplay.c"

#ifdef _cplusplus
}
#endif
