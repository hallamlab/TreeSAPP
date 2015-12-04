#ifdef _cplusplus
extern "C" {
#endif
#include "hitlist.h"

/* Function:  sort_HitList_by_score(hl)
 *
 * Descrip:    Sorts by score
 *
 *
 * Arg:        hl [UNKN ] Undocumented argument [HitList *]
 *
 */
# line 59 "hitlist.dy"
void sort_HitList_by_score(HitList * hl)
{
  sort_HitList(hl,compare_HitPair_score);
}


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
# line 68 "hitlist.dy"
int compare_HitPair_score(HitPair * one,HitPair * two) 
{
  return two->raw_score - one->raw_score;
}


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
# line 77 "hitlist.dy"
void apply_SearchStat_to_HitList(HitList * hspm,SearchStatInterface * ssi,int database_size)
{
  int i;
  int j;

  hspm->stat_attrib = stringalloc((*ssi->attribution)(ssi->data));

  for(i=0;i<hspm->len;i++) {
    hspm->pair[i]->bit_score = (*ssi->calc_bits)(ssi->data,hspm->pair[i]->query->len,hspm->pair[i]->target->len,hspm->pair[i]->raw_score);
    hspm->pair[i]->evalue = (*ssi->calc_evalue)(ssi->data,hspm->pair[i]->query,hspm->pair[i]->target,hspm->pair[i]->raw_score,database_size);
    
    for(j=0;j<hspm->pair[i]->len;j++ ) {
      hspm->pair[i]->aln[j]->bit_score = (*ssi->calc_bits)(ssi->data,hspm->pair[i]->query->len,hspm->pair[i]->target->len,hspm->pair[i]->aln[j]->raw_score);
      hspm->pair[i]->aln[j]->evalue = (*ssi->calc_evalue)(ssi->data,hspm->pair[i]->query,hspm->pair[i]->target,hspm->pair[i]->aln[j]->raw_score,database_size); 
    }
  }

}


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
# line 100 "hitlist.dy"
HitList * HitList_from_LinearHSPmanager(LinearHSPmanager * lm)
{
  HitList * out;
  HitPair * pair;
  int i;


  out = HitList_alloc_std();
  if( lm->mat != NULL ) 
    out->mat = hard_link_CompMat(lm->mat);

  for(i=0;i<lm->len;i++) {
    pair = HitPair_from_HSPset(lm->set[i],lm->mat);
    add_HitList(out,pair);
  }
  
  return out;
}



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
# line 125 "hitlist.dy"
HitPair * HitPair_from_HSPset(HSPset * set,CompMat * mat)
{
  HitPair * out;
  HitAln * aln;
  int i;

  out = HitPair_alloc_std();
  out->query  = hard_link_Sequence(set->hsp[0]->query);
  out->target = hard_link_Sequence(set->hsp[0]->target);

  out->raw_score = 0.0;

  for(i=0;i<set->len;i++) {
    aln = HitAln_alloc();
    aln->raw_score = set->hsp[i]->score;
    aln->bit_score = aln->raw_score/2.0;
    

    aln->alb = ungapped_AlnBlock_from_HSP(set->hsp[i],out->query,out->target,mat);
    add_HitPair(out,aln);
    out->raw_score += set->hsp[i]->score;
  }


  return out;
}

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
# line 156 "hitlist.dy"
AlnBlock * ungapped_AlnBlock_from_HSP(HSP * hsp,Sequence * q,Sequence * t,CompMat * mat)
{
  AlnBlock * alb;
  AlnColumn * alc;
  AlnColumn * prev = NULL;
  int i;

  alb = AlnBlock_alloc_len(2);
  
  add_AlnBlock(alb,AlnSequence_alloc());
  add_AlnBlock(alb,AlnSequence_alloc());
  
  for(i=0;i<hsp->length && hsp->query_start+i < q->len && hsp->target_start+i < t->len;i++) {
    alc = new_pairwise_AlnColumn();

    alc->alu[0]->start = hsp->query_start+i-1;
    alc->alu[0]->end   = hsp->query_start+i;
    alc->alu[0]->text_label = "SEQUENCE";

    alc->alu[1]->start = hsp->target_start+i-1;
    alc->alu[1]->end   = hsp->target_start+i;
    alc->alu[1]->text_label = "SEQUENCE";

    if( mat != NULL ) {
      alc->alu[0]->score[0] = alc->alu[1]->score[0] = mat->comp[toupper(q->seq[hsp->query_start+i])-'A'][toupper(t->seq[hsp->target_start+i])-'A'];
    }


    if( prev == NULL ) {
      alb->start = alc;
      prev = alc;
    } else {
      prev->next = alc;
      prev = alc;
    }
  }

  return alb;
}

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
# line 199 "hitlist.dy"
HitListOutputImpl * new_HitListOutputImpl_from_argv(int * argc,char ** argv)
{
  HitListOutputImpl * out;
  char * temp;

  out = HitListOutputImpl_alloc();

  if( strip_out_boolean_argument(argc,argv,"hithelp") == TRUE ) {
    fprintf(stdout,"Hitlist help\npseudoblast gives a format similar to BLAST output\nTab gives a tab delimited format one line foreach ungapped block with columns\n<bit_score> <query-id> <query-start> <query-end> <query-strand> <query-len> <target-id> <target_start> <target_end> <target_strand> <target-len> <alignment-group-id>\n");
    fprintf(stdout,"aln gives cumlative score align label dumping, good for debugging\n");
    exit(0);
  }

  if( (temp = strip_out_assigned_argument(argc,argv,"hitoutput")) != NULL ) {
    if( strcmp(temp,"pseudoblast") == 0 ) {
      out->type = HitListOutputFormatPseudoBlast;
    } 
    if( strcmp(temp,"xml") == 0 ) {
      out->type = HitListOutputFormatXML;
    } 
    if( strcmp(temp,"aln") == 0 ) {
      out->type = HitListAlnCumlative;
    }
    if( strcmp(temp,"tab") == 0 ) {
      out->type = HitListOutputFormatTab;
    }
  }
  return out;
}

/* Function:  show_help_HitListOutputImpl(ofp)
 *
 * Descrip:    Shows help for HitList output
 *
 *
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 232 "hitlist.dy"
void show_help_HitListOutputImpl(FILE * ofp)
{
  fprintf(ofp,"Hit list output options\n");
  fprintf(ofp,"   -hitoutput [pseudoblast/xml/tab] pseudoblast by default\n");
  fprintf(ofp,"   -hithelp   more detailed help on hitlist formats\n");
}


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
# line 243 "hitlist.dy"
void show_HitList_HitListOutputImpl(HitListOutputImpl * hloi,HitList * hl,FILE * ofp)
{
  switch(hloi->type) {
    
  case HitListOutputFormatPseudoBlast :
    write_pseudoblast_HitList(hl,ofp);
    break;
  case HitListOutputFormatXML :
    write_XML_HitList(hl,ofp);
    break;
  case HitListOutputFormatTab :
    write_tab_HitList(hl,ofp);
    break;
  case HitListAlnCumlative :
    write_alb_HitList(hl,ofp);
    break;
  default :
    error("No valid HitListOutputFormat!");
  }

}

/* Function:  write_alb_HitList(hl,ofp)
 *
 * Descrip:    Writes Alb output
 *
 *
 * Arg:         hl [UNKN ] Undocumented argument [HitList *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 268 "hitlist.dy"
void write_alb_HitList(HitList * hl,FILE * ofp)
{
  int i,j;
  for(i=0;i<hl->len;i++) {
    fprintf(stdout,"%s %s\n",hl->pair[i]->query->name,hl->pair[i]->target->name);
    
    for(j=0;j<hl->pair[i]->len;j++) {
      fprintf(ofp,"Alignment %d\n",j);
      if( hl->pair[i]->aln[j]->alb != NULL ) 
	mapped_ascii_AlnBlock(hl->pair[i]->aln[j]->alb,Score2Bits,1,ofp);
    }
  }
}


/* Function:  write_XML_HitList(hl,ofp)
 *
 * Descrip:    Writes XML output
 *
 *
 * Arg:         hl [UNKN ] Undocumented argument [HitList *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 286 "hitlist.dy"
void write_XML_HitList(HitList * hl,FILE * ofp)
{
  int i,j;
  btCanvas * btc;
  

  fprintf(ofp,"<?xml version=\"1.0\"?>\n");
  fprintf(ofp,"<!DOCTYPE SequenceHitList PUBLIC \"-//EBI//SequenceHitList\" \"SequenceHitList.dtd\">\n");
  fprintf(ofp,"<sequencehitlist>\n");
  for(i=0;i<hl->len;i++) {
    fprintf(ofp," <sequencehit>\n");
    fprintf(ofp,"   <hitrank>%d</hitrank>\n",i);
    fprintf(ofp,"   <sequence>\n");
    fprintf(ofp,"     <id>%s</id>\n",hl->pair[i]->target->name);
    if( hl->pair[i]->target->desc != NULL ) {
      fprintf(ofp,"     <desc>%s</desc>\n",hl->pair[i]->target->desc);
    }
    fprintf(ofp,"     <residues>%s</residues>\n",hl->pair[i]->target->seq);
    fprintf(ofp,"   </sequence>\n");
    fprintf(ofp,"   <similaritymeasure>\n");
    fprintf(ofp,"     <raw_score>%d</raw_score>\n",hl->pair[i]->raw_score);
    fprintf(ofp,"     <bits_score>%.2f</bits_score>\n",hl->pair[i]->bit_score);
    fprintf(ofp,"     <evalue>%g</evalue>\n",hl->pair[i]->evalue);
    fprintf(ofp,"   </similaritymeasure>\n");
    for(j=0;j<hl->pair[i]->len;j++) {
      auto HitAln * haln = hl->pair[i]->aln[j];

      fprintf(ofp,"   <hitalignment>\n");
      fprintf(ofp,"   <similaritymeasure>\n");
      fprintf(ofp,"     <raw_score>%d</raw_score>\n",haln->raw_score);
      fprintf(ofp,"     <bits_score>%.2f</bits_score>\n",haln->bit_score);
      fprintf(ofp,"     <evalue>%g</evalue>\n",haln->evalue);
      fprintf(ofp,"   </similaritymeasure>\n");
      fprintf(ofp,"   <formatted_alignment>\n");

      btc = new_Ascii_btCanvas(stdout,20,50,7,3);  
      write_pretty_str_blast_align_btc(hl->pair[i]->aln[j]->alb,"Query:",hl->pair[i]->query->seq,"Sbjct:",hl->pair[i]->target->seq,btc);
      free_btCanvas(btc);
      fprintf(ofp,"\n");
      
      fprintf(ofp,"   </formatted_alignment>\n");
      fprintf(ofp,"   </hitalignment>\n");
    }
    fprintf(ofp,"</sequencehit>\n");
  } 
  fprintf(ofp,"</sequencehitlist>\n");
  
}

/* Function:  write_tab_HitList(hl,ofp)
 *
 * Descrip:    Writes tab delimited tab like output
 *
 *
 * Arg:         hl [UNKN ] Undocumented argument [HitList *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 338 "hitlist.dy"
void write_tab_HitList(HitList * hl,FILE * ofp)
{
  int i;
  int j;
  AlnColumn * alc;
  int qstart;
  int qend = 0;
  int tstart;
  int tend = 0;
  int strand;
  int tempt;

  for(i=0;i<hl->len;i++) {
    for(j=0;j<hl->pair[i]->len;j++) {
      for(alc = hl->pair[i]->aln[j]->alb->start;alc != NULL && strcmp(alc->alu[1]->text_label,"END") != 0 ;) {
	/* start of a block - remember start coordinates */
	qstart = alc->alu[0]->start+1+1;
	tstart = alc->alu[1]->start+1+1;
	/* progress any number of cases where the progression is the same */
	for(;alc != NULL;alc = alc->next ) {
	  if( (alc->alu[0]->end - alc->alu[0]->start) != (alc->alu[1]->end - alc->alu[1]->start) ) {
	    break;
	  }
	  /* otherwise set end points */

	  qend = alc->alu[0]->end+1;
	  tend = alc->alu[1]->end+1;
	}

	if( hl->pair[i]->target_reversed == 1 ) {
	  strand = -1;
	  tempt = tstart;
	  tstart = hl->pair[i]->target->len - tend+1;
	  tend   = hl->pair[i]->target->len - tempt+1;
	} else {
	  strand = 1;
	}

	/* end of a block. Print line */
	fprintf(ofp,"%.2f\t%s\t%d\t%d\t%d\t%d\t%s\t%d\t%d\t%d\t%d\tgroup_%d_%d\n",hl->pair[i]->bit_score,hl->pair[i]->query->name,qstart,qend,1,hl->pair[i]->query->len,hl->pair[i]->target->name,tstart,tend,strand,hl->pair[i]->target->len,i,j);
	
	/* find start of next block */
	for(;alc != NULL && strcmp(alc->alu[1]->text_label,"END") != 0;alc = alc->next ) {
	  if( alc->alu[0]->end - alc->alu[0]->start == alc->alu[1]->end - alc->alu[1]->start ) {
	    break;
	  } 
	}
	/* top loop will break at alc == NULL */
      }
    }
  }
  
}


/* Function:  write_pseudoblast_HitList(hl,ofp)
 *
 * Descrip:    Writes pseudoblast output
 *
 *
 * Arg:         hl [UNKN ] Undocumented argument [HitList *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 396 "hitlist.dy"
void write_pseudoblast_HitList(HitList * hl,FILE * ofp)
{
  int i,j;
  btCanvas * btc;
  char buffer[MAXLINE];

  fprintf(ofp,"BLASTP 2.1.2\n");
  fprintf(ofp,"\n\nReference: Wise2 Package, Ewan Birney\n");
  fprintf(ofp,"BLAST like format to play well with existing parsers. Other options are available\nSee help on the program that generated this for other options\n");
  if( hl->stat_attrib != NULL ) {
    fprintf(ofp,"  Statistics from : %s\n",hl->stat_attrib);
  }

  fprintf(ofp,"\n\n");
  fprintf(ofp,"Query= Not specified\n");
  fprintf(ofp,"\n\nSearch.................................done\n\n");
  fprintf(ofp,"                                                                   Score     E\n");
  fprintf(ofp,"Sequences producing significant alignments:                        (bits)  Value\n");

  for(i=0;i<hl->len;i++) {
    if( hl->pair[i]->target->desc != NULL ) {
      strcpy(buffer,hl->pair[i]->target->desc);	
    } else {
      strcpy(buffer," ");
    }

    buffer[50] = '\0';
    fprintf(ofp,"%15s %50s %-.2f    %.2g\n",hl->pair[i]->target->name,buffer,hl->pair[i]->bit_score,hl->pair[i]->evalue);
  }

  fprintf(ofp,"\n");

  for(i=0;i<hl->len;i++) {
    fprintf(ofp,">%s %s\n          Length = %d Reversed %d\n\n",hl->pair[i]->target->name,(hl->pair[i]->target->desc != NULL ? hl->pair[i]->target->desc : " "),hl->pair[i]->target->len,hl->pair[i]->target_reversed);
    for(j=0;j<hl->pair[i]->len;j++) {
      
      fprintf(ofp," Score = %.1f bits (%d), Expect = %g\n",
	      hl->pair[i]->aln[j]->bit_score,hl->pair[i]->aln[j]->raw_score,hl->pair[i]->aln[j]->evalue);
      fprintf(ofp,"\n");

      btc = new_Ascii_btCanvas(stdout,20,50,5,3);
      if( hl->write_btc_func == NULL ) {	  
	write_pretty_str_blast_align_btc(hl->pair[i]->aln[j]->alb,"Query:",hl->pair[i]->query->seq,"Sbjct:",hl->pair[i]->target->seq,btc);
      } else {
	(*hl->write_btc_func)(hl->pair[i]->aln[j]->alb,hl->pair[i]->query,hl->pair[i]->target,btc);
      }
      free_btCanvas(btc);
      fprintf(ofp,"\n");
    }
  }
      
}


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
# line 453 "hitlist.dy"
boolean write_pretty_Seq_blast_align_btc(AlnBlock * alb,Sequence * one,Sequence * two,btCanvas * btc)
{
  return write_pretty_str_blast_align_btc(alb,one->name,one->seq,two->name,two->seq,btc);
}

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
# line 470 "hitlist.dy"
boolean write_pretty_str_blast_align_btc(AlnBlock * alb,char * qname,char * query,char * tname,char * target,btCanvas * btc)
{
  int finished = 0;
  AlnColumn * alc;
  AlnColumn * prev = NULL;
  AlnUnit * q;
  AlnUnit * t;
  char buffer[14];

  btPasteArea * btp;

  for(alc=alb->start;alc != NULL && finished == 0;) {

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

    for(;finished == 0 && alc != NULL &&  can_get_paste_area_btCanvas(btc,1) == TRUE;prev=alc,alc=alc->next) {
      
      q = alc->alu[0];
      t = alc->alu[1];

      /*
       * at the end, break
       */
      if( strcmp(q->text_label,"END") == 0 ) {
	finished = 1;
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
	paste_char_btPasteArea(btp,0,0,toupper((int)query[q->start+1]),0);
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


    if( prev != NULL ) {
      btp = get_reserved_right_btCanvas(btc);
      
      sprintf(buffer,"%d",prev->alu[0]->end+1);
      
      paste_string_btPasteArea(btp,0,0,buffer,BC_RIGHT,0);
      
      sprintf(buffer,"%d",prev->alu[1]->end+1);
      
      paste_string_btPasteArea(btp,0,2,buffer,BC_RIGHT,0);
      
      free_btPasteArea(btp);
    }

    advance_line_btCanvas(btc);
    if( alc->next != NULL && strcmp(alc->next->alu[1]->text_label,"END") == 0 ) {
      break;
    }
  } /* end of for the alignment */

  return TRUE; /* we never returned false. Ooops! */
}
# line 646 "hitlist.c"
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
HitAln * hard_link_HitAln(HitAln * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a HitAln object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  HitAln_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [HitAln *]
 *
 */
HitAln * HitAln_alloc(void) 
{
    HitAln * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(HitAln *) ckalloc (sizeof(HitAln))) == NULL)    {  
      warn("HitAln_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->raw_score = 0;  
    out->bit_score = 0;  
    out->evalue = 0; 
    out->alb = NULL; 


    return out;  
}    


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
HitAln * free_HitAln(HitAln * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a HitAln obj. Should be trappable");    
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
    if( obj->alb != NULL)    
      free_AlnBlock(obj->alb);   


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_HitPair(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_HitPair
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [HitAln **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_HitPair(HitAln ** list,int i,int j)  
{
    HitAln * temp;   
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_HitPair(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_HitPair which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [HitAln **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_HitPair(HitAln ** list,int left,int right,int (*comp)(HitAln * ,HitAln * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_HitPair(list,left,(left+right)/2);  
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_HitPair (list,++last,i);    
      }  
    swap_HitPair (list,left,last);   
    qsort_HitPair(list,left,last-1,comp);    
    qsort_HitPair(list,last+1,right,comp);   
}    


/* Function:  sort_HitPair(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_HitPair
 *
 *
 * Arg:         obj [UNKN ] Object containing list [HitPair *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_HitPair(HitPair * obj,int (*comp)(HitAln *, HitAln *)) 
{
    qsort_HitPair(obj->aln,0,obj->len-1,comp);   
    return;  
}    


/* Function:  expand_HitPair(obj,len)
 *
 * Descrip:    Really an internal function for add_HitPair
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [HitPair *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_HitPair(HitPair * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_HitPair called with no need");    
      return TRUE;   
      }  


    if( (obj->aln = (HitAln ** ) ckrealloc (obj->aln,sizeof(HitAln *)*len)) == NULL)     {  
      warn("ckrealloc failed for expand_HitPair, returning FALSE");  
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


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
/* will expand function if necessary */ 
boolean add_HitPair(HitPair * obj,HitAln * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_HitPair(obj,obj->len + HitPairLISTLENGTH) == FALSE) 
        return FALSE;    
      }  


    obj->aln[obj->len++]=add;    
    return TRUE; 
}    


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
int flush_HitPair(HitPair * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->aln[i] != NULL)   {  
        free_HitAln(obj->aln[i]);    
        obj->aln[i] = NULL;  
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  HitPair_alloc_std(void)
 *
 * Descrip:    Equivalent to HitPair_alloc_len(HitPairLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [HitPair *]
 *
 */
HitPair * HitPair_alloc_std(void) 
{
    return HitPair_alloc_len(HitPairLISTLENGTH); 
}    


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
HitPair * HitPair_alloc_len(int len) 
{
    HitPair * out;  /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = HitPair_alloc()) == NULL)  
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->aln = (HitAln ** ) ckcalloc (len,sizeof(HitAln *))) == NULL)    {  
      warn("Warning, ckcalloc failed in HitPair_alloc_len"); 
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


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
HitPair * hard_link_HitPair(HitPair * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a HitPair object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  HitPair_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [HitPair *]
 *
 */
HitPair * HitPair_alloc(void) 
{
    HitPair * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(HitPair *) ckalloc (sizeof(HitPair))) == NULL)  {  
      warn("HitPair_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->query = NULL;   
    out->target = NULL;  
    out->raw_score = 0;  
    out->bit_score = 0;  
    out->evalue = 0; 
    out->aln = NULL; 
    out->len = out->maxlen = 0;  
    out->target_reversed = FALSE;    


    return out;  
}    


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
HitPair * free_HitPair(HitPair * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a HitPair obj. Should be trappable");   
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
    if( obj->query != NULL)  
      free_Sequence(obj->query);     
    if( obj->target != NULL) 
      free_Sequence(obj->target);    
    if( obj->aln != NULL)    {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->aln[i] != NULL) 
          free_HitAln(obj->aln[i]);  
        }  
      ckfree(obj->aln);  
      }  


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_HitList(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_HitList
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [HitPair **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_HitList(HitPair ** list,int i,int j)  
{
    HitPair * temp;  
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_HitList(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_HitList which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [HitPair **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_HitList(HitPair ** list,int left,int right,int (*comp)(HitPair * ,HitPair * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_HitList(list,left,(left+right)/2);  
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_HitList (list,++last,i);    
      }  
    swap_HitList (list,left,last);   
    qsort_HitList(list,left,last-1,comp);    
    qsort_HitList(list,last+1,right,comp);   
}    


/* Function:  sort_HitList(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_HitList
 *
 *
 * Arg:         obj [UNKN ] Object containing list [HitList *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_HitList(HitList * obj,int (*comp)(HitPair *, HitPair *)) 
{
    qsort_HitList(obj->pair,0,obj->len-1,comp);  
    return;  
}    


/* Function:  expand_HitList(obj,len)
 *
 * Descrip:    Really an internal function for add_HitList
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [HitList *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_HitList(HitList * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_HitList called with no need");    
      return TRUE;   
      }  


    if( (obj->pair = (HitPair ** ) ckrealloc (obj->pair,sizeof(HitPair *)*len)) == NULL)     {  
      warn("ckrealloc failed for expand_HitList, returning FALSE");  
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


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
/* will expand function if necessary */ 
boolean add_HitList(HitList * obj,HitPair * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_HitList(obj,obj->len + HitListLISTLENGTH) == FALSE) 
        return FALSE;    
      }  


    obj->pair[obj->len++]=add;   
    return TRUE; 
}    


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
int flush_HitList(HitList * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->pair[i] != NULL)  {  
        free_HitPair(obj->pair[i]);  
        obj->pair[i] = NULL; 
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  HitList_alloc_std(void)
 *
 * Descrip:    Equivalent to HitList_alloc_len(HitListLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [HitList *]
 *
 */
HitList * HitList_alloc_std(void) 
{
    return HitList_alloc_len(HitListLISTLENGTH); 
}    


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
HitList * HitList_alloc_len(int len) 
{
    HitList * out;  /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = HitList_alloc()) == NULL)  
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->pair = (HitPair ** ) ckcalloc (len,sizeof(HitPair *))) == NULL) {  
      warn("Warning, ckcalloc failed in HitList_alloc_len"); 
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


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
HitList * hard_link_HitList(HitList * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a HitList object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  HitList_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [HitList *]
 *
 */
HitList * HitList_alloc(void) 
{
    HitList * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(HitList *) ckalloc (sizeof(HitList))) == NULL)  {  
      warn("HitList_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->pair = NULL;    
    out->len = out->maxlen = 0;  
    out->mat = NULL; 
    out->write_btc_func = NULL;  
    out->stat_attrib = NULL; 


    return out;  
}    


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
HitList * free_HitList(HitList * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a HitList obj. Should be trappable");   
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
    if( obj->pair != NULL)   {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->pair[i] != NULL)    
          free_HitPair(obj->pair[i]);    
        }  
      ckfree(obj->pair); 
      }  
    if( obj->mat != NULL)    
      free_CompMat(obj->mat);    
    /* obj->write_btc_func is a function pointer */ 
    if( obj->stat_attrib != NULL)    
      ckfree(obj->stat_attrib);  


    ckfree(obj); 
    return NULL; 
}    


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
HitListOutputImpl * hard_link_HitListOutputImpl(HitListOutputImpl * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a HitListOutputImpl object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  HitListOutputImpl_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [HitListOutputImpl *]
 *
 */
HitListOutputImpl * HitListOutputImpl_alloc(void) 
{
    HitListOutputImpl * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(HitListOutputImpl *) ckalloc (sizeof(HitListOutputImpl))) == NULL)  {  
      warn("HitListOutputImpl_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->type = HitListOutputFormatPseudoBlast;  


    return out;  
}    


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
HitListOutputImpl * free_HitListOutputImpl(HitListOutputImpl * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a HitListOutputImpl obj. Should be trappable"); 
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



#ifdef _cplusplus
}
#endif
