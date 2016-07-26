#ifdef _cplusplus
extern "C" {
#endif
#include "geneutil.h"


/* Function:  new_GenomicRegion_from_GeneWise(gen,pseudo,alb)
 *
 * Descrip:    Makes a new GenomicRegion with the genes
 *             predicted from this AlnBlock
 *
 *             Really a wrapper around the add_Genes_to_GenomicRegion_GeneWise
 *             and other functionality
 *
 *
 * Arg:           gen [UNKN ] genomic sequence to use [Genomic *]
 * Arg:        pseudo [UNKN ] If true, predicts frameshifted genes as pseudogenes [boolean]
 * Arg:           alb [UNKN ] genewise alignment to predict genes from [AlnBlock *]
 *
 * Return [UNKN ]  a newly allocated structure [GenomicRegion *]
 *
 */
# line 25 "geneutil.dy"
GenomicRegion * new_GenomicRegion_from_GeneWise(Genomic * gen,boolean pseudo,AlnBlock * alb)
{
  GenomicRegion * out;

  out = new_GenomicRegion(gen);
  
  add_Genes_to_GenomicRegion_GeneWise(out,gen->baseseq->offset,gen->baseseq->end,alb,NULL,pseudo,NULL);

  return out;
}

/* Function:  add_Genes_to_GenomicRegion_GeneWise(gr,org_start,org_end,alb,root,pseudo,make_name)
 *
 * Descrip:    Potential an Alnblock may have more
 *             than one gene due to looping models
 *
 *             This adds all the genes to gr
 *
 *
 *
 * Arg:               gr [UNKN ] genomic region to add genes to [GenomicRegion *]
 * Arg:        org_start [UNKN ] start point of the dna to which the alb was made from [int]
 * Arg:          org_end [UNKN ] end point of the dna to which the alb was made from [int]
 * Arg:              alb [UNKN ] logical label alignment [AlnBlock *]
 * Arg:             root [UNKN ] the second argument to make_name [char *]
 * Arg:           pseudo [UNKN ] If true, frameshifted genes are predicted as pseudo genes [boolean]
 * Arg:        make_name [FUNCP] a pointer to a function to actually make the name of the gene [char * (*make_name]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
# line 51 "geneutil.dy"
int add_Genes_to_GenomicRegion_GeneWise(GenomicRegion * gr,int org_start,int org_end,AlnBlock * alb,char * root,boolean pseudo,char * (*make_name)(Wise2_Genomic *,char *,int ,Wise2_Gene * ))
{
  int count = 0;
  Gene * gene;
  AlnColumn * alc;

  alc = alb->start;

  while( (gene = Gene_from_AlnColumn_GeneWise(alc,org_start,org_end,pseudo,&alc)) != NULL ) {
    if( make_name != NULL ) {
      gene->name = (*make_name)(gr->genomic,root,gr->len,gene);
    } 
    if( root != NULL ) {
      gene->seqname = stringalloc(root);
    }

    add_Gene_to_GenomicRegion(gr,gene);
    count++;
    if( alc == NULL ) 
      break;
  }

  return count;
}
 

/* Function:  Gene_from_AlnColumn_GeneWise(alc,org_start,org_end,predict_pseudo_for_frameshift,end)
 *
 * Descrip:    A wrap for making a gene structure from
 *             an AlnBlock derived from one of the genewise
 *             methods
 *
 *
 * Arg:                                  alc [UNKN ] Alignment column in an AlnBlock produced by genewise [AlnColumn *]
 * Arg:                            org_start [UNKN ] offset that the genewise alignment was to the coordinate system [int]
 * Arg:                              org_end [UNKN ] emd that the genewise alignment was to the coordinate system [int]
 * Arg:        predict_pseudo_for_frameshift [UNKN ] Undocumented argument [boolean]
 * Arg:                                  end [WRITE] pointer to a AlnColumn * to say when it has finished with this gene [AlnColumn **]
 *
 * Return [UNKN ]  Undocumented return value [Gene *]
 *
 */
# line 87 "geneutil.dy"
Gene * Gene_from_AlnColumn_GeneWise(AlnColumn * alc,int org_start,int org_end,boolean predict_pseudo_for_frameshift,AlnColumn ** end)
{
  Gene * out;
  Transcript * tr;
  Translation * ts;
  Exon * ex;
  AlnColumn * prev;
  SupportingFeature * sf;
  int score = 0; 
  int dosf = 0;
  int phase = 0;
  int frame_break = 0;

  while(alc != NULL && is_random_AlnColumn_genewise(alc) == TRUE )
    alc = alc->next;

  while (alc != NULL && strcmp(alc->alu[1]->text_label,"CODON") != 0 ) 
    alc = alc->next;

  if( alc == NULL)
    return NULL;

  out = Gene_alloc_len(1);
  tr = Transcript_alloc_std();
  add_Gene(out,tr);

  out->start = alc->alu[1]->start +1; 

  score += alc->alu[0]->score[0];

  for(;;) {

    /*    fprintf(stderr,"Top loop - alc at %s %s\n",alc->alu[1]->text_label,alc->alu[0]->text_label); */

    /*    fprintf(stderr,"2 Score is %.2f %s %d %.2f\n",Score2Bits(score),alc->alu[1]->text_label,alc->alu[0]->score[0],Score2Bits(alc->alu[0]->score[0]));*/
    ex = Exon_alloc_std();
    add_ex_Transcript(tr,ex);

    /*
     * this is always the start of an exon 
     */

    dosf = 0;
 
    if( strcmp(alc->alu[1]->text_label,"CODON") == 0 ) {      
      ex->start = alc->alu[1]->start+1 - out->start; /* coordinated in alignment coords */
      dosf = 1;
      phase = 0;
    } else if ( strcmp(alc->alu[1]->text_label,"3SS_PHASE_0") == 0  ) {
      ex->start = alc->alu[1]->start +4 - out->start; 
      dosf = 1;
      phase = 0;
    } else if ( strcmp(alc->alu[1]->text_label,"3SS_PHASE_1") == 0  ) {
      ex->start = alc->alu[1]->start +4 - out->start; 
      phase = 1;
    } else if ( strcmp(alc->alu[1]->text_label,"3SS_PHASE_2") == 0  ) {
      ex->start = alc->alu[1]->start +4 - out->start; 
      phase = 2;
    } else {
      ex->start = alc->alu[1]->start +1 - out->start; /* coordinated in alignment coords */
      phase = -1;
    }

    ex->phase = phase;

    /* 
     * Exons can start in INSERTs (yuk). In which case we don't 
     * make a supporting feature
     */


    if( dosf == 1 && strstartcmp(alc->alu[0]->text_label,"INSERT") != 0 ) {
      sf = SupportingFeature_alloc();
      /* we fill in start and end from the exon at the moment. Not pretty */
      sf->hstart  = alc->alu[0]->start+1;
      sf->hstrand = 1; /* currently only got proteins. Thank the lord! */
      sf->start = ex->start;
    } else {
      sf = NULL; /* make sure we don't generate a sf here */
    }
  


    for(prev=alc,alc=alc->next;alc != NULL; ) {
      /*      fprintf(stderr,"Exon loop - alc at %s %s %d\n",alc->alu[1]->text_label,alc->alu[0]->text_label,sf); */

      score += alc->alu[0]->score[0];
      /*      fprintf(stderr,"1 Score is %.2f %s\n",Score2Bits(score),alc->alu[1]->text_label); */
      if( is_frameshift_AlnColumn_genewise(alc) == TRUE && predict_pseudo_for_frameshift == TRUE ) {
	score += alc->alu[0]->score[0];
	fprintf(stderr,"Score is %.2f\n",Score2Bits(score));
	out->ispseudo = TRUE;
	alc = alc->next;
	continue;
      }

      if( is_random_AlnColumn_genewise(alc) == TRUE)
	break;
      if( strcmp(alc->alu[1]->text_label,"CODON") != 0 ) {
	if( strstartcmp(alc->alu[0]->text_label,"DELETE") == 0 ) {
	  /* must add sf and start a new one */
	  /*fprintf(stderr,"Looking at alc at %s %s %d %d %d\n",alc->alu[1]->text_label,alc->alu[0]->text_label,prev,out,sf);*/

	  if( sf != NULL ) {
	    sf->end  = prev->alu[1]->end+1 - out->start;
	    sf->hend = prev->alu[0]->end+1;
	    add_Exon(ex,sf);
	    sf = NULL;
	  }
	
	  /* 
	   * go the end of this delete run, which are residues in the query with no
	   * target info
	   */
	  while( alc->next != NULL && strstartcmp(alc->next->alu[0]->text_label,"DELETE") == 0 ) {
	    alc = alc->next;
	  }
	  if( alc->next != NULL && strcmp(alc->next->alu[1]->text_label,"CODON") == 0 ) {
	    /* the next position is the start of the new alignment */
	    sf = SupportingFeature_alloc();
	    sf->hstart = alc->next->alu[0]->start+1;
	    sf->start  = alc->next->alu[1]->start+1 - out->start;
	  } else {
	    sf = NULL;
	  }
	} else {
	  break;
	}
      
      } else {
	/* it is a codon match, but it could be an insert */
	if( strstartcmp(alc->alu[0]->text_label,"INSERT") == 0) {
	  /* break at this point, add this supporting feature */
	  if( sf != NULL ) {
	    sf->end  = prev->alu[1]->end+1 - out->start;
	    sf->hend = prev->alu[0]->end+1;
	    add_Exon(ex,sf);
	    sf = NULL;
	  }
	  
	  frame_break = 0;

	  /* go to the end of this insert run, watching for frameshifts */
	  while( alc->next != NULL && strstartcmp(alc->next->alu[0]->text_label,"INSERT") == 0 ) {

	    if( is_frameshift_AlnColumn_genewise(alc->next) == TRUE || is_random_AlnColumn_genewise(alc->next) == TRUE) {

	      if( is_frameshift_AlnColumn_genewise(alc->next) == TRUE && predict_pseudo_for_frameshift == TRUE ) {
		out->ispseudo = TRUE;
		alc = alc->next;
		continue;
	      } else {
		alc = alc->next;
		frame_break = 1;
		break;
	      }
	    }
	    alc = alc->next;
	  }
	  if( frame_break == 1 ) {
	    break; /* out of this gene */
	  }

	  if( alc->next != NULL && strcmp(alc->next->alu[1]->text_label,"CODON") == 0 ) {
	    /* the next position is the start of the new alignment */
	    sf = SupportingFeature_alloc();
	    /* do not understand why not having a +1 here is correct. Hmph */
	    sf->hstart = alc->next->alu[0]->start+1;
	    sf->start  = alc->next->alu[1]->start+1 - out->start;
	  } else {
	    sf = NULL;
	  }
	
	} else {
	  /* could be the start of a run from INSERT into match */
	  if( sf == NULL ) {
	    sf = SupportingFeature_alloc();
	    /* we fill in start and end from the exon at the moment. Not pretty */
	    sf->hstart  = alc->alu[0]->start+1;
	    sf->hstrand = 1; /* currently only got proteins. Thank the lord! */
	    sf->start = alc->alu[1]->start+1 - out->start;
	  }
	}
      } /* end of else it is a codon match */
      prev  = alc;
      alc = alc->next;
      
    }

    /*
     * The exon has ended. But why?
     */

    if( sf != NULL ) {
      sf->hend = prev->alu[0]->end+1;
      add_Exon(ex,sf);
    }

    if( alc == NULL ) {
      out->end = prev->alu[1]->end +1;
      ex->end = out->end - out->start;
      if( sf != NULL ) {
	sf->end = ex->end;
      }
      break;
    } 

    if( is_random_AlnColumn_genewise(alc) == TRUE) {
      out->end = prev->alu[1]->end +1;
      ex->end = out->end - out->start;
      if( sf != NULL ) {
	sf->end = ex->end;
      }
      break;
    }

   
   


    if( strcmp(alc->alu[1]->text_label,"5SS_PHASE_0") == 0 ) {
      out->end = alc->alu[1]->start+1;
      phase = 0;
    } else if ( strcmp(alc->alu[1]->text_label,"5SS_PHASE_1") == 0 ) {
      out->end = alc->alu[1]->start+2;
      phase = 1;
    } else if ( strcmp(alc->alu[1]->text_label,"5SS_PHASE_2") == 0 ) {
      out->end = alc->alu[1]->start+3;
      phase = 2;
    } else {
      phase = 0;
      out->end = prev->alu[1]->end +1;
      ex->end = out->end - out->start;
      if( sf != NULL ) {
	sf->end = ex->end;
      }
      break;
      
    }


    /** set end of exon to the correct size from here **/
    ex->end = out->end - out->start;


    /** sf is to the codon, not to the exon */
    if( sf != NULL ) {
      if( phase == 0 ) {
	sf->end = ex->end;
      } else if ( phase == 1 ) {
	sf->end = ex->end-1;
      } else {
	sf->end = ex->end-2;
      }
    }
      

    while( alc != NULL && strstartcmp(alc->alu[1]->text_label,"3SS") != 0 ) {
      alc = alc->next;
      score += alc->alu[0]->score[0];
    }

    if( alc == NULL ) {
      warn("Got to the end of an alignment inside an intron. Oooops!");
      break;
    }

  }/* back to for(;;) */

  
  if( end != NULL ) 
    *end = alc;

  if( org_start < org_end) {
    out->start += org_start-1;
    out->end   += org_start-1;
  } else {
    /*    fprintf(stderr,"was %d to %d\n",out->start,out->end); */
    out->start = org_start-1 - out->start;
    out->end   = org_start-1 - out->end;
    /*    fprintf(stderr,"now %d to %d (%d)\n",out->start,out->end,org_start); */
  }

  if( out->ispseudo == FALSE ) {
    ts = Translation_alloc();
    ts->start = 0;
    ts->end = length_Transcript(tr);
    ts->parent = tr;  
    add_Transcript(tr,ts);
  }

  tr->parent = out;
  
  out->bits = Score2Bits(score);
  return out;
}


/* Function:  is_frameshift_AlnColumn_genewise(alc)
 *
 * Descrip:    This function is to say what is a frameshift label
 *
 *
 * Arg:        alc [UNKN ] Undocumented argument [const AlnColumn *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 389 "geneutil.dy"
boolean is_frameshift_AlnColumn_genewise(const AlnColumn * alc)
{
  if( strcmp(alc->alu[1]->text_label,"SEQUENCE_INSERTION") == 0 ) {
    return TRUE;
  }
  if( strcmp(alc->alu[1]->text_label,"SEQUENCE_DELETION") == 0 ) {
    return TRUE;
  }
  return FALSE;
}

/* Function:  is_random_AlnColumn_genewise(alc)
 *
 * Descrip:    This function is to say where this should
 *             be skipped in alignment/gene prediction problems
 *
 *
 * Arg:        alc [UNKN ] Undocumented argument [const AlnColumn *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 405 "geneutil.dy"
boolean is_random_AlnColumn_genewise(const AlnColumn * alc)
{
  char * la;

  la  = alc->alu[1]->text_label;

  if( strcmp(la,"RANDOM_SEQUENCE") == 0 )
    return TRUE;
  if( strcmp(la,"END") == 0 )
    return TRUE;

  la  = alc->alu[0]->text_label;
  
  if( strstr(la,"_RND_") != NULL )
    return TRUE;

  return FALSE;
}


/* Function:  is_intron_AlnColumn_genewise(alc)
 *
 * Descrip:    This function is to say where things are introns
 *
 *
 * Arg:        alc [UNKN ] Undocumented argument [const AlnColumn *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 429 "geneutil.dy"
boolean is_intron_AlnColumn_genewise(const AlnColumn * alc)
{
  char * la;

  la  = alc->alu[1]->text_label;

  if( strcmp(la,"CENTRAL_INTRON") == 0 )
    return TRUE;
  if( strcmp(la,"PYRIMIDINE_TRACT") == 0 )
    return TRUE;
  if( strcmp(la,"SPACER") == 0 )
    return TRUE;

  return FALSE;
}


# line 466 "geneutil.c"

#ifdef _cplusplus
}
#endif
