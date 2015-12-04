#ifdef _cplusplus
extern "C" {
#endif
#include "genomicregion.h"

/* Function:  new_GenomicRegion_discard_short(gr,multiexon,singleexon)
 *
 * Descrip:    Makes a new genomic region with genes
 *             greater than XXX length going through
 *
 *
 * Arg:                gr [UNKN ] Undocumented argument [GenomicRegion *]
 * Arg:         multiexon [UNKN ] Undocumented argument [int]
 * Arg:        singleexon [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [GenomicRegion *]
 *
 */
# line 75 "genomicregion.dy"
GenomicRegion * new_GenomicRegion_discard_short(GenomicRegion * gr,int multiexon,int singleexon)
{
  GenomicRegion * out;
  int i;

  assert(gr);
  assert(gr->genomic);

  out = GenomicRegion_alloc_std();
  assert(out);

  if( gr->genomic != NULL ) {
    out->genomic = hard_link_Genomic(gr->genomic);
  }

  for(i=0;i<gr->len;i++) {
    if( gr->gene[i]->len == 0 ) {
      continue;
    }

    assert(gr->gene[i]->transcript[0]);

    if( gr->gene[i]->transcript[0]->ex_len > 1 ) {
      if( length_Transcript(gr->gene[i]->transcript[0]) > multiexon ) {
	add_GenomicRegion(out,hard_link_Gene(gr->gene[i]));
      }
    } else {
      if( length_Transcript(gr->gene[i]->transcript[0]) > singleexon ) {
	add_GenomicRegion(out,hard_link_Gene(gr->gene[i]));
      }
    }

  }

  return out;

}

/* Function:  show_GenomicRegionOptions(sgro,gr,ct,dividestr,ofp)
 *
 * Descrip:    Actually shows a genomic region wrt to the options
 *
 *
 * Arg:             sgro [UNKN ] Undocumented argument [ShowGenomicRegionOptions *]
 * Arg:               gr [UNKN ] Undocumented argument [GenomicRegion *]
 * Arg:               ct [UNKN ] Undocumented argument [CodonTable *]
 * Arg:        dividestr [UNKN ] Undocumented argument [char *]
 * Arg:              ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 116 "genomicregion.dy"
void show_GenomicRegionOptions(ShowGenomicRegionOptions * sgro,GenomicRegion * gr,CodonTable * ct,char * dividestr,FILE * ofp)
{
  if( sgro->show_raw == TRUE ) {
    show_GenomicRegion(gr,ofp);
    fprintf(ofp,"%s\n",dividestr);
  }

  if( sgro->show_trans == TRUE ) {
    dump_translations_GenomicRegion(gr,ct,ofp);
    fprintf(ofp,"%s\n",dividestr);
  }

  if( sgro->show_gene_str == TRUE ) {
    show_pretty_GenomicRegion(gr,0,ofp);
    fprintf(ofp,"%s\n",dividestr);
  }

  if( sgro->show_gene_supp == TRUE ) {
    show_pretty_GenomicRegion(gr,1,ofp);
    fprintf(ofp,"%s\n",dividestr);
  }
    
}

/* Function:  show_help_ShowGenomicRegionOptions(ofp)
 *
 * Descrip:    Help for ShowGenomicRegionOptions
 *
 *
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 143 "genomicregion.dy"
void show_help_ShowGenomicRegionOptions(FILE * ofp)
{
  fprintf(ofp,"Gene Structure Output options\n");
  fprintf(ofp,"   -trans        show translation\n");
  fprintf(ofp,"   -cdna         show virtual cdna\n");
  fprintf(ofp,"   -genes        show gene structure\n");
  fprintf(ofp,"   -genesf       show gene structure with supporting evidence\n");
  fprintf(ofp,"   -gener        show raw (offset) gene structure\n");
  fprintf(ofp,"   -ace          show ace file\n");
  
}


/* Function:  new_ShowGenomicRegionOptions_from_argv(argc,argv)
 *
 * Descrip:    Makes a ShowGenomicRegionOptions from command line
 *
 *
 * Arg:        argc [UNKN ] Undocumented argument [int *]
 * Arg:        argv [UNKN ] Undocumented argument [char **]
 *
 * Return [UNKN ]  Undocumented return value [ShowGenomicRegionOptions *]
 *
 */
# line 159 "genomicregion.dy"
ShowGenomicRegionOptions * new_ShowGenomicRegionOptions_from_argv(int * argc,char ** argv)
{
  ShowGenomicRegionOptions * out;

  out = ShowGenomicRegionOptions_alloc();

  out->show_trans    = strip_out_boolean_argument(argc,argv,"trans");
  out->show_raw      = strip_out_boolean_argument(argc,argv,"gener");
  out->show_gene_str = strip_out_boolean_argument(argc,argv,"genes");
  out->show_gene_supp = strip_out_boolean_argument(argc,argv,"genesf");
  out->show_cdna     = strip_out_boolean_argument(argc,argv,"cdna");
  out->show_ace     = strip_out_boolean_argument(argc,argv,"ace");
  out->show_ace     = strip_out_boolean_argument(argc,argv,"aceh");


  return out;
}

/* Function:  read_genes_GenomicRegion(ifp)
 *
 * Descrip:    read genes output
 *
 *
 * Arg:        ifp [UNKN ] Undocumented argument [FILE *]
 *
 * Return [UNKN ]  Undocumented return value [GenomicRegion *]
 *
 */
# line 180 "genomicregion.dy"
GenomicRegion * read_genes_GenomicRegion(FILE * ifp)
{
  char buffer[MAXLINE];
  GenomicRegion * out;
  Gene * g;
  Transcript * t;
  Translation * ts;
  Exon * e;
  int is_start = 1;


  out = GenomicRegion_alloc_std();

  
  while( fgets(buffer,MAXLINE,ifp) ) {
    if( strstartcmp(buffer,"//") == 0 ) {
      break;
    }
    
    if( strstartcmp(buffer,"Gene") == 0 ) {
      if( is_start == 1 ) {
	is_start =0;
	continue; /* next gene triggered */
      }

      /* Line is Gene xxx yyy here */
      
      g = Gene_alloc_std();

      add_GenomicRegion(out,g);

      sscanf(buffer,"Gene %d %d",&g->start,&g->end);
      g->start = g->start -1;
      t = Transcript_alloc_std();
      add_Gene(g,t);
      g->parent = out;

      while( fgets(buffer,MAXLINE,ifp) != NULL ) {
	if( strstartcmp(buffer,"Gene") == 0 ) {
	  break;
	}
	if( strstartcmp(buffer,"//") == 0 ) {
	  break;
	}
	e = Exon_alloc_std();
	add_ex_Transcript(t,e);

	sscanf(buffer," Exon %d %d",&e->start,&e->end);
	e->start = e->start - g->start -1;
	e->end   = e->end   - g->start;
      }

      ts = Translation_alloc();
      ts->start = 0;
      ts->end = length_Transcript(t);
      ts->parent = t;
      add_Transcript(t,ts);
      t->parent = g;
      
      if( strstartcmp(buffer,"//") == 0 ) {
	break;
      }
    }
  }
      
  
    return out;
}


/* Function:  show_GenomicOverlapResults(gor,ofp)
 *
 * Descrip:    shows overlap resuls vaguely humanely
 *
 *
 * Arg:        gor [UNKN ] Undocumented argument [GenomicOverlapResults *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 253 "genomicregion.dy"
void show_GenomicOverlapResults(GenomicOverlapResults * gor,FILE * ofp)
{
  int i;
  
  fprintf(ofp,"%d genes overlapped\n",gor->gene_overlap);
  for(i=0;i<gor->len;i++)
    show_GenomicOverlapGene(gor->gog[i],ofp);
}

/* Function:  show_GenomicOverlapGene(gog,ofp)
 *
 * Descrip:    shows overlap genes
 *
 *
 * Arg:        gog [UNKN ] Undocumented argument [GenomicOverlapGene *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 266 "genomicregion.dy"
void show_GenomicOverlapGene(GenomicOverlapGene * gog,FILE * ofp)
{
  fprintf(ofp,"Perfect      exons %d\n",gog->exon_perfect);
  fprintf(ofp,"Truncated    exons %d\n",gog->exon_truncated);
  fprintf(ofp,"Partial      exons %d\n",gog->exon_partial);
  fprintf(ofp,"Mispredicted exons %d\n",gog->exon_mispredicted);
  fprintf(ofp,"Missed(int)  exons %d\n",gog->exon_missed_internal);
  fprintf(ofp,"Missed(ext)  exons %d\n",gog->exon_missed_external);
}

/* Function:  Genomic_overlap(query,truth)
 *
 * Descrip:    Gives the overlap of query in target. It is reported
 *             back in the GenomicOverlapResults structure
 *
 *
 * Arg:        query [UNKN ] Undocumented argument [GenomicRegion *]
 * Arg:        truth [UNKN ] Undocumented argument [GenomicRegion *]
 *
 * Return [UNKN ]  Undocumented return value [GenomicOverlapResults *]
 *
 */
# line 280 "genomicregion.dy"
GenomicOverlapResults * Genomic_overlap(GenomicRegion * query,GenomicRegion * truth)
{
  GenomicOverlapResults * out;
  int i;
  int j;

  out = GenomicOverlapResults_alloc_std();

  for(i=0;i<query->len;i++) {
    auto Gene * gene;

    gene = query->gene[i];

    if( gene->start > gene->end ) {
      /* backward strand */
      for(j=0;j<truth->len;j++) {
	if( truth->gene[j]->start < truth->gene[j]->end )
	  continue; /* in forward strand - don't care */
	if( (gene->start < truth->gene[j]->start && gene->start > truth->gene[j]->end) || (gene->end > truth->gene[j]->end && gene->end < truth->gene[j]->start) ) {
	  out->gene_overlap++;
	  add_GenomicOverlapResults(out,Gene_overlap_backward(gene,truth->gene[j]));
	}
      }
    } else {
      for(j=0;j<truth->len;j++) {
	if( truth->gene[j]->start > truth->gene[j]->end )
	  continue; /* in backward strand - don't care */
	if( (gene->start > truth->gene[j]->start && gene->start < truth->gene[j]->end) || (gene->end < truth->gene[j]->end && gene->end > truth->gene[j]->start) ) {
	  out->gene_overlap++;
	  add_GenomicOverlapResults(out,Gene_overlap_forward(gene,truth->gene[j]));
	}
      }
    }

  }

  return out;
}


/* Function:  Gene_overlap_forward(test,truth)
 *
 * Descrip:    Works out a gene overlap for two forward genes
 *
 *
 * Arg:         test [UNKN ] Undocumented argument [Gene *]
 * Arg:        truth [UNKN ] Undocumented argument [Gene *]
 *
 * Return [UNKN ]  Undocumented return value [GenomicOverlapGene *]
 *
 */
# line 323 "genomicregion.dy"
GenomicOverlapGene * Gene_overlap_forward(Gene * test,Gene * truth)
{
  int i;
  int j;
  Transcript * tr;
  Transcript * te;
  int start_te;
  int end_te;
  int count;
  

  GenomicOverlapGene * out;


  out = GenomicOverlapGene_alloc();

  out->test = hard_link_Gene(test);
  out->truth = hard_link_Gene(truth);
  
  tr = truth->transcript[0];
  te = test->transcript[0];

  for(i=0;i<te->ex_len;i++) 
    te->exon[i]->used = FALSE;

  start_te = te->exon[0]->start + test->start;
  end_te = te->exon[te->ex_len-1]->end + test->start;

  for(i=0;i<tr->ex_len;i++) {
    /** exon is to the left of start_te **/
    if( tr->exon[i]->end + truth->start < start_te ) 
      out->exon_missed_external++;
    else if( tr->exon[i]->start + truth->start > end_te)
      out->exon_missed_external++;
    else {
      for(j=0;j<te->ex_len;j++) {
	if( (tr->exon[i]->start + truth->start == te->exon[j]->start + test->start) ) {
	  if( tr->exon[i]->end + truth->start == te->exon[j]->end + test->start) {
	    out->exon_perfect++;
	    te->exon[j]->used = TRUE;
	  }
	  else {
	    out->exon_truncated++;
	    te->exon[j]->used = TRUE;
	  }
	  break;
	}
	if( tr->exon[i]->end + truth->start == te->exon[j]->end + test->start ) {
	  out->exon_truncated++;
	  te->exon[j]->used = TRUE;
	  break;
	}
	
	if( ((tr->exon[i]->start + truth->start) > te->exon[j]->start + test->start) && ((tr->exon[i]->start + truth->start) < te->exon[j]->end + test->start) ) { 
	  out->exon_partial++;
	  te->exon[j]->used = TRUE;
	  break;
	}

	if( ((tr->exon[i]->end + truth->start) > te->exon[j]->start + test->start) && ((tr->exon[i]->end + truth->start) < te->exon[j]->end + test->start) ) { 
	  out->exon_partial++;
	  te->exon[j]->used = TRUE;
	  break;

	}


	if( ((test->start + te->exon[j]->start) > (truth->start + tr->exon[i]->start)) && (test->start + te->exon[j]->end) < (truth->start + tr->exon[i]->end)) {
	  if( j == 0 && j+1 == test->len ) {
	    /* single exon prediction */
	    out->exon_truncated++;
	    te->exon[j]->used = TRUE;
	  } else
	    out->exon_partial++;
	    te->exon[j]->used = TRUE;
	  break;
	}

      }
      if( j == te->ex_len ) {
	out->exon_missed_internal++;
      }
    }
  }

  for(i=0,count=0;i<te->ex_len;i++) {
    if( te->exon[i]->used != TRUE ) {
      count++;
    }
  }

  out->exon_mispredicted = count;
  
  return out;
}

/* Function:  Gene_overlap_backward(test,truth)
 *
 * Descrip:    Works out a gene overlap for two backward genes
 *
 *
 * Arg:         test [UNKN ] Undocumented argument [Gene *]
 * Arg:        truth [UNKN ] Undocumented argument [Gene *]
 *
 * Return [UNKN ]  Undocumented return value [GenomicOverlapGene *]
 *
 */
# line 422 "genomicregion.dy"
GenomicOverlapGene * Gene_overlap_backward(Gene * test,Gene * truth)
{
  int i;
  int j;
  Transcript * tr;
  Transcript * te;
  int start_te;
  int end_te;
  int count;

  GenomicOverlapGene * out;

  out = GenomicOverlapGene_alloc();


  out->test = hard_link_Gene(test);
  out->truth = hard_link_Gene(truth);
  
  tr = truth->transcript[0];
  te = test->transcript[0];

  for(i=0;i<te->ex_len;i++) 
    te->exon[i]->used = FALSE;

  start_te = test->start - te->exon[0]->start;
  end_te = test->start - te->exon[te->ex_len-1]->end;

  for(i=0;i<tr->ex_len;i++) {
    /** exon is to the left of start_te **/

    if( truth->start - tr->exon[i]->end < start_te ) 
      out->exon_missed_external++;
    else if( truth->start - tr->exon[i]->start  > end_te)
      out->exon_missed_external++;
    else {
      for(j=0;j<te->ex_len;j++) {
	if( (truth->start - tr->exon[i]->start) == (test->start - te->exon[j]->start) ) {
	  if( (truth->start - tr->exon[i]->end) == (test->start - te->exon[j]->end) ) {
	    out->exon_perfect++;
	    te->exon[j]->used = TRUE;
	  } 
	  else {
	    out->exon_truncated++;
	    te->exon[j]->used = TRUE;
	  }
	  break;
	}
	if( (truth->start-tr->exon[i]->end) == (test->start - te->exon[j]->end)) {
	  out->exon_truncated++;
	  te->exon[j]->used = TRUE;
	  break;
	}
	
	if( ((truth->start - tr->exon[i]->start) < test->start - te->exon[j]->start) && ((truth->start-tr->exon[i]->start) > test->start - te->exon[j]->end) ) { 
	  out->exon_partial++;
	  te->exon[j]->used = TRUE;
	  break;
	}

	if( ((truth->start - tr->exon[i]->end) < (test->start - te->exon[j]->start)) && ((truth->start - tr->exon[i]->end) > test->start - te->exon[j]->end) ) { 
	  out->exon_partial++;
	  te->exon[j]->used = TRUE;
	  break;
	}

	if( ((test->start - te->exon[j]->start) < (truth->start - tr->exon[i]->start)) && (test->start - te->exon[j]->end) > (truth->start - tr->exon[i]->end)) {
	  if( j == 0 && j+1 == test->len ) {
	    /* single exon prediction */
	    out->exon_truncated++;
	    te->exon[j]->used = TRUE;
	  } else
	    out->exon_partial++;
	    te->exon[j]->used = TRUE;
	  break;
	}


      }
      if( j == te->ex_len ) {
	out->exon_missed_internal++;
      }
    }
  }

  for(i=0,count=0;i<te->ex_len;i++) {
    if( te->exon[i]->used != TRUE ) {
      count++;
    }
  }

  out->exon_mispredicted = count;

  return out;
}
      
/* Function:  simple_merged_GenomicRegion(gr,bits_cutoff,max_ext)
 *
 * Descrip:    Makes a new genomic region from the given
 *             genomic region, trying to merge close
 *             gene predictions that can be made by
 *             extending open reading frames
 *
 *
 * Arg:                 gr [UNKN ] Undocumented argument [GenomicRegion *]
 * Arg:        bits_cutoff [UNKN ] Undocumented argument [double]
 * Arg:            max_ext [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [GenomicRegion *]
 *
 */
# line 523 "genomicregion.dy"
GenomicRegion * simple_merged_GenomicRegion(GenomicRegion * gr,double bits_cutoff,int max_ext)
{
  GenomicRegion * out;
  int i,j;
  Gene * gene;

  for(i=0;i<gr->len;i++) {
    if( is_simple_prediction_Gene(gr->gene[i]) == FALSE ) {
      warn("Sorry - can only merge simple predictions of genes");
      return FALSE;
    }
  }

  sort_GenomicRegion_absolute(gr);

  out = new_GenomicRegion(gr->genomic);

  /* start with a new gene from the same position 
     as the first gene in the list */

  for(i=0;i<gr->len;i++)
    if( gr->gene[i]->bits > bits_cutoff ) 
      break;

  if( i == gr->len )
    return out;

  while( i < gr->len ) {

    /* copy this gene, and add it */
    gene = copy_Gene(gr->gene[i]);
    add_GenomicRegion(out,gene);

    /* look at the next genes over the cutoff */

    for(j=i+1;j<gr->len;j++) {
      if( gr->gene[j]->bits < bits_cutoff ) 
	continue; /* look at the next gene */

      /* if it is on the other strand - we can forget it */

      if( reversed_Gene(gene) != reversed_Gene(gr->gene[j]) )
	break; 

      /* ok - see if we can merge it */

      
    }
  }


  return out;
}	
	

/* Function:  sort_GenomicRegion_absolute(gr)
 *
 * Descrip:    sorts the genomicregion by absolute start points.
 *
 *
 * Arg:        gr [UNKN ] Undocumented argument [GenomicRegion *]
 *
 */
# line 581 "genomicregion.dy"
void sort_GenomicRegion_absolute(GenomicRegion * gr)
{
  sort_GenomicRegion(gr,compare_Gene_absolute);
}

/* Function:  compare_Gene_absolute(one,two)
 *
 * Descrip:    internal sort function
 *
 *
 * Arg:        one [UNKN ] Undocumented argument [Gene *]
 * Arg:        two [UNKN ] Undocumented argument [Gene *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
# line 589 "genomicregion.dy"
int compare_Gene_absolute(Gene * one,Gene * two)
{
  int startone;
  int starttwo;

  if( reversed_Gene(one) == TRUE ) 
    startone = one->end;
  else
    startone = one->start;


  if( reversed_Gene(two) == TRUE ) 
    starttwo = two->end;
  else
    starttwo = two->start;

  return (startone-starttwo);
}


/* Function:  get_Genomic_from_GenomicRegion(gr)
 *
 * Descrip:    gives back genomic sequence from a genomic region. This is *soft
 *             linked* - ie, dont free it and use /hard_link_Genomic if you do want to...
 *
 *
 * Arg:        gr [UNKN ] genomic region input [GenomicRegion *]
 *
 * Return [SOFT ]  a Genomic sequence [Genomic *]
 *
 */
# line 616 "genomicregion.dy"
Genomic * get_Genomic_from_GenomicRegion(GenomicRegion * gr)
{
  return gr->genomic;
}

/* Function:  read_EMBL_GenomicRegion_efetch(efetch)
 *
 * Descrip:    Reads both feature table and sequence
 *
 *
 * Arg:        efetch [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [GenomicRegion *]
 *
 */
# line 624 "genomicregion.dy"
GenomicRegion * read_EMBL_GenomicRegion_efetch(char * efetch) 
{
  FILE * ifp;
  char buffer[MAXLINE];
  GenomicRegion * out;

  sprintf(buffer,"efetch -a  %s",efetch);

  ifp = popen(buffer,"r");

  out = GenomicRegion_alloc_std();

  read_EMBL_FT_into_GenomicRegion(buffer,MAXLINE,out,ifp);

  pclose(ifp);

  sprintf(buffer,"efetch -a  -f %s",efetch);

  ifp = popen(buffer,"r");

  out->genomic = read_fasta_Genomic(ifp,-1);

  pclose(ifp);

  return out;
}


/* Function:  read_EMBL_GenomicRegion_SRS(srsquery)
 *
 * Descrip:    Reads both feature table and sequence
 *
 *
 * Arg:        srsquery [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [GenomicRegion *]
 *
 */
# line 655 "genomicregion.dy"
GenomicRegion * read_EMBL_GenomicRegion_SRS(char * srsquery) 
{
  FILE * ifp;
  char buffer[MAXLINE];
  GenomicRegion * out;

  sprintf(buffer,"getz -td '[%s]'",srsquery);

  ifp = popen(buffer,"r");

  out = GenomicRegion_alloc_std();

  read_EMBL_FT_into_GenomicRegion(buffer,MAXLINE,out,ifp);

  out->genomic = read_fasta_Genomic(ifp,-1);

  pclose(ifp);

  return out;
}

/* Function:  read_EMBL_GenomicRegion_file(filename)
 *
 * Descrip:    Reads in both EMBL sequence and features 
 *
 *
 * Arg:        filename [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [GenomicRegion *]
 *
 */
# line 679 "genomicregion.dy"
GenomicRegion * read_EMBL_GenomicRegion_file(char * filename)
{
  FILE * ifp;
  GenomicRegion * out;

  ifp = openfile(filename,"r");
  if( ifp == NULL ) {
    warn("Could not open %s for EMBL reading",filename);
    return NULL;
  }

  out = read_EMBL_GenomicRegion(ifp);

  fclose(ifp);

  return out;
}

/* Function:  read_EMBL_GenomicRegion(ifp)
 *
 * Descrip:    Reads in both EMBL sequence and features 
 *
 *
 * Arg:        ifp [UNKN ] Undocumented argument [FILE *]
 *
 * Return [UNKN ]  Undocumented return value [GenomicRegion *]
 *
 */
# line 700 "genomicregion.dy"
GenomicRegion * read_EMBL_GenomicRegion(FILE * ifp)
{
  GenomicRegion * out;
  Sequence * seq;
  char buffer[MAXLINE];
  char * name = NULL;
  char * runner;

  out = GenomicRegion_alloc_std();
  
  while( fgets(buffer,MAXLINE,ifp) != NULL ) {
    if( strstartcmp(buffer,"ID") == 0 ) {
      if( (runner=strtok(buffer+3,spacestr)) == NULL ) {
	warn("Very weird. Got an EMBL ID line with no name...");
      } else {
	name = stringalloc(runner);
      }
    } else if ( strstartcmp(buffer,"FH") == 0 ) {
      break;
    }
  }


  if( read_EMBL_FT_into_GenomicRegion(buffer,MAXLINE,out,ifp) == FALSE ) {
    warn("Could not read EMBL feature table into GenomicRegion");
    return free_GenomicRegion(out);
  }

  seq = read_Sequence_EMBL_seq(buffer,MAXLINE,ifp);

  if( name != NULL ) {
    ckfree(seq->name);
    seq->name = name;
  }

  if( seq == NULL ) {
    warn("In reading EMBL file %s, could not read sequence",name == NULL ? "Null" : name);
    return free_GenomicRegion(out);
  }

  if( (out->genomic = Genomic_from_Sequence(seq)) == NULL ) {
    warn("In reading EMBL file %s, sequence was not DNA",name == NULL ? "Null" : name);
    return free_GenomicRegion(out);
  }

  return out;
}
  

/* Function:  read_EMBL_FT_into_GenomicRegion(buffer,maxlen,gr,ifp)
 *
 * Descrip:    Reads in EMBL *features*, not sequence.
 *
 *
 * Arg:        buffer [UNKN ] Undocumented argument [char *]
 * Arg:        maxlen [UNKN ] Undocumented argument [int]
 * Arg:            gr [UNKN ] Undocumented argument [GenomicRegion *]
 * Arg:           ifp [UNKN ] Undocumented argument [FILE *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 752 "genomicregion.dy"
boolean read_EMBL_FT_into_GenomicRegion(char * buffer,int maxlen,GenomicRegion * gr,FILE * ifp)
{
  Gene * gene;

  fgets(buffer,maxlen,ifp);
 
  for(;;)  {
    if( strstartcmp(buffer,"FT") == 0 && strstr(buffer," CDS ") ) {
      gene = read_EMBL_feature_Gene(buffer,maxlen,ifp);
      if( gene != NULL ) {
	gene->parent = gr;
	add_GenomicRegion(gr,gene);
      }
    } else if ( strstartcmp(buffer,"SQ") == 0 ) {
      break;
    } else {
      if( fgets(buffer,maxlen,ifp) == NULL ) {
	break;
      }
    }
  }
  return TRUE;
}

/* Function:  show_GenomicRegion(gr,ofp)
 *
 * Descrip:    dumps genomic region in vaguely human form
 *
 *
 * Arg:         gr [UNKN ] Undocumented argument [GenomicRegion *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 779 "genomicregion.dy"
void show_GenomicRegion(GenomicRegion * gr,FILE * ofp)
{
  int i;

  for(i=0;i<gr->len;i++) {
    fprintf(ofp,"Gene %d\n",i);
    show_Gene(gr->gene[i],ofp);
    fprintf(ofp,"\n");
  }

}

/* Function:  dump_translations_GenomicRegion(gr,ct,ofp)
 *
 * Descrip:    shows all the translations in this genomic region
 *
 *
 * Arg:         gr [UNKN ] Undocumented argument [GenomicRegion *]
 * Arg:         ct [UNKN ] Undocumented argument [CodonTable *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 794 "genomicregion.dy"
void dump_translations_GenomicRegion(GenomicRegion * gr,CodonTable * ct,FILE * ofp)
{
  int i,j,k;
  Protein * trans;

  for(i=0;i<gr->len;i++) 
    for(j=0;j<gr->gene[i]->len;j++)
      for(k=0;k<gr->gene[i]->transcript[j]->len;k++) {
	trans = get_Protein_from_Translation(gr->gene[i]->transcript[j]->translation[k],ct);
 	write_fasta_Sequence(trans->baseseq,ofp);
      } 
}

/* Function:  dump_transcripts_GenomicRegion(gr,ofp)
 *
 * Descrip:    shows all the transcripts in this genomic region
 *
 *
 * Arg:         gr [UNKN ] Undocumented argument [GenomicRegion *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 810 "genomicregion.dy"
void dump_transcripts_GenomicRegion(GenomicRegion * gr,FILE * ofp)
{
  int i,j;
  cDNA * cd;

  for(i=0;i<gr->len;i++) 
    for(j=0;j<gr->gene[i]->len;j++) {
	cd = get_cDNA_from_Transcript(gr->gene[i]->transcript[j]);
 	write_fasta_Sequence(cd->baseseq,ofp);
      } 
}


/* Function:  new_GenomicRegion(gen)
 *
 * Descrip:    makes a genomicregion from a genomic sequence
 *
 *
 * Arg:        gen [UNKN ] Undocumented argument [Genomic *]
 *
 * Return [UNKN ]  Undocumented return value [GenomicRegion *]
 *
 */
# line 826 "genomicregion.dy"
GenomicRegion * new_GenomicRegion(Genomic * gen)
{
  GenomicRegion * out;

  out = GenomicRegion_alloc_std();
  out->genomic = hard_link_Genomic(gen);
  
  return out;
}

/* Function:  add_Gene_to_GenomicRegion(gr,gene)
 *
 * Descrip:    adds a Gene to this GenomicRegion, making
 *             sure that it parent/son relationship is ok
 *
 *
 * Arg:          gr [UNKN ] GenomicRegion to be added to [GenomicRegion *]
 * Arg:        gene [UNKN ] Gene to be added [Gene *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 843 "genomicregion.dy"
boolean add_Gene_to_GenomicRegion(GenomicRegion * gr,Gene * gene)
{
  gene->parent = gr;
  return add_GenomicRegion(gr,gene);
}

/* Function:  write_Embl_FT_GenomicRegion(gr,ofp)
 *
 * Descrip:    Writes Embl feature table. Does assumme that
 *             there is only one transcript per gene and only
 *             cds exons are used
 *
 *             Output like
 *
 *                FT   CDS          join(100..200)
 *
 *
 * Arg:         gr [UNKN ] Undocumented argument [GenomicRegion *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 858 "genomicregion.dy"
void write_Embl_FT_GenomicRegion(GenomicRegion * gr,FILE * ofp)
{
  int i;
  
  for(i=0;i<gr->len;i++) {
    write_Embl_FT_Gene(gr->gene[i],"CDS",ofp);
    if( gr->gene[i]->seqname == NULL) {
      fprintf(ofp,"FT                   /note=\"Wise2 gene object\"\n");
    } else {
      fprintf(ofp,"FT                   /note=\"Match to %s\"\n",gr->gene[i]->seqname);
    }
    if( gr->gene[i]->ispseudo == TRUE ) {
      fprintf(ofp,"FT                   /note=Pseudogene\n");
    }
  }
}

/* Function:  write_Diana_FT_GenomicRegion(gr,ofp)
 *
 * Descrip:    Writes Embl feature table for diana use. Does assumme that
 *             there is only one transcript per gene and only
 *             cds exons are used
 *
 *             Output like
 *
 *                FT   misc_feature       join(100..200)
 *
 *
 * Arg:         gr [UNKN ] Undocumented argument [GenomicRegion *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 884 "genomicregion.dy"
void write_Diana_FT_GenomicRegion(GenomicRegion * gr,FILE * ofp)
{
  int i;
  
  for(i=0;i<gr->len;i++) {
    write_Embl_FT_Gene(gr->gene[i],"misc_feature",ofp);
    if( gr->gene[i]->seqname == NULL) {
      fprintf(ofp,"FT                   /note=\"Wise2 gene object\"\n");
    } else {
      fprintf(ofp,"FT                   /note=\"Match to %s Score %.2f\"\n",gr->gene[i]->seqname,gr->gene[i]->bits);
    }
    if( gr->gene[i]->ispseudo == TRUE ) {
      fprintf(ofp,"FT                   /note=Pseudogene\n");
    }
  }
}
  


/* Function:  show_ace_GenomicRegion(gr,seq_name,ofp)
 *
 * Descrip:    shows ACeDB subsequence source.
 *
 *             Assummes
 *               a only one transcript per gene
 *               b only cds exons are used
 *
 *
 * Arg:              gr [UNKN ] Undocumented argument [GenomicRegion *]
 * Arg:        seq_name [UNKN ] Undocumented argument [char *]
 * Arg:             ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 910 "genomicregion.dy"
void show_ace_GenomicRegion(GenomicRegion * gr,char * seq_name,FILE * ofp)
{
  int i,j;
  char buffer[64];
  
  for(i=0;i<gr->len;i++) {
    fprintf(ofp,"Sequence %s\n",seq_name);
    if ( gr->gene[i]->name != NULL ) {
      strcpy(buffer,gr->gene[i]->name);
    } else {
      sprintf(buffer,"%s.%d",seq_name,i+1);
    }

    if( gr->gene[i]->start < gr->gene[i]->end )  
      fprintf(ofp,"subsequence %s %d %d\n\n",buffer,gr->gene[i]->start+1,gr->gene[i]->end);
    else 
      fprintf(ofp,"subsequence %s %d %d\n\n",buffer,gr->gene[i]->start+1,gr->gene[i]->end+2);
    
    fprintf(ofp,"Sequence %s\n",buffer);
    if( gr->gene[i]->ispseudo == FALSE ) 
      fprintf(ofp,"CDS\n");
    else 
      fprintf(ofp,"Pseudogene\n");
    fprintf(ofp,"Start_not_found\n");
    fprintf(ofp,"End_not_found\n");
    fprintf(ofp,"CDS_predicted_by genewise %.2f\n",gr->gene[i]->bits);
    for(j=0;j<gr->gene[i]->transcript[0]->ex_len;j++)
      fprintf(ofp,"source_Exons %d %d\n",gr->gene[i]->transcript[0]->exon[j]->start+1,gr->gene[i]->transcript[0]->exon[j]->end);
    fprintf(ofp,"\n");
  
  }
}

/* Function:  show_halfwise_GenomicRegion(gr,seq_name,method,db,doweb,weblocation,ofp)
 *
 * Descrip:    shows ACeDB subsequence source for halfwise
 *             method.
 *
 *             Assummes
 *               a only one transcript per gene
 *               b only cds exons are used
 *
 *
 * Arg:                 gr [UNKN ] Undocumented argument [GenomicRegion *]
 * Arg:           seq_name [UNKN ] Undocumented argument [char *]
 * Arg:             method [UNKN ] Undocumented argument [char *]
 * Arg:                 db [UNKN ] Undocumented argument [char *]
 * Arg:              doweb [UNKN ] Undocumented argument [boolean]
 * Arg:        weblocation [UNKN ] Undocumented argument [char *]
 * Arg:                ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 951 "genomicregion.dy"
void show_halfwise_GenomicRegion(GenomicRegion * gr,char * seq_name,char * method,char * db,boolean doweb,char * weblocation,FILE * ofp)
{
  int i,j;
  char buffer[64];
  
  for(i=0;i<gr->len;i++) {
    fprintf(ofp,"Sequence %s\n",seq_name);
    
    sprintf(buffer,"%s.%s.%d",seq_name,gr->gene[i]->seqname,i+1);

    if( gr->gene[i]->start < gr->gene[i]->end )  
      fprintf(ofp,"subsequence %s %d %d\n\n",buffer,gr->gene[i]->start+1,gr->gene[i]->end);
    else 
      fprintf(ofp,"subsequence %s %d %d\n\n",buffer,gr->gene[i]->start+1,gr->gene[i]->end+2);
    
    fprintf(ofp,"Sequence %s\n",buffer);
    fprintf(ofp,"Method %s\n",method);
    fprintf(ofp,"Database %s %s\n",db,gr->gene[i]->seqname);
    if( doweb) {
      fprintf(ofp,"Web_Location %s\n",weblocation);
    }

    if( gr->gene[i]->ispseudo == FALSE ) 
      fprintf(ofp,"CDS\n");
    else 
      fprintf(ofp,"Pseudogene\n");
    fprintf(ofp,"Start_not_found\n");
    fprintf(ofp,"End_not_found\n");
    fprintf(ofp,"CDS_predicted_by genewise %.2f\n",gr->gene[i]->bits);
    for(j=0;j<gr->gene[i]->transcript[0]->ex_len;j++)
      fprintf(ofp,"Source_Exons %d %d\n",gr->gene[i]->transcript[0]->exon[j]->start+1,gr->gene[i]->transcript[0]->exon[j]->end);
    fprintf(ofp,"\n");
  
  }
}

/* Function:  show_GFF_GenomicRegion(gr,seq_name,source,ofp)
 *
 * Descrip:    shows GFF output
 *
 *             Assummes
 *               a only cds exons are used
 *
 *
 * Arg:              gr [UNKN ] Undocumented argument [GenomicRegion *]
 * Arg:        seq_name [UNKN ] Undocumented argument [char *]
 * Arg:          source [UNKN ] Undocumented argument [char *]
 * Arg:             ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 993 "genomicregion.dy"
void show_GFF_GenomicRegion(GenomicRegion * gr,char * seq_name,char * source,FILE * ofp)
{
  int i,j,k,phase,len;
  char pname[64];
  int count;
  
  if( seq_name == NULL ) {
    seq_name = "SEQ";
  }
  if( source == NULL ) {
    source = "Wise-Generated";
  }

  count=0;
  for(i=0;i<gr->len;i++) {
    auto Gene * ge;
    count++;
    ge = gr->gene[i];

    sprintf(pname,"%s-genewise-prediction-%d",seq_name == NULL ? "seq" : seq_name,count);

    if( ge->start < ge->end ) {
      fprintf(ofp,"%s\t%s\tmatch\t%d\t%d\t%.2f\t+\t.\t%s\n",seq_name,source,ge->start+1,ge->end,ge->bits,pname);
    } else {
      fprintf(ofp,"%s\t%s\tmatch\t%d\t%d\t%.2f\t-\t.\t%s\n",seq_name,source,ge->start+1,ge->end+2,ge->bits,pname);
    }
    
    for(j=0;j<ge->len;j++) {
      auto Transcript * tr;

      phase = 0;
      len = 0;
      tr = ge->transcript[j];

      if( ge->start < ge->end ) {
	for(k=0;k<tr->ex_len;k++) {
	  fprintf(ofp,"%s\t%s\tcds\t%d\t%d\t%.2f\t+\t%d\t%s\n",seq_name,source,ge->start+tr->exon[k]->start+1,ge->start+tr->exon[k]->end,tr->exon[k]->score,phase,pname);

	  if( k < tr->ex_len-1 )
	    fprintf(ofp,"%s\t%s\tintron\t%d\t%d\t%.2f\t+\t.\t%s\n",seq_name,source,ge->start+tr->exon[k]->end+1,ge->start+tr->exon[k+1]->start,0.0,pname);

	  len = len + (tr->exon[k]->end - tr->exon[k]->start);
	  phase = len%3;
	  if( phase == 2 ) 
	    phase = 1;
	  else if( phase == 1 )
	    phase = 2;
	  /* else 0 */
	  
	}
      } else {
	/* reverse strand */
	for(k=0;k<tr->ex_len;k++) {
	  fprintf(ofp,"%s\t%s\tcds\t%d\t%d\t%.2f\t-\t%d\t%s\n",seq_name,source,(ge->start+1 - tr->exon[k]->start),ge->start - tr->exon[k]->end+2,tr->exon[k]->score,phase,pname);

	  if( k < tr->ex_len-1 )
	    fprintf(ofp,"%s\t%s\tintron\t%d\t%d\t%.2f\t-\t.\t%s\n",seq_name,source,ge->start+1 - tr->exon[k]->end,(ge->start - tr->exon[k+1]->start)+2,0.0,pname);

	  len += (tr->exon[k]->end - tr->exon[k]->start);
	  phase = len%3;
	  if( phase == 2 ) 
	    phase = 1;
	  else if( phase == 1 )
	    phase = 2;
	  /* else 0 */
	  
	}
      } /* end of else */
    } /* end of over transcripts */
  }
}

/* Function:  show_pretty_GenomicRegion(gr,show_supporting,ofp)
 *
 * Descrip: No Description
 *
 * Arg:                     gr [UNKN ] Undocumented argument [GenomicRegion *]
 * Arg:        show_supporting [UNKN ] Undocumented argument [boolean]
 * Arg:                    ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 1069 "genomicregion.dy"
void show_pretty_GenomicRegion(GenomicRegion * gr,boolean show_supporting,FILE * ofp)
{
  int i;

  for(i=0;i<gr->len;i++) {
    fprintf(ofp,"Gene %d\n",i+1);
    show_pretty_Gene(gr->gene[i],show_supporting,ofp);
  }
}



# line 1196 "genomicregion.c"
/* Function:  swap_GenomicRegion(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_GenomicRegion
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [Gene **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_GenomicRegion(Gene ** list,int i,int j)  
{
    Gene * temp; 
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_GenomicRegion(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_GenomicRegion which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [Gene **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_GenomicRegion(Gene ** list,int left,int right,int (*comp)(Gene * ,Gene * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_GenomicRegion(list,left,(left+right)/2);    
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_GenomicRegion (list,++last,i);  
      }  
    swap_GenomicRegion (list,left,last); 
    qsort_GenomicRegion(list,left,last-1,comp);  
    qsort_GenomicRegion(list,last+1,right,comp); 
}    


/* Function:  sort_GenomicRegion(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_GenomicRegion
 *
 *
 * Arg:         obj [UNKN ] Object containing list [GenomicRegion *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_GenomicRegion(GenomicRegion * obj,int (*comp)(Gene *, Gene *)) 
{
    qsort_GenomicRegion(obj->gene,0,obj->len-1,comp);    
    return;  
}    


/* Function:  expand_GenomicRegion(obj,len)
 *
 * Descrip:    Really an internal function for add_GenomicRegion
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [GenomicRegion *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_GenomicRegion(GenomicRegion * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_GenomicRegion called with no need");  
      return TRUE;   
      }  


    if( (obj->gene = (Gene ** ) ckrealloc (obj->gene,sizeof(Gene *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_GenomicRegion, returning FALSE");    
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_GenomicRegion(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [GenomicRegion *]
 * Arg:        add [OWNER] Object to add to the list [Gene *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_GenomicRegion(GenomicRegion * obj,Gene * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_GenomicRegion(obj,obj->len + GenomicRegionLISTLENGTH) == FALSE) 
        return FALSE;    
      }  


    obj->gene[obj->len++]=add;   
    return TRUE; 
}    


/* Function:  flush_GenomicRegion(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [GenomicRegion *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_GenomicRegion(GenomicRegion * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->gene[i] != NULL)  {  
        free_Gene(obj->gene[i]); 
        obj->gene[i] = NULL; 
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  GenomicRegion_alloc_std(void)
 *
 * Descrip:    Equivalent to GenomicRegion_alloc_len(GenomicRegionLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GenomicRegion *]
 *
 */
GenomicRegion * GenomicRegion_alloc_std(void) 
{
    return GenomicRegion_alloc_len(GenomicRegionLISTLENGTH); 
}    


/* Function:  GenomicRegion_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [GenomicRegion *]
 *
 */
GenomicRegion * GenomicRegion_alloc_len(int len) 
{
    GenomicRegion * out;/* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = GenomicRegion_alloc()) == NULL)    
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->gene = (Gene ** ) ckcalloc (len,sizeof(Gene *))) == NULL)   {  
      warn("Warning, ckcalloc failed in GenomicRegion_alloc_len");   
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_GenomicRegion(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GenomicRegion *]
 *
 * Return [UNKN ]  Undocumented return value [GenomicRegion *]
 *
 */
GenomicRegion * hard_link_GenomicRegion(GenomicRegion * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a GenomicRegion object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  GenomicRegion_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GenomicRegion *]
 *
 */
GenomicRegion * GenomicRegion_alloc(void) 
{
    GenomicRegion * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(GenomicRegion *) ckalloc (sizeof(GenomicRegion))) == NULL)  {  
      warn("GenomicRegion_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->gene = NULL;    
    out->len = out->maxlen = 0;  
    out->genomic = NULL; 


    return out;  
}    


/* Function:  free_GenomicRegion(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GenomicRegion *]
 *
 * Return [UNKN ]  Undocumented return value [GenomicRegion *]
 *
 */
GenomicRegion * free_GenomicRegion(GenomicRegion * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a GenomicRegion obj. Should be trappable"); 
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
    if( obj->gene != NULL)   {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->gene[i] != NULL)    
          free_Gene(obj->gene[i]);   
        }  
      ckfree(obj->gene); 
      }  
    if( obj->genomic != NULL)    
      free_Genomic(obj->genomic);    


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_GenomicOverlapGene(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GenomicOverlapGene *]
 *
 * Return [UNKN ]  Undocumented return value [GenomicOverlapGene *]
 *
 */
GenomicOverlapGene * hard_link_GenomicOverlapGene(GenomicOverlapGene * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a GenomicOverlapGene object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  GenomicOverlapGene_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GenomicOverlapGene *]
 *
 */
GenomicOverlapGene * GenomicOverlapGene_alloc(void) 
{
    GenomicOverlapGene * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(GenomicOverlapGene *) ckalloc (sizeof(GenomicOverlapGene))) == NULL)    {  
      warn("GenomicOverlapGene_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->exon_perfect = 0;   
    out->exon_truncated = 0; 
    out->exon_partial = 0;   
    out->exon_missed_internal = 0;   
    out->exon_missed_external = 0;   
    out->exon_mispredicted = 0;  
    out->truth = NULL;   
    out->test = NULL;    


    return out;  
}    


/* Function:  free_GenomicOverlapGene(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GenomicOverlapGene *]
 *
 * Return [UNKN ]  Undocumented return value [GenomicOverlapGene *]
 *
 */
GenomicOverlapGene * free_GenomicOverlapGene(GenomicOverlapGene * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a GenomicOverlapGene obj. Should be trappable");    
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
    if( obj->truth != NULL)  
      free_Gene(obj->truth);     
    if( obj->test != NULL)   
      free_Gene(obj->test);  


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_GenomicOverlapResults(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_GenomicOverlapResults
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [GenomicOverlapGene **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_GenomicOverlapResults(GenomicOverlapGene ** list,int i,int j)  
{
    GenomicOverlapGene * temp;   
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_GenomicOverlapResults(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_GenomicOverlapResults which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [GenomicOverlapGene **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_GenomicOverlapResults(GenomicOverlapGene ** list,int left,int right,int (*comp)(GenomicOverlapGene * ,GenomicOverlapGene * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_GenomicOverlapResults(list,left,(left+right)/2);    
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_GenomicOverlapResults (list,++last,i);  
      }  
    swap_GenomicOverlapResults (list,left,last); 
    qsort_GenomicOverlapResults(list,left,last-1,comp);  
    qsort_GenomicOverlapResults(list,last+1,right,comp); 
}    


/* Function:  sort_GenomicOverlapResults(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_GenomicOverlapResults
 *
 *
 * Arg:         obj [UNKN ] Object containing list [GenomicOverlapResults *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_GenomicOverlapResults(GenomicOverlapResults * obj,int (*comp)(GenomicOverlapGene *, GenomicOverlapGene *)) 
{
    qsort_GenomicOverlapResults(obj->gog,0,obj->len-1,comp); 
    return;  
}    


/* Function:  expand_GenomicOverlapResults(obj,len)
 *
 * Descrip:    Really an internal function for add_GenomicOverlapResults
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [GenomicOverlapResults *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_GenomicOverlapResults(GenomicOverlapResults * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_GenomicOverlapResults called with no need");  
      return TRUE;   
      }  


    if( (obj->gog = (GenomicOverlapGene ** ) ckrealloc (obj->gog,sizeof(GenomicOverlapGene *)*len)) == NULL)     {  
      warn("ckrealloc failed for expand_GenomicOverlapResults, returning FALSE");    
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_GenomicOverlapResults(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [GenomicOverlapResults *]
 * Arg:        add [OWNER] Object to add to the list [GenomicOverlapGene *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_GenomicOverlapResults(GenomicOverlapResults * obj,GenomicOverlapGene * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_GenomicOverlapResults(obj,obj->len + GenomicOverlapResultsLISTLENGTH) == FALSE) 
        return FALSE;    
      }  


    obj->gog[obj->len++]=add;    
    return TRUE; 
}    


/* Function:  flush_GenomicOverlapResults(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [GenomicOverlapResults *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_GenomicOverlapResults(GenomicOverlapResults * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->gog[i] != NULL)   {  
        free_GenomicOverlapGene(obj->gog[i]);    
        obj->gog[i] = NULL;  
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  GenomicOverlapResults_alloc_std(void)
 *
 * Descrip:    Equivalent to GenomicOverlapResults_alloc_len(GenomicOverlapResultsLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GenomicOverlapResults *]
 *
 */
GenomicOverlapResults * GenomicOverlapResults_alloc_std(void) 
{
    return GenomicOverlapResults_alloc_len(GenomicOverlapResultsLISTLENGTH); 
}    


/* Function:  GenomicOverlapResults_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [GenomicOverlapResults *]
 *
 */
GenomicOverlapResults * GenomicOverlapResults_alloc_len(int len) 
{
    GenomicOverlapResults * out;/* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = GenomicOverlapResults_alloc()) == NULL)    
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->gog = (GenomicOverlapGene ** ) ckcalloc (len,sizeof(GenomicOverlapGene *))) == NULL)    {  
      warn("Warning, ckcalloc failed in GenomicOverlapResults_alloc_len");   
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_GenomicOverlapResults(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GenomicOverlapResults *]
 *
 * Return [UNKN ]  Undocumented return value [GenomicOverlapResults *]
 *
 */
GenomicOverlapResults * hard_link_GenomicOverlapResults(GenomicOverlapResults * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a GenomicOverlapResults object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  GenomicOverlapResults_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GenomicOverlapResults *]
 *
 */
GenomicOverlapResults * GenomicOverlapResults_alloc(void) 
{
    GenomicOverlapResults * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(GenomicOverlapResults *) ckalloc (sizeof(GenomicOverlapResults))) == NULL)  {  
      warn("GenomicOverlapResults_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->gene_overlap = 0;   
    out->gog = NULL; 
    out->len = out->maxlen = 0;  


    return out;  
}    


/* Function:  free_GenomicOverlapResults(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GenomicOverlapResults *]
 *
 * Return [UNKN ]  Undocumented return value [GenomicOverlapResults *]
 *
 */
GenomicOverlapResults * free_GenomicOverlapResults(GenomicOverlapResults * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a GenomicOverlapResults obj. Should be trappable"); 
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
    if( obj->gog != NULL)    {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->gog[i] != NULL) 
          free_GenomicOverlapGene(obj->gog[i]);  
        }  
      ckfree(obj->gog);  
      }  


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_ShowGenomicRegionOptions(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [ShowGenomicRegionOptions *]
 *
 * Return [UNKN ]  Undocumented return value [ShowGenomicRegionOptions *]
 *
 */
ShowGenomicRegionOptions * hard_link_ShowGenomicRegionOptions(ShowGenomicRegionOptions * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a ShowGenomicRegionOptions object: passed a NULL object");    
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  ShowGenomicRegionOptions_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ShowGenomicRegionOptions *]
 *
 */
ShowGenomicRegionOptions * ShowGenomicRegionOptions_alloc(void) 
{
    ShowGenomicRegionOptions * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(ShowGenomicRegionOptions *) ckalloc (sizeof(ShowGenomicRegionOptions))) == NULL)    {  
      warn("ShowGenomicRegionOptions_alloc failed ");    
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->show_trans = 0; 
    out->show_raw = 0;   
    out->show_cdna = 0;  
    out->show_ace = 0;   
    out->show_ace_halfwise = 0;  
    out->show_GFF = 0;   
    out->show_gene_str = 0;  
    out->show_gene_supp = 0; 


    return out;  
}    


/* Function:  free_ShowGenomicRegionOptions(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [ShowGenomicRegionOptions *]
 *
 * Return [UNKN ]  Undocumented return value [ShowGenomicRegionOptions *]
 *
 */
ShowGenomicRegionOptions * free_ShowGenomicRegionOptions(ShowGenomicRegionOptions * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a ShowGenomicRegionOptions obj. Should be trappable");  
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


/* Function:  access_gene_GenomicRegion(obj,i)
 *
 * Descrip:    Access members stored in the gene list
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the list [GenomicRegion *]
 * Arg:          i [UNKN ] Position in the list [int]
 *
 * Return [SOFT ]  Element of the list [Gene *]
 *
 */
Gene * access_gene_GenomicRegion(GenomicRegion * obj,int i) 
{
    if( obj == NULL)     {  
      warn("In accessor function gene for object GenomicRegion, got a NULL object"); 
      return NULL;   
      }  
    if( obj->len <= i )  {  
      warn("In accessor function gene for object GenomicRegion, index %%d is greater than list length %%d",i,obj->len);  
      return NULL;   
      }  
    return obj->gene[i];     
}    


/* Function:  length_gene_GenomicRegion(obj)
 *
 * Descrip:    discover the length of the list
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the list [GenomicRegion *]
 *
 * Return [UNKN ]  length of the list [int]
 *
 */
int length_gene_GenomicRegion(GenomicRegion * obj) 
{
    if( obj == NULL)     {  
      warn("In length function gene for object GenomicRegion, got a NULL object");   
      return -1;     
      }  
    return obj->len;     
}    


/* Function:  replace_genomic_GenomicRegion(obj,genomic)
 *
 * Descrip:    Replace member variable genomic
 *             For use principly by API functions
 *
 *
 * Arg:            obj [UNKN ] Object holding the variable [GenomicRegion *]
 * Arg:        genomic [OWNER] New value of the variable [Genomic *]
 *
 * Return [SOFT ]  member variable genomic [boolean]
 *
 */
boolean replace_genomic_GenomicRegion(GenomicRegion * obj,Genomic * genomic) 
{
    if( obj == NULL)     {  
      warn("In replacement function genomic for object GenomicRegion, got a NULL object");   
      return FALSE;  
      }  
    obj->genomic = genomic;  
    return TRUE; 
}    


/* Function:  access_genomic_GenomicRegion(obj)
 *
 * Descrip:    Access member variable genomic
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [GenomicRegion *]
 *
 * Return [SOFT ]  member variable genomic [Genomic *]
 *
 */
Genomic * access_genomic_GenomicRegion(GenomicRegion * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function genomic for object GenomicRegion, got a NULL object");  
      return NULL;   
      }  
    return obj->genomic;     
}    



#ifdef _cplusplus
}
#endif
