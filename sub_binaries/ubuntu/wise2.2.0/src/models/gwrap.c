#ifdef _cplusplus
extern "C" {
#endif
#include "gwrap.h"

/* Function:  show_PotentialGeneList(*pgl,ofp)
 *
 * Descrip:    Mainly for debugging. Shows a potential gene list
 *
 *
 * Arg:        *pgl [UNKN ] Undocumented argument [PotentialGeneList]
 * Arg:         ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 125 "gwrap.dy"
void show_PotentialGeneList(PotentialGeneList *pgl,FILE * ofp)
{
  int i;

  for(i=0;i<pgl->len;i++)
    fprintf(ofp,"%s from %d %d\n",pgl->pg[i]->name,pgl->pg[i]->guess_start,pgl->pg[i]->guess_end);
}

/* Function:  read_PotentialGene_file(filename)
 *
 * Descrip:    reads in potential gene format having open the
 *             file
 *
 *
 * Arg:        filename [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [PotentialGene *]
 *
 */
# line 137 "gwrap.dy"
PotentialGene * read_PotentialGene_file(char * filename)
{
  PotentialGene * out = NULL;
  FILE * ifp;
  char buffer[MAXLINE];

  ifp = openfile(filename,"r");
  if( ifp == NULL ) {
    warn("Could not open %s for reading a potential gene",filename);
    return NULL;
  }

  while( fgets(buffer,MAXLINE,ifp) != NULL ) {
    if( strstartcmp(buffer,"pgene") == 0 ) {
      out = read_PotentialGene(buffer,ifp);
      break;
    }
  }

  fclose(ifp);

  return out;
}
  
/* Function:  read_PotentialGene(*line,ifp)
 *
 * Descrip:    Reads in a potential gene format
 *
 *
 * Arg:        *line [UNKN ] Undocumented argument [char]
 * Arg:          ifp [UNKN ] Undocumented argument [FILE *]
 *
 * Return [UNKN ]  Undocumented return value [PotentialGene *]
 *
 */
# line 164 "gwrap.dy"
PotentialGene * read_PotentialGene(char *line,FILE * ifp)
{
  PotentialGene * out;
  char buffer[MAXLINE];
  PotentialTranscript * pet;

  out = PotentialGene_alloc_std();
  

  while( fgets(buffer,MAXLINE,ifp) != NULL ) {
    if( strstartcmp(buffer,"end") == 0 ) 
      break;
    if( strstartcmp(buffer,"ptrans") == 0 ) {
      pet = read_PotentialTranscript(buffer,ifp);
      if( pet == NULL ) {
	warn("Could not read a potential transcript. Continuing");
      } else {
	add_PotentialGene(out,pet);
      }
    } else {
      warn("Did not understand line [%s] in PotentialGene",buffer);
    }
  }

  return out;
}
	     
      
/* Function:  read_PotentialTranscript(line,ifp)
 *
 * Descrip:    read in format ->
 *                pgene
 *                ptrans
 *                pexon start end qstart qend
 *                endptrans
 *                endpgene
 *
 *
 * Arg:        line [UNKN ] Undocumented argument [char *]
 * Arg:         ifp [UNKN ] Undocumented argument [FILE *]
 *
 * Return [UNKN ]  Undocumented return value [PotentialTranscript *]
 *
 */
# line 200 "gwrap.dy"
PotentialTranscript *  read_PotentialTranscript(char * line,FILE * ifp)
{
  int i,tempi;
  char buffer[MAXLINE];
  char temp[MAXLINE];
  PotentialTranscript * out;
  char * runner;
  PotentialExon * pex;

  out = PotentialTranscript_alloc_std();

  while( fgets(buffer,MAXLINE,ifp) != NULL ) {
    if( strstartcmp(buffer,"end") == 0 ) 
      break;
    if( strstartcmp(buffer,"pexon") == 0 ) {
      strcpy(temp,buffer);
      pex = PotentialExon_alloc();
      runner = strtok(buffer,spacestr);
      for(i=0;i<4;i++) {
	runner = strtok(NULL,spacestr);
	if( runner == NULL || is_integer_string(runner,&tempi) == FALSE ) {
	  warn("Could not parse [%s] as a Potential Exon line",temp);
	  free_PotentialExon(pex);
	  break;
	}
	switch (i) {
	case 0 : pex->tstart = tempi; break; 
	case 1 : pex->tend   = tempi; break;
	case 2 : pex->qstart = tempi; break;
	case 3 : pex->qend = tempi; break;
	default : warn("Bollocks. V.bad bug"); break;
	}
      }
      if( i == 4 ) {
	add_PotentialTranscript(out,pex);
      } else {
	free_PotentialExon(pex);
      }
      /* back to while */
    } else {
      warn("Got an unparsable line %s in PotentialTranscript",buffer);
    }
    
  } 
	  
  return out;
}
  

/* Function:  DPEnvelope_from_PotentialGene(pg)
 *
 * Descrip:    Makes a DPEnv structure from a potential gene
 *             with Exons in it
 *
 *             If there are no potential transcripts, returns NULL,
 *             which is what you probably want
 *
 *
 * Arg:        pg [UNKN ] Undocumented argument [PotentialGene *]
 *
 * Return [UNKN ]  Undocumented return value [DPEnvelope *]
 *
 */
# line 256 "gwrap.dy"
DPEnvelope * DPEnvelope_from_PotentialGene(PotentialGene * pg)
{
  int i;
  DPEnvelope * out;

  if( pg->len == 0 ) {
    return NULL;
  }
  
  out = DPEnvelope_alloc_std();

  for(i=0;i<pg->len;i++) {
    add_PotentialTranscript_to_DPEnvelope(out,pg->pet[i],pg);
  }

  return out;
}

/* Function:  add_PotentialTranscript_to_DPEnvelope(dpenv,pet,pg)
 *
 * Descrip:    Adds the potential exons in the Potential transcript to
 *             the dpenvelope.
 *
 *
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:          pet [UNKN ] Undocumented argument [PotentialTranscript *]
 * Arg:           pg [UNKN ] Undocumented argument [PotentialGene *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 278 "gwrap.dy"
boolean add_PotentialTranscript_to_DPEnvelope(DPEnvelope * dpenv,PotentialTranscript * pet,PotentialGene * pg)
{
  int i;
  DPUnit * unit;
  boolean is_rev = FALSE;

  if( pg->guess_start > pg->guess_end ) {
    is_rev = TRUE;
  }

  if( is_rev ) 
    invsort_PotentialExons_by_start(pet);
  else
    sort_PotentialExons_by_start(pet);

  /* do links to start points */

  /* by definition we start at the start point of the Potential gene... */

  unit = DPUnit_alloc();
  unit->starti = 0; /* always */
  unit->startj = 0;
  unit->height = pet->pex[0]->qstart + 1 + pg->slop_query;
  if( is_rev ) 
    unit->length = pg->guess_start - pet->pex[0]->tstart +1 + pg->slop_target;
  else 
    unit->length = pet->pex[0]->tstart - pg->guess_start +1 + pg->slop_target;
  
  add_DPEnvelope(dpenv,unit);

     
  
  for(i=0;i<pet->len;i++) {
    unit = DPUnit_alloc();
    unit->starti = pet->pex[i]->qstart - pg->slop_query;
    unit->height = pet->pex[i]->qend - pet->pex[i]->qstart + 1 + 2*pg->slop_query;

    if( is_rev ) 
      unit->startj = pg->guess_start - pet->pex[i]->tstart - pg->slop_target;
    else 
      unit->startj = pet->pex[i]->tstart - pg->guess_start - pg->slop_target;

    unit->length = abs(pet->pex[i]->tend - pet->pex[i]->tstart) + 1 + 2*pg->slop_target;

    add_DPEnvelope(dpenv,unit);
    if( i == pet->len-1 ) 
      continue;

    /* connecting unit */
    unit = DPUnit_alloc();
    unit->starti = pet->pex[i]->qend - pg->slop_query;
    unit->height = pet->pex[i+1]->qstart - pet->pex[i]->qend +1 + 2*pg->slop_query;

    if( is_rev ) 
      unit->startj =  pg->guess_start - pet->pex[i]->tend - pg->slop_target;
    else 
      unit->startj = pet->pex[i]->tend - pg->guess_start - pg->slop_target;

    unit->length = abs(pet->pex[i+1]->tstart - pet->pex[i]->tend) +1 + 2*pg->slop_target;
    add_DPEnvelope(dpenv,unit);      
  }

  /* end point */

  i--;
  unit = DPUnit_alloc();
  unit->starti = pet->pex[i]->qend - pg->slop_query;
  unit->height = pg->query_length - pet->pex[i]->qend +1 + pg->slop_query;

  if( is_rev ) 
    unit->startj =  pg->guess_start - pet->pex[i]->tend - pg->slop_target;
  else 
    unit->startj = pet->pex[i]->tend - pg->guess_start - pg->slop_target;

  unit->length = abs(pg->guess_end - pet->pex[i]->tend) +1 + pg->slop_target;

  if( is_rev ) 
    unit->length = pg->guess_start - pet->pex[0]->tstart +1 + 2*pg->slop_target;
  else 
    unit->length = pet->pex[0]->tstart - pg->guess_start +1 + 2*pg->slop_target;
  
  add_DPEnvelope(dpenv,unit);

  return TRUE;;
}

/* Function:  sort_PotentialExons_by_start(pet)
 *
 * Descrip:    Sorts by start point
 *
 *
 * Arg:        pet [UNKN ] Undocumented argument [PotentialTranscript *]
 *
 */
# line 368 "gwrap.dy"
void sort_PotentialExons_by_start(PotentialTranscript * pet)
{
  sort_PotentialTranscript(pet,compare_PotentialExons_start);
}

/* Function:  compare_PotentialExons_start(one,two)
 *
 * Descrip:    compares by start point
 *
 *
 * Arg:        one [UNKN ] Undocumented argument [PotentialExon *]
 * Arg:        two [UNKN ] Undocumented argument [PotentialExon *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
# line 377 "gwrap.dy"
int compare_PotentialExons_start(PotentialExon * one,PotentialExon * two)
{
  return one->tstart - two->tstart;
}

/* Function:  invsort_PotentialExons_by_start(pet)
 *
 * Descrip:    Sorts by start point backwards
 *
 *
 * Arg:        pet [UNKN ] Undocumented argument [PotentialTranscript *]
 *
 */
# line 386 "gwrap.dy"
void invsort_PotentialExons_by_start(PotentialTranscript * pet)
{
  sort_PotentialTranscript(pet,invcompare_PotentialExons_start);
}

/* Function:  invcompare_PotentialExons_start(one,two)
 *
 * Descrip:    compares by start point but reversed order
 *
 *
 * Arg:        one [UNKN ] Undocumented argument [PotentialExon *]
 * Arg:        two [UNKN ] Undocumented argument [PotentialExon *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
# line 395 "gwrap.dy"
int invcompare_PotentialExons_start(PotentialExon * one,PotentialExon * two)
{
  return two->tstart - one->tstart;
}

/* Function:  PotentialGeneList_from_DnaSequenceHitList(dsl,window,wing_length,min_score)
 *
 * Descrip:    This makes a PotentialGeneList from a
 *             DnaSequenceHitList, which is a module which,
 *             for example, abstracts the MSP crunch output.
 *
 *             The three parameters are:
 *                 window - what window size to consider a potential gene in
 *                 wing_length - length of wing sequences to add onto the start/end points of a hit
 *                 min_score - minimum score to trigger a potential gene.
 *
 *             The potential genes are selected as follows:
 *
 *                 foreach window
 *                        Take the best scoring segment.
 *                        if( > min_score) 
 *                             if( there_is_a_segment which start/end + wing_length overlaps + the same name)
 *                                  extend that segment
 *                             else
 *                                  make a new potential gene
 *
 *
 *
 * Arg:                dsl [UNKN ] Undocumented argument [DnaSequenceHitList *]
 * Arg:             window [UNKN ] Undocumented argument [int]
 * Arg:        wing_length [UNKN ] Undocumented argument [int]
 * Arg:          min_score [UNKN ] Undocumented argument [double]
 *
 * Return [UNKN ]  Undocumented return value [PotentialGeneList *]
 *
 */
# line 421 "gwrap.dy"
PotentialGeneList * PotentialGeneList_from_DnaSequenceHitList(DnaSequenceHitList * dsl,int window,int wing_length,double min_score)
{
  int pos;
  int win_start;
  double best_score;
  SegmentHit * best_segment;
  PotentialGeneList * out;
  PotentialGene * pg;
  int i;
  int back_st_point;


  out = PotentialGeneList_alloc_std();

  /** assumme dsl is sorted **/

  pos = 0;
  win_start = 0;

  for(;pos < dsl->forward->len;) {
    best_score = -10000000;
    best_segment = NULL;

    for(;pos < dsl->forward->len && dsl->forward->seghit[pos]->qstart < win_start + window;pos++) {
      if( dsl->forward->seghit[pos]->score > best_score ) {
	best_score = dsl->forward->seghit[pos]->score;
	best_segment = dsl->forward->seghit[pos];
      }
    }

    if( best_segment == NULL || best_score < min_score ) {
      win_start += window;
      continue; /* no hits here! */
    }

    /* check to see if we already have used this name */

    for(i=0;i<out->len;i++) {
      if( strcmp(best_segment->name,out->pg[i]->name) == 0 ) {
	if( best_segment->qstart < out->pg[i]->guess_end ) {
	  out->pg[i]->guess_end = best_segment->qend + wing_length;
	  break;
	}
      }
    }

    if( i != out->len ) {
      /** means we have placed this into another hit **/
      win_start += window;
      continue;
    }


    /** ok - we have a new potential gene **/
    
    pg = PotentialGene_alloc();
    pg->name = stringalloc(best_segment->name);
    pg->guess_start = best_segment->qstart - wing_length;
    pg->guess_end   = best_segment->qend   + wing_length;
    
    add_PotentialGeneList(out,pg);
    win_start += window;

  }

  /** backward strands **/
  win_start = 0;
  back_st_point = out->len;

  for(pos=0;pos < dsl->backward->len;) {
    best_score = -10000;
    best_segment = NULL;

    for(;pos < dsl->backward->len && dsl->backward->seghit[pos]->qend < win_start + window;pos++) {
      if( dsl->backward->seghit[pos]->score > best_score ) {
	best_score = dsl->backward->seghit[pos]->score;
	best_segment = dsl->backward->seghit[pos];
      }
    }


    if( best_segment == NULL || best_score < min_score ) {
      win_start += window;
      continue; /* no hits here! */
    }


    /* check to see if we already have used this name Only need to check backwards ones*/

    for(i=back_st_point;i<out->len;i++) {
      if( strcmp(best_segment->name,out->pg[i]->name) == 0 ) {
	if( best_segment->qstart < out->pg[i]->guess_start ) {
	  out->pg[i]->guess_start = best_segment->qstart + wing_length;
	  break;
	}
      }
    }


    if( i != out->len ) {
      /** means we have placed this into another hit **/
      win_start += window;
      continue;
    }


    /** ok - we have a new potential gene **/
    
    pg = PotentialGene_alloc();
    pg->name = stringalloc(best_segment->name);

    /** it is reversed, hence the flip from end to start **/
    pg->guess_end     = best_segment->qend     - wing_length;
    pg->guess_start   = best_segment->qstart   + wing_length;
    
    add_PotentialGeneList(out,pg);
    win_start += window;

  }

  return out;
}

/* Function:  read_PotentialGeneList_pgfasta_file(filename)
 *
 * Descrip:    reads a file with lines like
 *
 *              >ROA1_HUMAN:1234:3445
 *              <sequence>
 *
 *             As potential gene with a guess start point at 1234 and end
 *             at 3445
 *
 *
 * Arg:        filename [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [PotentialGeneList *]
 *
 */
# line 553 "gwrap.dy"
PotentialGeneList * read_PotentialGeneList_pgfasta_file(char * filename)
{
  FILE * ifp;
  PotentialGeneList * out;

  ifp = openfile(filename,"r");
  if( ifp == NULL ) {
    warn("Could not open %s for reading Potential Gene List",filename);
    return NULL;
  }

  out = read_PotentialGeneList_pgfasta(ifp);
  fclose(ifp);

  return out;
}


/* Function:  read_PotentialGeneList_pgfasta(ifp)
 *
 * Descrip:    reads a file with lines like
 *
 *             >ROA1_HUMAN:1234:3445
 *             <sequence>
 *
 *             As potential gene with a guess start point at 1234 and end
 *             at 3445
 *
 *
 * Arg:        ifp [UNKN ] Undocumented argument [FILE *]
 *
 * Return [UNKN ]  Undocumented return value [PotentialGeneList *]
 *
 */
# line 580 "gwrap.dy"
PotentialGeneList * read_PotentialGeneList_pgfasta(FILE * ifp)
{
  Sequence * seq;
  PotentialGeneList * pgl;
  PotentialGene * pg;
  char * runner,*run2;

  pgl = PotentialGeneList_alloc_std();

  while( (seq = read_fasta_Sequence(ifp)) != NULL ) {
    pg = PotentialGene_alloc();
    pg->guess_start = -1;
    pg->guess_end   = -1;
    add_PotentialGeneList(pgl,pg);
    pg->homolog = Protein_from_Sequence(seq);
    if ( (runner=strchr(seq->name,':')) != NULL ) {
      *runner = '\0';
      ++runner;
      if( (run2 = strchr(runner,':')) != NULL ) {
	*run2 = '\0';
	++run2;
	is_integer_string(runner,&pg->guess_start);
	is_integer_string(run2,&pg->guess_end);
      }
    }
  } 

  return pgl;
}

    
      

   
/* Function:  resolve_PotentialGenes_on_GenomicRegion(gr,pgl,alg_protein,alg_hmm,prot_thr,hmm_thr,comp,gap,ext,gpara,rmd,inter,fetch_from_pipe,should_free,make_name,bit_cut_off,dpri)
 *
 * Descrip:    Takes the potential gene list, (made for example from
 *             MSP crunch file through DnaSequenceHitList object),
 *             and calculates genes on it.
 *
 *             This is the core functionality to postwise.
 *
 *             There are two basic modes to this routine that probably
 *             should be split up. In the first mode, the Potential Gene List
 *             has actual proteins/or HMMs (tsm - threestatemodels) in 
 *             memory in the structure. This then loops through the list
 *             calling genewise functions (going through the wrap functions
 *             in this module).
 *
 *             The other mode has no proteins in memory but a way of fetching
 *             proteins from a pipe using a string passed in. The string looks
 *             like a sprintf string where the name of the protein is the
 *             only thing to be substituted. For example
 *                 
 *                efetch -f %s
 *
 *             or
 *
 *                getz -d '[swissprot-id:%s]'
 *
 *             would be sensible examples. Remember that the command should produce
 *             things on stdout, and that you should (obviously) make sure that the
 *             indexer uses the same name as the name in, for example, the msp crunch
 *             file.
 *
 *             If should_free is true then it frees the protein sequence and any
 *             alignment after the calculation. Otherwise it stores both of these
 *             in the potentialgene structure.
 *
 *             For the business end of the algorithm, this function uses the
 *             /AlnBlock_from_protein_genewise_wrap and /AlnBlock_from_TSM_genewise_wrap
 *             functions in this module. 
 *                 
 *
 *
 * Arg:                     gr [UNKN ] genomic region to make genes on [GenomicRegion *]
 * Arg:                    pgl [UNKN ] potential gene list of homologs and start/end points [PotentialGeneList *]
 * Arg:            alg_protein [UNKN ] algorithm to use with protein comparisons [int]
 * Arg:                alg_hmm [UNKN ] algorithm to use with HMM comparisons [int]
 * Arg:               prot_thr [UNKN ] Threshold under which genes are not predicted [double]
 * Arg:                hmm_thr [UNKN ] Undocumented argument [double]
 * Arg:                   comp [UNKN ] Comparison matrix for protein sequences [CompMat *]
 * Arg:                    gap [UNKN ] Gap cost for protein sequences [int]
 * Arg:                    ext [UNKN ] Extension cost for protein sequences [int]
 * Arg:                  gpara [UNKN ] Gene parameters [GeneParameter21 *]
 * Arg:                    rmd [UNKN ] Random Model to compare against [RandomModelDNA *]
 * Arg:                  inter [UNKN ] Random Model for intergenic regions [RandomModelDNA *]
 * Arg:        fetch_from_pipe [UNKN ] For potential genes with just a name, a pipe to get sequences from [char *]
 * Arg:            should_free [UNKN ] Free the protein sequences after resolving it [boolean]
 * Arg:              make_name [FUNCP] a pointer to a function which makes the name of the gene [char *(*make_name]
 * Arg:            bit_cut_off [UNKN ] Undocumented argument [double]
 * Arg:                   dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
# line 668 "gwrap.dy"
int resolve_PotentialGenes_on_GenomicRegion(GenomicRegion * gr,PotentialGeneList * pgl,int alg_protein,int alg_hmm,double prot_thr,double hmm_thr,CompMat * comp,int gap,int ext,GeneParameter21 * gpara,RandomModelDNA * rmd,RandomModelDNA * inter,char * fetch_from_pipe,boolean should_free,char *(*make_name)(Wise2_Genomic * gen,char *,int,Wise2_Gene *),double bit_cut_off,DPRunImpl * dpri)
{
  int i;
  Genomic * gentemp;
  RandomModel * rm;
  AlnBlock * alb;
  int count = 0;
  char buffer[128];
  FILE * pipe;


  if( gr == NULL || pgl == NULL ) {
    warn("Look - you've passed in some NULL genomicregion or pgl's to resolve_PotentialGenes. You can't expect me to do something with them??");
    return -1;
  }

  if( gr->genomic == NULL ) {
    warn("Your genomic region has no DNA, so I can't put any genes onto it!");
    return -1;
  }

  rm = default_RandomModel();

  for(i=0;i<pgl->len;i++) {
    auto PotentialGene * pg;

    pg = pgl->pg[i];

    if( pg->tsm == NULL && pg->homolog == NULL && fetch_from_pipe == NULL) {
      warn("You have neither a HMM nor a homolog for a potential gene, or a fetchable name. Yikes!");
      continue;
    }

    if( pg->tsm == NULL && pg->homolog == NULL ) {
      /* lets fetch from pipe */

      sprintf(buffer,fetch_from_pipe,pg->name);

      pipe = popen(buffer,"r");

      pg->homolog = read_fasta_Protein(pipe);

      if( pg->homolog == NULL ) {
	warn("Could not read protein [%s] with pipe [%s]",pg->name,fetch_from_pipe);
	pclose(pipe);
	continue;
      }
      pclose(pipe);
    }
      
	


    /*
     * ok - use guess_start and end to pick out region to feed genewise.
     * if end < start then it truncates (clever eh!)
     */

    /*    fprintf(stderr,"Here we go %d and %d\n",pg->guess_start,pg->guess_end); */
    if( pg->guess_start == -1 || pg->guess_start < 0) 
      pg->guess_start = 0;
    if( pg->guess_end < 0) 
      pg->guess_end = gr->genomic->baseseq->len-1;


    if( pg->guess_end == -1 || pg->guess_end > gr->genomic->baseseq->len-1) 
      pg->guess_end = gr->genomic->baseseq->len-1;
    if( pg->guess_start > gr->genomic->baseseq->len-1) 
      pg->guess_start = gr->genomic->baseseq->len-1;


    /*fprintf(stderr,"And now we go %d and %d\n",pg->guess_start,pg->guess_end); */
    gentemp = truncate_Genomic(gr->genomic,pg->guess_start,pg->guess_end);

    if( gentemp == NULL ) {
      warn("Cannot make genomic truncation!");
      continue;
    }

    /*   gentemp->baseseq->name = stringalloc("TempGene");
	 write_fasta_Sequence(gentemp->baseseq,stdout);
    */

    /*
     * decide whether this is protein sequence or HMM
     */

    if( pg->tsm != NULL && pg->homolog != NULL ) {
      warn("Nope - you have a potential gene with both a threestatemodel and homolog thing. Taking the HMM!");
    }

    if( pg->tsm != NULL ) {
      
      log_full_error(INFO,0,"Using HMM  [%s] [%d/%d] in region %d %d",pg->tsm->name,i+1,pgl->len,pg->guess_start,pg->guess_end);
      alb = AlnBlock_from_TSM_genewise_wrap(pg->tsm,gentemp,gpara,rmd,inter,TRUE,alg_hmm,1.0,1,pg,dpri,NULL);
      
      if( alb == NULL ) {
	warn("In attempting to map a region of %s to %s, got no alignment. This seems like a bad error!",gr->genomic->baseseq->name,pg->tsm->name);
	continue;
      }

      /*** could be multiple genes ***/

      count += add_Genes_to_GenomicRegion_GeneWise(gr,pg->guess_start,pg->guess_end,alb,pg->tsm->name,FALSE,make_name);

    } else {
      
      log_full_error(INFO,0,"Using protein [%s] [%d/%d] in region %d %d",pg->homolog->baseseq->name,i+1,pgl->len,pg->guess_start,pg->guess_end);


      alb = AlnBlock_from_protein_genewise_wrap(pg->homolog,gentemp,comp,gap,ext,gpara,rmd,inter,alg_protein,pg->is_global,TRUE,rm,1.0,pg,dpri,NULL);
      
      
      if( alb == NULL ) {
	warn("In attempting to map a region of %s to %s, got no alignment. This seems like a bad error!",gr->genomic->baseseq->name,pg->homolog->baseseq->name);
	continue;
      }

      /*** check the score against the bit_cut_off ***/

      pg->bitscore = Score2Bits(alb->score); /* wont work well with 21:93 */

      if( pg->bitscore >= bit_cut_off) {
	/*** could be multiple genes ***/

	count += add_Genes_to_GenomicRegion_GeneWise(gr,gentemp->baseseq->offset,gentemp->baseseq->end,alb,pg->homolog->baseseq->name,FALSE,make_name);
      }

    }

    free_Genomic(gentemp);
    if( should_free == TRUE ) {
      pg->homolog = free_Protein(pg->homolog);
    }

    
    if( should_free == TRUE)
      alb = free_AlnBlock(alb);
    else pg->alb = alb;
  }
    
  free_RandomModel(rm);
  return count;

}
      
    
  
/* Function:  GeneParameter21_wrap(gf,subs_error,indel_error,rmd,use_modelled_codon,use_modelled_splice,tie_intron_prob,ct,rnd_loop,cds_loop,rnd_to_model,link_loop,link_to_model)
 *
 * Descrip:    A general wrap over the production of parameters for the
 *             GeneWise programs. The geneparameter21 holds all the parameters,
 *             and can be approximated for the 6:23 and 4:21 algorithms
 *
 *             This function is the best way to make a GeneParameter21 object
 *             as all the different options for how to make it or modify its
 *             contents are laid out as arguments to this function
 *
 *
 * Arg:                         gf [READ ] Gene Frequency data structure, holding counts for splice sites etc [GeneFrequency21 *]
 * Arg:                 subs_error [UNKN ] substitution error on the dna sequence [double]
 * Arg:                indel_error [UNKN ] rough estimate of the insertion/deletion per base error rate [double]
 * Arg:                        rmd [UNKN ] the random model of the DNA that is used [RandomModelDNA *]
 * Arg:         use_modelled_codon [UNKN ] if TRUE, model codon frequency [boolean]
 * Arg:        use_modelled_splice [UNKN ] if TRUE, make splice models from gf parameters [boolean]
 * Arg:            tie_intron_prob [UNKN ] Undocumented argument [boolean]
 * Arg:                         ct [UNKN ] codon table which is used for codon->aa mapping [CodonTable *]
 * Arg:                   rnd_loop [UNKN ] Undocumented argument [Probability]
 * Arg:                   cds_loop [UNKN ] Undocumented argument [Probability]
 * Arg:               rnd_to_model [UNKN ] Undocumented argument [Probability]
 * Arg:                  link_loop [UNKN ] Undocumented argument [Probability]
 * Arg:              link_to_model [UNKN ] Undocumented argument [Probability]
 *
 * Return [UNKN ]  A newly allocated structure [GeneParameter21 *]
 *
 */
# line 834 "gwrap.dy"
GeneParameter21 * GeneParameter21_wrap(GeneFrequency21 * gf,double subs_error,double indel_error,RandomModelDNA * rmd,boolean use_modelled_codon,boolean use_modelled_splice,boolean tie_intron_prob,CodonTable * ct,Probability rnd_loop,Probability cds_loop,Probability rnd_to_model,Probability link_loop,Probability link_to_model)
{
  GeneParameter21 * out;
  int i;

  out = GeneParameter21_from_GeneFrequency21(gf,ct,rmd,rnd_loop,cds_loop,rnd_to_model,link_loop,link_to_model);


  if( use_modelled_codon == FALSE ) {
    out->cm = free_CodonMapper(out->cm);
    out->cm = flat_CodonMapper(ct);
  }

  if( use_modelled_splice == FALSE )  {
    out->cses = free_ComplexSequenceEvalSet(out->cses);
    out->cses = default_genomic_ComplexSequenceEvalSet();
    out->modelled_splice = FALSE;
  }

  if( tie_intron_prob == TRUE ) {
    for(i=0;i<5;i++) 
      out->gp->central[i] = rmd->base[i];
  }

  /*** errors ***/

  sprinkle_errors_over_CodonMapper(out->cm,subs_error);

  add_flat_error_probabilities_GeneParser21(out->gp,indel_error);

  GeneParser21_fold_in_RandomModelDNA(out->gp,rmd);

  fold_in_RandomModelDNA_into_RandomCodon(out->rc,rmd);

  return out;
}

/* Function:  AlnBlock_from_protein_genewise_wrap(protein,dna,comp,gap,ext,gpara,rmd,intergenic,alg,is_global,use_syn,rm,allN,pg,dpri,pal)
 *
 * Descrip:    A function which aligns a Protein sequecne to a Genomic sequence
 *             under the Comparison matrix comp and the gene paras in gpara.
 *
 *             This is the best function for accessing GeneWise functionality
 *             for a protein to dna comparison, allowing for introns.
 *
 *             To make the protein object, you will first read in a generic
 *             sequence object using something like read_fasta_Sequence and
 *             then convert it to a protein object using new_Protein_from_Sequence
 *
 *             To make the genomic object, you will first read in a generic
 *             sequence object using something like read_fasta_Sequence and
 *             then convert it to a genomic object using new_Genomic_from_Sequence
 *
 *             To make a CompMat object you will use read_Blast_file_CompMat
 *             from the compmat module. It is likely, if the Wise2 enviroment
 *             has been set up correctly that read_Blast_file_CompMat("blosum62.bla")
 *             will be fine. You should at the moment only use halfbit matrices
 *             (blosum62 is one such matrix)
 *
 *             To make the necessary random modules use the default construtors
 *             in the randommodel module
 *
 *             To make the gene parameter object use the GeneParameter21_wrap
 *             function found in this module. It will need GeneFrequencies
 *             read in using the read_GeneFrequency21_file function in
 *             the genefrequency module.  Again if Wise2 has been set up
 *             correctly, read_GeneFrequency21_file("human.gf") should work
 *
 *             To again a valid algorithm type use gwrap_alg_type_from_string
 *             found in this module. gwrap_alg_type_from_string("623") would
 *             be a good choice
 *
 *
 *             This function basically makes a threestatemodel (standard HMM) from
 *             the protein and the comparison matrix with the *scary* assumption that
 *             the comparison matrix is in half bit form. It then calls 
 *              /AlnBlock_from_TSM_genewise_wrap to do the nasty stuff. 
 *
 *
 * Arg:           protein [UNKN ] protein sequence used in the comparison [Protein *]
 * Arg:               dna [UNKN ] genomic DNA sequence used  [Genomic *]
 * Arg:              comp [UNKN ] protein comparison matrix *in half bits* [CompMat *]
 * Arg:               gap [UNKN ] gap penalty (negative) [int]
 * Arg:               ext [UNKN ] extension penalty (negative) [int]
 * Arg:             gpara [UNKN ] Gene parameters. [GeneParameter21 *]
 * Arg:               rmd [UNKN ] models to be compared to [RandomModelDNA *]
 * Arg:        intergenic [UNKN ] model of random dna between genes [RandomModelDNA *]
 * Arg:               alg [UNKN ] algorithm type [int]
 * Arg:         is_global [UNKN ] has now become flag for local/global/end-biased switch [int]
 * Arg:           use_syn [UNKN ] Undocumented argument [boolean]
 * Arg:                rm [UNKN ] Undocumented argument [RandomModel *]
 * Arg:              allN [UNKN ] Undocumented argument [Probability]
 * Arg:                pg [READ ] Potential gene - could be NULL - if rough exon positions are known  [PotentialGene *]
 * Arg:              dpri [UNKN ] Undocumented argument [DPRunImpl *]
 * Arg:               pal [WRITE] Raw alginment to be saved if non-NULL [PackAln **]
 *
 * Return [UNKN ]  Undocumented return value [AlnBlock *]
 *
 */
# line 924 "gwrap.dy"
AlnBlock * AlnBlock_from_protein_genewise_wrap(Protein * protein,Genomic * dna,CompMat * comp,int gap,int ext,GeneParameter21 * gpara,RandomModelDNA * rmd,RandomModelDNA * intergenic,int alg,int is_global,boolean use_syn,RandomModel * rm,Probability allN,PotentialGene * pg,DPRunImpl * dpri,PackAln ** pal)
{
  ThreeStateModel * tsm;
  RandomModel * rm2;
  AlnBlock * out;

  if( protein == NULL || dna == NULL || comp == NULL || gpara == NULL || rmd == NULL ){
    warn("trappable error in PackAln from protein sequence, passed some NULL objects, Complain!");
    return NULL;
  }

  rm2 = default_RandomModel();
  
  if( is_global == 1) 
    tsm = global_ThreeStateModel_from_half_bit_Sequence(protein,comp,rm2,gap,ext);
  else {
    tsm = ThreeStateModel_from_half_bit_Sequence(protein,comp,rm2,gap,ext);
    if( is_global > 2 ) {
      set_startend_policy_ThreeStateModel(tsm,is_global,30,halfbit2Probability(-15));
      /* set end bias */
    }
	 
  }
  
  out = AlnBlock_from_TSM_genewise_wrap(tsm,dna,gpara,rmd,intergenic,use_syn,alg,allN,1,pg,dpri,pal);

  free_ThreeStateModel(tsm);
  free_RandomModel(rm2);

  return out;

}


/* Function:  AlnBlock_from_TSM_genewise_wrap(tsm,gen,gpara,rmd,intergenic,use_syn,alg,allN,flat_insert,pg,dpri,palpoi)
 *
 * Descrip:    A function which aligns a protein HMM (as found
 *             in my threestatemodel structure) to a genomic DNA 
 *             sequence. 
 *
 *             At the moment you are unlikely to be reading in the
 *             HMM structure yourself, so this is not something
 *             you will be doing.
 *
 *             The core algorithms for each method are found in
 *             genewise21/geneloop21 etc files. 
 *
 *
 *
 * Arg:                tsm [UNKN ] protein TSM to be used in the comparison [ThreeStateModel *]
 * Arg:                gen [UNKN ] genomic DNA sequence used  [Genomic *]
 * Arg:              gpara [UNKN ] Gene parameters. [GeneParameter21 *]
 * Arg:                rmd [UNKN ] models to be compared to [RandomModelDNA *]
 * Arg:         intergenic [UNKN ] model of random dna between genes [RandomModelDNA *]
 * Arg:            use_syn [UNKN ] use a synchronous null model [boolean]
 * Arg:                alg [UNKN ] algorithm type [int]
 * Arg:               allN [UNKN ] Undocumented argument [Probability]
 * Arg:        flat_insert [UNKN ] Undocumented argument [boolean]
 * Arg:                 pg [READ ] Potential gene - could be NULL - if rough exon positions are known [PotentialGene *]
 * Arg:               dpri [UNKN ] Undocumented argument [DPRunImpl *]
 * Arg:             palpoi [WRITE] Raw alginment to be saved if non-NULL [PackAln **]
 *
 * Return [UNKN ]  Undocumented return value [AlnBlock *]
 *
 */
# line 981 "gwrap.dy"
AlnBlock * AlnBlock_from_TSM_genewise_wrap(ThreeStateModel * tsm,Genomic * gen,GeneParameter21 * gpara,RandomModelDNA * rmd,RandomModelDNA * intergenic,boolean use_syn,int alg,Probability allN,boolean flat_insert,PotentialGene * pg,DPRunImpl * dpri,PackAln ** palpoi)
{
  AlnBlock * out = NULL;
  PackAln * pal = NULL;
  ComplexSequence * cs = NULL;
  GeneWise * gw = NULL;
  GeneWiseScore * gws = NULL;
  RandomCodonScore * rcs = NULL ;
  GeneParser21Score  * gps = NULL;
  GeneParser4Score * gp4s = NULL;
  RandomModelDNAScore * ids = NULL;
  DPEnvelope * dpenv;
  Sequence * dna;
  cDNAParserScore * cps = NULL; /* for estwise type algorithms */
  GwLite * gwl = NULL;
  GwLiteScore * gwls = NULL;
  ComplexSequenceEval * tempcse;
  ComplexSequenceEvalSet * cses;
  dna = gen->baseseq;
  

  assert(tsm);
  assert(gen);
  assert(gpara);
  assert(rmd);
  assert(gpara->rc);


  /*show_Genomic(gen,stderr);*/

  /*show_GeneParser21(gpara->gp,stderr); */



  if( tsm == NULL || dna == NULL || gpara == NULL || rmd == NULL){
    warn("trappable error in PackAln from TSM  sequence, passed some NULL objects, Complain!");
    return NULL;
  }

  /*** prepare cses ***/

  if( prepare_ComplexSequenceEvalSet(gpara->cses) == FALSE ) {
    warn("Unable to prepare complexsequenceevalset in TMS2DNA wrap");
    goto exit;
  }

  if( alg == GWWRAP_333 || alg == GWWRAP_333L ) {
    cses = default_cDNA_ComplexSequenceEvalSet();
    cs = new_ComplexSequence(gen->baseseq,cses);
    free_ComplexSequenceEvalSet(cses);
  } else if ( alg == GWWRAP_6LITE ) {
    /* yup. This is scary. */
    tempcse = gpara->cses->cse[1];
    gpara->cses->cse[1] = codon64_number_ComplexSequenceEval();
    cs = new_ComplexSequence(gen->baseseq,gpara->cses);
    free_ComplexSequenceEval(gpara->cses->cse[1]);
    gpara->cses->cse[1] = tempcse;
  } else {
    if( (cs=evaluate_ComplexSequence_Genomic(gen,gpara->cses,0,Probability2Score(0.01))) == FALSE ) {
      warn("Unable to make ComplexSequence in TMS2DNA wrap");
      goto exit;
    }
  }

  /*show_ComplexSequence(cs,stderr);*/


  if( (gw=GeneWise_from_ThreeStateModel(tsm,gpara->gp,gpara->cm,allN,gpara->gwcm)) == NULL) {
    warn("Unable to make GeneWise model");
    goto exit;
  }

  if( gpara->modelled_splice == FALSE) {
    flatten_balance_scores_GeneWise(gw);
  }
	

  if( pg == NULL || (pg->guess_start == -1) ) {
    dpenv = NULL ; 
  } else {
    info("Using DPEnvelope over matrix");
    pg->guess_start = dna->offset;
    pg->guess_end = dna->end;
    pg->query_length = tsm->len;
    dpenv = DPEnvelope_from_PotentialGene(pg);
    /* show_DPEnvelope(dpenv,stderr); */
  }

  /*  show_GeneWiseSegment(gw->seg[0],stderr); */
  
  if( use_syn == TRUE ) {
    if( tsm->rm == NULL ) {
      warn("Ugh - a threestatemodel without a random model. Not in this code matey");
      goto exit;
    }


    GeneWise_fold_in_synchronised_RandomModel(gw,tsm->rm,gpara->cm,gpara->ct,0.5);
    flatten_RandomCodon(gpara->rc);
  } else {
    GeneWise_fold_in_RandomModelDNA(gw,rmd);
  }

  if( alg == GWWRAP_6LITE ) {
    gwl = GwLite_from_GeneWise(gw);
    gwls = GwLiteScore_from_GwLite(gwl);
  }


  if( flat_insert == TRUE ) {
    check_flat_insert(gw,1,0,gpara->cm->ct);
  }

  if( (gws = GeneWiseScore_from_GeneWise(gw)) == NULL) {
    warn("Unable to make GeneWiseScore model");
    goto exit;
  }


  if( (gps = GeneParser21Score_from_GeneParser21(gpara->gp)) == NULL) {
    warn("Unable to make GeneParserScore model");
    goto exit;
  }


  if( (rcs = RandomCodonScore_from_RandomCodon(gpara->rc)) == NULL) {
    warn("Unable to make RandomCodonScore model");
    goto exit;
  }

  ids = folded_RandomModelDNAScore_from_2RMD(intergenic,rmd);

  gp4s = GeneParser4Score_from_GeneParser21Score(gps);


  if( alg == GWWRAP_333 || alg == GWWRAP_333L ) {
    cps = cDNAParserScore_from_GeneParser21Score(gps);
  }

  switch(alg) {
  case GWWRAP_2193 :

    pal = PackAln_bestmemory_GeneWise21(gws,cs,gps,rcs,ids,dpenv,dpri);

    out = convert_PackAln_to_AlnBlock_GeneWise21(pal);
    break;

  case GWWRAP_2193I :

    warn("Algorithm currently disabled! Sorry!");
    break;

    /**
    pal = PackAln_dc_build_GeneLinker21(gws,cs,gps,rcs,ids);
    if(pal == NULL)
      goto exit;
    out = convert_PackAln_to_AlnBlock_GeneLinker21(pal,NULL);
    break;
    **/

  case GWWRAP_2193L :

    pal = PackAln_bestmemory_GeneLoop21(gws,cs,gps,rcs,ids,dpenv,dpri);
    if(pal == NULL)
      goto exit;
    out = convert_PackAln_to_AlnBlock_GeneLoop21(pal);
    break;
  case GWWRAP_623L :

    pal = PackAln_bestmemory_GeneLoop6(gws,cs,gp4s,dpenv,dpri);
    if(pal == NULL)
      goto exit;
    out = convert_PackAln_to_AlnBlock_GeneLoop6(pal);
    break;
  case GWWRAP_623 :

    pal = PackAln_bestmemory_GeneWise6(gws,cs,gp4s,dpenv,dpri);
    if(pal == NULL)
      goto exit;
    out = convert_PackAln_to_AlnBlock_GeneWise6(pal);
    break;

  case GWWRAP_333 :
    cps = cDNAParserScore_from_GeneParser21Score(gps);
    pal = PackAln_bestmemory_EstWise3(gws,cs,cps,dpenv,dpri);
    if(pal == NULL)
      goto exit;
    out = convert_PackAln_to_AlnBlock_EstWise3(pal);
    break;

  case GWWRAP_333L :
    cps = cDNAParserScore_from_GeneParser21Score(gps);
    pal = PackAln_bestmemory_EstLoop3(gws,cs,cps,dpenv,dpri);
    if(pal == NULL)
      goto exit;
    out = convert_PackAln_to_AlnBlock_EstLoop3(pal);
    break;
    
  case GWWRAP_421 :
    
    pal = PackAln_bestmemory_GeneWise4(gws,cs,gp4s,dpenv,dpri);
    if(pal == NULL)
      goto exit;
    out = convert_PackAln_to_AlnBlock_GeneWise4(pal);
    break;
    
  case GWWRAP_6LITE :

    pal = PackAln_bestmemory_GeneLiteModel(gwls,cs,gp4s,dpenv,dpri);
    if(pal == NULL)
      goto exit;
    out = convert_PackAln_to_AlnBlock_GeneLiteModel(pal);
    GwLite_AlnBlock_surgery(out);
    break;

  default :
    warn("A major problem. No valid algorithm type passed in");
    goto exit;
  }

  map_phase0_codons_AlnBlock_GeneWise(out,gws,cs);

  if( palpoi != NULL ) {
    *palpoi = pal;
    pal = NULL;
  }

  goto exit;




  exit :

  
  if(pal != NULL)
    pal = free_PackAln(pal);
  if( cps != NULL )
    cps = free_cDNAParserScore(cps);
  if( ids != NULL )
    ids = free_RandomModelDNAScore(ids);
  if(cs != NULL )
    cs = free_ComplexSequence(cs);
  if(gw != NULL )
    gw = free_GeneWise(gw);
  if(gws != NULL )
    gws = free_GeneWiseScore(gws);
  if(gps != NULL )
    free_GeneParser21Score(gps);
  if(rcs != NULL )
    rcs = free_RandomCodonScore(rcs);
  if(gp4s != NULL)
    gp4s = free_GeneParser4Score(gp4s);

  return out;
}

/* Function:  Hscore_from_TSM_genewise(tdb,gdb,gpara,rmd,intergenic,use_syn,alg,bits_cutoff,allN,report_level,die_on_error,flat_insert,dbsi)
 *
 * Descrip:    Runs a database search of the genewise algorithm. 
 *
 *             This makes a high score object which you can then use 
 *             to retrieve enteries as well as print out the top score (!)
 *
 *
 * Arg:                 tdb [READ ] a database of profileHMMs  [ThreeStateDB *]
 * Arg:                 gdb [READ ] a database of genomic sequence [GenomicDB *]
 * Arg:               gpara [READ ] geneparameters [GeneParameter21 *]
 * Arg:                 rmd [READ ] random model to be compared with in non syn mode [RandomModelDNA *]
 * Arg:          intergenic [READ ] random model of intergenic DNA (usually the same as rmd) [RandomModelDNA *]
 * Arg:             use_syn [UNKN ] use synchronous random model [boolean]
 * Arg:                 alg [UNKN ] algorithm type [int]
 * Arg:         bits_cutoff [UNKN ] cutoff in bits of the scores to store [double]
 * Arg:                allN [UNKN ] Undocumented argument [Probability]
 * Arg:        report_level [UNKN ] stagger rate of reporting progress on stderr  -1 means never [int]
 * Arg:        die_on_error [UNKN ] if true, exits on error (not used at the moment) [boolean]
 * Arg:         flat_insert [UNKN ] Undocumented argument [boolean]
 * Arg:                dbsi [UNKN ] Undocumented argument [DBSearchImpl *]
 *
 * Return [UNKN ]  a new Hscore object of the entire db search [Hscore *]
 *
 */
# line 1256 "gwrap.dy"
Hscore * Hscore_from_TSM_genewise(ThreeStateDB * tdb,GenomicDB * gdb,GeneParameter21 * gpara,RandomModelDNA * rmd,RandomModelDNA * intergenic,boolean use_syn,int alg,double bits_cutoff,Probability allN,int report_level,boolean die_on_error,boolean flat_insert,DBSearchImpl * dbsi)
{
  Hscore * out = NULL;
  GeneWiseDB * gwdb = NULL;
  cDNADB * cdb = NULL;
  cDNAParserScore * cps= NULL;
  GeneParser21Score * gps = NULL;
  GeneParser4Score * gp4s = NULL;
  RandomCodonScore * rcs = NULL;
  RandomModelDNAScore * ids = NULL;
  cDNA * temp;
  Search_Return_Type ret;
  ComplexSequenceEval * tempcse;
  
  ret = SEARCH_ERROR;

  gwdb = new_GeneWiseDB(tdb,gpara,rmd,use_syn,allN);
  gwdb->flat_insert = flat_insert;
  if( gwdb == NULL ) {
    warn("Could not build a new GeneWiseDB from the objects provided. Exiting without completing the search");
    goto exit;
  }


  if( (gps = GeneParser21Score_from_GeneParser21(gpara->gp)) == NULL) {
    warn("Unable to make GeneParserScore model");
    goto exit;
  }
  

  if( (rcs = RandomCodonScore_from_RandomCodon(gpara->rc)) == NULL) {
    warn("Unable to make RandomCodonScore model");
    goto exit;
  }

  ids = folded_RandomModelDNAScore_from_2RMD(intergenic,rmd);

  gp4s = GeneParser4Score_from_GeneParser21Score(gps);

  if( alg == GWWRAP_333 || alg == GWWRAP_333L ) {
    /* could be a single dna sequence */
    if( gdb->is_single_seq == TRUE ) {
      temp = cDNA_from_Sequence(hard_link_Sequence(gdb->forw->seq));
      cdb = new_cDNADB_from_single_seq(temp);
      free_cDNA(temp); /* hard linked by database */
    } else {
      cdb = new_cDNADB(gdb->sdb);
    }
  }  

  /*** allocate Hscore structure ***/

  out = std_bits_Hscore(bits_cutoff,report_level);

  switch(alg) {
  case GWWRAP_2193 :

    ret = Wise2_search_GeneWise21(dbsi,out,gwdb,gdb,gps,rcs,ids);
    break;

  case GWWRAP_2193I :

    warn("Algorithm currently disabled! Sorry!");
    break;


  case GWWRAP_2193L :

    warn("Unable to do db looping searches (sort of - pointless - )!");
    break;
  case GWWRAP_623L :
    warn("Unable to do db looping searches (sort of - pointless - )!");
    break;

  case GWWRAP_623 :

    ret = Wise2_search_GeneWise6(dbsi,out,gwdb,gdb,gp4s);
    break;

  case GWWRAP_6LITE :

    ret = Wise2_search_GeneLiteModel(dbsi,out,gwdb,gdb,gp4s);
    break;

  case GWWRAP_333 :
    cps = cDNAParserScore_from_GeneParser21Score(gps);
    ret = search_EstWise3(dbsi,out,gwdb,cdb,cps);

    break;

  case GWWRAP_333L :
    warn("Unable to do db looping searches (sort of - pointless - )!");
    break;
    

  case GWWRAP_421 :

    ret = Wise2_search_GeneWise4(dbsi,out,gwdb,gdb,gp4s);
    break;


  default :
    warn("A major problem. No valid algorithm type passed in");
    goto exit;
  }

  goto exit;




  exit :

    /* for 6LITE leaking a tiny amount of memory. Oh well... */
  if( ids != NULL )
    ids = free_RandomModelDNAScore(ids);
  if( cps != NULL ) 
    free_cDNAParserScore(cps);
  if( cdb != NULL ) 
    free_cDNADB(cdb);
  if(gps != NULL )
    free_GeneParser21Score(gps);
  if(rcs != NULL )
    rcs = free_RandomCodonScore(rcs);
  if(gp4s != NULL)
    gp4s = free_GeneParser4Score(gp4s);
  if( gwdb != NULL ) {
    free_GeneWiseDB(gwdb);
  }

  if( die_on_error == TRUE  && ret == SEARCH_ERROR) {
    if( out != NULL ) {
      free_Hscore(out);
    } 
    return NULL;
  }


  return out;
}

/* Function:  cDNAParserScore_from_GeneParser21Score(gps)
 *
 * Descrip:    Makes a cdna parser from a genewise parser. Basically
 *             copies the indel penalties across.
 *
 *
 * Arg:        gps [UNKN ] Undocumented argument [GeneParser21Score *]
 *
 * Return [UNKN ]  Undocumented return value [cDNAParserScore *]
 *
 */
# line 1402 "gwrap.dy"
cDNAParserScore * cDNAParserScore_from_GeneParser21Score(GeneParser21Score * gps)
{
  cDNAParserScore * out;

  out = cDNAParserScore_alloc();

  out->trans[PCD_INSERT_2_BASE] = gps->transition[GP21_INSERT_2_BASE];
  out->trans[PCD_INSERT_1_BASE] = gps->transition[GP21_INSERT_1_BASE];
  out->trans[PCD_DELETE_2_BASE] = gps->transition[GP21_DELETE_2_BASE];
  out->trans[PCD_DELETE_1_BASE] = gps->transition[GP21_DELETE_1_BASE];


  return out;
}

/* Function:  gwrap_alg_type_from_string(str)
 *
 * Descrip:    Gives you the integer interpretation from
 *             the string, which is one of
 *             2193 2193L, 623, 623L, 421, 2193LINK
 *
 *             This integer can then be passed into routines
 *             like AlnBlock_from_protein_genewise_wrap
 *
 *
 *
 * Arg:        str [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
# line 1426 "gwrap.dy"
int gwrap_alg_type_from_string(char * str)
{
  int t;

  t = get_number_from_slashed_string(str,"2193/2193L/623/623L/421/2193LINK/333/333L/6LITE");

  switch (t) {
  case 0 : return GWWRAP_2193;
  case 1 : return GWWRAP_2193L;
  case 2 : return GWWRAP_623;
  case 3 : return GWWRAP_623L;
  case 4 : return GWWRAP_421;
  case 5 : return GWWRAP_2193I;
  case 6 : return GWWRAP_333;
  case 7 : return GWWRAP_333L;
  case 8 : return GWWRAP_6LITE;
  default : warn("Cannot convert string %s into a valid genewise algorithm type\n",str);
    return -1;
  }
}


# line 1457 "gwrap.c"
/* Function:  hard_link_PotentialExon(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [PotentialExon *]
 *
 * Return [UNKN ]  Undocumented return value [PotentialExon *]
 *
 */
PotentialExon * hard_link_PotentialExon(PotentialExon * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a PotentialExon object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  PotentialExon_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [PotentialExon *]
 *
 */
PotentialExon * PotentialExon_alloc(void) 
{
    PotentialExon * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(PotentialExon *) ckalloc (sizeof(PotentialExon))) == NULL)  {  
      warn("PotentialExon_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
    out->tstart = 0; 
    out->tend = 0;   
    out->qstart = 0; 
    out->qend = 0;   


    return out;  
}    


/* Function:  free_PotentialExon(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [PotentialExon *]
 *
 * Return [UNKN ]  Undocumented return value [PotentialExon *]
 *
 */
PotentialExon * free_PotentialExon(PotentialExon * obj) 
{


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a PotentialExon obj. Should be trappable"); 
      return NULL;   
      }  


    if( obj->dynamite_hard_link > 1)     {  
      obj->dynamite_hard_link--; 
      return NULL;   
      }  


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_PotentialTranscript(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_PotentialTranscript
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [PotentialExon **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_PotentialTranscript(PotentialExon ** list,int i,int j)  
{
    PotentialExon * temp;    
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_PotentialTranscript(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_PotentialTranscript which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [PotentialExon **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_PotentialTranscript(PotentialExon ** list,int left,int right,int (*comp)(PotentialExon * ,PotentialExon * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_PotentialTranscript(list,left,(left+right)/2);  
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_PotentialTranscript (list,++last,i);    
      }  
    swap_PotentialTranscript (list,left,last);   
    qsort_PotentialTranscript(list,left,last-1,comp);    
    qsort_PotentialTranscript(list,last+1,right,comp);   
}    


/* Function:  sort_PotentialTranscript(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_PotentialTranscript
 *
 *
 * Arg:         obj [UNKN ] Object containing list [PotentialTranscript *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_PotentialTranscript(PotentialTranscript * obj,int (*comp)(PotentialExon *, PotentialExon *)) 
{
    qsort_PotentialTranscript(obj->pex,0,obj->len-1,comp);   
    return;  
}    


/* Function:  expand_PotentialTranscript(obj,len)
 *
 * Descrip:    Really an internal function for add_PotentialTranscript
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [PotentialTranscript *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_PotentialTranscript(PotentialTranscript * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_PotentialTranscript called with no need");    
      return TRUE;   
      }  


    if( (obj->pex = (PotentialExon ** ) ckrealloc (obj->pex,sizeof(PotentialExon *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_PotentialTranscript, returning FALSE");  
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_PotentialTranscript(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [PotentialTranscript *]
 * Arg:        add [OWNER] Object to add to the list [PotentialExon *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_PotentialTranscript(PotentialTranscript * obj,PotentialExon * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_PotentialTranscript(obj,obj->len + PotentialTranscriptLISTLENGTH) == FALSE) 
        return FALSE;    
      }  


    obj->pex[obj->len++]=add;    
    return TRUE; 
}    


/* Function:  flush_PotentialTranscript(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [PotentialTranscript *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_PotentialTranscript(PotentialTranscript * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->pex[i] != NULL)   {  
        free_PotentialExon(obj->pex[i]); 
        obj->pex[i] = NULL;  
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  PotentialTranscript_alloc_std(void)
 *
 * Descrip:    Equivalent to PotentialTranscript_alloc_len(PotentialTranscriptLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [PotentialTranscript *]
 *
 */
PotentialTranscript * PotentialTranscript_alloc_std(void) 
{
    return PotentialTranscript_alloc_len(PotentialTranscriptLISTLENGTH); 
}    


/* Function:  PotentialTranscript_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [PotentialTranscript *]
 *
 */
PotentialTranscript * PotentialTranscript_alloc_len(int len) 
{
    PotentialTranscript * out;  /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = PotentialTranscript_alloc()) == NULL)  
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->pex = (PotentialExon ** ) ckcalloc (len,sizeof(PotentialExon *))) == NULL)  {  
      warn("Warning, ckcalloc failed in PotentialTranscript_alloc_len"); 
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_PotentialTranscript(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [PotentialTranscript *]
 *
 * Return [UNKN ]  Undocumented return value [PotentialTranscript *]
 *
 */
PotentialTranscript * hard_link_PotentialTranscript(PotentialTranscript * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a PotentialTranscript object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  PotentialTranscript_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [PotentialTranscript *]
 *
 */
PotentialTranscript * PotentialTranscript_alloc(void) 
{
    PotentialTranscript * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(PotentialTranscript *) ckalloc (sizeof(PotentialTranscript))) == NULL)  {  
      warn("PotentialTranscript_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
    out->pex = NULL; 
    out->len = out->maxlen = 0;  


    return out;  
}    


/* Function:  free_PotentialTranscript(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [PotentialTranscript *]
 *
 * Return [UNKN ]  Undocumented return value [PotentialTranscript *]
 *
 */
PotentialTranscript * free_PotentialTranscript(PotentialTranscript * obj) 
{
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a PotentialTranscript obj. Should be trappable");   
      return NULL;   
      }  


    if( obj->dynamite_hard_link > 1)     {  
      obj->dynamite_hard_link--; 
      return NULL;   
      }  
    if( obj->pex != NULL)    {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->pex[i] != NULL) 
          free_PotentialExon(obj->pex[i]);   
        }  
      ckfree(obj->pex);  
      }  


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_PotentialGene(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_PotentialGene
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [PotentialTranscript **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_PotentialGene(PotentialTranscript ** list,int i,int j)  
{
    PotentialTranscript * temp;  
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_PotentialGene(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_PotentialGene which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [PotentialTranscript **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_PotentialGene(PotentialTranscript ** list,int left,int right,int (*comp)(PotentialTranscript * ,PotentialTranscript * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_PotentialGene(list,left,(left+right)/2);    
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_PotentialGene (list,++last,i);  
      }  
    swap_PotentialGene (list,left,last); 
    qsort_PotentialGene(list,left,last-1,comp);  
    qsort_PotentialGene(list,last+1,right,comp); 
}    


/* Function:  sort_PotentialGene(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_PotentialGene
 *
 *
 * Arg:         obj [UNKN ] Object containing list [PotentialGene *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_PotentialGene(PotentialGene * obj,int (*comp)(PotentialTranscript *, PotentialTranscript *)) 
{
    qsort_PotentialGene(obj->pet,0,obj->len-1,comp); 
    return;  
}    


/* Function:  expand_PotentialGene(obj,len)
 *
 * Descrip:    Really an internal function for add_PotentialGene
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [PotentialGene *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_PotentialGene(PotentialGene * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_PotentialGene called with no need");  
      return TRUE;   
      }  


    if( (obj->pet = (PotentialTranscript ** ) ckrealloc (obj->pet,sizeof(PotentialTranscript *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_PotentialGene, returning FALSE");    
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_PotentialGene(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [PotentialGene *]
 * Arg:        add [OWNER] Object to add to the list [PotentialTranscript *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_PotentialGene(PotentialGene * obj,PotentialTranscript * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_PotentialGene(obj,obj->len + PotentialGeneLISTLENGTH) == FALSE) 
        return FALSE;    
      }  


    obj->pet[obj->len++]=add;    
    return TRUE; 
}    


/* Function:  flush_PotentialGene(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [PotentialGene *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_PotentialGene(PotentialGene * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->pet[i] != NULL)   {  
        free_PotentialTranscript(obj->pet[i]);   
        obj->pet[i] = NULL;  
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  PotentialGene_alloc_std(void)
 *
 * Descrip:    Equivalent to PotentialGene_alloc_len(PotentialGeneLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [PotentialGene *]
 *
 */
PotentialGene * PotentialGene_alloc_std(void) 
{
    return PotentialGene_alloc_len(PotentialGeneLISTLENGTH); 
}    


/* Function:  PotentialGene_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [PotentialGene *]
 *
 */
PotentialGene * PotentialGene_alloc_len(int len) 
{
    PotentialGene * out;/* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = PotentialGene_alloc()) == NULL)    
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->pet = (PotentialTranscript ** ) ckcalloc (len,sizeof(PotentialTranscript *))) == NULL)  {  
      warn("Warning, ckcalloc failed in PotentialGene_alloc_len");   
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_PotentialGene(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [PotentialGene *]
 *
 * Return [UNKN ]  Undocumented return value [PotentialGene *]
 *
 */
PotentialGene * hard_link_PotentialGene(PotentialGene * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a PotentialGene object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  PotentialGene_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [PotentialGene *]
 *
 */
PotentialGene * PotentialGene_alloc(void) 
{
    PotentialGene * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(PotentialGene *) ckalloc (sizeof(PotentialGene))) == NULL)  {  
      warn("PotentialGene_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
    out->guess_start = -1;   
    out->guess_end = -1; 
    out->pet = NULL; 
    out->len = out->maxlen = 0;  
    out->is_global = FALSE;  
    out->name = NULL;    
    out->homolog = NULL; 
    out->tsm = NULL; 
    out->alb = NULL; 
    out->bitscore = 0;   
    out->slop_query = 5; 
    out->slop_target = 25;   
    out->query_length = 0;   


    return out;  
}    


/* Function:  free_PotentialGene(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [PotentialGene *]
 *
 * Return [UNKN ]  Undocumented return value [PotentialGene *]
 *
 */
PotentialGene * free_PotentialGene(PotentialGene * obj) 
{
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a PotentialGene obj. Should be trappable"); 
      return NULL;   
      }  


    if( obj->dynamite_hard_link > 1)     {  
      obj->dynamite_hard_link--; 
      return NULL;   
      }  
    if( obj->pet != NULL)    {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->pet[i] != NULL) 
          free_PotentialTranscript(obj->pet[i]); 
        }  
      ckfree(obj->pet);  
      }  
    if( obj->name != NULL)   
      ckfree(obj->name);     
    if( obj->homolog != NULL)    
      free_Protein(obj->homolog);    
    if( obj->tsm != NULL)    
      free_ThreeStateModel(obj->tsm);    
    if( obj->alb != NULL)    
      free_AlnBlock(obj->alb);   


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_PotentialGeneList(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_PotentialGeneList
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [PotentialGene **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_PotentialGeneList(PotentialGene ** list,int i,int j)  
{
    PotentialGene * temp;    
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_PotentialGeneList(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_PotentialGeneList which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [PotentialGene **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_PotentialGeneList(PotentialGene ** list,int left,int right,int (*comp)(PotentialGene * ,PotentialGene * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_PotentialGeneList(list,left,(left+right)/2);    
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_PotentialGeneList (list,++last,i);  
      }  
    swap_PotentialGeneList (list,left,last); 
    qsort_PotentialGeneList(list,left,last-1,comp);  
    qsort_PotentialGeneList(list,last+1,right,comp); 
}    


/* Function:  sort_PotentialGeneList(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_PotentialGeneList
 *
 *
 * Arg:         obj [UNKN ] Object containing list [PotentialGeneList *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_PotentialGeneList(PotentialGeneList * obj,int (*comp)(PotentialGene *, PotentialGene *)) 
{
    qsort_PotentialGeneList(obj->pg,0,obj->len-1,comp);  
    return;  
}    


/* Function:  expand_PotentialGeneList(obj,len)
 *
 * Descrip:    Really an internal function for add_PotentialGeneList
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [PotentialGeneList *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_PotentialGeneList(PotentialGeneList * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_PotentialGeneList called with no need");  
      return TRUE;   
      }  


    if( (obj->pg = (PotentialGene ** ) ckrealloc (obj->pg,sizeof(PotentialGene *)*len)) == NULL)     {  
      warn("ckrealloc failed for expand_PotentialGeneList, returning FALSE");    
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_PotentialGeneList(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [PotentialGeneList *]
 * Arg:        add [OWNER] Object to add to the list [PotentialGene *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_PotentialGeneList(PotentialGeneList * obj,PotentialGene * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_PotentialGeneList(obj,obj->len + PotentialGeneListLISTLENGTH) == FALSE) 
        return FALSE;    
      }  


    obj->pg[obj->len++]=add; 
    return TRUE; 
}    


/* Function:  flush_PotentialGeneList(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [PotentialGeneList *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_PotentialGeneList(PotentialGeneList * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->pg[i] != NULL)    {  
        free_PotentialGene(obj->pg[i]);  
        obj->pg[i] = NULL;   
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  PotentialGeneList_alloc_std(void)
 *
 * Descrip:    Equivalent to PotentialGeneList_alloc_len(PotentialGeneListLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [PotentialGeneList *]
 *
 */
PotentialGeneList * PotentialGeneList_alloc_std(void) 
{
    return PotentialGeneList_alloc_len(PotentialGeneListLISTLENGTH); 
}    


/* Function:  PotentialGeneList_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [PotentialGeneList *]
 *
 */
PotentialGeneList * PotentialGeneList_alloc_len(int len) 
{
    PotentialGeneList * out;/* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = PotentialGeneList_alloc()) == NULL)    
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->pg = (PotentialGene ** ) ckcalloc (len,sizeof(PotentialGene *))) == NULL)   {  
      warn("Warning, ckcalloc failed in PotentialGeneList_alloc_len");   
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_PotentialGeneList(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [PotentialGeneList *]
 *
 * Return [UNKN ]  Undocumented return value [PotentialGeneList *]
 *
 */
PotentialGeneList * hard_link_PotentialGeneList(PotentialGeneList * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a PotentialGeneList object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  PotentialGeneList_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [PotentialGeneList *]
 *
 */
PotentialGeneList * PotentialGeneList_alloc(void) 
{
    PotentialGeneList * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(PotentialGeneList *) ckalloc (sizeof(PotentialGeneList))) == NULL)  {  
      warn("PotentialGeneList_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
    out->pg = NULL;  
    out->len = out->maxlen = 0;  


    return out;  
}    


/* Function:  free_PotentialGeneList(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [PotentialGeneList *]
 *
 * Return [UNKN ]  Undocumented return value [PotentialGeneList *]
 *
 */
PotentialGeneList * free_PotentialGeneList(PotentialGeneList * obj) 
{
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a PotentialGeneList obj. Should be trappable"); 
      return NULL;   
      }  


    if( obj->dynamite_hard_link > 1)     {  
      obj->dynamite_hard_link--; 
      return NULL;   
      }  
    if( obj->pg != NULL) {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->pg[i] != NULL)  
          free_PotentialGene(obj->pg[i]);    
        }  
      ckfree(obj->pg);   
      }  


    ckfree(obj); 
    return NULL; 
}    


/* Function:  access_pg_PotentialGeneList(obj,i)
 *
 * Descrip:    Access members stored in the pg list
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the list [PotentialGeneList *]
 * Arg:          i [UNKN ] Position in the list [int]
 *
 * Return [SOFT ]  Element of the list [PotentialGene *]
 *
 */
PotentialGene * access_pg_PotentialGeneList(PotentialGeneList * obj,int i) 
{
    if( obj == NULL)     {  
      warn("In accessor function pg for object PotentialGeneList, got a NULL object");   
      return NULL;   
      }  
    if( obj->len <= i )  {  
      warn("In accessor function pg for object PotentialGeneList, index %%d is greater than list length %%d",i,obj->len);    
      return NULL;   
      }  
    return obj->pg[i];   
}    


/* Function:  length_pg_PotentialGeneList(obj)
 *
 * Descrip:    discover the length of the list
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the list [PotentialGeneList *]
 *
 * Return [UNKN ]  length of the list [int]
 *
 */
int length_pg_PotentialGeneList(PotentialGeneList * obj) 
{
    if( obj == NULL)     {  
      warn("In length function pg for object PotentialGeneList, got a NULL object"); 
      return -1;     
      }  
    return obj->len;     
}    


/* Function:  replace_guess_start_PotentialGene(obj,guess_start)
 *
 * Descrip:    Replace member variable guess_start
 *             For use principly by API functions
 *
 *
 * Arg:                obj [UNKN ] Object holding the variable [PotentialGene *]
 * Arg:        guess_start [OWNER] New value of the variable [int]
 *
 * Return [SOFT ]  member variable guess_start [boolean]
 *
 */
boolean replace_guess_start_PotentialGene(PotentialGene * obj,int guess_start) 
{
    if( obj == NULL)     {  
      warn("In replacement function guess_start for object PotentialGene, got a NULL object");   
      return FALSE;  
      }  
    obj->guess_start = guess_start;  
    return TRUE; 
}    


/* Function:  access_guess_start_PotentialGene(obj)
 *
 * Descrip:    Access member variable guess_start
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [PotentialGene *]
 *
 * Return [SOFT ]  member variable guess_start [int]
 *
 */
int access_guess_start_PotentialGene(PotentialGene * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function guess_start for object PotentialGene, got a NULL object");  
      return 0;  
      }  
    return obj->guess_start;     
}    


/* Function:  replace_guess_end_PotentialGene(obj,guess_end)
 *
 * Descrip:    Replace member variable guess_end
 *             For use principly by API functions
 *
 *
 * Arg:              obj [UNKN ] Object holding the variable [PotentialGene *]
 * Arg:        guess_end [OWNER] New value of the variable [int]
 *
 * Return [SOFT ]  member variable guess_end [boolean]
 *
 */
boolean replace_guess_end_PotentialGene(PotentialGene * obj,int guess_end) 
{
    if( obj == NULL)     {  
      warn("In replacement function guess_end for object PotentialGene, got a NULL object"); 
      return FALSE;  
      }  
    obj->guess_end = guess_end;  
    return TRUE; 
}    


/* Function:  access_guess_end_PotentialGene(obj)
 *
 * Descrip:    Access member variable guess_end
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [PotentialGene *]
 *
 * Return [SOFT ]  member variable guess_end [int]
 *
 */
int access_guess_end_PotentialGene(PotentialGene * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function guess_end for object PotentialGene, got a NULL object");    
      return 0;  
      }  
    return obj->guess_end;   
}    


/* Function:  access_pet_PotentialGene(obj,i)
 *
 * Descrip:    Access members stored in the pet list
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the list [PotentialGene *]
 * Arg:          i [UNKN ] Position in the list [int]
 *
 * Return [SOFT ]  Element of the list [PotentialTranscript *]
 *
 */
PotentialTranscript * access_pet_PotentialGene(PotentialGene * obj,int i) 
{
    if( obj == NULL)     {  
      warn("In accessor function pet for object PotentialGene, got a NULL object");  
      return NULL;   
      }  
    if( obj->len <= i )  {  
      warn("In accessor function pet for object PotentialGene, index %%d is greater than list length %%d",i,obj->len);   
      return NULL;   
      }  
    return obj->pet[i];  
}    


/* Function:  length_pet_PotentialGene(obj)
 *
 * Descrip:    discover the length of the list
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the list [PotentialGene *]
 *
 * Return [UNKN ]  length of the list [int]
 *
 */
int length_pet_PotentialGene(PotentialGene * obj) 
{
    if( obj == NULL)     {  
      warn("In length function pet for object PotentialGene, got a NULL object");    
      return -1;     
      }  
    return obj->len;     
}    


/* Function:  replace_is_global_PotentialGene(obj,is_global)
 *
 * Descrip:    Replace member variable is_global
 *             For use principly by API functions
 *
 *
 * Arg:              obj [UNKN ] Object holding the variable [PotentialGene *]
 * Arg:        is_global [OWNER] New value of the variable [boolean]
 *
 * Return [SOFT ]  member variable is_global [boolean]
 *
 */
boolean replace_is_global_PotentialGene(PotentialGene * obj,boolean is_global) 
{
    if( obj == NULL)     {  
      warn("In replacement function is_global for object PotentialGene, got a NULL object"); 
      return FALSE;  
      }  
    obj->is_global = is_global;  
    return TRUE; 
}    


/* Function:  access_is_global_PotentialGene(obj)
 *
 * Descrip:    Access member variable is_global
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [PotentialGene *]
 *
 * Return [SOFT ]  member variable is_global [boolean]
 *
 */
boolean access_is_global_PotentialGene(PotentialGene * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function is_global for object PotentialGene, got a NULL object");    
      return FALSE;  
      }  
    return obj->is_global;   
}    


/* Function:  replace_name_PotentialGene(obj,name)
 *
 * Descrip:    Replace member variable name
 *             For use principly by API functions
 *
 *
 * Arg:         obj [UNKN ] Object holding the variable [PotentialGene *]
 * Arg:        name [OWNER] New value of the variable [char *]
 *
 * Return [SOFT ]  member variable name [boolean]
 *
 */
boolean replace_name_PotentialGene(PotentialGene * obj,char * name) 
{
    if( obj == NULL)     {  
      warn("In replacement function name for object PotentialGene, got a NULL object");  
      return FALSE;  
      }  
    obj->name = name;    
    return TRUE; 
}    


/* Function:  access_name_PotentialGene(obj)
 *
 * Descrip:    Access member variable name
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [PotentialGene *]
 *
 * Return [SOFT ]  member variable name [char *]
 *
 */
char * access_name_PotentialGene(PotentialGene * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function name for object PotentialGene, got a NULL object"); 
      return NULL;   
      }  
    return obj->name;    
}    


/* Function:  replace_homolog_PotentialGene(obj,homolog)
 *
 * Descrip:    Replace member variable homolog
 *             For use principly by API functions
 *
 *
 * Arg:            obj [UNKN ] Object holding the variable [PotentialGene *]
 * Arg:        homolog [OWNER] New value of the variable [Protein *]
 *
 * Return [SOFT ]  member variable homolog [boolean]
 *
 */
boolean replace_homolog_PotentialGene(PotentialGene * obj,Protein * homolog) 
{
    if( obj == NULL)     {  
      warn("In replacement function homolog for object PotentialGene, got a NULL object");   
      return FALSE;  
      }  
    obj->homolog = homolog;  
    return TRUE; 
}    


/* Function:  access_homolog_PotentialGene(obj)
 *
 * Descrip:    Access member variable homolog
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [PotentialGene *]
 *
 * Return [SOFT ]  member variable homolog [Protein *]
 *
 */
Protein * access_homolog_PotentialGene(PotentialGene * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function homolog for object PotentialGene, got a NULL object");  
      return NULL;   
      }  
    return obj->homolog;     
}    


/* Function:  replace_tsm_PotentialGene(obj,tsm)
 *
 * Descrip:    Replace member variable tsm
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [PotentialGene *]
 * Arg:        tsm [OWNER] New value of the variable [ThreeStateModel *]
 *
 * Return [SOFT ]  member variable tsm [boolean]
 *
 */
boolean replace_tsm_PotentialGene(PotentialGene * obj,ThreeStateModel * tsm) 
{
    if( obj == NULL)     {  
      warn("In replacement function tsm for object PotentialGene, got a NULL object");   
      return FALSE;  
      }  
    obj->tsm = tsm;  
    return TRUE; 
}    


/* Function:  access_tsm_PotentialGene(obj)
 *
 * Descrip:    Access member variable tsm
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [PotentialGene *]
 *
 * Return [SOFT ]  member variable tsm [ThreeStateModel *]
 *
 */
ThreeStateModel * access_tsm_PotentialGene(PotentialGene * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function tsm for object PotentialGene, got a NULL object");  
      return NULL;   
      }  
    return obj->tsm;     
}    


/* Function:  replace_alb_PotentialGene(obj,alb)
 *
 * Descrip:    Replace member variable alb
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [PotentialGene *]
 * Arg:        alb [OWNER] New value of the variable [AlnBlock *]
 *
 * Return [SOFT ]  member variable alb [boolean]
 *
 */
boolean replace_alb_PotentialGene(PotentialGene * obj,AlnBlock * alb) 
{
    if( obj == NULL)     {  
      warn("In replacement function alb for object PotentialGene, got a NULL object");   
      return FALSE;  
      }  
    obj->alb = alb;  
    return TRUE; 
}    


/* Function:  access_alb_PotentialGene(obj)
 *
 * Descrip:    Access member variable alb
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [PotentialGene *]
 *
 * Return [SOFT ]  member variable alb [AlnBlock *]
 *
 */
AlnBlock * access_alb_PotentialGene(PotentialGene * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function alb for object PotentialGene, got a NULL object");  
      return NULL;   
      }  
    return obj->alb;     
}    


/* Function:  replace_bitscore_PotentialGene(obj,bitscore)
 *
 * Descrip:    Replace member variable bitscore
 *             For use principly by API functions
 *
 *
 * Arg:             obj [UNKN ] Object holding the variable [PotentialGene *]
 * Arg:        bitscore [OWNER] New value of the variable [double]
 *
 * Return [SOFT ]  member variable bitscore [boolean]
 *
 */
boolean replace_bitscore_PotentialGene(PotentialGene * obj,double bitscore) 
{
    if( obj == NULL)     {  
      warn("In replacement function bitscore for object PotentialGene, got a NULL object");  
      return FALSE;  
      }  
    obj->bitscore = bitscore;    
    return TRUE; 
}    


/* Function:  access_bitscore_PotentialGene(obj)
 *
 * Descrip:    Access member variable bitscore
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [PotentialGene *]
 *
 * Return [SOFT ]  member variable bitscore [double]
 *
 */
double access_bitscore_PotentialGene(PotentialGene * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function bitscore for object PotentialGene, got a NULL object"); 
      return 0;  
      }  
    return obj->bitscore;    
}    


/* Function:  replace_slop_query_PotentialGene(obj,slop_query)
 *
 * Descrip:    Replace member variable slop_query
 *             For use principly by API functions
 *
 *
 * Arg:               obj [UNKN ] Object holding the variable [PotentialGene *]
 * Arg:        slop_query [OWNER] New value of the variable [int]
 *
 * Return [SOFT ]  member variable slop_query [boolean]
 *
 */
boolean replace_slop_query_PotentialGene(PotentialGene * obj,int slop_query) 
{
    if( obj == NULL)     {  
      warn("In replacement function slop_query for object PotentialGene, got a NULL object");    
      return FALSE;  
      }  
    obj->slop_query = slop_query;    
    return TRUE; 
}    


/* Function:  access_slop_query_PotentialGene(obj)
 *
 * Descrip:    Access member variable slop_query
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [PotentialGene *]
 *
 * Return [SOFT ]  member variable slop_query [int]
 *
 */
int access_slop_query_PotentialGene(PotentialGene * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function slop_query for object PotentialGene, got a NULL object");   
      return 0;  
      }  
    return obj->slop_query;  
}    


/* Function:  replace_slop_target_PotentialGene(obj,slop_target)
 *
 * Descrip:    Replace member variable slop_target
 *             For use principly by API functions
 *
 *
 * Arg:                obj [UNKN ] Object holding the variable [PotentialGene *]
 * Arg:        slop_target [OWNER] New value of the variable [int]
 *
 * Return [SOFT ]  member variable slop_target [boolean]
 *
 */
boolean replace_slop_target_PotentialGene(PotentialGene * obj,int slop_target) 
{
    if( obj == NULL)     {  
      warn("In replacement function slop_target for object PotentialGene, got a NULL object");   
      return FALSE;  
      }  
    obj->slop_target = slop_target;  
    return TRUE; 
}    


/* Function:  access_slop_target_PotentialGene(obj)
 *
 * Descrip:    Access member variable slop_target
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [PotentialGene *]
 *
 * Return [SOFT ]  member variable slop_target [int]
 *
 */
int access_slop_target_PotentialGene(PotentialGene * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function slop_target for object PotentialGene, got a NULL object");  
      return 0;  
      }  
    return obj->slop_target;     
}    


/* Function:  replace_query_length_PotentialGene(obj,query_length)
 *
 * Descrip:    Replace member variable query_length
 *             For use principly by API functions
 *
 *
 * Arg:                 obj [UNKN ] Object holding the variable [PotentialGene *]
 * Arg:        query_length [OWNER] New value of the variable [int]
 *
 * Return [SOFT ]  member variable query_length [boolean]
 *
 */
boolean replace_query_length_PotentialGene(PotentialGene * obj,int query_length) 
{
    if( obj == NULL)     {  
      warn("In replacement function query_length for object PotentialGene, got a NULL object");  
      return FALSE;  
      }  
    obj->query_length = query_length;    
    return TRUE; 
}    


/* Function:  access_query_length_PotentialGene(obj)
 *
 * Descrip:    Access member variable query_length
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [PotentialGene *]
 *
 * Return [SOFT ]  member variable query_length [int]
 *
 */
int access_query_length_PotentialGene(PotentialGene * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function query_length for object PotentialGene, got a NULL object"); 
      return 0;  
      }  
    return obj->query_length;    
}    


/* Function:  access_pex_PotentialTranscript(obj,i)
 *
 * Descrip:    Access members stored in the pex list
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the list [PotentialTranscript *]
 * Arg:          i [UNKN ] Position in the list [int]
 *
 * Return [SOFT ]  Element of the list [PotentialExon *]
 *
 */
PotentialExon * access_pex_PotentialTranscript(PotentialTranscript * obj,int i) 
{
    if( obj == NULL)     {  
      warn("In accessor function pex for object PotentialTranscript, got a NULL object");    
      return NULL;   
      }  
    if( obj->len <= i )  {  
      warn("In accessor function pex for object PotentialTranscript, index %%d is greater than list length %%d",i,obj->len); 
      return NULL;   
      }  
    return obj->pex[i];  
}    


/* Function:  length_pex_PotentialTranscript(obj)
 *
 * Descrip:    discover the length of the list
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the list [PotentialTranscript *]
 *
 * Return [UNKN ]  length of the list [int]
 *
 */
int length_pex_PotentialTranscript(PotentialTranscript * obj) 
{
    if( obj == NULL)     {  
      warn("In length function pex for object PotentialTranscript, got a NULL object");  
      return -1;     
      }  
    return obj->len;     
}    


/* Function:  replace_tstart_PotentialExon(obj,tstart)
 *
 * Descrip:    Replace member variable tstart
 *             For use principly by API functions
 *
 *
 * Arg:           obj [UNKN ] Object holding the variable [PotentialExon *]
 * Arg:        tstart [OWNER] New value of the variable [int]
 *
 * Return [SOFT ]  member variable tstart [boolean]
 *
 */
boolean replace_tstart_PotentialExon(PotentialExon * obj,int tstart) 
{
    if( obj == NULL)     {  
      warn("In replacement function tstart for object PotentialExon, got a NULL object");    
      return FALSE;  
      }  
    obj->tstart = tstart;    
    return TRUE; 
}    


/* Function:  access_tstart_PotentialExon(obj)
 *
 * Descrip:    Access member variable tstart
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [PotentialExon *]
 *
 * Return [SOFT ]  member variable tstart [int]
 *
 */
int access_tstart_PotentialExon(PotentialExon * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function tstart for object PotentialExon, got a NULL object");   
      return 0;  
      }  
    return obj->tstart;  
}    


/* Function:  replace_tend_PotentialExon(obj,tend)
 *
 * Descrip:    Replace member variable tend
 *             For use principly by API functions
 *
 *
 * Arg:         obj [UNKN ] Object holding the variable [PotentialExon *]
 * Arg:        tend [OWNER] New value of the variable [int]
 *
 * Return [SOFT ]  member variable tend [boolean]
 *
 */
boolean replace_tend_PotentialExon(PotentialExon * obj,int tend) 
{
    if( obj == NULL)     {  
      warn("In replacement function tend for object PotentialExon, got a NULL object");  
      return FALSE;  
      }  
    obj->tend = tend;    
    return TRUE; 
}    


/* Function:  access_tend_PotentialExon(obj)
 *
 * Descrip:    Access member variable tend
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [PotentialExon *]
 *
 * Return [SOFT ]  member variable tend [int]
 *
 */
int access_tend_PotentialExon(PotentialExon * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function tend for object PotentialExon, got a NULL object"); 
      return 0;  
      }  
    return obj->tend;    
}    


/* Function:  replace_qstart_PotentialExon(obj,qstart)
 *
 * Descrip:    Replace member variable qstart
 *             For use principly by API functions
 *
 *
 * Arg:           obj [UNKN ] Object holding the variable [PotentialExon *]
 * Arg:        qstart [OWNER] New value of the variable [int]
 *
 * Return [SOFT ]  member variable qstart [boolean]
 *
 */
boolean replace_qstart_PotentialExon(PotentialExon * obj,int qstart) 
{
    if( obj == NULL)     {  
      warn("In replacement function qstart for object PotentialExon, got a NULL object");    
      return FALSE;  
      }  
    obj->qstart = qstart;    
    return TRUE; 
}    


/* Function:  access_qstart_PotentialExon(obj)
 *
 * Descrip:    Access member variable qstart
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [PotentialExon *]
 *
 * Return [SOFT ]  member variable qstart [int]
 *
 */
int access_qstart_PotentialExon(PotentialExon * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function qstart for object PotentialExon, got a NULL object");   
      return 0;  
      }  
    return obj->qstart;  
}    


/* Function:  replace_qend_PotentialExon(obj,qend)
 *
 * Descrip:    Replace member variable qend
 *             For use principly by API functions
 *
 *
 * Arg:         obj [UNKN ] Object holding the variable [PotentialExon *]
 * Arg:        qend [OWNER] New value of the variable [int]
 *
 * Return [SOFT ]  member variable qend [boolean]
 *
 */
boolean replace_qend_PotentialExon(PotentialExon * obj,int qend) 
{
    if( obj == NULL)     {  
      warn("In replacement function qend for object PotentialExon, got a NULL object");  
      return FALSE;  
      }  
    obj->qend = qend;    
    return TRUE; 
}    


/* Function:  access_qend_PotentialExon(obj)
 *
 * Descrip:    Access member variable qend
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [PotentialExon *]
 *
 * Return [SOFT ]  member variable qend [int]
 *
 */
int access_qend_PotentialExon(PotentialExon * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function qend for object PotentialExon, got a NULL object"); 
      return 0;  
      }  
    return obj->qend;    
}    



#ifdef _cplusplus
}
#endif
