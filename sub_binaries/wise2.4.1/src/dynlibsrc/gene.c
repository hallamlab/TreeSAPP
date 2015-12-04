#ifdef _cplusplus
extern "C" {
#endif
#include "gene.h"


/* Function:  reversed_Gene(g)
 *
 * Descrip:    is this gene reversed?
 *
 *
 * Arg:        g [UNKN ] Undocumented argument [Gene *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 67 "gene.dy"
boolean reversed_Gene(Gene * g)
{
  if( g->start < g->end ) 
    return FALSE;
  return TRUE;
}

/* Function:  copy_Gene(g)
 *
 * Descrip:    Makes a completely fresh copy of a
 *             gene
 *
 *
 * Arg:        g [UNKN ] Undocumented argument [Gene *]
 *
 * Return [UNKN ]  Undocumented return value [Gene *]
 *
 */
# line 78 "gene.dy"
Gene * copy_Gene(Gene * g)
{
  int i;
  Gene * out;

  out = Gene_alloc();
  for(i=0;i<g->len;i++)
    add_Gene(out,copy_Transcript(g->transcript[i]));

  return out;
}

/* Function:  is_simple_prediction_Gene(g)
 *
 * Descrip:    Does this gene have 
 *             	a single transcript
 *             	that transcript with translation start/end 
 *             	at the ends
 *
 *
 * Arg:        g [UNKN ] Undocumented argument [Gene *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 96 "gene.dy"
boolean is_simple_prediction_Gene(Gene * g)
{
  if( g->len > 1 ) 
    return FALSE;
  if( g->transcript[0]->len > 1 ) 
    return FALSE;
  if( g->transcript[0]->translation[0]->start != 0 || g->transcript[0]->translation[0]->end != length_Transcript(g->transcript[0]) )
    return FALSE;
  
  return TRUE;
}


      
/* Function:  get_Genomic_from_Gene(gene)
 *
 * Descrip:    Gives back a Genomic sequence type
 *             from a gene.
 *
 *
 * Arg:        gene [READ ] gene to get Genomic from [Gene *]
 *
 * Return [SOFT ]  Genomic DNA data structure [Genomic *]
 *
 */
# line 117 "gene.dy"
Genomic * get_Genomic_from_Gene(Gene * gene)
{
  Genomic * gn;
  char buffer[64];

  /*  fprintf(stdout,"Getting genomic...\n"); */
  if( gene->genomic != NULL )
    return gene->genomic;

  if( gene->parent == NULL ) {
    warn("Cannot get Gene, as no parent genomic region!");
    return NULL;
  }


  gn = get_Genomic_from_GenomicRegion(gene->parent);

  if( gn == NULL) {
    warn("Cannot get Gene, as no sequence in genomic region!");
    return NULL;
  }

  if( gn->baseseq->offset < gn->baseseq->end) {
    if( gene->start > gene->end ) 
      gene->genomic = truncate_Genomic(gn,gene->start-gn->baseseq->offset+2,gene->end-gn->baseseq->offset+2);
    else 
      gene->genomic = truncate_Genomic(gn,gene->start-gn->baseseq->offset+1,gene->end-gn->baseseq->offset+1);
  }
  else {
    gene->genomic = truncate_Genomic(gn,gn->baseseq->offset-1 - gene->start,gn->baseseq->offset-1 - gene->end);
  }

  sprintf(buffer,"%s.[%d:%d]",Genomic_name(gn),gene->start+1,gene->end);
  ckfree(gene->genomic->baseseq->name);
  gene->genomic->baseseq->name = stringalloc(buffer);


  return gene->genomic;
}

#define MAX_EMBL_EXON_PARSE 128
/* Function:  read_EMBL_feature_Gene(buffer,maxlen,ifp)
 *
 * Descrip:    Reads in an EMBL feature table.
 *
 *             It expects to be passed a buffer with 'FT   CDS'.
 *             or 'FT   mRNA' in it. It will then 
 *             use the buffer to read successive lines of the Feature table
 *             until it comes to the next 'meta' feature line (ie, 3 place point).
 *
 *             It will use functions in /embl module for the reading.
 *
 *
 * Arg:        buffer [UNKN ] a string with FT  CDS line in it [char *]
 * Arg:        maxlen [UNKN ] length of the buffer [int]
 * Arg:           ifp [UNKN ] file stream with the rest of the feature table in it [FILE *]
 *
 * Return [UNKN ]  Undocumented return value [Gene *]
 *
 */
# line 172 "gene.dy"
Gene * read_EMBL_feature_Gene(char * buffer,int maxlen,FILE * ifp)
{
  Gene * gene;
  Transcript * tr;
  Translation * ts;
  Exon * exon;

  char * runner;
  char * base;
  char * next;
  int i;
  int exon_start[MAX_EMBL_EXON_PARSE];
  int exon_end[MAX_EMBL_EXON_PARSE];
  int number;
  int exon_no = 0;
  int isstart = 1;
  int is_complement = 0;
  int is_cds = 0;
  int break_at_end = 0;

  if( strstartcmp(buffer,"FT") != 0 ) {
    warn("passed in a bad line [%s] to be used for feature table parsing",buffer);
    return NULL;
  }

  if( (runner=strtok(buffer+2,spacestr)) == NULL ) {
    warn("Bad embl feature line [%s]",buffer);
    return NULL;
  }

  if( strcmp(runner,"CDS") != 0 && strcmp(runner,"mRNA") != 0 ) {
    warn("passed in a feature line to read_EMBL_feature_Gene with a %s tag. This only handles CDS and mRNA tags",runner);
    return NULL;
  }

  if( strcmp(runner,"CDS") == 0 ) {
    is_cds = TRUE;
  }

  runner = strtok(NULL,spacestr);

  if( runner == NULL ) {
    warn("Bad embl feature line [%s]",buffer);
    return NULL;
  }

  if( strstartcmp(runner,"complement") == 0 ) {
    runner = strchr(runner,'(');
    if( runner == NULL) {
      warn("Could not find bracket on EMBL feature complement line");
      return NULL;
    }
    is_complement = 1;
    runner++;
  }


  if( strstartcmp(runner,"join") == 0 ) {
    runner = strchr(runner,'(');
    runner++;
  } else if( isdigit((int)*runner)  || *runner == '<' ) {
    /** ok - starts with the numbers. We'll cope!**/
  } else {
    warn("Expecting a join statement, got a [%s]",runner);
    return NULL;
  }

  
  /*** ok, now the major number loop ***/

  for(;;) {
    base= runner;
    for(;*runner && *runner != ')' && *runner != '.' && *runner != ',' && *runner != '>' && !isspace((int)*runner);runner++) 
      ;

    /*fprintf(stderr,"Got a runner of %s\n ",runner); */
    if( *runner == '\0' )
      next = runner;
    else next = runner+1;

    if( *runner == ')' ) {
      break_at_end = TRUE; /* out of reading exons */
    }

    
    *runner='\0';
    if( strstartcmp(base,"complement(") == 0 ) {
      is_complement = TRUE;
      for(;*base != '(';base++) 
	;
      base++;
      break_at_end = FALSE; /* we found an bracket too early! */
    }

    if( is_integer_string(base,&number) == FALSE ) {
      warn("Got a non integer [%s] in the middle of a join statement in EMBL parsing",runner);
      return NULL;
    }

    /** put this number away **/

    if( isstart ) {
      exon_start[exon_no] = number;
      isstart = 0;
    } else {
      exon_end[exon_no++] = number;
      isstart = 1;
    }
    if( break_at_end == TRUE)
      break;

    for(runner=next;*runner && (*runner == '.' || isspace((int)*runner));runner++)
      ;

    if( *runner == '\0' ) {
      if( next_feature_tab_line(buffer,maxlen,ifp) == FALSE) {
	warn("In the middle of getting a join statement, got a [%s]. Yuk!",buffer);
	return NULL;
      }

      if( !isdigit((int)buffer[0]) && buffer[0] != '.' && buffer[0] != ',') {
	/*** ok - sometimes people very boring end things in here ***/
	/* warn("In the middle of getting a join statement, got a [%s]. Ugh!",buffer); */
	break;
      }

      runner = buffer;
    }

  }

  if( isstart == 0 ) {
    warn("I have read an uneven number of start-end points in the exon thing. Yuk!");
    return NULL;
  }

  /** runner should now be on bracket **/

  if( is_complement == 1 ) {
    /** ok . should be another bracket. Do we care? **/
  }

  gene = Gene_alloc_len(1);
  tr  = Transcript_alloc_len(exon_no);
  add_Gene(gene,tr);
  tr->parent = gene;

  if( is_complement == 1 ) {
    gene->start = exon_end[exon_no-1]-1;
    gene->end = exon_start[0] -1;

    for(i=exon_no -1;i >= 0;i--) {
      exon = Exon_alloc();
      exon->start = (gene->start+1) - exon_end[i];
      exon->end = (gene->start+1) - exon_start[i] +1;
      add_ex_Transcript(tr,exon);
    }
  } else {
    gene->start = exon_start[0] -1;
    gene->end = exon_end[exon_no-1] -1;

    for(i=0;i<exon_no;i++) {
      exon = Exon_alloc();
      exon->start = exon_start[i] - (gene->start+1);
      exon->end = exon_end[i] - (gene->start+1)+1;
      add_ex_Transcript(tr,exon);
    }
  }

  if( is_cds == TRUE ) {
    ts = Translation_alloc();
    ts->start = 0;
    ts->end = length_Transcript(tr);
    ts->parent = tr;
    add_Transcript(tr,ts);
  }

  /*** read the rest of this feature ***/

  while( next_feature_tab_line(buffer,maxlen,ifp) == TRUE)
    ;

  return gene;

}

/* Function:  write_Embl_FT_Gene(ge,key,ofp)
 *
 * Descrip:    shows a embl feature table part
 *
 *
 * Arg:         ge [UNKN ] Undocumented argument [Gene *]
 * Arg:        key [UNKN ] Undocumented argument [char *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 361 "gene.dy"
void write_Embl_FT_Gene(Gene * ge,char * key,FILE * ofp)
{
  int i;
  Transcript * tr;
  int j;


  if( ge->start > ge->end ) {
    for(i=0;i<ge->len;i++) {
      tr = ge->transcript[i];
      if( tr->ex_len > 1 ) {
	fprintf(ofp,"FT   %10s    complement(join(",key);
      } else {
	fprintf(ofp,"FT   %10s    complement(",key);
      }
      for(j=0;j<tr->ex_len;j++) {
	fprintf(ofp,"%d..%d%s",ge->start+1 - tr->exon[j]->start,ge->start - tr->exon[j]->end+2,j+1 == tr->ex_len ? "" : ",");
	if( j != 0 && j+1 != tr->ex_len && (j%3) == 0 ) {
	  fprintf(ofp,"\nFT                    ");
	}
	
      } /* end of for over exons */
      if( tr->ex_len > 1 ) {
	fprintf(ofp,"))\n");
      } else {
	fprintf(ofp,")\n");
      }
    }
  } else {
    for(i=0;i<ge->len;i++) {
      tr = ge->transcript[i];
      if( tr->ex_len > 1 ) {
	fprintf(ofp,"FT   %10s    join(",key);
      } else {
	fprintf(ofp,"FT   %10s    ",key);
      }
      for(j=0;j<tr->ex_len;j++) {
	  fprintf(ofp,"%d..%d%s",tr->exon[j]->start+ge->start+1,tr->exon[j]->end+ge->start,j+1 == tr->ex_len ? "" : ",");
	  if( j != 0 && j+1 != tr->ex_len && (j%3) == 0 ) {
	    fprintf(ofp,"\nFT                   ");
	  }
      }
    }
    if( tr->ex_len > 1 ) {
      fprintf(ofp,")\n");
    } else {
      fprintf(ofp,"\n");
    }
  }

}
  
/* Function:  show_pretty_Gene(ge,show_supporting,ofp)
 *
 * Descrip:    Shows a gene in the biologically accepted form
 *
 *
 * Arg:                     ge [UNKN ] Undocumented argument [Gene *]
 * Arg:        show_supporting [UNKN ] Undocumented argument [boolean]
 * Arg:                    ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 416 "gene.dy"
void show_pretty_Gene(Gene * ge,boolean show_supporting,FILE * ofp)
{
  int i;
  int j;
  int k;

  if( ge->start > ge->end ) {
    fprintf(ofp,"Gene %d %d%s\n",ge->start+1,ge->end+2,ge->ispseudo == TRUE ? " [pseudogene]" : "");
    for(i=0;i<ge->len;i++) {
      auto Transcript * tr;
      if( ge->len != 1 ) {
	fprintf(ofp," Transcript %d\n",i);
      }
      tr = ge->transcript[i];
      for(j=0;j<tr->ex_len;j++) {
	fprintf(ofp,"  Exon %d %d phase %d\n",ge->start+1 - tr->exon[j]->start,ge->start - tr->exon[j]->end+2,tr->exon[j]->phase);
	if( show_supporting ) {
	  auto Exon * ex;
	  ex = tr->exon[j];
	  for(k=0;k<tr->exon[j]->len;k++) {
	    fprintf(ofp,"     Supporting %d %d %d %d\n",ge->start+1 - ex->sf[k]->start,ge->start - ex->sf[k]->end+2,ex->sf[k]->hstart+1,ex->sf[k]->hend);
	  }
	}
      }
    }
  } else {
    fprintf(ofp,"Gene %d %d %s\n",ge->start+1,ge->end,ge->ispseudo == TRUE ? "[pseudogene]" : "");
    for(i=0;i<ge->len;i++) {
      auto Transcript * tr;
      if( ge->len != 1 ) {
	fprintf(ofp," Transcript %d\n",i);
      }
      tr = ge->transcript[i];
      for(j=0;j<tr->ex_len;j++) {
	fprintf(ofp,"  Exon %d %d phase %d\n",tr->exon[j]->start+ge->start+1,
		tr->exon[j]->end+ge->start,tr->exon[j]->phase);

	if( show_supporting ) {
	  auto Exon * ex;
	  ex = tr->exon[j];
	  for(k=0;k<tr->exon[j]->len;k++) {
	    fprintf(ofp,"     Supporting %d %d %d %d\n",ge->start+1+ex->sf[k]->start,ge->start+ex->sf[k]->end,ex->sf[k]->hstart+1,ex->sf[k]->hend);
	  }
	}
      }
    }
  } 

}

/* Function:  show_Gene(ge,ofp)
 *
 * Descrip:    shows a gene in a vaguely human readable form
 *
 *
 * Arg:         ge [UNKN ] Undocumented argument [Gene *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 469 "gene.dy"
void show_Gene(Gene * ge,FILE * ofp)
{
  int i;

  fprintf(ofp,"Gene %d - %d\n",ge->start,ge->end);
  for(i=0;i<ge->len;i++) {
    fprintf(ofp,"Transcript %d\n",i);
    show_Transcript(ge->transcript[i],ofp);
  }

}

# line 465 "gene.c"
/* Function:  swap_Gene(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_Gene
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [Transcript **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_Gene(Transcript ** list,int i,int j)  
{
    Transcript * temp;   
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_Gene(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_Gene which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [Transcript **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_Gene(Transcript ** list,int left,int right,int (*comp)(Transcript * ,Transcript * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_Gene(list,left,(left+right)/2); 
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_Gene (list,++last,i);   
      }  
    swap_Gene (list,left,last);  
    qsort_Gene(list,left,last-1,comp);   
    qsort_Gene(list,last+1,right,comp);  
}    


/* Function:  sort_Gene(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_Gene
 *
 *
 * Arg:         obj [UNKN ] Object containing list [Gene *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_Gene(Gene * obj,int (*comp)(Transcript *, Transcript *)) 
{
    qsort_Gene(obj->transcript,0,obj->len-1,comp);   
    return;  
}    


/* Function:  expand_Gene(obj,len)
 *
 * Descrip:    Really an internal function for add_Gene
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [Gene *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_Gene(Gene * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_Gene called with no need");   
      return TRUE;   
      }  


    if( (obj->transcript = (Transcript ** ) ckrealloc (obj->transcript,sizeof(Transcript *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_Gene, returning FALSE"); 
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_Gene(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [Gene *]
 * Arg:        add [OWNER] Object to add to the list [Transcript *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_Gene(Gene * obj,Transcript * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_Gene(obj,obj->len + GeneLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->transcript[obj->len++]=add; 
    return TRUE; 
}    


/* Function:  flush_Gene(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [Gene *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_Gene(Gene * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->transcript[i] != NULL)    {  
        free_Transcript(obj->transcript[i]); 
        obj->transcript[i] = NULL;   
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  Gene_alloc_std(void)
 *
 * Descrip:    Equivalent to Gene_alloc_len(GeneLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Gene *]
 *
 */
Gene * Gene_alloc_std(void) 
{
    return Gene_alloc_len(GeneLISTLENGTH);   
}    


/* Function:  Gene_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [Gene *]
 *
 */
Gene * Gene_alloc_len(int len) 
{
    Gene * out; /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = Gene_alloc()) == NULL) 
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->transcript = (Transcript ** ) ckcalloc (len,sizeof(Transcript *))) == NULL) {  
      warn("Warning, ckcalloc failed in Gene_alloc_len");    
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_Gene(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [Gene *]
 *
 * Return [UNKN ]  Undocumented return value [Gene *]
 *
 */
Gene * hard_link_Gene(Gene * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a Gene object: passed a NULL object");    
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  Gene_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Gene *]
 *
 */
Gene * Gene_alloc(void) 
{
    Gene * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(Gene *) ckalloc (sizeof(Gene))) == NULL)    {  
      warn("Gene_alloc failed ");    
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->start = 0;  
    out->end = 0;    
    out->genomic = NULL; 
    out->transcript = NULL;  
    out->len = out->maxlen = 0;  
    out->name = NULL;    
    out->bits = 0;   
    out->seqname = NULL; 
    out->ispseudo = FALSE;   


    return out;  
}    


/* Function:  free_Gene(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [Gene *]
 *
 * Return [UNKN ]  Undocumented return value [Gene *]
 *
 */
Gene * free_Gene(Gene * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a Gene obj. Should be trappable");  
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
    /* obj->parent is linked in */ 
    if( obj->genomic != NULL)    
      free_Genomic(obj->genomic);    
    if( obj->transcript != NULL) {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->transcript[i] != NULL)  
          free_Transcript(obj->transcript[i]);   
        }  
      ckfree(obj->transcript);   
      }  
    if( obj->name != NULL)   
      ckfree(obj->name);     
    if( obj->seqname != NULL)    
      ckfree(obj->seqname);  


    ckfree(obj); 
    return NULL; 
}    


/* Function:  replace_start_Gene(obj,start)
 *
 * Descrip:    Replace member variable start
 *             For use principly by API functions
 *
 *
 * Arg:          obj [UNKN ] Object holding the variable [Gene *]
 * Arg:        start [OWNER] New value of the variable [int]
 *
 * Return [SOFT ]  member variable start [boolean]
 *
 */
boolean replace_start_Gene(Gene * obj,int start) 
{
    if( obj == NULL)     {  
      warn("In replacement function start for object Gene, got a NULL object");  
      return FALSE;  
      }  
    obj->start = start;  
    return TRUE; 
}    


/* Function:  access_start_Gene(obj)
 *
 * Descrip:    Access member variable start
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [Gene *]
 *
 * Return [SOFT ]  member variable start [int]
 *
 */
int access_start_Gene(Gene * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function start for object Gene, got a NULL object"); 
      return 0;  
      }  
    return obj->start;   
}    


/* Function:  replace_end_Gene(obj,end)
 *
 * Descrip:    Replace member variable end
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [Gene *]
 * Arg:        end [OWNER] New value of the variable [int]
 *
 * Return [SOFT ]  member variable end [boolean]
 *
 */
boolean replace_end_Gene(Gene * obj,int end) 
{
    if( obj == NULL)     {  
      warn("In replacement function end for object Gene, got a NULL object");    
      return FALSE;  
      }  
    obj->end = end;  
    return TRUE; 
}    


/* Function:  access_end_Gene(obj)
 *
 * Descrip:    Access member variable end
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [Gene *]
 *
 * Return [SOFT ]  member variable end [int]
 *
 */
int access_end_Gene(Gene * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function end for object Gene, got a NULL object");   
      return 0;  
      }  
    return obj->end;     
}    


/* Function:  replace_parent_Gene(obj,parent)
 *
 * Descrip:    Replace member variable parent
 *             For use principly by API functions
 *
 *
 * Arg:           obj [UNKN ] Object holding the variable [Gene *]
 * Arg:        parent [OWNER] New value of the variable [GenomicRegion *]
 *
 * Return [SOFT ]  member variable parent [boolean]
 *
 */
boolean replace_parent_Gene(Gene * obj,GenomicRegion * parent) 
{
    if( obj == NULL)     {  
      warn("In replacement function parent for object Gene, got a NULL object"); 
      return FALSE;  
      }  
    obj->parent = parent;    
    return TRUE; 
}    


/* Function:  access_parent_Gene(obj)
 *
 * Descrip:    Access member variable parent
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [Gene *]
 *
 * Return [SOFT ]  member variable parent [GenomicRegion *]
 *
 */
GenomicRegion * access_parent_Gene(Gene * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function parent for object Gene, got a NULL object");    
      return NULL;   
      }  
    return obj->parent;  
}    


/* Function:  replace_genomic_Gene(obj,genomic)
 *
 * Descrip:    Replace member variable genomic
 *             For use principly by API functions
 *
 *
 * Arg:            obj [UNKN ] Object holding the variable [Gene *]
 * Arg:        genomic [OWNER] New value of the variable [Genomic *]
 *
 * Return [SOFT ]  member variable genomic [boolean]
 *
 */
boolean replace_genomic_Gene(Gene * obj,Genomic * genomic) 
{
    if( obj == NULL)     {  
      warn("In replacement function genomic for object Gene, got a NULL object");    
      return FALSE;  
      }  
    obj->genomic = genomic;  
    return TRUE; 
}    


/* Function:  access_genomic_Gene(obj)
 *
 * Descrip:    Access member variable genomic
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [Gene *]
 *
 * Return [SOFT ]  member variable genomic [Genomic *]
 *
 */
Genomic * access_genomic_Gene(Gene * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function genomic for object Gene, got a NULL object");   
      return NULL;   
      }  
    return obj->genomic;     
}    


/* Function:  access_transcript_Gene(obj,i)
 *
 * Descrip:    Access members stored in the transcript list
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the list [Gene *]
 * Arg:          i [UNKN ] Position in the list [int]
 *
 * Return [SOFT ]  Element of the list [Transcript *]
 *
 */
Transcript * access_transcript_Gene(Gene * obj,int i) 
{
    if( obj == NULL)     {  
      warn("In accessor function transcript for object Gene, got a NULL object");    
      return NULL;   
      }  
    if( obj->len <= i )  {  
      warn("In accessor function transcript for object Gene, index %%d is greater than list length %%d",i,obj->len); 
      return NULL;   
      }  
    return obj->transcript[i];   
}    


/* Function:  length_transcript_Gene(obj)
 *
 * Descrip:    discover the length of the list
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the list [Gene *]
 *
 * Return [UNKN ]  length of the list [int]
 *
 */
int length_transcript_Gene(Gene * obj) 
{
    if( obj == NULL)     {  
      warn("In length function transcript for object Gene, got a NULL object");  
      return -1;     
      }  
    return obj->len;     
}    


/* Function:  replace_name_Gene(obj,name)
 *
 * Descrip:    Replace member variable name
 *             For use principly by API functions
 *
 *
 * Arg:         obj [UNKN ] Object holding the variable [Gene *]
 * Arg:        name [OWNER] New value of the variable [char *]
 *
 * Return [SOFT ]  member variable name [boolean]
 *
 */
boolean replace_name_Gene(Gene * obj,char * name) 
{
    if( obj == NULL)     {  
      warn("In replacement function name for object Gene, got a NULL object");   
      return FALSE;  
      }  
    obj->name = name;    
    return TRUE; 
}    


/* Function:  access_name_Gene(obj)
 *
 * Descrip:    Access member variable name
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [Gene *]
 *
 * Return [SOFT ]  member variable name [char *]
 *
 */
char * access_name_Gene(Gene * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function name for object Gene, got a NULL object");  
      return NULL;   
      }  
    return obj->name;    
}    


/* Function:  replace_bits_Gene(obj,bits)
 *
 * Descrip:    Replace member variable bits
 *             For use principly by API functions
 *
 *
 * Arg:         obj [UNKN ] Object holding the variable [Gene *]
 * Arg:        bits [OWNER] New value of the variable [double]
 *
 * Return [SOFT ]  member variable bits [boolean]
 *
 */
boolean replace_bits_Gene(Gene * obj,double bits) 
{
    if( obj == NULL)     {  
      warn("In replacement function bits for object Gene, got a NULL object");   
      return FALSE;  
      }  
    obj->bits = bits;    
    return TRUE; 
}    


/* Function:  access_bits_Gene(obj)
 *
 * Descrip:    Access member variable bits
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [Gene *]
 *
 * Return [SOFT ]  member variable bits [double]
 *
 */
double access_bits_Gene(Gene * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function bits for object Gene, got a NULL object");  
      return 0;  
      }  
    return obj->bits;    
}    


/* Function:  replace_seqname_Gene(obj,seqname)
 *
 * Descrip:    Replace member variable seqname
 *             For use principly by API functions
 *
 *
 * Arg:            obj [UNKN ] Object holding the variable [Gene *]
 * Arg:        seqname [OWNER] New value of the variable [char *]
 *
 * Return [SOFT ]  member variable seqname [boolean]
 *
 */
boolean replace_seqname_Gene(Gene * obj,char * seqname) 
{
    if( obj == NULL)     {  
      warn("In replacement function seqname for object Gene, got a NULL object");    
      return FALSE;  
      }  
    obj->seqname = seqname;  
    return TRUE; 
}    


/* Function:  access_seqname_Gene(obj)
 *
 * Descrip:    Access member variable seqname
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [Gene *]
 *
 * Return [SOFT ]  member variable seqname [char *]
 *
 */
char * access_seqname_Gene(Gene * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function seqname for object Gene, got a NULL object");   
      return NULL;   
      }  
    return obj->seqname;     
}    


/* Function:  replace_ispseudo_Gene(obj,ispseudo)
 *
 * Descrip:    Replace member variable ispseudo
 *             For use principly by API functions
 *
 *
 * Arg:             obj [UNKN ] Object holding the variable [Gene *]
 * Arg:        ispseudo [OWNER] New value of the variable [boolean]
 *
 * Return [SOFT ]  member variable ispseudo [boolean]
 *
 */
boolean replace_ispseudo_Gene(Gene * obj,boolean ispseudo) 
{
    if( obj == NULL)     {  
      warn("In replacement function ispseudo for object Gene, got a NULL object");   
      return FALSE;  
      }  
    obj->ispseudo = ispseudo;    
    return TRUE; 
}    


/* Function:  access_ispseudo_Gene(obj)
 *
 * Descrip:    Access member variable ispseudo
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [Gene *]
 *
 * Return [SOFT ]  member variable ispseudo [boolean]
 *
 */
boolean access_ispseudo_Gene(Gene * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function ispseudo for object Gene, got a NULL object");  
      return FALSE;  
      }  
    return obj->ispseudo;    
}    



#ifdef _cplusplus
}
#endif
