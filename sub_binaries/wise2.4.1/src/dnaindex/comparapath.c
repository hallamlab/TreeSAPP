#ifdef _cplusplus
extern "C" {
#endif
#include "comparapath.h"


# line 88 "comparapath.dy"
ComparaLinkStart * free_ComparaLinkStart(ComparaLinkStart * cls)
{
  free(cls);
  return NULL;
}


# line 95 "comparapath.dy"
void show_SetofHSPset(SetofHSPset * set,FILE * ofp)
{
  int i;

  for(i=0;i<set->len;i++) {
    show_HSPset(set->hspset[i],ofp); 
  }

}


# line 106 "comparapath.dy"
SetofHSPset * SetofHSPset_from_ComparaIndex(ComparaIndex * ci,ComparaLinkStartSet * clss,FILE * logfp)
{
  int i;
  HSPset * hsp;
  SetofHSPset * out;

  assert(ci != NULL);
  assert(clss != NULL);

  out = SetofHSPset_alloc_std();

  for(i=0;i<clss->len;i++) {
    if( clss->cls[i]->start != NULL ) {
      fprintf(logfp,"Considering sequence %s\n",clss->cls[i]->seq->name);
      fflush(logfp);

      hsp = HSPset_from_ComparaIndex(ci,clss->cls[i],logfp);
      add_SetofHSPset(out,hsp);
    }
  }

  return out;
}

# line 130 "comparapath.dy"
HSPset * HSPset_from_ComparaIndex(ComparaIndex * ci,ComparaLinkStart * cls,FILE * logfp)
{
  HSPset * out;
  HSP * hsp;
  ComparaHead * current;
  long int target_pos;
  ComparaHead * current_head;
  ComparaHead * prev = NULL;

  SinglePosSequence * query = NULL;
  SinglePosSequence * target = NULL;

  long int start_pos;
  long int next_report_post;
  long int last_pos;
  long int unplaced = 0;

  out = HSPset_alloc_std();

  query = lookup_Sequence_SinglePosSpace(ci->sps,cls->start->position[0]);

  start_pos = cls->start->position[0];
  next_report_post = start_pos;

  for(current = cls->start;current != NULL;) {
    if( current->next_query == current ) {
      fprintf(logfp,"Impossible; loop in next_query list, causing massive problems. have to abort");
      fflush(logfp);

      fatal("loop in next query, cannot contiue");
    }

    last_pos = current->position[0];

    current_head = current;
    if( current_head == prev || (prev != NULL &&  current_head->position[0] == prev->position[0]) ) {
      fprintf(logfp,"Dangerous case of unmoving loop; artificalling moving on from position %ld\n",current->position[0]);
      fflush(logfp);
      current = current->next_query;
      current_head = current;
      if( current == NULL ) {
	break;
      }
    }
    
    prev = current;
      
    if( current->position[0] >= next_report_post ) {
      fprintf(logfp,"Seen %ld unplaced vs %ld positions\n",unplaced,current->position[0] - start_pos);

      fprintf(logfp,"Considering positions %ld,%ld with size %d and spline %ld, (current %ld, prev %ld next %ld)\n",current->position[0],current->position[1],current->size,current->spline,current,prev,current->next_query);
      next_report_post += 10000000;
    }


    if( COMPARA_IS_JOINT_FORWARD(current_head) ) {
      /*
      fprintf(logfp,"Found forward, with state %d (len %d)\n",current_head->state,current_head->size);
      fflush(logfp);
      */
      target_pos  = current_head->position[1];


      if( target == NULL || target->end > target_pos || target->start < target_pos ) {
	target = lookup_Sequence_SinglePosSpace(ci->sps,target_pos);
      }


      assert(target);
      assert(target->data);
/*
      fprintf(logfp,"  ...forward is %s %ld\n",((Sequence*)target->data)->name,target_pos - target->start);
      fflush(logfp);
*/

      hsp = new_dna_identical_HSP((Sequence*)query->data,(Sequence*)target->data,current->position[0] - query->start,target_pos - target->start,0);
      add_HSPset(out,hsp);


      for(;current != NULL && (current->position[0] - query->start) < hsp->query_start + hsp->length;current = current->next_query)
	;

    } else if ( COMPARA_IS_JOINT_REVERSE(current_head) ) {
      /*    
	    fprintf(logfp,"Found reverse, with state %d (len %d)\n",current_head->state,current_head->size);
	    fflush(logfp);
      */  

      target_pos = current_head->spline->position[0];

      if( target == NULL || target->end > target_pos || target->start < target_pos ) {
	target = lookup_Sequence_SinglePosSpace(ci->sps,target_pos);
      }

      /*
      fprintf(logfp,"  ...reverse is %s %ld\n",((Sequence *)target->data)->name,target_pos - target->start);
      fflush(logfp);
      */

      hsp = new_dna_identical_HSP((Sequence*)query->data,(Sequence*)target->data,current->position[0] - query->start,target_pos - target->start +  ci->kii->kmer_size -1,1);


      add_HSPset(out,hsp);

      
      for(;current != NULL && (current->position[0] - query->start) < hsp->query_start + hsp->length;current = current->next_query)
	;

      
    } else {
      unplaced++;
      current = current->next_query;
    }
  }


  fprintf(logfp,"PLACED %s %ld unplaced positions vs %ld considered\n",cls->seq->name,unplaced,last_pos - start_pos);
  fflush(logfp);
  return out;
  
}



# line 254 "comparapath.dy"
boolean is_joint_forward(ComparaHead * h)
{
   if( h->spline != NULL ) {
     return 0;
   } 
   if( (h->state & COMPARA_QUERY_UNIQUE) && (h->state & COMPARA_TARGET_UNIQUE) ) {
     return 1;
   }

   return 0;
}

# line 266 "comparapath.dy"
boolean is_joint_reverse(ComparaHead * h)
{
   if( h->spline == NULL ) {
     return 0; 
   } 
   
   if( (h->state & COMPARA_QUERY_UNIQUE) && (h->spline->state & COMPARA_TARGET_UNIQUE) && 
       !(h->spline->state & COMPARA_QUERY_UNIQUE) &&
       !(h->spline->state & COMPARA_QUERY_MULTIPLE) 
       ) {
     return 1;
   }

   /*
   if( (h->spline->state & COMPARA_QUERY_UNIQUE) && (h->state & COMPARA_TARGET_UNIQUE) &&
       !(h->spline->state & COMPARA_TARGET_UNIQUE) &&
       !(h->spline->state & COMPARA_TARGET_MULTIPLE)
       ) {
     return 1;
   }
   */

   return 0;
}


# line 292 "comparapath.dy"
void show_distrib_ComparaIndex(ComparaIndex * ci,ComparaLinkStart * cls,FILE * ofp)
{
  long len = 0;
  int count[4];
  int i;
  ComparaHead * current;
  ComparaHead * chead;

  for(i=0;i<4;i++) {
    count[i] =0;
  }

  for(current = cls->start;current != NULL;current = current->next_query) {

    chead = current;

    if( chead == NULL ) {
      continue;
    }

    count[(int)chead->size]++;
    len++;
  }

  fprintf(ofp,"Seen %ld links...\n",len);

  for(i=0;i<3;i++) {
    fprintf(ofp,"%5d  %d\n",i,count[i]);
  }

  
}

# line 325 "comparapath.dy"
void show_stats_ComparaIndex(ComparaIndex * ci,ComparaLinkStart * cls,FILE * ofp)
{
  int i;
  int t_unique = 0;
  int q_unique = 0;
  int joint = 0;
  int rev_joint  =0;

  int t_multiple = 0;
  int q_multiple = 0;
  ComparaHead *chead;
  
  ComparaHead * current;


  for(current = cls->start;current != NULL;current = current->next_query) {

    i = current->number;

    chead = current;

    if( chead == NULL ) {
      continue;
    }


    if( chead->state & COMPARA_TARGET_UNIQUE )
      t_unique++;
    if( chead->state & COMPARA_QUERY_UNIQUE )
      q_unique++;

    if( (chead->state & COMPARA_TARGET_UNIQUE) &&
	(chead->state & COMPARA_QUERY_UNIQUE) ) 
      joint++;




    if( chead->spline != NULL ) {
      /* currently only put in if spline is unique */
      if( ((chead->state & COMPARA_TARGET_UNIQUE) &&
	  (chead->spline->state & COMPARA_QUERY_UNIQUE)) ||
	  ((chead->state & COMPARA_QUERY_UNIQUE) &&
	   (chead->spline->state & COMPARA_TARGET_UNIQUE))
	  )
	
	rev_joint++;
    }

    if( chead->state & COMPARA_TARGET_MULTIPLE )
      t_multiple++;
    if( chead->state & COMPARA_QUERY_MULTIPLE )
      q_multiple++;
  }

  fprintf(ofp,"Target %d unique positions, %d (rev: %d) shared (%.2f %%)\n",t_unique,joint,rev_joint,((joint+rev_joint)*100.0)/(t_unique*1.0));
  fprintf(ofp,"Query  %d unique positions, %d (rev: %d) shared (%.2f %%)\n",q_unique,joint,rev_joint,((joint+rev_joint)*100.0)/(q_unique*1.0));


}


# line 387 "comparapath.dy"
long int insert_revcom_Splines_in_set(ComparaIndex * ci,ComparaLinkStartSet * clss,FILE * logfp)
{
  int i;
  long int out = 0;

  for(i=0;i<clss->len;i++) {
    if( clss->cls[i]->start != NULL ) {
      out += insert_revcom_Splines(ci,clss->cls[i],logfp);
    }
  }

  return out;
}


# line 402 "comparapath.dy"
long int insert_revcom_Splines(ComparaIndex * ci,ComparaLinkStart * cls,FILE * logfp)
{
  long int total = 0;
  kmer_t rev;
  ComparaHead *chead;
  ComparaHead * revhead;
  ComparaHead * current;
  long int count = 0;


  for(current = cls->start ; current != NULL ; current = current->next_query) {
    chead = current;

    if( chead == NULL ) {
      continue;
    }

    if( count % 1000000 == 0 ) {
      fprintf(logfp,"Splines for %s, considered %ld positions\n",cls->seq->name,count);
      fflush(logfp);
    }
    count++;

    rev = reverse_complement_dna_number(current->number,ci->kii->kmer_size);

    revhead = (ComparaHead*) (*ci->kii->retrieve_by_kmer)(ci->kii->handle,rev);

    if( revhead == NULL ) {
      continue;
    }

    /*
    reverse_map_dna_number(current->number,ci->kii->kmer_size,forward);
    reverse_map_dna_number(rev,ci->kii->kmer_size,reverse);

    forward[ci->kii->kmer_size] = reverse[ci->kii->kmer_size] = '\0';

    fprintf(stderr,"Inserting spline at %ld with %s vs %s\n",current->position,forward,reverse);
    */

    chead->spline = revhead;
    revhead->spline = chead;
    total++;
    
  }

  return total;

}


# line 453 "comparapath.dy"
ComparaLinkStartSet * add_Sequence_stream_ComparaIndex(ComparaIndex * ci,FILE * ifp,boolean is_target,int lognumber,int test_rev,FILE * logfp,char * tag)
{
  Sequence * input;
  ComparaLinkStartSet * clss;
  ComparaLinkStart * cls;
  char buffer[512];

  /**
101001010100101010010
1001101
  */


  int skip_query[]  = {1,0,1,0,0,1,0};
  int skip_target[] = {1,0,0,1,1,0,1};

  clss = ComparaLinkStartSet_alloc_std();
  
  while( (input = read_large_dna_Sequence(ifp,lognumber,logfp)) != NULL ) {
    if( tag != NULL ) {
      sprintf(buffer,"%s_%s",tag,input->name);
      free(input->name);
      input->name = stringalloc(buffer);
    }

    if( is_target ) {
      cls = add_Sequence_ComparaIndex(ci,input,is_target,lognumber,0,7,skip_target,test_rev,logfp);
    } else {
      cls = add_Sequence_ComparaIndex(ci,input,is_target,lognumber,0,7,skip_query,test_rev,logfp);
    }
    add_ComparaLinkStartSet(clss,cls);
  }

  return clss;
}


# line 490 "comparapath.dy"
ComparaLinkStart * add_Sequence_ComparaIndex(ComparaIndex * ci,Sequence * seq,boolean is_target,int lognumber,long truncate,int skipsize,int * skipflag,int test_rev,FILE * logfp)
{
  int i;
  ComparaLinkStart * cls;
  long int pos_start;
  kmer_t number;
  ComparaHead * prev  = NULL;
  ComparaHead *chead;
  int j;
  int s;
  char base;

  kmer_t rev_number;
  ComparaHead * revhead;

  char * base_numbers;
  long int base_nos[100];
  int k;
  
  kmer_t skipped = 0;
  kmer_t repeated = 0;
  kmer_t seen_repeats = 0;

  assert(ci != NULL);
  assert(seq != NULL);

  assert(ci != NULL);
  assert(ci->kii != NULL);

  s = 1;
  for(i=0;i<ci->kii->kmer_size;i++) {
    base_nos[i] = s;
    s = s*4;
  }

  cls= new_ComparaLinkStart(seq);
  cls->start = NULL;

  pos_start = add_Sequence_SinglePosSpace(ci->sps,seq->len,(void*)seq);

  fprintf(logfp,"Loading position starting at %ld\n",pos_start);

  base_numbers = map_to_basepair_numbers(seq->seq,seq->len);

  for(i=0;i<seq->len - ci->kii->kmer_size;) {
    for(s=0;s<skipsize && i < seq->len - ci->kii->kmer_size;s++,i++) {

      if( lognumber != 0 && i % lognumber == 0 ) {
	fprintf(logfp,"Loading position %d (%ld skipped, %ld repeated, %ld repeats) from %s\n",i,skipped,repeated,seen_repeats,seq->name);
	fflush(logfp);
      }
  
      if( skipsize > 1 ) {
	if( skipflag[s] == 0 ) {
	  continue;
	}
      }
    
      if( truncate != 0 && i > truncate ) {
	break;
      }
      
      
      for(j=0;j<ci->kii->kmer_size;j++) {
	/*      if( seq->seq[i+j] == 'N' || islower(seq->seq[i+j]) ) {*/
	base = base_from_char(toupper(seq->seq[i+j]));
	if( base >= 4  ) {
	  break;
	}
      }
      if( j < ci->kii->kmer_size ) {
	skipped++;
	continue;
      }
      
      
      /*      number = forward_dna_number_from_string(seq->seq+i,ci->kii->kmer_size); */
      number = 0;
      for(k=0;k<ci->kii->kmer_size;k++) {
	number = number + (base_numbers[i+k]*base_nos[k]);
      }
      
      if( number < 0 ) {
	fprintf(logfp,"Yikes. Got bad number %ld at position %d, %.*s\n",number,i,ci->kii->kmer_size,seq->seq+i);
      }
      
      assert(number>= 0);
      
      chead = (ComparaHead*) (*ci->kii->retrieve_by_kmer)(ci->kii->handle,number);
      
      if( chead == NULL && is_target == 1 && test_rev == 1 ) {
	/* have to test reverse strand, otherwise pointless to put in */
	rev_number = reverse_complement_dna_number(number,ci->kii->kmer_size);
	revhead = (ComparaHead*) (*ci->kii->retrieve_by_kmer)(ci->kii->handle,rev_number);
	if( revhead == NULL ) {
	  continue;
	}
      }
      
      if( chead != NULL ) {
	if( is_target == 0 ) {
	  seen_repeats++;
	  chead->state |= COMPARA_IS_REPEATED;
	  chead->state = chead->state | COMPARA_QUERY_MULTIPLE;
	  continue;
	}

	if( (chead->state & COMPARA_IS_REPEATED)  ) {
	  seen_repeats++;
	  continue;
	}
	
	if( chead->size >= 2 ) {
	  chead->state |= COMPARA_IS_REPEATED;
	  repeated++;
	  continue;
	}
      }
      
      if( chead == NULL ) {
	chead = new_ComparaHead(ci->blockalloc);
	if( chead == NULL ) {
	  warn("Unable to make new ComparaHead for position %lld\n",i);
	  assert(chead);
	}
	chead->number = number;
	(*ci->kii->insert_by_kmer)(ci->kii->handle,number,chead);
	
	if( is_target ) {
	  chead->state = COMPARA_TARGET_UNIQUE;
	} else {
	  chead->state = COMPARA_QUERY_UNIQUE;
	}
      } else {
	/* must change state... */
	if( is_target ) {
	  if( (chead->state & COMPARA_TARGET_UNIQUE) ) {
	    chead->state = chead->state & (~COMPARA_TARGET_UNIQUE);
	    chead->state = chead->state | COMPARA_TARGET_MULTIPLE;
	  } else if( !(chead->state & COMPARA_TARGET_MULTIPLE) ) {
	    chead->state = chead->state | COMPARA_TARGET_UNIQUE;
	  }
	} else {
	  if( (chead->state & COMPARA_QUERY_UNIQUE) ) {
	    chead->state = chead->state & (~COMPARA_QUERY_UNIQUE);
	    chead->state = chead->state | COMPARA_QUERY_MULTIPLE;
	  } else if ( !(chead->state & COMPARA_QUERY_MULTIPLE) ) {
	    chead->state = chead->state | COMPARA_QUERY_UNIQUE;
	  }
	}
      }
      
      new_position_in_ComparaHead(chead,pos_start+i);
      
      if( cls->start == NULL ) {
	if( chead == NULL ) {
	  warn("Very weird; NULL chead at start assignment");
	} else {
	  cls->start = chead;
	}
      }
      
      if( is_target == 0 && prev != NULL ) {
	prev->next_query = chead;
	chead->next_query = NULL;
      }
      prev = chead;
    }
  }

  free(base_numbers);
  
  return cls;
}



# line 667 "comparapath.dy"
ComparaIndex * new_ComparaIndex(KmerIndexInterface * kii)
{
  ComparaIndex * out;
  
  assert(kii);

  out = malloc(sizeof(ComparaIndex));
  out->kii = kii;

  
  out->linkstart = calloc(COMPARAINDEX_LINK_START,sizeof(ComparaLinkStart*));
  out->current_link = 0;
  out->link_len = COMPARAINDEX_LINK_START;
  out->sps = new_SinglePosSpace(1,50000);
  out->blockalloc = new_ComparaHeadBlockAllocator(COMPARAHEAD_BA_BLOCK_LENGTH,COMPARAHEAD_BA_UNIT_LENGTH);
  return out;
}

# line 685 "comparapath.dy"
ComparaLinkStart * new_ComparaLinkStart(Sequence * seq)
{
  ComparaLinkStart * out;

  out = malloc(sizeof(ComparaLinkStart));
  
  out->start = NULL;
  out->seq = seq;

  return out;
}

# line 697 "comparapath.dy"
void add_ComparaLinkStart_to_ComparaIndex(ComparaIndex * ci,ComparaLinkStart * cls)
{
  assert(ci);
  assert(cls);

  if( ci->current_link >= ci->link_len ) {
    if( ci->link_len > COMPARAINDEX_LINK_LINEAR ) {
      ci->linkstart = realloc(ci->linkstart,sizeof(ComparaLinkStart*)*(ci->link_len + COMPARAINDEX_LINK_LINEAR));
      ci->link_len = ci->link_len + COMPARAINDEX_LINK_LINEAR;
    } else {
      ci->linkstart = realloc(ci->linkstart,sizeof(ComparaLinkStart*)*(ci->link_len*2));
      ci->link_len = ci->link_len*2;
    }
  }


  ci->linkstart[ci->current_link++] = cls;

}


# line 718 "comparapath.dy"
ComparaHead * new_ComparaHead(ComparaHeadBlockAllocator * ba)
{
  ComparaHead * out;

  out = new_ComparaHead_from_ComparaHeadBlockAllocator(ba);
  assert(out != NULL);

  out->size = 0;
  out->position[0] = -1;
  out->position[1] = -1;

  out->state    = COMPARA_NOTHING;
  out->spline  = NULL;
  return out;
}

# line 734 "comparapath.dy"
void new_position_in_ComparaHead(ComparaHead * h,long int position)
{
  assert(h != NULL);
  if( h->size == 2 ) {
    h->state |= COMPARA_IS_REPEATED;
    return;
  }

  h->position[(int)h->size] = position;
  h->size++;
  
  return;
}

# line 748 "comparapath.dy"
ComparaHeadBlockAllocator * new_ComparaHeadBlockAllocator(int block_length,int unit_length)
{
  ComparaHeadBlockAllocator * out;

  out = (ComparaHeadBlockAllocator*) malloc(sizeof(ComparaHeadBlockAllocator));
  out->current_block = 0;
  out->current_unit = 0;
  out->unit_length = unit_length;
  out->block_length = block_length;

  out->block = calloc(block_length,sizeof(ComparaHead*));
  out->block[0] = calloc(unit_length,sizeof(ComparaHead));

  return out;

}

# line 765 "comparapath.dy"
ComparaHead * new_ComparaHead_from_ComparaHeadBlockAllocator(ComparaHeadBlockAllocator * chba)
{
  ComparaHead * ret;
  assert(chba != NULL);
  if( chba->current_unit >= chba->unit_length ) {
    if( chba->current_block >= chba->block_length ) {
      chba->block = realloc(chba->block,sizeof(ComparaHead*)*(2*chba->block_length));
      chba->block_length *= 2;
    }
    chba->current_block++;
    chba->block[chba->current_block] = calloc(chba->unit_length,sizeof(ComparaHead));
    chba->current_unit = 0;
  }
  
  
  ret =  chba->block[chba->current_block]+(chba->current_unit);
  chba->current_unit++;

  return ret;

}

# line 724 "comparapath.c"
/* Function:  swap_ComparaLinkStartSet(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_ComparaLinkStartSet
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [ComparaLinkStart **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_ComparaLinkStartSet(ComparaLinkStart ** list,int i,int j)  
{
    ComparaLinkStart * temp; 
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_ComparaLinkStartSet(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_ComparaLinkStartSet which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [ComparaLinkStart **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_ComparaLinkStartSet(ComparaLinkStart ** list,int left,int right,int (*comp)(ComparaLinkStart * ,ComparaLinkStart * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_ComparaLinkStartSet(list,left,(left+right)/2);  
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_ComparaLinkStartSet (list,++last,i);    
      }  
    swap_ComparaLinkStartSet (list,left,last);   
    qsort_ComparaLinkStartSet(list,left,last-1,comp);    
    qsort_ComparaLinkStartSet(list,last+1,right,comp);   
}    


/* Function:  sort_ComparaLinkStartSet(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_ComparaLinkStartSet
 *
 *
 * Arg:         obj [UNKN ] Object containing list [ComparaLinkStartSet *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_ComparaLinkStartSet(ComparaLinkStartSet * obj,int (*comp)(ComparaLinkStart *, ComparaLinkStart *)) 
{
    qsort_ComparaLinkStartSet(obj->cls,0,obj->len-1,comp);   
    return;  
}    


/* Function:  expand_ComparaLinkStartSet(obj,len)
 *
 * Descrip:    Really an internal function for add_ComparaLinkStartSet
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [ComparaLinkStartSet *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_ComparaLinkStartSet(ComparaLinkStartSet * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_ComparaLinkStartSet called with no need");    
      return TRUE;   
      }  


    if( (obj->cls = (ComparaLinkStart ** ) ckrealloc (obj->cls,sizeof(ComparaLinkStart *)*len)) == NULL)     {  
      warn("ckrealloc failed for expand_ComparaLinkStartSet, returning FALSE");  
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_ComparaLinkStartSet(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [ComparaLinkStartSet *]
 * Arg:        add [OWNER] Object to add to the list [ComparaLinkStart *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_ComparaLinkStartSet(ComparaLinkStartSet * obj,ComparaLinkStart * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_ComparaLinkStartSet(obj,obj->len + ComparaLinkStartSetLISTLENGTH) == FALSE) 
        return FALSE;    
      }  


    obj->cls[obj->len++]=add;    
    return TRUE; 
}    


/* Function:  flush_ComparaLinkStartSet(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [ComparaLinkStartSet *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_ComparaLinkStartSet(ComparaLinkStartSet * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->cls[i] != NULL)   {  
        free_ComparaLinkStart(obj->cls[i]);  
        obj->cls[i] = NULL;  
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  ComparaLinkStartSet_alloc_std(void)
 *
 * Descrip:    Equivalent to ComparaLinkStartSet_alloc_len(ComparaLinkStartSetLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ComparaLinkStartSet *]
 *
 */
ComparaLinkStartSet * ComparaLinkStartSet_alloc_std(void) 
{
    return ComparaLinkStartSet_alloc_len(ComparaLinkStartSetLISTLENGTH); 
}    


/* Function:  ComparaLinkStartSet_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [ComparaLinkStartSet *]
 *
 */
ComparaLinkStartSet * ComparaLinkStartSet_alloc_len(int len) 
{
    ComparaLinkStartSet * out;  /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = ComparaLinkStartSet_alloc()) == NULL)  
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->cls = (ComparaLinkStart ** ) ckcalloc (len,sizeof(ComparaLinkStart *))) == NULL)    {  
      warn("Warning, ckcalloc failed in ComparaLinkStartSet_alloc_len"); 
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_ComparaLinkStartSet(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [ComparaLinkStartSet *]
 *
 * Return [UNKN ]  Undocumented return value [ComparaLinkStartSet *]
 *
 */
ComparaLinkStartSet * hard_link_ComparaLinkStartSet(ComparaLinkStartSet * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a ComparaLinkStartSet object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  ComparaLinkStartSet_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ComparaLinkStartSet *]
 *
 */
ComparaLinkStartSet * ComparaLinkStartSet_alloc(void) 
{
    ComparaLinkStartSet * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(ComparaLinkStartSet *) ckalloc (sizeof(ComparaLinkStartSet))) == NULL)  {  
      warn("ComparaLinkStartSet_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->cls = NULL; 
    out->len = out->maxlen = 0;  


    return out;  
}    


/* Function:  free_ComparaLinkStartSet(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [ComparaLinkStartSet *]
 *
 * Return [UNKN ]  Undocumented return value [ComparaLinkStartSet *]
 *
 */
ComparaLinkStartSet * free_ComparaLinkStartSet(ComparaLinkStartSet * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a ComparaLinkStartSet obj. Should be trappable");   
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
    if( obj->cls != NULL)    {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->cls[i] != NULL) 
          free_ComparaLinkStart(obj->cls[i]);    
        }  
      ckfree(obj->cls);  
      }  


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_SetofHSPset(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_SetofHSPset
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [HSPset **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_SetofHSPset(HSPset ** list,int i,int j)  
{
    HSPset * temp;   
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_SetofHSPset(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_SetofHSPset which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [HSPset **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_SetofHSPset(HSPset ** list,int left,int right,int (*comp)(HSPset * ,HSPset * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_SetofHSPset(list,left,(left+right)/2);  
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_SetofHSPset (list,++last,i);    
      }  
    swap_SetofHSPset (list,left,last);   
    qsort_SetofHSPset(list,left,last-1,comp);    
    qsort_SetofHSPset(list,last+1,right,comp);   
}    


/* Function:  sort_SetofHSPset(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_SetofHSPset
 *
 *
 * Arg:         obj [UNKN ] Object containing list [SetofHSPset *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_SetofHSPset(SetofHSPset * obj,int (*comp)(HSPset *, HSPset *)) 
{
    qsort_SetofHSPset(obj->hspset,0,obj->len-1,comp);    
    return;  
}    


/* Function:  expand_SetofHSPset(obj,len)
 *
 * Descrip:    Really an internal function for add_SetofHSPset
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [SetofHSPset *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_SetofHSPset(SetofHSPset * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_SetofHSPset called with no need");    
      return TRUE;   
      }  


    if( (obj->hspset = (HSPset ** ) ckrealloc (obj->hspset,sizeof(HSPset *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_SetofHSPset, returning FALSE");  
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_SetofHSPset(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [SetofHSPset *]
 * Arg:        add [OWNER] Object to add to the list [HSPset *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_SetofHSPset(SetofHSPset * obj,HSPset * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_SetofHSPset(obj,obj->len + SetofHSPsetLISTLENGTH) == FALSE) 
        return FALSE;    
      }  


    obj->hspset[obj->len++]=add; 
    return TRUE; 
}    


/* Function:  flush_SetofHSPset(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [SetofHSPset *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_SetofHSPset(SetofHSPset * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->hspset[i] != NULL)    {  
        free_HSPset(obj->hspset[i]); 
        obj->hspset[i] = NULL;   
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  SetofHSPset_alloc_std(void)
 *
 * Descrip:    Equivalent to SetofHSPset_alloc_len(SetofHSPsetLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SetofHSPset *]
 *
 */
SetofHSPset * SetofHSPset_alloc_std(void) 
{
    return SetofHSPset_alloc_len(SetofHSPsetLISTLENGTH); 
}    


/* Function:  SetofHSPset_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [SetofHSPset *]
 *
 */
SetofHSPset * SetofHSPset_alloc_len(int len) 
{
    SetofHSPset * out;  /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = SetofHSPset_alloc()) == NULL)  
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->hspset = (HSPset ** ) ckcalloc (len,sizeof(HSPset *))) == NULL) {  
      warn("Warning, ckcalloc failed in SetofHSPset_alloc_len"); 
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_SetofHSPset(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [SetofHSPset *]
 *
 * Return [UNKN ]  Undocumented return value [SetofHSPset *]
 *
 */
SetofHSPset * hard_link_SetofHSPset(SetofHSPset * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a SetofHSPset object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  SetofHSPset_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SetofHSPset *]
 *
 */
SetofHSPset * SetofHSPset_alloc(void) 
{
    SetofHSPset * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(SetofHSPset *) ckalloc (sizeof(SetofHSPset))) == NULL)  {  
      warn("SetofHSPset_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->hspset = NULL;  
    out->len = out->maxlen = 0;  


    return out;  
}    


/* Function:  free_SetofHSPset(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [SetofHSPset *]
 *
 * Return [UNKN ]  Undocumented return value [SetofHSPset *]
 *
 */
SetofHSPset * free_SetofHSPset(SetofHSPset * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a SetofHSPset obj. Should be trappable");   
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
    if( obj->hspset != NULL) {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->hspset[i] != NULL)  
          free_HSPset(obj->hspset[i]);   
        }  
      ckfree(obj->hspset);   
      }  


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif
