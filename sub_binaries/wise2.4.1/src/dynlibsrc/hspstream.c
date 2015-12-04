#ifdef _cplusplus
extern "C" {
#endif
#include "hspstream.h"

/* Function:  untyped_LinearHSPmanager_to_Stream(lm,ws)
 *
 * Descrip:    untyped linear hsp manager to stream
 *
 *
 * Arg:        lm [UNKN ] Undocumented argument [void *]
 * Arg:        ws [UNKN ] Undocumented argument [Wise2WriteStreamInterface *]
 *
 */
# line 16 "hspstream.dy"
void untyped_LinearHSPmanager_to_Stream(void * lm,Wise2WriteStreamInterface * ws)
{
  return typed_LinearHSPmanager_to_Stream((void*)lm,ws);
}

/* Function:  typed_LinearHSPmanager_to_Stream(lm,ws)
 *
 * Descrip:    typed linear hsp manager to stream
 *
 *             writes out sequneces first, with query as first
 *             sequence and then each target in turn, and then all the
 *             hsps in turn, only writing down the target name
 *
 *
 * Arg:        lm [UNKN ] Undocumented argument [LinearHSPmanager *]
 * Arg:        ws [UNKN ] Undocumented argument [Wise2WriteStreamInterface *]
 *
 */
# line 28 "hspstream.dy"
void typed_LinearHSPmanager_to_Stream(LinearHSPmanager * lm,Wise2WriteStreamInterface * ws)
{
  int i,j;
  SequenceSet * out;
  char buffer[512];

  /* assumme that each HSPset has only one target sequence */

  out = SequenceSet_alloc_len(lm->len+1);  

  /* query goes in first */
  add_SequenceSet(out,hard_link_Sequence(lm->query));
  for(i=0;i<lm->len;i++) {
    if( lm->set[i]->len > 0 ) {
      if( lm->set[i]->hsp[0]->target == NULL) {
	warn("Unable to transport HSP at position %d, as no linked target sequence",i);
      } else {
	add_SequenceSet(out,hard_link_Sequence(lm->set[i]->hsp[0]->target));
      }
    }
  }

  typed_write_SequenceSet_to_Stream(out,ws);
  free_SequenceSet(out);

  /* now write the hsps */

  for(i=0;i<lm->len;i++) {
    if( lm->set[i]->len > 0 && lm->set[i]->hsp[0]->target != NULL ) {
      for(j=0;j<lm->set[i]->len;j++) {
	sprintf(buffer,"%s %d %d %d %d %d\n",
	      lm->set[i]->hsp[j]->target->name,
		lm->set[i]->hsp[j]->query_start,
		lm->set[i]->hsp[j]->target_start,
		lm->set[i]->hsp[j]->length,
		lm->set[i]->hsp[j]->score,
		lm->set[i]->hsp[j]->target_reverse);
	/*	fprintf(stderr,"%s %d %d %d %d %d\n",
		lm->set[i]->hsp[j]->target->name,
		lm->set[i]->hsp[j]->query_start,
		lm->set[i]->hsp[j]->target_start,
		lm->set[i]->hsp[j]->length,
		lm->set[i]->hsp[j]->score,
		lm->set[i]->hsp[j]->target_reverse);
	*/

	WISE2_WRITE_STRING(buffer,ws);
      }
    }
  }

  WISE2_WRITE_STRING("//\n",ws);
  
}


/* Function:  untyped_LinearHSPmanager_from_Stream(rs)
 *
 * Descrip:    untyped read
 *
 *
 * Arg:        rs [UNKN ] Undocumented argument [Wise2ReadStreamInterface *]
 *
 * Return [UNKN ]  Undocumented return value [void *]
 *
 */
# line 87 "hspstream.dy"
void * untyped_LinearHSPmanager_from_Stream(Wise2ReadStreamInterface * rs)
{
  return (void*) typed_LinearHSPmanager_from_Stream(rs);
}


/* Function:  typed_LinearHSPmanager_from_Stream(rs)
 *
 * Descrip:    typed linear hsp manager from stream
 *
 *             assummes sequences coming in first, with query as first sequence
 *             and then hsps one after the other
 *
 *
 * Arg:        rs [UNKN ] Undocumented argument [Wise2ReadStreamInterface *]
 *
 * Return [UNKN ]  Undocumented return value [LinearHSPmanager *]
 *
 */
# line 99 "hspstream.dy"
LinearHSPmanager * typed_LinearHSPmanager_from_Stream(Wise2ReadStreamInterface * rs)
{

  SequenceSet * set;
  LinearHSPmanager * lm;
  HSPset * curr;
  HSP * h;

  int seq_pos;
  int dummy;
  int i;

  char buffer[MAXLINE];

  set = typed_SequenceSet_from_Stream(rs);

  /*
  for(i=0;i<set->len;i++) {
    fprintf(stderr,"Got %d, %s with %s\n",i,set->set[i]->name,set->set[i]->desc);
  }
  */


  lm = LinearHSPmanager_alloc_std();
  lm->query = hard_link_Sequence(set->set[0]);

  if( set->len == 1 ) {
    WISE2_READ_BUFFER(buffer,MAXLINE,rs);
    /* should check it is // */

    return lm;
  }



  seq_pos = 1;
  curr = HSPset_alloc_std();
  add_LinearHSPmanager(lm,curr);

  while( WISE2_READ_BUFFER(buffer,MAXLINE,rs) != NULL ) {
    auto int k;
    auto int temp_len = strlen(set->set[seq_pos]->name);

    if( strncmp(buffer,"//",2) == 0 ) {
      break;
    }
    

    /* is space is because we can have partial matching
       at the start of the of an identifier, in particular with
       uniprot alternative splice -X names
    */

    for(k=0;k < MAXLINE;k++) {
      if( isspace(buffer[k]) ) {
	buffer[k] = '\0';
	break;
      }
    }

    /*fprintf(stderr,"Comparing %s to %s with %d (%d) %d %d\n",buffer,set->set[seq_pos]->name,strncmp(buffer,set->set[seq_pos]->name,temp_len),temp_len,isspace(buffer[temp_len]),(int)(buffer[temp_len]) );
     */

    if( strcmp(buffer,set->set[seq_pos]->name) != 0 ) {

      seq_pos++;
      if( seq_pos >= set->len ) {
	warn("Overrun set sequences\n");
	break;
      }
      curr = HSPset_alloc_std();
      add_LinearHSPmanager(lm,curr);

      if( strncmp(buffer,set->set[seq_pos]->name,temp_len) != 0 ) {
	warn("Unable to match hsp name of [%s] to sequence name of %s, skipping\n",buffer,set->set[seq_pos]->name);
	continue;
      }


    }

    /* parse this hsp */
    h = HSP_alloc();
    h->query = hard_link_Sequence(set->set[0]);
    h->target = hard_link_Sequence(set->set[seq_pos]);

    sscanf(buffer+strlen(set->set[seq_pos]->name)+1,"%d %d %d %d %d",&h->query_start,&h->target_start,&h->length,&h->score,&dummy);
    if( dummy == 0 ) {
      h->target_reverse = 0;
    } else {
      h->target_reverse = 1;
    }
    
    add_HSPset(curr,h);

  }



  free_SequenceSet(set);


  return lm;
}

# line 218 "hspstream.c"

#ifdef _cplusplus
}
#endif
