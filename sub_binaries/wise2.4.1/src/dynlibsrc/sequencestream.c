#ifdef _cplusplus
extern "C" {
#endif
#include "sequencestream.h"


/* Function:  untyped_SequenceSet_from_Stream(rs)
 *
 * Descrip:    untyped version of reading from Wise2ReadStreamInterface
 *
 *
 * Arg:        rs [UNKN ] Undocumented argument [Wise2ReadStreamInterface *]
 *
 * Return [UNKN ]  Undocumented return value [void *]
 *
 */
# line 17 "sequencestream.dy"
void * untyped_SequenceSet_from_Stream(Wise2ReadStreamInterface * rs)
{
  SequenceSet * set;

  set = typed_SequenceSet_from_Stream(rs);

  return (void*) set;
}

/* Function:  typed_SequenceSet_from_Stream(rs)
 *
 * Descrip:    Typed version of reading from Stream, making a SequenceSet
 *
 *
 * Arg:        rs [UNKN ] Undocumented argument [Wise2ReadStreamInterface *]
 *
 * Return [UNKN ]  Undocumented return value [SequenceSet *]
 *
 */
# line 29 "sequencestream.dy"
SequenceSet * typed_SequenceSet_from_Stream(Wise2ReadStreamInterface * rs)
{
  SequenceSet * out;
  char buffer[512];
  Sequence * curr;

  char * pos;
  char * desc;

  out = SequenceSet_alloc_std();
  curr = NULL;

  while( WISE2_READ_BUFFER(buffer,512,rs) != NULL ) {
    if( strncmp(buffer,"//",2) == 0 ) {

/*
      fprintf(stderr,"In reading sequence set, have finally"
	" a sequence with length %d and first sequence %s\n",
	out->len,out->set[0]->seq);
*/

      return out;
    }
    if( buffer[0] == '>' ) {
      /* make new sequence with this position */

      /* delimit name after '>' */
      for(pos=buffer+1;*pos;pos++) {
	if( isspace(*pos) ) {
	  *pos = '\0';
	  break;
	}
      }

      /* find next string if there, will be description */

      for(++pos;*pos && isspace(*pos);pos++) {
	;
      }

      if( *pos != '\0' && !isspace(*pos) ) {
	desc = pos;
	/* could be '\n' delimited or '\0' delimited */
	for(;*pos != '\0' && *pos != '\n';pos++) {
	  ;
	}
	*pos = '\0';
      } else {
	desc = NULL;
      }

      /* if there has been a previous sequence,
	 it is in curr. Terminate the string at curr->len */

      if( curr != NULL ) {
	curr->seq[curr->len] = '\0';
      }


      curr = Sequence_alloc();
      curr->name = stringalloc(buffer+1);
      curr->seq  = calloc(512,sizeof(char));
      curr->len  = 0;
      curr->maxlen = 512;

      if( desc != NULL ) {
	curr->desc = stringalloc(desc);
      }

      add_SequenceSet(out,curr);

	/* fprintf(stderr, "SEQUENCE: %s\n", curr->name); */
      continue;
    }
    /* this is a sequence line */
    for(pos = buffer;*pos;pos++) {
      if( isalpha(*pos) ) {
	curr->seq[curr->len++] = *pos;
      } 
      if( curr->len+1 >= curr->maxlen ) {
	void *tmpp;

	if( curr->maxlen >= 32768 ) {
	  curr->maxlen += 32768;
	} else {
	  curr->maxlen *= 2;
	}
	tmpp = realloc(curr->seq,sizeof(char)*curr->maxlen);
	if (tmpp == NULL) {
		fatal("Failed in realloc()");
	}
	curr->seq = tmpp;
      }
    }
  }


	/*
	fprintf(stderr,"In reading sequence set, "
		"have finally a sequence with length %d and "
		"first sequence %s\n",out->len,out->set[0]->seq);
	*/

  return out;
  
}

/* Function:  typed_write_SequenceSet_to_Stream(set,ws)
 *
 * Descrip:    typed version of writing a sequence set to a stream
 *
 *
 * Arg:        set [UNKN ] Undocumented argument [SequenceSet *]
 * Arg:         ws [UNKN ] Undocumented argument [Wise2WriteStreamInterface *]
 *
 */
# line 139 "sequencestream.dy"
void typed_write_SequenceSet_to_Stream(SequenceSet * set,Wise2WriteStreamInterface * ws)
{
  int i;
  
  for(i=0;i<set->len;i++) {
    typed_write_one_Sequence_to_Stream(set->set[i],ws);
  }

  WISE2_WRITE_STRING("//\n",ws);

}

/* Function:  untyped_write_SequenceSet_to_Stream(s,ws)
 *
 * Descrip:    untyped version for writing a sequence set to a stream
 *
 *
 * Arg:         s [UNKN ] Undocumented argument [void *]
 * Arg:        ws [UNKN ] Undocumented argument [Wise2WriteStreamInterface *]
 *
 */
# line 154 "sequencestream.dy"
void untyped_write_SequenceSet_to_Stream(void * s,Wise2WriteStreamInterface * ws)
{
  typed_write_SequenceSet_to_Stream((SequenceSet*)s,ws);
}

/* Function:  typed_write_one_Sequence_to_Stream(seq,ws)
 *
 * Descrip:    internal function for writing out one sequence
 *
 *
 * Arg:        seq [UNKN ] Undocumented argument [Sequence *]
 * Arg:         ws [UNKN ] Undocumented argument [Wise2WriteStreamInterface *]
 *
 */
# line 162 "sequencestream.dy"
void typed_write_one_Sequence_to_Stream(Sequence * seq,Wise2WriteStreamInterface * ws)
{
  char buffer[512];
  int i;

 
  if( seq->desc != NULL) {
    sprintf(buffer,">%s %-80s\n",seq->name,seq->desc);
  } else {
    sprintf(buffer,">%s\n",seq->name);
  }

  WISE2_WRITE_STRING(buffer,ws);

  for(i=0;i<seq->len;) {
    buffer[0] = '\0';
    if( i+80 < seq->len ) {
      strncpy(buffer,seq->seq+i,80);
      i = i+80;
      buffer[80] = '\0';
    } else {
      strcpy(buffer,seq->seq+i);
      i = seq->len;
    }
    strcat(buffer,"\n");
    WISE2_WRITE_STRING(buffer,ws);
  }

}


/* Function:  untyped_write_Sequence_to_Stream(seq,ws)
 *
 * Descrip:    untyped write one sequence stream
 *
 *
 * Arg:        seq [UNKN ] Undocumented argument [void *]
 * Arg:         ws [UNKN ] Undocumented argument [Wise2WriteStreamInterface *]
 *
 */
# line 196 "sequencestream.dy"
void untyped_write_Sequence_to_Stream(void * seq,Wise2WriteStreamInterface * ws)
{
  typed_write_one_Sequence_to_Stream((Sequence *)seq,ws);
  WISE2_WRITE_STRING("//\n",ws);
}

/* Function:  untyped_read_Sequence_from_Stream(rs)
 *
 * Descrip:    reading one sequence from stream
 *
 *
 * Arg:        rs [UNKN ] Undocumented argument [Wise2ReadStreamInterface *]
 *
 * Return [UNKN ]  Undocumented return value [void *]
 *
 */
# line 205 "sequencestream.dy"
void * untyped_read_Sequence_from_Stream(Wise2ReadStreamInterface * rs)
{
  SequenceSet * out;
  Sequence * seq;

  out = typed_SequenceSet_from_Stream(rs);

  if( out == NULL || out->len == 0 ) {
    warn("In reading sequences, no sequences found");
    return NULL;
  }

  seq = hard_link_Sequence(out->set[0]);
  free_SequenceSet(out);

  
  return (void*) seq;
}



# line 256 "sequencestream.c"

#ifdef _cplusplus
}
#endif
