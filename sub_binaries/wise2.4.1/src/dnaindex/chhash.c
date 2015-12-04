#include <stdio.h>
#include <stdlib.h>
#include "chhash.h"

void chhash_init_bucket(
CHBucket *b
){
  b->num_entries = 0;
  b->buffer_size = 0;
  b->entry = NULL;
  b->index = NULL;
}

void chhash_fini_bucket(
CHBucket *b
){
  free(b->index);
  free(b->entry);
}

void chhash_init_hash(
CHHash *hash
){
  int i;

  hash->bucket = (CHBucket*)malloc(sizeof(CHBucket)*CHHASH_ENTRIES);
  if(!hash->bucket){
    fprintf(stderr, "Cannot allocate buckets\n");
    exit(1);
  }

  for(i = 0; i < CHHASH_ENTRIES; ++i){
    chhash_init_bucket(&(hash->bucket[i]));
  }
}

void chhash_fini_hash(
CHHash *hash
){
  int i;

  for(i = 0; i < CHHASH_ENTRIES; ++i){
    chhash_fini_bucket(&(hash->bucket[i]));
  }

  free(hash->bucket);
}


void chhash_bucket_put(
CHBucket *b, 
ComparaHead *ch, 
long long pos
){
  int i;

  for(i = 0; i < b->num_entries; ++i){
    if(b->index[i] == pos){
      b->entry[i] = ch;
      return;
    }
  }

  if(b->num_entries >= b->buffer_size){
    b->buffer_size += 1024;
    b->entry = realloc(b->entry, sizeof(struct ComparaHead*)*b->buffer_size);
    if(!b->entry){
      exit(1);
    }
    b->index = realloc(b->index, sizeof(long long)*b->buffer_size);
    if(!b->index){
      exit(1);
    }
  }

  b->entry[b->num_entries] = ch;
  b->index[b->num_entries] = pos;
  b->num_entries++;
}

ComparaHead *chhash_bucket_get(
CHBucket *b, 
long long pos
){
  int i;

  for(i = 0; i < b->num_entries; ++i){
    if(b->index[i] == pos){
      return b->entry[i];
    }
  }

  return NULL;
}


void chhash_put(
CHHash *h,
ComparaHead *ch, 
long long pos
){
  long long bno = pos & CHHASH_MASK;

  chhash_bucket_put(&h->bucket[bno], ch, pos);
}

ComparaHead *chhash_get(
CHHash *h,
long long pos
){
  long long bno = pos & CHHASH_MASK;

  return chhash_bucket_get(&h->bucket[bno], pos);
}

void chhash_stats(
CHHash *hash
){
  int  i;
  int empty = 0;
  CHBucket *b;

  printf("Hash chains longer than 1\n");
  for(i = 0; i < CHHASH_ENTRIES; ++i){
    b = &(hash->bucket[i]);
    if(b->num_entries == 0){
      empty++;
      continue;
    }
    if(b->num_entries > 1){
      printf("%12d: %d   ", i, b->num_entries);
    }
  }
  printf("\n");

  printf("Empty Chains: %d/%d\n", empty, CHHASH_ENTRIES);
}
