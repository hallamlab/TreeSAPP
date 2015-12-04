#ifndef CHHASH_H
#define CHHASH_H

#define CHHASH_BITS    24
#define CHHASH_MASK    0x0000000000ffffffL
#define CHHASH_ENTRIES 16777216

typedef struct ComparaHead ComparaHead;

typedef struct chbucket {
  int num_entries;
  int buffer_size;
  struct ComparaHead **entry;
  long long *index;
} CHBucket;

typedef struct chhash {
  CHBucket *bucket;
} CHHash;

void chhash_init_bucket(CHBucket*);
void chhash_fini_bucket(CHBucket*);

void chhash_init_hash(CHHash*);
void chhash_fini_hash(CHHash*);

void chhash_bucket_put(CHBucket*, ComparaHead*, long long pos);
ComparaHead *chhash_bucket_get(CHBucket*, long long pos);

void chhash_put(CHHash*, ComparaHead*, long long pos);
ComparaHead *chhash_get(CHHash*, long long pos);

void chhash_stats(CHHash*);

#endif /* !CHHASH_H */
