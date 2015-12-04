/*  RAxML-HPC, a program for sequential and parallel estimation of phylogenetic trees 
 *  Copyright March 2006 by Alexandros Stamatakis
 *
 *  Partially derived from
 *  fastDNAml, a program for estimation of phylogenetic trees from sequences by Gary J. Olsen
 *  
 *  and 
 *
 *  Programs of the PHYLIP package by Joe Felsenstein.
 *
 *  This program is free software; you may redistribute it and/or modify its
 *  under the terms of the GNU General Public License as published by the Free
 *  Software Foundation; either version 2 of the License, or (at your option)
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 *  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 *  for more details.
 * 
 *
 *  For any other enquiries send an Email to Alexandros Stamatakis
 *  stamatak@ics.forth.gr
 *
 *  When publishing work that is based on the results from RAxML-VI-HPC please cite:
 *  
 *  Alexandros Stamatakis: "An Efficient Program for phylogenetic Inference Using Simulated Annealing". 
 *  Proceedings of IPDPS2005,  Denver, Colorado, April 2005.
 *  
 *  AND
 *
 *  Alexandros Stamatakis:"RAxML-VI-HPC: maximum likelihood-based phylogenetic analyses with thousands of taxa and mixed models". 
 *  Bioinformatics 2006; doi: 10.1093/bioinformatics/btl446
 */


#ifndef WIN32  
#include <sys/times.h>
#include <sys/types.h>
#include <sys/time.h>
#include <unistd.h>  
#endif

#include <limits.h>
#include <math.h>
#include <time.h> 
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdint.h>
#include "axml.h"

extern FILE *INFILE;
extern char run_id[128];
extern char workdir[1024];
extern char bootStrapFile[1024];
extern char tree_file[1024];
extern char infoFileName[1024];
extern char resultFileName[1024];
extern volatile double          *reductionBuffer;
extern volatile double          *reductionBufferTwo;
extern volatile double          *reductionBufferThree;
extern volatile int             NumberOfThreads;
extern char inverseMeaningBINARY[4];

extern double masterTime;

extern const int mask32[32];




      



static entry *initEntry(void)
{
  entry *e = (entry*)malloc(sizeof(entry));

  e->bitVector     = (unsigned int*)NULL;
  e->treeVector    = (unsigned int*)NULL;
  e->supportVector = (int*)NULL;
  e->bipNumber  = 0;
  e->bipNumber2 = 0;
  e->next       = (entry*)NULL;

  return e;
} 

hashtable *initHashTable(hashNumberType n, unsigned int numberOfTips)
{
  /*static const hashNumberType initTable[] = {53, 97, 193, 389, 769, 1543, 3079, 6151, 12289, 24593, 49157, 98317,
					  196613, 393241, 786433, 1572869, 3145739, 6291469, 12582917, 25165843,
					  50331653, 100663319, 201326611, 402653189, 805306457, 1610612741};*/

  static const  hashNumberType initTable[] = {64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384,
					      32768, 65536, 131072, 262144, 524288, 1048576, 2097152,
					      4194304, 8388608, 16777216, 33554432, 67108864, 134217728,
					      268435456, 536870912, 1073741824, 2147483648};
  /* powers of two */
  
  hashtable *h = (hashtable*)malloc(sizeof(hashtable));
  
  hashNumberType
    tableSize,
    i,
    primeTableLength = sizeof(initTable)/sizeof(initTable[0]),
    maxSize = (hashNumberType)-1;    
  long seed = 12345;

  assert(n <= maxSize);

  i = 0;

  while(initTable[i] < n && i < primeTableLength)
    i++;

  assert(i < primeTableLength);

  tableSize = initTable[i];

  printf("Hash table init with size %u\n", tableSize);

  h->table = (entry**)calloc(tableSize, sizeof(entry*));
  h->tableSize = tableSize;  
  h->entryCount = 0;
 
  h->randomNumbers = (hashNumberType*)malloc(sizeof(hashNumberType) * numberOfTips); 

  for(i = 0; i < numberOfTips; i++)
    {
      /* SOS caution here with unsigned 32 and 64 integers */
      h->randomNumbers[i] = (hashNumberType)(randum(&seed) * tableSize);
      assert(h->randomNumbers[i] < tableSize);
    }

  
  return h;
}

void freeHashTable(hashtable *h)
{
  hashNumberType
    i,
    entryCount = 0,
    maxChainLength = 0,
    collisionCount = 0;

  for(i = 0; i < h->tableSize; i++)
    {
      if(h->table[i] != NULL)
	{
	  entry *e = h->table[i];
	  entry *previous;
	  hashNumberType count = 0;

	  do
	    {
	      previous = e;
	      e = e->next;

	      if(previous->bitVector)
		free(previous->bitVector);
	      if(previous->treeVector)
		free(previous->treeVector);
	      if(previous->supportVector)
		free(previous->supportVector);
	      count++;
	      entryCount++;
	    }
	  while(e != NULL);

	  if(count > 1)
	    {
	      collisionCount += (count - 1);
	      if(count > maxChainLength)
		maxChainLength = count;
	    }
	}

    }
  assert(entryCount == h->entryCount);

  free(h->randomNumbers);
  free(h->table);

  printf("Hash Table Stats: \n max chain: %u \n collisions: %u \n entries: %u\n", maxChainLength, collisionCount, entryCount);

}




#define NN 312
#define MM 156
#define MATRIX_A 0xB5026F5AA96619E9LLU
#define UM 0xFFFFFFFF80000000LLU /* Most significant 33 bits */
#define LM 0x7FFFFFFFLLU /* Least significant 31 bits */

static uint64_t mt[NN]; 
static int mti=NN+1; 

/* initializes mt[NN] with a seed */
static void init_genrand64(uint64_t seed)
{
    mt[0] = seed;
    for (mti=1; mti<NN; mti++) 
        mt[mti] =  (6364136223846793005LLU * (mt[mti-1] ^ (mt[mti-1] >> 62)) + mti);
}

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */

static void init_by_array64(uint64_t init_key[],
			    uint64_t key_length)
{
    uint64_t i, j, k;
    init_genrand64(19650218LLU);
    i=1; j=0;
    k = (NN>key_length ? NN : key_length);
    for (; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 62)) * 3935559000370003845LLU))
          + init_key[j] + j; /* non linear */
        i++; j++;
        if (i>=NN) { mt[0] = mt[NN-1]; i=1; }
        if (j>=key_length) j=0;
    }
    for (k=NN-1; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 62)) * 2862933555777941757LLU))
          - i; /* non linear */
        i++;
        if (i>=NN) { mt[0] = mt[NN-1]; i=1; }
    }

    mt[0] = 1LLU << 63; /* MSB is 1; assuring non-zero initial array */ 
}

/* generates a random number on [0, 2^64-1]-interval */
static uint64_t genrand64_int64(void)
{
    int i;
    uint64_t x;
    static uint64_t mag01[2]={0LLU, MATRIX_A};

    if (mti >= NN) { /* generate NN words at one time */

        /* if init_genrand64() has not been called, */
        /* a default initial seed is used     */
        if (mti == NN+1) 
            init_genrand64(5489LLU); 

        for (i=0;i<NN-MM;i++) {
            x = (mt[i]&UM)|(mt[i+1]&LM);
            mt[i] = mt[i+MM] ^ (x>>1) ^ mag01[(int)(x&1LLU)];
        }
        for (;i<NN-1;i++) {
            x = (mt[i]&UM)|(mt[i+1]&LM);
            mt[i] = mt[i+(MM-NN)] ^ (x>>1) ^ mag01[(int)(x&1LLU)];
        }
        x = (mt[NN-1]&UM)|(mt[0]&LM);
        mt[NN-1] = mt[MM-1] ^ (x>>1) ^ mag01[(int)(x&1LLU)];

        mti = 0;
    }
  
    x = mt[mti++];

    x ^= (x >> 29) & 0x5555555555555555LLU;
    x ^= (x << 17) & 0x71D67FFFEDA60000LLU;
    x ^= (x << 37) & 0xFFF7EEE000000000LLU;
    x ^= (x >> 43);

    return x;
}

/* generates a random number on [0, 2^63-1]-interval */
/*long long genrand64_int63(void)
{
    return (long long)(genrand64_int64() >> 1);
    }*/




unsigned int **initBitVector(tree *tr, int *vectorLength)
{
  unsigned int **bitVectors = (unsigned int **)malloc(sizeof(unsigned int*) * 2 * tr->mxtips);
  int i;

  if(tr->mxtips % MASK_LENGTH == 0)
    *vectorLength = tr->mxtips / MASK_LENGTH;
  else
    *vectorLength = 1 + (tr->mxtips / MASK_LENGTH); 
  
  for(i = 1; i <= tr->mxtips; i++)
    {
      bitVectors[i] = (unsigned int *)calloc(*vectorLength, sizeof(unsigned int));
      bitVectors[i][(i - 1) / MASK_LENGTH] |= mask32[(i - 1) % MASK_LENGTH];
    }
  
  for(i = tr->mxtips + 1; i < 2 * tr->mxtips; i++)
    bitVectors[i] = (unsigned int *)malloc(sizeof(unsigned int) * *vectorLength);

  return bitVectors;
}

void freeBitVectors(unsigned int **v, int n)
{
  int i;

  for(i = 1; i < n; i++)
    free(v[i]);
}

static inline hashNumberType computeHash(unsigned int *bitVector, hashtable *h, int vectorLength)
{
  hashNumberType hashIndex = 0;
  hashNumberType *randomNumbers = h->randomNumbers;
  int i;

  for(i = 0; i < vectorLength; i++)         
    hashIndex += randomNumbers[i] * bitVector[i];    

  hashIndex = hashIndex % h->tableSize;

  assert(hashIndex < h->tableSize);

  return hashIndex;
}


#if !defined (get16bits)
#define get16bits(d) ((((uint32_t)(((const uint8_t *)(d))[1])) << 8)\
                       +(uint32_t)(((const uint8_t *)(d))[0]) )
#endif

static uint32_t computeHash6(const char * data, int len) 
{
uint32_t hash = len, tmp;
int rem;

    if (len <= 0 || data == NULL) return 0;

    rem = len & 3;
    len >>= 2;

    /* Main loop */
    for (;len > 0; len--) {
        hash  += get16bits (data);
        tmp    = (get16bits (data+2) << 11) ^ hash;
        hash   = (hash << 16) ^ tmp;
        data  += 2*sizeof (uint16_t);
        hash  += hash >> 11;
    }

    /* Handle end cases */
    switch (rem) {
        case 3: hash += get16bits (data);
                hash ^= hash << 16;
                hash ^= data[sizeof (uint16_t)] << 18;
                hash += hash >> 11;
                break;
        case 2: hash += get16bits (data);
                hash ^= hash << 11;
                hash += hash >> 17;
                break;
        case 1: hash += *data;
                hash ^= hash << 10;
                hash += hash >> 1;
    }

    /* Force "avalanching" of final 127 bits */
    hash ^= hash << 3;
    hash += hash >> 5;
    hash ^= hash << 4;
    hash += hash >> 17;
    hash ^= hash << 25;
    hash += hash >> 6;

    return hash;
}


/*ub4 one_at_a_time(char *key, ub4 len)
{
  ub4   hash, i;
  for (hash=0, i=0; i<len; ++i)
  {
    hash += key[i];
    hash += (hash << 10);
    hash ^= (hash >> 6);
  }
  hash += (hash << 3);
  hash ^= (hash >> 11);
  hash += (hash << 15);
  return (hash & mask);
  } */


static inline unsigned int computeHash4(unsigned int *key, unsigned int key_len)
{
  unsigned 
    int hash = 0,
    i,
    length = 4 * key_len;
  unsigned char *ckey = (unsigned char *)key;
 
  for (i = 0; i < length; ++i) 
    {
      hash += ckey[i];
      hash += (hash << 10);
      hash ^= (hash >> 6);
    }
  hash += (hash << 3);
  hash ^= (hash >> 11);
  hash += (hash << 15);
  
  return hash;
}

static inline unsigned int computeHash3(unsigned int *key, unsigned int key_len)
{
  unsigned int hash, i;
 
  for(i = 0, hash = 0; i < key_len; ++i) 
    {
      hash += key[i];
      hash += (hash << 10);
      hash ^= (hash >> 6);
    }
  
  hash += (hash << 3);
  hash ^= (hash >> 11);
  hash += (hash << 15);
  
  return hash;
}

static inline unsigned int computeHash2(const void *key, int len, unsigned int seed)
{
  const unsigned int m = 0x5bd1e995;
  const int r = 24;

  unsigned int h = seed ^ len;

  const unsigned char * data = (const unsigned char *)key;

  while(len >= 4)
    {
      unsigned int k;

      k  = data[0];
      k |= data[1] << 8;
      k |= data[2] << 16;
      k |= data[3] << 24;
      
      k *= m; 
      k ^= k >> r; 
      k *= m;
      
      h *= m;
      h ^= k;
      
      data += 4;
      len -= 4;
    }
	
  switch(len)
    {
    case 3: h ^= data[2] << 16;
    case 2: h ^= data[1] << 8;
    case 1: h ^= data[0];
      h *= m;
    }
  ;

  h ^= h >> 13;
  h *= m;
  h ^= h >> 15;
  
  return h;
} 








static void newviewBipartitions(unsigned int **bitVectors, nodeptr p, int numsp, int vectorLength)
{
  if(isTip(p->number, numsp))
    return;
  {
    nodeptr 
      q = p->next->back, 
      r = p->next->next->back;
    unsigned int       
      *vector = bitVectors[p->number],
      *left  = bitVectors[q->number],
      *right = bitVectors[r->number];
    int i;           

    while(!p->x)
      {	
	if(!p->x)
	  getxnode(p);
      }

    p->hash = q->hash ^ r->hash;

    if(isTip(q->number, numsp) && isTip(r->number, numsp))
      {       		
	for(i = 0; i < vectorLength; i++)
	  vector[i] = left[i] | right[i];	  	
      }
    else
      {	
	if(isTip(q->number, numsp) || isTip(r->number, numsp))
	  {
	    if(isTip(r->number, numsp))
	      {	
		nodeptr tmp = r;
		r = q;
		q = tmp;
	      }	   
	    	    
	    while(!r->x)
	      {
		if(!r->x)
		  newviewBipartitions(bitVectors, r, numsp, vectorLength);
	      }	   

	    for(i = 0; i < vectorLength; i++)
	      vector[i] = left[i] | right[i];	    	 
	  }
	else
	  {	    
	    while((!r->x) || (!q->x))
	      {
		if(!q->x)
		  newviewBipartitions(bitVectors, q, numsp, vectorLength);
		if(!r->x)
		  newviewBipartitions(bitVectors, r, numsp, vectorLength);
	      }	   	    	    	    	   

	    for(i = 0; i < vectorLength; i++)
	      vector[i] = left[i] | right[i];	    
	  }

      }     
  }     
}

static void insertHash(unsigned int *bitVector, hashtable *h, int vectorLength, int bipNumber, hashNumberType position)
{
  entry *e = initEntry();

  e->bipNumber = bipNumber; 
  e->bitVector = (unsigned int*)calloc(vectorLength, sizeof(unsigned int)); 
 
  memcpy(e->bitVector, bitVector, sizeof(unsigned int) * vectorLength);
  
  if(h->table[position] != NULL)
    {
      e->next = h->table[position];
      h->table[position] = e;           
    }
  else
    h->table[position] = e;

  h->entryCount =  h->entryCount + 1;
}




static int countHash(unsigned int *bitVector, hashtable *h, int vectorLength, hashNumberType position)
{
  if(h->table[position] == NULL)
    return -1;
  {
    entry *e = h->table[position];     

    do
      {	 
	int i;

	for(i = 0; i < vectorLength; i++)
	  if(bitVector[i] != e->bitVector[i])
	    break;

	if(i == vectorLength)
	  return (e->bipNumber);

	e = e->next;
      }
    while(e != (entry*)NULL); 
  
    return -1;
   
  }

}

static void insertHashAll(unsigned int *bitVector, hashtable *h, int vectorLength, int treeNumber,  hashNumberType position)
{    
  if(h->table[position] != NULL)
    {
      entry *e = h->table[position];     

      do
	{	 
	  int i;
	  
	  for(i = 0; i < vectorLength; i++)
	    if(bitVector[i] != e->bitVector[i])
	      break;
	  
	  if(i == vectorLength)
	    {
	      if(treeNumber == 0)
		e->bipNumber = 	e->bipNumber  + 1;
	      else
		e->bipNumber2 = e->bipNumber2 + 1;
	      return;
	    }
	  
	  e = e->next;	 
	}
      while(e != (entry*)NULL); 

      e = initEntry(); 
  
      e->bitVector  = (unsigned int*)calloc(vectorLength, sizeof(unsigned int)); 
      memcpy(e->bitVector, bitVector, sizeof(unsigned int) * vectorLength);

      if(treeNumber == 0)	
	e->bipNumber  = 1;       	
      else		 
	e->bipNumber2 = 1;
	
      e->next = h->table[position];
      h->table[position] = e;              
    }
  else
    {
      entry *e = initEntry(); 
  
      e->bitVector  = (unsigned int*)calloc(vectorLength, sizeof(unsigned int)); 
      memcpy(e->bitVector, bitVector, sizeof(unsigned int) * vectorLength);

      if(treeNumber == 0)	
	e->bipNumber  = 1;	  	
      else    
	e->bipNumber2 = 1;	

      h->table[position] = e;
    }

  h->entryCount =  h->entryCount + 1;
}




static void insertHashBootstop(unsigned int *bitVector, hashtable *h, int vectorLength, int treeNumber, int treeVectorLength, hashNumberType position)
{    
  if(h->table[position] != NULL)
    {
      entry *e = h->table[position];     

      do
	{	 
	  int i;
	  
	  for(i = 0; i < vectorLength; i++)
	    if(bitVector[i] != e->bitVector[i])
	      break;
	  
	  if(i == vectorLength)
	    {
	      e->treeVector[treeNumber / MASK_LENGTH] |= mask32[treeNumber % MASK_LENGTH];
	      return;
	    }
	  
	  e = e->next;
	}
      while(e != (entry*)NULL); 

      e = initEntry(); 
       
      e->bitVector  = (unsigned int*)calloc(vectorLength, sizeof(unsigned int));
      e->treeVector = (unsigned int*)calloc(treeVectorLength, sizeof(unsigned int));
      
      e->treeVector[treeNumber / MASK_LENGTH] |= mask32[treeNumber % MASK_LENGTH];
      memcpy(e->bitVector, bitVector, sizeof(unsigned int) * vectorLength);
     
      e->next = h->table[position];
      h->table[position] = e;          
    }
  else
    {
      entry *e = initEntry(); 
       
      e->bitVector  = (unsigned int*)calloc(vectorLength, sizeof(unsigned int)); 
      e->treeVector = (unsigned int*)calloc(treeVectorLength, sizeof(unsigned int));

      e->treeVector[treeNumber / MASK_LENGTH] |= mask32[treeNumber % MASK_LENGTH];
      memcpy(e->bitVector, bitVector, sizeof(unsigned int) * vectorLength);     

      h->table[position] = e;
    }

  h->entryCount =  h->entryCount + 1;
}

static void insertHashRF(unsigned int *bitVector, hashtable *h, int vectorLength, int treeNumber, int treeVectorLength, hashNumberType position, int support, 
			 boolean computeWRF)
{    
  if(h->table[position] != NULL)
    {
      entry *e = h->table[position];     

      do
	{	 
	  int i;
	  
	  for(i = 0; i < vectorLength; i++)
	    if(bitVector[i] != e->bitVector[i])
	      break;
	  
	  if(i == vectorLength)
	    {
	      e->treeVector[treeNumber / MASK_LENGTH] |= mask32[treeNumber % MASK_LENGTH];
	      if(computeWRF)
		e->supportVector[treeNumber] = support;
	      return;
	    }
	  
	  e = e->next;
	}
      while(e != (entry*)NULL); 

      e = initEntry(); 
       
      e->bitVector  = (unsigned int*)calloc(vectorLength, sizeof(unsigned int));
      e->treeVector = (unsigned int*)calloc(treeVectorLength, sizeof(unsigned int));
      if(computeWRF)
	e->supportVector = (int*)calloc(treeVectorLength * MASK_LENGTH, sizeof(int));

      e->treeVector[treeNumber / MASK_LENGTH] |= mask32[treeNumber % MASK_LENGTH];
      if(computeWRF)
	e->supportVector[treeNumber] = support;

      memcpy(e->bitVector, bitVector, sizeof(unsigned int) * vectorLength);
     
      e->next = h->table[position];
      h->table[position] = e;          
    }
  else
    {
      entry *e = initEntry(); 
       
      e->bitVector  = (unsigned int*)calloc(vectorLength, sizeof(unsigned int)); 
      e->treeVector = (unsigned int*)calloc(treeVectorLength, sizeof(unsigned int));
      if(computeWRF)
	e->supportVector = (int*)calloc(treeVectorLength * MASK_LENGTH, sizeof(int));

      e->treeVector[treeNumber / MASK_LENGTH] |= mask32[treeNumber % MASK_LENGTH];
      if(computeWRF)
	e->supportVector[treeNumber] = support;

      memcpy(e->bitVector, bitVector, sizeof(unsigned int) * vectorLength);     

      h->table[position] = e;
    }

  h->entryCount =  h->entryCount + 1;
}





static void bitVectorInitravSpecial(unsigned int **bitVectors, nodeptr p, int numsp, int vectorLength, hashtable *h, int treeNumber, int function, branchInfo *bInf, 
				    int *countBranches, int treeVectorLength, boolean traverseOnly, boolean computeWRF)
{
  if(isTip(p->number, numsp))
    return;
  else
    {
      nodeptr q = p->next;          

      do 
	{
	  bitVectorInitravSpecial(bitVectors, q->back, numsp, vectorLength, h, treeNumber, function, bInf, countBranches, treeVectorLength, traverseOnly, computeWRF);
	  q = q->next;
	}
      while(q != p);
      
      newviewBipartitions(bitVectors, p, numsp, vectorLength);
      assert(p->x);

      if(traverseOnly)
	{
	  if(!(isTip(p->back->number, numsp)))
	    *countBranches =  *countBranches + 1;
	  return;
	}

      if(!(isTip(p->back->number, numsp)))
	{
	  unsigned int *toInsert  = bitVectors[p->number];
	  hashNumberType position;

	  /* 
	     position = computeHash(toInsert, h, vectorLength);
	     position = p->hash % h->tableSize;
	     position = computeHash6((const char *)toInsert, vectorLength) % h->tableSize;
	     position = computeHash4(toInsert, (unsigned int)vectorLength) % h->tableSize;
	     position = computeHash3(toInsert, (unsigned int)vectorLength) % h->tableSize;
	  */

	  position = p->hash % h->tableSize;
	 
	  assert(!(toInsert[0] & 1));	 

	  switch(function)
	    {
	    case BIPARTITIONS_ALL:	      
	      insertHashAll(toInsert, h, vectorLength, treeNumber, position);
	      *countBranches =  *countBranches + 1;	
	      break;
	    case GET_BIPARTITIONS_BEST:	      
	      insertHash(toInsert, h, vectorLength, *countBranches, position);	     
	      p->bInf            = &bInf[*countBranches];
	      p->back->bInf      = &bInf[*countBranches];        
	      p->bInf->support   = 0;	  	 
	      p->bInf->oP = p->number;
	      p->bInf->oQ = p->back->number;	      
	      *countBranches =  *countBranches + 1;
	      break;
	    case DRAW_BIPARTITIONS_BEST:
	      {
		int found = countHash(toInsert, h, vectorLength, position);
		if(found >= 0)
		  bInf[found].support =  bInf[found].support + 1;
		*countBranches =  *countBranches + 1;
	      }
	      break;
	    case BIPARTITIONS_BOOTSTOP:	      
	      insertHashBootstop(toInsert, h, vectorLength, treeNumber, treeVectorLength, position);
	      *countBranches =  *countBranches + 1;
	      break;
	    case BIPARTITIONS_RF:
	      if(computeWRF)
		assert(p->support = p->back->support);
	      insertHashRF(toInsert, h, vectorLength, treeNumber, treeVectorLength, position, p->support, computeWRF);
	      *countBranches =  *countBranches + 1;
	      break;
	    default:
	      assert(0);
	    }	  	  
	}
      
    }
}

static void linkBipartitions(nodeptr p, tree *tr, branchInfo *bInf, int *counter, int numberOfTrees)
{
  if(isTip(p->number, tr->mxtips))    
    {
      assert(p->bInf == (branchInfo*) NULL && p->back->bInf == (branchInfo*) NULL);      
      return;
    }
  else
    {
      nodeptr q;          
      
      q = p->next;

      while(q != p)
	{
	  linkBipartitions(q->back, tr, bInf, counter, numberOfTrees);	
	  q = q->next;
	}
     
      if(!(isTip(p->back->number, tr->mxtips)))
	{
	  double support;

	  p->bInf       = &bInf[*counter];
	  p->back->bInf = &bInf[*counter]; 

	  support = ((double)(p->bInf->support)) / ((double) (numberOfTrees));
	  p->bInf->support = (int)(0.5 + support * 100.0);	 	       
	
	  assert(p->bInf->oP == p->number);
	  assert(p->bInf->oQ == p->back->number);
	  

	  *counter = *counter + 1;
	}


      return;
    }
}





void calcBipartitions(tree *tr, analdef *adef, char *bestTreeFileName, char *bootStrapFileName)
{
  int numberOfTrees = 0, i;
  int vLength;
  FILE *infoFile;
  branchInfo *bInf;
  int branchCounter = 0;
  int counter = 0;
  unsigned int **bitVectors = initBitVector(tr, &vLength);
  hashtable *h = initHashTable(tr->mxtips * 10, tr->mxtips);       
 
  INFILE = myfopen(bestTreeFileName, "r");
  treeReadTopologyOnly(INFILE, tr, FALSE, FALSE, FALSE); 
  assert(tr->ntips == tr->mxtips);
  fclose(INFILE);  

  bInf = (branchInfo*)malloc(sizeof(branchInfo) * (tr->mxtips - 3));

  bitVectorInitravSpecial(bitVectors, tr->nodep[1]->back, tr->mxtips, vLength, h, 0, GET_BIPARTITIONS_BEST, bInf, &branchCounter, 0, FALSE, FALSE);   
 
  assert((int)h->entryCount == (tr->mxtips - 3));  
  assert(branchCounter == (tr->mxtips - 3));
  
  INFILE = myfopen(bootStrapFileName, "r");       
  numberOfTrees = countTrees(INFILE); 
 
  if(!adef->allInOne)
    {
      infoFile = myfopen(infoFileName, "a");
      printf("\n\nFound %d trees in File %s\n\n", numberOfTrees, bootStrapFileName);
      fprintf(infoFile, "\n\nFound %d trees in File %s\n\n", numberOfTrees, bootStrapFileName);
      fclose(infoFile);
    }

  for(i = 0; i < numberOfTrees; i++)
    {                
      int bCount = 0;      
      treeReadTopologyOnly(INFILE, tr, FALSE, FALSE, FALSE);        
      assert(tr->ntips == tr->mxtips);
      bitVectorInitravSpecial(bitVectors, tr->nodep[1]->back, tr->mxtips, vLength, h, 0, DRAW_BIPARTITIONS_BEST, bInf, &bCount, 0, FALSE, FALSE);
      assert(bCount == tr->mxtips - 3);
      
    }

  fclose(INFILE);  
 
  INFILE = myfopen(bestTreeFileName, "r");
  treeReadTopologyOnly(INFILE, tr, TRUE, FALSE, FALSE);
  fclose(INFILE); 
   
  linkBipartitions(tr->nodep[1]->back, tr, bInf, &counter, numberOfTrees);

  assert(counter == branchCounter);

  printBipartitionResult(tr, adef, TRUE);    

  freeBitVectors(bitVectors, 2 * tr->mxtips);
  free(bitVectors);
  freeHashTable(h);
  free(h); 

  free(bInf);
}

/*************************************************************/

static double testFreq(double *vect1, double *vect2, int n, tree *tr);



void compareBips(tree *tr, char *bootStrapFileName)
{
  int 
    numberOfTreesAll = 0, 
    numberOfTreesStop = 0,
    i; 
  unsigned int k, entryCount;
  double *vect1, *vect2, p, avg1 = 0.0, avg2 = 0.0, scaleAll, scaleStop;
  int 
    bipAll = 0,
    bipStop = 0;
  char bipFileName[1024];
  FILE *outf; 
  int vLength;
  unsigned int **bitVectors = initBitVector(tr, &vLength);
  hashtable *h = initHashTable(tr->mxtips * 100, tr->mxtips);    
  unsigned long int c1 = 0;
  unsigned long int c2 = 0;  


 


  INFILE = myfopen(bootStrapFileName, "r"); 
  numberOfTreesAll = countTrees(INFILE);        

  printf("\n\nFound %d trees in File %s\n\n", numberOfTreesAll, bootStrapFileName);
              
  for(i = 0; i < numberOfTreesAll; i++)
    { 
      int bCounter = 0;
      
      treeReadTopologyOnly(INFILE, tr, FALSE, FALSE, FALSE);
      assert(tr->mxtips == tr->ntips); 
      bitVectorInitravSpecial(bitVectors, tr->nodep[1]->back, tr->mxtips, vLength, h, 0, BIPARTITIONS_ALL, (branchInfo*)NULL, &bCounter, 0, FALSE, FALSE);
      assert(bCounter == tr->mxtips - 3);      
    }
	  
  fclose(INFILE); 

  /* do BOOTSTOP ********************************************************************************************************/  

  INFILE = myfopen(tree_file, "r");              
  numberOfTreesStop = countTrees(INFILE);


  printf("\n\nFound %d trees in File %s\n\n", numberOfTreesStop, tree_file);
       
  for(i = 0; i < numberOfTreesStop; i++)
    {              
      int bCounter = 0;

      treeReadTopologyOnly(INFILE, tr, FALSE, FALSE, FALSE);
      assert(tr->mxtips == tr->ntips);
      bitVectorInitravSpecial(bitVectors, tr->nodep[1]->back, tr->mxtips, vLength, h, 1, BIPARTITIONS_ALL, (branchInfo*)NULL, &bCounter, 0, FALSE, FALSE);
      assert(bCounter == tr->mxtips - 3);     
    }
	  
  fclose(INFILE); 

  /***************************************************************************************************/
   
 
  
  vect1 = (double *)malloc(sizeof(double) * h->entryCount);
  vect2 = (double *)malloc(sizeof(double) * h->entryCount);

  strcpy(bipFileName,         workdir);  
  strcat(bipFileName,         "RAxML_bipartitionFrequencies.");
  strcat(bipFileName,         run_id);

  outf = myfopen(bipFileName, "w");


  scaleAll  = 1.0 / (double)numberOfTreesAll;
  scaleStop = 1.0 / (double)numberOfTreesStop;

  for(k = 0, entryCount = 0; k < h->tableSize; k++)	     
    {
      
      if(h->table[k] != NULL)
	{
	  entry *e = h->table[k];

	  do
	    {
	      c1 += e->bipNumber;
	      c2 += e->bipNumber2;
	      vect1[entryCount] = ((double)e->bipNumber) * scaleAll;
	      if(vect1[entryCount] > 0)
		bipAll++;
	      vect2[entryCount] = ((double)e->bipNumber2) * scaleStop;
	      if(vect2[entryCount] > 0)
		bipStop++;
	      fprintf(outf, "%f %f\n", vect1[entryCount], vect2[entryCount]);
	      entryCount++;
	      e = e->next;
	    }
	  while(e != NULL);
	}

     
    }
  
  printf("%ld %ld\n", c1, c2);

  assert(entryCount == h->entryCount);

  fclose(outf);

  p = testFreq(vect1, vect2, h->entryCount, tr);

  for(k = 0; k < h->entryCount; k++)
    {
      avg1 += vect1[k];
      avg2 += vect2[k];
    }

  avg1 /= ((double)h->entryCount);
  avg2 /= ((double)h->entryCount);
  
  
  printf("Average [%s]: %1.40f [%s]: %1.40f\n", bootStrapFileName, avg1, tree_file, avg2);
  printf("Pearson: %f Bipartitions in [%s]: %d Bipartitions in [%s]: %d Total Bipartitions: %d\n", p, bootStrapFileName, bipAll, tree_file, bipStop, h->entryCount);

  printf("\nFile containing pair-wise bipartition frequencies written to %s\n\n", bipFileName);
  
  freeBitVectors(bitVectors, 2 * tr->mxtips);
  free(bitVectors);
  freeHashTable(h);
  free(h);

  free(vect1);
  free(vect2);

  exit(0);
}

/*************************************************************************************************************/



/* 
   static unsigned int getMat(int i, int j, int n)
   {
   
   
   unsigned int index = i * (n - 1) - ((i * (i + 1)) / 2) + j - 1;
   
   
   
   return index;
   }
   
   static unsigned int getIndex(int i, int j, int n, unsigned int *mat)
   {
   
   
   int index = i * (n - 1) - ((i * (i + 1)) / 2) + j - 1;
   
   rfMat = (unsigned int*)calloc(((numberOfTrees * numberOfTrees - numberOfTrees) / 2), sizeof(unsigned int));
   return mat[index];
   }
   
   
   
   static void incMat(int i, int j, int n, unsigned int *mat)
   {
   int index = i * (n - 1) - ((i * (i + 1)) / 2) + j - 1;  
   
   mat[index] = mat[index] + 1;
   }
   
*/




void computeRF(tree *tr, char *bootStrapFileName)
{
  int  
    treeVectorLength, 
    numberOfTrees = 0, 
    i,
    j, 
    *rfMat,
    *wrfMat,
    *wrf2Mat,
    vLength; 

  unsigned int 
    k, 
    entryCount,    
    **bitVectors = initBitVector(tr, &vLength);
  
  char rfFileName[1024];

  boolean computeWRF = FALSE;

  double 
    maxRF, 
    avgRF,
    avgWRF,
    avgWRF2;

  FILE *outf;

  hashtable *h = initHashTable(tr->mxtips * 100, tr->mxtips);    

  INFILE = myfopen(bootStrapFileName, "r"); 
  numberOfTrees = countTrees(INFILE);        

  printf("\n\nFound %d trees in File %s\n\n", numberOfTrees, bootStrapFileName);
  
  if(numberOfTrees % MASK_LENGTH == 0)
    treeVectorLength = numberOfTrees / MASK_LENGTH;
  else
    treeVectorLength = 1 + (numberOfTrees / MASK_LENGTH);

  
  rfMat = (int*)calloc(numberOfTrees * numberOfTrees, sizeof(int));
  wrfMat = (int*)calloc(numberOfTrees * numberOfTrees, sizeof(int));
  wrf2Mat = (int*)calloc(numberOfTrees * numberOfTrees, sizeof(int));

  for(i = 0; i < numberOfTrees; i++)
    { 
      int bCounter = 0;
      
      int lcount = treeReadTopologyOnly(INFILE, tr, FALSE, FALSE, TRUE);
           
      assert(tr->mxtips == tr->ntips); 
      if(i == 0)
	{
	  assert(lcount == 0 || lcount == tr->mxtips - 3);
	  if(lcount == 0)
	    computeWRF = FALSE;
	  else
	    computeWRF = TRUE;
	}
      else
	{
	  if(computeWRF)
	    assert(lcount == tr->mxtips - 3);
	  else
	    assert(lcount == 0);
	}

      bitVectorInitravSpecial(bitVectors, tr->nodep[1]->back, tr->mxtips, vLength, h, i, BIPARTITIONS_RF, (branchInfo *)NULL,
			      &bCounter, treeVectorLength, FALSE, computeWRF);
     
      assert(bCounter == tr->mxtips - 3);      
    }
	  
  fclose(INFILE);   

  /*printf("Done 1 at %f\n", gettime() - masterTime);*/

  for(k = 0, entryCount = 0; k < h->tableSize; k++)	     
    {
      
      if(h->table[k] != NULL)
	{
	  entry *e = h->table[k];

	  do
	    {
	      unsigned int *vector = e->treeVector;
	      int *supportVector = e->supportVector;
	      
	      for(i = 0; i < numberOfTrees; i++)
		{
		  int
		    *r = &rfMat[i * numberOfTrees],
		    *w,
		    *w2,
		    val1 = ((vector[i / MASK_LENGTH] & mask32[i % MASK_LENGTH]) > 0),
		    s1;
		  
		  if(computeWRF)
		    {
		      
		      w = &wrfMat[i * numberOfTrees];
		      w2 = &wrf2Mat[i * numberOfTrees];
		      s1   = supportVector[i];
		    }

		  for(j = i + 1; j < numberOfTrees; j++)		  
		    {
		      if(computeWRF)
			{
			  int s2 = supportVector[j];
			  
			  if(s1 > 0 && s2 > 0)
			    w2[j] += ABS(s1 - s2);

			  if((s1 > 0 && s2 == 0) || (s2 > 0 && s1 == 0))
			    {
			      w2[j] += ABS(s1 - s2);
			      w[j] += ABS(s1 - s2);
			    }
			}

		      if(val1 + ((vector[j / MASK_LENGTH] & mask32[j % MASK_LENGTH]) > 0) == 1)		    	       		      
			r[j]++;
		    }
		}
	      
	      entryCount++;
	      e = e->next;
	    }
	  while(e != NULL);
	}

     
    }
  assert(entryCount == h->entryCount);  

  /*printf("Done 2 at %f\n", gettime() - masterTime);*/

  strcpy(rfFileName,         workdir);  
  strcat(rfFileName,         "RAxML_RF-Distances.");
  strcat(rfFileName,         run_id);

  outf = myfopen(rfFileName, "w");
  
  maxRF  = ((double)(2 * (tr->mxtips - 3)));
  avgRF  = 0.0;
  avgWRF = 0.0;
  avgWRF2 = 0.0;

  for(i = 0; i < numberOfTrees; i++)
    for(j = i + 1; j < numberOfTrees; j++)
      {
	int    rf = rfMat[i * numberOfTrees + j];
	double rrf = (double)rf / maxRF;
	if(computeWRF)
	  {
	    double wrf = wrfMat[i * numberOfTrees + j] / 100.0;
	    double rwrf = wrf / maxRF;
	    double wrf2 = wrf2Mat[i * numberOfTrees + j] / 100.0;
	    double rwrf2 = wrf2 / maxRF;
	
	    fprintf(outf, "%d %d: %d %f, %f %f, %f %f\n", i, j, rf, rrf, wrf, rwrf, wrf2, rwrf2);
	    avgWRF  += rwrf;
	    avgWRF2 += rwrf2; 
	  }
	else
	  fprintf(outf, "%d %d: %d %f\n", i, j, rf, rrf);
	/*printf("%d %d: %d %f\n", i, j, rf, rrf);*/
	avgRF += rrf;
      }
  
  fclose(outf);

  
  printf("\n\nAverage relative RF in this set: %f\n", avgRF / ((double)((numberOfTrees * numberOfTrees - numberOfTrees) / 2)));
  if(computeWRF)
    {
      printf("\n\nAverage relative WRF in this set: %f\n", avgWRF / ((double)((numberOfTrees * numberOfTrees - numberOfTrees) / 2)));
      printf("\n\nAverage relative WRF2 in this set: %f\n", avgWRF2 / ((double)((numberOfTrees * numberOfTrees - numberOfTrees) / 2)));
      printf("\nFile containing all %d pair-wise RF, WRF and WRF2 distances written to file %s\n\n", (numberOfTrees * numberOfTrees - numberOfTrees) / 2, rfFileName);
    }    
  else
    printf("\nFile containing all %d pair-wise RF distances written to file %s\n\n", (numberOfTrees * numberOfTrees - numberOfTrees) / 2, rfFileName);

  free(rfMat);
  free(wrfMat);
  free(wrf2Mat);
  freeBitVectors(bitVectors, 2 * tr->mxtips);
  free(bitVectors);
  freeHashTable(h);
  free(h);  

  exit(0);
}




/*************************************************************************************************************/




static void permute(unsigned int *perm, unsigned int n, long *seed)
{
  unsigned int  i, j, k;
 
  for (i = 0; i < n; i++) 
    {
      k =  (int)((double)(n - i) * randum(seed));
      j        = perm[i];    
      perm[i]     = perm[i + k];
      perm[i + k] = j; 
      /*assert(i + k < n);*/
    }
}


/* #ifdef _USE_PTHREADS
  
void threadComputeAverage(tree *tr, int tid)
{
  int
    i,
    n = tr->numberOfBipartitions,
    lower,
    upper;
  double
    *vect1 = tr->v1,
    *vect2 = tr->v2,
    avg1 = 0.0,
    avg2 = 0.0;

  if(tid == 0)
    {
      lower = 0;
      upper = n/2;
    }
  else
    {
      lower = n/2;
      upper = n;
    } 

  for(i = lower; i < upper; i++)
    {   
      avg1 += vect1[i];
      avg2 += vect2[i];	
    }
  reductionBuffer[tid]    = avg1;
  reductionBufferTwo[tid] = avg2;
}

void threadComputePearson(tree *tr, int tid)
{
  int
    i,
    n = tr->numberOfBipartitions,
    lower,
    upper;
  double
    avg1 = tr->avg1,
    avg2 = tr->avg2,
    *vect1 = tr->v1,
    *vect2 = tr->v2,
    sum_xy = 0.0,
    sum_x = 0.0,
    sum_y = 0.0;

   if(tid == 0)
    {
      lower = 0;
      upper = n/2;
    }
  else
    {
      lower = n/2;
      upper = n;
    } 

  for(i = lower; i < upper; i++)
    {      
      sum_xy += ((vect1[i] - avg1) * (vect2[i] - avg2));
      sum_x  += ((vect1[i] - avg1) * (vect1[i] - avg1));
      sum_y  += ((vect2[i] - avg2) * (vect2[i] - avg2));    
    }
  reductionBuffer[tid]      = sum_xy;
  reductionBufferTwo[tid]   = sum_x;
  reductionBufferThree[tid] = sum_y;
  }

  #endif*/


static double testFreq(double *vect1, double *vect2, int n, tree *tr)
{
  int i;
  double
    avg1 = 0.0, 
    avg2 = 0.0,
    sum_xy = 0.0, 
    sum_x  = 0.0, 
    sum_y  = 0.0,
    corr   = 0.0;

  /*#ifdef _USE_PTHREADS
  tr->numberOfBipartitions = n;
  tr->v1 = vect1;
  tr->v2 = vect2;
  
  masterBarrier(THREAD_COMPUTE_AVERAGE, tr);

  for(i = 0; i < NumberOfThreads; i++)
    {
      avg1 += reductionBuffer[i];
      avg2 += reductionBufferTwo[i]; 
    }

    #else*/
  for(i = 0; i < n; i++)
    {	     
      avg1 += vect1[i];
      avg2 += vect2[i];
    }
  /*#endif*/
      
  avg1 /= ((double)n);
  avg2 /= ((double)n); 
    

  /*#ifdef _USE_PTHREADS
  tr->avg1 = avg1;
  tr->avg2 = avg2;

  masterBarrier(THREAD_COMPUTE_PEARSON, tr);

  for(i = 0; i < NumberOfThreads; i++)
    {
      sum_xy += reductionBuffer[i];
      sum_x  += reductionBufferTwo[i]; 
      sum_y  += reductionBufferThree[i];
    }
    #else*/
  for(i = 0; i < n; i++)
    {
      sum_xy += ((vect1[i] - avg1) * (vect2[i] - avg2));
      sum_x  += ((vect1[i] - avg1) * (vect1[i] - avg1));
      sum_y  += ((vect2[i] - avg2) * (vect2[i] - avg2));	 
    }
  /*#endif*/

  corr = sum_xy / (sqrt(sum_x) * sqrt(sum_y));
   
#ifndef WIN32
  if(isnan(corr))
    {
      printf("Numerical Error pearson correlation is not a number\n");
      assert(0);
    }
#endif

  return corr;
}

















/*#ifdef _USE_PTHREADS
void threadMakeVector(tree *tr, int tid)
{
  double *vect1 = tr->v1;
  double *vect2 = tr->v2;
  BL *b         = tr->b;
  int *perm     = tr->perm;
  int numberOfTrees = tr->numberOfBootstopTrees;
  int lower, upper, l, j;
  unsigned char *set;

  if(tid == 0)
    {
      lower = 0;
      upper = b->count/2;
    }
  else
    {
      lower = b->count/2;
      upper = b->count;
    }   
	      
  for(j = lower; j < upper; j++)
    {		       
      set = b->b[j].isSet;	      
      
      for(l = 0; l < numberOfTrees; l++)
	{			     
	  if(get_bit(set,l))
	    {
	      if(perm[l] % 2 == 0)
		vect1[j] = vect1[j] + 1.0;
	      else			
		vect2[j] = vect2[j] + 1.0;
	    }			     
	}		
    }	
}
#endif*/

static double frequencyCriterion(int numberOfTrees, hashtable *h, int *countBetter, tree *tr)
{
  int 
    k, 
    l;
    
  long 
    seed = 12345;

  double 
    t,
    result, 
    avg = 0, 
    *vect1, 
    *vect2; 

  unsigned int
    *perm =  (unsigned int *)malloc(sizeof(unsigned int) * numberOfTrees),
    j;

  assert(*countBetter == 0);
	  
 
	  
  for(j = 0; j < (unsigned int)numberOfTrees; j++)
    perm[j] = j;
	  
  for(k = 0; k < BOOTSTOP_PERMUTATIONS; k++)
    {   	      		      	      
      unsigned int entryCount = 0;

      permute(perm, numberOfTrees, &seed);
      
      t = gettime();

      vect1 = (double *)calloc(h->entryCount, sizeof(double));
      vect2 = (double *)calloc(h->entryCount, sizeof(double));	     

      /*#ifdef _USE_PTHREADS
      tr->v1 = vect1;
      tr->v2 = vect2;
      tr->b = b;
      tr->numberOfBootstopTrees = numberOfTrees;
      tr->perm = perm;

      masterBarrier(THREAD_MAKE_VECTORS, tr);      
      #else	*/      
      
      for(j = 0; j < h->tableSize; j++)
	{		
	  if(h->table[j] != NULL)
	    {		
	      entry *e = h->table[j];
	      
	      do
		{
		  unsigned int *set = e->treeVector;       
		  
		  for(l = 0; l < numberOfTrees; l++)
		    {			     
		      if((set[l / MASK_LENGTH] != 0) && (set[l / MASK_LENGTH] & mask32[l % MASK_LENGTH]))
			{
			  if(perm[l] % 2 == 0)
			    vect1[entryCount] = vect1[entryCount] + 1.0;
			  else			
			    vect2[entryCount] = vect2[entryCount] + 1.0;
			}
		    }
		  entryCount++;
		  e = e->next;
		}
	      while(e != NULL);
	    }			     
	}		    	
      /*#endif*/
      
      
      
      assert(entryCount == h->entryCount);
      
      

      result = testFreq(vect1, vect2, entryCount, tr);
	  
          
      
      if(result >= FC_LOWER)
	*countBetter = *countBetter + 1;
	      
      avg += result;
	      
      free(vect1);		  
      free(vect2);		 
    }

  free(perm);
	  
  avg /= BOOTSTOP_PERMUTATIONS;
	
      

  return avg;
}




static double wcCriterion(int numberOfTrees, hashtable *h, int *countBetter, double *wrf_thresh_avg, double *wrf_avg, tree *tr)
{
  int 
    k, 
    l,   
    wrf,
    mr_thresh = ((double)numberOfTrees/4.0);
   
  unsigned int 
    *perm =  (unsigned int *)malloc(sizeof(unsigned int) * numberOfTrees),
    j;

  long seed = 12345;  
  
  double 
    wrf_thresh = 0.0,
    pct_avg = 0.0;

  assert(*countBetter == 0 && *wrf_thresh_avg == 0.0 && *wrf_avg == 0.0);
	   	  
  for(j = 0; j < (unsigned int)numberOfTrees; j++)
    perm[j] = j;
	  
  for(k = 0; k < BOOTSTOP_PERMUTATIONS; k++)
    {   	      		           
      int mcnt1 = 0;			  
      int mcnt2 = 0;
      unsigned int entryCount = 0;
      
      wrf = 0;
      	      
      permute(perm, numberOfTrees, &seed);      
    
      for(j = 0; j < h->tableSize; j++)
	{		
	  if(h->table[j] != NULL)
	    {
	      entry *e = h->table[j];

	        do
		  {
		    int cnt1 = 0;
		    int cnt2 = 0;

		    unsigned int *set = e->treeVector;
		    for(l = 0; l < numberOfTrees; l++)
		      {			     
			if((set[l / MASK_LENGTH] != 0) && (set[l / MASK_LENGTH] & mask32[l % MASK_LENGTH]))
			  {			    
			    if(perm[l] % 2 == 0)
			      cnt1++;
			    else			
			      cnt2++;
			  }			     
		      }
		    
		    if(cnt1 <= mr_thresh)			      
		      cnt1 = 0;
		       
		    if(cnt2 <= mr_thresh)	    
		      cnt2 = 0;

		    if(cnt1 > 0)			      
		      mcnt1++;

		    if(cnt2 > 0)			      
		      mcnt2++;

		    wrf += ((cnt1 > cnt2) ? cnt1 - cnt2 : cnt2 - cnt1);

		    entryCount++;
		    e = e->next;
		  }
		while(e != NULL);

	    }	  	  
	  	
	}	
      
      assert(entryCount == h->entryCount);

      /* 
	 wrf_thresh is the 'custom' threshold computed for this pair
	 of majority rules trees (i.e. one of the BS_PERMS splits),
	 and simply takes into account the resolution of the two trees
      */

      wrf_thresh = (tr->wcThreshold) * ( ((((double)numberOfTrees/2.0) * (double)mcnt1)) + ((((double)numberOfTrees/2.0) * (double)mcnt2)) );      
      
      /*
	we count this random split as 'succeeding' when
	 the wrf between maj rules trees is exceeded
	 by its custom threshold
      */

      if((double)wrf <= wrf_thresh)			        
	*countBetter = *countBetter + 1;

      /* 
	 here we accumulate outcomes and thresholds, because
	 we're not going to stop until the avg dist is less
	 than the avg threshold
      */
      pct_avg += ((double)wrf / ( ((double)numberOfTrees/2.0 * (double)mcnt1) + ((double)numberOfTrees/2.0 * (double)mcnt2) )) * 100.0;
      *wrf_avg += (double)wrf;
      *wrf_thresh_avg += wrf_thresh;
    }
 
  free(perm);

  pct_avg /= (double)BOOTSTOP_PERMUTATIONS; 
  *wrf_avg /= (double)BOOTSTOP_PERMUTATIONS; 
  *wrf_thresh_avg /= (double)BOOTSTOP_PERMUTATIONS;   
  /*printf("%d \t\t %f \t\t %d \t\t\t\t %f\n", numberOfTrees, *wrf_avg, *countBetter, *wrf_thresh_avg);	  	      */
  return pct_avg; 
}	  






void computeBootStopOnly(tree *tr, char *bootStrapFileName)
{
  int numberOfTrees = 0, i;
  boolean stop = FALSE;
  double avg;
  int checkEvery;
  int treesAdded = 0;
  hashtable *h = initHashTable(tr->mxtips * FC_INIT * 10, tr->mxtips); 
  int 
    treeVectorLength, 
    vectorLength;
  unsigned int **bitVectors = initBitVector(tr, &vectorLength);

  double tt = 0.0, tc = 0.0;

  assert((FC_SPACING % 2 == 0) && (FC_THRESHOLD < BOOTSTOP_PERMUTATIONS));

  INFILE = myfopen(bootStrapFileName, "r");       
  numberOfTrees = countTrees(INFILE); 
  
  printf("\n\nFound %d trees in File %s\n\n", numberOfTrees, bootStrapFileName);
  
  assert(sizeof(unsigned char) == 1);
  
  if(numberOfTrees % MASK_LENGTH == 0)
    treeVectorLength = numberOfTrees / MASK_LENGTH;
  else
    treeVectorLength = 1 + (numberOfTrees / MASK_LENGTH);  
 
  checkEvery = FC_SPACING;
        
  switch(tr->bootStopCriterion)
    {
    case FREQUENCY_STOP:
      printf("# Trees \t Average Pearson Coefficient \t # Permutations: pearson >= %f\n", 
	     FC_LOWER);
      break;
    case WC_STOP:
      printf("# Trees \t Avg WRF in %s \t # Perms: wrf <= %1.2f %s\n","%", 100.0 * tr->wcThreshold, "%");
      break;
    default:
      printf("%d \n", tr->bootStopCriterion);
      assert(0);
    }
  
  for(i = 1; i <= numberOfTrees && !stop; i++)
    {                  
      int bCount = 0;      
      double t = gettime();
      treeReadTopologyOnly(INFILE, tr, FALSE, FALSE, FALSE);	  	     
      assert(tr->mxtips == tr->ntips);
      
      bitVectorInitravSpecial(bitVectors, tr->nodep[1]->back, tr->mxtips, vectorLength, h, (i - 1), BIPARTITIONS_BOOTSTOP, (branchInfo *)NULL,
			      &bCount, treeVectorLength, FALSE, FALSE);
      tt += gettime() - t;
      assert(bCount == tr->mxtips - 3);
                 
      treesAdded++;	
            
      if(i > START_BSTOP_TEST && i % checkEvery == 0)
	{ 
	  int countBetter = 0;
	  
	  t = gettime();
	  switch(tr->bootStopCriterion)
	    {
	    case FREQUENCY_STOP:
	      avg = frequencyCriterion(i, h, &countBetter, tr);	  	  	  
	      printf("%d \t\t\t %f \t\t\t\t %d\n", i, avg, countBetter);
	  	  
	      stop = (countBetter >= FC_THRESHOLD && avg >= FC_LOWER);	  	 
	      break;
	    case WC_STOP:
	      {
		double 
		  wrf_thresh_avg = 0.0,
		  wrf_avg = 0.0;
		avg = wcCriterion(i, h, &countBetter, &wrf_thresh_avg, &wrf_avg, tr);
		printf("%d \t\t %1.2f \t\t\t %d\n", i, avg, countBetter);	       
		
		stop = (countBetter >= FC_THRESHOLD && wrf_avg <= wrf_thresh_avg);
	      }
	      break;
	    default:
	      assert(0);
	    }
	  tc += gettime() - t;
	}	 	   
      
    }
  
  /* 
     printf("Time traverse add %f compute %f\n", tt, tc);
     actual computation of respective statistics consumes the by far largest amount of time 
  */

  if(stop)              
    printf("Converged after %d replicates\n", treesAdded);           
  else    
    printf("Bootstopping test did not converge after %d trees\n", treesAdded);
     
  freeBitVectors(bitVectors, 2 * tr->mxtips);
  free(bitVectors);
  freeHashTable(h);
  free(h);
  
  fclose(INFILE);
  exit(0);
}

boolean bootStop(tree *tr, hashtable *h, int numberOfTrees, double *pearsonAverage, unsigned int **bitVectors, int treeVectorLength, int vectorLength)
{
  int 
    n = numberOfTrees + 1,
    bCount = 0;

  assert((FC_SPACING % 2 == 0) && (FC_THRESHOLD < BOOTSTOP_PERMUTATIONS));
  assert(tr->mxtips == tr->rdta->numsp);

  bitVectorInitravSpecial(bitVectors, tr->nodep[1]->back, tr->mxtips, vectorLength, h, numberOfTrees, BIPARTITIONS_BOOTSTOP, (branchInfo *)NULL,
			  &bCount, treeVectorLength, FALSE, FALSE);
  assert(bCount == tr->mxtips - 3); 

  if((n > START_BSTOP_TEST) && (n % FC_SPACING == 0))
    {     
      int countBetter = 0;

      switch(tr->bootStopCriterion)
	{
	case FREQUENCY_STOP:
	  *pearsonAverage = frequencyCriterion(n, h, &countBetter, tr);	  	        	       

	  printf("WC %d %f %d\n", n, *pearsonAverage, countBetter);    

	  if(countBetter >= FC_THRESHOLD)
	    return TRUE;
	  else
	    return FALSE;
	  break;
	case WC_STOP:
	 {	   
	   double 
	     wrf_thresh_avg = 0.0,
	     wrf_avg = 0.0;
	   
	   *pearsonAverage = wcCriterion(n, h, &countBetter, &wrf_thresh_avg, &wrf_avg, tr);
		  
	   printf("WC %d %f %d\n", n, *pearsonAverage, countBetter);	 
	  
	   if(countBetter >= FC_THRESHOLD && wrf_avg <= wrf_thresh_avg)
	     return TRUE;
	   else
	     return FALSE;
	 } 
	default:
	  assert(0);
	}
    }
  else
    return FALSE;
}


/******************************************************** MRP *******************************************************************************/

void getMRP_Bipartitions(nodeptr p, tree *tr, int *countBips, unsigned int **bitVectors, int vectorLength, int *countColumns)
{
  if(isTip(p->number, tr->mxtips))          
    return;
  else
    {
      nodeptr q = p->next;
                 
            

      if(!(isTip(p->back->number, tr->mxtips)))
	{
	  int i, c = 0;
	  unsigned int *left  = bitVectors[p->number];	       	 	  	 
	  unsigned int *right = bitVectors[p->back->number];
	  int totalLength = vectorLength * MASK_LENGTH;

	  if(!p->x)
	    newviewBipartitions(bitVectors, p, tr->mxtips, vectorLength);
	  assert(p->x);
	  
	  if(!p->back->x)
	    newviewBipartitions(bitVectors, p->back, tr->mxtips, vectorLength);
	  assert(p->back->x);

	  for(i = 0; i < totalLength; i++)
	    {
	      unsigned int le = left[i / MASK_LENGTH];
	      unsigned int ri = right[i / MASK_LENGTH];

	      if(le & mask32[i % MASK_LENGTH])
		{
		  tr->yVector[i + 1][*countColumns] = 1;
		  c++;
		}
	      if(ri & mask32[i % MASK_LENGTH])
		{
		   tr->yVector[i + 1][*countColumns] = 2;
		   c++;
		}
	    }
	 
	  assert(c == tr->ntips);
	  
	  *countColumns = *countColumns + 1;
	  *countBips = *countBips + 1;
	} 
      
      do 
	{
	  getMRP_Bipartitions(q->back, tr, countBips, bitVectors, vectorLength, countColumns);	  
	  q = q->next;
	}
      while(q != p); 
          
     
      return;
    }
}





void encodeMRP(tree *tr, rawdata *rdta)
{
  FILE *f = myfopen(bootStrapFile, "r");
  char 
    **nameList,
    mrpFileName[1024];
  int
    countColumns = 0,
    matrixLength = 0,
    tipCount = 0,
    tips,
    inter,
    taxaSize = 1000,
    taxaCount = 0,
    numberOfTrees,
    c,
    i,
    j,
    vLength;
  char buffer[nmlngth + 2];
  nodeptr p, q, p0;
  unsigned int **bitVectors;   
 

  nameList = (char**)malloc(sizeof(char*) * taxaSize);  

  printf("here\n");

  while((c = fgetc(f)) != EOF)
    {
      /*printf("%c\n", c);*/
      if(c == '(' || c == ',')
	{
	  c = fgetc(f);
	  if(c ==  '(' || c == ',')
	    ungetc(c, f);
	  else
	    {
	      int k = 0;	      
	      boolean found = FALSE;
	      /*ungetc(c, f);*/
	      do
		{
		  buffer[k++] = c;
		  c = fgetc(f);
		}
	      while(c != ':' && c != ')' && c != ',');
	      buffer[k] = '\0';

	      tipCount++;

	      for(i = 0; i < taxaCount && !found; i++)
		{
		  if(strcmp(buffer, nameList[i]) == 0)
		    found = TRUE;
		}

	      if(!found)
		{
		  if(taxaCount == taxaSize)
		    {
		      char **b = (char**)malloc(sizeof(char*) * taxaSize * 2);
		      for(i = 0; i < taxaSize; i++)
			b[i] = nameList[i];
		      taxaSize *= 2;
		      free(nameList);
		      nameList = b;
		    }
		  nameList[taxaCount] = (char*)malloc(sizeof(char) * (strlen(buffer) + 1));
		  strcpy(nameList[taxaCount], buffer);
		  printf("%s\n", nameList[taxaCount]);
		  taxaCount++;
		}

	      /*printf("%s\n", buffer);*/
	      ungetc(c, f);
	    }
	}
      else
	if(c == ';')
	  {
	    matrixLength += (tipCount - 3);
	    /*printf("Tree has %d tips\n", tipCount);*/
	    tipCount = 0;
	  }
    }
 
  printf("found %d taxa in tree collection total of %d MRP columns\n", taxaCount, matrixLength);

  tips = tr->mxtips = taxaCount;
  inter = tr->mxtips - 1;
  tr->rdta        = rdta;
  tr->rdta->numsp = tips;
  tr->start       = (node *) NULL;
  tr->ntips       = 0;
  tr->nextnode    = 0;  

  tr->yVector = (unsigned char**)malloc(sizeof(unsigned char*) * (tr->mxtips + 1));

  for(i = 1; i <= tr->mxtips; i++)
    {
      tr->yVector[i] = (unsigned char*)malloc(sizeof(unsigned char) * matrixLength);
      for(j = 0; j < matrixLength; j++)
	tr->yVector[i][j] = UNDETERMINED_BINARY; 
    }

 
  bitVectors = initBitVector(tr, &vLength);

  tr->nameList = (char **)malloc(sizeof(char *) * (tips + 1));  
  for(i = 1; i <= tips; i++)
    tr->nameList[i] = nameList[i - 1];
  
  free(nameList);

  if(!(p0 = (nodeptr) malloc((tips + 3*inter) * sizeof(node)))) 
    assert(0);
     

  if(!(tr->nodep = (nodeptr *) malloc((2*tr->mxtips) * sizeof(nodeptr)))) 
    assert(0);
    
  tr->nodep[0] = (node *) NULL;

  for (i = 1; i <= tips; i++) 
    {
      p = p0++;
     
      p->x      =  0;    
      p->number =  i;
      p->next   =  p;
      p->back   = (node *)NULL; 
      p->bInf   = (branchInfo *)NULL;
      tr->nodep[i] = p;
    }

  for (i = tips + 1; i <= tips + inter; i++) 
    {
      q = (node *) NULL;
      for (j = 1; j <= 3; j++) 
	{
	  p = p0++;
	  if(j == 1)
	    p->x      =  1;
	  else
	    p->x      =  0;	  
	  p->number = i;
	  p->next   = q;
	  p->bInf   = (branchInfo *)NULL;
	  p->back   = (node *) NULL;
	  q = p;
	}
      p->next->next->next = p;
      tr->nodep[i] = p;
    }

  tr->start       = (node *) NULL;
  tr->ntips       = 0;
  tr->nextnode    = 0;   

  tr->numBranches = 1; 
 
  
  rewind(f);
  numberOfTrees = countTrees(f);
  
  printf("Found %d trees in tree collection\n", numberOfTrees);

  for(i = 0; i < numberOfTrees; i++)
    {
      int countBips = 0;
      int countBranches = 0;

      treeReadTopologyOnly(f, tr, FALSE, TRUE, FALSE);
      printf("Succesfully read tree topol %d\n", i);
      tr->start = findAnyTip(tr->start, tr->ntips);
      assert(isTip(tr->start->number, tr->ntips));
      {
	int k, j;
	for(k = tr->mxtips + 1; k < 2 * tr->mxtips; k++)
	  {
	    for(j = 0; j < vLength; j++)	      
	      bitVectors[k][j] = 0;	      
	  }
      }
      
      bitVectorInitravSpecial(bitVectors, tr->start->back, tr->mxtips, vLength, (hashtable *)NULL, -1, 0, (branchInfo *)NULL, 
			      &countBranches, 0, TRUE, FALSE);
      assert(countBranches == (tr->ntips - 3));
      printf("Branches %d\n", countBranches);
      
      getMRP_Bipartitions(tr->start->back, tr, &countBips, bitVectors, vLength, &countColumns); 
      assert(countBips == (tr->ntips - 3));
    }
  fclose(f);

  assert(countColumns == matrixLength);

  strcpy(mrpFileName, workdir); 
  strcat(mrpFileName, "RAxML_mrpAlignment.");
  strcat(mrpFileName, run_id);    

  f = myfopen(mrpFileName, "w");

  fprintf(f, "%d %d\n", tr->mxtips, countColumns);

  for(i = 1; i <= tr->mxtips; i++)
    {
      fprintf(f, "%s ", tr->nameList[i]);
      for(j = 0; j < matrixLength; j++)
	fprintf(f, "%c", inverseMeaningBINARY[tr->yVector[i][j]]);
      fprintf(f, "\n");
    }
  
  fclose(f);
  exit(1);
}

