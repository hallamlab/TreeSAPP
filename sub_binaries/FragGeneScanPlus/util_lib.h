#ifndef __UTIL_LIB_H
#define __UTIL_LIB_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "hmm.h"

double **dmatrix(int num_col);
int **imatrix(int num_col);
int *ivector(int nh);

void free_dmatrix(double **m);
void free_imatrix(int **m);

int tr2int (char *nt);
int nt2int (char nt);
int nt2int_rc (char nt);

int trinucleotide (char a, char b, char c);
void get_protein(char *dna, char *protein, int strand);
void print_usage();

typedef struct q {
    struct q* next;
    thread_data* td;
    unsigned int buffer;
} QUEUE;

QUEUE* q_empty_head;
QUEUE* q_empty_tail;
QUEUE* q_done_head;
QUEUE* q_done_tail;


void printq(unsigned int which);
void enq(thread_data* td, unsigned int buffer, unsigned int which);
QUEUE* deq(unsigned int which);

void cutnpaste_q(QUEUE** dest, unsigned int which);

void stopMemset(char* ptr, int length);

#endif
