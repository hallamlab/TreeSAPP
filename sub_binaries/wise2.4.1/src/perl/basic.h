
#include <stdio.h>

typedef char aa;
typedef int base;
typedef double Probability;
typedef double Bits;
typedef int Score;
typedef int codon;
typedef int boolean;

#define WISE2_FATAL    1
#define WISE2_WARNING  2
#define WISE2_INFO     8
#define WISE2_REPORT   16


char * Wise2_stringalloc(char *);
void   Wise2_error_off(int type);
void   Wise2_error_on(int type);
