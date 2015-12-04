
#include "intallocator.h"

double a_random_0_to_1(void)
{
  double ret;
  long rnd;

  rnd = random();

  ret = (double)rnd / (double)LONG_MAX;

/*  fprintf(stderr," got... %D %f\n",rnd,ret);*/

  return ret;
}


int a_random_integer(int l)
{
  double rand;
  
  rand = a_random_0_to_1();
  
  return (int) (l*rand);
}


int main(int argc,char ** argv)
{
  int * alloc_array[4024];
  int i;
  int rnd;
  IntAllocator * inta;

  inta = new_IntAllocator(4);
  srandom(10);

  for(i=0;i<4024;i++) {
    alloc_array[i] = NULL;
  }

  for(i=0;i<7000;i++) {
    rnd = a_random_integer(4024);
    if( alloc_array[rnd] != NULL ) {
      fprintf(stderr,"Freeing position %d, (count %d)\n",rnd,i);
      free_intarray_IntAllocator(inta,alloc_array[rnd]);
      alloc_array[rnd] = NULL;
    }
    
    rnd = a_random_integer(4024);
    if( alloc_array[rnd] == NULL ) {
      fprintf(stderr,"allocating position %d, (count %d)\n",rnd,i);
      alloc_array[rnd] = alloc_intarray_IntAllocator(inta);
    }
  
  }

  show_allocator_status_IntAllocator(inta,stdout);

  return 0;


}
