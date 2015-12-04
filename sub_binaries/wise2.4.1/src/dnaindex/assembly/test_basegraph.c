#include "basegraph.h"
#include "kmer_glib_index.h"



int test_basic_data(void)
{
  int i;

  BaseLink * bl;

  BaseNode * bn_a;
  BaseNode * bn_b;
  BaseNode * st = NULL;
  
  bn_a = new_BaseNode(&st,456);
  bn_b = new_BaseNode(&bn_a->next_node,67);

  for(i=0;i<56;i++) {
    bl = new_BaseLink(bn_a,bn_b,1);
  }

  for(i=0;i<412;i++) {
    add_Sequence_BaseLink(bl,78);
  }
  
  assert(bn_a->next_node == bn_b);
  assert(bn_a->link_len == 56);
  assert(bn_b->link_len == 56);
  assert(bl->sequence_label_len == 412);
  assert(bl->sequence_label[77] == 78);

  fprintf(stdout,"basic data tests successful\n");

  return 0;

}

int test_graph(void)
{
  KmerIndexInterface * kii;
  BaseGraph * bg;
  BaseNode * bn;
  BaseLink * bl_1;
  BaseLink * bl_2;

  kii = new_interface_KmerGlibIndex(9);

  bg = new_BaseGraph(kii);

  insert_seqpos_BaseGraph(bg,1445,678,1,56);
  insert_seqpos_BaseGraph(bg,1445,678,1,57);
  insert_seqpos_BaseGraph(bg,678,43677,0,78);
  
  fprintf(stdout,"basic graph tests successful\n");


  bn = retrieve_BaseNode(bg,1445);
  assert(bn->link_len == 1);
  assert(bn->link[0]->sequence_label_len == 2);
  assert(bn->link[0]->twist == 1);

  bn = retrieve_BaseNode(bg,678);
  assert(bn->link_len == 2);

  bl_1 = bn->link[0];
  assert(next_link_if_unique_BaseLink(bl_1,bn,&bl_2) == 1);
  assert( (bl_1->sequence_label_len == 2 && bl_2->sequence_label_len == 1) ||
	  (bl_1->sequence_label_len == 1 && bl_2->sequence_label_len == 2) );


  fprintf(stdout,"link following tests successful\n");

  return 0;
}


int main(int argc,char ** argv)
{

  test_basic_data();

  test_graph();

  return 0;
}
