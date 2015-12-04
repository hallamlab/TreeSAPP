#ifndef DYNAMITEkmer_assemblyHEADERFILE
#define DYNAMITEkmer_assemblyHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "kmer_index_interface.h"
#include "singleseqspace.h"
#include "assembly.h"

#define KMER_ASSEMBLY_NEXT_TANGLED 2
#define KMER_ASSEMBLY_PREV_TANGLED 4


typedef struct KmerAssemblyLink {
  struct KmerAssemblyNode * prev;
  struct KmerAssemblyNode * next;
  char base;
  long * sequence_label;
  int sequence_label_len;
  int sequence_label_maxlen;
  char state;
} KmerAssemblyLink;


typedef struct KmerAssemblyNode {
  kmer_t number;
  struct KmerAssemblyLink ** prev;
  int prev_len;
  int prev_maxlen;
  struct KmerAssemblyLink ** next;
  int next_len;
  int next_maxlen;
  struct KmerAssemblyNode * node_chain;
} KmerAssemblyNode;

typedef struct KmerAssemblyIndex {
  KmerIndexInterface * kii;
  SinglePosSpace * sps;
} KmerAssemblyIndex;


#define KmerAssemblyNode_LINK_START  8
#define KmerAssemblyNode_LINK_LINEAR 64

#define KmerAssemblyLink_LABEL_START  2
#define KmerAssemblyLink_LABEL_LINEAR 8




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
void Wise2_show_extensive_stats_KmerAssemblyIndex(KmerAssemblyIndex * kai,FILE * ofp);
#define show_extensive_stats_KmerAssemblyIndex Wise2_show_extensive_stats_KmerAssemblyIndex
void Wise2_add_AssemblySequence_KmerAssemblyIndex(KmerAssemblyIndex * kai,AssemblySequence * aseq,int report);
#define add_AssemblySequence_KmerAssemblyIndex Wise2_add_AssemblySequence_KmerAssemblyIndex
boolean Wise2_store_KmerAssemblyLink_KmerAssemblyIndex(KmerAssemblyIndex * kai,kmer_t prev_number,kmer_t next_number,char base,long label);
#define store_KmerAssemblyLink_KmerAssemblyIndex Wise2_store_KmerAssemblyLink_KmerAssemblyIndex
void Wise2_show_KmerAssemblyIndex(KmerAssemblyIndex * kai,FILE * ofp);
#define show_KmerAssemblyIndex Wise2_show_KmerAssemblyIndex
void Wise2_show_KmerAssemblyNode(KmerAssemblyNode * node,int kmer_size,int level,FILE * ofp);
#define show_KmerAssemblyNode Wise2_show_KmerAssemblyNode
KmerAssemblyIndex * Wise2_new_KmerAssemblyIndex(KmerIndexInterface * kii,SinglePosSpace * sps);
#define new_KmerAssemblyIndex Wise2_new_KmerAssemblyIndex
KmerAssemblyLink * Wise2_new_KmerAssemblyLink(char base);
#define new_KmerAssemblyLink Wise2_new_KmerAssemblyLink
void Wise2_remove_sequence_label_KmerAssemblyLink(KmerAssemblyLink * kal,long label);
#define remove_sequence_label_KmerAssemblyLink Wise2_remove_sequence_label_KmerAssemblyLink
void Wise2_add_sequence_label_KmerAssemblyLink(KmerAssemblyLink * kal,long label);
#define add_sequence_label_KmerAssemblyLink Wise2_add_sequence_label_KmerAssemblyLink
void Wise2_detach_KmerAssemblyLink(KmerAssemblyIndex * kai,KmerAssemblyLink * link);
#define detach_KmerAssemblyLink Wise2_detach_KmerAssemblyLink
void Wise2_remove_next_KmerAssemblyNode(KmerAssemblyNode * node,KmerAssemblyLink * next);
#define remove_next_KmerAssemblyNode Wise2_remove_next_KmerAssemblyNode
void Wise2_add_next_KmerAssemblyNode(KmerAssemblyNode * kan,KmerAssemblyLink * kal);
#define add_next_KmerAssemblyNode Wise2_add_next_KmerAssemblyNode
void Wise2_remove_prev_KmerAssemblyNode(KmerAssemblyNode * node,KmerAssemblyLink * prev);
#define remove_prev_KmerAssemblyNode Wise2_remove_prev_KmerAssemblyNode
void Wise2_add_prev_KmerAssemblyNode(KmerAssemblyNode * kan,KmerAssemblyLink * kal);
#define add_prev_KmerAssemblyNode Wise2_add_prev_KmerAssemblyNode
KmerAssemblyNode * Wise2_new_KmerAssemblyNode(kmer_t number);
#define new_KmerAssemblyNode Wise2_new_KmerAssemblyNode


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/

#ifdef _cplusplus
}
#endif

#endif
