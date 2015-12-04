#ifndef DYNAMITEkmer_assembly_errorHEADERFILE
#define DYNAMITEkmer_assembly_errorHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "kmer_assembly.h"
#include "kmer_assembly_untangler.h"



    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
boolean Wise2_lift_indel_error_KmerAssembly(KmerAssemblyIndex * kai,KmerAssemblyPath * real,KmerAssemblyPath * error,long int * start_labels,int label_len,int error_position,int error_len);
#define lift_indel_error_KmerAssembly Wise2_lift_indel_error_KmerAssembly
boolean Wise2_lift_forward_error_KmerAssembly(KmerAssemblyIndex * kai,KmerAssemblyPath * real,KmerAssemblyPath * error,long int * start_labels,int label_len);
#define lift_forward_error_KmerAssembly Wise2_lift_forward_error_KmerAssembly
int Wise2_mark_tangles_KmerAssembly(KmerAssemblyIndex * kai);
#define mark_tangles_KmerAssembly Wise2_mark_tangles_KmerAssembly
int Wise2_resolve_forward_errors_KmerAssembly(KmerAssemblyIndex * kai,int depth,int verbose,int max_path_enum);
#define resolve_forward_errors_KmerAssembly Wise2_resolve_forward_errors_KmerAssembly
void Wise2_remove_errors_KmerAssemblyIndex(KmerAssemblyIndex * kai,int max_error_depth);
#define remove_errors_KmerAssemblyIndex Wise2_remove_errors_KmerAssemblyIndex
int Wise2_attempt_resolve_forward_error_KmerAssembly(KmerAssemblyIndex * kai,KmerAssemblyNode * node,int depth,int max_path_enum);
#define attempt_resolve_forward_error_KmerAssembly Wise2_attempt_resolve_forward_error_KmerAssembly
boolean Wise2_find_indel_path_KmerAssembly(KmerAssemblyLink *real,KmerAssemblyLink * error,int delete_length,int proposed_position,KmerAssemblyPath * final_real,KmerAssemblyPath * final_error,int max_search,int max_path_enum);
#define find_indel_path_KmerAssembly Wise2_find_indel_path_KmerAssembly
boolean Wise2_extend_indel_path_KmerAssembly(KmerAssemblyLink* real,KmerAssemblyLink * error,int real_pos,int error_pos,int delete_length,int proposed_position,KmerAssemblyPath * final_real,KmerAssemblyPath * final_error,int max_search,int *current_path,int max_path_enum);
#define extend_indel_path_KmerAssembly Wise2_extend_indel_path_KmerAssembly
boolean Wise2_find_error_path_KmerAssembly(KmerAssemblyLink * start_real,KmerAssemblyLink * start_error,char proposed_fix,int proposed_fix_position,KmerAssemblyPath * final_real,KmerAssemblyPath * final_error,int max_search,int max_path_enum);
#define find_error_path_KmerAssembly Wise2_find_error_path_KmerAssembly
boolean Wise2_extend_error_path_KmerAssembly(KmerAssemblyLink * real,KmerAssemblyLink * error,int current_pos,char proposed_fix,int proposed_fix_position,KmerAssemblyPath * final_real,KmerAssemblyPath * final_error,int max_search,int * current_path,int max_path_enum);
#define extend_error_path_KmerAssembly Wise2_extend_error_path_KmerAssembly


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/

#ifdef _cplusplus
}
#endif

#endif
