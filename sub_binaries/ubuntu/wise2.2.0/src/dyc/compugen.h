
#ifndef COMPUGEN_HEADER
#define COMPUGEN_HEADER

#include "dyna2.h"

#define CUGEN_MAX_ONEMODEL_TRANS 1024
#define CUGEN_MAX_STATE_NAMES    128
#define CUGEN_MAX_FUNCS          128
#define CUGEN_MAX_PRF_NUM_LEN    20
#define CUGEN_MININF             -1000000

typedef struct OneModelTrans {
  ExprTree * tprf;
  ExprTree * wprf;
  char * tseq;
  ExprTree * wseq_node_to_replace;
  char *wseq;
  int tpos_tprf;
  int tpos_wprf;
  int range;
  int dseq;
  int dwx;
  int dwy;
  int dprf;
  char * from_state;
  char * to_state;
  char * name;
  int    dx;
  int    dy;
  char * xlabel;
  char * ylabel;
  char * map_func;
} OneModelTrans;


typedef struct OneModel {
  char *name;
  char * datatype_str;
  char * statenames    [CUGEN_MAX_STATE_NAMES];
  char * semistatenames[CUGEN_MAX_STATE_NAMES];
  char * midstatenames [CUGEN_MAX_STATE_NAMES];
  char * seqfuncs      [CUGEN_MAX_FUNCS];
  OneModelTrans * trans[CUGEN_MAX_ONEMODEL_TRANS];
  Scope * sc;
  char  * arg_str;
  char  * chain_str;
  int     states_num;
  int     semistates_num;
  int     prof_line_len;
} OneModel;

typedef enum {
  CALC_PRF,
  CALC_SEQ,
  CALC_W,
  MAX_CALC_TYPES
}calc_type_t;

typedef enum {
  COMPUGEN_OK,
  COMPUGEN_TOO_MANY_CALC_FUNCS,
  COMPUGEN_ILLEGAL_SEQ_EXPR,
  COMPUGEN_NO_WSEQ
} compugen_rc_t;

typedef struct {
  char * dy_func_name;
  char * om_func_name;
  int    range;
  char * map_func;
} seq_func_t;

typedef struct {
  int          funcs_num;
  seq_func_t * func;
} seq_func_table_t;

enum {
  OM_TMPSTART,
  OM_TMPEND,
  OM_TMP_STATES_NUM
};

#define is_sum_tree(tree)								\
  ((tree->type == ETR_EXPRESSION) && (tree->nochild == 3) &&				\
   (tree->child[1]->type == ETR_OPERATOR) && (strcmp(tree->child[1]->word, "+") == 0))

OneModel * OneModel_from_GenericMatrix(GenericMatrix * gm);

boolean write_cugen_funcs(OneModel * om,MethodTypeSet * mts,DYNFILE * dfp);
#endif






