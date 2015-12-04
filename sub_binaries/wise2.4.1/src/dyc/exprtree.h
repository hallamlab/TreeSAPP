#ifndef DYNAMITEexprtreeHEADERFILE
#define DYNAMITEexprtreeHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "wisebase.h"

/*
 *

 Expr Tree is about the yacc parse tree. Each node in the
 tree is going to be ExprTree *. The children of each node
 are in the child array, the number of children are in the
 nochild field.

 Each node has a type, and this is the critical layout of
 the node. Each type will layout its children in a particular
 way:

 ETR_NUMBER - has no children. The number is word
 ETR_OPERATOR - has no children. The operator character is in word
 ETR_EXPRESSION - has potentially any number of children, but
    the yacc parser restricts it to (possible_tag) (operator) (possible_tag).
    Don't rely on this though. Generally allowed types in possible_tag
    are NUMBEr, TAG, METHOD or EXPRESSION

 ETR_STATEMENT - is the top level tag. It only has one child which must
     be an ETR_EXPRESSION
 
 ETR_NAME  - is for absolute variable names (or method names). It has no
     children. The actual name is in the word

 ETR_ARRAY - has 2 children. The first must be a "TAG" type and is what
     is indexed, the second must be an EXPRESSION and is what indexes it.

 ETR_TAG -   has any number of children to build up a TAG from NAMES, 
     ->,. (struct refs) or * (REFERNECES).

 ETR_STRUCTREF -  -> or . constructions. It has 3 children: 1st is the
     TAG to the left of the STRUCTREF, second holds either . or -> in word
     and third is the tag to the right.

 ETR_REFERENCE - * or & constructions. It has two children. 1st is * or & in 
     word, the second is the tag which is it is acting on.

 ETR_METHOD - function calls. It has one or two children. The first is a tag (hence
     could be a pointer to a function, or similar). If it is a void function
     it has no other children the second is a commalist.

 ETR_COMMALIST - can only be found in methods. Any number of children, each being
     an argument of the method

*/

    
 
	

	
enum types {
	ETR_NUMBER = 0,
	ETR_OPERATOR,
	ETR_EXPRESSION,
	ETR_STATEMENT,
	ETR_ARRAY,
	ETR_NAME,
	ETR_REFERENCE,
	ETR_STRUCTREF,
	ETR_METHOD,
	ETR_COMMALIST,
	ETR_DECL_VARIABLE,
	ETR_DECL_METHOD,
	ETR_DECL_LIST,
	ETR_TAG
	};

#define EXPRTREE_MAXCHILD 128

#define IS_DECLARED 2
#define IS_TOPLEVEL 4


struct ExprTree {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    struct ExprTree * child[EXPRTREE_MAXCHILD];  
    int nochild;     
    char token;  
    char * word;     
    int type;    
    int attrib;  
    struct ExprTree * parent;    
    int position_in_parent;  
    } ;  
/* ExprTree defined */ 
#ifndef DYNAMITE_DEFINED_ExprTree
typedef struct ExprTree ExprTree;
#define DYNAMITE_DEFINED_ExprTree
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  hard_link_ExprTree(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [ExprTree *]
 *
 * Return [UNKN ]  Undocumented return value [ExprTree *]
 *
 */
ExprTree * hard_link_ExprTree(ExprTree * obj);


/* Function:  ExprTree_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ExprTree *]
 *
 */
ExprTree * ExprTree_alloc(void);


/* Function:  free_ExprTree(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [ExprTree *]
 *
 * Return [UNKN ]  Undocumented return value [ExprTree *]
 *
 */
ExprTree * free_ExprTree(ExprTree * obj);


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
void parentfy_ExprTree(ExprTree * et);
void declared_ExprTree(ExprTree * et);
void find_toplevel_name(ExprTree * et);
void strcat_ExprTree(ExprTree * ExprTree,char * buffer);
void print_ExprTree(ExprTree * ExprTree);
ExprTree * new_ExprTree_decl_method(ExprTree * name,ExprTree * list);
ExprTree * new_ExprTree_decl_variable(ExprTree * type,ExprTree * name);
ExprTree * add_to_decl_list_ExprTree(ExprTree * list,ExprTree * add);
ExprTree * new_ExprTree_decl_list(ExprTree * start);
ExprTree * new_ExprTree_struct_ref(ExprTree * left, ExprTree * ref,ExprTree * right);
ExprTree * new_ExprTree_ref(char op,ExprTree * right);
ExprTree * add_to_commalist_ExprTree(ExprTree * list,ExprTree * add);
ExprTree * new_ExprTree_commalist(ExprTree * start);
ExprTree * new_ExprTree_method(ExprTree * one,ExprTree * other);
ExprTree * new_ExprTree_tag_from_name(ExprTree * name);
ExprTree * new_ExprTree_array(ExprTree * tag,ExprTree * expr) ;
ExprTree * new_ExprTree_binary_expr(ExprTree * left,char op,ExprTree * rgt);
ExprTree * new_ExprTree_token(char t);
boolean add_ExprTree(ExprTree * one,ExprTree * child);
ExprTree * new_ExprTree(void);


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/

#ifdef _cplusplus
}
#endif

#endif
