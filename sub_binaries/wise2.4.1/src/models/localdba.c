#ifdef _cplusplus
extern "C" {
#endif
#include "localdba.h"

# line 5 "localdba.c"


  /*****************   C functions  ****************/
  /*             Written using dynamite            */
  /*            Sat Sep  8 09:05:32 2007           */
  /*            email birney@sanger.ac.uk          */
  /* http://www.sanger.ac.uk/Users/birney/dynamite */
  /*************************************************/


  /* Please report any problems or bugs to         */
  /* Ewan Birney, birney@sanger.ac.uk              */


/* basic set of macros to map states to numbers */ 
#define MATCH65 0    
#define MATCH75 1    
#define MATCH85 2    
#define MATCH95 3    
#define UNMATCHED_QUERY 4    
#define UNMATCHED_TARGET 5   


#define START 0  
#define END 1    


#define LocalDnaMatchBlock_EXPL_MATRIX(this_matrix,i,j,STATE) this_matrix->basematrix->matrix[((j+1)*6)+STATE][i+1]  
#define LocalDnaMatchBlock_EXPL_SPECIAL(matrix,i,j,STATE) matrix->basematrix->specmatrix[STATE][j+1] 
#define LocalDnaMatchBlock_READ_OFF_ERROR -3
    


#define LocalDnaMatchBlock_VSMALL_MATRIX(mat,i,j,STATE) mat->basematrix->matrix[(j+2)%2][((i+1)*6)+STATE]    
#define LocalDnaMatchBlock_VSMALL_SPECIAL(mat,i,j,STATE) mat->basematrix->specmatrix[(j+2)%2][STATE] 




#define LocalDnaMatchBlock_SHATTER_SPECIAL(matrix,i,j,STATE) matrix->shatter->special[STATE][j]  
#define LocalDnaMatchBlock_SHATTER_MATRIX(matrix,i,j,STATE)  fetch_cell_value_ShatterMatrix(mat->shatter,i,j,STATE)  


/* Function:  PackAln_read_Shatter_LocalDnaMatchBlock(mat)
 *
 * Descrip:    Reads off PackAln from shatter matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [LocalDnaMatchBlock *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * PackAln_read_Shatter_LocalDnaMatchBlock(LocalDnaMatchBlock * mat) 
{
    LocalDnaMatchBlock_access_func_holder holder;    


    holder.access_main    = LocalDnaMatchBlock_shatter_access_main;  
    holder.access_special = LocalDnaMatchBlock_shatter_access_special;   
    assert(mat);     
    assert(mat->shatter);    
    return PackAln_read_generic_LocalDnaMatchBlock(mat,holder);  
}    


/* Function:  LocalDnaMatchBlock_shatter_access_main(mat,i,j,state)
 *
 * Descrip: No Description
 *
 * Arg:          mat [UNKN ] Undocumented argument [LocalDnaMatchBlock *]
 * Arg:            i [UNKN ] Undocumented argument [int]
 * Arg:            j [UNKN ] Undocumented argument [int]
 * Arg:        state [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int LocalDnaMatchBlock_shatter_access_main(LocalDnaMatchBlock * mat,int i,int j,int state) 
{
    return LocalDnaMatchBlock_SHATTER_MATRIX(mat,i,j,state); 
}    


/* Function:  LocalDnaMatchBlock_shatter_access_special(mat,i,j,state)
 *
 * Descrip: No Description
 *
 * Arg:          mat [UNKN ] Undocumented argument [LocalDnaMatchBlock *]
 * Arg:            i [UNKN ] Undocumented argument [int]
 * Arg:            j [UNKN ] Undocumented argument [int]
 * Arg:        state [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int LocalDnaMatchBlock_shatter_access_special(LocalDnaMatchBlock * mat,int i,int j,int state) 
{
    return LocalDnaMatchBlock_SHATTER_SPECIAL(mat,i,j,state);    
}    


/* Function:  calculate_shatter_LocalDnaMatchBlock(mat,dpenv)
 *
 * Descrip:    This function calculates the LocalDnaMatchBlock matrix when in shatter mode
 *
 *
 * Arg:          mat [UNKN ] (null) [LocalDnaMatchBlock *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean calculate_shatter_LocalDnaMatchBlock(LocalDnaMatchBlock * mat,DPEnvelope * dpenv) 
{
    int i;   
    int j;   
    int k;   
    int should_calc;     
    int leni;    
    int lenj;    
    int tot; 
    int num; 
    int starti;  
    int startj;  
    int endi;    
    int endj;    


    int * SIG_0_0;   
    int * SIG_1_1;   
    int * SIG_0_1;   
    int * SIG_1_0;   


    leni = mat->leni;    
    lenj = mat->lenj;    


    mat->shatter = new_ShatterMatrix(dpenv,6,lenj,2);    
    prepare_DPEnvelope(dpenv);   
    starti = dpenv->starti;  
    if( starti < 0 ) 
      starti = 0;    
    startj = dpenv->startj;  
    if( startj < 0 ) 
      startj = 0;    
    endi = dpenv->endi;  
    if( endi > mat->leni )   
      endi = mat->leni;  
    endj = dpenv->endj;  
    if( endj > mat->lenj )   
      endj = mat->lenj;  
    tot = (endi-starti) * (endj-startj); 
    num = 0; 


    start_reporting("LocalDnaMatchBlock Matrix calculation: ");  
    for(j=startj;j<endj;j++) {  
      auto int score;    
      auto int temp;     
      for(i=starti;i<endi;i++)   {  
        /* Check if is in envelope - code identical to is_in_DPEnvelope, but aggressively inlined here for speed */ 
        should_calc = 0; 
        for(k=0;k<dpenv->len;k++)    {  
          auto DPUnit * u;   
          u = dpenv->dpu[k]; 
          switch(u->type)    {  
            case DPENV_RECT :    
              if( i >= u->starti && j >= u->startj && i < (u->starti+u->height) && j < (u->startj+u->length))    
                should_calc = 1;     
              break; 
            case DPENV_DIAG :    
              if(  abs( (i-j) - (u->starti-u->startj)) <= u->height && i+j >= u->starti+u->startj && i+j+u->length >= u->starti+u->startj)   
                should_calc = 1;     
              break; 
            }  
          if( should_calc == 1 ) 
            break;   
          }  
        if( should_calc == 0)    
          continue;  


        SIG_0_0 = fetch_cell_from_ShatterMatrix(mat->shatter,i,j);   
        SIG_1_1 = fetch_cell_from_ShatterMatrix(mat->shatter,i-1,j-1);   
        SIG_0_1 = fetch_cell_from_ShatterMatrix(mat->shatter,i-0,j-1);   
        SIG_1_0 = fetch_cell_from_ShatterMatrix(mat->shatter,i-1,j-0);   




        /* For state MATCH65 */ 
        /* setting first movement to score */ 
        score = SIG_1_1[MATCH65] + (mat->lchs->comp65->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->s);   
        /* From state MATCH65 to state MATCH65 */ 
        temp = SIG_0_1[MATCH65] + mat->lchs->g;  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH65 to state MATCH65 */ 
        temp = SIG_1_0[MATCH65] + mat->lchs->g;  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state START to state MATCH65 */ 
        temp = LocalDnaMatchBlock_SHATTER_SPECIAL(mat,i-1,j-1,START) + (mat->lchs->comp65->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->v);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state UNMATCHED_TARGET to state MATCH65 */ 
        temp = SIG_1_1[UNMATCHED_TARGET] + (mat->lchs->comp65->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->v);   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for MATCH65 */ 
        /* Add any movement independant score and put away */ 
         SIG_0_0[MATCH65] = score;   


        /* state MATCH65 is a source for special END */ 
        temp = score + (0) + (0) ;   
        if( temp > LocalDnaMatchBlock_SHATTER_SPECIAL(mat,i,j,END) )     {  
          LocalDnaMatchBlock_SHATTER_SPECIAL(mat,i,j,END) = temp;    
          }  




        /* Finished calculating state MATCH65 */ 


        /* For state MATCH75 */ 
        /* setting first movement to score */ 
        score = SIG_1_1[MATCH75] + (mat->lchs->comp75->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->s);   
        /* From state MATCH75 to state MATCH75 */ 
        temp = SIG_0_1[MATCH75] + mat->lchs->g;  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH75 to state MATCH75 */ 
        temp = SIG_1_0[MATCH75] + mat->lchs->g;  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state START to state MATCH75 */ 
        temp = LocalDnaMatchBlock_SHATTER_SPECIAL(mat,i-1,j-1,START) + (mat->lchs->comp75->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->v);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state UNMATCHED_TARGET to state MATCH75 */ 
        temp = SIG_1_1[UNMATCHED_TARGET] + (mat->lchs->comp75->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->v);   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for MATCH75 */ 
        /* Add any movement independant score and put away */ 
         SIG_0_0[MATCH75] = score;   


        /* state MATCH75 is a source for special END */ 
        temp = score + (0) + (0) ;   
        if( temp > LocalDnaMatchBlock_SHATTER_SPECIAL(mat,i,j,END) )     {  
          LocalDnaMatchBlock_SHATTER_SPECIAL(mat,i,j,END) = temp;    
          }  




        /* Finished calculating state MATCH75 */ 


        /* For state MATCH85 */ 
        /* setting first movement to score */ 
        score = SIG_1_1[MATCH85] + (mat->lchs->comp85->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->s);   
        /* From state MATCH85 to state MATCH85 */ 
        temp = SIG_0_1[MATCH85] + mat->lchs->g;  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH85 to state MATCH85 */ 
        temp = SIG_1_0[MATCH85] + mat->lchs->g;  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state START to state MATCH85 */ 
        temp = LocalDnaMatchBlock_SHATTER_SPECIAL(mat,i-1,j-1,START) + (mat->lchs->comp85->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->v);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state UNMATCHED_TARGET to state MATCH85 */ 
        temp = SIG_1_1[UNMATCHED_TARGET] + (mat->lchs->comp85->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->v);   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for MATCH85 */ 
        /* Add any movement independant score and put away */ 
         SIG_0_0[MATCH85] = score;   


        /* state MATCH85 is a source for special END */ 
        temp = score + (0) + (0) ;   
        if( temp > LocalDnaMatchBlock_SHATTER_SPECIAL(mat,i,j,END) )     {  
          LocalDnaMatchBlock_SHATTER_SPECIAL(mat,i,j,END) = temp;    
          }  




        /* Finished calculating state MATCH85 */ 


        /* For state MATCH95 */ 
        /* setting first movement to score */ 
        score = SIG_1_1[MATCH95] + (mat->lchs->comp95->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->s);   
        /* From state MATCH95 to state MATCH95 */ 
        temp = SIG_0_1[MATCH95] + mat->lchs->g;  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH95 to state MATCH95 */ 
        temp = SIG_1_0[MATCH95] + mat->lchs->g;  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state START to state MATCH95 */ 
        temp = LocalDnaMatchBlock_SHATTER_SPECIAL(mat,i-1,j-1,START) + (mat->lchs->comp95->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->v);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state UNMATCHED_TARGET to state MATCH95 */ 
        temp = SIG_1_1[UNMATCHED_TARGET] + (mat->lchs->comp95->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->v);   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for MATCH95 */ 
        /* Add any movement independant score and put away */ 
         SIG_0_0[MATCH95] = score;   


        /* state MATCH95 is a source for special END */ 
        temp = score + (0) + (0) ;   
        if( temp > LocalDnaMatchBlock_SHATTER_SPECIAL(mat,i,j,END) )     {  
          LocalDnaMatchBlock_SHATTER_SPECIAL(mat,i,j,END) = temp;    
          }  




        /* Finished calculating state MATCH95 */ 


        /* For state UNMATCHED_QUERY */ 
        /* setting first movement to score */ 
        score = SIG_1_0[MATCH65] + mat->lchs->b;     
        /* From state MATCH75 to state UNMATCHED_QUERY */ 
        temp = SIG_1_0[MATCH75] + mat->lchs->b;  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH85 to state UNMATCHED_QUERY */ 
        temp = SIG_1_0[MATCH85] + mat->lchs->b;  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH95 to state UNMATCHED_QUERY */ 
        temp = SIG_1_0[MATCH95] + mat->lchs->b;  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state UNMATCHED_QUERY to state UNMATCHED_QUERY */ 
        temp = SIG_1_0[UNMATCHED_QUERY] + mat->lchs->u;  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for UNMATCHED_QUERY */ 
        /* Add any movement independant score and put away */ 
         SIG_0_0[UNMATCHED_QUERY] = score;   


        /* Finished calculating state UNMATCHED_QUERY */ 


        /* For state UNMATCHED_TARGET */ 
        /* setting first movement to score */ 
        score = SIG_0_1[UNMATCHED_QUERY] + mat->lchs->v;     
        /* From state UNMATCHED_TARGET to state UNMATCHED_TARGET */ 
        temp = SIG_0_1[UNMATCHED_TARGET] + mat->lchs->u;     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for UNMATCHED_TARGET */ 
        /* Add any movement independant score and put away */ 
         SIG_0_0[UNMATCHED_TARGET] = score;  


        /* Finished calculating state UNMATCHED_TARGET */ 
        }  


      /* Special state START has no special to special movements */ 


      /* Special state END has no special to special movements */ 
      }  
    stop_reporting();    
    return TRUE;     
}    


/* Function:  search_LocalDnaMatchBlock(dbsi,out,query,target,lchs)
 *
 * Descrip:    This function makes a database search of LocalDnaMatchBlock
 *             It uses the dbsi structure to choose which implementation
 *             to use of the database searching. This way at run time you
 *             can switch between single threaded/multi-threaded or hardware
 *
 *
 * Arg:          dbsi [UNKN ] Undocumented argument [DBSearchImpl *]
 * Arg:           out [UNKN ] Undocumented argument [Hscore *]
 * Arg:         query [UNKN ] Undocumented argument [ComplexSequence*]
 * Arg:        target [UNKN ] Undocumented argument [ComplexSequence*]
 * Arg:          lchs [UNKN ] Undocumented argument [LocalCisHitScore*]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type search_LocalDnaMatchBlock(DBSearchImpl * dbsi,Hscore * out,ComplexSequence* query,ComplexSequence* target ,LocalCisHitScore* lchs) 
{
    if( out == NULL )    {  
      warn("Passed in a null Hscore object into search_LocalDnaMatchBlock. Can't process results!"); 
      return SEARCH_ERROR;   
      }  
    if( dbsi == NULL )   {  
      warn("Passed in a null DBSearchImpl object into search_LocalDnaMatchBlock. Can't process results!");   
      return SEARCH_ERROR;   
      }  
    if( dbsi->trace_level > 5 )  
      warn("Asking for trace level of %d in database search for LocalDnaMatchBlock, but it was compiled with a trace level of -2139062144. Not all trace statements can be shown",dbsi->trace_level);    
    switch(dbsi->type)   { /*switch on implementation*/ 
      case DBSearchImpl_Serial : 
        return serial_search_LocalDnaMatchBlock(out,query,target ,lchs); 
      case DBSearchImpl_Pthreads :   
        warn("This matrix LocalDnaMatchBlock was not dyc compiled with thread support"); 
        return SEARCH_ERROR; 
      default :  
        warn("database search implementation %s was not provided in the compiled dynamite file from LocalDnaMatchBlock",impl_string_DBSearchImpl(dbsi)); 
        return SEARCH_ERROR; 
      } /* end of switch on implementation */ 


}    


/* Function:  serial_search_LocalDnaMatchBlock(out,query,target,lchs)
 *
 * Descrip:    This function makes a database search of LocalDnaMatchBlock
 *             It is a single processor implementation
 *
 *
 * Arg:           out [UNKN ] Undocumented argument [Hscore *]
 * Arg:         query [UNKN ] Undocumented argument [ComplexSequence*]
 * Arg:        target [UNKN ] Undocumented argument [ComplexSequence*]
 * Arg:          lchs [UNKN ] Undocumented argument [LocalCisHitScore*]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type serial_search_LocalDnaMatchBlock(Hscore * out,ComplexSequence* query,ComplexSequence* target ,LocalCisHitScore* lchs) 
{
    int db_status;   
    int score;   
    int query_pos = 0;   
    int target_pos = 0;  
    DataScore * ds;  


    push_errormsg_stack("Before any actual search in db searching"); 


    target_pos = 0;  




    /* No maximum length - allocated on-the-fly */ 
    score = score_only_LocalDnaMatchBlock(query, target , lchs);     
    if( should_store_Hscore(out,score) == TRUE )     { /*if storing datascore*/ 
      ds = new_DataScore_from_storage(out);  
      if( ds == NULL )   {  
        warn("LocalDnaMatchBlock search had a memory error in allocating a new_DataScore (?a leak somewhere - DataScore is a very small datastructure"); 
        return SEARCH_ERROR; 
        }  
      /* Now: add query/target information to the entry */ 
      ds->score = score;     
      add_Hscore(out,ds);    
      } /* end of if storing datascore */ 
    pop_errormsg_stack();    
    push_errormsg_stack("DB searching: just finished [Query Pos: %d] [Target Pos: %d]",query_pos,target_pos);    


    pop_errormsg_stack();    
    return SEARCH_OK;    
}    


/* Function:  score_only_LocalDnaMatchBlock(query,target,lchs)
 *
 * Descrip:    This function just calculates the score for the matrix
 *             I am pretty sure we can do this better, but hey, for the moment...
 *             It calls /allocate_LocalDnaMatchBlock_only
 *
 *
 * Arg:         query [UNKN ] query data structure [ComplexSequence*]
 * Arg:        target [UNKN ] target data structure [ComplexSequence*]
 * Arg:          lchs [UNKN ] Resource [LocalCisHitScore*]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int score_only_LocalDnaMatchBlock(ComplexSequence* query,ComplexSequence* target ,LocalCisHitScore* lchs) 
{
    int bestscore = NEGI;    
    int i;   
    int j;   
    int k;   
    LocalDnaMatchBlock * mat;    


    mat = allocate_LocalDnaMatchBlock_only(query, target , lchs);    
    if( mat == NULL )    {  
      warn("Memory allocation error in the db search - unable to communicate to calling function. this spells DIASTER!");    
      return NEGI;   
      }  
    if((mat->basematrix = BaseMatrix_alloc_matrix_and_specials(2,(mat->leni + 1) * 6,2,2)) == NULL)  {  
      warn("Score only matrix for LocalDnaMatchBlock cannot be allocated, (asking for 1  by %d  cells)",mat->leni*6);    
      mat = free_LocalDnaMatchBlock(mat);    
      return 0;  
      }  
    mat->basematrix->type = BASEMATRIX_TYPE_VERYSMALL;   


    /* Now, initiate matrix */ 
    for(j=0;j<3;j++) {  
      for(i=(-1);i<mat->leni;i++)    {  
        for(k=0;k<6;k++) 
          LocalDnaMatchBlock_VSMALL_MATRIX(mat,i,j,k) = NEGI;    
        }  
      LocalDnaMatchBlock_VSMALL_SPECIAL(mat,i,j,START) = 0;  
      LocalDnaMatchBlock_VSMALL_SPECIAL(mat,i,j,END) = NEGI; 
      }  


    /* Ok, lets do-o-o-o-o it */ 


    for(j=0;j<mat->lenj;j++) { /*for all target positions*/ 
      auto int score;    
      auto int temp;     
      for(i=0;i<mat->leni;i++)   { /*for all query positions*/ 


        /* For state MATCH65 */ 
        /* setting first movement to score */ 
        score = LocalDnaMatchBlock_VSMALL_MATRIX(mat,i-1,j-1,MATCH65) + (mat->lchs->comp65->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->s);  
        /* From state MATCH65 to state MATCH65 */ 
        temp = LocalDnaMatchBlock_VSMALL_MATRIX(mat,i-0,j-1,MATCH65) + mat->lchs->g;     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH65 to state MATCH65 */ 
        temp = LocalDnaMatchBlock_VSMALL_MATRIX(mat,i-1,j-0,MATCH65) + mat->lchs->g;     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state START to state MATCH65 */ 
        temp = LocalDnaMatchBlock_VSMALL_SPECIAL(mat,i-1,j-1,START) + (mat->lchs->comp65->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->v);    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state UNMATCHED_TARGET to state MATCH65 */ 
        temp = LocalDnaMatchBlock_VSMALL_MATRIX(mat,i-1,j-1,UNMATCHED_TARGET) + (mat->lchs->comp65->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->v);  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for MATCH65 */ 
        /* Add any movement independant score and put away */ 
         LocalDnaMatchBlock_VSMALL_MATRIX(mat,i,j,MATCH65) = score;  


        /* state MATCH65 is a source for special END */ 
        temp = score + (0) + (0) ;   
        if( temp > LocalDnaMatchBlock_VSMALL_SPECIAL(mat,i,j,END) )  {  
          LocalDnaMatchBlock_VSMALL_SPECIAL(mat,i,j,END) = temp;     
          }  




        /* Finished calculating state MATCH65 */ 


        /* For state MATCH75 */ 
        /* setting first movement to score */ 
        score = LocalDnaMatchBlock_VSMALL_MATRIX(mat,i-1,j-1,MATCH75) + (mat->lchs->comp75->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->s);  
        /* From state MATCH75 to state MATCH75 */ 
        temp = LocalDnaMatchBlock_VSMALL_MATRIX(mat,i-0,j-1,MATCH75) + mat->lchs->g;     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH75 to state MATCH75 */ 
        temp = LocalDnaMatchBlock_VSMALL_MATRIX(mat,i-1,j-0,MATCH75) + mat->lchs->g;     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state START to state MATCH75 */ 
        temp = LocalDnaMatchBlock_VSMALL_SPECIAL(mat,i-1,j-1,START) + (mat->lchs->comp75->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->v);    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state UNMATCHED_TARGET to state MATCH75 */ 
        temp = LocalDnaMatchBlock_VSMALL_MATRIX(mat,i-1,j-1,UNMATCHED_TARGET) + (mat->lchs->comp75->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->v);  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for MATCH75 */ 
        /* Add any movement independant score and put away */ 
         LocalDnaMatchBlock_VSMALL_MATRIX(mat,i,j,MATCH75) = score;  


        /* state MATCH75 is a source for special END */ 
        temp = score + (0) + (0) ;   
        if( temp > LocalDnaMatchBlock_VSMALL_SPECIAL(mat,i,j,END) )  {  
          LocalDnaMatchBlock_VSMALL_SPECIAL(mat,i,j,END) = temp;     
          }  




        /* Finished calculating state MATCH75 */ 


        /* For state MATCH85 */ 
        /* setting first movement to score */ 
        score = LocalDnaMatchBlock_VSMALL_MATRIX(mat,i-1,j-1,MATCH85) + (mat->lchs->comp85->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->s);  
        /* From state MATCH85 to state MATCH85 */ 
        temp = LocalDnaMatchBlock_VSMALL_MATRIX(mat,i-0,j-1,MATCH85) + mat->lchs->g;     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH85 to state MATCH85 */ 
        temp = LocalDnaMatchBlock_VSMALL_MATRIX(mat,i-1,j-0,MATCH85) + mat->lchs->g;     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state START to state MATCH85 */ 
        temp = LocalDnaMatchBlock_VSMALL_SPECIAL(mat,i-1,j-1,START) + (mat->lchs->comp85->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->v);    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state UNMATCHED_TARGET to state MATCH85 */ 
        temp = LocalDnaMatchBlock_VSMALL_MATRIX(mat,i-1,j-1,UNMATCHED_TARGET) + (mat->lchs->comp85->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->v);  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for MATCH85 */ 
        /* Add any movement independant score and put away */ 
         LocalDnaMatchBlock_VSMALL_MATRIX(mat,i,j,MATCH85) = score;  


        /* state MATCH85 is a source for special END */ 
        temp = score + (0) + (0) ;   
        if( temp > LocalDnaMatchBlock_VSMALL_SPECIAL(mat,i,j,END) )  {  
          LocalDnaMatchBlock_VSMALL_SPECIAL(mat,i,j,END) = temp;     
          }  




        /* Finished calculating state MATCH85 */ 


        /* For state MATCH95 */ 
        /* setting first movement to score */ 
        score = LocalDnaMatchBlock_VSMALL_MATRIX(mat,i-1,j-1,MATCH95) + (mat->lchs->comp95->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->s);  
        /* From state MATCH95 to state MATCH95 */ 
        temp = LocalDnaMatchBlock_VSMALL_MATRIX(mat,i-0,j-1,MATCH95) + mat->lchs->g;     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH95 to state MATCH95 */ 
        temp = LocalDnaMatchBlock_VSMALL_MATRIX(mat,i-1,j-0,MATCH95) + mat->lchs->g;     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state START to state MATCH95 */ 
        temp = LocalDnaMatchBlock_VSMALL_SPECIAL(mat,i-1,j-1,START) + (mat->lchs->comp95->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->v);    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state UNMATCHED_TARGET to state MATCH95 */ 
        temp = LocalDnaMatchBlock_VSMALL_MATRIX(mat,i-1,j-1,UNMATCHED_TARGET) + (mat->lchs->comp95->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->v);  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for MATCH95 */ 
        /* Add any movement independant score and put away */ 
         LocalDnaMatchBlock_VSMALL_MATRIX(mat,i,j,MATCH95) = score;  


        /* state MATCH95 is a source for special END */ 
        temp = score + (0) + (0) ;   
        if( temp > LocalDnaMatchBlock_VSMALL_SPECIAL(mat,i,j,END) )  {  
          LocalDnaMatchBlock_VSMALL_SPECIAL(mat,i,j,END) = temp;     
          }  




        /* Finished calculating state MATCH95 */ 


        /* For state UNMATCHED_QUERY */ 
        /* setting first movement to score */ 
        score = LocalDnaMatchBlock_VSMALL_MATRIX(mat,i-1,j-0,MATCH65) + mat->lchs->b;    
        /* From state MATCH75 to state UNMATCHED_QUERY */ 
        temp = LocalDnaMatchBlock_VSMALL_MATRIX(mat,i-1,j-0,MATCH75) + mat->lchs->b;     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH85 to state UNMATCHED_QUERY */ 
        temp = LocalDnaMatchBlock_VSMALL_MATRIX(mat,i-1,j-0,MATCH85) + mat->lchs->b;     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH95 to state UNMATCHED_QUERY */ 
        temp = LocalDnaMatchBlock_VSMALL_MATRIX(mat,i-1,j-0,MATCH95) + mat->lchs->b;     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state UNMATCHED_QUERY to state UNMATCHED_QUERY */ 
        temp = LocalDnaMatchBlock_VSMALL_MATRIX(mat,i-1,j-0,UNMATCHED_QUERY) + mat->lchs->u;     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for UNMATCHED_QUERY */ 
        /* Add any movement independant score and put away */ 
         LocalDnaMatchBlock_VSMALL_MATRIX(mat,i,j,UNMATCHED_QUERY) = score;  


        /* Finished calculating state UNMATCHED_QUERY */ 


        /* For state UNMATCHED_TARGET */ 
        /* setting first movement to score */ 
        score = LocalDnaMatchBlock_VSMALL_MATRIX(mat,i-0,j-1,UNMATCHED_QUERY) + mat->lchs->v;    
        /* From state UNMATCHED_TARGET to state UNMATCHED_TARGET */ 
        temp = LocalDnaMatchBlock_VSMALL_MATRIX(mat,i-0,j-1,UNMATCHED_TARGET) + mat->lchs->u;    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for UNMATCHED_TARGET */ 
        /* Add any movement independant score and put away */ 
         LocalDnaMatchBlock_VSMALL_MATRIX(mat,i,j,UNMATCHED_TARGET) = score; 


        /* Finished calculating state UNMATCHED_TARGET */ 
        } /* end of for all query positions */ 




      /* Special state START has no special to special movements */ 


      /* Special state END has no special to special movements */ 
      if( bestscore < LocalDnaMatchBlock_VSMALL_SPECIAL(mat,0,j,END) )   
        bestscore = LocalDnaMatchBlock_VSMALL_SPECIAL(mat,0,j,END);  
      } /* end of for all target positions */ 


    mat = free_LocalDnaMatchBlock(mat);  
    return bestscore;    
}    


/* Function:  PackAln_bestmemory_LocalDnaMatchBlock(query,target,lchs,dpenv,dpri)
 *
 * Descrip:    This function chooses the best memory set-up for the alignment
 *             using calls to basematrix, and then implements either a large
 *             or small memory model.
 *
 *             It is the best function to use if you just want an alignment
 *
 *             If you want a label alignment, you will need
 *             /convert_PackAln_to_AlnBlock_LocalDnaMatchBlock
 *
 *
 * Arg:         query [UNKN ] query data structure [ComplexSequence*]
 * Arg:        target [UNKN ] target data structure [ComplexSequence*]
 * Arg:          lchs [UNKN ] Resource [LocalCisHitScore*]
 * Arg:         dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:          dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * PackAln_bestmemory_LocalDnaMatchBlock(ComplexSequence* query,ComplexSequence* target ,LocalCisHitScore* lchs,DPEnvelope * dpenv,DPRunImpl * dpri) 
{
    long long total; 
    LocalDnaMatchBlock * mat;    
    PackAln * out;   
    DebugMatrix * de;    
    DPRunImplMemory strategy;    
    assert(dpri);    


    total = query->seq->len * target->seq->len;  
    if( dpri->memory == DPIM_Default )   {  
      if( (total * 6 * sizeof(int)) > 1000*dpri->kbyte_size) {  
        strategy = DPIM_Linear;  
        }  
      else   {  
        strategy = DPIM_Explicit;    
        }  
      }  
    else {  
      strategy = dpri->memory;   
      }  


    if( dpenv != NULL )  {  
      if( strategy == DPIM_Explicit) {  
        if( (mat=allocate_Expl_LocalDnaMatchBlock(query, target , lchs,dpri)) == NULL )  {  
          warn("Unable to allocate large LocalDnaMatchBlock version");   
          return NULL;   
          }  
        calculate_dpenv_LocalDnaMatchBlock(mat,dpenv);   
        out =  PackAln_read_Expl_LocalDnaMatchBlock(mat);    
        }  
      else   {  
        mat = allocate_LocalDnaMatchBlock_only(query, target , lchs);    
        calculate_shatter_LocalDnaMatchBlock(mat,dpenv);     
        out = PackAln_read_Shatter_LocalDnaMatchBlock(mat);  
        }  
      }  
    else {  
      if( strategy == DPIM_Linear )  {  
        /* use small implementation */ 
        if( (mat=allocate_Small_LocalDnaMatchBlock(query, target , lchs)) == NULL )  {  
          warn("Unable to allocate small LocalDnaMatchBlock version");   
          return NULL;   
          }  
        out = PackAln_calculate_Small_LocalDnaMatchBlock(mat,dpenv);     
        }  
      else   {  
        /* use Large implementation */ 
        if( (mat=allocate_Expl_LocalDnaMatchBlock(query, target , lchs,dpri)) == NULL )  {  
          warn("Unable to allocate large LocalDnaMatchBlock version");   
          return NULL;   
          }  
        if( dpri->debug == TRUE) {  
          fatal("Asked for dydebug, but dynamite file not compiled with -g. Need to recompile dynamite source"); 
          }  
        else {  
          calculate_LocalDnaMatchBlock(mat);     
          out =  PackAln_read_Expl_LocalDnaMatchBlock(mat);  
          }  
        }  
      }  


    mat = free_LocalDnaMatchBlock(mat);  
    return out;  
}    


/* Function:  allocate_LocalDnaMatchBlock_only(query,target,lchs)
 *
 * Descrip:    This function only allocates the LocalDnaMatchBlock structure
 *             checks types where possible and determines leni and lenj
 *             The basematrix area is delt with elsewhere
 *
 *
 * Arg:         query [UNKN ] query data structure [ComplexSequence*]
 * Arg:        target [UNKN ] target data structure [ComplexSequence*]
 * Arg:          lchs [UNKN ] Resource [LocalCisHitScore*]
 *
 * Return [UNKN ]  Undocumented return value [LocalDnaMatchBlock *]
 *
 */
LocalDnaMatchBlock * allocate_LocalDnaMatchBlock_only(ComplexSequence* query,ComplexSequence* target ,LocalCisHitScore* lchs) 
{
    LocalDnaMatchBlock * out;    


    if((out= LocalDnaMatchBlock_alloc()) == NULL)    {  
      warn("Allocation of basic LocalDnaMatchBlock structure failed...");    
      return NULL;   
      }  


    out->query = query;  
    out->target = target;    
    out->lchs = lchs;    
    out->leni = query->seq->len;     
    out->lenj = target->seq->len;    
    return out;  
}    


/* Function:  allocate_Expl_LocalDnaMatchBlock(query,target,lchs,dpri)
 *
 * Descrip:    This function allocates the LocalDnaMatchBlock structure
 *             and the basematrix area for explicit memory implementations
 *             It calls /allocate_LocalDnaMatchBlock_only
 *
 *
 * Arg:         query [UNKN ] query data structure [ComplexSequence*]
 * Arg:        target [UNKN ] target data structure [ComplexSequence*]
 * Arg:          lchs [UNKN ] Resource [LocalCisHitScore*]
 * Arg:          dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [LocalDnaMatchBlock *]
 *
 */
LocalDnaMatchBlock * allocate_Expl_LocalDnaMatchBlock(ComplexSequence* query,ComplexSequence* target ,LocalCisHitScore* lchs,DPRunImpl * dpri) 
{
    LocalDnaMatchBlock * out;    


    out = allocate_LocalDnaMatchBlock_only(query, target , lchs);    
    if( out == NULL )    
      return NULL;   
    if( dpri->should_cache == TRUE ) {  
      if( dpri->cache != NULL )  {  
        if( dpri->cache->maxleni >= (out->lenj+1)*6 && dpri->cache->maxlenj >= (out->leni+1))    
          out->basematrix = hard_link_BaseMatrix(dpri->cache);   
        else 
          dpri->cache = free_BaseMatrix(dpri->cache);    
        }  
      }  
    if( out->basematrix == NULL )    {  
      if( (out->basematrix = BaseMatrix_alloc_matrix_and_specials((out->lenj+1)*6,(out->leni+1),2,out->lenj+1)) == NULL) {  
        warn("Explicit matrix LocalDnaMatchBlock cannot be allocated, (asking for %d by %d main cells)",out->leni,out->lenj);    
        free_LocalDnaMatchBlock(out);    
        return NULL; 
        }  
      }  
    if( dpri->should_cache == TRUE && dpri->cache == NULL)   
      dpri->cache = hard_link_BaseMatrix(out->basematrix);   
    out->basematrix->type = BASEMATRIX_TYPE_EXPLICIT;    
    init_LocalDnaMatchBlock(out);    
    return out;  
}    


/* Function:  init_LocalDnaMatchBlock(mat)
 *
 * Descrip:    This function initates LocalDnaMatchBlock matrix when in explicit mode
 *             Called in /allocate_Expl_LocalDnaMatchBlock
 *
 *
 * Arg:        mat [UNKN ] LocalDnaMatchBlock which contains explicit basematrix memory [LocalDnaMatchBlock *]
 *
 */
void init_LocalDnaMatchBlock(LocalDnaMatchBlock * mat) 
{
    register int i;  
    register int j;  
    if( mat->basematrix->type != BASEMATRIX_TYPE_EXPLICIT)   {  
      warn("Cannot iniate matrix, is not an explicit memory type and you have assummed that");   
      return;    
      }  


    for(i= (-1);i<mat->query->seq->len;i++)  {  
      for(j= (-1);j<2;j++)   {  
        LocalDnaMatchBlock_EXPL_MATRIX(mat,i,j,MATCH65) = NEGI;  
        LocalDnaMatchBlock_EXPL_MATRIX(mat,i,j,MATCH75) = NEGI;  
        LocalDnaMatchBlock_EXPL_MATRIX(mat,i,j,MATCH85) = NEGI;  
        LocalDnaMatchBlock_EXPL_MATRIX(mat,i,j,MATCH95) = NEGI;  
        LocalDnaMatchBlock_EXPL_MATRIX(mat,i,j,UNMATCHED_QUERY) = NEGI;  
        LocalDnaMatchBlock_EXPL_MATRIX(mat,i,j,UNMATCHED_TARGET) = NEGI; 
        }  
      }  
    for(j= (-1);j<mat->target->seq->len;j++) {  
      for(i= (-1);i<2;i++)   {  
        LocalDnaMatchBlock_EXPL_MATRIX(mat,i,j,MATCH65) = NEGI;  
        LocalDnaMatchBlock_EXPL_MATRIX(mat,i,j,MATCH75) = NEGI;  
        LocalDnaMatchBlock_EXPL_MATRIX(mat,i,j,MATCH85) = NEGI;  
        LocalDnaMatchBlock_EXPL_MATRIX(mat,i,j,MATCH95) = NEGI;  
        LocalDnaMatchBlock_EXPL_MATRIX(mat,i,j,UNMATCHED_QUERY) = NEGI;  
        LocalDnaMatchBlock_EXPL_MATRIX(mat,i,j,UNMATCHED_TARGET) = NEGI; 
        }  
      LocalDnaMatchBlock_EXPL_SPECIAL(mat,i,j,START) = 0;    
      LocalDnaMatchBlock_EXPL_SPECIAL(mat,i,j,END) = NEGI;   
      }  
    return;  
}    


/* Function:  recalculate_PackAln_LocalDnaMatchBlock(pal,mat)
 *
 * Descrip:    This function recalculates the PackAln structure produced by LocalDnaMatchBlock
 *             For example, in linear space methods this is used to score them
 *
 *
 * Arg:        pal [UNKN ] Undocumented argument [PackAln *]
 * Arg:        mat [UNKN ] Undocumented argument [LocalDnaMatchBlock *]
 *
 */
void recalculate_PackAln_LocalDnaMatchBlock(PackAln * pal,LocalDnaMatchBlock * mat) 
{
    int i,j,k,offi,offj; 
    PackAlnUnit * prev;  
    PackAlnUnit * pau;   


    for(k=1,prev=pal->pau[0];k < pal->len;k++,prev=pau)  {  
      pau = pal->pau[k]; 
      i = pau->i;    
      j = pau->j;    
      offi = pau->i - prev->i;   
      offj = pau->j - prev->j;   
      switch(pau->state) {  
        case MATCH65 :   
          if( offi == 1 && offj == 1 && prev->state == MATCH65 ) {  
            pau->score = (mat->lchs->comp65->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->s) + (0);   
            continue;    
            }  
          if( offi == 0 && offj == 1 && prev->state == MATCH65 ) {  
            pau->score = mat->lchs->g + (0);     
            continue;    
            }  
          if( offi == 1 && offj == 0 && prev->state == MATCH65 ) {  
            pau->score = mat->lchs->g + (0);     
            continue;    
            }  
          if( offj == 1 && prev->state == (START+6) )    {  
            pau->score = (mat->lchs->comp65->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->v) + (0);   
            continue;    
            }  
          if( offi == 1 && offj == 1 && prev->state == UNMATCHED_TARGET )    {  
            pau->score = (mat->lchs->comp65->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->v) + (0);   
            continue;    
            }  
          warn("In recaluclating PackAln with state MATCH65, from [%d,%d,%d], got a bad source state. Error!",offi,offj,prev->state);    
          break; 
        case MATCH75 :   
          if( offi == 1 && offj == 1 && prev->state == MATCH75 ) {  
            pau->score = (mat->lchs->comp75->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->s) + (0);   
            continue;    
            }  
          if( offi == 0 && offj == 1 && prev->state == MATCH75 ) {  
            pau->score = mat->lchs->g + (0);     
            continue;    
            }  
          if( offi == 1 && offj == 0 && prev->state == MATCH75 ) {  
            pau->score = mat->lchs->g + (0);     
            continue;    
            }  
          if( offj == 1 && prev->state == (START+6) )    {  
            pau->score = (mat->lchs->comp75->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->v) + (0);   
            continue;    
            }  
          if( offi == 1 && offj == 1 && prev->state == UNMATCHED_TARGET )    {  
            pau->score = (mat->lchs->comp75->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->v) + (0);   
            continue;    
            }  
          warn("In recaluclating PackAln with state MATCH75, from [%d,%d,%d], got a bad source state. Error!",offi,offj,prev->state);    
          break; 
        case MATCH85 :   
          if( offi == 1 && offj == 1 && prev->state == MATCH85 ) {  
            pau->score = (mat->lchs->comp85->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->s) + (0);   
            continue;    
            }  
          if( offi == 0 && offj == 1 && prev->state == MATCH85 ) {  
            pau->score = mat->lchs->g + (0);     
            continue;    
            }  
          if( offi == 1 && offj == 0 && prev->state == MATCH85 ) {  
            pau->score = mat->lchs->g + (0);     
            continue;    
            }  
          if( offj == 1 && prev->state == (START+6) )    {  
            pau->score = (mat->lchs->comp85->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->v) + (0);   
            continue;    
            }  
          if( offi == 1 && offj == 1 && prev->state == UNMATCHED_TARGET )    {  
            pau->score = (mat->lchs->comp85->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->v) + (0);   
            continue;    
            }  
          warn("In recaluclating PackAln with state MATCH85, from [%d,%d,%d], got a bad source state. Error!",offi,offj,prev->state);    
          break; 
        case MATCH95 :   
          if( offi == 1 && offj == 1 && prev->state == MATCH95 ) {  
            pau->score = (mat->lchs->comp95->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->s) + (0);   
            continue;    
            }  
          if( offi == 0 && offj == 1 && prev->state == MATCH95 ) {  
            pau->score = mat->lchs->g + (0);     
            continue;    
            }  
          if( offi == 1 && offj == 0 && prev->state == MATCH95 ) {  
            pau->score = mat->lchs->g + (0);     
            continue;    
            }  
          if( offj == 1 && prev->state == (START+6) )    {  
            pau->score = (mat->lchs->comp95->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->v) + (0);   
            continue;    
            }  
          if( offi == 1 && offj == 1 && prev->state == UNMATCHED_TARGET )    {  
            pau->score = (mat->lchs->comp95->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->v) + (0);   
            continue;    
            }  
          warn("In recaluclating PackAln with state MATCH95, from [%d,%d,%d], got a bad source state. Error!",offi,offj,prev->state);    
          break; 
        case UNMATCHED_QUERY :   
          if( offi == 1 && offj == 0 && prev->state == MATCH65 ) {  
            pau->score = mat->lchs->b + (0);     
            continue;    
            }  
          if( offi == 1 && offj == 0 && prev->state == MATCH75 ) {  
            pau->score = mat->lchs->b + (0);     
            continue;    
            }  
          if( offi == 1 && offj == 0 && prev->state == MATCH85 ) {  
            pau->score = mat->lchs->b + (0);     
            continue;    
            }  
          if( offi == 1 && offj == 0 && prev->state == MATCH95 ) {  
            pau->score = mat->lchs->b + (0);     
            continue;    
            }  
          if( offi == 1 && offj == 0 && prev->state == UNMATCHED_QUERY ) {  
            pau->score = mat->lchs->u + (0);     
            continue;    
            }  
          warn("In recaluclating PackAln with state UNMATCHED_QUERY, from [%d,%d,%d], got a bad source state. Error!",offi,offj,prev->state);    
          break; 
        case UNMATCHED_TARGET :  
          if( offi == 0 && offj == 1 && prev->state == UNMATCHED_QUERY ) {  
            pau->score = mat->lchs->v + (0);     
            continue;    
            }  
          if( offi == 0 && offj == 1 && prev->state == UNMATCHED_TARGET )    {  
            pau->score = mat->lchs->u + (0);     
            continue;    
            }  
          warn("In recaluclating PackAln with state UNMATCHED_TARGET, from [%d,%d,%d], got a bad source state. Error!",offi,offj,prev->state);   
          break; 
        case (START+6) :     
          warn("In recaluclating PackAln with state START, got a bad source state. Error!"); 
          break; 
        case (END+6) :   
          if( offj == 0 && prev->state == MATCH65 )  {  
            /* i here comes from the previous state ;) - not the real one */ 
            i = prev->i; 
            pau->score = 0 + (0);    
            continue;    
            }  
          if( offj == 0 && prev->state == MATCH75 )  {  
            /* i here comes from the previous state ;) - not the real one */ 
            i = prev->i; 
            pau->score = 0 + (0);    
            continue;    
            }  
          if( offj == 0 && prev->state == MATCH85 )  {  
            /* i here comes from the previous state ;) - not the real one */ 
            i = prev->i; 
            pau->score = 0 + (0);    
            continue;    
            }  
          if( offj == 0 && prev->state == MATCH95 )  {  
            /* i here comes from the previous state ;) - not the real one */ 
            i = prev->i; 
            pau->score = 0 + (0);    
            continue;    
            }  
          warn("In recaluclating PackAln with state END, got a bad source state. Error!");   
          break; 
        default :    
          warn("In recaluclating PackAln got a bad recipient state. Error!");    
        }  
      prev = pau;    
      }  
    return;  
}    
/* divide and conquor macros are next */ 
#define LocalDnaMatchBlock_HIDDEN_MATRIX(thismatrix,i,j,state) (thismatrix->basematrix->matrix[(j-hiddenj+1)][(i+1)*6+state])    
#define LocalDnaMatchBlock_DC_SHADOW_MATRIX(thismatrix,i,j,state) (thismatrix->basematrix->matrix[((j+2)*8) % 16][(i+1)*6+state])    
#define LocalDnaMatchBlock_HIDDEN_SPECIAL(thismatrix,i,j,state) (thismatrix->basematrix->specmatrix[state][(j+1)])   
#define LocalDnaMatchBlock_DC_SHADOW_SPECIAL(thismatrix,i,j,state) (thismatrix->basematrix->specmatrix[state*8][(j+1)])  
#define LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(thismatrix,i,j,state,shadow) (thismatrix->basematrix->matrix[((((j+2)*8)+(shadow+1)) % 16)][(i+1)*6 + state]) 
#define LocalDnaMatchBlock_DC_SHADOW_SPECIAL_SP(thismatrix,i,j,state,shadow) (thismatrix->basematrix->specmatrix[state*8 +shadow+1][(j+1)])  
#define LocalDnaMatchBlock_DC_OPT_SHADOW_MATRIX(thismatrix,i,j,state) (score_pointers[(((j+1)% 1) * (leni+1) * 6) + ((i+1) * 6) + (state)])  
#define LocalDnaMatchBlock_DC_OPT_SHADOW_MATRIX_SP(thismatrix,i,j,state,shadow) (shadow_pointers[(((j+1)% 1) * (leni+1) * 48) + ((i+1) * 48) + (state * 8) + shadow+1])  
#define LocalDnaMatchBlock_DC_OPT_SHADOW_SPECIAL(thismatrix,i,j,state) (thismatrix->basematrix->specmatrix[state*8][(j+1)])  
/* Function:  allocate_Small_LocalDnaMatchBlock(query,target,lchs)
 *
 * Descrip:    This function allocates the LocalDnaMatchBlock structure
 *             and the basematrix area for a small memory implementations
 *             It calls /allocate_LocalDnaMatchBlock_only
 *
 *
 * Arg:         query [UNKN ] query data structure [ComplexSequence*]
 * Arg:        target [UNKN ] target data structure [ComplexSequence*]
 * Arg:          lchs [UNKN ] Resource [LocalCisHitScore*]
 *
 * Return [UNKN ]  Undocumented return value [LocalDnaMatchBlock *]
 *
 */
#define LocalDnaMatchBlock_DC_OPT_SHADOW_SPECIAL_SP(thismatrix,i,j,state,shadow) (thismatrix->basematrix->specmatrix[state*8 +shadow+1][(j+1)])  
LocalDnaMatchBlock * allocate_Small_LocalDnaMatchBlock(ComplexSequence* query,ComplexSequence* target ,LocalCisHitScore* lchs) 
{
    LocalDnaMatchBlock * out;    


    out = allocate_LocalDnaMatchBlock_only(query, target , lchs);    
    if( out == NULL )    
      return NULL;   
    out->basematrix = BaseMatrix_alloc_matrix_and_specials(16,(out->leni + 1) * 6,16,out->lenj+1);   
    if(out == NULL)  {  
      warn("Small shadow matrix LocalDnaMatchBlock cannot be allocated, (asking for 2 by %d main cells)",out->leni+2);   
      free_LocalDnaMatchBlock(out);  
      return NULL;   
      }  
    out->basematrix->type = BASEMATRIX_TYPE_SHADOW;  
    return out;  
}    


/* Function:  PackAln_calculate_Small_LocalDnaMatchBlock(mat,dpenv)
 *
 * Descrip:    This function calculates an alignment for LocalDnaMatchBlock structure in linear space
 *             If you want only the start/end points
 *             use /AlnRangeSet_calculate_Small_LocalDnaMatchBlock 
 *
 *             The function basically
 *               finds start/end points 
 *               foreach start/end point 
 *                 calls /full_dc_LocalDnaMatchBlock 
 *
 *
 * Arg:          mat [UNKN ] Undocumented argument [LocalDnaMatchBlock *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * PackAln_calculate_Small_LocalDnaMatchBlock(LocalDnaMatchBlock * mat,DPEnvelope * dpenv) 
{
    int endj;    
    int score;   
    PackAln * out;   
    PackAlnUnit * pau;   
    int starti;  
    int startj;  
    int startstate;  
    int stopi;   
    int stopj;   
    int stopstate;   
    int temp;    
    int donej;  /* This is for reporting, will be passed as a & arg in */ 
    int totalj; /* This also is for reporting, but as is not changed, can be passed by value */ 


    if( mat->basematrix->type != BASEMATRIX_TYPE_SHADOW )    {  
      warn("Could not calculate packaln small for LocalDnaMatchBlock due to wrong type of matrix");  
      return NULL;   
      }  


    out = PackAln_alloc_std();   


    start_reporting("Find start end points: ");  
    dc_optimised_start_end_calc_LocalDnaMatchBlock(mat,dpenv);   
    score = start_end_find_end_LocalDnaMatchBlock(mat,&endj);    
    out->score = score;  
    stopstate = END;
    
    /* No special to specials: one matrix alignment: simply remove and get */ 
    starti = LocalDnaMatchBlock_DC_SHADOW_SPECIAL_SP(mat,0,endj,END,0);  
    startj = LocalDnaMatchBlock_DC_SHADOW_SPECIAL_SP(mat,0,endj,END,1);  
    startstate = LocalDnaMatchBlock_DC_SHADOW_SPECIAL_SP(mat,0,endj,END,2);  
    stopi = LocalDnaMatchBlock_DC_SHADOW_SPECIAL_SP(mat,0,endj,END,3);   
    stopj = LocalDnaMatchBlock_DC_SHADOW_SPECIAL_SP(mat,0,endj,END,4);   
    stopstate = LocalDnaMatchBlock_DC_SHADOW_SPECIAL_SP(mat,0,endj,END,5);   
    temp = LocalDnaMatchBlock_DC_SHADOW_SPECIAL_SP(mat,0,endj,END,6);    
    log_full_error(REPORT,0,"[%d,%d][%d,%d] Score %d",starti,startj,stopi,stopj,score);  
    stop_reporting();    
    start_reporting("Recovering alignment: ");   


    /* Figuring how much j we have to align for reporting purposes */ 
    donej = 0;   
    totalj = stopj - startj; 
    full_dc_LocalDnaMatchBlock(mat,starti,startj,startstate,stopi,stopj,stopstate,out,&donej,totalj,dpenv);  


    /* Although we have no specials, need to get start. Better to check than assume */ 


    max_matrix_to_special_LocalDnaMatchBlock(mat,starti,startj,startstate,temp,&stopi,&stopj,&stopstate,&temp,NULL); 
    if( stopi == LocalDnaMatchBlock_READ_OFF_ERROR || stopstate != START )   {  
      warn("Problem in reading off special state system, hit a non start state (or an internal error) in a single alignment mode");  
      invert_PackAln(out);   
      recalculate_PackAln_LocalDnaMatchBlock(out,mat);   
      return out;    
      }  


    /* Ok. Put away start start... */ 
    pau = PackAlnUnit_alloc();   
    pau->i = stopi;  
    pau->j = stopj;  
    pau->state = stopstate + 6;  
    add_PackAln(out,pau);    


    log_full_error(REPORT,0,"Alignment recovered");  
    stop_reporting();    
    invert_PackAln(out); 
    recalculate_PackAln_LocalDnaMatchBlock(out,mat); 
    return out;  


}    


/* Function:  AlnRangeSet_calculate_Small_LocalDnaMatchBlock(mat)
 *
 * Descrip:    This function calculates an alignment for LocalDnaMatchBlock structure in linear space
 *             If you want the full alignment, use /PackAln_calculate_Small_LocalDnaMatchBlock 
 *             If you have already got the full alignment, but want the range set, use /AlnRangeSet_from_PackAln_LocalDnaMatchBlock
 *             If you have got the small matrix but not the alignment, use /AlnRangeSet_from_LocalDnaMatchBlock 
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [LocalDnaMatchBlock *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * AlnRangeSet_calculate_Small_LocalDnaMatchBlock(LocalDnaMatchBlock * mat) 
{
    AlnRangeSet * out;   


    start_reporting("Find start end points: ");  
    dc_optimised_start_end_calc_LocalDnaMatchBlock(mat,NULL);    
    log_full_error(REPORT,0,"Calculated");   


    out = AlnRangeSet_from_LocalDnaMatchBlock(mat);  
    return out;  
}    


/* Function:  AlnRangeSet_from_LocalDnaMatchBlock(mat)
 *
 * Descrip:    This function reads off a start/end structure
 *             for LocalDnaMatchBlock structure in linear space
 *             If you want the full alignment use
 *             /PackAln_calculate_Small_LocalDnaMatchBlock 
 *             If you have not calculated the matrix use
 *             /AlnRange_calculate_Small_LocalDnaMatchBlock
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [LocalDnaMatchBlock *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * AlnRangeSet_from_LocalDnaMatchBlock(LocalDnaMatchBlock * mat) 
{
    AlnRangeSet * out;   
    AlnRange * temp; 
    int jpos;    
    int state;   


    if( mat->basematrix->type != BASEMATRIX_TYPE_SHADOW) {  
      warn("Bad error! - non shadow matrix type in AlnRangeSet_from_LocalDnaMatchBlock");    
      return NULL;   
      }  


    out = AlnRangeSet_alloc_std();   
    /* Find the end position */ 
    out->score = start_end_find_end_LocalDnaMatchBlock(mat,&jpos);   
    state = END; 


    while( (temp = AlnRange_build_LocalDnaMatchBlock(mat,jpos,state,&jpos,&state)) != NULL)  
      add_AlnRangeSet(out,temp); 
    return out;  
}    


/* Function:  AlnRange_build_LocalDnaMatchBlock(mat,stopj,stopspecstate,startj,startspecstate)
 *
 * Descrip:    This function calculates a single start/end set in linear space
 *             Really a sub-routine for /AlnRangeSet_from_PackAln_LocalDnaMatchBlock
 *
 *
 * Arg:                   mat [UNKN ] Undocumented argument [LocalDnaMatchBlock *]
 * Arg:                 stopj [UNKN ] Undocumented argument [int]
 * Arg:         stopspecstate [UNKN ] Undocumented argument [int]
 * Arg:                startj [UNKN ] Undocumented argument [int *]
 * Arg:        startspecstate [UNKN ] Undocumented argument [int *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRange *]
 *
 */
AlnRange * AlnRange_build_LocalDnaMatchBlock(LocalDnaMatchBlock * mat,int stopj,int stopspecstate,int * startj,int * startspecstate) 
{
    AlnRange * out;  
    int jpos;    
    int state;   


    if( mat->basematrix->type != BASEMATRIX_TYPE_SHADOW) {  
      warn("Bad error! - non shadow matrix type in AlnRangeSet_from_LocalDnaMatchBlock");    
      return NULL;   
      }  


    /* Assumme that we have specials (we should!). Read back along the specials till we have the finish point */ 
    if( read_special_strip_LocalDnaMatchBlock(mat,0,stopj,stopspecstate,&jpos,&state,NULL) == FALSE) {  
      warn("In AlnRanger_build_LocalDnaMatchBlock alignment ending at %d, unable to read back specials. Will (evenutally) return a partial range set... BEWARE!",stopj); 
      return NULL;   
      }  
    if( state == START || jpos <= 0) 
      return NULL;   


    out = AlnRange_alloc();  


    out->starti = LocalDnaMatchBlock_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,0);   
    out->startj = LocalDnaMatchBlock_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,1);   
    out->startstate = LocalDnaMatchBlock_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,2);   
    out->stopi = LocalDnaMatchBlock_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,3);    
    out->stopj = LocalDnaMatchBlock_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,4);    
    out->stopstate = LocalDnaMatchBlock_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,5);    
    out->startscore = LocalDnaMatchBlock_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,6);   
    out->stopscore = LocalDnaMatchBlock_DC_SHADOW_SPECIAL(mat,0,jpos,state); 


    /* Now, we have to figure out where this state came from in the specials */ 
    max_matrix_to_special_LocalDnaMatchBlock(mat,out->starti,out->startj,out->startstate,out->startscore,&jpos,startj,startspecstate,&state,NULL);   
    if( jpos == LocalDnaMatchBlock_READ_OFF_ERROR)   {  
      warn("In AlnRange_build_LocalDnaMatchBlock alignment ending at %d, with aln range between %d-%d in j, unable to find source special, returning this range, but this could get tricky!",stopj,out->startj,out->stopj);  
      return out;    
      }  


    /* Put in the correct score for startstate, from the special */ 
    out->startscore = LocalDnaMatchBlock_DC_SHADOW_SPECIAL(mat,0,*startj,*startspecstate);   
    /* The correct j coords have been put into startj, startspecstate... so just return out */ 
    return out;  
}    


/* Function:  read_hidden_LocalDnaMatchBlock(mat,starti,startj,startstate,stopi,stopj,stopstate,out)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [LocalDnaMatchBlock *]
 * Arg:            starti [UNKN ] Undocumented argument [int]
 * Arg:            startj [UNKN ] Undocumented argument [int]
 * Arg:        startstate [UNKN ] Undocumented argument [int]
 * Arg:             stopi [UNKN ] Undocumented argument [int]
 * Arg:             stopj [UNKN ] Undocumented argument [int]
 * Arg:         stopstate [UNKN ] Undocumented argument [int]
 * Arg:               out [UNKN ] Undocumented argument [PackAln *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean read_hidden_LocalDnaMatchBlock(LocalDnaMatchBlock * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out) 
{
    int i;   
    int j;   
    int state;   
    int cellscore;   
    int isspecial;   
    /* We don't need hiddenj here, 'cause matrix access handled by max funcs */ 
    PackAlnUnit * pau;   


    /* stop position is on the path */ 
    i = stopi;   
    j = stopj;   
    state= stopstate;    
    isspecial = FALSE;   


    while( i >= starti && j >= startj)   {  
      /* Put away current i,j,state */ 
      pau = PackAlnUnit_alloc();/* Should deal with memory overflow */ 
      pau->i = i;    
      pau->j = j;    
      pau->state =  state;   
      add_PackAln(out,pau);  


      max_hidden_LocalDnaMatchBlock(mat,startj,i,j,state,isspecial,&i,&j,&state,&isspecial,&cellscore);  


      if( i == LocalDnaMatchBlock_READ_OFF_ERROR)    {  
        warn("In LocalDnaMatchBlock hidden read off, between %d:%d,%d:%d - at got bad read off. Problem!",starti,startj,stopi,stopj);    
        return FALSE;    
        }  


      if( i == starti && j == startj && state == startstate) {  
/* Put away final state (start of this block) */ 
        pau = PackAlnUnit_alloc();  /* Should deal with memory overflow */ 
        pau->i = i;  
        pau->j = j;  
        pau->state =  state; 
        add_PackAln(out,pau);    
          return TRUE;   
        }  
      if( i == starti && j == startj)    {  
        warn("In LocalDnaMatchBlock hidden read off, between %d:%d,%d:%d - hit start cell, but not in start state. Can't be good!.",starti,startj,stopi,stopj);  
        return FALSE;    
        }  
      }  
    warn("In LocalDnaMatchBlock hidden read off, between %d:%d,%d:%d - gone past start cell (now in %d,%d,%d), can't be good news!.",starti,startj,stopi,stopj,i,j,state);   
    return FALSE;    
}    


/* Function:  max_hidden_LocalDnaMatchBlock(mat,hiddenj,i,j,state,isspecial,reti,retj,retstate,retspecial,cellscore)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [LocalDnaMatchBlock *]
 * Arg:           hiddenj [UNKN ] Undocumented argument [int]
 * Arg:                 i [UNKN ] Undocumented argument [int]
 * Arg:                 j [UNKN ] Undocumented argument [int]
 * Arg:             state [UNKN ] Undocumented argument [int]
 * Arg:         isspecial [UNKN ] Undocumented argument [boolean]
 * Arg:              reti [UNKN ] Undocumented argument [int *]
 * Arg:              retj [UNKN ] Undocumented argument [int *]
 * Arg:          retstate [UNKN ] Undocumented argument [int *]
 * Arg:        retspecial [UNKN ] Undocumented argument [boolean *]
 * Arg:         cellscore [UNKN ] Undocumented argument [int *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int max_hidden_LocalDnaMatchBlock(LocalDnaMatchBlock * mat,int hiddenj,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore) 
{
    register int temp;   
    register int cscore; 


    *reti = (*retj) = (*retstate) = LocalDnaMatchBlock_READ_OFF_ERROR;   


    if( i < 0 || j < 0 || i > mat->query->seq->len || j > mat->target->seq->len) {  
      warn("In LocalDnaMatchBlock matrix special read off - out of bounds on matrix [i,j is %d,%d state %d in standard matrix]",i,j,state);  
      return -1; 
      }  


    /* Then you have to select the correct switch statement to figure out the readoff      */ 
    /* Somewhat odd - reverse the order of calculation and return as soon as it is correct */ 
    cscore = LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i,j,state);    
    switch(state)    { /*Switch state */ 
      case MATCH65 :     
        temp = cscore - ((mat->lchs->comp65->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->v)) -  (0); 
        if( temp == LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i - 1,j - 1,UNMATCHED_TARGET) ) {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = UNMATCHED_TARGET;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i-1,j-1,UNMATCHED_TARGET);    
            }  
          return LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i - 1,j - 1,UNMATCHED_TARGET);     
          }  
        /* Not allowing special sources.. skipping START */ 
        temp = cscore - (mat->lchs->g) -  (0);   
        if( temp == LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i - 1,j - 0,MATCH65) )  {  
          *reti = i - 1; 
          *retj = j - 0; 
          *retstate = MATCH65;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i-1,j-0,MATCH65); 
            }  
          return LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i - 1,j - 0,MATCH65);  
          }  
        temp = cscore - (mat->lchs->g) -  (0);   
        if( temp == LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i - 0,j - 1,MATCH65) )  {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = MATCH65;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i-0,j-1,MATCH65); 
            }  
          return LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i - 0,j - 1,MATCH65);  
          }  
        temp = cscore - ((mat->lchs->comp65->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->s)) -  (0); 
        if( temp == LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i - 1,j - 1,MATCH65) )  {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = MATCH65;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i-1,j-1,MATCH65); 
            }  
          return LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i - 1,j - 1,MATCH65);  
          }  
        warn("Major problem (!) - in LocalDnaMatchBlock read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      case MATCH75 :     
        temp = cscore - ((mat->lchs->comp75->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->v)) -  (0); 
        if( temp == LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i - 1,j - 1,UNMATCHED_TARGET) ) {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = UNMATCHED_TARGET;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i-1,j-1,UNMATCHED_TARGET);    
            }  
          return LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i - 1,j - 1,UNMATCHED_TARGET);     
          }  
        /* Not allowing special sources.. skipping START */ 
        temp = cscore - (mat->lchs->g) -  (0);   
        if( temp == LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i - 1,j - 0,MATCH75) )  {  
          *reti = i - 1; 
          *retj = j - 0; 
          *retstate = MATCH75;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i-1,j-0,MATCH75); 
            }  
          return LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i - 1,j - 0,MATCH75);  
          }  
        temp = cscore - (mat->lchs->g) -  (0);   
        if( temp == LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i - 0,j - 1,MATCH75) )  {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = MATCH75;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i-0,j-1,MATCH75); 
            }  
          return LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i - 0,j - 1,MATCH75);  
          }  
        temp = cscore - ((mat->lchs->comp75->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->s)) -  (0); 
        if( temp == LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i - 1,j - 1,MATCH75) )  {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = MATCH75;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i-1,j-1,MATCH75); 
            }  
          return LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i - 1,j - 1,MATCH75);  
          }  
        warn("Major problem (!) - in LocalDnaMatchBlock read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      case MATCH85 :     
        temp = cscore - ((mat->lchs->comp85->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->v)) -  (0); 
        if( temp == LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i - 1,j - 1,UNMATCHED_TARGET) ) {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = UNMATCHED_TARGET;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i-1,j-1,UNMATCHED_TARGET);    
            }  
          return LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i - 1,j - 1,UNMATCHED_TARGET);     
          }  
        /* Not allowing special sources.. skipping START */ 
        temp = cscore - (mat->lchs->g) -  (0);   
        if( temp == LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i - 1,j - 0,MATCH85) )  {  
          *reti = i - 1; 
          *retj = j - 0; 
          *retstate = MATCH85;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i-1,j-0,MATCH85); 
            }  
          return LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i - 1,j - 0,MATCH85);  
          }  
        temp = cscore - (mat->lchs->g) -  (0);   
        if( temp == LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i - 0,j - 1,MATCH85) )  {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = MATCH85;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i-0,j-1,MATCH85); 
            }  
          return LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i - 0,j - 1,MATCH85);  
          }  
        temp = cscore - ((mat->lchs->comp85->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->s)) -  (0); 
        if( temp == LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i - 1,j - 1,MATCH85) )  {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = MATCH85;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i-1,j-1,MATCH85); 
            }  
          return LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i - 1,j - 1,MATCH85);  
          }  
        warn("Major problem (!) - in LocalDnaMatchBlock read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      case MATCH95 :     
        temp = cscore - ((mat->lchs->comp95->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->v)) -  (0); 
        if( temp == LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i - 1,j - 1,UNMATCHED_TARGET) ) {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = UNMATCHED_TARGET;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i-1,j-1,UNMATCHED_TARGET);    
            }  
          return LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i - 1,j - 1,UNMATCHED_TARGET);     
          }  
        /* Not allowing special sources.. skipping START */ 
        temp = cscore - (mat->lchs->g) -  (0);   
        if( temp == LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i - 1,j - 0,MATCH95) )  {  
          *reti = i - 1; 
          *retj = j - 0; 
          *retstate = MATCH95;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i-1,j-0,MATCH95); 
            }  
          return LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i - 1,j - 0,MATCH95);  
          }  
        temp = cscore - (mat->lchs->g) -  (0);   
        if( temp == LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i - 0,j - 1,MATCH95) )  {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = MATCH95;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i-0,j-1,MATCH95); 
            }  
          return LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i - 0,j - 1,MATCH95);  
          }  
        temp = cscore - ((mat->lchs->comp95->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->s)) -  (0); 
        if( temp == LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i - 1,j - 1,MATCH95) )  {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = MATCH95;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i-1,j-1,MATCH95); 
            }  
          return LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i - 1,j - 1,MATCH95);  
          }  
        warn("Major problem (!) - in LocalDnaMatchBlock read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      case UNMATCHED_QUERY :     
        temp = cscore - (mat->lchs->u) -  (0);   
        if( temp == LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i - 1,j - 0,UNMATCHED_QUERY) )  {  
          *reti = i - 1; 
          *retj = j - 0; 
          *retstate = UNMATCHED_QUERY;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i-1,j-0,UNMATCHED_QUERY); 
            }  
          return LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i - 1,j - 0,UNMATCHED_QUERY);  
          }  
        temp = cscore - (mat->lchs->b) -  (0);   
        if( temp == LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i - 1,j - 0,MATCH95) )  {  
          *reti = i - 1; 
          *retj = j - 0; 
          *retstate = MATCH95;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i-1,j-0,MATCH95); 
            }  
          return LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i - 1,j - 0,MATCH95);  
          }  
        temp = cscore - (mat->lchs->b) -  (0);   
        if( temp == LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i - 1,j - 0,MATCH85) )  {  
          *reti = i - 1; 
          *retj = j - 0; 
          *retstate = MATCH85;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i-1,j-0,MATCH85); 
            }  
          return LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i - 1,j - 0,MATCH85);  
          }  
        temp = cscore - (mat->lchs->b) -  (0);   
        if( temp == LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i - 1,j - 0,MATCH75) )  {  
          *reti = i - 1; 
          *retj = j - 0; 
          *retstate = MATCH75;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i-1,j-0,MATCH75); 
            }  
          return LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i - 1,j - 0,MATCH75);  
          }  
        temp = cscore - (mat->lchs->b) -  (0);   
        if( temp == LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i - 1,j - 0,MATCH65) )  {  
          *reti = i - 1; 
          *retj = j - 0; 
          *retstate = MATCH65;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i-1,j-0,MATCH65); 
            }  
          return LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i - 1,j - 0,MATCH65);  
          }  
        warn("Major problem (!) - in LocalDnaMatchBlock read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      case UNMATCHED_TARGET :    
        temp = cscore - (mat->lchs->u) -  (0);   
        if( temp == LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i - 0,j - 1,UNMATCHED_TARGET) ) {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = UNMATCHED_TARGET;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i-0,j-1,UNMATCHED_TARGET);    
            }  
          return LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i - 0,j - 1,UNMATCHED_TARGET);     
          }  
        temp = cscore - (mat->lchs->v) -  (0);   
        if( temp == LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i - 0,j - 1,UNMATCHED_QUERY) )  {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = UNMATCHED_QUERY;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i-0,j-1,UNMATCHED_QUERY); 
            }  
          return LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i - 0,j - 1,UNMATCHED_QUERY);  
          }  
        warn("Major problem (!) - in LocalDnaMatchBlock read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      default:   
        warn("Major problem (!) - in LocalDnaMatchBlock read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      } /* end of Switch state  */ 
}    


/* Function:  read_special_strip_LocalDnaMatchBlock(mat,stopi,stopj,stopstate,startj,startstate,out)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [LocalDnaMatchBlock *]
 * Arg:             stopi [UNKN ] Undocumented argument [int]
 * Arg:             stopj [UNKN ] Undocumented argument [int]
 * Arg:         stopstate [UNKN ] Undocumented argument [int]
 * Arg:            startj [UNKN ] Undocumented argument [int *]
 * Arg:        startstate [UNKN ] Undocumented argument [int *]
 * Arg:               out [UNKN ] Undocumented argument [PackAln *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean read_special_strip_LocalDnaMatchBlock(LocalDnaMatchBlock * mat,int stopi,int stopj,int stopstate,int * startj,int * startstate,PackAln * out) 
{
    int i;   
    int j;   
    int state;   
    int cellscore;   
    int isspecial;   
    PackAlnUnit * pau;   


    /* stop position is on the path */ 
    i = stopi;   
    j = stopj;   
    state= stopstate;    
    isspecial = TRUE;    


    /* Loop until state has the same j as its stop in shadow pointers */ 
    /* This will be the state is came out from, OR it has hit !start */ 
    /* We may not want to get the alignment, in which case out will be NULL */ 
    while( j > LocalDnaMatchBlock_DC_SHADOW_SPECIAL_SP(mat,i,j,state,4) && state != START)   { /*while more specials to eat up*/ 
      /* Put away current state, if we should */ 
      if(out != NULL)    {  
        pau = PackAlnUnit_alloc();  /* Should deal with memory overflow */ 
        pau->i = i;  
        pau->j = j;  
        pau->state =  state + 6; 
        add_PackAln(out,pau);    
        }  


      max_special_strip_LocalDnaMatchBlock(mat,i,j,state,isspecial,&i,&j,&state,&isspecial,&cellscore);  
      if( i == LocalDnaMatchBlock_READ_OFF_ERROR)    {  
        warn("In special strip read LocalDnaMatchBlock, got a bad read off error. Sorry!");  
        return FALSE;    
        }  
      } /* end of while more specials to eat up */ 


    /* check to see we have not gone too far! */ 
    if( state != START && j < LocalDnaMatchBlock_DC_SHADOW_SPECIAL_SP(mat,i,j,state,4))  {  
      warn("In special strip read LocalDnaMatchBlock, at special [%d] state [%d] overshot!",j,state);    
      return FALSE;  
      }  
    /* Put away last state */ 
    if(out != NULL)  {  
      pau = PackAlnUnit_alloc();/* Should deal with memory overflow */ 
      pau->i = i;    
      pau->j = j;    
      pau->state =  state + 6;   
      add_PackAln(out,pau);  
      }  


    /* Put away where we are in startj and startstate */ 
    *startj = j; 
    *startstate = state; 
    return TRUE; 
}    


/* Function:  max_special_strip_LocalDnaMatchBlock(mat,i,j,state,isspecial,reti,retj,retstate,retspecial,cellscore)
 *
 * Descrip:    A pretty intense internal function. Deals with read-off only in specials
 *
 *
 * Arg:               mat [UNKN ] Undocumented argument [LocalDnaMatchBlock *]
 * Arg:                 i [UNKN ] Undocumented argument [int]
 * Arg:                 j [UNKN ] Undocumented argument [int]
 * Arg:             state [UNKN ] Undocumented argument [int]
 * Arg:         isspecial [UNKN ] Undocumented argument [boolean]
 * Arg:              reti [UNKN ] Undocumented argument [int *]
 * Arg:              retj [UNKN ] Undocumented argument [int *]
 * Arg:          retstate [UNKN ] Undocumented argument [int *]
 * Arg:        retspecial [UNKN ] Undocumented argument [boolean *]
 * Arg:         cellscore [UNKN ] Undocumented argument [int *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int max_special_strip_LocalDnaMatchBlock(LocalDnaMatchBlock * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore) 
{
    int temp;    
    int cscore;  


    *reti = (*retj) = (*retstate) = LocalDnaMatchBlock_READ_OFF_ERROR;   
    if( isspecial == FALSE ) {  
      warn("In special strip max function for LocalDnaMatchBlock, got a non special start point. Problem! (bad!)");  
      return (-1);   
      }  


    if( j < 0 || j > mat->target->seq->len)  {  
      warn("In LocalDnaMatchBlock matrix special read off - out of bounds on matrix [j is %d in special]",j);    
      return -1; 
      }  


    cscore = LocalDnaMatchBlock_DC_SHADOW_SPECIAL(mat,i,j,state);    
    switch(state)    { /*switch on special states*/ 
      case START :   
      case END :     
        /* Source MATCH95 is not a special */ 
        /* Source MATCH85 is not a special */ 
        /* Source MATCH75 is not a special */ 
        /* Source MATCH65 is not a special */ 
      default:   
        warn("Major problem (!) - in LocalDnaMatchBlock special strip read off, position %d,%d state %d no source found  dropped into default on source switch!",i,j,state); 
        return (-1); 
      } /* end of switch on special states */ 
}    


/* Function:  max_matrix_to_special_LocalDnaMatchBlock(mat,i,j,state,cscore,reti,retj,retstate,retspecial,cellscore)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [LocalDnaMatchBlock *]
 * Arg:                 i [UNKN ] Undocumented argument [int]
 * Arg:                 j [UNKN ] Undocumented argument [int]
 * Arg:             state [UNKN ] Undocumented argument [int]
 * Arg:            cscore [UNKN ] Undocumented argument [int]
 * Arg:              reti [UNKN ] Undocumented argument [int *]
 * Arg:              retj [UNKN ] Undocumented argument [int *]
 * Arg:          retstate [UNKN ] Undocumented argument [int *]
 * Arg:        retspecial [UNKN ] Undocumented argument [boolean *]
 * Arg:         cellscore [UNKN ] Undocumented argument [int *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int max_matrix_to_special_LocalDnaMatchBlock(LocalDnaMatchBlock * mat,int i,int j,int state,int cscore,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore) 
{
    int temp;    
    *reti = (*retj) = (*retstate) = LocalDnaMatchBlock_READ_OFF_ERROR;   


    if( j < 0 || j > mat->lenj)  {  
      warn("In LocalDnaMatchBlock matrix to special read off - out of bounds on matrix [j is %d in special]",j); 
      return -1; 
      }  


    switch(state)    { /*Switch state */ 
      case MATCH65 :     
        /* Source UNMATCHED_TARGET is not a special, should not get here! */ 
        temp = cscore - ((mat->lchs->comp65->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->v)) -  (0);     
        if( temp == LocalDnaMatchBlock_DC_SHADOW_SPECIAL(mat,i - 1,j - 1,START) )    {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = START; 
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - LocalDnaMatchBlock_DC_SHADOW_SPECIAL(mat,i-1,j-1,START);   
            }  
          return LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i - 1,j - 1,START) ;    
          }  
        /* Source MATCH65 is not a special, should not get here! */ 
        /* Source MATCH65 is not a special, should not get here! */ 
        /* Source MATCH65 is not a special, should not get here! */ 
        warn("Major problem (!) - in LocalDnaMatchBlock matrix to special read off, position %d,%d state %d no source found!",i,j,state);    
        return (-1); 
      case MATCH75 :     
        /* Source UNMATCHED_TARGET is not a special, should not get here! */ 
        temp = cscore - ((mat->lchs->comp75->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->v)) -  (0);     
        if( temp == LocalDnaMatchBlock_DC_SHADOW_SPECIAL(mat,i - 1,j - 1,START) )    {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = START; 
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - LocalDnaMatchBlock_DC_SHADOW_SPECIAL(mat,i-1,j-1,START);   
            }  
          return LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i - 1,j - 1,START) ;    
          }  
        /* Source MATCH75 is not a special, should not get here! */ 
        /* Source MATCH75 is not a special, should not get here! */ 
        /* Source MATCH75 is not a special, should not get here! */ 
        warn("Major problem (!) - in LocalDnaMatchBlock matrix to special read off, position %d,%d state %d no source found!",i,j,state);    
        return (-1); 
      case MATCH85 :     
        /* Source UNMATCHED_TARGET is not a special, should not get here! */ 
        temp = cscore - ((mat->lchs->comp85->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->v)) -  (0);     
        if( temp == LocalDnaMatchBlock_DC_SHADOW_SPECIAL(mat,i - 1,j - 1,START) )    {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = START; 
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - LocalDnaMatchBlock_DC_SHADOW_SPECIAL(mat,i-1,j-1,START);   
            }  
          return LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i - 1,j - 1,START) ;    
          }  
        /* Source MATCH85 is not a special, should not get here! */ 
        /* Source MATCH85 is not a special, should not get here! */ 
        /* Source MATCH85 is not a special, should not get here! */ 
        warn("Major problem (!) - in LocalDnaMatchBlock matrix to special read off, position %d,%d state %d no source found!",i,j,state);    
        return (-1); 
      case MATCH95 :     
        /* Source UNMATCHED_TARGET is not a special, should not get here! */ 
        temp = cscore - ((mat->lchs->comp95->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->v)) -  (0);     
        if( temp == LocalDnaMatchBlock_DC_SHADOW_SPECIAL(mat,i - 1,j - 1,START) )    {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = START; 
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - LocalDnaMatchBlock_DC_SHADOW_SPECIAL(mat,i-1,j-1,START);   
            }  
          return LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i - 1,j - 1,START) ;    
          }  
        /* Source MATCH95 is not a special, should not get here! */ 
        /* Source MATCH95 is not a special, should not get here! */ 
        /* Source MATCH95 is not a special, should not get here! */ 
        warn("Major problem (!) - in LocalDnaMatchBlock matrix to special read off, position %d,%d state %d no source found!",i,j,state);    
        return (-1); 
      case UNMATCHED_QUERY :     
        /* Source UNMATCHED_QUERY is not a special, should not get here! */ 
        /* Source MATCH95 is not a special, should not get here! */ 
        /* Source MATCH85 is not a special, should not get here! */ 
        /* Source MATCH75 is not a special, should not get here! */ 
        /* Source MATCH65 is not a special, should not get here! */ 
        warn("Major problem (!) - in LocalDnaMatchBlock matrix to special read off, position %d,%d state %d no source found!",i,j,state);    
        return (-1); 
      case UNMATCHED_TARGET :    
        /* Source UNMATCHED_TARGET is not a special, should not get here! */ 
        /* Source UNMATCHED_QUERY is not a special, should not get here! */ 
        warn("Major problem (!) - in LocalDnaMatchBlock matrix to special read off, position %d,%d state %d no source found!",i,j,state);    
        return (-1); 
      default:   
        warn("Major problem (!) - in LocalDnaMatchBlock read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      } /* end of Switch state  */ 


}    


/* Function:  calculate_hidden_LocalDnaMatchBlock(mat,starti,startj,startstate,stopi,stopj,stopstate,dpenv)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [LocalDnaMatchBlock *]
 * Arg:            starti [UNKN ] Undocumented argument [int]
 * Arg:            startj [UNKN ] Undocumented argument [int]
 * Arg:        startstate [UNKN ] Undocumented argument [int]
 * Arg:             stopi [UNKN ] Undocumented argument [int]
 * Arg:             stopj [UNKN ] Undocumented argument [int]
 * Arg:         stopstate [UNKN ] Undocumented argument [int]
 * Arg:             dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 */
void calculate_hidden_LocalDnaMatchBlock(LocalDnaMatchBlock * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv) 
{
    register int i;  
    register int j;  
    register int score;  
    register int temp;   
    register int hiddenj;    


    hiddenj = startj;    


    init_hidden_LocalDnaMatchBlock(mat,starti,startj,stopi,stopj);   


    LocalDnaMatchBlock_HIDDEN_MATRIX(mat,starti,startj,startstate) = 0;  


    for(j=startj;j<=stopj;j++)   {  
      for(i=starti;i<=stopi;i++) {  
        /* Should *not* do very first cell as this is the one set to zero in one state! */ 
        if( i == starti && j == startj ) 
          continue;  
        if( dpenv != NULL && is_in_DPEnvelope(dpenv,i,j) == FALSE )  { /*Is not in envelope*/ 
          LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i,j,MATCH65) = NEGI;  
          LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i,j,MATCH75) = NEGI;  
          LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i,j,MATCH85) = NEGI;  
          LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i,j,MATCH95) = NEGI;  
          LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i,j,UNMATCHED_QUERY) = NEGI;  
          LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i,j,UNMATCHED_TARGET) = NEGI;     
          continue;  
          } /* end of Is not in envelope */ 


        /* For state MATCH65 */ 
        /* setting first movement to score */ 
        score = LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i-1,j-1,MATCH65) + (mat->lchs->comp65->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->s);  
        /* From state MATCH65 to state MATCH65 */ 
        temp = LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i-0,j-1,MATCH65) + mat->lchs->g;     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH65 to state MATCH65 */ 
        temp = LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i-1,j-0,MATCH65) + mat->lchs->g;     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state UNMATCHED_TARGET to state MATCH65 */ 
        temp = LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i-1,j-1,UNMATCHED_TARGET) + (mat->lchs->comp65->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->v);  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for MATCH65 */ 
        /* Add any movement independant score and put away */ 
         LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i,j,MATCH65) = score;  
        /* Finished calculating state MATCH65 */ 


        /* For state MATCH75 */ 
        /* setting first movement to score */ 
        score = LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i-1,j-1,MATCH75) + (mat->lchs->comp75->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->s);  
        /* From state MATCH75 to state MATCH75 */ 
        temp = LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i-0,j-1,MATCH75) + mat->lchs->g;     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH75 to state MATCH75 */ 
        temp = LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i-1,j-0,MATCH75) + mat->lchs->g;     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state UNMATCHED_TARGET to state MATCH75 */ 
        temp = LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i-1,j-1,UNMATCHED_TARGET) + (mat->lchs->comp75->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->v);  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for MATCH75 */ 
        /* Add any movement independant score and put away */ 
         LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i,j,MATCH75) = score;  
        /* Finished calculating state MATCH75 */ 


        /* For state MATCH85 */ 
        /* setting first movement to score */ 
        score = LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i-1,j-1,MATCH85) + (mat->lchs->comp85->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->s);  
        /* From state MATCH85 to state MATCH85 */ 
        temp = LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i-0,j-1,MATCH85) + mat->lchs->g;     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH85 to state MATCH85 */ 
        temp = LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i-1,j-0,MATCH85) + mat->lchs->g;     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state UNMATCHED_TARGET to state MATCH85 */ 
        temp = LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i-1,j-1,UNMATCHED_TARGET) + (mat->lchs->comp85->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->v);  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for MATCH85 */ 
        /* Add any movement independant score and put away */ 
         LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i,j,MATCH85) = score;  
        /* Finished calculating state MATCH85 */ 


        /* For state MATCH95 */ 
        /* setting first movement to score */ 
        score = LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i-1,j-1,MATCH95) + (mat->lchs->comp95->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->s);  
        /* From state MATCH95 to state MATCH95 */ 
        temp = LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i-0,j-1,MATCH95) + mat->lchs->g;     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH95 to state MATCH95 */ 
        temp = LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i-1,j-0,MATCH95) + mat->lchs->g;     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state UNMATCHED_TARGET to state MATCH95 */ 
        temp = LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i-1,j-1,UNMATCHED_TARGET) + (mat->lchs->comp95->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->v);  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for MATCH95 */ 
        /* Add any movement independant score and put away */ 
         LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i,j,MATCH95) = score;  
        /* Finished calculating state MATCH95 */ 


        /* For state UNMATCHED_QUERY */ 
        /* setting first movement to score */ 
        score = LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i-1,j-0,MATCH65) + mat->lchs->b;    
        /* From state MATCH75 to state UNMATCHED_QUERY */ 
        temp = LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i-1,j-0,MATCH75) + mat->lchs->b;     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH85 to state UNMATCHED_QUERY */ 
        temp = LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i-1,j-0,MATCH85) + mat->lchs->b;     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH95 to state UNMATCHED_QUERY */ 
        temp = LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i-1,j-0,MATCH95) + mat->lchs->b;     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state UNMATCHED_QUERY to state UNMATCHED_QUERY */ 
        temp = LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i-1,j-0,UNMATCHED_QUERY) + mat->lchs->u;     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for UNMATCHED_QUERY */ 
        /* Add any movement independant score and put away */ 
         LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i,j,UNMATCHED_QUERY) = score;  
        /* Finished calculating state UNMATCHED_QUERY */ 


        /* For state UNMATCHED_TARGET */ 
        /* setting first movement to score */ 
        score = LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i-0,j-1,UNMATCHED_QUERY) + mat->lchs->v;    
        /* From state UNMATCHED_TARGET to state UNMATCHED_TARGET */ 
        temp = LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i-0,j-1,UNMATCHED_TARGET) + mat->lchs->u;    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for UNMATCHED_TARGET */ 
        /* Add any movement independant score and put away */ 
         LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i,j,UNMATCHED_TARGET) = score; 
        /* Finished calculating state UNMATCHED_TARGET */ 
        }  
      }  


    return;  
}    


/* Function:  init_hidden_LocalDnaMatchBlock(mat,starti,startj,stopi,stopj)
 *
 * Descrip: No Description
 *
 * Arg:           mat [UNKN ] Undocumented argument [LocalDnaMatchBlock *]
 * Arg:        starti [UNKN ] Undocumented argument [int]
 * Arg:        startj [UNKN ] Undocumented argument [int]
 * Arg:         stopi [UNKN ] Undocumented argument [int]
 * Arg:         stopj [UNKN ] Undocumented argument [int]
 *
 */
void init_hidden_LocalDnaMatchBlock(LocalDnaMatchBlock * mat,int starti,int startj,int stopi,int stopj) 
{
    register int i;  
    register int j;  
    register int hiddenj;    


    hiddenj = startj;    
    for(j=(startj-1);j<=stopj;j++)   {  
      for(i=(starti-1);i<=stopi;i++) {  
        LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i,j,MATCH65) = NEGI;
   
        LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i,j,MATCH75) = NEGI;
   
        LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i,j,MATCH85) = NEGI;
   
        LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i,j,MATCH95) = NEGI;
   
        LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i,j,UNMATCHED_QUERY) = NEGI;
   
        LocalDnaMatchBlock_HIDDEN_MATRIX(mat,i,j,UNMATCHED_TARGET) = NEGI;
  
        }  
      }  


    return;  
}    


/* Function:  full_dc_LocalDnaMatchBlock(mat,starti,startj,startstate,stopi,stopj,stopstate,out,donej,totalj,dpenv)
 *
 * Descrip:    The main divide-and-conquor routine. Basically, call /PackAln_calculate_small_LocalDnaMatchBlock
 *             Not this function, which is pretty hard core. 
 *             Function is given start/end points (in main matrix) for alignment
 *             It does some checks, decides whether start/end in j is small enough for explicit calc
 *               - if yes, calculates it, reads off into PackAln (out), adds the j distance to donej and returns TRUE
 *               - if no,  uses /do_dc_single_pass_LocalDnaMatchBlock to get mid-point
 *                          saves midpoint, and calls itself to do right portion then left portion
 *             right then left ensures PackAln is added the 'right' way, ie, back-to-front
 *             returns FALSE on any error, with a warning
 *
 *
 * Arg:               mat [UNKN ] Matrix with small memory implementation [LocalDnaMatchBlock *]
 * Arg:            starti [UNKN ] Start position in i [int]
 * Arg:            startj [UNKN ] Start position in j [int]
 * Arg:        startstate [UNKN ] Start position state number [int]
 * Arg:             stopi [UNKN ] Stop position in i [int]
 * Arg:             stopj [UNKN ] Stop position in j [int]
 * Arg:         stopstate [UNKN ] Stop position state number [int]
 * Arg:               out [UNKN ] PackAln structure to put alignment into [PackAln *]
 * Arg:             donej [UNKN ] pointer to a number with the amount of alignment done [int *]
 * Arg:            totalj [UNKN ] total amount of alignment to do (in j coordinates) [int]
 * Arg:             dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean full_dc_LocalDnaMatchBlock(LocalDnaMatchBlock * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out,int * donej,int totalj,DPEnvelope * dpenv) 
{
    int lstarti; 
    int lstartj; 
    int lstate;  


    if( mat->basematrix->type != BASEMATRIX_TYPE_SHADOW) {  
      warn("*Very* bad error! - non shadow matrix type in full_dc_LocalDnaMatchBlock");  
      return FALSE;  
      }  


    if( starti == -1 || startj == -1 || startstate == -1 || stopi == -1 || stopstate == -1)  {  
      warn("In full dc program, passed bad indices, indices passed were %d:%d[%d] to %d:%d[%d]\n",starti,startj,startstate,stopi,stopj,stopstate);   
      return FALSE;  
      }  


    if( stopj - startj < 5)  {  
      log_full_error(REPORT,0,"[%d,%d][%d,%d] Explicit read off",starti,startj,stopi,stopj);/* Build hidden explicit matrix */ 
      calculate_hidden_LocalDnaMatchBlock(mat,starti,startj,startstate,stopi,stopj,stopstate,dpenv);     
      *donej += (stopj - startj);   /* Now read it off into out */ 
      if( read_hidden_LocalDnaMatchBlock(mat,starti,startj,startstate,stopi,stopj,stopstate,out) == FALSE)   {  
        warn("In full dc, at %d:%d,%d:%d got a bad hidden explicit read off... ",starti,startj,stopi,stopj); 
        return FALSE;    
        }  
      return TRUE;   
      }  


/* In actual divide and conquor */ 
    if( do_dc_single_pass_LocalDnaMatchBlock(mat,starti,startj,startstate,stopi,stopj,stopstate,dpenv,(int)(*donej*100)/totalj) == FALSE)    {  
      warn("In divide and conquor for LocalDnaMatchBlock, at bound %d:%d to %d:%d, unable to calculate midpoint. Problem!",starti,startj,stopi,stopj);   
      return FALSE;  
      }  


/* Ok... now we have to call on each side of the matrix */ 
/* We have to retrieve left hand side positions, as they will be vapped by the time we call LHS */ 
    lstarti= LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,0);    
    lstartj= LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,1);    
    lstate = LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,2);    


/* Call on right hand side: this lets us do the correct read off */ 
    if( full_dc_LocalDnaMatchBlock(mat,LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,3),LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,4),LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,5),stopi,stopj,stopstate,out,donej,totalj,dpenv) == FALSE)   {  
/* Warning already issued, simply chained back up to top */ 
      return FALSE;  
      }  
/* Call on left hand side */ 
    if( full_dc_LocalDnaMatchBlock(mat,starti,startj,startstate,lstarti,lstartj,lstate,out,donej,totalj,dpenv) == FALSE) {  
/* Warning already issued, simply chained back up to top */ 
      return FALSE;  
      }  


    return TRUE;     
}    


/* Function:  do_dc_single_pass_LocalDnaMatchBlock(mat,starti,startj,startstate,stopi,stopj,stopstate,dpenv,perc_done)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [LocalDnaMatchBlock *]
 * Arg:            starti [UNKN ] Undocumented argument [int]
 * Arg:            startj [UNKN ] Undocumented argument [int]
 * Arg:        startstate [UNKN ] Undocumented argument [int]
 * Arg:             stopi [UNKN ] Undocumented argument [int]
 * Arg:             stopj [UNKN ] Undocumented argument [int]
 * Arg:         stopstate [UNKN ] Undocumented argument [int]
 * Arg:             dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:         perc_done [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean do_dc_single_pass_LocalDnaMatchBlock(LocalDnaMatchBlock * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv,int perc_done) 
{
    int halfj;   
    halfj = startj + ((stopj - startj)/2);   


    init_dc_LocalDnaMatchBlock(mat); 


    LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,starti,startj,startstate) = 0;   
    run_up_dc_LocalDnaMatchBlock(mat,starti,stopi,startj,halfj-1,dpenv,perc_done);   
    push_dc_at_merge_LocalDnaMatchBlock(mat,starti,stopi,halfj,&halfj,dpenv);    
    follow_on_dc_LocalDnaMatchBlock(mat,starti,stopi,halfj,stopj,dpenv,perc_done);   
    return TRUE; 
}    


/* Function:  push_dc_at_merge_LocalDnaMatchBlock(mat,starti,stopi,startj,stopj,dpenv)
 *
 * Descrip: No Description
 *
 * Arg:           mat [UNKN ] Undocumented argument [LocalDnaMatchBlock *]
 * Arg:        starti [UNKN ] Undocumented argument [int]
 * Arg:         stopi [UNKN ] Undocumented argument [int]
 * Arg:        startj [UNKN ] Undocumented argument [int]
 * Arg:         stopj [UNKN ] Undocumented argument [int *]
 * Arg:         dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 */
void push_dc_at_merge_LocalDnaMatchBlock(LocalDnaMatchBlock * mat,int starti,int stopi,int startj,int * stopj,DPEnvelope * dpenv) 
{
    register int i;  
    register int j;  
    register int k;  
    register int count;  
    register int mergej;/* Sources below this j will be stamped by triples */ 
    register int score;  
    register int temp;   


    mergej = startj -1;  
    for(count=0,j=startj;count<1;count++,j++)    {  
      for(i=starti;i<=stopi;i++) {  
        if( dpenv != NULL && is_in_DPEnvelope(dpenv,i,j) == FALSE )  { /*Is not in envelope*/ 
          LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i,j,MATCH65) = NEGI;   
          LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH65,0) = (-100);    
          LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH65,1) = (-100);    
          LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i,j,MATCH75) = NEGI;   
          LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH75,0) = (-100);    
          LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH75,1) = (-100);    
          LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i,j,MATCH85) = NEGI;   
          LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH85,0) = (-100);    
          LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH85,1) = (-100);    
          LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i,j,MATCH95) = NEGI;   
          LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH95,0) = (-100);    
          LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH95,1) = (-100);    
          LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i,j,UNMATCHED_QUERY) = NEGI;   
          LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED_QUERY,0) = (-100);    
          LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED_QUERY,1) = (-100);    
          LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i,j,UNMATCHED_TARGET) = NEGI;  
          LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED_TARGET,0) = (-100);   
          LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED_TARGET,1) = (-100);   
          continue;  
          } /* end of Is not in envelope */ 


        /* For state MATCH65, pushing when j - offj <= mergej */ 
        score = LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i-1,j-1,MATCH65) + (mat->lchs->comp65->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->s);   
        if( j - 1 <= mergej) {  
          LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH65,0) = i-1;   
          LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH65,1) = j-1;   
          LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH65,2) = MATCH65;   
          LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH65,3) = i; 
          LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH65,4) = j; 
          LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH65,5) = MATCH65;   
          }  
        else {  
          for(k=0;k<7;k++)   
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH65,k) = LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,MATCH65,k);   
          }  


        temp = LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i-0,j-1,MATCH65) + mat->lchs->g;  
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH65,0) = i-0; 
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH65,1) = j-1; 
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH65,2) = MATCH65; 
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH65,3) = i;   
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH65,4) = j;   
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH65,5) = MATCH65; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH65,k) = LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,MATCH65,k); 
            }  
          }  


        temp = LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i-1,j-0,MATCH65) + mat->lchs->g;  
        if( temp > score)    {  
          score = temp;  


          if( j - 0 <= mergej)   {  
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH65,0) = i-1; 
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH65,1) = j-0; 
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH65,2) = MATCH65; 
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH65,3) = i;   
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH65,4) = j;   
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH65,5) = MATCH65; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH65,k) = LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 0,MATCH65,k); 
            }  
          }  


        temp = LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i-1,j-1,UNMATCHED_TARGET) + (mat->lchs->comp65->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->v);   
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH65,0) = i-1; 
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH65,1) = j-1; 
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH65,2) = UNMATCHED_TARGET;    
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH65,3) = i;   
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH65,4) = j;   
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH65,5) = MATCH65; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH65,k) = LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,UNMATCHED_TARGET,k);    
            }  
          }  
        /* Add any movement independant score */ 
        LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i,j,MATCH65) = score;    
        /* Finished with state MATCH65 */ 


        /* For state MATCH75, pushing when j - offj <= mergej */ 
        score = LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i-1,j-1,MATCH75) + (mat->lchs->comp75->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->s);   
        if( j - 1 <= mergej) {  
          LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH75,0) = i-1;   
          LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH75,1) = j-1;   
          LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH75,2) = MATCH75;   
          LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH75,3) = i; 
          LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH75,4) = j; 
          LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH75,5) = MATCH75;   
          }  
        else {  
          for(k=0;k<7;k++)   
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH75,k) = LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,MATCH75,k);   
          }  


        temp = LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i-0,j-1,MATCH75) + mat->lchs->g;  
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH75,0) = i-0; 
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH75,1) = j-1; 
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH75,2) = MATCH75; 
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH75,3) = i;   
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH75,4) = j;   
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH75,5) = MATCH75; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH75,k) = LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,MATCH75,k); 
            }  
          }  


        temp = LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i-1,j-0,MATCH75) + mat->lchs->g;  
        if( temp > score)    {  
          score = temp;  


          if( j - 0 <= mergej)   {  
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH75,0) = i-1; 
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH75,1) = j-0; 
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH75,2) = MATCH75; 
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH75,3) = i;   
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH75,4) = j;   
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH75,5) = MATCH75; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH75,k) = LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 0,MATCH75,k); 
            }  
          }  


        temp = LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i-1,j-1,UNMATCHED_TARGET) + (mat->lchs->comp75->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->v);   
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH75,0) = i-1; 
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH75,1) = j-1; 
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH75,2) = UNMATCHED_TARGET;    
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH75,3) = i;   
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH75,4) = j;   
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH75,5) = MATCH75; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH75,k) = LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,UNMATCHED_TARGET,k);    
            }  
          }  
        /* Add any movement independant score */ 
        LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i,j,MATCH75) = score;    
        /* Finished with state MATCH75 */ 


        /* For state MATCH85, pushing when j - offj <= mergej */ 
        score = LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i-1,j-1,MATCH85) + (mat->lchs->comp85->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->s);   
        if( j - 1 <= mergej) {  
          LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH85,0) = i-1;   
          LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH85,1) = j-1;   
          LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH85,2) = MATCH85;   
          LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH85,3) = i; 
          LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH85,4) = j; 
          LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH85,5) = MATCH85;   
          }  
        else {  
          for(k=0;k<7;k++)   
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH85,k) = LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,MATCH85,k);   
          }  


        temp = LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i-0,j-1,MATCH85) + mat->lchs->g;  
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH85,0) = i-0; 
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH85,1) = j-1; 
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH85,2) = MATCH85; 
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH85,3) = i;   
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH85,4) = j;   
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH85,5) = MATCH85; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH85,k) = LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,MATCH85,k); 
            }  
          }  


        temp = LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i-1,j-0,MATCH85) + mat->lchs->g;  
        if( temp > score)    {  
          score = temp;  


          if( j - 0 <= mergej)   {  
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH85,0) = i-1; 
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH85,1) = j-0; 
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH85,2) = MATCH85; 
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH85,3) = i;   
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH85,4) = j;   
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH85,5) = MATCH85; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH85,k) = LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 0,MATCH85,k); 
            }  
          }  


        temp = LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i-1,j-1,UNMATCHED_TARGET) + (mat->lchs->comp85->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->v);   
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH85,0) = i-1; 
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH85,1) = j-1; 
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH85,2) = UNMATCHED_TARGET;    
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH85,3) = i;   
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH85,4) = j;   
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH85,5) = MATCH85; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH85,k) = LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,UNMATCHED_TARGET,k);    
            }  
          }  
        /* Add any movement independant score */ 
        LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i,j,MATCH85) = score;    
        /* Finished with state MATCH85 */ 


        /* For state MATCH95, pushing when j - offj <= mergej */ 
        score = LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i-1,j-1,MATCH95) + (mat->lchs->comp95->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->s);   
        if( j - 1 <= mergej) {  
          LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH95,0) = i-1;   
          LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH95,1) = j-1;   
          LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH95,2) = MATCH95;   
          LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH95,3) = i; 
          LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH95,4) = j; 
          LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH95,5) = MATCH95;   
          }  
        else {  
          for(k=0;k<7;k++)   
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH95,k) = LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,MATCH95,k);   
          }  


        temp = LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i-0,j-1,MATCH95) + mat->lchs->g;  
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH95,0) = i-0; 
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH95,1) = j-1; 
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH95,2) = MATCH95; 
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH95,3) = i;   
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH95,4) = j;   
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH95,5) = MATCH95; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH95,k) = LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,MATCH95,k); 
            }  
          }  


        temp = LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i-1,j-0,MATCH95) + mat->lchs->g;  
        if( temp > score)    {  
          score = temp;  


          if( j - 0 <= mergej)   {  
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH95,0) = i-1; 
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH95,1) = j-0; 
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH95,2) = MATCH95; 
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH95,3) = i;   
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH95,4) = j;   
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH95,5) = MATCH95; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH95,k) = LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 0,MATCH95,k); 
            }  
          }  


        temp = LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i-1,j-1,UNMATCHED_TARGET) + (mat->lchs->comp95->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->v);   
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH95,0) = i-1; 
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH95,1) = j-1; 
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH95,2) = UNMATCHED_TARGET;    
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH95,3) = i;   
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH95,4) = j;   
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH95,5) = MATCH95; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH95,k) = LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,UNMATCHED_TARGET,k);    
            }  
          }  
        /* Add any movement independant score */ 
        LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i,j,MATCH95) = score;    
        /* Finished with state MATCH95 */ 


        /* For state UNMATCHED_QUERY, pushing when j - offj <= mergej */ 
        score = LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i-1,j-0,MATCH65) + mat->lchs->b;     
        if( j - 0 <= mergej) {  
          LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED_QUERY,0) = i-1;   
          LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED_QUERY,1) = j-0;   
          LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED_QUERY,2) = MATCH65;   
          LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED_QUERY,3) = i; 
          LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED_QUERY,4) = j; 
          LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED_QUERY,5) = UNMATCHED_QUERY;   
          }  
        else {  
          for(k=0;k<7;k++)   
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED_QUERY,k) = LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 0,MATCH65,k);   
          }  


        temp = LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i-1,j-0,MATCH75) + mat->lchs->b;  
        if( temp > score)    {  
          score = temp;  


          if( j - 0 <= mergej)   {  
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED_QUERY,0) = i-1; 
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED_QUERY,1) = j-0; 
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED_QUERY,2) = MATCH75; 
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED_QUERY,3) = i;   
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED_QUERY,4) = j;   
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED_QUERY,5) = UNMATCHED_QUERY; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED_QUERY,k) = LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 0,MATCH75,k); 
            }  
          }  


        temp = LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i-1,j-0,MATCH85) + mat->lchs->b;  
        if( temp > score)    {  
          score = temp;  


          if( j - 0 <= mergej)   {  
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED_QUERY,0) = i-1; 
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED_QUERY,1) = j-0; 
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED_QUERY,2) = MATCH85; 
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED_QUERY,3) = i;   
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED_QUERY,4) = j;   
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED_QUERY,5) = UNMATCHED_QUERY; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED_QUERY,k) = LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 0,MATCH85,k); 
            }  
          }  


        temp = LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i-1,j-0,MATCH95) + mat->lchs->b;  
        if( temp > score)    {  
          score = temp;  


          if( j - 0 <= mergej)   {  
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED_QUERY,0) = i-1; 
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED_QUERY,1) = j-0; 
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED_QUERY,2) = MATCH95; 
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED_QUERY,3) = i;   
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED_QUERY,4) = j;   
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED_QUERY,5) = UNMATCHED_QUERY; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED_QUERY,k) = LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 0,MATCH95,k); 
            }  
          }  


        temp = LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i-1,j-0,UNMATCHED_QUERY) + mat->lchs->u;  
        if( temp > score)    {  
          score = temp;  


          if( j - 0 <= mergej)   {  
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED_QUERY,0) = i-1; 
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED_QUERY,1) = j-0; 
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED_QUERY,2) = UNMATCHED_QUERY; 
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED_QUERY,3) = i;   
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED_QUERY,4) = j;   
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED_QUERY,5) = UNMATCHED_QUERY; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED_QUERY,k) = LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 0,UNMATCHED_QUERY,k); 
            }  
          }  
        /* Add any movement independant score */ 
        LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i,j,UNMATCHED_QUERY) = score;    
        /* Finished with state UNMATCHED_QUERY */ 


        /* For state UNMATCHED_TARGET, pushing when j - offj <= mergej */ 
        score = LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i-0,j-1,UNMATCHED_QUERY) + mat->lchs->v;     
        if( j - 1 <= mergej) {  
          LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED_TARGET,0) = i-0;  
          LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED_TARGET,1) = j-1;  
          LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED_TARGET,2) = UNMATCHED_QUERY;  
          LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED_TARGET,3) = i;    
          LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED_TARGET,4) = j;    
          LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED_TARGET,5) = UNMATCHED_TARGET; 
          }  
        else {  
          for(k=0;k<7;k++)   
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED_TARGET,k) = LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,UNMATCHED_QUERY,k);  
          }  


        temp = LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i-0,j-1,UNMATCHED_TARGET) + mat->lchs->u;     
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED_TARGET,0) = i-0;    
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED_TARGET,1) = j-1;    
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED_TARGET,2) = UNMATCHED_TARGET;   
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED_TARGET,3) = i;  
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED_TARGET,4) = j;  
            LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED_TARGET,5) = UNMATCHED_TARGET;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED_TARGET,k) = LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,UNMATCHED_TARGET,k);   
            }  
          }  
        /* Add any movement independant score */ 
        LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i,j,UNMATCHED_TARGET) = score;   
        /* Finished with state UNMATCHED_TARGET */ 
        }  
      }  
    /* Put back j into * stop j so that calling function gets it correct */ 
    if( stopj == NULL)   
      warn("Bad news... NULL stopj pointer in push dc function. This means that calling function does not know how many cells I have done!");    
    else 
      *stopj = j;    


    return;  
}    


/* Function:  follow_on_dc_LocalDnaMatchBlock(mat,starti,stopi,startj,stopj,dpenv,perc_done)
 *
 * Descrip: No Description
 *
 * Arg:              mat [UNKN ] Undocumented argument [LocalDnaMatchBlock *]
 * Arg:           starti [UNKN ] Undocumented argument [int]
 * Arg:            stopi [UNKN ] Undocumented argument [int]
 * Arg:           startj [UNKN ] Undocumented argument [int]
 * Arg:            stopj [UNKN ] Undocumented argument [int]
 * Arg:            dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:        perc_done [UNKN ] Undocumented argument [int]
 *
 */
void follow_on_dc_LocalDnaMatchBlock(LocalDnaMatchBlock * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done) 
{
    int i;   
    int j;   
    int k;   
    int score;   
    int temp;    
    int localshadow[7];  
    long int total;  
    long int num;    


    total = (stopi - starti+1) * (stopj - startj+1); 
    num = 0;     


    for(j=startj;j<=stopj;j++)   { /*for each valid j column*/ 
      for(i=starti;i<=stopi;i++) { /*this is strip*/ 
        num++;   
        if( dpenv != NULL && is_in_DPEnvelope(dpenv,i,j) == FALSE )  { /*Is not in envelope*/ 
          LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i,j,MATCH65) = NEGI;   
          LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i,j,MATCH75) = NEGI;   
          LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i,j,MATCH85) = NEGI;   
          LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i,j,MATCH95) = NEGI;   
          LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i,j,UNMATCHED_QUERY) = NEGI;   
          LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i,j,UNMATCHED_TARGET) = NEGI;  
          continue;  
          } /* end of Is not in envelope */ 
        if( num % 1000 == 0 )    
          log_full_error(REPORT,0,"[%d%%%% done]After  mid-j %5d Cells done %d%%%%",perc_done,startj,(num*100)/total);   


        /* For state MATCH65 */ 
        /* setting first movement to score */ 
        score = LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i-1,j-1,MATCH65) + (mat->lchs->comp65->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->s);   
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,MATCH65,k);    
        /* From state MATCH65 to state MATCH65 */ 
        temp = LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i-0,j-1,MATCH65) + mat->lchs->g;  
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,MATCH65,k);  
          }  
        /* From state MATCH65 to state MATCH65 */ 
        temp = LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i-1,j-0,MATCH65) + mat->lchs->g;  
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 0,MATCH65,k);  
          }  
        /* From state UNMATCHED_TARGET to state MATCH65 */ 
        temp = LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i-1,j-1,UNMATCHED_TARGET) + (mat->lchs->comp65->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->v);   
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,UNMATCHED_TARGET,k); 
          }  


        /* Ok - finished max calculation for MATCH65 */ 
        /* Add any movement independant score and put away */ 
         LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i,j,MATCH65) = score;   
        for(k=0;k<7;k++) 
          LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH65,k) = localshadow[k];    
        /* Now figure out if any specials need this score */ 
        /* Finished calculating state MATCH65 */ 


        /* For state MATCH75 */ 
        /* setting first movement to score */ 
        score = LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i-1,j-1,MATCH75) + (mat->lchs->comp75->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->s);   
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,MATCH75,k);    
        /* From state MATCH75 to state MATCH75 */ 
        temp = LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i-0,j-1,MATCH75) + mat->lchs->g;  
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,MATCH75,k);  
          }  
        /* From state MATCH75 to state MATCH75 */ 
        temp = LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i-1,j-0,MATCH75) + mat->lchs->g;  
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 0,MATCH75,k);  
          }  
        /* From state UNMATCHED_TARGET to state MATCH75 */ 
        temp = LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i-1,j-1,UNMATCHED_TARGET) + (mat->lchs->comp75->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->v);   
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,UNMATCHED_TARGET,k); 
          }  


        /* Ok - finished max calculation for MATCH75 */ 
        /* Add any movement independant score and put away */ 
         LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i,j,MATCH75) = score;   
        for(k=0;k<7;k++) 
          LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH75,k) = localshadow[k];    
        /* Now figure out if any specials need this score */ 
        /* Finished calculating state MATCH75 */ 


        /* For state MATCH85 */ 
        /* setting first movement to score */ 
        score = LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i-1,j-1,MATCH85) + (mat->lchs->comp85->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->s);   
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,MATCH85,k);    
        /* From state MATCH85 to state MATCH85 */ 
        temp = LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i-0,j-1,MATCH85) + mat->lchs->g;  
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,MATCH85,k);  
          }  
        /* From state MATCH85 to state MATCH85 */ 
        temp = LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i-1,j-0,MATCH85) + mat->lchs->g;  
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 0,MATCH85,k);  
          }  
        /* From state UNMATCHED_TARGET to state MATCH85 */ 
        temp = LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i-1,j-1,UNMATCHED_TARGET) + (mat->lchs->comp85->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->v);   
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,UNMATCHED_TARGET,k); 
          }  


        /* Ok - finished max calculation for MATCH85 */ 
        /* Add any movement independant score and put away */ 
         LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i,j,MATCH85) = score;   
        for(k=0;k<7;k++) 
          LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH85,k) = localshadow[k];    
        /* Now figure out if any specials need this score */ 
        /* Finished calculating state MATCH85 */ 


        /* For state MATCH95 */ 
        /* setting first movement to score */ 
        score = LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i-1,j-1,MATCH95) + (mat->lchs->comp95->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->s);   
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,MATCH95,k);    
        /* From state MATCH95 to state MATCH95 */ 
        temp = LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i-0,j-1,MATCH95) + mat->lchs->g;  
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,MATCH95,k);  
          }  
        /* From state MATCH95 to state MATCH95 */ 
        temp = LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i-1,j-0,MATCH95) + mat->lchs->g;  
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 0,MATCH95,k);  
          }  
        /* From state UNMATCHED_TARGET to state MATCH95 */ 
        temp = LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i-1,j-1,UNMATCHED_TARGET) + (mat->lchs->comp95->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->v);   
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,UNMATCHED_TARGET,k); 
          }  


        /* Ok - finished max calculation for MATCH95 */ 
        /* Add any movement independant score and put away */ 
         LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i,j,MATCH95) = score;   
        for(k=0;k<7;k++) 
          LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH95,k) = localshadow[k];    
        /* Now figure out if any specials need this score */ 
        /* Finished calculating state MATCH95 */ 


        /* For state UNMATCHED_QUERY */ 
        /* setting first movement to score */ 
        score = LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i-1,j-0,MATCH65) + mat->lchs->b;     
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 0,MATCH65,k);    
        /* From state MATCH75 to state UNMATCHED_QUERY */ 
        temp = LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i-1,j-0,MATCH75) + mat->lchs->b;  
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 0,MATCH75,k);  
          }  
        /* From state MATCH85 to state UNMATCHED_QUERY */ 
        temp = LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i-1,j-0,MATCH85) + mat->lchs->b;  
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 0,MATCH85,k);  
          }  
        /* From state MATCH95 to state UNMATCHED_QUERY */ 
        temp = LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i-1,j-0,MATCH95) + mat->lchs->b;  
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 0,MATCH95,k);  
          }  
        /* From state UNMATCHED_QUERY to state UNMATCHED_QUERY */ 
        temp = LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i-1,j-0,UNMATCHED_QUERY) + mat->lchs->u;  
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 0,UNMATCHED_QUERY,k);  
          }  


        /* Ok - finished max calculation for UNMATCHED_QUERY */ 
        /* Add any movement independant score and put away */ 
         LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i,j,UNMATCHED_QUERY) = score;   
        for(k=0;k<7;k++) 
          LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED_QUERY,k) = localshadow[k];    
        /* Now figure out if any specials need this score */ 
        /* Finished calculating state UNMATCHED_QUERY */ 


        /* For state UNMATCHED_TARGET */ 
        /* setting first movement to score */ 
        score = LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i-0,j-1,UNMATCHED_QUERY) + mat->lchs->v;     
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,UNMATCHED_QUERY,k);    
        /* From state UNMATCHED_TARGET to state UNMATCHED_TARGET */ 
        temp = LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i-0,j-1,UNMATCHED_TARGET) + mat->lchs->u;     
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,UNMATCHED_TARGET,k); 
          }  


        /* Ok - finished max calculation for UNMATCHED_TARGET */ 
        /* Add any movement independant score and put away */ 
         LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i,j,UNMATCHED_TARGET) = score;  
        for(k=0;k<7;k++) 
          LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED_TARGET,k) = localshadow[k];   
        /* Now figure out if any specials need this score */ 
        /* Finished calculating state UNMATCHED_TARGET */ 
        } /* end of this is strip */ 
      } /* end of for each valid j column */ 


/* Function:  run_up_dc_LocalDnaMatchBlock(mat,starti,stopi,startj,stopj,dpenv,perc_done)
 *
 * Descrip: No Description
 *
 * Arg:              mat [UNKN ] Undocumented argument [LocalDnaMatchBlock *]
 * Arg:           starti [UNKN ] Undocumented argument [int]
 * Arg:            stopi [UNKN ] Undocumented argument [int]
 * Arg:           startj [UNKN ] Undocumented argument [int]
 * Arg:            stopj [UNKN ] Undocumented argument [int]
 * Arg:            dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:        perc_done [UNKN ] Undocumented argument [int]
 *
 */
}    
void run_up_dc_LocalDnaMatchBlock(LocalDnaMatchBlock * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done) 
{
    register int i;  
    register int j;  
    register int score;  
    register int temp;   
    long int total;  
    long int num;    


    total = (stopi - starti+1) * (stopj - startj+1); 
    if( total <= 0 ) 
      total = 1; 
    num = 0;     


    for(j=startj;j<=stopj;j++)   { /*for each valid j column*/ 
      for(i=starti;i<=stopi;i++) { /*this is strip*/ 
        if( j == startj && i == starti)  
          continue;  
        num++;   
        if( dpenv != NULL && is_in_DPEnvelope(dpenv,i,j) == FALSE )  { /*Is not in envelope*/ 
          LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i,j,MATCH65) = NEGI;   
          LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i,j,MATCH75) = NEGI;   
          LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i,j,MATCH85) = NEGI;   
          LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i,j,MATCH95) = NEGI;   
          LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i,j,UNMATCHED_QUERY) = NEGI;   
          LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i,j,UNMATCHED_TARGET) = NEGI;  
          continue;  
          } /* end of Is not in envelope */ 
        if( num % 1000 == 0 )    
          log_full_error(REPORT,0,"[%d%%%% done]Before mid-j %5d Cells done %d%%%%",perc_done,stopj,(num*100)/total);    


        /* For state MATCH65 */ 
        /* setting first movement to score */ 
        score = LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i-1,j-1,MATCH65) + (mat->lchs->comp65->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->s);   
        /* From state MATCH65 to state MATCH65 */ 
        temp = LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i-0,j-1,MATCH65) + mat->lchs->g;  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH65 to state MATCH65 */ 
        temp = LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i-1,j-0,MATCH65) + mat->lchs->g;  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state UNMATCHED_TARGET to state MATCH65 */ 
        temp = LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i-1,j-1,UNMATCHED_TARGET) + (mat->lchs->comp65->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->v);   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for MATCH65 */ 
        /* Add any movement independant score and put away */ 
         LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i,j,MATCH65) = score;   
        /* Finished calculating state MATCH65 */ 


        /* For state MATCH75 */ 
        /* setting first movement to score */ 
        score = LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i-1,j-1,MATCH75) + (mat->lchs->comp75->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->s);   
        /* From state MATCH75 to state MATCH75 */ 
        temp = LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i-0,j-1,MATCH75) + mat->lchs->g;  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH75 to state MATCH75 */ 
        temp = LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i-1,j-0,MATCH75) + mat->lchs->g;  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state UNMATCHED_TARGET to state MATCH75 */ 
        temp = LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i-1,j-1,UNMATCHED_TARGET) + (mat->lchs->comp75->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->v);   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for MATCH75 */ 
        /* Add any movement independant score and put away */ 
         LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i,j,MATCH75) = score;   
        /* Finished calculating state MATCH75 */ 


        /* For state MATCH85 */ 
        /* setting first movement to score */ 
        score = LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i-1,j-1,MATCH85) + (mat->lchs->comp85->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->s);   
        /* From state MATCH85 to state MATCH85 */ 
        temp = LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i-0,j-1,MATCH85) + mat->lchs->g;  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH85 to state MATCH85 */ 
        temp = LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i-1,j-0,MATCH85) + mat->lchs->g;  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state UNMATCHED_TARGET to state MATCH85 */ 
        temp = LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i-1,j-1,UNMATCHED_TARGET) + (mat->lchs->comp85->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->v);   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for MATCH85 */ 
        /* Add any movement independant score and put away */ 
         LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i,j,MATCH85) = score;   
        /* Finished calculating state MATCH85 */ 


        /* For state MATCH95 */ 
        /* setting first movement to score */ 
        score = LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i-1,j-1,MATCH95) + (mat->lchs->comp95->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->s);   
        /* From state MATCH95 to state MATCH95 */ 
        temp = LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i-0,j-1,MATCH95) + mat->lchs->g;  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH95 to state MATCH95 */ 
        temp = LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i-1,j-0,MATCH95) + mat->lchs->g;  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state UNMATCHED_TARGET to state MATCH95 */ 
        temp = LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i-1,j-1,UNMATCHED_TARGET) + (mat->lchs->comp95->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->v);   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for MATCH95 */ 
        /* Add any movement independant score and put away */ 
         LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i,j,MATCH95) = score;   
        /* Finished calculating state MATCH95 */ 


        /* For state UNMATCHED_QUERY */ 
        /* setting first movement to score */ 
        score = LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i-1,j-0,MATCH65) + mat->lchs->b;     
        /* From state MATCH75 to state UNMATCHED_QUERY */ 
        temp = LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i-1,j-0,MATCH75) + mat->lchs->b;  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH85 to state UNMATCHED_QUERY */ 
        temp = LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i-1,j-0,MATCH85) + mat->lchs->b;  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH95 to state UNMATCHED_QUERY */ 
        temp = LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i-1,j-0,MATCH95) + mat->lchs->b;  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state UNMATCHED_QUERY to state UNMATCHED_QUERY */ 
        temp = LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i-1,j-0,UNMATCHED_QUERY) + mat->lchs->u;  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for UNMATCHED_QUERY */ 
        /* Add any movement independant score and put away */ 
         LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i,j,UNMATCHED_QUERY) = score;   
        /* Finished calculating state UNMATCHED_QUERY */ 


        /* For state UNMATCHED_TARGET */ 
        /* setting first movement to score */ 
        score = LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i-0,j-1,UNMATCHED_QUERY) + mat->lchs->v;     
        /* From state UNMATCHED_TARGET to state UNMATCHED_TARGET */ 
        temp = LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i-0,j-1,UNMATCHED_TARGET) + mat->lchs->u;     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for UNMATCHED_TARGET */ 
        /* Add any movement independant score and put away */ 
         LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i,j,UNMATCHED_TARGET) = score;  
        /* Finished calculating state UNMATCHED_TARGET */ 
        } /* end of this is strip */ 
      } /* end of for each valid j column */ 


/* Function:  init_dc_LocalDnaMatchBlock(mat)
 *
 * Descrip: No Description
 *
 * Arg:        mat [UNKN ] Undocumented argument [LocalDnaMatchBlock *]
 *
 */
}    
void init_dc_LocalDnaMatchBlock(LocalDnaMatchBlock * mat) 
{
    register int i;  
    register int j;  
    register int k;  


    for(j=0;j<3;j++) {  
      for(i=(-1);i<mat->query->seq->len;i++) {  
        LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i,j,MATCH65) = NEGI; 
        LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i,j,MATCH75) = NEGI; 
        LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i,j,MATCH85) = NEGI; 
        LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i,j,MATCH95) = NEGI; 
        LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i,j,UNMATCHED_QUERY) = NEGI; 
        LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i,j,UNMATCHED_TARGET) = NEGI;    
        for(k=0;k<7;k++) {  
          LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH65,k) = (-1);  
          LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH75,k) = (-1);  
          LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH85,k) = (-1);  
          LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH95,k) = (-1);  
          LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED_QUERY,k) = (-1);  
          LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED_TARGET,k) = (-1); 
          }  
        }  
      }  


    return;  
}    


/* Function:  start_end_find_end_LocalDnaMatchBlock(mat,endj)
 *
 * Descrip:    First function used to find end of the best path in the special state !end
 *
 *
 * Arg:         mat [UNKN ] Matrix in small mode [LocalDnaMatchBlock *]
 * Arg:        endj [WRITE] position of end in j (meaningless in i) [int *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int start_end_find_end_LocalDnaMatchBlock(LocalDnaMatchBlock * mat,int * endj) 
{
    register int j;  
    register int max;    
    register int maxj;   


    max = LocalDnaMatchBlock_DC_SHADOW_SPECIAL(mat,0,mat->target->seq->len-1,END);   
    maxj = mat->target->seq->len-1;  
    for(j= mat->target->seq->len-2 ;j >= 0 ;j--) {  
      if( LocalDnaMatchBlock_DC_SHADOW_SPECIAL(mat,0,j,END) > max )  {  
        max = LocalDnaMatchBlock_DC_SHADOW_SPECIAL(mat,0,j,END); 
        maxj = j;    
        }  
      }  


    if( endj != NULL)    
      *endj = maxj;  


    return max;  
}    


/* Function:  dc_optimised_start_end_calc_LocalDnaMatchBlock(*mat,dpenv)
 *
 * Descrip:    Calculates special strip, leaving start/end/score points in shadow matrix
 *             Works off specially laid out memory from steve searle
 *
 *
 * Arg:         *mat [UNKN ] Undocumented argument [LocalDnaMatchBlock]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean dc_optimised_start_end_calc_LocalDnaMatchBlock(LocalDnaMatchBlock *mat,DPEnvelope * dpenv) 
{
    int i;   
    int j;   
    int k;   
    int score;   
    int temp;    
    int leni;    
    int lenj;    
    int localshadow[7];  
    long int total;  
    long int num=0;  
    int * score_pointers;    
    int * shadow_pointers;   
    int * localsp;   
    leni = mat->query->seq->len; 
    lenj = mat->target->seq->len;    
    total = leni * lenj; 


    score_pointers = (int *) calloc (1 * (leni + 1) * 6,sizeof(int));    
    shadow_pointers = (int *) calloc (1 * (leni + 1) * 6 * 8,sizeof(int));   


    for(j=0;j<lenj;j++)  { /*for each j strip*/ 
      for(i=0;i<leni;i++)    { /*for each i position in strip*/ 
        num++;   
        if( dpenv != NULL && is_in_DPEnvelope(dpenv,i,j) == FALSE )  { /*Is not in envelope*/ 
          LocalDnaMatchBlock_DC_OPT_SHADOW_MATRIX(mat,i,j,MATCH65) = NEGI;   
          LocalDnaMatchBlock_DC_OPT_SHADOW_MATRIX(mat,i,j,MATCH75) = NEGI;   
          LocalDnaMatchBlock_DC_OPT_SHADOW_MATRIX(mat,i,j,MATCH85) = NEGI;   
          LocalDnaMatchBlock_DC_OPT_SHADOW_MATRIX(mat,i,j,MATCH95) = NEGI;   
          LocalDnaMatchBlock_DC_OPT_SHADOW_MATRIX(mat,i,j,UNMATCHED_QUERY) = NEGI;   
          LocalDnaMatchBlock_DC_OPT_SHADOW_MATRIX(mat,i,j,UNMATCHED_TARGET) = NEGI;  
          continue;  
          } /* end of Is not in envelope */ 
        if( num%1000 == 0)   
          log_full_error(REPORT,0,"%6d Cells done [%2d%%%%]",num,num*100/total); 




        /* For state MATCH65 */ 
        /* setting first movement to score */ 
        score = LocalDnaMatchBlock_DC_OPT_SHADOW_MATRIX(mat,i-1,j-1,MATCH65) + (mat->lchs->comp65->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->s) + (0);     
        /* assign local shadown pointer */ 
        localsp = &(LocalDnaMatchBlock_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 1,MATCH65,0));  
        /* From state MATCH65 to state MATCH65 */ 
        temp = LocalDnaMatchBlock_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,MATCH65) + mat->lchs->g +(0);     
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(LocalDnaMatchBlock_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,MATCH65,0));    
          }  
        /* From state MATCH65 to state MATCH65 */ 
        temp = LocalDnaMatchBlock_DC_OPT_SHADOW_MATRIX(mat,i-1,j-0,MATCH65) + mat->lchs->g +(0);     
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(LocalDnaMatchBlock_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 0,MATCH65,0));    
          }  
        /* From state START to state MATCH65 */ 
        temp = LocalDnaMatchBlock_DC_OPT_SHADOW_SPECIAL(mat,i-1,j-1,START) + (mat->lchs->comp65->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->v) + (0);   
        if( temp  > score )  {  
          score = temp;  
          /* This state [START] is a special for MATCH65... push top shadow pointers here */ 
          localshadow[0]= i; 
          localshadow[1]= j; 
          localshadow[2]= MATCH65;   
          localshadow[3]= (-1);  
          localshadow[4]= (-1);  
          localshadow[5]= (-1);  
          localshadow[6]= score; 
          localsp = localshadow; 
          }  
        /* From state UNMATCHED_TARGET to state MATCH65 */ 
        temp = LocalDnaMatchBlock_DC_OPT_SHADOW_MATRIX(mat,i-1,j-1,UNMATCHED_TARGET) + (mat->lchs->comp65->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->v) +(0);  
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(LocalDnaMatchBlock_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 1,UNMATCHED_TARGET,0));   
          }  


        /* Ok - finished max calculation for MATCH65 */ 
        /* Add any movement independant score and put away */ 
        /* Actually, already done inside scores */ 
         LocalDnaMatchBlock_DC_OPT_SHADOW_MATRIX(mat,i,j,MATCH65) = score;   
        for(k=0;k<7;k++) 
          LocalDnaMatchBlock_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,MATCH65,k) = localsp[k];    
        /* Now figure out if any specials need this score */ 


        /* state MATCH65 is a source for special END */ 
        temp = score + (0) + (0) ;   
        if( temp > LocalDnaMatchBlock_DC_OPT_SHADOW_SPECIAL(mat,i,j,END) )   {  
          LocalDnaMatchBlock_DC_OPT_SHADOW_SPECIAL(mat,i,j,END) = temp;  
          /* Have to push only bottem half of system here */ 
          for(k=0;k<3;k++)   
            LocalDnaMatchBlock_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,k) = LocalDnaMatchBlock_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,MATCH65,k);  
          LocalDnaMatchBlock_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,6) = LocalDnaMatchBlock_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,MATCH65,6);    
          LocalDnaMatchBlock_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,3) = i;    
          LocalDnaMatchBlock_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,4) = j;    
          LocalDnaMatchBlock_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,5) = MATCH65;  
          }  




        /* Finished calculating state MATCH65 */ 


        /* For state MATCH75 */ 
        /* setting first movement to score */ 
        score = LocalDnaMatchBlock_DC_OPT_SHADOW_MATRIX(mat,i-1,j-1,MATCH75) + (mat->lchs->comp75->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->s) + (0);     
        /* assign local shadown pointer */ 
        localsp = &(LocalDnaMatchBlock_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 1,MATCH75,0));  
        /* From state MATCH75 to state MATCH75 */ 
        temp = LocalDnaMatchBlock_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,MATCH75) + mat->lchs->g +(0);     
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(LocalDnaMatchBlock_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,MATCH75,0));    
          }  
        /* From state MATCH75 to state MATCH75 */ 
        temp = LocalDnaMatchBlock_DC_OPT_SHADOW_MATRIX(mat,i-1,j-0,MATCH75) + mat->lchs->g +(0);     
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(LocalDnaMatchBlock_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 0,MATCH75,0));    
          }  
        /* From state START to state MATCH75 */ 
        temp = LocalDnaMatchBlock_DC_OPT_SHADOW_SPECIAL(mat,i-1,j-1,START) + (mat->lchs->comp75->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->v) + (0);   
        if( temp  > score )  {  
          score = temp;  
          /* This state [START] is a special for MATCH75... push top shadow pointers here */ 
          localshadow[0]= i; 
          localshadow[1]= j; 
          localshadow[2]= MATCH75;   
          localshadow[3]= (-1);  
          localshadow[4]= (-1);  
          localshadow[5]= (-1);  
          localshadow[6]= score; 
          localsp = localshadow; 
          }  
        /* From state UNMATCHED_TARGET to state MATCH75 */ 
        temp = LocalDnaMatchBlock_DC_OPT_SHADOW_MATRIX(mat,i-1,j-1,UNMATCHED_TARGET) + (mat->lchs->comp75->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->v) +(0);  
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(LocalDnaMatchBlock_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 1,UNMATCHED_TARGET,0));   
          }  


        /* Ok - finished max calculation for MATCH75 */ 
        /* Add any movement independant score and put away */ 
        /* Actually, already done inside scores */ 
         LocalDnaMatchBlock_DC_OPT_SHADOW_MATRIX(mat,i,j,MATCH75) = score;   
        for(k=0;k<7;k++) 
          LocalDnaMatchBlock_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,MATCH75,k) = localsp[k];    
        /* Now figure out if any specials need this score */ 


        /* state MATCH75 is a source for special END */ 
        temp = score + (0) + (0) ;   
        if( temp > LocalDnaMatchBlock_DC_OPT_SHADOW_SPECIAL(mat,i,j,END) )   {  
          LocalDnaMatchBlock_DC_OPT_SHADOW_SPECIAL(mat,i,j,END) = temp;  
          /* Have to push only bottem half of system here */ 
          for(k=0;k<3;k++)   
            LocalDnaMatchBlock_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,k) = LocalDnaMatchBlock_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,MATCH75,k);  
          LocalDnaMatchBlock_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,6) = LocalDnaMatchBlock_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,MATCH75,6);    
          LocalDnaMatchBlock_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,3) = i;    
          LocalDnaMatchBlock_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,4) = j;    
          LocalDnaMatchBlock_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,5) = MATCH75;  
          }  




        /* Finished calculating state MATCH75 */ 


        /* For state MATCH85 */ 
        /* setting first movement to score */ 
        score = LocalDnaMatchBlock_DC_OPT_SHADOW_MATRIX(mat,i-1,j-1,MATCH85) + (mat->lchs->comp85->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->s) + (0);     
        /* assign local shadown pointer */ 
        localsp = &(LocalDnaMatchBlock_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 1,MATCH85,0));  
        /* From state MATCH85 to state MATCH85 */ 
        temp = LocalDnaMatchBlock_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,MATCH85) + mat->lchs->g +(0);     
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(LocalDnaMatchBlock_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,MATCH85,0));    
          }  
        /* From state MATCH85 to state MATCH85 */ 
        temp = LocalDnaMatchBlock_DC_OPT_SHADOW_MATRIX(mat,i-1,j-0,MATCH85) + mat->lchs->g +(0);     
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(LocalDnaMatchBlock_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 0,MATCH85,0));    
          }  
        /* From state START to state MATCH85 */ 
        temp = LocalDnaMatchBlock_DC_OPT_SHADOW_SPECIAL(mat,i-1,j-1,START) + (mat->lchs->comp85->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->v) + (0);   
        if( temp  > score )  {  
          score = temp;  
          /* This state [START] is a special for MATCH85... push top shadow pointers here */ 
          localshadow[0]= i; 
          localshadow[1]= j; 
          localshadow[2]= MATCH85;   
          localshadow[3]= (-1);  
          localshadow[4]= (-1);  
          localshadow[5]= (-1);  
          localshadow[6]= score; 
          localsp = localshadow; 
          }  
        /* From state UNMATCHED_TARGET to state MATCH85 */ 
        temp = LocalDnaMatchBlock_DC_OPT_SHADOW_MATRIX(mat,i-1,j-1,UNMATCHED_TARGET) + (mat->lchs->comp85->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->v) +(0);  
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(LocalDnaMatchBlock_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 1,UNMATCHED_TARGET,0));   
          }  


        /* Ok - finished max calculation for MATCH85 */ 
        /* Add any movement independant score and put away */ 
        /* Actually, already done inside scores */ 
         LocalDnaMatchBlock_DC_OPT_SHADOW_MATRIX(mat,i,j,MATCH85) = score;   
        for(k=0;k<7;k++) 
          LocalDnaMatchBlock_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,MATCH85,k) = localsp[k];    
        /* Now figure out if any specials need this score */ 


        /* state MATCH85 is a source for special END */ 
        temp = score + (0) + (0) ;   
        if( temp > LocalDnaMatchBlock_DC_OPT_SHADOW_SPECIAL(mat,i,j,END) )   {  
          LocalDnaMatchBlock_DC_OPT_SHADOW_SPECIAL(mat,i,j,END) = temp;  
          /* Have to push only bottem half of system here */ 
          for(k=0;k<3;k++)   
            LocalDnaMatchBlock_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,k) = LocalDnaMatchBlock_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,MATCH85,k);  
          LocalDnaMatchBlock_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,6) = LocalDnaMatchBlock_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,MATCH85,6);    
          LocalDnaMatchBlock_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,3) = i;    
          LocalDnaMatchBlock_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,4) = j;    
          LocalDnaMatchBlock_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,5) = MATCH85;  
          }  




        /* Finished calculating state MATCH85 */ 


        /* For state MATCH95 */ 
        /* setting first movement to score */ 
        score = LocalDnaMatchBlock_DC_OPT_SHADOW_MATRIX(mat,i-1,j-1,MATCH95) + (mat->lchs->comp95->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->s) + (0);     
        /* assign local shadown pointer */ 
        localsp = &(LocalDnaMatchBlock_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 1,MATCH95,0));  
        /* From state MATCH95 to state MATCH95 */ 
        temp = LocalDnaMatchBlock_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,MATCH95) + mat->lchs->g +(0);     
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(LocalDnaMatchBlock_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,MATCH95,0));    
          }  
        /* From state MATCH95 to state MATCH95 */ 
        temp = LocalDnaMatchBlock_DC_OPT_SHADOW_MATRIX(mat,i-1,j-0,MATCH95) + mat->lchs->g +(0);     
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(LocalDnaMatchBlock_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 0,MATCH95,0));    
          }  
        /* From state START to state MATCH95 */ 
        temp = LocalDnaMatchBlock_DC_OPT_SHADOW_SPECIAL(mat,i-1,j-1,START) + (mat->lchs->comp95->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->v) + (0);   
        if( temp  > score )  {  
          score = temp;  
          /* This state [START] is a special for MATCH95... push top shadow pointers here */ 
          localshadow[0]= i; 
          localshadow[1]= j; 
          localshadow[2]= MATCH95;   
          localshadow[3]= (-1);  
          localshadow[4]= (-1);  
          localshadow[5]= (-1);  
          localshadow[6]= score; 
          localsp = localshadow; 
          }  
        /* From state UNMATCHED_TARGET to state MATCH95 */ 
        temp = LocalDnaMatchBlock_DC_OPT_SHADOW_MATRIX(mat,i-1,j-1,UNMATCHED_TARGET) + (mat->lchs->comp95->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->v) +(0);  
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(LocalDnaMatchBlock_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 1,UNMATCHED_TARGET,0));   
          }  


        /* Ok - finished max calculation for MATCH95 */ 
        /* Add any movement independant score and put away */ 
        /* Actually, already done inside scores */ 
         LocalDnaMatchBlock_DC_OPT_SHADOW_MATRIX(mat,i,j,MATCH95) = score;   
        for(k=0;k<7;k++) 
          LocalDnaMatchBlock_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,MATCH95,k) = localsp[k];    
        /* Now figure out if any specials need this score */ 


        /* state MATCH95 is a source for special END */ 
        temp = score + (0) + (0) ;   
        if( temp > LocalDnaMatchBlock_DC_OPT_SHADOW_SPECIAL(mat,i,j,END) )   {  
          LocalDnaMatchBlock_DC_OPT_SHADOW_SPECIAL(mat,i,j,END) = temp;  
          /* Have to push only bottem half of system here */ 
          for(k=0;k<3;k++)   
            LocalDnaMatchBlock_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,k) = LocalDnaMatchBlock_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,MATCH95,k);  
          LocalDnaMatchBlock_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,6) = LocalDnaMatchBlock_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,MATCH95,6);    
          LocalDnaMatchBlock_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,3) = i;    
          LocalDnaMatchBlock_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,4) = j;    
          LocalDnaMatchBlock_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,5) = MATCH95;  
          }  




        /* Finished calculating state MATCH95 */ 


        /* For state UNMATCHED_QUERY */ 
        /* setting first movement to score */ 
        score = LocalDnaMatchBlock_DC_OPT_SHADOW_MATRIX(mat,i-1,j-0,MATCH65) + mat->lchs->b + (0);   
        /* assign local shadown pointer */ 
        localsp = &(LocalDnaMatchBlock_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 0,MATCH65,0));  
        /* From state MATCH75 to state UNMATCHED_QUERY */ 
        temp = LocalDnaMatchBlock_DC_OPT_SHADOW_MATRIX(mat,i-1,j-0,MATCH75) + mat->lchs->b +(0);     
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(LocalDnaMatchBlock_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 0,MATCH75,0));    
          }  
        /* From state MATCH85 to state UNMATCHED_QUERY */ 
        temp = LocalDnaMatchBlock_DC_OPT_SHADOW_MATRIX(mat,i-1,j-0,MATCH85) + mat->lchs->b +(0);     
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(LocalDnaMatchBlock_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 0,MATCH85,0));    
          }  
        /* From state MATCH95 to state UNMATCHED_QUERY */ 
        temp = LocalDnaMatchBlock_DC_OPT_SHADOW_MATRIX(mat,i-1,j-0,MATCH95) + mat->lchs->b +(0);     
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(LocalDnaMatchBlock_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 0,MATCH95,0));    
          }  
        /* From state UNMATCHED_QUERY to state UNMATCHED_QUERY */ 
        temp = LocalDnaMatchBlock_DC_OPT_SHADOW_MATRIX(mat,i-1,j-0,UNMATCHED_QUERY) + mat->lchs->u +(0);     
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(LocalDnaMatchBlock_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 0,UNMATCHED_QUERY,0));    
          }  


        /* Ok - finished max calculation for UNMATCHED_QUERY */ 
        /* Add any movement independant score and put away */ 
        /* Actually, already done inside scores */ 
         LocalDnaMatchBlock_DC_OPT_SHADOW_MATRIX(mat,i,j,UNMATCHED_QUERY) = score;   
        for(k=0;k<7;k++) 
          LocalDnaMatchBlock_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED_QUERY,k) = localsp[k];    
        /* Now figure out if any specials need this score */ 


        /* Finished calculating state UNMATCHED_QUERY */ 


        /* For state UNMATCHED_TARGET */ 
        /* setting first movement to score */ 
        score = LocalDnaMatchBlock_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,UNMATCHED_QUERY) + mat->lchs->v + (0);   
        /* assign local shadown pointer */ 
        localsp = &(LocalDnaMatchBlock_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,UNMATCHED_QUERY,0));  
        /* From state UNMATCHED_TARGET to state UNMATCHED_TARGET */ 
        temp = LocalDnaMatchBlock_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,UNMATCHED_TARGET) + mat->lchs->u +(0);    
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(LocalDnaMatchBlock_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,UNMATCHED_TARGET,0));   
          }  


        /* Ok - finished max calculation for UNMATCHED_TARGET */ 
        /* Add any movement independant score and put away */ 
        /* Actually, already done inside scores */ 
         LocalDnaMatchBlock_DC_OPT_SHADOW_MATRIX(mat,i,j,UNMATCHED_TARGET) = score;  
        for(k=0;k<7;k++) 
          LocalDnaMatchBlock_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED_TARGET,k) = localsp[k];   
        /* Now figure out if any specials need this score */ 


        /* Finished calculating state UNMATCHED_TARGET */ 


        } /* end of for each i position in strip */ 
      } /* end of for each j strip */ 
    free(score_pointers);    
    free(shadow_pointers);   
    return TRUE;     
}    


/* Function:  init_start_end_linear_LocalDnaMatchBlock(mat)
 *
 * Descrip: No Description
 *
 * Arg:        mat [UNKN ] Undocumented argument [LocalDnaMatchBlock *]
 *
 */
void init_start_end_linear_LocalDnaMatchBlock(LocalDnaMatchBlock * mat) 
{
    register int i;  
    register int j;  
    for(j=0;j<3;j++) {  
      for(i=(-1);i<mat->query->seq->len;i++) {  
        LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i,j,MATCH65) = NEGI; 
        LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH65,0) = (-1);    
        LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i,j,MATCH75) = NEGI; 
        LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH75,0) = (-1);    
        LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i,j,MATCH85) = NEGI; 
        LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH85,0) = (-1);    
        LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i,j,MATCH95) = NEGI; 
        LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH95,0) = (-1);    
        LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i,j,UNMATCHED_QUERY) = NEGI; 
        LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED_QUERY,0) = (-1);    
        LocalDnaMatchBlock_DC_SHADOW_MATRIX(mat,i,j,UNMATCHED_TARGET) = NEGI;    
        LocalDnaMatchBlock_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED_TARGET,0) = (-1);   
        }  
      }  


    for(j=(-1);j<mat->target->seq->len;j++)  {  
      LocalDnaMatchBlock_DC_SHADOW_SPECIAL(mat,0,j,START) = 0;   
      LocalDnaMatchBlock_DC_SHADOW_SPECIAL_SP(mat,0,j,START,0) = j;  
      LocalDnaMatchBlock_DC_SHADOW_SPECIAL(mat,0,j,END) = NEGI;  
      LocalDnaMatchBlock_DC_SHADOW_SPECIAL_SP(mat,0,j,END,0) = (-1); 
      }  


    return;  
}    


/* Function:  convert_PackAln_to_AlnBlock_LocalDnaMatchBlock(pal)
 *
 * Descrip:    Converts a path alignment to a label alignment
 *             The label alignment is probably much more useful than the path
 *
 *
 * Arg:        pal [UNKN ] Undocumented argument [PackAln *]
 *
 * Return [UNKN ]  Undocumented return value [AlnBlock *]
 *
 */
AlnBlock * convert_PackAln_to_AlnBlock_LocalDnaMatchBlock(PackAln * pal) 
{
    AlnConvertSet * acs; 
    AlnBlock * alb;  


    acs = AlnConvertSet_LocalDnaMatchBlock();    
    alb = AlnBlock_from_PackAln(acs,pal);    
    free_AlnConvertSet(acs); 
    return alb;  
}    


 static char * query_label[] = { "MM65","MI65","MM75","MI75","MM85","MI85","MM95","MI95","UM","UI","END" };  
/* Function:  AlnConvertSet_LocalDnaMatchBlock(void)
 *
 * Descrip: No Description
 *
 *
 * Return [UNKN ]  Undocumented return value [AlnConvertSet *]
 *
 */
 static char * target_label[] = { "MM65","MI65","MM75","MI75","MM85","MI85","MM95","MI95","UI","UM","END" }; 
AlnConvertSet * AlnConvertSet_LocalDnaMatchBlock(void) 
{
    AlnConvertUnit * acu;    
    AlnConvertSet  * out;    


    out = AlnConvertSet_alloc_std(); 


    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = MATCH65;   
    acu->state2 = MATCH65;   
    acu->offi = 1;   
    acu->offj = 1;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[0];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = MATCH65;   
    acu->state2 = MATCH65;   
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[1];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = MATCH65;   
    acu->state2 = MATCH65;   
    acu->offi = 1;   
    acu->offj = 0;   
    acu->label1 = query_label[1];    
    acu->label2 = target_label[0];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = START + 6; 
    acu->is_from_special = TRUE; 
    acu->state2 = MATCH65;   
    acu->offi = (-1);    
    acu->offj = 1;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[0];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = UNMATCHED_TARGET;  
    acu->state2 = MATCH65;   
    acu->offi = 1;   
    acu->offj = 1;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[0];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = MATCH75;   
    acu->state2 = MATCH75;   
    acu->offi = 1;   
    acu->offj = 1;   
    acu->label1 = query_label[2];    
    acu->label2 = target_label[2];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = MATCH75;   
    acu->state2 = MATCH75;   
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[2];    
    acu->label2 = target_label[3];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = MATCH75;   
    acu->state2 = MATCH75;   
    acu->offi = 1;   
    acu->offj = 0;   
    acu->label1 = query_label[3];    
    acu->label2 = target_label[2];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = START + 6; 
    acu->is_from_special = TRUE; 
    acu->state2 = MATCH75;   
    acu->offi = (-1);    
    acu->offj = 1;   
    acu->label1 = query_label[2];    
    acu->label2 = target_label[2];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = UNMATCHED_TARGET;  
    acu->state2 = MATCH75;   
    acu->offi = 1;   
    acu->offj = 1;   
    acu->label1 = query_label[2];    
    acu->label2 = target_label[2];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = MATCH85;   
    acu->state2 = MATCH85;   
    acu->offi = 1;   
    acu->offj = 1;   
    acu->label1 = query_label[4];    
    acu->label2 = target_label[4];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = MATCH85;   
    acu->state2 = MATCH85;   
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[4];    
    acu->label2 = target_label[5];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = MATCH85;   
    acu->state2 = MATCH85;   
    acu->offi = 1;   
    acu->offj = 0;   
    acu->label1 = query_label[5];    
    acu->label2 = target_label[4];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = START + 6; 
    acu->is_from_special = TRUE; 
    acu->state2 = MATCH85;   
    acu->offi = (-1);    
    acu->offj = 1;   
    acu->label1 = query_label[4];    
    acu->label2 = target_label[4];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = UNMATCHED_TARGET;  
    acu->state2 = MATCH85;   
    acu->offi = 1;   
    acu->offj = 1;   
    acu->label1 = query_label[4];    
    acu->label2 = target_label[4];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = MATCH95;   
    acu->state2 = MATCH95;   
    acu->offi = 1;   
    acu->offj = 1;   
    acu->label1 = query_label[6];    
    acu->label2 = target_label[6];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = MATCH95;   
    acu->state2 = MATCH95;   
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[6];    
    acu->label2 = target_label[7];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = MATCH95;   
    acu->state2 = MATCH95;   
    acu->offi = 1;   
    acu->offj = 0;   
    acu->label1 = query_label[7];    
    acu->label2 = target_label[6];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = START + 6; 
    acu->is_from_special = TRUE; 
    acu->state2 = MATCH95;   
    acu->offi = (-1);    
    acu->offj = 1;   
    acu->label1 = query_label[6];    
    acu->label2 = target_label[6];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = UNMATCHED_TARGET;  
    acu->state2 = MATCH95;   
    acu->offi = 1;   
    acu->offj = 1;   
    acu->label1 = query_label[6];    
    acu->label2 = target_label[6];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = MATCH65;   
    acu->state2 = UNMATCHED_QUERY;   
    acu->offi = 1;   
    acu->offj = 0;   
    acu->label1 = query_label[8];    
    acu->label2 = target_label[8];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = MATCH75;   
    acu->state2 = UNMATCHED_QUERY;   
    acu->offi = 1;   
    acu->offj = 0;   
    acu->label1 = query_label[8];    
    acu->label2 = target_label[8];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = MATCH85;   
    acu->state2 = UNMATCHED_QUERY;   
    acu->offi = 1;   
    acu->offj = 0;   
    acu->label1 = query_label[8];    
    acu->label2 = target_label[8];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = MATCH95;   
    acu->state2 = UNMATCHED_QUERY;   
    acu->offi = 1;   
    acu->offj = 0;   
    acu->label1 = query_label[8];    
    acu->label2 = target_label[8];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = UNMATCHED_QUERY;   
    acu->state2 = UNMATCHED_QUERY;   
    acu->offi = 1;   
    acu->offj = 0;   
    acu->label1 = query_label[8];    
    acu->label2 = target_label[8];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = UNMATCHED_QUERY;   
    acu->state2 = UNMATCHED_TARGET;  
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[9];    
    acu->label2 = target_label[9];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = UNMATCHED_TARGET;  
    acu->state2 = UNMATCHED_TARGET;  
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[9];    
    acu->label2 = target_label[9];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = MATCH65;   
    acu->state2 = END + 6;   
    acu->offi = (-1);    
    acu->offj = 0;   
    acu->label1 = query_label[10];   
    acu->label2 = target_label[10];  
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = MATCH75;   
    acu->state2 = END + 6;   
    acu->offi = (-1);    
    acu->offj = 0;   
    acu->label1 = query_label[10];   
    acu->label2 = target_label[10];  
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = MATCH85;   
    acu->state2 = END + 6;   
    acu->offi = (-1);    
    acu->offj = 0;   
    acu->label1 = query_label[10];   
    acu->label2 = target_label[10];  
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = MATCH95;   
    acu->state2 = END + 6;   
    acu->offi = (-1);    
    acu->offj = 0;   
    acu->label1 = query_label[10];   
    acu->label2 = target_label[10];  
    return out;  
}    


/* Function:  PackAln_read_Expl_LocalDnaMatchBlock(mat)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [LocalDnaMatchBlock *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * PackAln_read_Expl_LocalDnaMatchBlock(LocalDnaMatchBlock * mat) 
{
    LocalDnaMatchBlock_access_func_holder holder;    


    holder.access_main    = LocalDnaMatchBlock_explicit_access_main; 
    holder.access_special = LocalDnaMatchBlock_explicit_access_special;  
    return PackAln_read_generic_LocalDnaMatchBlock(mat,holder);  
}    


/* Function:  LocalDnaMatchBlock_explicit_access_main(mat,i,j,state)
 *
 * Descrip: No Description
 *
 * Arg:          mat [UNKN ] Undocumented argument [LocalDnaMatchBlock *]
 * Arg:            i [UNKN ] Undocumented argument [int]
 * Arg:            j [UNKN ] Undocumented argument [int]
 * Arg:        state [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int LocalDnaMatchBlock_explicit_access_main(LocalDnaMatchBlock * mat,int i,int j,int state) 
{
    return LocalDnaMatchBlock_EXPL_MATRIX(mat,i,j,state);    
}    


/* Function:  LocalDnaMatchBlock_explicit_access_special(mat,i,j,state)
 *
 * Descrip: No Description
 *
 * Arg:          mat [UNKN ] Undocumented argument [LocalDnaMatchBlock *]
 * Arg:            i [UNKN ] Undocumented argument [int]
 * Arg:            j [UNKN ] Undocumented argument [int]
 * Arg:        state [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int LocalDnaMatchBlock_explicit_access_special(LocalDnaMatchBlock * mat,int i,int j,int state) 
{
    return LocalDnaMatchBlock_EXPL_SPECIAL(mat,i,j,state);   
}    


/* Function:  PackAln_read_generic_LocalDnaMatchBlock(mat,h)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [LocalDnaMatchBlock *]
 * Arg:          h [UNKN ] Undocumented argument [LocalDnaMatchBlock_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * PackAln_read_generic_LocalDnaMatchBlock(LocalDnaMatchBlock * mat,LocalDnaMatchBlock_access_func_holder h) 
{
    register PackAln * out;  
    int i;   
    int j;   
    int state;   
    int cellscore = (-1);    
    boolean isspecial;   
    PackAlnUnit * pau = NULL;    
    PackAlnUnit * prev = NULL;   


    assert(mat);     
    assert(h.access_main);   
    assert(h.access_special);    


    out = PackAln_alloc_std();   
    if( out == NULL )    
      return NULL;   


    out->score =  find_end_LocalDnaMatchBlock(mat,&i,&j,&state,&isspecial,h);    


    /* Add final end transition (at the moment we have not got the score! */ 
    if( (pau= PackAlnUnit_alloc()) == NULL  || add_PackAln(out,pau) == FALSE )   {  
      warn("Failed the first PackAlnUnit alloc, %d length of Alignment in LocalDnaMatchBlock_basic_read, returning a mess.(Sorry!)",out->len);   
      return out;    
      }  


    /* Put in positions for end trans. Remember that coordinates in C style */ 
    pau->i = i;  
    pau->j = j;  
    if( isspecial != TRUE)   
      pau->state = state;    
    else pau->state = state + 6;     
    prev=pau;    
    while( state != START || isspecial != TRUE)  { /*while state != START*/ 


      if( isspecial == TRUE )    
        max_calc_special_LocalDnaMatchBlock(mat,i,j,state,isspecial,&i,&j,&state,&isspecial,&cellscore,h);   
      else   
        max_calc_LocalDnaMatchBlock(mat,i,j,state,isspecial,&i,&j,&state,&isspecial,&cellscore,h);   
      if(i == LocalDnaMatchBlock_READ_OFF_ERROR || j == LocalDnaMatchBlock_READ_OFF_ERROR || state == LocalDnaMatchBlock_READ_OFF_ERROR )    {  
        warn("Problem - hit bad read off system, exiting now");  
        break;   
        }  
      if( (pau= PackAlnUnit_alloc()) == NULL  || add_PackAln(out,pau) == FALSE ) {  
        warn("Failed a PackAlnUnit alloc, %d length of Alignment in LocalDnaMatchBlock_basic_read, returning partial alignment",out->len);   
        break;   
        }  


      /* Put in positions for block. Remember that coordinates in C style */ 
      pau->i = i;    
      pau->j = j;    
      if( isspecial != TRUE)     
        pau->state = state;  
      else pau->state = state + 6;   
      prev->score = cellscore;   
      prev = pau;    
      } /* end of while state != START */ 


    invert_PackAln(out); 
    return out;  
}    


/* Function:  find_end_LocalDnaMatchBlock(mat,ri,rj,state,isspecial,h)
 *
 * Descrip: No Description
 *
 * Arg:              mat [UNKN ] Undocumented argument [LocalDnaMatchBlock *]
 * Arg:               ri [UNKN ] Undocumented argument [int *]
 * Arg:               rj [UNKN ] Undocumented argument [int *]
 * Arg:            state [UNKN ] Undocumented argument [int *]
 * Arg:        isspecial [UNKN ] Undocumented argument [boolean *]
 * Arg:                h [UNKN ] Undocumented argument [LocalDnaMatchBlock_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int find_end_LocalDnaMatchBlock(LocalDnaMatchBlock * mat,int * ri,int * rj,int * state,boolean * isspecial,LocalDnaMatchBlock_access_func_holder h) 
{
    int j;   
    int max; 
    int maxj;    
    int temp;    


    max = (*h.access_special)(mat,0,mat->target->seq->len-1,END);    
    maxj = mat->target->seq->len-1;  
    for(j= mat->target->seq->len-2 ;j >= 0 ;j--) {  
      if( (temp =(*h.access_special)(mat,0,j,END)) > max )   {  
        max = temp;  
        maxj = j;    
        }  
      }  


    if( ri != NULL)  
       *ri = 0;  
    if( rj != NULL)  
       *rj = maxj;   
    if( state != NULL)   
       *state = END; 
    if( isspecial != NULL)   
       *isspecial = TRUE;    


    return max;  
}    


/* Function:  LocalDnaMatchBlock_debug_show_matrix(mat,starti,stopi,startj,stopj,ofp)
 *
 * Descrip: No Description
 *
 * Arg:           mat [UNKN ] Undocumented argument [LocalDnaMatchBlock *]
 * Arg:        starti [UNKN ] Undocumented argument [int]
 * Arg:         stopi [UNKN ] Undocumented argument [int]
 * Arg:        startj [UNKN ] Undocumented argument [int]
 * Arg:         stopj [UNKN ] Undocumented argument [int]
 * Arg:           ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
void LocalDnaMatchBlock_debug_show_matrix(LocalDnaMatchBlock * mat,int starti,int stopi,int startj,int stopj,FILE * ofp) 
{
    register int i;  
    register int j;  


    for(i=starti;i<stopi && i < mat->query->seq->len;i++)    {  
      for(j=startj;j<stopj && j < mat->target->seq->len;j++) {  
        fprintf(ofp,"Cell [%d - %d]\n",i,j);     
        fprintf(ofp,"State MATCH65 %d\n",LocalDnaMatchBlock_EXPL_MATRIX(mat,i,j,MATCH65));   
        fprintf(ofp,"State MATCH75 %d\n",LocalDnaMatchBlock_EXPL_MATRIX(mat,i,j,MATCH75));   
        fprintf(ofp,"State MATCH85 %d\n",LocalDnaMatchBlock_EXPL_MATRIX(mat,i,j,MATCH85));   
        fprintf(ofp,"State MATCH95 %d\n",LocalDnaMatchBlock_EXPL_MATRIX(mat,i,j,MATCH95));   
        fprintf(ofp,"State UNMATCHED_QUERY %d\n",LocalDnaMatchBlock_EXPL_MATRIX(mat,i,j,UNMATCHED_QUERY));   
        fprintf(ofp,"State UNMATCHED_TARGET %d\n",LocalDnaMatchBlock_EXPL_MATRIX(mat,i,j,UNMATCHED_TARGET)); 
        fprintf(ofp,"\n\n"); 
        }  
      }  


}    


/* Function:  max_calc_LocalDnaMatchBlock(mat,i,j,state,isspecial,reti,retj,retstate,retspecial,cellscore,h)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [LocalDnaMatchBlock *]
 * Arg:                 i [UNKN ] Undocumented argument [int]
 * Arg:                 j [UNKN ] Undocumented argument [int]
 * Arg:             state [UNKN ] Undocumented argument [int]
 * Arg:         isspecial [UNKN ] Undocumented argument [boolean]
 * Arg:              reti [UNKN ] Undocumented argument [int *]
 * Arg:              retj [UNKN ] Undocumented argument [int *]
 * Arg:          retstate [UNKN ] Undocumented argument [int *]
 * Arg:        retspecial [UNKN ] Undocumented argument [boolean *]
 * Arg:         cellscore [UNKN ] Undocumented argument [int *]
 * Arg:                 h [UNKN ] Undocumented argument [LocalDnaMatchBlock_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int max_calc_LocalDnaMatchBlock(LocalDnaMatchBlock * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,LocalDnaMatchBlock_access_func_holder h) 
{
    register int temp;   
    register int cscore; 


    *reti = (*retj) = (*retstate) = LocalDnaMatchBlock_READ_OFF_ERROR;   


    if( i < 0 || j < 0 || i > mat->query->seq->len || j > mat->target->seq->len) {  
      warn("In LocalDnaMatchBlock matrix special read off - out of bounds on matrix [i,j is %d,%d state %d in standard matrix]",i,j,state);  
      return -1;     
      }  


    /* Then you have to select the correct switch statement to figure out the readoff      */ 
    /* Somewhat odd - reverse the order of calculation and return as soon as it is correct */ 
    cscore = (*h.access_main)(mat,i,j,state);    
    switch(state)    { /*Switch state */ 
      case MATCH65 :     
        temp = cscore - ((mat->lchs->comp65->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->v)) -  (0); 
        if( temp == (*h.access_main)(mat,i - 1,j - 1,UNMATCHED_TARGET) ) {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = UNMATCHED_TARGET;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-1,UNMATCHED_TARGET);    
            }  
          return (*h.access_main)(mat,i - 1,j - 1,UNMATCHED_TARGET);     
          }  
        temp = cscore - ((mat->lchs->comp65->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->v)) -  (0); 
        if( temp == (*h.access_special)(mat,i - 1,j - 1,START) ) {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = START; 
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_special)(mat,i-1,j-1,START);    
            }  
          return (*h.access_main)(mat,i - 1,j - 1,START);    
          }  
        temp = cscore - (mat->lchs->g) -  (0);   
        if( temp == (*h.access_main)(mat,i - 1,j - 0,MATCH65) )  {  
          *reti = i - 1; 
          *retj = j - 0; 
          *retstate = MATCH65;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-0,MATCH65); 
            }  
          return (*h.access_main)(mat,i - 1,j - 0,MATCH65);  
          }  
        temp = cscore - (mat->lchs->g) -  (0);   
        if( temp == (*h.access_main)(mat,i - 0,j - 1,MATCH65) )  {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = MATCH65;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,MATCH65); 
            }  
          return (*h.access_main)(mat,i - 0,j - 1,MATCH65);  
          }  
        temp = cscore - ((mat->lchs->comp65->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->s)) -  (0); 
        if( temp == (*h.access_main)(mat,i - 1,j - 1,MATCH65) )  {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = MATCH65;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-1,MATCH65); 
            }  
          return (*h.access_main)(mat,i - 1,j - 1,MATCH65);  
          }  
        warn("Major problem (!) - in LocalDnaMatchBlock read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      case MATCH75 :     
        temp = cscore - ((mat->lchs->comp75->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->v)) -  (0); 
        if( temp == (*h.access_main)(mat,i - 1,j - 1,UNMATCHED_TARGET) ) {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = UNMATCHED_TARGET;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-1,UNMATCHED_TARGET);    
            }  
          return (*h.access_main)(mat,i - 1,j - 1,UNMATCHED_TARGET);     
          }  
        temp = cscore - ((mat->lchs->comp75->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->v)) -  (0); 
        if( temp == (*h.access_special)(mat,i - 1,j - 1,START) ) {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = START; 
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_special)(mat,i-1,j-1,START);    
            }  
          return (*h.access_main)(mat,i - 1,j - 1,START);    
          }  
        temp = cscore - (mat->lchs->g) -  (0);   
        if( temp == (*h.access_main)(mat,i - 1,j - 0,MATCH75) )  {  
          *reti = i - 1; 
          *retj = j - 0; 
          *retstate = MATCH75;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-0,MATCH75); 
            }  
          return (*h.access_main)(mat,i - 1,j - 0,MATCH75);  
          }  
        temp = cscore - (mat->lchs->g) -  (0);   
        if( temp == (*h.access_main)(mat,i - 0,j - 1,MATCH75) )  {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = MATCH75;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,MATCH75); 
            }  
          return (*h.access_main)(mat,i - 0,j - 1,MATCH75);  
          }  
        temp = cscore - ((mat->lchs->comp75->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->s)) -  (0); 
        if( temp == (*h.access_main)(mat,i - 1,j - 1,MATCH75) )  {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = MATCH75;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-1,MATCH75); 
            }  
          return (*h.access_main)(mat,i - 1,j - 1,MATCH75);  
          }  
        warn("Major problem (!) - in LocalDnaMatchBlock read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      case MATCH85 :     
        temp = cscore - ((mat->lchs->comp85->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->v)) -  (0); 
        if( temp == (*h.access_main)(mat,i - 1,j - 1,UNMATCHED_TARGET) ) {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = UNMATCHED_TARGET;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-1,UNMATCHED_TARGET);    
            }  
          return (*h.access_main)(mat,i - 1,j - 1,UNMATCHED_TARGET);     
          }  
        temp = cscore - ((mat->lchs->comp85->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->v)) -  (0); 
        if( temp == (*h.access_special)(mat,i - 1,j - 1,START) ) {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = START; 
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_special)(mat,i-1,j-1,START);    
            }  
          return (*h.access_main)(mat,i - 1,j - 1,START);    
          }  
        temp = cscore - (mat->lchs->g) -  (0);   
        if( temp == (*h.access_main)(mat,i - 1,j - 0,MATCH85) )  {  
          *reti = i - 1; 
          *retj = j - 0; 
          *retstate = MATCH85;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-0,MATCH85); 
            }  
          return (*h.access_main)(mat,i - 1,j - 0,MATCH85);  
          }  
        temp = cscore - (mat->lchs->g) -  (0);   
        if( temp == (*h.access_main)(mat,i - 0,j - 1,MATCH85) )  {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = MATCH85;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,MATCH85); 
            }  
          return (*h.access_main)(mat,i - 0,j - 1,MATCH85);  
          }  
        temp = cscore - ((mat->lchs->comp85->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->s)) -  (0); 
        if( temp == (*h.access_main)(mat,i - 1,j - 1,MATCH85) )  {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = MATCH85;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-1,MATCH85); 
            }  
          return (*h.access_main)(mat,i - 1,j - 1,MATCH85);  
          }  
        warn("Major problem (!) - in LocalDnaMatchBlock read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      case MATCH95 :     
        temp = cscore - ((mat->lchs->comp95->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->v)) -  (0); 
        if( temp == (*h.access_main)(mat,i - 1,j - 1,UNMATCHED_TARGET) ) {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = UNMATCHED_TARGET;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-1,UNMATCHED_TARGET);    
            }  
          return (*h.access_main)(mat,i - 1,j - 1,UNMATCHED_TARGET);     
          }  
        temp = cscore - ((mat->lchs->comp95->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->v)) -  (0); 
        if( temp == (*h.access_special)(mat,i - 1,j - 1,START) ) {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = START; 
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_special)(mat,i-1,j-1,START);    
            }  
          return (*h.access_main)(mat,i - 1,j - 1,START);    
          }  
        temp = cscore - (mat->lchs->g) -  (0);   
        if( temp == (*h.access_main)(mat,i - 1,j - 0,MATCH95) )  {  
          *reti = i - 1; 
          *retj = j - 0; 
          *retstate = MATCH95;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-0,MATCH95); 
            }  
          return (*h.access_main)(mat,i - 1,j - 0,MATCH95);  
          }  
        temp = cscore - (mat->lchs->g) -  (0);   
        if( temp == (*h.access_main)(mat,i - 0,j - 1,MATCH95) )  {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = MATCH95;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,MATCH95); 
            }  
          return (*h.access_main)(mat,i - 0,j - 1,MATCH95);  
          }  
        temp = cscore - ((mat->lchs->comp95->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->s)) -  (0); 
        if( temp == (*h.access_main)(mat,i - 1,j - 1,MATCH95) )  {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = MATCH95;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-1,MATCH95); 
            }  
          return (*h.access_main)(mat,i - 1,j - 1,MATCH95);  
          }  
        warn("Major problem (!) - in LocalDnaMatchBlock read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      case UNMATCHED_QUERY :     
        temp = cscore - (mat->lchs->u) -  (0);   
        if( temp == (*h.access_main)(mat,i - 1,j - 0,UNMATCHED_QUERY) )  {  
          *reti = i - 1; 
          *retj = j - 0; 
          *retstate = UNMATCHED_QUERY;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-0,UNMATCHED_QUERY); 
            }  
          return (*h.access_main)(mat,i - 1,j - 0,UNMATCHED_QUERY);  
          }  
        temp = cscore - (mat->lchs->b) -  (0);   
        if( temp == (*h.access_main)(mat,i - 1,j - 0,MATCH95) )  {  
          *reti = i - 1; 
          *retj = j - 0; 
          *retstate = MATCH95;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-0,MATCH95); 
            }  
          return (*h.access_main)(mat,i - 1,j - 0,MATCH95);  
          }  
        temp = cscore - (mat->lchs->b) -  (0);   
        if( temp == (*h.access_main)(mat,i - 1,j - 0,MATCH85) )  {  
          *reti = i - 1; 
          *retj = j - 0; 
          *retstate = MATCH85;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-0,MATCH85); 
            }  
          return (*h.access_main)(mat,i - 1,j - 0,MATCH85);  
          }  
        temp = cscore - (mat->lchs->b) -  (0);   
        if( temp == (*h.access_main)(mat,i - 1,j - 0,MATCH75) )  {  
          *reti = i - 1; 
          *retj = j - 0; 
          *retstate = MATCH75;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-0,MATCH75); 
            }  
          return (*h.access_main)(mat,i - 1,j - 0,MATCH75);  
          }  
        temp = cscore - (mat->lchs->b) -  (0);   
        if( temp == (*h.access_main)(mat,i - 1,j - 0,MATCH65) )  {  
          *reti = i - 1; 
          *retj = j - 0; 
          *retstate = MATCH65;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-0,MATCH65); 
            }  
          return (*h.access_main)(mat,i - 1,j - 0,MATCH65);  
          }  
        warn("Major problem (!) - in LocalDnaMatchBlock read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      case UNMATCHED_TARGET :    
        temp = cscore - (mat->lchs->u) -  (0);   
        if( temp == (*h.access_main)(mat,i - 0,j - 1,UNMATCHED_TARGET) ) {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = UNMATCHED_TARGET;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,UNMATCHED_TARGET);    
            }  
          return (*h.access_main)(mat,i - 0,j - 1,UNMATCHED_TARGET);     
          }  
        temp = cscore - (mat->lchs->v) -  (0);   
        if( temp == (*h.access_main)(mat,i - 0,j - 1,UNMATCHED_QUERY) )  {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = UNMATCHED_QUERY;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,UNMATCHED_QUERY); 
            }  
          return (*h.access_main)(mat,i - 0,j - 1,UNMATCHED_QUERY);  
          }  
        warn("Major problem (!) - in LocalDnaMatchBlock read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      default:   
        warn("Major problem (!) - in LocalDnaMatchBlock read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      } /* end of Switch state  */ 
}    


/* Function:  max_calc_special_LocalDnaMatchBlock(mat,i,j,state,isspecial,reti,retj,retstate,retspecial,cellscore,h)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [LocalDnaMatchBlock *]
 * Arg:                 i [UNKN ] Undocumented argument [int]
 * Arg:                 j [UNKN ] Undocumented argument [int]
 * Arg:             state [UNKN ] Undocumented argument [int]
 * Arg:         isspecial [UNKN ] Undocumented argument [boolean]
 * Arg:              reti [UNKN ] Undocumented argument [int *]
 * Arg:              retj [UNKN ] Undocumented argument [int *]
 * Arg:          retstate [UNKN ] Undocumented argument [int *]
 * Arg:        retspecial [UNKN ] Undocumented argument [boolean *]
 * Arg:         cellscore [UNKN ] Undocumented argument [int *]
 * Arg:                 h [UNKN ] Undocumented argument [LocalDnaMatchBlock_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int max_calc_special_LocalDnaMatchBlock(LocalDnaMatchBlock * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,LocalDnaMatchBlock_access_func_holder h) 
{
    register int temp;   
    register int cscore; 


    *reti = (*retj) = (*retstate) = LocalDnaMatchBlock_READ_OFF_ERROR;   


    if( j < 0 || j > mat->target->seq->len)  {  
      warn("In LocalDnaMatchBlock matrix special read off - out of bounds on matrix [j is %d in special]",j);    
      return -1;     
      }  


    cscore = (*h.access_special)(mat,i,j,state); 
    switch(state)    { /*switch on special states*/ 
      case START :   
      case END :     
        /* source MATCH95 is from main matrix */ 
        for(i= mat->query->seq->len-1;i >= 0 ;i--)   { /*for i >= 0*/ 
          temp = cscore - (0) - (0);     
          if( temp == (*h.access_main)(mat,i - 0,j - 0,MATCH95) )    {  
            *reti = i - 0;   
            *retj = j - 0;   
            *retstate = MATCH95; 
            *retspecial = FALSE; 
            if( cellscore != NULL)   {  
              *cellscore = cscore - (*h.access_main)(mat,i-0,j-0,MATCH95);   
              }  
            return (*h.access_main)(mat,i - 0,j - 0,MATCH95) ;   
            }  
          } /* end of for i >= 0 */ 
        /* source MATCH85 is from main matrix */ 
        for(i= mat->query->seq->len-1;i >= 0 ;i--)   { /*for i >= 0*/ 
          temp = cscore - (0) - (0);     
          if( temp == (*h.access_main)(mat,i - 0,j - 0,MATCH85) )    {  
            *reti = i - 0;   
            *retj = j - 0;   
            *retstate = MATCH85; 
            *retspecial = FALSE; 
            if( cellscore != NULL)   {  
              *cellscore = cscore - (*h.access_main)(mat,i-0,j-0,MATCH85);   
              }  
            return (*h.access_main)(mat,i - 0,j - 0,MATCH85) ;   
            }  
          } /* end of for i >= 0 */ 
        /* source MATCH75 is from main matrix */ 
        for(i= mat->query->seq->len-1;i >= 0 ;i--)   { /*for i >= 0*/ 
          temp = cscore - (0) - (0);     
          if( temp == (*h.access_main)(mat,i - 0,j - 0,MATCH75) )    {  
            *reti = i - 0;   
            *retj = j - 0;   
            *retstate = MATCH75; 
            *retspecial = FALSE; 
            if( cellscore != NULL)   {  
              *cellscore = cscore - (*h.access_main)(mat,i-0,j-0,MATCH75);   
              }  
            return (*h.access_main)(mat,i - 0,j - 0,MATCH75) ;   
            }  
          } /* end of for i >= 0 */ 
        /* source MATCH65 is from main matrix */ 
        for(i= mat->query->seq->len-1;i >= 0 ;i--)   { /*for i >= 0*/ 
          temp = cscore - (0) - (0);     
          if( temp == (*h.access_main)(mat,i - 0,j - 0,MATCH65) )    {  
            *reti = i - 0;   
            *retj = j - 0;   
            *retstate = MATCH65; 
            *retspecial = FALSE; 
            if( cellscore != NULL)   {  
              *cellscore = cscore - (*h.access_main)(mat,i-0,j-0,MATCH65);   
              }  
            return (*h.access_main)(mat,i - 0,j - 0,MATCH65) ;   
            }  
          } /* end of for i >= 0 */ 
      default:   
        warn("Major problem (!) - in LocalDnaMatchBlock read off, position %d,%d state %d no source found  dropped into default on source switch!",i,j,state);   
        return (-1); 
      } /* end of switch on special states */ 
}    


/* Function:  calculate_LocalDnaMatchBlock(mat)
 *
 * Descrip:    This function calculates the LocalDnaMatchBlock matrix when in explicit mode
 *             To allocate the matrix use /allocate_Expl_LocalDnaMatchBlock
 *
 *
 * Arg:        mat [UNKN ] LocalDnaMatchBlock which contains explicit basematrix memory [LocalDnaMatchBlock *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean calculate_LocalDnaMatchBlock(LocalDnaMatchBlock * mat) 
{
    int i;   
    int j;   
    int leni;    
    int lenj;    
    int tot; 
    int num; 


    if( mat->basematrix->type != BASEMATRIX_TYPE_EXPLICIT )  {  
      warn("in calculate_LocalDnaMatchBlock, passed a non Explicit matrix type, cannot calculate!"); 
      return FALSE;  
      }  


    leni = mat->leni;    
    lenj = mat->lenj;    
    tot = leni * lenj;   
    num = 0; 


    start_reporting("LocalDnaMatchBlock Matrix calculation: ");  
    for(j=0;j<lenj;j++)  {  
      auto int score;    
      auto int temp;     
      for(i=0;i<leni;i++)    {  
        if( num%1000 == 0 )  
          log_full_error(REPORT,0,"[%7d] Cells %2d%%%%",num,num*100/tot);    
        num++;   


        /* For state MATCH65 */ 
        /* setting first movement to score */ 
        score = LocalDnaMatchBlock_EXPL_MATRIX(mat,i-1,j-1,MATCH65) + (mat->lchs->comp65->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->s);    
        /* From state MATCH65 to state MATCH65 */ 
        temp = LocalDnaMatchBlock_EXPL_MATRIX(mat,i-0,j-1,MATCH65) + mat->lchs->g;   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH65 to state MATCH65 */ 
        temp = LocalDnaMatchBlock_EXPL_MATRIX(mat,i-1,j-0,MATCH65) + mat->lchs->g;   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state START to state MATCH65 */ 
        temp = LocalDnaMatchBlock_EXPL_SPECIAL(mat,i-1,j-1,START) + (mat->lchs->comp65->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->v);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state UNMATCHED_TARGET to state MATCH65 */ 
        temp = LocalDnaMatchBlock_EXPL_MATRIX(mat,i-1,j-1,UNMATCHED_TARGET) + (mat->lchs->comp65->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->v);    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for MATCH65 */ 
        /* Add any movement independant score and put away */ 
         LocalDnaMatchBlock_EXPL_MATRIX(mat,i,j,MATCH65) = score;    


        /* state MATCH65 is a source for special END */ 
        temp = score + (0) + (0) ;   
        if( temp > LocalDnaMatchBlock_EXPL_SPECIAL(mat,i,j,END) )    {  
          LocalDnaMatchBlock_EXPL_SPECIAL(mat,i,j,END) = temp;   
          }  




        /* Finished calculating state MATCH65 */ 


        /* For state MATCH75 */ 
        /* setting first movement to score */ 
        score = LocalDnaMatchBlock_EXPL_MATRIX(mat,i-1,j-1,MATCH75) + (mat->lchs->comp75->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->s);    
        /* From state MATCH75 to state MATCH75 */ 
        temp = LocalDnaMatchBlock_EXPL_MATRIX(mat,i-0,j-1,MATCH75) + mat->lchs->g;   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH75 to state MATCH75 */ 
        temp = LocalDnaMatchBlock_EXPL_MATRIX(mat,i-1,j-0,MATCH75) + mat->lchs->g;   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state START to state MATCH75 */ 
        temp = LocalDnaMatchBlock_EXPL_SPECIAL(mat,i-1,j-1,START) + (mat->lchs->comp75->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->v);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state UNMATCHED_TARGET to state MATCH75 */ 
        temp = LocalDnaMatchBlock_EXPL_MATRIX(mat,i-1,j-1,UNMATCHED_TARGET) + (mat->lchs->comp75->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->v);    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for MATCH75 */ 
        /* Add any movement independant score and put away */ 
         LocalDnaMatchBlock_EXPL_MATRIX(mat,i,j,MATCH75) = score;    


        /* state MATCH75 is a source for special END */ 
        temp = score + (0) + (0) ;   
        if( temp > LocalDnaMatchBlock_EXPL_SPECIAL(mat,i,j,END) )    {  
          LocalDnaMatchBlock_EXPL_SPECIAL(mat,i,j,END) = temp;   
          }  




        /* Finished calculating state MATCH75 */ 


        /* For state MATCH85 */ 
        /* setting first movement to score */ 
        score = LocalDnaMatchBlock_EXPL_MATRIX(mat,i-1,j-1,MATCH85) + (mat->lchs->comp85->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->s);    
        /* From state MATCH85 to state MATCH85 */ 
        temp = LocalDnaMatchBlock_EXPL_MATRIX(mat,i-0,j-1,MATCH85) + mat->lchs->g;   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH85 to state MATCH85 */ 
        temp = LocalDnaMatchBlock_EXPL_MATRIX(mat,i-1,j-0,MATCH85) + mat->lchs->g;   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state START to state MATCH85 */ 
        temp = LocalDnaMatchBlock_EXPL_SPECIAL(mat,i-1,j-1,START) + (mat->lchs->comp85->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->v);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state UNMATCHED_TARGET to state MATCH85 */ 
        temp = LocalDnaMatchBlock_EXPL_MATRIX(mat,i-1,j-1,UNMATCHED_TARGET) + (mat->lchs->comp85->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->v);    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for MATCH85 */ 
        /* Add any movement independant score and put away */ 
         LocalDnaMatchBlock_EXPL_MATRIX(mat,i,j,MATCH85) = score;    


        /* state MATCH85 is a source for special END */ 
        temp = score + (0) + (0) ;   
        if( temp > LocalDnaMatchBlock_EXPL_SPECIAL(mat,i,j,END) )    {  
          LocalDnaMatchBlock_EXPL_SPECIAL(mat,i,j,END) = temp;   
          }  




        /* Finished calculating state MATCH85 */ 


        /* For state MATCH95 */ 
        /* setting first movement to score */ 
        score = LocalDnaMatchBlock_EXPL_MATRIX(mat,i-1,j-1,MATCH95) + (mat->lchs->comp95->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->s);    
        /* From state MATCH95 to state MATCH95 */ 
        temp = LocalDnaMatchBlock_EXPL_MATRIX(mat,i-0,j-1,MATCH95) + mat->lchs->g;   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH95 to state MATCH95 */ 
        temp = LocalDnaMatchBlock_EXPL_MATRIX(mat,i-1,j-0,MATCH95) + mat->lchs->g;   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state START to state MATCH95 */ 
        temp = LocalDnaMatchBlock_EXPL_SPECIAL(mat,i-1,j-1,START) + (mat->lchs->comp95->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->v);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state UNMATCHED_TARGET to state MATCH95 */ 
        temp = LocalDnaMatchBlock_EXPL_MATRIX(mat,i-1,j-1,UNMATCHED_TARGET) + (mat->lchs->comp95->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->v);    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for MATCH95 */ 
        /* Add any movement independant score and put away */ 
         LocalDnaMatchBlock_EXPL_MATRIX(mat,i,j,MATCH95) = score;    


        /* state MATCH95 is a source for special END */ 
        temp = score + (0) + (0) ;   
        if( temp > LocalDnaMatchBlock_EXPL_SPECIAL(mat,i,j,END) )    {  
          LocalDnaMatchBlock_EXPL_SPECIAL(mat,i,j,END) = temp;   
          }  




        /* Finished calculating state MATCH95 */ 


        /* For state UNMATCHED_QUERY */ 
        /* setting first movement to score */ 
        score = LocalDnaMatchBlock_EXPL_MATRIX(mat,i-1,j-0,MATCH65) + mat->lchs->b;  
        /* From state MATCH75 to state UNMATCHED_QUERY */ 
        temp = LocalDnaMatchBlock_EXPL_MATRIX(mat,i-1,j-0,MATCH75) + mat->lchs->b;   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH85 to state UNMATCHED_QUERY */ 
        temp = LocalDnaMatchBlock_EXPL_MATRIX(mat,i-1,j-0,MATCH85) + mat->lchs->b;   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH95 to state UNMATCHED_QUERY */ 
        temp = LocalDnaMatchBlock_EXPL_MATRIX(mat,i-1,j-0,MATCH95) + mat->lchs->b;   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state UNMATCHED_QUERY to state UNMATCHED_QUERY */ 
        temp = LocalDnaMatchBlock_EXPL_MATRIX(mat,i-1,j-0,UNMATCHED_QUERY) + mat->lchs->u;   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for UNMATCHED_QUERY */ 
        /* Add any movement independant score and put away */ 
         LocalDnaMatchBlock_EXPL_MATRIX(mat,i,j,UNMATCHED_QUERY) = score;    


        /* Finished calculating state UNMATCHED_QUERY */ 


        /* For state UNMATCHED_TARGET */ 
        /* setting first movement to score */ 
        score = LocalDnaMatchBlock_EXPL_MATRIX(mat,i-0,j-1,UNMATCHED_QUERY) + mat->lchs->v;  
        /* From state UNMATCHED_TARGET to state UNMATCHED_TARGET */ 
        temp = LocalDnaMatchBlock_EXPL_MATRIX(mat,i-0,j-1,UNMATCHED_TARGET) + mat->lchs->u;  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for UNMATCHED_TARGET */ 
        /* Add any movement independant score and put away */ 
         LocalDnaMatchBlock_EXPL_MATRIX(mat,i,j,UNMATCHED_TARGET) = score;   


        /* Finished calculating state UNMATCHED_TARGET */ 
        }  


      /* Special state START has no special to special movements */ 


      /* Special state END has no special to special movements */ 
      }  
    stop_reporting();    
    return TRUE;     
}    


/* Function:  calculate_dpenv_LocalDnaMatchBlock(mat,dpenv)
 *
 * Descrip:    This function calculates the LocalDnaMatchBlock matrix when in explicit mode, subject to the envelope
 *
 *
 * Arg:          mat [UNKN ] LocalDnaMatchBlock which contains explicit basematrix memory [LocalDnaMatchBlock *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean calculate_dpenv_LocalDnaMatchBlock(LocalDnaMatchBlock * mat,DPEnvelope * dpenv) 
{
    int i;   
    int j;   
    int k;   
    int starti;  
    int startj;  
    int endi;    
    int endj;    
    int tot; 
    int num; 
    int should_calc; 


    if( mat->basematrix->type != BASEMATRIX_TYPE_EXPLICIT )  {  
      warn("in calculate_LocalDnaMatchBlock, passed a non Explicit matrix type, cannot calculate!"); 
      return FALSE;  
      }  


    prepare_DPEnvelope(dpenv);   
    starti = dpenv->starti;  
    if( starti < 0 ) 
      starti = 0;    
    startj = dpenv->startj;  
    if( startj < 0 ) 
      startj = 0;    
    endi = dpenv->endi;  
    if( endi > mat->leni )   
      endi = mat->leni;  
    endj = dpenv->endj;  
    if( endj > mat->lenj )   
      endj = mat->lenj;  
    tot = (endi-starti) * (endj-startj); 
    num = 0; 


    for(j=startj-1;j<endj;j++)   {  
      for(i=1;i<mat->leni;i++)   {  
        LocalDnaMatchBlock_EXPL_MATRIX(mat,i,j,MATCH65) = NEGI;  
        LocalDnaMatchBlock_EXPL_MATRIX(mat,i,j,MATCH75) = NEGI;  
        LocalDnaMatchBlock_EXPL_MATRIX(mat,i,j,MATCH85) = NEGI;  
        LocalDnaMatchBlock_EXPL_MATRIX(mat,i,j,MATCH95) = NEGI;  
        LocalDnaMatchBlock_EXPL_MATRIX(mat,i,j,UNMATCHED_QUERY) = NEGI;  
        LocalDnaMatchBlock_EXPL_MATRIX(mat,i,j,UNMATCHED_TARGET) = NEGI; 
        }  
      }  
    for(j=-1;j<mat->lenj;j++)    {  
      LocalDnaMatchBlock_EXPL_SPECIAL(mat,i,j,START) = 0;    
      LocalDnaMatchBlock_EXPL_SPECIAL(mat,i,j,END) = NEGI;   
      }  


    start_reporting("LocalDnaMatchBlock Matrix calculation: ");  
    for(j=startj;j<endj;j++) {  
      auto int score;    
      auto int temp;     
      for(i=starti;i<endi;i++)   {  
        /* Check if is in envelope - code identical to is_in_DPEnvelope, but aggressively inlined here for speed */ 
        should_calc = 0; 
        for(k=0;k<dpenv->len;k++)    {  
          auto DPUnit * u;   
          u = dpenv->dpu[k]; 
          switch(u->type)    {  
            case DPENV_RECT :    
              if( i >= u->starti && j >= u->startj && i <= (u->starti+u->height) && j <= (u->startj+u->length))  
                should_calc = 1;     
              break; 
            case DPENV_DIAG :    
              if(  abs( (i-j) - (u->starti-u->startj)) <= u->height && i+j >= u->starti+u->startj && i+j+u->length >= u->starti+u->startj)   
                should_calc = 1;     
              break; 
            }  
          if( should_calc == 1 ) 
            break;   
          }  
        if( should_calc == 0)    {  
          LocalDnaMatchBlock_EXPL_MATRIX(mat,i,j,MATCH65) = NEGI;    
          LocalDnaMatchBlock_EXPL_MATRIX(mat,i,j,MATCH75) = NEGI;    
          LocalDnaMatchBlock_EXPL_MATRIX(mat,i,j,MATCH85) = NEGI;    
          LocalDnaMatchBlock_EXPL_MATRIX(mat,i,j,MATCH95) = NEGI;    
          LocalDnaMatchBlock_EXPL_MATRIX(mat,i,j,UNMATCHED_QUERY) = NEGI;    
          LocalDnaMatchBlock_EXPL_MATRIX(mat,i,j,UNMATCHED_TARGET) = NEGI;   
          continue;  
          }  


        if( num%1000 == 0 )  
          log_full_error(REPORT,0,"[%7d] Cells %2d%%%%",num,num*100/tot);    
        num++;   


        /* For state MATCH65 */ 
        /* setting first movement to score */ 
        score = LocalDnaMatchBlock_EXPL_MATRIX(mat,i-1,j-1,MATCH65) + (mat->lchs->comp65->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->s);    
        /* From state MATCH65 to state MATCH65 */ 
        temp = LocalDnaMatchBlock_EXPL_MATRIX(mat,i-0,j-1,MATCH65) + mat->lchs->g;   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH65 to state MATCH65 */ 
        temp = LocalDnaMatchBlock_EXPL_MATRIX(mat,i-1,j-0,MATCH65) + mat->lchs->g;   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state START to state MATCH65 */ 
        temp = LocalDnaMatchBlock_EXPL_SPECIAL(mat,i-1,j-1,START) + (mat->lchs->comp65->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->v);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state UNMATCHED_TARGET to state MATCH65 */ 
        temp = LocalDnaMatchBlock_EXPL_MATRIX(mat,i-1,j-1,UNMATCHED_TARGET) + (mat->lchs->comp65->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->v);    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for MATCH65 */ 
        /* Add any movement independant score and put away */ 
         LocalDnaMatchBlock_EXPL_MATRIX(mat,i,j,MATCH65) = score;    


        /* state MATCH65 is a source for special END */ 
        temp = score + (0) + (0) ;   
        if( temp > LocalDnaMatchBlock_EXPL_SPECIAL(mat,i,j,END) )    {  
          LocalDnaMatchBlock_EXPL_SPECIAL(mat,i,j,END) = temp;   
          }  




        /* Finished calculating state MATCH65 */ 


        /* For state MATCH75 */ 
        /* setting first movement to score */ 
        score = LocalDnaMatchBlock_EXPL_MATRIX(mat,i-1,j-1,MATCH75) + (mat->lchs->comp75->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->s);    
        /* From state MATCH75 to state MATCH75 */ 
        temp = LocalDnaMatchBlock_EXPL_MATRIX(mat,i-0,j-1,MATCH75) + mat->lchs->g;   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH75 to state MATCH75 */ 
        temp = LocalDnaMatchBlock_EXPL_MATRIX(mat,i-1,j-0,MATCH75) + mat->lchs->g;   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state START to state MATCH75 */ 
        temp = LocalDnaMatchBlock_EXPL_SPECIAL(mat,i-1,j-1,START) + (mat->lchs->comp75->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->v);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state UNMATCHED_TARGET to state MATCH75 */ 
        temp = LocalDnaMatchBlock_EXPL_MATRIX(mat,i-1,j-1,UNMATCHED_TARGET) + (mat->lchs->comp75->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->v);    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for MATCH75 */ 
        /* Add any movement independant score and put away */ 
         LocalDnaMatchBlock_EXPL_MATRIX(mat,i,j,MATCH75) = score;    


        /* state MATCH75 is a source for special END */ 
        temp = score + (0) + (0) ;   
        if( temp > LocalDnaMatchBlock_EXPL_SPECIAL(mat,i,j,END) )    {  
          LocalDnaMatchBlock_EXPL_SPECIAL(mat,i,j,END) = temp;   
          }  




        /* Finished calculating state MATCH75 */ 


        /* For state MATCH85 */ 
        /* setting first movement to score */ 
        score = LocalDnaMatchBlock_EXPL_MATRIX(mat,i-1,j-1,MATCH85) + (mat->lchs->comp85->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->s);    
        /* From state MATCH85 to state MATCH85 */ 
        temp = LocalDnaMatchBlock_EXPL_MATRIX(mat,i-0,j-1,MATCH85) + mat->lchs->g;   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH85 to state MATCH85 */ 
        temp = LocalDnaMatchBlock_EXPL_MATRIX(mat,i-1,j-0,MATCH85) + mat->lchs->g;   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state START to state MATCH85 */ 
        temp = LocalDnaMatchBlock_EXPL_SPECIAL(mat,i-1,j-1,START) + (mat->lchs->comp85->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->v);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state UNMATCHED_TARGET to state MATCH85 */ 
        temp = LocalDnaMatchBlock_EXPL_MATRIX(mat,i-1,j-1,UNMATCHED_TARGET) + (mat->lchs->comp85->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->v);    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for MATCH85 */ 
        /* Add any movement independant score and put away */ 
         LocalDnaMatchBlock_EXPL_MATRIX(mat,i,j,MATCH85) = score;    


        /* state MATCH85 is a source for special END */ 
        temp = score + (0) + (0) ;   
        if( temp > LocalDnaMatchBlock_EXPL_SPECIAL(mat,i,j,END) )    {  
          LocalDnaMatchBlock_EXPL_SPECIAL(mat,i,j,END) = temp;   
          }  




        /* Finished calculating state MATCH85 */ 


        /* For state MATCH95 */ 
        /* setting first movement to score */ 
        score = LocalDnaMatchBlock_EXPL_MATRIX(mat,i-1,j-1,MATCH95) + (mat->lchs->comp95->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->s);    
        /* From state MATCH95 to state MATCH95 */ 
        temp = LocalDnaMatchBlock_EXPL_MATRIX(mat,i-0,j-1,MATCH95) + mat->lchs->g;   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH95 to state MATCH95 */ 
        temp = LocalDnaMatchBlock_EXPL_MATRIX(mat,i-1,j-0,MATCH95) + mat->lchs->g;   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state START to state MATCH95 */ 
        temp = LocalDnaMatchBlock_EXPL_SPECIAL(mat,i-1,j-1,START) + (mat->lchs->comp95->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->v);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state UNMATCHED_TARGET to state MATCH95 */ 
        temp = LocalDnaMatchBlock_EXPL_MATRIX(mat,i-1,j-1,UNMATCHED_TARGET) + (mat->lchs->comp95->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->lchs->v);    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for MATCH95 */ 
        /* Add any movement independant score and put away */ 
         LocalDnaMatchBlock_EXPL_MATRIX(mat,i,j,MATCH95) = score;    


        /* state MATCH95 is a source for special END */ 
        temp = score + (0) + (0) ;   
        if( temp > LocalDnaMatchBlock_EXPL_SPECIAL(mat,i,j,END) )    {  
          LocalDnaMatchBlock_EXPL_SPECIAL(mat,i,j,END) = temp;   
          }  




        /* Finished calculating state MATCH95 */ 


        /* For state UNMATCHED_QUERY */ 
        /* setting first movement to score */ 
        score = LocalDnaMatchBlock_EXPL_MATRIX(mat,i-1,j-0,MATCH65) + mat->lchs->b;  
        /* From state MATCH75 to state UNMATCHED_QUERY */ 
        temp = LocalDnaMatchBlock_EXPL_MATRIX(mat,i-1,j-0,MATCH75) + mat->lchs->b;   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH85 to state UNMATCHED_QUERY */ 
        temp = LocalDnaMatchBlock_EXPL_MATRIX(mat,i-1,j-0,MATCH85) + mat->lchs->b;   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH95 to state UNMATCHED_QUERY */ 
        temp = LocalDnaMatchBlock_EXPL_MATRIX(mat,i-1,j-0,MATCH95) + mat->lchs->b;   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state UNMATCHED_QUERY to state UNMATCHED_QUERY */ 
        temp = LocalDnaMatchBlock_EXPL_MATRIX(mat,i-1,j-0,UNMATCHED_QUERY) + mat->lchs->u;   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for UNMATCHED_QUERY */ 
        /* Add any movement independant score and put away */ 
         LocalDnaMatchBlock_EXPL_MATRIX(mat,i,j,UNMATCHED_QUERY) = score;    


        /* Finished calculating state UNMATCHED_QUERY */ 


        /* For state UNMATCHED_TARGET */ 
        /* setting first movement to score */ 
        score = LocalDnaMatchBlock_EXPL_MATRIX(mat,i-0,j-1,UNMATCHED_QUERY) + mat->lchs->v;  
        /* From state UNMATCHED_TARGET to state UNMATCHED_TARGET */ 
        temp = LocalDnaMatchBlock_EXPL_MATRIX(mat,i-0,j-1,UNMATCHED_TARGET) + mat->lchs->u;  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for UNMATCHED_TARGET */ 
        /* Add any movement independant score and put away */ 
         LocalDnaMatchBlock_EXPL_MATRIX(mat,i,j,UNMATCHED_TARGET) = score;   


        /* Finished calculating state UNMATCHED_TARGET */ 
        }  


      /* Special state START has no special to special movements */ 


      /* Special state END has no special to special movements */ 
      }  
    stop_reporting();    
    return TRUE;     
}    


/* Function:  LocalDnaMatchBlock_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [LocalDnaMatchBlock *]
 *
 */
LocalDnaMatchBlock * LocalDnaMatchBlock_alloc(void) 
{
    LocalDnaMatchBlock * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(LocalDnaMatchBlock *) ckalloc (sizeof(LocalDnaMatchBlock))) == NULL)    {  
      warn("LocalDnaMatchBlock_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->basematrix = NULL;  
    out->shatter = NULL; 
    out->leni = 0;   
    out->lenj = 0;   


    return out;  
}    


/* Function:  free_LocalDnaMatchBlock(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [LocalDnaMatchBlock *]
 *
 * Return [UNKN ]  Undocumented return value [LocalDnaMatchBlock *]
 *
 */
LocalDnaMatchBlock * free_LocalDnaMatchBlock(LocalDnaMatchBlock * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a LocalDnaMatchBlock obj. Should be trappable");    
      return NULL;   
      }  


#ifdef PTHREAD   
    assert(pthread_mutex_lock(&(obj->dynamite_mutex)) == 0); 
#endif   
    if( obj->dynamite_hard_link > 1)     {  
      return_early = 1;  
      obj->dynamite_hard_link--; 
      }  
#ifdef PTHREAD   
    assert(pthread_mutex_unlock(&(obj->dynamite_mutex)) == 0);   
#endif   
    if( return_early == 1)   
      return NULL;   
    if( obj->basematrix != NULL) 
      free_BaseMatrix(obj->basematrix);  
    if( obj->shatter != NULL)    
      free_ShatterMatrix(obj->shatter);  
    /* obj->query is linked in */ 
    /* obj->target is linked in */ 
    /* obj->lchs is linked in */ 


    ckfree(obj); 
    return NULL; 
}    





#ifdef _cplusplus
}
#endif
