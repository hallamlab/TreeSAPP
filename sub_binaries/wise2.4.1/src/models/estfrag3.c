#ifdef _cplusplus
extern "C" {
#endif
#include "estfrag3.h"

# line 5 "estfrag3.c"


  /*****************   C functions  ****************/
  /*             Written using dynamite            */
  /*            Sat Sep  8 09:05:31 2007           */
  /*            email birney@sanger.ac.uk          */
  /* http://www.sanger.ac.uk/Users/birney/dynamite */
  /*************************************************/


  /* Please report any problems or bugs to         */
  /* Ewan Birney, birney@sanger.ac.uk              */


/* basic set of macros to map states to numbers */ 
#define MATCH 0  
#define INSERT 1 
#define DELETE 2 


#define START 0  
#define END 1    


#define EstFrag3_EXPL_MATRIX(this_matrix,i,j,STATE) this_matrix->basematrix->matrix[((j+5)*3)+STATE][i+1]    
#define EstFrag3_EXPL_SPECIAL(matrix,i,j,STATE) matrix->basematrix->specmatrix[STATE][j+5]   
#define EstFrag3_READ_OFF_ERROR -7
  


#define EstFrag3_VSMALL_MATRIX(mat,i,j,STATE) internal_matrix[(j+6)%6][((i+1)*3)+STATE]  
#define EstFrag3_VSMALL_SPECIAL(mat,i,j,STATE) internal_specials[(j+6)%6][STATE] 




#define EstFrag3_SHATTER_SPECIAL(matrix,i,j,STATE) matrix->shatter->special[STATE][j]    
#define EstFrag3_SHATTER_MATRIX(matrix,i,j,STATE)  fetch_cell_value_ShatterMatrix(mat->shatter,i,j,STATE)    


/* Function:  PackAln_read_Shatter_EstFrag3(mat)
 *
 * Descrip:    Reads off PackAln from shatter matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [EstFrag3 *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * PackAln_read_Shatter_EstFrag3(EstFrag3 * mat) 
{
    EstFrag3_access_func_holder holder;  


    holder.access_main    = EstFrag3_shatter_access_main;    
    holder.access_special = EstFrag3_shatter_access_special; 
    assert(mat);     
    assert(mat->shatter);    
    return PackAln_read_generic_EstFrag3(mat,holder);    
}    


/* Function:  EstFrag3_shatter_access_main(mat,i,j,state)
 *
 * Descrip: No Description
 *
 * Arg:          mat [UNKN ] Undocumented argument [EstFrag3 *]
 * Arg:            i [UNKN ] Undocumented argument [int]
 * Arg:            j [UNKN ] Undocumented argument [int]
 * Arg:        state [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int EstFrag3_shatter_access_main(EstFrag3 * mat,int i,int j,int state) 
{
    return EstFrag3_SHATTER_MATRIX(mat,i,j,state);   
}    


/* Function:  EstFrag3_shatter_access_special(mat,i,j,state)
 *
 * Descrip: No Description
 *
 * Arg:          mat [UNKN ] Undocumented argument [EstFrag3 *]
 * Arg:            i [UNKN ] Undocumented argument [int]
 * Arg:            j [UNKN ] Undocumented argument [int]
 * Arg:        state [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int EstFrag3_shatter_access_special(EstFrag3 * mat,int i,int j,int state) 
{
    return EstFrag3_SHATTER_SPECIAL(mat,i,j,state);  
}    


/* Function:  calculate_shatter_EstFrag3(mat,dpenv)
 *
 * Descrip:    This function calculates the EstFrag3 matrix when in shatter mode
 *
 *
 * Arg:          mat [UNKN ] (null) [EstFrag3 *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean calculate_shatter_EstFrag3(EstFrag3 * mat,DPEnvelope * dpenv) 
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
    int * SIG_1_3;   
    int * SIG_1_2;   
    int * SIG_1_1;   
    int * SIG_1_4;   
    int * SIG_1_5;   
    int * SIG_0_3;   
    int * SIG_0_2;   
    int * SIG_0_1;   
    int * SIG_0_4;   
    int * SIG_0_5;   
    int * SIG_1_0;   


    leni = mat->leni;    
    lenj = mat->lenj;    


    mat->shatter = new_ShatterMatrix(dpenv,3,lenj,2);    
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


    start_reporting("EstFrag3 Matrix calculation: ");    
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
        SIG_1_3 = fetch_cell_from_ShatterMatrix(mat->shatter,i-1,j-3);   
        SIG_1_2 = fetch_cell_from_ShatterMatrix(mat->shatter,i-1,j-2);   
        SIG_1_1 = fetch_cell_from_ShatterMatrix(mat->shatter,i-1,j-1);   
        SIG_1_4 = fetch_cell_from_ShatterMatrix(mat->shatter,i-1,j-4);   
        SIG_1_5 = fetch_cell_from_ShatterMatrix(mat->shatter,i-1,j-5);   
        SIG_0_3 = fetch_cell_from_ShatterMatrix(mat->shatter,i-0,j-3);   
        SIG_0_2 = fetch_cell_from_ShatterMatrix(mat->shatter,i-0,j-2);   
        SIG_0_1 = fetch_cell_from_ShatterMatrix(mat->shatter,i-0,j-1);   
        SIG_0_4 = fetch_cell_from_ShatterMatrix(mat->shatter,i-0,j-4);   
        SIG_0_5 = fetch_cell_from_ShatterMatrix(mat->shatter,i-0,j-5);   
        SIG_1_0 = fetch_cell_from_ShatterMatrix(mat->shatter,i-1,j-0);   




        /* For state MATCH */ 
        /* setting first movement to score */ 
        score = SIG_1_3[MATCH] + (mat->query->seg[i]->transition[GW_MATCH2MATCH]+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)]);     
        /* From state INSERT to state MATCH */ 
        temp = SIG_1_3[INSERT] + (mat->query->seg[i]->transition[GW_INSERT2MATCH]+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)]);    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state MATCH */ 
        temp = SIG_1_3[DELETE] + (mat->query->seg[i]->transition[GW_DELETE2MATCH]+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)]);    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state START to state MATCH */ 
        temp = EstFrag3_SHATTER_SPECIAL(mat,i-1,j-3,START) + (mat->query->seg[i]->transition[GW_START2MATCH]+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)]);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* Has restricted position */ 
        if( (j-3) == 0  )    {  
          /* From state START to state MATCH */ 
          temp = EstFrag3_SHATTER_SPECIAL(mat,i-1,j-3,START) + (mat->start_frag+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)]);  
          if( temp  > score )    {  
            score = temp;    
            }  
          }  
        /* From state MATCH to state MATCH */ 
        temp = SIG_1_2[MATCH] + mat->cp->trans[PCD_DELETE_1_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = SIG_1_1[MATCH] + mat->cp->trans[PCD_DELETE_2_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = SIG_1_4[MATCH] + mat->cp->trans[PCD_INSERT_1_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = SIG_1_5[MATCH] + mat->cp->trans[PCD_INSERT_2_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state MATCH */ 
        temp = SIG_1_2[INSERT] + mat->cp->trans[PCD_DELETE_1_BASE];  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state MATCH */ 
        temp = SIG_1_1[INSERT] + mat->cp->trans[PCD_DELETE_2_BASE];  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state MATCH */ 
        temp = SIG_1_4[INSERT] + mat->cp->trans[PCD_INSERT_1_BASE];  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state MATCH */ 
        temp = SIG_1_5[INSERT] + mat->cp->trans[PCD_INSERT_2_BASE];  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state MATCH */ 
        temp = SIG_1_2[DELETE] + mat->cp->trans[PCD_DELETE_1_BASE];  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state MATCH */ 
        temp = SIG_1_1[DELETE] + mat->cp->trans[PCD_DELETE_2_BASE];  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state MATCH */ 
        temp = SIG_1_4[DELETE] + mat->cp->trans[PCD_INSERT_1_BASE];  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state MATCH */ 
        temp = SIG_1_5[DELETE] + mat->cp->trans[PCD_INSERT_2_BASE];  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for MATCH */ 
        /* Add any movement independant score and put away */ 
         SIG_0_0[MATCH] = score; 


        /* state MATCH is a source for special END */ 
        temp = score + (mat->query->seg[i]->transition[GW_MATCH2END]) + (0) ;    
        if( temp > EstFrag3_SHATTER_SPECIAL(mat,i,j,END) )   {  
          EstFrag3_SHATTER_SPECIAL(mat,i,j,END) = temp;  
          }  




        /* state MATCH is a source for special END */ 
        /* Has restricted position */ 
        if( j == mat->lenj-1 )   {  
          temp = score + (mat->end_frag) + (0) ;     
          if( temp > EstFrag3_SHATTER_SPECIAL(mat,i,j,END) )     {  
            EstFrag3_SHATTER_SPECIAL(mat,i,j,END) = temp;    
            }  


          }  


        /* Finished calculating state MATCH */ 


        /* For state INSERT */ 
        /* setting first movement to score */ 
        score = SIG_0_3[MATCH] + (mat->query->seg[i]->transition[GW_MATCH2INSERT]+mat->query->seg[i]->insert[CSEQ_CDNA_CODON(mat->target,j)]);   
        /* From state INSERT to state INSERT */ 
        temp = SIG_0_3[INSERT] + (mat->query->seg[i]->transition[GW_INSERT2INSERT]+mat->query->seg[i]->insert[CSEQ_CDNA_CODON(mat->target,j)]);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state INSERT */ 
        temp = SIG_0_3[DELETE] + (mat->query->seg[i]->transition[GW_DELETE2INSERT]+mat->query->seg[i]->insert[CSEQ_CDNA_CODON(mat->target,j)]);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state START to state INSERT */ 
        temp = EstFrag3_SHATTER_SPECIAL(mat,i-0,j-3,START) + (mat->query->seg[i]->transition[GW_START2INSERT]+mat->query->seg[i]->insert[CSEQ_CDNA_CODON(mat->target,j)]);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state INSERT */ 
        temp = SIG_0_2[INSERT] + mat->cp->trans[PCD_DELETE_1_BASE];  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state INSERT */ 
        temp = SIG_0_1[INSERT] + mat->cp->trans[PCD_DELETE_2_BASE];  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state INSERT */ 
        temp = SIG_0_4[INSERT] + mat->cp->trans[PCD_INSERT_1_BASE];  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state INSERT */ 
        temp = SIG_0_5[INSERT] + mat->cp->trans[PCD_INSERT_2_BASE];  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state INSERT */ 
        temp = SIG_0_2[MATCH] + mat->cp->trans[PCD_DELETE_1_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state INSERT */ 
        temp = SIG_0_1[MATCH] + mat->cp->trans[PCD_DELETE_2_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state INSERT */ 
        temp = SIG_0_4[MATCH] + mat->cp->trans[PCD_INSERT_1_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state INSERT */ 
        temp = SIG_0_5[MATCH] + mat->cp->trans[PCD_INSERT_2_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state INSERT */ 
        temp = SIG_0_2[DELETE] + mat->cp->trans[PCD_DELETE_1_BASE];  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state INSERT */ 
        temp = SIG_0_1[DELETE] + mat->cp->trans[PCD_DELETE_2_BASE];  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state INSERT */ 
        temp = SIG_0_4[DELETE] + mat->cp->trans[PCD_INSERT_1_BASE];  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state INSERT */ 
        temp = SIG_0_5[DELETE] + mat->cp->trans[PCD_INSERT_2_BASE];  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for INSERT */ 
        /* Add any movement independant score and put away */ 
         SIG_0_0[INSERT] = score;    


        /* state INSERT is a source for special END */ 
        temp = score + (mat->query->seg[i]->transition[GW_INSERT2END]) + (0) ;   
        if( temp > EstFrag3_SHATTER_SPECIAL(mat,i,j,END) )   {  
          EstFrag3_SHATTER_SPECIAL(mat,i,j,END) = temp;  
          }  




        /* Finished calculating state INSERT */ 


        /* For state DELETE */ 
        /* setting first movement to score */ 
        score = SIG_1_0[MATCH] + mat->query->seg[i]->transition[GW_MATCH2DELETE];    
        /* From state INSERT to state DELETE */ 
        temp = SIG_1_0[INSERT] + mat->query->seg[i]->transition[GW_INSERT2DELETE];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state DELETE */ 
        temp = SIG_1_0[DELETE] + mat->query->seg[i]->transition[GW_DELETE2DELETE];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state START to state DELETE */ 
        temp = EstFrag3_SHATTER_SPECIAL(mat,i-1,j-0,START) + mat->query->seg[i]->transition[GW_START2DELETE];    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for DELETE */ 
        /* Add any movement independant score and put away */ 
         SIG_0_0[DELETE] = score;    


        /* state DELETE is a source for special END */ 
        temp = score + (mat->query->seg[i]->transition[GW_DELETE2END]) + (0) ;   
        if( temp > EstFrag3_SHATTER_SPECIAL(mat,i,j,END) )   {  
          EstFrag3_SHATTER_SPECIAL(mat,i,j,END) = temp;  
          }  




        /* Finished calculating state DELETE */ 
        }  


      /* Special state START has no special to special movements */ 


      /* Special state END has no special to special movements */ 
      }  
    stop_reporting();    
    return TRUE;     
}    


/* Function:  search_EstFrag3(dbsi,out,querydb,targetdb,cp,start_frag,end_frag)
 *
 * Descrip:    This function makes a database search of EstFrag3
 *             It uses the dbsi structure to choose which implementation
 *             to use of the database searching. This way at run time you
 *             can switch between single threaded/multi-threaded or hardware
 *
 *
 * Arg:              dbsi [UNKN ] Undocumented argument [DBSearchImpl *]
 * Arg:               out [UNKN ] Undocumented argument [Hscore *]
 * Arg:           querydb [UNKN ] Undocumented argument [GeneWiseDB*]
 * Arg:          targetdb [UNKN ] Undocumented argument [cDNADB*]
 * Arg:                cp [UNKN ] Undocumented argument [cDNAParserScore *]
 * Arg:        start_frag [UNKN ] Undocumented argument [Score]
 * Arg:          end_frag [UNKN ] Undocumented argument [Score]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type search_EstFrag3(DBSearchImpl * dbsi,Hscore * out,GeneWiseDB* querydb,cDNADB* targetdb ,cDNAParserScore * cp,Score start_frag,Score end_frag) 
{
    if( out == NULL )    {  
      warn("Passed in a null Hscore object into search_EstFrag3. Can't process results!");   
      return SEARCH_ERROR;   
      }  
    if( dbsi == NULL )   {  
      warn("Passed in a null DBSearchImpl object into search_EstFrag3. Can't process results!"); 
      return SEARCH_ERROR;   
      }  
    if( dbsi->trace_level > 5 )  
      warn("Asking for trace level of %d in database search for EstFrag3, but it was compiled with a trace level of -2139062144. Not all trace statements can be shown",dbsi->trace_level);  
    switch(dbsi->type)   { /*switch on implementation*/ 
      case DBSearchImpl_Serial : 
        return serial_search_EstFrag3(out,querydb, targetdb ,cp,start_frag,end_frag);    
      case DBSearchImpl_Pthreads :   
        warn("This matrix EstFrag3 was not dyc compiled with thread support");   
        return SEARCH_ERROR; 
      default :  
        warn("database search implementation %s was not provided in the compiled dynamite file from EstFrag3",impl_string_DBSearchImpl(dbsi));   
        return SEARCH_ERROR; 
      } /* end of switch on implementation */ 


}    


/* Function:  serial_search_EstFrag3(out,querydb,targetdb,cp,start_frag,end_frag)
 *
 * Descrip:    This function makes a database search of EstFrag3
 *             It is a single processor implementation
 *
 *
 * Arg:               out [UNKN ] Undocumented argument [Hscore *]
 * Arg:           querydb [UNKN ] Undocumented argument [GeneWiseDB*]
 * Arg:          targetdb [UNKN ] Undocumented argument [cDNADB*]
 * Arg:                cp [UNKN ] Undocumented argument [cDNAParserScore *]
 * Arg:        start_frag [UNKN ] Undocumented argument [Score]
 * Arg:          end_frag [UNKN ] Undocumented argument [Score]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type serial_search_EstFrag3(Hscore * out,GeneWiseDB* querydb,cDNADB* targetdb ,cDNAParserScore * cp,Score start_frag,Score end_frag) 
{
    GeneWiseScore* query;    
    ComplexSequence* target;     
    int db_status;   
    int score;   
    int query_pos = 0;   
    int target_pos = 0;  
    DataScore * ds;  


    push_errormsg_stack("Before any actual search in db searching"); 
    query = init_GeneWiseDB(querydb,&db_status); 
    if( db_status == DB_RETURN_ERROR )   {  
      warn("In searching EstFrag3, got a database reload error on the query [query] database");  
      return SEARCH_ERROR;   
      }  
    for(;;)  { /*For all query entries*/ 


      target_pos = 0;    


      target = init_cDNADB(targetdb,&db_status); 
      if( db_status == DB_RETURN_ERROR )     {  
        warn("In searching EstFrag3, got a database init error on the target [target] database");    
        return SEARCH_ERROR; 
        }  
      for(;;)    { /*For all target entries*/ 


        /* Maximum length to search - should check */ 
        if( query->len > 1500 )  { /*if over length*/ 
          warn("A query over the length when maxlen provided. Problem!");    
          score = -10000;
   
          } /* end of if over length */ 
        else 
          score = score_only_EstFrag3(query, target , cp, start_frag, end_frag);     
        if( should_store_Hscore(out,score) == TRUE )     { /*if storing datascore*/ 
          ds = new_DataScore_from_storage(out);  
          if( ds == NULL )   {  
            warn("EstFrag3 search had a memory error in allocating a new_DataScore (?a leak somewhere - DataScore is a very small datastructure");   
            return SEARCH_ERROR; 
            }  
          /* Now: add query/target information to the entry */ 
          dataentry_add_GeneWiseDB(ds->query,query,querydb);     
          dataentry_add_cDNADB(ds->target,target,targetdb);  
          ds->score = score;     
          add_Hscore(out,ds);    
          } /* end of if storing datascore */ 
        pop_errormsg_stack();    
        push_errormsg_stack("DB searching: just finished [Query Pos: %d] [Target Pos: %d]",query_pos,target_pos);    


         target = reload_cDNADB(target,targetdb,&db_status);     
        if( db_status == DB_RETURN_ERROR )   {  
          warn("In searching EstFrag3, Reload error on database target, position %d,%d",query_pos,target_pos);   
          return SEARCH_ERROR;   
          }  
        if( db_status == DB_RETURN_END ) 
          break;/* Out of target loop */ 
        target_pos++;    
        } /* end of For all target entries */ 
      close_cDNADB(target,targetdb);     
       query = reload_GeneWiseDB(query,querydb,&db_status);  
      if( db_status == DB_RETURN_ERROR)  {  
        warn("In searching EstFrag3, Reload error on database query, position %d,%d",query_pos,target_pos);  
        return SEARCH_ERROR; 
        }  
      if( db_status == DB_RETURN_END)    
        break;  /* Out of query loop */ 
      query_pos++;   
      } /* end of For all query entries */ 
    close_GeneWiseDB(query,querydb);     
    pop_errormsg_stack();    
    return SEARCH_OK;    
}    


/* Function:  score_only_EstFrag3(query,target,cp,start_frag,end_frag)
 *
 * Descrip:    This function just calculates the score for the matrix
 *             I am pretty sure we can do this better, but hey, for the moment...
 *             It calls /allocate_EstFrag3_only
 *
 *
 * Arg:             query [UNKN ] query data structure [GeneWiseScore*]
 * Arg:            target [UNKN ] target data structure [ComplexSequence*]
 * Arg:                cp [UNKN ] Resource [cDNAParserScore *]
 * Arg:        start_frag [UNKN ] Resource [Score]
 * Arg:          end_frag [UNKN ] Resource [Score]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int score_only_EstFrag3(GeneWiseScore* query,ComplexSequence* target ,cDNAParserScore * cp,Score start_frag,Score end_frag) 
{
    int bestscore = NEGI;    
    int i;   
    int j;   
    int k;   
    EstFrag3 * mat;  
    int internal_matrix[6][(1500+1) * 3];    
    int internal_specials[6][2]; 


    mat = allocate_EstFrag3_only(query, target , cp, start_frag, end_frag);  
    if( mat == NULL )    {  
      warn("Memory allocation error in the db search - unable to communicate to calling function. this spells DIASTER!");    
      return NEGI;   
      }  
    /* Ok,don't need to allocate matrix as it is internal, because we have a max length */ 


    /* Now, initiate matrix */ 
    for(j=0;j<7;j++) {  
      for(i=(-1);i<mat->leni;i++)    {  
        for(k=0;k<3;k++) 
          EstFrag3_VSMALL_MATRIX(mat,i,j,k) = NEGI;  
        }  
      EstFrag3_VSMALL_SPECIAL(mat,i,j,START) = 0;    
      EstFrag3_VSMALL_SPECIAL(mat,i,j,END) = NEGI;   
      }  


    /* Ok, lets do-o-o-o-o it */ 


    for(j=0;j<mat->lenj;j++) { /*for all target positions*/ 
      auto int score;    
      auto int temp;     
      for(i=0;i<mat->leni;i++)   { /*for all query positions*/ 


        /* For state MATCH */ 
        /* setting first movement to score */ 
        score = EstFrag3_VSMALL_MATRIX(mat,i-1,j-3,MATCH) + (mat->query->seg[i]->transition[GW_MATCH2MATCH]+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)]);  
        /* From state INSERT to state MATCH */ 
        temp = EstFrag3_VSMALL_MATRIX(mat,i-1,j-3,INSERT) + (mat->query->seg[i]->transition[GW_INSERT2MATCH]+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)]);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state MATCH */ 
        temp = EstFrag3_VSMALL_MATRIX(mat,i-1,j-3,DELETE) + (mat->query->seg[i]->transition[GW_DELETE2MATCH]+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)]);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state START to state MATCH */ 
        temp = EstFrag3_VSMALL_SPECIAL(mat,i-1,j-3,START) + (mat->query->seg[i]->transition[GW_START2MATCH]+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)]);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* Has restricted position */ 
        if( (j-3) == 0  )    {  
          /* From state START to state MATCH */ 
          temp = EstFrag3_VSMALL_SPECIAL(mat,i-1,j-3,START) + (mat->start_frag+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)]);   
          if( temp  > score )    {  
            score = temp;    
            }  
          }  
        /* From state MATCH to state MATCH */ 
        temp = EstFrag3_VSMALL_MATRIX(mat,i-1,j-2,MATCH) + mat->cp->trans[PCD_DELETE_1_BASE];    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = EstFrag3_VSMALL_MATRIX(mat,i-1,j-1,MATCH) + mat->cp->trans[PCD_DELETE_2_BASE];    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = EstFrag3_VSMALL_MATRIX(mat,i-1,j-4,MATCH) + mat->cp->trans[PCD_INSERT_1_BASE];    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = EstFrag3_VSMALL_MATRIX(mat,i-1,j-5,MATCH) + mat->cp->trans[PCD_INSERT_2_BASE];    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state MATCH */ 
        temp = EstFrag3_VSMALL_MATRIX(mat,i-1,j-2,INSERT) + mat->cp->trans[PCD_DELETE_1_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state MATCH */ 
        temp = EstFrag3_VSMALL_MATRIX(mat,i-1,j-1,INSERT) + mat->cp->trans[PCD_DELETE_2_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state MATCH */ 
        temp = EstFrag3_VSMALL_MATRIX(mat,i-1,j-4,INSERT) + mat->cp->trans[PCD_INSERT_1_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state MATCH */ 
        temp = EstFrag3_VSMALL_MATRIX(mat,i-1,j-5,INSERT) + mat->cp->trans[PCD_INSERT_2_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state MATCH */ 
        temp = EstFrag3_VSMALL_MATRIX(mat,i-1,j-2,DELETE) + mat->cp->trans[PCD_DELETE_1_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state MATCH */ 
        temp = EstFrag3_VSMALL_MATRIX(mat,i-1,j-1,DELETE) + mat->cp->trans[PCD_DELETE_2_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state MATCH */ 
        temp = EstFrag3_VSMALL_MATRIX(mat,i-1,j-4,DELETE) + mat->cp->trans[PCD_INSERT_1_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state MATCH */ 
        temp = EstFrag3_VSMALL_MATRIX(mat,i-1,j-5,DELETE) + mat->cp->trans[PCD_INSERT_2_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for MATCH */ 
        /* Add any movement independant score and put away */ 
         EstFrag3_VSMALL_MATRIX(mat,i,j,MATCH) = score;  


        /* state MATCH is a source for special END */ 
        temp = score + (mat->query->seg[i]->transition[GW_MATCH2END]) + (0) ;    
        if( temp > EstFrag3_VSMALL_SPECIAL(mat,i,j,END) )    {  
          EstFrag3_VSMALL_SPECIAL(mat,i,j,END) = temp;   
          }  




        /* state MATCH is a source for special END */ 
        /* Has restricted position */ 
        if( j == mat->lenj-1 )   {  
          temp = score + (mat->end_frag) + (0) ;     
          if( temp > EstFrag3_VSMALL_SPECIAL(mat,i,j,END) )  {  
            EstFrag3_VSMALL_SPECIAL(mat,i,j,END) = temp;     
            }  


          }  


        /* Finished calculating state MATCH */ 


        /* For state INSERT */ 
        /* setting first movement to score */ 
        score = EstFrag3_VSMALL_MATRIX(mat,i-0,j-3,MATCH) + (mat->query->seg[i]->transition[GW_MATCH2INSERT]+mat->query->seg[i]->insert[CSEQ_CDNA_CODON(mat->target,j)]);    
        /* From state INSERT to state INSERT */ 
        temp = EstFrag3_VSMALL_MATRIX(mat,i-0,j-3,INSERT) + (mat->query->seg[i]->transition[GW_INSERT2INSERT]+mat->query->seg[i]->insert[CSEQ_CDNA_CODON(mat->target,j)]);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state INSERT */ 
        temp = EstFrag3_VSMALL_MATRIX(mat,i-0,j-3,DELETE) + (mat->query->seg[i]->transition[GW_DELETE2INSERT]+mat->query->seg[i]->insert[CSEQ_CDNA_CODON(mat->target,j)]);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state START to state INSERT */ 
        temp = EstFrag3_VSMALL_SPECIAL(mat,i-0,j-3,START) + (mat->query->seg[i]->transition[GW_START2INSERT]+mat->query->seg[i]->insert[CSEQ_CDNA_CODON(mat->target,j)]);    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state INSERT */ 
        temp = EstFrag3_VSMALL_MATRIX(mat,i-0,j-2,INSERT) + mat->cp->trans[PCD_DELETE_1_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state INSERT */ 
        temp = EstFrag3_VSMALL_MATRIX(mat,i-0,j-1,INSERT) + mat->cp->trans[PCD_DELETE_2_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state INSERT */ 
        temp = EstFrag3_VSMALL_MATRIX(mat,i-0,j-4,INSERT) + mat->cp->trans[PCD_INSERT_1_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state INSERT */ 
        temp = EstFrag3_VSMALL_MATRIX(mat,i-0,j-5,INSERT) + mat->cp->trans[PCD_INSERT_2_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state INSERT */ 
        temp = EstFrag3_VSMALL_MATRIX(mat,i-0,j-2,MATCH) + mat->cp->trans[PCD_DELETE_1_BASE];    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state INSERT */ 
        temp = EstFrag3_VSMALL_MATRIX(mat,i-0,j-1,MATCH) + mat->cp->trans[PCD_DELETE_2_BASE];    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state INSERT */ 
        temp = EstFrag3_VSMALL_MATRIX(mat,i-0,j-4,MATCH) + mat->cp->trans[PCD_INSERT_1_BASE];    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state INSERT */ 
        temp = EstFrag3_VSMALL_MATRIX(mat,i-0,j-5,MATCH) + mat->cp->trans[PCD_INSERT_2_BASE];    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state INSERT */ 
        temp = EstFrag3_VSMALL_MATRIX(mat,i-0,j-2,DELETE) + mat->cp->trans[PCD_DELETE_1_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state INSERT */ 
        temp = EstFrag3_VSMALL_MATRIX(mat,i-0,j-1,DELETE) + mat->cp->trans[PCD_DELETE_2_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state INSERT */ 
        temp = EstFrag3_VSMALL_MATRIX(mat,i-0,j-4,DELETE) + mat->cp->trans[PCD_INSERT_1_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state INSERT */ 
        temp = EstFrag3_VSMALL_MATRIX(mat,i-0,j-5,DELETE) + mat->cp->trans[PCD_INSERT_2_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for INSERT */ 
        /* Add any movement independant score and put away */ 
         EstFrag3_VSMALL_MATRIX(mat,i,j,INSERT) = score; 


        /* state INSERT is a source for special END */ 
        temp = score + (mat->query->seg[i]->transition[GW_INSERT2END]) + (0) ;   
        if( temp > EstFrag3_VSMALL_SPECIAL(mat,i,j,END) )    {  
          EstFrag3_VSMALL_SPECIAL(mat,i,j,END) = temp;   
          }  




        /* Finished calculating state INSERT */ 


        /* For state DELETE */ 
        /* setting first movement to score */ 
        score = EstFrag3_VSMALL_MATRIX(mat,i-1,j-0,MATCH) + mat->query->seg[i]->transition[GW_MATCH2DELETE];     
        /* From state INSERT to state DELETE */ 
        temp = EstFrag3_VSMALL_MATRIX(mat,i-1,j-0,INSERT) + mat->query->seg[i]->transition[GW_INSERT2DELETE];    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state DELETE */ 
        temp = EstFrag3_VSMALL_MATRIX(mat,i-1,j-0,DELETE) + mat->query->seg[i]->transition[GW_DELETE2DELETE];    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state START to state DELETE */ 
        temp = EstFrag3_VSMALL_SPECIAL(mat,i-1,j-0,START) + mat->query->seg[i]->transition[GW_START2DELETE];     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for DELETE */ 
        /* Add any movement independant score and put away */ 
         EstFrag3_VSMALL_MATRIX(mat,i,j,DELETE) = score; 


        /* state DELETE is a source for special END */ 
        temp = score + (mat->query->seg[i]->transition[GW_DELETE2END]) + (0) ;   
        if( temp > EstFrag3_VSMALL_SPECIAL(mat,i,j,END) )    {  
          EstFrag3_VSMALL_SPECIAL(mat,i,j,END) = temp;   
          }  




        /* Finished calculating state DELETE */ 
        } /* end of for all query positions */ 




      /* Special state START has no special to special movements */ 


      /* Special state END has no special to special movements */ 
      if( bestscore < EstFrag3_VSMALL_SPECIAL(mat,0,j,END) ) 
        bestscore = EstFrag3_VSMALL_SPECIAL(mat,0,j,END);    
      } /* end of for all target positions */ 


    mat = free_EstFrag3(mat);    
    return bestscore;    
}    


/* Function:  PackAln_bestmemory_EstFrag3(query,target,cp,start_frag,end_frag,dpenv,dpri)
 *
 * Descrip:    This function chooses the best memory set-up for the alignment
 *             using calls to basematrix, and then implements either a large
 *             or small memory model.
 *
 *             It is the best function to use if you just want an alignment
 *
 *             If you want a label alignment, you will need
 *             /convert_PackAln_to_AlnBlock_EstFrag3
 *
 *
 * Arg:             query [UNKN ] query data structure [GeneWiseScore*]
 * Arg:            target [UNKN ] target data structure [ComplexSequence*]
 * Arg:                cp [UNKN ] Resource [cDNAParserScore *]
 * Arg:        start_frag [UNKN ] Resource [Score]
 * Arg:          end_frag [UNKN ] Resource [Score]
 * Arg:             dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:              dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * PackAln_bestmemory_EstFrag3(GeneWiseScore* query,ComplexSequence* target ,cDNAParserScore * cp,Score start_frag,Score end_frag,DPEnvelope * dpenv,DPRunImpl * dpri) 
{
    long long total; 
    EstFrag3 * mat;  
    PackAln * out;   
    DebugMatrix * de;    
    DPRunImplMemory strategy;    
    assert(dpri);    


    total = query->len * target->seq->len;   
    if( dpri->memory == DPIM_Default )   {  
      if( (total * 3 * sizeof(int)) > 1000*dpri->kbyte_size) {  
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
        if( (mat=allocate_Expl_EstFrag3(query, target , cp, start_frag, end_frag,dpri)) == NULL )    {  
          warn("Unable to allocate large EstFrag3 version"); 
          return NULL;   
          }  
        calculate_dpenv_EstFrag3(mat,dpenv);     
        out =  PackAln_read_Expl_EstFrag3(mat);  
        }  
      else   {  
        mat = allocate_EstFrag3_only(query, target , cp, start_frag, end_frag);  
        calculate_shatter_EstFrag3(mat,dpenv);   
        out = PackAln_read_Shatter_EstFrag3(mat);    
        }  
      }  
    else {  
      if( strategy == DPIM_Linear )  {  
        /* use small implementation */ 
        if( (mat=allocate_Small_EstFrag3(query, target , cp, start_frag, end_frag)) == NULL )    {  
          warn("Unable to allocate small EstFrag3 version"); 
          return NULL;   
          }  
        out = PackAln_calculate_Small_EstFrag3(mat,dpenv);   
        }  
      else   {  
        /* use Large implementation */ 
        if( (mat=allocate_Expl_EstFrag3(query, target , cp, start_frag, end_frag,dpri)) == NULL )    {  
          warn("Unable to allocate large EstFrag3 version"); 
          return NULL;   
          }  
        if( dpri->debug == TRUE) {  
          fatal("Asked for dydebug, but dynamite file not compiled with -g. Need to recompile dynamite source"); 
          }  
        else {  
          calculate_EstFrag3(mat);   
          out =  PackAln_read_Expl_EstFrag3(mat);    
          }  
        }  
      }  


    mat = free_EstFrag3(mat);    
    return out;  
}    


/* Function:  allocate_EstFrag3_only(query,target,cp,start_frag,end_frag)
 *
 * Descrip:    This function only allocates the EstFrag3 structure
 *             checks types where possible and determines leni and lenj
 *             The basematrix area is delt with elsewhere
 *
 *
 * Arg:             query [UNKN ] query data structure [GeneWiseScore*]
 * Arg:            target [UNKN ] target data structure [ComplexSequence*]
 * Arg:                cp [UNKN ] Resource [cDNAParserScore *]
 * Arg:        start_frag [UNKN ] Resource [Score]
 * Arg:          end_frag [UNKN ] Resource [Score]
 *
 * Return [UNKN ]  Undocumented return value [EstFrag3 *]
 *
 */
EstFrag3 * allocate_EstFrag3_only(GeneWiseScore* query,ComplexSequence* target ,cDNAParserScore * cp,Score start_frag,Score end_frag) 
{
    EstFrag3 * out;  


    if((out= EstFrag3_alloc()) == NULL)  {  
      warn("Allocation of basic EstFrag3 structure failed...");  
      return NULL;   
      }  


    out->query = query;  
    out->target = target;    
    out->cp = cp;    
    out->start_frag = start_frag;    
    out->end_frag = end_frag;    
    out->leni = query->len;  
    out->lenj = target->seq->len;    
    return out;  
}    


/* Function:  allocate_Expl_EstFrag3(query,target,cp,start_frag,end_frag,dpri)
 *
 * Descrip:    This function allocates the EstFrag3 structure
 *             and the basematrix area for explicit memory implementations
 *             It calls /allocate_EstFrag3_only
 *
 *
 * Arg:             query [UNKN ] query data structure [GeneWiseScore*]
 * Arg:            target [UNKN ] target data structure [ComplexSequence*]
 * Arg:                cp [UNKN ] Resource [cDNAParserScore *]
 * Arg:        start_frag [UNKN ] Resource [Score]
 * Arg:          end_frag [UNKN ] Resource [Score]
 * Arg:              dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [EstFrag3 *]
 *
 */
EstFrag3 * allocate_Expl_EstFrag3(GeneWiseScore* query,ComplexSequence* target ,cDNAParserScore * cp,Score start_frag,Score end_frag,DPRunImpl * dpri) 
{
    EstFrag3 * out;  


    out = allocate_EstFrag3_only(query, target , cp, start_frag, end_frag);  
    if( out == NULL )    
      return NULL;   
    if( dpri->should_cache == TRUE ) {  
      if( dpri->cache != NULL )  {  
        if( dpri->cache->maxleni >= (out->lenj+5)*3 && dpri->cache->maxlenj >= (out->leni+1))    
          out->basematrix = hard_link_BaseMatrix(dpri->cache);   
        else 
          dpri->cache = free_BaseMatrix(dpri->cache);    
        }  
      }  
    if( out->basematrix == NULL )    {  
      if( (out->basematrix = BaseMatrix_alloc_matrix_and_specials((out->lenj+5)*3,(out->leni+1),2,out->lenj+5)) == NULL) {  
        warn("Explicit matrix EstFrag3 cannot be allocated, (asking for %d by %d main cells)",out->leni,out->lenj);  
        free_EstFrag3(out);  
        return NULL; 
        }  
      }  
    if( dpri->should_cache == TRUE && dpri->cache == NULL)   
      dpri->cache = hard_link_BaseMatrix(out->basematrix);   
    out->basematrix->type = BASEMATRIX_TYPE_EXPLICIT;    
    init_EstFrag3(out);  
    return out;  
}    


/* Function:  init_EstFrag3(mat)
 *
 * Descrip:    This function initates EstFrag3 matrix when in explicit mode
 *             Called in /allocate_Expl_EstFrag3
 *
 *
 * Arg:        mat [UNKN ] EstFrag3 which contains explicit basematrix memory [EstFrag3 *]
 *
 */
void init_EstFrag3(EstFrag3 * mat) 
{
    register int i;  
    register int j;  
    if( mat->basematrix->type != BASEMATRIX_TYPE_EXPLICIT)   {  
      warn("Cannot iniate matrix, is not an explicit memory type and you have assummed that");   
      return;    
      }  


    for(i= (-1);i<mat->query->len;i++)   {  
      for(j= (-5);j<6;j++)   {  
        EstFrag3_EXPL_MATRIX(mat,i,j,MATCH) = NEGI;  
        EstFrag3_EXPL_MATRIX(mat,i,j,INSERT) = NEGI; 
        EstFrag3_EXPL_MATRIX(mat,i,j,DELETE) = NEGI; 
        }  
      }  
    for(j= (-5);j<mat->target->seq->len;j++) {  
      for(i= (-1);i<2;i++)   {  
        EstFrag3_EXPL_MATRIX(mat,i,j,MATCH) = NEGI;  
        EstFrag3_EXPL_MATRIX(mat,i,j,INSERT) = NEGI; 
        EstFrag3_EXPL_MATRIX(mat,i,j,DELETE) = NEGI; 
        }  
      EstFrag3_EXPL_SPECIAL(mat,i,j,START) = 0;  
      EstFrag3_EXPL_SPECIAL(mat,i,j,END) = NEGI; 
      }  
    return;  
}    


/* Function:  recalculate_PackAln_EstFrag3(pal,mat)
 *
 * Descrip:    This function recalculates the PackAln structure produced by EstFrag3
 *             For example, in linear space methods this is used to score them
 *
 *
 * Arg:        pal [UNKN ] Undocumented argument [PackAln *]
 * Arg:        mat [UNKN ] Undocumented argument [EstFrag3 *]
 *
 */
void recalculate_PackAln_EstFrag3(PackAln * pal,EstFrag3 * mat) 
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
        case MATCH :     
          if( offi == 1 && offj == 3 && prev->state == MATCH )   {  
            pau->score = (mat->query->seg[i]->transition[GW_MATCH2MATCH]+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)]) + (0);   
            continue;    
            }  
          if( offi == 1 && offj == 3 && prev->state == INSERT )  {  
            pau->score = (mat->query->seg[i]->transition[GW_INSERT2MATCH]+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)]) + (0);  
            continue;    
            }  
          if( offi == 1 && offj == 3 && prev->state == DELETE )  {  
            pau->score = (mat->query->seg[i]->transition[GW_DELETE2MATCH]+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)]) + (0);  
            continue;    
            }  
          if( offj == 3 && prev->state == (START+3) )    {  
            pau->score = (mat->query->seg[i]->transition[GW_START2MATCH]+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)]) + (0);   
            continue;    
            }  
          if( offj == 3 && prev->state == (START+3) )    {  
            pau->score = (mat->start_frag+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)]) + (0);  
            continue;    
            }  
          if( offi == 1 && offj == 2 && prev->state == MATCH )   {  
            pau->score = mat->cp->trans[PCD_DELETE_1_BASE] + (0);    
            continue;    
            }  
          if( offi == 1 && offj == 1 && prev->state == MATCH )   {  
            pau->score = mat->cp->trans[PCD_DELETE_2_BASE] + (0);    
            continue;    
            }  
          if( offi == 1 && offj == 4 && prev->state == MATCH )   {  
            pau->score = mat->cp->trans[PCD_INSERT_1_BASE] + (0);    
            continue;    
            }  
          if( offi == 1 && offj == 5 && prev->state == MATCH )   {  
            pau->score = mat->cp->trans[PCD_INSERT_2_BASE] + (0);    
            continue;    
            }  
          if( offi == 1 && offj == 2 && prev->state == INSERT )  {  
            pau->score = mat->cp->trans[PCD_DELETE_1_BASE] + (0);    
            continue;    
            }  
          if( offi == 1 && offj == 1 && prev->state == INSERT )  {  
            pau->score = mat->cp->trans[PCD_DELETE_2_BASE] + (0);    
            continue;    
            }  
          if( offi == 1 && offj == 4 && prev->state == INSERT )  {  
            pau->score = mat->cp->trans[PCD_INSERT_1_BASE] + (0);    
            continue;    
            }  
          if( offi == 1 && offj == 5 && prev->state == INSERT )  {  
            pau->score = mat->cp->trans[PCD_INSERT_2_BASE] + (0);    
            continue;    
            }  
          if( offi == 1 && offj == 2 && prev->state == DELETE )  {  
            pau->score = mat->cp->trans[PCD_DELETE_1_BASE] + (0);    
            continue;    
            }  
          if( offi == 1 && offj == 1 && prev->state == DELETE )  {  
            pau->score = mat->cp->trans[PCD_DELETE_2_BASE] + (0);    
            continue;    
            }  
          if( offi == 1 && offj == 4 && prev->state == DELETE )  {  
            pau->score = mat->cp->trans[PCD_INSERT_1_BASE] + (0);    
            continue;    
            }  
          if( offi == 1 && offj == 5 && prev->state == DELETE )  {  
            pau->score = mat->cp->trans[PCD_INSERT_2_BASE] + (0);    
            continue;    
            }  
          warn("In recaluclating PackAln with state MATCH, from [%d,%d,%d], got a bad source state. Error!",offi,offj,prev->state);  
          break; 
        case INSERT :    
          if( offi == 0 && offj == 3 && prev->state == MATCH )   {  
            pau->score = (mat->query->seg[i]->transition[GW_MATCH2INSERT]+mat->query->seg[i]->insert[CSEQ_CDNA_CODON(mat->target,j)]) + (0);     
            continue;    
            }  
          if( offi == 0 && offj == 3 && prev->state == INSERT )  {  
            pau->score = (mat->query->seg[i]->transition[GW_INSERT2INSERT]+mat->query->seg[i]->insert[CSEQ_CDNA_CODON(mat->target,j)]) + (0);    
            continue;    
            }  
          if( offi == 0 && offj == 3 && prev->state == DELETE )  {  
            pau->score = (mat->query->seg[i]->transition[GW_DELETE2INSERT]+mat->query->seg[i]->insert[CSEQ_CDNA_CODON(mat->target,j)]) + (0);    
            continue;    
            }  
          if( offj == 3 && prev->state == (START+3) )    {  
            pau->score = (mat->query->seg[i]->transition[GW_START2INSERT]+mat->query->seg[i]->insert[CSEQ_CDNA_CODON(mat->target,j)]) + (0);     
            continue;    
            }  
          if( offi == 0 && offj == 2 && prev->state == INSERT )  {  
            pau->score = mat->cp->trans[PCD_DELETE_1_BASE] + (0);    
            continue;    
            }  
          if( offi == 0 && offj == 1 && prev->state == INSERT )  {  
            pau->score = mat->cp->trans[PCD_DELETE_2_BASE] + (0);    
            continue;    
            }  
          if( offi == 0 && offj == 4 && prev->state == INSERT )  {  
            pau->score = mat->cp->trans[PCD_INSERT_1_BASE] + (0);    
            continue;    
            }  
          if( offi == 0 && offj == 5 && prev->state == INSERT )  {  
            pau->score = mat->cp->trans[PCD_INSERT_2_BASE] + (0);    
            continue;    
            }  
          if( offi == 0 && offj == 2 && prev->state == MATCH )   {  
            pau->score = mat->cp->trans[PCD_DELETE_1_BASE] + (0);    
            continue;    
            }  
          if( offi == 0 && offj == 1 && prev->state == MATCH )   {  
            pau->score = mat->cp->trans[PCD_DELETE_2_BASE] + (0);    
            continue;    
            }  
          if( offi == 0 && offj == 4 && prev->state == MATCH )   {  
            pau->score = mat->cp->trans[PCD_INSERT_1_BASE] + (0);    
            continue;    
            }  
          if( offi == 0 && offj == 5 && prev->state == MATCH )   {  
            pau->score = mat->cp->trans[PCD_INSERT_2_BASE] + (0);    
            continue;    
            }  
          if( offi == 0 && offj == 2 && prev->state == DELETE )  {  
            pau->score = mat->cp->trans[PCD_DELETE_1_BASE] + (0);    
            continue;    
            }  
          if( offi == 0 && offj == 1 && prev->state == DELETE )  {  
            pau->score = mat->cp->trans[PCD_DELETE_2_BASE] + (0);    
            continue;    
            }  
          if( offi == 0 && offj == 4 && prev->state == DELETE )  {  
            pau->score = mat->cp->trans[PCD_INSERT_1_BASE] + (0);    
            continue;    
            }  
          if( offi == 0 && offj == 5 && prev->state == DELETE )  {  
            pau->score = mat->cp->trans[PCD_INSERT_2_BASE] + (0);    
            continue;    
            }  
          warn("In recaluclating PackAln with state INSERT, from [%d,%d,%d], got a bad source state. Error!",offi,offj,prev->state); 
          break; 
        case DELETE :    
          if( offi == 1 && offj == 0 && prev->state == MATCH )   {  
            pau->score = mat->query->seg[i]->transition[GW_MATCH2DELETE] + (0);  
            continue;    
            }  
          if( offi == 1 && offj == 0 && prev->state == INSERT )  {  
            pau->score = mat->query->seg[i]->transition[GW_INSERT2DELETE] + (0);     
            continue;    
            }  
          if( offi == 1 && offj == 0 && prev->state == DELETE )  {  
            pau->score = mat->query->seg[i]->transition[GW_DELETE2DELETE] + (0);     
            continue;    
            }  
          if( offj == 0 && prev->state == (START+3) )    {  
            pau->score = mat->query->seg[i]->transition[GW_START2DELETE] + (0);  
            continue;    
            }  
          warn("In recaluclating PackAln with state DELETE, from [%d,%d,%d], got a bad source state. Error!",offi,offj,prev->state); 
          break; 
        case (START+3) :     
          warn("In recaluclating PackAln with state START, got a bad source state. Error!"); 
          break; 
        case (END+3) :   
          if( offj == 0 && prev->state == MATCH )    {  
            /* i here comes from the previous state ;) - not the real one */ 
            i = prev->i; 
            pau->score = mat->query->seg[i]->transition[GW_MATCH2END] + (0);     
            continue;    
            }  
          if( offj == 0 && prev->state == MATCH )    {  
            /* i here comes from the previous state ;) - not the real one */ 
            i = prev->i; 
            pau->score = mat->end_frag + (0);    
            continue;    
            }  
          if( offj == 0 && prev->state == INSERT )   {  
            /* i here comes from the previous state ;) - not the real one */ 
            i = prev->i; 
            pau->score = mat->query->seg[i]->transition[GW_INSERT2END] + (0);    
            continue;    
            }  
          if( offj == 0 && prev->state == DELETE )   {  
            /* i here comes from the previous state ;) - not the real one */ 
            i = prev->i; 
            pau->score = mat->query->seg[i]->transition[GW_DELETE2END] + (0);    
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
#define EstFrag3_HIDDEN_MATRIX(thismatrix,i,j,state) (thismatrix->basematrix->matrix[(j-hiddenj+5)][(i+1)*3+state])  
#define EstFrag3_DC_SHADOW_MATRIX(thismatrix,i,j,state) (thismatrix->basematrix->matrix[((j+6)*8) % 48][(i+1)*3+state])  
#define EstFrag3_HIDDEN_SPECIAL(thismatrix,i,j,state) (thismatrix->basematrix->specmatrix[state][(j+5)]) 
#define EstFrag3_DC_SHADOW_SPECIAL(thismatrix,i,j,state) (thismatrix->basematrix->specmatrix[state*8][(j+5)])    
#define EstFrag3_DC_SHADOW_MATRIX_SP(thismatrix,i,j,state,shadow) (thismatrix->basematrix->matrix[((((j+6)*8)+(shadow+1)) % 48)][(i+1)*3 + state])   
#define EstFrag3_DC_SHADOW_SPECIAL_SP(thismatrix,i,j,state,shadow) (thismatrix->basematrix->specmatrix[state*8 +shadow+1][(j+5)])    
#define EstFrag3_DC_OPT_SHADOW_MATRIX(thismatrix,i,j,state) (score_pointers[(((j+5)% 5) * (leni+1) * 3) + ((i+1) * 3) + (state)])    
#define EstFrag3_DC_OPT_SHADOW_MATRIX_SP(thismatrix,i,j,state,shadow) (shadow_pointers[(((j+5)% 5) * (leni+1) * 24) + ((i+1) * 24) + (state * 8) + shadow+1])    
#define EstFrag3_DC_OPT_SHADOW_SPECIAL(thismatrix,i,j,state) (thismatrix->basematrix->specmatrix[state*8][(j+5)])    
/* Function:  allocate_Small_EstFrag3(query,target,cp,start_frag,end_frag)
 *
 * Descrip:    This function allocates the EstFrag3 structure
 *             and the basematrix area for a small memory implementations
 *             It calls /allocate_EstFrag3_only
 *
 *
 * Arg:             query [UNKN ] query data structure [GeneWiseScore*]
 * Arg:            target [UNKN ] target data structure [ComplexSequence*]
 * Arg:                cp [UNKN ] Resource [cDNAParserScore *]
 * Arg:        start_frag [UNKN ] Resource [Score]
 * Arg:          end_frag [UNKN ] Resource [Score]
 *
 * Return [UNKN ]  Undocumented return value [EstFrag3 *]
 *
 */
#define EstFrag3_DC_OPT_SHADOW_SPECIAL_SP(thismatrix,i,j,state,shadow) (thismatrix->basematrix->specmatrix[state*8 +shadow+1][(j+5)])    
EstFrag3 * allocate_Small_EstFrag3(GeneWiseScore* query,ComplexSequence* target ,cDNAParserScore * cp,Score start_frag,Score end_frag) 
{
    EstFrag3 * out;  


    out = allocate_EstFrag3_only(query, target , cp, start_frag, end_frag);  
    if( out == NULL )    
      return NULL;   
    out->basematrix = BaseMatrix_alloc_matrix_and_specials(48,(out->leni + 1) * 3,16,out->lenj+5);   
    if(out == NULL)  {  
      warn("Small shadow matrix EstFrag3 cannot be allocated, (asking for 6 by %d main cells)",out->leni+2); 
      free_EstFrag3(out);    
      return NULL;   
      }  
    out->basematrix->type = BASEMATRIX_TYPE_SHADOW;  
    return out;  
}    


/* Function:  PackAln_calculate_Small_EstFrag3(mat,dpenv)
 *
 * Descrip:    This function calculates an alignment for EstFrag3 structure in linear space
 *             If you want only the start/end points
 *             use /AlnRangeSet_calculate_Small_EstFrag3 
 *
 *             The function basically
 *               finds start/end points 
 *               foreach start/end point 
 *                 calls /full_dc_EstFrag3 
 *
 *
 * Arg:          mat [UNKN ] Undocumented argument [EstFrag3 *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * PackAln_calculate_Small_EstFrag3(EstFrag3 * mat,DPEnvelope * dpenv) 
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
      warn("Could not calculate packaln small for EstFrag3 due to wrong type of matrix");    
      return NULL;   
      }  


    out = PackAln_alloc_std();   


    start_reporting("Find start end points: ");  
    dc_optimised_start_end_calc_EstFrag3(mat,dpenv); 
    score = start_end_find_end_EstFrag3(mat,&endj);  
    out->score = score;  
    stopstate = END;
    
    /* No special to specials: one matrix alignment: simply remove and get */ 
    starti = EstFrag3_DC_SHADOW_SPECIAL_SP(mat,0,endj,END,0);    
    startj = EstFrag3_DC_SHADOW_SPECIAL_SP(mat,0,endj,END,1);    
    startstate = EstFrag3_DC_SHADOW_SPECIAL_SP(mat,0,endj,END,2);    
    stopi = EstFrag3_DC_SHADOW_SPECIAL_SP(mat,0,endj,END,3); 
    stopj = EstFrag3_DC_SHADOW_SPECIAL_SP(mat,0,endj,END,4); 
    stopstate = EstFrag3_DC_SHADOW_SPECIAL_SP(mat,0,endj,END,5); 
    temp = EstFrag3_DC_SHADOW_SPECIAL_SP(mat,0,endj,END,6);  
    log_full_error(REPORT,0,"[%d,%d][%d,%d] Score %d",starti,startj,stopi,stopj,score);  
    stop_reporting();    
    start_reporting("Recovering alignment: ");   


    /* Figuring how much j we have to align for reporting purposes */ 
    donej = 0;   
    totalj = stopj - startj; 
    full_dc_EstFrag3(mat,starti,startj,startstate,stopi,stopj,stopstate,out,&donej,totalj,dpenv);    


    /* Although we have no specials, need to get start. Better to check than assume */ 


    max_matrix_to_special_EstFrag3(mat,starti,startj,startstate,temp,&stopi,&stopj,&stopstate,&temp,NULL);   
    if( stopi == EstFrag3_READ_OFF_ERROR || stopstate != START ) {  
      warn("Problem in reading off special state system, hit a non start state (or an internal error) in a single alignment mode");  
      invert_PackAln(out);   
      recalculate_PackAln_EstFrag3(out,mat); 
      return out;    
      }  


    /* Ok. Put away start start... */ 
    pau = PackAlnUnit_alloc();   
    pau->i = stopi;  
    pau->j = stopj;  
    pau->state = stopstate + 3;  
    add_PackAln(out,pau);    


    log_full_error(REPORT,0,"Alignment recovered");  
    stop_reporting();    
    invert_PackAln(out); 
    recalculate_PackAln_EstFrag3(out,mat);   
    return out;  


}    


/* Function:  AlnRangeSet_calculate_Small_EstFrag3(mat)
 *
 * Descrip:    This function calculates an alignment for EstFrag3 structure in linear space
 *             If you want the full alignment, use /PackAln_calculate_Small_EstFrag3 
 *             If you have already got the full alignment, but want the range set, use /AlnRangeSet_from_PackAln_EstFrag3
 *             If you have got the small matrix but not the alignment, use /AlnRangeSet_from_EstFrag3 
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [EstFrag3 *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * AlnRangeSet_calculate_Small_EstFrag3(EstFrag3 * mat) 
{
    AlnRangeSet * out;   


    start_reporting("Find start end points: ");  
    dc_optimised_start_end_calc_EstFrag3(mat,NULL);  
    log_full_error(REPORT,0,"Calculated");   


    out = AlnRangeSet_from_EstFrag3(mat);    
    return out;  
}    


/* Function:  AlnRangeSet_from_EstFrag3(mat)
 *
 * Descrip:    This function reads off a start/end structure
 *             for EstFrag3 structure in linear space
 *             If you want the full alignment use
 *             /PackAln_calculate_Small_EstFrag3 
 *             If you have not calculated the matrix use
 *             /AlnRange_calculate_Small_EstFrag3
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [EstFrag3 *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * AlnRangeSet_from_EstFrag3(EstFrag3 * mat) 
{
    AlnRangeSet * out;   
    AlnRange * temp; 
    int jpos;    
    int state;   


    if( mat->basematrix->type != BASEMATRIX_TYPE_SHADOW) {  
      warn("Bad error! - non shadow matrix type in AlnRangeSet_from_EstFrag3");  
      return NULL;   
      }  


    out = AlnRangeSet_alloc_std();   
    /* Find the end position */ 
    out->score = start_end_find_end_EstFrag3(mat,&jpos); 
    state = END; 


    while( (temp = AlnRange_build_EstFrag3(mat,jpos,state,&jpos,&state)) != NULL)    
      add_AlnRangeSet(out,temp); 
    return out;  
}    


/* Function:  AlnRange_build_EstFrag3(mat,stopj,stopspecstate,startj,startspecstate)
 *
 * Descrip:    This function calculates a single start/end set in linear space
 *             Really a sub-routine for /AlnRangeSet_from_PackAln_EstFrag3
 *
 *
 * Arg:                   mat [UNKN ] Undocumented argument [EstFrag3 *]
 * Arg:                 stopj [UNKN ] Undocumented argument [int]
 * Arg:         stopspecstate [UNKN ] Undocumented argument [int]
 * Arg:                startj [UNKN ] Undocumented argument [int *]
 * Arg:        startspecstate [UNKN ] Undocumented argument [int *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRange *]
 *
 */
AlnRange * AlnRange_build_EstFrag3(EstFrag3 * mat,int stopj,int stopspecstate,int * startj,int * startspecstate) 
{
    AlnRange * out;  
    int jpos;    
    int state;   


    if( mat->basematrix->type != BASEMATRIX_TYPE_SHADOW) {  
      warn("Bad error! - non shadow matrix type in AlnRangeSet_from_EstFrag3");  
      return NULL;   
      }  


    /* Assumme that we have specials (we should!). Read back along the specials till we have the finish point */ 
    if( read_special_strip_EstFrag3(mat,0,stopj,stopspecstate,&jpos,&state,NULL) == FALSE)   {  
      warn("In AlnRanger_build_EstFrag3 alignment ending at %d, unable to read back specials. Will (evenutally) return a partial range set... BEWARE!",stopj);   
      return NULL;   
      }  
    if( state == START || jpos <= 0) 
      return NULL;   


    out = AlnRange_alloc();  


    out->starti = EstFrag3_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,0); 
    out->startj = EstFrag3_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,1); 
    out->startstate = EstFrag3_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,2); 
    out->stopi = EstFrag3_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,3);  
    out->stopj = EstFrag3_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,4);  
    out->stopstate = EstFrag3_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,5);  
    out->startscore = EstFrag3_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,6); 
    out->stopscore = EstFrag3_DC_SHADOW_SPECIAL(mat,0,jpos,state);   


    /* Now, we have to figure out where this state came from in the specials */ 
    max_matrix_to_special_EstFrag3(mat,out->starti,out->startj,out->startstate,out->startscore,&jpos,startj,startspecstate,&state,NULL); 
    if( jpos == EstFrag3_READ_OFF_ERROR) {  
      warn("In AlnRange_build_EstFrag3 alignment ending at %d, with aln range between %d-%d in j, unable to find source special, returning this range, but this could get tricky!",stopj,out->startj,out->stopj);    
      return out;    
      }  


    /* Put in the correct score for startstate, from the special */ 
    out->startscore = EstFrag3_DC_SHADOW_SPECIAL(mat,0,*startj,*startspecstate); 
    /* The correct j coords have been put into startj, startspecstate... so just return out */ 
    return out;  
}    


/* Function:  read_hidden_EstFrag3(mat,starti,startj,startstate,stopi,stopj,stopstate,out)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [EstFrag3 *]
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
boolean read_hidden_EstFrag3(EstFrag3 * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out) 
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


      max_hidden_EstFrag3(mat,startj,i,j,state,isspecial,&i,&j,&state,&isspecial,&cellscore);    


      if( i == EstFrag3_READ_OFF_ERROR)  {  
        warn("In EstFrag3 hidden read off, between %d:%d,%d:%d - at got bad read off. Problem!",starti,startj,stopi,stopj);  
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
        warn("In EstFrag3 hidden read off, between %d:%d,%d:%d - hit start cell, but not in start state. Can't be good!.",starti,startj,stopi,stopj);    
        return FALSE;    
        }  
      }  
    warn("In EstFrag3 hidden read off, between %d:%d,%d:%d - gone past start cell (now in %d,%d,%d), can't be good news!.",starti,startj,stopi,stopj,i,j,state); 
    return FALSE;    
}    


/* Function:  max_hidden_EstFrag3(mat,hiddenj,i,j,state,isspecial,reti,retj,retstate,retspecial,cellscore)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [EstFrag3 *]
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
int max_hidden_EstFrag3(EstFrag3 * mat,int hiddenj,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore) 
{
    register int temp;   
    register int cscore; 


    *reti = (*retj) = (*retstate) = EstFrag3_READ_OFF_ERROR; 


    if( i < 0 || j < 0 || i > mat->query->len || j > mat->target->seq->len)  {  
      warn("In EstFrag3 matrix special read off - out of bounds on matrix [i,j is %d,%d state %d in standard matrix]",i,j,state);    
      return -1; 
      }  


    /* Then you have to select the correct switch statement to figure out the readoff      */ 
    /* Somewhat odd - reverse the order of calculation and return as soon as it is correct */ 
    cscore = EstFrag3_HIDDEN_MATRIX(mat,i,j,state);  
    switch(state)    { /*Switch state */ 
      case MATCH :   
        temp = cscore - (mat->cp->trans[PCD_INSERT_2_BASE]) -  (0);  
        if( temp == EstFrag3_HIDDEN_MATRIX(mat,i - 1,j - 5,DELETE) ) {  
          *reti = i - 1; 
          *retj = j - 5; 
          *retstate = DELETE;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - EstFrag3_HIDDEN_MATRIX(mat,i-1,j-5,DELETE);    
            }  
          return EstFrag3_HIDDEN_MATRIX(mat,i - 1,j - 5,DELETE);     
          }  
        temp = cscore - (mat->cp->trans[PCD_INSERT_1_BASE]) -  (0);  
        if( temp == EstFrag3_HIDDEN_MATRIX(mat,i - 1,j - 4,DELETE) ) {  
          *reti = i - 1; 
          *retj = j - 4; 
          *retstate = DELETE;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - EstFrag3_HIDDEN_MATRIX(mat,i-1,j-4,DELETE);    
            }  
          return EstFrag3_HIDDEN_MATRIX(mat,i - 1,j - 4,DELETE);     
          }  
        temp = cscore - (mat->cp->trans[PCD_DELETE_2_BASE]) -  (0);  
        if( temp == EstFrag3_HIDDEN_MATRIX(mat,i - 1,j - 1,DELETE) ) {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = DELETE;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - EstFrag3_HIDDEN_MATRIX(mat,i-1,j-1,DELETE);    
            }  
          return EstFrag3_HIDDEN_MATRIX(mat,i - 1,j - 1,DELETE);     
          }  
        temp = cscore - (mat->cp->trans[PCD_DELETE_1_BASE]) -  (0);  
        if( temp == EstFrag3_HIDDEN_MATRIX(mat,i - 1,j - 2,DELETE) ) {  
          *reti = i - 1; 
          *retj = j - 2; 
          *retstate = DELETE;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - EstFrag3_HIDDEN_MATRIX(mat,i-1,j-2,DELETE);    
            }  
          return EstFrag3_HIDDEN_MATRIX(mat,i - 1,j - 2,DELETE);     
          }  
        temp = cscore - (mat->cp->trans[PCD_INSERT_2_BASE]) -  (0);  
        if( temp == EstFrag3_HIDDEN_MATRIX(mat,i - 1,j - 5,INSERT) ) {  
          *reti = i - 1; 
          *retj = j - 5; 
          *retstate = INSERT;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - EstFrag3_HIDDEN_MATRIX(mat,i-1,j-5,INSERT);    
            }  
          return EstFrag3_HIDDEN_MATRIX(mat,i - 1,j - 5,INSERT);     
          }  
        temp = cscore - (mat->cp->trans[PCD_INSERT_1_BASE]) -  (0);  
        if( temp == EstFrag3_HIDDEN_MATRIX(mat,i - 1,j - 4,INSERT) ) {  
          *reti = i - 1; 
          *retj = j - 4; 
          *retstate = INSERT;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - EstFrag3_HIDDEN_MATRIX(mat,i-1,j-4,INSERT);    
            }  
          return EstFrag3_HIDDEN_MATRIX(mat,i - 1,j - 4,INSERT);     
          }  
        temp = cscore - (mat->cp->trans[PCD_DELETE_2_BASE]) -  (0);  
        if( temp == EstFrag3_HIDDEN_MATRIX(mat,i - 1,j - 1,INSERT) ) {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = INSERT;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - EstFrag3_HIDDEN_MATRIX(mat,i-1,j-1,INSERT);    
            }  
          return EstFrag3_HIDDEN_MATRIX(mat,i - 1,j - 1,INSERT);     
          }  
        temp = cscore - (mat->cp->trans[PCD_DELETE_1_BASE]) -  (0);  
        if( temp == EstFrag3_HIDDEN_MATRIX(mat,i - 1,j - 2,INSERT) ) {  
          *reti = i - 1; 
          *retj = j - 2; 
          *retstate = INSERT;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - EstFrag3_HIDDEN_MATRIX(mat,i-1,j-2,INSERT);    
            }  
          return EstFrag3_HIDDEN_MATRIX(mat,i - 1,j - 2,INSERT);     
          }  
        temp = cscore - (mat->cp->trans[PCD_INSERT_2_BASE]) -  (0);  
        if( temp == EstFrag3_HIDDEN_MATRIX(mat,i - 1,j - 5,MATCH) )  {  
          *reti = i - 1; 
          *retj = j - 5; 
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - EstFrag3_HIDDEN_MATRIX(mat,i-1,j-5,MATCH); 
            }  
          return EstFrag3_HIDDEN_MATRIX(mat,i - 1,j - 5,MATCH);  
          }  
        temp = cscore - (mat->cp->trans[PCD_INSERT_1_BASE]) -  (0);  
        if( temp == EstFrag3_HIDDEN_MATRIX(mat,i - 1,j - 4,MATCH) )  {  
          *reti = i - 1; 
          *retj = j - 4; 
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - EstFrag3_HIDDEN_MATRIX(mat,i-1,j-4,MATCH); 
            }  
          return EstFrag3_HIDDEN_MATRIX(mat,i - 1,j - 4,MATCH);  
          }  
        temp = cscore - (mat->cp->trans[PCD_DELETE_2_BASE]) -  (0);  
        if( temp == EstFrag3_HIDDEN_MATRIX(mat,i - 1,j - 1,MATCH) )  {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - EstFrag3_HIDDEN_MATRIX(mat,i-1,j-1,MATCH); 
            }  
          return EstFrag3_HIDDEN_MATRIX(mat,i - 1,j - 1,MATCH);  
          }  
        temp = cscore - (mat->cp->trans[PCD_DELETE_1_BASE]) -  (0);  
        if( temp == EstFrag3_HIDDEN_MATRIX(mat,i - 1,j - 2,MATCH) )  {  
          *reti = i - 1; 
          *retj = j - 2; 
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - EstFrag3_HIDDEN_MATRIX(mat,i-1,j-2,MATCH); 
            }  
          return EstFrag3_HIDDEN_MATRIX(mat,i - 1,j - 2,MATCH);  
          }  
        /* Not allowing special sources.. skipping START */ 
        /* Not allowing special sources.. skipping START */ 
        temp = cscore - ((mat->query->seg[i]->transition[GW_DELETE2MATCH]+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)])) -  (0);    
        if( temp == EstFrag3_HIDDEN_MATRIX(mat,i - 1,j - 3,DELETE) ) {  
          *reti = i - 1; 
          *retj = j - 3; 
          *retstate = DELETE;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - EstFrag3_HIDDEN_MATRIX(mat,i-1,j-3,DELETE);    
            }  
          return EstFrag3_HIDDEN_MATRIX(mat,i - 1,j - 3,DELETE);     
          }  
        temp = cscore - ((mat->query->seg[i]->transition[GW_INSERT2MATCH]+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)])) -  (0);    
        if( temp == EstFrag3_HIDDEN_MATRIX(mat,i - 1,j - 3,INSERT) ) {  
          *reti = i - 1; 
          *retj = j - 3; 
          *retstate = INSERT;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - EstFrag3_HIDDEN_MATRIX(mat,i-1,j-3,INSERT);    
            }  
          return EstFrag3_HIDDEN_MATRIX(mat,i - 1,j - 3,INSERT);     
          }  
        temp = cscore - ((mat->query->seg[i]->transition[GW_MATCH2MATCH]+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)])) -  (0); 
        if( temp == EstFrag3_HIDDEN_MATRIX(mat,i - 1,j - 3,MATCH) )  {  
          *reti = i - 1; 
          *retj = j - 3; 
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - EstFrag3_HIDDEN_MATRIX(mat,i-1,j-3,MATCH); 
            }  
          return EstFrag3_HIDDEN_MATRIX(mat,i - 1,j - 3,MATCH);  
          }  
        warn("Major problem (!) - in EstFrag3 read off, position %d,%d state %d no source found!",i,j,state);    
        return (-1); 
      case INSERT :  
        temp = cscore - (mat->cp->trans[PCD_INSERT_2_BASE]) -  (0);  
        if( temp == EstFrag3_HIDDEN_MATRIX(mat,i - 0,j - 5,DELETE) ) {  
          *reti = i - 0; 
          *retj = j - 5; 
          *retstate = DELETE;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - EstFrag3_HIDDEN_MATRIX(mat,i-0,j-5,DELETE);    
            }  
          return EstFrag3_HIDDEN_MATRIX(mat,i - 0,j - 5,DELETE);     
          }  
        temp = cscore - (mat->cp->trans[PCD_INSERT_1_BASE]) -  (0);  
        if( temp == EstFrag3_HIDDEN_MATRIX(mat,i - 0,j - 4,DELETE) ) {  
          *reti = i - 0; 
          *retj = j - 4; 
          *retstate = DELETE;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - EstFrag3_HIDDEN_MATRIX(mat,i-0,j-4,DELETE);    
            }  
          return EstFrag3_HIDDEN_MATRIX(mat,i - 0,j - 4,DELETE);     
          }  
        temp = cscore - (mat->cp->trans[PCD_DELETE_2_BASE]) -  (0);  
        if( temp == EstFrag3_HIDDEN_MATRIX(mat,i - 0,j - 1,DELETE) ) {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = DELETE;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - EstFrag3_HIDDEN_MATRIX(mat,i-0,j-1,DELETE);    
            }  
          return EstFrag3_HIDDEN_MATRIX(mat,i - 0,j - 1,DELETE);     
          }  
        temp = cscore - (mat->cp->trans[PCD_DELETE_1_BASE]) -  (0);  
        if( temp == EstFrag3_HIDDEN_MATRIX(mat,i - 0,j - 2,DELETE) ) {  
          *reti = i - 0; 
          *retj = j - 2; 
          *retstate = DELETE;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - EstFrag3_HIDDEN_MATRIX(mat,i-0,j-2,DELETE);    
            }  
          return EstFrag3_HIDDEN_MATRIX(mat,i - 0,j - 2,DELETE);     
          }  
        temp = cscore - (mat->cp->trans[PCD_INSERT_2_BASE]) -  (0);  
        if( temp == EstFrag3_HIDDEN_MATRIX(mat,i - 0,j - 5,MATCH) )  {  
          *reti = i - 0; 
          *retj = j - 5; 
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - EstFrag3_HIDDEN_MATRIX(mat,i-0,j-5,MATCH); 
            }  
          return EstFrag3_HIDDEN_MATRIX(mat,i - 0,j - 5,MATCH);  
          }  
        temp = cscore - (mat->cp->trans[PCD_INSERT_1_BASE]) -  (0);  
        if( temp == EstFrag3_HIDDEN_MATRIX(mat,i - 0,j - 4,MATCH) )  {  
          *reti = i - 0; 
          *retj = j - 4; 
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - EstFrag3_HIDDEN_MATRIX(mat,i-0,j-4,MATCH); 
            }  
          return EstFrag3_HIDDEN_MATRIX(mat,i - 0,j - 4,MATCH);  
          }  
        temp = cscore - (mat->cp->trans[PCD_DELETE_2_BASE]) -  (0);  
        if( temp == EstFrag3_HIDDEN_MATRIX(mat,i - 0,j - 1,MATCH) )  {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - EstFrag3_HIDDEN_MATRIX(mat,i-0,j-1,MATCH); 
            }  
          return EstFrag3_HIDDEN_MATRIX(mat,i - 0,j - 1,MATCH);  
          }  
        temp = cscore - (mat->cp->trans[PCD_DELETE_1_BASE]) -  (0);  
        if( temp == EstFrag3_HIDDEN_MATRIX(mat,i - 0,j - 2,MATCH) )  {  
          *reti = i - 0; 
          *retj = j - 2; 
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - EstFrag3_HIDDEN_MATRIX(mat,i-0,j-2,MATCH); 
            }  
          return EstFrag3_HIDDEN_MATRIX(mat,i - 0,j - 2,MATCH);  
          }  
        temp = cscore - (mat->cp->trans[PCD_INSERT_2_BASE]) -  (0);  
        if( temp == EstFrag3_HIDDEN_MATRIX(mat,i - 0,j - 5,INSERT) ) {  
          *reti = i - 0; 
          *retj = j - 5; 
          *retstate = INSERT;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - EstFrag3_HIDDEN_MATRIX(mat,i-0,j-5,INSERT);    
            }  
          return EstFrag3_HIDDEN_MATRIX(mat,i - 0,j - 5,INSERT);     
          }  
        temp = cscore - (mat->cp->trans[PCD_INSERT_1_BASE]) -  (0);  
        if( temp == EstFrag3_HIDDEN_MATRIX(mat,i - 0,j - 4,INSERT) ) {  
          *reti = i - 0; 
          *retj = j - 4; 
          *retstate = INSERT;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - EstFrag3_HIDDEN_MATRIX(mat,i-0,j-4,INSERT);    
            }  
          return EstFrag3_HIDDEN_MATRIX(mat,i - 0,j - 4,INSERT);     
          }  
        temp = cscore - (mat->cp->trans[PCD_DELETE_2_BASE]) -  (0);  
        if( temp == EstFrag3_HIDDEN_MATRIX(mat,i - 0,j - 1,INSERT) ) {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = INSERT;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - EstFrag3_HIDDEN_MATRIX(mat,i-0,j-1,INSERT);    
            }  
          return EstFrag3_HIDDEN_MATRIX(mat,i - 0,j - 1,INSERT);     
          }  
        temp = cscore - (mat->cp->trans[PCD_DELETE_1_BASE]) -  (0);  
        if( temp == EstFrag3_HIDDEN_MATRIX(mat,i - 0,j - 2,INSERT) ) {  
          *reti = i - 0; 
          *retj = j - 2; 
          *retstate = INSERT;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - EstFrag3_HIDDEN_MATRIX(mat,i-0,j-2,INSERT);    
            }  
          return EstFrag3_HIDDEN_MATRIX(mat,i - 0,j - 2,INSERT);     
          }  
        /* Not allowing special sources.. skipping START */ 
        temp = cscore - ((mat->query->seg[i]->transition[GW_DELETE2INSERT]+mat->query->seg[i]->insert[CSEQ_CDNA_CODON(mat->target,j)])) -  (0);  
        if( temp == EstFrag3_HIDDEN_MATRIX(mat,i - 0,j - 3,DELETE) ) {  
          *reti = i - 0; 
          *retj = j - 3; 
          *retstate = DELETE;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - EstFrag3_HIDDEN_MATRIX(mat,i-0,j-3,DELETE);    
            }  
          return EstFrag3_HIDDEN_MATRIX(mat,i - 0,j - 3,DELETE);     
          }  
        temp = cscore - ((mat->query->seg[i]->transition[GW_INSERT2INSERT]+mat->query->seg[i]->insert[CSEQ_CDNA_CODON(mat->target,j)])) -  (0);  
        if( temp == EstFrag3_HIDDEN_MATRIX(mat,i - 0,j - 3,INSERT) ) {  
          *reti = i - 0; 
          *retj = j - 3; 
          *retstate = INSERT;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - EstFrag3_HIDDEN_MATRIX(mat,i-0,j-3,INSERT);    
            }  
          return EstFrag3_HIDDEN_MATRIX(mat,i - 0,j - 3,INSERT);     
          }  
        temp = cscore - ((mat->query->seg[i]->transition[GW_MATCH2INSERT]+mat->query->seg[i]->insert[CSEQ_CDNA_CODON(mat->target,j)])) -  (0);   
        if( temp == EstFrag3_HIDDEN_MATRIX(mat,i - 0,j - 3,MATCH) )  {  
          *reti = i - 0; 
          *retj = j - 3; 
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - EstFrag3_HIDDEN_MATRIX(mat,i-0,j-3,MATCH); 
            }  
          return EstFrag3_HIDDEN_MATRIX(mat,i - 0,j - 3,MATCH);  
          }  
        warn("Major problem (!) - in EstFrag3 read off, position %d,%d state %d no source found!",i,j,state);    
        return (-1); 
      case DELETE :  
        /* Not allowing special sources.. skipping START */ 
        temp = cscore - (mat->query->seg[i]->transition[GW_DELETE2DELETE]) -  (0);   
        if( temp == EstFrag3_HIDDEN_MATRIX(mat,i - 1,j - 0,DELETE) ) {  
          *reti = i - 1; 
          *retj = j - 0; 
          *retstate = DELETE;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - EstFrag3_HIDDEN_MATRIX(mat,i-1,j-0,DELETE);    
            }  
          return EstFrag3_HIDDEN_MATRIX(mat,i - 1,j - 0,DELETE);     
          }  
        temp = cscore - (mat->query->seg[i]->transition[GW_INSERT2DELETE]) -  (0);   
        if( temp == EstFrag3_HIDDEN_MATRIX(mat,i - 1,j - 0,INSERT) ) {  
          *reti = i - 1; 
          *retj = j - 0; 
          *retstate = INSERT;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - EstFrag3_HIDDEN_MATRIX(mat,i-1,j-0,INSERT);    
            }  
          return EstFrag3_HIDDEN_MATRIX(mat,i - 1,j - 0,INSERT);     
          }  
        temp = cscore - (mat->query->seg[i]->transition[GW_MATCH2DELETE]) -  (0);    
        if( temp == EstFrag3_HIDDEN_MATRIX(mat,i - 1,j - 0,MATCH) )  {  
          *reti = i - 1; 
          *retj = j - 0; 
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - EstFrag3_HIDDEN_MATRIX(mat,i-1,j-0,MATCH); 
            }  
          return EstFrag3_HIDDEN_MATRIX(mat,i - 1,j - 0,MATCH);  
          }  
        warn("Major problem (!) - in EstFrag3 read off, position %d,%d state %d no source found!",i,j,state);    
        return (-1); 
      default:   
        warn("Major problem (!) - in EstFrag3 read off, position %d,%d state %d no source found!",i,j,state);    
        return (-1); 
      } /* end of Switch state  */ 
}    


/* Function:  read_special_strip_EstFrag3(mat,stopi,stopj,stopstate,startj,startstate,out)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [EstFrag3 *]
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
boolean read_special_strip_EstFrag3(EstFrag3 * mat,int stopi,int stopj,int stopstate,int * startj,int * startstate,PackAln * out) 
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
    while( j > EstFrag3_DC_SHADOW_SPECIAL_SP(mat,i,j,state,4) && state != START) { /*while more specials to eat up*/ 
      /* Put away current state, if we should */ 
      if(out != NULL)    {  
        pau = PackAlnUnit_alloc();  /* Should deal with memory overflow */ 
        pau->i = i;  
        pau->j = j;  
        pau->state =  state + 3; 
        add_PackAln(out,pau);    
        }  


      max_special_strip_EstFrag3(mat,i,j,state,isspecial,&i,&j,&state,&isspecial,&cellscore);    
      if( i == EstFrag3_READ_OFF_ERROR)  {  
        warn("In special strip read EstFrag3, got a bad read off error. Sorry!");    
        return FALSE;    
        }  
      } /* end of while more specials to eat up */ 


    /* check to see we have not gone too far! */ 
    if( state != START && j < EstFrag3_DC_SHADOW_SPECIAL_SP(mat,i,j,state,4))    {  
      warn("In special strip read EstFrag3, at special [%d] state [%d] overshot!",j,state);  
      return FALSE;  
      }  
    /* Put away last state */ 
    if(out != NULL)  {  
      pau = PackAlnUnit_alloc();/* Should deal with memory overflow */ 
      pau->i = i;    
      pau->j = j;    
      pau->state =  state + 3;   
      add_PackAln(out,pau);  
      }  


    /* Put away where we are in startj and startstate */ 
    *startj = j; 
    *startstate = state; 
    return TRUE; 
}    


/* Function:  max_special_strip_EstFrag3(mat,i,j,state,isspecial,reti,retj,retstate,retspecial,cellscore)
 *
 * Descrip:    A pretty intense internal function. Deals with read-off only in specials
 *
 *
 * Arg:               mat [UNKN ] Undocumented argument [EstFrag3 *]
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
int max_special_strip_EstFrag3(EstFrag3 * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore) 
{
    int temp;    
    int cscore;  


    *reti = (*retj) = (*retstate) = EstFrag3_READ_OFF_ERROR; 
    if( isspecial == FALSE ) {  
      warn("In special strip max function for EstFrag3, got a non special start point. Problem! (bad!)");    
      return (-1);   
      }  


    if( j < 0 || j > mat->target->seq->len)  {  
      warn("In EstFrag3 matrix special read off - out of bounds on matrix [j is %d in special]",j);  
      return -1; 
      }  


    cscore = EstFrag3_DC_SHADOW_SPECIAL(mat,i,j,state);  
    switch(state)    { /*switch on special states*/ 
      case START :   
      case END :     
        /* Source DELETE is not a special */ 
        /* Source INSERT is not a special */ 
        /* Source MATCH is not a special */ 
        /* Source MATCH is not a special */ 
      default:   
        warn("Major problem (!) - in EstFrag3 special strip read off, position %d,%d state %d no source found  dropped into default on source switch!",i,j,state);   
        return (-1); 
      } /* end of switch on special states */ 
}    


/* Function:  max_matrix_to_special_EstFrag3(mat,i,j,state,cscore,reti,retj,retstate,retspecial,cellscore)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [EstFrag3 *]
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
int max_matrix_to_special_EstFrag3(EstFrag3 * mat,int i,int j,int state,int cscore,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore) 
{
    int temp;    
    *reti = (*retj) = (*retstate) = EstFrag3_READ_OFF_ERROR; 


    if( j < 0 || j > mat->lenj)  {  
      warn("In EstFrag3 matrix to special read off - out of bounds on matrix [j is %d in special]",j);   
      return -1; 
      }  


    switch(state)    { /*Switch state */ 
      case MATCH :   
        /* Source DELETE is not a special, should not get here! */ 
        /* Source DELETE is not a special, should not get here! */ 
        /* Source DELETE is not a special, should not get here! */ 
        /* Source DELETE is not a special, should not get here! */ 
        /* Source INSERT is not a special, should not get here! */ 
        /* Source INSERT is not a special, should not get here! */ 
        /* Source INSERT is not a special, should not get here! */ 
        /* Source INSERT is not a special, should not get here! */ 
        /* Source MATCH is not a special, should not get here! */ 
        /* Source MATCH is not a special, should not get here! */ 
        /* Source MATCH is not a special, should not get here! */ 
        /* Source MATCH is not a special, should not get here! */ 
        temp = cscore - ((mat->start_frag+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)])) -  (0);    
        if( temp == EstFrag3_DC_SHADOW_SPECIAL(mat,i - 1,j - 3,START) )  {  
          *reti = i - 1; 
          *retj = j - 3; 
          *retstate = START; 
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - EstFrag3_DC_SHADOW_SPECIAL(mat,i-1,j-3,START);     
            }  
          return EstFrag3_DC_SHADOW_MATRIX(mat,i - 1,j - 3,START) ;  
          }  
        temp = cscore - ((mat->query->seg[i]->transition[GW_START2MATCH]+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)])) -  (0);     
        if( temp == EstFrag3_DC_SHADOW_SPECIAL(mat,i - 1,j - 3,START) )  {  
          *reti = i - 1; 
          *retj = j - 3; 
          *retstate = START; 
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - EstFrag3_DC_SHADOW_SPECIAL(mat,i-1,j-3,START);     
            }  
          return EstFrag3_DC_SHADOW_MATRIX(mat,i - 1,j - 3,START) ;  
          }  
        /* Source DELETE is not a special, should not get here! */ 
        /* Source INSERT is not a special, should not get here! */ 
        /* Source MATCH is not a special, should not get here! */ 
        warn("Major problem (!) - in EstFrag3 matrix to special read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      case INSERT :  
        /* Source DELETE is not a special, should not get here! */ 
        /* Source DELETE is not a special, should not get here! */ 
        /* Source DELETE is not a special, should not get here! */ 
        /* Source DELETE is not a special, should not get here! */ 
        /* Source MATCH is not a special, should not get here! */ 
        /* Source MATCH is not a special, should not get here! */ 
        /* Source MATCH is not a special, should not get here! */ 
        /* Source MATCH is not a special, should not get here! */ 
        /* Source INSERT is not a special, should not get here! */ 
        /* Source INSERT is not a special, should not get here! */ 
        /* Source INSERT is not a special, should not get here! */ 
        /* Source INSERT is not a special, should not get here! */ 
        temp = cscore - ((mat->query->seg[i]->transition[GW_START2INSERT]+mat->query->seg[i]->insert[CSEQ_CDNA_CODON(mat->target,j)])) -  (0);   
        if( temp == EstFrag3_DC_SHADOW_SPECIAL(mat,i - 0,j - 3,START) )  {  
          *reti = i - 0; 
          *retj = j - 3; 
          *retstate = START; 
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - EstFrag3_DC_SHADOW_SPECIAL(mat,i-0,j-3,START);     
            }  
          return EstFrag3_DC_SHADOW_MATRIX(mat,i - 0,j - 3,START) ;  
          }  
        /* Source DELETE is not a special, should not get here! */ 
        /* Source INSERT is not a special, should not get here! */ 
        /* Source MATCH is not a special, should not get here! */ 
        warn("Major problem (!) - in EstFrag3 matrix to special read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      case DELETE :  
        temp = cscore - (mat->query->seg[i]->transition[GW_START2DELETE]) -  (0);    
        if( temp == EstFrag3_DC_SHADOW_SPECIAL(mat,i - 1,j - 0,START) )  {  
          *reti = i - 1; 
          *retj = j - 0; 
          *retstate = START; 
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - EstFrag3_DC_SHADOW_SPECIAL(mat,i-1,j-0,START);     
            }  
          return EstFrag3_DC_SHADOW_MATRIX(mat,i - 1,j - 0,START) ;  
          }  
        /* Source DELETE is not a special, should not get here! */ 
        /* Source INSERT is not a special, should not get here! */ 
        /* Source MATCH is not a special, should not get here! */ 
        warn("Major problem (!) - in EstFrag3 matrix to special read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      default:   
        warn("Major problem (!) - in EstFrag3 read off, position %d,%d state %d no source found!",i,j,state);    
        return (-1); 
      } /* end of Switch state  */ 


}    


/* Function:  calculate_hidden_EstFrag3(mat,starti,startj,startstate,stopi,stopj,stopstate,dpenv)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [EstFrag3 *]
 * Arg:            starti [UNKN ] Undocumented argument [int]
 * Arg:            startj [UNKN ] Undocumented argument [int]
 * Arg:        startstate [UNKN ] Undocumented argument [int]
 * Arg:             stopi [UNKN ] Undocumented argument [int]
 * Arg:             stopj [UNKN ] Undocumented argument [int]
 * Arg:         stopstate [UNKN ] Undocumented argument [int]
 * Arg:             dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 */
void calculate_hidden_EstFrag3(EstFrag3 * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv) 
{
    register int i;  
    register int j;  
    register int score;  
    register int temp;   
    register int hiddenj;    


    hiddenj = startj;    


    init_hidden_EstFrag3(mat,starti,startj,stopi,stopj);     


    EstFrag3_HIDDEN_MATRIX(mat,starti,startj,startstate) = 0;    


    for(j=startj;j<=stopj;j++)   {  
      for(i=starti;i<=stopi;i++) {  
        /* Should *not* do very first cell as this is the one set to zero in one state! */ 
        if( i == starti && j == startj ) 
          continue;  
        if( dpenv != NULL && is_in_DPEnvelope(dpenv,i,j) == FALSE )  { /*Is not in envelope*/ 
          EstFrag3_HIDDEN_MATRIX(mat,i,j,MATCH) = NEGI;  
          EstFrag3_HIDDEN_MATRIX(mat,i,j,INSERT) = NEGI;     
          EstFrag3_HIDDEN_MATRIX(mat,i,j,DELETE) = NEGI;     
          continue;  
          } /* end of Is not in envelope */ 


        /* For state MATCH */ 
        /* setting first movement to score */ 
        score = EstFrag3_HIDDEN_MATRIX(mat,i-1,j-3,MATCH) + (mat->query->seg[i]->transition[GW_MATCH2MATCH]+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)]);  
        /* From state INSERT to state MATCH */ 
        temp = EstFrag3_HIDDEN_MATRIX(mat,i-1,j-3,INSERT) + (mat->query->seg[i]->transition[GW_INSERT2MATCH]+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)]);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state MATCH */ 
        temp = EstFrag3_HIDDEN_MATRIX(mat,i-1,j-3,DELETE) + (mat->query->seg[i]->transition[GW_DELETE2MATCH]+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)]);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = EstFrag3_HIDDEN_MATRIX(mat,i-1,j-2,MATCH) + mat->cp->trans[PCD_DELETE_1_BASE];    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = EstFrag3_HIDDEN_MATRIX(mat,i-1,j-1,MATCH) + mat->cp->trans[PCD_DELETE_2_BASE];    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = EstFrag3_HIDDEN_MATRIX(mat,i-1,j-4,MATCH) + mat->cp->trans[PCD_INSERT_1_BASE];    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = EstFrag3_HIDDEN_MATRIX(mat,i-1,j-5,MATCH) + mat->cp->trans[PCD_INSERT_2_BASE];    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state MATCH */ 
        temp = EstFrag3_HIDDEN_MATRIX(mat,i-1,j-2,INSERT) + mat->cp->trans[PCD_DELETE_1_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state MATCH */ 
        temp = EstFrag3_HIDDEN_MATRIX(mat,i-1,j-1,INSERT) + mat->cp->trans[PCD_DELETE_2_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state MATCH */ 
        temp = EstFrag3_HIDDEN_MATRIX(mat,i-1,j-4,INSERT) + mat->cp->trans[PCD_INSERT_1_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state MATCH */ 
        temp = EstFrag3_HIDDEN_MATRIX(mat,i-1,j-5,INSERT) + mat->cp->trans[PCD_INSERT_2_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state MATCH */ 
        temp = EstFrag3_HIDDEN_MATRIX(mat,i-1,j-2,DELETE) + mat->cp->trans[PCD_DELETE_1_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state MATCH */ 
        temp = EstFrag3_HIDDEN_MATRIX(mat,i-1,j-1,DELETE) + mat->cp->trans[PCD_DELETE_2_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state MATCH */ 
        temp = EstFrag3_HIDDEN_MATRIX(mat,i-1,j-4,DELETE) + mat->cp->trans[PCD_INSERT_1_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state MATCH */ 
        temp = EstFrag3_HIDDEN_MATRIX(mat,i-1,j-5,DELETE) + mat->cp->trans[PCD_INSERT_2_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for MATCH */ 
        /* Add any movement independant score and put away */ 
         EstFrag3_HIDDEN_MATRIX(mat,i,j,MATCH) = score;  
        /* Finished calculating state MATCH */ 


        /* For state INSERT */ 
        /* setting first movement to score */ 
        score = EstFrag3_HIDDEN_MATRIX(mat,i-0,j-3,MATCH) + (mat->query->seg[i]->transition[GW_MATCH2INSERT]+mat->query->seg[i]->insert[CSEQ_CDNA_CODON(mat->target,j)]);    
        /* From state INSERT to state INSERT */ 
        temp = EstFrag3_HIDDEN_MATRIX(mat,i-0,j-3,INSERT) + (mat->query->seg[i]->transition[GW_INSERT2INSERT]+mat->query->seg[i]->insert[CSEQ_CDNA_CODON(mat->target,j)]);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state INSERT */ 
        temp = EstFrag3_HIDDEN_MATRIX(mat,i-0,j-3,DELETE) + (mat->query->seg[i]->transition[GW_DELETE2INSERT]+mat->query->seg[i]->insert[CSEQ_CDNA_CODON(mat->target,j)]);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state INSERT */ 
        temp = EstFrag3_HIDDEN_MATRIX(mat,i-0,j-2,INSERT) + mat->cp->trans[PCD_DELETE_1_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state INSERT */ 
        temp = EstFrag3_HIDDEN_MATRIX(mat,i-0,j-1,INSERT) + mat->cp->trans[PCD_DELETE_2_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state INSERT */ 
        temp = EstFrag3_HIDDEN_MATRIX(mat,i-0,j-4,INSERT) + mat->cp->trans[PCD_INSERT_1_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state INSERT */ 
        temp = EstFrag3_HIDDEN_MATRIX(mat,i-0,j-5,INSERT) + mat->cp->trans[PCD_INSERT_2_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state INSERT */ 
        temp = EstFrag3_HIDDEN_MATRIX(mat,i-0,j-2,MATCH) + mat->cp->trans[PCD_DELETE_1_BASE];    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state INSERT */ 
        temp = EstFrag3_HIDDEN_MATRIX(mat,i-0,j-1,MATCH) + mat->cp->trans[PCD_DELETE_2_BASE];    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state INSERT */ 
        temp = EstFrag3_HIDDEN_MATRIX(mat,i-0,j-4,MATCH) + mat->cp->trans[PCD_INSERT_1_BASE];    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state INSERT */ 
        temp = EstFrag3_HIDDEN_MATRIX(mat,i-0,j-5,MATCH) + mat->cp->trans[PCD_INSERT_2_BASE];    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state INSERT */ 
        temp = EstFrag3_HIDDEN_MATRIX(mat,i-0,j-2,DELETE) + mat->cp->trans[PCD_DELETE_1_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state INSERT */ 
        temp = EstFrag3_HIDDEN_MATRIX(mat,i-0,j-1,DELETE) + mat->cp->trans[PCD_DELETE_2_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state INSERT */ 
        temp = EstFrag3_HIDDEN_MATRIX(mat,i-0,j-4,DELETE) + mat->cp->trans[PCD_INSERT_1_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state INSERT */ 
        temp = EstFrag3_HIDDEN_MATRIX(mat,i-0,j-5,DELETE) + mat->cp->trans[PCD_INSERT_2_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for INSERT */ 
        /* Add any movement independant score and put away */ 
         EstFrag3_HIDDEN_MATRIX(mat,i,j,INSERT) = score; 
        /* Finished calculating state INSERT */ 


        /* For state DELETE */ 
        /* setting first movement to score */ 
        score = EstFrag3_HIDDEN_MATRIX(mat,i-1,j-0,MATCH) + mat->query->seg[i]->transition[GW_MATCH2DELETE];     
        /* From state INSERT to state DELETE */ 
        temp = EstFrag3_HIDDEN_MATRIX(mat,i-1,j-0,INSERT) + mat->query->seg[i]->transition[GW_INSERT2DELETE];    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state DELETE */ 
        temp = EstFrag3_HIDDEN_MATRIX(mat,i-1,j-0,DELETE) + mat->query->seg[i]->transition[GW_DELETE2DELETE];    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for DELETE */ 
        /* Add any movement independant score and put away */ 
         EstFrag3_HIDDEN_MATRIX(mat,i,j,DELETE) = score; 
        /* Finished calculating state DELETE */ 
        }  
      }  


    return;  
}    


/* Function:  init_hidden_EstFrag3(mat,starti,startj,stopi,stopj)
 *
 * Descrip: No Description
 *
 * Arg:           mat [UNKN ] Undocumented argument [EstFrag3 *]
 * Arg:        starti [UNKN ] Undocumented argument [int]
 * Arg:        startj [UNKN ] Undocumented argument [int]
 * Arg:         stopi [UNKN ] Undocumented argument [int]
 * Arg:         stopj [UNKN ] Undocumented argument [int]
 *
 */
void init_hidden_EstFrag3(EstFrag3 * mat,int starti,int startj,int stopi,int stopj) 
{
    register int i;  
    register int j;  
    register int hiddenj;    


    hiddenj = startj;    
    for(j=(startj-5);j<=stopj;j++)   {  
      for(i=(starti-1);i<=stopi;i++) {  
        EstFrag3_HIDDEN_MATRIX(mat,i,j,MATCH) = NEGI;
   
        EstFrag3_HIDDEN_MATRIX(mat,i,j,INSERT) = NEGI;
  
        EstFrag3_HIDDEN_MATRIX(mat,i,j,DELETE) = NEGI;
  
        }  
      }  


    return;  
}    


/* Function:  full_dc_EstFrag3(mat,starti,startj,startstate,stopi,stopj,stopstate,out,donej,totalj,dpenv)
 *
 * Descrip:    The main divide-and-conquor routine. Basically, call /PackAln_calculate_small_EstFrag3
 *             Not this function, which is pretty hard core. 
 *             Function is given start/end points (in main matrix) for alignment
 *             It does some checks, decides whether start/end in j is small enough for explicit calc
 *               - if yes, calculates it, reads off into PackAln (out), adds the j distance to donej and returns TRUE
 *               - if no,  uses /do_dc_single_pass_EstFrag3 to get mid-point
 *                          saves midpoint, and calls itself to do right portion then left portion
 *             right then left ensures PackAln is added the 'right' way, ie, back-to-front
 *             returns FALSE on any error, with a warning
 *
 *
 * Arg:               mat [UNKN ] Matrix with small memory implementation [EstFrag3 *]
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
boolean full_dc_EstFrag3(EstFrag3 * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out,int * donej,int totalj,DPEnvelope * dpenv) 
{
    int lstarti; 
    int lstartj; 
    int lstate;  


    if( mat->basematrix->type != BASEMATRIX_TYPE_SHADOW) {  
      warn("*Very* bad error! - non shadow matrix type in full_dc_EstFrag3");    
      return FALSE;  
      }  


    if( starti == -1 || startj == -1 || startstate == -1 || stopi == -1 || stopstate == -1)  {  
      warn("In full dc program, passed bad indices, indices passed were %d:%d[%d] to %d:%d[%d]\n",starti,startj,startstate,stopi,stopj,stopstate);   
      return FALSE;  
      }  


    if( stopj - startj < 25) {  
      log_full_error(REPORT,0,"[%d,%d][%d,%d] Explicit read off",starti,startj,stopi,stopj);/* Build hidden explicit matrix */ 
      calculate_hidden_EstFrag3(mat,starti,startj,startstate,stopi,stopj,stopstate,dpenv);   
      *donej += (stopj - startj);   /* Now read it off into out */ 
      if( read_hidden_EstFrag3(mat,starti,startj,startstate,stopi,stopj,stopstate,out) == FALSE) {  
        warn("In full dc, at %d:%d,%d:%d got a bad hidden explicit read off... ",starti,startj,stopi,stopj); 
        return FALSE;    
        }  
      return TRUE;   
      }  


/* In actual divide and conquor */ 
    if( do_dc_single_pass_EstFrag3(mat,starti,startj,startstate,stopi,stopj,stopstate,dpenv,(int)(*donej*100)/totalj) == FALSE)  {  
      warn("In divide and conquor for EstFrag3, at bound %d:%d to %d:%d, unable to calculate midpoint. Problem!",starti,startj,stopi,stopj); 
      return FALSE;  
      }  


/* Ok... now we have to call on each side of the matrix */ 
/* We have to retrieve left hand side positions, as they will be vapped by the time we call LHS */ 
    lstarti= EstFrag3_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,0);  
    lstartj= EstFrag3_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,1);  
    lstate = EstFrag3_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,2);  


/* Call on right hand side: this lets us do the correct read off */ 
    if( full_dc_EstFrag3(mat,EstFrag3_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,3),EstFrag3_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,4),EstFrag3_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,5),stopi,stopj,stopstate,out,donej,totalj,dpenv) == FALSE)   {  
/* Warning already issued, simply chained back up to top */ 
      return FALSE;  
      }  
/* Call on left hand side */ 
    if( full_dc_EstFrag3(mat,starti,startj,startstate,lstarti,lstartj,lstate,out,donej,totalj,dpenv) == FALSE)   {  
/* Warning already issued, simply chained back up to top */ 
      return FALSE;  
      }  


    return TRUE;     
}    


/* Function:  do_dc_single_pass_EstFrag3(mat,starti,startj,startstate,stopi,stopj,stopstate,dpenv,perc_done)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [EstFrag3 *]
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
boolean do_dc_single_pass_EstFrag3(EstFrag3 * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv,int perc_done) 
{
    int halfj;   
    halfj = startj + ((stopj - startj)/2);   


    init_dc_EstFrag3(mat);   


    EstFrag3_DC_SHADOW_MATRIX(mat,starti,startj,startstate) = 0; 
    run_up_dc_EstFrag3(mat,starti,stopi,startj,halfj-1,dpenv,perc_done);     
    push_dc_at_merge_EstFrag3(mat,starti,stopi,halfj,&halfj,dpenv);  
    follow_on_dc_EstFrag3(mat,starti,stopi,halfj,stopj,dpenv,perc_done);     
    return TRUE; 
}    


/* Function:  push_dc_at_merge_EstFrag3(mat,starti,stopi,startj,stopj,dpenv)
 *
 * Descrip: No Description
 *
 * Arg:           mat [UNKN ] Undocumented argument [EstFrag3 *]
 * Arg:        starti [UNKN ] Undocumented argument [int]
 * Arg:         stopi [UNKN ] Undocumented argument [int]
 * Arg:        startj [UNKN ] Undocumented argument [int]
 * Arg:         stopj [UNKN ] Undocumented argument [int *]
 * Arg:         dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 */
void push_dc_at_merge_EstFrag3(EstFrag3 * mat,int starti,int stopi,int startj,int * stopj,DPEnvelope * dpenv) 
{
    register int i;  
    register int j;  
    register int k;  
    register int count;  
    register int mergej;/* Sources below this j will be stamped by triples */ 
    register int score;  
    register int temp;   


    mergej = startj -1;  
    for(count=0,j=startj;count<5;count++,j++)    {  
      for(i=starti;i<=stopi;i++) {  
        if( dpenv != NULL && is_in_DPEnvelope(dpenv,i,j) == FALSE )  { /*Is not in envelope*/ 
          EstFrag3_DC_SHADOW_MATRIX(mat,i,j,MATCH) = NEGI;   
          EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,0) = (-100);    
          EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,1) = (-100);    
          EstFrag3_DC_SHADOW_MATRIX(mat,i,j,INSERT) = NEGI;  
          EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,0) = (-100);   
          EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,1) = (-100);   
          EstFrag3_DC_SHADOW_MATRIX(mat,i,j,DELETE) = NEGI;  
          EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,0) = (-100);   
          EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,1) = (-100);   
          continue;  
          } /* end of Is not in envelope */ 


        /* For state MATCH, pushing when j - offj <= mergej */ 
        score = EstFrag3_DC_SHADOW_MATRIX(mat,i-1,j-3,MATCH) + (mat->query->seg[i]->transition[GW_MATCH2MATCH]+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)]);   
        if( j - 3 <= mergej) {  
          EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,0) = i-1;   
          EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,1) = j-3;   
          EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,2) = MATCH; 
          EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,3) = i; 
          EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,4) = j; 
          EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,5) = MATCH; 
          }  
        else {  
          for(k=0;k<7;k++)   
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = EstFrag3_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 3,MATCH,k);   
          }  


        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-1,j-3,INSERT) + (mat->query->seg[i]->transition[GW_INSERT2MATCH]+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)]);  
        if( temp > score)    {  
          score = temp;  


          if( j - 3 <= mergej)   {  
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,0) = i-1; 
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,1) = j-3; 
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,2) = INSERT;  
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,3) = i;   
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,4) = j;   
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,5) = MATCH;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = EstFrag3_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 3,INSERT,k);    
            }  
          }  


        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-1,j-3,DELETE) + (mat->query->seg[i]->transition[GW_DELETE2MATCH]+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)]);  
        if( temp > score)    {  
          score = temp;  


          if( j - 3 <= mergej)   {  
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,0) = i-1; 
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,1) = j-3; 
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,2) = DELETE;  
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,3) = i;   
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,4) = j;   
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,5) = MATCH;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = EstFrag3_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 3,DELETE,k);    
            }  
          }  


        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-1,j-2,MATCH) + mat->cp->trans[PCD_DELETE_1_BASE];     
        if( temp > score)    {  
          score = temp;  


          if( j - 2 <= mergej)   {  
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,0) = i-1; 
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,1) = j-2; 
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,2) = MATCH;   
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,3) = i;   
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,4) = j;   
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,5) = MATCH;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = EstFrag3_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 2,MATCH,k); 
            }  
          }  


        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-1,j-1,MATCH) + mat->cp->trans[PCD_DELETE_2_BASE];     
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,0) = i-1; 
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,1) = j-1; 
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,2) = MATCH;   
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,3) = i;   
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,4) = j;   
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,5) = MATCH;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = EstFrag3_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,MATCH,k); 
            }  
          }  


        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-1,j-4,MATCH) + mat->cp->trans[PCD_INSERT_1_BASE];     
        if( temp > score)    {  
          score = temp;  


          if( j - 4 <= mergej)   {  
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,0) = i-1; 
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,1) = j-4; 
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,2) = MATCH;   
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,3) = i;   
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,4) = j;   
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,5) = MATCH;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = EstFrag3_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 4,MATCH,k); 
            }  
          }  


        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-1,j-5,MATCH) + mat->cp->trans[PCD_INSERT_2_BASE];     
        if( temp > score)    {  
          score = temp;  


          if( j - 5 <= mergej)   {  
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,0) = i-1; 
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,1) = j-5; 
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,2) = MATCH;   
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,3) = i;   
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,4) = j;   
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,5) = MATCH;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = EstFrag3_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 5,MATCH,k); 
            }  
          }  


        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-1,j-2,INSERT) + mat->cp->trans[PCD_DELETE_1_BASE];    
        if( temp > score)    {  
          score = temp;  


          if( j - 2 <= mergej)   {  
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,0) = i-1; 
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,1) = j-2; 
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,2) = INSERT;  
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,3) = i;   
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,4) = j;   
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,5) = MATCH;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = EstFrag3_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 2,INSERT,k);    
            }  
          }  


        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-1,j-1,INSERT) + mat->cp->trans[PCD_DELETE_2_BASE];    
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,0) = i-1; 
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,1) = j-1; 
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,2) = INSERT;  
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,3) = i;   
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,4) = j;   
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,5) = MATCH;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = EstFrag3_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,INSERT,k);    
            }  
          }  


        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-1,j-4,INSERT) + mat->cp->trans[PCD_INSERT_1_BASE];    
        if( temp > score)    {  
          score = temp;  


          if( j - 4 <= mergej)   {  
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,0) = i-1; 
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,1) = j-4; 
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,2) = INSERT;  
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,3) = i;   
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,4) = j;   
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,5) = MATCH;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = EstFrag3_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 4,INSERT,k);    
            }  
          }  


        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-1,j-5,INSERT) + mat->cp->trans[PCD_INSERT_2_BASE];    
        if( temp > score)    {  
          score = temp;  


          if( j - 5 <= mergej)   {  
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,0) = i-1; 
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,1) = j-5; 
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,2) = INSERT;  
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,3) = i;   
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,4) = j;   
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,5) = MATCH;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = EstFrag3_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 5,INSERT,k);    
            }  
          }  


        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-1,j-2,DELETE) + mat->cp->trans[PCD_DELETE_1_BASE];    
        if( temp > score)    {  
          score = temp;  


          if( j - 2 <= mergej)   {  
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,0) = i-1; 
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,1) = j-2; 
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,2) = DELETE;  
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,3) = i;   
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,4) = j;   
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,5) = MATCH;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = EstFrag3_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 2,DELETE,k);    
            }  
          }  


        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-1,j-1,DELETE) + mat->cp->trans[PCD_DELETE_2_BASE];    
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,0) = i-1; 
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,1) = j-1; 
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,2) = DELETE;  
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,3) = i;   
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,4) = j;   
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,5) = MATCH;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = EstFrag3_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,DELETE,k);    
            }  
          }  


        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-1,j-4,DELETE) + mat->cp->trans[PCD_INSERT_1_BASE];    
        if( temp > score)    {  
          score = temp;  


          if( j - 4 <= mergej)   {  
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,0) = i-1; 
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,1) = j-4; 
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,2) = DELETE;  
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,3) = i;   
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,4) = j;   
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,5) = MATCH;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = EstFrag3_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 4,DELETE,k);    
            }  
          }  


        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-1,j-5,DELETE) + mat->cp->trans[PCD_INSERT_2_BASE];    
        if( temp > score)    {  
          score = temp;  


          if( j - 5 <= mergej)   {  
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,0) = i-1; 
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,1) = j-5; 
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,2) = DELETE;  
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,3) = i;   
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,4) = j;   
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,5) = MATCH;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = EstFrag3_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 5,DELETE,k);    
            }  
          }  
        /* Add any movement independant score */ 
        EstFrag3_DC_SHADOW_MATRIX(mat,i,j,MATCH) = score;    
        /* Finished with state MATCH */ 


        /* For state INSERT, pushing when j - offj <= mergej */ 
        score = EstFrag3_DC_SHADOW_MATRIX(mat,i-0,j-3,MATCH) + (mat->query->seg[i]->transition[GW_MATCH2INSERT]+mat->query->seg[i]->insert[CSEQ_CDNA_CODON(mat->target,j)]);     
        if( j - 3 <= mergej) {  
          EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,0) = i-0;  
          EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,1) = j-3;  
          EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,2) = MATCH;    
          EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,3) = i;    
          EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,4) = j;    
          EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,5) = INSERT;   
          }  
        else {  
          for(k=0;k<7;k++)   
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,k) = EstFrag3_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 3,MATCH,k);  
          }  


        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-0,j-3,INSERT) + (mat->query->seg[i]->transition[GW_INSERT2INSERT]+mat->query->seg[i]->insert[CSEQ_CDNA_CODON(mat->target,j)]);    
        if( temp > score)    {  
          score = temp;  


          if( j - 3 <= mergej)   {  
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,0) = i-0;    
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,1) = j-3;    
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,2) = INSERT; 
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,3) = i;  
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,4) = j;  
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,5) = INSERT; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,k) = EstFrag3_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 3,INSERT,k);   
            }  
          }  


        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-0,j-3,DELETE) + (mat->query->seg[i]->transition[GW_DELETE2INSERT]+mat->query->seg[i]->insert[CSEQ_CDNA_CODON(mat->target,j)]);    
        if( temp > score)    {  
          score = temp;  


          if( j - 3 <= mergej)   {  
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,0) = i-0;    
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,1) = j-3;    
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,2) = DELETE; 
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,3) = i;  
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,4) = j;  
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,5) = INSERT; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,k) = EstFrag3_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 3,DELETE,k);   
            }  
          }  


        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-0,j-2,INSERT) + mat->cp->trans[PCD_DELETE_1_BASE];    
        if( temp > score)    {  
          score = temp;  


          if( j - 2 <= mergej)   {  
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,0) = i-0;    
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,1) = j-2;    
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,2) = INSERT; 
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,3) = i;  
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,4) = j;  
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,5) = INSERT; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,k) = EstFrag3_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 2,INSERT,k);   
            }  
          }  


        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-0,j-1,INSERT) + mat->cp->trans[PCD_DELETE_2_BASE];    
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,0) = i-0;    
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,1) = j-1;    
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,2) = INSERT; 
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,3) = i;  
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,4) = j;  
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,5) = INSERT; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,k) = EstFrag3_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,INSERT,k);   
            }  
          }  


        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-0,j-4,INSERT) + mat->cp->trans[PCD_INSERT_1_BASE];    
        if( temp > score)    {  
          score = temp;  


          if( j - 4 <= mergej)   {  
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,0) = i-0;    
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,1) = j-4;    
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,2) = INSERT; 
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,3) = i;  
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,4) = j;  
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,5) = INSERT; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,k) = EstFrag3_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 4,INSERT,k);   
            }  
          }  


        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-0,j-5,INSERT) + mat->cp->trans[PCD_INSERT_2_BASE];    
        if( temp > score)    {  
          score = temp;  


          if( j - 5 <= mergej)   {  
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,0) = i-0;    
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,1) = j-5;    
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,2) = INSERT; 
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,3) = i;  
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,4) = j;  
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,5) = INSERT; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,k) = EstFrag3_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 5,INSERT,k);   
            }  
          }  


        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-0,j-2,MATCH) + mat->cp->trans[PCD_DELETE_1_BASE];     
        if( temp > score)    {  
          score = temp;  


          if( j - 2 <= mergej)   {  
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,0) = i-0;    
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,1) = j-2;    
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,2) = MATCH;  
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,3) = i;  
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,4) = j;  
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,5) = INSERT; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,k) = EstFrag3_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 2,MATCH,k);    
            }  
          }  


        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-0,j-1,MATCH) + mat->cp->trans[PCD_DELETE_2_BASE];     
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,0) = i-0;    
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,1) = j-1;    
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,2) = MATCH;  
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,3) = i;  
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,4) = j;  
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,5) = INSERT; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,k) = EstFrag3_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,MATCH,k);    
            }  
          }  


        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-0,j-4,MATCH) + mat->cp->trans[PCD_INSERT_1_BASE];     
        if( temp > score)    {  
          score = temp;  


          if( j - 4 <= mergej)   {  
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,0) = i-0;    
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,1) = j-4;    
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,2) = MATCH;  
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,3) = i;  
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,4) = j;  
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,5) = INSERT; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,k) = EstFrag3_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 4,MATCH,k);    
            }  
          }  


        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-0,j-5,MATCH) + mat->cp->trans[PCD_INSERT_2_BASE];     
        if( temp > score)    {  
          score = temp;  


          if( j - 5 <= mergej)   {  
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,0) = i-0;    
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,1) = j-5;    
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,2) = MATCH;  
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,3) = i;  
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,4) = j;  
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,5) = INSERT; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,k) = EstFrag3_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 5,MATCH,k);    
            }  
          }  


        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-0,j-2,DELETE) + mat->cp->trans[PCD_DELETE_1_BASE];    
        if( temp > score)    {  
          score = temp;  


          if( j - 2 <= mergej)   {  
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,0) = i-0;    
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,1) = j-2;    
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,2) = DELETE; 
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,3) = i;  
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,4) = j;  
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,5) = INSERT; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,k) = EstFrag3_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 2,DELETE,k);   
            }  
          }  


        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-0,j-1,DELETE) + mat->cp->trans[PCD_DELETE_2_BASE];    
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,0) = i-0;    
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,1) = j-1;    
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,2) = DELETE; 
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,3) = i;  
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,4) = j;  
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,5) = INSERT; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,k) = EstFrag3_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,DELETE,k);   
            }  
          }  


        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-0,j-4,DELETE) + mat->cp->trans[PCD_INSERT_1_BASE];    
        if( temp > score)    {  
          score = temp;  


          if( j - 4 <= mergej)   {  
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,0) = i-0;    
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,1) = j-4;    
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,2) = DELETE; 
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,3) = i;  
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,4) = j;  
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,5) = INSERT; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,k) = EstFrag3_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 4,DELETE,k);   
            }  
          }  


        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-0,j-5,DELETE) + mat->cp->trans[PCD_INSERT_2_BASE];    
        if( temp > score)    {  
          score = temp;  


          if( j - 5 <= mergej)   {  
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,0) = i-0;    
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,1) = j-5;    
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,2) = DELETE; 
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,3) = i;  
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,4) = j;  
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,5) = INSERT; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,k) = EstFrag3_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 5,DELETE,k);   
            }  
          }  
        /* Add any movement independant score */ 
        EstFrag3_DC_SHADOW_MATRIX(mat,i,j,INSERT) = score;   
        /* Finished with state INSERT */ 


        /* For state DELETE, pushing when j - offj <= mergej */ 
        score = EstFrag3_DC_SHADOW_MATRIX(mat,i-1,j-0,MATCH) + mat->query->seg[i]->transition[GW_MATCH2DELETE];  
        if( j - 0 <= mergej) {  
          EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,0) = i-1;  
          EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,1) = j-0;  
          EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,2) = MATCH;    
          EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,3) = i;    
          EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,4) = j;    
          EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,5) = DELETE;   
          }  
        else {  
          for(k=0;k<7;k++)   
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,k) = EstFrag3_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 0,MATCH,k);  
          }  


        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-1,j-0,INSERT) + mat->query->seg[i]->transition[GW_INSERT2DELETE];     
        if( temp > score)    {  
          score = temp;  


          if( j - 0 <= mergej)   {  
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,0) = i-1;    
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,1) = j-0;    
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,2) = INSERT; 
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,3) = i;  
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,4) = j;  
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,5) = DELETE; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,k) = EstFrag3_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 0,INSERT,k);   
            }  
          }  


        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-1,j-0,DELETE) + mat->query->seg[i]->transition[GW_DELETE2DELETE];     
        if( temp > score)    {  
          score = temp;  


          if( j - 0 <= mergej)   {  
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,0) = i-1;    
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,1) = j-0;    
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,2) = DELETE; 
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,3) = i;  
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,4) = j;  
            EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,5) = DELETE; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,k) = EstFrag3_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 0,DELETE,k);   
            }  
          }  
        /* Add any movement independant score */ 
        EstFrag3_DC_SHADOW_MATRIX(mat,i,j,DELETE) = score;   
        /* Finished with state DELETE */ 
        }  
      }  
    /* Put back j into * stop j so that calling function gets it correct */ 
    if( stopj == NULL)   
      warn("Bad news... NULL stopj pointer in push dc function. This means that calling function does not know how many cells I have done!");    
    else 
      *stopj = j;    


    return;  
}    


/* Function:  follow_on_dc_EstFrag3(mat,starti,stopi,startj,stopj,dpenv,perc_done)
 *
 * Descrip: No Description
 *
 * Arg:              mat [UNKN ] Undocumented argument [EstFrag3 *]
 * Arg:           starti [UNKN ] Undocumented argument [int]
 * Arg:            stopi [UNKN ] Undocumented argument [int]
 * Arg:           startj [UNKN ] Undocumented argument [int]
 * Arg:            stopj [UNKN ] Undocumented argument [int]
 * Arg:            dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:        perc_done [UNKN ] Undocumented argument [int]
 *
 */
void follow_on_dc_EstFrag3(EstFrag3 * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done) 
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
          EstFrag3_DC_SHADOW_MATRIX(mat,i,j,MATCH) = NEGI;   
          EstFrag3_DC_SHADOW_MATRIX(mat,i,j,INSERT) = NEGI;  
          EstFrag3_DC_SHADOW_MATRIX(mat,i,j,DELETE) = NEGI;  
          continue;  
          } /* end of Is not in envelope */ 
        if( num % 1000 == 0 )    
          log_full_error(REPORT,0,"[%d%%%% done]After  mid-j %5d Cells done %d%%%%",perc_done,startj,(num*100)/total);   


        /* For state MATCH */ 
        /* setting first movement to score */ 
        score = EstFrag3_DC_SHADOW_MATRIX(mat,i-1,j-3,MATCH) + (mat->query->seg[i]->transition[GW_MATCH2MATCH]+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)]);   
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = EstFrag3_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 3,MATCH,k);    
        /* From state INSERT to state MATCH */ 
        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-1,j-3,INSERT) + (mat->query->seg[i]->transition[GW_INSERT2MATCH]+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)]);  
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = EstFrag3_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 3,INSERT,k); 
          }  
        /* From state DELETE to state MATCH */ 
        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-1,j-3,DELETE) + (mat->query->seg[i]->transition[GW_DELETE2MATCH]+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)]);  
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = EstFrag3_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 3,DELETE,k); 
          }  
        /* From state MATCH to state MATCH */ 
        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-1,j-2,MATCH) + mat->cp->trans[PCD_DELETE_1_BASE];     
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = EstFrag3_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 2,MATCH,k);  
          }  
        /* From state MATCH to state MATCH */ 
        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-1,j-1,MATCH) + mat->cp->trans[PCD_DELETE_2_BASE];     
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = EstFrag3_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,MATCH,k);  
          }  
        /* From state MATCH to state MATCH */ 
        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-1,j-4,MATCH) + mat->cp->trans[PCD_INSERT_1_BASE];     
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = EstFrag3_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 4,MATCH,k);  
          }  
        /* From state MATCH to state MATCH */ 
        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-1,j-5,MATCH) + mat->cp->trans[PCD_INSERT_2_BASE];     
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = EstFrag3_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 5,MATCH,k);  
          }  
        /* From state INSERT to state MATCH */ 
        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-1,j-2,INSERT) + mat->cp->trans[PCD_DELETE_1_BASE];    
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = EstFrag3_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 2,INSERT,k); 
          }  
        /* From state INSERT to state MATCH */ 
        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-1,j-1,INSERT) + mat->cp->trans[PCD_DELETE_2_BASE];    
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = EstFrag3_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,INSERT,k); 
          }  
        /* From state INSERT to state MATCH */ 
        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-1,j-4,INSERT) + mat->cp->trans[PCD_INSERT_1_BASE];    
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = EstFrag3_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 4,INSERT,k); 
          }  
        /* From state INSERT to state MATCH */ 
        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-1,j-5,INSERT) + mat->cp->trans[PCD_INSERT_2_BASE];    
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = EstFrag3_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 5,INSERT,k); 
          }  
        /* From state DELETE to state MATCH */ 
        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-1,j-2,DELETE) + mat->cp->trans[PCD_DELETE_1_BASE];    
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = EstFrag3_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 2,DELETE,k); 
          }  
        /* From state DELETE to state MATCH */ 
        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-1,j-1,DELETE) + mat->cp->trans[PCD_DELETE_2_BASE];    
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = EstFrag3_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,DELETE,k); 
          }  
        /* From state DELETE to state MATCH */ 
        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-1,j-4,DELETE) + mat->cp->trans[PCD_INSERT_1_BASE];    
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = EstFrag3_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 4,DELETE,k); 
          }  
        /* From state DELETE to state MATCH */ 
        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-1,j-5,DELETE) + mat->cp->trans[PCD_INSERT_2_BASE];    
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = EstFrag3_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 5,DELETE,k); 
          }  


        /* Ok - finished max calculation for MATCH */ 
        /* Add any movement independant score and put away */ 
         EstFrag3_DC_SHADOW_MATRIX(mat,i,j,MATCH) = score;   
        for(k=0;k<7;k++) 
          EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = localshadow[k];    
        /* Now figure out if any specials need this score */ 
        /* Finished calculating state MATCH */ 


        /* For state INSERT */ 
        /* setting first movement to score */ 
        score = EstFrag3_DC_SHADOW_MATRIX(mat,i-0,j-3,MATCH) + (mat->query->seg[i]->transition[GW_MATCH2INSERT]+mat->query->seg[i]->insert[CSEQ_CDNA_CODON(mat->target,j)]);     
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = EstFrag3_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 3,MATCH,k);    
        /* From state INSERT to state INSERT */ 
        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-0,j-3,INSERT) + (mat->query->seg[i]->transition[GW_INSERT2INSERT]+mat->query->seg[i]->insert[CSEQ_CDNA_CODON(mat->target,j)]);    
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = EstFrag3_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 3,INSERT,k); 
          }  
        /* From state DELETE to state INSERT */ 
        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-0,j-3,DELETE) + (mat->query->seg[i]->transition[GW_DELETE2INSERT]+mat->query->seg[i]->insert[CSEQ_CDNA_CODON(mat->target,j)]);    
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = EstFrag3_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 3,DELETE,k); 
          }  
        /* From state INSERT to state INSERT */ 
        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-0,j-2,INSERT) + mat->cp->trans[PCD_DELETE_1_BASE];    
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = EstFrag3_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 2,INSERT,k); 
          }  
        /* From state INSERT to state INSERT */ 
        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-0,j-1,INSERT) + mat->cp->trans[PCD_DELETE_2_BASE];    
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = EstFrag3_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,INSERT,k); 
          }  
        /* From state INSERT to state INSERT */ 
        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-0,j-4,INSERT) + mat->cp->trans[PCD_INSERT_1_BASE];    
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = EstFrag3_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 4,INSERT,k); 
          }  
        /* From state INSERT to state INSERT */ 
        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-0,j-5,INSERT) + mat->cp->trans[PCD_INSERT_2_BASE];    
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = EstFrag3_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 5,INSERT,k); 
          }  
        /* From state MATCH to state INSERT */ 
        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-0,j-2,MATCH) + mat->cp->trans[PCD_DELETE_1_BASE];     
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = EstFrag3_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 2,MATCH,k);  
          }  
        /* From state MATCH to state INSERT */ 
        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-0,j-1,MATCH) + mat->cp->trans[PCD_DELETE_2_BASE];     
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = EstFrag3_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,MATCH,k);  
          }  
        /* From state MATCH to state INSERT */ 
        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-0,j-4,MATCH) + mat->cp->trans[PCD_INSERT_1_BASE];     
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = EstFrag3_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 4,MATCH,k);  
          }  
        /* From state MATCH to state INSERT */ 
        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-0,j-5,MATCH) + mat->cp->trans[PCD_INSERT_2_BASE];     
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = EstFrag3_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 5,MATCH,k);  
          }  
        /* From state DELETE to state INSERT */ 
        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-0,j-2,DELETE) + mat->cp->trans[PCD_DELETE_1_BASE];    
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = EstFrag3_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 2,DELETE,k); 
          }  
        /* From state DELETE to state INSERT */ 
        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-0,j-1,DELETE) + mat->cp->trans[PCD_DELETE_2_BASE];    
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = EstFrag3_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,DELETE,k); 
          }  
        /* From state DELETE to state INSERT */ 
        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-0,j-4,DELETE) + mat->cp->trans[PCD_INSERT_1_BASE];    
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = EstFrag3_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 4,DELETE,k); 
          }  
        /* From state DELETE to state INSERT */ 
        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-0,j-5,DELETE) + mat->cp->trans[PCD_INSERT_2_BASE];    
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = EstFrag3_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 5,DELETE,k); 
          }  


        /* Ok - finished max calculation for INSERT */ 
        /* Add any movement independant score and put away */ 
         EstFrag3_DC_SHADOW_MATRIX(mat,i,j,INSERT) = score;  
        for(k=0;k<7;k++) 
          EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,k) = localshadow[k];   
        /* Now figure out if any specials need this score */ 
        /* Finished calculating state INSERT */ 


        /* For state DELETE */ 
        /* setting first movement to score */ 
        score = EstFrag3_DC_SHADOW_MATRIX(mat,i-1,j-0,MATCH) + mat->query->seg[i]->transition[GW_MATCH2DELETE];  
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = EstFrag3_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 0,MATCH,k);    
        /* From state INSERT to state DELETE */ 
        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-1,j-0,INSERT) + mat->query->seg[i]->transition[GW_INSERT2DELETE];     
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = EstFrag3_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 0,INSERT,k); 
          }  
        /* From state DELETE to state DELETE */ 
        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-1,j-0,DELETE) + mat->query->seg[i]->transition[GW_DELETE2DELETE];     
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = EstFrag3_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 0,DELETE,k); 
          }  


        /* Ok - finished max calculation for DELETE */ 
        /* Add any movement independant score and put away */ 
         EstFrag3_DC_SHADOW_MATRIX(mat,i,j,DELETE) = score;  
        for(k=0;k<7;k++) 
          EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,k) = localshadow[k];   
        /* Now figure out if any specials need this score */ 
        /* Finished calculating state DELETE */ 
        } /* end of this is strip */ 
      } /* end of for each valid j column */ 


/* Function:  run_up_dc_EstFrag3(mat,starti,stopi,startj,stopj,dpenv,perc_done)
 *
 * Descrip: No Description
 *
 * Arg:              mat [UNKN ] Undocumented argument [EstFrag3 *]
 * Arg:           starti [UNKN ] Undocumented argument [int]
 * Arg:            stopi [UNKN ] Undocumented argument [int]
 * Arg:           startj [UNKN ] Undocumented argument [int]
 * Arg:            stopj [UNKN ] Undocumented argument [int]
 * Arg:            dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:        perc_done [UNKN ] Undocumented argument [int]
 *
 */
}    
void run_up_dc_EstFrag3(EstFrag3 * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done) 
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
          EstFrag3_DC_SHADOW_MATRIX(mat,i,j,MATCH) = NEGI;   
          EstFrag3_DC_SHADOW_MATRIX(mat,i,j,INSERT) = NEGI;  
          EstFrag3_DC_SHADOW_MATRIX(mat,i,j,DELETE) = NEGI;  
          continue;  
          } /* end of Is not in envelope */ 
        if( num % 1000 == 0 )    
          log_full_error(REPORT,0,"[%d%%%% done]Before mid-j %5d Cells done %d%%%%",perc_done,stopj,(num*100)/total);    


        /* For state MATCH */ 
        /* setting first movement to score */ 
        score = EstFrag3_DC_SHADOW_MATRIX(mat,i-1,j-3,MATCH) + (mat->query->seg[i]->transition[GW_MATCH2MATCH]+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)]);   
        /* From state INSERT to state MATCH */ 
        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-1,j-3,INSERT) + (mat->query->seg[i]->transition[GW_INSERT2MATCH]+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)]);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state MATCH */ 
        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-1,j-3,DELETE) + (mat->query->seg[i]->transition[GW_DELETE2MATCH]+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)]);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-1,j-2,MATCH) + mat->cp->trans[PCD_DELETE_1_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-1,j-1,MATCH) + mat->cp->trans[PCD_DELETE_2_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-1,j-4,MATCH) + mat->cp->trans[PCD_INSERT_1_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-1,j-5,MATCH) + mat->cp->trans[PCD_INSERT_2_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state MATCH */ 
        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-1,j-2,INSERT) + mat->cp->trans[PCD_DELETE_1_BASE];    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state MATCH */ 
        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-1,j-1,INSERT) + mat->cp->trans[PCD_DELETE_2_BASE];    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state MATCH */ 
        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-1,j-4,INSERT) + mat->cp->trans[PCD_INSERT_1_BASE];    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state MATCH */ 
        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-1,j-5,INSERT) + mat->cp->trans[PCD_INSERT_2_BASE];    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state MATCH */ 
        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-1,j-2,DELETE) + mat->cp->trans[PCD_DELETE_1_BASE];    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state MATCH */ 
        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-1,j-1,DELETE) + mat->cp->trans[PCD_DELETE_2_BASE];    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state MATCH */ 
        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-1,j-4,DELETE) + mat->cp->trans[PCD_INSERT_1_BASE];    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state MATCH */ 
        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-1,j-5,DELETE) + mat->cp->trans[PCD_INSERT_2_BASE];    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for MATCH */ 
        /* Add any movement independant score and put away */ 
         EstFrag3_DC_SHADOW_MATRIX(mat,i,j,MATCH) = score;   
        /* Finished calculating state MATCH */ 


        /* For state INSERT */ 
        /* setting first movement to score */ 
        score = EstFrag3_DC_SHADOW_MATRIX(mat,i-0,j-3,MATCH) + (mat->query->seg[i]->transition[GW_MATCH2INSERT]+mat->query->seg[i]->insert[CSEQ_CDNA_CODON(mat->target,j)]);     
        /* From state INSERT to state INSERT */ 
        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-0,j-3,INSERT) + (mat->query->seg[i]->transition[GW_INSERT2INSERT]+mat->query->seg[i]->insert[CSEQ_CDNA_CODON(mat->target,j)]);    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state INSERT */ 
        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-0,j-3,DELETE) + (mat->query->seg[i]->transition[GW_DELETE2INSERT]+mat->query->seg[i]->insert[CSEQ_CDNA_CODON(mat->target,j)]);    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state INSERT */ 
        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-0,j-2,INSERT) + mat->cp->trans[PCD_DELETE_1_BASE];    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state INSERT */ 
        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-0,j-1,INSERT) + mat->cp->trans[PCD_DELETE_2_BASE];    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state INSERT */ 
        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-0,j-4,INSERT) + mat->cp->trans[PCD_INSERT_1_BASE];    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state INSERT */ 
        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-0,j-5,INSERT) + mat->cp->trans[PCD_INSERT_2_BASE];    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state INSERT */ 
        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-0,j-2,MATCH) + mat->cp->trans[PCD_DELETE_1_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state INSERT */ 
        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-0,j-1,MATCH) + mat->cp->trans[PCD_DELETE_2_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state INSERT */ 
        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-0,j-4,MATCH) + mat->cp->trans[PCD_INSERT_1_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state INSERT */ 
        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-0,j-5,MATCH) + mat->cp->trans[PCD_INSERT_2_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state INSERT */ 
        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-0,j-2,DELETE) + mat->cp->trans[PCD_DELETE_1_BASE];    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state INSERT */ 
        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-0,j-1,DELETE) + mat->cp->trans[PCD_DELETE_2_BASE];    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state INSERT */ 
        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-0,j-4,DELETE) + mat->cp->trans[PCD_INSERT_1_BASE];    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state INSERT */ 
        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-0,j-5,DELETE) + mat->cp->trans[PCD_INSERT_2_BASE];    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for INSERT */ 
        /* Add any movement independant score and put away */ 
         EstFrag3_DC_SHADOW_MATRIX(mat,i,j,INSERT) = score;  
        /* Finished calculating state INSERT */ 


        /* For state DELETE */ 
        /* setting first movement to score */ 
        score = EstFrag3_DC_SHADOW_MATRIX(mat,i-1,j-0,MATCH) + mat->query->seg[i]->transition[GW_MATCH2DELETE];  
        /* From state INSERT to state DELETE */ 
        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-1,j-0,INSERT) + mat->query->seg[i]->transition[GW_INSERT2DELETE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state DELETE */ 
        temp = EstFrag3_DC_SHADOW_MATRIX(mat,i-1,j-0,DELETE) + mat->query->seg[i]->transition[GW_DELETE2DELETE];     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for DELETE */ 
        /* Add any movement independant score and put away */ 
         EstFrag3_DC_SHADOW_MATRIX(mat,i,j,DELETE) = score;  
        /* Finished calculating state DELETE */ 
        } /* end of this is strip */ 
      } /* end of for each valid j column */ 


/* Function:  init_dc_EstFrag3(mat)
 *
 * Descrip: No Description
 *
 * Arg:        mat [UNKN ] Undocumented argument [EstFrag3 *]
 *
 */
}    
void init_dc_EstFrag3(EstFrag3 * mat) 
{
    register int i;  
    register int j;  
    register int k;  


    for(j=0;j<7;j++) {  
      for(i=(-1);i<mat->query->len;i++)  {  
        EstFrag3_DC_SHADOW_MATRIX(mat,i,j,MATCH) = NEGI; 
        EstFrag3_DC_SHADOW_MATRIX(mat,i,j,INSERT) = NEGI;    
        EstFrag3_DC_SHADOW_MATRIX(mat,i,j,DELETE) = NEGI;    
        for(k=0;k<7;k++) {  
          EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = (-1);  
          EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,k) = (-1); 
          EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,k) = (-1); 
          }  
        }  
      }  


    return;  
}    


/* Function:  start_end_find_end_EstFrag3(mat,endj)
 *
 * Descrip:    First function used to find end of the best path in the special state !end
 *
 *
 * Arg:         mat [UNKN ] Matrix in small mode [EstFrag3 *]
 * Arg:        endj [WRITE] position of end in j (meaningless in i) [int *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int start_end_find_end_EstFrag3(EstFrag3 * mat,int * endj) 
{
    register int j;  
    register int max;    
    register int maxj;   


    max = EstFrag3_DC_SHADOW_SPECIAL(mat,0,mat->target->seq->len-1,END); 
    maxj = mat->target->seq->len-1;  
    for(j= mat->target->seq->len-2 ;j >= 0 ;j--) {  
      if( EstFrag3_DC_SHADOW_SPECIAL(mat,0,j,END) > max )    {  
        max = EstFrag3_DC_SHADOW_SPECIAL(mat,0,j,END);   
        maxj = j;    
        }  
      }  


    if( endj != NULL)    
      *endj = maxj;  


    return max;  
}    


/* Function:  dc_optimised_start_end_calc_EstFrag3(*mat,dpenv)
 *
 * Descrip:    Calculates special strip, leaving start/end/score points in shadow matrix
 *             Works off specially laid out memory from steve searle
 *
 *
 * Arg:         *mat [UNKN ] Undocumented argument [EstFrag3]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean dc_optimised_start_end_calc_EstFrag3(EstFrag3 *mat,DPEnvelope * dpenv) 
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
    leni = mat->query->len;  
    lenj = mat->target->seq->len;    
    total = leni * lenj; 


    score_pointers = (int *) calloc (5 * (leni + 1) * 3,sizeof(int));    
    shadow_pointers = (int *) calloc (5 * (leni + 1) * 3 * 8,sizeof(int));   


    for(j=0;j<lenj;j++)  { /*for each j strip*/ 
      for(i=0;i<leni;i++)    { /*for each i position in strip*/ 
        num++;   
        if( dpenv != NULL && is_in_DPEnvelope(dpenv,i,j) == FALSE )  { /*Is not in envelope*/ 
          EstFrag3_DC_OPT_SHADOW_MATRIX(mat,i,j,MATCH) = NEGI;   
          EstFrag3_DC_OPT_SHADOW_MATRIX(mat,i,j,INSERT) = NEGI;  
          EstFrag3_DC_OPT_SHADOW_MATRIX(mat,i,j,DELETE) = NEGI;  
          continue;  
          } /* end of Is not in envelope */ 
        if( num%1000 == 0)   
          log_full_error(REPORT,0,"%6d Cells done [%2d%%%%]",num,num*100/total); 




        /* For state MATCH */ 
        /* setting first movement to score */ 
        score = EstFrag3_DC_OPT_SHADOW_MATRIX(mat,i-1,j-3,MATCH) + (mat->query->seg[i]->transition[GW_MATCH2MATCH]+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)]) + (0);     
        /* assign local shadown pointer */ 
        localsp = &(EstFrag3_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 3,MATCH,0));  
        /* From state INSERT to state MATCH */ 
        temp = EstFrag3_DC_OPT_SHADOW_MATRIX(mat,i-1,j-3,INSERT) + (mat->query->seg[i]->transition[GW_INSERT2MATCH]+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)]) +(0);     
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(EstFrag3_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 3,INSERT,0));   
          }  
        /* From state DELETE to state MATCH */ 
        temp = EstFrag3_DC_OPT_SHADOW_MATRIX(mat,i-1,j-3,DELETE) + (mat->query->seg[i]->transition[GW_DELETE2MATCH]+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)]) +(0);     
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(EstFrag3_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 3,DELETE,0));   
          }  
        /* From state START to state MATCH */ 
        temp = EstFrag3_DC_OPT_SHADOW_SPECIAL(mat,i-1,j-3,START) + (mat->query->seg[i]->transition[GW_START2MATCH]+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)]) + (0);     
        if( temp  > score )  {  
          score = temp;  
          /* This state [START] is a special for MATCH... push top shadow pointers here */ 
          localshadow[0]= i; 
          localshadow[1]= j; 
          localshadow[2]= MATCH; 
          localshadow[3]= (-1);  
          localshadow[4]= (-1);  
          localshadow[5]= (-1);  
          localshadow[6]= score; 
          localsp = localshadow; 
          }  
        /* From state START to state MATCH */ 
        temp = EstFrag3_DC_OPT_SHADOW_SPECIAL(mat,i-1,j-3,START) + (mat->start_frag+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)]) + (0);    
        if( temp  > score )  {  
          score = temp;  
          /* This state [START] is a special for MATCH... push top shadow pointers here */ 
          localshadow[0]= i; 
          localshadow[1]= j; 
          localshadow[2]= MATCH; 
          localshadow[3]= (-1);  
          localshadow[4]= (-1);  
          localshadow[5]= (-1);  
          localshadow[6]= score; 
          localsp = localshadow; 
          }  
        /* From state MATCH to state MATCH */ 
        temp = EstFrag3_DC_OPT_SHADOW_MATRIX(mat,i-1,j-2,MATCH) + mat->cp->trans[PCD_DELETE_1_BASE] +(0);    
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(EstFrag3_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 2,MATCH,0));    
          }  
        /* From state MATCH to state MATCH */ 
        temp = EstFrag3_DC_OPT_SHADOW_MATRIX(mat,i-1,j-1,MATCH) + mat->cp->trans[PCD_DELETE_2_BASE] +(0);    
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(EstFrag3_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 1,MATCH,0));    
          }  
        /* From state MATCH to state MATCH */ 
        temp = EstFrag3_DC_OPT_SHADOW_MATRIX(mat,i-1,j-4,MATCH) + mat->cp->trans[PCD_INSERT_1_BASE] +(0);    
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(EstFrag3_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 4,MATCH,0));    
          }  
        /* From state MATCH to state MATCH */ 
        temp = EstFrag3_DC_OPT_SHADOW_MATRIX(mat,i-1,j-5,MATCH) + mat->cp->trans[PCD_INSERT_2_BASE] +(0);    
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(EstFrag3_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 5,MATCH,0));    
          }  
        /* From state INSERT to state MATCH */ 
        temp = EstFrag3_DC_OPT_SHADOW_MATRIX(mat,i-1,j-2,INSERT) + mat->cp->trans[PCD_DELETE_1_BASE] +(0);   
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(EstFrag3_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 2,INSERT,0));   
          }  
        /* From state INSERT to state MATCH */ 
        temp = EstFrag3_DC_OPT_SHADOW_MATRIX(mat,i-1,j-1,INSERT) + mat->cp->trans[PCD_DELETE_2_BASE] +(0);   
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(EstFrag3_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 1,INSERT,0));   
          }  
        /* From state INSERT to state MATCH */ 
        temp = EstFrag3_DC_OPT_SHADOW_MATRIX(mat,i-1,j-4,INSERT) + mat->cp->trans[PCD_INSERT_1_BASE] +(0);   
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(EstFrag3_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 4,INSERT,0));   
          }  
        /* From state INSERT to state MATCH */ 
        temp = EstFrag3_DC_OPT_SHADOW_MATRIX(mat,i-1,j-5,INSERT) + mat->cp->trans[PCD_INSERT_2_BASE] +(0);   
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(EstFrag3_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 5,INSERT,0));   
          }  
        /* From state DELETE to state MATCH */ 
        temp = EstFrag3_DC_OPT_SHADOW_MATRIX(mat,i-1,j-2,DELETE) + mat->cp->trans[PCD_DELETE_1_BASE] +(0);   
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(EstFrag3_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 2,DELETE,0));   
          }  
        /* From state DELETE to state MATCH */ 
        temp = EstFrag3_DC_OPT_SHADOW_MATRIX(mat,i-1,j-1,DELETE) + mat->cp->trans[PCD_DELETE_2_BASE] +(0);   
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(EstFrag3_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 1,DELETE,0));   
          }  
        /* From state DELETE to state MATCH */ 
        temp = EstFrag3_DC_OPT_SHADOW_MATRIX(mat,i-1,j-4,DELETE) + mat->cp->trans[PCD_INSERT_1_BASE] +(0);   
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(EstFrag3_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 4,DELETE,0));   
          }  
        /* From state DELETE to state MATCH */ 
        temp = EstFrag3_DC_OPT_SHADOW_MATRIX(mat,i-1,j-5,DELETE) + mat->cp->trans[PCD_INSERT_2_BASE] +(0);   
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(EstFrag3_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 5,DELETE,0));   
          }  


        /* Ok - finished max calculation for MATCH */ 
        /* Add any movement independant score and put away */ 
        /* Actually, already done inside scores */ 
         EstFrag3_DC_OPT_SHADOW_MATRIX(mat,i,j,MATCH) = score;   
        for(k=0;k<7;k++) 
          EstFrag3_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = localsp[k];    
        /* Now figure out if any specials need this score */ 


        /* state MATCH is a source for special END */ 
        temp = score + (mat->query->seg[i]->transition[GW_MATCH2END]) + (0) ;    
        if( temp > EstFrag3_DC_OPT_SHADOW_SPECIAL(mat,i,j,END) )     {  
          EstFrag3_DC_OPT_SHADOW_SPECIAL(mat,i,j,END) = temp;    
          /* Have to push only bottem half of system here */ 
          for(k=0;k<3;k++)   
            EstFrag3_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,k) = EstFrag3_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,MATCH,k);    
          EstFrag3_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,6) = EstFrag3_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,MATCH,6);  
          EstFrag3_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,3) = i;  
          EstFrag3_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,4) = j;  
          EstFrag3_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,5) = MATCH;  
          }  




        /* state MATCH is a source for special END */ 
        temp = score + (mat->end_frag) + (0) ;   
        if( temp > EstFrag3_DC_OPT_SHADOW_SPECIAL(mat,i,j,END) )     {  
          EstFrag3_DC_OPT_SHADOW_SPECIAL(mat,i,j,END) = temp;    
          /* Have to push only bottem half of system here */ 
          for(k=0;k<3;k++)   
            EstFrag3_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,k) = EstFrag3_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,MATCH,k);    
          EstFrag3_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,6) = EstFrag3_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,MATCH,6);  
          EstFrag3_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,3) = i;  
          EstFrag3_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,4) = j;  
          EstFrag3_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,5) = MATCH;  
          }  




        /* Finished calculating state MATCH */ 


        /* For state INSERT */ 
        /* setting first movement to score */ 
        score = EstFrag3_DC_OPT_SHADOW_MATRIX(mat,i-0,j-3,MATCH) + (mat->query->seg[i]->transition[GW_MATCH2INSERT]+mat->query->seg[i]->insert[CSEQ_CDNA_CODON(mat->target,j)]) + (0);   
        /* assign local shadown pointer */ 
        localsp = &(EstFrag3_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 3,MATCH,0));  
        /* From state INSERT to state INSERT */ 
        temp = EstFrag3_DC_OPT_SHADOW_MATRIX(mat,i-0,j-3,INSERT) + (mat->query->seg[i]->transition[GW_INSERT2INSERT]+mat->query->seg[i]->insert[CSEQ_CDNA_CODON(mat->target,j)]) +(0);   
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(EstFrag3_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 3,INSERT,0));   
          }  
        /* From state DELETE to state INSERT */ 
        temp = EstFrag3_DC_OPT_SHADOW_MATRIX(mat,i-0,j-3,DELETE) + (mat->query->seg[i]->transition[GW_DELETE2INSERT]+mat->query->seg[i]->insert[CSEQ_CDNA_CODON(mat->target,j)]) +(0);   
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(EstFrag3_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 3,DELETE,0));   
          }  
        /* From state START to state INSERT */ 
        temp = EstFrag3_DC_OPT_SHADOW_SPECIAL(mat,i-0,j-3,START) + (mat->query->seg[i]->transition[GW_START2INSERT]+mat->query->seg[i]->insert[CSEQ_CDNA_CODON(mat->target,j)]) + (0);   
        if( temp  > score )  {  
          score = temp;  
          /* This state [START] is a special for INSERT... push top shadow pointers here */ 
          localshadow[0]= i; 
          localshadow[1]= j; 
          localshadow[2]= INSERT;    
          localshadow[3]= (-1);  
          localshadow[4]= (-1);  
          localshadow[5]= (-1);  
          localshadow[6]= score; 
          localsp = localshadow; 
          }  
        /* From state INSERT to state INSERT */ 
        temp = EstFrag3_DC_OPT_SHADOW_MATRIX(mat,i-0,j-2,INSERT) + mat->cp->trans[PCD_DELETE_1_BASE] +(0);   
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(EstFrag3_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 2,INSERT,0));   
          }  
        /* From state INSERT to state INSERT */ 
        temp = EstFrag3_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,INSERT) + mat->cp->trans[PCD_DELETE_2_BASE] +(0);   
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(EstFrag3_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,INSERT,0));   
          }  
        /* From state INSERT to state INSERT */ 
        temp = EstFrag3_DC_OPT_SHADOW_MATRIX(mat,i-0,j-4,INSERT) + mat->cp->trans[PCD_INSERT_1_BASE] +(0);   
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(EstFrag3_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 4,INSERT,0));   
          }  
        /* From state INSERT to state INSERT */ 
        temp = EstFrag3_DC_OPT_SHADOW_MATRIX(mat,i-0,j-5,INSERT) + mat->cp->trans[PCD_INSERT_2_BASE] +(0);   
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(EstFrag3_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 5,INSERT,0));   
          }  
        /* From state MATCH to state INSERT */ 
        temp = EstFrag3_DC_OPT_SHADOW_MATRIX(mat,i-0,j-2,MATCH) + mat->cp->trans[PCD_DELETE_1_BASE] +(0);    
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(EstFrag3_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 2,MATCH,0));    
          }  
        /* From state MATCH to state INSERT */ 
        temp = EstFrag3_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,MATCH) + mat->cp->trans[PCD_DELETE_2_BASE] +(0);    
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(EstFrag3_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,MATCH,0));    
          }  
        /* From state MATCH to state INSERT */ 
        temp = EstFrag3_DC_OPT_SHADOW_MATRIX(mat,i-0,j-4,MATCH) + mat->cp->trans[PCD_INSERT_1_BASE] +(0);    
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(EstFrag3_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 4,MATCH,0));    
          }  
        /* From state MATCH to state INSERT */ 
        temp = EstFrag3_DC_OPT_SHADOW_MATRIX(mat,i-0,j-5,MATCH) + mat->cp->trans[PCD_INSERT_2_BASE] +(0);    
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(EstFrag3_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 5,MATCH,0));    
          }  
        /* From state DELETE to state INSERT */ 
        temp = EstFrag3_DC_OPT_SHADOW_MATRIX(mat,i-0,j-2,DELETE) + mat->cp->trans[PCD_DELETE_1_BASE] +(0);   
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(EstFrag3_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 2,DELETE,0));   
          }  
        /* From state DELETE to state INSERT */ 
        temp = EstFrag3_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,DELETE) + mat->cp->trans[PCD_DELETE_2_BASE] +(0);   
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(EstFrag3_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,DELETE,0));   
          }  
        /* From state DELETE to state INSERT */ 
        temp = EstFrag3_DC_OPT_SHADOW_MATRIX(mat,i-0,j-4,DELETE) + mat->cp->trans[PCD_INSERT_1_BASE] +(0);   
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(EstFrag3_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 4,DELETE,0));   
          }  
        /* From state DELETE to state INSERT */ 
        temp = EstFrag3_DC_OPT_SHADOW_MATRIX(mat,i-0,j-5,DELETE) + mat->cp->trans[PCD_INSERT_2_BASE] +(0);   
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(EstFrag3_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 5,DELETE,0));   
          }  


        /* Ok - finished max calculation for INSERT */ 
        /* Add any movement independant score and put away */ 
        /* Actually, already done inside scores */ 
         EstFrag3_DC_OPT_SHADOW_MATRIX(mat,i,j,INSERT) = score;  
        for(k=0;k<7;k++) 
          EstFrag3_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,INSERT,k) = localsp[k];   
        /* Now figure out if any specials need this score */ 


        /* state INSERT is a source for special END */ 
        temp = score + (mat->query->seg[i]->transition[GW_INSERT2END]) + (0) ;   
        if( temp > EstFrag3_DC_OPT_SHADOW_SPECIAL(mat,i,j,END) )     {  
          EstFrag3_DC_OPT_SHADOW_SPECIAL(mat,i,j,END) = temp;    
          /* Have to push only bottem half of system here */ 
          for(k=0;k<3;k++)   
            EstFrag3_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,k) = EstFrag3_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,INSERT,k);   
          EstFrag3_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,6) = EstFrag3_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,INSERT,6); 
          EstFrag3_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,3) = i;  
          EstFrag3_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,4) = j;  
          EstFrag3_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,5) = INSERT; 
          }  




        /* Finished calculating state INSERT */ 


        /* For state DELETE */ 
        /* setting first movement to score */ 
        score = EstFrag3_DC_OPT_SHADOW_MATRIX(mat,i-1,j-0,MATCH) + mat->query->seg[i]->transition[GW_MATCH2DELETE] + (0);    
        /* assign local shadown pointer */ 
        localsp = &(EstFrag3_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 0,MATCH,0));  
        /* From state INSERT to state DELETE */ 
        temp = EstFrag3_DC_OPT_SHADOW_MATRIX(mat,i-1,j-0,INSERT) + mat->query->seg[i]->transition[GW_INSERT2DELETE] +(0);    
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(EstFrag3_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 0,INSERT,0));   
          }  
        /* From state DELETE to state DELETE */ 
        temp = EstFrag3_DC_OPT_SHADOW_MATRIX(mat,i-1,j-0,DELETE) + mat->query->seg[i]->transition[GW_DELETE2DELETE] +(0);    
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(EstFrag3_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 0,DELETE,0));   
          }  
        /* From state START to state DELETE */ 
        temp = EstFrag3_DC_OPT_SHADOW_SPECIAL(mat,i-1,j-0,START) + mat->query->seg[i]->transition[GW_START2DELETE] + (0);    
        if( temp  > score )  {  
          score = temp;  
          /* This state [START] is a special for DELETE... push top shadow pointers here */ 
          localshadow[0]= i; 
          localshadow[1]= j; 
          localshadow[2]= DELETE;    
          localshadow[3]= (-1);  
          localshadow[4]= (-1);  
          localshadow[5]= (-1);  
          localshadow[6]= score; 
          localsp = localshadow; 
          }  


        /* Ok - finished max calculation for DELETE */ 
        /* Add any movement independant score and put away */ 
        /* Actually, already done inside scores */ 
         EstFrag3_DC_OPT_SHADOW_MATRIX(mat,i,j,DELETE) = score;  
        for(k=0;k<7;k++) 
          EstFrag3_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,DELETE,k) = localsp[k];   
        /* Now figure out if any specials need this score */ 


        /* state DELETE is a source for special END */ 
        temp = score + (mat->query->seg[i]->transition[GW_DELETE2END]) + (0) ;   
        if( temp > EstFrag3_DC_OPT_SHADOW_SPECIAL(mat,i,j,END) )     {  
          EstFrag3_DC_OPT_SHADOW_SPECIAL(mat,i,j,END) = temp;    
          /* Have to push only bottem half of system here */ 
          for(k=0;k<3;k++)   
            EstFrag3_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,k) = EstFrag3_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,DELETE,k);   
          EstFrag3_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,6) = EstFrag3_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,DELETE,6); 
          EstFrag3_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,3) = i;  
          EstFrag3_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,4) = j;  
          EstFrag3_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,5) = DELETE; 
          }  




        /* Finished calculating state DELETE */ 


        } /* end of for each i position in strip */ 
      } /* end of for each j strip */ 
    free(score_pointers);    
    free(shadow_pointers);   
    return TRUE;     
}    


/* Function:  init_start_end_linear_EstFrag3(mat)
 *
 * Descrip: No Description
 *
 * Arg:        mat [UNKN ] Undocumented argument [EstFrag3 *]
 *
 */
void init_start_end_linear_EstFrag3(EstFrag3 * mat) 
{
    register int i;  
    register int j;  
    for(j=0;j<7;j++) {  
      for(i=(-1);i<mat->query->len;i++)  {  
        EstFrag3_DC_SHADOW_MATRIX(mat,i,j,MATCH) = NEGI; 
        EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,0) = (-1);    
        EstFrag3_DC_SHADOW_MATRIX(mat,i,j,INSERT) = NEGI;    
        EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,0) = (-1);   
        EstFrag3_DC_SHADOW_MATRIX(mat,i,j,DELETE) = NEGI;    
        EstFrag3_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,0) = (-1);   
        }  
      }  


    for(j=(-5);j<mat->target->seq->len;j++)  {  
      EstFrag3_DC_SHADOW_SPECIAL(mat,0,j,START) = 0; 
      EstFrag3_DC_SHADOW_SPECIAL_SP(mat,0,j,START,0) = j;    
      EstFrag3_DC_SHADOW_SPECIAL(mat,0,j,END) = NEGI;    
      EstFrag3_DC_SHADOW_SPECIAL_SP(mat,0,j,END,0) = (-1);   
      }  


    return;  
}    


/* Function:  convert_PackAln_to_AlnBlock_EstFrag3(pal)
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
AlnBlock * convert_PackAln_to_AlnBlock_EstFrag3(PackAln * pal) 
{
    AlnConvertSet * acs; 
    AlnBlock * alb;  


    acs = AlnConvertSet_EstFrag3();  
    alb = AlnBlock_from_PackAln(acs,pal);    
    free_AlnConvertSet(acs); 
    return alb;  
}    


 static char * query_label[] = { "MATCH_STATE","INSERT_STATE","DELETE_STATE","END" };    
/* Function:  AlnConvertSet_EstFrag3(void)
 *
 * Descrip: No Description
 *
 *
 * Return [UNKN ]  Undocumented return value [AlnConvertSet *]
 *
 */
 static char * target_label[] = { "CODON","SEQUENCE_DELETION","SEQUENCE_INSERTION","INSERT","END" }; 
AlnConvertSet * AlnConvertSet_EstFrag3(void) 
{
    AlnConvertUnit * acu;    
    AlnConvertSet  * out;    


    out = AlnConvertSet_alloc_std(); 


    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = MATCH; 
    acu->state2 = MATCH;     
    acu->offi = 1;   
    acu->offj = 3;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[0];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = INSERT;    
    acu->state2 = MATCH;     
    acu->offi = 1;   
    acu->offj = 3;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[0];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = DELETE;    
    acu->state2 = MATCH;     
    acu->offi = 1;   
    acu->offj = 3;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[0];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = START + 3; 
    acu->is_from_special = TRUE; 
    acu->state2 = MATCH;     
    acu->offi = (-1);    
    acu->offj = 3;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[0];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = START + 3; 
    acu->is_from_special = TRUE; 
    acu->state2 = MATCH;     
    acu->offi = (-1);    
    acu->offj = 3;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[0];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = MATCH; 
    acu->state2 = MATCH;     
    acu->offi = 1;   
    acu->offj = 2;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[1];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = MATCH; 
    acu->state2 = MATCH;     
    acu->offi = 1;   
    acu->offj = 1;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[1];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = MATCH; 
    acu->state2 = MATCH;     
    acu->offi = 1;   
    acu->offj = 4;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[2];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = MATCH; 
    acu->state2 = MATCH;     
    acu->offi = 1;   
    acu->offj = 5;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[2];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = INSERT;    
    acu->state2 = MATCH;     
    acu->offi = 1;   
    acu->offj = 2;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[1];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = INSERT;    
    acu->state2 = MATCH;     
    acu->offi = 1;   
    acu->offj = 1;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[1];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = INSERT;    
    acu->state2 = MATCH;     
    acu->offi = 1;   
    acu->offj = 4;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[2];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = INSERT;    
    acu->state2 = MATCH;     
    acu->offi = 1;   
    acu->offj = 5;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[2];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = DELETE;    
    acu->state2 = MATCH;     
    acu->offi = 1;   
    acu->offj = 2;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[1];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = DELETE;    
    acu->state2 = MATCH;     
    acu->offi = 1;   
    acu->offj = 1;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[1];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = DELETE;    
    acu->state2 = MATCH;     
    acu->offi = 1;   
    acu->offj = 4;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[2];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = DELETE;    
    acu->state2 = MATCH;     
    acu->offi = 1;   
    acu->offj = 5;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[2];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = MATCH; 
    acu->state2 = INSERT;    
    acu->offi = 0;   
    acu->offj = 3;   
    acu->label1 = query_label[1];    
    acu->label2 = target_label[0];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = INSERT;    
    acu->state2 = INSERT;    
    acu->offi = 0;   
    acu->offj = 3;   
    acu->label1 = query_label[1];    
    acu->label2 = target_label[0];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = DELETE;    
    acu->state2 = INSERT;    
    acu->offi = 0;   
    acu->offj = 3;   
    acu->label1 = query_label[1];    
    acu->label2 = target_label[0];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = START + 3; 
    acu->is_from_special = TRUE; 
    acu->state2 = INSERT;    
    acu->offi = (-1);    
    acu->offj = 3;   
    acu->label1 = query_label[1];    
    acu->label2 = target_label[0];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = INSERT;    
    acu->state2 = INSERT;    
    acu->offi = 0;   
    acu->offj = 2;   
    acu->label1 = query_label[1];    
    acu->label2 = target_label[1];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = INSERT;    
    acu->state2 = INSERT;    
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[1];    
    acu->label2 = target_label[1];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = INSERT;    
    acu->state2 = INSERT;    
    acu->offi = 0;   
    acu->offj = 4;   
    acu->label1 = query_label[1];    
    acu->label2 = target_label[2];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = INSERT;    
    acu->state2 = INSERT;    
    acu->offi = 0;   
    acu->offj = 5;   
    acu->label1 = query_label[1];    
    acu->label2 = target_label[2];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = MATCH; 
    acu->state2 = INSERT;    
    acu->offi = 0;   
    acu->offj = 2;   
    acu->label1 = query_label[1];    
    acu->label2 = target_label[1];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = MATCH; 
    acu->state2 = INSERT;    
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[1];    
    acu->label2 = target_label[1];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = MATCH; 
    acu->state2 = INSERT;    
    acu->offi = 0;   
    acu->offj = 4;   
    acu->label1 = query_label[1];    
    acu->label2 = target_label[2];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = MATCH; 
    acu->state2 = INSERT;    
    acu->offi = 0;   
    acu->offj = 5;   
    acu->label1 = query_label[1];    
    acu->label2 = target_label[2];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = DELETE;    
    acu->state2 = INSERT;    
    acu->offi = 0;   
    acu->offj = 2;   
    acu->label1 = query_label[1];    
    acu->label2 = target_label[1];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = DELETE;    
    acu->state2 = INSERT;    
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[1];    
    acu->label2 = target_label[1];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = DELETE;    
    acu->state2 = INSERT;    
    acu->offi = 0;   
    acu->offj = 4;   
    acu->label1 = query_label[1];    
    acu->label2 = target_label[2];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = DELETE;    
    acu->state2 = INSERT;    
    acu->offi = 0;   
    acu->offj = 5;   
    acu->label1 = query_label[1];    
    acu->label2 = target_label[2];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = MATCH; 
    acu->state2 = DELETE;    
    acu->offi = 1;   
    acu->offj = 0;   
    acu->label1 = query_label[2];    
    acu->label2 = target_label[3];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = INSERT;    
    acu->state2 = DELETE;    
    acu->offi = 1;   
    acu->offj = 0;   
    acu->label1 = query_label[2];    
    acu->label2 = target_label[3];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = DELETE;    
    acu->state2 = DELETE;    
    acu->offi = 1;   
    acu->offj = 0;   
    acu->label1 = query_label[2];    
    acu->label2 = target_label[3];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = START + 3; 
    acu->is_from_special = TRUE; 
    acu->state2 = DELETE;    
    acu->offi = (-1);    
    acu->offj = 0;   
    acu->label1 = query_label[2];    
    acu->label2 = target_label[3];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = MATCH; 
    acu->state2 = END + 3;   
    acu->offi = (-1);    
    acu->offj = 0;   
    acu->label1 = query_label[3];    
    acu->label2 = target_label[4];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = MATCH; 
    acu->state2 = END + 3;   
    acu->offi = (-1);    
    acu->offj = 0;   
    acu->label1 = query_label[3];    
    acu->label2 = target_label[4];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = INSERT;    
    acu->state2 = END + 3;   
    acu->offi = (-1);    
    acu->offj = 0;   
    acu->label1 = query_label[3];    
    acu->label2 = target_label[4];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = DELETE;    
    acu->state2 = END + 3;   
    acu->offi = (-1);    
    acu->offj = 0;   
    acu->label1 = query_label[3];    
    acu->label2 = target_label[4];   
    return out;  
}    


/* Function:  PackAln_read_Expl_EstFrag3(mat)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [EstFrag3 *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * PackAln_read_Expl_EstFrag3(EstFrag3 * mat) 
{
    EstFrag3_access_func_holder holder;  


    holder.access_main    = EstFrag3_explicit_access_main;   
    holder.access_special = EstFrag3_explicit_access_special;    
    return PackAln_read_generic_EstFrag3(mat,holder);    
}    


/* Function:  EstFrag3_explicit_access_main(mat,i,j,state)
 *
 * Descrip: No Description
 *
 * Arg:          mat [UNKN ] Undocumented argument [EstFrag3 *]
 * Arg:            i [UNKN ] Undocumented argument [int]
 * Arg:            j [UNKN ] Undocumented argument [int]
 * Arg:        state [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int EstFrag3_explicit_access_main(EstFrag3 * mat,int i,int j,int state) 
{
    return EstFrag3_EXPL_MATRIX(mat,i,j,state);  
}    


/* Function:  EstFrag3_explicit_access_special(mat,i,j,state)
 *
 * Descrip: No Description
 *
 * Arg:          mat [UNKN ] Undocumented argument [EstFrag3 *]
 * Arg:            i [UNKN ] Undocumented argument [int]
 * Arg:            j [UNKN ] Undocumented argument [int]
 * Arg:        state [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int EstFrag3_explicit_access_special(EstFrag3 * mat,int i,int j,int state) 
{
    return EstFrag3_EXPL_SPECIAL(mat,i,j,state); 
}    


/* Function:  PackAln_read_generic_EstFrag3(mat,h)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [EstFrag3 *]
 * Arg:          h [UNKN ] Undocumented argument [EstFrag3_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * PackAln_read_generic_EstFrag3(EstFrag3 * mat,EstFrag3_access_func_holder h) 
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


    out->score =  find_end_EstFrag3(mat,&i,&j,&state,&isspecial,h);  


    /* Add final end transition (at the moment we have not got the score! */ 
    if( (pau= PackAlnUnit_alloc()) == NULL  || add_PackAln(out,pau) == FALSE )   {  
      warn("Failed the first PackAlnUnit alloc, %d length of Alignment in EstFrag3_basic_read, returning a mess.(Sorry!)",out->len); 
      return out;    
      }  


    /* Put in positions for end trans. Remember that coordinates in C style */ 
    pau->i = i;  
    pau->j = j;  
    if( isspecial != TRUE)   
      pau->state = state;    
    else pau->state = state + 3;     
    prev=pau;    
    while( state != START || isspecial != TRUE)  { /*while state != START*/ 


      if( isspecial == TRUE )    
        max_calc_special_EstFrag3(mat,i,j,state,isspecial,&i,&j,&state,&isspecial,&cellscore,h);     
      else   
        max_calc_EstFrag3(mat,i,j,state,isspecial,&i,&j,&state,&isspecial,&cellscore,h);     
      if(i == EstFrag3_READ_OFF_ERROR || j == EstFrag3_READ_OFF_ERROR || state == EstFrag3_READ_OFF_ERROR )  {  
        warn("Problem - hit bad read off system, exiting now");  
        break;   
        }  
      if( (pau= PackAlnUnit_alloc()) == NULL  || add_PackAln(out,pau) == FALSE ) {  
        warn("Failed a PackAlnUnit alloc, %d length of Alignment in EstFrag3_basic_read, returning partial alignment",out->len); 
        break;   
        }  


      /* Put in positions for block. Remember that coordinates in C style */ 
      pau->i = i;    
      pau->j = j;    
      if( isspecial != TRUE)     
        pau->state = state;  
      else pau->state = state + 3;   
      prev->score = cellscore;   
      prev = pau;    
      } /* end of while state != START */ 


    invert_PackAln(out); 
    return out;  
}    


/* Function:  find_end_EstFrag3(mat,ri,rj,state,isspecial,h)
 *
 * Descrip: No Description
 *
 * Arg:              mat [UNKN ] Undocumented argument [EstFrag3 *]
 * Arg:               ri [UNKN ] Undocumented argument [int *]
 * Arg:               rj [UNKN ] Undocumented argument [int *]
 * Arg:            state [UNKN ] Undocumented argument [int *]
 * Arg:        isspecial [UNKN ] Undocumented argument [boolean *]
 * Arg:                h [UNKN ] Undocumented argument [EstFrag3_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int find_end_EstFrag3(EstFrag3 * mat,int * ri,int * rj,int * state,boolean * isspecial,EstFrag3_access_func_holder h) 
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


/* Function:  EstFrag3_debug_show_matrix(mat,starti,stopi,startj,stopj,ofp)
 *
 * Descrip: No Description
 *
 * Arg:           mat [UNKN ] Undocumented argument [EstFrag3 *]
 * Arg:        starti [UNKN ] Undocumented argument [int]
 * Arg:         stopi [UNKN ] Undocumented argument [int]
 * Arg:        startj [UNKN ] Undocumented argument [int]
 * Arg:         stopj [UNKN ] Undocumented argument [int]
 * Arg:           ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
void EstFrag3_debug_show_matrix(EstFrag3 * mat,int starti,int stopi,int startj,int stopj,FILE * ofp) 
{
    register int i;  
    register int j;  


    for(i=starti;i<stopi && i < mat->query->len;i++) {  
      for(j=startj;j<stopj && j < mat->target->seq->len;j++) {  
        fprintf(ofp,"Cell [%d - %d]\n",i,j);     
        fprintf(ofp,"State MATCH %d\n",EstFrag3_EXPL_MATRIX(mat,i,j,MATCH)); 
        fprintf(ofp,"State INSERT %d\n",EstFrag3_EXPL_MATRIX(mat,i,j,INSERT));   
        fprintf(ofp,"State DELETE %d\n",EstFrag3_EXPL_MATRIX(mat,i,j,DELETE));   
        fprintf(ofp,"\n\n"); 
        }  
      }  


}    


/* Function:  max_calc_EstFrag3(mat,i,j,state,isspecial,reti,retj,retstate,retspecial,cellscore,h)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [EstFrag3 *]
 * Arg:                 i [UNKN ] Undocumented argument [int]
 * Arg:                 j [UNKN ] Undocumented argument [int]
 * Arg:             state [UNKN ] Undocumented argument [int]
 * Arg:         isspecial [UNKN ] Undocumented argument [boolean]
 * Arg:              reti [UNKN ] Undocumented argument [int *]
 * Arg:              retj [UNKN ] Undocumented argument [int *]
 * Arg:          retstate [UNKN ] Undocumented argument [int *]
 * Arg:        retspecial [UNKN ] Undocumented argument [boolean *]
 * Arg:         cellscore [UNKN ] Undocumented argument [int *]
 * Arg:                 h [UNKN ] Undocumented argument [EstFrag3_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int max_calc_EstFrag3(EstFrag3 * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,EstFrag3_access_func_holder h) 
{
    register int temp;   
    register int cscore; 


    *reti = (*retj) = (*retstate) = EstFrag3_READ_OFF_ERROR; 


    if( i < 0 || j < 0 || i > mat->query->len || j > mat->target->seq->len)  {  
      warn("In EstFrag3 matrix special read off - out of bounds on matrix [i,j is %d,%d state %d in standard matrix]",i,j,state);    
      return -1;     
      }  


    /* Then you have to select the correct switch statement to figure out the readoff      */ 
    /* Somewhat odd - reverse the order of calculation and return as soon as it is correct */ 
    cscore = (*h.access_main)(mat,i,j,state);    
    switch(state)    { /*Switch state */ 
      case MATCH :   
        temp = cscore - (mat->cp->trans[PCD_INSERT_2_BASE]) -  (0);  
        if( temp == (*h.access_main)(mat,i - 1,j - 5,DELETE) )   {  
          *reti = i - 1; 
          *retj = j - 5; 
          *retstate = DELETE;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-5,DELETE);  
            }  
          return (*h.access_main)(mat,i - 1,j - 5,DELETE);   
          }  
        temp = cscore - (mat->cp->trans[PCD_INSERT_1_BASE]) -  (0);  
        if( temp == (*h.access_main)(mat,i - 1,j - 4,DELETE) )   {  
          *reti = i - 1; 
          *retj = j - 4; 
          *retstate = DELETE;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-4,DELETE);  
            }  
          return (*h.access_main)(mat,i - 1,j - 4,DELETE);   
          }  
        temp = cscore - (mat->cp->trans[PCD_DELETE_2_BASE]) -  (0);  
        if( temp == (*h.access_main)(mat,i - 1,j - 1,DELETE) )   {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = DELETE;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-1,DELETE);  
            }  
          return (*h.access_main)(mat,i - 1,j - 1,DELETE);   
          }  
        temp = cscore - (mat->cp->trans[PCD_DELETE_1_BASE]) -  (0);  
        if( temp == (*h.access_main)(mat,i - 1,j - 2,DELETE) )   {  
          *reti = i - 1; 
          *retj = j - 2; 
          *retstate = DELETE;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-2,DELETE);  
            }  
          return (*h.access_main)(mat,i - 1,j - 2,DELETE);   
          }  
        temp = cscore - (mat->cp->trans[PCD_INSERT_2_BASE]) -  (0);  
        if( temp == (*h.access_main)(mat,i - 1,j - 5,INSERT) )   {  
          *reti = i - 1; 
          *retj = j - 5; 
          *retstate = INSERT;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-5,INSERT);  
            }  
          return (*h.access_main)(mat,i - 1,j - 5,INSERT);   
          }  
        temp = cscore - (mat->cp->trans[PCD_INSERT_1_BASE]) -  (0);  
        if( temp == (*h.access_main)(mat,i - 1,j - 4,INSERT) )   {  
          *reti = i - 1; 
          *retj = j - 4; 
          *retstate = INSERT;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-4,INSERT);  
            }  
          return (*h.access_main)(mat,i - 1,j - 4,INSERT);   
          }  
        temp = cscore - (mat->cp->trans[PCD_DELETE_2_BASE]) -  (0);  
        if( temp == (*h.access_main)(mat,i - 1,j - 1,INSERT) )   {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = INSERT;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-1,INSERT);  
            }  
          return (*h.access_main)(mat,i - 1,j - 1,INSERT);   
          }  
        temp = cscore - (mat->cp->trans[PCD_DELETE_1_BASE]) -  (0);  
        if( temp == (*h.access_main)(mat,i - 1,j - 2,INSERT) )   {  
          *reti = i - 1; 
          *retj = j - 2; 
          *retstate = INSERT;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-2,INSERT);  
            }  
          return (*h.access_main)(mat,i - 1,j - 2,INSERT);   
          }  
        temp = cscore - (mat->cp->trans[PCD_INSERT_2_BASE]) -  (0);  
        if( temp == (*h.access_main)(mat,i - 1,j - 5,MATCH) )    {  
          *reti = i - 1; 
          *retj = j - 5; 
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-5,MATCH);   
            }  
          return (*h.access_main)(mat,i - 1,j - 5,MATCH);    
          }  
        temp = cscore - (mat->cp->trans[PCD_INSERT_1_BASE]) -  (0);  
        if( temp == (*h.access_main)(mat,i - 1,j - 4,MATCH) )    {  
          *reti = i - 1; 
          *retj = j - 4; 
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-4,MATCH);   
            }  
          return (*h.access_main)(mat,i - 1,j - 4,MATCH);    
          }  
        temp = cscore - (mat->cp->trans[PCD_DELETE_2_BASE]) -  (0);  
        if( temp == (*h.access_main)(mat,i - 1,j - 1,MATCH) )    {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-1,MATCH);   
            }  
          return (*h.access_main)(mat,i - 1,j - 1,MATCH);    
          }  
        temp = cscore - (mat->cp->trans[PCD_DELETE_1_BASE]) -  (0);  
        if( temp == (*h.access_main)(mat,i - 1,j - 2,MATCH) )    {  
          *reti = i - 1; 
          *retj = j - 2; 
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-2,MATCH);   
            }  
          return (*h.access_main)(mat,i - 1,j - 2,MATCH);    
          }  
        /* Has restricted position */ 
        if( (j-3) == 0  )    {  
          temp = cscore - ((mat->start_frag+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)])) -  (0);  
          if( temp == (*h.access_special)(mat,i - 1,j - 3,START) )   {  
            *reti = i - 1;   
            *retj = j - 3;   
            *retstate = START;   
            *retspecial = TRUE;  
            if( cellscore != NULL)   {  
              *cellscore = cscore - (*h.access_special)(mat,i-1,j-3,START);  
              }  
            return (*h.access_main)(mat,i - 1,j - 3,START);  
            }  
          }  
        temp = cscore - ((mat->query->seg[i]->transition[GW_START2MATCH]+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)])) -  (0); 
        if( temp == (*h.access_special)(mat,i - 1,j - 3,START) ) {  
          *reti = i - 1; 
          *retj = j - 3; 
          *retstate = START; 
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_special)(mat,i-1,j-3,START);    
            }  
          return (*h.access_main)(mat,i - 1,j - 3,START);    
          }  
        temp = cscore - ((mat->query->seg[i]->transition[GW_DELETE2MATCH]+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)])) -  (0);    
        if( temp == (*h.access_main)(mat,i - 1,j - 3,DELETE) )   {  
          *reti = i - 1; 
          *retj = j - 3; 
          *retstate = DELETE;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-3,DELETE);  
            }  
          return (*h.access_main)(mat,i - 1,j - 3,DELETE);   
          }  
        temp = cscore - ((mat->query->seg[i]->transition[GW_INSERT2MATCH]+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)])) -  (0);    
        if( temp == (*h.access_main)(mat,i - 1,j - 3,INSERT) )   {  
          *reti = i - 1; 
          *retj = j - 3; 
          *retstate = INSERT;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-3,INSERT);  
            }  
          return (*h.access_main)(mat,i - 1,j - 3,INSERT);   
          }  
        temp = cscore - ((mat->query->seg[i]->transition[GW_MATCH2MATCH]+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)])) -  (0); 
        if( temp == (*h.access_main)(mat,i - 1,j - 3,MATCH) )    {  
          *reti = i - 1; 
          *retj = j - 3; 
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-3,MATCH);   
            }  
          return (*h.access_main)(mat,i - 1,j - 3,MATCH);    
          }  
        warn("Major problem (!) - in EstFrag3 read off, position %d,%d state %d no source found!",i,j,state);    
        return (-1); 
      case INSERT :  
        temp = cscore - (mat->cp->trans[PCD_INSERT_2_BASE]) -  (0);  
        if( temp == (*h.access_main)(mat,i - 0,j - 5,DELETE) )   {  
          *reti = i - 0; 
          *retj = j - 5; 
          *retstate = DELETE;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-5,DELETE);  
            }  
          return (*h.access_main)(mat,i - 0,j - 5,DELETE);   
          }  
        temp = cscore - (mat->cp->trans[PCD_INSERT_1_BASE]) -  (0);  
        if( temp == (*h.access_main)(mat,i - 0,j - 4,DELETE) )   {  
          *reti = i - 0; 
          *retj = j - 4; 
          *retstate = DELETE;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-4,DELETE);  
            }  
          return (*h.access_main)(mat,i - 0,j - 4,DELETE);   
          }  
        temp = cscore - (mat->cp->trans[PCD_DELETE_2_BASE]) -  (0);  
        if( temp == (*h.access_main)(mat,i - 0,j - 1,DELETE) )   {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = DELETE;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,DELETE);  
            }  
          return (*h.access_main)(mat,i - 0,j - 1,DELETE);   
          }  
        temp = cscore - (mat->cp->trans[PCD_DELETE_1_BASE]) -  (0);  
        if( temp == (*h.access_main)(mat,i - 0,j - 2,DELETE) )   {  
          *reti = i - 0; 
          *retj = j - 2; 
          *retstate = DELETE;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-2,DELETE);  
            }  
          return (*h.access_main)(mat,i - 0,j - 2,DELETE);   
          }  
        temp = cscore - (mat->cp->trans[PCD_INSERT_2_BASE]) -  (0);  
        if( temp == (*h.access_main)(mat,i - 0,j - 5,MATCH) )    {  
          *reti = i - 0; 
          *retj = j - 5; 
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-5,MATCH);   
            }  
          return (*h.access_main)(mat,i - 0,j - 5,MATCH);    
          }  
        temp = cscore - (mat->cp->trans[PCD_INSERT_1_BASE]) -  (0);  
        if( temp == (*h.access_main)(mat,i - 0,j - 4,MATCH) )    {  
          *reti = i - 0; 
          *retj = j - 4; 
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-4,MATCH);   
            }  
          return (*h.access_main)(mat,i - 0,j - 4,MATCH);    
          }  
        temp = cscore - (mat->cp->trans[PCD_DELETE_2_BASE]) -  (0);  
        if( temp == (*h.access_main)(mat,i - 0,j - 1,MATCH) )    {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,MATCH);   
            }  
          return (*h.access_main)(mat,i - 0,j - 1,MATCH);    
          }  
        temp = cscore - (mat->cp->trans[PCD_DELETE_1_BASE]) -  (0);  
        if( temp == (*h.access_main)(mat,i - 0,j - 2,MATCH) )    {  
          *reti = i - 0; 
          *retj = j - 2; 
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-2,MATCH);   
            }  
          return (*h.access_main)(mat,i - 0,j - 2,MATCH);    
          }  
        temp = cscore - (mat->cp->trans[PCD_INSERT_2_BASE]) -  (0);  
        if( temp == (*h.access_main)(mat,i - 0,j - 5,INSERT) )   {  
          *reti = i - 0; 
          *retj = j - 5; 
          *retstate = INSERT;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-5,INSERT);  
            }  
          return (*h.access_main)(mat,i - 0,j - 5,INSERT);   
          }  
        temp = cscore - (mat->cp->trans[PCD_INSERT_1_BASE]) -  (0);  
        if( temp == (*h.access_main)(mat,i - 0,j - 4,INSERT) )   {  
          *reti = i - 0; 
          *retj = j - 4; 
          *retstate = INSERT;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-4,INSERT);  
            }  
          return (*h.access_main)(mat,i - 0,j - 4,INSERT);   
          }  
        temp = cscore - (mat->cp->trans[PCD_DELETE_2_BASE]) -  (0);  
        if( temp == (*h.access_main)(mat,i - 0,j - 1,INSERT) )   {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = INSERT;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,INSERT);  
            }  
          return (*h.access_main)(mat,i - 0,j - 1,INSERT);   
          }  
        temp = cscore - (mat->cp->trans[PCD_DELETE_1_BASE]) -  (0);  
        if( temp == (*h.access_main)(mat,i - 0,j - 2,INSERT) )   {  
          *reti = i - 0; 
          *retj = j - 2; 
          *retstate = INSERT;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-2,INSERT);  
            }  
          return (*h.access_main)(mat,i - 0,j - 2,INSERT);   
          }  
        temp = cscore - ((mat->query->seg[i]->transition[GW_START2INSERT]+mat->query->seg[i]->insert[CSEQ_CDNA_CODON(mat->target,j)])) -  (0);   
        if( temp == (*h.access_special)(mat,i - 0,j - 3,START) ) {  
          *reti = i - 0; 
          *retj = j - 3; 
          *retstate = START; 
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_special)(mat,i-0,j-3,START);    
            }  
          return (*h.access_main)(mat,i - 0,j - 3,START);    
          }  
        temp = cscore - ((mat->query->seg[i]->transition[GW_DELETE2INSERT]+mat->query->seg[i]->insert[CSEQ_CDNA_CODON(mat->target,j)])) -  (0);  
        if( temp == (*h.access_main)(mat,i - 0,j - 3,DELETE) )   {  
          *reti = i - 0; 
          *retj = j - 3; 
          *retstate = DELETE;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-3,DELETE);  
            }  
          return (*h.access_main)(mat,i - 0,j - 3,DELETE);   
          }  
        temp = cscore - ((mat->query->seg[i]->transition[GW_INSERT2INSERT]+mat->query->seg[i]->insert[CSEQ_CDNA_CODON(mat->target,j)])) -  (0);  
        if( temp == (*h.access_main)(mat,i - 0,j - 3,INSERT) )   {  
          *reti = i - 0; 
          *retj = j - 3; 
          *retstate = INSERT;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-3,INSERT);  
            }  
          return (*h.access_main)(mat,i - 0,j - 3,INSERT);   
          }  
        temp = cscore - ((mat->query->seg[i]->transition[GW_MATCH2INSERT]+mat->query->seg[i]->insert[CSEQ_CDNA_CODON(mat->target,j)])) -  (0);   
        if( temp == (*h.access_main)(mat,i - 0,j - 3,MATCH) )    {  
          *reti = i - 0; 
          *retj = j - 3; 
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-3,MATCH);   
            }  
          return (*h.access_main)(mat,i - 0,j - 3,MATCH);    
          }  
        warn("Major problem (!) - in EstFrag3 read off, position %d,%d state %d no source found!",i,j,state);    
        return (-1); 
      case DELETE :  
        temp = cscore - (mat->query->seg[i]->transition[GW_START2DELETE]) -  (0);    
        if( temp == (*h.access_special)(mat,i - 1,j - 0,START) ) {  
          *reti = i - 1; 
          *retj = j - 0; 
          *retstate = START; 
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_special)(mat,i-1,j-0,START);    
            }  
          return (*h.access_main)(mat,i - 1,j - 0,START);    
          }  
        temp = cscore - (mat->query->seg[i]->transition[GW_DELETE2DELETE]) -  (0);   
        if( temp == (*h.access_main)(mat,i - 1,j - 0,DELETE) )   {  
          *reti = i - 1; 
          *retj = j - 0; 
          *retstate = DELETE;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-0,DELETE);  
            }  
          return (*h.access_main)(mat,i - 1,j - 0,DELETE);   
          }  
        temp = cscore - (mat->query->seg[i]->transition[GW_INSERT2DELETE]) -  (0);   
        if( temp == (*h.access_main)(mat,i - 1,j - 0,INSERT) )   {  
          *reti = i - 1; 
          *retj = j - 0; 
          *retstate = INSERT;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-0,INSERT);  
            }  
          return (*h.access_main)(mat,i - 1,j - 0,INSERT);   
          }  
        temp = cscore - (mat->query->seg[i]->transition[GW_MATCH2DELETE]) -  (0);    
        if( temp == (*h.access_main)(mat,i - 1,j - 0,MATCH) )    {  
          *reti = i - 1; 
          *retj = j - 0; 
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-0,MATCH);   
            }  
          return (*h.access_main)(mat,i - 1,j - 0,MATCH);    
          }  
        warn("Major problem (!) - in EstFrag3 read off, position %d,%d state %d no source found!",i,j,state);    
        return (-1); 
      default:   
        warn("Major problem (!) - in EstFrag3 read off, position %d,%d state %d no source found!",i,j,state);    
        return (-1); 
      } /* end of Switch state  */ 
}    


/* Function:  max_calc_special_EstFrag3(mat,i,j,state,isspecial,reti,retj,retstate,retspecial,cellscore,h)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [EstFrag3 *]
 * Arg:                 i [UNKN ] Undocumented argument [int]
 * Arg:                 j [UNKN ] Undocumented argument [int]
 * Arg:             state [UNKN ] Undocumented argument [int]
 * Arg:         isspecial [UNKN ] Undocumented argument [boolean]
 * Arg:              reti [UNKN ] Undocumented argument [int *]
 * Arg:              retj [UNKN ] Undocumented argument [int *]
 * Arg:          retstate [UNKN ] Undocumented argument [int *]
 * Arg:        retspecial [UNKN ] Undocumented argument [boolean *]
 * Arg:         cellscore [UNKN ] Undocumented argument [int *]
 * Arg:                 h [UNKN ] Undocumented argument [EstFrag3_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int max_calc_special_EstFrag3(EstFrag3 * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,EstFrag3_access_func_holder h) 
{
    register int temp;   
    register int cscore; 


    *reti = (*retj) = (*retstate) = EstFrag3_READ_OFF_ERROR; 


    if( j < 0 || j > mat->target->seq->len)  {  
      warn("In EstFrag3 matrix special read off - out of bounds on matrix [j is %d in special]",j);  
      return -1;     
      }  


    cscore = (*h.access_special)(mat,i,j,state); 
    switch(state)    { /*switch on special states*/ 
      case START :   
      case END :     
        /* source DELETE is from main matrix */ 
        for(i= mat->query->len-1;i >= 0 ;i--)    { /*for i >= 0*/ 
          temp = cscore - (mat->query->seg[i]->transition[GW_DELETE2END]) - (0);     
          if( temp == (*h.access_main)(mat,i - 0,j - 0,DELETE) ) {  
            *reti = i - 0;   
            *retj = j - 0;   
            *retstate = DELETE;  
            *retspecial = FALSE; 
            if( cellscore != NULL)   {  
              *cellscore = cscore - (*h.access_main)(mat,i-0,j-0,DELETE);    
              }  
            return (*h.access_main)(mat,i - 0,j - 0,DELETE) ;    
            }  
          } /* end of for i >= 0 */ 
        /* source INSERT is from main matrix */ 
        for(i= mat->query->len-1;i >= 0 ;i--)    { /*for i >= 0*/ 
          temp = cscore - (mat->query->seg[i]->transition[GW_INSERT2END]) - (0);     
          if( temp == (*h.access_main)(mat,i - 0,j - 0,INSERT) ) {  
            *reti = i - 0;   
            *retj = j - 0;   
            *retstate = INSERT;  
            *retspecial = FALSE; 
            if( cellscore != NULL)   {  
              *cellscore = cscore - (*h.access_main)(mat,i-0,j-0,INSERT);    
              }  
            return (*h.access_main)(mat,i - 0,j - 0,INSERT) ;    
            }  
          } /* end of for i >= 0 */ 
        /* source MATCH is from main matrix */ 
        for(i= mat->query->len-1;i >= 0 ;i--)    { /*for i >= 0*/ 
          temp = cscore - (mat->end_frag) - (0);     
          if( temp == (*h.access_main)(mat,i - 0,j - 0,MATCH) )  {  
            *reti = i - 0;   
            *retj = j - 0;   
            *retstate = MATCH;   
            *retspecial = FALSE; 
            if( cellscore != NULL)   {  
              *cellscore = cscore - (*h.access_main)(mat,i-0,j-0,MATCH);     
              }  
            return (*h.access_main)(mat,i - 0,j - 0,MATCH) ;     
            }  
          } /* end of for i >= 0 */ 
        /* source MATCH is from main matrix */ 
        for(i= mat->query->len-1;i >= 0 ;i--)    { /*for i >= 0*/ 
          temp = cscore - (mat->query->seg[i]->transition[GW_MATCH2END]) - (0);  
          if( temp == (*h.access_main)(mat,i - 0,j - 0,MATCH) )  {  
            *reti = i - 0;   
            *retj = j - 0;   
            *retstate = MATCH;   
            *retspecial = FALSE; 
            if( cellscore != NULL)   {  
              *cellscore = cscore - (*h.access_main)(mat,i-0,j-0,MATCH);     
              }  
            return (*h.access_main)(mat,i - 0,j - 0,MATCH) ;     
            }  
          } /* end of for i >= 0 */ 
      default:   
        warn("Major problem (!) - in EstFrag3 read off, position %d,%d state %d no source found  dropped into default on source switch!",i,j,state); 
        return (-1); 
      } /* end of switch on special states */ 
}    


/* Function:  calculate_EstFrag3(mat)
 *
 * Descrip:    This function calculates the EstFrag3 matrix when in explicit mode
 *             To allocate the matrix use /allocate_Expl_EstFrag3
 *
 *
 * Arg:        mat [UNKN ] EstFrag3 which contains explicit basematrix memory [EstFrag3 *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean calculate_EstFrag3(EstFrag3 * mat) 
{
    int i;   
    int j;   
    int leni;    
    int lenj;    
    int tot; 
    int num; 


    if( mat->basematrix->type != BASEMATRIX_TYPE_EXPLICIT )  {  
      warn("in calculate_EstFrag3, passed a non Explicit matrix type, cannot calculate!");   
      return FALSE;  
      }  


    leni = mat->leni;    
    lenj = mat->lenj;    
    tot = leni * lenj;   
    num = 0; 


    start_reporting("EstFrag3 Matrix calculation: ");    
    for(j=0;j<lenj;j++)  {  
      auto int score;    
      auto int temp;     
      for(i=0;i<leni;i++)    {  
        if( num%1000 == 0 )  
          log_full_error(REPORT,0,"[%7d] Cells %2d%%%%",num,num*100/tot);    
        num++;   


        /* For state MATCH */ 
        /* setting first movement to score */ 
        score = EstFrag3_EXPL_MATRIX(mat,i-1,j-3,MATCH) + (mat->query->seg[i]->transition[GW_MATCH2MATCH]+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)]);    
        /* From state INSERT to state MATCH */ 
        temp = EstFrag3_EXPL_MATRIX(mat,i-1,j-3,INSERT) + (mat->query->seg[i]->transition[GW_INSERT2MATCH]+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)]);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state MATCH */ 
        temp = EstFrag3_EXPL_MATRIX(mat,i-1,j-3,DELETE) + (mat->query->seg[i]->transition[GW_DELETE2MATCH]+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)]);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state START to state MATCH */ 
        temp = EstFrag3_EXPL_SPECIAL(mat,i-1,j-3,START) + (mat->query->seg[i]->transition[GW_START2MATCH]+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)]);    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* Has restricted position */ 
        if( (j-3) == 0  )    {  
          /* From state START to state MATCH */ 
          temp = EstFrag3_EXPL_SPECIAL(mat,i-1,j-3,START) + (mat->start_frag+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)]);     
          if( temp  > score )    {  
            score = temp;    
            }  
          }  
        /* From state MATCH to state MATCH */ 
        temp = EstFrag3_EXPL_MATRIX(mat,i-1,j-2,MATCH) + mat->cp->trans[PCD_DELETE_1_BASE];  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = EstFrag3_EXPL_MATRIX(mat,i-1,j-1,MATCH) + mat->cp->trans[PCD_DELETE_2_BASE];  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = EstFrag3_EXPL_MATRIX(mat,i-1,j-4,MATCH) + mat->cp->trans[PCD_INSERT_1_BASE];  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = EstFrag3_EXPL_MATRIX(mat,i-1,j-5,MATCH) + mat->cp->trans[PCD_INSERT_2_BASE];  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state MATCH */ 
        temp = EstFrag3_EXPL_MATRIX(mat,i-1,j-2,INSERT) + mat->cp->trans[PCD_DELETE_1_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state MATCH */ 
        temp = EstFrag3_EXPL_MATRIX(mat,i-1,j-1,INSERT) + mat->cp->trans[PCD_DELETE_2_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state MATCH */ 
        temp = EstFrag3_EXPL_MATRIX(mat,i-1,j-4,INSERT) + mat->cp->trans[PCD_INSERT_1_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state MATCH */ 
        temp = EstFrag3_EXPL_MATRIX(mat,i-1,j-5,INSERT) + mat->cp->trans[PCD_INSERT_2_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state MATCH */ 
        temp = EstFrag3_EXPL_MATRIX(mat,i-1,j-2,DELETE) + mat->cp->trans[PCD_DELETE_1_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state MATCH */ 
        temp = EstFrag3_EXPL_MATRIX(mat,i-1,j-1,DELETE) + mat->cp->trans[PCD_DELETE_2_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state MATCH */ 
        temp = EstFrag3_EXPL_MATRIX(mat,i-1,j-4,DELETE) + mat->cp->trans[PCD_INSERT_1_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state MATCH */ 
        temp = EstFrag3_EXPL_MATRIX(mat,i-1,j-5,DELETE) + mat->cp->trans[PCD_INSERT_2_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for MATCH */ 
        /* Add any movement independant score and put away */ 
         EstFrag3_EXPL_MATRIX(mat,i,j,MATCH) = score;    


        /* state MATCH is a source for special END */ 
        temp = score + (mat->query->seg[i]->transition[GW_MATCH2END]) + (0) ;    
        if( temp > EstFrag3_EXPL_SPECIAL(mat,i,j,END) )  {  
          EstFrag3_EXPL_SPECIAL(mat,i,j,END) = temp;     
          }  




        /* state MATCH is a source for special END */ 
        /* Has restricted position */ 
        if( j == mat->lenj-1 )   {  
          temp = score + (mat->end_frag) + (0) ;     
          if( temp > EstFrag3_EXPL_SPECIAL(mat,i,j,END) )    {  
            EstFrag3_EXPL_SPECIAL(mat,i,j,END) = temp;   
            }  


          }  


        /* Finished calculating state MATCH */ 


        /* For state INSERT */ 
        /* setting first movement to score */ 
        score = EstFrag3_EXPL_MATRIX(mat,i-0,j-3,MATCH) + (mat->query->seg[i]->transition[GW_MATCH2INSERT]+mat->query->seg[i]->insert[CSEQ_CDNA_CODON(mat->target,j)]);  
        /* From state INSERT to state INSERT */ 
        temp = EstFrag3_EXPL_MATRIX(mat,i-0,j-3,INSERT) + (mat->query->seg[i]->transition[GW_INSERT2INSERT]+mat->query->seg[i]->insert[CSEQ_CDNA_CODON(mat->target,j)]);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state INSERT */ 
        temp = EstFrag3_EXPL_MATRIX(mat,i-0,j-3,DELETE) + (mat->query->seg[i]->transition[GW_DELETE2INSERT]+mat->query->seg[i]->insert[CSEQ_CDNA_CODON(mat->target,j)]);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state START to state INSERT */ 
        temp = EstFrag3_EXPL_SPECIAL(mat,i-0,j-3,START) + (mat->query->seg[i]->transition[GW_START2INSERT]+mat->query->seg[i]->insert[CSEQ_CDNA_CODON(mat->target,j)]);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state INSERT */ 
        temp = EstFrag3_EXPL_MATRIX(mat,i-0,j-2,INSERT) + mat->cp->trans[PCD_DELETE_1_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state INSERT */ 
        temp = EstFrag3_EXPL_MATRIX(mat,i-0,j-1,INSERT) + mat->cp->trans[PCD_DELETE_2_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state INSERT */ 
        temp = EstFrag3_EXPL_MATRIX(mat,i-0,j-4,INSERT) + mat->cp->trans[PCD_INSERT_1_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state INSERT */ 
        temp = EstFrag3_EXPL_MATRIX(mat,i-0,j-5,INSERT) + mat->cp->trans[PCD_INSERT_2_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state INSERT */ 
        temp = EstFrag3_EXPL_MATRIX(mat,i-0,j-2,MATCH) + mat->cp->trans[PCD_DELETE_1_BASE];  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state INSERT */ 
        temp = EstFrag3_EXPL_MATRIX(mat,i-0,j-1,MATCH) + mat->cp->trans[PCD_DELETE_2_BASE];  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state INSERT */ 
        temp = EstFrag3_EXPL_MATRIX(mat,i-0,j-4,MATCH) + mat->cp->trans[PCD_INSERT_1_BASE];  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state INSERT */ 
        temp = EstFrag3_EXPL_MATRIX(mat,i-0,j-5,MATCH) + mat->cp->trans[PCD_INSERT_2_BASE];  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state INSERT */ 
        temp = EstFrag3_EXPL_MATRIX(mat,i-0,j-2,DELETE) + mat->cp->trans[PCD_DELETE_1_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state INSERT */ 
        temp = EstFrag3_EXPL_MATRIX(mat,i-0,j-1,DELETE) + mat->cp->trans[PCD_DELETE_2_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state INSERT */ 
        temp = EstFrag3_EXPL_MATRIX(mat,i-0,j-4,DELETE) + mat->cp->trans[PCD_INSERT_1_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state INSERT */ 
        temp = EstFrag3_EXPL_MATRIX(mat,i-0,j-5,DELETE) + mat->cp->trans[PCD_INSERT_2_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for INSERT */ 
        /* Add any movement independant score and put away */ 
         EstFrag3_EXPL_MATRIX(mat,i,j,INSERT) = score;   


        /* state INSERT is a source for special END */ 
        temp = score + (mat->query->seg[i]->transition[GW_INSERT2END]) + (0) ;   
        if( temp > EstFrag3_EXPL_SPECIAL(mat,i,j,END) )  {  
          EstFrag3_EXPL_SPECIAL(mat,i,j,END) = temp;     
          }  




        /* Finished calculating state INSERT */ 


        /* For state DELETE */ 
        /* setting first movement to score */ 
        score = EstFrag3_EXPL_MATRIX(mat,i-1,j-0,MATCH) + mat->query->seg[i]->transition[GW_MATCH2DELETE];   
        /* From state INSERT to state DELETE */ 
        temp = EstFrag3_EXPL_MATRIX(mat,i-1,j-0,INSERT) + mat->query->seg[i]->transition[GW_INSERT2DELETE];  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state DELETE */ 
        temp = EstFrag3_EXPL_MATRIX(mat,i-1,j-0,DELETE) + mat->query->seg[i]->transition[GW_DELETE2DELETE];  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state START to state DELETE */ 
        temp = EstFrag3_EXPL_SPECIAL(mat,i-1,j-0,START) + mat->query->seg[i]->transition[GW_START2DELETE];   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for DELETE */ 
        /* Add any movement independant score and put away */ 
         EstFrag3_EXPL_MATRIX(mat,i,j,DELETE) = score;   


        /* state DELETE is a source for special END */ 
        temp = score + (mat->query->seg[i]->transition[GW_DELETE2END]) + (0) ;   
        if( temp > EstFrag3_EXPL_SPECIAL(mat,i,j,END) )  {  
          EstFrag3_EXPL_SPECIAL(mat,i,j,END) = temp;     
          }  




        /* Finished calculating state DELETE */ 
        }  


      /* Special state START has no special to special movements */ 


      /* Special state END has no special to special movements */ 
      }  
    stop_reporting();    
    return TRUE;     
}    


/* Function:  calculate_dpenv_EstFrag3(mat,dpenv)
 *
 * Descrip:    This function calculates the EstFrag3 matrix when in explicit mode, subject to the envelope
 *
 *
 * Arg:          mat [UNKN ] EstFrag3 which contains explicit basematrix memory [EstFrag3 *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean calculate_dpenv_EstFrag3(EstFrag3 * mat,DPEnvelope * dpenv) 
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
      warn("in calculate_EstFrag3, passed a non Explicit matrix type, cannot calculate!");   
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


    for(j=startj-5;j<endj;j++)   {  
      for(i=1;i<mat->leni;i++)   {  
        EstFrag3_EXPL_MATRIX(mat,i,j,MATCH) = NEGI;  
        EstFrag3_EXPL_MATRIX(mat,i,j,INSERT) = NEGI; 
        EstFrag3_EXPL_MATRIX(mat,i,j,DELETE) = NEGI; 
        }  
      }  
    for(j=-5;j<mat->lenj;j++)    {  
      EstFrag3_EXPL_SPECIAL(mat,i,j,START) = 0;  
      EstFrag3_EXPL_SPECIAL(mat,i,j,END) = NEGI; 
      }  


    start_reporting("EstFrag3 Matrix calculation: ");    
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
          EstFrag3_EXPL_MATRIX(mat,i,j,MATCH) = NEGI;    
          EstFrag3_EXPL_MATRIX(mat,i,j,INSERT) = NEGI;   
          EstFrag3_EXPL_MATRIX(mat,i,j,DELETE) = NEGI;   
          continue;  
          }  


        if( num%1000 == 0 )  
          log_full_error(REPORT,0,"[%7d] Cells %2d%%%%",num,num*100/tot);    
        num++;   


        /* For state MATCH */ 
        /* setting first movement to score */ 
        score = EstFrag3_EXPL_MATRIX(mat,i-1,j-3,MATCH) + (mat->query->seg[i]->transition[GW_MATCH2MATCH]+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)]);    
        /* From state INSERT to state MATCH */ 
        temp = EstFrag3_EXPL_MATRIX(mat,i-1,j-3,INSERT) + (mat->query->seg[i]->transition[GW_INSERT2MATCH]+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)]);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state MATCH */ 
        temp = EstFrag3_EXPL_MATRIX(mat,i-1,j-3,DELETE) + (mat->query->seg[i]->transition[GW_DELETE2MATCH]+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)]);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state START to state MATCH */ 
        temp = EstFrag3_EXPL_SPECIAL(mat,i-1,j-3,START) + (mat->query->seg[i]->transition[GW_START2MATCH]+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)]);    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* Has restricted position */ 
        if( (j-3) == 0  )    {  
          /* From state START to state MATCH */ 
          temp = EstFrag3_EXPL_SPECIAL(mat,i-1,j-3,START) + (mat->start_frag+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)]);     
          if( temp  > score )    {  
            score = temp;    
            }  
          }  
        /* From state MATCH to state MATCH */ 
        temp = EstFrag3_EXPL_MATRIX(mat,i-1,j-2,MATCH) + mat->cp->trans[PCD_DELETE_1_BASE];  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = EstFrag3_EXPL_MATRIX(mat,i-1,j-1,MATCH) + mat->cp->trans[PCD_DELETE_2_BASE];  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = EstFrag3_EXPL_MATRIX(mat,i-1,j-4,MATCH) + mat->cp->trans[PCD_INSERT_1_BASE];  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = EstFrag3_EXPL_MATRIX(mat,i-1,j-5,MATCH) + mat->cp->trans[PCD_INSERT_2_BASE];  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state MATCH */ 
        temp = EstFrag3_EXPL_MATRIX(mat,i-1,j-2,INSERT) + mat->cp->trans[PCD_DELETE_1_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state MATCH */ 
        temp = EstFrag3_EXPL_MATRIX(mat,i-1,j-1,INSERT) + mat->cp->trans[PCD_DELETE_2_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state MATCH */ 
        temp = EstFrag3_EXPL_MATRIX(mat,i-1,j-4,INSERT) + mat->cp->trans[PCD_INSERT_1_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state MATCH */ 
        temp = EstFrag3_EXPL_MATRIX(mat,i-1,j-5,INSERT) + mat->cp->trans[PCD_INSERT_2_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state MATCH */ 
        temp = EstFrag3_EXPL_MATRIX(mat,i-1,j-2,DELETE) + mat->cp->trans[PCD_DELETE_1_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state MATCH */ 
        temp = EstFrag3_EXPL_MATRIX(mat,i-1,j-1,DELETE) + mat->cp->trans[PCD_DELETE_2_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state MATCH */ 
        temp = EstFrag3_EXPL_MATRIX(mat,i-1,j-4,DELETE) + mat->cp->trans[PCD_INSERT_1_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state MATCH */ 
        temp = EstFrag3_EXPL_MATRIX(mat,i-1,j-5,DELETE) + mat->cp->trans[PCD_INSERT_2_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for MATCH */ 
        /* Add any movement independant score and put away */ 
         EstFrag3_EXPL_MATRIX(mat,i,j,MATCH) = score;    


        /* state MATCH is a source for special END */ 
        temp = score + (mat->query->seg[i]->transition[GW_MATCH2END]) + (0) ;    
        if( temp > EstFrag3_EXPL_SPECIAL(mat,i,j,END) )  {  
          EstFrag3_EXPL_SPECIAL(mat,i,j,END) = temp;     
          }  




        /* state MATCH is a source for special END */ 
        /* Has restricted position */ 
        if( j == mat->lenj-1 )   {  
          temp = score + (mat->end_frag) + (0) ;     
          if( temp > EstFrag3_EXPL_SPECIAL(mat,i,j,END) )    {  
            EstFrag3_EXPL_SPECIAL(mat,i,j,END) = temp;   
            }  


          }  


        /* Finished calculating state MATCH */ 


        /* For state INSERT */ 
        /* setting first movement to score */ 
        score = EstFrag3_EXPL_MATRIX(mat,i-0,j-3,MATCH) + (mat->query->seg[i]->transition[GW_MATCH2INSERT]+mat->query->seg[i]->insert[CSEQ_CDNA_CODON(mat->target,j)]);  
        /* From state INSERT to state INSERT */ 
        temp = EstFrag3_EXPL_MATRIX(mat,i-0,j-3,INSERT) + (mat->query->seg[i]->transition[GW_INSERT2INSERT]+mat->query->seg[i]->insert[CSEQ_CDNA_CODON(mat->target,j)]);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state INSERT */ 
        temp = EstFrag3_EXPL_MATRIX(mat,i-0,j-3,DELETE) + (mat->query->seg[i]->transition[GW_DELETE2INSERT]+mat->query->seg[i]->insert[CSEQ_CDNA_CODON(mat->target,j)]);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state START to state INSERT */ 
        temp = EstFrag3_EXPL_SPECIAL(mat,i-0,j-3,START) + (mat->query->seg[i]->transition[GW_START2INSERT]+mat->query->seg[i]->insert[CSEQ_CDNA_CODON(mat->target,j)]);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state INSERT */ 
        temp = EstFrag3_EXPL_MATRIX(mat,i-0,j-2,INSERT) + mat->cp->trans[PCD_DELETE_1_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state INSERT */ 
        temp = EstFrag3_EXPL_MATRIX(mat,i-0,j-1,INSERT) + mat->cp->trans[PCD_DELETE_2_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state INSERT */ 
        temp = EstFrag3_EXPL_MATRIX(mat,i-0,j-4,INSERT) + mat->cp->trans[PCD_INSERT_1_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state INSERT */ 
        temp = EstFrag3_EXPL_MATRIX(mat,i-0,j-5,INSERT) + mat->cp->trans[PCD_INSERT_2_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state INSERT */ 
        temp = EstFrag3_EXPL_MATRIX(mat,i-0,j-2,MATCH) + mat->cp->trans[PCD_DELETE_1_BASE];  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state INSERT */ 
        temp = EstFrag3_EXPL_MATRIX(mat,i-0,j-1,MATCH) + mat->cp->trans[PCD_DELETE_2_BASE];  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state INSERT */ 
        temp = EstFrag3_EXPL_MATRIX(mat,i-0,j-4,MATCH) + mat->cp->trans[PCD_INSERT_1_BASE];  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state INSERT */ 
        temp = EstFrag3_EXPL_MATRIX(mat,i-0,j-5,MATCH) + mat->cp->trans[PCD_INSERT_2_BASE];  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state INSERT */ 
        temp = EstFrag3_EXPL_MATRIX(mat,i-0,j-2,DELETE) + mat->cp->trans[PCD_DELETE_1_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state INSERT */ 
        temp = EstFrag3_EXPL_MATRIX(mat,i-0,j-1,DELETE) + mat->cp->trans[PCD_DELETE_2_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state INSERT */ 
        temp = EstFrag3_EXPL_MATRIX(mat,i-0,j-4,DELETE) + mat->cp->trans[PCD_INSERT_1_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state INSERT */ 
        temp = EstFrag3_EXPL_MATRIX(mat,i-0,j-5,DELETE) + mat->cp->trans[PCD_INSERT_2_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for INSERT */ 
        /* Add any movement independant score and put away */ 
         EstFrag3_EXPL_MATRIX(mat,i,j,INSERT) = score;   


        /* state INSERT is a source for special END */ 
        temp = score + (mat->query->seg[i]->transition[GW_INSERT2END]) + (0) ;   
        if( temp > EstFrag3_EXPL_SPECIAL(mat,i,j,END) )  {  
          EstFrag3_EXPL_SPECIAL(mat,i,j,END) = temp;     
          }  




        /* Finished calculating state INSERT */ 


        /* For state DELETE */ 
        /* setting first movement to score */ 
        score = EstFrag3_EXPL_MATRIX(mat,i-1,j-0,MATCH) + mat->query->seg[i]->transition[GW_MATCH2DELETE];   
        /* From state INSERT to state DELETE */ 
        temp = EstFrag3_EXPL_MATRIX(mat,i-1,j-0,INSERT) + mat->query->seg[i]->transition[GW_INSERT2DELETE];  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state DELETE */ 
        temp = EstFrag3_EXPL_MATRIX(mat,i-1,j-0,DELETE) + mat->query->seg[i]->transition[GW_DELETE2DELETE];  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state START to state DELETE */ 
        temp = EstFrag3_EXPL_SPECIAL(mat,i-1,j-0,START) + mat->query->seg[i]->transition[GW_START2DELETE];   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for DELETE */ 
        /* Add any movement independant score and put away */ 
         EstFrag3_EXPL_MATRIX(mat,i,j,DELETE) = score;   


        /* state DELETE is a source for special END */ 
        temp = score + (mat->query->seg[i]->transition[GW_DELETE2END]) + (0) ;   
        if( temp > EstFrag3_EXPL_SPECIAL(mat,i,j,END) )  {  
          EstFrag3_EXPL_SPECIAL(mat,i,j,END) = temp;     
          }  




        /* Finished calculating state DELETE */ 
        }  


      /* Special state START has no special to special movements */ 


      /* Special state END has no special to special movements */ 
      }  
    stop_reporting();    
    return TRUE;     
}    


/* Function:  EstFrag3_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [EstFrag3 *]
 *
 */
EstFrag3 * EstFrag3_alloc(void) 
{
    EstFrag3 * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(EstFrag3 *) ckalloc (sizeof(EstFrag3))) == NULL)    {  
      warn("EstFrag3_alloc failed ");    
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


/* Function:  free_EstFrag3(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [EstFrag3 *]
 *
 * Return [UNKN ]  Undocumented return value [EstFrag3 *]
 *
 */
EstFrag3 * free_EstFrag3(EstFrag3 * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a EstFrag3 obj. Should be trappable");  
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
    /* obj->cp is linked in */ 
    /* obj->start_frag is linked in */ 
    /* obj->end_frag is linked in */ 


    ckfree(obj); 
    return NULL; 
}    





#ifdef _cplusplus
}
#endif
