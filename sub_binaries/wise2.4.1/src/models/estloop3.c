#ifdef _cplusplus
extern "C" {
#endif
#include "estloop3.h"

# line 5 "estloop3.c"


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


#define LOOP 0   
#define START 1  
#define END 2    


#define EstLoop3_EXPL_MATRIX(this_matrix,i,j,STATE) this_matrix->basematrix->matrix[((j+5)*3)+STATE][i+1]    
#define EstLoop3_EXPL_SPECIAL(matrix,i,j,STATE) matrix->basematrix->specmatrix[STATE][j+5]   
#define EstLoop3_READ_OFF_ERROR -7
  


#define EstLoop3_VSMALL_MATRIX(mat,i,j,STATE) mat->basematrix->matrix[(j+6)%6][((i+1)*3)+STATE]  
#define EstLoop3_VSMALL_SPECIAL(mat,i,j,STATE) mat->basematrix->specmatrix[(j+6)%6][STATE]   




#define EstLoop3_SHATTER_SPECIAL(matrix,i,j,STATE) matrix->shatter->special[STATE][j]    
#define EstLoop3_SHATTER_MATRIX(matrix,i,j,STATE)  fetch_cell_value_ShatterMatrix(mat->shatter,i,j,STATE)    


/* Function:  PackAln_read_Shatter_EstLoop3(mat)
 *
 * Descrip:    Reads off PackAln from shatter matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [EstLoop3 *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * PackAln_read_Shatter_EstLoop3(EstLoop3 * mat) 
{
    EstLoop3_access_func_holder holder;  


    holder.access_main    = EstLoop3_shatter_access_main;    
    holder.access_special = EstLoop3_shatter_access_special; 
    assert(mat);     
    assert(mat->shatter);    
    return PackAln_read_generic_EstLoop3(mat,holder);    
}    


/* Function:  EstLoop3_shatter_access_main(mat,i,j,state)
 *
 * Descrip: No Description
 *
 * Arg:          mat [UNKN ] Undocumented argument [EstLoop3 *]
 * Arg:            i [UNKN ] Undocumented argument [int]
 * Arg:            j [UNKN ] Undocumented argument [int]
 * Arg:        state [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int EstLoop3_shatter_access_main(EstLoop3 * mat,int i,int j,int state) 
{
    return EstLoop3_SHATTER_MATRIX(mat,i,j,state);   
}    


/* Function:  EstLoop3_shatter_access_special(mat,i,j,state)
 *
 * Descrip: No Description
 *
 * Arg:          mat [UNKN ] Undocumented argument [EstLoop3 *]
 * Arg:            i [UNKN ] Undocumented argument [int]
 * Arg:            j [UNKN ] Undocumented argument [int]
 * Arg:        state [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int EstLoop3_shatter_access_special(EstLoop3 * mat,int i,int j,int state) 
{
    return EstLoop3_SHATTER_SPECIAL(mat,i,j,state);  
}    


/* Function:  calculate_shatter_EstLoop3(mat,dpenv)
 *
 * Descrip:    This function calculates the EstLoop3 matrix when in shatter mode
 *
 *
 * Arg:          mat [UNKN ] (null) [EstLoop3 *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean calculate_shatter_EstLoop3(EstLoop3 * mat,DPEnvelope * dpenv) 
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


    mat->shatter = new_ShatterMatrix(dpenv,3,lenj,3);    
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


    start_reporting("EstLoop3 Matrix calculation: ");    
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
        /* From state LOOP to state MATCH */ 
        temp = EstLoop3_SHATTER_SPECIAL(mat,i-1,j-3,LOOP) + (mat->query->seg[i]->transition[GW_START2MATCH]+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)]);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state START to state MATCH */ 
        temp = EstLoop3_SHATTER_SPECIAL(mat,i-1,j-3,START) + (mat->query->seg[i]->transition[GW_START2MATCH]+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)]);     
        if( temp  > score )  {  
          score = temp;  
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


        /* state MATCH is a source for special LOOP */ 
        temp = score + (mat->query->seg[i]->transition[GW_MATCH2END]) + (0) ;    
        if( temp > EstLoop3_SHATTER_SPECIAL(mat,i,j,LOOP) )  {  
          EstLoop3_SHATTER_SPECIAL(mat,i,j,LOOP) = temp;     
          }  




        /* state MATCH is a source for special END */ 
        temp = score + (mat->query->seg[i]->transition[GW_MATCH2END]) + (0) ;    
        if( temp > EstLoop3_SHATTER_SPECIAL(mat,i,j,END) )   {  
          EstLoop3_SHATTER_SPECIAL(mat,i,j,END) = temp;  
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
        /* From state LOOP to state INSERT */ 
        temp = EstLoop3_SHATTER_SPECIAL(mat,i-0,j-3,LOOP) + (mat->query->seg[i]->transition[GW_START2INSERT]+mat->query->seg[i]->insert[CSEQ_CDNA_CODON(mat->target,j)]);    
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


        /* state INSERT is a source for special LOOP */ 
        temp = score + (mat->query->seg[i]->transition[GW_INSERT2END]) + (0) ;   
        if( temp > EstLoop3_SHATTER_SPECIAL(mat,i,j,LOOP) )  {  
          EstLoop3_SHATTER_SPECIAL(mat,i,j,LOOP) = temp;     
          }  




        /* state INSERT is a source for special END */ 
        temp = score + (mat->query->seg[i]->transition[GW_INSERT2END]) + (0) ;   
        if( temp > EstLoop3_SHATTER_SPECIAL(mat,i,j,END) )   {  
          EstLoop3_SHATTER_SPECIAL(mat,i,j,END) = temp;  
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
        /* From state LOOP to state DELETE */ 
        temp = EstLoop3_SHATTER_SPECIAL(mat,i-1,j-0,LOOP) + mat->query->seg[i]->transition[GW_START2DELETE];     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for DELETE */ 
        /* Add any movement independant score and put away */ 
         SIG_0_0[DELETE] = score;    


        /* state DELETE is a source for special LOOP */ 
        temp = score + (mat->query->seg[i]->transition[GW_DELETE2END]) + (0) ;   
        if( temp > EstLoop3_SHATTER_SPECIAL(mat,i,j,LOOP) )  {  
          EstLoop3_SHATTER_SPECIAL(mat,i,j,LOOP) = temp;     
          }  




        /* state DELETE is a source for special END */ 
        temp = score + (mat->query->seg[i]->transition[GW_DELETE2END]) + (0) ;   
        if( temp > EstLoop3_SHATTER_SPECIAL(mat,i,j,END) )   {  
          EstLoop3_SHATTER_SPECIAL(mat,i,j,END) = temp;  
          }  




        /* Finished calculating state DELETE */ 
        }  


      /* Special state LOOP has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = EstLoop3_SHATTER_SPECIAL(mat,0,j,LOOP);    


      /* Source START is a special source for LOOP */ 
      /* Has restricted position */ 
      if( (j-1) == 0  )  {  
        temp = EstLoop3_SHATTER_SPECIAL(mat,0,j - 1,START) + (0) + (0);  
        if( temp > score )   
          score = temp;  
        }  


      /* Source MATCH for state LOOP is not special... already calculated */ 
      /* Source INSERT for state LOOP is not special... already calculated */ 
      /* Source DELETE for state LOOP is not special... already calculated */ 
      /* Source LOOP is a special source for LOOP */ 
      temp = EstLoop3_SHATTER_SPECIAL(mat,0,j - 1,LOOP) + (0) + (0);     
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      EstLoop3_SHATTER_SPECIAL(mat,0,j,LOOP) = score;    
      /* Finished updating state LOOP */ 




      /* Special state START has no special to special movements */ 


      /* Special state END has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = EstLoop3_SHATTER_SPECIAL(mat,0,j,END); 


      /* Source LOOP is a special source for END */ 
      /* Has restricted position */ 
      if( j == mat->lenj-1 ) {  
        temp = EstLoop3_SHATTER_SPECIAL(mat,0,j - 1,LOOP) + (0) + (0);   
        if( temp > score )   
          score = temp;  
        }  


      /* Source MATCH for state END is not special... already calculated */ 
      /* Source INSERT for state END is not special... already calculated */ 
      /* Source DELETE for state END is not special... already calculated */ 
      /* Put back score... (now updated!) */ 
      EstLoop3_SHATTER_SPECIAL(mat,0,j,END) = score; 
      /* Finished updating state END */ 


      }  
    stop_reporting();    
    return TRUE;     
}    


/* Function:  search_EstLoop3(dbsi,out,querydb,targetdb,cp)
 *
 * Descrip:    This function makes a database search of EstLoop3
 *             It uses the dbsi structure to choose which implementation
 *             to use of the database searching. This way at run time you
 *             can switch between single threaded/multi-threaded or hardware
 *
 *
 * Arg:            dbsi [UNKN ] Undocumented argument [DBSearchImpl *]
 * Arg:             out [UNKN ] Undocumented argument [Hscore *]
 * Arg:         querydb [UNKN ] Undocumented argument [GeneWiseDB*]
 * Arg:        targetdb [UNKN ] Undocumented argument [cDNADB*]
 * Arg:              cp [UNKN ] Undocumented argument [cDNAParserScore *]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type search_EstLoop3(DBSearchImpl * dbsi,Hscore * out,GeneWiseDB* querydb,cDNADB* targetdb ,cDNAParserScore * cp) 
{
#ifdef PTHREAD   
    int i;   
    int thr_no;  
    pthread_attr_t pat;  
    struct thread_pool_holder_EstLoop3 * holder;     
#endif   
    if( out == NULL )    {  
      warn("Passed in a null Hscore object into search_EstLoop3. Can't process results!");   
      return SEARCH_ERROR;   
      }  
    if( dbsi == NULL )   {  
      warn("Passed in a null DBSearchImpl object into search_EstLoop3. Can't process results!"); 
      return SEARCH_ERROR;   
      }  
    if( dbsi->trace_level > 5 )  
      warn("Asking for trace level of %d in database search for EstLoop3, but it was compiled with a trace level of -2139062144. Not all trace statements can be shown",dbsi->trace_level);  
    switch(dbsi->type)   { /*switch on implementation*/ 
      case DBSearchImpl_Serial : 
        return serial_search_EstLoop3(out,querydb, targetdb ,cp);    
      case DBSearchImpl_Pthreads :   
#ifdef PTHREAD   
        holder = (struct thread_pool_holder_EstLoop3 *) ckalloc(sizeof(struct thread_pool_holder_EstLoop3)); 
        if( holder == NULL )     {  
          warn("Unable to allocated thread pool datastructure...");  
          return SEARCH_ERROR;   
          }  
        holder->out = out;   
        holder->dbsi = dbsi; 
        holder->querydb = querydb;   
        holder->targetdb = targetdb; 
        holder->cp = cp; 
        if( pthread_mutex_init(&(holder->input_lock),NULL) != 0 )    
        fatal("Unable to iniated input mutex lock"); 
        if( pthread_mutex_init(&(holder->output_lock),NULL) != 0 )   
        fatal("Unable to iniated output mutex lock");    
        /* Let us rock! */ 
        thr_no = number_of_threads_DBSearchImpl(dbsi);   
        holder->pool = ckcalloc (thr_no,sizeof(pthread_t));  
        if( holder->pool == NULL )   {  
          warn("Unable to allocated thread pools");  
          return SEARCH_ERROR;   
          }  
        /* Build a thread attribute to make sure we get the most out of SMP boxes */ 
        pthread_attr_init(&pat);     
        /* Give thread libraries a hint that threads should be kernel threads */ 
#ifndef __sgi /* SGI can't set system scope ... */   
#ifdef  HAS_PTHREAD_SETSCOPE 
        pthread_attr_setscope(&pat, PTHREAD_SCOPE_SYSTEM);   
#endif /* set scope */   
#endif /* sgi */ 
        /* Give thread libraries a hint that there are num of threads to run */ 
#ifdef HAS_PTHREAD_SETCONCURRENCY    
        pthread_setconcurrency(thr_no+1);    
#endif /* set concurrency */ 
        for(i=0;i<thr_no;i++)    {  
          if( pthread_create(holder->pool+i,&pat,thread_loop_EstLoop3,(void *)holder) )  
            fatal("Unable to create a thread!"); 
          }  
        /* Now - wait for all the threads to exit */ 
        for(i=0;i<thr_no;i++)    {  
          if( pthread_join(holder->pool[i],NULL) != 0 )  
            fatal("Unable to join a thread!");   
          }  
        /* Deallocate the thread structures */ 
        ckfree(holder->pool);    
        ckfree(holder);  
        return SEARCH_OK;    
#else /* not compiled with threads */    
        warn("You did not specifiy the PTHREAD compile when compiled the C code for EstLoop3");  
#endif /* finished threads */    
      default :  
        warn("database search implementation %s was not provided in the compiled dynamite file from EstLoop3",impl_string_DBSearchImpl(dbsi));   
        return SEARCH_ERROR; 
      } /* end of switch on implementation */ 


}    


/* Function:  thread_loop_EstLoop3(ptr)
 *
 * Descrip:    Infinite loop code foreach thread for EstLoop3
 *
 *
 * Arg:        ptr [UNKN ] Undocumented argument [void *]
 *
 * Return [UNKN ]  Undocumented return value [void *]
 *
 */
#ifdef PTHREAD   
void * thread_loop_EstLoop3(void * ptr) 
{
    struct thread_pool_holder_EstLoop3 * holder;     
    int db_status;   
    int score;   
    DataScore * ds;  
    GeneWiseScore* query;    
    ComplexSequence* target;     


    holder = (struct thread_pool_holder_EstLoop3 *) ptr; 
    if ( holder->dbsi->trace_level >= 1 )    
      fprintf(holder->dbsi->trace_file,"Entering infinite loop for thread...\n");    


    while(1) { /*Infinite loop over all models*/ 
      /* Get input lock */ 


      if ( holder->dbsi->trace_level >= 2 )  
        fprintf(holder->dbsi->trace_file,"About to get input lock for main reload\n");   
      if( pthread_mutex_lock(&(holder->input_lock))!= 0 )    
        fatal("Error on getting input lock for EstLoop3");   
      if ( holder->dbsi->trace_level >= 2 )  
        fprintf(holder->dbsi->trace_file,"Got input lock for main reload\n");    


      if( holder->search_has_ended == TRUE ) {  
        if ( holder->dbsi->trace_level >= 2 )    
          fprintf(holder->dbsi->trace_file,"Database search finished for me!...\n"); 
        if( pthread_mutex_unlock(&(holder->input_lock))!= 0 )    
          fatal("Error in releasing input lock for EstLoop3");   
        if ( holder->dbsi->trace_level >= 2 )    
          fprintf(holder->dbsi->trace_file,"Released lock and broken out of loop\n");    
        break;   
        }  


      /* Get storage space now, as we have to read in the info for the db now */ 
      if ( holder->dbsi->trace_level >= 3 )  
        fprintf(holder->dbsi->trace_file,"Getting new DataScore from storage...\n"); 
      ds = new_DataScore();  


      /* We need to get our query object */ 
      if ( holder->dbsi->trace_level >= 3 )  
        fprintf(holder->dbsi->trace_file,"Starting query database...\n");    


      if( holder->query_init == FALSE)   {  
        holder->query = init_GeneWiseDB(holder->querydb,&db_status); 
        holder->query_init = TRUE;   
        if( db_status == DB_RETURN_ERROR )   
          fatal("Unable to initalise query database in EstLoop3 search");    
        }  
      query = hard_link_GeneWiseScore(holder->query);    
      /* get query information into datascore */ 
      dataentry_add_GeneWiseDB(ds->query,query,holder->querydb); 


      if( holder->target_init == FALSE ) { /*if the db has not been init'd*/ 
        if ( holder->dbsi->trace_level >= 3 )    
          fprintf(holder->dbsi->trace_file,"Starting target database...\n"); 
        target = init_cDNADB(holder->targetdb,&db_status);   
        holder->target_init = TRUE;  
        } /* end of if the db has not been init'd */ 
      else   { /*Normal reload*/ 
         target = reload_cDNADB(NULL,holder->targetdb,&db_status);   
        } /* end of Normal reload */ 


      /* Check to see what the reload is like */ 


      if( db_status == DB_RETURN_ERROR ) {  
        fatal("In searching EstLoop3, Reload error on database target, in threads"); 
        }  


      if( db_status == DB_RETURN_END)    { /*End of target database*/ 
        /* close target database and schedule it for initalisation by next thread */ 
        close_cDNADB(NULL,holder->targetdb); 
        holder->target_init = FALSE; 
        if ( holder->dbsi->trace_level >= 2 )    
          fprintf(holder->dbsi->trace_file,"Target Database to be reloaded...\n");   


        /* free'ing the query object */ 
        free_GeneWiseScore(holder->query);   
        /* get the next query object for the next thread */ 
        holder->query = reload_GeneWiseDB(NULL,holder->querydb,&db_status);  
        if( db_status == DB_RETURN_ERROR )   
          fatal("In searching EstLoop3, reload error on database query, in threads");    
        if( db_status == DB_RETURN_END ) { /*last load!*/ 
          /* End of target and query database - finished search! */ 
          close_GeneWiseDB(NULL,holder->querydb);    
          holder->search_has_ended = TRUE;   
          } /* end of last load! */ 


        /* release input mutex */ 
        if ( holder->dbsi->trace_level >= 2 )    
          fprintf(holder->dbsi->trace_file,"Releasing input lock after end of target\n");    
        if( pthread_mutex_unlock(&(holder->input_lock))!= 0 )    
          fatal("Error in releasing input lock for EstLoop3");   
        continue;    
        } /* end of End of target database */ 
      else   { /*Normal reload*/ 
        if ( holder->dbsi->trace_level >= 2 )    
          fprintf(holder->dbsi->trace_file,"Releasing input lock for normal reload\n");  
        if( pthread_mutex_unlock(&(holder->input_lock))!= 0 )    
          fatal("Error in releasing input lock for EstLoop3");   
        } /* end of Normal reload */ 
      /* get target information into datascore */ 
      dataentry_add_cDNADB(ds->target,target,holder->targetdb);  


      /* Now there is a new query/target pair ready for comparison */ 
      if ( holder->dbsi->trace_level >= 1 )  
        fprintf(holder->dbsi->trace_file,"A new pair to be compared...\n");  
      score = score_only_EstLoop3(query, target ,holder->cp);    


      if ( holder->dbsi->trace_level >= 2 )  
        fprintf(holder->dbsi->trace_file,"Getting output lock\n");   
      /* Getting lock on output */ 
      if( pthread_mutex_lock(&(holder->output_lock))!= 0 )   
        fatal("Error on getting output lock for EstLoop3");  
      /* If the score is less than cutoff, schedule the datascore for reuse */ 
      if( should_store_Hscore(holder->out,score) != TRUE)    {  
        free_DataScore(ds);  
        }  
      else   { /*storing score*/ 
        ds->score = score;   
        add_Hscore(holder->out,ds);  
        } /* end of storing score */ 
      if( pthread_mutex_unlock(&(holder->output_lock))!= 0 ) 
        fatal("Error on releasing output lock for EstLoop3");    
      if ( holder->dbsi->trace_level >= 2 )  
        fprintf(holder->dbsi->trace_file,"Released output lock\n");  


      /* Now free database objects */ 
      if ( holder->dbsi->trace_level >= 2 )  
        fprintf(holder->dbsi->trace_file,"About to get input lock for free func\n"); 
      if( pthread_mutex_lock(&(holder->input_lock))!= 0 )    
        fatal("Error on getting input lock for EstLoop3");   
      if ( holder->dbsi->trace_level >= 2 )  
        fprintf(holder->dbsi->trace_file,"Got input lock for free func\n");  
      free_GeneWiseScore(query); 
      free_ComplexSequence(target);  
      if ( holder->dbsi->trace_level >= 2 )  
        fprintf(holder->dbsi->trace_file,"Releasing input lock after free'ing\n");   
      if( pthread_mutex_unlock(&(holder->input_lock))!= 0 )  
        fatal("Error in releasing input lock for EstLoop3"); 
      } /* end of Infinite loop over all models */ 


    if ( holder->dbsi->trace_level >= 1 )    
      fprintf(holder->dbsi->trace_file,"Exiting forever loop\n");    
    return NULL; 
}    


/* Function:  serial_search_EstLoop3(out,querydb,targetdb,cp)
 *
 * Descrip:    This function makes a database search of EstLoop3
 *             It is a single processor implementation
 *
 *
 * Arg:             out [UNKN ] Undocumented argument [Hscore *]
 * Arg:         querydb [UNKN ] Undocumented argument [GeneWiseDB*]
 * Arg:        targetdb [UNKN ] Undocumented argument [cDNADB*]
 * Arg:              cp [UNKN ] Undocumented argument [cDNAParserScore *]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
#endif /* PTHREAD */ 
Search_Return_Type serial_search_EstLoop3(Hscore * out,GeneWiseDB* querydb,cDNADB* targetdb ,cDNAParserScore * cp) 
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
      warn("In searching EstLoop3, got a database reload error on the query [query] database");  
      return SEARCH_ERROR;   
      }  
    for(;;)  { /*For all query entries*/ 


      target_pos = 0;    


      target = init_cDNADB(targetdb,&db_status); 
      if( db_status == DB_RETURN_ERROR )     {  
        warn("In searching EstLoop3, got a database init error on the target [target] database");    
        return SEARCH_ERROR; 
        }  
      for(;;)    { /*For all target entries*/ 


        /* No maximum length - allocated on-the-fly */ 
        score = score_only_EstLoop3(query, target , cp);     
        if( should_store_Hscore(out,score) == TRUE )     { /*if storing datascore*/ 
          ds = new_DataScore_from_storage(out);  
          if( ds == NULL )   {  
            warn("EstLoop3 search had a memory error in allocating a new_DataScore (?a leak somewhere - DataScore is a very small datastructure");   
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
          warn("In searching EstLoop3, Reload error on database target, position %d,%d",query_pos,target_pos);   
          return SEARCH_ERROR;   
          }  
        if( db_status == DB_RETURN_END ) 
          break;/* Out of target loop */ 
        target_pos++;    
        } /* end of For all target entries */ 
      close_cDNADB(target,targetdb);     
       query = reload_GeneWiseDB(query,querydb,&db_status);  
      if( db_status == DB_RETURN_ERROR)  {  
        warn("In searching EstLoop3, Reload error on database query, position %d,%d",query_pos,target_pos);  
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


/* Function:  score_only_EstLoop3(query,target,cp)
 *
 * Descrip:    This function just calculates the score for the matrix
 *             I am pretty sure we can do this better, but hey, for the moment...
 *             It calls /allocate_EstLoop3_only
 *
 *
 * Arg:         query [UNKN ] query data structure [GeneWiseScore*]
 * Arg:        target [UNKN ] target data structure [ComplexSequence*]
 * Arg:            cp [UNKN ] Resource [cDNAParserScore *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int score_only_EstLoop3(GeneWiseScore* query,ComplexSequence* target ,cDNAParserScore * cp) 
{
    int bestscore = NEGI;    
    int i;   
    int j;   
    int k;   
    EstLoop3 * mat;  


    mat = allocate_EstLoop3_only(query, target , cp);    
    if( mat == NULL )    {  
      warn("Memory allocation error in the db search - unable to communicate to calling function. this spells DIASTER!");    
      return NEGI;   
      }  
    if((mat->basematrix = BaseMatrix_alloc_matrix_and_specials(6,(mat->leni + 1) * 3,6,3)) == NULL)  {  
      warn("Score only matrix for EstLoop3 cannot be allocated, (asking for 5  by %d  cells)",mat->leni*3);  
      mat = free_EstLoop3(mat);  
      return 0;  
      }  
    mat->basematrix->type = BASEMATRIX_TYPE_VERYSMALL;   


    /* Now, initiate matrix */ 
    for(j=0;j<7;j++) {  
      for(i=(-1);i<mat->leni;i++)    {  
        for(k=0;k<3;k++) 
          EstLoop3_VSMALL_MATRIX(mat,i,j,k) = NEGI;  
        }  
      EstLoop3_VSMALL_SPECIAL(mat,i,j,LOOP) = NEGI;  
      EstLoop3_VSMALL_SPECIAL(mat,i,j,START) = 0;    
      EstLoop3_VSMALL_SPECIAL(mat,i,j,END) = NEGI;   
      }  


    /* Ok, lets do-o-o-o-o it */ 


    for(j=0;j<mat->lenj;j++) { /*for all target positions*/ 
      auto int score;    
      auto int temp;     
      for(i=0;i<mat->leni;i++)   { /*for all query positions*/ 


        /* For state MATCH */ 
        /* setting first movement to score */ 
        score = EstLoop3_VSMALL_MATRIX(mat,i-1,j-3,MATCH) + (mat->query->seg[i]->transition[GW_MATCH2MATCH]+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)]);  
        /* From state INSERT to state MATCH */ 
        temp = EstLoop3_VSMALL_MATRIX(mat,i-1,j-3,INSERT) + (mat->query->seg[i]->transition[GW_INSERT2MATCH]+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)]);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state MATCH */ 
        temp = EstLoop3_VSMALL_MATRIX(mat,i-1,j-3,DELETE) + (mat->query->seg[i]->transition[GW_DELETE2MATCH]+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)]);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state LOOP to state MATCH */ 
        temp = EstLoop3_VSMALL_SPECIAL(mat,i-1,j-3,LOOP) + (mat->query->seg[i]->transition[GW_START2MATCH]+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)]);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state START to state MATCH */ 
        temp = EstLoop3_VSMALL_SPECIAL(mat,i-1,j-3,START) + (mat->query->seg[i]->transition[GW_START2MATCH]+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)]);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = EstLoop3_VSMALL_MATRIX(mat,i-1,j-2,MATCH) + mat->cp->trans[PCD_DELETE_1_BASE];    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = EstLoop3_VSMALL_MATRIX(mat,i-1,j-1,MATCH) + mat->cp->trans[PCD_DELETE_2_BASE];    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = EstLoop3_VSMALL_MATRIX(mat,i-1,j-4,MATCH) + mat->cp->trans[PCD_INSERT_1_BASE];    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = EstLoop3_VSMALL_MATRIX(mat,i-1,j-5,MATCH) + mat->cp->trans[PCD_INSERT_2_BASE];    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state MATCH */ 
        temp = EstLoop3_VSMALL_MATRIX(mat,i-1,j-2,INSERT) + mat->cp->trans[PCD_DELETE_1_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state MATCH */ 
        temp = EstLoop3_VSMALL_MATRIX(mat,i-1,j-1,INSERT) + mat->cp->trans[PCD_DELETE_2_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state MATCH */ 
        temp = EstLoop3_VSMALL_MATRIX(mat,i-1,j-4,INSERT) + mat->cp->trans[PCD_INSERT_1_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state MATCH */ 
        temp = EstLoop3_VSMALL_MATRIX(mat,i-1,j-5,INSERT) + mat->cp->trans[PCD_INSERT_2_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state MATCH */ 
        temp = EstLoop3_VSMALL_MATRIX(mat,i-1,j-2,DELETE) + mat->cp->trans[PCD_DELETE_1_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state MATCH */ 
        temp = EstLoop3_VSMALL_MATRIX(mat,i-1,j-1,DELETE) + mat->cp->trans[PCD_DELETE_2_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state MATCH */ 
        temp = EstLoop3_VSMALL_MATRIX(mat,i-1,j-4,DELETE) + mat->cp->trans[PCD_INSERT_1_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state MATCH */ 
        temp = EstLoop3_VSMALL_MATRIX(mat,i-1,j-5,DELETE) + mat->cp->trans[PCD_INSERT_2_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for MATCH */ 
        /* Add any movement independant score and put away */ 
         EstLoop3_VSMALL_MATRIX(mat,i,j,MATCH) = score;  


        /* state MATCH is a source for special LOOP */ 
        temp = score + (mat->query->seg[i]->transition[GW_MATCH2END]) + (0) ;    
        if( temp > EstLoop3_VSMALL_SPECIAL(mat,i,j,LOOP) )   {  
          EstLoop3_VSMALL_SPECIAL(mat,i,j,LOOP) = temp;  
          }  




        /* state MATCH is a source for special END */ 
        temp = score + (mat->query->seg[i]->transition[GW_MATCH2END]) + (0) ;    
        if( temp > EstLoop3_VSMALL_SPECIAL(mat,i,j,END) )    {  
          EstLoop3_VSMALL_SPECIAL(mat,i,j,END) = temp;   
          }  




        /* Finished calculating state MATCH */ 


        /* For state INSERT */ 
        /* setting first movement to score */ 
        score = EstLoop3_VSMALL_MATRIX(mat,i-0,j-3,MATCH) + (mat->query->seg[i]->transition[GW_MATCH2INSERT]+mat->query->seg[i]->insert[CSEQ_CDNA_CODON(mat->target,j)]);    
        /* From state INSERT to state INSERT */ 
        temp = EstLoop3_VSMALL_MATRIX(mat,i-0,j-3,INSERT) + (mat->query->seg[i]->transition[GW_INSERT2INSERT]+mat->query->seg[i]->insert[CSEQ_CDNA_CODON(mat->target,j)]);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state INSERT */ 
        temp = EstLoop3_VSMALL_MATRIX(mat,i-0,j-3,DELETE) + (mat->query->seg[i]->transition[GW_DELETE2INSERT]+mat->query->seg[i]->insert[CSEQ_CDNA_CODON(mat->target,j)]);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state LOOP to state INSERT */ 
        temp = EstLoop3_VSMALL_SPECIAL(mat,i-0,j-3,LOOP) + (mat->query->seg[i]->transition[GW_START2INSERT]+mat->query->seg[i]->insert[CSEQ_CDNA_CODON(mat->target,j)]);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state INSERT */ 
        temp = EstLoop3_VSMALL_MATRIX(mat,i-0,j-2,INSERT) + mat->cp->trans[PCD_DELETE_1_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state INSERT */ 
        temp = EstLoop3_VSMALL_MATRIX(mat,i-0,j-1,INSERT) + mat->cp->trans[PCD_DELETE_2_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state INSERT */ 
        temp = EstLoop3_VSMALL_MATRIX(mat,i-0,j-4,INSERT) + mat->cp->trans[PCD_INSERT_1_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state INSERT */ 
        temp = EstLoop3_VSMALL_MATRIX(mat,i-0,j-5,INSERT) + mat->cp->trans[PCD_INSERT_2_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state INSERT */ 
        temp = EstLoop3_VSMALL_MATRIX(mat,i-0,j-2,MATCH) + mat->cp->trans[PCD_DELETE_1_BASE];    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state INSERT */ 
        temp = EstLoop3_VSMALL_MATRIX(mat,i-0,j-1,MATCH) + mat->cp->trans[PCD_DELETE_2_BASE];    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state INSERT */ 
        temp = EstLoop3_VSMALL_MATRIX(mat,i-0,j-4,MATCH) + mat->cp->trans[PCD_INSERT_1_BASE];    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state INSERT */ 
        temp = EstLoop3_VSMALL_MATRIX(mat,i-0,j-5,MATCH) + mat->cp->trans[PCD_INSERT_2_BASE];    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state INSERT */ 
        temp = EstLoop3_VSMALL_MATRIX(mat,i-0,j-2,DELETE) + mat->cp->trans[PCD_DELETE_1_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state INSERT */ 
        temp = EstLoop3_VSMALL_MATRIX(mat,i-0,j-1,DELETE) + mat->cp->trans[PCD_DELETE_2_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state INSERT */ 
        temp = EstLoop3_VSMALL_MATRIX(mat,i-0,j-4,DELETE) + mat->cp->trans[PCD_INSERT_1_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state INSERT */ 
        temp = EstLoop3_VSMALL_MATRIX(mat,i-0,j-5,DELETE) + mat->cp->trans[PCD_INSERT_2_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for INSERT */ 
        /* Add any movement independant score and put away */ 
         EstLoop3_VSMALL_MATRIX(mat,i,j,INSERT) = score; 


        /* state INSERT is a source for special LOOP */ 
        temp = score + (mat->query->seg[i]->transition[GW_INSERT2END]) + (0) ;   
        if( temp > EstLoop3_VSMALL_SPECIAL(mat,i,j,LOOP) )   {  
          EstLoop3_VSMALL_SPECIAL(mat,i,j,LOOP) = temp;  
          }  




        /* state INSERT is a source for special END */ 
        temp = score + (mat->query->seg[i]->transition[GW_INSERT2END]) + (0) ;   
        if( temp > EstLoop3_VSMALL_SPECIAL(mat,i,j,END) )    {  
          EstLoop3_VSMALL_SPECIAL(mat,i,j,END) = temp;   
          }  




        /* Finished calculating state INSERT */ 


        /* For state DELETE */ 
        /* setting first movement to score */ 
        score = EstLoop3_VSMALL_MATRIX(mat,i-1,j-0,MATCH) + mat->query->seg[i]->transition[GW_MATCH2DELETE];     
        /* From state INSERT to state DELETE */ 
        temp = EstLoop3_VSMALL_MATRIX(mat,i-1,j-0,INSERT) + mat->query->seg[i]->transition[GW_INSERT2DELETE];    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state DELETE */ 
        temp = EstLoop3_VSMALL_MATRIX(mat,i-1,j-0,DELETE) + mat->query->seg[i]->transition[GW_DELETE2DELETE];    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state LOOP to state DELETE */ 
        temp = EstLoop3_VSMALL_SPECIAL(mat,i-1,j-0,LOOP) + mat->query->seg[i]->transition[GW_START2DELETE];  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for DELETE */ 
        /* Add any movement independant score and put away */ 
         EstLoop3_VSMALL_MATRIX(mat,i,j,DELETE) = score; 


        /* state DELETE is a source for special LOOP */ 
        temp = score + (mat->query->seg[i]->transition[GW_DELETE2END]) + (0) ;   
        if( temp > EstLoop3_VSMALL_SPECIAL(mat,i,j,LOOP) )   {  
          EstLoop3_VSMALL_SPECIAL(mat,i,j,LOOP) = temp;  
          }  




        /* state DELETE is a source for special END */ 
        temp = score + (mat->query->seg[i]->transition[GW_DELETE2END]) + (0) ;   
        if( temp > EstLoop3_VSMALL_SPECIAL(mat,i,j,END) )    {  
          EstLoop3_VSMALL_SPECIAL(mat,i,j,END) = temp;   
          }  




        /* Finished calculating state DELETE */ 
        } /* end of for all query positions */ 




      /* Special state LOOP has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = EstLoop3_VSMALL_SPECIAL(mat,0,j,LOOP); 


      /* Source START is a special source for LOOP */ 
      /* Has restricted position */ 
      if( (j-1) == 0  )  {  
        temp = EstLoop3_VSMALL_SPECIAL(mat,0,j - 1,START) + (0) + (0);   
        if( temp > score )   
          score = temp;  
        }  


      /* Source MATCH for state LOOP is not special... already calculated */ 
      /* Source INSERT for state LOOP is not special... already calculated */ 
      /* Source DELETE for state LOOP is not special... already calculated */ 
      /* Source LOOP is a special source for LOOP */ 
      temp = EstLoop3_VSMALL_SPECIAL(mat,0,j - 1,LOOP) + (0) + (0);  
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      EstLoop3_VSMALL_SPECIAL(mat,0,j,LOOP) = score; 
      /* Finished updating state LOOP */ 




      /* Special state START has no special to special movements */ 


      /* Special state END has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = EstLoop3_VSMALL_SPECIAL(mat,0,j,END);  


      /* Source LOOP is a special source for END */ 
      /* Has restricted position */ 
      if( j == mat->lenj-1 ) {  
        temp = EstLoop3_VSMALL_SPECIAL(mat,0,j - 1,LOOP) + (0) + (0);    
        if( temp > score )   
          score = temp;  
        }  


      /* Source MATCH for state END is not special... already calculated */ 
      /* Source INSERT for state END is not special... already calculated */ 
      /* Source DELETE for state END is not special... already calculated */ 
      /* Put back score... (now updated!) */ 
      EstLoop3_VSMALL_SPECIAL(mat,0,j,END) = score;  
      /* Finished updating state END */ 


      if( bestscore < EstLoop3_VSMALL_SPECIAL(mat,0,j,END) ) 
        bestscore = EstLoop3_VSMALL_SPECIAL(mat,0,j,END);    
      } /* end of for all target positions */ 


    mat = free_EstLoop3(mat);    
    return bestscore;    
}    


/* Function:  PackAln_bestmemory_EstLoop3(query,target,cp,dpenv,dpri)
 *
 * Descrip:    This function chooses the best memory set-up for the alignment
 *             using calls to basematrix, and then implements either a large
 *             or small memory model.
 *
 *             It is the best function to use if you just want an alignment
 *
 *             If you want a label alignment, you will need
 *             /convert_PackAln_to_AlnBlock_EstLoop3
 *
 *
 * Arg:         query [UNKN ] query data structure [GeneWiseScore*]
 * Arg:        target [UNKN ] target data structure [ComplexSequence*]
 * Arg:            cp [UNKN ] Resource [cDNAParserScore *]
 * Arg:         dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:          dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * PackAln_bestmemory_EstLoop3(GeneWiseScore* query,ComplexSequence* target ,cDNAParserScore * cp,DPEnvelope * dpenv,DPRunImpl * dpri) 
{
    long long total; 
    EstLoop3 * mat;  
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
        if( (mat=allocate_Expl_EstLoop3(query, target , cp,dpri)) == NULL )  {  
          warn("Unable to allocate large EstLoop3 version"); 
          return NULL;   
          }  
        calculate_dpenv_EstLoop3(mat,dpenv);     
        out =  PackAln_read_Expl_EstLoop3(mat);  
        }  
      else   {  
        mat = allocate_EstLoop3_only(query, target , cp);    
        calculate_shatter_EstLoop3(mat,dpenv);   
        out = PackAln_read_Shatter_EstLoop3(mat);    
        }  
      }  
    else {  
      if( strategy == DPIM_Linear )  {  
        /* use small implementation */ 
        if( (mat=allocate_Small_EstLoop3(query, target , cp)) == NULL )  {  
          warn("Unable to allocate small EstLoop3 version"); 
          return NULL;   
          }  
        out = PackAln_calculate_Small_EstLoop3(mat,dpenv);   
        }  
      else   {  
        /* use Large implementation */ 
        if( (mat=allocate_Expl_EstLoop3(query, target , cp,dpri)) == NULL )  {  
          warn("Unable to allocate large EstLoop3 version"); 
          return NULL;   
          }  
        if( dpri->debug == TRUE) {  
          fatal("Asked for dydebug, but dynamite file not compiled with -g. Need to recompile dynamite source"); 
          }  
        else {  
          calculate_EstLoop3(mat);   
          out =  PackAln_read_Expl_EstLoop3(mat);    
          }  
        }  
      }  


    mat = free_EstLoop3(mat);    
    return out;  
}    


/* Function:  allocate_EstLoop3_only(query,target,cp)
 *
 * Descrip:    This function only allocates the EstLoop3 structure
 *             checks types where possible and determines leni and lenj
 *             The basematrix area is delt with elsewhere
 *
 *
 * Arg:         query [UNKN ] query data structure [GeneWiseScore*]
 * Arg:        target [UNKN ] target data structure [ComplexSequence*]
 * Arg:            cp [UNKN ] Resource [cDNAParserScore *]
 *
 * Return [UNKN ]  Undocumented return value [EstLoop3 *]
 *
 */
EstLoop3 * allocate_EstLoop3_only(GeneWiseScore* query,ComplexSequence* target ,cDNAParserScore * cp) 
{
    EstLoop3 * out;  


    if((out= EstLoop3_alloc()) == NULL)  {  
      warn("Allocation of basic EstLoop3 structure failed...");  
      return NULL;   
      }  


    out->query = query;  
    out->target = target;    
    out->cp = cp;    
    out->leni = query->len;  
    out->lenj = target->seq->len;    
    return out;  
}    


/* Function:  allocate_Expl_EstLoop3(query,target,cp,dpri)
 *
 * Descrip:    This function allocates the EstLoop3 structure
 *             and the basematrix area for explicit memory implementations
 *             It calls /allocate_EstLoop3_only
 *
 *
 * Arg:         query [UNKN ] query data structure [GeneWiseScore*]
 * Arg:        target [UNKN ] target data structure [ComplexSequence*]
 * Arg:            cp [UNKN ] Resource [cDNAParserScore *]
 * Arg:          dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [EstLoop3 *]
 *
 */
EstLoop3 * allocate_Expl_EstLoop3(GeneWiseScore* query,ComplexSequence* target ,cDNAParserScore * cp,DPRunImpl * dpri) 
{
    EstLoop3 * out;  


    out = allocate_EstLoop3_only(query, target , cp);    
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
      if( (out->basematrix = BaseMatrix_alloc_matrix_and_specials((out->lenj+5)*3,(out->leni+1),3,out->lenj+5)) == NULL) {  
        warn("Explicit matrix EstLoop3 cannot be allocated, (asking for %d by %d main cells)",out->leni,out->lenj);  
        free_EstLoop3(out);  
        return NULL; 
        }  
      }  
    if( dpri->should_cache == TRUE && dpri->cache == NULL)   
      dpri->cache = hard_link_BaseMatrix(out->basematrix);   
    out->basematrix->type = BASEMATRIX_TYPE_EXPLICIT;    
    init_EstLoop3(out);  
    return out;  
}    


/* Function:  init_EstLoop3(mat)
 *
 * Descrip:    This function initates EstLoop3 matrix when in explicit mode
 *             Called in /allocate_Expl_EstLoop3
 *
 *
 * Arg:        mat [UNKN ] EstLoop3 which contains explicit basematrix memory [EstLoop3 *]
 *
 */
void init_EstLoop3(EstLoop3 * mat) 
{
    register int i;  
    register int j;  
    if( mat->basematrix->type != BASEMATRIX_TYPE_EXPLICIT)   {  
      warn("Cannot iniate matrix, is not an explicit memory type and you have assummed that");   
      return;    
      }  


    for(i= (-1);i<mat->query->len;i++)   {  
      for(j= (-5);j<6;j++)   {  
        EstLoop3_EXPL_MATRIX(mat,i,j,MATCH) = NEGI;  
        EstLoop3_EXPL_MATRIX(mat,i,j,INSERT) = NEGI; 
        EstLoop3_EXPL_MATRIX(mat,i,j,DELETE) = NEGI; 
        }  
      }  
    for(j= (-5);j<mat->target->seq->len;j++) {  
      for(i= (-1);i<2;i++)   {  
        EstLoop3_EXPL_MATRIX(mat,i,j,MATCH) = NEGI;  
        EstLoop3_EXPL_MATRIX(mat,i,j,INSERT) = NEGI; 
        EstLoop3_EXPL_MATRIX(mat,i,j,DELETE) = NEGI; 
        }  
      EstLoop3_EXPL_SPECIAL(mat,i,j,LOOP) = NEGI;    
      EstLoop3_EXPL_SPECIAL(mat,i,j,START) = 0;  
      EstLoop3_EXPL_SPECIAL(mat,i,j,END) = NEGI; 
      }  
    return;  
}    


/* Function:  recalculate_PackAln_EstLoop3(pal,mat)
 *
 * Descrip:    This function recalculates the PackAln structure produced by EstLoop3
 *             For example, in linear space methods this is used to score them
 *
 *
 * Arg:        pal [UNKN ] Undocumented argument [PackAln *]
 * Arg:        mat [UNKN ] Undocumented argument [EstLoop3 *]
 *
 */
void recalculate_PackAln_EstLoop3(PackAln * pal,EstLoop3 * mat) 
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
          if( offj == 3 && prev->state == (LOOP+3) ) {  
            pau->score = (mat->query->seg[i]->transition[GW_START2MATCH]+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)]) + (0);   
            continue;    
            }  
          if( offj == 3 && prev->state == (START+3) )    {  
            pau->score = (mat->query->seg[i]->transition[GW_START2MATCH]+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)]) + (0);   
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
          if( offj == 3 && prev->state == (LOOP+3) ) {  
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
          if( offj == 0 && prev->state == (LOOP+3) ) {  
            pau->score = mat->query->seg[i]->transition[GW_START2DELETE] + (0);  
            continue;    
            }  
          warn("In recaluclating PackAln with state DELETE, from [%d,%d,%d], got a bad source state. Error!",offi,offj,prev->state); 
          break; 
        case (LOOP+3) :  
          if( offj == 1 && prev->state == (START+3) )    {  
            pau->score = 0 + (0);    
            continue;    
            }  
          if( offj == 0 && prev->state == MATCH )    {  
            /* i here comes from the previous state ;) - not the real one */ 
            i = prev->i; 
            pau->score = mat->query->seg[i]->transition[GW_MATCH2END] + (0);     
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
          if( offj == 1 && prev->state == (LOOP+3) ) {  
            pau->score = 0 + (0);    
            continue;    
            }  
          warn("In recaluclating PackAln with state LOOP, got a bad source state. Error!");  
          break; 
        case (START+3) :     
          warn("In recaluclating PackAln with state START, got a bad source state. Error!"); 
          break; 
        case (END+3) :   
          if( offj == 1 && prev->state == (LOOP+3) ) {  
            pau->score = 0 + (0);    
            continue;    
            }  
          if( offj == 0 && prev->state == MATCH )    {  
            /* i here comes from the previous state ;) - not the real one */ 
            i = prev->i; 
            pau->score = mat->query->seg[i]->transition[GW_MATCH2END] + (0);     
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
#define EstLoop3_HIDDEN_MATRIX(thismatrix,i,j,state) (thismatrix->basematrix->matrix[(j-hiddenj+5)][(i+1)*3+state])  
#define EstLoop3_DC_SHADOW_MATRIX(thismatrix,i,j,state) (thismatrix->basematrix->matrix[((j+6)*8) % 48][(i+1)*3+state])  
#define EstLoop3_HIDDEN_SPECIAL(thismatrix,i,j,state) (thismatrix->basematrix->specmatrix[state][(j+5)]) 
#define EstLoop3_DC_SHADOW_SPECIAL(thismatrix,i,j,state) (thismatrix->basematrix->specmatrix[state*8][(j+5)])    
#define EstLoop3_DC_SHADOW_MATRIX_SP(thismatrix,i,j,state,shadow) (thismatrix->basematrix->matrix[((((j+6)*8)+(shadow+1)) % 48)][(i+1)*3 + state])   
#define EstLoop3_DC_SHADOW_SPECIAL_SP(thismatrix,i,j,state,shadow) (thismatrix->basematrix->specmatrix[state*8 +shadow+1][(j+5)])    
#define EstLoop3_DC_OPT_SHADOW_MATRIX(thismatrix,i,j,state) (score_pointers[(((j+5)% 5) * (leni+1) * 3) + ((i+1) * 3) + (state)])    
#define EstLoop3_DC_OPT_SHADOW_MATRIX_SP(thismatrix,i,j,state,shadow) (shadow_pointers[(((j+5)% 5) * (leni+1) * 24) + ((i+1) * 24) + (state * 8) + shadow+1])    
#define EstLoop3_DC_OPT_SHADOW_SPECIAL(thismatrix,i,j,state) (thismatrix->basematrix->specmatrix[state*8][(j+5)])    
/* Function:  allocate_Small_EstLoop3(query,target,cp)
 *
 * Descrip:    This function allocates the EstLoop3 structure
 *             and the basematrix area for a small memory implementations
 *             It calls /allocate_EstLoop3_only
 *
 *
 * Arg:         query [UNKN ] query data structure [GeneWiseScore*]
 * Arg:        target [UNKN ] target data structure [ComplexSequence*]
 * Arg:            cp [UNKN ] Resource [cDNAParserScore *]
 *
 * Return [UNKN ]  Undocumented return value [EstLoop3 *]
 *
 */
#define EstLoop3_DC_OPT_SHADOW_SPECIAL_SP(thismatrix,i,j,state,shadow) (thismatrix->basematrix->specmatrix[state*8 +shadow+1][(j+5)])    
EstLoop3 * allocate_Small_EstLoop3(GeneWiseScore* query,ComplexSequence* target ,cDNAParserScore * cp) 
{
    EstLoop3 * out;  


    out = allocate_EstLoop3_only(query, target , cp);    
    if( out == NULL )    
      return NULL;   
    out->basematrix = BaseMatrix_alloc_matrix_and_specials(48,(out->leni + 1) * 3,24,out->lenj+5);   
    if(out == NULL)  {  
      warn("Small shadow matrix EstLoop3 cannot be allocated, (asking for 6 by %d main cells)",out->leni+2); 
      free_EstLoop3(out);    
      return NULL;   
      }  
    out->basematrix->type = BASEMATRIX_TYPE_SHADOW;  
    return out;  
}    


/* Function:  PackAln_calculate_Small_EstLoop3(mat,dpenv)
 *
 * Descrip:    This function calculates an alignment for EstLoop3 structure in linear space
 *             If you want only the start/end points
 *             use /AlnRangeSet_calculate_Small_EstLoop3 
 *
 *             The function basically
 *               finds start/end points 
 *               foreach start/end point 
 *                 calls /full_dc_EstLoop3 
 *
 *
 * Arg:          mat [UNKN ] Undocumented argument [EstLoop3 *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * PackAln_calculate_Small_EstLoop3(EstLoop3 * mat,DPEnvelope * dpenv) 
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
      warn("Could not calculate packaln small for EstLoop3 due to wrong type of matrix");    
      return NULL;   
      }  


    out = PackAln_alloc_std();   


    start_reporting("Find start end points: ");  
    dc_optimised_start_end_calc_EstLoop3(mat,dpenv); 
    score = start_end_find_end_EstLoop3(mat,&endj);  
    out->score = score;  
    stopstate = END;
    
    /* Special to specials: have to eat up in strip and then drop back to full_dc for intervening bits */ 
    log_full_error(REPORT,0,"End at %d Score %d",endj,score);    
    stop_reporting();    
    for(;;)  { /*while there are more special bits to recover*/ 
      start_reporting("Special cell aln end   %d:",endj);    
      if( read_special_strip_EstLoop3(mat,0,endj,stopstate,&endj,&startstate,out) == FALSE ) {  
        warn("Problem in reading off special state system... going to return partial alignment");    
        break;   
        }  
      if( startstate == START || endj <= 0)  {  
        log_full_error(REPORT,0,"Recovered complete alignment"); 
        stop_reporting();    
        break;   
        }  


      log_full_error(REPORT,0,"Finished to %d",endj);    
      stop_reporting();  


      /* Ok... have to eat up another piece of matrix <sigh> */ 
      temp = startstate; 
      starti = EstLoop3_DC_SHADOW_SPECIAL_SP(mat,0,endj,temp,0); 
      startj = EstLoop3_DC_SHADOW_SPECIAL_SP(mat,0,endj,temp,1); 
      startstate = EstLoop3_DC_SHADOW_SPECIAL_SP(mat,0,endj,temp,2); 
      stopi = EstLoop3_DC_SHADOW_SPECIAL_SP(mat,0,endj,temp,3);  
      stopj = EstLoop3_DC_SHADOW_SPECIAL_SP(mat,0,endj,temp,4);  
      stopstate = EstLoop3_DC_SHADOW_SPECIAL_SP(mat,0,endj,temp,5);  


      /* Get out the score of this block. V. important! */ 
      temp = EstLoop3_DC_SHADOW_SPECIAL_SP(mat,0,endj,temp,6);   
      totalj = stopj - startj;   
      donej  = 0;    
      start_reporting("Main matrix  aln [%d,%d]:",startj,stopj);     
      if(full_dc_EstLoop3(mat,starti,startj,startstate,stopi,stopj,stopstate,out,&donej,totalj,dpenv) == FALSE)  {  
        warn("In the alignment EstLoop3 [%d,%d][%d,%d], got a problem. Please report bug ... giving you back a partial alignment",starti,startj,stopi,stopj);    
        return out;  
        }  


      /* now have to figure out which special we came from... yikes */ 
      max_matrix_to_special_EstLoop3(mat,starti,startj,startstate,temp,&stopi,&stopj,&stopstate,&temp,NULL); 
      if( stopi == EstLoop3_READ_OFF_ERROR)  {  
        warn("In EstLoop3 read off ending at %d ... got a bad matrix to special read off... returning partial alignment",startj);    
        invert_PackAln(out); 
        recalculate_PackAln_EstLoop3(out,mat);   
        return out;  
        }  
      /* if at start, break, otherwise, back to eat another strip */ 
      if( stopstate == START)    {  
        log_full_error(REPORT,0,"Recovered complete alignment      ");   
        stop_reporting();    
        break;   
        }  
      log_full_error(REPORT,0,"Finished  alignment to %d           ",startj);    
      stop_reporting();  
      endj = stopj;  
      /* stopstate is correct as it is */ 
      } /* end of while there are more special bits to recover */ 
    invert_PackAln(out); 
    recalculate_PackAln_EstLoop3(out,mat);   
    return out;  


}    


/* Function:  AlnRangeSet_calculate_Small_EstLoop3(mat)
 *
 * Descrip:    This function calculates an alignment for EstLoop3 structure in linear space
 *             If you want the full alignment, use /PackAln_calculate_Small_EstLoop3 
 *             If you have already got the full alignment, but want the range set, use /AlnRangeSet_from_PackAln_EstLoop3
 *             If you have got the small matrix but not the alignment, use /AlnRangeSet_from_EstLoop3 
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [EstLoop3 *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * AlnRangeSet_calculate_Small_EstLoop3(EstLoop3 * mat) 
{
    AlnRangeSet * out;   


    start_reporting("Find start end points: ");  
    dc_optimised_start_end_calc_EstLoop3(mat,NULL);  
    log_full_error(REPORT,0,"Calculated");   


    out = AlnRangeSet_from_EstLoop3(mat);    
    return out;  
}    


/* Function:  AlnRangeSet_from_EstLoop3(mat)
 *
 * Descrip:    This function reads off a start/end structure
 *             for EstLoop3 structure in linear space
 *             If you want the full alignment use
 *             /PackAln_calculate_Small_EstLoop3 
 *             If you have not calculated the matrix use
 *             /AlnRange_calculate_Small_EstLoop3
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [EstLoop3 *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * AlnRangeSet_from_EstLoop3(EstLoop3 * mat) 
{
    AlnRangeSet * out;   
    AlnRange * temp; 
    int jpos;    
    int state;   


    if( mat->basematrix->type != BASEMATRIX_TYPE_SHADOW) {  
      warn("Bad error! - non shadow matrix type in AlnRangeSet_from_EstLoop3");  
      return NULL;   
      }  


    out = AlnRangeSet_alloc_std();   
    /* Find the end position */ 
    out->score = start_end_find_end_EstLoop3(mat,&jpos); 
    state = END; 


    while( (temp = AlnRange_build_EstLoop3(mat,jpos,state,&jpos,&state)) != NULL)    
      add_AlnRangeSet(out,temp); 
    return out;  
}    


/* Function:  AlnRange_build_EstLoop3(mat,stopj,stopspecstate,startj,startspecstate)
 *
 * Descrip:    This function calculates a single start/end set in linear space
 *             Really a sub-routine for /AlnRangeSet_from_PackAln_EstLoop3
 *
 *
 * Arg:                   mat [UNKN ] Undocumented argument [EstLoop3 *]
 * Arg:                 stopj [UNKN ] Undocumented argument [int]
 * Arg:         stopspecstate [UNKN ] Undocumented argument [int]
 * Arg:                startj [UNKN ] Undocumented argument [int *]
 * Arg:        startspecstate [UNKN ] Undocumented argument [int *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRange *]
 *
 */
AlnRange * AlnRange_build_EstLoop3(EstLoop3 * mat,int stopj,int stopspecstate,int * startj,int * startspecstate) 
{
    AlnRange * out;  
    int jpos;    
    int state;   


    if( mat->basematrix->type != BASEMATRIX_TYPE_SHADOW) {  
      warn("Bad error! - non shadow matrix type in AlnRangeSet_from_EstLoop3");  
      return NULL;   
      }  


    /* Assumme that we have specials (we should!). Read back along the specials till we have the finish point */ 
    if( read_special_strip_EstLoop3(mat,0,stopj,stopspecstate,&jpos,&state,NULL) == FALSE)   {  
      warn("In AlnRanger_build_EstLoop3 alignment ending at %d, unable to read back specials. Will (evenutally) return a partial range set... BEWARE!",stopj);   
      return NULL;   
      }  
    if( state == START || jpos <= 0) 
      return NULL;   


    out = AlnRange_alloc();  


    out->starti = EstLoop3_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,0); 
    out->startj = EstLoop3_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,1); 
    out->startstate = EstLoop3_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,2); 
    out->stopi = EstLoop3_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,3);  
    out->stopj = EstLoop3_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,4);  
    out->stopstate = EstLoop3_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,5);  
    out->startscore = EstLoop3_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,6); 
    out->stopscore = EstLoop3_DC_SHADOW_SPECIAL(mat,0,jpos,state);   


    /* Now, we have to figure out where this state came from in the specials */ 
    max_matrix_to_special_EstLoop3(mat,out->starti,out->startj,out->startstate,out->startscore,&jpos,startj,startspecstate,&state,NULL); 
    if( jpos == EstLoop3_READ_OFF_ERROR) {  
      warn("In AlnRange_build_EstLoop3 alignment ending at %d, with aln range between %d-%d in j, unable to find source special, returning this range, but this could get tricky!",stopj,out->startj,out->stopj);    
      return out;    
      }  


    /* Put in the correct score for startstate, from the special */ 
    out->startscore = EstLoop3_DC_SHADOW_SPECIAL(mat,0,*startj,*startspecstate); 
    /* The correct j coords have been put into startj, startspecstate... so just return out */ 
    return out;  
}    


/* Function:  read_hidden_EstLoop3(mat,starti,startj,startstate,stopi,stopj,stopstate,out)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [EstLoop3 *]
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
boolean read_hidden_EstLoop3(EstLoop3 * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out) 
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


      max_hidden_EstLoop3(mat,startj,i,j,state,isspecial,&i,&j,&state,&isspecial,&cellscore);    


      if( i == EstLoop3_READ_OFF_ERROR)  {  
        warn("In EstLoop3 hidden read off, between %d:%d,%d:%d - at got bad read off. Problem!",starti,startj,stopi,stopj);  
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
        warn("In EstLoop3 hidden read off, between %d:%d,%d:%d - hit start cell, but not in start state. Can't be good!.",starti,startj,stopi,stopj);    
        return FALSE;    
        }  
      }  
    warn("In EstLoop3 hidden read off, between %d:%d,%d:%d - gone past start cell (now in %d,%d,%d), can't be good news!.",starti,startj,stopi,stopj,i,j,state); 
    return FALSE;    
}    


/* Function:  max_hidden_EstLoop3(mat,hiddenj,i,j,state,isspecial,reti,retj,retstate,retspecial,cellscore)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [EstLoop3 *]
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
int max_hidden_EstLoop3(EstLoop3 * mat,int hiddenj,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore) 
{
    register int temp;   
    register int cscore; 


    *reti = (*retj) = (*retstate) = EstLoop3_READ_OFF_ERROR; 


    if( i < 0 || j < 0 || i > mat->query->len || j > mat->target->seq->len)  {  
      warn("In EstLoop3 matrix special read off - out of bounds on matrix [i,j is %d,%d state %d in standard matrix]",i,j,state);    
      return -1; 
      }  


    /* Then you have to select the correct switch statement to figure out the readoff      */ 
    /* Somewhat odd - reverse the order of calculation and return as soon as it is correct */ 
    cscore = EstLoop3_HIDDEN_MATRIX(mat,i,j,state);  
    switch(state)    { /*Switch state */ 
      case MATCH :   
        temp = cscore - (mat->cp->trans[PCD_INSERT_2_BASE]) -  (0);  
        if( temp == EstLoop3_HIDDEN_MATRIX(mat,i - 1,j - 5,DELETE) ) {  
          *reti = i - 1; 
          *retj = j - 5; 
          *retstate = DELETE;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - EstLoop3_HIDDEN_MATRIX(mat,i-1,j-5,DELETE);    
            }  
          return EstLoop3_HIDDEN_MATRIX(mat,i - 1,j - 5,DELETE);     
          }  
        temp = cscore - (mat->cp->trans[PCD_INSERT_1_BASE]) -  (0);  
        if( temp == EstLoop3_HIDDEN_MATRIX(mat,i - 1,j - 4,DELETE) ) {  
          *reti = i - 1; 
          *retj = j - 4; 
          *retstate = DELETE;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - EstLoop3_HIDDEN_MATRIX(mat,i-1,j-4,DELETE);    
            }  
          return EstLoop3_HIDDEN_MATRIX(mat,i - 1,j - 4,DELETE);     
          }  
        temp = cscore - (mat->cp->trans[PCD_DELETE_2_BASE]) -  (0);  
        if( temp == EstLoop3_HIDDEN_MATRIX(mat,i - 1,j - 1,DELETE) ) {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = DELETE;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - EstLoop3_HIDDEN_MATRIX(mat,i-1,j-1,DELETE);    
            }  
          return EstLoop3_HIDDEN_MATRIX(mat,i - 1,j - 1,DELETE);     
          }  
        temp = cscore - (mat->cp->trans[PCD_DELETE_1_BASE]) -  (0);  
        if( temp == EstLoop3_HIDDEN_MATRIX(mat,i - 1,j - 2,DELETE) ) {  
          *reti = i - 1; 
          *retj = j - 2; 
          *retstate = DELETE;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - EstLoop3_HIDDEN_MATRIX(mat,i-1,j-2,DELETE);    
            }  
          return EstLoop3_HIDDEN_MATRIX(mat,i - 1,j - 2,DELETE);     
          }  
        temp = cscore - (mat->cp->trans[PCD_INSERT_2_BASE]) -  (0);  
        if( temp == EstLoop3_HIDDEN_MATRIX(mat,i - 1,j - 5,INSERT) ) {  
          *reti = i - 1; 
          *retj = j - 5; 
          *retstate = INSERT;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - EstLoop3_HIDDEN_MATRIX(mat,i-1,j-5,INSERT);    
            }  
          return EstLoop3_HIDDEN_MATRIX(mat,i - 1,j - 5,INSERT);     
          }  
        temp = cscore - (mat->cp->trans[PCD_INSERT_1_BASE]) -  (0);  
        if( temp == EstLoop3_HIDDEN_MATRIX(mat,i - 1,j - 4,INSERT) ) {  
          *reti = i - 1; 
          *retj = j - 4; 
          *retstate = INSERT;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - EstLoop3_HIDDEN_MATRIX(mat,i-1,j-4,INSERT);    
            }  
          return EstLoop3_HIDDEN_MATRIX(mat,i - 1,j - 4,INSERT);     
          }  
        temp = cscore - (mat->cp->trans[PCD_DELETE_2_BASE]) -  (0);  
        if( temp == EstLoop3_HIDDEN_MATRIX(mat,i - 1,j - 1,INSERT) ) {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = INSERT;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - EstLoop3_HIDDEN_MATRIX(mat,i-1,j-1,INSERT);    
            }  
          return EstLoop3_HIDDEN_MATRIX(mat,i - 1,j - 1,INSERT);     
          }  
        temp = cscore - (mat->cp->trans[PCD_DELETE_1_BASE]) -  (0);  
        if( temp == EstLoop3_HIDDEN_MATRIX(mat,i - 1,j - 2,INSERT) ) {  
          *reti = i - 1; 
          *retj = j - 2; 
          *retstate = INSERT;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - EstLoop3_HIDDEN_MATRIX(mat,i-1,j-2,INSERT);    
            }  
          return EstLoop3_HIDDEN_MATRIX(mat,i - 1,j - 2,INSERT);     
          }  
        temp = cscore - (mat->cp->trans[PCD_INSERT_2_BASE]) -  (0);  
        if( temp == EstLoop3_HIDDEN_MATRIX(mat,i - 1,j - 5,MATCH) )  {  
          *reti = i - 1; 
          *retj = j - 5; 
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - EstLoop3_HIDDEN_MATRIX(mat,i-1,j-5,MATCH); 
            }  
          return EstLoop3_HIDDEN_MATRIX(mat,i - 1,j - 5,MATCH);  
          }  
        temp = cscore - (mat->cp->trans[PCD_INSERT_1_BASE]) -  (0);  
        if( temp == EstLoop3_HIDDEN_MATRIX(mat,i - 1,j - 4,MATCH) )  {  
          *reti = i - 1; 
          *retj = j - 4; 
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - EstLoop3_HIDDEN_MATRIX(mat,i-1,j-4,MATCH); 
            }  
          return EstLoop3_HIDDEN_MATRIX(mat,i - 1,j - 4,MATCH);  
          }  
        temp = cscore - (mat->cp->trans[PCD_DELETE_2_BASE]) -  (0);  
        if( temp == EstLoop3_HIDDEN_MATRIX(mat,i - 1,j - 1,MATCH) )  {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - EstLoop3_HIDDEN_MATRIX(mat,i-1,j-1,MATCH); 
            }  
          return EstLoop3_HIDDEN_MATRIX(mat,i - 1,j - 1,MATCH);  
          }  
        temp = cscore - (mat->cp->trans[PCD_DELETE_1_BASE]) -  (0);  
        if( temp == EstLoop3_HIDDEN_MATRIX(mat,i - 1,j - 2,MATCH) )  {  
          *reti = i - 1; 
          *retj = j - 2; 
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - EstLoop3_HIDDEN_MATRIX(mat,i-1,j-2,MATCH); 
            }  
          return EstLoop3_HIDDEN_MATRIX(mat,i - 1,j - 2,MATCH);  
          }  
        /* Not allowing special sources.. skipping START */ 
        /* Not allowing special sources.. skipping LOOP */ 
        temp = cscore - ((mat->query->seg[i]->transition[GW_DELETE2MATCH]+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)])) -  (0);    
        if( temp == EstLoop3_HIDDEN_MATRIX(mat,i - 1,j - 3,DELETE) ) {  
          *reti = i - 1; 
          *retj = j - 3; 
          *retstate = DELETE;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - EstLoop3_HIDDEN_MATRIX(mat,i-1,j-3,DELETE);    
            }  
          return EstLoop3_HIDDEN_MATRIX(mat,i - 1,j - 3,DELETE);     
          }  
        temp = cscore - ((mat->query->seg[i]->transition[GW_INSERT2MATCH]+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)])) -  (0);    
        if( temp == EstLoop3_HIDDEN_MATRIX(mat,i - 1,j - 3,INSERT) ) {  
          *reti = i - 1; 
          *retj = j - 3; 
          *retstate = INSERT;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - EstLoop3_HIDDEN_MATRIX(mat,i-1,j-3,INSERT);    
            }  
          return EstLoop3_HIDDEN_MATRIX(mat,i - 1,j - 3,INSERT);     
          }  
        temp = cscore - ((mat->query->seg[i]->transition[GW_MATCH2MATCH]+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)])) -  (0); 
        if( temp == EstLoop3_HIDDEN_MATRIX(mat,i - 1,j - 3,MATCH) )  {  
          *reti = i - 1; 
          *retj = j - 3; 
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - EstLoop3_HIDDEN_MATRIX(mat,i-1,j-3,MATCH); 
            }  
          return EstLoop3_HIDDEN_MATRIX(mat,i - 1,j - 3,MATCH);  
          }  
        warn("Major problem (!) - in EstLoop3 read off, position %d,%d state %d no source found!",i,j,state);    
        return (-1); 
      case INSERT :  
        temp = cscore - (mat->cp->trans[PCD_INSERT_2_BASE]) -  (0);  
        if( temp == EstLoop3_HIDDEN_MATRIX(mat,i - 0,j - 5,DELETE) ) {  
          *reti = i - 0; 
          *retj = j - 5; 
          *retstate = DELETE;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - EstLoop3_HIDDEN_MATRIX(mat,i-0,j-5,DELETE);    
            }  
          return EstLoop3_HIDDEN_MATRIX(mat,i - 0,j - 5,DELETE);     
          }  
        temp = cscore - (mat->cp->trans[PCD_INSERT_1_BASE]) -  (0);  
        if( temp == EstLoop3_HIDDEN_MATRIX(mat,i - 0,j - 4,DELETE) ) {  
          *reti = i - 0; 
          *retj = j - 4; 
          *retstate = DELETE;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - EstLoop3_HIDDEN_MATRIX(mat,i-0,j-4,DELETE);    
            }  
          return EstLoop3_HIDDEN_MATRIX(mat,i - 0,j - 4,DELETE);     
          }  
        temp = cscore - (mat->cp->trans[PCD_DELETE_2_BASE]) -  (0);  
        if( temp == EstLoop3_HIDDEN_MATRIX(mat,i - 0,j - 1,DELETE) ) {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = DELETE;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - EstLoop3_HIDDEN_MATRIX(mat,i-0,j-1,DELETE);    
            }  
          return EstLoop3_HIDDEN_MATRIX(mat,i - 0,j - 1,DELETE);     
          }  
        temp = cscore - (mat->cp->trans[PCD_DELETE_1_BASE]) -  (0);  
        if( temp == EstLoop3_HIDDEN_MATRIX(mat,i - 0,j - 2,DELETE) ) {  
          *reti = i - 0; 
          *retj = j - 2; 
          *retstate = DELETE;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - EstLoop3_HIDDEN_MATRIX(mat,i-0,j-2,DELETE);    
            }  
          return EstLoop3_HIDDEN_MATRIX(mat,i - 0,j - 2,DELETE);     
          }  
        temp = cscore - (mat->cp->trans[PCD_INSERT_2_BASE]) -  (0);  
        if( temp == EstLoop3_HIDDEN_MATRIX(mat,i - 0,j - 5,MATCH) )  {  
          *reti = i - 0; 
          *retj = j - 5; 
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - EstLoop3_HIDDEN_MATRIX(mat,i-0,j-5,MATCH); 
            }  
          return EstLoop3_HIDDEN_MATRIX(mat,i - 0,j - 5,MATCH);  
          }  
        temp = cscore - (mat->cp->trans[PCD_INSERT_1_BASE]) -  (0);  
        if( temp == EstLoop3_HIDDEN_MATRIX(mat,i - 0,j - 4,MATCH) )  {  
          *reti = i - 0; 
          *retj = j - 4; 
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - EstLoop3_HIDDEN_MATRIX(mat,i-0,j-4,MATCH); 
            }  
          return EstLoop3_HIDDEN_MATRIX(mat,i - 0,j - 4,MATCH);  
          }  
        temp = cscore - (mat->cp->trans[PCD_DELETE_2_BASE]) -  (0);  
        if( temp == EstLoop3_HIDDEN_MATRIX(mat,i - 0,j - 1,MATCH) )  {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - EstLoop3_HIDDEN_MATRIX(mat,i-0,j-1,MATCH); 
            }  
          return EstLoop3_HIDDEN_MATRIX(mat,i - 0,j - 1,MATCH);  
          }  
        temp = cscore - (mat->cp->trans[PCD_DELETE_1_BASE]) -  (0);  
        if( temp == EstLoop3_HIDDEN_MATRIX(mat,i - 0,j - 2,MATCH) )  {  
          *reti = i - 0; 
          *retj = j - 2; 
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - EstLoop3_HIDDEN_MATRIX(mat,i-0,j-2,MATCH); 
            }  
          return EstLoop3_HIDDEN_MATRIX(mat,i - 0,j - 2,MATCH);  
          }  
        temp = cscore - (mat->cp->trans[PCD_INSERT_2_BASE]) -  (0);  
        if( temp == EstLoop3_HIDDEN_MATRIX(mat,i - 0,j - 5,INSERT) ) {  
          *reti = i - 0; 
          *retj = j - 5; 
          *retstate = INSERT;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - EstLoop3_HIDDEN_MATRIX(mat,i-0,j-5,INSERT);    
            }  
          return EstLoop3_HIDDEN_MATRIX(mat,i - 0,j - 5,INSERT);     
          }  
        temp = cscore - (mat->cp->trans[PCD_INSERT_1_BASE]) -  (0);  
        if( temp == EstLoop3_HIDDEN_MATRIX(mat,i - 0,j - 4,INSERT) ) {  
          *reti = i - 0; 
          *retj = j - 4; 
          *retstate = INSERT;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - EstLoop3_HIDDEN_MATRIX(mat,i-0,j-4,INSERT);    
            }  
          return EstLoop3_HIDDEN_MATRIX(mat,i - 0,j - 4,INSERT);     
          }  
        temp = cscore - (mat->cp->trans[PCD_DELETE_2_BASE]) -  (0);  
        if( temp == EstLoop3_HIDDEN_MATRIX(mat,i - 0,j - 1,INSERT) ) {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = INSERT;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - EstLoop3_HIDDEN_MATRIX(mat,i-0,j-1,INSERT);    
            }  
          return EstLoop3_HIDDEN_MATRIX(mat,i - 0,j - 1,INSERT);     
          }  
        temp = cscore - (mat->cp->trans[PCD_DELETE_1_BASE]) -  (0);  
        if( temp == EstLoop3_HIDDEN_MATRIX(mat,i - 0,j - 2,INSERT) ) {  
          *reti = i - 0; 
          *retj = j - 2; 
          *retstate = INSERT;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - EstLoop3_HIDDEN_MATRIX(mat,i-0,j-2,INSERT);    
            }  
          return EstLoop3_HIDDEN_MATRIX(mat,i - 0,j - 2,INSERT);     
          }  
        /* Not allowing special sources.. skipping LOOP */ 
        temp = cscore - ((mat->query->seg[i]->transition[GW_DELETE2INSERT]+mat->query->seg[i]->insert[CSEQ_CDNA_CODON(mat->target,j)])) -  (0);  
        if( temp == EstLoop3_HIDDEN_MATRIX(mat,i - 0,j - 3,DELETE) ) {  
          *reti = i - 0; 
          *retj = j - 3; 
          *retstate = DELETE;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - EstLoop3_HIDDEN_MATRIX(mat,i-0,j-3,DELETE);    
            }  
          return EstLoop3_HIDDEN_MATRIX(mat,i - 0,j - 3,DELETE);     
          }  
        temp = cscore - ((mat->query->seg[i]->transition[GW_INSERT2INSERT]+mat->query->seg[i]->insert[CSEQ_CDNA_CODON(mat->target,j)])) -  (0);  
        if( temp == EstLoop3_HIDDEN_MATRIX(mat,i - 0,j - 3,INSERT) ) {  
          *reti = i - 0; 
          *retj = j - 3; 
          *retstate = INSERT;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - EstLoop3_HIDDEN_MATRIX(mat,i-0,j-3,INSERT);    
            }  
          return EstLoop3_HIDDEN_MATRIX(mat,i - 0,j - 3,INSERT);     
          }  
        temp = cscore - ((mat->query->seg[i]->transition[GW_MATCH2INSERT]+mat->query->seg[i]->insert[CSEQ_CDNA_CODON(mat->target,j)])) -  (0);   
        if( temp == EstLoop3_HIDDEN_MATRIX(mat,i - 0,j - 3,MATCH) )  {  
          *reti = i - 0; 
          *retj = j - 3; 
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - EstLoop3_HIDDEN_MATRIX(mat,i-0,j-3,MATCH); 
            }  
          return EstLoop3_HIDDEN_MATRIX(mat,i - 0,j - 3,MATCH);  
          }  
        warn("Major problem (!) - in EstLoop3 read off, position %d,%d state %d no source found!",i,j,state);    
        return (-1); 
      case DELETE :  
        /* Not allowing special sources.. skipping LOOP */ 
        temp = cscore - (mat->query->seg[i]->transition[GW_DELETE2DELETE]) -  (0);   
        if( temp == EstLoop3_HIDDEN_MATRIX(mat,i - 1,j - 0,DELETE) ) {  
          *reti = i - 1; 
          *retj = j - 0; 
          *retstate = DELETE;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - EstLoop3_HIDDEN_MATRIX(mat,i-1,j-0,DELETE);    
            }  
          return EstLoop3_HIDDEN_MATRIX(mat,i - 1,j - 0,DELETE);     
          }  
        temp = cscore - (mat->query->seg[i]->transition[GW_INSERT2DELETE]) -  (0);   
        if( temp == EstLoop3_HIDDEN_MATRIX(mat,i - 1,j - 0,INSERT) ) {  
          *reti = i - 1; 
          *retj = j - 0; 
          *retstate = INSERT;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - EstLoop3_HIDDEN_MATRIX(mat,i-1,j-0,INSERT);    
            }  
          return EstLoop3_HIDDEN_MATRIX(mat,i - 1,j - 0,INSERT);     
          }  
        temp = cscore - (mat->query->seg[i]->transition[GW_MATCH2DELETE]) -  (0);    
        if( temp == EstLoop3_HIDDEN_MATRIX(mat,i - 1,j - 0,MATCH) )  {  
          *reti = i - 1; 
          *retj = j - 0; 
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - EstLoop3_HIDDEN_MATRIX(mat,i-1,j-0,MATCH); 
            }  
          return EstLoop3_HIDDEN_MATRIX(mat,i - 1,j - 0,MATCH);  
          }  
        warn("Major problem (!) - in EstLoop3 read off, position %d,%d state %d no source found!",i,j,state);    
        return (-1); 
      default:   
        warn("Major problem (!) - in EstLoop3 read off, position %d,%d state %d no source found!",i,j,state);    
        return (-1); 
      } /* end of Switch state  */ 
}    


/* Function:  read_special_strip_EstLoop3(mat,stopi,stopj,stopstate,startj,startstate,out)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [EstLoop3 *]
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
boolean read_special_strip_EstLoop3(EstLoop3 * mat,int stopi,int stopj,int stopstate,int * startj,int * startstate,PackAln * out) 
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
    while( j > EstLoop3_DC_SHADOW_SPECIAL_SP(mat,i,j,state,4) && state != START) { /*while more specials to eat up*/ 
      /* Put away current state, if we should */ 
      if(out != NULL)    {  
        pau = PackAlnUnit_alloc();  /* Should deal with memory overflow */ 
        pau->i = i;  
        pau->j = j;  
        pau->state =  state + 3; 
        add_PackAln(out,pau);    
        }  


      max_special_strip_EstLoop3(mat,i,j,state,isspecial,&i,&j,&state,&isspecial,&cellscore);    
      if( i == EstLoop3_READ_OFF_ERROR)  {  
        warn("In special strip read EstLoop3, got a bad read off error. Sorry!");    
        return FALSE;    
        }  
      } /* end of while more specials to eat up */ 


    /* check to see we have not gone too far! */ 
    if( state != START && j < EstLoop3_DC_SHADOW_SPECIAL_SP(mat,i,j,state,4))    {  
      warn("In special strip read EstLoop3, at special [%d] state [%d] overshot!",j,state);  
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


/* Function:  max_special_strip_EstLoop3(mat,i,j,state,isspecial,reti,retj,retstate,retspecial,cellscore)
 *
 * Descrip:    A pretty intense internal function. Deals with read-off only in specials
 *
 *
 * Arg:               mat [UNKN ] Undocumented argument [EstLoop3 *]
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
int max_special_strip_EstLoop3(EstLoop3 * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore) 
{
    int temp;    
    int cscore;  


    *reti = (*retj) = (*retstate) = EstLoop3_READ_OFF_ERROR; 
    if( isspecial == FALSE ) {  
      warn("In special strip max function for EstLoop3, got a non special start point. Problem! (bad!)");    
      return (-1);   
      }  


    if( j < 0 || j > mat->target->seq->len)  {  
      warn("In EstLoop3 matrix special read off - out of bounds on matrix [j is %d in special]",j);  
      return -1; 
      }  


    cscore = EstLoop3_DC_SHADOW_SPECIAL(mat,i,j,state);  
    switch(state)    { /*switch on special states*/ 
      case LOOP :    
        /* source LOOP is a special */ 
        temp = cscore - (0) - (0);   
        if( temp == EstLoop3_DC_SHADOW_SPECIAL(mat,i - 0,j - 1,LOOP) )   {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = LOOP;  
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - EstLoop3_DC_SHADOW_SPECIAL(mat,i-0,j-1,LOOP);  
            }  
          return EstLoop3_DC_SHADOW_MATRIX(mat,i - 0,j - 1,LOOP) ;   
          }  
        /* Source DELETE is not a special */ 
        /* Source INSERT is not a special */ 
        /* Source MATCH is not a special */ 
        /* source START is a special */ 
        temp = cscore - (0) - (0);   
        if( temp == EstLoop3_DC_SHADOW_SPECIAL(mat,i - 0,j - 1,START) )  {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = START; 
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - EstLoop3_DC_SHADOW_SPECIAL(mat,i-0,j-1,START);     
            }  
          return EstLoop3_DC_SHADOW_MATRIX(mat,i - 0,j - 1,START) ;  
          }  
      case START :   
      case END :     
        /* Source DELETE is not a special */ 
        /* Source INSERT is not a special */ 
        /* Source MATCH is not a special */ 
        /* source LOOP is a special */ 
        temp = cscore - (0) - (0);   
        if( temp == EstLoop3_DC_SHADOW_SPECIAL(mat,i - 0,j - 1,LOOP) )   {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = LOOP;  
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - EstLoop3_DC_SHADOW_SPECIAL(mat,i-0,j-1,LOOP);  
            }  
          return EstLoop3_DC_SHADOW_MATRIX(mat,i - 0,j - 1,LOOP) ;   
          }  
      default:   
        warn("Major problem (!) - in EstLoop3 special strip read off, position %d,%d state %d no source found  dropped into default on source switch!",i,j,state);   
        return (-1); 
      } /* end of switch on special states */ 
}    


/* Function:  max_matrix_to_special_EstLoop3(mat,i,j,state,cscore,reti,retj,retstate,retspecial,cellscore)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [EstLoop3 *]
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
int max_matrix_to_special_EstLoop3(EstLoop3 * mat,int i,int j,int state,int cscore,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore) 
{
    int temp;    
    *reti = (*retj) = (*retstate) = EstLoop3_READ_OFF_ERROR; 


    if( j < 0 || j > mat->lenj)  {  
      warn("In EstLoop3 matrix to special read off - out of bounds on matrix [j is %d in special]",j);   
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
        temp = cscore - ((mat->query->seg[i]->transition[GW_START2MATCH]+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)])) -  (0);     
        if( temp == EstLoop3_DC_SHADOW_SPECIAL(mat,i - 1,j - 3,START) )  {  
          *reti = i - 1; 
          *retj = j - 3; 
          *retstate = START; 
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - EstLoop3_DC_SHADOW_SPECIAL(mat,i-1,j-3,START);     
            }  
          return EstLoop3_DC_SHADOW_MATRIX(mat,i - 1,j - 3,START) ;  
          }  
        temp = cscore - ((mat->query->seg[i]->transition[GW_START2MATCH]+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)])) -  (0);     
        if( temp == EstLoop3_DC_SHADOW_SPECIAL(mat,i - 1,j - 3,LOOP) )   {  
          *reti = i - 1; 
          *retj = j - 3; 
          *retstate = LOOP;  
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - EstLoop3_DC_SHADOW_SPECIAL(mat,i-1,j-3,LOOP);  
            }  
          return EstLoop3_DC_SHADOW_MATRIX(mat,i - 1,j - 3,LOOP) ;   
          }  
        /* Source DELETE is not a special, should not get here! */ 
        /* Source INSERT is not a special, should not get here! */ 
        /* Source MATCH is not a special, should not get here! */ 
        warn("Major problem (!) - in EstLoop3 matrix to special read off, position %d,%d state %d no source found!",i,j,state);  
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
        if( temp == EstLoop3_DC_SHADOW_SPECIAL(mat,i - 0,j - 3,LOOP) )   {  
          *reti = i - 0; 
          *retj = j - 3; 
          *retstate = LOOP;  
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - EstLoop3_DC_SHADOW_SPECIAL(mat,i-0,j-3,LOOP);  
            }  
          return EstLoop3_DC_SHADOW_MATRIX(mat,i - 0,j - 3,LOOP) ;   
          }  
        /* Source DELETE is not a special, should not get here! */ 
        /* Source INSERT is not a special, should not get here! */ 
        /* Source MATCH is not a special, should not get here! */ 
        warn("Major problem (!) - in EstLoop3 matrix to special read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      case DELETE :  
        temp = cscore - (mat->query->seg[i]->transition[GW_START2DELETE]) -  (0);    
        if( temp == EstLoop3_DC_SHADOW_SPECIAL(mat,i - 1,j - 0,LOOP) )   {  
          *reti = i - 1; 
          *retj = j - 0; 
          *retstate = LOOP;  
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - EstLoop3_DC_SHADOW_SPECIAL(mat,i-1,j-0,LOOP);  
            }  
          return EstLoop3_DC_SHADOW_MATRIX(mat,i - 1,j - 0,LOOP) ;   
          }  
        /* Source DELETE is not a special, should not get here! */ 
        /* Source INSERT is not a special, should not get here! */ 
        /* Source MATCH is not a special, should not get here! */ 
        warn("Major problem (!) - in EstLoop3 matrix to special read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      default:   
        warn("Major problem (!) - in EstLoop3 read off, position %d,%d state %d no source found!",i,j,state);    
        return (-1); 
      } /* end of Switch state  */ 


}    


/* Function:  calculate_hidden_EstLoop3(mat,starti,startj,startstate,stopi,stopj,stopstate,dpenv)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [EstLoop3 *]
 * Arg:            starti [UNKN ] Undocumented argument [int]
 * Arg:            startj [UNKN ] Undocumented argument [int]
 * Arg:        startstate [UNKN ] Undocumented argument [int]
 * Arg:             stopi [UNKN ] Undocumented argument [int]
 * Arg:             stopj [UNKN ] Undocumented argument [int]
 * Arg:         stopstate [UNKN ] Undocumented argument [int]
 * Arg:             dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 */
void calculate_hidden_EstLoop3(EstLoop3 * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv) 
{
    register int i;  
    register int j;  
    register int score;  
    register int temp;   
    register int hiddenj;    


    hiddenj = startj;    


    init_hidden_EstLoop3(mat,starti,startj,stopi,stopj);     


    EstLoop3_HIDDEN_MATRIX(mat,starti,startj,startstate) = 0;    


    for(j=startj;j<=stopj;j++)   {  
      for(i=starti;i<=stopi;i++) {  
        /* Should *not* do very first cell as this is the one set to zero in one state! */ 
        if( i == starti && j == startj ) 
          continue;  
        if( dpenv != NULL && is_in_DPEnvelope(dpenv,i,j) == FALSE )  { /*Is not in envelope*/ 
          EstLoop3_HIDDEN_MATRIX(mat,i,j,MATCH) = NEGI;  
          EstLoop3_HIDDEN_MATRIX(mat,i,j,INSERT) = NEGI;     
          EstLoop3_HIDDEN_MATRIX(mat,i,j,DELETE) = NEGI;     
          continue;  
          } /* end of Is not in envelope */ 


        /* For state MATCH */ 
        /* setting first movement to score */ 
        score = EstLoop3_HIDDEN_MATRIX(mat,i-1,j-3,MATCH) + (mat->query->seg[i]->transition[GW_MATCH2MATCH]+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)]);  
        /* From state INSERT to state MATCH */ 
        temp = EstLoop3_HIDDEN_MATRIX(mat,i-1,j-3,INSERT) + (mat->query->seg[i]->transition[GW_INSERT2MATCH]+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)]);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state MATCH */ 
        temp = EstLoop3_HIDDEN_MATRIX(mat,i-1,j-3,DELETE) + (mat->query->seg[i]->transition[GW_DELETE2MATCH]+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)]);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = EstLoop3_HIDDEN_MATRIX(mat,i-1,j-2,MATCH) + mat->cp->trans[PCD_DELETE_1_BASE];    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = EstLoop3_HIDDEN_MATRIX(mat,i-1,j-1,MATCH) + mat->cp->trans[PCD_DELETE_2_BASE];    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = EstLoop3_HIDDEN_MATRIX(mat,i-1,j-4,MATCH) + mat->cp->trans[PCD_INSERT_1_BASE];    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = EstLoop3_HIDDEN_MATRIX(mat,i-1,j-5,MATCH) + mat->cp->trans[PCD_INSERT_2_BASE];    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state MATCH */ 
        temp = EstLoop3_HIDDEN_MATRIX(mat,i-1,j-2,INSERT) + mat->cp->trans[PCD_DELETE_1_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state MATCH */ 
        temp = EstLoop3_HIDDEN_MATRIX(mat,i-1,j-1,INSERT) + mat->cp->trans[PCD_DELETE_2_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state MATCH */ 
        temp = EstLoop3_HIDDEN_MATRIX(mat,i-1,j-4,INSERT) + mat->cp->trans[PCD_INSERT_1_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state MATCH */ 
        temp = EstLoop3_HIDDEN_MATRIX(mat,i-1,j-5,INSERT) + mat->cp->trans[PCD_INSERT_2_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state MATCH */ 
        temp = EstLoop3_HIDDEN_MATRIX(mat,i-1,j-2,DELETE) + mat->cp->trans[PCD_DELETE_1_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state MATCH */ 
        temp = EstLoop3_HIDDEN_MATRIX(mat,i-1,j-1,DELETE) + mat->cp->trans[PCD_DELETE_2_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state MATCH */ 
        temp = EstLoop3_HIDDEN_MATRIX(mat,i-1,j-4,DELETE) + mat->cp->trans[PCD_INSERT_1_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state MATCH */ 
        temp = EstLoop3_HIDDEN_MATRIX(mat,i-1,j-5,DELETE) + mat->cp->trans[PCD_INSERT_2_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for MATCH */ 
        /* Add any movement independant score and put away */ 
         EstLoop3_HIDDEN_MATRIX(mat,i,j,MATCH) = score;  
        /* Finished calculating state MATCH */ 


        /* For state INSERT */ 
        /* setting first movement to score */ 
        score = EstLoop3_HIDDEN_MATRIX(mat,i-0,j-3,MATCH) + (mat->query->seg[i]->transition[GW_MATCH2INSERT]+mat->query->seg[i]->insert[CSEQ_CDNA_CODON(mat->target,j)]);    
        /* From state INSERT to state INSERT */ 
        temp = EstLoop3_HIDDEN_MATRIX(mat,i-0,j-3,INSERT) + (mat->query->seg[i]->transition[GW_INSERT2INSERT]+mat->query->seg[i]->insert[CSEQ_CDNA_CODON(mat->target,j)]);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state INSERT */ 
        temp = EstLoop3_HIDDEN_MATRIX(mat,i-0,j-3,DELETE) + (mat->query->seg[i]->transition[GW_DELETE2INSERT]+mat->query->seg[i]->insert[CSEQ_CDNA_CODON(mat->target,j)]);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state INSERT */ 
        temp = EstLoop3_HIDDEN_MATRIX(mat,i-0,j-2,INSERT) + mat->cp->trans[PCD_DELETE_1_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state INSERT */ 
        temp = EstLoop3_HIDDEN_MATRIX(mat,i-0,j-1,INSERT) + mat->cp->trans[PCD_DELETE_2_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state INSERT */ 
        temp = EstLoop3_HIDDEN_MATRIX(mat,i-0,j-4,INSERT) + mat->cp->trans[PCD_INSERT_1_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state INSERT */ 
        temp = EstLoop3_HIDDEN_MATRIX(mat,i-0,j-5,INSERT) + mat->cp->trans[PCD_INSERT_2_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state INSERT */ 
        temp = EstLoop3_HIDDEN_MATRIX(mat,i-0,j-2,MATCH) + mat->cp->trans[PCD_DELETE_1_BASE];    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state INSERT */ 
        temp = EstLoop3_HIDDEN_MATRIX(mat,i-0,j-1,MATCH) + mat->cp->trans[PCD_DELETE_2_BASE];    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state INSERT */ 
        temp = EstLoop3_HIDDEN_MATRIX(mat,i-0,j-4,MATCH) + mat->cp->trans[PCD_INSERT_1_BASE];    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state INSERT */ 
        temp = EstLoop3_HIDDEN_MATRIX(mat,i-0,j-5,MATCH) + mat->cp->trans[PCD_INSERT_2_BASE];    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state INSERT */ 
        temp = EstLoop3_HIDDEN_MATRIX(mat,i-0,j-2,DELETE) + mat->cp->trans[PCD_DELETE_1_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state INSERT */ 
        temp = EstLoop3_HIDDEN_MATRIX(mat,i-0,j-1,DELETE) + mat->cp->trans[PCD_DELETE_2_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state INSERT */ 
        temp = EstLoop3_HIDDEN_MATRIX(mat,i-0,j-4,DELETE) + mat->cp->trans[PCD_INSERT_1_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state INSERT */ 
        temp = EstLoop3_HIDDEN_MATRIX(mat,i-0,j-5,DELETE) + mat->cp->trans[PCD_INSERT_2_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for INSERT */ 
        /* Add any movement independant score and put away */ 
         EstLoop3_HIDDEN_MATRIX(mat,i,j,INSERT) = score; 
        /* Finished calculating state INSERT */ 


        /* For state DELETE */ 
        /* setting first movement to score */ 
        score = EstLoop3_HIDDEN_MATRIX(mat,i-1,j-0,MATCH) + mat->query->seg[i]->transition[GW_MATCH2DELETE];     
        /* From state INSERT to state DELETE */ 
        temp = EstLoop3_HIDDEN_MATRIX(mat,i-1,j-0,INSERT) + mat->query->seg[i]->transition[GW_INSERT2DELETE];    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state DELETE */ 
        temp = EstLoop3_HIDDEN_MATRIX(mat,i-1,j-0,DELETE) + mat->query->seg[i]->transition[GW_DELETE2DELETE];    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for DELETE */ 
        /* Add any movement independant score and put away */ 
         EstLoop3_HIDDEN_MATRIX(mat,i,j,DELETE) = score; 
        /* Finished calculating state DELETE */ 
        }  
      }  


    return;  
}    


/* Function:  init_hidden_EstLoop3(mat,starti,startj,stopi,stopj)
 *
 * Descrip: No Description
 *
 * Arg:           mat [UNKN ] Undocumented argument [EstLoop3 *]
 * Arg:        starti [UNKN ] Undocumented argument [int]
 * Arg:        startj [UNKN ] Undocumented argument [int]
 * Arg:         stopi [UNKN ] Undocumented argument [int]
 * Arg:         stopj [UNKN ] Undocumented argument [int]
 *
 */
void init_hidden_EstLoop3(EstLoop3 * mat,int starti,int startj,int stopi,int stopj) 
{
    register int i;  
    register int j;  
    register int hiddenj;    


    hiddenj = startj;    
    for(j=(startj-5);j<=stopj;j++)   {  
      for(i=(starti-1);i<=stopi;i++) {  
        EstLoop3_HIDDEN_MATRIX(mat,i,j,MATCH) = NEGI;
   
        EstLoop3_HIDDEN_MATRIX(mat,i,j,INSERT) = NEGI;
  
        EstLoop3_HIDDEN_MATRIX(mat,i,j,DELETE) = NEGI;
  
        }  
      }  


    return;  
}    


/* Function:  full_dc_EstLoop3(mat,starti,startj,startstate,stopi,stopj,stopstate,out,donej,totalj,dpenv)
 *
 * Descrip:    The main divide-and-conquor routine. Basically, call /PackAln_calculate_small_EstLoop3
 *             Not this function, which is pretty hard core. 
 *             Function is given start/end points (in main matrix) for alignment
 *             It does some checks, decides whether start/end in j is small enough for explicit calc
 *               - if yes, calculates it, reads off into PackAln (out), adds the j distance to donej and returns TRUE
 *               - if no,  uses /do_dc_single_pass_EstLoop3 to get mid-point
 *                          saves midpoint, and calls itself to do right portion then left portion
 *             right then left ensures PackAln is added the 'right' way, ie, back-to-front
 *             returns FALSE on any error, with a warning
 *
 *
 * Arg:               mat [UNKN ] Matrix with small memory implementation [EstLoop3 *]
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
boolean full_dc_EstLoop3(EstLoop3 * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out,int * donej,int totalj,DPEnvelope * dpenv) 
{
    int lstarti; 
    int lstartj; 
    int lstate;  


    if( mat->basematrix->type != BASEMATRIX_TYPE_SHADOW) {  
      warn("*Very* bad error! - non shadow matrix type in full_dc_EstLoop3");    
      return FALSE;  
      }  


    if( starti == -1 || startj == -1 || startstate == -1 || stopi == -1 || stopstate == -1)  {  
      warn("In full dc program, passed bad indices, indices passed were %d:%d[%d] to %d:%d[%d]\n",starti,startj,startstate,stopi,stopj,stopstate);   
      return FALSE;  
      }  


    if( stopj - startj < 25) {  
      log_full_error(REPORT,0,"[%d,%d][%d,%d] Explicit read off",starti,startj,stopi,stopj);/* Build hidden explicit matrix */ 
      calculate_hidden_EstLoop3(mat,starti,startj,startstate,stopi,stopj,stopstate,dpenv);   
      *donej += (stopj - startj);   /* Now read it off into out */ 
      if( read_hidden_EstLoop3(mat,starti,startj,startstate,stopi,stopj,stopstate,out) == FALSE) {  
        warn("In full dc, at %d:%d,%d:%d got a bad hidden explicit read off... ",starti,startj,stopi,stopj); 
        return FALSE;    
        }  
      return TRUE;   
      }  


/* In actual divide and conquor */ 
    if( do_dc_single_pass_EstLoop3(mat,starti,startj,startstate,stopi,stopj,stopstate,dpenv,(int)(*donej*100)/totalj) == FALSE)  {  
      warn("In divide and conquor for EstLoop3, at bound %d:%d to %d:%d, unable to calculate midpoint. Problem!",starti,startj,stopi,stopj); 
      return FALSE;  
      }  


/* Ok... now we have to call on each side of the matrix */ 
/* We have to retrieve left hand side positions, as they will be vapped by the time we call LHS */ 
    lstarti= EstLoop3_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,0);  
    lstartj= EstLoop3_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,1);  
    lstate = EstLoop3_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,2);  


/* Call on right hand side: this lets us do the correct read off */ 
    if( full_dc_EstLoop3(mat,EstLoop3_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,3),EstLoop3_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,4),EstLoop3_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,5),stopi,stopj,stopstate,out,donej,totalj,dpenv) == FALSE)   {  
/* Warning already issued, simply chained back up to top */ 
      return FALSE;  
      }  
/* Call on left hand side */ 
    if( full_dc_EstLoop3(mat,starti,startj,startstate,lstarti,lstartj,lstate,out,donej,totalj,dpenv) == FALSE)   {  
/* Warning already issued, simply chained back up to top */ 
      return FALSE;  
      }  


    return TRUE;     
}    


/* Function:  do_dc_single_pass_EstLoop3(mat,starti,startj,startstate,stopi,stopj,stopstate,dpenv,perc_done)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [EstLoop3 *]
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
boolean do_dc_single_pass_EstLoop3(EstLoop3 * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv,int perc_done) 
{
    int halfj;   
    halfj = startj + ((stopj - startj)/2);   


    init_dc_EstLoop3(mat);   


    EstLoop3_DC_SHADOW_MATRIX(mat,starti,startj,startstate) = 0; 
    run_up_dc_EstLoop3(mat,starti,stopi,startj,halfj-1,dpenv,perc_done);     
    push_dc_at_merge_EstLoop3(mat,starti,stopi,halfj,&halfj,dpenv);  
    follow_on_dc_EstLoop3(mat,starti,stopi,halfj,stopj,dpenv,perc_done);     
    return TRUE; 
}    


/* Function:  push_dc_at_merge_EstLoop3(mat,starti,stopi,startj,stopj,dpenv)
 *
 * Descrip: No Description
 *
 * Arg:           mat [UNKN ] Undocumented argument [EstLoop3 *]
 * Arg:        starti [UNKN ] Undocumented argument [int]
 * Arg:         stopi [UNKN ] Undocumented argument [int]
 * Arg:        startj [UNKN ] Undocumented argument [int]
 * Arg:         stopj [UNKN ] Undocumented argument [int *]
 * Arg:         dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 */
void push_dc_at_merge_EstLoop3(EstLoop3 * mat,int starti,int stopi,int startj,int * stopj,DPEnvelope * dpenv) 
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
          EstLoop3_DC_SHADOW_MATRIX(mat,i,j,MATCH) = NEGI;   
          EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,0) = (-100);    
          EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,1) = (-100);    
          EstLoop3_DC_SHADOW_MATRIX(mat,i,j,INSERT) = NEGI;  
          EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,0) = (-100);   
          EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,1) = (-100);   
          EstLoop3_DC_SHADOW_MATRIX(mat,i,j,DELETE) = NEGI;  
          EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,0) = (-100);   
          EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,1) = (-100);   
          continue;  
          } /* end of Is not in envelope */ 


        /* For state MATCH, pushing when j - offj <= mergej */ 
        score = EstLoop3_DC_SHADOW_MATRIX(mat,i-1,j-3,MATCH) + (mat->query->seg[i]->transition[GW_MATCH2MATCH]+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)]);   
        if( j - 3 <= mergej) {  
          EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,0) = i-1;   
          EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,1) = j-3;   
          EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,2) = MATCH; 
          EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,3) = i; 
          EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,4) = j; 
          EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,5) = MATCH; 
          }  
        else {  
          for(k=0;k<7;k++)   
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = EstLoop3_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 3,MATCH,k);   
          }  


        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-1,j-3,INSERT) + (mat->query->seg[i]->transition[GW_INSERT2MATCH]+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)]);  
        if( temp > score)    {  
          score = temp;  


          if( j - 3 <= mergej)   {  
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,0) = i-1; 
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,1) = j-3; 
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,2) = INSERT;  
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,3) = i;   
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,4) = j;   
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,5) = MATCH;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = EstLoop3_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 3,INSERT,k);    
            }  
          }  


        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-1,j-3,DELETE) + (mat->query->seg[i]->transition[GW_DELETE2MATCH]+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)]);  
        if( temp > score)    {  
          score = temp;  


          if( j - 3 <= mergej)   {  
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,0) = i-1; 
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,1) = j-3; 
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,2) = DELETE;  
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,3) = i;   
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,4) = j;   
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,5) = MATCH;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = EstLoop3_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 3,DELETE,k);    
            }  
          }  


        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-1,j-2,MATCH) + mat->cp->trans[PCD_DELETE_1_BASE];     
        if( temp > score)    {  
          score = temp;  


          if( j - 2 <= mergej)   {  
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,0) = i-1; 
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,1) = j-2; 
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,2) = MATCH;   
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,3) = i;   
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,4) = j;   
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,5) = MATCH;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = EstLoop3_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 2,MATCH,k); 
            }  
          }  


        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-1,j-1,MATCH) + mat->cp->trans[PCD_DELETE_2_BASE];     
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,0) = i-1; 
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,1) = j-1; 
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,2) = MATCH;   
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,3) = i;   
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,4) = j;   
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,5) = MATCH;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = EstLoop3_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,MATCH,k); 
            }  
          }  


        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-1,j-4,MATCH) + mat->cp->trans[PCD_INSERT_1_BASE];     
        if( temp > score)    {  
          score = temp;  


          if( j - 4 <= mergej)   {  
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,0) = i-1; 
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,1) = j-4; 
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,2) = MATCH;   
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,3) = i;   
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,4) = j;   
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,5) = MATCH;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = EstLoop3_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 4,MATCH,k); 
            }  
          }  


        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-1,j-5,MATCH) + mat->cp->trans[PCD_INSERT_2_BASE];     
        if( temp > score)    {  
          score = temp;  


          if( j - 5 <= mergej)   {  
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,0) = i-1; 
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,1) = j-5; 
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,2) = MATCH;   
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,3) = i;   
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,4) = j;   
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,5) = MATCH;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = EstLoop3_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 5,MATCH,k); 
            }  
          }  


        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-1,j-2,INSERT) + mat->cp->trans[PCD_DELETE_1_BASE];    
        if( temp > score)    {  
          score = temp;  


          if( j - 2 <= mergej)   {  
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,0) = i-1; 
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,1) = j-2; 
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,2) = INSERT;  
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,3) = i;   
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,4) = j;   
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,5) = MATCH;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = EstLoop3_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 2,INSERT,k);    
            }  
          }  


        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-1,j-1,INSERT) + mat->cp->trans[PCD_DELETE_2_BASE];    
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,0) = i-1; 
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,1) = j-1; 
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,2) = INSERT;  
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,3) = i;   
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,4) = j;   
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,5) = MATCH;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = EstLoop3_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,INSERT,k);    
            }  
          }  


        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-1,j-4,INSERT) + mat->cp->trans[PCD_INSERT_1_BASE];    
        if( temp > score)    {  
          score = temp;  


          if( j - 4 <= mergej)   {  
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,0) = i-1; 
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,1) = j-4; 
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,2) = INSERT;  
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,3) = i;   
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,4) = j;   
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,5) = MATCH;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = EstLoop3_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 4,INSERT,k);    
            }  
          }  


        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-1,j-5,INSERT) + mat->cp->trans[PCD_INSERT_2_BASE];    
        if( temp > score)    {  
          score = temp;  


          if( j - 5 <= mergej)   {  
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,0) = i-1; 
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,1) = j-5; 
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,2) = INSERT;  
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,3) = i;   
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,4) = j;   
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,5) = MATCH;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = EstLoop3_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 5,INSERT,k);    
            }  
          }  


        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-1,j-2,DELETE) + mat->cp->trans[PCD_DELETE_1_BASE];    
        if( temp > score)    {  
          score = temp;  


          if( j - 2 <= mergej)   {  
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,0) = i-1; 
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,1) = j-2; 
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,2) = DELETE;  
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,3) = i;   
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,4) = j;   
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,5) = MATCH;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = EstLoop3_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 2,DELETE,k);    
            }  
          }  


        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-1,j-1,DELETE) + mat->cp->trans[PCD_DELETE_2_BASE];    
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,0) = i-1; 
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,1) = j-1; 
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,2) = DELETE;  
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,3) = i;   
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,4) = j;   
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,5) = MATCH;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = EstLoop3_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,DELETE,k);    
            }  
          }  


        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-1,j-4,DELETE) + mat->cp->trans[PCD_INSERT_1_BASE];    
        if( temp > score)    {  
          score = temp;  


          if( j - 4 <= mergej)   {  
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,0) = i-1; 
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,1) = j-4; 
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,2) = DELETE;  
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,3) = i;   
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,4) = j;   
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,5) = MATCH;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = EstLoop3_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 4,DELETE,k);    
            }  
          }  


        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-1,j-5,DELETE) + mat->cp->trans[PCD_INSERT_2_BASE];    
        if( temp > score)    {  
          score = temp;  


          if( j - 5 <= mergej)   {  
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,0) = i-1; 
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,1) = j-5; 
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,2) = DELETE;  
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,3) = i;   
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,4) = j;   
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,5) = MATCH;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = EstLoop3_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 5,DELETE,k);    
            }  
          }  
        /* Add any movement independant score */ 
        EstLoop3_DC_SHADOW_MATRIX(mat,i,j,MATCH) = score;    
        /* Finished with state MATCH */ 


        /* For state INSERT, pushing when j - offj <= mergej */ 
        score = EstLoop3_DC_SHADOW_MATRIX(mat,i-0,j-3,MATCH) + (mat->query->seg[i]->transition[GW_MATCH2INSERT]+mat->query->seg[i]->insert[CSEQ_CDNA_CODON(mat->target,j)]);     
        if( j - 3 <= mergej) {  
          EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,0) = i-0;  
          EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,1) = j-3;  
          EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,2) = MATCH;    
          EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,3) = i;    
          EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,4) = j;    
          EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,5) = INSERT;   
          }  
        else {  
          for(k=0;k<7;k++)   
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,k) = EstLoop3_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 3,MATCH,k);  
          }  


        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-0,j-3,INSERT) + (mat->query->seg[i]->transition[GW_INSERT2INSERT]+mat->query->seg[i]->insert[CSEQ_CDNA_CODON(mat->target,j)]);    
        if( temp > score)    {  
          score = temp;  


          if( j - 3 <= mergej)   {  
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,0) = i-0;    
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,1) = j-3;    
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,2) = INSERT; 
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,3) = i;  
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,4) = j;  
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,5) = INSERT; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,k) = EstLoop3_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 3,INSERT,k);   
            }  
          }  


        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-0,j-3,DELETE) + (mat->query->seg[i]->transition[GW_DELETE2INSERT]+mat->query->seg[i]->insert[CSEQ_CDNA_CODON(mat->target,j)]);    
        if( temp > score)    {  
          score = temp;  


          if( j - 3 <= mergej)   {  
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,0) = i-0;    
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,1) = j-3;    
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,2) = DELETE; 
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,3) = i;  
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,4) = j;  
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,5) = INSERT; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,k) = EstLoop3_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 3,DELETE,k);   
            }  
          }  


        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-0,j-2,INSERT) + mat->cp->trans[PCD_DELETE_1_BASE];    
        if( temp > score)    {  
          score = temp;  


          if( j - 2 <= mergej)   {  
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,0) = i-0;    
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,1) = j-2;    
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,2) = INSERT; 
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,3) = i;  
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,4) = j;  
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,5) = INSERT; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,k) = EstLoop3_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 2,INSERT,k);   
            }  
          }  


        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-0,j-1,INSERT) + mat->cp->trans[PCD_DELETE_2_BASE];    
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,0) = i-0;    
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,1) = j-1;    
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,2) = INSERT; 
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,3) = i;  
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,4) = j;  
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,5) = INSERT; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,k) = EstLoop3_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,INSERT,k);   
            }  
          }  


        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-0,j-4,INSERT) + mat->cp->trans[PCD_INSERT_1_BASE];    
        if( temp > score)    {  
          score = temp;  


          if( j - 4 <= mergej)   {  
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,0) = i-0;    
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,1) = j-4;    
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,2) = INSERT; 
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,3) = i;  
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,4) = j;  
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,5) = INSERT; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,k) = EstLoop3_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 4,INSERT,k);   
            }  
          }  


        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-0,j-5,INSERT) + mat->cp->trans[PCD_INSERT_2_BASE];    
        if( temp > score)    {  
          score = temp;  


          if( j - 5 <= mergej)   {  
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,0) = i-0;    
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,1) = j-5;    
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,2) = INSERT; 
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,3) = i;  
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,4) = j;  
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,5) = INSERT; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,k) = EstLoop3_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 5,INSERT,k);   
            }  
          }  


        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-0,j-2,MATCH) + mat->cp->trans[PCD_DELETE_1_BASE];     
        if( temp > score)    {  
          score = temp;  


          if( j - 2 <= mergej)   {  
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,0) = i-0;    
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,1) = j-2;    
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,2) = MATCH;  
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,3) = i;  
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,4) = j;  
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,5) = INSERT; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,k) = EstLoop3_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 2,MATCH,k);    
            }  
          }  


        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-0,j-1,MATCH) + mat->cp->trans[PCD_DELETE_2_BASE];     
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,0) = i-0;    
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,1) = j-1;    
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,2) = MATCH;  
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,3) = i;  
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,4) = j;  
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,5) = INSERT; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,k) = EstLoop3_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,MATCH,k);    
            }  
          }  


        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-0,j-4,MATCH) + mat->cp->trans[PCD_INSERT_1_BASE];     
        if( temp > score)    {  
          score = temp;  


          if( j - 4 <= mergej)   {  
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,0) = i-0;    
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,1) = j-4;    
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,2) = MATCH;  
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,3) = i;  
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,4) = j;  
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,5) = INSERT; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,k) = EstLoop3_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 4,MATCH,k);    
            }  
          }  


        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-0,j-5,MATCH) + mat->cp->trans[PCD_INSERT_2_BASE];     
        if( temp > score)    {  
          score = temp;  


          if( j - 5 <= mergej)   {  
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,0) = i-0;    
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,1) = j-5;    
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,2) = MATCH;  
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,3) = i;  
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,4) = j;  
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,5) = INSERT; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,k) = EstLoop3_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 5,MATCH,k);    
            }  
          }  


        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-0,j-2,DELETE) + mat->cp->trans[PCD_DELETE_1_BASE];    
        if( temp > score)    {  
          score = temp;  


          if( j - 2 <= mergej)   {  
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,0) = i-0;    
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,1) = j-2;    
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,2) = DELETE; 
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,3) = i;  
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,4) = j;  
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,5) = INSERT; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,k) = EstLoop3_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 2,DELETE,k);   
            }  
          }  


        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-0,j-1,DELETE) + mat->cp->trans[PCD_DELETE_2_BASE];    
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,0) = i-0;    
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,1) = j-1;    
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,2) = DELETE; 
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,3) = i;  
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,4) = j;  
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,5) = INSERT; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,k) = EstLoop3_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,DELETE,k);   
            }  
          }  


        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-0,j-4,DELETE) + mat->cp->trans[PCD_INSERT_1_BASE];    
        if( temp > score)    {  
          score = temp;  


          if( j - 4 <= mergej)   {  
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,0) = i-0;    
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,1) = j-4;    
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,2) = DELETE; 
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,3) = i;  
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,4) = j;  
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,5) = INSERT; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,k) = EstLoop3_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 4,DELETE,k);   
            }  
          }  


        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-0,j-5,DELETE) + mat->cp->trans[PCD_INSERT_2_BASE];    
        if( temp > score)    {  
          score = temp;  


          if( j - 5 <= mergej)   {  
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,0) = i-0;    
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,1) = j-5;    
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,2) = DELETE; 
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,3) = i;  
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,4) = j;  
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,5) = INSERT; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,k) = EstLoop3_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 5,DELETE,k);   
            }  
          }  
        /* Add any movement independant score */ 
        EstLoop3_DC_SHADOW_MATRIX(mat,i,j,INSERT) = score;   
        /* Finished with state INSERT */ 


        /* For state DELETE, pushing when j - offj <= mergej */ 
        score = EstLoop3_DC_SHADOW_MATRIX(mat,i-1,j-0,MATCH) + mat->query->seg[i]->transition[GW_MATCH2DELETE];  
        if( j - 0 <= mergej) {  
          EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,0) = i-1;  
          EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,1) = j-0;  
          EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,2) = MATCH;    
          EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,3) = i;    
          EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,4) = j;    
          EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,5) = DELETE;   
          }  
        else {  
          for(k=0;k<7;k++)   
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,k) = EstLoop3_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 0,MATCH,k);  
          }  


        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-1,j-0,INSERT) + mat->query->seg[i]->transition[GW_INSERT2DELETE];     
        if( temp > score)    {  
          score = temp;  


          if( j - 0 <= mergej)   {  
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,0) = i-1;    
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,1) = j-0;    
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,2) = INSERT; 
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,3) = i;  
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,4) = j;  
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,5) = DELETE; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,k) = EstLoop3_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 0,INSERT,k);   
            }  
          }  


        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-1,j-0,DELETE) + mat->query->seg[i]->transition[GW_DELETE2DELETE];     
        if( temp > score)    {  
          score = temp;  


          if( j - 0 <= mergej)   {  
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,0) = i-1;    
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,1) = j-0;    
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,2) = DELETE; 
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,3) = i;  
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,4) = j;  
            EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,5) = DELETE; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,k) = EstLoop3_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 0,DELETE,k);   
            }  
          }  
        /* Add any movement independant score */ 
        EstLoop3_DC_SHADOW_MATRIX(mat,i,j,DELETE) = score;   
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


/* Function:  follow_on_dc_EstLoop3(mat,starti,stopi,startj,stopj,dpenv,perc_done)
 *
 * Descrip: No Description
 *
 * Arg:              mat [UNKN ] Undocumented argument [EstLoop3 *]
 * Arg:           starti [UNKN ] Undocumented argument [int]
 * Arg:            stopi [UNKN ] Undocumented argument [int]
 * Arg:           startj [UNKN ] Undocumented argument [int]
 * Arg:            stopj [UNKN ] Undocumented argument [int]
 * Arg:            dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:        perc_done [UNKN ] Undocumented argument [int]
 *
 */
void follow_on_dc_EstLoop3(EstLoop3 * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done) 
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
          EstLoop3_DC_SHADOW_MATRIX(mat,i,j,MATCH) = NEGI;   
          EstLoop3_DC_SHADOW_MATRIX(mat,i,j,INSERT) = NEGI;  
          EstLoop3_DC_SHADOW_MATRIX(mat,i,j,DELETE) = NEGI;  
          continue;  
          } /* end of Is not in envelope */ 
        if( num % 1000 == 0 )    
          log_full_error(REPORT,0,"[%d%%%% done]After  mid-j %5d Cells done %d%%%%",perc_done,startj,(num*100)/total);   


        /* For state MATCH */ 
        /* setting first movement to score */ 
        score = EstLoop3_DC_SHADOW_MATRIX(mat,i-1,j-3,MATCH) + (mat->query->seg[i]->transition[GW_MATCH2MATCH]+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)]);   
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = EstLoop3_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 3,MATCH,k);    
        /* From state INSERT to state MATCH */ 
        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-1,j-3,INSERT) + (mat->query->seg[i]->transition[GW_INSERT2MATCH]+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)]);  
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = EstLoop3_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 3,INSERT,k); 
          }  
        /* From state DELETE to state MATCH */ 
        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-1,j-3,DELETE) + (mat->query->seg[i]->transition[GW_DELETE2MATCH]+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)]);  
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = EstLoop3_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 3,DELETE,k); 
          }  
        /* From state MATCH to state MATCH */ 
        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-1,j-2,MATCH) + mat->cp->trans[PCD_DELETE_1_BASE];     
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = EstLoop3_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 2,MATCH,k);  
          }  
        /* From state MATCH to state MATCH */ 
        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-1,j-1,MATCH) + mat->cp->trans[PCD_DELETE_2_BASE];     
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = EstLoop3_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,MATCH,k);  
          }  
        /* From state MATCH to state MATCH */ 
        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-1,j-4,MATCH) + mat->cp->trans[PCD_INSERT_1_BASE];     
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = EstLoop3_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 4,MATCH,k);  
          }  
        /* From state MATCH to state MATCH */ 
        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-1,j-5,MATCH) + mat->cp->trans[PCD_INSERT_2_BASE];     
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = EstLoop3_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 5,MATCH,k);  
          }  
        /* From state INSERT to state MATCH */ 
        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-1,j-2,INSERT) + mat->cp->trans[PCD_DELETE_1_BASE];    
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = EstLoop3_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 2,INSERT,k); 
          }  
        /* From state INSERT to state MATCH */ 
        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-1,j-1,INSERT) + mat->cp->trans[PCD_DELETE_2_BASE];    
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = EstLoop3_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,INSERT,k); 
          }  
        /* From state INSERT to state MATCH */ 
        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-1,j-4,INSERT) + mat->cp->trans[PCD_INSERT_1_BASE];    
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = EstLoop3_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 4,INSERT,k); 
          }  
        /* From state INSERT to state MATCH */ 
        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-1,j-5,INSERT) + mat->cp->trans[PCD_INSERT_2_BASE];    
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = EstLoop3_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 5,INSERT,k); 
          }  
        /* From state DELETE to state MATCH */ 
        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-1,j-2,DELETE) + mat->cp->trans[PCD_DELETE_1_BASE];    
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = EstLoop3_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 2,DELETE,k); 
          }  
        /* From state DELETE to state MATCH */ 
        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-1,j-1,DELETE) + mat->cp->trans[PCD_DELETE_2_BASE];    
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = EstLoop3_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,DELETE,k); 
          }  
        /* From state DELETE to state MATCH */ 
        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-1,j-4,DELETE) + mat->cp->trans[PCD_INSERT_1_BASE];    
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = EstLoop3_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 4,DELETE,k); 
          }  
        /* From state DELETE to state MATCH */ 
        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-1,j-5,DELETE) + mat->cp->trans[PCD_INSERT_2_BASE];    
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = EstLoop3_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 5,DELETE,k); 
          }  


        /* Ok - finished max calculation for MATCH */ 
        /* Add any movement independant score and put away */ 
         EstLoop3_DC_SHADOW_MATRIX(mat,i,j,MATCH) = score;   
        for(k=0;k<7;k++) 
          EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = localshadow[k];    
        /* Now figure out if any specials need this score */ 
        /* Finished calculating state MATCH */ 


        /* For state INSERT */ 
        /* setting first movement to score */ 
        score = EstLoop3_DC_SHADOW_MATRIX(mat,i-0,j-3,MATCH) + (mat->query->seg[i]->transition[GW_MATCH2INSERT]+mat->query->seg[i]->insert[CSEQ_CDNA_CODON(mat->target,j)]);     
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = EstLoop3_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 3,MATCH,k);    
        /* From state INSERT to state INSERT */ 
        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-0,j-3,INSERT) + (mat->query->seg[i]->transition[GW_INSERT2INSERT]+mat->query->seg[i]->insert[CSEQ_CDNA_CODON(mat->target,j)]);    
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = EstLoop3_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 3,INSERT,k); 
          }  
        /* From state DELETE to state INSERT */ 
        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-0,j-3,DELETE) + (mat->query->seg[i]->transition[GW_DELETE2INSERT]+mat->query->seg[i]->insert[CSEQ_CDNA_CODON(mat->target,j)]);    
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = EstLoop3_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 3,DELETE,k); 
          }  
        /* From state INSERT to state INSERT */ 
        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-0,j-2,INSERT) + mat->cp->trans[PCD_DELETE_1_BASE];    
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = EstLoop3_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 2,INSERT,k); 
          }  
        /* From state INSERT to state INSERT */ 
        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-0,j-1,INSERT) + mat->cp->trans[PCD_DELETE_2_BASE];    
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = EstLoop3_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,INSERT,k); 
          }  
        /* From state INSERT to state INSERT */ 
        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-0,j-4,INSERT) + mat->cp->trans[PCD_INSERT_1_BASE];    
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = EstLoop3_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 4,INSERT,k); 
          }  
        /* From state INSERT to state INSERT */ 
        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-0,j-5,INSERT) + mat->cp->trans[PCD_INSERT_2_BASE];    
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = EstLoop3_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 5,INSERT,k); 
          }  
        /* From state MATCH to state INSERT */ 
        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-0,j-2,MATCH) + mat->cp->trans[PCD_DELETE_1_BASE];     
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = EstLoop3_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 2,MATCH,k);  
          }  
        /* From state MATCH to state INSERT */ 
        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-0,j-1,MATCH) + mat->cp->trans[PCD_DELETE_2_BASE];     
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = EstLoop3_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,MATCH,k);  
          }  
        /* From state MATCH to state INSERT */ 
        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-0,j-4,MATCH) + mat->cp->trans[PCD_INSERT_1_BASE];     
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = EstLoop3_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 4,MATCH,k);  
          }  
        /* From state MATCH to state INSERT */ 
        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-0,j-5,MATCH) + mat->cp->trans[PCD_INSERT_2_BASE];     
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = EstLoop3_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 5,MATCH,k);  
          }  
        /* From state DELETE to state INSERT */ 
        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-0,j-2,DELETE) + mat->cp->trans[PCD_DELETE_1_BASE];    
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = EstLoop3_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 2,DELETE,k); 
          }  
        /* From state DELETE to state INSERT */ 
        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-0,j-1,DELETE) + mat->cp->trans[PCD_DELETE_2_BASE];    
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = EstLoop3_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,DELETE,k); 
          }  
        /* From state DELETE to state INSERT */ 
        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-0,j-4,DELETE) + mat->cp->trans[PCD_INSERT_1_BASE];    
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = EstLoop3_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 4,DELETE,k); 
          }  
        /* From state DELETE to state INSERT */ 
        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-0,j-5,DELETE) + mat->cp->trans[PCD_INSERT_2_BASE];    
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = EstLoop3_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 5,DELETE,k); 
          }  


        /* Ok - finished max calculation for INSERT */ 
        /* Add any movement independant score and put away */ 
         EstLoop3_DC_SHADOW_MATRIX(mat,i,j,INSERT) = score;  
        for(k=0;k<7;k++) 
          EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,k) = localshadow[k];   
        /* Now figure out if any specials need this score */ 
        /* Finished calculating state INSERT */ 


        /* For state DELETE */ 
        /* setting first movement to score */ 
        score = EstLoop3_DC_SHADOW_MATRIX(mat,i-1,j-0,MATCH) + mat->query->seg[i]->transition[GW_MATCH2DELETE];  
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = EstLoop3_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 0,MATCH,k);    
        /* From state INSERT to state DELETE */ 
        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-1,j-0,INSERT) + mat->query->seg[i]->transition[GW_INSERT2DELETE];     
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = EstLoop3_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 0,INSERT,k); 
          }  
        /* From state DELETE to state DELETE */ 
        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-1,j-0,DELETE) + mat->query->seg[i]->transition[GW_DELETE2DELETE];     
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = EstLoop3_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 0,DELETE,k); 
          }  


        /* Ok - finished max calculation for DELETE */ 
        /* Add any movement independant score and put away */ 
         EstLoop3_DC_SHADOW_MATRIX(mat,i,j,DELETE) = score;  
        for(k=0;k<7;k++) 
          EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,k) = localshadow[k];   
        /* Now figure out if any specials need this score */ 
        /* Finished calculating state DELETE */ 
        } /* end of this is strip */ 
      } /* end of for each valid j column */ 


/* Function:  run_up_dc_EstLoop3(mat,starti,stopi,startj,stopj,dpenv,perc_done)
 *
 * Descrip: No Description
 *
 * Arg:              mat [UNKN ] Undocumented argument [EstLoop3 *]
 * Arg:           starti [UNKN ] Undocumented argument [int]
 * Arg:            stopi [UNKN ] Undocumented argument [int]
 * Arg:           startj [UNKN ] Undocumented argument [int]
 * Arg:            stopj [UNKN ] Undocumented argument [int]
 * Arg:            dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:        perc_done [UNKN ] Undocumented argument [int]
 *
 */
}    
void run_up_dc_EstLoop3(EstLoop3 * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done) 
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
          EstLoop3_DC_SHADOW_MATRIX(mat,i,j,MATCH) = NEGI;   
          EstLoop3_DC_SHADOW_MATRIX(mat,i,j,INSERT) = NEGI;  
          EstLoop3_DC_SHADOW_MATRIX(mat,i,j,DELETE) = NEGI;  
          continue;  
          } /* end of Is not in envelope */ 
        if( num % 1000 == 0 )    
          log_full_error(REPORT,0,"[%d%%%% done]Before mid-j %5d Cells done %d%%%%",perc_done,stopj,(num*100)/total);    


        /* For state MATCH */ 
        /* setting first movement to score */ 
        score = EstLoop3_DC_SHADOW_MATRIX(mat,i-1,j-3,MATCH) + (mat->query->seg[i]->transition[GW_MATCH2MATCH]+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)]);   
        /* From state INSERT to state MATCH */ 
        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-1,j-3,INSERT) + (mat->query->seg[i]->transition[GW_INSERT2MATCH]+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)]);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state MATCH */ 
        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-1,j-3,DELETE) + (mat->query->seg[i]->transition[GW_DELETE2MATCH]+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)]);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-1,j-2,MATCH) + mat->cp->trans[PCD_DELETE_1_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-1,j-1,MATCH) + mat->cp->trans[PCD_DELETE_2_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-1,j-4,MATCH) + mat->cp->trans[PCD_INSERT_1_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-1,j-5,MATCH) + mat->cp->trans[PCD_INSERT_2_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state MATCH */ 
        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-1,j-2,INSERT) + mat->cp->trans[PCD_DELETE_1_BASE];    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state MATCH */ 
        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-1,j-1,INSERT) + mat->cp->trans[PCD_DELETE_2_BASE];    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state MATCH */ 
        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-1,j-4,INSERT) + mat->cp->trans[PCD_INSERT_1_BASE];    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state MATCH */ 
        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-1,j-5,INSERT) + mat->cp->trans[PCD_INSERT_2_BASE];    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state MATCH */ 
        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-1,j-2,DELETE) + mat->cp->trans[PCD_DELETE_1_BASE];    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state MATCH */ 
        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-1,j-1,DELETE) + mat->cp->trans[PCD_DELETE_2_BASE];    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state MATCH */ 
        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-1,j-4,DELETE) + mat->cp->trans[PCD_INSERT_1_BASE];    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state MATCH */ 
        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-1,j-5,DELETE) + mat->cp->trans[PCD_INSERT_2_BASE];    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for MATCH */ 
        /* Add any movement independant score and put away */ 
         EstLoop3_DC_SHADOW_MATRIX(mat,i,j,MATCH) = score;   
        /* Finished calculating state MATCH */ 


        /* For state INSERT */ 
        /* setting first movement to score */ 
        score = EstLoop3_DC_SHADOW_MATRIX(mat,i-0,j-3,MATCH) + (mat->query->seg[i]->transition[GW_MATCH2INSERT]+mat->query->seg[i]->insert[CSEQ_CDNA_CODON(mat->target,j)]);     
        /* From state INSERT to state INSERT */ 
        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-0,j-3,INSERT) + (mat->query->seg[i]->transition[GW_INSERT2INSERT]+mat->query->seg[i]->insert[CSEQ_CDNA_CODON(mat->target,j)]);    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state INSERT */ 
        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-0,j-3,DELETE) + (mat->query->seg[i]->transition[GW_DELETE2INSERT]+mat->query->seg[i]->insert[CSEQ_CDNA_CODON(mat->target,j)]);    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state INSERT */ 
        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-0,j-2,INSERT) + mat->cp->trans[PCD_DELETE_1_BASE];    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state INSERT */ 
        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-0,j-1,INSERT) + mat->cp->trans[PCD_DELETE_2_BASE];    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state INSERT */ 
        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-0,j-4,INSERT) + mat->cp->trans[PCD_INSERT_1_BASE];    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state INSERT */ 
        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-0,j-5,INSERT) + mat->cp->trans[PCD_INSERT_2_BASE];    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state INSERT */ 
        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-0,j-2,MATCH) + mat->cp->trans[PCD_DELETE_1_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state INSERT */ 
        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-0,j-1,MATCH) + mat->cp->trans[PCD_DELETE_2_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state INSERT */ 
        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-0,j-4,MATCH) + mat->cp->trans[PCD_INSERT_1_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state INSERT */ 
        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-0,j-5,MATCH) + mat->cp->trans[PCD_INSERT_2_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state INSERT */ 
        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-0,j-2,DELETE) + mat->cp->trans[PCD_DELETE_1_BASE];    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state INSERT */ 
        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-0,j-1,DELETE) + mat->cp->trans[PCD_DELETE_2_BASE];    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state INSERT */ 
        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-0,j-4,DELETE) + mat->cp->trans[PCD_INSERT_1_BASE];    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state INSERT */ 
        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-0,j-5,DELETE) + mat->cp->trans[PCD_INSERT_2_BASE];    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for INSERT */ 
        /* Add any movement independant score and put away */ 
         EstLoop3_DC_SHADOW_MATRIX(mat,i,j,INSERT) = score;  
        /* Finished calculating state INSERT */ 


        /* For state DELETE */ 
        /* setting first movement to score */ 
        score = EstLoop3_DC_SHADOW_MATRIX(mat,i-1,j-0,MATCH) + mat->query->seg[i]->transition[GW_MATCH2DELETE];  
        /* From state INSERT to state DELETE */ 
        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-1,j-0,INSERT) + mat->query->seg[i]->transition[GW_INSERT2DELETE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state DELETE */ 
        temp = EstLoop3_DC_SHADOW_MATRIX(mat,i-1,j-0,DELETE) + mat->query->seg[i]->transition[GW_DELETE2DELETE];     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for DELETE */ 
        /* Add any movement independant score and put away */ 
         EstLoop3_DC_SHADOW_MATRIX(mat,i,j,DELETE) = score;  
        /* Finished calculating state DELETE */ 
        } /* end of this is strip */ 
      } /* end of for each valid j column */ 


/* Function:  init_dc_EstLoop3(mat)
 *
 * Descrip: No Description
 *
 * Arg:        mat [UNKN ] Undocumented argument [EstLoop3 *]
 *
 */
}    
void init_dc_EstLoop3(EstLoop3 * mat) 
{
    register int i;  
    register int j;  
    register int k;  


    for(j=0;j<7;j++) {  
      for(i=(-1);i<mat->query->len;i++)  {  
        EstLoop3_DC_SHADOW_MATRIX(mat,i,j,MATCH) = NEGI; 
        EstLoop3_DC_SHADOW_MATRIX(mat,i,j,INSERT) = NEGI;    
        EstLoop3_DC_SHADOW_MATRIX(mat,i,j,DELETE) = NEGI;    
        for(k=0;k<7;k++) {  
          EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = (-1);  
          EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,k) = (-1); 
          EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,k) = (-1); 
          }  
        }  
      }  


    return;  
}    


/* Function:  start_end_find_end_EstLoop3(mat,endj)
 *
 * Descrip:    First function used to find end of the best path in the special state !end
 *
 *
 * Arg:         mat [UNKN ] Matrix in small mode [EstLoop3 *]
 * Arg:        endj [WRITE] position of end in j (meaningless in i) [int *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int start_end_find_end_EstLoop3(EstLoop3 * mat,int * endj) 
{
    register int j;  
    register int max;    
    register int maxj;   


    max = EstLoop3_DC_SHADOW_SPECIAL(mat,0,mat->target->seq->len-1,END); 
    maxj = mat->target->seq->len-1;  
    for(j= mat->target->seq->len-2 ;j >= 0 ;j--) {  
      if( EstLoop3_DC_SHADOW_SPECIAL(mat,0,j,END) > max )    {  
        max = EstLoop3_DC_SHADOW_SPECIAL(mat,0,j,END);   
        maxj = j;    
        }  
      }  


    if( endj != NULL)    
      *endj = maxj;  


    return max;  
}    


/* Function:  dc_optimised_start_end_calc_EstLoop3(*mat,dpenv)
 *
 * Descrip:    Calculates special strip, leaving start/end/score points in shadow matrix
 *             Works off specially laid out memory from steve searle
 *
 *
 * Arg:         *mat [UNKN ] Undocumented argument [EstLoop3]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean dc_optimised_start_end_calc_EstLoop3(EstLoop3 *mat,DPEnvelope * dpenv) 
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
          EstLoop3_DC_OPT_SHADOW_MATRIX(mat,i,j,MATCH) = NEGI;   
          EstLoop3_DC_OPT_SHADOW_MATRIX(mat,i,j,INSERT) = NEGI;  
          EstLoop3_DC_OPT_SHADOW_MATRIX(mat,i,j,DELETE) = NEGI;  
          continue;  
          } /* end of Is not in envelope */ 
        if( num%1000 == 0)   
          log_full_error(REPORT,0,"%6d Cells done [%2d%%%%]",num,num*100/total); 




        /* For state MATCH */ 
        /* setting first movement to score */ 
        score = EstLoop3_DC_OPT_SHADOW_MATRIX(mat,i-1,j-3,MATCH) + (mat->query->seg[i]->transition[GW_MATCH2MATCH]+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)]) + (0);     
        /* assign local shadown pointer */ 
        localsp = &(EstLoop3_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 3,MATCH,0));  
        /* From state INSERT to state MATCH */ 
        temp = EstLoop3_DC_OPT_SHADOW_MATRIX(mat,i-1,j-3,INSERT) + (mat->query->seg[i]->transition[GW_INSERT2MATCH]+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)]) +(0);     
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(EstLoop3_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 3,INSERT,0));   
          }  
        /* From state DELETE to state MATCH */ 
        temp = EstLoop3_DC_OPT_SHADOW_MATRIX(mat,i-1,j-3,DELETE) + (mat->query->seg[i]->transition[GW_DELETE2MATCH]+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)]) +(0);     
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(EstLoop3_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 3,DELETE,0));   
          }  
        /* From state LOOP to state MATCH */ 
        temp = EstLoop3_DC_OPT_SHADOW_SPECIAL(mat,i-1,j-3,LOOP) + (mat->query->seg[i]->transition[GW_START2MATCH]+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)]) + (0);  
        if( temp  > score )  {  
          score = temp;  
          /* This state [LOOP] is a special for MATCH... push top shadow pointers here */ 
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
        temp = EstLoop3_DC_OPT_SHADOW_SPECIAL(mat,i-1,j-3,START) + (mat->query->seg[i]->transition[GW_START2MATCH]+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)]) + (0);     
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
        temp = EstLoop3_DC_OPT_SHADOW_MATRIX(mat,i-1,j-2,MATCH) + mat->cp->trans[PCD_DELETE_1_BASE] +(0);    
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(EstLoop3_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 2,MATCH,0));    
          }  
        /* From state MATCH to state MATCH */ 
        temp = EstLoop3_DC_OPT_SHADOW_MATRIX(mat,i-1,j-1,MATCH) + mat->cp->trans[PCD_DELETE_2_BASE] +(0);    
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(EstLoop3_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 1,MATCH,0));    
          }  
        /* From state MATCH to state MATCH */ 
        temp = EstLoop3_DC_OPT_SHADOW_MATRIX(mat,i-1,j-4,MATCH) + mat->cp->trans[PCD_INSERT_1_BASE] +(0);    
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(EstLoop3_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 4,MATCH,0));    
          }  
        /* From state MATCH to state MATCH */ 
        temp = EstLoop3_DC_OPT_SHADOW_MATRIX(mat,i-1,j-5,MATCH) + mat->cp->trans[PCD_INSERT_2_BASE] +(0);    
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(EstLoop3_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 5,MATCH,0));    
          }  
        /* From state INSERT to state MATCH */ 
        temp = EstLoop3_DC_OPT_SHADOW_MATRIX(mat,i-1,j-2,INSERT) + mat->cp->trans[PCD_DELETE_1_BASE] +(0);   
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(EstLoop3_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 2,INSERT,0));   
          }  
        /* From state INSERT to state MATCH */ 
        temp = EstLoop3_DC_OPT_SHADOW_MATRIX(mat,i-1,j-1,INSERT) + mat->cp->trans[PCD_DELETE_2_BASE] +(0);   
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(EstLoop3_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 1,INSERT,0));   
          }  
        /* From state INSERT to state MATCH */ 
        temp = EstLoop3_DC_OPT_SHADOW_MATRIX(mat,i-1,j-4,INSERT) + mat->cp->trans[PCD_INSERT_1_BASE] +(0);   
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(EstLoop3_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 4,INSERT,0));   
          }  
        /* From state INSERT to state MATCH */ 
        temp = EstLoop3_DC_OPT_SHADOW_MATRIX(mat,i-1,j-5,INSERT) + mat->cp->trans[PCD_INSERT_2_BASE] +(0);   
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(EstLoop3_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 5,INSERT,0));   
          }  
        /* From state DELETE to state MATCH */ 
        temp = EstLoop3_DC_OPT_SHADOW_MATRIX(mat,i-1,j-2,DELETE) + mat->cp->trans[PCD_DELETE_1_BASE] +(0);   
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(EstLoop3_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 2,DELETE,0));   
          }  
        /* From state DELETE to state MATCH */ 
        temp = EstLoop3_DC_OPT_SHADOW_MATRIX(mat,i-1,j-1,DELETE) + mat->cp->trans[PCD_DELETE_2_BASE] +(0);   
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(EstLoop3_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 1,DELETE,0));   
          }  
        /* From state DELETE to state MATCH */ 
        temp = EstLoop3_DC_OPT_SHADOW_MATRIX(mat,i-1,j-4,DELETE) + mat->cp->trans[PCD_INSERT_1_BASE] +(0);   
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(EstLoop3_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 4,DELETE,0));   
          }  
        /* From state DELETE to state MATCH */ 
        temp = EstLoop3_DC_OPT_SHADOW_MATRIX(mat,i-1,j-5,DELETE) + mat->cp->trans[PCD_INSERT_2_BASE] +(0);   
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(EstLoop3_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 5,DELETE,0));   
          }  


        /* Ok - finished max calculation for MATCH */ 
        /* Add any movement independant score and put away */ 
        /* Actually, already done inside scores */ 
         EstLoop3_DC_OPT_SHADOW_MATRIX(mat,i,j,MATCH) = score;   
        for(k=0;k<7;k++) 
          EstLoop3_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = localsp[k];    
        /* Now figure out if any specials need this score */ 


        /* state MATCH is a source for special LOOP */ 
        temp = score + (mat->query->seg[i]->transition[GW_MATCH2END]) + (0) ;    
        if( temp > EstLoop3_DC_OPT_SHADOW_SPECIAL(mat,i,j,LOOP) )    {  
          EstLoop3_DC_OPT_SHADOW_SPECIAL(mat,i,j,LOOP) = temp;   
          /* Have to push only bottem half of system here */ 
          for(k=0;k<3;k++)   
            EstLoop3_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,LOOP,k) = EstLoop3_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,MATCH,k);   
          EstLoop3_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,LOOP,6) = EstLoop3_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,MATCH,6); 
          EstLoop3_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,LOOP,3) = i; 
          EstLoop3_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,LOOP,4) = j; 
          EstLoop3_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,LOOP,5) = MATCH; 
          }  




        /* state MATCH is a source for special END */ 
        temp = score + (mat->query->seg[i]->transition[GW_MATCH2END]) + (0) ;    
        if( temp > EstLoop3_DC_OPT_SHADOW_SPECIAL(mat,i,j,END) )     {  
          EstLoop3_DC_OPT_SHADOW_SPECIAL(mat,i,j,END) = temp;    
          /* Have to push only bottem half of system here */ 
          for(k=0;k<3;k++)   
            EstLoop3_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,k) = EstLoop3_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,MATCH,k);    
          EstLoop3_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,6) = EstLoop3_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,MATCH,6);  
          EstLoop3_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,3) = i;  
          EstLoop3_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,4) = j;  
          EstLoop3_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,5) = MATCH;  
          }  




        /* Finished calculating state MATCH */ 


        /* For state INSERT */ 
        /* setting first movement to score */ 
        score = EstLoop3_DC_OPT_SHADOW_MATRIX(mat,i-0,j-3,MATCH) + (mat->query->seg[i]->transition[GW_MATCH2INSERT]+mat->query->seg[i]->insert[CSEQ_CDNA_CODON(mat->target,j)]) + (0);   
        /* assign local shadown pointer */ 
        localsp = &(EstLoop3_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 3,MATCH,0));  
        /* From state INSERT to state INSERT */ 
        temp = EstLoop3_DC_OPT_SHADOW_MATRIX(mat,i-0,j-3,INSERT) + (mat->query->seg[i]->transition[GW_INSERT2INSERT]+mat->query->seg[i]->insert[CSEQ_CDNA_CODON(mat->target,j)]) +(0);   
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(EstLoop3_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 3,INSERT,0));   
          }  
        /* From state DELETE to state INSERT */ 
        temp = EstLoop3_DC_OPT_SHADOW_MATRIX(mat,i-0,j-3,DELETE) + (mat->query->seg[i]->transition[GW_DELETE2INSERT]+mat->query->seg[i]->insert[CSEQ_CDNA_CODON(mat->target,j)]) +(0);   
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(EstLoop3_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 3,DELETE,0));   
          }  
        /* From state LOOP to state INSERT */ 
        temp = EstLoop3_DC_OPT_SHADOW_SPECIAL(mat,i-0,j-3,LOOP) + (mat->query->seg[i]->transition[GW_START2INSERT]+mat->query->seg[i]->insert[CSEQ_CDNA_CODON(mat->target,j)]) + (0);    
        if( temp  > score )  {  
          score = temp;  
          /* This state [LOOP] is a special for INSERT... push top shadow pointers here */ 
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
        temp = EstLoop3_DC_OPT_SHADOW_MATRIX(mat,i-0,j-2,INSERT) + mat->cp->trans[PCD_DELETE_1_BASE] +(0);   
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(EstLoop3_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 2,INSERT,0));   
          }  
        /* From state INSERT to state INSERT */ 
        temp = EstLoop3_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,INSERT) + mat->cp->trans[PCD_DELETE_2_BASE] +(0);   
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(EstLoop3_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,INSERT,0));   
          }  
        /* From state INSERT to state INSERT */ 
        temp = EstLoop3_DC_OPT_SHADOW_MATRIX(mat,i-0,j-4,INSERT) + mat->cp->trans[PCD_INSERT_1_BASE] +(0);   
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(EstLoop3_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 4,INSERT,0));   
          }  
        /* From state INSERT to state INSERT */ 
        temp = EstLoop3_DC_OPT_SHADOW_MATRIX(mat,i-0,j-5,INSERT) + mat->cp->trans[PCD_INSERT_2_BASE] +(0);   
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(EstLoop3_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 5,INSERT,0));   
          }  
        /* From state MATCH to state INSERT */ 
        temp = EstLoop3_DC_OPT_SHADOW_MATRIX(mat,i-0,j-2,MATCH) + mat->cp->trans[PCD_DELETE_1_BASE] +(0);    
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(EstLoop3_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 2,MATCH,0));    
          }  
        /* From state MATCH to state INSERT */ 
        temp = EstLoop3_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,MATCH) + mat->cp->trans[PCD_DELETE_2_BASE] +(0);    
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(EstLoop3_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,MATCH,0));    
          }  
        /* From state MATCH to state INSERT */ 
        temp = EstLoop3_DC_OPT_SHADOW_MATRIX(mat,i-0,j-4,MATCH) + mat->cp->trans[PCD_INSERT_1_BASE] +(0);    
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(EstLoop3_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 4,MATCH,0));    
          }  
        /* From state MATCH to state INSERT */ 
        temp = EstLoop3_DC_OPT_SHADOW_MATRIX(mat,i-0,j-5,MATCH) + mat->cp->trans[PCD_INSERT_2_BASE] +(0);    
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(EstLoop3_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 5,MATCH,0));    
          }  
        /* From state DELETE to state INSERT */ 
        temp = EstLoop3_DC_OPT_SHADOW_MATRIX(mat,i-0,j-2,DELETE) + mat->cp->trans[PCD_DELETE_1_BASE] +(0);   
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(EstLoop3_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 2,DELETE,0));   
          }  
        /* From state DELETE to state INSERT */ 
        temp = EstLoop3_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,DELETE) + mat->cp->trans[PCD_DELETE_2_BASE] +(0);   
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(EstLoop3_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,DELETE,0));   
          }  
        /* From state DELETE to state INSERT */ 
        temp = EstLoop3_DC_OPT_SHADOW_MATRIX(mat,i-0,j-4,DELETE) + mat->cp->trans[PCD_INSERT_1_BASE] +(0);   
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(EstLoop3_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 4,DELETE,0));   
          }  
        /* From state DELETE to state INSERT */ 
        temp = EstLoop3_DC_OPT_SHADOW_MATRIX(mat,i-0,j-5,DELETE) + mat->cp->trans[PCD_INSERT_2_BASE] +(0);   
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(EstLoop3_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 5,DELETE,0));   
          }  


        /* Ok - finished max calculation for INSERT */ 
        /* Add any movement independant score and put away */ 
        /* Actually, already done inside scores */ 
         EstLoop3_DC_OPT_SHADOW_MATRIX(mat,i,j,INSERT) = score;  
        for(k=0;k<7;k++) 
          EstLoop3_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,INSERT,k) = localsp[k];   
        /* Now figure out if any specials need this score */ 


        /* state INSERT is a source for special LOOP */ 
        temp = score + (mat->query->seg[i]->transition[GW_INSERT2END]) + (0) ;   
        if( temp > EstLoop3_DC_OPT_SHADOW_SPECIAL(mat,i,j,LOOP) )    {  
          EstLoop3_DC_OPT_SHADOW_SPECIAL(mat,i,j,LOOP) = temp;   
          /* Have to push only bottem half of system here */ 
          for(k=0;k<3;k++)   
            EstLoop3_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,LOOP,k) = EstLoop3_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,INSERT,k);  
          EstLoop3_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,LOOP,6) = EstLoop3_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,INSERT,6);    
          EstLoop3_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,LOOP,3) = i; 
          EstLoop3_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,LOOP,4) = j; 
          EstLoop3_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,LOOP,5) = INSERT;    
          }  




        /* state INSERT is a source for special END */ 
        temp = score + (mat->query->seg[i]->transition[GW_INSERT2END]) + (0) ;   
        if( temp > EstLoop3_DC_OPT_SHADOW_SPECIAL(mat,i,j,END) )     {  
          EstLoop3_DC_OPT_SHADOW_SPECIAL(mat,i,j,END) = temp;    
          /* Have to push only bottem half of system here */ 
          for(k=0;k<3;k++)   
            EstLoop3_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,k) = EstLoop3_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,INSERT,k);   
          EstLoop3_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,6) = EstLoop3_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,INSERT,6); 
          EstLoop3_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,3) = i;  
          EstLoop3_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,4) = j;  
          EstLoop3_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,5) = INSERT; 
          }  




        /* Finished calculating state INSERT */ 


        /* For state DELETE */ 
        /* setting first movement to score */ 
        score = EstLoop3_DC_OPT_SHADOW_MATRIX(mat,i-1,j-0,MATCH) + mat->query->seg[i]->transition[GW_MATCH2DELETE] + (0);    
        /* assign local shadown pointer */ 
        localsp = &(EstLoop3_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 0,MATCH,0));  
        /* From state INSERT to state DELETE */ 
        temp = EstLoop3_DC_OPT_SHADOW_MATRIX(mat,i-1,j-0,INSERT) + mat->query->seg[i]->transition[GW_INSERT2DELETE] +(0);    
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(EstLoop3_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 0,INSERT,0));   
          }  
        /* From state DELETE to state DELETE */ 
        temp = EstLoop3_DC_OPT_SHADOW_MATRIX(mat,i-1,j-0,DELETE) + mat->query->seg[i]->transition[GW_DELETE2DELETE] +(0);    
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(EstLoop3_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 0,DELETE,0));   
          }  
        /* From state LOOP to state DELETE */ 
        temp = EstLoop3_DC_OPT_SHADOW_SPECIAL(mat,i-1,j-0,LOOP) + mat->query->seg[i]->transition[GW_START2DELETE] + (0);     
        if( temp  > score )  {  
          score = temp;  
          /* This state [LOOP] is a special for DELETE... push top shadow pointers here */ 
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
         EstLoop3_DC_OPT_SHADOW_MATRIX(mat,i,j,DELETE) = score;  
        for(k=0;k<7;k++) 
          EstLoop3_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,DELETE,k) = localsp[k];   
        /* Now figure out if any specials need this score */ 


        /* state DELETE is a source for special LOOP */ 
        temp = score + (mat->query->seg[i]->transition[GW_DELETE2END]) + (0) ;   
        if( temp > EstLoop3_DC_OPT_SHADOW_SPECIAL(mat,i,j,LOOP) )    {  
          EstLoop3_DC_OPT_SHADOW_SPECIAL(mat,i,j,LOOP) = temp;   
          /* Have to push only bottem half of system here */ 
          for(k=0;k<3;k++)   
            EstLoop3_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,LOOP,k) = EstLoop3_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,DELETE,k);  
          EstLoop3_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,LOOP,6) = EstLoop3_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,DELETE,6);    
          EstLoop3_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,LOOP,3) = i; 
          EstLoop3_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,LOOP,4) = j; 
          EstLoop3_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,LOOP,5) = DELETE;    
          }  




        /* state DELETE is a source for special END */ 
        temp = score + (mat->query->seg[i]->transition[GW_DELETE2END]) + (0) ;   
        if( temp > EstLoop3_DC_OPT_SHADOW_SPECIAL(mat,i,j,END) )     {  
          EstLoop3_DC_OPT_SHADOW_SPECIAL(mat,i,j,END) = temp;    
          /* Have to push only bottem half of system here */ 
          for(k=0;k<3;k++)   
            EstLoop3_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,k) = EstLoop3_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,DELETE,k);   
          EstLoop3_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,6) = EstLoop3_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,DELETE,6); 
          EstLoop3_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,3) = i;  
          EstLoop3_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,4) = j;  
          EstLoop3_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,5) = DELETE; 
          }  




        /* Finished calculating state DELETE */ 


        } /* end of for each i position in strip */ 


      /* Special state LOOP has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = EstLoop3_DC_OPT_SHADOW_SPECIAL(mat,0,j,LOOP);  


      /* Source START is a special source for LOOP */ 
      temp = EstLoop3_DC_OPT_SHADOW_SPECIAL(mat,0,j - 1,START) + (0) + (0);  
      if( temp > score ) {  
        score = temp;    
        /* Also got to propagate shadows  */ 
        for(k=0;k<7;k++) 
          EstLoop3_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,LOOP,k) = EstLoop3_DC_OPT_SHADOW_SPECIAL_SP(mat,i - 0,j - 1,START,k);    
        }  


      /* Source MATCH for state LOOP is not special... already calculated */ 
      /* Source INSERT for state LOOP is not special... already calculated */ 
      /* Source DELETE for state LOOP is not special... already calculated */ 
      /* Source LOOP is a special source for LOOP */ 
      temp = EstLoop3_DC_OPT_SHADOW_SPECIAL(mat,0,j - 1,LOOP) + (0) + (0);   
      if( temp > score ) {  
        score = temp;    
        /* Also got to propagate shadows  */ 
        for(k=0;k<7;k++) 
          EstLoop3_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,LOOP,k) = EstLoop3_DC_OPT_SHADOW_SPECIAL_SP(mat,i - 0,j - 1,LOOP,k); 
        }  


      /* Put back score... (now updated!) */ 
      EstLoop3_DC_OPT_SHADOW_SPECIAL(mat,0,j,LOOP) = score;  
      /* Finished updating state LOOP */ 




      /* Special state START has no special to special movements */ 


      /* Special state END has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = EstLoop3_DC_OPT_SHADOW_SPECIAL(mat,0,j,END);   


      /* Source LOOP is a special source for END */ 
      temp = EstLoop3_DC_OPT_SHADOW_SPECIAL(mat,0,j - 1,LOOP) + (0) + (0);   
      if( temp > score ) {  
        score = temp;    
        /* Also got to propagate shadows  */ 
        for(k=0;k<7;k++) 
          EstLoop3_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,k) = EstLoop3_DC_OPT_SHADOW_SPECIAL_SP(mat,i - 0,j - 1,LOOP,k);  
        }  


      /* Source MATCH for state END is not special... already calculated */ 
      /* Source INSERT for state END is not special... already calculated */ 
      /* Source DELETE for state END is not special... already calculated */ 
      /* Put back score... (now updated!) */ 
      EstLoop3_DC_OPT_SHADOW_SPECIAL(mat,0,j,END) = score;   
      /* Finished updating state END */ 


      } /* end of for each j strip */ 
    free(score_pointers);    
    free(shadow_pointers);   
    return TRUE;     
}    


/* Function:  init_start_end_linear_EstLoop3(mat)
 *
 * Descrip: No Description
 *
 * Arg:        mat [UNKN ] Undocumented argument [EstLoop3 *]
 *
 */
void init_start_end_linear_EstLoop3(EstLoop3 * mat) 
{
    register int i;  
    register int j;  
    for(j=0;j<7;j++) {  
      for(i=(-1);i<mat->query->len;i++)  {  
        EstLoop3_DC_SHADOW_MATRIX(mat,i,j,MATCH) = NEGI; 
        EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,0) = (-1);    
        EstLoop3_DC_SHADOW_MATRIX(mat,i,j,INSERT) = NEGI;    
        EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,0) = (-1);   
        EstLoop3_DC_SHADOW_MATRIX(mat,i,j,DELETE) = NEGI;    
        EstLoop3_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,0) = (-1);   
        }  
      }  


    for(j=(-5);j<mat->target->seq->len;j++)  {  
      EstLoop3_DC_SHADOW_SPECIAL(mat,0,j,LOOP) = NEGI;   
      EstLoop3_DC_SHADOW_SPECIAL_SP(mat,0,j,LOOP,0) = (-1);  
      EstLoop3_DC_SHADOW_SPECIAL(mat,0,j,START) = 0; 
      EstLoop3_DC_SHADOW_SPECIAL_SP(mat,0,j,START,0) = j;    
      EstLoop3_DC_SHADOW_SPECIAL(mat,0,j,END) = NEGI;    
      EstLoop3_DC_SHADOW_SPECIAL_SP(mat,0,j,END,0) = (-1);   
      }  


    return;  
}    


/* Function:  convert_PackAln_to_AlnBlock_EstLoop3(pal)
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
AlnBlock * convert_PackAln_to_AlnBlock_EstLoop3(PackAln * pal) 
{
    AlnConvertSet * acs; 
    AlnBlock * alb;  


    acs = AlnConvertSet_EstLoop3();  
    alb = AlnBlock_from_PackAln(acs,pal);    
    free_AlnConvertSet(acs); 
    return alb;  
}    


 static char * query_label[] = { "MATCH_STATE","INSERT_STATE","DELETE_STATE","LOOP","END" }; 
/* Function:  AlnConvertSet_EstLoop3(void)
 *
 * Descrip: No Description
 *
 *
 * Return [UNKN ]  Undocumented return value [AlnConvertSet *]
 *
 */
 static char * target_label[] = { "CODON","SEQUENCE_DELETION","SEQUENCE_INSERTION","INSERT","RANDOM_SEQUENCE","END" };   
AlnConvertSet * AlnConvertSet_EstLoop3(void) 
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
    acu->state1 = LOOP + 3;  
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
    acu->state1 = LOOP + 3;  
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
    acu->state1 = LOOP + 3;  
    acu->is_from_special = TRUE; 
    acu->state2 = DELETE;    
    acu->offi = (-1);    
    acu->offj = 0;   
    acu->label1 = query_label[2];    
    acu->label2 = target_label[3];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = START + 3; 
    acu->state2 = LOOP + 3;  
    acu->offi = (-1);    
    acu->offj = 1;   
    acu->label1 = query_label[3];    
    acu->label2 = target_label[4];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = MATCH; 
    acu->state2 = LOOP + 3;  
    acu->offi = (-1);    
    acu->offj = 0;   
    acu->label1 = query_label[3];    
    acu->label2 = target_label[4];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = INSERT;    
    acu->state2 = LOOP + 3;  
    acu->offi = (-1);    
    acu->offj = 0;   
    acu->label1 = query_label[3];    
    acu->label2 = target_label[4];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = DELETE;    
    acu->state2 = LOOP + 3;  
    acu->offi = (-1);    
    acu->offj = 0;   
    acu->label1 = query_label[3];    
    acu->label2 = target_label[4];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = LOOP + 3;  
    acu->state2 = LOOP + 3;  
    acu->offi = (-1);    
    acu->offj = 1;   
    acu->label1 = query_label[3];    
    acu->label2 = target_label[4];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = LOOP + 3;  
    acu->state2 = END + 3;   
    acu->offi = (-1);    
    acu->offj = 1;   
    acu->label1 = query_label[4];    
    acu->label2 = target_label[5];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = MATCH; 
    acu->state2 = END + 3;   
    acu->offi = (-1);    
    acu->offj = 0;   
    acu->label1 = query_label[4];    
    acu->label2 = target_label[5];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = INSERT;    
    acu->state2 = END + 3;   
    acu->offi = (-1);    
    acu->offj = 0;   
    acu->label1 = query_label[4];    
    acu->label2 = target_label[5];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = DELETE;    
    acu->state2 = END + 3;   
    acu->offi = (-1);    
    acu->offj = 0;   
    acu->label1 = query_label[4];    
    acu->label2 = target_label[5];   
    add_collapse_label_AlnConvertSet(out,"LOOP","RANDOM_SEQUENCE");  
    return out;  
}    


/* Function:  PackAln_read_Expl_EstLoop3(mat)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [EstLoop3 *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * PackAln_read_Expl_EstLoop3(EstLoop3 * mat) 
{
    EstLoop3_access_func_holder holder;  


    holder.access_main    = EstLoop3_explicit_access_main;   
    holder.access_special = EstLoop3_explicit_access_special;    
    return PackAln_read_generic_EstLoop3(mat,holder);    
}    


/* Function:  EstLoop3_explicit_access_main(mat,i,j,state)
 *
 * Descrip: No Description
 *
 * Arg:          mat [UNKN ] Undocumented argument [EstLoop3 *]
 * Arg:            i [UNKN ] Undocumented argument [int]
 * Arg:            j [UNKN ] Undocumented argument [int]
 * Arg:        state [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int EstLoop3_explicit_access_main(EstLoop3 * mat,int i,int j,int state) 
{
    return EstLoop3_EXPL_MATRIX(mat,i,j,state);  
}    


/* Function:  EstLoop3_explicit_access_special(mat,i,j,state)
 *
 * Descrip: No Description
 *
 * Arg:          mat [UNKN ] Undocumented argument [EstLoop3 *]
 * Arg:            i [UNKN ] Undocumented argument [int]
 * Arg:            j [UNKN ] Undocumented argument [int]
 * Arg:        state [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int EstLoop3_explicit_access_special(EstLoop3 * mat,int i,int j,int state) 
{
    return EstLoop3_EXPL_SPECIAL(mat,i,j,state); 
}    


/* Function:  PackAln_read_generic_EstLoop3(mat,h)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [EstLoop3 *]
 * Arg:          h [UNKN ] Undocumented argument [EstLoop3_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * PackAln_read_generic_EstLoop3(EstLoop3 * mat,EstLoop3_access_func_holder h) 
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


    out->score =  find_end_EstLoop3(mat,&i,&j,&state,&isspecial,h);  


    /* Add final end transition (at the moment we have not got the score! */ 
    if( (pau= PackAlnUnit_alloc()) == NULL  || add_PackAln(out,pau) == FALSE )   {  
      warn("Failed the first PackAlnUnit alloc, %d length of Alignment in EstLoop3_basic_read, returning a mess.(Sorry!)",out->len); 
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
        max_calc_special_EstLoop3(mat,i,j,state,isspecial,&i,&j,&state,&isspecial,&cellscore,h);     
      else   
        max_calc_EstLoop3(mat,i,j,state,isspecial,&i,&j,&state,&isspecial,&cellscore,h);     
      if(i == EstLoop3_READ_OFF_ERROR || j == EstLoop3_READ_OFF_ERROR || state == EstLoop3_READ_OFF_ERROR )  {  
        warn("Problem - hit bad read off system, exiting now");  
        break;   
        }  
      if( (pau= PackAlnUnit_alloc()) == NULL  || add_PackAln(out,pau) == FALSE ) {  
        warn("Failed a PackAlnUnit alloc, %d length of Alignment in EstLoop3_basic_read, returning partial alignment",out->len); 
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


/* Function:  find_end_EstLoop3(mat,ri,rj,state,isspecial,h)
 *
 * Descrip: No Description
 *
 * Arg:              mat [UNKN ] Undocumented argument [EstLoop3 *]
 * Arg:               ri [UNKN ] Undocumented argument [int *]
 * Arg:               rj [UNKN ] Undocumented argument [int *]
 * Arg:            state [UNKN ] Undocumented argument [int *]
 * Arg:        isspecial [UNKN ] Undocumented argument [boolean *]
 * Arg:                h [UNKN ] Undocumented argument [EstLoop3_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int find_end_EstLoop3(EstLoop3 * mat,int * ri,int * rj,int * state,boolean * isspecial,EstLoop3_access_func_holder h) 
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


/* Function:  EstLoop3_debug_show_matrix(mat,starti,stopi,startj,stopj,ofp)
 *
 * Descrip: No Description
 *
 * Arg:           mat [UNKN ] Undocumented argument [EstLoop3 *]
 * Arg:        starti [UNKN ] Undocumented argument [int]
 * Arg:         stopi [UNKN ] Undocumented argument [int]
 * Arg:        startj [UNKN ] Undocumented argument [int]
 * Arg:         stopj [UNKN ] Undocumented argument [int]
 * Arg:           ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
void EstLoop3_debug_show_matrix(EstLoop3 * mat,int starti,int stopi,int startj,int stopj,FILE * ofp) 
{
    register int i;  
    register int j;  


    for(i=starti;i<stopi && i < mat->query->len;i++) {  
      for(j=startj;j<stopj && j < mat->target->seq->len;j++) {  
        fprintf(ofp,"Cell [%d - %d]\n",i,j);     
        fprintf(ofp,"State MATCH %d\n",EstLoop3_EXPL_MATRIX(mat,i,j,MATCH)); 
        fprintf(ofp,"State INSERT %d\n",EstLoop3_EXPL_MATRIX(mat,i,j,INSERT));   
        fprintf(ofp,"State DELETE %d\n",EstLoop3_EXPL_MATRIX(mat,i,j,DELETE));   
        fprintf(ofp,"\n\n"); 
        }  
      }  


}    


/* Function:  max_calc_EstLoop3(mat,i,j,state,isspecial,reti,retj,retstate,retspecial,cellscore,h)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [EstLoop3 *]
 * Arg:                 i [UNKN ] Undocumented argument [int]
 * Arg:                 j [UNKN ] Undocumented argument [int]
 * Arg:             state [UNKN ] Undocumented argument [int]
 * Arg:         isspecial [UNKN ] Undocumented argument [boolean]
 * Arg:              reti [UNKN ] Undocumented argument [int *]
 * Arg:              retj [UNKN ] Undocumented argument [int *]
 * Arg:          retstate [UNKN ] Undocumented argument [int *]
 * Arg:        retspecial [UNKN ] Undocumented argument [boolean *]
 * Arg:         cellscore [UNKN ] Undocumented argument [int *]
 * Arg:                 h [UNKN ] Undocumented argument [EstLoop3_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int max_calc_EstLoop3(EstLoop3 * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,EstLoop3_access_func_holder h) 
{
    register int temp;   
    register int cscore; 


    *reti = (*retj) = (*retstate) = EstLoop3_READ_OFF_ERROR; 


    if( i < 0 || j < 0 || i > mat->query->len || j > mat->target->seq->len)  {  
      warn("In EstLoop3 matrix special read off - out of bounds on matrix [i,j is %d,%d state %d in standard matrix]",i,j,state);    
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
        temp = cscore - ((mat->query->seg[i]->transition[GW_START2MATCH]+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)])) -  (0); 
        if( temp == (*h.access_special)(mat,i - 1,j - 3,LOOP) )  {  
          *reti = i - 1; 
          *retj = j - 3; 
          *retstate = LOOP;  
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_special)(mat,i-1,j-3,LOOP); 
            }  
          return (*h.access_main)(mat,i - 1,j - 3,LOOP);     
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
        warn("Major problem (!) - in EstLoop3 read off, position %d,%d state %d no source found!",i,j,state);    
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
        if( temp == (*h.access_special)(mat,i - 0,j - 3,LOOP) )  {  
          *reti = i - 0; 
          *retj = j - 3; 
          *retstate = LOOP;  
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_special)(mat,i-0,j-3,LOOP); 
            }  
          return (*h.access_main)(mat,i - 0,j - 3,LOOP);     
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
        warn("Major problem (!) - in EstLoop3 read off, position %d,%d state %d no source found!",i,j,state);    
        return (-1); 
      case DELETE :  
        temp = cscore - (mat->query->seg[i]->transition[GW_START2DELETE]) -  (0);    
        if( temp == (*h.access_special)(mat,i - 1,j - 0,LOOP) )  {  
          *reti = i - 1; 
          *retj = j - 0; 
          *retstate = LOOP;  
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_special)(mat,i-1,j-0,LOOP); 
            }  
          return (*h.access_main)(mat,i - 1,j - 0,LOOP);     
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
        warn("Major problem (!) - in EstLoop3 read off, position %d,%d state %d no source found!",i,j,state);    
        return (-1); 
      default:   
        warn("Major problem (!) - in EstLoop3 read off, position %d,%d state %d no source found!",i,j,state);    
        return (-1); 
      } /* end of Switch state  */ 
}    


/* Function:  max_calc_special_EstLoop3(mat,i,j,state,isspecial,reti,retj,retstate,retspecial,cellscore,h)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [EstLoop3 *]
 * Arg:                 i [UNKN ] Undocumented argument [int]
 * Arg:                 j [UNKN ] Undocumented argument [int]
 * Arg:             state [UNKN ] Undocumented argument [int]
 * Arg:         isspecial [UNKN ] Undocumented argument [boolean]
 * Arg:              reti [UNKN ] Undocumented argument [int *]
 * Arg:              retj [UNKN ] Undocumented argument [int *]
 * Arg:          retstate [UNKN ] Undocumented argument [int *]
 * Arg:        retspecial [UNKN ] Undocumented argument [boolean *]
 * Arg:         cellscore [UNKN ] Undocumented argument [int *]
 * Arg:                 h [UNKN ] Undocumented argument [EstLoop3_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int max_calc_special_EstLoop3(EstLoop3 * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,EstLoop3_access_func_holder h) 
{
    register int temp;   
    register int cscore; 


    *reti = (*retj) = (*retstate) = EstLoop3_READ_OFF_ERROR; 


    if( j < 0 || j > mat->target->seq->len)  {  
      warn("In EstLoop3 matrix special read off - out of bounds on matrix [j is %d in special]",j);  
      return -1;     
      }  


    cscore = (*h.access_special)(mat,i,j,state); 
    switch(state)    { /*switch on special states*/ 
      case LOOP :    
        /* source LOOP is a special */ 
        temp = cscore - (0) - (0);   
        if( temp == (*h.access_special)(mat,i - 0,j - 1,LOOP) )  {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = LOOP;  
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_special)(mat,i-0,j-1,LOOP);     
            }  
          return (*h.access_special)(mat,i - 0,j - 1,LOOP) ;     
          }  
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
        /* source START is a special */ 
        temp = cscore - (0) - (0);   
        if( temp == (*h.access_special)(mat,i - 0,j - 1,START) ) {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = START; 
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_special)(mat,i-0,j-1,START);    
            }  
          return (*h.access_special)(mat,i - 0,j - 1,START) ;    
          }  
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
        /* source LOOP is a special */ 
        temp = cscore - (0) - (0);   
        if( temp == (*h.access_special)(mat,i - 0,j - 1,LOOP) )  {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = LOOP;  
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_special)(mat,i-0,j-1,LOOP);     
            }  
          return (*h.access_special)(mat,i - 0,j - 1,LOOP) ;     
          }  
      default:   
        warn("Major problem (!) - in EstLoop3 read off, position %d,%d state %d no source found  dropped into default on source switch!",i,j,state); 
        return (-1); 
      } /* end of switch on special states */ 
}    


/* Function:  calculate_EstLoop3(mat)
 *
 * Descrip:    This function calculates the EstLoop3 matrix when in explicit mode
 *             To allocate the matrix use /allocate_Expl_EstLoop3
 *
 *
 * Arg:        mat [UNKN ] EstLoop3 which contains explicit basematrix memory [EstLoop3 *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean calculate_EstLoop3(EstLoop3 * mat) 
{
    int i;   
    int j;   
    int leni;    
    int lenj;    
    int tot; 
    int num; 


    if( mat->basematrix->type != BASEMATRIX_TYPE_EXPLICIT )  {  
      warn("in calculate_EstLoop3, passed a non Explicit matrix type, cannot calculate!");   
      return FALSE;  
      }  


    leni = mat->leni;    
    lenj = mat->lenj;    
    tot = leni * lenj;   
    num = 0; 


    start_reporting("EstLoop3 Matrix calculation: ");    
    for(j=0;j<lenj;j++)  {  
      auto int score;    
      auto int temp;     
      for(i=0;i<leni;i++)    {  
        if( num%1000 == 0 )  
          log_full_error(REPORT,0,"[%7d] Cells %2d%%%%",num,num*100/tot);    
        num++;   


        /* For state MATCH */ 
        /* setting first movement to score */ 
        score = EstLoop3_EXPL_MATRIX(mat,i-1,j-3,MATCH) + (mat->query->seg[i]->transition[GW_MATCH2MATCH]+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)]);    
        /* From state INSERT to state MATCH */ 
        temp = EstLoop3_EXPL_MATRIX(mat,i-1,j-3,INSERT) + (mat->query->seg[i]->transition[GW_INSERT2MATCH]+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)]);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state MATCH */ 
        temp = EstLoop3_EXPL_MATRIX(mat,i-1,j-3,DELETE) + (mat->query->seg[i]->transition[GW_DELETE2MATCH]+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)]);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state LOOP to state MATCH */ 
        temp = EstLoop3_EXPL_SPECIAL(mat,i-1,j-3,LOOP) + (mat->query->seg[i]->transition[GW_START2MATCH]+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)]);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state START to state MATCH */ 
        temp = EstLoop3_EXPL_SPECIAL(mat,i-1,j-3,START) + (mat->query->seg[i]->transition[GW_START2MATCH]+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)]);    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = EstLoop3_EXPL_MATRIX(mat,i-1,j-2,MATCH) + mat->cp->trans[PCD_DELETE_1_BASE];  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = EstLoop3_EXPL_MATRIX(mat,i-1,j-1,MATCH) + mat->cp->trans[PCD_DELETE_2_BASE];  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = EstLoop3_EXPL_MATRIX(mat,i-1,j-4,MATCH) + mat->cp->trans[PCD_INSERT_1_BASE];  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = EstLoop3_EXPL_MATRIX(mat,i-1,j-5,MATCH) + mat->cp->trans[PCD_INSERT_2_BASE];  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state MATCH */ 
        temp = EstLoop3_EXPL_MATRIX(mat,i-1,j-2,INSERT) + mat->cp->trans[PCD_DELETE_1_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state MATCH */ 
        temp = EstLoop3_EXPL_MATRIX(mat,i-1,j-1,INSERT) + mat->cp->trans[PCD_DELETE_2_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state MATCH */ 
        temp = EstLoop3_EXPL_MATRIX(mat,i-1,j-4,INSERT) + mat->cp->trans[PCD_INSERT_1_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state MATCH */ 
        temp = EstLoop3_EXPL_MATRIX(mat,i-1,j-5,INSERT) + mat->cp->trans[PCD_INSERT_2_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state MATCH */ 
        temp = EstLoop3_EXPL_MATRIX(mat,i-1,j-2,DELETE) + mat->cp->trans[PCD_DELETE_1_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state MATCH */ 
        temp = EstLoop3_EXPL_MATRIX(mat,i-1,j-1,DELETE) + mat->cp->trans[PCD_DELETE_2_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state MATCH */ 
        temp = EstLoop3_EXPL_MATRIX(mat,i-1,j-4,DELETE) + mat->cp->trans[PCD_INSERT_1_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state MATCH */ 
        temp = EstLoop3_EXPL_MATRIX(mat,i-1,j-5,DELETE) + mat->cp->trans[PCD_INSERT_2_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for MATCH */ 
        /* Add any movement independant score and put away */ 
         EstLoop3_EXPL_MATRIX(mat,i,j,MATCH) = score;    


        /* state MATCH is a source for special LOOP */ 
        temp = score + (mat->query->seg[i]->transition[GW_MATCH2END]) + (0) ;    
        if( temp > EstLoop3_EXPL_SPECIAL(mat,i,j,LOOP) )     {  
          EstLoop3_EXPL_SPECIAL(mat,i,j,LOOP) = temp;    
          }  




        /* state MATCH is a source for special END */ 
        temp = score + (mat->query->seg[i]->transition[GW_MATCH2END]) + (0) ;    
        if( temp > EstLoop3_EXPL_SPECIAL(mat,i,j,END) )  {  
          EstLoop3_EXPL_SPECIAL(mat,i,j,END) = temp;     
          }  




        /* Finished calculating state MATCH */ 


        /* For state INSERT */ 
        /* setting first movement to score */ 
        score = EstLoop3_EXPL_MATRIX(mat,i-0,j-3,MATCH) + (mat->query->seg[i]->transition[GW_MATCH2INSERT]+mat->query->seg[i]->insert[CSEQ_CDNA_CODON(mat->target,j)]);  
        /* From state INSERT to state INSERT */ 
        temp = EstLoop3_EXPL_MATRIX(mat,i-0,j-3,INSERT) + (mat->query->seg[i]->transition[GW_INSERT2INSERT]+mat->query->seg[i]->insert[CSEQ_CDNA_CODON(mat->target,j)]);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state INSERT */ 
        temp = EstLoop3_EXPL_MATRIX(mat,i-0,j-3,DELETE) + (mat->query->seg[i]->transition[GW_DELETE2INSERT]+mat->query->seg[i]->insert[CSEQ_CDNA_CODON(mat->target,j)]);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state LOOP to state INSERT */ 
        temp = EstLoop3_EXPL_SPECIAL(mat,i-0,j-3,LOOP) + (mat->query->seg[i]->transition[GW_START2INSERT]+mat->query->seg[i]->insert[CSEQ_CDNA_CODON(mat->target,j)]);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state INSERT */ 
        temp = EstLoop3_EXPL_MATRIX(mat,i-0,j-2,INSERT) + mat->cp->trans[PCD_DELETE_1_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state INSERT */ 
        temp = EstLoop3_EXPL_MATRIX(mat,i-0,j-1,INSERT) + mat->cp->trans[PCD_DELETE_2_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state INSERT */ 
        temp = EstLoop3_EXPL_MATRIX(mat,i-0,j-4,INSERT) + mat->cp->trans[PCD_INSERT_1_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state INSERT */ 
        temp = EstLoop3_EXPL_MATRIX(mat,i-0,j-5,INSERT) + mat->cp->trans[PCD_INSERT_2_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state INSERT */ 
        temp = EstLoop3_EXPL_MATRIX(mat,i-0,j-2,MATCH) + mat->cp->trans[PCD_DELETE_1_BASE];  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state INSERT */ 
        temp = EstLoop3_EXPL_MATRIX(mat,i-0,j-1,MATCH) + mat->cp->trans[PCD_DELETE_2_BASE];  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state INSERT */ 
        temp = EstLoop3_EXPL_MATRIX(mat,i-0,j-4,MATCH) + mat->cp->trans[PCD_INSERT_1_BASE];  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state INSERT */ 
        temp = EstLoop3_EXPL_MATRIX(mat,i-0,j-5,MATCH) + mat->cp->trans[PCD_INSERT_2_BASE];  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state INSERT */ 
        temp = EstLoop3_EXPL_MATRIX(mat,i-0,j-2,DELETE) + mat->cp->trans[PCD_DELETE_1_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state INSERT */ 
        temp = EstLoop3_EXPL_MATRIX(mat,i-0,j-1,DELETE) + mat->cp->trans[PCD_DELETE_2_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state INSERT */ 
        temp = EstLoop3_EXPL_MATRIX(mat,i-0,j-4,DELETE) + mat->cp->trans[PCD_INSERT_1_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state INSERT */ 
        temp = EstLoop3_EXPL_MATRIX(mat,i-0,j-5,DELETE) + mat->cp->trans[PCD_INSERT_2_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for INSERT */ 
        /* Add any movement independant score and put away */ 
         EstLoop3_EXPL_MATRIX(mat,i,j,INSERT) = score;   


        /* state INSERT is a source for special LOOP */ 
        temp = score + (mat->query->seg[i]->transition[GW_INSERT2END]) + (0) ;   
        if( temp > EstLoop3_EXPL_SPECIAL(mat,i,j,LOOP) )     {  
          EstLoop3_EXPL_SPECIAL(mat,i,j,LOOP) = temp;    
          }  




        /* state INSERT is a source for special END */ 
        temp = score + (mat->query->seg[i]->transition[GW_INSERT2END]) + (0) ;   
        if( temp > EstLoop3_EXPL_SPECIAL(mat,i,j,END) )  {  
          EstLoop3_EXPL_SPECIAL(mat,i,j,END) = temp;     
          }  




        /* Finished calculating state INSERT */ 


        /* For state DELETE */ 
        /* setting first movement to score */ 
        score = EstLoop3_EXPL_MATRIX(mat,i-1,j-0,MATCH) + mat->query->seg[i]->transition[GW_MATCH2DELETE];   
        /* From state INSERT to state DELETE */ 
        temp = EstLoop3_EXPL_MATRIX(mat,i-1,j-0,INSERT) + mat->query->seg[i]->transition[GW_INSERT2DELETE];  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state DELETE */ 
        temp = EstLoop3_EXPL_MATRIX(mat,i-1,j-0,DELETE) + mat->query->seg[i]->transition[GW_DELETE2DELETE];  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state LOOP to state DELETE */ 
        temp = EstLoop3_EXPL_SPECIAL(mat,i-1,j-0,LOOP) + mat->query->seg[i]->transition[GW_START2DELETE];    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for DELETE */ 
        /* Add any movement independant score and put away */ 
         EstLoop3_EXPL_MATRIX(mat,i,j,DELETE) = score;   


        /* state DELETE is a source for special LOOP */ 
        temp = score + (mat->query->seg[i]->transition[GW_DELETE2END]) + (0) ;   
        if( temp > EstLoop3_EXPL_SPECIAL(mat,i,j,LOOP) )     {  
          EstLoop3_EXPL_SPECIAL(mat,i,j,LOOP) = temp;    
          }  




        /* state DELETE is a source for special END */ 
        temp = score + (mat->query->seg[i]->transition[GW_DELETE2END]) + (0) ;   
        if( temp > EstLoop3_EXPL_SPECIAL(mat,i,j,END) )  {  
          EstLoop3_EXPL_SPECIAL(mat,i,j,END) = temp;     
          }  




        /* Finished calculating state DELETE */ 
        }  


      /* Special state LOOP has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = EstLoop3_EXPL_SPECIAL(mat,0,j,LOOP);   


      /* Source START is a special source for LOOP */ 
      /* Has restricted position */ 
      if( (j-1) == 0  )  {  
        temp = EstLoop3_EXPL_SPECIAL(mat,0,j - 1,START) + (0) + (0);     
        if( temp > score )   
          score = temp;  
        }  


      /* Source MATCH for state LOOP is not special... already calculated */ 
      /* Source INSERT for state LOOP is not special... already calculated */ 
      /* Source DELETE for state LOOP is not special... already calculated */ 
      /* Source LOOP is a special source for LOOP */ 
      temp = EstLoop3_EXPL_SPECIAL(mat,0,j - 1,LOOP) + (0) + (0);    
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      EstLoop3_EXPL_SPECIAL(mat,0,j,LOOP) = score;   
      /* Finished updating state LOOP */ 




      /* Special state START has no special to special movements */ 


      /* Special state END has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = EstLoop3_EXPL_SPECIAL(mat,0,j,END);    


      /* Source LOOP is a special source for END */ 
      /* Has restricted position */ 
      if( j == mat->lenj-1 ) {  
        temp = EstLoop3_EXPL_SPECIAL(mat,0,j - 1,LOOP) + (0) + (0);  
        if( temp > score )   
          score = temp;  
        }  


      /* Source MATCH for state END is not special... already calculated */ 
      /* Source INSERT for state END is not special... already calculated */ 
      /* Source DELETE for state END is not special... already calculated */ 
      /* Put back score... (now updated!) */ 
      EstLoop3_EXPL_SPECIAL(mat,0,j,END) = score;    
      /* Finished updating state END */ 


      }  
    stop_reporting();    
    return TRUE;     
}    


/* Function:  calculate_dpenv_EstLoop3(mat,dpenv)
 *
 * Descrip:    This function calculates the EstLoop3 matrix when in explicit mode, subject to the envelope
 *
 *
 * Arg:          mat [UNKN ] EstLoop3 which contains explicit basematrix memory [EstLoop3 *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean calculate_dpenv_EstLoop3(EstLoop3 * mat,DPEnvelope * dpenv) 
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
      warn("in calculate_EstLoop3, passed a non Explicit matrix type, cannot calculate!");   
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
        EstLoop3_EXPL_MATRIX(mat,i,j,MATCH) = NEGI;  
        EstLoop3_EXPL_MATRIX(mat,i,j,INSERT) = NEGI; 
        EstLoop3_EXPL_MATRIX(mat,i,j,DELETE) = NEGI; 
        }  
      }  
    for(j=-5;j<mat->lenj;j++)    {  
      EstLoop3_EXPL_SPECIAL(mat,i,j,LOOP) = NEGI;    
      EstLoop3_EXPL_SPECIAL(mat,i,j,START) = 0;  
      EstLoop3_EXPL_SPECIAL(mat,i,j,END) = NEGI; 
      }  


    start_reporting("EstLoop3 Matrix calculation: ");    
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
          EstLoop3_EXPL_MATRIX(mat,i,j,MATCH) = NEGI;    
          EstLoop3_EXPL_MATRIX(mat,i,j,INSERT) = NEGI;   
          EstLoop3_EXPL_MATRIX(mat,i,j,DELETE) = NEGI;   
          continue;  
          }  


        if( num%1000 == 0 )  
          log_full_error(REPORT,0,"[%7d] Cells %2d%%%%",num,num*100/tot);    
        num++;   


        /* For state MATCH */ 
        /* setting first movement to score */ 
        score = EstLoop3_EXPL_MATRIX(mat,i-1,j-3,MATCH) + (mat->query->seg[i]->transition[GW_MATCH2MATCH]+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)]);    
        /* From state INSERT to state MATCH */ 
        temp = EstLoop3_EXPL_MATRIX(mat,i-1,j-3,INSERT) + (mat->query->seg[i]->transition[GW_INSERT2MATCH]+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)]);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state MATCH */ 
        temp = EstLoop3_EXPL_MATRIX(mat,i-1,j-3,DELETE) + (mat->query->seg[i]->transition[GW_DELETE2MATCH]+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)]);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state LOOP to state MATCH */ 
        temp = EstLoop3_EXPL_SPECIAL(mat,i-1,j-3,LOOP) + (mat->query->seg[i]->transition[GW_START2MATCH]+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)]);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state START to state MATCH */ 
        temp = EstLoop3_EXPL_SPECIAL(mat,i-1,j-3,START) + (mat->query->seg[i]->transition[GW_START2MATCH]+mat->query->seg[i]->match[CSEQ_CDNA_CODON(mat->target,j)]);    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = EstLoop3_EXPL_MATRIX(mat,i-1,j-2,MATCH) + mat->cp->trans[PCD_DELETE_1_BASE];  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = EstLoop3_EXPL_MATRIX(mat,i-1,j-1,MATCH) + mat->cp->trans[PCD_DELETE_2_BASE];  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = EstLoop3_EXPL_MATRIX(mat,i-1,j-4,MATCH) + mat->cp->trans[PCD_INSERT_1_BASE];  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = EstLoop3_EXPL_MATRIX(mat,i-1,j-5,MATCH) + mat->cp->trans[PCD_INSERT_2_BASE];  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state MATCH */ 
        temp = EstLoop3_EXPL_MATRIX(mat,i-1,j-2,INSERT) + mat->cp->trans[PCD_DELETE_1_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state MATCH */ 
        temp = EstLoop3_EXPL_MATRIX(mat,i-1,j-1,INSERT) + mat->cp->trans[PCD_DELETE_2_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state MATCH */ 
        temp = EstLoop3_EXPL_MATRIX(mat,i-1,j-4,INSERT) + mat->cp->trans[PCD_INSERT_1_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state MATCH */ 
        temp = EstLoop3_EXPL_MATRIX(mat,i-1,j-5,INSERT) + mat->cp->trans[PCD_INSERT_2_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state MATCH */ 
        temp = EstLoop3_EXPL_MATRIX(mat,i-1,j-2,DELETE) + mat->cp->trans[PCD_DELETE_1_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state MATCH */ 
        temp = EstLoop3_EXPL_MATRIX(mat,i-1,j-1,DELETE) + mat->cp->trans[PCD_DELETE_2_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state MATCH */ 
        temp = EstLoop3_EXPL_MATRIX(mat,i-1,j-4,DELETE) + mat->cp->trans[PCD_INSERT_1_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state MATCH */ 
        temp = EstLoop3_EXPL_MATRIX(mat,i-1,j-5,DELETE) + mat->cp->trans[PCD_INSERT_2_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for MATCH */ 
        /* Add any movement independant score and put away */ 
         EstLoop3_EXPL_MATRIX(mat,i,j,MATCH) = score;    


        /* state MATCH is a source for special LOOP */ 
        temp = score + (mat->query->seg[i]->transition[GW_MATCH2END]) + (0) ;    
        if( temp > EstLoop3_EXPL_SPECIAL(mat,i,j,LOOP) )     {  
          EstLoop3_EXPL_SPECIAL(mat,i,j,LOOP) = temp;    
          }  




        /* state MATCH is a source for special END */ 
        temp = score + (mat->query->seg[i]->transition[GW_MATCH2END]) + (0) ;    
        if( temp > EstLoop3_EXPL_SPECIAL(mat,i,j,END) )  {  
          EstLoop3_EXPL_SPECIAL(mat,i,j,END) = temp;     
          }  




        /* Finished calculating state MATCH */ 


        /* For state INSERT */ 
        /* setting first movement to score */ 
        score = EstLoop3_EXPL_MATRIX(mat,i-0,j-3,MATCH) + (mat->query->seg[i]->transition[GW_MATCH2INSERT]+mat->query->seg[i]->insert[CSEQ_CDNA_CODON(mat->target,j)]);  
        /* From state INSERT to state INSERT */ 
        temp = EstLoop3_EXPL_MATRIX(mat,i-0,j-3,INSERT) + (mat->query->seg[i]->transition[GW_INSERT2INSERT]+mat->query->seg[i]->insert[CSEQ_CDNA_CODON(mat->target,j)]);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state INSERT */ 
        temp = EstLoop3_EXPL_MATRIX(mat,i-0,j-3,DELETE) + (mat->query->seg[i]->transition[GW_DELETE2INSERT]+mat->query->seg[i]->insert[CSEQ_CDNA_CODON(mat->target,j)]);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state LOOP to state INSERT */ 
        temp = EstLoop3_EXPL_SPECIAL(mat,i-0,j-3,LOOP) + (mat->query->seg[i]->transition[GW_START2INSERT]+mat->query->seg[i]->insert[CSEQ_CDNA_CODON(mat->target,j)]);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state INSERT */ 
        temp = EstLoop3_EXPL_MATRIX(mat,i-0,j-2,INSERT) + mat->cp->trans[PCD_DELETE_1_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state INSERT */ 
        temp = EstLoop3_EXPL_MATRIX(mat,i-0,j-1,INSERT) + mat->cp->trans[PCD_DELETE_2_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state INSERT */ 
        temp = EstLoop3_EXPL_MATRIX(mat,i-0,j-4,INSERT) + mat->cp->trans[PCD_INSERT_1_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state INSERT */ 
        temp = EstLoop3_EXPL_MATRIX(mat,i-0,j-5,INSERT) + mat->cp->trans[PCD_INSERT_2_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state INSERT */ 
        temp = EstLoop3_EXPL_MATRIX(mat,i-0,j-2,MATCH) + mat->cp->trans[PCD_DELETE_1_BASE];  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state INSERT */ 
        temp = EstLoop3_EXPL_MATRIX(mat,i-0,j-1,MATCH) + mat->cp->trans[PCD_DELETE_2_BASE];  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state INSERT */ 
        temp = EstLoop3_EXPL_MATRIX(mat,i-0,j-4,MATCH) + mat->cp->trans[PCD_INSERT_1_BASE];  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state INSERT */ 
        temp = EstLoop3_EXPL_MATRIX(mat,i-0,j-5,MATCH) + mat->cp->trans[PCD_INSERT_2_BASE];  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state INSERT */ 
        temp = EstLoop3_EXPL_MATRIX(mat,i-0,j-2,DELETE) + mat->cp->trans[PCD_DELETE_1_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state INSERT */ 
        temp = EstLoop3_EXPL_MATRIX(mat,i-0,j-1,DELETE) + mat->cp->trans[PCD_DELETE_2_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state INSERT */ 
        temp = EstLoop3_EXPL_MATRIX(mat,i-0,j-4,DELETE) + mat->cp->trans[PCD_INSERT_1_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state INSERT */ 
        temp = EstLoop3_EXPL_MATRIX(mat,i-0,j-5,DELETE) + mat->cp->trans[PCD_INSERT_2_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for INSERT */ 
        /* Add any movement independant score and put away */ 
         EstLoop3_EXPL_MATRIX(mat,i,j,INSERT) = score;   


        /* state INSERT is a source for special LOOP */ 
        temp = score + (mat->query->seg[i]->transition[GW_INSERT2END]) + (0) ;   
        if( temp > EstLoop3_EXPL_SPECIAL(mat,i,j,LOOP) )     {  
          EstLoop3_EXPL_SPECIAL(mat,i,j,LOOP) = temp;    
          }  




        /* state INSERT is a source for special END */ 
        temp = score + (mat->query->seg[i]->transition[GW_INSERT2END]) + (0) ;   
        if( temp > EstLoop3_EXPL_SPECIAL(mat,i,j,END) )  {  
          EstLoop3_EXPL_SPECIAL(mat,i,j,END) = temp;     
          }  




        /* Finished calculating state INSERT */ 


        /* For state DELETE */ 
        /* setting first movement to score */ 
        score = EstLoop3_EXPL_MATRIX(mat,i-1,j-0,MATCH) + mat->query->seg[i]->transition[GW_MATCH2DELETE];   
        /* From state INSERT to state DELETE */ 
        temp = EstLoop3_EXPL_MATRIX(mat,i-1,j-0,INSERT) + mat->query->seg[i]->transition[GW_INSERT2DELETE];  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state DELETE */ 
        temp = EstLoop3_EXPL_MATRIX(mat,i-1,j-0,DELETE) + mat->query->seg[i]->transition[GW_DELETE2DELETE];  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state LOOP to state DELETE */ 
        temp = EstLoop3_EXPL_SPECIAL(mat,i-1,j-0,LOOP) + mat->query->seg[i]->transition[GW_START2DELETE];    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for DELETE */ 
        /* Add any movement independant score and put away */ 
         EstLoop3_EXPL_MATRIX(mat,i,j,DELETE) = score;   


        /* state DELETE is a source for special LOOP */ 
        temp = score + (mat->query->seg[i]->transition[GW_DELETE2END]) + (0) ;   
        if( temp > EstLoop3_EXPL_SPECIAL(mat,i,j,LOOP) )     {  
          EstLoop3_EXPL_SPECIAL(mat,i,j,LOOP) = temp;    
          }  




        /* state DELETE is a source for special END */ 
        temp = score + (mat->query->seg[i]->transition[GW_DELETE2END]) + (0) ;   
        if( temp > EstLoop3_EXPL_SPECIAL(mat,i,j,END) )  {  
          EstLoop3_EXPL_SPECIAL(mat,i,j,END) = temp;     
          }  




        /* Finished calculating state DELETE */ 
        }  


      /* Special state LOOP has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = EstLoop3_EXPL_SPECIAL(mat,0,j,LOOP);   


      /* Source START is a special source for LOOP */ 
      /* Has restricted position */ 
      if( (j-1) == 0  )  {  
        temp = EstLoop3_EXPL_SPECIAL(mat,0,j - 1,START) + (0) + (0);     
        if( temp > score )   
          score = temp;  
        }  


      /* Source MATCH for state LOOP is not special... already calculated */ 
      /* Source INSERT for state LOOP is not special... already calculated */ 
      /* Source DELETE for state LOOP is not special... already calculated */ 
      /* Source LOOP is a special source for LOOP */ 
      temp = EstLoop3_EXPL_SPECIAL(mat,0,j - 1,LOOP) + (0) + (0);    
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      EstLoop3_EXPL_SPECIAL(mat,0,j,LOOP) = score;   
      /* Finished updating state LOOP */ 




      /* Special state START has no special to special movements */ 


      /* Special state END has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = EstLoop3_EXPL_SPECIAL(mat,0,j,END);    


      /* Source LOOP is a special source for END */ 
      /* Has restricted position */ 
      if( j == mat->lenj-1 ) {  
        temp = EstLoop3_EXPL_SPECIAL(mat,0,j - 1,LOOP) + (0) + (0);  
        if( temp > score )   
          score = temp;  
        }  


      /* Source MATCH for state END is not special... already calculated */ 
      /* Source INSERT for state END is not special... already calculated */ 
      /* Source DELETE for state END is not special... already calculated */ 
      /* Put back score... (now updated!) */ 
      EstLoop3_EXPL_SPECIAL(mat,0,j,END) = score;    
      /* Finished updating state END */ 


      }  
    stop_reporting();    
    return TRUE;     
}    


/* Function:  EstLoop3_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [EstLoop3 *]
 *
 */
EstLoop3 * EstLoop3_alloc(void) 
{
    EstLoop3 * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(EstLoop3 *) ckalloc (sizeof(EstLoop3))) == NULL)    {  
      warn("EstLoop3_alloc failed ");    
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


/* Function:  free_EstLoop3(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [EstLoop3 *]
 *
 * Return [UNKN ]  Undocumented return value [EstLoop3 *]
 *
 */
EstLoop3 * free_EstLoop3(EstLoop3 * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a EstLoop3 obj. Should be trappable");  
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


    ckfree(obj); 
    return NULL; 
}    





#ifdef _cplusplus
}
#endif
