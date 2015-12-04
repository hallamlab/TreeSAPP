#ifdef _cplusplus
extern "C" {
#endif
#include "pba.h"


# line 6 "pba.c"


  /*****************   C functions  ****************/
  /*             Written using dynamite            */
  /*            Sat Sep  8 09:05:32 2007           */
  /*            email birney@sanger.ac.uk          */
  /* http://www.sanger.ac.uk/Users/birney/dynamite */
  /*************************************************/


  /* Please report any problems or bugs to         */
  /* Ewan Birney, birney@sanger.ac.uk              */


/* basic set of macros to map states to numbers */ 
#define BLOCK_1 0    
#define BLOCK_2 1    
#define BLOCK_3 2    
#define UNALIGNED 3  


#define START 0  
#define END 1    


#define ProteinBlockAligner_EXPL_MATRIX(this_matrix,i,j,STATE) this_matrix->basematrix->matrix[((j+1)*4)+STATE][i+1] 
#define ProteinBlockAligner_EXPL_SPECIAL(matrix,i,j,STATE) matrix->basematrix->specmatrix[STATE][j+1]    
#define ProteinBlockAligner_READ_OFF_ERROR -3
   


#define ProteinBlockAligner_VSMALL_MATRIX(mat,i,j,STATE) mat->basematrix->matrix[(j+2)%2][((i+1)*4)+STATE]   
#define ProteinBlockAligner_VSMALL_SPECIAL(mat,i,j,STATE) mat->basematrix->specmatrix[(j+2)%2][STATE]    




#define ProteinBlockAligner_SHATTER_SPECIAL(matrix,i,j,STATE) matrix->shatter->special[STATE][j] 
#define ProteinBlockAligner_SHATTER_MATRIX(matrix,i,j,STATE)  fetch_cell_value_ShatterMatrix(mat->shatter,i,j,STATE) 


/* Function:  PackAln_read_Shatter_ProteinBlockAligner(mat)
 *
 * Descrip:    Reads off PackAln from shatter matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [ProteinBlockAligner *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * PackAln_read_Shatter_ProteinBlockAligner(ProteinBlockAligner * mat) 
{
    ProteinBlockAligner_access_func_holder holder;   


    holder.access_main    = ProteinBlockAligner_shatter_access_main; 
    holder.access_special = ProteinBlockAligner_shatter_access_special;  
    assert(mat);     
    assert(mat->shatter);    
    return PackAln_read_generic_ProteinBlockAligner(mat,holder); 
}    


/* Function:  ProteinBlockAligner_shatter_access_main(mat,i,j,state)
 *
 * Descrip: No Description
 *
 * Arg:          mat [UNKN ] Undocumented argument [ProteinBlockAligner *]
 * Arg:            i [UNKN ] Undocumented argument [int]
 * Arg:            j [UNKN ] Undocumented argument [int]
 * Arg:        state [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int ProteinBlockAligner_shatter_access_main(ProteinBlockAligner * mat,int i,int j,int state) 
{
    return ProteinBlockAligner_SHATTER_MATRIX(mat,i,j,state);    
}    


/* Function:  ProteinBlockAligner_shatter_access_special(mat,i,j,state)
 *
 * Descrip: No Description
 *
 * Arg:          mat [UNKN ] Undocumented argument [ProteinBlockAligner *]
 * Arg:            i [UNKN ] Undocumented argument [int]
 * Arg:            j [UNKN ] Undocumented argument [int]
 * Arg:        state [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int ProteinBlockAligner_shatter_access_special(ProteinBlockAligner * mat,int i,int j,int state) 
{
    return ProteinBlockAligner_SHATTER_SPECIAL(mat,i,j,state);   
}    


/* Function:  calculate_shatter_ProteinBlockAligner(mat,dpenv)
 *
 * Descrip:    This function calculates the ProteinBlockAligner matrix when in shatter mode
 *
 *
 * Arg:          mat [UNKN ] (null) [ProteinBlockAligner *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean calculate_shatter_ProteinBlockAligner(ProteinBlockAligner * mat,DPEnvelope * dpenv) 
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


    mat->shatter = new_ShatterMatrix(dpenv,4,lenj,2);    
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


    start_reporting("ProteinBlockAligner Matrix calculation: "); 
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




        /* For state BLOCK_1 */ 
        /* setting first movement to score */ 
        score = SIG_1_1[BLOCK_1] + mat->b_self_trans;    
        /* From state START to state BLOCK_1 */ 
        temp = ProteinBlockAligner_SHATTER_SPECIAL(mat,i-1,j-1,START) + 0;   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state UNALIGNED to state BLOCK_1 */ 
        temp = SIG_1_1[UNALIGNED] + mat->bentry;     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for BLOCK_1 */ 
        /* Add any movement independant score and put away */ 
         score += CompMat_AAMATCH(mat->m,CSEQ_PROTEIN_AMINOACID(mat->q,i),CSEQ_PROTEIN_AMINOACID(mat->t,j));     
         SIG_0_0[BLOCK_1] = score;   


        /* state BLOCK_1 is a source for special END */ 
        temp = score + (0) + (0) ;   
        if( temp > ProteinBlockAligner_SHATTER_SPECIAL(mat,i,j,END) )    {  
          ProteinBlockAligner_SHATTER_SPECIAL(mat,i,j,END) = temp;   
          }  




        /* Finished calculating state BLOCK_1 */ 


        /* For state BLOCK_2 */ 
        /* setting first movement to score */ 
        score = SIG_1_1[BLOCK_2] + mat->b_self_trans;    
        /* From state BLOCK_1 to state BLOCK_2 */ 
        temp = SIG_1_1[BLOCK_1] + mat->bfor_trans;   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for BLOCK_2 */ 
        /* Add any movement independant score and put away */ 
         score += CompMat_AAMATCH(mat->m,CSEQ_PROTEIN_AMINOACID(mat->q,i),CSEQ_PROTEIN_AMINOACID(mat->t,j));     
         SIG_0_0[BLOCK_2] = score;   


        /* Finished calculating state BLOCK_2 */ 


        /* For state BLOCK_3 */ 
        /* setting first movement to score */ 
        score = SIG_1_1[BLOCK_3] + mat->b_self_trans;    
        /* From state BLOCK_2 to state BLOCK_3 */ 
        temp = SIG_1_1[BLOCK_2] + mat->bfor_trans;   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for BLOCK_3 */ 
        /* Add any movement independant score and put away */ 
         score += CompMat_AAMATCH(mat->m,CSEQ_PROTEIN_AMINOACID(mat->q,i),CSEQ_PROTEIN_AMINOACID(mat->t,j));     
         SIG_0_0[BLOCK_3] = score;   


        /* Finished calculating state BLOCK_3 */ 


        /* For state UNALIGNED */ 
        /* setting first movement to score */ 
        score = SIG_0_1[BLOCK_1] + mat->bexit;   
        /* From state BLOCK_1 to state UNALIGNED */ 
        temp = SIG_1_0[BLOCK_1] + mat->bexit;    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state BLOCK_2 to state UNALIGNED */ 
        temp = SIG_0_1[BLOCK_2] + mat->b3exit;   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state BLOCK_2 to state UNALIGNED */ 
        temp = SIG_1_0[BLOCK_2] + mat->b3exit;   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state UNALIGNED to state UNALIGNED */ 
        temp = SIG_0_1[UNALIGNED] + 0;   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state UNALIGNED to state UNALIGNED */ 
        temp = SIG_1_0[UNALIGNED] + 0;   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for UNALIGNED */ 
        /* Add any movement independant score and put away */ 
         SIG_0_0[UNALIGNED] = score; 


        /* state UNALIGNED is a source for special END */ 
        temp = score + (0) + (0) ;   
        if( temp > ProteinBlockAligner_SHATTER_SPECIAL(mat,i,j,END) )    {  
          ProteinBlockAligner_SHATTER_SPECIAL(mat,i,j,END) = temp;   
          }  




        /* Finished calculating state UNALIGNED */ 
        }  


      /* Special state START has no special to special movements */ 


      /* Special state END has no special to special movements */ 
      }  
    stop_reporting();    
    return TRUE;     
}    


/* Function:  search_ProteinBlockAligner(dbsi,out,querydb,targetdb,m,bentry,bexit,bfor_trans,b_self_trans,b3exit)
 *
 * Descrip:    This function makes a database search of ProteinBlockAligner
 *             It uses the dbsi structure to choose which implementation
 *             to use of the database searching. This way at run time you
 *             can switch between single threaded/multi-threaded or hardware
 *
 *
 * Arg:                dbsi [UNKN ] Undocumented argument [DBSearchImpl *]
 * Arg:                 out [UNKN ] Undocumented argument [Hscore *]
 * Arg:             querydb [UNKN ] Undocumented argument [ProteinDB*]
 * Arg:            targetdb [UNKN ] Undocumented argument [ProteinDB*]
 * Arg:                   m [UNKN ] Undocumented argument [CompMat*]
 * Arg:              bentry [UNKN ] Undocumented argument [Score]
 * Arg:               bexit [UNKN ] Undocumented argument [Score]
 * Arg:          bfor_trans [UNKN ] Undocumented argument [Score]
 * Arg:        b_self_trans [UNKN ] Undocumented argument [Score]
 * Arg:              b3exit [UNKN ] Undocumented argument [Score]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type search_ProteinBlockAligner(DBSearchImpl * dbsi,Hscore * out,ProteinDB* querydb,ProteinDB* targetdb ,CompMat* m,Score bentry,Score bexit,Score bfor_trans,Score b_self_trans,Score b3exit) 
{
#ifdef PTHREAD   
    int i;   
    int thr_no;  
    pthread_attr_t pat;  
    struct thread_pool_holder_ProteinBlockAligner * holder;  
#endif   
    if( out == NULL )    {  
      warn("Passed in a null Hscore object into search_ProteinBlockAligner. Can't process results!");    
      return SEARCH_ERROR;   
      }  
    if( dbsi == NULL )   {  
      warn("Passed in a null DBSearchImpl object into search_ProteinBlockAligner. Can't process results!");  
      return SEARCH_ERROR;   
      }  
    if( dbsi->trace_level > 5 )  
      warn("Asking for trace level of %d in database search for ProteinBlockAligner, but it was compiled with a trace level of -2139062144. Not all trace statements can be shown",dbsi->trace_level);   
    switch(dbsi->type)   { /*switch on implementation*/ 
      case DBSearchImpl_Serial : 
        return serial_search_ProteinBlockAligner(out,querydb, targetdb ,m,bentry,bexit,bfor_trans,b_self_trans,b3exit);  
      case DBSearchImpl_Pthreads :   
#ifdef PTHREAD   
        holder = (struct thread_pool_holder_ProteinBlockAligner *) ckalloc(sizeof(struct thread_pool_holder_ProteinBlockAligner));   
        if( holder == NULL )     {  
          warn("Unable to allocated thread pool datastructure...");  
          return SEARCH_ERROR;   
          }  
        holder->out = out;   
        holder->dbsi = dbsi; 
        holder->querydb = querydb;   
        holder->targetdb = targetdb; 
        holder->m = m;   
        holder->bentry = bentry; 
        holder->bexit = bexit;   
        holder->bfor_trans = bfor_trans; 
        holder->b_self_trans = b_self_trans; 
        holder->b3exit = b3exit; 
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
          if( pthread_create(holder->pool+i,&pat,thread_loop_ProteinBlockAligner,(void *)holder) )   
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
        warn("You did not specifiy the PTHREAD compile when compiled the C code for ProteinBlockAligner");   
#endif /* finished threads */    
      default :  
        warn("database search implementation %s was not provided in the compiled dynamite file from ProteinBlockAligner",impl_string_DBSearchImpl(dbsi));    
        return SEARCH_ERROR; 
      } /* end of switch on implementation */ 


}    


/* Function:  thread_loop_ProteinBlockAligner(ptr)
 *
 * Descrip:    Infinite loop code foreach thread for ProteinBlockAligner
 *
 *
 * Arg:        ptr [UNKN ] Undocumented argument [void *]
 *
 * Return [UNKN ]  Undocumented return value [void *]
 *
 */
#ifdef PTHREAD   
void * thread_loop_ProteinBlockAligner(void * ptr) 
{
    struct thread_pool_holder_ProteinBlockAligner * holder;  
    int db_status;   
    int score;   
    DataScore * ds;  
    ComplexSequence* q;  
    ComplexSequence* t;  


    holder = (struct thread_pool_holder_ProteinBlockAligner *) ptr;  
    if ( holder->dbsi->trace_level >= 1 )    
      fprintf(holder->dbsi->trace_file,"Entering infinite loop for thread...\n");    


    while(1) { /*Infinite loop over all models*/ 
      /* Get input lock */ 


      if ( holder->dbsi->trace_level >= 2 )  
        fprintf(holder->dbsi->trace_file,"About to get input lock for main reload\n");   
      if( pthread_mutex_lock(&(holder->input_lock))!= 0 )    
        fatal("Error on getting input lock for ProteinBlockAligner");    
      if ( holder->dbsi->trace_level >= 2 )  
        fprintf(holder->dbsi->trace_file,"Got input lock for main reload\n");    


      if( holder->search_has_ended == TRUE ) {  
        if ( holder->dbsi->trace_level >= 2 )    
          fprintf(holder->dbsi->trace_file,"Database search finished for me!...\n"); 
        if( pthread_mutex_unlock(&(holder->input_lock))!= 0 )    
          fatal("Error in releasing input lock for ProteinBlockAligner");    
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
        holder->q = init_ProteinDB(holder->querydb,&db_status);  
        holder->query_init = TRUE;   
        if( db_status == DB_RETURN_ERROR )   
          fatal("Unable to initalise query database in ProteinBlockAligner search"); 
        }  
      q = hard_link_ComplexSequence(holder->q);  
      /* get query information into datascore */ 
      dataentry_add_ProteinDB(ds->query,q,holder->querydb);  


      if( holder->target_init == FALSE ) { /*if the db has not been init'd*/ 
        if ( holder->dbsi->trace_level >= 3 )    
          fprintf(holder->dbsi->trace_file,"Starting target database...\n"); 
        t = init_ProteinDB(holder->targetdb,&db_status); 
        holder->target_init = TRUE;  
        } /* end of if the db has not been init'd */ 
      else   { /*Normal reload*/ 
         t = reload_ProteinDB(NULL,holder->targetdb,&db_status);     
        } /* end of Normal reload */ 


      /* Check to see what the reload is like */ 


      if( db_status == DB_RETURN_ERROR ) {  
        fatal("In searching ProteinBlockAligner, Reload error on database t, in threads");   
        }  


      if( db_status == DB_RETURN_END)    { /*End of target database*/ 
        /* close target database and schedule it for initalisation by next thread */ 
        close_ProteinDB(NULL,holder->targetdb);  
        holder->target_init = FALSE; 
        if ( holder->dbsi->trace_level >= 2 )    
          fprintf(holder->dbsi->trace_file,"Target Database to be reloaded...\n");   


        /* free'ing the query object */ 
        free_ComplexSequence(holder->q); 
        /* get the next query object for the next thread */ 
        holder->q = reload_ProteinDB(NULL,holder->querydb,&db_status);   
        if( db_status == DB_RETURN_ERROR )   
          fatal("In searching ProteinBlockAligner, reload error on database q, in threads"); 
        if( db_status == DB_RETURN_END ) { /*last load!*/ 
          /* End of target and query database - finished search! */ 
          close_ProteinDB(NULL,holder->querydb); 
          holder->search_has_ended = TRUE;   
          } /* end of last load! */ 


        /* release input mutex */ 
        if ( holder->dbsi->trace_level >= 2 )    
          fprintf(holder->dbsi->trace_file,"Releasing input lock after end of target\n");    
        if( pthread_mutex_unlock(&(holder->input_lock))!= 0 )    
          fatal("Error in releasing input lock for ProteinBlockAligner");    
        continue;    
        } /* end of End of target database */ 
      else   { /*Normal reload*/ 
        if ( holder->dbsi->trace_level >= 2 )    
          fprintf(holder->dbsi->trace_file,"Releasing input lock for normal reload\n");  
        if( pthread_mutex_unlock(&(holder->input_lock))!= 0 )    
          fatal("Error in releasing input lock for ProteinBlockAligner");    
        } /* end of Normal reload */ 
      /* get target information into datascore */ 
      dataentry_add_ProteinDB(ds->target,t,holder->targetdb);    


      /* Now there is a new query/target pair ready for comparison */ 
      if ( holder->dbsi->trace_level >= 1 )  
        fprintf(holder->dbsi->trace_file,"A new pair to be compared...\n");  
      score = score_only_ProteinBlockAligner(q, t ,holder->m,holder->bentry,holder->bexit,holder->bfor_trans,holder->b_self_trans,holder->b3exit);   


      if ( holder->dbsi->trace_level >= 2 )  
        fprintf(holder->dbsi->trace_file,"Getting output lock\n");   
      /* Getting lock on output */ 
      if( pthread_mutex_lock(&(holder->output_lock))!= 0 )   
        fatal("Error on getting output lock for ProteinBlockAligner");   
      /* If the score is less than cutoff, schedule the datascore for reuse */ 
      if( should_store_Hscore(holder->out,score) != TRUE)    {  
        free_DataScore(ds);  
        }  
      else   { /*storing score*/ 
        ds->score = score;   
        add_Hscore(holder->out,ds);  
        } /* end of storing score */ 
      if( pthread_mutex_unlock(&(holder->output_lock))!= 0 ) 
        fatal("Error on releasing output lock for ProteinBlockAligner"); 
      if ( holder->dbsi->trace_level >= 2 )  
        fprintf(holder->dbsi->trace_file,"Released output lock\n");  


      /* Now free database objects */ 
      if ( holder->dbsi->trace_level >= 2 )  
        fprintf(holder->dbsi->trace_file,"About to get input lock for free func\n"); 
      if( pthread_mutex_lock(&(holder->input_lock))!= 0 )    
        fatal("Error on getting input lock for ProteinBlockAligner");    
      if ( holder->dbsi->trace_level >= 2 )  
        fprintf(holder->dbsi->trace_file,"Got input lock for free func\n");  
      free_ComplexSequence(q);   
      free_ComplexSequence(t);   
      if ( holder->dbsi->trace_level >= 2 )  
        fprintf(holder->dbsi->trace_file,"Releasing input lock after free'ing\n");   
      if( pthread_mutex_unlock(&(holder->input_lock))!= 0 )  
        fatal("Error in releasing input lock for ProteinBlockAligner");  
      } /* end of Infinite loop over all models */ 


    if ( holder->dbsi->trace_level >= 1 )    
      fprintf(holder->dbsi->trace_file,"Exiting forever loop\n");    
    return NULL; 
}    


/* Function:  serial_search_ProteinBlockAligner(out,querydb,targetdb,m,bentry,bexit,bfor_trans,b_self_trans,b3exit)
 *
 * Descrip:    This function makes a database search of ProteinBlockAligner
 *             It is a single processor implementation
 *
 *
 * Arg:                 out [UNKN ] Undocumented argument [Hscore *]
 * Arg:             querydb [UNKN ] Undocumented argument [ProteinDB*]
 * Arg:            targetdb [UNKN ] Undocumented argument [ProteinDB*]
 * Arg:                   m [UNKN ] Undocumented argument [CompMat*]
 * Arg:              bentry [UNKN ] Undocumented argument [Score]
 * Arg:               bexit [UNKN ] Undocumented argument [Score]
 * Arg:          bfor_trans [UNKN ] Undocumented argument [Score]
 * Arg:        b_self_trans [UNKN ] Undocumented argument [Score]
 * Arg:              b3exit [UNKN ] Undocumented argument [Score]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
#endif /* PTHREAD */ 
Search_Return_Type serial_search_ProteinBlockAligner(Hscore * out,ProteinDB* querydb,ProteinDB* targetdb ,CompMat* m,Score bentry,Score bexit,Score bfor_trans,Score b_self_trans,Score b3exit) 
{
    ComplexSequence* q;  
    ComplexSequence* t;  
    int db_status;   
    int score;   
    int query_pos = 0;   
    int target_pos = 0;  
    DataScore * ds;  


    push_errormsg_stack("Before any actual search in db searching"); 
    q = init_ProteinDB(querydb,&db_status);  
    if( db_status == DB_RETURN_ERROR )   {  
      warn("In searching ProteinBlockAligner, got a database reload error on the query [q] database");   
      return SEARCH_ERROR;   
      }  
    for(;;)  { /*For all query entries*/ 


      target_pos = 0;    


      t = init_ProteinDB(targetdb,&db_status);   
      if( db_status == DB_RETURN_ERROR )     {  
        warn("In searching ProteinBlockAligner, got a database init error on the target [t] database");  
        return SEARCH_ERROR; 
        }  
      for(;;)    { /*For all target entries*/ 


        /* No maximum length - allocated on-the-fly */ 
        score = score_only_ProteinBlockAligner(q, t , m, bentry, bexit, bfor_trans, b_self_trans, b3exit);   
        if( should_store_Hscore(out,score) == TRUE )     { /*if storing datascore*/ 
          ds = new_DataScore_from_storage(out);  
          if( ds == NULL )   {  
            warn("ProteinBlockAligner search had a memory error in allocating a new_DataScore (?a leak somewhere - DataScore is a very small datastructure");    
            return SEARCH_ERROR; 
            }  
          /* Now: add query/target information to the entry */ 
          dataentry_add_ProteinDB(ds->query,q,querydb);  
          dataentry_add_ProteinDB(ds->target,t,targetdb);    
          ds->score = score;     
          add_Hscore(out,ds);    
          } /* end of if storing datascore */ 
        pop_errormsg_stack();    
        push_errormsg_stack("DB searching: just finished [Query Pos: %d] [Target Pos: %d]",query_pos,target_pos);    


         t = reload_ProteinDB(t,targetdb,&db_status);    
        if( db_status == DB_RETURN_ERROR )   {  
          warn("In searching ProteinBlockAligner, Reload error on database t, position %d,%d",query_pos,target_pos); 
          return SEARCH_ERROR;   
          }  
        if( db_status == DB_RETURN_END ) 
          break;/* Out of target loop */ 
        target_pos++;    
        } /* end of For all target entries */ 
      close_ProteinDB(t,targetdb);   
       q = reload_ProteinDB(q,querydb,&db_status);   
      if( db_status == DB_RETURN_ERROR)  {  
        warn("In searching ProteinBlockAligner, Reload error on database q, position %d,%d",query_pos,target_pos);   
        return SEARCH_ERROR; 
        }  
      if( db_status == DB_RETURN_END)    
        break;  /* Out of query loop */ 
      query_pos++;   
      } /* end of For all query entries */ 
    close_ProteinDB(q,querydb);  
    pop_errormsg_stack();    
    return SEARCH_OK;    
}    


/* Function:  score_only_ProteinBlockAligner(q,t,m,bentry,bexit,bfor_trans,b_self_trans,b3exit)
 *
 * Descrip:    This function just calculates the score for the matrix
 *             I am pretty sure we can do this better, but hey, for the moment...
 *             It calls /allocate_ProteinBlockAligner_only
 *
 *
 * Arg:                   q [UNKN ] query data structure [ComplexSequence*]
 * Arg:                   t [UNKN ] target data structure [ComplexSequence*]
 * Arg:                   m [UNKN ] Resource [CompMat*]
 * Arg:              bentry [UNKN ] Resource [Score]
 * Arg:               bexit [UNKN ] Resource [Score]
 * Arg:          bfor_trans [UNKN ] Resource [Score]
 * Arg:        b_self_trans [UNKN ] Resource [Score]
 * Arg:              b3exit [UNKN ] Resource [Score]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int score_only_ProteinBlockAligner(ComplexSequence* q,ComplexSequence* t ,CompMat* m,Score bentry,Score bexit,Score bfor_trans,Score b_self_trans,Score b3exit) 
{
    int bestscore = NEGI;    
    int i;   
    int j;   
    int k;   
    ProteinBlockAligner * mat;   


    mat = allocate_ProteinBlockAligner_only(q, t , m, bentry, bexit, bfor_trans, b_self_trans, b3exit);  
    if( mat == NULL )    {  
      warn("Memory allocation error in the db search - unable to communicate to calling function. this spells DIASTER!");    
      return NEGI;   
      }  
    if((mat->basematrix = BaseMatrix_alloc_matrix_and_specials(2,(mat->leni + 1) * 4,2,2)) == NULL)  {  
      warn("Score only matrix for ProteinBlockAligner cannot be allocated, (asking for 1  by %d  cells)",mat->leni*4);   
      mat = free_ProteinBlockAligner(mat);   
      return 0;  
      }  
    mat->basematrix->type = BASEMATRIX_TYPE_VERYSMALL;   


    /* Now, initiate matrix */ 
    for(j=0;j<3;j++) {  
      for(i=(-1);i<mat->leni;i++)    {  
        for(k=0;k<4;k++) 
          ProteinBlockAligner_VSMALL_MATRIX(mat,i,j,k) = NEGI;   
        }  
      ProteinBlockAligner_VSMALL_SPECIAL(mat,i,j,START) = 0; 
      ProteinBlockAligner_VSMALL_SPECIAL(mat,i,j,END) = NEGI;    
      }  


    /* Ok, lets do-o-o-o-o it */ 


    for(j=0;j<mat->lenj;j++) { /*for all target positions*/ 
      auto int score;    
      auto int temp;     
      for(i=0;i<mat->leni;i++)   { /*for all query positions*/ 


        /* For state BLOCK_1 */ 
        /* setting first movement to score */ 
        score = ProteinBlockAligner_VSMALL_MATRIX(mat,i-1,j-1,BLOCK_1) + mat->b_self_trans;  
        /* From state START to state BLOCK_1 */ 
        temp = ProteinBlockAligner_VSMALL_SPECIAL(mat,i-1,j-1,START) + 0;    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state UNALIGNED to state BLOCK_1 */ 
        temp = ProteinBlockAligner_VSMALL_MATRIX(mat,i-1,j-1,UNALIGNED) + mat->bentry;   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for BLOCK_1 */ 
        /* Add any movement independant score and put away */ 
         score += CompMat_AAMATCH(mat->m,CSEQ_PROTEIN_AMINOACID(mat->q,i),CSEQ_PROTEIN_AMINOACID(mat->t,j));     
         ProteinBlockAligner_VSMALL_MATRIX(mat,i,j,BLOCK_1) = score; 


        /* state BLOCK_1 is a source for special END */ 
        temp = score + (0) + (0) ;   
        if( temp > ProteinBlockAligner_VSMALL_SPECIAL(mat,i,j,END) )     {  
          ProteinBlockAligner_VSMALL_SPECIAL(mat,i,j,END) = temp;    
          }  




        /* Finished calculating state BLOCK_1 */ 


        /* For state BLOCK_2 */ 
        /* setting first movement to score */ 
        score = ProteinBlockAligner_VSMALL_MATRIX(mat,i-1,j-1,BLOCK_2) + mat->b_self_trans;  
        /* From state BLOCK_1 to state BLOCK_2 */ 
        temp = ProteinBlockAligner_VSMALL_MATRIX(mat,i-1,j-1,BLOCK_1) + mat->bfor_trans;     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for BLOCK_2 */ 
        /* Add any movement independant score and put away */ 
         score += CompMat_AAMATCH(mat->m,CSEQ_PROTEIN_AMINOACID(mat->q,i),CSEQ_PROTEIN_AMINOACID(mat->t,j));     
         ProteinBlockAligner_VSMALL_MATRIX(mat,i,j,BLOCK_2) = score; 


        /* Finished calculating state BLOCK_2 */ 


        /* For state BLOCK_3 */ 
        /* setting first movement to score */ 
        score = ProteinBlockAligner_VSMALL_MATRIX(mat,i-1,j-1,BLOCK_3) + mat->b_self_trans;  
        /* From state BLOCK_2 to state BLOCK_3 */ 
        temp = ProteinBlockAligner_VSMALL_MATRIX(mat,i-1,j-1,BLOCK_2) + mat->bfor_trans;     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for BLOCK_3 */ 
        /* Add any movement independant score and put away */ 
         score += CompMat_AAMATCH(mat->m,CSEQ_PROTEIN_AMINOACID(mat->q,i),CSEQ_PROTEIN_AMINOACID(mat->t,j));     
         ProteinBlockAligner_VSMALL_MATRIX(mat,i,j,BLOCK_3) = score; 


        /* Finished calculating state BLOCK_3 */ 


        /* For state UNALIGNED */ 
        /* setting first movement to score */ 
        score = ProteinBlockAligner_VSMALL_MATRIX(mat,i-0,j-1,BLOCK_1) + mat->bexit;     
        /* From state BLOCK_1 to state UNALIGNED */ 
        temp = ProteinBlockAligner_VSMALL_MATRIX(mat,i-1,j-0,BLOCK_1) + mat->bexit;  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state BLOCK_2 to state UNALIGNED */ 
        temp = ProteinBlockAligner_VSMALL_MATRIX(mat,i-0,j-1,BLOCK_2) + mat->b3exit;     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state BLOCK_2 to state UNALIGNED */ 
        temp = ProteinBlockAligner_VSMALL_MATRIX(mat,i-1,j-0,BLOCK_2) + mat->b3exit;     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state UNALIGNED to state UNALIGNED */ 
        temp = ProteinBlockAligner_VSMALL_MATRIX(mat,i-0,j-1,UNALIGNED) + 0;     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state UNALIGNED to state UNALIGNED */ 
        temp = ProteinBlockAligner_VSMALL_MATRIX(mat,i-1,j-0,UNALIGNED) + 0;     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for UNALIGNED */ 
        /* Add any movement independant score and put away */ 
         ProteinBlockAligner_VSMALL_MATRIX(mat,i,j,UNALIGNED) = score;   


        /* state UNALIGNED is a source for special END */ 
        temp = score + (0) + (0) ;   
        if( temp > ProteinBlockAligner_VSMALL_SPECIAL(mat,i,j,END) )     {  
          ProteinBlockAligner_VSMALL_SPECIAL(mat,i,j,END) = temp;    
          }  




        /* Finished calculating state UNALIGNED */ 
        } /* end of for all query positions */ 




      /* Special state START has no special to special movements */ 


      /* Special state END has no special to special movements */ 
      if( bestscore < ProteinBlockAligner_VSMALL_SPECIAL(mat,0,j,END) )  
        bestscore = ProteinBlockAligner_VSMALL_SPECIAL(mat,0,j,END);     
      } /* end of for all target positions */ 


    mat = free_ProteinBlockAligner(mat);     
    return bestscore;    
}    


/* Function:  PackAln_bestmemory_ProteinBlockAligner(q,t,m,bentry,bexit,bfor_trans,b_self_trans,b3exit,dpenv,dpri)
 *
 * Descrip:    This function chooses the best memory set-up for the alignment
 *             using calls to basematrix, and then implements either a large
 *             or small memory model.
 *
 *             It is the best function to use if you just want an alignment
 *
 *             If you want a label alignment, you will need
 *             /convert_PackAln_to_AlnBlock_ProteinBlockAligner
 *
 *
 * Arg:                   q [UNKN ] query data structure [ComplexSequence*]
 * Arg:                   t [UNKN ] target data structure [ComplexSequence*]
 * Arg:                   m [UNKN ] Resource [CompMat*]
 * Arg:              bentry [UNKN ] Resource [Score]
 * Arg:               bexit [UNKN ] Resource [Score]
 * Arg:          bfor_trans [UNKN ] Resource [Score]
 * Arg:        b_self_trans [UNKN ] Resource [Score]
 * Arg:              b3exit [UNKN ] Resource [Score]
 * Arg:               dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:                dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * PackAln_bestmemory_ProteinBlockAligner(ComplexSequence* q,ComplexSequence* t ,CompMat* m,Score bentry,Score bexit,Score bfor_trans,Score b_self_trans,Score b3exit,DPEnvelope * dpenv,DPRunImpl * dpri) 
{
    long long total; 
    ProteinBlockAligner * mat;   
    PackAln * out;   
    DebugMatrix * de;    
    DPRunImplMemory strategy;    
    assert(dpri);    


    total = q->seq->len * t->seq->len;   
    if( dpri->memory == DPIM_Default )   {  
      if( (total * 4 * sizeof(int)) > 1000*dpri->kbyte_size) {  
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
        if( (mat=allocate_Expl_ProteinBlockAligner(q, t , m, bentry, bexit, bfor_trans, b_self_trans, b3exit,dpri)) == NULL )    {  
          warn("Unable to allocate large ProteinBlockAligner version");  
          return NULL;   
          }  
        calculate_dpenv_ProteinBlockAligner(mat,dpenv);  
        out =  PackAln_read_Expl_ProteinBlockAligner(mat);   
        }  
      else   {  
        mat = allocate_ProteinBlockAligner_only(q, t , m, bentry, bexit, bfor_trans, b_self_trans, b3exit);  
        calculate_shatter_ProteinBlockAligner(mat,dpenv);    
        out = PackAln_read_Shatter_ProteinBlockAligner(mat);     
        }  
      }  
    else {  
      if( strategy == DPIM_Linear )  {  
        /* use small implementation */ 
        if( (mat=allocate_Small_ProteinBlockAligner(q, t , m, bentry, bexit, bfor_trans, b_self_trans, b3exit)) == NULL )    {  
          warn("Unable to allocate small ProteinBlockAligner version");  
          return NULL;   
          }  
        out = PackAln_calculate_Small_ProteinBlockAligner(mat,dpenv);    
        }  
      else   {  
        /* use Large implementation */ 
        if( (mat=allocate_Expl_ProteinBlockAligner(q, t , m, bentry, bexit, bfor_trans, b_self_trans, b3exit,dpri)) == NULL )    {  
          warn("Unable to allocate large ProteinBlockAligner version");  
          return NULL;   
          }  
        if( dpri->debug == TRUE) {  
          fatal("Asked for dydebug, but dynamite file not compiled with -g. Need to recompile dynamite source"); 
          }  
        else {  
          calculate_ProteinBlockAligner(mat);    
          out =  PackAln_read_Expl_ProteinBlockAligner(mat); 
          }  
        }  
      }  


    mat = free_ProteinBlockAligner(mat);     
    return out;  
}    


/* Function:  allocate_ProteinBlockAligner_only(q,t,m,bentry,bexit,bfor_trans,b_self_trans,b3exit)
 *
 * Descrip:    This function only allocates the ProteinBlockAligner structure
 *             checks types where possible and determines leni and lenj
 *             The basematrix area is delt with elsewhere
 *
 *
 * Arg:                   q [UNKN ] query data structure [ComplexSequence*]
 * Arg:                   t [UNKN ] target data structure [ComplexSequence*]
 * Arg:                   m [UNKN ] Resource [CompMat*]
 * Arg:              bentry [UNKN ] Resource [Score]
 * Arg:               bexit [UNKN ] Resource [Score]
 * Arg:          bfor_trans [UNKN ] Resource [Score]
 * Arg:        b_self_trans [UNKN ] Resource [Score]
 * Arg:              b3exit [UNKN ] Resource [Score]
 *
 * Return [UNKN ]  Undocumented return value [ProteinBlockAligner *]
 *
 */
ProteinBlockAligner * allocate_ProteinBlockAligner_only(ComplexSequence* q,ComplexSequence* t ,CompMat* m,Score bentry,Score bexit,Score bfor_trans,Score b_self_trans,Score b3exit) 
{
    ProteinBlockAligner * out;   


    if((out= ProteinBlockAligner_alloc()) == NULL)   {  
      warn("Allocation of basic ProteinBlockAligner structure failed...");   
      return NULL;   
      }  


    out->q = q;  
    out->t = t;  
    out->m = m;  
    out->bentry = bentry;    
    out->bexit = bexit;  
    out->bfor_trans = bfor_trans;    
    out->b_self_trans = b_self_trans;    
    out->b3exit = b3exit;    
    out->leni = q->seq->len;     
    out->lenj = t->seq->len;     
    return out;  
}    


/* Function:  allocate_Expl_ProteinBlockAligner(q,t,m,bentry,bexit,bfor_trans,b_self_trans,b3exit,dpri)
 *
 * Descrip:    This function allocates the ProteinBlockAligner structure
 *             and the basematrix area for explicit memory implementations
 *             It calls /allocate_ProteinBlockAligner_only
 *
 *
 * Arg:                   q [UNKN ] query data structure [ComplexSequence*]
 * Arg:                   t [UNKN ] target data structure [ComplexSequence*]
 * Arg:                   m [UNKN ] Resource [CompMat*]
 * Arg:              bentry [UNKN ] Resource [Score]
 * Arg:               bexit [UNKN ] Resource [Score]
 * Arg:          bfor_trans [UNKN ] Resource [Score]
 * Arg:        b_self_trans [UNKN ] Resource [Score]
 * Arg:              b3exit [UNKN ] Resource [Score]
 * Arg:                dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [ProteinBlockAligner *]
 *
 */
ProteinBlockAligner * allocate_Expl_ProteinBlockAligner(ComplexSequence* q,ComplexSequence* t ,CompMat* m,Score bentry,Score bexit,Score bfor_trans,Score b_self_trans,Score b3exit,DPRunImpl * dpri) 
{
    ProteinBlockAligner * out;   


    out = allocate_ProteinBlockAligner_only(q, t , m, bentry, bexit, bfor_trans, b_self_trans, b3exit);  
    if( out == NULL )    
      return NULL;   
    if( dpri->should_cache == TRUE ) {  
      if( dpri->cache != NULL )  {  
        if( dpri->cache->maxleni >= (out->lenj+1)*4 && dpri->cache->maxlenj >= (out->leni+1))    
          out->basematrix = hard_link_BaseMatrix(dpri->cache);   
        else 
          dpri->cache = free_BaseMatrix(dpri->cache);    
        }  
      }  
    if( out->basematrix == NULL )    {  
      if( (out->basematrix = BaseMatrix_alloc_matrix_and_specials((out->lenj+1)*4,(out->leni+1),2,out->lenj+1)) == NULL) {  
        warn("Explicit matrix ProteinBlockAligner cannot be allocated, (asking for %d by %d main cells)",out->leni,out->lenj);   
        free_ProteinBlockAligner(out);   
        return NULL; 
        }  
      }  
    if( dpri->should_cache == TRUE && dpri->cache == NULL)   
      dpri->cache = hard_link_BaseMatrix(out->basematrix);   
    out->basematrix->type = BASEMATRIX_TYPE_EXPLICIT;    
    init_ProteinBlockAligner(out);   
    return out;  
}    


/* Function:  init_ProteinBlockAligner(mat)
 *
 * Descrip:    This function initates ProteinBlockAligner matrix when in explicit mode
 *             Called in /allocate_Expl_ProteinBlockAligner
 *
 *
 * Arg:        mat [UNKN ] ProteinBlockAligner which contains explicit basematrix memory [ProteinBlockAligner *]
 *
 */
void init_ProteinBlockAligner(ProteinBlockAligner * mat) 
{
    register int i;  
    register int j;  
    if( mat->basematrix->type != BASEMATRIX_TYPE_EXPLICIT)   {  
      warn("Cannot iniate matrix, is not an explicit memory type and you have assummed that");   
      return;    
      }  


    for(i= (-1);i<mat->q->seq->len;i++)  {  
      for(j= (-1);j<2;j++)   {  
        ProteinBlockAligner_EXPL_MATRIX(mat,i,j,BLOCK_1) = NEGI; 
        ProteinBlockAligner_EXPL_MATRIX(mat,i,j,BLOCK_2) = NEGI; 
        ProteinBlockAligner_EXPL_MATRIX(mat,i,j,BLOCK_3) = NEGI; 
        ProteinBlockAligner_EXPL_MATRIX(mat,i,j,UNALIGNED) = NEGI;   
        }  
      }  
    for(j= (-1);j<mat->t->seq->len;j++)  {  
      for(i= (-1);i<2;i++)   {  
        ProteinBlockAligner_EXPL_MATRIX(mat,i,j,BLOCK_1) = NEGI; 
        ProteinBlockAligner_EXPL_MATRIX(mat,i,j,BLOCK_2) = NEGI; 
        ProteinBlockAligner_EXPL_MATRIX(mat,i,j,BLOCK_3) = NEGI; 
        ProteinBlockAligner_EXPL_MATRIX(mat,i,j,UNALIGNED) = NEGI;   
        }  
      ProteinBlockAligner_EXPL_SPECIAL(mat,i,j,START) = 0;   
      ProteinBlockAligner_EXPL_SPECIAL(mat,i,j,END) = NEGI;  
      }  
    return;  
}    


/* Function:  recalculate_PackAln_ProteinBlockAligner(pal,mat)
 *
 * Descrip:    This function recalculates the PackAln structure produced by ProteinBlockAligner
 *             For example, in linear space methods this is used to score them
 *
 *
 * Arg:        pal [UNKN ] Undocumented argument [PackAln *]
 * Arg:        mat [UNKN ] Undocumented argument [ProteinBlockAligner *]
 *
 */
void recalculate_PackAln_ProteinBlockAligner(PackAln * pal,ProteinBlockAligner * mat) 
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
        case BLOCK_1 :   
          if( offi == 1 && offj == 1 && prev->state == BLOCK_1 ) {  
            pau->score = mat->b_self_trans + (CompMat_AAMATCH(mat->m,CSEQ_PROTEIN_AMINOACID(mat->q,i),CSEQ_PROTEIN_AMINOACID(mat->t,j)));    
            continue;    
            }  
          if( offj == 1 && prev->state == (START+4) )    {  
            pau->score = 0 + (CompMat_AAMATCH(mat->m,CSEQ_PROTEIN_AMINOACID(mat->q,i),CSEQ_PROTEIN_AMINOACID(mat->t,j)));    
            continue;    
            }  
          if( offi == 1 && offj == 1 && prev->state == UNALIGNED )   {  
            pau->score = mat->bentry + (CompMat_AAMATCH(mat->m,CSEQ_PROTEIN_AMINOACID(mat->q,i),CSEQ_PROTEIN_AMINOACID(mat->t,j)));  
            continue;    
            }  
          warn("In recaluclating PackAln with state BLOCK_1, from [%d,%d,%d], got a bad source state. Error!",offi,offj,prev->state);    
          break; 
        case BLOCK_2 :   
          if( offi == 1 && offj == 1 && prev->state == BLOCK_2 ) {  
            pau->score = mat->b_self_trans + (CompMat_AAMATCH(mat->m,CSEQ_PROTEIN_AMINOACID(mat->q,i),CSEQ_PROTEIN_AMINOACID(mat->t,j)));    
            continue;    
            }  
          if( offi == 1 && offj == 1 && prev->state == BLOCK_1 ) {  
            pau->score = mat->bfor_trans + (CompMat_AAMATCH(mat->m,CSEQ_PROTEIN_AMINOACID(mat->q,i),CSEQ_PROTEIN_AMINOACID(mat->t,j)));  
            continue;    
            }  
          warn("In recaluclating PackAln with state BLOCK_2, from [%d,%d,%d], got a bad source state. Error!",offi,offj,prev->state);    
          break; 
        case BLOCK_3 :   
          if( offi == 1 && offj == 1 && prev->state == BLOCK_3 ) {  
            pau->score = mat->b_self_trans + (CompMat_AAMATCH(mat->m,CSEQ_PROTEIN_AMINOACID(mat->q,i),CSEQ_PROTEIN_AMINOACID(mat->t,j)));    
            continue;    
            }  
          if( offi == 1 && offj == 1 && prev->state == BLOCK_2 ) {  
            pau->score = mat->bfor_trans + (CompMat_AAMATCH(mat->m,CSEQ_PROTEIN_AMINOACID(mat->q,i),CSEQ_PROTEIN_AMINOACID(mat->t,j)));  
            continue;    
            }  
          warn("In recaluclating PackAln with state BLOCK_3, from [%d,%d,%d], got a bad source state. Error!",offi,offj,prev->state);    
          break; 
        case UNALIGNED :     
          if( offi == 0 && offj == 1 && prev->state == BLOCK_1 ) {  
            pau->score = mat->bexit + (0);   
            continue;    
            }  
          if( offi == 1 && offj == 0 && prev->state == BLOCK_1 ) {  
            pau->score = mat->bexit + (0);   
            continue;    
            }  
          if( offi == 0 && offj == 1 && prev->state == BLOCK_2 ) {  
            pau->score = mat->b3exit + (0);  
            continue;    
            }  
          if( offi == 1 && offj == 0 && prev->state == BLOCK_2 ) {  
            pau->score = mat->b3exit + (0);  
            continue;    
            }  
          if( offi == 0 && offj == 1 && prev->state == UNALIGNED )   {  
            pau->score = 0 + (0);    
            continue;    
            }  
          if( offi == 1 && offj == 0 && prev->state == UNALIGNED )   {  
            pau->score = 0 + (0);    
            continue;    
            }  
          warn("In recaluclating PackAln with state UNALIGNED, from [%d,%d,%d], got a bad source state. Error!",offi,offj,prev->state);  
          break; 
        case (START+4) :     
          warn("In recaluclating PackAln with state START, got a bad source state. Error!"); 
          break; 
        case (END+4) :   
          if( offj == 0 && prev->state == UNALIGNED )    {  
            /* i here comes from the previous state ;) - not the real one */ 
            i = prev->i; 
            pau->score = 0 + (0);    
            continue;    
            }  
          if( offj == 0 && prev->state == BLOCK_1 )  {  
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
#define ProteinBlockAligner_HIDDEN_MATRIX(thismatrix,i,j,state) (thismatrix->basematrix->matrix[(j-hiddenj+1)][(i+1)*4+state])   
#define ProteinBlockAligner_DC_SHADOW_MATRIX(thismatrix,i,j,state) (thismatrix->basematrix->matrix[((j+2)*8) % 16][(i+1)*4+state])   
#define ProteinBlockAligner_HIDDEN_SPECIAL(thismatrix,i,j,state) (thismatrix->basematrix->specmatrix[state][(j+1)])  
#define ProteinBlockAligner_DC_SHADOW_SPECIAL(thismatrix,i,j,state) (thismatrix->basematrix->specmatrix[state*8][(j+1)]) 
#define ProteinBlockAligner_DC_SHADOW_MATRIX_SP(thismatrix,i,j,state,shadow) (thismatrix->basematrix->matrix[((((j+2)*8)+(shadow+1)) % 16)][(i+1)*4 + state])    
#define ProteinBlockAligner_DC_SHADOW_SPECIAL_SP(thismatrix,i,j,state,shadow) (thismatrix->basematrix->specmatrix[state*8 +shadow+1][(j+1)]) 
#define ProteinBlockAligner_DC_OPT_SHADOW_MATRIX(thismatrix,i,j,state) (score_pointers[(((j+1)% 1) * (leni+1) * 4) + ((i+1) * 4) + (state)]) 
#define ProteinBlockAligner_DC_OPT_SHADOW_MATRIX_SP(thismatrix,i,j,state,shadow) (shadow_pointers[(((j+1)% 1) * (leni+1) * 32) + ((i+1) * 32) + (state * 8) + shadow+1]) 
#define ProteinBlockAligner_DC_OPT_SHADOW_SPECIAL(thismatrix,i,j,state) (thismatrix->basematrix->specmatrix[state*8][(j+1)]) 
/* Function:  allocate_Small_ProteinBlockAligner(q,t,m,bentry,bexit,bfor_trans,b_self_trans,b3exit)
 *
 * Descrip:    This function allocates the ProteinBlockAligner structure
 *             and the basematrix area for a small memory implementations
 *             It calls /allocate_ProteinBlockAligner_only
 *
 *
 * Arg:                   q [UNKN ] query data structure [ComplexSequence*]
 * Arg:                   t [UNKN ] target data structure [ComplexSequence*]
 * Arg:                   m [UNKN ] Resource [CompMat*]
 * Arg:              bentry [UNKN ] Resource [Score]
 * Arg:               bexit [UNKN ] Resource [Score]
 * Arg:          bfor_trans [UNKN ] Resource [Score]
 * Arg:        b_self_trans [UNKN ] Resource [Score]
 * Arg:              b3exit [UNKN ] Resource [Score]
 *
 * Return [UNKN ]  Undocumented return value [ProteinBlockAligner *]
 *
 */
#define ProteinBlockAligner_DC_OPT_SHADOW_SPECIAL_SP(thismatrix,i,j,state,shadow) (thismatrix->basematrix->specmatrix[state*8 +shadow+1][(j+1)]) 
ProteinBlockAligner * allocate_Small_ProteinBlockAligner(ComplexSequence* q,ComplexSequence* t ,CompMat* m,Score bentry,Score bexit,Score bfor_trans,Score b_self_trans,Score b3exit) 
{
    ProteinBlockAligner * out;   


    out = allocate_ProteinBlockAligner_only(q, t , m, bentry, bexit, bfor_trans, b_self_trans, b3exit);  
    if( out == NULL )    
      return NULL;   
    out->basematrix = BaseMatrix_alloc_matrix_and_specials(16,(out->leni + 1) * 4,16,out->lenj+1);   
    if(out == NULL)  {  
      warn("Small shadow matrix ProteinBlockAligner cannot be allocated, (asking for 2 by %d main cells)",out->leni+2);  
      free_ProteinBlockAligner(out);     
      return NULL;   
      }  
    out->basematrix->type = BASEMATRIX_TYPE_SHADOW;  
    return out;  
}    


/* Function:  PackAln_calculate_Small_ProteinBlockAligner(mat,dpenv)
 *
 * Descrip:    This function calculates an alignment for ProteinBlockAligner structure in linear space
 *             If you want only the start/end points
 *             use /AlnRangeSet_calculate_Small_ProteinBlockAligner 
 *
 *             The function basically
 *               finds start/end points 
 *               foreach start/end point 
 *                 calls /full_dc_ProteinBlockAligner 
 *
 *
 * Arg:          mat [UNKN ] Undocumented argument [ProteinBlockAligner *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * PackAln_calculate_Small_ProteinBlockAligner(ProteinBlockAligner * mat,DPEnvelope * dpenv) 
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
      warn("Could not calculate packaln small for ProteinBlockAligner due to wrong type of matrix"); 
      return NULL;   
      }  


    out = PackAln_alloc_std();   


    start_reporting("Find start end points: ");  
    dc_optimised_start_end_calc_ProteinBlockAligner(mat,dpenv);  
    score = start_end_find_end_ProteinBlockAligner(mat,&endj);   
    out->score = score;  
    stopstate = END;
    
    /* No special to specials: one matrix alignment: simply remove and get */ 
    starti = ProteinBlockAligner_DC_SHADOW_SPECIAL_SP(mat,0,endj,END,0); 
    startj = ProteinBlockAligner_DC_SHADOW_SPECIAL_SP(mat,0,endj,END,1); 
    startstate = ProteinBlockAligner_DC_SHADOW_SPECIAL_SP(mat,0,endj,END,2); 
    stopi = ProteinBlockAligner_DC_SHADOW_SPECIAL_SP(mat,0,endj,END,3);  
    stopj = ProteinBlockAligner_DC_SHADOW_SPECIAL_SP(mat,0,endj,END,4);  
    stopstate = ProteinBlockAligner_DC_SHADOW_SPECIAL_SP(mat,0,endj,END,5);  
    temp = ProteinBlockAligner_DC_SHADOW_SPECIAL_SP(mat,0,endj,END,6);   
    log_full_error(REPORT,0,"[%d,%d][%d,%d] Score %d",starti,startj,stopi,stopj,score);  
    stop_reporting();    
    start_reporting("Recovering alignment: ");   


    /* Figuring how much j we have to align for reporting purposes */ 
    donej = 0;   
    totalj = stopj - startj; 
    full_dc_ProteinBlockAligner(mat,starti,startj,startstate,stopi,stopj,stopstate,out,&donej,totalj,dpenv); 


    /* Although we have no specials, need to get start. Better to check than assume */ 


    max_matrix_to_special_ProteinBlockAligner(mat,starti,startj,startstate,temp,&stopi,&stopj,&stopstate,&temp,NULL);    
    if( stopi == ProteinBlockAligner_READ_OFF_ERROR || stopstate != START )  {  
      warn("Problem in reading off special state system, hit a non start state (or an internal error) in a single alignment mode");  
      invert_PackAln(out);   
      recalculate_PackAln_ProteinBlockAligner(out,mat);  
      return out;    
      }  


    /* Ok. Put away start start... */ 
    pau = PackAlnUnit_alloc();   
    pau->i = stopi;  
    pau->j = stopj;  
    pau->state = stopstate + 4;  
    add_PackAln(out,pau);    


    log_full_error(REPORT,0,"Alignment recovered");  
    stop_reporting();    
    invert_PackAln(out); 
    recalculate_PackAln_ProteinBlockAligner(out,mat);    
    return out;  


}    


/* Function:  AlnRangeSet_calculate_Small_ProteinBlockAligner(mat)
 *
 * Descrip:    This function calculates an alignment for ProteinBlockAligner structure in linear space
 *             If you want the full alignment, use /PackAln_calculate_Small_ProteinBlockAligner 
 *             If you have already got the full alignment, but want the range set, use /AlnRangeSet_from_PackAln_ProteinBlockAligner
 *             If you have got the small matrix but not the alignment, use /AlnRangeSet_from_ProteinBlockAligner 
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [ProteinBlockAligner *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * AlnRangeSet_calculate_Small_ProteinBlockAligner(ProteinBlockAligner * mat) 
{
    AlnRangeSet * out;   


    start_reporting("Find start end points: ");  
    dc_optimised_start_end_calc_ProteinBlockAligner(mat,NULL);   
    log_full_error(REPORT,0,"Calculated");   


    out = AlnRangeSet_from_ProteinBlockAligner(mat); 
    return out;  
}    


/* Function:  AlnRangeSet_from_ProteinBlockAligner(mat)
 *
 * Descrip:    This function reads off a start/end structure
 *             for ProteinBlockAligner structure in linear space
 *             If you want the full alignment use
 *             /PackAln_calculate_Small_ProteinBlockAligner 
 *             If you have not calculated the matrix use
 *             /AlnRange_calculate_Small_ProteinBlockAligner
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [ProteinBlockAligner *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * AlnRangeSet_from_ProteinBlockAligner(ProteinBlockAligner * mat) 
{
    AlnRangeSet * out;   
    AlnRange * temp; 
    int jpos;    
    int state;   


    if( mat->basematrix->type != BASEMATRIX_TYPE_SHADOW) {  
      warn("Bad error! - non shadow matrix type in AlnRangeSet_from_ProteinBlockAligner");   
      return NULL;   
      }  


    out = AlnRangeSet_alloc_std();   
    /* Find the end position */ 
    out->score = start_end_find_end_ProteinBlockAligner(mat,&jpos);  
    state = END; 


    while( (temp = AlnRange_build_ProteinBlockAligner(mat,jpos,state,&jpos,&state)) != NULL) 
      add_AlnRangeSet(out,temp); 
    return out;  
}    


/* Function:  AlnRange_build_ProteinBlockAligner(mat,stopj,stopspecstate,startj,startspecstate)
 *
 * Descrip:    This function calculates a single start/end set in linear space
 *             Really a sub-routine for /AlnRangeSet_from_PackAln_ProteinBlockAligner
 *
 *
 * Arg:                   mat [UNKN ] Undocumented argument [ProteinBlockAligner *]
 * Arg:                 stopj [UNKN ] Undocumented argument [int]
 * Arg:         stopspecstate [UNKN ] Undocumented argument [int]
 * Arg:                startj [UNKN ] Undocumented argument [int *]
 * Arg:        startspecstate [UNKN ] Undocumented argument [int *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRange *]
 *
 */
AlnRange * AlnRange_build_ProteinBlockAligner(ProteinBlockAligner * mat,int stopj,int stopspecstate,int * startj,int * startspecstate) 
{
    AlnRange * out;  
    int jpos;    
    int state;   


    if( mat->basematrix->type != BASEMATRIX_TYPE_SHADOW) {  
      warn("Bad error! - non shadow matrix type in AlnRangeSet_from_ProteinBlockAligner");   
      return NULL;   
      }  


    /* Assumme that we have specials (we should!). Read back along the specials till we have the finish point */ 
    if( read_special_strip_ProteinBlockAligner(mat,0,stopj,stopspecstate,&jpos,&state,NULL) == FALSE)    {  
      warn("In AlnRanger_build_ProteinBlockAligner alignment ending at %d, unable to read back specials. Will (evenutally) return a partial range set... BEWARE!",stopj);    
      return NULL;   
      }  
    if( state == START || jpos <= 0) 
      return NULL;   


    out = AlnRange_alloc();  


    out->starti = ProteinBlockAligner_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,0);  
    out->startj = ProteinBlockAligner_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,1);  
    out->startstate = ProteinBlockAligner_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,2);  
    out->stopi = ProteinBlockAligner_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,3);   
    out->stopj = ProteinBlockAligner_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,4);   
    out->stopstate = ProteinBlockAligner_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,5);   
    out->startscore = ProteinBlockAligner_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,6);  
    out->stopscore = ProteinBlockAligner_DC_SHADOW_SPECIAL(mat,0,jpos,state);    


    /* Now, we have to figure out where this state came from in the specials */ 
    max_matrix_to_special_ProteinBlockAligner(mat,out->starti,out->startj,out->startstate,out->startscore,&jpos,startj,startspecstate,&state,NULL);  
    if( jpos == ProteinBlockAligner_READ_OFF_ERROR)  {  
      warn("In AlnRange_build_ProteinBlockAligner alignment ending at %d, with aln range between %d-%d in j, unable to find source special, returning this range, but this could get tricky!",stopj,out->startj,out->stopj); 
      return out;    
      }  


    /* Put in the correct score for startstate, from the special */ 
    out->startscore = ProteinBlockAligner_DC_SHADOW_SPECIAL(mat,0,*startj,*startspecstate);  
    /* The correct j coords have been put into startj, startspecstate... so just return out */ 
    return out;  
}    


/* Function:  read_hidden_ProteinBlockAligner(mat,starti,startj,startstate,stopi,stopj,stopstate,out)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [ProteinBlockAligner *]
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
boolean read_hidden_ProteinBlockAligner(ProteinBlockAligner * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out) 
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


      max_hidden_ProteinBlockAligner(mat,startj,i,j,state,isspecial,&i,&j,&state,&isspecial,&cellscore); 


      if( i == ProteinBlockAligner_READ_OFF_ERROR)   {  
        warn("In ProteinBlockAligner hidden read off, between %d:%d,%d:%d - at got bad read off. Problem!",starti,startj,stopi,stopj);   
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
        warn("In ProteinBlockAligner hidden read off, between %d:%d,%d:%d - hit start cell, but not in start state. Can't be good!.",starti,startj,stopi,stopj); 
        return FALSE;    
        }  
      }  
    warn("In ProteinBlockAligner hidden read off, between %d:%d,%d:%d - gone past start cell (now in %d,%d,%d), can't be good news!.",starti,startj,stopi,stopj,i,j,state);  
    return FALSE;    
}    


/* Function:  max_hidden_ProteinBlockAligner(mat,hiddenj,i,j,state,isspecial,reti,retj,retstate,retspecial,cellscore)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [ProteinBlockAligner *]
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
int max_hidden_ProteinBlockAligner(ProteinBlockAligner * mat,int hiddenj,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore) 
{
    register int temp;   
    register int cscore; 


    *reti = (*retj) = (*retstate) = ProteinBlockAligner_READ_OFF_ERROR;  


    if( i < 0 || j < 0 || i > mat->q->seq->len || j > mat->t->seq->len)  {  
      warn("In ProteinBlockAligner matrix special read off - out of bounds on matrix [i,j is %d,%d state %d in standard matrix]",i,j,state); 
      return -1; 
      }  


    /* Then you have to select the correct switch statement to figure out the readoff      */ 
    /* Somewhat odd - reverse the order of calculation and return as soon as it is correct */ 
    cscore = ProteinBlockAligner_HIDDEN_MATRIX(mat,i,j,state);   
    switch(state)    { /*Switch state */ 
      case BLOCK_1 :     
        temp = cscore - (mat->bentry) -  (CompMat_AAMATCH(mat->m,CSEQ_PROTEIN_AMINOACID(mat->q,i),CSEQ_PROTEIN_AMINOACID(mat->t,j)));    
        if( temp == ProteinBlockAligner_HIDDEN_MATRIX(mat,i - 1,j - 1,UNALIGNED) )   {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = UNALIGNED; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - ProteinBlockAligner_HIDDEN_MATRIX(mat,i-1,j-1,UNALIGNED);  
            }  
          return ProteinBlockAligner_HIDDEN_MATRIX(mat,i - 1,j - 1,UNALIGNED);   
          }  
        /* Not allowing special sources.. skipping START */ 
        temp = cscore - (mat->b_self_trans) -  (CompMat_AAMATCH(mat->m,CSEQ_PROTEIN_AMINOACID(mat->q,i),CSEQ_PROTEIN_AMINOACID(mat->t,j)));  
        if( temp == ProteinBlockAligner_HIDDEN_MATRIX(mat,i - 1,j - 1,BLOCK_1) ) {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = BLOCK_1;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - ProteinBlockAligner_HIDDEN_MATRIX(mat,i-1,j-1,BLOCK_1);    
            }  
          return ProteinBlockAligner_HIDDEN_MATRIX(mat,i - 1,j - 1,BLOCK_1);     
          }  
        warn("Major problem (!) - in ProteinBlockAligner read off, position %d,%d state %d no source found!",i,j,state); 
        return (-1); 
      case BLOCK_2 :     
        temp = cscore - (mat->bfor_trans) -  (CompMat_AAMATCH(mat->m,CSEQ_PROTEIN_AMINOACID(mat->q,i),CSEQ_PROTEIN_AMINOACID(mat->t,j)));    
        if( temp == ProteinBlockAligner_HIDDEN_MATRIX(mat,i - 1,j - 1,BLOCK_1) ) {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = BLOCK_1;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - ProteinBlockAligner_HIDDEN_MATRIX(mat,i-1,j-1,BLOCK_1);    
            }  
          return ProteinBlockAligner_HIDDEN_MATRIX(mat,i - 1,j - 1,BLOCK_1);     
          }  
        temp = cscore - (mat->b_self_trans) -  (CompMat_AAMATCH(mat->m,CSEQ_PROTEIN_AMINOACID(mat->q,i),CSEQ_PROTEIN_AMINOACID(mat->t,j)));  
        if( temp == ProteinBlockAligner_HIDDEN_MATRIX(mat,i - 1,j - 1,BLOCK_2) ) {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = BLOCK_2;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - ProteinBlockAligner_HIDDEN_MATRIX(mat,i-1,j-1,BLOCK_2);    
            }  
          return ProteinBlockAligner_HIDDEN_MATRIX(mat,i - 1,j - 1,BLOCK_2);     
          }  
        warn("Major problem (!) - in ProteinBlockAligner read off, position %d,%d state %d no source found!",i,j,state); 
        return (-1); 
      case BLOCK_3 :     
        temp = cscore - (mat->bfor_trans) -  (CompMat_AAMATCH(mat->m,CSEQ_PROTEIN_AMINOACID(mat->q,i),CSEQ_PROTEIN_AMINOACID(mat->t,j)));    
        if( temp == ProteinBlockAligner_HIDDEN_MATRIX(mat,i - 1,j - 1,BLOCK_2) ) {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = BLOCK_2;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - ProteinBlockAligner_HIDDEN_MATRIX(mat,i-1,j-1,BLOCK_2);    
            }  
          return ProteinBlockAligner_HIDDEN_MATRIX(mat,i - 1,j - 1,BLOCK_2);     
          }  
        temp = cscore - (mat->b_self_trans) -  (CompMat_AAMATCH(mat->m,CSEQ_PROTEIN_AMINOACID(mat->q,i),CSEQ_PROTEIN_AMINOACID(mat->t,j)));  
        if( temp == ProteinBlockAligner_HIDDEN_MATRIX(mat,i - 1,j - 1,BLOCK_3) ) {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = BLOCK_3;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - ProteinBlockAligner_HIDDEN_MATRIX(mat,i-1,j-1,BLOCK_3);    
            }  
          return ProteinBlockAligner_HIDDEN_MATRIX(mat,i - 1,j - 1,BLOCK_3);     
          }  
        warn("Major problem (!) - in ProteinBlockAligner read off, position %d,%d state %d no source found!",i,j,state); 
        return (-1); 
      case UNALIGNED :   
        temp = cscore - (0) -  (0);  
        if( temp == ProteinBlockAligner_HIDDEN_MATRIX(mat,i - 1,j - 0,UNALIGNED) )   {  
          *reti = i - 1; 
          *retj = j - 0; 
          *retstate = UNALIGNED; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - ProteinBlockAligner_HIDDEN_MATRIX(mat,i-1,j-0,UNALIGNED);  
            }  
          return ProteinBlockAligner_HIDDEN_MATRIX(mat,i - 1,j - 0,UNALIGNED);   
          }  
        temp = cscore - (0) -  (0);  
        if( temp == ProteinBlockAligner_HIDDEN_MATRIX(mat,i - 0,j - 1,UNALIGNED) )   {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = UNALIGNED; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - ProteinBlockAligner_HIDDEN_MATRIX(mat,i-0,j-1,UNALIGNED);  
            }  
          return ProteinBlockAligner_HIDDEN_MATRIX(mat,i - 0,j - 1,UNALIGNED);   
          }  
        temp = cscore - (mat->b3exit) -  (0);    
        if( temp == ProteinBlockAligner_HIDDEN_MATRIX(mat,i - 1,j - 0,BLOCK_2) ) {  
          *reti = i - 1; 
          *retj = j - 0; 
          *retstate = BLOCK_2;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - ProteinBlockAligner_HIDDEN_MATRIX(mat,i-1,j-0,BLOCK_2);    
            }  
          return ProteinBlockAligner_HIDDEN_MATRIX(mat,i - 1,j - 0,BLOCK_2);     
          }  
        temp = cscore - (mat->b3exit) -  (0);    
        if( temp == ProteinBlockAligner_HIDDEN_MATRIX(mat,i - 0,j - 1,BLOCK_2) ) {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = BLOCK_2;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - ProteinBlockAligner_HIDDEN_MATRIX(mat,i-0,j-1,BLOCK_2);    
            }  
          return ProteinBlockAligner_HIDDEN_MATRIX(mat,i - 0,j - 1,BLOCK_2);     
          }  
        temp = cscore - (mat->bexit) -  (0); 
        if( temp == ProteinBlockAligner_HIDDEN_MATRIX(mat,i - 1,j - 0,BLOCK_1) ) {  
          *reti = i - 1; 
          *retj = j - 0; 
          *retstate = BLOCK_1;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - ProteinBlockAligner_HIDDEN_MATRIX(mat,i-1,j-0,BLOCK_1);    
            }  
          return ProteinBlockAligner_HIDDEN_MATRIX(mat,i - 1,j - 0,BLOCK_1);     
          }  
        temp = cscore - (mat->bexit) -  (0); 
        if( temp == ProteinBlockAligner_HIDDEN_MATRIX(mat,i - 0,j - 1,BLOCK_1) ) {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = BLOCK_1;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - ProteinBlockAligner_HIDDEN_MATRIX(mat,i-0,j-1,BLOCK_1);    
            }  
          return ProteinBlockAligner_HIDDEN_MATRIX(mat,i - 0,j - 1,BLOCK_1);     
          }  
        warn("Major problem (!) - in ProteinBlockAligner read off, position %d,%d state %d no source found!",i,j,state); 
        return (-1); 
      default:   
        warn("Major problem (!) - in ProteinBlockAligner read off, position %d,%d state %d no source found!",i,j,state); 
        return (-1); 
      } /* end of Switch state  */ 
}    


/* Function:  read_special_strip_ProteinBlockAligner(mat,stopi,stopj,stopstate,startj,startstate,out)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [ProteinBlockAligner *]
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
boolean read_special_strip_ProteinBlockAligner(ProteinBlockAligner * mat,int stopi,int stopj,int stopstate,int * startj,int * startstate,PackAln * out) 
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
    while( j > ProteinBlockAligner_DC_SHADOW_SPECIAL_SP(mat,i,j,state,4) && state != START)  { /*while more specials to eat up*/ 
      /* Put away current state, if we should */ 
      if(out != NULL)    {  
        pau = PackAlnUnit_alloc();  /* Should deal with memory overflow */ 
        pau->i = i;  
        pau->j = j;  
        pau->state =  state + 4; 
        add_PackAln(out,pau);    
        }  


      max_special_strip_ProteinBlockAligner(mat,i,j,state,isspecial,&i,&j,&state,&isspecial,&cellscore); 
      if( i == ProteinBlockAligner_READ_OFF_ERROR)   {  
        warn("In special strip read ProteinBlockAligner, got a bad read off error. Sorry!"); 
        return FALSE;    
        }  
      } /* end of while more specials to eat up */ 


    /* check to see we have not gone too far! */ 
    if( state != START && j < ProteinBlockAligner_DC_SHADOW_SPECIAL_SP(mat,i,j,state,4)) {  
      warn("In special strip read ProteinBlockAligner, at special [%d] state [%d] overshot!",j,state);   
      return FALSE;  
      }  
    /* Put away last state */ 
    if(out != NULL)  {  
      pau = PackAlnUnit_alloc();/* Should deal with memory overflow */ 
      pau->i = i;    
      pau->j = j;    
      pau->state =  state + 4;   
      add_PackAln(out,pau);  
      }  


    /* Put away where we are in startj and startstate */ 
    *startj = j; 
    *startstate = state; 
    return TRUE; 
}    


/* Function:  max_special_strip_ProteinBlockAligner(mat,i,j,state,isspecial,reti,retj,retstate,retspecial,cellscore)
 *
 * Descrip:    A pretty intense internal function. Deals with read-off only in specials
 *
 *
 * Arg:               mat [UNKN ] Undocumented argument [ProteinBlockAligner *]
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
int max_special_strip_ProteinBlockAligner(ProteinBlockAligner * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore) 
{
    int temp;    
    int cscore;  


    *reti = (*retj) = (*retstate) = ProteinBlockAligner_READ_OFF_ERROR;  
    if( isspecial == FALSE ) {  
      warn("In special strip max function for ProteinBlockAligner, got a non special start point. Problem! (bad!)"); 
      return (-1);   
      }  


    if( j < 0 || j > mat->t->seq->len)   {  
      warn("In ProteinBlockAligner matrix special read off - out of bounds on matrix [j is %d in special]",j);   
      return -1; 
      }  


    cscore = ProteinBlockAligner_DC_SHADOW_SPECIAL(mat,i,j,state);   
    switch(state)    { /*switch on special states*/ 
      case START :   
      case END :     
        /* Source BLOCK_1 is not a special */ 
        /* Source UNALIGNED is not a special */ 
      default:   
        warn("Major problem (!) - in ProteinBlockAligner special strip read off, position %d,%d state %d no source found  dropped into default on source switch!",i,j,state);    
        return (-1); 
      } /* end of switch on special states */ 
}    


/* Function:  max_matrix_to_special_ProteinBlockAligner(mat,i,j,state,cscore,reti,retj,retstate,retspecial,cellscore)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [ProteinBlockAligner *]
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
int max_matrix_to_special_ProteinBlockAligner(ProteinBlockAligner * mat,int i,int j,int state,int cscore,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore) 
{
    int temp;    
    *reti = (*retj) = (*retstate) = ProteinBlockAligner_READ_OFF_ERROR;  


    if( j < 0 || j > mat->lenj)  {  
      warn("In ProteinBlockAligner matrix to special read off - out of bounds on matrix [j is %d in special]",j);    
      return -1; 
      }  


    switch(state)    { /*Switch state */ 
      case BLOCK_1 :     
        /* Source UNALIGNED is not a special, should not get here! */ 
        temp = cscore - (0) -  (CompMat_AAMATCH(mat->m,CSEQ_PROTEIN_AMINOACID(mat->q,i),CSEQ_PROTEIN_AMINOACID(mat->t,j)));  
        if( temp == ProteinBlockAligner_DC_SHADOW_SPECIAL(mat,i - 1,j - 1,START) )   {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = START; 
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - ProteinBlockAligner_DC_SHADOW_SPECIAL(mat,i-1,j-1,START);  
            }  
          return ProteinBlockAligner_DC_SHADOW_MATRIX(mat,i - 1,j - 1,START) ;   
          }  
        /* Source BLOCK_1 is not a special, should not get here! */ 
        warn("Major problem (!) - in ProteinBlockAligner matrix to special read off, position %d,%d state %d no source found!",i,j,state);   
        return (-1); 
      case BLOCK_2 :     
        /* Source BLOCK_1 is not a special, should not get here! */ 
        /* Source BLOCK_2 is not a special, should not get here! */ 
        warn("Major problem (!) - in ProteinBlockAligner matrix to special read off, position %d,%d state %d no source found!",i,j,state);   
        return (-1); 
      case BLOCK_3 :     
        /* Source BLOCK_2 is not a special, should not get here! */ 
        /* Source BLOCK_3 is not a special, should not get here! */ 
        warn("Major problem (!) - in ProteinBlockAligner matrix to special read off, position %d,%d state %d no source found!",i,j,state);   
        return (-1); 
      case UNALIGNED :   
        /* Source UNALIGNED is not a special, should not get here! */ 
        /* Source UNALIGNED is not a special, should not get here! */ 
        /* Source BLOCK_2 is not a special, should not get here! */ 
        /* Source BLOCK_2 is not a special, should not get here! */ 
        /* Source BLOCK_1 is not a special, should not get here! */ 
        /* Source BLOCK_1 is not a special, should not get here! */ 
        warn("Major problem (!) - in ProteinBlockAligner matrix to special read off, position %d,%d state %d no source found!",i,j,state);   
        return (-1); 
      default:   
        warn("Major problem (!) - in ProteinBlockAligner read off, position %d,%d state %d no source found!",i,j,state); 
        return (-1); 
      } /* end of Switch state  */ 


}    


/* Function:  calculate_hidden_ProteinBlockAligner(mat,starti,startj,startstate,stopi,stopj,stopstate,dpenv)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [ProteinBlockAligner *]
 * Arg:            starti [UNKN ] Undocumented argument [int]
 * Arg:            startj [UNKN ] Undocumented argument [int]
 * Arg:        startstate [UNKN ] Undocumented argument [int]
 * Arg:             stopi [UNKN ] Undocumented argument [int]
 * Arg:             stopj [UNKN ] Undocumented argument [int]
 * Arg:         stopstate [UNKN ] Undocumented argument [int]
 * Arg:             dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 */
void calculate_hidden_ProteinBlockAligner(ProteinBlockAligner * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv) 
{
    register int i;  
    register int j;  
    register int score;  
    register int temp;   
    register int hiddenj;    


    hiddenj = startj;    


    init_hidden_ProteinBlockAligner(mat,starti,startj,stopi,stopj);  


    ProteinBlockAligner_HIDDEN_MATRIX(mat,starti,startj,startstate) = 0; 


    for(j=startj;j<=stopj;j++)   {  
      for(i=starti;i<=stopi;i++) {  
        /* Should *not* do very first cell as this is the one set to zero in one state! */ 
        if( i == starti && j == startj ) 
          continue;  
        if( dpenv != NULL && is_in_DPEnvelope(dpenv,i,j) == FALSE )  { /*Is not in envelope*/ 
          ProteinBlockAligner_HIDDEN_MATRIX(mat,i,j,BLOCK_1) = NEGI;     
          ProteinBlockAligner_HIDDEN_MATRIX(mat,i,j,BLOCK_2) = NEGI;     
          ProteinBlockAligner_HIDDEN_MATRIX(mat,i,j,BLOCK_3) = NEGI;     
          ProteinBlockAligner_HIDDEN_MATRIX(mat,i,j,UNALIGNED) = NEGI;   
          continue;  
          } /* end of Is not in envelope */ 


        /* For state BLOCK_1 */ 
        /* setting first movement to score */ 
        score = ProteinBlockAligner_HIDDEN_MATRIX(mat,i-1,j-1,BLOCK_1) + mat->b_self_trans;  
        /* From state UNALIGNED to state BLOCK_1 */ 
        temp = ProteinBlockAligner_HIDDEN_MATRIX(mat,i-1,j-1,UNALIGNED) + mat->bentry;   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for BLOCK_1 */ 
        /* Add any movement independant score and put away */ 
         score += CompMat_AAMATCH(mat->m,CSEQ_PROTEIN_AMINOACID(mat->q,i),CSEQ_PROTEIN_AMINOACID(mat->t,j));     
         ProteinBlockAligner_HIDDEN_MATRIX(mat,i,j,BLOCK_1) = score; 
        /* Finished calculating state BLOCK_1 */ 


        /* For state BLOCK_2 */ 
        /* setting first movement to score */ 
        score = ProteinBlockAligner_HIDDEN_MATRIX(mat,i-1,j-1,BLOCK_2) + mat->b_self_trans;  
        /* From state BLOCK_1 to state BLOCK_2 */ 
        temp = ProteinBlockAligner_HIDDEN_MATRIX(mat,i-1,j-1,BLOCK_1) + mat->bfor_trans;     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for BLOCK_2 */ 
        /* Add any movement independant score and put away */ 
         score += CompMat_AAMATCH(mat->m,CSEQ_PROTEIN_AMINOACID(mat->q,i),CSEQ_PROTEIN_AMINOACID(mat->t,j));     
         ProteinBlockAligner_HIDDEN_MATRIX(mat,i,j,BLOCK_2) = score; 
        /* Finished calculating state BLOCK_2 */ 


        /* For state BLOCK_3 */ 
        /* setting first movement to score */ 
        score = ProteinBlockAligner_HIDDEN_MATRIX(mat,i-1,j-1,BLOCK_3) + mat->b_self_trans;  
        /* From state BLOCK_2 to state BLOCK_3 */ 
        temp = ProteinBlockAligner_HIDDEN_MATRIX(mat,i-1,j-1,BLOCK_2) + mat->bfor_trans;     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for BLOCK_3 */ 
        /* Add any movement independant score and put away */ 
         score += CompMat_AAMATCH(mat->m,CSEQ_PROTEIN_AMINOACID(mat->q,i),CSEQ_PROTEIN_AMINOACID(mat->t,j));     
         ProteinBlockAligner_HIDDEN_MATRIX(mat,i,j,BLOCK_3) = score; 
        /* Finished calculating state BLOCK_3 */ 


        /* For state UNALIGNED */ 
        /* setting first movement to score */ 
        score = ProteinBlockAligner_HIDDEN_MATRIX(mat,i-0,j-1,BLOCK_1) + mat->bexit;     
        /* From state BLOCK_1 to state UNALIGNED */ 
        temp = ProteinBlockAligner_HIDDEN_MATRIX(mat,i-1,j-0,BLOCK_1) + mat->bexit;  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state BLOCK_2 to state UNALIGNED */ 
        temp = ProteinBlockAligner_HIDDEN_MATRIX(mat,i-0,j-1,BLOCK_2) + mat->b3exit;     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state BLOCK_2 to state UNALIGNED */ 
        temp = ProteinBlockAligner_HIDDEN_MATRIX(mat,i-1,j-0,BLOCK_2) + mat->b3exit;     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state UNALIGNED to state UNALIGNED */ 
        temp = ProteinBlockAligner_HIDDEN_MATRIX(mat,i-0,j-1,UNALIGNED) + 0;     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state UNALIGNED to state UNALIGNED */ 
        temp = ProteinBlockAligner_HIDDEN_MATRIX(mat,i-1,j-0,UNALIGNED) + 0;     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for UNALIGNED */ 
        /* Add any movement independant score and put away */ 
         ProteinBlockAligner_HIDDEN_MATRIX(mat,i,j,UNALIGNED) = score;   
        /* Finished calculating state UNALIGNED */ 
        }  
      }  


    return;  
}    


/* Function:  init_hidden_ProteinBlockAligner(mat,starti,startj,stopi,stopj)
 *
 * Descrip: No Description
 *
 * Arg:           mat [UNKN ] Undocumented argument [ProteinBlockAligner *]
 * Arg:        starti [UNKN ] Undocumented argument [int]
 * Arg:        startj [UNKN ] Undocumented argument [int]
 * Arg:         stopi [UNKN ] Undocumented argument [int]
 * Arg:         stopj [UNKN ] Undocumented argument [int]
 *
 */
void init_hidden_ProteinBlockAligner(ProteinBlockAligner * mat,int starti,int startj,int stopi,int stopj) 
{
    register int i;  
    register int j;  
    register int hiddenj;    


    hiddenj = startj;    
    for(j=(startj-1);j<=stopj;j++)   {  
      for(i=(starti-1);i<=stopi;i++) {  
        ProteinBlockAligner_HIDDEN_MATRIX(mat,i,j,BLOCK_1) = NEGI;
  
        ProteinBlockAligner_HIDDEN_MATRIX(mat,i,j,BLOCK_2) = NEGI;
  
        ProteinBlockAligner_HIDDEN_MATRIX(mat,i,j,BLOCK_3) = NEGI;
  
        ProteinBlockAligner_HIDDEN_MATRIX(mat,i,j,UNALIGNED) = NEGI;
    
        }  
      }  


    return;  
}    


/* Function:  full_dc_ProteinBlockAligner(mat,starti,startj,startstate,stopi,stopj,stopstate,out,donej,totalj,dpenv)
 *
 * Descrip:    The main divide-and-conquor routine. Basically, call /PackAln_calculate_small_ProteinBlockAligner
 *             Not this function, which is pretty hard core. 
 *             Function is given start/end points (in main matrix) for alignment
 *             It does some checks, decides whether start/end in j is small enough for explicit calc
 *               - if yes, calculates it, reads off into PackAln (out), adds the j distance to donej and returns TRUE
 *               - if no,  uses /do_dc_single_pass_ProteinBlockAligner to get mid-point
 *                          saves midpoint, and calls itself to do right portion then left portion
 *             right then left ensures PackAln is added the 'right' way, ie, back-to-front
 *             returns FALSE on any error, with a warning
 *
 *
 * Arg:               mat [UNKN ] Matrix with small memory implementation [ProteinBlockAligner *]
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
boolean full_dc_ProteinBlockAligner(ProteinBlockAligner * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out,int * donej,int totalj,DPEnvelope * dpenv) 
{
    int lstarti; 
    int lstartj; 
    int lstate;  


    if( mat->basematrix->type != BASEMATRIX_TYPE_SHADOW) {  
      warn("*Very* bad error! - non shadow matrix type in full_dc_ProteinBlockAligner"); 
      return FALSE;  
      }  


    if( starti == -1 || startj == -1 || startstate == -1 || stopi == -1 || stopstate == -1)  {  
      warn("In full dc program, passed bad indices, indices passed were %d:%d[%d] to %d:%d[%d]\n",starti,startj,startstate,stopi,stopj,stopstate);   
      return FALSE;  
      }  


    if( stopj - startj < 5)  {  
      log_full_error(REPORT,0,"[%d,%d][%d,%d] Explicit read off",starti,startj,stopi,stopj);/* Build hidden explicit matrix */ 
      calculate_hidden_ProteinBlockAligner(mat,starti,startj,startstate,stopi,stopj,stopstate,dpenv);    
      *donej += (stopj - startj);   /* Now read it off into out */ 
      if( read_hidden_ProteinBlockAligner(mat,starti,startj,startstate,stopi,stopj,stopstate,out) == FALSE)  {  
        warn("In full dc, at %d:%d,%d:%d got a bad hidden explicit read off... ",starti,startj,stopi,stopj); 
        return FALSE;    
        }  
      return TRUE;   
      }  


/* In actual divide and conquor */ 
    if( do_dc_single_pass_ProteinBlockAligner(mat,starti,startj,startstate,stopi,stopj,stopstate,dpenv,(int)(*donej*100)/totalj) == FALSE)   {  
      warn("In divide and conquor for ProteinBlockAligner, at bound %d:%d to %d:%d, unable to calculate midpoint. Problem!",starti,startj,stopi,stopj);  
      return FALSE;  
      }  


/* Ok... now we have to call on each side of the matrix */ 
/* We have to retrieve left hand side positions, as they will be vapped by the time we call LHS */ 
    lstarti= ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,0);   
    lstartj= ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,1);   
    lstate = ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,2);   


/* Call on right hand side: this lets us do the correct read off */ 
    if( full_dc_ProteinBlockAligner(mat,ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,3),ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,4),ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,5),stopi,stopj,stopstate,out,donej,totalj,dpenv) == FALSE)   {  
/* Warning already issued, simply chained back up to top */ 
      return FALSE;  
      }  
/* Call on left hand side */ 
    if( full_dc_ProteinBlockAligner(mat,starti,startj,startstate,lstarti,lstartj,lstate,out,donej,totalj,dpenv) == FALSE)    {  
/* Warning already issued, simply chained back up to top */ 
      return FALSE;  
      }  


    return TRUE;     
}    


/* Function:  do_dc_single_pass_ProteinBlockAligner(mat,starti,startj,startstate,stopi,stopj,stopstate,dpenv,perc_done)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [ProteinBlockAligner *]
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
boolean do_dc_single_pass_ProteinBlockAligner(ProteinBlockAligner * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv,int perc_done) 
{
    int halfj;   
    halfj = startj + ((stopj - startj)/2);   


    init_dc_ProteinBlockAligner(mat);    


    ProteinBlockAligner_DC_SHADOW_MATRIX(mat,starti,startj,startstate) = 0;  
    run_up_dc_ProteinBlockAligner(mat,starti,stopi,startj,halfj-1,dpenv,perc_done);  
    push_dc_at_merge_ProteinBlockAligner(mat,starti,stopi,halfj,&halfj,dpenv);   
    follow_on_dc_ProteinBlockAligner(mat,starti,stopi,halfj,stopj,dpenv,perc_done);  
    return TRUE; 
}    


/* Function:  push_dc_at_merge_ProteinBlockAligner(mat,starti,stopi,startj,stopj,dpenv)
 *
 * Descrip: No Description
 *
 * Arg:           mat [UNKN ] Undocumented argument [ProteinBlockAligner *]
 * Arg:        starti [UNKN ] Undocumented argument [int]
 * Arg:         stopi [UNKN ] Undocumented argument [int]
 * Arg:        startj [UNKN ] Undocumented argument [int]
 * Arg:         stopj [UNKN ] Undocumented argument [int *]
 * Arg:         dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 */
void push_dc_at_merge_ProteinBlockAligner(ProteinBlockAligner * mat,int starti,int stopi,int startj,int * stopj,DPEnvelope * dpenv) 
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
          ProteinBlockAligner_DC_SHADOW_MATRIX(mat,i,j,BLOCK_1) = NEGI;  
          ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,BLOCK_1,0) = (-100);   
          ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,BLOCK_1,1) = (-100);   
          ProteinBlockAligner_DC_SHADOW_MATRIX(mat,i,j,BLOCK_2) = NEGI;  
          ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,BLOCK_2,0) = (-100);   
          ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,BLOCK_2,1) = (-100);   
          ProteinBlockAligner_DC_SHADOW_MATRIX(mat,i,j,BLOCK_3) = NEGI;  
          ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,BLOCK_3,0) = (-100);   
          ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,BLOCK_3,1) = (-100);   
          ProteinBlockAligner_DC_SHADOW_MATRIX(mat,i,j,UNALIGNED) = NEGI;    
          ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,UNALIGNED,0) = (-100); 
          ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,UNALIGNED,1) = (-100); 
          continue;  
          } /* end of Is not in envelope */ 


        /* For state BLOCK_1, pushing when j - offj <= mergej */ 
        score = ProteinBlockAligner_DC_SHADOW_MATRIX(mat,i-1,j-1,BLOCK_1) + mat->b_self_trans;   
        if( j - 1 <= mergej) {  
          ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,BLOCK_1,0) = i-1;  
          ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,BLOCK_1,1) = j-1;  
          ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,BLOCK_1,2) = BLOCK_1;  
          ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,BLOCK_1,3) = i;    
          ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,BLOCK_1,4) = j;    
          ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,BLOCK_1,5) = BLOCK_1;  
          }  
        else {  
          for(k=0;k<7;k++)   
            ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,BLOCK_1,k) = ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,BLOCK_1,k); 
          }  


        temp = ProteinBlockAligner_DC_SHADOW_MATRIX(mat,i-1,j-1,UNALIGNED) + mat->bentry;    
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,BLOCK_1,0) = i-1;    
            ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,BLOCK_1,1) = j-1;    
            ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,BLOCK_1,2) = UNALIGNED;  
            ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,BLOCK_1,3) = i;  
            ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,BLOCK_1,4) = j;  
            ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,BLOCK_1,5) = BLOCK_1;    
            }  
          else   {  
            for(k=0;k<7;k++) 
              ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,BLOCK_1,k) = ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,UNALIGNED,k); 
            }  
          }  
        /* Add any movement independant score */ 
        score += CompMat_AAMATCH(mat->m,CSEQ_PROTEIN_AMINOACID(mat->q,i),CSEQ_PROTEIN_AMINOACID(mat->t,j));  
        ProteinBlockAligner_DC_SHADOW_MATRIX(mat,i,j,BLOCK_1) = score;   
        /* Finished with state BLOCK_1 */ 


        /* For state BLOCK_2, pushing when j - offj <= mergej */ 
        score = ProteinBlockAligner_DC_SHADOW_MATRIX(mat,i-1,j-1,BLOCK_2) + mat->b_self_trans;   
        if( j - 1 <= mergej) {  
          ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,BLOCK_2,0) = i-1;  
          ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,BLOCK_2,1) = j-1;  
          ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,BLOCK_2,2) = BLOCK_2;  
          ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,BLOCK_2,3) = i;    
          ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,BLOCK_2,4) = j;    
          ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,BLOCK_2,5) = BLOCK_2;  
          }  
        else {  
          for(k=0;k<7;k++)   
            ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,BLOCK_2,k) = ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,BLOCK_2,k); 
          }  


        temp = ProteinBlockAligner_DC_SHADOW_MATRIX(mat,i-1,j-1,BLOCK_1) + mat->bfor_trans;  
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,BLOCK_2,0) = i-1;    
            ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,BLOCK_2,1) = j-1;    
            ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,BLOCK_2,2) = BLOCK_1;    
            ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,BLOCK_2,3) = i;  
            ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,BLOCK_2,4) = j;  
            ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,BLOCK_2,5) = BLOCK_2;    
            }  
          else   {  
            for(k=0;k<7;k++) 
              ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,BLOCK_2,k) = ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,BLOCK_1,k);   
            }  
          }  
        /* Add any movement independant score */ 
        score += CompMat_AAMATCH(mat->m,CSEQ_PROTEIN_AMINOACID(mat->q,i),CSEQ_PROTEIN_AMINOACID(mat->t,j));  
        ProteinBlockAligner_DC_SHADOW_MATRIX(mat,i,j,BLOCK_2) = score;   
        /* Finished with state BLOCK_2 */ 


        /* For state BLOCK_3, pushing when j - offj <= mergej */ 
        score = ProteinBlockAligner_DC_SHADOW_MATRIX(mat,i-1,j-1,BLOCK_3) + mat->b_self_trans;   
        if( j - 1 <= mergej) {  
          ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,BLOCK_3,0) = i-1;  
          ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,BLOCK_3,1) = j-1;  
          ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,BLOCK_3,2) = BLOCK_3;  
          ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,BLOCK_3,3) = i;    
          ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,BLOCK_3,4) = j;    
          ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,BLOCK_3,5) = BLOCK_3;  
          }  
        else {  
          for(k=0;k<7;k++)   
            ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,BLOCK_3,k) = ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,BLOCK_3,k); 
          }  


        temp = ProteinBlockAligner_DC_SHADOW_MATRIX(mat,i-1,j-1,BLOCK_2) + mat->bfor_trans;  
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,BLOCK_3,0) = i-1;    
            ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,BLOCK_3,1) = j-1;    
            ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,BLOCK_3,2) = BLOCK_2;    
            ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,BLOCK_3,3) = i;  
            ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,BLOCK_3,4) = j;  
            ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,BLOCK_3,5) = BLOCK_3;    
            }  
          else   {  
            for(k=0;k<7;k++) 
              ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,BLOCK_3,k) = ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,BLOCK_2,k);   
            }  
          }  
        /* Add any movement independant score */ 
        score += CompMat_AAMATCH(mat->m,CSEQ_PROTEIN_AMINOACID(mat->q,i),CSEQ_PROTEIN_AMINOACID(mat->t,j));  
        ProteinBlockAligner_DC_SHADOW_MATRIX(mat,i,j,BLOCK_3) = score;   
        /* Finished with state BLOCK_3 */ 


        /* For state UNALIGNED, pushing when j - offj <= mergej */ 
        score = ProteinBlockAligner_DC_SHADOW_MATRIX(mat,i-0,j-1,BLOCK_1) + mat->bexit;  
        if( j - 1 <= mergej) {  
          ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,UNALIGNED,0) = i-0;    
          ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,UNALIGNED,1) = j-1;    
          ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,UNALIGNED,2) = BLOCK_1;    
          ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,UNALIGNED,3) = i;  
          ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,UNALIGNED,4) = j;  
          ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,UNALIGNED,5) = UNALIGNED;  
          }  
        else {  
          for(k=0;k<7;k++)   
            ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,UNALIGNED,k) = ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,BLOCK_1,k);   
          }  


        temp = ProteinBlockAligner_DC_SHADOW_MATRIX(mat,i-1,j-0,BLOCK_1) + mat->bexit;   
        if( temp > score)    {  
          score = temp;  


          if( j - 0 <= mergej)   {  
            ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,UNALIGNED,0) = i-1;  
            ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,UNALIGNED,1) = j-0;  
            ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,UNALIGNED,2) = BLOCK_1;  
            ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,UNALIGNED,3) = i;    
            ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,UNALIGNED,4) = j;    
            ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,UNALIGNED,5) = UNALIGNED;    
            }  
          else   {  
            for(k=0;k<7;k++) 
              ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,UNALIGNED,k) = ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 0,BLOCK_1,k); 
            }  
          }  


        temp = ProteinBlockAligner_DC_SHADOW_MATRIX(mat,i-0,j-1,BLOCK_2) + mat->b3exit;  
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,UNALIGNED,0) = i-0;  
            ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,UNALIGNED,1) = j-1;  
            ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,UNALIGNED,2) = BLOCK_2;  
            ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,UNALIGNED,3) = i;    
            ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,UNALIGNED,4) = j;    
            ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,UNALIGNED,5) = UNALIGNED;    
            }  
          else   {  
            for(k=0;k<7;k++) 
              ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,UNALIGNED,k) = ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,BLOCK_2,k); 
            }  
          }  


        temp = ProteinBlockAligner_DC_SHADOW_MATRIX(mat,i-1,j-0,BLOCK_2) + mat->b3exit;  
        if( temp > score)    {  
          score = temp;  


          if( j - 0 <= mergej)   {  
            ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,UNALIGNED,0) = i-1;  
            ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,UNALIGNED,1) = j-0;  
            ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,UNALIGNED,2) = BLOCK_2;  
            ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,UNALIGNED,3) = i;    
            ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,UNALIGNED,4) = j;    
            ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,UNALIGNED,5) = UNALIGNED;    
            }  
          else   {  
            for(k=0;k<7;k++) 
              ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,UNALIGNED,k) = ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 0,BLOCK_2,k); 
            }  
          }  


        temp = ProteinBlockAligner_DC_SHADOW_MATRIX(mat,i-0,j-1,UNALIGNED) + 0;  
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,UNALIGNED,0) = i-0;  
            ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,UNALIGNED,1) = j-1;  
            ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,UNALIGNED,2) = UNALIGNED;    
            ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,UNALIGNED,3) = i;    
            ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,UNALIGNED,4) = j;    
            ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,UNALIGNED,5) = UNALIGNED;    
            }  
          else   {  
            for(k=0;k<7;k++) 
              ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,UNALIGNED,k) = ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,UNALIGNED,k);   
            }  
          }  


        temp = ProteinBlockAligner_DC_SHADOW_MATRIX(mat,i-1,j-0,UNALIGNED) + 0;  
        if( temp > score)    {  
          score = temp;  


          if( j - 0 <= mergej)   {  
            ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,UNALIGNED,0) = i-1;  
            ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,UNALIGNED,1) = j-0;  
            ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,UNALIGNED,2) = UNALIGNED;    
            ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,UNALIGNED,3) = i;    
            ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,UNALIGNED,4) = j;    
            ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,UNALIGNED,5) = UNALIGNED;    
            }  
          else   {  
            for(k=0;k<7;k++) 
              ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,UNALIGNED,k) = ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 0,UNALIGNED,k);   
            }  
          }  
        /* Add any movement independant score */ 
        ProteinBlockAligner_DC_SHADOW_MATRIX(mat,i,j,UNALIGNED) = score;     
        /* Finished with state UNALIGNED */ 
        }  
      }  
    /* Put back j into * stop j so that calling function gets it correct */ 
    if( stopj == NULL)   
      warn("Bad news... NULL stopj pointer in push dc function. This means that calling function does not know how many cells I have done!");    
    else 
      *stopj = j;    


    return;  
}    


/* Function:  follow_on_dc_ProteinBlockAligner(mat,starti,stopi,startj,stopj,dpenv,perc_done)
 *
 * Descrip: No Description
 *
 * Arg:              mat [UNKN ] Undocumented argument [ProteinBlockAligner *]
 * Arg:           starti [UNKN ] Undocumented argument [int]
 * Arg:            stopi [UNKN ] Undocumented argument [int]
 * Arg:           startj [UNKN ] Undocumented argument [int]
 * Arg:            stopj [UNKN ] Undocumented argument [int]
 * Arg:            dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:        perc_done [UNKN ] Undocumented argument [int]
 *
 */
void follow_on_dc_ProteinBlockAligner(ProteinBlockAligner * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done) 
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
          ProteinBlockAligner_DC_SHADOW_MATRIX(mat,i,j,BLOCK_1) = NEGI;  
          ProteinBlockAligner_DC_SHADOW_MATRIX(mat,i,j,BLOCK_2) = NEGI;  
          ProteinBlockAligner_DC_SHADOW_MATRIX(mat,i,j,BLOCK_3) = NEGI;  
          ProteinBlockAligner_DC_SHADOW_MATRIX(mat,i,j,UNALIGNED) = NEGI;    
          continue;  
          } /* end of Is not in envelope */ 
        if( num % 1000 == 0 )    
          log_full_error(REPORT,0,"[%d%%%% done]After  mid-j %5d Cells done %d%%%%",perc_done,startj,(num*100)/total);   


        /* For state BLOCK_1 */ 
        /* setting first movement to score */ 
        score = ProteinBlockAligner_DC_SHADOW_MATRIX(mat,i-1,j-1,BLOCK_1) + mat->b_self_trans;   
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,BLOCK_1,k);   
        /* From state UNALIGNED to state BLOCK_1 */ 
        temp = ProteinBlockAligner_DC_SHADOW_MATRIX(mat,i-1,j-1,UNALIGNED) + mat->bentry;    
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,UNALIGNED,k);   
          }  


        /* Ok - finished max calculation for BLOCK_1 */ 
        /* Add any movement independant score and put away */ 
         score += CompMat_AAMATCH(mat->m,CSEQ_PROTEIN_AMINOACID(mat->q,i),CSEQ_PROTEIN_AMINOACID(mat->t,j));     
         ProteinBlockAligner_DC_SHADOW_MATRIX(mat,i,j,BLOCK_1) = score;  
        for(k=0;k<7;k++) 
          ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,BLOCK_1,k) = localshadow[k];   
        /* Now figure out if any specials need this score */ 
        /* Finished calculating state BLOCK_1 */ 


        /* For state BLOCK_2 */ 
        /* setting first movement to score */ 
        score = ProteinBlockAligner_DC_SHADOW_MATRIX(mat,i-1,j-1,BLOCK_2) + mat->b_self_trans;   
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,BLOCK_2,k);   
        /* From state BLOCK_1 to state BLOCK_2 */ 
        temp = ProteinBlockAligner_DC_SHADOW_MATRIX(mat,i-1,j-1,BLOCK_1) + mat->bfor_trans;  
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,BLOCK_1,k); 
          }  


        /* Ok - finished max calculation for BLOCK_2 */ 
        /* Add any movement independant score and put away */ 
         score += CompMat_AAMATCH(mat->m,CSEQ_PROTEIN_AMINOACID(mat->q,i),CSEQ_PROTEIN_AMINOACID(mat->t,j));     
         ProteinBlockAligner_DC_SHADOW_MATRIX(mat,i,j,BLOCK_2) = score;  
        for(k=0;k<7;k++) 
          ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,BLOCK_2,k) = localshadow[k];   
        /* Now figure out if any specials need this score */ 
        /* Finished calculating state BLOCK_2 */ 


        /* For state BLOCK_3 */ 
        /* setting first movement to score */ 
        score = ProteinBlockAligner_DC_SHADOW_MATRIX(mat,i-1,j-1,BLOCK_3) + mat->b_self_trans;   
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,BLOCK_3,k);   
        /* From state BLOCK_2 to state BLOCK_3 */ 
        temp = ProteinBlockAligner_DC_SHADOW_MATRIX(mat,i-1,j-1,BLOCK_2) + mat->bfor_trans;  
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,BLOCK_2,k); 
          }  


        /* Ok - finished max calculation for BLOCK_3 */ 
        /* Add any movement independant score and put away */ 
         score += CompMat_AAMATCH(mat->m,CSEQ_PROTEIN_AMINOACID(mat->q,i),CSEQ_PROTEIN_AMINOACID(mat->t,j));     
         ProteinBlockAligner_DC_SHADOW_MATRIX(mat,i,j,BLOCK_3) = score;  
        for(k=0;k<7;k++) 
          ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,BLOCK_3,k) = localshadow[k];   
        /* Now figure out if any specials need this score */ 
        /* Finished calculating state BLOCK_3 */ 


        /* For state UNALIGNED */ 
        /* setting first movement to score */ 
        score = ProteinBlockAligner_DC_SHADOW_MATRIX(mat,i-0,j-1,BLOCK_1) + mat->bexit;  
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,BLOCK_1,k);   
        /* From state BLOCK_1 to state UNALIGNED */ 
        temp = ProteinBlockAligner_DC_SHADOW_MATRIX(mat,i-1,j-0,BLOCK_1) + mat->bexit;   
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 0,BLOCK_1,k); 
          }  
        /* From state BLOCK_2 to state UNALIGNED */ 
        temp = ProteinBlockAligner_DC_SHADOW_MATRIX(mat,i-0,j-1,BLOCK_2) + mat->b3exit;  
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,BLOCK_2,k); 
          }  
        /* From state BLOCK_2 to state UNALIGNED */ 
        temp = ProteinBlockAligner_DC_SHADOW_MATRIX(mat,i-1,j-0,BLOCK_2) + mat->b3exit;  
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 0,BLOCK_2,k); 
          }  
        /* From state UNALIGNED to state UNALIGNED */ 
        temp = ProteinBlockAligner_DC_SHADOW_MATRIX(mat,i-0,j-1,UNALIGNED) + 0;  
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,UNALIGNED,k);   
          }  
        /* From state UNALIGNED to state UNALIGNED */ 
        temp = ProteinBlockAligner_DC_SHADOW_MATRIX(mat,i-1,j-0,UNALIGNED) + 0;  
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 0,UNALIGNED,k);   
          }  


        /* Ok - finished max calculation for UNALIGNED */ 
        /* Add any movement independant score and put away */ 
         ProteinBlockAligner_DC_SHADOW_MATRIX(mat,i,j,UNALIGNED) = score;    
        for(k=0;k<7;k++) 
          ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,UNALIGNED,k) = localshadow[k]; 
        /* Now figure out if any specials need this score */ 
        /* Finished calculating state UNALIGNED */ 
        } /* end of this is strip */ 
      } /* end of for each valid j column */ 


/* Function:  run_up_dc_ProteinBlockAligner(mat,starti,stopi,startj,stopj,dpenv,perc_done)
 *
 * Descrip: No Description
 *
 * Arg:              mat [UNKN ] Undocumented argument [ProteinBlockAligner *]
 * Arg:           starti [UNKN ] Undocumented argument [int]
 * Arg:            stopi [UNKN ] Undocumented argument [int]
 * Arg:           startj [UNKN ] Undocumented argument [int]
 * Arg:            stopj [UNKN ] Undocumented argument [int]
 * Arg:            dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:        perc_done [UNKN ] Undocumented argument [int]
 *
 */
}    
void run_up_dc_ProteinBlockAligner(ProteinBlockAligner * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done) 
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
          ProteinBlockAligner_DC_SHADOW_MATRIX(mat,i,j,BLOCK_1) = NEGI;  
          ProteinBlockAligner_DC_SHADOW_MATRIX(mat,i,j,BLOCK_2) = NEGI;  
          ProteinBlockAligner_DC_SHADOW_MATRIX(mat,i,j,BLOCK_3) = NEGI;  
          ProteinBlockAligner_DC_SHADOW_MATRIX(mat,i,j,UNALIGNED) = NEGI;    
          continue;  
          } /* end of Is not in envelope */ 
        if( num % 1000 == 0 )    
          log_full_error(REPORT,0,"[%d%%%% done]Before mid-j %5d Cells done %d%%%%",perc_done,stopj,(num*100)/total);    


        /* For state BLOCK_1 */ 
        /* setting first movement to score */ 
        score = ProteinBlockAligner_DC_SHADOW_MATRIX(mat,i-1,j-1,BLOCK_1) + mat->b_self_trans;   
        /* From state UNALIGNED to state BLOCK_1 */ 
        temp = ProteinBlockAligner_DC_SHADOW_MATRIX(mat,i-1,j-1,UNALIGNED) + mat->bentry;    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for BLOCK_1 */ 
        /* Add any movement independant score and put away */ 
         score += CompMat_AAMATCH(mat->m,CSEQ_PROTEIN_AMINOACID(mat->q,i),CSEQ_PROTEIN_AMINOACID(mat->t,j));     
         ProteinBlockAligner_DC_SHADOW_MATRIX(mat,i,j,BLOCK_1) = score;  
        /* Finished calculating state BLOCK_1 */ 


        /* For state BLOCK_2 */ 
        /* setting first movement to score */ 
        score = ProteinBlockAligner_DC_SHADOW_MATRIX(mat,i-1,j-1,BLOCK_2) + mat->b_self_trans;   
        /* From state BLOCK_1 to state BLOCK_2 */ 
        temp = ProteinBlockAligner_DC_SHADOW_MATRIX(mat,i-1,j-1,BLOCK_1) + mat->bfor_trans;  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for BLOCK_2 */ 
        /* Add any movement independant score and put away */ 
         score += CompMat_AAMATCH(mat->m,CSEQ_PROTEIN_AMINOACID(mat->q,i),CSEQ_PROTEIN_AMINOACID(mat->t,j));     
         ProteinBlockAligner_DC_SHADOW_MATRIX(mat,i,j,BLOCK_2) = score;  
        /* Finished calculating state BLOCK_2 */ 


        /* For state BLOCK_3 */ 
        /* setting first movement to score */ 
        score = ProteinBlockAligner_DC_SHADOW_MATRIX(mat,i-1,j-1,BLOCK_3) + mat->b_self_trans;   
        /* From state BLOCK_2 to state BLOCK_3 */ 
        temp = ProteinBlockAligner_DC_SHADOW_MATRIX(mat,i-1,j-1,BLOCK_2) + mat->bfor_trans;  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for BLOCK_3 */ 
        /* Add any movement independant score and put away */ 
         score += CompMat_AAMATCH(mat->m,CSEQ_PROTEIN_AMINOACID(mat->q,i),CSEQ_PROTEIN_AMINOACID(mat->t,j));     
         ProteinBlockAligner_DC_SHADOW_MATRIX(mat,i,j,BLOCK_3) = score;  
        /* Finished calculating state BLOCK_3 */ 


        /* For state UNALIGNED */ 
        /* setting first movement to score */ 
        score = ProteinBlockAligner_DC_SHADOW_MATRIX(mat,i-0,j-1,BLOCK_1) + mat->bexit;  
        /* From state BLOCK_1 to state UNALIGNED */ 
        temp = ProteinBlockAligner_DC_SHADOW_MATRIX(mat,i-1,j-0,BLOCK_1) + mat->bexit;   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state BLOCK_2 to state UNALIGNED */ 
        temp = ProteinBlockAligner_DC_SHADOW_MATRIX(mat,i-0,j-1,BLOCK_2) + mat->b3exit;  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state BLOCK_2 to state UNALIGNED */ 
        temp = ProteinBlockAligner_DC_SHADOW_MATRIX(mat,i-1,j-0,BLOCK_2) + mat->b3exit;  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state UNALIGNED to state UNALIGNED */ 
        temp = ProteinBlockAligner_DC_SHADOW_MATRIX(mat,i-0,j-1,UNALIGNED) + 0;  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state UNALIGNED to state UNALIGNED */ 
        temp = ProteinBlockAligner_DC_SHADOW_MATRIX(mat,i-1,j-0,UNALIGNED) + 0;  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for UNALIGNED */ 
        /* Add any movement independant score and put away */ 
         ProteinBlockAligner_DC_SHADOW_MATRIX(mat,i,j,UNALIGNED) = score;    
        /* Finished calculating state UNALIGNED */ 
        } /* end of this is strip */ 
      } /* end of for each valid j column */ 


/* Function:  init_dc_ProteinBlockAligner(mat)
 *
 * Descrip: No Description
 *
 * Arg:        mat [UNKN ] Undocumented argument [ProteinBlockAligner *]
 *
 */
}    
void init_dc_ProteinBlockAligner(ProteinBlockAligner * mat) 
{
    register int i;  
    register int j;  
    register int k;  


    for(j=0;j<3;j++) {  
      for(i=(-1);i<mat->q->seq->len;i++) {  
        ProteinBlockAligner_DC_SHADOW_MATRIX(mat,i,j,BLOCK_1) = NEGI;    
        ProteinBlockAligner_DC_SHADOW_MATRIX(mat,i,j,BLOCK_2) = NEGI;    
        ProteinBlockAligner_DC_SHADOW_MATRIX(mat,i,j,BLOCK_3) = NEGI;    
        ProteinBlockAligner_DC_SHADOW_MATRIX(mat,i,j,UNALIGNED) = NEGI;  
        for(k=0;k<7;k++) {  
          ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,BLOCK_1,k) = (-1); 
          ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,BLOCK_2,k) = (-1); 
          ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,BLOCK_3,k) = (-1); 
          ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,UNALIGNED,k) = (-1);   
          }  
        }  
      }  


    return;  
}    


/* Function:  start_end_find_end_ProteinBlockAligner(mat,endj)
 *
 * Descrip:    First function used to find end of the best path in the special state !end
 *
 *
 * Arg:         mat [UNKN ] Matrix in small mode [ProteinBlockAligner *]
 * Arg:        endj [WRITE] position of end in j (meaningless in i) [int *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int start_end_find_end_ProteinBlockAligner(ProteinBlockAligner * mat,int * endj) 
{
    register int j;  
    register int max;    
    register int maxj;   


    max = ProteinBlockAligner_DC_SHADOW_SPECIAL(mat,0,mat->t->seq->len-1,END);   
    maxj = mat->t->seq->len-1;   
    for(j= mat->t->seq->len-2 ;j >= 0 ;j--)  {  
      if( ProteinBlockAligner_DC_SHADOW_SPECIAL(mat,0,j,END) > max ) {  
        max = ProteinBlockAligner_DC_SHADOW_SPECIAL(mat,0,j,END);    
        maxj = j;    
        }  
      }  


    if( endj != NULL)    
      *endj = maxj;  


    return max;  
}    


/* Function:  dc_optimised_start_end_calc_ProteinBlockAligner(*mat,dpenv)
 *
 * Descrip:    Calculates special strip, leaving start/end/score points in shadow matrix
 *             Works off specially laid out memory from steve searle
 *
 *
 * Arg:         *mat [UNKN ] Undocumented argument [ProteinBlockAligner]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean dc_optimised_start_end_calc_ProteinBlockAligner(ProteinBlockAligner *mat,DPEnvelope * dpenv) 
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
    leni = mat->q->seq->len; 
    lenj = mat->t->seq->len; 
    total = leni * lenj; 


    score_pointers = (int *) calloc (1 * (leni + 1) * 4,sizeof(int));    
    shadow_pointers = (int *) calloc (1 * (leni + 1) * 4 * 8,sizeof(int));   


    for(j=0;j<lenj;j++)  { /*for each j strip*/ 
      for(i=0;i<leni;i++)    { /*for each i position in strip*/ 
        num++;   
        if( dpenv != NULL && is_in_DPEnvelope(dpenv,i,j) == FALSE )  { /*Is not in envelope*/ 
          ProteinBlockAligner_DC_OPT_SHADOW_MATRIX(mat,i,j,BLOCK_1) = NEGI;  
          ProteinBlockAligner_DC_OPT_SHADOW_MATRIX(mat,i,j,BLOCK_2) = NEGI;  
          ProteinBlockAligner_DC_OPT_SHADOW_MATRIX(mat,i,j,BLOCK_3) = NEGI;  
          ProteinBlockAligner_DC_OPT_SHADOW_MATRIX(mat,i,j,UNALIGNED) = NEGI;    
          continue;  
          } /* end of Is not in envelope */ 
        if( num%1000 == 0)   
          log_full_error(REPORT,0,"%6d Cells done [%2d%%%%]",num,num*100/total); 




        /* For state BLOCK_1 */ 
        /* setting first movement to score */ 
        score = ProteinBlockAligner_DC_OPT_SHADOW_MATRIX(mat,i-1,j-1,BLOCK_1) + mat->b_self_trans + (CompMat_AAMATCH(mat->m,CSEQ_PROTEIN_AMINOACID(mat->q,i),CSEQ_PROTEIN_AMINOACID(mat->t,j)));     
        /* assign local shadown pointer */ 
        localsp = &(ProteinBlockAligner_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 1,BLOCK_1,0)); 
        /* From state START to state BLOCK_1 */ 
        temp = ProteinBlockAligner_DC_OPT_SHADOW_SPECIAL(mat,i-1,j-1,START) + 0 + (CompMat_AAMATCH(mat->m,CSEQ_PROTEIN_AMINOACID(mat->q,i),CSEQ_PROTEIN_AMINOACID(mat->t,j)));   
        if( temp  > score )  {  
          score = temp;  
          /* This state [START] is a special for BLOCK_1... push top shadow pointers here */ 
          localshadow[0]= i; 
          localshadow[1]= j; 
          localshadow[2]= BLOCK_1;   
          localshadow[3]= (-1);  
          localshadow[4]= (-1);  
          localshadow[5]= (-1);  
          localshadow[6]= score; 
          localsp = localshadow; 
          }  
        /* From state UNALIGNED to state BLOCK_1 */ 
        temp = ProteinBlockAligner_DC_OPT_SHADOW_MATRIX(mat,i-1,j-1,UNALIGNED) + mat->bentry +(CompMat_AAMATCH(mat->m,CSEQ_PROTEIN_AMINOACID(mat->q,i),CSEQ_PROTEIN_AMINOACID(mat->t,j)));   
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(ProteinBlockAligner_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 1,UNALIGNED,0)); 
          }  


        /* Ok - finished max calculation for BLOCK_1 */ 
        /* Add any movement independant score and put away */ 
        /* Actually, already done inside scores */ 
         ProteinBlockAligner_DC_OPT_SHADOW_MATRIX(mat,i,j,BLOCK_1) = score;  
        for(k=0;k<7;k++) 
          ProteinBlockAligner_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,BLOCK_1,k) = localsp[k];   
        /* Now figure out if any specials need this score */ 


        /* state BLOCK_1 is a source for special END */ 
        temp = score + (0) + (0) ;   
        if( temp > ProteinBlockAligner_DC_OPT_SHADOW_SPECIAL(mat,i,j,END) )  {  
          ProteinBlockAligner_DC_OPT_SHADOW_SPECIAL(mat,i,j,END) = temp;     
          /* Have to push only bottem half of system here */ 
          for(k=0;k<3;k++)   
            ProteinBlockAligner_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,k) = ProteinBlockAligner_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,BLOCK_1,k);    
          ProteinBlockAligner_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,6) = ProteinBlockAligner_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,BLOCK_1,6);  
          ProteinBlockAligner_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,3) = i;   
          ProteinBlockAligner_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,4) = j;   
          ProteinBlockAligner_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,5) = BLOCK_1; 
          }  




        /* Finished calculating state BLOCK_1 */ 


        /* For state BLOCK_2 */ 
        /* setting first movement to score */ 
        score = ProteinBlockAligner_DC_OPT_SHADOW_MATRIX(mat,i-1,j-1,BLOCK_2) + mat->b_self_trans + (CompMat_AAMATCH(mat->m,CSEQ_PROTEIN_AMINOACID(mat->q,i),CSEQ_PROTEIN_AMINOACID(mat->t,j)));     
        /* assign local shadown pointer */ 
        localsp = &(ProteinBlockAligner_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 1,BLOCK_2,0)); 
        /* From state BLOCK_1 to state BLOCK_2 */ 
        temp = ProteinBlockAligner_DC_OPT_SHADOW_MATRIX(mat,i-1,j-1,BLOCK_1) + mat->bfor_trans +(CompMat_AAMATCH(mat->m,CSEQ_PROTEIN_AMINOACID(mat->q,i),CSEQ_PROTEIN_AMINOACID(mat->t,j)));     
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(ProteinBlockAligner_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 1,BLOCK_1,0));   
          }  


        /* Ok - finished max calculation for BLOCK_2 */ 
        /* Add any movement independant score and put away */ 
        /* Actually, already done inside scores */ 
         ProteinBlockAligner_DC_OPT_SHADOW_MATRIX(mat,i,j,BLOCK_2) = score;  
        for(k=0;k<7;k++) 
          ProteinBlockAligner_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,BLOCK_2,k) = localsp[k];   
        /* Now figure out if any specials need this score */ 


        /* Finished calculating state BLOCK_2 */ 


        /* For state BLOCK_3 */ 
        /* setting first movement to score */ 
        score = ProteinBlockAligner_DC_OPT_SHADOW_MATRIX(mat,i-1,j-1,BLOCK_3) + mat->b_self_trans + (CompMat_AAMATCH(mat->m,CSEQ_PROTEIN_AMINOACID(mat->q,i),CSEQ_PROTEIN_AMINOACID(mat->t,j)));     
        /* assign local shadown pointer */ 
        localsp = &(ProteinBlockAligner_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 1,BLOCK_3,0)); 
        /* From state BLOCK_2 to state BLOCK_3 */ 
        temp = ProteinBlockAligner_DC_OPT_SHADOW_MATRIX(mat,i-1,j-1,BLOCK_2) + mat->bfor_trans +(CompMat_AAMATCH(mat->m,CSEQ_PROTEIN_AMINOACID(mat->q,i),CSEQ_PROTEIN_AMINOACID(mat->t,j)));     
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(ProteinBlockAligner_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 1,BLOCK_2,0));   
          }  


        /* Ok - finished max calculation for BLOCK_3 */ 
        /* Add any movement independant score and put away */ 
        /* Actually, already done inside scores */ 
         ProteinBlockAligner_DC_OPT_SHADOW_MATRIX(mat,i,j,BLOCK_3) = score;  
        for(k=0;k<7;k++) 
          ProteinBlockAligner_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,BLOCK_3,k) = localsp[k];   
        /* Now figure out if any specials need this score */ 


        /* Finished calculating state BLOCK_3 */ 


        /* For state UNALIGNED */ 
        /* setting first movement to score */ 
        score = ProteinBlockAligner_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,BLOCK_1) + mat->bexit + (0);    
        /* assign local shadown pointer */ 
        localsp = &(ProteinBlockAligner_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,BLOCK_1,0)); 
        /* From state BLOCK_1 to state UNALIGNED */ 
        temp = ProteinBlockAligner_DC_OPT_SHADOW_MATRIX(mat,i-1,j-0,BLOCK_1) + mat->bexit +(0);  
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(ProteinBlockAligner_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 0,BLOCK_1,0));   
          }  
        /* From state BLOCK_2 to state UNALIGNED */ 
        temp = ProteinBlockAligner_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,BLOCK_2) + mat->b3exit +(0);     
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(ProteinBlockAligner_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,BLOCK_2,0));   
          }  
        /* From state BLOCK_2 to state UNALIGNED */ 
        temp = ProteinBlockAligner_DC_OPT_SHADOW_MATRIX(mat,i-1,j-0,BLOCK_2) + mat->b3exit +(0);     
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(ProteinBlockAligner_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 0,BLOCK_2,0));   
          }  
        /* From state UNALIGNED to state UNALIGNED */ 
        temp = ProteinBlockAligner_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,UNALIGNED) + 0 +(0);     
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(ProteinBlockAligner_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,UNALIGNED,0)); 
          }  
        /* From state UNALIGNED to state UNALIGNED */ 
        temp = ProteinBlockAligner_DC_OPT_SHADOW_MATRIX(mat,i-1,j-0,UNALIGNED) + 0 +(0);     
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(ProteinBlockAligner_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 0,UNALIGNED,0)); 
          }  


        /* Ok - finished max calculation for UNALIGNED */ 
        /* Add any movement independant score and put away */ 
        /* Actually, already done inside scores */ 
         ProteinBlockAligner_DC_OPT_SHADOW_MATRIX(mat,i,j,UNALIGNED) = score;    
        for(k=0;k<7;k++) 
          ProteinBlockAligner_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,UNALIGNED,k) = localsp[k]; 
        /* Now figure out if any specials need this score */ 


        /* state UNALIGNED is a source for special END */ 
        temp = score + (0) + (0) ;   
        if( temp > ProteinBlockAligner_DC_OPT_SHADOW_SPECIAL(mat,i,j,END) )  {  
          ProteinBlockAligner_DC_OPT_SHADOW_SPECIAL(mat,i,j,END) = temp;     
          /* Have to push only bottem half of system here */ 
          for(k=0;k<3;k++)   
            ProteinBlockAligner_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,k) = ProteinBlockAligner_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,UNALIGNED,k);  
          ProteinBlockAligner_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,6) = ProteinBlockAligner_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,UNALIGNED,6);    
          ProteinBlockAligner_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,3) = i;   
          ProteinBlockAligner_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,4) = j;   
          ProteinBlockAligner_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,5) = UNALIGNED;   
          }  




        /* Finished calculating state UNALIGNED */ 


        } /* end of for each i position in strip */ 
      } /* end of for each j strip */ 
    free(score_pointers);    
    free(shadow_pointers);   
    return TRUE;     
}    


/* Function:  init_start_end_linear_ProteinBlockAligner(mat)
 *
 * Descrip: No Description
 *
 * Arg:        mat [UNKN ] Undocumented argument [ProteinBlockAligner *]
 *
 */
void init_start_end_linear_ProteinBlockAligner(ProteinBlockAligner * mat) 
{
    register int i;  
    register int j;  
    for(j=0;j<3;j++) {  
      for(i=(-1);i<mat->q->seq->len;i++) {  
        ProteinBlockAligner_DC_SHADOW_MATRIX(mat,i,j,BLOCK_1) = NEGI;    
        ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,BLOCK_1,0) = (-1);   
        ProteinBlockAligner_DC_SHADOW_MATRIX(mat,i,j,BLOCK_2) = NEGI;    
        ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,BLOCK_2,0) = (-1);   
        ProteinBlockAligner_DC_SHADOW_MATRIX(mat,i,j,BLOCK_3) = NEGI;    
        ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,BLOCK_3,0) = (-1);   
        ProteinBlockAligner_DC_SHADOW_MATRIX(mat,i,j,UNALIGNED) = NEGI;  
        ProteinBlockAligner_DC_SHADOW_MATRIX_SP(mat,i,j,UNALIGNED,0) = (-1); 
        }  
      }  


    for(j=(-1);j<mat->t->seq->len;j++)   {  
      ProteinBlockAligner_DC_SHADOW_SPECIAL(mat,0,j,START) = 0;  
      ProteinBlockAligner_DC_SHADOW_SPECIAL_SP(mat,0,j,START,0) = j; 
      ProteinBlockAligner_DC_SHADOW_SPECIAL(mat,0,j,END) = NEGI; 
      ProteinBlockAligner_DC_SHADOW_SPECIAL_SP(mat,0,j,END,0) = (-1);    
      }  


    return;  
}    


/* Function:  convert_PackAln_to_AlnBlock_ProteinBlockAligner(pal)
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
AlnBlock * convert_PackAln_to_AlnBlock_ProteinBlockAligner(PackAln * pal) 
{
    AlnConvertSet * acs; 
    AlnBlock * alb;  


    acs = AlnConvertSet_ProteinBlockAligner();   
    alb = AlnBlock_from_PackAln(acs,pal);    
    free_AlnConvertSet(acs); 
    return alb;  
}    


 static char * query_label[] = { "SEQUENCE_BLOCK_1","SEQUENCE_BLOCK_2","UNALIGNED","END" };  
/* Function:  AlnConvertSet_ProteinBlockAligner(void)
 *
 * Descrip: No Description
 *
 *
 * Return [UNKN ]  Undocumented return value [AlnConvertSet *]
 *
 */
 static char * target_label[] = { "SEQUENCE_BLOCK_1","SEQUENCE_BLOCK_2","UNALIGNED","END" }; 
AlnConvertSet * AlnConvertSet_ProteinBlockAligner(void) 
{
    AlnConvertUnit * acu;    
    AlnConvertSet  * out;    


    out = AlnConvertSet_alloc_std(); 


    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = BLOCK_1;   
    acu->state2 = BLOCK_1;   
    acu->offi = 1;   
    acu->offj = 1;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[0];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = START + 4; 
    acu->is_from_special = TRUE; 
    acu->state2 = BLOCK_1;   
    acu->offi = (-1);    
    acu->offj = 1;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[0];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = UNALIGNED; 
    acu->state2 = BLOCK_1;   
    acu->offi = 1;   
    acu->offj = 1;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[0];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = BLOCK_2;   
    acu->state2 = BLOCK_2;   
    acu->offi = 1;   
    acu->offj = 1;   
    acu->label1 = query_label[1];    
    acu->label2 = target_label[1];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = BLOCK_1;   
    acu->state2 = BLOCK_2;   
    acu->offi = 1;   
    acu->offj = 1;   
    acu->label1 = query_label[1];    
    acu->label2 = target_label[1];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = BLOCK_3;   
    acu->state2 = BLOCK_3;   
    acu->offi = 1;   
    acu->offj = 1;   
    acu->label1 = query_label[1];    
    acu->label2 = target_label[1];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = BLOCK_2;   
    acu->state2 = BLOCK_3;   
    acu->offi = 1;   
    acu->offj = 1;   
    acu->label1 = query_label[1];    
    acu->label2 = target_label[1];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = BLOCK_1;   
    acu->state2 = UNALIGNED;     
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[2];    
    acu->label2 = target_label[2];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = BLOCK_1;   
    acu->state2 = UNALIGNED;     
    acu->offi = 1;   
    acu->offj = 0;   
    acu->label1 = query_label[2];    
    acu->label2 = target_label[2];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = BLOCK_2;   
    acu->state2 = UNALIGNED;     
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[2];    
    acu->label2 = target_label[2];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = BLOCK_2;   
    acu->state2 = UNALIGNED;     
    acu->offi = 1;   
    acu->offj = 0;   
    acu->label1 = query_label[2];    
    acu->label2 = target_label[2];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = UNALIGNED; 
    acu->state2 = UNALIGNED;     
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[2];    
    acu->label2 = target_label[2];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = UNALIGNED; 
    acu->state2 = UNALIGNED;     
    acu->offi = 1;   
    acu->offj = 0;   
    acu->label1 = query_label[2];    
    acu->label2 = target_label[2];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = UNALIGNED; 
    acu->state2 = END + 4;   
    acu->offi = (-1);    
    acu->offj = 0;   
    acu->label1 = query_label[3];    
    acu->label2 = target_label[3];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = BLOCK_1;   
    acu->state2 = END + 4;   
    acu->offi = (-1);    
    acu->offj = 0;   
    acu->label1 = query_label[3];    
    acu->label2 = target_label[3];   
    return out;  
}    


/* Function:  PackAln_read_Expl_ProteinBlockAligner(mat)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [ProteinBlockAligner *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * PackAln_read_Expl_ProteinBlockAligner(ProteinBlockAligner * mat) 
{
    ProteinBlockAligner_access_func_holder holder;   


    holder.access_main    = ProteinBlockAligner_explicit_access_main;    
    holder.access_special = ProteinBlockAligner_explicit_access_special; 
    return PackAln_read_generic_ProteinBlockAligner(mat,holder); 
}    


/* Function:  ProteinBlockAligner_explicit_access_main(mat,i,j,state)
 *
 * Descrip: No Description
 *
 * Arg:          mat [UNKN ] Undocumented argument [ProteinBlockAligner *]
 * Arg:            i [UNKN ] Undocumented argument [int]
 * Arg:            j [UNKN ] Undocumented argument [int]
 * Arg:        state [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int ProteinBlockAligner_explicit_access_main(ProteinBlockAligner * mat,int i,int j,int state) 
{
    return ProteinBlockAligner_EXPL_MATRIX(mat,i,j,state);   
}    


/* Function:  ProteinBlockAligner_explicit_access_special(mat,i,j,state)
 *
 * Descrip: No Description
 *
 * Arg:          mat [UNKN ] Undocumented argument [ProteinBlockAligner *]
 * Arg:            i [UNKN ] Undocumented argument [int]
 * Arg:            j [UNKN ] Undocumented argument [int]
 * Arg:        state [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int ProteinBlockAligner_explicit_access_special(ProteinBlockAligner * mat,int i,int j,int state) 
{
    return ProteinBlockAligner_EXPL_SPECIAL(mat,i,j,state);  
}    


/* Function:  PackAln_read_generic_ProteinBlockAligner(mat,h)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [ProteinBlockAligner *]
 * Arg:          h [UNKN ] Undocumented argument [ProteinBlockAligner_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * PackAln_read_generic_ProteinBlockAligner(ProteinBlockAligner * mat,ProteinBlockAligner_access_func_holder h) 
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


    out->score =  find_end_ProteinBlockAligner(mat,&i,&j,&state,&isspecial,h);   


    /* Add final end transition (at the moment we have not got the score! */ 
    if( (pau= PackAlnUnit_alloc()) == NULL  || add_PackAln(out,pau) == FALSE )   {  
      warn("Failed the first PackAlnUnit alloc, %d length of Alignment in ProteinBlockAligner_basic_read, returning a mess.(Sorry!)",out->len);  
      return out;    
      }  


    /* Put in positions for end trans. Remember that coordinates in C style */ 
    pau->i = i;  
    pau->j = j;  
    if( isspecial != TRUE)   
      pau->state = state;    
    else pau->state = state + 4;     
    prev=pau;    
    while( state != START || isspecial != TRUE)  { /*while state != START*/ 


      if( isspecial == TRUE )    
        max_calc_special_ProteinBlockAligner(mat,i,j,state,isspecial,&i,&j,&state,&isspecial,&cellscore,h);  
      else   
        max_calc_ProteinBlockAligner(mat,i,j,state,isspecial,&i,&j,&state,&isspecial,&cellscore,h);  
      if(i == ProteinBlockAligner_READ_OFF_ERROR || j == ProteinBlockAligner_READ_OFF_ERROR || state == ProteinBlockAligner_READ_OFF_ERROR ) {  
        warn("Problem - hit bad read off system, exiting now");  
        break;   
        }  
      if( (pau= PackAlnUnit_alloc()) == NULL  || add_PackAln(out,pau) == FALSE ) {  
        warn("Failed a PackAlnUnit alloc, %d length of Alignment in ProteinBlockAligner_basic_read, returning partial alignment",out->len);  
        break;   
        }  


      /* Put in positions for block. Remember that coordinates in C style */ 
      pau->i = i;    
      pau->j = j;    
      if( isspecial != TRUE)     
        pau->state = state;  
      else pau->state = state + 4;   
      prev->score = cellscore;   
      prev = pau;    
      } /* end of while state != START */ 


    invert_PackAln(out); 
    return out;  
}    


/* Function:  find_end_ProteinBlockAligner(mat,ri,rj,state,isspecial,h)
 *
 * Descrip: No Description
 *
 * Arg:              mat [UNKN ] Undocumented argument [ProteinBlockAligner *]
 * Arg:               ri [UNKN ] Undocumented argument [int *]
 * Arg:               rj [UNKN ] Undocumented argument [int *]
 * Arg:            state [UNKN ] Undocumented argument [int *]
 * Arg:        isspecial [UNKN ] Undocumented argument [boolean *]
 * Arg:                h [UNKN ] Undocumented argument [ProteinBlockAligner_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int find_end_ProteinBlockAligner(ProteinBlockAligner * mat,int * ri,int * rj,int * state,boolean * isspecial,ProteinBlockAligner_access_func_holder h) 
{
    int j;   
    int max; 
    int maxj;    
    int temp;    


    max = (*h.access_special)(mat,0,mat->t->seq->len-1,END); 
    maxj = mat->t->seq->len-1;   
    for(j= mat->t->seq->len-2 ;j >= 0 ;j--)  {  
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


/* Function:  ProteinBlockAligner_debug_show_matrix(mat,starti,stopi,startj,stopj,ofp)
 *
 * Descrip: No Description
 *
 * Arg:           mat [UNKN ] Undocumented argument [ProteinBlockAligner *]
 * Arg:        starti [UNKN ] Undocumented argument [int]
 * Arg:         stopi [UNKN ] Undocumented argument [int]
 * Arg:        startj [UNKN ] Undocumented argument [int]
 * Arg:         stopj [UNKN ] Undocumented argument [int]
 * Arg:           ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
void ProteinBlockAligner_debug_show_matrix(ProteinBlockAligner * mat,int starti,int stopi,int startj,int stopj,FILE * ofp) 
{
    register int i;  
    register int j;  


    for(i=starti;i<stopi && i < mat->q->seq->len;i++)    {  
      for(j=startj;j<stopj && j < mat->t->seq->len;j++)  {  
        fprintf(ofp,"Cell [%d - %d]\n",i,j);     
        fprintf(ofp,"State BLOCK_1 %d\n",ProteinBlockAligner_EXPL_MATRIX(mat,i,j,BLOCK_1));  
        fprintf(ofp,"State BLOCK_2 %d\n",ProteinBlockAligner_EXPL_MATRIX(mat,i,j,BLOCK_2));  
        fprintf(ofp,"State BLOCK_3 %d\n",ProteinBlockAligner_EXPL_MATRIX(mat,i,j,BLOCK_3));  
        fprintf(ofp,"State UNALIGNED %d\n",ProteinBlockAligner_EXPL_MATRIX(mat,i,j,UNALIGNED));  
        fprintf(ofp,"\n\n"); 
        }  
      }  


}    


/* Function:  max_calc_ProteinBlockAligner(mat,i,j,state,isspecial,reti,retj,retstate,retspecial,cellscore,h)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [ProteinBlockAligner *]
 * Arg:                 i [UNKN ] Undocumented argument [int]
 * Arg:                 j [UNKN ] Undocumented argument [int]
 * Arg:             state [UNKN ] Undocumented argument [int]
 * Arg:         isspecial [UNKN ] Undocumented argument [boolean]
 * Arg:              reti [UNKN ] Undocumented argument [int *]
 * Arg:              retj [UNKN ] Undocumented argument [int *]
 * Arg:          retstate [UNKN ] Undocumented argument [int *]
 * Arg:        retspecial [UNKN ] Undocumented argument [boolean *]
 * Arg:         cellscore [UNKN ] Undocumented argument [int *]
 * Arg:                 h [UNKN ] Undocumented argument [ProteinBlockAligner_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int max_calc_ProteinBlockAligner(ProteinBlockAligner * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,ProteinBlockAligner_access_func_holder h) 
{
    register int temp;   
    register int cscore; 


    *reti = (*retj) = (*retstate) = ProteinBlockAligner_READ_OFF_ERROR;  


    if( i < 0 || j < 0 || i > mat->q->seq->len || j > mat->t->seq->len)  {  
      warn("In ProteinBlockAligner matrix special read off - out of bounds on matrix [i,j is %d,%d state %d in standard matrix]",i,j,state); 
      return -1;     
      }  


    /* Then you have to select the correct switch statement to figure out the readoff      */ 
    /* Somewhat odd - reverse the order of calculation and return as soon as it is correct */ 
    cscore = (*h.access_main)(mat,i,j,state);    
    switch(state)    { /*Switch state */ 
      case BLOCK_1 :     
        temp = cscore - (mat->bentry) -  (CompMat_AAMATCH(mat->m,CSEQ_PROTEIN_AMINOACID(mat->q,i),CSEQ_PROTEIN_AMINOACID(mat->t,j)));    
        if( temp == (*h.access_main)(mat,i - 1,j - 1,UNALIGNED) )    {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = UNALIGNED; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-1,UNALIGNED);   
            }  
          return (*h.access_main)(mat,i - 1,j - 1,UNALIGNED);    
          }  
        temp = cscore - (0) -  (CompMat_AAMATCH(mat->m,CSEQ_PROTEIN_AMINOACID(mat->q,i),CSEQ_PROTEIN_AMINOACID(mat->t,j)));  
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
        temp = cscore - (mat->b_self_trans) -  (CompMat_AAMATCH(mat->m,CSEQ_PROTEIN_AMINOACID(mat->q,i),CSEQ_PROTEIN_AMINOACID(mat->t,j)));  
        if( temp == (*h.access_main)(mat,i - 1,j - 1,BLOCK_1) )  {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = BLOCK_1;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-1,BLOCK_1); 
            }  
          return (*h.access_main)(mat,i - 1,j - 1,BLOCK_1);  
          }  
        warn("Major problem (!) - in ProteinBlockAligner read off, position %d,%d state %d no source found!",i,j,state); 
        return (-1); 
      case BLOCK_2 :     
        temp = cscore - (mat->bfor_trans) -  (CompMat_AAMATCH(mat->m,CSEQ_PROTEIN_AMINOACID(mat->q,i),CSEQ_PROTEIN_AMINOACID(mat->t,j)));    
        if( temp == (*h.access_main)(mat,i - 1,j - 1,BLOCK_1) )  {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = BLOCK_1;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-1,BLOCK_1); 
            }  
          return (*h.access_main)(mat,i - 1,j - 1,BLOCK_1);  
          }  
        temp = cscore - (mat->b_self_trans) -  (CompMat_AAMATCH(mat->m,CSEQ_PROTEIN_AMINOACID(mat->q,i),CSEQ_PROTEIN_AMINOACID(mat->t,j)));  
        if( temp == (*h.access_main)(mat,i - 1,j - 1,BLOCK_2) )  {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = BLOCK_2;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-1,BLOCK_2); 
            }  
          return (*h.access_main)(mat,i - 1,j - 1,BLOCK_2);  
          }  
        warn("Major problem (!) - in ProteinBlockAligner read off, position %d,%d state %d no source found!",i,j,state); 
        return (-1); 
      case BLOCK_3 :     
        temp = cscore - (mat->bfor_trans) -  (CompMat_AAMATCH(mat->m,CSEQ_PROTEIN_AMINOACID(mat->q,i),CSEQ_PROTEIN_AMINOACID(mat->t,j)));    
        if( temp == (*h.access_main)(mat,i - 1,j - 1,BLOCK_2) )  {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = BLOCK_2;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-1,BLOCK_2); 
            }  
          return (*h.access_main)(mat,i - 1,j - 1,BLOCK_2);  
          }  
        temp = cscore - (mat->b_self_trans) -  (CompMat_AAMATCH(mat->m,CSEQ_PROTEIN_AMINOACID(mat->q,i),CSEQ_PROTEIN_AMINOACID(mat->t,j)));  
        if( temp == (*h.access_main)(mat,i - 1,j - 1,BLOCK_3) )  {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = BLOCK_3;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-1,BLOCK_3); 
            }  
          return (*h.access_main)(mat,i - 1,j - 1,BLOCK_3);  
          }  
        warn("Major problem (!) - in ProteinBlockAligner read off, position %d,%d state %d no source found!",i,j,state); 
        return (-1); 
      case UNALIGNED :   
        temp = cscore - (0) -  (0);  
        if( temp == (*h.access_main)(mat,i - 1,j - 0,UNALIGNED) )    {  
          *reti = i - 1; 
          *retj = j - 0; 
          *retstate = UNALIGNED; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-0,UNALIGNED);   
            }  
          return (*h.access_main)(mat,i - 1,j - 0,UNALIGNED);    
          }  
        temp = cscore - (0) -  (0);  
        if( temp == (*h.access_main)(mat,i - 0,j - 1,UNALIGNED) )    {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = UNALIGNED; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,UNALIGNED);   
            }  
          return (*h.access_main)(mat,i - 0,j - 1,UNALIGNED);    
          }  
        temp = cscore - (mat->b3exit) -  (0);    
        if( temp == (*h.access_main)(mat,i - 1,j - 0,BLOCK_2) )  {  
          *reti = i - 1; 
          *retj = j - 0; 
          *retstate = BLOCK_2;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-0,BLOCK_2); 
            }  
          return (*h.access_main)(mat,i - 1,j - 0,BLOCK_2);  
          }  
        temp = cscore - (mat->b3exit) -  (0);    
        if( temp == (*h.access_main)(mat,i - 0,j - 1,BLOCK_2) )  {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = BLOCK_2;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,BLOCK_2); 
            }  
          return (*h.access_main)(mat,i - 0,j - 1,BLOCK_2);  
          }  
        temp = cscore - (mat->bexit) -  (0); 
        if( temp == (*h.access_main)(mat,i - 1,j - 0,BLOCK_1) )  {  
          *reti = i - 1; 
          *retj = j - 0; 
          *retstate = BLOCK_1;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-0,BLOCK_1); 
            }  
          return (*h.access_main)(mat,i - 1,j - 0,BLOCK_1);  
          }  
        temp = cscore - (mat->bexit) -  (0); 
        if( temp == (*h.access_main)(mat,i - 0,j - 1,BLOCK_1) )  {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = BLOCK_1;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,BLOCK_1); 
            }  
          return (*h.access_main)(mat,i - 0,j - 1,BLOCK_1);  
          }  
        warn("Major problem (!) - in ProteinBlockAligner read off, position %d,%d state %d no source found!",i,j,state); 
        return (-1); 
      default:   
        warn("Major problem (!) - in ProteinBlockAligner read off, position %d,%d state %d no source found!",i,j,state); 
        return (-1); 
      } /* end of Switch state  */ 
}    


/* Function:  max_calc_special_ProteinBlockAligner(mat,i,j,state,isspecial,reti,retj,retstate,retspecial,cellscore,h)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [ProteinBlockAligner *]
 * Arg:                 i [UNKN ] Undocumented argument [int]
 * Arg:                 j [UNKN ] Undocumented argument [int]
 * Arg:             state [UNKN ] Undocumented argument [int]
 * Arg:         isspecial [UNKN ] Undocumented argument [boolean]
 * Arg:              reti [UNKN ] Undocumented argument [int *]
 * Arg:              retj [UNKN ] Undocumented argument [int *]
 * Arg:          retstate [UNKN ] Undocumented argument [int *]
 * Arg:        retspecial [UNKN ] Undocumented argument [boolean *]
 * Arg:         cellscore [UNKN ] Undocumented argument [int *]
 * Arg:                 h [UNKN ] Undocumented argument [ProteinBlockAligner_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int max_calc_special_ProteinBlockAligner(ProteinBlockAligner * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,ProteinBlockAligner_access_func_holder h) 
{
    register int temp;   
    register int cscore; 


    *reti = (*retj) = (*retstate) = ProteinBlockAligner_READ_OFF_ERROR;  


    if( j < 0 || j > mat->t->seq->len)   {  
      warn("In ProteinBlockAligner matrix special read off - out of bounds on matrix [j is %d in special]",j);   
      return -1;     
      }  


    cscore = (*h.access_special)(mat,i,j,state); 
    switch(state)    { /*switch on special states*/ 
      case START :   
      case END :     
        /* source BLOCK_1 is from main matrix */ 
        for(i= mat->q->seq->len-1;i >= 0 ;i--)   { /*for i >= 0*/ 
          temp = cscore - (0) - (0);     
          if( temp == (*h.access_main)(mat,i - 0,j - 0,BLOCK_1) )    {  
            *reti = i - 0;   
            *retj = j - 0;   
            *retstate = BLOCK_1; 
            *retspecial = FALSE; 
            if( cellscore != NULL)   {  
              *cellscore = cscore - (*h.access_main)(mat,i-0,j-0,BLOCK_1);   
              }  
            return (*h.access_main)(mat,i - 0,j - 0,BLOCK_1) ;   
            }  
          } /* end of for i >= 0 */ 
        /* source UNALIGNED is from main matrix */ 
        for(i= mat->q->seq->len-1;i >= 0 ;i--)   { /*for i >= 0*/ 
          temp = cscore - (0) - (0);     
          if( temp == (*h.access_main)(mat,i - 0,j - 0,UNALIGNED) )  {  
            *reti = i - 0;   
            *retj = j - 0;   
            *retstate = UNALIGNED;   
            *retspecial = FALSE; 
            if( cellscore != NULL)   {  
              *cellscore = cscore - (*h.access_main)(mat,i-0,j-0,UNALIGNED);     
              }  
            return (*h.access_main)(mat,i - 0,j - 0,UNALIGNED) ;     
            }  
          } /* end of for i >= 0 */ 
      default:   
        warn("Major problem (!) - in ProteinBlockAligner read off, position %d,%d state %d no source found  dropped into default on source switch!",i,j,state);  
        return (-1); 
      } /* end of switch on special states */ 
}    


/* Function:  calculate_ProteinBlockAligner(mat)
 *
 * Descrip:    This function calculates the ProteinBlockAligner matrix when in explicit mode
 *             To allocate the matrix use /allocate_Expl_ProteinBlockAligner
 *
 *
 * Arg:        mat [UNKN ] ProteinBlockAligner which contains explicit basematrix memory [ProteinBlockAligner *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean calculate_ProteinBlockAligner(ProteinBlockAligner * mat) 
{
    int i;   
    int j;   
    int leni;    
    int lenj;    
    int tot; 
    int num; 


    if( mat->basematrix->type != BASEMATRIX_TYPE_EXPLICIT )  {  
      warn("in calculate_ProteinBlockAligner, passed a non Explicit matrix type, cannot calculate!");    
      return FALSE;  
      }  


    leni = mat->leni;    
    lenj = mat->lenj;    
    tot = leni * lenj;   
    num = 0; 


    start_reporting("ProteinBlockAligner Matrix calculation: "); 
    for(j=0;j<lenj;j++)  {  
      auto int score;    
      auto int temp;     
      for(i=0;i<leni;i++)    {  
        if( num%1000 == 0 )  
          log_full_error(REPORT,0,"[%7d] Cells %2d%%%%",num,num*100/tot);    
        num++;   


        /* For state BLOCK_1 */ 
        /* setting first movement to score */ 
        score = ProteinBlockAligner_EXPL_MATRIX(mat,i-1,j-1,BLOCK_1) + mat->b_self_trans;    
        /* From state START to state BLOCK_1 */ 
        temp = ProteinBlockAligner_EXPL_SPECIAL(mat,i-1,j-1,START) + 0;  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state UNALIGNED to state BLOCK_1 */ 
        temp = ProteinBlockAligner_EXPL_MATRIX(mat,i-1,j-1,UNALIGNED) + mat->bentry;     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for BLOCK_1 */ 
        /* Add any movement independant score and put away */ 
         score += CompMat_AAMATCH(mat->m,CSEQ_PROTEIN_AMINOACID(mat->q,i),CSEQ_PROTEIN_AMINOACID(mat->t,j));     
         ProteinBlockAligner_EXPL_MATRIX(mat,i,j,BLOCK_1) = score;   


        /* state BLOCK_1 is a source for special END */ 
        temp = score + (0) + (0) ;   
        if( temp > ProteinBlockAligner_EXPL_SPECIAL(mat,i,j,END) )   {  
          ProteinBlockAligner_EXPL_SPECIAL(mat,i,j,END) = temp;  
          }  




        /* Finished calculating state BLOCK_1 */ 


        /* For state BLOCK_2 */ 
        /* setting first movement to score */ 
        score = ProteinBlockAligner_EXPL_MATRIX(mat,i-1,j-1,BLOCK_2) + mat->b_self_trans;    
        /* From state BLOCK_1 to state BLOCK_2 */ 
        temp = ProteinBlockAligner_EXPL_MATRIX(mat,i-1,j-1,BLOCK_1) + mat->bfor_trans;   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for BLOCK_2 */ 
        /* Add any movement independant score and put away */ 
         score += CompMat_AAMATCH(mat->m,CSEQ_PROTEIN_AMINOACID(mat->q,i),CSEQ_PROTEIN_AMINOACID(mat->t,j));     
         ProteinBlockAligner_EXPL_MATRIX(mat,i,j,BLOCK_2) = score;   


        /* Finished calculating state BLOCK_2 */ 


        /* For state BLOCK_3 */ 
        /* setting first movement to score */ 
        score = ProteinBlockAligner_EXPL_MATRIX(mat,i-1,j-1,BLOCK_3) + mat->b_self_trans;    
        /* From state BLOCK_2 to state BLOCK_3 */ 
        temp = ProteinBlockAligner_EXPL_MATRIX(mat,i-1,j-1,BLOCK_2) + mat->bfor_trans;   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for BLOCK_3 */ 
        /* Add any movement independant score and put away */ 
         score += CompMat_AAMATCH(mat->m,CSEQ_PROTEIN_AMINOACID(mat->q,i),CSEQ_PROTEIN_AMINOACID(mat->t,j));     
         ProteinBlockAligner_EXPL_MATRIX(mat,i,j,BLOCK_3) = score;   


        /* Finished calculating state BLOCK_3 */ 


        /* For state UNALIGNED */ 
        /* setting first movement to score */ 
        score = ProteinBlockAligner_EXPL_MATRIX(mat,i-0,j-1,BLOCK_1) + mat->bexit;   
        /* From state BLOCK_1 to state UNALIGNED */ 
        temp = ProteinBlockAligner_EXPL_MATRIX(mat,i-1,j-0,BLOCK_1) + mat->bexit;    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state BLOCK_2 to state UNALIGNED */ 
        temp = ProteinBlockAligner_EXPL_MATRIX(mat,i-0,j-1,BLOCK_2) + mat->b3exit;   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state BLOCK_2 to state UNALIGNED */ 
        temp = ProteinBlockAligner_EXPL_MATRIX(mat,i-1,j-0,BLOCK_2) + mat->b3exit;   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state UNALIGNED to state UNALIGNED */ 
        temp = ProteinBlockAligner_EXPL_MATRIX(mat,i-0,j-1,UNALIGNED) + 0;   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state UNALIGNED to state UNALIGNED */ 
        temp = ProteinBlockAligner_EXPL_MATRIX(mat,i-1,j-0,UNALIGNED) + 0;   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for UNALIGNED */ 
        /* Add any movement independant score and put away */ 
         ProteinBlockAligner_EXPL_MATRIX(mat,i,j,UNALIGNED) = score; 


        /* state UNALIGNED is a source for special END */ 
        temp = score + (0) + (0) ;   
        if( temp > ProteinBlockAligner_EXPL_SPECIAL(mat,i,j,END) )   {  
          ProteinBlockAligner_EXPL_SPECIAL(mat,i,j,END) = temp;  
          }  




        /* Finished calculating state UNALIGNED */ 
        }  


      /* Special state START has no special to special movements */ 


      /* Special state END has no special to special movements */ 
      }  
    stop_reporting();    
    return TRUE;     
}    


/* Function:  calculate_dpenv_ProteinBlockAligner(mat,dpenv)
 *
 * Descrip:    This function calculates the ProteinBlockAligner matrix when in explicit mode, subject to the envelope
 *
 *
 * Arg:          mat [UNKN ] ProteinBlockAligner which contains explicit basematrix memory [ProteinBlockAligner *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean calculate_dpenv_ProteinBlockAligner(ProteinBlockAligner * mat,DPEnvelope * dpenv) 
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
      warn("in calculate_ProteinBlockAligner, passed a non Explicit matrix type, cannot calculate!");    
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
        ProteinBlockAligner_EXPL_MATRIX(mat,i,j,BLOCK_1) = NEGI; 
        ProteinBlockAligner_EXPL_MATRIX(mat,i,j,BLOCK_2) = NEGI; 
        ProteinBlockAligner_EXPL_MATRIX(mat,i,j,BLOCK_3) = NEGI; 
        ProteinBlockAligner_EXPL_MATRIX(mat,i,j,UNALIGNED) = NEGI;   
        }  
      }  
    for(j=-1;j<mat->lenj;j++)    {  
      ProteinBlockAligner_EXPL_SPECIAL(mat,i,j,START) = 0;   
      ProteinBlockAligner_EXPL_SPECIAL(mat,i,j,END) = NEGI;  
      }  


    start_reporting("ProteinBlockAligner Matrix calculation: "); 
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
          ProteinBlockAligner_EXPL_MATRIX(mat,i,j,BLOCK_1) = NEGI;   
          ProteinBlockAligner_EXPL_MATRIX(mat,i,j,BLOCK_2) = NEGI;   
          ProteinBlockAligner_EXPL_MATRIX(mat,i,j,BLOCK_3) = NEGI;   
          ProteinBlockAligner_EXPL_MATRIX(mat,i,j,UNALIGNED) = NEGI; 
          continue;  
          }  


        if( num%1000 == 0 )  
          log_full_error(REPORT,0,"[%7d] Cells %2d%%%%",num,num*100/tot);    
        num++;   


        /* For state BLOCK_1 */ 
        /* setting first movement to score */ 
        score = ProteinBlockAligner_EXPL_MATRIX(mat,i-1,j-1,BLOCK_1) + mat->b_self_trans;    
        /* From state START to state BLOCK_1 */ 
        temp = ProteinBlockAligner_EXPL_SPECIAL(mat,i-1,j-1,START) + 0;  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state UNALIGNED to state BLOCK_1 */ 
        temp = ProteinBlockAligner_EXPL_MATRIX(mat,i-1,j-1,UNALIGNED) + mat->bentry;     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for BLOCK_1 */ 
        /* Add any movement independant score and put away */ 
         score += CompMat_AAMATCH(mat->m,CSEQ_PROTEIN_AMINOACID(mat->q,i),CSEQ_PROTEIN_AMINOACID(mat->t,j));     
         ProteinBlockAligner_EXPL_MATRIX(mat,i,j,BLOCK_1) = score;   


        /* state BLOCK_1 is a source for special END */ 
        temp = score + (0) + (0) ;   
        if( temp > ProteinBlockAligner_EXPL_SPECIAL(mat,i,j,END) )   {  
          ProteinBlockAligner_EXPL_SPECIAL(mat,i,j,END) = temp;  
          }  




        /* Finished calculating state BLOCK_1 */ 


        /* For state BLOCK_2 */ 
        /* setting first movement to score */ 
        score = ProteinBlockAligner_EXPL_MATRIX(mat,i-1,j-1,BLOCK_2) + mat->b_self_trans;    
        /* From state BLOCK_1 to state BLOCK_2 */ 
        temp = ProteinBlockAligner_EXPL_MATRIX(mat,i-1,j-1,BLOCK_1) + mat->bfor_trans;   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for BLOCK_2 */ 
        /* Add any movement independant score and put away */ 
         score += CompMat_AAMATCH(mat->m,CSEQ_PROTEIN_AMINOACID(mat->q,i),CSEQ_PROTEIN_AMINOACID(mat->t,j));     
         ProteinBlockAligner_EXPL_MATRIX(mat,i,j,BLOCK_2) = score;   


        /* Finished calculating state BLOCK_2 */ 


        /* For state BLOCK_3 */ 
        /* setting first movement to score */ 
        score = ProteinBlockAligner_EXPL_MATRIX(mat,i-1,j-1,BLOCK_3) + mat->b_self_trans;    
        /* From state BLOCK_2 to state BLOCK_3 */ 
        temp = ProteinBlockAligner_EXPL_MATRIX(mat,i-1,j-1,BLOCK_2) + mat->bfor_trans;   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for BLOCK_3 */ 
        /* Add any movement independant score and put away */ 
         score += CompMat_AAMATCH(mat->m,CSEQ_PROTEIN_AMINOACID(mat->q,i),CSEQ_PROTEIN_AMINOACID(mat->t,j));     
         ProteinBlockAligner_EXPL_MATRIX(mat,i,j,BLOCK_3) = score;   


        /* Finished calculating state BLOCK_3 */ 


        /* For state UNALIGNED */ 
        /* setting first movement to score */ 
        score = ProteinBlockAligner_EXPL_MATRIX(mat,i-0,j-1,BLOCK_1) + mat->bexit;   
        /* From state BLOCK_1 to state UNALIGNED */ 
        temp = ProteinBlockAligner_EXPL_MATRIX(mat,i-1,j-0,BLOCK_1) + mat->bexit;    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state BLOCK_2 to state UNALIGNED */ 
        temp = ProteinBlockAligner_EXPL_MATRIX(mat,i-0,j-1,BLOCK_2) + mat->b3exit;   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state BLOCK_2 to state UNALIGNED */ 
        temp = ProteinBlockAligner_EXPL_MATRIX(mat,i-1,j-0,BLOCK_2) + mat->b3exit;   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state UNALIGNED to state UNALIGNED */ 
        temp = ProteinBlockAligner_EXPL_MATRIX(mat,i-0,j-1,UNALIGNED) + 0;   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state UNALIGNED to state UNALIGNED */ 
        temp = ProteinBlockAligner_EXPL_MATRIX(mat,i-1,j-0,UNALIGNED) + 0;   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for UNALIGNED */ 
        /* Add any movement independant score and put away */ 
         ProteinBlockAligner_EXPL_MATRIX(mat,i,j,UNALIGNED) = score; 


        /* state UNALIGNED is a source for special END */ 
        temp = score + (0) + (0) ;   
        if( temp > ProteinBlockAligner_EXPL_SPECIAL(mat,i,j,END) )   {  
          ProteinBlockAligner_EXPL_SPECIAL(mat,i,j,END) = temp;  
          }  




        /* Finished calculating state UNALIGNED */ 
        }  


      /* Special state START has no special to special movements */ 


      /* Special state END has no special to special movements */ 
      }  
    stop_reporting();    
    return TRUE;     
}    


/* Function:  ProteinBlockAligner_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ProteinBlockAligner *]
 *
 */
ProteinBlockAligner * ProteinBlockAligner_alloc(void) 
{
    ProteinBlockAligner * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(ProteinBlockAligner *) ckalloc (sizeof(ProteinBlockAligner))) == NULL)  {  
      warn("ProteinBlockAligner_alloc failed "); 
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


/* Function:  free_ProteinBlockAligner(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [ProteinBlockAligner *]
 *
 * Return [UNKN ]  Undocumented return value [ProteinBlockAligner *]
 *
 */
ProteinBlockAligner * free_ProteinBlockAligner(ProteinBlockAligner * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a ProteinBlockAligner obj. Should be trappable");   
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
    /* obj->q is linked in */ 
    /* obj->t is linked in */ 
    /* obj->m is linked in */ 
    /* obj->bentry is linked in */ 
    /* obj->bexit is linked in */ 
    /* obj->bfor_trans is linked in */ 
    /* obj->b_self_trans is linked in */ 
    /* obj->b3exit is linked in */ 


    ckfree(obj); 
    return NULL; 
}    





#ifdef _cplusplus
}
#endif
