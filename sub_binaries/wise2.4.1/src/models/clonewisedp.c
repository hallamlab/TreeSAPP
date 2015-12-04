#ifdef _cplusplus
extern "C" {
#endif
#include "clonewisedp.h"

# line 5 "clonewisedp.c"


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
#define SKIP_QUERY 1 
#define SKIP_TARGET 2    


#define START 0  
#define TRUSTED_SPECIAL 1    
#define END 2    


#define CloneWise_EXPL_MATRIX(this_matrix,i,j,STATE) this_matrix->basematrix->matrix[((j+1)*3)+STATE][i+1]   
#define CloneWise_EXPL_SPECIAL(matrix,i,j,STATE) matrix->basematrix->specmatrix[STATE][j+1]  
#define CloneWise_READ_OFF_ERROR -3
 


#define CloneWise_VSMALL_MATRIX(mat,i,j,STATE) mat->basematrix->matrix[(j+2)%2][((i+1)*3)+STATE] 
#define CloneWise_VSMALL_SPECIAL(mat,i,j,STATE) mat->basematrix->specmatrix[(j+2)%2][STATE]  




#define CloneWise_SHATTER_SPECIAL(matrix,i,j,STATE) matrix->shatter->special[STATE][j]   
#define CloneWise_SHATTER_MATRIX(matrix,i,j,STATE)  fetch_cell_value_ShatterMatrix(mat->shatter,i,j,STATE)   


/* Function:  PackAln_read_Shatter_CloneWise(mat)
 *
 * Descrip:    Reads off PackAln from shatter matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [CloneWise *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * PackAln_read_Shatter_CloneWise(CloneWise * mat) 
{
    CloneWise_access_func_holder holder;     


    holder.access_main    = CloneWise_shatter_access_main;   
    holder.access_special = CloneWise_shatter_access_special;    
    assert(mat);     
    assert(mat->shatter);    
    return PackAln_read_generic_CloneWise(mat,holder);   
}    


/* Function:  CloneWise_shatter_access_main(mat,i,j,state)
 *
 * Descrip: No Description
 *
 * Arg:          mat [UNKN ] Undocumented argument [CloneWise *]
 * Arg:            i [UNKN ] Undocumented argument [int]
 * Arg:            j [UNKN ] Undocumented argument [int]
 * Arg:        state [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int CloneWise_shatter_access_main(CloneWise * mat,int i,int j,int state) 
{
    return CloneWise_SHATTER_MATRIX(mat,i,j,state);  
}    


/* Function:  CloneWise_shatter_access_special(mat,i,j,state)
 *
 * Descrip: No Description
 *
 * Arg:          mat [UNKN ] Undocumented argument [CloneWise *]
 * Arg:            i [UNKN ] Undocumented argument [int]
 * Arg:            j [UNKN ] Undocumented argument [int]
 * Arg:        state [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int CloneWise_shatter_access_special(CloneWise * mat,int i,int j,int state) 
{
    return CloneWise_SHATTER_SPECIAL(mat,i,j,state); 
}    


/* Function:  calculate_shatter_CloneWise(mat,dpenv)
 *
 * Descrip:    This function calculates the CloneWise matrix when in shatter mode
 *
 *
 * Arg:          mat [UNKN ] (null) [CloneWise *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean calculate_shatter_CloneWise(CloneWise * mat,DPEnvelope * dpenv) 
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


    start_reporting("CloneWise Matrix calculation: ");   
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




        /* For state MATCH */ 
        /* setting first movement to score */ 
        score = SIG_1_1[MATCH] + (mat->match->matrix[i][j]+1);   
        /* From state MATCH to state MATCH */ 
        temp = SIG_0_1[MATCH] + 0;   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = SIG_1_0[MATCH] + 0;   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SKIP_QUERY to state MATCH */ 
        temp = SIG_1_1[SKIP_QUERY] + 0;  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SKIP_TARGET to state MATCH */ 
        temp = SIG_1_1[SKIP_TARGET] + 0;     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* Has restricted position */ 
        if( (i-1) == 0 && (j-1) == 0 )   {  
          /* From state START to state MATCH */ 
          temp = CloneWise_SHATTER_SPECIAL(mat,i-1,j-1,START) + 0;   
          if( temp  > score )    {  
            score = temp;    
            }  
          }  


        /* Ok - finished max calculation for MATCH */ 
        /* Add any movement independant score and put away */ 
         score += mat->match->matrix[i][j];  
         SIG_0_0[MATCH] = score; 


        /* state MATCH is a source for special TRUSTED_SPECIAL */ 
        temp = score + (mat->target_special_s) + (0) ;   
        if( temp > CloneWise_SHATTER_SPECIAL(mat,i,j,TRUSTED_SPECIAL) )  {  
          CloneWise_SHATTER_SPECIAL(mat,i,j,TRUSTED_SPECIAL) = temp;     
          }  




        /* Finished calculating state MATCH */ 


        /* For state SKIP_QUERY */ 
        /* setting first movement to score */ 
        score = SIG_1_0[MATCH] + mat->query_skip_start;  
        /* From state SKIP_TARGET to state SKIP_QUERY */ 
        temp = SIG_1_0[SKIP_TARGET] + mat->query_skip_start;     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SKIP_QUERY to state SKIP_QUERY */ 
        temp = SIG_1_0[SKIP_QUERY] + 0;  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* Has restricted position */ 
        if( (i-1) == 0 && (j-0) == 0 )   {  
          /* From state START to state SKIP_QUERY */ 
          temp = CloneWise_SHATTER_SPECIAL(mat,i-1,j-0,START) + mat->query_skip_start;   
          if( temp  > score )    {  
            score = temp;    
            }  
          }  


        /* Ok - finished max calculation for SKIP_QUERY */ 
        /* Add any movement independant score and put away */ 
         score += mat->match->skip_iset[i];  
         SIG_0_0[SKIP_QUERY] = score;    


        /* state SKIP_QUERY is a source for special TRUSTED_SPECIAL */ 
        temp = score + (mat->target_special_s) + (0) ;   
        if( temp > CloneWise_SHATTER_SPECIAL(mat,i,j,TRUSTED_SPECIAL) )  {  
          CloneWise_SHATTER_SPECIAL(mat,i,j,TRUSTED_SPECIAL) = temp;     
          }  




        /* Finished calculating state SKIP_QUERY */ 


        /* For state SKIP_TARGET */ 
        /* setting first movement to score */ 
        score = SIG_0_1[MATCH] + mat->target_skip_start;     
        /* From state SKIP_QUERY to state SKIP_TARGET */ 
        temp = SIG_0_1[SKIP_QUERY] + mat->target_skip_start;     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SKIP_TARGET to state SKIP_TARGET */ 
        temp = SIG_0_1[SKIP_TARGET] + 0;     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* Has restricted position */ 
        if( (i-0) == 0 && (j-1) == 0 )   {  
          /* From state START to state SKIP_TARGET */ 
          temp = CloneWise_SHATTER_SPECIAL(mat,i-0,j-1,START) + mat->target_skip_start;  
          if( temp  > score )    {  
            score = temp;    
            }  
          }  


        /* Ok - finished max calculation for SKIP_TARGET */ 
        /* Add any movement independant score and put away */ 
         score += mat->match->skip_jset[j];  
         SIG_0_0[SKIP_TARGET] = score;   


        /* state SKIP_TARGET is a source for special TRUSTED_SPECIAL */ 
        temp = score + (mat->target_special_s) + (0) ;   
        if( temp > CloneWise_SHATTER_SPECIAL(mat,i,j,TRUSTED_SPECIAL) )  {  
          CloneWise_SHATTER_SPECIAL(mat,i,j,TRUSTED_SPECIAL) = temp;     
          }  




        /* Finished calculating state SKIP_TARGET */ 
        }  


      /* Special state START has no special to special movements */ 


      /* Special state TRUSTED_SPECIAL has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = CloneWise_SHATTER_SPECIAL(mat,0,j,TRUSTED_SPECIAL);    


      /* Source START is a special source for TRUSTED_SPECIAL */ 
      /* Has restricted position */ 
      if( (j-1) == 0  )  {  
        temp = CloneWise_SHATTER_SPECIAL(mat,0,j - 1,START) + (0) + (0);     
        if( temp > score )   
          score = temp;  
        }  


      /* Source TRUSTED_SPECIAL is a special source for TRUSTED_SPECIAL */ 
      temp = CloneWise_SHATTER_SPECIAL(mat,0,j - 1,TRUSTED_SPECIAL) + ((mat->match->skip_jset[j]-1)) + (0);  
      if( temp > score ) 
        score = temp;    


      /* Source MATCH for state TRUSTED_SPECIAL is not special... already calculated */ 
      /* Source SKIP_QUERY for state TRUSTED_SPECIAL is not special... already calculated */ 
      /* Source SKIP_TARGET for state TRUSTED_SPECIAL is not special... already calculated */ 
      /* Put back score... (now updated!) */ 
      CloneWise_SHATTER_SPECIAL(mat,0,j,TRUSTED_SPECIAL) = score;    
      /* Finished updating state TRUSTED_SPECIAL */ 




      /* Special state END has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = CloneWise_SHATTER_SPECIAL(mat,0,j,END);    


      /* Source TRUSTED_SPECIAL is a special source for END */ 
      /* Has restricted position */ 
      if( j == mat->lenj-1 ) {  
        temp = CloneWise_SHATTER_SPECIAL(mat,0,j - 1,TRUSTED_SPECIAL) + (0) + (0);   
        if( temp > score )   
          score = temp;  
        }  


      /* Put back score... (now updated!) */ 
      CloneWise_SHATTER_SPECIAL(mat,0,j,END) = score;    
      /* Finished updating state END */ 


      }  
    stop_reporting();    
    return TRUE;     
}    


/* Function:  search_CloneWise(dbsi,out,q,t,match,target_skip_start,target_skip,query_skip_start,query_skip,spread,target_special_s)
 *
 * Descrip:    This function makes a database search of CloneWise
 *             It uses the dbsi structure to choose which implementation
 *             to use of the database searching. This way at run time you
 *             can switch between single threaded/multi-threaded or hardware
 *
 *
 * Arg:                     dbsi [UNKN ] Undocumented argument [DBSearchImpl *]
 * Arg:                      out [UNKN ] Undocumented argument [Hscore *]
 * Arg:                        q [UNKN ] Undocumented argument [MappedCloneSet *]
 * Arg:                        t [UNKN ] Undocumented argument [MappedCloneSet *]
 * Arg:                    match [UNKN ] Undocumented argument [MappedCloneMatch*]
 * Arg:        target_skip_start [UNKN ] Undocumented argument [Score]
 * Arg:              target_skip [UNKN ] Undocumented argument [Score]
 * Arg:         query_skip_start [UNKN ] Undocumented argument [Score]
 * Arg:               query_skip [UNKN ] Undocumented argument [Score]
 * Arg:                   spread [UNKN ] Undocumented argument [int]
 * Arg:         target_special_s [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type search_CloneWise(DBSearchImpl * dbsi,Hscore * out,MappedCloneSet * q,MappedCloneSet * t ,MappedCloneMatch* match,Score target_skip_start,Score target_skip,Score query_skip_start,Score query_skip,int spread,int target_special_s) 
{
#ifdef PTHREAD   
    int i;   
    int thr_no;  
    pthread_attr_t pat;  
    struct thread_pool_holder_CloneWise * holder;    
#endif   
    if( out == NULL )    {  
      warn("Passed in a null Hscore object into search_CloneWise. Can't process results!");  
      return SEARCH_ERROR;   
      }  
    if( dbsi == NULL )   {  
      warn("Passed in a null DBSearchImpl object into search_CloneWise. Can't process results!");    
      return SEARCH_ERROR;   
      }  
    if( dbsi->trace_level > 5 )  
      warn("Asking for trace level of %d in database search for CloneWise, but it was compiled with a trace level of -2139062144. Not all trace statements can be shown",dbsi->trace_level); 
    switch(dbsi->type)   { /*switch on implementation*/ 
      case DBSearchImpl_Serial : 
        return serial_search_CloneWise(out,q,t ,match,target_skip_start,target_skip,query_skip_start,query_skip,spread,target_special_s);    
      case DBSearchImpl_Pthreads :   
#ifdef PTHREAD   
        holder = (struct thread_pool_holder_CloneWise *) ckalloc(sizeof(struct thread_pool_holder_CloneWise));   
        if( holder == NULL )     {  
          warn("Unable to allocated thread pool datastructure...");  
          return SEARCH_ERROR;   
          }  
        holder->out = out;   
        holder->dbsi = dbsi; 
        holder->q = q;   
        holder->t = t;   
        holder->match = match;   
        holder->target_skip_start = target_skip_start;   
        holder->target_skip = target_skip;   
        holder->query_skip_start = query_skip_start; 
        holder->query_skip = query_skip; 
        holder->spread = spread; 
        holder->target_special_s = target_special_s; 
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
          if( pthread_create(holder->pool+i,&pat,thread_loop_CloneWise,(void *)holder) ) 
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
        warn("You did not specifiy the PTHREAD compile when compiled the C code for CloneWise"); 
#endif /* finished threads */    
      default :  
        warn("database search implementation %s was not provided in the compiled dynamite file from CloneWise",impl_string_DBSearchImpl(dbsi));  
        return SEARCH_ERROR; 
      } /* end of switch on implementation */ 


}    


/* Function:  thread_loop_CloneWise(ptr)
 *
 * Descrip:    dummy loop code foreach thread for CloneWise
 *
 *
 * Arg:        ptr [UNKN ] Undocumented argument [void *]
 *
 * Return [UNKN ]  Undocumented return value [void *]
 *
 */
void * thread_loop_CloneWise(void * ptr) 
{
    fatal("dummy thread loop function"); 
}    


/* Function:  serial_search_CloneWise(out,q,t,match,target_skip_start,target_skip,query_skip_start,query_skip,spread,target_special_s)
 *
 * Descrip:    This function makes a database search of CloneWise
 *             It is a single processor implementation
 *
 *
 * Arg:                      out [UNKN ] Undocumented argument [Hscore *]
 * Arg:                        q [UNKN ] Undocumented argument [MappedCloneSet *]
 * Arg:                        t [UNKN ] Undocumented argument [MappedCloneSet *]
 * Arg:                    match [UNKN ] Undocumented argument [MappedCloneMatch*]
 * Arg:        target_skip_start [UNKN ] Undocumented argument [Score]
 * Arg:              target_skip [UNKN ] Undocumented argument [Score]
 * Arg:         query_skip_start [UNKN ] Undocumented argument [Score]
 * Arg:               query_skip [UNKN ] Undocumented argument [Score]
 * Arg:                   spread [UNKN ] Undocumented argument [int]
 * Arg:         target_special_s [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type serial_search_CloneWise(Hscore * out,MappedCloneSet * q,MappedCloneSet * t ,MappedCloneMatch* match,Score target_skip_start,Score target_skip,Score query_skip_start,Score query_skip,int spread,int target_special_s) 
{
    int db_status;   
    int score;   
    int query_pos = 0;   
    int target_pos = 0;  
    DataScore * ds;  


    push_errormsg_stack("Before any actual search in db searching"); 


    target_pos = 0;  




    /* No maximum length - allocated on-the-fly */ 
    score = score_only_CloneWise(q, t , match, target_skip_start, target_skip, query_skip_start, query_skip, spread, target_special_s);  
    if( should_store_Hscore(out,score) == TRUE )     { /*if storing datascore*/ 
      ds = new_DataScore_from_storage(out);  
      if( ds == NULL )   {  
        warn("CloneWise search had a memory error in allocating a new_DataScore (?a leak somewhere - DataScore is a very small datastructure");  
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


/* Function:  score_only_CloneWise(q,t,match,target_skip_start,target_skip,query_skip_start,query_skip,spread,target_special_s)
 *
 * Descrip:    This function just calculates the score for the matrix
 *             I am pretty sure we can do this better, but hey, for the moment...
 *             It calls /allocate_CloneWise_only
 *
 *
 * Arg:                        q [UNKN ] query data structure [MappedCloneSet *]
 * Arg:                        t [UNKN ] target data structure [MappedCloneSet *]
 * Arg:                    match [UNKN ] Resource [MappedCloneMatch*]
 * Arg:        target_skip_start [UNKN ] Resource [Score]
 * Arg:              target_skip [UNKN ] Resource [Score]
 * Arg:         query_skip_start [UNKN ] Resource [Score]
 * Arg:               query_skip [UNKN ] Resource [Score]
 * Arg:                   spread [UNKN ] Resource [int]
 * Arg:         target_special_s [UNKN ] Resource [int]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int score_only_CloneWise(MappedCloneSet * q,MappedCloneSet * t ,MappedCloneMatch* match,Score target_skip_start,Score target_skip,Score query_skip_start,Score query_skip,int spread,int target_special_s) 
{
    int bestscore = NEGI;    
    int i;   
    int j;   
    int k;   
    CloneWise * mat;     


    mat = allocate_CloneWise_only(q, t , match, target_skip_start, target_skip, query_skip_start, query_skip, spread, target_special_s); 
    if( mat == NULL )    {  
      warn("Memory allocation error in the db search - unable to communicate to calling function. this spells DIASTER!");    
      return NEGI;   
      }  
    if((mat->basematrix = BaseMatrix_alloc_matrix_and_specials(2,(mat->leni + 1) * 3,2,3)) == NULL)  {  
      warn("Score only matrix for CloneWise cannot be allocated, (asking for 1  by %d  cells)",mat->leni*3); 
      mat = free_CloneWise(mat);     
      return 0;  
      }  
    mat->basematrix->type = BASEMATRIX_TYPE_VERYSMALL;   


    /* Now, initiate matrix */ 
    for(j=0;j<3;j++) {  
      for(i=(-1);i<mat->leni;i++)    {  
        for(k=0;k<3;k++) 
          CloneWise_VSMALL_MATRIX(mat,i,j,k) = NEGI; 
        }  
      CloneWise_VSMALL_SPECIAL(mat,i,j,START) = 0;   
      CloneWise_VSMALL_SPECIAL(mat,i,j,TRUSTED_SPECIAL) = NEGI;  
      CloneWise_VSMALL_SPECIAL(mat,i,j,END) = NEGI;  
      }  


    /* Ok, lets do-o-o-o-o it */ 


    for(j=0;j<mat->lenj;j++) { /*for all target positions*/ 
      auto int score;    
      auto int temp;     
      for(i=0;i<mat->leni;i++)   { /*for all query positions*/ 


        /* For state MATCH */ 
        /* setting first movement to score */ 
        score = CloneWise_VSMALL_MATRIX(mat,i-1,j-1,MATCH) + (mat->match->matrix[i][j]+1);   
        /* From state MATCH to state MATCH */ 
        temp = CloneWise_VSMALL_MATRIX(mat,i-0,j-1,MATCH) + 0;   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = CloneWise_VSMALL_MATRIX(mat,i-1,j-0,MATCH) + 0;   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SKIP_QUERY to state MATCH */ 
        temp = CloneWise_VSMALL_MATRIX(mat,i-1,j-1,SKIP_QUERY) + 0;  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SKIP_TARGET to state MATCH */ 
        temp = CloneWise_VSMALL_MATRIX(mat,i-1,j-1,SKIP_TARGET) + 0;     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* Has restricted position */ 
        if( (i-1) == 0 && (j-1) == 0 )   {  
          /* From state START to state MATCH */ 
          temp = CloneWise_VSMALL_SPECIAL(mat,i-1,j-1,START) + 0;    
          if( temp  > score )    {  
            score = temp;    
            }  
          }  


        /* Ok - finished max calculation for MATCH */ 
        /* Add any movement independant score and put away */ 
         score += mat->match->matrix[i][j];  
         CloneWise_VSMALL_MATRIX(mat,i,j,MATCH) = score; 


        /* state MATCH is a source for special TRUSTED_SPECIAL */ 
        temp = score + (mat->target_special_s) + (0) ;   
        if( temp > CloneWise_VSMALL_SPECIAL(mat,i,j,TRUSTED_SPECIAL) )   {  
          CloneWise_VSMALL_SPECIAL(mat,i,j,TRUSTED_SPECIAL) = temp;  
          }  




        /* Finished calculating state MATCH */ 


        /* For state SKIP_QUERY */ 
        /* setting first movement to score */ 
        score = CloneWise_VSMALL_MATRIX(mat,i-1,j-0,MATCH) + mat->query_skip_start;  
        /* From state SKIP_TARGET to state SKIP_QUERY */ 
        temp = CloneWise_VSMALL_MATRIX(mat,i-1,j-0,SKIP_TARGET) + mat->query_skip_start;     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SKIP_QUERY to state SKIP_QUERY */ 
        temp = CloneWise_VSMALL_MATRIX(mat,i-1,j-0,SKIP_QUERY) + 0;  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* Has restricted position */ 
        if( (i-1) == 0 && (j-0) == 0 )   {  
          /* From state START to state SKIP_QUERY */ 
          temp = CloneWise_VSMALL_SPECIAL(mat,i-1,j-0,START) + mat->query_skip_start;    
          if( temp  > score )    {  
            score = temp;    
            }  
          }  


        /* Ok - finished max calculation for SKIP_QUERY */ 
        /* Add any movement independant score and put away */ 
         score += mat->match->skip_iset[i];  
         CloneWise_VSMALL_MATRIX(mat,i,j,SKIP_QUERY) = score;    


        /* state SKIP_QUERY is a source for special TRUSTED_SPECIAL */ 
        temp = score + (mat->target_special_s) + (0) ;   
        if( temp > CloneWise_VSMALL_SPECIAL(mat,i,j,TRUSTED_SPECIAL) )   {  
          CloneWise_VSMALL_SPECIAL(mat,i,j,TRUSTED_SPECIAL) = temp;  
          }  




        /* Finished calculating state SKIP_QUERY */ 


        /* For state SKIP_TARGET */ 
        /* setting first movement to score */ 
        score = CloneWise_VSMALL_MATRIX(mat,i-0,j-1,MATCH) + mat->target_skip_start;     
        /* From state SKIP_QUERY to state SKIP_TARGET */ 
        temp = CloneWise_VSMALL_MATRIX(mat,i-0,j-1,SKIP_QUERY) + mat->target_skip_start;     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SKIP_TARGET to state SKIP_TARGET */ 
        temp = CloneWise_VSMALL_MATRIX(mat,i-0,j-1,SKIP_TARGET) + 0;     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* Has restricted position */ 
        if( (i-0) == 0 && (j-1) == 0 )   {  
          /* From state START to state SKIP_TARGET */ 
          temp = CloneWise_VSMALL_SPECIAL(mat,i-0,j-1,START) + mat->target_skip_start;   
          if( temp  > score )    {  
            score = temp;    
            }  
          }  


        /* Ok - finished max calculation for SKIP_TARGET */ 
        /* Add any movement independant score and put away */ 
         score += mat->match->skip_jset[j];  
         CloneWise_VSMALL_MATRIX(mat,i,j,SKIP_TARGET) = score;   


        /* state SKIP_TARGET is a source for special TRUSTED_SPECIAL */ 
        temp = score + (mat->target_special_s) + (0) ;   
        if( temp > CloneWise_VSMALL_SPECIAL(mat,i,j,TRUSTED_SPECIAL) )   {  
          CloneWise_VSMALL_SPECIAL(mat,i,j,TRUSTED_SPECIAL) = temp;  
          }  




        /* Finished calculating state SKIP_TARGET */ 
        } /* end of for all query positions */ 




      /* Special state START has no special to special movements */ 


      /* Special state TRUSTED_SPECIAL has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = CloneWise_VSMALL_SPECIAL(mat,0,j,TRUSTED_SPECIAL); 


      /* Source START is a special source for TRUSTED_SPECIAL */ 
      /* Has restricted position */ 
      if( (j-1) == 0  )  {  
        temp = CloneWise_VSMALL_SPECIAL(mat,0,j - 1,START) + (0) + (0);  
        if( temp > score )   
          score = temp;  
        }  


      /* Source TRUSTED_SPECIAL is a special source for TRUSTED_SPECIAL */ 
      temp = CloneWise_VSMALL_SPECIAL(mat,0,j - 1,TRUSTED_SPECIAL) + ((mat->match->skip_jset[j]-1)) + (0);   
      if( temp > score ) 
        score = temp;    


      /* Source MATCH for state TRUSTED_SPECIAL is not special... already calculated */ 
      /* Source SKIP_QUERY for state TRUSTED_SPECIAL is not special... already calculated */ 
      /* Source SKIP_TARGET for state TRUSTED_SPECIAL is not special... already calculated */ 
      /* Put back score... (now updated!) */ 
      CloneWise_VSMALL_SPECIAL(mat,0,j,TRUSTED_SPECIAL) = score; 
      /* Finished updating state TRUSTED_SPECIAL */ 




      /* Special state END has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = CloneWise_VSMALL_SPECIAL(mat,0,j,END); 


      /* Source TRUSTED_SPECIAL is a special source for END */ 
      /* Has restricted position */ 
      if( j == mat->lenj-1 ) {  
        temp = CloneWise_VSMALL_SPECIAL(mat,0,j - 1,TRUSTED_SPECIAL) + (0) + (0);    
        if( temp > score )   
          score = temp;  
        }  


      /* Put back score... (now updated!) */ 
      CloneWise_VSMALL_SPECIAL(mat,0,j,END) = score; 
      /* Finished updating state END */ 


      if( bestscore < CloneWise_VSMALL_SPECIAL(mat,0,j,END) )    
        bestscore = CloneWise_VSMALL_SPECIAL(mat,0,j,END);   
      } /* end of for all target positions */ 


    mat = free_CloneWise(mat);   
    return bestscore;    
}    


/* Function:  PackAln_bestmemory_CloneWise(q,t,match,target_skip_start,target_skip,query_skip_start,query_skip,spread,target_special_s,dpenv,dpri)
 *
 * Descrip:    This function chooses the best memory set-up for the alignment
 *             using calls to basematrix, and then implements either a large
 *             or small memory model.
 *
 *             It is the best function to use if you just want an alignment
 *
 *             If you want a label alignment, you will need
 *             /convert_PackAln_to_AlnBlock_CloneWise
 *
 *
 * Arg:                        q [UNKN ] query data structure [MappedCloneSet *]
 * Arg:                        t [UNKN ] target data structure [MappedCloneSet *]
 * Arg:                    match [UNKN ] Resource [MappedCloneMatch*]
 * Arg:        target_skip_start [UNKN ] Resource [Score]
 * Arg:              target_skip [UNKN ] Resource [Score]
 * Arg:         query_skip_start [UNKN ] Resource [Score]
 * Arg:               query_skip [UNKN ] Resource [Score]
 * Arg:                   spread [UNKN ] Resource [int]
 * Arg:         target_special_s [UNKN ] Resource [int]
 * Arg:                    dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:                     dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * PackAln_bestmemory_CloneWise(MappedCloneSet * q,MappedCloneSet * t ,MappedCloneMatch* match,Score target_skip_start,Score target_skip,Score query_skip_start,Score query_skip,int spread,int target_special_s,DPEnvelope * dpenv,DPRunImpl * dpri) 
{
    long long total; 
    CloneWise * mat; 
    PackAln * out;   
    DebugMatrix * de;    
    DPRunImplMemory strategy;    
    assert(dpri);    


    total = q->length * t->length;   
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
        if( (mat=allocate_Expl_CloneWise(q, t , match, target_skip_start, target_skip, query_skip_start, query_skip, spread, target_special_s,dpri)) == NULL )   {  
          warn("Unable to allocate large CloneWise version");    
          return NULL;   
          }  
        calculate_dpenv_CloneWise(mat,dpenv);    
        out =  PackAln_read_Expl_CloneWise(mat); 
        }  
      else   {  
        mat = allocate_CloneWise_only(q, t , match, target_skip_start, target_skip, query_skip_start, query_skip, spread, target_special_s);     
        calculate_shatter_CloneWise(mat,dpenv);  
        out = PackAln_read_Shatter_CloneWise(mat);   
        }  
      }  
    else {  
      if( strategy == DPIM_Linear )  {  
        /* use small implementation */ 
        if( (mat=allocate_Small_CloneWise(q, t , match, target_skip_start, target_skip, query_skip_start, query_skip, spread, target_special_s)) == NULL )   {  
          warn("Unable to allocate small CloneWise version");    
          return NULL;   
          }  
        out = PackAln_calculate_Small_CloneWise(mat,dpenv);  
        }  
      else   {  
        /* use Large implementation */ 
        if( (mat=allocate_Expl_CloneWise(q, t , match, target_skip_start, target_skip, query_skip_start, query_skip, spread, target_special_s,dpri)) == NULL )   {  
          warn("Unable to allocate large CloneWise version");    
          return NULL;   
          }  
        if( dpri->debug == TRUE) {  
          fatal("Asked for dydebug, but dynamite file not compiled with -g. Need to recompile dynamite source"); 
          }  
        else {  
          calculate_CloneWise(mat);  
          out =  PackAln_read_Expl_CloneWise(mat);   
          }  
        }  
      }  


    mat = free_CloneWise(mat);   
    return out;  
}    


/* Function:  allocate_CloneWise_only(q,t,match,target_skip_start,target_skip,query_skip_start,query_skip,spread,target_special_s)
 *
 * Descrip:    This function only allocates the CloneWise structure
 *             checks types where possible and determines leni and lenj
 *             The basematrix area is delt with elsewhere
 *
 *
 * Arg:                        q [UNKN ] query data structure [MappedCloneSet *]
 * Arg:                        t [UNKN ] target data structure [MappedCloneSet *]
 * Arg:                    match [UNKN ] Resource [MappedCloneMatch*]
 * Arg:        target_skip_start [UNKN ] Resource [Score]
 * Arg:              target_skip [UNKN ] Resource [Score]
 * Arg:         query_skip_start [UNKN ] Resource [Score]
 * Arg:               query_skip [UNKN ] Resource [Score]
 * Arg:                   spread [UNKN ] Resource [int]
 * Arg:         target_special_s [UNKN ] Resource [int]
 *
 * Return [UNKN ]  Undocumented return value [CloneWise *]
 *
 */
CloneWise * allocate_CloneWise_only(MappedCloneSet * q,MappedCloneSet * t ,MappedCloneMatch* match,Score target_skip_start,Score target_skip,Score query_skip_start,Score query_skip,int spread,int target_special_s) 
{
    CloneWise * out;     


    if((out= CloneWise_alloc()) == NULL) {  
      warn("Allocation of basic CloneWise structure failed..."); 
      return NULL;   
      }  


    out->q = q;  
    out->t = t;  
    out->match = match;  
    out->target_skip_start = target_skip_start;  
    out->target_skip = target_skip;  
    out->query_skip_start = query_skip_start;    
    out->query_skip = query_skip;    
    out->spread = spread;    
    out->target_special_s = target_special_s;    
    out->leni = q->length;   
    out->lenj = t->length;   
    return out;  
}    


/* Function:  allocate_Expl_CloneWise(q,t,match,target_skip_start,target_skip,query_skip_start,query_skip,spread,target_special_s,dpri)
 *
 * Descrip:    This function allocates the CloneWise structure
 *             and the basematrix area for explicit memory implementations
 *             It calls /allocate_CloneWise_only
 *
 *
 * Arg:                        q [UNKN ] query data structure [MappedCloneSet *]
 * Arg:                        t [UNKN ] target data structure [MappedCloneSet *]
 * Arg:                    match [UNKN ] Resource [MappedCloneMatch*]
 * Arg:        target_skip_start [UNKN ] Resource [Score]
 * Arg:              target_skip [UNKN ] Resource [Score]
 * Arg:         query_skip_start [UNKN ] Resource [Score]
 * Arg:               query_skip [UNKN ] Resource [Score]
 * Arg:                   spread [UNKN ] Resource [int]
 * Arg:         target_special_s [UNKN ] Resource [int]
 * Arg:                     dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [CloneWise *]
 *
 */
CloneWise * allocate_Expl_CloneWise(MappedCloneSet * q,MappedCloneSet * t ,MappedCloneMatch* match,Score target_skip_start,Score target_skip,Score query_skip_start,Score query_skip,int spread,int target_special_s,DPRunImpl * dpri) 
{
    CloneWise * out; 


    out = allocate_CloneWise_only(q, t , match, target_skip_start, target_skip, query_skip_start, query_skip, spread, target_special_s); 
    if( out == NULL )    
      return NULL;   
    if( dpri->should_cache == TRUE ) {  
      if( dpri->cache != NULL )  {  
        if( dpri->cache->maxleni >= (out->lenj+1)*3 && dpri->cache->maxlenj >= (out->leni+1))    
          out->basematrix = hard_link_BaseMatrix(dpri->cache);   
        else 
          dpri->cache = free_BaseMatrix(dpri->cache);    
        }  
      }  
    if( out->basematrix == NULL )    {  
      if( (out->basematrix = BaseMatrix_alloc_matrix_and_specials((out->lenj+1)*3,(out->leni+1),3,out->lenj+1)) == NULL) {  
        warn("Explicit matrix CloneWise cannot be allocated, (asking for %d by %d main cells)",out->leni,out->lenj); 
        free_CloneWise(out);     
        return NULL; 
        }  
      }  
    if( dpri->should_cache == TRUE && dpri->cache == NULL)   
      dpri->cache = hard_link_BaseMatrix(out->basematrix);   
    out->basematrix->type = BASEMATRIX_TYPE_EXPLICIT;    
    init_CloneWise(out);     
    return out;  
}    


/* Function:  init_CloneWise(mat)
 *
 * Descrip:    This function initates CloneWise matrix when in explicit mode
 *             Called in /allocate_Expl_CloneWise
 *
 *
 * Arg:        mat [UNKN ] CloneWise which contains explicit basematrix memory [CloneWise *]
 *
 */
void init_CloneWise(CloneWise * mat) 
{
    register int i;  
    register int j;  
    if( mat->basematrix->type != BASEMATRIX_TYPE_EXPLICIT)   {  
      warn("Cannot iniate matrix, is not an explicit memory type and you have assummed that");   
      return;    
      }  


    for(i= (-1);i<mat->q->length;i++)    {  
      for(j= (-1);j<2;j++)   {  
        CloneWise_EXPL_MATRIX(mat,i,j,MATCH) = NEGI; 
        CloneWise_EXPL_MATRIX(mat,i,j,SKIP_QUERY) = NEGI;    
        CloneWise_EXPL_MATRIX(mat,i,j,SKIP_TARGET) = NEGI;   
        }  
      }  
    for(j= (-1);j<mat->t->length;j++)    {  
      for(i= (-1);i<2;i++)   {  
        CloneWise_EXPL_MATRIX(mat,i,j,MATCH) = NEGI; 
        CloneWise_EXPL_MATRIX(mat,i,j,SKIP_QUERY) = NEGI;    
        CloneWise_EXPL_MATRIX(mat,i,j,SKIP_TARGET) = NEGI;   
        }  
      CloneWise_EXPL_SPECIAL(mat,i,j,START) = 0; 
      CloneWise_EXPL_SPECIAL(mat,i,j,TRUSTED_SPECIAL) = NEGI;    
      CloneWise_EXPL_SPECIAL(mat,i,j,END) = NEGI;    
      }  
    return;  
}    


/* Function:  recalculate_PackAln_CloneWise(pal,mat)
 *
 * Descrip:    This function recalculates the PackAln structure produced by CloneWise
 *             For example, in linear space methods this is used to score them
 *
 *
 * Arg:        pal [UNKN ] Undocumented argument [PackAln *]
 * Arg:        mat [UNKN ] Undocumented argument [CloneWise *]
 *
 */
void recalculate_PackAln_CloneWise(PackAln * pal,CloneWise * mat) 
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
          if( offi == 1 && offj == 1 && prev->state == MATCH )   {  
            pau->score = (mat->match->matrix[i][j]+1) + (mat->match->matrix[i][j]);  
            continue;    
            }  
          if( offi == 0 && offj == 1 && prev->state == MATCH )   {  
            pau->score = 0 + (mat->match->matrix[i][j]);     
            continue;    
            }  
          if( offi == 1 && offj == 0 && prev->state == MATCH )   {  
            pau->score = 0 + (mat->match->matrix[i][j]);     
            continue;    
            }  
          if( offi == 1 && offj == 1 && prev->state == SKIP_QUERY )  {  
            pau->score = 0 + (mat->match->matrix[i][j]);     
            continue;    
            }  
          if( offi == 1 && offj == 1 && prev->state == SKIP_TARGET ) {  
            pau->score = 0 + (mat->match->matrix[i][j]);     
            continue;    
            }  
          if( offj == 1 && prev->state == (START+3) )    {  
            pau->score = 0 + (mat->match->matrix[i][j]);     
            continue;    
            }  
          warn("In recaluclating PackAln with state MATCH, from [%d,%d,%d], got a bad source state. Error!",offi,offj,prev->state);  
          break; 
        case SKIP_QUERY :    
          if( offi == 1 && offj == 0 && prev->state == MATCH )   {  
            pau->score = mat->query_skip_start + (mat->match->skip_iset[i]);     
            continue;    
            }  
          if( offi == 1 && offj == 0 && prev->state == SKIP_TARGET ) {  
            pau->score = mat->query_skip_start + (mat->match->skip_iset[i]);     
            continue;    
            }  
          if( offi == 1 && offj == 0 && prev->state == SKIP_QUERY )  {  
            pau->score = 0 + (mat->match->skip_iset[i]);     
            continue;    
            }  
          if( offj == 0 && prev->state == (START+3) )    {  
            pau->score = mat->query_skip_start + (mat->match->skip_iset[i]);     
            continue;    
            }  
          warn("In recaluclating PackAln with state SKIP_QUERY, from [%d,%d,%d], got a bad source state. Error!",offi,offj,prev->state); 
          break; 
        case SKIP_TARGET :   
          if( offi == 0 && offj == 1 && prev->state == MATCH )   {  
            pau->score = mat->target_skip_start + (mat->match->skip_jset[j]);    
            continue;    
            }  
          if( offi == 0 && offj == 1 && prev->state == SKIP_QUERY )  {  
            pau->score = mat->target_skip_start + (mat->match->skip_jset[j]);    
            continue;    
            }  
          if( offi == 0 && offj == 1 && prev->state == SKIP_TARGET ) {  
            pau->score = 0 + (mat->match->skip_jset[j]);     
            continue;    
            }  
          if( offj == 1 && prev->state == (START+3) )    {  
            pau->score = mat->target_skip_start + (mat->match->skip_jset[j]);    
            continue;    
            }  
          warn("In recaluclating PackAln with state SKIP_TARGET, from [%d,%d,%d], got a bad source state. Error!",offi,offj,prev->state);    
          break; 
        case (START+3) :     
          warn("In recaluclating PackAln with state START, got a bad source state. Error!"); 
          break; 
        case (TRUSTED_SPECIAL+3) :   
          if( offj == 1 && prev->state == (START+3) )    {  
            pau->score = 0 + (0);    
            continue;    
            }  
          if( offj == 1 && prev->state == (TRUSTED_SPECIAL+3) )  {  
            pau->score = (mat->match->skip_jset[j]-1) + (0);     
            continue;    
            }  
          if( offj == 0 && prev->state == MATCH )    {  
            /* i here comes from the previous state ;) - not the real one */ 
            i = prev->i; 
            pau->score = mat->target_special_s + (0);    
            continue;    
            }  
          if( offj == 0 && prev->state == SKIP_QUERY )   {  
            /* i here comes from the previous state ;) - not the real one */ 
            i = prev->i; 
            pau->score = mat->target_special_s + (0);    
            continue;    
            }  
          if( offj == 0 && prev->state == SKIP_TARGET )  {  
            /* i here comes from the previous state ;) - not the real one */ 
            i = prev->i; 
            pau->score = mat->target_special_s + (0);    
            continue;    
            }  
          warn("In recaluclating PackAln with state TRUSTED_SPECIAL, got a bad source state. Error!");   
          break; 
        case (END+3) :   
          if( offj == 1 && prev->state == (TRUSTED_SPECIAL+3) )  {  
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
#define CloneWise_HIDDEN_MATRIX(thismatrix,i,j,state) (thismatrix->basematrix->matrix[(j-hiddenj+1)][(i+1)*3+state]) 
#define CloneWise_DC_SHADOW_MATRIX(thismatrix,i,j,state) (thismatrix->basematrix->matrix[((j+2)*8) % 16][(i+1)*3+state]) 
#define CloneWise_HIDDEN_SPECIAL(thismatrix,i,j,state) (thismatrix->basematrix->specmatrix[state][(j+1)])    
#define CloneWise_DC_SHADOW_SPECIAL(thismatrix,i,j,state) (thismatrix->basematrix->specmatrix[state*8][(j+1)])   
#define CloneWise_DC_SHADOW_MATRIX_SP(thismatrix,i,j,state,shadow) (thismatrix->basematrix->matrix[((((j+2)*8)+(shadow+1)) % 16)][(i+1)*3 + state])  
#define CloneWise_DC_SHADOW_SPECIAL_SP(thismatrix,i,j,state,shadow) (thismatrix->basematrix->specmatrix[state*8 +shadow+1][(j+1)])   
#define CloneWise_DC_OPT_SHADOW_MATRIX(thismatrix,i,j,state) (score_pointers[(((j+1)% 1) * (leni+1) * 3) + ((i+1) * 3) + (state)])   
#define CloneWise_DC_OPT_SHADOW_MATRIX_SP(thismatrix,i,j,state,shadow) (shadow_pointers[(((j+1)% 1) * (leni+1) * 24) + ((i+1) * 24) + (state * 8) + shadow+1])   
#define CloneWise_DC_OPT_SHADOW_SPECIAL(thismatrix,i,j,state) (thismatrix->basematrix->specmatrix[state*8][(j+1)])   
/* Function:  allocate_Small_CloneWise(q,t,match,target_skip_start,target_skip,query_skip_start,query_skip,spread,target_special_s)
 *
 * Descrip:    This function allocates the CloneWise structure
 *             and the basematrix area for a small memory implementations
 *             It calls /allocate_CloneWise_only
 *
 *
 * Arg:                        q [UNKN ] query data structure [MappedCloneSet *]
 * Arg:                        t [UNKN ] target data structure [MappedCloneSet *]
 * Arg:                    match [UNKN ] Resource [MappedCloneMatch*]
 * Arg:        target_skip_start [UNKN ] Resource [Score]
 * Arg:              target_skip [UNKN ] Resource [Score]
 * Arg:         query_skip_start [UNKN ] Resource [Score]
 * Arg:               query_skip [UNKN ] Resource [Score]
 * Arg:                   spread [UNKN ] Resource [int]
 * Arg:         target_special_s [UNKN ] Resource [int]
 *
 * Return [UNKN ]  Undocumented return value [CloneWise *]
 *
 */
#define CloneWise_DC_OPT_SHADOW_SPECIAL_SP(thismatrix,i,j,state,shadow) (thismatrix->basematrix->specmatrix[state*8 +shadow+1][(j+1)])   
CloneWise * allocate_Small_CloneWise(MappedCloneSet * q,MappedCloneSet * t ,MappedCloneMatch* match,Score target_skip_start,Score target_skip,Score query_skip_start,Score query_skip,int spread,int target_special_s) 
{
    CloneWise * out; 


    out = allocate_CloneWise_only(q, t , match, target_skip_start, target_skip, query_skip_start, query_skip, spread, target_special_s); 
    if( out == NULL )    
      return NULL;   
    out->basematrix = BaseMatrix_alloc_matrix_and_specials(16,(out->leni + 1) * 3,24,out->lenj+1);   
    if(out == NULL)  {  
      warn("Small shadow matrix CloneWise cannot be allocated, (asking for 2 by %d main cells)",out->leni+2);    
      free_CloneWise(out);   
      return NULL;   
      }  
    out->basematrix->type = BASEMATRIX_TYPE_SHADOW;  
    return out;  
}    


/* Function:  PackAln_calculate_Small_CloneWise(mat,dpenv)
 *
 * Descrip:    This function calculates an alignment for CloneWise structure in linear space
 *             If you want only the start/end points
 *             use /AlnRangeSet_calculate_Small_CloneWise 
 *
 *             The function basically
 *               finds start/end points 
 *               foreach start/end point 
 *                 calls /full_dc_CloneWise 
 *
 *
 * Arg:          mat [UNKN ] Undocumented argument [CloneWise *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * PackAln_calculate_Small_CloneWise(CloneWise * mat,DPEnvelope * dpenv) 
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
      warn("Could not calculate packaln small for CloneWise due to wrong type of matrix");   
      return NULL;   
      }  


    out = PackAln_alloc_std();   


    start_reporting("Find start end points: ");  
    dc_optimised_start_end_calc_CloneWise(mat,dpenv);    
    score = start_end_find_end_CloneWise(mat,&endj); 
    out->score = score;  
    stopstate = END;
    
    /* Special to specials: have to eat up in strip and then drop back to full_dc for intervening bits */ 
    log_full_error(REPORT,0,"End at %d Score %d",endj,score);    
    stop_reporting();    
    for(;;)  { /*while there are more special bits to recover*/ 
      start_reporting("Special cell aln end   %d:",endj);    
      if( read_special_strip_CloneWise(mat,0,endj,stopstate,&endj,&startstate,out) == FALSE )    {  
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
      starti = CloneWise_DC_SHADOW_SPECIAL_SP(mat,0,endj,temp,0);    
      startj = CloneWise_DC_SHADOW_SPECIAL_SP(mat,0,endj,temp,1);    
      startstate = CloneWise_DC_SHADOW_SPECIAL_SP(mat,0,endj,temp,2);    
      stopi = CloneWise_DC_SHADOW_SPECIAL_SP(mat,0,endj,temp,3); 
      stopj = CloneWise_DC_SHADOW_SPECIAL_SP(mat,0,endj,temp,4); 
      stopstate = CloneWise_DC_SHADOW_SPECIAL_SP(mat,0,endj,temp,5); 


      /* Get out the score of this block. V. important! */ 
      temp = CloneWise_DC_SHADOW_SPECIAL_SP(mat,0,endj,temp,6);  
      totalj = stopj - startj;   
      donej  = 0;    
      start_reporting("Main matrix  aln [%d,%d]:",startj,stopj);     
      if(full_dc_CloneWise(mat,starti,startj,startstate,stopi,stopj,stopstate,out,&donej,totalj,dpenv) == FALSE) {  
        warn("In the alignment CloneWise [%d,%d][%d,%d], got a problem. Please report bug ... giving you back a partial alignment",starti,startj,stopi,stopj);   
        return out;  
        }  


      /* now have to figure out which special we came from... yikes */ 
      max_matrix_to_special_CloneWise(mat,starti,startj,startstate,temp,&stopi,&stopj,&stopstate,&temp,NULL);    
      if( stopi == CloneWise_READ_OFF_ERROR) {  
        warn("In CloneWise read off ending at %d ... got a bad matrix to special read off... returning partial alignment",startj);   
        invert_PackAln(out); 
        recalculate_PackAln_CloneWise(out,mat);  
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
    recalculate_PackAln_CloneWise(out,mat);  
    return out;  


}    


/* Function:  AlnRangeSet_calculate_Small_CloneWise(mat)
 *
 * Descrip:    This function calculates an alignment for CloneWise structure in linear space
 *             If you want the full alignment, use /PackAln_calculate_Small_CloneWise 
 *             If you have already got the full alignment, but want the range set, use /AlnRangeSet_from_PackAln_CloneWise
 *             If you have got the small matrix but not the alignment, use /AlnRangeSet_from_CloneWise 
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [CloneWise *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * AlnRangeSet_calculate_Small_CloneWise(CloneWise * mat) 
{
    AlnRangeSet * out;   


    start_reporting("Find start end points: ");  
    dc_optimised_start_end_calc_CloneWise(mat,NULL); 
    log_full_error(REPORT,0,"Calculated");   


    out = AlnRangeSet_from_CloneWise(mat);   
    return out;  
}    


/* Function:  AlnRangeSet_from_CloneWise(mat)
 *
 * Descrip:    This function reads off a start/end structure
 *             for CloneWise structure in linear space
 *             If you want the full alignment use
 *             /PackAln_calculate_Small_CloneWise 
 *             If you have not calculated the matrix use
 *             /AlnRange_calculate_Small_CloneWise
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [CloneWise *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * AlnRangeSet_from_CloneWise(CloneWise * mat) 
{
    AlnRangeSet * out;   
    AlnRange * temp; 
    int jpos;    
    int state;   


    if( mat->basematrix->type != BASEMATRIX_TYPE_SHADOW) {  
      warn("Bad error! - non shadow matrix type in AlnRangeSet_from_CloneWise"); 
      return NULL;   
      }  


    out = AlnRangeSet_alloc_std();   
    /* Find the end position */ 
    out->score = start_end_find_end_CloneWise(mat,&jpos);    
    state = END; 


    while( (temp = AlnRange_build_CloneWise(mat,jpos,state,&jpos,&state)) != NULL)   
      add_AlnRangeSet(out,temp); 
    return out;  
}    


/* Function:  AlnRange_build_CloneWise(mat,stopj,stopspecstate,startj,startspecstate)
 *
 * Descrip:    This function calculates a single start/end set in linear space
 *             Really a sub-routine for /AlnRangeSet_from_PackAln_CloneWise
 *
 *
 * Arg:                   mat [UNKN ] Undocumented argument [CloneWise *]
 * Arg:                 stopj [UNKN ] Undocumented argument [int]
 * Arg:         stopspecstate [UNKN ] Undocumented argument [int]
 * Arg:                startj [UNKN ] Undocumented argument [int *]
 * Arg:        startspecstate [UNKN ] Undocumented argument [int *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRange *]
 *
 */
AlnRange * AlnRange_build_CloneWise(CloneWise * mat,int stopj,int stopspecstate,int * startj,int * startspecstate) 
{
    AlnRange * out;  
    int jpos;    
    int state;   


    if( mat->basematrix->type != BASEMATRIX_TYPE_SHADOW) {  
      warn("Bad error! - non shadow matrix type in AlnRangeSet_from_CloneWise"); 
      return NULL;   
      }  


    /* Assumme that we have specials (we should!). Read back along the specials till we have the finish point */ 
    if( read_special_strip_CloneWise(mat,0,stopj,stopspecstate,&jpos,&state,NULL) == FALSE)  {  
      warn("In AlnRanger_build_CloneWise alignment ending at %d, unable to read back specials. Will (evenutally) return a partial range set... BEWARE!",stopj);  
      return NULL;   
      }  
    if( state == START || jpos <= 0) 
      return NULL;   


    out = AlnRange_alloc();  


    out->starti = CloneWise_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,0);    
    out->startj = CloneWise_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,1);    
    out->startstate = CloneWise_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,2);    
    out->stopi = CloneWise_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,3); 
    out->stopj = CloneWise_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,4); 
    out->stopstate = CloneWise_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,5); 
    out->startscore = CloneWise_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,6);    
    out->stopscore = CloneWise_DC_SHADOW_SPECIAL(mat,0,jpos,state);  


    /* Now, we have to figure out where this state came from in the specials */ 
    max_matrix_to_special_CloneWise(mat,out->starti,out->startj,out->startstate,out->startscore,&jpos,startj,startspecstate,&state,NULL);    
    if( jpos == CloneWise_READ_OFF_ERROR)    {  
      warn("In AlnRange_build_CloneWise alignment ending at %d, with aln range between %d-%d in j, unable to find source special, returning this range, but this could get tricky!",stopj,out->startj,out->stopj);   
      return out;    
      }  


    /* Put in the correct score for startstate, from the special */ 
    out->startscore = CloneWise_DC_SHADOW_SPECIAL(mat,0,*startj,*startspecstate);    
    /* The correct j coords have been put into startj, startspecstate... so just return out */ 
    return out;  
}    


/* Function:  read_hidden_CloneWise(mat,starti,startj,startstate,stopi,stopj,stopstate,out)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [CloneWise *]
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
boolean read_hidden_CloneWise(CloneWise * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out) 
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


      max_hidden_CloneWise(mat,startj,i,j,state,isspecial,&i,&j,&state,&isspecial,&cellscore);   


      if( i == CloneWise_READ_OFF_ERROR) {  
        warn("In CloneWise hidden read off, between %d:%d,%d:%d - at got bad read off. Problem!",starti,startj,stopi,stopj); 
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
        warn("In CloneWise hidden read off, between %d:%d,%d:%d - hit start cell, but not in start state. Can't be good!.",starti,startj,stopi,stopj);   
        return FALSE;    
        }  
      }  
    warn("In CloneWise hidden read off, between %d:%d,%d:%d - gone past start cell (now in %d,%d,%d), can't be good news!.",starti,startj,stopi,stopj,i,j,state);    
    return FALSE;    
}    


/* Function:  max_hidden_CloneWise(mat,hiddenj,i,j,state,isspecial,reti,retj,retstate,retspecial,cellscore)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [CloneWise *]
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
int max_hidden_CloneWise(CloneWise * mat,int hiddenj,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore) 
{
    register int temp;   
    register int cscore; 


    *reti = (*retj) = (*retstate) = CloneWise_READ_OFF_ERROR;    


    if( i < 0 || j < 0 || i > mat->q->length || j > mat->t->length)  {  
      warn("In CloneWise matrix special read off - out of bounds on matrix [i,j is %d,%d state %d in standard matrix]",i,j,state);   
      return -1; 
      }  


    /* Then you have to select the correct switch statement to figure out the readoff      */ 
    /* Somewhat odd - reverse the order of calculation and return as soon as it is correct */ 
    cscore = CloneWise_HIDDEN_MATRIX(mat,i,j,state); 
    switch(state)    { /*Switch state */ 
      case MATCH :   
        /* Not allowing special sources.. skipping START */ 
        temp = cscore - (0) -  (mat->match->matrix[i][j]);   
        if( temp == CloneWise_HIDDEN_MATRIX(mat,i - 1,j - 1,SKIP_TARGET) )   {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = SKIP_TARGET;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - CloneWise_HIDDEN_MATRIX(mat,i-1,j-1,SKIP_TARGET);  
            }  
          return CloneWise_HIDDEN_MATRIX(mat,i - 1,j - 1,SKIP_TARGET);   
          }  
        temp = cscore - (0) -  (mat->match->matrix[i][j]);   
        if( temp == CloneWise_HIDDEN_MATRIX(mat,i - 1,j - 1,SKIP_QUERY) )    {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = SKIP_QUERY;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - CloneWise_HIDDEN_MATRIX(mat,i-1,j-1,SKIP_QUERY);   
            }  
          return CloneWise_HIDDEN_MATRIX(mat,i - 1,j - 1,SKIP_QUERY);    
          }  
        temp = cscore - (0) -  (mat->match->matrix[i][j]);   
        if( temp == CloneWise_HIDDEN_MATRIX(mat,i - 1,j - 0,MATCH) ) {  
          *reti = i - 1; 
          *retj = j - 0; 
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - CloneWise_HIDDEN_MATRIX(mat,i-1,j-0,MATCH);    
            }  
          return CloneWise_HIDDEN_MATRIX(mat,i - 1,j - 0,MATCH);     
          }  
        temp = cscore - (0) -  (mat->match->matrix[i][j]);   
        if( temp == CloneWise_HIDDEN_MATRIX(mat,i - 0,j - 1,MATCH) ) {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - CloneWise_HIDDEN_MATRIX(mat,i-0,j-1,MATCH);    
            }  
          return CloneWise_HIDDEN_MATRIX(mat,i - 0,j - 1,MATCH);     
          }  
        temp = cscore - ((mat->match->matrix[i][j]+1)) -  (mat->match->matrix[i][j]);    
        if( temp == CloneWise_HIDDEN_MATRIX(mat,i - 1,j - 1,MATCH) ) {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - CloneWise_HIDDEN_MATRIX(mat,i-1,j-1,MATCH);    
            }  
          return CloneWise_HIDDEN_MATRIX(mat,i - 1,j - 1,MATCH);     
          }  
        warn("Major problem (!) - in CloneWise read off, position %d,%d state %d no source found!",i,j,state);   
        return (-1); 
      case SKIP_QUERY :  
        /* Not allowing special sources.. skipping START */ 
        temp = cscore - (0) -  (mat->match->skip_iset[i]);   
        if( temp == CloneWise_HIDDEN_MATRIX(mat,i - 1,j - 0,SKIP_QUERY) )    {  
          *reti = i - 1; 
          *retj = j - 0; 
          *retstate = SKIP_QUERY;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - CloneWise_HIDDEN_MATRIX(mat,i-1,j-0,SKIP_QUERY);   
            }  
          return CloneWise_HIDDEN_MATRIX(mat,i - 1,j - 0,SKIP_QUERY);    
          }  
        temp = cscore - (mat->query_skip_start) -  (mat->match->skip_iset[i]);   
        if( temp == CloneWise_HIDDEN_MATRIX(mat,i - 1,j - 0,SKIP_TARGET) )   {  
          *reti = i - 1; 
          *retj = j - 0; 
          *retstate = SKIP_TARGET;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - CloneWise_HIDDEN_MATRIX(mat,i-1,j-0,SKIP_TARGET);  
            }  
          return CloneWise_HIDDEN_MATRIX(mat,i - 1,j - 0,SKIP_TARGET);   
          }  
        temp = cscore - (mat->query_skip_start) -  (mat->match->skip_iset[i]);   
        if( temp == CloneWise_HIDDEN_MATRIX(mat,i - 1,j - 0,MATCH) ) {  
          *reti = i - 1; 
          *retj = j - 0; 
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - CloneWise_HIDDEN_MATRIX(mat,i-1,j-0,MATCH);    
            }  
          return CloneWise_HIDDEN_MATRIX(mat,i - 1,j - 0,MATCH);     
          }  
        warn("Major problem (!) - in CloneWise read off, position %d,%d state %d no source found!",i,j,state);   
        return (-1); 
      case SKIP_TARGET :     
        /* Not allowing special sources.. skipping START */ 
        temp = cscore - (0) -  (mat->match->skip_jset[j]);   
        if( temp == CloneWise_HIDDEN_MATRIX(mat,i - 0,j - 1,SKIP_TARGET) )   {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = SKIP_TARGET;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - CloneWise_HIDDEN_MATRIX(mat,i-0,j-1,SKIP_TARGET);  
            }  
          return CloneWise_HIDDEN_MATRIX(mat,i - 0,j - 1,SKIP_TARGET);   
          }  
        temp = cscore - (mat->target_skip_start) -  (mat->match->skip_jset[j]);  
        if( temp == CloneWise_HIDDEN_MATRIX(mat,i - 0,j - 1,SKIP_QUERY) )    {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = SKIP_QUERY;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - CloneWise_HIDDEN_MATRIX(mat,i-0,j-1,SKIP_QUERY);   
            }  
          return CloneWise_HIDDEN_MATRIX(mat,i - 0,j - 1,SKIP_QUERY);    
          }  
        temp = cscore - (mat->target_skip_start) -  (mat->match->skip_jset[j]);  
        if( temp == CloneWise_HIDDEN_MATRIX(mat,i - 0,j - 1,MATCH) ) {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - CloneWise_HIDDEN_MATRIX(mat,i-0,j-1,MATCH);    
            }  
          return CloneWise_HIDDEN_MATRIX(mat,i - 0,j - 1,MATCH);     
          }  
        warn("Major problem (!) - in CloneWise read off, position %d,%d state %d no source found!",i,j,state);   
        return (-1); 
      default:   
        warn("Major problem (!) - in CloneWise read off, position %d,%d state %d no source found!",i,j,state);   
        return (-1); 
      } /* end of Switch state  */ 
}    


/* Function:  read_special_strip_CloneWise(mat,stopi,stopj,stopstate,startj,startstate,out)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [CloneWise *]
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
boolean read_special_strip_CloneWise(CloneWise * mat,int stopi,int stopj,int stopstate,int * startj,int * startstate,PackAln * out) 
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
    while( j > CloneWise_DC_SHADOW_SPECIAL_SP(mat,i,j,state,4) && state != START)    { /*while more specials to eat up*/ 
      /* Put away current state, if we should */ 
      if(out != NULL)    {  
        pau = PackAlnUnit_alloc();  /* Should deal with memory overflow */ 
        pau->i = i;  
        pau->j = j;  
        pau->state =  state + 3; 
        add_PackAln(out,pau);    
        }  


      max_special_strip_CloneWise(mat,i,j,state,isspecial,&i,&j,&state,&isspecial,&cellscore);   
      if( i == CloneWise_READ_OFF_ERROR) {  
        warn("In special strip read CloneWise, got a bad read off error. Sorry!");   
        return FALSE;    
        }  
      } /* end of while more specials to eat up */ 


    /* check to see we have not gone too far! */ 
    if( state != START && j < CloneWise_DC_SHADOW_SPECIAL_SP(mat,i,j,state,4))   {  
      warn("In special strip read CloneWise, at special [%d] state [%d] overshot!",j,state); 
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


/* Function:  max_special_strip_CloneWise(mat,i,j,state,isspecial,reti,retj,retstate,retspecial,cellscore)
 *
 * Descrip:    A pretty intense internal function. Deals with read-off only in specials
 *
 *
 * Arg:               mat [UNKN ] Undocumented argument [CloneWise *]
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
int max_special_strip_CloneWise(CloneWise * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore) 
{
    int temp;    
    int cscore;  


    *reti = (*retj) = (*retstate) = CloneWise_READ_OFF_ERROR;    
    if( isspecial == FALSE ) {  
      warn("In special strip max function for CloneWise, got a non special start point. Problem! (bad!)");   
      return (-1);   
      }  


    if( j < 0 || j > mat->t->length) {  
      warn("In CloneWise matrix special read off - out of bounds on matrix [j is %d in special]",j); 
      return -1; 
      }  


    cscore = CloneWise_DC_SHADOW_SPECIAL(mat,i,j,state); 
    switch(state)    { /*switch on special states*/ 
      case START :   
      case TRUSTED_SPECIAL :     
        /* Source SKIP_TARGET is not a special */ 
        /* Source SKIP_QUERY is not a special */ 
        /* Source MATCH is not a special */ 
        /* source TRUSTED_SPECIAL is a special */ 
        temp = cscore - ((mat->match->skip_jset[j]-1)) - (0);    
        if( temp == CloneWise_DC_SHADOW_SPECIAL(mat,i - 0,j - 1,TRUSTED_SPECIAL) )   {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = TRUSTED_SPECIAL;   
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - CloneWise_DC_SHADOW_SPECIAL(mat,i-0,j-1,TRUSTED_SPECIAL);  
            }  
          return CloneWise_DC_SHADOW_MATRIX(mat,i - 0,j - 1,TRUSTED_SPECIAL) ;   
          }  
        /* source START is a special */ 
        temp = cscore - (0) - (0);   
        if( temp == CloneWise_DC_SHADOW_SPECIAL(mat,i - 0,j - 1,START) ) {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = START; 
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - CloneWise_DC_SHADOW_SPECIAL(mat,i-0,j-1,START);    
            }  
          return CloneWise_DC_SHADOW_MATRIX(mat,i - 0,j - 1,START) ;     
          }  
      case END :     
        /* source TRUSTED_SPECIAL is a special */ 
        temp = cscore - (0) - (0);   
        if( temp == CloneWise_DC_SHADOW_SPECIAL(mat,i - 0,j - 1,TRUSTED_SPECIAL) )   {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = TRUSTED_SPECIAL;   
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - CloneWise_DC_SHADOW_SPECIAL(mat,i-0,j-1,TRUSTED_SPECIAL);  
            }  
          return CloneWise_DC_SHADOW_MATRIX(mat,i - 0,j - 1,TRUSTED_SPECIAL) ;   
          }  
      default:   
        warn("Major problem (!) - in CloneWise special strip read off, position %d,%d state %d no source found  dropped into default on source switch!",i,j,state);  
        return (-1); 
      } /* end of switch on special states */ 
}    


/* Function:  max_matrix_to_special_CloneWise(mat,i,j,state,cscore,reti,retj,retstate,retspecial,cellscore)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [CloneWise *]
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
int max_matrix_to_special_CloneWise(CloneWise * mat,int i,int j,int state,int cscore,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore) 
{
    int temp;    
    *reti = (*retj) = (*retstate) = CloneWise_READ_OFF_ERROR;    


    if( j < 0 || j > mat->lenj)  {  
      warn("In CloneWise matrix to special read off - out of bounds on matrix [j is %d in special]",j);  
      return -1; 
      }  


    switch(state)    { /*Switch state */ 
      case MATCH :   
        temp = cscore - (0) -  (mat->match->matrix[i][j]);   
        if( temp == CloneWise_DC_SHADOW_SPECIAL(mat,i - 1,j - 1,START) ) {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = START; 
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - CloneWise_DC_SHADOW_SPECIAL(mat,i-1,j-1,START);    
            }  
          return CloneWise_DC_SHADOW_MATRIX(mat,i - 1,j - 1,START) ;     
          }  
        /* Source SKIP_TARGET is not a special, should not get here! */ 
        /* Source SKIP_QUERY is not a special, should not get here! */ 
        /* Source MATCH is not a special, should not get here! */ 
        /* Source MATCH is not a special, should not get here! */ 
        /* Source MATCH is not a special, should not get here! */ 
        warn("Major problem (!) - in CloneWise matrix to special read off, position %d,%d state %d no source found!",i,j,state); 
        return (-1); 
      case SKIP_QUERY :  
        temp = cscore - (mat->query_skip_start) -  (mat->match->skip_iset[i]);   
        if( temp == CloneWise_DC_SHADOW_SPECIAL(mat,i - 1,j - 0,START) ) {  
          *reti = i - 1; 
          *retj = j - 0; 
          *retstate = START; 
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - CloneWise_DC_SHADOW_SPECIAL(mat,i-1,j-0,START);    
            }  
          return CloneWise_DC_SHADOW_MATRIX(mat,i - 1,j - 0,START) ;     
          }  
        /* Source SKIP_QUERY is not a special, should not get here! */ 
        /* Source SKIP_TARGET is not a special, should not get here! */ 
        /* Source MATCH is not a special, should not get here! */ 
        warn("Major problem (!) - in CloneWise matrix to special read off, position %d,%d state %d no source found!",i,j,state); 
        return (-1); 
      case SKIP_TARGET :     
        temp = cscore - (mat->target_skip_start) -  (mat->match->skip_jset[j]);  
        if( temp == CloneWise_DC_SHADOW_SPECIAL(mat,i - 0,j - 1,START) ) {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = START; 
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - CloneWise_DC_SHADOW_SPECIAL(mat,i-0,j-1,START);    
            }  
          return CloneWise_DC_SHADOW_MATRIX(mat,i - 0,j - 1,START) ;     
          }  
        /* Source SKIP_TARGET is not a special, should not get here! */ 
        /* Source SKIP_QUERY is not a special, should not get here! */ 
        /* Source MATCH is not a special, should not get here! */ 
        warn("Major problem (!) - in CloneWise matrix to special read off, position %d,%d state %d no source found!",i,j,state); 
        return (-1); 
      default:   
        warn("Major problem (!) - in CloneWise read off, position %d,%d state %d no source found!",i,j,state);   
        return (-1); 
      } /* end of Switch state  */ 


}    


/* Function:  calculate_hidden_CloneWise(mat,starti,startj,startstate,stopi,stopj,stopstate,dpenv)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [CloneWise *]
 * Arg:            starti [UNKN ] Undocumented argument [int]
 * Arg:            startj [UNKN ] Undocumented argument [int]
 * Arg:        startstate [UNKN ] Undocumented argument [int]
 * Arg:             stopi [UNKN ] Undocumented argument [int]
 * Arg:             stopj [UNKN ] Undocumented argument [int]
 * Arg:         stopstate [UNKN ] Undocumented argument [int]
 * Arg:             dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 */
void calculate_hidden_CloneWise(CloneWise * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv) 
{
    register int i;  
    register int j;  
    register int score;  
    register int temp;   
    register int hiddenj;    


    hiddenj = startj;    


    init_hidden_CloneWise(mat,starti,startj,stopi,stopj);    


    CloneWise_HIDDEN_MATRIX(mat,starti,startj,startstate) = 0;   


    for(j=startj;j<=stopj;j++)   {  
      for(i=starti;i<=stopi;i++) {  
        /* Should *not* do very first cell as this is the one set to zero in one state! */ 
        if( i == starti && j == startj ) 
          continue;  
        if( dpenv != NULL && is_in_DPEnvelope(dpenv,i,j) == FALSE )  { /*Is not in envelope*/ 
          CloneWise_HIDDEN_MATRIX(mat,i,j,MATCH) = NEGI;     
          CloneWise_HIDDEN_MATRIX(mat,i,j,SKIP_QUERY) = NEGI;    
          CloneWise_HIDDEN_MATRIX(mat,i,j,SKIP_TARGET) = NEGI;   
          continue;  
          } /* end of Is not in envelope */ 


        /* For state MATCH */ 
        /* setting first movement to score */ 
        score = CloneWise_HIDDEN_MATRIX(mat,i-1,j-1,MATCH) + (mat->match->matrix[i][j]+1);   
        /* From state MATCH to state MATCH */ 
        temp = CloneWise_HIDDEN_MATRIX(mat,i-0,j-1,MATCH) + 0;   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = CloneWise_HIDDEN_MATRIX(mat,i-1,j-0,MATCH) + 0;   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SKIP_QUERY to state MATCH */ 
        temp = CloneWise_HIDDEN_MATRIX(mat,i-1,j-1,SKIP_QUERY) + 0;  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SKIP_TARGET to state MATCH */ 
        temp = CloneWise_HIDDEN_MATRIX(mat,i-1,j-1,SKIP_TARGET) + 0;     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for MATCH */ 
        /* Add any movement independant score and put away */ 
         score += mat->match->matrix[i][j];  
         CloneWise_HIDDEN_MATRIX(mat,i,j,MATCH) = score; 
        /* Finished calculating state MATCH */ 


        /* For state SKIP_QUERY */ 
        /* setting first movement to score */ 
        score = CloneWise_HIDDEN_MATRIX(mat,i-1,j-0,MATCH) + mat->query_skip_start;  
        /* From state SKIP_TARGET to state SKIP_QUERY */ 
        temp = CloneWise_HIDDEN_MATRIX(mat,i-1,j-0,SKIP_TARGET) + mat->query_skip_start;     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SKIP_QUERY to state SKIP_QUERY */ 
        temp = CloneWise_HIDDEN_MATRIX(mat,i-1,j-0,SKIP_QUERY) + 0;  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for SKIP_QUERY */ 
        /* Add any movement independant score and put away */ 
         score += mat->match->skip_iset[i];  
         CloneWise_HIDDEN_MATRIX(mat,i,j,SKIP_QUERY) = score;    
        /* Finished calculating state SKIP_QUERY */ 


        /* For state SKIP_TARGET */ 
        /* setting first movement to score */ 
        score = CloneWise_HIDDEN_MATRIX(mat,i-0,j-1,MATCH) + mat->target_skip_start;     
        /* From state SKIP_QUERY to state SKIP_TARGET */ 
        temp = CloneWise_HIDDEN_MATRIX(mat,i-0,j-1,SKIP_QUERY) + mat->target_skip_start;     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SKIP_TARGET to state SKIP_TARGET */ 
        temp = CloneWise_HIDDEN_MATRIX(mat,i-0,j-1,SKIP_TARGET) + 0;     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for SKIP_TARGET */ 
        /* Add any movement independant score and put away */ 
         score += mat->match->skip_jset[j];  
         CloneWise_HIDDEN_MATRIX(mat,i,j,SKIP_TARGET) = score;   
        /* Finished calculating state SKIP_TARGET */ 
        }  
      }  


    return;  
}    


/* Function:  init_hidden_CloneWise(mat,starti,startj,stopi,stopj)
 *
 * Descrip: No Description
 *
 * Arg:           mat [UNKN ] Undocumented argument [CloneWise *]
 * Arg:        starti [UNKN ] Undocumented argument [int]
 * Arg:        startj [UNKN ] Undocumented argument [int]
 * Arg:         stopi [UNKN ] Undocumented argument [int]
 * Arg:         stopj [UNKN ] Undocumented argument [int]
 *
 */
void init_hidden_CloneWise(CloneWise * mat,int starti,int startj,int stopi,int stopj) 
{
    register int i;  
    register int j;  
    register int hiddenj;    


    hiddenj = startj;    
    for(j=(startj-1);j<=stopj;j++)   {  
      for(i=(starti-1);i<=stopi;i++) {  
        CloneWise_HIDDEN_MATRIX(mat,i,j,MATCH) = NEGI;
  
        CloneWise_HIDDEN_MATRIX(mat,i,j,SKIP_QUERY) = NEGI;
 
        CloneWise_HIDDEN_MATRIX(mat,i,j,SKIP_TARGET) = NEGI;
    
        }  
      }  


    return;  
}    


/* Function:  full_dc_CloneWise(mat,starti,startj,startstate,stopi,stopj,stopstate,out,donej,totalj,dpenv)
 *
 * Descrip:    The main divide-and-conquor routine. Basically, call /PackAln_calculate_small_CloneWise
 *             Not this function, which is pretty hard core. 
 *             Function is given start/end points (in main matrix) for alignment
 *             It does some checks, decides whether start/end in j is small enough for explicit calc
 *               - if yes, calculates it, reads off into PackAln (out), adds the j distance to donej and returns TRUE
 *               - if no,  uses /do_dc_single_pass_CloneWise to get mid-point
 *                          saves midpoint, and calls itself to do right portion then left portion
 *             right then left ensures PackAln is added the 'right' way, ie, back-to-front
 *             returns FALSE on any error, with a warning
 *
 *
 * Arg:               mat [UNKN ] Matrix with small memory implementation [CloneWise *]
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
boolean full_dc_CloneWise(CloneWise * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out,int * donej,int totalj,DPEnvelope * dpenv) 
{
    int lstarti; 
    int lstartj; 
    int lstate;  


    if( mat->basematrix->type != BASEMATRIX_TYPE_SHADOW) {  
      warn("*Very* bad error! - non shadow matrix type in full_dc_CloneWise");   
      return FALSE;  
      }  


    if( starti == -1 || startj == -1 || startstate == -1 || stopi == -1 || stopstate == -1)  {  
      warn("In full dc program, passed bad indices, indices passed were %d:%d[%d] to %d:%d[%d]\n",starti,startj,startstate,stopi,stopj,stopstate);   
      return FALSE;  
      }  


    if( stopj - startj < 5)  {  
      log_full_error(REPORT,0,"[%d,%d][%d,%d] Explicit read off",starti,startj,stopi,stopj);/* Build hidden explicit matrix */ 
      calculate_hidden_CloneWise(mat,starti,startj,startstate,stopi,stopj,stopstate,dpenv);  
      *donej += (stopj - startj);   /* Now read it off into out */ 
      if( read_hidden_CloneWise(mat,starti,startj,startstate,stopi,stopj,stopstate,out) == FALSE)    {  
        warn("In full dc, at %d:%d,%d:%d got a bad hidden explicit read off... ",starti,startj,stopi,stopj); 
        return FALSE;    
        }  
      return TRUE;   
      }  


/* In actual divide and conquor */ 
    if( do_dc_single_pass_CloneWise(mat,starti,startj,startstate,stopi,stopj,stopstate,dpenv,(int)(*donej*100)/totalj) == FALSE) {  
      warn("In divide and conquor for CloneWise, at bound %d:%d to %d:%d, unable to calculate midpoint. Problem!",starti,startj,stopi,stopj);    
      return FALSE;  
      }  


/* Ok... now we have to call on each side of the matrix */ 
/* We have to retrieve left hand side positions, as they will be vapped by the time we call LHS */ 
    lstarti= CloneWise_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,0);     
    lstartj= CloneWise_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,1);     
    lstate = CloneWise_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,2);     


/* Call on right hand side: this lets us do the correct read off */ 
    if( full_dc_CloneWise(mat,CloneWise_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,3),CloneWise_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,4),CloneWise_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,5),stopi,stopj,stopstate,out,donej,totalj,dpenv) == FALSE)   {  
/* Warning already issued, simply chained back up to top */ 
      return FALSE;  
      }  
/* Call on left hand side */ 
    if( full_dc_CloneWise(mat,starti,startj,startstate,lstarti,lstartj,lstate,out,donej,totalj,dpenv) == FALSE)  {  
/* Warning already issued, simply chained back up to top */ 
      return FALSE;  
      }  


    return TRUE;     
}    


/* Function:  do_dc_single_pass_CloneWise(mat,starti,startj,startstate,stopi,stopj,stopstate,dpenv,perc_done)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [CloneWise *]
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
boolean do_dc_single_pass_CloneWise(CloneWise * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv,int perc_done) 
{
    int halfj;   
    halfj = startj + ((stopj - startj)/2);   


    init_dc_CloneWise(mat);  


    CloneWise_DC_SHADOW_MATRIX(mat,starti,startj,startstate) = 0;    
    run_up_dc_CloneWise(mat,starti,stopi,startj,halfj-1,dpenv,perc_done);    
    push_dc_at_merge_CloneWise(mat,starti,stopi,halfj,&halfj,dpenv);     
    follow_on_dc_CloneWise(mat,starti,stopi,halfj,stopj,dpenv,perc_done);    
    return TRUE; 
}    


/* Function:  push_dc_at_merge_CloneWise(mat,starti,stopi,startj,stopj,dpenv)
 *
 * Descrip: No Description
 *
 * Arg:           mat [UNKN ] Undocumented argument [CloneWise *]
 * Arg:        starti [UNKN ] Undocumented argument [int]
 * Arg:         stopi [UNKN ] Undocumented argument [int]
 * Arg:        startj [UNKN ] Undocumented argument [int]
 * Arg:         stopj [UNKN ] Undocumented argument [int *]
 * Arg:         dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 */
void push_dc_at_merge_CloneWise(CloneWise * mat,int starti,int stopi,int startj,int * stopj,DPEnvelope * dpenv) 
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
          CloneWise_DC_SHADOW_MATRIX(mat,i,j,MATCH) = NEGI;  
          CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,0) = (-100);   
          CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,1) = (-100);   
          CloneWise_DC_SHADOW_MATRIX(mat,i,j,SKIP_QUERY) = NEGI;     
          CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,SKIP_QUERY,0) = (-100);  
          CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,SKIP_QUERY,1) = (-100);  
          CloneWise_DC_SHADOW_MATRIX(mat,i,j,SKIP_TARGET) = NEGI;    
          CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,SKIP_TARGET,0) = (-100); 
          CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,SKIP_TARGET,1) = (-100); 
          continue;  
          } /* end of Is not in envelope */ 


        /* For state MATCH, pushing when j - offj <= mergej */ 
        score = CloneWise_DC_SHADOW_MATRIX(mat,i-1,j-1,MATCH) + (mat->match->matrix[i][j]+1);    
        if( j - 1 <= mergej) {  
          CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,0) = i-1;  
          CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,1) = j-1;  
          CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,2) = MATCH;    
          CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,3) = i;    
          CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,4) = j;    
          CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,5) = MATCH;    
          }  
        else {  
          for(k=0;k<7;k++)   
            CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = CloneWise_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,MATCH,k); 
          }  


        temp = CloneWise_DC_SHADOW_MATRIX(mat,i-0,j-1,MATCH) + 0;    
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,0) = i-0;    
            CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,1) = j-1;    
            CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,2) = MATCH;  
            CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,3) = i;  
            CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,4) = j;  
            CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,5) = MATCH;  
            }  
          else   {  
            for(k=0;k<7;k++) 
              CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = CloneWise_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,MATCH,k);   
            }  
          }  


        temp = CloneWise_DC_SHADOW_MATRIX(mat,i-1,j-0,MATCH) + 0;    
        if( temp > score)    {  
          score = temp;  


          if( j - 0 <= mergej)   {  
            CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,0) = i-1;    
            CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,1) = j-0;    
            CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,2) = MATCH;  
            CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,3) = i;  
            CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,4) = j;  
            CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,5) = MATCH;  
            }  
          else   {  
            for(k=0;k<7;k++) 
              CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = CloneWise_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 0,MATCH,k);   
            }  
          }  


        temp = CloneWise_DC_SHADOW_MATRIX(mat,i-1,j-1,SKIP_QUERY) + 0;   
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,0) = i-1;    
            CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,1) = j-1;    
            CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,2) = SKIP_QUERY; 
            CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,3) = i;  
            CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,4) = j;  
            CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,5) = MATCH;  
            }  
          else   {  
            for(k=0;k<7;k++) 
              CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = CloneWise_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,SKIP_QUERY,k);  
            }  
          }  


        temp = CloneWise_DC_SHADOW_MATRIX(mat,i-1,j-1,SKIP_TARGET) + 0;  
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,0) = i-1;    
            CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,1) = j-1;    
            CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,2) = SKIP_TARGET;    
            CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,3) = i;  
            CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,4) = j;  
            CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,5) = MATCH;  
            }  
          else   {  
            for(k=0;k<7;k++) 
              CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = CloneWise_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,SKIP_TARGET,k); 
            }  
          }  
        /* Add any movement independant score */ 
        score += mat->match->matrix[i][j];   
        CloneWise_DC_SHADOW_MATRIX(mat,i,j,MATCH) = score;   
        /* Finished with state MATCH */ 


        /* For state SKIP_QUERY, pushing when j - offj <= mergej */ 
        score = CloneWise_DC_SHADOW_MATRIX(mat,i-1,j-0,MATCH) + mat->query_skip_start;   
        if( j - 0 <= mergej) {  
          CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,SKIP_QUERY,0) = i-1; 
          CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,SKIP_QUERY,1) = j-0; 
          CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,SKIP_QUERY,2) = MATCH;   
          CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,SKIP_QUERY,3) = i;   
          CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,SKIP_QUERY,4) = j;   
          CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,SKIP_QUERY,5) = SKIP_QUERY;  
          }  
        else {  
          for(k=0;k<7;k++)   
            CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,SKIP_QUERY,k) = CloneWise_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 0,MATCH,k);    
          }  


        temp = CloneWise_DC_SHADOW_MATRIX(mat,i-1,j-0,SKIP_TARGET) + mat->query_skip_start;  
        if( temp > score)    {  
          score = temp;  


          if( j - 0 <= mergej)   {  
            CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,SKIP_QUERY,0) = i-1;   
            CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,SKIP_QUERY,1) = j-0;   
            CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,SKIP_QUERY,2) = SKIP_TARGET;   
            CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,SKIP_QUERY,3) = i; 
            CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,SKIP_QUERY,4) = j; 
            CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,SKIP_QUERY,5) = SKIP_QUERY;    
            }  
          else   {  
            for(k=0;k<7;k++) 
              CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,SKIP_QUERY,k) = CloneWise_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 0,SKIP_TARGET,k);    
            }  
          }  


        temp = CloneWise_DC_SHADOW_MATRIX(mat,i-1,j-0,SKIP_QUERY) + 0;   
        if( temp > score)    {  
          score = temp;  


          if( j - 0 <= mergej)   {  
            CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,SKIP_QUERY,0) = i-1;   
            CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,SKIP_QUERY,1) = j-0;   
            CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,SKIP_QUERY,2) = SKIP_QUERY;    
            CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,SKIP_QUERY,3) = i; 
            CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,SKIP_QUERY,4) = j; 
            CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,SKIP_QUERY,5) = SKIP_QUERY;    
            }  
          else   {  
            for(k=0;k<7;k++) 
              CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,SKIP_QUERY,k) = CloneWise_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 0,SKIP_QUERY,k); 
            }  
          }  
        /* Add any movement independant score */ 
        score += mat->match->skip_iset[i];   
        CloneWise_DC_SHADOW_MATRIX(mat,i,j,SKIP_QUERY) = score;  
        /* Finished with state SKIP_QUERY */ 


        /* For state SKIP_TARGET, pushing when j - offj <= mergej */ 
        score = CloneWise_DC_SHADOW_MATRIX(mat,i-0,j-1,MATCH) + mat->target_skip_start;  
        if( j - 1 <= mergej) {  
          CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,SKIP_TARGET,0) = i-0;    
          CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,SKIP_TARGET,1) = j-1;    
          CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,SKIP_TARGET,2) = MATCH;  
          CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,SKIP_TARGET,3) = i;  
          CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,SKIP_TARGET,4) = j;  
          CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,SKIP_TARGET,5) = SKIP_TARGET;    
          }  
        else {  
          for(k=0;k<7;k++)   
            CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,SKIP_TARGET,k) = CloneWise_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,MATCH,k);   
          }  


        temp = CloneWise_DC_SHADOW_MATRIX(mat,i-0,j-1,SKIP_QUERY) + mat->target_skip_start;  
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,SKIP_TARGET,0) = i-0;  
            CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,SKIP_TARGET,1) = j-1;  
            CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,SKIP_TARGET,2) = SKIP_QUERY;   
            CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,SKIP_TARGET,3) = i;    
            CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,SKIP_TARGET,4) = j;    
            CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,SKIP_TARGET,5) = SKIP_TARGET;  
            }  
          else   {  
            for(k=0;k<7;k++) 
              CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,SKIP_TARGET,k) = CloneWise_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,SKIP_QUERY,k);    
            }  
          }  


        temp = CloneWise_DC_SHADOW_MATRIX(mat,i-0,j-1,SKIP_TARGET) + 0;  
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,SKIP_TARGET,0) = i-0;  
            CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,SKIP_TARGET,1) = j-1;  
            CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,SKIP_TARGET,2) = SKIP_TARGET;  
            CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,SKIP_TARGET,3) = i;    
            CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,SKIP_TARGET,4) = j;    
            CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,SKIP_TARGET,5) = SKIP_TARGET;  
            }  
          else   {  
            for(k=0;k<7;k++) 
              CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,SKIP_TARGET,k) = CloneWise_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,SKIP_TARGET,k);   
            }  
          }  
        /* Add any movement independant score */ 
        score += mat->match->skip_jset[j];   
        CloneWise_DC_SHADOW_MATRIX(mat,i,j,SKIP_TARGET) = score;     
        /* Finished with state SKIP_TARGET */ 
        }  
      }  
    /* Put back j into * stop j so that calling function gets it correct */ 
    if( stopj == NULL)   
      warn("Bad news... NULL stopj pointer in push dc function. This means that calling function does not know how many cells I have done!");    
    else 
      *stopj = j;    


    return;  
}    


/* Function:  follow_on_dc_CloneWise(mat,starti,stopi,startj,stopj,dpenv,perc_done)
 *
 * Descrip: No Description
 *
 * Arg:              mat [UNKN ] Undocumented argument [CloneWise *]
 * Arg:           starti [UNKN ] Undocumented argument [int]
 * Arg:            stopi [UNKN ] Undocumented argument [int]
 * Arg:           startj [UNKN ] Undocumented argument [int]
 * Arg:            stopj [UNKN ] Undocumented argument [int]
 * Arg:            dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:        perc_done [UNKN ] Undocumented argument [int]
 *
 */
void follow_on_dc_CloneWise(CloneWise * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done) 
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
          CloneWise_DC_SHADOW_MATRIX(mat,i,j,MATCH) = NEGI;  
          CloneWise_DC_SHADOW_MATRIX(mat,i,j,SKIP_QUERY) = NEGI;     
          CloneWise_DC_SHADOW_MATRIX(mat,i,j,SKIP_TARGET) = NEGI;    
          continue;  
          } /* end of Is not in envelope */ 
        if( num % 1000 == 0 )    
          log_full_error(REPORT,0,"[%d%%%% done]After  mid-j %5d Cells done %d%%%%",perc_done,startj,(num*100)/total);   


        /* For state MATCH */ 
        /* setting first movement to score */ 
        score = CloneWise_DC_SHADOW_MATRIX(mat,i-1,j-1,MATCH) + (mat->match->matrix[i][j]+1);    
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = CloneWise_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,MATCH,k);   
        /* From state MATCH to state MATCH */ 
        temp = CloneWise_DC_SHADOW_MATRIX(mat,i-0,j-1,MATCH) + 0;    
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = CloneWise_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,MATCH,k); 
          }  
        /* From state MATCH to state MATCH */ 
        temp = CloneWise_DC_SHADOW_MATRIX(mat,i-1,j-0,MATCH) + 0;    
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = CloneWise_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 0,MATCH,k); 
          }  
        /* From state SKIP_QUERY to state MATCH */ 
        temp = CloneWise_DC_SHADOW_MATRIX(mat,i-1,j-1,SKIP_QUERY) + 0;   
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = CloneWise_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,SKIP_QUERY,k);    
          }  
        /* From state SKIP_TARGET to state MATCH */ 
        temp = CloneWise_DC_SHADOW_MATRIX(mat,i-1,j-1,SKIP_TARGET) + 0;  
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = CloneWise_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,SKIP_TARGET,k);   
          }  


        /* Ok - finished max calculation for MATCH */ 
        /* Add any movement independant score and put away */ 
         score += mat->match->matrix[i][j];  
         CloneWise_DC_SHADOW_MATRIX(mat,i,j,MATCH) = score;  
        for(k=0;k<7;k++) 
          CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = localshadow[k];   
        /* Now figure out if any specials need this score */ 
        /* Finished calculating state MATCH */ 


        /* For state SKIP_QUERY */ 
        /* setting first movement to score */ 
        score = CloneWise_DC_SHADOW_MATRIX(mat,i-1,j-0,MATCH) + mat->query_skip_start;   
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = CloneWise_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 0,MATCH,k);   
        /* From state SKIP_TARGET to state SKIP_QUERY */ 
        temp = CloneWise_DC_SHADOW_MATRIX(mat,i-1,j-0,SKIP_TARGET) + mat->query_skip_start;  
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = CloneWise_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 0,SKIP_TARGET,k);   
          }  
        /* From state SKIP_QUERY to state SKIP_QUERY */ 
        temp = CloneWise_DC_SHADOW_MATRIX(mat,i-1,j-0,SKIP_QUERY) + 0;   
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = CloneWise_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 0,SKIP_QUERY,k);    
          }  


        /* Ok - finished max calculation for SKIP_QUERY */ 
        /* Add any movement independant score and put away */ 
         score += mat->match->skip_iset[i];  
         CloneWise_DC_SHADOW_MATRIX(mat,i,j,SKIP_QUERY) = score; 
        for(k=0;k<7;k++) 
          CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,SKIP_QUERY,k) = localshadow[k];  
        /* Now figure out if any specials need this score */ 
        /* Finished calculating state SKIP_QUERY */ 


        /* For state SKIP_TARGET */ 
        /* setting first movement to score */ 
        score = CloneWise_DC_SHADOW_MATRIX(mat,i-0,j-1,MATCH) + mat->target_skip_start;  
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = CloneWise_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,MATCH,k);   
        /* From state SKIP_QUERY to state SKIP_TARGET */ 
        temp = CloneWise_DC_SHADOW_MATRIX(mat,i-0,j-1,SKIP_QUERY) + mat->target_skip_start;  
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = CloneWise_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,SKIP_QUERY,k);    
          }  
        /* From state SKIP_TARGET to state SKIP_TARGET */ 
        temp = CloneWise_DC_SHADOW_MATRIX(mat,i-0,j-1,SKIP_TARGET) + 0;  
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = CloneWise_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,SKIP_TARGET,k);   
          }  


        /* Ok - finished max calculation for SKIP_TARGET */ 
        /* Add any movement independant score and put away */ 
         score += mat->match->skip_jset[j];  
         CloneWise_DC_SHADOW_MATRIX(mat,i,j,SKIP_TARGET) = score;    
        for(k=0;k<7;k++) 
          CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,SKIP_TARGET,k) = localshadow[k]; 
        /* Now figure out if any specials need this score */ 
        /* Finished calculating state SKIP_TARGET */ 
        } /* end of this is strip */ 
      } /* end of for each valid j column */ 


/* Function:  run_up_dc_CloneWise(mat,starti,stopi,startj,stopj,dpenv,perc_done)
 *
 * Descrip: No Description
 *
 * Arg:              mat [UNKN ] Undocumented argument [CloneWise *]
 * Arg:           starti [UNKN ] Undocumented argument [int]
 * Arg:            stopi [UNKN ] Undocumented argument [int]
 * Arg:           startj [UNKN ] Undocumented argument [int]
 * Arg:            stopj [UNKN ] Undocumented argument [int]
 * Arg:            dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:        perc_done [UNKN ] Undocumented argument [int]
 *
 */
}    
void run_up_dc_CloneWise(CloneWise * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done) 
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
          CloneWise_DC_SHADOW_MATRIX(mat,i,j,MATCH) = NEGI;  
          CloneWise_DC_SHADOW_MATRIX(mat,i,j,SKIP_QUERY) = NEGI;     
          CloneWise_DC_SHADOW_MATRIX(mat,i,j,SKIP_TARGET) = NEGI;    
          continue;  
          } /* end of Is not in envelope */ 
        if( num % 1000 == 0 )    
          log_full_error(REPORT,0,"[%d%%%% done]Before mid-j %5d Cells done %d%%%%",perc_done,stopj,(num*100)/total);    


        /* For state MATCH */ 
        /* setting first movement to score */ 
        score = CloneWise_DC_SHADOW_MATRIX(mat,i-1,j-1,MATCH) + (mat->match->matrix[i][j]+1);    
        /* From state MATCH to state MATCH */ 
        temp = CloneWise_DC_SHADOW_MATRIX(mat,i-0,j-1,MATCH) + 0;    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = CloneWise_DC_SHADOW_MATRIX(mat,i-1,j-0,MATCH) + 0;    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SKIP_QUERY to state MATCH */ 
        temp = CloneWise_DC_SHADOW_MATRIX(mat,i-1,j-1,SKIP_QUERY) + 0;   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SKIP_TARGET to state MATCH */ 
        temp = CloneWise_DC_SHADOW_MATRIX(mat,i-1,j-1,SKIP_TARGET) + 0;  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for MATCH */ 
        /* Add any movement independant score and put away */ 
         score += mat->match->matrix[i][j];  
         CloneWise_DC_SHADOW_MATRIX(mat,i,j,MATCH) = score;  
        /* Finished calculating state MATCH */ 


        /* For state SKIP_QUERY */ 
        /* setting first movement to score */ 
        score = CloneWise_DC_SHADOW_MATRIX(mat,i-1,j-0,MATCH) + mat->query_skip_start;   
        /* From state SKIP_TARGET to state SKIP_QUERY */ 
        temp = CloneWise_DC_SHADOW_MATRIX(mat,i-1,j-0,SKIP_TARGET) + mat->query_skip_start;  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SKIP_QUERY to state SKIP_QUERY */ 
        temp = CloneWise_DC_SHADOW_MATRIX(mat,i-1,j-0,SKIP_QUERY) + 0;   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for SKIP_QUERY */ 
        /* Add any movement independant score and put away */ 
         score += mat->match->skip_iset[i];  
         CloneWise_DC_SHADOW_MATRIX(mat,i,j,SKIP_QUERY) = score; 
        /* Finished calculating state SKIP_QUERY */ 


        /* For state SKIP_TARGET */ 
        /* setting first movement to score */ 
        score = CloneWise_DC_SHADOW_MATRIX(mat,i-0,j-1,MATCH) + mat->target_skip_start;  
        /* From state SKIP_QUERY to state SKIP_TARGET */ 
        temp = CloneWise_DC_SHADOW_MATRIX(mat,i-0,j-1,SKIP_QUERY) + mat->target_skip_start;  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SKIP_TARGET to state SKIP_TARGET */ 
        temp = CloneWise_DC_SHADOW_MATRIX(mat,i-0,j-1,SKIP_TARGET) + 0;  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for SKIP_TARGET */ 
        /* Add any movement independant score and put away */ 
         score += mat->match->skip_jset[j];  
         CloneWise_DC_SHADOW_MATRIX(mat,i,j,SKIP_TARGET) = score;    
        /* Finished calculating state SKIP_TARGET */ 
        } /* end of this is strip */ 
      } /* end of for each valid j column */ 


/* Function:  init_dc_CloneWise(mat)
 *
 * Descrip: No Description
 *
 * Arg:        mat [UNKN ] Undocumented argument [CloneWise *]
 *
 */
}    
void init_dc_CloneWise(CloneWise * mat) 
{
    register int i;  
    register int j;  
    register int k;  


    for(j=0;j<3;j++) {  
      for(i=(-1);i<mat->q->length;i++)   {  
        CloneWise_DC_SHADOW_MATRIX(mat,i,j,MATCH) = NEGI;    
        CloneWise_DC_SHADOW_MATRIX(mat,i,j,SKIP_QUERY) = NEGI;   
        CloneWise_DC_SHADOW_MATRIX(mat,i,j,SKIP_TARGET) = NEGI;  
        for(k=0;k<7;k++) {  
          CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = (-1); 
          CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,SKIP_QUERY,k) = (-1);    
          CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,SKIP_TARGET,k) = (-1);   
          }  
        }  
      }  


    return;  
}    


/* Function:  start_end_find_end_CloneWise(mat,endj)
 *
 * Descrip:    First function used to find end of the best path in the special state !end
 *
 *
 * Arg:         mat [UNKN ] Matrix in small mode [CloneWise *]
 * Arg:        endj [WRITE] position of end in j (meaningless in i) [int *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int start_end_find_end_CloneWise(CloneWise * mat,int * endj) 
{
    register int j;  
    register int max;    
    register int maxj;   


    max = CloneWise_DC_SHADOW_SPECIAL(mat,0,mat->t->length-1,END);   
    maxj = mat->t->length-1;     
    for(j= mat->t->length-2 ;j >= 0 ;j--)    {  
      if( CloneWise_DC_SHADOW_SPECIAL(mat,0,j,END) > max )   {  
        max = CloneWise_DC_SHADOW_SPECIAL(mat,0,j,END);  
        maxj = j;    
        }  
      }  


    if( endj != NULL)    
      *endj = maxj;  


    return max;  
}    


/* Function:  dc_optimised_start_end_calc_CloneWise(*mat,dpenv)
 *
 * Descrip:    Calculates special strip, leaving start/end/score points in shadow matrix
 *             Works off specially laid out memory from steve searle
 *
 *
 * Arg:         *mat [UNKN ] Undocumented argument [CloneWise]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean dc_optimised_start_end_calc_CloneWise(CloneWise *mat,DPEnvelope * dpenv) 
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
    leni = mat->q->length;   
    lenj = mat->t->length;   
    total = leni * lenj; 


    score_pointers = (int *) calloc (1 * (leni + 1) * 3,sizeof(int));    
    shadow_pointers = (int *) calloc (1 * (leni + 1) * 3 * 8,sizeof(int));   


    for(j=0;j<lenj;j++)  { /*for each j strip*/ 
      for(i=0;i<leni;i++)    { /*for each i position in strip*/ 
        num++;   
        if( dpenv != NULL && is_in_DPEnvelope(dpenv,i,j) == FALSE )  { /*Is not in envelope*/ 
          CloneWise_DC_OPT_SHADOW_MATRIX(mat,i,j,MATCH) = NEGI;  
          CloneWise_DC_OPT_SHADOW_MATRIX(mat,i,j,SKIP_QUERY) = NEGI;     
          CloneWise_DC_OPT_SHADOW_MATRIX(mat,i,j,SKIP_TARGET) = NEGI;    
          continue;  
          } /* end of Is not in envelope */ 
        if( num%1000 == 0)   
          log_full_error(REPORT,0,"%6d Cells done [%2d%%%%]",num,num*100/total); 




        /* For state MATCH */ 
        /* setting first movement to score */ 
        score = CloneWise_DC_OPT_SHADOW_MATRIX(mat,i-1,j-1,MATCH) + (mat->match->matrix[i][j]+1) + (mat->match->matrix[i][j]);   
        /* assign local shadown pointer */ 
        localsp = &(CloneWise_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 1,MATCH,0)); 
        /* From state MATCH to state MATCH */ 
        temp = CloneWise_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,MATCH) + 0 +(mat->match->matrix[i][j]);    
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(CloneWise_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,MATCH,0));   
          }  
        /* From state MATCH to state MATCH */ 
        temp = CloneWise_DC_OPT_SHADOW_MATRIX(mat,i-1,j-0,MATCH) + 0 +(mat->match->matrix[i][j]);    
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(CloneWise_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 0,MATCH,0));   
          }  
        /* From state SKIP_QUERY to state MATCH */ 
        temp = CloneWise_DC_OPT_SHADOW_MATRIX(mat,i-1,j-1,SKIP_QUERY) + 0 +(mat->match->matrix[i][j]);   
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(CloneWise_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 1,SKIP_QUERY,0));  
          }  
        /* From state SKIP_TARGET to state MATCH */ 
        temp = CloneWise_DC_OPT_SHADOW_MATRIX(mat,i-1,j-1,SKIP_TARGET) + 0 +(mat->match->matrix[i][j]);  
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(CloneWise_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 1,SKIP_TARGET,0)); 
          }  
        /* From state START to state MATCH */ 
        temp = CloneWise_DC_OPT_SHADOW_SPECIAL(mat,i-1,j-1,START) + 0 + (mat->match->matrix[i][j]);  
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


        /* Ok - finished max calculation for MATCH */ 
        /* Add any movement independant score and put away */ 
        /* Actually, already done inside scores */ 
         CloneWise_DC_OPT_SHADOW_MATRIX(mat,i,j,MATCH) = score;  
        for(k=0;k<7;k++) 
          CloneWise_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = localsp[k];   
        /* Now figure out if any specials need this score */ 


        /* state MATCH is a source for special TRUSTED_SPECIAL */ 
        temp = score + (mat->target_special_s) + (0) ;   
        if( temp > CloneWise_DC_OPT_SHADOW_SPECIAL(mat,i,j,TRUSTED_SPECIAL) )    {  
          CloneWise_DC_OPT_SHADOW_SPECIAL(mat,i,j,TRUSTED_SPECIAL) = temp;   
          /* Have to push only bottem half of system here */ 
          for(k=0;k<3;k++)   
            CloneWise_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,TRUSTED_SPECIAL,k) = CloneWise_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,MATCH,k);  
          CloneWise_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,TRUSTED_SPECIAL,6) = CloneWise_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,MATCH,6);    
          CloneWise_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,TRUSTED_SPECIAL,3) = i; 
          CloneWise_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,TRUSTED_SPECIAL,4) = j; 
          CloneWise_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,TRUSTED_SPECIAL,5) = MATCH; 
          }  




        /* Finished calculating state MATCH */ 


        /* For state SKIP_QUERY */ 
        /* setting first movement to score */ 
        score = CloneWise_DC_OPT_SHADOW_MATRIX(mat,i-1,j-0,MATCH) + mat->query_skip_start + (mat->match->skip_iset[i]);  
        /* assign local shadown pointer */ 
        localsp = &(CloneWise_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 0,MATCH,0)); 
        /* From state SKIP_TARGET to state SKIP_QUERY */ 
        temp = CloneWise_DC_OPT_SHADOW_MATRIX(mat,i-1,j-0,SKIP_TARGET) + mat->query_skip_start +(mat->match->skip_iset[i]);  
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(CloneWise_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 0,SKIP_TARGET,0)); 
          }  
        /* From state SKIP_QUERY to state SKIP_QUERY */ 
        temp = CloneWise_DC_OPT_SHADOW_MATRIX(mat,i-1,j-0,SKIP_QUERY) + 0 +(mat->match->skip_iset[i]);   
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(CloneWise_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 0,SKIP_QUERY,0));  
          }  
        /* From state START to state SKIP_QUERY */ 
        temp = CloneWise_DC_OPT_SHADOW_SPECIAL(mat,i-1,j-0,START) + mat->query_skip_start + (mat->match->skip_iset[i]);  
        if( temp  > score )  {  
          score = temp;  
          /* This state [START] is a special for SKIP_QUERY... push top shadow pointers here */ 
          localshadow[0]= i; 
          localshadow[1]= j; 
          localshadow[2]= SKIP_QUERY;    
          localshadow[3]= (-1);  
          localshadow[4]= (-1);  
          localshadow[5]= (-1);  
          localshadow[6]= score; 
          localsp = localshadow; 
          }  


        /* Ok - finished max calculation for SKIP_QUERY */ 
        /* Add any movement independant score and put away */ 
        /* Actually, already done inside scores */ 
         CloneWise_DC_OPT_SHADOW_MATRIX(mat,i,j,SKIP_QUERY) = score; 
        for(k=0;k<7;k++) 
          CloneWise_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,SKIP_QUERY,k) = localsp[k];  
        /* Now figure out if any specials need this score */ 


        /* state SKIP_QUERY is a source for special TRUSTED_SPECIAL */ 
        temp = score + (mat->target_special_s) + (0) ;   
        if( temp > CloneWise_DC_OPT_SHADOW_SPECIAL(mat,i,j,TRUSTED_SPECIAL) )    {  
          CloneWise_DC_OPT_SHADOW_SPECIAL(mat,i,j,TRUSTED_SPECIAL) = temp;   
          /* Have to push only bottem half of system here */ 
          for(k=0;k<3;k++)   
            CloneWise_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,TRUSTED_SPECIAL,k) = CloneWise_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,SKIP_QUERY,k); 
          CloneWise_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,TRUSTED_SPECIAL,6) = CloneWise_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,SKIP_QUERY,6);   
          CloneWise_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,TRUSTED_SPECIAL,3) = i; 
          CloneWise_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,TRUSTED_SPECIAL,4) = j; 
          CloneWise_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,TRUSTED_SPECIAL,5) = SKIP_QUERY;    
          }  




        /* Finished calculating state SKIP_QUERY */ 


        /* For state SKIP_TARGET */ 
        /* setting first movement to score */ 
        score = CloneWise_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,MATCH) + mat->target_skip_start + (mat->match->skip_jset[j]);     
        /* assign local shadown pointer */ 
        localsp = &(CloneWise_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,MATCH,0)); 
        /* From state SKIP_QUERY to state SKIP_TARGET */ 
        temp = CloneWise_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,SKIP_QUERY) + mat->target_skip_start +(mat->match->skip_jset[j]);  
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(CloneWise_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,SKIP_QUERY,0));  
          }  
        /* From state SKIP_TARGET to state SKIP_TARGET */ 
        temp = CloneWise_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,SKIP_TARGET) + 0 +(mat->match->skip_jset[j]);  
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(CloneWise_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,SKIP_TARGET,0)); 
          }  
        /* From state START to state SKIP_TARGET */ 
        temp = CloneWise_DC_OPT_SHADOW_SPECIAL(mat,i-0,j-1,START) + mat->target_skip_start + (mat->match->skip_jset[j]);     
        if( temp  > score )  {  
          score = temp;  
          /* This state [START] is a special for SKIP_TARGET... push top shadow pointers here */ 
          localshadow[0]= i; 
          localshadow[1]= j; 
          localshadow[2]= SKIP_TARGET;   
          localshadow[3]= (-1);  
          localshadow[4]= (-1);  
          localshadow[5]= (-1);  
          localshadow[6]= score; 
          localsp = localshadow; 
          }  


        /* Ok - finished max calculation for SKIP_TARGET */ 
        /* Add any movement independant score and put away */ 
        /* Actually, already done inside scores */ 
         CloneWise_DC_OPT_SHADOW_MATRIX(mat,i,j,SKIP_TARGET) = score;    
        for(k=0;k<7;k++) 
          CloneWise_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,SKIP_TARGET,k) = localsp[k]; 
        /* Now figure out if any specials need this score */ 


        /* state SKIP_TARGET is a source for special TRUSTED_SPECIAL */ 
        temp = score + (mat->target_special_s) + (0) ;   
        if( temp > CloneWise_DC_OPT_SHADOW_SPECIAL(mat,i,j,TRUSTED_SPECIAL) )    {  
          CloneWise_DC_OPT_SHADOW_SPECIAL(mat,i,j,TRUSTED_SPECIAL) = temp;   
          /* Have to push only bottem half of system here */ 
          for(k=0;k<3;k++)   
            CloneWise_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,TRUSTED_SPECIAL,k) = CloneWise_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,SKIP_TARGET,k);    
          CloneWise_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,TRUSTED_SPECIAL,6) = CloneWise_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,SKIP_TARGET,6);  
          CloneWise_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,TRUSTED_SPECIAL,3) = i; 
          CloneWise_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,TRUSTED_SPECIAL,4) = j; 
          CloneWise_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,TRUSTED_SPECIAL,5) = SKIP_TARGET;   
          }  




        /* Finished calculating state SKIP_TARGET */ 


        } /* end of for each i position in strip */ 


      /* Special state START has no special to special movements */ 


      /* Special state TRUSTED_SPECIAL has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = CloneWise_DC_OPT_SHADOW_SPECIAL(mat,0,j,TRUSTED_SPECIAL);  


      /* Source START is a special source for TRUSTED_SPECIAL */ 
      temp = CloneWise_DC_OPT_SHADOW_SPECIAL(mat,0,j - 1,START) + (0) + (0);     
      if( temp > score ) {  
        score = temp;    
        /* Also got to propagate shadows  */ 
        for(k=0;k<7;k++) 
          CloneWise_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,TRUSTED_SPECIAL,k) = CloneWise_DC_OPT_SHADOW_SPECIAL_SP(mat,i - 0,j - 1,START,k);   
        }  


      /* Source TRUSTED_SPECIAL is a special source for TRUSTED_SPECIAL */ 
      temp = CloneWise_DC_OPT_SHADOW_SPECIAL(mat,0,j - 1,TRUSTED_SPECIAL) + ((mat->match->skip_jset[j]-1)) + (0);    
      if( temp > score ) {  
        score = temp;    
        /* Also got to propagate shadows  */ 
        for(k=0;k<7;k++) 
          CloneWise_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,TRUSTED_SPECIAL,k) = CloneWise_DC_OPT_SHADOW_SPECIAL_SP(mat,i - 0,j - 1,TRUSTED_SPECIAL,k); 
        }  


      /* Source MATCH for state TRUSTED_SPECIAL is not special... already calculated */ 
      /* Source SKIP_QUERY for state TRUSTED_SPECIAL is not special... already calculated */ 
      /* Source SKIP_TARGET for state TRUSTED_SPECIAL is not special... already calculated */ 
      /* Put back score... (now updated!) */ 
      CloneWise_DC_OPT_SHADOW_SPECIAL(mat,0,j,TRUSTED_SPECIAL) = score;  
      /* Finished updating state TRUSTED_SPECIAL */ 




      /* Special state END has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = CloneWise_DC_OPT_SHADOW_SPECIAL(mat,0,j,END);  


      /* Source TRUSTED_SPECIAL is a special source for END */ 
      temp = CloneWise_DC_OPT_SHADOW_SPECIAL(mat,0,j - 1,TRUSTED_SPECIAL) + (0) + (0);   
      if( temp > score ) {  
        score = temp;    
        /* Also got to propagate shadows  */ 
        for(k=0;k<7;k++) 
          CloneWise_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,k) = CloneWise_DC_OPT_SHADOW_SPECIAL_SP(mat,i - 0,j - 1,TRUSTED_SPECIAL,k); 
        }  


      /* Put back score... (now updated!) */ 
      CloneWise_DC_OPT_SHADOW_SPECIAL(mat,0,j,END) = score;  
      /* Finished updating state END */ 


      } /* end of for each j strip */ 
    free(score_pointers);    
    free(shadow_pointers);   
    return TRUE;     
}    


/* Function:  init_start_end_linear_CloneWise(mat)
 *
 * Descrip: No Description
 *
 * Arg:        mat [UNKN ] Undocumented argument [CloneWise *]
 *
 */
void init_start_end_linear_CloneWise(CloneWise * mat) 
{
    register int i;  
    register int j;  
    for(j=0;j<3;j++) {  
      for(i=(-1);i<mat->q->length;i++)   {  
        CloneWise_DC_SHADOW_MATRIX(mat,i,j,MATCH) = NEGI;    
        CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,0) = (-1);   
        CloneWise_DC_SHADOW_MATRIX(mat,i,j,SKIP_QUERY) = NEGI;   
        CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,SKIP_QUERY,0) = (-1);  
        CloneWise_DC_SHADOW_MATRIX(mat,i,j,SKIP_TARGET) = NEGI;  
        CloneWise_DC_SHADOW_MATRIX_SP(mat,i,j,SKIP_TARGET,0) = (-1); 
        }  
      }  


    for(j=(-1);j<mat->t->length;j++) {  
      CloneWise_DC_SHADOW_SPECIAL(mat,0,j,START) = 0;    
      CloneWise_DC_SHADOW_SPECIAL_SP(mat,0,j,START,0) = j;   
      CloneWise_DC_SHADOW_SPECIAL(mat,0,j,TRUSTED_SPECIAL) = NEGI;   
      CloneWise_DC_SHADOW_SPECIAL_SP(mat,0,j,TRUSTED_SPECIAL,0) = (-1);  
      CloneWise_DC_SHADOW_SPECIAL(mat,0,j,END) = NEGI;   
      CloneWise_DC_SHADOW_SPECIAL_SP(mat,0,j,END,0) = (-1);  
      }  


    return;  
}    


/* Function:  convert_PackAln_to_AlnBlock_CloneWise(pal)
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
AlnBlock * convert_PackAln_to_AlnBlock_CloneWise(PackAln * pal) 
{
    AlnConvertSet * acs; 
    AlnBlock * alb;  


    acs = AlnConvertSet_CloneWise(); 
    alb = AlnBlock_from_PackAln(acs,pal);    
    free_AlnConvertSet(acs); 
    return alb;  
}    


 static char * query_label[] = { "QUERY_MATCH","QUERY_MATCH_PAUSE","QUERY_SKIP","QUERY_PAUSE","NO_QUERY_PAUSE","END" };  
/* Function:  AlnConvertSet_CloneWise(void)
 *
 * Descrip: No Description
 *
 *
 * Return [UNKN ]  Undocumented return value [AlnConvertSet *]
 *
 */
 static char * target_label[] = { "TARGET_MATCH","TARGET_MATCH_PAUSE","TARGET_PAUSE","TARGET_SKIP","END" };  
AlnConvertSet * AlnConvertSet_CloneWise(void) 
{
    AlnConvertUnit * acu;    
    AlnConvertSet  * out;    


    out = AlnConvertSet_alloc_std(); 


    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = MATCH; 
    acu->state2 = MATCH;     
    acu->offi = 1;   
    acu->offj = 1;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[0];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = MATCH; 
    acu->state2 = MATCH;     
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[1];    
    acu->label2 = target_label[0];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = MATCH; 
    acu->state2 = MATCH;     
    acu->offi = 1;   
    acu->offj = 0;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[1];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = SKIP_QUERY;    
    acu->state2 = MATCH;     
    acu->offi = 1;   
    acu->offj = 1;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[0];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = SKIP_TARGET;   
    acu->state2 = MATCH;     
    acu->offi = 1;   
    acu->offj = 1;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[0];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = START + 3; 
    acu->is_from_special = TRUE; 
    acu->state2 = MATCH;     
    acu->offi = (-1);    
    acu->offj = 1;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[0];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = MATCH; 
    acu->state2 = SKIP_QUERY;    
    acu->offi = 1;   
    acu->offj = 0;   
    acu->label1 = query_label[2];    
    acu->label2 = target_label[2];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = SKIP_TARGET;   
    acu->state2 = SKIP_QUERY;    
    acu->offi = 1;   
    acu->offj = 0;   
    acu->label1 = query_label[2];    
    acu->label2 = target_label[2];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = SKIP_QUERY;    
    acu->state2 = SKIP_QUERY;    
    acu->offi = 1;   
    acu->offj = 0;   
    acu->label1 = query_label[2];    
    acu->label2 = target_label[2];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = START + 3; 
    acu->is_from_special = TRUE; 
    acu->state2 = SKIP_QUERY;    
    acu->offi = (-1);    
    acu->offj = 0;   
    acu->label1 = query_label[2];    
    acu->label2 = target_label[2];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = MATCH; 
    acu->state2 = SKIP_TARGET;   
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[3];    
    acu->label2 = target_label[3];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = SKIP_QUERY;    
    acu->state2 = SKIP_TARGET;   
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[3];    
    acu->label2 = target_label[3];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = SKIP_TARGET;   
    acu->state2 = SKIP_TARGET;   
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[3];    
    acu->label2 = target_label[3];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = START + 3; 
    acu->is_from_special = TRUE; 
    acu->state2 = SKIP_TARGET;   
    acu->offi = (-1);    
    acu->offj = 1;   
    acu->label1 = query_label[3];    
    acu->label2 = target_label[3];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = START + 3; 
    acu->state2 = TRUSTED_SPECIAL + 3;   
    acu->offi = (-1);    
    acu->offj = 1;   
    acu->label1 = query_label[4];    
    acu->label2 = target_label[3];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = TRUSTED_SPECIAL + 3;   
    acu->state2 = TRUSTED_SPECIAL + 3;   
    acu->offi = (-1);    
    acu->offj = 1;   
    acu->label1 = query_label[4];    
    acu->label2 = target_label[3];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = MATCH; 
    acu->state2 = TRUSTED_SPECIAL + 3;   
    acu->offi = (-1);    
    acu->offj = 0;   
    acu->label1 = query_label[4];    
    acu->label2 = target_label[3];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = SKIP_QUERY;    
    acu->state2 = TRUSTED_SPECIAL + 3;   
    acu->offi = (-1);    
    acu->offj = 0;   
    acu->label1 = query_label[4];    
    acu->label2 = target_label[3];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = SKIP_TARGET;   
    acu->state2 = TRUSTED_SPECIAL + 3;   
    acu->offi = (-1);    
    acu->offj = 0;   
    acu->label1 = query_label[4];    
    acu->label2 = target_label[3];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = TRUSTED_SPECIAL + 3;   
    acu->state2 = END + 3;   
    acu->offi = (-1);    
    acu->offj = 1;   
    acu->label1 = query_label[5];    
    acu->label2 = target_label[4];   
    return out;  
}    


/* Function:  PackAln_read_Expl_CloneWise(mat)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [CloneWise *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * PackAln_read_Expl_CloneWise(CloneWise * mat) 
{
    CloneWise_access_func_holder holder;     


    holder.access_main    = CloneWise_explicit_access_main;  
    holder.access_special = CloneWise_explicit_access_special;   
    return PackAln_read_generic_CloneWise(mat,holder);   
}    


/* Function:  CloneWise_explicit_access_main(mat,i,j,state)
 *
 * Descrip: No Description
 *
 * Arg:          mat [UNKN ] Undocumented argument [CloneWise *]
 * Arg:            i [UNKN ] Undocumented argument [int]
 * Arg:            j [UNKN ] Undocumented argument [int]
 * Arg:        state [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int CloneWise_explicit_access_main(CloneWise * mat,int i,int j,int state) 
{
    return CloneWise_EXPL_MATRIX(mat,i,j,state); 
}    


/* Function:  CloneWise_explicit_access_special(mat,i,j,state)
 *
 * Descrip: No Description
 *
 * Arg:          mat [UNKN ] Undocumented argument [CloneWise *]
 * Arg:            i [UNKN ] Undocumented argument [int]
 * Arg:            j [UNKN ] Undocumented argument [int]
 * Arg:        state [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int CloneWise_explicit_access_special(CloneWise * mat,int i,int j,int state) 
{
    return CloneWise_EXPL_SPECIAL(mat,i,j,state);    
}    


/* Function:  PackAln_read_generic_CloneWise(mat,h)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [CloneWise *]
 * Arg:          h [UNKN ] Undocumented argument [CloneWise_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * PackAln_read_generic_CloneWise(CloneWise * mat,CloneWise_access_func_holder h) 
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


    out->score =  find_end_CloneWise(mat,&i,&j,&state,&isspecial,h); 


    /* Add final end transition (at the moment we have not got the score! */ 
    if( (pau= PackAlnUnit_alloc()) == NULL  || add_PackAln(out,pau) == FALSE )   {  
      warn("Failed the first PackAlnUnit alloc, %d length of Alignment in CloneWise_basic_read, returning a mess.(Sorry!)",out->len);    
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
        max_calc_special_CloneWise(mat,i,j,state,isspecial,&i,&j,&state,&isspecial,&cellscore,h);    
      else   
        max_calc_CloneWise(mat,i,j,state,isspecial,&i,&j,&state,&isspecial,&cellscore,h);    
      if(i == CloneWise_READ_OFF_ERROR || j == CloneWise_READ_OFF_ERROR || state == CloneWise_READ_OFF_ERROR )   {  
        warn("Problem - hit bad read off system, exiting now");  
        break;   
        }  
      if( (pau= PackAlnUnit_alloc()) == NULL  || add_PackAln(out,pau) == FALSE ) {  
        warn("Failed a PackAlnUnit alloc, %d length of Alignment in CloneWise_basic_read, returning partial alignment",out->len);    
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


/* Function:  find_end_CloneWise(mat,ri,rj,state,isspecial,h)
 *
 * Descrip: No Description
 *
 * Arg:              mat [UNKN ] Undocumented argument [CloneWise *]
 * Arg:               ri [UNKN ] Undocumented argument [int *]
 * Arg:               rj [UNKN ] Undocumented argument [int *]
 * Arg:            state [UNKN ] Undocumented argument [int *]
 * Arg:        isspecial [UNKN ] Undocumented argument [boolean *]
 * Arg:                h [UNKN ] Undocumented argument [CloneWise_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int find_end_CloneWise(CloneWise * mat,int * ri,int * rj,int * state,boolean * isspecial,CloneWise_access_func_holder h) 
{
    int j;   
    int max; 
    int maxj;    
    int temp;    


    max = (*h.access_special)(mat,0,mat->t->length-1,END);   
    maxj = mat->t->length-1;     
    for(j= mat->t->length-2 ;j >= 0 ;j--)    {  
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


/* Function:  CloneWise_debug_show_matrix(mat,starti,stopi,startj,stopj,ofp)
 *
 * Descrip: No Description
 *
 * Arg:           mat [UNKN ] Undocumented argument [CloneWise *]
 * Arg:        starti [UNKN ] Undocumented argument [int]
 * Arg:         stopi [UNKN ] Undocumented argument [int]
 * Arg:        startj [UNKN ] Undocumented argument [int]
 * Arg:         stopj [UNKN ] Undocumented argument [int]
 * Arg:           ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
void CloneWise_debug_show_matrix(CloneWise * mat,int starti,int stopi,int startj,int stopj,FILE * ofp) 
{
    register int i;  
    register int j;  


    for(i=starti;i<stopi && i < mat->q->length;i++)  {  
      for(j=startj;j<stopj && j < mat->t->length;j++)    {  
        fprintf(ofp,"Cell [%d - %d]\n",i,j);     
        fprintf(ofp,"State MATCH %d\n",CloneWise_EXPL_MATRIX(mat,i,j,MATCH));    
        fprintf(ofp,"State SKIP_QUERY %d\n",CloneWise_EXPL_MATRIX(mat,i,j,SKIP_QUERY));  
        fprintf(ofp,"State SKIP_TARGET %d\n",CloneWise_EXPL_MATRIX(mat,i,j,SKIP_TARGET));    
        fprintf(ofp,"\n\n"); 
        }  
      }  


}    


/* Function:  max_calc_CloneWise(mat,i,j,state,isspecial,reti,retj,retstate,retspecial,cellscore,h)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [CloneWise *]
 * Arg:                 i [UNKN ] Undocumented argument [int]
 * Arg:                 j [UNKN ] Undocumented argument [int]
 * Arg:             state [UNKN ] Undocumented argument [int]
 * Arg:         isspecial [UNKN ] Undocumented argument [boolean]
 * Arg:              reti [UNKN ] Undocumented argument [int *]
 * Arg:              retj [UNKN ] Undocumented argument [int *]
 * Arg:          retstate [UNKN ] Undocumented argument [int *]
 * Arg:        retspecial [UNKN ] Undocumented argument [boolean *]
 * Arg:         cellscore [UNKN ] Undocumented argument [int *]
 * Arg:                 h [UNKN ] Undocumented argument [CloneWise_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int max_calc_CloneWise(CloneWise * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,CloneWise_access_func_holder h) 
{
    register int temp;   
    register int cscore; 


    *reti = (*retj) = (*retstate) = CloneWise_READ_OFF_ERROR;    


    if( i < 0 || j < 0 || i > mat->q->length || j > mat->t->length)  {  
      warn("In CloneWise matrix special read off - out of bounds on matrix [i,j is %d,%d state %d in standard matrix]",i,j,state);   
      return -1;     
      }  


    /* Then you have to select the correct switch statement to figure out the readoff      */ 
    /* Somewhat odd - reverse the order of calculation and return as soon as it is correct */ 
    cscore = (*h.access_main)(mat,i,j,state);    
    switch(state)    { /*Switch state */ 
      case MATCH :   
        /* Has restricted position */ 
        if( (i-1) == 0 && (j-1) == 0 )   {  
          temp = cscore - (0) -  (mat->match->matrix[i][j]); 
          if( temp == (*h.access_special)(mat,i - 1,j - 1,START) )   {  
            *reti = i - 1;   
            *retj = j - 1;   
            *retstate = START;   
            *retspecial = TRUE;  
            if( cellscore != NULL)   {  
              *cellscore = cscore - (*h.access_special)(mat,i-1,j-1,START);  
              }  
            return (*h.access_main)(mat,i - 1,j - 1,START);  
            }  
          }  
        temp = cscore - (0) -  (mat->match->matrix[i][j]);   
        if( temp == (*h.access_main)(mat,i - 1,j - 1,SKIP_TARGET) )  {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = SKIP_TARGET;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-1,SKIP_TARGET); 
            }  
          return (*h.access_main)(mat,i - 1,j - 1,SKIP_TARGET);  
          }  
        temp = cscore - (0) -  (mat->match->matrix[i][j]);   
        if( temp == (*h.access_main)(mat,i - 1,j - 1,SKIP_QUERY) )   {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = SKIP_QUERY;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-1,SKIP_QUERY);  
            }  
          return (*h.access_main)(mat,i - 1,j - 1,SKIP_QUERY);   
          }  
        temp = cscore - (0) -  (mat->match->matrix[i][j]);   
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
        temp = cscore - (0) -  (mat->match->matrix[i][j]);   
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
        temp = cscore - ((mat->match->matrix[i][j]+1)) -  (mat->match->matrix[i][j]);    
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
        warn("Major problem (!) - in CloneWise read off, position %d,%d state %d no source found!",i,j,state);   
        return (-1); 
      case SKIP_QUERY :  
        /* Has restricted position */ 
        if( (i-1) == 0 && (j-0) == 0 )   {  
          temp = cscore - (mat->query_skip_start) -  (mat->match->skip_iset[i]); 
          if( temp == (*h.access_special)(mat,i - 1,j - 0,START) )   {  
            *reti = i - 1;   
            *retj = j - 0;   
            *retstate = START;   
            *retspecial = TRUE;  
            if( cellscore != NULL)   {  
              *cellscore = cscore - (*h.access_special)(mat,i-1,j-0,START);  
              }  
            return (*h.access_main)(mat,i - 1,j - 0,START);  
            }  
          }  
        temp = cscore - (0) -  (mat->match->skip_iset[i]);   
        if( temp == (*h.access_main)(mat,i - 1,j - 0,SKIP_QUERY) )   {  
          *reti = i - 1; 
          *retj = j - 0; 
          *retstate = SKIP_QUERY;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-0,SKIP_QUERY);  
            }  
          return (*h.access_main)(mat,i - 1,j - 0,SKIP_QUERY);   
          }  
        temp = cscore - (mat->query_skip_start) -  (mat->match->skip_iset[i]);   
        if( temp == (*h.access_main)(mat,i - 1,j - 0,SKIP_TARGET) )  {  
          *reti = i - 1; 
          *retj = j - 0; 
          *retstate = SKIP_TARGET;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-0,SKIP_TARGET); 
            }  
          return (*h.access_main)(mat,i - 1,j - 0,SKIP_TARGET);  
          }  
        temp = cscore - (mat->query_skip_start) -  (mat->match->skip_iset[i]);   
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
        warn("Major problem (!) - in CloneWise read off, position %d,%d state %d no source found!",i,j,state);   
        return (-1); 
      case SKIP_TARGET :     
        /* Has restricted position */ 
        if( (i-0) == 0 && (j-1) == 0 )   {  
          temp = cscore - (mat->target_skip_start) -  (mat->match->skip_jset[j]);    
          if( temp == (*h.access_special)(mat,i - 0,j - 1,START) )   {  
            *reti = i - 0;   
            *retj = j - 1;   
            *retstate = START;   
            *retspecial = TRUE;  
            if( cellscore != NULL)   {  
              *cellscore = cscore - (*h.access_special)(mat,i-0,j-1,START);  
              }  
            return (*h.access_main)(mat,i - 0,j - 1,START);  
            }  
          }  
        temp = cscore - (0) -  (mat->match->skip_jset[j]);   
        if( temp == (*h.access_main)(mat,i - 0,j - 1,SKIP_TARGET) )  {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = SKIP_TARGET;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,SKIP_TARGET); 
            }  
          return (*h.access_main)(mat,i - 0,j - 1,SKIP_TARGET);  
          }  
        temp = cscore - (mat->target_skip_start) -  (mat->match->skip_jset[j]);  
        if( temp == (*h.access_main)(mat,i - 0,j - 1,SKIP_QUERY) )   {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = SKIP_QUERY;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,SKIP_QUERY);  
            }  
          return (*h.access_main)(mat,i - 0,j - 1,SKIP_QUERY);   
          }  
        temp = cscore - (mat->target_skip_start) -  (mat->match->skip_jset[j]);  
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
        warn("Major problem (!) - in CloneWise read off, position %d,%d state %d no source found!",i,j,state);   
        return (-1); 
      default:   
        warn("Major problem (!) - in CloneWise read off, position %d,%d state %d no source found!",i,j,state);   
        return (-1); 
      } /* end of Switch state  */ 
}    


/* Function:  max_calc_special_CloneWise(mat,i,j,state,isspecial,reti,retj,retstate,retspecial,cellscore,h)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [CloneWise *]
 * Arg:                 i [UNKN ] Undocumented argument [int]
 * Arg:                 j [UNKN ] Undocumented argument [int]
 * Arg:             state [UNKN ] Undocumented argument [int]
 * Arg:         isspecial [UNKN ] Undocumented argument [boolean]
 * Arg:              reti [UNKN ] Undocumented argument [int *]
 * Arg:              retj [UNKN ] Undocumented argument [int *]
 * Arg:          retstate [UNKN ] Undocumented argument [int *]
 * Arg:        retspecial [UNKN ] Undocumented argument [boolean *]
 * Arg:         cellscore [UNKN ] Undocumented argument [int *]
 * Arg:                 h [UNKN ] Undocumented argument [CloneWise_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int max_calc_special_CloneWise(CloneWise * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,CloneWise_access_func_holder h) 
{
    register int temp;   
    register int cscore; 


    *reti = (*retj) = (*retstate) = CloneWise_READ_OFF_ERROR;    


    if( j < 0 || j > mat->t->length) {  
      warn("In CloneWise matrix special read off - out of bounds on matrix [j is %d in special]",j); 
      return -1;     
      }  


    cscore = (*h.access_special)(mat,i,j,state); 
    switch(state)    { /*switch on special states*/ 
      case START :   
      case TRUSTED_SPECIAL :     
        /* source SKIP_TARGET is from main matrix */ 
        for(i= mat->q->length-1;i >= 0 ;i--) { /*for i >= 0*/ 
          temp = cscore - (mat->target_special_s) - (0);     
          if( temp == (*h.access_main)(mat,i - 0,j - 0,SKIP_TARGET) )    {  
            *reti = i - 0;   
            *retj = j - 0;   
            *retstate = SKIP_TARGET; 
            *retspecial = FALSE; 
            if( cellscore != NULL)   {  
              *cellscore = cscore - (*h.access_main)(mat,i-0,j-0,SKIP_TARGET);   
              }  
            return (*h.access_main)(mat,i - 0,j - 0,SKIP_TARGET) ;   
            }  
          } /* end of for i >= 0 */ 
        /* source SKIP_QUERY is from main matrix */ 
        for(i= mat->q->length-1;i >= 0 ;i--) { /*for i >= 0*/ 
          temp = cscore - (mat->target_special_s) - (0);     
          if( temp == (*h.access_main)(mat,i - 0,j - 0,SKIP_QUERY) ) {  
            *reti = i - 0;   
            *retj = j - 0;   
            *retstate = SKIP_QUERY;  
            *retspecial = FALSE; 
            if( cellscore != NULL)   {  
              *cellscore = cscore - (*h.access_main)(mat,i-0,j-0,SKIP_QUERY);    
              }  
            return (*h.access_main)(mat,i - 0,j - 0,SKIP_QUERY) ;    
            }  
          } /* end of for i >= 0 */ 
        /* source MATCH is from main matrix */ 
        for(i= mat->q->length-1;i >= 0 ;i--) { /*for i >= 0*/ 
          temp = cscore - (mat->target_special_s) - (0);     
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
        /* source TRUSTED_SPECIAL is a special */ 
        temp = cscore - ((mat->match->skip_jset[j]-1)) - (0);    
        if( temp == (*h.access_special)(mat,i - 0,j - 1,TRUSTED_SPECIAL) )   {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = TRUSTED_SPECIAL;   
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_special)(mat,i-0,j-1,TRUSTED_SPECIAL);  
            }  
          return (*h.access_special)(mat,i - 0,j - 1,TRUSTED_SPECIAL) ;  
          }  
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
      case END :     
        /* source TRUSTED_SPECIAL is a special */ 
        temp = cscore - (0) - (0);   
        if( temp == (*h.access_special)(mat,i - 0,j - 1,TRUSTED_SPECIAL) )   {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = TRUSTED_SPECIAL;   
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_special)(mat,i-0,j-1,TRUSTED_SPECIAL);  
            }  
          return (*h.access_special)(mat,i - 0,j - 1,TRUSTED_SPECIAL) ;  
          }  
      default:   
        warn("Major problem (!) - in CloneWise read off, position %d,%d state %d no source found  dropped into default on source switch!",i,j,state);    
        return (-1); 
      } /* end of switch on special states */ 
}    


/* Function:  calculate_CloneWise(mat)
 *
 * Descrip:    This function calculates the CloneWise matrix when in explicit mode
 *             To allocate the matrix use /allocate_Expl_CloneWise
 *
 *
 * Arg:        mat [UNKN ] CloneWise which contains explicit basematrix memory [CloneWise *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean calculate_CloneWise(CloneWise * mat) 
{
    int i;   
    int j;   
    int leni;    
    int lenj;    
    int tot; 
    int num; 


    if( mat->basematrix->type != BASEMATRIX_TYPE_EXPLICIT )  {  
      warn("in calculate_CloneWise, passed a non Explicit matrix type, cannot calculate!");  
      return FALSE;  
      }  


    leni = mat->leni;    
    lenj = mat->lenj;    
    tot = leni * lenj;   
    num = 0; 


    start_reporting("CloneWise Matrix calculation: ");   
    for(j=0;j<lenj;j++)  {  
      auto int score;    
      auto int temp;     
      for(i=0;i<leni;i++)    {  
        if( num%1000 == 0 )  
          log_full_error(REPORT,0,"[%7d] Cells %2d%%%%",num,num*100/tot);    
        num++;   


        /* For state MATCH */ 
        /* setting first movement to score */ 
        score = CloneWise_EXPL_MATRIX(mat,i-1,j-1,MATCH) + (mat->match->matrix[i][j]+1);     
        /* From state MATCH to state MATCH */ 
        temp = CloneWise_EXPL_MATRIX(mat,i-0,j-1,MATCH) + 0;     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = CloneWise_EXPL_MATRIX(mat,i-1,j-0,MATCH) + 0;     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SKIP_QUERY to state MATCH */ 
        temp = CloneWise_EXPL_MATRIX(mat,i-1,j-1,SKIP_QUERY) + 0;    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SKIP_TARGET to state MATCH */ 
        temp = CloneWise_EXPL_MATRIX(mat,i-1,j-1,SKIP_TARGET) + 0;   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* Has restricted position */ 
        if( (i-1) == 0 && (j-1) == 0 )   {  
          /* From state START to state MATCH */ 
          temp = CloneWise_EXPL_SPECIAL(mat,i-1,j-1,START) + 0;  
          if( temp  > score )    {  
            score = temp;    
            }  
          }  


        /* Ok - finished max calculation for MATCH */ 
        /* Add any movement independant score and put away */ 
         score += mat->match->matrix[i][j];  
         CloneWise_EXPL_MATRIX(mat,i,j,MATCH) = score;   


        /* state MATCH is a source for special TRUSTED_SPECIAL */ 
        temp = score + (mat->target_special_s) + (0) ;   
        if( temp > CloneWise_EXPL_SPECIAL(mat,i,j,TRUSTED_SPECIAL) )     {  
          CloneWise_EXPL_SPECIAL(mat,i,j,TRUSTED_SPECIAL) = temp;    
          }  




        /* Finished calculating state MATCH */ 


        /* For state SKIP_QUERY */ 
        /* setting first movement to score */ 
        score = CloneWise_EXPL_MATRIX(mat,i-1,j-0,MATCH) + mat->query_skip_start;    
        /* From state SKIP_TARGET to state SKIP_QUERY */ 
        temp = CloneWise_EXPL_MATRIX(mat,i-1,j-0,SKIP_TARGET) + mat->query_skip_start;   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SKIP_QUERY to state SKIP_QUERY */ 
        temp = CloneWise_EXPL_MATRIX(mat,i-1,j-0,SKIP_QUERY) + 0;    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* Has restricted position */ 
        if( (i-1) == 0 && (j-0) == 0 )   {  
          /* From state START to state SKIP_QUERY */ 
          temp = CloneWise_EXPL_SPECIAL(mat,i-1,j-0,START) + mat->query_skip_start;  
          if( temp  > score )    {  
            score = temp;    
            }  
          }  


        /* Ok - finished max calculation for SKIP_QUERY */ 
        /* Add any movement independant score and put away */ 
         score += mat->match->skip_iset[i];  
         CloneWise_EXPL_MATRIX(mat,i,j,SKIP_QUERY) = score;  


        /* state SKIP_QUERY is a source for special TRUSTED_SPECIAL */ 
        temp = score + (mat->target_special_s) + (0) ;   
        if( temp > CloneWise_EXPL_SPECIAL(mat,i,j,TRUSTED_SPECIAL) )     {  
          CloneWise_EXPL_SPECIAL(mat,i,j,TRUSTED_SPECIAL) = temp;    
          }  




        /* Finished calculating state SKIP_QUERY */ 


        /* For state SKIP_TARGET */ 
        /* setting first movement to score */ 
        score = CloneWise_EXPL_MATRIX(mat,i-0,j-1,MATCH) + mat->target_skip_start;   
        /* From state SKIP_QUERY to state SKIP_TARGET */ 
        temp = CloneWise_EXPL_MATRIX(mat,i-0,j-1,SKIP_QUERY) + mat->target_skip_start;   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SKIP_TARGET to state SKIP_TARGET */ 
        temp = CloneWise_EXPL_MATRIX(mat,i-0,j-1,SKIP_TARGET) + 0;   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* Has restricted position */ 
        if( (i-0) == 0 && (j-1) == 0 )   {  
          /* From state START to state SKIP_TARGET */ 
          temp = CloneWise_EXPL_SPECIAL(mat,i-0,j-1,START) + mat->target_skip_start;     
          if( temp  > score )    {  
            score = temp;    
            }  
          }  


        /* Ok - finished max calculation for SKIP_TARGET */ 
        /* Add any movement independant score and put away */ 
         score += mat->match->skip_jset[j];  
         CloneWise_EXPL_MATRIX(mat,i,j,SKIP_TARGET) = score; 


        /* state SKIP_TARGET is a source for special TRUSTED_SPECIAL */ 
        temp = score + (mat->target_special_s) + (0) ;   
        if( temp > CloneWise_EXPL_SPECIAL(mat,i,j,TRUSTED_SPECIAL) )     {  
          CloneWise_EXPL_SPECIAL(mat,i,j,TRUSTED_SPECIAL) = temp;    
          }  




        /* Finished calculating state SKIP_TARGET */ 
        }  


      /* Special state START has no special to special movements */ 


      /* Special state TRUSTED_SPECIAL has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = CloneWise_EXPL_SPECIAL(mat,0,j,TRUSTED_SPECIAL);   


      /* Source START is a special source for TRUSTED_SPECIAL */ 
      /* Has restricted position */ 
      if( (j-1) == 0  )  {  
        temp = CloneWise_EXPL_SPECIAL(mat,0,j - 1,START) + (0) + (0);    
        if( temp > score )   
          score = temp;  
        }  


      /* Source TRUSTED_SPECIAL is a special source for TRUSTED_SPECIAL */ 
      temp = CloneWise_EXPL_SPECIAL(mat,0,j - 1,TRUSTED_SPECIAL) + ((mat->match->skip_jset[j]-1)) + (0);     
      if( temp > score ) 
        score = temp;    


      /* Source MATCH for state TRUSTED_SPECIAL is not special... already calculated */ 
      /* Source SKIP_QUERY for state TRUSTED_SPECIAL is not special... already calculated */ 
      /* Source SKIP_TARGET for state TRUSTED_SPECIAL is not special... already calculated */ 
      /* Put back score... (now updated!) */ 
      CloneWise_EXPL_SPECIAL(mat,0,j,TRUSTED_SPECIAL) = score;   
      /* Finished updating state TRUSTED_SPECIAL */ 




      /* Special state END has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = CloneWise_EXPL_SPECIAL(mat,0,j,END);   


      /* Source TRUSTED_SPECIAL is a special source for END */ 
      /* Has restricted position */ 
      if( j == mat->lenj-1 ) {  
        temp = CloneWise_EXPL_SPECIAL(mat,0,j - 1,TRUSTED_SPECIAL) + (0) + (0);  
        if( temp > score )   
          score = temp;  
        }  


      /* Put back score... (now updated!) */ 
      CloneWise_EXPL_SPECIAL(mat,0,j,END) = score;   
      /* Finished updating state END */ 


      }  
    stop_reporting();    
    return TRUE;     
}    


/* Function:  calculate_dpenv_CloneWise(mat,dpenv)
 *
 * Descrip:    This function calculates the CloneWise matrix when in explicit mode, subject to the envelope
 *
 *
 * Arg:          mat [UNKN ] CloneWise which contains explicit basematrix memory [CloneWise *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean calculate_dpenv_CloneWise(CloneWise * mat,DPEnvelope * dpenv) 
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
      warn("in calculate_CloneWise, passed a non Explicit matrix type, cannot calculate!");  
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
        CloneWise_EXPL_MATRIX(mat,i,j,MATCH) = NEGI; 
        CloneWise_EXPL_MATRIX(mat,i,j,SKIP_QUERY) = NEGI;    
        CloneWise_EXPL_MATRIX(mat,i,j,SKIP_TARGET) = NEGI;   
        }  
      }  
    for(j=-1;j<mat->lenj;j++)    {  
      CloneWise_EXPL_SPECIAL(mat,i,j,START) = 0; 
      CloneWise_EXPL_SPECIAL(mat,i,j,TRUSTED_SPECIAL) = NEGI;    
      CloneWise_EXPL_SPECIAL(mat,i,j,END) = NEGI;    
      }  


    start_reporting("CloneWise Matrix calculation: ");   
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
          CloneWise_EXPL_MATRIX(mat,i,j,MATCH) = NEGI;   
          CloneWise_EXPL_MATRIX(mat,i,j,SKIP_QUERY) = NEGI;  
          CloneWise_EXPL_MATRIX(mat,i,j,SKIP_TARGET) = NEGI; 
          continue;  
          }  


        if( num%1000 == 0 )  
          log_full_error(REPORT,0,"[%7d] Cells %2d%%%%",num,num*100/tot);    
        num++;   


        /* For state MATCH */ 
        /* setting first movement to score */ 
        score = CloneWise_EXPL_MATRIX(mat,i-1,j-1,MATCH) + (mat->match->matrix[i][j]+1);     
        /* From state MATCH to state MATCH */ 
        temp = CloneWise_EXPL_MATRIX(mat,i-0,j-1,MATCH) + 0;     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = CloneWise_EXPL_MATRIX(mat,i-1,j-0,MATCH) + 0;     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SKIP_QUERY to state MATCH */ 
        temp = CloneWise_EXPL_MATRIX(mat,i-1,j-1,SKIP_QUERY) + 0;    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SKIP_TARGET to state MATCH */ 
        temp = CloneWise_EXPL_MATRIX(mat,i-1,j-1,SKIP_TARGET) + 0;   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* Has restricted position */ 
        if( (i-1) == 0 && (j-1) == 0 )   {  
          /* From state START to state MATCH */ 
          temp = CloneWise_EXPL_SPECIAL(mat,i-1,j-1,START) + 0;  
          if( temp  > score )    {  
            score = temp;    
            }  
          }  


        /* Ok - finished max calculation for MATCH */ 
        /* Add any movement independant score and put away */ 
         score += mat->match->matrix[i][j];  
         CloneWise_EXPL_MATRIX(mat,i,j,MATCH) = score;   


        /* state MATCH is a source for special TRUSTED_SPECIAL */ 
        temp = score + (mat->target_special_s) + (0) ;   
        if( temp > CloneWise_EXPL_SPECIAL(mat,i,j,TRUSTED_SPECIAL) )     {  
          CloneWise_EXPL_SPECIAL(mat,i,j,TRUSTED_SPECIAL) = temp;    
          }  




        /* Finished calculating state MATCH */ 


        /* For state SKIP_QUERY */ 
        /* setting first movement to score */ 
        score = CloneWise_EXPL_MATRIX(mat,i-1,j-0,MATCH) + mat->query_skip_start;    
        /* From state SKIP_TARGET to state SKIP_QUERY */ 
        temp = CloneWise_EXPL_MATRIX(mat,i-1,j-0,SKIP_TARGET) + mat->query_skip_start;   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SKIP_QUERY to state SKIP_QUERY */ 
        temp = CloneWise_EXPL_MATRIX(mat,i-1,j-0,SKIP_QUERY) + 0;    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* Has restricted position */ 
        if( (i-1) == 0 && (j-0) == 0 )   {  
          /* From state START to state SKIP_QUERY */ 
          temp = CloneWise_EXPL_SPECIAL(mat,i-1,j-0,START) + mat->query_skip_start;  
          if( temp  > score )    {  
            score = temp;    
            }  
          }  


        /* Ok - finished max calculation for SKIP_QUERY */ 
        /* Add any movement independant score and put away */ 
         score += mat->match->skip_iset[i];  
         CloneWise_EXPL_MATRIX(mat,i,j,SKIP_QUERY) = score;  


        /* state SKIP_QUERY is a source for special TRUSTED_SPECIAL */ 
        temp = score + (mat->target_special_s) + (0) ;   
        if( temp > CloneWise_EXPL_SPECIAL(mat,i,j,TRUSTED_SPECIAL) )     {  
          CloneWise_EXPL_SPECIAL(mat,i,j,TRUSTED_SPECIAL) = temp;    
          }  




        /* Finished calculating state SKIP_QUERY */ 


        /* For state SKIP_TARGET */ 
        /* setting first movement to score */ 
        score = CloneWise_EXPL_MATRIX(mat,i-0,j-1,MATCH) + mat->target_skip_start;   
        /* From state SKIP_QUERY to state SKIP_TARGET */ 
        temp = CloneWise_EXPL_MATRIX(mat,i-0,j-1,SKIP_QUERY) + mat->target_skip_start;   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SKIP_TARGET to state SKIP_TARGET */ 
        temp = CloneWise_EXPL_MATRIX(mat,i-0,j-1,SKIP_TARGET) + 0;   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* Has restricted position */ 
        if( (i-0) == 0 && (j-1) == 0 )   {  
          /* From state START to state SKIP_TARGET */ 
          temp = CloneWise_EXPL_SPECIAL(mat,i-0,j-1,START) + mat->target_skip_start;     
          if( temp  > score )    {  
            score = temp;    
            }  
          }  


        /* Ok - finished max calculation for SKIP_TARGET */ 
        /* Add any movement independant score and put away */ 
         score += mat->match->skip_jset[j];  
         CloneWise_EXPL_MATRIX(mat,i,j,SKIP_TARGET) = score; 


        /* state SKIP_TARGET is a source for special TRUSTED_SPECIAL */ 
        temp = score + (mat->target_special_s) + (0) ;   
        if( temp > CloneWise_EXPL_SPECIAL(mat,i,j,TRUSTED_SPECIAL) )     {  
          CloneWise_EXPL_SPECIAL(mat,i,j,TRUSTED_SPECIAL) = temp;    
          }  




        /* Finished calculating state SKIP_TARGET */ 
        }  


      /* Special state START has no special to special movements */ 


      /* Special state TRUSTED_SPECIAL has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = CloneWise_EXPL_SPECIAL(mat,0,j,TRUSTED_SPECIAL);   


      /* Source START is a special source for TRUSTED_SPECIAL */ 
      /* Has restricted position */ 
      if( (j-1) == 0  )  {  
        temp = CloneWise_EXPL_SPECIAL(mat,0,j - 1,START) + (0) + (0);    
        if( temp > score )   
          score = temp;  
        }  


      /* Source TRUSTED_SPECIAL is a special source for TRUSTED_SPECIAL */ 
      temp = CloneWise_EXPL_SPECIAL(mat,0,j - 1,TRUSTED_SPECIAL) + ((mat->match->skip_jset[j]-1)) + (0);     
      if( temp > score ) 
        score = temp;    


      /* Source MATCH for state TRUSTED_SPECIAL is not special... already calculated */ 
      /* Source SKIP_QUERY for state TRUSTED_SPECIAL is not special... already calculated */ 
      /* Source SKIP_TARGET for state TRUSTED_SPECIAL is not special... already calculated */ 
      /* Put back score... (now updated!) */ 
      CloneWise_EXPL_SPECIAL(mat,0,j,TRUSTED_SPECIAL) = score;   
      /* Finished updating state TRUSTED_SPECIAL */ 




      /* Special state END has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = CloneWise_EXPL_SPECIAL(mat,0,j,END);   


      /* Source TRUSTED_SPECIAL is a special source for END */ 
      /* Has restricted position */ 
      if( j == mat->lenj-1 ) {  
        temp = CloneWise_EXPL_SPECIAL(mat,0,j - 1,TRUSTED_SPECIAL) + (0) + (0);  
        if( temp > score )   
          score = temp;  
        }  


      /* Put back score... (now updated!) */ 
      CloneWise_EXPL_SPECIAL(mat,0,j,END) = score;   
      /* Finished updating state END */ 


      }  
    stop_reporting();    
    return TRUE;     
}    


/* Function:  CloneWise_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [CloneWise *]
 *
 */
CloneWise * CloneWise_alloc(void) 
{
    CloneWise * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(CloneWise *) ckalloc (sizeof(CloneWise))) == NULL)  {  
      warn("CloneWise_alloc failed ");   
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


/* Function:  free_CloneWise(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [CloneWise *]
 *
 * Return [UNKN ]  Undocumented return value [CloneWise *]
 *
 */
CloneWise * free_CloneWise(CloneWise * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a CloneWise obj. Should be trappable"); 
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
    /* obj->match is linked in */ 
    /* obj->target_skip_start is linked in */ 
    /* obj->target_skip is linked in */ 
    /* obj->query_skip_start is linked in */ 
    /* obj->query_skip is linked in */ 
    /* obj->spread is linked in */ 
    /* obj->target_special_s is linked in */ 


    ckfree(obj); 
    return NULL; 
}    





#ifdef _cplusplus
}
#endif
