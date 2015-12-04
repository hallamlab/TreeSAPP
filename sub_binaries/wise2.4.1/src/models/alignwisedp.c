#ifdef _cplusplus
extern "C" {
#endif
#include "alignwisedp.h"



# line 7 "alignwisedp.c"
/* Function:  hard_link_AlignGeneModelFrame(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [AlignGeneModelFrame *]
 *
 * Return [UNKN ]  Undocumented return value [AlignGeneModelFrame *]
 *
 */
AlignGeneModelFrame * hard_link_AlignGeneModelFrame(AlignGeneModelFrame * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a AlignGeneModelFrame object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  AlignGeneModelFrame_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [AlignGeneModelFrame *]
 *
 */
AlignGeneModelFrame * AlignGeneModelFrame_alloc(void) 
{
    AlignGeneModelFrame * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(AlignGeneModelFrame *) ckalloc (sizeof(AlignGeneModelFrame))) == NULL)  {  
      warn("AlignGeneModelFrame_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->len = 0;    


    return out;  
}    


/* Function:  free_AlignGeneModelFrame(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [AlignGeneModelFrame *]
 *
 * Return [UNKN ]  Undocumented return value [AlignGeneModelFrame *]
 *
 */
AlignGeneModelFrame * free_AlignGeneModelFrame(AlignGeneModelFrame * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a AlignGeneModelFrame obj. Should be trappable");   
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


    ckfree(obj); 
    return NULL; 
}    




  /*****************   C functions  ****************/
  /*             Written using dynamite            */
  /*            Sat Sep  8 09:05:31 2007           */
  /*            email birney@sanger.ac.uk          */
  /* http://www.sanger.ac.uk/Users/birney/dynamite */
  /*************************************************/


  /* Please report any problems or bugs to         */
  /* Ewan Birney, birney@sanger.ac.uk              */


/* basic set of macros to map states to numbers */ 
#define CODON 0  
#define INTRON_0 1   
#define INTRON_1 2   
#define INTRON_2 3   
#define REVERSE_CODON 4  
#define REVERSE_INTRON_0 5   
#define REVERSE_INTRON_1 6   
#define REVERSE_INTRON_2 7   


#define START 0  
#define INTERGENIC 1 
#define END 2    


#define AlignWise_EXPL_MATRIX(this_matrix,i,j,STATE) this_matrix->basematrix->matrix[((j+10)*8)+STATE][i+0]  
#define AlignWise_EXPL_SPECIAL(matrix,i,j,STATE) matrix->basematrix->specmatrix[STATE][j+10] 
#define AlignWise_READ_OFF_ERROR -11
    


#define AlignWise_VSMALL_MATRIX(mat,i,j,STATE) mat->basematrix->matrix[(j+11)%11][((i+0)*8)+STATE]   
#define AlignWise_VSMALL_SPECIAL(mat,i,j,STATE) mat->basematrix->specmatrix[(j+11)%11][STATE]    




#define AlignWise_SHATTER_SPECIAL(matrix,i,j,STATE) matrix->shatter->special[STATE][j]   
#define AlignWise_SHATTER_MATRIX(matrix,i,j,STATE)  fetch_cell_value_ShatterMatrix(mat->shatter,i,j,STATE)   


/* Function:  PackAln_read_Shatter_AlignWise(mat)
 *
 * Descrip:    Reads off PackAln from shatter matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [AlignWise *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * PackAln_read_Shatter_AlignWise(AlignWise * mat) 
{
    AlignWise_access_func_holder holder;     


    holder.access_main    = AlignWise_shatter_access_main;   
    holder.access_special = AlignWise_shatter_access_special;    
    assert(mat);     
    assert(mat->shatter);    
    return PackAln_read_generic_AlignWise(mat,holder);   
}    


/* Function:  AlignWise_shatter_access_main(mat,i,j,state)
 *
 * Descrip: No Description
 *
 * Arg:          mat [UNKN ] Undocumented argument [AlignWise *]
 * Arg:            i [UNKN ] Undocumented argument [int]
 * Arg:            j [UNKN ] Undocumented argument [int]
 * Arg:        state [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int AlignWise_shatter_access_main(AlignWise * mat,int i,int j,int state) 
{
    return AlignWise_SHATTER_MATRIX(mat,i,j,state);  
}    


/* Function:  AlignWise_shatter_access_special(mat,i,j,state)
 *
 * Descrip: No Description
 *
 * Arg:          mat [UNKN ] Undocumented argument [AlignWise *]
 * Arg:            i [UNKN ] Undocumented argument [int]
 * Arg:            j [UNKN ] Undocumented argument [int]
 * Arg:        state [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int AlignWise_shatter_access_special(AlignWise * mat,int i,int j,int state) 
{
    return AlignWise_SHATTER_SPECIAL(mat,i,j,state); 
}    


/* Function:  calculate_shatter_AlignWise(mat,dpenv)
 *
 * Descrip:    This function calculates the AlignWise matrix when in shatter mode
 *
 *
 * Arg:          mat [UNKN ] (null) [AlignWise *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean calculate_shatter_AlignWise(AlignWise * mat,DPEnvelope * dpenv) 
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
    int * SIG_0_3;   
    int * SIG_0_6;   
    int * SIG_0_5;   
    int * SIG_0_4;   
    int * SIG_0_8;   
    int * SIG_0_1;   
    int * SIG_0_9;   
    int * SIG_0_10;  


    leni = mat->leni;    
    lenj = mat->lenj;    


    mat->shatter = new_ShatterMatrix(dpenv,8,lenj,3);    
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


    start_reporting("AlignWise Matrix calculation: ");   
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
        SIG_0_3 = fetch_cell_from_ShatterMatrix(mat->shatter,i-0,j-3);   
        SIG_0_6 = fetch_cell_from_ShatterMatrix(mat->shatter,i-0,j-6);   
        SIG_0_5 = fetch_cell_from_ShatterMatrix(mat->shatter,i-0,j-5);   
        SIG_0_4 = fetch_cell_from_ShatterMatrix(mat->shatter,i-0,j-4);   
        SIG_0_8 = fetch_cell_from_ShatterMatrix(mat->shatter,i-0,j-8);   
        SIG_0_1 = fetch_cell_from_ShatterMatrix(mat->shatter,i-0,j-1);   
        SIG_0_9 = fetch_cell_from_ShatterMatrix(mat->shatter,i-0,j-9);   
        SIG_0_10 = fetch_cell_from_ShatterMatrix(mat->shatter,i-0,j-10); 




        /* For state CODON */ 
        /* setting first movement to score */ 
        score = SIG_0_3[CODON] + mat->align->forward_coding[j];  
        /* From state INTERGENIC to state CODON */ 
        temp = AlignWise_SHATTER_SPECIAL(mat,i-0,j-3,INTERGENIC) + mat->align->forward_coding[j];    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTRON_0 to state CODON */ 
        temp = SIG_0_6[INTRON_0] + (mat->align->splice3_forward[(j-3)]+mat->align->forward_coding[j]);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTRON_1 to state CODON */ 
        temp = SIG_0_5[INTRON_1] + mat->align->splice3_forward[(j-2)];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTRON_2 to state CODON */ 
        temp = SIG_0_4[INTRON_2] + mat->align->splice3_forward[(j-1)];   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for CODON */ 
        /* Add any movement independant score and put away */ 
         SIG_0_0[CODON] = score; 


        /* state CODON is a source for special INTERGENIC */ 
        temp = score + (mat->geneopen) + (0) ;   
        if( temp > AlignWise_SHATTER_SPECIAL(mat,i,j,INTERGENIC) )   {  
          AlignWise_SHATTER_SPECIAL(mat,i,j,INTERGENIC) = temp;  
          }  




        /* Finished calculating state CODON */ 


        /* For state INTRON_0 */ 
        /* setting first movement to score */ 
        score = SIG_0_8[CODON] + (mat->align->splice5_forward[(j-7)]+mat->intronopen);   
        /* From state INTRON_0 to state INTRON_0 */ 
        temp = SIG_0_1[INTRON_0] + 0;    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for INTRON_0 */ 
        /* Add any movement independant score and put away */ 
         SIG_0_0[INTRON_0] = score;  


        /* Finished calculating state INTRON_0 */ 


        /* For state INTRON_1 */ 
        /* setting first movement to score */ 
        score = SIG_0_9[CODON] + (mat->align->splice5_forward[(j-7)]+mat->intronopen);   
        /* From state INTRON_1 to state INTRON_1 */ 
        temp = SIG_0_1[INTRON_1] + 0;    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for INTRON_1 */ 
        /* Add any movement independant score and put away */ 
         SIG_0_0[INTRON_1] = score;  


        /* Finished calculating state INTRON_1 */ 


        /* For state INTRON_2 */ 
        /* setting first movement to score */ 
        score = SIG_0_10[CODON] + (mat->align->splice5_forward[(j-7)]+mat->intronopen);  
        /* From state INTRON_2 to state INTRON_2 */ 
        temp = SIG_0_1[INTRON_2] + 0;    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for INTRON_2 */ 
        /* Add any movement independant score and put away */ 
         SIG_0_0[INTRON_2] = score;  


        /* Finished calculating state INTRON_2 */ 


        /* For state REVERSE_CODON */ 
        /* setting first movement to score */ 
        score = SIG_0_3[REVERSE_CODON] + mat->align->reverse_coding[j];  
        /* From state INTERGENIC to state REVERSE_CODON */ 
        temp = AlignWise_SHATTER_SPECIAL(mat,i-0,j-3,INTERGENIC) + mat->align->reverse_coding[j];    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state REVERSE_INTRON_0 to state REVERSE_CODON */ 
        temp = SIG_0_6[REVERSE_INTRON_0] + mat->align->splice5_reverse[(j-3)];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state REVERSE_INTRON_1 to state REVERSE_CODON */ 
        temp = SIG_0_5[REVERSE_INTRON_1] + mat->align->splice5_reverse[(j-2)];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state REVERSE_INTRON_2 to state REVERSE_CODON */ 
        temp = SIG_0_4[REVERSE_INTRON_2] + mat->align->splice5_reverse[(j-1)];   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for REVERSE_CODON */ 
        /* Add any movement independant score and put away */ 
         SIG_0_0[REVERSE_CODON] = score; 


        /* state REVERSE_CODON is a source for special INTERGENIC */ 
        temp = score + (mat->geneopen) + (0) ;   
        if( temp > AlignWise_SHATTER_SPECIAL(mat,i,j,INTERGENIC) )   {  
          AlignWise_SHATTER_SPECIAL(mat,i,j,INTERGENIC) = temp;  
          }  




        /* Finished calculating state REVERSE_CODON */ 


        /* For state REVERSE_INTRON_0 */ 
        /* setting first movement to score */ 
        score = SIG_0_8[REVERSE_CODON] + (mat->align->splice3_reverse[(j-7)]+mat->intronopen);   
        /* From state REVERSE_INTRON_0 to state REVERSE_INTRON_0 */ 
        temp = SIG_0_1[REVERSE_INTRON_0] + 0;    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for REVERSE_INTRON_0 */ 
        /* Add any movement independant score and put away */ 
         SIG_0_0[REVERSE_INTRON_0] = score;  


        /* Finished calculating state REVERSE_INTRON_0 */ 


        /* For state REVERSE_INTRON_1 */ 
        /* setting first movement to score */ 
        score = SIG_0_9[REVERSE_CODON] + (mat->align->splice3_reverse[(j-7)]+mat->intronopen);   
        /* From state REVERSE_INTRON_1 to state REVERSE_INTRON_1 */ 
        temp = SIG_0_1[REVERSE_INTRON_1] + 0;    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for REVERSE_INTRON_1 */ 
        /* Add any movement independant score and put away */ 
         SIG_0_0[REVERSE_INTRON_1] = score;  


        /* Finished calculating state REVERSE_INTRON_1 */ 


        /* For state REVERSE_INTRON_2 */ 
        /* setting first movement to score */ 
        score = SIG_0_10[REVERSE_CODON] + (mat->align->splice3_reverse[(j-7)]+mat->intronopen);  
        /* From state REVERSE_INTRON_2 to state REVERSE_INTRON_2 */ 
        temp = SIG_0_1[REVERSE_INTRON_2] + 0;    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for REVERSE_INTRON_2 */ 
        /* Add any movement independant score and put away */ 
         SIG_0_0[REVERSE_INTRON_2] = score;  


        /* Finished calculating state REVERSE_INTRON_2 */ 
        }  


      /* Special state START has no special to special movements */ 


      /* Special state INTERGENIC has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = AlignWise_SHATTER_SPECIAL(mat,0,j,INTERGENIC); 


      /* Source START is a special source for INTERGENIC */ 
      temp = AlignWise_SHATTER_SPECIAL(mat,0,j - 1,START) + (0) + (0);   
      if( temp > score ) 
        score = temp;    


      /* Source INTERGENIC is a special source for INTERGENIC */ 
      temp = AlignWise_SHATTER_SPECIAL(mat,0,j - 1,INTERGENIC) + (0) + (0);  
      if( temp > score ) 
        score = temp;    


      /* Source CODON for state INTERGENIC is not special... already calculated */ 
      /* Source REVERSE_CODON for state INTERGENIC is not special... already calculated */ 
      /* Put back score... (now updated!) */ 
      AlignWise_SHATTER_SPECIAL(mat,0,j,INTERGENIC) = score; 
      /* Finished updating state INTERGENIC */ 




      /* Special state END has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = AlignWise_SHATTER_SPECIAL(mat,0,j,END);    


      /* Source INTERGENIC is a special source for END */ 
      temp = AlignWise_SHATTER_SPECIAL(mat,0,j - 1,INTERGENIC) + (0) + (0);  
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      AlignWise_SHATTER_SPECIAL(mat,0,j,END) = score;    
      /* Finished updating state END */ 


      }  
    stop_reporting();    
    return TRUE;     
}    


/* Function:  search_AlignWise(dbsi,out,model,align,intronopen,geneopen)
 *
 * Descrip:    This function makes a database search of AlignWise
 *             It uses the dbsi structure to choose which implementation
 *             to use of the database searching. This way at run time you
 *             can switch between single threaded/multi-threaded or hardware
 *
 *
 * Arg:              dbsi [UNKN ] Undocumented argument [DBSearchImpl *]
 * Arg:               out [UNKN ] Undocumented argument [Hscore *]
 * Arg:             model [UNKN ] Undocumented argument [AlignGeneModelFrame*]
 * Arg:             align [UNKN ] Undocumented argument [AlignGeneModelScore*]
 * Arg:        intronopen [UNKN ] Undocumented argument [Score]
 * Arg:          geneopen [UNKN ] Undocumented argument [Score]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type search_AlignWise(DBSearchImpl * dbsi,Hscore * out,AlignGeneModelFrame* model,AlignGeneModelScore* align ,Score intronopen,Score geneopen) 
{
#ifdef PTHREAD   
    int i;   
    int thr_no;  
    pthread_attr_t pat;  
    struct thread_pool_holder_AlignWise * holder;    
#endif   
    if( out == NULL )    {  
      warn("Passed in a null Hscore object into search_AlignWise. Can't process results!");  
      return SEARCH_ERROR;   
      }  
    if( dbsi == NULL )   {  
      warn("Passed in a null DBSearchImpl object into search_AlignWise. Can't process results!");    
      return SEARCH_ERROR;   
      }  
    if( dbsi->trace_level > 5 )  
      warn("Asking for trace level of %d in database search for AlignWise, but it was compiled with a trace level of -2139062144. Not all trace statements can be shown",dbsi->trace_level); 
    switch(dbsi->type)   { /*switch on implementation*/ 
      case DBSearchImpl_Serial : 
        return serial_search_AlignWise(out,model,align ,intronopen,geneopen);    
      case DBSearchImpl_Pthreads :   
#ifdef PTHREAD   
        holder = (struct thread_pool_holder_AlignWise *) ckalloc(sizeof(struct thread_pool_holder_AlignWise));   
        if( holder == NULL )     {  
          warn("Unable to allocated thread pool datastructure...");  
          return SEARCH_ERROR;   
          }  
        holder->out = out;   
        holder->dbsi = dbsi; 
        holder->model = model;   
        holder->align = align;   
        holder->intronopen = intronopen; 
        holder->geneopen = geneopen; 
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
          if( pthread_create(holder->pool+i,&pat,thread_loop_AlignWise,(void *)holder) ) 
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
        warn("You did not specifiy the PTHREAD compile when compiled the C code for AlignWise"); 
#endif /* finished threads */    
      default :  
        warn("database search implementation %s was not provided in the compiled dynamite file from AlignWise",impl_string_DBSearchImpl(dbsi));  
        return SEARCH_ERROR; 
      } /* end of switch on implementation */ 


}    


/* Function:  thread_loop_AlignWise(ptr)
 *
 * Descrip:    dummy loop code foreach thread for AlignWise
 *
 *
 * Arg:        ptr [UNKN ] Undocumented argument [void *]
 *
 * Return [UNKN ]  Undocumented return value [void *]
 *
 */
void * thread_loop_AlignWise(void * ptr) 
{
    fatal("dummy thread loop function"); 
}    


/* Function:  serial_search_AlignWise(out,model,align,intronopen,geneopen)
 *
 * Descrip:    This function makes a database search of AlignWise
 *             It is a single processor implementation
 *
 *
 * Arg:               out [UNKN ] Undocumented argument [Hscore *]
 * Arg:             model [UNKN ] Undocumented argument [AlignGeneModelFrame*]
 * Arg:             align [UNKN ] Undocumented argument [AlignGeneModelScore*]
 * Arg:        intronopen [UNKN ] Undocumented argument [Score]
 * Arg:          geneopen [UNKN ] Undocumented argument [Score]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type serial_search_AlignWise(Hscore * out,AlignGeneModelFrame* model,AlignGeneModelScore* align ,Score intronopen,Score geneopen) 
{
    int db_status;   
    int score;   
    int query_pos = 0;   
    int target_pos = 0;  
    DataScore * ds;  


    push_errormsg_stack("Before any actual search in db searching"); 


    target_pos = 0;  




    /* No maximum length - allocated on-the-fly */ 
    score = score_only_AlignWise(model, align , intronopen, geneopen);   
    if( should_store_Hscore(out,score) == TRUE )     { /*if storing datascore*/ 
      ds = new_DataScore_from_storage(out);  
      if( ds == NULL )   {  
        warn("AlignWise search had a memory error in allocating a new_DataScore (?a leak somewhere - DataScore is a very small datastructure");  
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


/* Function:  score_only_AlignWise(model,align,intronopen,geneopen)
 *
 * Descrip:    This function just calculates the score for the matrix
 *             I am pretty sure we can do this better, but hey, for the moment...
 *             It calls /allocate_AlignWise_only
 *
 *
 * Arg:             model [UNKN ] query data structure [AlignGeneModelFrame*]
 * Arg:             align [UNKN ] target data structure [AlignGeneModelScore*]
 * Arg:        intronopen [UNKN ] Resource [Score]
 * Arg:          geneopen [UNKN ] Resource [Score]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int score_only_AlignWise(AlignGeneModelFrame* model,AlignGeneModelScore* align ,Score intronopen,Score geneopen) 
{
    int bestscore = NEGI;    
    int i;   
    int j;   
    int k;   
    AlignWise * mat;     


    mat = allocate_AlignWise_only(model, align , intronopen, geneopen);  
    if( mat == NULL )    {  
      warn("Memory allocation error in the db search - unable to communicate to calling function. this spells DIASTER!");    
      return NEGI;   
      }  
    if((mat->basematrix = BaseMatrix_alloc_matrix_and_specials(11,(mat->leni + 0) * 8,11,3)) == NULL)    {  
      warn("Score only matrix for AlignWise cannot be allocated, (asking for 10  by %d  cells)",mat->leni*8);    
      mat = free_AlignWise(mat);     
      return 0;  
      }  
    mat->basematrix->type = BASEMATRIX_TYPE_VERYSMALL;   


    /* Now, initiate matrix */ 
    for(j=0;j<12;j++)    {  
      for(i=(-0);i<mat->leni;i++)    {  
        for(k=0;k<8;k++) 
          AlignWise_VSMALL_MATRIX(mat,i,j,k) = NEGI; 
        }  
      AlignWise_VSMALL_SPECIAL(mat,i,j,START) = 0;   
      AlignWise_VSMALL_SPECIAL(mat,i,j,INTERGENIC) = NEGI;   
      AlignWise_VSMALL_SPECIAL(mat,i,j,END) = NEGI;  
      }  


    /* Ok, lets do-o-o-o-o it */ 


    for(j=0;j<mat->lenj;j++) { /*for all target positions*/ 
      auto int score;    
      auto int temp;     
      for(i=0;i<mat->leni;i++)   { /*for all query positions*/ 


        /* For state CODON */ 
        /* setting first movement to score */ 
        score = AlignWise_VSMALL_MATRIX(mat,i-0,j-3,CODON) + mat->align->forward_coding[j];  
        /* From state INTERGENIC to state CODON */ 
        temp = AlignWise_VSMALL_SPECIAL(mat,i-0,j-3,INTERGENIC) + mat->align->forward_coding[j];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTRON_0 to state CODON */ 
        temp = AlignWise_VSMALL_MATRIX(mat,i-0,j-6,INTRON_0) + (mat->align->splice3_forward[(j-3)]+mat->align->forward_coding[j]);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTRON_1 to state CODON */ 
        temp = AlignWise_VSMALL_MATRIX(mat,i-0,j-5,INTRON_1) + mat->align->splice3_forward[(j-2)];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTRON_2 to state CODON */ 
        temp = AlignWise_VSMALL_MATRIX(mat,i-0,j-4,INTRON_2) + mat->align->splice3_forward[(j-1)];   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for CODON */ 
        /* Add any movement independant score and put away */ 
         AlignWise_VSMALL_MATRIX(mat,i,j,CODON) = score; 


        /* state CODON is a source for special INTERGENIC */ 
        temp = score + (mat->geneopen) + (0) ;   
        if( temp > AlignWise_VSMALL_SPECIAL(mat,i,j,INTERGENIC) )    {  
          AlignWise_VSMALL_SPECIAL(mat,i,j,INTERGENIC) = temp;   
          }  




        /* Finished calculating state CODON */ 


        /* For state INTRON_0 */ 
        /* setting first movement to score */ 
        score = AlignWise_VSMALL_MATRIX(mat,i-0,j-8,CODON) + (mat->align->splice5_forward[(j-7)]+mat->intronopen);   
        /* From state INTRON_0 to state INTRON_0 */ 
        temp = AlignWise_VSMALL_MATRIX(mat,i-0,j-1,INTRON_0) + 0;    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for INTRON_0 */ 
        /* Add any movement independant score and put away */ 
         AlignWise_VSMALL_MATRIX(mat,i,j,INTRON_0) = score;  


        /* Finished calculating state INTRON_0 */ 


        /* For state INTRON_1 */ 
        /* setting first movement to score */ 
        score = AlignWise_VSMALL_MATRIX(mat,i-0,j-9,CODON) + (mat->align->splice5_forward[(j-7)]+mat->intronopen);   
        /* From state INTRON_1 to state INTRON_1 */ 
        temp = AlignWise_VSMALL_MATRIX(mat,i-0,j-1,INTRON_1) + 0;    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for INTRON_1 */ 
        /* Add any movement independant score and put away */ 
         AlignWise_VSMALL_MATRIX(mat,i,j,INTRON_1) = score;  


        /* Finished calculating state INTRON_1 */ 


        /* For state INTRON_2 */ 
        /* setting first movement to score */ 
        score = AlignWise_VSMALL_MATRIX(mat,i-0,j-10,CODON) + (mat->align->splice5_forward[(j-7)]+mat->intronopen);  
        /* From state INTRON_2 to state INTRON_2 */ 
        temp = AlignWise_VSMALL_MATRIX(mat,i-0,j-1,INTRON_2) + 0;    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for INTRON_2 */ 
        /* Add any movement independant score and put away */ 
         AlignWise_VSMALL_MATRIX(mat,i,j,INTRON_2) = score;  


        /* Finished calculating state INTRON_2 */ 


        /* For state REVERSE_CODON */ 
        /* setting first movement to score */ 
        score = AlignWise_VSMALL_MATRIX(mat,i-0,j-3,REVERSE_CODON) + mat->align->reverse_coding[j];  
        /* From state INTERGENIC to state REVERSE_CODON */ 
        temp = AlignWise_VSMALL_SPECIAL(mat,i-0,j-3,INTERGENIC) + mat->align->reverse_coding[j];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state REVERSE_INTRON_0 to state REVERSE_CODON */ 
        temp = AlignWise_VSMALL_MATRIX(mat,i-0,j-6,REVERSE_INTRON_0) + mat->align->splice5_reverse[(j-3)];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state REVERSE_INTRON_1 to state REVERSE_CODON */ 
        temp = AlignWise_VSMALL_MATRIX(mat,i-0,j-5,REVERSE_INTRON_1) + mat->align->splice5_reverse[(j-2)];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state REVERSE_INTRON_2 to state REVERSE_CODON */ 
        temp = AlignWise_VSMALL_MATRIX(mat,i-0,j-4,REVERSE_INTRON_2) + mat->align->splice5_reverse[(j-1)];   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for REVERSE_CODON */ 
        /* Add any movement independant score and put away */ 
         AlignWise_VSMALL_MATRIX(mat,i,j,REVERSE_CODON) = score; 


        /* state REVERSE_CODON is a source for special INTERGENIC */ 
        temp = score + (mat->geneopen) + (0) ;   
        if( temp > AlignWise_VSMALL_SPECIAL(mat,i,j,INTERGENIC) )    {  
          AlignWise_VSMALL_SPECIAL(mat,i,j,INTERGENIC) = temp;   
          }  




        /* Finished calculating state REVERSE_CODON */ 


        /* For state REVERSE_INTRON_0 */ 
        /* setting first movement to score */ 
        score = AlignWise_VSMALL_MATRIX(mat,i-0,j-8,REVERSE_CODON) + (mat->align->splice3_reverse[(j-7)]+mat->intronopen);   
        /* From state REVERSE_INTRON_0 to state REVERSE_INTRON_0 */ 
        temp = AlignWise_VSMALL_MATRIX(mat,i-0,j-1,REVERSE_INTRON_0) + 0;    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for REVERSE_INTRON_0 */ 
        /* Add any movement independant score and put away */ 
         AlignWise_VSMALL_MATRIX(mat,i,j,REVERSE_INTRON_0) = score;  


        /* Finished calculating state REVERSE_INTRON_0 */ 


        /* For state REVERSE_INTRON_1 */ 
        /* setting first movement to score */ 
        score = AlignWise_VSMALL_MATRIX(mat,i-0,j-9,REVERSE_CODON) + (mat->align->splice3_reverse[(j-7)]+mat->intronopen);   
        /* From state REVERSE_INTRON_1 to state REVERSE_INTRON_1 */ 
        temp = AlignWise_VSMALL_MATRIX(mat,i-0,j-1,REVERSE_INTRON_1) + 0;    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for REVERSE_INTRON_1 */ 
        /* Add any movement independant score and put away */ 
         AlignWise_VSMALL_MATRIX(mat,i,j,REVERSE_INTRON_1) = score;  


        /* Finished calculating state REVERSE_INTRON_1 */ 


        /* For state REVERSE_INTRON_2 */ 
        /* setting first movement to score */ 
        score = AlignWise_VSMALL_MATRIX(mat,i-0,j-10,REVERSE_CODON) + (mat->align->splice3_reverse[(j-7)]+mat->intronopen);  
        /* From state REVERSE_INTRON_2 to state REVERSE_INTRON_2 */ 
        temp = AlignWise_VSMALL_MATRIX(mat,i-0,j-1,REVERSE_INTRON_2) + 0;    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for REVERSE_INTRON_2 */ 
        /* Add any movement independant score and put away */ 
         AlignWise_VSMALL_MATRIX(mat,i,j,REVERSE_INTRON_2) = score;  


        /* Finished calculating state REVERSE_INTRON_2 */ 
        } /* end of for all query positions */ 




      /* Special state START has no special to special movements */ 


      /* Special state INTERGENIC has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = AlignWise_VSMALL_SPECIAL(mat,0,j,INTERGENIC);  


      /* Source START is a special source for INTERGENIC */ 
      temp = AlignWise_VSMALL_SPECIAL(mat,0,j - 1,START) + (0) + (0);    
      if( temp > score ) 
        score = temp;    


      /* Source INTERGENIC is a special source for INTERGENIC */ 
      temp = AlignWise_VSMALL_SPECIAL(mat,0,j - 1,INTERGENIC) + (0) + (0);   
      if( temp > score ) 
        score = temp;    


      /* Source CODON for state INTERGENIC is not special... already calculated */ 
      /* Source REVERSE_CODON for state INTERGENIC is not special... already calculated */ 
      /* Put back score... (now updated!) */ 
      AlignWise_VSMALL_SPECIAL(mat,0,j,INTERGENIC) = score;  
      /* Finished updating state INTERGENIC */ 




      /* Special state END has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = AlignWise_VSMALL_SPECIAL(mat,0,j,END); 


      /* Source INTERGENIC is a special source for END */ 
      temp = AlignWise_VSMALL_SPECIAL(mat,0,j - 1,INTERGENIC) + (0) + (0);   
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      AlignWise_VSMALL_SPECIAL(mat,0,j,END) = score; 
      /* Finished updating state END */ 


      if( bestscore < AlignWise_VSMALL_SPECIAL(mat,0,j,END) )    
        bestscore = AlignWise_VSMALL_SPECIAL(mat,0,j,END);   
      } /* end of for all target positions */ 


    mat = free_AlignWise(mat);   
    return bestscore;    
}    


/* Function:  PackAln_bestmemory_AlignWise(model,align,intronopen,geneopen,dpenv,dpri)
 *
 * Descrip:    This function chooses the best memory set-up for the alignment
 *             using calls to basematrix, and then implements either a large
 *             or small memory model.
 *
 *             It is the best function to use if you just want an alignment
 *
 *             If you want a label alignment, you will need
 *             /convert_PackAln_to_AlnBlock_AlignWise
 *
 *
 * Arg:             model [UNKN ] query data structure [AlignGeneModelFrame*]
 * Arg:             align [UNKN ] target data structure [AlignGeneModelScore*]
 * Arg:        intronopen [UNKN ] Resource [Score]
 * Arg:          geneopen [UNKN ] Resource [Score]
 * Arg:             dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:              dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * PackAln_bestmemory_AlignWise(AlignGeneModelFrame* model,AlignGeneModelScore* align ,Score intronopen,Score geneopen,DPEnvelope * dpenv,DPRunImpl * dpri) 
{
    long long total; 
    AlignWise * mat; 
    PackAln * out;   
    DebugMatrix * de;    
    DPRunImplMemory strategy;    
    assert(dpri);    


    total = model->len * align->len;     
    if( dpri->memory == DPIM_Default )   {  
      if( (total * 8 * sizeof(int)) > 1000*dpri->kbyte_size) {  
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
        if( (mat=allocate_Expl_AlignWise(model, align , intronopen, geneopen,dpri)) == NULL )    {  
          warn("Unable to allocate large AlignWise version");    
          return NULL;   
          }  
        calculate_dpenv_AlignWise(mat,dpenv);    
        out =  PackAln_read_Expl_AlignWise(mat); 
        }  
      else   {  
        mat = allocate_AlignWise_only(model, align , intronopen, geneopen);  
        calculate_shatter_AlignWise(mat,dpenv);  
        out = PackAln_read_Shatter_AlignWise(mat);   
        }  
      }  
    else {  
      if( strategy == DPIM_Linear )  {  
        /* use small implementation */ 
        if( (mat=allocate_Small_AlignWise(model, align , intronopen, geneopen)) == NULL )    {  
          warn("Unable to allocate small AlignWise version");    
          return NULL;   
          }  
        out = PackAln_calculate_Small_AlignWise(mat,dpenv);  
        }  
      else   {  
        /* use Large implementation */ 
        if( (mat=allocate_Expl_AlignWise(model, align , intronopen, geneopen,dpri)) == NULL )    {  
          warn("Unable to allocate large AlignWise version");    
          return NULL;   
          }  
        if( dpri->debug == TRUE) {  
          fatal("Asked for dydebug, but dynamite file not compiled with -g. Need to recompile dynamite source"); 
          }  
        else {  
          calculate_AlignWise(mat);  
          out =  PackAln_read_Expl_AlignWise(mat);   
          }  
        }  
      }  


    mat = free_AlignWise(mat);   
    return out;  
}    


/* Function:  allocate_AlignWise_only(model,align,intronopen,geneopen)
 *
 * Descrip:    This function only allocates the AlignWise structure
 *             checks types where possible and determines leni and lenj
 *             The basematrix area is delt with elsewhere
 *
 *
 * Arg:             model [UNKN ] query data structure [AlignGeneModelFrame*]
 * Arg:             align [UNKN ] target data structure [AlignGeneModelScore*]
 * Arg:        intronopen [UNKN ] Resource [Score]
 * Arg:          geneopen [UNKN ] Resource [Score]
 *
 * Return [UNKN ]  Undocumented return value [AlignWise *]
 *
 */
AlignWise * allocate_AlignWise_only(AlignGeneModelFrame* model,AlignGeneModelScore* align ,Score intronopen,Score geneopen) 
{
    AlignWise * out;     


    if((out= AlignWise_alloc()) == NULL) {  
      warn("Allocation of basic AlignWise structure failed..."); 
      return NULL;   
      }  


    out->model = model;  
    out->align = align;  
    out->intronopen = intronopen;    
    out->geneopen = geneopen;    
    out->leni = model->len;  
    out->lenj = align->len;  
    return out;  
}    


/* Function:  allocate_Expl_AlignWise(model,align,intronopen,geneopen,dpri)
 *
 * Descrip:    This function allocates the AlignWise structure
 *             and the basematrix area for explicit memory implementations
 *             It calls /allocate_AlignWise_only
 *
 *
 * Arg:             model [UNKN ] query data structure [AlignGeneModelFrame*]
 * Arg:             align [UNKN ] target data structure [AlignGeneModelScore*]
 * Arg:        intronopen [UNKN ] Resource [Score]
 * Arg:          geneopen [UNKN ] Resource [Score]
 * Arg:              dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [AlignWise *]
 *
 */
AlignWise * allocate_Expl_AlignWise(AlignGeneModelFrame* model,AlignGeneModelScore* align ,Score intronopen,Score geneopen,DPRunImpl * dpri) 
{
    AlignWise * out; 


    out = allocate_AlignWise_only(model, align , intronopen, geneopen);  
    if( out == NULL )    
      return NULL;   
    if( dpri->should_cache == TRUE ) {  
      if( dpri->cache != NULL )  {  
        if( dpri->cache->maxleni >= (out->lenj+10)*8 && dpri->cache->maxlenj >= (out->leni+0))   
          out->basematrix = hard_link_BaseMatrix(dpri->cache);   
        else 
          dpri->cache = free_BaseMatrix(dpri->cache);    
        }  
      }  
    if( out->basematrix == NULL )    {  
      if( (out->basematrix = BaseMatrix_alloc_matrix_and_specials((out->lenj+10)*8,(out->leni+0),3,out->lenj+10)) == NULL)   {  
        warn("Explicit matrix AlignWise cannot be allocated, (asking for %d by %d main cells)",out->leni,out->lenj); 
        free_AlignWise(out);     
        return NULL; 
        }  
      }  
    if( dpri->should_cache == TRUE && dpri->cache == NULL)   
      dpri->cache = hard_link_BaseMatrix(out->basematrix);   
    out->basematrix->type = BASEMATRIX_TYPE_EXPLICIT;    
    init_AlignWise(out);     
    return out;  
}    


/* Function:  init_AlignWise(mat)
 *
 * Descrip:    This function initates AlignWise matrix when in explicit mode
 *             Called in /allocate_Expl_AlignWise
 *
 *
 * Arg:        mat [UNKN ] AlignWise which contains explicit basematrix memory [AlignWise *]
 *
 */
void init_AlignWise(AlignWise * mat) 
{
    register int i;  
    register int j;  
    if( mat->basematrix->type != BASEMATRIX_TYPE_EXPLICIT)   {  
      warn("Cannot iniate matrix, is not an explicit memory type and you have assummed that");   
      return;    
      }  


    for(i= (-0);i<mat->model->len;i++)   {  
      for(j= (-10);j<11;j++) {  
        AlignWise_EXPL_MATRIX(mat,i,j,CODON) = NEGI; 
        AlignWise_EXPL_MATRIX(mat,i,j,INTRON_0) = NEGI;  
        AlignWise_EXPL_MATRIX(mat,i,j,INTRON_1) = NEGI;  
        AlignWise_EXPL_MATRIX(mat,i,j,INTRON_2) = NEGI;  
        AlignWise_EXPL_MATRIX(mat,i,j,REVERSE_CODON) = NEGI; 
        AlignWise_EXPL_MATRIX(mat,i,j,REVERSE_INTRON_0) = NEGI;  
        AlignWise_EXPL_MATRIX(mat,i,j,REVERSE_INTRON_1) = NEGI;  
        AlignWise_EXPL_MATRIX(mat,i,j,REVERSE_INTRON_2) = NEGI;  
        }  
      }  
    for(j= (-10);j<mat->align->len;j++)  {  
      for(i= (-0);i<1;i++)   {  
        AlignWise_EXPL_MATRIX(mat,i,j,CODON) = NEGI; 
        AlignWise_EXPL_MATRIX(mat,i,j,INTRON_0) = NEGI;  
        AlignWise_EXPL_MATRIX(mat,i,j,INTRON_1) = NEGI;  
        AlignWise_EXPL_MATRIX(mat,i,j,INTRON_2) = NEGI;  
        AlignWise_EXPL_MATRIX(mat,i,j,REVERSE_CODON) = NEGI; 
        AlignWise_EXPL_MATRIX(mat,i,j,REVERSE_INTRON_0) = NEGI;  
        AlignWise_EXPL_MATRIX(mat,i,j,REVERSE_INTRON_1) = NEGI;  
        AlignWise_EXPL_MATRIX(mat,i,j,REVERSE_INTRON_2) = NEGI;  
        }  
      AlignWise_EXPL_SPECIAL(mat,i,j,START) = 0; 
      AlignWise_EXPL_SPECIAL(mat,i,j,INTERGENIC) = NEGI; 
      AlignWise_EXPL_SPECIAL(mat,i,j,END) = NEGI;    
      }  
    return;  
}    


/* Function:  recalculate_PackAln_AlignWise(pal,mat)
 *
 * Descrip:    This function recalculates the PackAln structure produced by AlignWise
 *             For example, in linear space methods this is used to score them
 *
 *
 * Arg:        pal [UNKN ] Undocumented argument [PackAln *]
 * Arg:        mat [UNKN ] Undocumented argument [AlignWise *]
 *
 */
void recalculate_PackAln_AlignWise(PackAln * pal,AlignWise * mat) 
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
        case CODON :     
          if( offi == 0 && offj == 3 && prev->state == CODON )   {  
            pau->score = mat->align->forward_coding[j] + (0);    
            continue;    
            }  
          if( offj == 3 && prev->state == (INTERGENIC+8) )   {  
            pau->score = mat->align->forward_coding[j] + (0);    
            continue;    
            }  
          if( offi == 0 && offj == 6 && prev->state == INTRON_0 )    {  
            pau->score = (mat->align->splice3_forward[(j-3)]+mat->align->forward_coding[j]) + (0);   
            continue;    
            }  
          if( offi == 0 && offj == 5 && prev->state == INTRON_1 )    {  
            pau->score = mat->align->splice3_forward[(j-2)] + (0);   
            continue;    
            }  
          if( offi == 0 && offj == 4 && prev->state == INTRON_2 )    {  
            pau->score = mat->align->splice3_forward[(j-1)] + (0);   
            continue;    
            }  
          warn("In recaluclating PackAln with state CODON, from [%d,%d,%d], got a bad source state. Error!",offi,offj,prev->state);  
          break; 
        case INTRON_0 :  
          if( offi == 0 && offj == 8 && prev->state == CODON )   {  
            pau->score = (mat->align->splice5_forward[(j-7)]+mat->intronopen) + (0);     
            continue;    
            }  
          if( offi == 0 && offj == 1 && prev->state == INTRON_0 )    {  
            pau->score = 0 + (0);    
            continue;    
            }  
          warn("In recaluclating PackAln with state INTRON_0, from [%d,%d,%d], got a bad source state. Error!",offi,offj,prev->state);   
          break; 
        case INTRON_1 :  
          if( offi == 0 && offj == 9 && prev->state == CODON )   {  
            pau->score = (mat->align->splice5_forward[(j-7)]+mat->intronopen) + (0);     
            continue;    
            }  
          if( offi == 0 && offj == 1 && prev->state == INTRON_1 )    {  
            pau->score = 0 + (0);    
            continue;    
            }  
          warn("In recaluclating PackAln with state INTRON_1, from [%d,%d,%d], got a bad source state. Error!",offi,offj,prev->state);   
          break; 
        case INTRON_2 :  
          if( offi == 0 && offj == 10 && prev->state == CODON )  {  
            pau->score = (mat->align->splice5_forward[(j-7)]+mat->intronopen) + (0);     
            continue;    
            }  
          if( offi == 0 && offj == 1 && prev->state == INTRON_2 )    {  
            pau->score = 0 + (0);    
            continue;    
            }  
          warn("In recaluclating PackAln with state INTRON_2, from [%d,%d,%d], got a bad source state. Error!",offi,offj,prev->state);   
          break; 
        case REVERSE_CODON :     
          if( offi == 0 && offj == 3 && prev->state == REVERSE_CODON )   {  
            pau->score = mat->align->reverse_coding[j] + (0);    
            continue;    
            }  
          if( offj == 3 && prev->state == (INTERGENIC+8) )   {  
            pau->score = mat->align->reverse_coding[j] + (0);    
            continue;    
            }  
          if( offi == 0 && offj == 6 && prev->state == REVERSE_INTRON_0 )    {  
            pau->score = mat->align->splice5_reverse[(j-3)] + (0);   
            continue;    
            }  
          if( offi == 0 && offj == 5 && prev->state == REVERSE_INTRON_1 )    {  
            pau->score = mat->align->splice5_reverse[(j-2)] + (0);   
            continue;    
            }  
          if( offi == 0 && offj == 4 && prev->state == REVERSE_INTRON_2 )    {  
            pau->score = mat->align->splice5_reverse[(j-1)] + (0);   
            continue;    
            }  
          warn("In recaluclating PackAln with state REVERSE_CODON, from [%d,%d,%d], got a bad source state. Error!",offi,offj,prev->state);  
          break; 
        case REVERSE_INTRON_0 :  
          if( offi == 0 && offj == 8 && prev->state == REVERSE_CODON )   {  
            pau->score = (mat->align->splice3_reverse[(j-7)]+mat->intronopen) + (0);     
            continue;    
            }  
          if( offi == 0 && offj == 1 && prev->state == REVERSE_INTRON_0 )    {  
            pau->score = 0 + (0);    
            continue;    
            }  
          warn("In recaluclating PackAln with state REVERSE_INTRON_0, from [%d,%d,%d], got a bad source state. Error!",offi,offj,prev->state);   
          break; 
        case REVERSE_INTRON_1 :  
          if( offi == 0 && offj == 9 && prev->state == REVERSE_CODON )   {  
            pau->score = (mat->align->splice3_reverse[(j-7)]+mat->intronopen) + (0);     
            continue;    
            }  
          if( offi == 0 && offj == 1 && prev->state == REVERSE_INTRON_1 )    {  
            pau->score = 0 + (0);    
            continue;    
            }  
          warn("In recaluclating PackAln with state REVERSE_INTRON_1, from [%d,%d,%d], got a bad source state. Error!",offi,offj,prev->state);   
          break; 
        case REVERSE_INTRON_2 :  
          if( offi == 0 && offj == 10 && prev->state == REVERSE_CODON )  {  
            pau->score = (mat->align->splice3_reverse[(j-7)]+mat->intronopen) + (0);     
            continue;    
            }  
          if( offi == 0 && offj == 1 && prev->state == REVERSE_INTRON_2 )    {  
            pau->score = 0 + (0);    
            continue;    
            }  
          warn("In recaluclating PackAln with state REVERSE_INTRON_2, from [%d,%d,%d], got a bad source state. Error!",offi,offj,prev->state);   
          break; 
        case (START+8) :     
          warn("In recaluclating PackAln with state START, got a bad source state. Error!"); 
          break; 
        case (INTERGENIC+8) :    
          if( offj == 1 && prev->state == (START+8) )    {  
            pau->score = 0 + (0);    
            continue;    
            }  
          if( offj == 1 && prev->state == (INTERGENIC+8) )   {  
            pau->score = 0 + (0);    
            continue;    
            }  
          if( offj == 0 && prev->state == CODON )    {  
            /* i here comes from the previous state ;) - not the real one */ 
            i = prev->i; 
            pau->score = mat->geneopen + (0);    
            continue;    
            }  
          if( offj == 0 && prev->state == REVERSE_CODON )    {  
            /* i here comes from the previous state ;) - not the real one */ 
            i = prev->i; 
            pau->score = mat->geneopen + (0);    
            continue;    
            }  
          warn("In recaluclating PackAln with state INTERGENIC, got a bad source state. Error!");    
          break; 
        case (END+8) :   
          if( offj == 1 && prev->state == (INTERGENIC+8) )   {  
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
#define AlignWise_HIDDEN_MATRIX(thismatrix,i,j,state) (thismatrix->basematrix->matrix[(j-hiddenj+10)][(i+0)*8+state])    
#define AlignWise_DC_SHADOW_MATRIX(thismatrix,i,j,state) (thismatrix->basematrix->matrix[((j+11)*8) % 88][(i+0)*8+state])    
#define AlignWise_HIDDEN_SPECIAL(thismatrix,i,j,state) (thismatrix->basematrix->specmatrix[state][(j+10)])   
#define AlignWise_DC_SHADOW_SPECIAL(thismatrix,i,j,state) (thismatrix->basematrix->specmatrix[state*8][(j+10)])  
#define AlignWise_DC_SHADOW_MATRIX_SP(thismatrix,i,j,state,shadow) (thismatrix->basematrix->matrix[((((j+11)*8)+(shadow+1)) % 88)][(i+0)*8 + state]) 
#define AlignWise_DC_SHADOW_SPECIAL_SP(thismatrix,i,j,state,shadow) (thismatrix->basematrix->specmatrix[state*8 +shadow+1][(j+10)])  
#define AlignWise_DC_OPT_SHADOW_MATRIX(thismatrix,i,j,state) (score_pointers[(((j+10)% 10) * (leni+1) * 8) + ((i+0) * 8) + (state)]) 
#define AlignWise_DC_OPT_SHADOW_MATRIX_SP(thismatrix,i,j,state,shadow) (shadow_pointers[(((j+10)% 10) * (leni+1) * 64) + ((i+0) * 64) + (state * 8) + shadow+1]) 
#define AlignWise_DC_OPT_SHADOW_SPECIAL(thismatrix,i,j,state) (thismatrix->basematrix->specmatrix[state*8][(j+10)])  
/* Function:  allocate_Small_AlignWise(model,align,intronopen,geneopen)
 *
 * Descrip:    This function allocates the AlignWise structure
 *             and the basematrix area for a small memory implementations
 *             It calls /allocate_AlignWise_only
 *
 *
 * Arg:             model [UNKN ] query data structure [AlignGeneModelFrame*]
 * Arg:             align [UNKN ] target data structure [AlignGeneModelScore*]
 * Arg:        intronopen [UNKN ] Resource [Score]
 * Arg:          geneopen [UNKN ] Resource [Score]
 *
 * Return [UNKN ]  Undocumented return value [AlignWise *]
 *
 */
#define AlignWise_DC_OPT_SHADOW_SPECIAL_SP(thismatrix,i,j,state,shadow) (thismatrix->basematrix->specmatrix[state*8 +shadow+1][(j+10)])  
AlignWise * allocate_Small_AlignWise(AlignGeneModelFrame* model,AlignGeneModelScore* align ,Score intronopen,Score geneopen) 
{
    AlignWise * out; 


    out = allocate_AlignWise_only(model, align , intronopen, geneopen);  
    if( out == NULL )    
      return NULL;   
    out->basematrix = BaseMatrix_alloc_matrix_and_specials(88,(out->leni + 0) * 8,24,out->lenj+10);  
    if(out == NULL)  {  
      warn("Small shadow matrix AlignWise cannot be allocated, (asking for 11 by %d main cells)",out->leni+1);   
      free_AlignWise(out);   
      return NULL;   
      }  
    out->basematrix->type = BASEMATRIX_TYPE_SHADOW;  
    return out;  
}    


/* Function:  PackAln_calculate_Small_AlignWise(mat,dpenv)
 *
 * Descrip:    This function calculates an alignment for AlignWise structure in linear space
 *             If you want only the start/end points
 *             use /AlnRangeSet_calculate_Small_AlignWise 
 *
 *             The function basically
 *               finds start/end points 
 *               foreach start/end point 
 *                 calls /full_dc_AlignWise 
 *
 *
 * Arg:          mat [UNKN ] Undocumented argument [AlignWise *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * PackAln_calculate_Small_AlignWise(AlignWise * mat,DPEnvelope * dpenv) 
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
      warn("Could not calculate packaln small for AlignWise due to wrong type of matrix");   
      return NULL;   
      }  


    out = PackAln_alloc_std();   


    start_reporting("Find start end points: ");  
    dc_optimised_start_end_calc_AlignWise(mat,dpenv);    
    score = start_end_find_end_AlignWise(mat,&endj); 
    out->score = score;  
    stopstate = END;
    
    /* Special to specials: have to eat up in strip and then drop back to full_dc for intervening bits */ 
    log_full_error(REPORT,0,"End at %d Score %d",endj,score);    
    stop_reporting();    
    for(;;)  { /*while there are more special bits to recover*/ 
      start_reporting("Special cell aln end   %d:",endj);    
      if( read_special_strip_AlignWise(mat,0,endj,stopstate,&endj,&startstate,out) == FALSE )    {  
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
      starti = AlignWise_DC_SHADOW_SPECIAL_SP(mat,0,endj,temp,0);    
      startj = AlignWise_DC_SHADOW_SPECIAL_SP(mat,0,endj,temp,1);    
      startstate = AlignWise_DC_SHADOW_SPECIAL_SP(mat,0,endj,temp,2);    
      stopi = AlignWise_DC_SHADOW_SPECIAL_SP(mat,0,endj,temp,3); 
      stopj = AlignWise_DC_SHADOW_SPECIAL_SP(mat,0,endj,temp,4); 
      stopstate = AlignWise_DC_SHADOW_SPECIAL_SP(mat,0,endj,temp,5); 


      /* Get out the score of this block. V. important! */ 
      temp = AlignWise_DC_SHADOW_SPECIAL_SP(mat,0,endj,temp,6);  
      totalj = stopj - startj;   
      donej  = 0;    
      start_reporting("Main matrix  aln [%d,%d]:",startj,stopj);     
      if(full_dc_AlignWise(mat,starti,startj,startstate,stopi,stopj,stopstate,out,&donej,totalj,dpenv) == FALSE) {  
        warn("In the alignment AlignWise [%d,%d][%d,%d], got a problem. Please report bug ... giving you back a partial alignment",starti,startj,stopi,stopj);   
        return out;  
        }  


      /* now have to figure out which special we came from... yikes */ 
      max_matrix_to_special_AlignWise(mat,starti,startj,startstate,temp,&stopi,&stopj,&stopstate,&temp,NULL);    
      if( stopi == AlignWise_READ_OFF_ERROR) {  
        warn("In AlignWise read off ending at %d ... got a bad matrix to special read off... returning partial alignment",startj);   
        invert_PackAln(out); 
        recalculate_PackAln_AlignWise(out,mat);  
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
    recalculate_PackAln_AlignWise(out,mat);  
    return out;  


}    


/* Function:  AlnRangeSet_calculate_Small_AlignWise(mat)
 *
 * Descrip:    This function calculates an alignment for AlignWise structure in linear space
 *             If you want the full alignment, use /PackAln_calculate_Small_AlignWise 
 *             If you have already got the full alignment, but want the range set, use /AlnRangeSet_from_PackAln_AlignWise
 *             If you have got the small matrix but not the alignment, use /AlnRangeSet_from_AlignWise 
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [AlignWise *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * AlnRangeSet_calculate_Small_AlignWise(AlignWise * mat) 
{
    AlnRangeSet * out;   


    start_reporting("Find start end points: ");  
    dc_optimised_start_end_calc_AlignWise(mat,NULL); 
    log_full_error(REPORT,0,"Calculated");   


    out = AlnRangeSet_from_AlignWise(mat);   
    return out;  
}    


/* Function:  AlnRangeSet_from_AlignWise(mat)
 *
 * Descrip:    This function reads off a start/end structure
 *             for AlignWise structure in linear space
 *             If you want the full alignment use
 *             /PackAln_calculate_Small_AlignWise 
 *             If you have not calculated the matrix use
 *             /AlnRange_calculate_Small_AlignWise
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [AlignWise *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * AlnRangeSet_from_AlignWise(AlignWise * mat) 
{
    AlnRangeSet * out;   
    AlnRange * temp; 
    int jpos;    
    int state;   


    if( mat->basematrix->type != BASEMATRIX_TYPE_SHADOW) {  
      warn("Bad error! - non shadow matrix type in AlnRangeSet_from_AlignWise"); 
      return NULL;   
      }  


    out = AlnRangeSet_alloc_std();   
    /* Find the end position */ 
    out->score = start_end_find_end_AlignWise(mat,&jpos);    
    state = END; 


    while( (temp = AlnRange_build_AlignWise(mat,jpos,state,&jpos,&state)) != NULL)   
      add_AlnRangeSet(out,temp); 
    return out;  
}    


/* Function:  AlnRange_build_AlignWise(mat,stopj,stopspecstate,startj,startspecstate)
 *
 * Descrip:    This function calculates a single start/end set in linear space
 *             Really a sub-routine for /AlnRangeSet_from_PackAln_AlignWise
 *
 *
 * Arg:                   mat [UNKN ] Undocumented argument [AlignWise *]
 * Arg:                 stopj [UNKN ] Undocumented argument [int]
 * Arg:         stopspecstate [UNKN ] Undocumented argument [int]
 * Arg:                startj [UNKN ] Undocumented argument [int *]
 * Arg:        startspecstate [UNKN ] Undocumented argument [int *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRange *]
 *
 */
AlnRange * AlnRange_build_AlignWise(AlignWise * mat,int stopj,int stopspecstate,int * startj,int * startspecstate) 
{
    AlnRange * out;  
    int jpos;    
    int state;   


    if( mat->basematrix->type != BASEMATRIX_TYPE_SHADOW) {  
      warn("Bad error! - non shadow matrix type in AlnRangeSet_from_AlignWise"); 
      return NULL;   
      }  


    /* Assumme that we have specials (we should!). Read back along the specials till we have the finish point */ 
    if( read_special_strip_AlignWise(mat,0,stopj,stopspecstate,&jpos,&state,NULL) == FALSE)  {  
      warn("In AlnRanger_build_AlignWise alignment ending at %d, unable to read back specials. Will (evenutally) return a partial range set... BEWARE!",stopj);  
      return NULL;   
      }  
    if( state == START || jpos <= 0) 
      return NULL;   


    out = AlnRange_alloc();  


    out->starti = AlignWise_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,0);    
    out->startj = AlignWise_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,1);    
    out->startstate = AlignWise_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,2);    
    out->stopi = AlignWise_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,3); 
    out->stopj = AlignWise_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,4); 
    out->stopstate = AlignWise_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,5); 
    out->startscore = AlignWise_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,6);    
    out->stopscore = AlignWise_DC_SHADOW_SPECIAL(mat,0,jpos,state);  


    /* Now, we have to figure out where this state came from in the specials */ 
    max_matrix_to_special_AlignWise(mat,out->starti,out->startj,out->startstate,out->startscore,&jpos,startj,startspecstate,&state,NULL);    
    if( jpos == AlignWise_READ_OFF_ERROR)    {  
      warn("In AlnRange_build_AlignWise alignment ending at %d, with aln range between %d-%d in j, unable to find source special, returning this range, but this could get tricky!",stopj,out->startj,out->stopj);   
      return out;    
      }  


    /* Put in the correct score for startstate, from the special */ 
    out->startscore = AlignWise_DC_SHADOW_SPECIAL(mat,0,*startj,*startspecstate);    
    /* The correct j coords have been put into startj, startspecstate... so just return out */ 
    return out;  
}    


/* Function:  read_hidden_AlignWise(mat,starti,startj,startstate,stopi,stopj,stopstate,out)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [AlignWise *]
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
boolean read_hidden_AlignWise(AlignWise * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out) 
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


      max_hidden_AlignWise(mat,startj,i,j,state,isspecial,&i,&j,&state,&isspecial,&cellscore);   


      if( i == AlignWise_READ_OFF_ERROR) {  
        warn("In AlignWise hidden read off, between %d:%d,%d:%d - at got bad read off. Problem!",starti,startj,stopi,stopj); 
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
        warn("In AlignWise hidden read off, between %d:%d,%d:%d - hit start cell, but not in start state. Can't be good!.",starti,startj,stopi,stopj);   
        return FALSE;    
        }  
      }  
    warn("In AlignWise hidden read off, between %d:%d,%d:%d - gone past start cell (now in %d,%d,%d), can't be good news!.",starti,startj,stopi,stopj,i,j,state);    
    return FALSE;    
}    


/* Function:  max_hidden_AlignWise(mat,hiddenj,i,j,state,isspecial,reti,retj,retstate,retspecial,cellscore)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [AlignWise *]
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
int max_hidden_AlignWise(AlignWise * mat,int hiddenj,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore) 
{
    register int temp;   
    register int cscore; 


    *reti = (*retj) = (*retstate) = AlignWise_READ_OFF_ERROR;    


    if( i < 0 || j < 0 || i > mat->model->len || j > mat->align->len)    {  
      warn("In AlignWise matrix special read off - out of bounds on matrix [i,j is %d,%d state %d in standard matrix]",i,j,state);   
      return -1; 
      }  


    /* Then you have to select the correct switch statement to figure out the readoff      */ 
    /* Somewhat odd - reverse the order of calculation and return as soon as it is correct */ 
    cscore = AlignWise_HIDDEN_MATRIX(mat,i,j,state); 
    switch(state)    { /*Switch state */ 
      case CODON :   
        temp = cscore - (mat->align->splice3_forward[(j-1)]) -  (0); 
        if( temp == AlignWise_HIDDEN_MATRIX(mat,i - 0,j - 4,INTRON_2) )  {  
          *reti = i - 0; 
          *retj = j - 4; 
          *retstate = INTRON_2;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - AlignWise_HIDDEN_MATRIX(mat,i-0,j-4,INTRON_2); 
            }  
          return AlignWise_HIDDEN_MATRIX(mat,i - 0,j - 4,INTRON_2);  
          }  
        temp = cscore - (mat->align->splice3_forward[(j-2)]) -  (0); 
        if( temp == AlignWise_HIDDEN_MATRIX(mat,i - 0,j - 5,INTRON_1) )  {  
          *reti = i - 0; 
          *retj = j - 5; 
          *retstate = INTRON_1;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - AlignWise_HIDDEN_MATRIX(mat,i-0,j-5,INTRON_1); 
            }  
          return AlignWise_HIDDEN_MATRIX(mat,i - 0,j - 5,INTRON_1);  
          }  
        temp = cscore - ((mat->align->splice3_forward[(j-3)]+mat->align->forward_coding[j])) -  (0); 
        if( temp == AlignWise_HIDDEN_MATRIX(mat,i - 0,j - 6,INTRON_0) )  {  
          *reti = i - 0; 
          *retj = j - 6; 
          *retstate = INTRON_0;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - AlignWise_HIDDEN_MATRIX(mat,i-0,j-6,INTRON_0); 
            }  
          return AlignWise_HIDDEN_MATRIX(mat,i - 0,j - 6,INTRON_0);  
          }  
        /* Not allowing special sources.. skipping INTERGENIC */ 
        temp = cscore - (mat->align->forward_coding[j]) -  (0);  
        if( temp == AlignWise_HIDDEN_MATRIX(mat,i - 0,j - 3,CODON) ) {  
          *reti = i - 0; 
          *retj = j - 3; 
          *retstate = CODON; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - AlignWise_HIDDEN_MATRIX(mat,i-0,j-3,CODON);    
            }  
          return AlignWise_HIDDEN_MATRIX(mat,i - 0,j - 3,CODON);     
          }  
        warn("Major problem (!) - in AlignWise read off, position %d,%d state %d no source found!",i,j,state);   
        return (-1); 
      case INTRON_0 :    
        temp = cscore - (0) -  (0);  
        if( temp == AlignWise_HIDDEN_MATRIX(mat,i - 0,j - 1,INTRON_0) )  {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = INTRON_0;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - AlignWise_HIDDEN_MATRIX(mat,i-0,j-1,INTRON_0); 
            }  
          return AlignWise_HIDDEN_MATRIX(mat,i - 0,j - 1,INTRON_0);  
          }  
        temp = cscore - ((mat->align->splice5_forward[(j-7)]+mat->intronopen)) -  (0);   
        if( temp == AlignWise_HIDDEN_MATRIX(mat,i - 0,j - 8,CODON) ) {  
          *reti = i - 0; 
          *retj = j - 8; 
          *retstate = CODON; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - AlignWise_HIDDEN_MATRIX(mat,i-0,j-8,CODON);    
            }  
          return AlignWise_HIDDEN_MATRIX(mat,i - 0,j - 8,CODON);     
          }  
        warn("Major problem (!) - in AlignWise read off, position %d,%d state %d no source found!",i,j,state);   
        return (-1); 
      case INTRON_1 :    
        temp = cscore - (0) -  (0);  
        if( temp == AlignWise_HIDDEN_MATRIX(mat,i - 0,j - 1,INTRON_1) )  {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = INTRON_1;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - AlignWise_HIDDEN_MATRIX(mat,i-0,j-1,INTRON_1); 
            }  
          return AlignWise_HIDDEN_MATRIX(mat,i - 0,j - 1,INTRON_1);  
          }  
        temp = cscore - ((mat->align->splice5_forward[(j-7)]+mat->intronopen)) -  (0);   
        if( temp == AlignWise_HIDDEN_MATRIX(mat,i - 0,j - 9,CODON) ) {  
          *reti = i - 0; 
          *retj = j - 9; 
          *retstate = CODON; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - AlignWise_HIDDEN_MATRIX(mat,i-0,j-9,CODON);    
            }  
          return AlignWise_HIDDEN_MATRIX(mat,i - 0,j - 9,CODON);     
          }  
        warn("Major problem (!) - in AlignWise read off, position %d,%d state %d no source found!",i,j,state);   
        return (-1); 
      case INTRON_2 :    
        temp = cscore - (0) -  (0);  
        if( temp == AlignWise_HIDDEN_MATRIX(mat,i - 0,j - 1,INTRON_2) )  {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = INTRON_2;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - AlignWise_HIDDEN_MATRIX(mat,i-0,j-1,INTRON_2); 
            }  
          return AlignWise_HIDDEN_MATRIX(mat,i - 0,j - 1,INTRON_2);  
          }  
        temp = cscore - ((mat->align->splice5_forward[(j-7)]+mat->intronopen)) -  (0);   
        if( temp == AlignWise_HIDDEN_MATRIX(mat,i - 0,j - 10,CODON) )    {  
          *reti = i - 0; 
          *retj = j - 10;    
          *retstate = CODON; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - AlignWise_HIDDEN_MATRIX(mat,i-0,j-10,CODON);   
            }  
          return AlignWise_HIDDEN_MATRIX(mat,i - 0,j - 10,CODON);    
          }  
        warn("Major problem (!) - in AlignWise read off, position %d,%d state %d no source found!",i,j,state);   
        return (-1); 
      case REVERSE_CODON :   
        temp = cscore - (mat->align->splice5_reverse[(j-1)]) -  (0); 
        if( temp == AlignWise_HIDDEN_MATRIX(mat,i - 0,j - 4,REVERSE_INTRON_2) )  {  
          *reti = i - 0; 
          *retj = j - 4; 
          *retstate = REVERSE_INTRON_2;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - AlignWise_HIDDEN_MATRIX(mat,i-0,j-4,REVERSE_INTRON_2); 
            }  
          return AlignWise_HIDDEN_MATRIX(mat,i - 0,j - 4,REVERSE_INTRON_2);  
          }  
        temp = cscore - (mat->align->splice5_reverse[(j-2)]) -  (0); 
        if( temp == AlignWise_HIDDEN_MATRIX(mat,i - 0,j - 5,REVERSE_INTRON_1) )  {  
          *reti = i - 0; 
          *retj = j - 5; 
          *retstate = REVERSE_INTRON_1;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - AlignWise_HIDDEN_MATRIX(mat,i-0,j-5,REVERSE_INTRON_1); 
            }  
          return AlignWise_HIDDEN_MATRIX(mat,i - 0,j - 5,REVERSE_INTRON_1);  
          }  
        temp = cscore - (mat->align->splice5_reverse[(j-3)]) -  (0); 
        if( temp == AlignWise_HIDDEN_MATRIX(mat,i - 0,j - 6,REVERSE_INTRON_0) )  {  
          *reti = i - 0; 
          *retj = j - 6; 
          *retstate = REVERSE_INTRON_0;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - AlignWise_HIDDEN_MATRIX(mat,i-0,j-6,REVERSE_INTRON_0); 
            }  
          return AlignWise_HIDDEN_MATRIX(mat,i - 0,j - 6,REVERSE_INTRON_0);  
          }  
        /* Not allowing special sources.. skipping INTERGENIC */ 
        temp = cscore - (mat->align->reverse_coding[j]) -  (0);  
        if( temp == AlignWise_HIDDEN_MATRIX(mat,i - 0,j - 3,REVERSE_CODON) ) {  
          *reti = i - 0; 
          *retj = j - 3; 
          *retstate = REVERSE_CODON; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - AlignWise_HIDDEN_MATRIX(mat,i-0,j-3,REVERSE_CODON);    
            }  
          return AlignWise_HIDDEN_MATRIX(mat,i - 0,j - 3,REVERSE_CODON);     
          }  
        warn("Major problem (!) - in AlignWise read off, position %d,%d state %d no source found!",i,j,state);   
        return (-1); 
      case REVERSE_INTRON_0 :    
        temp = cscore - (0) -  (0);  
        if( temp == AlignWise_HIDDEN_MATRIX(mat,i - 0,j - 1,REVERSE_INTRON_0) )  {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = REVERSE_INTRON_0;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - AlignWise_HIDDEN_MATRIX(mat,i-0,j-1,REVERSE_INTRON_0); 
            }  
          return AlignWise_HIDDEN_MATRIX(mat,i - 0,j - 1,REVERSE_INTRON_0);  
          }  
        temp = cscore - ((mat->align->splice3_reverse[(j-7)]+mat->intronopen)) -  (0);   
        if( temp == AlignWise_HIDDEN_MATRIX(mat,i - 0,j - 8,REVERSE_CODON) ) {  
          *reti = i - 0; 
          *retj = j - 8; 
          *retstate = REVERSE_CODON; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - AlignWise_HIDDEN_MATRIX(mat,i-0,j-8,REVERSE_CODON);    
            }  
          return AlignWise_HIDDEN_MATRIX(mat,i - 0,j - 8,REVERSE_CODON);     
          }  
        warn("Major problem (!) - in AlignWise read off, position %d,%d state %d no source found!",i,j,state);   
        return (-1); 
      case REVERSE_INTRON_1 :    
        temp = cscore - (0) -  (0);  
        if( temp == AlignWise_HIDDEN_MATRIX(mat,i - 0,j - 1,REVERSE_INTRON_1) )  {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = REVERSE_INTRON_1;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - AlignWise_HIDDEN_MATRIX(mat,i-0,j-1,REVERSE_INTRON_1); 
            }  
          return AlignWise_HIDDEN_MATRIX(mat,i - 0,j - 1,REVERSE_INTRON_1);  
          }  
        temp = cscore - ((mat->align->splice3_reverse[(j-7)]+mat->intronopen)) -  (0);   
        if( temp == AlignWise_HIDDEN_MATRIX(mat,i - 0,j - 9,REVERSE_CODON) ) {  
          *reti = i - 0; 
          *retj = j - 9; 
          *retstate = REVERSE_CODON; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - AlignWise_HIDDEN_MATRIX(mat,i-0,j-9,REVERSE_CODON);    
            }  
          return AlignWise_HIDDEN_MATRIX(mat,i - 0,j - 9,REVERSE_CODON);     
          }  
        warn("Major problem (!) - in AlignWise read off, position %d,%d state %d no source found!",i,j,state);   
        return (-1); 
      case REVERSE_INTRON_2 :    
        temp = cscore - (0) -  (0);  
        if( temp == AlignWise_HIDDEN_MATRIX(mat,i - 0,j - 1,REVERSE_INTRON_2) )  {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = REVERSE_INTRON_2;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - AlignWise_HIDDEN_MATRIX(mat,i-0,j-1,REVERSE_INTRON_2); 
            }  
          return AlignWise_HIDDEN_MATRIX(mat,i - 0,j - 1,REVERSE_INTRON_2);  
          }  
        temp = cscore - ((mat->align->splice3_reverse[(j-7)]+mat->intronopen)) -  (0);   
        if( temp == AlignWise_HIDDEN_MATRIX(mat,i - 0,j - 10,REVERSE_CODON) )    {  
          *reti = i - 0; 
          *retj = j - 10;    
          *retstate = REVERSE_CODON; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - AlignWise_HIDDEN_MATRIX(mat,i-0,j-10,REVERSE_CODON);   
            }  
          return AlignWise_HIDDEN_MATRIX(mat,i - 0,j - 10,REVERSE_CODON);    
          }  
        warn("Major problem (!) - in AlignWise read off, position %d,%d state %d no source found!",i,j,state);   
        return (-1); 
      default:   
        warn("Major problem (!) - in AlignWise read off, position %d,%d state %d no source found!",i,j,state);   
        return (-1); 
      } /* end of Switch state  */ 
}    


/* Function:  read_special_strip_AlignWise(mat,stopi,stopj,stopstate,startj,startstate,out)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [AlignWise *]
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
boolean read_special_strip_AlignWise(AlignWise * mat,int stopi,int stopj,int stopstate,int * startj,int * startstate,PackAln * out) 
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
    while( j > AlignWise_DC_SHADOW_SPECIAL_SP(mat,i,j,state,4) && state != START)    { /*while more specials to eat up*/ 
      /* Put away current state, if we should */ 
      if(out != NULL)    {  
        pau = PackAlnUnit_alloc();  /* Should deal with memory overflow */ 
        pau->i = i;  
        pau->j = j;  
        pau->state =  state + 8; 
        add_PackAln(out,pau);    
        }  


      max_special_strip_AlignWise(mat,i,j,state,isspecial,&i,&j,&state,&isspecial,&cellscore);   
      if( i == AlignWise_READ_OFF_ERROR) {  
        warn("In special strip read AlignWise, got a bad read off error. Sorry!");   
        return FALSE;    
        }  
      } /* end of while more specials to eat up */ 


    /* check to see we have not gone too far! */ 
    if( state != START && j < AlignWise_DC_SHADOW_SPECIAL_SP(mat,i,j,state,4))   {  
      warn("In special strip read AlignWise, at special [%d] state [%d] overshot!",j,state); 
      return FALSE;  
      }  
    /* Put away last state */ 
    if(out != NULL)  {  
      pau = PackAlnUnit_alloc();/* Should deal with memory overflow */ 
      pau->i = i;    
      pau->j = j;    
      pau->state =  state + 8;   
      add_PackAln(out,pau);  
      }  


    /* Put away where we are in startj and startstate */ 
    *startj = j; 
    *startstate = state; 
    return TRUE; 
}    


/* Function:  max_special_strip_AlignWise(mat,i,j,state,isspecial,reti,retj,retstate,retspecial,cellscore)
 *
 * Descrip:    A pretty intense internal function. Deals with read-off only in specials
 *
 *
 * Arg:               mat [UNKN ] Undocumented argument [AlignWise *]
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
int max_special_strip_AlignWise(AlignWise * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore) 
{
    int temp;    
    int cscore;  


    *reti = (*retj) = (*retstate) = AlignWise_READ_OFF_ERROR;    
    if( isspecial == FALSE ) {  
      warn("In special strip max function for AlignWise, got a non special start point. Problem! (bad!)");   
      return (-1);   
      }  


    if( j < 0 || j > mat->align->len)    {  
      warn("In AlignWise matrix special read off - out of bounds on matrix [j is %d in special]",j); 
      return -1; 
      }  


    cscore = AlignWise_DC_SHADOW_SPECIAL(mat,i,j,state); 
    switch(state)    { /*switch on special states*/ 
      case START :   
      case INTERGENIC :  
        /* Source REVERSE_CODON is not a special */ 
        /* Source CODON is not a special */ 
        /* source INTERGENIC is a special */ 
        temp = cscore - (0) - (0);   
        if( temp == AlignWise_DC_SHADOW_SPECIAL(mat,i - 0,j - 1,INTERGENIC) )    {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = INTERGENIC;    
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - AlignWise_DC_SHADOW_SPECIAL(mat,i-0,j-1,INTERGENIC);   
            }  
          return AlignWise_DC_SHADOW_MATRIX(mat,i - 0,j - 1,INTERGENIC) ;    
          }  
        /* source START is a special */ 
        temp = cscore - (0) - (0);   
        if( temp == AlignWise_DC_SHADOW_SPECIAL(mat,i - 0,j - 1,START) ) {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = START; 
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - AlignWise_DC_SHADOW_SPECIAL(mat,i-0,j-1,START);    
            }  
          return AlignWise_DC_SHADOW_MATRIX(mat,i - 0,j - 1,START) ;     
          }  
      case END :     
        /* source INTERGENIC is a special */ 
        temp = cscore - (0) - (0);   
        if( temp == AlignWise_DC_SHADOW_SPECIAL(mat,i - 0,j - 1,INTERGENIC) )    {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = INTERGENIC;    
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - AlignWise_DC_SHADOW_SPECIAL(mat,i-0,j-1,INTERGENIC);   
            }  
          return AlignWise_DC_SHADOW_MATRIX(mat,i - 0,j - 1,INTERGENIC) ;    
          }  
      default:   
        warn("Major problem (!) - in AlignWise special strip read off, position %d,%d state %d no source found  dropped into default on source switch!",i,j,state);  
        return (-1); 
      } /* end of switch on special states */ 
}    


/* Function:  max_matrix_to_special_AlignWise(mat,i,j,state,cscore,reti,retj,retstate,retspecial,cellscore)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [AlignWise *]
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
int max_matrix_to_special_AlignWise(AlignWise * mat,int i,int j,int state,int cscore,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore) 
{
    int temp;    
    *reti = (*retj) = (*retstate) = AlignWise_READ_OFF_ERROR;    


    if( j < 0 || j > mat->lenj)  {  
      warn("In AlignWise matrix to special read off - out of bounds on matrix [j is %d in special]",j);  
      return -1; 
      }  


    switch(state)    { /*Switch state */ 
      case CODON :   
        /* Source INTRON_2 is not a special, should not get here! */ 
        /* Source INTRON_1 is not a special, should not get here! */ 
        /* Source INTRON_0 is not a special, should not get here! */ 
        temp = cscore - (mat->align->forward_coding[j]) -  (0);  
        if( temp == AlignWise_DC_SHADOW_SPECIAL(mat,i - 0,j - 3,INTERGENIC) )    {  
          *reti = i - 0; 
          *retj = j - 3; 
          *retstate = INTERGENIC;    
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - AlignWise_DC_SHADOW_SPECIAL(mat,i-0,j-3,INTERGENIC);   
            }  
          return AlignWise_DC_SHADOW_MATRIX(mat,i - 0,j - 3,INTERGENIC) ;    
          }  
        /* Source CODON is not a special, should not get here! */ 
        warn("Major problem (!) - in AlignWise matrix to special read off, position %d,%d state %d no source found!",i,j,state); 
        return (-1); 
      case INTRON_0 :    
        /* Source INTRON_0 is not a special, should not get here! */ 
        /* Source CODON is not a special, should not get here! */ 
        warn("Major problem (!) - in AlignWise matrix to special read off, position %d,%d state %d no source found!",i,j,state); 
        return (-1); 
      case INTRON_1 :    
        /* Source INTRON_1 is not a special, should not get here! */ 
        /* Source CODON is not a special, should not get here! */ 
        warn("Major problem (!) - in AlignWise matrix to special read off, position %d,%d state %d no source found!",i,j,state); 
        return (-1); 
      case INTRON_2 :    
        /* Source INTRON_2 is not a special, should not get here! */ 
        /* Source CODON is not a special, should not get here! */ 
        warn("Major problem (!) - in AlignWise matrix to special read off, position %d,%d state %d no source found!",i,j,state); 
        return (-1); 
      case REVERSE_CODON :   
        /* Source REVERSE_INTRON_2 is not a special, should not get here! */ 
        /* Source REVERSE_INTRON_1 is not a special, should not get here! */ 
        /* Source REVERSE_INTRON_0 is not a special, should not get here! */ 
        temp = cscore - (mat->align->reverse_coding[j]) -  (0);  
        if( temp == AlignWise_DC_SHADOW_SPECIAL(mat,i - 0,j - 3,INTERGENIC) )    {  
          *reti = i - 0; 
          *retj = j - 3; 
          *retstate = INTERGENIC;    
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - AlignWise_DC_SHADOW_SPECIAL(mat,i-0,j-3,INTERGENIC);   
            }  
          return AlignWise_DC_SHADOW_MATRIX(mat,i - 0,j - 3,INTERGENIC) ;    
          }  
        /* Source REVERSE_CODON is not a special, should not get here! */ 
        warn("Major problem (!) - in AlignWise matrix to special read off, position %d,%d state %d no source found!",i,j,state); 
        return (-1); 
      case REVERSE_INTRON_0 :    
        /* Source REVERSE_INTRON_0 is not a special, should not get here! */ 
        /* Source REVERSE_CODON is not a special, should not get here! */ 
        warn("Major problem (!) - in AlignWise matrix to special read off, position %d,%d state %d no source found!",i,j,state); 
        return (-1); 
      case REVERSE_INTRON_1 :    
        /* Source REVERSE_INTRON_1 is not a special, should not get here! */ 
        /* Source REVERSE_CODON is not a special, should not get here! */ 
        warn("Major problem (!) - in AlignWise matrix to special read off, position %d,%d state %d no source found!",i,j,state); 
        return (-1); 
      case REVERSE_INTRON_2 :    
        /* Source REVERSE_INTRON_2 is not a special, should not get here! */ 
        /* Source REVERSE_CODON is not a special, should not get here! */ 
        warn("Major problem (!) - in AlignWise matrix to special read off, position %d,%d state %d no source found!",i,j,state); 
        return (-1); 
      default:   
        warn("Major problem (!) - in AlignWise read off, position %d,%d state %d no source found!",i,j,state);   
        return (-1); 
      } /* end of Switch state  */ 


}    


/* Function:  calculate_hidden_AlignWise(mat,starti,startj,startstate,stopi,stopj,stopstate,dpenv)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [AlignWise *]
 * Arg:            starti [UNKN ] Undocumented argument [int]
 * Arg:            startj [UNKN ] Undocumented argument [int]
 * Arg:        startstate [UNKN ] Undocumented argument [int]
 * Arg:             stopi [UNKN ] Undocumented argument [int]
 * Arg:             stopj [UNKN ] Undocumented argument [int]
 * Arg:         stopstate [UNKN ] Undocumented argument [int]
 * Arg:             dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 */
void calculate_hidden_AlignWise(AlignWise * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv) 
{
    register int i;  
    register int j;  
    register int score;  
    register int temp;   
    register int hiddenj;    


    hiddenj = startj;    


    init_hidden_AlignWise(mat,starti,startj,stopi,stopj);    


    AlignWise_HIDDEN_MATRIX(mat,starti,startj,startstate) = 0;   


    for(j=startj;j<=stopj;j++)   {  
      for(i=starti;i<=stopi;i++) {  
        /* Should *not* do very first cell as this is the one set to zero in one state! */ 
        if( i == starti && j == startj ) 
          continue;  
        if( dpenv != NULL && is_in_DPEnvelope(dpenv,i,j) == FALSE )  { /*Is not in envelope*/ 
          AlignWise_HIDDEN_MATRIX(mat,i,j,CODON) = NEGI;     
          AlignWise_HIDDEN_MATRIX(mat,i,j,INTRON_0) = NEGI;  
          AlignWise_HIDDEN_MATRIX(mat,i,j,INTRON_1) = NEGI;  
          AlignWise_HIDDEN_MATRIX(mat,i,j,INTRON_2) = NEGI;  
          AlignWise_HIDDEN_MATRIX(mat,i,j,REVERSE_CODON) = NEGI;     
          AlignWise_HIDDEN_MATRIX(mat,i,j,REVERSE_INTRON_0) = NEGI;  
          AlignWise_HIDDEN_MATRIX(mat,i,j,REVERSE_INTRON_1) = NEGI;  
          AlignWise_HIDDEN_MATRIX(mat,i,j,REVERSE_INTRON_2) = NEGI;  
          continue;  
          } /* end of Is not in envelope */ 


        /* For state CODON */ 
        /* setting first movement to score */ 
        score = AlignWise_HIDDEN_MATRIX(mat,i-0,j-3,CODON) + mat->align->forward_coding[j];  
        /* From state INTRON_0 to state CODON */ 
        temp = AlignWise_HIDDEN_MATRIX(mat,i-0,j-6,INTRON_0) + (mat->align->splice3_forward[(j-3)]+mat->align->forward_coding[j]);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTRON_1 to state CODON */ 
        temp = AlignWise_HIDDEN_MATRIX(mat,i-0,j-5,INTRON_1) + mat->align->splice3_forward[(j-2)];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTRON_2 to state CODON */ 
        temp = AlignWise_HIDDEN_MATRIX(mat,i-0,j-4,INTRON_2) + mat->align->splice3_forward[(j-1)];   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for CODON */ 
        /* Add any movement independant score and put away */ 
         AlignWise_HIDDEN_MATRIX(mat,i,j,CODON) = score; 
        /* Finished calculating state CODON */ 


        /* For state INTRON_0 */ 
        /* setting first movement to score */ 
        score = AlignWise_HIDDEN_MATRIX(mat,i-0,j-8,CODON) + (mat->align->splice5_forward[(j-7)]+mat->intronopen);   
        /* From state INTRON_0 to state INTRON_0 */ 
        temp = AlignWise_HIDDEN_MATRIX(mat,i-0,j-1,INTRON_0) + 0;    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for INTRON_0 */ 
        /* Add any movement independant score and put away */ 
         AlignWise_HIDDEN_MATRIX(mat,i,j,INTRON_0) = score;  
        /* Finished calculating state INTRON_0 */ 


        /* For state INTRON_1 */ 
        /* setting first movement to score */ 
        score = AlignWise_HIDDEN_MATRIX(mat,i-0,j-9,CODON) + (mat->align->splice5_forward[(j-7)]+mat->intronopen);   
        /* From state INTRON_1 to state INTRON_1 */ 
        temp = AlignWise_HIDDEN_MATRIX(mat,i-0,j-1,INTRON_1) + 0;    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for INTRON_1 */ 
        /* Add any movement independant score and put away */ 
         AlignWise_HIDDEN_MATRIX(mat,i,j,INTRON_1) = score;  
        /* Finished calculating state INTRON_1 */ 


        /* For state INTRON_2 */ 
        /* setting first movement to score */ 
        score = AlignWise_HIDDEN_MATRIX(mat,i-0,j-10,CODON) + (mat->align->splice5_forward[(j-7)]+mat->intronopen);  
        /* From state INTRON_2 to state INTRON_2 */ 
        temp = AlignWise_HIDDEN_MATRIX(mat,i-0,j-1,INTRON_2) + 0;    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for INTRON_2 */ 
        /* Add any movement independant score and put away */ 
         AlignWise_HIDDEN_MATRIX(mat,i,j,INTRON_2) = score;  
        /* Finished calculating state INTRON_2 */ 


        /* For state REVERSE_CODON */ 
        /* setting first movement to score */ 
        score = AlignWise_HIDDEN_MATRIX(mat,i-0,j-3,REVERSE_CODON) + mat->align->reverse_coding[j];  
        /* From state REVERSE_INTRON_0 to state REVERSE_CODON */ 
        temp = AlignWise_HIDDEN_MATRIX(mat,i-0,j-6,REVERSE_INTRON_0) + mat->align->splice5_reverse[(j-3)];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state REVERSE_INTRON_1 to state REVERSE_CODON */ 
        temp = AlignWise_HIDDEN_MATRIX(mat,i-0,j-5,REVERSE_INTRON_1) + mat->align->splice5_reverse[(j-2)];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state REVERSE_INTRON_2 to state REVERSE_CODON */ 
        temp = AlignWise_HIDDEN_MATRIX(mat,i-0,j-4,REVERSE_INTRON_2) + mat->align->splice5_reverse[(j-1)];   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for REVERSE_CODON */ 
        /* Add any movement independant score and put away */ 
         AlignWise_HIDDEN_MATRIX(mat,i,j,REVERSE_CODON) = score; 
        /* Finished calculating state REVERSE_CODON */ 


        /* For state REVERSE_INTRON_0 */ 
        /* setting first movement to score */ 
        score = AlignWise_HIDDEN_MATRIX(mat,i-0,j-8,REVERSE_CODON) + (mat->align->splice3_reverse[(j-7)]+mat->intronopen);   
        /* From state REVERSE_INTRON_0 to state REVERSE_INTRON_0 */ 
        temp = AlignWise_HIDDEN_MATRIX(mat,i-0,j-1,REVERSE_INTRON_0) + 0;    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for REVERSE_INTRON_0 */ 
        /* Add any movement independant score and put away */ 
         AlignWise_HIDDEN_MATRIX(mat,i,j,REVERSE_INTRON_0) = score;  
        /* Finished calculating state REVERSE_INTRON_0 */ 


        /* For state REVERSE_INTRON_1 */ 
        /* setting first movement to score */ 
        score = AlignWise_HIDDEN_MATRIX(mat,i-0,j-9,REVERSE_CODON) + (mat->align->splice3_reverse[(j-7)]+mat->intronopen);   
        /* From state REVERSE_INTRON_1 to state REVERSE_INTRON_1 */ 
        temp = AlignWise_HIDDEN_MATRIX(mat,i-0,j-1,REVERSE_INTRON_1) + 0;    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for REVERSE_INTRON_1 */ 
        /* Add any movement independant score and put away */ 
         AlignWise_HIDDEN_MATRIX(mat,i,j,REVERSE_INTRON_1) = score;  
        /* Finished calculating state REVERSE_INTRON_1 */ 


        /* For state REVERSE_INTRON_2 */ 
        /* setting first movement to score */ 
        score = AlignWise_HIDDEN_MATRIX(mat,i-0,j-10,REVERSE_CODON) + (mat->align->splice3_reverse[(j-7)]+mat->intronopen);  
        /* From state REVERSE_INTRON_2 to state REVERSE_INTRON_2 */ 
        temp = AlignWise_HIDDEN_MATRIX(mat,i-0,j-1,REVERSE_INTRON_2) + 0;    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for REVERSE_INTRON_2 */ 
        /* Add any movement independant score and put away */ 
         AlignWise_HIDDEN_MATRIX(mat,i,j,REVERSE_INTRON_2) = score;  
        /* Finished calculating state REVERSE_INTRON_2 */ 
        }  
      }  


    return;  
}    


/* Function:  init_hidden_AlignWise(mat,starti,startj,stopi,stopj)
 *
 * Descrip: No Description
 *
 * Arg:           mat [UNKN ] Undocumented argument [AlignWise *]
 * Arg:        starti [UNKN ] Undocumented argument [int]
 * Arg:        startj [UNKN ] Undocumented argument [int]
 * Arg:         stopi [UNKN ] Undocumented argument [int]
 * Arg:         stopj [UNKN ] Undocumented argument [int]
 *
 */
void init_hidden_AlignWise(AlignWise * mat,int starti,int startj,int stopi,int stopj) 
{
    register int i;  
    register int j;  
    register int hiddenj;    


    hiddenj = startj;    
    for(j=(startj-10);j<=stopj;j++)  {  
      for(i=(starti-0);i<=stopi;i++) {  
        AlignWise_HIDDEN_MATRIX(mat,i,j,CODON) = NEGI;
  
        AlignWise_HIDDEN_MATRIX(mat,i,j,INTRON_0) = NEGI;
   
        AlignWise_HIDDEN_MATRIX(mat,i,j,INTRON_1) = NEGI;
   
        AlignWise_HIDDEN_MATRIX(mat,i,j,INTRON_2) = NEGI;
   
        AlignWise_HIDDEN_MATRIX(mat,i,j,REVERSE_CODON) = NEGI;
  
        AlignWise_HIDDEN_MATRIX(mat,i,j,REVERSE_INTRON_0) = NEGI;
   
        AlignWise_HIDDEN_MATRIX(mat,i,j,REVERSE_INTRON_1) = NEGI;
   
        AlignWise_HIDDEN_MATRIX(mat,i,j,REVERSE_INTRON_2) = NEGI;
   
        }  
      }  


    return;  
}    


/* Function:  full_dc_AlignWise(mat,starti,startj,startstate,stopi,stopj,stopstate,out,donej,totalj,dpenv)
 *
 * Descrip:    The main divide-and-conquor routine. Basically, call /PackAln_calculate_small_AlignWise
 *             Not this function, which is pretty hard core. 
 *             Function is given start/end points (in main matrix) for alignment
 *             It does some checks, decides whether start/end in j is small enough for explicit calc
 *               - if yes, calculates it, reads off into PackAln (out), adds the j distance to donej and returns TRUE
 *               - if no,  uses /do_dc_single_pass_AlignWise to get mid-point
 *                          saves midpoint, and calls itself to do right portion then left portion
 *             right then left ensures PackAln is added the 'right' way, ie, back-to-front
 *             returns FALSE on any error, with a warning
 *
 *
 * Arg:               mat [UNKN ] Matrix with small memory implementation [AlignWise *]
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
boolean full_dc_AlignWise(AlignWise * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out,int * donej,int totalj,DPEnvelope * dpenv) 
{
    int lstarti; 
    int lstartj; 
    int lstate;  


    if( mat->basematrix->type != BASEMATRIX_TYPE_SHADOW) {  
      warn("*Very* bad error! - non shadow matrix type in full_dc_AlignWise");   
      return FALSE;  
      }  


    if( starti == -1 || startj == -1 || startstate == -1 || stopi == -1 || stopstate == -1)  {  
      warn("In full dc program, passed bad indices, indices passed were %d:%d[%d] to %d:%d[%d]\n",starti,startj,startstate,stopi,stopj,stopstate);   
      return FALSE;  
      }  


    if( stopj - startj < 50) {  
      log_full_error(REPORT,0,"[%d,%d][%d,%d] Explicit read off",starti,startj,stopi,stopj);/* Build hidden explicit matrix */ 
      calculate_hidden_AlignWise(mat,starti,startj,startstate,stopi,stopj,stopstate,dpenv);  
      *donej += (stopj - startj);   /* Now read it off into out */ 
      if( read_hidden_AlignWise(mat,starti,startj,startstate,stopi,stopj,stopstate,out) == FALSE)    {  
        warn("In full dc, at %d:%d,%d:%d got a bad hidden explicit read off... ",starti,startj,stopi,stopj); 
        return FALSE;    
        }  
      return TRUE;   
      }  


/* In actual divide and conquor */ 
    if( do_dc_single_pass_AlignWise(mat,starti,startj,startstate,stopi,stopj,stopstate,dpenv,(int)(*donej*100)/totalj) == FALSE) {  
      warn("In divide and conquor for AlignWise, at bound %d:%d to %d:%d, unable to calculate midpoint. Problem!",starti,startj,stopi,stopj);    
      return FALSE;  
      }  


/* Ok... now we have to call on each side of the matrix */ 
/* We have to retrieve left hand side positions, as they will be vapped by the time we call LHS */ 
    lstarti= AlignWise_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,0);     
    lstartj= AlignWise_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,1);     
    lstate = AlignWise_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,2);     


/* Call on right hand side: this lets us do the correct read off */ 
    if( full_dc_AlignWise(mat,AlignWise_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,3),AlignWise_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,4),AlignWise_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,5),stopi,stopj,stopstate,out,donej,totalj,dpenv) == FALSE)   {  
/* Warning already issued, simply chained back up to top */ 
      return FALSE;  
      }  
/* Call on left hand side */ 
    if( full_dc_AlignWise(mat,starti,startj,startstate,lstarti,lstartj,lstate,out,donej,totalj,dpenv) == FALSE)  {  
/* Warning already issued, simply chained back up to top */ 
      return FALSE;  
      }  


    return TRUE;     
}    


/* Function:  do_dc_single_pass_AlignWise(mat,starti,startj,startstate,stopi,stopj,stopstate,dpenv,perc_done)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [AlignWise *]
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
boolean do_dc_single_pass_AlignWise(AlignWise * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv,int perc_done) 
{
    int halfj;   
    halfj = startj + ((stopj - startj)/2);   


    init_dc_AlignWise(mat);  


    AlignWise_DC_SHADOW_MATRIX(mat,starti,startj,startstate) = 0;    
    run_up_dc_AlignWise(mat,starti,stopi,startj,halfj-1,dpenv,perc_done);    
    push_dc_at_merge_AlignWise(mat,starti,stopi,halfj,&halfj,dpenv);     
    follow_on_dc_AlignWise(mat,starti,stopi,halfj,stopj,dpenv,perc_done);    
    return TRUE; 
}    


/* Function:  push_dc_at_merge_AlignWise(mat,starti,stopi,startj,stopj,dpenv)
 *
 * Descrip: No Description
 *
 * Arg:           mat [UNKN ] Undocumented argument [AlignWise *]
 * Arg:        starti [UNKN ] Undocumented argument [int]
 * Arg:         stopi [UNKN ] Undocumented argument [int]
 * Arg:        startj [UNKN ] Undocumented argument [int]
 * Arg:         stopj [UNKN ] Undocumented argument [int *]
 * Arg:         dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 */
void push_dc_at_merge_AlignWise(AlignWise * mat,int starti,int stopi,int startj,int * stopj,DPEnvelope * dpenv) 
{
    register int i;  
    register int j;  
    register int k;  
    register int count;  
    register int mergej;/* Sources below this j will be stamped by triples */ 
    register int score;  
    register int temp;   


    mergej = startj -1;  
    for(count=0,j=startj;count<10;count++,j++)   {  
      for(i=starti;i<=stopi;i++) {  
        if( dpenv != NULL && is_in_DPEnvelope(dpenv,i,j) == FALSE )  { /*Is not in envelope*/ 
          AlignWise_DC_SHADOW_MATRIX(mat,i,j,CODON) = NEGI;  
          AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,CODON,0) = (-100);   
          AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,CODON,1) = (-100);   
          AlignWise_DC_SHADOW_MATRIX(mat,i,j,INTRON_0) = NEGI;   
          AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_0,0) = (-100);    
          AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_0,1) = (-100);    
          AlignWise_DC_SHADOW_MATRIX(mat,i,j,INTRON_1) = NEGI;   
          AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_1,0) = (-100);    
          AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_1,1) = (-100);    
          AlignWise_DC_SHADOW_MATRIX(mat,i,j,INTRON_2) = NEGI;   
          AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_2,0) = (-100);    
          AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_2,1) = (-100);    
          AlignWise_DC_SHADOW_MATRIX(mat,i,j,REVERSE_CODON) = NEGI;  
          AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_CODON,0) = (-100);   
          AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_CODON,1) = (-100);   
          AlignWise_DC_SHADOW_MATRIX(mat,i,j,REVERSE_INTRON_0) = NEGI;   
          AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_INTRON_0,0) = (-100);    
          AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_INTRON_0,1) = (-100);    
          AlignWise_DC_SHADOW_MATRIX(mat,i,j,REVERSE_INTRON_1) = NEGI;   
          AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_INTRON_1,0) = (-100);    
          AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_INTRON_1,1) = (-100);    
          AlignWise_DC_SHADOW_MATRIX(mat,i,j,REVERSE_INTRON_2) = NEGI;   
          AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_INTRON_2,0) = (-100);    
          AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_INTRON_2,1) = (-100);    
          continue;  
          } /* end of Is not in envelope */ 


        /* For state CODON, pushing when j - offj <= mergej */ 
        score = AlignWise_DC_SHADOW_MATRIX(mat,i-0,j-3,CODON) + mat->align->forward_coding[j];   
        if( j - 3 <= mergej) {  
          AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,CODON,0) = i-0;  
          AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,CODON,1) = j-3;  
          AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,CODON,2) = CODON;    
          AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,CODON,3) = i;    
          AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,CODON,4) = j;    
          AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,CODON,5) = CODON;    
          }  
        else {  
          for(k=0;k<7;k++)   
            AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,CODON,k) = AlignWise_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 3,CODON,k); 
          }  


        temp = AlignWise_DC_SHADOW_MATRIX(mat,i-0,j-6,INTRON_0) + (mat->align->splice3_forward[(j-3)]+mat->align->forward_coding[j]);    
        if( temp > score)    {  
          score = temp;  


          if( j - 6 <= mergej)   {  
            AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,CODON,0) = i-0;    
            AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,CODON,1) = j-6;    
            AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,CODON,2) = INTRON_0;   
            AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,CODON,3) = i;  
            AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,CODON,4) = j;  
            AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,CODON,5) = CODON;  
            }  
          else   {  
            for(k=0;k<7;k++) 
              AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,CODON,k) = AlignWise_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 6,INTRON_0,k);    
            }  
          }  


        temp = AlignWise_DC_SHADOW_MATRIX(mat,i-0,j-5,INTRON_1) + mat->align->splice3_forward[(j-2)];    
        if( temp > score)    {  
          score = temp;  


          if( j - 5 <= mergej)   {  
            AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,CODON,0) = i-0;    
            AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,CODON,1) = j-5;    
            AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,CODON,2) = INTRON_1;   
            AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,CODON,3) = i;  
            AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,CODON,4) = j;  
            AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,CODON,5) = CODON;  
            }  
          else   {  
            for(k=0;k<7;k++) 
              AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,CODON,k) = AlignWise_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 5,INTRON_1,k);    
            }  
          }  


        temp = AlignWise_DC_SHADOW_MATRIX(mat,i-0,j-4,INTRON_2) + mat->align->splice3_forward[(j-1)];    
        if( temp > score)    {  
          score = temp;  


          if( j - 4 <= mergej)   {  
            AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,CODON,0) = i-0;    
            AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,CODON,1) = j-4;    
            AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,CODON,2) = INTRON_2;   
            AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,CODON,3) = i;  
            AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,CODON,4) = j;  
            AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,CODON,5) = CODON;  
            }  
          else   {  
            for(k=0;k<7;k++) 
              AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,CODON,k) = AlignWise_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 4,INTRON_2,k);    
            }  
          }  
        /* Add any movement independant score */ 
        AlignWise_DC_SHADOW_MATRIX(mat,i,j,CODON) = score;   
        /* Finished with state CODON */ 


        /* For state INTRON_0, pushing when j - offj <= mergej */ 
        score = AlignWise_DC_SHADOW_MATRIX(mat,i-0,j-8,CODON) + (mat->align->splice5_forward[(j-7)]+mat->intronopen);    
        if( j - 8 <= mergej) {  
          AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_0,0) = i-0;   
          AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_0,1) = j-8;   
          AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_0,2) = CODON; 
          AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_0,3) = i; 
          AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_0,4) = j; 
          AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_0,5) = INTRON_0;  
          }  
        else {  
          for(k=0;k<7;k++)   
            AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_0,k) = AlignWise_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 8,CODON,k);  
          }  


        temp = AlignWise_DC_SHADOW_MATRIX(mat,i-0,j-1,INTRON_0) + 0;     
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_0,0) = i-0; 
            AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_0,1) = j-1; 
            AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_0,2) = INTRON_0;    
            AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_0,3) = i;   
            AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_0,4) = j;   
            AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_0,5) = INTRON_0;    
            }  
          else   {  
            for(k=0;k<7;k++) 
              AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_0,k) = AlignWise_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,INTRON_0,k); 
            }  
          }  
        /* Add any movement independant score */ 
        AlignWise_DC_SHADOW_MATRIX(mat,i,j,INTRON_0) = score;    
        /* Finished with state INTRON_0 */ 


        /* For state INTRON_1, pushing when j - offj <= mergej */ 
        score = AlignWise_DC_SHADOW_MATRIX(mat,i-0,j-9,CODON) + (mat->align->splice5_forward[(j-7)]+mat->intronopen);    
        if( j - 9 <= mergej) {  
          AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_1,0) = i-0;   
          AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_1,1) = j-9;   
          AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_1,2) = CODON; 
          AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_1,3) = i; 
          AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_1,4) = j; 
          AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_1,5) = INTRON_1;  
          }  
        else {  
          for(k=0;k<7;k++)   
            AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_1,k) = AlignWise_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 9,CODON,k);  
          }  


        temp = AlignWise_DC_SHADOW_MATRIX(mat,i-0,j-1,INTRON_1) + 0;     
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_1,0) = i-0; 
            AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_1,1) = j-1; 
            AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_1,2) = INTRON_1;    
            AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_1,3) = i;   
            AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_1,4) = j;   
            AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_1,5) = INTRON_1;    
            }  
          else   {  
            for(k=0;k<7;k++) 
              AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_1,k) = AlignWise_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,INTRON_1,k); 
            }  
          }  
        /* Add any movement independant score */ 
        AlignWise_DC_SHADOW_MATRIX(mat,i,j,INTRON_1) = score;    
        /* Finished with state INTRON_1 */ 


        /* For state INTRON_2, pushing when j - offj <= mergej */ 
        score = AlignWise_DC_SHADOW_MATRIX(mat,i-0,j-10,CODON) + (mat->align->splice5_forward[(j-7)]+mat->intronopen);   
        if( j - 10 <= mergej)    {  
          AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_2,0) = i-0;   
          AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_2,1) = j-10;  
          AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_2,2) = CODON; 
          AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_2,3) = i; 
          AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_2,4) = j; 
          AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_2,5) = INTRON_2;  
          }  
        else {  
          for(k=0;k<7;k++)   
            AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_2,k) = AlignWise_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 10,CODON,k); 
          }  


        temp = AlignWise_DC_SHADOW_MATRIX(mat,i-0,j-1,INTRON_2) + 0;     
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_2,0) = i-0; 
            AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_2,1) = j-1; 
            AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_2,2) = INTRON_2;    
            AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_2,3) = i;   
            AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_2,4) = j;   
            AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_2,5) = INTRON_2;    
            }  
          else   {  
            for(k=0;k<7;k++) 
              AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_2,k) = AlignWise_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,INTRON_2,k); 
            }  
          }  
        /* Add any movement independant score */ 
        AlignWise_DC_SHADOW_MATRIX(mat,i,j,INTRON_2) = score;    
        /* Finished with state INTRON_2 */ 


        /* For state REVERSE_CODON, pushing when j - offj <= mergej */ 
        score = AlignWise_DC_SHADOW_MATRIX(mat,i-0,j-3,REVERSE_CODON) + mat->align->reverse_coding[j];   
        if( j - 3 <= mergej) {  
          AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_CODON,0) = i-0;  
          AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_CODON,1) = j-3;  
          AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_CODON,2) = REVERSE_CODON;    
          AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_CODON,3) = i;    
          AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_CODON,4) = j;    
          AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_CODON,5) = REVERSE_CODON;    
          }  
        else {  
          for(k=0;k<7;k++)   
            AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_CODON,k) = AlignWise_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 3,REVERSE_CODON,k); 
          }  


        temp = AlignWise_DC_SHADOW_MATRIX(mat,i-0,j-6,REVERSE_INTRON_0) + mat->align->splice5_reverse[(j-3)];    
        if( temp > score)    {  
          score = temp;  


          if( j - 6 <= mergej)   {  
            AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_CODON,0) = i-0;    
            AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_CODON,1) = j-6;    
            AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_CODON,2) = REVERSE_INTRON_0;   
            AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_CODON,3) = i;  
            AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_CODON,4) = j;  
            AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_CODON,5) = REVERSE_CODON;  
            }  
          else   {  
            for(k=0;k<7;k++) 
              AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_CODON,k) = AlignWise_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 6,REVERSE_INTRON_0,k);    
            }  
          }  


        temp = AlignWise_DC_SHADOW_MATRIX(mat,i-0,j-5,REVERSE_INTRON_1) + mat->align->splice5_reverse[(j-2)];    
        if( temp > score)    {  
          score = temp;  


          if( j - 5 <= mergej)   {  
            AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_CODON,0) = i-0;    
            AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_CODON,1) = j-5;    
            AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_CODON,2) = REVERSE_INTRON_1;   
            AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_CODON,3) = i;  
            AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_CODON,4) = j;  
            AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_CODON,5) = REVERSE_CODON;  
            }  
          else   {  
            for(k=0;k<7;k++) 
              AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_CODON,k) = AlignWise_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 5,REVERSE_INTRON_1,k);    
            }  
          }  


        temp = AlignWise_DC_SHADOW_MATRIX(mat,i-0,j-4,REVERSE_INTRON_2) + mat->align->splice5_reverse[(j-1)];    
        if( temp > score)    {  
          score = temp;  


          if( j - 4 <= mergej)   {  
            AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_CODON,0) = i-0;    
            AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_CODON,1) = j-4;    
            AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_CODON,2) = REVERSE_INTRON_2;   
            AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_CODON,3) = i;  
            AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_CODON,4) = j;  
            AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_CODON,5) = REVERSE_CODON;  
            }  
          else   {  
            for(k=0;k<7;k++) 
              AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_CODON,k) = AlignWise_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 4,REVERSE_INTRON_2,k);    
            }  
          }  
        /* Add any movement independant score */ 
        AlignWise_DC_SHADOW_MATRIX(mat,i,j,REVERSE_CODON) = score;   
        /* Finished with state REVERSE_CODON */ 


        /* For state REVERSE_INTRON_0, pushing when j - offj <= mergej */ 
        score = AlignWise_DC_SHADOW_MATRIX(mat,i-0,j-8,REVERSE_CODON) + (mat->align->splice3_reverse[(j-7)]+mat->intronopen);    
        if( j - 8 <= mergej) {  
          AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_INTRON_0,0) = i-0;   
          AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_INTRON_0,1) = j-8;   
          AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_INTRON_0,2) = REVERSE_CODON; 
          AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_INTRON_0,3) = i; 
          AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_INTRON_0,4) = j; 
          AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_INTRON_0,5) = REVERSE_INTRON_0;  
          }  
        else {  
          for(k=0;k<7;k++)   
            AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_INTRON_0,k) = AlignWise_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 8,REVERSE_CODON,k);  
          }  


        temp = AlignWise_DC_SHADOW_MATRIX(mat,i-0,j-1,REVERSE_INTRON_0) + 0;     
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_INTRON_0,0) = i-0; 
            AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_INTRON_0,1) = j-1; 
            AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_INTRON_0,2) = REVERSE_INTRON_0;    
            AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_INTRON_0,3) = i;   
            AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_INTRON_0,4) = j;   
            AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_INTRON_0,5) = REVERSE_INTRON_0;    
            }  
          else   {  
            for(k=0;k<7;k++) 
              AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_INTRON_0,k) = AlignWise_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,REVERSE_INTRON_0,k); 
            }  
          }  
        /* Add any movement independant score */ 
        AlignWise_DC_SHADOW_MATRIX(mat,i,j,REVERSE_INTRON_0) = score;    
        /* Finished with state REVERSE_INTRON_0 */ 


        /* For state REVERSE_INTRON_1, pushing when j - offj <= mergej */ 
        score = AlignWise_DC_SHADOW_MATRIX(mat,i-0,j-9,REVERSE_CODON) + (mat->align->splice3_reverse[(j-7)]+mat->intronopen);    
        if( j - 9 <= mergej) {  
          AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_INTRON_1,0) = i-0;   
          AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_INTRON_1,1) = j-9;   
          AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_INTRON_1,2) = REVERSE_CODON; 
          AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_INTRON_1,3) = i; 
          AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_INTRON_1,4) = j; 
          AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_INTRON_1,5) = REVERSE_INTRON_1;  
          }  
        else {  
          for(k=0;k<7;k++)   
            AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_INTRON_1,k) = AlignWise_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 9,REVERSE_CODON,k);  
          }  


        temp = AlignWise_DC_SHADOW_MATRIX(mat,i-0,j-1,REVERSE_INTRON_1) + 0;     
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_INTRON_1,0) = i-0; 
            AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_INTRON_1,1) = j-1; 
            AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_INTRON_1,2) = REVERSE_INTRON_1;    
            AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_INTRON_1,3) = i;   
            AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_INTRON_1,4) = j;   
            AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_INTRON_1,5) = REVERSE_INTRON_1;    
            }  
          else   {  
            for(k=0;k<7;k++) 
              AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_INTRON_1,k) = AlignWise_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,REVERSE_INTRON_1,k); 
            }  
          }  
        /* Add any movement independant score */ 
        AlignWise_DC_SHADOW_MATRIX(mat,i,j,REVERSE_INTRON_1) = score;    
        /* Finished with state REVERSE_INTRON_1 */ 


        /* For state REVERSE_INTRON_2, pushing when j - offj <= mergej */ 
        score = AlignWise_DC_SHADOW_MATRIX(mat,i-0,j-10,REVERSE_CODON) + (mat->align->splice3_reverse[(j-7)]+mat->intronopen);   
        if( j - 10 <= mergej)    {  
          AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_INTRON_2,0) = i-0;   
          AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_INTRON_2,1) = j-10;  
          AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_INTRON_2,2) = REVERSE_CODON; 
          AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_INTRON_2,3) = i; 
          AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_INTRON_2,4) = j; 
          AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_INTRON_2,5) = REVERSE_INTRON_2;  
          }  
        else {  
          for(k=0;k<7;k++)   
            AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_INTRON_2,k) = AlignWise_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 10,REVERSE_CODON,k); 
          }  


        temp = AlignWise_DC_SHADOW_MATRIX(mat,i-0,j-1,REVERSE_INTRON_2) + 0;     
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_INTRON_2,0) = i-0; 
            AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_INTRON_2,1) = j-1; 
            AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_INTRON_2,2) = REVERSE_INTRON_2;    
            AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_INTRON_2,3) = i;   
            AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_INTRON_2,4) = j;   
            AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_INTRON_2,5) = REVERSE_INTRON_2;    
            }  
          else   {  
            for(k=0;k<7;k++) 
              AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_INTRON_2,k) = AlignWise_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,REVERSE_INTRON_2,k); 
            }  
          }  
        /* Add any movement independant score */ 
        AlignWise_DC_SHADOW_MATRIX(mat,i,j,REVERSE_INTRON_2) = score;    
        /* Finished with state REVERSE_INTRON_2 */ 
        }  
      }  
    /* Put back j into * stop j so that calling function gets it correct */ 
    if( stopj == NULL)   
      warn("Bad news... NULL stopj pointer in push dc function. This means that calling function does not know how many cells I have done!");    
    else 
      *stopj = j;    


    return;  
}    


/* Function:  follow_on_dc_AlignWise(mat,starti,stopi,startj,stopj,dpenv,perc_done)
 *
 * Descrip: No Description
 *
 * Arg:              mat [UNKN ] Undocumented argument [AlignWise *]
 * Arg:           starti [UNKN ] Undocumented argument [int]
 * Arg:            stopi [UNKN ] Undocumented argument [int]
 * Arg:           startj [UNKN ] Undocumented argument [int]
 * Arg:            stopj [UNKN ] Undocumented argument [int]
 * Arg:            dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:        perc_done [UNKN ] Undocumented argument [int]
 *
 */
void follow_on_dc_AlignWise(AlignWise * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done) 
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
          AlignWise_DC_SHADOW_MATRIX(mat,i,j,CODON) = NEGI;  
          AlignWise_DC_SHADOW_MATRIX(mat,i,j,INTRON_0) = NEGI;   
          AlignWise_DC_SHADOW_MATRIX(mat,i,j,INTRON_1) = NEGI;   
          AlignWise_DC_SHADOW_MATRIX(mat,i,j,INTRON_2) = NEGI;   
          AlignWise_DC_SHADOW_MATRIX(mat,i,j,REVERSE_CODON) = NEGI;  
          AlignWise_DC_SHADOW_MATRIX(mat,i,j,REVERSE_INTRON_0) = NEGI;   
          AlignWise_DC_SHADOW_MATRIX(mat,i,j,REVERSE_INTRON_1) = NEGI;   
          AlignWise_DC_SHADOW_MATRIX(mat,i,j,REVERSE_INTRON_2) = NEGI;   
          continue;  
          } /* end of Is not in envelope */ 
        if( num % 1000 == 0 )    
          log_full_error(REPORT,0,"[%d%%%% done]After  mid-j %5d Cells done %d%%%%",perc_done,startj,(num*100)/total);   


        /* For state CODON */ 
        /* setting first movement to score */ 
        score = AlignWise_DC_SHADOW_MATRIX(mat,i-0,j-3,CODON) + mat->align->forward_coding[j];   
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = AlignWise_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 3,CODON,k);   
        /* From state INTRON_0 to state CODON */ 
        temp = AlignWise_DC_SHADOW_MATRIX(mat,i-0,j-6,INTRON_0) + (mat->align->splice3_forward[(j-3)]+mat->align->forward_coding[j]);    
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = AlignWise_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 6,INTRON_0,k);  
          }  
        /* From state INTRON_1 to state CODON */ 
        temp = AlignWise_DC_SHADOW_MATRIX(mat,i-0,j-5,INTRON_1) + mat->align->splice3_forward[(j-2)];    
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = AlignWise_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 5,INTRON_1,k);  
          }  
        /* From state INTRON_2 to state CODON */ 
        temp = AlignWise_DC_SHADOW_MATRIX(mat,i-0,j-4,INTRON_2) + mat->align->splice3_forward[(j-1)];    
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = AlignWise_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 4,INTRON_2,k);  
          }  


        /* Ok - finished max calculation for CODON */ 
        /* Add any movement independant score and put away */ 
         AlignWise_DC_SHADOW_MATRIX(mat,i,j,CODON) = score;  
        for(k=0;k<7;k++) 
          AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,CODON,k) = localshadow[k];   
        /* Now figure out if any specials need this score */ 
        /* Finished calculating state CODON */ 


        /* For state INTRON_0 */ 
        /* setting first movement to score */ 
        score = AlignWise_DC_SHADOW_MATRIX(mat,i-0,j-8,CODON) + (mat->align->splice5_forward[(j-7)]+mat->intronopen);    
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = AlignWise_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 8,CODON,k);   
        /* From state INTRON_0 to state INTRON_0 */ 
        temp = AlignWise_DC_SHADOW_MATRIX(mat,i-0,j-1,INTRON_0) + 0;     
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = AlignWise_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,INTRON_0,k);  
          }  


        /* Ok - finished max calculation for INTRON_0 */ 
        /* Add any movement independant score and put away */ 
         AlignWise_DC_SHADOW_MATRIX(mat,i,j,INTRON_0) = score;   
        for(k=0;k<7;k++) 
          AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_0,k) = localshadow[k];    
        /* Now figure out if any specials need this score */ 
        /* Finished calculating state INTRON_0 */ 


        /* For state INTRON_1 */ 
        /* setting first movement to score */ 
        score = AlignWise_DC_SHADOW_MATRIX(mat,i-0,j-9,CODON) + (mat->align->splice5_forward[(j-7)]+mat->intronopen);    
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = AlignWise_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 9,CODON,k);   
        /* From state INTRON_1 to state INTRON_1 */ 
        temp = AlignWise_DC_SHADOW_MATRIX(mat,i-0,j-1,INTRON_1) + 0;     
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = AlignWise_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,INTRON_1,k);  
          }  


        /* Ok - finished max calculation for INTRON_1 */ 
        /* Add any movement independant score and put away */ 
         AlignWise_DC_SHADOW_MATRIX(mat,i,j,INTRON_1) = score;   
        for(k=0;k<7;k++) 
          AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_1,k) = localshadow[k];    
        /* Now figure out if any specials need this score */ 
        /* Finished calculating state INTRON_1 */ 


        /* For state INTRON_2 */ 
        /* setting first movement to score */ 
        score = AlignWise_DC_SHADOW_MATRIX(mat,i-0,j-10,CODON) + (mat->align->splice5_forward[(j-7)]+mat->intronopen);   
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = AlignWise_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 10,CODON,k);  
        /* From state INTRON_2 to state INTRON_2 */ 
        temp = AlignWise_DC_SHADOW_MATRIX(mat,i-0,j-1,INTRON_2) + 0;     
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = AlignWise_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,INTRON_2,k);  
          }  


        /* Ok - finished max calculation for INTRON_2 */ 
        /* Add any movement independant score and put away */ 
         AlignWise_DC_SHADOW_MATRIX(mat,i,j,INTRON_2) = score;   
        for(k=0;k<7;k++) 
          AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_2,k) = localshadow[k];    
        /* Now figure out if any specials need this score */ 
        /* Finished calculating state INTRON_2 */ 


        /* For state REVERSE_CODON */ 
        /* setting first movement to score */ 
        score = AlignWise_DC_SHADOW_MATRIX(mat,i-0,j-3,REVERSE_CODON) + mat->align->reverse_coding[j];   
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = AlignWise_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 3,REVERSE_CODON,k);   
        /* From state REVERSE_INTRON_0 to state REVERSE_CODON */ 
        temp = AlignWise_DC_SHADOW_MATRIX(mat,i-0,j-6,REVERSE_INTRON_0) + mat->align->splice5_reverse[(j-3)];    
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = AlignWise_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 6,REVERSE_INTRON_0,k);  
          }  
        /* From state REVERSE_INTRON_1 to state REVERSE_CODON */ 
        temp = AlignWise_DC_SHADOW_MATRIX(mat,i-0,j-5,REVERSE_INTRON_1) + mat->align->splice5_reverse[(j-2)];    
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = AlignWise_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 5,REVERSE_INTRON_1,k);  
          }  
        /* From state REVERSE_INTRON_2 to state REVERSE_CODON */ 
        temp = AlignWise_DC_SHADOW_MATRIX(mat,i-0,j-4,REVERSE_INTRON_2) + mat->align->splice5_reverse[(j-1)];    
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = AlignWise_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 4,REVERSE_INTRON_2,k);  
          }  


        /* Ok - finished max calculation for REVERSE_CODON */ 
        /* Add any movement independant score and put away */ 
         AlignWise_DC_SHADOW_MATRIX(mat,i,j,REVERSE_CODON) = score;  
        for(k=0;k<7;k++) 
          AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_CODON,k) = localshadow[k];   
        /* Now figure out if any specials need this score */ 
        /* Finished calculating state REVERSE_CODON */ 


        /* For state REVERSE_INTRON_0 */ 
        /* setting first movement to score */ 
        score = AlignWise_DC_SHADOW_MATRIX(mat,i-0,j-8,REVERSE_CODON) + (mat->align->splice3_reverse[(j-7)]+mat->intronopen);    
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = AlignWise_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 8,REVERSE_CODON,k);   
        /* From state REVERSE_INTRON_0 to state REVERSE_INTRON_0 */ 
        temp = AlignWise_DC_SHADOW_MATRIX(mat,i-0,j-1,REVERSE_INTRON_0) + 0;     
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = AlignWise_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,REVERSE_INTRON_0,k);  
          }  


        /* Ok - finished max calculation for REVERSE_INTRON_0 */ 
        /* Add any movement independant score and put away */ 
         AlignWise_DC_SHADOW_MATRIX(mat,i,j,REVERSE_INTRON_0) = score;   
        for(k=0;k<7;k++) 
          AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_INTRON_0,k) = localshadow[k];    
        /* Now figure out if any specials need this score */ 
        /* Finished calculating state REVERSE_INTRON_0 */ 


        /* For state REVERSE_INTRON_1 */ 
        /* setting first movement to score */ 
        score = AlignWise_DC_SHADOW_MATRIX(mat,i-0,j-9,REVERSE_CODON) + (mat->align->splice3_reverse[(j-7)]+mat->intronopen);    
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = AlignWise_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 9,REVERSE_CODON,k);   
        /* From state REVERSE_INTRON_1 to state REVERSE_INTRON_1 */ 
        temp = AlignWise_DC_SHADOW_MATRIX(mat,i-0,j-1,REVERSE_INTRON_1) + 0;     
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = AlignWise_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,REVERSE_INTRON_1,k);  
          }  


        /* Ok - finished max calculation for REVERSE_INTRON_1 */ 
        /* Add any movement independant score and put away */ 
         AlignWise_DC_SHADOW_MATRIX(mat,i,j,REVERSE_INTRON_1) = score;   
        for(k=0;k<7;k++) 
          AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_INTRON_1,k) = localshadow[k];    
        /* Now figure out if any specials need this score */ 
        /* Finished calculating state REVERSE_INTRON_1 */ 


        /* For state REVERSE_INTRON_2 */ 
        /* setting first movement to score */ 
        score = AlignWise_DC_SHADOW_MATRIX(mat,i-0,j-10,REVERSE_CODON) + (mat->align->splice3_reverse[(j-7)]+mat->intronopen);   
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = AlignWise_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 10,REVERSE_CODON,k);  
        /* From state REVERSE_INTRON_2 to state REVERSE_INTRON_2 */ 
        temp = AlignWise_DC_SHADOW_MATRIX(mat,i-0,j-1,REVERSE_INTRON_2) + 0;     
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = AlignWise_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,REVERSE_INTRON_2,k);  
          }  


        /* Ok - finished max calculation for REVERSE_INTRON_2 */ 
        /* Add any movement independant score and put away */ 
         AlignWise_DC_SHADOW_MATRIX(mat,i,j,REVERSE_INTRON_2) = score;   
        for(k=0;k<7;k++) 
          AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_INTRON_2,k) = localshadow[k];    
        /* Now figure out if any specials need this score */ 
        /* Finished calculating state REVERSE_INTRON_2 */ 
        } /* end of this is strip */ 
      } /* end of for each valid j column */ 


/* Function:  run_up_dc_AlignWise(mat,starti,stopi,startj,stopj,dpenv,perc_done)
 *
 * Descrip: No Description
 *
 * Arg:              mat [UNKN ] Undocumented argument [AlignWise *]
 * Arg:           starti [UNKN ] Undocumented argument [int]
 * Arg:            stopi [UNKN ] Undocumented argument [int]
 * Arg:           startj [UNKN ] Undocumented argument [int]
 * Arg:            stopj [UNKN ] Undocumented argument [int]
 * Arg:            dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:        perc_done [UNKN ] Undocumented argument [int]
 *
 */
}    
void run_up_dc_AlignWise(AlignWise * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done) 
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
          AlignWise_DC_SHADOW_MATRIX(mat,i,j,CODON) = NEGI;  
          AlignWise_DC_SHADOW_MATRIX(mat,i,j,INTRON_0) = NEGI;   
          AlignWise_DC_SHADOW_MATRIX(mat,i,j,INTRON_1) = NEGI;   
          AlignWise_DC_SHADOW_MATRIX(mat,i,j,INTRON_2) = NEGI;   
          AlignWise_DC_SHADOW_MATRIX(mat,i,j,REVERSE_CODON) = NEGI;  
          AlignWise_DC_SHADOW_MATRIX(mat,i,j,REVERSE_INTRON_0) = NEGI;   
          AlignWise_DC_SHADOW_MATRIX(mat,i,j,REVERSE_INTRON_1) = NEGI;   
          AlignWise_DC_SHADOW_MATRIX(mat,i,j,REVERSE_INTRON_2) = NEGI;   
          continue;  
          } /* end of Is not in envelope */ 
        if( num % 1000 == 0 )    
          log_full_error(REPORT,0,"[%d%%%% done]Before mid-j %5d Cells done %d%%%%",perc_done,stopj,(num*100)/total);    


        /* For state CODON */ 
        /* setting first movement to score */ 
        score = AlignWise_DC_SHADOW_MATRIX(mat,i-0,j-3,CODON) + mat->align->forward_coding[j];   
        /* From state INTRON_0 to state CODON */ 
        temp = AlignWise_DC_SHADOW_MATRIX(mat,i-0,j-6,INTRON_0) + (mat->align->splice3_forward[(j-3)]+mat->align->forward_coding[j]);    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTRON_1 to state CODON */ 
        temp = AlignWise_DC_SHADOW_MATRIX(mat,i-0,j-5,INTRON_1) + mat->align->splice3_forward[(j-2)];    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTRON_2 to state CODON */ 
        temp = AlignWise_DC_SHADOW_MATRIX(mat,i-0,j-4,INTRON_2) + mat->align->splice3_forward[(j-1)];    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for CODON */ 
        /* Add any movement independant score and put away */ 
         AlignWise_DC_SHADOW_MATRIX(mat,i,j,CODON) = score;  
        /* Finished calculating state CODON */ 


        /* For state INTRON_0 */ 
        /* setting first movement to score */ 
        score = AlignWise_DC_SHADOW_MATRIX(mat,i-0,j-8,CODON) + (mat->align->splice5_forward[(j-7)]+mat->intronopen);    
        /* From state INTRON_0 to state INTRON_0 */ 
        temp = AlignWise_DC_SHADOW_MATRIX(mat,i-0,j-1,INTRON_0) + 0;     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for INTRON_0 */ 
        /* Add any movement independant score and put away */ 
         AlignWise_DC_SHADOW_MATRIX(mat,i,j,INTRON_0) = score;   
        /* Finished calculating state INTRON_0 */ 


        /* For state INTRON_1 */ 
        /* setting first movement to score */ 
        score = AlignWise_DC_SHADOW_MATRIX(mat,i-0,j-9,CODON) + (mat->align->splice5_forward[(j-7)]+mat->intronopen);    
        /* From state INTRON_1 to state INTRON_1 */ 
        temp = AlignWise_DC_SHADOW_MATRIX(mat,i-0,j-1,INTRON_1) + 0;     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for INTRON_1 */ 
        /* Add any movement independant score and put away */ 
         AlignWise_DC_SHADOW_MATRIX(mat,i,j,INTRON_1) = score;   
        /* Finished calculating state INTRON_1 */ 


        /* For state INTRON_2 */ 
        /* setting first movement to score */ 
        score = AlignWise_DC_SHADOW_MATRIX(mat,i-0,j-10,CODON) + (mat->align->splice5_forward[(j-7)]+mat->intronopen);   
        /* From state INTRON_2 to state INTRON_2 */ 
        temp = AlignWise_DC_SHADOW_MATRIX(mat,i-0,j-1,INTRON_2) + 0;     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for INTRON_2 */ 
        /* Add any movement independant score and put away */ 
         AlignWise_DC_SHADOW_MATRIX(mat,i,j,INTRON_2) = score;   
        /* Finished calculating state INTRON_2 */ 


        /* For state REVERSE_CODON */ 
        /* setting first movement to score */ 
        score = AlignWise_DC_SHADOW_MATRIX(mat,i-0,j-3,REVERSE_CODON) + mat->align->reverse_coding[j];   
        /* From state REVERSE_INTRON_0 to state REVERSE_CODON */ 
        temp = AlignWise_DC_SHADOW_MATRIX(mat,i-0,j-6,REVERSE_INTRON_0) + mat->align->splice5_reverse[(j-3)];    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state REVERSE_INTRON_1 to state REVERSE_CODON */ 
        temp = AlignWise_DC_SHADOW_MATRIX(mat,i-0,j-5,REVERSE_INTRON_1) + mat->align->splice5_reverse[(j-2)];    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state REVERSE_INTRON_2 to state REVERSE_CODON */ 
        temp = AlignWise_DC_SHADOW_MATRIX(mat,i-0,j-4,REVERSE_INTRON_2) + mat->align->splice5_reverse[(j-1)];    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for REVERSE_CODON */ 
        /* Add any movement independant score and put away */ 
         AlignWise_DC_SHADOW_MATRIX(mat,i,j,REVERSE_CODON) = score;  
        /* Finished calculating state REVERSE_CODON */ 


        /* For state REVERSE_INTRON_0 */ 
        /* setting first movement to score */ 
        score = AlignWise_DC_SHADOW_MATRIX(mat,i-0,j-8,REVERSE_CODON) + (mat->align->splice3_reverse[(j-7)]+mat->intronopen);    
        /* From state REVERSE_INTRON_0 to state REVERSE_INTRON_0 */ 
        temp = AlignWise_DC_SHADOW_MATRIX(mat,i-0,j-1,REVERSE_INTRON_0) + 0;     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for REVERSE_INTRON_0 */ 
        /* Add any movement independant score and put away */ 
         AlignWise_DC_SHADOW_MATRIX(mat,i,j,REVERSE_INTRON_0) = score;   
        /* Finished calculating state REVERSE_INTRON_0 */ 


        /* For state REVERSE_INTRON_1 */ 
        /* setting first movement to score */ 
        score = AlignWise_DC_SHADOW_MATRIX(mat,i-0,j-9,REVERSE_CODON) + (mat->align->splice3_reverse[(j-7)]+mat->intronopen);    
        /* From state REVERSE_INTRON_1 to state REVERSE_INTRON_1 */ 
        temp = AlignWise_DC_SHADOW_MATRIX(mat,i-0,j-1,REVERSE_INTRON_1) + 0;     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for REVERSE_INTRON_1 */ 
        /* Add any movement independant score and put away */ 
         AlignWise_DC_SHADOW_MATRIX(mat,i,j,REVERSE_INTRON_1) = score;   
        /* Finished calculating state REVERSE_INTRON_1 */ 


        /* For state REVERSE_INTRON_2 */ 
        /* setting first movement to score */ 
        score = AlignWise_DC_SHADOW_MATRIX(mat,i-0,j-10,REVERSE_CODON) + (mat->align->splice3_reverse[(j-7)]+mat->intronopen);   
        /* From state REVERSE_INTRON_2 to state REVERSE_INTRON_2 */ 
        temp = AlignWise_DC_SHADOW_MATRIX(mat,i-0,j-1,REVERSE_INTRON_2) + 0;     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for REVERSE_INTRON_2 */ 
        /* Add any movement independant score and put away */ 
         AlignWise_DC_SHADOW_MATRIX(mat,i,j,REVERSE_INTRON_2) = score;   
        /* Finished calculating state REVERSE_INTRON_2 */ 
        } /* end of this is strip */ 
      } /* end of for each valid j column */ 


/* Function:  init_dc_AlignWise(mat)
 *
 * Descrip: No Description
 *
 * Arg:        mat [UNKN ] Undocumented argument [AlignWise *]
 *
 */
}    
void init_dc_AlignWise(AlignWise * mat) 
{
    register int i;  
    register int j;  
    register int k;  


    for(j=0;j<12;j++)    {  
      for(i=(-0);i<mat->model->len;i++)  {  
        AlignWise_DC_SHADOW_MATRIX(mat,i,j,CODON) = NEGI;    
        AlignWise_DC_SHADOW_MATRIX(mat,i,j,INTRON_0) = NEGI; 
        AlignWise_DC_SHADOW_MATRIX(mat,i,j,INTRON_1) = NEGI; 
        AlignWise_DC_SHADOW_MATRIX(mat,i,j,INTRON_2) = NEGI; 
        AlignWise_DC_SHADOW_MATRIX(mat,i,j,REVERSE_CODON) = NEGI;    
        AlignWise_DC_SHADOW_MATRIX(mat,i,j,REVERSE_INTRON_0) = NEGI; 
        AlignWise_DC_SHADOW_MATRIX(mat,i,j,REVERSE_INTRON_1) = NEGI; 
        AlignWise_DC_SHADOW_MATRIX(mat,i,j,REVERSE_INTRON_2) = NEGI; 
        for(k=0;k<7;k++) {  
          AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,CODON,k) = (-1); 
          AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_0,k) = (-1);  
          AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_1,k) = (-1);  
          AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_2,k) = (-1);  
          AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_CODON,k) = (-1); 
          AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_INTRON_0,k) = (-1);  
          AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_INTRON_1,k) = (-1);  
          AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_INTRON_2,k) = (-1);  
          }  
        }  
      }  


    return;  
}    


/* Function:  start_end_find_end_AlignWise(mat,endj)
 *
 * Descrip:    First function used to find end of the best path in the special state !end
 *
 *
 * Arg:         mat [UNKN ] Matrix in small mode [AlignWise *]
 * Arg:        endj [WRITE] position of end in j (meaningless in i) [int *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int start_end_find_end_AlignWise(AlignWise * mat,int * endj) 
{
    register int j;  
    register int max;    
    register int maxj;   


    max = AlignWise_DC_SHADOW_SPECIAL(mat,0,mat->align->len-1,END);  
    maxj = mat->align->len-1;    
    for(j= mat->align->len-2 ;j >= 0 ;j--)   {  
      if( AlignWise_DC_SHADOW_SPECIAL(mat,0,j,END) > max )   {  
        max = AlignWise_DC_SHADOW_SPECIAL(mat,0,j,END);  
        maxj = j;    
        }  
      }  


    if( endj != NULL)    
      *endj = maxj;  


    return max;  
}    


/* Function:  dc_optimised_start_end_calc_AlignWise(*mat,dpenv)
 *
 * Descrip:    Calculates special strip, leaving start/end/score points in shadow matrix
 *             Works off specially laid out memory from steve searle
 *
 *
 * Arg:         *mat [UNKN ] Undocumented argument [AlignWise]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean dc_optimised_start_end_calc_AlignWise(AlignWise *mat,DPEnvelope * dpenv) 
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
    leni = mat->model->len;  
    lenj = mat->align->len;  
    total = leni * lenj; 


    score_pointers = (int *) calloc (10 * (leni + 0) * 8,sizeof(int));   
    shadow_pointers = (int *) calloc (10 * (leni + 0) * 8 * 8,sizeof(int));  


    for(j=0;j<lenj;j++)  { /*for each j strip*/ 
      for(i=0;i<leni;i++)    { /*for each i position in strip*/ 
        num++;   
        if( dpenv != NULL && is_in_DPEnvelope(dpenv,i,j) == FALSE )  { /*Is not in envelope*/ 
          AlignWise_DC_OPT_SHADOW_MATRIX(mat,i,j,CODON) = NEGI;  
          AlignWise_DC_OPT_SHADOW_MATRIX(mat,i,j,INTRON_0) = NEGI;   
          AlignWise_DC_OPT_SHADOW_MATRIX(mat,i,j,INTRON_1) = NEGI;   
          AlignWise_DC_OPT_SHADOW_MATRIX(mat,i,j,INTRON_2) = NEGI;   
          AlignWise_DC_OPT_SHADOW_MATRIX(mat,i,j,REVERSE_CODON) = NEGI;  
          AlignWise_DC_OPT_SHADOW_MATRIX(mat,i,j,REVERSE_INTRON_0) = NEGI;   
          AlignWise_DC_OPT_SHADOW_MATRIX(mat,i,j,REVERSE_INTRON_1) = NEGI;   
          AlignWise_DC_OPT_SHADOW_MATRIX(mat,i,j,REVERSE_INTRON_2) = NEGI;   
          continue;  
          } /* end of Is not in envelope */ 
        if( num%1000 == 0)   
          log_full_error(REPORT,0,"%6d Cells done [%2d%%%%]",num,num*100/total); 




        /* For state CODON */ 
        /* setting first movement to score */ 
        score = AlignWise_DC_OPT_SHADOW_MATRIX(mat,i-0,j-3,CODON) + mat->align->forward_coding[j] + (0);     
        /* assign local shadown pointer */ 
        localsp = &(AlignWise_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 3,CODON,0)); 
        /* From state INTERGENIC to state CODON */ 
        temp = AlignWise_DC_OPT_SHADOW_SPECIAL(mat,i-0,j-3,INTERGENIC) + mat->align->forward_coding[j] + (0);    
        if( temp  > score )  {  
          score = temp;  
          /* This state [INTERGENIC] is a special for CODON... push top shadow pointers here */ 
          localshadow[0]= i; 
          localshadow[1]= j; 
          localshadow[2]= CODON; 
          localshadow[3]= (-1);  
          localshadow[4]= (-1);  
          localshadow[5]= (-1);  
          localshadow[6]= score; 
          localsp = localshadow; 
          }  
        /* From state INTRON_0 to state CODON */ 
        temp = AlignWise_DC_OPT_SHADOW_MATRIX(mat,i-0,j-6,INTRON_0) + (mat->align->splice3_forward[(j-3)]+mat->align->forward_coding[j]) +(0);   
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(AlignWise_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 6,INTRON_0,0));    
          }  
        /* From state INTRON_1 to state CODON */ 
        temp = AlignWise_DC_OPT_SHADOW_MATRIX(mat,i-0,j-5,INTRON_1) + mat->align->splice3_forward[(j-2)] +(0);   
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(AlignWise_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 5,INTRON_1,0));    
          }  
        /* From state INTRON_2 to state CODON */ 
        temp = AlignWise_DC_OPT_SHADOW_MATRIX(mat,i-0,j-4,INTRON_2) + mat->align->splice3_forward[(j-1)] +(0);   
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(AlignWise_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 4,INTRON_2,0));    
          }  


        /* Ok - finished max calculation for CODON */ 
        /* Add any movement independant score and put away */ 
        /* Actually, already done inside scores */ 
         AlignWise_DC_OPT_SHADOW_MATRIX(mat,i,j,CODON) = score;  
        for(k=0;k<7;k++) 
          AlignWise_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,CODON,k) = localsp[k];   
        /* Now figure out if any specials need this score */ 


        /* state CODON is a source for special INTERGENIC */ 
        temp = score + (mat->geneopen) + (0) ;   
        if( temp > AlignWise_DC_OPT_SHADOW_SPECIAL(mat,i,j,INTERGENIC) )     {  
          AlignWise_DC_OPT_SHADOW_SPECIAL(mat,i,j,INTERGENIC) = temp;    
          /* Have to push only bottem half of system here */ 
          for(k=0;k<3;k++)   
            AlignWise_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,INTERGENIC,k) = AlignWise_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,CODON,k);   
          AlignWise_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,INTERGENIC,6) = AlignWise_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,CODON,6); 
          AlignWise_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,INTERGENIC,3) = i;  
          AlignWise_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,INTERGENIC,4) = j;  
          AlignWise_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,INTERGENIC,5) = CODON;  
          }  




        /* Finished calculating state CODON */ 


        /* For state INTRON_0 */ 
        /* setting first movement to score */ 
        score = AlignWise_DC_OPT_SHADOW_MATRIX(mat,i-0,j-8,CODON) + (mat->align->splice5_forward[(j-7)]+mat->intronopen) + (0);  
        /* assign local shadown pointer */ 
        localsp = &(AlignWise_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 8,CODON,0)); 
        /* From state INTRON_0 to state INTRON_0 */ 
        temp = AlignWise_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,INTRON_0) + 0 +(0);    
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(AlignWise_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,INTRON_0,0));    
          }  


        /* Ok - finished max calculation for INTRON_0 */ 
        /* Add any movement independant score and put away */ 
        /* Actually, already done inside scores */ 
         AlignWise_DC_OPT_SHADOW_MATRIX(mat,i,j,INTRON_0) = score;   
        for(k=0;k<7;k++) 
          AlignWise_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,INTRON_0,k) = localsp[k];    
        /* Now figure out if any specials need this score */ 


        /* Finished calculating state INTRON_0 */ 


        /* For state INTRON_1 */ 
        /* setting first movement to score */ 
        score = AlignWise_DC_OPT_SHADOW_MATRIX(mat,i-0,j-9,CODON) + (mat->align->splice5_forward[(j-7)]+mat->intronopen) + (0);  
        /* assign local shadown pointer */ 
        localsp = &(AlignWise_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 9,CODON,0)); 
        /* From state INTRON_1 to state INTRON_1 */ 
        temp = AlignWise_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,INTRON_1) + 0 +(0);    
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(AlignWise_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,INTRON_1,0));    
          }  


        /* Ok - finished max calculation for INTRON_1 */ 
        /* Add any movement independant score and put away */ 
        /* Actually, already done inside scores */ 
         AlignWise_DC_OPT_SHADOW_MATRIX(mat,i,j,INTRON_1) = score;   
        for(k=0;k<7;k++) 
          AlignWise_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,INTRON_1,k) = localsp[k];    
        /* Now figure out if any specials need this score */ 


        /* Finished calculating state INTRON_1 */ 


        /* For state INTRON_2 */ 
        /* setting first movement to score */ 
        score = AlignWise_DC_OPT_SHADOW_MATRIX(mat,i-0,j-10,CODON) + (mat->align->splice5_forward[(j-7)]+mat->intronopen) + (0);     
        /* assign local shadown pointer */ 
        localsp = &(AlignWise_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 10,CODON,0));    
        /* From state INTRON_2 to state INTRON_2 */ 
        temp = AlignWise_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,INTRON_2) + 0 +(0);    
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(AlignWise_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,INTRON_2,0));    
          }  


        /* Ok - finished max calculation for INTRON_2 */ 
        /* Add any movement independant score and put away */ 
        /* Actually, already done inside scores */ 
         AlignWise_DC_OPT_SHADOW_MATRIX(mat,i,j,INTRON_2) = score;   
        for(k=0;k<7;k++) 
          AlignWise_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,INTRON_2,k) = localsp[k];    
        /* Now figure out if any specials need this score */ 


        /* Finished calculating state INTRON_2 */ 


        /* For state REVERSE_CODON */ 
        /* setting first movement to score */ 
        score = AlignWise_DC_OPT_SHADOW_MATRIX(mat,i-0,j-3,REVERSE_CODON) + mat->align->reverse_coding[j] + (0);     
        /* assign local shadown pointer */ 
        localsp = &(AlignWise_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 3,REVERSE_CODON,0)); 
        /* From state INTERGENIC to state REVERSE_CODON */ 
        temp = AlignWise_DC_OPT_SHADOW_SPECIAL(mat,i-0,j-3,INTERGENIC) + mat->align->reverse_coding[j] + (0);    
        if( temp  > score )  {  
          score = temp;  
          /* This state [INTERGENIC] is a special for REVERSE_CODON... push top shadow pointers here */ 
          localshadow[0]= i; 
          localshadow[1]= j; 
          localshadow[2]= REVERSE_CODON; 
          localshadow[3]= (-1);  
          localshadow[4]= (-1);  
          localshadow[5]= (-1);  
          localshadow[6]= score; 
          localsp = localshadow; 
          }  
        /* From state REVERSE_INTRON_0 to state REVERSE_CODON */ 
        temp = AlignWise_DC_OPT_SHADOW_MATRIX(mat,i-0,j-6,REVERSE_INTRON_0) + mat->align->splice5_reverse[(j-3)] +(0);   
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(AlignWise_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 6,REVERSE_INTRON_0,0));    
          }  
        /* From state REVERSE_INTRON_1 to state REVERSE_CODON */ 
        temp = AlignWise_DC_OPT_SHADOW_MATRIX(mat,i-0,j-5,REVERSE_INTRON_1) + mat->align->splice5_reverse[(j-2)] +(0);   
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(AlignWise_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 5,REVERSE_INTRON_1,0));    
          }  
        /* From state REVERSE_INTRON_2 to state REVERSE_CODON */ 
        temp = AlignWise_DC_OPT_SHADOW_MATRIX(mat,i-0,j-4,REVERSE_INTRON_2) + mat->align->splice5_reverse[(j-1)] +(0);   
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(AlignWise_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 4,REVERSE_INTRON_2,0));    
          }  


        /* Ok - finished max calculation for REVERSE_CODON */ 
        /* Add any movement independant score and put away */ 
        /* Actually, already done inside scores */ 
         AlignWise_DC_OPT_SHADOW_MATRIX(mat,i,j,REVERSE_CODON) = score;  
        for(k=0;k<7;k++) 
          AlignWise_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,REVERSE_CODON,k) = localsp[k];   
        /* Now figure out if any specials need this score */ 


        /* state REVERSE_CODON is a source for special INTERGENIC */ 
        temp = score + (mat->geneopen) + (0) ;   
        if( temp > AlignWise_DC_OPT_SHADOW_SPECIAL(mat,i,j,INTERGENIC) )     {  
          AlignWise_DC_OPT_SHADOW_SPECIAL(mat,i,j,INTERGENIC) = temp;    
          /* Have to push only bottem half of system here */ 
          for(k=0;k<3;k++)   
            AlignWise_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,INTERGENIC,k) = AlignWise_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,REVERSE_CODON,k);   
          AlignWise_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,INTERGENIC,6) = AlignWise_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,REVERSE_CODON,6); 
          AlignWise_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,INTERGENIC,3) = i;  
          AlignWise_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,INTERGENIC,4) = j;  
          AlignWise_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,INTERGENIC,5) = REVERSE_CODON;  
          }  




        /* Finished calculating state REVERSE_CODON */ 


        /* For state REVERSE_INTRON_0 */ 
        /* setting first movement to score */ 
        score = AlignWise_DC_OPT_SHADOW_MATRIX(mat,i-0,j-8,REVERSE_CODON) + (mat->align->splice3_reverse[(j-7)]+mat->intronopen) + (0);  
        /* assign local shadown pointer */ 
        localsp = &(AlignWise_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 8,REVERSE_CODON,0)); 
        /* From state REVERSE_INTRON_0 to state REVERSE_INTRON_0 */ 
        temp = AlignWise_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,REVERSE_INTRON_0) + 0 +(0);    
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(AlignWise_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,REVERSE_INTRON_0,0));    
          }  


        /* Ok - finished max calculation for REVERSE_INTRON_0 */ 
        /* Add any movement independant score and put away */ 
        /* Actually, already done inside scores */ 
         AlignWise_DC_OPT_SHADOW_MATRIX(mat,i,j,REVERSE_INTRON_0) = score;   
        for(k=0;k<7;k++) 
          AlignWise_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,REVERSE_INTRON_0,k) = localsp[k];    
        /* Now figure out if any specials need this score */ 


        /* Finished calculating state REVERSE_INTRON_0 */ 


        /* For state REVERSE_INTRON_1 */ 
        /* setting first movement to score */ 
        score = AlignWise_DC_OPT_SHADOW_MATRIX(mat,i-0,j-9,REVERSE_CODON) + (mat->align->splice3_reverse[(j-7)]+mat->intronopen) + (0);  
        /* assign local shadown pointer */ 
        localsp = &(AlignWise_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 9,REVERSE_CODON,0)); 
        /* From state REVERSE_INTRON_1 to state REVERSE_INTRON_1 */ 
        temp = AlignWise_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,REVERSE_INTRON_1) + 0 +(0);    
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(AlignWise_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,REVERSE_INTRON_1,0));    
          }  


        /* Ok - finished max calculation for REVERSE_INTRON_1 */ 
        /* Add any movement independant score and put away */ 
        /* Actually, already done inside scores */ 
         AlignWise_DC_OPT_SHADOW_MATRIX(mat,i,j,REVERSE_INTRON_1) = score;   
        for(k=0;k<7;k++) 
          AlignWise_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,REVERSE_INTRON_1,k) = localsp[k];    
        /* Now figure out if any specials need this score */ 


        /* Finished calculating state REVERSE_INTRON_1 */ 


        /* For state REVERSE_INTRON_2 */ 
        /* setting first movement to score */ 
        score = AlignWise_DC_OPT_SHADOW_MATRIX(mat,i-0,j-10,REVERSE_CODON) + (mat->align->splice3_reverse[(j-7)]+mat->intronopen) + (0);     
        /* assign local shadown pointer */ 
        localsp = &(AlignWise_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 10,REVERSE_CODON,0));    
        /* From state REVERSE_INTRON_2 to state REVERSE_INTRON_2 */ 
        temp = AlignWise_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,REVERSE_INTRON_2) + 0 +(0);    
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(AlignWise_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,REVERSE_INTRON_2,0));    
          }  


        /* Ok - finished max calculation for REVERSE_INTRON_2 */ 
        /* Add any movement independant score and put away */ 
        /* Actually, already done inside scores */ 
         AlignWise_DC_OPT_SHADOW_MATRIX(mat,i,j,REVERSE_INTRON_2) = score;   
        for(k=0;k<7;k++) 
          AlignWise_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,REVERSE_INTRON_2,k) = localsp[k];    
        /* Now figure out if any specials need this score */ 


        /* Finished calculating state REVERSE_INTRON_2 */ 


        } /* end of for each i position in strip */ 


      /* Special state START has no special to special movements */ 


      /* Special state INTERGENIC has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = AlignWise_DC_OPT_SHADOW_SPECIAL(mat,0,j,INTERGENIC);   


      /* Source START is a special source for INTERGENIC */ 
      temp = AlignWise_DC_OPT_SHADOW_SPECIAL(mat,0,j - 1,START) + (0) + (0);     
      if( temp > score ) {  
        score = temp;    
        /* Also got to propagate shadows  */ 
        for(k=0;k<7;k++) 
          AlignWise_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,INTERGENIC,k) = AlignWise_DC_OPT_SHADOW_SPECIAL_SP(mat,i - 0,j - 1,START,k);    
        }  


      /* Source INTERGENIC is a special source for INTERGENIC */ 
      temp = AlignWise_DC_OPT_SHADOW_SPECIAL(mat,0,j - 1,INTERGENIC) + (0) + (0);    
      if( temp > score ) {  
        score = temp;    
        /* Also got to propagate shadows  */ 
        for(k=0;k<7;k++) 
          AlignWise_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,INTERGENIC,k) = AlignWise_DC_OPT_SHADOW_SPECIAL_SP(mat,i - 0,j - 1,INTERGENIC,k);   
        }  


      /* Source CODON for state INTERGENIC is not special... already calculated */ 
      /* Source REVERSE_CODON for state INTERGENIC is not special... already calculated */ 
      /* Put back score... (now updated!) */ 
      AlignWise_DC_OPT_SHADOW_SPECIAL(mat,0,j,INTERGENIC) = score;   
      /* Finished updating state INTERGENIC */ 




      /* Special state END has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = AlignWise_DC_OPT_SHADOW_SPECIAL(mat,0,j,END);  


      /* Source INTERGENIC is a special source for END */ 
      temp = AlignWise_DC_OPT_SHADOW_SPECIAL(mat,0,j - 1,INTERGENIC) + (0) + (0);    
      if( temp > score ) {  
        score = temp;    
        /* Also got to propagate shadows  */ 
        for(k=0;k<7;k++) 
          AlignWise_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,k) = AlignWise_DC_OPT_SHADOW_SPECIAL_SP(mat,i - 0,j - 1,INTERGENIC,k);  
        }  


      /* Put back score... (now updated!) */ 
      AlignWise_DC_OPT_SHADOW_SPECIAL(mat,0,j,END) = score;  
      /* Finished updating state END */ 


      } /* end of for each j strip */ 
    free(score_pointers);    
    free(shadow_pointers);   
    return TRUE;     
}    


/* Function:  init_start_end_linear_AlignWise(mat)
 *
 * Descrip: No Description
 *
 * Arg:        mat [UNKN ] Undocumented argument [AlignWise *]
 *
 */
void init_start_end_linear_AlignWise(AlignWise * mat) 
{
    register int i;  
    register int j;  
    for(j=0;j<12;j++)    {  
      for(i=(-0);i<mat->model->len;i++)  {  
        AlignWise_DC_SHADOW_MATRIX(mat,i,j,CODON) = NEGI;    
        AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,CODON,0) = (-1);   
        AlignWise_DC_SHADOW_MATRIX(mat,i,j,INTRON_0) = NEGI; 
        AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_0,0) = (-1);    
        AlignWise_DC_SHADOW_MATRIX(mat,i,j,INTRON_1) = NEGI; 
        AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_1,0) = (-1);    
        AlignWise_DC_SHADOW_MATRIX(mat,i,j,INTRON_2) = NEGI; 
        AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_2,0) = (-1);    
        AlignWise_DC_SHADOW_MATRIX(mat,i,j,REVERSE_CODON) = NEGI;    
        AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_CODON,0) = (-1);   
        AlignWise_DC_SHADOW_MATRIX(mat,i,j,REVERSE_INTRON_0) = NEGI; 
        AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_INTRON_0,0) = (-1);    
        AlignWise_DC_SHADOW_MATRIX(mat,i,j,REVERSE_INTRON_1) = NEGI; 
        AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_INTRON_1,0) = (-1);    
        AlignWise_DC_SHADOW_MATRIX(mat,i,j,REVERSE_INTRON_2) = NEGI; 
        AlignWise_DC_SHADOW_MATRIX_SP(mat,i,j,REVERSE_INTRON_2,0) = (-1);    
        }  
      }  


    for(j=(-10);j<mat->align->len;j++)   {  
      AlignWise_DC_SHADOW_SPECIAL(mat,0,j,START) = 0;    
      AlignWise_DC_SHADOW_SPECIAL_SP(mat,0,j,START,0) = j;   
      AlignWise_DC_SHADOW_SPECIAL(mat,0,j,INTERGENIC) = NEGI;    
      AlignWise_DC_SHADOW_SPECIAL_SP(mat,0,j,INTERGENIC,0) = (-1);   
      AlignWise_DC_SHADOW_SPECIAL(mat,0,j,END) = NEGI;   
      AlignWise_DC_SHADOW_SPECIAL_SP(mat,0,j,END,0) = (-1);  
      }  


    return;  
}    


/* Function:  convert_PackAln_to_AlnBlock_AlignWise(pal)
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
AlnBlock * convert_PackAln_to_AlnBlock_AlignWise(PackAln * pal) 
{
    AlnConvertSet * acs; 
    AlnBlock * alb;  


    acs = AlnConvertSet_AlignWise(); 
    alb = AlnBlock_from_PackAln(acs,pal);    
    free_AlnConvertSet(acs); 
    return alb;  
}    


 static char * query_label[] = { "EXON_STATE","INTRON_STATE","REV_EXON_STATE","REV_INTRON_STATE","INTERGENIC","END" };   
/* Function:  AlnConvertSet_AlignWise(void)
 *
 * Descrip: No Description
 *
 *
 * Return [UNKN ]  Undocumented return value [AlnConvertSet *]
 *
 */
 static char * target_label[] = { "CODON","3SS_PHASE_0","3SS_PHASE_1","3SS_PHASE_2","5SS_PHASE_0","INTRON","5SS_PHASE_1","5SS_PHASE_2","REV_CODON","REV_5SS_PHASE_0","REV_5SS_PHASE_1","REV_5SS_PHASE_2","REV_3SS_PHASE_0","REV_INTRON","REV_3SS_PHASE_1","REV_3SS_PHASE_2","SEQUENCE","STOP","END" };   
AlnConvertSet * AlnConvertSet_AlignWise(void) 
{
    AlnConvertUnit * acu;    
    AlnConvertSet  * out;    


    out = AlnConvertSet_alloc_std(); 


    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = CODON; 
    acu->state2 = CODON;     
    acu->offi = 0;   
    acu->offj = 3;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[0];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = INTERGENIC + 8;    
    acu->is_from_special = TRUE; 
    acu->state2 = CODON;     
    acu->offi = (-1);    
    acu->offj = 3;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[0];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = INTRON_0;  
    acu->state2 = CODON;     
    acu->offi = 0;   
    acu->offj = 6;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[1];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = INTRON_1;  
    acu->state2 = CODON;     
    acu->offi = 0;   
    acu->offj = 5;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[2];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = INTRON_2;  
    acu->state2 = CODON;     
    acu->offi = 0;   
    acu->offj = 4;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[3];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = CODON; 
    acu->state2 = INTRON_0;  
    acu->offi = 0;   
    acu->offj = 8;   
    acu->label1 = query_label[1];    
    acu->label2 = target_label[4];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = INTRON_0;  
    acu->state2 = INTRON_0;  
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[1];    
    acu->label2 = target_label[5];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = CODON; 
    acu->state2 = INTRON_1;  
    acu->offi = 0;   
    acu->offj = 9;   
    acu->label1 = query_label[1];    
    acu->label2 = target_label[6];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = INTRON_1;  
    acu->state2 = INTRON_1;  
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[1];    
    acu->label2 = target_label[5];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = CODON; 
    acu->state2 = INTRON_2;  
    acu->offi = 0;   
    acu->offj = 10;  
    acu->label1 = query_label[1];    
    acu->label2 = target_label[7];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = INTRON_2;  
    acu->state2 = INTRON_2;  
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[1];    
    acu->label2 = target_label[5];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = REVERSE_CODON; 
    acu->state2 = REVERSE_CODON;     
    acu->offi = 0;   
    acu->offj = 3;   
    acu->label1 = query_label[2];    
    acu->label2 = target_label[8];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = INTERGENIC + 8;    
    acu->is_from_special = TRUE; 
    acu->state2 = REVERSE_CODON;     
    acu->offi = (-1);    
    acu->offj = 3;   
    acu->label1 = query_label[2];    
    acu->label2 = target_label[8];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = REVERSE_INTRON_0;  
    acu->state2 = REVERSE_CODON;     
    acu->offi = 0;   
    acu->offj = 6;   
    acu->label1 = query_label[2];    
    acu->label2 = target_label[9];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = REVERSE_INTRON_1;  
    acu->state2 = REVERSE_CODON;     
    acu->offi = 0;   
    acu->offj = 5;   
    acu->label1 = query_label[2];    
    acu->label2 = target_label[10];  
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = REVERSE_INTRON_2;  
    acu->state2 = REVERSE_CODON;     
    acu->offi = 0;   
    acu->offj = 4;   
    acu->label1 = query_label[2];    
    acu->label2 = target_label[11];  
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = REVERSE_CODON; 
    acu->state2 = REVERSE_INTRON_0;  
    acu->offi = 0;   
    acu->offj = 8;   
    acu->label1 = query_label[3];    
    acu->label2 = target_label[12];  
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = REVERSE_INTRON_0;  
    acu->state2 = REVERSE_INTRON_0;  
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[3];    
    acu->label2 = target_label[13];  
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = REVERSE_CODON; 
    acu->state2 = REVERSE_INTRON_1;  
    acu->offi = 0;   
    acu->offj = 9;   
    acu->label1 = query_label[3];    
    acu->label2 = target_label[14];  
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = REVERSE_INTRON_1;  
    acu->state2 = REVERSE_INTRON_1;  
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[3];    
    acu->label2 = target_label[13];  
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = REVERSE_CODON; 
    acu->state2 = REVERSE_INTRON_2;  
    acu->offi = 0;   
    acu->offj = 10;  
    acu->label1 = query_label[3];    
    acu->label2 = target_label[15];  
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = REVERSE_INTRON_2;  
    acu->state2 = REVERSE_INTRON_2;  
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[3];    
    acu->label2 = target_label[13];  
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = START + 8; 
    acu->state2 = INTERGENIC + 8;    
    acu->offi = (-1);    
    acu->offj = 1;   
    acu->label1 = query_label[4];    
    acu->label2 = target_label[16];  
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = INTERGENIC + 8;    
    acu->state2 = INTERGENIC + 8;    
    acu->offi = (-1);    
    acu->offj = 1;   
    acu->label1 = query_label[4];    
    acu->label2 = target_label[16];  
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = CODON; 
    acu->state2 = INTERGENIC + 8;    
    acu->offi = (-1);    
    acu->offj = 0;   
    acu->label1 = query_label[4];    
    acu->label2 = target_label[17];  
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = REVERSE_CODON; 
    acu->state2 = INTERGENIC + 8;    
    acu->offi = (-1);    
    acu->offj = 0;   
    acu->label1 = query_label[4];    
    acu->label2 = target_label[8];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = INTERGENIC + 8;    
    acu->state2 = END + 8;   
    acu->offi = (-1);    
    acu->offj = 1;   
    acu->label1 = query_label[5];    
    acu->label2 = target_label[18];  
    add_collapse_label_AlnConvertSet(out,"INTRON_STATE","INTRON");   
    add_collapse_label_AlnConvertSet(out,"REV_INTRON_STATE","REV_INTRON");   
    add_collapse_label_AlnConvertSet(out,"INTERGENIC","SEQUENCE");   
    return out;  
}    


/* Function:  PackAln_read_Expl_AlignWise(mat)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [AlignWise *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * PackAln_read_Expl_AlignWise(AlignWise * mat) 
{
    AlignWise_access_func_holder holder;     


    holder.access_main    = AlignWise_explicit_access_main;  
    holder.access_special = AlignWise_explicit_access_special;   
    return PackAln_read_generic_AlignWise(mat,holder);   
}    


/* Function:  AlignWise_explicit_access_main(mat,i,j,state)
 *
 * Descrip: No Description
 *
 * Arg:          mat [UNKN ] Undocumented argument [AlignWise *]
 * Arg:            i [UNKN ] Undocumented argument [int]
 * Arg:            j [UNKN ] Undocumented argument [int]
 * Arg:        state [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int AlignWise_explicit_access_main(AlignWise * mat,int i,int j,int state) 
{
    return AlignWise_EXPL_MATRIX(mat,i,j,state); 
}    


/* Function:  AlignWise_explicit_access_special(mat,i,j,state)
 *
 * Descrip: No Description
 *
 * Arg:          mat [UNKN ] Undocumented argument [AlignWise *]
 * Arg:            i [UNKN ] Undocumented argument [int]
 * Arg:            j [UNKN ] Undocumented argument [int]
 * Arg:        state [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int AlignWise_explicit_access_special(AlignWise * mat,int i,int j,int state) 
{
    return AlignWise_EXPL_SPECIAL(mat,i,j,state);    
}    


/* Function:  PackAln_read_generic_AlignWise(mat,h)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [AlignWise *]
 * Arg:          h [UNKN ] Undocumented argument [AlignWise_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * PackAln_read_generic_AlignWise(AlignWise * mat,AlignWise_access_func_holder h) 
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


    out->score =  find_end_AlignWise(mat,&i,&j,&state,&isspecial,h); 


    /* Add final end transition (at the moment we have not got the score! */ 
    if( (pau= PackAlnUnit_alloc()) == NULL  || add_PackAln(out,pau) == FALSE )   {  
      warn("Failed the first PackAlnUnit alloc, %d length of Alignment in AlignWise_basic_read, returning a mess.(Sorry!)",out->len);    
      return out;    
      }  


    /* Put in positions for end trans. Remember that coordinates in C style */ 
    pau->i = i;  
    pau->j = j;  
    if( isspecial != TRUE)   
      pau->state = state;    
    else pau->state = state + 8;     
    prev=pau;    
    while( state != START || isspecial != TRUE)  { /*while state != START*/ 


      if( isspecial == TRUE )    
        max_calc_special_AlignWise(mat,i,j,state,isspecial,&i,&j,&state,&isspecial,&cellscore,h);    
      else   
        max_calc_AlignWise(mat,i,j,state,isspecial,&i,&j,&state,&isspecial,&cellscore,h);    
      if(i == AlignWise_READ_OFF_ERROR || j == AlignWise_READ_OFF_ERROR || state == AlignWise_READ_OFF_ERROR )   {  
        warn("Problem - hit bad read off system, exiting now");  
        break;   
        }  
      if( (pau= PackAlnUnit_alloc()) == NULL  || add_PackAln(out,pau) == FALSE ) {  
        warn("Failed a PackAlnUnit alloc, %d length of Alignment in AlignWise_basic_read, returning partial alignment",out->len);    
        break;   
        }  


      /* Put in positions for block. Remember that coordinates in C style */ 
      pau->i = i;    
      pau->j = j;    
      if( isspecial != TRUE)     
        pau->state = state;  
      else pau->state = state + 8;   
      prev->score = cellscore;   
      prev = pau;    
      } /* end of while state != START */ 


    invert_PackAln(out); 
    return out;  
}    


/* Function:  find_end_AlignWise(mat,ri,rj,state,isspecial,h)
 *
 * Descrip: No Description
 *
 * Arg:              mat [UNKN ] Undocumented argument [AlignWise *]
 * Arg:               ri [UNKN ] Undocumented argument [int *]
 * Arg:               rj [UNKN ] Undocumented argument [int *]
 * Arg:            state [UNKN ] Undocumented argument [int *]
 * Arg:        isspecial [UNKN ] Undocumented argument [boolean *]
 * Arg:                h [UNKN ] Undocumented argument [AlignWise_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int find_end_AlignWise(AlignWise * mat,int * ri,int * rj,int * state,boolean * isspecial,AlignWise_access_func_holder h) 
{
    int j;   
    int max; 
    int maxj;    
    int temp;    


    max = (*h.access_special)(mat,0,mat->align->len-1,END);  
    maxj = mat->align->len-1;    
    for(j= mat->align->len-2 ;j >= 0 ;j--)   {  
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


/* Function:  AlignWise_debug_show_matrix(mat,starti,stopi,startj,stopj,ofp)
 *
 * Descrip: No Description
 *
 * Arg:           mat [UNKN ] Undocumented argument [AlignWise *]
 * Arg:        starti [UNKN ] Undocumented argument [int]
 * Arg:         stopi [UNKN ] Undocumented argument [int]
 * Arg:        startj [UNKN ] Undocumented argument [int]
 * Arg:         stopj [UNKN ] Undocumented argument [int]
 * Arg:           ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
void AlignWise_debug_show_matrix(AlignWise * mat,int starti,int stopi,int startj,int stopj,FILE * ofp) 
{
    register int i;  
    register int j;  


    for(i=starti;i<stopi && i < mat->model->len;i++) {  
      for(j=startj;j<stopj && j < mat->align->len;j++)   {  
        fprintf(ofp,"Cell [%d - %d]\n",i,j);     
        fprintf(ofp,"State CODON %d\n",AlignWise_EXPL_MATRIX(mat,i,j,CODON));    
        fprintf(ofp,"State INTRON_0 %d\n",AlignWise_EXPL_MATRIX(mat,i,j,INTRON_0));  
        fprintf(ofp,"State INTRON_1 %d\n",AlignWise_EXPL_MATRIX(mat,i,j,INTRON_1));  
        fprintf(ofp,"State INTRON_2 %d\n",AlignWise_EXPL_MATRIX(mat,i,j,INTRON_2));  
        fprintf(ofp,"State REVERSE_CODON %d\n",AlignWise_EXPL_MATRIX(mat,i,j,REVERSE_CODON));    
        fprintf(ofp,"State REVERSE_INTRON_0 %d\n",AlignWise_EXPL_MATRIX(mat,i,j,REVERSE_INTRON_0));  
        fprintf(ofp,"State REVERSE_INTRON_1 %d\n",AlignWise_EXPL_MATRIX(mat,i,j,REVERSE_INTRON_1));  
        fprintf(ofp,"State REVERSE_INTRON_2 %d\n",AlignWise_EXPL_MATRIX(mat,i,j,REVERSE_INTRON_2));  
        fprintf(ofp,"\n\n"); 
        }  
      }  


}    


/* Function:  max_calc_AlignWise(mat,i,j,state,isspecial,reti,retj,retstate,retspecial,cellscore,h)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [AlignWise *]
 * Arg:                 i [UNKN ] Undocumented argument [int]
 * Arg:                 j [UNKN ] Undocumented argument [int]
 * Arg:             state [UNKN ] Undocumented argument [int]
 * Arg:         isspecial [UNKN ] Undocumented argument [boolean]
 * Arg:              reti [UNKN ] Undocumented argument [int *]
 * Arg:              retj [UNKN ] Undocumented argument [int *]
 * Arg:          retstate [UNKN ] Undocumented argument [int *]
 * Arg:        retspecial [UNKN ] Undocumented argument [boolean *]
 * Arg:         cellscore [UNKN ] Undocumented argument [int *]
 * Arg:                 h [UNKN ] Undocumented argument [AlignWise_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int max_calc_AlignWise(AlignWise * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,AlignWise_access_func_holder h) 
{
    register int temp;   
    register int cscore; 


    *reti = (*retj) = (*retstate) = AlignWise_READ_OFF_ERROR;    


    if( i < 0 || j < 0 || i > mat->model->len || j > mat->align->len)    {  
      warn("In AlignWise matrix special read off - out of bounds on matrix [i,j is %d,%d state %d in standard matrix]",i,j,state);   
      return -1;     
      }  


    /* Then you have to select the correct switch statement to figure out the readoff      */ 
    /* Somewhat odd - reverse the order of calculation and return as soon as it is correct */ 
    cscore = (*h.access_main)(mat,i,j,state);    
    switch(state)    { /*Switch state */ 
      case CODON :   
        temp = cscore - (mat->align->splice3_forward[(j-1)]) -  (0); 
        if( temp == (*h.access_main)(mat,i - 0,j - 4,INTRON_2) ) {  
          *reti = i - 0; 
          *retj = j - 4; 
          *retstate = INTRON_2;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-4,INTRON_2);    
            }  
          return (*h.access_main)(mat,i - 0,j - 4,INTRON_2);     
          }  
        temp = cscore - (mat->align->splice3_forward[(j-2)]) -  (0); 
        if( temp == (*h.access_main)(mat,i - 0,j - 5,INTRON_1) ) {  
          *reti = i - 0; 
          *retj = j - 5; 
          *retstate = INTRON_1;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-5,INTRON_1);    
            }  
          return (*h.access_main)(mat,i - 0,j - 5,INTRON_1);     
          }  
        temp = cscore - ((mat->align->splice3_forward[(j-3)]+mat->align->forward_coding[j])) -  (0); 
        if( temp == (*h.access_main)(mat,i - 0,j - 6,INTRON_0) ) {  
          *reti = i - 0; 
          *retj = j - 6; 
          *retstate = INTRON_0;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-6,INTRON_0);    
            }  
          return (*h.access_main)(mat,i - 0,j - 6,INTRON_0);     
          }  
        temp = cscore - (mat->align->forward_coding[j]) -  (0);  
        if( temp == (*h.access_special)(mat,i - 0,j - 3,INTERGENIC) )    {  
          *reti = i - 0; 
          *retj = j - 3; 
          *retstate = INTERGENIC;    
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_special)(mat,i-0,j-3,INTERGENIC);   
            }  
          return (*h.access_main)(mat,i - 0,j - 3,INTERGENIC);   
          }  
        temp = cscore - (mat->align->forward_coding[j]) -  (0);  
        if( temp == (*h.access_main)(mat,i - 0,j - 3,CODON) )    {  
          *reti = i - 0; 
          *retj = j - 3; 
          *retstate = CODON; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-3,CODON);   
            }  
          return (*h.access_main)(mat,i - 0,j - 3,CODON);    
          }  
        warn("Major problem (!) - in AlignWise read off, position %d,%d state %d no source found!",i,j,state);   
        return (-1); 
      case INTRON_0 :    
        temp = cscore - (0) -  (0);  
        if( temp == (*h.access_main)(mat,i - 0,j - 1,INTRON_0) ) {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = INTRON_0;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,INTRON_0);    
            }  
          return (*h.access_main)(mat,i - 0,j - 1,INTRON_0);     
          }  
        temp = cscore - ((mat->align->splice5_forward[(j-7)]+mat->intronopen)) -  (0);   
        if( temp == (*h.access_main)(mat,i - 0,j - 8,CODON) )    {  
          *reti = i - 0; 
          *retj = j - 8; 
          *retstate = CODON; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-8,CODON);   
            }  
          return (*h.access_main)(mat,i - 0,j - 8,CODON);    
          }  
        warn("Major problem (!) - in AlignWise read off, position %d,%d state %d no source found!",i,j,state);   
        return (-1); 
      case INTRON_1 :    
        temp = cscore - (0) -  (0);  
        if( temp == (*h.access_main)(mat,i - 0,j - 1,INTRON_1) ) {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = INTRON_1;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,INTRON_1);    
            }  
          return (*h.access_main)(mat,i - 0,j - 1,INTRON_1);     
          }  
        temp = cscore - ((mat->align->splice5_forward[(j-7)]+mat->intronopen)) -  (0);   
        if( temp == (*h.access_main)(mat,i - 0,j - 9,CODON) )    {  
          *reti = i - 0; 
          *retj = j - 9; 
          *retstate = CODON; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-9,CODON);   
            }  
          return (*h.access_main)(mat,i - 0,j - 9,CODON);    
          }  
        warn("Major problem (!) - in AlignWise read off, position %d,%d state %d no source found!",i,j,state);   
        return (-1); 
      case INTRON_2 :    
        temp = cscore - (0) -  (0);  
        if( temp == (*h.access_main)(mat,i - 0,j - 1,INTRON_2) ) {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = INTRON_2;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,INTRON_2);    
            }  
          return (*h.access_main)(mat,i - 0,j - 1,INTRON_2);     
          }  
        temp = cscore - ((mat->align->splice5_forward[(j-7)]+mat->intronopen)) -  (0);   
        if( temp == (*h.access_main)(mat,i - 0,j - 10,CODON) )   {  
          *reti = i - 0; 
          *retj = j - 10;    
          *retstate = CODON; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-10,CODON);  
            }  
          return (*h.access_main)(mat,i - 0,j - 10,CODON);   
          }  
        warn("Major problem (!) - in AlignWise read off, position %d,%d state %d no source found!",i,j,state);   
        return (-1); 
      case REVERSE_CODON :   
        temp = cscore - (mat->align->splice5_reverse[(j-1)]) -  (0); 
        if( temp == (*h.access_main)(mat,i - 0,j - 4,REVERSE_INTRON_2) ) {  
          *reti = i - 0; 
          *retj = j - 4; 
          *retstate = REVERSE_INTRON_2;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-4,REVERSE_INTRON_2);    
            }  
          return (*h.access_main)(mat,i - 0,j - 4,REVERSE_INTRON_2);     
          }  
        temp = cscore - (mat->align->splice5_reverse[(j-2)]) -  (0); 
        if( temp == (*h.access_main)(mat,i - 0,j - 5,REVERSE_INTRON_1) ) {  
          *reti = i - 0; 
          *retj = j - 5; 
          *retstate = REVERSE_INTRON_1;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-5,REVERSE_INTRON_1);    
            }  
          return (*h.access_main)(mat,i - 0,j - 5,REVERSE_INTRON_1);     
          }  
        temp = cscore - (mat->align->splice5_reverse[(j-3)]) -  (0); 
        if( temp == (*h.access_main)(mat,i - 0,j - 6,REVERSE_INTRON_0) ) {  
          *reti = i - 0; 
          *retj = j - 6; 
          *retstate = REVERSE_INTRON_0;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-6,REVERSE_INTRON_0);    
            }  
          return (*h.access_main)(mat,i - 0,j - 6,REVERSE_INTRON_0);     
          }  
        temp = cscore - (mat->align->reverse_coding[j]) -  (0);  
        if( temp == (*h.access_special)(mat,i - 0,j - 3,INTERGENIC) )    {  
          *reti = i - 0; 
          *retj = j - 3; 
          *retstate = INTERGENIC;    
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_special)(mat,i-0,j-3,INTERGENIC);   
            }  
          return (*h.access_main)(mat,i - 0,j - 3,INTERGENIC);   
          }  
        temp = cscore - (mat->align->reverse_coding[j]) -  (0);  
        if( temp == (*h.access_main)(mat,i - 0,j - 3,REVERSE_CODON) )    {  
          *reti = i - 0; 
          *retj = j - 3; 
          *retstate = REVERSE_CODON; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-3,REVERSE_CODON);   
            }  
          return (*h.access_main)(mat,i - 0,j - 3,REVERSE_CODON);    
          }  
        warn("Major problem (!) - in AlignWise read off, position %d,%d state %d no source found!",i,j,state);   
        return (-1); 
      case REVERSE_INTRON_0 :    
        temp = cscore - (0) -  (0);  
        if( temp == (*h.access_main)(mat,i - 0,j - 1,REVERSE_INTRON_0) ) {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = REVERSE_INTRON_0;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,REVERSE_INTRON_0);    
            }  
          return (*h.access_main)(mat,i - 0,j - 1,REVERSE_INTRON_0);     
          }  
        temp = cscore - ((mat->align->splice3_reverse[(j-7)]+mat->intronopen)) -  (0);   
        if( temp == (*h.access_main)(mat,i - 0,j - 8,REVERSE_CODON) )    {  
          *reti = i - 0; 
          *retj = j - 8; 
          *retstate = REVERSE_CODON; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-8,REVERSE_CODON);   
            }  
          return (*h.access_main)(mat,i - 0,j - 8,REVERSE_CODON);    
          }  
        warn("Major problem (!) - in AlignWise read off, position %d,%d state %d no source found!",i,j,state);   
        return (-1); 
      case REVERSE_INTRON_1 :    
        temp = cscore - (0) -  (0);  
        if( temp == (*h.access_main)(mat,i - 0,j - 1,REVERSE_INTRON_1) ) {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = REVERSE_INTRON_1;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,REVERSE_INTRON_1);    
            }  
          return (*h.access_main)(mat,i - 0,j - 1,REVERSE_INTRON_1);     
          }  
        temp = cscore - ((mat->align->splice3_reverse[(j-7)]+mat->intronopen)) -  (0);   
        if( temp == (*h.access_main)(mat,i - 0,j - 9,REVERSE_CODON) )    {  
          *reti = i - 0; 
          *retj = j - 9; 
          *retstate = REVERSE_CODON; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-9,REVERSE_CODON);   
            }  
          return (*h.access_main)(mat,i - 0,j - 9,REVERSE_CODON);    
          }  
        warn("Major problem (!) - in AlignWise read off, position %d,%d state %d no source found!",i,j,state);   
        return (-1); 
      case REVERSE_INTRON_2 :    
        temp = cscore - (0) -  (0);  
        if( temp == (*h.access_main)(mat,i - 0,j - 1,REVERSE_INTRON_2) ) {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = REVERSE_INTRON_2;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,REVERSE_INTRON_2);    
            }  
          return (*h.access_main)(mat,i - 0,j - 1,REVERSE_INTRON_2);     
          }  
        temp = cscore - ((mat->align->splice3_reverse[(j-7)]+mat->intronopen)) -  (0);   
        if( temp == (*h.access_main)(mat,i - 0,j - 10,REVERSE_CODON) )   {  
          *reti = i - 0; 
          *retj = j - 10;    
          *retstate = REVERSE_CODON; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-10,REVERSE_CODON);  
            }  
          return (*h.access_main)(mat,i - 0,j - 10,REVERSE_CODON);   
          }  
        warn("Major problem (!) - in AlignWise read off, position %d,%d state %d no source found!",i,j,state);   
        return (-1); 
      default:   
        warn("Major problem (!) - in AlignWise read off, position %d,%d state %d no source found!",i,j,state);   
        return (-1); 
      } /* end of Switch state  */ 
}    


/* Function:  max_calc_special_AlignWise(mat,i,j,state,isspecial,reti,retj,retstate,retspecial,cellscore,h)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [AlignWise *]
 * Arg:                 i [UNKN ] Undocumented argument [int]
 * Arg:                 j [UNKN ] Undocumented argument [int]
 * Arg:             state [UNKN ] Undocumented argument [int]
 * Arg:         isspecial [UNKN ] Undocumented argument [boolean]
 * Arg:              reti [UNKN ] Undocumented argument [int *]
 * Arg:              retj [UNKN ] Undocumented argument [int *]
 * Arg:          retstate [UNKN ] Undocumented argument [int *]
 * Arg:        retspecial [UNKN ] Undocumented argument [boolean *]
 * Arg:         cellscore [UNKN ] Undocumented argument [int *]
 * Arg:                 h [UNKN ] Undocumented argument [AlignWise_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int max_calc_special_AlignWise(AlignWise * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,AlignWise_access_func_holder h) 
{
    register int temp;   
    register int cscore; 


    *reti = (*retj) = (*retstate) = AlignWise_READ_OFF_ERROR;    


    if( j < 0 || j > mat->align->len)    {  
      warn("In AlignWise matrix special read off - out of bounds on matrix [j is %d in special]",j); 
      return -1;     
      }  


    cscore = (*h.access_special)(mat,i,j,state); 
    switch(state)    { /*switch on special states*/ 
      case START :   
      case INTERGENIC :  
        /* source REVERSE_CODON is from main matrix */ 
        for(i= mat->model->len-1;i >= 0 ;i--)    { /*for i >= 0*/ 
          temp = cscore - (mat->geneopen) - (0);     
          if( temp == (*h.access_main)(mat,i - 0,j - 0,REVERSE_CODON) )  {  
            *reti = i - 0;   
            *retj = j - 0;   
            *retstate = REVERSE_CODON;   
            *retspecial = FALSE; 
            if( cellscore != NULL)   {  
              *cellscore = cscore - (*h.access_main)(mat,i-0,j-0,REVERSE_CODON);     
              }  
            return (*h.access_main)(mat,i - 0,j - 0,REVERSE_CODON) ;     
            }  
          } /* end of for i >= 0 */ 
        /* source CODON is from main matrix */ 
        for(i= mat->model->len-1;i >= 0 ;i--)    { /*for i >= 0*/ 
          temp = cscore - (mat->geneopen) - (0);     
          if( temp == (*h.access_main)(mat,i - 0,j - 0,CODON) )  {  
            *reti = i - 0;   
            *retj = j - 0;   
            *retstate = CODON;   
            *retspecial = FALSE; 
            if( cellscore != NULL)   {  
              *cellscore = cscore - (*h.access_main)(mat,i-0,j-0,CODON);     
              }  
            return (*h.access_main)(mat,i - 0,j - 0,CODON) ;     
            }  
          } /* end of for i >= 0 */ 
        /* source INTERGENIC is a special */ 
        temp = cscore - (0) - (0);   
        if( temp == (*h.access_special)(mat,i - 0,j - 1,INTERGENIC) )    {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = INTERGENIC;    
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_special)(mat,i-0,j-1,INTERGENIC);   
            }  
          return (*h.access_special)(mat,i - 0,j - 1,INTERGENIC) ;   
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
        /* source INTERGENIC is a special */ 
        temp = cscore - (0) - (0);   
        if( temp == (*h.access_special)(mat,i - 0,j - 1,INTERGENIC) )    {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = INTERGENIC;    
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_special)(mat,i-0,j-1,INTERGENIC);   
            }  
          return (*h.access_special)(mat,i - 0,j - 1,INTERGENIC) ;   
          }  
      default:   
        warn("Major problem (!) - in AlignWise read off, position %d,%d state %d no source found  dropped into default on source switch!",i,j,state);    
        return (-1); 
      } /* end of switch on special states */ 
}    


/* Function:  calculate_AlignWise(mat)
 *
 * Descrip:    This function calculates the AlignWise matrix when in explicit mode
 *             To allocate the matrix use /allocate_Expl_AlignWise
 *
 *
 * Arg:        mat [UNKN ] AlignWise which contains explicit basematrix memory [AlignWise *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean calculate_AlignWise(AlignWise * mat) 
{
    int i;   
    int j;   
    int leni;    
    int lenj;    
    int tot; 
    int num; 


    if( mat->basematrix->type != BASEMATRIX_TYPE_EXPLICIT )  {  
      warn("in calculate_AlignWise, passed a non Explicit matrix type, cannot calculate!");  
      return FALSE;  
      }  


    leni = mat->leni;    
    lenj = mat->lenj;    
    tot = leni * lenj;   
    num = 0; 


    start_reporting("AlignWise Matrix calculation: ");   
    for(j=0;j<lenj;j++)  {  
      auto int score;    
      auto int temp;     
      for(i=0;i<leni;i++)    {  
        if( num%1000 == 0 )  
          log_full_error(REPORT,0,"[%7d] Cells %2d%%%%",num,num*100/tot);    
        num++;   


        /* For state CODON */ 
        /* setting first movement to score */ 
        score = AlignWise_EXPL_MATRIX(mat,i-0,j-3,CODON) + mat->align->forward_coding[j];    
        /* From state INTERGENIC to state CODON */ 
        temp = AlignWise_EXPL_SPECIAL(mat,i-0,j-3,INTERGENIC) + mat->align->forward_coding[j];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTRON_0 to state CODON */ 
        temp = AlignWise_EXPL_MATRIX(mat,i-0,j-6,INTRON_0) + (mat->align->splice3_forward[(j-3)]+mat->align->forward_coding[j]);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTRON_1 to state CODON */ 
        temp = AlignWise_EXPL_MATRIX(mat,i-0,j-5,INTRON_1) + mat->align->splice3_forward[(j-2)];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTRON_2 to state CODON */ 
        temp = AlignWise_EXPL_MATRIX(mat,i-0,j-4,INTRON_2) + mat->align->splice3_forward[(j-1)];     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for CODON */ 
        /* Add any movement independant score and put away */ 
         AlignWise_EXPL_MATRIX(mat,i,j,CODON) = score;   


        /* state CODON is a source for special INTERGENIC */ 
        temp = score + (mat->geneopen) + (0) ;   
        if( temp > AlignWise_EXPL_SPECIAL(mat,i,j,INTERGENIC) )  {  
          AlignWise_EXPL_SPECIAL(mat,i,j,INTERGENIC) = temp;     
          }  




        /* Finished calculating state CODON */ 


        /* For state INTRON_0 */ 
        /* setting first movement to score */ 
        score = AlignWise_EXPL_MATRIX(mat,i-0,j-8,CODON) + (mat->align->splice5_forward[(j-7)]+mat->intronopen);     
        /* From state INTRON_0 to state INTRON_0 */ 
        temp = AlignWise_EXPL_MATRIX(mat,i-0,j-1,INTRON_0) + 0;  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for INTRON_0 */ 
        /* Add any movement independant score and put away */ 
         AlignWise_EXPL_MATRIX(mat,i,j,INTRON_0) = score;    


        /* Finished calculating state INTRON_0 */ 


        /* For state INTRON_1 */ 
        /* setting first movement to score */ 
        score = AlignWise_EXPL_MATRIX(mat,i-0,j-9,CODON) + (mat->align->splice5_forward[(j-7)]+mat->intronopen);     
        /* From state INTRON_1 to state INTRON_1 */ 
        temp = AlignWise_EXPL_MATRIX(mat,i-0,j-1,INTRON_1) + 0;  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for INTRON_1 */ 
        /* Add any movement independant score and put away */ 
         AlignWise_EXPL_MATRIX(mat,i,j,INTRON_1) = score;    


        /* Finished calculating state INTRON_1 */ 


        /* For state INTRON_2 */ 
        /* setting first movement to score */ 
        score = AlignWise_EXPL_MATRIX(mat,i-0,j-10,CODON) + (mat->align->splice5_forward[(j-7)]+mat->intronopen);    
        /* From state INTRON_2 to state INTRON_2 */ 
        temp = AlignWise_EXPL_MATRIX(mat,i-0,j-1,INTRON_2) + 0;  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for INTRON_2 */ 
        /* Add any movement independant score and put away */ 
         AlignWise_EXPL_MATRIX(mat,i,j,INTRON_2) = score;    


        /* Finished calculating state INTRON_2 */ 


        /* For state REVERSE_CODON */ 
        /* setting first movement to score */ 
        score = AlignWise_EXPL_MATRIX(mat,i-0,j-3,REVERSE_CODON) + mat->align->reverse_coding[j];    
        /* From state INTERGENIC to state REVERSE_CODON */ 
        temp = AlignWise_EXPL_SPECIAL(mat,i-0,j-3,INTERGENIC) + mat->align->reverse_coding[j];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state REVERSE_INTRON_0 to state REVERSE_CODON */ 
        temp = AlignWise_EXPL_MATRIX(mat,i-0,j-6,REVERSE_INTRON_0) + mat->align->splice5_reverse[(j-3)];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state REVERSE_INTRON_1 to state REVERSE_CODON */ 
        temp = AlignWise_EXPL_MATRIX(mat,i-0,j-5,REVERSE_INTRON_1) + mat->align->splice5_reverse[(j-2)];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state REVERSE_INTRON_2 to state REVERSE_CODON */ 
        temp = AlignWise_EXPL_MATRIX(mat,i-0,j-4,REVERSE_INTRON_2) + mat->align->splice5_reverse[(j-1)];     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for REVERSE_CODON */ 
        /* Add any movement independant score and put away */ 
         AlignWise_EXPL_MATRIX(mat,i,j,REVERSE_CODON) = score;   


        /* state REVERSE_CODON is a source for special INTERGENIC */ 
        temp = score + (mat->geneopen) + (0) ;   
        if( temp > AlignWise_EXPL_SPECIAL(mat,i,j,INTERGENIC) )  {  
          AlignWise_EXPL_SPECIAL(mat,i,j,INTERGENIC) = temp;     
          }  




        /* Finished calculating state REVERSE_CODON */ 


        /* For state REVERSE_INTRON_0 */ 
        /* setting first movement to score */ 
        score = AlignWise_EXPL_MATRIX(mat,i-0,j-8,REVERSE_CODON) + (mat->align->splice3_reverse[(j-7)]+mat->intronopen);     
        /* From state REVERSE_INTRON_0 to state REVERSE_INTRON_0 */ 
        temp = AlignWise_EXPL_MATRIX(mat,i-0,j-1,REVERSE_INTRON_0) + 0;  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for REVERSE_INTRON_0 */ 
        /* Add any movement independant score and put away */ 
         AlignWise_EXPL_MATRIX(mat,i,j,REVERSE_INTRON_0) = score;    


        /* Finished calculating state REVERSE_INTRON_0 */ 


        /* For state REVERSE_INTRON_1 */ 
        /* setting first movement to score */ 
        score = AlignWise_EXPL_MATRIX(mat,i-0,j-9,REVERSE_CODON) + (mat->align->splice3_reverse[(j-7)]+mat->intronopen);     
        /* From state REVERSE_INTRON_1 to state REVERSE_INTRON_1 */ 
        temp = AlignWise_EXPL_MATRIX(mat,i-0,j-1,REVERSE_INTRON_1) + 0;  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for REVERSE_INTRON_1 */ 
        /* Add any movement independant score and put away */ 
         AlignWise_EXPL_MATRIX(mat,i,j,REVERSE_INTRON_1) = score;    


        /* Finished calculating state REVERSE_INTRON_1 */ 


        /* For state REVERSE_INTRON_2 */ 
        /* setting first movement to score */ 
        score = AlignWise_EXPL_MATRIX(mat,i-0,j-10,REVERSE_CODON) + (mat->align->splice3_reverse[(j-7)]+mat->intronopen);    
        /* From state REVERSE_INTRON_2 to state REVERSE_INTRON_2 */ 
        temp = AlignWise_EXPL_MATRIX(mat,i-0,j-1,REVERSE_INTRON_2) + 0;  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for REVERSE_INTRON_2 */ 
        /* Add any movement independant score and put away */ 
         AlignWise_EXPL_MATRIX(mat,i,j,REVERSE_INTRON_2) = score;    


        /* Finished calculating state REVERSE_INTRON_2 */ 
        }  


      /* Special state START has no special to special movements */ 


      /* Special state INTERGENIC has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = AlignWise_EXPL_SPECIAL(mat,0,j,INTERGENIC);    


      /* Source START is a special source for INTERGENIC */ 
      temp = AlignWise_EXPL_SPECIAL(mat,0,j - 1,START) + (0) + (0);  
      if( temp > score ) 
        score = temp;    


      /* Source INTERGENIC is a special source for INTERGENIC */ 
      temp = AlignWise_EXPL_SPECIAL(mat,0,j - 1,INTERGENIC) + (0) + (0);     
      if( temp > score ) 
        score = temp;    


      /* Source CODON for state INTERGENIC is not special... already calculated */ 
      /* Source REVERSE_CODON for state INTERGENIC is not special... already calculated */ 
      /* Put back score... (now updated!) */ 
      AlignWise_EXPL_SPECIAL(mat,0,j,INTERGENIC) = score;    
      /* Finished updating state INTERGENIC */ 




      /* Special state END has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = AlignWise_EXPL_SPECIAL(mat,0,j,END);   


      /* Source INTERGENIC is a special source for END */ 
      temp = AlignWise_EXPL_SPECIAL(mat,0,j - 1,INTERGENIC) + (0) + (0);     
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      AlignWise_EXPL_SPECIAL(mat,0,j,END) = score;   
      /* Finished updating state END */ 


      }  
    stop_reporting();    
    return TRUE;     
}    


/* Function:  calculate_dpenv_AlignWise(mat,dpenv)
 *
 * Descrip:    This function calculates the AlignWise matrix when in explicit mode, subject to the envelope
 *
 *
 * Arg:          mat [UNKN ] AlignWise which contains explicit basematrix memory [AlignWise *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean calculate_dpenv_AlignWise(AlignWise * mat,DPEnvelope * dpenv) 
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
      warn("in calculate_AlignWise, passed a non Explicit matrix type, cannot calculate!");  
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


    for(j=startj-10;j<endj;j++)  {  
      for(i=0;i<mat->leni;i++)   {  
        AlignWise_EXPL_MATRIX(mat,i,j,CODON) = NEGI; 
        AlignWise_EXPL_MATRIX(mat,i,j,INTRON_0) = NEGI;  
        AlignWise_EXPL_MATRIX(mat,i,j,INTRON_1) = NEGI;  
        AlignWise_EXPL_MATRIX(mat,i,j,INTRON_2) = NEGI;  
        AlignWise_EXPL_MATRIX(mat,i,j,REVERSE_CODON) = NEGI; 
        AlignWise_EXPL_MATRIX(mat,i,j,REVERSE_INTRON_0) = NEGI;  
        AlignWise_EXPL_MATRIX(mat,i,j,REVERSE_INTRON_1) = NEGI;  
        AlignWise_EXPL_MATRIX(mat,i,j,REVERSE_INTRON_2) = NEGI;  
        }  
      }  
    for(j=-10;j<mat->lenj;j++)   {  
      AlignWise_EXPL_SPECIAL(mat,i,j,START) = 0; 
      AlignWise_EXPL_SPECIAL(mat,i,j,INTERGENIC) = NEGI; 
      AlignWise_EXPL_SPECIAL(mat,i,j,END) = NEGI;    
      }  


    start_reporting("AlignWise Matrix calculation: ");   
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
          AlignWise_EXPL_MATRIX(mat,i,j,CODON) = NEGI;   
          AlignWise_EXPL_MATRIX(mat,i,j,INTRON_0) = NEGI;    
          AlignWise_EXPL_MATRIX(mat,i,j,INTRON_1) = NEGI;    
          AlignWise_EXPL_MATRIX(mat,i,j,INTRON_2) = NEGI;    
          AlignWise_EXPL_MATRIX(mat,i,j,REVERSE_CODON) = NEGI;   
          AlignWise_EXPL_MATRIX(mat,i,j,REVERSE_INTRON_0) = NEGI;    
          AlignWise_EXPL_MATRIX(mat,i,j,REVERSE_INTRON_1) = NEGI;    
          AlignWise_EXPL_MATRIX(mat,i,j,REVERSE_INTRON_2) = NEGI;    
          continue;  
          }  


        if( num%1000 == 0 )  
          log_full_error(REPORT,0,"[%7d] Cells %2d%%%%",num,num*100/tot);    
        num++;   


        /* For state CODON */ 
        /* setting first movement to score */ 
        score = AlignWise_EXPL_MATRIX(mat,i-0,j-3,CODON) + mat->align->forward_coding[j];    
        /* From state INTERGENIC to state CODON */ 
        temp = AlignWise_EXPL_SPECIAL(mat,i-0,j-3,INTERGENIC) + mat->align->forward_coding[j];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTRON_0 to state CODON */ 
        temp = AlignWise_EXPL_MATRIX(mat,i-0,j-6,INTRON_0) + (mat->align->splice3_forward[(j-3)]+mat->align->forward_coding[j]);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTRON_1 to state CODON */ 
        temp = AlignWise_EXPL_MATRIX(mat,i-0,j-5,INTRON_1) + mat->align->splice3_forward[(j-2)];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTRON_2 to state CODON */ 
        temp = AlignWise_EXPL_MATRIX(mat,i-0,j-4,INTRON_2) + mat->align->splice3_forward[(j-1)];     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for CODON */ 
        /* Add any movement independant score and put away */ 
         AlignWise_EXPL_MATRIX(mat,i,j,CODON) = score;   


        /* state CODON is a source for special INTERGENIC */ 
        temp = score + (mat->geneopen) + (0) ;   
        if( temp > AlignWise_EXPL_SPECIAL(mat,i,j,INTERGENIC) )  {  
          AlignWise_EXPL_SPECIAL(mat,i,j,INTERGENIC) = temp;     
          }  




        /* Finished calculating state CODON */ 


        /* For state INTRON_0 */ 
        /* setting first movement to score */ 
        score = AlignWise_EXPL_MATRIX(mat,i-0,j-8,CODON) + (mat->align->splice5_forward[(j-7)]+mat->intronopen);     
        /* From state INTRON_0 to state INTRON_0 */ 
        temp = AlignWise_EXPL_MATRIX(mat,i-0,j-1,INTRON_0) + 0;  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for INTRON_0 */ 
        /* Add any movement independant score and put away */ 
         AlignWise_EXPL_MATRIX(mat,i,j,INTRON_0) = score;    


        /* Finished calculating state INTRON_0 */ 


        /* For state INTRON_1 */ 
        /* setting first movement to score */ 
        score = AlignWise_EXPL_MATRIX(mat,i-0,j-9,CODON) + (mat->align->splice5_forward[(j-7)]+mat->intronopen);     
        /* From state INTRON_1 to state INTRON_1 */ 
        temp = AlignWise_EXPL_MATRIX(mat,i-0,j-1,INTRON_1) + 0;  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for INTRON_1 */ 
        /* Add any movement independant score and put away */ 
         AlignWise_EXPL_MATRIX(mat,i,j,INTRON_1) = score;    


        /* Finished calculating state INTRON_1 */ 


        /* For state INTRON_2 */ 
        /* setting first movement to score */ 
        score = AlignWise_EXPL_MATRIX(mat,i-0,j-10,CODON) + (mat->align->splice5_forward[(j-7)]+mat->intronopen);    
        /* From state INTRON_2 to state INTRON_2 */ 
        temp = AlignWise_EXPL_MATRIX(mat,i-0,j-1,INTRON_2) + 0;  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for INTRON_2 */ 
        /* Add any movement independant score and put away */ 
         AlignWise_EXPL_MATRIX(mat,i,j,INTRON_2) = score;    


        /* Finished calculating state INTRON_2 */ 


        /* For state REVERSE_CODON */ 
        /* setting first movement to score */ 
        score = AlignWise_EXPL_MATRIX(mat,i-0,j-3,REVERSE_CODON) + mat->align->reverse_coding[j];    
        /* From state INTERGENIC to state REVERSE_CODON */ 
        temp = AlignWise_EXPL_SPECIAL(mat,i-0,j-3,INTERGENIC) + mat->align->reverse_coding[j];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state REVERSE_INTRON_0 to state REVERSE_CODON */ 
        temp = AlignWise_EXPL_MATRIX(mat,i-0,j-6,REVERSE_INTRON_0) + mat->align->splice5_reverse[(j-3)];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state REVERSE_INTRON_1 to state REVERSE_CODON */ 
        temp = AlignWise_EXPL_MATRIX(mat,i-0,j-5,REVERSE_INTRON_1) + mat->align->splice5_reverse[(j-2)];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state REVERSE_INTRON_2 to state REVERSE_CODON */ 
        temp = AlignWise_EXPL_MATRIX(mat,i-0,j-4,REVERSE_INTRON_2) + mat->align->splice5_reverse[(j-1)];     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for REVERSE_CODON */ 
        /* Add any movement independant score and put away */ 
         AlignWise_EXPL_MATRIX(mat,i,j,REVERSE_CODON) = score;   


        /* state REVERSE_CODON is a source for special INTERGENIC */ 
        temp = score + (mat->geneopen) + (0) ;   
        if( temp > AlignWise_EXPL_SPECIAL(mat,i,j,INTERGENIC) )  {  
          AlignWise_EXPL_SPECIAL(mat,i,j,INTERGENIC) = temp;     
          }  




        /* Finished calculating state REVERSE_CODON */ 


        /* For state REVERSE_INTRON_0 */ 
        /* setting first movement to score */ 
        score = AlignWise_EXPL_MATRIX(mat,i-0,j-8,REVERSE_CODON) + (mat->align->splice3_reverse[(j-7)]+mat->intronopen);     
        /* From state REVERSE_INTRON_0 to state REVERSE_INTRON_0 */ 
        temp = AlignWise_EXPL_MATRIX(mat,i-0,j-1,REVERSE_INTRON_0) + 0;  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for REVERSE_INTRON_0 */ 
        /* Add any movement independant score and put away */ 
         AlignWise_EXPL_MATRIX(mat,i,j,REVERSE_INTRON_0) = score;    


        /* Finished calculating state REVERSE_INTRON_0 */ 


        /* For state REVERSE_INTRON_1 */ 
        /* setting first movement to score */ 
        score = AlignWise_EXPL_MATRIX(mat,i-0,j-9,REVERSE_CODON) + (mat->align->splice3_reverse[(j-7)]+mat->intronopen);     
        /* From state REVERSE_INTRON_1 to state REVERSE_INTRON_1 */ 
        temp = AlignWise_EXPL_MATRIX(mat,i-0,j-1,REVERSE_INTRON_1) + 0;  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for REVERSE_INTRON_1 */ 
        /* Add any movement independant score and put away */ 
         AlignWise_EXPL_MATRIX(mat,i,j,REVERSE_INTRON_1) = score;    


        /* Finished calculating state REVERSE_INTRON_1 */ 


        /* For state REVERSE_INTRON_2 */ 
        /* setting first movement to score */ 
        score = AlignWise_EXPL_MATRIX(mat,i-0,j-10,REVERSE_CODON) + (mat->align->splice3_reverse[(j-7)]+mat->intronopen);    
        /* From state REVERSE_INTRON_2 to state REVERSE_INTRON_2 */ 
        temp = AlignWise_EXPL_MATRIX(mat,i-0,j-1,REVERSE_INTRON_2) + 0;  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for REVERSE_INTRON_2 */ 
        /* Add any movement independant score and put away */ 
         AlignWise_EXPL_MATRIX(mat,i,j,REVERSE_INTRON_2) = score;    


        /* Finished calculating state REVERSE_INTRON_2 */ 
        }  


      /* Special state START has no special to special movements */ 


      /* Special state INTERGENIC has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = AlignWise_EXPL_SPECIAL(mat,0,j,INTERGENIC);    


      /* Source START is a special source for INTERGENIC */ 
      temp = AlignWise_EXPL_SPECIAL(mat,0,j - 1,START) + (0) + (0);  
      if( temp > score ) 
        score = temp;    


      /* Source INTERGENIC is a special source for INTERGENIC */ 
      temp = AlignWise_EXPL_SPECIAL(mat,0,j - 1,INTERGENIC) + (0) + (0);     
      if( temp > score ) 
        score = temp;    


      /* Source CODON for state INTERGENIC is not special... already calculated */ 
      /* Source REVERSE_CODON for state INTERGENIC is not special... already calculated */ 
      /* Put back score... (now updated!) */ 
      AlignWise_EXPL_SPECIAL(mat,0,j,INTERGENIC) = score;    
      /* Finished updating state INTERGENIC */ 




      /* Special state END has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = AlignWise_EXPL_SPECIAL(mat,0,j,END);   


      /* Source INTERGENIC is a special source for END */ 
      temp = AlignWise_EXPL_SPECIAL(mat,0,j - 1,INTERGENIC) + (0) + (0);     
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      AlignWise_EXPL_SPECIAL(mat,0,j,END) = score;   
      /* Finished updating state END */ 


      }  
    stop_reporting();    
    return TRUE;     
}    


/* Function:  AlignWise_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [AlignWise *]
 *
 */
AlignWise * AlignWise_alloc(void) 
{
    AlignWise * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(AlignWise *) ckalloc (sizeof(AlignWise))) == NULL)  {  
      warn("AlignWise_alloc failed ");   
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


/* Function:  free_AlignWise(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [AlignWise *]
 *
 * Return [UNKN ]  Undocumented return value [AlignWise *]
 *
 */
AlignWise * free_AlignWise(AlignWise * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a AlignWise obj. Should be trappable"); 
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
    /* obj->model is linked in */ 
    /* obj->align is linked in */ 
    /* obj->intronopen is linked in */ 
    /* obj->geneopen is linked in */ 


    ckfree(obj); 
    return NULL; 
}    





#ifdef _cplusplus
}
#endif
