#ifdef _cplusplus
extern "C" {
#endif
#include "statwise10.h"

# line 5 "statwise10.c"


  /*****************   C functions  ****************/
  /*             Written using dynamite            */
  /*            Tue Apr  2 10:59:22 2002           */
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


#define RND_SEQ 0    
#define START 1  
#define END 2    


#define StatWise10_EXPL_MATRIX(this_matrix,i,j,STATE) this_matrix->basematrix->matrix[((j+3)*4)+STATE][i+1]  
#define StatWise10_EXPL_SPECIAL(matrix,i,j,STATE) matrix->basematrix->specmatrix[STATE][j+3] 
#define StatWise10_READ_OFF_ERROR -5
    


#define StatWise10_VSMALL_MATRIX(mat,i,j,STATE) mat->basematrix->matrix[(j+4)%4][((i+1)*4)+STATE]    
#define StatWise10_VSMALL_SPECIAL(mat,i,j,STATE) mat->basematrix->specmatrix[(j+4)%4][STATE] 




/* Function:  search_StatWise10(dbsi,out,exonmodel,seq,codon,intron_open,gene_open)
 *
 * Descrip:    This function makes a database search of StatWise10
 *             It uses the dbsi structure to choose which implementation
 *             to use of the database searching. This way at run time you
 *             can switch between single threaded/multi-threaded or hardware
 *
 *
 * Arg:               dbsi [UNKN ] Undocumented argument [DBSearchImpl *]
 * Arg:                out [UNKN ] Undocumented argument [Hscore *]
 * Arg:          exonmodel [UNKN ] Undocumented argument [SyExonScore*]
 * Arg:                seq [UNKN ] Undocumented argument [ComplexSequence*]
 * Arg:              codon [UNKN ] Undocumented argument [RandomCodonScore*]
 * Arg:        intron_open [UNKN ] Undocumented argument [Score]
 * Arg:          gene_open [UNKN ] Undocumented argument [Score]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type search_StatWise10(DBSearchImpl * dbsi,Hscore * out,SyExonScore* exonmodel,ComplexSequence* seq ,RandomCodonScore* codon,Score intron_open,Score gene_open) 
{
#ifdef PTHREAD   
    int i;   
    int thr_no;  
    pthread_attr_t pat;  
    struct thread_pool_holder_StatWise10 * holder;   
#endif   
    if( out == NULL )    {  
      warn("Passed in a null Hscore object into search_StatWise10. Can't process results!"); 
      return SEARCH_ERROR;   
      }  
    if( dbsi == NULL )   {  
      warn("Passed in a null DBSearchImpl object into search_StatWise10. Can't process results!");   
      return SEARCH_ERROR;   
      }  
    if( dbsi->trace_level > 5 )  
      warn("Asking for trace level of %d in database search for StatWise10, but it was compiled with a trace level of 1702062423. Not all trace statements can be shown",dbsi->trace_level); 
    switch(dbsi->type)   { /*switch on implementation*/ 
      case DBSearchImpl_Serial : 
        return serial_search_StatWise10(out,exonmodel,seq ,codon,intron_open,gene_open); 
      case DBSearchImpl_Pthreads :   
#ifdef PTHREAD   
        holder = (struct thread_pool_holder_StatWise10 *) ckalloc(sizeof(struct thread_pool_holder_StatWise10)); 
        if( holder == NULL )     {  
          warn("Unable to allocated thread pool datastructure...");  
          return SEARCH_ERROR;   
          }  
        holder->out = out;   
        holder->dbsi = dbsi; 
        holder->exonmodel = exonmodel;   
        holder->seq = seq;   
        holder->codon = codon;   
        holder->intron_open = intron_open;   
        holder->gene_open = gene_open;   
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
          if( pthread_create(holder->pool+i,&pat,thread_loop_StatWise10,(void *)holder) )    
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
        warn("You did not specifiy the PTHREAD compile when compiled the C code for StatWise10");    
#endif /* finished threads */    
      default :  
        warn("database search implementation %s was not provided in the compiled dynamite file from StatWise10",impl_string_DBSearchImpl(dbsi)); 
        return SEARCH_ERROR; 
      } /* end of switch on implementation */ 


}    


/* Function:  serial_search_StatWise10(out,exonmodel,seq,codon,intron_open,gene_open)
 *
 * Descrip:    This function makes a database search of StatWise10
 *             It is a single processor implementation
 *
 *
 * Arg:                out [UNKN ] Undocumented argument [Hscore *]
 * Arg:          exonmodel [UNKN ] Undocumented argument [SyExonScore*]
 * Arg:                seq [UNKN ] Undocumented argument [ComplexSequence*]
 * Arg:              codon [UNKN ] Undocumented argument [RandomCodonScore*]
 * Arg:        intron_open [UNKN ] Undocumented argument [Score]
 * Arg:          gene_open [UNKN ] Undocumented argument [Score]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type serial_search_StatWise10(Hscore * out,SyExonScore* exonmodel,ComplexSequence* seq ,RandomCodonScore* codon,Score intron_open,Score gene_open) 
{
    int db_status;   
    int score;   
    int query_pos = 0;   
    int target_pos = 0;  
    DataScore * ds;  


    push_errormsg_stack("Before any actual search in db searching"); 


    target_pos = 0;  




    /* No maximum length - allocated on-the-fly */ 
    score = score_only_StatWise10(exonmodel, seq , codon, intron_open, gene_open);   
    if( should_store_Hscore(out,score) == TRUE )     { /*if storing datascore*/ 
      ds = new_DataScore_from_storage(out);  
      if( ds == NULL )   {  
        warn("StatWise10 search had a memory error in allocating a new_DataScore (?a leak somewhere - DataScore is a very small datastructure"); 
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


/* Function:  score_only_StatWise10(exonmodel,seq,codon,intron_open,gene_open)
 *
 * Descrip:    This function just calculates the score for the matrix
 *             I am pretty sure we can do this better, but hey, for the moment...
 *             It calls /allocate_StatWise10_only
 *
 *
 * Arg:          exonmodel [UNKN ] query data structure [SyExonScore*]
 * Arg:                seq [UNKN ] target data structure [ComplexSequence*]
 * Arg:              codon [UNKN ] Resource [RandomCodonScore*]
 * Arg:        intron_open [UNKN ] Resource [Score]
 * Arg:          gene_open [UNKN ] Resource [Score]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int score_only_StatWise10(SyExonScore* exonmodel,ComplexSequence* seq ,RandomCodonScore* codon,Score intron_open,Score gene_open) 
{
    int bestscore = NEGI;    
    int i;   
    int j;   
    int k;   
    StatWise10 * mat;    


    mat = allocate_StatWise10_only(exonmodel, seq , codon, intron_open, gene_open);  
    if( mat == NULL )    {  
      warn("Memory allocation error in the db search - unable to communicate to calling function. this spells DIASTER!");    
      return NEGI;   
      }  
    if((mat->basematrix = BaseMatrix_alloc_matrix_and_specials(4,(mat->leni + 1) * 4,4,3)) == NULL)  {  
      warn("Score only matrix for StatWise10 cannot be allocated, (asking for 3  by %d  cells)",mat->leni*4);    
      mat = free_StatWise10(mat);    
      return 0;  
      }  
    mat->basematrix->type = BASEMATRIX_TYPE_VERYSMALL;   


    /* Now, initiate matrix */ 
    for(j=0;j<5;j++) {  
      for(i=(-1);i<mat->leni;i++)    {  
        for(k=0;k<4;k++) 
          StatWise10_VSMALL_MATRIX(mat,i,j,k) = NEGI;    
        }  
      StatWise10_VSMALL_SPECIAL(mat,i,j,RND_SEQ) = NEGI; 
      StatWise10_VSMALL_SPECIAL(mat,i,j,START) = 0;  
      StatWise10_VSMALL_SPECIAL(mat,i,j,END) = NEGI; 
      }  


    /* Ok, lets do-o-o-o-o it */ 


    for(j=0;j<mat->lenj;j++) { /*for all target positions*/ 
      auto int score;    
      auto int temp;     
      for(i=0;i<mat->leni;i++)   { /*for all query positions*/ 


        /* For state CODON */ 
        /* setting first movement to score */ 
        score = StatWise10_VSMALL_MATRIX(mat,i-1,j-3,CODON) + mat->codon->codon[CSEQ_GENOMIC_CODON(mat->seq,j)];     
        /* Has restricted position */ 
        if( (i-1) == 0  )    {  
          /* From state RND_SEQ to state CODON */ 
          temp = StatWise10_VSMALL_SPECIAL(mat,i-1,j-3,RND_SEQ) + (mat->codon->codon[CSEQ_GENOMIC_CODON(mat->seq,j)]+mat->gene_open);    
          if( temp  > score )    {  
            score = temp;    
            }  
          }  
        /* From state CODON to state CODON */ 
        temp = StatWise10_VSMALL_MATRIX(mat,i-0,j-3,CODON) + (mat->codon->codon[CSEQ_GENOMIC_CODON(mat->seq,j)]+mat->exonmodel->exon[i]->stay_score);    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTRON_0 to state CODON */ 
        temp = StatWise10_VSMALL_MATRIX(mat,i-1,j-3,INTRON_0) + CSEQ_GENOMIC_3SS(mat->seq,(j-3));    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTRON_1 to state CODON */ 
        temp = StatWise10_VSMALL_MATRIX(mat,i-1,j-2,INTRON_1) + CSEQ_GENOMIC_3SS(mat->seq,(j-2));    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTRON_2 to state CODON */ 
        temp = StatWise10_VSMALL_MATRIX(mat,i-1,j-1,INTRON_2) + CSEQ_GENOMIC_3SS(mat->seq,(j-1));    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for CODON */ 
        /* Add any movement independant score and put away */ 
         StatWise10_VSMALL_MATRIX(mat,i,j,CODON) = score;    


        /* state CODON is a source for special RND_SEQ */ 
        temp = score + (mat->exonmodel->exon[i]->exit_score) + (0) ;     
        if( temp > StatWise10_VSMALL_SPECIAL(mat,i,j,RND_SEQ) )  {  
          StatWise10_VSMALL_SPECIAL(mat,i,j,RND_SEQ) = temp;     
          }  




        /* Finished calculating state CODON */ 


        /* For state INTRON_0 */ 
        /* setting first movement to score */ 
        score = StatWise10_VSMALL_MATRIX(mat,i-0,j-1,CODON) + (CSEQ_GENOMIC_5SS(mat->seq,j)+mat->intron_open);   
        /* From state INTRON_0 to state INTRON_0 */ 
        temp = StatWise10_VSMALL_MATRIX(mat,i-0,j-1,INTRON_0) + 0;   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for INTRON_0 */ 
        /* Add any movement independant score and put away */ 
         StatWise10_VSMALL_MATRIX(mat,i,j,INTRON_0) = score; 


        /* Finished calculating state INTRON_0 */ 


        /* For state INTRON_1 */ 
        /* setting first movement to score */ 
        score = StatWise10_VSMALL_MATRIX(mat,i-0,j-2,CODON) + (CSEQ_GENOMIC_5SS(mat->seq,j)+mat->intron_open);   
        /* From state INTRON_1 to state INTRON_1 */ 
        temp = StatWise10_VSMALL_MATRIX(mat,i-0,j-1,INTRON_1) + 0;   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for INTRON_1 */ 
        /* Add any movement independant score and put away */ 
         StatWise10_VSMALL_MATRIX(mat,i,j,INTRON_1) = score; 


        /* Finished calculating state INTRON_1 */ 


        /* For state INTRON_2 */ 
        /* setting first movement to score */ 
        score = StatWise10_VSMALL_MATRIX(mat,i-0,j-3,CODON) + (CSEQ_GENOMIC_5SS(mat->seq,j)+mat->intron_open);   
        /* From state INTRON_2 to state INTRON_2 */ 
        temp = StatWise10_VSMALL_MATRIX(mat,i-0,j-1,INTRON_2) + 0;   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for INTRON_2 */ 
        /* Add any movement independant score and put away */ 
         StatWise10_VSMALL_MATRIX(mat,i,j,INTRON_2) = score; 


        /* Finished calculating state INTRON_2 */ 
        } /* end of for all query positions */ 




      /* Special state RND_SEQ has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = StatWise10_VSMALL_SPECIAL(mat,0,j,RND_SEQ);    


      /* Source CODON for state RND_SEQ is not special... already calculated */ 
      /* Source RND_SEQ is a special source for RND_SEQ */ 
      temp = StatWise10_VSMALL_SPECIAL(mat,0,j - 1,RND_SEQ) + (0) + (0);     
      if( temp > score ) 
        score = temp;    


      /* Source START is a special source for RND_SEQ */ 
      temp = StatWise10_VSMALL_SPECIAL(mat,0,j - 1,START) + (0) + (0);   
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      StatWise10_VSMALL_SPECIAL(mat,0,j,RND_SEQ) = score;    
      /* Finished updating state RND_SEQ */ 




      /* Special state START has no special to special movements */ 


      /* Special state END has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = StatWise10_VSMALL_SPECIAL(mat,0,j,END);    


      /* Source RND_SEQ is a special source for END */ 
      temp = StatWise10_VSMALL_SPECIAL(mat,0,j - 1,RND_SEQ) + (0) + (0);     
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      StatWise10_VSMALL_SPECIAL(mat,0,j,END) = score;    
      /* Finished updating state END */ 


      if( bestscore < StatWise10_VSMALL_SPECIAL(mat,0,j,END) )   
        bestscore = StatWise10_VSMALL_SPECIAL(mat,0,j,END);  
      } /* end of for all target positions */ 


    mat = free_StatWise10(mat);  
    return bestscore;    
}    


/* Function:  PackAln_bestmemory_StatWise10(exonmodel,seq,codon,intron_open,gene_open,dpenv,dpri)
 *
 * Descrip:    This function chooses the best memory set-up for the alignment
 *             using calls to basematrix, and then implements either a large
 *             or small memory model.
 *
 *             It is the best function to use if you just want an alignment
 *
 *             If you want a label alignment, you will need
 *             /convert_PackAln_to_AlnBlock_StatWise10
 *
 *
 * Arg:          exonmodel [UNKN ] query data structure [SyExonScore*]
 * Arg:                seq [UNKN ] target data structure [ComplexSequence*]
 * Arg:              codon [UNKN ] Resource [RandomCodonScore*]
 * Arg:        intron_open [UNKN ] Resource [Score]
 * Arg:          gene_open [UNKN ] Resource [Score]
 * Arg:              dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:               dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * PackAln_bestmemory_StatWise10(SyExonScore* exonmodel,ComplexSequence* seq ,RandomCodonScore* codon,Score intron_open,Score gene_open,DPEnvelope * dpenv,DPRunImpl * dpri) 
{
    int total;   
    StatWise10 * mat;    
    PackAln * out;   
    DebugMatrix * de;    
    DPRunImplMemory strategy;    
    assert(dpri);    


    total = exonmodel->len * seq->length;    
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


    if( strategy == DPIM_Linear )    {  
      /* use small implementation */ 
      if( (mat=allocate_Small_StatWise10(exonmodel, seq , codon, intron_open, gene_open)) == NULL )  {  
        warn("Unable to allocate small StatWise10 version"); 
        return NULL; 
        }  
      out = PackAln_calculate_Small_StatWise10(mat,dpenv);   
      }  
    else {  
      /* use Large implementation */ 
      if( (mat=allocate_Expl_StatWise10(exonmodel, seq , codon, intron_open, gene_open)) == NULL )   {  
        warn("Unable to allocate large StatWise10 version"); 
        return NULL; 
        }  
      if( dpri->debug == TRUE)   {  
        fatal("Asked for dydebug, but dynamite file not compiled with -g. Need to recompile dynamite source");   
        }  
      else   
        calculate_StatWise10(mat);   
      out =  PackAln_read_Expl_StatWise10(mat);  
      }  


    mat = free_StatWise10(mat);  
    return out;  
}    


/* Function:  allocate_StatWise10_only(exonmodel,seq,codon,intron_open,gene_open)
 *
 * Descrip:    This function only allocates the StatWise10 structure
 *             checks types where possible and determines leni and lenj
 *             The basematrix area is delt with elsewhere
 *
 *
 * Arg:          exonmodel [UNKN ] query data structure [SyExonScore*]
 * Arg:                seq [UNKN ] target data structure [ComplexSequence*]
 * Arg:              codon [UNKN ] Resource [RandomCodonScore*]
 * Arg:        intron_open [UNKN ] Resource [Score]
 * Arg:          gene_open [UNKN ] Resource [Score]
 *
 * Return [UNKN ]  Undocumented return value [StatWise10 *]
 *
 */
StatWise10 * allocate_StatWise10_only(SyExonScore* exonmodel,ComplexSequence* seq ,RandomCodonScore* codon,Score intron_open,Score gene_open) 
{
    StatWise10 * out;    


    if((out= StatWise10_alloc()) == NULL)    {  
      warn("Allocation of basic StatWise10 structure failed...");    
      return NULL;   
      }  


    out->exonmodel = exonmodel;  
    out->seq = seq;  
    out->codon = codon;  
    out->intron_open = intron_open;  
    out->gene_open = gene_open;  
    out->leni = exonmodel->len;  
    out->lenj = seq->length;     
    return out;  
}    


/* Function:  allocate_Expl_StatWise10(exonmodel,seq,codon,intron_open,gene_open)
 *
 * Descrip:    This function allocates the StatWise10 structure
 *             and the basematrix area for explicit memory implementations
 *             It calls /allocate_StatWise10_only
 *
 *
 * Arg:          exonmodel [UNKN ] query data structure [SyExonScore*]
 * Arg:                seq [UNKN ] target data structure [ComplexSequence*]
 * Arg:              codon [UNKN ] Resource [RandomCodonScore*]
 * Arg:        intron_open [UNKN ] Resource [Score]
 * Arg:          gene_open [UNKN ] Resource [Score]
 *
 * Return [UNKN ]  Undocumented return value [StatWise10 *]
 *
 */
StatWise10 * allocate_Expl_StatWise10(SyExonScore* exonmodel,ComplexSequence* seq ,RandomCodonScore* codon,Score intron_open,Score gene_open) 
{
    StatWise10 * out;    


    out = allocate_StatWise10_only(exonmodel, seq , codon, intron_open, gene_open);  
    if( out == NULL )    
      return NULL;   
    if( (out->basematrix = BaseMatrix_alloc_matrix_and_specials((out->lenj+3)*4,(out->leni+1),3,out->lenj+3)) == NULL)   {  
      warn("Explicit matrix StatWise10 cannot be allocated, (asking for %d by %d main cells)",out->leni,out->lenj);  
      free_StatWise10(out);  
      return NULL;   
      }  
    out->basematrix->type = BASEMATRIX_TYPE_EXPLICIT;    
    init_StatWise10(out);    
    return out;  
}    


/* Function:  init_StatWise10(mat)
 *
 * Descrip:    This function initates StatWise10 matrix when in explicit mode
 *             Called in /allocate_Expl_StatWise10
 *
 *
 * Arg:        mat [UNKN ] StatWise10 which contains explicit basematrix memory [StatWise10 *]
 *
 */
void init_StatWise10(StatWise10 * mat) 
{
    register int i;  
    register int j;  
    if( mat->basematrix->type != BASEMATRIX_TYPE_EXPLICIT)   {  
      warn("Cannot iniate matrix, is not an explicit memory type and you have assummed that");   
      return;    
      }  


    for(i= (-1);i<mat->exonmodel->len;i++)   {  
      for(j= (-3);j<4;j++)   {  
        StatWise10_EXPL_MATRIX(mat,i,j,CODON) = NEGI;    
        StatWise10_EXPL_MATRIX(mat,i,j,INTRON_0) = NEGI; 
        StatWise10_EXPL_MATRIX(mat,i,j,INTRON_1) = NEGI; 
        StatWise10_EXPL_MATRIX(mat,i,j,INTRON_2) = NEGI; 
        }  
      }  
    for(j= (-3);j<mat->seq->length;j++)  {  
      for(i= (-1);i<2;i++)   {  
        StatWise10_EXPL_MATRIX(mat,i,j,CODON) = NEGI;    
        StatWise10_EXPL_MATRIX(mat,i,j,INTRON_0) = NEGI; 
        StatWise10_EXPL_MATRIX(mat,i,j,INTRON_1) = NEGI; 
        StatWise10_EXPL_MATRIX(mat,i,j,INTRON_2) = NEGI; 
        }  
      StatWise10_EXPL_SPECIAL(mat,i,j,RND_SEQ) = NEGI;   
      StatWise10_EXPL_SPECIAL(mat,i,j,START) = 0;    
      StatWise10_EXPL_SPECIAL(mat,i,j,END) = NEGI;   
      }  
    return;  
}    


/* Function:  recalculate_PackAln_StatWise10(pal,mat)
 *
 * Descrip:    This function recalculates the PackAln structure produced by StatWise10
 *             For example, in linear space methods this is used to score them
 *
 *
 * Arg:        pal [UNKN ] Undocumented argument [PackAln *]
 * Arg:        mat [UNKN ] Undocumented argument [StatWise10 *]
 *
 */
void recalculate_PackAln_StatWise10(PackAln * pal,StatWise10 * mat) 
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
          if( offi == 1 && offj == 3 && prev->state == CODON )   {  
            pau->score = mat->codon->codon[CSEQ_GENOMIC_CODON(mat->seq,j)] + (0);    
            continue;    
            }  
          if( offj == 3 && prev->state == (RND_SEQ+4) )  {  
            pau->score = (mat->codon->codon[CSEQ_GENOMIC_CODON(mat->seq,j)]+mat->gene_open) + (0);   
            continue;    
            }  
          if( offi == 0 && offj == 3 && prev->state == CODON )   {  
            pau->score = (mat->codon->codon[CSEQ_GENOMIC_CODON(mat->seq,j)]+mat->exonmodel->exon[i]->stay_score) + (0);  
            continue;    
            }  
          if( offi == 1 && offj == 3 && prev->state == INTRON_0 )    {  
            pau->score = CSEQ_GENOMIC_3SS(mat->seq,(j-3)) + (0);     
            continue;    
            }  
          if( offi == 1 && offj == 2 && prev->state == INTRON_1 )    {  
            pau->score = CSEQ_GENOMIC_3SS(mat->seq,(j-2)) + (0);     
            continue;    
            }  
          if( offi == 1 && offj == 1 && prev->state == INTRON_2 )    {  
            pau->score = CSEQ_GENOMIC_3SS(mat->seq,(j-1)) + (0);     
            continue;    
            }  
          warn("In recaluclating PackAln with state CODON, from [%d,%d,%d], got a bad source state. Error!",offi,offj,prev->state);  
          break; 
        case INTRON_0 :  
          if( offi == 0 && offj == 1 && prev->state == CODON )   {  
            pau->score = (CSEQ_GENOMIC_5SS(mat->seq,j)+mat->intron_open) + (0);  
            continue;    
            }  
          if( offi == 0 && offj == 1 && prev->state == INTRON_0 )    {  
            pau->score = 0 + (0);    
            continue;    
            }  
          warn("In recaluclating PackAln with state INTRON_0, from [%d,%d,%d], got a bad source state. Error!",offi,offj,prev->state);   
          break; 
        case INTRON_1 :  
          if( offi == 0 && offj == 2 && prev->state == CODON )   {  
            pau->score = (CSEQ_GENOMIC_5SS(mat->seq,j)+mat->intron_open) + (0);  
            continue;    
            }  
          if( offi == 0 && offj == 1 && prev->state == INTRON_1 )    {  
            pau->score = 0 + (0);    
            continue;    
            }  
          warn("In recaluclating PackAln with state INTRON_1, from [%d,%d,%d], got a bad source state. Error!",offi,offj,prev->state);   
          break; 
        case INTRON_2 :  
          if( offi == 0 && offj == 3 && prev->state == CODON )   {  
            pau->score = (CSEQ_GENOMIC_5SS(mat->seq,j)+mat->intron_open) + (0);  
            continue;    
            }  
          if( offi == 0 && offj == 1 && prev->state == INTRON_2 )    {  
            pau->score = 0 + (0);    
            continue;    
            }  
          warn("In recaluclating PackAln with state INTRON_2, from [%d,%d,%d], got a bad source state. Error!",offi,offj,prev->state);   
          break; 
        case (RND_SEQ+4) :   
          if( offj == 0 && prev->state == CODON )    {  
            /* i here comes from the previous state ;) - not the real one */ 
            i = prev->i; 
            pau->score = mat->exonmodel->exon[i]->exit_score + (0);  
            continue;    
            }  
          if( offj == 1 && prev->state == (RND_SEQ+4) )  {  
            pau->score = 0 + (0);    
            continue;    
            }  
          if( offj == 1 && prev->state == (START+4) )    {  
            pau->score = 0 + (0);    
            continue;    
            }  
          warn("In recaluclating PackAln with state RND_SEQ, got a bad source state. Error!");   
          break; 
        case (START+4) :     
          warn("In recaluclating PackAln with state START, got a bad source state. Error!"); 
          break; 
        case (END+4) :   
          if( offj == 1 && prev->state == (RND_SEQ+4) )  {  
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
#define StatWise10_HIDDEN_MATRIX(thismatrix,i,j,state) (thismatrix->basematrix->matrix[(j-hiddenj+3)][(i+1)*4+state])    
#define StatWise10_DC_SHADOW_MATRIX(thismatrix,i,j,state) (thismatrix->basematrix->matrix[((j+4)*8) % 32][(i+1)*4+state])    
#define StatWise10_HIDDEN_SPECIAL(thismatrix,i,j,state) (thismatrix->basematrix->specmatrix[state][(j+3)])   
#define StatWise10_DC_SHADOW_SPECIAL(thismatrix,i,j,state) (thismatrix->basematrix->specmatrix[state*8][(j+3)])  
#define StatWise10_DC_SHADOW_MATRIX_SP(thismatrix,i,j,state,shadow) (thismatrix->basematrix->matrix[((((j+4)*8)+(shadow+1)) % 32)][(i+1)*4 + state]) 
/* Function:  allocate_Small_StatWise10(exonmodel,seq,codon,intron_open,gene_open)
 *
 * Descrip:    This function allocates the StatWise10 structure
 *             and the basematrix area for a small memory implementations
 *             It calls /allocate_StatWise10_only
 *
 *
 * Arg:          exonmodel [UNKN ] query data structure [SyExonScore*]
 * Arg:                seq [UNKN ] target data structure [ComplexSequence*]
 * Arg:              codon [UNKN ] Resource [RandomCodonScore*]
 * Arg:        intron_open [UNKN ] Resource [Score]
 * Arg:          gene_open [UNKN ] Resource [Score]
 *
 * Return [UNKN ]  Undocumented return value [StatWise10 *]
 *
 */
#define StatWise10_DC_SHADOW_SPECIAL_SP(thismatrix,i,j,state,shadow) (thismatrix->basematrix->specmatrix[state*8 +shadow+1][(j+3)])  
StatWise10 * allocate_Small_StatWise10(SyExonScore* exonmodel,ComplexSequence* seq ,RandomCodonScore* codon,Score intron_open,Score gene_open) 
{
    StatWise10 * out;    


    out = allocate_StatWise10_only(exonmodel, seq , codon, intron_open, gene_open);  
    if( out == NULL )    
      return NULL;   
    out->basematrix = BaseMatrix_alloc_matrix_and_specials(32,(out->leni + 1) * 4,24,out->lenj+3);   
    if(out == NULL)  {  
      warn("Small shadow matrix StatWise10 cannot be allocated, (asking for 4 by %d main cells)",out->leni+2);   
      free_StatWise10(out);  
      return NULL;   
      }  
    out->basematrix->type = BASEMATRIX_TYPE_SHADOW;  
    return out;  
}    


/* Function:  PackAln_calculate_Small_StatWise10(mat,dpenv)
 *
 * Descrip:    This function calculates an alignment for StatWise10 structure in linear space
 *             If you want only the start/end points
 *             use /AlnRangeSet_calculate_Small_StatWise10 
 *
 *             The function basically
 *               finds start/end points 
 *               foreach start/end point 
 *                 calls /full_dc_StatWise10 
 *
 *
 * Arg:          mat [UNKN ] Undocumented argument [StatWise10 *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * PackAln_calculate_Small_StatWise10(StatWise10 * mat,DPEnvelope * dpenv) 
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
      warn("Could not calculate packaln small for StatWise10 due to wrong type of matrix");  
      return NULL;   
      }  


    out = PackAln_alloc_std();   


    start_reporting("Find start end points: ");  
    dc_start_end_calculate_StatWise10(mat,dpenv);    
    score = start_end_find_end_StatWise10(mat,&endj);    
    out->score = score;  
    stopstate = END;
    
    /* Special to specials: have to eat up in strip and then drop back to full_dc for intervening bits */ 
    log_full_error(REPORT,0,"End at %d Score %d",endj,score);    
    stop_reporting();    
    for(;;)  { /*while there are more special bits to recover*/ 
      start_reporting("Special cell aln end   %d:",endj);    
      if( read_special_strip_StatWise10(mat,0,endj,stopstate,&endj,&startstate,out) == FALSE )   {  
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
      starti = StatWise10_DC_SHADOW_SPECIAL_SP(mat,0,endj,temp,0);   
      startj = StatWise10_DC_SHADOW_SPECIAL_SP(mat,0,endj,temp,1);   
      startstate = StatWise10_DC_SHADOW_SPECIAL_SP(mat,0,endj,temp,2);   
      stopi = StatWise10_DC_SHADOW_SPECIAL_SP(mat,0,endj,temp,3);    
      stopj = StatWise10_DC_SHADOW_SPECIAL_SP(mat,0,endj,temp,4);    
      stopstate = StatWise10_DC_SHADOW_SPECIAL_SP(mat,0,endj,temp,5);    


      /* Get out the score of this block. V. important! */ 
      temp = StatWise10_DC_SHADOW_SPECIAL_SP(mat,0,endj,temp,6); 
      totalj = stopj - startj;   
      donej  = 0;    
      start_reporting("Main matrix  aln [%d,%d]:",startj,stopj);     
      if(full_dc_StatWise10(mat,starti,startj,startstate,stopi,stopj,stopstate,out,&donej,totalj,dpenv) == FALSE)    {  
        warn("In the alignment StatWise10 [%d,%d][%d,%d], got a problem. Please report bug ... giving you back a partial alignment",starti,startj,stopi,stopj);  
        return out;  
        }  


      /* now have to figure out which special we came from... yikes */ 
      max_matrix_to_special_StatWise10(mat,starti,startj,startstate,temp,&stopi,&stopj,&stopstate,&temp,NULL);   
      if( stopi == StatWise10_READ_OFF_ERROR)    {  
        warn("In StatWise10 read off ending at %d ... got a bad matrix to special read off... returning partial alignment",startj);  
        invert_PackAln(out); 
        recalculate_PackAln_StatWise10(out,mat); 
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
    recalculate_PackAln_StatWise10(out,mat); 
    return out;  


}    


/* Function:  AlnRangeSet_calculate_Small_StatWise10(mat)
 *
 * Descrip:    This function calculates an alignment for StatWise10 structure in linear space
 *             If you want the full alignment, use /PackAln_calculate_Small_StatWise10 
 *             If you have already got the full alignment, but want the range set, use /AlnRangeSet_from_PackAln_StatWise10
 *             If you have got the small matrix but not the alignment, use /AlnRangeSet_from_StatWise10 
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [StatWise10 *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * AlnRangeSet_calculate_Small_StatWise10(StatWise10 * mat) 
{
    AlnRangeSet * out;   


    start_reporting("Find start end points: ");  
    dc_start_end_calculate_StatWise10(mat,NULL); 
    log_full_error(REPORT,0,"Calculated");   


    out = AlnRangeSet_from_StatWise10(mat);  
    return out;  
}    


/* Function:  AlnRangeSet_from_StatWise10(mat)
 *
 * Descrip:    This function reads off a start/end structure
 *             for StatWise10 structure in linear space
 *             If you want the full alignment use
 *             /PackAln_calculate_Small_StatWise10 
 *             If you have not calculated the matrix use
 *             /AlnRange_calculate_Small_StatWise10
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [StatWise10 *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * AlnRangeSet_from_StatWise10(StatWise10 * mat) 
{
    AlnRangeSet * out;   
    AlnRange * temp; 
    int jpos;    
    int state;   


    if( mat->basematrix->type != BASEMATRIX_TYPE_SHADOW) {  
      warn("Bad error! - non shadow matrix type in AlnRangeSet_from_StatWise10");    
      return NULL;   
      }  


    out = AlnRangeSet_alloc_std();   
    /* Find the end position */ 
    out->score = start_end_find_end_StatWise10(mat,&jpos);   
    state = END; 


    while( (temp = AlnRange_build_StatWise10(mat,jpos,state,&jpos,&state)) != NULL)  
      add_AlnRangeSet(out,temp); 
    return out;  
}    


/* Function:  AlnRange_build_StatWise10(mat,stopj,stopspecstate,startj,startspecstate)
 *
 * Descrip:    This function calculates a single start/end set in linear space
 *             Really a sub-routine for /AlnRangeSet_from_PackAln_StatWise10
 *
 *
 * Arg:                   mat [UNKN ] Undocumented argument [StatWise10 *]
 * Arg:                 stopj [UNKN ] Undocumented argument [int]
 * Arg:         stopspecstate [UNKN ] Undocumented argument [int]
 * Arg:                startj [UNKN ] Undocumented argument [int *]
 * Arg:        startspecstate [UNKN ] Undocumented argument [int *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRange *]
 *
 */
AlnRange * AlnRange_build_StatWise10(StatWise10 * mat,int stopj,int stopspecstate,int * startj,int * startspecstate) 
{
    AlnRange * out;  
    int jpos;    
    int state;   


    if( mat->basematrix->type != BASEMATRIX_TYPE_SHADOW) {  
      warn("Bad error! - non shadow matrix type in AlnRangeSet_from_StatWise10");    
      return NULL;   
      }  


    /* Assumme that we have specials (we should!). Read back along the specials till we have the finish point */ 
    if( read_special_strip_StatWise10(mat,0,stopj,stopspecstate,&jpos,&state,NULL) == FALSE) {  
      warn("In AlnRanger_build_StatWise10 alignment ending at %d, unable to read back specials. Will (evenutally) return a partial range set... BEWARE!",stopj); 
      return NULL;   
      }  
    if( state == START || jpos <= 0) 
      return NULL;   


    out = AlnRange_alloc();  


    out->starti = StatWise10_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,0);   
    out->startj = StatWise10_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,1);   
    out->startstate = StatWise10_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,2);   
    out->stopi = StatWise10_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,3);    
    out->stopj = StatWise10_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,4);    
    out->stopstate = StatWise10_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,5);    
    out->startscore = StatWise10_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,6);   
    out->stopscore = StatWise10_DC_SHADOW_SPECIAL(mat,0,jpos,state); 


    /* Now, we have to figure out where this state came from in the specials */ 
    max_matrix_to_special_StatWise10(mat,out->starti,out->startj,out->startstate,out->startscore,&jpos,startj,startspecstate,&state,NULL);   
    if( jpos == StatWise10_READ_OFF_ERROR)   {  
      warn("In AlnRange_build_StatWise10 alignment ending at %d, with aln range between %d-%d in j, unable to find source special, returning this range, but this could get tricky!",stopj,out->startj,out->stopj);  
      return out;    
      }  


    /* Put in the correct score for startstate, from the special */ 
    out->startscore = StatWise10_DC_SHADOW_SPECIAL(mat,0,*startj,*startspecstate);   
    /* The correct j coords have been put into startj, startspecstate... so just return out */ 
    return out;  
}    


/* Function:  read_hidden_StatWise10(mat,starti,startj,startstate,stopi,stopj,stopstate,out)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [StatWise10 *]
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
boolean read_hidden_StatWise10(StatWise10 * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out) 
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


      max_hidden_StatWise10(mat,startj,i,j,state,isspecial,&i,&j,&state,&isspecial,&cellscore);  


      if( i == StatWise10_READ_OFF_ERROR)    {  
        warn("In StatWise10 hidden read off, between %d:%d,%d:%d - at got bad read off. Problem!",starti,startj,stopi,stopj);    
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
        warn("In StatWise10 hidden read off, between %d:%d,%d:%d - hit start cell, but not in start state. Can't be good!.",starti,startj,stopi,stopj);  
        return FALSE;    
        }  
      }  
    warn("In StatWise10 hidden read off, between %d:%d,%d:%d - gone past start cell (now in %d,%d,%d), can't be good news!.",starti,startj,stopi,stopj,i,j,state);   
    return FALSE;    
}    


/* Function:  max_hidden_StatWise10(mat,hiddenj,i,j,state,isspecial,reti,retj,retstate,retspecial,cellscore)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [StatWise10 *]
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
int max_hidden_StatWise10(StatWise10 * mat,int hiddenj,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore) 
{
    register int temp;   
    register int cscore; 


    *reti = (*retj) = (*retstate) = StatWise10_READ_OFF_ERROR;   


    if( i < 0 || j < 0 || i > mat->exonmodel->len || j > mat->seq->length)   {  
      warn("In StatWise10 matrix special read off - out of bounds on matrix [i,j is %d,%d state %d in standard matrix]",i,j,state);  
      return -1; 
      }  


    /* Then you have to select the correct switch statement to figure out the readoff      */ 
    /* Somewhat odd - reverse the order of calculation and return as soon as it is correct */ 
    cscore = StatWise10_HIDDEN_MATRIX(mat,i,j,state);    
    switch(state)    { /*Switch state */ 
      case CODON :   
        temp = cscore - (CSEQ_GENOMIC_3SS(mat->seq,(j-1))) -  (0);   
        if( temp == StatWise10_HIDDEN_MATRIX(mat,i - 1,j - 1,INTRON_2) ) {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = INTRON_2;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - StatWise10_HIDDEN_MATRIX(mat,i-1,j-1,INTRON_2);    
            }  
          return StatWise10_HIDDEN_MATRIX(mat,i - 1,j - 1,INTRON_2);     
          }  
        temp = cscore - (CSEQ_GENOMIC_3SS(mat->seq,(j-2))) -  (0);   
        if( temp == StatWise10_HIDDEN_MATRIX(mat,i - 1,j - 2,INTRON_1) ) {  
          *reti = i - 1; 
          *retj = j - 2; 
          *retstate = INTRON_1;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - StatWise10_HIDDEN_MATRIX(mat,i-1,j-2,INTRON_1);    
            }  
          return StatWise10_HIDDEN_MATRIX(mat,i - 1,j - 2,INTRON_1);     
          }  
        temp = cscore - (CSEQ_GENOMIC_3SS(mat->seq,(j-3))) -  (0);   
        if( temp == StatWise10_HIDDEN_MATRIX(mat,i - 1,j - 3,INTRON_0) ) {  
          *reti = i - 1; 
          *retj = j - 3; 
          *retstate = INTRON_0;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - StatWise10_HIDDEN_MATRIX(mat,i-1,j-3,INTRON_0);    
            }  
          return StatWise10_HIDDEN_MATRIX(mat,i - 1,j - 3,INTRON_0);     
          }  
        temp = cscore - ((mat->codon->codon[CSEQ_GENOMIC_CODON(mat->seq,j)]+mat->exonmodel->exon[i]->stay_score)) -  (0);    
        if( temp == StatWise10_HIDDEN_MATRIX(mat,i - 0,j - 3,CODON) )    {  
          *reti = i - 0; 
          *retj = j - 3; 
          *retstate = CODON; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - StatWise10_HIDDEN_MATRIX(mat,i-0,j-3,CODON);   
            }  
          return StatWise10_HIDDEN_MATRIX(mat,i - 0,j - 3,CODON);    
          }  
        /* Not allowing special sources.. skipping RND_SEQ */ 
        temp = cscore - (mat->codon->codon[CSEQ_GENOMIC_CODON(mat->seq,j)]) -  (0);  
        if( temp == StatWise10_HIDDEN_MATRIX(mat,i - 1,j - 3,CODON) )    {  
          *reti = i - 1; 
          *retj = j - 3; 
          *retstate = CODON; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - StatWise10_HIDDEN_MATRIX(mat,i-1,j-3,CODON);   
            }  
          return StatWise10_HIDDEN_MATRIX(mat,i - 1,j - 3,CODON);    
          }  
        warn("Major problem (!) - in StatWise10 read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      case INTRON_0 :    
        temp = cscore - (0) -  (0);  
        if( temp == StatWise10_HIDDEN_MATRIX(mat,i - 0,j - 1,INTRON_0) ) {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = INTRON_0;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - StatWise10_HIDDEN_MATRIX(mat,i-0,j-1,INTRON_0);    
            }  
          return StatWise10_HIDDEN_MATRIX(mat,i - 0,j - 1,INTRON_0);     
          }  
        temp = cscore - ((CSEQ_GENOMIC_5SS(mat->seq,j)+mat->intron_open)) -  (0);    
        if( temp == StatWise10_HIDDEN_MATRIX(mat,i - 0,j - 1,CODON) )    {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = CODON; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - StatWise10_HIDDEN_MATRIX(mat,i-0,j-1,CODON);   
            }  
          return StatWise10_HIDDEN_MATRIX(mat,i - 0,j - 1,CODON);    
          }  
        warn("Major problem (!) - in StatWise10 read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      case INTRON_1 :    
        temp = cscore - (0) -  (0);  
        if( temp == StatWise10_HIDDEN_MATRIX(mat,i - 0,j - 1,INTRON_1) ) {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = INTRON_1;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - StatWise10_HIDDEN_MATRIX(mat,i-0,j-1,INTRON_1);    
            }  
          return StatWise10_HIDDEN_MATRIX(mat,i - 0,j - 1,INTRON_1);     
          }  
        temp = cscore - ((CSEQ_GENOMIC_5SS(mat->seq,j)+mat->intron_open)) -  (0);    
        if( temp == StatWise10_HIDDEN_MATRIX(mat,i - 0,j - 2,CODON) )    {  
          *reti = i - 0; 
          *retj = j - 2; 
          *retstate = CODON; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - StatWise10_HIDDEN_MATRIX(mat,i-0,j-2,CODON);   
            }  
          return StatWise10_HIDDEN_MATRIX(mat,i - 0,j - 2,CODON);    
          }  
        warn("Major problem (!) - in StatWise10 read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      case INTRON_2 :    
        temp = cscore - (0) -  (0);  
        if( temp == StatWise10_HIDDEN_MATRIX(mat,i - 0,j - 1,INTRON_2) ) {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = INTRON_2;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - StatWise10_HIDDEN_MATRIX(mat,i-0,j-1,INTRON_2);    
            }  
          return StatWise10_HIDDEN_MATRIX(mat,i - 0,j - 1,INTRON_2);     
          }  
        temp = cscore - ((CSEQ_GENOMIC_5SS(mat->seq,j)+mat->intron_open)) -  (0);    
        if( temp == StatWise10_HIDDEN_MATRIX(mat,i - 0,j - 3,CODON) )    {  
          *reti = i - 0; 
          *retj = j - 3; 
          *retstate = CODON; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - StatWise10_HIDDEN_MATRIX(mat,i-0,j-3,CODON);   
            }  
          return StatWise10_HIDDEN_MATRIX(mat,i - 0,j - 3,CODON);    
          }  
        warn("Major problem (!) - in StatWise10 read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      default:   
        warn("Major problem (!) - in StatWise10 read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      } /* end of Switch state  */ 
}    


/* Function:  read_special_strip_StatWise10(mat,stopi,stopj,stopstate,startj,startstate,out)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [StatWise10 *]
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
boolean read_special_strip_StatWise10(StatWise10 * mat,int stopi,int stopj,int stopstate,int * startj,int * startstate,PackAln * out) 
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
    while( j > StatWise10_DC_SHADOW_SPECIAL_SP(mat,i,j,state,4) && state != START)   { /*while more specials to eat up*/ 
      /* Put away current state, if we should */ 
      if(out != NULL)    {  
        pau = PackAlnUnit_alloc();  /* Should deal with memory overflow */ 
        pau->i = i;  
        pau->j = j;  
        pau->state =  state + 4; 
        add_PackAln(out,pau);    
        }  


      max_special_strip_StatWise10(mat,i,j,state,isspecial,&i,&j,&state,&isspecial,&cellscore);  
      if( i == StatWise10_READ_OFF_ERROR)    {  
        warn("In special strip read StatWise10, got a bad read off error. Sorry!");  
        return FALSE;    
        }  
      } /* end of while more specials to eat up */ 


    /* check to see we have not gone too far! */ 
    if( state != START && j < StatWise10_DC_SHADOW_SPECIAL_SP(mat,i,j,state,4))  {  
      warn("In special strip read StatWise10, at special [%d] state [%d] overshot!",j,state);    
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


/* Function:  max_special_strip_StatWise10(mat,i,j,state,isspecial,reti,retj,retstate,retspecial,cellscore)
 *
 * Descrip:    A pretty intense internal function. Deals with read-off only in specials
 *
 *
 * Arg:               mat [UNKN ] Undocumented argument [StatWise10 *]
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
int max_special_strip_StatWise10(StatWise10 * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore) 
{
    int temp;    
    int cscore;  


    *reti = (*retj) = (*retstate) = StatWise10_READ_OFF_ERROR;   
    if( isspecial == FALSE ) {  
      warn("In special strip max function for StatWise10, got a non special start point. Problem! (bad!)");  
      return (-1);   
      }  


    if( j < 0 || j > mat->seq->length)   {  
      warn("In StatWise10 matrix special read off - out of bounds on matrix [j is %d in special]",j);    
      return -1; 
      }  


    cscore = StatWise10_DC_SHADOW_SPECIAL(mat,i,j,state);    
    switch(state)    { /*switch on special states*/ 
      case RND_SEQ :     
        /* source START is a special */ 
        temp = cscore - (0) - (0);   
        if( temp == StatWise10_DC_SHADOW_SPECIAL(mat,i - 0,j - 1,START) )    {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = START; 
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - StatWise10_DC_SHADOW_SPECIAL(mat,i-0,j-1,START);   
            }  
          return StatWise10_DC_SHADOW_MATRIX(mat,i - 0,j - 1,START) ;    
          }  
        /* source RND_SEQ is a special */ 
        temp = cscore - (0) - (0);   
        if( temp == StatWise10_DC_SHADOW_SPECIAL(mat,i - 0,j - 1,RND_SEQ) )  {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = RND_SEQ;   
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - StatWise10_DC_SHADOW_SPECIAL(mat,i-0,j-1,RND_SEQ);     
            }  
          return StatWise10_DC_SHADOW_MATRIX(mat,i - 0,j - 1,RND_SEQ) ;  
          }  
        /* Source CODON is not a special */ 
      case START :   
      case END :     
        /* source RND_SEQ is a special */ 
        temp = cscore - (0) - (0);   
        if( temp == StatWise10_DC_SHADOW_SPECIAL(mat,i - 0,j - 1,RND_SEQ) )  {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = RND_SEQ;   
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - StatWise10_DC_SHADOW_SPECIAL(mat,i-0,j-1,RND_SEQ);     
            }  
          return StatWise10_DC_SHADOW_MATRIX(mat,i - 0,j - 1,RND_SEQ) ;  
          }  
      default:   
        warn("Major problem (!) - in StatWise10 special strip read off, position %d,%d state %d no source found  dropped into default on source switch!",i,j,state); 
        return (-1); 
      } /* end of switch on special states */ 
}    


/* Function:  max_matrix_to_special_StatWise10(mat,i,j,state,cscore,reti,retj,retstate,retspecial,cellscore)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [StatWise10 *]
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
int max_matrix_to_special_StatWise10(StatWise10 * mat,int i,int j,int state,int cscore,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore) 
{
    int temp;    
    *reti = (*retj) = (*retstate) = StatWise10_READ_OFF_ERROR;   


    if( j < 0 || j > mat->lenj)  {  
      warn("In StatWise10 matrix to special read off - out of bounds on matrix [j is %d in special]",j); 
      return -1; 
      }  


    switch(state)    { /*Switch state */ 
      case CODON :   
        /* Source INTRON_2 is not a special, should not get here! */ 
        /* Source INTRON_1 is not a special, should not get here! */ 
        /* Source INTRON_0 is not a special, should not get here! */ 
        /* Source CODON is not a special, should not get here! */ 
        temp = cscore - ((mat->codon->codon[CSEQ_GENOMIC_CODON(mat->seq,j)]+mat->gene_open)) -  (0);     
        if( temp == StatWise10_DC_SHADOW_SPECIAL(mat,i - 1,j - 3,RND_SEQ) )  {  
          *reti = i - 1; 
          *retj = j - 3; 
          *retstate = RND_SEQ;   
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - StatWise10_DC_SHADOW_SPECIAL(mat,i-1,j-3,RND_SEQ);     
            }  
          return StatWise10_DC_SHADOW_MATRIX(mat,i - 1,j - 3,RND_SEQ) ;  
          }  
        /* Source CODON is not a special, should not get here! */ 
        warn("Major problem (!) - in StatWise10 matrix to special read off, position %d,%d state %d no source found!",i,j,state);    
        return (-1); 
      case INTRON_0 :    
        /* Source INTRON_0 is not a special, should not get here! */ 
        /* Source CODON is not a special, should not get here! */ 
        warn("Major problem (!) - in StatWise10 matrix to special read off, position %d,%d state %d no source found!",i,j,state);    
        return (-1); 
      case INTRON_1 :    
        /* Source INTRON_1 is not a special, should not get here! */ 
        /* Source CODON is not a special, should not get here! */ 
        warn("Major problem (!) - in StatWise10 matrix to special read off, position %d,%d state %d no source found!",i,j,state);    
        return (-1); 
      case INTRON_2 :    
        /* Source INTRON_2 is not a special, should not get here! */ 
        /* Source CODON is not a special, should not get here! */ 
        warn("Major problem (!) - in StatWise10 matrix to special read off, position %d,%d state %d no source found!",i,j,state);    
        return (-1); 
      default:   
        warn("Major problem (!) - in StatWise10 read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      } /* end of Switch state  */ 


}    


/* Function:  calculate_hidden_StatWise10(mat,starti,startj,startstate,stopi,stopj,stopstate,dpenv)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [StatWise10 *]
 * Arg:            starti [UNKN ] Undocumented argument [int]
 * Arg:            startj [UNKN ] Undocumented argument [int]
 * Arg:        startstate [UNKN ] Undocumented argument [int]
 * Arg:             stopi [UNKN ] Undocumented argument [int]
 * Arg:             stopj [UNKN ] Undocumented argument [int]
 * Arg:         stopstate [UNKN ] Undocumented argument [int]
 * Arg:             dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 */
void calculate_hidden_StatWise10(StatWise10 * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv) 
{
    register int i;  
    register int j;  
    register int score;  
    register int temp;   
    register int hiddenj;    


    hiddenj = startj;    


    init_hidden_StatWise10(mat,starti,startj,stopi,stopj);   


    StatWise10_HIDDEN_MATRIX(mat,starti,startj,startstate) = 0;  


    for(j=startj;j<=stopj;j++)   {  
      for(i=starti;i<=stopi;i++) {  
        /* Should *not* do very first cell as this is the one set to zero in one state! */ 
        if( i == starti && j == startj ) 
          continue;  
        if( dpenv != NULL && is_in_DPEnvelope(dpenv,i,j) == FALSE )  { /*Is not in envelope*/ 
          StatWise10_HIDDEN_MATRIX(mat,i,j,CODON) = NEGI;    
          StatWise10_HIDDEN_MATRIX(mat,i,j,INTRON_0) = NEGI;     
          StatWise10_HIDDEN_MATRIX(mat,i,j,INTRON_1) = NEGI;     
          StatWise10_HIDDEN_MATRIX(mat,i,j,INTRON_2) = NEGI;     
          continue;  
          } /* end of Is not in envelope */ 


        /* For state CODON */ 
        /* setting first movement to score */ 
        score = StatWise10_HIDDEN_MATRIX(mat,i-1,j-3,CODON) + mat->codon->codon[CSEQ_GENOMIC_CODON(mat->seq,j)];     
        /* From state CODON to state CODON */ 
        temp = StatWise10_HIDDEN_MATRIX(mat,i-0,j-3,CODON) + (mat->codon->codon[CSEQ_GENOMIC_CODON(mat->seq,j)]+mat->exonmodel->exon[i]->stay_score);    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTRON_0 to state CODON */ 
        temp = StatWise10_HIDDEN_MATRIX(mat,i-1,j-3,INTRON_0) + CSEQ_GENOMIC_3SS(mat->seq,(j-3));    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTRON_1 to state CODON */ 
        temp = StatWise10_HIDDEN_MATRIX(mat,i-1,j-2,INTRON_1) + CSEQ_GENOMIC_3SS(mat->seq,(j-2));    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTRON_2 to state CODON */ 
        temp = StatWise10_HIDDEN_MATRIX(mat,i-1,j-1,INTRON_2) + CSEQ_GENOMIC_3SS(mat->seq,(j-1));    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for CODON */ 
        /* Add any movement independant score and put away */ 
         StatWise10_HIDDEN_MATRIX(mat,i,j,CODON) = score;    
        /* Finished calculating state CODON */ 


        /* For state INTRON_0 */ 
        /* setting first movement to score */ 
        score = StatWise10_HIDDEN_MATRIX(mat,i-0,j-1,CODON) + (CSEQ_GENOMIC_5SS(mat->seq,j)+mat->intron_open);   
        /* From state INTRON_0 to state INTRON_0 */ 
        temp = StatWise10_HIDDEN_MATRIX(mat,i-0,j-1,INTRON_0) + 0;   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for INTRON_0 */ 
        /* Add any movement independant score and put away */ 
         StatWise10_HIDDEN_MATRIX(mat,i,j,INTRON_0) = score; 
        /* Finished calculating state INTRON_0 */ 


        /* For state INTRON_1 */ 
        /* setting first movement to score */ 
        score = StatWise10_HIDDEN_MATRIX(mat,i-0,j-2,CODON) + (CSEQ_GENOMIC_5SS(mat->seq,j)+mat->intron_open);   
        /* From state INTRON_1 to state INTRON_1 */ 
        temp = StatWise10_HIDDEN_MATRIX(mat,i-0,j-1,INTRON_1) + 0;   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for INTRON_1 */ 
        /* Add any movement independant score and put away */ 
         StatWise10_HIDDEN_MATRIX(mat,i,j,INTRON_1) = score; 
        /* Finished calculating state INTRON_1 */ 


        /* For state INTRON_2 */ 
        /* setting first movement to score */ 
        score = StatWise10_HIDDEN_MATRIX(mat,i-0,j-3,CODON) + (CSEQ_GENOMIC_5SS(mat->seq,j)+mat->intron_open);   
        /* From state INTRON_2 to state INTRON_2 */ 
        temp = StatWise10_HIDDEN_MATRIX(mat,i-0,j-1,INTRON_2) + 0;   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for INTRON_2 */ 
        /* Add any movement independant score and put away */ 
         StatWise10_HIDDEN_MATRIX(mat,i,j,INTRON_2) = score; 
        /* Finished calculating state INTRON_2 */ 
        }  
      }  


    return;  
}    


/* Function:  init_hidden_StatWise10(mat,starti,startj,stopi,stopj)
 *
 * Descrip: No Description
 *
 * Arg:           mat [UNKN ] Undocumented argument [StatWise10 *]
 * Arg:        starti [UNKN ] Undocumented argument [int]
 * Arg:        startj [UNKN ] Undocumented argument [int]
 * Arg:         stopi [UNKN ] Undocumented argument [int]
 * Arg:         stopj [UNKN ] Undocumented argument [int]
 *
 */
void init_hidden_StatWise10(StatWise10 * mat,int starti,int startj,int stopi,int stopj) 
{
    register int i;  
    register int j;  
    register int hiddenj;    


    hiddenj = startj;    
    for(j=(startj-3);j<=stopj;j++)   {  
      for(i=(starti-1);i<=stopi;i++) {  
        StatWise10_HIDDEN_MATRIX(mat,i,j,CODON) = NEGI;
 
        StatWise10_HIDDEN_MATRIX(mat,i,j,INTRON_0) = NEGI;
  
        StatWise10_HIDDEN_MATRIX(mat,i,j,INTRON_1) = NEGI;
  
        StatWise10_HIDDEN_MATRIX(mat,i,j,INTRON_2) = NEGI;
  
        }  
      }  


    return;  
}    


/* Function:  full_dc_StatWise10(mat,starti,startj,startstate,stopi,stopj,stopstate,out,donej,totalj,dpenv)
 *
 * Descrip:    The main divide-and-conquor routine. Basically, call /PackAln_calculate_small_StatWise10
 *             Not this function, which is pretty hard core. 
 *             Function is given start/end points (in main matrix) for alignment
 *             It does some checks, decides whether start/end in j is small enough for explicit calc
 *               - if yes, calculates it, reads off into PackAln (out), adds the j distance to donej and returns TRUE
 *               - if no,  uses /do_dc_single_pass_StatWise10 to get mid-point
 *                          saves midpoint, and calls itself to do right portion then left portion
 *             right then left ensures PackAln is added the 'right' way, ie, back-to-front
 *             returns FALSE on any error, with a warning
 *
 *
 * Arg:               mat [UNKN ] Matrix with small memory implementation [StatWise10 *]
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
boolean full_dc_StatWise10(StatWise10 * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out,int * donej,int totalj,DPEnvelope * dpenv) 
{
    int lstarti; 
    int lstartj; 
    int lstate;  


    if( mat->basematrix->type != BASEMATRIX_TYPE_SHADOW) {  
      warn("*Very* bad error! - non shadow matrix type in full_dc_StatWise10");  
      return FALSE;  
      }  


    if( starti == -1 || startj == -1 || startstate == -1 || stopi == -1 || stopstate == -1)  {  
      warn("In full dc program, passed bad indices, indices passed were %d:%d[%d] to %d:%d[%d]\n",starti,startj,startstate,stopi,stopj,stopstate);   
      return FALSE;  
      }  


    if( stopj - startj < 15) {  
      log_full_error(REPORT,0,"[%d,%d][%d,%d] Explicit read off",starti,startj,stopi,stopj);/* Build hidden explicit matrix */ 
      calculate_hidden_StatWise10(mat,starti,startj,startstate,stopi,stopj,stopstate,dpenv);     
      *donej += (stopj - startj);   /* Now read it off into out */ 
      if( read_hidden_StatWise10(mat,starti,startj,startstate,stopi,stopj,stopstate,out) == FALSE)   {  
        warn("In full dc, at %d:%d,%d:%d got a bad hidden explicit read off... ",starti,startj,stopi,stopj); 
        return FALSE;    
        }  
      return TRUE;   
      }  


/* In actual divide and conquor */ 
    if( do_dc_single_pass_StatWise10(mat,starti,startj,startstate,stopi,stopj,stopstate,dpenv,(int)(*donej*100)/totalj) == FALSE)    {  
      warn("In divide and conquor for StatWise10, at bound %d:%d to %d:%d, unable to calculate midpoint. Problem!",starti,startj,stopi,stopj);   
      return FALSE;  
      }  


/* Ok... now we have to call on each side of the matrix */ 
/* We have to retrieve left hand side positions, as they will be vapped by the time we call LHS */ 
    lstarti= StatWise10_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,0);    
    lstartj= StatWise10_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,1);    
    lstate = StatWise10_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,2);    


/* Call on right hand side: this lets us do the correct read off */ 
    if( full_dc_StatWise10(mat,StatWise10_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,3),StatWise10_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,4),StatWise10_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,5),stopi,stopj,stopstate,out,donej,totalj,dpenv) == FALSE)   {  
/* Warning already issued, simply chained back up to top */ 
      return FALSE;  
      }  
/* Call on left hand side */ 
    if( full_dc_StatWise10(mat,starti,startj,startstate,lstarti,lstartj,lstate,out,donej,totalj,dpenv) == FALSE) {  
/* Warning already issued, simply chained back up to top */ 
      return FALSE;  
      }  


    return TRUE;     
}    


/* Function:  do_dc_single_pass_StatWise10(mat,starti,startj,startstate,stopi,stopj,stopstate,dpenv,perc_done)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [StatWise10 *]
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
boolean do_dc_single_pass_StatWise10(StatWise10 * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv,int perc_done) 
{
    int halfj;   
    halfj = startj + ((stopj - startj)/2);   


    init_dc_StatWise10(mat); 


    StatWise10_DC_SHADOW_MATRIX(mat,starti,startj,startstate) = 0;   
    run_up_dc_StatWise10(mat,starti,stopi,startj,halfj-1,dpenv,perc_done);   
    push_dc_at_merge_StatWise10(mat,starti,stopi,halfj,&halfj,dpenv);    
    follow_on_dc_StatWise10(mat,starti,stopi,halfj,stopj,dpenv,perc_done);   
    return TRUE; 
}    


/* Function:  push_dc_at_merge_StatWise10(mat,starti,stopi,startj,stopj,dpenv)
 *
 * Descrip: No Description
 *
 * Arg:           mat [UNKN ] Undocumented argument [StatWise10 *]
 * Arg:        starti [UNKN ] Undocumented argument [int]
 * Arg:         stopi [UNKN ] Undocumented argument [int]
 * Arg:        startj [UNKN ] Undocumented argument [int]
 * Arg:         stopj [UNKN ] Undocumented argument [int *]
 * Arg:         dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 */
void push_dc_at_merge_StatWise10(StatWise10 * mat,int starti,int stopi,int startj,int * stopj,DPEnvelope * dpenv) 
{
    register int i;  
    register int j;  
    register int k;  
    register int count;  
    register int mergej;/* Sources below this j will be stamped by triples */ 
    register int score;  
    register int temp;   


    mergej = startj -1;  
    for(count=0,j=startj;count<3;count++,j++)    {  
      for(i=starti;i<=stopi;i++) {  
        if( dpenv != NULL && is_in_DPEnvelope(dpenv,i,j) == FALSE )  { /*Is not in envelope*/ 
          StatWise10_DC_SHADOW_MATRIX(mat,i,j,CODON) = NEGI;     
          StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,CODON,0) = (-100);  
          StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,CODON,1) = (-100);  
          StatWise10_DC_SHADOW_MATRIX(mat,i,j,INTRON_0) = NEGI;  
          StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_0,0) = (-100);   
          StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_0,1) = (-100);   
          StatWise10_DC_SHADOW_MATRIX(mat,i,j,INTRON_1) = NEGI;  
          StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_1,0) = (-100);   
          StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_1,1) = (-100);   
          StatWise10_DC_SHADOW_MATRIX(mat,i,j,INTRON_2) = NEGI;  
          StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_2,0) = (-100);   
          StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_2,1) = (-100);   
          continue;  
          } /* end of Is not in envelope */ 


        /* For state CODON, pushing when j - offj <= mergej */ 
        score = StatWise10_DC_SHADOW_MATRIX(mat,i-1,j-3,CODON) + mat->codon->codon[CSEQ_GENOMIC_CODON(mat->seq,j)];  
        if( j - 3 <= mergej) {  
          StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,CODON,0) = i-1; 
          StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,CODON,1) = j-3; 
          StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,CODON,2) = CODON;   
          StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,CODON,3) = i;   
          StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,CODON,4) = j;   
          StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,CODON,5) = CODON;   
          }  
        else {  
          for(k=0;k<7;k++)   
            StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,CODON,k) = StatWise10_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 3,CODON,k);   
          }  


        temp = StatWise10_DC_SHADOW_MATRIX(mat,i-0,j-3,CODON) + (mat->codon->codon[CSEQ_GENOMIC_CODON(mat->seq,j)]+mat->exonmodel->exon[i]->stay_score);     
        if( temp > score)    {  
          score = temp;  


          if( j - 3 <= mergej)   {  
            StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,CODON,0) = i-0;   
            StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,CODON,1) = j-3;   
            StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,CODON,2) = CODON; 
            StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,CODON,3) = i; 
            StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,CODON,4) = j; 
            StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,CODON,5) = CODON; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,CODON,k) = StatWise10_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 3,CODON,k); 
            }  
          }  


        temp = StatWise10_DC_SHADOW_MATRIX(mat,i-1,j-3,INTRON_0) + CSEQ_GENOMIC_3SS(mat->seq,(j-3));     
        if( temp > score)    {  
          score = temp;  


          if( j - 3 <= mergej)   {  
            StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,CODON,0) = i-1;   
            StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,CODON,1) = j-3;   
            StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,CODON,2) = INTRON_0;  
            StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,CODON,3) = i; 
            StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,CODON,4) = j; 
            StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,CODON,5) = CODON; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,CODON,k) = StatWise10_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 3,INTRON_0,k);  
            }  
          }  


        temp = StatWise10_DC_SHADOW_MATRIX(mat,i-1,j-2,INTRON_1) + CSEQ_GENOMIC_3SS(mat->seq,(j-2));     
        if( temp > score)    {  
          score = temp;  


          if( j - 2 <= mergej)   {  
            StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,CODON,0) = i-1;   
            StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,CODON,1) = j-2;   
            StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,CODON,2) = INTRON_1;  
            StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,CODON,3) = i; 
            StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,CODON,4) = j; 
            StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,CODON,5) = CODON; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,CODON,k) = StatWise10_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 2,INTRON_1,k);  
            }  
          }  


        temp = StatWise10_DC_SHADOW_MATRIX(mat,i-1,j-1,INTRON_2) + CSEQ_GENOMIC_3SS(mat->seq,(j-1));     
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,CODON,0) = i-1;   
            StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,CODON,1) = j-1;   
            StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,CODON,2) = INTRON_2;  
            StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,CODON,3) = i; 
            StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,CODON,4) = j; 
            StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,CODON,5) = CODON; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,CODON,k) = StatWise10_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,INTRON_2,k);  
            }  
          }  
        /* Add any movement independant score */ 
        StatWise10_DC_SHADOW_MATRIX(mat,i,j,CODON) = score;  
        /* Finished with state CODON */ 


        /* For state INTRON_0, pushing when j - offj <= mergej */ 
        score = StatWise10_DC_SHADOW_MATRIX(mat,i-0,j-1,CODON) + (CSEQ_GENOMIC_5SS(mat->seq,j)+mat->intron_open);    
        if( j - 1 <= mergej) {  
          StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_0,0) = i-0;  
          StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_0,1) = j-1;  
          StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_0,2) = CODON;    
          StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_0,3) = i;    
          StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_0,4) = j;    
          StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_0,5) = INTRON_0; 
          }  
        else {  
          for(k=0;k<7;k++)   
            StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_0,k) = StatWise10_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,CODON,k);    
          }  


        temp = StatWise10_DC_SHADOW_MATRIX(mat,i-0,j-1,INTRON_0) + 0;    
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_0,0) = i-0;    
            StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_0,1) = j-1;    
            StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_0,2) = INTRON_0;   
            StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_0,3) = i;  
            StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_0,4) = j;  
            StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_0,5) = INTRON_0;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_0,k) = StatWise10_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,INTRON_0,k);   
            }  
          }  
        /* Add any movement independant score */ 
        StatWise10_DC_SHADOW_MATRIX(mat,i,j,INTRON_0) = score;   
        /* Finished with state INTRON_0 */ 


        /* For state INTRON_1, pushing when j - offj <= mergej */ 
        score = StatWise10_DC_SHADOW_MATRIX(mat,i-0,j-2,CODON) + (CSEQ_GENOMIC_5SS(mat->seq,j)+mat->intron_open);    
        if( j - 2 <= mergej) {  
          StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_1,0) = i-0;  
          StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_1,1) = j-2;  
          StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_1,2) = CODON;    
          StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_1,3) = i;    
          StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_1,4) = j;    
          StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_1,5) = INTRON_1; 
          }  
        else {  
          for(k=0;k<7;k++)   
            StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_1,k) = StatWise10_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 2,CODON,k);    
          }  


        temp = StatWise10_DC_SHADOW_MATRIX(mat,i-0,j-1,INTRON_1) + 0;    
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_1,0) = i-0;    
            StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_1,1) = j-1;    
            StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_1,2) = INTRON_1;   
            StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_1,3) = i;  
            StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_1,4) = j;  
            StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_1,5) = INTRON_1;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_1,k) = StatWise10_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,INTRON_1,k);   
            }  
          }  
        /* Add any movement independant score */ 
        StatWise10_DC_SHADOW_MATRIX(mat,i,j,INTRON_1) = score;   
        /* Finished with state INTRON_1 */ 


        /* For state INTRON_2, pushing when j - offj <= mergej */ 
        score = StatWise10_DC_SHADOW_MATRIX(mat,i-0,j-3,CODON) + (CSEQ_GENOMIC_5SS(mat->seq,j)+mat->intron_open);    
        if( j - 3 <= mergej) {  
          StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_2,0) = i-0;  
          StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_2,1) = j-3;  
          StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_2,2) = CODON;    
          StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_2,3) = i;    
          StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_2,4) = j;    
          StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_2,5) = INTRON_2; 
          }  
        else {  
          for(k=0;k<7;k++)   
            StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_2,k) = StatWise10_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 3,CODON,k);    
          }  


        temp = StatWise10_DC_SHADOW_MATRIX(mat,i-0,j-1,INTRON_2) + 0;    
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_2,0) = i-0;    
            StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_2,1) = j-1;    
            StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_2,2) = INTRON_2;   
            StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_2,3) = i;  
            StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_2,4) = j;  
            StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_2,5) = INTRON_2;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_2,k) = StatWise10_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,INTRON_2,k);   
            }  
          }  
        /* Add any movement independant score */ 
        StatWise10_DC_SHADOW_MATRIX(mat,i,j,INTRON_2) = score;   
        /* Finished with state INTRON_2 */ 
        }  
      }  
    /* Put back j into * stop j so that calling function gets it correct */ 
    if( stopj == NULL)   
      warn("Bad news... NULL stopj pointer in push dc function. This means that calling function does not know how many cells I have done!");    
    else 
      *stopj = j;    


    return;  
}    


/* Function:  follow_on_dc_StatWise10(mat,starti,stopi,startj,stopj,dpenv,perc_done)
 *
 * Descrip: No Description
 *
 * Arg:              mat [UNKN ] Undocumented argument [StatWise10 *]
 * Arg:           starti [UNKN ] Undocumented argument [int]
 * Arg:            stopi [UNKN ] Undocumented argument [int]
 * Arg:           startj [UNKN ] Undocumented argument [int]
 * Arg:            stopj [UNKN ] Undocumented argument [int]
 * Arg:            dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:        perc_done [UNKN ] Undocumented argument [int]
 *
 */
void follow_on_dc_StatWise10(StatWise10 * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done) 
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
          StatWise10_DC_SHADOW_MATRIX(mat,i,j,CODON) = NEGI;     
          StatWise10_DC_SHADOW_MATRIX(mat,i,j,INTRON_0) = NEGI;  
          StatWise10_DC_SHADOW_MATRIX(mat,i,j,INTRON_1) = NEGI;  
          StatWise10_DC_SHADOW_MATRIX(mat,i,j,INTRON_2) = NEGI;  
          continue;  
          } /* end of Is not in envelope */ 
        if( num % 1000 == 0 )    
          log_full_error(REPORT,0,"[%d%%%% done]After  mid-j %5d Cells done %d%%%%",perc_done,startj,(num*100)/total);   


        /* For state CODON */ 
        /* setting first movement to score */ 
        score = StatWise10_DC_SHADOW_MATRIX(mat,i-1,j-3,CODON) + mat->codon->codon[CSEQ_GENOMIC_CODON(mat->seq,j)];  
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = StatWise10_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 3,CODON,k);  
        /* From state CODON to state CODON */ 
        temp = StatWise10_DC_SHADOW_MATRIX(mat,i-0,j-3,CODON) + (mat->codon->codon[CSEQ_GENOMIC_CODON(mat->seq,j)]+mat->exonmodel->exon[i]->stay_score);     
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = StatWise10_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 3,CODON,k);    
          }  
        /* From state INTRON_0 to state CODON */ 
        temp = StatWise10_DC_SHADOW_MATRIX(mat,i-1,j-3,INTRON_0) + CSEQ_GENOMIC_3SS(mat->seq,(j-3));     
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = StatWise10_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 3,INTRON_0,k); 
          }  
        /* From state INTRON_1 to state CODON */ 
        temp = StatWise10_DC_SHADOW_MATRIX(mat,i-1,j-2,INTRON_1) + CSEQ_GENOMIC_3SS(mat->seq,(j-2));     
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = StatWise10_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 2,INTRON_1,k); 
          }  
        /* From state INTRON_2 to state CODON */ 
        temp = StatWise10_DC_SHADOW_MATRIX(mat,i-1,j-1,INTRON_2) + CSEQ_GENOMIC_3SS(mat->seq,(j-1));     
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = StatWise10_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,INTRON_2,k); 
          }  


        /* Ok - finished max calculation for CODON */ 
        /* Add any movement independant score and put away */ 
         StatWise10_DC_SHADOW_MATRIX(mat,i,j,CODON) = score; 
        for(k=0;k<7;k++) 
          StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,CODON,k) = localshadow[k];  
        /* Now figure out if any specials need this score */ 
        /* Finished calculating state CODON */ 


        /* For state INTRON_0 */ 
        /* setting first movement to score */ 
        score = StatWise10_DC_SHADOW_MATRIX(mat,i-0,j-1,CODON) + (CSEQ_GENOMIC_5SS(mat->seq,j)+mat->intron_open);    
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = StatWise10_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,CODON,k);  
        /* From state INTRON_0 to state INTRON_0 */ 
        temp = StatWise10_DC_SHADOW_MATRIX(mat,i-0,j-1,INTRON_0) + 0;    
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = StatWise10_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,INTRON_0,k); 
          }  


        /* Ok - finished max calculation for INTRON_0 */ 
        /* Add any movement independant score and put away */ 
         StatWise10_DC_SHADOW_MATRIX(mat,i,j,INTRON_0) = score;  
        for(k=0;k<7;k++) 
          StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_0,k) = localshadow[k];   
        /* Now figure out if any specials need this score */ 
        /* Finished calculating state INTRON_0 */ 


        /* For state INTRON_1 */ 
        /* setting first movement to score */ 
        score = StatWise10_DC_SHADOW_MATRIX(mat,i-0,j-2,CODON) + (CSEQ_GENOMIC_5SS(mat->seq,j)+mat->intron_open);    
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = StatWise10_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 2,CODON,k);  
        /* From state INTRON_1 to state INTRON_1 */ 
        temp = StatWise10_DC_SHADOW_MATRIX(mat,i-0,j-1,INTRON_1) + 0;    
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = StatWise10_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,INTRON_1,k); 
          }  


        /* Ok - finished max calculation for INTRON_1 */ 
        /* Add any movement independant score and put away */ 
         StatWise10_DC_SHADOW_MATRIX(mat,i,j,INTRON_1) = score;  
        for(k=0;k<7;k++) 
          StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_1,k) = localshadow[k];   
        /* Now figure out if any specials need this score */ 
        /* Finished calculating state INTRON_1 */ 


        /* For state INTRON_2 */ 
        /* setting first movement to score */ 
        score = StatWise10_DC_SHADOW_MATRIX(mat,i-0,j-3,CODON) + (CSEQ_GENOMIC_5SS(mat->seq,j)+mat->intron_open);    
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = StatWise10_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 3,CODON,k);  
        /* From state INTRON_2 to state INTRON_2 */ 
        temp = StatWise10_DC_SHADOW_MATRIX(mat,i-0,j-1,INTRON_2) + 0;    
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = StatWise10_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,INTRON_2,k); 
          }  


        /* Ok - finished max calculation for INTRON_2 */ 
        /* Add any movement independant score and put away */ 
         StatWise10_DC_SHADOW_MATRIX(mat,i,j,INTRON_2) = score;  
        for(k=0;k<7;k++) 
          StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_2,k) = localshadow[k];   
        /* Now figure out if any specials need this score */ 
        /* Finished calculating state INTRON_2 */ 
        } /* end of this is strip */ 
      } /* end of for each valid j column */ 


/* Function:  run_up_dc_StatWise10(mat,starti,stopi,startj,stopj,dpenv,perc_done)
 *
 * Descrip: No Description
 *
 * Arg:              mat [UNKN ] Undocumented argument [StatWise10 *]
 * Arg:           starti [UNKN ] Undocumented argument [int]
 * Arg:            stopi [UNKN ] Undocumented argument [int]
 * Arg:           startj [UNKN ] Undocumented argument [int]
 * Arg:            stopj [UNKN ] Undocumented argument [int]
 * Arg:            dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:        perc_done [UNKN ] Undocumented argument [int]
 *
 */
}    
void run_up_dc_StatWise10(StatWise10 * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done) 
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
          StatWise10_DC_SHADOW_MATRIX(mat,i,j,CODON) = NEGI;     
          StatWise10_DC_SHADOW_MATRIX(mat,i,j,INTRON_0) = NEGI;  
          StatWise10_DC_SHADOW_MATRIX(mat,i,j,INTRON_1) = NEGI;  
          StatWise10_DC_SHADOW_MATRIX(mat,i,j,INTRON_2) = NEGI;  
          continue;  
          } /* end of Is not in envelope */ 
        if( num % 1000 == 0 )    
          log_full_error(REPORT,0,"[%d%%%% done]Before mid-j %5d Cells done %d%%%%",perc_done,stopj,(num*100)/total);    


        /* For state CODON */ 
        /* setting first movement to score */ 
        score = StatWise10_DC_SHADOW_MATRIX(mat,i-1,j-3,CODON) + mat->codon->codon[CSEQ_GENOMIC_CODON(mat->seq,j)];  
        /* From state CODON to state CODON */ 
        temp = StatWise10_DC_SHADOW_MATRIX(mat,i-0,j-3,CODON) + (mat->codon->codon[CSEQ_GENOMIC_CODON(mat->seq,j)]+mat->exonmodel->exon[i]->stay_score);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTRON_0 to state CODON */ 
        temp = StatWise10_DC_SHADOW_MATRIX(mat,i-1,j-3,INTRON_0) + CSEQ_GENOMIC_3SS(mat->seq,(j-3));     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTRON_1 to state CODON */ 
        temp = StatWise10_DC_SHADOW_MATRIX(mat,i-1,j-2,INTRON_1) + CSEQ_GENOMIC_3SS(mat->seq,(j-2));     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTRON_2 to state CODON */ 
        temp = StatWise10_DC_SHADOW_MATRIX(mat,i-1,j-1,INTRON_2) + CSEQ_GENOMIC_3SS(mat->seq,(j-1));     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for CODON */ 
        /* Add any movement independant score and put away */ 
         StatWise10_DC_SHADOW_MATRIX(mat,i,j,CODON) = score; 
        /* Finished calculating state CODON */ 


        /* For state INTRON_0 */ 
        /* setting first movement to score */ 
        score = StatWise10_DC_SHADOW_MATRIX(mat,i-0,j-1,CODON) + (CSEQ_GENOMIC_5SS(mat->seq,j)+mat->intron_open);    
        /* From state INTRON_0 to state INTRON_0 */ 
        temp = StatWise10_DC_SHADOW_MATRIX(mat,i-0,j-1,INTRON_0) + 0;    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for INTRON_0 */ 
        /* Add any movement independant score and put away */ 
         StatWise10_DC_SHADOW_MATRIX(mat,i,j,INTRON_0) = score;  
        /* Finished calculating state INTRON_0 */ 


        /* For state INTRON_1 */ 
        /* setting first movement to score */ 
        score = StatWise10_DC_SHADOW_MATRIX(mat,i-0,j-2,CODON) + (CSEQ_GENOMIC_5SS(mat->seq,j)+mat->intron_open);    
        /* From state INTRON_1 to state INTRON_1 */ 
        temp = StatWise10_DC_SHADOW_MATRIX(mat,i-0,j-1,INTRON_1) + 0;    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for INTRON_1 */ 
        /* Add any movement independant score and put away */ 
         StatWise10_DC_SHADOW_MATRIX(mat,i,j,INTRON_1) = score;  
        /* Finished calculating state INTRON_1 */ 


        /* For state INTRON_2 */ 
        /* setting first movement to score */ 
        score = StatWise10_DC_SHADOW_MATRIX(mat,i-0,j-3,CODON) + (CSEQ_GENOMIC_5SS(mat->seq,j)+mat->intron_open);    
        /* From state INTRON_2 to state INTRON_2 */ 
        temp = StatWise10_DC_SHADOW_MATRIX(mat,i-0,j-1,INTRON_2) + 0;    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for INTRON_2 */ 
        /* Add any movement independant score and put away */ 
         StatWise10_DC_SHADOW_MATRIX(mat,i,j,INTRON_2) = score;  
        /* Finished calculating state INTRON_2 */ 
        } /* end of this is strip */ 
      } /* end of for each valid j column */ 


/* Function:  init_dc_StatWise10(mat)
 *
 * Descrip: No Description
 *
 * Arg:        mat [UNKN ] Undocumented argument [StatWise10 *]
 *
 */
}    
void init_dc_StatWise10(StatWise10 * mat) 
{
    register int i;  
    register int j;  
    register int k;  


    for(j=0;j<5;j++) {  
      for(i=(-1);i<mat->exonmodel->len;i++)  {  
        StatWise10_DC_SHADOW_MATRIX(mat,i,j,CODON) = NEGI;   
        StatWise10_DC_SHADOW_MATRIX(mat,i,j,INTRON_0) = NEGI;    
        StatWise10_DC_SHADOW_MATRIX(mat,i,j,INTRON_1) = NEGI;    
        StatWise10_DC_SHADOW_MATRIX(mat,i,j,INTRON_2) = NEGI;    
        for(k=0;k<7;k++) {  
          StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,CODON,k) = (-1);    
          StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_0,k) = (-1); 
          StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_1,k) = (-1); 
          StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_2,k) = (-1); 
          }  
        }  
      }  


    return;  
}    


/* Function:  dc_start_end_calculate_StatWise10(mat,dpenv)
 *
 * Descrip:    Calculates special strip, leaving start/end/score points in the shadow matrix 
 *             One tricky thing is that we need to add score-independent calcs in the states
 *             As we have to evaluate them then. This is not ideally implemented therefore 
 *             In fact it is *definitely* not ideal. Will have to do for now
 *
 *
 * Arg:          mat [UNKN ] Undocumented argument [StatWise10 *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean dc_start_end_calculate_StatWise10(StatWise10 * mat,DPEnvelope * dpenv) 
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


    init_start_end_linear_StatWise10(mat);   


    leni = mat->exonmodel->len;  
    lenj = mat->seq->length; 
    total = leni * lenj; 


    for(j=0;j<lenj;j++)  { /*for each j strip*/ 
      for(i=0;i<leni;i++)    { /*for each i position in strip*/ 
        num++;   
        if( dpenv != NULL && is_in_DPEnvelope(dpenv,i,j) == FALSE )  { /*Is not in envelope*/ 
          StatWise10_DC_SHADOW_MATRIX(mat,i,j,CODON) = NEGI;     
          StatWise10_DC_SHADOW_MATRIX(mat,i,j,INTRON_0) = NEGI;  
          StatWise10_DC_SHADOW_MATRIX(mat,i,j,INTRON_1) = NEGI;  
          StatWise10_DC_SHADOW_MATRIX(mat,i,j,INTRON_2) = NEGI;  
          continue;  
          } /* end of Is not in envelope */ 
        if( num%1000 == 0)   
          log_full_error(REPORT,0,"%6d Cells done [%2d%%%%]",num,num*100/total); 




        /* For state CODON */ 
        /* setting first movement to score */ 
        score = StatWise10_DC_SHADOW_MATRIX(mat,i-1,j-3,CODON) + mat->codon->codon[CSEQ_GENOMIC_CODON(mat->seq,j)] + (0);    
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = StatWise10_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 3,CODON,k);  
        /* From state RND_SEQ to state CODON */ 
        temp = StatWise10_DC_SHADOW_SPECIAL(mat,i-1,j-3,RND_SEQ) + (mat->codon->codon[CSEQ_GENOMIC_CODON(mat->seq,j)]+mat->gene_open) + (0);     
        if( temp  > score )  {  
          score = temp;  
          /* This state [RND_SEQ] is a special for CODON... push top shadow pointers here */ 
          localshadow[0]= i; 
          localshadow[1]= j; 
          localshadow[2]= CODON; 
          localshadow[3]= (-1);  
          localshadow[4]= (-1);  
          localshadow[5]= (-1);  
          localshadow[6]= score; 
          }  
        /* From state CODON to state CODON */ 
        temp = StatWise10_DC_SHADOW_MATRIX(mat,i-0,j-3,CODON) + (mat->codon->codon[CSEQ_GENOMIC_CODON(mat->seq,j)]+mat->exonmodel->exon[i]->stay_score) +(0);    
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = StatWise10_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 3,CODON,k);    
          }  
        /* From state INTRON_0 to state CODON */ 
        temp = StatWise10_DC_SHADOW_MATRIX(mat,i-1,j-3,INTRON_0) + CSEQ_GENOMIC_3SS(mat->seq,(j-3)) +(0);    
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = StatWise10_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 3,INTRON_0,k); 
          }  
        /* From state INTRON_1 to state CODON */ 
        temp = StatWise10_DC_SHADOW_MATRIX(mat,i-1,j-2,INTRON_1) + CSEQ_GENOMIC_3SS(mat->seq,(j-2)) +(0);    
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = StatWise10_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 2,INTRON_1,k); 
          }  
        /* From state INTRON_2 to state CODON */ 
        temp = StatWise10_DC_SHADOW_MATRIX(mat,i-1,j-1,INTRON_2) + CSEQ_GENOMIC_3SS(mat->seq,(j-1)) +(0);    
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = StatWise10_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,INTRON_2,k); 
          }  


        /* Ok - finished max calculation for CODON */ 
        /* Add any movement independant score and put away */ 
        /* Actually, already done inside scores */ 
         StatWise10_DC_SHADOW_MATRIX(mat,i,j,CODON) = score; 
        for(k=0;k<7;k++) 
          StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,CODON,k) = localshadow[k];  
        /* Now figure out if any specials need this score */ 


        /* state CODON is a source for special RND_SEQ */ 
        temp = score + (mat->exonmodel->exon[i]->exit_score) + (0) ;     
        if( temp > StatWise10_DC_SHADOW_SPECIAL(mat,i,j,RND_SEQ) )   {  
          StatWise10_DC_SHADOW_SPECIAL(mat,i,j,RND_SEQ) = temp;  
          /* Have to push only bottem half of system here */ 
          for(k=0;k<3;k++)   
            StatWise10_DC_SHADOW_SPECIAL_SP(mat,i,j,RND_SEQ,k) = StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,CODON,k);    
          StatWise10_DC_SHADOW_SPECIAL_SP(mat,i,j,RND_SEQ,6) = StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,CODON,6);  
          StatWise10_DC_SHADOW_SPECIAL_SP(mat,i,j,RND_SEQ,3) = i;    
          StatWise10_DC_SHADOW_SPECIAL_SP(mat,i,j,RND_SEQ,4) = j;    
          StatWise10_DC_SHADOW_SPECIAL_SP(mat,i,j,RND_SEQ,5) = CODON;    
          }  




        /* Finished calculating state CODON */ 


        /* For state INTRON_0 */ 
        /* setting first movement to score */ 
        score = StatWise10_DC_SHADOW_MATRIX(mat,i-0,j-1,CODON) + (CSEQ_GENOMIC_5SS(mat->seq,j)+mat->intron_open) + (0);  
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = StatWise10_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,CODON,k);  
        /* From state INTRON_0 to state INTRON_0 */ 
        temp = StatWise10_DC_SHADOW_MATRIX(mat,i-0,j-1,INTRON_0) + 0 +(0);   
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = StatWise10_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,INTRON_0,k); 
          }  


        /* Ok - finished max calculation for INTRON_0 */ 
        /* Add any movement independant score and put away */ 
        /* Actually, already done inside scores */ 
         StatWise10_DC_SHADOW_MATRIX(mat,i,j,INTRON_0) = score;  
        for(k=0;k<7;k++) 
          StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_0,k) = localshadow[k];   
        /* Now figure out if any specials need this score */ 


        /* Finished calculating state INTRON_0 */ 


        /* For state INTRON_1 */ 
        /* setting first movement to score */ 
        score = StatWise10_DC_SHADOW_MATRIX(mat,i-0,j-2,CODON) + (CSEQ_GENOMIC_5SS(mat->seq,j)+mat->intron_open) + (0);  
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = StatWise10_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 2,CODON,k);  
        /* From state INTRON_1 to state INTRON_1 */ 
        temp = StatWise10_DC_SHADOW_MATRIX(mat,i-0,j-1,INTRON_1) + 0 +(0);   
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = StatWise10_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,INTRON_1,k); 
          }  


        /* Ok - finished max calculation for INTRON_1 */ 
        /* Add any movement independant score and put away */ 
        /* Actually, already done inside scores */ 
         StatWise10_DC_SHADOW_MATRIX(mat,i,j,INTRON_1) = score;  
        for(k=0;k<7;k++) 
          StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_1,k) = localshadow[k];   
        /* Now figure out if any specials need this score */ 


        /* Finished calculating state INTRON_1 */ 


        /* For state INTRON_2 */ 
        /* setting first movement to score */ 
        score = StatWise10_DC_SHADOW_MATRIX(mat,i-0,j-3,CODON) + (CSEQ_GENOMIC_5SS(mat->seq,j)+mat->intron_open) + (0);  
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = StatWise10_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 3,CODON,k);  
        /* From state INTRON_2 to state INTRON_2 */ 
        temp = StatWise10_DC_SHADOW_MATRIX(mat,i-0,j-1,INTRON_2) + 0 +(0);   
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = StatWise10_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,INTRON_2,k); 
          }  


        /* Ok - finished max calculation for INTRON_2 */ 
        /* Add any movement independant score and put away */ 
        /* Actually, already done inside scores */ 
         StatWise10_DC_SHADOW_MATRIX(mat,i,j,INTRON_2) = score;  
        for(k=0;k<7;k++) 
          StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_2,k) = localshadow[k];   
        /* Now figure out if any specials need this score */ 


        /* Finished calculating state INTRON_2 */ 


        } /* end of for each i position in strip */ 


      /* Special state RND_SEQ has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = StatWise10_DC_SHADOW_SPECIAL(mat,0,j,RND_SEQ); 


      /* Source CODON for state RND_SEQ is not special... already calculated */ 
      /* Source RND_SEQ is a special source for RND_SEQ */ 
      temp = StatWise10_DC_SHADOW_SPECIAL(mat,0,j - 1,RND_SEQ) + (0) + (0);  
      if( temp > score ) {  
        score = temp;    
        /* Also got to propagate shadows  */ 
        for(k=0;k<7;k++) 
          StatWise10_DC_SHADOW_SPECIAL_SP(mat,i,j,RND_SEQ,k) = StatWise10_DC_SHADOW_SPECIAL_SP(mat,i - 0,j - 1,RND_SEQ,k);   
        }  


      /* Source START is a special source for RND_SEQ */ 
      temp = StatWise10_DC_SHADOW_SPECIAL(mat,0,j - 1,START) + (0) + (0);    
      if( temp > score ) {  
        score = temp;    
        /* Also got to propagate shadows  */ 
        for(k=0;k<7;k++) 
          StatWise10_DC_SHADOW_SPECIAL_SP(mat,i,j,RND_SEQ,k) = StatWise10_DC_SHADOW_SPECIAL_SP(mat,i - 0,j - 1,START,k); 
        }  


      /* Put back score... (now updated!) */ 
      StatWise10_DC_SHADOW_SPECIAL(mat,0,j,RND_SEQ) = score; 
      /* Finished updating state RND_SEQ */ 




      /* Special state START has no special to special movements */ 


      /* Special state END has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = StatWise10_DC_SHADOW_SPECIAL(mat,0,j,END); 


      /* Source RND_SEQ is a special source for END */ 
      temp = StatWise10_DC_SHADOW_SPECIAL(mat,0,j - 1,RND_SEQ) + (0) + (0);  
      if( temp > score ) {  
        score = temp;    
        /* Also got to propagate shadows  */ 
        for(k=0;k<7;k++) 
          StatWise10_DC_SHADOW_SPECIAL_SP(mat,i,j,END,k) = StatWise10_DC_SHADOW_SPECIAL_SP(mat,i - 0,j - 1,RND_SEQ,k);   
        }  


      /* Put back score... (now updated!) */ 
      StatWise10_DC_SHADOW_SPECIAL(mat,0,j,END) = score; 
      /* Finished updating state END */ 


      } /* end of for each j strip */ 
    return TRUE;     
}    


/* Function:  start_end_find_end_StatWise10(mat,endj)
 *
 * Descrip:    First function used to find end of the best path in the special state !end
 *
 *
 * Arg:         mat [UNKN ] Matrix in small mode [StatWise10 *]
 * Arg:        endj [WRITE] position of end in j (meaningless in i) [int *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int start_end_find_end_StatWise10(StatWise10 * mat,int * endj) 
{
    register int j;  
    register int max;    
    register int maxj;   


    max = StatWise10_DC_SHADOW_SPECIAL(mat,0,mat->seq->length-1,END);    
    maxj = mat->seq->length-1;   
    for(j= mat->seq->length-2 ;j >= 0 ;j--)  {  
      if( StatWise10_DC_SHADOW_SPECIAL(mat,0,j,END) > max )  {  
        max = StatWise10_DC_SHADOW_SPECIAL(mat,0,j,END); 
        maxj = j;    
        }  
      }  


    if( endj != NULL)    
      *endj = maxj;  


    return max;  
}    


/* Function:  init_start_end_linear_StatWise10(mat)
 *
 * Descrip: No Description
 *
 * Arg:        mat [UNKN ] Undocumented argument [StatWise10 *]
 *
 */
void init_start_end_linear_StatWise10(StatWise10 * mat) 
{
    register int i;  
    register int j;  
    for(j=0;j<5;j++) {  
      for(i=(-1);i<mat->exonmodel->len;i++)  {  
        StatWise10_DC_SHADOW_MATRIX(mat,i,j,CODON) = NEGI;   
        StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,CODON,0) = (-1);  
        StatWise10_DC_SHADOW_MATRIX(mat,i,j,INTRON_0) = NEGI;    
        StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_0,0) = (-1);   
        StatWise10_DC_SHADOW_MATRIX(mat,i,j,INTRON_1) = NEGI;    
        StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_1,0) = (-1);   
        StatWise10_DC_SHADOW_MATRIX(mat,i,j,INTRON_2) = NEGI;    
        StatWise10_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_2,0) = (-1);   
        }  
      }  


    for(j=(-3);j<mat->seq->length;j++)   {  
      StatWise10_DC_SHADOW_SPECIAL(mat,0,j,RND_SEQ) = NEGI;  
      StatWise10_DC_SHADOW_SPECIAL_SP(mat,0,j,RND_SEQ,0) = (-1); 
      StatWise10_DC_SHADOW_SPECIAL(mat,0,j,START) = 0;   
      StatWise10_DC_SHADOW_SPECIAL_SP(mat,0,j,START,0) = j;  
      StatWise10_DC_SHADOW_SPECIAL(mat,0,j,END) = NEGI;  
      StatWise10_DC_SHADOW_SPECIAL_SP(mat,0,j,END,0) = (-1); 
      }  


    return;  
}    


/* Function:  convert_PackAln_to_AlnBlock_StatWise10(pal)
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
AlnBlock * convert_PackAln_to_AlnBlock_StatWise10(PackAln * pal) 
{
    AlnConvertSet * acs; 
    AlnBlock * alb;  


    acs = AlnConvertSet_StatWise10();    
    alb = AlnBlock_from_PackAln(acs,pal);    
    free_AlnConvertSet(acs); 
    return alb;  
}    


 static char * query_label[] = { "EXON_STATE","INTRON_STATE","LOOP_STATE","END" };   
/* Function:  AlnConvertSet_StatWise10(void)
 *
 * Descrip: No Description
 *
 *
 * Return [UNKN ]  Undocumented return value [AlnConvertSet *]
 *
 */
 static char * target_label[] = { "CODON","3SS_PHASE_0","3SS_PHASE_1","3SS_PHASE_2","5SS_PHASE_0","CENTRAL_INTRON","5SS_PHASE_1","5SS_PHASE_2","RANDOM_SEQUENCE","END" };    
AlnConvertSet * AlnConvertSet_StatWise10(void) 
{
    AlnConvertUnit * acu;    
    AlnConvertSet  * out;    


    out = AlnConvertSet_alloc_std(); 


    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = CODON; 
    acu->state2 = CODON;     
    acu->offi = 1;   
    acu->offj = 3;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[0];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = RND_SEQ + 4;   
    acu->is_from_special = TRUE; 
    acu->state2 = CODON;     
    acu->offi = (-1);    
    acu->offj = 3;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[0];   
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
    acu->state1 = INTRON_0;  
    acu->state2 = CODON;     
    acu->offi = 1;   
    acu->offj = 3;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[1];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = INTRON_1;  
    acu->state2 = CODON;     
    acu->offi = 1;   
    acu->offj = 2;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[2];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = INTRON_2;  
    acu->state2 = CODON;     
    acu->offi = 1;   
    acu->offj = 1;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[3];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = CODON; 
    acu->state2 = INTRON_0;  
    acu->offi = 0;   
    acu->offj = 1;   
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
    acu->offj = 2;   
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
    acu->offj = 3;   
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
    acu->state1 = CODON; 
    acu->state2 = RND_SEQ + 4;   
    acu->offi = (-1);    
    acu->offj = 0;   
    acu->label1 = query_label[2];    
    acu->label2 = target_label[8];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = RND_SEQ + 4;   
    acu->state2 = RND_SEQ + 4;   
    acu->offi = (-1);    
    acu->offj = 1;   
    acu->label1 = query_label[2];    
    acu->label2 = target_label[8];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = START + 4; 
    acu->state2 = RND_SEQ + 4;   
    acu->offi = (-1);    
    acu->offj = 1;   
    acu->label1 = query_label[2];    
    acu->label2 = target_label[8];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = RND_SEQ + 4;   
    acu->state2 = END + 4;   
    acu->offi = (-1);    
    acu->offj = 1;   
    acu->label1 = query_label[3];    
    acu->label2 = target_label[9];   
    add_collapse_label_AlnConvertSet(out,"LOOP_STATE","RANDOM_SEQUENCE");    
    add_collapse_label_AlnConvertSet(out,"NON_CDS_CONSERVED","RANDOM_SEQUENCE"); 
    add_collapse_label_AlnConvertSet(out,"INTRON_STATE","CENTRAL_INTRON");   
    return out;  
}    


/* Function:  PackAln_read_Expl_StatWise10(mat)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [StatWise10 *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * PackAln_read_Expl_StatWise10(StatWise10 * mat) 
{
    register PackAln * out;  
    int i;   
    int j;   
    int state;   
    int cellscore = (-1);    
    boolean isspecial;   
    PackAlnUnit * pau = NULL;    
    PackAlnUnit * prev = NULL;   


    if( mat->basematrix->type != BASEMATRIX_TYPE_EXPLICIT)   {  
      warn("In StatWise10_basic_read You have asked for an alignment from a non-explicit matrix: c'est impossible [current type is %d - %s]", mat->basematrix->type,basematrix_type_to_string(mat->basematrix->type));   
      return NULL;   
      }  


    out = PackAln_alloc_std();   
    if( out == NULL )    
      return NULL;   


    out->score =  find_end_StatWise10(mat,&i,&j,&state,&isspecial);  


    /* Add final end transition (at the moment we have not got the score! */ 
    if( (pau= PackAlnUnit_alloc()) == NULL  || add_PackAln(out,pau) == FALSE )   {  
      warn("Failed the first PackAlnUnit alloc, %d length of Alignment in StatWise10_basic_read, returning a mess.(Sorry!)",out->len);   
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
        max_calc_special_StatWise10(mat,i,j,state,isspecial,&i,&j,&state,&isspecial,&cellscore);     
      else   
        max_calc_StatWise10(mat,i,j,state,isspecial,&i,&j,&state,&isspecial,&cellscore);     
      if(i == StatWise10_READ_OFF_ERROR || j == StatWise10_READ_OFF_ERROR || state == StatWise10_READ_OFF_ERROR )    {  
        warn("Problem - hit bad read off system, exiting now");  
        break;   
        }  
      if( (pau= PackAlnUnit_alloc()) == NULL  || add_PackAln(out,pau) == FALSE ) {  
        warn("Failed a PackAlnUnit alloc, %d length of Alignment in StatWise10_basic_read, returning partial alignment",out->len);   
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


/* Function:  find_end_StatWise10(mat,ri,rj,state,isspecial)
 *
 * Descrip: No Description
 *
 * Arg:              mat [UNKN ] Undocumented argument [StatWise10 *]
 * Arg:               ri [UNKN ] Undocumented argument [int *]
 * Arg:               rj [UNKN ] Undocumented argument [int *]
 * Arg:            state [UNKN ] Undocumented argument [int *]
 * Arg:        isspecial [UNKN ] Undocumented argument [boolean *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int find_end_StatWise10(StatWise10 * mat,int * ri,int * rj,int * state,boolean * isspecial) 
{
    register int j;  
    register int max;    
    register int maxj;   


    max = StatWise10_EXPL_SPECIAL(mat,0,mat->seq->length-1,END); 
    maxj = mat->seq->length-1;   
    for(j= mat->seq->length-2 ;j >= 0 ;j--)  {  
      if( StatWise10_EXPL_SPECIAL(mat,0,j,END) > max )   {  
        max = StatWise10_EXPL_SPECIAL(mat,0,j,END);  
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


/* Function:  StatWise10_debug_show_matrix(mat,starti,stopi,startj,stopj,ofp)
 *
 * Descrip: No Description
 *
 * Arg:           mat [UNKN ] Undocumented argument [StatWise10 *]
 * Arg:        starti [UNKN ] Undocumented argument [int]
 * Arg:         stopi [UNKN ] Undocumented argument [int]
 * Arg:        startj [UNKN ] Undocumented argument [int]
 * Arg:         stopj [UNKN ] Undocumented argument [int]
 * Arg:           ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
void StatWise10_debug_show_matrix(StatWise10 * mat,int starti,int stopi,int startj,int stopj,FILE * ofp) 
{
    register int i;  
    register int j;  


    for(i=starti;i<stopi && i < mat->exonmodel->len;i++) {  
      for(j=startj;j<stopj && j < mat->seq->length;j++)  {  
        fprintf(ofp,"Cell [%d - %d]\n",i,j);     
        fprintf(ofp,"State CODON %d\n",StatWise10_EXPL_MATRIX(mat,i,j,CODON));   
        fprintf(ofp,"State INTRON_0 %d\n",StatWise10_EXPL_MATRIX(mat,i,j,INTRON_0)); 
        fprintf(ofp,"State INTRON_1 %d\n",StatWise10_EXPL_MATRIX(mat,i,j,INTRON_1)); 
        fprintf(ofp,"State INTRON_2 %d\n",StatWise10_EXPL_MATRIX(mat,i,j,INTRON_2)); 
        fprintf(ofp,"\n\n"); 
        }  
      }  


}    


/* Function:  max_calc_StatWise10(mat,i,j,state,isspecial,reti,retj,retstate,retspecial,cellscore)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [StatWise10 *]
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
int max_calc_StatWise10(StatWise10 * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore) 
{
    register int temp;   
    register int cscore; 


    *reti = (*retj) = (*retstate) = StatWise10_READ_OFF_ERROR;   


    if( i < 0 || j < 0 || i > mat->exonmodel->len || j > mat->seq->length)   {  
      warn("In StatWise10 matrix special read off - out of bounds on matrix [i,j is %d,%d state %d in standard matrix]",i,j,state);  
      return -1;     
      }  


    /* Then you have to select the correct switch statement to figure out the readoff      */ 
    /* Somewhat odd - reverse the order of calculation and return as soon as it is correct */ 
    cscore = StatWise10_EXPL_MATRIX(mat,i,j,state);  
    switch(state)    { /*Switch state */ 
      case CODON :   
        temp = cscore - (CSEQ_GENOMIC_3SS(mat->seq,(j-1))) -  (0);   
        if( temp == StatWise10_EXPL_MATRIX(mat,i - 1,j - 1,INTRON_2) )   {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = INTRON_2;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - StatWise10_EXPL_MATRIX(mat,i-1,j-1,INTRON_2);  
            }  
          return StatWise10_EXPL_MATRIX(mat,i - 1,j - 1,INTRON_2);   
          }  
        temp = cscore - (CSEQ_GENOMIC_3SS(mat->seq,(j-2))) -  (0);   
        if( temp == StatWise10_EXPL_MATRIX(mat,i - 1,j - 2,INTRON_1) )   {  
          *reti = i - 1; 
          *retj = j - 2; 
          *retstate = INTRON_1;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - StatWise10_EXPL_MATRIX(mat,i-1,j-2,INTRON_1);  
            }  
          return StatWise10_EXPL_MATRIX(mat,i - 1,j - 2,INTRON_1);   
          }  
        temp = cscore - (CSEQ_GENOMIC_3SS(mat->seq,(j-3))) -  (0);   
        if( temp == StatWise10_EXPL_MATRIX(mat,i - 1,j - 3,INTRON_0) )   {  
          *reti = i - 1; 
          *retj = j - 3; 
          *retstate = INTRON_0;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - StatWise10_EXPL_MATRIX(mat,i-1,j-3,INTRON_0);  
            }  
          return StatWise10_EXPL_MATRIX(mat,i - 1,j - 3,INTRON_0);   
          }  
        temp = cscore - ((mat->codon->codon[CSEQ_GENOMIC_CODON(mat->seq,j)]+mat->exonmodel->exon[i]->stay_score)) -  (0);    
        if( temp == StatWise10_EXPL_MATRIX(mat,i - 0,j - 3,CODON) )  {  
          *reti = i - 0; 
          *retj = j - 3; 
          *retstate = CODON; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - StatWise10_EXPL_MATRIX(mat,i-0,j-3,CODON); 
            }  
          return StatWise10_EXPL_MATRIX(mat,i - 0,j - 3,CODON);  
          }  
        /* Has restricted position */ 
        if( (i-1) == 0  )    {  
          temp = cscore - ((mat->codon->codon[CSEQ_GENOMIC_CODON(mat->seq,j)]+mat->gene_open)) -  (0);   
          if( temp == StatWise10_EXPL_SPECIAL(mat,i - 1,j - 3,RND_SEQ) ) {  
            *reti = i - 1;   
            *retj = j - 3;   
            *retstate = RND_SEQ; 
            *retspecial = TRUE;  
            if( cellscore != NULL)   {  
              *cellscore = cscore - StatWise10_EXPL_SPECIAL(mat,i-1,j-3,RND_SEQ);    
              }  
            return StatWise10_EXPL_MATRIX(mat,i - 1,j - 3,RND_SEQ);  
            }  
          }  
        temp = cscore - (mat->codon->codon[CSEQ_GENOMIC_CODON(mat->seq,j)]) -  (0);  
        if( temp == StatWise10_EXPL_MATRIX(mat,i - 1,j - 3,CODON) )  {  
          *reti = i - 1; 
          *retj = j - 3; 
          *retstate = CODON; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - StatWise10_EXPL_MATRIX(mat,i-1,j-3,CODON); 
            }  
          return StatWise10_EXPL_MATRIX(mat,i - 1,j - 3,CODON);  
          }  
        warn("Major problem (!) - in StatWise10 read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      case INTRON_0 :    
        temp = cscore - (0) -  (0);  
        if( temp == StatWise10_EXPL_MATRIX(mat,i - 0,j - 1,INTRON_0) )   {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = INTRON_0;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - StatWise10_EXPL_MATRIX(mat,i-0,j-1,INTRON_0);  
            }  
          return StatWise10_EXPL_MATRIX(mat,i - 0,j - 1,INTRON_0);   
          }  
        temp = cscore - ((CSEQ_GENOMIC_5SS(mat->seq,j)+mat->intron_open)) -  (0);    
        if( temp == StatWise10_EXPL_MATRIX(mat,i - 0,j - 1,CODON) )  {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = CODON; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - StatWise10_EXPL_MATRIX(mat,i-0,j-1,CODON); 
            }  
          return StatWise10_EXPL_MATRIX(mat,i - 0,j - 1,CODON);  
          }  
        warn("Major problem (!) - in StatWise10 read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      case INTRON_1 :    
        temp = cscore - (0) -  (0);  
        if( temp == StatWise10_EXPL_MATRIX(mat,i - 0,j - 1,INTRON_1) )   {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = INTRON_1;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - StatWise10_EXPL_MATRIX(mat,i-0,j-1,INTRON_1);  
            }  
          return StatWise10_EXPL_MATRIX(mat,i - 0,j - 1,INTRON_1);   
          }  
        temp = cscore - ((CSEQ_GENOMIC_5SS(mat->seq,j)+mat->intron_open)) -  (0);    
        if( temp == StatWise10_EXPL_MATRIX(mat,i - 0,j - 2,CODON) )  {  
          *reti = i - 0; 
          *retj = j - 2; 
          *retstate = CODON; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - StatWise10_EXPL_MATRIX(mat,i-0,j-2,CODON); 
            }  
          return StatWise10_EXPL_MATRIX(mat,i - 0,j - 2,CODON);  
          }  
        warn("Major problem (!) - in StatWise10 read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      case INTRON_2 :    
        temp = cscore - (0) -  (0);  
        if( temp == StatWise10_EXPL_MATRIX(mat,i - 0,j - 1,INTRON_2) )   {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = INTRON_2;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - StatWise10_EXPL_MATRIX(mat,i-0,j-1,INTRON_2);  
            }  
          return StatWise10_EXPL_MATRIX(mat,i - 0,j - 1,INTRON_2);   
          }  
        temp = cscore - ((CSEQ_GENOMIC_5SS(mat->seq,j)+mat->intron_open)) -  (0);    
        if( temp == StatWise10_EXPL_MATRIX(mat,i - 0,j - 3,CODON) )  {  
          *reti = i - 0; 
          *retj = j - 3; 
          *retstate = CODON; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - StatWise10_EXPL_MATRIX(mat,i-0,j-3,CODON); 
            }  
          return StatWise10_EXPL_MATRIX(mat,i - 0,j - 3,CODON);  
          }  
        warn("Major problem (!) - in StatWise10 read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      default:   
        warn("Major problem (!) - in StatWise10 read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      } /* end of Switch state  */ 
}    


/* Function:  max_calc_special_StatWise10(mat,i,j,state,isspecial,reti,retj,retstate,retspecial,cellscore)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [StatWise10 *]
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
int max_calc_special_StatWise10(StatWise10 * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore) 
{
    register int temp;   
    register int cscore; 


    *reti = (*retj) = (*retstate) = StatWise10_READ_OFF_ERROR;   


    if( j < 0 || j > mat->seq->length)   {  
      warn("In StatWise10 matrix special read off - out of bounds on matrix [j is %d in special]",j);    
      return -1;     
      }  


    cscore = StatWise10_EXPL_SPECIAL(mat,i,j,state); 
    switch(state)    { /*switch on special states*/ 
      case RND_SEQ :     
        /* source START is a special */ 
        temp = cscore - (0) - (0);   
        if( temp == StatWise10_EXPL_SPECIAL(mat,i - 0,j - 1,START) ) {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = START; 
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - StatWise10_EXPL_SPECIAL(mat,i-0,j-1,START);    
            }  
          return StatWise10_EXPL_MATRIX(mat,i - 0,j - 1,START) ;     
          }  
        /* source RND_SEQ is a special */ 
        temp = cscore - (0) - (0);   
        if( temp == StatWise10_EXPL_SPECIAL(mat,i - 0,j - 1,RND_SEQ) )   {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = RND_SEQ;   
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - StatWise10_EXPL_SPECIAL(mat,i-0,j-1,RND_SEQ);  
            }  
          return StatWise10_EXPL_MATRIX(mat,i - 0,j - 1,RND_SEQ) ;   
          }  
        /* source CODON is from main matrix */ 
        for(i= mat->exonmodel->len-1;i >= 0 ;i--)    { /*for i >= 0*/ 
          temp = cscore - (mat->exonmodel->exon[i]->exit_score) - (0);   
          if( temp == StatWise10_EXPL_MATRIX(mat,i - 0,j - 0,CODON) )    {  
            *reti = i - 0;   
            *retj = j - 0;   
            *retstate = CODON;   
            *retspecial = FALSE; 
            if( cellscore != NULL)   {  
              *cellscore = cscore - StatWise10_EXPL_MATRIX(mat,i-0,j-0,CODON);   
              }  
            return StatWise10_EXPL_MATRIX(mat,i - 0,j - 0,CODON) ;   
            }  
          } /* end of for i >= 0 */ 
      case START :   
      case END :     
        /* source RND_SEQ is a special */ 
        temp = cscore - (0) - (0);   
        if( temp == StatWise10_EXPL_SPECIAL(mat,i - 0,j - 1,RND_SEQ) )   {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = RND_SEQ;   
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - StatWise10_EXPL_SPECIAL(mat,i-0,j-1,RND_SEQ);  
            }  
          return StatWise10_EXPL_MATRIX(mat,i - 0,j - 1,RND_SEQ) ;   
          }  
      default:   
        warn("Major problem (!) - in StatWise10 read off, position %d,%d state %d no source found  dropped into default on source switch!",i,j,state);   
        return (-1); 
      } /* end of switch on special states */ 
}    


/* Function:  calculate_StatWise10(mat)
 *
 * Descrip:    This function calculates the StatWise10 matrix when in explicit mode
 *             To allocate the matrix use /allocate_Expl_StatWise10
 *
 *
 * Arg:        mat [UNKN ] StatWise10 which contains explicit basematrix memory [StatWise10 *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean calculate_StatWise10(StatWise10 * mat) 
{
    int i;   
    int j;   
    int leni;    
    int lenj;    
    int tot; 
    int num; 


    if( mat->basematrix->type != BASEMATRIX_TYPE_EXPLICIT )  {  
      warn("in calculate_StatWise10, passed a non Explicit matrix type, cannot calculate!"); 
      return FALSE;  
      }  


    leni = mat->leni;    
    lenj = mat->lenj;    
    tot = leni * lenj;   
    num = 0; 


    start_reporting("StatWise10 Matrix calculation: ");  
    for(j=0;j<lenj;j++)  {  
      auto int score;    
      auto int temp;     
      for(i=0;i<leni;i++)    {  
        if( num%1000 == 0 )  
          log_full_error(REPORT,0,"[%7d] Cells %2d%%%%",num,num*100/tot);    
        num++;   


        /* For state CODON */ 
        /* setting first movement to score */ 
        score = StatWise10_EXPL_MATRIX(mat,i-1,j-3,CODON) + mat->codon->codon[CSEQ_GENOMIC_CODON(mat->seq,j)];   
        /* Has restricted position */ 
        if( (i-1) == 0  )    {  
          /* From state RND_SEQ to state CODON */ 
          temp = StatWise10_EXPL_SPECIAL(mat,i-1,j-3,RND_SEQ) + (mat->codon->codon[CSEQ_GENOMIC_CODON(mat->seq,j)]+mat->gene_open);  
          if( temp  > score )    {  
            score = temp;    
            }  
          }  
        /* From state CODON to state CODON */ 
        temp = StatWise10_EXPL_MATRIX(mat,i-0,j-3,CODON) + (mat->codon->codon[CSEQ_GENOMIC_CODON(mat->seq,j)]+mat->exonmodel->exon[i]->stay_score);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTRON_0 to state CODON */ 
        temp = StatWise10_EXPL_MATRIX(mat,i-1,j-3,INTRON_0) + CSEQ_GENOMIC_3SS(mat->seq,(j-3));  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTRON_1 to state CODON */ 
        temp = StatWise10_EXPL_MATRIX(mat,i-1,j-2,INTRON_1) + CSEQ_GENOMIC_3SS(mat->seq,(j-2));  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTRON_2 to state CODON */ 
        temp = StatWise10_EXPL_MATRIX(mat,i-1,j-1,INTRON_2) + CSEQ_GENOMIC_3SS(mat->seq,(j-1));  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for CODON */ 
        /* Add any movement independant score and put away */ 
         StatWise10_EXPL_MATRIX(mat,i,j,CODON) = score;  


        /* state CODON is a source for special RND_SEQ */ 
        temp = score + (mat->exonmodel->exon[i]->exit_score) + (0) ;     
        if( temp > StatWise10_EXPL_SPECIAL(mat,i,j,RND_SEQ) )    {  
          StatWise10_EXPL_SPECIAL(mat,i,j,RND_SEQ) = temp;   
          }  




        /* Finished calculating state CODON */ 


        /* For state INTRON_0 */ 
        /* setting first movement to score */ 
        score = StatWise10_EXPL_MATRIX(mat,i-0,j-1,CODON) + (CSEQ_GENOMIC_5SS(mat->seq,j)+mat->intron_open);     
        /* From state INTRON_0 to state INTRON_0 */ 
        temp = StatWise10_EXPL_MATRIX(mat,i-0,j-1,INTRON_0) + 0;     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for INTRON_0 */ 
        /* Add any movement independant score and put away */ 
         StatWise10_EXPL_MATRIX(mat,i,j,INTRON_0) = score;   


        /* Finished calculating state INTRON_0 */ 


        /* For state INTRON_1 */ 
        /* setting first movement to score */ 
        score = StatWise10_EXPL_MATRIX(mat,i-0,j-2,CODON) + (CSEQ_GENOMIC_5SS(mat->seq,j)+mat->intron_open);     
        /* From state INTRON_1 to state INTRON_1 */ 
        temp = StatWise10_EXPL_MATRIX(mat,i-0,j-1,INTRON_1) + 0;     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for INTRON_1 */ 
        /* Add any movement independant score and put away */ 
         StatWise10_EXPL_MATRIX(mat,i,j,INTRON_1) = score;   


        /* Finished calculating state INTRON_1 */ 


        /* For state INTRON_2 */ 
        /* setting first movement to score */ 
        score = StatWise10_EXPL_MATRIX(mat,i-0,j-3,CODON) + (CSEQ_GENOMIC_5SS(mat->seq,j)+mat->intron_open);     
        /* From state INTRON_2 to state INTRON_2 */ 
        temp = StatWise10_EXPL_MATRIX(mat,i-0,j-1,INTRON_2) + 0;     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for INTRON_2 */ 
        /* Add any movement independant score and put away */ 
         StatWise10_EXPL_MATRIX(mat,i,j,INTRON_2) = score;   


        /* Finished calculating state INTRON_2 */ 
        }  


      /* Special state RND_SEQ has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = StatWise10_EXPL_SPECIAL(mat,0,j,RND_SEQ);  


      /* Source CODON for state RND_SEQ is not special... already calculated */ 
      /* Source RND_SEQ is a special source for RND_SEQ */ 
      temp = StatWise10_EXPL_SPECIAL(mat,0,j - 1,RND_SEQ) + (0) + (0);   
      if( temp > score ) 
        score = temp;    


      /* Source START is a special source for RND_SEQ */ 
      temp = StatWise10_EXPL_SPECIAL(mat,0,j - 1,START) + (0) + (0);     
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      StatWise10_EXPL_SPECIAL(mat,0,j,RND_SEQ) = score;  
      /* Finished updating state RND_SEQ */ 




      /* Special state START has no special to special movements */ 


      /* Special state END has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = StatWise10_EXPL_SPECIAL(mat,0,j,END);  


      /* Source RND_SEQ is a special source for END */ 
      temp = StatWise10_EXPL_SPECIAL(mat,0,j - 1,RND_SEQ) + (0) + (0);   
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      StatWise10_EXPL_SPECIAL(mat,0,j,END) = score;  
      /* Finished updating state END */ 


      }  
    stop_reporting();    
    return TRUE;     
}    


/* Function:  StatWise10_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [StatWise10 *]
 *
 */
StatWise10 * StatWise10_alloc(void) 
{
    StatWise10 * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(StatWise10 *) ckalloc (sizeof(StatWise10))) == NULL)    {  
      warn("StatWise10_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
    out->basematrix = NULL;  
    out->leni = 0;   
    out->lenj = 0;   


    return out;  
}    


/* Function:  free_StatWise10(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [StatWise10 *]
 *
 * Return [UNKN ]  Undocumented return value [StatWise10 *]
 *
 */
StatWise10 * free_StatWise10(StatWise10 * obj) 
{


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a StatWise10 obj. Should be trappable");    
      return NULL;   
      }  


    if( obj->dynamite_hard_link > 1)     {  
      obj->dynamite_hard_link--; 
      return NULL;   
      }  
    if( obj->basematrix != NULL) 
      free_BaseMatrix(obj->basematrix);  
    /* obj->exonmodel is linked in */ 
    /* obj->seq is linked in */ 
    /* obj->codon is linked in */ 
    /* obj->intron_open is linked in */ 
    /* obj->gene_open is linked in */ 


    ckfree(obj); 
    return NULL; 
}    





#ifdef _cplusplus
}
#endif
