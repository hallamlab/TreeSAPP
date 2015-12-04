#ifdef _cplusplus
extern "C" {
#endif
#include "dnaprofiledp.h"


# line 94 "dnaprofiledp.dy"
DnaProfileMatchScore * new_ALLR_DnaProfileMatchScore(DnaProfile * q,DnaProfile * t)
{
  DnaProfileMatchScore * out;
  int i;
  int j;
  int k;
  int length_q;
  int length_t;
  double score;
  double total;


  assert(q != NULL);
  assert(t != NULL);
  assert(q->folded_random == TRUE);
  assert(t->folded_random == TRUE);


  length_q = q->sa->len;
  length_t = t->sa->len;
  out = DnaProfileMatchScore_alloc_matrix(q->len,t->len);

  assert(out != NULL);

  for(i=0;i<q->len;i++) {
    for(j=0;j<t->len;j++) {      
      /*
	Sum_all_bases( count_bt.ln(fr_bq/rnd_b) + count_bq.ln(fr_bt/rnd_b) ) /
	Sum_all_bases( count_bt + count_bq )
      */
      score = 0.0;
      total = 0.0;
      for(k=0;k<4;k++) {
	score += (length_t * t->col[j]->emit[k]) * (log(q->col[i]->emit[k]));
	score += (length_q * q->col[i]->emit[k]) * (log(t->col[j]->emit[k]));
	total += (length_t * t->col[j]->emit[k]) + (length_q * q->col[i]->emit[k]);
      }

      score = score / total;
      
      out->score[i][j] = (int)(score * (double)INTEGER_FACTOR);
      
    }
  }

  return out;
}

# line 142 "dnaprofiledp.dy"
DnaProfileMatchScore * new_DnaProfileMatchScore(DnaProfileScore * q,DnaProfileScore * t)
{
  int i,j;
  int k;
  DnaProfileMatchScore * out;
  int score;

  assert(q != NULL);
  assert(t != NULL);

  out = DnaProfileMatchScore_alloc_matrix(q->len,t->len);

  
  for(i=0;i<q->len;i++) {
    for(j=0;j<t->len;j++) {
      score  =0;
      for(k=0;k<4;k++) {
	if( q->col[i]->emit[k] > 0 || t->col[j]->emit[k] > 0 ) {
	  score += ((q->col[i]->emit[k] + t->col[j]->emit[k]) / 4);
	}
      }
      out->score[i][j] = score;
    }
  }

  return out;
}


# line 85 "dnaprofiledp.c"
/* Function:  DnaProfileMatchScore_alloc_matrix(leni,lenj)
 *
 * Descrip:    Allocates structure and matrix
 *
 *
 * Arg:        leni [UNKN ] Length of first dimension of matrix [int]
 * Arg:        lenj [UNKN ] Length of second dimension of matrix [int]
 *
 * Return [UNKN ]  Undocumented return value [DnaProfileMatchScore *]
 *
 */
DnaProfileMatchScore * DnaProfileMatchScore_alloc_matrix(int leni,int lenj) 
{
    DnaProfileMatchScore * out; /* out is exported          */ 
    register int i; /* for stepping down matrix */ 
    register int j; /* for stepping across matrix */ 


    /* Call alloc function, return NULL if NULL */ 
    if((out = DnaProfileMatchScore_alloc()) == NULL) 
      return NULL;   


    /* Allocate memory for score  */ 
    if((out->score = (Score **) ckcalloc (leni,sizeof(Score *))) == NULL)    {  
      warn("Memory allocation problem in matrix for DnaProfileMatchScore score, first pointer set"); 
      ckfree(out);   
      return NULL;   
      }  


    /* Add NULL to all matrix pointers so free can be called */ 
    for(i=0;i<leni;i++)  
      out->score[i] = NULL;  


    /* Allocate each matrix row */ 
    for(i=0;i<leni;i++)  {  
      out->score[i] = (Score *) ckcalloc (lenj,sizeof(Score ));  
      if( out->score[i] == NULL) {  
        warn("Failed alloc on %d, calling free and returning NULL",i);   
        free_DnaProfileMatchScore(out);  
        return NULL;     
        }  
      }  


    for(i=0;i<leni;i++)  {  
      for(j=0;j<lenj;j++)    
        out->score[i][j] = NULL;     
      }  


    out->leni=out->maxleni=leni;     
    out->lenj=out->maxlenj=lenj;     


    return out;  
}    


/* Function:  expand_DnaProfileMatchScore(obj,leni,lenj)
 *
 * Descrip:    Expands matrix. Rarely used
 *
 *
 * Arg:         obj [UNKN ] Undocumented argument [DnaProfileMatchScore *]
 * Arg:        leni [UNKN ] Undocumented argument [int]
 * Arg:        lenj [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_DnaProfileMatchScore(DnaProfileMatchScore * obj,int leni,int lenj) 
{
    int i;   
    int actualj;     


    if( obj == NULL)     {  
      warn("Trying to expand a DnaProfileMatchScore but is NULL!");  
      return FALSE;  
      }  


    if( leni <= obj->maxleni && lenj <= obj->maxlenj)    
      return TRUE;   


    if( obj->maxleni < leni )    {  
      if( (obj->score=(Score **) ckrealloc (obj->score,sizeof(Score *)*leni)) == NULL)   
        return FALSE;    
      obj->maxleni=obj->leni=leni;   
      }  
    if( lenj > obj->maxlenj )    
      actualj = lenj;    
    else actualj = obj->maxlenj; 
    for(i=0;i<obj->leni;i++) {  
      if((obj->score[i] = (Score *) realloc (obj->score[i],sizeof(Score ) * actualj)) == NULL)   
        return FALSE;    
      }  
    return TRUE;     
}    


/* Function:  hard_link_DnaProfileMatchScore(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [DnaProfileMatchScore *]
 *
 * Return [UNKN ]  Undocumented return value [DnaProfileMatchScore *]
 *
 */
DnaProfileMatchScore * hard_link_DnaProfileMatchScore(DnaProfileMatchScore * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a DnaProfileMatchScore object: passed a NULL object");    
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  DnaProfileMatchScore_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DnaProfileMatchScore *]
 *
 */
DnaProfileMatchScore * DnaProfileMatchScore_alloc(void) 
{
    DnaProfileMatchScore * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(DnaProfileMatchScore *) ckalloc (sizeof(DnaProfileMatchScore))) == NULL)    {  
      warn("DnaProfileMatchScore_alloc failed ");    
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->score = NULL;   
    out->leni=out->maxleni=0;    
    out->lenj=out->maxlenj=0;    


    return out;  
}    


/* Function:  free_DnaProfileMatchScore(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [DnaProfileMatchScore *]
 *
 * Return [UNKN ]  Undocumented return value [DnaProfileMatchScore *]
 *
 */
DnaProfileMatchScore * free_DnaProfileMatchScore(DnaProfileMatchScore * obj) 
{
    int return_early = 0;    
    int i;   
    int j;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a DnaProfileMatchScore obj. Should be trappable");  
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
    if( obj->score != NULL)  {  
      for(i=0;i<obj->leni;i++)   {  
        if( obj->score[i] != NULL)   
          ckfree(obj->score[i]); 
        }  
      ckfree(obj->score);    
      }  


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
#define MATCH 0  
#define UNMATCHED 1  
#define INSERT 2 
#define DELETE 3 


#define START 0  
#define END 1    


#define DnaProfileMat_EXPL_MATRIX(this_matrix,i,j,STATE) this_matrix->basematrix->matrix[((j+1)*4)+STATE][i+1]   
#define DnaProfileMat_EXPL_SPECIAL(matrix,i,j,STATE) matrix->basematrix->specmatrix[STATE][j+1]  
#define DnaProfileMat_READ_OFF_ERROR -3
 


#define DnaProfileMat_VSMALL_MATRIX(mat,i,j,STATE) mat->basematrix->matrix[(j+2)%2][((i+1)*4)+STATE] 
#define DnaProfileMat_VSMALL_SPECIAL(mat,i,j,STATE) mat->basematrix->specmatrix[(j+2)%2][STATE]  




#define DnaProfileMat_SHATTER_SPECIAL(matrix,i,j,STATE) matrix->shatter->special[STATE][j]   
#define DnaProfileMat_SHATTER_MATRIX(matrix,i,j,STATE)  fetch_cell_value_ShatterMatrix(mat->shatter,i,j,STATE)   


/* Function:  PackAln_read_Shatter_DnaProfileMat(mat)
 *
 * Descrip:    Reads off PackAln from shatter matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [DnaProfileMat *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * PackAln_read_Shatter_DnaProfileMat(DnaProfileMat * mat) 
{
    DnaProfileMat_access_func_holder holder;     


    holder.access_main    = DnaProfileMat_shatter_access_main;   
    holder.access_special = DnaProfileMat_shatter_access_special;    
    assert(mat);     
    assert(mat->shatter);    
    return PackAln_read_generic_DnaProfileMat(mat,holder);   
}    


/* Function:  DnaProfileMat_shatter_access_main(mat,i,j,state)
 *
 * Descrip: No Description
 *
 * Arg:          mat [UNKN ] Undocumented argument [DnaProfileMat *]
 * Arg:            i [UNKN ] Undocumented argument [int]
 * Arg:            j [UNKN ] Undocumented argument [int]
 * Arg:        state [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int DnaProfileMat_shatter_access_main(DnaProfileMat * mat,int i,int j,int state) 
{
    return DnaProfileMat_SHATTER_MATRIX(mat,i,j,state);  
}    


/* Function:  DnaProfileMat_shatter_access_special(mat,i,j,state)
 *
 * Descrip: No Description
 *
 * Arg:          mat [UNKN ] Undocumented argument [DnaProfileMat *]
 * Arg:            i [UNKN ] Undocumented argument [int]
 * Arg:            j [UNKN ] Undocumented argument [int]
 * Arg:        state [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int DnaProfileMat_shatter_access_special(DnaProfileMat * mat,int i,int j,int state) 
{
    return DnaProfileMat_SHATTER_SPECIAL(mat,i,j,state); 
}    


/* Function:  calculate_shatter_DnaProfileMat(mat,dpenv)
 *
 * Descrip:    This function calculates the DnaProfileMat matrix when in shatter mode
 *
 *
 * Arg:          mat [UNKN ] (null) [DnaProfileMat *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean calculate_shatter_DnaProfileMat(DnaProfileMat * mat,DPEnvelope * dpenv) 
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


    start_reporting("DnaProfileMat Matrix calculation: ");   
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
        score = SIG_1_1[MATCH] + 0;  
        /* From state UNMATCHED to state MATCH */ 
        temp = SIG_1_1[UNMATCHED] + 0;   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state MATCH */ 
        temp = SIG_1_1[INSERT] + 0;  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state MATCH */ 
        temp = SIG_1_1[DELETE] + 0;  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state START to state MATCH */ 
        temp = DnaProfileMat_SHATTER_SPECIAL(mat,i-1,j-1,START) + 0;     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for MATCH */ 
        /* Add any movement independant score and put away */ 
         score += mat->m->score[i][j];   
         SIG_0_0[MATCH] = score; 


        /* state MATCH is a source for special END */ 
        temp = score + (0) + (0) ;   
        if( temp > DnaProfileMat_SHATTER_SPECIAL(mat,i,j,END) )  {  
          DnaProfileMat_SHATTER_SPECIAL(mat,i,j,END) = temp;     
          }  




        /* Finished calculating state MATCH */ 


        /* For state UNMATCHED */ 
        /* setting first movement to score */ 
        score = SIG_1_1[MATCH] + mat->open_unmatched;    
        /* From state UNMATCHED to state UNMATCHED */ 
        temp = SIG_1_1[UNMATCHED] + mat->ext_unmatched;  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state UNMATCHED to state UNMATCHED */ 
        temp = SIG_0_1[UNMATCHED] + (mat->ext_unmatched+mat->gap_unmatched);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state UNMATCHED to state UNMATCHED */ 
        temp = SIG_1_0[UNMATCHED] + (mat->ext_unmatched+mat->gap_unmatched);     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for UNMATCHED */ 
        /* Add any movement independant score and put away */ 
         SIG_0_0[UNMATCHED] = score; 


        /* Finished calculating state UNMATCHED */ 


        /* For state INSERT */ 
        /* setting first movement to score */ 
        score = SIG_0_1[MATCH] + (mat->q->col[i]->trans[DnaProfile_M2I]+mat->t->col[j]->trans[DnaProfile_M2D]);  
        /* From state INSERT to state INSERT */ 
        temp = SIG_0_1[INSERT] + (mat->q->col[i]->trans[DnaProfile_I2I]+mat->t->col[j]->trans[DnaProfile_D2D]);  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for INSERT */ 
        /* Add any movement independant score and put away */ 
         SIG_0_0[INSERT] = score;    


        /* Finished calculating state INSERT */ 


        /* For state DELETE */ 
        /* setting first movement to score */ 
        score = SIG_1_0[MATCH] + (mat->t->col[j]->trans[DnaProfile_M2I]+mat->q->col[i]->trans[DnaProfile_M2D]);  
        /* From state INSERT to state DELETE */ 
        temp = SIG_1_0[INSERT] + (mat->t->col[j]->trans[DnaProfile_I2I]+mat->q->col[i]->trans[DnaProfile_D2D]);  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for DELETE */ 
        /* Add any movement independant score and put away */ 
         SIG_0_0[DELETE] = score;    


        /* Finished calculating state DELETE */ 
        }  


      /* Special state START has no special to special movements */ 


      /* Special state END has no special to special movements */ 
      }  
    stop_reporting();    
    return TRUE;     
}    


/* Function:  search_DnaProfileMat(dbsi,out,q,t,m,open_unmatched,ext_unmatched,gap_unmatched)
 *
 * Descrip:    This function makes a database search of DnaProfileMat
 *             It uses the dbsi structure to choose which implementation
 *             to use of the database searching. This way at run time you
 *             can switch between single threaded/multi-threaded or hardware
 *
 *
 * Arg:                  dbsi [UNKN ] Undocumented argument [DBSearchImpl *]
 * Arg:                   out [UNKN ] Undocumented argument [Hscore *]
 * Arg:                     q [UNKN ] Undocumented argument [DnaProfileScore *]
 * Arg:                     t [UNKN ] Undocumented argument [DnaProfileScore *]
 * Arg:                     m [UNKN ] Undocumented argument [DnaProfileMatchScore*]
 * Arg:        open_unmatched [UNKN ] Undocumented argument [Score]
 * Arg:         ext_unmatched [UNKN ] Undocumented argument [Score]
 * Arg:         gap_unmatched [UNKN ] Undocumented argument [Score]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type search_DnaProfileMat(DBSearchImpl * dbsi,Hscore * out,DnaProfileScore * q,DnaProfileScore * t ,DnaProfileMatchScore* m,Score open_unmatched,Score ext_unmatched,Score gap_unmatched) 
{
#ifdef PTHREAD   
    int i;   
    int thr_no;  
    pthread_attr_t pat;  
    struct thread_pool_holder_DnaProfileMat * holder;    
#endif   
    if( out == NULL )    {  
      warn("Passed in a null Hscore object into search_DnaProfileMat. Can't process results!");  
      return SEARCH_ERROR;   
      }  
    if( dbsi == NULL )   {  
      warn("Passed in a null DBSearchImpl object into search_DnaProfileMat. Can't process results!");    
      return SEARCH_ERROR;   
      }  
    if( dbsi->trace_level > 5 )  
      warn("Asking for trace level of %d in database search for DnaProfileMat, but it was compiled with a trace level of -2139062144. Not all trace statements can be shown",dbsi->trace_level); 
    switch(dbsi->type)   { /*switch on implementation*/ 
      case DBSearchImpl_Serial : 
        return serial_search_DnaProfileMat(out,q,t ,m,open_unmatched,ext_unmatched,gap_unmatched);   
      case DBSearchImpl_Pthreads :   
#ifdef PTHREAD   
        holder = (struct thread_pool_holder_DnaProfileMat *) ckalloc(sizeof(struct thread_pool_holder_DnaProfileMat));   
        if( holder == NULL )     {  
          warn("Unable to allocated thread pool datastructure...");  
          return SEARCH_ERROR;   
          }  
        holder->out = out;   
        holder->dbsi = dbsi; 
        holder->q = q;   
        holder->t = t;   
        holder->m = m;   
        holder->open_unmatched = open_unmatched; 
        holder->ext_unmatched = ext_unmatched;   
        holder->gap_unmatched = gap_unmatched;   
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
          if( pthread_create(holder->pool+i,&pat,thread_loop_DnaProfileMat,(void *)holder) ) 
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
        warn("You did not specifiy the PTHREAD compile when compiled the C code for DnaProfileMat"); 
#endif /* finished threads */    
      default :  
        warn("database search implementation %s was not provided in the compiled dynamite file from DnaProfileMat",impl_string_DBSearchImpl(dbsi));  
        return SEARCH_ERROR; 
      } /* end of switch on implementation */ 


}    


/* Function:  thread_loop_DnaProfileMat(ptr)
 *
 * Descrip:    dummy loop code foreach thread for DnaProfileMat
 *
 *
 * Arg:        ptr [UNKN ] Undocumented argument [void *]
 *
 * Return [UNKN ]  Undocumented return value [void *]
 *
 */
void * thread_loop_DnaProfileMat(void * ptr) 
{
    fatal("dummy thread loop function"); 
}    


/* Function:  serial_search_DnaProfileMat(out,q,t,m,open_unmatched,ext_unmatched,gap_unmatched)
 *
 * Descrip:    This function makes a database search of DnaProfileMat
 *             It is a single processor implementation
 *
 *
 * Arg:                   out [UNKN ] Undocumented argument [Hscore *]
 * Arg:                     q [UNKN ] Undocumented argument [DnaProfileScore *]
 * Arg:                     t [UNKN ] Undocumented argument [DnaProfileScore *]
 * Arg:                     m [UNKN ] Undocumented argument [DnaProfileMatchScore*]
 * Arg:        open_unmatched [UNKN ] Undocumented argument [Score]
 * Arg:         ext_unmatched [UNKN ] Undocumented argument [Score]
 * Arg:         gap_unmatched [UNKN ] Undocumented argument [Score]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type serial_search_DnaProfileMat(Hscore * out,DnaProfileScore * q,DnaProfileScore * t ,DnaProfileMatchScore* m,Score open_unmatched,Score ext_unmatched,Score gap_unmatched) 
{
    int db_status;   
    int score;   
    int query_pos = 0;   
    int target_pos = 0;  
    DataScore * ds;  


    push_errormsg_stack("Before any actual search in db searching"); 


    target_pos = 0;  




    /* No maximum length - allocated on-the-fly */ 
    score = score_only_DnaProfileMat(q, t , m, open_unmatched, ext_unmatched, gap_unmatched);    
    if( should_store_Hscore(out,score) == TRUE )     { /*if storing datascore*/ 
      ds = new_DataScore_from_storage(out);  
      if( ds == NULL )   {  
        warn("DnaProfileMat search had a memory error in allocating a new_DataScore (?a leak somewhere - DataScore is a very small datastructure");  
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


/* Function:  score_only_DnaProfileMat(q,t,m,open_unmatched,ext_unmatched,gap_unmatched)
 *
 * Descrip:    This function just calculates the score for the matrix
 *             I am pretty sure we can do this better, but hey, for the moment...
 *             It calls /allocate_DnaProfileMat_only
 *
 *
 * Arg:                     q [UNKN ] query data structure [DnaProfileScore *]
 * Arg:                     t [UNKN ] target data structure [DnaProfileScore *]
 * Arg:                     m [UNKN ] Resource [DnaProfileMatchScore*]
 * Arg:        open_unmatched [UNKN ] Resource [Score]
 * Arg:         ext_unmatched [UNKN ] Resource [Score]
 * Arg:         gap_unmatched [UNKN ] Resource [Score]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int score_only_DnaProfileMat(DnaProfileScore * q,DnaProfileScore * t ,DnaProfileMatchScore* m,Score open_unmatched,Score ext_unmatched,Score gap_unmatched) 
{
    int bestscore = NEGI;    
    int i;   
    int j;   
    int k;   
    DnaProfileMat * mat;     


    mat = allocate_DnaProfileMat_only(q, t , m, open_unmatched, ext_unmatched, gap_unmatched);   
    if( mat == NULL )    {  
      warn("Memory allocation error in the db search - unable to communicate to calling function. this spells DIASTER!");    
      return NEGI;   
      }  
    if((mat->basematrix = BaseMatrix_alloc_matrix_and_specials(2,(mat->leni + 1) * 4,2,2)) == NULL)  {  
      warn("Score only matrix for DnaProfileMat cannot be allocated, (asking for 1  by %d  cells)",mat->leni*4); 
      mat = free_DnaProfileMat(mat);     
      return 0;  
      }  
    mat->basematrix->type = BASEMATRIX_TYPE_VERYSMALL;   


    /* Now, initiate matrix */ 
    for(j=0;j<3;j++) {  
      for(i=(-1);i<mat->leni;i++)    {  
        for(k=0;k<4;k++) 
          DnaProfileMat_VSMALL_MATRIX(mat,i,j,k) = NEGI; 
        }  
      DnaProfileMat_VSMALL_SPECIAL(mat,i,j,START) = 0;   
      DnaProfileMat_VSMALL_SPECIAL(mat,i,j,END) = NEGI;  
      }  


    /* Ok, lets do-o-o-o-o it */ 


    for(j=0;j<mat->lenj;j++) { /*for all target positions*/ 
      auto int score;    
      auto int temp;     
      for(i=0;i<mat->leni;i++)   { /*for all query positions*/ 


        /* For state MATCH */ 
        /* setting first movement to score */ 
        score = DnaProfileMat_VSMALL_MATRIX(mat,i-1,j-1,MATCH) + 0;  
        /* From state UNMATCHED to state MATCH */ 
        temp = DnaProfileMat_VSMALL_MATRIX(mat,i-1,j-1,UNMATCHED) + 0;   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state MATCH */ 
        temp = DnaProfileMat_VSMALL_MATRIX(mat,i-1,j-1,INSERT) + 0;  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state MATCH */ 
        temp = DnaProfileMat_VSMALL_MATRIX(mat,i-1,j-1,DELETE) + 0;  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state START to state MATCH */ 
        temp = DnaProfileMat_VSMALL_SPECIAL(mat,i-1,j-1,START) + 0;  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for MATCH */ 
        /* Add any movement independant score and put away */ 
         score += mat->m->score[i][j];   
         DnaProfileMat_VSMALL_MATRIX(mat,i,j,MATCH) = score; 


        /* state MATCH is a source for special END */ 
        temp = score + (0) + (0) ;   
        if( temp > DnaProfileMat_VSMALL_SPECIAL(mat,i,j,END) )   {  
          DnaProfileMat_VSMALL_SPECIAL(mat,i,j,END) = temp;  
          }  




        /* Finished calculating state MATCH */ 


        /* For state UNMATCHED */ 
        /* setting first movement to score */ 
        score = DnaProfileMat_VSMALL_MATRIX(mat,i-1,j-1,MATCH) + mat->open_unmatched;    
        /* From state UNMATCHED to state UNMATCHED */ 
        temp = DnaProfileMat_VSMALL_MATRIX(mat,i-1,j-1,UNMATCHED) + mat->ext_unmatched;  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state UNMATCHED to state UNMATCHED */ 
        temp = DnaProfileMat_VSMALL_MATRIX(mat,i-0,j-1,UNMATCHED) + (mat->ext_unmatched+mat->gap_unmatched);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state UNMATCHED to state UNMATCHED */ 
        temp = DnaProfileMat_VSMALL_MATRIX(mat,i-1,j-0,UNMATCHED) + (mat->ext_unmatched+mat->gap_unmatched);     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for UNMATCHED */ 
        /* Add any movement independant score and put away */ 
         DnaProfileMat_VSMALL_MATRIX(mat,i,j,UNMATCHED) = score; 


        /* Finished calculating state UNMATCHED */ 


        /* For state INSERT */ 
        /* setting first movement to score */ 
        score = DnaProfileMat_VSMALL_MATRIX(mat,i-0,j-1,MATCH) + (mat->q->col[i]->trans[DnaProfile_M2I]+mat->t->col[j]->trans[DnaProfile_M2D]);  
        /* From state INSERT to state INSERT */ 
        temp = DnaProfileMat_VSMALL_MATRIX(mat,i-0,j-1,INSERT) + (mat->q->col[i]->trans[DnaProfile_I2I]+mat->t->col[j]->trans[DnaProfile_D2D]);  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for INSERT */ 
        /* Add any movement independant score and put away */ 
         DnaProfileMat_VSMALL_MATRIX(mat,i,j,INSERT) = score;    


        /* Finished calculating state INSERT */ 


        /* For state DELETE */ 
        /* setting first movement to score */ 
        score = DnaProfileMat_VSMALL_MATRIX(mat,i-1,j-0,MATCH) + (mat->t->col[j]->trans[DnaProfile_M2I]+mat->q->col[i]->trans[DnaProfile_M2D]);  
        /* From state INSERT to state DELETE */ 
        temp = DnaProfileMat_VSMALL_MATRIX(mat,i-1,j-0,INSERT) + (mat->t->col[j]->trans[DnaProfile_I2I]+mat->q->col[i]->trans[DnaProfile_D2D]);  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for DELETE */ 
        /* Add any movement independant score and put away */ 
         DnaProfileMat_VSMALL_MATRIX(mat,i,j,DELETE) = score;    


        /* Finished calculating state DELETE */ 
        } /* end of for all query positions */ 




      /* Special state START has no special to special movements */ 


      /* Special state END has no special to special movements */ 
      if( bestscore < DnaProfileMat_VSMALL_SPECIAL(mat,0,j,END) )    
        bestscore = DnaProfileMat_VSMALL_SPECIAL(mat,0,j,END);   
      } /* end of for all target positions */ 


    mat = free_DnaProfileMat(mat);   
    return bestscore;    
}    


/* Function:  PackAln_bestmemory_DnaProfileMat(q,t,m,open_unmatched,ext_unmatched,gap_unmatched,dpenv,dpri)
 *
 * Descrip:    This function chooses the best memory set-up for the alignment
 *             using calls to basematrix, and then implements either a large
 *             or small memory model.
 *
 *             It is the best function to use if you just want an alignment
 *
 *             If you want a label alignment, you will need
 *             /convert_PackAln_to_AlnBlock_DnaProfileMat
 *
 *
 * Arg:                     q [UNKN ] query data structure [DnaProfileScore *]
 * Arg:                     t [UNKN ] target data structure [DnaProfileScore *]
 * Arg:                     m [UNKN ] Resource [DnaProfileMatchScore*]
 * Arg:        open_unmatched [UNKN ] Resource [Score]
 * Arg:         ext_unmatched [UNKN ] Resource [Score]
 * Arg:         gap_unmatched [UNKN ] Resource [Score]
 * Arg:                 dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:                  dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * PackAln_bestmemory_DnaProfileMat(DnaProfileScore * q,DnaProfileScore * t ,DnaProfileMatchScore* m,Score open_unmatched,Score ext_unmatched,Score gap_unmatched,DPEnvelope * dpenv,DPRunImpl * dpri) 
{
    long long total; 
    DnaProfileMat * mat; 
    PackAln * out;   
    DebugMatrix * de;    
    DPRunImplMemory strategy;    
    assert(dpri);    


    total = q->len * t->len;     
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
        if( (mat=allocate_Expl_DnaProfileMat(q, t , m, open_unmatched, ext_unmatched, gap_unmatched,dpri)) == NULL ) {  
          warn("Unable to allocate large DnaProfileMat version");    
          return NULL;   
          }  
        calculate_dpenv_DnaProfileMat(mat,dpenv);    
        out =  PackAln_read_Expl_DnaProfileMat(mat); 
        }  
      else   {  
        mat = allocate_DnaProfileMat_only(q, t , m, open_unmatched, ext_unmatched, gap_unmatched);   
        calculate_shatter_DnaProfileMat(mat,dpenv);  
        out = PackAln_read_Shatter_DnaProfileMat(mat);   
        }  
      }  
    else {  
      if( strategy == DPIM_Linear )  {  
        /* use small implementation */ 
        if( (mat=allocate_Small_DnaProfileMat(q, t , m, open_unmatched, ext_unmatched, gap_unmatched)) == NULL ) {  
          warn("Unable to allocate small DnaProfileMat version");    
          return NULL;   
          }  
        out = PackAln_calculate_Small_DnaProfileMat(mat,dpenv);  
        }  
      else   {  
        /* use Large implementation */ 
        if( (mat=allocate_Expl_DnaProfileMat(q, t , m, open_unmatched, ext_unmatched, gap_unmatched,dpri)) == NULL ) {  
          warn("Unable to allocate large DnaProfileMat version");    
          return NULL;   
          }  
        if( dpri->debug == TRUE) {  
          fatal("Asked for dydebug, but dynamite file not compiled with -g. Need to recompile dynamite source"); 
          }  
        else {  
          calculate_DnaProfileMat(mat);  
          out =  PackAln_read_Expl_DnaProfileMat(mat);   
          }  
        }  
      }  


    mat = free_DnaProfileMat(mat);   
    return out;  
}    


/* Function:  allocate_DnaProfileMat_only(q,t,m,open_unmatched,ext_unmatched,gap_unmatched)
 *
 * Descrip:    This function only allocates the DnaProfileMat structure
 *             checks types where possible and determines leni and lenj
 *             The basematrix area is delt with elsewhere
 *
 *
 * Arg:                     q [UNKN ] query data structure [DnaProfileScore *]
 * Arg:                     t [UNKN ] target data structure [DnaProfileScore *]
 * Arg:                     m [UNKN ] Resource [DnaProfileMatchScore*]
 * Arg:        open_unmatched [UNKN ] Resource [Score]
 * Arg:         ext_unmatched [UNKN ] Resource [Score]
 * Arg:         gap_unmatched [UNKN ] Resource [Score]
 *
 * Return [UNKN ]  Undocumented return value [DnaProfileMat *]
 *
 */
DnaProfileMat * allocate_DnaProfileMat_only(DnaProfileScore * q,DnaProfileScore * t ,DnaProfileMatchScore* m,Score open_unmatched,Score ext_unmatched,Score gap_unmatched) 
{
    DnaProfileMat * out;     


    if((out= DnaProfileMat_alloc()) == NULL) {  
      warn("Allocation of basic DnaProfileMat structure failed..."); 
      return NULL;   
      }  


    out->q = q;  
    out->t = t;  
    out->m = m;  
    out->open_unmatched = open_unmatched;    
    out->ext_unmatched = ext_unmatched;  
    out->gap_unmatched = gap_unmatched;  
    out->leni = q->len;  
    out->lenj = t->len;  
    return out;  
}    


/* Function:  allocate_Expl_DnaProfileMat(q,t,m,open_unmatched,ext_unmatched,gap_unmatched,dpri)
 *
 * Descrip:    This function allocates the DnaProfileMat structure
 *             and the basematrix area for explicit memory implementations
 *             It calls /allocate_DnaProfileMat_only
 *
 *
 * Arg:                     q [UNKN ] query data structure [DnaProfileScore *]
 * Arg:                     t [UNKN ] target data structure [DnaProfileScore *]
 * Arg:                     m [UNKN ] Resource [DnaProfileMatchScore*]
 * Arg:        open_unmatched [UNKN ] Resource [Score]
 * Arg:         ext_unmatched [UNKN ] Resource [Score]
 * Arg:         gap_unmatched [UNKN ] Resource [Score]
 * Arg:                  dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [DnaProfileMat *]
 *
 */
DnaProfileMat * allocate_Expl_DnaProfileMat(DnaProfileScore * q,DnaProfileScore * t ,DnaProfileMatchScore* m,Score open_unmatched,Score ext_unmatched,Score gap_unmatched,DPRunImpl * dpri) 
{
    DnaProfileMat * out; 


    out = allocate_DnaProfileMat_only(q, t , m, open_unmatched, ext_unmatched, gap_unmatched);   
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
        warn("Explicit matrix DnaProfileMat cannot be allocated, (asking for %d by %d main cells)",out->leni,out->lenj); 
        free_DnaProfileMat(out);     
        return NULL; 
        }  
      }  
    if( dpri->should_cache == TRUE && dpri->cache == NULL)   
      dpri->cache = hard_link_BaseMatrix(out->basematrix);   
    out->basematrix->type = BASEMATRIX_TYPE_EXPLICIT;    
    init_DnaProfileMat(out);     
    return out;  
}    


/* Function:  init_DnaProfileMat(mat)
 *
 * Descrip:    This function initates DnaProfileMat matrix when in explicit mode
 *             Called in /allocate_Expl_DnaProfileMat
 *
 *
 * Arg:        mat [UNKN ] DnaProfileMat which contains explicit basematrix memory [DnaProfileMat *]
 *
 */
void init_DnaProfileMat(DnaProfileMat * mat) 
{
    register int i;  
    register int j;  
    if( mat->basematrix->type != BASEMATRIX_TYPE_EXPLICIT)   {  
      warn("Cannot iniate matrix, is not an explicit memory type and you have assummed that");   
      return;    
      }  


    for(i= (-1);i<mat->q->len;i++)   {  
      for(j= (-1);j<2;j++)   {  
        DnaProfileMat_EXPL_MATRIX(mat,i,j,MATCH) = NEGI; 
        DnaProfileMat_EXPL_MATRIX(mat,i,j,UNMATCHED) = NEGI; 
        DnaProfileMat_EXPL_MATRIX(mat,i,j,INSERT) = NEGI;    
        DnaProfileMat_EXPL_MATRIX(mat,i,j,DELETE) = NEGI;    
        }  
      }  
    for(j= (-1);j<mat->t->len;j++)   {  
      for(i= (-1);i<2;i++)   {  
        DnaProfileMat_EXPL_MATRIX(mat,i,j,MATCH) = NEGI; 
        DnaProfileMat_EXPL_MATRIX(mat,i,j,UNMATCHED) = NEGI; 
        DnaProfileMat_EXPL_MATRIX(mat,i,j,INSERT) = NEGI;    
        DnaProfileMat_EXPL_MATRIX(mat,i,j,DELETE) = NEGI;    
        }  
      DnaProfileMat_EXPL_SPECIAL(mat,i,j,START) = 0; 
      DnaProfileMat_EXPL_SPECIAL(mat,i,j,END) = NEGI;    
      }  
    return;  
}    


/* Function:  recalculate_PackAln_DnaProfileMat(pal,mat)
 *
 * Descrip:    This function recalculates the PackAln structure produced by DnaProfileMat
 *             For example, in linear space methods this is used to score them
 *
 *
 * Arg:        pal [UNKN ] Undocumented argument [PackAln *]
 * Arg:        mat [UNKN ] Undocumented argument [DnaProfileMat *]
 *
 */
void recalculate_PackAln_DnaProfileMat(PackAln * pal,DnaProfileMat * mat) 
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
            pau->score = 0 + (mat->m->score[i][j]);  
            continue;    
            }  
          if( offi == 1 && offj == 1 && prev->state == UNMATCHED )   {  
            pau->score = 0 + (mat->m->score[i][j]);  
            continue;    
            }  
          if( offi == 1 && offj == 1 && prev->state == INSERT )  {  
            pau->score = 0 + (mat->m->score[i][j]);  
            continue;    
            }  
          if( offi == 1 && offj == 1 && prev->state == DELETE )  {  
            pau->score = 0 + (mat->m->score[i][j]);  
            continue;    
            }  
          if( offj == 1 && prev->state == (START+4) )    {  
            pau->score = 0 + (mat->m->score[i][j]);  
            continue;    
            }  
          warn("In recaluclating PackAln with state MATCH, from [%d,%d,%d], got a bad source state. Error!",offi,offj,prev->state);  
          break; 
        case UNMATCHED :     
          if( offi == 1 && offj == 1 && prev->state == MATCH )   {  
            pau->score = mat->open_unmatched + (0);  
            continue;    
            }  
          if( offi == 1 && offj == 1 && prev->state == UNMATCHED )   {  
            pau->score = mat->ext_unmatched + (0);   
            continue;    
            }  
          if( offi == 0 && offj == 1 && prev->state == UNMATCHED )   {  
            pau->score = (mat->ext_unmatched+mat->gap_unmatched) + (0);  
            continue;    
            }  
          if( offi == 1 && offj == 0 && prev->state == UNMATCHED )   {  
            pau->score = (mat->ext_unmatched+mat->gap_unmatched) + (0);  
            continue;    
            }  
          warn("In recaluclating PackAln with state UNMATCHED, from [%d,%d,%d], got a bad source state. Error!",offi,offj,prev->state);  
          break; 
        case INSERT :    
          if( offi == 0 && offj == 1 && prev->state == MATCH )   {  
            pau->score = (mat->q->col[i]->trans[DnaProfile_M2I]+mat->t->col[j]->trans[DnaProfile_M2D]) + (0);    
            continue;    
            }  
          if( offi == 0 && offj == 1 && prev->state == INSERT )  {  
            pau->score = (mat->q->col[i]->trans[DnaProfile_I2I]+mat->t->col[j]->trans[DnaProfile_D2D]) + (0);    
            continue;    
            }  
          warn("In recaluclating PackAln with state INSERT, from [%d,%d,%d], got a bad source state. Error!",offi,offj,prev->state); 
          break; 
        case DELETE :    
          if( offi == 1 && offj == 0 && prev->state == MATCH )   {  
            pau->score = (mat->t->col[j]->trans[DnaProfile_M2I]+mat->q->col[i]->trans[DnaProfile_M2D]) + (0);    
            continue;    
            }  
          if( offi == 1 && offj == 0 && prev->state == INSERT )  {  
            pau->score = (mat->t->col[j]->trans[DnaProfile_I2I]+mat->q->col[i]->trans[DnaProfile_D2D]) + (0);    
            continue;    
            }  
          warn("In recaluclating PackAln with state DELETE, from [%d,%d,%d], got a bad source state. Error!",offi,offj,prev->state); 
          break; 
        case (START+4) :     
          warn("In recaluclating PackAln with state START, got a bad source state. Error!"); 
          break; 
        case (END+4) :   
          if( offj == 0 && prev->state == MATCH )    {  
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
#define DnaProfileMat_HIDDEN_MATRIX(thismatrix,i,j,state) (thismatrix->basematrix->matrix[(j-hiddenj+1)][(i+1)*4+state]) 
#define DnaProfileMat_DC_SHADOW_MATRIX(thismatrix,i,j,state) (thismatrix->basematrix->matrix[((j+2)*8) % 16][(i+1)*4+state]) 
#define DnaProfileMat_HIDDEN_SPECIAL(thismatrix,i,j,state) (thismatrix->basematrix->specmatrix[state][(j+1)])    
#define DnaProfileMat_DC_SHADOW_SPECIAL(thismatrix,i,j,state) (thismatrix->basematrix->specmatrix[state*8][(j+1)])   
#define DnaProfileMat_DC_SHADOW_MATRIX_SP(thismatrix,i,j,state,shadow) (thismatrix->basematrix->matrix[((((j+2)*8)+(shadow+1)) % 16)][(i+1)*4 + state])  
#define DnaProfileMat_DC_SHADOW_SPECIAL_SP(thismatrix,i,j,state,shadow) (thismatrix->basematrix->specmatrix[state*8 +shadow+1][(j+1)])   
#define DnaProfileMat_DC_OPT_SHADOW_MATRIX(thismatrix,i,j,state) (score_pointers[(((j+1)% 1) * (leni+1) * 4) + ((i+1) * 4) + (state)])   
#define DnaProfileMat_DC_OPT_SHADOW_MATRIX_SP(thismatrix,i,j,state,shadow) (shadow_pointers[(((j+1)% 1) * (leni+1) * 32) + ((i+1) * 32) + (state * 8) + shadow+1])   
#define DnaProfileMat_DC_OPT_SHADOW_SPECIAL(thismatrix,i,j,state) (thismatrix->basematrix->specmatrix[state*8][(j+1)])   
/* Function:  allocate_Small_DnaProfileMat(q,t,m,open_unmatched,ext_unmatched,gap_unmatched)
 *
 * Descrip:    This function allocates the DnaProfileMat structure
 *             and the basematrix area for a small memory implementations
 *             It calls /allocate_DnaProfileMat_only
 *
 *
 * Arg:                     q [UNKN ] query data structure [DnaProfileScore *]
 * Arg:                     t [UNKN ] target data structure [DnaProfileScore *]
 * Arg:                     m [UNKN ] Resource [DnaProfileMatchScore*]
 * Arg:        open_unmatched [UNKN ] Resource [Score]
 * Arg:         ext_unmatched [UNKN ] Resource [Score]
 * Arg:         gap_unmatched [UNKN ] Resource [Score]
 *
 * Return [UNKN ]  Undocumented return value [DnaProfileMat *]
 *
 */
#define DnaProfileMat_DC_OPT_SHADOW_SPECIAL_SP(thismatrix,i,j,state,shadow) (thismatrix->basematrix->specmatrix[state*8 +shadow+1][(j+1)])   
DnaProfileMat * allocate_Small_DnaProfileMat(DnaProfileScore * q,DnaProfileScore * t ,DnaProfileMatchScore* m,Score open_unmatched,Score ext_unmatched,Score gap_unmatched) 
{
    DnaProfileMat * out; 


    out = allocate_DnaProfileMat_only(q, t , m, open_unmatched, ext_unmatched, gap_unmatched);   
    if( out == NULL )    
      return NULL;   
    out->basematrix = BaseMatrix_alloc_matrix_and_specials(16,(out->leni + 1) * 4,16,out->lenj+1);   
    if(out == NULL)  {  
      warn("Small shadow matrix DnaProfileMat cannot be allocated, (asking for 2 by %d main cells)",out->leni+2);    
      free_DnaProfileMat(out);   
      return NULL;   
      }  
    out->basematrix->type = BASEMATRIX_TYPE_SHADOW;  
    return out;  
}    


/* Function:  PackAln_calculate_Small_DnaProfileMat(mat,dpenv)
 *
 * Descrip:    This function calculates an alignment for DnaProfileMat structure in linear space
 *             If you want only the start/end points
 *             use /AlnRangeSet_calculate_Small_DnaProfileMat 
 *
 *             The function basically
 *               finds start/end points 
 *               foreach start/end point 
 *                 calls /full_dc_DnaProfileMat 
 *
 *
 * Arg:          mat [UNKN ] Undocumented argument [DnaProfileMat *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * PackAln_calculate_Small_DnaProfileMat(DnaProfileMat * mat,DPEnvelope * dpenv) 
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
      warn("Could not calculate packaln small for DnaProfileMat due to wrong type of matrix");   
      return NULL;   
      }  


    out = PackAln_alloc_std();   


    start_reporting("Find start end points: ");  
    dc_optimised_start_end_calc_DnaProfileMat(mat,dpenv);    
    score = start_end_find_end_DnaProfileMat(mat,&endj); 
    out->score = score;  
    stopstate = END;
    
    /* No special to specials: one matrix alignment: simply remove and get */ 
    starti = DnaProfileMat_DC_SHADOW_SPECIAL_SP(mat,0,endj,END,0);   
    startj = DnaProfileMat_DC_SHADOW_SPECIAL_SP(mat,0,endj,END,1);   
    startstate = DnaProfileMat_DC_SHADOW_SPECIAL_SP(mat,0,endj,END,2);   
    stopi = DnaProfileMat_DC_SHADOW_SPECIAL_SP(mat,0,endj,END,3);    
    stopj = DnaProfileMat_DC_SHADOW_SPECIAL_SP(mat,0,endj,END,4);    
    stopstate = DnaProfileMat_DC_SHADOW_SPECIAL_SP(mat,0,endj,END,5);    
    temp = DnaProfileMat_DC_SHADOW_SPECIAL_SP(mat,0,endj,END,6); 
    log_full_error(REPORT,0,"[%d,%d][%d,%d] Score %d",starti,startj,stopi,stopj,score);  
    stop_reporting();    
    start_reporting("Recovering alignment: ");   


    /* Figuring how much j we have to align for reporting purposes */ 
    donej = 0;   
    totalj = stopj - startj; 
    full_dc_DnaProfileMat(mat,starti,startj,startstate,stopi,stopj,stopstate,out,&donej,totalj,dpenv);   


    /* Although we have no specials, need to get start. Better to check than assume */ 


    max_matrix_to_special_DnaProfileMat(mat,starti,startj,startstate,temp,&stopi,&stopj,&stopstate,&temp,NULL);  
    if( stopi == DnaProfileMat_READ_OFF_ERROR || stopstate != START )    {  
      warn("Problem in reading off special state system, hit a non start state (or an internal error) in a single alignment mode");  
      invert_PackAln(out);   
      recalculate_PackAln_DnaProfileMat(out,mat);    
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
    recalculate_PackAln_DnaProfileMat(out,mat);  
    return out;  


}    


/* Function:  AlnRangeSet_calculate_Small_DnaProfileMat(mat)
 *
 * Descrip:    This function calculates an alignment for DnaProfileMat structure in linear space
 *             If you want the full alignment, use /PackAln_calculate_Small_DnaProfileMat 
 *             If you have already got the full alignment, but want the range set, use /AlnRangeSet_from_PackAln_DnaProfileMat
 *             If you have got the small matrix but not the alignment, use /AlnRangeSet_from_DnaProfileMat 
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [DnaProfileMat *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * AlnRangeSet_calculate_Small_DnaProfileMat(DnaProfileMat * mat) 
{
    AlnRangeSet * out;   


    start_reporting("Find start end points: ");  
    dc_optimised_start_end_calc_DnaProfileMat(mat,NULL); 
    log_full_error(REPORT,0,"Calculated");   


    out = AlnRangeSet_from_DnaProfileMat(mat);   
    return out;  
}    


/* Function:  AlnRangeSet_from_DnaProfileMat(mat)
 *
 * Descrip:    This function reads off a start/end structure
 *             for DnaProfileMat structure in linear space
 *             If you want the full alignment use
 *             /PackAln_calculate_Small_DnaProfileMat 
 *             If you have not calculated the matrix use
 *             /AlnRange_calculate_Small_DnaProfileMat
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [DnaProfileMat *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * AlnRangeSet_from_DnaProfileMat(DnaProfileMat * mat) 
{
    AlnRangeSet * out;   
    AlnRange * temp; 
    int jpos;    
    int state;   


    if( mat->basematrix->type != BASEMATRIX_TYPE_SHADOW) {  
      warn("Bad error! - non shadow matrix type in AlnRangeSet_from_DnaProfileMat"); 
      return NULL;   
      }  


    out = AlnRangeSet_alloc_std();   
    /* Find the end position */ 
    out->score = start_end_find_end_DnaProfileMat(mat,&jpos);    
    state = END; 


    while( (temp = AlnRange_build_DnaProfileMat(mat,jpos,state,&jpos,&state)) != NULL)   
      add_AlnRangeSet(out,temp); 
    return out;  
}    


/* Function:  AlnRange_build_DnaProfileMat(mat,stopj,stopspecstate,startj,startspecstate)
 *
 * Descrip:    This function calculates a single start/end set in linear space
 *             Really a sub-routine for /AlnRangeSet_from_PackAln_DnaProfileMat
 *
 *
 * Arg:                   mat [UNKN ] Undocumented argument [DnaProfileMat *]
 * Arg:                 stopj [UNKN ] Undocumented argument [int]
 * Arg:         stopspecstate [UNKN ] Undocumented argument [int]
 * Arg:                startj [UNKN ] Undocumented argument [int *]
 * Arg:        startspecstate [UNKN ] Undocumented argument [int *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRange *]
 *
 */
AlnRange * AlnRange_build_DnaProfileMat(DnaProfileMat * mat,int stopj,int stopspecstate,int * startj,int * startspecstate) 
{
    AlnRange * out;  
    int jpos;    
    int state;   


    if( mat->basematrix->type != BASEMATRIX_TYPE_SHADOW) {  
      warn("Bad error! - non shadow matrix type in AlnRangeSet_from_DnaProfileMat"); 
      return NULL;   
      }  


    /* Assumme that we have specials (we should!). Read back along the specials till we have the finish point */ 
    if( read_special_strip_DnaProfileMat(mat,0,stopj,stopspecstate,&jpos,&state,NULL) == FALSE)  {  
      warn("In AlnRanger_build_DnaProfileMat alignment ending at %d, unable to read back specials. Will (evenutally) return a partial range set... BEWARE!",stopj);  
      return NULL;   
      }  
    if( state == START || jpos <= 0) 
      return NULL;   


    out = AlnRange_alloc();  


    out->starti = DnaProfileMat_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,0);    
    out->startj = DnaProfileMat_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,1);    
    out->startstate = DnaProfileMat_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,2);    
    out->stopi = DnaProfileMat_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,3); 
    out->stopj = DnaProfileMat_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,4); 
    out->stopstate = DnaProfileMat_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,5); 
    out->startscore = DnaProfileMat_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,6);    
    out->stopscore = DnaProfileMat_DC_SHADOW_SPECIAL(mat,0,jpos,state);  


    /* Now, we have to figure out where this state came from in the specials */ 
    max_matrix_to_special_DnaProfileMat(mat,out->starti,out->startj,out->startstate,out->startscore,&jpos,startj,startspecstate,&state,NULL);    
    if( jpos == DnaProfileMat_READ_OFF_ERROR)    {  
      warn("In AlnRange_build_DnaProfileMat alignment ending at %d, with aln range between %d-%d in j, unable to find source special, returning this range, but this could get tricky!",stopj,out->startj,out->stopj);   
      return out;    
      }  


    /* Put in the correct score for startstate, from the special */ 
    out->startscore = DnaProfileMat_DC_SHADOW_SPECIAL(mat,0,*startj,*startspecstate);    
    /* The correct j coords have been put into startj, startspecstate... so just return out */ 
    return out;  
}    


/* Function:  read_hidden_DnaProfileMat(mat,starti,startj,startstate,stopi,stopj,stopstate,out)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [DnaProfileMat *]
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
boolean read_hidden_DnaProfileMat(DnaProfileMat * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out) 
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


      max_hidden_DnaProfileMat(mat,startj,i,j,state,isspecial,&i,&j,&state,&isspecial,&cellscore);   


      if( i == DnaProfileMat_READ_OFF_ERROR) {  
        warn("In DnaProfileMat hidden read off, between %d:%d,%d:%d - at got bad read off. Problem!",starti,startj,stopi,stopj); 
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
        warn("In DnaProfileMat hidden read off, between %d:%d,%d:%d - hit start cell, but not in start state. Can't be good!.",starti,startj,stopi,stopj);   
        return FALSE;    
        }  
      }  
    warn("In DnaProfileMat hidden read off, between %d:%d,%d:%d - gone past start cell (now in %d,%d,%d), can't be good news!.",starti,startj,stopi,stopj,i,j,state);    
    return FALSE;    
}    


/* Function:  max_hidden_DnaProfileMat(mat,hiddenj,i,j,state,isspecial,reti,retj,retstate,retspecial,cellscore)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [DnaProfileMat *]
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
int max_hidden_DnaProfileMat(DnaProfileMat * mat,int hiddenj,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore) 
{
    register int temp;   
    register int cscore; 


    *reti = (*retj) = (*retstate) = DnaProfileMat_READ_OFF_ERROR;    


    if( i < 0 || j < 0 || i > mat->q->len || j > mat->t->len)    {  
      warn("In DnaProfileMat matrix special read off - out of bounds on matrix [i,j is %d,%d state %d in standard matrix]",i,j,state);   
      return -1; 
      }  


    /* Then you have to select the correct switch statement to figure out the readoff      */ 
    /* Somewhat odd - reverse the order of calculation and return as soon as it is correct */ 
    cscore = DnaProfileMat_HIDDEN_MATRIX(mat,i,j,state); 
    switch(state)    { /*Switch state */ 
      case MATCH :   
        /* Not allowing special sources.. skipping START */ 
        temp = cscore - (0) -  (mat->m->score[i][j]);    
        if( temp == DnaProfileMat_HIDDEN_MATRIX(mat,i - 1,j - 1,DELETE) )    {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = DELETE;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - DnaProfileMat_HIDDEN_MATRIX(mat,i-1,j-1,DELETE);   
            }  
          return DnaProfileMat_HIDDEN_MATRIX(mat,i - 1,j - 1,DELETE);    
          }  
        temp = cscore - (0) -  (mat->m->score[i][j]);    
        if( temp == DnaProfileMat_HIDDEN_MATRIX(mat,i - 1,j - 1,INSERT) )    {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = INSERT;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - DnaProfileMat_HIDDEN_MATRIX(mat,i-1,j-1,INSERT);   
            }  
          return DnaProfileMat_HIDDEN_MATRIX(mat,i - 1,j - 1,INSERT);    
          }  
        temp = cscore - (0) -  (mat->m->score[i][j]);    
        if( temp == DnaProfileMat_HIDDEN_MATRIX(mat,i - 1,j - 1,UNMATCHED) ) {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = UNMATCHED; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - DnaProfileMat_HIDDEN_MATRIX(mat,i-1,j-1,UNMATCHED);    
            }  
          return DnaProfileMat_HIDDEN_MATRIX(mat,i - 1,j - 1,UNMATCHED);     
          }  
        temp = cscore - (0) -  (mat->m->score[i][j]);    
        if( temp == DnaProfileMat_HIDDEN_MATRIX(mat,i - 1,j - 1,MATCH) ) {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - DnaProfileMat_HIDDEN_MATRIX(mat,i-1,j-1,MATCH);    
            }  
          return DnaProfileMat_HIDDEN_MATRIX(mat,i - 1,j - 1,MATCH);     
          }  
        warn("Major problem (!) - in DnaProfileMat read off, position %d,%d state %d no source found!",i,j,state);   
        return (-1); 
      case UNMATCHED :   
        temp = cscore - ((mat->ext_unmatched+mat->gap_unmatched)) -  (0);    
        if( temp == DnaProfileMat_HIDDEN_MATRIX(mat,i - 1,j - 0,UNMATCHED) ) {  
          *reti = i - 1; 
          *retj = j - 0; 
          *retstate = UNMATCHED; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - DnaProfileMat_HIDDEN_MATRIX(mat,i-1,j-0,UNMATCHED);    
            }  
          return DnaProfileMat_HIDDEN_MATRIX(mat,i - 1,j - 0,UNMATCHED);     
          }  
        temp = cscore - ((mat->ext_unmatched+mat->gap_unmatched)) -  (0);    
        if( temp == DnaProfileMat_HIDDEN_MATRIX(mat,i - 0,j - 1,UNMATCHED) ) {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = UNMATCHED; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - DnaProfileMat_HIDDEN_MATRIX(mat,i-0,j-1,UNMATCHED);    
            }  
          return DnaProfileMat_HIDDEN_MATRIX(mat,i - 0,j - 1,UNMATCHED);     
          }  
        temp = cscore - (mat->ext_unmatched) -  (0); 
        if( temp == DnaProfileMat_HIDDEN_MATRIX(mat,i - 1,j - 1,UNMATCHED) ) {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = UNMATCHED; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - DnaProfileMat_HIDDEN_MATRIX(mat,i-1,j-1,UNMATCHED);    
            }  
          return DnaProfileMat_HIDDEN_MATRIX(mat,i - 1,j - 1,UNMATCHED);     
          }  
        temp = cscore - (mat->open_unmatched) -  (0);    
        if( temp == DnaProfileMat_HIDDEN_MATRIX(mat,i - 1,j - 1,MATCH) ) {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - DnaProfileMat_HIDDEN_MATRIX(mat,i-1,j-1,MATCH);    
            }  
          return DnaProfileMat_HIDDEN_MATRIX(mat,i - 1,j - 1,MATCH);     
          }  
        warn("Major problem (!) - in DnaProfileMat read off, position %d,%d state %d no source found!",i,j,state);   
        return (-1); 
      case INSERT :  
        temp = cscore - ((mat->q->col[i]->trans[DnaProfile_I2I]+mat->t->col[j]->trans[DnaProfile_D2D])) -  (0);  
        if( temp == DnaProfileMat_HIDDEN_MATRIX(mat,i - 0,j - 1,INSERT) )    {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = INSERT;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - DnaProfileMat_HIDDEN_MATRIX(mat,i-0,j-1,INSERT);   
            }  
          return DnaProfileMat_HIDDEN_MATRIX(mat,i - 0,j - 1,INSERT);    
          }  
        temp = cscore - ((mat->q->col[i]->trans[DnaProfile_M2I]+mat->t->col[j]->trans[DnaProfile_M2D])) -  (0);  
        if( temp == DnaProfileMat_HIDDEN_MATRIX(mat,i - 0,j - 1,MATCH) ) {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - DnaProfileMat_HIDDEN_MATRIX(mat,i-0,j-1,MATCH);    
            }  
          return DnaProfileMat_HIDDEN_MATRIX(mat,i - 0,j - 1,MATCH);     
          }  
        warn("Major problem (!) - in DnaProfileMat read off, position %d,%d state %d no source found!",i,j,state);   
        return (-1); 
      case DELETE :  
        temp = cscore - ((mat->t->col[j]->trans[DnaProfile_I2I]+mat->q->col[i]->trans[DnaProfile_D2D])) -  (0);  
        if( temp == DnaProfileMat_HIDDEN_MATRIX(mat,i - 1,j - 0,INSERT) )    {  
          *reti = i - 1; 
          *retj = j - 0; 
          *retstate = INSERT;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - DnaProfileMat_HIDDEN_MATRIX(mat,i-1,j-0,INSERT);   
            }  
          return DnaProfileMat_HIDDEN_MATRIX(mat,i - 1,j - 0,INSERT);    
          }  
        temp = cscore - ((mat->t->col[j]->trans[DnaProfile_M2I]+mat->q->col[i]->trans[DnaProfile_M2D])) -  (0);  
        if( temp == DnaProfileMat_HIDDEN_MATRIX(mat,i - 1,j - 0,MATCH) ) {  
          *reti = i - 1; 
          *retj = j - 0; 
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - DnaProfileMat_HIDDEN_MATRIX(mat,i-1,j-0,MATCH);    
            }  
          return DnaProfileMat_HIDDEN_MATRIX(mat,i - 1,j - 0,MATCH);     
          }  
        warn("Major problem (!) - in DnaProfileMat read off, position %d,%d state %d no source found!",i,j,state);   
        return (-1); 
      default:   
        warn("Major problem (!) - in DnaProfileMat read off, position %d,%d state %d no source found!",i,j,state);   
        return (-1); 
      } /* end of Switch state  */ 
}    


/* Function:  read_special_strip_DnaProfileMat(mat,stopi,stopj,stopstate,startj,startstate,out)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [DnaProfileMat *]
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
boolean read_special_strip_DnaProfileMat(DnaProfileMat * mat,int stopi,int stopj,int stopstate,int * startj,int * startstate,PackAln * out) 
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
    while( j > DnaProfileMat_DC_SHADOW_SPECIAL_SP(mat,i,j,state,4) && state != START)    { /*while more specials to eat up*/ 
      /* Put away current state, if we should */ 
      if(out != NULL)    {  
        pau = PackAlnUnit_alloc();  /* Should deal with memory overflow */ 
        pau->i = i;  
        pau->j = j;  
        pau->state =  state + 4; 
        add_PackAln(out,pau);    
        }  


      max_special_strip_DnaProfileMat(mat,i,j,state,isspecial,&i,&j,&state,&isspecial,&cellscore);   
      if( i == DnaProfileMat_READ_OFF_ERROR) {  
        warn("In special strip read DnaProfileMat, got a bad read off error. Sorry!");   
        return FALSE;    
        }  
      } /* end of while more specials to eat up */ 


    /* check to see we have not gone too far! */ 
    if( state != START && j < DnaProfileMat_DC_SHADOW_SPECIAL_SP(mat,i,j,state,4))   {  
      warn("In special strip read DnaProfileMat, at special [%d] state [%d] overshot!",j,state); 
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


/* Function:  max_special_strip_DnaProfileMat(mat,i,j,state,isspecial,reti,retj,retstate,retspecial,cellscore)
 *
 * Descrip:    A pretty intense internal function. Deals with read-off only in specials
 *
 *
 * Arg:               mat [UNKN ] Undocumented argument [DnaProfileMat *]
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
int max_special_strip_DnaProfileMat(DnaProfileMat * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore) 
{
    int temp;    
    int cscore;  


    *reti = (*retj) = (*retstate) = DnaProfileMat_READ_OFF_ERROR;    
    if( isspecial == FALSE ) {  
      warn("In special strip max function for DnaProfileMat, got a non special start point. Problem! (bad!)");   
      return (-1);   
      }  


    if( j < 0 || j > mat->t->len)    {  
      warn("In DnaProfileMat matrix special read off - out of bounds on matrix [j is %d in special]",j); 
      return -1; 
      }  


    cscore = DnaProfileMat_DC_SHADOW_SPECIAL(mat,i,j,state); 
    switch(state)    { /*switch on special states*/ 
      case START :   
      case END :     
        /* Source MATCH is not a special */ 
      default:   
        warn("Major problem (!) - in DnaProfileMat special strip read off, position %d,%d state %d no source found  dropped into default on source switch!",i,j,state);  
        return (-1); 
      } /* end of switch on special states */ 
}    


/* Function:  max_matrix_to_special_DnaProfileMat(mat,i,j,state,cscore,reti,retj,retstate,retspecial,cellscore)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [DnaProfileMat *]
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
int max_matrix_to_special_DnaProfileMat(DnaProfileMat * mat,int i,int j,int state,int cscore,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore) 
{
    int temp;    
    *reti = (*retj) = (*retstate) = DnaProfileMat_READ_OFF_ERROR;    


    if( j < 0 || j > mat->lenj)  {  
      warn("In DnaProfileMat matrix to special read off - out of bounds on matrix [j is %d in special]",j);  
      return -1; 
      }  


    switch(state)    { /*Switch state */ 
      case MATCH :   
        temp = cscore - (0) -  (mat->m->score[i][j]);    
        if( temp == DnaProfileMat_DC_SHADOW_SPECIAL(mat,i - 1,j - 1,START) ) {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = START; 
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - DnaProfileMat_DC_SHADOW_SPECIAL(mat,i-1,j-1,START);    
            }  
          return DnaProfileMat_DC_SHADOW_MATRIX(mat,i - 1,j - 1,START) ;     
          }  
        /* Source DELETE is not a special, should not get here! */ 
        /* Source INSERT is not a special, should not get here! */ 
        /* Source UNMATCHED is not a special, should not get here! */ 
        /* Source MATCH is not a special, should not get here! */ 
        warn("Major problem (!) - in DnaProfileMat matrix to special read off, position %d,%d state %d no source found!",i,j,state); 
        return (-1); 
      case UNMATCHED :   
        /* Source UNMATCHED is not a special, should not get here! */ 
        /* Source UNMATCHED is not a special, should not get here! */ 
        /* Source UNMATCHED is not a special, should not get here! */ 
        /* Source MATCH is not a special, should not get here! */ 
        warn("Major problem (!) - in DnaProfileMat matrix to special read off, position %d,%d state %d no source found!",i,j,state); 
        return (-1); 
      case INSERT :  
        /* Source INSERT is not a special, should not get here! */ 
        /* Source MATCH is not a special, should not get here! */ 
        warn("Major problem (!) - in DnaProfileMat matrix to special read off, position %d,%d state %d no source found!",i,j,state); 
        return (-1); 
      case DELETE :  
        /* Source INSERT is not a special, should not get here! */ 
        /* Source MATCH is not a special, should not get here! */ 
        warn("Major problem (!) - in DnaProfileMat matrix to special read off, position %d,%d state %d no source found!",i,j,state); 
        return (-1); 
      default:   
        warn("Major problem (!) - in DnaProfileMat read off, position %d,%d state %d no source found!",i,j,state);   
        return (-1); 
      } /* end of Switch state  */ 


}    


/* Function:  calculate_hidden_DnaProfileMat(mat,starti,startj,startstate,stopi,stopj,stopstate,dpenv)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [DnaProfileMat *]
 * Arg:            starti [UNKN ] Undocumented argument [int]
 * Arg:            startj [UNKN ] Undocumented argument [int]
 * Arg:        startstate [UNKN ] Undocumented argument [int]
 * Arg:             stopi [UNKN ] Undocumented argument [int]
 * Arg:             stopj [UNKN ] Undocumented argument [int]
 * Arg:         stopstate [UNKN ] Undocumented argument [int]
 * Arg:             dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 */
void calculate_hidden_DnaProfileMat(DnaProfileMat * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv) 
{
    register int i;  
    register int j;  
    register int score;  
    register int temp;   
    register int hiddenj;    


    hiddenj = startj;    


    init_hidden_DnaProfileMat(mat,starti,startj,stopi,stopj);    


    DnaProfileMat_HIDDEN_MATRIX(mat,starti,startj,startstate) = 0;   


    for(j=startj;j<=stopj;j++)   {  
      for(i=starti;i<=stopi;i++) {  
        /* Should *not* do very first cell as this is the one set to zero in one state! */ 
        if( i == starti && j == startj ) 
          continue;  
        if( dpenv != NULL && is_in_DPEnvelope(dpenv,i,j) == FALSE )  { /*Is not in envelope*/ 
          DnaProfileMat_HIDDEN_MATRIX(mat,i,j,MATCH) = NEGI;     
          DnaProfileMat_HIDDEN_MATRIX(mat,i,j,UNMATCHED) = NEGI;     
          DnaProfileMat_HIDDEN_MATRIX(mat,i,j,INSERT) = NEGI;    
          DnaProfileMat_HIDDEN_MATRIX(mat,i,j,DELETE) = NEGI;    
          continue;  
          } /* end of Is not in envelope */ 


        /* For state MATCH */ 
        /* setting first movement to score */ 
        score = DnaProfileMat_HIDDEN_MATRIX(mat,i-1,j-1,MATCH) + 0;  
        /* From state UNMATCHED to state MATCH */ 
        temp = DnaProfileMat_HIDDEN_MATRIX(mat,i-1,j-1,UNMATCHED) + 0;   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state MATCH */ 
        temp = DnaProfileMat_HIDDEN_MATRIX(mat,i-1,j-1,INSERT) + 0;  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state MATCH */ 
        temp = DnaProfileMat_HIDDEN_MATRIX(mat,i-1,j-1,DELETE) + 0;  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for MATCH */ 
        /* Add any movement independant score and put away */ 
         score += mat->m->score[i][j];   
         DnaProfileMat_HIDDEN_MATRIX(mat,i,j,MATCH) = score; 
        /* Finished calculating state MATCH */ 


        /* For state UNMATCHED */ 
        /* setting first movement to score */ 
        score = DnaProfileMat_HIDDEN_MATRIX(mat,i-1,j-1,MATCH) + mat->open_unmatched;    
        /* From state UNMATCHED to state UNMATCHED */ 
        temp = DnaProfileMat_HIDDEN_MATRIX(mat,i-1,j-1,UNMATCHED) + mat->ext_unmatched;  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state UNMATCHED to state UNMATCHED */ 
        temp = DnaProfileMat_HIDDEN_MATRIX(mat,i-0,j-1,UNMATCHED) + (mat->ext_unmatched+mat->gap_unmatched);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state UNMATCHED to state UNMATCHED */ 
        temp = DnaProfileMat_HIDDEN_MATRIX(mat,i-1,j-0,UNMATCHED) + (mat->ext_unmatched+mat->gap_unmatched);     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for UNMATCHED */ 
        /* Add any movement independant score and put away */ 
         DnaProfileMat_HIDDEN_MATRIX(mat,i,j,UNMATCHED) = score; 
        /* Finished calculating state UNMATCHED */ 


        /* For state INSERT */ 
        /* setting first movement to score */ 
        score = DnaProfileMat_HIDDEN_MATRIX(mat,i-0,j-1,MATCH) + (mat->q->col[i]->trans[DnaProfile_M2I]+mat->t->col[j]->trans[DnaProfile_M2D]);  
        /* From state INSERT to state INSERT */ 
        temp = DnaProfileMat_HIDDEN_MATRIX(mat,i-0,j-1,INSERT) + (mat->q->col[i]->trans[DnaProfile_I2I]+mat->t->col[j]->trans[DnaProfile_D2D]);  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for INSERT */ 
        /* Add any movement independant score and put away */ 
         DnaProfileMat_HIDDEN_MATRIX(mat,i,j,INSERT) = score;    
        /* Finished calculating state INSERT */ 


        /* For state DELETE */ 
        /* setting first movement to score */ 
        score = DnaProfileMat_HIDDEN_MATRIX(mat,i-1,j-0,MATCH) + (mat->t->col[j]->trans[DnaProfile_M2I]+mat->q->col[i]->trans[DnaProfile_M2D]);  
        /* From state INSERT to state DELETE */ 
        temp = DnaProfileMat_HIDDEN_MATRIX(mat,i-1,j-0,INSERT) + (mat->t->col[j]->trans[DnaProfile_I2I]+mat->q->col[i]->trans[DnaProfile_D2D]);  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for DELETE */ 
        /* Add any movement independant score and put away */ 
         DnaProfileMat_HIDDEN_MATRIX(mat,i,j,DELETE) = score;    
        /* Finished calculating state DELETE */ 
        }  
      }  


    return;  
}    


/* Function:  init_hidden_DnaProfileMat(mat,starti,startj,stopi,stopj)
 *
 * Descrip: No Description
 *
 * Arg:           mat [UNKN ] Undocumented argument [DnaProfileMat *]
 * Arg:        starti [UNKN ] Undocumented argument [int]
 * Arg:        startj [UNKN ] Undocumented argument [int]
 * Arg:         stopi [UNKN ] Undocumented argument [int]
 * Arg:         stopj [UNKN ] Undocumented argument [int]
 *
 */
void init_hidden_DnaProfileMat(DnaProfileMat * mat,int starti,int startj,int stopi,int stopj) 
{
    register int i;  
    register int j;  
    register int hiddenj;    


    hiddenj = startj;    
    for(j=(startj-1);j<=stopj;j++)   {  
      for(i=(starti-1);i<=stopi;i++) {  
        DnaProfileMat_HIDDEN_MATRIX(mat,i,j,MATCH) = NEGI;
  
        DnaProfileMat_HIDDEN_MATRIX(mat,i,j,UNMATCHED) = NEGI;
  
        DnaProfileMat_HIDDEN_MATRIX(mat,i,j,INSERT) = NEGI;
 
        DnaProfileMat_HIDDEN_MATRIX(mat,i,j,DELETE) = NEGI;
 
        }  
      }  


    return;  
}    


/* Function:  full_dc_DnaProfileMat(mat,starti,startj,startstate,stopi,stopj,stopstate,out,donej,totalj,dpenv)
 *
 * Descrip:    The main divide-and-conquor routine. Basically, call /PackAln_calculate_small_DnaProfileMat
 *             Not this function, which is pretty hard core. 
 *             Function is given start/end points (in main matrix) for alignment
 *             It does some checks, decides whether start/end in j is small enough for explicit calc
 *               - if yes, calculates it, reads off into PackAln (out), adds the j distance to donej and returns TRUE
 *               - if no,  uses /do_dc_single_pass_DnaProfileMat to get mid-point
 *                          saves midpoint, and calls itself to do right portion then left portion
 *             right then left ensures PackAln is added the 'right' way, ie, back-to-front
 *             returns FALSE on any error, with a warning
 *
 *
 * Arg:               mat [UNKN ] Matrix with small memory implementation [DnaProfileMat *]
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
boolean full_dc_DnaProfileMat(DnaProfileMat * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out,int * donej,int totalj,DPEnvelope * dpenv) 
{
    int lstarti; 
    int lstartj; 
    int lstate;  


    if( mat->basematrix->type != BASEMATRIX_TYPE_SHADOW) {  
      warn("*Very* bad error! - non shadow matrix type in full_dc_DnaProfileMat");   
      return FALSE;  
      }  


    if( starti == -1 || startj == -1 || startstate == -1 || stopi == -1 || stopstate == -1)  {  
      warn("In full dc program, passed bad indices, indices passed were %d:%d[%d] to %d:%d[%d]\n",starti,startj,startstate,stopi,stopj,stopstate);   
      return FALSE;  
      }  


    if( stopj - startj < 5)  {  
      log_full_error(REPORT,0,"[%d,%d][%d,%d] Explicit read off",starti,startj,stopi,stopj);/* Build hidden explicit matrix */ 
      calculate_hidden_DnaProfileMat(mat,starti,startj,startstate,stopi,stopj,stopstate,dpenv);  
      *donej += (stopj - startj);   /* Now read it off into out */ 
      if( read_hidden_DnaProfileMat(mat,starti,startj,startstate,stopi,stopj,stopstate,out) == FALSE)    {  
        warn("In full dc, at %d:%d,%d:%d got a bad hidden explicit read off... ",starti,startj,stopi,stopj); 
        return FALSE;    
        }  
      return TRUE;   
      }  


/* In actual divide and conquor */ 
    if( do_dc_single_pass_DnaProfileMat(mat,starti,startj,startstate,stopi,stopj,stopstate,dpenv,(int)(*donej*100)/totalj) == FALSE) {  
      warn("In divide and conquor for DnaProfileMat, at bound %d:%d to %d:%d, unable to calculate midpoint. Problem!",starti,startj,stopi,stopj);    
      return FALSE;  
      }  


/* Ok... now we have to call on each side of the matrix */ 
/* We have to retrieve left hand side positions, as they will be vapped by the time we call LHS */ 
    lstarti= DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,0);     
    lstartj= DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,1);     
    lstate = DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,2);     


/* Call on right hand side: this lets us do the correct read off */ 
    if( full_dc_DnaProfileMat(mat,DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,3),DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,4),DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,5),stopi,stopj,stopstate,out,donej,totalj,dpenv) == FALSE)   {  
/* Warning already issued, simply chained back up to top */ 
      return FALSE;  
      }  
/* Call on left hand side */ 
    if( full_dc_DnaProfileMat(mat,starti,startj,startstate,lstarti,lstartj,lstate,out,donej,totalj,dpenv) == FALSE)  {  
/* Warning already issued, simply chained back up to top */ 
      return FALSE;  
      }  


    return TRUE;     
}    


/* Function:  do_dc_single_pass_DnaProfileMat(mat,starti,startj,startstate,stopi,stopj,stopstate,dpenv,perc_done)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [DnaProfileMat *]
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
boolean do_dc_single_pass_DnaProfileMat(DnaProfileMat * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv,int perc_done) 
{
    int halfj;   
    halfj = startj + ((stopj - startj)/2);   


    init_dc_DnaProfileMat(mat);  


    DnaProfileMat_DC_SHADOW_MATRIX(mat,starti,startj,startstate) = 0;    
    run_up_dc_DnaProfileMat(mat,starti,stopi,startj,halfj-1,dpenv,perc_done);    
    push_dc_at_merge_DnaProfileMat(mat,starti,stopi,halfj,&halfj,dpenv);     
    follow_on_dc_DnaProfileMat(mat,starti,stopi,halfj,stopj,dpenv,perc_done);    
    return TRUE; 
}    


/* Function:  push_dc_at_merge_DnaProfileMat(mat,starti,stopi,startj,stopj,dpenv)
 *
 * Descrip: No Description
 *
 * Arg:           mat [UNKN ] Undocumented argument [DnaProfileMat *]
 * Arg:        starti [UNKN ] Undocumented argument [int]
 * Arg:         stopi [UNKN ] Undocumented argument [int]
 * Arg:        startj [UNKN ] Undocumented argument [int]
 * Arg:         stopj [UNKN ] Undocumented argument [int *]
 * Arg:         dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 */
void push_dc_at_merge_DnaProfileMat(DnaProfileMat * mat,int starti,int stopi,int startj,int * stopj,DPEnvelope * dpenv) 
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
          DnaProfileMat_DC_SHADOW_MATRIX(mat,i,j,MATCH) = NEGI;  
          DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,0) = (-100);   
          DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,1) = (-100);   
          DnaProfileMat_DC_SHADOW_MATRIX(mat,i,j,UNMATCHED) = NEGI;  
          DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED,0) = (-100);   
          DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED,1) = (-100);   
          DnaProfileMat_DC_SHADOW_MATRIX(mat,i,j,INSERT) = NEGI;     
          DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,0) = (-100);  
          DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,1) = (-100);  
          DnaProfileMat_DC_SHADOW_MATRIX(mat,i,j,DELETE) = NEGI;     
          DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,0) = (-100);  
          DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,1) = (-100);  
          continue;  
          } /* end of Is not in envelope */ 


        /* For state MATCH, pushing when j - offj <= mergej */ 
        score = DnaProfileMat_DC_SHADOW_MATRIX(mat,i-1,j-1,MATCH) + 0;   
        if( j - 1 <= mergej) {  
          DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,0) = i-1;  
          DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,1) = j-1;  
          DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,2) = MATCH;    
          DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,3) = i;    
          DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,4) = j;    
          DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,5) = MATCH;    
          }  
        else {  
          for(k=0;k<7;k++)   
            DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,MATCH,k); 
          }  


        temp = DnaProfileMat_DC_SHADOW_MATRIX(mat,i-1,j-1,UNMATCHED) + 0;    
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,0) = i-1;    
            DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,1) = j-1;    
            DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,2) = UNMATCHED;  
            DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,3) = i;  
            DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,4) = j;  
            DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,5) = MATCH;  
            }  
          else   {  
            for(k=0;k<7;k++) 
              DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,UNMATCHED,k);   
            }  
          }  


        temp = DnaProfileMat_DC_SHADOW_MATRIX(mat,i-1,j-1,INSERT) + 0;   
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,0) = i-1;    
            DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,1) = j-1;    
            DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,2) = INSERT; 
            DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,3) = i;  
            DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,4) = j;  
            DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,5) = MATCH;  
            }  
          else   {  
            for(k=0;k<7;k++) 
              DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,INSERT,k);  
            }  
          }  


        temp = DnaProfileMat_DC_SHADOW_MATRIX(mat,i-1,j-1,DELETE) + 0;   
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,0) = i-1;    
            DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,1) = j-1;    
            DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,2) = DELETE; 
            DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,3) = i;  
            DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,4) = j;  
            DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,5) = MATCH;  
            }  
          else   {  
            for(k=0;k<7;k++) 
              DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,DELETE,k);  
            }  
          }  
        /* Add any movement independant score */ 
        score += mat->m->score[i][j];    
        DnaProfileMat_DC_SHADOW_MATRIX(mat,i,j,MATCH) = score;   
        /* Finished with state MATCH */ 


        /* For state UNMATCHED, pushing when j - offj <= mergej */ 
        score = DnaProfileMat_DC_SHADOW_MATRIX(mat,i-1,j-1,MATCH) + mat->open_unmatched;     
        if( j - 1 <= mergej) {  
          DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED,0) = i-1;  
          DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED,1) = j-1;  
          DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED,2) = MATCH;    
          DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED,3) = i;    
          DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED,4) = j;    
          DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED,5) = UNMATCHED;    
          }  
        else {  
          for(k=0;k<7;k++)   
            DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED,k) = DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,MATCH,k); 
          }  


        temp = DnaProfileMat_DC_SHADOW_MATRIX(mat,i-1,j-1,UNMATCHED) + mat->ext_unmatched;   
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED,0) = i-1;    
            DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED,1) = j-1;    
            DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED,2) = UNMATCHED;  
            DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED,3) = i;  
            DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED,4) = j;  
            DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED,5) = UNMATCHED;  
            }  
          else   {  
            for(k=0;k<7;k++) 
              DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED,k) = DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,UNMATCHED,k);   
            }  
          }  


        temp = DnaProfileMat_DC_SHADOW_MATRIX(mat,i-0,j-1,UNMATCHED) + (mat->ext_unmatched+mat->gap_unmatched);  
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED,0) = i-0;    
            DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED,1) = j-1;    
            DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED,2) = UNMATCHED;  
            DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED,3) = i;  
            DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED,4) = j;  
            DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED,5) = UNMATCHED;  
            }  
          else   {  
            for(k=0;k<7;k++) 
              DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED,k) = DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,UNMATCHED,k);   
            }  
          }  


        temp = DnaProfileMat_DC_SHADOW_MATRIX(mat,i-1,j-0,UNMATCHED) + (mat->ext_unmatched+mat->gap_unmatched);  
        if( temp > score)    {  
          score = temp;  


          if( j - 0 <= mergej)   {  
            DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED,0) = i-1;    
            DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED,1) = j-0;    
            DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED,2) = UNMATCHED;  
            DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED,3) = i;  
            DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED,4) = j;  
            DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED,5) = UNMATCHED;  
            }  
          else   {  
            for(k=0;k<7;k++) 
              DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED,k) = DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 0,UNMATCHED,k);   
            }  
          }  
        /* Add any movement independant score */ 
        DnaProfileMat_DC_SHADOW_MATRIX(mat,i,j,UNMATCHED) = score;   
        /* Finished with state UNMATCHED */ 


        /* For state INSERT, pushing when j - offj <= mergej */ 
        score = DnaProfileMat_DC_SHADOW_MATRIX(mat,i-0,j-1,MATCH) + (mat->q->col[i]->trans[DnaProfile_M2I]+mat->t->col[j]->trans[DnaProfile_M2D]);   
        if( j - 1 <= mergej) {  
          DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,0) = i-0; 
          DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,1) = j-1; 
          DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,2) = MATCH;   
          DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,3) = i;   
          DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,4) = j;   
          DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,5) = INSERT;  
          }  
        else {  
          for(k=0;k<7;k++)   
            DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,k) = DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,MATCH,k);    
          }  


        temp = DnaProfileMat_DC_SHADOW_MATRIX(mat,i-0,j-1,INSERT) + (mat->q->col[i]->trans[DnaProfile_I2I]+mat->t->col[j]->trans[DnaProfile_D2D]);   
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,0) = i-0;   
            DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,1) = j-1;   
            DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,2) = INSERT;    
            DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,3) = i; 
            DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,4) = j; 
            DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,5) = INSERT;    
            }  
          else   {  
            for(k=0;k<7;k++) 
              DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,k) = DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,INSERT,k); 
            }  
          }  
        /* Add any movement independant score */ 
        DnaProfileMat_DC_SHADOW_MATRIX(mat,i,j,INSERT) = score;  
        /* Finished with state INSERT */ 


        /* For state DELETE, pushing when j - offj <= mergej */ 
        score = DnaProfileMat_DC_SHADOW_MATRIX(mat,i-1,j-0,MATCH) + (mat->t->col[j]->trans[DnaProfile_M2I]+mat->q->col[i]->trans[DnaProfile_M2D]);   
        if( j - 0 <= mergej) {  
          DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,0) = i-1; 
          DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,1) = j-0; 
          DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,2) = MATCH;   
          DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,3) = i;   
          DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,4) = j;   
          DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,5) = DELETE;  
          }  
        else {  
          for(k=0;k<7;k++)   
            DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,k) = DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 0,MATCH,k);    
          }  


        temp = DnaProfileMat_DC_SHADOW_MATRIX(mat,i-1,j-0,INSERT) + (mat->t->col[j]->trans[DnaProfile_I2I]+mat->q->col[i]->trans[DnaProfile_D2D]);   
        if( temp > score)    {  
          score = temp;  


          if( j - 0 <= mergej)   {  
            DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,0) = i-1;   
            DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,1) = j-0;   
            DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,2) = INSERT;    
            DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,3) = i; 
            DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,4) = j; 
            DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,5) = DELETE;    
            }  
          else   {  
            for(k=0;k<7;k++) 
              DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,k) = DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 0,INSERT,k); 
            }  
          }  
        /* Add any movement independant score */ 
        DnaProfileMat_DC_SHADOW_MATRIX(mat,i,j,DELETE) = score;  
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


/* Function:  follow_on_dc_DnaProfileMat(mat,starti,stopi,startj,stopj,dpenv,perc_done)
 *
 * Descrip: No Description
 *
 * Arg:              mat [UNKN ] Undocumented argument [DnaProfileMat *]
 * Arg:           starti [UNKN ] Undocumented argument [int]
 * Arg:            stopi [UNKN ] Undocumented argument [int]
 * Arg:           startj [UNKN ] Undocumented argument [int]
 * Arg:            stopj [UNKN ] Undocumented argument [int]
 * Arg:            dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:        perc_done [UNKN ] Undocumented argument [int]
 *
 */
void follow_on_dc_DnaProfileMat(DnaProfileMat * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done) 
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
          DnaProfileMat_DC_SHADOW_MATRIX(mat,i,j,MATCH) = NEGI;  
          DnaProfileMat_DC_SHADOW_MATRIX(mat,i,j,UNMATCHED) = NEGI;  
          DnaProfileMat_DC_SHADOW_MATRIX(mat,i,j,INSERT) = NEGI;     
          DnaProfileMat_DC_SHADOW_MATRIX(mat,i,j,DELETE) = NEGI;     
          continue;  
          } /* end of Is not in envelope */ 
        if( num % 1000 == 0 )    
          log_full_error(REPORT,0,"[%d%%%% done]After  mid-j %5d Cells done %d%%%%",perc_done,startj,(num*100)/total);   


        /* For state MATCH */ 
        /* setting first movement to score */ 
        score = DnaProfileMat_DC_SHADOW_MATRIX(mat,i-1,j-1,MATCH) + 0;   
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,MATCH,k);   
        /* From state UNMATCHED to state MATCH */ 
        temp = DnaProfileMat_DC_SHADOW_MATRIX(mat,i-1,j-1,UNMATCHED) + 0;    
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,UNMATCHED,k); 
          }  
        /* From state INSERT to state MATCH */ 
        temp = DnaProfileMat_DC_SHADOW_MATRIX(mat,i-1,j-1,INSERT) + 0;   
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,INSERT,k);    
          }  
        /* From state DELETE to state MATCH */ 
        temp = DnaProfileMat_DC_SHADOW_MATRIX(mat,i-1,j-1,DELETE) + 0;   
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,DELETE,k);    
          }  


        /* Ok - finished max calculation for MATCH */ 
        /* Add any movement independant score and put away */ 
         score += mat->m->score[i][j];   
         DnaProfileMat_DC_SHADOW_MATRIX(mat,i,j,MATCH) = score;  
        for(k=0;k<7;k++) 
          DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = localshadow[k];   
        /* Now figure out if any specials need this score */ 
        /* Finished calculating state MATCH */ 


        /* For state UNMATCHED */ 
        /* setting first movement to score */ 
        score = DnaProfileMat_DC_SHADOW_MATRIX(mat,i-1,j-1,MATCH) + mat->open_unmatched;     
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,MATCH,k);   
        /* From state UNMATCHED to state UNMATCHED */ 
        temp = DnaProfileMat_DC_SHADOW_MATRIX(mat,i-1,j-1,UNMATCHED) + mat->ext_unmatched;   
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,UNMATCHED,k); 
          }  
        /* From state UNMATCHED to state UNMATCHED */ 
        temp = DnaProfileMat_DC_SHADOW_MATRIX(mat,i-0,j-1,UNMATCHED) + (mat->ext_unmatched+mat->gap_unmatched);  
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,UNMATCHED,k); 
          }  
        /* From state UNMATCHED to state UNMATCHED */ 
        temp = DnaProfileMat_DC_SHADOW_MATRIX(mat,i-1,j-0,UNMATCHED) + (mat->ext_unmatched+mat->gap_unmatched);  
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 0,UNMATCHED,k); 
          }  


        /* Ok - finished max calculation for UNMATCHED */ 
        /* Add any movement independant score and put away */ 
         DnaProfileMat_DC_SHADOW_MATRIX(mat,i,j,UNMATCHED) = score;  
        for(k=0;k<7;k++) 
          DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED,k) = localshadow[k];   
        /* Now figure out if any specials need this score */ 
        /* Finished calculating state UNMATCHED */ 


        /* For state INSERT */ 
        /* setting first movement to score */ 
        score = DnaProfileMat_DC_SHADOW_MATRIX(mat,i-0,j-1,MATCH) + (mat->q->col[i]->trans[DnaProfile_M2I]+mat->t->col[j]->trans[DnaProfile_M2D]);   
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,MATCH,k);   
        /* From state INSERT to state INSERT */ 
        temp = DnaProfileMat_DC_SHADOW_MATRIX(mat,i-0,j-1,INSERT) + (mat->q->col[i]->trans[DnaProfile_I2I]+mat->t->col[j]->trans[DnaProfile_D2D]);   
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,INSERT,k);    
          }  


        /* Ok - finished max calculation for INSERT */ 
        /* Add any movement independant score and put away */ 
         DnaProfileMat_DC_SHADOW_MATRIX(mat,i,j,INSERT) = score; 
        for(k=0;k<7;k++) 
          DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,k) = localshadow[k];  
        /* Now figure out if any specials need this score */ 
        /* Finished calculating state INSERT */ 


        /* For state DELETE */ 
        /* setting first movement to score */ 
        score = DnaProfileMat_DC_SHADOW_MATRIX(mat,i-1,j-0,MATCH) + (mat->t->col[j]->trans[DnaProfile_M2I]+mat->q->col[i]->trans[DnaProfile_M2D]);   
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 0,MATCH,k);   
        /* From state INSERT to state DELETE */ 
        temp = DnaProfileMat_DC_SHADOW_MATRIX(mat,i-1,j-0,INSERT) + (mat->t->col[j]->trans[DnaProfile_I2I]+mat->q->col[i]->trans[DnaProfile_D2D]);   
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 0,INSERT,k);    
          }  


        /* Ok - finished max calculation for DELETE */ 
        /* Add any movement independant score and put away */ 
         DnaProfileMat_DC_SHADOW_MATRIX(mat,i,j,DELETE) = score; 
        for(k=0;k<7;k++) 
          DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,k) = localshadow[k];  
        /* Now figure out if any specials need this score */ 
        /* Finished calculating state DELETE */ 
        } /* end of this is strip */ 
      } /* end of for each valid j column */ 


/* Function:  run_up_dc_DnaProfileMat(mat,starti,stopi,startj,stopj,dpenv,perc_done)
 *
 * Descrip: No Description
 *
 * Arg:              mat [UNKN ] Undocumented argument [DnaProfileMat *]
 * Arg:           starti [UNKN ] Undocumented argument [int]
 * Arg:            stopi [UNKN ] Undocumented argument [int]
 * Arg:           startj [UNKN ] Undocumented argument [int]
 * Arg:            stopj [UNKN ] Undocumented argument [int]
 * Arg:            dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:        perc_done [UNKN ] Undocumented argument [int]
 *
 */
}    
void run_up_dc_DnaProfileMat(DnaProfileMat * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done) 
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
          DnaProfileMat_DC_SHADOW_MATRIX(mat,i,j,MATCH) = NEGI;  
          DnaProfileMat_DC_SHADOW_MATRIX(mat,i,j,UNMATCHED) = NEGI;  
          DnaProfileMat_DC_SHADOW_MATRIX(mat,i,j,INSERT) = NEGI;     
          DnaProfileMat_DC_SHADOW_MATRIX(mat,i,j,DELETE) = NEGI;     
          continue;  
          } /* end of Is not in envelope */ 
        if( num % 1000 == 0 )    
          log_full_error(REPORT,0,"[%d%%%% done]Before mid-j %5d Cells done %d%%%%",perc_done,stopj,(num*100)/total);    


        /* For state MATCH */ 
        /* setting first movement to score */ 
        score = DnaProfileMat_DC_SHADOW_MATRIX(mat,i-1,j-1,MATCH) + 0;   
        /* From state UNMATCHED to state MATCH */ 
        temp = DnaProfileMat_DC_SHADOW_MATRIX(mat,i-1,j-1,UNMATCHED) + 0;    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state MATCH */ 
        temp = DnaProfileMat_DC_SHADOW_MATRIX(mat,i-1,j-1,INSERT) + 0;   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state MATCH */ 
        temp = DnaProfileMat_DC_SHADOW_MATRIX(mat,i-1,j-1,DELETE) + 0;   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for MATCH */ 
        /* Add any movement independant score and put away */ 
         score += mat->m->score[i][j];   
         DnaProfileMat_DC_SHADOW_MATRIX(mat,i,j,MATCH) = score;  
        /* Finished calculating state MATCH */ 


        /* For state UNMATCHED */ 
        /* setting first movement to score */ 
        score = DnaProfileMat_DC_SHADOW_MATRIX(mat,i-1,j-1,MATCH) + mat->open_unmatched;     
        /* From state UNMATCHED to state UNMATCHED */ 
        temp = DnaProfileMat_DC_SHADOW_MATRIX(mat,i-1,j-1,UNMATCHED) + mat->ext_unmatched;   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state UNMATCHED to state UNMATCHED */ 
        temp = DnaProfileMat_DC_SHADOW_MATRIX(mat,i-0,j-1,UNMATCHED) + (mat->ext_unmatched+mat->gap_unmatched);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state UNMATCHED to state UNMATCHED */ 
        temp = DnaProfileMat_DC_SHADOW_MATRIX(mat,i-1,j-0,UNMATCHED) + (mat->ext_unmatched+mat->gap_unmatched);  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for UNMATCHED */ 
        /* Add any movement independant score and put away */ 
         DnaProfileMat_DC_SHADOW_MATRIX(mat,i,j,UNMATCHED) = score;  
        /* Finished calculating state UNMATCHED */ 


        /* For state INSERT */ 
        /* setting first movement to score */ 
        score = DnaProfileMat_DC_SHADOW_MATRIX(mat,i-0,j-1,MATCH) + (mat->q->col[i]->trans[DnaProfile_M2I]+mat->t->col[j]->trans[DnaProfile_M2D]);   
        /* From state INSERT to state INSERT */ 
        temp = DnaProfileMat_DC_SHADOW_MATRIX(mat,i-0,j-1,INSERT) + (mat->q->col[i]->trans[DnaProfile_I2I]+mat->t->col[j]->trans[DnaProfile_D2D]);   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for INSERT */ 
        /* Add any movement independant score and put away */ 
         DnaProfileMat_DC_SHADOW_MATRIX(mat,i,j,INSERT) = score; 
        /* Finished calculating state INSERT */ 


        /* For state DELETE */ 
        /* setting first movement to score */ 
        score = DnaProfileMat_DC_SHADOW_MATRIX(mat,i-1,j-0,MATCH) + (mat->t->col[j]->trans[DnaProfile_M2I]+mat->q->col[i]->trans[DnaProfile_M2D]);   
        /* From state INSERT to state DELETE */ 
        temp = DnaProfileMat_DC_SHADOW_MATRIX(mat,i-1,j-0,INSERT) + (mat->t->col[j]->trans[DnaProfile_I2I]+mat->q->col[i]->trans[DnaProfile_D2D]);   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for DELETE */ 
        /* Add any movement independant score and put away */ 
         DnaProfileMat_DC_SHADOW_MATRIX(mat,i,j,DELETE) = score; 
        /* Finished calculating state DELETE */ 
        } /* end of this is strip */ 
      } /* end of for each valid j column */ 


/* Function:  init_dc_DnaProfileMat(mat)
 *
 * Descrip: No Description
 *
 * Arg:        mat [UNKN ] Undocumented argument [DnaProfileMat *]
 *
 */
}    
void init_dc_DnaProfileMat(DnaProfileMat * mat) 
{
    register int i;  
    register int j;  
    register int k;  


    for(j=0;j<3;j++) {  
      for(i=(-1);i<mat->q->len;i++)  {  
        DnaProfileMat_DC_SHADOW_MATRIX(mat,i,j,MATCH) = NEGI;    
        DnaProfileMat_DC_SHADOW_MATRIX(mat,i,j,UNMATCHED) = NEGI;    
        DnaProfileMat_DC_SHADOW_MATRIX(mat,i,j,INSERT) = NEGI;   
        DnaProfileMat_DC_SHADOW_MATRIX(mat,i,j,DELETE) = NEGI;   
        for(k=0;k<7;k++) {  
          DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = (-1); 
          DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED,k) = (-1); 
          DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,k) = (-1);    
          DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,k) = (-1);    
          }  
        }  
      }  


    return;  
}    


/* Function:  start_end_find_end_DnaProfileMat(mat,endj)
 *
 * Descrip:    First function used to find end of the best path in the special state !end
 *
 *
 * Arg:         mat [UNKN ] Matrix in small mode [DnaProfileMat *]
 * Arg:        endj [WRITE] position of end in j (meaningless in i) [int *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int start_end_find_end_DnaProfileMat(DnaProfileMat * mat,int * endj) 
{
    register int j;  
    register int max;    
    register int maxj;   


    max = DnaProfileMat_DC_SHADOW_SPECIAL(mat,0,mat->t->len-1,END);  
    maxj = mat->t->len-1;    
    for(j= mat->t->len-2 ;j >= 0 ;j--)   {  
      if( DnaProfileMat_DC_SHADOW_SPECIAL(mat,0,j,END) > max )   {  
        max = DnaProfileMat_DC_SHADOW_SPECIAL(mat,0,j,END);  
        maxj = j;    
        }  
      }  


    if( endj != NULL)    
      *endj = maxj;  


    return max;  
}    


/* Function:  dc_optimised_start_end_calc_DnaProfileMat(*mat,dpenv)
 *
 * Descrip:    Calculates special strip, leaving start/end/score points in shadow matrix
 *             Works off specially laid out memory from steve searle
 *
 *
 * Arg:         *mat [UNKN ] Undocumented argument [DnaProfileMat]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean dc_optimised_start_end_calc_DnaProfileMat(DnaProfileMat *mat,DPEnvelope * dpenv) 
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
    leni = mat->q->len;  
    lenj = mat->t->len;  
    total = leni * lenj; 


    score_pointers = (int *) calloc (1 * (leni + 1) * 4,sizeof(int));    
    shadow_pointers = (int *) calloc (1 * (leni + 1) * 4 * 8,sizeof(int));   


    for(j=0;j<lenj;j++)  { /*for each j strip*/ 
      for(i=0;i<leni;i++)    { /*for each i position in strip*/ 
        num++;   
        if( dpenv != NULL && is_in_DPEnvelope(dpenv,i,j) == FALSE )  { /*Is not in envelope*/ 
          DnaProfileMat_DC_OPT_SHADOW_MATRIX(mat,i,j,MATCH) = NEGI;  
          DnaProfileMat_DC_OPT_SHADOW_MATRIX(mat,i,j,UNMATCHED) = NEGI;  
          DnaProfileMat_DC_OPT_SHADOW_MATRIX(mat,i,j,INSERT) = NEGI;     
          DnaProfileMat_DC_OPT_SHADOW_MATRIX(mat,i,j,DELETE) = NEGI;     
          continue;  
          } /* end of Is not in envelope */ 
        if( num%1000 == 0)   
          log_full_error(REPORT,0,"%6d Cells done [%2d%%%%]",num,num*100/total); 




        /* For state MATCH */ 
        /* setting first movement to score */ 
        score = DnaProfileMat_DC_OPT_SHADOW_MATRIX(mat,i-1,j-1,MATCH) + 0 + (mat->m->score[i][j]);   
        /* assign local shadown pointer */ 
        localsp = &(DnaProfileMat_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 1,MATCH,0)); 
        /* From state UNMATCHED to state MATCH */ 
        temp = DnaProfileMat_DC_OPT_SHADOW_MATRIX(mat,i-1,j-1,UNMATCHED) + 0 +(mat->m->score[i][j]);     
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(DnaProfileMat_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 1,UNMATCHED,0));   
          }  
        /* From state INSERT to state MATCH */ 
        temp = DnaProfileMat_DC_OPT_SHADOW_MATRIX(mat,i-1,j-1,INSERT) + 0 +(mat->m->score[i][j]);    
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(DnaProfileMat_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 1,INSERT,0));  
          }  
        /* From state DELETE to state MATCH */ 
        temp = DnaProfileMat_DC_OPT_SHADOW_MATRIX(mat,i-1,j-1,DELETE) + 0 +(mat->m->score[i][j]);    
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(DnaProfileMat_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 1,DELETE,0));  
          }  
        /* From state START to state MATCH */ 
        temp = DnaProfileMat_DC_OPT_SHADOW_SPECIAL(mat,i-1,j-1,START) + 0 + (mat->m->score[i][j]);   
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
         DnaProfileMat_DC_OPT_SHADOW_MATRIX(mat,i,j,MATCH) = score;  
        for(k=0;k<7;k++) 
          DnaProfileMat_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = localsp[k];   
        /* Now figure out if any specials need this score */ 


        /* state MATCH is a source for special END */ 
        temp = score + (0) + (0) ;   
        if( temp > DnaProfileMat_DC_OPT_SHADOW_SPECIAL(mat,i,j,END) )    {  
          DnaProfileMat_DC_OPT_SHADOW_SPECIAL(mat,i,j,END) = temp;   
          /* Have to push only bottem half of system here */ 
          for(k=0;k<3;k++)   
            DnaProfileMat_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,k) = DnaProfileMat_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,MATCH,k);  
          DnaProfileMat_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,6) = DnaProfileMat_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,MATCH,6);    
          DnaProfileMat_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,3) = i; 
          DnaProfileMat_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,4) = j; 
          DnaProfileMat_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,5) = MATCH; 
          }  




        /* Finished calculating state MATCH */ 


        /* For state UNMATCHED */ 
        /* setting first movement to score */ 
        score = DnaProfileMat_DC_OPT_SHADOW_MATRIX(mat,i-1,j-1,MATCH) + mat->open_unmatched + (0);   
        /* assign local shadown pointer */ 
        localsp = &(DnaProfileMat_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 1,MATCH,0)); 
        /* From state UNMATCHED to state UNMATCHED */ 
        temp = DnaProfileMat_DC_OPT_SHADOW_MATRIX(mat,i-1,j-1,UNMATCHED) + mat->ext_unmatched +(0);  
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(DnaProfileMat_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 1,UNMATCHED,0));   
          }  
        /* From state UNMATCHED to state UNMATCHED */ 
        temp = DnaProfileMat_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,UNMATCHED) + (mat->ext_unmatched+mat->gap_unmatched) +(0);     
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(DnaProfileMat_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,UNMATCHED,0));   
          }  
        /* From state UNMATCHED to state UNMATCHED */ 
        temp = DnaProfileMat_DC_OPT_SHADOW_MATRIX(mat,i-1,j-0,UNMATCHED) + (mat->ext_unmatched+mat->gap_unmatched) +(0);     
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(DnaProfileMat_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 0,UNMATCHED,0));   
          }  


        /* Ok - finished max calculation for UNMATCHED */ 
        /* Add any movement independant score and put away */ 
        /* Actually, already done inside scores */ 
         DnaProfileMat_DC_OPT_SHADOW_MATRIX(mat,i,j,UNMATCHED) = score;  
        for(k=0;k<7;k++) 
          DnaProfileMat_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED,k) = localsp[k];   
        /* Now figure out if any specials need this score */ 


        /* Finished calculating state UNMATCHED */ 


        /* For state INSERT */ 
        /* setting first movement to score */ 
        score = DnaProfileMat_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,MATCH) + (mat->q->col[i]->trans[DnaProfile_M2I]+mat->t->col[j]->trans[DnaProfile_M2D]) + (0);     
        /* assign local shadown pointer */ 
        localsp = &(DnaProfileMat_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,MATCH,0)); 
        /* From state INSERT to state INSERT */ 
        temp = DnaProfileMat_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,INSERT) + (mat->q->col[i]->trans[DnaProfile_I2I]+mat->t->col[j]->trans[DnaProfile_D2D]) +(0);  
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(DnaProfileMat_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,INSERT,0));  
          }  


        /* Ok - finished max calculation for INSERT */ 
        /* Add any movement independant score and put away */ 
        /* Actually, already done inside scores */ 
         DnaProfileMat_DC_OPT_SHADOW_MATRIX(mat,i,j,INSERT) = score; 
        for(k=0;k<7;k++) 
          DnaProfileMat_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,INSERT,k) = localsp[k];  
        /* Now figure out if any specials need this score */ 


        /* Finished calculating state INSERT */ 


        /* For state DELETE */ 
        /* setting first movement to score */ 
        score = DnaProfileMat_DC_OPT_SHADOW_MATRIX(mat,i-1,j-0,MATCH) + (mat->t->col[j]->trans[DnaProfile_M2I]+mat->q->col[i]->trans[DnaProfile_M2D]) + (0);     
        /* assign local shadown pointer */ 
        localsp = &(DnaProfileMat_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 0,MATCH,0)); 
        /* From state INSERT to state DELETE */ 
        temp = DnaProfileMat_DC_OPT_SHADOW_MATRIX(mat,i-1,j-0,INSERT) + (mat->t->col[j]->trans[DnaProfile_I2I]+mat->q->col[i]->trans[DnaProfile_D2D]) +(0);  
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(DnaProfileMat_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 0,INSERT,0));  
          }  


        /* Ok - finished max calculation for DELETE */ 
        /* Add any movement independant score and put away */ 
        /* Actually, already done inside scores */ 
         DnaProfileMat_DC_OPT_SHADOW_MATRIX(mat,i,j,DELETE) = score; 
        for(k=0;k<7;k++) 
          DnaProfileMat_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,DELETE,k) = localsp[k];  
        /* Now figure out if any specials need this score */ 


        /* Finished calculating state DELETE */ 


        } /* end of for each i position in strip */ 
      } /* end of for each j strip */ 
    free(score_pointers);    
    free(shadow_pointers);   
    return TRUE;     
}    


/* Function:  init_start_end_linear_DnaProfileMat(mat)
 *
 * Descrip: No Description
 *
 * Arg:        mat [UNKN ] Undocumented argument [DnaProfileMat *]
 *
 */
void init_start_end_linear_DnaProfileMat(DnaProfileMat * mat) 
{
    register int i;  
    register int j;  
    for(j=0;j<3;j++) {  
      for(i=(-1);i<mat->q->len;i++)  {  
        DnaProfileMat_DC_SHADOW_MATRIX(mat,i,j,MATCH) = NEGI;    
        DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,0) = (-1);   
        DnaProfileMat_DC_SHADOW_MATRIX(mat,i,j,UNMATCHED) = NEGI;    
        DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,UNMATCHED,0) = (-1);   
        DnaProfileMat_DC_SHADOW_MATRIX(mat,i,j,INSERT) = NEGI;   
        DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,0) = (-1);  
        DnaProfileMat_DC_SHADOW_MATRIX(mat,i,j,DELETE) = NEGI;   
        DnaProfileMat_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,0) = (-1);  
        }  
      }  


    for(j=(-1);j<mat->t->len;j++)    {  
      DnaProfileMat_DC_SHADOW_SPECIAL(mat,0,j,START) = 0;    
      DnaProfileMat_DC_SHADOW_SPECIAL_SP(mat,0,j,START,0) = j;   
      DnaProfileMat_DC_SHADOW_SPECIAL(mat,0,j,END) = NEGI;   
      DnaProfileMat_DC_SHADOW_SPECIAL_SP(mat,0,j,END,0) = (-1);  
      }  


    return;  
}    


/* Function:  convert_PackAln_to_AlnBlock_DnaProfileMat(pal)
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
AlnBlock * convert_PackAln_to_AlnBlock_DnaProfileMat(PackAln * pal) 
{
    AlnConvertSet * acs; 
    AlnBlock * alb;  


    acs = AlnConvertSet_DnaProfileMat(); 
    alb = AlnBlock_from_PackAln(acs,pal);    
    free_AlnConvertSet(acs); 
    return alb;  
}    


 static char * query_label[] = { "MATCH","SEQ_UNMATCHED","PAUSED_UNMATCHED","PAUSE","INSERT","END" };    
/* Function:  AlnConvertSet_DnaProfileMat(void)
 *
 * Descrip: No Description
 *
 *
 * Return [UNKN ]  Undocumented return value [AlnConvertSet *]
 *
 */
 static char * target_label[] = { "MATCH","SEQ_UNMATCHED","PAUSED_UNMATCHED","INSERT","PAUSE","END" };   
AlnConvertSet * AlnConvertSet_DnaProfileMat(void) 
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
    acu->state1 = UNMATCHED; 
    acu->state2 = MATCH;     
    acu->offi = 1;   
    acu->offj = 1;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[0];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = INSERT;    
    acu->state2 = MATCH;     
    acu->offi = 1;   
    acu->offj = 1;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[0];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = DELETE;    
    acu->state2 = MATCH;     
    acu->offi = 1;   
    acu->offj = 1;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[0];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = START + 4; 
    acu->is_from_special = TRUE; 
    acu->state2 = MATCH;     
    acu->offi = (-1);    
    acu->offj = 1;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[0];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = MATCH; 
    acu->state2 = UNMATCHED;     
    acu->offi = 1;   
    acu->offj = 1;   
    acu->label1 = query_label[1];    
    acu->label2 = target_label[1];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = UNMATCHED; 
    acu->state2 = UNMATCHED;     
    acu->offi = 1;   
    acu->offj = 1;   
    acu->label1 = query_label[1];    
    acu->label2 = target_label[1];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = UNMATCHED; 
    acu->state2 = UNMATCHED;     
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[2];    
    acu->label2 = target_label[1];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = UNMATCHED; 
    acu->state2 = UNMATCHED;     
    acu->offi = 1;   
    acu->offj = 0;   
    acu->label1 = query_label[1];    
    acu->label2 = target_label[2];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = MATCH; 
    acu->state2 = INSERT;    
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[3];    
    acu->label2 = target_label[3];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = INSERT;    
    acu->state2 = INSERT;    
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[3];    
    acu->label2 = target_label[3];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = MATCH; 
    acu->state2 = DELETE;    
    acu->offi = 1;   
    acu->offj = 0;   
    acu->label1 = query_label[4];    
    acu->label2 = target_label[4];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = INSERT;    
    acu->state2 = DELETE;    
    acu->offi = 1;   
    acu->offj = 0;   
    acu->label1 = query_label[4];    
    acu->label2 = target_label[4];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = MATCH; 
    acu->state2 = END + 4;   
    acu->offi = (-1);    
    acu->offj = 0;   
    acu->label1 = query_label[5];    
    acu->label2 = target_label[5];   
    return out;  
}    


/* Function:  PackAln_read_Expl_DnaProfileMat(mat)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [DnaProfileMat *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * PackAln_read_Expl_DnaProfileMat(DnaProfileMat * mat) 
{
    DnaProfileMat_access_func_holder holder;     


    holder.access_main    = DnaProfileMat_explicit_access_main;  
    holder.access_special = DnaProfileMat_explicit_access_special;   
    return PackAln_read_generic_DnaProfileMat(mat,holder);   
}    


/* Function:  DnaProfileMat_explicit_access_main(mat,i,j,state)
 *
 * Descrip: No Description
 *
 * Arg:          mat [UNKN ] Undocumented argument [DnaProfileMat *]
 * Arg:            i [UNKN ] Undocumented argument [int]
 * Arg:            j [UNKN ] Undocumented argument [int]
 * Arg:        state [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int DnaProfileMat_explicit_access_main(DnaProfileMat * mat,int i,int j,int state) 
{
    return DnaProfileMat_EXPL_MATRIX(mat,i,j,state); 
}    


/* Function:  DnaProfileMat_explicit_access_special(mat,i,j,state)
 *
 * Descrip: No Description
 *
 * Arg:          mat [UNKN ] Undocumented argument [DnaProfileMat *]
 * Arg:            i [UNKN ] Undocumented argument [int]
 * Arg:            j [UNKN ] Undocumented argument [int]
 * Arg:        state [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int DnaProfileMat_explicit_access_special(DnaProfileMat * mat,int i,int j,int state) 
{
    return DnaProfileMat_EXPL_SPECIAL(mat,i,j,state);    
}    


/* Function:  PackAln_read_generic_DnaProfileMat(mat,h)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [DnaProfileMat *]
 * Arg:          h [UNKN ] Undocumented argument [DnaProfileMat_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * PackAln_read_generic_DnaProfileMat(DnaProfileMat * mat,DnaProfileMat_access_func_holder h) 
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


    out->score =  find_end_DnaProfileMat(mat,&i,&j,&state,&isspecial,h); 


    /* Add final end transition (at the moment we have not got the score! */ 
    if( (pau= PackAlnUnit_alloc()) == NULL  || add_PackAln(out,pau) == FALSE )   {  
      warn("Failed the first PackAlnUnit alloc, %d length of Alignment in DnaProfileMat_basic_read, returning a mess.(Sorry!)",out->len);    
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
        max_calc_special_DnaProfileMat(mat,i,j,state,isspecial,&i,&j,&state,&isspecial,&cellscore,h);    
      else   
        max_calc_DnaProfileMat(mat,i,j,state,isspecial,&i,&j,&state,&isspecial,&cellscore,h);    
      if(i == DnaProfileMat_READ_OFF_ERROR || j == DnaProfileMat_READ_OFF_ERROR || state == DnaProfileMat_READ_OFF_ERROR )   {  
        warn("Problem - hit bad read off system, exiting now");  
        break;   
        }  
      if( (pau= PackAlnUnit_alloc()) == NULL  || add_PackAln(out,pau) == FALSE ) {  
        warn("Failed a PackAlnUnit alloc, %d length of Alignment in DnaProfileMat_basic_read, returning partial alignment",out->len);    
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


/* Function:  find_end_DnaProfileMat(mat,ri,rj,state,isspecial,h)
 *
 * Descrip: No Description
 *
 * Arg:              mat [UNKN ] Undocumented argument [DnaProfileMat *]
 * Arg:               ri [UNKN ] Undocumented argument [int *]
 * Arg:               rj [UNKN ] Undocumented argument [int *]
 * Arg:            state [UNKN ] Undocumented argument [int *]
 * Arg:        isspecial [UNKN ] Undocumented argument [boolean *]
 * Arg:                h [UNKN ] Undocumented argument [DnaProfileMat_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int find_end_DnaProfileMat(DnaProfileMat * mat,int * ri,int * rj,int * state,boolean * isspecial,DnaProfileMat_access_func_holder h) 
{
    int j;   
    int max; 
    int maxj;    
    int temp;    


    max = (*h.access_special)(mat,0,mat->t->len-1,END);  
    maxj = mat->t->len-1;    
    for(j= mat->t->len-2 ;j >= 0 ;j--)   {  
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


/* Function:  DnaProfileMat_debug_show_matrix(mat,starti,stopi,startj,stopj,ofp)
 *
 * Descrip: No Description
 *
 * Arg:           mat [UNKN ] Undocumented argument [DnaProfileMat *]
 * Arg:        starti [UNKN ] Undocumented argument [int]
 * Arg:         stopi [UNKN ] Undocumented argument [int]
 * Arg:        startj [UNKN ] Undocumented argument [int]
 * Arg:         stopj [UNKN ] Undocumented argument [int]
 * Arg:           ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
void DnaProfileMat_debug_show_matrix(DnaProfileMat * mat,int starti,int stopi,int startj,int stopj,FILE * ofp) 
{
    register int i;  
    register int j;  


    for(i=starti;i<stopi && i < mat->q->len;i++) {  
      for(j=startj;j<stopj && j < mat->t->len;j++)   {  
        fprintf(ofp,"Cell [%d - %d]\n",i,j);     
        fprintf(ofp,"State MATCH %d\n",DnaProfileMat_EXPL_MATRIX(mat,i,j,MATCH));    
        fprintf(ofp,"State UNMATCHED %d\n",DnaProfileMat_EXPL_MATRIX(mat,i,j,UNMATCHED));    
        fprintf(ofp,"State INSERT %d\n",DnaProfileMat_EXPL_MATRIX(mat,i,j,INSERT));  
        fprintf(ofp,"State DELETE %d\n",DnaProfileMat_EXPL_MATRIX(mat,i,j,DELETE));  
        fprintf(ofp,"\n\n"); 
        }  
      }  


}    


/* Function:  max_calc_DnaProfileMat(mat,i,j,state,isspecial,reti,retj,retstate,retspecial,cellscore,h)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [DnaProfileMat *]
 * Arg:                 i [UNKN ] Undocumented argument [int]
 * Arg:                 j [UNKN ] Undocumented argument [int]
 * Arg:             state [UNKN ] Undocumented argument [int]
 * Arg:         isspecial [UNKN ] Undocumented argument [boolean]
 * Arg:              reti [UNKN ] Undocumented argument [int *]
 * Arg:              retj [UNKN ] Undocumented argument [int *]
 * Arg:          retstate [UNKN ] Undocumented argument [int *]
 * Arg:        retspecial [UNKN ] Undocumented argument [boolean *]
 * Arg:         cellscore [UNKN ] Undocumented argument [int *]
 * Arg:                 h [UNKN ] Undocumented argument [DnaProfileMat_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int max_calc_DnaProfileMat(DnaProfileMat * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,DnaProfileMat_access_func_holder h) 
{
    register int temp;   
    register int cscore; 


    *reti = (*retj) = (*retstate) = DnaProfileMat_READ_OFF_ERROR;    


    if( i < 0 || j < 0 || i > mat->q->len || j > mat->t->len)    {  
      warn("In DnaProfileMat matrix special read off - out of bounds on matrix [i,j is %d,%d state %d in standard matrix]",i,j,state);   
      return -1;     
      }  


    /* Then you have to select the correct switch statement to figure out the readoff      */ 
    /* Somewhat odd - reverse the order of calculation and return as soon as it is correct */ 
    cscore = (*h.access_main)(mat,i,j,state);    
    switch(state)    { /*Switch state */ 
      case MATCH :   
        temp = cscore - (0) -  (mat->m->score[i][j]);    
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
        temp = cscore - (0) -  (mat->m->score[i][j]);    
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
        temp = cscore - (0) -  (mat->m->score[i][j]);    
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
        temp = cscore - (0) -  (mat->m->score[i][j]);    
        if( temp == (*h.access_main)(mat,i - 1,j - 1,UNMATCHED) )    {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = UNMATCHED; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-1,UNMATCHED);   
            }  
          return (*h.access_main)(mat,i - 1,j - 1,UNMATCHED);    
          }  
        temp = cscore - (0) -  (mat->m->score[i][j]);    
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
        warn("Major problem (!) - in DnaProfileMat read off, position %d,%d state %d no source found!",i,j,state);   
        return (-1); 
      case UNMATCHED :   
        temp = cscore - ((mat->ext_unmatched+mat->gap_unmatched)) -  (0);    
        if( temp == (*h.access_main)(mat,i - 1,j - 0,UNMATCHED) )    {  
          *reti = i - 1; 
          *retj = j - 0; 
          *retstate = UNMATCHED; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-0,UNMATCHED);   
            }  
          return (*h.access_main)(mat,i - 1,j - 0,UNMATCHED);    
          }  
        temp = cscore - ((mat->ext_unmatched+mat->gap_unmatched)) -  (0);    
        if( temp == (*h.access_main)(mat,i - 0,j - 1,UNMATCHED) )    {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = UNMATCHED; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,UNMATCHED);   
            }  
          return (*h.access_main)(mat,i - 0,j - 1,UNMATCHED);    
          }  
        temp = cscore - (mat->ext_unmatched) -  (0); 
        if( temp == (*h.access_main)(mat,i - 1,j - 1,UNMATCHED) )    {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = UNMATCHED; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-1,UNMATCHED);   
            }  
          return (*h.access_main)(mat,i - 1,j - 1,UNMATCHED);    
          }  
        temp = cscore - (mat->open_unmatched) -  (0);    
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
        warn("Major problem (!) - in DnaProfileMat read off, position %d,%d state %d no source found!",i,j,state);   
        return (-1); 
      case INSERT :  
        temp = cscore - ((mat->q->col[i]->trans[DnaProfile_I2I]+mat->t->col[j]->trans[DnaProfile_D2D])) -  (0);  
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
        temp = cscore - ((mat->q->col[i]->trans[DnaProfile_M2I]+mat->t->col[j]->trans[DnaProfile_M2D])) -  (0);  
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
        warn("Major problem (!) - in DnaProfileMat read off, position %d,%d state %d no source found!",i,j,state);   
        return (-1); 
      case DELETE :  
        temp = cscore - ((mat->t->col[j]->trans[DnaProfile_I2I]+mat->q->col[i]->trans[DnaProfile_D2D])) -  (0);  
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
        temp = cscore - ((mat->t->col[j]->trans[DnaProfile_M2I]+mat->q->col[i]->trans[DnaProfile_M2D])) -  (0);  
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
        warn("Major problem (!) - in DnaProfileMat read off, position %d,%d state %d no source found!",i,j,state);   
        return (-1); 
      default:   
        warn("Major problem (!) - in DnaProfileMat read off, position %d,%d state %d no source found!",i,j,state);   
        return (-1); 
      } /* end of Switch state  */ 
}    


/* Function:  max_calc_special_DnaProfileMat(mat,i,j,state,isspecial,reti,retj,retstate,retspecial,cellscore,h)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [DnaProfileMat *]
 * Arg:                 i [UNKN ] Undocumented argument [int]
 * Arg:                 j [UNKN ] Undocumented argument [int]
 * Arg:             state [UNKN ] Undocumented argument [int]
 * Arg:         isspecial [UNKN ] Undocumented argument [boolean]
 * Arg:              reti [UNKN ] Undocumented argument [int *]
 * Arg:              retj [UNKN ] Undocumented argument [int *]
 * Arg:          retstate [UNKN ] Undocumented argument [int *]
 * Arg:        retspecial [UNKN ] Undocumented argument [boolean *]
 * Arg:         cellscore [UNKN ] Undocumented argument [int *]
 * Arg:                 h [UNKN ] Undocumented argument [DnaProfileMat_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int max_calc_special_DnaProfileMat(DnaProfileMat * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,DnaProfileMat_access_func_holder h) 
{
    register int temp;   
    register int cscore; 


    *reti = (*retj) = (*retstate) = DnaProfileMat_READ_OFF_ERROR;    


    if( j < 0 || j > mat->t->len)    {  
      warn("In DnaProfileMat matrix special read off - out of bounds on matrix [j is %d in special]",j); 
      return -1;     
      }  


    cscore = (*h.access_special)(mat,i,j,state); 
    switch(state)    { /*switch on special states*/ 
      case START :   
      case END :     
        /* source MATCH is from main matrix */ 
        for(i= mat->q->len-1;i >= 0 ;i--)    { /*for i >= 0*/ 
          temp = cscore - (0) - (0);     
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
        warn("Major problem (!) - in DnaProfileMat read off, position %d,%d state %d no source found  dropped into default on source switch!",i,j,state);    
        return (-1); 
      } /* end of switch on special states */ 
}    


/* Function:  calculate_DnaProfileMat(mat)
 *
 * Descrip:    This function calculates the DnaProfileMat matrix when in explicit mode
 *             To allocate the matrix use /allocate_Expl_DnaProfileMat
 *
 *
 * Arg:        mat [UNKN ] DnaProfileMat which contains explicit basematrix memory [DnaProfileMat *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean calculate_DnaProfileMat(DnaProfileMat * mat) 
{
    int i;   
    int j;   
    int leni;    
    int lenj;    
    int tot; 
    int num; 


    if( mat->basematrix->type != BASEMATRIX_TYPE_EXPLICIT )  {  
      warn("in calculate_DnaProfileMat, passed a non Explicit matrix type, cannot calculate!");  
      return FALSE;  
      }  


    leni = mat->leni;    
    lenj = mat->lenj;    
    tot = leni * lenj;   
    num = 0; 


    start_reporting("DnaProfileMat Matrix calculation: ");   
    for(j=0;j<lenj;j++)  {  
      auto int score;    
      auto int temp;     
      for(i=0;i<leni;i++)    {  
        if( num%1000 == 0 )  
          log_full_error(REPORT,0,"[%7d] Cells %2d%%%%",num,num*100/tot);    
        num++;   


        /* For state MATCH */ 
        /* setting first movement to score */ 
        score = DnaProfileMat_EXPL_MATRIX(mat,i-1,j-1,MATCH) + 0;    
        /* From state UNMATCHED to state MATCH */ 
        temp = DnaProfileMat_EXPL_MATRIX(mat,i-1,j-1,UNMATCHED) + 0;     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state MATCH */ 
        temp = DnaProfileMat_EXPL_MATRIX(mat,i-1,j-1,INSERT) + 0;    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state MATCH */ 
        temp = DnaProfileMat_EXPL_MATRIX(mat,i-1,j-1,DELETE) + 0;    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state START to state MATCH */ 
        temp = DnaProfileMat_EXPL_SPECIAL(mat,i-1,j-1,START) + 0;    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for MATCH */ 
        /* Add any movement independant score and put away */ 
         score += mat->m->score[i][j];   
         DnaProfileMat_EXPL_MATRIX(mat,i,j,MATCH) = score;   


        /* state MATCH is a source for special END */ 
        temp = score + (0) + (0) ;   
        if( temp > DnaProfileMat_EXPL_SPECIAL(mat,i,j,END) )     {  
          DnaProfileMat_EXPL_SPECIAL(mat,i,j,END) = temp;    
          }  




        /* Finished calculating state MATCH */ 


        /* For state UNMATCHED */ 
        /* setting first movement to score */ 
        score = DnaProfileMat_EXPL_MATRIX(mat,i-1,j-1,MATCH) + mat->open_unmatched;  
        /* From state UNMATCHED to state UNMATCHED */ 
        temp = DnaProfileMat_EXPL_MATRIX(mat,i-1,j-1,UNMATCHED) + mat->ext_unmatched;    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state UNMATCHED to state UNMATCHED */ 
        temp = DnaProfileMat_EXPL_MATRIX(mat,i-0,j-1,UNMATCHED) + (mat->ext_unmatched+mat->gap_unmatched);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state UNMATCHED to state UNMATCHED */ 
        temp = DnaProfileMat_EXPL_MATRIX(mat,i-1,j-0,UNMATCHED) + (mat->ext_unmatched+mat->gap_unmatched);   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for UNMATCHED */ 
        /* Add any movement independant score and put away */ 
         DnaProfileMat_EXPL_MATRIX(mat,i,j,UNMATCHED) = score;   


        /* Finished calculating state UNMATCHED */ 


        /* For state INSERT */ 
        /* setting first movement to score */ 
        score = DnaProfileMat_EXPL_MATRIX(mat,i-0,j-1,MATCH) + (mat->q->col[i]->trans[DnaProfile_M2I]+mat->t->col[j]->trans[DnaProfile_M2D]);    
        /* From state INSERT to state INSERT */ 
        temp = DnaProfileMat_EXPL_MATRIX(mat,i-0,j-1,INSERT) + (mat->q->col[i]->trans[DnaProfile_I2I]+mat->t->col[j]->trans[DnaProfile_D2D]);    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for INSERT */ 
        /* Add any movement independant score and put away */ 
         DnaProfileMat_EXPL_MATRIX(mat,i,j,INSERT) = score;  


        /* Finished calculating state INSERT */ 


        /* For state DELETE */ 
        /* setting first movement to score */ 
        score = DnaProfileMat_EXPL_MATRIX(mat,i-1,j-0,MATCH) + (mat->t->col[j]->trans[DnaProfile_M2I]+mat->q->col[i]->trans[DnaProfile_M2D]);    
        /* From state INSERT to state DELETE */ 
        temp = DnaProfileMat_EXPL_MATRIX(mat,i-1,j-0,INSERT) + (mat->t->col[j]->trans[DnaProfile_I2I]+mat->q->col[i]->trans[DnaProfile_D2D]);    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for DELETE */ 
        /* Add any movement independant score and put away */ 
         DnaProfileMat_EXPL_MATRIX(mat,i,j,DELETE) = score;  


        /* Finished calculating state DELETE */ 
        }  


      /* Special state START has no special to special movements */ 


      /* Special state END has no special to special movements */ 
      }  
    stop_reporting();    
    return TRUE;     
}    


/* Function:  calculate_dpenv_DnaProfileMat(mat,dpenv)
 *
 * Descrip:    This function calculates the DnaProfileMat matrix when in explicit mode, subject to the envelope
 *
 *
 * Arg:          mat [UNKN ] DnaProfileMat which contains explicit basematrix memory [DnaProfileMat *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean calculate_dpenv_DnaProfileMat(DnaProfileMat * mat,DPEnvelope * dpenv) 
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
      warn("in calculate_DnaProfileMat, passed a non Explicit matrix type, cannot calculate!");  
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
        DnaProfileMat_EXPL_MATRIX(mat,i,j,MATCH) = NEGI; 
        DnaProfileMat_EXPL_MATRIX(mat,i,j,UNMATCHED) = NEGI; 
        DnaProfileMat_EXPL_MATRIX(mat,i,j,INSERT) = NEGI;    
        DnaProfileMat_EXPL_MATRIX(mat,i,j,DELETE) = NEGI;    
        }  
      }  
    for(j=-1;j<mat->lenj;j++)    {  
      DnaProfileMat_EXPL_SPECIAL(mat,i,j,START) = 0; 
      DnaProfileMat_EXPL_SPECIAL(mat,i,j,END) = NEGI;    
      }  


    start_reporting("DnaProfileMat Matrix calculation: ");   
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
          DnaProfileMat_EXPL_MATRIX(mat,i,j,MATCH) = NEGI;   
          DnaProfileMat_EXPL_MATRIX(mat,i,j,UNMATCHED) = NEGI;   
          DnaProfileMat_EXPL_MATRIX(mat,i,j,INSERT) = NEGI;  
          DnaProfileMat_EXPL_MATRIX(mat,i,j,DELETE) = NEGI;  
          continue;  
          }  


        if( num%1000 == 0 )  
          log_full_error(REPORT,0,"[%7d] Cells %2d%%%%",num,num*100/tot);    
        num++;   


        /* For state MATCH */ 
        /* setting first movement to score */ 
        score = DnaProfileMat_EXPL_MATRIX(mat,i-1,j-1,MATCH) + 0;    
        /* From state UNMATCHED to state MATCH */ 
        temp = DnaProfileMat_EXPL_MATRIX(mat,i-1,j-1,UNMATCHED) + 0;     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state MATCH */ 
        temp = DnaProfileMat_EXPL_MATRIX(mat,i-1,j-1,INSERT) + 0;    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state MATCH */ 
        temp = DnaProfileMat_EXPL_MATRIX(mat,i-1,j-1,DELETE) + 0;    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state START to state MATCH */ 
        temp = DnaProfileMat_EXPL_SPECIAL(mat,i-1,j-1,START) + 0;    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for MATCH */ 
        /* Add any movement independant score and put away */ 
         score += mat->m->score[i][j];   
         DnaProfileMat_EXPL_MATRIX(mat,i,j,MATCH) = score;   


        /* state MATCH is a source for special END */ 
        temp = score + (0) + (0) ;   
        if( temp > DnaProfileMat_EXPL_SPECIAL(mat,i,j,END) )     {  
          DnaProfileMat_EXPL_SPECIAL(mat,i,j,END) = temp;    
          }  




        /* Finished calculating state MATCH */ 


        /* For state UNMATCHED */ 
        /* setting first movement to score */ 
        score = DnaProfileMat_EXPL_MATRIX(mat,i-1,j-1,MATCH) + mat->open_unmatched;  
        /* From state UNMATCHED to state UNMATCHED */ 
        temp = DnaProfileMat_EXPL_MATRIX(mat,i-1,j-1,UNMATCHED) + mat->ext_unmatched;    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state UNMATCHED to state UNMATCHED */ 
        temp = DnaProfileMat_EXPL_MATRIX(mat,i-0,j-1,UNMATCHED) + (mat->ext_unmatched+mat->gap_unmatched);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state UNMATCHED to state UNMATCHED */ 
        temp = DnaProfileMat_EXPL_MATRIX(mat,i-1,j-0,UNMATCHED) + (mat->ext_unmatched+mat->gap_unmatched);   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for UNMATCHED */ 
        /* Add any movement independant score and put away */ 
         DnaProfileMat_EXPL_MATRIX(mat,i,j,UNMATCHED) = score;   


        /* Finished calculating state UNMATCHED */ 


        /* For state INSERT */ 
        /* setting first movement to score */ 
        score = DnaProfileMat_EXPL_MATRIX(mat,i-0,j-1,MATCH) + (mat->q->col[i]->trans[DnaProfile_M2I]+mat->t->col[j]->trans[DnaProfile_M2D]);    
        /* From state INSERT to state INSERT */ 
        temp = DnaProfileMat_EXPL_MATRIX(mat,i-0,j-1,INSERT) + (mat->q->col[i]->trans[DnaProfile_I2I]+mat->t->col[j]->trans[DnaProfile_D2D]);    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for INSERT */ 
        /* Add any movement independant score and put away */ 
         DnaProfileMat_EXPL_MATRIX(mat,i,j,INSERT) = score;  


        /* Finished calculating state INSERT */ 


        /* For state DELETE */ 
        /* setting first movement to score */ 
        score = DnaProfileMat_EXPL_MATRIX(mat,i-1,j-0,MATCH) + (mat->t->col[j]->trans[DnaProfile_M2I]+mat->q->col[i]->trans[DnaProfile_M2D]);    
        /* From state INSERT to state DELETE */ 
        temp = DnaProfileMat_EXPL_MATRIX(mat,i-1,j-0,INSERT) + (mat->t->col[j]->trans[DnaProfile_I2I]+mat->q->col[i]->trans[DnaProfile_D2D]);    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for DELETE */ 
        /* Add any movement independant score and put away */ 
         DnaProfileMat_EXPL_MATRIX(mat,i,j,DELETE) = score;  


        /* Finished calculating state DELETE */ 
        }  


      /* Special state START has no special to special movements */ 


      /* Special state END has no special to special movements */ 
      }  
    stop_reporting();    
    return TRUE;     
}    


/* Function:  DnaProfileMat_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DnaProfileMat *]
 *
 */
DnaProfileMat * DnaProfileMat_alloc(void) 
{
    DnaProfileMat * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(DnaProfileMat *) ckalloc (sizeof(DnaProfileMat))) == NULL)  {  
      warn("DnaProfileMat_alloc failed ");   
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


/* Function:  free_DnaProfileMat(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [DnaProfileMat *]
 *
 * Return [UNKN ]  Undocumented return value [DnaProfileMat *]
 *
 */
DnaProfileMat * free_DnaProfileMat(DnaProfileMat * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a DnaProfileMat obj. Should be trappable"); 
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
    /* obj->open_unmatched is linked in */ 
    /* obj->ext_unmatched is linked in */ 
    /* obj->gap_unmatched is linked in */ 


    ckfree(obj); 
    return NULL; 
}    





#ifdef _cplusplus
}
#endif
