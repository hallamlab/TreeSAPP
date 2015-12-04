#ifdef _cplusplus
extern "C" {
#endif
#include "genomewise9.h"

# line 5 "genomewise9.c"


  /*****************   C functions  ****************/
  /*             Written using dynamite            */
  /*            Sat Sep  8 09:05:32 2007           */
  /*            email birney@sanger.ac.uk          */
  /* http://www.sanger.ac.uk/Users/birney/dynamite */
  /*************************************************/


  /* Please report any problems or bugs to         */
  /* Ewan Birney, birney@sanger.ac.uk              */


/* basic set of macros to map states to numbers */ 
#define UTR5 0   
#define UTR5_INTRON 1    
#define START_CODON 2    
#define CDS 3    
#define CDS_INTRON_0 4   
#define CDS_INTRON_1 5   
#define CDS_INTRON_2 6   
#define STOP_CODON 7 
#define UTR3 8   
#define UTR3_INTRON 9    


#define PREGENE_INTERGENIC 0 
#define POSTGENE_INTERGENIC 1    
#define INTERGENIC 2 
#define SPECIAL_UTR5 3   
#define SPECIAL_UTR3 4   
#define SPECIAL_CDS 5    
#define START 6  
#define END 7    


#define GenomeWise9_EXPL_MATRIX(this_matrix,i,j,STATE) this_matrix->basematrix->matrix[((j+10)*10)+STATE][i+0]   
#define GenomeWise9_EXPL_SPECIAL(matrix,i,j,STATE) matrix->basematrix->specmatrix[STATE][j+10]   
#define GenomeWise9_READ_OFF_ERROR -11
  


#define GenomeWise9_VSMALL_MATRIX(mat,i,j,STATE) mat->basematrix->matrix[(j+11)%11][((i+0)*10)+STATE]    
#define GenomeWise9_VSMALL_SPECIAL(mat,i,j,STATE) mat->basematrix->specmatrix[(j+11)%11][STATE]  




#define GenomeWise9_SHATTER_SPECIAL(matrix,i,j,STATE) matrix->shatter->special[STATE][j] 
#define GenomeWise9_SHATTER_MATRIX(matrix,i,j,STATE)  fetch_cell_value_ShatterMatrix(mat->shatter,i,j,STATE) 


/* Function:  PackAln_read_Shatter_GenomeWise9(mat)
 *
 * Descrip:    Reads off PackAln from shatter matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [GenomeWise9 *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * PackAln_read_Shatter_GenomeWise9(GenomeWise9 * mat) 
{
    GenomeWise9_access_func_holder holder;   


    holder.access_main    = GenomeWise9_shatter_access_main; 
    holder.access_special = GenomeWise9_shatter_access_special;  
    assert(mat);     
    assert(mat->shatter);    
    return PackAln_read_generic_GenomeWise9(mat,holder); 
}    


/* Function:  GenomeWise9_shatter_access_main(mat,i,j,state)
 *
 * Descrip: No Description
 *
 * Arg:          mat [UNKN ] Undocumented argument [GenomeWise9 *]
 * Arg:            i [UNKN ] Undocumented argument [int]
 * Arg:            j [UNKN ] Undocumented argument [int]
 * Arg:        state [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int GenomeWise9_shatter_access_main(GenomeWise9 * mat,int i,int j,int state) 
{
    return GenomeWise9_SHATTER_MATRIX(mat,i,j,state);    
}    


/* Function:  GenomeWise9_shatter_access_special(mat,i,j,state)
 *
 * Descrip: No Description
 *
 * Arg:          mat [UNKN ] Undocumented argument [GenomeWise9 *]
 * Arg:            i [UNKN ] Undocumented argument [int]
 * Arg:            j [UNKN ] Undocumented argument [int]
 * Arg:        state [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int GenomeWise9_shatter_access_special(GenomeWise9 * mat,int i,int j,int state) 
{
    return GenomeWise9_SHATTER_SPECIAL(mat,i,j,state);   
}    


/* Function:  calculate_shatter_GenomeWise9(mat,dpenv)
 *
 * Descrip:    This function calculates the GenomeWise9 matrix when in shatter mode
 *
 *
 * Arg:          mat [UNKN ] (null) [GenomeWise9 *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean calculate_shatter_GenomeWise9(GenomeWise9 * mat,DPEnvelope * dpenv) 
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
    int * SIG_0_1;   
    int * SIG_0_3;   
    int * SIG_0_6;   
    int * SIG_0_5;   
    int * SIG_0_4;   
    int * SIG_0_2;   
    int * SIG_0_8;   
    int * SIG_0_9;   
    int * SIG_0_10;  


    leni = mat->leni;    
    lenj = mat->lenj;    


    mat->shatter = new_ShatterMatrix(dpenv,10,lenj,8);   
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


    start_reporting("GenomeWise9 Matrix calculation: "); 
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
        SIG_0_1 = fetch_cell_from_ShatterMatrix(mat->shatter,i-0,j-1);   
        SIG_0_3 = fetch_cell_from_ShatterMatrix(mat->shatter,i-0,j-3);   
        SIG_0_6 = fetch_cell_from_ShatterMatrix(mat->shatter,i-0,j-6);   
        SIG_0_5 = fetch_cell_from_ShatterMatrix(mat->shatter,i-0,j-5);   
        SIG_0_4 = fetch_cell_from_ShatterMatrix(mat->shatter,i-0,j-4);   
        SIG_0_2 = fetch_cell_from_ShatterMatrix(mat->shatter,i-0,j-2);   
        SIG_0_8 = fetch_cell_from_ShatterMatrix(mat->shatter,i-0,j-8);   
        SIG_0_9 = fetch_cell_from_ShatterMatrix(mat->shatter,i-0,j-9);   
        SIG_0_10 = fetch_cell_from_ShatterMatrix(mat->shatter,i-0,j-10); 




        /* For state UTR5 */ 
        /* setting first movement to score */ 
        score = SIG_0_1[UTR5] + GNE_UTR(mat->evi,i,mat->gen,j);  
        /* From state UTR5_INTRON to state UTR5 */ 
        temp = SIG_0_1[UTR5_INTRON] + GNE_UTR_3SS(mat->evi,i,mat->gen,j);    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SPECIAL_UTR5 to state UTR5 */ 
        temp = GenomeWise9_SHATTER_SPECIAL(mat,i-0,j-1,SPECIAL_UTR5) + GNE_UTR(mat->evi,i,mat->gen,j);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state START to state UTR5 */ 
        temp = GenomeWise9_SHATTER_SPECIAL(mat,i-0,j-1,START) + GNE_UTR5_START(mat->evi,i,mat->gen,j);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTERGENIC to state UTR5 */ 
        temp = GenomeWise9_SHATTER_SPECIAL(mat,i-0,j-1,INTERGENIC) + (GNE_UTR(mat->evi,i,mat->gen,j)+GNE_UTR5_START(mat->evi,i,mat->gen,j));     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state PREGENE_INTERGENIC to state UTR5 */ 
        temp = GenomeWise9_SHATTER_SPECIAL(mat,i-0,j-1,PREGENE_INTERGENIC) + GNE_UTR(mat->evi,i,mat->gen,j);     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for UTR5 */ 
        /* Add any movement independant score and put away */ 
         SIG_0_0[UTR5] = score;  


        /* state UTR5 is a source for special SPECIAL_UTR5 */ 
        temp = score + (mat->switchcost) + (0) ;     
        if( temp > GenomeWise9_SHATTER_SPECIAL(mat,i,j,SPECIAL_UTR5) )   {  
          GenomeWise9_SHATTER_SPECIAL(mat,i,j,SPECIAL_UTR5) = temp;  
          }  




        /* Finished calculating state UTR5 */ 


        /* For state UTR5_INTRON */ 
        /* setting first movement to score */ 
        score = SIG_0_1[UTR5_INTRON] + GNE_UTR_INTRON(mat->evi,i,mat->gen,j);    
        /* From state UTR5 to state UTR5_INTRON */ 
        temp = SIG_0_1[UTR5] + GNE_UTR_5SS(mat->evi,i,mat->gen,j);   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for UTR5_INTRON */ 
        /* Add any movement independant score and put away */ 
         SIG_0_0[UTR5_INTRON] = score;   


        /* Finished calculating state UTR5_INTRON */ 


        /* For state START_CODON */ 
        /* setting first movement to score */ 
        score = SIG_0_3[UTR5_INTRON] + GNE_START_CODON(mat->evi,i,mat->gen,j);   
        /* From state UTR5 to state START_CODON */ 
        temp = SIG_0_3[UTR5] + GNE_START_CODON(mat->evi,i,mat->gen,j);   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for START_CODON */ 
        /* Add any movement independant score and put away */ 
         SIG_0_0[START_CODON] = score;   


        /* Finished calculating state START_CODON */ 


        /* For state CDS */ 
        /* setting first movement to score */ 
        score = SIG_0_3[CDS] + GNE_CDS(mat->evi,i,mat->gen,j);   
        /* From state CDS_INTRON_0 to state CDS */ 
        temp = SIG_0_6[CDS_INTRON_0] + GNE_CDS_3SS(mat->evi,i,mat->gen,(j-3),0);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state CDS_INTRON_1 to state CDS */ 
        temp = SIG_0_5[CDS_INTRON_1] + GNE_CDS_3SS(mat->evi,i,mat->gen,(j-2),1);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state CDS_INTRON_2 to state CDS */ 
        temp = SIG_0_4[CDS_INTRON_2] + GNE_CDS_3SS(mat->evi,i,mat->gen,(j-1),2);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state CDS to state CDS */ 
        temp = SIG_0_2[CDS] + GNE_CDS_FRAMESHIFT(mat->evi,i,mat->gen,j,2);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state CDS to state CDS */ 
        temp = SIG_0_4[CDS] + GNE_CDS_FRAMESHIFT(mat->evi,i,mat->gen,j,4);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state UTR5 to state CDS */ 
        temp = SIG_0_3[UTR5] + (GNE_CDS(mat->evi,i,mat->gen,j)+mat->non_start_codon);    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state START_CODON to state CDS */ 
        temp = SIG_0_3[START_CODON] + GNE_CDS(mat->evi,i,mat->gen,j);    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SPECIAL_CDS to state CDS */ 
        temp = GenomeWise9_SHATTER_SPECIAL(mat,i-0,j-3,SPECIAL_CDS) + GNE_CDS(mat->evi,i,mat->gen,j);    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTERGENIC to state CDS */ 
        temp = GenomeWise9_SHATTER_SPECIAL(mat,i-0,j-3,INTERGENIC) + (GNE_CDS(mat->evi,i,mat->gen,j)+GNE_UTR5_START(mat->evi,i,mat->gen,j));     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for CDS */ 
        /* Add any movement independant score and put away */ 
         SIG_0_0[CDS] = score;   


        /* state CDS is a source for special INTERGENIC */ 
        temp = score + ((mat->newgenecost+GNE_UTR3_END(mat->evi,i,mat->gen,j))) + (0) ;  
        if( temp > GenomeWise9_SHATTER_SPECIAL(mat,i,j,INTERGENIC) )     {  
          GenomeWise9_SHATTER_SPECIAL(mat,i,j,INTERGENIC) = temp;    
          }  




        /* state CDS is a source for special SPECIAL_CDS */ 
        temp = score + ((mat->switchcost+mat->rndcodon->codon[CSEQ_GENOMIC_CODON(mat->gen,j)])) + (0) ;  
        if( temp > GenomeWise9_SHATTER_SPECIAL(mat,i,j,SPECIAL_CDS) )    {  
          GenomeWise9_SHATTER_SPECIAL(mat,i,j,SPECIAL_CDS) = temp;   
          }  




        /* Finished calculating state CDS */ 


        /* For state CDS_INTRON_0 */ 
        /* setting first movement to score */ 
        score = SIG_0_1[CDS_INTRON_0] + GNE_CDS_INTRON(mat->evi,i,mat->gen,j);   
        /* From state CDS to state CDS_INTRON_0 */ 
        temp = SIG_0_8[CDS] + GNE_CDS_5SS(mat->evi,i,mat->gen,(j-7),0);  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for CDS_INTRON_0 */ 
        /* Add any movement independant score and put away */ 
         SIG_0_0[CDS_INTRON_0] = score;  


        /* Finished calculating state CDS_INTRON_0 */ 


        /* For state CDS_INTRON_1 */ 
        /* setting first movement to score */ 
        score = SIG_0_1[CDS_INTRON_1] + GNE_CDS_INTRON(mat->evi,i,mat->gen,j);   
        /* From state CDS to state CDS_INTRON_1 */ 
        temp = SIG_0_9[CDS] + GNE_CDS_5SS(mat->evi,i,mat->gen,(j-7),1);  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for CDS_INTRON_1 */ 
        /* Add any movement independant score and put away */ 
         SIG_0_0[CDS_INTRON_1] = score;  


        /* Finished calculating state CDS_INTRON_1 */ 


        /* For state CDS_INTRON_2 */ 
        /* setting first movement to score */ 
        score = SIG_0_1[CDS_INTRON_2] + GNE_CDS_INTRON(mat->evi,i,mat->gen,j);   
        /* From state CDS to state CDS_INTRON_2 */ 
        temp = SIG_0_10[CDS] + GNE_CDS_5SS(mat->evi,i,mat->gen,(j-7),2);     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for CDS_INTRON_2 */ 
        /* Add any movement independant score and put away */ 
         SIG_0_0[CDS_INTRON_2] = score;  


        /* Finished calculating state CDS_INTRON_2 */ 


        /* For state STOP_CODON */ 
        /* setting first movement to score */ 
        score = SIG_0_3[CDS] + GNE_STOP_CODON(mat->evi,i,mat->gen,j);    


        /* Ok - finished max calculation for STOP_CODON */ 
        /* Add any movement independant score and put away */ 
         SIG_0_0[STOP_CODON] = score;    


        /* Finished calculating state STOP_CODON */ 


        /* For state UTR3 */ 
        /* setting first movement to score */ 
        score = SIG_0_1[UTR3] + GNE_UTR(mat->evi,i,mat->gen,j);  
        /* From state CDS to state UTR3 */ 
        temp = SIG_0_1[CDS] + (GNE_UTR(mat->evi,i,mat->gen,j)+mat->non_stop_codon);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state STOP_CODON to state UTR3 */ 
        temp = SIG_0_1[STOP_CODON] + GNE_UTR(mat->evi,i,mat->gen,j);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state UTR3_INTRON to state UTR3 */ 
        temp = SIG_0_1[UTR3_INTRON] + GNE_UTR_3SS(mat->evi,i,mat->gen,j);    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* Has restricted position */ 
        if( (j-1) == 0  )    {  
          /* From state INTERGENIC to state UTR3 */ 
          temp = GenomeWise9_SHATTER_SPECIAL(mat,i-0,j-1,INTERGENIC) + GNE_UTR(mat->evi,i,mat->gen,j);   
          if( temp  > score )    {  
            score = temp;    
            }  
          }  


        /* Ok - finished max calculation for UTR3 */ 
        /* Add any movement independant score and put away */ 
         SIG_0_0[UTR3] = score;  


        /* state UTR3 is a source for special POSTGENE_INTERGENIC */ 
        temp = score + ((mat->newgenecost+GNE_UTR3_END(mat->evi,i,mat->gen,j))) + (0) ;  
        if( temp > GenomeWise9_SHATTER_SPECIAL(mat,i,j,POSTGENE_INTERGENIC) )    {  
          GenomeWise9_SHATTER_SPECIAL(mat,i,j,POSTGENE_INTERGENIC) = temp;   
          }  




        /* state UTR3 is a source for special INTERGENIC */ 
        temp = score + ((mat->newgenecost+GNE_UTR3_END(mat->evi,i,mat->gen,j))) + (0) ;  
        if( temp > GenomeWise9_SHATTER_SPECIAL(mat,i,j,INTERGENIC) )     {  
          GenomeWise9_SHATTER_SPECIAL(mat,i,j,INTERGENIC) = temp;    
          }  




        /* state UTR3 is a source for special SPECIAL_UTR3 */ 
        temp = score + (mat->switchcost) + (0) ;     
        if( temp > GenomeWise9_SHATTER_SPECIAL(mat,i,j,SPECIAL_UTR3) )   {  
          GenomeWise9_SHATTER_SPECIAL(mat,i,j,SPECIAL_UTR3) = temp;  
          }  




        /* Finished calculating state UTR3 */ 


        /* For state UTR3_INTRON */ 
        /* setting first movement to score */ 
        score = SIG_0_1[UTR3_INTRON] + GNE_UTR_INTRON(mat->evi,i,mat->gen,j);    
        /* From state UTR3 to state UTR3_INTRON */ 
        temp = SIG_0_1[UTR3] + GNE_UTR_5SS(mat->evi,i,mat->gen,j);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SPECIAL_UTR3 to state UTR3_INTRON */ 
        temp = GenomeWise9_SHATTER_SPECIAL(mat,i-0,j-1,SPECIAL_UTR3) + GNE_UTR(mat->evi,i,mat->gen,j);   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for UTR3_INTRON */ 
        /* Add any movement independant score and put away */ 
         SIG_0_0[UTR3_INTRON] = score;   


        /* Finished calculating state UTR3_INTRON */ 
        }  


      /* Special state PREGENE_INTERGENIC has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GenomeWise9_SHATTER_SPECIAL(mat,0,j,PREGENE_INTERGENIC);   


      /* Source START is a special source for PREGENE_INTERGENIC */ 
      temp = GenomeWise9_SHATTER_SPECIAL(mat,0,j - 1,START) + (0) + (0);     
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      GenomeWise9_SHATTER_SPECIAL(mat,0,j,PREGENE_INTERGENIC) = score;   
      /* Finished updating state PREGENE_INTERGENIC */ 




      /* Special state POSTGENE_INTERGENIC has no special to special movements */ 


      /* Special state INTERGENIC has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GenomeWise9_SHATTER_SPECIAL(mat,0,j,INTERGENIC);   


      /* Source INTERGENIC is a special source for INTERGENIC */ 
      temp = GenomeWise9_SHATTER_SPECIAL(mat,0,j - 1,INTERGENIC) + (0) + (0);    
      if( temp > score ) 
        score = temp;    


      /* Source CDS for state INTERGENIC is not special... already calculated */ 
      /* Source UTR3 for state INTERGENIC is not special... already calculated */ 
      /* Put back score... (now updated!) */ 
      GenomeWise9_SHATTER_SPECIAL(mat,0,j,INTERGENIC) = score;   
      /* Finished updating state INTERGENIC */ 




      /* Special state SPECIAL_UTR5 has no special to special movements */ 


      /* Special state SPECIAL_UTR3 has no special to special movements */ 


      /* Special state SPECIAL_CDS has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GenomeWise9_SHATTER_SPECIAL(mat,0,j,SPECIAL_CDS);  


      /* Source CDS for state SPECIAL_CDS is not special... already calculated */ 
      /* Source SPECIAL_CDS is a special source for SPECIAL_CDS */ 
      temp = GenomeWise9_SHATTER_SPECIAL(mat,0,j - 3,SPECIAL_CDS) + (mat->rndcodon->codon[CSEQ_GENOMIC_CODON(mat->gen,j)]) + (0);    
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      GenomeWise9_SHATTER_SPECIAL(mat,0,j,SPECIAL_CDS) = score;  
      /* Finished updating state SPECIAL_CDS */ 




      /* Special state START has no special to special movements */ 


      /* Special state END has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GenomeWise9_SHATTER_SPECIAL(mat,0,j,END);  


      /* Source INTERGENIC is a special source for END */ 
      /* Has restricted position */ 
      if( j == mat->lenj-1 ) {  
        temp = GenomeWise9_SHATTER_SPECIAL(mat,0,j - 1,INTERGENIC) + (0) + (0);  
        if( temp > score )   
          score = temp;  
        }  


      /* Put back score... (now updated!) */ 
      GenomeWise9_SHATTER_SPECIAL(mat,0,j,END) = score;  
      /* Finished updating state END */ 


      }  
    stop_reporting();    
    return TRUE;     
}    


/* Function:  search_GenomeWise9(dbsi,out,evi,targetdb,switchcost,newgenecost,non_start_codon,non_stop_codon,rndcodon)
 *
 * Descrip:    This function makes a database search of GenomeWise9
 *             It uses the dbsi structure to choose which implementation
 *             to use of the database searching. This way at run time you
 *             can switch between single threaded/multi-threaded or hardware
 *
 *
 * Arg:                   dbsi [UNKN ] Undocumented argument [DBSearchImpl *]
 * Arg:                    out [UNKN ] Undocumented argument [Hscore *]
 * Arg:                    evi [UNKN ] Undocumented argument [GenomeEvidenceSet*]
 * Arg:               targetdb [UNKN ] Undocumented argument [GenomicDB*]
 * Arg:             switchcost [UNKN ] Undocumented argument [int]
 * Arg:            newgenecost [UNKN ] Undocumented argument [int]
 * Arg:        non_start_codon [UNKN ] Undocumented argument [int]
 * Arg:         non_stop_codon [UNKN ] Undocumented argument [int]
 * Arg:               rndcodon [UNKN ] Undocumented argument [RandomCodonScore *]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type search_GenomeWise9(DBSearchImpl * dbsi,Hscore * out,GenomeEvidenceSet* evi,GenomicDB* targetdb ,int switchcost,int newgenecost,int non_start_codon,int non_stop_codon,RandomCodonScore * rndcodon) 
{
#ifdef PTHREAD   
    int i;   
    int thr_no;  
    pthread_attr_t pat;  
    struct thread_pool_holder_GenomeWise9 * holder;  
#endif   
    if( out == NULL )    {  
      warn("Passed in a null Hscore object into search_GenomeWise9. Can't process results!");    
      return SEARCH_ERROR;   
      }  
    if( dbsi == NULL )   {  
      warn("Passed in a null DBSearchImpl object into search_GenomeWise9. Can't process results!");  
      return SEARCH_ERROR;   
      }  
    if( dbsi->trace_level > 5 )  
      warn("Asking for trace level of %d in database search for GenomeWise9, but it was compiled with a trace level of -2139062144. Not all trace statements can be shown",dbsi->trace_level);   
    switch(dbsi->type)   { /*switch on implementation*/ 
      case DBSearchImpl_Serial : 
        return serial_search_GenomeWise9(out,evi, targetdb ,switchcost,newgenecost,non_start_codon,non_stop_codon,rndcodon); 
      case DBSearchImpl_Pthreads :   
#ifdef PTHREAD   
        holder = (struct thread_pool_holder_GenomeWise9 *) ckalloc(sizeof(struct thread_pool_holder_GenomeWise9));   
        if( holder == NULL )     {  
          warn("Unable to allocated thread pool datastructure...");  
          return SEARCH_ERROR;   
          }  
        holder->out = out;   
        holder->dbsi = dbsi; 
        holder->evi = evi;   
        holder->targetdb = targetdb; 
        holder->switchcost = switchcost; 
        holder->newgenecost = newgenecost;   
        holder->non_start_codon = non_start_codon;   
        holder->non_stop_codon = non_stop_codon; 
        holder->rndcodon = rndcodon; 
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
          if( pthread_create(holder->pool+i,&pat,thread_loop_GenomeWise9,(void *)holder) )   
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
        warn("You did not specifiy the PTHREAD compile when compiled the C code for GenomeWise9");   
#endif /* finished threads */    
      default :  
        warn("database search implementation %s was not provided in the compiled dynamite file from GenomeWise9",impl_string_DBSearchImpl(dbsi));    
        return SEARCH_ERROR; 
      } /* end of switch on implementation */ 


}    


/* Function:  thread_loop_GenomeWise9(ptr)
 *
 * Descrip:    dummy loop code foreach thread for GenomeWise9
 *
 *
 * Arg:        ptr [UNKN ] Undocumented argument [void *]
 *
 * Return [UNKN ]  Undocumented return value [void *]
 *
 */
void * thread_loop_GenomeWise9(void * ptr) 
{
    fatal("dummy thread loop function"); 
}    


/* Function:  serial_search_GenomeWise9(out,evi,targetdb,switchcost,newgenecost,non_start_codon,non_stop_codon,rndcodon)
 *
 * Descrip:    This function makes a database search of GenomeWise9
 *             It is a single processor implementation
 *
 *
 * Arg:                    out [UNKN ] Undocumented argument [Hscore *]
 * Arg:                    evi [UNKN ] Undocumented argument [GenomeEvidenceSet*]
 * Arg:               targetdb [UNKN ] Undocumented argument [GenomicDB*]
 * Arg:             switchcost [UNKN ] Undocumented argument [int]
 * Arg:            newgenecost [UNKN ] Undocumented argument [int]
 * Arg:        non_start_codon [UNKN ] Undocumented argument [int]
 * Arg:         non_stop_codon [UNKN ] Undocumented argument [int]
 * Arg:               rndcodon [UNKN ] Undocumented argument [RandomCodonScore *]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type serial_search_GenomeWise9(Hscore * out,GenomeEvidenceSet* evi,GenomicDB* targetdb ,int switchcost,int newgenecost,int non_start_codon,int non_stop_codon,RandomCodonScore * rndcodon) 
{
    ComplexSequence* gen;    
    int db_status;   
    int score;   
    int query_pos = 0;   
    int target_pos = 0;  
    DataScore * ds;  


    push_errormsg_stack("Before any actual search in db searching"); 


    target_pos = 0;  


    gen = init_GenomicDB(targetdb,&db_status);   
    if( db_status == DB_RETURN_ERROR )   {  
      warn("In searching GenomeWise9, got a database init error on the target [gen] database");  
      return SEARCH_ERROR;   
      }  
    for(;;)  { /*For all target entries*/ 


      /* No maximum length - allocated on-the-fly */ 
      score = score_only_GenomeWise9(evi, gen , switchcost, newgenecost, non_start_codon, non_stop_codon, rndcodon);     
      if( should_store_Hscore(out,score) == TRUE )   { /*if storing datascore*/ 
        ds = new_DataScore_from_storage(out);    
        if( ds == NULL )     {  
          warn("GenomeWise9 search had a memory error in allocating a new_DataScore (?a leak somewhere - DataScore is a very small datastructure");  
          return SEARCH_ERROR;   
          }  
        /* Now: add query/target information to the entry */ 
        dataentry_add_GenomicDB(ds->target,gen,targetdb);    
        ds->score = score;   
        add_Hscore(out,ds);  
        } /* end of if storing datascore */ 
      pop_errormsg_stack();  
      push_errormsg_stack("DB searching: just finished [Query Pos: %d] [Target Pos: %d]",query_pos,target_pos);  


       gen = reload_GenomicDB(gen,targetdb,&db_status);  
      if( db_status == DB_RETURN_ERROR ) {  
        warn("In searching GenomeWise9, Reload error on database gen, position %d,%d",query_pos,target_pos); 
        return SEARCH_ERROR; 
        }  
      if( db_status == DB_RETURN_END )   
        break;  /* Out of target loop */ 
      target_pos++;  
      } /* end of For all target entries */ 
    close_GenomicDB(gen,targetdb);   
    pop_errormsg_stack();    
    return SEARCH_OK;    
}    


/* Function:  score_only_GenomeWise9(evi,gen,switchcost,newgenecost,non_start_codon,non_stop_codon,rndcodon)
 *
 * Descrip:    This function just calculates the score for the matrix
 *             I am pretty sure we can do this better, but hey, for the moment...
 *             It calls /allocate_GenomeWise9_only
 *
 *
 * Arg:                    evi [UNKN ] query data structure [GenomeEvidenceSet*]
 * Arg:                    gen [UNKN ] target data structure [ComplexSequence*]
 * Arg:             switchcost [UNKN ] Resource [int]
 * Arg:            newgenecost [UNKN ] Resource [int]
 * Arg:        non_start_codon [UNKN ] Resource [int]
 * Arg:         non_stop_codon [UNKN ] Resource [int]
 * Arg:               rndcodon [UNKN ] Resource [RandomCodonScore *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int score_only_GenomeWise9(GenomeEvidenceSet* evi,ComplexSequence* gen ,int switchcost,int newgenecost,int non_start_codon,int non_stop_codon,RandomCodonScore * rndcodon) 
{
    int bestscore = NEGI;    
    int i;   
    int j;   
    int k;   
    GenomeWise9 * mat;   


    mat = allocate_GenomeWise9_only(evi, gen , switchcost, newgenecost, non_start_codon, non_stop_codon, rndcodon);  
    if( mat == NULL )    {  
      warn("Memory allocation error in the db search - unable to communicate to calling function. this spells DIASTER!");    
      return NEGI;   
      }  
    if((mat->basematrix = BaseMatrix_alloc_matrix_and_specials(11,(mat->leni + 0) * 10,11,8)) == NULL)   {  
      warn("Score only matrix for GenomeWise9 cannot be allocated, (asking for 10  by %d  cells)",mat->leni*10); 
      mat = free_GenomeWise9(mat);   
      return 0;  
      }  
    mat->basematrix->type = BASEMATRIX_TYPE_VERYSMALL;   


    /* Now, initiate matrix */ 
    for(j=0;j<12;j++)    {  
      for(i=(-0);i<mat->leni;i++)    {  
        for(k=0;k<10;k++)    
          GenomeWise9_VSMALL_MATRIX(mat,i,j,k) = NEGI;   
        }  
      GenomeWise9_VSMALL_SPECIAL(mat,i,j,PREGENE_INTERGENIC) = NEGI; 
      GenomeWise9_VSMALL_SPECIAL(mat,i,j,POSTGENE_INTERGENIC) = NEGI;    
      GenomeWise9_VSMALL_SPECIAL(mat,i,j,INTERGENIC) = NEGI; 
      GenomeWise9_VSMALL_SPECIAL(mat,i,j,SPECIAL_UTR5) = NEGI;   
      GenomeWise9_VSMALL_SPECIAL(mat,i,j,SPECIAL_UTR3) = NEGI;   
      GenomeWise9_VSMALL_SPECIAL(mat,i,j,SPECIAL_CDS) = NEGI;    
      GenomeWise9_VSMALL_SPECIAL(mat,i,j,START) = 0; 
      GenomeWise9_VSMALL_SPECIAL(mat,i,j,END) = NEGI;    
      }  


    /* Ok, lets do-o-o-o-o it */ 


    for(j=0;j<mat->lenj;j++) { /*for all target positions*/ 
      auto int score;    
      auto int temp;     
      for(i=0;i<mat->leni;i++)   { /*for all query positions*/ 


        /* For state UTR5 */ 
        /* setting first movement to score */ 
        score = GenomeWise9_VSMALL_MATRIX(mat,i-0,j-1,UTR5) + GNE_UTR(mat->evi,i,mat->gen,j);    
        /* From state UTR5_INTRON to state UTR5 */ 
        temp = GenomeWise9_VSMALL_MATRIX(mat,i-0,j-1,UTR5_INTRON) + GNE_UTR_3SS(mat->evi,i,mat->gen,j);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SPECIAL_UTR5 to state UTR5 */ 
        temp = GenomeWise9_VSMALL_SPECIAL(mat,i-0,j-1,SPECIAL_UTR5) + GNE_UTR(mat->evi,i,mat->gen,j);    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state START to state UTR5 */ 
        temp = GenomeWise9_VSMALL_SPECIAL(mat,i-0,j-1,START) + GNE_UTR5_START(mat->evi,i,mat->gen,j);    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTERGENIC to state UTR5 */ 
        temp = GenomeWise9_VSMALL_SPECIAL(mat,i-0,j-1,INTERGENIC) + (GNE_UTR(mat->evi,i,mat->gen,j)+GNE_UTR5_START(mat->evi,i,mat->gen,j));  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state PREGENE_INTERGENIC to state UTR5 */ 
        temp = GenomeWise9_VSMALL_SPECIAL(mat,i-0,j-1,PREGENE_INTERGENIC) + GNE_UTR(mat->evi,i,mat->gen,j);  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for UTR5 */ 
        /* Add any movement independant score and put away */ 
         GenomeWise9_VSMALL_MATRIX(mat,i,j,UTR5) = score;    


        /* state UTR5 is a source for special SPECIAL_UTR5 */ 
        temp = score + (mat->switchcost) + (0) ;     
        if( temp > GenomeWise9_VSMALL_SPECIAL(mat,i,j,SPECIAL_UTR5) )    {  
          GenomeWise9_VSMALL_SPECIAL(mat,i,j,SPECIAL_UTR5) = temp;   
          }  




        /* Finished calculating state UTR5 */ 


        /* For state UTR5_INTRON */ 
        /* setting first movement to score */ 
        score = GenomeWise9_VSMALL_MATRIX(mat,i-0,j-1,UTR5_INTRON) + GNE_UTR_INTRON(mat->evi,i,mat->gen,j);  
        /* From state UTR5 to state UTR5_INTRON */ 
        temp = GenomeWise9_VSMALL_MATRIX(mat,i-0,j-1,UTR5) + GNE_UTR_5SS(mat->evi,i,mat->gen,j);     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for UTR5_INTRON */ 
        /* Add any movement independant score and put away */ 
         GenomeWise9_VSMALL_MATRIX(mat,i,j,UTR5_INTRON) = score; 


        /* Finished calculating state UTR5_INTRON */ 


        /* For state START_CODON */ 
        /* setting first movement to score */ 
        score = GenomeWise9_VSMALL_MATRIX(mat,i-0,j-3,UTR5_INTRON) + GNE_START_CODON(mat->evi,i,mat->gen,j);     
        /* From state UTR5 to state START_CODON */ 
        temp = GenomeWise9_VSMALL_MATRIX(mat,i-0,j-3,UTR5) + GNE_START_CODON(mat->evi,i,mat->gen,j);     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for START_CODON */ 
        /* Add any movement independant score and put away */ 
         GenomeWise9_VSMALL_MATRIX(mat,i,j,START_CODON) = score; 


        /* Finished calculating state START_CODON */ 


        /* For state CDS */ 
        /* setting first movement to score */ 
        score = GenomeWise9_VSMALL_MATRIX(mat,i-0,j-3,CDS) + GNE_CDS(mat->evi,i,mat->gen,j);     
        /* From state CDS_INTRON_0 to state CDS */ 
        temp = GenomeWise9_VSMALL_MATRIX(mat,i-0,j-6,CDS_INTRON_0) + GNE_CDS_3SS(mat->evi,i,mat->gen,(j-3),0);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state CDS_INTRON_1 to state CDS */ 
        temp = GenomeWise9_VSMALL_MATRIX(mat,i-0,j-5,CDS_INTRON_1) + GNE_CDS_3SS(mat->evi,i,mat->gen,(j-2),1);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state CDS_INTRON_2 to state CDS */ 
        temp = GenomeWise9_VSMALL_MATRIX(mat,i-0,j-4,CDS_INTRON_2) + GNE_CDS_3SS(mat->evi,i,mat->gen,(j-1),2);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state CDS to state CDS */ 
        temp = GenomeWise9_VSMALL_MATRIX(mat,i-0,j-2,CDS) + GNE_CDS_FRAMESHIFT(mat->evi,i,mat->gen,j,2);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state CDS to state CDS */ 
        temp = GenomeWise9_VSMALL_MATRIX(mat,i-0,j-4,CDS) + GNE_CDS_FRAMESHIFT(mat->evi,i,mat->gen,j,4);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state UTR5 to state CDS */ 
        temp = GenomeWise9_VSMALL_MATRIX(mat,i-0,j-3,UTR5) + (GNE_CDS(mat->evi,i,mat->gen,j)+mat->non_start_codon);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state START_CODON to state CDS */ 
        temp = GenomeWise9_VSMALL_MATRIX(mat,i-0,j-3,START_CODON) + GNE_CDS(mat->evi,i,mat->gen,j);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SPECIAL_CDS to state CDS */ 
        temp = GenomeWise9_VSMALL_SPECIAL(mat,i-0,j-3,SPECIAL_CDS) + GNE_CDS(mat->evi,i,mat->gen,j);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTERGENIC to state CDS */ 
        temp = GenomeWise9_VSMALL_SPECIAL(mat,i-0,j-3,INTERGENIC) + (GNE_CDS(mat->evi,i,mat->gen,j)+GNE_UTR5_START(mat->evi,i,mat->gen,j));  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for CDS */ 
        /* Add any movement independant score and put away */ 
         GenomeWise9_VSMALL_MATRIX(mat,i,j,CDS) = score; 


        /* state CDS is a source for special INTERGENIC */ 
        temp = score + ((mat->newgenecost+GNE_UTR3_END(mat->evi,i,mat->gen,j))) + (0) ;  
        if( temp > GenomeWise9_VSMALL_SPECIAL(mat,i,j,INTERGENIC) )  {  
          GenomeWise9_VSMALL_SPECIAL(mat,i,j,INTERGENIC) = temp;     
          }  




        /* state CDS is a source for special SPECIAL_CDS */ 
        temp = score + ((mat->switchcost+mat->rndcodon->codon[CSEQ_GENOMIC_CODON(mat->gen,j)])) + (0) ;  
        if( temp > GenomeWise9_VSMALL_SPECIAL(mat,i,j,SPECIAL_CDS) )     {  
          GenomeWise9_VSMALL_SPECIAL(mat,i,j,SPECIAL_CDS) = temp;    
          }  




        /* Finished calculating state CDS */ 


        /* For state CDS_INTRON_0 */ 
        /* setting first movement to score */ 
        score = GenomeWise9_VSMALL_MATRIX(mat,i-0,j-1,CDS_INTRON_0) + GNE_CDS_INTRON(mat->evi,i,mat->gen,j);     
        /* From state CDS to state CDS_INTRON_0 */ 
        temp = GenomeWise9_VSMALL_MATRIX(mat,i-0,j-8,CDS) + GNE_CDS_5SS(mat->evi,i,mat->gen,(j-7),0);    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for CDS_INTRON_0 */ 
        /* Add any movement independant score and put away */ 
         GenomeWise9_VSMALL_MATRIX(mat,i,j,CDS_INTRON_0) = score;    


        /* Finished calculating state CDS_INTRON_0 */ 


        /* For state CDS_INTRON_1 */ 
        /* setting first movement to score */ 
        score = GenomeWise9_VSMALL_MATRIX(mat,i-0,j-1,CDS_INTRON_1) + GNE_CDS_INTRON(mat->evi,i,mat->gen,j);     
        /* From state CDS to state CDS_INTRON_1 */ 
        temp = GenomeWise9_VSMALL_MATRIX(mat,i-0,j-9,CDS) + GNE_CDS_5SS(mat->evi,i,mat->gen,(j-7),1);    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for CDS_INTRON_1 */ 
        /* Add any movement independant score and put away */ 
         GenomeWise9_VSMALL_MATRIX(mat,i,j,CDS_INTRON_1) = score;    


        /* Finished calculating state CDS_INTRON_1 */ 


        /* For state CDS_INTRON_2 */ 
        /* setting first movement to score */ 
        score = GenomeWise9_VSMALL_MATRIX(mat,i-0,j-1,CDS_INTRON_2) + GNE_CDS_INTRON(mat->evi,i,mat->gen,j);     
        /* From state CDS to state CDS_INTRON_2 */ 
        temp = GenomeWise9_VSMALL_MATRIX(mat,i-0,j-10,CDS) + GNE_CDS_5SS(mat->evi,i,mat->gen,(j-7),2);   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for CDS_INTRON_2 */ 
        /* Add any movement independant score and put away */ 
         GenomeWise9_VSMALL_MATRIX(mat,i,j,CDS_INTRON_2) = score;    


        /* Finished calculating state CDS_INTRON_2 */ 


        /* For state STOP_CODON */ 
        /* setting first movement to score */ 
        score = GenomeWise9_VSMALL_MATRIX(mat,i-0,j-3,CDS) + GNE_STOP_CODON(mat->evi,i,mat->gen,j);  


        /* Ok - finished max calculation for STOP_CODON */ 
        /* Add any movement independant score and put away */ 
         GenomeWise9_VSMALL_MATRIX(mat,i,j,STOP_CODON) = score;  


        /* Finished calculating state STOP_CODON */ 


        /* For state UTR3 */ 
        /* setting first movement to score */ 
        score = GenomeWise9_VSMALL_MATRIX(mat,i-0,j-1,UTR3) + GNE_UTR(mat->evi,i,mat->gen,j);    
        /* From state CDS to state UTR3 */ 
        temp = GenomeWise9_VSMALL_MATRIX(mat,i-0,j-1,CDS) + (GNE_UTR(mat->evi,i,mat->gen,j)+mat->non_stop_codon);    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state STOP_CODON to state UTR3 */ 
        temp = GenomeWise9_VSMALL_MATRIX(mat,i-0,j-1,STOP_CODON) + GNE_UTR(mat->evi,i,mat->gen,j);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state UTR3_INTRON to state UTR3 */ 
        temp = GenomeWise9_VSMALL_MATRIX(mat,i-0,j-1,UTR3_INTRON) + GNE_UTR_3SS(mat->evi,i,mat->gen,j);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* Has restricted position */ 
        if( (j-1) == 0  )    {  
          /* From state INTERGENIC to state UTR3 */ 
          temp = GenomeWise9_VSMALL_SPECIAL(mat,i-0,j-1,INTERGENIC) + GNE_UTR(mat->evi,i,mat->gen,j);    
          if( temp  > score )    {  
            score = temp;    
            }  
          }  


        /* Ok - finished max calculation for UTR3 */ 
        /* Add any movement independant score and put away */ 
         GenomeWise9_VSMALL_MATRIX(mat,i,j,UTR3) = score;    


        /* state UTR3 is a source for special POSTGENE_INTERGENIC */ 
        temp = score + ((mat->newgenecost+GNE_UTR3_END(mat->evi,i,mat->gen,j))) + (0) ;  
        if( temp > GenomeWise9_VSMALL_SPECIAL(mat,i,j,POSTGENE_INTERGENIC) )     {  
          GenomeWise9_VSMALL_SPECIAL(mat,i,j,POSTGENE_INTERGENIC) = temp;    
          }  




        /* state UTR3 is a source for special INTERGENIC */ 
        temp = score + ((mat->newgenecost+GNE_UTR3_END(mat->evi,i,mat->gen,j))) + (0) ;  
        if( temp > GenomeWise9_VSMALL_SPECIAL(mat,i,j,INTERGENIC) )  {  
          GenomeWise9_VSMALL_SPECIAL(mat,i,j,INTERGENIC) = temp;     
          }  




        /* state UTR3 is a source for special SPECIAL_UTR3 */ 
        temp = score + (mat->switchcost) + (0) ;     
        if( temp > GenomeWise9_VSMALL_SPECIAL(mat,i,j,SPECIAL_UTR3) )    {  
          GenomeWise9_VSMALL_SPECIAL(mat,i,j,SPECIAL_UTR3) = temp;   
          }  




        /* Finished calculating state UTR3 */ 


        /* For state UTR3_INTRON */ 
        /* setting first movement to score */ 
        score = GenomeWise9_VSMALL_MATRIX(mat,i-0,j-1,UTR3_INTRON) + GNE_UTR_INTRON(mat->evi,i,mat->gen,j);  
        /* From state UTR3 to state UTR3_INTRON */ 
        temp = GenomeWise9_VSMALL_MATRIX(mat,i-0,j-1,UTR3) + GNE_UTR_5SS(mat->evi,i,mat->gen,j);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SPECIAL_UTR3 to state UTR3_INTRON */ 
        temp = GenomeWise9_VSMALL_SPECIAL(mat,i-0,j-1,SPECIAL_UTR3) + GNE_UTR(mat->evi,i,mat->gen,j);    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for UTR3_INTRON */ 
        /* Add any movement independant score and put away */ 
         GenomeWise9_VSMALL_MATRIX(mat,i,j,UTR3_INTRON) = score; 


        /* Finished calculating state UTR3_INTRON */ 
        } /* end of for all query positions */ 




      /* Special state PREGENE_INTERGENIC has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GenomeWise9_VSMALL_SPECIAL(mat,0,j,PREGENE_INTERGENIC);    


      /* Source START is a special source for PREGENE_INTERGENIC */ 
      temp = GenomeWise9_VSMALL_SPECIAL(mat,0,j - 1,START) + (0) + (0);  
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      GenomeWise9_VSMALL_SPECIAL(mat,0,j,PREGENE_INTERGENIC) = score;    
      /* Finished updating state PREGENE_INTERGENIC */ 




      /* Special state POSTGENE_INTERGENIC has no special to special movements */ 


      /* Special state INTERGENIC has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GenomeWise9_VSMALL_SPECIAL(mat,0,j,INTERGENIC);    


      /* Source INTERGENIC is a special source for INTERGENIC */ 
      temp = GenomeWise9_VSMALL_SPECIAL(mat,0,j - 1,INTERGENIC) + (0) + (0);     
      if( temp > score ) 
        score = temp;    


      /* Source CDS for state INTERGENIC is not special... already calculated */ 
      /* Source UTR3 for state INTERGENIC is not special... already calculated */ 
      /* Put back score... (now updated!) */ 
      GenomeWise9_VSMALL_SPECIAL(mat,0,j,INTERGENIC) = score;    
      /* Finished updating state INTERGENIC */ 




      /* Special state SPECIAL_UTR5 has no special to special movements */ 


      /* Special state SPECIAL_UTR3 has no special to special movements */ 


      /* Special state SPECIAL_CDS has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GenomeWise9_VSMALL_SPECIAL(mat,0,j,SPECIAL_CDS);   


      /* Source CDS for state SPECIAL_CDS is not special... already calculated */ 
      /* Source SPECIAL_CDS is a special source for SPECIAL_CDS */ 
      temp = GenomeWise9_VSMALL_SPECIAL(mat,0,j - 3,SPECIAL_CDS) + (mat->rndcodon->codon[CSEQ_GENOMIC_CODON(mat->gen,j)]) + (0);     
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      GenomeWise9_VSMALL_SPECIAL(mat,0,j,SPECIAL_CDS) = score;   
      /* Finished updating state SPECIAL_CDS */ 




      /* Special state START has no special to special movements */ 


      /* Special state END has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GenomeWise9_VSMALL_SPECIAL(mat,0,j,END);   


      /* Source INTERGENIC is a special source for END */ 
      /* Has restricted position */ 
      if( j == mat->lenj-1 ) {  
        temp = GenomeWise9_VSMALL_SPECIAL(mat,0,j - 1,INTERGENIC) + (0) + (0);   
        if( temp > score )   
          score = temp;  
        }  


      /* Put back score... (now updated!) */ 
      GenomeWise9_VSMALL_SPECIAL(mat,0,j,END) = score;   
      /* Finished updating state END */ 


      if( bestscore < GenomeWise9_VSMALL_SPECIAL(mat,0,j,END) )  
        bestscore = GenomeWise9_VSMALL_SPECIAL(mat,0,j,END);     
      } /* end of for all target positions */ 


    mat = free_GenomeWise9(mat);     
    return bestscore;    
}    


/* Function:  PackAln_bestmemory_GenomeWise9(evi,gen,switchcost,newgenecost,non_start_codon,non_stop_codon,rndcodon,dpenv,dpri)
 *
 * Descrip:    This function chooses the best memory set-up for the alignment
 *             using calls to basematrix, and then implements either a large
 *             or small memory model.
 *
 *             It is the best function to use if you just want an alignment
 *
 *             If you want a label alignment, you will need
 *             /convert_PackAln_to_AlnBlock_GenomeWise9
 *
 *
 * Arg:                    evi [UNKN ] query data structure [GenomeEvidenceSet*]
 * Arg:                    gen [UNKN ] target data structure [ComplexSequence*]
 * Arg:             switchcost [UNKN ] Resource [int]
 * Arg:            newgenecost [UNKN ] Resource [int]
 * Arg:        non_start_codon [UNKN ] Resource [int]
 * Arg:         non_stop_codon [UNKN ] Resource [int]
 * Arg:               rndcodon [UNKN ] Resource [RandomCodonScore *]
 * Arg:                  dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:                   dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * PackAln_bestmemory_GenomeWise9(GenomeEvidenceSet* evi,ComplexSequence* gen ,int switchcost,int newgenecost,int non_start_codon,int non_stop_codon,RandomCodonScore * rndcodon,DPEnvelope * dpenv,DPRunImpl * dpri) 
{
    long long total; 
    GenomeWise9 * mat;   
    PackAln * out;   
    DebugMatrix * de;    
    DPRunImplMemory strategy;    
    assert(dpri);    


    total = evi->len * gen->seq->len;    
    if( dpri->memory == DPIM_Default )   {  
      if( (total * 10 * sizeof(int)) > 1000*dpri->kbyte_size)    {  
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
        if( (mat=allocate_Expl_GenomeWise9(evi, gen , switchcost, newgenecost, non_start_codon, non_stop_codon, rndcodon,dpri)) == NULL )    {  
          warn("Unable to allocate large GenomeWise9 version");  
          return NULL;   
          }  
        calculate_dpenv_GenomeWise9(mat,dpenv);  
        out =  PackAln_read_Expl_GenomeWise9(mat);   
        }  
      else   {  
        mat = allocate_GenomeWise9_only(evi, gen , switchcost, newgenecost, non_start_codon, non_stop_codon, rndcodon);  
        calculate_shatter_GenomeWise9(mat,dpenv);    
        out = PackAln_read_Shatter_GenomeWise9(mat);     
        }  
      }  
    else {  
      if( strategy == DPIM_Linear )  {  
        /* use small implementation */ 
        if( (mat=allocate_Small_GenomeWise9(evi, gen , switchcost, newgenecost, non_start_codon, non_stop_codon, rndcodon)) == NULL )    {  
          warn("Unable to allocate small GenomeWise9 version");  
          return NULL;   
          }  
        out = PackAln_calculate_Small_GenomeWise9(mat,dpenv);    
        }  
      else   {  
        /* use Large implementation */ 
        if( (mat=allocate_Expl_GenomeWise9(evi, gen , switchcost, newgenecost, non_start_codon, non_stop_codon, rndcodon,dpri)) == NULL )    {  
          warn("Unable to allocate large GenomeWise9 version");  
          return NULL;   
          }  
        if( dpri->debug == TRUE) {  
          fatal("Asked for dydebug, but dynamite file not compiled with -g. Need to recompile dynamite source"); 
          }  
        else {  
          calculate_GenomeWise9(mat);    
          out =  PackAln_read_Expl_GenomeWise9(mat); 
          }  
        }  
      }  


    mat = free_GenomeWise9(mat);     
    return out;  
}    


/* Function:  allocate_GenomeWise9_only(evi,gen,switchcost,newgenecost,non_start_codon,non_stop_codon,rndcodon)
 *
 * Descrip:    This function only allocates the GenomeWise9 structure
 *             checks types where possible and determines leni and lenj
 *             The basematrix area is delt with elsewhere
 *
 *
 * Arg:                    evi [UNKN ] query data structure [GenomeEvidenceSet*]
 * Arg:                    gen [UNKN ] target data structure [ComplexSequence*]
 * Arg:             switchcost [UNKN ] Resource [int]
 * Arg:            newgenecost [UNKN ] Resource [int]
 * Arg:        non_start_codon [UNKN ] Resource [int]
 * Arg:         non_stop_codon [UNKN ] Resource [int]
 * Arg:               rndcodon [UNKN ] Resource [RandomCodonScore *]
 *
 * Return [UNKN ]  Undocumented return value [GenomeWise9 *]
 *
 */
GenomeWise9 * allocate_GenomeWise9_only(GenomeEvidenceSet* evi,ComplexSequence* gen ,int switchcost,int newgenecost,int non_start_codon,int non_stop_codon,RandomCodonScore * rndcodon) 
{
    GenomeWise9 * out;   


    if((out= GenomeWise9_alloc()) == NULL)   {  
      warn("Allocation of basic GenomeWise9 structure failed...");   
      return NULL;   
      }  


    out->evi = evi;  
    out->gen = gen;  
    out->switchcost = switchcost;    
    out->newgenecost = newgenecost;  
    out->non_start_codon = non_start_codon;  
    out->non_stop_codon = non_stop_codon;    
    out->rndcodon = rndcodon;    
    out->leni = evi->len;    
    out->lenj = gen->seq->len;   
    return out;  
}    


/* Function:  allocate_Expl_GenomeWise9(evi,gen,switchcost,newgenecost,non_start_codon,non_stop_codon,rndcodon,dpri)
 *
 * Descrip:    This function allocates the GenomeWise9 structure
 *             and the basematrix area for explicit memory implementations
 *             It calls /allocate_GenomeWise9_only
 *
 *
 * Arg:                    evi [UNKN ] query data structure [GenomeEvidenceSet*]
 * Arg:                    gen [UNKN ] target data structure [ComplexSequence*]
 * Arg:             switchcost [UNKN ] Resource [int]
 * Arg:            newgenecost [UNKN ] Resource [int]
 * Arg:        non_start_codon [UNKN ] Resource [int]
 * Arg:         non_stop_codon [UNKN ] Resource [int]
 * Arg:               rndcodon [UNKN ] Resource [RandomCodonScore *]
 * Arg:                   dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [GenomeWise9 *]
 *
 */
GenomeWise9 * allocate_Expl_GenomeWise9(GenomeEvidenceSet* evi,ComplexSequence* gen ,int switchcost,int newgenecost,int non_start_codon,int non_stop_codon,RandomCodonScore * rndcodon,DPRunImpl * dpri) 
{
    GenomeWise9 * out;   


    out = allocate_GenomeWise9_only(evi, gen , switchcost, newgenecost, non_start_codon, non_stop_codon, rndcodon);  
    if( out == NULL )    
      return NULL;   
    if( dpri->should_cache == TRUE ) {  
      if( dpri->cache != NULL )  {  
        if( dpri->cache->maxleni >= (out->lenj+10)*10 && dpri->cache->maxlenj >= (out->leni+0))  
          out->basematrix = hard_link_BaseMatrix(dpri->cache);   
        else 
          dpri->cache = free_BaseMatrix(dpri->cache);    
        }  
      }  
    if( out->basematrix == NULL )    {  
      if( (out->basematrix = BaseMatrix_alloc_matrix_and_specials((out->lenj+10)*10,(out->leni+0),8,out->lenj+10)) == NULL)  {  
        warn("Explicit matrix GenomeWise9 cannot be allocated, (asking for %d by %d main cells)",out->leni,out->lenj);   
        free_GenomeWise9(out);   
        return NULL; 
        }  
      }  
    if( dpri->should_cache == TRUE && dpri->cache == NULL)   
      dpri->cache = hard_link_BaseMatrix(out->basematrix);   
    out->basematrix->type = BASEMATRIX_TYPE_EXPLICIT;    
    init_GenomeWise9(out);   
    return out;  
}    


/* Function:  init_GenomeWise9(mat)
 *
 * Descrip:    This function initates GenomeWise9 matrix when in explicit mode
 *             Called in /allocate_Expl_GenomeWise9
 *
 *
 * Arg:        mat [UNKN ] GenomeWise9 which contains explicit basematrix memory [GenomeWise9 *]
 *
 */
void init_GenomeWise9(GenomeWise9 * mat) 
{
    register int i;  
    register int j;  
    if( mat->basematrix->type != BASEMATRIX_TYPE_EXPLICIT)   {  
      warn("Cannot iniate matrix, is not an explicit memory type and you have assummed that");   
      return;    
      }  


    for(i= (-0);i<mat->evi->len;i++) {  
      for(j= (-10);j<11;j++) {  
        GenomeWise9_EXPL_MATRIX(mat,i,j,UTR5) = NEGI;    
        GenomeWise9_EXPL_MATRIX(mat,i,j,UTR5_INTRON) = NEGI; 
        GenomeWise9_EXPL_MATRIX(mat,i,j,START_CODON) = NEGI; 
        GenomeWise9_EXPL_MATRIX(mat,i,j,CDS) = NEGI; 
        GenomeWise9_EXPL_MATRIX(mat,i,j,CDS_INTRON_0) = NEGI;    
        GenomeWise9_EXPL_MATRIX(mat,i,j,CDS_INTRON_1) = NEGI;    
        GenomeWise9_EXPL_MATRIX(mat,i,j,CDS_INTRON_2) = NEGI;    
        GenomeWise9_EXPL_MATRIX(mat,i,j,STOP_CODON) = NEGI;  
        GenomeWise9_EXPL_MATRIX(mat,i,j,UTR3) = NEGI;    
        GenomeWise9_EXPL_MATRIX(mat,i,j,UTR3_INTRON) = NEGI; 
        }  
      }  
    for(j= (-10);j<mat->gen->seq->len;j++)   {  
      for(i= (-0);i<1;i++)   {  
        GenomeWise9_EXPL_MATRIX(mat,i,j,UTR5) = NEGI;    
        GenomeWise9_EXPL_MATRIX(mat,i,j,UTR5_INTRON) = NEGI; 
        GenomeWise9_EXPL_MATRIX(mat,i,j,START_CODON) = NEGI; 
        GenomeWise9_EXPL_MATRIX(mat,i,j,CDS) = NEGI; 
        GenomeWise9_EXPL_MATRIX(mat,i,j,CDS_INTRON_0) = NEGI;    
        GenomeWise9_EXPL_MATRIX(mat,i,j,CDS_INTRON_1) = NEGI;    
        GenomeWise9_EXPL_MATRIX(mat,i,j,CDS_INTRON_2) = NEGI;    
        GenomeWise9_EXPL_MATRIX(mat,i,j,STOP_CODON) = NEGI;  
        GenomeWise9_EXPL_MATRIX(mat,i,j,UTR3) = NEGI;    
        GenomeWise9_EXPL_MATRIX(mat,i,j,UTR3_INTRON) = NEGI; 
        }  
      GenomeWise9_EXPL_SPECIAL(mat,i,j,PREGENE_INTERGENIC) = NEGI;   
      GenomeWise9_EXPL_SPECIAL(mat,i,j,POSTGENE_INTERGENIC) = NEGI;  
      GenomeWise9_EXPL_SPECIAL(mat,i,j,INTERGENIC) = NEGI;   
      GenomeWise9_EXPL_SPECIAL(mat,i,j,SPECIAL_UTR5) = NEGI; 
      GenomeWise9_EXPL_SPECIAL(mat,i,j,SPECIAL_UTR3) = NEGI; 
      GenomeWise9_EXPL_SPECIAL(mat,i,j,SPECIAL_CDS) = NEGI;  
      GenomeWise9_EXPL_SPECIAL(mat,i,j,START) = 0;   
      GenomeWise9_EXPL_SPECIAL(mat,i,j,END) = NEGI;  
      }  
    return;  
}    


/* Function:  recalculate_PackAln_GenomeWise9(pal,mat)
 *
 * Descrip:    This function recalculates the PackAln structure produced by GenomeWise9
 *             For example, in linear space methods this is used to score them
 *
 *
 * Arg:        pal [UNKN ] Undocumented argument [PackAln *]
 * Arg:        mat [UNKN ] Undocumented argument [GenomeWise9 *]
 *
 */
void recalculate_PackAln_GenomeWise9(PackAln * pal,GenomeWise9 * mat) 
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
        case UTR5 :  
          if( offi == 0 && offj == 1 && prev->state == UTR5 )    {  
            pau->score = GNE_UTR(mat->evi,i,mat->gen,j) + (0);   
            continue;    
            }  
          if( offi == 0 && offj == 1 && prev->state == UTR5_INTRON ) {  
            pau->score = GNE_UTR_3SS(mat->evi,i,mat->gen,j) + (0);   
            continue;    
            }  
          if( offj == 1 && prev->state == (SPECIAL_UTR5+10) )    {  
            pau->score = GNE_UTR(mat->evi,i,mat->gen,j) + (0);   
            continue;    
            }  
          if( offj == 1 && prev->state == (START+10) )   {  
            pau->score = GNE_UTR5_START(mat->evi,i,mat->gen,j) + (0);    
            continue;    
            }  
          if( offj == 1 && prev->state == (INTERGENIC+10) )  {  
            pau->score = (GNE_UTR(mat->evi,i,mat->gen,j)+GNE_UTR5_START(mat->evi,i,mat->gen,j)) + (0);   
            continue;    
            }  
          if( offj == 1 && prev->state == (PREGENE_INTERGENIC+10) )  {  
            pau->score = GNE_UTR(mat->evi,i,mat->gen,j) + (0);   
            continue;    
            }  
          warn("In recaluclating PackAln with state UTR5, from [%d,%d,%d], got a bad source state. Error!",offi,offj,prev->state);   
          break; 
        case UTR5_INTRON :   
          if( offi == 0 && offj == 1 && prev->state == UTR5_INTRON ) {  
            pau->score = GNE_UTR_INTRON(mat->evi,i,mat->gen,j) + (0);    
            continue;    
            }  
          if( offi == 0 && offj == 1 && prev->state == UTR5 )    {  
            pau->score = GNE_UTR_5SS(mat->evi,i,mat->gen,j) + (0);   
            continue;    
            }  
          warn("In recaluclating PackAln with state UTR5_INTRON, from [%d,%d,%d], got a bad source state. Error!",offi,offj,prev->state);    
          break; 
        case START_CODON :   
          if( offi == 0 && offj == 3 && prev->state == UTR5_INTRON ) {  
            pau->score = GNE_START_CODON(mat->evi,i,mat->gen,j) + (0);   
            continue;    
            }  
          if( offi == 0 && offj == 3 && prev->state == UTR5 )    {  
            pau->score = GNE_START_CODON(mat->evi,i,mat->gen,j) + (0);   
            continue;    
            }  
          warn("In recaluclating PackAln with state START_CODON, from [%d,%d,%d], got a bad source state. Error!",offi,offj,prev->state);    
          break; 
        case CDS :   
          if( offi == 0 && offj == 3 && prev->state == CDS ) {  
            pau->score = GNE_CDS(mat->evi,i,mat->gen,j) + (0);   
            continue;    
            }  
          if( offi == 0 && offj == 6 && prev->state == CDS_INTRON_0 )    {  
            pau->score = GNE_CDS_3SS(mat->evi,i,mat->gen,(j-3),0) + (0);     
            continue;    
            }  
          if( offi == 0 && offj == 5 && prev->state == CDS_INTRON_1 )    {  
            pau->score = GNE_CDS_3SS(mat->evi,i,mat->gen,(j-2),1) + (0);     
            continue;    
            }  
          if( offi == 0 && offj == 4 && prev->state == CDS_INTRON_2 )    {  
            pau->score = GNE_CDS_3SS(mat->evi,i,mat->gen,(j-1),2) + (0);     
            continue;    
            }  
          if( offi == 0 && offj == 2 && prev->state == CDS ) {  
            pau->score = GNE_CDS_FRAMESHIFT(mat->evi,i,mat->gen,j,2) + (0);  
            continue;    
            }  
          if( offi == 0 && offj == 4 && prev->state == CDS ) {  
            pau->score = GNE_CDS_FRAMESHIFT(mat->evi,i,mat->gen,j,4) + (0);  
            continue;    
            }  
          if( offi == 0 && offj == 3 && prev->state == UTR5 )    {  
            pau->score = (GNE_CDS(mat->evi,i,mat->gen,j)+mat->non_start_codon) + (0);    
            continue;    
            }  
          if( offi == 0 && offj == 3 && prev->state == START_CODON ) {  
            pau->score = GNE_CDS(mat->evi,i,mat->gen,j) + (0);   
            continue;    
            }  
          if( offj == 3 && prev->state == (SPECIAL_CDS+10) ) {  
            pau->score = GNE_CDS(mat->evi,i,mat->gen,j) + (0);   
            continue;    
            }  
          if( offj == 3 && prev->state == (INTERGENIC+10) )  {  
            pau->score = (GNE_CDS(mat->evi,i,mat->gen,j)+GNE_UTR5_START(mat->evi,i,mat->gen,j)) + (0);   
            continue;    
            }  
          warn("In recaluclating PackAln with state CDS, from [%d,%d,%d], got a bad source state. Error!",offi,offj,prev->state);    
          break; 
        case CDS_INTRON_0 :  
          if( offi == 0 && offj == 1 && prev->state == CDS_INTRON_0 )    {  
            pau->score = GNE_CDS_INTRON(mat->evi,i,mat->gen,j) + (0);    
            continue;    
            }  
          if( offi == 0 && offj == 8 && prev->state == CDS ) {  
            pau->score = GNE_CDS_5SS(mat->evi,i,mat->gen,(j-7),0) + (0);     
            continue;    
            }  
          warn("In recaluclating PackAln with state CDS_INTRON_0, from [%d,%d,%d], got a bad source state. Error!",offi,offj,prev->state);   
          break; 
        case CDS_INTRON_1 :  
          if( offi == 0 && offj == 1 && prev->state == CDS_INTRON_1 )    {  
            pau->score = GNE_CDS_INTRON(mat->evi,i,mat->gen,j) + (0);    
            continue;    
            }  
          if( offi == 0 && offj == 9 && prev->state == CDS ) {  
            pau->score = GNE_CDS_5SS(mat->evi,i,mat->gen,(j-7),1) + (0);     
            continue;    
            }  
          warn("In recaluclating PackAln with state CDS_INTRON_1, from [%d,%d,%d], got a bad source state. Error!",offi,offj,prev->state);   
          break; 
        case CDS_INTRON_2 :  
          if( offi == 0 && offj == 1 && prev->state == CDS_INTRON_2 )    {  
            pau->score = GNE_CDS_INTRON(mat->evi,i,mat->gen,j) + (0);    
            continue;    
            }  
          if( offi == 0 && offj == 10 && prev->state == CDS )    {  
            pau->score = GNE_CDS_5SS(mat->evi,i,mat->gen,(j-7),2) + (0);     
            continue;    
            }  
          warn("In recaluclating PackAln with state CDS_INTRON_2, from [%d,%d,%d], got a bad source state. Error!",offi,offj,prev->state);   
          break; 
        case STOP_CODON :    
          if( offi == 0 && offj == 3 && prev->state == CDS ) {  
            pau->score = GNE_STOP_CODON(mat->evi,i,mat->gen,j) + (0);    
            continue;    
            }  
          warn("In recaluclating PackAln with state STOP_CODON, from [%d,%d,%d], got a bad source state. Error!",offi,offj,prev->state); 
          break; 
        case UTR3 :  
          if( offi == 0 && offj == 1 && prev->state == UTR3 )    {  
            pau->score = GNE_UTR(mat->evi,i,mat->gen,j) + (0);   
            continue;    
            }  
          if( offi == 0 && offj == 1 && prev->state == CDS ) {  
            pau->score = (GNE_UTR(mat->evi,i,mat->gen,j)+mat->non_stop_codon) + (0);     
            continue;    
            }  
          if( offi == 0 && offj == 1 && prev->state == STOP_CODON )  {  
            pau->score = GNE_UTR(mat->evi,i,mat->gen,j) + (0);   
            continue;    
            }  
          if( offi == 0 && offj == 1 && prev->state == UTR3_INTRON ) {  
            pau->score = GNE_UTR_3SS(mat->evi,i,mat->gen,j) + (0);   
            continue;    
            }  
          if( offj == 1 && prev->state == (INTERGENIC+10) )  {  
            pau->score = GNE_UTR(mat->evi,i,mat->gen,j) + (0);   
            continue;    
            }  
          warn("In recaluclating PackAln with state UTR3, from [%d,%d,%d], got a bad source state. Error!",offi,offj,prev->state);   
          break; 
        case UTR3_INTRON :   
          if( offi == 0 && offj == 1 && prev->state == UTR3_INTRON ) {  
            pau->score = GNE_UTR_INTRON(mat->evi,i,mat->gen,j) + (0);    
            continue;    
            }  
          if( offi == 0 && offj == 1 && prev->state == UTR3 )    {  
            pau->score = GNE_UTR_5SS(mat->evi,i,mat->gen,j) + (0);   
            continue;    
            }  
          if( offj == 1 && prev->state == (SPECIAL_UTR3+10) )    {  
            pau->score = GNE_UTR(mat->evi,i,mat->gen,j) + (0);   
            continue;    
            }  
          warn("In recaluclating PackAln with state UTR3_INTRON, from [%d,%d,%d], got a bad source state. Error!",offi,offj,prev->state);    
          break; 
        case (PREGENE_INTERGENIC+10) :   
          if( offj == 1 && prev->state == (START+10) )   {  
            pau->score = 0 + (0);    
            continue;    
            }  
          warn("In recaluclating PackAln with state PREGENE_INTERGENIC, got a bad source state. Error!");    
          break; 
        case (POSTGENE_INTERGENIC+10) :  
          if( offj == 0 && prev->state == UTR3 ) {  
            /* i here comes from the previous state ;) - not the real one */ 
            i = prev->i; 
            pau->score = (mat->newgenecost+GNE_UTR3_END(mat->evi,i,mat->gen,j)) + (0);   
            continue;    
            }  
          warn("In recaluclating PackAln with state POSTGENE_INTERGENIC, got a bad source state. Error!");   
          break; 
        case (INTERGENIC+10) :   
          if( offj == 1 && prev->state == (INTERGENIC+10) )  {  
            pau->score = 0 + (0);    
            continue;    
            }  
          if( offj == 0 && prev->state == CDS )  {  
            /* i here comes from the previous state ;) - not the real one */ 
            i = prev->i; 
            pau->score = (mat->newgenecost+GNE_UTR3_END(mat->evi,i,mat->gen,j)) + (0);   
            continue;    
            }  
          if( offj == 0 && prev->state == UTR3 ) {  
            /* i here comes from the previous state ;) - not the real one */ 
            i = prev->i; 
            pau->score = (mat->newgenecost+GNE_UTR3_END(mat->evi,i,mat->gen,j)) + (0);   
            continue;    
            }  
          warn("In recaluclating PackAln with state INTERGENIC, got a bad source state. Error!");    
          break; 
        case (SPECIAL_UTR5+10) :     
          if( offj == 0 && prev->state == UTR5 ) {  
            /* i here comes from the previous state ;) - not the real one */ 
            i = prev->i; 
            pau->score = mat->switchcost + (0);  
            continue;    
            }  
          warn("In recaluclating PackAln with state SPECIAL_UTR5, got a bad source state. Error!");  
          break; 
        case (SPECIAL_UTR3+10) :     
          if( offj == 0 && prev->state == UTR3 ) {  
            /* i here comes from the previous state ;) - not the real one */ 
            i = prev->i; 
            pau->score = mat->switchcost + (0);  
            continue;    
            }  
          warn("In recaluclating PackAln with state SPECIAL_UTR3, got a bad source state. Error!");  
          break; 
        case (SPECIAL_CDS+10) :  
          if( offj == 0 && prev->state == CDS )  {  
            /* i here comes from the previous state ;) - not the real one */ 
            i = prev->i; 
            pau->score = (mat->switchcost+mat->rndcodon->codon[CSEQ_GENOMIC_CODON(mat->gen,j)]) + (0);   
            continue;    
            }  
          if( offj == 3 && prev->state == (SPECIAL_CDS+10) ) {  
            pau->score = mat->rndcodon->codon[CSEQ_GENOMIC_CODON(mat->gen,j)] + (0);     
            continue;    
            }  
          warn("In recaluclating PackAln with state SPECIAL_CDS, got a bad source state. Error!");   
          break; 
        case (START+10) :    
          warn("In recaluclating PackAln with state START, got a bad source state. Error!"); 
          break; 
        case (END+10) :  
          if( offj == 1 && prev->state == (INTERGENIC+10) )  {  
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
#define GenomeWise9_HIDDEN_MATRIX(thismatrix,i,j,state) (thismatrix->basematrix->matrix[(j-hiddenj+10)][(i+0)*10+state]) 
#define GenomeWise9_DC_SHADOW_MATRIX(thismatrix,i,j,state) (thismatrix->basematrix->matrix[((j+11)*8) % 88][(i+0)*10+state]) 
#define GenomeWise9_HIDDEN_SPECIAL(thismatrix,i,j,state) (thismatrix->basematrix->specmatrix[state][(j+10)]) 
#define GenomeWise9_DC_SHADOW_SPECIAL(thismatrix,i,j,state) (thismatrix->basematrix->specmatrix[state*8][(j+10)])    
#define GenomeWise9_DC_SHADOW_MATRIX_SP(thismatrix,i,j,state,shadow) (thismatrix->basematrix->matrix[((((j+11)*8)+(shadow+1)) % 88)][(i+0)*10 + state])  
#define GenomeWise9_DC_SHADOW_SPECIAL_SP(thismatrix,i,j,state,shadow) (thismatrix->basematrix->specmatrix[state*8 +shadow+1][(j+10)])    
#define GenomeWise9_DC_OPT_SHADOW_MATRIX(thismatrix,i,j,state) (score_pointers[(((j+10)% 10) * (leni+1) * 10) + ((i+0) * 10) + (state)]) 
#define GenomeWise9_DC_OPT_SHADOW_MATRIX_SP(thismatrix,i,j,state,shadow) (shadow_pointers[(((j+10)% 10) * (leni+1) * 80) + ((i+0) * 80) + (state * 8) + shadow+1])   
#define GenomeWise9_DC_OPT_SHADOW_SPECIAL(thismatrix,i,j,state) (thismatrix->basematrix->specmatrix[state*8][(j+10)])    
/* Function:  allocate_Small_GenomeWise9(evi,gen,switchcost,newgenecost,non_start_codon,non_stop_codon,rndcodon)
 *
 * Descrip:    This function allocates the GenomeWise9 structure
 *             and the basematrix area for a small memory implementations
 *             It calls /allocate_GenomeWise9_only
 *
 *
 * Arg:                    evi [UNKN ] query data structure [GenomeEvidenceSet*]
 * Arg:                    gen [UNKN ] target data structure [ComplexSequence*]
 * Arg:             switchcost [UNKN ] Resource [int]
 * Arg:            newgenecost [UNKN ] Resource [int]
 * Arg:        non_start_codon [UNKN ] Resource [int]
 * Arg:         non_stop_codon [UNKN ] Resource [int]
 * Arg:               rndcodon [UNKN ] Resource [RandomCodonScore *]
 *
 * Return [UNKN ]  Undocumented return value [GenomeWise9 *]
 *
 */
#define GenomeWise9_DC_OPT_SHADOW_SPECIAL_SP(thismatrix,i,j,state,shadow) (thismatrix->basematrix->specmatrix[state*8 +shadow+1][(j+10)])    
GenomeWise9 * allocate_Small_GenomeWise9(GenomeEvidenceSet* evi,ComplexSequence* gen ,int switchcost,int newgenecost,int non_start_codon,int non_stop_codon,RandomCodonScore * rndcodon) 
{
    GenomeWise9 * out;   


    out = allocate_GenomeWise9_only(evi, gen , switchcost, newgenecost, non_start_codon, non_stop_codon, rndcodon);  
    if( out == NULL )    
      return NULL;   
    out->basematrix = BaseMatrix_alloc_matrix_and_specials(88,(out->leni + 0) * 10,64,out->lenj+10);     
    if(out == NULL)  {  
      warn("Small shadow matrix GenomeWise9 cannot be allocated, (asking for 11 by %d main cells)",out->leni+1); 
      free_GenomeWise9(out);     
      return NULL;   
      }  
    out->basematrix->type = BASEMATRIX_TYPE_SHADOW;  
    return out;  
}    


/* Function:  PackAln_calculate_Small_GenomeWise9(mat,dpenv)
 *
 * Descrip:    This function calculates an alignment for GenomeWise9 structure in linear space
 *             If you want only the start/end points
 *             use /AlnRangeSet_calculate_Small_GenomeWise9 
 *
 *             The function basically
 *               finds start/end points 
 *               foreach start/end point 
 *                 calls /full_dc_GenomeWise9 
 *
 *
 * Arg:          mat [UNKN ] Undocumented argument [GenomeWise9 *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * PackAln_calculate_Small_GenomeWise9(GenomeWise9 * mat,DPEnvelope * dpenv) 
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
      warn("Could not calculate packaln small for GenomeWise9 due to wrong type of matrix"); 
      return NULL;   
      }  


    out = PackAln_alloc_std();   


    start_reporting("Find start end points: ");  
    dc_optimised_start_end_calc_GenomeWise9(mat,dpenv);  
    score = start_end_find_end_GenomeWise9(mat,&endj);   
    out->score = score;  
    stopstate = END;
    
    /* Special to specials: have to eat up in strip and then drop back to full_dc for intervening bits */ 
    log_full_error(REPORT,0,"End at %d Score %d",endj,score);    
    stop_reporting();    
    for(;;)  { /*while there are more special bits to recover*/ 
      start_reporting("Special cell aln end   %d:",endj);    
      if( read_special_strip_GenomeWise9(mat,0,endj,stopstate,&endj,&startstate,out) == FALSE )  {  
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
      starti = GenomeWise9_DC_SHADOW_SPECIAL_SP(mat,0,endj,temp,0);  
      startj = GenomeWise9_DC_SHADOW_SPECIAL_SP(mat,0,endj,temp,1);  
      startstate = GenomeWise9_DC_SHADOW_SPECIAL_SP(mat,0,endj,temp,2);  
      stopi = GenomeWise9_DC_SHADOW_SPECIAL_SP(mat,0,endj,temp,3);   
      stopj = GenomeWise9_DC_SHADOW_SPECIAL_SP(mat,0,endj,temp,4);   
      stopstate = GenomeWise9_DC_SHADOW_SPECIAL_SP(mat,0,endj,temp,5);   


      /* Get out the score of this block. V. important! */ 
      temp = GenomeWise9_DC_SHADOW_SPECIAL_SP(mat,0,endj,temp,6);    
      totalj = stopj - startj;   
      donej  = 0;    
      start_reporting("Main matrix  aln [%d,%d]:",startj,stopj);     
      if(full_dc_GenomeWise9(mat,starti,startj,startstate,stopi,stopj,stopstate,out,&donej,totalj,dpenv) == FALSE)   {  
        warn("In the alignment GenomeWise9 [%d,%d][%d,%d], got a problem. Please report bug ... giving you back a partial alignment",starti,startj,stopi,stopj); 
        return out;  
        }  


      /* now have to figure out which special we came from... yikes */ 
      max_matrix_to_special_GenomeWise9(mat,starti,startj,startstate,temp,&stopi,&stopj,&stopstate,&temp,NULL);  
      if( stopi == GenomeWise9_READ_OFF_ERROR)   {  
        warn("In GenomeWise9 read off ending at %d ... got a bad matrix to special read off... returning partial alignment",startj); 
        invert_PackAln(out); 
        recalculate_PackAln_GenomeWise9(out,mat);    
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
    recalculate_PackAln_GenomeWise9(out,mat);    
    return out;  


}    


/* Function:  AlnRangeSet_calculate_Small_GenomeWise9(mat)
 *
 * Descrip:    This function calculates an alignment for GenomeWise9 structure in linear space
 *             If you want the full alignment, use /PackAln_calculate_Small_GenomeWise9 
 *             If you have already got the full alignment, but want the range set, use /AlnRangeSet_from_PackAln_GenomeWise9
 *             If you have got the small matrix but not the alignment, use /AlnRangeSet_from_GenomeWise9 
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [GenomeWise9 *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * AlnRangeSet_calculate_Small_GenomeWise9(GenomeWise9 * mat) 
{
    AlnRangeSet * out;   


    start_reporting("Find start end points: ");  
    dc_optimised_start_end_calc_GenomeWise9(mat,NULL);   
    log_full_error(REPORT,0,"Calculated");   


    out = AlnRangeSet_from_GenomeWise9(mat); 
    return out;  
}    


/* Function:  AlnRangeSet_from_GenomeWise9(mat)
 *
 * Descrip:    This function reads off a start/end structure
 *             for GenomeWise9 structure in linear space
 *             If you want the full alignment use
 *             /PackAln_calculate_Small_GenomeWise9 
 *             If you have not calculated the matrix use
 *             /AlnRange_calculate_Small_GenomeWise9
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [GenomeWise9 *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * AlnRangeSet_from_GenomeWise9(GenomeWise9 * mat) 
{
    AlnRangeSet * out;   
    AlnRange * temp; 
    int jpos;    
    int state;   


    if( mat->basematrix->type != BASEMATRIX_TYPE_SHADOW) {  
      warn("Bad error! - non shadow matrix type in AlnRangeSet_from_GenomeWise9");   
      return NULL;   
      }  


    out = AlnRangeSet_alloc_std();   
    /* Find the end position */ 
    out->score = start_end_find_end_GenomeWise9(mat,&jpos);  
    state = END; 


    while( (temp = AlnRange_build_GenomeWise9(mat,jpos,state,&jpos,&state)) != NULL) 
      add_AlnRangeSet(out,temp); 
    return out;  
}    


/* Function:  AlnRange_build_GenomeWise9(mat,stopj,stopspecstate,startj,startspecstate)
 *
 * Descrip:    This function calculates a single start/end set in linear space
 *             Really a sub-routine for /AlnRangeSet_from_PackAln_GenomeWise9
 *
 *
 * Arg:                   mat [UNKN ] Undocumented argument [GenomeWise9 *]
 * Arg:                 stopj [UNKN ] Undocumented argument [int]
 * Arg:         stopspecstate [UNKN ] Undocumented argument [int]
 * Arg:                startj [UNKN ] Undocumented argument [int *]
 * Arg:        startspecstate [UNKN ] Undocumented argument [int *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRange *]
 *
 */
AlnRange * AlnRange_build_GenomeWise9(GenomeWise9 * mat,int stopj,int stopspecstate,int * startj,int * startspecstate) 
{
    AlnRange * out;  
    int jpos;    
    int state;   


    if( mat->basematrix->type != BASEMATRIX_TYPE_SHADOW) {  
      warn("Bad error! - non shadow matrix type in AlnRangeSet_from_GenomeWise9");   
      return NULL;   
      }  


    /* Assumme that we have specials (we should!). Read back along the specials till we have the finish point */ 
    if( read_special_strip_GenomeWise9(mat,0,stopj,stopspecstate,&jpos,&state,NULL) == FALSE)    {  
      warn("In AlnRanger_build_GenomeWise9 alignment ending at %d, unable to read back specials. Will (evenutally) return a partial range set... BEWARE!",stopj);    
      return NULL;   
      }  
    if( state == START || jpos <= 0) 
      return NULL;   


    out = AlnRange_alloc();  


    out->starti = GenomeWise9_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,0);  
    out->startj = GenomeWise9_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,1);  
    out->startstate = GenomeWise9_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,2);  
    out->stopi = GenomeWise9_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,3);   
    out->stopj = GenomeWise9_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,4);   
    out->stopstate = GenomeWise9_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,5);   
    out->startscore = GenomeWise9_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,6);  
    out->stopscore = GenomeWise9_DC_SHADOW_SPECIAL(mat,0,jpos,state);    


    /* Now, we have to figure out where this state came from in the specials */ 
    max_matrix_to_special_GenomeWise9(mat,out->starti,out->startj,out->startstate,out->startscore,&jpos,startj,startspecstate,&state,NULL);  
    if( jpos == GenomeWise9_READ_OFF_ERROR)  {  
      warn("In AlnRange_build_GenomeWise9 alignment ending at %d, with aln range between %d-%d in j, unable to find source special, returning this range, but this could get tricky!",stopj,out->startj,out->stopj); 
      return out;    
      }  


    /* Put in the correct score for startstate, from the special */ 
    out->startscore = GenomeWise9_DC_SHADOW_SPECIAL(mat,0,*startj,*startspecstate);  
    /* The correct j coords have been put into startj, startspecstate... so just return out */ 
    return out;  
}    


/* Function:  read_hidden_GenomeWise9(mat,starti,startj,startstate,stopi,stopj,stopstate,out)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [GenomeWise9 *]
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
boolean read_hidden_GenomeWise9(GenomeWise9 * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out) 
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


      max_hidden_GenomeWise9(mat,startj,i,j,state,isspecial,&i,&j,&state,&isspecial,&cellscore); 


      if( i == GenomeWise9_READ_OFF_ERROR)   {  
        warn("In GenomeWise9 hidden read off, between %d:%d,%d:%d - at got bad read off. Problem!",starti,startj,stopi,stopj);   
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
        warn("In GenomeWise9 hidden read off, between %d:%d,%d:%d - hit start cell, but not in start state. Can't be good!.",starti,startj,stopi,stopj); 
        return FALSE;    
        }  
      }  
    warn("In GenomeWise9 hidden read off, between %d:%d,%d:%d - gone past start cell (now in %d,%d,%d), can't be good news!.",starti,startj,stopi,stopj,i,j,state);  
    return FALSE;    
}    


/* Function:  max_hidden_GenomeWise9(mat,hiddenj,i,j,state,isspecial,reti,retj,retstate,retspecial,cellscore)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [GenomeWise9 *]
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
int max_hidden_GenomeWise9(GenomeWise9 * mat,int hiddenj,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore) 
{
    register int temp;   
    register int cscore; 


    *reti = (*retj) = (*retstate) = GenomeWise9_READ_OFF_ERROR;  


    if( i < 0 || j < 0 || i > mat->evi->len || j > mat->gen->seq->len)   {  
      warn("In GenomeWise9 matrix special read off - out of bounds on matrix [i,j is %d,%d state %d in standard matrix]",i,j,state); 
      return -1; 
      }  


    /* Then you have to select the correct switch statement to figure out the readoff      */ 
    /* Somewhat odd - reverse the order of calculation and return as soon as it is correct */ 
    cscore = GenomeWise9_HIDDEN_MATRIX(mat,i,j,state);   
    switch(state)    { /*Switch state */ 
      case UTR5 :    
        /* Not allowing special sources.. skipping PREGENE_INTERGENIC */ 
        /* Not allowing special sources.. skipping INTERGENIC */ 
        /* Not allowing special sources.. skipping START */ 
        /* Not allowing special sources.. skipping SPECIAL_UTR5 */ 
        temp = cscore - (GNE_UTR_3SS(mat->evi,i,mat->gen,j)) -  (0); 
        if( temp == GenomeWise9_HIDDEN_MATRIX(mat,i - 0,j - 1,UTR5_INTRON) ) {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = UTR5_INTRON;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GenomeWise9_HIDDEN_MATRIX(mat,i-0,j-1,UTR5_INTRON);    
            }  
          return GenomeWise9_HIDDEN_MATRIX(mat,i - 0,j - 1,UTR5_INTRON);     
          }  
        temp = cscore - (GNE_UTR(mat->evi,i,mat->gen,j)) -  (0); 
        if( temp == GenomeWise9_HIDDEN_MATRIX(mat,i - 0,j - 1,UTR5) )    {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = UTR5;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GenomeWise9_HIDDEN_MATRIX(mat,i-0,j-1,UTR5);   
            }  
          return GenomeWise9_HIDDEN_MATRIX(mat,i - 0,j - 1,UTR5);    
          }  
        warn("Major problem (!) - in GenomeWise9 read off, position %d,%d state %d no source found!",i,j,state); 
        return (-1); 
      case UTR5_INTRON :     
        temp = cscore - (GNE_UTR_5SS(mat->evi,i,mat->gen,j)) -  (0); 
        if( temp == GenomeWise9_HIDDEN_MATRIX(mat,i - 0,j - 1,UTR5) )    {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = UTR5;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GenomeWise9_HIDDEN_MATRIX(mat,i-0,j-1,UTR5);   
            }  
          return GenomeWise9_HIDDEN_MATRIX(mat,i - 0,j - 1,UTR5);    
          }  
        temp = cscore - (GNE_UTR_INTRON(mat->evi,i,mat->gen,j)) -  (0);  
        if( temp == GenomeWise9_HIDDEN_MATRIX(mat,i - 0,j - 1,UTR5_INTRON) ) {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = UTR5_INTRON;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GenomeWise9_HIDDEN_MATRIX(mat,i-0,j-1,UTR5_INTRON);    
            }  
          return GenomeWise9_HIDDEN_MATRIX(mat,i - 0,j - 1,UTR5_INTRON);     
          }  
        warn("Major problem (!) - in GenomeWise9 read off, position %d,%d state %d no source found!",i,j,state); 
        return (-1); 
      case START_CODON :     
        temp = cscore - (GNE_START_CODON(mat->evi,i,mat->gen,j)) -  (0); 
        if( temp == GenomeWise9_HIDDEN_MATRIX(mat,i - 0,j - 3,UTR5) )    {  
          *reti = i - 0; 
          *retj = j - 3; 
          *retstate = UTR5;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GenomeWise9_HIDDEN_MATRIX(mat,i-0,j-3,UTR5);   
            }  
          return GenomeWise9_HIDDEN_MATRIX(mat,i - 0,j - 3,UTR5);    
          }  
        temp = cscore - (GNE_START_CODON(mat->evi,i,mat->gen,j)) -  (0); 
        if( temp == GenomeWise9_HIDDEN_MATRIX(mat,i - 0,j - 3,UTR5_INTRON) ) {  
          *reti = i - 0; 
          *retj = j - 3; 
          *retstate = UTR5_INTRON;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GenomeWise9_HIDDEN_MATRIX(mat,i-0,j-3,UTR5_INTRON);    
            }  
          return GenomeWise9_HIDDEN_MATRIX(mat,i - 0,j - 3,UTR5_INTRON);     
          }  
        warn("Major problem (!) - in GenomeWise9 read off, position %d,%d state %d no source found!",i,j,state); 
        return (-1); 
      case CDS :     
        /* Not allowing special sources.. skipping INTERGENIC */ 
        /* Not allowing special sources.. skipping SPECIAL_CDS */ 
        temp = cscore - (GNE_CDS(mat->evi,i,mat->gen,j)) -  (0); 
        if( temp == GenomeWise9_HIDDEN_MATRIX(mat,i - 0,j - 3,START_CODON) ) {  
          *reti = i - 0; 
          *retj = j - 3; 
          *retstate = START_CODON;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GenomeWise9_HIDDEN_MATRIX(mat,i-0,j-3,START_CODON);    
            }  
          return GenomeWise9_HIDDEN_MATRIX(mat,i - 0,j - 3,START_CODON);     
          }  
        temp = cscore - ((GNE_CDS(mat->evi,i,mat->gen,j)+mat->non_start_codon)) -  (0);  
        if( temp == GenomeWise9_HIDDEN_MATRIX(mat,i - 0,j - 3,UTR5) )    {  
          *reti = i - 0; 
          *retj = j - 3; 
          *retstate = UTR5;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GenomeWise9_HIDDEN_MATRIX(mat,i-0,j-3,UTR5);   
            }  
          return GenomeWise9_HIDDEN_MATRIX(mat,i - 0,j - 3,UTR5);    
          }  
        temp = cscore - (GNE_CDS_FRAMESHIFT(mat->evi,i,mat->gen,j,4)) -  (0);    
        if( temp == GenomeWise9_HIDDEN_MATRIX(mat,i - 0,j - 4,CDS) ) {  
          *reti = i - 0; 
          *retj = j - 4; 
          *retstate = CDS;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GenomeWise9_HIDDEN_MATRIX(mat,i-0,j-4,CDS);    
            }  
          return GenomeWise9_HIDDEN_MATRIX(mat,i - 0,j - 4,CDS);     
          }  
        temp = cscore - (GNE_CDS_FRAMESHIFT(mat->evi,i,mat->gen,j,2)) -  (0);    
        if( temp == GenomeWise9_HIDDEN_MATRIX(mat,i - 0,j - 2,CDS) ) {  
          *reti = i - 0; 
          *retj = j - 2; 
          *retstate = CDS;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GenomeWise9_HIDDEN_MATRIX(mat,i-0,j-2,CDS);    
            }  
          return GenomeWise9_HIDDEN_MATRIX(mat,i - 0,j - 2,CDS);     
          }  
        temp = cscore - (GNE_CDS_3SS(mat->evi,i,mat->gen,(j-1),2)) -  (0);   
        if( temp == GenomeWise9_HIDDEN_MATRIX(mat,i - 0,j - 4,CDS_INTRON_2) )    {  
          *reti = i - 0; 
          *retj = j - 4; 
          *retstate = CDS_INTRON_2;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GenomeWise9_HIDDEN_MATRIX(mat,i-0,j-4,CDS_INTRON_2);   
            }  
          return GenomeWise9_HIDDEN_MATRIX(mat,i - 0,j - 4,CDS_INTRON_2);    
          }  
        temp = cscore - (GNE_CDS_3SS(mat->evi,i,mat->gen,(j-2),1)) -  (0);   
        if( temp == GenomeWise9_HIDDEN_MATRIX(mat,i - 0,j - 5,CDS_INTRON_1) )    {  
          *reti = i - 0; 
          *retj = j - 5; 
          *retstate = CDS_INTRON_1;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GenomeWise9_HIDDEN_MATRIX(mat,i-0,j-5,CDS_INTRON_1);   
            }  
          return GenomeWise9_HIDDEN_MATRIX(mat,i - 0,j - 5,CDS_INTRON_1);    
          }  
        temp = cscore - (GNE_CDS_3SS(mat->evi,i,mat->gen,(j-3),0)) -  (0);   
        if( temp == GenomeWise9_HIDDEN_MATRIX(mat,i - 0,j - 6,CDS_INTRON_0) )    {  
          *reti = i - 0; 
          *retj = j - 6; 
          *retstate = CDS_INTRON_0;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GenomeWise9_HIDDEN_MATRIX(mat,i-0,j-6,CDS_INTRON_0);   
            }  
          return GenomeWise9_HIDDEN_MATRIX(mat,i - 0,j - 6,CDS_INTRON_0);    
          }  
        temp = cscore - (GNE_CDS(mat->evi,i,mat->gen,j)) -  (0); 
        if( temp == GenomeWise9_HIDDEN_MATRIX(mat,i - 0,j - 3,CDS) ) {  
          *reti = i - 0; 
          *retj = j - 3; 
          *retstate = CDS;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GenomeWise9_HIDDEN_MATRIX(mat,i-0,j-3,CDS);    
            }  
          return GenomeWise9_HIDDEN_MATRIX(mat,i - 0,j - 3,CDS);     
          }  
        warn("Major problem (!) - in GenomeWise9 read off, position %d,%d state %d no source found!",i,j,state); 
        return (-1); 
      case CDS_INTRON_0 :    
        temp = cscore - (GNE_CDS_5SS(mat->evi,i,mat->gen,(j-7),0)) -  (0);   
        if( temp == GenomeWise9_HIDDEN_MATRIX(mat,i - 0,j - 8,CDS) ) {  
          *reti = i - 0; 
          *retj = j - 8; 
          *retstate = CDS;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GenomeWise9_HIDDEN_MATRIX(mat,i-0,j-8,CDS);    
            }  
          return GenomeWise9_HIDDEN_MATRIX(mat,i - 0,j - 8,CDS);     
          }  
        temp = cscore - (GNE_CDS_INTRON(mat->evi,i,mat->gen,j)) -  (0);  
        if( temp == GenomeWise9_HIDDEN_MATRIX(mat,i - 0,j - 1,CDS_INTRON_0) )    {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = CDS_INTRON_0;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GenomeWise9_HIDDEN_MATRIX(mat,i-0,j-1,CDS_INTRON_0);   
            }  
          return GenomeWise9_HIDDEN_MATRIX(mat,i - 0,j - 1,CDS_INTRON_0);    
          }  
        warn("Major problem (!) - in GenomeWise9 read off, position %d,%d state %d no source found!",i,j,state); 
        return (-1); 
      case CDS_INTRON_1 :    
        temp = cscore - (GNE_CDS_5SS(mat->evi,i,mat->gen,(j-7),1)) -  (0);   
        if( temp == GenomeWise9_HIDDEN_MATRIX(mat,i - 0,j - 9,CDS) ) {  
          *reti = i - 0; 
          *retj = j - 9; 
          *retstate = CDS;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GenomeWise9_HIDDEN_MATRIX(mat,i-0,j-9,CDS);    
            }  
          return GenomeWise9_HIDDEN_MATRIX(mat,i - 0,j - 9,CDS);     
          }  
        temp = cscore - (GNE_CDS_INTRON(mat->evi,i,mat->gen,j)) -  (0);  
        if( temp == GenomeWise9_HIDDEN_MATRIX(mat,i - 0,j - 1,CDS_INTRON_1) )    {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = CDS_INTRON_1;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GenomeWise9_HIDDEN_MATRIX(mat,i-0,j-1,CDS_INTRON_1);   
            }  
          return GenomeWise9_HIDDEN_MATRIX(mat,i - 0,j - 1,CDS_INTRON_1);    
          }  
        warn("Major problem (!) - in GenomeWise9 read off, position %d,%d state %d no source found!",i,j,state); 
        return (-1); 
      case CDS_INTRON_2 :    
        temp = cscore - (GNE_CDS_5SS(mat->evi,i,mat->gen,(j-7),2)) -  (0);   
        if( temp == GenomeWise9_HIDDEN_MATRIX(mat,i - 0,j - 10,CDS) )    {  
          *reti = i - 0; 
          *retj = j - 10;    
          *retstate = CDS;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GenomeWise9_HIDDEN_MATRIX(mat,i-0,j-10,CDS);   
            }  
          return GenomeWise9_HIDDEN_MATRIX(mat,i - 0,j - 10,CDS);    
          }  
        temp = cscore - (GNE_CDS_INTRON(mat->evi,i,mat->gen,j)) -  (0);  
        if( temp == GenomeWise9_HIDDEN_MATRIX(mat,i - 0,j - 1,CDS_INTRON_2) )    {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = CDS_INTRON_2;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GenomeWise9_HIDDEN_MATRIX(mat,i-0,j-1,CDS_INTRON_2);   
            }  
          return GenomeWise9_HIDDEN_MATRIX(mat,i - 0,j - 1,CDS_INTRON_2);    
          }  
        warn("Major problem (!) - in GenomeWise9 read off, position %d,%d state %d no source found!",i,j,state); 
        return (-1); 
      case STOP_CODON :  
        temp = cscore - (GNE_STOP_CODON(mat->evi,i,mat->gen,j)) -  (0);  
        if( temp == GenomeWise9_HIDDEN_MATRIX(mat,i - 0,j - 3,CDS) ) {  
          *reti = i - 0; 
          *retj = j - 3; 
          *retstate = CDS;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GenomeWise9_HIDDEN_MATRIX(mat,i-0,j-3,CDS);    
            }  
          return GenomeWise9_HIDDEN_MATRIX(mat,i - 0,j - 3,CDS);     
          }  
        warn("Major problem (!) - in GenomeWise9 read off, position %d,%d state %d no source found!",i,j,state); 
        return (-1); 
      case UTR3 :    
        /* Not allowing special sources.. skipping INTERGENIC */ 
        temp = cscore - (GNE_UTR_3SS(mat->evi,i,mat->gen,j)) -  (0); 
        if( temp == GenomeWise9_HIDDEN_MATRIX(mat,i - 0,j - 1,UTR3_INTRON) ) {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = UTR3_INTRON;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GenomeWise9_HIDDEN_MATRIX(mat,i-0,j-1,UTR3_INTRON);    
            }  
          return GenomeWise9_HIDDEN_MATRIX(mat,i - 0,j - 1,UTR3_INTRON);     
          }  
        temp = cscore - (GNE_UTR(mat->evi,i,mat->gen,j)) -  (0); 
        if( temp == GenomeWise9_HIDDEN_MATRIX(mat,i - 0,j - 1,STOP_CODON) )  {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = STOP_CODON;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GenomeWise9_HIDDEN_MATRIX(mat,i-0,j-1,STOP_CODON); 
            }  
          return GenomeWise9_HIDDEN_MATRIX(mat,i - 0,j - 1,STOP_CODON);  
          }  
        temp = cscore - ((GNE_UTR(mat->evi,i,mat->gen,j)+mat->non_stop_codon)) -  (0);   
        if( temp == GenomeWise9_HIDDEN_MATRIX(mat,i - 0,j - 1,CDS) ) {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = CDS;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GenomeWise9_HIDDEN_MATRIX(mat,i-0,j-1,CDS);    
            }  
          return GenomeWise9_HIDDEN_MATRIX(mat,i - 0,j - 1,CDS);     
          }  
        temp = cscore - (GNE_UTR(mat->evi,i,mat->gen,j)) -  (0); 
        if( temp == GenomeWise9_HIDDEN_MATRIX(mat,i - 0,j - 1,UTR3) )    {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = UTR3;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GenomeWise9_HIDDEN_MATRIX(mat,i-0,j-1,UTR3);   
            }  
          return GenomeWise9_HIDDEN_MATRIX(mat,i - 0,j - 1,UTR3);    
          }  
        warn("Major problem (!) - in GenomeWise9 read off, position %d,%d state %d no source found!",i,j,state); 
        return (-1); 
      case UTR3_INTRON :     
        /* Not allowing special sources.. skipping SPECIAL_UTR3 */ 
        temp = cscore - (GNE_UTR_5SS(mat->evi,i,mat->gen,j)) -  (0); 
        if( temp == GenomeWise9_HIDDEN_MATRIX(mat,i - 0,j - 1,UTR3) )    {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = UTR3;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GenomeWise9_HIDDEN_MATRIX(mat,i-0,j-1,UTR3);   
            }  
          return GenomeWise9_HIDDEN_MATRIX(mat,i - 0,j - 1,UTR3);    
          }  
        temp = cscore - (GNE_UTR_INTRON(mat->evi,i,mat->gen,j)) -  (0);  
        if( temp == GenomeWise9_HIDDEN_MATRIX(mat,i - 0,j - 1,UTR3_INTRON) ) {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = UTR3_INTRON;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GenomeWise9_HIDDEN_MATRIX(mat,i-0,j-1,UTR3_INTRON);    
            }  
          return GenomeWise9_HIDDEN_MATRIX(mat,i - 0,j - 1,UTR3_INTRON);     
          }  
        warn("Major problem (!) - in GenomeWise9 read off, position %d,%d state %d no source found!",i,j,state); 
        return (-1); 
      default:   
        warn("Major problem (!) - in GenomeWise9 read off, position %d,%d state %d no source found!",i,j,state); 
        return (-1); 
      } /* end of Switch state  */ 
}    


/* Function:  read_special_strip_GenomeWise9(mat,stopi,stopj,stopstate,startj,startstate,out)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [GenomeWise9 *]
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
boolean read_special_strip_GenomeWise9(GenomeWise9 * mat,int stopi,int stopj,int stopstate,int * startj,int * startstate,PackAln * out) 
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
    while( j > GenomeWise9_DC_SHADOW_SPECIAL_SP(mat,i,j,state,4) && state != START)  { /*while more specials to eat up*/ 
      /* Put away current state, if we should */ 
      if(out != NULL)    {  
        pau = PackAlnUnit_alloc();  /* Should deal with memory overflow */ 
        pau->i = i;  
        pau->j = j;  
        pau->state =  state + 10;    
        add_PackAln(out,pau);    
        }  


      max_special_strip_GenomeWise9(mat,i,j,state,isspecial,&i,&j,&state,&isspecial,&cellscore); 
      if( i == GenomeWise9_READ_OFF_ERROR)   {  
        warn("In special strip read GenomeWise9, got a bad read off error. Sorry!"); 
        return FALSE;    
        }  
      } /* end of while more specials to eat up */ 


    /* check to see we have not gone too far! */ 
    if( state != START && j < GenomeWise9_DC_SHADOW_SPECIAL_SP(mat,i,j,state,4)) {  
      warn("In special strip read GenomeWise9, at special [%d] state [%d] overshot!",j,state);   
      return FALSE;  
      }  
    /* Put away last state */ 
    if(out != NULL)  {  
      pau = PackAlnUnit_alloc();/* Should deal with memory overflow */ 
      pau->i = i;    
      pau->j = j;    
      pau->state =  state + 10;  
      add_PackAln(out,pau);  
      }  


    /* Put away where we are in startj and startstate */ 
    *startj = j; 
    *startstate = state; 
    return TRUE; 
}    


/* Function:  max_special_strip_GenomeWise9(mat,i,j,state,isspecial,reti,retj,retstate,retspecial,cellscore)
 *
 * Descrip:    A pretty intense internal function. Deals with read-off only in specials
 *
 *
 * Arg:               mat [UNKN ] Undocumented argument [GenomeWise9 *]
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
int max_special_strip_GenomeWise9(GenomeWise9 * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore) 
{
    int temp;    
    int cscore;  


    *reti = (*retj) = (*retstate) = GenomeWise9_READ_OFF_ERROR;  
    if( isspecial == FALSE ) {  
      warn("In special strip max function for GenomeWise9, got a non special start point. Problem! (bad!)"); 
      return (-1);   
      }  


    if( j < 0 || j > mat->gen->seq->len) {  
      warn("In GenomeWise9 matrix special read off - out of bounds on matrix [j is %d in special]",j);   
      return -1; 
      }  


    cscore = GenomeWise9_DC_SHADOW_SPECIAL(mat,i,j,state);   
    switch(state)    { /*switch on special states*/ 
      case PREGENE_INTERGENIC :  
        /* source START is a special */ 
        temp = cscore - (0) - (0);   
        if( temp == GenomeWise9_DC_SHADOW_SPECIAL(mat,i - 0,j - 1,START) )   {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = START; 
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - GenomeWise9_DC_SHADOW_SPECIAL(mat,i-0,j-1,START);  
            }  
          return GenomeWise9_DC_SHADOW_MATRIX(mat,i - 0,j - 1,START) ;   
          }  
      case POSTGENE_INTERGENIC :     
        /* Source UTR3 is not a special */ 
      case INTERGENIC :  
        /* Source UTR3 is not a special */ 
        /* Source CDS is not a special */ 
        /* source INTERGENIC is a special */ 
        temp = cscore - (0) - (0);   
        if( temp == GenomeWise9_DC_SHADOW_SPECIAL(mat,i - 0,j - 1,INTERGENIC) )  {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = INTERGENIC;    
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - GenomeWise9_DC_SHADOW_SPECIAL(mat,i-0,j-1,INTERGENIC);     
            }  
          return GenomeWise9_DC_SHADOW_MATRIX(mat,i - 0,j - 1,INTERGENIC) ;  
          }  
      case SPECIAL_UTR5 :    
        /* Source UTR5 is not a special */ 
      case SPECIAL_UTR3 :    
        /* Source UTR3 is not a special */ 
      case SPECIAL_CDS :     
        /* source SPECIAL_CDS is a special */ 
        temp = cscore - (mat->rndcodon->codon[CSEQ_GENOMIC_CODON(mat->gen,j)]) - (0);    
        if( temp == GenomeWise9_DC_SHADOW_SPECIAL(mat,i - 0,j - 3,SPECIAL_CDS) ) {  
          *reti = i - 0; 
          *retj = j - 3; 
          *retstate = SPECIAL_CDS;   
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - GenomeWise9_DC_SHADOW_SPECIAL(mat,i-0,j-3,SPECIAL_CDS);    
            }  
          return GenomeWise9_DC_SHADOW_MATRIX(mat,i - 0,j - 3,SPECIAL_CDS) ;     
          }  
        /* Source CDS is not a special */ 
      case START :   
      case END :     
        /* source INTERGENIC is a special */ 
        temp = cscore - (0) - (0);   
        if( temp == GenomeWise9_DC_SHADOW_SPECIAL(mat,i - 0,j - 1,INTERGENIC) )  {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = INTERGENIC;    
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - GenomeWise9_DC_SHADOW_SPECIAL(mat,i-0,j-1,INTERGENIC);     
            }  
          return GenomeWise9_DC_SHADOW_MATRIX(mat,i - 0,j - 1,INTERGENIC) ;  
          }  
      default:   
        warn("Major problem (!) - in GenomeWise9 special strip read off, position %d,%d state %d no source found  dropped into default on source switch!",i,j,state);    
        return (-1); 
      } /* end of switch on special states */ 
}    


/* Function:  max_matrix_to_special_GenomeWise9(mat,i,j,state,cscore,reti,retj,retstate,retspecial,cellscore)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [GenomeWise9 *]
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
int max_matrix_to_special_GenomeWise9(GenomeWise9 * mat,int i,int j,int state,int cscore,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore) 
{
    int temp;    
    *reti = (*retj) = (*retstate) = GenomeWise9_READ_OFF_ERROR;  


    if( j < 0 || j > mat->lenj)  {  
      warn("In GenomeWise9 matrix to special read off - out of bounds on matrix [j is %d in special]",j);    
      return -1; 
      }  


    switch(state)    { /*Switch state */ 
      case UTR5 :    
        temp = cscore - (GNE_UTR(mat->evi,i,mat->gen,j)) -  (0);     
        if( temp == GenomeWise9_DC_SHADOW_SPECIAL(mat,i - 0,j - 1,PREGENE_INTERGENIC) )  {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = PREGENE_INTERGENIC;    
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - GenomeWise9_DC_SHADOW_SPECIAL(mat,i-0,j-1,PREGENE_INTERGENIC);     
            }  
          return GenomeWise9_DC_SHADOW_MATRIX(mat,i - 0,j - 1,PREGENE_INTERGENIC) ;  
          }  
        temp = cscore - ((GNE_UTR(mat->evi,i,mat->gen,j)+GNE_UTR5_START(mat->evi,i,mat->gen,j))) -  (0);     
        if( temp == GenomeWise9_DC_SHADOW_SPECIAL(mat,i - 0,j - 1,INTERGENIC) )  {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = INTERGENIC;    
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - GenomeWise9_DC_SHADOW_SPECIAL(mat,i-0,j-1,INTERGENIC);     
            }  
          return GenomeWise9_DC_SHADOW_MATRIX(mat,i - 0,j - 1,INTERGENIC) ;  
          }  
        temp = cscore - (GNE_UTR5_START(mat->evi,i,mat->gen,j)) -  (0);  
        if( temp == GenomeWise9_DC_SHADOW_SPECIAL(mat,i - 0,j - 1,START) )   {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = START; 
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - GenomeWise9_DC_SHADOW_SPECIAL(mat,i-0,j-1,START);  
            }  
          return GenomeWise9_DC_SHADOW_MATRIX(mat,i - 0,j - 1,START) ;   
          }  
        temp = cscore - (GNE_UTR(mat->evi,i,mat->gen,j)) -  (0);     
        if( temp == GenomeWise9_DC_SHADOW_SPECIAL(mat,i - 0,j - 1,SPECIAL_UTR5) )    {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = SPECIAL_UTR5;  
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - GenomeWise9_DC_SHADOW_SPECIAL(mat,i-0,j-1,SPECIAL_UTR5);   
            }  
          return GenomeWise9_DC_SHADOW_MATRIX(mat,i - 0,j - 1,SPECIAL_UTR5) ;    
          }  
        /* Source UTR5_INTRON is not a special, should not get here! */ 
        /* Source UTR5 is not a special, should not get here! */ 
        warn("Major problem (!) - in GenomeWise9 matrix to special read off, position %d,%d state %d no source found!",i,j,state);   
        return (-1); 
      case UTR5_INTRON :     
        /* Source UTR5 is not a special, should not get here! */ 
        /* Source UTR5_INTRON is not a special, should not get here! */ 
        warn("Major problem (!) - in GenomeWise9 matrix to special read off, position %d,%d state %d no source found!",i,j,state);   
        return (-1); 
      case START_CODON :     
        /* Source UTR5 is not a special, should not get here! */ 
        /* Source UTR5_INTRON is not a special, should not get here! */ 
        warn("Major problem (!) - in GenomeWise9 matrix to special read off, position %d,%d state %d no source found!",i,j,state);   
        return (-1); 
      case CDS :     
        temp = cscore - ((GNE_CDS(mat->evi,i,mat->gen,j)+GNE_UTR5_START(mat->evi,i,mat->gen,j))) -  (0);     
        if( temp == GenomeWise9_DC_SHADOW_SPECIAL(mat,i - 0,j - 3,INTERGENIC) )  {  
          *reti = i - 0; 
          *retj = j - 3; 
          *retstate = INTERGENIC;    
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - GenomeWise9_DC_SHADOW_SPECIAL(mat,i-0,j-3,INTERGENIC);     
            }  
          return GenomeWise9_DC_SHADOW_MATRIX(mat,i - 0,j - 3,INTERGENIC) ;  
          }  
        temp = cscore - (GNE_CDS(mat->evi,i,mat->gen,j)) -  (0);     
        if( temp == GenomeWise9_DC_SHADOW_SPECIAL(mat,i - 0,j - 3,SPECIAL_CDS) ) {  
          *reti = i - 0; 
          *retj = j - 3; 
          *retstate = SPECIAL_CDS;   
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - GenomeWise9_DC_SHADOW_SPECIAL(mat,i-0,j-3,SPECIAL_CDS);    
            }  
          return GenomeWise9_DC_SHADOW_MATRIX(mat,i - 0,j - 3,SPECIAL_CDS) ;     
          }  
        /* Source START_CODON is not a special, should not get here! */ 
        /* Source UTR5 is not a special, should not get here! */ 
        /* Source CDS is not a special, should not get here! */ 
        /* Source CDS is not a special, should not get here! */ 
        /* Source CDS_INTRON_2 is not a special, should not get here! */ 
        /* Source CDS_INTRON_1 is not a special, should not get here! */ 
        /* Source CDS_INTRON_0 is not a special, should not get here! */ 
        /* Source CDS is not a special, should not get here! */ 
        warn("Major problem (!) - in GenomeWise9 matrix to special read off, position %d,%d state %d no source found!",i,j,state);   
        return (-1); 
      case CDS_INTRON_0 :    
        /* Source CDS is not a special, should not get here! */ 
        /* Source CDS_INTRON_0 is not a special, should not get here! */ 
        warn("Major problem (!) - in GenomeWise9 matrix to special read off, position %d,%d state %d no source found!",i,j,state);   
        return (-1); 
      case CDS_INTRON_1 :    
        /* Source CDS is not a special, should not get here! */ 
        /* Source CDS_INTRON_1 is not a special, should not get here! */ 
        warn("Major problem (!) - in GenomeWise9 matrix to special read off, position %d,%d state %d no source found!",i,j,state);   
        return (-1); 
      case CDS_INTRON_2 :    
        /* Source CDS is not a special, should not get here! */ 
        /* Source CDS_INTRON_2 is not a special, should not get here! */ 
        warn("Major problem (!) - in GenomeWise9 matrix to special read off, position %d,%d state %d no source found!",i,j,state);   
        return (-1); 
      case STOP_CODON :  
        /* Source CDS is not a special, should not get here! */ 
        warn("Major problem (!) - in GenomeWise9 matrix to special read off, position %d,%d state %d no source found!",i,j,state);   
        return (-1); 
      case UTR3 :    
        temp = cscore - (GNE_UTR(mat->evi,i,mat->gen,j)) -  (0);     
        if( temp == GenomeWise9_DC_SHADOW_SPECIAL(mat,i - 0,j - 1,INTERGENIC) )  {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = INTERGENIC;    
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - GenomeWise9_DC_SHADOW_SPECIAL(mat,i-0,j-1,INTERGENIC);     
            }  
          return GenomeWise9_DC_SHADOW_MATRIX(mat,i - 0,j - 1,INTERGENIC) ;  
          }  
        /* Source UTR3_INTRON is not a special, should not get here! */ 
        /* Source STOP_CODON is not a special, should not get here! */ 
        /* Source CDS is not a special, should not get here! */ 
        /* Source UTR3 is not a special, should not get here! */ 
        warn("Major problem (!) - in GenomeWise9 matrix to special read off, position %d,%d state %d no source found!",i,j,state);   
        return (-1); 
      case UTR3_INTRON :     
        temp = cscore - (GNE_UTR(mat->evi,i,mat->gen,j)) -  (0);     
        if( temp == GenomeWise9_DC_SHADOW_SPECIAL(mat,i - 0,j - 1,SPECIAL_UTR3) )    {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = SPECIAL_UTR3;  
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - GenomeWise9_DC_SHADOW_SPECIAL(mat,i-0,j-1,SPECIAL_UTR3);   
            }  
          return GenomeWise9_DC_SHADOW_MATRIX(mat,i - 0,j - 1,SPECIAL_UTR3) ;    
          }  
        /* Source UTR3 is not a special, should not get here! */ 
        /* Source UTR3_INTRON is not a special, should not get here! */ 
        warn("Major problem (!) - in GenomeWise9 matrix to special read off, position %d,%d state %d no source found!",i,j,state);   
        return (-1); 
      default:   
        warn("Major problem (!) - in GenomeWise9 read off, position %d,%d state %d no source found!",i,j,state); 
        return (-1); 
      } /* end of Switch state  */ 


}    


/* Function:  calculate_hidden_GenomeWise9(mat,starti,startj,startstate,stopi,stopj,stopstate,dpenv)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [GenomeWise9 *]
 * Arg:            starti [UNKN ] Undocumented argument [int]
 * Arg:            startj [UNKN ] Undocumented argument [int]
 * Arg:        startstate [UNKN ] Undocumented argument [int]
 * Arg:             stopi [UNKN ] Undocumented argument [int]
 * Arg:             stopj [UNKN ] Undocumented argument [int]
 * Arg:         stopstate [UNKN ] Undocumented argument [int]
 * Arg:             dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 */
void calculate_hidden_GenomeWise9(GenomeWise9 * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv) 
{
    register int i;  
    register int j;  
    register int score;  
    register int temp;   
    register int hiddenj;    


    hiddenj = startj;    


    init_hidden_GenomeWise9(mat,starti,startj,stopi,stopj);  


    GenomeWise9_HIDDEN_MATRIX(mat,starti,startj,startstate) = 0; 


    for(j=startj;j<=stopj;j++)   {  
      for(i=starti;i<=stopi;i++) {  
        /* Should *not* do very first cell as this is the one set to zero in one state! */ 
        if( i == starti && j == startj ) 
          continue;  
        if( dpenv != NULL && is_in_DPEnvelope(dpenv,i,j) == FALSE )  { /*Is not in envelope*/ 
          GenomeWise9_HIDDEN_MATRIX(mat,i,j,UTR5) = NEGI;    
          GenomeWise9_HIDDEN_MATRIX(mat,i,j,UTR5_INTRON) = NEGI;     
          GenomeWise9_HIDDEN_MATRIX(mat,i,j,START_CODON) = NEGI;     
          GenomeWise9_HIDDEN_MATRIX(mat,i,j,CDS) = NEGI;     
          GenomeWise9_HIDDEN_MATRIX(mat,i,j,CDS_INTRON_0) = NEGI;    
          GenomeWise9_HIDDEN_MATRIX(mat,i,j,CDS_INTRON_1) = NEGI;    
          GenomeWise9_HIDDEN_MATRIX(mat,i,j,CDS_INTRON_2) = NEGI;    
          GenomeWise9_HIDDEN_MATRIX(mat,i,j,STOP_CODON) = NEGI;  
          GenomeWise9_HIDDEN_MATRIX(mat,i,j,UTR3) = NEGI;    
          GenomeWise9_HIDDEN_MATRIX(mat,i,j,UTR3_INTRON) = NEGI;     
          continue;  
          } /* end of Is not in envelope */ 


        /* For state UTR5 */ 
        /* setting first movement to score */ 
        score = GenomeWise9_HIDDEN_MATRIX(mat,i-0,j-1,UTR5) + GNE_UTR(mat->evi,i,mat->gen,j);    
        /* From state UTR5_INTRON to state UTR5 */ 
        temp = GenomeWise9_HIDDEN_MATRIX(mat,i-0,j-1,UTR5_INTRON) + GNE_UTR_3SS(mat->evi,i,mat->gen,j);  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for UTR5 */ 
        /* Add any movement independant score and put away */ 
         GenomeWise9_HIDDEN_MATRIX(mat,i,j,UTR5) = score;    
        /* Finished calculating state UTR5 */ 


        /* For state UTR5_INTRON */ 
        /* setting first movement to score */ 
        score = GenomeWise9_HIDDEN_MATRIX(mat,i-0,j-1,UTR5_INTRON) + GNE_UTR_INTRON(mat->evi,i,mat->gen,j);  
        /* From state UTR5 to state UTR5_INTRON */ 
        temp = GenomeWise9_HIDDEN_MATRIX(mat,i-0,j-1,UTR5) + GNE_UTR_5SS(mat->evi,i,mat->gen,j);     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for UTR5_INTRON */ 
        /* Add any movement independant score and put away */ 
         GenomeWise9_HIDDEN_MATRIX(mat,i,j,UTR5_INTRON) = score; 
        /* Finished calculating state UTR5_INTRON */ 


        /* For state START_CODON */ 
        /* setting first movement to score */ 
        score = GenomeWise9_HIDDEN_MATRIX(mat,i-0,j-3,UTR5_INTRON) + GNE_START_CODON(mat->evi,i,mat->gen,j);     
        /* From state UTR5 to state START_CODON */ 
        temp = GenomeWise9_HIDDEN_MATRIX(mat,i-0,j-3,UTR5) + GNE_START_CODON(mat->evi,i,mat->gen,j);     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for START_CODON */ 
        /* Add any movement independant score and put away */ 
         GenomeWise9_HIDDEN_MATRIX(mat,i,j,START_CODON) = score; 
        /* Finished calculating state START_CODON */ 


        /* For state CDS */ 
        /* setting first movement to score */ 
        score = GenomeWise9_HIDDEN_MATRIX(mat,i-0,j-3,CDS) + GNE_CDS(mat->evi,i,mat->gen,j);     
        /* From state CDS_INTRON_0 to state CDS */ 
        temp = GenomeWise9_HIDDEN_MATRIX(mat,i-0,j-6,CDS_INTRON_0) + GNE_CDS_3SS(mat->evi,i,mat->gen,(j-3),0);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state CDS_INTRON_1 to state CDS */ 
        temp = GenomeWise9_HIDDEN_MATRIX(mat,i-0,j-5,CDS_INTRON_1) + GNE_CDS_3SS(mat->evi,i,mat->gen,(j-2),1);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state CDS_INTRON_2 to state CDS */ 
        temp = GenomeWise9_HIDDEN_MATRIX(mat,i-0,j-4,CDS_INTRON_2) + GNE_CDS_3SS(mat->evi,i,mat->gen,(j-1),2);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state CDS to state CDS */ 
        temp = GenomeWise9_HIDDEN_MATRIX(mat,i-0,j-2,CDS) + GNE_CDS_FRAMESHIFT(mat->evi,i,mat->gen,j,2);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state CDS to state CDS */ 
        temp = GenomeWise9_HIDDEN_MATRIX(mat,i-0,j-4,CDS) + GNE_CDS_FRAMESHIFT(mat->evi,i,mat->gen,j,4);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state UTR5 to state CDS */ 
        temp = GenomeWise9_HIDDEN_MATRIX(mat,i-0,j-3,UTR5) + (GNE_CDS(mat->evi,i,mat->gen,j)+mat->non_start_codon);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state START_CODON to state CDS */ 
        temp = GenomeWise9_HIDDEN_MATRIX(mat,i-0,j-3,START_CODON) + GNE_CDS(mat->evi,i,mat->gen,j);  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for CDS */ 
        /* Add any movement independant score and put away */ 
         GenomeWise9_HIDDEN_MATRIX(mat,i,j,CDS) = score; 
        /* Finished calculating state CDS */ 


        /* For state CDS_INTRON_0 */ 
        /* setting first movement to score */ 
        score = GenomeWise9_HIDDEN_MATRIX(mat,i-0,j-1,CDS_INTRON_0) + GNE_CDS_INTRON(mat->evi,i,mat->gen,j);     
        /* From state CDS to state CDS_INTRON_0 */ 
        temp = GenomeWise9_HIDDEN_MATRIX(mat,i-0,j-8,CDS) + GNE_CDS_5SS(mat->evi,i,mat->gen,(j-7),0);    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for CDS_INTRON_0 */ 
        /* Add any movement independant score and put away */ 
         GenomeWise9_HIDDEN_MATRIX(mat,i,j,CDS_INTRON_0) = score;    
        /* Finished calculating state CDS_INTRON_0 */ 


        /* For state CDS_INTRON_1 */ 
        /* setting first movement to score */ 
        score = GenomeWise9_HIDDEN_MATRIX(mat,i-0,j-1,CDS_INTRON_1) + GNE_CDS_INTRON(mat->evi,i,mat->gen,j);     
        /* From state CDS to state CDS_INTRON_1 */ 
        temp = GenomeWise9_HIDDEN_MATRIX(mat,i-0,j-9,CDS) + GNE_CDS_5SS(mat->evi,i,mat->gen,(j-7),1);    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for CDS_INTRON_1 */ 
        /* Add any movement independant score and put away */ 
         GenomeWise9_HIDDEN_MATRIX(mat,i,j,CDS_INTRON_1) = score;    
        /* Finished calculating state CDS_INTRON_1 */ 


        /* For state CDS_INTRON_2 */ 
        /* setting first movement to score */ 
        score = GenomeWise9_HIDDEN_MATRIX(mat,i-0,j-1,CDS_INTRON_2) + GNE_CDS_INTRON(mat->evi,i,mat->gen,j);     
        /* From state CDS to state CDS_INTRON_2 */ 
        temp = GenomeWise9_HIDDEN_MATRIX(mat,i-0,j-10,CDS) + GNE_CDS_5SS(mat->evi,i,mat->gen,(j-7),2);   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for CDS_INTRON_2 */ 
        /* Add any movement independant score and put away */ 
         GenomeWise9_HIDDEN_MATRIX(mat,i,j,CDS_INTRON_2) = score;    
        /* Finished calculating state CDS_INTRON_2 */ 


        /* For state STOP_CODON */ 
        /* setting first movement to score */ 
        score = GenomeWise9_HIDDEN_MATRIX(mat,i-0,j-3,CDS) + GNE_STOP_CODON(mat->evi,i,mat->gen,j);  


        /* Ok - finished max calculation for STOP_CODON */ 
        /* Add any movement independant score and put away */ 
         GenomeWise9_HIDDEN_MATRIX(mat,i,j,STOP_CODON) = score;  
        /* Finished calculating state STOP_CODON */ 


        /* For state UTR3 */ 
        /* setting first movement to score */ 
        score = GenomeWise9_HIDDEN_MATRIX(mat,i-0,j-1,UTR3) + GNE_UTR(mat->evi,i,mat->gen,j);    
        /* From state CDS to state UTR3 */ 
        temp = GenomeWise9_HIDDEN_MATRIX(mat,i-0,j-1,CDS) + (GNE_UTR(mat->evi,i,mat->gen,j)+mat->non_stop_codon);    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state STOP_CODON to state UTR3 */ 
        temp = GenomeWise9_HIDDEN_MATRIX(mat,i-0,j-1,STOP_CODON) + GNE_UTR(mat->evi,i,mat->gen,j);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state UTR3_INTRON to state UTR3 */ 
        temp = GenomeWise9_HIDDEN_MATRIX(mat,i-0,j-1,UTR3_INTRON) + GNE_UTR_3SS(mat->evi,i,mat->gen,j);  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for UTR3 */ 
        /* Add any movement independant score and put away */ 
         GenomeWise9_HIDDEN_MATRIX(mat,i,j,UTR3) = score;    
        /* Finished calculating state UTR3 */ 


        /* For state UTR3_INTRON */ 
        /* setting first movement to score */ 
        score = GenomeWise9_HIDDEN_MATRIX(mat,i-0,j-1,UTR3_INTRON) + GNE_UTR_INTRON(mat->evi,i,mat->gen,j);  
        /* From state UTR3 to state UTR3_INTRON */ 
        temp = GenomeWise9_HIDDEN_MATRIX(mat,i-0,j-1,UTR3) + GNE_UTR_5SS(mat->evi,i,mat->gen,j);     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for UTR3_INTRON */ 
        /* Add any movement independant score and put away */ 
         GenomeWise9_HIDDEN_MATRIX(mat,i,j,UTR3_INTRON) = score; 
        /* Finished calculating state UTR3_INTRON */ 
        }  
      }  


    return;  
}    


/* Function:  init_hidden_GenomeWise9(mat,starti,startj,stopi,stopj)
 *
 * Descrip: No Description
 *
 * Arg:           mat [UNKN ] Undocumented argument [GenomeWise9 *]
 * Arg:        starti [UNKN ] Undocumented argument [int]
 * Arg:        startj [UNKN ] Undocumented argument [int]
 * Arg:         stopi [UNKN ] Undocumented argument [int]
 * Arg:         stopj [UNKN ] Undocumented argument [int]
 *
 */
void init_hidden_GenomeWise9(GenomeWise9 * mat,int starti,int startj,int stopi,int stopj) 
{
    register int i;  
    register int j;  
    register int hiddenj;    


    hiddenj = startj;    
    for(j=(startj-10);j<=stopj;j++)  {  
      for(i=(starti-0);i<=stopi;i++) {  
        GenomeWise9_HIDDEN_MATRIX(mat,i,j,UTR5) = NEGI;
 
        GenomeWise9_HIDDEN_MATRIX(mat,i,j,UTR5_INTRON) = NEGI;
  
        GenomeWise9_HIDDEN_MATRIX(mat,i,j,START_CODON) = NEGI;
  
        GenomeWise9_HIDDEN_MATRIX(mat,i,j,CDS) = NEGI;
  
        GenomeWise9_HIDDEN_MATRIX(mat,i,j,CDS_INTRON_0) = NEGI;
 
        GenomeWise9_HIDDEN_MATRIX(mat,i,j,CDS_INTRON_1) = NEGI;
 
        GenomeWise9_HIDDEN_MATRIX(mat,i,j,CDS_INTRON_2) = NEGI;
 
        GenomeWise9_HIDDEN_MATRIX(mat,i,j,STOP_CODON) = NEGI;
   
        GenomeWise9_HIDDEN_MATRIX(mat,i,j,UTR3) = NEGI;
 
        GenomeWise9_HIDDEN_MATRIX(mat,i,j,UTR3_INTRON) = NEGI;
  
        }  
      }  


    return;  
}    


/* Function:  full_dc_GenomeWise9(mat,starti,startj,startstate,stopi,stopj,stopstate,out,donej,totalj,dpenv)
 *
 * Descrip:    The main divide-and-conquor routine. Basically, call /PackAln_calculate_small_GenomeWise9
 *             Not this function, which is pretty hard core. 
 *             Function is given start/end points (in main matrix) for alignment
 *             It does some checks, decides whether start/end in j is small enough for explicit calc
 *               - if yes, calculates it, reads off into PackAln (out), adds the j distance to donej and returns TRUE
 *               - if no,  uses /do_dc_single_pass_GenomeWise9 to get mid-point
 *                          saves midpoint, and calls itself to do right portion then left portion
 *             right then left ensures PackAln is added the 'right' way, ie, back-to-front
 *             returns FALSE on any error, with a warning
 *
 *
 * Arg:               mat [UNKN ] Matrix with small memory implementation [GenomeWise9 *]
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
boolean full_dc_GenomeWise9(GenomeWise9 * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out,int * donej,int totalj,DPEnvelope * dpenv) 
{
    int lstarti; 
    int lstartj; 
    int lstate;  


    if( mat->basematrix->type != BASEMATRIX_TYPE_SHADOW) {  
      warn("*Very* bad error! - non shadow matrix type in full_dc_GenomeWise9"); 
      return FALSE;  
      }  


    if( starti == -1 || startj == -1 || startstate == -1 || stopi == -1 || stopstate == -1)  {  
      warn("In full dc program, passed bad indices, indices passed were %d:%d[%d] to %d:%d[%d]\n",starti,startj,startstate,stopi,stopj,stopstate);   
      return FALSE;  
      }  


    if( stopj - startj < 50) {  
      log_full_error(REPORT,0,"[%d,%d][%d,%d] Explicit read off",starti,startj,stopi,stopj);/* Build hidden explicit matrix */ 
      calculate_hidden_GenomeWise9(mat,starti,startj,startstate,stopi,stopj,stopstate,dpenv);    
      *donej += (stopj - startj);   /* Now read it off into out */ 
      if( read_hidden_GenomeWise9(mat,starti,startj,startstate,stopi,stopj,stopstate,out) == FALSE)  {  
        warn("In full dc, at %d:%d,%d:%d got a bad hidden explicit read off... ",starti,startj,stopi,stopj); 
        return FALSE;    
        }  
      return TRUE;   
      }  


/* In actual divide and conquor */ 
    if( do_dc_single_pass_GenomeWise9(mat,starti,startj,startstate,stopi,stopj,stopstate,dpenv,(int)(*donej*100)/totalj) == FALSE)   {  
      warn("In divide and conquor for GenomeWise9, at bound %d:%d to %d:%d, unable to calculate midpoint. Problem!",starti,startj,stopi,stopj);  
      return FALSE;  
      }  


/* Ok... now we have to call on each side of the matrix */ 
/* We have to retrieve left hand side positions, as they will be vapped by the time we call LHS */ 
    lstarti= GenomeWise9_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,0);   
    lstartj= GenomeWise9_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,1);   
    lstate = GenomeWise9_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,2);   


/* Call on right hand side: this lets us do the correct read off */ 
    if( full_dc_GenomeWise9(mat,GenomeWise9_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,3),GenomeWise9_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,4),GenomeWise9_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,5),stopi,stopj,stopstate,out,donej,totalj,dpenv) == FALSE)   {  
/* Warning already issued, simply chained back up to top */ 
      return FALSE;  
      }  
/* Call on left hand side */ 
    if( full_dc_GenomeWise9(mat,starti,startj,startstate,lstarti,lstartj,lstate,out,donej,totalj,dpenv) == FALSE)    {  
/* Warning already issued, simply chained back up to top */ 
      return FALSE;  
      }  


    return TRUE;     
}    


/* Function:  do_dc_single_pass_GenomeWise9(mat,starti,startj,startstate,stopi,stopj,stopstate,dpenv,perc_done)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [GenomeWise9 *]
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
boolean do_dc_single_pass_GenomeWise9(GenomeWise9 * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv,int perc_done) 
{
    int halfj;   
    halfj = startj + ((stopj - startj)/2);   


    init_dc_GenomeWise9(mat);    


    GenomeWise9_DC_SHADOW_MATRIX(mat,starti,startj,startstate) = 0;  
    run_up_dc_GenomeWise9(mat,starti,stopi,startj,halfj-1,dpenv,perc_done);  
    push_dc_at_merge_GenomeWise9(mat,starti,stopi,halfj,&halfj,dpenv);   
    follow_on_dc_GenomeWise9(mat,starti,stopi,halfj,stopj,dpenv,perc_done);  
    return TRUE; 
}    


/* Function:  push_dc_at_merge_GenomeWise9(mat,starti,stopi,startj,stopj,dpenv)
 *
 * Descrip: No Description
 *
 * Arg:           mat [UNKN ] Undocumented argument [GenomeWise9 *]
 * Arg:        starti [UNKN ] Undocumented argument [int]
 * Arg:         stopi [UNKN ] Undocumented argument [int]
 * Arg:        startj [UNKN ] Undocumented argument [int]
 * Arg:         stopj [UNKN ] Undocumented argument [int *]
 * Arg:         dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 */
void push_dc_at_merge_GenomeWise9(GenomeWise9 * mat,int starti,int stopi,int startj,int * stopj,DPEnvelope * dpenv) 
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
          GenomeWise9_DC_SHADOW_MATRIX(mat,i,j,UTR5) = NEGI;     
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR5,0) = (-100);  
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR5,1) = (-100);  
          GenomeWise9_DC_SHADOW_MATRIX(mat,i,j,UTR5_INTRON) = NEGI;  
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR5_INTRON,0) = (-100);   
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR5_INTRON,1) = (-100);   
          GenomeWise9_DC_SHADOW_MATRIX(mat,i,j,START_CODON) = NEGI;  
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,START_CODON,0) = (-100);   
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,START_CODON,1) = (-100);   
          GenomeWise9_DC_SHADOW_MATRIX(mat,i,j,CDS) = NEGI;  
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS,0) = (-100);   
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS,1) = (-100);   
          GenomeWise9_DC_SHADOW_MATRIX(mat,i,j,CDS_INTRON_0) = NEGI;     
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS_INTRON_0,0) = (-100);  
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS_INTRON_0,1) = (-100);  
          GenomeWise9_DC_SHADOW_MATRIX(mat,i,j,CDS_INTRON_1) = NEGI;     
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS_INTRON_1,0) = (-100);  
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS_INTRON_1,1) = (-100);  
          GenomeWise9_DC_SHADOW_MATRIX(mat,i,j,CDS_INTRON_2) = NEGI;     
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS_INTRON_2,0) = (-100);  
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS_INTRON_2,1) = (-100);  
          GenomeWise9_DC_SHADOW_MATRIX(mat,i,j,STOP_CODON) = NEGI;   
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,STOP_CODON,0) = (-100);    
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,STOP_CODON,1) = (-100);    
          GenomeWise9_DC_SHADOW_MATRIX(mat,i,j,UTR3) = NEGI;     
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR3,0) = (-100);  
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR3,1) = (-100);  
          GenomeWise9_DC_SHADOW_MATRIX(mat,i,j,UTR3_INTRON) = NEGI;  
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR3_INTRON,0) = (-100);   
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR3_INTRON,1) = (-100);   
          continue;  
          } /* end of Is not in envelope */ 


        /* For state UTR5, pushing when j - offj <= mergej */ 
        score = GenomeWise9_DC_SHADOW_MATRIX(mat,i-0,j-1,UTR5) + GNE_UTR(mat->evi,i,mat->gen,j);     
        if( j - 1 <= mergej) {  
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR5,0) = i-0; 
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR5,1) = j-1; 
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR5,2) = UTR5;    
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR5,3) = i;   
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR5,4) = j;   
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR5,5) = UTR5;    
          }  
        else {  
          for(k=0;k<7;k++)   
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR5,k) = GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,UTR5,k);   
          }  


        temp = GenomeWise9_DC_SHADOW_MATRIX(mat,i-0,j-1,UTR5_INTRON) + GNE_UTR_3SS(mat->evi,i,mat->gen,j);   
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR5,0) = i-0;   
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR5,1) = j-1;   
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR5,2) = UTR5_INTRON;   
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR5,3) = i; 
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR5,4) = j; 
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR5,5) = UTR5;  
            }  
          else   {  
            for(k=0;k<7;k++) 
              GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR5,k) = GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,UTR5_INTRON,k);  
            }  
          }  
        /* Add any movement independant score */ 
        GenomeWise9_DC_SHADOW_MATRIX(mat,i,j,UTR5) = score;  
        /* Finished with state UTR5 */ 


        /* For state UTR5_INTRON, pushing when j - offj <= mergej */ 
        score = GenomeWise9_DC_SHADOW_MATRIX(mat,i-0,j-1,UTR5_INTRON) + GNE_UTR_INTRON(mat->evi,i,mat->gen,j);   
        if( j - 1 <= mergej) {  
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR5_INTRON,0) = i-0;  
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR5_INTRON,1) = j-1;  
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR5_INTRON,2) = UTR5_INTRON;  
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR5_INTRON,3) = i;    
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR5_INTRON,4) = j;    
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR5_INTRON,5) = UTR5_INTRON;  
          }  
        else {  
          for(k=0;k<7;k++)   
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR5_INTRON,k) = GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,UTR5_INTRON,k); 
          }  


        temp = GenomeWise9_DC_SHADOW_MATRIX(mat,i-0,j-1,UTR5) + GNE_UTR_5SS(mat->evi,i,mat->gen,j);  
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR5_INTRON,0) = i-0;    
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR5_INTRON,1) = j-1;    
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR5_INTRON,2) = UTR5;   
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR5_INTRON,3) = i;  
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR5_INTRON,4) = j;  
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR5_INTRON,5) = UTR5_INTRON;    
            }  
          else   {  
            for(k=0;k<7;k++) 
              GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR5_INTRON,k) = GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,UTR5,k);  
            }  
          }  
        /* Add any movement independant score */ 
        GenomeWise9_DC_SHADOW_MATRIX(mat,i,j,UTR5_INTRON) = score;   
        /* Finished with state UTR5_INTRON */ 


        /* For state START_CODON, pushing when j - offj <= mergej */ 
        score = GenomeWise9_DC_SHADOW_MATRIX(mat,i-0,j-3,UTR5_INTRON) + GNE_START_CODON(mat->evi,i,mat->gen,j);  
        if( j - 3 <= mergej) {  
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,START_CODON,0) = i-0;  
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,START_CODON,1) = j-3;  
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,START_CODON,2) = UTR5_INTRON;  
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,START_CODON,3) = i;    
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,START_CODON,4) = j;    
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,START_CODON,5) = START_CODON;  
          }  
        else {  
          for(k=0;k<7;k++)   
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,START_CODON,k) = GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 3,UTR5_INTRON,k); 
          }  


        temp = GenomeWise9_DC_SHADOW_MATRIX(mat,i-0,j-3,UTR5) + GNE_START_CODON(mat->evi,i,mat->gen,j);  
        if( temp > score)    {  
          score = temp;  


          if( j - 3 <= mergej)   {  
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,START_CODON,0) = i-0;    
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,START_CODON,1) = j-3;    
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,START_CODON,2) = UTR5;   
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,START_CODON,3) = i;  
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,START_CODON,4) = j;  
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,START_CODON,5) = START_CODON;    
            }  
          else   {  
            for(k=0;k<7;k++) 
              GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,START_CODON,k) = GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 3,UTR5,k);  
            }  
          }  
        /* Add any movement independant score */ 
        GenomeWise9_DC_SHADOW_MATRIX(mat,i,j,START_CODON) = score;   
        /* Finished with state START_CODON */ 


        /* For state CDS, pushing when j - offj <= mergej */ 
        score = GenomeWise9_DC_SHADOW_MATRIX(mat,i-0,j-3,CDS) + GNE_CDS(mat->evi,i,mat->gen,j);  
        if( j - 3 <= mergej) {  
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS,0) = i-0;  
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS,1) = j-3;  
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS,2) = CDS;  
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS,3) = i;    
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS,4) = j;    
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS,5) = CDS;  
          }  
        else {  
          for(k=0;k<7;k++)   
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS,k) = GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 3,CDS,k); 
          }  


        temp = GenomeWise9_DC_SHADOW_MATRIX(mat,i-0,j-6,CDS_INTRON_0) + GNE_CDS_3SS(mat->evi,i,mat->gen,(j-3),0);    
        if( temp > score)    {  
          score = temp;  


          if( j - 6 <= mergej)   {  
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS,0) = i-0;    
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS,1) = j-6;    
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS,2) = CDS_INTRON_0;   
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS,3) = i;  
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS,4) = j;  
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS,5) = CDS;    
            }  
          else   {  
            for(k=0;k<7;k++) 
              GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS,k) = GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 6,CDS_INTRON_0,k);  
            }  
          }  


        temp = GenomeWise9_DC_SHADOW_MATRIX(mat,i-0,j-5,CDS_INTRON_1) + GNE_CDS_3SS(mat->evi,i,mat->gen,(j-2),1);    
        if( temp > score)    {  
          score = temp;  


          if( j - 5 <= mergej)   {  
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS,0) = i-0;    
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS,1) = j-5;    
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS,2) = CDS_INTRON_1;   
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS,3) = i;  
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS,4) = j;  
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS,5) = CDS;    
            }  
          else   {  
            for(k=0;k<7;k++) 
              GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS,k) = GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 5,CDS_INTRON_1,k);  
            }  
          }  


        temp = GenomeWise9_DC_SHADOW_MATRIX(mat,i-0,j-4,CDS_INTRON_2) + GNE_CDS_3SS(mat->evi,i,mat->gen,(j-1),2);    
        if( temp > score)    {  
          score = temp;  


          if( j - 4 <= mergej)   {  
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS,0) = i-0;    
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS,1) = j-4;    
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS,2) = CDS_INTRON_2;   
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS,3) = i;  
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS,4) = j;  
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS,5) = CDS;    
            }  
          else   {  
            for(k=0;k<7;k++) 
              GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS,k) = GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 4,CDS_INTRON_2,k);  
            }  
          }  


        temp = GenomeWise9_DC_SHADOW_MATRIX(mat,i-0,j-2,CDS) + GNE_CDS_FRAMESHIFT(mat->evi,i,mat->gen,j,2);  
        if( temp > score)    {  
          score = temp;  


          if( j - 2 <= mergej)   {  
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS,0) = i-0;    
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS,1) = j-2;    
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS,2) = CDS;    
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS,3) = i;  
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS,4) = j;  
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS,5) = CDS;    
            }  
          else   {  
            for(k=0;k<7;k++) 
              GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS,k) = GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 2,CDS,k);   
            }  
          }  


        temp = GenomeWise9_DC_SHADOW_MATRIX(mat,i-0,j-4,CDS) + GNE_CDS_FRAMESHIFT(mat->evi,i,mat->gen,j,4);  
        if( temp > score)    {  
          score = temp;  


          if( j - 4 <= mergej)   {  
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS,0) = i-0;    
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS,1) = j-4;    
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS,2) = CDS;    
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS,3) = i;  
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS,4) = j;  
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS,5) = CDS;    
            }  
          else   {  
            for(k=0;k<7;k++) 
              GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS,k) = GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 4,CDS,k);   
            }  
          }  


        temp = GenomeWise9_DC_SHADOW_MATRIX(mat,i-0,j-3,UTR5) + (GNE_CDS(mat->evi,i,mat->gen,j)+mat->non_start_codon);   
        if( temp > score)    {  
          score = temp;  


          if( j - 3 <= mergej)   {  
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS,0) = i-0;    
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS,1) = j-3;    
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS,2) = UTR5;   
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS,3) = i;  
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS,4) = j;  
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS,5) = CDS;    
            }  
          else   {  
            for(k=0;k<7;k++) 
              GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS,k) = GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 3,UTR5,k);  
            }  
          }  


        temp = GenomeWise9_DC_SHADOW_MATRIX(mat,i-0,j-3,START_CODON) + GNE_CDS(mat->evi,i,mat->gen,j);   
        if( temp > score)    {  
          score = temp;  


          if( j - 3 <= mergej)   {  
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS,0) = i-0;    
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS,1) = j-3;    
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS,2) = START_CODON;    
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS,3) = i;  
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS,4) = j;  
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS,5) = CDS;    
            }  
          else   {  
            for(k=0;k<7;k++) 
              GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS,k) = GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 3,START_CODON,k);   
            }  
          }  
        /* Add any movement independant score */ 
        GenomeWise9_DC_SHADOW_MATRIX(mat,i,j,CDS) = score;   
        /* Finished with state CDS */ 


        /* For state CDS_INTRON_0, pushing when j - offj <= mergej */ 
        score = GenomeWise9_DC_SHADOW_MATRIX(mat,i-0,j-1,CDS_INTRON_0) + GNE_CDS_INTRON(mat->evi,i,mat->gen,j);  
        if( j - 1 <= mergej) {  
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS_INTRON_0,0) = i-0; 
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS_INTRON_0,1) = j-1; 
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS_INTRON_0,2) = CDS_INTRON_0;    
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS_INTRON_0,3) = i;   
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS_INTRON_0,4) = j;   
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS_INTRON_0,5) = CDS_INTRON_0;    
          }  
        else {  
          for(k=0;k<7;k++)   
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS_INTRON_0,k) = GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,CDS_INTRON_0,k);   
          }  


        temp = GenomeWise9_DC_SHADOW_MATRIX(mat,i-0,j-8,CDS) + GNE_CDS_5SS(mat->evi,i,mat->gen,(j-7),0);     
        if( temp > score)    {  
          score = temp;  


          if( j - 8 <= mergej)   {  
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS_INTRON_0,0) = i-0;   
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS_INTRON_0,1) = j-8;   
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS_INTRON_0,2) = CDS;   
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS_INTRON_0,3) = i; 
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS_INTRON_0,4) = j; 
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS_INTRON_0,5) = CDS_INTRON_0;  
            }  
          else   {  
            for(k=0;k<7;k++) 
              GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS_INTRON_0,k) = GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 8,CDS,k);  
            }  
          }  
        /* Add any movement independant score */ 
        GenomeWise9_DC_SHADOW_MATRIX(mat,i,j,CDS_INTRON_0) = score;  
        /* Finished with state CDS_INTRON_0 */ 


        /* For state CDS_INTRON_1, pushing when j - offj <= mergej */ 
        score = GenomeWise9_DC_SHADOW_MATRIX(mat,i-0,j-1,CDS_INTRON_1) + GNE_CDS_INTRON(mat->evi,i,mat->gen,j);  
        if( j - 1 <= mergej) {  
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS_INTRON_1,0) = i-0; 
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS_INTRON_1,1) = j-1; 
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS_INTRON_1,2) = CDS_INTRON_1;    
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS_INTRON_1,3) = i;   
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS_INTRON_1,4) = j;   
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS_INTRON_1,5) = CDS_INTRON_1;    
          }  
        else {  
          for(k=0;k<7;k++)   
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS_INTRON_1,k) = GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,CDS_INTRON_1,k);   
          }  


        temp = GenomeWise9_DC_SHADOW_MATRIX(mat,i-0,j-9,CDS) + GNE_CDS_5SS(mat->evi,i,mat->gen,(j-7),1);     
        if( temp > score)    {  
          score = temp;  


          if( j - 9 <= mergej)   {  
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS_INTRON_1,0) = i-0;   
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS_INTRON_1,1) = j-9;   
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS_INTRON_1,2) = CDS;   
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS_INTRON_1,3) = i; 
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS_INTRON_1,4) = j; 
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS_INTRON_1,5) = CDS_INTRON_1;  
            }  
          else   {  
            for(k=0;k<7;k++) 
              GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS_INTRON_1,k) = GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 9,CDS,k);  
            }  
          }  
        /* Add any movement independant score */ 
        GenomeWise9_DC_SHADOW_MATRIX(mat,i,j,CDS_INTRON_1) = score;  
        /* Finished with state CDS_INTRON_1 */ 


        /* For state CDS_INTRON_2, pushing when j - offj <= mergej */ 
        score = GenomeWise9_DC_SHADOW_MATRIX(mat,i-0,j-1,CDS_INTRON_2) + GNE_CDS_INTRON(mat->evi,i,mat->gen,j);  
        if( j - 1 <= mergej) {  
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS_INTRON_2,0) = i-0; 
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS_INTRON_2,1) = j-1; 
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS_INTRON_2,2) = CDS_INTRON_2;    
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS_INTRON_2,3) = i;   
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS_INTRON_2,4) = j;   
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS_INTRON_2,5) = CDS_INTRON_2;    
          }  
        else {  
          for(k=0;k<7;k++)   
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS_INTRON_2,k) = GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,CDS_INTRON_2,k);   
          }  


        temp = GenomeWise9_DC_SHADOW_MATRIX(mat,i-0,j-10,CDS) + GNE_CDS_5SS(mat->evi,i,mat->gen,(j-7),2);    
        if( temp > score)    {  
          score = temp;  


          if( j - 10 <= mergej)  {  
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS_INTRON_2,0) = i-0;   
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS_INTRON_2,1) = j-10;  
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS_INTRON_2,2) = CDS;   
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS_INTRON_2,3) = i; 
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS_INTRON_2,4) = j; 
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS_INTRON_2,5) = CDS_INTRON_2;  
            }  
          else   {  
            for(k=0;k<7;k++) 
              GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS_INTRON_2,k) = GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 10,CDS,k); 
            }  
          }  
        /* Add any movement independant score */ 
        GenomeWise9_DC_SHADOW_MATRIX(mat,i,j,CDS_INTRON_2) = score;  
        /* Finished with state CDS_INTRON_2 */ 


        /* For state STOP_CODON, pushing when j - offj <= mergej */ 
        score = GenomeWise9_DC_SHADOW_MATRIX(mat,i-0,j-3,CDS) + GNE_STOP_CODON(mat->evi,i,mat->gen,j);   
        if( j - 3 <= mergej) {  
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,STOP_CODON,0) = i-0;   
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,STOP_CODON,1) = j-3;   
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,STOP_CODON,2) = CDS;   
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,STOP_CODON,3) = i; 
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,STOP_CODON,4) = j; 
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,STOP_CODON,5) = STOP_CODON;    
          }  
        else {  
          for(k=0;k<7;k++)   
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,STOP_CODON,k) = GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 3,CDS,k);  
          }  
        /* Add any movement independant score */ 
        GenomeWise9_DC_SHADOW_MATRIX(mat,i,j,STOP_CODON) = score;    
        /* Finished with state STOP_CODON */ 


        /* For state UTR3, pushing when j - offj <= mergej */ 
        score = GenomeWise9_DC_SHADOW_MATRIX(mat,i-0,j-1,UTR3) + GNE_UTR(mat->evi,i,mat->gen,j);     
        if( j - 1 <= mergej) {  
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR3,0) = i-0; 
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR3,1) = j-1; 
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR3,2) = UTR3;    
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR3,3) = i;   
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR3,4) = j;   
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR3,5) = UTR3;    
          }  
        else {  
          for(k=0;k<7;k++)   
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR3,k) = GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,UTR3,k);   
          }  


        temp = GenomeWise9_DC_SHADOW_MATRIX(mat,i-0,j-1,CDS) + (GNE_UTR(mat->evi,i,mat->gen,j)+mat->non_stop_codon);     
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR3,0) = i-0;   
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR3,1) = j-1;   
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR3,2) = CDS;   
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR3,3) = i; 
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR3,4) = j; 
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR3,5) = UTR3;  
            }  
          else   {  
            for(k=0;k<7;k++) 
              GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR3,k) = GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,CDS,k);  
            }  
          }  


        temp = GenomeWise9_DC_SHADOW_MATRIX(mat,i-0,j-1,STOP_CODON) + GNE_UTR(mat->evi,i,mat->gen,j);    
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR3,0) = i-0;   
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR3,1) = j-1;   
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR3,2) = STOP_CODON;    
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR3,3) = i; 
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR3,4) = j; 
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR3,5) = UTR3;  
            }  
          else   {  
            for(k=0;k<7;k++) 
              GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR3,k) = GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,STOP_CODON,k);   
            }  
          }  


        temp = GenomeWise9_DC_SHADOW_MATRIX(mat,i-0,j-1,UTR3_INTRON) + GNE_UTR_3SS(mat->evi,i,mat->gen,j);   
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR3,0) = i-0;   
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR3,1) = j-1;   
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR3,2) = UTR3_INTRON;   
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR3,3) = i; 
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR3,4) = j; 
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR3,5) = UTR3;  
            }  
          else   {  
            for(k=0;k<7;k++) 
              GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR3,k) = GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,UTR3_INTRON,k);  
            }  
          }  
        /* Add any movement independant score */ 
        GenomeWise9_DC_SHADOW_MATRIX(mat,i,j,UTR3) = score;  
        /* Finished with state UTR3 */ 


        /* For state UTR3_INTRON, pushing when j - offj <= mergej */ 
        score = GenomeWise9_DC_SHADOW_MATRIX(mat,i-0,j-1,UTR3_INTRON) + GNE_UTR_INTRON(mat->evi,i,mat->gen,j);   
        if( j - 1 <= mergej) {  
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR3_INTRON,0) = i-0;  
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR3_INTRON,1) = j-1;  
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR3_INTRON,2) = UTR3_INTRON;  
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR3_INTRON,3) = i;    
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR3_INTRON,4) = j;    
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR3_INTRON,5) = UTR3_INTRON;  
          }  
        else {  
          for(k=0;k<7;k++)   
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR3_INTRON,k) = GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,UTR3_INTRON,k); 
          }  


        temp = GenomeWise9_DC_SHADOW_MATRIX(mat,i-0,j-1,UTR3) + GNE_UTR_5SS(mat->evi,i,mat->gen,j);  
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR3_INTRON,0) = i-0;    
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR3_INTRON,1) = j-1;    
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR3_INTRON,2) = UTR3;   
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR3_INTRON,3) = i;  
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR3_INTRON,4) = j;  
            GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR3_INTRON,5) = UTR3_INTRON;    
            }  
          else   {  
            for(k=0;k<7;k++) 
              GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR3_INTRON,k) = GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,UTR3,k);  
            }  
          }  
        /* Add any movement independant score */ 
        GenomeWise9_DC_SHADOW_MATRIX(mat,i,j,UTR3_INTRON) = score;   
        /* Finished with state UTR3_INTRON */ 
        }  
      }  
    /* Put back j into * stop j so that calling function gets it correct */ 
    if( stopj == NULL)   
      warn("Bad news... NULL stopj pointer in push dc function. This means that calling function does not know how many cells I have done!");    
    else 
      *stopj = j;    


    return;  
}    


/* Function:  follow_on_dc_GenomeWise9(mat,starti,stopi,startj,stopj,dpenv,perc_done)
 *
 * Descrip: No Description
 *
 * Arg:              mat [UNKN ] Undocumented argument [GenomeWise9 *]
 * Arg:           starti [UNKN ] Undocumented argument [int]
 * Arg:            stopi [UNKN ] Undocumented argument [int]
 * Arg:           startj [UNKN ] Undocumented argument [int]
 * Arg:            stopj [UNKN ] Undocumented argument [int]
 * Arg:            dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:        perc_done [UNKN ] Undocumented argument [int]
 *
 */
void follow_on_dc_GenomeWise9(GenomeWise9 * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done) 
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
          GenomeWise9_DC_SHADOW_MATRIX(mat,i,j,UTR5) = NEGI;     
          GenomeWise9_DC_SHADOW_MATRIX(mat,i,j,UTR5_INTRON) = NEGI;  
          GenomeWise9_DC_SHADOW_MATRIX(mat,i,j,START_CODON) = NEGI;  
          GenomeWise9_DC_SHADOW_MATRIX(mat,i,j,CDS) = NEGI;  
          GenomeWise9_DC_SHADOW_MATRIX(mat,i,j,CDS_INTRON_0) = NEGI;     
          GenomeWise9_DC_SHADOW_MATRIX(mat,i,j,CDS_INTRON_1) = NEGI;     
          GenomeWise9_DC_SHADOW_MATRIX(mat,i,j,CDS_INTRON_2) = NEGI;     
          GenomeWise9_DC_SHADOW_MATRIX(mat,i,j,STOP_CODON) = NEGI;   
          GenomeWise9_DC_SHADOW_MATRIX(mat,i,j,UTR3) = NEGI;     
          GenomeWise9_DC_SHADOW_MATRIX(mat,i,j,UTR3_INTRON) = NEGI;  
          continue;  
          } /* end of Is not in envelope */ 
        if( num % 1000 == 0 )    
          log_full_error(REPORT,0,"[%d%%%% done]After  mid-j %5d Cells done %d%%%%",perc_done,startj,(num*100)/total);   


        /* For state UTR5 */ 
        /* setting first movement to score */ 
        score = GenomeWise9_DC_SHADOW_MATRIX(mat,i-0,j-1,UTR5) + GNE_UTR(mat->evi,i,mat->gen,j);     
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,UTR5,k);  
        /* From state UTR5_INTRON to state UTR5 */ 
        temp = GenomeWise9_DC_SHADOW_MATRIX(mat,i-0,j-1,UTR5_INTRON) + GNE_UTR_3SS(mat->evi,i,mat->gen,j);   
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,UTR5_INTRON,k); 
          }  


        /* Ok - finished max calculation for UTR5 */ 
        /* Add any movement independant score and put away */ 
         GenomeWise9_DC_SHADOW_MATRIX(mat,i,j,UTR5) = score; 
        for(k=0;k<7;k++) 
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR5,k) = localshadow[k];  
        /* Now figure out if any specials need this score */ 
        /* Finished calculating state UTR5 */ 


        /* For state UTR5_INTRON */ 
        /* setting first movement to score */ 
        score = GenomeWise9_DC_SHADOW_MATRIX(mat,i-0,j-1,UTR5_INTRON) + GNE_UTR_INTRON(mat->evi,i,mat->gen,j);   
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,UTR5_INTRON,k);   
        /* From state UTR5 to state UTR5_INTRON */ 
        temp = GenomeWise9_DC_SHADOW_MATRIX(mat,i-0,j-1,UTR5) + GNE_UTR_5SS(mat->evi,i,mat->gen,j);  
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,UTR5,k);    
          }  


        /* Ok - finished max calculation for UTR5_INTRON */ 
        /* Add any movement independant score and put away */ 
         GenomeWise9_DC_SHADOW_MATRIX(mat,i,j,UTR5_INTRON) = score;  
        for(k=0;k<7;k++) 
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR5_INTRON,k) = localshadow[k];   
        /* Now figure out if any specials need this score */ 
        /* Finished calculating state UTR5_INTRON */ 


        /* For state START_CODON */ 
        /* setting first movement to score */ 
        score = GenomeWise9_DC_SHADOW_MATRIX(mat,i-0,j-3,UTR5_INTRON) + GNE_START_CODON(mat->evi,i,mat->gen,j);  
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 3,UTR5_INTRON,k);   
        /* From state UTR5 to state START_CODON */ 
        temp = GenomeWise9_DC_SHADOW_MATRIX(mat,i-0,j-3,UTR5) + GNE_START_CODON(mat->evi,i,mat->gen,j);  
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 3,UTR5,k);    
          }  


        /* Ok - finished max calculation for START_CODON */ 
        /* Add any movement independant score and put away */ 
         GenomeWise9_DC_SHADOW_MATRIX(mat,i,j,START_CODON) = score;  
        for(k=0;k<7;k++) 
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,START_CODON,k) = localshadow[k];   
        /* Now figure out if any specials need this score */ 
        /* Finished calculating state START_CODON */ 


        /* For state CDS */ 
        /* setting first movement to score */ 
        score = GenomeWise9_DC_SHADOW_MATRIX(mat,i-0,j-3,CDS) + GNE_CDS(mat->evi,i,mat->gen,j);  
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 3,CDS,k);   
        /* From state CDS_INTRON_0 to state CDS */ 
        temp = GenomeWise9_DC_SHADOW_MATRIX(mat,i-0,j-6,CDS_INTRON_0) + GNE_CDS_3SS(mat->evi,i,mat->gen,(j-3),0);    
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 6,CDS_INTRON_0,k);    
          }  
        /* From state CDS_INTRON_1 to state CDS */ 
        temp = GenomeWise9_DC_SHADOW_MATRIX(mat,i-0,j-5,CDS_INTRON_1) + GNE_CDS_3SS(mat->evi,i,mat->gen,(j-2),1);    
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 5,CDS_INTRON_1,k);    
          }  
        /* From state CDS_INTRON_2 to state CDS */ 
        temp = GenomeWise9_DC_SHADOW_MATRIX(mat,i-0,j-4,CDS_INTRON_2) + GNE_CDS_3SS(mat->evi,i,mat->gen,(j-1),2);    
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 4,CDS_INTRON_2,k);    
          }  
        /* From state CDS to state CDS */ 
        temp = GenomeWise9_DC_SHADOW_MATRIX(mat,i-0,j-2,CDS) + GNE_CDS_FRAMESHIFT(mat->evi,i,mat->gen,j,2);  
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 2,CDS,k); 
          }  
        /* From state CDS to state CDS */ 
        temp = GenomeWise9_DC_SHADOW_MATRIX(mat,i-0,j-4,CDS) + GNE_CDS_FRAMESHIFT(mat->evi,i,mat->gen,j,4);  
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 4,CDS,k); 
          }  
        /* From state UTR5 to state CDS */ 
        temp = GenomeWise9_DC_SHADOW_MATRIX(mat,i-0,j-3,UTR5) + (GNE_CDS(mat->evi,i,mat->gen,j)+mat->non_start_codon);   
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 3,UTR5,k);    
          }  
        /* From state START_CODON to state CDS */ 
        temp = GenomeWise9_DC_SHADOW_MATRIX(mat,i-0,j-3,START_CODON) + GNE_CDS(mat->evi,i,mat->gen,j);   
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 3,START_CODON,k); 
          }  


        /* Ok - finished max calculation for CDS */ 
        /* Add any movement independant score and put away */ 
         GenomeWise9_DC_SHADOW_MATRIX(mat,i,j,CDS) = score;  
        for(k=0;k<7;k++) 
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS,k) = localshadow[k];   
        /* Now figure out if any specials need this score */ 
        /* Finished calculating state CDS */ 


        /* For state CDS_INTRON_0 */ 
        /* setting first movement to score */ 
        score = GenomeWise9_DC_SHADOW_MATRIX(mat,i-0,j-1,CDS_INTRON_0) + GNE_CDS_INTRON(mat->evi,i,mat->gen,j);  
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,CDS_INTRON_0,k);  
        /* From state CDS to state CDS_INTRON_0 */ 
        temp = GenomeWise9_DC_SHADOW_MATRIX(mat,i-0,j-8,CDS) + GNE_CDS_5SS(mat->evi,i,mat->gen,(j-7),0);     
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 8,CDS,k); 
          }  


        /* Ok - finished max calculation for CDS_INTRON_0 */ 
        /* Add any movement independant score and put away */ 
         GenomeWise9_DC_SHADOW_MATRIX(mat,i,j,CDS_INTRON_0) = score; 
        for(k=0;k<7;k++) 
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS_INTRON_0,k) = localshadow[k];  
        /* Now figure out if any specials need this score */ 
        /* Finished calculating state CDS_INTRON_0 */ 


        /* For state CDS_INTRON_1 */ 
        /* setting first movement to score */ 
        score = GenomeWise9_DC_SHADOW_MATRIX(mat,i-0,j-1,CDS_INTRON_1) + GNE_CDS_INTRON(mat->evi,i,mat->gen,j);  
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,CDS_INTRON_1,k);  
        /* From state CDS to state CDS_INTRON_1 */ 
        temp = GenomeWise9_DC_SHADOW_MATRIX(mat,i-0,j-9,CDS) + GNE_CDS_5SS(mat->evi,i,mat->gen,(j-7),1);     
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 9,CDS,k); 
          }  


        /* Ok - finished max calculation for CDS_INTRON_1 */ 
        /* Add any movement independant score and put away */ 
         GenomeWise9_DC_SHADOW_MATRIX(mat,i,j,CDS_INTRON_1) = score; 
        for(k=0;k<7;k++) 
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS_INTRON_1,k) = localshadow[k];  
        /* Now figure out if any specials need this score */ 
        /* Finished calculating state CDS_INTRON_1 */ 


        /* For state CDS_INTRON_2 */ 
        /* setting first movement to score */ 
        score = GenomeWise9_DC_SHADOW_MATRIX(mat,i-0,j-1,CDS_INTRON_2) + GNE_CDS_INTRON(mat->evi,i,mat->gen,j);  
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,CDS_INTRON_2,k);  
        /* From state CDS to state CDS_INTRON_2 */ 
        temp = GenomeWise9_DC_SHADOW_MATRIX(mat,i-0,j-10,CDS) + GNE_CDS_5SS(mat->evi,i,mat->gen,(j-7),2);    
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 10,CDS,k);    
          }  


        /* Ok - finished max calculation for CDS_INTRON_2 */ 
        /* Add any movement independant score and put away */ 
         GenomeWise9_DC_SHADOW_MATRIX(mat,i,j,CDS_INTRON_2) = score; 
        for(k=0;k<7;k++) 
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS_INTRON_2,k) = localshadow[k];  
        /* Now figure out if any specials need this score */ 
        /* Finished calculating state CDS_INTRON_2 */ 


        /* For state STOP_CODON */ 
        /* setting first movement to score */ 
        score = GenomeWise9_DC_SHADOW_MATRIX(mat,i-0,j-3,CDS) + GNE_STOP_CODON(mat->evi,i,mat->gen,j);   
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 3,CDS,k);   


        /* Ok - finished max calculation for STOP_CODON */ 
        /* Add any movement independant score and put away */ 
         GenomeWise9_DC_SHADOW_MATRIX(mat,i,j,STOP_CODON) = score;   
        for(k=0;k<7;k++) 
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,STOP_CODON,k) = localshadow[k];    
        /* Now figure out if any specials need this score */ 
        /* Finished calculating state STOP_CODON */ 


        /* For state UTR3 */ 
        /* setting first movement to score */ 
        score = GenomeWise9_DC_SHADOW_MATRIX(mat,i-0,j-1,UTR3) + GNE_UTR(mat->evi,i,mat->gen,j);     
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,UTR3,k);  
        /* From state CDS to state UTR3 */ 
        temp = GenomeWise9_DC_SHADOW_MATRIX(mat,i-0,j-1,CDS) + (GNE_UTR(mat->evi,i,mat->gen,j)+mat->non_stop_codon);     
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,CDS,k); 
          }  
        /* From state STOP_CODON to state UTR3 */ 
        temp = GenomeWise9_DC_SHADOW_MATRIX(mat,i-0,j-1,STOP_CODON) + GNE_UTR(mat->evi,i,mat->gen,j);    
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,STOP_CODON,k);  
          }  
        /* From state UTR3_INTRON to state UTR3 */ 
        temp = GenomeWise9_DC_SHADOW_MATRIX(mat,i-0,j-1,UTR3_INTRON) + GNE_UTR_3SS(mat->evi,i,mat->gen,j);   
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,UTR3_INTRON,k); 
          }  


        /* Ok - finished max calculation for UTR3 */ 
        /* Add any movement independant score and put away */ 
         GenomeWise9_DC_SHADOW_MATRIX(mat,i,j,UTR3) = score; 
        for(k=0;k<7;k++) 
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR3,k) = localshadow[k];  
        /* Now figure out if any specials need this score */ 
        /* Finished calculating state UTR3 */ 


        /* For state UTR3_INTRON */ 
        /* setting first movement to score */ 
        score = GenomeWise9_DC_SHADOW_MATRIX(mat,i-0,j-1,UTR3_INTRON) + GNE_UTR_INTRON(mat->evi,i,mat->gen,j);   
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,UTR3_INTRON,k);   
        /* From state UTR3 to state UTR3_INTRON */ 
        temp = GenomeWise9_DC_SHADOW_MATRIX(mat,i-0,j-1,UTR3) + GNE_UTR_5SS(mat->evi,i,mat->gen,j);  
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,UTR3,k);    
          }  


        /* Ok - finished max calculation for UTR3_INTRON */ 
        /* Add any movement independant score and put away */ 
         GenomeWise9_DC_SHADOW_MATRIX(mat,i,j,UTR3_INTRON) = score;  
        for(k=0;k<7;k++) 
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR3_INTRON,k) = localshadow[k];   
        /* Now figure out if any specials need this score */ 
        /* Finished calculating state UTR3_INTRON */ 
        } /* end of this is strip */ 
      } /* end of for each valid j column */ 


/* Function:  run_up_dc_GenomeWise9(mat,starti,stopi,startj,stopj,dpenv,perc_done)
 *
 * Descrip: No Description
 *
 * Arg:              mat [UNKN ] Undocumented argument [GenomeWise9 *]
 * Arg:           starti [UNKN ] Undocumented argument [int]
 * Arg:            stopi [UNKN ] Undocumented argument [int]
 * Arg:           startj [UNKN ] Undocumented argument [int]
 * Arg:            stopj [UNKN ] Undocumented argument [int]
 * Arg:            dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:        perc_done [UNKN ] Undocumented argument [int]
 *
 */
}    
void run_up_dc_GenomeWise9(GenomeWise9 * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done) 
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
          GenomeWise9_DC_SHADOW_MATRIX(mat,i,j,UTR5) = NEGI;     
          GenomeWise9_DC_SHADOW_MATRIX(mat,i,j,UTR5_INTRON) = NEGI;  
          GenomeWise9_DC_SHADOW_MATRIX(mat,i,j,START_CODON) = NEGI;  
          GenomeWise9_DC_SHADOW_MATRIX(mat,i,j,CDS) = NEGI;  
          GenomeWise9_DC_SHADOW_MATRIX(mat,i,j,CDS_INTRON_0) = NEGI;     
          GenomeWise9_DC_SHADOW_MATRIX(mat,i,j,CDS_INTRON_1) = NEGI;     
          GenomeWise9_DC_SHADOW_MATRIX(mat,i,j,CDS_INTRON_2) = NEGI;     
          GenomeWise9_DC_SHADOW_MATRIX(mat,i,j,STOP_CODON) = NEGI;   
          GenomeWise9_DC_SHADOW_MATRIX(mat,i,j,UTR3) = NEGI;     
          GenomeWise9_DC_SHADOW_MATRIX(mat,i,j,UTR3_INTRON) = NEGI;  
          continue;  
          } /* end of Is not in envelope */ 
        if( num % 1000 == 0 )    
          log_full_error(REPORT,0,"[%d%%%% done]Before mid-j %5d Cells done %d%%%%",perc_done,stopj,(num*100)/total);    


        /* For state UTR5 */ 
        /* setting first movement to score */ 
        score = GenomeWise9_DC_SHADOW_MATRIX(mat,i-0,j-1,UTR5) + GNE_UTR(mat->evi,i,mat->gen,j);     
        /* From state UTR5_INTRON to state UTR5 */ 
        temp = GenomeWise9_DC_SHADOW_MATRIX(mat,i-0,j-1,UTR5_INTRON) + GNE_UTR_3SS(mat->evi,i,mat->gen,j);   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for UTR5 */ 
        /* Add any movement independant score and put away */ 
         GenomeWise9_DC_SHADOW_MATRIX(mat,i,j,UTR5) = score; 
        /* Finished calculating state UTR5 */ 


        /* For state UTR5_INTRON */ 
        /* setting first movement to score */ 
        score = GenomeWise9_DC_SHADOW_MATRIX(mat,i-0,j-1,UTR5_INTRON) + GNE_UTR_INTRON(mat->evi,i,mat->gen,j);   
        /* From state UTR5 to state UTR5_INTRON */ 
        temp = GenomeWise9_DC_SHADOW_MATRIX(mat,i-0,j-1,UTR5) + GNE_UTR_5SS(mat->evi,i,mat->gen,j);  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for UTR5_INTRON */ 
        /* Add any movement independant score and put away */ 
         GenomeWise9_DC_SHADOW_MATRIX(mat,i,j,UTR5_INTRON) = score;  
        /* Finished calculating state UTR5_INTRON */ 


        /* For state START_CODON */ 
        /* setting first movement to score */ 
        score = GenomeWise9_DC_SHADOW_MATRIX(mat,i-0,j-3,UTR5_INTRON) + GNE_START_CODON(mat->evi,i,mat->gen,j);  
        /* From state UTR5 to state START_CODON */ 
        temp = GenomeWise9_DC_SHADOW_MATRIX(mat,i-0,j-3,UTR5) + GNE_START_CODON(mat->evi,i,mat->gen,j);  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for START_CODON */ 
        /* Add any movement independant score and put away */ 
         GenomeWise9_DC_SHADOW_MATRIX(mat,i,j,START_CODON) = score;  
        /* Finished calculating state START_CODON */ 


        /* For state CDS */ 
        /* setting first movement to score */ 
        score = GenomeWise9_DC_SHADOW_MATRIX(mat,i-0,j-3,CDS) + GNE_CDS(mat->evi,i,mat->gen,j);  
        /* From state CDS_INTRON_0 to state CDS */ 
        temp = GenomeWise9_DC_SHADOW_MATRIX(mat,i-0,j-6,CDS_INTRON_0) + GNE_CDS_3SS(mat->evi,i,mat->gen,(j-3),0);    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state CDS_INTRON_1 to state CDS */ 
        temp = GenomeWise9_DC_SHADOW_MATRIX(mat,i-0,j-5,CDS_INTRON_1) + GNE_CDS_3SS(mat->evi,i,mat->gen,(j-2),1);    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state CDS_INTRON_2 to state CDS */ 
        temp = GenomeWise9_DC_SHADOW_MATRIX(mat,i-0,j-4,CDS_INTRON_2) + GNE_CDS_3SS(mat->evi,i,mat->gen,(j-1),2);    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state CDS to state CDS */ 
        temp = GenomeWise9_DC_SHADOW_MATRIX(mat,i-0,j-2,CDS) + GNE_CDS_FRAMESHIFT(mat->evi,i,mat->gen,j,2);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state CDS to state CDS */ 
        temp = GenomeWise9_DC_SHADOW_MATRIX(mat,i-0,j-4,CDS) + GNE_CDS_FRAMESHIFT(mat->evi,i,mat->gen,j,4);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state UTR5 to state CDS */ 
        temp = GenomeWise9_DC_SHADOW_MATRIX(mat,i-0,j-3,UTR5) + (GNE_CDS(mat->evi,i,mat->gen,j)+mat->non_start_codon);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state START_CODON to state CDS */ 
        temp = GenomeWise9_DC_SHADOW_MATRIX(mat,i-0,j-3,START_CODON) + GNE_CDS(mat->evi,i,mat->gen,j);   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for CDS */ 
        /* Add any movement independant score and put away */ 
         GenomeWise9_DC_SHADOW_MATRIX(mat,i,j,CDS) = score;  
        /* Finished calculating state CDS */ 


        /* For state CDS_INTRON_0 */ 
        /* setting first movement to score */ 
        score = GenomeWise9_DC_SHADOW_MATRIX(mat,i-0,j-1,CDS_INTRON_0) + GNE_CDS_INTRON(mat->evi,i,mat->gen,j);  
        /* From state CDS to state CDS_INTRON_0 */ 
        temp = GenomeWise9_DC_SHADOW_MATRIX(mat,i-0,j-8,CDS) + GNE_CDS_5SS(mat->evi,i,mat->gen,(j-7),0);     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for CDS_INTRON_0 */ 
        /* Add any movement independant score and put away */ 
         GenomeWise9_DC_SHADOW_MATRIX(mat,i,j,CDS_INTRON_0) = score; 
        /* Finished calculating state CDS_INTRON_0 */ 


        /* For state CDS_INTRON_1 */ 
        /* setting first movement to score */ 
        score = GenomeWise9_DC_SHADOW_MATRIX(mat,i-0,j-1,CDS_INTRON_1) + GNE_CDS_INTRON(mat->evi,i,mat->gen,j);  
        /* From state CDS to state CDS_INTRON_1 */ 
        temp = GenomeWise9_DC_SHADOW_MATRIX(mat,i-0,j-9,CDS) + GNE_CDS_5SS(mat->evi,i,mat->gen,(j-7),1);     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for CDS_INTRON_1 */ 
        /* Add any movement independant score and put away */ 
         GenomeWise9_DC_SHADOW_MATRIX(mat,i,j,CDS_INTRON_1) = score; 
        /* Finished calculating state CDS_INTRON_1 */ 


        /* For state CDS_INTRON_2 */ 
        /* setting first movement to score */ 
        score = GenomeWise9_DC_SHADOW_MATRIX(mat,i-0,j-1,CDS_INTRON_2) + GNE_CDS_INTRON(mat->evi,i,mat->gen,j);  
        /* From state CDS to state CDS_INTRON_2 */ 
        temp = GenomeWise9_DC_SHADOW_MATRIX(mat,i-0,j-10,CDS) + GNE_CDS_5SS(mat->evi,i,mat->gen,(j-7),2);    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for CDS_INTRON_2 */ 
        /* Add any movement independant score and put away */ 
         GenomeWise9_DC_SHADOW_MATRIX(mat,i,j,CDS_INTRON_2) = score; 
        /* Finished calculating state CDS_INTRON_2 */ 


        /* For state STOP_CODON */ 
        /* setting first movement to score */ 
        score = GenomeWise9_DC_SHADOW_MATRIX(mat,i-0,j-3,CDS) + GNE_STOP_CODON(mat->evi,i,mat->gen,j);   


        /* Ok - finished max calculation for STOP_CODON */ 
        /* Add any movement independant score and put away */ 
         GenomeWise9_DC_SHADOW_MATRIX(mat,i,j,STOP_CODON) = score;   
        /* Finished calculating state STOP_CODON */ 


        /* For state UTR3 */ 
        /* setting first movement to score */ 
        score = GenomeWise9_DC_SHADOW_MATRIX(mat,i-0,j-1,UTR3) + GNE_UTR(mat->evi,i,mat->gen,j);     
        /* From state CDS to state UTR3 */ 
        temp = GenomeWise9_DC_SHADOW_MATRIX(mat,i-0,j-1,CDS) + (GNE_UTR(mat->evi,i,mat->gen,j)+mat->non_stop_codon);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state STOP_CODON to state UTR3 */ 
        temp = GenomeWise9_DC_SHADOW_MATRIX(mat,i-0,j-1,STOP_CODON) + GNE_UTR(mat->evi,i,mat->gen,j);    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state UTR3_INTRON to state UTR3 */ 
        temp = GenomeWise9_DC_SHADOW_MATRIX(mat,i-0,j-1,UTR3_INTRON) + GNE_UTR_3SS(mat->evi,i,mat->gen,j);   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for UTR3 */ 
        /* Add any movement independant score and put away */ 
         GenomeWise9_DC_SHADOW_MATRIX(mat,i,j,UTR3) = score; 
        /* Finished calculating state UTR3 */ 


        /* For state UTR3_INTRON */ 
        /* setting first movement to score */ 
        score = GenomeWise9_DC_SHADOW_MATRIX(mat,i-0,j-1,UTR3_INTRON) + GNE_UTR_INTRON(mat->evi,i,mat->gen,j);   
        /* From state UTR3 to state UTR3_INTRON */ 
        temp = GenomeWise9_DC_SHADOW_MATRIX(mat,i-0,j-1,UTR3) + GNE_UTR_5SS(mat->evi,i,mat->gen,j);  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for UTR3_INTRON */ 
        /* Add any movement independant score and put away */ 
         GenomeWise9_DC_SHADOW_MATRIX(mat,i,j,UTR3_INTRON) = score;  
        /* Finished calculating state UTR3_INTRON */ 
        } /* end of this is strip */ 
      } /* end of for each valid j column */ 


/* Function:  init_dc_GenomeWise9(mat)
 *
 * Descrip: No Description
 *
 * Arg:        mat [UNKN ] Undocumented argument [GenomeWise9 *]
 *
 */
}    
void init_dc_GenomeWise9(GenomeWise9 * mat) 
{
    register int i;  
    register int j;  
    register int k;  


    for(j=0;j<12;j++)    {  
      for(i=(-0);i<mat->evi->len;i++)    {  
        GenomeWise9_DC_SHADOW_MATRIX(mat,i,j,UTR5) = NEGI;   
        GenomeWise9_DC_SHADOW_MATRIX(mat,i,j,UTR5_INTRON) = NEGI;    
        GenomeWise9_DC_SHADOW_MATRIX(mat,i,j,START_CODON) = NEGI;    
        GenomeWise9_DC_SHADOW_MATRIX(mat,i,j,CDS) = NEGI;    
        GenomeWise9_DC_SHADOW_MATRIX(mat,i,j,CDS_INTRON_0) = NEGI;   
        GenomeWise9_DC_SHADOW_MATRIX(mat,i,j,CDS_INTRON_1) = NEGI;   
        GenomeWise9_DC_SHADOW_MATRIX(mat,i,j,CDS_INTRON_2) = NEGI;   
        GenomeWise9_DC_SHADOW_MATRIX(mat,i,j,STOP_CODON) = NEGI; 
        GenomeWise9_DC_SHADOW_MATRIX(mat,i,j,UTR3) = NEGI;   
        GenomeWise9_DC_SHADOW_MATRIX(mat,i,j,UTR3_INTRON) = NEGI;    
        for(k=0;k<7;k++) {  
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR5,k) = (-1);    
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR5_INTRON,k) = (-1); 
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,START_CODON,k) = (-1); 
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS,k) = (-1); 
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS_INTRON_0,k) = (-1);    
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS_INTRON_1,k) = (-1);    
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS_INTRON_2,k) = (-1);    
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,STOP_CODON,k) = (-1);  
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR3,k) = (-1);    
          GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR3_INTRON,k) = (-1); 
          }  
        }  
      }  


    return;  
}    


/* Function:  start_end_find_end_GenomeWise9(mat,endj)
 *
 * Descrip:    First function used to find end of the best path in the special state !end
 *
 *
 * Arg:         mat [UNKN ] Matrix in small mode [GenomeWise9 *]
 * Arg:        endj [WRITE] position of end in j (meaningless in i) [int *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int start_end_find_end_GenomeWise9(GenomeWise9 * mat,int * endj) 
{
    register int j;  
    register int max;    
    register int maxj;   


    max = GenomeWise9_DC_SHADOW_SPECIAL(mat,0,mat->gen->seq->len-1,END); 
    maxj = mat->gen->seq->len-1;     
    for(j= mat->gen->seq->len-2 ;j >= 0 ;j--)    {  
      if( GenomeWise9_DC_SHADOW_SPECIAL(mat,0,j,END) > max ) {  
        max = GenomeWise9_DC_SHADOW_SPECIAL(mat,0,j,END);    
        maxj = j;    
        }  
      }  


    if( endj != NULL)    
      *endj = maxj;  


    return max;  
}    


/* Function:  dc_optimised_start_end_calc_GenomeWise9(*mat,dpenv)
 *
 * Descrip:    Calculates special strip, leaving start/end/score points in shadow matrix
 *             Works off specially laid out memory from steve searle
 *
 *
 * Arg:         *mat [UNKN ] Undocumented argument [GenomeWise9]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean dc_optimised_start_end_calc_GenomeWise9(GenomeWise9 *mat,DPEnvelope * dpenv) 
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
    leni = mat->evi->len;    
    lenj = mat->gen->seq->len;   
    total = leni * lenj; 


    score_pointers = (int *) calloc (10 * (leni + 0) * 10,sizeof(int));  
    shadow_pointers = (int *) calloc (10 * (leni + 0) * 10 * 8,sizeof(int)); 


    for(j=0;j<lenj;j++)  { /*for each j strip*/ 
      for(i=0;i<leni;i++)    { /*for each i position in strip*/ 
        num++;   
        if( dpenv != NULL && is_in_DPEnvelope(dpenv,i,j) == FALSE )  { /*Is not in envelope*/ 
          GenomeWise9_DC_OPT_SHADOW_MATRIX(mat,i,j,UTR5) = NEGI;     
          GenomeWise9_DC_OPT_SHADOW_MATRIX(mat,i,j,UTR5_INTRON) = NEGI;  
          GenomeWise9_DC_OPT_SHADOW_MATRIX(mat,i,j,START_CODON) = NEGI;  
          GenomeWise9_DC_OPT_SHADOW_MATRIX(mat,i,j,CDS) = NEGI;  
          GenomeWise9_DC_OPT_SHADOW_MATRIX(mat,i,j,CDS_INTRON_0) = NEGI;     
          GenomeWise9_DC_OPT_SHADOW_MATRIX(mat,i,j,CDS_INTRON_1) = NEGI;     
          GenomeWise9_DC_OPT_SHADOW_MATRIX(mat,i,j,CDS_INTRON_2) = NEGI;     
          GenomeWise9_DC_OPT_SHADOW_MATRIX(mat,i,j,STOP_CODON) = NEGI;   
          GenomeWise9_DC_OPT_SHADOW_MATRIX(mat,i,j,UTR3) = NEGI;     
          GenomeWise9_DC_OPT_SHADOW_MATRIX(mat,i,j,UTR3_INTRON) = NEGI;  
          continue;  
          } /* end of Is not in envelope */ 
        if( num%1000 == 0)   
          log_full_error(REPORT,0,"%6d Cells done [%2d%%%%]",num,num*100/total); 




        /* For state UTR5 */ 
        /* setting first movement to score */ 
        score = GenomeWise9_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,UTR5) + GNE_UTR(mat->evi,i,mat->gen,j) + (0);   
        /* assign local shadown pointer */ 
        localsp = &(GenomeWise9_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,UTR5,0));    
        /* From state UTR5_INTRON to state UTR5 */ 
        temp = GenomeWise9_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,UTR5_INTRON) + GNE_UTR_3SS(mat->evi,i,mat->gen,j) +(0);  
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GenomeWise9_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,UTR5_INTRON,0));   
          }  
        /* From state SPECIAL_UTR5 to state UTR5 */ 
        temp = GenomeWise9_DC_OPT_SHADOW_SPECIAL(mat,i-0,j-1,SPECIAL_UTR5) + GNE_UTR(mat->evi,i,mat->gen,j) + (0);   
        if( temp  > score )  {  
          score = temp;  
          /* This state [SPECIAL_UTR5] is a special for UTR5... push top shadow pointers here */ 
          localshadow[0]= i; 
          localshadow[1]= j; 
          localshadow[2]= UTR5;  
          localshadow[3]= (-1);  
          localshadow[4]= (-1);  
          localshadow[5]= (-1);  
          localshadow[6]= score; 
          localsp = localshadow; 
          }  
        /* From state START to state UTR5 */ 
        temp = GenomeWise9_DC_OPT_SHADOW_SPECIAL(mat,i-0,j-1,START) + GNE_UTR5_START(mat->evi,i,mat->gen,j) + (0);   
        if( temp  > score )  {  
          score = temp;  
          /* This state [START] is a special for UTR5... push top shadow pointers here */ 
          localshadow[0]= i; 
          localshadow[1]= j; 
          localshadow[2]= UTR5;  
          localshadow[3]= (-1);  
          localshadow[4]= (-1);  
          localshadow[5]= (-1);  
          localshadow[6]= score; 
          localsp = localshadow; 
          }  
        /* From state INTERGENIC to state UTR5 */ 
        temp = GenomeWise9_DC_OPT_SHADOW_SPECIAL(mat,i-0,j-1,INTERGENIC) + (GNE_UTR(mat->evi,i,mat->gen,j)+GNE_UTR5_START(mat->evi,i,mat->gen,j)) + (0);     
        if( temp  > score )  {  
          score = temp;  
          /* This state [INTERGENIC] is a special for UTR5... push top shadow pointers here */ 
          localshadow[0]= i; 
          localshadow[1]= j; 
          localshadow[2]= UTR5;  
          localshadow[3]= (-1);  
          localshadow[4]= (-1);  
          localshadow[5]= (-1);  
          localshadow[6]= score; 
          localsp = localshadow; 
          }  
        /* From state PREGENE_INTERGENIC to state UTR5 */ 
        temp = GenomeWise9_DC_OPT_SHADOW_SPECIAL(mat,i-0,j-1,PREGENE_INTERGENIC) + GNE_UTR(mat->evi,i,mat->gen,j) + (0);     
        if( temp  > score )  {  
          score = temp;  
          /* This state [PREGENE_INTERGENIC] is a special for UTR5... push top shadow pointers here */ 
          localshadow[0]= i; 
          localshadow[1]= j; 
          localshadow[2]= UTR5;  
          localshadow[3]= (-1);  
          localshadow[4]= (-1);  
          localshadow[5]= (-1);  
          localshadow[6]= score; 
          localsp = localshadow; 
          }  


        /* Ok - finished max calculation for UTR5 */ 
        /* Add any movement independant score and put away */ 
        /* Actually, already done inside scores */ 
         GenomeWise9_DC_OPT_SHADOW_MATRIX(mat,i,j,UTR5) = score; 
        for(k=0;k<7;k++) 
          GenomeWise9_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,UTR5,k) = localsp[k];  
        /* Now figure out if any specials need this score */ 


        /* state UTR5 is a source for special SPECIAL_UTR5 */ 
        temp = score + (mat->switchcost) + (0) ;     
        if( temp > GenomeWise9_DC_OPT_SHADOW_SPECIAL(mat,i,j,SPECIAL_UTR5) )     {  
          GenomeWise9_DC_OPT_SHADOW_SPECIAL(mat,i,j,SPECIAL_UTR5) = temp;    
          /* Have to push only bottem half of system here */ 
          for(k=0;k<3;k++)   
            GenomeWise9_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,SPECIAL_UTR5,k) = GenomeWise9_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,UTR5,k);  
          GenomeWise9_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,SPECIAL_UTR5,6) = GenomeWise9_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,UTR5,6);    
          GenomeWise9_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,SPECIAL_UTR5,3) = i;  
          GenomeWise9_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,SPECIAL_UTR5,4) = j;  
          GenomeWise9_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,SPECIAL_UTR5,5) = UTR5;   
          }  




        /* Finished calculating state UTR5 */ 


        /* For state UTR5_INTRON */ 
        /* setting first movement to score */ 
        score = GenomeWise9_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,UTR5_INTRON) + GNE_UTR_INTRON(mat->evi,i,mat->gen,j) + (0);     
        /* assign local shadown pointer */ 
        localsp = &(GenomeWise9_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,UTR5_INTRON,0)); 
        /* From state UTR5 to state UTR5_INTRON */ 
        temp = GenomeWise9_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,UTR5) + GNE_UTR_5SS(mat->evi,i,mat->gen,j) +(0);     
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GenomeWise9_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,UTR5,0));  
          }  


        /* Ok - finished max calculation for UTR5_INTRON */ 
        /* Add any movement independant score and put away */ 
        /* Actually, already done inside scores */ 
         GenomeWise9_DC_OPT_SHADOW_MATRIX(mat,i,j,UTR5_INTRON) = score;  
        for(k=0;k<7;k++) 
          GenomeWise9_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,UTR5_INTRON,k) = localsp[k];   
        /* Now figure out if any specials need this score */ 


        /* Finished calculating state UTR5_INTRON */ 


        /* For state START_CODON */ 
        /* setting first movement to score */ 
        score = GenomeWise9_DC_OPT_SHADOW_MATRIX(mat,i-0,j-3,UTR5_INTRON) + GNE_START_CODON(mat->evi,i,mat->gen,j) + (0);    
        /* assign local shadown pointer */ 
        localsp = &(GenomeWise9_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 3,UTR5_INTRON,0)); 
        /* From state UTR5 to state START_CODON */ 
        temp = GenomeWise9_DC_OPT_SHADOW_MATRIX(mat,i-0,j-3,UTR5) + GNE_START_CODON(mat->evi,i,mat->gen,j) +(0);     
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GenomeWise9_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 3,UTR5,0));  
          }  


        /* Ok - finished max calculation for START_CODON */ 
        /* Add any movement independant score and put away */ 
        /* Actually, already done inside scores */ 
         GenomeWise9_DC_OPT_SHADOW_MATRIX(mat,i,j,START_CODON) = score;  
        for(k=0;k<7;k++) 
          GenomeWise9_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,START_CODON,k) = localsp[k];   
        /* Now figure out if any specials need this score */ 


        /* Finished calculating state START_CODON */ 


        /* For state CDS */ 
        /* setting first movement to score */ 
        score = GenomeWise9_DC_OPT_SHADOW_MATRIX(mat,i-0,j-3,CDS) + GNE_CDS(mat->evi,i,mat->gen,j) + (0);    
        /* assign local shadown pointer */ 
        localsp = &(GenomeWise9_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 3,CDS,0)); 
        /* From state CDS_INTRON_0 to state CDS */ 
        temp = GenomeWise9_DC_OPT_SHADOW_MATRIX(mat,i-0,j-6,CDS_INTRON_0) + GNE_CDS_3SS(mat->evi,i,mat->gen,(j-3),0) +(0);   
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GenomeWise9_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 6,CDS_INTRON_0,0));  
          }  
        /* From state CDS_INTRON_1 to state CDS */ 
        temp = GenomeWise9_DC_OPT_SHADOW_MATRIX(mat,i-0,j-5,CDS_INTRON_1) + GNE_CDS_3SS(mat->evi,i,mat->gen,(j-2),1) +(0);   
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GenomeWise9_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 5,CDS_INTRON_1,0));  
          }  
        /* From state CDS_INTRON_2 to state CDS */ 
        temp = GenomeWise9_DC_OPT_SHADOW_MATRIX(mat,i-0,j-4,CDS_INTRON_2) + GNE_CDS_3SS(mat->evi,i,mat->gen,(j-1),2) +(0);   
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GenomeWise9_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 4,CDS_INTRON_2,0));  
          }  
        /* From state CDS to state CDS */ 
        temp = GenomeWise9_DC_OPT_SHADOW_MATRIX(mat,i-0,j-2,CDS) + GNE_CDS_FRAMESHIFT(mat->evi,i,mat->gen,j,2) +(0);     
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GenomeWise9_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 2,CDS,0));   
          }  
        /* From state CDS to state CDS */ 
        temp = GenomeWise9_DC_OPT_SHADOW_MATRIX(mat,i-0,j-4,CDS) + GNE_CDS_FRAMESHIFT(mat->evi,i,mat->gen,j,4) +(0);     
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GenomeWise9_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 4,CDS,0));   
          }  
        /* From state UTR5 to state CDS */ 
        temp = GenomeWise9_DC_OPT_SHADOW_MATRIX(mat,i-0,j-3,UTR5) + (GNE_CDS(mat->evi,i,mat->gen,j)+mat->non_start_codon) +(0);  
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GenomeWise9_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 3,UTR5,0));  
          }  
        /* From state START_CODON to state CDS */ 
        temp = GenomeWise9_DC_OPT_SHADOW_MATRIX(mat,i-0,j-3,START_CODON) + GNE_CDS(mat->evi,i,mat->gen,j) +(0);  
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GenomeWise9_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 3,START_CODON,0));   
          }  
        /* From state SPECIAL_CDS to state CDS */ 
        temp = GenomeWise9_DC_OPT_SHADOW_SPECIAL(mat,i-0,j-3,SPECIAL_CDS) + GNE_CDS(mat->evi,i,mat->gen,j) + (0);    
        if( temp  > score )  {  
          score = temp;  
          /* This state [SPECIAL_CDS] is a special for CDS... push top shadow pointers here */ 
          localshadow[0]= i; 
          localshadow[1]= j; 
          localshadow[2]= CDS;   
          localshadow[3]= (-1);  
          localshadow[4]= (-1);  
          localshadow[5]= (-1);  
          localshadow[6]= score; 
          localsp = localshadow; 
          }  
        /* From state INTERGENIC to state CDS */ 
        temp = GenomeWise9_DC_OPT_SHADOW_SPECIAL(mat,i-0,j-3,INTERGENIC) + (GNE_CDS(mat->evi,i,mat->gen,j)+GNE_UTR5_START(mat->evi,i,mat->gen,j)) + (0);     
        if( temp  > score )  {  
          score = temp;  
          /* This state [INTERGENIC] is a special for CDS... push top shadow pointers here */ 
          localshadow[0]= i; 
          localshadow[1]= j; 
          localshadow[2]= CDS;   
          localshadow[3]= (-1);  
          localshadow[4]= (-1);  
          localshadow[5]= (-1);  
          localshadow[6]= score; 
          localsp = localshadow; 
          }  


        /* Ok - finished max calculation for CDS */ 
        /* Add any movement independant score and put away */ 
        /* Actually, already done inside scores */ 
         GenomeWise9_DC_OPT_SHADOW_MATRIX(mat,i,j,CDS) = score;  
        for(k=0;k<7;k++) 
          GenomeWise9_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,CDS,k) = localsp[k];   
        /* Now figure out if any specials need this score */ 


        /* state CDS is a source for special INTERGENIC */ 
        temp = score + ((mat->newgenecost+GNE_UTR3_END(mat->evi,i,mat->gen,j))) + (0) ;  
        if( temp > GenomeWise9_DC_OPT_SHADOW_SPECIAL(mat,i,j,INTERGENIC) )   {  
          GenomeWise9_DC_OPT_SHADOW_SPECIAL(mat,i,j,INTERGENIC) = temp;  
          /* Have to push only bottem half of system here */ 
          for(k=0;k<3;k++)   
            GenomeWise9_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,INTERGENIC,k) = GenomeWise9_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,CDS,k); 
          GenomeWise9_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,INTERGENIC,6) = GenomeWise9_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,CDS,6);   
          GenomeWise9_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,INTERGENIC,3) = i;    
          GenomeWise9_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,INTERGENIC,4) = j;    
          GenomeWise9_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,INTERGENIC,5) = CDS;  
          }  




        /* state CDS is a source for special SPECIAL_CDS */ 
        temp = score + ((mat->switchcost+mat->rndcodon->codon[CSEQ_GENOMIC_CODON(mat->gen,j)])) + (0) ;  
        if( temp > GenomeWise9_DC_OPT_SHADOW_SPECIAL(mat,i,j,SPECIAL_CDS) )  {  
          GenomeWise9_DC_OPT_SHADOW_SPECIAL(mat,i,j,SPECIAL_CDS) = temp;     
          /* Have to push only bottem half of system here */ 
          for(k=0;k<3;k++)   
            GenomeWise9_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,SPECIAL_CDS,k) = GenomeWise9_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,CDS,k);    
          GenomeWise9_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,SPECIAL_CDS,6) = GenomeWise9_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,CDS,6);  
          GenomeWise9_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,SPECIAL_CDS,3) = i;   
          GenomeWise9_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,SPECIAL_CDS,4) = j;   
          GenomeWise9_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,SPECIAL_CDS,5) = CDS; 
          }  




        /* Finished calculating state CDS */ 


        /* For state CDS_INTRON_0 */ 
        /* setting first movement to score */ 
        score = GenomeWise9_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,CDS_INTRON_0) + GNE_CDS_INTRON(mat->evi,i,mat->gen,j) + (0);    
        /* assign local shadown pointer */ 
        localsp = &(GenomeWise9_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,CDS_INTRON_0,0));    
        /* From state CDS to state CDS_INTRON_0 */ 
        temp = GenomeWise9_DC_OPT_SHADOW_MATRIX(mat,i-0,j-8,CDS) + GNE_CDS_5SS(mat->evi,i,mat->gen,(j-7),0) +(0);    
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GenomeWise9_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 8,CDS,0));   
          }  


        /* Ok - finished max calculation for CDS_INTRON_0 */ 
        /* Add any movement independant score and put away */ 
        /* Actually, already done inside scores */ 
         GenomeWise9_DC_OPT_SHADOW_MATRIX(mat,i,j,CDS_INTRON_0) = score; 
        for(k=0;k<7;k++) 
          GenomeWise9_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,CDS_INTRON_0,k) = localsp[k];  
        /* Now figure out if any specials need this score */ 


        /* Finished calculating state CDS_INTRON_0 */ 


        /* For state CDS_INTRON_1 */ 
        /* setting first movement to score */ 
        score = GenomeWise9_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,CDS_INTRON_1) + GNE_CDS_INTRON(mat->evi,i,mat->gen,j) + (0);    
        /* assign local shadown pointer */ 
        localsp = &(GenomeWise9_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,CDS_INTRON_1,0));    
        /* From state CDS to state CDS_INTRON_1 */ 
        temp = GenomeWise9_DC_OPT_SHADOW_MATRIX(mat,i-0,j-9,CDS) + GNE_CDS_5SS(mat->evi,i,mat->gen,(j-7),1) +(0);    
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GenomeWise9_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 9,CDS,0));   
          }  


        /* Ok - finished max calculation for CDS_INTRON_1 */ 
        /* Add any movement independant score and put away */ 
        /* Actually, already done inside scores */ 
         GenomeWise9_DC_OPT_SHADOW_MATRIX(mat,i,j,CDS_INTRON_1) = score; 
        for(k=0;k<7;k++) 
          GenomeWise9_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,CDS_INTRON_1,k) = localsp[k];  
        /* Now figure out if any specials need this score */ 


        /* Finished calculating state CDS_INTRON_1 */ 


        /* For state CDS_INTRON_2 */ 
        /* setting first movement to score */ 
        score = GenomeWise9_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,CDS_INTRON_2) + GNE_CDS_INTRON(mat->evi,i,mat->gen,j) + (0);    
        /* assign local shadown pointer */ 
        localsp = &(GenomeWise9_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,CDS_INTRON_2,0));    
        /* From state CDS to state CDS_INTRON_2 */ 
        temp = GenomeWise9_DC_OPT_SHADOW_MATRIX(mat,i-0,j-10,CDS) + GNE_CDS_5SS(mat->evi,i,mat->gen,(j-7),2) +(0);   
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GenomeWise9_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 10,CDS,0));  
          }  


        /* Ok - finished max calculation for CDS_INTRON_2 */ 
        /* Add any movement independant score and put away */ 
        /* Actually, already done inside scores */ 
         GenomeWise9_DC_OPT_SHADOW_MATRIX(mat,i,j,CDS_INTRON_2) = score; 
        for(k=0;k<7;k++) 
          GenomeWise9_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,CDS_INTRON_2,k) = localsp[k];  
        /* Now figure out if any specials need this score */ 


        /* Finished calculating state CDS_INTRON_2 */ 


        /* For state STOP_CODON */ 
        /* setting first movement to score */ 
        score = GenomeWise9_DC_OPT_SHADOW_MATRIX(mat,i-0,j-3,CDS) + GNE_STOP_CODON(mat->evi,i,mat->gen,j) + (0);     
        /* assign local shadown pointer */ 
        localsp = &(GenomeWise9_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 3,CDS,0)); 


        /* Ok - finished max calculation for STOP_CODON */ 
        /* Add any movement independant score and put away */ 
        /* Actually, already done inside scores */ 
         GenomeWise9_DC_OPT_SHADOW_MATRIX(mat,i,j,STOP_CODON) = score;   
        for(k=0;k<7;k++) 
          GenomeWise9_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,STOP_CODON,k) = localsp[k];    
        /* Now figure out if any specials need this score */ 


        /* Finished calculating state STOP_CODON */ 


        /* For state UTR3 */ 
        /* setting first movement to score */ 
        score = GenomeWise9_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,UTR3) + GNE_UTR(mat->evi,i,mat->gen,j) + (0);   
        /* assign local shadown pointer */ 
        localsp = &(GenomeWise9_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,UTR3,0));    
        /* From state CDS to state UTR3 */ 
        temp = GenomeWise9_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,CDS) + (GNE_UTR(mat->evi,i,mat->gen,j)+mat->non_stop_codon) +(0);    
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GenomeWise9_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,CDS,0));   
          }  
        /* From state STOP_CODON to state UTR3 */ 
        temp = GenomeWise9_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,STOP_CODON) + GNE_UTR(mat->evi,i,mat->gen,j) +(0);   
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GenomeWise9_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,STOP_CODON,0));    
          }  
        /* From state UTR3_INTRON to state UTR3 */ 
        temp = GenomeWise9_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,UTR3_INTRON) + GNE_UTR_3SS(mat->evi,i,mat->gen,j) +(0);  
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GenomeWise9_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,UTR3_INTRON,0));   
          }  
        /* From state INTERGENIC to state UTR3 */ 
        temp = GenomeWise9_DC_OPT_SHADOW_SPECIAL(mat,i-0,j-1,INTERGENIC) + GNE_UTR(mat->evi,i,mat->gen,j) + (0);     
        if( temp  > score )  {  
          score = temp;  
          /* This state [INTERGENIC] is a special for UTR3... push top shadow pointers here */ 
          localshadow[0]= i; 
          localshadow[1]= j; 
          localshadow[2]= UTR3;  
          localshadow[3]= (-1);  
          localshadow[4]= (-1);  
          localshadow[5]= (-1);  
          localshadow[6]= score; 
          localsp = localshadow; 
          }  


        /* Ok - finished max calculation for UTR3 */ 
        /* Add any movement independant score and put away */ 
        /* Actually, already done inside scores */ 
         GenomeWise9_DC_OPT_SHADOW_MATRIX(mat,i,j,UTR3) = score; 
        for(k=0;k<7;k++) 
          GenomeWise9_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,UTR3,k) = localsp[k];  
        /* Now figure out if any specials need this score */ 


        /* state UTR3 is a source for special POSTGENE_INTERGENIC */ 
        temp = score + ((mat->newgenecost+GNE_UTR3_END(mat->evi,i,mat->gen,j))) + (0) ;  
        if( temp > GenomeWise9_DC_OPT_SHADOW_SPECIAL(mat,i,j,POSTGENE_INTERGENIC) )  {  
          GenomeWise9_DC_OPT_SHADOW_SPECIAL(mat,i,j,POSTGENE_INTERGENIC) = temp;     
          /* Have to push only bottem half of system here */ 
          for(k=0;k<3;k++)   
            GenomeWise9_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,POSTGENE_INTERGENIC,k) = GenomeWise9_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,UTR3,k);   
          GenomeWise9_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,POSTGENE_INTERGENIC,6) = GenomeWise9_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,UTR3,6); 
          GenomeWise9_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,POSTGENE_INTERGENIC,3) = i;   
          GenomeWise9_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,POSTGENE_INTERGENIC,4) = j;   
          GenomeWise9_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,POSTGENE_INTERGENIC,5) = UTR3;    
          }  




        /* state UTR3 is a source for special INTERGENIC */ 
        temp = score + ((mat->newgenecost+GNE_UTR3_END(mat->evi,i,mat->gen,j))) + (0) ;  
        if( temp > GenomeWise9_DC_OPT_SHADOW_SPECIAL(mat,i,j,INTERGENIC) )   {  
          GenomeWise9_DC_OPT_SHADOW_SPECIAL(mat,i,j,INTERGENIC) = temp;  
          /* Have to push only bottem half of system here */ 
          for(k=0;k<3;k++)   
            GenomeWise9_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,INTERGENIC,k) = GenomeWise9_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,UTR3,k);    
          GenomeWise9_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,INTERGENIC,6) = GenomeWise9_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,UTR3,6);  
          GenomeWise9_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,INTERGENIC,3) = i;    
          GenomeWise9_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,INTERGENIC,4) = j;    
          GenomeWise9_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,INTERGENIC,5) = UTR3; 
          }  




        /* state UTR3 is a source for special SPECIAL_UTR3 */ 
        temp = score + (mat->switchcost) + (0) ;     
        if( temp > GenomeWise9_DC_OPT_SHADOW_SPECIAL(mat,i,j,SPECIAL_UTR3) )     {  
          GenomeWise9_DC_OPT_SHADOW_SPECIAL(mat,i,j,SPECIAL_UTR3) = temp;    
          /* Have to push only bottem half of system here */ 
          for(k=0;k<3;k++)   
            GenomeWise9_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,SPECIAL_UTR3,k) = GenomeWise9_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,UTR3,k);  
          GenomeWise9_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,SPECIAL_UTR3,6) = GenomeWise9_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,UTR3,6);    
          GenomeWise9_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,SPECIAL_UTR3,3) = i;  
          GenomeWise9_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,SPECIAL_UTR3,4) = j;  
          GenomeWise9_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,SPECIAL_UTR3,5) = UTR3;   
          }  




        /* Finished calculating state UTR3 */ 


        /* For state UTR3_INTRON */ 
        /* setting first movement to score */ 
        score = GenomeWise9_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,UTR3_INTRON) + GNE_UTR_INTRON(mat->evi,i,mat->gen,j) + (0);     
        /* assign local shadown pointer */ 
        localsp = &(GenomeWise9_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,UTR3_INTRON,0)); 
        /* From state UTR3 to state UTR3_INTRON */ 
        temp = GenomeWise9_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,UTR3) + GNE_UTR_5SS(mat->evi,i,mat->gen,j) +(0);     
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GenomeWise9_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,UTR3,0));  
          }  
        /* From state SPECIAL_UTR3 to state UTR3_INTRON */ 
        temp = GenomeWise9_DC_OPT_SHADOW_SPECIAL(mat,i-0,j-1,SPECIAL_UTR3) + GNE_UTR(mat->evi,i,mat->gen,j) + (0);   
        if( temp  > score )  {  
          score = temp;  
          /* This state [SPECIAL_UTR3] is a special for UTR3_INTRON... push top shadow pointers here */ 
          localshadow[0]= i; 
          localshadow[1]= j; 
          localshadow[2]= UTR3_INTRON;   
          localshadow[3]= (-1);  
          localshadow[4]= (-1);  
          localshadow[5]= (-1);  
          localshadow[6]= score; 
          localsp = localshadow; 
          }  


        /* Ok - finished max calculation for UTR3_INTRON */ 
        /* Add any movement independant score and put away */ 
        /* Actually, already done inside scores */ 
         GenomeWise9_DC_OPT_SHADOW_MATRIX(mat,i,j,UTR3_INTRON) = score;  
        for(k=0;k<7;k++) 
          GenomeWise9_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,UTR3_INTRON,k) = localsp[k];   
        /* Now figure out if any specials need this score */ 


        /* Finished calculating state UTR3_INTRON */ 


        } /* end of for each i position in strip */ 


      /* Special state PREGENE_INTERGENIC has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GenomeWise9_DC_OPT_SHADOW_SPECIAL(mat,0,j,PREGENE_INTERGENIC); 


      /* Source START is a special source for PREGENE_INTERGENIC */ 
      temp = GenomeWise9_DC_OPT_SHADOW_SPECIAL(mat,0,j - 1,START) + (0) + (0);   
      if( temp > score ) {  
        score = temp;    
        /* Also got to propagate shadows  */ 
        for(k=0;k<7;k++) 
          GenomeWise9_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,PREGENE_INTERGENIC,k) = GenomeWise9_DC_OPT_SHADOW_SPECIAL_SP(mat,i - 0,j - 1,START,k);    
        }  


      /* Put back score... (now updated!) */ 
      GenomeWise9_DC_OPT_SHADOW_SPECIAL(mat,0,j,PREGENE_INTERGENIC) = score; 
      /* Finished updating state PREGENE_INTERGENIC */ 




      /* Special state POSTGENE_INTERGENIC has no special to special movements */ 


      /* Special state INTERGENIC has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GenomeWise9_DC_OPT_SHADOW_SPECIAL(mat,0,j,INTERGENIC); 


      /* Source INTERGENIC is a special source for INTERGENIC */ 
      temp = GenomeWise9_DC_OPT_SHADOW_SPECIAL(mat,0,j - 1,INTERGENIC) + (0) + (0);  
      if( temp > score ) {  
        score = temp;    
        /* Also got to propagate shadows  */ 
        for(k=0;k<7;k++) 
          GenomeWise9_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,INTERGENIC,k) = GenomeWise9_DC_OPT_SHADOW_SPECIAL_SP(mat,i - 0,j - 1,INTERGENIC,k);   
        }  


      /* Source CDS for state INTERGENIC is not special... already calculated */ 
      /* Source UTR3 for state INTERGENIC is not special... already calculated */ 
      /* Put back score... (now updated!) */ 
      GenomeWise9_DC_OPT_SHADOW_SPECIAL(mat,0,j,INTERGENIC) = score; 
      /* Finished updating state INTERGENIC */ 




      /* Special state SPECIAL_UTR5 has no special to special movements */ 


      /* Special state SPECIAL_UTR3 has no special to special movements */ 


      /* Special state SPECIAL_CDS has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GenomeWise9_DC_OPT_SHADOW_SPECIAL(mat,0,j,SPECIAL_CDS);    


      /* Source CDS for state SPECIAL_CDS is not special... already calculated */ 
      /* Source SPECIAL_CDS is a special source for SPECIAL_CDS */ 
      temp = GenomeWise9_DC_OPT_SHADOW_SPECIAL(mat,0,j - 3,SPECIAL_CDS) + (mat->rndcodon->codon[CSEQ_GENOMIC_CODON(mat->gen,j)]) + (0);  
      if( temp > score ) {  
        score = temp;    
        /* Also got to propagate shadows  */ 
        for(k=0;k<7;k++) 
          GenomeWise9_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,SPECIAL_CDS,k) = GenomeWise9_DC_OPT_SHADOW_SPECIAL_SP(mat,i - 0,j - 3,SPECIAL_CDS,k); 
        }  


      /* Put back score... (now updated!) */ 
      GenomeWise9_DC_OPT_SHADOW_SPECIAL(mat,0,j,SPECIAL_CDS) = score;    
      /* Finished updating state SPECIAL_CDS */ 




      /* Special state START has no special to special movements */ 


      /* Special state END has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GenomeWise9_DC_OPT_SHADOW_SPECIAL(mat,0,j,END);    


      /* Source INTERGENIC is a special source for END */ 
      temp = GenomeWise9_DC_OPT_SHADOW_SPECIAL(mat,0,j - 1,INTERGENIC) + (0) + (0);  
      if( temp > score ) {  
        score = temp;    
        /* Also got to propagate shadows  */ 
        for(k=0;k<7;k++) 
          GenomeWise9_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,k) = GenomeWise9_DC_OPT_SHADOW_SPECIAL_SP(mat,i - 0,j - 1,INTERGENIC,k);  
        }  


      /* Put back score... (now updated!) */ 
      GenomeWise9_DC_OPT_SHADOW_SPECIAL(mat,0,j,END) = score;    
      /* Finished updating state END */ 


      } /* end of for each j strip */ 
    free(score_pointers);    
    free(shadow_pointers);   
    return TRUE;     
}    


/* Function:  init_start_end_linear_GenomeWise9(mat)
 *
 * Descrip: No Description
 *
 * Arg:        mat [UNKN ] Undocumented argument [GenomeWise9 *]
 *
 */
void init_start_end_linear_GenomeWise9(GenomeWise9 * mat) 
{
    register int i;  
    register int j;  
    for(j=0;j<12;j++)    {  
      for(i=(-0);i<mat->evi->len;i++)    {  
        GenomeWise9_DC_SHADOW_MATRIX(mat,i,j,UTR5) = NEGI;   
        GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR5,0) = (-1);  
        GenomeWise9_DC_SHADOW_MATRIX(mat,i,j,UTR5_INTRON) = NEGI;    
        GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR5_INTRON,0) = (-1);   
        GenomeWise9_DC_SHADOW_MATRIX(mat,i,j,START_CODON) = NEGI;    
        GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,START_CODON,0) = (-1);   
        GenomeWise9_DC_SHADOW_MATRIX(mat,i,j,CDS) = NEGI;    
        GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS,0) = (-1);   
        GenomeWise9_DC_SHADOW_MATRIX(mat,i,j,CDS_INTRON_0) = NEGI;   
        GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS_INTRON_0,0) = (-1);  
        GenomeWise9_DC_SHADOW_MATRIX(mat,i,j,CDS_INTRON_1) = NEGI;   
        GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS_INTRON_1,0) = (-1);  
        GenomeWise9_DC_SHADOW_MATRIX(mat,i,j,CDS_INTRON_2) = NEGI;   
        GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,CDS_INTRON_2,0) = (-1);  
        GenomeWise9_DC_SHADOW_MATRIX(mat,i,j,STOP_CODON) = NEGI; 
        GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,STOP_CODON,0) = (-1);    
        GenomeWise9_DC_SHADOW_MATRIX(mat,i,j,UTR3) = NEGI;   
        GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR3,0) = (-1);  
        GenomeWise9_DC_SHADOW_MATRIX(mat,i,j,UTR3_INTRON) = NEGI;    
        GenomeWise9_DC_SHADOW_MATRIX_SP(mat,i,j,UTR3_INTRON,0) = (-1);   
        }  
      }  


    for(j=(-10);j<mat->gen->seq->len;j++)    {  
      GenomeWise9_DC_SHADOW_SPECIAL(mat,0,j,PREGENE_INTERGENIC) = NEGI;  
      GenomeWise9_DC_SHADOW_SPECIAL_SP(mat,0,j,PREGENE_INTERGENIC,0) = (-1); 
      GenomeWise9_DC_SHADOW_SPECIAL(mat,0,j,POSTGENE_INTERGENIC) = NEGI; 
      GenomeWise9_DC_SHADOW_SPECIAL_SP(mat,0,j,POSTGENE_INTERGENIC,0) = (-1);    
      GenomeWise9_DC_SHADOW_SPECIAL(mat,0,j,INTERGENIC) = NEGI;  
      GenomeWise9_DC_SHADOW_SPECIAL_SP(mat,0,j,INTERGENIC,0) = (-1); 
      GenomeWise9_DC_SHADOW_SPECIAL(mat,0,j,SPECIAL_UTR5) = NEGI;    
      GenomeWise9_DC_SHADOW_SPECIAL_SP(mat,0,j,SPECIAL_UTR5,0) = (-1);   
      GenomeWise9_DC_SHADOW_SPECIAL(mat,0,j,SPECIAL_UTR3) = NEGI;    
      GenomeWise9_DC_SHADOW_SPECIAL_SP(mat,0,j,SPECIAL_UTR3,0) = (-1);   
      GenomeWise9_DC_SHADOW_SPECIAL(mat,0,j,SPECIAL_CDS) = NEGI; 
      GenomeWise9_DC_SHADOW_SPECIAL_SP(mat,0,j,SPECIAL_CDS,0) = (-1);    
      GenomeWise9_DC_SHADOW_SPECIAL(mat,0,j,START) = 0;  
      GenomeWise9_DC_SHADOW_SPECIAL_SP(mat,0,j,START,0) = j; 
      GenomeWise9_DC_SHADOW_SPECIAL(mat,0,j,END) = NEGI; 
      GenomeWise9_DC_SHADOW_SPECIAL_SP(mat,0,j,END,0) = (-1);    
      }  


    return;  
}    


/* Function:  convert_PackAln_to_AlnBlock_GenomeWise9(pal)
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
AlnBlock * convert_PackAln_to_AlnBlock_GenomeWise9(PackAln * pal) 
{
    AlnConvertSet * acs; 
    AlnBlock * alb;  


    acs = AlnConvertSet_GenomeWise9();   
    alb = AlnBlock_from_PackAln(acs,pal);    
    free_AlnConvertSet(acs); 
    return alb;  
}    


 static char * query_label[] = { "UTR5","START_CODON","CDS","CDS_INTRON","STOP_CODON","UTR3","INTERGENIC","SPECIAL","END" }; 
/* Function:  AlnConvertSet_GenomeWise9(void)
 *
 * Descrip: No Description
 *
 *
 * Return [UNKN ]  Undocumented return value [AlnConvertSet *]
 *
 */
 static char * target_label[] = { "UTR5","UTR5_INTRON","CODON","3SS_PHASE_0","3SS_PHASE_1","3SS_PHASE_2","SEQUENCE_DELETION","CDS_INTRON","5SS_PHASE_0","5SS_PHASE_1","5SS_PHASE_2","STOP_CODON","UTR3","UTR3_INTRON","RANDOM_DNA","END" };  
AlnConvertSet * AlnConvertSet_GenomeWise9(void) 
{
    AlnConvertUnit * acu;    
    AlnConvertSet  * out;    


    out = AlnConvertSet_alloc_std(); 


    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = UTR5;  
    acu->state2 = UTR5;  
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[0];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = UTR5_INTRON;   
    acu->state2 = UTR5;  
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[0];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = SPECIAL_UTR5 + 10; 
    acu->is_from_special = TRUE; 
    acu->state2 = UTR5;  
    acu->offi = (-1);    
    acu->offj = 1;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[0];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = START + 10;    
    acu->is_from_special = TRUE; 
    acu->state2 = UTR5;  
    acu->offi = (-1);    
    acu->offj = 1;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[0];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = INTERGENIC + 10;   
    acu->is_from_special = TRUE; 
    acu->state2 = UTR5;  
    acu->offi = (-1);    
    acu->offj = 1;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[0];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = PREGENE_INTERGENIC + 10;   
    acu->is_from_special = TRUE; 
    acu->state2 = UTR5;  
    acu->offi = (-1);    
    acu->offj = 1;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[0];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = UTR5_INTRON;   
    acu->state2 = UTR5_INTRON;   
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[1];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = UTR5;  
    acu->state2 = UTR5_INTRON;   
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[1];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = UTR5_INTRON;   
    acu->state2 = START_CODON;   
    acu->offi = 0;   
    acu->offj = 3;   
    acu->label1 = query_label[1];    
    acu->label2 = target_label[2];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = UTR5;  
    acu->state2 = START_CODON;   
    acu->offi = 0;   
    acu->offj = 3;   
    acu->label1 = query_label[1];    
    acu->label2 = target_label[2];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = CDS;   
    acu->state2 = CDS;   
    acu->offi = 0;   
    acu->offj = 3;   
    acu->label1 = query_label[2];    
    acu->label2 = target_label[2];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = CDS_INTRON_0;  
    acu->state2 = CDS;   
    acu->offi = 0;   
    acu->offj = 6;   
    acu->label1 = query_label[2];    
    acu->label2 = target_label[3];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = CDS_INTRON_1;  
    acu->state2 = CDS;   
    acu->offi = 0;   
    acu->offj = 5;   
    acu->label1 = query_label[2];    
    acu->label2 = target_label[4];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = CDS_INTRON_2;  
    acu->state2 = CDS;   
    acu->offi = 0;   
    acu->offj = 4;   
    acu->label1 = query_label[2];    
    acu->label2 = target_label[5];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = CDS;   
    acu->state2 = CDS;   
    acu->offi = 0;   
    acu->offj = 2;   
    acu->label1 = query_label[2];    
    acu->label2 = target_label[6];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = CDS;   
    acu->state2 = CDS;   
    acu->offi = 0;   
    acu->offj = 4;   
    acu->label1 = query_label[2];    
    acu->label2 = target_label[6];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = UTR5;  
    acu->state2 = CDS;   
    acu->offi = 0;   
    acu->offj = 3;   
    acu->label1 = query_label[2];    
    acu->label2 = target_label[2];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = START_CODON;   
    acu->state2 = CDS;   
    acu->offi = 0;   
    acu->offj = 3;   
    acu->label1 = query_label[2];    
    acu->label2 = target_label[2];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = SPECIAL_CDS + 10;  
    acu->is_from_special = TRUE; 
    acu->state2 = CDS;   
    acu->offi = (-1);    
    acu->offj = 3;   
    acu->label1 = query_label[2];    
    acu->label2 = target_label[2];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = INTERGENIC + 10;   
    acu->is_from_special = TRUE; 
    acu->state2 = CDS;   
    acu->offi = (-1);    
    acu->offj = 3;   
    acu->label1 = query_label[2];    
    acu->label2 = target_label[2];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = CDS_INTRON_0;  
    acu->state2 = CDS_INTRON_0;  
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[3];    
    acu->label2 = target_label[7];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = CDS;   
    acu->state2 = CDS_INTRON_0;  
    acu->offi = 0;   
    acu->offj = 8;   
    acu->label1 = query_label[3];    
    acu->label2 = target_label[8];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = CDS_INTRON_1;  
    acu->state2 = CDS_INTRON_1;  
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[3];    
    acu->label2 = target_label[7];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = CDS;   
    acu->state2 = CDS_INTRON_1;  
    acu->offi = 0;   
    acu->offj = 9;   
    acu->label1 = query_label[3];    
    acu->label2 = target_label[9];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = CDS_INTRON_2;  
    acu->state2 = CDS_INTRON_2;  
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[3];    
    acu->label2 = target_label[7];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = CDS;   
    acu->state2 = CDS_INTRON_2;  
    acu->offi = 0;   
    acu->offj = 10;  
    acu->label1 = query_label[3];    
    acu->label2 = target_label[10];  
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = CDS;   
    acu->state2 = STOP_CODON;    
    acu->offi = 0;   
    acu->offj = 3;   
    acu->label1 = query_label[4];    
    acu->label2 = target_label[11];  
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = UTR3;  
    acu->state2 = UTR3;  
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[5];    
    acu->label2 = target_label[12];  
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = CDS;   
    acu->state2 = UTR3;  
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[5];    
    acu->label2 = target_label[12];  
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = STOP_CODON;    
    acu->state2 = UTR3;  
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[5];    
    acu->label2 = target_label[12];  
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = UTR3_INTRON;   
    acu->state2 = UTR3;  
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[5];    
    acu->label2 = target_label[12];  
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = INTERGENIC + 10;   
    acu->is_from_special = TRUE; 
    acu->state2 = UTR3;  
    acu->offi = (-1);    
    acu->offj = 1;   
    acu->label1 = query_label[5];    
    acu->label2 = target_label[12];  
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = UTR3_INTRON;   
    acu->state2 = UTR3_INTRON;   
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[5];    
    acu->label2 = target_label[13];  
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = UTR3;  
    acu->state2 = UTR3_INTRON;   
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[5];    
    acu->label2 = target_label[13];  
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = SPECIAL_UTR3 + 10; 
    acu->is_from_special = TRUE; 
    acu->state2 = UTR3_INTRON;   
    acu->offi = (-1);    
    acu->offj = 1;   
    acu->label1 = query_label[5];    
    acu->label2 = target_label[13];  
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = START + 10;    
    acu->state2 = PREGENE_INTERGENIC + 10;   
    acu->offi = (-1);    
    acu->offj = 1;   
    acu->label1 = query_label[6];    
    acu->label2 = target_label[14];  
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = UTR3;  
    acu->state2 = POSTGENE_INTERGENIC + 10;  
    acu->offi = (-1);    
    acu->offj = 0;   
    acu->label1 = query_label[6];    
    acu->label2 = target_label[14];  
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = INTERGENIC + 10;   
    acu->state2 = INTERGENIC + 10;   
    acu->offi = (-1);    
    acu->offj = 1;   
    acu->label1 = query_label[6];    
    acu->label2 = target_label[14];  
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = CDS;   
    acu->state2 = INTERGENIC + 10;   
    acu->offi = (-1);    
    acu->offj = 0;   
    acu->label1 = query_label[6];    
    acu->label2 = target_label[14];  
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = UTR3;  
    acu->state2 = INTERGENIC + 10;   
    acu->offi = (-1);    
    acu->offj = 0;   
    acu->label1 = query_label[6];    
    acu->label2 = target_label[14];  
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = UTR5;  
    acu->state2 = SPECIAL_UTR5 + 10;     
    acu->offi = (-1);    
    acu->offj = 0;   
    acu->label1 = query_label[7];    
    acu->label2 = target_label[0];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = UTR3;  
    acu->state2 = SPECIAL_UTR3 + 10;     
    acu->offi = (-1);    
    acu->offj = 0;   
    acu->label1 = query_label[7];    
    acu->label2 = target_label[12];  
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = CDS;   
    acu->state2 = SPECIAL_CDS + 10;  
    acu->offi = (-1);    
    acu->offj = 0;   
    acu->label1 = query_label[7];    
    acu->label2 = target_label[2];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = SPECIAL_CDS + 10;  
    acu->state2 = SPECIAL_CDS + 10;  
    acu->offi = (-1);    
    acu->offj = 3;   
    acu->label1 = query_label[7];    
    acu->label2 = target_label[2];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = INTERGENIC + 10;   
    acu->state2 = END + 10;  
    acu->offi = (-1);    
    acu->offj = 1;   
    acu->label1 = query_label[8];    
    acu->label2 = target_label[15];  
    add_collapse_label_AlnConvertSet(out,"UTR5","UTR5"); 
    add_collapse_label_AlnConvertSet(out,"UTR3","UTR3"); 
    add_collapse_label_AlnConvertSet(out,"UTR5","UTR5_INTRON");  
    add_collapse_label_AlnConvertSet(out,"UTR3","UTR3_INTRON");  
    add_collapse_label_AlnConvertSet(out,"INTERGENIC","RANDOM_DNA"); 
    add_collapse_label_AlnConvertSet(out,"CDS_INTRON","CDS_INTRON"); 
    return out;  
}    


/* Function:  PackAln_read_Expl_GenomeWise9(mat)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [GenomeWise9 *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * PackAln_read_Expl_GenomeWise9(GenomeWise9 * mat) 
{
    GenomeWise9_access_func_holder holder;   


    holder.access_main    = GenomeWise9_explicit_access_main;    
    holder.access_special = GenomeWise9_explicit_access_special; 
    return PackAln_read_generic_GenomeWise9(mat,holder); 
}    


/* Function:  GenomeWise9_explicit_access_main(mat,i,j,state)
 *
 * Descrip: No Description
 *
 * Arg:          mat [UNKN ] Undocumented argument [GenomeWise9 *]
 * Arg:            i [UNKN ] Undocumented argument [int]
 * Arg:            j [UNKN ] Undocumented argument [int]
 * Arg:        state [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int GenomeWise9_explicit_access_main(GenomeWise9 * mat,int i,int j,int state) 
{
    return GenomeWise9_EXPL_MATRIX(mat,i,j,state);   
}    


/* Function:  GenomeWise9_explicit_access_special(mat,i,j,state)
 *
 * Descrip: No Description
 *
 * Arg:          mat [UNKN ] Undocumented argument [GenomeWise9 *]
 * Arg:            i [UNKN ] Undocumented argument [int]
 * Arg:            j [UNKN ] Undocumented argument [int]
 * Arg:        state [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int GenomeWise9_explicit_access_special(GenomeWise9 * mat,int i,int j,int state) 
{
    return GenomeWise9_EXPL_SPECIAL(mat,i,j,state);  
}    


/* Function:  PackAln_read_generic_GenomeWise9(mat,h)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [GenomeWise9 *]
 * Arg:          h [UNKN ] Undocumented argument [GenomeWise9_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * PackAln_read_generic_GenomeWise9(GenomeWise9 * mat,GenomeWise9_access_func_holder h) 
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


    out->score =  find_end_GenomeWise9(mat,&i,&j,&state,&isspecial,h);   


    /* Add final end transition (at the moment we have not got the score! */ 
    if( (pau= PackAlnUnit_alloc()) == NULL  || add_PackAln(out,pau) == FALSE )   {  
      warn("Failed the first PackAlnUnit alloc, %d length of Alignment in GenomeWise9_basic_read, returning a mess.(Sorry!)",out->len);  
      return out;    
      }  


    /* Put in positions for end trans. Remember that coordinates in C style */ 
    pau->i = i;  
    pau->j = j;  
    if( isspecial != TRUE)   
      pau->state = state;    
    else pau->state = state + 10;    
    prev=pau;    
    while( state != START || isspecial != TRUE)  { /*while state != START*/ 


      if( isspecial == TRUE )    
        max_calc_special_GenomeWise9(mat,i,j,state,isspecial,&i,&j,&state,&isspecial,&cellscore,h);  
      else   
        max_calc_GenomeWise9(mat,i,j,state,isspecial,&i,&j,&state,&isspecial,&cellscore,h);  
      if(i == GenomeWise9_READ_OFF_ERROR || j == GenomeWise9_READ_OFF_ERROR || state == GenomeWise9_READ_OFF_ERROR ) {  
        warn("Problem - hit bad read off system, exiting now");  
        break;   
        }  
      if( (pau= PackAlnUnit_alloc()) == NULL  || add_PackAln(out,pau) == FALSE ) {  
        warn("Failed a PackAlnUnit alloc, %d length of Alignment in GenomeWise9_basic_read, returning partial alignment",out->len);  
        break;   
        }  


      /* Put in positions for block. Remember that coordinates in C style */ 
      pau->i = i;    
      pau->j = j;    
      if( isspecial != TRUE)     
        pau->state = state;  
      else pau->state = state + 10;  
      prev->score = cellscore;   
      prev = pau;    
      } /* end of while state != START */ 


    invert_PackAln(out); 
    return out;  
}    


/* Function:  find_end_GenomeWise9(mat,ri,rj,state,isspecial,h)
 *
 * Descrip: No Description
 *
 * Arg:              mat [UNKN ] Undocumented argument [GenomeWise9 *]
 * Arg:               ri [UNKN ] Undocumented argument [int *]
 * Arg:               rj [UNKN ] Undocumented argument [int *]
 * Arg:            state [UNKN ] Undocumented argument [int *]
 * Arg:        isspecial [UNKN ] Undocumented argument [boolean *]
 * Arg:                h [UNKN ] Undocumented argument [GenomeWise9_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int find_end_GenomeWise9(GenomeWise9 * mat,int * ri,int * rj,int * state,boolean * isspecial,GenomeWise9_access_func_holder h) 
{
    int j;   
    int max; 
    int maxj;    
    int temp;    


    max = (*h.access_special)(mat,0,mat->gen->seq->len-1,END);   
    maxj = mat->gen->seq->len-1;     
    for(j= mat->gen->seq->len-2 ;j >= 0 ;j--)    {  
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


/* Function:  GenomeWise9_debug_show_matrix(mat,starti,stopi,startj,stopj,ofp)
 *
 * Descrip: No Description
 *
 * Arg:           mat [UNKN ] Undocumented argument [GenomeWise9 *]
 * Arg:        starti [UNKN ] Undocumented argument [int]
 * Arg:         stopi [UNKN ] Undocumented argument [int]
 * Arg:        startj [UNKN ] Undocumented argument [int]
 * Arg:         stopj [UNKN ] Undocumented argument [int]
 * Arg:           ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
void GenomeWise9_debug_show_matrix(GenomeWise9 * mat,int starti,int stopi,int startj,int stopj,FILE * ofp) 
{
    register int i;  
    register int j;  


    for(i=starti;i<stopi && i < mat->evi->len;i++)   {  
      for(j=startj;j<stopj && j < mat->gen->seq->len;j++)    {  
        fprintf(ofp,"Cell [%d - %d]\n",i,j);     
        fprintf(ofp,"State UTR5 %d\n",GenomeWise9_EXPL_MATRIX(mat,i,j,UTR5));    
        fprintf(ofp,"State UTR5_INTRON %d\n",GenomeWise9_EXPL_MATRIX(mat,i,j,UTR5_INTRON));  
        fprintf(ofp,"State START_CODON %d\n",GenomeWise9_EXPL_MATRIX(mat,i,j,START_CODON));  
        fprintf(ofp,"State CDS %d\n",GenomeWise9_EXPL_MATRIX(mat,i,j,CDS));  
        fprintf(ofp,"State CDS_INTRON_0 %d\n",GenomeWise9_EXPL_MATRIX(mat,i,j,CDS_INTRON_0));    
        fprintf(ofp,"State CDS_INTRON_1 %d\n",GenomeWise9_EXPL_MATRIX(mat,i,j,CDS_INTRON_1));    
        fprintf(ofp,"State CDS_INTRON_2 %d\n",GenomeWise9_EXPL_MATRIX(mat,i,j,CDS_INTRON_2));    
        fprintf(ofp,"State STOP_CODON %d\n",GenomeWise9_EXPL_MATRIX(mat,i,j,STOP_CODON));    
        fprintf(ofp,"State UTR3 %d\n",GenomeWise9_EXPL_MATRIX(mat,i,j,UTR3));    
        fprintf(ofp,"State UTR3_INTRON %d\n",GenomeWise9_EXPL_MATRIX(mat,i,j,UTR3_INTRON));  
        fprintf(ofp,"\n\n"); 
        }  
      }  


}    


/* Function:  max_calc_GenomeWise9(mat,i,j,state,isspecial,reti,retj,retstate,retspecial,cellscore,h)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [GenomeWise9 *]
 * Arg:                 i [UNKN ] Undocumented argument [int]
 * Arg:                 j [UNKN ] Undocumented argument [int]
 * Arg:             state [UNKN ] Undocumented argument [int]
 * Arg:         isspecial [UNKN ] Undocumented argument [boolean]
 * Arg:              reti [UNKN ] Undocumented argument [int *]
 * Arg:              retj [UNKN ] Undocumented argument [int *]
 * Arg:          retstate [UNKN ] Undocumented argument [int *]
 * Arg:        retspecial [UNKN ] Undocumented argument [boolean *]
 * Arg:         cellscore [UNKN ] Undocumented argument [int *]
 * Arg:                 h [UNKN ] Undocumented argument [GenomeWise9_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int max_calc_GenomeWise9(GenomeWise9 * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,GenomeWise9_access_func_holder h) 
{
    register int temp;   
    register int cscore; 


    *reti = (*retj) = (*retstate) = GenomeWise9_READ_OFF_ERROR;  


    if( i < 0 || j < 0 || i > mat->evi->len || j > mat->gen->seq->len)   {  
      warn("In GenomeWise9 matrix special read off - out of bounds on matrix [i,j is %d,%d state %d in standard matrix]",i,j,state); 
      return -1;     
      }  


    /* Then you have to select the correct switch statement to figure out the readoff      */ 
    /* Somewhat odd - reverse the order of calculation and return as soon as it is correct */ 
    cscore = (*h.access_main)(mat,i,j,state);    
    switch(state)    { /*Switch state */ 
      case UTR5 :    
        temp = cscore - (GNE_UTR(mat->evi,i,mat->gen,j)) -  (0); 
        if( temp == (*h.access_special)(mat,i - 0,j - 1,PREGENE_INTERGENIC) )    {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = PREGENE_INTERGENIC;    
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_special)(mat,i-0,j-1,PREGENE_INTERGENIC);   
            }  
          return (*h.access_main)(mat,i - 0,j - 1,PREGENE_INTERGENIC);   
          }  
        temp = cscore - ((GNE_UTR(mat->evi,i,mat->gen,j)+GNE_UTR5_START(mat->evi,i,mat->gen,j))) -  (0); 
        if( temp == (*h.access_special)(mat,i - 0,j - 1,INTERGENIC) )    {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = INTERGENIC;    
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_special)(mat,i-0,j-1,INTERGENIC);   
            }  
          return (*h.access_main)(mat,i - 0,j - 1,INTERGENIC);   
          }  
        temp = cscore - (GNE_UTR5_START(mat->evi,i,mat->gen,j)) -  (0);  
        if( temp == (*h.access_special)(mat,i - 0,j - 1,START) ) {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = START; 
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_special)(mat,i-0,j-1,START);    
            }  
          return (*h.access_main)(mat,i - 0,j - 1,START);    
          }  
        temp = cscore - (GNE_UTR(mat->evi,i,mat->gen,j)) -  (0); 
        if( temp == (*h.access_special)(mat,i - 0,j - 1,SPECIAL_UTR5) )  {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = SPECIAL_UTR5;  
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_special)(mat,i-0,j-1,SPECIAL_UTR5); 
            }  
          return (*h.access_main)(mat,i - 0,j - 1,SPECIAL_UTR5);     
          }  
        temp = cscore - (GNE_UTR_3SS(mat->evi,i,mat->gen,j)) -  (0); 
        if( temp == (*h.access_main)(mat,i - 0,j - 1,UTR5_INTRON) )  {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = UTR5_INTRON;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,UTR5_INTRON); 
            }  
          return (*h.access_main)(mat,i - 0,j - 1,UTR5_INTRON);  
          }  
        temp = cscore - (GNE_UTR(mat->evi,i,mat->gen,j)) -  (0); 
        if( temp == (*h.access_main)(mat,i - 0,j - 1,UTR5) ) {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = UTR5;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,UTR5);    
            }  
          return (*h.access_main)(mat,i - 0,j - 1,UTR5);     
          }  
        warn("Major problem (!) - in GenomeWise9 read off, position %d,%d state %d no source found!",i,j,state); 
        return (-1); 
      case UTR5_INTRON :     
        temp = cscore - (GNE_UTR_5SS(mat->evi,i,mat->gen,j)) -  (0); 
        if( temp == (*h.access_main)(mat,i - 0,j - 1,UTR5) ) {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = UTR5;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,UTR5);    
            }  
          return (*h.access_main)(mat,i - 0,j - 1,UTR5);     
          }  
        temp = cscore - (GNE_UTR_INTRON(mat->evi,i,mat->gen,j)) -  (0);  
        if( temp == (*h.access_main)(mat,i - 0,j - 1,UTR5_INTRON) )  {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = UTR5_INTRON;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,UTR5_INTRON); 
            }  
          return (*h.access_main)(mat,i - 0,j - 1,UTR5_INTRON);  
          }  
        warn("Major problem (!) - in GenomeWise9 read off, position %d,%d state %d no source found!",i,j,state); 
        return (-1); 
      case START_CODON :     
        temp = cscore - (GNE_START_CODON(mat->evi,i,mat->gen,j)) -  (0); 
        if( temp == (*h.access_main)(mat,i - 0,j - 3,UTR5) ) {  
          *reti = i - 0; 
          *retj = j - 3; 
          *retstate = UTR5;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-3,UTR5);    
            }  
          return (*h.access_main)(mat,i - 0,j - 3,UTR5);     
          }  
        temp = cscore - (GNE_START_CODON(mat->evi,i,mat->gen,j)) -  (0); 
        if( temp == (*h.access_main)(mat,i - 0,j - 3,UTR5_INTRON) )  {  
          *reti = i - 0; 
          *retj = j - 3; 
          *retstate = UTR5_INTRON;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-3,UTR5_INTRON); 
            }  
          return (*h.access_main)(mat,i - 0,j - 3,UTR5_INTRON);  
          }  
        warn("Major problem (!) - in GenomeWise9 read off, position %d,%d state %d no source found!",i,j,state); 
        return (-1); 
      case CDS :     
        temp = cscore - ((GNE_CDS(mat->evi,i,mat->gen,j)+GNE_UTR5_START(mat->evi,i,mat->gen,j))) -  (0); 
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
        temp = cscore - (GNE_CDS(mat->evi,i,mat->gen,j)) -  (0); 
        if( temp == (*h.access_special)(mat,i - 0,j - 3,SPECIAL_CDS) )   {  
          *reti = i - 0; 
          *retj = j - 3; 
          *retstate = SPECIAL_CDS;   
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_special)(mat,i-0,j-3,SPECIAL_CDS);  
            }  
          return (*h.access_main)(mat,i - 0,j - 3,SPECIAL_CDS);  
          }  
        temp = cscore - (GNE_CDS(mat->evi,i,mat->gen,j)) -  (0); 
        if( temp == (*h.access_main)(mat,i - 0,j - 3,START_CODON) )  {  
          *reti = i - 0; 
          *retj = j - 3; 
          *retstate = START_CODON;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-3,START_CODON); 
            }  
          return (*h.access_main)(mat,i - 0,j - 3,START_CODON);  
          }  
        temp = cscore - ((GNE_CDS(mat->evi,i,mat->gen,j)+mat->non_start_codon)) -  (0);  
        if( temp == (*h.access_main)(mat,i - 0,j - 3,UTR5) ) {  
          *reti = i - 0; 
          *retj = j - 3; 
          *retstate = UTR5;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-3,UTR5);    
            }  
          return (*h.access_main)(mat,i - 0,j - 3,UTR5);     
          }  
        temp = cscore - (GNE_CDS_FRAMESHIFT(mat->evi,i,mat->gen,j,4)) -  (0);    
        if( temp == (*h.access_main)(mat,i - 0,j - 4,CDS) )  {  
          *reti = i - 0; 
          *retj = j - 4; 
          *retstate = CDS;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-4,CDS); 
            }  
          return (*h.access_main)(mat,i - 0,j - 4,CDS);  
          }  
        temp = cscore - (GNE_CDS_FRAMESHIFT(mat->evi,i,mat->gen,j,2)) -  (0);    
        if( temp == (*h.access_main)(mat,i - 0,j - 2,CDS) )  {  
          *reti = i - 0; 
          *retj = j - 2; 
          *retstate = CDS;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-2,CDS); 
            }  
          return (*h.access_main)(mat,i - 0,j - 2,CDS);  
          }  
        temp = cscore - (GNE_CDS_3SS(mat->evi,i,mat->gen,(j-1),2)) -  (0);   
        if( temp == (*h.access_main)(mat,i - 0,j - 4,CDS_INTRON_2) ) {  
          *reti = i - 0; 
          *retj = j - 4; 
          *retstate = CDS_INTRON_2;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-4,CDS_INTRON_2);    
            }  
          return (*h.access_main)(mat,i - 0,j - 4,CDS_INTRON_2);     
          }  
        temp = cscore - (GNE_CDS_3SS(mat->evi,i,mat->gen,(j-2),1)) -  (0);   
        if( temp == (*h.access_main)(mat,i - 0,j - 5,CDS_INTRON_1) ) {  
          *reti = i - 0; 
          *retj = j - 5; 
          *retstate = CDS_INTRON_1;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-5,CDS_INTRON_1);    
            }  
          return (*h.access_main)(mat,i - 0,j - 5,CDS_INTRON_1);     
          }  
        temp = cscore - (GNE_CDS_3SS(mat->evi,i,mat->gen,(j-3),0)) -  (0);   
        if( temp == (*h.access_main)(mat,i - 0,j - 6,CDS_INTRON_0) ) {  
          *reti = i - 0; 
          *retj = j - 6; 
          *retstate = CDS_INTRON_0;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-6,CDS_INTRON_0);    
            }  
          return (*h.access_main)(mat,i - 0,j - 6,CDS_INTRON_0);     
          }  
        temp = cscore - (GNE_CDS(mat->evi,i,mat->gen,j)) -  (0); 
        if( temp == (*h.access_main)(mat,i - 0,j - 3,CDS) )  {  
          *reti = i - 0; 
          *retj = j - 3; 
          *retstate = CDS;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-3,CDS); 
            }  
          return (*h.access_main)(mat,i - 0,j - 3,CDS);  
          }  
        warn("Major problem (!) - in GenomeWise9 read off, position %d,%d state %d no source found!",i,j,state); 
        return (-1); 
      case CDS_INTRON_0 :    
        temp = cscore - (GNE_CDS_5SS(mat->evi,i,mat->gen,(j-7),0)) -  (0);   
        if( temp == (*h.access_main)(mat,i - 0,j - 8,CDS) )  {  
          *reti = i - 0; 
          *retj = j - 8; 
          *retstate = CDS;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-8,CDS); 
            }  
          return (*h.access_main)(mat,i - 0,j - 8,CDS);  
          }  
        temp = cscore - (GNE_CDS_INTRON(mat->evi,i,mat->gen,j)) -  (0);  
        if( temp == (*h.access_main)(mat,i - 0,j - 1,CDS_INTRON_0) ) {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = CDS_INTRON_0;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,CDS_INTRON_0);    
            }  
          return (*h.access_main)(mat,i - 0,j - 1,CDS_INTRON_0);     
          }  
        warn("Major problem (!) - in GenomeWise9 read off, position %d,%d state %d no source found!",i,j,state); 
        return (-1); 
      case CDS_INTRON_1 :    
        temp = cscore - (GNE_CDS_5SS(mat->evi,i,mat->gen,(j-7),1)) -  (0);   
        if( temp == (*h.access_main)(mat,i - 0,j - 9,CDS) )  {  
          *reti = i - 0; 
          *retj = j - 9; 
          *retstate = CDS;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-9,CDS); 
            }  
          return (*h.access_main)(mat,i - 0,j - 9,CDS);  
          }  
        temp = cscore - (GNE_CDS_INTRON(mat->evi,i,mat->gen,j)) -  (0);  
        if( temp == (*h.access_main)(mat,i - 0,j - 1,CDS_INTRON_1) ) {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = CDS_INTRON_1;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,CDS_INTRON_1);    
            }  
          return (*h.access_main)(mat,i - 0,j - 1,CDS_INTRON_1);     
          }  
        warn("Major problem (!) - in GenomeWise9 read off, position %d,%d state %d no source found!",i,j,state); 
        return (-1); 
      case CDS_INTRON_2 :    
        temp = cscore - (GNE_CDS_5SS(mat->evi,i,mat->gen,(j-7),2)) -  (0);   
        if( temp == (*h.access_main)(mat,i - 0,j - 10,CDS) ) {  
          *reti = i - 0; 
          *retj = j - 10;    
          *retstate = CDS;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-10,CDS);    
            }  
          return (*h.access_main)(mat,i - 0,j - 10,CDS);     
          }  
        temp = cscore - (GNE_CDS_INTRON(mat->evi,i,mat->gen,j)) -  (0);  
        if( temp == (*h.access_main)(mat,i - 0,j - 1,CDS_INTRON_2) ) {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = CDS_INTRON_2;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,CDS_INTRON_2);    
            }  
          return (*h.access_main)(mat,i - 0,j - 1,CDS_INTRON_2);     
          }  
        warn("Major problem (!) - in GenomeWise9 read off, position %d,%d state %d no source found!",i,j,state); 
        return (-1); 
      case STOP_CODON :  
        temp = cscore - (GNE_STOP_CODON(mat->evi,i,mat->gen,j)) -  (0);  
        if( temp == (*h.access_main)(mat,i - 0,j - 3,CDS) )  {  
          *reti = i - 0; 
          *retj = j - 3; 
          *retstate = CDS;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-3,CDS); 
            }  
          return (*h.access_main)(mat,i - 0,j - 3,CDS);  
          }  
        warn("Major problem (!) - in GenomeWise9 read off, position %d,%d state %d no source found!",i,j,state); 
        return (-1); 
      case UTR3 :    
        /* Has restricted position */ 
        if( (j-1) == 0  )    {  
          temp = cscore - (GNE_UTR(mat->evi,i,mat->gen,j)) -  (0);   
          if( temp == (*h.access_special)(mat,i - 0,j - 1,INTERGENIC) )  {  
            *reti = i - 0;   
            *retj = j - 1;   
            *retstate = INTERGENIC;  
            *retspecial = TRUE;  
            if( cellscore != NULL)   {  
              *cellscore = cscore - (*h.access_special)(mat,i-0,j-1,INTERGENIC); 
              }  
            return (*h.access_main)(mat,i - 0,j - 1,INTERGENIC);     
            }  
          }  
        temp = cscore - (GNE_UTR_3SS(mat->evi,i,mat->gen,j)) -  (0); 
        if( temp == (*h.access_main)(mat,i - 0,j - 1,UTR3_INTRON) )  {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = UTR3_INTRON;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,UTR3_INTRON); 
            }  
          return (*h.access_main)(mat,i - 0,j - 1,UTR3_INTRON);  
          }  
        temp = cscore - (GNE_UTR(mat->evi,i,mat->gen,j)) -  (0); 
        if( temp == (*h.access_main)(mat,i - 0,j - 1,STOP_CODON) )   {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = STOP_CODON;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,STOP_CODON);  
            }  
          return (*h.access_main)(mat,i - 0,j - 1,STOP_CODON);   
          }  
        temp = cscore - ((GNE_UTR(mat->evi,i,mat->gen,j)+mat->non_stop_codon)) -  (0);   
        if( temp == (*h.access_main)(mat,i - 0,j - 1,CDS) )  {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = CDS;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,CDS); 
            }  
          return (*h.access_main)(mat,i - 0,j - 1,CDS);  
          }  
        temp = cscore - (GNE_UTR(mat->evi,i,mat->gen,j)) -  (0); 
        if( temp == (*h.access_main)(mat,i - 0,j - 1,UTR3) ) {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = UTR3;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,UTR3);    
            }  
          return (*h.access_main)(mat,i - 0,j - 1,UTR3);     
          }  
        warn("Major problem (!) - in GenomeWise9 read off, position %d,%d state %d no source found!",i,j,state); 
        return (-1); 
      case UTR3_INTRON :     
        temp = cscore - (GNE_UTR(mat->evi,i,mat->gen,j)) -  (0); 
        if( temp == (*h.access_special)(mat,i - 0,j - 1,SPECIAL_UTR3) )  {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = SPECIAL_UTR3;  
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_special)(mat,i-0,j-1,SPECIAL_UTR3); 
            }  
          return (*h.access_main)(mat,i - 0,j - 1,SPECIAL_UTR3);     
          }  
        temp = cscore - (GNE_UTR_5SS(mat->evi,i,mat->gen,j)) -  (0); 
        if( temp == (*h.access_main)(mat,i - 0,j - 1,UTR3) ) {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = UTR3;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,UTR3);    
            }  
          return (*h.access_main)(mat,i - 0,j - 1,UTR3);     
          }  
        temp = cscore - (GNE_UTR_INTRON(mat->evi,i,mat->gen,j)) -  (0);  
        if( temp == (*h.access_main)(mat,i - 0,j - 1,UTR3_INTRON) )  {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = UTR3_INTRON;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,UTR3_INTRON); 
            }  
          return (*h.access_main)(mat,i - 0,j - 1,UTR3_INTRON);  
          }  
        warn("Major problem (!) - in GenomeWise9 read off, position %d,%d state %d no source found!",i,j,state); 
        return (-1); 
      default:   
        warn("Major problem (!) - in GenomeWise9 read off, position %d,%d state %d no source found!",i,j,state); 
        return (-1); 
      } /* end of Switch state  */ 
}    


/* Function:  max_calc_special_GenomeWise9(mat,i,j,state,isspecial,reti,retj,retstate,retspecial,cellscore,h)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [GenomeWise9 *]
 * Arg:                 i [UNKN ] Undocumented argument [int]
 * Arg:                 j [UNKN ] Undocumented argument [int]
 * Arg:             state [UNKN ] Undocumented argument [int]
 * Arg:         isspecial [UNKN ] Undocumented argument [boolean]
 * Arg:              reti [UNKN ] Undocumented argument [int *]
 * Arg:              retj [UNKN ] Undocumented argument [int *]
 * Arg:          retstate [UNKN ] Undocumented argument [int *]
 * Arg:        retspecial [UNKN ] Undocumented argument [boolean *]
 * Arg:         cellscore [UNKN ] Undocumented argument [int *]
 * Arg:                 h [UNKN ] Undocumented argument [GenomeWise9_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int max_calc_special_GenomeWise9(GenomeWise9 * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,GenomeWise9_access_func_holder h) 
{
    register int temp;   
    register int cscore; 


    *reti = (*retj) = (*retstate) = GenomeWise9_READ_OFF_ERROR;  


    if( j < 0 || j > mat->gen->seq->len) {  
      warn("In GenomeWise9 matrix special read off - out of bounds on matrix [j is %d in special]",j);   
      return -1;     
      }  


    cscore = (*h.access_special)(mat,i,j,state); 
    switch(state)    { /*switch on special states*/ 
      case PREGENE_INTERGENIC :  
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
      case POSTGENE_INTERGENIC :     
        /* source UTR3 is from main matrix */ 
        for(i= mat->evi->len-1;i >= 0 ;i--)  { /*for i >= 0*/ 
          temp = cscore - ((mat->newgenecost+GNE_UTR3_END(mat->evi,i,mat->gen,j))) - (0);    
          if( temp == (*h.access_main)(mat,i - 0,j - 0,UTR3) )   {  
            *reti = i - 0;   
            *retj = j - 0;   
            *retstate = UTR3;    
            *retspecial = FALSE; 
            if( cellscore != NULL)   {  
              *cellscore = cscore - (*h.access_main)(mat,i-0,j-0,UTR3);  
              }  
            return (*h.access_main)(mat,i - 0,j - 0,UTR3) ;  
            }  
          } /* end of for i >= 0 */ 
      case INTERGENIC :  
        /* source UTR3 is from main matrix */ 
        for(i= mat->evi->len-1;i >= 0 ;i--)  { /*for i >= 0*/ 
          temp = cscore - ((mat->newgenecost+GNE_UTR3_END(mat->evi,i,mat->gen,j))) - (0);    
          if( temp == (*h.access_main)(mat,i - 0,j - 0,UTR3) )   {  
            *reti = i - 0;   
            *retj = j - 0;   
            *retstate = UTR3;    
            *retspecial = FALSE; 
            if( cellscore != NULL)   {  
              *cellscore = cscore - (*h.access_main)(mat,i-0,j-0,UTR3);  
              }  
            return (*h.access_main)(mat,i - 0,j - 0,UTR3) ;  
            }  
          } /* end of for i >= 0 */ 
        /* source CDS is from main matrix */ 
        for(i= mat->evi->len-1;i >= 0 ;i--)  { /*for i >= 0*/ 
          temp = cscore - ((mat->newgenecost+GNE_UTR3_END(mat->evi,i,mat->gen,j))) - (0);    
          if( temp == (*h.access_main)(mat,i - 0,j - 0,CDS) )    {  
            *reti = i - 0;   
            *retj = j - 0;   
            *retstate = CDS; 
            *retspecial = FALSE; 
            if( cellscore != NULL)   {  
              *cellscore = cscore - (*h.access_main)(mat,i-0,j-0,CDS);   
              }  
            return (*h.access_main)(mat,i - 0,j - 0,CDS) ;   
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
      case SPECIAL_UTR5 :    
        /* source UTR5 is from main matrix */ 
        for(i= mat->evi->len-1;i >= 0 ;i--)  { /*for i >= 0*/ 
          temp = cscore - (mat->switchcost) - (0);   
          if( temp == (*h.access_main)(mat,i - 0,j - 0,UTR5) )   {  
            *reti = i - 0;   
            *retj = j - 0;   
            *retstate = UTR5;    
            *retspecial = FALSE; 
            if( cellscore != NULL)   {  
              *cellscore = cscore - (*h.access_main)(mat,i-0,j-0,UTR5);  
              }  
            return (*h.access_main)(mat,i - 0,j - 0,UTR5) ;  
            }  
          } /* end of for i >= 0 */ 
      case SPECIAL_UTR3 :    
        /* source UTR3 is from main matrix */ 
        for(i= mat->evi->len-1;i >= 0 ;i--)  { /*for i >= 0*/ 
          temp = cscore - (mat->switchcost) - (0);   
          if( temp == (*h.access_main)(mat,i - 0,j - 0,UTR3) )   {  
            *reti = i - 0;   
            *retj = j - 0;   
            *retstate = UTR3;    
            *retspecial = FALSE; 
            if( cellscore != NULL)   {  
              *cellscore = cscore - (*h.access_main)(mat,i-0,j-0,UTR3);  
              }  
            return (*h.access_main)(mat,i - 0,j - 0,UTR3) ;  
            }  
          } /* end of for i >= 0 */ 
      case SPECIAL_CDS :     
        /* source SPECIAL_CDS is a special */ 
        temp = cscore - (mat->rndcodon->codon[CSEQ_GENOMIC_CODON(mat->gen,j)]) - (0);    
        if( temp == (*h.access_special)(mat,i - 0,j - 3,SPECIAL_CDS) )   {  
          *reti = i - 0; 
          *retj = j - 3; 
          *retstate = SPECIAL_CDS;   
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_special)(mat,i-0,j-3,SPECIAL_CDS);  
            }  
          return (*h.access_special)(mat,i - 0,j - 3,SPECIAL_CDS) ;  
          }  
        /* source CDS is from main matrix */ 
        for(i= mat->evi->len-1;i >= 0 ;i--)  { /*for i >= 0*/ 
          temp = cscore - ((mat->switchcost+mat->rndcodon->codon[CSEQ_GENOMIC_CODON(mat->gen,j)])) - (0);    
          if( temp == (*h.access_main)(mat,i - 0,j - 0,CDS) )    {  
            *reti = i - 0;   
            *retj = j - 0;   
            *retstate = CDS; 
            *retspecial = FALSE; 
            if( cellscore != NULL)   {  
              *cellscore = cscore - (*h.access_main)(mat,i-0,j-0,CDS);   
              }  
            return (*h.access_main)(mat,i - 0,j - 0,CDS) ;   
            }  
          } /* end of for i >= 0 */ 
      case START :   
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
        warn("Major problem (!) - in GenomeWise9 read off, position %d,%d state %d no source found  dropped into default on source switch!",i,j,state);  
        return (-1); 
      } /* end of switch on special states */ 
}    


/* Function:  calculate_GenomeWise9(mat)
 *
 * Descrip:    This function calculates the GenomeWise9 matrix when in explicit mode
 *             To allocate the matrix use /allocate_Expl_GenomeWise9
 *
 *
 * Arg:        mat [UNKN ] GenomeWise9 which contains explicit basematrix memory [GenomeWise9 *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean calculate_GenomeWise9(GenomeWise9 * mat) 
{
    int i;   
    int j;   
    int leni;    
    int lenj;    
    int tot; 
    int num; 


    if( mat->basematrix->type != BASEMATRIX_TYPE_EXPLICIT )  {  
      warn("in calculate_GenomeWise9, passed a non Explicit matrix type, cannot calculate!");    
      return FALSE;  
      }  


    leni = mat->leni;    
    lenj = mat->lenj;    
    tot = leni * lenj;   
    num = 0; 


    start_reporting("GenomeWise9 Matrix calculation: "); 
    for(j=0;j<lenj;j++)  {  
      auto int score;    
      auto int temp;     
      for(i=0;i<leni;i++)    {  
        if( num%1000 == 0 )  
          log_full_error(REPORT,0,"[%7d] Cells %2d%%%%",num,num*100/tot);    
        num++;   


        /* For state UTR5 */ 
        /* setting first movement to score */ 
        score = GenomeWise9_EXPL_MATRIX(mat,i-0,j-1,UTR5) + GNE_UTR(mat->evi,i,mat->gen,j);  
        /* From state UTR5_INTRON to state UTR5 */ 
        temp = GenomeWise9_EXPL_MATRIX(mat,i-0,j-1,UTR5_INTRON) + GNE_UTR_3SS(mat->evi,i,mat->gen,j);    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SPECIAL_UTR5 to state UTR5 */ 
        temp = GenomeWise9_EXPL_SPECIAL(mat,i-0,j-1,SPECIAL_UTR5) + GNE_UTR(mat->evi,i,mat->gen,j);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state START to state UTR5 */ 
        temp = GenomeWise9_EXPL_SPECIAL(mat,i-0,j-1,START) + GNE_UTR5_START(mat->evi,i,mat->gen,j);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTERGENIC to state UTR5 */ 
        temp = GenomeWise9_EXPL_SPECIAL(mat,i-0,j-1,INTERGENIC) + (GNE_UTR(mat->evi,i,mat->gen,j)+GNE_UTR5_START(mat->evi,i,mat->gen,j));    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state PREGENE_INTERGENIC to state UTR5 */ 
        temp = GenomeWise9_EXPL_SPECIAL(mat,i-0,j-1,PREGENE_INTERGENIC) + GNE_UTR(mat->evi,i,mat->gen,j);    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for UTR5 */ 
        /* Add any movement independant score and put away */ 
         GenomeWise9_EXPL_MATRIX(mat,i,j,UTR5) = score;  


        /* state UTR5 is a source for special SPECIAL_UTR5 */ 
        temp = score + (mat->switchcost) + (0) ;     
        if( temp > GenomeWise9_EXPL_SPECIAL(mat,i,j,SPECIAL_UTR5) )  {  
          GenomeWise9_EXPL_SPECIAL(mat,i,j,SPECIAL_UTR5) = temp;     
          }  




        /* Finished calculating state UTR5 */ 


        /* For state UTR5_INTRON */ 
        /* setting first movement to score */ 
        score = GenomeWise9_EXPL_MATRIX(mat,i-0,j-1,UTR5_INTRON) + GNE_UTR_INTRON(mat->evi,i,mat->gen,j);    
        /* From state UTR5 to state UTR5_INTRON */ 
        temp = GenomeWise9_EXPL_MATRIX(mat,i-0,j-1,UTR5) + GNE_UTR_5SS(mat->evi,i,mat->gen,j);   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for UTR5_INTRON */ 
        /* Add any movement independant score and put away */ 
         GenomeWise9_EXPL_MATRIX(mat,i,j,UTR5_INTRON) = score;   


        /* Finished calculating state UTR5_INTRON */ 


        /* For state START_CODON */ 
        /* setting first movement to score */ 
        score = GenomeWise9_EXPL_MATRIX(mat,i-0,j-3,UTR5_INTRON) + GNE_START_CODON(mat->evi,i,mat->gen,j);   
        /* From state UTR5 to state START_CODON */ 
        temp = GenomeWise9_EXPL_MATRIX(mat,i-0,j-3,UTR5) + GNE_START_CODON(mat->evi,i,mat->gen,j);   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for START_CODON */ 
        /* Add any movement independant score and put away */ 
         GenomeWise9_EXPL_MATRIX(mat,i,j,START_CODON) = score;   


        /* Finished calculating state START_CODON */ 


        /* For state CDS */ 
        /* setting first movement to score */ 
        score = GenomeWise9_EXPL_MATRIX(mat,i-0,j-3,CDS) + GNE_CDS(mat->evi,i,mat->gen,j);   
        /* From state CDS_INTRON_0 to state CDS */ 
        temp = GenomeWise9_EXPL_MATRIX(mat,i-0,j-6,CDS_INTRON_0) + GNE_CDS_3SS(mat->evi,i,mat->gen,(j-3),0);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state CDS_INTRON_1 to state CDS */ 
        temp = GenomeWise9_EXPL_MATRIX(mat,i-0,j-5,CDS_INTRON_1) + GNE_CDS_3SS(mat->evi,i,mat->gen,(j-2),1);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state CDS_INTRON_2 to state CDS */ 
        temp = GenomeWise9_EXPL_MATRIX(mat,i-0,j-4,CDS_INTRON_2) + GNE_CDS_3SS(mat->evi,i,mat->gen,(j-1),2);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state CDS to state CDS */ 
        temp = GenomeWise9_EXPL_MATRIX(mat,i-0,j-2,CDS) + GNE_CDS_FRAMESHIFT(mat->evi,i,mat->gen,j,2);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state CDS to state CDS */ 
        temp = GenomeWise9_EXPL_MATRIX(mat,i-0,j-4,CDS) + GNE_CDS_FRAMESHIFT(mat->evi,i,mat->gen,j,4);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state UTR5 to state CDS */ 
        temp = GenomeWise9_EXPL_MATRIX(mat,i-0,j-3,UTR5) + (GNE_CDS(mat->evi,i,mat->gen,j)+mat->non_start_codon);    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state START_CODON to state CDS */ 
        temp = GenomeWise9_EXPL_MATRIX(mat,i-0,j-3,START_CODON) + GNE_CDS(mat->evi,i,mat->gen,j);    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SPECIAL_CDS to state CDS */ 
        temp = GenomeWise9_EXPL_SPECIAL(mat,i-0,j-3,SPECIAL_CDS) + GNE_CDS(mat->evi,i,mat->gen,j);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTERGENIC to state CDS */ 
        temp = GenomeWise9_EXPL_SPECIAL(mat,i-0,j-3,INTERGENIC) + (GNE_CDS(mat->evi,i,mat->gen,j)+GNE_UTR5_START(mat->evi,i,mat->gen,j));    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for CDS */ 
        /* Add any movement independant score and put away */ 
         GenomeWise9_EXPL_MATRIX(mat,i,j,CDS) = score;   


        /* state CDS is a source for special INTERGENIC */ 
        temp = score + ((mat->newgenecost+GNE_UTR3_END(mat->evi,i,mat->gen,j))) + (0) ;  
        if( temp > GenomeWise9_EXPL_SPECIAL(mat,i,j,INTERGENIC) )    {  
          GenomeWise9_EXPL_SPECIAL(mat,i,j,INTERGENIC) = temp;   
          }  




        /* state CDS is a source for special SPECIAL_CDS */ 
        temp = score + ((mat->switchcost+mat->rndcodon->codon[CSEQ_GENOMIC_CODON(mat->gen,j)])) + (0) ;  
        if( temp > GenomeWise9_EXPL_SPECIAL(mat,i,j,SPECIAL_CDS) )   {  
          GenomeWise9_EXPL_SPECIAL(mat,i,j,SPECIAL_CDS) = temp;  
          }  




        /* Finished calculating state CDS */ 


        /* For state CDS_INTRON_0 */ 
        /* setting first movement to score */ 
        score = GenomeWise9_EXPL_MATRIX(mat,i-0,j-1,CDS_INTRON_0) + GNE_CDS_INTRON(mat->evi,i,mat->gen,j);   
        /* From state CDS to state CDS_INTRON_0 */ 
        temp = GenomeWise9_EXPL_MATRIX(mat,i-0,j-8,CDS) + GNE_CDS_5SS(mat->evi,i,mat->gen,(j-7),0);  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for CDS_INTRON_0 */ 
        /* Add any movement independant score and put away */ 
         GenomeWise9_EXPL_MATRIX(mat,i,j,CDS_INTRON_0) = score;  


        /* Finished calculating state CDS_INTRON_0 */ 


        /* For state CDS_INTRON_1 */ 
        /* setting first movement to score */ 
        score = GenomeWise9_EXPL_MATRIX(mat,i-0,j-1,CDS_INTRON_1) + GNE_CDS_INTRON(mat->evi,i,mat->gen,j);   
        /* From state CDS to state CDS_INTRON_1 */ 
        temp = GenomeWise9_EXPL_MATRIX(mat,i-0,j-9,CDS) + GNE_CDS_5SS(mat->evi,i,mat->gen,(j-7),1);  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for CDS_INTRON_1 */ 
        /* Add any movement independant score and put away */ 
         GenomeWise9_EXPL_MATRIX(mat,i,j,CDS_INTRON_1) = score;  


        /* Finished calculating state CDS_INTRON_1 */ 


        /* For state CDS_INTRON_2 */ 
        /* setting first movement to score */ 
        score = GenomeWise9_EXPL_MATRIX(mat,i-0,j-1,CDS_INTRON_2) + GNE_CDS_INTRON(mat->evi,i,mat->gen,j);   
        /* From state CDS to state CDS_INTRON_2 */ 
        temp = GenomeWise9_EXPL_MATRIX(mat,i-0,j-10,CDS) + GNE_CDS_5SS(mat->evi,i,mat->gen,(j-7),2);     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for CDS_INTRON_2 */ 
        /* Add any movement independant score and put away */ 
         GenomeWise9_EXPL_MATRIX(mat,i,j,CDS_INTRON_2) = score;  


        /* Finished calculating state CDS_INTRON_2 */ 


        /* For state STOP_CODON */ 
        /* setting first movement to score */ 
        score = GenomeWise9_EXPL_MATRIX(mat,i-0,j-3,CDS) + GNE_STOP_CODON(mat->evi,i,mat->gen,j);    


        /* Ok - finished max calculation for STOP_CODON */ 
        /* Add any movement independant score and put away */ 
         GenomeWise9_EXPL_MATRIX(mat,i,j,STOP_CODON) = score;    


        /* Finished calculating state STOP_CODON */ 


        /* For state UTR3 */ 
        /* setting first movement to score */ 
        score = GenomeWise9_EXPL_MATRIX(mat,i-0,j-1,UTR3) + GNE_UTR(mat->evi,i,mat->gen,j);  
        /* From state CDS to state UTR3 */ 
        temp = GenomeWise9_EXPL_MATRIX(mat,i-0,j-1,CDS) + (GNE_UTR(mat->evi,i,mat->gen,j)+mat->non_stop_codon);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state STOP_CODON to state UTR3 */ 
        temp = GenomeWise9_EXPL_MATRIX(mat,i-0,j-1,STOP_CODON) + GNE_UTR(mat->evi,i,mat->gen,j);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state UTR3_INTRON to state UTR3 */ 
        temp = GenomeWise9_EXPL_MATRIX(mat,i-0,j-1,UTR3_INTRON) + GNE_UTR_3SS(mat->evi,i,mat->gen,j);    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* Has restricted position */ 
        if( (j-1) == 0  )    {  
          /* From state INTERGENIC to state UTR3 */ 
          temp = GenomeWise9_EXPL_SPECIAL(mat,i-0,j-1,INTERGENIC) + GNE_UTR(mat->evi,i,mat->gen,j);  
          if( temp  > score )    {  
            score = temp;    
            }  
          }  


        /* Ok - finished max calculation for UTR3 */ 
        /* Add any movement independant score and put away */ 
         GenomeWise9_EXPL_MATRIX(mat,i,j,UTR3) = score;  


        /* state UTR3 is a source for special POSTGENE_INTERGENIC */ 
        temp = score + ((mat->newgenecost+GNE_UTR3_END(mat->evi,i,mat->gen,j))) + (0) ;  
        if( temp > GenomeWise9_EXPL_SPECIAL(mat,i,j,POSTGENE_INTERGENIC) )   {  
          GenomeWise9_EXPL_SPECIAL(mat,i,j,POSTGENE_INTERGENIC) = temp;  
          }  




        /* state UTR3 is a source for special INTERGENIC */ 
        temp = score + ((mat->newgenecost+GNE_UTR3_END(mat->evi,i,mat->gen,j))) + (0) ;  
        if( temp > GenomeWise9_EXPL_SPECIAL(mat,i,j,INTERGENIC) )    {  
          GenomeWise9_EXPL_SPECIAL(mat,i,j,INTERGENIC) = temp;   
          }  




        /* state UTR3 is a source for special SPECIAL_UTR3 */ 
        temp = score + (mat->switchcost) + (0) ;     
        if( temp > GenomeWise9_EXPL_SPECIAL(mat,i,j,SPECIAL_UTR3) )  {  
          GenomeWise9_EXPL_SPECIAL(mat,i,j,SPECIAL_UTR3) = temp;     
          }  




        /* Finished calculating state UTR3 */ 


        /* For state UTR3_INTRON */ 
        /* setting first movement to score */ 
        score = GenomeWise9_EXPL_MATRIX(mat,i-0,j-1,UTR3_INTRON) + GNE_UTR_INTRON(mat->evi,i,mat->gen,j);    
        /* From state UTR3 to state UTR3_INTRON */ 
        temp = GenomeWise9_EXPL_MATRIX(mat,i-0,j-1,UTR3) + GNE_UTR_5SS(mat->evi,i,mat->gen,j);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SPECIAL_UTR3 to state UTR3_INTRON */ 
        temp = GenomeWise9_EXPL_SPECIAL(mat,i-0,j-1,SPECIAL_UTR3) + GNE_UTR(mat->evi,i,mat->gen,j);  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for UTR3_INTRON */ 
        /* Add any movement independant score and put away */ 
         GenomeWise9_EXPL_MATRIX(mat,i,j,UTR3_INTRON) = score;   


        /* Finished calculating state UTR3_INTRON */ 
        }  


      /* Special state PREGENE_INTERGENIC has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GenomeWise9_EXPL_SPECIAL(mat,0,j,PREGENE_INTERGENIC);  


      /* Source START is a special source for PREGENE_INTERGENIC */ 
      temp = GenomeWise9_EXPL_SPECIAL(mat,0,j - 1,START) + (0) + (0);    
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      GenomeWise9_EXPL_SPECIAL(mat,0,j,PREGENE_INTERGENIC) = score;  
      /* Finished updating state PREGENE_INTERGENIC */ 




      /* Special state POSTGENE_INTERGENIC has no special to special movements */ 


      /* Special state INTERGENIC has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GenomeWise9_EXPL_SPECIAL(mat,0,j,INTERGENIC);  


      /* Source INTERGENIC is a special source for INTERGENIC */ 
      temp = GenomeWise9_EXPL_SPECIAL(mat,0,j - 1,INTERGENIC) + (0) + (0);   
      if( temp > score ) 
        score = temp;    


      /* Source CDS for state INTERGENIC is not special... already calculated */ 
      /* Source UTR3 for state INTERGENIC is not special... already calculated */ 
      /* Put back score... (now updated!) */ 
      GenomeWise9_EXPL_SPECIAL(mat,0,j,INTERGENIC) = score;  
      /* Finished updating state INTERGENIC */ 




      /* Special state SPECIAL_UTR5 has no special to special movements */ 


      /* Special state SPECIAL_UTR3 has no special to special movements */ 


      /* Special state SPECIAL_CDS has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GenomeWise9_EXPL_SPECIAL(mat,0,j,SPECIAL_CDS); 


      /* Source CDS for state SPECIAL_CDS is not special... already calculated */ 
      /* Source SPECIAL_CDS is a special source for SPECIAL_CDS */ 
      temp = GenomeWise9_EXPL_SPECIAL(mat,0,j - 3,SPECIAL_CDS) + (mat->rndcodon->codon[CSEQ_GENOMIC_CODON(mat->gen,j)]) + (0);   
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      GenomeWise9_EXPL_SPECIAL(mat,0,j,SPECIAL_CDS) = score; 
      /* Finished updating state SPECIAL_CDS */ 




      /* Special state START has no special to special movements */ 


      /* Special state END has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GenomeWise9_EXPL_SPECIAL(mat,0,j,END); 


      /* Source INTERGENIC is a special source for END */ 
      /* Has restricted position */ 
      if( j == mat->lenj-1 ) {  
        temp = GenomeWise9_EXPL_SPECIAL(mat,0,j - 1,INTERGENIC) + (0) + (0);     
        if( temp > score )   
          score = temp;  
        }  


      /* Put back score... (now updated!) */ 
      GenomeWise9_EXPL_SPECIAL(mat,0,j,END) = score; 
      /* Finished updating state END */ 


      }  
    stop_reporting();    
    return TRUE;     
}    


/* Function:  calculate_dpenv_GenomeWise9(mat,dpenv)
 *
 * Descrip:    This function calculates the GenomeWise9 matrix when in explicit mode, subject to the envelope
 *
 *
 * Arg:          mat [UNKN ] GenomeWise9 which contains explicit basematrix memory [GenomeWise9 *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean calculate_dpenv_GenomeWise9(GenomeWise9 * mat,DPEnvelope * dpenv) 
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
      warn("in calculate_GenomeWise9, passed a non Explicit matrix type, cannot calculate!");    
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
        GenomeWise9_EXPL_MATRIX(mat,i,j,UTR5) = NEGI;    
        GenomeWise9_EXPL_MATRIX(mat,i,j,UTR5_INTRON) = NEGI; 
        GenomeWise9_EXPL_MATRIX(mat,i,j,START_CODON) = NEGI; 
        GenomeWise9_EXPL_MATRIX(mat,i,j,CDS) = NEGI; 
        GenomeWise9_EXPL_MATRIX(mat,i,j,CDS_INTRON_0) = NEGI;    
        GenomeWise9_EXPL_MATRIX(mat,i,j,CDS_INTRON_1) = NEGI;    
        GenomeWise9_EXPL_MATRIX(mat,i,j,CDS_INTRON_2) = NEGI;    
        GenomeWise9_EXPL_MATRIX(mat,i,j,STOP_CODON) = NEGI;  
        GenomeWise9_EXPL_MATRIX(mat,i,j,UTR3) = NEGI;    
        GenomeWise9_EXPL_MATRIX(mat,i,j,UTR3_INTRON) = NEGI; 
        }  
      }  
    for(j=-10;j<mat->lenj;j++)   {  
      GenomeWise9_EXPL_SPECIAL(mat,i,j,PREGENE_INTERGENIC) = NEGI;   
      GenomeWise9_EXPL_SPECIAL(mat,i,j,POSTGENE_INTERGENIC) = NEGI;  
      GenomeWise9_EXPL_SPECIAL(mat,i,j,INTERGENIC) = NEGI;   
      GenomeWise9_EXPL_SPECIAL(mat,i,j,SPECIAL_UTR5) = NEGI; 
      GenomeWise9_EXPL_SPECIAL(mat,i,j,SPECIAL_UTR3) = NEGI; 
      GenomeWise9_EXPL_SPECIAL(mat,i,j,SPECIAL_CDS) = NEGI;  
      GenomeWise9_EXPL_SPECIAL(mat,i,j,START) = 0;   
      GenomeWise9_EXPL_SPECIAL(mat,i,j,END) = NEGI;  
      }  


    start_reporting("GenomeWise9 Matrix calculation: "); 
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
          GenomeWise9_EXPL_MATRIX(mat,i,j,UTR5) = NEGI;  
          GenomeWise9_EXPL_MATRIX(mat,i,j,UTR5_INTRON) = NEGI;   
          GenomeWise9_EXPL_MATRIX(mat,i,j,START_CODON) = NEGI;   
          GenomeWise9_EXPL_MATRIX(mat,i,j,CDS) = NEGI;   
          GenomeWise9_EXPL_MATRIX(mat,i,j,CDS_INTRON_0) = NEGI;  
          GenomeWise9_EXPL_MATRIX(mat,i,j,CDS_INTRON_1) = NEGI;  
          GenomeWise9_EXPL_MATRIX(mat,i,j,CDS_INTRON_2) = NEGI;  
          GenomeWise9_EXPL_MATRIX(mat,i,j,STOP_CODON) = NEGI;    
          GenomeWise9_EXPL_MATRIX(mat,i,j,UTR3) = NEGI;  
          GenomeWise9_EXPL_MATRIX(mat,i,j,UTR3_INTRON) = NEGI;   
          continue;  
          }  


        if( num%1000 == 0 )  
          log_full_error(REPORT,0,"[%7d] Cells %2d%%%%",num,num*100/tot);    
        num++;   


        /* For state UTR5 */ 
        /* setting first movement to score */ 
        score = GenomeWise9_EXPL_MATRIX(mat,i-0,j-1,UTR5) + GNE_UTR(mat->evi,i,mat->gen,j);  
        /* From state UTR5_INTRON to state UTR5 */ 
        temp = GenomeWise9_EXPL_MATRIX(mat,i-0,j-1,UTR5_INTRON) + GNE_UTR_3SS(mat->evi,i,mat->gen,j);    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SPECIAL_UTR5 to state UTR5 */ 
        temp = GenomeWise9_EXPL_SPECIAL(mat,i-0,j-1,SPECIAL_UTR5) + GNE_UTR(mat->evi,i,mat->gen,j);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state START to state UTR5 */ 
        temp = GenomeWise9_EXPL_SPECIAL(mat,i-0,j-1,START) + GNE_UTR5_START(mat->evi,i,mat->gen,j);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTERGENIC to state UTR5 */ 
        temp = GenomeWise9_EXPL_SPECIAL(mat,i-0,j-1,INTERGENIC) + (GNE_UTR(mat->evi,i,mat->gen,j)+GNE_UTR5_START(mat->evi,i,mat->gen,j));    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state PREGENE_INTERGENIC to state UTR5 */ 
        temp = GenomeWise9_EXPL_SPECIAL(mat,i-0,j-1,PREGENE_INTERGENIC) + GNE_UTR(mat->evi,i,mat->gen,j);    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for UTR5 */ 
        /* Add any movement independant score and put away */ 
         GenomeWise9_EXPL_MATRIX(mat,i,j,UTR5) = score;  


        /* state UTR5 is a source for special SPECIAL_UTR5 */ 
        temp = score + (mat->switchcost) + (0) ;     
        if( temp > GenomeWise9_EXPL_SPECIAL(mat,i,j,SPECIAL_UTR5) )  {  
          GenomeWise9_EXPL_SPECIAL(mat,i,j,SPECIAL_UTR5) = temp;     
          }  




        /* Finished calculating state UTR5 */ 


        /* For state UTR5_INTRON */ 
        /* setting first movement to score */ 
        score = GenomeWise9_EXPL_MATRIX(mat,i-0,j-1,UTR5_INTRON) + GNE_UTR_INTRON(mat->evi,i,mat->gen,j);    
        /* From state UTR5 to state UTR5_INTRON */ 
        temp = GenomeWise9_EXPL_MATRIX(mat,i-0,j-1,UTR5) + GNE_UTR_5SS(mat->evi,i,mat->gen,j);   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for UTR5_INTRON */ 
        /* Add any movement independant score and put away */ 
         GenomeWise9_EXPL_MATRIX(mat,i,j,UTR5_INTRON) = score;   


        /* Finished calculating state UTR5_INTRON */ 


        /* For state START_CODON */ 
        /* setting first movement to score */ 
        score = GenomeWise9_EXPL_MATRIX(mat,i-0,j-3,UTR5_INTRON) + GNE_START_CODON(mat->evi,i,mat->gen,j);   
        /* From state UTR5 to state START_CODON */ 
        temp = GenomeWise9_EXPL_MATRIX(mat,i-0,j-3,UTR5) + GNE_START_CODON(mat->evi,i,mat->gen,j);   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for START_CODON */ 
        /* Add any movement independant score and put away */ 
         GenomeWise9_EXPL_MATRIX(mat,i,j,START_CODON) = score;   


        /* Finished calculating state START_CODON */ 


        /* For state CDS */ 
        /* setting first movement to score */ 
        score = GenomeWise9_EXPL_MATRIX(mat,i-0,j-3,CDS) + GNE_CDS(mat->evi,i,mat->gen,j);   
        /* From state CDS_INTRON_0 to state CDS */ 
        temp = GenomeWise9_EXPL_MATRIX(mat,i-0,j-6,CDS_INTRON_0) + GNE_CDS_3SS(mat->evi,i,mat->gen,(j-3),0);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state CDS_INTRON_1 to state CDS */ 
        temp = GenomeWise9_EXPL_MATRIX(mat,i-0,j-5,CDS_INTRON_1) + GNE_CDS_3SS(mat->evi,i,mat->gen,(j-2),1);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state CDS_INTRON_2 to state CDS */ 
        temp = GenomeWise9_EXPL_MATRIX(mat,i-0,j-4,CDS_INTRON_2) + GNE_CDS_3SS(mat->evi,i,mat->gen,(j-1),2);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state CDS to state CDS */ 
        temp = GenomeWise9_EXPL_MATRIX(mat,i-0,j-2,CDS) + GNE_CDS_FRAMESHIFT(mat->evi,i,mat->gen,j,2);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state CDS to state CDS */ 
        temp = GenomeWise9_EXPL_MATRIX(mat,i-0,j-4,CDS) + GNE_CDS_FRAMESHIFT(mat->evi,i,mat->gen,j,4);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state UTR5 to state CDS */ 
        temp = GenomeWise9_EXPL_MATRIX(mat,i-0,j-3,UTR5) + (GNE_CDS(mat->evi,i,mat->gen,j)+mat->non_start_codon);    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state START_CODON to state CDS */ 
        temp = GenomeWise9_EXPL_MATRIX(mat,i-0,j-3,START_CODON) + GNE_CDS(mat->evi,i,mat->gen,j);    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SPECIAL_CDS to state CDS */ 
        temp = GenomeWise9_EXPL_SPECIAL(mat,i-0,j-3,SPECIAL_CDS) + GNE_CDS(mat->evi,i,mat->gen,j);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTERGENIC to state CDS */ 
        temp = GenomeWise9_EXPL_SPECIAL(mat,i-0,j-3,INTERGENIC) + (GNE_CDS(mat->evi,i,mat->gen,j)+GNE_UTR5_START(mat->evi,i,mat->gen,j));    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for CDS */ 
        /* Add any movement independant score and put away */ 
         GenomeWise9_EXPL_MATRIX(mat,i,j,CDS) = score;   


        /* state CDS is a source for special INTERGENIC */ 
        temp = score + ((mat->newgenecost+GNE_UTR3_END(mat->evi,i,mat->gen,j))) + (0) ;  
        if( temp > GenomeWise9_EXPL_SPECIAL(mat,i,j,INTERGENIC) )    {  
          GenomeWise9_EXPL_SPECIAL(mat,i,j,INTERGENIC) = temp;   
          }  




        /* state CDS is a source for special SPECIAL_CDS */ 
        temp = score + ((mat->switchcost+mat->rndcodon->codon[CSEQ_GENOMIC_CODON(mat->gen,j)])) + (0) ;  
        if( temp > GenomeWise9_EXPL_SPECIAL(mat,i,j,SPECIAL_CDS) )   {  
          GenomeWise9_EXPL_SPECIAL(mat,i,j,SPECIAL_CDS) = temp;  
          }  




        /* Finished calculating state CDS */ 


        /* For state CDS_INTRON_0 */ 
        /* setting first movement to score */ 
        score = GenomeWise9_EXPL_MATRIX(mat,i-0,j-1,CDS_INTRON_0) + GNE_CDS_INTRON(mat->evi,i,mat->gen,j);   
        /* From state CDS to state CDS_INTRON_0 */ 
        temp = GenomeWise9_EXPL_MATRIX(mat,i-0,j-8,CDS) + GNE_CDS_5SS(mat->evi,i,mat->gen,(j-7),0);  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for CDS_INTRON_0 */ 
        /* Add any movement independant score and put away */ 
         GenomeWise9_EXPL_MATRIX(mat,i,j,CDS_INTRON_0) = score;  


        /* Finished calculating state CDS_INTRON_0 */ 


        /* For state CDS_INTRON_1 */ 
        /* setting first movement to score */ 
        score = GenomeWise9_EXPL_MATRIX(mat,i-0,j-1,CDS_INTRON_1) + GNE_CDS_INTRON(mat->evi,i,mat->gen,j);   
        /* From state CDS to state CDS_INTRON_1 */ 
        temp = GenomeWise9_EXPL_MATRIX(mat,i-0,j-9,CDS) + GNE_CDS_5SS(mat->evi,i,mat->gen,(j-7),1);  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for CDS_INTRON_1 */ 
        /* Add any movement independant score and put away */ 
         GenomeWise9_EXPL_MATRIX(mat,i,j,CDS_INTRON_1) = score;  


        /* Finished calculating state CDS_INTRON_1 */ 


        /* For state CDS_INTRON_2 */ 
        /* setting first movement to score */ 
        score = GenomeWise9_EXPL_MATRIX(mat,i-0,j-1,CDS_INTRON_2) + GNE_CDS_INTRON(mat->evi,i,mat->gen,j);   
        /* From state CDS to state CDS_INTRON_2 */ 
        temp = GenomeWise9_EXPL_MATRIX(mat,i-0,j-10,CDS) + GNE_CDS_5SS(mat->evi,i,mat->gen,(j-7),2);     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for CDS_INTRON_2 */ 
        /* Add any movement independant score and put away */ 
         GenomeWise9_EXPL_MATRIX(mat,i,j,CDS_INTRON_2) = score;  


        /* Finished calculating state CDS_INTRON_2 */ 


        /* For state STOP_CODON */ 
        /* setting first movement to score */ 
        score = GenomeWise9_EXPL_MATRIX(mat,i-0,j-3,CDS) + GNE_STOP_CODON(mat->evi,i,mat->gen,j);    


        /* Ok - finished max calculation for STOP_CODON */ 
        /* Add any movement independant score and put away */ 
         GenomeWise9_EXPL_MATRIX(mat,i,j,STOP_CODON) = score;    


        /* Finished calculating state STOP_CODON */ 


        /* For state UTR3 */ 
        /* setting first movement to score */ 
        score = GenomeWise9_EXPL_MATRIX(mat,i-0,j-1,UTR3) + GNE_UTR(mat->evi,i,mat->gen,j);  
        /* From state CDS to state UTR3 */ 
        temp = GenomeWise9_EXPL_MATRIX(mat,i-0,j-1,CDS) + (GNE_UTR(mat->evi,i,mat->gen,j)+mat->non_stop_codon);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state STOP_CODON to state UTR3 */ 
        temp = GenomeWise9_EXPL_MATRIX(mat,i-0,j-1,STOP_CODON) + GNE_UTR(mat->evi,i,mat->gen,j);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state UTR3_INTRON to state UTR3 */ 
        temp = GenomeWise9_EXPL_MATRIX(mat,i-0,j-1,UTR3_INTRON) + GNE_UTR_3SS(mat->evi,i,mat->gen,j);    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* Has restricted position */ 
        if( (j-1) == 0  )    {  
          /* From state INTERGENIC to state UTR3 */ 
          temp = GenomeWise9_EXPL_SPECIAL(mat,i-0,j-1,INTERGENIC) + GNE_UTR(mat->evi,i,mat->gen,j);  
          if( temp  > score )    {  
            score = temp;    
            }  
          }  


        /* Ok - finished max calculation for UTR3 */ 
        /* Add any movement independant score and put away */ 
         GenomeWise9_EXPL_MATRIX(mat,i,j,UTR3) = score;  


        /* state UTR3 is a source for special POSTGENE_INTERGENIC */ 
        temp = score + ((mat->newgenecost+GNE_UTR3_END(mat->evi,i,mat->gen,j))) + (0) ;  
        if( temp > GenomeWise9_EXPL_SPECIAL(mat,i,j,POSTGENE_INTERGENIC) )   {  
          GenomeWise9_EXPL_SPECIAL(mat,i,j,POSTGENE_INTERGENIC) = temp;  
          }  




        /* state UTR3 is a source for special INTERGENIC */ 
        temp = score + ((mat->newgenecost+GNE_UTR3_END(mat->evi,i,mat->gen,j))) + (0) ;  
        if( temp > GenomeWise9_EXPL_SPECIAL(mat,i,j,INTERGENIC) )    {  
          GenomeWise9_EXPL_SPECIAL(mat,i,j,INTERGENIC) = temp;   
          }  




        /* state UTR3 is a source for special SPECIAL_UTR3 */ 
        temp = score + (mat->switchcost) + (0) ;     
        if( temp > GenomeWise9_EXPL_SPECIAL(mat,i,j,SPECIAL_UTR3) )  {  
          GenomeWise9_EXPL_SPECIAL(mat,i,j,SPECIAL_UTR3) = temp;     
          }  




        /* Finished calculating state UTR3 */ 


        /* For state UTR3_INTRON */ 
        /* setting first movement to score */ 
        score = GenomeWise9_EXPL_MATRIX(mat,i-0,j-1,UTR3_INTRON) + GNE_UTR_INTRON(mat->evi,i,mat->gen,j);    
        /* From state UTR3 to state UTR3_INTRON */ 
        temp = GenomeWise9_EXPL_MATRIX(mat,i-0,j-1,UTR3) + GNE_UTR_5SS(mat->evi,i,mat->gen,j);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SPECIAL_UTR3 to state UTR3_INTRON */ 
        temp = GenomeWise9_EXPL_SPECIAL(mat,i-0,j-1,SPECIAL_UTR3) + GNE_UTR(mat->evi,i,mat->gen,j);  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for UTR3_INTRON */ 
        /* Add any movement independant score and put away */ 
         GenomeWise9_EXPL_MATRIX(mat,i,j,UTR3_INTRON) = score;   


        /* Finished calculating state UTR3_INTRON */ 
        }  


      /* Special state PREGENE_INTERGENIC has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GenomeWise9_EXPL_SPECIAL(mat,0,j,PREGENE_INTERGENIC);  


      /* Source START is a special source for PREGENE_INTERGENIC */ 
      temp = GenomeWise9_EXPL_SPECIAL(mat,0,j - 1,START) + (0) + (0);    
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      GenomeWise9_EXPL_SPECIAL(mat,0,j,PREGENE_INTERGENIC) = score;  
      /* Finished updating state PREGENE_INTERGENIC */ 




      /* Special state POSTGENE_INTERGENIC has no special to special movements */ 


      /* Special state INTERGENIC has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GenomeWise9_EXPL_SPECIAL(mat,0,j,INTERGENIC);  


      /* Source INTERGENIC is a special source for INTERGENIC */ 
      temp = GenomeWise9_EXPL_SPECIAL(mat,0,j - 1,INTERGENIC) + (0) + (0);   
      if( temp > score ) 
        score = temp;    


      /* Source CDS for state INTERGENIC is not special... already calculated */ 
      /* Source UTR3 for state INTERGENIC is not special... already calculated */ 
      /* Put back score... (now updated!) */ 
      GenomeWise9_EXPL_SPECIAL(mat,0,j,INTERGENIC) = score;  
      /* Finished updating state INTERGENIC */ 




      /* Special state SPECIAL_UTR5 has no special to special movements */ 


      /* Special state SPECIAL_UTR3 has no special to special movements */ 


      /* Special state SPECIAL_CDS has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GenomeWise9_EXPL_SPECIAL(mat,0,j,SPECIAL_CDS); 


      /* Source CDS for state SPECIAL_CDS is not special... already calculated */ 
      /* Source SPECIAL_CDS is a special source for SPECIAL_CDS */ 
      temp = GenomeWise9_EXPL_SPECIAL(mat,0,j - 3,SPECIAL_CDS) + (mat->rndcodon->codon[CSEQ_GENOMIC_CODON(mat->gen,j)]) + (0);   
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      GenomeWise9_EXPL_SPECIAL(mat,0,j,SPECIAL_CDS) = score; 
      /* Finished updating state SPECIAL_CDS */ 




      /* Special state START has no special to special movements */ 


      /* Special state END has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GenomeWise9_EXPL_SPECIAL(mat,0,j,END); 


      /* Source INTERGENIC is a special source for END */ 
      /* Has restricted position */ 
      if( j == mat->lenj-1 ) {  
        temp = GenomeWise9_EXPL_SPECIAL(mat,0,j - 1,INTERGENIC) + (0) + (0);     
        if( temp > score )   
          score = temp;  
        }  


      /* Put back score... (now updated!) */ 
      GenomeWise9_EXPL_SPECIAL(mat,0,j,END) = score; 
      /* Finished updating state END */ 


      }  
    stop_reporting();    
    return TRUE;     
}    


/* Function:  GenomeWise9_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GenomeWise9 *]
 *
 */
GenomeWise9 * GenomeWise9_alloc(void) 
{
    GenomeWise9 * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(GenomeWise9 *) ckalloc (sizeof(GenomeWise9))) == NULL)  {  
      warn("GenomeWise9_alloc failed "); 
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


/* Function:  free_GenomeWise9(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GenomeWise9 *]
 *
 * Return [UNKN ]  Undocumented return value [GenomeWise9 *]
 *
 */
GenomeWise9 * free_GenomeWise9(GenomeWise9 * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a GenomeWise9 obj. Should be trappable");   
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
    /* obj->evi is linked in */ 
    /* obj->gen is linked in */ 
    /* obj->switchcost is linked in */ 
    /* obj->newgenecost is linked in */ 
    /* obj->non_start_codon is linked in */ 
    /* obj->non_stop_codon is linked in */ 
    /* obj->rndcodon is linked in */ 


    ckfree(obj); 
    return NULL; 
}    





#ifdef _cplusplus
}
#endif
