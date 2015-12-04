#ifdef _cplusplus
extern "C" {
#endif
#include "genewise21.h"



# line 7 "genewise21.c"


  /*****************   C functions  ****************/
  /*             Written using dynamite            */
  /*            Sat Sep  8 09:05:32 2007           */
  /*            email birney@sanger.ac.uk          */
  /* http://www.sanger.ac.uk/Users/birney/dynamite */
  /*************************************************/


  /* Please report any problems or bugs to         */
  /* Ewan Birney, birney@sanger.ac.uk              */


/* basic set of macros to map states to numbers */ 
#define MATCH 0  
#define INSERT 1 
#define DELETE 2 
#define CENTRAL_0M 3 
#define PY_0M 4  
#define SPACER_0M 5  
#define CENTRAL_1M 6 
#define PY_1M 7  
#define SPACER_1M 8  
#define CENTRAL_2M 9 
#define PY_2M 10 
#define SPACER_2M 11 
#define CENTRAL_0I 12    
#define PY_0I 13 
#define SPACER_0I 14 
#define CENTRAL_1I 15    
#define PY_1I 16 
#define SPACER_1I 17 
#define CENTRAL_2I 18    
#define PY_2I 19 
#define SPACER_2I 20 


#define START 0  
#define END 1    
#define BEFORE_RND 2 
#define BEFORE_CODON 3   
#define BEFORE_CENTRAL 4 
#define BEFORE_PY_TRACT 5    
#define BEFORE_SPACER 6  
#define AFTER_RND 7  
#define AFTER_CODON 8    
#define AFTER_CENTRAL 9  
#define AFTER_PY_TRACT 10    
#define AFTER_SPACER 11  


#define GeneWise21_EXPL_MATRIX(this_matrix,i,j,STATE) this_matrix->basematrix->matrix[((j+10)*21)+STATE][i+1]    
#define GeneWise21_EXPL_SPECIAL(matrix,i,j,STATE) matrix->basematrix->specmatrix[STATE][j+10]    
#define GeneWise21_READ_OFF_ERROR -12
   


#define GeneWise21_VSMALL_MATRIX(mat,i,j,STATE) mat->basematrix->matrix[(j+11)%11][((i+1)*21)+STATE] 
#define GeneWise21_VSMALL_SPECIAL(mat,i,j,STATE) mat->basematrix->specmatrix[(j+11)%11][STATE]   




#define GeneWise21_SHATTER_SPECIAL(matrix,i,j,STATE) matrix->shatter->special[STATE][j]  
#define GeneWise21_SHATTER_MATRIX(matrix,i,j,STATE)  fetch_cell_value_ShatterMatrix(mat->shatter,i,j,STATE)  


/* Function:  PackAln_read_Shatter_GeneWise21(mat)
 *
 * Descrip:    Reads off PackAln from shatter matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [GeneWise21 *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * PackAln_read_Shatter_GeneWise21(GeneWise21 * mat) 
{
    GeneWise21_access_func_holder holder;    


    holder.access_main    = GeneWise21_shatter_access_main;  
    holder.access_special = GeneWise21_shatter_access_special;   
    assert(mat);     
    assert(mat->shatter);    
    return PackAln_read_generic_GeneWise21(mat,holder);  
}    


/* Function:  GeneWise21_shatter_access_main(mat,i,j,state)
 *
 * Descrip: No Description
 *
 * Arg:          mat [UNKN ] Undocumented argument [GeneWise21 *]
 * Arg:            i [UNKN ] Undocumented argument [int]
 * Arg:            j [UNKN ] Undocumented argument [int]
 * Arg:        state [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int GeneWise21_shatter_access_main(GeneWise21 * mat,int i,int j,int state) 
{
    return GeneWise21_SHATTER_MATRIX(mat,i,j,state); 
}    


/* Function:  GeneWise21_shatter_access_special(mat,i,j,state)
 *
 * Descrip: No Description
 *
 * Arg:          mat [UNKN ] Undocumented argument [GeneWise21 *]
 * Arg:            i [UNKN ] Undocumented argument [int]
 * Arg:            j [UNKN ] Undocumented argument [int]
 * Arg:        state [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int GeneWise21_shatter_access_special(GeneWise21 * mat,int i,int j,int state) 
{
    return GeneWise21_SHATTER_SPECIAL(mat,i,j,state);    
}    


/* Function:  calculate_shatter_GeneWise21(mat,dpenv)
 *
 * Descrip:    This function calculates the GeneWise21 matrix when in shatter mode
 *
 *
 * Arg:          mat [UNKN ] (null) [GeneWise21 *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean calculate_shatter_GeneWise21(GeneWise21 * mat,DPEnvelope * dpenv) 
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
    int * SIG_1_6;   
    int * SIG_1_5;   
    int * SIG_1_4;   
    int * SIG_1_2;   
    int * SIG_1_1;   
    int * SIG_0_3;   
    int * SIG_1_0;   
    int * SIG_0_8;   
    int * SIG_0_1;   
    int * SIG_0_9;   
    int * SIG_0_10;  


    leni = mat->leni;    
    lenj = mat->lenj;    


    mat->shatter = new_ShatterMatrix(dpenv,21,lenj,12);  
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


    start_reporting("GeneWise21 Matrix calculation: ");  
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
        SIG_1_6 = fetch_cell_from_ShatterMatrix(mat->shatter,i-1,j-6);   
        SIG_1_5 = fetch_cell_from_ShatterMatrix(mat->shatter,i-1,j-5);   
        SIG_1_4 = fetch_cell_from_ShatterMatrix(mat->shatter,i-1,j-4);   
        SIG_1_2 = fetch_cell_from_ShatterMatrix(mat->shatter,i-1,j-2);   
        SIG_1_1 = fetch_cell_from_ShatterMatrix(mat->shatter,i-1,j-1);   
        SIG_0_3 = fetch_cell_from_ShatterMatrix(mat->shatter,i-0,j-3);   
        SIG_1_0 = fetch_cell_from_ShatterMatrix(mat->shatter,i-1,j-0);   
        SIG_0_8 = fetch_cell_from_ShatterMatrix(mat->shatter,i-0,j-8);   
        SIG_0_1 = fetch_cell_from_ShatterMatrix(mat->shatter,i-0,j-1);   
        SIG_0_9 = fetch_cell_from_ShatterMatrix(mat->shatter,i-0,j-9);   
        SIG_0_10 = fetch_cell_from_ShatterMatrix(mat->shatter,i-0,j-10); 




        /* For state MATCH */ 
        /* setting first movement to score */ 
        score = SIG_1_3[MATCH] + (mat->query->seg[i]->transition[GW_MATCH2MATCH]+mat->query->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]);  
        /* From state INSERT to state MATCH */ 
        temp = SIG_1_3[INSERT] + (mat->query->seg[i]->transition[GW_INSERT2MATCH]+mat->query->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state MATCH */ 
        temp = SIG_1_3[DELETE] + (mat->query->seg[i]->transition[GW_DELETE2MATCH]+mat->query->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state BEFORE_CODON to state MATCH */ 
        temp = GeneWise21_SHATTER_SPECIAL(mat,i-1,j-3,BEFORE_CODON) + (mat->query->seg[i]->transition[GW_START2MATCH]+mat->query->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SPACER_0M to state MATCH */ 
        temp = SIG_1_6[SPACER_0M] + (((mat->gp21->transition[GP21_SPACER2CDS]+mat->query->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)])+CSEQ_GENOMIC_3SS(mat->target,(j-3)))+mat->query->seg[i]->transition[GW_MATCH_BALANCE_3SS]);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state PY_0M to state MATCH */ 
        temp = SIG_1_6[PY_0M] + (((mat->gp21->transition[GP21_PY2CDS]+mat->query->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)])+CSEQ_GENOMIC_3SS(mat->target,(j-3)))+mat->query->seg[i]->transition[GW_MATCH_BALANCE_3SS]);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SPACER_1M to state MATCH */ 
        temp = SIG_1_5[SPACER_1M] + (mat->gp21->transition[GP21_SPACER2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-2)));    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state PY_1M to state MATCH */ 
        temp = SIG_1_5[PY_1M] + (mat->gp21->transition[GP21_PY2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-2)));    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SPACER_2M to state MATCH */ 
        temp = SIG_1_4[SPACER_2M] + (mat->gp21->transition[GP21_SPACER2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-1)));    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state PY_2M to state MATCH */ 
        temp = SIG_1_4[PY_2M] + (mat->gp21->transition[GP21_PY2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-1)));    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = SIG_1_2[MATCH] + mat->gp21->transition[GP21_DELETE_1_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = SIG_1_1[MATCH] + mat->gp21->transition[GP21_DELETE_2_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = SIG_1_4[MATCH] + mat->gp21->transition[GP21_INSERT_1_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = SIG_1_5[MATCH] + mat->gp21->transition[GP21_INSERT_2_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for MATCH */ 
        /* Add any movement independant score and put away */ 
         score += CSEQ_GENOMIC_CDSPOT(mat->target,j);    
         SIG_0_0[MATCH] = score; 


        /* state MATCH is a source for special AFTER_CODON */ 
        temp = score + (mat->query->seg[i]->transition[GW_MATCH2END]) + (0) ;    
        if( temp > GeneWise21_SHATTER_SPECIAL(mat,i,j,AFTER_CODON) )     {  
          GeneWise21_SHATTER_SPECIAL(mat,i,j,AFTER_CODON) = temp;    
          }  




        /* Finished calculating state MATCH */ 


        /* For state INSERT */ 
        /* setting first movement to score */ 
        score = SIG_0_3[MATCH] + (mat->query->seg[i]->transition[GW_MATCH2INSERT]+mat->query->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)]);    
        /* From state INSERT to state INSERT */ 
        temp = SIG_0_3[INSERT] + (mat->query->seg[i]->transition[GW_INSERT2INSERT]+mat->query->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)]);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state INSERT */ 
        temp = SIG_0_3[DELETE] + (mat->query->seg[i]->transition[GW_DELETE2INSERT]+mat->query->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)]);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state BEFORE_CODON to state INSERT */ 
        temp = GeneWise21_SHATTER_SPECIAL(mat,i-0,j-3,BEFORE_CODON) + (mat->query->seg[i]->transition[GW_START2INSERT]+mat->query->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)]);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SPACER_0I to state INSERT */ 
        temp = SIG_1_6[SPACER_0I] + (((mat->gp21->transition[GP21_SPACER2CDS]+mat->query->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)])+CSEQ_GENOMIC_3SS(mat->target,(j-3)))+mat->query->seg[i]->transition[GW_MATCH_BALANCE_3SS]);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SPACER_1I to state INSERT */ 
        temp = SIG_1_5[SPACER_1I] + (mat->gp21->transition[GP21_SPACER2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-2)));    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SPACER_2I to state INSERT */ 
        temp = SIG_1_4[SPACER_2I] + (mat->gp21->transition[GP21_SPACER2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-1)));    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state PY_0I to state INSERT */ 
        temp = SIG_1_6[PY_0I] + (((mat->gp21->transition[GP21_PY2CDS]+mat->query->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)])+CSEQ_GENOMIC_3SS(mat->target,(j-3)))+mat->query->seg[i]->transition[GW_MATCH_BALANCE_3SS]);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state PY_1I to state INSERT */ 
        temp = SIG_1_5[PY_1I] + (mat->gp21->transition[GP21_PY2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-2)));    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state PY_2I to state INSERT */ 
        temp = SIG_1_4[PY_2I] + (mat->gp21->transition[GP21_PY2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-1)));    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state INSERT */ 
        temp = SIG_1_2[INSERT] + mat->gp21->transition[GP21_DELETE_1_BASE];  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state INSERT */ 
        temp = SIG_1_1[INSERT] + mat->gp21->transition[GP21_DELETE_2_BASE];  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state INSERT */ 
        temp = SIG_1_4[INSERT] + mat->gp21->transition[GP21_INSERT_1_BASE];  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state INSERT */ 
        temp = SIG_1_5[INSERT] + mat->gp21->transition[GP21_INSERT_2_BASE];  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for INSERT */ 
        /* Add any movement independant score and put away */ 
         score += CSEQ_GENOMIC_CDSPOT(mat->target,j);    
         SIG_0_0[INSERT] = score;    


        /* state INSERT is a source for special AFTER_CODON */ 
        temp = score + (mat->query->seg[i]->transition[GW_INSERT2END]) + (0) ;   
        if( temp > GeneWise21_SHATTER_SPECIAL(mat,i,j,AFTER_CODON) )     {  
          GeneWise21_SHATTER_SPECIAL(mat,i,j,AFTER_CODON) = temp;    
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
        /* From state BEFORE_CODON to state DELETE */ 
        temp = GeneWise21_SHATTER_SPECIAL(mat,i-1,j-0,BEFORE_CODON) + mat->query->seg[i]->transition[GW_START2DELETE];   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for DELETE */ 
        /* Add any movement independant score and put away */ 
         SIG_0_0[DELETE] = score;    


        /* state DELETE is a source for special AFTER_CODON */ 
        temp = score + (mat->query->seg[i]->transition[GW_DELETE2END]) + (0) ;   
        if( temp > GeneWise21_SHATTER_SPECIAL(mat,i,j,AFTER_CODON) )     {  
          GeneWise21_SHATTER_SPECIAL(mat,i,j,AFTER_CODON) = temp;    
          }  




        /* Finished calculating state DELETE */ 


        /* For state CENTRAL_0M */ 
        /* setting first movement to score */ 
        score = SIG_0_8[MATCH] + (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_MATCH2MATCH])+mat->query->seg[i]->transition[GW_MATCH_BALANCE_5SS]);     
        /* From state INSERT to state CENTRAL_0M */ 
        temp = SIG_0_8[INSERT] + (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_INSERT2MATCH])+mat->query->seg[i]->transition[GW_INSERT_BALANCE_5SS]);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state CENTRAL_0M */ 
        temp = SIG_0_8[DELETE] + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_DELETE2MATCH]);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state CENTRAL_0M to state CENTRAL_0M */ 
        temp = SIG_0_1[CENTRAL_0M] + (mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CENTRAL2CENTRAL]);     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for CENTRAL_0M */ 
        /* Add any movement independant score and put away */ 
         SIG_0_0[CENTRAL_0M] = score;    


        /* Finished calculating state CENTRAL_0M */ 


        /* For state PY_0M */ 
        /* setting first movement to score */ 
        score = SIG_0_1[CENTRAL_0M] + mat->gp21->transition[GP21_CENTRAL2PY];    
        /* From state PY_0M to state PY_0M */ 
        temp = SIG_0_1[PY_0M] + mat->gp21->transition[GP21_PY2PY];   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for PY_0M */ 
        /* Add any movement independant score and put away */ 
         score += mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)];   
         SIG_0_0[PY_0M] = score; 


        /* Finished calculating state PY_0M */ 


        /* For state SPACER_0M */ 
        /* setting first movement to score */ 
        score = SIG_0_1[PY_0M] + mat->gp21->transition[GP21_PY2SPACER];  
        /* From state SPACER_0M to state SPACER_0M */ 
        temp = SIG_0_1[SPACER_0M] + mat->gp21->transition[GP21_SPACER2SPACER];   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for SPACER_0M */ 
        /* Add any movement independant score and put away */ 
         score += mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)];   
         SIG_0_0[SPACER_0M] = score; 


        /* Finished calculating state SPACER_0M */ 


        /* For state CENTRAL_1M */ 
        /* setting first movement to score */ 
        score = SIG_0_9[MATCH] + (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_MATCH2MATCH])+mat->query->seg[i]->transition[GW_MATCH_BALANCE_5SS]);     
        /* From state INSERT to state CENTRAL_1M */ 
        temp = SIG_0_9[INSERT] + (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_INSERT2MATCH])+mat->query->seg[i]->transition[GW_MATCH_BALANCE_5SS]);    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state CENTRAL_1M */ 
        temp = SIG_0_9[DELETE] + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_DELETE2MATCH]);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state CENTRAL_1M to state CENTRAL_1M */ 
        temp = SIG_0_1[CENTRAL_1M] + (mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CENTRAL2CENTRAL]);     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for CENTRAL_1M */ 
        /* Add any movement independant score and put away */ 
         SIG_0_0[CENTRAL_1M] = score;    


        /* Finished calculating state CENTRAL_1M */ 


        /* For state PY_1M */ 
        /* setting first movement to score */ 
        score = SIG_0_1[CENTRAL_1M] + mat->gp21->transition[GP21_CENTRAL2PY];    
        /* From state PY_1M to state PY_1M */ 
        temp = SIG_0_1[PY_1M] + mat->gp21->transition[GP21_PY2PY];   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for PY_1M */ 
        /* Add any movement independant score and put away */ 
         score += mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)];   
         SIG_0_0[PY_1M] = score; 


        /* Finished calculating state PY_1M */ 


        /* For state SPACER_1M */ 
        /* setting first movement to score */ 
        score = SIG_0_1[PY_1M] + mat->gp21->transition[GP21_PY2SPACER];  
        /* From state SPACER_1M to state SPACER_1M */ 
        temp = SIG_0_1[SPACER_1M] + mat->gp21->transition[GP21_SPACER2SPACER];   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for SPACER_1M */ 
        /* Add any movement independant score and put away */ 
         score += mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)];   
         SIG_0_0[SPACER_1M] = score; 


        /* Finished calculating state SPACER_1M */ 


        /* For state CENTRAL_2M */ 
        /* setting first movement to score */ 
        score = SIG_0_10[MATCH] + (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_MATCH2MATCH])+mat->query->seg[i]->transition[GW_MATCH_BALANCE_5SS]);    
        /* From state INSERT to state CENTRAL_2M */ 
        temp = SIG_0_10[INSERT] + (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_INSERT2MATCH])+mat->query->seg[i]->transition[GW_MATCH_BALANCE_5SS]);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state CENTRAL_2M */ 
        temp = SIG_0_10[DELETE] + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_DELETE2MATCH]);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state CENTRAL_2M to state CENTRAL_2M */ 
        temp = SIG_0_1[CENTRAL_2M] + (mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CENTRAL2CENTRAL]);     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for CENTRAL_2M */ 
        /* Add any movement independant score and put away */ 
         SIG_0_0[CENTRAL_2M] = score;    


        /* Finished calculating state CENTRAL_2M */ 


        /* For state PY_2M */ 
        /* setting first movement to score */ 
        score = SIG_0_1[CENTRAL_2M] + mat->gp21->transition[GP21_CENTRAL2PY];    
        /* From state PY_2M to state PY_2M */ 
        temp = SIG_0_1[PY_2M] + mat->gp21->transition[GP21_PY2PY];   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for PY_2M */ 
        /* Add any movement independant score and put away */ 
         score += mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)];   
         SIG_0_0[PY_2M] = score; 


        /* Finished calculating state PY_2M */ 


        /* For state SPACER_2M */ 
        /* setting first movement to score */ 
        score = SIG_0_1[PY_2M] + mat->gp21->transition[GP21_PY2SPACER];  
        /* From state SPACER_2M to state SPACER_2M */ 
        temp = SIG_0_1[SPACER_2M] + mat->gp21->transition[GP21_SPACER2SPACER];   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for SPACER_2M */ 
        /* Add any movement independant score and put away */ 
         score += mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)];   
         SIG_0_0[SPACER_2M] = score; 


        /* Finished calculating state SPACER_2M */ 


        /* For state CENTRAL_0I */ 
        /* setting first movement to score */ 
        score = SIG_0_8[MATCH] + (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_MATCH2INSERT])+mat->query->seg[i]->transition[GW_MATCH_BALANCE_5SS]);    
        /* From state INSERT to state CENTRAL_0I */ 
        temp = SIG_0_8[INSERT] + (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_INSERT2INSERT])+mat->query->seg[i]->transition[GW_INSERT_BALANCE_5SS]);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state CENTRAL_0I */ 
        temp = SIG_0_8[DELETE] + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_DELETE2INSERT]);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state CENTRAL_0I to state CENTRAL_0I */ 
        temp = SIG_0_1[CENTRAL_0I] + (mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CENTRAL2CENTRAL]);     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for CENTRAL_0I */ 
        /* Add any movement independant score and put away */ 
         SIG_0_0[CENTRAL_0I] = score;    


        /* Finished calculating state CENTRAL_0I */ 


        /* For state PY_0I */ 
        /* setting first movement to score */ 
        score = SIG_0_1[CENTRAL_0I] + mat->gp21->transition[GP21_CENTRAL2PY];    
        /* From state PY_0I to state PY_0I */ 
        temp = SIG_0_1[PY_0I] + mat->gp21->transition[GP21_PY2PY];   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for PY_0I */ 
        /* Add any movement independant score and put away */ 
         score += mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)];   
         SIG_0_0[PY_0I] = score; 


        /* Finished calculating state PY_0I */ 


        /* For state SPACER_0I */ 
        /* setting first movement to score */ 
        score = SIG_0_1[PY_0I] + mat->gp21->transition[GP21_PY2SPACER];  
        /* From state SPACER_0I to state SPACER_0I */ 
        temp = SIG_0_1[SPACER_0I] + mat->gp21->transition[GP21_SPACER2SPACER];   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for SPACER_0I */ 
        /* Add any movement independant score and put away */ 
         score += mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)];   
         SIG_0_0[SPACER_0I] = score; 


        /* Finished calculating state SPACER_0I */ 


        /* For state CENTRAL_1I */ 
        /* setting first movement to score */ 
        score = SIG_0_9[MATCH] + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_MATCH2INSERT]);   
        /* From state INSERT to state CENTRAL_1I */ 
        temp = SIG_0_9[INSERT] + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_INSERT2INSERT]);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state CENTRAL_1I */ 
        temp = SIG_0_9[DELETE] + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_DELETE2INSERT]);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state CENTRAL_1I to state CENTRAL_1I */ 
        temp = SIG_0_1[CENTRAL_1I] + (mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CENTRAL2CENTRAL]);     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for CENTRAL_1I */ 
        /* Add any movement independant score and put away */ 
         SIG_0_0[CENTRAL_1I] = score;    


        /* Finished calculating state CENTRAL_1I */ 


        /* For state PY_1I */ 
        /* setting first movement to score */ 
        score = SIG_0_1[CENTRAL_1I] + mat->gp21->transition[GP21_CENTRAL2PY];    
        /* From state PY_1I to state PY_1I */ 
        temp = SIG_0_1[PY_1I] + mat->gp21->transition[GP21_PY2PY];   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for PY_1I */ 
        /* Add any movement independant score and put away */ 
         score += mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)];   
         SIG_0_0[PY_1I] = score; 


        /* Finished calculating state PY_1I */ 


        /* For state SPACER_1I */ 
        /* setting first movement to score */ 
        score = SIG_0_1[PY_1I] + mat->gp21->transition[GP21_PY2SPACER];  
        /* From state SPACER_1I to state SPACER_1I */ 
        temp = SIG_0_1[SPACER_1I] + mat->gp21->transition[GP21_SPACER2SPACER];   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for SPACER_1I */ 
        /* Add any movement independant score and put away */ 
         score += mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)];   
         SIG_0_0[SPACER_1I] = score; 


        /* Finished calculating state SPACER_1I */ 


        /* For state CENTRAL_2I */ 
        /* setting first movement to score */ 
        score = SIG_0_10[MATCH] + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_MATCH2INSERT]);  
        /* From state INSERT to state CENTRAL_2I */ 
        temp = SIG_0_10[INSERT] + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_INSERT2INSERT]);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state CENTRAL_2I */ 
        temp = SIG_0_10[INSERT] + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_DELETE2INSERT]);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state CENTRAL_2I to state CENTRAL_2I */ 
        temp = SIG_0_1[CENTRAL_2I] + (mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CENTRAL2CENTRAL]);     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for CENTRAL_2I */ 
        /* Add any movement independant score and put away */ 
         SIG_0_0[CENTRAL_2I] = score;    


        /* Finished calculating state CENTRAL_2I */ 


        /* For state PY_2I */ 
        /* setting first movement to score */ 
        score = SIG_0_1[CENTRAL_2I] + mat->gp21->transition[GP21_CENTRAL2PY];    
        /* From state PY_2I to state PY_2I */ 
        temp = SIG_0_1[PY_2I] + mat->gp21->transition[GP21_PY2PY];   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for PY_2I */ 
        /* Add any movement independant score and put away */ 
         score += mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)];   
         SIG_0_0[PY_2I] = score; 


        /* Finished calculating state PY_2I */ 


        /* For state SPACER_2I */ 
        /* setting first movement to score */ 
        score = SIG_0_1[PY_2I] + mat->gp21->transition[GP21_PY2SPACER];  
        /* From state SPACER_2I to state SPACER_2I */ 
        temp = SIG_0_1[SPACER_2I] + mat->gp21->transition[GP21_SPACER2SPACER];   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for SPACER_2I */ 
        /* Add any movement independant score and put away */ 
         score += mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)];   
         SIG_0_0[SPACER_2I] = score; 


        /* Finished calculating state SPACER_2I */ 
        }  


      /* Special state START has no special to special movements */ 


      /* Special state END has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GeneWise21_SHATTER_SPECIAL(mat,0,j,END);   


      /* Source AFTER_RND is a special source for END */ 
      /* Has restricted position */ 
      if( j == mat->lenj-1 ) {  
        temp = GeneWise21_SHATTER_SPECIAL(mat,0,j - 1,AFTER_RND) + (0) + (0);    
        if( temp > score )   
          score = temp;  
        }  


      /* Put back score... (now updated!) */ 
      GeneWise21_SHATTER_SPECIAL(mat,0,j,END) = score;   
      /* Finished updating state END */ 




      /* Special state BEFORE_RND has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GeneWise21_SHATTER_SPECIAL(mat,0,j,BEFORE_RND);    


      /* Source START is a special source for BEFORE_RND */ 
      temp = GeneWise21_SHATTER_SPECIAL(mat,0,j - 1,START) + (mat->rndbase->base[CSEQ_GENOMIC_BASE(mat->target,j)]) + (0);   
      if( temp > score ) 
        score = temp;    


      /* Source BEFORE_RND is a special source for BEFORE_RND */ 
      temp = GeneWise21_SHATTER_SPECIAL(mat,0,j - 1,BEFORE_RND) + ((mat->rndbase->base[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_RND2RND])) + (0);    
      if( temp > score ) 
        score = temp;    


      /* Source BEFORE_CODON is a special source for BEFORE_RND */ 
      temp = GeneWise21_SHATTER_SPECIAL(mat,0,j - 1,BEFORE_CODON) + ((mat->rndbase->base[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CDS2RND])) + (0);  
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      GeneWise21_SHATTER_SPECIAL(mat,0,j,BEFORE_RND) = score;    
      /* Finished updating state BEFORE_RND */ 




      /* Special state BEFORE_CODON has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GeneWise21_SHATTER_SPECIAL(mat,0,j,BEFORE_CODON);  


      /* Source BEFORE_RND is a special source for BEFORE_CODON */ 
      temp = GeneWise21_SHATTER_SPECIAL(mat,0,j - 3,BEFORE_RND) + ((mat->rndcodon->codon[CSEQ_GENOMIC_CODON(mat->target,j)]+mat->gp21->transition[GP21_RND2CDS])) + (0);     
      if( temp > score ) 
        score = temp;    


      /* Source BEFORE_CODON is a special source for BEFORE_CODON */ 
      temp = GeneWise21_SHATTER_SPECIAL(mat,0,j - 3,BEFORE_CODON) + ((mat->rndcodon->codon[CSEQ_GENOMIC_CODON(mat->target,j)]+mat->gp21->transition[GP21_CDS2CDS])) + (0);   
      if( temp > score ) 
        score = temp;    


      /* Source BEFORE_SPACER is a special source for BEFORE_CODON */ 
      temp = GeneWise21_SHATTER_SPECIAL(mat,0,j - 5,BEFORE_SPACER) + (((mat->gp21->transition[GP21_SPACER2CDS]+mat->rndcodon->codon[CSEQ_GENOMIC_CODON(mat->target,j)])+CSEQ_GENOMIC_3SS(mat->target,(j-3)))) + (0);     
      if( temp > score ) 
        score = temp;    


      /* Source BEFORE_SPACER is a special source for BEFORE_CODON */ 
      temp = GeneWise21_SHATTER_SPECIAL(mat,0,j - 4,BEFORE_SPACER) + ((mat->gp21->transition[GP21_SPACER2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-2)))) + (0);   
      if( temp > score ) 
        score = temp;    


      /* Source BEFORE_SPACER is a special source for BEFORE_CODON */ 
      temp = GeneWise21_SHATTER_SPECIAL(mat,0,j - 3,BEFORE_SPACER) + ((mat->gp21->transition[GP21_SPACER2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-1)))) + (0);   
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      GeneWise21_SHATTER_SPECIAL(mat,0,j,BEFORE_CODON) = score;  
      /* Finished updating state BEFORE_CODON */ 




      /* Special state BEFORE_CENTRAL has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GeneWise21_SHATTER_SPECIAL(mat,0,j,BEFORE_CENTRAL);    


      /* Source BEFORE_CENTRAL is a special source for BEFORE_CENTRAL */ 
      temp = GeneWise21_SHATTER_SPECIAL(mat,0,j - 1,BEFORE_CENTRAL) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CENTRAL2CENTRAL])) + (0);    
      if( temp > score ) 
        score = temp;    


      /* Source BEFORE_CODON is a special source for BEFORE_CENTRAL */ 
      temp = GeneWise21_SHATTER_SPECIAL(mat,0,j - 8,BEFORE_CODON) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))) + (0);  
      if( temp > score ) 
        score = temp;    


      /* Source BEFORE_CODON is a special source for BEFORE_CENTRAL */ 
      temp = GeneWise21_SHATTER_SPECIAL(mat,0,j - 9,BEFORE_CODON) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))) + (0);  
      if( temp > score ) 
        score = temp;    


      /* Source BEFORE_CODON is a special source for BEFORE_CENTRAL */ 
      temp = GeneWise21_SHATTER_SPECIAL(mat,0,j - 10,BEFORE_CODON) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))) + (0);     
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      GeneWise21_SHATTER_SPECIAL(mat,0,j,BEFORE_CENTRAL) = score;    
      /* Finished updating state BEFORE_CENTRAL */ 




      /* Special state BEFORE_PY_TRACT has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GeneWise21_SHATTER_SPECIAL(mat,0,j,BEFORE_PY_TRACT);   


      /* Source BEFORE_CENTRAL is a special source for BEFORE_PY_TRACT */ 
      temp = GeneWise21_SHATTER_SPECIAL(mat,0,j - 1,BEFORE_CENTRAL) + (mat->gp21->transition[GP21_CENTRAL2PY]) + (mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)]);  
      if( temp > score ) 
        score = temp;    


      /* Source BEFORE_PY_TRACT is a special source for BEFORE_PY_TRACT */ 
      temp = GeneWise21_SHATTER_SPECIAL(mat,0,j - 1,BEFORE_PY_TRACT) + (mat->gp21->transition[GP21_PY2PY]) + (mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)]);  
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      GeneWise21_SHATTER_SPECIAL(mat,0,j,BEFORE_PY_TRACT) = score;   
      /* Finished updating state BEFORE_PY_TRACT */ 




      /* Special state BEFORE_SPACER has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GeneWise21_SHATTER_SPECIAL(mat,0,j,BEFORE_SPACER); 


      /* Source BEFORE_PY_TRACT is a special source for BEFORE_SPACER */ 
      temp = GeneWise21_SHATTER_SPECIAL(mat,0,j - 1,BEFORE_PY_TRACT) + (mat->gp21->transition[GP21_PY2SPACER]) + (mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)]);  
      if( temp > score ) 
        score = temp;    


      /* Source BEFORE_SPACER is a special source for BEFORE_SPACER */ 
      temp = GeneWise21_SHATTER_SPECIAL(mat,0,j - 1,BEFORE_SPACER) + (mat->gp21->transition[GP21_SPACER2SPACER]) + (mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)]);    
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      GeneWise21_SHATTER_SPECIAL(mat,0,j,BEFORE_SPACER) = score; 
      /* Finished updating state BEFORE_SPACER */ 




      /* Special state AFTER_RND has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GeneWise21_SHATTER_SPECIAL(mat,0,j,AFTER_RND); 


      /* Source AFTER_RND is a special source for AFTER_RND */ 
      temp = GeneWise21_SHATTER_SPECIAL(mat,0,j - 1,AFTER_RND) + ((mat->rndbase->base[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_RND2RND])) + (0);     
      if( temp > score ) 
        score = temp;    


      /* Source AFTER_CODON is a special source for AFTER_RND */ 
      temp = GeneWise21_SHATTER_SPECIAL(mat,0,j - 1,AFTER_CODON) + ((mat->rndbase->base[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CDS2RND])) + (0);   
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      GeneWise21_SHATTER_SPECIAL(mat,0,j,AFTER_RND) = score; 
      /* Finished updating state AFTER_RND */ 




      /* Special state AFTER_CODON has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GeneWise21_SHATTER_SPECIAL(mat,0,j,AFTER_CODON);   


      /* Source AFTER_RND is a special source for AFTER_CODON */ 
      temp = GeneWise21_SHATTER_SPECIAL(mat,0,j - 3,AFTER_RND) + ((mat->rndcodon->codon[CSEQ_GENOMIC_CODON(mat->target,j)]+mat->gp21->transition[GP21_RND2CDS])) + (0);  
      if( temp > score ) 
        score = temp;    


      /* Source MATCH for state AFTER_CODON is not special... already calculated */ 
      /* Source INSERT for state AFTER_CODON is not special... already calculated */ 
      /* Source DELETE for state AFTER_CODON is not special... already calculated */ 
      /* Source AFTER_CODON is a special source for AFTER_CODON */ 
      temp = GeneWise21_SHATTER_SPECIAL(mat,0,j - 3,AFTER_CODON) + ((mat->rndcodon->codon[CSEQ_GENOMIC_CODON(mat->target,j)]+mat->gp21->transition[GP21_CDS2CDS])) + (0);    
      if( temp > score ) 
        score = temp;    


      /* Source AFTER_SPACER is a special source for AFTER_CODON */ 
      temp = GeneWise21_SHATTER_SPECIAL(mat,0,j - 5,AFTER_SPACER) + (((mat->gp21->transition[GP21_SPACER2CDS]+mat->rndcodon->codon[CSEQ_GENOMIC_CODON(mat->target,j)])+CSEQ_GENOMIC_3SS(mat->target,(j-3)))) + (0);  
      if( temp > score ) 
        score = temp;    


      /* Source AFTER_SPACER is a special source for AFTER_CODON */ 
      temp = GeneWise21_SHATTER_SPECIAL(mat,0,j - 4,AFTER_SPACER) + ((mat->gp21->transition[GP21_SPACER2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-2)))) + (0);    
      if( temp > score ) 
        score = temp;    


      /* Source AFTER_SPACER is a special source for AFTER_CODON */ 
      temp = GeneWise21_SHATTER_SPECIAL(mat,0,j - 3,AFTER_SPACER) + ((mat->gp21->transition[GP21_SPACER2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-1)))) + (0);    
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      GeneWise21_SHATTER_SPECIAL(mat,0,j,AFTER_CODON) = score;   
      /* Finished updating state AFTER_CODON */ 




      /* Special state AFTER_CENTRAL has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GeneWise21_SHATTER_SPECIAL(mat,0,j,AFTER_CENTRAL); 


      /* Source AFTER_CENTRAL is a special source for AFTER_CENTRAL */ 
      temp = GeneWise21_SHATTER_SPECIAL(mat,0,j - 1,AFTER_CENTRAL) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CENTRAL2CENTRAL])) + (0);     
      if( temp > score ) 
        score = temp;    


      /* Source AFTER_CODON is a special source for AFTER_CENTRAL */ 
      temp = GeneWise21_SHATTER_SPECIAL(mat,0,j - 8,AFTER_CODON) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))) + (0);   
      if( temp > score ) 
        score = temp;    


      /* Source AFTER_CODON is a special source for AFTER_CENTRAL */ 
      temp = GeneWise21_SHATTER_SPECIAL(mat,0,j - 9,AFTER_CODON) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))) + (0);   
      if( temp > score ) 
        score = temp;    


      /* Source AFTER_CODON is a special source for AFTER_CENTRAL */ 
      temp = GeneWise21_SHATTER_SPECIAL(mat,0,j - 10,AFTER_CODON) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))) + (0);  
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      GeneWise21_SHATTER_SPECIAL(mat,0,j,AFTER_CENTRAL) = score; 
      /* Finished updating state AFTER_CENTRAL */ 




      /* Special state AFTER_PY_TRACT has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GeneWise21_SHATTER_SPECIAL(mat,0,j,AFTER_PY_TRACT);    


      /* Source AFTER_CENTRAL is a special source for AFTER_PY_TRACT */ 
      temp = GeneWise21_SHATTER_SPECIAL(mat,0,j - 1,AFTER_CENTRAL) + (mat->gp21->transition[GP21_CENTRAL2PY]) + (mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)]);   
      if( temp > score ) 
        score = temp;    


      /* Source AFTER_PY_TRACT is a special source for AFTER_PY_TRACT */ 
      temp = GeneWise21_SHATTER_SPECIAL(mat,0,j - 1,AFTER_PY_TRACT) + (mat->gp21->transition[GP21_PY2PY]) + (mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)]);   
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      GeneWise21_SHATTER_SPECIAL(mat,0,j,AFTER_PY_TRACT) = score;    
      /* Finished updating state AFTER_PY_TRACT */ 




      /* Special state AFTER_SPACER has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GeneWise21_SHATTER_SPECIAL(mat,0,j,AFTER_SPACER);  


      /* Source AFTER_PY_TRACT is a special source for AFTER_SPACER */ 
      temp = GeneWise21_SHATTER_SPECIAL(mat,0,j - 1,AFTER_PY_TRACT) + (mat->gp21->transition[GP21_PY2SPACER]) + (mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)]);   
      if( temp > score ) 
        score = temp;    


      /* Source AFTER_SPACER is a special source for AFTER_SPACER */ 
      temp = GeneWise21_SHATTER_SPECIAL(mat,0,j - 1,AFTER_SPACER) + (mat->gp21->transition[GP21_SPACER2SPACER]) + (mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)]);     
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      GeneWise21_SHATTER_SPECIAL(mat,0,j,AFTER_SPACER) = score;  
      /* Finished updating state AFTER_SPACER */ 


      }  
    stop_reporting();    
    return TRUE;     
}    


/* Function:  search_GeneWise21(dbsi,out,querydb,targetdb,gp21,rndcodon,rndbase)
 *
 * Descrip:    This function makes a database search of GeneWise21
 *             It uses the dbsi structure to choose which implementation
 *             to use of the database searching. This way at run time you
 *             can switch between single threaded/multi-threaded or hardware
 *
 *
 * Arg:            dbsi [UNKN ] Undocumented argument [DBSearchImpl *]
 * Arg:             out [UNKN ] Undocumented argument [Hscore *]
 * Arg:         querydb [UNKN ] Undocumented argument [GeneWiseDB*]
 * Arg:        targetdb [UNKN ] Undocumented argument [GenomicDB*]
 * Arg:            gp21 [UNKN ] Undocumented argument [GeneParser21Score *]
 * Arg:        rndcodon [UNKN ] Undocumented argument [RandomCodonScore *]
 * Arg:         rndbase [UNKN ] Undocumented argument [RandomModelDNAScore *]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type search_GeneWise21(DBSearchImpl * dbsi,Hscore * out,GeneWiseDB* querydb,GenomicDB* targetdb ,GeneParser21Score * gp21,RandomCodonScore * rndcodon,RandomModelDNAScore * rndbase) 
{
#ifdef PTHREAD   
    int i;   
    int thr_no;  
    pthread_attr_t pat;  
    struct thread_pool_holder_GeneWise21 * holder;   
#endif   
    if( out == NULL )    {  
      warn("Passed in a null Hscore object into search_GeneWise21. Can't process results!"); 
      return SEARCH_ERROR;   
      }  
    if( dbsi == NULL )   {  
      warn("Passed in a null DBSearchImpl object into search_GeneWise21. Can't process results!");   
      return SEARCH_ERROR;   
      }  
    if( dbsi->trace_level > 5 )  
      warn("Asking for trace level of %d in database search for GeneWise21, but it was compiled with a trace level of -2139062144. Not all trace statements can be shown",dbsi->trace_level);    
    switch(dbsi->type)   { /*switch on implementation*/ 
      case DBSearchImpl_Serial : 
        return serial_search_GeneWise21(out,querydb, targetdb ,gp21,rndcodon,rndbase);   
      case DBSearchImpl_Pthreads :   
#ifdef PTHREAD   
        holder = (struct thread_pool_holder_GeneWise21 *) ckalloc(sizeof(struct thread_pool_holder_GeneWise21)); 
        if( holder == NULL )     {  
          warn("Unable to allocated thread pool datastructure...");  
          return SEARCH_ERROR;   
          }  
        holder->out = out;   
        holder->dbsi = dbsi; 
        holder->querydb = querydb;   
        holder->targetdb = targetdb; 
        holder->gp21 = gp21; 
        holder->rndcodon = rndcodon; 
        holder->rndbase = rndbase;   
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
          if( pthread_create(holder->pool+i,&pat,thread_loop_GeneWise21,(void *)holder) )    
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
        warn("You did not specifiy the PTHREAD compile when compiled the C code for GeneWise21");    
#endif /* finished threads */    
      default :  
        warn("database search implementation %s was not provided in the compiled dynamite file from GeneWise21",impl_string_DBSearchImpl(dbsi)); 
        return SEARCH_ERROR; 
      } /* end of switch on implementation */ 


}    


/* Function:  thread_loop_GeneWise21(ptr)
 *
 * Descrip:    Infinite loop code foreach thread for GeneWise21
 *
 *
 * Arg:        ptr [UNKN ] Undocumented argument [void *]
 *
 * Return [UNKN ]  Undocumented return value [void *]
 *
 */
#ifdef PTHREAD   
void * thread_loop_GeneWise21(void * ptr) 
{
    struct thread_pool_holder_GeneWise21 * holder;   
    int db_status;   
    int score;   
    DataScore * ds;  
    GeneWiseScore* query;    
    ComplexSequence* target;     


    holder = (struct thread_pool_holder_GeneWise21 *) ptr;   
    if ( holder->dbsi->trace_level >= 1 )    
      fprintf(holder->dbsi->trace_file,"Entering infinite loop for thread...\n");    


    while(1) { /*Infinite loop over all models*/ 
      /* Get input lock */ 


      if ( holder->dbsi->trace_level >= 2 )  
        fprintf(holder->dbsi->trace_file,"About to get input lock for main reload\n");   
      if( pthread_mutex_lock(&(holder->input_lock))!= 0 )    
        fatal("Error on getting input lock for GeneWise21"); 
      if ( holder->dbsi->trace_level >= 2 )  
        fprintf(holder->dbsi->trace_file,"Got input lock for main reload\n");    


      if( holder->search_has_ended == TRUE ) {  
        if ( holder->dbsi->trace_level >= 2 )    
          fprintf(holder->dbsi->trace_file,"Database search finished for me!...\n"); 
        if( pthread_mutex_unlock(&(holder->input_lock))!= 0 )    
          fatal("Error in releasing input lock for GeneWise21"); 
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
          fatal("Unable to initalise query database in GeneWise21 search");  
        }  
      query = hard_link_GeneWiseScore(holder->query);    
      /* get query information into datascore */ 
      dataentry_add_GeneWiseDB(ds->query,query,holder->querydb); 


      if( holder->target_init == FALSE ) { /*if the db has not been init'd*/ 
        if ( holder->dbsi->trace_level >= 3 )    
          fprintf(holder->dbsi->trace_file,"Starting target database...\n"); 
        target = init_GenomicDB(holder->targetdb,&db_status);    
        holder->target_init = TRUE;  
        } /* end of if the db has not been init'd */ 
      else   { /*Normal reload*/ 
         target = reload_GenomicDB(NULL,holder->targetdb,&db_status);    
        } /* end of Normal reload */ 


      /* Check to see what the reload is like */ 


      if( db_status == DB_RETURN_ERROR ) {  
        fatal("In searching GeneWise21, Reload error on database target, in threads");   
        }  


      if( db_status == DB_RETURN_END)    { /*End of target database*/ 
        /* close target database and schedule it for initalisation by next thread */ 
        close_GenomicDB(NULL,holder->targetdb);  
        holder->target_init = FALSE; 
        if ( holder->dbsi->trace_level >= 2 )    
          fprintf(holder->dbsi->trace_file,"Target Database to be reloaded...\n");   


        /* free'ing the query object */ 
        free_GeneWiseScore(holder->query);   
        /* get the next query object for the next thread */ 
        holder->query = reload_GeneWiseDB(NULL,holder->querydb,&db_status);  
        if( db_status == DB_RETURN_ERROR )   
          fatal("In searching GeneWise21, reload error on database query, in threads");  
        if( db_status == DB_RETURN_END ) { /*last load!*/ 
          /* End of target and query database - finished search! */ 
          close_GeneWiseDB(NULL,holder->querydb);    
          holder->search_has_ended = TRUE;   
          } /* end of last load! */ 


        /* release input mutex */ 
        if ( holder->dbsi->trace_level >= 2 )    
          fprintf(holder->dbsi->trace_file,"Releasing input lock after end of target\n");    
        if( pthread_mutex_unlock(&(holder->input_lock))!= 0 )    
          fatal("Error in releasing input lock for GeneWise21"); 
        continue;    
        } /* end of End of target database */ 
      else   { /*Normal reload*/ 
        if ( holder->dbsi->trace_level >= 2 )    
          fprintf(holder->dbsi->trace_file,"Releasing input lock for normal reload\n");  
        if( pthread_mutex_unlock(&(holder->input_lock))!= 0 )    
          fatal("Error in releasing input lock for GeneWise21"); 
        } /* end of Normal reload */ 
      /* get target information into datascore */ 
      dataentry_add_GenomicDB(ds->target,target,holder->targetdb);   


      /* Now there is a new query/target pair ready for comparison */ 
      if ( holder->dbsi->trace_level >= 1 )  
        fprintf(holder->dbsi->trace_file,"A new pair to be compared...\n");  
      score = score_only_GeneWise21(query, target ,holder->gp21,holder->rndcodon,holder->rndbase);   


      if ( holder->dbsi->trace_level >= 2 )  
        fprintf(holder->dbsi->trace_file,"Getting output lock\n");   
      /* Getting lock on output */ 
      if( pthread_mutex_lock(&(holder->output_lock))!= 0 )   
        fatal("Error on getting output lock for GeneWise21");    
      /* If the score is less than cutoff, schedule the datascore for reuse */ 
      if( should_store_Hscore(holder->out,score) != TRUE)    {  
        free_DataScore(ds);  
        }  
      else   { /*storing score*/ 
        ds->score = score;   
        add_Hscore(holder->out,ds);  
        } /* end of storing score */ 
      if( pthread_mutex_unlock(&(holder->output_lock))!= 0 ) 
        fatal("Error on releasing output lock for GeneWise21");  
      if ( holder->dbsi->trace_level >= 2 )  
        fprintf(holder->dbsi->trace_file,"Released output lock\n");  


      /* Now free database objects */ 
      if ( holder->dbsi->trace_level >= 2 )  
        fprintf(holder->dbsi->trace_file,"About to get input lock for free func\n"); 
      if( pthread_mutex_lock(&(holder->input_lock))!= 0 )    
        fatal("Error on getting input lock for GeneWise21"); 
      if ( holder->dbsi->trace_level >= 2 )  
        fprintf(holder->dbsi->trace_file,"Got input lock for free func\n");  
      free_GeneWiseScore(query); 
      free_ComplexSequence(target);  
      if ( holder->dbsi->trace_level >= 2 )  
        fprintf(holder->dbsi->trace_file,"Releasing input lock after free'ing\n");   
      if( pthread_mutex_unlock(&(holder->input_lock))!= 0 )  
        fatal("Error in releasing input lock for GeneWise21");   
      } /* end of Infinite loop over all models */ 


    if ( holder->dbsi->trace_level >= 1 )    
      fprintf(holder->dbsi->trace_file,"Exiting forever loop\n");    
    return NULL; 
}    


/* Function:  serial_search_GeneWise21(out,querydb,targetdb,gp21,rndcodon,rndbase)
 *
 * Descrip:    This function makes a database search of GeneWise21
 *             It is a single processor implementation
 *
 *
 * Arg:             out [UNKN ] Undocumented argument [Hscore *]
 * Arg:         querydb [UNKN ] Undocumented argument [GeneWiseDB*]
 * Arg:        targetdb [UNKN ] Undocumented argument [GenomicDB*]
 * Arg:            gp21 [UNKN ] Undocumented argument [GeneParser21Score *]
 * Arg:        rndcodon [UNKN ] Undocumented argument [RandomCodonScore *]
 * Arg:         rndbase [UNKN ] Undocumented argument [RandomModelDNAScore *]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
#endif /* PTHREAD */ 
Search_Return_Type serial_search_GeneWise21(Hscore * out,GeneWiseDB* querydb,GenomicDB* targetdb ,GeneParser21Score * gp21,RandomCodonScore * rndcodon,RandomModelDNAScore * rndbase) 
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
      warn("In searching GeneWise21, got a database reload error on the query [query] database");    
      return SEARCH_ERROR;   
      }  
    for(;;)  { /*For all query entries*/ 


      target_pos = 0;    


      target = init_GenomicDB(targetdb,&db_status);  
      if( db_status == DB_RETURN_ERROR )     {  
        warn("In searching GeneWise21, got a database init error on the target [target] database");  
        return SEARCH_ERROR; 
        }  
      for(;;)    { /*For all target entries*/ 


        /* No maximum length - allocated on-the-fly */ 
        score = score_only_GeneWise21(query, target , gp21, rndcodon, rndbase);  
        if( should_store_Hscore(out,score) == TRUE )     { /*if storing datascore*/ 
          ds = new_DataScore_from_storage(out);  
          if( ds == NULL )   {  
            warn("GeneWise21 search had a memory error in allocating a new_DataScore (?a leak somewhere - DataScore is a very small datastructure"); 
            return SEARCH_ERROR; 
            }  
          /* Now: add query/target information to the entry */ 
          dataentry_add_GeneWiseDB(ds->query,query,querydb);     
          dataentry_add_GenomicDB(ds->target,target,targetdb);   
          ds->score = score;     
          add_Hscore(out,ds);    
          } /* end of if storing datascore */ 
        pop_errormsg_stack();    
        push_errormsg_stack("DB searching: just finished [Query Pos: %d] [Target Pos: %d]",query_pos,target_pos);    


         target = reload_GenomicDB(target,targetdb,&db_status);  
        if( db_status == DB_RETURN_ERROR )   {  
          warn("In searching GeneWise21, Reload error on database target, position %d,%d",query_pos,target_pos); 
          return SEARCH_ERROR;   
          }  
        if( db_status == DB_RETURN_END ) 
          break;/* Out of target loop */ 
        target_pos++;    
        } /* end of For all target entries */ 
      close_GenomicDB(target,targetdb);  
       query = reload_GeneWiseDB(query,querydb,&db_status);  
      if( db_status == DB_RETURN_ERROR)  {  
        warn("In searching GeneWise21, Reload error on database query, position %d,%d",query_pos,target_pos);    
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


/* Function:  score_only_GeneWise21(query,target,gp21,rndcodon,rndbase)
 *
 * Descrip:    This function just calculates the score for the matrix
 *             I am pretty sure we can do this better, but hey, for the moment...
 *             It calls /allocate_GeneWise21_only
 *
 *
 * Arg:           query [UNKN ] query data structure [GeneWiseScore*]
 * Arg:          target [UNKN ] target data structure [ComplexSequence*]
 * Arg:            gp21 [UNKN ] Resource [GeneParser21Score *]
 * Arg:        rndcodon [UNKN ] Resource [RandomCodonScore *]
 * Arg:         rndbase [UNKN ] Resource [RandomModelDNAScore *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int score_only_GeneWise21(GeneWiseScore* query,ComplexSequence* target ,GeneParser21Score * gp21,RandomCodonScore * rndcodon,RandomModelDNAScore * rndbase) 
{
    int bestscore = NEGI;    
    int i;   
    int j;   
    int k;   
    GeneWise21 * mat;    


    mat = allocate_GeneWise21_only(query, target , gp21, rndcodon, rndbase); 
    if( mat == NULL )    {  
      warn("Memory allocation error in the db search - unable to communicate to calling function. this spells DIASTER!");    
      return NEGI;   
      }  
    if((mat->basematrix = BaseMatrix_alloc_matrix_and_specials(11,(mat->leni + 1) * 21,11,12)) == NULL)  {  
      warn("Score only matrix for GeneWise21 cannot be allocated, (asking for 10  by %d  cells)",mat->leni*21);  
      mat = free_GeneWise21(mat);    
      return 0;  
      }  
    mat->basematrix->type = BASEMATRIX_TYPE_VERYSMALL;   


    /* Now, initiate matrix */ 
    for(j=0;j<12;j++)    {  
      for(i=(-1);i<mat->leni;i++)    {  
        for(k=0;k<21;k++)    
          GeneWise21_VSMALL_MATRIX(mat,i,j,k) = NEGI;    
        }  
      GeneWise21_VSMALL_SPECIAL(mat,i,j,START) = 0;  
      GeneWise21_VSMALL_SPECIAL(mat,i,j,END) = NEGI; 
      GeneWise21_VSMALL_SPECIAL(mat,i,j,BEFORE_RND) = NEGI;  
      GeneWise21_VSMALL_SPECIAL(mat,i,j,BEFORE_CODON) = NEGI;    
      GeneWise21_VSMALL_SPECIAL(mat,i,j,BEFORE_CENTRAL) = NEGI;  
      GeneWise21_VSMALL_SPECIAL(mat,i,j,BEFORE_PY_TRACT) = NEGI; 
      GeneWise21_VSMALL_SPECIAL(mat,i,j,BEFORE_SPACER) = NEGI;   
      GeneWise21_VSMALL_SPECIAL(mat,i,j,AFTER_RND) = NEGI;   
      GeneWise21_VSMALL_SPECIAL(mat,i,j,AFTER_CODON) = NEGI; 
      GeneWise21_VSMALL_SPECIAL(mat,i,j,AFTER_CENTRAL) = NEGI;   
      GeneWise21_VSMALL_SPECIAL(mat,i,j,AFTER_PY_TRACT) = NEGI;  
      GeneWise21_VSMALL_SPECIAL(mat,i,j,AFTER_SPACER) = NEGI;    
      }  


    /* Ok, lets do-o-o-o-o it */ 


    for(j=0;j<mat->lenj;j++) { /*for all target positions*/ 
      auto int score;    
      auto int temp;     
      for(i=0;i<mat->leni;i++)   { /*for all query positions*/ 


        /* For state MATCH */ 
        /* setting first movement to score */ 
        score = GeneWise21_VSMALL_MATRIX(mat,i-1,j-3,MATCH) + (mat->query->seg[i]->transition[GW_MATCH2MATCH]+mat->query->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]);     
        /* From state INSERT to state MATCH */ 
        temp = GeneWise21_VSMALL_MATRIX(mat,i-1,j-3,INSERT) + (mat->query->seg[i]->transition[GW_INSERT2MATCH]+mat->query->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]);    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state MATCH */ 
        temp = GeneWise21_VSMALL_MATRIX(mat,i-1,j-3,DELETE) + (mat->query->seg[i]->transition[GW_DELETE2MATCH]+mat->query->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]);    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state BEFORE_CODON to state MATCH */ 
        temp = GeneWise21_VSMALL_SPECIAL(mat,i-1,j-3,BEFORE_CODON) + (mat->query->seg[i]->transition[GW_START2MATCH]+mat->query->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SPACER_0M to state MATCH */ 
        temp = GeneWise21_VSMALL_MATRIX(mat,i-1,j-6,SPACER_0M) + (((mat->gp21->transition[GP21_SPACER2CDS]+mat->query->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)])+CSEQ_GENOMIC_3SS(mat->target,(j-3)))+mat->query->seg[i]->transition[GW_MATCH_BALANCE_3SS]);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state PY_0M to state MATCH */ 
        temp = GeneWise21_VSMALL_MATRIX(mat,i-1,j-6,PY_0M) + (((mat->gp21->transition[GP21_PY2CDS]+mat->query->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)])+CSEQ_GENOMIC_3SS(mat->target,(j-3)))+mat->query->seg[i]->transition[GW_MATCH_BALANCE_3SS]);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SPACER_1M to state MATCH */ 
        temp = GeneWise21_VSMALL_MATRIX(mat,i-1,j-5,SPACER_1M) + (mat->gp21->transition[GP21_SPACER2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-2)));   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state PY_1M to state MATCH */ 
        temp = GeneWise21_VSMALL_MATRIX(mat,i-1,j-5,PY_1M) + (mat->gp21->transition[GP21_PY2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-2)));   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SPACER_2M to state MATCH */ 
        temp = GeneWise21_VSMALL_MATRIX(mat,i-1,j-4,SPACER_2M) + (mat->gp21->transition[GP21_SPACER2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-1)));   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state PY_2M to state MATCH */ 
        temp = GeneWise21_VSMALL_MATRIX(mat,i-1,j-4,PY_2M) + (mat->gp21->transition[GP21_PY2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-1)));   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = GeneWise21_VSMALL_MATRIX(mat,i-1,j-2,MATCH) + mat->gp21->transition[GP21_DELETE_1_BASE];  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = GeneWise21_VSMALL_MATRIX(mat,i-1,j-1,MATCH) + mat->gp21->transition[GP21_DELETE_2_BASE];  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = GeneWise21_VSMALL_MATRIX(mat,i-1,j-4,MATCH) + mat->gp21->transition[GP21_INSERT_1_BASE];  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = GeneWise21_VSMALL_MATRIX(mat,i-1,j-5,MATCH) + mat->gp21->transition[GP21_INSERT_2_BASE];  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for MATCH */ 
        /* Add any movement independant score and put away */ 
         score += CSEQ_GENOMIC_CDSPOT(mat->target,j);    
         GeneWise21_VSMALL_MATRIX(mat,i,j,MATCH) = score;    


        /* state MATCH is a source for special AFTER_CODON */ 
        temp = score + (mat->query->seg[i]->transition[GW_MATCH2END]) + (0) ;    
        if( temp > GeneWise21_VSMALL_SPECIAL(mat,i,j,AFTER_CODON) )  {  
          GeneWise21_VSMALL_SPECIAL(mat,i,j,AFTER_CODON) = temp;     
          }  




        /* Finished calculating state MATCH */ 


        /* For state INSERT */ 
        /* setting first movement to score */ 
        score = GeneWise21_VSMALL_MATRIX(mat,i-0,j-3,MATCH) + (mat->query->seg[i]->transition[GW_MATCH2INSERT]+mat->query->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)]);   
        /* From state INSERT to state INSERT */ 
        temp = GeneWise21_VSMALL_MATRIX(mat,i-0,j-3,INSERT) + (mat->query->seg[i]->transition[GW_INSERT2INSERT]+mat->query->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)]);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state INSERT */ 
        temp = GeneWise21_VSMALL_MATRIX(mat,i-0,j-3,DELETE) + (mat->query->seg[i]->transition[GW_DELETE2INSERT]+mat->query->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)]);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state BEFORE_CODON to state INSERT */ 
        temp = GeneWise21_VSMALL_SPECIAL(mat,i-0,j-3,BEFORE_CODON) + (mat->query->seg[i]->transition[GW_START2INSERT]+mat->query->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)]);    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SPACER_0I to state INSERT */ 
        temp = GeneWise21_VSMALL_MATRIX(mat,i-1,j-6,SPACER_0I) + (((mat->gp21->transition[GP21_SPACER2CDS]+mat->query->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)])+CSEQ_GENOMIC_3SS(mat->target,(j-3)))+mat->query->seg[i]->transition[GW_MATCH_BALANCE_3SS]);    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SPACER_1I to state INSERT */ 
        temp = GeneWise21_VSMALL_MATRIX(mat,i-1,j-5,SPACER_1I) + (mat->gp21->transition[GP21_SPACER2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-2)));   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SPACER_2I to state INSERT */ 
        temp = GeneWise21_VSMALL_MATRIX(mat,i-1,j-4,SPACER_2I) + (mat->gp21->transition[GP21_SPACER2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-1)));   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state PY_0I to state INSERT */ 
        temp = GeneWise21_VSMALL_MATRIX(mat,i-1,j-6,PY_0I) + (((mat->gp21->transition[GP21_PY2CDS]+mat->query->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)])+CSEQ_GENOMIC_3SS(mat->target,(j-3)))+mat->query->seg[i]->transition[GW_MATCH_BALANCE_3SS]);    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state PY_1I to state INSERT */ 
        temp = GeneWise21_VSMALL_MATRIX(mat,i-1,j-5,PY_1I) + (mat->gp21->transition[GP21_PY2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-2)));   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state PY_2I to state INSERT */ 
        temp = GeneWise21_VSMALL_MATRIX(mat,i-1,j-4,PY_2I) + (mat->gp21->transition[GP21_PY2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-1)));   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state INSERT */ 
        temp = GeneWise21_VSMALL_MATRIX(mat,i-1,j-2,INSERT) + mat->gp21->transition[GP21_DELETE_1_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state INSERT */ 
        temp = GeneWise21_VSMALL_MATRIX(mat,i-1,j-1,INSERT) + mat->gp21->transition[GP21_DELETE_2_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state INSERT */ 
        temp = GeneWise21_VSMALL_MATRIX(mat,i-1,j-4,INSERT) + mat->gp21->transition[GP21_INSERT_1_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state INSERT */ 
        temp = GeneWise21_VSMALL_MATRIX(mat,i-1,j-5,INSERT) + mat->gp21->transition[GP21_INSERT_2_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for INSERT */ 
        /* Add any movement independant score and put away */ 
         score += CSEQ_GENOMIC_CDSPOT(mat->target,j);    
         GeneWise21_VSMALL_MATRIX(mat,i,j,INSERT) = score;   


        /* state INSERT is a source for special AFTER_CODON */ 
        temp = score + (mat->query->seg[i]->transition[GW_INSERT2END]) + (0) ;   
        if( temp > GeneWise21_VSMALL_SPECIAL(mat,i,j,AFTER_CODON) )  {  
          GeneWise21_VSMALL_SPECIAL(mat,i,j,AFTER_CODON) = temp;     
          }  




        /* Finished calculating state INSERT */ 


        /* For state DELETE */ 
        /* setting first movement to score */ 
        score = GeneWise21_VSMALL_MATRIX(mat,i-1,j-0,MATCH) + mat->query->seg[i]->transition[GW_MATCH2DELETE];   
        /* From state INSERT to state DELETE */ 
        temp = GeneWise21_VSMALL_MATRIX(mat,i-1,j-0,INSERT) + mat->query->seg[i]->transition[GW_INSERT2DELETE];  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state DELETE */ 
        temp = GeneWise21_VSMALL_MATRIX(mat,i-1,j-0,DELETE) + mat->query->seg[i]->transition[GW_DELETE2DELETE];  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state BEFORE_CODON to state DELETE */ 
        temp = GeneWise21_VSMALL_SPECIAL(mat,i-1,j-0,BEFORE_CODON) + mat->query->seg[i]->transition[GW_START2DELETE];    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for DELETE */ 
        /* Add any movement independant score and put away */ 
         GeneWise21_VSMALL_MATRIX(mat,i,j,DELETE) = score;   


        /* state DELETE is a source for special AFTER_CODON */ 
        temp = score + (mat->query->seg[i]->transition[GW_DELETE2END]) + (0) ;   
        if( temp > GeneWise21_VSMALL_SPECIAL(mat,i,j,AFTER_CODON) )  {  
          GeneWise21_VSMALL_SPECIAL(mat,i,j,AFTER_CODON) = temp;     
          }  




        /* Finished calculating state DELETE */ 


        /* For state CENTRAL_0M */ 
        /* setting first movement to score */ 
        score = GeneWise21_VSMALL_MATRIX(mat,i-0,j-8,MATCH) + (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_MATCH2MATCH])+mat->query->seg[i]->transition[GW_MATCH_BALANCE_5SS]);    
        /* From state INSERT to state CENTRAL_0M */ 
        temp = GeneWise21_VSMALL_MATRIX(mat,i-0,j-8,INSERT) + (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_INSERT2MATCH])+mat->query->seg[i]->transition[GW_INSERT_BALANCE_5SS]);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state CENTRAL_0M */ 
        temp = GeneWise21_VSMALL_MATRIX(mat,i-0,j-8,DELETE) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_DELETE2MATCH]);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state CENTRAL_0M to state CENTRAL_0M */ 
        temp = GeneWise21_VSMALL_MATRIX(mat,i-0,j-1,CENTRAL_0M) + (mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CENTRAL2CENTRAL]);    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for CENTRAL_0M */ 
        /* Add any movement independant score and put away */ 
         GeneWise21_VSMALL_MATRIX(mat,i,j,CENTRAL_0M) = score;   


        /* Finished calculating state CENTRAL_0M */ 


        /* For state PY_0M */ 
        /* setting first movement to score */ 
        score = GeneWise21_VSMALL_MATRIX(mat,i-0,j-1,CENTRAL_0M) + mat->gp21->transition[GP21_CENTRAL2PY];   
        /* From state PY_0M to state PY_0M */ 
        temp = GeneWise21_VSMALL_MATRIX(mat,i-0,j-1,PY_0M) + mat->gp21->transition[GP21_PY2PY];  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for PY_0M */ 
        /* Add any movement independant score and put away */ 
         score += mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)];   
         GeneWise21_VSMALL_MATRIX(mat,i,j,PY_0M) = score;    


        /* Finished calculating state PY_0M */ 


        /* For state SPACER_0M */ 
        /* setting first movement to score */ 
        score = GeneWise21_VSMALL_MATRIX(mat,i-0,j-1,PY_0M) + mat->gp21->transition[GP21_PY2SPACER];     
        /* From state SPACER_0M to state SPACER_0M */ 
        temp = GeneWise21_VSMALL_MATRIX(mat,i-0,j-1,SPACER_0M) + mat->gp21->transition[GP21_SPACER2SPACER];  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for SPACER_0M */ 
        /* Add any movement independant score and put away */ 
         score += mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)];   
         GeneWise21_VSMALL_MATRIX(mat,i,j,SPACER_0M) = score;    


        /* Finished calculating state SPACER_0M */ 


        /* For state CENTRAL_1M */ 
        /* setting first movement to score */ 
        score = GeneWise21_VSMALL_MATRIX(mat,i-0,j-9,MATCH) + (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_MATCH2MATCH])+mat->query->seg[i]->transition[GW_MATCH_BALANCE_5SS]);    
        /* From state INSERT to state CENTRAL_1M */ 
        temp = GeneWise21_VSMALL_MATRIX(mat,i-0,j-9,INSERT) + (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_INSERT2MATCH])+mat->query->seg[i]->transition[GW_MATCH_BALANCE_5SS]);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state CENTRAL_1M */ 
        temp = GeneWise21_VSMALL_MATRIX(mat,i-0,j-9,DELETE) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_DELETE2MATCH]);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state CENTRAL_1M to state CENTRAL_1M */ 
        temp = GeneWise21_VSMALL_MATRIX(mat,i-0,j-1,CENTRAL_1M) + (mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CENTRAL2CENTRAL]);    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for CENTRAL_1M */ 
        /* Add any movement independant score and put away */ 
         GeneWise21_VSMALL_MATRIX(mat,i,j,CENTRAL_1M) = score;   


        /* Finished calculating state CENTRAL_1M */ 


        /* For state PY_1M */ 
        /* setting first movement to score */ 
        score = GeneWise21_VSMALL_MATRIX(mat,i-0,j-1,CENTRAL_1M) + mat->gp21->transition[GP21_CENTRAL2PY];   
        /* From state PY_1M to state PY_1M */ 
        temp = GeneWise21_VSMALL_MATRIX(mat,i-0,j-1,PY_1M) + mat->gp21->transition[GP21_PY2PY];  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for PY_1M */ 
        /* Add any movement independant score and put away */ 
         score += mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)];   
         GeneWise21_VSMALL_MATRIX(mat,i,j,PY_1M) = score;    


        /* Finished calculating state PY_1M */ 


        /* For state SPACER_1M */ 
        /* setting first movement to score */ 
        score = GeneWise21_VSMALL_MATRIX(mat,i-0,j-1,PY_1M) + mat->gp21->transition[GP21_PY2SPACER];     
        /* From state SPACER_1M to state SPACER_1M */ 
        temp = GeneWise21_VSMALL_MATRIX(mat,i-0,j-1,SPACER_1M) + mat->gp21->transition[GP21_SPACER2SPACER];  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for SPACER_1M */ 
        /* Add any movement independant score and put away */ 
         score += mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)];   
         GeneWise21_VSMALL_MATRIX(mat,i,j,SPACER_1M) = score;    


        /* Finished calculating state SPACER_1M */ 


        /* For state CENTRAL_2M */ 
        /* setting first movement to score */ 
        score = GeneWise21_VSMALL_MATRIX(mat,i-0,j-10,MATCH) + (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_MATCH2MATCH])+mat->query->seg[i]->transition[GW_MATCH_BALANCE_5SS]);   
        /* From state INSERT to state CENTRAL_2M */ 
        temp = GeneWise21_VSMALL_MATRIX(mat,i-0,j-10,INSERT) + (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_INSERT2MATCH])+mat->query->seg[i]->transition[GW_MATCH_BALANCE_5SS]);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state CENTRAL_2M */ 
        temp = GeneWise21_VSMALL_MATRIX(mat,i-0,j-10,DELETE) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_DELETE2MATCH]);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state CENTRAL_2M to state CENTRAL_2M */ 
        temp = GeneWise21_VSMALL_MATRIX(mat,i-0,j-1,CENTRAL_2M) + (mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CENTRAL2CENTRAL]);    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for CENTRAL_2M */ 
        /* Add any movement independant score and put away */ 
         GeneWise21_VSMALL_MATRIX(mat,i,j,CENTRAL_2M) = score;   


        /* Finished calculating state CENTRAL_2M */ 


        /* For state PY_2M */ 
        /* setting first movement to score */ 
        score = GeneWise21_VSMALL_MATRIX(mat,i-0,j-1,CENTRAL_2M) + mat->gp21->transition[GP21_CENTRAL2PY];   
        /* From state PY_2M to state PY_2M */ 
        temp = GeneWise21_VSMALL_MATRIX(mat,i-0,j-1,PY_2M) + mat->gp21->transition[GP21_PY2PY];  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for PY_2M */ 
        /* Add any movement independant score and put away */ 
         score += mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)];   
         GeneWise21_VSMALL_MATRIX(mat,i,j,PY_2M) = score;    


        /* Finished calculating state PY_2M */ 


        /* For state SPACER_2M */ 
        /* setting first movement to score */ 
        score = GeneWise21_VSMALL_MATRIX(mat,i-0,j-1,PY_2M) + mat->gp21->transition[GP21_PY2SPACER];     
        /* From state SPACER_2M to state SPACER_2M */ 
        temp = GeneWise21_VSMALL_MATRIX(mat,i-0,j-1,SPACER_2M) + mat->gp21->transition[GP21_SPACER2SPACER];  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for SPACER_2M */ 
        /* Add any movement independant score and put away */ 
         score += mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)];   
         GeneWise21_VSMALL_MATRIX(mat,i,j,SPACER_2M) = score;    


        /* Finished calculating state SPACER_2M */ 


        /* For state CENTRAL_0I */ 
        /* setting first movement to score */ 
        score = GeneWise21_VSMALL_MATRIX(mat,i-0,j-8,MATCH) + (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_MATCH2INSERT])+mat->query->seg[i]->transition[GW_MATCH_BALANCE_5SS]);   
        /* From state INSERT to state CENTRAL_0I */ 
        temp = GeneWise21_VSMALL_MATRIX(mat,i-0,j-8,INSERT) + (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_INSERT2INSERT])+mat->query->seg[i]->transition[GW_INSERT_BALANCE_5SS]);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state CENTRAL_0I */ 
        temp = GeneWise21_VSMALL_MATRIX(mat,i-0,j-8,DELETE) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_DELETE2INSERT]);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state CENTRAL_0I to state CENTRAL_0I */ 
        temp = GeneWise21_VSMALL_MATRIX(mat,i-0,j-1,CENTRAL_0I) + (mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CENTRAL2CENTRAL]);    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for CENTRAL_0I */ 
        /* Add any movement independant score and put away */ 
         GeneWise21_VSMALL_MATRIX(mat,i,j,CENTRAL_0I) = score;   


        /* Finished calculating state CENTRAL_0I */ 


        /* For state PY_0I */ 
        /* setting first movement to score */ 
        score = GeneWise21_VSMALL_MATRIX(mat,i-0,j-1,CENTRAL_0I) + mat->gp21->transition[GP21_CENTRAL2PY];   
        /* From state PY_0I to state PY_0I */ 
        temp = GeneWise21_VSMALL_MATRIX(mat,i-0,j-1,PY_0I) + mat->gp21->transition[GP21_PY2PY];  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for PY_0I */ 
        /* Add any movement independant score and put away */ 
         score += mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)];   
         GeneWise21_VSMALL_MATRIX(mat,i,j,PY_0I) = score;    


        /* Finished calculating state PY_0I */ 


        /* For state SPACER_0I */ 
        /* setting first movement to score */ 
        score = GeneWise21_VSMALL_MATRIX(mat,i-0,j-1,PY_0I) + mat->gp21->transition[GP21_PY2SPACER];     
        /* From state SPACER_0I to state SPACER_0I */ 
        temp = GeneWise21_VSMALL_MATRIX(mat,i-0,j-1,SPACER_0I) + mat->gp21->transition[GP21_SPACER2SPACER];  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for SPACER_0I */ 
        /* Add any movement independant score and put away */ 
         score += mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)];   
         GeneWise21_VSMALL_MATRIX(mat,i,j,SPACER_0I) = score;    


        /* Finished calculating state SPACER_0I */ 


        /* For state CENTRAL_1I */ 
        /* setting first movement to score */ 
        score = GeneWise21_VSMALL_MATRIX(mat,i-0,j-9,MATCH) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_MATCH2INSERT]);  
        /* From state INSERT to state CENTRAL_1I */ 
        temp = GeneWise21_VSMALL_MATRIX(mat,i-0,j-9,INSERT) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_INSERT2INSERT]);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state CENTRAL_1I */ 
        temp = GeneWise21_VSMALL_MATRIX(mat,i-0,j-9,DELETE) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_DELETE2INSERT]);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state CENTRAL_1I to state CENTRAL_1I */ 
        temp = GeneWise21_VSMALL_MATRIX(mat,i-0,j-1,CENTRAL_1I) + (mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CENTRAL2CENTRAL]);    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for CENTRAL_1I */ 
        /* Add any movement independant score and put away */ 
         GeneWise21_VSMALL_MATRIX(mat,i,j,CENTRAL_1I) = score;   


        /* Finished calculating state CENTRAL_1I */ 


        /* For state PY_1I */ 
        /* setting first movement to score */ 
        score = GeneWise21_VSMALL_MATRIX(mat,i-0,j-1,CENTRAL_1I) + mat->gp21->transition[GP21_CENTRAL2PY];   
        /* From state PY_1I to state PY_1I */ 
        temp = GeneWise21_VSMALL_MATRIX(mat,i-0,j-1,PY_1I) + mat->gp21->transition[GP21_PY2PY];  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for PY_1I */ 
        /* Add any movement independant score and put away */ 
         score += mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)];   
         GeneWise21_VSMALL_MATRIX(mat,i,j,PY_1I) = score;    


        /* Finished calculating state PY_1I */ 


        /* For state SPACER_1I */ 
        /* setting first movement to score */ 
        score = GeneWise21_VSMALL_MATRIX(mat,i-0,j-1,PY_1I) + mat->gp21->transition[GP21_PY2SPACER];     
        /* From state SPACER_1I to state SPACER_1I */ 
        temp = GeneWise21_VSMALL_MATRIX(mat,i-0,j-1,SPACER_1I) + mat->gp21->transition[GP21_SPACER2SPACER];  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for SPACER_1I */ 
        /* Add any movement independant score and put away */ 
         score += mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)];   
         GeneWise21_VSMALL_MATRIX(mat,i,j,SPACER_1I) = score;    


        /* Finished calculating state SPACER_1I */ 


        /* For state CENTRAL_2I */ 
        /* setting first movement to score */ 
        score = GeneWise21_VSMALL_MATRIX(mat,i-0,j-10,MATCH) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_MATCH2INSERT]);     
        /* From state INSERT to state CENTRAL_2I */ 
        temp = GeneWise21_VSMALL_MATRIX(mat,i-0,j-10,INSERT) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_INSERT2INSERT]);    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state CENTRAL_2I */ 
        temp = GeneWise21_VSMALL_MATRIX(mat,i-0,j-10,INSERT) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_DELETE2INSERT]);    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state CENTRAL_2I to state CENTRAL_2I */ 
        temp = GeneWise21_VSMALL_MATRIX(mat,i-0,j-1,CENTRAL_2I) + (mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CENTRAL2CENTRAL]);    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for CENTRAL_2I */ 
        /* Add any movement independant score and put away */ 
         GeneWise21_VSMALL_MATRIX(mat,i,j,CENTRAL_2I) = score;   


        /* Finished calculating state CENTRAL_2I */ 


        /* For state PY_2I */ 
        /* setting first movement to score */ 
        score = GeneWise21_VSMALL_MATRIX(mat,i-0,j-1,CENTRAL_2I) + mat->gp21->transition[GP21_CENTRAL2PY];   
        /* From state PY_2I to state PY_2I */ 
        temp = GeneWise21_VSMALL_MATRIX(mat,i-0,j-1,PY_2I) + mat->gp21->transition[GP21_PY2PY];  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for PY_2I */ 
        /* Add any movement independant score and put away */ 
         score += mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)];   
         GeneWise21_VSMALL_MATRIX(mat,i,j,PY_2I) = score;    


        /* Finished calculating state PY_2I */ 


        /* For state SPACER_2I */ 
        /* setting first movement to score */ 
        score = GeneWise21_VSMALL_MATRIX(mat,i-0,j-1,PY_2I) + mat->gp21->transition[GP21_PY2SPACER];     
        /* From state SPACER_2I to state SPACER_2I */ 
        temp = GeneWise21_VSMALL_MATRIX(mat,i-0,j-1,SPACER_2I) + mat->gp21->transition[GP21_SPACER2SPACER];  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for SPACER_2I */ 
        /* Add any movement independant score and put away */ 
         score += mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)];   
         GeneWise21_VSMALL_MATRIX(mat,i,j,SPACER_2I) = score;    


        /* Finished calculating state SPACER_2I */ 
        } /* end of for all query positions */ 




      /* Special state START has no special to special movements */ 


      /* Special state END has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GeneWise21_VSMALL_SPECIAL(mat,0,j,END);    


      /* Source AFTER_RND is a special source for END */ 
      /* Has restricted position */ 
      if( j == mat->lenj-1 ) {  
        temp = GeneWise21_VSMALL_SPECIAL(mat,0,j - 1,AFTER_RND) + (0) + (0);     
        if( temp > score )   
          score = temp;  
        }  


      /* Put back score... (now updated!) */ 
      GeneWise21_VSMALL_SPECIAL(mat,0,j,END) = score;    
      /* Finished updating state END */ 




      /* Special state BEFORE_RND has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GeneWise21_VSMALL_SPECIAL(mat,0,j,BEFORE_RND); 


      /* Source START is a special source for BEFORE_RND */ 
      temp = GeneWise21_VSMALL_SPECIAL(mat,0,j - 1,START) + (mat->rndbase->base[CSEQ_GENOMIC_BASE(mat->target,j)]) + (0);    
      if( temp > score ) 
        score = temp;    


      /* Source BEFORE_RND is a special source for BEFORE_RND */ 
      temp = GeneWise21_VSMALL_SPECIAL(mat,0,j - 1,BEFORE_RND) + ((mat->rndbase->base[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_RND2RND])) + (0);     
      if( temp > score ) 
        score = temp;    


      /* Source BEFORE_CODON is a special source for BEFORE_RND */ 
      temp = GeneWise21_VSMALL_SPECIAL(mat,0,j - 1,BEFORE_CODON) + ((mat->rndbase->base[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CDS2RND])) + (0);   
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      GeneWise21_VSMALL_SPECIAL(mat,0,j,BEFORE_RND) = score; 
      /* Finished updating state BEFORE_RND */ 




      /* Special state BEFORE_CODON has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GeneWise21_VSMALL_SPECIAL(mat,0,j,BEFORE_CODON);   


      /* Source BEFORE_RND is a special source for BEFORE_CODON */ 
      temp = GeneWise21_VSMALL_SPECIAL(mat,0,j - 3,BEFORE_RND) + ((mat->rndcodon->codon[CSEQ_GENOMIC_CODON(mat->target,j)]+mat->gp21->transition[GP21_RND2CDS])) + (0);  
      if( temp > score ) 
        score = temp;    


      /* Source BEFORE_CODON is a special source for BEFORE_CODON */ 
      temp = GeneWise21_VSMALL_SPECIAL(mat,0,j - 3,BEFORE_CODON) + ((mat->rndcodon->codon[CSEQ_GENOMIC_CODON(mat->target,j)]+mat->gp21->transition[GP21_CDS2CDS])) + (0);    
      if( temp > score ) 
        score = temp;    


      /* Source BEFORE_SPACER is a special source for BEFORE_CODON */ 
      temp = GeneWise21_VSMALL_SPECIAL(mat,0,j - 5,BEFORE_SPACER) + (((mat->gp21->transition[GP21_SPACER2CDS]+mat->rndcodon->codon[CSEQ_GENOMIC_CODON(mat->target,j)])+CSEQ_GENOMIC_3SS(mat->target,(j-3)))) + (0);  
      if( temp > score ) 
        score = temp;    


      /* Source BEFORE_SPACER is a special source for BEFORE_CODON */ 
      temp = GeneWise21_VSMALL_SPECIAL(mat,0,j - 4,BEFORE_SPACER) + ((mat->gp21->transition[GP21_SPACER2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-2)))) + (0);    
      if( temp > score ) 
        score = temp;    


      /* Source BEFORE_SPACER is a special source for BEFORE_CODON */ 
      temp = GeneWise21_VSMALL_SPECIAL(mat,0,j - 3,BEFORE_SPACER) + ((mat->gp21->transition[GP21_SPACER2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-1)))) + (0);    
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      GeneWise21_VSMALL_SPECIAL(mat,0,j,BEFORE_CODON) = score;   
      /* Finished updating state BEFORE_CODON */ 




      /* Special state BEFORE_CENTRAL has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GeneWise21_VSMALL_SPECIAL(mat,0,j,BEFORE_CENTRAL); 


      /* Source BEFORE_CENTRAL is a special source for BEFORE_CENTRAL */ 
      temp = GeneWise21_VSMALL_SPECIAL(mat,0,j - 1,BEFORE_CENTRAL) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CENTRAL2CENTRAL])) + (0);     
      if( temp > score ) 
        score = temp;    


      /* Source BEFORE_CODON is a special source for BEFORE_CENTRAL */ 
      temp = GeneWise21_VSMALL_SPECIAL(mat,0,j - 8,BEFORE_CODON) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))) + (0);   
      if( temp > score ) 
        score = temp;    


      /* Source BEFORE_CODON is a special source for BEFORE_CENTRAL */ 
      temp = GeneWise21_VSMALL_SPECIAL(mat,0,j - 9,BEFORE_CODON) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))) + (0);   
      if( temp > score ) 
        score = temp;    


      /* Source BEFORE_CODON is a special source for BEFORE_CENTRAL */ 
      temp = GeneWise21_VSMALL_SPECIAL(mat,0,j - 10,BEFORE_CODON) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))) + (0);  
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      GeneWise21_VSMALL_SPECIAL(mat,0,j,BEFORE_CENTRAL) = score; 
      /* Finished updating state BEFORE_CENTRAL */ 




      /* Special state BEFORE_PY_TRACT has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GeneWise21_VSMALL_SPECIAL(mat,0,j,BEFORE_PY_TRACT);    


      /* Source BEFORE_CENTRAL is a special source for BEFORE_PY_TRACT */ 
      temp = GeneWise21_VSMALL_SPECIAL(mat,0,j - 1,BEFORE_CENTRAL) + (mat->gp21->transition[GP21_CENTRAL2PY]) + (mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)]);   
      if( temp > score ) 
        score = temp;    


      /* Source BEFORE_PY_TRACT is a special source for BEFORE_PY_TRACT */ 
      temp = GeneWise21_VSMALL_SPECIAL(mat,0,j - 1,BEFORE_PY_TRACT) + (mat->gp21->transition[GP21_PY2PY]) + (mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)]);   
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      GeneWise21_VSMALL_SPECIAL(mat,0,j,BEFORE_PY_TRACT) = score;    
      /* Finished updating state BEFORE_PY_TRACT */ 




      /* Special state BEFORE_SPACER has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GeneWise21_VSMALL_SPECIAL(mat,0,j,BEFORE_SPACER);  


      /* Source BEFORE_PY_TRACT is a special source for BEFORE_SPACER */ 
      temp = GeneWise21_VSMALL_SPECIAL(mat,0,j - 1,BEFORE_PY_TRACT) + (mat->gp21->transition[GP21_PY2SPACER]) + (mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)]);   
      if( temp > score ) 
        score = temp;    


      /* Source BEFORE_SPACER is a special source for BEFORE_SPACER */ 
      temp = GeneWise21_VSMALL_SPECIAL(mat,0,j - 1,BEFORE_SPACER) + (mat->gp21->transition[GP21_SPACER2SPACER]) + (mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)]);     
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      GeneWise21_VSMALL_SPECIAL(mat,0,j,BEFORE_SPACER) = score;  
      /* Finished updating state BEFORE_SPACER */ 




      /* Special state AFTER_RND has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GeneWise21_VSMALL_SPECIAL(mat,0,j,AFTER_RND);  


      /* Source AFTER_RND is a special source for AFTER_RND */ 
      temp = GeneWise21_VSMALL_SPECIAL(mat,0,j - 1,AFTER_RND) + ((mat->rndbase->base[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_RND2RND])) + (0);  
      if( temp > score ) 
        score = temp;    


      /* Source AFTER_CODON is a special source for AFTER_RND */ 
      temp = GeneWise21_VSMALL_SPECIAL(mat,0,j - 1,AFTER_CODON) + ((mat->rndbase->base[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CDS2RND])) + (0);    
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      GeneWise21_VSMALL_SPECIAL(mat,0,j,AFTER_RND) = score;  
      /* Finished updating state AFTER_RND */ 




      /* Special state AFTER_CODON has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GeneWise21_VSMALL_SPECIAL(mat,0,j,AFTER_CODON);    


      /* Source AFTER_RND is a special source for AFTER_CODON */ 
      temp = GeneWise21_VSMALL_SPECIAL(mat,0,j - 3,AFTER_RND) + ((mat->rndcodon->codon[CSEQ_GENOMIC_CODON(mat->target,j)]+mat->gp21->transition[GP21_RND2CDS])) + (0);   
      if( temp > score ) 
        score = temp;    


      /* Source MATCH for state AFTER_CODON is not special... already calculated */ 
      /* Source INSERT for state AFTER_CODON is not special... already calculated */ 
      /* Source DELETE for state AFTER_CODON is not special... already calculated */ 
      /* Source AFTER_CODON is a special source for AFTER_CODON */ 
      temp = GeneWise21_VSMALL_SPECIAL(mat,0,j - 3,AFTER_CODON) + ((mat->rndcodon->codon[CSEQ_GENOMIC_CODON(mat->target,j)]+mat->gp21->transition[GP21_CDS2CDS])) + (0);     
      if( temp > score ) 
        score = temp;    


      /* Source AFTER_SPACER is a special source for AFTER_CODON */ 
      temp = GeneWise21_VSMALL_SPECIAL(mat,0,j - 5,AFTER_SPACER) + (((mat->gp21->transition[GP21_SPACER2CDS]+mat->rndcodon->codon[CSEQ_GENOMIC_CODON(mat->target,j)])+CSEQ_GENOMIC_3SS(mat->target,(j-3)))) + (0);   
      if( temp > score ) 
        score = temp;    


      /* Source AFTER_SPACER is a special source for AFTER_CODON */ 
      temp = GeneWise21_VSMALL_SPECIAL(mat,0,j - 4,AFTER_SPACER) + ((mat->gp21->transition[GP21_SPACER2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-2)))) + (0);     
      if( temp > score ) 
        score = temp;    


      /* Source AFTER_SPACER is a special source for AFTER_CODON */ 
      temp = GeneWise21_VSMALL_SPECIAL(mat,0,j - 3,AFTER_SPACER) + ((mat->gp21->transition[GP21_SPACER2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-1)))) + (0);     
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      GeneWise21_VSMALL_SPECIAL(mat,0,j,AFTER_CODON) = score;    
      /* Finished updating state AFTER_CODON */ 




      /* Special state AFTER_CENTRAL has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GeneWise21_VSMALL_SPECIAL(mat,0,j,AFTER_CENTRAL);  


      /* Source AFTER_CENTRAL is a special source for AFTER_CENTRAL */ 
      temp = GeneWise21_VSMALL_SPECIAL(mat,0,j - 1,AFTER_CENTRAL) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CENTRAL2CENTRAL])) + (0);  
      if( temp > score ) 
        score = temp;    


      /* Source AFTER_CODON is a special source for AFTER_CENTRAL */ 
      temp = GeneWise21_VSMALL_SPECIAL(mat,0,j - 8,AFTER_CODON) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))) + (0);    
      if( temp > score ) 
        score = temp;    


      /* Source AFTER_CODON is a special source for AFTER_CENTRAL */ 
      temp = GeneWise21_VSMALL_SPECIAL(mat,0,j - 9,AFTER_CODON) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))) + (0);    
      if( temp > score ) 
        score = temp;    


      /* Source AFTER_CODON is a special source for AFTER_CENTRAL */ 
      temp = GeneWise21_VSMALL_SPECIAL(mat,0,j - 10,AFTER_CODON) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))) + (0);   
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      GeneWise21_VSMALL_SPECIAL(mat,0,j,AFTER_CENTRAL) = score;  
      /* Finished updating state AFTER_CENTRAL */ 




      /* Special state AFTER_PY_TRACT has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GeneWise21_VSMALL_SPECIAL(mat,0,j,AFTER_PY_TRACT); 


      /* Source AFTER_CENTRAL is a special source for AFTER_PY_TRACT */ 
      temp = GeneWise21_VSMALL_SPECIAL(mat,0,j - 1,AFTER_CENTRAL) + (mat->gp21->transition[GP21_CENTRAL2PY]) + (mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)]);    
      if( temp > score ) 
        score = temp;    


      /* Source AFTER_PY_TRACT is a special source for AFTER_PY_TRACT */ 
      temp = GeneWise21_VSMALL_SPECIAL(mat,0,j - 1,AFTER_PY_TRACT) + (mat->gp21->transition[GP21_PY2PY]) + (mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)]);    
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      GeneWise21_VSMALL_SPECIAL(mat,0,j,AFTER_PY_TRACT) = score; 
      /* Finished updating state AFTER_PY_TRACT */ 




      /* Special state AFTER_SPACER has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GeneWise21_VSMALL_SPECIAL(mat,0,j,AFTER_SPACER);   


      /* Source AFTER_PY_TRACT is a special source for AFTER_SPACER */ 
      temp = GeneWise21_VSMALL_SPECIAL(mat,0,j - 1,AFTER_PY_TRACT) + (mat->gp21->transition[GP21_PY2SPACER]) + (mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)]);    
      if( temp > score ) 
        score = temp;    


      /* Source AFTER_SPACER is a special source for AFTER_SPACER */ 
      temp = GeneWise21_VSMALL_SPECIAL(mat,0,j - 1,AFTER_SPACER) + (mat->gp21->transition[GP21_SPACER2SPACER]) + (mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)]);  
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      GeneWise21_VSMALL_SPECIAL(mat,0,j,AFTER_SPACER) = score;   
      /* Finished updating state AFTER_SPACER */ 


      if( bestscore < GeneWise21_VSMALL_SPECIAL(mat,0,j,END) )   
        bestscore = GeneWise21_VSMALL_SPECIAL(mat,0,j,END);  
      } /* end of for all target positions */ 


    mat = free_GeneWise21(mat);  
    return bestscore;    
}    


/* Function:  PackAln_bestmemory_GeneWise21(query,target,gp21,rndcodon,rndbase,dpenv,dpri)
 *
 * Descrip:    This function chooses the best memory set-up for the alignment
 *             using calls to basematrix, and then implements either a large
 *             or small memory model.
 *
 *             It is the best function to use if you just want an alignment
 *
 *             If you want a label alignment, you will need
 *             /convert_PackAln_to_AlnBlock_GeneWise21
 *
 *
 * Arg:           query [UNKN ] query data structure [GeneWiseScore*]
 * Arg:          target [UNKN ] target data structure [ComplexSequence*]
 * Arg:            gp21 [UNKN ] Resource [GeneParser21Score *]
 * Arg:        rndcodon [UNKN ] Resource [RandomCodonScore *]
 * Arg:         rndbase [UNKN ] Resource [RandomModelDNAScore *]
 * Arg:           dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:            dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * PackAln_bestmemory_GeneWise21(GeneWiseScore* query,ComplexSequence* target ,GeneParser21Score * gp21,RandomCodonScore * rndcodon,RandomModelDNAScore * rndbase,DPEnvelope * dpenv,DPRunImpl * dpri) 
{
    long long total; 
    GeneWise21 * mat;    
    PackAln * out;   
    DebugMatrix * de;    
    DPRunImplMemory strategy;    
    assert(dpri);    


    total = query->len * target->seq->len;   
    if( dpri->memory == DPIM_Default )   {  
      if( (total * 21 * sizeof(int)) > 1000*dpri->kbyte_size)    {  
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
        if( (mat=allocate_Expl_GeneWise21(query, target , gp21, rndcodon, rndbase,dpri)) == NULL )   {  
          warn("Unable to allocate large GeneWise21 version");   
          return NULL;   
          }  
        calculate_dpenv_GeneWise21(mat,dpenv);   
        out =  PackAln_read_Expl_GeneWise21(mat);    
        }  
      else   {  
        mat = allocate_GeneWise21_only(query, target , gp21, rndcodon, rndbase);     
        calculate_shatter_GeneWise21(mat,dpenv);     
        out = PackAln_read_Shatter_GeneWise21(mat);  
        }  
      }  
    else {  
      if( strategy == DPIM_Linear )  {  
        /* use small implementation */ 
        if( (mat=allocate_Small_GeneWise21(query, target , gp21, rndcodon, rndbase)) == NULL )   {  
          warn("Unable to allocate small GeneWise21 version");   
          return NULL;   
          }  
        out = PackAln_calculate_Small_GeneWise21(mat,dpenv);     
        }  
      else   {  
        /* use Large implementation */ 
        if( (mat=allocate_Expl_GeneWise21(query, target , gp21, rndcodon, rndbase,dpri)) == NULL )   {  
          warn("Unable to allocate large GeneWise21 version");   
          return NULL;   
          }  
        if( dpri->debug == TRUE) {  
          fatal("Asked for dydebug, but dynamite file not compiled with -g. Need to recompile dynamite source"); 
          }  
        else {  
          calculate_GeneWise21(mat);     
          out =  PackAln_read_Expl_GeneWise21(mat);  
          }  
        }  
      }  


    mat = free_GeneWise21(mat);  
    return out;  
}    


/* Function:  allocate_GeneWise21_only(query,target,gp21,rndcodon,rndbase)
 *
 * Descrip:    This function only allocates the GeneWise21 structure
 *             checks types where possible and determines leni and lenj
 *             The basematrix area is delt with elsewhere
 *
 *
 * Arg:           query [UNKN ] query data structure [GeneWiseScore*]
 * Arg:          target [UNKN ] target data structure [ComplexSequence*]
 * Arg:            gp21 [UNKN ] Resource [GeneParser21Score *]
 * Arg:        rndcodon [UNKN ] Resource [RandomCodonScore *]
 * Arg:         rndbase [UNKN ] Resource [RandomModelDNAScore *]
 *
 * Return [UNKN ]  Undocumented return value [GeneWise21 *]
 *
 */
GeneWise21 * allocate_GeneWise21_only(GeneWiseScore* query,ComplexSequence* target ,GeneParser21Score * gp21,RandomCodonScore * rndcodon,RandomModelDNAScore * rndbase) 
{
    GeneWise21 * out;    


    if((out= GeneWise21_alloc()) == NULL)    {  
      warn("Allocation of basic GeneWise21 structure failed...");    
      return NULL;   
      }  


    out->query = query;  
    out->target = target;    
    out->gp21 = gp21;    
    out->rndcodon = rndcodon;    
    out->rndbase = rndbase;  
    out->leni = query->len;  
    out->lenj = target->seq->len;    
    return out;  
}    


/* Function:  allocate_Expl_GeneWise21(query,target,gp21,rndcodon,rndbase,dpri)
 *
 * Descrip:    This function allocates the GeneWise21 structure
 *             and the basematrix area for explicit memory implementations
 *             It calls /allocate_GeneWise21_only
 *
 *
 * Arg:           query [UNKN ] query data structure [GeneWiseScore*]
 * Arg:          target [UNKN ] target data structure [ComplexSequence*]
 * Arg:            gp21 [UNKN ] Resource [GeneParser21Score *]
 * Arg:        rndcodon [UNKN ] Resource [RandomCodonScore *]
 * Arg:         rndbase [UNKN ] Resource [RandomModelDNAScore *]
 * Arg:            dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [GeneWise21 *]
 *
 */
GeneWise21 * allocate_Expl_GeneWise21(GeneWiseScore* query,ComplexSequence* target ,GeneParser21Score * gp21,RandomCodonScore * rndcodon,RandomModelDNAScore * rndbase,DPRunImpl * dpri) 
{
    GeneWise21 * out;    


    out = allocate_GeneWise21_only(query, target , gp21, rndcodon, rndbase); 
    if( out == NULL )    
      return NULL;   
    if( dpri->should_cache == TRUE ) {  
      if( dpri->cache != NULL )  {  
        if( dpri->cache->maxleni >= (out->lenj+10)*21 && dpri->cache->maxlenj >= (out->leni+1))  
          out->basematrix = hard_link_BaseMatrix(dpri->cache);   
        else 
          dpri->cache = free_BaseMatrix(dpri->cache);    
        }  
      }  
    if( out->basematrix == NULL )    {  
      if( (out->basematrix = BaseMatrix_alloc_matrix_and_specials((out->lenj+10)*21,(out->leni+1),12,out->lenj+10)) == NULL) {  
        warn("Explicit matrix GeneWise21 cannot be allocated, (asking for %d by %d main cells)",out->leni,out->lenj);    
        free_GeneWise21(out);    
        return NULL; 
        }  
      }  
    if( dpri->should_cache == TRUE && dpri->cache == NULL)   
      dpri->cache = hard_link_BaseMatrix(out->basematrix);   
    out->basematrix->type = BASEMATRIX_TYPE_EXPLICIT;    
    init_GeneWise21(out);    
    return out;  
}    


/* Function:  init_GeneWise21(mat)
 *
 * Descrip:    This function initates GeneWise21 matrix when in explicit mode
 *             Called in /allocate_Expl_GeneWise21
 *
 *
 * Arg:        mat [UNKN ] GeneWise21 which contains explicit basematrix memory [GeneWise21 *]
 *
 */
void init_GeneWise21(GeneWise21 * mat) 
{
    register int i;  
    register int j;  
    if( mat->basematrix->type != BASEMATRIX_TYPE_EXPLICIT)   {  
      warn("Cannot iniate matrix, is not an explicit memory type and you have assummed that");   
      return;    
      }  


    for(i= (-1);i<mat->query->len;i++)   {  
      for(j= (-10);j<11;j++) {  
        GeneWise21_EXPL_MATRIX(mat,i,j,MATCH) = NEGI;    
        GeneWise21_EXPL_MATRIX(mat,i,j,INSERT) = NEGI;   
        GeneWise21_EXPL_MATRIX(mat,i,j,DELETE) = NEGI;   
        GeneWise21_EXPL_MATRIX(mat,i,j,CENTRAL_0M) = NEGI;   
        GeneWise21_EXPL_MATRIX(mat,i,j,PY_0M) = NEGI;    
        GeneWise21_EXPL_MATRIX(mat,i,j,SPACER_0M) = NEGI;    
        GeneWise21_EXPL_MATRIX(mat,i,j,CENTRAL_1M) = NEGI;   
        GeneWise21_EXPL_MATRIX(mat,i,j,PY_1M) = NEGI;    
        GeneWise21_EXPL_MATRIX(mat,i,j,SPACER_1M) = NEGI;    
        GeneWise21_EXPL_MATRIX(mat,i,j,CENTRAL_2M) = NEGI;   
        GeneWise21_EXPL_MATRIX(mat,i,j,PY_2M) = NEGI;    
        GeneWise21_EXPL_MATRIX(mat,i,j,SPACER_2M) = NEGI;    
        GeneWise21_EXPL_MATRIX(mat,i,j,CENTRAL_0I) = NEGI;   
        GeneWise21_EXPL_MATRIX(mat,i,j,PY_0I) = NEGI;    
        GeneWise21_EXPL_MATRIX(mat,i,j,SPACER_0I) = NEGI;    
        GeneWise21_EXPL_MATRIX(mat,i,j,CENTRAL_1I) = NEGI;   
        GeneWise21_EXPL_MATRIX(mat,i,j,PY_1I) = NEGI;    
        GeneWise21_EXPL_MATRIX(mat,i,j,SPACER_1I) = NEGI;    
        GeneWise21_EXPL_MATRIX(mat,i,j,CENTRAL_2I) = NEGI;   
        GeneWise21_EXPL_MATRIX(mat,i,j,PY_2I) = NEGI;    
        GeneWise21_EXPL_MATRIX(mat,i,j,SPACER_2I) = NEGI;    
        }  
      }  
    for(j= (-10);j<mat->target->seq->len;j++)    {  
      for(i= (-1);i<2;i++)   {  
        GeneWise21_EXPL_MATRIX(mat,i,j,MATCH) = NEGI;    
        GeneWise21_EXPL_MATRIX(mat,i,j,INSERT) = NEGI;   
        GeneWise21_EXPL_MATRIX(mat,i,j,DELETE) = NEGI;   
        GeneWise21_EXPL_MATRIX(mat,i,j,CENTRAL_0M) = NEGI;   
        GeneWise21_EXPL_MATRIX(mat,i,j,PY_0M) = NEGI;    
        GeneWise21_EXPL_MATRIX(mat,i,j,SPACER_0M) = NEGI;    
        GeneWise21_EXPL_MATRIX(mat,i,j,CENTRAL_1M) = NEGI;   
        GeneWise21_EXPL_MATRIX(mat,i,j,PY_1M) = NEGI;    
        GeneWise21_EXPL_MATRIX(mat,i,j,SPACER_1M) = NEGI;    
        GeneWise21_EXPL_MATRIX(mat,i,j,CENTRAL_2M) = NEGI;   
        GeneWise21_EXPL_MATRIX(mat,i,j,PY_2M) = NEGI;    
        GeneWise21_EXPL_MATRIX(mat,i,j,SPACER_2M) = NEGI;    
        GeneWise21_EXPL_MATRIX(mat,i,j,CENTRAL_0I) = NEGI;   
        GeneWise21_EXPL_MATRIX(mat,i,j,PY_0I) = NEGI;    
        GeneWise21_EXPL_MATRIX(mat,i,j,SPACER_0I) = NEGI;    
        GeneWise21_EXPL_MATRIX(mat,i,j,CENTRAL_1I) = NEGI;   
        GeneWise21_EXPL_MATRIX(mat,i,j,PY_1I) = NEGI;    
        GeneWise21_EXPL_MATRIX(mat,i,j,SPACER_1I) = NEGI;    
        GeneWise21_EXPL_MATRIX(mat,i,j,CENTRAL_2I) = NEGI;   
        GeneWise21_EXPL_MATRIX(mat,i,j,PY_2I) = NEGI;    
        GeneWise21_EXPL_MATRIX(mat,i,j,SPACER_2I) = NEGI;    
        }  
      GeneWise21_EXPL_SPECIAL(mat,i,j,START) = 0;    
      GeneWise21_EXPL_SPECIAL(mat,i,j,END) = NEGI;   
      GeneWise21_EXPL_SPECIAL(mat,i,j,BEFORE_RND) = NEGI;    
      GeneWise21_EXPL_SPECIAL(mat,i,j,BEFORE_CODON) = NEGI;  
      GeneWise21_EXPL_SPECIAL(mat,i,j,BEFORE_CENTRAL) = NEGI;    
      GeneWise21_EXPL_SPECIAL(mat,i,j,BEFORE_PY_TRACT) = NEGI;   
      GeneWise21_EXPL_SPECIAL(mat,i,j,BEFORE_SPACER) = NEGI; 
      GeneWise21_EXPL_SPECIAL(mat,i,j,AFTER_RND) = NEGI; 
      GeneWise21_EXPL_SPECIAL(mat,i,j,AFTER_CODON) = NEGI;   
      GeneWise21_EXPL_SPECIAL(mat,i,j,AFTER_CENTRAL) = NEGI; 
      GeneWise21_EXPL_SPECIAL(mat,i,j,AFTER_PY_TRACT) = NEGI;    
      GeneWise21_EXPL_SPECIAL(mat,i,j,AFTER_SPACER) = NEGI;  
      }  
    return;  
}    


/* Function:  recalculate_PackAln_GeneWise21(pal,mat)
 *
 * Descrip:    This function recalculates the PackAln structure produced by GeneWise21
 *             For example, in linear space methods this is used to score them
 *
 *
 * Arg:        pal [UNKN ] Undocumented argument [PackAln *]
 * Arg:        mat [UNKN ] Undocumented argument [GeneWise21 *]
 *
 */
void recalculate_PackAln_GeneWise21(PackAln * pal,GeneWise21 * mat) 
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
            pau->score = (mat->query->seg[i]->transition[GW_MATCH2MATCH]+mat->query->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]) + (CSEQ_GENOMIC_CDSPOT(mat->target,j));   
            continue;    
            }  
          if( offi == 1 && offj == 3 && prev->state == INSERT )  {  
            pau->score = (mat->query->seg[i]->transition[GW_INSERT2MATCH]+mat->query->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]) + (CSEQ_GENOMIC_CDSPOT(mat->target,j));  
            continue;    
            }  
          if( offi == 1 && offj == 3 && prev->state == DELETE )  {  
            pau->score = (mat->query->seg[i]->transition[GW_DELETE2MATCH]+mat->query->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]) + (CSEQ_GENOMIC_CDSPOT(mat->target,j));  
            continue;    
            }  
          if( offj == 3 && prev->state == (BEFORE_CODON+21) )    {  
            pau->score = (mat->query->seg[i]->transition[GW_START2MATCH]+mat->query->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]) + (CSEQ_GENOMIC_CDSPOT(mat->target,j));   
            continue;    
            }  
          if( offi == 1 && offj == 6 && prev->state == SPACER_0M )   {  
            pau->score = (((mat->gp21->transition[GP21_SPACER2CDS]+mat->query->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)])+CSEQ_GENOMIC_3SS(mat->target,(j-3)))+mat->query->seg[i]->transition[GW_MATCH_BALANCE_3SS]) + (CSEQ_GENOMIC_CDSPOT(mat->target,j));  
            continue;    
            }  
          if( offi == 1 && offj == 6 && prev->state == PY_0M )   {  
            pau->score = (((mat->gp21->transition[GP21_PY2CDS]+mat->query->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)])+CSEQ_GENOMIC_3SS(mat->target,(j-3)))+mat->query->seg[i]->transition[GW_MATCH_BALANCE_3SS]) + (CSEQ_GENOMIC_CDSPOT(mat->target,j));  
            continue;    
            }  
          if( offi == 1 && offj == 5 && prev->state == SPACER_1M )   {  
            pau->score = (mat->gp21->transition[GP21_SPACER2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-2))) + (CSEQ_GENOMIC_CDSPOT(mat->target,j));    
            continue;    
            }  
          if( offi == 1 && offj == 5 && prev->state == PY_1M )   {  
            pau->score = (mat->gp21->transition[GP21_PY2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-2))) + (CSEQ_GENOMIC_CDSPOT(mat->target,j));    
            continue;    
            }  
          if( offi == 1 && offj == 4 && prev->state == SPACER_2M )   {  
            pau->score = (mat->gp21->transition[GP21_SPACER2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-1))) + (CSEQ_GENOMIC_CDSPOT(mat->target,j));    
            continue;    
            }  
          if( offi == 1 && offj == 4 && prev->state == PY_2M )   {  
            pau->score = (mat->gp21->transition[GP21_PY2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-1))) + (CSEQ_GENOMIC_CDSPOT(mat->target,j));    
            continue;    
            }  
          if( offi == 1 && offj == 2 && prev->state == MATCH )   {  
            pau->score = mat->gp21->transition[GP21_DELETE_1_BASE] + (CSEQ_GENOMIC_CDSPOT(mat->target,j));   
            continue;    
            }  
          if( offi == 1 && offj == 1 && prev->state == MATCH )   {  
            pau->score = mat->gp21->transition[GP21_DELETE_2_BASE] + (CSEQ_GENOMIC_CDSPOT(mat->target,j));   
            continue;    
            }  
          if( offi == 1 && offj == 4 && prev->state == MATCH )   {  
            pau->score = mat->gp21->transition[GP21_INSERT_1_BASE] + (CSEQ_GENOMIC_CDSPOT(mat->target,j));   
            continue;    
            }  
          if( offi == 1 && offj == 5 && prev->state == MATCH )   {  
            pau->score = mat->gp21->transition[GP21_INSERT_2_BASE] + (CSEQ_GENOMIC_CDSPOT(mat->target,j));   
            continue;    
            }  
          warn("In recaluclating PackAln with state MATCH, from [%d,%d,%d], got a bad source state. Error!",offi,offj,prev->state);  
          break; 
        case INSERT :    
          if( offi == 0 && offj == 3 && prev->state == MATCH )   {  
            pau->score = (mat->query->seg[i]->transition[GW_MATCH2INSERT]+mat->query->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)]) + (CSEQ_GENOMIC_CDSPOT(mat->target,j));     
            continue;    
            }  
          if( offi == 0 && offj == 3 && prev->state == INSERT )  {  
            pau->score = (mat->query->seg[i]->transition[GW_INSERT2INSERT]+mat->query->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)]) + (CSEQ_GENOMIC_CDSPOT(mat->target,j));    
            continue;    
            }  
          if( offi == 0 && offj == 3 && prev->state == DELETE )  {  
            pau->score = (mat->query->seg[i]->transition[GW_DELETE2INSERT]+mat->query->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)]) + (CSEQ_GENOMIC_CDSPOT(mat->target,j));    
            continue;    
            }  
          if( offj == 3 && prev->state == (BEFORE_CODON+21) )    {  
            pau->score = (mat->query->seg[i]->transition[GW_START2INSERT]+mat->query->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)]) + (CSEQ_GENOMIC_CDSPOT(mat->target,j));     
            continue;    
            }  
          if( offi == 1 && offj == 6 && prev->state == SPACER_0I )   {  
            pau->score = (((mat->gp21->transition[GP21_SPACER2CDS]+mat->query->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)])+CSEQ_GENOMIC_3SS(mat->target,(j-3)))+mat->query->seg[i]->transition[GW_MATCH_BALANCE_3SS]) + (CSEQ_GENOMIC_CDSPOT(mat->target,j));     
            continue;    
            }  
          if( offi == 1 && offj == 5 && prev->state == SPACER_1I )   {  
            pau->score = (mat->gp21->transition[GP21_SPACER2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-2))) + (CSEQ_GENOMIC_CDSPOT(mat->target,j));    
            continue;    
            }  
          if( offi == 1 && offj == 4 && prev->state == SPACER_2I )   {  
            pau->score = (mat->gp21->transition[GP21_SPACER2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-1))) + (CSEQ_GENOMIC_CDSPOT(mat->target,j));    
            continue;    
            }  
          if( offi == 1 && offj == 6 && prev->state == PY_0I )   {  
            pau->score = (((mat->gp21->transition[GP21_PY2CDS]+mat->query->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)])+CSEQ_GENOMIC_3SS(mat->target,(j-3)))+mat->query->seg[i]->transition[GW_MATCH_BALANCE_3SS]) + (CSEQ_GENOMIC_CDSPOT(mat->target,j));     
            continue;    
            }  
          if( offi == 1 && offj == 5 && prev->state == PY_1I )   {  
            pau->score = (mat->gp21->transition[GP21_PY2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-2))) + (CSEQ_GENOMIC_CDSPOT(mat->target,j));    
            continue;    
            }  
          if( offi == 1 && offj == 4 && prev->state == PY_2I )   {  
            pau->score = (mat->gp21->transition[GP21_PY2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-1))) + (CSEQ_GENOMIC_CDSPOT(mat->target,j));    
            continue;    
            }  
          if( offi == 1 && offj == 2 && prev->state == INSERT )  {  
            pau->score = mat->gp21->transition[GP21_DELETE_1_BASE] + (CSEQ_GENOMIC_CDSPOT(mat->target,j));   
            continue;    
            }  
          if( offi == 1 && offj == 1 && prev->state == INSERT )  {  
            pau->score = mat->gp21->transition[GP21_DELETE_2_BASE] + (CSEQ_GENOMIC_CDSPOT(mat->target,j));   
            continue;    
            }  
          if( offi == 1 && offj == 4 && prev->state == INSERT )  {  
            pau->score = mat->gp21->transition[GP21_INSERT_1_BASE] + (CSEQ_GENOMIC_CDSPOT(mat->target,j));   
            continue;    
            }  
          if( offi == 1 && offj == 5 && prev->state == INSERT )  {  
            pau->score = mat->gp21->transition[GP21_INSERT_2_BASE] + (CSEQ_GENOMIC_CDSPOT(mat->target,j));   
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
          if( offj == 0 && prev->state == (BEFORE_CODON+21) )    {  
            pau->score = mat->query->seg[i]->transition[GW_START2DELETE] + (0);  
            continue;    
            }  
          warn("In recaluclating PackAln with state DELETE, from [%d,%d,%d], got a bad source state. Error!",offi,offj,prev->state); 
          break; 
        case CENTRAL_0M :    
          if( offi == 0 && offj == 8 && prev->state == MATCH )   {  
            pau->score = (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_MATCH2MATCH])+mat->query->seg[i]->transition[GW_MATCH_BALANCE_5SS]) + (0);   
            continue;    
            }  
          if( offi == 0 && offj == 8 && prev->state == INSERT )  {  
            pau->score = (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_INSERT2MATCH])+mat->query->seg[i]->transition[GW_INSERT_BALANCE_5SS]) + (0);     
            continue;    
            }  
          if( offi == 0 && offj == 8 && prev->state == DELETE )  {  
            pau->score = ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_DELETE2MATCH]) + (0);     
            continue;    
            }  
          if( offi == 0 && offj == 1 && prev->state == CENTRAL_0M )  {  
            pau->score = (mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CENTRAL2CENTRAL]) + (0);   
            continue;    
            }  
          warn("In recaluclating PackAln with state CENTRAL_0M, from [%d,%d,%d], got a bad source state. Error!",offi,offj,prev->state); 
          break; 
        case PY_0M :     
          if( offi == 0 && offj == 1 && prev->state == CENTRAL_0M )  {  
            pau->score = mat->gp21->transition[GP21_CENTRAL2PY] + (mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)]);     
            continue;    
            }  
          if( offi == 0 && offj == 1 && prev->state == PY_0M )   {  
            pau->score = mat->gp21->transition[GP21_PY2PY] + (mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)]);  
            continue;    
            }  
          warn("In recaluclating PackAln with state PY_0M, from [%d,%d,%d], got a bad source state. Error!",offi,offj,prev->state);  
          break; 
        case SPACER_0M :     
          if( offi == 0 && offj == 1 && prev->state == PY_0M )   {  
            pau->score = mat->gp21->transition[GP21_PY2SPACER] + (mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)]);  
            continue;    
            }  
          if( offi == 0 && offj == 1 && prev->state == SPACER_0M )   {  
            pau->score = mat->gp21->transition[GP21_SPACER2SPACER] + (mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)]);  
            continue;    
            }  
          warn("In recaluclating PackAln with state SPACER_0M, from [%d,%d,%d], got a bad source state. Error!",offi,offj,prev->state);  
          break; 
        case CENTRAL_1M :    
          if( offi == 0 && offj == 9 && prev->state == MATCH )   {  
            pau->score = (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_MATCH2MATCH])+mat->query->seg[i]->transition[GW_MATCH_BALANCE_5SS]) + (0);   
            continue;    
            }  
          if( offi == 0 && offj == 9 && prev->state == INSERT )  {  
            pau->score = (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_INSERT2MATCH])+mat->query->seg[i]->transition[GW_MATCH_BALANCE_5SS]) + (0);  
            continue;    
            }  
          if( offi == 0 && offj == 9 && prev->state == DELETE )  {  
            pau->score = ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_DELETE2MATCH]) + (0);     
            continue;    
            }  
          if( offi == 0 && offj == 1 && prev->state == CENTRAL_1M )  {  
            pau->score = (mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CENTRAL2CENTRAL]) + (0);   
            continue;    
            }  
          warn("In recaluclating PackAln with state CENTRAL_1M, from [%d,%d,%d], got a bad source state. Error!",offi,offj,prev->state); 
          break; 
        case PY_1M :     
          if( offi == 0 && offj == 1 && prev->state == CENTRAL_1M )  {  
            pau->score = mat->gp21->transition[GP21_CENTRAL2PY] + (mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)]);     
            continue;    
            }  
          if( offi == 0 && offj == 1 && prev->state == PY_1M )   {  
            pau->score = mat->gp21->transition[GP21_PY2PY] + (mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)]);  
            continue;    
            }  
          warn("In recaluclating PackAln with state PY_1M, from [%d,%d,%d], got a bad source state. Error!",offi,offj,prev->state);  
          break; 
        case SPACER_1M :     
          if( offi == 0 && offj == 1 && prev->state == PY_1M )   {  
            pau->score = mat->gp21->transition[GP21_PY2SPACER] + (mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)]);  
            continue;    
            }  
          if( offi == 0 && offj == 1 && prev->state == SPACER_1M )   {  
            pau->score = mat->gp21->transition[GP21_SPACER2SPACER] + (mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)]);  
            continue;    
            }  
          warn("In recaluclating PackAln with state SPACER_1M, from [%d,%d,%d], got a bad source state. Error!",offi,offj,prev->state);  
          break; 
        case CENTRAL_2M :    
          if( offi == 0 && offj == 10 && prev->state == MATCH )  {  
            pau->score = (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_MATCH2MATCH])+mat->query->seg[i]->transition[GW_MATCH_BALANCE_5SS]) + (0);   
            continue;    
            }  
          if( offi == 0 && offj == 10 && prev->state == INSERT ) {  
            pau->score = (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_INSERT2MATCH])+mat->query->seg[i]->transition[GW_MATCH_BALANCE_5SS]) + (0);  
            continue;    
            }  
          if( offi == 0 && offj == 10 && prev->state == DELETE ) {  
            pau->score = ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_DELETE2MATCH]) + (0);     
            continue;    
            }  
          if( offi == 0 && offj == 1 && prev->state == CENTRAL_2M )  {  
            pau->score = (mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CENTRAL2CENTRAL]) + (0);   
            continue;    
            }  
          warn("In recaluclating PackAln with state CENTRAL_2M, from [%d,%d,%d], got a bad source state. Error!",offi,offj,prev->state); 
          break; 
        case PY_2M :     
          if( offi == 0 && offj == 1 && prev->state == CENTRAL_2M )  {  
            pau->score = mat->gp21->transition[GP21_CENTRAL2PY] + (mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)]);     
            continue;    
            }  
          if( offi == 0 && offj == 1 && prev->state == PY_2M )   {  
            pau->score = mat->gp21->transition[GP21_PY2PY] + (mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)]);  
            continue;    
            }  
          warn("In recaluclating PackAln with state PY_2M, from [%d,%d,%d], got a bad source state. Error!",offi,offj,prev->state);  
          break; 
        case SPACER_2M :     
          if( offi == 0 && offj == 1 && prev->state == PY_2M )   {  
            pau->score = mat->gp21->transition[GP21_PY2SPACER] + (mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)]);  
            continue;    
            }  
          if( offi == 0 && offj == 1 && prev->state == SPACER_2M )   {  
            pau->score = mat->gp21->transition[GP21_SPACER2SPACER] + (mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)]);  
            continue;    
            }  
          warn("In recaluclating PackAln with state SPACER_2M, from [%d,%d,%d], got a bad source state. Error!",offi,offj,prev->state);  
          break; 
        case CENTRAL_0I :    
          if( offi == 0 && offj == 8 && prev->state == MATCH )   {  
            pau->score = (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_MATCH2INSERT])+mat->query->seg[i]->transition[GW_MATCH_BALANCE_5SS]) + (0);  
            continue;    
            }  
          if( offi == 0 && offj == 8 && prev->state == INSERT )  {  
            pau->score = (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_INSERT2INSERT])+mat->query->seg[i]->transition[GW_INSERT_BALANCE_5SS]) + (0);    
            continue;    
            }  
          if( offi == 0 && offj == 8 && prev->state == DELETE )  {  
            pau->score = ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_DELETE2INSERT]) + (0);    
            continue;    
            }  
          if( offi == 0 && offj == 1 && prev->state == CENTRAL_0I )  {  
            pau->score = (mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CENTRAL2CENTRAL]) + (0);   
            continue;    
            }  
          warn("In recaluclating PackAln with state CENTRAL_0I, from [%d,%d,%d], got a bad source state. Error!",offi,offj,prev->state); 
          break; 
        case PY_0I :     
          if( offi == 0 && offj == 1 && prev->state == CENTRAL_0I )  {  
            pau->score = mat->gp21->transition[GP21_CENTRAL2PY] + (mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)]);     
            continue;    
            }  
          if( offi == 0 && offj == 1 && prev->state == PY_0I )   {  
            pau->score = mat->gp21->transition[GP21_PY2PY] + (mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)]);  
            continue;    
            }  
          warn("In recaluclating PackAln with state PY_0I, from [%d,%d,%d], got a bad source state. Error!",offi,offj,prev->state);  
          break; 
        case SPACER_0I :     
          if( offi == 0 && offj == 1 && prev->state == PY_0I )   {  
            pau->score = mat->gp21->transition[GP21_PY2SPACER] + (mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)]);  
            continue;    
            }  
          if( offi == 0 && offj == 1 && prev->state == SPACER_0I )   {  
            pau->score = mat->gp21->transition[GP21_SPACER2SPACER] + (mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)]);  
            continue;    
            }  
          warn("In recaluclating PackAln with state SPACER_0I, from [%d,%d,%d], got a bad source state. Error!",offi,offj,prev->state);  
          break; 
        case CENTRAL_1I :    
          if( offi == 0 && offj == 9 && prev->state == MATCH )   {  
            pau->score = ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_MATCH2INSERT]) + (0);     
            continue;    
            }  
          if( offi == 0 && offj == 9 && prev->state == INSERT )  {  
            pau->score = ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_INSERT2INSERT]) + (0);    
            continue;    
            }  
          if( offi == 0 && offj == 9 && prev->state == DELETE )  {  
            pau->score = ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_DELETE2INSERT]) + (0);    
            continue;    
            }  
          if( offi == 0 && offj == 1 && prev->state == CENTRAL_1I )  {  
            pau->score = (mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CENTRAL2CENTRAL]) + (0);   
            continue;    
            }  
          warn("In recaluclating PackAln with state CENTRAL_1I, from [%d,%d,%d], got a bad source state. Error!",offi,offj,prev->state); 
          break; 
        case PY_1I :     
          if( offi == 0 && offj == 1 && prev->state == CENTRAL_1I )  {  
            pau->score = mat->gp21->transition[GP21_CENTRAL2PY] + (mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)]);     
            continue;    
            }  
          if( offi == 0 && offj == 1 && prev->state == PY_1I )   {  
            pau->score = mat->gp21->transition[GP21_PY2PY] + (mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)]);  
            continue;    
            }  
          warn("In recaluclating PackAln with state PY_1I, from [%d,%d,%d], got a bad source state. Error!",offi,offj,prev->state);  
          break; 
        case SPACER_1I :     
          if( offi == 0 && offj == 1 && prev->state == PY_1I )   {  
            pau->score = mat->gp21->transition[GP21_PY2SPACER] + (mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)]);  
            continue;    
            }  
          if( offi == 0 && offj == 1 && prev->state == SPACER_1I )   {  
            pau->score = mat->gp21->transition[GP21_SPACER2SPACER] + (mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)]);  
            continue;    
            }  
          warn("In recaluclating PackAln with state SPACER_1I, from [%d,%d,%d], got a bad source state. Error!",offi,offj,prev->state);  
          break; 
        case CENTRAL_2I :    
          if( offi == 0 && offj == 10 && prev->state == MATCH )  {  
            pau->score = ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_MATCH2INSERT]) + (0);     
            continue;    
            }  
          if( offi == 0 && offj == 10 && prev->state == INSERT ) {  
            pau->score = ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_INSERT2INSERT]) + (0);    
            continue;    
            }  
          if( offi == 0 && offj == 10 && prev->state == INSERT ) {  
            pau->score = ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_DELETE2INSERT]) + (0);    
            continue;    
            }  
          if( offi == 0 && offj == 1 && prev->state == CENTRAL_2I )  {  
            pau->score = (mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CENTRAL2CENTRAL]) + (0);   
            continue;    
            }  
          warn("In recaluclating PackAln with state CENTRAL_2I, from [%d,%d,%d], got a bad source state. Error!",offi,offj,prev->state); 
          break; 
        case PY_2I :     
          if( offi == 0 && offj == 1 && prev->state == CENTRAL_2I )  {  
            pau->score = mat->gp21->transition[GP21_CENTRAL2PY] + (mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)]);     
            continue;    
            }  
          if( offi == 0 && offj == 1 && prev->state == PY_2I )   {  
            pau->score = mat->gp21->transition[GP21_PY2PY] + (mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)]);  
            continue;    
            }  
          warn("In recaluclating PackAln with state PY_2I, from [%d,%d,%d], got a bad source state. Error!",offi,offj,prev->state);  
          break; 
        case SPACER_2I :     
          if( offi == 0 && offj == 1 && prev->state == PY_2I )   {  
            pau->score = mat->gp21->transition[GP21_PY2SPACER] + (mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)]);  
            continue;    
            }  
          if( offi == 0 && offj == 1 && prev->state == SPACER_2I )   {  
            pau->score = mat->gp21->transition[GP21_SPACER2SPACER] + (mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)]);  
            continue;    
            }  
          warn("In recaluclating PackAln with state SPACER_2I, from [%d,%d,%d], got a bad source state. Error!",offi,offj,prev->state);  
          break; 
        case (START+21) :    
          warn("In recaluclating PackAln with state START, got a bad source state. Error!"); 
          break; 
        case (END+21) :  
          if( offj == 1 && prev->state == (AFTER_RND+21) )   {  
            pau->score = 0 + (0);    
            continue;    
            }  
          warn("In recaluclating PackAln with state END, got a bad source state. Error!");   
          break; 
        case (BEFORE_RND+21) :   
          if( offj == 1 && prev->state == (START+21) )   {  
            pau->score = mat->rndbase->base[CSEQ_GENOMIC_BASE(mat->target,j)] + (0);     
            continue;    
            }  
          if( offj == 1 && prev->state == (BEFORE_RND+21) )  {  
            pau->score = (mat->rndbase->base[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_RND2RND]) + (0);   
            continue;    
            }  
          if( offj == 1 && prev->state == (BEFORE_CODON+21) )    {  
            pau->score = (mat->rndbase->base[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CDS2RND]) + (0);   
            continue;    
            }  
          warn("In recaluclating PackAln with state BEFORE_RND, got a bad source state. Error!");    
          break; 
        case (BEFORE_CODON+21) :     
          if( offj == 3 && prev->state == (BEFORE_RND+21) )  {  
            pau->score = (mat->rndcodon->codon[CSEQ_GENOMIC_CODON(mat->target,j)]+mat->gp21->transition[GP21_RND2CDS]) + (0);    
            continue;    
            }  
          if( offj == 3 && prev->state == (BEFORE_CODON+21) )    {  
            pau->score = (mat->rndcodon->codon[CSEQ_GENOMIC_CODON(mat->target,j)]+mat->gp21->transition[GP21_CDS2CDS]) + (0);    
            continue;    
            }  
          if( offj == 5 && prev->state == (BEFORE_SPACER+21) )   {  
            pau->score = ((mat->gp21->transition[GP21_SPACER2CDS]+mat->rndcodon->codon[CSEQ_GENOMIC_CODON(mat->target,j)])+CSEQ_GENOMIC_3SS(mat->target,(j-3))) + (0);   
            continue;    
            }  
          if( offj == 4 && prev->state == (BEFORE_SPACER+21) )   {  
            pau->score = (mat->gp21->transition[GP21_SPACER2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-2))) + (0);     
            continue;    
            }  
          if( offj == 3 && prev->state == (BEFORE_SPACER+21) )   {  
            pau->score = (mat->gp21->transition[GP21_SPACER2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-1))) + (0);     
            continue;    
            }  
          warn("In recaluclating PackAln with state BEFORE_CODON, got a bad source state. Error!");  
          break; 
        case (BEFORE_CENTRAL+21) :   
          if( offj == 1 && prev->state == (BEFORE_CENTRAL+21) )  {  
            pau->score = (mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CENTRAL2CENTRAL]) + (0);   
            continue;    
            }  
          if( offj == 8 && prev->state == (BEFORE_CODON+21) )    {  
            pau->score = (mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7))) + (0);   
            continue;    
            }  
          if( offj == 9 && prev->state == (BEFORE_CODON+21) )    {  
            pau->score = (mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7))) + (0);   
            continue;    
            }  
          if( offj == 10 && prev->state == (BEFORE_CODON+21) )   {  
            pau->score = (mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7))) + (0);   
            continue;    
            }  
          warn("In recaluclating PackAln with state BEFORE_CENTRAL, got a bad source state. Error!");    
          break; 
        case (BEFORE_PY_TRACT+21) :  
          if( offj == 1 && prev->state == (BEFORE_CENTRAL+21) )  {  
            pau->score = mat->gp21->transition[GP21_CENTRAL2PY] + (mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)]);     
            continue;    
            }  
          if( offj == 1 && prev->state == (BEFORE_PY_TRACT+21) ) {  
            pau->score = mat->gp21->transition[GP21_PY2PY] + (mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)]);  
            continue;    
            }  
          warn("In recaluclating PackAln with state BEFORE_PY_TRACT, got a bad source state. Error!");   
          break; 
        case (BEFORE_SPACER+21) :    
          if( offj == 1 && prev->state == (BEFORE_PY_TRACT+21) ) {  
            pau->score = mat->gp21->transition[GP21_PY2SPACER] + (mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)]);  
            continue;    
            }  
          if( offj == 1 && prev->state == (BEFORE_SPACER+21) )   {  
            pau->score = mat->gp21->transition[GP21_SPACER2SPACER] + (mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)]);  
            continue;    
            }  
          warn("In recaluclating PackAln with state BEFORE_SPACER, got a bad source state. Error!"); 
          break; 
        case (AFTER_RND+21) :    
          if( offj == 1 && prev->state == (AFTER_RND+21) )   {  
            pau->score = (mat->rndbase->base[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_RND2RND]) + (0);   
            continue;    
            }  
          if( offj == 1 && prev->state == (AFTER_CODON+21) ) {  
            pau->score = (mat->rndbase->base[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CDS2RND]) + (0);   
            continue;    
            }  
          warn("In recaluclating PackAln with state AFTER_RND, got a bad source state. Error!"); 
          break; 
        case (AFTER_CODON+21) :  
          if( offj == 3 && prev->state == (AFTER_RND+21) )   {  
            pau->score = (mat->rndcodon->codon[CSEQ_GENOMIC_CODON(mat->target,j)]+mat->gp21->transition[GP21_RND2CDS]) + (0);    
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
          if( offj == 3 && prev->state == (AFTER_CODON+21) ) {  
            pau->score = (mat->rndcodon->codon[CSEQ_GENOMIC_CODON(mat->target,j)]+mat->gp21->transition[GP21_CDS2CDS]) + (0);    
            continue;    
            }  
          if( offj == 5 && prev->state == (AFTER_SPACER+21) )    {  
            pau->score = ((mat->gp21->transition[GP21_SPACER2CDS]+mat->rndcodon->codon[CSEQ_GENOMIC_CODON(mat->target,j)])+CSEQ_GENOMIC_3SS(mat->target,(j-3))) + (0);   
            continue;    
            }  
          if( offj == 4 && prev->state == (AFTER_SPACER+21) )    {  
            pau->score = (mat->gp21->transition[GP21_SPACER2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-2))) + (0);     
            continue;    
            }  
          if( offj == 3 && prev->state == (AFTER_SPACER+21) )    {  
            pau->score = (mat->gp21->transition[GP21_SPACER2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-1))) + (0);     
            continue;    
            }  
          warn("In recaluclating PackAln with state AFTER_CODON, got a bad source state. Error!");   
          break; 
        case (AFTER_CENTRAL+21) :    
          if( offj == 1 && prev->state == (AFTER_CENTRAL+21) )   {  
            pau->score = (mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CENTRAL2CENTRAL]) + (0);   
            continue;    
            }  
          if( offj == 8 && prev->state == (AFTER_CODON+21) ) {  
            pau->score = (mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7))) + (0);   
            continue;    
            }  
          if( offj == 9 && prev->state == (AFTER_CODON+21) ) {  
            pau->score = (mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7))) + (0);   
            continue;    
            }  
          if( offj == 10 && prev->state == (AFTER_CODON+21) )    {  
            pau->score = (mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7))) + (0);   
            continue;    
            }  
          warn("In recaluclating PackAln with state AFTER_CENTRAL, got a bad source state. Error!"); 
          break; 
        case (AFTER_PY_TRACT+21) :   
          if( offj == 1 && prev->state == (AFTER_CENTRAL+21) )   {  
            pau->score = mat->gp21->transition[GP21_CENTRAL2PY] + (mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)]);     
            continue;    
            }  
          if( offj == 1 && prev->state == (AFTER_PY_TRACT+21) )  {  
            pau->score = mat->gp21->transition[GP21_PY2PY] + (mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)]);  
            continue;    
            }  
          warn("In recaluclating PackAln with state AFTER_PY_TRACT, got a bad source state. Error!");    
          break; 
        case (AFTER_SPACER+21) :     
          if( offj == 1 && prev->state == (AFTER_PY_TRACT+21) )  {  
            pau->score = mat->gp21->transition[GP21_PY2SPACER] + (mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)]);  
            continue;    
            }  
          if( offj == 1 && prev->state == (AFTER_SPACER+21) )    {  
            pau->score = mat->gp21->transition[GP21_SPACER2SPACER] + (mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)]);  
            continue;    
            }  
          warn("In recaluclating PackAln with state AFTER_SPACER, got a bad source state. Error!");  
          break; 
        default :    
          warn("In recaluclating PackAln got a bad recipient state. Error!");    
        }  
      prev = pau;    
      }  
    return;  
}    
/* divide and conquor macros are next */ 
#define GeneWise21_HIDDEN_MATRIX(thismatrix,i,j,state) (thismatrix->basematrix->matrix[(j-hiddenj+10)][(i+1)*21+state])  
#define GeneWise21_DC_SHADOW_MATRIX(thismatrix,i,j,state) (thismatrix->basematrix->matrix[((j+11)*8) % 88][(i+1)*21+state])  
#define GeneWise21_HIDDEN_SPECIAL(thismatrix,i,j,state) (thismatrix->basematrix->specmatrix[state][(j+10)])  
#define GeneWise21_DC_SHADOW_SPECIAL(thismatrix,i,j,state) (thismatrix->basematrix->specmatrix[state*8][(j+10)]) 
#define GeneWise21_DC_SHADOW_MATRIX_SP(thismatrix,i,j,state,shadow) (thismatrix->basematrix->matrix[((((j+11)*8)+(shadow+1)) % 88)][(i+1)*21 + state])   
#define GeneWise21_DC_SHADOW_SPECIAL_SP(thismatrix,i,j,state,shadow) (thismatrix->basematrix->specmatrix[state*8 +shadow+1][(j+10)]) 
#define GeneWise21_DC_OPT_SHADOW_MATRIX(thismatrix,i,j,state) (score_pointers[(((j+10)% 10) * (leni+1) * 21) + ((i+1) * 21) + (state)])  
#define GeneWise21_DC_OPT_SHADOW_MATRIX_SP(thismatrix,i,j,state,shadow) (shadow_pointers[(((j+10)% 10) * (leni+1) * 168) + ((i+1) * 168) + (state * 8) + shadow+1])  
#define GeneWise21_DC_OPT_SHADOW_SPECIAL(thismatrix,i,j,state) (thismatrix->basematrix->specmatrix[state*8][(j+10)]) 
/* Function:  allocate_Small_GeneWise21(query,target,gp21,rndcodon,rndbase)
 *
 * Descrip:    This function allocates the GeneWise21 structure
 *             and the basematrix area for a small memory implementations
 *             It calls /allocate_GeneWise21_only
 *
 *
 * Arg:           query [UNKN ] query data structure [GeneWiseScore*]
 * Arg:          target [UNKN ] target data structure [ComplexSequence*]
 * Arg:            gp21 [UNKN ] Resource [GeneParser21Score *]
 * Arg:        rndcodon [UNKN ] Resource [RandomCodonScore *]
 * Arg:         rndbase [UNKN ] Resource [RandomModelDNAScore *]
 *
 * Return [UNKN ]  Undocumented return value [GeneWise21 *]
 *
 */
#define GeneWise21_DC_OPT_SHADOW_SPECIAL_SP(thismatrix,i,j,state,shadow) (thismatrix->basematrix->specmatrix[state*8 +shadow+1][(j+10)]) 
GeneWise21 * allocate_Small_GeneWise21(GeneWiseScore* query,ComplexSequence* target ,GeneParser21Score * gp21,RandomCodonScore * rndcodon,RandomModelDNAScore * rndbase) 
{
    GeneWise21 * out;    


    out = allocate_GeneWise21_only(query, target , gp21, rndcodon, rndbase); 
    if( out == NULL )    
      return NULL;   
    out->basematrix = BaseMatrix_alloc_matrix_and_specials(88,(out->leni + 1) * 21,96,out->lenj+10);     
    if(out == NULL)  {  
      warn("Small shadow matrix GeneWise21 cannot be allocated, (asking for 11 by %d main cells)",out->leni+2);  
      free_GeneWise21(out);  
      return NULL;   
      }  
    out->basematrix->type = BASEMATRIX_TYPE_SHADOW;  
    return out;  
}    


/* Function:  PackAln_calculate_Small_GeneWise21(mat,dpenv)
 *
 * Descrip:    This function calculates an alignment for GeneWise21 structure in linear space
 *             If you want only the start/end points
 *             use /AlnRangeSet_calculate_Small_GeneWise21 
 *
 *             The function basically
 *               finds start/end points 
 *               foreach start/end point 
 *                 calls /full_dc_GeneWise21 
 *
 *
 * Arg:          mat [UNKN ] Undocumented argument [GeneWise21 *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * PackAln_calculate_Small_GeneWise21(GeneWise21 * mat,DPEnvelope * dpenv) 
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
      warn("Could not calculate packaln small for GeneWise21 due to wrong type of matrix");  
      return NULL;   
      }  


    out = PackAln_alloc_std();   


    start_reporting("Find start end points: ");  
    dc_optimised_start_end_calc_GeneWise21(mat,dpenv);   
    score = start_end_find_end_GeneWise21(mat,&endj);    
    out->score = score;  
    stopstate = END;
    
    /* Special to specials: have to eat up in strip and then drop back to full_dc for intervening bits */ 
    log_full_error(REPORT,0,"End at %d Score %d",endj,score);    
    stop_reporting();    
    for(;;)  { /*while there are more special bits to recover*/ 
      start_reporting("Special cell aln end   %d:",endj);    
      if( read_special_strip_GeneWise21(mat,0,endj,stopstate,&endj,&startstate,out) == FALSE )   {  
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
      starti = GeneWise21_DC_SHADOW_SPECIAL_SP(mat,0,endj,temp,0);   
      startj = GeneWise21_DC_SHADOW_SPECIAL_SP(mat,0,endj,temp,1);   
      startstate = GeneWise21_DC_SHADOW_SPECIAL_SP(mat,0,endj,temp,2);   
      stopi = GeneWise21_DC_SHADOW_SPECIAL_SP(mat,0,endj,temp,3);    
      stopj = GeneWise21_DC_SHADOW_SPECIAL_SP(mat,0,endj,temp,4);    
      stopstate = GeneWise21_DC_SHADOW_SPECIAL_SP(mat,0,endj,temp,5);    


      /* Get out the score of this block. V. important! */ 
      temp = GeneWise21_DC_SHADOW_SPECIAL_SP(mat,0,endj,temp,6); 
      totalj = stopj - startj;   
      donej  = 0;    
      start_reporting("Main matrix  aln [%d,%d]:",startj,stopj);     
      if(full_dc_GeneWise21(mat,starti,startj,startstate,stopi,stopj,stopstate,out,&donej,totalj,dpenv) == FALSE)    {  
        warn("In the alignment GeneWise21 [%d,%d][%d,%d], got a problem. Please report bug ... giving you back a partial alignment",starti,startj,stopi,stopj);  
        return out;  
        }  


      /* now have to figure out which special we came from... yikes */ 
      max_matrix_to_special_GeneWise21(mat,starti,startj,startstate,temp,&stopi,&stopj,&stopstate,&temp,NULL);   
      if( stopi == GeneWise21_READ_OFF_ERROR)    {  
        warn("In GeneWise21 read off ending at %d ... got a bad matrix to special read off... returning partial alignment",startj);  
        invert_PackAln(out); 
        recalculate_PackAln_GeneWise21(out,mat); 
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
    recalculate_PackAln_GeneWise21(out,mat); 
    return out;  


}    


/* Function:  AlnRangeSet_calculate_Small_GeneWise21(mat)
 *
 * Descrip:    This function calculates an alignment for GeneWise21 structure in linear space
 *             If you want the full alignment, use /PackAln_calculate_Small_GeneWise21 
 *             If you have already got the full alignment, but want the range set, use /AlnRangeSet_from_PackAln_GeneWise21
 *             If you have got the small matrix but not the alignment, use /AlnRangeSet_from_GeneWise21 
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [GeneWise21 *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * AlnRangeSet_calculate_Small_GeneWise21(GeneWise21 * mat) 
{
    AlnRangeSet * out;   


    start_reporting("Find start end points: ");  
    dc_optimised_start_end_calc_GeneWise21(mat,NULL);    
    log_full_error(REPORT,0,"Calculated");   


    out = AlnRangeSet_from_GeneWise21(mat);  
    return out;  
}    


/* Function:  AlnRangeSet_from_GeneWise21(mat)
 *
 * Descrip:    This function reads off a start/end structure
 *             for GeneWise21 structure in linear space
 *             If you want the full alignment use
 *             /PackAln_calculate_Small_GeneWise21 
 *             If you have not calculated the matrix use
 *             /AlnRange_calculate_Small_GeneWise21
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [GeneWise21 *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * AlnRangeSet_from_GeneWise21(GeneWise21 * mat) 
{
    AlnRangeSet * out;   
    AlnRange * temp; 
    int jpos;    
    int state;   


    if( mat->basematrix->type != BASEMATRIX_TYPE_SHADOW) {  
      warn("Bad error! - non shadow matrix type in AlnRangeSet_from_GeneWise21");    
      return NULL;   
      }  


    out = AlnRangeSet_alloc_std();   
    /* Find the end position */ 
    out->score = start_end_find_end_GeneWise21(mat,&jpos);   
    state = END; 


    while( (temp = AlnRange_build_GeneWise21(mat,jpos,state,&jpos,&state)) != NULL)  
      add_AlnRangeSet(out,temp); 
    return out;  
}    


/* Function:  AlnRange_build_GeneWise21(mat,stopj,stopspecstate,startj,startspecstate)
 *
 * Descrip:    This function calculates a single start/end set in linear space
 *             Really a sub-routine for /AlnRangeSet_from_PackAln_GeneWise21
 *
 *
 * Arg:                   mat [UNKN ] Undocumented argument [GeneWise21 *]
 * Arg:                 stopj [UNKN ] Undocumented argument [int]
 * Arg:         stopspecstate [UNKN ] Undocumented argument [int]
 * Arg:                startj [UNKN ] Undocumented argument [int *]
 * Arg:        startspecstate [UNKN ] Undocumented argument [int *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRange *]
 *
 */
AlnRange * AlnRange_build_GeneWise21(GeneWise21 * mat,int stopj,int stopspecstate,int * startj,int * startspecstate) 
{
    AlnRange * out;  
    int jpos;    
    int state;   


    if( mat->basematrix->type != BASEMATRIX_TYPE_SHADOW) {  
      warn("Bad error! - non shadow matrix type in AlnRangeSet_from_GeneWise21");    
      return NULL;   
      }  


    /* Assumme that we have specials (we should!). Read back along the specials till we have the finish point */ 
    if( read_special_strip_GeneWise21(mat,0,stopj,stopspecstate,&jpos,&state,NULL) == FALSE) {  
      warn("In AlnRanger_build_GeneWise21 alignment ending at %d, unable to read back specials. Will (evenutally) return a partial range set... BEWARE!",stopj); 
      return NULL;   
      }  
    if( state == START || jpos <= 0) 
      return NULL;   


    out = AlnRange_alloc();  


    out->starti = GeneWise21_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,0);   
    out->startj = GeneWise21_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,1);   
    out->startstate = GeneWise21_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,2);   
    out->stopi = GeneWise21_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,3);    
    out->stopj = GeneWise21_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,4);    
    out->stopstate = GeneWise21_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,5);    
    out->startscore = GeneWise21_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,6);   
    out->stopscore = GeneWise21_DC_SHADOW_SPECIAL(mat,0,jpos,state); 


    /* Now, we have to figure out where this state came from in the specials */ 
    max_matrix_to_special_GeneWise21(mat,out->starti,out->startj,out->startstate,out->startscore,&jpos,startj,startspecstate,&state,NULL);   
    if( jpos == GeneWise21_READ_OFF_ERROR)   {  
      warn("In AlnRange_build_GeneWise21 alignment ending at %d, with aln range between %d-%d in j, unable to find source special, returning this range, but this could get tricky!",stopj,out->startj,out->stopj);  
      return out;    
      }  


    /* Put in the correct score for startstate, from the special */ 
    out->startscore = GeneWise21_DC_SHADOW_SPECIAL(mat,0,*startj,*startspecstate);   
    /* The correct j coords have been put into startj, startspecstate... so just return out */ 
    return out;  
}    


/* Function:  read_hidden_GeneWise21(mat,starti,startj,startstate,stopi,stopj,stopstate,out)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [GeneWise21 *]
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
boolean read_hidden_GeneWise21(GeneWise21 * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out) 
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


      max_hidden_GeneWise21(mat,startj,i,j,state,isspecial,&i,&j,&state,&isspecial,&cellscore);  


      if( i == GeneWise21_READ_OFF_ERROR)    {  
        warn("In GeneWise21 hidden read off, between %d:%d,%d:%d - at got bad read off. Problem!",starti,startj,stopi,stopj);    
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
        warn("In GeneWise21 hidden read off, between %d:%d,%d:%d - hit start cell, but not in start state. Can't be good!.",starti,startj,stopi,stopj);  
        return FALSE;    
        }  
      }  
    warn("In GeneWise21 hidden read off, between %d:%d,%d:%d - gone past start cell (now in %d,%d,%d), can't be good news!.",starti,startj,stopi,stopj,i,j,state);   
    return FALSE;    
}    


/* Function:  max_hidden_GeneWise21(mat,hiddenj,i,j,state,isspecial,reti,retj,retstate,retspecial,cellscore)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [GeneWise21 *]
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
int max_hidden_GeneWise21(GeneWise21 * mat,int hiddenj,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore) 
{
    register int temp;   
    register int cscore; 


    *reti = (*retj) = (*retstate) = GeneWise21_READ_OFF_ERROR;   


    if( i < 0 || j < 0 || i > mat->query->len || j > mat->target->seq->len)  {  
      warn("In GeneWise21 matrix special read off - out of bounds on matrix [i,j is %d,%d state %d in standard matrix]",i,j,state);  
      return -1; 
      }  


    /* Then you have to select the correct switch statement to figure out the readoff      */ 
    /* Somewhat odd - reverse the order of calculation and return as soon as it is correct */ 
    cscore = GeneWise21_HIDDEN_MATRIX(mat,i,j,state);    
    switch(state)    { /*Switch state */ 
      case MATCH :   
        temp = cscore - (mat->gp21->transition[GP21_INSERT_2_BASE]) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j)); 
        if( temp == GeneWise21_HIDDEN_MATRIX(mat,i - 1,j - 5,MATCH) )    {  
          *reti = i - 1; 
          *retj = j - 5; 
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_HIDDEN_MATRIX(mat,i-1,j-5,MATCH);   
            }  
          return GeneWise21_HIDDEN_MATRIX(mat,i - 1,j - 5,MATCH);    
          }  
        temp = cscore - (mat->gp21->transition[GP21_INSERT_1_BASE]) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j)); 
        if( temp == GeneWise21_HIDDEN_MATRIX(mat,i - 1,j - 4,MATCH) )    {  
          *reti = i - 1; 
          *retj = j - 4; 
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_HIDDEN_MATRIX(mat,i-1,j-4,MATCH);   
            }  
          return GeneWise21_HIDDEN_MATRIX(mat,i - 1,j - 4,MATCH);    
          }  
        temp = cscore - (mat->gp21->transition[GP21_DELETE_2_BASE]) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j)); 
        if( temp == GeneWise21_HIDDEN_MATRIX(mat,i - 1,j - 1,MATCH) )    {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_HIDDEN_MATRIX(mat,i-1,j-1,MATCH);   
            }  
          return GeneWise21_HIDDEN_MATRIX(mat,i - 1,j - 1,MATCH);    
          }  
        temp = cscore - (mat->gp21->transition[GP21_DELETE_1_BASE]) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j)); 
        if( temp == GeneWise21_HIDDEN_MATRIX(mat,i - 1,j - 2,MATCH) )    {  
          *reti = i - 1; 
          *retj = j - 2; 
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_HIDDEN_MATRIX(mat,i-1,j-2,MATCH);   
            }  
          return GeneWise21_HIDDEN_MATRIX(mat,i - 1,j - 2,MATCH);    
          }  
        temp = cscore - ((mat->gp21->transition[GP21_PY2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-1)))) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j));  
        if( temp == GeneWise21_HIDDEN_MATRIX(mat,i - 1,j - 4,PY_2M) )    {  
          *reti = i - 1; 
          *retj = j - 4; 
          *retstate = PY_2M; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_HIDDEN_MATRIX(mat,i-1,j-4,PY_2M);   
            }  
          return GeneWise21_HIDDEN_MATRIX(mat,i - 1,j - 4,PY_2M);    
          }  
        temp = cscore - ((mat->gp21->transition[GP21_SPACER2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-1)))) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j));  
        if( temp == GeneWise21_HIDDEN_MATRIX(mat,i - 1,j - 4,SPACER_2M) )    {  
          *reti = i - 1; 
          *retj = j - 4; 
          *retstate = SPACER_2M; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_HIDDEN_MATRIX(mat,i-1,j-4,SPACER_2M);   
            }  
          return GeneWise21_HIDDEN_MATRIX(mat,i - 1,j - 4,SPACER_2M);    
          }  
        temp = cscore - ((mat->gp21->transition[GP21_PY2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-2)))) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j));  
        if( temp == GeneWise21_HIDDEN_MATRIX(mat,i - 1,j - 5,PY_1M) )    {  
          *reti = i - 1; 
          *retj = j - 5; 
          *retstate = PY_1M; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_HIDDEN_MATRIX(mat,i-1,j-5,PY_1M);   
            }  
          return GeneWise21_HIDDEN_MATRIX(mat,i - 1,j - 5,PY_1M);    
          }  
        temp = cscore - ((mat->gp21->transition[GP21_SPACER2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-2)))) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j));  
        if( temp == GeneWise21_HIDDEN_MATRIX(mat,i - 1,j - 5,SPACER_1M) )    {  
          *reti = i - 1; 
          *retj = j - 5; 
          *retstate = SPACER_1M; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_HIDDEN_MATRIX(mat,i-1,j-5,SPACER_1M);   
            }  
          return GeneWise21_HIDDEN_MATRIX(mat,i - 1,j - 5,SPACER_1M);    
          }  
        temp = cscore - ((((mat->gp21->transition[GP21_PY2CDS]+mat->query->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)])+CSEQ_GENOMIC_3SS(mat->target,(j-3)))+mat->query->seg[i]->transition[GW_MATCH_BALANCE_3SS])) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j));    
        if( temp == GeneWise21_HIDDEN_MATRIX(mat,i - 1,j - 6,PY_0M) )    {  
          *reti = i - 1; 
          *retj = j - 6; 
          *retstate = PY_0M; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_HIDDEN_MATRIX(mat,i-1,j-6,PY_0M);   
            }  
          return GeneWise21_HIDDEN_MATRIX(mat,i - 1,j - 6,PY_0M);    
          }  
        temp = cscore - ((((mat->gp21->transition[GP21_SPACER2CDS]+mat->query->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)])+CSEQ_GENOMIC_3SS(mat->target,(j-3)))+mat->query->seg[i]->transition[GW_MATCH_BALANCE_3SS])) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j));    
        if( temp == GeneWise21_HIDDEN_MATRIX(mat,i - 1,j - 6,SPACER_0M) )    {  
          *reti = i - 1; 
          *retj = j - 6; 
          *retstate = SPACER_0M; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_HIDDEN_MATRIX(mat,i-1,j-6,SPACER_0M);   
            }  
          return GeneWise21_HIDDEN_MATRIX(mat,i - 1,j - 6,SPACER_0M);    
          }  
        /* Not allowing special sources.. skipping BEFORE_CODON */ 
        temp = cscore - ((mat->query->seg[i]->transition[GW_DELETE2MATCH]+mat->query->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)])) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j));    
        if( temp == GeneWise21_HIDDEN_MATRIX(mat,i - 1,j - 3,DELETE) )   {  
          *reti = i - 1; 
          *retj = j - 3; 
          *retstate = DELETE;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_HIDDEN_MATRIX(mat,i-1,j-3,DELETE);  
            }  
          return GeneWise21_HIDDEN_MATRIX(mat,i - 1,j - 3,DELETE);   
          }  
        temp = cscore - ((mat->query->seg[i]->transition[GW_INSERT2MATCH]+mat->query->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)])) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j));    
        if( temp == GeneWise21_HIDDEN_MATRIX(mat,i - 1,j - 3,INSERT) )   {  
          *reti = i - 1; 
          *retj = j - 3; 
          *retstate = INSERT;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_HIDDEN_MATRIX(mat,i-1,j-3,INSERT);  
            }  
          return GeneWise21_HIDDEN_MATRIX(mat,i - 1,j - 3,INSERT);   
          }  
        temp = cscore - ((mat->query->seg[i]->transition[GW_MATCH2MATCH]+mat->query->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)])) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j)); 
        if( temp == GeneWise21_HIDDEN_MATRIX(mat,i - 1,j - 3,MATCH) )    {  
          *reti = i - 1; 
          *retj = j - 3; 
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_HIDDEN_MATRIX(mat,i-1,j-3,MATCH);   
            }  
          return GeneWise21_HIDDEN_MATRIX(mat,i - 1,j - 3,MATCH);    
          }  
        warn("Major problem (!) - in GeneWise21 read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      case INSERT :  
        temp = cscore - (mat->gp21->transition[GP21_INSERT_2_BASE]) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j)); 
        if( temp == GeneWise21_HIDDEN_MATRIX(mat,i - 1,j - 5,INSERT) )   {  
          *reti = i - 1; 
          *retj = j - 5; 
          *retstate = INSERT;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_HIDDEN_MATRIX(mat,i-1,j-5,INSERT);  
            }  
          return GeneWise21_HIDDEN_MATRIX(mat,i - 1,j - 5,INSERT);   
          }  
        temp = cscore - (mat->gp21->transition[GP21_INSERT_1_BASE]) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j)); 
        if( temp == GeneWise21_HIDDEN_MATRIX(mat,i - 1,j - 4,INSERT) )   {  
          *reti = i - 1; 
          *retj = j - 4; 
          *retstate = INSERT;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_HIDDEN_MATRIX(mat,i-1,j-4,INSERT);  
            }  
          return GeneWise21_HIDDEN_MATRIX(mat,i - 1,j - 4,INSERT);   
          }  
        temp = cscore - (mat->gp21->transition[GP21_DELETE_2_BASE]) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j)); 
        if( temp == GeneWise21_HIDDEN_MATRIX(mat,i - 1,j - 1,INSERT) )   {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = INSERT;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_HIDDEN_MATRIX(mat,i-1,j-1,INSERT);  
            }  
          return GeneWise21_HIDDEN_MATRIX(mat,i - 1,j - 1,INSERT);   
          }  
        temp = cscore - (mat->gp21->transition[GP21_DELETE_1_BASE]) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j)); 
        if( temp == GeneWise21_HIDDEN_MATRIX(mat,i - 1,j - 2,INSERT) )   {  
          *reti = i - 1; 
          *retj = j - 2; 
          *retstate = INSERT;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_HIDDEN_MATRIX(mat,i-1,j-2,INSERT);  
            }  
          return GeneWise21_HIDDEN_MATRIX(mat,i - 1,j - 2,INSERT);   
          }  
        temp = cscore - ((mat->gp21->transition[GP21_PY2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-1)))) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j));  
        if( temp == GeneWise21_HIDDEN_MATRIX(mat,i - 1,j - 4,PY_2I) )    {  
          *reti = i - 1; 
          *retj = j - 4; 
          *retstate = PY_2I; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_HIDDEN_MATRIX(mat,i-1,j-4,PY_2I);   
            }  
          return GeneWise21_HIDDEN_MATRIX(mat,i - 1,j - 4,PY_2I);    
          }  
        temp = cscore - ((mat->gp21->transition[GP21_PY2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-2)))) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j));  
        if( temp == GeneWise21_HIDDEN_MATRIX(mat,i - 1,j - 5,PY_1I) )    {  
          *reti = i - 1; 
          *retj = j - 5; 
          *retstate = PY_1I; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_HIDDEN_MATRIX(mat,i-1,j-5,PY_1I);   
            }  
          return GeneWise21_HIDDEN_MATRIX(mat,i - 1,j - 5,PY_1I);    
          }  
        temp = cscore - ((((mat->gp21->transition[GP21_PY2CDS]+mat->query->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)])+CSEQ_GENOMIC_3SS(mat->target,(j-3)))+mat->query->seg[i]->transition[GW_MATCH_BALANCE_3SS])) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j));   
        if( temp == GeneWise21_HIDDEN_MATRIX(mat,i - 1,j - 6,PY_0I) )    {  
          *reti = i - 1; 
          *retj = j - 6; 
          *retstate = PY_0I; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_HIDDEN_MATRIX(mat,i-1,j-6,PY_0I);   
            }  
          return GeneWise21_HIDDEN_MATRIX(mat,i - 1,j - 6,PY_0I);    
          }  
        temp = cscore - ((mat->gp21->transition[GP21_SPACER2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-1)))) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j));  
        if( temp == GeneWise21_HIDDEN_MATRIX(mat,i - 1,j - 4,SPACER_2I) )    {  
          *reti = i - 1; 
          *retj = j - 4; 
          *retstate = SPACER_2I; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_HIDDEN_MATRIX(mat,i-1,j-4,SPACER_2I);   
            }  
          return GeneWise21_HIDDEN_MATRIX(mat,i - 1,j - 4,SPACER_2I);    
          }  
        temp = cscore - ((mat->gp21->transition[GP21_SPACER2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-2)))) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j));  
        if( temp == GeneWise21_HIDDEN_MATRIX(mat,i - 1,j - 5,SPACER_1I) )    {  
          *reti = i - 1; 
          *retj = j - 5; 
          *retstate = SPACER_1I; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_HIDDEN_MATRIX(mat,i-1,j-5,SPACER_1I);   
            }  
          return GeneWise21_HIDDEN_MATRIX(mat,i - 1,j - 5,SPACER_1I);    
          }  
        temp = cscore - ((((mat->gp21->transition[GP21_SPACER2CDS]+mat->query->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)])+CSEQ_GENOMIC_3SS(mat->target,(j-3)))+mat->query->seg[i]->transition[GW_MATCH_BALANCE_3SS])) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j));   
        if( temp == GeneWise21_HIDDEN_MATRIX(mat,i - 1,j - 6,SPACER_0I) )    {  
          *reti = i - 1; 
          *retj = j - 6; 
          *retstate = SPACER_0I; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_HIDDEN_MATRIX(mat,i-1,j-6,SPACER_0I);   
            }  
          return GeneWise21_HIDDEN_MATRIX(mat,i - 1,j - 6,SPACER_0I);    
          }  
        /* Not allowing special sources.. skipping BEFORE_CODON */ 
        temp = cscore - ((mat->query->seg[i]->transition[GW_DELETE2INSERT]+mat->query->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)])) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j));  
        if( temp == GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 3,DELETE) )   {  
          *reti = i - 0; 
          *retj = j - 3; 
          *retstate = DELETE;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_HIDDEN_MATRIX(mat,i-0,j-3,DELETE);  
            }  
          return GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 3,DELETE);   
          }  
        temp = cscore - ((mat->query->seg[i]->transition[GW_INSERT2INSERT]+mat->query->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)])) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j));  
        if( temp == GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 3,INSERT) )   {  
          *reti = i - 0; 
          *retj = j - 3; 
          *retstate = INSERT;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_HIDDEN_MATRIX(mat,i-0,j-3,INSERT);  
            }  
          return GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 3,INSERT);   
          }  
        temp = cscore - ((mat->query->seg[i]->transition[GW_MATCH2INSERT]+mat->query->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)])) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j));   
        if( temp == GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 3,MATCH) )    {  
          *reti = i - 0; 
          *retj = j - 3; 
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_HIDDEN_MATRIX(mat,i-0,j-3,MATCH);   
            }  
          return GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 3,MATCH);    
          }  
        warn("Major problem (!) - in GeneWise21 read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      case DELETE :  
        /* Not allowing special sources.. skipping BEFORE_CODON */ 
        temp = cscore - (mat->query->seg[i]->transition[GW_DELETE2DELETE]) -  (0);   
        if( temp == GeneWise21_HIDDEN_MATRIX(mat,i - 1,j - 0,DELETE) )   {  
          *reti = i - 1; 
          *retj = j - 0; 
          *retstate = DELETE;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_HIDDEN_MATRIX(mat,i-1,j-0,DELETE);  
            }  
          return GeneWise21_HIDDEN_MATRIX(mat,i - 1,j - 0,DELETE);   
          }  
        temp = cscore - (mat->query->seg[i]->transition[GW_INSERT2DELETE]) -  (0);   
        if( temp == GeneWise21_HIDDEN_MATRIX(mat,i - 1,j - 0,INSERT) )   {  
          *reti = i - 1; 
          *retj = j - 0; 
          *retstate = INSERT;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_HIDDEN_MATRIX(mat,i-1,j-0,INSERT);  
            }  
          return GeneWise21_HIDDEN_MATRIX(mat,i - 1,j - 0,INSERT);   
          }  
        temp = cscore - (mat->query->seg[i]->transition[GW_MATCH2DELETE]) -  (0);    
        if( temp == GeneWise21_HIDDEN_MATRIX(mat,i - 1,j - 0,MATCH) )    {  
          *reti = i - 1; 
          *retj = j - 0; 
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_HIDDEN_MATRIX(mat,i-1,j-0,MATCH);   
            }  
          return GeneWise21_HIDDEN_MATRIX(mat,i - 1,j - 0,MATCH);    
          }  
        warn("Major problem (!) - in GeneWise21 read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      case CENTRAL_0M :  
        temp = cscore - ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CENTRAL2CENTRAL])) -  (0); 
        if( temp == GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 1,CENTRAL_0M) )   {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = CENTRAL_0M;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_HIDDEN_MATRIX(mat,i-0,j-1,CENTRAL_0M);  
            }  
          return GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 1,CENTRAL_0M);   
          }  
        temp = cscore - (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_DELETE2MATCH])) -  (0);   
        if( temp == GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 8,DELETE) )   {  
          *reti = i - 0; 
          *retj = j - 8; 
          *retstate = DELETE;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_HIDDEN_MATRIX(mat,i-0,j-8,DELETE);  
            }  
          return GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 8,DELETE);   
          }  
        temp = cscore - ((((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_INSERT2MATCH])+mat->query->seg[i]->transition[GW_INSERT_BALANCE_5SS])) -  (0);   
        if( temp == GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 8,INSERT) )   {  
          *reti = i - 0; 
          *retj = j - 8; 
          *retstate = INSERT;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_HIDDEN_MATRIX(mat,i-0,j-8,INSERT);  
            }  
          return GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 8,INSERT);   
          }  
        temp = cscore - ((((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_MATCH2MATCH])+mat->query->seg[i]->transition[GW_MATCH_BALANCE_5SS])) -  (0); 
        if( temp == GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 8,MATCH) )    {  
          *reti = i - 0; 
          *retj = j - 8; 
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_HIDDEN_MATRIX(mat,i-0,j-8,MATCH);   
            }  
          return GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 8,MATCH);    
          }  
        warn("Major problem (!) - in GeneWise21 read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      case PY_0M :   
        temp = cscore - (mat->gp21->transition[GP21_PY2PY]) -  (mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)]);    
        if( temp == GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 1,PY_0M) )    {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = PY_0M; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_HIDDEN_MATRIX(mat,i-0,j-1,PY_0M);   
            }  
          return GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 1,PY_0M);    
          }  
        temp = cscore - (mat->gp21->transition[GP21_CENTRAL2PY]) -  (mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)]);   
        if( temp == GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 1,CENTRAL_0M) )   {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = CENTRAL_0M;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_HIDDEN_MATRIX(mat,i-0,j-1,CENTRAL_0M);  
            }  
          return GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 1,CENTRAL_0M);   
          }  
        warn("Major problem (!) - in GeneWise21 read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      case SPACER_0M :   
        temp = cscore - (mat->gp21->transition[GP21_SPACER2SPACER]) -  (mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)]);    
        if( temp == GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 1,SPACER_0M) )    {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = SPACER_0M; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_HIDDEN_MATRIX(mat,i-0,j-1,SPACER_0M);   
            }  
          return GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 1,SPACER_0M);    
          }  
        temp = cscore - (mat->gp21->transition[GP21_PY2SPACER]) -  (mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)]);    
        if( temp == GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 1,PY_0M) )    {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = PY_0M; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_HIDDEN_MATRIX(mat,i-0,j-1,PY_0M);   
            }  
          return GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 1,PY_0M);    
          }  
        warn("Major problem (!) - in GeneWise21 read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      case CENTRAL_1M :  
        temp = cscore - ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CENTRAL2CENTRAL])) -  (0); 
        if( temp == GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 1,CENTRAL_1M) )   {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = CENTRAL_1M;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_HIDDEN_MATRIX(mat,i-0,j-1,CENTRAL_1M);  
            }  
          return GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 1,CENTRAL_1M);   
          }  
        temp = cscore - (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_DELETE2MATCH])) -  (0);   
        if( temp == GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 9,DELETE) )   {  
          *reti = i - 0; 
          *retj = j - 9; 
          *retstate = DELETE;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_HIDDEN_MATRIX(mat,i-0,j-9,DELETE);  
            }  
          return GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 9,DELETE);   
          }  
        temp = cscore - ((((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_INSERT2MATCH])+mat->query->seg[i]->transition[GW_MATCH_BALANCE_5SS])) -  (0);    
        if( temp == GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 9,INSERT) )   {  
          *reti = i - 0; 
          *retj = j - 9; 
          *retstate = INSERT;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_HIDDEN_MATRIX(mat,i-0,j-9,INSERT);  
            }  
          return GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 9,INSERT);   
          }  
        temp = cscore - ((((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_MATCH2MATCH])+mat->query->seg[i]->transition[GW_MATCH_BALANCE_5SS])) -  (0); 
        if( temp == GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 9,MATCH) )    {  
          *reti = i - 0; 
          *retj = j - 9; 
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_HIDDEN_MATRIX(mat,i-0,j-9,MATCH);   
            }  
          return GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 9,MATCH);    
          }  
        warn("Major problem (!) - in GeneWise21 read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      case PY_1M :   
        temp = cscore - (mat->gp21->transition[GP21_PY2PY]) -  (mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)]);    
        if( temp == GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 1,PY_1M) )    {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = PY_1M; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_HIDDEN_MATRIX(mat,i-0,j-1,PY_1M);   
            }  
          return GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 1,PY_1M);    
          }  
        temp = cscore - (mat->gp21->transition[GP21_CENTRAL2PY]) -  (mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)]);   
        if( temp == GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 1,CENTRAL_1M) )   {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = CENTRAL_1M;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_HIDDEN_MATRIX(mat,i-0,j-1,CENTRAL_1M);  
            }  
          return GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 1,CENTRAL_1M);   
          }  
        warn("Major problem (!) - in GeneWise21 read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      case SPACER_1M :   
        temp = cscore - (mat->gp21->transition[GP21_SPACER2SPACER]) -  (mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)]);    
        if( temp == GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 1,SPACER_1M) )    {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = SPACER_1M; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_HIDDEN_MATRIX(mat,i-0,j-1,SPACER_1M);   
            }  
          return GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 1,SPACER_1M);    
          }  
        temp = cscore - (mat->gp21->transition[GP21_PY2SPACER]) -  (mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)]);    
        if( temp == GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 1,PY_1M) )    {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = PY_1M; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_HIDDEN_MATRIX(mat,i-0,j-1,PY_1M);   
            }  
          return GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 1,PY_1M);    
          }  
        warn("Major problem (!) - in GeneWise21 read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      case CENTRAL_2M :  
        temp = cscore - ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CENTRAL2CENTRAL])) -  (0); 
        if( temp == GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 1,CENTRAL_2M) )   {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = CENTRAL_2M;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_HIDDEN_MATRIX(mat,i-0,j-1,CENTRAL_2M);  
            }  
          return GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 1,CENTRAL_2M);   
          }  
        temp = cscore - (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_DELETE2MATCH])) -  (0);   
        if( temp == GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 10,DELETE) )  {  
          *reti = i - 0; 
          *retj = j - 10;    
          *retstate = DELETE;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_HIDDEN_MATRIX(mat,i-0,j-10,DELETE); 
            }  
          return GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 10,DELETE);  
          }  
        temp = cscore - ((((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_INSERT2MATCH])+mat->query->seg[i]->transition[GW_MATCH_BALANCE_5SS])) -  (0);    
        if( temp == GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 10,INSERT) )  {  
          *reti = i - 0; 
          *retj = j - 10;    
          *retstate = INSERT;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_HIDDEN_MATRIX(mat,i-0,j-10,INSERT); 
            }  
          return GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 10,INSERT);  
          }  
        temp = cscore - ((((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_MATCH2MATCH])+mat->query->seg[i]->transition[GW_MATCH_BALANCE_5SS])) -  (0); 
        if( temp == GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 10,MATCH) )   {  
          *reti = i - 0; 
          *retj = j - 10;    
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_HIDDEN_MATRIX(mat,i-0,j-10,MATCH);  
            }  
          return GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 10,MATCH);   
          }  
        warn("Major problem (!) - in GeneWise21 read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      case PY_2M :   
        temp = cscore - (mat->gp21->transition[GP21_PY2PY]) -  (mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)]);    
        if( temp == GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 1,PY_2M) )    {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = PY_2M; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_HIDDEN_MATRIX(mat,i-0,j-1,PY_2M);   
            }  
          return GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 1,PY_2M);    
          }  
        temp = cscore - (mat->gp21->transition[GP21_CENTRAL2PY]) -  (mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)]);   
        if( temp == GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 1,CENTRAL_2M) )   {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = CENTRAL_2M;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_HIDDEN_MATRIX(mat,i-0,j-1,CENTRAL_2M);  
            }  
          return GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 1,CENTRAL_2M);   
          }  
        warn("Major problem (!) - in GeneWise21 read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      case SPACER_2M :   
        temp = cscore - (mat->gp21->transition[GP21_SPACER2SPACER]) -  (mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)]);    
        if( temp == GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 1,SPACER_2M) )    {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = SPACER_2M; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_HIDDEN_MATRIX(mat,i-0,j-1,SPACER_2M);   
            }  
          return GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 1,SPACER_2M);    
          }  
        temp = cscore - (mat->gp21->transition[GP21_PY2SPACER]) -  (mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)]);    
        if( temp == GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 1,PY_2M) )    {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = PY_2M; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_HIDDEN_MATRIX(mat,i-0,j-1,PY_2M);   
            }  
          return GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 1,PY_2M);    
          }  
        warn("Major problem (!) - in GeneWise21 read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      case CENTRAL_0I :  
        temp = cscore - ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CENTRAL2CENTRAL])) -  (0); 
        if( temp == GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 1,CENTRAL_0I) )   {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = CENTRAL_0I;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_HIDDEN_MATRIX(mat,i-0,j-1,CENTRAL_0I);  
            }  
          return GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 1,CENTRAL_0I);   
          }  
        temp = cscore - (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_DELETE2INSERT])) -  (0);  
        if( temp == GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 8,DELETE) )   {  
          *reti = i - 0; 
          *retj = j - 8; 
          *retstate = DELETE;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_HIDDEN_MATRIX(mat,i-0,j-8,DELETE);  
            }  
          return GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 8,DELETE);   
          }  
        temp = cscore - ((((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_INSERT2INSERT])+mat->query->seg[i]->transition[GW_INSERT_BALANCE_5SS])) -  (0);  
        if( temp == GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 8,INSERT) )   {  
          *reti = i - 0; 
          *retj = j - 8; 
          *retstate = INSERT;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_HIDDEN_MATRIX(mat,i-0,j-8,INSERT);  
            }  
          return GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 8,INSERT);   
          }  
        temp = cscore - ((((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_MATCH2INSERT])+mat->query->seg[i]->transition[GW_MATCH_BALANCE_5SS])) -  (0);    
        if( temp == GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 8,MATCH) )    {  
          *reti = i - 0; 
          *retj = j - 8; 
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_HIDDEN_MATRIX(mat,i-0,j-8,MATCH);   
            }  
          return GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 8,MATCH);    
          }  
        warn("Major problem (!) - in GeneWise21 read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      case PY_0I :   
        temp = cscore - (mat->gp21->transition[GP21_PY2PY]) -  (mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)]);    
        if( temp == GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 1,PY_0I) )    {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = PY_0I; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_HIDDEN_MATRIX(mat,i-0,j-1,PY_0I);   
            }  
          return GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 1,PY_0I);    
          }  
        temp = cscore - (mat->gp21->transition[GP21_CENTRAL2PY]) -  (mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)]);   
        if( temp == GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 1,CENTRAL_0I) )   {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = CENTRAL_0I;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_HIDDEN_MATRIX(mat,i-0,j-1,CENTRAL_0I);  
            }  
          return GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 1,CENTRAL_0I);   
          }  
        warn("Major problem (!) - in GeneWise21 read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      case SPACER_0I :   
        temp = cscore - (mat->gp21->transition[GP21_SPACER2SPACER]) -  (mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)]);    
        if( temp == GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 1,SPACER_0I) )    {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = SPACER_0I; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_HIDDEN_MATRIX(mat,i-0,j-1,SPACER_0I);   
            }  
          return GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 1,SPACER_0I);    
          }  
        temp = cscore - (mat->gp21->transition[GP21_PY2SPACER]) -  (mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)]);    
        if( temp == GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 1,PY_0I) )    {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = PY_0I; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_HIDDEN_MATRIX(mat,i-0,j-1,PY_0I);   
            }  
          return GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 1,PY_0I);    
          }  
        warn("Major problem (!) - in GeneWise21 read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      case CENTRAL_1I :  
        temp = cscore - ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CENTRAL2CENTRAL])) -  (0); 
        if( temp == GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 1,CENTRAL_1I) )   {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = CENTRAL_1I;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_HIDDEN_MATRIX(mat,i-0,j-1,CENTRAL_1I);  
            }  
          return GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 1,CENTRAL_1I);   
          }  
        temp = cscore - (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_DELETE2INSERT])) -  (0);  
        if( temp == GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 9,DELETE) )   {  
          *reti = i - 0; 
          *retj = j - 9; 
          *retstate = DELETE;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_HIDDEN_MATRIX(mat,i-0,j-9,DELETE);  
            }  
          return GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 9,DELETE);   
          }  
        temp = cscore - (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_INSERT2INSERT])) -  (0);  
        if( temp == GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 9,INSERT) )   {  
          *reti = i - 0; 
          *retj = j - 9; 
          *retstate = INSERT;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_HIDDEN_MATRIX(mat,i-0,j-9,INSERT);  
            }  
          return GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 9,INSERT);   
          }  
        temp = cscore - (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_MATCH2INSERT])) -  (0);   
        if( temp == GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 9,MATCH) )    {  
          *reti = i - 0; 
          *retj = j - 9; 
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_HIDDEN_MATRIX(mat,i-0,j-9,MATCH);   
            }  
          return GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 9,MATCH);    
          }  
        warn("Major problem (!) - in GeneWise21 read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      case PY_1I :   
        temp = cscore - (mat->gp21->transition[GP21_PY2PY]) -  (mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)]);    
        if( temp == GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 1,PY_1I) )    {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = PY_1I; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_HIDDEN_MATRIX(mat,i-0,j-1,PY_1I);   
            }  
          return GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 1,PY_1I);    
          }  
        temp = cscore - (mat->gp21->transition[GP21_CENTRAL2PY]) -  (mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)]);   
        if( temp == GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 1,CENTRAL_1I) )   {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = CENTRAL_1I;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_HIDDEN_MATRIX(mat,i-0,j-1,CENTRAL_1I);  
            }  
          return GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 1,CENTRAL_1I);   
          }  
        warn("Major problem (!) - in GeneWise21 read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      case SPACER_1I :   
        temp = cscore - (mat->gp21->transition[GP21_SPACER2SPACER]) -  (mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)]);    
        if( temp == GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 1,SPACER_1I) )    {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = SPACER_1I; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_HIDDEN_MATRIX(mat,i-0,j-1,SPACER_1I);   
            }  
          return GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 1,SPACER_1I);    
          }  
        temp = cscore - (mat->gp21->transition[GP21_PY2SPACER]) -  (mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)]);    
        if( temp == GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 1,PY_1I) )    {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = PY_1I; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_HIDDEN_MATRIX(mat,i-0,j-1,PY_1I);   
            }  
          return GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 1,PY_1I);    
          }  
        warn("Major problem (!) - in GeneWise21 read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      case CENTRAL_2I :  
        temp = cscore - ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CENTRAL2CENTRAL])) -  (0); 
        if( temp == GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 1,CENTRAL_2I) )   {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = CENTRAL_2I;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_HIDDEN_MATRIX(mat,i-0,j-1,CENTRAL_2I);  
            }  
          return GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 1,CENTRAL_2I);   
          }  
        temp = cscore - (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_DELETE2INSERT])) -  (0);  
        if( temp == GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 10,INSERT) )  {  
          *reti = i - 0; 
          *retj = j - 10;    
          *retstate = INSERT;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_HIDDEN_MATRIX(mat,i-0,j-10,INSERT); 
            }  
          return GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 10,INSERT);  
          }  
        temp = cscore - (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_INSERT2INSERT])) -  (0);  
        if( temp == GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 10,INSERT) )  {  
          *reti = i - 0; 
          *retj = j - 10;    
          *retstate = INSERT;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_HIDDEN_MATRIX(mat,i-0,j-10,INSERT); 
            }  
          return GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 10,INSERT);  
          }  
        temp = cscore - (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_MATCH2INSERT])) -  (0);   
        if( temp == GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 10,MATCH) )   {  
          *reti = i - 0; 
          *retj = j - 10;    
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_HIDDEN_MATRIX(mat,i-0,j-10,MATCH);  
            }  
          return GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 10,MATCH);   
          }  
        warn("Major problem (!) - in GeneWise21 read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      case PY_2I :   
        temp = cscore - (mat->gp21->transition[GP21_PY2PY]) -  (mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)]);    
        if( temp == GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 1,PY_2I) )    {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = PY_2I; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_HIDDEN_MATRIX(mat,i-0,j-1,PY_2I);   
            }  
          return GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 1,PY_2I);    
          }  
        temp = cscore - (mat->gp21->transition[GP21_CENTRAL2PY]) -  (mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)]);   
        if( temp == GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 1,CENTRAL_2I) )   {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = CENTRAL_2I;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_HIDDEN_MATRIX(mat,i-0,j-1,CENTRAL_2I);  
            }  
          return GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 1,CENTRAL_2I);   
          }  
        warn("Major problem (!) - in GeneWise21 read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      case SPACER_2I :   
        temp = cscore - (mat->gp21->transition[GP21_SPACER2SPACER]) -  (mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)]);    
        if( temp == GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 1,SPACER_2I) )    {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = SPACER_2I; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_HIDDEN_MATRIX(mat,i-0,j-1,SPACER_2I);   
            }  
          return GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 1,SPACER_2I);    
          }  
        temp = cscore - (mat->gp21->transition[GP21_PY2SPACER]) -  (mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)]);    
        if( temp == GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 1,PY_2I) )    {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = PY_2I; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_HIDDEN_MATRIX(mat,i-0,j-1,PY_2I);   
            }  
          return GeneWise21_HIDDEN_MATRIX(mat,i - 0,j - 1,PY_2I);    
          }  
        warn("Major problem (!) - in GeneWise21 read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      default:   
        warn("Major problem (!) - in GeneWise21 read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      } /* end of Switch state  */ 
}    


/* Function:  read_special_strip_GeneWise21(mat,stopi,stopj,stopstate,startj,startstate,out)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [GeneWise21 *]
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
boolean read_special_strip_GeneWise21(GeneWise21 * mat,int stopi,int stopj,int stopstate,int * startj,int * startstate,PackAln * out) 
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
    while( j > GeneWise21_DC_SHADOW_SPECIAL_SP(mat,i,j,state,4) && state != START)   { /*while more specials to eat up*/ 
      /* Put away current state, if we should */ 
      if(out != NULL)    {  
        pau = PackAlnUnit_alloc();  /* Should deal with memory overflow */ 
        pau->i = i;  
        pau->j = j;  
        pau->state =  state + 21;    
        add_PackAln(out,pau);    
        }  


      max_special_strip_GeneWise21(mat,i,j,state,isspecial,&i,&j,&state,&isspecial,&cellscore);  
      if( i == GeneWise21_READ_OFF_ERROR)    {  
        warn("In special strip read GeneWise21, got a bad read off error. Sorry!");  
        return FALSE;    
        }  
      } /* end of while more specials to eat up */ 


    /* check to see we have not gone too far! */ 
    if( state != START && j < GeneWise21_DC_SHADOW_SPECIAL_SP(mat,i,j,state,4))  {  
      warn("In special strip read GeneWise21, at special [%d] state [%d] overshot!",j,state);    
      return FALSE;  
      }  
    /* Put away last state */ 
    if(out != NULL)  {  
      pau = PackAlnUnit_alloc();/* Should deal with memory overflow */ 
      pau->i = i;    
      pau->j = j;    
      pau->state =  state + 21;  
      add_PackAln(out,pau);  
      }  


    /* Put away where we are in startj and startstate */ 
    *startj = j; 
    *startstate = state; 
    return TRUE; 
}    


/* Function:  max_special_strip_GeneWise21(mat,i,j,state,isspecial,reti,retj,retstate,retspecial,cellscore)
 *
 * Descrip:    A pretty intense internal function. Deals with read-off only in specials
 *
 *
 * Arg:               mat [UNKN ] Undocumented argument [GeneWise21 *]
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
int max_special_strip_GeneWise21(GeneWise21 * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore) 
{
    int temp;    
    int cscore;  


    *reti = (*retj) = (*retstate) = GeneWise21_READ_OFF_ERROR;   
    if( isspecial == FALSE ) {  
      warn("In special strip max function for GeneWise21, got a non special start point. Problem! (bad!)");  
      return (-1);   
      }  


    if( j < 0 || j > mat->target->seq->len)  {  
      warn("In GeneWise21 matrix special read off - out of bounds on matrix [j is %d in special]",j);    
      return -1; 
      }  


    cscore = GeneWise21_DC_SHADOW_SPECIAL(mat,i,j,state);    
    switch(state)    { /*switch on special states*/ 
      case START :   
      case END :     
        /* source AFTER_RND is a special */ 
        temp = cscore - (0) - (0);   
        if( temp == GeneWise21_DC_SHADOW_SPECIAL(mat,i - 0,j - 1,AFTER_RND) )    {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = AFTER_RND; 
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_DC_SHADOW_SPECIAL(mat,i-0,j-1,AFTER_RND);   
            }  
          return GeneWise21_DC_SHADOW_MATRIX(mat,i - 0,j - 1,AFTER_RND) ;    
          }  
      case BEFORE_RND :  
        /* source BEFORE_CODON is a special */ 
        temp = cscore - ((mat->rndbase->base[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CDS2RND])) - (0);  
        if( temp == GeneWise21_DC_SHADOW_SPECIAL(mat,i - 0,j - 1,BEFORE_CODON) ) {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = BEFORE_CODON;  
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_DC_SHADOW_SPECIAL(mat,i-0,j-1,BEFORE_CODON);    
            }  
          return GeneWise21_DC_SHADOW_MATRIX(mat,i - 0,j - 1,BEFORE_CODON) ;     
          }  
        /* source BEFORE_RND is a special */ 
        temp = cscore - ((mat->rndbase->base[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_RND2RND])) - (0);  
        if( temp == GeneWise21_DC_SHADOW_SPECIAL(mat,i - 0,j - 1,BEFORE_RND) )   {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = BEFORE_RND;    
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_DC_SHADOW_SPECIAL(mat,i-0,j-1,BEFORE_RND);  
            }  
          return GeneWise21_DC_SHADOW_MATRIX(mat,i - 0,j - 1,BEFORE_RND) ;   
          }  
        /* source START is a special */ 
        temp = cscore - (mat->rndbase->base[CSEQ_GENOMIC_BASE(mat->target,j)]) - (0);    
        if( temp == GeneWise21_DC_SHADOW_SPECIAL(mat,i - 0,j - 1,START) )    {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = START; 
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_DC_SHADOW_SPECIAL(mat,i-0,j-1,START);   
            }  
          return GeneWise21_DC_SHADOW_MATRIX(mat,i - 0,j - 1,START) ;    
          }  
      case BEFORE_CODON :    
        /* source BEFORE_SPACER is a special */ 
        temp = cscore - ((mat->gp21->transition[GP21_SPACER2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-1)))) - (0);    
        if( temp == GeneWise21_DC_SHADOW_SPECIAL(mat,i - 0,j - 3,BEFORE_SPACER) )    {  
          *reti = i - 0; 
          *retj = j - 3; 
          *retstate = BEFORE_SPACER; 
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_DC_SHADOW_SPECIAL(mat,i-0,j-3,BEFORE_SPACER);   
            }  
          return GeneWise21_DC_SHADOW_MATRIX(mat,i - 0,j - 3,BEFORE_SPACER) ;    
          }  
        /* source BEFORE_SPACER is a special */ 
        temp = cscore - ((mat->gp21->transition[GP21_SPACER2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-2)))) - (0);    
        if( temp == GeneWise21_DC_SHADOW_SPECIAL(mat,i - 0,j - 4,BEFORE_SPACER) )    {  
          *reti = i - 0; 
          *retj = j - 4; 
          *retstate = BEFORE_SPACER; 
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_DC_SHADOW_SPECIAL(mat,i-0,j-4,BEFORE_SPACER);   
            }  
          return GeneWise21_DC_SHADOW_MATRIX(mat,i - 0,j - 4,BEFORE_SPACER) ;    
          }  
        /* source BEFORE_SPACER is a special */ 
        temp = cscore - (((mat->gp21->transition[GP21_SPACER2CDS]+mat->rndcodon->codon[CSEQ_GENOMIC_CODON(mat->target,j)])+CSEQ_GENOMIC_3SS(mat->target,(j-3)))) - (0);  
        if( temp == GeneWise21_DC_SHADOW_SPECIAL(mat,i - 0,j - 5,BEFORE_SPACER) )    {  
          *reti = i - 0; 
          *retj = j - 5; 
          *retstate = BEFORE_SPACER; 
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_DC_SHADOW_SPECIAL(mat,i-0,j-5,BEFORE_SPACER);   
            }  
          return GeneWise21_DC_SHADOW_MATRIX(mat,i - 0,j - 5,BEFORE_SPACER) ;    
          }  
        /* source BEFORE_CODON is a special */ 
        temp = cscore - ((mat->rndcodon->codon[CSEQ_GENOMIC_CODON(mat->target,j)]+mat->gp21->transition[GP21_CDS2CDS])) - (0);   
        if( temp == GeneWise21_DC_SHADOW_SPECIAL(mat,i - 0,j - 3,BEFORE_CODON) ) {  
          *reti = i - 0; 
          *retj = j - 3; 
          *retstate = BEFORE_CODON;  
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_DC_SHADOW_SPECIAL(mat,i-0,j-3,BEFORE_CODON);    
            }  
          return GeneWise21_DC_SHADOW_MATRIX(mat,i - 0,j - 3,BEFORE_CODON) ;     
          }  
        /* source BEFORE_RND is a special */ 
        temp = cscore - ((mat->rndcodon->codon[CSEQ_GENOMIC_CODON(mat->target,j)]+mat->gp21->transition[GP21_RND2CDS])) - (0);   
        if( temp == GeneWise21_DC_SHADOW_SPECIAL(mat,i - 0,j - 3,BEFORE_RND) )   {  
          *reti = i - 0; 
          *retj = j - 3; 
          *retstate = BEFORE_RND;    
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_DC_SHADOW_SPECIAL(mat,i-0,j-3,BEFORE_RND);  
            }  
          return GeneWise21_DC_SHADOW_MATRIX(mat,i - 0,j - 3,BEFORE_RND) ;   
          }  
      case BEFORE_CENTRAL :  
        /* source BEFORE_CODON is a special */ 
        temp = cscore - ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))) - (0);  
        if( temp == GeneWise21_DC_SHADOW_SPECIAL(mat,i - 0,j - 10,BEFORE_CODON) )    {  
          *reti = i - 0; 
          *retj = j - 10;    
          *retstate = BEFORE_CODON;  
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_DC_SHADOW_SPECIAL(mat,i-0,j-10,BEFORE_CODON);   
            }  
          return GeneWise21_DC_SHADOW_MATRIX(mat,i - 0,j - 10,BEFORE_CODON) ;    
          }  
        /* source BEFORE_CODON is a special */ 
        temp = cscore - ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))) - (0);  
        if( temp == GeneWise21_DC_SHADOW_SPECIAL(mat,i - 0,j - 9,BEFORE_CODON) ) {  
          *reti = i - 0; 
          *retj = j - 9; 
          *retstate = BEFORE_CODON;  
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_DC_SHADOW_SPECIAL(mat,i-0,j-9,BEFORE_CODON);    
            }  
          return GeneWise21_DC_SHADOW_MATRIX(mat,i - 0,j - 9,BEFORE_CODON) ;     
          }  
        /* source BEFORE_CODON is a special */ 
        temp = cscore - ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))) - (0);  
        if( temp == GeneWise21_DC_SHADOW_SPECIAL(mat,i - 0,j - 8,BEFORE_CODON) ) {  
          *reti = i - 0; 
          *retj = j - 8; 
          *retstate = BEFORE_CODON;  
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_DC_SHADOW_SPECIAL(mat,i-0,j-8,BEFORE_CODON);    
            }  
          return GeneWise21_DC_SHADOW_MATRIX(mat,i - 0,j - 8,BEFORE_CODON) ;     
          }  
        /* source BEFORE_CENTRAL is a special */ 
        temp = cscore - ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CENTRAL2CENTRAL])) - (0);  
        if( temp == GeneWise21_DC_SHADOW_SPECIAL(mat,i - 0,j - 1,BEFORE_CENTRAL) )   {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = BEFORE_CENTRAL;    
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_DC_SHADOW_SPECIAL(mat,i-0,j-1,BEFORE_CENTRAL);  
            }  
          return GeneWise21_DC_SHADOW_MATRIX(mat,i - 0,j - 1,BEFORE_CENTRAL) ;   
          }  
      case BEFORE_PY_TRACT :     
        /* source BEFORE_PY_TRACT is a special */ 
        temp = cscore - (mat->gp21->transition[GP21_PY2PY]) - (mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)]);     
        if( temp == GeneWise21_DC_SHADOW_SPECIAL(mat,i - 0,j - 1,BEFORE_PY_TRACT) )  {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = BEFORE_PY_TRACT;   
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_DC_SHADOW_SPECIAL(mat,i-0,j-1,BEFORE_PY_TRACT);     
            }  
          return GeneWise21_DC_SHADOW_MATRIX(mat,i - 0,j - 1,BEFORE_PY_TRACT) ;  
          }  
        /* source BEFORE_CENTRAL is a special */ 
        temp = cscore - (mat->gp21->transition[GP21_CENTRAL2PY]) - (mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)]);    
        if( temp == GeneWise21_DC_SHADOW_SPECIAL(mat,i - 0,j - 1,BEFORE_CENTRAL) )   {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = BEFORE_CENTRAL;    
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_DC_SHADOW_SPECIAL(mat,i-0,j-1,BEFORE_CENTRAL);  
            }  
          return GeneWise21_DC_SHADOW_MATRIX(mat,i - 0,j - 1,BEFORE_CENTRAL) ;   
          }  
      case BEFORE_SPACER :   
        /* source BEFORE_SPACER is a special */ 
        temp = cscore - (mat->gp21->transition[GP21_SPACER2SPACER]) - (mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)]);     
        if( temp == GeneWise21_DC_SHADOW_SPECIAL(mat,i - 0,j - 1,BEFORE_SPACER) )    {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = BEFORE_SPACER; 
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_DC_SHADOW_SPECIAL(mat,i-0,j-1,BEFORE_SPACER);   
            }  
          return GeneWise21_DC_SHADOW_MATRIX(mat,i - 0,j - 1,BEFORE_SPACER) ;    
          }  
        /* source BEFORE_PY_TRACT is a special */ 
        temp = cscore - (mat->gp21->transition[GP21_PY2SPACER]) - (mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)]);     
        if( temp == GeneWise21_DC_SHADOW_SPECIAL(mat,i - 0,j - 1,BEFORE_PY_TRACT) )  {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = BEFORE_PY_TRACT;   
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_DC_SHADOW_SPECIAL(mat,i-0,j-1,BEFORE_PY_TRACT);     
            }  
          return GeneWise21_DC_SHADOW_MATRIX(mat,i - 0,j - 1,BEFORE_PY_TRACT) ;  
          }  
      case AFTER_RND :   
        /* source AFTER_CODON is a special */ 
        temp = cscore - ((mat->rndbase->base[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CDS2RND])) - (0);  
        if( temp == GeneWise21_DC_SHADOW_SPECIAL(mat,i - 0,j - 1,AFTER_CODON) )  {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = AFTER_CODON;   
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_DC_SHADOW_SPECIAL(mat,i-0,j-1,AFTER_CODON);     
            }  
          return GeneWise21_DC_SHADOW_MATRIX(mat,i - 0,j - 1,AFTER_CODON) ;  
          }  
        /* source AFTER_RND is a special */ 
        temp = cscore - ((mat->rndbase->base[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_RND2RND])) - (0);  
        if( temp == GeneWise21_DC_SHADOW_SPECIAL(mat,i - 0,j - 1,AFTER_RND) )    {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = AFTER_RND; 
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_DC_SHADOW_SPECIAL(mat,i-0,j-1,AFTER_RND);   
            }  
          return GeneWise21_DC_SHADOW_MATRIX(mat,i - 0,j - 1,AFTER_RND) ;    
          }  
      case AFTER_CODON :     
        /* source AFTER_SPACER is a special */ 
        temp = cscore - ((mat->gp21->transition[GP21_SPACER2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-1)))) - (0);    
        if( temp == GeneWise21_DC_SHADOW_SPECIAL(mat,i - 0,j - 3,AFTER_SPACER) ) {  
          *reti = i - 0; 
          *retj = j - 3; 
          *retstate = AFTER_SPACER;  
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_DC_SHADOW_SPECIAL(mat,i-0,j-3,AFTER_SPACER);    
            }  
          return GeneWise21_DC_SHADOW_MATRIX(mat,i - 0,j - 3,AFTER_SPACER) ;     
          }  
        /* source AFTER_SPACER is a special */ 
        temp = cscore - ((mat->gp21->transition[GP21_SPACER2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-2)))) - (0);    
        if( temp == GeneWise21_DC_SHADOW_SPECIAL(mat,i - 0,j - 4,AFTER_SPACER) ) {  
          *reti = i - 0; 
          *retj = j - 4; 
          *retstate = AFTER_SPACER;  
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_DC_SHADOW_SPECIAL(mat,i-0,j-4,AFTER_SPACER);    
            }  
          return GeneWise21_DC_SHADOW_MATRIX(mat,i - 0,j - 4,AFTER_SPACER) ;     
          }  
        /* source AFTER_SPACER is a special */ 
        temp = cscore - (((mat->gp21->transition[GP21_SPACER2CDS]+mat->rndcodon->codon[CSEQ_GENOMIC_CODON(mat->target,j)])+CSEQ_GENOMIC_3SS(mat->target,(j-3)))) - (0);  
        if( temp == GeneWise21_DC_SHADOW_SPECIAL(mat,i - 0,j - 5,AFTER_SPACER) ) {  
          *reti = i - 0; 
          *retj = j - 5; 
          *retstate = AFTER_SPACER;  
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_DC_SHADOW_SPECIAL(mat,i-0,j-5,AFTER_SPACER);    
            }  
          return GeneWise21_DC_SHADOW_MATRIX(mat,i - 0,j - 5,AFTER_SPACER) ;     
          }  
        /* source AFTER_CODON is a special */ 
        temp = cscore - ((mat->rndcodon->codon[CSEQ_GENOMIC_CODON(mat->target,j)]+mat->gp21->transition[GP21_CDS2CDS])) - (0);   
        if( temp == GeneWise21_DC_SHADOW_SPECIAL(mat,i - 0,j - 3,AFTER_CODON) )  {  
          *reti = i - 0; 
          *retj = j - 3; 
          *retstate = AFTER_CODON;   
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_DC_SHADOW_SPECIAL(mat,i-0,j-3,AFTER_CODON);     
            }  
          return GeneWise21_DC_SHADOW_MATRIX(mat,i - 0,j - 3,AFTER_CODON) ;  
          }  
        /* Source DELETE is not a special */ 
        /* Source INSERT is not a special */ 
        /* Source MATCH is not a special */ 
        /* source AFTER_RND is a special */ 
        temp = cscore - ((mat->rndcodon->codon[CSEQ_GENOMIC_CODON(mat->target,j)]+mat->gp21->transition[GP21_RND2CDS])) - (0);   
        if( temp == GeneWise21_DC_SHADOW_SPECIAL(mat,i - 0,j - 3,AFTER_RND) )    {  
          *reti = i - 0; 
          *retj = j - 3; 
          *retstate = AFTER_RND; 
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_DC_SHADOW_SPECIAL(mat,i-0,j-3,AFTER_RND);   
            }  
          return GeneWise21_DC_SHADOW_MATRIX(mat,i - 0,j - 3,AFTER_RND) ;    
          }  
      case AFTER_CENTRAL :   
        /* source AFTER_CODON is a special */ 
        temp = cscore - ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))) - (0);  
        if( temp == GeneWise21_DC_SHADOW_SPECIAL(mat,i - 0,j - 10,AFTER_CODON) ) {  
          *reti = i - 0; 
          *retj = j - 10;    
          *retstate = AFTER_CODON;   
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_DC_SHADOW_SPECIAL(mat,i-0,j-10,AFTER_CODON);    
            }  
          return GeneWise21_DC_SHADOW_MATRIX(mat,i - 0,j - 10,AFTER_CODON) ;     
          }  
        /* source AFTER_CODON is a special */ 
        temp = cscore - ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))) - (0);  
        if( temp == GeneWise21_DC_SHADOW_SPECIAL(mat,i - 0,j - 9,AFTER_CODON) )  {  
          *reti = i - 0; 
          *retj = j - 9; 
          *retstate = AFTER_CODON;   
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_DC_SHADOW_SPECIAL(mat,i-0,j-9,AFTER_CODON);     
            }  
          return GeneWise21_DC_SHADOW_MATRIX(mat,i - 0,j - 9,AFTER_CODON) ;  
          }  
        /* source AFTER_CODON is a special */ 
        temp = cscore - ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))) - (0);  
        if( temp == GeneWise21_DC_SHADOW_SPECIAL(mat,i - 0,j - 8,AFTER_CODON) )  {  
          *reti = i - 0; 
          *retj = j - 8; 
          *retstate = AFTER_CODON;   
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_DC_SHADOW_SPECIAL(mat,i-0,j-8,AFTER_CODON);     
            }  
          return GeneWise21_DC_SHADOW_MATRIX(mat,i - 0,j - 8,AFTER_CODON) ;  
          }  
        /* source AFTER_CENTRAL is a special */ 
        temp = cscore - ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CENTRAL2CENTRAL])) - (0);  
        if( temp == GeneWise21_DC_SHADOW_SPECIAL(mat,i - 0,j - 1,AFTER_CENTRAL) )    {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = AFTER_CENTRAL; 
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_DC_SHADOW_SPECIAL(mat,i-0,j-1,AFTER_CENTRAL);   
            }  
          return GeneWise21_DC_SHADOW_MATRIX(mat,i - 0,j - 1,AFTER_CENTRAL) ;    
          }  
      case AFTER_PY_TRACT :  
        /* source AFTER_PY_TRACT is a special */ 
        temp = cscore - (mat->gp21->transition[GP21_PY2PY]) - (mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)]);     
        if( temp == GeneWise21_DC_SHADOW_SPECIAL(mat,i - 0,j - 1,AFTER_PY_TRACT) )   {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = AFTER_PY_TRACT;    
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_DC_SHADOW_SPECIAL(mat,i-0,j-1,AFTER_PY_TRACT);  
            }  
          return GeneWise21_DC_SHADOW_MATRIX(mat,i - 0,j - 1,AFTER_PY_TRACT) ;   
          }  
        /* source AFTER_CENTRAL is a special */ 
        temp = cscore - (mat->gp21->transition[GP21_CENTRAL2PY]) - (mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)]);    
        if( temp == GeneWise21_DC_SHADOW_SPECIAL(mat,i - 0,j - 1,AFTER_CENTRAL) )    {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = AFTER_CENTRAL; 
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_DC_SHADOW_SPECIAL(mat,i-0,j-1,AFTER_CENTRAL);   
            }  
          return GeneWise21_DC_SHADOW_MATRIX(mat,i - 0,j - 1,AFTER_CENTRAL) ;    
          }  
      case AFTER_SPACER :    
        /* source AFTER_SPACER is a special */ 
        temp = cscore - (mat->gp21->transition[GP21_SPACER2SPACER]) - (mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)]);     
        if( temp == GeneWise21_DC_SHADOW_SPECIAL(mat,i - 0,j - 1,AFTER_SPACER) ) {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = AFTER_SPACER;  
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_DC_SHADOW_SPECIAL(mat,i-0,j-1,AFTER_SPACER);    
            }  
          return GeneWise21_DC_SHADOW_MATRIX(mat,i - 0,j - 1,AFTER_SPACER) ;     
          }  
        /* source AFTER_PY_TRACT is a special */ 
        temp = cscore - (mat->gp21->transition[GP21_PY2SPACER]) - (mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)]);     
        if( temp == GeneWise21_DC_SHADOW_SPECIAL(mat,i - 0,j - 1,AFTER_PY_TRACT) )   {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = AFTER_PY_TRACT;    
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_DC_SHADOW_SPECIAL(mat,i-0,j-1,AFTER_PY_TRACT);  
            }  
          return GeneWise21_DC_SHADOW_MATRIX(mat,i - 0,j - 1,AFTER_PY_TRACT) ;   
          }  
      default:   
        warn("Major problem (!) - in GeneWise21 special strip read off, position %d,%d state %d no source found  dropped into default on source switch!",i,j,state); 
        return (-1); 
      } /* end of switch on special states */ 
}    


/* Function:  max_matrix_to_special_GeneWise21(mat,i,j,state,cscore,reti,retj,retstate,retspecial,cellscore)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [GeneWise21 *]
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
int max_matrix_to_special_GeneWise21(GeneWise21 * mat,int i,int j,int state,int cscore,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore) 
{
    int temp;    
    *reti = (*retj) = (*retstate) = GeneWise21_READ_OFF_ERROR;   


    if( j < 0 || j > mat->lenj)  {  
      warn("In GeneWise21 matrix to special read off - out of bounds on matrix [j is %d in special]",j); 
      return -1; 
      }  


    switch(state)    { /*Switch state */ 
      case MATCH :   
        /* Source MATCH is not a special, should not get here! */ 
        /* Source MATCH is not a special, should not get here! */ 
        /* Source MATCH is not a special, should not get here! */ 
        /* Source MATCH is not a special, should not get here! */ 
        /* Source PY_2M is not a special, should not get here! */ 
        /* Source SPACER_2M is not a special, should not get here! */ 
        /* Source PY_1M is not a special, should not get here! */ 
        /* Source SPACER_1M is not a special, should not get here! */ 
        /* Source PY_0M is not a special, should not get here! */ 
        /* Source SPACER_0M is not a special, should not get here! */ 
        temp = cscore - ((mat->query->seg[i]->transition[GW_START2MATCH]+mat->query->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)])) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j));     
        if( temp == GeneWise21_DC_SHADOW_SPECIAL(mat,i - 1,j - 3,BEFORE_CODON) ) {  
          *reti = i - 1; 
          *retj = j - 3; 
          *retstate = BEFORE_CODON;  
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_DC_SHADOW_SPECIAL(mat,i-1,j-3,BEFORE_CODON);    
            }  
          return GeneWise21_DC_SHADOW_MATRIX(mat,i - 1,j - 3,BEFORE_CODON) ;     
          }  
        /* Source DELETE is not a special, should not get here! */ 
        /* Source INSERT is not a special, should not get here! */ 
        /* Source MATCH is not a special, should not get here! */ 
        warn("Major problem (!) - in GeneWise21 matrix to special read off, position %d,%d state %d no source found!",i,j,state);    
        return (-1); 
      case INSERT :  
        /* Source INSERT is not a special, should not get here! */ 
        /* Source INSERT is not a special, should not get here! */ 
        /* Source INSERT is not a special, should not get here! */ 
        /* Source INSERT is not a special, should not get here! */ 
        /* Source PY_2I is not a special, should not get here! */ 
        /* Source PY_1I is not a special, should not get here! */ 
        /* Source PY_0I is not a special, should not get here! */ 
        /* Source SPACER_2I is not a special, should not get here! */ 
        /* Source SPACER_1I is not a special, should not get here! */ 
        /* Source SPACER_0I is not a special, should not get here! */ 
        temp = cscore - ((mat->query->seg[i]->transition[GW_START2INSERT]+mat->query->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)])) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j));   
        if( temp == GeneWise21_DC_SHADOW_SPECIAL(mat,i - 0,j - 3,BEFORE_CODON) ) {  
          *reti = i - 0; 
          *retj = j - 3; 
          *retstate = BEFORE_CODON;  
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_DC_SHADOW_SPECIAL(mat,i-0,j-3,BEFORE_CODON);    
            }  
          return GeneWise21_DC_SHADOW_MATRIX(mat,i - 0,j - 3,BEFORE_CODON) ;     
          }  
        /* Source DELETE is not a special, should not get here! */ 
        /* Source INSERT is not a special, should not get here! */ 
        /* Source MATCH is not a special, should not get here! */ 
        warn("Major problem (!) - in GeneWise21 matrix to special read off, position %d,%d state %d no source found!",i,j,state);    
        return (-1); 
      case DELETE :  
        temp = cscore - (mat->query->seg[i]->transition[GW_START2DELETE]) -  (0);    
        if( temp == GeneWise21_DC_SHADOW_SPECIAL(mat,i - 1,j - 0,BEFORE_CODON) ) {  
          *reti = i - 1; 
          *retj = j - 0; 
          *retstate = BEFORE_CODON;  
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - GeneWise21_DC_SHADOW_SPECIAL(mat,i-1,j-0,BEFORE_CODON);    
            }  
          return GeneWise21_DC_SHADOW_MATRIX(mat,i - 1,j - 0,BEFORE_CODON) ;     
          }  
        /* Source DELETE is not a special, should not get here! */ 
        /* Source INSERT is not a special, should not get here! */ 
        /* Source MATCH is not a special, should not get here! */ 
        warn("Major problem (!) - in GeneWise21 matrix to special read off, position %d,%d state %d no source found!",i,j,state);    
        return (-1); 
      case CENTRAL_0M :  
        /* Source CENTRAL_0M is not a special, should not get here! */ 
        /* Source DELETE is not a special, should not get here! */ 
        /* Source INSERT is not a special, should not get here! */ 
        /* Source MATCH is not a special, should not get here! */ 
        warn("Major problem (!) - in GeneWise21 matrix to special read off, position %d,%d state %d no source found!",i,j,state);    
        return (-1); 
      case PY_0M :   
        /* Source PY_0M is not a special, should not get here! */ 
        /* Source CENTRAL_0M is not a special, should not get here! */ 
        warn("Major problem (!) - in GeneWise21 matrix to special read off, position %d,%d state %d no source found!",i,j,state);    
        return (-1); 
      case SPACER_0M :   
        /* Source SPACER_0M is not a special, should not get here! */ 
        /* Source PY_0M is not a special, should not get here! */ 
        warn("Major problem (!) - in GeneWise21 matrix to special read off, position %d,%d state %d no source found!",i,j,state);    
        return (-1); 
      case CENTRAL_1M :  
        /* Source CENTRAL_1M is not a special, should not get here! */ 
        /* Source DELETE is not a special, should not get here! */ 
        /* Source INSERT is not a special, should not get here! */ 
        /* Source MATCH is not a special, should not get here! */ 
        warn("Major problem (!) - in GeneWise21 matrix to special read off, position %d,%d state %d no source found!",i,j,state);    
        return (-1); 
      case PY_1M :   
        /* Source PY_1M is not a special, should not get here! */ 
        /* Source CENTRAL_1M is not a special, should not get here! */ 
        warn("Major problem (!) - in GeneWise21 matrix to special read off, position %d,%d state %d no source found!",i,j,state);    
        return (-1); 
      case SPACER_1M :   
        /* Source SPACER_1M is not a special, should not get here! */ 
        /* Source PY_1M is not a special, should not get here! */ 
        warn("Major problem (!) - in GeneWise21 matrix to special read off, position %d,%d state %d no source found!",i,j,state);    
        return (-1); 
      case CENTRAL_2M :  
        /* Source CENTRAL_2M is not a special, should not get here! */ 
        /* Source DELETE is not a special, should not get here! */ 
        /* Source INSERT is not a special, should not get here! */ 
        /* Source MATCH is not a special, should not get here! */ 
        warn("Major problem (!) - in GeneWise21 matrix to special read off, position %d,%d state %d no source found!",i,j,state);    
        return (-1); 
      case PY_2M :   
        /* Source PY_2M is not a special, should not get here! */ 
        /* Source CENTRAL_2M is not a special, should not get here! */ 
        warn("Major problem (!) - in GeneWise21 matrix to special read off, position %d,%d state %d no source found!",i,j,state);    
        return (-1); 
      case SPACER_2M :   
        /* Source SPACER_2M is not a special, should not get here! */ 
        /* Source PY_2M is not a special, should not get here! */ 
        warn("Major problem (!) - in GeneWise21 matrix to special read off, position %d,%d state %d no source found!",i,j,state);    
        return (-1); 
      case CENTRAL_0I :  
        /* Source CENTRAL_0I is not a special, should not get here! */ 
        /* Source DELETE is not a special, should not get here! */ 
        /* Source INSERT is not a special, should not get here! */ 
        /* Source MATCH is not a special, should not get here! */ 
        warn("Major problem (!) - in GeneWise21 matrix to special read off, position %d,%d state %d no source found!",i,j,state);    
        return (-1); 
      case PY_0I :   
        /* Source PY_0I is not a special, should not get here! */ 
        /* Source CENTRAL_0I is not a special, should not get here! */ 
        warn("Major problem (!) - in GeneWise21 matrix to special read off, position %d,%d state %d no source found!",i,j,state);    
        return (-1); 
      case SPACER_0I :   
        /* Source SPACER_0I is not a special, should not get here! */ 
        /* Source PY_0I is not a special, should not get here! */ 
        warn("Major problem (!) - in GeneWise21 matrix to special read off, position %d,%d state %d no source found!",i,j,state);    
        return (-1); 
      case CENTRAL_1I :  
        /* Source CENTRAL_1I is not a special, should not get here! */ 
        /* Source DELETE is not a special, should not get here! */ 
        /* Source INSERT is not a special, should not get here! */ 
        /* Source MATCH is not a special, should not get here! */ 
        warn("Major problem (!) - in GeneWise21 matrix to special read off, position %d,%d state %d no source found!",i,j,state);    
        return (-1); 
      case PY_1I :   
        /* Source PY_1I is not a special, should not get here! */ 
        /* Source CENTRAL_1I is not a special, should not get here! */ 
        warn("Major problem (!) - in GeneWise21 matrix to special read off, position %d,%d state %d no source found!",i,j,state);    
        return (-1); 
      case SPACER_1I :   
        /* Source SPACER_1I is not a special, should not get here! */ 
        /* Source PY_1I is not a special, should not get here! */ 
        warn("Major problem (!) - in GeneWise21 matrix to special read off, position %d,%d state %d no source found!",i,j,state);    
        return (-1); 
      case CENTRAL_2I :  
        /* Source CENTRAL_2I is not a special, should not get here! */ 
        /* Source INSERT is not a special, should not get here! */ 
        /* Source INSERT is not a special, should not get here! */ 
        /* Source MATCH is not a special, should not get here! */ 
        warn("Major problem (!) - in GeneWise21 matrix to special read off, position %d,%d state %d no source found!",i,j,state);    
        return (-1); 
      case PY_2I :   
        /* Source PY_2I is not a special, should not get here! */ 
        /* Source CENTRAL_2I is not a special, should not get here! */ 
        warn("Major problem (!) - in GeneWise21 matrix to special read off, position %d,%d state %d no source found!",i,j,state);    
        return (-1); 
      case SPACER_2I :   
        /* Source SPACER_2I is not a special, should not get here! */ 
        /* Source PY_2I is not a special, should not get here! */ 
        warn("Major problem (!) - in GeneWise21 matrix to special read off, position %d,%d state %d no source found!",i,j,state);    
        return (-1); 
      default:   
        warn("Major problem (!) - in GeneWise21 read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      } /* end of Switch state  */ 


}    


/* Function:  calculate_hidden_GeneWise21(mat,starti,startj,startstate,stopi,stopj,stopstate,dpenv)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [GeneWise21 *]
 * Arg:            starti [UNKN ] Undocumented argument [int]
 * Arg:            startj [UNKN ] Undocumented argument [int]
 * Arg:        startstate [UNKN ] Undocumented argument [int]
 * Arg:             stopi [UNKN ] Undocumented argument [int]
 * Arg:             stopj [UNKN ] Undocumented argument [int]
 * Arg:         stopstate [UNKN ] Undocumented argument [int]
 * Arg:             dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 */
void calculate_hidden_GeneWise21(GeneWise21 * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv) 
{
    register int i;  
    register int j;  
    register int score;  
    register int temp;   
    register int hiddenj;    


    hiddenj = startj;    


    init_hidden_GeneWise21(mat,starti,startj,stopi,stopj);   


    GeneWise21_HIDDEN_MATRIX(mat,starti,startj,startstate) = 0;  


    for(j=startj;j<=stopj;j++)   {  
      for(i=starti;i<=stopi;i++) {  
        /* Should *not* do very first cell as this is the one set to zero in one state! */ 
        if( i == starti && j == startj ) 
          continue;  
        if( dpenv != NULL && is_in_DPEnvelope(dpenv,i,j) == FALSE )  { /*Is not in envelope*/ 
          GeneWise21_HIDDEN_MATRIX(mat,i,j,MATCH) = NEGI;    
          GeneWise21_HIDDEN_MATRIX(mat,i,j,INSERT) = NEGI;   
          GeneWise21_HIDDEN_MATRIX(mat,i,j,DELETE) = NEGI;   
          GeneWise21_HIDDEN_MATRIX(mat,i,j,CENTRAL_0M) = NEGI;   
          GeneWise21_HIDDEN_MATRIX(mat,i,j,PY_0M) = NEGI;    
          GeneWise21_HIDDEN_MATRIX(mat,i,j,SPACER_0M) = NEGI;    
          GeneWise21_HIDDEN_MATRIX(mat,i,j,CENTRAL_1M) = NEGI;   
          GeneWise21_HIDDEN_MATRIX(mat,i,j,PY_1M) = NEGI;    
          GeneWise21_HIDDEN_MATRIX(mat,i,j,SPACER_1M) = NEGI;    
          GeneWise21_HIDDEN_MATRIX(mat,i,j,CENTRAL_2M) = NEGI;   
          GeneWise21_HIDDEN_MATRIX(mat,i,j,PY_2M) = NEGI;    
          GeneWise21_HIDDEN_MATRIX(mat,i,j,SPACER_2M) = NEGI;    
          GeneWise21_HIDDEN_MATRIX(mat,i,j,CENTRAL_0I) = NEGI;   
          GeneWise21_HIDDEN_MATRIX(mat,i,j,PY_0I) = NEGI;    
          GeneWise21_HIDDEN_MATRIX(mat,i,j,SPACER_0I) = NEGI;    
          GeneWise21_HIDDEN_MATRIX(mat,i,j,CENTRAL_1I) = NEGI;   
          GeneWise21_HIDDEN_MATRIX(mat,i,j,PY_1I) = NEGI;    
          GeneWise21_HIDDEN_MATRIX(mat,i,j,SPACER_1I) = NEGI;    
          GeneWise21_HIDDEN_MATRIX(mat,i,j,CENTRAL_2I) = NEGI;   
          GeneWise21_HIDDEN_MATRIX(mat,i,j,PY_2I) = NEGI;    
          GeneWise21_HIDDEN_MATRIX(mat,i,j,SPACER_2I) = NEGI;    
          continue;  
          } /* end of Is not in envelope */ 


        /* For state MATCH */ 
        /* setting first movement to score */ 
        score = GeneWise21_HIDDEN_MATRIX(mat,i-1,j-3,MATCH) + (mat->query->seg[i]->transition[GW_MATCH2MATCH]+mat->query->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]);     
        /* From state INSERT to state MATCH */ 
        temp = GeneWise21_HIDDEN_MATRIX(mat,i-1,j-3,INSERT) + (mat->query->seg[i]->transition[GW_INSERT2MATCH]+mat->query->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]);    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state MATCH */ 
        temp = GeneWise21_HIDDEN_MATRIX(mat,i-1,j-3,DELETE) + (mat->query->seg[i]->transition[GW_DELETE2MATCH]+mat->query->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]);    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SPACER_0M to state MATCH */ 
        temp = GeneWise21_HIDDEN_MATRIX(mat,i-1,j-6,SPACER_0M) + (((mat->gp21->transition[GP21_SPACER2CDS]+mat->query->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)])+CSEQ_GENOMIC_3SS(mat->target,(j-3)))+mat->query->seg[i]->transition[GW_MATCH_BALANCE_3SS]);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state PY_0M to state MATCH */ 
        temp = GeneWise21_HIDDEN_MATRIX(mat,i-1,j-6,PY_0M) + (((mat->gp21->transition[GP21_PY2CDS]+mat->query->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)])+CSEQ_GENOMIC_3SS(mat->target,(j-3)))+mat->query->seg[i]->transition[GW_MATCH_BALANCE_3SS]);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SPACER_1M to state MATCH */ 
        temp = GeneWise21_HIDDEN_MATRIX(mat,i-1,j-5,SPACER_1M) + (mat->gp21->transition[GP21_SPACER2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-2)));   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state PY_1M to state MATCH */ 
        temp = GeneWise21_HIDDEN_MATRIX(mat,i-1,j-5,PY_1M) + (mat->gp21->transition[GP21_PY2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-2)));   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SPACER_2M to state MATCH */ 
        temp = GeneWise21_HIDDEN_MATRIX(mat,i-1,j-4,SPACER_2M) + (mat->gp21->transition[GP21_SPACER2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-1)));   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state PY_2M to state MATCH */ 
        temp = GeneWise21_HIDDEN_MATRIX(mat,i-1,j-4,PY_2M) + (mat->gp21->transition[GP21_PY2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-1)));   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = GeneWise21_HIDDEN_MATRIX(mat,i-1,j-2,MATCH) + mat->gp21->transition[GP21_DELETE_1_BASE];  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = GeneWise21_HIDDEN_MATRIX(mat,i-1,j-1,MATCH) + mat->gp21->transition[GP21_DELETE_2_BASE];  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = GeneWise21_HIDDEN_MATRIX(mat,i-1,j-4,MATCH) + mat->gp21->transition[GP21_INSERT_1_BASE];  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = GeneWise21_HIDDEN_MATRIX(mat,i-1,j-5,MATCH) + mat->gp21->transition[GP21_INSERT_2_BASE];  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for MATCH */ 
        /* Add any movement independant score and put away */ 
         score += CSEQ_GENOMIC_CDSPOT(mat->target,j);    
         GeneWise21_HIDDEN_MATRIX(mat,i,j,MATCH) = score;    
        /* Finished calculating state MATCH */ 


        /* For state INSERT */ 
        /* setting first movement to score */ 
        score = GeneWise21_HIDDEN_MATRIX(mat,i-0,j-3,MATCH) + (mat->query->seg[i]->transition[GW_MATCH2INSERT]+mat->query->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)]);   
        /* From state INSERT to state INSERT */ 
        temp = GeneWise21_HIDDEN_MATRIX(mat,i-0,j-3,INSERT) + (mat->query->seg[i]->transition[GW_INSERT2INSERT]+mat->query->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)]);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state INSERT */ 
        temp = GeneWise21_HIDDEN_MATRIX(mat,i-0,j-3,DELETE) + (mat->query->seg[i]->transition[GW_DELETE2INSERT]+mat->query->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)]);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SPACER_0I to state INSERT */ 
        temp = GeneWise21_HIDDEN_MATRIX(mat,i-1,j-6,SPACER_0I) + (((mat->gp21->transition[GP21_SPACER2CDS]+mat->query->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)])+CSEQ_GENOMIC_3SS(mat->target,(j-3)))+mat->query->seg[i]->transition[GW_MATCH_BALANCE_3SS]);    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SPACER_1I to state INSERT */ 
        temp = GeneWise21_HIDDEN_MATRIX(mat,i-1,j-5,SPACER_1I) + (mat->gp21->transition[GP21_SPACER2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-2)));   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SPACER_2I to state INSERT */ 
        temp = GeneWise21_HIDDEN_MATRIX(mat,i-1,j-4,SPACER_2I) + (mat->gp21->transition[GP21_SPACER2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-1)));   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state PY_0I to state INSERT */ 
        temp = GeneWise21_HIDDEN_MATRIX(mat,i-1,j-6,PY_0I) + (((mat->gp21->transition[GP21_PY2CDS]+mat->query->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)])+CSEQ_GENOMIC_3SS(mat->target,(j-3)))+mat->query->seg[i]->transition[GW_MATCH_BALANCE_3SS]);    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state PY_1I to state INSERT */ 
        temp = GeneWise21_HIDDEN_MATRIX(mat,i-1,j-5,PY_1I) + (mat->gp21->transition[GP21_PY2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-2)));   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state PY_2I to state INSERT */ 
        temp = GeneWise21_HIDDEN_MATRIX(mat,i-1,j-4,PY_2I) + (mat->gp21->transition[GP21_PY2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-1)));   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state INSERT */ 
        temp = GeneWise21_HIDDEN_MATRIX(mat,i-1,j-2,INSERT) + mat->gp21->transition[GP21_DELETE_1_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state INSERT */ 
        temp = GeneWise21_HIDDEN_MATRIX(mat,i-1,j-1,INSERT) + mat->gp21->transition[GP21_DELETE_2_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state INSERT */ 
        temp = GeneWise21_HIDDEN_MATRIX(mat,i-1,j-4,INSERT) + mat->gp21->transition[GP21_INSERT_1_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state INSERT */ 
        temp = GeneWise21_HIDDEN_MATRIX(mat,i-1,j-5,INSERT) + mat->gp21->transition[GP21_INSERT_2_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for INSERT */ 
        /* Add any movement independant score and put away */ 
         score += CSEQ_GENOMIC_CDSPOT(mat->target,j);    
         GeneWise21_HIDDEN_MATRIX(mat,i,j,INSERT) = score;   
        /* Finished calculating state INSERT */ 


        /* For state DELETE */ 
        /* setting first movement to score */ 
        score = GeneWise21_HIDDEN_MATRIX(mat,i-1,j-0,MATCH) + mat->query->seg[i]->transition[GW_MATCH2DELETE];   
        /* From state INSERT to state DELETE */ 
        temp = GeneWise21_HIDDEN_MATRIX(mat,i-1,j-0,INSERT) + mat->query->seg[i]->transition[GW_INSERT2DELETE];  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state DELETE */ 
        temp = GeneWise21_HIDDEN_MATRIX(mat,i-1,j-0,DELETE) + mat->query->seg[i]->transition[GW_DELETE2DELETE];  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for DELETE */ 
        /* Add any movement independant score and put away */ 
         GeneWise21_HIDDEN_MATRIX(mat,i,j,DELETE) = score;   
        /* Finished calculating state DELETE */ 


        /* For state CENTRAL_0M */ 
        /* setting first movement to score */ 
        score = GeneWise21_HIDDEN_MATRIX(mat,i-0,j-8,MATCH) + (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_MATCH2MATCH])+mat->query->seg[i]->transition[GW_MATCH_BALANCE_5SS]);    
        /* From state INSERT to state CENTRAL_0M */ 
        temp = GeneWise21_HIDDEN_MATRIX(mat,i-0,j-8,INSERT) + (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_INSERT2MATCH])+mat->query->seg[i]->transition[GW_INSERT_BALANCE_5SS]);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state CENTRAL_0M */ 
        temp = GeneWise21_HIDDEN_MATRIX(mat,i-0,j-8,DELETE) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_DELETE2MATCH]);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state CENTRAL_0M to state CENTRAL_0M */ 
        temp = GeneWise21_HIDDEN_MATRIX(mat,i-0,j-1,CENTRAL_0M) + (mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CENTRAL2CENTRAL]);    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for CENTRAL_0M */ 
        /* Add any movement independant score and put away */ 
         GeneWise21_HIDDEN_MATRIX(mat,i,j,CENTRAL_0M) = score;   
        /* Finished calculating state CENTRAL_0M */ 


        /* For state PY_0M */ 
        /* setting first movement to score */ 
        score = GeneWise21_HIDDEN_MATRIX(mat,i-0,j-1,CENTRAL_0M) + mat->gp21->transition[GP21_CENTRAL2PY];   
        /* From state PY_0M to state PY_0M */ 
        temp = GeneWise21_HIDDEN_MATRIX(mat,i-0,j-1,PY_0M) + mat->gp21->transition[GP21_PY2PY];  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for PY_0M */ 
        /* Add any movement independant score and put away */ 
         score += mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)];   
         GeneWise21_HIDDEN_MATRIX(mat,i,j,PY_0M) = score;    
        /* Finished calculating state PY_0M */ 


        /* For state SPACER_0M */ 
        /* setting first movement to score */ 
        score = GeneWise21_HIDDEN_MATRIX(mat,i-0,j-1,PY_0M) + mat->gp21->transition[GP21_PY2SPACER];     
        /* From state SPACER_0M to state SPACER_0M */ 
        temp = GeneWise21_HIDDEN_MATRIX(mat,i-0,j-1,SPACER_0M) + mat->gp21->transition[GP21_SPACER2SPACER];  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for SPACER_0M */ 
        /* Add any movement independant score and put away */ 
         score += mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)];   
         GeneWise21_HIDDEN_MATRIX(mat,i,j,SPACER_0M) = score;    
        /* Finished calculating state SPACER_0M */ 


        /* For state CENTRAL_1M */ 
        /* setting first movement to score */ 
        score = GeneWise21_HIDDEN_MATRIX(mat,i-0,j-9,MATCH) + (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_MATCH2MATCH])+mat->query->seg[i]->transition[GW_MATCH_BALANCE_5SS]);    
        /* From state INSERT to state CENTRAL_1M */ 
        temp = GeneWise21_HIDDEN_MATRIX(mat,i-0,j-9,INSERT) + (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_INSERT2MATCH])+mat->query->seg[i]->transition[GW_MATCH_BALANCE_5SS]);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state CENTRAL_1M */ 
        temp = GeneWise21_HIDDEN_MATRIX(mat,i-0,j-9,DELETE) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_DELETE2MATCH]);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state CENTRAL_1M to state CENTRAL_1M */ 
        temp = GeneWise21_HIDDEN_MATRIX(mat,i-0,j-1,CENTRAL_1M) + (mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CENTRAL2CENTRAL]);    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for CENTRAL_1M */ 
        /* Add any movement independant score and put away */ 
         GeneWise21_HIDDEN_MATRIX(mat,i,j,CENTRAL_1M) = score;   
        /* Finished calculating state CENTRAL_1M */ 


        /* For state PY_1M */ 
        /* setting first movement to score */ 
        score = GeneWise21_HIDDEN_MATRIX(mat,i-0,j-1,CENTRAL_1M) + mat->gp21->transition[GP21_CENTRAL2PY];   
        /* From state PY_1M to state PY_1M */ 
        temp = GeneWise21_HIDDEN_MATRIX(mat,i-0,j-1,PY_1M) + mat->gp21->transition[GP21_PY2PY];  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for PY_1M */ 
        /* Add any movement independant score and put away */ 
         score += mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)];   
         GeneWise21_HIDDEN_MATRIX(mat,i,j,PY_1M) = score;    
        /* Finished calculating state PY_1M */ 


        /* For state SPACER_1M */ 
        /* setting first movement to score */ 
        score = GeneWise21_HIDDEN_MATRIX(mat,i-0,j-1,PY_1M) + mat->gp21->transition[GP21_PY2SPACER];     
        /* From state SPACER_1M to state SPACER_1M */ 
        temp = GeneWise21_HIDDEN_MATRIX(mat,i-0,j-1,SPACER_1M) + mat->gp21->transition[GP21_SPACER2SPACER];  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for SPACER_1M */ 
        /* Add any movement independant score and put away */ 
         score += mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)];   
         GeneWise21_HIDDEN_MATRIX(mat,i,j,SPACER_1M) = score;    
        /* Finished calculating state SPACER_1M */ 


        /* For state CENTRAL_2M */ 
        /* setting first movement to score */ 
        score = GeneWise21_HIDDEN_MATRIX(mat,i-0,j-10,MATCH) + (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_MATCH2MATCH])+mat->query->seg[i]->transition[GW_MATCH_BALANCE_5SS]);   
        /* From state INSERT to state CENTRAL_2M */ 
        temp = GeneWise21_HIDDEN_MATRIX(mat,i-0,j-10,INSERT) + (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_INSERT2MATCH])+mat->query->seg[i]->transition[GW_MATCH_BALANCE_5SS]);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state CENTRAL_2M */ 
        temp = GeneWise21_HIDDEN_MATRIX(mat,i-0,j-10,DELETE) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_DELETE2MATCH]);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state CENTRAL_2M to state CENTRAL_2M */ 
        temp = GeneWise21_HIDDEN_MATRIX(mat,i-0,j-1,CENTRAL_2M) + (mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CENTRAL2CENTRAL]);    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for CENTRAL_2M */ 
        /* Add any movement independant score and put away */ 
         GeneWise21_HIDDEN_MATRIX(mat,i,j,CENTRAL_2M) = score;   
        /* Finished calculating state CENTRAL_2M */ 


        /* For state PY_2M */ 
        /* setting first movement to score */ 
        score = GeneWise21_HIDDEN_MATRIX(mat,i-0,j-1,CENTRAL_2M) + mat->gp21->transition[GP21_CENTRAL2PY];   
        /* From state PY_2M to state PY_2M */ 
        temp = GeneWise21_HIDDEN_MATRIX(mat,i-0,j-1,PY_2M) + mat->gp21->transition[GP21_PY2PY];  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for PY_2M */ 
        /* Add any movement independant score and put away */ 
         score += mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)];   
         GeneWise21_HIDDEN_MATRIX(mat,i,j,PY_2M) = score;    
        /* Finished calculating state PY_2M */ 


        /* For state SPACER_2M */ 
        /* setting first movement to score */ 
        score = GeneWise21_HIDDEN_MATRIX(mat,i-0,j-1,PY_2M) + mat->gp21->transition[GP21_PY2SPACER];     
        /* From state SPACER_2M to state SPACER_2M */ 
        temp = GeneWise21_HIDDEN_MATRIX(mat,i-0,j-1,SPACER_2M) + mat->gp21->transition[GP21_SPACER2SPACER];  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for SPACER_2M */ 
        /* Add any movement independant score and put away */ 
         score += mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)];   
         GeneWise21_HIDDEN_MATRIX(mat,i,j,SPACER_2M) = score;    
        /* Finished calculating state SPACER_2M */ 


        /* For state CENTRAL_0I */ 
        /* setting first movement to score */ 
        score = GeneWise21_HIDDEN_MATRIX(mat,i-0,j-8,MATCH) + (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_MATCH2INSERT])+mat->query->seg[i]->transition[GW_MATCH_BALANCE_5SS]);   
        /* From state INSERT to state CENTRAL_0I */ 
        temp = GeneWise21_HIDDEN_MATRIX(mat,i-0,j-8,INSERT) + (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_INSERT2INSERT])+mat->query->seg[i]->transition[GW_INSERT_BALANCE_5SS]);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state CENTRAL_0I */ 
        temp = GeneWise21_HIDDEN_MATRIX(mat,i-0,j-8,DELETE) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_DELETE2INSERT]);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state CENTRAL_0I to state CENTRAL_0I */ 
        temp = GeneWise21_HIDDEN_MATRIX(mat,i-0,j-1,CENTRAL_0I) + (mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CENTRAL2CENTRAL]);    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for CENTRAL_0I */ 
        /* Add any movement independant score and put away */ 
         GeneWise21_HIDDEN_MATRIX(mat,i,j,CENTRAL_0I) = score;   
        /* Finished calculating state CENTRAL_0I */ 


        /* For state PY_0I */ 
        /* setting first movement to score */ 
        score = GeneWise21_HIDDEN_MATRIX(mat,i-0,j-1,CENTRAL_0I) + mat->gp21->transition[GP21_CENTRAL2PY];   
        /* From state PY_0I to state PY_0I */ 
        temp = GeneWise21_HIDDEN_MATRIX(mat,i-0,j-1,PY_0I) + mat->gp21->transition[GP21_PY2PY];  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for PY_0I */ 
        /* Add any movement independant score and put away */ 
         score += mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)];   
         GeneWise21_HIDDEN_MATRIX(mat,i,j,PY_0I) = score;    
        /* Finished calculating state PY_0I */ 


        /* For state SPACER_0I */ 
        /* setting first movement to score */ 
        score = GeneWise21_HIDDEN_MATRIX(mat,i-0,j-1,PY_0I) + mat->gp21->transition[GP21_PY2SPACER];     
        /* From state SPACER_0I to state SPACER_0I */ 
        temp = GeneWise21_HIDDEN_MATRIX(mat,i-0,j-1,SPACER_0I) + mat->gp21->transition[GP21_SPACER2SPACER];  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for SPACER_0I */ 
        /* Add any movement independant score and put away */ 
         score += mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)];   
         GeneWise21_HIDDEN_MATRIX(mat,i,j,SPACER_0I) = score;    
        /* Finished calculating state SPACER_0I */ 


        /* For state CENTRAL_1I */ 
        /* setting first movement to score */ 
        score = GeneWise21_HIDDEN_MATRIX(mat,i-0,j-9,MATCH) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_MATCH2INSERT]);  
        /* From state INSERT to state CENTRAL_1I */ 
        temp = GeneWise21_HIDDEN_MATRIX(mat,i-0,j-9,INSERT) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_INSERT2INSERT]);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state CENTRAL_1I */ 
        temp = GeneWise21_HIDDEN_MATRIX(mat,i-0,j-9,DELETE) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_DELETE2INSERT]);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state CENTRAL_1I to state CENTRAL_1I */ 
        temp = GeneWise21_HIDDEN_MATRIX(mat,i-0,j-1,CENTRAL_1I) + (mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CENTRAL2CENTRAL]);    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for CENTRAL_1I */ 
        /* Add any movement independant score and put away */ 
         GeneWise21_HIDDEN_MATRIX(mat,i,j,CENTRAL_1I) = score;   
        /* Finished calculating state CENTRAL_1I */ 


        /* For state PY_1I */ 
        /* setting first movement to score */ 
        score = GeneWise21_HIDDEN_MATRIX(mat,i-0,j-1,CENTRAL_1I) + mat->gp21->transition[GP21_CENTRAL2PY];   
        /* From state PY_1I to state PY_1I */ 
        temp = GeneWise21_HIDDEN_MATRIX(mat,i-0,j-1,PY_1I) + mat->gp21->transition[GP21_PY2PY];  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for PY_1I */ 
        /* Add any movement independant score and put away */ 
         score += mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)];   
         GeneWise21_HIDDEN_MATRIX(mat,i,j,PY_1I) = score;    
        /* Finished calculating state PY_1I */ 


        /* For state SPACER_1I */ 
        /* setting first movement to score */ 
        score = GeneWise21_HIDDEN_MATRIX(mat,i-0,j-1,PY_1I) + mat->gp21->transition[GP21_PY2SPACER];     
        /* From state SPACER_1I to state SPACER_1I */ 
        temp = GeneWise21_HIDDEN_MATRIX(mat,i-0,j-1,SPACER_1I) + mat->gp21->transition[GP21_SPACER2SPACER];  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for SPACER_1I */ 
        /* Add any movement independant score and put away */ 
         score += mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)];   
         GeneWise21_HIDDEN_MATRIX(mat,i,j,SPACER_1I) = score;    
        /* Finished calculating state SPACER_1I */ 


        /* For state CENTRAL_2I */ 
        /* setting first movement to score */ 
        score = GeneWise21_HIDDEN_MATRIX(mat,i-0,j-10,MATCH) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_MATCH2INSERT]);     
        /* From state INSERT to state CENTRAL_2I */ 
        temp = GeneWise21_HIDDEN_MATRIX(mat,i-0,j-10,INSERT) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_INSERT2INSERT]);    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state CENTRAL_2I */ 
        temp = GeneWise21_HIDDEN_MATRIX(mat,i-0,j-10,INSERT) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_DELETE2INSERT]);    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state CENTRAL_2I to state CENTRAL_2I */ 
        temp = GeneWise21_HIDDEN_MATRIX(mat,i-0,j-1,CENTRAL_2I) + (mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CENTRAL2CENTRAL]);    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for CENTRAL_2I */ 
        /* Add any movement independant score and put away */ 
         GeneWise21_HIDDEN_MATRIX(mat,i,j,CENTRAL_2I) = score;   
        /* Finished calculating state CENTRAL_2I */ 


        /* For state PY_2I */ 
        /* setting first movement to score */ 
        score = GeneWise21_HIDDEN_MATRIX(mat,i-0,j-1,CENTRAL_2I) + mat->gp21->transition[GP21_CENTRAL2PY];   
        /* From state PY_2I to state PY_2I */ 
        temp = GeneWise21_HIDDEN_MATRIX(mat,i-0,j-1,PY_2I) + mat->gp21->transition[GP21_PY2PY];  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for PY_2I */ 
        /* Add any movement independant score and put away */ 
         score += mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)];   
         GeneWise21_HIDDEN_MATRIX(mat,i,j,PY_2I) = score;    
        /* Finished calculating state PY_2I */ 


        /* For state SPACER_2I */ 
        /* setting first movement to score */ 
        score = GeneWise21_HIDDEN_MATRIX(mat,i-0,j-1,PY_2I) + mat->gp21->transition[GP21_PY2SPACER];     
        /* From state SPACER_2I to state SPACER_2I */ 
        temp = GeneWise21_HIDDEN_MATRIX(mat,i-0,j-1,SPACER_2I) + mat->gp21->transition[GP21_SPACER2SPACER];  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for SPACER_2I */ 
        /* Add any movement independant score and put away */ 
         score += mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)];   
         GeneWise21_HIDDEN_MATRIX(mat,i,j,SPACER_2I) = score;    
        /* Finished calculating state SPACER_2I */ 
        }  
      }  


    return;  
}    


/* Function:  init_hidden_GeneWise21(mat,starti,startj,stopi,stopj)
 *
 * Descrip: No Description
 *
 * Arg:           mat [UNKN ] Undocumented argument [GeneWise21 *]
 * Arg:        starti [UNKN ] Undocumented argument [int]
 * Arg:        startj [UNKN ] Undocumented argument [int]
 * Arg:         stopi [UNKN ] Undocumented argument [int]
 * Arg:         stopj [UNKN ] Undocumented argument [int]
 *
 */
void init_hidden_GeneWise21(GeneWise21 * mat,int starti,int startj,int stopi,int stopj) 
{
    register int i;  
    register int j;  
    register int hiddenj;    


    hiddenj = startj;    
    for(j=(startj-10);j<=stopj;j++)  {  
      for(i=(starti-1);i<=stopi;i++) {  
        GeneWise21_HIDDEN_MATRIX(mat,i,j,MATCH) = NEGI;
 
        GeneWise21_HIDDEN_MATRIX(mat,i,j,INSERT) = NEGI;
    
        GeneWise21_HIDDEN_MATRIX(mat,i,j,DELETE) = NEGI;
    
        GeneWise21_HIDDEN_MATRIX(mat,i,j,CENTRAL_0M) = NEGI;
    
        GeneWise21_HIDDEN_MATRIX(mat,i,j,PY_0M) = NEGI;
 
        GeneWise21_HIDDEN_MATRIX(mat,i,j,SPACER_0M) = NEGI;
 
        GeneWise21_HIDDEN_MATRIX(mat,i,j,CENTRAL_1M) = NEGI;
    
        GeneWise21_HIDDEN_MATRIX(mat,i,j,PY_1M) = NEGI;
 
        GeneWise21_HIDDEN_MATRIX(mat,i,j,SPACER_1M) = NEGI;
 
        GeneWise21_HIDDEN_MATRIX(mat,i,j,CENTRAL_2M) = NEGI;
    
        GeneWise21_HIDDEN_MATRIX(mat,i,j,PY_2M) = NEGI;
 
        GeneWise21_HIDDEN_MATRIX(mat,i,j,SPACER_2M) = NEGI;
 
        GeneWise21_HIDDEN_MATRIX(mat,i,j,CENTRAL_0I) = NEGI;
    
        GeneWise21_HIDDEN_MATRIX(mat,i,j,PY_0I) = NEGI;
 
        GeneWise21_HIDDEN_MATRIX(mat,i,j,SPACER_0I) = NEGI;
 
        GeneWise21_HIDDEN_MATRIX(mat,i,j,CENTRAL_1I) = NEGI;
    
        GeneWise21_HIDDEN_MATRIX(mat,i,j,PY_1I) = NEGI;
 
        GeneWise21_HIDDEN_MATRIX(mat,i,j,SPACER_1I) = NEGI;
 
        GeneWise21_HIDDEN_MATRIX(mat,i,j,CENTRAL_2I) = NEGI;
    
        GeneWise21_HIDDEN_MATRIX(mat,i,j,PY_2I) = NEGI;
 
        GeneWise21_HIDDEN_MATRIX(mat,i,j,SPACER_2I) = NEGI;
 
        }  
      }  


    return;  
}    


/* Function:  full_dc_GeneWise21(mat,starti,startj,startstate,stopi,stopj,stopstate,out,donej,totalj,dpenv)
 *
 * Descrip:    The main divide-and-conquor routine. Basically, call /PackAln_calculate_small_GeneWise21
 *             Not this function, which is pretty hard core. 
 *             Function is given start/end points (in main matrix) for alignment
 *             It does some checks, decides whether start/end in j is small enough for explicit calc
 *               - if yes, calculates it, reads off into PackAln (out), adds the j distance to donej and returns TRUE
 *               - if no,  uses /do_dc_single_pass_GeneWise21 to get mid-point
 *                          saves midpoint, and calls itself to do right portion then left portion
 *             right then left ensures PackAln is added the 'right' way, ie, back-to-front
 *             returns FALSE on any error, with a warning
 *
 *
 * Arg:               mat [UNKN ] Matrix with small memory implementation [GeneWise21 *]
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
boolean full_dc_GeneWise21(GeneWise21 * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out,int * donej,int totalj,DPEnvelope * dpenv) 
{
    int lstarti; 
    int lstartj; 
    int lstate;  


    if( mat->basematrix->type != BASEMATRIX_TYPE_SHADOW) {  
      warn("*Very* bad error! - non shadow matrix type in full_dc_GeneWise21");  
      return FALSE;  
      }  


    if( starti == -1 || startj == -1 || startstate == -1 || stopi == -1 || stopstate == -1)  {  
      warn("In full dc program, passed bad indices, indices passed were %d:%d[%d] to %d:%d[%d]\n",starti,startj,startstate,stopi,stopj,stopstate);   
      return FALSE;  
      }  


    if( stopj - startj < 50) {  
      log_full_error(REPORT,0,"[%d,%d][%d,%d] Explicit read off",starti,startj,stopi,stopj);/* Build hidden explicit matrix */ 
      calculate_hidden_GeneWise21(mat,starti,startj,startstate,stopi,stopj,stopstate,dpenv);     
      *donej += (stopj - startj);   /* Now read it off into out */ 
      if( read_hidden_GeneWise21(mat,starti,startj,startstate,stopi,stopj,stopstate,out) == FALSE)   {  
        warn("In full dc, at %d:%d,%d:%d got a bad hidden explicit read off... ",starti,startj,stopi,stopj); 
        return FALSE;    
        }  
      return TRUE;   
      }  


/* In actual divide and conquor */ 
    if( do_dc_single_pass_GeneWise21(mat,starti,startj,startstate,stopi,stopj,stopstate,dpenv,(int)(*donej*100)/totalj) == FALSE)    {  
      warn("In divide and conquor for GeneWise21, at bound %d:%d to %d:%d, unable to calculate midpoint. Problem!",starti,startj,stopi,stopj);   
      return FALSE;  
      }  


/* Ok... now we have to call on each side of the matrix */ 
/* We have to retrieve left hand side positions, as they will be vapped by the time we call LHS */ 
    lstarti= GeneWise21_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,0);    
    lstartj= GeneWise21_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,1);    
    lstate = GeneWise21_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,2);    


/* Call on right hand side: this lets us do the correct read off */ 
    if( full_dc_GeneWise21(mat,GeneWise21_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,3),GeneWise21_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,4),GeneWise21_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,5),stopi,stopj,stopstate,out,donej,totalj,dpenv) == FALSE)   {  
/* Warning already issued, simply chained back up to top */ 
      return FALSE;  
      }  
/* Call on left hand side */ 
    if( full_dc_GeneWise21(mat,starti,startj,startstate,lstarti,lstartj,lstate,out,donej,totalj,dpenv) == FALSE) {  
/* Warning already issued, simply chained back up to top */ 
      return FALSE;  
      }  


    return TRUE;     
}    


/* Function:  do_dc_single_pass_GeneWise21(mat,starti,startj,startstate,stopi,stopj,stopstate,dpenv,perc_done)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [GeneWise21 *]
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
boolean do_dc_single_pass_GeneWise21(GeneWise21 * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv,int perc_done) 
{
    int halfj;   
    halfj = startj + ((stopj - startj)/2);   


    init_dc_GeneWise21(mat); 


    GeneWise21_DC_SHADOW_MATRIX(mat,starti,startj,startstate) = 0;   
    run_up_dc_GeneWise21(mat,starti,stopi,startj,halfj-1,dpenv,perc_done);   
    push_dc_at_merge_GeneWise21(mat,starti,stopi,halfj,&halfj,dpenv);    
    follow_on_dc_GeneWise21(mat,starti,stopi,halfj,stopj,dpenv,perc_done);   
    return TRUE; 
}    


/* Function:  push_dc_at_merge_GeneWise21(mat,starti,stopi,startj,stopj,dpenv)
 *
 * Descrip: No Description
 *
 * Arg:           mat [UNKN ] Undocumented argument [GeneWise21 *]
 * Arg:        starti [UNKN ] Undocumented argument [int]
 * Arg:         stopi [UNKN ] Undocumented argument [int]
 * Arg:        startj [UNKN ] Undocumented argument [int]
 * Arg:         stopj [UNKN ] Undocumented argument [int *]
 * Arg:         dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 */
void push_dc_at_merge_GeneWise21(GeneWise21 * mat,int starti,int stopi,int startj,int * stopj,DPEnvelope * dpenv) 
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
          GeneWise21_DC_SHADOW_MATRIX(mat,i,j,MATCH) = NEGI;     
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,0) = (-100);  
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,1) = (-100);  
          GeneWise21_DC_SHADOW_MATRIX(mat,i,j,INSERT) = NEGI;    
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,0) = (-100); 
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,1) = (-100); 
          GeneWise21_DC_SHADOW_MATRIX(mat,i,j,DELETE) = NEGI;    
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,0) = (-100); 
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,1) = (-100); 
          GeneWise21_DC_SHADOW_MATRIX(mat,i,j,CENTRAL_0M) = NEGI;    
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_0M,0) = (-100); 
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_0M,1) = (-100); 
          GeneWise21_DC_SHADOW_MATRIX(mat,i,j,PY_0M) = NEGI;     
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_0M,0) = (-100);  
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_0M,1) = (-100);  
          GeneWise21_DC_SHADOW_MATRIX(mat,i,j,SPACER_0M) = NEGI;     
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_0M,0) = (-100);  
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_0M,1) = (-100);  
          GeneWise21_DC_SHADOW_MATRIX(mat,i,j,CENTRAL_1M) = NEGI;    
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_1M,0) = (-100); 
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_1M,1) = (-100); 
          GeneWise21_DC_SHADOW_MATRIX(mat,i,j,PY_1M) = NEGI;     
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_1M,0) = (-100);  
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_1M,1) = (-100);  
          GeneWise21_DC_SHADOW_MATRIX(mat,i,j,SPACER_1M) = NEGI;     
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_1M,0) = (-100);  
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_1M,1) = (-100);  
          GeneWise21_DC_SHADOW_MATRIX(mat,i,j,CENTRAL_2M) = NEGI;    
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_2M,0) = (-100); 
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_2M,1) = (-100); 
          GeneWise21_DC_SHADOW_MATRIX(mat,i,j,PY_2M) = NEGI;     
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_2M,0) = (-100);  
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_2M,1) = (-100);  
          GeneWise21_DC_SHADOW_MATRIX(mat,i,j,SPACER_2M) = NEGI;     
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_2M,0) = (-100);  
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_2M,1) = (-100);  
          GeneWise21_DC_SHADOW_MATRIX(mat,i,j,CENTRAL_0I) = NEGI;    
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_0I,0) = (-100); 
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_0I,1) = (-100); 
          GeneWise21_DC_SHADOW_MATRIX(mat,i,j,PY_0I) = NEGI;     
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_0I,0) = (-100);  
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_0I,1) = (-100);  
          GeneWise21_DC_SHADOW_MATRIX(mat,i,j,SPACER_0I) = NEGI;     
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_0I,0) = (-100);  
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_0I,1) = (-100);  
          GeneWise21_DC_SHADOW_MATRIX(mat,i,j,CENTRAL_1I) = NEGI;    
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_1I,0) = (-100); 
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_1I,1) = (-100); 
          GeneWise21_DC_SHADOW_MATRIX(mat,i,j,PY_1I) = NEGI;     
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_1I,0) = (-100);  
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_1I,1) = (-100);  
          GeneWise21_DC_SHADOW_MATRIX(mat,i,j,SPACER_1I) = NEGI;     
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_1I,0) = (-100);  
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_1I,1) = (-100);  
          GeneWise21_DC_SHADOW_MATRIX(mat,i,j,CENTRAL_2I) = NEGI;    
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_2I,0) = (-100); 
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_2I,1) = (-100); 
          GeneWise21_DC_SHADOW_MATRIX(mat,i,j,PY_2I) = NEGI;     
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_2I,0) = (-100);  
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_2I,1) = (-100);  
          GeneWise21_DC_SHADOW_MATRIX(mat,i,j,SPACER_2I) = NEGI;     
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_2I,0) = (-100);  
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_2I,1) = (-100);  
          continue;  
          } /* end of Is not in envelope */ 


        /* For state MATCH, pushing when j - offj <= mergej */ 
        score = GeneWise21_DC_SHADOW_MATRIX(mat,i-1,j-3,MATCH) + (mat->query->seg[i]->transition[GW_MATCH2MATCH]+mat->query->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]);  
        if( j - 3 <= mergej) {  
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,0) = i-1; 
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,1) = j-3; 
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,2) = MATCH;   
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,3) = i;   
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,4) = j;   
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,5) = MATCH;   
          }  
        else {  
          for(k=0;k<7;k++)   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 3,MATCH,k);   
          }  


        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-1,j-3,INSERT) + (mat->query->seg[i]->transition[GW_INSERT2MATCH]+mat->query->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]);     
        if( temp > score)    {  
          score = temp;  


          if( j - 3 <= mergej)   {  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,0) = i-1;   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,1) = j-3;   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,2) = INSERT;    
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,3) = i; 
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,4) = j; 
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,5) = MATCH; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 3,INSERT,k);    
            }  
          }  


        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-1,j-3,DELETE) + (mat->query->seg[i]->transition[GW_DELETE2MATCH]+mat->query->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]);     
        if( temp > score)    {  
          score = temp;  


          if( j - 3 <= mergej)   {  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,0) = i-1;   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,1) = j-3;   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,2) = DELETE;    
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,3) = i; 
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,4) = j; 
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,5) = MATCH; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 3,DELETE,k);    
            }  
          }  


        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-1,j-6,SPACER_0M) + (((mat->gp21->transition[GP21_SPACER2CDS]+mat->query->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)])+CSEQ_GENOMIC_3SS(mat->target,(j-3)))+mat->query->seg[i]->transition[GW_MATCH_BALANCE_3SS]);  
        if( temp > score)    {  
          score = temp;  


          if( j - 6 <= mergej)   {  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,0) = i-1;   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,1) = j-6;   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,2) = SPACER_0M; 
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,3) = i; 
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,4) = j; 
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,5) = MATCH; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 6,SPACER_0M,k); 
            }  
          }  


        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-1,j-6,PY_0M) + (((mat->gp21->transition[GP21_PY2CDS]+mat->query->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)])+CSEQ_GENOMIC_3SS(mat->target,(j-3)))+mat->query->seg[i]->transition[GW_MATCH_BALANCE_3SS]);  
        if( temp > score)    {  
          score = temp;  


          if( j - 6 <= mergej)   {  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,0) = i-1;   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,1) = j-6;   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,2) = PY_0M; 
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,3) = i; 
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,4) = j; 
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,5) = MATCH; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 6,PY_0M,k); 
            }  
          }  


        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-1,j-5,SPACER_1M) + (mat->gp21->transition[GP21_SPACER2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-2)));    
        if( temp > score)    {  
          score = temp;  


          if( j - 5 <= mergej)   {  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,0) = i-1;   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,1) = j-5;   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,2) = SPACER_1M; 
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,3) = i; 
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,4) = j; 
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,5) = MATCH; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 5,SPACER_1M,k); 
            }  
          }  


        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-1,j-5,PY_1M) + (mat->gp21->transition[GP21_PY2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-2)));    
        if( temp > score)    {  
          score = temp;  


          if( j - 5 <= mergej)   {  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,0) = i-1;   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,1) = j-5;   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,2) = PY_1M; 
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,3) = i; 
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,4) = j; 
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,5) = MATCH; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 5,PY_1M,k); 
            }  
          }  


        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-1,j-4,SPACER_2M) + (mat->gp21->transition[GP21_SPACER2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-1)));    
        if( temp > score)    {  
          score = temp;  


          if( j - 4 <= mergej)   {  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,0) = i-1;   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,1) = j-4;   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,2) = SPACER_2M; 
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,3) = i; 
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,4) = j; 
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,5) = MATCH; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 4,SPACER_2M,k); 
            }  
          }  


        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-1,j-4,PY_2M) + (mat->gp21->transition[GP21_PY2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-1)));    
        if( temp > score)    {  
          score = temp;  


          if( j - 4 <= mergej)   {  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,0) = i-1;   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,1) = j-4;   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,2) = PY_2M; 
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,3) = i; 
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,4) = j; 
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,5) = MATCH; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 4,PY_2M,k); 
            }  
          }  


        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-1,j-2,MATCH) + mat->gp21->transition[GP21_DELETE_1_BASE];   
        if( temp > score)    {  
          score = temp;  


          if( j - 2 <= mergej)   {  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,0) = i-1;   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,1) = j-2;   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,2) = MATCH; 
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,3) = i; 
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,4) = j; 
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,5) = MATCH; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 2,MATCH,k); 
            }  
          }  


        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-1,j-1,MATCH) + mat->gp21->transition[GP21_DELETE_2_BASE];   
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,0) = i-1;   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,1) = j-1;   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,2) = MATCH; 
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,3) = i; 
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,4) = j; 
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,5) = MATCH; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,MATCH,k); 
            }  
          }  


        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-1,j-4,MATCH) + mat->gp21->transition[GP21_INSERT_1_BASE];   
        if( temp > score)    {  
          score = temp;  


          if( j - 4 <= mergej)   {  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,0) = i-1;   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,1) = j-4;   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,2) = MATCH; 
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,3) = i; 
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,4) = j; 
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,5) = MATCH; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 4,MATCH,k); 
            }  
          }  


        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-1,j-5,MATCH) + mat->gp21->transition[GP21_INSERT_2_BASE];   
        if( temp > score)    {  
          score = temp;  


          if( j - 5 <= mergej)   {  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,0) = i-1;   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,1) = j-5;   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,2) = MATCH; 
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,3) = i; 
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,4) = j; 
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,5) = MATCH; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 5,MATCH,k); 
            }  
          }  
        /* Add any movement independant score */ 
        score += CSEQ_GENOMIC_CDSPOT(mat->target,j);     
        GeneWise21_DC_SHADOW_MATRIX(mat,i,j,MATCH) = score;  
        /* Finished with state MATCH */ 


        /* For state INSERT, pushing when j - offj <= mergej */ 
        score = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-3,MATCH) + (mat->query->seg[i]->transition[GW_MATCH2INSERT]+mat->query->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)]);    
        if( j - 3 <= mergej) {  
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,0) = i-0;    
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,1) = j-3;    
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,2) = MATCH;  
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,3) = i;  
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,4) = j;  
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,5) = INSERT; 
          }  
        else {  
          for(k=0;k<7;k++)   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,k) = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 3,MATCH,k);  
          }  


        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-3,INSERT) + (mat->query->seg[i]->transition[GW_INSERT2INSERT]+mat->query->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)]);   
        if( temp > score)    {  
          score = temp;  


          if( j - 3 <= mergej)   {  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,0) = i-0;  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,1) = j-3;  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,2) = INSERT;   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,3) = i;    
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,4) = j;    
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,5) = INSERT;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,k) = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 3,INSERT,k);   
            }  
          }  


        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-3,DELETE) + (mat->query->seg[i]->transition[GW_DELETE2INSERT]+mat->query->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)]);   
        if( temp > score)    {  
          score = temp;  


          if( j - 3 <= mergej)   {  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,0) = i-0;  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,1) = j-3;  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,2) = DELETE;   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,3) = i;    
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,4) = j;    
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,5) = INSERT;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,k) = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 3,DELETE,k);   
            }  
          }  


        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-1,j-6,SPACER_0I) + (((mat->gp21->transition[GP21_SPACER2CDS]+mat->query->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)])+CSEQ_GENOMIC_3SS(mat->target,(j-3)))+mat->query->seg[i]->transition[GW_MATCH_BALANCE_3SS]);     
        if( temp > score)    {  
          score = temp;  


          if( j - 6 <= mergej)   {  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,0) = i-1;  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,1) = j-6;  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,2) = SPACER_0I;    
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,3) = i;    
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,4) = j;    
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,5) = INSERT;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,k) = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 6,SPACER_0I,k);    
            }  
          }  


        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-1,j-5,SPACER_1I) + (mat->gp21->transition[GP21_SPACER2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-2)));    
        if( temp > score)    {  
          score = temp;  


          if( j - 5 <= mergej)   {  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,0) = i-1;  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,1) = j-5;  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,2) = SPACER_1I;    
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,3) = i;    
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,4) = j;    
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,5) = INSERT;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,k) = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 5,SPACER_1I,k);    
            }  
          }  


        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-1,j-4,SPACER_2I) + (mat->gp21->transition[GP21_SPACER2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-1)));    
        if( temp > score)    {  
          score = temp;  


          if( j - 4 <= mergej)   {  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,0) = i-1;  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,1) = j-4;  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,2) = SPACER_2I;    
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,3) = i;    
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,4) = j;    
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,5) = INSERT;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,k) = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 4,SPACER_2I,k);    
            }  
          }  


        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-1,j-6,PY_0I) + (((mat->gp21->transition[GP21_PY2CDS]+mat->query->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)])+CSEQ_GENOMIC_3SS(mat->target,(j-3)))+mat->query->seg[i]->transition[GW_MATCH_BALANCE_3SS]);     
        if( temp > score)    {  
          score = temp;  


          if( j - 6 <= mergej)   {  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,0) = i-1;  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,1) = j-6;  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,2) = PY_0I;    
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,3) = i;    
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,4) = j;    
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,5) = INSERT;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,k) = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 6,PY_0I,k);    
            }  
          }  


        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-1,j-5,PY_1I) + (mat->gp21->transition[GP21_PY2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-2)));    
        if( temp > score)    {  
          score = temp;  


          if( j - 5 <= mergej)   {  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,0) = i-1;  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,1) = j-5;  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,2) = PY_1I;    
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,3) = i;    
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,4) = j;    
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,5) = INSERT;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,k) = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 5,PY_1I,k);    
            }  
          }  


        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-1,j-4,PY_2I) + (mat->gp21->transition[GP21_PY2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-1)));    
        if( temp > score)    {  
          score = temp;  


          if( j - 4 <= mergej)   {  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,0) = i-1;  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,1) = j-4;  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,2) = PY_2I;    
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,3) = i;    
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,4) = j;    
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,5) = INSERT;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,k) = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 4,PY_2I,k);    
            }  
          }  


        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-1,j-2,INSERT) + mat->gp21->transition[GP21_DELETE_1_BASE];  
        if( temp > score)    {  
          score = temp;  


          if( j - 2 <= mergej)   {  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,0) = i-1;  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,1) = j-2;  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,2) = INSERT;   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,3) = i;    
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,4) = j;    
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,5) = INSERT;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,k) = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 2,INSERT,k);   
            }  
          }  


        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-1,j-1,INSERT) + mat->gp21->transition[GP21_DELETE_2_BASE];  
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,0) = i-1;  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,1) = j-1;  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,2) = INSERT;   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,3) = i;    
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,4) = j;    
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,5) = INSERT;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,k) = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,INSERT,k);   
            }  
          }  


        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-1,j-4,INSERT) + mat->gp21->transition[GP21_INSERT_1_BASE];  
        if( temp > score)    {  
          score = temp;  


          if( j - 4 <= mergej)   {  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,0) = i-1;  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,1) = j-4;  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,2) = INSERT;   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,3) = i;    
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,4) = j;    
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,5) = INSERT;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,k) = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 4,INSERT,k);   
            }  
          }  


        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-1,j-5,INSERT) + mat->gp21->transition[GP21_INSERT_2_BASE];  
        if( temp > score)    {  
          score = temp;  


          if( j - 5 <= mergej)   {  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,0) = i-1;  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,1) = j-5;  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,2) = INSERT;   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,3) = i;    
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,4) = j;    
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,5) = INSERT;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,k) = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 5,INSERT,k);   
            }  
          }  
        /* Add any movement independant score */ 
        score += CSEQ_GENOMIC_CDSPOT(mat->target,j);     
        GeneWise21_DC_SHADOW_MATRIX(mat,i,j,INSERT) = score;     
        /* Finished with state INSERT */ 


        /* For state DELETE, pushing when j - offj <= mergej */ 
        score = GeneWise21_DC_SHADOW_MATRIX(mat,i-1,j-0,MATCH) + mat->query->seg[i]->transition[GW_MATCH2DELETE];    
        if( j - 0 <= mergej) {  
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,0) = i-1;    
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,1) = j-0;    
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,2) = MATCH;  
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,3) = i;  
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,4) = j;  
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,5) = DELETE; 
          }  
        else {  
          for(k=0;k<7;k++)   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,k) = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 0,MATCH,k);  
          }  


        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-1,j-0,INSERT) + mat->query->seg[i]->transition[GW_INSERT2DELETE];   
        if( temp > score)    {  
          score = temp;  


          if( j - 0 <= mergej)   {  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,0) = i-1;  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,1) = j-0;  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,2) = INSERT;   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,3) = i;    
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,4) = j;    
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,5) = DELETE;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,k) = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 0,INSERT,k);   
            }  
          }  


        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-1,j-0,DELETE) + mat->query->seg[i]->transition[GW_DELETE2DELETE];   
        if( temp > score)    {  
          score = temp;  


          if( j - 0 <= mergej)   {  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,0) = i-1;  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,1) = j-0;  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,2) = DELETE;   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,3) = i;    
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,4) = j;    
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,5) = DELETE;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,k) = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 0,DELETE,k);   
            }  
          }  
        /* Add any movement independant score */ 
        GeneWise21_DC_SHADOW_MATRIX(mat,i,j,DELETE) = score;     
        /* Finished with state DELETE */ 


        /* For state CENTRAL_0M, pushing when j - offj <= mergej */ 
        score = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-8,MATCH) + (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_MATCH2MATCH])+mat->query->seg[i]->transition[GW_MATCH_BALANCE_5SS]);     
        if( j - 8 <= mergej) {  
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_0M,0) = i-0;    
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_0M,1) = j-8;    
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_0M,2) = MATCH;  
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_0M,3) = i;  
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_0M,4) = j;  
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_0M,5) = CENTRAL_0M; 
          }  
        else {  
          for(k=0;k<7;k++)   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_0M,k) = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 8,MATCH,k);  
          }  


        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-8,INSERT) + (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_INSERT2MATCH])+mat->query->seg[i]->transition[GW_INSERT_BALANCE_5SS]);   
        if( temp > score)    {  
          score = temp;  


          if( j - 8 <= mergej)   {  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_0M,0) = i-0;  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_0M,1) = j-8;  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_0M,2) = INSERT;   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_0M,3) = i;    
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_0M,4) = j;    
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_0M,5) = CENTRAL_0M;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_0M,k) = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 8,INSERT,k);   
            }  
          }  


        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-8,DELETE) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_DELETE2MATCH]);   
        if( temp > score)    {  
          score = temp;  


          if( j - 8 <= mergej)   {  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_0M,0) = i-0;  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_0M,1) = j-8;  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_0M,2) = DELETE;   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_0M,3) = i;    
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_0M,4) = j;    
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_0M,5) = CENTRAL_0M;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_0M,k) = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 8,DELETE,k);   
            }  
          }  


        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,CENTRAL_0M) + (mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CENTRAL2CENTRAL]);     
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_0M,0) = i-0;  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_0M,1) = j-1;  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_0M,2) = CENTRAL_0M;   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_0M,3) = i;    
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_0M,4) = j;    
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_0M,5) = CENTRAL_0M;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_0M,k) = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,CENTRAL_0M,k);   
            }  
          }  
        /* Add any movement independant score */ 
        GeneWise21_DC_SHADOW_MATRIX(mat,i,j,CENTRAL_0M) = score;     
        /* Finished with state CENTRAL_0M */ 


        /* For state PY_0M, pushing when j - offj <= mergej */ 
        score = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,CENTRAL_0M) + mat->gp21->transition[GP21_CENTRAL2PY];    
        if( j - 1 <= mergej) {  
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_0M,0) = i-0; 
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_0M,1) = j-1; 
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_0M,2) = CENTRAL_0M;  
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_0M,3) = i;   
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_0M,4) = j;   
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_0M,5) = PY_0M;   
          }  
        else {  
          for(k=0;k<7;k++)   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_0M,k) = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,CENTRAL_0M,k);  
          }  


        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,PY_0M) + mat->gp21->transition[GP21_PY2PY];   
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_0M,0) = i-0;   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_0M,1) = j-1;   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_0M,2) = PY_0M; 
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_0M,3) = i; 
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_0M,4) = j; 
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_0M,5) = PY_0M; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_0M,k) = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,PY_0M,k); 
            }  
          }  
        /* Add any movement independant score */ 
        score += mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)];    
        GeneWise21_DC_SHADOW_MATRIX(mat,i,j,PY_0M) = score;  
        /* Finished with state PY_0M */ 


        /* For state SPACER_0M, pushing when j - offj <= mergej */ 
        score = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,PY_0M) + mat->gp21->transition[GP21_PY2SPACER];  
        if( j - 1 <= mergej) {  
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_0M,0) = i-0; 
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_0M,1) = j-1; 
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_0M,2) = PY_0M;   
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_0M,3) = i;   
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_0M,4) = j;   
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_0M,5) = SPACER_0M;   
          }  
        else {  
          for(k=0;k<7;k++)   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_0M,k) = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,PY_0M,k);   
          }  


        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,SPACER_0M) + mat->gp21->transition[GP21_SPACER2SPACER];   
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_0M,0) = i-0;   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_0M,1) = j-1;   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_0M,2) = SPACER_0M; 
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_0M,3) = i; 
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_0M,4) = j; 
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_0M,5) = SPACER_0M; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_0M,k) = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,SPACER_0M,k); 
            }  
          }  
        /* Add any movement independant score */ 
        score += mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)];    
        GeneWise21_DC_SHADOW_MATRIX(mat,i,j,SPACER_0M) = score;  
        /* Finished with state SPACER_0M */ 


        /* For state CENTRAL_1M, pushing when j - offj <= mergej */ 
        score = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-9,MATCH) + (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_MATCH2MATCH])+mat->query->seg[i]->transition[GW_MATCH_BALANCE_5SS]);     
        if( j - 9 <= mergej) {  
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_1M,0) = i-0;    
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_1M,1) = j-9;    
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_1M,2) = MATCH;  
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_1M,3) = i;  
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_1M,4) = j;  
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_1M,5) = CENTRAL_1M; 
          }  
        else {  
          for(k=0;k<7;k++)   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_1M,k) = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 9,MATCH,k);  
          }  


        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-9,INSERT) + (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_INSERT2MATCH])+mat->query->seg[i]->transition[GW_MATCH_BALANCE_5SS]);    
        if( temp > score)    {  
          score = temp;  


          if( j - 9 <= mergej)   {  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_1M,0) = i-0;  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_1M,1) = j-9;  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_1M,2) = INSERT;   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_1M,3) = i;    
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_1M,4) = j;    
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_1M,5) = CENTRAL_1M;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_1M,k) = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 9,INSERT,k);   
            }  
          }  


        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-9,DELETE) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_DELETE2MATCH]);   
        if( temp > score)    {  
          score = temp;  


          if( j - 9 <= mergej)   {  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_1M,0) = i-0;  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_1M,1) = j-9;  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_1M,2) = DELETE;   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_1M,3) = i;    
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_1M,4) = j;    
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_1M,5) = CENTRAL_1M;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_1M,k) = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 9,DELETE,k);   
            }  
          }  


        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,CENTRAL_1M) + (mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CENTRAL2CENTRAL]);     
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_1M,0) = i-0;  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_1M,1) = j-1;  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_1M,2) = CENTRAL_1M;   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_1M,3) = i;    
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_1M,4) = j;    
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_1M,5) = CENTRAL_1M;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_1M,k) = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,CENTRAL_1M,k);   
            }  
          }  
        /* Add any movement independant score */ 
        GeneWise21_DC_SHADOW_MATRIX(mat,i,j,CENTRAL_1M) = score;     
        /* Finished with state CENTRAL_1M */ 


        /* For state PY_1M, pushing when j - offj <= mergej */ 
        score = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,CENTRAL_1M) + mat->gp21->transition[GP21_CENTRAL2PY];    
        if( j - 1 <= mergej) {  
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_1M,0) = i-0; 
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_1M,1) = j-1; 
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_1M,2) = CENTRAL_1M;  
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_1M,3) = i;   
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_1M,4) = j;   
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_1M,5) = PY_1M;   
          }  
        else {  
          for(k=0;k<7;k++)   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_1M,k) = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,CENTRAL_1M,k);  
          }  


        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,PY_1M) + mat->gp21->transition[GP21_PY2PY];   
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_1M,0) = i-0;   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_1M,1) = j-1;   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_1M,2) = PY_1M; 
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_1M,3) = i; 
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_1M,4) = j; 
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_1M,5) = PY_1M; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_1M,k) = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,PY_1M,k); 
            }  
          }  
        /* Add any movement independant score */ 
        score += mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)];    
        GeneWise21_DC_SHADOW_MATRIX(mat,i,j,PY_1M) = score;  
        /* Finished with state PY_1M */ 


        /* For state SPACER_1M, pushing when j - offj <= mergej */ 
        score = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,PY_1M) + mat->gp21->transition[GP21_PY2SPACER];  
        if( j - 1 <= mergej) {  
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_1M,0) = i-0; 
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_1M,1) = j-1; 
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_1M,2) = PY_1M;   
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_1M,3) = i;   
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_1M,4) = j;   
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_1M,5) = SPACER_1M;   
          }  
        else {  
          for(k=0;k<7;k++)   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_1M,k) = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,PY_1M,k);   
          }  


        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,SPACER_1M) + mat->gp21->transition[GP21_SPACER2SPACER];   
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_1M,0) = i-0;   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_1M,1) = j-1;   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_1M,2) = SPACER_1M; 
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_1M,3) = i; 
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_1M,4) = j; 
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_1M,5) = SPACER_1M; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_1M,k) = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,SPACER_1M,k); 
            }  
          }  
        /* Add any movement independant score */ 
        score += mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)];    
        GeneWise21_DC_SHADOW_MATRIX(mat,i,j,SPACER_1M) = score;  
        /* Finished with state SPACER_1M */ 


        /* For state CENTRAL_2M, pushing when j - offj <= mergej */ 
        score = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-10,MATCH) + (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_MATCH2MATCH])+mat->query->seg[i]->transition[GW_MATCH_BALANCE_5SS]);    
        if( j - 10 <= mergej)    {  
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_2M,0) = i-0;    
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_2M,1) = j-10;   
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_2M,2) = MATCH;  
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_2M,3) = i;  
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_2M,4) = j;  
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_2M,5) = CENTRAL_2M; 
          }  
        else {  
          for(k=0;k<7;k++)   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_2M,k) = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 10,MATCH,k); 
          }  


        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-10,INSERT) + (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_INSERT2MATCH])+mat->query->seg[i]->transition[GW_MATCH_BALANCE_5SS]);   
        if( temp > score)    {  
          score = temp;  


          if( j - 10 <= mergej)  {  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_2M,0) = i-0;  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_2M,1) = j-10; 
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_2M,2) = INSERT;   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_2M,3) = i;    
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_2M,4) = j;    
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_2M,5) = CENTRAL_2M;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_2M,k) = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 10,INSERT,k);  
            }  
          }  


        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-10,DELETE) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_DELETE2MATCH]);  
        if( temp > score)    {  
          score = temp;  


          if( j - 10 <= mergej)  {  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_2M,0) = i-0;  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_2M,1) = j-10; 
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_2M,2) = DELETE;   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_2M,3) = i;    
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_2M,4) = j;    
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_2M,5) = CENTRAL_2M;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_2M,k) = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 10,DELETE,k);  
            }  
          }  


        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,CENTRAL_2M) + (mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CENTRAL2CENTRAL]);     
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_2M,0) = i-0;  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_2M,1) = j-1;  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_2M,2) = CENTRAL_2M;   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_2M,3) = i;    
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_2M,4) = j;    
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_2M,5) = CENTRAL_2M;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_2M,k) = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,CENTRAL_2M,k);   
            }  
          }  
        /* Add any movement independant score */ 
        GeneWise21_DC_SHADOW_MATRIX(mat,i,j,CENTRAL_2M) = score;     
        /* Finished with state CENTRAL_2M */ 


        /* For state PY_2M, pushing when j - offj <= mergej */ 
        score = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,CENTRAL_2M) + mat->gp21->transition[GP21_CENTRAL2PY];    
        if( j - 1 <= mergej) {  
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_2M,0) = i-0; 
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_2M,1) = j-1; 
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_2M,2) = CENTRAL_2M;  
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_2M,3) = i;   
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_2M,4) = j;   
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_2M,5) = PY_2M;   
          }  
        else {  
          for(k=0;k<7;k++)   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_2M,k) = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,CENTRAL_2M,k);  
          }  


        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,PY_2M) + mat->gp21->transition[GP21_PY2PY];   
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_2M,0) = i-0;   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_2M,1) = j-1;   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_2M,2) = PY_2M; 
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_2M,3) = i; 
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_2M,4) = j; 
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_2M,5) = PY_2M; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_2M,k) = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,PY_2M,k); 
            }  
          }  
        /* Add any movement independant score */ 
        score += mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)];    
        GeneWise21_DC_SHADOW_MATRIX(mat,i,j,PY_2M) = score;  
        /* Finished with state PY_2M */ 


        /* For state SPACER_2M, pushing when j - offj <= mergej */ 
        score = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,PY_2M) + mat->gp21->transition[GP21_PY2SPACER];  
        if( j - 1 <= mergej) {  
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_2M,0) = i-0; 
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_2M,1) = j-1; 
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_2M,2) = PY_2M;   
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_2M,3) = i;   
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_2M,4) = j;   
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_2M,5) = SPACER_2M;   
          }  
        else {  
          for(k=0;k<7;k++)   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_2M,k) = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,PY_2M,k);   
          }  


        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,SPACER_2M) + mat->gp21->transition[GP21_SPACER2SPACER];   
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_2M,0) = i-0;   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_2M,1) = j-1;   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_2M,2) = SPACER_2M; 
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_2M,3) = i; 
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_2M,4) = j; 
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_2M,5) = SPACER_2M; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_2M,k) = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,SPACER_2M,k); 
            }  
          }  
        /* Add any movement independant score */ 
        score += mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)];    
        GeneWise21_DC_SHADOW_MATRIX(mat,i,j,SPACER_2M) = score;  
        /* Finished with state SPACER_2M */ 


        /* For state CENTRAL_0I, pushing when j - offj <= mergej */ 
        score = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-8,MATCH) + (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_MATCH2INSERT])+mat->query->seg[i]->transition[GW_MATCH_BALANCE_5SS]);    
        if( j - 8 <= mergej) {  
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_0I,0) = i-0;    
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_0I,1) = j-8;    
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_0I,2) = MATCH;  
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_0I,3) = i;  
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_0I,4) = j;  
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_0I,5) = CENTRAL_0I; 
          }  
        else {  
          for(k=0;k<7;k++)   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_0I,k) = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 8,MATCH,k);  
          }  


        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-8,INSERT) + (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_INSERT2INSERT])+mat->query->seg[i]->transition[GW_INSERT_BALANCE_5SS]);  
        if( temp > score)    {  
          score = temp;  


          if( j - 8 <= mergej)   {  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_0I,0) = i-0;  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_0I,1) = j-8;  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_0I,2) = INSERT;   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_0I,3) = i;    
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_0I,4) = j;    
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_0I,5) = CENTRAL_0I;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_0I,k) = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 8,INSERT,k);   
            }  
          }  


        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-8,DELETE) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_DELETE2INSERT]);  
        if( temp > score)    {  
          score = temp;  


          if( j - 8 <= mergej)   {  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_0I,0) = i-0;  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_0I,1) = j-8;  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_0I,2) = DELETE;   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_0I,3) = i;    
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_0I,4) = j;    
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_0I,5) = CENTRAL_0I;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_0I,k) = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 8,DELETE,k);   
            }  
          }  


        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,CENTRAL_0I) + (mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CENTRAL2CENTRAL]);     
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_0I,0) = i-0;  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_0I,1) = j-1;  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_0I,2) = CENTRAL_0I;   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_0I,3) = i;    
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_0I,4) = j;    
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_0I,5) = CENTRAL_0I;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_0I,k) = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,CENTRAL_0I,k);   
            }  
          }  
        /* Add any movement independant score */ 
        GeneWise21_DC_SHADOW_MATRIX(mat,i,j,CENTRAL_0I) = score;     
        /* Finished with state CENTRAL_0I */ 


        /* For state PY_0I, pushing when j - offj <= mergej */ 
        score = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,CENTRAL_0I) + mat->gp21->transition[GP21_CENTRAL2PY];    
        if( j - 1 <= mergej) {  
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_0I,0) = i-0; 
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_0I,1) = j-1; 
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_0I,2) = CENTRAL_0I;  
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_0I,3) = i;   
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_0I,4) = j;   
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_0I,5) = PY_0I;   
          }  
        else {  
          for(k=0;k<7;k++)   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_0I,k) = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,CENTRAL_0I,k);  
          }  


        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,PY_0I) + mat->gp21->transition[GP21_PY2PY];   
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_0I,0) = i-0;   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_0I,1) = j-1;   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_0I,2) = PY_0I; 
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_0I,3) = i; 
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_0I,4) = j; 
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_0I,5) = PY_0I; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_0I,k) = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,PY_0I,k); 
            }  
          }  
        /* Add any movement independant score */ 
        score += mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)];    
        GeneWise21_DC_SHADOW_MATRIX(mat,i,j,PY_0I) = score;  
        /* Finished with state PY_0I */ 


        /* For state SPACER_0I, pushing when j - offj <= mergej */ 
        score = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,PY_0I) + mat->gp21->transition[GP21_PY2SPACER];  
        if( j - 1 <= mergej) {  
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_0I,0) = i-0; 
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_0I,1) = j-1; 
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_0I,2) = PY_0I;   
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_0I,3) = i;   
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_0I,4) = j;   
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_0I,5) = SPACER_0I;   
          }  
        else {  
          for(k=0;k<7;k++)   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_0I,k) = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,PY_0I,k);   
          }  


        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,SPACER_0I) + mat->gp21->transition[GP21_SPACER2SPACER];   
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_0I,0) = i-0;   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_0I,1) = j-1;   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_0I,2) = SPACER_0I; 
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_0I,3) = i; 
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_0I,4) = j; 
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_0I,5) = SPACER_0I; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_0I,k) = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,SPACER_0I,k); 
            }  
          }  
        /* Add any movement independant score */ 
        score += mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)];    
        GeneWise21_DC_SHADOW_MATRIX(mat,i,j,SPACER_0I) = score;  
        /* Finished with state SPACER_0I */ 


        /* For state CENTRAL_1I, pushing when j - offj <= mergej */ 
        score = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-9,MATCH) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_MATCH2INSERT]);   
        if( j - 9 <= mergej) {  
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_1I,0) = i-0;    
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_1I,1) = j-9;    
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_1I,2) = MATCH;  
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_1I,3) = i;  
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_1I,4) = j;  
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_1I,5) = CENTRAL_1I; 
          }  
        else {  
          for(k=0;k<7;k++)   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_1I,k) = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 9,MATCH,k);  
          }  


        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-9,INSERT) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_INSERT2INSERT]);  
        if( temp > score)    {  
          score = temp;  


          if( j - 9 <= mergej)   {  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_1I,0) = i-0;  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_1I,1) = j-9;  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_1I,2) = INSERT;   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_1I,3) = i;    
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_1I,4) = j;    
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_1I,5) = CENTRAL_1I;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_1I,k) = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 9,INSERT,k);   
            }  
          }  


        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-9,DELETE) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_DELETE2INSERT]);  
        if( temp > score)    {  
          score = temp;  


          if( j - 9 <= mergej)   {  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_1I,0) = i-0;  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_1I,1) = j-9;  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_1I,2) = DELETE;   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_1I,3) = i;    
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_1I,4) = j;    
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_1I,5) = CENTRAL_1I;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_1I,k) = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 9,DELETE,k);   
            }  
          }  


        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,CENTRAL_1I) + (mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CENTRAL2CENTRAL]);     
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_1I,0) = i-0;  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_1I,1) = j-1;  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_1I,2) = CENTRAL_1I;   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_1I,3) = i;    
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_1I,4) = j;    
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_1I,5) = CENTRAL_1I;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_1I,k) = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,CENTRAL_1I,k);   
            }  
          }  
        /* Add any movement independant score */ 
        GeneWise21_DC_SHADOW_MATRIX(mat,i,j,CENTRAL_1I) = score;     
        /* Finished with state CENTRAL_1I */ 


        /* For state PY_1I, pushing when j - offj <= mergej */ 
        score = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,CENTRAL_1I) + mat->gp21->transition[GP21_CENTRAL2PY];    
        if( j - 1 <= mergej) {  
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_1I,0) = i-0; 
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_1I,1) = j-1; 
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_1I,2) = CENTRAL_1I;  
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_1I,3) = i;   
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_1I,4) = j;   
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_1I,5) = PY_1I;   
          }  
        else {  
          for(k=0;k<7;k++)   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_1I,k) = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,CENTRAL_1I,k);  
          }  


        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,PY_1I) + mat->gp21->transition[GP21_PY2PY];   
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_1I,0) = i-0;   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_1I,1) = j-1;   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_1I,2) = PY_1I; 
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_1I,3) = i; 
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_1I,4) = j; 
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_1I,5) = PY_1I; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_1I,k) = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,PY_1I,k); 
            }  
          }  
        /* Add any movement independant score */ 
        score += mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)];    
        GeneWise21_DC_SHADOW_MATRIX(mat,i,j,PY_1I) = score;  
        /* Finished with state PY_1I */ 


        /* For state SPACER_1I, pushing when j - offj <= mergej */ 
        score = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,PY_1I) + mat->gp21->transition[GP21_PY2SPACER];  
        if( j - 1 <= mergej) {  
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_1I,0) = i-0; 
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_1I,1) = j-1; 
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_1I,2) = PY_1I;   
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_1I,3) = i;   
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_1I,4) = j;   
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_1I,5) = SPACER_1I;   
          }  
        else {  
          for(k=0;k<7;k++)   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_1I,k) = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,PY_1I,k);   
          }  


        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,SPACER_1I) + mat->gp21->transition[GP21_SPACER2SPACER];   
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_1I,0) = i-0;   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_1I,1) = j-1;   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_1I,2) = SPACER_1I; 
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_1I,3) = i; 
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_1I,4) = j; 
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_1I,5) = SPACER_1I; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_1I,k) = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,SPACER_1I,k); 
            }  
          }  
        /* Add any movement independant score */ 
        score += mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)];    
        GeneWise21_DC_SHADOW_MATRIX(mat,i,j,SPACER_1I) = score;  
        /* Finished with state SPACER_1I */ 


        /* For state CENTRAL_2I, pushing when j - offj <= mergej */ 
        score = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-10,MATCH) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_MATCH2INSERT]);  
        if( j - 10 <= mergej)    {  
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_2I,0) = i-0;    
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_2I,1) = j-10;   
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_2I,2) = MATCH;  
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_2I,3) = i;  
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_2I,4) = j;  
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_2I,5) = CENTRAL_2I; 
          }  
        else {  
          for(k=0;k<7;k++)   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_2I,k) = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 10,MATCH,k); 
          }  


        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-10,INSERT) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_INSERT2INSERT]);     
        if( temp > score)    {  
          score = temp;  


          if( j - 10 <= mergej)  {  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_2I,0) = i-0;  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_2I,1) = j-10; 
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_2I,2) = INSERT;   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_2I,3) = i;    
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_2I,4) = j;    
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_2I,5) = CENTRAL_2I;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_2I,k) = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 10,INSERT,k);  
            }  
          }  


        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-10,INSERT) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_DELETE2INSERT]);     
        if( temp > score)    {  
          score = temp;  


          if( j - 10 <= mergej)  {  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_2I,0) = i-0;  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_2I,1) = j-10; 
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_2I,2) = INSERT;   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_2I,3) = i;    
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_2I,4) = j;    
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_2I,5) = CENTRAL_2I;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_2I,k) = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 10,INSERT,k);  
            }  
          }  


        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,CENTRAL_2I) + (mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CENTRAL2CENTRAL]);     
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_2I,0) = i-0;  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_2I,1) = j-1;  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_2I,2) = CENTRAL_2I;   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_2I,3) = i;    
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_2I,4) = j;    
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_2I,5) = CENTRAL_2I;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_2I,k) = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,CENTRAL_2I,k);   
            }  
          }  
        /* Add any movement independant score */ 
        GeneWise21_DC_SHADOW_MATRIX(mat,i,j,CENTRAL_2I) = score;     
        /* Finished with state CENTRAL_2I */ 


        /* For state PY_2I, pushing when j - offj <= mergej */ 
        score = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,CENTRAL_2I) + mat->gp21->transition[GP21_CENTRAL2PY];    
        if( j - 1 <= mergej) {  
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_2I,0) = i-0; 
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_2I,1) = j-1; 
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_2I,2) = CENTRAL_2I;  
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_2I,3) = i;   
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_2I,4) = j;   
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_2I,5) = PY_2I;   
          }  
        else {  
          for(k=0;k<7;k++)   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_2I,k) = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,CENTRAL_2I,k);  
          }  


        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,PY_2I) + mat->gp21->transition[GP21_PY2PY];   
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_2I,0) = i-0;   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_2I,1) = j-1;   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_2I,2) = PY_2I; 
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_2I,3) = i; 
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_2I,4) = j; 
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_2I,5) = PY_2I; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_2I,k) = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,PY_2I,k); 
            }  
          }  
        /* Add any movement independant score */ 
        score += mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)];    
        GeneWise21_DC_SHADOW_MATRIX(mat,i,j,PY_2I) = score;  
        /* Finished with state PY_2I */ 


        /* For state SPACER_2I, pushing when j - offj <= mergej */ 
        score = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,PY_2I) + mat->gp21->transition[GP21_PY2SPACER];  
        if( j - 1 <= mergej) {  
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_2I,0) = i-0; 
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_2I,1) = j-1; 
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_2I,2) = PY_2I;   
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_2I,3) = i;   
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_2I,4) = j;   
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_2I,5) = SPACER_2I;   
          }  
        else {  
          for(k=0;k<7;k++)   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_2I,k) = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,PY_2I,k);   
          }  


        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,SPACER_2I) + mat->gp21->transition[GP21_SPACER2SPACER];   
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_2I,0) = i-0;   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_2I,1) = j-1;   
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_2I,2) = SPACER_2I; 
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_2I,3) = i; 
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_2I,4) = j; 
            GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_2I,5) = SPACER_2I; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_2I,k) = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,SPACER_2I,k); 
            }  
          }  
        /* Add any movement independant score */ 
        score += mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)];    
        GeneWise21_DC_SHADOW_MATRIX(mat,i,j,SPACER_2I) = score;  
        /* Finished with state SPACER_2I */ 
        }  
      }  
    /* Put back j into * stop j so that calling function gets it correct */ 
    if( stopj == NULL)   
      warn("Bad news... NULL stopj pointer in push dc function. This means that calling function does not know how many cells I have done!");    
    else 
      *stopj = j;    


    return;  
}    


/* Function:  follow_on_dc_GeneWise21(mat,starti,stopi,startj,stopj,dpenv,perc_done)
 *
 * Descrip: No Description
 *
 * Arg:              mat [UNKN ] Undocumented argument [GeneWise21 *]
 * Arg:           starti [UNKN ] Undocumented argument [int]
 * Arg:            stopi [UNKN ] Undocumented argument [int]
 * Arg:           startj [UNKN ] Undocumented argument [int]
 * Arg:            stopj [UNKN ] Undocumented argument [int]
 * Arg:            dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:        perc_done [UNKN ] Undocumented argument [int]
 *
 */
void follow_on_dc_GeneWise21(GeneWise21 * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done) 
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
          GeneWise21_DC_SHADOW_MATRIX(mat,i,j,MATCH) = NEGI;     
          GeneWise21_DC_SHADOW_MATRIX(mat,i,j,INSERT) = NEGI;    
          GeneWise21_DC_SHADOW_MATRIX(mat,i,j,DELETE) = NEGI;    
          GeneWise21_DC_SHADOW_MATRIX(mat,i,j,CENTRAL_0M) = NEGI;    
          GeneWise21_DC_SHADOW_MATRIX(mat,i,j,PY_0M) = NEGI;     
          GeneWise21_DC_SHADOW_MATRIX(mat,i,j,SPACER_0M) = NEGI;     
          GeneWise21_DC_SHADOW_MATRIX(mat,i,j,CENTRAL_1M) = NEGI;    
          GeneWise21_DC_SHADOW_MATRIX(mat,i,j,PY_1M) = NEGI;     
          GeneWise21_DC_SHADOW_MATRIX(mat,i,j,SPACER_1M) = NEGI;     
          GeneWise21_DC_SHADOW_MATRIX(mat,i,j,CENTRAL_2M) = NEGI;    
          GeneWise21_DC_SHADOW_MATRIX(mat,i,j,PY_2M) = NEGI;     
          GeneWise21_DC_SHADOW_MATRIX(mat,i,j,SPACER_2M) = NEGI;     
          GeneWise21_DC_SHADOW_MATRIX(mat,i,j,CENTRAL_0I) = NEGI;    
          GeneWise21_DC_SHADOW_MATRIX(mat,i,j,PY_0I) = NEGI;     
          GeneWise21_DC_SHADOW_MATRIX(mat,i,j,SPACER_0I) = NEGI;     
          GeneWise21_DC_SHADOW_MATRIX(mat,i,j,CENTRAL_1I) = NEGI;    
          GeneWise21_DC_SHADOW_MATRIX(mat,i,j,PY_1I) = NEGI;     
          GeneWise21_DC_SHADOW_MATRIX(mat,i,j,SPACER_1I) = NEGI;     
          GeneWise21_DC_SHADOW_MATRIX(mat,i,j,CENTRAL_2I) = NEGI;    
          GeneWise21_DC_SHADOW_MATRIX(mat,i,j,PY_2I) = NEGI;     
          GeneWise21_DC_SHADOW_MATRIX(mat,i,j,SPACER_2I) = NEGI;     
          continue;  
          } /* end of Is not in envelope */ 
        if( num % 1000 == 0 )    
          log_full_error(REPORT,0,"[%d%%%% done]After  mid-j %5d Cells done %d%%%%",perc_done,startj,(num*100)/total);   


        /* For state MATCH */ 
        /* setting first movement to score */ 
        score = GeneWise21_DC_SHADOW_MATRIX(mat,i-1,j-3,MATCH) + (mat->query->seg[i]->transition[GW_MATCH2MATCH]+mat->query->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]);  
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 3,MATCH,k);  
        /* From state INSERT to state MATCH */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-1,j-3,INSERT) + (mat->query->seg[i]->transition[GW_INSERT2MATCH]+mat->query->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]);     
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 3,INSERT,k);   
          }  
        /* From state DELETE to state MATCH */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-1,j-3,DELETE) + (mat->query->seg[i]->transition[GW_DELETE2MATCH]+mat->query->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]);     
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 3,DELETE,k);   
          }  
        /* From state SPACER_0M to state MATCH */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-1,j-6,SPACER_0M) + (((mat->gp21->transition[GP21_SPACER2CDS]+mat->query->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)])+CSEQ_GENOMIC_3SS(mat->target,(j-3)))+mat->query->seg[i]->transition[GW_MATCH_BALANCE_3SS]);  
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 6,SPACER_0M,k);    
          }  
        /* From state PY_0M to state MATCH */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-1,j-6,PY_0M) + (((mat->gp21->transition[GP21_PY2CDS]+mat->query->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)])+CSEQ_GENOMIC_3SS(mat->target,(j-3)))+mat->query->seg[i]->transition[GW_MATCH_BALANCE_3SS]);  
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 6,PY_0M,k);    
          }  
        /* From state SPACER_1M to state MATCH */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-1,j-5,SPACER_1M) + (mat->gp21->transition[GP21_SPACER2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-2)));    
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 5,SPACER_1M,k);    
          }  
        /* From state PY_1M to state MATCH */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-1,j-5,PY_1M) + (mat->gp21->transition[GP21_PY2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-2)));    
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 5,PY_1M,k);    
          }  
        /* From state SPACER_2M to state MATCH */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-1,j-4,SPACER_2M) + (mat->gp21->transition[GP21_SPACER2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-1)));    
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 4,SPACER_2M,k);    
          }  
        /* From state PY_2M to state MATCH */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-1,j-4,PY_2M) + (mat->gp21->transition[GP21_PY2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-1)));    
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 4,PY_2M,k);    
          }  
        /* From state MATCH to state MATCH */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-1,j-2,MATCH) + mat->gp21->transition[GP21_DELETE_1_BASE];   
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 2,MATCH,k);    
          }  
        /* From state MATCH to state MATCH */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-1,j-1,MATCH) + mat->gp21->transition[GP21_DELETE_2_BASE];   
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,MATCH,k);    
          }  
        /* From state MATCH to state MATCH */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-1,j-4,MATCH) + mat->gp21->transition[GP21_INSERT_1_BASE];   
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 4,MATCH,k);    
          }  
        /* From state MATCH to state MATCH */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-1,j-5,MATCH) + mat->gp21->transition[GP21_INSERT_2_BASE];   
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 5,MATCH,k);    
          }  


        /* Ok - finished max calculation for MATCH */ 
        /* Add any movement independant score and put away */ 
         score += CSEQ_GENOMIC_CDSPOT(mat->target,j);    
         GeneWise21_DC_SHADOW_MATRIX(mat,i,j,MATCH) = score; 
        for(k=0;k<7;k++) 
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = localshadow[k];  
        /* Now figure out if any specials need this score */ 
        /* Finished calculating state MATCH */ 


        /* For state INSERT */ 
        /* setting first movement to score */ 
        score = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-3,MATCH) + (mat->query->seg[i]->transition[GW_MATCH2INSERT]+mat->query->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)]);    
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 3,MATCH,k);  
        /* From state INSERT to state INSERT */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-3,INSERT) + (mat->query->seg[i]->transition[GW_INSERT2INSERT]+mat->query->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)]);   
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 3,INSERT,k);   
          }  
        /* From state DELETE to state INSERT */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-3,DELETE) + (mat->query->seg[i]->transition[GW_DELETE2INSERT]+mat->query->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)]);   
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 3,DELETE,k);   
          }  
        /* From state SPACER_0I to state INSERT */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-1,j-6,SPACER_0I) + (((mat->gp21->transition[GP21_SPACER2CDS]+mat->query->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)])+CSEQ_GENOMIC_3SS(mat->target,(j-3)))+mat->query->seg[i]->transition[GW_MATCH_BALANCE_3SS]);     
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 6,SPACER_0I,k);    
          }  
        /* From state SPACER_1I to state INSERT */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-1,j-5,SPACER_1I) + (mat->gp21->transition[GP21_SPACER2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-2)));    
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 5,SPACER_1I,k);    
          }  
        /* From state SPACER_2I to state INSERT */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-1,j-4,SPACER_2I) + (mat->gp21->transition[GP21_SPACER2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-1)));    
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 4,SPACER_2I,k);    
          }  
        /* From state PY_0I to state INSERT */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-1,j-6,PY_0I) + (((mat->gp21->transition[GP21_PY2CDS]+mat->query->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)])+CSEQ_GENOMIC_3SS(mat->target,(j-3)))+mat->query->seg[i]->transition[GW_MATCH_BALANCE_3SS]);     
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 6,PY_0I,k);    
          }  
        /* From state PY_1I to state INSERT */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-1,j-5,PY_1I) + (mat->gp21->transition[GP21_PY2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-2)));    
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 5,PY_1I,k);    
          }  
        /* From state PY_2I to state INSERT */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-1,j-4,PY_2I) + (mat->gp21->transition[GP21_PY2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-1)));    
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 4,PY_2I,k);    
          }  
        /* From state INSERT to state INSERT */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-1,j-2,INSERT) + mat->gp21->transition[GP21_DELETE_1_BASE];  
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 2,INSERT,k);   
          }  
        /* From state INSERT to state INSERT */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-1,j-1,INSERT) + mat->gp21->transition[GP21_DELETE_2_BASE];  
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,INSERT,k);   
          }  
        /* From state INSERT to state INSERT */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-1,j-4,INSERT) + mat->gp21->transition[GP21_INSERT_1_BASE];  
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 4,INSERT,k);   
          }  
        /* From state INSERT to state INSERT */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-1,j-5,INSERT) + mat->gp21->transition[GP21_INSERT_2_BASE];  
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 5,INSERT,k);   
          }  


        /* Ok - finished max calculation for INSERT */ 
        /* Add any movement independant score and put away */ 
         score += CSEQ_GENOMIC_CDSPOT(mat->target,j);    
         GeneWise21_DC_SHADOW_MATRIX(mat,i,j,INSERT) = score;    
        for(k=0;k<7;k++) 
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,k) = localshadow[k]; 
        /* Now figure out if any specials need this score */ 
        /* Finished calculating state INSERT */ 


        /* For state DELETE */ 
        /* setting first movement to score */ 
        score = GeneWise21_DC_SHADOW_MATRIX(mat,i-1,j-0,MATCH) + mat->query->seg[i]->transition[GW_MATCH2DELETE];    
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 0,MATCH,k);  
        /* From state INSERT to state DELETE */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-1,j-0,INSERT) + mat->query->seg[i]->transition[GW_INSERT2DELETE];   
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 0,INSERT,k);   
          }  
        /* From state DELETE to state DELETE */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-1,j-0,DELETE) + mat->query->seg[i]->transition[GW_DELETE2DELETE];   
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 0,DELETE,k);   
          }  


        /* Ok - finished max calculation for DELETE */ 
        /* Add any movement independant score and put away */ 
         GeneWise21_DC_SHADOW_MATRIX(mat,i,j,DELETE) = score;    
        for(k=0;k<7;k++) 
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,k) = localshadow[k]; 
        /* Now figure out if any specials need this score */ 
        /* Finished calculating state DELETE */ 


        /* For state CENTRAL_0M */ 
        /* setting first movement to score */ 
        score = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-8,MATCH) + (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_MATCH2MATCH])+mat->query->seg[i]->transition[GW_MATCH_BALANCE_5SS]);     
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 8,MATCH,k);  
        /* From state INSERT to state CENTRAL_0M */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-8,INSERT) + (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_INSERT2MATCH])+mat->query->seg[i]->transition[GW_INSERT_BALANCE_5SS]);   
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 8,INSERT,k);   
          }  
        /* From state DELETE to state CENTRAL_0M */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-8,DELETE) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_DELETE2MATCH]);   
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 8,DELETE,k);   
          }  
        /* From state CENTRAL_0M to state CENTRAL_0M */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,CENTRAL_0M) + (mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CENTRAL2CENTRAL]);     
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,CENTRAL_0M,k);   
          }  


        /* Ok - finished max calculation for CENTRAL_0M */ 
        /* Add any movement independant score and put away */ 
         GeneWise21_DC_SHADOW_MATRIX(mat,i,j,CENTRAL_0M) = score;    
        for(k=0;k<7;k++) 
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_0M,k) = localshadow[k]; 
        /* Now figure out if any specials need this score */ 
        /* Finished calculating state CENTRAL_0M */ 


        /* For state PY_0M */ 
        /* setting first movement to score */ 
        score = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,CENTRAL_0M) + mat->gp21->transition[GP21_CENTRAL2PY];    
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,CENTRAL_0M,k); 
        /* From state PY_0M to state PY_0M */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,PY_0M) + mat->gp21->transition[GP21_PY2PY];   
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,PY_0M,k);    
          }  


        /* Ok - finished max calculation for PY_0M */ 
        /* Add any movement independant score and put away */ 
         score += mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)];   
         GeneWise21_DC_SHADOW_MATRIX(mat,i,j,PY_0M) = score; 
        for(k=0;k<7;k++) 
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_0M,k) = localshadow[k];  
        /* Now figure out if any specials need this score */ 
        /* Finished calculating state PY_0M */ 


        /* For state SPACER_0M */ 
        /* setting first movement to score */ 
        score = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,PY_0M) + mat->gp21->transition[GP21_PY2SPACER];  
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,PY_0M,k);  
        /* From state SPACER_0M to state SPACER_0M */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,SPACER_0M) + mat->gp21->transition[GP21_SPACER2SPACER];   
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,SPACER_0M,k);    
          }  


        /* Ok - finished max calculation for SPACER_0M */ 
        /* Add any movement independant score and put away */ 
         score += mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)];   
         GeneWise21_DC_SHADOW_MATRIX(mat,i,j,SPACER_0M) = score; 
        for(k=0;k<7;k++) 
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_0M,k) = localshadow[k];  
        /* Now figure out if any specials need this score */ 
        /* Finished calculating state SPACER_0M */ 


        /* For state CENTRAL_1M */ 
        /* setting first movement to score */ 
        score = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-9,MATCH) + (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_MATCH2MATCH])+mat->query->seg[i]->transition[GW_MATCH_BALANCE_5SS]);     
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 9,MATCH,k);  
        /* From state INSERT to state CENTRAL_1M */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-9,INSERT) + (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_INSERT2MATCH])+mat->query->seg[i]->transition[GW_MATCH_BALANCE_5SS]);    
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 9,INSERT,k);   
          }  
        /* From state DELETE to state CENTRAL_1M */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-9,DELETE) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_DELETE2MATCH]);   
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 9,DELETE,k);   
          }  
        /* From state CENTRAL_1M to state CENTRAL_1M */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,CENTRAL_1M) + (mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CENTRAL2CENTRAL]);     
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,CENTRAL_1M,k);   
          }  


        /* Ok - finished max calculation for CENTRAL_1M */ 
        /* Add any movement independant score and put away */ 
         GeneWise21_DC_SHADOW_MATRIX(mat,i,j,CENTRAL_1M) = score;    
        for(k=0;k<7;k++) 
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_1M,k) = localshadow[k]; 
        /* Now figure out if any specials need this score */ 
        /* Finished calculating state CENTRAL_1M */ 


        /* For state PY_1M */ 
        /* setting first movement to score */ 
        score = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,CENTRAL_1M) + mat->gp21->transition[GP21_CENTRAL2PY];    
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,CENTRAL_1M,k); 
        /* From state PY_1M to state PY_1M */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,PY_1M) + mat->gp21->transition[GP21_PY2PY];   
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,PY_1M,k);    
          }  


        /* Ok - finished max calculation for PY_1M */ 
        /* Add any movement independant score and put away */ 
         score += mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)];   
         GeneWise21_DC_SHADOW_MATRIX(mat,i,j,PY_1M) = score; 
        for(k=0;k<7;k++) 
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_1M,k) = localshadow[k];  
        /* Now figure out if any specials need this score */ 
        /* Finished calculating state PY_1M */ 


        /* For state SPACER_1M */ 
        /* setting first movement to score */ 
        score = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,PY_1M) + mat->gp21->transition[GP21_PY2SPACER];  
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,PY_1M,k);  
        /* From state SPACER_1M to state SPACER_1M */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,SPACER_1M) + mat->gp21->transition[GP21_SPACER2SPACER];   
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,SPACER_1M,k);    
          }  


        /* Ok - finished max calculation for SPACER_1M */ 
        /* Add any movement independant score and put away */ 
         score += mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)];   
         GeneWise21_DC_SHADOW_MATRIX(mat,i,j,SPACER_1M) = score; 
        for(k=0;k<7;k++) 
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_1M,k) = localshadow[k];  
        /* Now figure out if any specials need this score */ 
        /* Finished calculating state SPACER_1M */ 


        /* For state CENTRAL_2M */ 
        /* setting first movement to score */ 
        score = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-10,MATCH) + (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_MATCH2MATCH])+mat->query->seg[i]->transition[GW_MATCH_BALANCE_5SS]);    
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 10,MATCH,k); 
        /* From state INSERT to state CENTRAL_2M */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-10,INSERT) + (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_INSERT2MATCH])+mat->query->seg[i]->transition[GW_MATCH_BALANCE_5SS]);   
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 10,INSERT,k);  
          }  
        /* From state DELETE to state CENTRAL_2M */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-10,DELETE) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_DELETE2MATCH]);  
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 10,DELETE,k);  
          }  
        /* From state CENTRAL_2M to state CENTRAL_2M */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,CENTRAL_2M) + (mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CENTRAL2CENTRAL]);     
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,CENTRAL_2M,k);   
          }  


        /* Ok - finished max calculation for CENTRAL_2M */ 
        /* Add any movement independant score and put away */ 
         GeneWise21_DC_SHADOW_MATRIX(mat,i,j,CENTRAL_2M) = score;    
        for(k=0;k<7;k++) 
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_2M,k) = localshadow[k]; 
        /* Now figure out if any specials need this score */ 
        /* Finished calculating state CENTRAL_2M */ 


        /* For state PY_2M */ 
        /* setting first movement to score */ 
        score = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,CENTRAL_2M) + mat->gp21->transition[GP21_CENTRAL2PY];    
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,CENTRAL_2M,k); 
        /* From state PY_2M to state PY_2M */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,PY_2M) + mat->gp21->transition[GP21_PY2PY];   
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,PY_2M,k);    
          }  


        /* Ok - finished max calculation for PY_2M */ 
        /* Add any movement independant score and put away */ 
         score += mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)];   
         GeneWise21_DC_SHADOW_MATRIX(mat,i,j,PY_2M) = score; 
        for(k=0;k<7;k++) 
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_2M,k) = localshadow[k];  
        /* Now figure out if any specials need this score */ 
        /* Finished calculating state PY_2M */ 


        /* For state SPACER_2M */ 
        /* setting first movement to score */ 
        score = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,PY_2M) + mat->gp21->transition[GP21_PY2SPACER];  
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,PY_2M,k);  
        /* From state SPACER_2M to state SPACER_2M */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,SPACER_2M) + mat->gp21->transition[GP21_SPACER2SPACER];   
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,SPACER_2M,k);    
          }  


        /* Ok - finished max calculation for SPACER_2M */ 
        /* Add any movement independant score and put away */ 
         score += mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)];   
         GeneWise21_DC_SHADOW_MATRIX(mat,i,j,SPACER_2M) = score; 
        for(k=0;k<7;k++) 
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_2M,k) = localshadow[k];  
        /* Now figure out if any specials need this score */ 
        /* Finished calculating state SPACER_2M */ 


        /* For state CENTRAL_0I */ 
        /* setting first movement to score */ 
        score = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-8,MATCH) + (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_MATCH2INSERT])+mat->query->seg[i]->transition[GW_MATCH_BALANCE_5SS]);    
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 8,MATCH,k);  
        /* From state INSERT to state CENTRAL_0I */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-8,INSERT) + (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_INSERT2INSERT])+mat->query->seg[i]->transition[GW_INSERT_BALANCE_5SS]);  
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 8,INSERT,k);   
          }  
        /* From state DELETE to state CENTRAL_0I */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-8,DELETE) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_DELETE2INSERT]);  
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 8,DELETE,k);   
          }  
        /* From state CENTRAL_0I to state CENTRAL_0I */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,CENTRAL_0I) + (mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CENTRAL2CENTRAL]);     
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,CENTRAL_0I,k);   
          }  


        /* Ok - finished max calculation for CENTRAL_0I */ 
        /* Add any movement independant score and put away */ 
         GeneWise21_DC_SHADOW_MATRIX(mat,i,j,CENTRAL_0I) = score;    
        for(k=0;k<7;k++) 
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_0I,k) = localshadow[k]; 
        /* Now figure out if any specials need this score */ 
        /* Finished calculating state CENTRAL_0I */ 


        /* For state PY_0I */ 
        /* setting first movement to score */ 
        score = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,CENTRAL_0I) + mat->gp21->transition[GP21_CENTRAL2PY];    
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,CENTRAL_0I,k); 
        /* From state PY_0I to state PY_0I */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,PY_0I) + mat->gp21->transition[GP21_PY2PY];   
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,PY_0I,k);    
          }  


        /* Ok - finished max calculation for PY_0I */ 
        /* Add any movement independant score and put away */ 
         score += mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)];   
         GeneWise21_DC_SHADOW_MATRIX(mat,i,j,PY_0I) = score; 
        for(k=0;k<7;k++) 
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_0I,k) = localshadow[k];  
        /* Now figure out if any specials need this score */ 
        /* Finished calculating state PY_0I */ 


        /* For state SPACER_0I */ 
        /* setting first movement to score */ 
        score = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,PY_0I) + mat->gp21->transition[GP21_PY2SPACER];  
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,PY_0I,k);  
        /* From state SPACER_0I to state SPACER_0I */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,SPACER_0I) + mat->gp21->transition[GP21_SPACER2SPACER];   
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,SPACER_0I,k);    
          }  


        /* Ok - finished max calculation for SPACER_0I */ 
        /* Add any movement independant score and put away */ 
         score += mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)];   
         GeneWise21_DC_SHADOW_MATRIX(mat,i,j,SPACER_0I) = score; 
        for(k=0;k<7;k++) 
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_0I,k) = localshadow[k];  
        /* Now figure out if any specials need this score */ 
        /* Finished calculating state SPACER_0I */ 


        /* For state CENTRAL_1I */ 
        /* setting first movement to score */ 
        score = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-9,MATCH) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_MATCH2INSERT]);   
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 9,MATCH,k);  
        /* From state INSERT to state CENTRAL_1I */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-9,INSERT) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_INSERT2INSERT]);  
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 9,INSERT,k);   
          }  
        /* From state DELETE to state CENTRAL_1I */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-9,DELETE) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_DELETE2INSERT]);  
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 9,DELETE,k);   
          }  
        /* From state CENTRAL_1I to state CENTRAL_1I */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,CENTRAL_1I) + (mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CENTRAL2CENTRAL]);     
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,CENTRAL_1I,k);   
          }  


        /* Ok - finished max calculation for CENTRAL_1I */ 
        /* Add any movement independant score and put away */ 
         GeneWise21_DC_SHADOW_MATRIX(mat,i,j,CENTRAL_1I) = score;    
        for(k=0;k<7;k++) 
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_1I,k) = localshadow[k]; 
        /* Now figure out if any specials need this score */ 
        /* Finished calculating state CENTRAL_1I */ 


        /* For state PY_1I */ 
        /* setting first movement to score */ 
        score = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,CENTRAL_1I) + mat->gp21->transition[GP21_CENTRAL2PY];    
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,CENTRAL_1I,k); 
        /* From state PY_1I to state PY_1I */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,PY_1I) + mat->gp21->transition[GP21_PY2PY];   
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,PY_1I,k);    
          }  


        /* Ok - finished max calculation for PY_1I */ 
        /* Add any movement independant score and put away */ 
         score += mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)];   
         GeneWise21_DC_SHADOW_MATRIX(mat,i,j,PY_1I) = score; 
        for(k=0;k<7;k++) 
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_1I,k) = localshadow[k];  
        /* Now figure out if any specials need this score */ 
        /* Finished calculating state PY_1I */ 


        /* For state SPACER_1I */ 
        /* setting first movement to score */ 
        score = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,PY_1I) + mat->gp21->transition[GP21_PY2SPACER];  
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,PY_1I,k);  
        /* From state SPACER_1I to state SPACER_1I */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,SPACER_1I) + mat->gp21->transition[GP21_SPACER2SPACER];   
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,SPACER_1I,k);    
          }  


        /* Ok - finished max calculation for SPACER_1I */ 
        /* Add any movement independant score and put away */ 
         score += mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)];   
         GeneWise21_DC_SHADOW_MATRIX(mat,i,j,SPACER_1I) = score; 
        for(k=0;k<7;k++) 
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_1I,k) = localshadow[k];  
        /* Now figure out if any specials need this score */ 
        /* Finished calculating state SPACER_1I */ 


        /* For state CENTRAL_2I */ 
        /* setting first movement to score */ 
        score = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-10,MATCH) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_MATCH2INSERT]);  
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 10,MATCH,k); 
        /* From state INSERT to state CENTRAL_2I */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-10,INSERT) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_INSERT2INSERT]);     
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 10,INSERT,k);  
          }  
        /* From state INSERT to state CENTRAL_2I */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-10,INSERT) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_DELETE2INSERT]);     
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 10,INSERT,k);  
          }  
        /* From state CENTRAL_2I to state CENTRAL_2I */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,CENTRAL_2I) + (mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CENTRAL2CENTRAL]);     
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,CENTRAL_2I,k);   
          }  


        /* Ok - finished max calculation for CENTRAL_2I */ 
        /* Add any movement independant score and put away */ 
         GeneWise21_DC_SHADOW_MATRIX(mat,i,j,CENTRAL_2I) = score;    
        for(k=0;k<7;k++) 
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_2I,k) = localshadow[k]; 
        /* Now figure out if any specials need this score */ 
        /* Finished calculating state CENTRAL_2I */ 


        /* For state PY_2I */ 
        /* setting first movement to score */ 
        score = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,CENTRAL_2I) + mat->gp21->transition[GP21_CENTRAL2PY];    
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,CENTRAL_2I,k); 
        /* From state PY_2I to state PY_2I */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,PY_2I) + mat->gp21->transition[GP21_PY2PY];   
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,PY_2I,k);    
          }  


        /* Ok - finished max calculation for PY_2I */ 
        /* Add any movement independant score and put away */ 
         score += mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)];   
         GeneWise21_DC_SHADOW_MATRIX(mat,i,j,PY_2I) = score; 
        for(k=0;k<7;k++) 
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_2I,k) = localshadow[k];  
        /* Now figure out if any specials need this score */ 
        /* Finished calculating state PY_2I */ 


        /* For state SPACER_2I */ 
        /* setting first movement to score */ 
        score = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,PY_2I) + mat->gp21->transition[GP21_PY2SPACER];  
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,PY_2I,k);  
        /* From state SPACER_2I to state SPACER_2I */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,SPACER_2I) + mat->gp21->transition[GP21_SPACER2SPACER];   
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GeneWise21_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,SPACER_2I,k);    
          }  


        /* Ok - finished max calculation for SPACER_2I */ 
        /* Add any movement independant score and put away */ 
         score += mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)];   
         GeneWise21_DC_SHADOW_MATRIX(mat,i,j,SPACER_2I) = score; 
        for(k=0;k<7;k++) 
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_2I,k) = localshadow[k];  
        /* Now figure out if any specials need this score */ 
        /* Finished calculating state SPACER_2I */ 
        } /* end of this is strip */ 
      } /* end of for each valid j column */ 


/* Function:  run_up_dc_GeneWise21(mat,starti,stopi,startj,stopj,dpenv,perc_done)
 *
 * Descrip: No Description
 *
 * Arg:              mat [UNKN ] Undocumented argument [GeneWise21 *]
 * Arg:           starti [UNKN ] Undocumented argument [int]
 * Arg:            stopi [UNKN ] Undocumented argument [int]
 * Arg:           startj [UNKN ] Undocumented argument [int]
 * Arg:            stopj [UNKN ] Undocumented argument [int]
 * Arg:            dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:        perc_done [UNKN ] Undocumented argument [int]
 *
 */
}    
void run_up_dc_GeneWise21(GeneWise21 * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done) 
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
          GeneWise21_DC_SHADOW_MATRIX(mat,i,j,MATCH) = NEGI;     
          GeneWise21_DC_SHADOW_MATRIX(mat,i,j,INSERT) = NEGI;    
          GeneWise21_DC_SHADOW_MATRIX(mat,i,j,DELETE) = NEGI;    
          GeneWise21_DC_SHADOW_MATRIX(mat,i,j,CENTRAL_0M) = NEGI;    
          GeneWise21_DC_SHADOW_MATRIX(mat,i,j,PY_0M) = NEGI;     
          GeneWise21_DC_SHADOW_MATRIX(mat,i,j,SPACER_0M) = NEGI;     
          GeneWise21_DC_SHADOW_MATRIX(mat,i,j,CENTRAL_1M) = NEGI;    
          GeneWise21_DC_SHADOW_MATRIX(mat,i,j,PY_1M) = NEGI;     
          GeneWise21_DC_SHADOW_MATRIX(mat,i,j,SPACER_1M) = NEGI;     
          GeneWise21_DC_SHADOW_MATRIX(mat,i,j,CENTRAL_2M) = NEGI;    
          GeneWise21_DC_SHADOW_MATRIX(mat,i,j,PY_2M) = NEGI;     
          GeneWise21_DC_SHADOW_MATRIX(mat,i,j,SPACER_2M) = NEGI;     
          GeneWise21_DC_SHADOW_MATRIX(mat,i,j,CENTRAL_0I) = NEGI;    
          GeneWise21_DC_SHADOW_MATRIX(mat,i,j,PY_0I) = NEGI;     
          GeneWise21_DC_SHADOW_MATRIX(mat,i,j,SPACER_0I) = NEGI;     
          GeneWise21_DC_SHADOW_MATRIX(mat,i,j,CENTRAL_1I) = NEGI;    
          GeneWise21_DC_SHADOW_MATRIX(mat,i,j,PY_1I) = NEGI;     
          GeneWise21_DC_SHADOW_MATRIX(mat,i,j,SPACER_1I) = NEGI;     
          GeneWise21_DC_SHADOW_MATRIX(mat,i,j,CENTRAL_2I) = NEGI;    
          GeneWise21_DC_SHADOW_MATRIX(mat,i,j,PY_2I) = NEGI;     
          GeneWise21_DC_SHADOW_MATRIX(mat,i,j,SPACER_2I) = NEGI;     
          continue;  
          } /* end of Is not in envelope */ 
        if( num % 1000 == 0 )    
          log_full_error(REPORT,0,"[%d%%%% done]Before mid-j %5d Cells done %d%%%%",perc_done,stopj,(num*100)/total);    


        /* For state MATCH */ 
        /* setting first movement to score */ 
        score = GeneWise21_DC_SHADOW_MATRIX(mat,i-1,j-3,MATCH) + (mat->query->seg[i]->transition[GW_MATCH2MATCH]+mat->query->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]);  
        /* From state INSERT to state MATCH */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-1,j-3,INSERT) + (mat->query->seg[i]->transition[GW_INSERT2MATCH]+mat->query->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state MATCH */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-1,j-3,DELETE) + (mat->query->seg[i]->transition[GW_DELETE2MATCH]+mat->query->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SPACER_0M to state MATCH */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-1,j-6,SPACER_0M) + (((mat->gp21->transition[GP21_SPACER2CDS]+mat->query->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)])+CSEQ_GENOMIC_3SS(mat->target,(j-3)))+mat->query->seg[i]->transition[GW_MATCH_BALANCE_3SS]);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state PY_0M to state MATCH */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-1,j-6,PY_0M) + (((mat->gp21->transition[GP21_PY2CDS]+mat->query->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)])+CSEQ_GENOMIC_3SS(mat->target,(j-3)))+mat->query->seg[i]->transition[GW_MATCH_BALANCE_3SS]);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SPACER_1M to state MATCH */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-1,j-5,SPACER_1M) + (mat->gp21->transition[GP21_SPACER2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-2)));    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state PY_1M to state MATCH */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-1,j-5,PY_1M) + (mat->gp21->transition[GP21_PY2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-2)));    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SPACER_2M to state MATCH */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-1,j-4,SPACER_2M) + (mat->gp21->transition[GP21_SPACER2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-1)));    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state PY_2M to state MATCH */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-1,j-4,PY_2M) + (mat->gp21->transition[GP21_PY2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-1)));    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-1,j-2,MATCH) + mat->gp21->transition[GP21_DELETE_1_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-1,j-1,MATCH) + mat->gp21->transition[GP21_DELETE_2_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-1,j-4,MATCH) + mat->gp21->transition[GP21_INSERT_1_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-1,j-5,MATCH) + mat->gp21->transition[GP21_INSERT_2_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for MATCH */ 
        /* Add any movement independant score and put away */ 
         score += CSEQ_GENOMIC_CDSPOT(mat->target,j);    
         GeneWise21_DC_SHADOW_MATRIX(mat,i,j,MATCH) = score; 
        /* Finished calculating state MATCH */ 


        /* For state INSERT */ 
        /* setting first movement to score */ 
        score = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-3,MATCH) + (mat->query->seg[i]->transition[GW_MATCH2INSERT]+mat->query->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)]);    
        /* From state INSERT to state INSERT */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-3,INSERT) + (mat->query->seg[i]->transition[GW_INSERT2INSERT]+mat->query->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)]);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state INSERT */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-3,DELETE) + (mat->query->seg[i]->transition[GW_DELETE2INSERT]+mat->query->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)]);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SPACER_0I to state INSERT */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-1,j-6,SPACER_0I) + (((mat->gp21->transition[GP21_SPACER2CDS]+mat->query->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)])+CSEQ_GENOMIC_3SS(mat->target,(j-3)))+mat->query->seg[i]->transition[GW_MATCH_BALANCE_3SS]);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SPACER_1I to state INSERT */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-1,j-5,SPACER_1I) + (mat->gp21->transition[GP21_SPACER2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-2)));    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SPACER_2I to state INSERT */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-1,j-4,SPACER_2I) + (mat->gp21->transition[GP21_SPACER2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-1)));    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state PY_0I to state INSERT */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-1,j-6,PY_0I) + (((mat->gp21->transition[GP21_PY2CDS]+mat->query->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)])+CSEQ_GENOMIC_3SS(mat->target,(j-3)))+mat->query->seg[i]->transition[GW_MATCH_BALANCE_3SS]);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state PY_1I to state INSERT */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-1,j-5,PY_1I) + (mat->gp21->transition[GP21_PY2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-2)));    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state PY_2I to state INSERT */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-1,j-4,PY_2I) + (mat->gp21->transition[GP21_PY2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-1)));    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state INSERT */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-1,j-2,INSERT) + mat->gp21->transition[GP21_DELETE_1_BASE];  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state INSERT */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-1,j-1,INSERT) + mat->gp21->transition[GP21_DELETE_2_BASE];  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state INSERT */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-1,j-4,INSERT) + mat->gp21->transition[GP21_INSERT_1_BASE];  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state INSERT */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-1,j-5,INSERT) + mat->gp21->transition[GP21_INSERT_2_BASE];  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for INSERT */ 
        /* Add any movement independant score and put away */ 
         score += CSEQ_GENOMIC_CDSPOT(mat->target,j);    
         GeneWise21_DC_SHADOW_MATRIX(mat,i,j,INSERT) = score;    
        /* Finished calculating state INSERT */ 


        /* For state DELETE */ 
        /* setting first movement to score */ 
        score = GeneWise21_DC_SHADOW_MATRIX(mat,i-1,j-0,MATCH) + mat->query->seg[i]->transition[GW_MATCH2DELETE];    
        /* From state INSERT to state DELETE */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-1,j-0,INSERT) + mat->query->seg[i]->transition[GW_INSERT2DELETE];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state DELETE */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-1,j-0,DELETE) + mat->query->seg[i]->transition[GW_DELETE2DELETE];   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for DELETE */ 
        /* Add any movement independant score and put away */ 
         GeneWise21_DC_SHADOW_MATRIX(mat,i,j,DELETE) = score;    
        /* Finished calculating state DELETE */ 


        /* For state CENTRAL_0M */ 
        /* setting first movement to score */ 
        score = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-8,MATCH) + (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_MATCH2MATCH])+mat->query->seg[i]->transition[GW_MATCH_BALANCE_5SS]);     
        /* From state INSERT to state CENTRAL_0M */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-8,INSERT) + (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_INSERT2MATCH])+mat->query->seg[i]->transition[GW_INSERT_BALANCE_5SS]);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state CENTRAL_0M */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-8,DELETE) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_DELETE2MATCH]);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state CENTRAL_0M to state CENTRAL_0M */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,CENTRAL_0M) + (mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CENTRAL2CENTRAL]);     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for CENTRAL_0M */ 
        /* Add any movement independant score and put away */ 
         GeneWise21_DC_SHADOW_MATRIX(mat,i,j,CENTRAL_0M) = score;    
        /* Finished calculating state CENTRAL_0M */ 


        /* For state PY_0M */ 
        /* setting first movement to score */ 
        score = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,CENTRAL_0M) + mat->gp21->transition[GP21_CENTRAL2PY];    
        /* From state PY_0M to state PY_0M */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,PY_0M) + mat->gp21->transition[GP21_PY2PY];   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for PY_0M */ 
        /* Add any movement independant score and put away */ 
         score += mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)];   
         GeneWise21_DC_SHADOW_MATRIX(mat,i,j,PY_0M) = score; 
        /* Finished calculating state PY_0M */ 


        /* For state SPACER_0M */ 
        /* setting first movement to score */ 
        score = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,PY_0M) + mat->gp21->transition[GP21_PY2SPACER];  
        /* From state SPACER_0M to state SPACER_0M */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,SPACER_0M) + mat->gp21->transition[GP21_SPACER2SPACER];   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for SPACER_0M */ 
        /* Add any movement independant score and put away */ 
         score += mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)];   
         GeneWise21_DC_SHADOW_MATRIX(mat,i,j,SPACER_0M) = score; 
        /* Finished calculating state SPACER_0M */ 


        /* For state CENTRAL_1M */ 
        /* setting first movement to score */ 
        score = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-9,MATCH) + (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_MATCH2MATCH])+mat->query->seg[i]->transition[GW_MATCH_BALANCE_5SS]);     
        /* From state INSERT to state CENTRAL_1M */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-9,INSERT) + (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_INSERT2MATCH])+mat->query->seg[i]->transition[GW_MATCH_BALANCE_5SS]);    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state CENTRAL_1M */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-9,DELETE) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_DELETE2MATCH]);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state CENTRAL_1M to state CENTRAL_1M */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,CENTRAL_1M) + (mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CENTRAL2CENTRAL]);     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for CENTRAL_1M */ 
        /* Add any movement independant score and put away */ 
         GeneWise21_DC_SHADOW_MATRIX(mat,i,j,CENTRAL_1M) = score;    
        /* Finished calculating state CENTRAL_1M */ 


        /* For state PY_1M */ 
        /* setting first movement to score */ 
        score = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,CENTRAL_1M) + mat->gp21->transition[GP21_CENTRAL2PY];    
        /* From state PY_1M to state PY_1M */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,PY_1M) + mat->gp21->transition[GP21_PY2PY];   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for PY_1M */ 
        /* Add any movement independant score and put away */ 
         score += mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)];   
         GeneWise21_DC_SHADOW_MATRIX(mat,i,j,PY_1M) = score; 
        /* Finished calculating state PY_1M */ 


        /* For state SPACER_1M */ 
        /* setting first movement to score */ 
        score = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,PY_1M) + mat->gp21->transition[GP21_PY2SPACER];  
        /* From state SPACER_1M to state SPACER_1M */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,SPACER_1M) + mat->gp21->transition[GP21_SPACER2SPACER];   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for SPACER_1M */ 
        /* Add any movement independant score and put away */ 
         score += mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)];   
         GeneWise21_DC_SHADOW_MATRIX(mat,i,j,SPACER_1M) = score; 
        /* Finished calculating state SPACER_1M */ 


        /* For state CENTRAL_2M */ 
        /* setting first movement to score */ 
        score = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-10,MATCH) + (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_MATCH2MATCH])+mat->query->seg[i]->transition[GW_MATCH_BALANCE_5SS]);    
        /* From state INSERT to state CENTRAL_2M */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-10,INSERT) + (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_INSERT2MATCH])+mat->query->seg[i]->transition[GW_MATCH_BALANCE_5SS]);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state CENTRAL_2M */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-10,DELETE) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_DELETE2MATCH]);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state CENTRAL_2M to state CENTRAL_2M */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,CENTRAL_2M) + (mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CENTRAL2CENTRAL]);     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for CENTRAL_2M */ 
        /* Add any movement independant score and put away */ 
         GeneWise21_DC_SHADOW_MATRIX(mat,i,j,CENTRAL_2M) = score;    
        /* Finished calculating state CENTRAL_2M */ 


        /* For state PY_2M */ 
        /* setting first movement to score */ 
        score = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,CENTRAL_2M) + mat->gp21->transition[GP21_CENTRAL2PY];    
        /* From state PY_2M to state PY_2M */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,PY_2M) + mat->gp21->transition[GP21_PY2PY];   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for PY_2M */ 
        /* Add any movement independant score and put away */ 
         score += mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)];   
         GeneWise21_DC_SHADOW_MATRIX(mat,i,j,PY_2M) = score; 
        /* Finished calculating state PY_2M */ 


        /* For state SPACER_2M */ 
        /* setting first movement to score */ 
        score = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,PY_2M) + mat->gp21->transition[GP21_PY2SPACER];  
        /* From state SPACER_2M to state SPACER_2M */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,SPACER_2M) + mat->gp21->transition[GP21_SPACER2SPACER];   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for SPACER_2M */ 
        /* Add any movement independant score and put away */ 
         score += mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)];   
         GeneWise21_DC_SHADOW_MATRIX(mat,i,j,SPACER_2M) = score; 
        /* Finished calculating state SPACER_2M */ 


        /* For state CENTRAL_0I */ 
        /* setting first movement to score */ 
        score = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-8,MATCH) + (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_MATCH2INSERT])+mat->query->seg[i]->transition[GW_MATCH_BALANCE_5SS]);    
        /* From state INSERT to state CENTRAL_0I */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-8,INSERT) + (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_INSERT2INSERT])+mat->query->seg[i]->transition[GW_INSERT_BALANCE_5SS]);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state CENTRAL_0I */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-8,DELETE) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_DELETE2INSERT]);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state CENTRAL_0I to state CENTRAL_0I */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,CENTRAL_0I) + (mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CENTRAL2CENTRAL]);     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for CENTRAL_0I */ 
        /* Add any movement independant score and put away */ 
         GeneWise21_DC_SHADOW_MATRIX(mat,i,j,CENTRAL_0I) = score;    
        /* Finished calculating state CENTRAL_0I */ 


        /* For state PY_0I */ 
        /* setting first movement to score */ 
        score = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,CENTRAL_0I) + mat->gp21->transition[GP21_CENTRAL2PY];    
        /* From state PY_0I to state PY_0I */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,PY_0I) + mat->gp21->transition[GP21_PY2PY];   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for PY_0I */ 
        /* Add any movement independant score and put away */ 
         score += mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)];   
         GeneWise21_DC_SHADOW_MATRIX(mat,i,j,PY_0I) = score; 
        /* Finished calculating state PY_0I */ 


        /* For state SPACER_0I */ 
        /* setting first movement to score */ 
        score = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,PY_0I) + mat->gp21->transition[GP21_PY2SPACER];  
        /* From state SPACER_0I to state SPACER_0I */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,SPACER_0I) + mat->gp21->transition[GP21_SPACER2SPACER];   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for SPACER_0I */ 
        /* Add any movement independant score and put away */ 
         score += mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)];   
         GeneWise21_DC_SHADOW_MATRIX(mat,i,j,SPACER_0I) = score; 
        /* Finished calculating state SPACER_0I */ 


        /* For state CENTRAL_1I */ 
        /* setting first movement to score */ 
        score = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-9,MATCH) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_MATCH2INSERT]);   
        /* From state INSERT to state CENTRAL_1I */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-9,INSERT) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_INSERT2INSERT]);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state CENTRAL_1I */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-9,DELETE) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_DELETE2INSERT]);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state CENTRAL_1I to state CENTRAL_1I */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,CENTRAL_1I) + (mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CENTRAL2CENTRAL]);     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for CENTRAL_1I */ 
        /* Add any movement independant score and put away */ 
         GeneWise21_DC_SHADOW_MATRIX(mat,i,j,CENTRAL_1I) = score;    
        /* Finished calculating state CENTRAL_1I */ 


        /* For state PY_1I */ 
        /* setting first movement to score */ 
        score = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,CENTRAL_1I) + mat->gp21->transition[GP21_CENTRAL2PY];    
        /* From state PY_1I to state PY_1I */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,PY_1I) + mat->gp21->transition[GP21_PY2PY];   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for PY_1I */ 
        /* Add any movement independant score and put away */ 
         score += mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)];   
         GeneWise21_DC_SHADOW_MATRIX(mat,i,j,PY_1I) = score; 
        /* Finished calculating state PY_1I */ 


        /* For state SPACER_1I */ 
        /* setting first movement to score */ 
        score = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,PY_1I) + mat->gp21->transition[GP21_PY2SPACER];  
        /* From state SPACER_1I to state SPACER_1I */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,SPACER_1I) + mat->gp21->transition[GP21_SPACER2SPACER];   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for SPACER_1I */ 
        /* Add any movement independant score and put away */ 
         score += mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)];   
         GeneWise21_DC_SHADOW_MATRIX(mat,i,j,SPACER_1I) = score; 
        /* Finished calculating state SPACER_1I */ 


        /* For state CENTRAL_2I */ 
        /* setting first movement to score */ 
        score = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-10,MATCH) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_MATCH2INSERT]);  
        /* From state INSERT to state CENTRAL_2I */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-10,INSERT) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_INSERT2INSERT]);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state CENTRAL_2I */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-10,INSERT) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_DELETE2INSERT]);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state CENTRAL_2I to state CENTRAL_2I */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,CENTRAL_2I) + (mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CENTRAL2CENTRAL]);     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for CENTRAL_2I */ 
        /* Add any movement independant score and put away */ 
         GeneWise21_DC_SHADOW_MATRIX(mat,i,j,CENTRAL_2I) = score;    
        /* Finished calculating state CENTRAL_2I */ 


        /* For state PY_2I */ 
        /* setting first movement to score */ 
        score = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,CENTRAL_2I) + mat->gp21->transition[GP21_CENTRAL2PY];    
        /* From state PY_2I to state PY_2I */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,PY_2I) + mat->gp21->transition[GP21_PY2PY];   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for PY_2I */ 
        /* Add any movement independant score and put away */ 
         score += mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)];   
         GeneWise21_DC_SHADOW_MATRIX(mat,i,j,PY_2I) = score; 
        /* Finished calculating state PY_2I */ 


        /* For state SPACER_2I */ 
        /* setting first movement to score */ 
        score = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,PY_2I) + mat->gp21->transition[GP21_PY2SPACER];  
        /* From state SPACER_2I to state SPACER_2I */ 
        temp = GeneWise21_DC_SHADOW_MATRIX(mat,i-0,j-1,SPACER_2I) + mat->gp21->transition[GP21_SPACER2SPACER];   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for SPACER_2I */ 
        /* Add any movement independant score and put away */ 
         score += mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)];   
         GeneWise21_DC_SHADOW_MATRIX(mat,i,j,SPACER_2I) = score; 
        /* Finished calculating state SPACER_2I */ 
        } /* end of this is strip */ 
      } /* end of for each valid j column */ 


/* Function:  init_dc_GeneWise21(mat)
 *
 * Descrip: No Description
 *
 * Arg:        mat [UNKN ] Undocumented argument [GeneWise21 *]
 *
 */
}    
void init_dc_GeneWise21(GeneWise21 * mat) 
{
    register int i;  
    register int j;  
    register int k;  


    for(j=0;j<12;j++)    {  
      for(i=(-1);i<mat->query->len;i++)  {  
        GeneWise21_DC_SHADOW_MATRIX(mat,i,j,MATCH) = NEGI;   
        GeneWise21_DC_SHADOW_MATRIX(mat,i,j,INSERT) = NEGI;  
        GeneWise21_DC_SHADOW_MATRIX(mat,i,j,DELETE) = NEGI;  
        GeneWise21_DC_SHADOW_MATRIX(mat,i,j,CENTRAL_0M) = NEGI;  
        GeneWise21_DC_SHADOW_MATRIX(mat,i,j,PY_0M) = NEGI;   
        GeneWise21_DC_SHADOW_MATRIX(mat,i,j,SPACER_0M) = NEGI;   
        GeneWise21_DC_SHADOW_MATRIX(mat,i,j,CENTRAL_1M) = NEGI;  
        GeneWise21_DC_SHADOW_MATRIX(mat,i,j,PY_1M) = NEGI;   
        GeneWise21_DC_SHADOW_MATRIX(mat,i,j,SPACER_1M) = NEGI;   
        GeneWise21_DC_SHADOW_MATRIX(mat,i,j,CENTRAL_2M) = NEGI;  
        GeneWise21_DC_SHADOW_MATRIX(mat,i,j,PY_2M) = NEGI;   
        GeneWise21_DC_SHADOW_MATRIX(mat,i,j,SPACER_2M) = NEGI;   
        GeneWise21_DC_SHADOW_MATRIX(mat,i,j,CENTRAL_0I) = NEGI;  
        GeneWise21_DC_SHADOW_MATRIX(mat,i,j,PY_0I) = NEGI;   
        GeneWise21_DC_SHADOW_MATRIX(mat,i,j,SPACER_0I) = NEGI;   
        GeneWise21_DC_SHADOW_MATRIX(mat,i,j,CENTRAL_1I) = NEGI;  
        GeneWise21_DC_SHADOW_MATRIX(mat,i,j,PY_1I) = NEGI;   
        GeneWise21_DC_SHADOW_MATRIX(mat,i,j,SPACER_1I) = NEGI;   
        GeneWise21_DC_SHADOW_MATRIX(mat,i,j,CENTRAL_2I) = NEGI;  
        GeneWise21_DC_SHADOW_MATRIX(mat,i,j,PY_2I) = NEGI;   
        GeneWise21_DC_SHADOW_MATRIX(mat,i,j,SPACER_2I) = NEGI;   
        for(k=0;k<7;k++) {  
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = (-1);    
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,k) = (-1);   
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,k) = (-1);   
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_0M,k) = (-1);   
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_0M,k) = (-1);    
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_0M,k) = (-1);    
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_1M,k) = (-1);   
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_1M,k) = (-1);    
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_1M,k) = (-1);    
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_2M,k) = (-1);   
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_2M,k) = (-1);    
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_2M,k) = (-1);    
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_0I,k) = (-1);   
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_0I,k) = (-1);    
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_0I,k) = (-1);    
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_1I,k) = (-1);   
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_1I,k) = (-1);    
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_1I,k) = (-1);    
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_2I,k) = (-1);   
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_2I,k) = (-1);    
          GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_2I,k) = (-1);    
          }  
        }  
      }  


    return;  
}    


/* Function:  start_end_find_end_GeneWise21(mat,endj)
 *
 * Descrip:    First function used to find end of the best path in the special state !end
 *
 *
 * Arg:         mat [UNKN ] Matrix in small mode [GeneWise21 *]
 * Arg:        endj [WRITE] position of end in j (meaningless in i) [int *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int start_end_find_end_GeneWise21(GeneWise21 * mat,int * endj) 
{
    register int j;  
    register int max;    
    register int maxj;   


    max = GeneWise21_DC_SHADOW_SPECIAL(mat,0,mat->target->seq->len-1,END);   
    maxj = mat->target->seq->len-1;  
    for(j= mat->target->seq->len-2 ;j >= 0 ;j--) {  
      if( GeneWise21_DC_SHADOW_SPECIAL(mat,0,j,END) > max )  {  
        max = GeneWise21_DC_SHADOW_SPECIAL(mat,0,j,END); 
        maxj = j;    
        }  
      }  


    if( endj != NULL)    
      *endj = maxj;  


    return max;  
}    


/* Function:  dc_optimised_start_end_calc_GeneWise21(*mat,dpenv)
 *
 * Descrip:    Calculates special strip, leaving start/end/score points in shadow matrix
 *             Works off specially laid out memory from steve searle
 *
 *
 * Arg:         *mat [UNKN ] Undocumented argument [GeneWise21]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean dc_optimised_start_end_calc_GeneWise21(GeneWise21 *mat,DPEnvelope * dpenv) 
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


    score_pointers = (int *) calloc (10 * (leni + 1) * 21,sizeof(int));  
    shadow_pointers = (int *) calloc (10 * (leni + 1) * 21 * 8,sizeof(int)); 


    for(j=0;j<lenj;j++)  { /*for each j strip*/ 
      for(i=0;i<leni;i++)    { /*for each i position in strip*/ 
        num++;   
        if( dpenv != NULL && is_in_DPEnvelope(dpenv,i,j) == FALSE )  { /*Is not in envelope*/ 
          GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i,j,MATCH) = NEGI;     
          GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i,j,INSERT) = NEGI;    
          GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i,j,DELETE) = NEGI;    
          GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i,j,CENTRAL_0M) = NEGI;    
          GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i,j,PY_0M) = NEGI;     
          GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i,j,SPACER_0M) = NEGI;     
          GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i,j,CENTRAL_1M) = NEGI;    
          GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i,j,PY_1M) = NEGI;     
          GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i,j,SPACER_1M) = NEGI;     
          GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i,j,CENTRAL_2M) = NEGI;    
          GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i,j,PY_2M) = NEGI;     
          GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i,j,SPACER_2M) = NEGI;     
          GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i,j,CENTRAL_0I) = NEGI;    
          GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i,j,PY_0I) = NEGI;     
          GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i,j,SPACER_0I) = NEGI;     
          GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i,j,CENTRAL_1I) = NEGI;    
          GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i,j,PY_1I) = NEGI;     
          GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i,j,SPACER_1I) = NEGI;     
          GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i,j,CENTRAL_2I) = NEGI;    
          GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i,j,PY_2I) = NEGI;     
          GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i,j,SPACER_2I) = NEGI;     
          continue;  
          } /* end of Is not in envelope */ 
        if( num%1000 == 0)   
          log_full_error(REPORT,0,"%6d Cells done [%2d%%%%]",num,num*100/total); 




        /* For state MATCH */ 
        /* setting first movement to score */ 
        score = GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i-1,j-3,MATCH) + (mat->query->seg[i]->transition[GW_MATCH2MATCH]+mat->query->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]) + (CSEQ_GENOMIC_CDSPOT(mat->target,j));   
        /* assign local shadown pointer */ 
        localsp = &(GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 3,MATCH,0));    
        /* From state INSERT to state MATCH */ 
        temp = GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i-1,j-3,INSERT) + (mat->query->seg[i]->transition[GW_INSERT2MATCH]+mat->query->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]) +(CSEQ_GENOMIC_CDSPOT(mat->target,j));   
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 3,INSERT,0)); 
          }  
        /* From state DELETE to state MATCH */ 
        temp = GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i-1,j-3,DELETE) + (mat->query->seg[i]->transition[GW_DELETE2MATCH]+mat->query->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]) +(CSEQ_GENOMIC_CDSPOT(mat->target,j));   
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 3,DELETE,0)); 
          }  
        /* From state BEFORE_CODON to state MATCH */ 
        temp = GeneWise21_DC_OPT_SHADOW_SPECIAL(mat,i-1,j-3,BEFORE_CODON) + (mat->query->seg[i]->transition[GW_START2MATCH]+mat->query->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]) + (CSEQ_GENOMIC_CDSPOT(mat->target,j));    
        if( temp  > score )  {  
          score = temp;  
          /* This state [BEFORE_CODON] is a special for MATCH... push top shadow pointers here */ 
          localshadow[0]= i; 
          localshadow[1]= j; 
          localshadow[2]= MATCH; 
          localshadow[3]= (-1);  
          localshadow[4]= (-1);  
          localshadow[5]= (-1);  
          localshadow[6]= score; 
          localsp = localshadow; 
          }  
        /* From state SPACER_0M to state MATCH */ 
        temp = GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i-1,j-6,SPACER_0M) + (((mat->gp21->transition[GP21_SPACER2CDS]+mat->query->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)])+CSEQ_GENOMIC_3SS(mat->target,(j-3)))+mat->query->seg[i]->transition[GW_MATCH_BALANCE_3SS]) +(CSEQ_GENOMIC_CDSPOT(mat->target,j));    
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 6,SPACER_0M,0));  
          }  
        /* From state PY_0M to state MATCH */ 
        temp = GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i-1,j-6,PY_0M) + (((mat->gp21->transition[GP21_PY2CDS]+mat->query->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)])+CSEQ_GENOMIC_3SS(mat->target,(j-3)))+mat->query->seg[i]->transition[GW_MATCH_BALANCE_3SS]) +(CSEQ_GENOMIC_CDSPOT(mat->target,j));    
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 6,PY_0M,0));  
          }  
        /* From state SPACER_1M to state MATCH */ 
        temp = GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i-1,j-5,SPACER_1M) + (mat->gp21->transition[GP21_SPACER2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-2))) +(CSEQ_GENOMIC_CDSPOT(mat->target,j));  
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 5,SPACER_1M,0));  
          }  
        /* From state PY_1M to state MATCH */ 
        temp = GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i-1,j-5,PY_1M) + (mat->gp21->transition[GP21_PY2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-2))) +(CSEQ_GENOMIC_CDSPOT(mat->target,j));  
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 5,PY_1M,0));  
          }  
        /* From state SPACER_2M to state MATCH */ 
        temp = GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i-1,j-4,SPACER_2M) + (mat->gp21->transition[GP21_SPACER2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-1))) +(CSEQ_GENOMIC_CDSPOT(mat->target,j));  
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 4,SPACER_2M,0));  
          }  
        /* From state PY_2M to state MATCH */ 
        temp = GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i-1,j-4,PY_2M) + (mat->gp21->transition[GP21_PY2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-1))) +(CSEQ_GENOMIC_CDSPOT(mat->target,j));  
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 4,PY_2M,0));  
          }  
        /* From state MATCH to state MATCH */ 
        temp = GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i-1,j-2,MATCH) + mat->gp21->transition[GP21_DELETE_1_BASE] +(CSEQ_GENOMIC_CDSPOT(mat->target,j));     
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 2,MATCH,0));  
          }  
        /* From state MATCH to state MATCH */ 
        temp = GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i-1,j-1,MATCH) + mat->gp21->transition[GP21_DELETE_2_BASE] +(CSEQ_GENOMIC_CDSPOT(mat->target,j));     
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 1,MATCH,0));  
          }  
        /* From state MATCH to state MATCH */ 
        temp = GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i-1,j-4,MATCH) + mat->gp21->transition[GP21_INSERT_1_BASE] +(CSEQ_GENOMIC_CDSPOT(mat->target,j));     
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 4,MATCH,0));  
          }  
        /* From state MATCH to state MATCH */ 
        temp = GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i-1,j-5,MATCH) + mat->gp21->transition[GP21_INSERT_2_BASE] +(CSEQ_GENOMIC_CDSPOT(mat->target,j));     
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 5,MATCH,0));  
          }  


        /* Ok - finished max calculation for MATCH */ 
        /* Add any movement independant score and put away */ 
        /* Actually, already done inside scores */ 
         GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i,j,MATCH) = score; 
        for(k=0;k<7;k++) 
          GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = localsp[k];  
        /* Now figure out if any specials need this score */ 


        /* state MATCH is a source for special AFTER_CODON */ 
        temp = score + (mat->query->seg[i]->transition[GW_MATCH2END]) + (0) ;    
        if( temp > GeneWise21_DC_OPT_SHADOW_SPECIAL(mat,i,j,AFTER_CODON) )   {  
          GeneWise21_DC_OPT_SHADOW_SPECIAL(mat,i,j,AFTER_CODON) = temp;  
          /* Have to push only bottem half of system here */ 
          for(k=0;k<3;k++)   
            GeneWise21_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,AFTER_CODON,k) = GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,MATCH,k);    
          GeneWise21_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,AFTER_CODON,6) = GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,MATCH,6);  
          GeneWise21_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,AFTER_CODON,3) = i;    
          GeneWise21_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,AFTER_CODON,4) = j;    
          GeneWise21_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,AFTER_CODON,5) = MATCH;    
          }  




        /* Finished calculating state MATCH */ 


        /* For state INSERT */ 
        /* setting first movement to score */ 
        score = GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i-0,j-3,MATCH) + (mat->query->seg[i]->transition[GW_MATCH2INSERT]+mat->query->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)]) + (CSEQ_GENOMIC_CDSPOT(mat->target,j));     
        /* assign local shadown pointer */ 
        localsp = &(GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 3,MATCH,0));    
        /* From state INSERT to state INSERT */ 
        temp = GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i-0,j-3,INSERT) + (mat->query->seg[i]->transition[GW_INSERT2INSERT]+mat->query->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)]) +(CSEQ_GENOMIC_CDSPOT(mat->target,j));     
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 3,INSERT,0)); 
          }  
        /* From state DELETE to state INSERT */ 
        temp = GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i-0,j-3,DELETE) + (mat->query->seg[i]->transition[GW_DELETE2INSERT]+mat->query->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)]) +(CSEQ_GENOMIC_CDSPOT(mat->target,j));     
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 3,DELETE,0)); 
          }  
        /* From state BEFORE_CODON to state INSERT */ 
        temp = GeneWise21_DC_OPT_SHADOW_SPECIAL(mat,i-0,j-3,BEFORE_CODON) + (mat->query->seg[i]->transition[GW_START2INSERT]+mat->query->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)]) + (CSEQ_GENOMIC_CDSPOT(mat->target,j));  
        if( temp  > score )  {  
          score = temp;  
          /* This state [BEFORE_CODON] is a special for INSERT... push top shadow pointers here */ 
          localshadow[0]= i; 
          localshadow[1]= j; 
          localshadow[2]= INSERT;    
          localshadow[3]= (-1);  
          localshadow[4]= (-1);  
          localshadow[5]= (-1);  
          localshadow[6]= score; 
          localsp = localshadow; 
          }  
        /* From state SPACER_0I to state INSERT */ 
        temp = GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i-1,j-6,SPACER_0I) + (((mat->gp21->transition[GP21_SPACER2CDS]+mat->query->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)])+CSEQ_GENOMIC_3SS(mat->target,(j-3)))+mat->query->seg[i]->transition[GW_MATCH_BALANCE_3SS]) +(CSEQ_GENOMIC_CDSPOT(mat->target,j));   
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 6,SPACER_0I,0));  
          }  
        /* From state SPACER_1I to state INSERT */ 
        temp = GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i-1,j-5,SPACER_1I) + (mat->gp21->transition[GP21_SPACER2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-2))) +(CSEQ_GENOMIC_CDSPOT(mat->target,j));  
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 5,SPACER_1I,0));  
          }  
        /* From state SPACER_2I to state INSERT */ 
        temp = GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i-1,j-4,SPACER_2I) + (mat->gp21->transition[GP21_SPACER2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-1))) +(CSEQ_GENOMIC_CDSPOT(mat->target,j));  
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 4,SPACER_2I,0));  
          }  
        /* From state PY_0I to state INSERT */ 
        temp = GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i-1,j-6,PY_0I) + (((mat->gp21->transition[GP21_PY2CDS]+mat->query->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)])+CSEQ_GENOMIC_3SS(mat->target,(j-3)))+mat->query->seg[i]->transition[GW_MATCH_BALANCE_3SS]) +(CSEQ_GENOMIC_CDSPOT(mat->target,j));   
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 6,PY_0I,0));  
          }  
        /* From state PY_1I to state INSERT */ 
        temp = GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i-1,j-5,PY_1I) + (mat->gp21->transition[GP21_PY2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-2))) +(CSEQ_GENOMIC_CDSPOT(mat->target,j));  
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 5,PY_1I,0));  
          }  
        /* From state PY_2I to state INSERT */ 
        temp = GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i-1,j-4,PY_2I) + (mat->gp21->transition[GP21_PY2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-1))) +(CSEQ_GENOMIC_CDSPOT(mat->target,j));  
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 4,PY_2I,0));  
          }  
        /* From state INSERT to state INSERT */ 
        temp = GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i-1,j-2,INSERT) + mat->gp21->transition[GP21_DELETE_1_BASE] +(CSEQ_GENOMIC_CDSPOT(mat->target,j));    
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 2,INSERT,0)); 
          }  
        /* From state INSERT to state INSERT */ 
        temp = GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i-1,j-1,INSERT) + mat->gp21->transition[GP21_DELETE_2_BASE] +(CSEQ_GENOMIC_CDSPOT(mat->target,j));    
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 1,INSERT,0)); 
          }  
        /* From state INSERT to state INSERT */ 
        temp = GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i-1,j-4,INSERT) + mat->gp21->transition[GP21_INSERT_1_BASE] +(CSEQ_GENOMIC_CDSPOT(mat->target,j));    
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 4,INSERT,0)); 
          }  
        /* From state INSERT to state INSERT */ 
        temp = GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i-1,j-5,INSERT) + mat->gp21->transition[GP21_INSERT_2_BASE] +(CSEQ_GENOMIC_CDSPOT(mat->target,j));    
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 5,INSERT,0)); 
          }  


        /* Ok - finished max calculation for INSERT */ 
        /* Add any movement independant score and put away */ 
        /* Actually, already done inside scores */ 
         GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i,j,INSERT) = score;    
        for(k=0;k<7;k++) 
          GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,INSERT,k) = localsp[k]; 
        /* Now figure out if any specials need this score */ 


        /* state INSERT is a source for special AFTER_CODON */ 
        temp = score + (mat->query->seg[i]->transition[GW_INSERT2END]) + (0) ;   
        if( temp > GeneWise21_DC_OPT_SHADOW_SPECIAL(mat,i,j,AFTER_CODON) )   {  
          GeneWise21_DC_OPT_SHADOW_SPECIAL(mat,i,j,AFTER_CODON) = temp;  
          /* Have to push only bottem half of system here */ 
          for(k=0;k<3;k++)   
            GeneWise21_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,AFTER_CODON,k) = GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,INSERT,k);   
          GeneWise21_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,AFTER_CODON,6) = GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,INSERT,6); 
          GeneWise21_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,AFTER_CODON,3) = i;    
          GeneWise21_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,AFTER_CODON,4) = j;    
          GeneWise21_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,AFTER_CODON,5) = INSERT;   
          }  




        /* Finished calculating state INSERT */ 


        /* For state DELETE */ 
        /* setting first movement to score */ 
        score = GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i-1,j-0,MATCH) + mat->query->seg[i]->transition[GW_MATCH2DELETE] + (0);  
        /* assign local shadown pointer */ 
        localsp = &(GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 0,MATCH,0));    
        /* From state INSERT to state DELETE */ 
        temp = GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i-1,j-0,INSERT) + mat->query->seg[i]->transition[GW_INSERT2DELETE] +(0);  
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 0,INSERT,0)); 
          }  
        /* From state DELETE to state DELETE */ 
        temp = GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i-1,j-0,DELETE) + mat->query->seg[i]->transition[GW_DELETE2DELETE] +(0);  
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 0,DELETE,0)); 
          }  
        /* From state BEFORE_CODON to state DELETE */ 
        temp = GeneWise21_DC_OPT_SHADOW_SPECIAL(mat,i-1,j-0,BEFORE_CODON) + mat->query->seg[i]->transition[GW_START2DELETE] + (0);   
        if( temp  > score )  {  
          score = temp;  
          /* This state [BEFORE_CODON] is a special for DELETE... push top shadow pointers here */ 
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
         GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i,j,DELETE) = score;    
        for(k=0;k<7;k++) 
          GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,DELETE,k) = localsp[k]; 
        /* Now figure out if any specials need this score */ 


        /* state DELETE is a source for special AFTER_CODON */ 
        temp = score + (mat->query->seg[i]->transition[GW_DELETE2END]) + (0) ;   
        if( temp > GeneWise21_DC_OPT_SHADOW_SPECIAL(mat,i,j,AFTER_CODON) )   {  
          GeneWise21_DC_OPT_SHADOW_SPECIAL(mat,i,j,AFTER_CODON) = temp;  
          /* Have to push only bottem half of system here */ 
          for(k=0;k<3;k++)   
            GeneWise21_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,AFTER_CODON,k) = GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,DELETE,k);   
          GeneWise21_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,AFTER_CODON,6) = GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,DELETE,6); 
          GeneWise21_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,AFTER_CODON,3) = i;    
          GeneWise21_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,AFTER_CODON,4) = j;    
          GeneWise21_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,AFTER_CODON,5) = DELETE;   
          }  




        /* Finished calculating state DELETE */ 


        /* For state CENTRAL_0M */ 
        /* setting first movement to score */ 
        score = GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i-0,j-8,MATCH) + (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_MATCH2MATCH])+mat->query->seg[i]->transition[GW_MATCH_BALANCE_5SS]) + (0);   
        /* assign local shadown pointer */ 
        localsp = &(GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 8,MATCH,0));    
        /* From state INSERT to state CENTRAL_0M */ 
        temp = GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i-0,j-8,INSERT) + (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_INSERT2MATCH])+mat->query->seg[i]->transition[GW_INSERT_BALANCE_5SS]) +(0);  
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 8,INSERT,0)); 
          }  
        /* From state DELETE to state CENTRAL_0M */ 
        temp = GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i-0,j-8,DELETE) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_DELETE2MATCH]) +(0);  
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 8,DELETE,0)); 
          }  
        /* From state CENTRAL_0M to state CENTRAL_0M */ 
        temp = GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,CENTRAL_0M) + (mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CENTRAL2CENTRAL]) +(0);    
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,CENTRAL_0M,0)); 
          }  


        /* Ok - finished max calculation for CENTRAL_0M */ 
        /* Add any movement independant score and put away */ 
        /* Actually, already done inside scores */ 
         GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i,j,CENTRAL_0M) = score;    
        for(k=0;k<7;k++) 
          GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_0M,k) = localsp[k]; 
        /* Now figure out if any specials need this score */ 


        /* Finished calculating state CENTRAL_0M */ 


        /* For state PY_0M */ 
        /* setting first movement to score */ 
        score = GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,CENTRAL_0M) + mat->gp21->transition[GP21_CENTRAL2PY] + (mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)]);    
        /* assign local shadown pointer */ 
        localsp = &(GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,CENTRAL_0M,0));   
        /* From state PY_0M to state PY_0M */ 
        temp = GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,PY_0M) + mat->gp21->transition[GP21_PY2PY] +(mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)]);    
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,PY_0M,0));  
          }  


        /* Ok - finished max calculation for PY_0M */ 
        /* Add any movement independant score and put away */ 
        /* Actually, already done inside scores */ 
         GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i,j,PY_0M) = score; 
        for(k=0;k<7;k++) 
          GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,PY_0M,k) = localsp[k];  
        /* Now figure out if any specials need this score */ 


        /* Finished calculating state PY_0M */ 


        /* For state SPACER_0M */ 
        /* setting first movement to score */ 
        score = GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,PY_0M) + mat->gp21->transition[GP21_PY2SPACER] + (mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)]);  
        /* assign local shadown pointer */ 
        localsp = &(GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,PY_0M,0));    
        /* From state SPACER_0M to state SPACER_0M */ 
        temp = GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,SPACER_0M) + mat->gp21->transition[GP21_SPACER2SPACER] +(mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)]);    
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,SPACER_0M,0));  
          }  


        /* Ok - finished max calculation for SPACER_0M */ 
        /* Add any movement independant score and put away */ 
        /* Actually, already done inside scores */ 
         GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i,j,SPACER_0M) = score; 
        for(k=0;k<7;k++) 
          GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,SPACER_0M,k) = localsp[k];  
        /* Now figure out if any specials need this score */ 


        /* Finished calculating state SPACER_0M */ 


        /* For state CENTRAL_1M */ 
        /* setting first movement to score */ 
        score = GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i-0,j-9,MATCH) + (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_MATCH2MATCH])+mat->query->seg[i]->transition[GW_MATCH_BALANCE_5SS]) + (0);   
        /* assign local shadown pointer */ 
        localsp = &(GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 9,MATCH,0));    
        /* From state INSERT to state CENTRAL_1M */ 
        temp = GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i-0,j-9,INSERT) + (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_INSERT2MATCH])+mat->query->seg[i]->transition[GW_MATCH_BALANCE_5SS]) +(0);   
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 9,INSERT,0)); 
          }  
        /* From state DELETE to state CENTRAL_1M */ 
        temp = GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i-0,j-9,DELETE) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_DELETE2MATCH]) +(0);  
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 9,DELETE,0)); 
          }  
        /* From state CENTRAL_1M to state CENTRAL_1M */ 
        temp = GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,CENTRAL_1M) + (mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CENTRAL2CENTRAL]) +(0);    
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,CENTRAL_1M,0)); 
          }  


        /* Ok - finished max calculation for CENTRAL_1M */ 
        /* Add any movement independant score and put away */ 
        /* Actually, already done inside scores */ 
         GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i,j,CENTRAL_1M) = score;    
        for(k=0;k<7;k++) 
          GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_1M,k) = localsp[k]; 
        /* Now figure out if any specials need this score */ 


        /* Finished calculating state CENTRAL_1M */ 


        /* For state PY_1M */ 
        /* setting first movement to score */ 
        score = GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,CENTRAL_1M) + mat->gp21->transition[GP21_CENTRAL2PY] + (mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)]);    
        /* assign local shadown pointer */ 
        localsp = &(GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,CENTRAL_1M,0));   
        /* From state PY_1M to state PY_1M */ 
        temp = GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,PY_1M) + mat->gp21->transition[GP21_PY2PY] +(mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)]);    
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,PY_1M,0));  
          }  


        /* Ok - finished max calculation for PY_1M */ 
        /* Add any movement independant score and put away */ 
        /* Actually, already done inside scores */ 
         GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i,j,PY_1M) = score; 
        for(k=0;k<7;k++) 
          GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,PY_1M,k) = localsp[k];  
        /* Now figure out if any specials need this score */ 


        /* Finished calculating state PY_1M */ 


        /* For state SPACER_1M */ 
        /* setting first movement to score */ 
        score = GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,PY_1M) + mat->gp21->transition[GP21_PY2SPACER] + (mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)]);  
        /* assign local shadown pointer */ 
        localsp = &(GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,PY_1M,0));    
        /* From state SPACER_1M to state SPACER_1M */ 
        temp = GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,SPACER_1M) + mat->gp21->transition[GP21_SPACER2SPACER] +(mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)]);    
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,SPACER_1M,0));  
          }  


        /* Ok - finished max calculation for SPACER_1M */ 
        /* Add any movement independant score and put away */ 
        /* Actually, already done inside scores */ 
         GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i,j,SPACER_1M) = score; 
        for(k=0;k<7;k++) 
          GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,SPACER_1M,k) = localsp[k];  
        /* Now figure out if any specials need this score */ 


        /* Finished calculating state SPACER_1M */ 


        /* For state CENTRAL_2M */ 
        /* setting first movement to score */ 
        score = GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i-0,j-10,MATCH) + (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_MATCH2MATCH])+mat->query->seg[i]->transition[GW_MATCH_BALANCE_5SS]) + (0);  
        /* assign local shadown pointer */ 
        localsp = &(GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 10,MATCH,0));   
        /* From state INSERT to state CENTRAL_2M */ 
        temp = GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i-0,j-10,INSERT) + (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_INSERT2MATCH])+mat->query->seg[i]->transition[GW_MATCH_BALANCE_5SS]) +(0);  
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 10,INSERT,0));    
          }  
        /* From state DELETE to state CENTRAL_2M */ 
        temp = GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i-0,j-10,DELETE) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_DELETE2MATCH]) +(0);     
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 10,DELETE,0));    
          }  
        /* From state CENTRAL_2M to state CENTRAL_2M */ 
        temp = GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,CENTRAL_2M) + (mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CENTRAL2CENTRAL]) +(0);    
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,CENTRAL_2M,0)); 
          }  


        /* Ok - finished max calculation for CENTRAL_2M */ 
        /* Add any movement independant score and put away */ 
        /* Actually, already done inside scores */ 
         GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i,j,CENTRAL_2M) = score;    
        for(k=0;k<7;k++) 
          GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_2M,k) = localsp[k]; 
        /* Now figure out if any specials need this score */ 


        /* Finished calculating state CENTRAL_2M */ 


        /* For state PY_2M */ 
        /* setting first movement to score */ 
        score = GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,CENTRAL_2M) + mat->gp21->transition[GP21_CENTRAL2PY] + (mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)]);    
        /* assign local shadown pointer */ 
        localsp = &(GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,CENTRAL_2M,0));   
        /* From state PY_2M to state PY_2M */ 
        temp = GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,PY_2M) + mat->gp21->transition[GP21_PY2PY] +(mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)]);    
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,PY_2M,0));  
          }  


        /* Ok - finished max calculation for PY_2M */ 
        /* Add any movement independant score and put away */ 
        /* Actually, already done inside scores */ 
         GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i,j,PY_2M) = score; 
        for(k=0;k<7;k++) 
          GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,PY_2M,k) = localsp[k];  
        /* Now figure out if any specials need this score */ 


        /* Finished calculating state PY_2M */ 


        /* For state SPACER_2M */ 
        /* setting first movement to score */ 
        score = GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,PY_2M) + mat->gp21->transition[GP21_PY2SPACER] + (mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)]);  
        /* assign local shadown pointer */ 
        localsp = &(GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,PY_2M,0));    
        /* From state SPACER_2M to state SPACER_2M */ 
        temp = GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,SPACER_2M) + mat->gp21->transition[GP21_SPACER2SPACER] +(mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)]);    
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,SPACER_2M,0));  
          }  


        /* Ok - finished max calculation for SPACER_2M */ 
        /* Add any movement independant score and put away */ 
        /* Actually, already done inside scores */ 
         GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i,j,SPACER_2M) = score; 
        for(k=0;k<7;k++) 
          GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,SPACER_2M,k) = localsp[k];  
        /* Now figure out if any specials need this score */ 


        /* Finished calculating state SPACER_2M */ 


        /* For state CENTRAL_0I */ 
        /* setting first movement to score */ 
        score = GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i-0,j-8,MATCH) + (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_MATCH2INSERT])+mat->query->seg[i]->transition[GW_MATCH_BALANCE_5SS]) + (0);  
        /* assign local shadown pointer */ 
        localsp = &(GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 8,MATCH,0));    
        /* From state INSERT to state CENTRAL_0I */ 
        temp = GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i-0,j-8,INSERT) + (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_INSERT2INSERT])+mat->query->seg[i]->transition[GW_INSERT_BALANCE_5SS]) +(0);     
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 8,INSERT,0)); 
          }  
        /* From state DELETE to state CENTRAL_0I */ 
        temp = GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i-0,j-8,DELETE) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_DELETE2INSERT]) +(0);     
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 8,DELETE,0)); 
          }  
        /* From state CENTRAL_0I to state CENTRAL_0I */ 
        temp = GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,CENTRAL_0I) + (mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CENTRAL2CENTRAL]) +(0);    
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,CENTRAL_0I,0)); 
          }  


        /* Ok - finished max calculation for CENTRAL_0I */ 
        /* Add any movement independant score and put away */ 
        /* Actually, already done inside scores */ 
         GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i,j,CENTRAL_0I) = score;    
        for(k=0;k<7;k++) 
          GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_0I,k) = localsp[k]; 
        /* Now figure out if any specials need this score */ 


        /* Finished calculating state CENTRAL_0I */ 


        /* For state PY_0I */ 
        /* setting first movement to score */ 
        score = GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,CENTRAL_0I) + mat->gp21->transition[GP21_CENTRAL2PY] + (mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)]);    
        /* assign local shadown pointer */ 
        localsp = &(GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,CENTRAL_0I,0));   
        /* From state PY_0I to state PY_0I */ 
        temp = GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,PY_0I) + mat->gp21->transition[GP21_PY2PY] +(mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)]);    
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,PY_0I,0));  
          }  


        /* Ok - finished max calculation for PY_0I */ 
        /* Add any movement independant score and put away */ 
        /* Actually, already done inside scores */ 
         GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i,j,PY_0I) = score; 
        for(k=0;k<7;k++) 
          GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,PY_0I,k) = localsp[k];  
        /* Now figure out if any specials need this score */ 


        /* Finished calculating state PY_0I */ 


        /* For state SPACER_0I */ 
        /* setting first movement to score */ 
        score = GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,PY_0I) + mat->gp21->transition[GP21_PY2SPACER] + (mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)]);  
        /* assign local shadown pointer */ 
        localsp = &(GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,PY_0I,0));    
        /* From state SPACER_0I to state SPACER_0I */ 
        temp = GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,SPACER_0I) + mat->gp21->transition[GP21_SPACER2SPACER] +(mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)]);    
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,SPACER_0I,0));  
          }  


        /* Ok - finished max calculation for SPACER_0I */ 
        /* Add any movement independant score and put away */ 
        /* Actually, already done inside scores */ 
         GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i,j,SPACER_0I) = score; 
        for(k=0;k<7;k++) 
          GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,SPACER_0I,k) = localsp[k];  
        /* Now figure out if any specials need this score */ 


        /* Finished calculating state SPACER_0I */ 


        /* For state CENTRAL_1I */ 
        /* setting first movement to score */ 
        score = GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i-0,j-9,MATCH) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_MATCH2INSERT]) + (0);     
        /* assign local shadown pointer */ 
        localsp = &(GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 9,MATCH,0));    
        /* From state INSERT to state CENTRAL_1I */ 
        temp = GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i-0,j-9,INSERT) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_INSERT2INSERT]) +(0);     
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 9,INSERT,0)); 
          }  
        /* From state DELETE to state CENTRAL_1I */ 
        temp = GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i-0,j-9,DELETE) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_DELETE2INSERT]) +(0);     
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 9,DELETE,0)); 
          }  
        /* From state CENTRAL_1I to state CENTRAL_1I */ 
        temp = GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,CENTRAL_1I) + (mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CENTRAL2CENTRAL]) +(0);    
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,CENTRAL_1I,0)); 
          }  


        /* Ok - finished max calculation for CENTRAL_1I */ 
        /* Add any movement independant score and put away */ 
        /* Actually, already done inside scores */ 
         GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i,j,CENTRAL_1I) = score;    
        for(k=0;k<7;k++) 
          GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_1I,k) = localsp[k]; 
        /* Now figure out if any specials need this score */ 


        /* Finished calculating state CENTRAL_1I */ 


        /* For state PY_1I */ 
        /* setting first movement to score */ 
        score = GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,CENTRAL_1I) + mat->gp21->transition[GP21_CENTRAL2PY] + (mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)]);    
        /* assign local shadown pointer */ 
        localsp = &(GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,CENTRAL_1I,0));   
        /* From state PY_1I to state PY_1I */ 
        temp = GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,PY_1I) + mat->gp21->transition[GP21_PY2PY] +(mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)]);    
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,PY_1I,0));  
          }  


        /* Ok - finished max calculation for PY_1I */ 
        /* Add any movement independant score and put away */ 
        /* Actually, already done inside scores */ 
         GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i,j,PY_1I) = score; 
        for(k=0;k<7;k++) 
          GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,PY_1I,k) = localsp[k];  
        /* Now figure out if any specials need this score */ 


        /* Finished calculating state PY_1I */ 


        /* For state SPACER_1I */ 
        /* setting first movement to score */ 
        score = GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,PY_1I) + mat->gp21->transition[GP21_PY2SPACER] + (mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)]);  
        /* assign local shadown pointer */ 
        localsp = &(GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,PY_1I,0));    
        /* From state SPACER_1I to state SPACER_1I */ 
        temp = GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,SPACER_1I) + mat->gp21->transition[GP21_SPACER2SPACER] +(mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)]);    
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,SPACER_1I,0));  
          }  


        /* Ok - finished max calculation for SPACER_1I */ 
        /* Add any movement independant score and put away */ 
        /* Actually, already done inside scores */ 
         GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i,j,SPACER_1I) = score; 
        for(k=0;k<7;k++) 
          GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,SPACER_1I,k) = localsp[k];  
        /* Now figure out if any specials need this score */ 


        /* Finished calculating state SPACER_1I */ 


        /* For state CENTRAL_2I */ 
        /* setting first movement to score */ 
        score = GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i-0,j-10,MATCH) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_MATCH2INSERT]) + (0);    
        /* assign local shadown pointer */ 
        localsp = &(GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 10,MATCH,0));   
        /* From state INSERT to state CENTRAL_2I */ 
        temp = GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i-0,j-10,INSERT) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_INSERT2INSERT]) +(0);    
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 10,INSERT,0));    
          }  
        /* From state INSERT to state CENTRAL_2I */ 
        temp = GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i-0,j-10,INSERT) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_DELETE2INSERT]) +(0);    
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 10,INSERT,0));    
          }  
        /* From state CENTRAL_2I to state CENTRAL_2I */ 
        temp = GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,CENTRAL_2I) + (mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CENTRAL2CENTRAL]) +(0);    
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,CENTRAL_2I,0)); 
          }  


        /* Ok - finished max calculation for CENTRAL_2I */ 
        /* Add any movement independant score and put away */ 
        /* Actually, already done inside scores */ 
         GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i,j,CENTRAL_2I) = score;    
        for(k=0;k<7;k++) 
          GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_2I,k) = localsp[k]; 
        /* Now figure out if any specials need this score */ 


        /* Finished calculating state CENTRAL_2I */ 


        /* For state PY_2I */ 
        /* setting first movement to score */ 
        score = GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,CENTRAL_2I) + mat->gp21->transition[GP21_CENTRAL2PY] + (mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)]);    
        /* assign local shadown pointer */ 
        localsp = &(GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,CENTRAL_2I,0));   
        /* From state PY_2I to state PY_2I */ 
        temp = GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,PY_2I) + mat->gp21->transition[GP21_PY2PY] +(mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)]);    
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,PY_2I,0));  
          }  


        /* Ok - finished max calculation for PY_2I */ 
        /* Add any movement independant score and put away */ 
        /* Actually, already done inside scores */ 
         GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i,j,PY_2I) = score; 
        for(k=0;k<7;k++) 
          GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,PY_2I,k) = localsp[k];  
        /* Now figure out if any specials need this score */ 


        /* Finished calculating state PY_2I */ 


        /* For state SPACER_2I */ 
        /* setting first movement to score */ 
        score = GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,PY_2I) + mat->gp21->transition[GP21_PY2SPACER] + (mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)]);  
        /* assign local shadown pointer */ 
        localsp = &(GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,PY_2I,0));    
        /* From state SPACER_2I to state SPACER_2I */ 
        temp = GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,SPACER_2I) + mat->gp21->transition[GP21_SPACER2SPACER] +(mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)]);    
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,SPACER_2I,0));  
          }  


        /* Ok - finished max calculation for SPACER_2I */ 
        /* Add any movement independant score and put away */ 
        /* Actually, already done inside scores */ 
         GeneWise21_DC_OPT_SHADOW_MATRIX(mat,i,j,SPACER_2I) = score; 
        for(k=0;k<7;k++) 
          GeneWise21_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,SPACER_2I,k) = localsp[k];  
        /* Now figure out if any specials need this score */ 


        /* Finished calculating state SPACER_2I */ 


        } /* end of for each i position in strip */ 


      /* Special state START has no special to special movements */ 


      /* Special state END has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GeneWise21_DC_OPT_SHADOW_SPECIAL(mat,0,j,END); 


      /* Source AFTER_RND is a special source for END */ 
      temp = GeneWise21_DC_OPT_SHADOW_SPECIAL(mat,0,j - 1,AFTER_RND) + (0) + (0);    
      if( temp > score ) {  
        score = temp;    
        /* Also got to propagate shadows  */ 
        for(k=0;k<7;k++) 
          GeneWise21_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,k) = GeneWise21_DC_OPT_SHADOW_SPECIAL_SP(mat,i - 0,j - 1,AFTER_RND,k); 
        }  


      /* Put back score... (now updated!) */ 
      GeneWise21_DC_OPT_SHADOW_SPECIAL(mat,0,j,END) = score; 
      /* Finished updating state END */ 




      /* Special state BEFORE_RND has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GeneWise21_DC_OPT_SHADOW_SPECIAL(mat,0,j,BEFORE_RND);  


      /* Source START is a special source for BEFORE_RND */ 
      temp = GeneWise21_DC_OPT_SHADOW_SPECIAL(mat,0,j - 1,START) + (mat->rndbase->base[CSEQ_GENOMIC_BASE(mat->target,j)]) + (0);     
      if( temp > score ) {  
        score = temp;    
        /* Also got to propagate shadows  */ 
        for(k=0;k<7;k++) 
          GeneWise21_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,BEFORE_RND,k) = GeneWise21_DC_OPT_SHADOW_SPECIAL_SP(mat,i - 0,j - 1,START,k);  
        }  


      /* Source BEFORE_RND is a special source for BEFORE_RND */ 
      temp = GeneWise21_DC_OPT_SHADOW_SPECIAL(mat,0,j - 1,BEFORE_RND) + ((mat->rndbase->base[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_RND2RND])) + (0);  
      if( temp > score ) {  
        score = temp;    
        /* Also got to propagate shadows  */ 
        for(k=0;k<7;k++) 
          GeneWise21_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,BEFORE_RND,k) = GeneWise21_DC_OPT_SHADOW_SPECIAL_SP(mat,i - 0,j - 1,BEFORE_RND,k); 
        }  


      /* Source BEFORE_CODON is a special source for BEFORE_RND */ 
      temp = GeneWise21_DC_OPT_SHADOW_SPECIAL(mat,0,j - 1,BEFORE_CODON) + ((mat->rndbase->base[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CDS2RND])) + (0);    
      if( temp > score ) {  
        score = temp;    
        /* Also got to propagate shadows  */ 
        for(k=0;k<7;k++) 
          GeneWise21_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,BEFORE_RND,k) = GeneWise21_DC_OPT_SHADOW_SPECIAL_SP(mat,i - 0,j - 1,BEFORE_CODON,k);   
        }  


      /* Put back score... (now updated!) */ 
      GeneWise21_DC_OPT_SHADOW_SPECIAL(mat,0,j,BEFORE_RND) = score;  
      /* Finished updating state BEFORE_RND */ 




      /* Special state BEFORE_CODON has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GeneWise21_DC_OPT_SHADOW_SPECIAL(mat,0,j,BEFORE_CODON);    


      /* Source BEFORE_RND is a special source for BEFORE_CODON */ 
      temp = GeneWise21_DC_OPT_SHADOW_SPECIAL(mat,0,j - 3,BEFORE_RND) + ((mat->rndcodon->codon[CSEQ_GENOMIC_CODON(mat->target,j)]+mat->gp21->transition[GP21_RND2CDS])) + (0);   
      if( temp > score ) {  
        score = temp;    
        /* Also got to propagate shadows  */ 
        for(k=0;k<7;k++) 
          GeneWise21_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,BEFORE_CODON,k) = GeneWise21_DC_OPT_SHADOW_SPECIAL_SP(mat,i - 0,j - 3,BEFORE_RND,k);   
        }  


      /* Source BEFORE_CODON is a special source for BEFORE_CODON */ 
      temp = GeneWise21_DC_OPT_SHADOW_SPECIAL(mat,0,j - 3,BEFORE_CODON) + ((mat->rndcodon->codon[CSEQ_GENOMIC_CODON(mat->target,j)]+mat->gp21->transition[GP21_CDS2CDS])) + (0);     
      if( temp > score ) {  
        score = temp;    
        /* Also got to propagate shadows  */ 
        for(k=0;k<7;k++) 
          GeneWise21_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,BEFORE_CODON,k) = GeneWise21_DC_OPT_SHADOW_SPECIAL_SP(mat,i - 0,j - 3,BEFORE_CODON,k); 
        }  


      /* Source BEFORE_SPACER is a special source for BEFORE_CODON */ 
      temp = GeneWise21_DC_OPT_SHADOW_SPECIAL(mat,0,j - 5,BEFORE_SPACER) + (((mat->gp21->transition[GP21_SPACER2CDS]+mat->rndcodon->codon[CSEQ_GENOMIC_CODON(mat->target,j)])+CSEQ_GENOMIC_3SS(mat->target,(j-3)))) + (0);   
      if( temp > score ) {  
        score = temp;    
        /* Also got to propagate shadows  */ 
        for(k=0;k<7;k++) 
          GeneWise21_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,BEFORE_CODON,k) = GeneWise21_DC_OPT_SHADOW_SPECIAL_SP(mat,i - 0,j - 5,BEFORE_SPACER,k);    
        }  


      /* Source BEFORE_SPACER is a special source for BEFORE_CODON */ 
      temp = GeneWise21_DC_OPT_SHADOW_SPECIAL(mat,0,j - 4,BEFORE_SPACER) + ((mat->gp21->transition[GP21_SPACER2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-2)))) + (0);     
      if( temp > score ) {  
        score = temp;    
        /* Also got to propagate shadows  */ 
        for(k=0;k<7;k++) 
          GeneWise21_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,BEFORE_CODON,k) = GeneWise21_DC_OPT_SHADOW_SPECIAL_SP(mat,i - 0,j - 4,BEFORE_SPACER,k);    
        }  


      /* Source BEFORE_SPACER is a special source for BEFORE_CODON */ 
      temp = GeneWise21_DC_OPT_SHADOW_SPECIAL(mat,0,j - 3,BEFORE_SPACER) + ((mat->gp21->transition[GP21_SPACER2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-1)))) + (0);     
      if( temp > score ) {  
        score = temp;    
        /* Also got to propagate shadows  */ 
        for(k=0;k<7;k++) 
          GeneWise21_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,BEFORE_CODON,k) = GeneWise21_DC_OPT_SHADOW_SPECIAL_SP(mat,i - 0,j - 3,BEFORE_SPACER,k);    
        }  


      /* Put back score... (now updated!) */ 
      GeneWise21_DC_OPT_SHADOW_SPECIAL(mat,0,j,BEFORE_CODON) = score;    
      /* Finished updating state BEFORE_CODON */ 




      /* Special state BEFORE_CENTRAL has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GeneWise21_DC_OPT_SHADOW_SPECIAL(mat,0,j,BEFORE_CENTRAL);  


      /* Source BEFORE_CENTRAL is a special source for BEFORE_CENTRAL */ 
      temp = GeneWise21_DC_OPT_SHADOW_SPECIAL(mat,0,j - 1,BEFORE_CENTRAL) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CENTRAL2CENTRAL])) + (0);  
      if( temp > score ) {  
        score = temp;    
        /* Also got to propagate shadows  */ 
        for(k=0;k<7;k++) 
          GeneWise21_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,BEFORE_CENTRAL,k) = GeneWise21_DC_OPT_SHADOW_SPECIAL_SP(mat,i - 0,j - 1,BEFORE_CENTRAL,k); 
        }  


      /* Source BEFORE_CODON is a special source for BEFORE_CENTRAL */ 
      temp = GeneWise21_DC_OPT_SHADOW_SPECIAL(mat,0,j - 8,BEFORE_CODON) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))) + (0);    
      if( temp > score ) {  
        score = temp;    
        /* Also got to propagate shadows  */ 
        for(k=0;k<7;k++) 
          GeneWise21_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,BEFORE_CENTRAL,k) = GeneWise21_DC_OPT_SHADOW_SPECIAL_SP(mat,i - 0,j - 8,BEFORE_CODON,k);   
        }  


      /* Source BEFORE_CODON is a special source for BEFORE_CENTRAL */ 
      temp = GeneWise21_DC_OPT_SHADOW_SPECIAL(mat,0,j - 9,BEFORE_CODON) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))) + (0);    
      if( temp > score ) {  
        score = temp;    
        /* Also got to propagate shadows  */ 
        for(k=0;k<7;k++) 
          GeneWise21_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,BEFORE_CENTRAL,k) = GeneWise21_DC_OPT_SHADOW_SPECIAL_SP(mat,i - 0,j - 9,BEFORE_CODON,k);   
        }  


      /* Source BEFORE_CODON is a special source for BEFORE_CENTRAL */ 
      temp = GeneWise21_DC_OPT_SHADOW_SPECIAL(mat,0,j - 10,BEFORE_CODON) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))) + (0);   
      if( temp > score ) {  
        score = temp;    
        /* Also got to propagate shadows  */ 
        for(k=0;k<7;k++) 
          GeneWise21_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,BEFORE_CENTRAL,k) = GeneWise21_DC_OPT_SHADOW_SPECIAL_SP(mat,i - 0,j - 10,BEFORE_CODON,k);  
        }  


      /* Put back score... (now updated!) */ 
      GeneWise21_DC_OPT_SHADOW_SPECIAL(mat,0,j,BEFORE_CENTRAL) = score;  
      /* Finished updating state BEFORE_CENTRAL */ 




      /* Special state BEFORE_PY_TRACT has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GeneWise21_DC_OPT_SHADOW_SPECIAL(mat,0,j,BEFORE_PY_TRACT); 


      /* Source BEFORE_CENTRAL is a special source for BEFORE_PY_TRACT */ 
      temp = GeneWise21_DC_OPT_SHADOW_SPECIAL(mat,0,j - 1,BEFORE_CENTRAL) + (mat->gp21->transition[GP21_CENTRAL2PY]) + (mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)]);    
      if( temp > score ) {  
        score = temp;    
        /* Also got to propagate shadows  */ 
        for(k=0;k<7;k++) 
          GeneWise21_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,BEFORE_PY_TRACT,k) = GeneWise21_DC_OPT_SHADOW_SPECIAL_SP(mat,i - 0,j - 1,BEFORE_CENTRAL,k);    
        }  


      /* Source BEFORE_PY_TRACT is a special source for BEFORE_PY_TRACT */ 
      temp = GeneWise21_DC_OPT_SHADOW_SPECIAL(mat,0,j - 1,BEFORE_PY_TRACT) + (mat->gp21->transition[GP21_PY2PY]) + (mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)]);    
      if( temp > score ) {  
        score = temp;    
        /* Also got to propagate shadows  */ 
        for(k=0;k<7;k++) 
          GeneWise21_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,BEFORE_PY_TRACT,k) = GeneWise21_DC_OPT_SHADOW_SPECIAL_SP(mat,i - 0,j - 1,BEFORE_PY_TRACT,k);   
        }  


      /* Put back score... (now updated!) */ 
      GeneWise21_DC_OPT_SHADOW_SPECIAL(mat,0,j,BEFORE_PY_TRACT) = score; 
      /* Finished updating state BEFORE_PY_TRACT */ 




      /* Special state BEFORE_SPACER has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GeneWise21_DC_OPT_SHADOW_SPECIAL(mat,0,j,BEFORE_SPACER);   


      /* Source BEFORE_PY_TRACT is a special source for BEFORE_SPACER */ 
      temp = GeneWise21_DC_OPT_SHADOW_SPECIAL(mat,0,j - 1,BEFORE_PY_TRACT) + (mat->gp21->transition[GP21_PY2SPACER]) + (mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)]);    
      if( temp > score ) {  
        score = temp;    
        /* Also got to propagate shadows  */ 
        for(k=0;k<7;k++) 
          GeneWise21_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,BEFORE_SPACER,k) = GeneWise21_DC_OPT_SHADOW_SPECIAL_SP(mat,i - 0,j - 1,BEFORE_PY_TRACT,k); 
        }  


      /* Source BEFORE_SPACER is a special source for BEFORE_SPACER */ 
      temp = GeneWise21_DC_OPT_SHADOW_SPECIAL(mat,0,j - 1,BEFORE_SPACER) + (mat->gp21->transition[GP21_SPACER2SPACER]) + (mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)]);  
      if( temp > score ) {  
        score = temp;    
        /* Also got to propagate shadows  */ 
        for(k=0;k<7;k++) 
          GeneWise21_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,BEFORE_SPACER,k) = GeneWise21_DC_OPT_SHADOW_SPECIAL_SP(mat,i - 0,j - 1,BEFORE_SPACER,k);   
        }  


      /* Put back score... (now updated!) */ 
      GeneWise21_DC_OPT_SHADOW_SPECIAL(mat,0,j,BEFORE_SPACER) = score;   
      /* Finished updating state BEFORE_SPACER */ 




      /* Special state AFTER_RND has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GeneWise21_DC_OPT_SHADOW_SPECIAL(mat,0,j,AFTER_RND);   


      /* Source AFTER_RND is a special source for AFTER_RND */ 
      temp = GeneWise21_DC_OPT_SHADOW_SPECIAL(mat,0,j - 1,AFTER_RND) + ((mat->rndbase->base[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_RND2RND])) + (0);   
      if( temp > score ) {  
        score = temp;    
        /* Also got to propagate shadows  */ 
        for(k=0;k<7;k++) 
          GeneWise21_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,AFTER_RND,k) = GeneWise21_DC_OPT_SHADOW_SPECIAL_SP(mat,i - 0,j - 1,AFTER_RND,k);   
        }  


      /* Source AFTER_CODON is a special source for AFTER_RND */ 
      temp = GeneWise21_DC_OPT_SHADOW_SPECIAL(mat,0,j - 1,AFTER_CODON) + ((mat->rndbase->base[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CDS2RND])) + (0);     
      if( temp > score ) {  
        score = temp;    
        /* Also got to propagate shadows  */ 
        for(k=0;k<7;k++) 
          GeneWise21_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,AFTER_RND,k) = GeneWise21_DC_OPT_SHADOW_SPECIAL_SP(mat,i - 0,j - 1,AFTER_CODON,k); 
        }  


      /* Put back score... (now updated!) */ 
      GeneWise21_DC_OPT_SHADOW_SPECIAL(mat,0,j,AFTER_RND) = score;   
      /* Finished updating state AFTER_RND */ 




      /* Special state AFTER_CODON has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GeneWise21_DC_OPT_SHADOW_SPECIAL(mat,0,j,AFTER_CODON); 


      /* Source AFTER_RND is a special source for AFTER_CODON */ 
      temp = GeneWise21_DC_OPT_SHADOW_SPECIAL(mat,0,j - 3,AFTER_RND) + ((mat->rndcodon->codon[CSEQ_GENOMIC_CODON(mat->target,j)]+mat->gp21->transition[GP21_RND2CDS])) + (0);    
      if( temp > score ) {  
        score = temp;    
        /* Also got to propagate shadows  */ 
        for(k=0;k<7;k++) 
          GeneWise21_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,AFTER_CODON,k) = GeneWise21_DC_OPT_SHADOW_SPECIAL_SP(mat,i - 0,j - 3,AFTER_RND,k); 
        }  


      /* Source MATCH for state AFTER_CODON is not special... already calculated */ 
      /* Source INSERT for state AFTER_CODON is not special... already calculated */ 
      /* Source DELETE for state AFTER_CODON is not special... already calculated */ 
      /* Source AFTER_CODON is a special source for AFTER_CODON */ 
      temp = GeneWise21_DC_OPT_SHADOW_SPECIAL(mat,0,j - 3,AFTER_CODON) + ((mat->rndcodon->codon[CSEQ_GENOMIC_CODON(mat->target,j)]+mat->gp21->transition[GP21_CDS2CDS])) + (0);  
      if( temp > score ) {  
        score = temp;    
        /* Also got to propagate shadows  */ 
        for(k=0;k<7;k++) 
          GeneWise21_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,AFTER_CODON,k) = GeneWise21_DC_OPT_SHADOW_SPECIAL_SP(mat,i - 0,j - 3,AFTER_CODON,k);   
        }  


      /* Source AFTER_SPACER is a special source for AFTER_CODON */ 
      temp = GeneWise21_DC_OPT_SHADOW_SPECIAL(mat,0,j - 5,AFTER_SPACER) + (((mat->gp21->transition[GP21_SPACER2CDS]+mat->rndcodon->codon[CSEQ_GENOMIC_CODON(mat->target,j)])+CSEQ_GENOMIC_3SS(mat->target,(j-3)))) + (0);    
      if( temp > score ) {  
        score = temp;    
        /* Also got to propagate shadows  */ 
        for(k=0;k<7;k++) 
          GeneWise21_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,AFTER_CODON,k) = GeneWise21_DC_OPT_SHADOW_SPECIAL_SP(mat,i - 0,j - 5,AFTER_SPACER,k);  
        }  


      /* Source AFTER_SPACER is a special source for AFTER_CODON */ 
      temp = GeneWise21_DC_OPT_SHADOW_SPECIAL(mat,0,j - 4,AFTER_SPACER) + ((mat->gp21->transition[GP21_SPACER2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-2)))) + (0);  
      if( temp > score ) {  
        score = temp;    
        /* Also got to propagate shadows  */ 
        for(k=0;k<7;k++) 
          GeneWise21_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,AFTER_CODON,k) = GeneWise21_DC_OPT_SHADOW_SPECIAL_SP(mat,i - 0,j - 4,AFTER_SPACER,k);  
        }  


      /* Source AFTER_SPACER is a special source for AFTER_CODON */ 
      temp = GeneWise21_DC_OPT_SHADOW_SPECIAL(mat,0,j - 3,AFTER_SPACER) + ((mat->gp21->transition[GP21_SPACER2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-1)))) + (0);  
      if( temp > score ) {  
        score = temp;    
        /* Also got to propagate shadows  */ 
        for(k=0;k<7;k++) 
          GeneWise21_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,AFTER_CODON,k) = GeneWise21_DC_OPT_SHADOW_SPECIAL_SP(mat,i - 0,j - 3,AFTER_SPACER,k);  
        }  


      /* Put back score... (now updated!) */ 
      GeneWise21_DC_OPT_SHADOW_SPECIAL(mat,0,j,AFTER_CODON) = score; 
      /* Finished updating state AFTER_CODON */ 




      /* Special state AFTER_CENTRAL has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GeneWise21_DC_OPT_SHADOW_SPECIAL(mat,0,j,AFTER_CENTRAL);   


      /* Source AFTER_CENTRAL is a special source for AFTER_CENTRAL */ 
      temp = GeneWise21_DC_OPT_SHADOW_SPECIAL(mat,0,j - 1,AFTER_CENTRAL) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CENTRAL2CENTRAL])) + (0);   
      if( temp > score ) {  
        score = temp;    
        /* Also got to propagate shadows  */ 
        for(k=0;k<7;k++) 
          GeneWise21_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,AFTER_CENTRAL,k) = GeneWise21_DC_OPT_SHADOW_SPECIAL_SP(mat,i - 0,j - 1,AFTER_CENTRAL,k);   
        }  


      /* Source AFTER_CODON is a special source for AFTER_CENTRAL */ 
      temp = GeneWise21_DC_OPT_SHADOW_SPECIAL(mat,0,j - 8,AFTER_CODON) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))) + (0);     
      if( temp > score ) {  
        score = temp;    
        /* Also got to propagate shadows  */ 
        for(k=0;k<7;k++) 
          GeneWise21_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,AFTER_CENTRAL,k) = GeneWise21_DC_OPT_SHADOW_SPECIAL_SP(mat,i - 0,j - 8,AFTER_CODON,k); 
        }  


      /* Source AFTER_CODON is a special source for AFTER_CENTRAL */ 
      temp = GeneWise21_DC_OPT_SHADOW_SPECIAL(mat,0,j - 9,AFTER_CODON) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))) + (0);     
      if( temp > score ) {  
        score = temp;    
        /* Also got to propagate shadows  */ 
        for(k=0;k<7;k++) 
          GeneWise21_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,AFTER_CENTRAL,k) = GeneWise21_DC_OPT_SHADOW_SPECIAL_SP(mat,i - 0,j - 9,AFTER_CODON,k); 
        }  


      /* Source AFTER_CODON is a special source for AFTER_CENTRAL */ 
      temp = GeneWise21_DC_OPT_SHADOW_SPECIAL(mat,0,j - 10,AFTER_CODON) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))) + (0);    
      if( temp > score ) {  
        score = temp;    
        /* Also got to propagate shadows  */ 
        for(k=0;k<7;k++) 
          GeneWise21_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,AFTER_CENTRAL,k) = GeneWise21_DC_OPT_SHADOW_SPECIAL_SP(mat,i - 0,j - 10,AFTER_CODON,k);    
        }  


      /* Put back score... (now updated!) */ 
      GeneWise21_DC_OPT_SHADOW_SPECIAL(mat,0,j,AFTER_CENTRAL) = score;   
      /* Finished updating state AFTER_CENTRAL */ 




      /* Special state AFTER_PY_TRACT has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GeneWise21_DC_OPT_SHADOW_SPECIAL(mat,0,j,AFTER_PY_TRACT);  


      /* Source AFTER_CENTRAL is a special source for AFTER_PY_TRACT */ 
      temp = GeneWise21_DC_OPT_SHADOW_SPECIAL(mat,0,j - 1,AFTER_CENTRAL) + (mat->gp21->transition[GP21_CENTRAL2PY]) + (mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)]);     
      if( temp > score ) {  
        score = temp;    
        /* Also got to propagate shadows  */ 
        for(k=0;k<7;k++) 
          GeneWise21_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,AFTER_PY_TRACT,k) = GeneWise21_DC_OPT_SHADOW_SPECIAL_SP(mat,i - 0,j - 1,AFTER_CENTRAL,k);  
        }  


      /* Source AFTER_PY_TRACT is a special source for AFTER_PY_TRACT */ 
      temp = GeneWise21_DC_OPT_SHADOW_SPECIAL(mat,0,j - 1,AFTER_PY_TRACT) + (mat->gp21->transition[GP21_PY2PY]) + (mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)]);     
      if( temp > score ) {  
        score = temp;    
        /* Also got to propagate shadows  */ 
        for(k=0;k<7;k++) 
          GeneWise21_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,AFTER_PY_TRACT,k) = GeneWise21_DC_OPT_SHADOW_SPECIAL_SP(mat,i - 0,j - 1,AFTER_PY_TRACT,k); 
        }  


      /* Put back score... (now updated!) */ 
      GeneWise21_DC_OPT_SHADOW_SPECIAL(mat,0,j,AFTER_PY_TRACT) = score;  
      /* Finished updating state AFTER_PY_TRACT */ 




      /* Special state AFTER_SPACER has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GeneWise21_DC_OPT_SHADOW_SPECIAL(mat,0,j,AFTER_SPACER);    


      /* Source AFTER_PY_TRACT is a special source for AFTER_SPACER */ 
      temp = GeneWise21_DC_OPT_SHADOW_SPECIAL(mat,0,j - 1,AFTER_PY_TRACT) + (mat->gp21->transition[GP21_PY2SPACER]) + (mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)]);     
      if( temp > score ) {  
        score = temp;    
        /* Also got to propagate shadows  */ 
        for(k=0;k<7;k++) 
          GeneWise21_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,AFTER_SPACER,k) = GeneWise21_DC_OPT_SHADOW_SPECIAL_SP(mat,i - 0,j - 1,AFTER_PY_TRACT,k);   
        }  


      /* Source AFTER_SPACER is a special source for AFTER_SPACER */ 
      temp = GeneWise21_DC_OPT_SHADOW_SPECIAL(mat,0,j - 1,AFTER_SPACER) + (mat->gp21->transition[GP21_SPACER2SPACER]) + (mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)]);   
      if( temp > score ) {  
        score = temp;    
        /* Also got to propagate shadows  */ 
        for(k=0;k<7;k++) 
          GeneWise21_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,AFTER_SPACER,k) = GeneWise21_DC_OPT_SHADOW_SPECIAL_SP(mat,i - 0,j - 1,AFTER_SPACER,k); 
        }  


      /* Put back score... (now updated!) */ 
      GeneWise21_DC_OPT_SHADOW_SPECIAL(mat,0,j,AFTER_SPACER) = score;    
      /* Finished updating state AFTER_SPACER */ 


      } /* end of for each j strip */ 
    free(score_pointers);    
    free(shadow_pointers);   
    return TRUE;     
}    


/* Function:  init_start_end_linear_GeneWise21(mat)
 *
 * Descrip: No Description
 *
 * Arg:        mat [UNKN ] Undocumented argument [GeneWise21 *]
 *
 */
void init_start_end_linear_GeneWise21(GeneWise21 * mat) 
{
    register int i;  
    register int j;  
    for(j=0;j<12;j++)    {  
      for(i=(-1);i<mat->query->len;i++)  {  
        GeneWise21_DC_SHADOW_MATRIX(mat,i,j,MATCH) = NEGI;   
        GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,0) = (-1);  
        GeneWise21_DC_SHADOW_MATRIX(mat,i,j,INSERT) = NEGI;  
        GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,0) = (-1); 
        GeneWise21_DC_SHADOW_MATRIX(mat,i,j,DELETE) = NEGI;  
        GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,0) = (-1); 
        GeneWise21_DC_SHADOW_MATRIX(mat,i,j,CENTRAL_0M) = NEGI;  
        GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_0M,0) = (-1); 
        GeneWise21_DC_SHADOW_MATRIX(mat,i,j,PY_0M) = NEGI;   
        GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_0M,0) = (-1);  
        GeneWise21_DC_SHADOW_MATRIX(mat,i,j,SPACER_0M) = NEGI;   
        GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_0M,0) = (-1);  
        GeneWise21_DC_SHADOW_MATRIX(mat,i,j,CENTRAL_1M) = NEGI;  
        GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_1M,0) = (-1); 
        GeneWise21_DC_SHADOW_MATRIX(mat,i,j,PY_1M) = NEGI;   
        GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_1M,0) = (-1);  
        GeneWise21_DC_SHADOW_MATRIX(mat,i,j,SPACER_1M) = NEGI;   
        GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_1M,0) = (-1);  
        GeneWise21_DC_SHADOW_MATRIX(mat,i,j,CENTRAL_2M) = NEGI;  
        GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_2M,0) = (-1); 
        GeneWise21_DC_SHADOW_MATRIX(mat,i,j,PY_2M) = NEGI;   
        GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_2M,0) = (-1);  
        GeneWise21_DC_SHADOW_MATRIX(mat,i,j,SPACER_2M) = NEGI;   
        GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_2M,0) = (-1);  
        GeneWise21_DC_SHADOW_MATRIX(mat,i,j,CENTRAL_0I) = NEGI;  
        GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_0I,0) = (-1); 
        GeneWise21_DC_SHADOW_MATRIX(mat,i,j,PY_0I) = NEGI;   
        GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_0I,0) = (-1);  
        GeneWise21_DC_SHADOW_MATRIX(mat,i,j,SPACER_0I) = NEGI;   
        GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_0I,0) = (-1);  
        GeneWise21_DC_SHADOW_MATRIX(mat,i,j,CENTRAL_1I) = NEGI;  
        GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_1I,0) = (-1); 
        GeneWise21_DC_SHADOW_MATRIX(mat,i,j,PY_1I) = NEGI;   
        GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_1I,0) = (-1);  
        GeneWise21_DC_SHADOW_MATRIX(mat,i,j,SPACER_1I) = NEGI;   
        GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_1I,0) = (-1);  
        GeneWise21_DC_SHADOW_MATRIX(mat,i,j,CENTRAL_2I) = NEGI;  
        GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,CENTRAL_2I,0) = (-1); 
        GeneWise21_DC_SHADOW_MATRIX(mat,i,j,PY_2I) = NEGI;   
        GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,PY_2I,0) = (-1);  
        GeneWise21_DC_SHADOW_MATRIX(mat,i,j,SPACER_2I) = NEGI;   
        GeneWise21_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_2I,0) = (-1);  
        }  
      }  


    for(j=(-10);j<mat->target->seq->len;j++) {  
      GeneWise21_DC_SHADOW_SPECIAL(mat,0,j,START) = 0;   
      GeneWise21_DC_SHADOW_SPECIAL_SP(mat,0,j,START,0) = j;  
      GeneWise21_DC_SHADOW_SPECIAL(mat,0,j,END) = NEGI;  
      GeneWise21_DC_SHADOW_SPECIAL_SP(mat,0,j,END,0) = (-1); 
      GeneWise21_DC_SHADOW_SPECIAL(mat,0,j,BEFORE_RND) = NEGI;   
      GeneWise21_DC_SHADOW_SPECIAL_SP(mat,0,j,BEFORE_RND,0) = (-1);  
      GeneWise21_DC_SHADOW_SPECIAL(mat,0,j,BEFORE_CODON) = NEGI; 
      GeneWise21_DC_SHADOW_SPECIAL_SP(mat,0,j,BEFORE_CODON,0) = (-1);    
      GeneWise21_DC_SHADOW_SPECIAL(mat,0,j,BEFORE_CENTRAL) = NEGI;   
      GeneWise21_DC_SHADOW_SPECIAL_SP(mat,0,j,BEFORE_CENTRAL,0) = (-1);  
      GeneWise21_DC_SHADOW_SPECIAL(mat,0,j,BEFORE_PY_TRACT) = NEGI;  
      GeneWise21_DC_SHADOW_SPECIAL_SP(mat,0,j,BEFORE_PY_TRACT,0) = (-1); 
      GeneWise21_DC_SHADOW_SPECIAL(mat,0,j,BEFORE_SPACER) = NEGI;    
      GeneWise21_DC_SHADOW_SPECIAL_SP(mat,0,j,BEFORE_SPACER,0) = (-1);   
      GeneWise21_DC_SHADOW_SPECIAL(mat,0,j,AFTER_RND) = NEGI;    
      GeneWise21_DC_SHADOW_SPECIAL_SP(mat,0,j,AFTER_RND,0) = (-1);   
      GeneWise21_DC_SHADOW_SPECIAL(mat,0,j,AFTER_CODON) = NEGI;  
      GeneWise21_DC_SHADOW_SPECIAL_SP(mat,0,j,AFTER_CODON,0) = (-1); 
      GeneWise21_DC_SHADOW_SPECIAL(mat,0,j,AFTER_CENTRAL) = NEGI;    
      GeneWise21_DC_SHADOW_SPECIAL_SP(mat,0,j,AFTER_CENTRAL,0) = (-1);   
      GeneWise21_DC_SHADOW_SPECIAL(mat,0,j,AFTER_PY_TRACT) = NEGI;   
      GeneWise21_DC_SHADOW_SPECIAL_SP(mat,0,j,AFTER_PY_TRACT,0) = (-1);  
      GeneWise21_DC_SHADOW_SPECIAL(mat,0,j,AFTER_SPACER) = NEGI; 
      GeneWise21_DC_SHADOW_SPECIAL_SP(mat,0,j,AFTER_SPACER,0) = (-1);    
      }  


    return;  
}    


/* Function:  convert_PackAln_to_AlnBlock_GeneWise21(pal)
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
AlnBlock * convert_PackAln_to_AlnBlock_GeneWise21(PackAln * pal) 
{
    AlnConvertSet * acs; 
    AlnBlock * alb;  


    acs = AlnConvertSet_GeneWise21();    
    alb = AlnBlock_from_PackAln(acs,pal);    
    free_AlnConvertSet(acs); 
    return alb;  
}    


 static char * query_label[] = { "MATCH_STATE","INSERT_STATE","DELETE_STATE","INTRON_MATCH_0","INTRON_MATCH_1","INTRON_MATCH_2","INTRON_INSERT_0","INTRON_INSERT_1","INTRON_INSERT_2","END","BEFORE_RND_STATE","AFTER_RND_STATE" };  
/* Function:  AlnConvertSet_GeneWise21(void)
 *
 * Descrip: No Description
 *
 *
 * Return [UNKN ]  Undocumented return value [AlnConvertSet *]
 *
 */
 static char * target_label[] = { "CODON","3SS_PHASE_0","3SS_PHASE_1","3SS_PHASE_2","SEQUENCE_DELETION","SEQUENCE_INSERTION","INSERT","5SS_PHASE_0","CENTRAL_INTRON","PYRIMIDINE_TRACT","SPACER","5SS_PHASE_1","5SS_PHASE_2","END","RANDOM_SEQUENCE" };  
AlnConvertSet * AlnConvertSet_GeneWise21(void) 
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
    acu->state1 = BEFORE_CODON + 21; 
    acu->is_from_special = TRUE; 
    acu->state2 = MATCH;     
    acu->offi = (-1);    
    acu->offj = 3;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[0];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = SPACER_0M; 
    acu->state2 = MATCH;     
    acu->offi = 1;   
    acu->offj = 6;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[1];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = PY_0M; 
    acu->state2 = MATCH;     
    acu->offi = 1;   
    acu->offj = 6;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[1];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = SPACER_1M; 
    acu->state2 = MATCH;     
    acu->offi = 1;   
    acu->offj = 5;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[2];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = PY_1M; 
    acu->state2 = MATCH;     
    acu->offi = 1;   
    acu->offj = 5;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[2];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = SPACER_2M; 
    acu->state2 = MATCH;     
    acu->offi = 1;   
    acu->offj = 4;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[3];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = PY_2M; 
    acu->state2 = MATCH;     
    acu->offi = 1;   
    acu->offj = 4;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[3];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = MATCH; 
    acu->state2 = MATCH;     
    acu->offi = 1;   
    acu->offj = 2;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[4];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = MATCH; 
    acu->state2 = MATCH;     
    acu->offi = 1;   
    acu->offj = 1;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[4];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = MATCH; 
    acu->state2 = MATCH;     
    acu->offi = 1;   
    acu->offj = 4;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[5];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = MATCH; 
    acu->state2 = MATCH;     
    acu->offi = 1;   
    acu->offj = 5;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[5];   
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
    acu->state1 = BEFORE_CODON + 21; 
    acu->is_from_special = TRUE; 
    acu->state2 = INSERT;    
    acu->offi = (-1);    
    acu->offj = 3;   
    acu->label1 = query_label[1];    
    acu->label2 = target_label[0];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = SPACER_0I; 
    acu->state2 = INSERT;    
    acu->offi = 1;   
    acu->offj = 6;   
    acu->label1 = query_label[1];    
    acu->label2 = target_label[1];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = SPACER_1I; 
    acu->state2 = INSERT;    
    acu->offi = 1;   
    acu->offj = 5;   
    acu->label1 = query_label[1];    
    acu->label2 = target_label[2];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = SPACER_2I; 
    acu->state2 = INSERT;    
    acu->offi = 1;   
    acu->offj = 4;   
    acu->label1 = query_label[1];    
    acu->label2 = target_label[3];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = PY_0I; 
    acu->state2 = INSERT;    
    acu->offi = 1;   
    acu->offj = 6;   
    acu->label1 = query_label[1];    
    acu->label2 = target_label[1];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = PY_1I; 
    acu->state2 = INSERT;    
    acu->offi = 1;   
    acu->offj = 5;   
    acu->label1 = query_label[1];    
    acu->label2 = target_label[2];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = PY_2I; 
    acu->state2 = INSERT;    
    acu->offi = 1;   
    acu->offj = 4;   
    acu->label1 = query_label[1];    
    acu->label2 = target_label[3];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = INSERT;    
    acu->state2 = INSERT;    
    acu->offi = 1;   
    acu->offj = 2;   
    acu->label1 = query_label[1];    
    acu->label2 = target_label[4];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = INSERT;    
    acu->state2 = INSERT;    
    acu->offi = 1;   
    acu->offj = 1;   
    acu->label1 = query_label[1];    
    acu->label2 = target_label[4];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = INSERT;    
    acu->state2 = INSERT;    
    acu->offi = 1;   
    acu->offj = 4;   
    acu->label1 = query_label[1];    
    acu->label2 = target_label[5];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = INSERT;    
    acu->state2 = INSERT;    
    acu->offi = 1;   
    acu->offj = 5;   
    acu->label1 = query_label[1];    
    acu->label2 = target_label[5];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = MATCH; 
    acu->state2 = DELETE;    
    acu->offi = 1;   
    acu->offj = 0;   
    acu->label1 = query_label[2];    
    acu->label2 = target_label[6];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = INSERT;    
    acu->state2 = DELETE;    
    acu->offi = 1;   
    acu->offj = 0;   
    acu->label1 = query_label[2];    
    acu->label2 = target_label[6];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = DELETE;    
    acu->state2 = DELETE;    
    acu->offi = 1;   
    acu->offj = 0;   
    acu->label1 = query_label[2];    
    acu->label2 = target_label[6];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = BEFORE_CODON + 21; 
    acu->is_from_special = TRUE; 
    acu->state2 = DELETE;    
    acu->offi = (-1);    
    acu->offj = 0;   
    acu->label1 = query_label[2];    
    acu->label2 = target_label[6];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = MATCH; 
    acu->state2 = CENTRAL_0M;    
    acu->offi = 0;   
    acu->offj = 8;   
    acu->label1 = query_label[3];    
    acu->label2 = target_label[7];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = INSERT;    
    acu->state2 = CENTRAL_0M;    
    acu->offi = 0;   
    acu->offj = 8;   
    acu->label1 = query_label[3];    
    acu->label2 = target_label[7];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = DELETE;    
    acu->state2 = CENTRAL_0M;    
    acu->offi = 0;   
    acu->offj = 8;   
    acu->label1 = query_label[3];    
    acu->label2 = target_label[7];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = CENTRAL_0M;    
    acu->state2 = CENTRAL_0M;    
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[3];    
    acu->label2 = target_label[8];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = CENTRAL_0M;    
    acu->state2 = PY_0M;     
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[3];    
    acu->label2 = target_label[9];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = PY_0M; 
    acu->state2 = PY_0M;     
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[3];    
    acu->label2 = target_label[9];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = PY_0M; 
    acu->state2 = SPACER_0M;     
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[3];    
    acu->label2 = target_label[10];  
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = SPACER_0M; 
    acu->state2 = SPACER_0M;     
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[3];    
    acu->label2 = target_label[10];  
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = MATCH; 
    acu->state2 = CENTRAL_1M;    
    acu->offi = 0;   
    acu->offj = 9;   
    acu->label1 = query_label[4];    
    acu->label2 = target_label[11];  
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = INSERT;    
    acu->state2 = CENTRAL_1M;    
    acu->offi = 0;   
    acu->offj = 9;   
    acu->label1 = query_label[4];    
    acu->label2 = target_label[11];  
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = DELETE;    
    acu->state2 = CENTRAL_1M;    
    acu->offi = 0;   
    acu->offj = 9;   
    acu->label1 = query_label[4];    
    acu->label2 = target_label[11];  
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = CENTRAL_1M;    
    acu->state2 = CENTRAL_1M;    
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[4];    
    acu->label2 = target_label[8];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = CENTRAL_1M;    
    acu->state2 = PY_1M;     
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[4];    
    acu->label2 = target_label[9];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = PY_1M; 
    acu->state2 = PY_1M;     
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[4];    
    acu->label2 = target_label[9];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = PY_1M; 
    acu->state2 = SPACER_1M;     
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[4];    
    acu->label2 = target_label[10];  
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = SPACER_1M; 
    acu->state2 = SPACER_1M;     
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[4];    
    acu->label2 = target_label[10];  
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = MATCH; 
    acu->state2 = CENTRAL_2M;    
    acu->offi = 0;   
    acu->offj = 10;  
    acu->label1 = query_label[5];    
    acu->label2 = target_label[12];  
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = INSERT;    
    acu->state2 = CENTRAL_2M;    
    acu->offi = 0;   
    acu->offj = 10;  
    acu->label1 = query_label[5];    
    acu->label2 = target_label[12];  
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = DELETE;    
    acu->state2 = CENTRAL_2M;    
    acu->offi = 0;   
    acu->offj = 10;  
    acu->label1 = query_label[5];    
    acu->label2 = target_label[12];  
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = CENTRAL_2M;    
    acu->state2 = CENTRAL_2M;    
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[5];    
    acu->label2 = target_label[8];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = CENTRAL_2M;    
    acu->state2 = PY_2M;     
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[5];    
    acu->label2 = target_label[9];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = PY_2M; 
    acu->state2 = PY_2M;     
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[5];    
    acu->label2 = target_label[9];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = PY_2M; 
    acu->state2 = SPACER_2M;     
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[5];    
    acu->label2 = target_label[10];  
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = SPACER_2M; 
    acu->state2 = SPACER_2M;     
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[5];    
    acu->label2 = target_label[10];  
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = MATCH; 
    acu->state2 = CENTRAL_0I;    
    acu->offi = 0;   
    acu->offj = 8;   
    acu->label1 = query_label[6];    
    acu->label2 = target_label[7];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = INSERT;    
    acu->state2 = CENTRAL_0I;    
    acu->offi = 0;   
    acu->offj = 8;   
    acu->label1 = query_label[6];    
    acu->label2 = target_label[7];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = DELETE;    
    acu->state2 = CENTRAL_0I;    
    acu->offi = 0;   
    acu->offj = 8;   
    acu->label1 = query_label[6];    
    acu->label2 = target_label[7];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = CENTRAL_0I;    
    acu->state2 = CENTRAL_0I;    
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[6];    
    acu->label2 = target_label[8];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = CENTRAL_0I;    
    acu->state2 = PY_0I;     
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[6];    
    acu->label2 = target_label[9];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = PY_0I; 
    acu->state2 = PY_0I;     
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[6];    
    acu->label2 = target_label[9];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = PY_0I; 
    acu->state2 = SPACER_0I;     
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[6];    
    acu->label2 = target_label[10];  
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = SPACER_0I; 
    acu->state2 = SPACER_0I;     
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[6];    
    acu->label2 = target_label[10];  
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = MATCH; 
    acu->state2 = CENTRAL_1I;    
    acu->offi = 0;   
    acu->offj = 9;   
    acu->label1 = query_label[7];    
    acu->label2 = target_label[11];  
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = INSERT;    
    acu->state2 = CENTRAL_1I;    
    acu->offi = 0;   
    acu->offj = 9;   
    acu->label1 = query_label[7];    
    acu->label2 = target_label[11];  
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = DELETE;    
    acu->state2 = CENTRAL_1I;    
    acu->offi = 0;   
    acu->offj = 9;   
    acu->label1 = query_label[7];    
    acu->label2 = target_label[11];  
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = CENTRAL_1I;    
    acu->state2 = CENTRAL_1I;    
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[7];    
    acu->label2 = target_label[8];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = CENTRAL_1I;    
    acu->state2 = PY_1I;     
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[7];    
    acu->label2 = target_label[9];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = PY_1I; 
    acu->state2 = PY_1I;     
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[7];    
    acu->label2 = target_label[9];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = PY_1I; 
    acu->state2 = SPACER_1I;     
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[7];    
    acu->label2 = target_label[10];  
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = SPACER_1I; 
    acu->state2 = SPACER_1I;     
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[7];    
    acu->label2 = target_label[10];  
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = MATCH; 
    acu->state2 = CENTRAL_2I;    
    acu->offi = 0;   
    acu->offj = 10;  
    acu->label1 = query_label[8];    
    acu->label2 = target_label[12];  
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = INSERT;    
    acu->state2 = CENTRAL_2I;    
    acu->offi = 0;   
    acu->offj = 10;  
    acu->label1 = query_label[8];    
    acu->label2 = target_label[12];  
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = INSERT;    
    acu->state2 = CENTRAL_2I;    
    acu->offi = 0;   
    acu->offj = 10;  
    acu->label1 = query_label[8];    
    acu->label2 = target_label[12];  
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = CENTRAL_2I;    
    acu->state2 = CENTRAL_2I;    
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[8];    
    acu->label2 = target_label[8];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = CENTRAL_2I;    
    acu->state2 = PY_2I;     
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[8];    
    acu->label2 = target_label[9];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = PY_2I; 
    acu->state2 = PY_2I;     
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[8];    
    acu->label2 = target_label[9];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = PY_2I; 
    acu->state2 = SPACER_2I;     
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[8];    
    acu->label2 = target_label[10];  
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = SPACER_2I; 
    acu->state2 = SPACER_2I;     
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[8];    
    acu->label2 = target_label[10];  
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = AFTER_RND + 21;    
    acu->state2 = END + 21;  
    acu->offi = (-1);    
    acu->offj = 1;   
    acu->label1 = query_label[9];    
    acu->label2 = target_label[13];  
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = START + 21;    
    acu->state2 = BEFORE_RND + 21;   
    acu->offi = (-1);    
    acu->offj = 1;   
    acu->label1 = query_label[10];   
    acu->label2 = target_label[14];  
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = BEFORE_RND + 21;   
    acu->state2 = BEFORE_RND + 21;   
    acu->offi = (-1);    
    acu->offj = 1;   
    acu->label1 = query_label[10];   
    acu->label2 = target_label[14];  
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = BEFORE_CODON + 21; 
    acu->state2 = BEFORE_RND + 21;   
    acu->offi = (-1);    
    acu->offj = 1;   
    acu->label1 = query_label[10];   
    acu->label2 = target_label[14];  
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = BEFORE_RND + 21;   
    acu->state2 = BEFORE_CODON + 21;     
    acu->offi = (-1);    
    acu->offj = 3;   
    acu->label1 = query_label[10];   
    acu->label2 = target_label[0];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = BEFORE_CODON + 21; 
    acu->state2 = BEFORE_CODON + 21;     
    acu->offi = (-1);    
    acu->offj = 3;   
    acu->label1 = query_label[10];   
    acu->label2 = target_label[0];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = BEFORE_SPACER + 21;    
    acu->state2 = BEFORE_CODON + 21;     
    acu->offi = (-1);    
    acu->offj = 5;   
    acu->label1 = query_label[10];   
    acu->label2 = target_label[1];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = BEFORE_SPACER + 21;    
    acu->state2 = BEFORE_CODON + 21;     
    acu->offi = (-1);    
    acu->offj = 4;   
    acu->label1 = query_label[10];   
    acu->label2 = target_label[2];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = BEFORE_SPACER + 21;    
    acu->state2 = BEFORE_CODON + 21;     
    acu->offi = (-1);    
    acu->offj = 3;   
    acu->label1 = query_label[10];   
    acu->label2 = target_label[3];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = BEFORE_CENTRAL + 21;   
    acu->state2 = BEFORE_CENTRAL + 21;   
    acu->offi = (-1);    
    acu->offj = 1;   
    acu->label1 = query_label[10];   
    acu->label2 = target_label[8];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = BEFORE_CODON + 21; 
    acu->state2 = BEFORE_CENTRAL + 21;   
    acu->offi = (-1);    
    acu->offj = 8;   
    acu->label1 = query_label[10];   
    acu->label2 = target_label[7];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = BEFORE_CODON + 21; 
    acu->state2 = BEFORE_CENTRAL + 21;   
    acu->offi = (-1);    
    acu->offj = 9;   
    acu->label1 = query_label[10];   
    acu->label2 = target_label[11];  
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = BEFORE_CODON + 21; 
    acu->state2 = BEFORE_CENTRAL + 21;   
    acu->offi = (-1);    
    acu->offj = 10;  
    acu->label1 = query_label[10];   
    acu->label2 = target_label[12];  
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = BEFORE_CENTRAL + 21;   
    acu->state2 = BEFORE_PY_TRACT + 21;  
    acu->offi = (-1);    
    acu->offj = 1;   
    acu->label1 = query_label[10];   
    acu->label2 = target_label[9];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = BEFORE_PY_TRACT + 21;  
    acu->state2 = BEFORE_PY_TRACT + 21;  
    acu->offi = (-1);    
    acu->offj = 1;   
    acu->label1 = query_label[10];   
    acu->label2 = target_label[9];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = BEFORE_PY_TRACT + 21;  
    acu->state2 = BEFORE_SPACER + 21;    
    acu->offi = (-1);    
    acu->offj = 1;   
    acu->label1 = query_label[10];   
    acu->label2 = target_label[10];  
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = BEFORE_SPACER + 21;    
    acu->state2 = BEFORE_SPACER + 21;    
    acu->offi = (-1);    
    acu->offj = 1;   
    acu->label1 = query_label[10];   
    acu->label2 = target_label[10];  
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = AFTER_RND + 21;    
    acu->state2 = AFTER_RND + 21;    
    acu->offi = (-1);    
    acu->offj = 1;   
    acu->label1 = query_label[11];   
    acu->label2 = target_label[14];  
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = AFTER_CODON + 21;  
    acu->state2 = AFTER_RND + 21;    
    acu->offi = (-1);    
    acu->offj = 1;   
    acu->label1 = query_label[11];   
    acu->label2 = target_label[14];  
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = AFTER_RND + 21;    
    acu->state2 = AFTER_CODON + 21;  
    acu->offi = (-1);    
    acu->offj = 3;   
    acu->label1 = query_label[11];   
    acu->label2 = target_label[0];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = MATCH; 
    acu->state2 = AFTER_CODON + 21;  
    acu->offi = (-1);    
    acu->offj = 0;   
    acu->label1 = query_label[11];   
    acu->label2 = target_label[0];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = INSERT;    
    acu->state2 = AFTER_CODON + 21;  
    acu->offi = (-1);    
    acu->offj = 0;   
    acu->label1 = query_label[11];   
    acu->label2 = target_label[0];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = DELETE;    
    acu->state2 = AFTER_CODON + 21;  
    acu->offi = (-1);    
    acu->offj = 0;   
    acu->label1 = query_label[11];   
    acu->label2 = target_label[0];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = AFTER_CODON + 21;  
    acu->state2 = AFTER_CODON + 21;  
    acu->offi = (-1);    
    acu->offj = 3;   
    acu->label1 = query_label[11];   
    acu->label2 = target_label[0];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = AFTER_SPACER + 21; 
    acu->state2 = AFTER_CODON + 21;  
    acu->offi = (-1);    
    acu->offj = 5;   
    acu->label1 = query_label[11];   
    acu->label2 = target_label[1];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = AFTER_SPACER + 21; 
    acu->state2 = AFTER_CODON + 21;  
    acu->offi = (-1);    
    acu->offj = 4;   
    acu->label1 = query_label[11];   
    acu->label2 = target_label[2];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = AFTER_SPACER + 21; 
    acu->state2 = AFTER_CODON + 21;  
    acu->offi = (-1);    
    acu->offj = 3;   
    acu->label1 = query_label[11];   
    acu->label2 = target_label[3];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = AFTER_CENTRAL + 21;    
    acu->state2 = AFTER_CENTRAL + 21;    
    acu->offi = (-1);    
    acu->offj = 1;   
    acu->label1 = query_label[11];   
    acu->label2 = target_label[8];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = AFTER_CODON + 21;  
    acu->state2 = AFTER_CENTRAL + 21;    
    acu->offi = (-1);    
    acu->offj = 8;   
    acu->label1 = query_label[11];   
    acu->label2 = target_label[7];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = AFTER_CODON + 21;  
    acu->state2 = AFTER_CENTRAL + 21;    
    acu->offi = (-1);    
    acu->offj = 9;   
    acu->label1 = query_label[11];   
    acu->label2 = target_label[11];  
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = AFTER_CODON + 21;  
    acu->state2 = AFTER_CENTRAL + 21;    
    acu->offi = (-1);    
    acu->offj = 10;  
    acu->label1 = query_label[11];   
    acu->label2 = target_label[12];  
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = AFTER_CENTRAL + 21;    
    acu->state2 = AFTER_PY_TRACT + 21;   
    acu->offi = (-1);    
    acu->offj = 1;   
    acu->label1 = query_label[11];   
    acu->label2 = target_label[9];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = AFTER_PY_TRACT + 21;   
    acu->state2 = AFTER_PY_TRACT + 21;   
    acu->offi = (-1);    
    acu->offj = 1;   
    acu->label1 = query_label[11];   
    acu->label2 = target_label[9];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = AFTER_PY_TRACT + 21;   
    acu->state2 = AFTER_SPACER + 21;     
    acu->offi = (-1);    
    acu->offj = 1;   
    acu->label1 = query_label[11];   
    acu->label2 = target_label[10];  
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = AFTER_SPACER + 21; 
    acu->state2 = AFTER_SPACER + 21;     
    acu->offi = (-1);    
    acu->offj = 1;   
    acu->label1 = query_label[11];   
    acu->label2 = target_label[10];  
    add_collapse_label_AlnConvertSet(out,"BEFORE_RND_STATE","RANDOM_SEQUENCE");  
    add_collapse_label_AlnConvertSet(out,"AFTER_RND_STATE","RANDOM_SEQUENCE");   
    add_collapse_label_AlnConvertSet(out,"BEFORE_RND_STATE","CENTRAL_INTRON");   
    add_collapse_label_AlnConvertSet(out,"AFTER_RND_STATE","CENTRAL_INTRON");    
    add_collapse_label_AlnConvertSet(out,"INTRON_MATCH_0","CENTRAL_INTRON"); 
    add_collapse_label_AlnConvertSet(out,"INTRON_MATCH_1","CENTRAL_INTRON"); 
    add_collapse_label_AlnConvertSet(out,"INTRON_MATCH_2","CENTRAL_INTRON"); 
    add_collapse_label_AlnConvertSet(out,"INTRON_INSERT_0","CENTRAL_INTRON");    
    add_collapse_label_AlnConvertSet(out,"INTRON_INSERT_1","CENTRAL_INTRON");    
    add_collapse_label_AlnConvertSet(out,"INTRON_INSERT_2","CENTRAL_INTRON");    
    return out;  
}    


/* Function:  PackAln_read_Expl_GeneWise21(mat)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [GeneWise21 *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * PackAln_read_Expl_GeneWise21(GeneWise21 * mat) 
{
    GeneWise21_access_func_holder holder;    


    holder.access_main    = GeneWise21_explicit_access_main; 
    holder.access_special = GeneWise21_explicit_access_special;  
    return PackAln_read_generic_GeneWise21(mat,holder);  
}    


/* Function:  GeneWise21_explicit_access_main(mat,i,j,state)
 *
 * Descrip: No Description
 *
 * Arg:          mat [UNKN ] Undocumented argument [GeneWise21 *]
 * Arg:            i [UNKN ] Undocumented argument [int]
 * Arg:            j [UNKN ] Undocumented argument [int]
 * Arg:        state [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int GeneWise21_explicit_access_main(GeneWise21 * mat,int i,int j,int state) 
{
    return GeneWise21_EXPL_MATRIX(mat,i,j,state);    
}    


/* Function:  GeneWise21_explicit_access_special(mat,i,j,state)
 *
 * Descrip: No Description
 *
 * Arg:          mat [UNKN ] Undocumented argument [GeneWise21 *]
 * Arg:            i [UNKN ] Undocumented argument [int]
 * Arg:            j [UNKN ] Undocumented argument [int]
 * Arg:        state [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int GeneWise21_explicit_access_special(GeneWise21 * mat,int i,int j,int state) 
{
    return GeneWise21_EXPL_SPECIAL(mat,i,j,state);   
}    


/* Function:  PackAln_read_generic_GeneWise21(mat,h)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [GeneWise21 *]
 * Arg:          h [UNKN ] Undocumented argument [GeneWise21_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * PackAln_read_generic_GeneWise21(GeneWise21 * mat,GeneWise21_access_func_holder h) 
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


    out->score =  find_end_GeneWise21(mat,&i,&j,&state,&isspecial,h);    


    /* Add final end transition (at the moment we have not got the score! */ 
    if( (pau= PackAlnUnit_alloc()) == NULL  || add_PackAln(out,pau) == FALSE )   {  
      warn("Failed the first PackAlnUnit alloc, %d length of Alignment in GeneWise21_basic_read, returning a mess.(Sorry!)",out->len);   
      return out;    
      }  


    /* Put in positions for end trans. Remember that coordinates in C style */ 
    pau->i = i;  
    pau->j = j;  
    if( isspecial != TRUE)   
      pau->state = state;    
    else pau->state = state + 21;    
    prev=pau;    
    while( state != START || isspecial != TRUE)  { /*while state != START*/ 


      if( isspecial == TRUE )    
        max_calc_special_GeneWise21(mat,i,j,state,isspecial,&i,&j,&state,&isspecial,&cellscore,h);   
      else   
        max_calc_GeneWise21(mat,i,j,state,isspecial,&i,&j,&state,&isspecial,&cellscore,h);   
      if(i == GeneWise21_READ_OFF_ERROR || j == GeneWise21_READ_OFF_ERROR || state == GeneWise21_READ_OFF_ERROR )    {  
        warn("Problem - hit bad read off system, exiting now");  
        break;   
        }  
      if( (pau= PackAlnUnit_alloc()) == NULL  || add_PackAln(out,pau) == FALSE ) {  
        warn("Failed a PackAlnUnit alloc, %d length of Alignment in GeneWise21_basic_read, returning partial alignment",out->len);   
        break;   
        }  


      /* Put in positions for block. Remember that coordinates in C style */ 
      pau->i = i;    
      pau->j = j;    
      if( isspecial != TRUE)     
        pau->state = state;  
      else pau->state = state + 21;  
      prev->score = cellscore;   
      prev = pau;    
      } /* end of while state != START */ 


    invert_PackAln(out); 
    return out;  
}    


/* Function:  find_end_GeneWise21(mat,ri,rj,state,isspecial,h)
 *
 * Descrip: No Description
 *
 * Arg:              mat [UNKN ] Undocumented argument [GeneWise21 *]
 * Arg:               ri [UNKN ] Undocumented argument [int *]
 * Arg:               rj [UNKN ] Undocumented argument [int *]
 * Arg:            state [UNKN ] Undocumented argument [int *]
 * Arg:        isspecial [UNKN ] Undocumented argument [boolean *]
 * Arg:                h [UNKN ] Undocumented argument [GeneWise21_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int find_end_GeneWise21(GeneWise21 * mat,int * ri,int * rj,int * state,boolean * isspecial,GeneWise21_access_func_holder h) 
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


/* Function:  GeneWise21_debug_show_matrix(mat,starti,stopi,startj,stopj,ofp)
 *
 * Descrip: No Description
 *
 * Arg:           mat [UNKN ] Undocumented argument [GeneWise21 *]
 * Arg:        starti [UNKN ] Undocumented argument [int]
 * Arg:         stopi [UNKN ] Undocumented argument [int]
 * Arg:        startj [UNKN ] Undocumented argument [int]
 * Arg:         stopj [UNKN ] Undocumented argument [int]
 * Arg:           ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
void GeneWise21_debug_show_matrix(GeneWise21 * mat,int starti,int stopi,int startj,int stopj,FILE * ofp) 
{
    register int i;  
    register int j;  


    for(i=starti;i<stopi && i < mat->query->len;i++) {  
      for(j=startj;j<stopj && j < mat->target->seq->len;j++) {  
        fprintf(ofp,"Cell [%d - %d]\n",i,j);     
        fprintf(ofp,"State MATCH %d\n",GeneWise21_EXPL_MATRIX(mat,i,j,MATCH));   
        fprintf(ofp,"State INSERT %d\n",GeneWise21_EXPL_MATRIX(mat,i,j,INSERT)); 
        fprintf(ofp,"State DELETE %d\n",GeneWise21_EXPL_MATRIX(mat,i,j,DELETE)); 
        fprintf(ofp,"State CENTRAL_0M %d\n",GeneWise21_EXPL_MATRIX(mat,i,j,CENTRAL_0M)); 
        fprintf(ofp,"State PY_0M %d\n",GeneWise21_EXPL_MATRIX(mat,i,j,PY_0M));   
        fprintf(ofp,"State SPACER_0M %d\n",GeneWise21_EXPL_MATRIX(mat,i,j,SPACER_0M));   
        fprintf(ofp,"State CENTRAL_1M %d\n",GeneWise21_EXPL_MATRIX(mat,i,j,CENTRAL_1M)); 
        fprintf(ofp,"State PY_1M %d\n",GeneWise21_EXPL_MATRIX(mat,i,j,PY_1M));   
        fprintf(ofp,"State SPACER_1M %d\n",GeneWise21_EXPL_MATRIX(mat,i,j,SPACER_1M));   
        fprintf(ofp,"State CENTRAL_2M %d\n",GeneWise21_EXPL_MATRIX(mat,i,j,CENTRAL_2M)); 
        fprintf(ofp,"State PY_2M %d\n",GeneWise21_EXPL_MATRIX(mat,i,j,PY_2M));   
        fprintf(ofp,"State SPACER_2M %d\n",GeneWise21_EXPL_MATRIX(mat,i,j,SPACER_2M));   
        fprintf(ofp,"State CENTRAL_0I %d\n",GeneWise21_EXPL_MATRIX(mat,i,j,CENTRAL_0I)); 
        fprintf(ofp,"State PY_0I %d\n",GeneWise21_EXPL_MATRIX(mat,i,j,PY_0I));   
        fprintf(ofp,"State SPACER_0I %d\n",GeneWise21_EXPL_MATRIX(mat,i,j,SPACER_0I));   
        fprintf(ofp,"State CENTRAL_1I %d\n",GeneWise21_EXPL_MATRIX(mat,i,j,CENTRAL_1I)); 
        fprintf(ofp,"State PY_1I %d\n",GeneWise21_EXPL_MATRIX(mat,i,j,PY_1I));   
        fprintf(ofp,"State SPACER_1I %d\n",GeneWise21_EXPL_MATRIX(mat,i,j,SPACER_1I));   
        fprintf(ofp,"State CENTRAL_2I %d\n",GeneWise21_EXPL_MATRIX(mat,i,j,CENTRAL_2I)); 
        fprintf(ofp,"State PY_2I %d\n",GeneWise21_EXPL_MATRIX(mat,i,j,PY_2I));   
        fprintf(ofp,"State SPACER_2I %d\n",GeneWise21_EXPL_MATRIX(mat,i,j,SPACER_2I));   
        fprintf(ofp,"\n\n"); 
        }  
      }  


}    


/* Function:  max_calc_GeneWise21(mat,i,j,state,isspecial,reti,retj,retstate,retspecial,cellscore,h)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [GeneWise21 *]
 * Arg:                 i [UNKN ] Undocumented argument [int]
 * Arg:                 j [UNKN ] Undocumented argument [int]
 * Arg:             state [UNKN ] Undocumented argument [int]
 * Arg:         isspecial [UNKN ] Undocumented argument [boolean]
 * Arg:              reti [UNKN ] Undocumented argument [int *]
 * Arg:              retj [UNKN ] Undocumented argument [int *]
 * Arg:          retstate [UNKN ] Undocumented argument [int *]
 * Arg:        retspecial [UNKN ] Undocumented argument [boolean *]
 * Arg:         cellscore [UNKN ] Undocumented argument [int *]
 * Arg:                 h [UNKN ] Undocumented argument [GeneWise21_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int max_calc_GeneWise21(GeneWise21 * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,GeneWise21_access_func_holder h) 
{
    register int temp;   
    register int cscore; 


    *reti = (*retj) = (*retstate) = GeneWise21_READ_OFF_ERROR;   


    if( i < 0 || j < 0 || i > mat->query->len || j > mat->target->seq->len)  {  
      warn("In GeneWise21 matrix special read off - out of bounds on matrix [i,j is %d,%d state %d in standard matrix]",i,j,state);  
      return -1;     
      }  


    /* Then you have to select the correct switch statement to figure out the readoff      */ 
    /* Somewhat odd - reverse the order of calculation and return as soon as it is correct */ 
    cscore = (*h.access_main)(mat,i,j,state);    
    switch(state)    { /*Switch state */ 
      case MATCH :   
        temp = cscore - (mat->gp21->transition[GP21_INSERT_2_BASE]) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j)); 
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
        temp = cscore - (mat->gp21->transition[GP21_INSERT_1_BASE]) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j)); 
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
        temp = cscore - (mat->gp21->transition[GP21_DELETE_2_BASE]) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j)); 
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
        temp = cscore - (mat->gp21->transition[GP21_DELETE_1_BASE]) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j)); 
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
        temp = cscore - ((mat->gp21->transition[GP21_PY2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-1)))) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j));  
        if( temp == (*h.access_main)(mat,i - 1,j - 4,PY_2M) )    {  
          *reti = i - 1; 
          *retj = j - 4; 
          *retstate = PY_2M; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-4,PY_2M);   
            }  
          return (*h.access_main)(mat,i - 1,j - 4,PY_2M);    
          }  
        temp = cscore - ((mat->gp21->transition[GP21_SPACER2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-1)))) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j));  
        if( temp == (*h.access_main)(mat,i - 1,j - 4,SPACER_2M) )    {  
          *reti = i - 1; 
          *retj = j - 4; 
          *retstate = SPACER_2M; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-4,SPACER_2M);   
            }  
          return (*h.access_main)(mat,i - 1,j - 4,SPACER_2M);    
          }  
        temp = cscore - ((mat->gp21->transition[GP21_PY2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-2)))) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j));  
        if( temp == (*h.access_main)(mat,i - 1,j - 5,PY_1M) )    {  
          *reti = i - 1; 
          *retj = j - 5; 
          *retstate = PY_1M; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-5,PY_1M);   
            }  
          return (*h.access_main)(mat,i - 1,j - 5,PY_1M);    
          }  
        temp = cscore - ((mat->gp21->transition[GP21_SPACER2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-2)))) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j));  
        if( temp == (*h.access_main)(mat,i - 1,j - 5,SPACER_1M) )    {  
          *reti = i - 1; 
          *retj = j - 5; 
          *retstate = SPACER_1M; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-5,SPACER_1M);   
            }  
          return (*h.access_main)(mat,i - 1,j - 5,SPACER_1M);    
          }  
        temp = cscore - ((((mat->gp21->transition[GP21_PY2CDS]+mat->query->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)])+CSEQ_GENOMIC_3SS(mat->target,(j-3)))+mat->query->seg[i]->transition[GW_MATCH_BALANCE_3SS])) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j));    
        if( temp == (*h.access_main)(mat,i - 1,j - 6,PY_0M) )    {  
          *reti = i - 1; 
          *retj = j - 6; 
          *retstate = PY_0M; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-6,PY_0M);   
            }  
          return (*h.access_main)(mat,i - 1,j - 6,PY_0M);    
          }  
        temp = cscore - ((((mat->gp21->transition[GP21_SPACER2CDS]+mat->query->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)])+CSEQ_GENOMIC_3SS(mat->target,(j-3)))+mat->query->seg[i]->transition[GW_MATCH_BALANCE_3SS])) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j));    
        if( temp == (*h.access_main)(mat,i - 1,j - 6,SPACER_0M) )    {  
          *reti = i - 1; 
          *retj = j - 6; 
          *retstate = SPACER_0M; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-6,SPACER_0M);   
            }  
          return (*h.access_main)(mat,i - 1,j - 6,SPACER_0M);    
          }  
        temp = cscore - ((mat->query->seg[i]->transition[GW_START2MATCH]+mat->query->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)])) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j)); 
        if( temp == (*h.access_special)(mat,i - 1,j - 3,BEFORE_CODON) )  {  
          *reti = i - 1; 
          *retj = j - 3; 
          *retstate = BEFORE_CODON;  
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_special)(mat,i-1,j-3,BEFORE_CODON); 
            }  
          return (*h.access_main)(mat,i - 1,j - 3,BEFORE_CODON);     
          }  
        temp = cscore - ((mat->query->seg[i]->transition[GW_DELETE2MATCH]+mat->query->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)])) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j));    
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
        temp = cscore - ((mat->query->seg[i]->transition[GW_INSERT2MATCH]+mat->query->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)])) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j));    
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
        temp = cscore - ((mat->query->seg[i]->transition[GW_MATCH2MATCH]+mat->query->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)])) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j)); 
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
        warn("Major problem (!) - in GeneWise21 read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      case INSERT :  
        temp = cscore - (mat->gp21->transition[GP21_INSERT_2_BASE]) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j)); 
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
        temp = cscore - (mat->gp21->transition[GP21_INSERT_1_BASE]) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j)); 
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
        temp = cscore - (mat->gp21->transition[GP21_DELETE_2_BASE]) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j)); 
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
        temp = cscore - (mat->gp21->transition[GP21_DELETE_1_BASE]) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j)); 
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
        temp = cscore - ((mat->gp21->transition[GP21_PY2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-1)))) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j));  
        if( temp == (*h.access_main)(mat,i - 1,j - 4,PY_2I) )    {  
          *reti = i - 1; 
          *retj = j - 4; 
          *retstate = PY_2I; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-4,PY_2I);   
            }  
          return (*h.access_main)(mat,i - 1,j - 4,PY_2I);    
          }  
        temp = cscore - ((mat->gp21->transition[GP21_PY2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-2)))) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j));  
        if( temp == (*h.access_main)(mat,i - 1,j - 5,PY_1I) )    {  
          *reti = i - 1; 
          *retj = j - 5; 
          *retstate = PY_1I; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-5,PY_1I);   
            }  
          return (*h.access_main)(mat,i - 1,j - 5,PY_1I);    
          }  
        temp = cscore - ((((mat->gp21->transition[GP21_PY2CDS]+mat->query->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)])+CSEQ_GENOMIC_3SS(mat->target,(j-3)))+mat->query->seg[i]->transition[GW_MATCH_BALANCE_3SS])) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j));   
        if( temp == (*h.access_main)(mat,i - 1,j - 6,PY_0I) )    {  
          *reti = i - 1; 
          *retj = j - 6; 
          *retstate = PY_0I; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-6,PY_0I);   
            }  
          return (*h.access_main)(mat,i - 1,j - 6,PY_0I);    
          }  
        temp = cscore - ((mat->gp21->transition[GP21_SPACER2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-1)))) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j));  
        if( temp == (*h.access_main)(mat,i - 1,j - 4,SPACER_2I) )    {  
          *reti = i - 1; 
          *retj = j - 4; 
          *retstate = SPACER_2I; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-4,SPACER_2I);   
            }  
          return (*h.access_main)(mat,i - 1,j - 4,SPACER_2I);    
          }  
        temp = cscore - ((mat->gp21->transition[GP21_SPACER2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-2)))) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j));  
        if( temp == (*h.access_main)(mat,i - 1,j - 5,SPACER_1I) )    {  
          *reti = i - 1; 
          *retj = j - 5; 
          *retstate = SPACER_1I; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-5,SPACER_1I);   
            }  
          return (*h.access_main)(mat,i - 1,j - 5,SPACER_1I);    
          }  
        temp = cscore - ((((mat->gp21->transition[GP21_SPACER2CDS]+mat->query->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)])+CSEQ_GENOMIC_3SS(mat->target,(j-3)))+mat->query->seg[i]->transition[GW_MATCH_BALANCE_3SS])) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j));   
        if( temp == (*h.access_main)(mat,i - 1,j - 6,SPACER_0I) )    {  
          *reti = i - 1; 
          *retj = j - 6; 
          *retstate = SPACER_0I; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-6,SPACER_0I);   
            }  
          return (*h.access_main)(mat,i - 1,j - 6,SPACER_0I);    
          }  
        temp = cscore - ((mat->query->seg[i]->transition[GW_START2INSERT]+mat->query->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)])) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j));   
        if( temp == (*h.access_special)(mat,i - 0,j - 3,BEFORE_CODON) )  {  
          *reti = i - 0; 
          *retj = j - 3; 
          *retstate = BEFORE_CODON;  
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_special)(mat,i-0,j-3,BEFORE_CODON); 
            }  
          return (*h.access_main)(mat,i - 0,j - 3,BEFORE_CODON);     
          }  
        temp = cscore - ((mat->query->seg[i]->transition[GW_DELETE2INSERT]+mat->query->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)])) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j));  
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
        temp = cscore - ((mat->query->seg[i]->transition[GW_INSERT2INSERT]+mat->query->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)])) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j));  
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
        temp = cscore - ((mat->query->seg[i]->transition[GW_MATCH2INSERT]+mat->query->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)])) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j));   
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
        warn("Major problem (!) - in GeneWise21 read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      case DELETE :  
        temp = cscore - (mat->query->seg[i]->transition[GW_START2DELETE]) -  (0);    
        if( temp == (*h.access_special)(mat,i - 1,j - 0,BEFORE_CODON) )  {  
          *reti = i - 1; 
          *retj = j - 0; 
          *retstate = BEFORE_CODON;  
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_special)(mat,i-1,j-0,BEFORE_CODON); 
            }  
          return (*h.access_main)(mat,i - 1,j - 0,BEFORE_CODON);     
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
        warn("Major problem (!) - in GeneWise21 read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      case CENTRAL_0M :  
        temp = cscore - ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CENTRAL2CENTRAL])) -  (0); 
        if( temp == (*h.access_main)(mat,i - 0,j - 1,CENTRAL_0M) )   {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = CENTRAL_0M;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,CENTRAL_0M);  
            }  
          return (*h.access_main)(mat,i - 0,j - 1,CENTRAL_0M);   
          }  
        temp = cscore - (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_DELETE2MATCH])) -  (0);   
        if( temp == (*h.access_main)(mat,i - 0,j - 8,DELETE) )   {  
          *reti = i - 0; 
          *retj = j - 8; 
          *retstate = DELETE;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-8,DELETE);  
            }  
          return (*h.access_main)(mat,i - 0,j - 8,DELETE);   
          }  
        temp = cscore - ((((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_INSERT2MATCH])+mat->query->seg[i]->transition[GW_INSERT_BALANCE_5SS])) -  (0);   
        if( temp == (*h.access_main)(mat,i - 0,j - 8,INSERT) )   {  
          *reti = i - 0; 
          *retj = j - 8; 
          *retstate = INSERT;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-8,INSERT);  
            }  
          return (*h.access_main)(mat,i - 0,j - 8,INSERT);   
          }  
        temp = cscore - ((((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_MATCH2MATCH])+mat->query->seg[i]->transition[GW_MATCH_BALANCE_5SS])) -  (0); 
        if( temp == (*h.access_main)(mat,i - 0,j - 8,MATCH) )    {  
          *reti = i - 0; 
          *retj = j - 8; 
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-8,MATCH);   
            }  
          return (*h.access_main)(mat,i - 0,j - 8,MATCH);    
          }  
        warn("Major problem (!) - in GeneWise21 read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      case PY_0M :   
        temp = cscore - (mat->gp21->transition[GP21_PY2PY]) -  (mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)]);    
        if( temp == (*h.access_main)(mat,i - 0,j - 1,PY_0M) )    {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = PY_0M; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,PY_0M);   
            }  
          return (*h.access_main)(mat,i - 0,j - 1,PY_0M);    
          }  
        temp = cscore - (mat->gp21->transition[GP21_CENTRAL2PY]) -  (mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)]);   
        if( temp == (*h.access_main)(mat,i - 0,j - 1,CENTRAL_0M) )   {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = CENTRAL_0M;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,CENTRAL_0M);  
            }  
          return (*h.access_main)(mat,i - 0,j - 1,CENTRAL_0M);   
          }  
        warn("Major problem (!) - in GeneWise21 read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      case SPACER_0M :   
        temp = cscore - (mat->gp21->transition[GP21_SPACER2SPACER]) -  (mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)]);    
        if( temp == (*h.access_main)(mat,i - 0,j - 1,SPACER_0M) )    {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = SPACER_0M; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,SPACER_0M);   
            }  
          return (*h.access_main)(mat,i - 0,j - 1,SPACER_0M);    
          }  
        temp = cscore - (mat->gp21->transition[GP21_PY2SPACER]) -  (mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)]);    
        if( temp == (*h.access_main)(mat,i - 0,j - 1,PY_0M) )    {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = PY_0M; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,PY_0M);   
            }  
          return (*h.access_main)(mat,i - 0,j - 1,PY_0M);    
          }  
        warn("Major problem (!) - in GeneWise21 read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      case CENTRAL_1M :  
        temp = cscore - ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CENTRAL2CENTRAL])) -  (0); 
        if( temp == (*h.access_main)(mat,i - 0,j - 1,CENTRAL_1M) )   {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = CENTRAL_1M;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,CENTRAL_1M);  
            }  
          return (*h.access_main)(mat,i - 0,j - 1,CENTRAL_1M);   
          }  
        temp = cscore - (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_DELETE2MATCH])) -  (0);   
        if( temp == (*h.access_main)(mat,i - 0,j - 9,DELETE) )   {  
          *reti = i - 0; 
          *retj = j - 9; 
          *retstate = DELETE;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-9,DELETE);  
            }  
          return (*h.access_main)(mat,i - 0,j - 9,DELETE);   
          }  
        temp = cscore - ((((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_INSERT2MATCH])+mat->query->seg[i]->transition[GW_MATCH_BALANCE_5SS])) -  (0);    
        if( temp == (*h.access_main)(mat,i - 0,j - 9,INSERT) )   {  
          *reti = i - 0; 
          *retj = j - 9; 
          *retstate = INSERT;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-9,INSERT);  
            }  
          return (*h.access_main)(mat,i - 0,j - 9,INSERT);   
          }  
        temp = cscore - ((((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_MATCH2MATCH])+mat->query->seg[i]->transition[GW_MATCH_BALANCE_5SS])) -  (0); 
        if( temp == (*h.access_main)(mat,i - 0,j - 9,MATCH) )    {  
          *reti = i - 0; 
          *retj = j - 9; 
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-9,MATCH);   
            }  
          return (*h.access_main)(mat,i - 0,j - 9,MATCH);    
          }  
        warn("Major problem (!) - in GeneWise21 read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      case PY_1M :   
        temp = cscore - (mat->gp21->transition[GP21_PY2PY]) -  (mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)]);    
        if( temp == (*h.access_main)(mat,i - 0,j - 1,PY_1M) )    {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = PY_1M; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,PY_1M);   
            }  
          return (*h.access_main)(mat,i - 0,j - 1,PY_1M);    
          }  
        temp = cscore - (mat->gp21->transition[GP21_CENTRAL2PY]) -  (mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)]);   
        if( temp == (*h.access_main)(mat,i - 0,j - 1,CENTRAL_1M) )   {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = CENTRAL_1M;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,CENTRAL_1M);  
            }  
          return (*h.access_main)(mat,i - 0,j - 1,CENTRAL_1M);   
          }  
        warn("Major problem (!) - in GeneWise21 read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      case SPACER_1M :   
        temp = cscore - (mat->gp21->transition[GP21_SPACER2SPACER]) -  (mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)]);    
        if( temp == (*h.access_main)(mat,i - 0,j - 1,SPACER_1M) )    {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = SPACER_1M; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,SPACER_1M);   
            }  
          return (*h.access_main)(mat,i - 0,j - 1,SPACER_1M);    
          }  
        temp = cscore - (mat->gp21->transition[GP21_PY2SPACER]) -  (mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)]);    
        if( temp == (*h.access_main)(mat,i - 0,j - 1,PY_1M) )    {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = PY_1M; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,PY_1M);   
            }  
          return (*h.access_main)(mat,i - 0,j - 1,PY_1M);    
          }  
        warn("Major problem (!) - in GeneWise21 read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      case CENTRAL_2M :  
        temp = cscore - ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CENTRAL2CENTRAL])) -  (0); 
        if( temp == (*h.access_main)(mat,i - 0,j - 1,CENTRAL_2M) )   {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = CENTRAL_2M;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,CENTRAL_2M);  
            }  
          return (*h.access_main)(mat,i - 0,j - 1,CENTRAL_2M);   
          }  
        temp = cscore - (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_DELETE2MATCH])) -  (0);   
        if( temp == (*h.access_main)(mat,i - 0,j - 10,DELETE) )  {  
          *reti = i - 0; 
          *retj = j - 10;    
          *retstate = DELETE;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-10,DELETE); 
            }  
          return (*h.access_main)(mat,i - 0,j - 10,DELETE);  
          }  
        temp = cscore - ((((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_INSERT2MATCH])+mat->query->seg[i]->transition[GW_MATCH_BALANCE_5SS])) -  (0);    
        if( temp == (*h.access_main)(mat,i - 0,j - 10,INSERT) )  {  
          *reti = i - 0; 
          *retj = j - 10;    
          *retstate = INSERT;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-10,INSERT); 
            }  
          return (*h.access_main)(mat,i - 0,j - 10,INSERT);  
          }  
        temp = cscore - ((((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_MATCH2MATCH])+mat->query->seg[i]->transition[GW_MATCH_BALANCE_5SS])) -  (0); 
        if( temp == (*h.access_main)(mat,i - 0,j - 10,MATCH) )   {  
          *reti = i - 0; 
          *retj = j - 10;    
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-10,MATCH);  
            }  
          return (*h.access_main)(mat,i - 0,j - 10,MATCH);   
          }  
        warn("Major problem (!) - in GeneWise21 read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      case PY_2M :   
        temp = cscore - (mat->gp21->transition[GP21_PY2PY]) -  (mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)]);    
        if( temp == (*h.access_main)(mat,i - 0,j - 1,PY_2M) )    {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = PY_2M; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,PY_2M);   
            }  
          return (*h.access_main)(mat,i - 0,j - 1,PY_2M);    
          }  
        temp = cscore - (mat->gp21->transition[GP21_CENTRAL2PY]) -  (mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)]);   
        if( temp == (*h.access_main)(mat,i - 0,j - 1,CENTRAL_2M) )   {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = CENTRAL_2M;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,CENTRAL_2M);  
            }  
          return (*h.access_main)(mat,i - 0,j - 1,CENTRAL_2M);   
          }  
        warn("Major problem (!) - in GeneWise21 read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      case SPACER_2M :   
        temp = cscore - (mat->gp21->transition[GP21_SPACER2SPACER]) -  (mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)]);    
        if( temp == (*h.access_main)(mat,i - 0,j - 1,SPACER_2M) )    {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = SPACER_2M; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,SPACER_2M);   
            }  
          return (*h.access_main)(mat,i - 0,j - 1,SPACER_2M);    
          }  
        temp = cscore - (mat->gp21->transition[GP21_PY2SPACER]) -  (mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)]);    
        if( temp == (*h.access_main)(mat,i - 0,j - 1,PY_2M) )    {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = PY_2M; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,PY_2M);   
            }  
          return (*h.access_main)(mat,i - 0,j - 1,PY_2M);    
          }  
        warn("Major problem (!) - in GeneWise21 read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      case CENTRAL_0I :  
        temp = cscore - ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CENTRAL2CENTRAL])) -  (0); 
        if( temp == (*h.access_main)(mat,i - 0,j - 1,CENTRAL_0I) )   {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = CENTRAL_0I;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,CENTRAL_0I);  
            }  
          return (*h.access_main)(mat,i - 0,j - 1,CENTRAL_0I);   
          }  
        temp = cscore - (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_DELETE2INSERT])) -  (0);  
        if( temp == (*h.access_main)(mat,i - 0,j - 8,DELETE) )   {  
          *reti = i - 0; 
          *retj = j - 8; 
          *retstate = DELETE;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-8,DELETE);  
            }  
          return (*h.access_main)(mat,i - 0,j - 8,DELETE);   
          }  
        temp = cscore - ((((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_INSERT2INSERT])+mat->query->seg[i]->transition[GW_INSERT_BALANCE_5SS])) -  (0);  
        if( temp == (*h.access_main)(mat,i - 0,j - 8,INSERT) )   {  
          *reti = i - 0; 
          *retj = j - 8; 
          *retstate = INSERT;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-8,INSERT);  
            }  
          return (*h.access_main)(mat,i - 0,j - 8,INSERT);   
          }  
        temp = cscore - ((((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_MATCH2INSERT])+mat->query->seg[i]->transition[GW_MATCH_BALANCE_5SS])) -  (0);    
        if( temp == (*h.access_main)(mat,i - 0,j - 8,MATCH) )    {  
          *reti = i - 0; 
          *retj = j - 8; 
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-8,MATCH);   
            }  
          return (*h.access_main)(mat,i - 0,j - 8,MATCH);    
          }  
        warn("Major problem (!) - in GeneWise21 read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      case PY_0I :   
        temp = cscore - (mat->gp21->transition[GP21_PY2PY]) -  (mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)]);    
        if( temp == (*h.access_main)(mat,i - 0,j - 1,PY_0I) )    {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = PY_0I; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,PY_0I);   
            }  
          return (*h.access_main)(mat,i - 0,j - 1,PY_0I);    
          }  
        temp = cscore - (mat->gp21->transition[GP21_CENTRAL2PY]) -  (mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)]);   
        if( temp == (*h.access_main)(mat,i - 0,j - 1,CENTRAL_0I) )   {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = CENTRAL_0I;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,CENTRAL_0I);  
            }  
          return (*h.access_main)(mat,i - 0,j - 1,CENTRAL_0I);   
          }  
        warn("Major problem (!) - in GeneWise21 read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      case SPACER_0I :   
        temp = cscore - (mat->gp21->transition[GP21_SPACER2SPACER]) -  (mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)]);    
        if( temp == (*h.access_main)(mat,i - 0,j - 1,SPACER_0I) )    {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = SPACER_0I; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,SPACER_0I);   
            }  
          return (*h.access_main)(mat,i - 0,j - 1,SPACER_0I);    
          }  
        temp = cscore - (mat->gp21->transition[GP21_PY2SPACER]) -  (mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)]);    
        if( temp == (*h.access_main)(mat,i - 0,j - 1,PY_0I) )    {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = PY_0I; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,PY_0I);   
            }  
          return (*h.access_main)(mat,i - 0,j - 1,PY_0I);    
          }  
        warn("Major problem (!) - in GeneWise21 read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      case CENTRAL_1I :  
        temp = cscore - ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CENTRAL2CENTRAL])) -  (0); 
        if( temp == (*h.access_main)(mat,i - 0,j - 1,CENTRAL_1I) )   {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = CENTRAL_1I;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,CENTRAL_1I);  
            }  
          return (*h.access_main)(mat,i - 0,j - 1,CENTRAL_1I);   
          }  
        temp = cscore - (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_DELETE2INSERT])) -  (0);  
        if( temp == (*h.access_main)(mat,i - 0,j - 9,DELETE) )   {  
          *reti = i - 0; 
          *retj = j - 9; 
          *retstate = DELETE;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-9,DELETE);  
            }  
          return (*h.access_main)(mat,i - 0,j - 9,DELETE);   
          }  
        temp = cscore - (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_INSERT2INSERT])) -  (0);  
        if( temp == (*h.access_main)(mat,i - 0,j - 9,INSERT) )   {  
          *reti = i - 0; 
          *retj = j - 9; 
          *retstate = INSERT;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-9,INSERT);  
            }  
          return (*h.access_main)(mat,i - 0,j - 9,INSERT);   
          }  
        temp = cscore - (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_MATCH2INSERT])) -  (0);   
        if( temp == (*h.access_main)(mat,i - 0,j - 9,MATCH) )    {  
          *reti = i - 0; 
          *retj = j - 9; 
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-9,MATCH);   
            }  
          return (*h.access_main)(mat,i - 0,j - 9,MATCH);    
          }  
        warn("Major problem (!) - in GeneWise21 read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      case PY_1I :   
        temp = cscore - (mat->gp21->transition[GP21_PY2PY]) -  (mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)]);    
        if( temp == (*h.access_main)(mat,i - 0,j - 1,PY_1I) )    {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = PY_1I; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,PY_1I);   
            }  
          return (*h.access_main)(mat,i - 0,j - 1,PY_1I);    
          }  
        temp = cscore - (mat->gp21->transition[GP21_CENTRAL2PY]) -  (mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)]);   
        if( temp == (*h.access_main)(mat,i - 0,j - 1,CENTRAL_1I) )   {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = CENTRAL_1I;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,CENTRAL_1I);  
            }  
          return (*h.access_main)(mat,i - 0,j - 1,CENTRAL_1I);   
          }  
        warn("Major problem (!) - in GeneWise21 read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      case SPACER_1I :   
        temp = cscore - (mat->gp21->transition[GP21_SPACER2SPACER]) -  (mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)]);    
        if( temp == (*h.access_main)(mat,i - 0,j - 1,SPACER_1I) )    {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = SPACER_1I; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,SPACER_1I);   
            }  
          return (*h.access_main)(mat,i - 0,j - 1,SPACER_1I);    
          }  
        temp = cscore - (mat->gp21->transition[GP21_PY2SPACER]) -  (mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)]);    
        if( temp == (*h.access_main)(mat,i - 0,j - 1,PY_1I) )    {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = PY_1I; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,PY_1I);   
            }  
          return (*h.access_main)(mat,i - 0,j - 1,PY_1I);    
          }  
        warn("Major problem (!) - in GeneWise21 read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      case CENTRAL_2I :  
        temp = cscore - ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CENTRAL2CENTRAL])) -  (0); 
        if( temp == (*h.access_main)(mat,i - 0,j - 1,CENTRAL_2I) )   {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = CENTRAL_2I;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,CENTRAL_2I);  
            }  
          return (*h.access_main)(mat,i - 0,j - 1,CENTRAL_2I);   
          }  
        temp = cscore - (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_DELETE2INSERT])) -  (0);  
        if( temp == (*h.access_main)(mat,i - 0,j - 10,INSERT) )  {  
          *reti = i - 0; 
          *retj = j - 10;    
          *retstate = INSERT;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-10,INSERT); 
            }  
          return (*h.access_main)(mat,i - 0,j - 10,INSERT);  
          }  
        temp = cscore - (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_INSERT2INSERT])) -  (0);  
        if( temp == (*h.access_main)(mat,i - 0,j - 10,INSERT) )  {  
          *reti = i - 0; 
          *retj = j - 10;    
          *retstate = INSERT;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-10,INSERT); 
            }  
          return (*h.access_main)(mat,i - 0,j - 10,INSERT);  
          }  
        temp = cscore - (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_MATCH2INSERT])) -  (0);   
        if( temp == (*h.access_main)(mat,i - 0,j - 10,MATCH) )   {  
          *reti = i - 0; 
          *retj = j - 10;    
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-10,MATCH);  
            }  
          return (*h.access_main)(mat,i - 0,j - 10,MATCH);   
          }  
        warn("Major problem (!) - in GeneWise21 read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      case PY_2I :   
        temp = cscore - (mat->gp21->transition[GP21_PY2PY]) -  (mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)]);    
        if( temp == (*h.access_main)(mat,i - 0,j - 1,PY_2I) )    {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = PY_2I; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,PY_2I);   
            }  
          return (*h.access_main)(mat,i - 0,j - 1,PY_2I);    
          }  
        temp = cscore - (mat->gp21->transition[GP21_CENTRAL2PY]) -  (mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)]);   
        if( temp == (*h.access_main)(mat,i - 0,j - 1,CENTRAL_2I) )   {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = CENTRAL_2I;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,CENTRAL_2I);  
            }  
          return (*h.access_main)(mat,i - 0,j - 1,CENTRAL_2I);   
          }  
        warn("Major problem (!) - in GeneWise21 read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      case SPACER_2I :   
        temp = cscore - (mat->gp21->transition[GP21_SPACER2SPACER]) -  (mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)]);    
        if( temp == (*h.access_main)(mat,i - 0,j - 1,SPACER_2I) )    {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = SPACER_2I; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,SPACER_2I);   
            }  
          return (*h.access_main)(mat,i - 0,j - 1,SPACER_2I);    
          }  
        temp = cscore - (mat->gp21->transition[GP21_PY2SPACER]) -  (mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)]);    
        if( temp == (*h.access_main)(mat,i - 0,j - 1,PY_2I) )    {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = PY_2I; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,PY_2I);   
            }  
          return (*h.access_main)(mat,i - 0,j - 1,PY_2I);    
          }  
        warn("Major problem (!) - in GeneWise21 read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      default:   
        warn("Major problem (!) - in GeneWise21 read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      } /* end of Switch state  */ 
}    


/* Function:  max_calc_special_GeneWise21(mat,i,j,state,isspecial,reti,retj,retstate,retspecial,cellscore,h)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [GeneWise21 *]
 * Arg:                 i [UNKN ] Undocumented argument [int]
 * Arg:                 j [UNKN ] Undocumented argument [int]
 * Arg:             state [UNKN ] Undocumented argument [int]
 * Arg:         isspecial [UNKN ] Undocumented argument [boolean]
 * Arg:              reti [UNKN ] Undocumented argument [int *]
 * Arg:              retj [UNKN ] Undocumented argument [int *]
 * Arg:          retstate [UNKN ] Undocumented argument [int *]
 * Arg:        retspecial [UNKN ] Undocumented argument [boolean *]
 * Arg:         cellscore [UNKN ] Undocumented argument [int *]
 * Arg:                 h [UNKN ] Undocumented argument [GeneWise21_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int max_calc_special_GeneWise21(GeneWise21 * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,GeneWise21_access_func_holder h) 
{
    register int temp;   
    register int cscore; 


    *reti = (*retj) = (*retstate) = GeneWise21_READ_OFF_ERROR;   


    if( j < 0 || j > mat->target->seq->len)  {  
      warn("In GeneWise21 matrix special read off - out of bounds on matrix [j is %d in special]",j);    
      return -1;     
      }  


    cscore = (*h.access_special)(mat,i,j,state); 
    switch(state)    { /*switch on special states*/ 
      case START :   
      case END :     
        /* source AFTER_RND is a special */ 
        temp = cscore - (0) - (0);   
        if( temp == (*h.access_special)(mat,i - 0,j - 1,AFTER_RND) ) {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = AFTER_RND; 
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_special)(mat,i-0,j-1,AFTER_RND);    
            }  
          return (*h.access_special)(mat,i - 0,j - 1,AFTER_RND) ;    
          }  
      case BEFORE_RND :  
        /* source BEFORE_CODON is a special */ 
        temp = cscore - ((mat->rndbase->base[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CDS2RND])) - (0);  
        if( temp == (*h.access_special)(mat,i - 0,j - 1,BEFORE_CODON) )  {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = BEFORE_CODON;  
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_special)(mat,i-0,j-1,BEFORE_CODON);     
            }  
          return (*h.access_special)(mat,i - 0,j - 1,BEFORE_CODON) ;     
          }  
        /* source BEFORE_RND is a special */ 
        temp = cscore - ((mat->rndbase->base[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_RND2RND])) - (0);  
        if( temp == (*h.access_special)(mat,i - 0,j - 1,BEFORE_RND) )    {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = BEFORE_RND;    
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_special)(mat,i-0,j-1,BEFORE_RND);   
            }  
          return (*h.access_special)(mat,i - 0,j - 1,BEFORE_RND) ;   
          }  
        /* source START is a special */ 
        temp = cscore - (mat->rndbase->base[CSEQ_GENOMIC_BASE(mat->target,j)]) - (0);    
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
      case BEFORE_CODON :    
        /* source BEFORE_SPACER is a special */ 
        temp = cscore - ((mat->gp21->transition[GP21_SPACER2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-1)))) - (0);    
        if( temp == (*h.access_special)(mat,i - 0,j - 3,BEFORE_SPACER) ) {  
          *reti = i - 0; 
          *retj = j - 3; 
          *retstate = BEFORE_SPACER; 
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_special)(mat,i-0,j-3,BEFORE_SPACER);    
            }  
          return (*h.access_special)(mat,i - 0,j - 3,BEFORE_SPACER) ;    
          }  
        /* source BEFORE_SPACER is a special */ 
        temp = cscore - ((mat->gp21->transition[GP21_SPACER2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-2)))) - (0);    
        if( temp == (*h.access_special)(mat,i - 0,j - 4,BEFORE_SPACER) ) {  
          *reti = i - 0; 
          *retj = j - 4; 
          *retstate = BEFORE_SPACER; 
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_special)(mat,i-0,j-4,BEFORE_SPACER);    
            }  
          return (*h.access_special)(mat,i - 0,j - 4,BEFORE_SPACER) ;    
          }  
        /* source BEFORE_SPACER is a special */ 
        temp = cscore - (((mat->gp21->transition[GP21_SPACER2CDS]+mat->rndcodon->codon[CSEQ_GENOMIC_CODON(mat->target,j)])+CSEQ_GENOMIC_3SS(mat->target,(j-3)))) - (0);  
        if( temp == (*h.access_special)(mat,i - 0,j - 5,BEFORE_SPACER) ) {  
          *reti = i - 0; 
          *retj = j - 5; 
          *retstate = BEFORE_SPACER; 
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_special)(mat,i-0,j-5,BEFORE_SPACER);    
            }  
          return (*h.access_special)(mat,i - 0,j - 5,BEFORE_SPACER) ;    
          }  
        /* source BEFORE_CODON is a special */ 
        temp = cscore - ((mat->rndcodon->codon[CSEQ_GENOMIC_CODON(mat->target,j)]+mat->gp21->transition[GP21_CDS2CDS])) - (0);   
        if( temp == (*h.access_special)(mat,i - 0,j - 3,BEFORE_CODON) )  {  
          *reti = i - 0; 
          *retj = j - 3; 
          *retstate = BEFORE_CODON;  
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_special)(mat,i-0,j-3,BEFORE_CODON);     
            }  
          return (*h.access_special)(mat,i - 0,j - 3,BEFORE_CODON) ;     
          }  
        /* source BEFORE_RND is a special */ 
        temp = cscore - ((mat->rndcodon->codon[CSEQ_GENOMIC_CODON(mat->target,j)]+mat->gp21->transition[GP21_RND2CDS])) - (0);   
        if( temp == (*h.access_special)(mat,i - 0,j - 3,BEFORE_RND) )    {  
          *reti = i - 0; 
          *retj = j - 3; 
          *retstate = BEFORE_RND;    
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_special)(mat,i-0,j-3,BEFORE_RND);   
            }  
          return (*h.access_special)(mat,i - 0,j - 3,BEFORE_RND) ;   
          }  
      case BEFORE_CENTRAL :  
        /* source BEFORE_CODON is a special */ 
        temp = cscore - ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))) - (0);  
        if( temp == (*h.access_special)(mat,i - 0,j - 10,BEFORE_CODON) ) {  
          *reti = i - 0; 
          *retj = j - 10;    
          *retstate = BEFORE_CODON;  
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_special)(mat,i-0,j-10,BEFORE_CODON);    
            }  
          return (*h.access_special)(mat,i - 0,j - 10,BEFORE_CODON) ;    
          }  
        /* source BEFORE_CODON is a special */ 
        temp = cscore - ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))) - (0);  
        if( temp == (*h.access_special)(mat,i - 0,j - 9,BEFORE_CODON) )  {  
          *reti = i - 0; 
          *retj = j - 9; 
          *retstate = BEFORE_CODON;  
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_special)(mat,i-0,j-9,BEFORE_CODON);     
            }  
          return (*h.access_special)(mat,i - 0,j - 9,BEFORE_CODON) ;     
          }  
        /* source BEFORE_CODON is a special */ 
        temp = cscore - ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))) - (0);  
        if( temp == (*h.access_special)(mat,i - 0,j - 8,BEFORE_CODON) )  {  
          *reti = i - 0; 
          *retj = j - 8; 
          *retstate = BEFORE_CODON;  
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_special)(mat,i-0,j-8,BEFORE_CODON);     
            }  
          return (*h.access_special)(mat,i - 0,j - 8,BEFORE_CODON) ;     
          }  
        /* source BEFORE_CENTRAL is a special */ 
        temp = cscore - ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CENTRAL2CENTRAL])) - (0);  
        if( temp == (*h.access_special)(mat,i - 0,j - 1,BEFORE_CENTRAL) )    {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = BEFORE_CENTRAL;    
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_special)(mat,i-0,j-1,BEFORE_CENTRAL);   
            }  
          return (*h.access_special)(mat,i - 0,j - 1,BEFORE_CENTRAL) ;   
          }  
      case BEFORE_PY_TRACT :     
        /* source BEFORE_PY_TRACT is a special */ 
        temp = cscore - (mat->gp21->transition[GP21_PY2PY]) - (mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)]);     
        if( temp == (*h.access_special)(mat,i - 0,j - 1,BEFORE_PY_TRACT) )   {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = BEFORE_PY_TRACT;   
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_special)(mat,i-0,j-1,BEFORE_PY_TRACT);  
            }  
          return (*h.access_special)(mat,i - 0,j - 1,BEFORE_PY_TRACT) ;  
          }  
        /* source BEFORE_CENTRAL is a special */ 
        temp = cscore - (mat->gp21->transition[GP21_CENTRAL2PY]) - (mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)]);    
        if( temp == (*h.access_special)(mat,i - 0,j - 1,BEFORE_CENTRAL) )    {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = BEFORE_CENTRAL;    
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_special)(mat,i-0,j-1,BEFORE_CENTRAL);   
            }  
          return (*h.access_special)(mat,i - 0,j - 1,BEFORE_CENTRAL) ;   
          }  
      case BEFORE_SPACER :   
        /* source BEFORE_SPACER is a special */ 
        temp = cscore - (mat->gp21->transition[GP21_SPACER2SPACER]) - (mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)]);     
        if( temp == (*h.access_special)(mat,i - 0,j - 1,BEFORE_SPACER) ) {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = BEFORE_SPACER; 
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_special)(mat,i-0,j-1,BEFORE_SPACER);    
            }  
          return (*h.access_special)(mat,i - 0,j - 1,BEFORE_SPACER) ;    
          }  
        /* source BEFORE_PY_TRACT is a special */ 
        temp = cscore - (mat->gp21->transition[GP21_PY2SPACER]) - (mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)]);     
        if( temp == (*h.access_special)(mat,i - 0,j - 1,BEFORE_PY_TRACT) )   {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = BEFORE_PY_TRACT;   
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_special)(mat,i-0,j-1,BEFORE_PY_TRACT);  
            }  
          return (*h.access_special)(mat,i - 0,j - 1,BEFORE_PY_TRACT) ;  
          }  
      case AFTER_RND :   
        /* source AFTER_CODON is a special */ 
        temp = cscore - ((mat->rndbase->base[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CDS2RND])) - (0);  
        if( temp == (*h.access_special)(mat,i - 0,j - 1,AFTER_CODON) )   {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = AFTER_CODON;   
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_special)(mat,i-0,j-1,AFTER_CODON);  
            }  
          return (*h.access_special)(mat,i - 0,j - 1,AFTER_CODON) ;  
          }  
        /* source AFTER_RND is a special */ 
        temp = cscore - ((mat->rndbase->base[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_RND2RND])) - (0);  
        if( temp == (*h.access_special)(mat,i - 0,j - 1,AFTER_RND) ) {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = AFTER_RND; 
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_special)(mat,i-0,j-1,AFTER_RND);    
            }  
          return (*h.access_special)(mat,i - 0,j - 1,AFTER_RND) ;    
          }  
      case AFTER_CODON :     
        /* source AFTER_SPACER is a special */ 
        temp = cscore - ((mat->gp21->transition[GP21_SPACER2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-1)))) - (0);    
        if( temp == (*h.access_special)(mat,i - 0,j - 3,AFTER_SPACER) )  {  
          *reti = i - 0; 
          *retj = j - 3; 
          *retstate = AFTER_SPACER;  
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_special)(mat,i-0,j-3,AFTER_SPACER);     
            }  
          return (*h.access_special)(mat,i - 0,j - 3,AFTER_SPACER) ;     
          }  
        /* source AFTER_SPACER is a special */ 
        temp = cscore - ((mat->gp21->transition[GP21_SPACER2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-2)))) - (0);    
        if( temp == (*h.access_special)(mat,i - 0,j - 4,AFTER_SPACER) )  {  
          *reti = i - 0; 
          *retj = j - 4; 
          *retstate = AFTER_SPACER;  
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_special)(mat,i-0,j-4,AFTER_SPACER);     
            }  
          return (*h.access_special)(mat,i - 0,j - 4,AFTER_SPACER) ;     
          }  
        /* source AFTER_SPACER is a special */ 
        temp = cscore - (((mat->gp21->transition[GP21_SPACER2CDS]+mat->rndcodon->codon[CSEQ_GENOMIC_CODON(mat->target,j)])+CSEQ_GENOMIC_3SS(mat->target,(j-3)))) - (0);  
        if( temp == (*h.access_special)(mat,i - 0,j - 5,AFTER_SPACER) )  {  
          *reti = i - 0; 
          *retj = j - 5; 
          *retstate = AFTER_SPACER;  
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_special)(mat,i-0,j-5,AFTER_SPACER);     
            }  
          return (*h.access_special)(mat,i - 0,j - 5,AFTER_SPACER) ;     
          }  
        /* source AFTER_CODON is a special */ 
        temp = cscore - ((mat->rndcodon->codon[CSEQ_GENOMIC_CODON(mat->target,j)]+mat->gp21->transition[GP21_CDS2CDS])) - (0);   
        if( temp == (*h.access_special)(mat,i - 0,j - 3,AFTER_CODON) )   {  
          *reti = i - 0; 
          *retj = j - 3; 
          *retstate = AFTER_CODON;   
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_special)(mat,i-0,j-3,AFTER_CODON);  
            }  
          return (*h.access_special)(mat,i - 0,j - 3,AFTER_CODON) ;  
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
        /* source AFTER_RND is a special */ 
        temp = cscore - ((mat->rndcodon->codon[CSEQ_GENOMIC_CODON(mat->target,j)]+mat->gp21->transition[GP21_RND2CDS])) - (0);   
        if( temp == (*h.access_special)(mat,i - 0,j - 3,AFTER_RND) ) {  
          *reti = i - 0; 
          *retj = j - 3; 
          *retstate = AFTER_RND; 
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_special)(mat,i-0,j-3,AFTER_RND);    
            }  
          return (*h.access_special)(mat,i - 0,j - 3,AFTER_RND) ;    
          }  
      case AFTER_CENTRAL :   
        /* source AFTER_CODON is a special */ 
        temp = cscore - ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))) - (0);  
        if( temp == (*h.access_special)(mat,i - 0,j - 10,AFTER_CODON) )  {  
          *reti = i - 0; 
          *retj = j - 10;    
          *retstate = AFTER_CODON;   
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_special)(mat,i-0,j-10,AFTER_CODON);     
            }  
          return (*h.access_special)(mat,i - 0,j - 10,AFTER_CODON) ;     
          }  
        /* source AFTER_CODON is a special */ 
        temp = cscore - ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))) - (0);  
        if( temp == (*h.access_special)(mat,i - 0,j - 9,AFTER_CODON) )   {  
          *reti = i - 0; 
          *retj = j - 9; 
          *retstate = AFTER_CODON;   
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_special)(mat,i-0,j-9,AFTER_CODON);  
            }  
          return (*h.access_special)(mat,i - 0,j - 9,AFTER_CODON) ;  
          }  
        /* source AFTER_CODON is a special */ 
        temp = cscore - ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))) - (0);  
        if( temp == (*h.access_special)(mat,i - 0,j - 8,AFTER_CODON) )   {  
          *reti = i - 0; 
          *retj = j - 8; 
          *retstate = AFTER_CODON;   
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_special)(mat,i-0,j-8,AFTER_CODON);  
            }  
          return (*h.access_special)(mat,i - 0,j - 8,AFTER_CODON) ;  
          }  
        /* source AFTER_CENTRAL is a special */ 
        temp = cscore - ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CENTRAL2CENTRAL])) - (0);  
        if( temp == (*h.access_special)(mat,i - 0,j - 1,AFTER_CENTRAL) ) {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = AFTER_CENTRAL; 
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_special)(mat,i-0,j-1,AFTER_CENTRAL);    
            }  
          return (*h.access_special)(mat,i - 0,j - 1,AFTER_CENTRAL) ;    
          }  
      case AFTER_PY_TRACT :  
        /* source AFTER_PY_TRACT is a special */ 
        temp = cscore - (mat->gp21->transition[GP21_PY2PY]) - (mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)]);     
        if( temp == (*h.access_special)(mat,i - 0,j - 1,AFTER_PY_TRACT) )    {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = AFTER_PY_TRACT;    
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_special)(mat,i-0,j-1,AFTER_PY_TRACT);   
            }  
          return (*h.access_special)(mat,i - 0,j - 1,AFTER_PY_TRACT) ;   
          }  
        /* source AFTER_CENTRAL is a special */ 
        temp = cscore - (mat->gp21->transition[GP21_CENTRAL2PY]) - (mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)]);    
        if( temp == (*h.access_special)(mat,i - 0,j - 1,AFTER_CENTRAL) ) {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = AFTER_CENTRAL; 
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_special)(mat,i-0,j-1,AFTER_CENTRAL);    
            }  
          return (*h.access_special)(mat,i - 0,j - 1,AFTER_CENTRAL) ;    
          }  
      case AFTER_SPACER :    
        /* source AFTER_SPACER is a special */ 
        temp = cscore - (mat->gp21->transition[GP21_SPACER2SPACER]) - (mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)]);     
        if( temp == (*h.access_special)(mat,i - 0,j - 1,AFTER_SPACER) )  {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = AFTER_SPACER;  
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_special)(mat,i-0,j-1,AFTER_SPACER);     
            }  
          return (*h.access_special)(mat,i - 0,j - 1,AFTER_SPACER) ;     
          }  
        /* source AFTER_PY_TRACT is a special */ 
        temp = cscore - (mat->gp21->transition[GP21_PY2SPACER]) - (mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)]);     
        if( temp == (*h.access_special)(mat,i - 0,j - 1,AFTER_PY_TRACT) )    {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = AFTER_PY_TRACT;    
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_special)(mat,i-0,j-1,AFTER_PY_TRACT);   
            }  
          return (*h.access_special)(mat,i - 0,j - 1,AFTER_PY_TRACT) ;   
          }  
      default:   
        warn("Major problem (!) - in GeneWise21 read off, position %d,%d state %d no source found  dropped into default on source switch!",i,j,state);   
        return (-1); 
      } /* end of switch on special states */ 
}    


/* Function:  calculate_GeneWise21(mat)
 *
 * Descrip:    This function calculates the GeneWise21 matrix when in explicit mode
 *             To allocate the matrix use /allocate_Expl_GeneWise21
 *
 *
 * Arg:        mat [UNKN ] GeneWise21 which contains explicit basematrix memory [GeneWise21 *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean calculate_GeneWise21(GeneWise21 * mat) 
{
    int i;   
    int j;   
    int leni;    
    int lenj;    
    int tot; 
    int num; 


    if( mat->basematrix->type != BASEMATRIX_TYPE_EXPLICIT )  {  
      warn("in calculate_GeneWise21, passed a non Explicit matrix type, cannot calculate!"); 
      return FALSE;  
      }  


    leni = mat->leni;    
    lenj = mat->lenj;    
    tot = leni * lenj;   
    num = 0; 


    start_reporting("GeneWise21 Matrix calculation: ");  
    for(j=0;j<lenj;j++)  {  
      auto int score;    
      auto int temp;     
      for(i=0;i<leni;i++)    {  
        if( num%1000 == 0 )  
          log_full_error(REPORT,0,"[%7d] Cells %2d%%%%",num,num*100/tot);    
        num++;   


        /* For state MATCH */ 
        /* setting first movement to score */ 
        score = GeneWise21_EXPL_MATRIX(mat,i-1,j-3,MATCH) + (mat->query->seg[i]->transition[GW_MATCH2MATCH]+mat->query->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]);   
        /* From state INSERT to state MATCH */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-1,j-3,INSERT) + (mat->query->seg[i]->transition[GW_INSERT2MATCH]+mat->query->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state MATCH */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-1,j-3,DELETE) + (mat->query->seg[i]->transition[GW_DELETE2MATCH]+mat->query->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state BEFORE_CODON to state MATCH */ 
        temp = GeneWise21_EXPL_SPECIAL(mat,i-1,j-3,BEFORE_CODON) + (mat->query->seg[i]->transition[GW_START2MATCH]+mat->query->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]);    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SPACER_0M to state MATCH */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-1,j-6,SPACER_0M) + (((mat->gp21->transition[GP21_SPACER2CDS]+mat->query->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)])+CSEQ_GENOMIC_3SS(mat->target,(j-3)))+mat->query->seg[i]->transition[GW_MATCH_BALANCE_3SS]);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state PY_0M to state MATCH */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-1,j-6,PY_0M) + (((mat->gp21->transition[GP21_PY2CDS]+mat->query->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)])+CSEQ_GENOMIC_3SS(mat->target,(j-3)))+mat->query->seg[i]->transition[GW_MATCH_BALANCE_3SS]);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SPACER_1M to state MATCH */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-1,j-5,SPACER_1M) + (mat->gp21->transition[GP21_SPACER2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-2)));     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state PY_1M to state MATCH */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-1,j-5,PY_1M) + (mat->gp21->transition[GP21_PY2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-2)));     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SPACER_2M to state MATCH */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-1,j-4,SPACER_2M) + (mat->gp21->transition[GP21_SPACER2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-1)));     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state PY_2M to state MATCH */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-1,j-4,PY_2M) + (mat->gp21->transition[GP21_PY2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-1)));     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-1,j-2,MATCH) + mat->gp21->transition[GP21_DELETE_1_BASE];    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-1,j-1,MATCH) + mat->gp21->transition[GP21_DELETE_2_BASE];    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-1,j-4,MATCH) + mat->gp21->transition[GP21_INSERT_1_BASE];    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-1,j-5,MATCH) + mat->gp21->transition[GP21_INSERT_2_BASE];    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for MATCH */ 
        /* Add any movement independant score and put away */ 
         score += CSEQ_GENOMIC_CDSPOT(mat->target,j);    
         GeneWise21_EXPL_MATRIX(mat,i,j,MATCH) = score;  


        /* state MATCH is a source for special AFTER_CODON */ 
        temp = score + (mat->query->seg[i]->transition[GW_MATCH2END]) + (0) ;    
        if( temp > GeneWise21_EXPL_SPECIAL(mat,i,j,AFTER_CODON) )    {  
          GeneWise21_EXPL_SPECIAL(mat,i,j,AFTER_CODON) = temp;   
          }  




        /* Finished calculating state MATCH */ 


        /* For state INSERT */ 
        /* setting first movement to score */ 
        score = GeneWise21_EXPL_MATRIX(mat,i-0,j-3,MATCH) + (mat->query->seg[i]->transition[GW_MATCH2INSERT]+mat->query->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)]);     
        /* From state INSERT to state INSERT */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-0,j-3,INSERT) + (mat->query->seg[i]->transition[GW_INSERT2INSERT]+mat->query->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)]);    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state INSERT */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-0,j-3,DELETE) + (mat->query->seg[i]->transition[GW_DELETE2INSERT]+mat->query->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)]);    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state BEFORE_CODON to state INSERT */ 
        temp = GeneWise21_EXPL_SPECIAL(mat,i-0,j-3,BEFORE_CODON) + (mat->query->seg[i]->transition[GW_START2INSERT]+mat->query->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)]);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SPACER_0I to state INSERT */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-1,j-6,SPACER_0I) + (((mat->gp21->transition[GP21_SPACER2CDS]+mat->query->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)])+CSEQ_GENOMIC_3SS(mat->target,(j-3)))+mat->query->seg[i]->transition[GW_MATCH_BALANCE_3SS]);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SPACER_1I to state INSERT */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-1,j-5,SPACER_1I) + (mat->gp21->transition[GP21_SPACER2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-2)));     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SPACER_2I to state INSERT */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-1,j-4,SPACER_2I) + (mat->gp21->transition[GP21_SPACER2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-1)));     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state PY_0I to state INSERT */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-1,j-6,PY_0I) + (((mat->gp21->transition[GP21_PY2CDS]+mat->query->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)])+CSEQ_GENOMIC_3SS(mat->target,(j-3)))+mat->query->seg[i]->transition[GW_MATCH_BALANCE_3SS]);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state PY_1I to state INSERT */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-1,j-5,PY_1I) + (mat->gp21->transition[GP21_PY2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-2)));     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state PY_2I to state INSERT */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-1,j-4,PY_2I) + (mat->gp21->transition[GP21_PY2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-1)));     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state INSERT */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-1,j-2,INSERT) + mat->gp21->transition[GP21_DELETE_1_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state INSERT */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-1,j-1,INSERT) + mat->gp21->transition[GP21_DELETE_2_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state INSERT */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-1,j-4,INSERT) + mat->gp21->transition[GP21_INSERT_1_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state INSERT */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-1,j-5,INSERT) + mat->gp21->transition[GP21_INSERT_2_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for INSERT */ 
        /* Add any movement independant score and put away */ 
         score += CSEQ_GENOMIC_CDSPOT(mat->target,j);    
         GeneWise21_EXPL_MATRIX(mat,i,j,INSERT) = score; 


        /* state INSERT is a source for special AFTER_CODON */ 
        temp = score + (mat->query->seg[i]->transition[GW_INSERT2END]) + (0) ;   
        if( temp > GeneWise21_EXPL_SPECIAL(mat,i,j,AFTER_CODON) )    {  
          GeneWise21_EXPL_SPECIAL(mat,i,j,AFTER_CODON) = temp;   
          }  




        /* Finished calculating state INSERT */ 


        /* For state DELETE */ 
        /* setting first movement to score */ 
        score = GeneWise21_EXPL_MATRIX(mat,i-1,j-0,MATCH) + mat->query->seg[i]->transition[GW_MATCH2DELETE];     
        /* From state INSERT to state DELETE */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-1,j-0,INSERT) + mat->query->seg[i]->transition[GW_INSERT2DELETE];    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state DELETE */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-1,j-0,DELETE) + mat->query->seg[i]->transition[GW_DELETE2DELETE];    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state BEFORE_CODON to state DELETE */ 
        temp = GeneWise21_EXPL_SPECIAL(mat,i-1,j-0,BEFORE_CODON) + mat->query->seg[i]->transition[GW_START2DELETE];  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for DELETE */ 
        /* Add any movement independant score and put away */ 
         GeneWise21_EXPL_MATRIX(mat,i,j,DELETE) = score; 


        /* state DELETE is a source for special AFTER_CODON */ 
        temp = score + (mat->query->seg[i]->transition[GW_DELETE2END]) + (0) ;   
        if( temp > GeneWise21_EXPL_SPECIAL(mat,i,j,AFTER_CODON) )    {  
          GeneWise21_EXPL_SPECIAL(mat,i,j,AFTER_CODON) = temp;   
          }  




        /* Finished calculating state DELETE */ 


        /* For state CENTRAL_0M */ 
        /* setting first movement to score */ 
        score = GeneWise21_EXPL_MATRIX(mat,i-0,j-8,MATCH) + (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_MATCH2MATCH])+mat->query->seg[i]->transition[GW_MATCH_BALANCE_5SS]);  
        /* From state INSERT to state CENTRAL_0M */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-0,j-8,INSERT) + (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_INSERT2MATCH])+mat->query->seg[i]->transition[GW_INSERT_BALANCE_5SS]);    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state CENTRAL_0M */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-0,j-8,DELETE) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_DELETE2MATCH]);    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state CENTRAL_0M to state CENTRAL_0M */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-0,j-1,CENTRAL_0M) + (mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CENTRAL2CENTRAL]);  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for CENTRAL_0M */ 
        /* Add any movement independant score and put away */ 
         GeneWise21_EXPL_MATRIX(mat,i,j,CENTRAL_0M) = score; 


        /* Finished calculating state CENTRAL_0M */ 


        /* For state PY_0M */ 
        /* setting first movement to score */ 
        score = GeneWise21_EXPL_MATRIX(mat,i-0,j-1,CENTRAL_0M) + mat->gp21->transition[GP21_CENTRAL2PY];     
        /* From state PY_0M to state PY_0M */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-0,j-1,PY_0M) + mat->gp21->transition[GP21_PY2PY];    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for PY_0M */ 
        /* Add any movement independant score and put away */ 
         score += mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)];   
         GeneWise21_EXPL_MATRIX(mat,i,j,PY_0M) = score;  


        /* Finished calculating state PY_0M */ 


        /* For state SPACER_0M */ 
        /* setting first movement to score */ 
        score = GeneWise21_EXPL_MATRIX(mat,i-0,j-1,PY_0M) + mat->gp21->transition[GP21_PY2SPACER];   
        /* From state SPACER_0M to state SPACER_0M */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-0,j-1,SPACER_0M) + mat->gp21->transition[GP21_SPACER2SPACER];    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for SPACER_0M */ 
        /* Add any movement independant score and put away */ 
         score += mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)];   
         GeneWise21_EXPL_MATRIX(mat,i,j,SPACER_0M) = score;  


        /* Finished calculating state SPACER_0M */ 


        /* For state CENTRAL_1M */ 
        /* setting first movement to score */ 
        score = GeneWise21_EXPL_MATRIX(mat,i-0,j-9,MATCH) + (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_MATCH2MATCH])+mat->query->seg[i]->transition[GW_MATCH_BALANCE_5SS]);  
        /* From state INSERT to state CENTRAL_1M */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-0,j-9,INSERT) + (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_INSERT2MATCH])+mat->query->seg[i]->transition[GW_MATCH_BALANCE_5SS]);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state CENTRAL_1M */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-0,j-9,DELETE) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_DELETE2MATCH]);    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state CENTRAL_1M to state CENTRAL_1M */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-0,j-1,CENTRAL_1M) + (mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CENTRAL2CENTRAL]);  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for CENTRAL_1M */ 
        /* Add any movement independant score and put away */ 
         GeneWise21_EXPL_MATRIX(mat,i,j,CENTRAL_1M) = score; 


        /* Finished calculating state CENTRAL_1M */ 


        /* For state PY_1M */ 
        /* setting first movement to score */ 
        score = GeneWise21_EXPL_MATRIX(mat,i-0,j-1,CENTRAL_1M) + mat->gp21->transition[GP21_CENTRAL2PY];     
        /* From state PY_1M to state PY_1M */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-0,j-1,PY_1M) + mat->gp21->transition[GP21_PY2PY];    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for PY_1M */ 
        /* Add any movement independant score and put away */ 
         score += mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)];   
         GeneWise21_EXPL_MATRIX(mat,i,j,PY_1M) = score;  


        /* Finished calculating state PY_1M */ 


        /* For state SPACER_1M */ 
        /* setting first movement to score */ 
        score = GeneWise21_EXPL_MATRIX(mat,i-0,j-1,PY_1M) + mat->gp21->transition[GP21_PY2SPACER];   
        /* From state SPACER_1M to state SPACER_1M */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-0,j-1,SPACER_1M) + mat->gp21->transition[GP21_SPACER2SPACER];    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for SPACER_1M */ 
        /* Add any movement independant score and put away */ 
         score += mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)];   
         GeneWise21_EXPL_MATRIX(mat,i,j,SPACER_1M) = score;  


        /* Finished calculating state SPACER_1M */ 


        /* For state CENTRAL_2M */ 
        /* setting first movement to score */ 
        score = GeneWise21_EXPL_MATRIX(mat,i-0,j-10,MATCH) + (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_MATCH2MATCH])+mat->query->seg[i]->transition[GW_MATCH_BALANCE_5SS]);     
        /* From state INSERT to state CENTRAL_2M */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-0,j-10,INSERT) + (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_INSERT2MATCH])+mat->query->seg[i]->transition[GW_MATCH_BALANCE_5SS]);    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state CENTRAL_2M */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-0,j-10,DELETE) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_DELETE2MATCH]);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state CENTRAL_2M to state CENTRAL_2M */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-0,j-1,CENTRAL_2M) + (mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CENTRAL2CENTRAL]);  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for CENTRAL_2M */ 
        /* Add any movement independant score and put away */ 
         GeneWise21_EXPL_MATRIX(mat,i,j,CENTRAL_2M) = score; 


        /* Finished calculating state CENTRAL_2M */ 


        /* For state PY_2M */ 
        /* setting first movement to score */ 
        score = GeneWise21_EXPL_MATRIX(mat,i-0,j-1,CENTRAL_2M) + mat->gp21->transition[GP21_CENTRAL2PY];     
        /* From state PY_2M to state PY_2M */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-0,j-1,PY_2M) + mat->gp21->transition[GP21_PY2PY];    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for PY_2M */ 
        /* Add any movement independant score and put away */ 
         score += mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)];   
         GeneWise21_EXPL_MATRIX(mat,i,j,PY_2M) = score;  


        /* Finished calculating state PY_2M */ 


        /* For state SPACER_2M */ 
        /* setting first movement to score */ 
        score = GeneWise21_EXPL_MATRIX(mat,i-0,j-1,PY_2M) + mat->gp21->transition[GP21_PY2SPACER];   
        /* From state SPACER_2M to state SPACER_2M */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-0,j-1,SPACER_2M) + mat->gp21->transition[GP21_SPACER2SPACER];    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for SPACER_2M */ 
        /* Add any movement independant score and put away */ 
         score += mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)];   
         GeneWise21_EXPL_MATRIX(mat,i,j,SPACER_2M) = score;  


        /* Finished calculating state SPACER_2M */ 


        /* For state CENTRAL_0I */ 
        /* setting first movement to score */ 
        score = GeneWise21_EXPL_MATRIX(mat,i-0,j-8,MATCH) + (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_MATCH2INSERT])+mat->query->seg[i]->transition[GW_MATCH_BALANCE_5SS]);     
        /* From state INSERT to state CENTRAL_0I */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-0,j-8,INSERT) + (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_INSERT2INSERT])+mat->query->seg[i]->transition[GW_INSERT_BALANCE_5SS]);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state CENTRAL_0I */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-0,j-8,DELETE) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_DELETE2INSERT]);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state CENTRAL_0I to state CENTRAL_0I */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-0,j-1,CENTRAL_0I) + (mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CENTRAL2CENTRAL]);  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for CENTRAL_0I */ 
        /* Add any movement independant score and put away */ 
         GeneWise21_EXPL_MATRIX(mat,i,j,CENTRAL_0I) = score; 


        /* Finished calculating state CENTRAL_0I */ 


        /* For state PY_0I */ 
        /* setting first movement to score */ 
        score = GeneWise21_EXPL_MATRIX(mat,i-0,j-1,CENTRAL_0I) + mat->gp21->transition[GP21_CENTRAL2PY];     
        /* From state PY_0I to state PY_0I */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-0,j-1,PY_0I) + mat->gp21->transition[GP21_PY2PY];    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for PY_0I */ 
        /* Add any movement independant score and put away */ 
         score += mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)];   
         GeneWise21_EXPL_MATRIX(mat,i,j,PY_0I) = score;  


        /* Finished calculating state PY_0I */ 


        /* For state SPACER_0I */ 
        /* setting first movement to score */ 
        score = GeneWise21_EXPL_MATRIX(mat,i-0,j-1,PY_0I) + mat->gp21->transition[GP21_PY2SPACER];   
        /* From state SPACER_0I to state SPACER_0I */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-0,j-1,SPACER_0I) + mat->gp21->transition[GP21_SPACER2SPACER];    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for SPACER_0I */ 
        /* Add any movement independant score and put away */ 
         score += mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)];   
         GeneWise21_EXPL_MATRIX(mat,i,j,SPACER_0I) = score;  


        /* Finished calculating state SPACER_0I */ 


        /* For state CENTRAL_1I */ 
        /* setting first movement to score */ 
        score = GeneWise21_EXPL_MATRIX(mat,i-0,j-9,MATCH) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_MATCH2INSERT]);    
        /* From state INSERT to state CENTRAL_1I */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-0,j-9,INSERT) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_INSERT2INSERT]);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state CENTRAL_1I */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-0,j-9,DELETE) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_DELETE2INSERT]);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state CENTRAL_1I to state CENTRAL_1I */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-0,j-1,CENTRAL_1I) + (mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CENTRAL2CENTRAL]);  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for CENTRAL_1I */ 
        /* Add any movement independant score and put away */ 
         GeneWise21_EXPL_MATRIX(mat,i,j,CENTRAL_1I) = score; 


        /* Finished calculating state CENTRAL_1I */ 


        /* For state PY_1I */ 
        /* setting first movement to score */ 
        score = GeneWise21_EXPL_MATRIX(mat,i-0,j-1,CENTRAL_1I) + mat->gp21->transition[GP21_CENTRAL2PY];     
        /* From state PY_1I to state PY_1I */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-0,j-1,PY_1I) + mat->gp21->transition[GP21_PY2PY];    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for PY_1I */ 
        /* Add any movement independant score and put away */ 
         score += mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)];   
         GeneWise21_EXPL_MATRIX(mat,i,j,PY_1I) = score;  


        /* Finished calculating state PY_1I */ 


        /* For state SPACER_1I */ 
        /* setting first movement to score */ 
        score = GeneWise21_EXPL_MATRIX(mat,i-0,j-1,PY_1I) + mat->gp21->transition[GP21_PY2SPACER];   
        /* From state SPACER_1I to state SPACER_1I */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-0,j-1,SPACER_1I) + mat->gp21->transition[GP21_SPACER2SPACER];    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for SPACER_1I */ 
        /* Add any movement independant score and put away */ 
         score += mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)];   
         GeneWise21_EXPL_MATRIX(mat,i,j,SPACER_1I) = score;  


        /* Finished calculating state SPACER_1I */ 


        /* For state CENTRAL_2I */ 
        /* setting first movement to score */ 
        score = GeneWise21_EXPL_MATRIX(mat,i-0,j-10,MATCH) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_MATCH2INSERT]);   
        /* From state INSERT to state CENTRAL_2I */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-0,j-10,INSERT) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_INSERT2INSERT]);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state CENTRAL_2I */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-0,j-10,INSERT) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_DELETE2INSERT]);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state CENTRAL_2I to state CENTRAL_2I */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-0,j-1,CENTRAL_2I) + (mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CENTRAL2CENTRAL]);  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for CENTRAL_2I */ 
        /* Add any movement independant score and put away */ 
         GeneWise21_EXPL_MATRIX(mat,i,j,CENTRAL_2I) = score; 


        /* Finished calculating state CENTRAL_2I */ 


        /* For state PY_2I */ 
        /* setting first movement to score */ 
        score = GeneWise21_EXPL_MATRIX(mat,i-0,j-1,CENTRAL_2I) + mat->gp21->transition[GP21_CENTRAL2PY];     
        /* From state PY_2I to state PY_2I */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-0,j-1,PY_2I) + mat->gp21->transition[GP21_PY2PY];    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for PY_2I */ 
        /* Add any movement independant score and put away */ 
         score += mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)];   
         GeneWise21_EXPL_MATRIX(mat,i,j,PY_2I) = score;  


        /* Finished calculating state PY_2I */ 


        /* For state SPACER_2I */ 
        /* setting first movement to score */ 
        score = GeneWise21_EXPL_MATRIX(mat,i-0,j-1,PY_2I) + mat->gp21->transition[GP21_PY2SPACER];   
        /* From state SPACER_2I to state SPACER_2I */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-0,j-1,SPACER_2I) + mat->gp21->transition[GP21_SPACER2SPACER];    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for SPACER_2I */ 
        /* Add any movement independant score and put away */ 
         score += mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)];   
         GeneWise21_EXPL_MATRIX(mat,i,j,SPACER_2I) = score;  


        /* Finished calculating state SPACER_2I */ 
        }  


      /* Special state START has no special to special movements */ 


      /* Special state END has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GeneWise21_EXPL_SPECIAL(mat,0,j,END);  


      /* Source AFTER_RND is a special source for END */ 
      /* Has restricted position */ 
      if( j == mat->lenj-1 ) {  
        temp = GeneWise21_EXPL_SPECIAL(mat,0,j - 1,AFTER_RND) + (0) + (0);   
        if( temp > score )   
          score = temp;  
        }  


      /* Put back score... (now updated!) */ 
      GeneWise21_EXPL_SPECIAL(mat,0,j,END) = score;  
      /* Finished updating state END */ 




      /* Special state BEFORE_RND has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GeneWise21_EXPL_SPECIAL(mat,0,j,BEFORE_RND);   


      /* Source START is a special source for BEFORE_RND */ 
      temp = GeneWise21_EXPL_SPECIAL(mat,0,j - 1,START) + (mat->rndbase->base[CSEQ_GENOMIC_BASE(mat->target,j)]) + (0);  
      if( temp > score ) 
        score = temp;    


      /* Source BEFORE_RND is a special source for BEFORE_RND */ 
      temp = GeneWise21_EXPL_SPECIAL(mat,0,j - 1,BEFORE_RND) + ((mat->rndbase->base[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_RND2RND])) + (0);   
      if( temp > score ) 
        score = temp;    


      /* Source BEFORE_CODON is a special source for BEFORE_RND */ 
      temp = GeneWise21_EXPL_SPECIAL(mat,0,j - 1,BEFORE_CODON) + ((mat->rndbase->base[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CDS2RND])) + (0);     
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      GeneWise21_EXPL_SPECIAL(mat,0,j,BEFORE_RND) = score;   
      /* Finished updating state BEFORE_RND */ 




      /* Special state BEFORE_CODON has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GeneWise21_EXPL_SPECIAL(mat,0,j,BEFORE_CODON); 


      /* Source BEFORE_RND is a special source for BEFORE_CODON */ 
      temp = GeneWise21_EXPL_SPECIAL(mat,0,j - 3,BEFORE_RND) + ((mat->rndcodon->codon[CSEQ_GENOMIC_CODON(mat->target,j)]+mat->gp21->transition[GP21_RND2CDS])) + (0);    
      if( temp > score ) 
        score = temp;    


      /* Source BEFORE_CODON is a special source for BEFORE_CODON */ 
      temp = GeneWise21_EXPL_SPECIAL(mat,0,j - 3,BEFORE_CODON) + ((mat->rndcodon->codon[CSEQ_GENOMIC_CODON(mat->target,j)]+mat->gp21->transition[GP21_CDS2CDS])) + (0);  
      if( temp > score ) 
        score = temp;    


      /* Source BEFORE_SPACER is a special source for BEFORE_CODON */ 
      temp = GeneWise21_EXPL_SPECIAL(mat,0,j - 5,BEFORE_SPACER) + (((mat->gp21->transition[GP21_SPACER2CDS]+mat->rndcodon->codon[CSEQ_GENOMIC_CODON(mat->target,j)])+CSEQ_GENOMIC_3SS(mat->target,(j-3)))) + (0);    
      if( temp > score ) 
        score = temp;    


      /* Source BEFORE_SPACER is a special source for BEFORE_CODON */ 
      temp = GeneWise21_EXPL_SPECIAL(mat,0,j - 4,BEFORE_SPACER) + ((mat->gp21->transition[GP21_SPACER2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-2)))) + (0);  
      if( temp > score ) 
        score = temp;    


      /* Source BEFORE_SPACER is a special source for BEFORE_CODON */ 
      temp = GeneWise21_EXPL_SPECIAL(mat,0,j - 3,BEFORE_SPACER) + ((mat->gp21->transition[GP21_SPACER2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-1)))) + (0);  
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      GeneWise21_EXPL_SPECIAL(mat,0,j,BEFORE_CODON) = score; 
      /* Finished updating state BEFORE_CODON */ 




      /* Special state BEFORE_CENTRAL has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GeneWise21_EXPL_SPECIAL(mat,0,j,BEFORE_CENTRAL);   


      /* Source BEFORE_CENTRAL is a special source for BEFORE_CENTRAL */ 
      temp = GeneWise21_EXPL_SPECIAL(mat,0,j - 1,BEFORE_CENTRAL) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CENTRAL2CENTRAL])) + (0);   
      if( temp > score ) 
        score = temp;    


      /* Source BEFORE_CODON is a special source for BEFORE_CENTRAL */ 
      temp = GeneWise21_EXPL_SPECIAL(mat,0,j - 8,BEFORE_CODON) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))) + (0);     
      if( temp > score ) 
        score = temp;    


      /* Source BEFORE_CODON is a special source for BEFORE_CENTRAL */ 
      temp = GeneWise21_EXPL_SPECIAL(mat,0,j - 9,BEFORE_CODON) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))) + (0);     
      if( temp > score ) 
        score = temp;    


      /* Source BEFORE_CODON is a special source for BEFORE_CENTRAL */ 
      temp = GeneWise21_EXPL_SPECIAL(mat,0,j - 10,BEFORE_CODON) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))) + (0);    
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      GeneWise21_EXPL_SPECIAL(mat,0,j,BEFORE_CENTRAL) = score;   
      /* Finished updating state BEFORE_CENTRAL */ 




      /* Special state BEFORE_PY_TRACT has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GeneWise21_EXPL_SPECIAL(mat,0,j,BEFORE_PY_TRACT);  


      /* Source BEFORE_CENTRAL is a special source for BEFORE_PY_TRACT */ 
      temp = GeneWise21_EXPL_SPECIAL(mat,0,j - 1,BEFORE_CENTRAL) + (mat->gp21->transition[GP21_CENTRAL2PY]) + (mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)]);     
      if( temp > score ) 
        score = temp;    


      /* Source BEFORE_PY_TRACT is a special source for BEFORE_PY_TRACT */ 
      temp = GeneWise21_EXPL_SPECIAL(mat,0,j - 1,BEFORE_PY_TRACT) + (mat->gp21->transition[GP21_PY2PY]) + (mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)]);     
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      GeneWise21_EXPL_SPECIAL(mat,0,j,BEFORE_PY_TRACT) = score;  
      /* Finished updating state BEFORE_PY_TRACT */ 




      /* Special state BEFORE_SPACER has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GeneWise21_EXPL_SPECIAL(mat,0,j,BEFORE_SPACER);    


      /* Source BEFORE_PY_TRACT is a special source for BEFORE_SPACER */ 
      temp = GeneWise21_EXPL_SPECIAL(mat,0,j - 1,BEFORE_PY_TRACT) + (mat->gp21->transition[GP21_PY2SPACER]) + (mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)]);     
      if( temp > score ) 
        score = temp;    


      /* Source BEFORE_SPACER is a special source for BEFORE_SPACER */ 
      temp = GeneWise21_EXPL_SPECIAL(mat,0,j - 1,BEFORE_SPACER) + (mat->gp21->transition[GP21_SPACER2SPACER]) + (mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)]);   
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      GeneWise21_EXPL_SPECIAL(mat,0,j,BEFORE_SPACER) = score;    
      /* Finished updating state BEFORE_SPACER */ 




      /* Special state AFTER_RND has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GeneWise21_EXPL_SPECIAL(mat,0,j,AFTER_RND);    


      /* Source AFTER_RND is a special source for AFTER_RND */ 
      temp = GeneWise21_EXPL_SPECIAL(mat,0,j - 1,AFTER_RND) + ((mat->rndbase->base[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_RND2RND])) + (0);    
      if( temp > score ) 
        score = temp;    


      /* Source AFTER_CODON is a special source for AFTER_RND */ 
      temp = GeneWise21_EXPL_SPECIAL(mat,0,j - 1,AFTER_CODON) + ((mat->rndbase->base[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CDS2RND])) + (0);  
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      GeneWise21_EXPL_SPECIAL(mat,0,j,AFTER_RND) = score;    
      /* Finished updating state AFTER_RND */ 




      /* Special state AFTER_CODON has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GeneWise21_EXPL_SPECIAL(mat,0,j,AFTER_CODON);  


      /* Source AFTER_RND is a special source for AFTER_CODON */ 
      temp = GeneWise21_EXPL_SPECIAL(mat,0,j - 3,AFTER_RND) + ((mat->rndcodon->codon[CSEQ_GENOMIC_CODON(mat->target,j)]+mat->gp21->transition[GP21_RND2CDS])) + (0);     
      if( temp > score ) 
        score = temp;    


      /* Source MATCH for state AFTER_CODON is not special... already calculated */ 
      /* Source INSERT for state AFTER_CODON is not special... already calculated */ 
      /* Source DELETE for state AFTER_CODON is not special... already calculated */ 
      /* Source AFTER_CODON is a special source for AFTER_CODON */ 
      temp = GeneWise21_EXPL_SPECIAL(mat,0,j - 3,AFTER_CODON) + ((mat->rndcodon->codon[CSEQ_GENOMIC_CODON(mat->target,j)]+mat->gp21->transition[GP21_CDS2CDS])) + (0);   
      if( temp > score ) 
        score = temp;    


      /* Source AFTER_SPACER is a special source for AFTER_CODON */ 
      temp = GeneWise21_EXPL_SPECIAL(mat,0,j - 5,AFTER_SPACER) + (((mat->gp21->transition[GP21_SPACER2CDS]+mat->rndcodon->codon[CSEQ_GENOMIC_CODON(mat->target,j)])+CSEQ_GENOMIC_3SS(mat->target,(j-3)))) + (0);     
      if( temp > score ) 
        score = temp;    


      /* Source AFTER_SPACER is a special source for AFTER_CODON */ 
      temp = GeneWise21_EXPL_SPECIAL(mat,0,j - 4,AFTER_SPACER) + ((mat->gp21->transition[GP21_SPACER2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-2)))) + (0);   
      if( temp > score ) 
        score = temp;    


      /* Source AFTER_SPACER is a special source for AFTER_CODON */ 
      temp = GeneWise21_EXPL_SPECIAL(mat,0,j - 3,AFTER_SPACER) + ((mat->gp21->transition[GP21_SPACER2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-1)))) + (0);   
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      GeneWise21_EXPL_SPECIAL(mat,0,j,AFTER_CODON) = score;  
      /* Finished updating state AFTER_CODON */ 




      /* Special state AFTER_CENTRAL has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GeneWise21_EXPL_SPECIAL(mat,0,j,AFTER_CENTRAL);    


      /* Source AFTER_CENTRAL is a special source for AFTER_CENTRAL */ 
      temp = GeneWise21_EXPL_SPECIAL(mat,0,j - 1,AFTER_CENTRAL) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CENTRAL2CENTRAL])) + (0);    
      if( temp > score ) 
        score = temp;    


      /* Source AFTER_CODON is a special source for AFTER_CENTRAL */ 
      temp = GeneWise21_EXPL_SPECIAL(mat,0,j - 8,AFTER_CODON) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))) + (0);  
      if( temp > score ) 
        score = temp;    


      /* Source AFTER_CODON is a special source for AFTER_CENTRAL */ 
      temp = GeneWise21_EXPL_SPECIAL(mat,0,j - 9,AFTER_CODON) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))) + (0);  
      if( temp > score ) 
        score = temp;    


      /* Source AFTER_CODON is a special source for AFTER_CENTRAL */ 
      temp = GeneWise21_EXPL_SPECIAL(mat,0,j - 10,AFTER_CODON) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))) + (0);     
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      GeneWise21_EXPL_SPECIAL(mat,0,j,AFTER_CENTRAL) = score;    
      /* Finished updating state AFTER_CENTRAL */ 




      /* Special state AFTER_PY_TRACT has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GeneWise21_EXPL_SPECIAL(mat,0,j,AFTER_PY_TRACT);   


      /* Source AFTER_CENTRAL is a special source for AFTER_PY_TRACT */ 
      temp = GeneWise21_EXPL_SPECIAL(mat,0,j - 1,AFTER_CENTRAL) + (mat->gp21->transition[GP21_CENTRAL2PY]) + (mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)]);  
      if( temp > score ) 
        score = temp;    


      /* Source AFTER_PY_TRACT is a special source for AFTER_PY_TRACT */ 
      temp = GeneWise21_EXPL_SPECIAL(mat,0,j - 1,AFTER_PY_TRACT) + (mat->gp21->transition[GP21_PY2PY]) + (mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)]);  
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      GeneWise21_EXPL_SPECIAL(mat,0,j,AFTER_PY_TRACT) = score;   
      /* Finished updating state AFTER_PY_TRACT */ 




      /* Special state AFTER_SPACER has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GeneWise21_EXPL_SPECIAL(mat,0,j,AFTER_SPACER); 


      /* Source AFTER_PY_TRACT is a special source for AFTER_SPACER */ 
      temp = GeneWise21_EXPL_SPECIAL(mat,0,j - 1,AFTER_PY_TRACT) + (mat->gp21->transition[GP21_PY2SPACER]) + (mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)]);  
      if( temp > score ) 
        score = temp;    


      /* Source AFTER_SPACER is a special source for AFTER_SPACER */ 
      temp = GeneWise21_EXPL_SPECIAL(mat,0,j - 1,AFTER_SPACER) + (mat->gp21->transition[GP21_SPACER2SPACER]) + (mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)]);    
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      GeneWise21_EXPL_SPECIAL(mat,0,j,AFTER_SPACER) = score; 
      /* Finished updating state AFTER_SPACER */ 


      }  
    stop_reporting();    
    return TRUE;     
}    


/* Function:  calculate_dpenv_GeneWise21(mat,dpenv)
 *
 * Descrip:    This function calculates the GeneWise21 matrix when in explicit mode, subject to the envelope
 *
 *
 * Arg:          mat [UNKN ] GeneWise21 which contains explicit basematrix memory [GeneWise21 *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean calculate_dpenv_GeneWise21(GeneWise21 * mat,DPEnvelope * dpenv) 
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
      warn("in calculate_GeneWise21, passed a non Explicit matrix type, cannot calculate!"); 
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
      for(i=1;i<mat->leni;i++)   {  
        GeneWise21_EXPL_MATRIX(mat,i,j,MATCH) = NEGI;    
        GeneWise21_EXPL_MATRIX(mat,i,j,INSERT) = NEGI;   
        GeneWise21_EXPL_MATRIX(mat,i,j,DELETE) = NEGI;   
        GeneWise21_EXPL_MATRIX(mat,i,j,CENTRAL_0M) = NEGI;   
        GeneWise21_EXPL_MATRIX(mat,i,j,PY_0M) = NEGI;    
        GeneWise21_EXPL_MATRIX(mat,i,j,SPACER_0M) = NEGI;    
        GeneWise21_EXPL_MATRIX(mat,i,j,CENTRAL_1M) = NEGI;   
        GeneWise21_EXPL_MATRIX(mat,i,j,PY_1M) = NEGI;    
        GeneWise21_EXPL_MATRIX(mat,i,j,SPACER_1M) = NEGI;    
        GeneWise21_EXPL_MATRIX(mat,i,j,CENTRAL_2M) = NEGI;   
        GeneWise21_EXPL_MATRIX(mat,i,j,PY_2M) = NEGI;    
        GeneWise21_EXPL_MATRIX(mat,i,j,SPACER_2M) = NEGI;    
        GeneWise21_EXPL_MATRIX(mat,i,j,CENTRAL_0I) = NEGI;   
        GeneWise21_EXPL_MATRIX(mat,i,j,PY_0I) = NEGI;    
        GeneWise21_EXPL_MATRIX(mat,i,j,SPACER_0I) = NEGI;    
        GeneWise21_EXPL_MATRIX(mat,i,j,CENTRAL_1I) = NEGI;   
        GeneWise21_EXPL_MATRIX(mat,i,j,PY_1I) = NEGI;    
        GeneWise21_EXPL_MATRIX(mat,i,j,SPACER_1I) = NEGI;    
        GeneWise21_EXPL_MATRIX(mat,i,j,CENTRAL_2I) = NEGI;   
        GeneWise21_EXPL_MATRIX(mat,i,j,PY_2I) = NEGI;    
        GeneWise21_EXPL_MATRIX(mat,i,j,SPACER_2I) = NEGI;    
        }  
      }  
    for(j=-10;j<mat->lenj;j++)   {  
      GeneWise21_EXPL_SPECIAL(mat,i,j,START) = 0;    
      GeneWise21_EXPL_SPECIAL(mat,i,j,END) = NEGI;   
      GeneWise21_EXPL_SPECIAL(mat,i,j,BEFORE_RND) = NEGI;    
      GeneWise21_EXPL_SPECIAL(mat,i,j,BEFORE_CODON) = NEGI;  
      GeneWise21_EXPL_SPECIAL(mat,i,j,BEFORE_CENTRAL) = NEGI;    
      GeneWise21_EXPL_SPECIAL(mat,i,j,BEFORE_PY_TRACT) = NEGI;   
      GeneWise21_EXPL_SPECIAL(mat,i,j,BEFORE_SPACER) = NEGI; 
      GeneWise21_EXPL_SPECIAL(mat,i,j,AFTER_RND) = NEGI; 
      GeneWise21_EXPL_SPECIAL(mat,i,j,AFTER_CODON) = NEGI;   
      GeneWise21_EXPL_SPECIAL(mat,i,j,AFTER_CENTRAL) = NEGI; 
      GeneWise21_EXPL_SPECIAL(mat,i,j,AFTER_PY_TRACT) = NEGI;    
      GeneWise21_EXPL_SPECIAL(mat,i,j,AFTER_SPACER) = NEGI;  
      }  


    start_reporting("GeneWise21 Matrix calculation: ");  
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
          GeneWise21_EXPL_MATRIX(mat,i,j,MATCH) = NEGI;  
          GeneWise21_EXPL_MATRIX(mat,i,j,INSERT) = NEGI; 
          GeneWise21_EXPL_MATRIX(mat,i,j,DELETE) = NEGI; 
          GeneWise21_EXPL_MATRIX(mat,i,j,CENTRAL_0M) = NEGI; 
          GeneWise21_EXPL_MATRIX(mat,i,j,PY_0M) = NEGI;  
          GeneWise21_EXPL_MATRIX(mat,i,j,SPACER_0M) = NEGI;  
          GeneWise21_EXPL_MATRIX(mat,i,j,CENTRAL_1M) = NEGI; 
          GeneWise21_EXPL_MATRIX(mat,i,j,PY_1M) = NEGI;  
          GeneWise21_EXPL_MATRIX(mat,i,j,SPACER_1M) = NEGI;  
          GeneWise21_EXPL_MATRIX(mat,i,j,CENTRAL_2M) = NEGI; 
          GeneWise21_EXPL_MATRIX(mat,i,j,PY_2M) = NEGI;  
          GeneWise21_EXPL_MATRIX(mat,i,j,SPACER_2M) = NEGI;  
          GeneWise21_EXPL_MATRIX(mat,i,j,CENTRAL_0I) = NEGI; 
          GeneWise21_EXPL_MATRIX(mat,i,j,PY_0I) = NEGI;  
          GeneWise21_EXPL_MATRIX(mat,i,j,SPACER_0I) = NEGI;  
          GeneWise21_EXPL_MATRIX(mat,i,j,CENTRAL_1I) = NEGI; 
          GeneWise21_EXPL_MATRIX(mat,i,j,PY_1I) = NEGI;  
          GeneWise21_EXPL_MATRIX(mat,i,j,SPACER_1I) = NEGI;  
          GeneWise21_EXPL_MATRIX(mat,i,j,CENTRAL_2I) = NEGI; 
          GeneWise21_EXPL_MATRIX(mat,i,j,PY_2I) = NEGI;  
          GeneWise21_EXPL_MATRIX(mat,i,j,SPACER_2I) = NEGI;  
          continue;  
          }  


        if( num%1000 == 0 )  
          log_full_error(REPORT,0,"[%7d] Cells %2d%%%%",num,num*100/tot);    
        num++;   


        /* For state MATCH */ 
        /* setting first movement to score */ 
        score = GeneWise21_EXPL_MATRIX(mat,i-1,j-3,MATCH) + (mat->query->seg[i]->transition[GW_MATCH2MATCH]+mat->query->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]);   
        /* From state INSERT to state MATCH */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-1,j-3,INSERT) + (mat->query->seg[i]->transition[GW_INSERT2MATCH]+mat->query->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state MATCH */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-1,j-3,DELETE) + (mat->query->seg[i]->transition[GW_DELETE2MATCH]+mat->query->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state BEFORE_CODON to state MATCH */ 
        temp = GeneWise21_EXPL_SPECIAL(mat,i-1,j-3,BEFORE_CODON) + (mat->query->seg[i]->transition[GW_START2MATCH]+mat->query->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]);    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SPACER_0M to state MATCH */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-1,j-6,SPACER_0M) + (((mat->gp21->transition[GP21_SPACER2CDS]+mat->query->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)])+CSEQ_GENOMIC_3SS(mat->target,(j-3)))+mat->query->seg[i]->transition[GW_MATCH_BALANCE_3SS]);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state PY_0M to state MATCH */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-1,j-6,PY_0M) + (((mat->gp21->transition[GP21_PY2CDS]+mat->query->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)])+CSEQ_GENOMIC_3SS(mat->target,(j-3)))+mat->query->seg[i]->transition[GW_MATCH_BALANCE_3SS]);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SPACER_1M to state MATCH */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-1,j-5,SPACER_1M) + (mat->gp21->transition[GP21_SPACER2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-2)));     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state PY_1M to state MATCH */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-1,j-5,PY_1M) + (mat->gp21->transition[GP21_PY2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-2)));     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SPACER_2M to state MATCH */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-1,j-4,SPACER_2M) + (mat->gp21->transition[GP21_SPACER2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-1)));     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state PY_2M to state MATCH */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-1,j-4,PY_2M) + (mat->gp21->transition[GP21_PY2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-1)));     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-1,j-2,MATCH) + mat->gp21->transition[GP21_DELETE_1_BASE];    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-1,j-1,MATCH) + mat->gp21->transition[GP21_DELETE_2_BASE];    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-1,j-4,MATCH) + mat->gp21->transition[GP21_INSERT_1_BASE];    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-1,j-5,MATCH) + mat->gp21->transition[GP21_INSERT_2_BASE];    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for MATCH */ 
        /* Add any movement independant score and put away */ 
         score += CSEQ_GENOMIC_CDSPOT(mat->target,j);    
         GeneWise21_EXPL_MATRIX(mat,i,j,MATCH) = score;  


        /* state MATCH is a source for special AFTER_CODON */ 
        temp = score + (mat->query->seg[i]->transition[GW_MATCH2END]) + (0) ;    
        if( temp > GeneWise21_EXPL_SPECIAL(mat,i,j,AFTER_CODON) )    {  
          GeneWise21_EXPL_SPECIAL(mat,i,j,AFTER_CODON) = temp;   
          }  




        /* Finished calculating state MATCH */ 


        /* For state INSERT */ 
        /* setting first movement to score */ 
        score = GeneWise21_EXPL_MATRIX(mat,i-0,j-3,MATCH) + (mat->query->seg[i]->transition[GW_MATCH2INSERT]+mat->query->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)]);     
        /* From state INSERT to state INSERT */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-0,j-3,INSERT) + (mat->query->seg[i]->transition[GW_INSERT2INSERT]+mat->query->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)]);    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state INSERT */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-0,j-3,DELETE) + (mat->query->seg[i]->transition[GW_DELETE2INSERT]+mat->query->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)]);    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state BEFORE_CODON to state INSERT */ 
        temp = GeneWise21_EXPL_SPECIAL(mat,i-0,j-3,BEFORE_CODON) + (mat->query->seg[i]->transition[GW_START2INSERT]+mat->query->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)]);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SPACER_0I to state INSERT */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-1,j-6,SPACER_0I) + (((mat->gp21->transition[GP21_SPACER2CDS]+mat->query->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)])+CSEQ_GENOMIC_3SS(mat->target,(j-3)))+mat->query->seg[i]->transition[GW_MATCH_BALANCE_3SS]);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SPACER_1I to state INSERT */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-1,j-5,SPACER_1I) + (mat->gp21->transition[GP21_SPACER2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-2)));     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SPACER_2I to state INSERT */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-1,j-4,SPACER_2I) + (mat->gp21->transition[GP21_SPACER2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-1)));     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state PY_0I to state INSERT */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-1,j-6,PY_0I) + (((mat->gp21->transition[GP21_PY2CDS]+mat->query->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)])+CSEQ_GENOMIC_3SS(mat->target,(j-3)))+mat->query->seg[i]->transition[GW_MATCH_BALANCE_3SS]);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state PY_1I to state INSERT */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-1,j-5,PY_1I) + (mat->gp21->transition[GP21_PY2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-2)));     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state PY_2I to state INSERT */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-1,j-4,PY_2I) + (mat->gp21->transition[GP21_PY2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-1)));     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state INSERT */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-1,j-2,INSERT) + mat->gp21->transition[GP21_DELETE_1_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state INSERT */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-1,j-1,INSERT) + mat->gp21->transition[GP21_DELETE_2_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state INSERT */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-1,j-4,INSERT) + mat->gp21->transition[GP21_INSERT_1_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state INSERT */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-1,j-5,INSERT) + mat->gp21->transition[GP21_INSERT_2_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for INSERT */ 
        /* Add any movement independant score and put away */ 
         score += CSEQ_GENOMIC_CDSPOT(mat->target,j);    
         GeneWise21_EXPL_MATRIX(mat,i,j,INSERT) = score; 


        /* state INSERT is a source for special AFTER_CODON */ 
        temp = score + (mat->query->seg[i]->transition[GW_INSERT2END]) + (0) ;   
        if( temp > GeneWise21_EXPL_SPECIAL(mat,i,j,AFTER_CODON) )    {  
          GeneWise21_EXPL_SPECIAL(mat,i,j,AFTER_CODON) = temp;   
          }  




        /* Finished calculating state INSERT */ 


        /* For state DELETE */ 
        /* setting first movement to score */ 
        score = GeneWise21_EXPL_MATRIX(mat,i-1,j-0,MATCH) + mat->query->seg[i]->transition[GW_MATCH2DELETE];     
        /* From state INSERT to state DELETE */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-1,j-0,INSERT) + mat->query->seg[i]->transition[GW_INSERT2DELETE];    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state DELETE */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-1,j-0,DELETE) + mat->query->seg[i]->transition[GW_DELETE2DELETE];    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state BEFORE_CODON to state DELETE */ 
        temp = GeneWise21_EXPL_SPECIAL(mat,i-1,j-0,BEFORE_CODON) + mat->query->seg[i]->transition[GW_START2DELETE];  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for DELETE */ 
        /* Add any movement independant score and put away */ 
         GeneWise21_EXPL_MATRIX(mat,i,j,DELETE) = score; 


        /* state DELETE is a source for special AFTER_CODON */ 
        temp = score + (mat->query->seg[i]->transition[GW_DELETE2END]) + (0) ;   
        if( temp > GeneWise21_EXPL_SPECIAL(mat,i,j,AFTER_CODON) )    {  
          GeneWise21_EXPL_SPECIAL(mat,i,j,AFTER_CODON) = temp;   
          }  




        /* Finished calculating state DELETE */ 


        /* For state CENTRAL_0M */ 
        /* setting first movement to score */ 
        score = GeneWise21_EXPL_MATRIX(mat,i-0,j-8,MATCH) + (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_MATCH2MATCH])+mat->query->seg[i]->transition[GW_MATCH_BALANCE_5SS]);  
        /* From state INSERT to state CENTRAL_0M */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-0,j-8,INSERT) + (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_INSERT2MATCH])+mat->query->seg[i]->transition[GW_INSERT_BALANCE_5SS]);    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state CENTRAL_0M */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-0,j-8,DELETE) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_DELETE2MATCH]);    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state CENTRAL_0M to state CENTRAL_0M */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-0,j-1,CENTRAL_0M) + (mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CENTRAL2CENTRAL]);  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for CENTRAL_0M */ 
        /* Add any movement independant score and put away */ 
         GeneWise21_EXPL_MATRIX(mat,i,j,CENTRAL_0M) = score; 


        /* Finished calculating state CENTRAL_0M */ 


        /* For state PY_0M */ 
        /* setting first movement to score */ 
        score = GeneWise21_EXPL_MATRIX(mat,i-0,j-1,CENTRAL_0M) + mat->gp21->transition[GP21_CENTRAL2PY];     
        /* From state PY_0M to state PY_0M */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-0,j-1,PY_0M) + mat->gp21->transition[GP21_PY2PY];    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for PY_0M */ 
        /* Add any movement independant score and put away */ 
         score += mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)];   
         GeneWise21_EXPL_MATRIX(mat,i,j,PY_0M) = score;  


        /* Finished calculating state PY_0M */ 


        /* For state SPACER_0M */ 
        /* setting first movement to score */ 
        score = GeneWise21_EXPL_MATRIX(mat,i-0,j-1,PY_0M) + mat->gp21->transition[GP21_PY2SPACER];   
        /* From state SPACER_0M to state SPACER_0M */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-0,j-1,SPACER_0M) + mat->gp21->transition[GP21_SPACER2SPACER];    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for SPACER_0M */ 
        /* Add any movement independant score and put away */ 
         score += mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)];   
         GeneWise21_EXPL_MATRIX(mat,i,j,SPACER_0M) = score;  


        /* Finished calculating state SPACER_0M */ 


        /* For state CENTRAL_1M */ 
        /* setting first movement to score */ 
        score = GeneWise21_EXPL_MATRIX(mat,i-0,j-9,MATCH) + (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_MATCH2MATCH])+mat->query->seg[i]->transition[GW_MATCH_BALANCE_5SS]);  
        /* From state INSERT to state CENTRAL_1M */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-0,j-9,INSERT) + (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_INSERT2MATCH])+mat->query->seg[i]->transition[GW_MATCH_BALANCE_5SS]);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state CENTRAL_1M */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-0,j-9,DELETE) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_DELETE2MATCH]);    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state CENTRAL_1M to state CENTRAL_1M */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-0,j-1,CENTRAL_1M) + (mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CENTRAL2CENTRAL]);  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for CENTRAL_1M */ 
        /* Add any movement independant score and put away */ 
         GeneWise21_EXPL_MATRIX(mat,i,j,CENTRAL_1M) = score; 


        /* Finished calculating state CENTRAL_1M */ 


        /* For state PY_1M */ 
        /* setting first movement to score */ 
        score = GeneWise21_EXPL_MATRIX(mat,i-0,j-1,CENTRAL_1M) + mat->gp21->transition[GP21_CENTRAL2PY];     
        /* From state PY_1M to state PY_1M */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-0,j-1,PY_1M) + mat->gp21->transition[GP21_PY2PY];    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for PY_1M */ 
        /* Add any movement independant score and put away */ 
         score += mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)];   
         GeneWise21_EXPL_MATRIX(mat,i,j,PY_1M) = score;  


        /* Finished calculating state PY_1M */ 


        /* For state SPACER_1M */ 
        /* setting first movement to score */ 
        score = GeneWise21_EXPL_MATRIX(mat,i-0,j-1,PY_1M) + mat->gp21->transition[GP21_PY2SPACER];   
        /* From state SPACER_1M to state SPACER_1M */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-0,j-1,SPACER_1M) + mat->gp21->transition[GP21_SPACER2SPACER];    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for SPACER_1M */ 
        /* Add any movement independant score and put away */ 
         score += mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)];   
         GeneWise21_EXPL_MATRIX(mat,i,j,SPACER_1M) = score;  


        /* Finished calculating state SPACER_1M */ 


        /* For state CENTRAL_2M */ 
        /* setting first movement to score */ 
        score = GeneWise21_EXPL_MATRIX(mat,i-0,j-10,MATCH) + (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_MATCH2MATCH])+mat->query->seg[i]->transition[GW_MATCH_BALANCE_5SS]);     
        /* From state INSERT to state CENTRAL_2M */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-0,j-10,INSERT) + (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_INSERT2MATCH])+mat->query->seg[i]->transition[GW_MATCH_BALANCE_5SS]);    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state CENTRAL_2M */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-0,j-10,DELETE) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_DELETE2MATCH]);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state CENTRAL_2M to state CENTRAL_2M */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-0,j-1,CENTRAL_2M) + (mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CENTRAL2CENTRAL]);  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for CENTRAL_2M */ 
        /* Add any movement independant score and put away */ 
         GeneWise21_EXPL_MATRIX(mat,i,j,CENTRAL_2M) = score; 


        /* Finished calculating state CENTRAL_2M */ 


        /* For state PY_2M */ 
        /* setting first movement to score */ 
        score = GeneWise21_EXPL_MATRIX(mat,i-0,j-1,CENTRAL_2M) + mat->gp21->transition[GP21_CENTRAL2PY];     
        /* From state PY_2M to state PY_2M */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-0,j-1,PY_2M) + mat->gp21->transition[GP21_PY2PY];    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for PY_2M */ 
        /* Add any movement independant score and put away */ 
         score += mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)];   
         GeneWise21_EXPL_MATRIX(mat,i,j,PY_2M) = score;  


        /* Finished calculating state PY_2M */ 


        /* For state SPACER_2M */ 
        /* setting first movement to score */ 
        score = GeneWise21_EXPL_MATRIX(mat,i-0,j-1,PY_2M) + mat->gp21->transition[GP21_PY2SPACER];   
        /* From state SPACER_2M to state SPACER_2M */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-0,j-1,SPACER_2M) + mat->gp21->transition[GP21_SPACER2SPACER];    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for SPACER_2M */ 
        /* Add any movement independant score and put away */ 
         score += mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)];   
         GeneWise21_EXPL_MATRIX(mat,i,j,SPACER_2M) = score;  


        /* Finished calculating state SPACER_2M */ 


        /* For state CENTRAL_0I */ 
        /* setting first movement to score */ 
        score = GeneWise21_EXPL_MATRIX(mat,i-0,j-8,MATCH) + (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_MATCH2INSERT])+mat->query->seg[i]->transition[GW_MATCH_BALANCE_5SS]);     
        /* From state INSERT to state CENTRAL_0I */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-0,j-8,INSERT) + (((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_INSERT2INSERT])+mat->query->seg[i]->transition[GW_INSERT_BALANCE_5SS]);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state CENTRAL_0I */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-0,j-8,DELETE) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_DELETE2INSERT]);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state CENTRAL_0I to state CENTRAL_0I */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-0,j-1,CENTRAL_0I) + (mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CENTRAL2CENTRAL]);  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for CENTRAL_0I */ 
        /* Add any movement independant score and put away */ 
         GeneWise21_EXPL_MATRIX(mat,i,j,CENTRAL_0I) = score; 


        /* Finished calculating state CENTRAL_0I */ 


        /* For state PY_0I */ 
        /* setting first movement to score */ 
        score = GeneWise21_EXPL_MATRIX(mat,i-0,j-1,CENTRAL_0I) + mat->gp21->transition[GP21_CENTRAL2PY];     
        /* From state PY_0I to state PY_0I */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-0,j-1,PY_0I) + mat->gp21->transition[GP21_PY2PY];    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for PY_0I */ 
        /* Add any movement independant score and put away */ 
         score += mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)];   
         GeneWise21_EXPL_MATRIX(mat,i,j,PY_0I) = score;  


        /* Finished calculating state PY_0I */ 


        /* For state SPACER_0I */ 
        /* setting first movement to score */ 
        score = GeneWise21_EXPL_MATRIX(mat,i-0,j-1,PY_0I) + mat->gp21->transition[GP21_PY2SPACER];   
        /* From state SPACER_0I to state SPACER_0I */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-0,j-1,SPACER_0I) + mat->gp21->transition[GP21_SPACER2SPACER];    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for SPACER_0I */ 
        /* Add any movement independant score and put away */ 
         score += mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)];   
         GeneWise21_EXPL_MATRIX(mat,i,j,SPACER_0I) = score;  


        /* Finished calculating state SPACER_0I */ 


        /* For state CENTRAL_1I */ 
        /* setting first movement to score */ 
        score = GeneWise21_EXPL_MATRIX(mat,i-0,j-9,MATCH) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_MATCH2INSERT]);    
        /* From state INSERT to state CENTRAL_1I */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-0,j-9,INSERT) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_INSERT2INSERT]);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state CENTRAL_1I */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-0,j-9,DELETE) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_DELETE2INSERT]);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state CENTRAL_1I to state CENTRAL_1I */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-0,j-1,CENTRAL_1I) + (mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CENTRAL2CENTRAL]);  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for CENTRAL_1I */ 
        /* Add any movement independant score and put away */ 
         GeneWise21_EXPL_MATRIX(mat,i,j,CENTRAL_1I) = score; 


        /* Finished calculating state CENTRAL_1I */ 


        /* For state PY_1I */ 
        /* setting first movement to score */ 
        score = GeneWise21_EXPL_MATRIX(mat,i-0,j-1,CENTRAL_1I) + mat->gp21->transition[GP21_CENTRAL2PY];     
        /* From state PY_1I to state PY_1I */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-0,j-1,PY_1I) + mat->gp21->transition[GP21_PY2PY];    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for PY_1I */ 
        /* Add any movement independant score and put away */ 
         score += mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)];   
         GeneWise21_EXPL_MATRIX(mat,i,j,PY_1I) = score;  


        /* Finished calculating state PY_1I */ 


        /* For state SPACER_1I */ 
        /* setting first movement to score */ 
        score = GeneWise21_EXPL_MATRIX(mat,i-0,j-1,PY_1I) + mat->gp21->transition[GP21_PY2SPACER];   
        /* From state SPACER_1I to state SPACER_1I */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-0,j-1,SPACER_1I) + mat->gp21->transition[GP21_SPACER2SPACER];    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for SPACER_1I */ 
        /* Add any movement independant score and put away */ 
         score += mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)];   
         GeneWise21_EXPL_MATRIX(mat,i,j,SPACER_1I) = score;  


        /* Finished calculating state SPACER_1I */ 


        /* For state CENTRAL_2I */ 
        /* setting first movement to score */ 
        score = GeneWise21_EXPL_MATRIX(mat,i-0,j-10,MATCH) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_MATCH2INSERT]);   
        /* From state INSERT to state CENTRAL_2I */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-0,j-10,INSERT) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_INSERT2INSERT]);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state CENTRAL_2I */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-0,j-10,INSERT) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->seg[i]->transition[GW_DELETE2INSERT]);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state CENTRAL_2I to state CENTRAL_2I */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-0,j-1,CENTRAL_2I) + (mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CENTRAL2CENTRAL]);  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for CENTRAL_2I */ 
        /* Add any movement independant score and put away */ 
         GeneWise21_EXPL_MATRIX(mat,i,j,CENTRAL_2I) = score; 


        /* Finished calculating state CENTRAL_2I */ 


        /* For state PY_2I */ 
        /* setting first movement to score */ 
        score = GeneWise21_EXPL_MATRIX(mat,i-0,j-1,CENTRAL_2I) + mat->gp21->transition[GP21_CENTRAL2PY];     
        /* From state PY_2I to state PY_2I */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-0,j-1,PY_2I) + mat->gp21->transition[GP21_PY2PY];    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for PY_2I */ 
        /* Add any movement independant score and put away */ 
         score += mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)];   
         GeneWise21_EXPL_MATRIX(mat,i,j,PY_2I) = score;  


        /* Finished calculating state PY_2I */ 


        /* For state SPACER_2I */ 
        /* setting first movement to score */ 
        score = GeneWise21_EXPL_MATRIX(mat,i-0,j-1,PY_2I) + mat->gp21->transition[GP21_PY2SPACER];   
        /* From state SPACER_2I to state SPACER_2I */ 
        temp = GeneWise21_EXPL_MATRIX(mat,i-0,j-1,SPACER_2I) + mat->gp21->transition[GP21_SPACER2SPACER];    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for SPACER_2I */ 
        /* Add any movement independant score and put away */ 
         score += mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)];   
         GeneWise21_EXPL_MATRIX(mat,i,j,SPACER_2I) = score;  


        /* Finished calculating state SPACER_2I */ 
        }  


      /* Special state START has no special to special movements */ 


      /* Special state END has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GeneWise21_EXPL_SPECIAL(mat,0,j,END);  


      /* Source AFTER_RND is a special source for END */ 
      /* Has restricted position */ 
      if( j == mat->lenj-1 ) {  
        temp = GeneWise21_EXPL_SPECIAL(mat,0,j - 1,AFTER_RND) + (0) + (0);   
        if( temp > score )   
          score = temp;  
        }  


      /* Put back score... (now updated!) */ 
      GeneWise21_EXPL_SPECIAL(mat,0,j,END) = score;  
      /* Finished updating state END */ 




      /* Special state BEFORE_RND has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GeneWise21_EXPL_SPECIAL(mat,0,j,BEFORE_RND);   


      /* Source START is a special source for BEFORE_RND */ 
      temp = GeneWise21_EXPL_SPECIAL(mat,0,j - 1,START) + (mat->rndbase->base[CSEQ_GENOMIC_BASE(mat->target,j)]) + (0);  
      if( temp > score ) 
        score = temp;    


      /* Source BEFORE_RND is a special source for BEFORE_RND */ 
      temp = GeneWise21_EXPL_SPECIAL(mat,0,j - 1,BEFORE_RND) + ((mat->rndbase->base[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_RND2RND])) + (0);   
      if( temp > score ) 
        score = temp;    


      /* Source BEFORE_CODON is a special source for BEFORE_RND */ 
      temp = GeneWise21_EXPL_SPECIAL(mat,0,j - 1,BEFORE_CODON) + ((mat->rndbase->base[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CDS2RND])) + (0);     
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      GeneWise21_EXPL_SPECIAL(mat,0,j,BEFORE_RND) = score;   
      /* Finished updating state BEFORE_RND */ 




      /* Special state BEFORE_CODON has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GeneWise21_EXPL_SPECIAL(mat,0,j,BEFORE_CODON); 


      /* Source BEFORE_RND is a special source for BEFORE_CODON */ 
      temp = GeneWise21_EXPL_SPECIAL(mat,0,j - 3,BEFORE_RND) + ((mat->rndcodon->codon[CSEQ_GENOMIC_CODON(mat->target,j)]+mat->gp21->transition[GP21_RND2CDS])) + (0);    
      if( temp > score ) 
        score = temp;    


      /* Source BEFORE_CODON is a special source for BEFORE_CODON */ 
      temp = GeneWise21_EXPL_SPECIAL(mat,0,j - 3,BEFORE_CODON) + ((mat->rndcodon->codon[CSEQ_GENOMIC_CODON(mat->target,j)]+mat->gp21->transition[GP21_CDS2CDS])) + (0);  
      if( temp > score ) 
        score = temp;    


      /* Source BEFORE_SPACER is a special source for BEFORE_CODON */ 
      temp = GeneWise21_EXPL_SPECIAL(mat,0,j - 5,BEFORE_SPACER) + (((mat->gp21->transition[GP21_SPACER2CDS]+mat->rndcodon->codon[CSEQ_GENOMIC_CODON(mat->target,j)])+CSEQ_GENOMIC_3SS(mat->target,(j-3)))) + (0);    
      if( temp > score ) 
        score = temp;    


      /* Source BEFORE_SPACER is a special source for BEFORE_CODON */ 
      temp = GeneWise21_EXPL_SPECIAL(mat,0,j - 4,BEFORE_SPACER) + ((mat->gp21->transition[GP21_SPACER2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-2)))) + (0);  
      if( temp > score ) 
        score = temp;    


      /* Source BEFORE_SPACER is a special source for BEFORE_CODON */ 
      temp = GeneWise21_EXPL_SPECIAL(mat,0,j - 3,BEFORE_SPACER) + ((mat->gp21->transition[GP21_SPACER2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-1)))) + (0);  
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      GeneWise21_EXPL_SPECIAL(mat,0,j,BEFORE_CODON) = score; 
      /* Finished updating state BEFORE_CODON */ 




      /* Special state BEFORE_CENTRAL has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GeneWise21_EXPL_SPECIAL(mat,0,j,BEFORE_CENTRAL);   


      /* Source BEFORE_CENTRAL is a special source for BEFORE_CENTRAL */ 
      temp = GeneWise21_EXPL_SPECIAL(mat,0,j - 1,BEFORE_CENTRAL) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CENTRAL2CENTRAL])) + (0);   
      if( temp > score ) 
        score = temp;    


      /* Source BEFORE_CODON is a special source for BEFORE_CENTRAL */ 
      temp = GeneWise21_EXPL_SPECIAL(mat,0,j - 8,BEFORE_CODON) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))) + (0);     
      if( temp > score ) 
        score = temp;    


      /* Source BEFORE_CODON is a special source for BEFORE_CENTRAL */ 
      temp = GeneWise21_EXPL_SPECIAL(mat,0,j - 9,BEFORE_CODON) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))) + (0);     
      if( temp > score ) 
        score = temp;    


      /* Source BEFORE_CODON is a special source for BEFORE_CENTRAL */ 
      temp = GeneWise21_EXPL_SPECIAL(mat,0,j - 10,BEFORE_CODON) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))) + (0);    
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      GeneWise21_EXPL_SPECIAL(mat,0,j,BEFORE_CENTRAL) = score;   
      /* Finished updating state BEFORE_CENTRAL */ 




      /* Special state BEFORE_PY_TRACT has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GeneWise21_EXPL_SPECIAL(mat,0,j,BEFORE_PY_TRACT);  


      /* Source BEFORE_CENTRAL is a special source for BEFORE_PY_TRACT */ 
      temp = GeneWise21_EXPL_SPECIAL(mat,0,j - 1,BEFORE_CENTRAL) + (mat->gp21->transition[GP21_CENTRAL2PY]) + (mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)]);     
      if( temp > score ) 
        score = temp;    


      /* Source BEFORE_PY_TRACT is a special source for BEFORE_PY_TRACT */ 
      temp = GeneWise21_EXPL_SPECIAL(mat,0,j - 1,BEFORE_PY_TRACT) + (mat->gp21->transition[GP21_PY2PY]) + (mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)]);     
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      GeneWise21_EXPL_SPECIAL(mat,0,j,BEFORE_PY_TRACT) = score;  
      /* Finished updating state BEFORE_PY_TRACT */ 




      /* Special state BEFORE_SPACER has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GeneWise21_EXPL_SPECIAL(mat,0,j,BEFORE_SPACER);    


      /* Source BEFORE_PY_TRACT is a special source for BEFORE_SPACER */ 
      temp = GeneWise21_EXPL_SPECIAL(mat,0,j - 1,BEFORE_PY_TRACT) + (mat->gp21->transition[GP21_PY2SPACER]) + (mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)]);     
      if( temp > score ) 
        score = temp;    


      /* Source BEFORE_SPACER is a special source for BEFORE_SPACER */ 
      temp = GeneWise21_EXPL_SPECIAL(mat,0,j - 1,BEFORE_SPACER) + (mat->gp21->transition[GP21_SPACER2SPACER]) + (mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)]);   
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      GeneWise21_EXPL_SPECIAL(mat,0,j,BEFORE_SPACER) = score;    
      /* Finished updating state BEFORE_SPACER */ 




      /* Special state AFTER_RND has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GeneWise21_EXPL_SPECIAL(mat,0,j,AFTER_RND);    


      /* Source AFTER_RND is a special source for AFTER_RND */ 
      temp = GeneWise21_EXPL_SPECIAL(mat,0,j - 1,AFTER_RND) + ((mat->rndbase->base[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_RND2RND])) + (0);    
      if( temp > score ) 
        score = temp;    


      /* Source AFTER_CODON is a special source for AFTER_RND */ 
      temp = GeneWise21_EXPL_SPECIAL(mat,0,j - 1,AFTER_CODON) + ((mat->rndbase->base[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CDS2RND])) + (0);  
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      GeneWise21_EXPL_SPECIAL(mat,0,j,AFTER_RND) = score;    
      /* Finished updating state AFTER_RND */ 




      /* Special state AFTER_CODON has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GeneWise21_EXPL_SPECIAL(mat,0,j,AFTER_CODON);  


      /* Source AFTER_RND is a special source for AFTER_CODON */ 
      temp = GeneWise21_EXPL_SPECIAL(mat,0,j - 3,AFTER_RND) + ((mat->rndcodon->codon[CSEQ_GENOMIC_CODON(mat->target,j)]+mat->gp21->transition[GP21_RND2CDS])) + (0);     
      if( temp > score ) 
        score = temp;    


      /* Source MATCH for state AFTER_CODON is not special... already calculated */ 
      /* Source INSERT for state AFTER_CODON is not special... already calculated */ 
      /* Source DELETE for state AFTER_CODON is not special... already calculated */ 
      /* Source AFTER_CODON is a special source for AFTER_CODON */ 
      temp = GeneWise21_EXPL_SPECIAL(mat,0,j - 3,AFTER_CODON) + ((mat->rndcodon->codon[CSEQ_GENOMIC_CODON(mat->target,j)]+mat->gp21->transition[GP21_CDS2CDS])) + (0);   
      if( temp > score ) 
        score = temp;    


      /* Source AFTER_SPACER is a special source for AFTER_CODON */ 
      temp = GeneWise21_EXPL_SPECIAL(mat,0,j - 5,AFTER_SPACER) + (((mat->gp21->transition[GP21_SPACER2CDS]+mat->rndcodon->codon[CSEQ_GENOMIC_CODON(mat->target,j)])+CSEQ_GENOMIC_3SS(mat->target,(j-3)))) + (0);     
      if( temp > score ) 
        score = temp;    


      /* Source AFTER_SPACER is a special source for AFTER_CODON */ 
      temp = GeneWise21_EXPL_SPECIAL(mat,0,j - 4,AFTER_SPACER) + ((mat->gp21->transition[GP21_SPACER2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-2)))) + (0);   
      if( temp > score ) 
        score = temp;    


      /* Source AFTER_SPACER is a special source for AFTER_CODON */ 
      temp = GeneWise21_EXPL_SPECIAL(mat,0,j - 3,AFTER_SPACER) + ((mat->gp21->transition[GP21_SPACER2CDS]+CSEQ_GENOMIC_3SS(mat->target,(j-1)))) + (0);   
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      GeneWise21_EXPL_SPECIAL(mat,0,j,AFTER_CODON) = score;  
      /* Finished updating state AFTER_CODON */ 




      /* Special state AFTER_CENTRAL has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GeneWise21_EXPL_SPECIAL(mat,0,j,AFTER_CENTRAL);    


      /* Source AFTER_CENTRAL is a special source for AFTER_CENTRAL */ 
      temp = GeneWise21_EXPL_SPECIAL(mat,0,j - 1,AFTER_CENTRAL) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp21->transition[GP21_CENTRAL2CENTRAL])) + (0);    
      if( temp > score ) 
        score = temp;    


      /* Source AFTER_CODON is a special source for AFTER_CENTRAL */ 
      temp = GeneWise21_EXPL_SPECIAL(mat,0,j - 8,AFTER_CODON) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))) + (0);  
      if( temp > score ) 
        score = temp;    


      /* Source AFTER_CODON is a special source for AFTER_CENTRAL */ 
      temp = GeneWise21_EXPL_SPECIAL(mat,0,j - 9,AFTER_CODON) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))) + (0);  
      if( temp > score ) 
        score = temp;    


      /* Source AFTER_CODON is a special source for AFTER_CENTRAL */ 
      temp = GeneWise21_EXPL_SPECIAL(mat,0,j - 10,AFTER_CODON) + ((mat->gp21->central[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))) + (0);     
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      GeneWise21_EXPL_SPECIAL(mat,0,j,AFTER_CENTRAL) = score;    
      /* Finished updating state AFTER_CENTRAL */ 




      /* Special state AFTER_PY_TRACT has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GeneWise21_EXPL_SPECIAL(mat,0,j,AFTER_PY_TRACT);   


      /* Source AFTER_CENTRAL is a special source for AFTER_PY_TRACT */ 
      temp = GeneWise21_EXPL_SPECIAL(mat,0,j - 1,AFTER_CENTRAL) + (mat->gp21->transition[GP21_CENTRAL2PY]) + (mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)]);  
      if( temp > score ) 
        score = temp;    


      /* Source AFTER_PY_TRACT is a special source for AFTER_PY_TRACT */ 
      temp = GeneWise21_EXPL_SPECIAL(mat,0,j - 1,AFTER_PY_TRACT) + (mat->gp21->transition[GP21_PY2PY]) + (mat->gp21->py[CSEQ_GENOMIC_BASE(mat->target,j)]);  
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      GeneWise21_EXPL_SPECIAL(mat,0,j,AFTER_PY_TRACT) = score;   
      /* Finished updating state AFTER_PY_TRACT */ 




      /* Special state AFTER_SPACER has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GeneWise21_EXPL_SPECIAL(mat,0,j,AFTER_SPACER); 


      /* Source AFTER_PY_TRACT is a special source for AFTER_SPACER */ 
      temp = GeneWise21_EXPL_SPECIAL(mat,0,j - 1,AFTER_PY_TRACT) + (mat->gp21->transition[GP21_PY2SPACER]) + (mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)]);  
      if( temp > score ) 
        score = temp;    


      /* Source AFTER_SPACER is a special source for AFTER_SPACER */ 
      temp = GeneWise21_EXPL_SPECIAL(mat,0,j - 1,AFTER_SPACER) + (mat->gp21->transition[GP21_SPACER2SPACER]) + (mat->gp21->spacer[CSEQ_GENOMIC_BASE(mat->target,j)]);    
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      GeneWise21_EXPL_SPECIAL(mat,0,j,AFTER_SPACER) = score; 
      /* Finished updating state AFTER_SPACER */ 


      }  
    stop_reporting();    
    return TRUE;     
}    


/* Function:  GeneWise21_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GeneWise21 *]
 *
 */
GeneWise21 * GeneWise21_alloc(void) 
{
    GeneWise21 * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(GeneWise21 *) ckalloc (sizeof(GeneWise21))) == NULL)    {  
      warn("GeneWise21_alloc failed ");  
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


/* Function:  free_GeneWise21(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GeneWise21 *]
 *
 * Return [UNKN ]  Undocumented return value [GeneWise21 *]
 *
 */
GeneWise21 * free_GeneWise21(GeneWise21 * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a GeneWise21 obj. Should be trappable");    
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
    /* obj->gp21 is linked in */ 
    /* obj->rndcodon is linked in */ 
    /* obj->rndbase is linked in */ 


    ckfree(obj); 
    return NULL; 
}    





#ifdef _cplusplus
}
#endif
