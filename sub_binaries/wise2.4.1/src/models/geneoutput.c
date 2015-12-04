#ifdef _cplusplus
extern "C" {
#endif
#include "geneoutput.h"


# line 35 "geneoutput.dy"
void show_GeneOutput(GeneOutputData * data,GeneOutputPara * para,FILE * ofp)
{
  Protein * trans;
  cDNA * cdna;
  int i;

  assert(data);
  assert(data->ct);
  assert(data->gr);
  assert(data->gen);


  if( para->show_genes ) {
    show_pretty_GenomicRegion(data->gr,0,ofp);
    fprintf(stdout,"%s\n",para->divide_string);
  }

  if( para->show_gff ) {
    show_GFF_GenomicRegion(data->gr,data->gen->baseseq->name,"genomwise",stdout);
    fprintf(stdout,"%s\n",para->divide_string);
  }

  if( para->show_trans ) {
    for(i=0;i<data->gr->len;i++) {
      if( data->gr->gene[i]->ispseudo == TRUE ) {
	fprintf(stdout,"#Gene %d is a pseudo gene - no translation possible\n",i);
      } else {
	trans = get_Protein_from_Translation(data->gr->gene[i]->transcript[0]->translation[0],data->ct);
	write_fasta_Sequence(trans->baseseq,ofp);
      }
    } 
    fprintf(stdout,"%s\n",para->divide_string);
  }

  if( para->show_cdna ) {
    for(i=0;i<data->gr->len;i++) {
      cdna = get_cDNA_from_Transcript(data->gr->gene[i]->transcript[0]);
      write_fasta_Sequence(cdna->baseseq,ofp);
    } 
    fprintf(stdout,"%s\n",para->divide_string);
  }

  if( para->show_geneutr ) {
    show_utr_exon_genomewise(data->alb,ofp);
    fprintf(stdout,"%s\n",para->divide_string);
  }



}

# line 86 "geneoutput.dy"
double id_map_func(int i)
{
  return (double)i;
}

# line 91 "geneoutput.dy"
GeneOutputPara * new_GeneOutputPara_from_argv(int * argc,char ** argv)
{
  GeneOutputPara * out;

  out = GeneOutputPara_alloc();
  out->show_genes = 1;
  out->show_trans = 1;
  out->divide_string = stringalloc("//");

  strip_out_boolean_def_argument(argc,argv,"geneutr",&out->show_geneutr);
  strip_out_boolean_def_argument(argc,argv,"genes",&out->show_genes);
  strip_out_boolean_def_argument(argc,argv,"trans",&out->show_trans);
  strip_out_boolean_def_argument(argc,argv,"gff",&out->show_gff);
  strip_out_boolean_def_argument(argc,argv,"cdna",&out->show_cdna);
  strip_out_boolean_def_argument(argc,argv,"genedebug",&out->show_debug);

  return out;
}


# line 111 "geneoutput.dy"
void show_help_GeneOutputPara(FILE * ofp)
{

  fprintf(ofp,"Gene Output\n");
  fprintf(ofp,"   -[no]genes       show gene structure (default yes)\n");
  fprintf(ofp,"   -[no]geneutr     show gene structure with utrs (default no)\n");
  fprintf(ofp,"   -[no]trans       show protein translation (default yes)\n");
  fprintf(ofp,"   -[no]gff         show gff (default no)\n");
  fprintf(ofp,"   -[no]genedebug   show gene debug\n");

}





# line 102 "geneoutput.c"
/* Function:  hard_link_GeneOutputData(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GeneOutputData *]
 *
 * Return [UNKN ]  Undocumented return value [GeneOutputData *]
 *
 */
GeneOutputData * hard_link_GeneOutputData(GeneOutputData * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a GeneOutputData object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  GeneOutputData_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GeneOutputData *]
 *
 */
GeneOutputData * GeneOutputData_alloc(void) 
{
    GeneOutputData * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(GeneOutputData *) ckalloc (sizeof(GeneOutputData))) == NULL)    {  
      warn("GeneOutputData_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->alb = NULL; 
    out->pal = NULL; 
    out->gr = NULL;  
    out->gen = NULL; 
    out->ct = NULL;  


    return out;  
}    


/* Function:  free_GeneOutputData(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GeneOutputData *]
 *
 * Return [UNKN ]  Undocumented return value [GeneOutputData *]
 *
 */
GeneOutputData * free_GeneOutputData(GeneOutputData * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a GeneOutputData obj. Should be trappable");    
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
    if( obj->alb != NULL)    
      free_AlnBlock(obj->alb);   
    if( obj->pal != NULL)    
      free_PackAln(obj->pal);    
    if( obj->gr != NULL) 
      free_GenomicRegion(obj->gr);   
    if( obj->gen != NULL)    
      free_Genomic(obj->gen);    
    if( obj->ct != NULL) 
      free_CodonTable(obj->ct);  


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_GeneOutputPara(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GeneOutputPara *]
 *
 * Return [UNKN ]  Undocumented return value [GeneOutputPara *]
 *
 */
GeneOutputPara * hard_link_GeneOutputPara(GeneOutputPara * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a GeneOutputPara object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  GeneOutputPara_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GeneOutputPara *]
 *
 */
GeneOutputPara * GeneOutputPara_alloc(void) 
{
    GeneOutputPara * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(GeneOutputPara *) ckalloc (sizeof(GeneOutputPara))) == NULL)    {  
      warn("GeneOutputPara_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->show_genes = 0; 
    out->show_gff = 0;   
    out->show_trans = 0; 
    out->show_cdna = 0;  
    out->show_geneutr = 0;   
    out->show_alb = 0;   
    out->show_pal = 0;   
    out->show_debug = 0; 
    out->divide_string = NULL;   


    return out;  
}    


/* Function:  free_GeneOutputPara(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GeneOutputPara *]
 *
 * Return [UNKN ]  Undocumented return value [GeneOutputPara *]
 *
 */
GeneOutputPara * free_GeneOutputPara(GeneOutputPara * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a GeneOutputPara obj. Should be trappable");    
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
    if( obj->divide_string != NULL)  
      ckfree(obj->divide_string);    


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif
