#ifdef _cplusplus
extern "C" {
#endif
#include "hspscaninterface.h"

/* Function:  untyped_read_HSPScanInterfacePara_from_Stream(rs)
 *
 * Descrip:    untyped version of read
 *
 *
 * Arg:        rs [UNKN ] Undocumented argument [Wise2ReadStreamInterface *]
 *
 * Return [UNKN ]  Undocumented return value [void *]
 *
 */
# line 43 "hspscaninterface.dy"
void * untyped_read_HSPScanInterfacePara_from_Stream(Wise2ReadStreamInterface * rs)
{
  return (void*) typed_HSPScanInterfacePara_from_Stream(rs);
}

/* Function:  typed_HSPScanInterfacePara_from_Stream(rs)
 *
 * Descrip:    reads a HSPscan interface para from a stream
 *
 *
 * Arg:        rs [UNKN ] Undocumented argument [Wise2ReadStreamInterface *]
 *
 * Return [UNKN ]  Undocumented return value [HSPScanInterfacePara *]
 *
 */
# line 51 "hspscaninterface.dy"
HSPScanInterfacePara * typed_HSPScanInterfacePara_from_Stream(Wise2ReadStreamInterface * rs)
{
  char buffer[MAXLINE];
  HSPScanInterfacePara * out;

  out = HSPScanInterfacePara_std();


  while( WISE2_READ_BUFFER(buffer,MAXLINE,rs) != NULL ) {
    if( buffer[0] == '/' && buffer[1] == '/' ) {
      break;
    }

    if( strstartcmp(buffer,"numb") == 0 ) {
      out->numb_level = strtol(buffer+5,NULL,0);
    } else if ( strstartcmp(buffer,"maxresults") == 0 ) {
      out->max_results = strtol(buffer+11,NULL,0);
    } else if ( strstartcmp(buffer,"worddepth") == 0 ) {
      out->word_depth = strtol(buffer+10,NULL,0);
    } else if ( strstartcmp(buffer,"minword") == 0 ) {
      out->min_word_score = strtol(buffer+8,NULL,0);
    } else if ( strstartcmp(buffer,"minhsp") == 0 ) {
      out->min_hsp_score = strtol(buffer+7,NULL,0);
    } else if ( strstartcmp(buffer,"impl") == 0 ) {
      out->implementation = strtol(buffer+5,NULL,0);
    } else if ( strstartcmp(buffer,"linkwidth") == 0 ) {
      out->hsp_link_width = strtol(buffer+10,NULL,0);
    } else if ( strstartcmp(buffer,"linklength") == 0 ) {
      out->hsp_link_length = strtol(buffer+11,NULL,0);
    } else if ( strstartcmp(buffer,"impl") == 0 ) {
      out->implementation = strtol(buffer+5,NULL,0);
    } else if ( strstartcmp(buffer,"verbose") == 0 ) {
      out->verbosity = strtol(buffer+8,NULL,0);
    } else if ( strstartcmp(buffer,"lownumb") == 0 ) {
      out->low_numb = strtol(buffer+8,NULL,0);
    } else if ( strstartcmp(buffer,"havgext") == 0 ) {
      out->hsp_avg_ext = strtol(buffer+8,NULL,0);
    } else {
      warn("In reading HspScanInterfacePara, unknown line %s",buffer);
    }
  }


  return out;
}

/* Function:  untyped_HSPScanInterfacePara_to_Stream(p,ws)
 *
 * Descrip:    writes out a standard HSPscan interface
 *
 *
 * Arg:         p [UNKN ] Undocumented argument [void *]
 * Arg:        ws [UNKN ] Undocumented argument [Wise2WriteStreamInterface *]
 *
 */
# line 100 "hspscaninterface.dy"
void untyped_HSPScanInterfacePara_to_Stream(void * p,Wise2WriteStreamInterface * ws)
{
  HSPScanInterfacePara * h = (HSPScanInterfacePara *) p;
  char buffer[MAXLINE];

  sprintf(buffer,"numb %d\n",h->numb_level);
  WISE2_WRITE_STRING(buffer,ws);

  sprintf(buffer,"maxresults %d\n",h->max_results);
  WISE2_WRITE_STRING(buffer,ws);

  sprintf(buffer,"worddepth %d\n",h->word_depth);
  WISE2_WRITE_STRING(buffer,ws);

  sprintf(buffer,"minword %d\n",h->min_word_score);
  WISE2_WRITE_STRING(buffer,ws);

  sprintf(buffer,"minhsp %d\n",h->min_hsp_score);
  WISE2_WRITE_STRING(buffer,ws);

  sprintf(buffer,"linkwidth %d\n",h->hsp_link_width);
  WISE2_WRITE_STRING(buffer,ws);

  sprintf(buffer,"linklength %d\n",h->hsp_link_length);
  WISE2_WRITE_STRING(buffer,ws);

  sprintf(buffer,"impl %d\n",h->implementation);
  WISE2_WRITE_STRING(buffer,ws);

  sprintf(buffer,"verbose %d\n",h->verbosity);
  WISE2_WRITE_STRING(buffer,ws);

  sprintf(buffer,"lownumb %d\n",h->low_numb);
  WISE2_WRITE_STRING(buffer,ws);
  

  WISE2_WRITE_STRING("//\n",ws);

}

/* Function:  untyped_HSPScanInterfacePara_free(p)
 *
 * Descrip:    untyped free function for AnonymousObjects
 *
 *
 * Arg:        p [UNKN ] Undocumented argument [void *]
 *
 */
# line 143 "hspscaninterface.dy"
void untyped_HSPScanInterfacePara_free(void * p)
{
  HSPScanInterfacePara * h;

  h = (HSPScanInterfacePara *) p;
  
  free_HSPScanInterfacePara(h);

}


/* Function:  show_help_HSPScanInterfacePara(ofp)
 *
 * Descrip:    help function for hspscan interface para
 *
 *
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 157 "hspscaninterface.dy"
void show_help_HSPScanInterfacePara(FILE * ofp)
{
  fprintf(ofp,"Parameters for word scan\n");
  fprintf(ofp,"   -hspscan_impl [vanilla/threaded/twohit] implementation to use if compliled for runtime\n");
  fprintf(ofp,"   -hspscan_maxres    [100] maximum results returned by scan\n");
  fprintf(ofp,"   -hspscan_numb     [1000] word count to numb word (for low complexity)\n");
  fprintf(ofp,"   -hspscan_worddepth   [2] maximum offset from word - [0,1,2]\n");
  fprintf(ofp,"   -hspscan_minword    [14] minimum word score\n");
  fprintf(ofp,"   -hspscan_minhsp     [22] minimum hsp score\n");
  fprintf(ofp,"   -hspscan_link_width [30] max width (gap) of scored HSP chain\n");
  fprintf(ofp,"   -hspscan_link_length [150] max length of scored HSP chain\n");
  fprintf(ofp,"   -hspscan_verbosity  [0] verbosity level of server, if server is trace compiled\n");
  fprintf(ofp,"   -hspscan_lownumb    [0] low complexity numb level, 0 means no special low complexity scores\n");
  fprintf(ofp,"   -hspscan_avgext     [-6] average extension minimum in hsp extension\n");

}

/* Function:  new_HSPScanInterfacePara_from_argv(argc,argv)
 *
 * Descrip:    makes a hspscan interface from the command line
 *
 *
 * Arg:        argc [UNKN ] Undocumented argument [int *]
 * Arg:        argv [UNKN ] Undocumented argument [char **]
 *
 * Return [UNKN ]  Undocumented return value [HSPScanInterfacePara *]
 *
 */
# line 177 "hspscaninterface.dy"
HSPScanInterfacePara * new_HSPScanInterfacePara_from_argv(int * argc,char ** argv)
{
  char * temp;
  HSPScanInterfacePara * out;

  out = HSPScanInterfacePara_std();

  if( (temp = strip_out_assigned_argument(argc,argv,"hspscan_impl")) != NULL ) {
    if( strcmp(temp,"vanilla") == 0 ) {
      out->implementation = HSPSCAN_IMPLEMENTATION_VANILLA;
    } else if ( strcmp(temp,"threaded") == 0 ) {
      out->implementation = HSPSCAN_IMPLEMENTATION_THREADED;
    } else if ( strcmp(temp,"twohit") == 0 ) {
      out->implementation = HSPSCAN_IMPLEMENTATION_TWOHIT;
    } else {
      fatal("Cannot parse %s as an implementation, exiting early before this gets snarled up",temp);
    }
  }



  strip_out_integer_argument(argc,argv,"hspscan_maxres",&out->max_results);
  strip_out_integer_argument(argc,argv,"hspscan_numb",&out->numb_level);
  strip_out_integer_argument(argc,argv,"hspscan_worddepth",&out->word_depth);
  strip_out_integer_argument(argc,argv,"hspscan_minword",&out->min_word_score);
  strip_out_integer_argument(argc,argv,"hspscan_minhsp",&out->min_hsp_score);
  strip_out_integer_argument(argc,argv,"hspscan_link_width",&out->hsp_link_width);
  strip_out_integer_argument(argc,argv,"hspscan_link_length",&out->hsp_link_length);
  strip_out_integer_argument(argc,argv,"hspscan_verbosity",&out->verbosity);
  strip_out_integer_argument(argc,argv,"hspscan_lownumb",&out->low_numb);
  strip_out_integer_argument(argc,argv,"hspscan_avgext",&out->hsp_avg_ext);

  return out;
}


/* Function:  HSPScanInterfacePara_std(void)
 *
 * Descrip:    makes a standard hsp scan interface para
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [HSPScanInterfacePara *]
 *
 */
# line 216 "hspscaninterface.dy"
HSPScanInterfacePara * HSPScanInterfacePara_std(void)
{
  HSPScanInterfacePara * out;

  out = HSPScanInterfacePara_alloc();
  out->min_score = 20;
  out->max_results = 100;
  out->use_protein_heuristic = TRUE;
  out->numb_level = 1000;
  out->word_depth = 2;
  out->min_word_score = 14;
  out->min_hsp_score = 22;
  out->hsp_link_width = 30;
  out->hsp_link_length = 150;
  out->hsp_avg_ext = -6;
  out->implementation = HSPSCAN_IMPLEMENTATION_VANILLA;
  out->verbosity = 0;
  out->low_numb = 0;


  return out;
}


/* Function:  free_HSPScanInterface(hsi)
 *
 * Descrip:    Frees overrides dynamite default
 *
 *
 * Arg:        hsi [UNKN ] Undocumented argument [HSPScanInterface *]
 *
 * Return [UNKN ]  Undocumented return value [HSPScanInterface *]
 *
 */
# line 244 "hspscaninterface.dy"
HSPScanInterface * free_HSPScanInterface(HSPScanInterface * hsi)
{
  if( hsi == NULL ) {
    return NULL;
  }

  (*hsi->free_data)(hsi->data);
  ckfree(hsi);

  return NULL;

}




# line 267 "hspscaninterface.c"
/* Function:  hard_link_HSPScanInterfacePara(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [HSPScanInterfacePara *]
 *
 * Return [UNKN ]  Undocumented return value [HSPScanInterfacePara *]
 *
 */
HSPScanInterfacePara * hard_link_HSPScanInterfacePara(HSPScanInterfacePara * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a HSPScanInterfacePara object: passed a NULL object");    
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  HSPScanInterfacePara_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [HSPScanInterfacePara *]
 *
 */
HSPScanInterfacePara * HSPScanInterfacePara_alloc(void) 
{
    HSPScanInterfacePara * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(HSPScanInterfacePara *) ckalloc (sizeof(HSPScanInterfacePara))) == NULL)    {  
      warn("HSPScanInterfacePara_alloc failed ");    
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->min_score = 0;  
    out->max_results = 0;    
    out->flags = 0;  
    out->use_protein_heuristic = TRUE;   
    out->numb_level = 1000;  
    out->word_depth = 2; 
    out->min_word_score = 14;    
    out->min_hsp_score = 22; 
    out->implementation = HSPSCAN_IMPLEMENTATION_VANILLA;    
    out->hsp_link_width = 30;    
    out->hsp_link_length = 150;  
    out->verbosity = 0;  
    out->low_numb = 0;   
    out->hsp_avg_ext = -6;   


    return out;  
}    


/* Function:  free_HSPScanInterfacePara(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [HSPScanInterfacePara *]
 *
 * Return [UNKN ]  Undocumented return value [HSPScanInterfacePara *]
 *
 */
HSPScanInterfacePara * free_HSPScanInterfacePara(HSPScanInterfacePara * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a HSPScanInterfacePara obj. Should be trappable");  
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


/* Function:  hard_link_HSPScanInterface(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [HSPScanInterface *]
 *
 * Return [UNKN ]  Undocumented return value [HSPScanInterface *]
 *
 */
HSPScanInterface * hard_link_HSPScanInterface(HSPScanInterface * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a HSPScanInterface object: passed a NULL object");    
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  HSPScanInterface_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [HSPScanInterface *]
 *
 */
HSPScanInterface * HSPScanInterface_alloc(void) 
{
    HSPScanInterface * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(HSPScanInterface *) ckalloc (sizeof(HSPScanInterface))) == NULL)    {  
      warn("HSPScanInterface_alloc failed ");    
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->scan_query = NULL;  
    out->free_data = NULL;   


    return out;  
}    



#ifdef _cplusplus
}
#endif
