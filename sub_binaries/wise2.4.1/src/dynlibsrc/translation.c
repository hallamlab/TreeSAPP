#ifdef _cplusplus
extern "C" {
#endif
#include "translation.h"

/* Function:  copy_Translation(t)
 *
 * Descrip:    Makes a complete clean copy of the translation
 *
 *
 * Arg:        t [UNKN ] Undocumented argument [Translation *]
 *
 * Return [UNKN ]  Undocumented return value [Translation *]
 *
 */
# line 50 "translation.dy"
Translation * copy_Translation(Translation * t)
{
  Translation * out;
  out = Translation_alloc();
  out->start = t->start;
  out->end   = t->end;

  return out;

}

/* Function:  get_Protein_from_Translation(ts,ct)
 *
 * Descrip:    Gets the protein
 *
 *
 * Arg:        ts [UNKN ] translation [Translation *]
 * Arg:        ct [UNKN ] codon table to use [CodonTable *]
 *
 * Return [SOFT ]  Protein sequence [Protein *]
 *
 */
# line 68 "translation.dy"
Protein * get_Protein_from_Translation(Translation * ts,CodonTable * ct)
{
  cDNA * cd;
  int i,j;
  Sequence * seq;
  char buffer[64];

  assert(ts);
  assert(ct);

  /*  fprintf(stderr,"Codon table is %d\n",ct);*/

  if( ts->protein != NULL)
    return ts->protein;

  if( ts->parent == NULL ) {
    warn("Cannot get Protein from translation as no parent!");
    return NULL;
  }


  cd = get_cDNA_from_Transcript(ts->parent);

  if( cd == NULL ) {
    warn("Cannot make translation as can't get transcript!");
    return NULL;
  }

  if( cd->baseseq == NULL ) {
    warn("A bad error - a non NULL cDNA with a null sequence object. No translation here!");
    return NULL;
  }
  if( cd->baseseq->len == 0 ) {
    warn("Attempting to translate a zero length cDNA. Yikes!");
    return NULL;
  }

  seq = Sequence_alloc();
  sprintf(buffer,"%s.tr",cDNA_name(cd));
  seq->name = stringalloc(buffer);
  seq->seq = ckcalloc((cd->baseseq->len/3) + 2,sizeof(char));
  seq->type = SEQUENCE_PROTEIN;

  if( cd->baseseq->len%3 != 0 ) {
    warn("Problem in making translation, cDNA is not mod3! - length is %d - transcript id %s",cd->baseseq->len,seq->name);
  }


  for(i=0,j=0;i<cd->baseseq->len;i+=3,j++) {
    if( is_stop_codon(codon_from_seq(cd->baseseq->seq+i),ct) == TRUE ) {
      if( i+3 >= cd->baseseq->len ) 
	break;
      else {
	warn("Got a stop codon in the middle of a translation at postion [%d]. Yuk!",i);
	seq->seq[j] = '*';
      }
    } else {
      seq->seq[j] = aminoacid_from_seq(ct,cd->baseseq->seq+i);

    }
  }
  seq->seq[j]='\0';
  make_len_type_Sequence(seq);

  /*write_fasta_Sequence(seq,stdout);*/
  seq->type = SEQUENCE_PROTEIN;
  ts->protein = Protein_from_Sequence(seq);

  return ts->protein;

}


/* Function:  show_Translation(*ts,ofp)
 *
 * Descrip:    shows a translation in vaguely human form
 *
 *
 * Arg:        *ts [UNKN ] Undocumented argument [Translation]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 144 "translation.dy"
void show_Translation(Translation *ts,FILE * ofp)
{
  fprintf(ofp,"Translation %d - %d\n",ts->start,ts->end);
}


# line 122 "translation.c"
/* Function:  hard_link_Translation(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [Translation *]
 *
 * Return [UNKN ]  Undocumented return value [Translation *]
 *
 */
Translation * hard_link_Translation(Translation * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a Translation object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  Translation_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Translation *]
 *
 */
Translation * Translation_alloc(void) 
{
    Translation * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(Translation *) ckalloc (sizeof(Translation))) == NULL)  {  
      warn("Translation_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->start = 0;  
    out->end = 0;    
    out->protein = NULL; 


    return out;  
}    


/* Function:  free_Translation(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [Translation *]
 *
 * Return [UNKN ]  Undocumented return value [Translation *]
 *
 */
Translation * free_Translation(Translation * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a Translation obj. Should be trappable");   
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
    /* obj->parent is linked in */ 
    if( obj->protein != NULL)    
      free_Protein(obj->protein);    


    ckfree(obj); 
    return NULL; 
}    


/* Function:  replace_start_Translation(obj,start)
 *
 * Descrip:    Replace member variable start
 *             For use principly by API functions
 *
 *
 * Arg:          obj [UNKN ] Object holding the variable [Translation *]
 * Arg:        start [OWNER] New value of the variable [int]
 *
 * Return [SOFT ]  member variable start [boolean]
 *
 */
boolean replace_start_Translation(Translation * obj,int start) 
{
    if( obj == NULL)     {  
      warn("In replacement function start for object Translation, got a NULL object");   
      return FALSE;  
      }  
    obj->start = start;  
    return TRUE; 
}    


/* Function:  access_start_Translation(obj)
 *
 * Descrip:    Access member variable start
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [Translation *]
 *
 * Return [SOFT ]  member variable start [int]
 *
 */
int access_start_Translation(Translation * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function start for object Translation, got a NULL object");  
      return 0;  
      }  
    return obj->start;   
}    


/* Function:  replace_end_Translation(obj,end)
 *
 * Descrip:    Replace member variable end
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [Translation *]
 * Arg:        end [OWNER] New value of the variable [int]
 *
 * Return [SOFT ]  member variable end [boolean]
 *
 */
boolean replace_end_Translation(Translation * obj,int end) 
{
    if( obj == NULL)     {  
      warn("In replacement function end for object Translation, got a NULL object"); 
      return FALSE;  
      }  
    obj->end = end;  
    return TRUE; 
}    


/* Function:  access_end_Translation(obj)
 *
 * Descrip:    Access member variable end
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [Translation *]
 *
 * Return [SOFT ]  member variable end [int]
 *
 */
int access_end_Translation(Translation * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function end for object Translation, got a NULL object");    
      return 0;  
      }  
    return obj->end;     
}    


/* Function:  replace_parent_Translation(obj,parent)
 *
 * Descrip:    Replace member variable parent
 *             For use principly by API functions
 *
 *
 * Arg:           obj [UNKN ] Object holding the variable [Translation *]
 * Arg:        parent [OWNER] New value of the variable [Transcript *]
 *
 * Return [SOFT ]  member variable parent [boolean]
 *
 */
boolean replace_parent_Translation(Translation * obj,Transcript * parent) 
{
    if( obj == NULL)     {  
      warn("In replacement function parent for object Translation, got a NULL object");  
      return FALSE;  
      }  
    obj->parent = parent;    
    return TRUE; 
}    


/* Function:  access_parent_Translation(obj)
 *
 * Descrip:    Access member variable parent
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [Translation *]
 *
 * Return [SOFT ]  member variable parent [Transcript *]
 *
 */
Transcript * access_parent_Translation(Translation * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function parent for object Translation, got a NULL object"); 
      return NULL;   
      }  
    return obj->parent;  
}    


/* Function:  replace_protein_Translation(obj,protein)
 *
 * Descrip:    Replace member variable protein
 *             For use principly by API functions
 *
 *
 * Arg:            obj [UNKN ] Object holding the variable [Translation *]
 * Arg:        protein [OWNER] New value of the variable [Protein *]
 *
 * Return [SOFT ]  member variable protein [boolean]
 *
 */
boolean replace_protein_Translation(Translation * obj,Protein * protein) 
{
    if( obj == NULL)     {  
      warn("In replacement function protein for object Translation, got a NULL object"); 
      return FALSE;  
      }  
    obj->protein = protein;  
    return TRUE; 
}    


/* Function:  access_protein_Translation(obj)
 *
 * Descrip:    Access member variable protein
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [Translation *]
 *
 * Return [SOFT ]  member variable protein [Protein *]
 *
 */
Protein * access_protein_Translation(Translation * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function protein for object Translation, got a NULL object");    
      return NULL;   
      }  
    return obj->protein;     
}    



#ifdef _cplusplus
}
#endif
