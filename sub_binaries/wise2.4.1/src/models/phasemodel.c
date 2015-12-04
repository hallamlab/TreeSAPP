#ifdef _cplusplus
extern "C" {
#endif
#include "phasemodel.h"

# line 59 "phasemodel.dy"
ProteinIntronList * read_ProteinIntronList(FILE * ifp)
{
  ProteinIntronList * out;
  ProteinIntron * in;
  char line[MAXLINE];
  char * run1;
  char * run2;
  int pos;
  int phase;

  out = ProteinIntronList_alloc_std();

  while( fgets(line,MAXLINE,ifp) != NULL ) {
    if( line[0] == '#' ) {
      continue;
    }
    if( !isnumber(line[0]) ) {
      warn("Bad looking line in intron file, %s",line);
      continue;
    }
    if( sscanf(line,"%d %d",&pos,&phase) != 2 ) {
      warn("Unable to parse line in intron file %s",line);
      continue;
    }
    if( phase > 2 || phase < 0 ) {
      warn("Got phase of %d, when must be between 0 and 2",phase);
      continue;
    }

    in = ProteinIntron_alloc();
    pos--;
    in->position = pos;
    in->phase = phase;
    add_ProteinIntronList(out,in);
  }

  return out;
}

# line 98 "phasemodel.dy"
ProteinIntronList * read_ProteinIntronList_from_filename(char * f)
{
  ProteinIntronList * pl;
  FILE * ifp;

  ifp = openfile(f,"r");
  if( ifp == NULL ) {
    warn("Cannot open file %s",f);
    return NULL;
  }

  pl = read_ProteinIntronList(ifp);

  fclose(ifp);

  return pl;
}


# line 117 "phasemodel.dy"
PhasedProteinPara * new_PhasedProteinPara_from_argv(int * argc,char ** argv)
{
  PhasedProteinPara * out;
  char * temp;

  out = PhasedProteinPara_alloc();
  out->marked_intron   = 0.95;
  out->unmarked_intron = 0.00001;
  out->use_phase = 0;

  strip_out_float_argument(argc,argv,"phase_marked",&out->marked_intron);
  strip_out_float_argument(argc,argv,"phase_unmarked",&out->unmarked_intron);
  /*  strip_out_boolean_def_argument(argc,argv,"phase_model",&out->use_phase); */
  if( (temp = strip_out_assigned_argument(argc,argv,"phase_file")) != NULL ) {
    out->intron_file = stringalloc(temp);
  }

  if( strip_out_boolean_argument(argc,argv,"phase_help") == TRUE ) {
    fprintf(stdout,"Phased marks provide the ability to restrict the position of introns\n");
    fprintf(stdout,"relative to the protein sequence; ie, assuming conserved introns. This\n");
    fprintf(stdout,"is most useful for fast evolving genes inside of relatively consistent\n");
    fprintf(stdout,"clades, eg for fast evolving genes, such as cytokines, in vertebrates\n");
    fprintf(stdout,"As moving between clades - say between Human and Drosophila - the intron\n");
    fprintf(stdout,"positions change, using these options would actively hinder good gene prediction\n");
    fprintf(stdout,"\n");
    fprintf(stdout,"This option can be used for either HMMs or proteins, although it is harder\n");
    fprintf(stdout,"to coordinate the HMM intron position than the protein positions.\n");
    fprintf(stdout,"Two things need to occur to use the phase information\n");
    fprintf(stdout,"    provide a phase mark file as -phase_file <xxxxxx>\n");
    fprintf(stdout,"    use the algorithm type 623P (6 states, 23 transitions, phased introns)\n");
    fprintf(stdout,"\n");
    fprintf(stdout,"The phase model attempts to make a ATG to STOP gene, even if the protein match\n");
    fprintf(stdout,"is not present across the entire gene. One major headache in this are introns in first\n");
    fprintf(stdout,"ATG, which is not handled at the moment\n\n");
    fprintf(stdout,"Genewise uses the protein position, in 1 coordinates, (first amino acid is 1)\n");
    fprintf(stdout,"for the definition of the intron. For phase 0 introns, it should be labeled as\n");
    fprintf(stdout,"the amino acid before the intron. For phase 1 and 2 introns, this is on the intron\n\n");
    fprintf(stdout,"We suggest using a small spread of positions to cope with intron positioning errors\n");
    fprintf(stdout,"  eg, defining an intron at position 4, phase 0, make postions 3,4 and 5 with position 0\n\n");
    fprintf(stdout,"The phase file format is\n");
    fprintf(stdout,"# lines starting with hash are comments\n");
    fprintf(stdout,"# three tab delimited columns\n");
    fprintf(stdout,"# <protein-position>  <phase>\n");
    fprintf(stdout,"# eg\n");
    fprintf(stdout,"4 0\n");
    exit(0);
  }

  return out;

}


# line 170 "phasemodel.dy"
void show_help_PhasedProteinPara(FILE * ofp)
{
  fprintf(ofp,"Phased Protein/HMM parameters (separate from other options)\n");
  fprintf(ofp,"  -phase_marked    [0.95]        Probability of marked introns\n");
  fprintf(ofp,"  -phase_unmarked  [0.00001]     Probability of unmarked introns\n");
  fprintf(ofp,"  -phase_file      [filename]    Intron positions and phases\n");
  fprintf(ofp,"  -phase_help      prints longer help on phase file and exits\n");
}


# line 180 "phasemodel.dy"
void write_fasta_PhasedProtein(PhasedProtein * pp,FILE * ofp)
{
  int i;
  int j;
  int line;
  assert(pp != NULL);
  
  write_fasta_Sequence(pp->protein,ofp);

  fprintf(ofp,">%s\n",pp->protein->name);

  for(i=0,j=0,line = 0;i<pp->protein->len;i++) {


    if( line != 0 && line % 50 == 0 ) {
      fputc('\n',ofp);
    }



    fputc(pp->protein->seq[i],ofp);
    line++;

    if( j < pp->list->len && pp->list->intron[j]->position == i ) {
      fputc('0'+pp->list->intron[j]->phase,ofp);
      line++;
      j++;
    }
  }
}


# line 212 "phasemodel.dy"
PhasedProtein * read_fasta_PhasedProtein_file(char * file)
{
  PhasedProtein * pp;
  FILE * ifp;

  ifp = openfile(file,"r");

  pp = read_fasta_PhasedProtein(ifp);

  fclose(ifp);

  return pp;

}

# line 227 "phasemodel.dy"
PhasedProtein * read_fasta_PhasedProtein(FILE * ifp)
{
  PhasedProtein * pp;
  ProteinIntron * intron;
  Sequence * input;
  int i;
  int j;
  char seqbuffer[10000];
  char name[2000];
  char c;

  pp = PhasedProtein_alloc();
  pp->list = ProteinIntronList_alloc_std();

  fgets(name,10000,ifp);
  assert(name[0] == '>');
  for(i=1; !isspace(name[i]);i++) {
    ;
  }
  name[i] = '\0';

  i = 0;

  while( (c = fgetc(ifp)) != EOF ) {
    if( c == '>' ) {
      ungetc('>',ifp);
      break;
    }
    
    if( isalpha(c) ) {
      seqbuffer[i++] = c;
    } else if( c == '0' || c == '1' || c == '2' ) {
      intron = ProteinIntron_alloc();
      intron->position = i-1;
      intron->phase    = c - '0';
      add_ProteinIntronList(pp->list,intron);
    } 
  }

  seqbuffer[i] = '\0';

  pp->protein = Sequence_from_static_memory(name+1,seqbuffer);

  return pp;
}



# line 275 "phasemodel.dy"
GenePhaseModel * GenePhaseModel_from_ThreeStateModel(ThreeStateModel * tsm,CodonMapper * cm,RandomModel * rm,CompMat * mat,PhasedProteinPara * ppp)
{
  GenePhaseModel * out;
  ProteinIntronList * pl;
  int i;
  int j;

  assert(ppp != NULL);
  assert(tsm != NULL);
  assert(rm != NULL);
  assert(mat != NULL);
  

  pl = read_ProteinIntronList_from_filename(ppp->intron_file);
  assert(pl != NULL);


  out = GenePhaseModel_alloc_std();


  /* current set this at global */

  for(i=0;i<tsm->len;i++) {
    tsm->unit[i]->transition[TSM_START2MATCH] = 0.0;
    tsm->unit[i]->transition[TSM_MATCH2END]   = 0.0;
  }

  tsm->unit[0]->transition[TSM_START2MATCH] = 1.0;
  tsm->unit[tsm->len-1]->transition[TSM_MATCH2END] = 1.0;

  
  out->gw = GeneWise_from_ThreeStateModel_setfactor(tsm,0.95,cm,0.1);


  for(i=0,j=0;i<tsm->len;i++) {
    /* set to default first */
    add_GenePhaseModel(out,GenePhaseSeg_alloc());
 
    out->phase[i]->intron_0 = ppp->unmarked_intron;
    out->phase[i]->intron_1 = ppp->unmarked_intron;
    out->phase[i]->intron_2 = ppp->unmarked_intron;
    out->phase[i]->insert_intron = ppp->unmarked_intron;

    if( j < pl->len && i == pl->intron[j]->position ) {
      if( pl->intron[j]->phase == 0 ) {
	out->phase[i]->intron_0 = ppp->marked_intron;
      } 
      if( pl->intron[j]->phase == 1 ) {
	out->phase[i]->intron_1 = ppp->marked_intron;
      } 
      if( pl->intron[j]->phase == 2 ) {
	out->phase[i]->intron_2 = ppp->marked_intron;
      } 

      out->gw->seg[i]->transition[GW_MATCH2MATCH]   = (1.0 - ppp->marked_intron)*out->gw->seg[i]->transition[GW_MATCH2MATCH];
      out->gw->seg[i]->transition[GW_MATCH2INSERT]  = (1.0 - ppp->marked_intron)*out->gw->seg[i]->transition[GW_MATCH2INSERT];
      out->gw->seg[i]->transition[GW_MATCH2DELETE]  = (1.0 - ppp->marked_intron)*out->gw->seg[i]->transition[GW_MATCH2DELETE];
      out->gw->seg[i]->transition[GW_DELETE2DELETE] = (1.0 - ppp->marked_intron)*out->gw->seg[i]->transition[GW_DELETE2DELETE];
      out->gw->seg[i]->transition[GW_DELETE2MATCH]  = (1.0 - ppp->marked_intron)*out->gw->seg[i]->transition[GW_DELETE2MATCH];

      j++;
    }
  }

  return out;
}




# line 345 "phasemodel.dy"
GenePhaseScore * GenePhaseScore_from_GenePhaseModel(GenePhaseModel * gpm)
{
  int i;
  GenePhaseScore * out;

  assert(gpm != NULL);
  assert(gpm->gw != NULL);
  assert(gpm->gw->len == gpm->len);

  out = GenePhaseScore_alloc_len(gpm->len);
  out->gws = GeneWiseScore_from_GeneWise(gpm->gw);

  for(i=0;i<gpm->len;i++) {
    add_GenePhaseScore(out,GenePhaseSegScore_alloc());
    out->phase[i]->intron_0  = Probability2Score(gpm->phase[i]->intron_0);
    out->phase[i]->intron_1  = Probability2Score(gpm->phase[i]->intron_1);
    out->phase[i]->intron_2  = Probability2Score(gpm->phase[i]->intron_2);
  }

  return out;

}



# line 325 "phasemodel.c"
/* Function:  hard_link_ProteinIntron(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [ProteinIntron *]
 *
 * Return [UNKN ]  Undocumented return value [ProteinIntron *]
 *
 */
ProteinIntron * hard_link_ProteinIntron(ProteinIntron * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a ProteinIntron object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  ProteinIntron_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ProteinIntron *]
 *
 */
ProteinIntron * ProteinIntron_alloc(void) 
{
    ProteinIntron * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(ProteinIntron *) ckalloc (sizeof(ProteinIntron))) == NULL)  {  
      warn("ProteinIntron_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->position = 0;   
    out->phase = 0;  


    return out;  
}    


/* Function:  free_ProteinIntron(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [ProteinIntron *]
 *
 * Return [UNKN ]  Undocumented return value [ProteinIntron *]
 *
 */
ProteinIntron * free_ProteinIntron(ProteinIntron * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a ProteinIntron obj. Should be trappable"); 
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


/* Function:  swap_ProteinIntronList(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_ProteinIntronList
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [ProteinIntron **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_ProteinIntronList(ProteinIntron ** list,int i,int j)  
{
    ProteinIntron * temp;    
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_ProteinIntronList(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_ProteinIntronList which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [ProteinIntron **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_ProteinIntronList(ProteinIntron ** list,int left,int right,int (*comp)(ProteinIntron * ,ProteinIntron * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_ProteinIntronList(list,left,(left+right)/2);    
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_ProteinIntronList (list,++last,i);  
      }  
    swap_ProteinIntronList (list,left,last); 
    qsort_ProteinIntronList(list,left,last-1,comp);  
    qsort_ProteinIntronList(list,last+1,right,comp); 
}    


/* Function:  sort_ProteinIntronList(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_ProteinIntronList
 *
 *
 * Arg:         obj [UNKN ] Object containing list [ProteinIntronList *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_ProteinIntronList(ProteinIntronList * obj,int (*comp)(ProteinIntron *, ProteinIntron *)) 
{
    qsort_ProteinIntronList(obj->intron,0,obj->len-1,comp);  
    return;  
}    


/* Function:  expand_ProteinIntronList(obj,len)
 *
 * Descrip:    Really an internal function for add_ProteinIntronList
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [ProteinIntronList *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_ProteinIntronList(ProteinIntronList * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_ProteinIntronList called with no need");  
      return TRUE;   
      }  


    if( (obj->intron = (ProteinIntron ** ) ckrealloc (obj->intron,sizeof(ProteinIntron *)*len)) == NULL)     {  
      warn("ckrealloc failed for expand_ProteinIntronList, returning FALSE");    
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_ProteinIntronList(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [ProteinIntronList *]
 * Arg:        add [OWNER] Object to add to the list [ProteinIntron *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_ProteinIntronList(ProteinIntronList * obj,ProteinIntron * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_ProteinIntronList(obj,obj->len + ProteinIntronListLISTLENGTH) == FALSE) 
        return FALSE;    
      }  


    obj->intron[obj->len++]=add; 
    return TRUE; 
}    


/* Function:  flush_ProteinIntronList(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [ProteinIntronList *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_ProteinIntronList(ProteinIntronList * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->intron[i] != NULL)    {  
        free_ProteinIntron(obj->intron[i]);  
        obj->intron[i] = NULL;   
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  ProteinIntronList_alloc_std(void)
 *
 * Descrip:    Equivalent to ProteinIntronList_alloc_len(ProteinIntronListLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ProteinIntronList *]
 *
 */
ProteinIntronList * ProteinIntronList_alloc_std(void) 
{
    return ProteinIntronList_alloc_len(ProteinIntronListLISTLENGTH); 
}    


/* Function:  ProteinIntronList_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [ProteinIntronList *]
 *
 */
ProteinIntronList * ProteinIntronList_alloc_len(int len) 
{
    ProteinIntronList * out;/* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = ProteinIntronList_alloc()) == NULL)    
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->intron = (ProteinIntron ** ) ckcalloc (len,sizeof(ProteinIntron *))) == NULL)   {  
      warn("Warning, ckcalloc failed in ProteinIntronList_alloc_len");   
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_ProteinIntronList(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [ProteinIntronList *]
 *
 * Return [UNKN ]  Undocumented return value [ProteinIntronList *]
 *
 */
ProteinIntronList * hard_link_ProteinIntronList(ProteinIntronList * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a ProteinIntronList object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  ProteinIntronList_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ProteinIntronList *]
 *
 */
ProteinIntronList * ProteinIntronList_alloc(void) 
{
    ProteinIntronList * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(ProteinIntronList *) ckalloc (sizeof(ProteinIntronList))) == NULL)  {  
      warn("ProteinIntronList_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->intron = NULL;  
    out->len = out->maxlen = 0;  


    return out;  
}    


/* Function:  free_ProteinIntronList(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [ProteinIntronList *]
 *
 * Return [UNKN ]  Undocumented return value [ProteinIntronList *]
 *
 */
ProteinIntronList * free_ProteinIntronList(ProteinIntronList * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a ProteinIntronList obj. Should be trappable"); 
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
    if( obj->intron != NULL) {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->intron[i] != NULL)  
          free_ProteinIntron(obj->intron[i]);    
        }  
      ckfree(obj->intron);   
      }  


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_PhasedProtein(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [PhasedProtein *]
 *
 * Return [UNKN ]  Undocumented return value [PhasedProtein *]
 *
 */
PhasedProtein * hard_link_PhasedProtein(PhasedProtein * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a PhasedProtein object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  PhasedProtein_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [PhasedProtein *]
 *
 */
PhasedProtein * PhasedProtein_alloc(void) 
{
    PhasedProtein * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(PhasedProtein *) ckalloc (sizeof(PhasedProtein))) == NULL)  {  
      warn("PhasedProtein_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->protein = NULL; 
    out->list = NULL;    


    return out;  
}    


/* Function:  free_PhasedProtein(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [PhasedProtein *]
 *
 * Return [UNKN ]  Undocumented return value [PhasedProtein *]
 *
 */
PhasedProtein * free_PhasedProtein(PhasedProtein * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a PhasedProtein obj. Should be trappable"); 
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
    if( obj->protein != NULL)    
      free_Sequence(obj->protein);   
    if( obj->list != NULL)   
      free_ProteinIntronList(obj->list);     


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_PhasedProteinPara(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [PhasedProteinPara *]
 *
 * Return [UNKN ]  Undocumented return value [PhasedProteinPara *]
 *
 */
PhasedProteinPara * hard_link_PhasedProteinPara(PhasedProteinPara * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a PhasedProteinPara object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  PhasedProteinPara_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [PhasedProteinPara *]
 *
 */
PhasedProteinPara * PhasedProteinPara_alloc(void) 
{
    PhasedProteinPara * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(PhasedProteinPara *) ckalloc (sizeof(PhasedProteinPara))) == NULL)  {  
      warn("PhasedProteinPara_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->marked_intron = 0.95;   
    out->unmarked_intron = 0.000001; 
    out->gap = -12;  
    out->ext = -2;   
    out->use_phase = FALSE;  
    out->intron_file = NULL; 


    return out;  
}    


/* Function:  free_PhasedProteinPara(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [PhasedProteinPara *]
 *
 * Return [UNKN ]  Undocumented return value [PhasedProteinPara *]
 *
 */
PhasedProteinPara * free_PhasedProteinPara(PhasedProteinPara * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a PhasedProteinPara obj. Should be trappable"); 
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
    if( obj->intron_file != NULL)    
      ckfree(obj->intron_file);  


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_PhasedHMM(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [PhasedHMM *]
 *
 * Return [UNKN ]  Undocumented return value [PhasedHMM *]
 *
 */
PhasedHMM * hard_link_PhasedHMM(PhasedHMM * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a PhasedHMM object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  PhasedHMM_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [PhasedHMM *]
 *
 */
PhasedHMM * PhasedHMM_alloc(void) 
{
    PhasedHMM * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(PhasedHMM *) ckalloc (sizeof(PhasedHMM))) == NULL)  {  
      warn("PhasedHMM_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->tsm = NULL; 
    out->list = NULL;    


    return out;  
}    


/* Function:  free_PhasedHMM(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [PhasedHMM *]
 *
 * Return [UNKN ]  Undocumented return value [PhasedHMM *]
 *
 */
PhasedHMM * free_PhasedHMM(PhasedHMM * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a PhasedHMM obj. Should be trappable"); 
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
    if( obj->tsm != NULL)    
      free_ThreeStateModel(obj->tsm);    
    if( obj->list != NULL)   
      free_ProteinIntronList(obj->list);     


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_GenePhaseSeg(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GenePhaseSeg *]
 *
 * Return [UNKN ]  Undocumented return value [GenePhaseSeg *]
 *
 */
GenePhaseSeg * hard_link_GenePhaseSeg(GenePhaseSeg * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a GenePhaseSeg object: passed a NULL object");    
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  GenePhaseSeg_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GenePhaseSeg *]
 *
 */
GenePhaseSeg * GenePhaseSeg_alloc(void) 
{
    GenePhaseSeg * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(GenePhaseSeg *) ckalloc (sizeof(GenePhaseSeg))) == NULL)    {  
      warn("GenePhaseSeg_alloc failed ");    
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->intron_0 = 0.0; 
    out->intron_1 = 0.0; 
    out->intron_2 = 0.0; 
    out->insert_intron = 0.0;    


    return out;  
}    


/* Function:  free_GenePhaseSeg(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GenePhaseSeg *]
 *
 * Return [UNKN ]  Undocumented return value [GenePhaseSeg *]
 *
 */
GenePhaseSeg * free_GenePhaseSeg(GenePhaseSeg * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a GenePhaseSeg obj. Should be trappable");  
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


/* Function:  swap_GenePhaseModel(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_GenePhaseModel
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [GenePhaseSeg **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_GenePhaseModel(GenePhaseSeg ** list,int i,int j)  
{
    GenePhaseSeg * temp; 
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_GenePhaseModel(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_GenePhaseModel which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [GenePhaseSeg **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_GenePhaseModel(GenePhaseSeg ** list,int left,int right,int (*comp)(GenePhaseSeg * ,GenePhaseSeg * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_GenePhaseModel(list,left,(left+right)/2);   
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_GenePhaseModel (list,++last,i); 
      }  
    swap_GenePhaseModel (list,left,last);    
    qsort_GenePhaseModel(list,left,last-1,comp); 
    qsort_GenePhaseModel(list,last+1,right,comp);    
}    


/* Function:  sort_GenePhaseModel(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_GenePhaseModel
 *
 *
 * Arg:         obj [UNKN ] Object containing list [GenePhaseModel *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_GenePhaseModel(GenePhaseModel * obj,int (*comp)(GenePhaseSeg *, GenePhaseSeg *)) 
{
    qsort_GenePhaseModel(obj->phase,0,obj->len-1,comp);  
    return;  
}    


/* Function:  expand_GenePhaseModel(obj,len)
 *
 * Descrip:    Really an internal function for add_GenePhaseModel
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [GenePhaseModel *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_GenePhaseModel(GenePhaseModel * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_GenePhaseModel called with no need"); 
      return TRUE;   
      }  


    if( (obj->phase = (GenePhaseSeg ** ) ckrealloc (obj->phase,sizeof(GenePhaseSeg *)*len)) == NULL)     {  
      warn("ckrealloc failed for expand_GenePhaseModel, returning FALSE");   
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_GenePhaseModel(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [GenePhaseModel *]
 * Arg:        add [OWNER] Object to add to the list [GenePhaseSeg *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_GenePhaseModel(GenePhaseModel * obj,GenePhaseSeg * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_GenePhaseModel(obj,obj->len + GenePhaseModelLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->phase[obj->len++]=add;  
    return TRUE; 
}    


/* Function:  flush_GenePhaseModel(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [GenePhaseModel *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_GenePhaseModel(GenePhaseModel * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->phase[i] != NULL) {  
        free_GenePhaseSeg(obj->phase[i]);    
        obj->phase[i] = NULL;    
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  GenePhaseModel_alloc_std(void)
 *
 * Descrip:    Equivalent to GenePhaseModel_alloc_len(GenePhaseModelLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GenePhaseModel *]
 *
 */
GenePhaseModel * GenePhaseModel_alloc_std(void) 
{
    return GenePhaseModel_alloc_len(GenePhaseModelLISTLENGTH);   
}    


/* Function:  GenePhaseModel_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [GenePhaseModel *]
 *
 */
GenePhaseModel * GenePhaseModel_alloc_len(int len) 
{
    GenePhaseModel * out;   /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = GenePhaseModel_alloc()) == NULL)   
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->phase = (GenePhaseSeg ** ) ckcalloc (len,sizeof(GenePhaseSeg *))) == NULL)  {  
      warn("Warning, ckcalloc failed in GenePhaseModel_alloc_len");  
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_GenePhaseModel(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GenePhaseModel *]
 *
 * Return [UNKN ]  Undocumented return value [GenePhaseModel *]
 *
 */
GenePhaseModel * hard_link_GenePhaseModel(GenePhaseModel * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a GenePhaseModel object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  GenePhaseModel_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GenePhaseModel *]
 *
 */
GenePhaseModel * GenePhaseModel_alloc(void) 
{
    GenePhaseModel * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(GenePhaseModel *) ckalloc (sizeof(GenePhaseModel))) == NULL)    {  
      warn("GenePhaseModel_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->phase = NULL;   
    out->len = out->maxlen = 0;  
    out->gw = NULL;  


    return out;  
}    


/* Function:  free_GenePhaseModel(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GenePhaseModel *]
 *
 * Return [UNKN ]  Undocumented return value [GenePhaseModel *]
 *
 */
GenePhaseModel * free_GenePhaseModel(GenePhaseModel * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a GenePhaseModel obj. Should be trappable");    
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
    if( obj->phase != NULL)  {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->phase[i] != NULL)   
          free_GenePhaseSeg(obj->phase[i]);  
        }  
      ckfree(obj->phase);    
      }  
    if( obj->gw != NULL) 
      free_GeneWise(obj->gw);    


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_GenePhaseSegScore(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GenePhaseSegScore *]
 *
 * Return [UNKN ]  Undocumented return value [GenePhaseSegScore *]
 *
 */
GenePhaseSegScore * hard_link_GenePhaseSegScore(GenePhaseSegScore * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a GenePhaseSegScore object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  GenePhaseSegScore_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GenePhaseSegScore *]
 *
 */
GenePhaseSegScore * GenePhaseSegScore_alloc(void) 
{
    GenePhaseSegScore * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(GenePhaseSegScore *) ckalloc (sizeof(GenePhaseSegScore))) == NULL)  {  
      warn("GenePhaseSegScore_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->intron_0 = 0;   
    out->intron_1 = 0;   
    out->intron_2 = 0;   
    out->insert_intron = 0;  


    return out;  
}    


/* Function:  free_GenePhaseSegScore(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GenePhaseSegScore *]
 *
 * Return [UNKN ]  Undocumented return value [GenePhaseSegScore *]
 *
 */
GenePhaseSegScore * free_GenePhaseSegScore(GenePhaseSegScore * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a GenePhaseSegScore obj. Should be trappable"); 
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


/* Function:  swap_GenePhaseScore(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_GenePhaseScore
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [GenePhaseSegScore **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_GenePhaseScore(GenePhaseSegScore ** list,int i,int j)  
{
    GenePhaseSegScore * temp;    
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_GenePhaseScore(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_GenePhaseScore which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [GenePhaseSegScore **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_GenePhaseScore(GenePhaseSegScore ** list,int left,int right,int (*comp)(GenePhaseSegScore * ,GenePhaseSegScore * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_GenePhaseScore(list,left,(left+right)/2);   
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_GenePhaseScore (list,++last,i); 
      }  
    swap_GenePhaseScore (list,left,last);    
    qsort_GenePhaseScore(list,left,last-1,comp); 
    qsort_GenePhaseScore(list,last+1,right,comp);    
}    


/* Function:  sort_GenePhaseScore(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_GenePhaseScore
 *
 *
 * Arg:         obj [UNKN ] Object containing list [GenePhaseScore *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_GenePhaseScore(GenePhaseScore * obj,int (*comp)(GenePhaseSegScore *, GenePhaseSegScore *)) 
{
    qsort_GenePhaseScore(obj->phase,0,obj->len-1,comp);  
    return;  
}    


/* Function:  expand_GenePhaseScore(obj,len)
 *
 * Descrip:    Really an internal function for add_GenePhaseScore
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [GenePhaseScore *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_GenePhaseScore(GenePhaseScore * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_GenePhaseScore called with no need"); 
      return TRUE;   
      }  


    if( (obj->phase = (GenePhaseSegScore ** ) ckrealloc (obj->phase,sizeof(GenePhaseSegScore *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_GenePhaseScore, returning FALSE");   
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_GenePhaseScore(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [GenePhaseScore *]
 * Arg:        add [OWNER] Object to add to the list [GenePhaseSegScore *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_GenePhaseScore(GenePhaseScore * obj,GenePhaseSegScore * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_GenePhaseScore(obj,obj->len + GenePhaseScoreLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->phase[obj->len++]=add;  
    return TRUE; 
}    


/* Function:  flush_GenePhaseScore(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [GenePhaseScore *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_GenePhaseScore(GenePhaseScore * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->phase[i] != NULL) {  
        free_GenePhaseSegScore(obj->phase[i]);   
        obj->phase[i] = NULL;    
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  GenePhaseScore_alloc_std(void)
 *
 * Descrip:    Equivalent to GenePhaseScore_alloc_len(GenePhaseScoreLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GenePhaseScore *]
 *
 */
GenePhaseScore * GenePhaseScore_alloc_std(void) 
{
    return GenePhaseScore_alloc_len(GenePhaseScoreLISTLENGTH);   
}    


/* Function:  GenePhaseScore_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [GenePhaseScore *]
 *
 */
GenePhaseScore * GenePhaseScore_alloc_len(int len) 
{
    GenePhaseScore * out;   /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = GenePhaseScore_alloc()) == NULL)   
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->phase = (GenePhaseSegScore ** ) ckcalloc (len,sizeof(GenePhaseSegScore *))) == NULL)    {  
      warn("Warning, ckcalloc failed in GenePhaseScore_alloc_len");  
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_GenePhaseScore(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GenePhaseScore *]
 *
 * Return [UNKN ]  Undocumented return value [GenePhaseScore *]
 *
 */
GenePhaseScore * hard_link_GenePhaseScore(GenePhaseScore * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a GenePhaseScore object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  GenePhaseScore_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GenePhaseScore *]
 *
 */
GenePhaseScore * GenePhaseScore_alloc(void) 
{
    GenePhaseScore * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(GenePhaseScore *) ckalloc (sizeof(GenePhaseScore))) == NULL)    {  
      warn("GenePhaseScore_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->phase = NULL;   
    out->len = out->maxlen = 0;  
    out->gws = NULL; 


    return out;  
}    


/* Function:  free_GenePhaseScore(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GenePhaseScore *]
 *
 * Return [UNKN ]  Undocumented return value [GenePhaseScore *]
 *
 */
GenePhaseScore * free_GenePhaseScore(GenePhaseScore * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a GenePhaseScore obj. Should be trappable");    
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
    if( obj->phase != NULL)  {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->phase[i] != NULL)   
          free_GenePhaseSegScore(obj->phase[i]); 
        }  
      ckfree(obj->phase);    
      }  
    if( obj->gws != NULL)    
      free_GeneWiseScore(obj->gws);  


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif
