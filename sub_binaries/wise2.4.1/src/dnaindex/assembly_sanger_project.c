#ifdef _cplusplus
extern "C" {
#endif
#include "assembly_sanger_project.h"


# line 23 "assembly_sanger_project.dy"
AssemblySequence * next_AssemblySequence_sanger_impl(void * h)
{
  char buffer[MAXLINE];
  AssemblySequence * aseq;
  char * run;

  FILE * ifp;
  SangerProjectDirectory * s = (SangerProjectDirectory *) h;
  struct dirent * dirp;
  
  while( (dirp = readdir(s->dir)) != NULL ) {
    if( (run = strstr(dirp->d_name,s->extension)) != NULL ) {
      /* should check it is at the end */
      if( strcmp(run,s->extension) != 0 ) {
	/* not at end */
	continue;
      }

      sprintf(buffer,"%s/%s",s->directory,dirp->d_name);
      ifp = fopen(buffer,"r");
      if( ifp == NULL ) {
	warn("This is crazy! name is in directory but cannot open file. Yikes");
	return NULL;
      }
      
      aseq = read_clipped_sanger_AssemblySequence(ifp);
      fclose(ifp);

      return aseq;
    }
  }


  return NULL;
}



# line 61 "assembly_sanger_project.dy"
AssemblySequence * read_clipped_sanger_AssemblySequence(FILE * ifp)
{
  AssemblySequence * out;
  Sequence * in;
  char * seqstr;
  char * run;
  int maxlen;
  int curr_pos;
  int in_seq = 0;
  int clip_start = -1;
  int clip_end = -1;

  int temp;

  char buffer[512];
  char name[60];

  seqstr = malloc(sizeof(char)*1028);
  maxlen = 1028;
  curr_pos = 0;

  while( fgets(buffer,512,ifp) != NULL ) {
    if( strstartcmp(buffer,"ID") == 0 ) {
      sscanf(buffer,"ID %s",name);
      continue;
    }
    if( strstartcmp(buffer,"QL") == 0 ) {
      sscanf(buffer,"QL %d",&clip_start);
      continue;
    } 
    if( strstartcmp(buffer,"SL") == 0 ) {
      sscanf(buffer,"SL %d",&temp);
      if( temp > clip_start ) {
	clip_start = temp;
      }
      continue;
    } 
    if( strstartcmp(buffer,"QR") == 0 ) {
      sscanf(buffer,"QR %d",&clip_end);
      continue;
    } 

    if( strstartcmp(buffer,"SQ") == 0 ) {
      in_seq = 1;
      continue;
    }

    if( strstartcmp(buffer,"//") == 0 ) {
      in_seq = 0;
      continue;
    }

    if( in_seq == 1 ) {
      for(run=buffer;*run;run++) {
	if( !isspace(*run) ) {
	  if( *run == '-' ) {
	    seqstr[curr_pos++] = 'N';
	  } else {
	    seqstr[curr_pos++] = *run;
	  }
	}
	if( curr_pos >= maxlen ) {
	  maxlen *= 2;
	  seqstr = realloc(seqstr,sizeof(char)*(maxlen));
	}
	 
      }
    }
  }
  seqstr[curr_pos] = '\0';
   
  if( clip_start == -1 ) {
    warn("No start clipping point found");
    clip_start = 1;
  } 
  if( clip_end == -1 ) {
    warn("No end clipping point found");
    clip_end = strlen(seqstr);
  }
  /* seq to C coords*/
  clip_start--;

  seqstr[clip_end] = '\0';

  in = Sequence_from_static_memory(name,seqstr+clip_start);
  in->type = SEQUENCE_DNA;

  out = AssemblySequence_alloc_std();
  out->seq = in;

  return out;

}


# line 156 "assembly_sanger_project.dy"
void free_handle_sanger_impl(void * h)
{
  SangerProjectDirectory * s = (SangerProjectDirectory *) h;

  free_SangerProjectDirectory(s);
}

# line 163 "assembly_sanger_project.dy"
AssemblySequenceStream * new_sanger_project_AssemblySequenceStream(char * dir_name,char * extension)
{
  AssemblySequenceStream * out;
  SangerProjectDirectory * spd;

  spd = new_SangerProjectDirectory(dir_name,extension);

  out = AssemblySequenceStream_alloc();
  out->handle = (void*) spd;

  out->next_AssemblySequence = next_AssemblySequence_sanger_impl;
  out->free_handle = free_handle_sanger_impl;

  return out;
}

# line 179 "assembly_sanger_project.dy"
SangerProjectDirectory * new_SangerProjectDirectory(char * dir_name,char * extension)
{
  DIR * d;
  SangerProjectDirectory * out;

  d = opendir(dir_name);
  if( d == NULL ) {
    warn("Could not open directory %s",dir_name);
    return NULL;
  }

  out = SangerProjectDirectory_alloc();

  out->dir = d;
  out->extension = stringalloc(extension);
  out->directory = stringalloc(dir_name);
  return out;

}

/* Function:  free_SangerProjectDirectory(obj)
 *
 * Descrip:    deconstructor
 *
 *
 * Arg:        obj [UNKN ] Undocumented argument [SangerProjectDirectory *]
 *
 * Return [UNKN ]  Undocumented return value [SangerProjectDirectory *]
 *
 */
# line 203 "assembly_sanger_project.dy"
SangerProjectDirectory * free_SangerProjectDirectory(SangerProjectDirectory * obj)
{
  int return_early = 0;    
  
  
  if( obj == NULL) {  
    warn("Attempting to free a NULL pointer to a Sanger Project obj. Should be trappable");    
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


  closedir(obj->dir);

  if( obj->extension != NULL ) {
    free(obj->extension);
  }

  free(obj);
  
  return NULL;
}


# line 232 "assembly_sanger_project.c"
/* Function:  hard_link_SangerProjectDirectory(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [SangerProjectDirectory *]
 *
 * Return [UNKN ]  Undocumented return value [SangerProjectDirectory *]
 *
 */
SangerProjectDirectory * hard_link_SangerProjectDirectory(SangerProjectDirectory * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a SangerProjectDirectory object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  SangerProjectDirectory_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SangerProjectDirectory *]
 *
 */
SangerProjectDirectory * SangerProjectDirectory_alloc(void) 
{
    SangerProjectDirectory * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(SangerProjectDirectory *) ckalloc (sizeof(SangerProjectDirectory))) == NULL)    {  
      warn("SangerProjectDirectory_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->dir = NULL; 
    out->extension = NULL;   
    out->directory = NULL;   


    return out;  
}    



#ifdef _cplusplus
}
#endif
