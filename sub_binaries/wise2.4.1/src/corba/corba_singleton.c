
#include "corba_singleton.h"

#include <stdio.h>

static Wise2Corba_Singleton singleton;

static char * static_init_string = NULL;


Wise2Corba_Singleton * get_Wise2Corba_Singleton(int * argc,char ** argv,char * init_string)
{

  if( static_init_string != NULL ) {
    if( strcmp(static_init_string,init_string) == 0 ) {
      return &singleton;
    } else {
      fprintf(stderr,"Trouble! Asked for orbs of different types %s and %s",init_string,static_init_string);
      return NULL;
    }
  }

  /* else, initiate orb */

  singleton.ev = g_new0(CORBA_Environment,1);

  CORBA_exception_init(singleton.ev);
  singleton.orb = CORBA_ORB_init(argc,argv,init_string,singleton.ev);
  /* FIXME: Exception check! */
  singleton.poa = (PortableServer_POA)CORBA_ORB_resolve_initial_references(singleton.orb,
									   "RootPOA",singleton.ev);

  static_init_string = CORBA_string_dup(init_string);

  return &singleton;
}
