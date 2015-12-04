
#ifndef CORBA_SINGLETON_HEADER
#define CORBA_SINGLETON_HEADER

#include <orb/orbit.h>


typedef struct Wise2Corba_Singleton {
  CORBA_Environment * ev;
  CORBA_ORB orb;
  PortableServer_POA poa;
} Wise2Corba_Singleton;


Wise2Corba_Singleton * get_Wise2Corba_Singleton(int * argc,char ** argv,char * init_string);

#endif
