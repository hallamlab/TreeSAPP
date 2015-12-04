

#ifndef HSPSCAN_SERVER_IMPL_HEADER
#define HSPSCAN_SERVER_IMPL_HEADER

#include "hspscaninterface.h"
#include "hspscan_corba.h"


Wise2HSP_HSPmanagerFactory new_Wise2HSP_HSPmanagerFactory(PortableServer_POA poa, Wise2_HSPScanInterface * hsi,CORBA_Environment * ev);



#endif
