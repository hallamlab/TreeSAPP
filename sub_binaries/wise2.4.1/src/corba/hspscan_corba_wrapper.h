

#ifndef HSPSCAN_CORBA_WRAPPER_HEADER
#define HSPSCAN_CORBA_WRAPPER_HEADER

#include "hspscan_corba.h"
#include "corba_singleton.h"
#include "hspscaninterface.h"


/*
 * This method provides a hspscaninterface (ie, the Wise2 interface)
 * into the CORBA layer, abstracting out the calling convetion
 * of the ORB and remapping the results into the standard Wise2
 * structures (the over the wire system is more compact)
 */


Wise2_HSPScanInterface * new_corba_HSPScan(Wise2Corba_Singleton * sorb,char * iorfile,Wise2_CompMat * client_side_compmat);



#endif
