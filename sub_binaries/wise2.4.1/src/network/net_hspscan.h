#ifndef DYNAMITEnet_hspscanHEADERFILE
#define DYNAMITEnet_hspscanHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "hspscaninterface.h"
#include "functionserver.h"
#include "functionclient.h"

#include "sequencestream.h"
#include "hspstream.h"




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
HSPScanInterface * Wise2_new_wise_transfer_HSPScanInterface(char * host,int port);
#define new_wise_transfer_HSPScanInterface Wise2_new_wise_transfer_HSPScanInterface
LinearHSPmanager * Wise2_dispatch_hspscan_FunctionProxyCoordinator(void * h,Sequence * seq,HSPScanInterfacePara * para);
#define dispatch_hspscan_FunctionProxyCoordinator Wise2_dispatch_hspscan_FunctionProxyCoordinator
void Wise2_free_hspscan_FunctionProxyCoordinator(void * h);
#define free_hspscan_FunctionProxyCoordinator Wise2_free_hspscan_FunctionProxyCoordinator
FunctionProxyCoordinator * Wise2_new_just_hspscan_FunctionProxyCoordinator(char * host,int port);
#define new_just_hspscan_FunctionProxyCoordinator Wise2_new_just_hspscan_FunctionProxyCoordinator
TransferedObjectMarshaller * Wise2_HSPScanInterfacePara_TransferedObjectMarshaller(void);
#define HSPScanInterfacePara_TransferedObjectMarshaller Wise2_HSPScanInterfacePara_TransferedObjectMarshaller
TransferedObjectMarshaller * Wise2_Sequence_TransferedObjectMarshaller(void);
#define Sequence_TransferedObjectMarshaller Wise2_Sequence_TransferedObjectMarshaller
TransferedObjectMarshaller * Wise2_LinearHSPmanager_TransferedObjectMarshaller(void);
#define LinearHSPmanager_TransferedObjectMarshaller Wise2_LinearHSPmanager_TransferedObjectMarshaller
TransferedFunctionCall * Wise2_new_hspscan_protein_TransferedFunctionCall(void);
#define new_hspscan_protein_TransferedFunctionCall Wise2_new_hspscan_protein_TransferedFunctionCall
FunctionImplementation * Wise2_new_hspscan_protein_FunctionImplementation(HSPScanInterface * hspi);
#define new_hspscan_protein_FunctionImplementation Wise2_new_hspscan_protein_FunctionImplementation
AnonymousObject * Wise2_hspscan_protein_simple_impl(void * h,AnonymousObjectList * aol);
#define hspscan_protein_simple_impl Wise2_hspscan_protein_simple_impl
void Wise2_untyped_free_Sequence(void * h);
#define untyped_free_Sequence Wise2_untyped_free_Sequence
void Wise2_untyped_free_LinearHSPmanager(void * h);
#define untyped_free_LinearHSPmanager Wise2_untyped_free_LinearHSPmanager


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/

#ifdef _cplusplus
}
#endif

#endif
