#ifdef _cplusplus
extern "C" {
#endif
#include "net_hspscan.h"


# line 20 "net_hspscan.dy"
HSPScanInterface * new_wise_transfer_HSPScanInterface(char * host,int port)
{
  HSPScanInterface * out;
  FunctionProxyCoordinator * fpc;

  out = HSPScanInterface_alloc();
  
  fpc = new_just_hspscan_FunctionProxyCoordinator(host,port);

  out->data = (void*) fpc;
  out->scan_query = dispatch_hspscan_FunctionProxyCoordinator;
  out->free_data  = free_hspscan_FunctionProxyCoordinator;
  
  return out;
}

# line 36 "net_hspscan.dy"
LinearHSPmanager * dispatch_hspscan_FunctionProxyCoordinator(void * h,Sequence * seq,HSPScanInterfacePara * para)
{
  FunctionProxyCoordinator * fpc = (FunctionProxyCoordinator*) h;
  LinearHSPmanager * lm;
  AnonymousObject * ao;
  AnonymousObjectList * aol;

  aol = AnonymousObjectList_alloc_len(2);
  ao = AnonymousObject_alloc();
  ao->obj = hard_link_Sequence(seq);
  ao->free_object = untyped_free_Sequence;

  add_AnonymousObjectList(aol,ao);

  ao = AnonymousObject_alloc();
  ao->obj = hard_link_HSPScanInterfacePara(para);
  ao->free_object = untyped_HSPScanInterfacePara_free;

  add_AnonymousObjectList(aol,ao);

  ao = dispatch_FunctionProxy(fpc,"hspscan_protein",aol);

  free_AnonymousObjectList(aol);

  /* transfer memory ownership to typed return value */
  lm = ao->obj;
  ao->obj = NULL;
  free_AnonymousObject(ao);


  return lm;

}

# line 70 "net_hspscan.dy"
void free_hspscan_FunctionProxyCoordinator(void * h)
{
  FunctionProxyCoordinator * fpc = (FunctionProxyCoordinator*) h;

  free_FunctionProxyCoordinator(fpc);
}


# line 78 "net_hspscan.dy"
FunctionProxyCoordinator * new_just_hspscan_FunctionProxyCoordinator(char * host,int port)
{
  FunctionProxyCoordinator * out;

  out = new_FunctionProxyCoordinator(host,port);

  add_FunctionProxyCoordinator(out,new_FunctionProxy(new_hspscan_protein_TransferedFunctionCall()));
  
  return out;
}

# line 89 "net_hspscan.dy"
TransferedObjectMarshaller * HSPScanInterfacePara_TransferedObjectMarshaller(void)
{
  TransferedObjectMarshaller * out;

  out = TransferedObjectMarshaller_alloc();
  out->object_to_stream = untyped_HSPScanInterfacePara_to_Stream;
  out->stream_to_object = untyped_read_HSPScanInterfacePara_from_Stream;
  out->object_type = "HSPScanInterfacePara";
  out->untyped_free_object = untyped_HSPScanInterfacePara_free;

  return out;
}
  

# line 103 "net_hspscan.dy"
TransferedObjectMarshaller * Sequence_TransferedObjectMarshaller(void)
{
  TransferedObjectMarshaller * out; 


  out = TransferedObjectMarshaller_alloc();
  out->object_to_stream = untyped_write_Sequence_to_Stream;
  out->stream_to_object = untyped_read_Sequence_from_Stream;
  out->object_type = "Sequence";
  out->untyped_free_object = untyped_free_Sequence;

  return out;
}


# line 118 "net_hspscan.dy"
TransferedObjectMarshaller * LinearHSPmanager_TransferedObjectMarshaller(void)
{
  TransferedObjectMarshaller * out; 

  out = TransferedObjectMarshaller_alloc();
  out->object_to_stream = untyped_LinearHSPmanager_to_Stream;

  out->stream_to_object = untyped_LinearHSPmanager_from_Stream;
  out->object_type = "LinearHSPmanager";
  out->untyped_free_object = untyped_free_LinearHSPmanager;


  return out;
}




# line 136 "net_hspscan.dy"
TransferedFunctionCall * new_hspscan_protein_TransferedFunctionCall(void)
{
  TransferedFunctionCall * out;

  out = TransferedFunctionCall_alloc_std();

  out->name = stringalloc("hspscan_protein");
  out->returned_type = LinearHSPmanager_TransferedObjectMarshaller();
  add_TransferedFunctionCall(out,Sequence_TransferedObjectMarshaller());
  add_TransferedFunctionCall(out,HSPScanInterfacePara_TransferedObjectMarshaller());

  return out;
}


# line 151 "net_hspscan.dy"
FunctionImplementation * new_hspscan_protein_FunctionImplementation(HSPScanInterface * hspi)
{
  FunctionImplementation * out;

  out = FunctionImplementation_alloc();

  out->transfer = new_hspscan_protein_TransferedFunctionCall();
  out->implementation = hspscan_protein_simple_impl;
  out->handle = hspi;

  return out;
}

# line 164 "net_hspscan.dy"
AnonymousObject * hspscan_protein_simple_impl(void * h,AnonymousObjectList * aol)
{
  HSPScanInterface * hsi = (HSPScanInterface*) h;
  Sequence * seq ;
  LinearHSPmanager * lm;
  AnonymousObject * ao;
  HSPScanInterfacePara * para;

  assert(aol != NULL);
  assert(aol->len == 2);

  seq  = (Sequence*)aol->anon[0]->obj;
  para = (HSPScanInterfacePara *)aol->anon[1]->obj;

  assert(seq != NULL);
  assert(para != NULL);
  

  lm = (*hsi->scan_query)(hsi->data,seq,para);
  lm->query = hard_link_Sequence(seq);


  ao = AnonymousObject_alloc();
  ao->obj = lm;
  ao->type = "LinearHSPmanager";
  ao->free_object = untyped_free_LinearHSPmanager;

  return ao;
}

# line 194 "net_hspscan.dy"
void untyped_free_Sequence(void * h)
{ 
  Sequence * s = (Sequence *) h;

  free_Sequence(s);
}

# line 201 "net_hspscan.dy"
void untyped_free_LinearHSPmanager(void * h)
{
  LinearHSPmanager * lm = (LinearHSPmanager *)h;

  free_LinearHSPmanager(lm);
}


# line 207 "net_hspscan.c"

#ifdef _cplusplus
}
#endif
