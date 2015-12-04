#include "hspscan_corba.h"

#include "hspscaninterface.h"
/*** App-specific servant structures ***/


typedef struct {
  POA_Wise2HSP_HSP_ts_iterator servant;
  PortableServer_POA poa;
  LinearHSPmanager * lm;
  int current_pos;
  int max_size;
} impl_POA_Wise2HSP_HSP_ts_iterator;


typedef struct {
  POA_Wise2HSP_HSPmanagerFactory servant;
  PortableServer_POA poa;
  HSPScanInterface * hsi;
  int is_verbose;
} impl_POA_Wise2HSP_HSPmanagerFactory;




/*** Implementation stub prototypes ***/


static void impl_Wise2HSP_HSP_ts_iterator__destroy(impl_POA_Wise2HSP_HSP_ts_iterator * servant,
						   CORBA_Environment * ev);
static Wise2HSP_HSPtargetset *
 impl_Wise2HSP_HSP_ts_iterator_next(impl_POA_Wise2HSP_HSP_ts_iterator * servant,
				    CORBA_Environment * ev);

static CORBA_boolean
 impl_Wise2HSP_HSP_ts_iterator_is_more(impl_POA_Wise2HSP_HSP_ts_iterator * servant,
				       CORBA_Environment * ev);





static void impl_Wise2HSP_HSPmanagerFactory__destroy(impl_POA_Wise2HSP_HSPmanagerFactory *servant,
CORBA_Environment *ev);
static Wise2HSP_HSPmanager* impl_Wise2HSP_HSPmanagerFactory_scan_query(impl_POA_Wise2HSP_HSPmanagerFactory *servant,
Wise2HSP_Sequence* query,
CORBA_Environment *ev);

static Wise2HSP_HSP_ts_iterator
 impl_Wise2HSP_HSPmanagerFactory_scan_query_iterator(impl_POA_Wise2HSP_HSPmanagerFactory * servant,
						     Wise2HSP_Sequence * query,CORBA_short max_size,
						     CORBA_Environment * ev);

/*** epv structures ***/

static PortableServer_ServantBase__epv impl_Wise2HSP_HSP_ts_iterator_base_epv =
{
   NULL,			/* _private data */
   NULL,			/* finalize routine */
   NULL,			/* default_POA routine */
};
static POA_Wise2HSP_HSP_ts_iterator__epv impl_Wise2HSP_HSP_ts_iterator_epv =
{
   NULL,			/* _private */
   (gpointer) & impl_Wise2HSP_HSP_ts_iterator_next,

   (gpointer) & impl_Wise2HSP_HSP_ts_iterator_is_more,

};
static PortableServer_ServantBase__epv impl_Wise2HSP_HSPmanagerFactory_base_epv =
{
   NULL,			/* _private data */
   NULL,			/* finalize routine */
   NULL,			/* default_POA routine */
};
static POA_Wise2HSP_HSPmanagerFactory__epv impl_Wise2HSP_HSPmanagerFactory_epv =
{
   NULL,			/* _private */
   (gpointer) & impl_Wise2HSP_HSPmanagerFactory_scan_query,

   (gpointer) & impl_Wise2HSP_HSPmanagerFactory_scan_query_iterator,

};

/*** vepv structures ***/

static POA_Wise2HSP_HSP_ts_iterator__vepv impl_Wise2HSP_HSP_ts_iterator_vepv =
{
   &impl_Wise2HSP_HSP_ts_iterator_base_epv,
   &impl_Wise2HSP_HSP_ts_iterator_epv,
};
static POA_Wise2HSP_HSPmanagerFactory__vepv impl_Wise2HSP_HSPmanagerFactory_vepv =
{
   &impl_Wise2HSP_HSPmanagerFactory_base_epv,
   &impl_Wise2HSP_HSPmanagerFactory_epv,
};


/*** iterator implementations ****/


/*** Stub implementations ***/


static Wise2HSP_HSP_ts_iterator 
impl_Wise2HSP_HSP_ts_iterator__create(PortableServer_POA poa, LinearHSPmanager * lm,int max_size,CORBA_Environment * ev)
{
   Wise2HSP_HSP_ts_iterator retval;
   impl_POA_Wise2HSP_HSP_ts_iterator *newservant;
   PortableServer_ObjectId *objid;

   newservant = g_new0(impl_POA_Wise2HSP_HSP_ts_iterator, 1);
   newservant->servant.vepv = &impl_Wise2HSP_HSP_ts_iterator_vepv;
   newservant->poa = poa;
   newservant->lm  = lm;
   newservant->max_size = max_size;
   newservant->current_pos = 0;

   POA_Wise2HSP_HSP_ts_iterator__init((PortableServer_Servant) newservant, ev);
   objid = PortableServer_POA_activate_object(poa, newservant, ev);
   CORBA_free(objid);
   retval = PortableServer_POA_servant_to_reference(poa, newservant, ev);

   return retval;
}

static void
impl_Wise2HSP_HSP_ts_iterator__destroy(impl_POA_Wise2HSP_HSP_ts_iterator * servant, CORBA_Environment * ev)
{
   PortableServer_ObjectId *objid;

   objid = PortableServer_POA_servant_to_id(servant->poa, servant, ev);
   PortableServer_POA_deactivate_object(servant->poa, objid, ev);
   CORBA_free(objid);

   POA_Wise2HSP_HSP_ts_iterator__fini((PortableServer_Servant) servant, ev);
   g_free(servant);
}

static Wise2HSP_HSPtargetset *
impl_Wise2HSP_HSP_ts_iterator_next(impl_POA_Wise2HSP_HSP_ts_iterator * servant,
				   CORBA_Environment * ev)
{
   Wise2HSP_HSPtargetset *retval;
   HSPset * s;
   int i;


   s = servant->lm->set[servant->current_pos];

   /*   fprintf(stderr,"Going to pass to client %d block %s\n",servant->current_pos,s->hsp[0]->target->name);*/

   retval = Wise2HSP_HSPtargetset__alloc();

   /* do we need to worry about string dup's ... probably!*/
   retval->target.name = CORBA_string_dup(s->hsp[0]->target->name);
   retval->target.seq  = CORBA_string_dup(s->hsp[0]->target->seq);
   
   retval->hsp._buffer =  CORBA_sequence_Wise2HSP_HSPinfo_allocbuf(s->len);
   retval->hsp._maximum = retval->hsp._length = s->len;
   
  for(i=0;i<s->len;i++) {
    retval->hsp._buffer[i].query_start   = s->hsp[i]->query_start;
    retval->hsp._buffer[i].target_start  = s->hsp[i]->target_start;
    retval->hsp._buffer[i].length        = s->hsp[i]->length;
    retval->hsp._buffer[i].score         = s->hsp[i]->score;
  }

   servant->current_pos++;

   /*   fprintf(stderr,"Going to return with score %d length %d\n",s->score,s->len); */

   return retval;
}

static CORBA_boolean
impl_Wise2HSP_HSP_ts_iterator_is_more(impl_POA_Wise2HSP_HSP_ts_iterator * servant,
				      CORBA_Environment * ev)
{
   CORBA_boolean retval;


   if( servant->current_pos > servant->max_size ) {
     return FALSE;
   }

   if( servant->current_pos < servant->lm->len ) {
     retval = TRUE;
   } else {
     retval = FALSE;
   }

   return retval;
}


/*** factory implementations ****/



static Wise2HSP_HSPmanagerFactory impl_Wise2HSP_HSPmanagerFactory__create(PortableServer_POA poa, CORBA_Environment *ev)
{
  Wise2HSP_HSPmanagerFactory retval;
  impl_POA_Wise2HSP_HSPmanagerFactory *newservant;
  PortableServer_ObjectId *objid;
  
  newservant = g_new0(impl_POA_Wise2HSP_HSPmanagerFactory, 1);
  newservant->servant.vepv = &impl_Wise2HSP_HSPmanagerFactory_vepv;
  newservant->poa = poa;
  POA_Wise2HSP_HSPmanagerFactory__init((PortableServer_Servant)newservant, ev);
  objid = PortableServer_POA_activate_object(poa, newservant, ev);
  CORBA_free(objid);
  retval = PortableServer_POA_servant_to_reference(poa, newservant, ev);
  
  return retval;
}


Wise2HSP_HSPmanagerFactory new_Wise2HSP_HSPmanagerFactory(PortableServer_POA poa, Wise2_HSPScanInterface * hsi,CORBA_Environment * ev)
{
  Wise2HSP_HSPmanagerFactory retval;
  impl_POA_Wise2HSP_HSPmanagerFactory *newservant;
  PortableServer_ObjectId *objid;
  
  newservant = g_new0(impl_POA_Wise2HSP_HSPmanagerFactory, 1);
  newservant->servant.vepv = &impl_Wise2HSP_HSPmanagerFactory_vepv;
  newservant->poa = poa;
  newservant->hsi = hsi;
  newservant->is_verbose = 1;

  POA_Wise2HSP_HSPmanagerFactory__init((PortableServer_Servant)newservant, ev);
  objid = PortableServer_POA_activate_object(poa, newservant, ev);
  CORBA_free(objid);
  retval = PortableServer_POA_servant_to_reference(poa, newservant, ev);
  
  return retval;
}



static void
impl_Wise2HSP_HSPmanagerFactory__destroy(impl_POA_Wise2HSP_HSPmanagerFactory *servant, CORBA_Environment *ev)
{
  PortableServer_ObjectId *objid;
  
  objid = PortableServer_POA_servant_to_id(servant->poa, servant, ev);
  PortableServer_POA_deactivate_object(servant->poa, objid, ev);
  CORBA_free(objid);
  
  POA_Wise2HSP_HSPmanagerFactory__fini((PortableServer_Servant)servant, ev);
  g_free(servant);
}


static void remap_hsp(gpointer key,gpointer data,gpointer user_data)
{
  HSPset * set = (HSPset*)data;
  Wise2HSP_HSPmanager* manager = (Wise2HSP_HSPmanager*) user_data;
  int pos;

  int i;

  pos = manager->target._length++;

  /* first off, allocate and store this target */

  manager->target._buffer[pos].target.name = CORBA_string_dup(set->hsp[0]->target->name);
  manager->target._buffer[pos].target.seq  = CORBA_string_dup(set->hsp[0]->target->seq);
  
  manager->target._buffer[pos].hsp._buffer  = CORBA_sequence_Wise2HSP_HSPinfo_allocbuf(set->len);
  manager->target._buffer[pos].hsp._maximum =set->len;
  manager->target._buffer[pos].hsp._length  =set->len;
  
  for(i=0;i<set->len;i++) {
    manager->target._buffer[pos].hsp._buffer[i].query_start   = set->hsp[i]->query_start;
    manager->target._buffer[pos].hsp._buffer[i].target_start  = set->hsp[i]->target_start;
    manager->target._buffer[pos].hsp._buffer[i].length        = set->hsp[i]->length;
    manager->target._buffer[pos].hsp._buffer[i].score         = set->hsp[i]->score;
  }
    
  return;
}


static Wise2HSP_HSPmanager * corba_manager_from_wise_manager(LinearHSPmanager * lm)
{
  Wise2HSP_HSPmanager* retval;
  int size;
  int i;
  int pos;

  size = lm->len;

  /* ok, we need to build the return structure now */

  retval = Wise2HSP_HSPmanager__alloc();
  retval->target._buffer  = CORBA_sequence_Wise2HSP_HSPtargetset_allocbuf(size);
  retval->target._length  = 0;
  retval->target._maximum = size;


  for(pos=0;pos<lm->len;pos++) {
    auto HSPset * set = lm->set[pos];

    retval->target._length++;
    
    /* first off, allocate and store this target */

    retval->target._buffer[pos].target.name = CORBA_string_dup(set->hsp[0]->target->name);
    retval->target._buffer[pos].target.seq  = CORBA_string_dup(set->hsp[0]->target->seq);
  
    retval->target._buffer[pos].hsp._buffer  = CORBA_sequence_Wise2HSP_HSPinfo_allocbuf(set->len);
    retval->target._buffer[pos].hsp._maximum =set->len;
    retval->target._buffer[pos].hsp._length  =set->len;
    
    for(i=0;i<set->len;i++) {
      retval->target._buffer[pos].hsp._buffer[i].query_start   = set->hsp[i]->query_start;
      retval->target._buffer[pos].hsp._buffer[i].target_start  = set->hsp[i]->target_start;
      retval->target._buffer[pos].hsp._buffer[i].length        = set->hsp[i]->length;
      retval->target._buffer[pos].hsp._buffer[i].score         = set->hsp[i]->score;
    }
    
  }

  return retval;
}

static Wise2HSP_HSPmanager*
impl_Wise2HSP_HSPmanagerFactory_scan_query(impl_POA_Wise2HSP_HSPmanagerFactory *servant,
Wise2HSP_Sequence* query,
					   CORBA_Environment *ev)
{
  Wise2HSP_HSPmanager* retval;
  Wise2_Sequence * wise2_seq;
  LinearHSPmanager * lm;
  HSPScanInterfacePara para;

  int size;
  
  if( servant->is_verbose == 1 ) {
    fprintf(stdout,"going to do %s\n",query->name);
  }

  wise2_seq = new_Sequence_from_strings(query->name,query->seq);

  para.max_results = 250;

  lm = (*servant->hsi->scan_query)(servant->hsi->data,wise2_seq,&para);

  retval = corba_manager_from_wise_manager(lm);

  if( servant->is_verbose == 1 ) {
    fprintf(stdout," ... remapped HSP manager\n");
  }

  Wise2_free_LinearHSPmanager(lm);

  if( servant->is_verbose == 1 ) {
    fprintf(stdout," ... freed data, about to return\n");
  }



  return retval;
}


static Wise2HSP_HSP_ts_iterator
impl_Wise2HSP_HSPmanagerFactory_scan_query_iterator(impl_POA_Wise2HSP_HSPmanagerFactory * servant,
						    Wise2HSP_Sequence * query,CORBA_short max_size,
						    CORBA_Environment * ev)
{
   Wise2HSP_HSP_ts_iterator retval;

   Wise2HSP_HSPmanager* corba_hspm;
   Wise2_Sequence * wise2_seq;

   LinearHSPmanager * lm;
   HSPScanInterfacePara para;
   int size;
   int i;
   
  if( servant->is_verbose == 1 ) {
    fprintf(stdout,"going to do %s\n",query->name);
  }

  wise2_seq = new_Sequence_from_strings(query->name,query->seq);

  para.max_results = max_size;

  lm = (*servant->hsi->scan_query)(servant->hsi->data,wise2_seq,&para);


  /*
  fprintf(stderr,"Sanity check in corba\n");
  for(i=0;i<lm->len;i++) {
    fprintf(stderr,"Got %s %d\n",lm->set[i]->hsp[0]->target->name,i);
  }
  */

  if( servant->is_verbose == 1 ) {
    fprintf(stdout," ... retrieved HSP manager size %d\n",lm->len);
  }

  if( servant->is_verbose == 1 ) {
    fprintf(stdout," ... remapped HSP manager\n");
  }
  
  if( servant->is_verbose == 1 ) {
    fprintf(stdout," ... freed data, about to return iterator\n");
  }

  retval = impl_Wise2HSP_HSP_ts_iterator__create(servant->poa,lm,max_size,ev);
  
  return retval;
}


