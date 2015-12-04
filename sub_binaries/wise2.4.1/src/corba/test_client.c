#include <string.h>
#include "hspscan_corba.h"
#include <stdio.h>
#include "sequence.h"

Wise2HSP_HSPmanagerFactory hspf;

Wise2HSP_HSPmanager * hspm;
Wise2HSP_Sequence query;


int
main (int argc, char *argv[])
{
    CORBA_Environment ev;
    CORBA_ORB orb;
    Wise2_Sequence * seq;
    int i;
    int j;

    FILE * ifp;
    char * ior;
    char filebuffer[1024];

    /*
     * Standard initalisation of the orb. Notice that
     * ORB_init 'eats' stuff off the command line
     */

    CORBA_exception_init(&ev);
    orb = CORBA_ORB_init(&argc, argv, "orbit-local-orb", &ev);

    /*
     * Get the IOR (object reference). It should be written out
     * by the echo-server into the file echo.ior. So - if you
     * are running the server in the same place as the client,
     * this should be fine!
     */

    seq = read_fasta_file_Sequence(argv[1]);

    ifp = fopen("hsp.ior","r");
    if( ifp == NULL ) {
      fprintf(stderr,"No echo.ior file!");
      exit(-1);
    }

    fgets(filebuffer,1023,ifp);
    ior = g_strdup(filebuffer);

    fclose(ifp);
    /*
     * Actually get the object. So easy!
     */

    hspf = CORBA_ORB_string_to_object(orb, ior, &ev);
    if (!hspf) {
	printf("Cannot bind to %s\n", ior);
	return 1;
    }

    query.name = CORBA_string_dup(seq->name);
    query.seq  = CORBA_string_dup(seq->seq);

    hspm = Wise2HSP_HSPmanagerFactory_scan_query(hspf,&query,&ev);
    

    for(i=0;i<hspm->target._length;i++) {
      fprintf(stdout,"Sequence %s\n",hspm->target._buffer[i].target.name);
      for(j=0;j<hspm->target._buffer[i].hsp._length;j++) {
	fprintf(stdout,"  hsp %d %d  %d\n",hspm->target._buffer[i].hsp._buffer[j].query_start,
		hspm->target._buffer[i].hsp._buffer[j].target_start,
		hspm->target._buffer[i].hsp._buffer[j].length);
      }
    }

    /* Clean up */
    CORBA_Object_release(hspf, &ev);
    CORBA_Object_release((CORBA_Object)orb, &ev);

    return 0;
}

