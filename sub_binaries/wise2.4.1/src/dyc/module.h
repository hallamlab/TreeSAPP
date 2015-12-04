/*  Last edited: Apr 29 10:15 1997 (birney) */
/********************************************************/
/* Wisetools version 3 file                             */
/*                                                      */
/* this file is copyright (c) Ewan Birney 1996          */
/* this file is part of the wisetools sequence analysis */
/* package                                              */
/********************************************************/

/********************************************************/
/* Header file for Dynamite structure - a wrapper around*/
/* Dynamite structures                                  */
/*                                                      */
/*                                                      */
/********************************************************/

/***** RCS Info *****************************************/
/*
   $Id: module.h,v 1.1.1.1 2001/06/18 13:59:56 birney Exp $

   $Log: module.h,v $
   Revision 1.1.1.1  2001/06/18 13:59:56  birney
   moved wise2 to ensembl cvs repository

   Revision 1.3  1999/05/03 14:52:59  birney
   added in alot of fixes to dyanmite.

   added in def boolean's in commandline

   added in seqalign/pwmdna stuff and linked it to hmmer bridge.

   Revision 1.2  1998/11/16 15:25:00  birney
   started the threads .... ;)

 * Revision 1.1.1.1  1998/08/28  09:30:58  birney
 * Wise2
 *
 * Revision 1.7  1998/07/27  11:39:36  birney
 * before compugen merge
 *
 * Revision 1.6  1998/01/09  15:12:45  birney
 * Added api handling, namespace protection and pdoc production. Lots of bugs
 * out from the main matrix construction including start point errors on
 * packaln
 *
 * Revision 1.5  1997/07/24  15:40:48  birney
 * added MethodTypeSet support for on-the-fly methods
 *
 * Revision 1.4  1997/06/10  13:55:37  birney
 * Lots of things: moved to teh dynfile system of keeping files around, addedfunction information, started to remodelled the matrix part
 *
 * Revision 1.3  1996/10/14  10:38:44  birney
 * added display handling
 *
 * Revision 1.2  1996/10/04  16:42:00  birney
 * added cons/decons control for basic structures
 *
 * Revision 1.1  1996/03/04  12:15:52  birney
 * Initial revision
 *

*/
/********************************************************/

#ifndef MODULEHEADER
#define MODULEHEADER

#define WRITECHEADER

#include "wisec.h"
#include "dyna2.h"
#include "display.h"
#include "dynafunc.h"
#include "friend.h"
#include "api.h"

typedef struct {  
	void * data;             
	int type;                
	}  DynamiteHolder;      /* DynamiteHolder defined */ 


typedef struct {  
	DynamiteHolder ** list;  
	int len;             /* len for above list  */ 
	int maxlen;          /* maxlen for above list */ 
	char * name;          
	}  Dynamite;         /* Dynamite defined */ 

#define DynamiteLISTLENGTH 32

#define DYNAMITETYPEGENERICMATRIX 56
#define DYNAMITETYPESTRUCT        57
#define DYNAMITETYPEFRONTEND      58
#define DYNAMITETYPEPARSER        59
#define DYNAMITETYPEDISPLAY       60
#define DYNAMITETYPEFRIEND        61
#define DYNAMITETYPEAPI           62




   /* prototypes */ 

boolean  prepare_Dynamite(Dynamite * dyn,MethodTypeSet * mts,DycWarning * dycw,ParseError fail);

dynAPI * promote_all_callable_dynAPI(Dynamite * dyn,DYNFILE * dfp);
boolean crosslink_to_StructHolder_dynAPI(dynAPI * api,Dynamite * dyn);

boolean prepare_dynAPI(dynAPI * api,DYNFILE * dfp,Dynamite * dyn);

void write_Dynamite_func(Dynamite * dyn,ModuleFunctionList * mfl,DPImplementation * dpi,DYNFILE * dfp);
void write_Dynamite_header(Dynamite * dyn,DPImplementation * dpi,DYNFILE * dfp);
Dynamite * read_Dynamite(FILE * ifp,char * fend,boolean * end_of_file,MethodTypeSet * mts);
DynamiteHolder * wrap_Aln2Display(Aln2Display * dis);
DynamiteHolder * wrap_Friend(Friend * fr);
DynamiteHolder * wrap_GenericMatrix(GenericMatrix * wm);
DynamiteHolder * wrap_StructHolder(StructHolder * wm);
DynamiteHolder * wrap_dynAPI(dynAPI * api);
DynamiteHolder * DynamiteHolder_alloc(void) ;
DynamiteHolder * free_DynamiteHolder(DynamiteHolder * obj) ;
void swap_Dynamite(DynamiteHolder ** list,int i,int j)  ;
void qsort_Dynamite(DynamiteHolder ** list,int left,int right,int (*comp)(DynamiteHolder * ,DynamiteHolder * )) ;
void sort_Dynamite(Dynamite * obj,int (*comp)(DynamiteHolder *, DynamiteHolder *)) ;
boolean add_Dynamite(Dynamite * obj,DynamiteHolder * add) ;
boolean expand_Dynamite(Dynamite * obj,int len) ;
Dynamite * Dynamite_alloc_std(void) ;
Dynamite * Dynamite_alloc_len(int len) ;
Dynamite * Dynamite_alloc(void) ;
Dynamite * free_Dynamite(Dynamite * obj) ;


#endif /* module loaded */






