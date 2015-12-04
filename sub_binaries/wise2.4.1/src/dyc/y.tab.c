
# line 3 "calc.y"
#include <stdio.h>
#include "exprtree.h"

extern ExprTree * root;
extern char * calc_lex_string;
extern int stringpos;


# line 12 "calc.y"
typedef union  {
	double dval;
	struct ExprTree *  tr;
} YYSTYPE;
# define NAME 257
# define NUMBER 258
# define STRUCTREF 259
# define DEREFERENCE 260
#define yyclearin yychar = -1
#define yyerrok yyerrflag = 0
extern int yychar;
extern int yyerrflag;
#ifndef YYMAXDEPTH
#define YYMAXDEPTH 150
#endif
YYSTYPE yylval, yyval;
typedef int yytabelem;
# define YYERRCODE 256

# line 72 "calc.y"


#ifdef YYWRAP_NEEDED
void yywrap(void) 
{
  fprintf(stderr,"In yywrap... going to exit badly");
  exit(1);
}
#endif

void yyerror(char * s) 
{
  /*** stringpos is position along the string ***/

  int i;

  warn("Calc line parser error: [%s]",s);
  fprintf(stderr,"Occured at:\n");
  fprintf(stderr,"%s\n",calc_lex_string);
  for(i=0;i<stringpos-1;i++) {
    fputc('-',stderr);
  }
  fputc('^',stderr);
  fputc('\n',stderr);

  root = NULL; /*** fuck - memory !!!! ***/

}



static const yytabelem yyexca[] ={
-1, 1,
	0, -1,
	-2, 0,
	};
# define YYNPROD 19
# define YYLAST 230

static const yytabelem yyact[]={

     8,    17,     3,    26,     9,     8,    12,     3,    15,     9,
     8,    13,    12,    10,     9,    11,    25,    13,    24,    12,
    10,     6,    11,     2,    13,    12,    10,    14,    11,    30,
    13,     1,    31,     5,    20,    21,    22,    23,     0,    27,
     0,    29,    18,    19,     0,     0,     0,     0,     0,     0,
    28,     0,     0,     0,     0,    33,     0,     0,     0,    17,
     0,     0,     0,    32,     0,     0,     0,     0,     0,     0,
     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     0,     0,     0,     0,     0,     0,     0,     0,     0,    16,
     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     0,     0,     0,     0,     0,     0,     0,     0,     0,     7,
     4,     0,     0,     0,     7,     4,     0,    16,     0,     7 };

static const yytabelem yypact[]={

   -33, -1000,   -17,   -33, -1000,   -32, -1000, -1000,   -28,   -28,
   -33,   -33,   -33,   -33,   -23,   -38,   -28,   -33,   -90,   -90,
   -36,   -36, -1000, -1000, -1000,   -12, -1000,   -17, -1000,   -30,
 -1000,   -33, -1000,   -17 };

static const yytabelem yypgo[]={

     0,    23,    31,    33,    21,    16 };

static const yytabelem yyr1[]={

     0,     2,     4,     4,     1,     1,     1,     1,     1,     1,
     1,     1,     5,     5,     3,     3,     3,     3,     3 };

static const yytabelem yyr2[]={

     0,     3,     9,     7,     7,     7,     7,     7,     7,     3,
     3,     3,     7,     3,     3,     7,     9,     5,     5 };

static const yytabelem yychk[]={

 -1000,    -2,    -1,    40,   258,    -3,    -4,   257,    38,    42,
    43,    45,    42,    47,    -1,    40,   259,    91,    -3,    -3,
    -1,    -1,    -1,    -1,    41,    -5,    41,    -1,    -3,    -1,
    41,    44,    93,    -1 };

static const yytabelem yydef[]={

     0,    -2,     1,     0,     9,    10,    11,    14,     0,     0,
     0,     0,     0,     0,     0,     0,     0,     0,    17,    18,
     4,     5,     6,     7,     8,     0,     3,    13,    15,     0,
     2,     0,    16,    12 };
typedef struct { char *t_name; int t_val; } yytoktype;
#ifndef YYDEBUG
#	define YYDEBUG	1	/* allow debugging */
#endif

#if YYDEBUG

yytoktype yytoks[] =
{
	"NAME",	257,
	"NUMBER",	258,
	"STRUCTREF",	259,
	"-",	45,
	"+",	43,
	"*",	42,
	"/",	47,
	"&",	38,
	"DEREFERENCE",	260,
	"[",	91,
	"]",	93,
	"-unknown-",	-1	/* ends search */
};

char * yyreds[] =
{
	"-no such reduction-",
      "statement : expression",
      "method : tag '(' commalist ')'",
      "method : tag '(' ')'",
      "expression : expression '+' expression",
      "expression : expression '-' expression",
      "expression : expression '*' expression",
      "expression : expression '/' expression",
      "expression : '(' expression ')'",
      "expression : NUMBER",
      "expression : tag",
      "expression : method",
      "commalist : commalist ',' expression",
      "commalist : expression",
      "tag : NAME",
      "tag : tag STRUCTREF tag",
      "tag : tag '[' expression ']'",
      "tag : '&' tag",
      "tag : '*' tag",
};
#endif /* YYDEBUG */
/*
 * *****************************************************************
 * *                                                               *
 * *   Copyright 2002 Compaq Information Technologies Group, L.P.  *
 * *                                                               *
 * *   The software contained on this media  is  proprietary  to   *
 * *   and  embodies  the  confidential  technology  of  Compaq    *
 * *   Computer Corporation.  Possession, use,  duplication  or    *
 * *   dissemination of the software and media is authorized only  *
 * *   pursuant to a valid written license from Compaq Computer    *
 * *   Corporation.                                                *
 * *                                                               *
 * *   RESTRICTED RIGHTS LEGEND   Use, duplication, or disclosure  *
 * *   by the U.S. Government is subject to restrictions  as  set  *
 * *   forth in Subparagraph (c)(1)(ii)  of  DFARS  252.227-7013,  *
 * *   or  in  FAR 52.227-19, as applicable.                       *
 * *                                                               *
 * *****************************************************************
 */
/*
 * HISTORY
 */
/*
 * @(#)$RCSfile: y.tab.c,v $ $Revision: 1.11 $ (DEC) $Date: 2006/04/02 08:53:40 $
 */
/*
** Skeleton parser driver for yacc output
*/

/* external references for c++ and ANSI C
 * Define YY_NOPROTO to suppress the prototype declarations
 * GNUC and DECC define __STDC__ differently
 */
#ifdef __GNUC__
#if !__STDC__
#define YY_NOPROTO
#endif /* __STDC__ */
#elif !defined(__STDC__) &&  !defined (__cplusplus)
#define YY_NOPROTO
#endif /* __STDC__ */

/* Disable array out of bounds info messages */
#if defined (__DECC)
#pragma message disable (badsubscript,subscrbounds,unreach)
#endif

#ifndef YY_NOPROTO
#if defined (__cplusplus)
extern "C" {
extern void yyerror(char *); 
extern int yylex();
#else /* __cplusplus */
extern int yylex(void);
#endif /* __cplusplus */
#if defined (__cplusplus)
}
#endif /* __cplusplus */
#endif /* YY_NOPROTO */
/*
** yacc user known macros and defines
*/
#ifdef YYSPLIT
#   define YYERROR	return(-2)
#else
#   define YYERROR	goto yyerrlab
#endif

#define YYACCEPT	return(0)
#define YYABORT		return(1)
#define YYBACKUP( newtoken, newvalue )\
{\
	if ( yychar >= 0 || ( yyr2[ yytmp ] >> 1 ) != 1 )\
	{\
		yyerror( "syntax error - cannot backup" );\
		goto yyerrlab;\
	}\
	yychar = newtoken;\
	yystate = *yyps;\
	yylval = newvalue;\
	goto yynewstate;\
}
#define YYRECOVERING()	(!!yyerrflag)
#ifndef YYDEBUG
#	define YYDEBUG	1	/* make debugging available */
#endif

/*
** user known globals
*/
int yydebug;			/* set to 1 to get debugging */

/*
** driver internal defines
*/
#define YYFLAG		(-1000)

#ifdef YYSPLIT
#   define YYSCODE { \
			extern int (*yyf[])(); \
			register int yyret; \
			if (yyf[yytmp]) \
			    if ((yyret=(*yyf[yytmp])()) == -2) \
				    goto yyerrlab; \
				else if (yyret>=0) return(yyret); \
		   }
#endif

/*
** local variables used by the parser
 * these should be static in order to support
 * multiple parsers in a single executable program. POSIX 1003.2-1993
 */
static YYSTYPE yyv[ YYMAXDEPTH ];	/* value stack */
static int yys[ YYMAXDEPTH ];		/* state stack */

static YYSTYPE *yypv;			/* top of value stack */
static YYSTYPE *yypvt;			/* top of value stack for $vars */
static int *yyps;			/* top of state stack */

static int yystate;			/* current state */
static int yytmp;			/* extra var (lasts between blocks) */

/*
** global variables used by the parser - renamed as a result of -p
*/
int yynerrs;			/* number of errors */
int yyerrflag;			/* error recovery flag */
int yychar;			/* current input token number */

/*
** yyparse - return 0 if worked, 1 if syntax error not recovered from
*/
int
yyparse()
{
	/*
	** Initialize externals - yyparse may be called more than once
	*/
	yypv = &yyv[-1];
	yyps = &yys[-1];
	yystate = 0;
	yytmp = 0;
	yynerrs = 0;
	yyerrflag = 0;
	yychar = -1;

	goto yystack;
	{
		register YYSTYPE *yy_pv;	/* top of value stack */
		register int *yy_ps;		/* top of state stack */
		register int yy_state;		/* current state */
		register int  yy_n;		/* internal state number info */

		/*
		** get globals into registers.
		** branch to here only if YYBACKUP was called.
		*/
	yynewstate:
		yy_pv = yypv;
		yy_ps = yyps;
		yy_state = yystate;
		goto yy_newstate;

		/*
		** get globals into registers.
		** either we just started, or we just finished a reduction
		*/
	yystack:
		yy_pv = yypv;
		yy_ps = yyps;
		yy_state = yystate;

		/*
		** top of for (;;) loop while no reductions done
		*/
	yy_stack:
		/*
		** put a state and value onto the stacks
		*/
#if YYDEBUG
		/*
		** if debugging, look up token value in list of value vs.
		** name pairs.  0 and negative (-1) are special values.
		** Note: linear search is used since time is not a real
		** consideration while debugging.
		*/
		if ( yydebug )
		{
			register int yy_i;

			printf( "State %d, token ", yy_state );
			if ( yychar == 0 )
				printf( "end-of-file\n" );
			else if ( yychar < 0 )
				printf( "-none-\n" );
			else
			{
				for ( yy_i = 0; yytoks[yy_i].t_val >= 0;
					yy_i++ )
				{
					if ( yytoks[yy_i].t_val == yychar )
						break;
				}
				printf( "%s\n", yytoks[yy_i].t_name );
			}
		}
#endif /* YYDEBUG */
		if ( ++yy_ps >= &yys[ YYMAXDEPTH ] )	/* room on stack? */
		{
			yyerror( "yacc stack overflow" );
			YYABORT;
		}
		*yy_ps = yy_state;
		*++yy_pv = yyval;

		/*
		** we have a new state - find out what to do
		*/
	yy_newstate:
		if ( ( yy_n = yypact[ yy_state ] ) <= YYFLAG )
			goto yydefault;		/* simple state */
#if YYDEBUG
		/*
		** if debugging, need to mark whether new token grabbed
		*/
		yytmp = yychar < 0;
#endif
		if ( ( yychar < 0 ) && ( ( yychar = yylex() ) < 0 ) )
			yychar = 0;		/* reached EOF */
#if YYDEBUG
		if ( yydebug && yytmp )
		{
			register int yy_i;

			printf( "Received token " );
			if ( yychar == 0 )
				printf( "end-of-file\n" );
			else if ( yychar < 0 )
				printf( "-none-\n" );
			else
			{
				for ( yy_i = 0; yytoks[yy_i].t_val >= 0;
					yy_i++ )
				{
					if ( yytoks[yy_i].t_val == yychar )
						break;
				}
				printf( "%s\n", yytoks[yy_i].t_name );
			}
		}
#endif /* YYDEBUG */
		if ( ( ( yy_n += yychar ) < 0 ) || 
               ( yy_n >= YYLAST )       || 
               (yyact[yy_n ] < 0))
               	goto yydefault;

		if ( yychk[ yy_n = yyact[ yy_n ] ] == yychar )	/*valid shift*/
		{
			yychar = -1;
			yyval = yylval;
			yy_state = yy_n;
			if ( yyerrflag > 0 )
				yyerrflag--;
			goto yy_stack;
		}

	yydefault:
		if ( ( yy_n = yydef[ yy_state ] ) == -2 )
		{
#if YYDEBUG
			yytmp = yychar < 0;
#endif
			if ( ( yychar < 0 ) && ( ( yychar = yylex() ) < 0 ) )
				yychar = 0;		/* reached EOF */
#if YYDEBUG
			if ( yydebug && yytmp )
			{
				register int yy_i;

				printf( "Received token " );
				if ( yychar == 0 )
					printf( "end-of-file\n" );
				else if ( yychar < 0 )
					printf( "-none-\n" );
				else
				{
					for ( yy_i = 0;
						yytoks[yy_i].t_val >= 0;
						yy_i++ )
					{
						if ( yytoks[yy_i].t_val
							== yychar )
						{
							break;
						}
					}
					printf( "%s\n", yytoks[yy_i].t_name );
				}
			}
#endif /* YYDEBUG */
			/*
			** look through exception table
			*/
			{
				register const int *yyxi = yyexca;

				while ( ( *yyxi != -1 ) ||
					( yyxi[1] != yy_state ) )
				{
					yyxi += 2;
				}
				while ( ( *(yyxi += 2) >= 0 ) &&
					( *yyxi != yychar ) )
					;
				if ( ( yy_n = yyxi[1] ) < 0 )
					YYACCEPT;
			}
		}

		/*
		** check for syntax error
		*/
		if ( yy_n == 0 )	/* have an error */
		{
			/* no worry about speed here! */
			switch ( yyerrflag )
			{
			case 0:		/* new error */
				yyerror( "syntax error" );
				goto skip_init;
			yyerrlab:
				/*
				** get globals into registers.
				** we have a user generated syntax type error
				*/
				yy_pv = yypv;
				yy_ps = yyps;
				yy_state = yystate;
				yynerrs++;
			skip_init:
			case 1:
			case 2:		/* incompletely recovered error */
					/* try again... */
				yyerrflag = 3;
				/*
				** find state where "error" is a legal
				** shift action
				*/
				while ( yy_ps >= yys )
				{
					yy_n = yypact[ *yy_ps ] + YYERRCODE;
					if ( yy_n >= 0 && yy_n < YYLAST &&
						yychk[yyact[yy_n]] == YYERRCODE)					{
						/*
						** simulate shift of "error"
						*/
						yy_state = yyact[ yy_n ];
						goto yy_stack;
					}
					/*
					** current state has no shift on
					** "error", pop stack
					*/
#if YYDEBUG
#	define _POP_ "Error recovery pops state %d, uncovers state %d\n"
					if ( yydebug )
						printf( _POP_, *yy_ps,
							yy_ps[-1] );
#	undef _POP_
#endif
					yy_ps--;
					yy_pv--;
				}
				/*
				** there is no state on stack with "error" as
				** a valid shift.  give up.
				*/
				YYABORT;
			case 3:		/* no shift yet; eat a token */
#if YYDEBUG
				/*
				** if debugging, look up token in list of
				** pairs.  0 and negative shouldn't occur,
				** but since timing doesn't matter when
				** debugging, it doesn't hurt to leave the
				** tests here.
				*/
				if ( yydebug )
				{
					register int yy_i;

					printf( "Error recovery discards " );
					if ( yychar == 0 )
						printf( "token end-of-file\n" );
					else if ( yychar < 0 )
						printf( "token -none-\n" );
					else
					{
						for ( yy_i = 0;
							yytoks[yy_i].t_val >= 0;
							yy_i++ )
						{
							if ( yytoks[yy_i].t_val
								== yychar )
							{
								break;
							}
						}
						printf( "token %s\n",
							yytoks[yy_i].t_name );
					}
				}
#endif /* YYDEBUG */
				if ( yychar == 0 )	/* reached EOF. quit */
					YYABORT;
				yychar = -1;
				goto yy_newstate;
			}
		}/* end if ( yy_n == 0 ) */
		/*
		** reduction by production yy_n
		** put stack tops, etc. so things right after switch
		*/
#if YYDEBUG
		/*
		** if debugging, print the string that is the user's
		** specification of the reduction which is just about
		** to be done.
		*/
		if ( yydebug )
			printf( "Reduce by (%d) \"%s\"\n",
				yy_n, yyreds[ yy_n ] );
#endif
		yytmp = yy_n;			/* value to switch over */
		yypvt = yy_pv;			/* $vars top of value stack */
		/*
		** Look in goto table for next state
		** Sorry about using yy_state here as temporary
		** register variable, but why not, if it works...
		** If yyr2[ yy_n ] doesn't have the low order bit
		** set, then there is no action to be done for
		** this reduction.  So, no saving & unsaving of
		** registers done.  The only difference between the
		** code just after the if and the body of the if is
		** the goto yy_stack in the body.  This way the test
		** can be made before the choice of what to do is needed.
		*/
		{
			/* length of production doubled with extra bit */
			register int yy_len = yyr2[ yy_n ];

			if ( !( yy_len & 01 ) )
			{
				yy_len >>= 1;
				yyval = ( yy_pv -= yy_len )[1];	/* $$ = $1 */
				yy_state = yypgo[ yy_n = yyr1[ yy_n ] ] +
					*( yy_ps -= yy_len ) + 1;
				if ( yy_state >= YYLAST ||
					yychk[ yy_state =
					yyact[ yy_state ] ] != -yy_n )
				{
					yy_state = yyact[ yypgo[ yy_n ] ];
				}
				goto yy_stack;
			}
			yy_len >>= 1;
			yyval = ( yy_pv -= yy_len )[1];	/* $$ = $1 */
			yy_state = yypgo[ yy_n = yyr1[ yy_n ] ] +
				*( yy_ps -= yy_len ) + 1;
			if ( yy_state >= YYLAST ||
				yychk[ yy_state = yyact[ yy_state ] ] != -yy_n )
			{
				yy_state = yyact[ yypgo[ yy_n ] ];
			}
		}
					/* save until reenter driver code */
		yystate = yy_state;
		yyps = yy_ps;
		yypv = yy_pv;
	}
	/*
	** code supplied by user is placed in this switch
	*/

		switch(yytmp){

case 1:	/* statement : expression */
# line 25 "calc.y"
{ 
	yyval.tr = new_ExprTree();
	yyval.tr->type = ETR_STATEMENT;
	
	yyval.tr->child[0] = yypvt[-0].tr;
	yyval.tr->nochild=1;
	root = yyval.tr;
	parentfy_ExprTree(root);
	} /*NOTREACHED*/ break;
case 2:	/* method : tag '(' commalist ')' */
# line 36 "calc.y"
{
	yyval.tr = new_ExprTree_method(yypvt[-3].tr,yypvt[-1].tr);
	} /*NOTREACHED*/ break;
case 3:	/* method : tag '(' ')' */
# line 39 "calc.y"
{
	yyval.tr = new_ExprTree_method(yypvt[-2].tr,NULL);
	} /*NOTREACHED*/ break;
case 4:	/* expression : expression '+' expression */
# line 44 "calc.y"
{ 
	yyval.tr = new_ExprTree_binary_expr(yypvt[-2].tr,'+',yypvt[-0].tr);
	} /*NOTREACHED*/ break;
case 5:	/* expression : expression '-' expression */
# line 47 "calc.y"
{
	yyval.tr = new_ExprTree_binary_expr(yypvt[-2].tr,'-',yypvt[-0].tr);
	} /*NOTREACHED*/ break;
case 6:	/* expression : expression '*' expression */
# line 50 "calc.y"
{
	yyval.tr = new_ExprTree_binary_expr(yypvt[-2].tr,'*',yypvt[-0].tr);
	} /*NOTREACHED*/ break;
case 7:	/* expression : expression '/' expression */
# line 53 "calc.y"
{
	yyval.tr = new_ExprTree_binary_expr(yypvt[-2].tr,'/',yypvt[-0].tr);
	} /*NOTREACHED*/ break;
case 8:	/* expression : '(' expression ')' */
# line 56 "calc.y"
{ yyval.tr = yypvt[-1].tr; } /*NOTREACHED*/ break;
case 9:	/* expression : NUMBER */
# line 57 "calc.y"
{ yyval.tr = yypvt[-0].tr; } /*NOTREACHED*/ break;
case 10:	/* expression : tag */
# line 58 "calc.y"
{ yyval.tr = yypvt[-0].tr; } /*NOTREACHED*/ break;
case 11:	/* expression : method */
# line 59 "calc.y"
{ yyval.tr = yypvt[-0].tr; } /*NOTREACHED*/ break;
case 12:	/* commalist : commalist ',' expression */
# line 62 "calc.y"
{ yyval.tr = add_to_commalist_ExprTree(yypvt[-2].tr,yypvt[-0].tr); } /*NOTREACHED*/ break;
case 13:	/* commalist : expression */
# line 63 "calc.y"
{ yyval.tr = new_ExprTree_commalist(yypvt[-0].tr); } /*NOTREACHED*/ break;
case 14:	/* tag : NAME */
# line 66 "calc.y"
{ yyval.tr = new_ExprTree_tag_from_name(yypvt[-0].tr); } /*NOTREACHED*/ break;
case 15:	/* tag : tag STRUCTREF tag */
# line 67 "calc.y"
{ yyval.tr = new_ExprTree_struct_ref(yypvt[-2].tr,yypvt[-1].tr,yypvt[-0].tr); } /*NOTREACHED*/ break;
case 16:	/* tag : tag '[' expression ']' */
# line 68 "calc.y"
{ yyval.tr = new_ExprTree_array(yypvt[-3].tr,yypvt[-1].tr); } /*NOTREACHED*/ break;
case 17:	/* tag : '&' tag */
# line 69 "calc.y"
{ yyval.tr = new_ExprTree_ref('*',yypvt[-0].tr); } /*NOTREACHED*/ break;
case 18:	/* tag : '*' tag */
# line 70 "calc.y"
{ yyval.tr = new_ExprTree_ref('*',yypvt[-0].tr); } /*NOTREACHED*/ break;
}


	goto yystack;		/* reset registers in driver code */
}
