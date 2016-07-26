#ifndef lint
static char yysccsid[] = "@(#)yaccpar	1.9 (Berkeley) 02/21/93";
#endif
#define YYBYACC 1
#define YYMAJOR 1
#define YYMINOR 9
#define yyclearin (yychar=(-1))
#define yyerrok (yyerrflag=0)
#define YYRECOVERING (yyerrflag!=0)
#define YYPREFIX "yy"
#line 3 "calc.y"
#include <stdio.h>
#include "exprtree.h"

extern ExprTree * root;
extern char * calc_lex_string;
extern int stringpos;

#line 12 "calc.y"
typedef union {
	double dval;
	struct ExprTree *  tr;
} YYSTYPE;
#line 25 "y.tab.c"
#define NAME 257
#define NUMBER 258
#define STRUCTREF 259
#define DEREFERENCE 260
#define YYERRCODE 256
short yylhs[] = {                                        -1,
    0,    3,    3,    1,    1,    1,    1,    1,    1,    1,
    1,    4,    4,    2,    2,    2,    2,    2,
};
short yylen[] = {                                         2,
    1,    4,    3,    3,    3,    3,    3,    3,    1,    1,
    1,    3,    1,    1,    3,    4,    2,    2,
};
short yydefred[] = {                                      0,
   14,    9,    0,    0,    0,    0,    0,    0,   11,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    8,
    0,    0,    6,    7,   15,    0,    3,    0,    0,   16,
    2,    0,    0,
};
short yydgoto[] = {                                       6,
    7,    8,    9,   29,
};
short yysindex[] = {                                    -33,
    0,    0,  -30,  -30,  -33,    0,   -2,  -29,    0,  -90,
  -90,   37,  -33,  -33,  -33,  -33,  -30,  -33,  -38,    0,
   -5,   -5,    0,    0,    0,  -19,    0,   -2,  -31,    0,
    0,  -33,   -2,
};
short yyrindex[] = {                                      0,
    0,    0,    0,    0,    0,    0,   19,   22,    0,    6,
   14,    0,    0,    0,    0,    0,    0,    0,    0,    0,
   27,   32,    0,    0,    0,    0,    0,  -24,    0,    0,
    0,    0,  -23,
};
short yygindex[] = {                                      0,
   20,   12,    0,    0,
};
#define YYTABLESIZE 230
short yytable[] = {                                       4,
   18,    5,   27,    3,    4,   18,    5,    4,    3,   31,
   19,    3,   32,   17,   10,   11,   13,   12,    1,   13,
   12,   10,   15,   14,   12,   13,    5,   16,   25,    0,
    0,    4,   21,   22,   23,   24,   15,   26,   28,   15,
   14,   16,   13,    0,   16,   18,   18,   18,   18,   18,
   18,   33,   18,   17,   17,   17,   17,   17,   17,    0,
   17,   18,   10,   10,   10,   10,   10,    5,   10,    5,
    5,    5,    4,   30,    4,    4,    4,   20,   15,   14,
    0,   13,    0,   16,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,   18,    0,
    0,    0,    0,    0,    0,    0,   17,    0,    0,    0,
    0,    0,    0,    0,   10,    0,    0,    0,    0,    5,
    0,    0,    0,    0,    4,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,   17,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    1,    2,
    0,    0,    0,    1,    2,    0,    1,    0,    0,   17,
};
short yycheck[] = {                                      38,
   91,   40,   41,   42,   38,    0,   40,   38,   42,   41,
   40,   42,   44,    0,    3,    4,   41,   41,    0,   44,
   44,    0,   42,   43,    5,   45,    0,   47,   17,   -1,
   -1,    0,   13,   14,   15,   16,   42,   18,   19,   42,
   43,   47,   45,   -1,   47,   40,   41,   42,   43,   44,
   45,   32,   47,   40,   41,   42,   43,   44,   45,   -1,
   47,   91,   41,   42,   43,   44,   45,   41,   47,   43,
   44,   45,   41,   93,   43,   44,   45,   41,   42,   43,
   -1,   45,   -1,   47,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   93,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   93,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   93,   -1,   -1,   -1,   -1,   93,
   -1,   -1,   -1,   -1,   93,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,  259,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,  257,  258,
   -1,   -1,   -1,  257,  258,   -1,  257,   -1,   -1,  259,
};
#define YYFINAL 6
#ifndef YYDEBUG
#define YYDEBUG 1
#endif
#define YYMAXTOKEN 260
#if YYDEBUG
char *yyname[] = {
"end-of-file",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,"'&'",0,"'('","')'","'*'","'+'","','","'-'",0,"'/'",0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"'['",0,"']'",
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,"NAME","NUMBER","STRUCTREF","DEREFERENCE",
};
char *yyrule[] = {
"$accept : statement",
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
#endif
#ifdef YYSTACKSIZE
#undef YYMAXDEPTH
#define YYMAXDEPTH YYSTACKSIZE
#else
#ifdef YYMAXDEPTH
#define YYSTACKSIZE YYMAXDEPTH
#else
#define YYSTACKSIZE 500
#define YYMAXDEPTH 500
#endif
#endif
int yydebug;
int yynerrs;
int yyerrflag;
int yychar;
short *yyssp;
YYSTYPE *yyvsp;
YYSTYPE yyval;
YYSTYPE yylval;
short yyss[YYSTACKSIZE];
YYSTYPE yyvs[YYSTACKSIZE];
#define yystacksize YYSTACKSIZE
#line 73 "calc.y"

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



#line 205 "y.tab.c"
#define YYABORT goto yyabort
#define YYREJECT goto yyabort
#define YYACCEPT goto yyaccept
#define YYERROR goto yyerrlab
int
yyparse()
{
    register int yym, yyn, yystate;
#if YYDEBUG
    register char *yys;
    extern char *getenv();

    if (yys = getenv("YYDEBUG"))
    {
        yyn = *yys;
        if (yyn >= '0' && yyn <= '9')
            yydebug = yyn - '0';
    }
#endif

    yynerrs = 0;
    yyerrflag = 0;
    yychar = (-1);

    yyssp = yyss;
    yyvsp = yyvs;
    *yyssp = yystate = 0;

yyloop:
    if (yyn = yydefred[yystate]) goto yyreduce;
    if (yychar < 0)
    {
        if ((yychar = yylex()) < 0) yychar = 0;
#if YYDEBUG
        if (yydebug)
        {
            yys = 0;
            if (yychar <= YYMAXTOKEN) yys = yyname[yychar];
            if (!yys) yys = "illegal-symbol";
            printf("%sdebug: state %d, reading %d (%s)\n",
                    YYPREFIX, yystate, yychar, yys);
        }
#endif
    }
    if ((yyn = yysindex[yystate]) && (yyn += yychar) >= 0 &&
            yyn <= YYTABLESIZE && yycheck[yyn] == yychar)
    {
#if YYDEBUG
        if (yydebug)
            printf("%sdebug: state %d, shifting to state %d\n",
                    YYPREFIX, yystate, yytable[yyn]);
#endif
        if (yyssp >= yyss + yystacksize - 1)
        {
            goto yyoverflow;
        }
        *++yyssp = yystate = yytable[yyn];
        *++yyvsp = yylval;
        yychar = (-1);
        if (yyerrflag > 0)  --yyerrflag;
        goto yyloop;
    }
    if ((yyn = yyrindex[yystate]) && (yyn += yychar) >= 0 &&
            yyn <= YYTABLESIZE && yycheck[yyn] == yychar)
    {
        yyn = yytable[yyn];
        goto yyreduce;
    }
    if (yyerrflag) goto yyinrecovery;
#ifdef lint
    goto yynewerror;
#endif
yynewerror:
    yyerror("syntax error");
#ifdef lint
    goto yyerrlab;
#endif
yyerrlab:
    ++yynerrs;
yyinrecovery:
    if (yyerrflag < 3)
    {
        yyerrflag = 3;
        for (;;)
        {
            if ((yyn = yysindex[*yyssp]) && (yyn += YYERRCODE) >= 0 &&
                    yyn <= YYTABLESIZE && yycheck[yyn] == YYERRCODE)
            {
#if YYDEBUG
                if (yydebug)
                    printf("%sdebug: state %d, error recovery shifting\
 to state %d\n", YYPREFIX, *yyssp, yytable[yyn]);
#endif
                if (yyssp >= yyss + yystacksize - 1)
                {
                    goto yyoverflow;
                }
                *++yyssp = yystate = yytable[yyn];
                *++yyvsp = yylval;
                goto yyloop;
            }
            else
            {
#if YYDEBUG
                if (yydebug)
                    printf("%sdebug: error recovery discarding state %d\n",
                            YYPREFIX, *yyssp);
#endif
                if (yyssp <= yyss) goto yyabort;
                --yyssp;
                --yyvsp;
            }
        }
    }
    else
    {
        if (yychar == 0) goto yyabort;
#if YYDEBUG
        if (yydebug)
        {
            yys = 0;
            if (yychar <= YYMAXTOKEN) yys = yyname[yychar];
            if (!yys) yys = "illegal-symbol";
            printf("%sdebug: state %d, error recovery discards token %d (%s)\n",
                    YYPREFIX, yystate, yychar, yys);
        }
#endif
        yychar = (-1);
        goto yyloop;
    }
yyreduce:
#if YYDEBUG
    if (yydebug)
        printf("%sdebug: state %d, reducing by rule %d (%s)\n",
                YYPREFIX, yystate, yyn, yyrule[yyn]);
#endif
    yym = yylen[yyn];
    yyval = yyvsp[1-yym];
    switch (yyn)
    {
case 1:
#line 25 "calc.y"
{ 
	yyval.tr = new_ExprTree();
	yyval.tr->type = ETR_STATEMENT;
	
	yyval.tr->child[0] = yyvsp[0].tr;
	yyval.tr->nochild=1;
	root = yyval.tr;
	parentfy_ExprTree(root);
	}
break;
case 2:
#line 36 "calc.y"
{
	yyval.tr = new_ExprTree_method(yyvsp[-3].tr,yyvsp[-1].tr);
	}
break;
case 3:
#line 39 "calc.y"
{
	yyval.tr = new_ExprTree_method(yyvsp[-2].tr,NULL);
	}
break;
case 4:
#line 44 "calc.y"
{ 
	yyval.tr = new_ExprTree_binary_expr(yyvsp[-2].tr,'+',yyvsp[0].tr);
	}
break;
case 5:
#line 47 "calc.y"
{
	yyval.tr = new_ExprTree_binary_expr(yyvsp[-2].tr,'-',yyvsp[0].tr);
	}
break;
case 6:
#line 50 "calc.y"
{
	yyval.tr = new_ExprTree_binary_expr(yyvsp[-2].tr,'*',yyvsp[0].tr);
	}
break;
case 7:
#line 53 "calc.y"
{
	yyval.tr = new_ExprTree_binary_expr(yyvsp[-2].tr,'/',yyvsp[0].tr);
	}
break;
case 8:
#line 56 "calc.y"
{ yyval.tr = yyvsp[-1].tr; }
break;
case 9:
#line 57 "calc.y"
{ yyval.tr = yyvsp[0].tr; }
break;
case 10:
#line 58 "calc.y"
{ yyval.tr = yyvsp[0].tr; }
break;
case 11:
#line 59 "calc.y"
{ yyval.tr = yyvsp[0].tr; }
break;
case 12:
#line 62 "calc.y"
{ yyval.tr = add_to_commalist_ExprTree(yyvsp[-2].tr,yyvsp[0].tr); }
break;
case 13:
#line 63 "calc.y"
{ yyval.tr = new_ExprTree_commalist(yyvsp[0].tr); }
break;
case 14:
#line 66 "calc.y"
{ yyval.tr = new_ExprTree_tag_from_name(yyvsp[0].tr); }
break;
case 15:
#line 67 "calc.y"
{ yyval.tr = new_ExprTree_struct_ref(yyvsp[-2].tr,yyvsp[-1].tr,yyvsp[0].tr); }
break;
case 16:
#line 68 "calc.y"
{ yyval.tr = new_ExprTree_array(yyvsp[-3].tr,yyvsp[-1].tr); }
break;
case 17:
#line 69 "calc.y"
{ yyval.tr = new_ExprTree_ref('*',yyvsp[0].tr); }
break;
case 18:
#line 70 "calc.y"
{ yyval.tr = new_ExprTree_ref('*',yyvsp[0].tr); }
break;
#line 438 "y.tab.c"
    }
    yyssp -= yym;
    yystate = *yyssp;
    yyvsp -= yym;
    yym = yylhs[yyn];
    if (yystate == 0 && yym == 0)
    {
#if YYDEBUG
        if (yydebug)
            printf("%sdebug: after reduction, shifting from state 0 to\
 state %d\n", YYPREFIX, YYFINAL);
#endif
        yystate = YYFINAL;
        *++yyssp = YYFINAL;
        *++yyvsp = yyval;
        if (yychar < 0)
        {
            if ((yychar = yylex()) < 0) yychar = 0;
#if YYDEBUG
            if (yydebug)
            {
                yys = 0;
                if (yychar <= YYMAXTOKEN) yys = yyname[yychar];
                if (!yys) yys = "illegal-symbol";
                printf("%sdebug: state %d, reading %d (%s)\n",
                        YYPREFIX, YYFINAL, yychar, yys);
            }
#endif
        }
        if (yychar == 0) goto yyaccept;
        goto yyloop;
    }
    if ((yyn = yygindex[yym]) && (yyn += yystate) >= 0 &&
            yyn <= YYTABLESIZE && yycheck[yyn] == yystate)
        yystate = yytable[yyn];
    else
        yystate = yydgoto[yym];
#if YYDEBUG
    if (yydebug)
        printf("%sdebug: after reduction, shifting from state %d \
to state %d\n", YYPREFIX, *yyssp, yystate);
#endif
    if (yyssp >= yyss + yystacksize - 1)
    {
        goto yyoverflow;
    }
    *++yyssp = yystate;
    *++yyvsp = yyval;
    goto yyloop;
yyoverflow:
    yyerror("yacc stack overflow");
yyabort:
    return (1);
yyaccept:
    return (0);
}
