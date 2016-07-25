
%{
#include "y.tab.h"
#include "exprtree.h"

#ifdef LINUX
#undef YY_INPUT
#define YY_INPUT(b, r, ms) (r = my_yyinput(b, ms))
#undef yywrap
#endif

extern char * calc_lex_string;
extern int stringpos;

#ifndef LINUX
#undef  input
#undef  unput
#endif

%}



%%
[ \t] ;
[0-9]+ { yylval.tr = new_ExprTree(); 
	yylval.tr->type = ETR_NUMBER;
	yylval.tr->word = stringalloc((char * )yytext);
	/*	printf("Got a number %s\n",yytext); */

	return NUMBER; }
[A-Za-z][A-Za-z0-9_]* { yylval.tr = new_ExprTree();
	yylval.tr->type = ETR_NAME;
	yylval.tr->word = stringalloc((char * )yytext);
	/*	printf("Got a name %s \n",yytext);*/
	return NAME; }

\->  { yylval.tr = new_ExprTree();
	yylval.tr->word = stringalloc((char * )yytext);
	yylval.tr->type = ETR_OPERATOR;
	return STRUCTREF;
	}
\.  { yylval.tr = new_ExprTree();
	yylval.tr->word = stringalloc((char * )yytext);
	yylval.tr->type = ETR_OPERATOR;
	return STRUCTREF;
	}

.  { 
  /* printf("operator\n"); */
  return yytext[0]; }

\n return 0;
%%

#ifdef LINUX
int my_yyinput(char * buf,int max_size)
{
  int i;

  /* fprintf(stderr,"[%s] Got into yyinput, calling %d and with %d in string\n",calc_lex_string,max_size,stringpos); */
  
  if( calc_lex_string[stringpos] == '\0') {
    /*     fprintf(stderr,"Going to return 0\n"); */
     buf[0] = EOF;
     return 0;
  }

  for(i=0;i< max_size && calc_lex_string[stringpos] != '\0';i++,stringpos++) {
    buf[i] = calc_lex_string[stringpos];
  }
  buf[i] = EOF; 
  /*  fprintf(stderr,"Returning %d\n",i); */

  
  return i;
}

int yywrap(void) 
{
  return 1;
}

#else

char input(void) {
  return calc_lex_string[stringpos++];
}

void unput(char c) {
  calc_lex_string[--stringpos] = c;
}
#endif

