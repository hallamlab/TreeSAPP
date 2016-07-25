
%{
#include <stdio.h>
#include "exprtree.h"

extern ExprTree * root;
extern char * calc_lex_string;
extern int stringpos;

%}

%union {
	double dval;
	struct ExprTree *  tr;
}

%token <tr> NAME NUMBER STRUCTREF
%type  <tr> expression statement tag method commalist
%left '-' '+'
%left '*' '/'
%left '&' DEREFERENCE
%left STRUCTREF '[' ']'
%%

statement : expression { 
	$$ = new_ExprTree();
	$$->type = ETR_STATEMENT;
	
	$$->child[0] = $1;
	$$->nochild=1;
	root = $$;
	parentfy_ExprTree(root);
	}
	;

method : tag '(' commalist ')' {
	$$ = new_ExprTree_method($1,$3);
	}
	| tag '(' ')' {
	$$ = new_ExprTree_method($1,NULL);
	}
	;

expression : expression '+'  expression { 
	$$ = new_ExprTree_binary_expr($1,'+',$3);
	}	
	| expression '-' expression {
	$$ = new_ExprTree_binary_expr($1,'-',$3);
	}
	| expression '*' expression {
	$$ = new_ExprTree_binary_expr($1,'*',$3);
	}	
	| expression '/' expression {
	$$ = new_ExprTree_binary_expr($1,'/',$3);
	}
	| '(' expression ')' { $$ = $2; }
	| NUMBER { $$ = $1; }
	| tag { $$ = $1; }
	| method { $$ = $1; }
	;

commalist : commalist ',' expression { $$ = add_to_commalist_ExprTree($1,$3); }
	| expression { $$ = new_ExprTree_commalist($1); }
	;

tag : NAME { $$ = new_ExprTree_tag_from_name($1); }
	| tag STRUCTREF tag { $$ = new_ExprTree_struct_ref($1,$2,$3); }
	| tag '[' expression ']' { $$ = new_ExprTree_array($1,$3); }
	| '&' tag { $$ = new_ExprTree_ref('*',$2); }
	| '*' tag %prec DEREFERENCE { $$ = new_ExprTree_ref('*',$2); }
	;	
%%

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



