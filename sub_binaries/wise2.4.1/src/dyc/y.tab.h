
typedef union  {
	double dval;
	struct ExprTree *  tr;
} YYSTYPE;
extern YYSTYPE yylval;
# define NAME 257
# define NUMBER 258
# define STRUCTREF 259
# define DEREFERENCE 260
