
.SUFFIXES :  .dy



WISECOBJ = dyna2.o dynfile.o wisec.o dynafunc.o module.o type.o method.o dynadb.o friend.o \
	 inputfile.o variable.o modulefunc.o api.o display.o dynashadow.o labelmaster.o ftext.o funcinfo.o \
	objectinfo.o exprtree.o compugen.o docugen.o input.o dpimpl.o dbthread.o probal.o telegraph.o dynadebug.o dyshatter.o




test : calc.y calc.l exprtree.o type.o
	yacc -d -t calc.y
	lex calc.l
	cc $(CFLAGS) y.tab.c
	cc $(CFLAGS) lex.yy.c
	cc -o test y.tab.o lex.yy.o type.o exprtree.o -ll -ly -lwisebase $(LFLAGS)

# module.o sometime!
# frontend.o taken out

# CFLAGS = -c -g -DUNIX -I ~/local/prog/hmmer-1.8.1/


CC = cc
CFLAGS = -pthread -c -g  -DUNIX  -I../base/
#CFLAGS = -c -g  -DUNIX -DNOERROR -DLINUX -I../base/
#CFLAGS = -c -DNOERROR
INCFLAGS = -I../base/
LFLAGS =  -g  -lm -L../base/ -lpthread
AR_OPTIONS = ru

RANLIB_NEEDED = 0

#
# Base library build - call wisebase
#

linux : dyc.o $(WISECOBJ)
	flex calc.l
	$(CC) -DLINUX -c -I../base/  y.tab.c
	$(CC) -DLINUX -c -I../base/  lex.yy.c
	$(CC) -o dyc dyc.o $(WISECOBJ) y.tab.o lex.yy.o  -lwisebase $(LFLAGS)

dyc  : dyc.o $(WISECOBJ) y.tab.o lex.yy.o 
	$(CC) -o dyc dyc.o $(WISECOBJ) y.tab.o lex.yy.o  -ll -lwisebase $(LFLAGS)

yacc :
	yacc -d -t calc.y
	lex calc.l

#
# Suffix rules.
#

DFLAGS = -l

.c.o :
	$(CC) $(CFLAGS) $(INCFLAGS) $?

.dy.o :
	dyc $(DFLAGS) $<	
	$(CC) $(CFLAGS) $*.c

.dy.c : 
	dyc $(DFLAGS) $<	


init : 
	dyc $(DFLAGS) *.dy

#
# pro2DNAmatrix.c pro2DNAmatrix.h : pro2DNAmatrix.dy
#	dy pro2DNAmatrix.dy







