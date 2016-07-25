

.SUFFIXES :  .dy


NEWG = geneparser21.o \
	geneparameter.o \
	genestats.o \
	geneutil.o \
	threestatemodel.o \
	genefrequency.o \
	splicesitemodeler.o \
	genewise4.o \
	genewise6.o \
	genewise21.o \
	geneloop21.o \
	geneloop6.o \
        gwlite.o \
        gwlitemodel.o \
	gwrap.o \
	matchsum.o \
	estwrap.o \
	genewisemodel.o \
	cdparser.o \
	genedisplay.o \
	estwise3.o \
	estslim3.o \
	estloop3.o \
	estfrag3.o \
	estslimloop.o \
	gwquickdb.o \
	threestatedb.o \
	pfamhmmer1db.o \
	pwmdna.o \
	wise2xhmmer2.o \
	genewisemodeldb.o \
	seqhit.o \
	geneparser4.o \
	sw_wrap.o \
	abc.o \
	pba.o \
	seqaligndisplay.o \
	dbadisplay.o \
	proteinsw.o \
	estquick3.o 

all : dnal psw pswdb dba estwise genewise genewisedb estwisedb genomewise libmodel.a

libmodel.a : $(NEWG) 
	ar $(AR_OPTIONS) libmodel.a $(NEWG)


wise2xhmmer2.o : wise2xhmmer2.c
	$(CC) $(CFLAGS) -I$(HMMER_INCLUDE) -D$(HMMER_DEFINE) $(INCFLAGS) wise2xhmmer2.c

#	$(CC) $(CFLAGS) -I../HMMer2/ $(INCFLAGS) wise2xhmmer2.c

testgenestat : testgenestat.o genestats.o pwmdna.o 
	$(CC) -o testgenestat testgenestat.o genestats.o pwmdna.o -ldyna -lwisebase $(LFLAGS) -lm

makepwm : makepwm.o pwmdna.o 
	$(CC) -o makepwm makepwm.o pwmdna.o -ldyna -lwisebase $(LFLAGS) -lm

testwise2xhmmer2.o : testwise2xhmmer2.c
	$(CC) $(CFLAGS) -I../HMMer2/ $(INCFLAGS) testwise2xhmmer2.c

testwise2xhmmer2 : testwise2xhmmer2.o threestatemodel.o threestatedb.o pfamhmmer1db.o 
	$(CC) -o testwise2xhmmer2 wise2xhmmer2.o testwise2xhmmer2.o threestatemodel.o threestatedb.o pfamhmmer1db.o -ldyna -lwisebase -lhmmer $(LFLAGS) -lm

plan7test.o : plan7test.c
	$(CC) $(CFLAGS) -I../HMMer2/ plan7test.c

plan7test : plan7test.o plan7.o wise2xhmmer2.o threestatemodel.o threestatedb.o  pfamhmmer1db.o
	$(CC) -o plan7test plan7test.o plan7.o threestatemodel.o threestatedb.o pfamhmmer1db.o wise2xhmmer2.o -ldyna -lwisebase -lhmmer $(LFLAGS) -lm

seqstat : seqstat.o
	$(CC) -o seqstat seqstat.o -ldyna -lwisebase $(LFLAGS)

testpwm : testpwm.o pwmdna.o
	$(CC) -o testpwm testpwm.o pwmdna.o -ldyna -lwisebase $(LFLAGS) -lm -lefence


threestate : $(MODELOBJ)

gw21 : $(MODELOBJ) $(GENEWISEMODEL) genewisemodel.o

pwise : $(MODELOBJ) threestatemat.o threestatedisplay.o

test : test.o proteinsw.o protprotdis.o
	$(CC) -o test test.o proteinsw.o protprotdis.o  -ldyna -lwisebase ../dynlibsrc/matrix.o $(LFLAGS)

testanc : testanc.o proteinsw.o abc.o pba.o sw_wrap.o
	$(CC) -o testanc testanc.o proteinsw.o abc.o pba.o sw_wrap.o $(LFLAGS) -lhmmer -ldyna -lwisebase -lm $(EXTRALIBS)

cdtest : cdtest.o cdna2genomic.o d2d_display.o
	$(CC) -o cdtest cdtest.o cdna2genomic.o d2d_display.o  -ldyna -lwisebase $(LFLAGS)

prob2bits : prob2bits.o
	$(CC) -o prob2bits  $(LFLAGS) -ldyna -lwisebase -lm

quick :
	cc -O -c -I../dynlibsrc/ -I../base/ genewise6.c
	cc -O -c -I../dynlibsrc/ -I../base/ geneloop6.c
	cc -O -c -I../dynlibsrc/ -I../base/ genewise21.c
	cc -O -c -I../dynlibsrc/ -I../base/ geneloop21.c
	cc -O -c -I../dynlibsrc/ -I../base/ estwise3.c
	cc -O -c -I../dynlibsrc/ -I../base/ estloop3.c

quicke :
	cc -O -c -I../dynlibsrc/ -I../base/ estwise3.c
	cc -O -c -I../dynlibsrc/ -I../base/ estloop3.c
	cc -O -c -I../dynlibsrc/ -I../base/ estslim3.c

quicksw :
	cc -O -c -I../dynlibsrc/ -I../base/ proteinsw.c
	cc -O -c -I../dynlibsrc/ -I../base/ abc.c

genewisedb.o : genewisedb.c
	$(CC) $(CFLAGS) genewisedb.c -I../base/ -I../dynlibsrc/ -I../HMMer2/ 

estwisedb.o : estwisedb.c
	$(CC) $(CFLAGS) estwisedb.c -I../base/ -I../dynlibsrc/ -I../HMMer2/

estwise.o : estwise.c
	$(CC) $(CFLAGS) estwise.c -I../base/ -I../dynlibsrc/ -I../HMMer2/ -I.

genewise.o : genewise.c
	$(CC) $(CFLAGS) genewise.c -I../base/ -I../dynlibsrc/ -I../HMMer2/ -I.

oldpostwise : oldpostwise.o $(NEWG) 
	$(CC) -g -o oldpostwise oldpostwise.o $(NEWG) $(LFLAGS) -lhmmer -ldyna -lwisebase -lm 

testdna : testdna.o dnaalign.o seqaligndisplay.o 
	$(CC) -o testdna testdna.o dnaalign.o seqaligndisplay.o $(LFLAGS) -ldyna -lwisebase -lm $(EXTRALIBS)

dnal : dnal.o dnaalign.o seqaligndisplay.o 
	$(CC) -o dnal dnal.o dnaalign.o seqaligndisplay.o $(LFLAGS) -ldyna -lwisebase -lm $(EXTRALIBS)

psw : psw.o seqaligndisplay.o proteinsw.o sw_wrap.o abc.o pba.o
	$(CC) -o psw psw.o sw_wrap.o seqaligndisplay.o proteinsw.o abc.o pba.o   $(LFLAGS) -ldyna -lwisebase -lm $(EXTRALIBS)

pswdb : pswdb.o seqaligndisplay.o proteinsw.o sw_wrap.o abc.o pba.o 
	$(CC) -o pswdb pswdb.o sw_wrap.o seqaligndisplay.o proteinsw.o abc.o pba.o $(LFLAGS) -ldyna -lwisebase -lm $(EXTRALIBS)

efpswdb : pswdb.o seqaligndisplay.o proteinsw.o sw_wrap.o abc.o pba.o 
	$(CC) -o pswdb pswdb.o sw_wrap.o seqaligndisplay.o proteinsw.o abc.o pba.o $(LFLAGS) -ldyna -lwisebase -lm -lefence $(EXTRALIBS)

dba : dbac.o dba.o dbadisplay.o
	$(CC) -o dba dbac.o dba.o dbadisplay.o $(LFLAGS) -ldyna -lwisebase -lm $(EXTRALIBS)

genewise : genewise.o $(NEWG) 
	$(CC) -g -o genewise genewise.o $(NEWG) $(LFLAGS) -lhmmer -ldyna -lwisebase -lm $(EXTRALIBS)

cdnawise : cdnawise.o $(NEWG) cdnawise10.o 
	$(CC) -g -o cdnawise cdnawise.o cdnawise10.o $(NEWG) $(LFLAGS) -lhmmer -ldyna -lwisebase -lm $(EXTRALIBS)

efgw   : genewise.o $(NEWG) 
	$(CC) -g -o genewise genewise.o  $(NEWG)  $(LFLAGS) -lhmmer -ldyna -lwisebase -lm $(EXTRALIBS) -lefence


genewisedb : genewisedb.o $(NEWG) 
	$(CC) -g -o genewisedb genewisedb.o  $(NEWG)  $(LFLAGS) -lhmmer -ldyna -lwisebase -lm $(EXTRALIBS) 

estwisedb : estwisedb.o $(NEWG) 
	$(CC) -g -o estwisedb estwisedb.o $(NEWG)  $(LFLAGS) -lhmmer -ldyna -lwisebase -lm $(EXTRALIBS)

edb : estwisedb.o $(NEWG) 
	$(CC) -g -o estwisedb estwisedb.o $(NEWG)  $(LFLAGS) -lhmmer -ldyna -lwisebase -lm $(EXTRALIBS) -lefence

estwise : estwise.o $(NEWG) 
	$(CC) -g -o estwise estwise.o  $(NEWG) $(LFLAGS) -lhmmer -ldyna -lwisebase -lm $(EXTRALIBS)

genomewise : genomewise.o genomewise9.o genome_evidence.o est_evidence.o geneutil.o
	$(CC) -g -o genomewise genomewise.o genomewise9.o genome_evidence.o est_evidence.o geneutil.o $(LFLAGS) -ldyna -lwisebase -lm $(EXTRALIBS)

clonewise : clonewise.o clonewisedp.o localclonewisedp.o mapstruct.o
	$(CC) -g -o clonewise clonewise.o clonewisedp.o localclonewisedp.o mapstruct.o $(LFLAGS) -ldyna -lwisebase -lm $(EXTRALIBS)


HMMER_LIBS = ../HMMer2/
HMMER_INCLUDE = ../HMMer2/
HMMER_DEFINE = 


CC = cc
#CFLAGS = -c -O  -I../base/ -I../dynlibsrc/
CFLAGS = -c -g -DUNIX  -I../base/ -I../dynlibsrc/ -I.
#CFLAGS = -c -g -O2 -Wall -DUNIX -I../base/ -I../dynlibsrc/
#LFLAGS =  -L../base/ -L../dynlibsrc/ -L../HMMer2/  

LFLAGS =  -g -L../base/ -L../dynlibsrc/ -L$(HMMER_LIBS)
AR_OPTIONS = ru

RANLIB_NEEDED = 0

#
# Base library build - call wisebase
#


#test : test.o libwisedb.a 
#	$(CC) -o test  test.o libwisedb.a  -lprob -ldyn  -lwisebase $(LFLAGS)


#
#
#

.c.o :
	$(CC) $(CFLAGS) $(INCFLAGS) $?

INCFLAGS = -I../base/ -I../dynlibsrc/
DFLAGS = -l -D -n Wise2_ -a _api.h -b _api.t -latex -perl -pthreads -dbtrace 5 -nocwarn
#DFLAGS = -l -O -D -n Wise2_ -a _api.h -b _api.t -latex -perl -pthreads -dbtrace 5
#DFLAGS = -l -n Wise2_ 

.dy.c : 
	dyc  $(DFLAGS) $<	

.dy.o :
	dyc $(DFLAGS) $<	
	$(CC) $(CFLAGS) $*.c

init : 
	dyc $(DFLAGS) *.dy

#
# pro2DNAmatrix.c pro2DNAmatrix.h : pro2DNAmatrix.dy
#	dy pro2DNAmatrix.dy







