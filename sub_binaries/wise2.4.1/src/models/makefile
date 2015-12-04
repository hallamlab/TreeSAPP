

.SUFFIXES :  .dy


NEWG = geneparser21.o \
	geneparameter.o \
	genestats.o \
	genewisehsp.o \
	geneutil.o \
	geneoutput.o \
	threestatemodel.o \
	genefrequency.o \
	splicesitemodeler.o \
	genewise4.o \
	genewise6.o \
	genestretch6.o \
	genewise21.o \
	geneloop21.o \
	geneloop6.o \
	genephase6.o \
        gwlite.o \
        gwlitemodel.o \
	gwrap.o \
	matchsum.o \
	estwrap.o \
	genewisemodel.o \
	phasemodel.o \
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
	standardout.o \
	geneparser4.o \
	sw_wrap.o \
	abc.o \
	pba.o \
	seqaligndisplay.o \
	dbadisplay.o \
	proteinsw.o \
	estquick3.o 

all : dnal psw pswdb dba estwise genewise genewisedb estwisedb genomewise sywise pseudowise promoterwise scanwise libmodel.a

libmodel.a : $(NEWG) 
	ar $(AR_OPTIONS) libmodel.a $(NEWG)


wise2xhmmer2.o : wise2xhmmer2.c
	$(CC) $(CFLAGS) -I$(HMMER_INCLUDE) -D$(HMMER_DEFINE) $(INCFLAGS) wise2xhmmer2.c

#	$(CC) $(CFLAGS) -I../HMMer2/ $(INCFLAGS) wise2xhmmer2.c

testgenestat : testgenestat.o genestats.o pwmdna.o 
	$(CC) -o testgenestat testgenestat.o genestats.o pwmdna.o -ldyna -lwisebase $(LFLAGS) 

test_aligng : test_aligng.o aligngenemodel.o pwmdna.o 
	$(CC) -o test_aligng test_aligng.o aligngenemodel.o pwmdna.o  -ldyna -lwisebase $(LFLAGS) 

makepwm : makepwm.o pwmdna.o 
	$(CC) -o makepwm makepwm.o pwmdna.o -ldyna -lwisebase $(LFLAGS) 

testwise2xhmmer2.o : testwise2xhmmer2.c
	$(CC) $(CFLAGS) -I../HMMer2/ $(INCFLAGS) testwise2xhmmer2.c

testwise2xhmmer2 : testwise2xhmmer2.o threestatemodel.o threestatedb.o pfamhmmer1db.o 
	$(CC) -o testwise2xhmmer2 wise2xhmmer2.o testwise2xhmmer2.o threestatemodel.o threestatedb.o pfamhmmer1db.o -ldyna -lwisebase -lhmmer $(LFLAGS) 

plan7test.o : plan7test.c
	$(CC) $(CFLAGS) -I../HMMer2/ plan7test.c

plan7test : plan7test.o plan7.o wise2xhmmer2.o threestatemodel.o threestatedb.o  pfamhmmer1db.o
	$(CC) -o plan7test plan7test.o plan7.o threestatemodel.o threestatedb.o pfamhmmer1db.o wise2xhmmer2.o -ldyna -lwisebase -lhmmer $(LFLAGS) 

seqstat : seqstat.o
	$(CC) -o seqstat seqstat.o -ldyna -lwisebase $(LFLAGS)

mutatedna : mutatedna.o
	$(CC) -o mutatedna mutatedna.o -ldyna -lwisebase $(LFLAGS) 

revcomp : revcomp.o
	$(CC) -o revcomp revcomp.o -ldyna -lwisebase $(LFLAGS) 

stupid_motif_scan : stupid_motif_scan.o
	$(CC) -o stupid_motif_scan stupid_motif_scan.o -ldyna -lwisebase $(LFLAGS) 

seqtrunc : seqtrunc.o
	$(CC) -o seqtrunc seqtrunc.o -ldyna -lwisebase $(LFLAGS) 

seqlist : seqlist.o
	$(CC) -o seqlist seqlist.o -ldyna -lwisebase $(LFLAGS) 

testpwm : testpwm.o pwmdna.o
	$(CC) -o testpwm testpwm.o pwmdna.o -ldyna -lwisebase $(LFLAGS)  -lefence


threestate : $(MODELOBJ)

gw21 : $(MODELOBJ) $(GENEWISEMODEL) genewisemodel.o

pwise : $(MODELOBJ) threestatemat.o threestatedisplay.o

test : test.o proteinsw.o protprotdis.o
	$(CC) -o test test.o proteinsw.o protprotdis.o  -ldyna_glib -ldyna -lwisebase ../dynlibsrc/matrix.o $(LFLAGS)

testanc : testanc.o proteinsw.o abc.o pba.o sw_wrap.o
	$(CC) -o testanc testanc.o proteinsw.o abc.o pba.o sw_wrap.o $(LFLAGS) -lhmmer -ldyna_glib -ldyna -lwisebase  $(EXTRALIBS)

test_tf : test_transfactor.o transfactor.o pwmdna.o
	$(CC) -o test_tf test_transfactor.o transfactor.o pwmdna.o $(LFLAGS) -ldyna -lwisebase 

motifwise : motifwise.o transfactor.o transregion.o pwmdna.o transregiondp.o 
	$(CC) -o motifwise motifwise.o transfactor.o transregion.o pwmdna.o transregiondp.o $(LFLAGS) -ldyna -lwisebase 

motifcluster : motifcluster.o transfactor.o transregion.o pwmdna.o transregiondp.o 
	$(CC) -o motifcluster motifcluster.o transfactor.o transregion.o pwmdna.o transregiondp.o $(LFLAGS) -ldyna -lwisebase 

motifdiff : motifdiff.o transfactor.o transregion.o pwmdna.o transregiondp.o 
	$(CC) -o motifdiff motifdiff.o transfactor.o transregion.o pwmdna.o transregiondp.o $(LFLAGS) -ldyna -lwisebase 


cdtest : cdtest.o cdna2genomic.o d2d_display.o
	$(CC) -o cdtest cdtest.o cdna2genomic.o d2d_display.o  -ldyna_glib -ldyna -lwisebase $(LFLAGS)

prob2bits : prob2bits.o
	$(CC) -o prob2bits  $(LFLAGS) -ldyna_glib -ldyna -lwisebase 

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

evopairwise.o : evopairwise.c
	$(CC) $(CFLAGS) evopairwise.c -I../base/ -I../dynlibsrc/ -I../HMMer2/ -I.

test_threestatehsp.o : test_threestatehsp.c
	$(CC) $(CFLAGS) test_threestatehsp.c -I../base/ -I../dynlibsrc/ -I../HMMer2/ -I.

fivestar.o : fivestar.c
	$(CC) $(CFLAGS) fivestar.c -I../base/ -I../dynlibsrc/ -I../HMMer2/ -I.

fivestarsearch.o : fivestarsearch.c
	$(CC) $(CFLAGS) fivestarsearch.c -I../base/ -I../dynlibsrc/ -I../HMMer2/ -I.

fivestarscan.o : fivestarscan.c
	$(CC) $(CFLAGS) fivestarscan.c -I../base/ -I../dynlibsrc/ -I../HMMer2/ -I.

oldpostwise : oldpostwise.o $(NEWG) 
	$(CC) -g -o oldpostwise oldpostwise.o $(NEWG) $(LFLAGS) -lhmmer -ldyna_glib -ldyna -lwisebase  

testdna : testdna.o dnaalign.o seqaligndisplay.o 
	$(CC) -o testdna testdna.o dnaalign.o seqaligndisplay.o $(LFLAGS) -ldyna_glib -ldyna -lwisebase  $(EXTRALIBS)

dnal : dnal.o dnaalign.o seqaligndisplay.o 
	$(CC) -o dnal dnal.o dnaalign.o seqaligndisplay.o $(LFLAGS) -ldyna_glib -ldyna -lwisebase  $(EXTRALIBS)

psw : psw.o seqaligndisplay.o proteinsw.o sw_wrap.o abc.o pba.o
	$(CC) -o psw psw.o sw_wrap.o seqaligndisplay.o proteinsw.o abc.o pba.o   $(LFLAGS) -ldyna_glib -ldyna -lwisebase  $(EXTRALIBS)

scanwisep.o : scanwisep.c
	$(CC) $(CFLAGS) $(INCFLAGS) -I../external/mott `glib-config --cflags` $?

hsp2aln_sw.o : hsp2aln_sw.c
	$(CC) $(CFLAGS) $(INCFLAGS) `glib-config --cflags` $?

scanwisep_mysql.o : scanwisep.c
	$(CC)  $(CFLAGS) -DSCAN_MYSQL -I../corba -I../external/mott -I../mysql/protein_index/ -I/usr/local/mysql/include scanwisep.c 

scanwisep_mysql : scanwisep_mysql.o seqaligndisplay.o proteinsw.o sw_wrap.o abc.o pba.o hsp2aln_sw.o
	$(CC) -o scanwisep scanwisep.o sw_wrap.o seqaligndisplay.o proteinsw.o abc.o pba.o hsp2aln_sw.o ../mysql/protein_index/mysql_protein_index.o  $(LFLAGS) -L../external/mott -L/usr/local/mysql/lib -ldyna_glib -ldyna ott -lwisebase -lmysqlclient -lm $(EXTRALIBS) -lpthread -lZ


scanwisep_corba.o : scanwisep.c
	$(CC)  $(CFLAGS) -DSCAN_CORBA -I../corba -I../external/mott `orbit-config --cflags server` scanwisep.c 

scanwisep_corba : scanwisep_corba.o seqaligndisplay.o proteinsw.o sw_wrap.o abc.o pba.o hsp2aln_sw.o
	$(CC) -o scanwisep scanwisep.o sw_wrap.o seqaligndisplay.o proteinsw.o abc.o pba.o hsp2aln_sw.o   ../corba/hspscan_corba-common.o ../corba/hspscan_corba-stubs.o   ../corba/corba_singleton.o ../corba/hspscan_corba_wrapper.o $(LFLAGS) -L../external/mott -ldyna_glib -ldyna -lmott -lwisebase  $(EXTRALIBS) `orbit-config --libs server` -lpthread 


scanwisep_wiseserver.o : scanwisep.c
	$(CC)  $(CFLAGS) $(INCFLAGS)  -o scanwisep_wiseserver.o -DSCAN_WISESERVER -I../network -I../socket -I../external/mott scanwisep.c 

scanwise : scanwisep_wiseserver.o seqaligndisplay.o proteinsw.o sw_wrap.o abc.o pba.o hsp2aln_sw.o
	$(CC) -o scanwise scanwisep_wiseserver.o sw_wrap.o seqaligndisplay.o proteinsw.o abc.o pba.o hsp2aln_sw.o   ../network/net_hspscan.o ../network/client_multihspscan.o  $(LFLAGS) -L../external/mott -L../socket -lmott -ldyna_glib -ldyna -lwisesocket -lwisebase  $(EXTRALIBS) -lpthread 

scanwisep_compress.o : scanwisep.c
	$(CC)  $(CFLAGS) -DSCAN_COMPRESS -I../dnaindex -I../network -I../socket -I../external/mott scanwisep.c 

scanwisep_compress : scanwisep_compress.o seqaligndisplay.o proteinsw.o sw_wrap.o abc.o pba.o hsp2aln_sw.o
	$(CC) -o scanwisep scanwisep.o sw_wrap.o seqaligndisplay.o proteinsw.o abc.o pba.o hsp2aln_sw.o   ../dnaindex/compressed_protein_index.o ../dnaindex/singleseqspace.o  ../dnaindex/kmer_direct.o  ../dnaindex/kmer_index_interface.o  $(LFLAGS) -L../external/mott -L../socket -lmott -ldyna_glib -ldyna -lwisesocket -lwisebase  $(EXTRALIBS) -lpthread 



scanwisep : scanwisep.o seqaligndisplay.o proteinsw.o sw_wrap.o abc.o pba.o hsp2aln_sw.o
	$(CC) -o scanwisep scanwisep.o sw_wrap.o seqaligndisplay.o proteinsw.o abc.o pba.o hsp2aln_sw.o   $(LFLAGS)  -L../external/mott -lmott -ldyna_glib -ldyna  -lwisebase  $(EXTRALIBS) `glib-config --libs` -lpthread

test_gwhsp : test_gwhsp.o genewisehsp.o
	$(CC) -o test_gwhsp test_gwhsp.o genewisehsp.o $(LFLAGS) -ldyna_glib -ldyna_glib -ldyna_glib -ldyna -lwisebase  $(EXTRALIBS) `glib-config --libs` -lpthread

localcishit.o : localcishit.c
	$(CC) $(CFLAGS) $(INCFLAGS) `glib-config --cflags` $?

promoterwise.o : promoterwise.c  
	$(CC) $(CFLAGS) $(INCFLAGS) `glib-config --cflags` $?

promoterwise : promoterwise.o localdba.o localcishit.o localcispara.o dbadisplay.o motifmatrix.o motifmatrixdp.o transfactor.o pwmdna.o pairwiseshortdna.o
	$(CC) -g -o promoterwise promoterwise.o localdba.o  localcishit.o localcispara.o dbadisplay.o  motifmatrix.o motifmatrixdp.o transfactor.o pwmdna.o  pairwiseshortdna.o $(LFLAGS) -ldyna_glib -ldyna -lwisebase  $(EXTRALIBS) `glib-config --libs` -lpthread

test_deng : test_deng.o localdba.o localcishit.o localcispara.o dbadisplay.o motifmatrix.o motifmatrixdp.o transfactor.o pwmdna.o pairwiseshortdna.o dnaprofile.o dnaprofiledp.o dnaprofileengine.o
	$(CC) -o test_deng test_deng.o  localdba.o localcishit.o localcispara.o dbadisplay.o motifmatrix.o motifmatrixdp.o transfactor.o pwmdna.o pairwiseshortdna.o dnaprofile.o dnaprofiledp.o dnaprofileengine.o $(LFLAGS) -ldyna_glib -ldyna -lwisebase  $(EXTRALIBS) `glib-config --libs` -lpthread

promotercluster : promotercluster.o localdba.o localcishit.o localcispara.o dbadisplay.o motifmatrix.o motifmatrixdp.o transfactor.o pwmdna.o pairwiseshortdna.o dnaprofile.o dnaprofiledp.o dnaprofileengine.o
	$(CC) -o promotercluster promotercluster.o  localdba.o localcishit.o localcispara.o dbadisplay.o motifmatrix.o motifmatrixdp.o transfactor.o pwmdna.o pairwiseshortdna.o dnaprofile.o dnaprofiledp.o dnaprofileengine.o $(LFLAGS) -ldyna_glib -ldyna -lwisebase  $(EXTRALIBS) `glib-config --libs` -lpthread

dnawise : dnawise.o dnahmm.o dnahmmdp.o seqaligndisplay.o 
	$(CC) -o dnawise dnawise.o dnahmm.o dnahmmdp.o  seqaligndisplay.o $(LFLAGS) -ldyna_glib -ldyna -lwisebase  -lpthread

pswdb : pswdb.o seqaligndisplay.o proteinsw.o sw_wrap.o abc.o pba.o 
	$(CC) -o pswdb pswdb.o sw_wrap.o seqaligndisplay.o proteinsw.o abc.o pba.o $(LFLAGS) -ldyna_glib -ldyna -lwisebase  $(EXTRALIBS)

efpswdb : pswdb.o seqaligndisplay.o proteinsw.o sw_wrap.o abc.o pba.o 
	$(CC) -o pswdb pswdb.o sw_wrap.o seqaligndisplay.o proteinsw.o abc.o pba.o $(LFLAGS) -ldyna_glib -ldyna -lwisebase  -lefence $(EXTRALIBS)

dba : dbac.o dba.o slimdba.o bigdba.o dbadisplay.o
	$(CC) -o dba dbac.o dba.o slimdba.o bigdba.o dbadisplay.o $(LFLAGS) -ldyna_glib -ldyna -lwisebase  $(EXTRALIBS)

lba : lba.o largeblockdp.o seqaligndisplay.o 
	$(CC) -o lba lba.o largeblockdp.o seqaligndisplay.o $(LFLAGS) -ldyna_glib -ldyna -lwisebase  $(EXTRALIBS)


fastdba : fastdba.o dba.o slimdba.o bigdba.o dbadisplay.o
	$(CC) -o fastdba fastdba.o dba.o slimdba.o bigdba.o dbadisplay.o $(LFLAGS) -ldyna_glib -ldyna -lwisebase  $(EXTRALIBS)

evopairwise : evopairwise.o threestatemodel.o threestateloop.o threestatedb.o pfamhmmer1db.o wise2xhmmer2.o seqaligndisplay.o 
	$(CC) -o evopairwise evopairwise.o threestatemodel.o threestateloop.o threestatedb.o pfamhmmer1db.o  wise2xhmmer2.o seqaligndisplay.o $(LFLAGS) -lhmmer -ldyna_glib -ldyna -lwisebase  $(EXTRALIBS)

test_threestatehsp : test_threestatehsp.o threestatemodel.o threestateloop.o threestatedb.o pfamhmmer1db.o wise2xhmmer2.o seqaligndisplay.o threestatehsp.o 
	$(CC) -o test_threestatehsp test_threestatehsp.o threestatemodel.o threestateloop.o threestatedb.o pfamhmmer1db.o  seqaligndisplay.o threestatehsp.o  wise2xhmmer2.o $(LFLAGS) -lhmmer -ldyna_glib -ldyna -lwisebase  $(EXTRALIBS)

makerandomdb : makerandomdb.o 
	$(CC) -o makerandomdb makerandomdb.o $(LFLAGS) -ldyna_glib -ldyna -lwisebase  $(EXTRALIBS)

genewise : genewise.o $(NEWG) 
	$(CC) -o genewise genewise.o $(NEWG) $(LFLAGS) -lhmmer -ldyna_glib -ldyna_glib -ldyna_glib -ldyna -lwisebase  $(EXTRALIBS)

cdnawise : cdnawise.o $(NEWG) cdnawise10.o 
	$(CC) -g -o cdnawise cdnawise.o cdnawise10.o $(NEWG) $(LFLAGS) -lhmmer -ldyna_glib -ldyna -lwisebase  $(EXTRALIBS)

fivestar : fivestar.o fivestatemodel.o fivestate.o threestatemodel.o threestatedb.o pfamhmmer1db.o wise2xhmmer2.o seqaligndisplay.o standardout.o threestatedp.o
	$(CC) -g -o fivestar fivestar.o fivestatemodel.o fivestate.o threestatedb.o threestatemodel.o pfamhmmer1db.o wise2xhmmer2.o seqaligndisplay.o standardout.o  threestatedp.o $(LFLAGS) -lhmmer -ldyna_glib -ldyna -lwisebase  $(EXTRALIBS)


fivestarsearch : fivestarsearch.o fivestatemodel.o fivestate.o threestatemodel.o threestatedb.o pfamhmmer1db.o wise2xhmmer2.o seqaligndisplay.o standardout.o threestatedp.o
	$(CC) -g -o fivestarsearch fivestarsearch.o fivestatemodel.o fivestate.o threestatedb.o threestatemodel.o pfamhmmer1db.o wise2xhmmer2.o seqaligndisplay.o standardout.o  threestatedp.o $(LFLAGS) -lhmmer -ldyna_glib -ldyna -lwisebase  $(EXTRALIBS)

fivestarscan : fivestarscan.o fivestatemodel.o fivestate.o threestatemodel.o threestatedb.o pfamhmmer1db.o wise2xhmmer2.o seqaligndisplay.o standardout.o threestatedp.o
	$(CC) -g -o fivestarscan fivestarscan.o fivestatemodel.o fivestate.o threestatedb.o threestatemodel.o pfamhmmer1db.o wise2xhmmer2.o seqaligndisplay.o standardout.o  threestatedp.o $(LFLAGS) -lhmmer -ldyna_glib -ldyna -lwisebase  $(EXTRALIBS)


editdist : editdist.o editdistdp.o standardout.o seqaligndisplay.o 
	$(CC) -g -o editdist editdist.o editdistdp.o standardout.o seqaligndisplay.o $(LFLAGS) -ldyna_glib -ldyna -lwisebase 


efgw   : genewise.o $(NEWG) 
	$(CC) -g -o genewise genewise.o  $(NEWG)  $(LFLAGS) -lhmmer -ldyna_glib -ldyna -lwisebase  $(EXTRALIBS) -lefence


genewisedb : genewisedb.o $(NEWG) 
	$(CC) -g -o genewisedb genewisedb.o  $(NEWG)  $(LFLAGS) -lhmmer -ldyna_glib -ldyna -lwisebase  $(EXTRALIBS) 

estwisedb : estwisedb.o $(NEWG) 
	$(CC) -g -o estwisedb estwisedb.o $(NEWG)  $(LFLAGS) -lhmmer -ldyna_glib -ldyna -lwisebase  $(EXTRALIBS)

edb : estwisedb.o $(NEWG) 
	$(CC) -g -o estwisedb estwisedb.o $(NEWG)  $(LFLAGS) -lhmmer -ldyna_glib -ldyna -lwisebase  $(EXTRALIBS) -lefence

estwise : estwise.o $(NEWG) 
	$(CC) -g -o estwise estwise.o  $(NEWG) $(LFLAGS) -lhmmer -ldyna_glib -ldyna -lwisebase  $(EXTRALIBS)


amplimer_resolver : amplimer_resolver.o dnaalign.o dnamatcher.o 
	$(CC) -g -o amplimer_resolver amplimer_resolver.o dnamatcher.o dnaalign.o $(LFLAGS) -ldyna_glib -ldyna -lwisebase  $(EXTRALIBS)

genomewise : genomewise.o genomewise9.o genome_evidence.o est_evidence.o geneutil.o geneoutput.o
	$(CC) -g -o genomewise genomewise.o genomewise9.o genome_evidence.o est_evidence.o geneoutput.o geneutil.o $(LFLAGS) -ldyna_glib -ldyna -lwisebase  $(EXTRALIBS)

sywise : sywise.o sywise20.o syexonmodel.o genestats.o pwmdna.o geneutil.o standardout.o
	$(CC) -g -o sywise sywise.o sywise20.o syexonmodel.o genestats.o pwmdna.o  standardout.o geneutil.o $(LFLAGS) -ldyna_glib -ldyna -lwisebase  $(EXTRALIBS)

alignwise : alignwise.o alignwisedp.o aligngenemodel.o genestats.o pwmdna.o geneutil.o standardout.o geneoutput.o
	$(CC) -g -o alignwise alignwise.o alignwisedp.o aligngenemodel.o genestats.o pwmdna.o geneutil.o standardout.o geneoutput.o $(LFLAGS) -ldyna_glib -ldyna -lwisebase  $(EXTRALIBS)

test_splice : test_splice.o  aligngenemodel.o genestats.o pwmdna.o geneutil.o standardout.o
	$(CC) -g -o test_splice test_splice.o aligngenemodel.o genestats.o pwmdna.o geneutil.o standardout.o $(LFLAGS) -ldyna_glib -ldyna -lwisebase  $(EXTRALIBS)

statwise : statwise.o statwise10.o syexonmodel.o genestats.o pwmdna.o geneutil.o
	$(CC) -g -o statwise statwise.o statwise10.o syexonmodel.o genestats.o pwmdna.o geneutil.o $(LFLAGS) -ldyna_glib -ldyna -lwisebase  $(EXTRALIBS)

pseudowise : pseudowise.o pseudowise7.o $(NEWG) 
	$(CC) -g -o pseudowise pseudowise.o pseudowise7.o $(NEWG) $(LFLAGS) -lhmmer -ldyna_glib -ldyna -lwisebase  $(EXTRALIBS)

clonewise : clonewise.o clonewisedp.o localclonewisedp.o mapstruct.o
	$(CC) -g -o clonewise clonewise.o clonewisedp.o localclonewisedp.o mapstruct.o $(LFLAGS) -ldyna_glib -ldyna -lwisebase  $(EXTRALIBS)


HMMER_LIBS = ../HMMer2/
HMMER_INCLUDE = ../HMMer2/
HMMER_DEFINE = 


CC = cc
#CFLAGS = -c -O2  -pg -I../base/ -I../dynlibsrc/
#CFLAGS = -c -O2 -DPTHREAD -DHAS_PTHREAD_SETSCOPE  -DUNIX  -I../base/ -I../dynlibsrc/ -I.
CFLAGS = -c -g3 -O2  -DUNIX  -I../base/ -I../dynlibsrc/ -I. -pthread `glib-config --cflags`
#CFLAGS = -c  -g -DHAS_PTHREAD_SETSCOPE -DUNIX -I../base/ -I../dynlibsrc/ `glib-config --cflags`
#LFLAGS =  -L../base/ -L../dynlibsrc/ -L../HMMer2/  

LFLAGS =  -L../base/ -L../dynlibsrc/ -L$(HMMER_LIBS) -lpthread `glib-config --libs`
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







