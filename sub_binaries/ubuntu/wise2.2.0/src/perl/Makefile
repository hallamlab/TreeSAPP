
int : 
	cat top.xs ../dynlibsrc/*.xs ../models/*.xs  > Wise2/Wise2.xs
	cat basic.typemap ../dynlibsrc/*.typemap ../models/*.typemap  >  Wise2/typemap
	cat basic.h ../dynlibsrc/*api.t ../models/*api.t ../dynlibsrc/*api.h ../models/*api.h > Wise2/libs/wise2.h
	cp ../dynlibsrc/*.pod Wise2
	cp ../models/*.pod Wise2


libs :  
	perl ./makelib.pl
	cd Wise2/libs ; make libwise2.a



