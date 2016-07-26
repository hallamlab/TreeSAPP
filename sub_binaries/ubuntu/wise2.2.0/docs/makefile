
.SUFFIXES : .tex .dvi .ps

dynamite : dynamite.ps

dynamite.ps : dynamite.dvi
	dvips dynamite.dvi -o



api :
	cat ../models/*.tex ../dynlibsrc/*.tex | perl gettex.pl > temp.tex
	cat wise2api.tex temp.tex apiend.tex > api.tex
	latex api.tex
	latex api.tex
	dvips api.dvi -o

