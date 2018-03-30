SHELL=/bin/bash
OPEN=xdg-open
-include ~/.Makefile.local

MAIN=p

BC_PATH=/home/ra/svn/as/p/single-cell/20171130-BCXX/OUTPUT/

default:
	make i
	make l || make b || make l || exit 0
	make l
	make p

i:
	p="$$(pwd)/../img/"; cd $(BC_PATH) && cp --parents \
		E*/*.pdf \
		E*/*.txt \
		$$p;

p:
	cd ../../; cp --parents */tex/*.pdf ~/public_html/files/s/

e:
	$(OPEN) $(MAIN).tex 1> /dev/null 2> /dev/null &

b:
	bibtex $(MAIN)

l:
	pdflatex -shell-escape -halt-on-error -interaction nonstopmode $(MAIN).tex

eb:
	$(OPEN) --nosplash *.bib 2> /dev/null &

up:
	svn up

co:
	svn commit -m ""

o: 
	$(OPEN) $(MAIN)*.pdf 2>/dev/null &

clean:
	@rm -fv *.bak *.log *.aux *.out *.blg *.dvi *.ps
	@rm -fv *.backup
	@rm -fv *~ 
	@rm -fv *.bib~*
	@rm -fv *.synctex.gz


