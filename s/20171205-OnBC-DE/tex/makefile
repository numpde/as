SHELL=/bin/bash
OPEN=xdg-open
-include ~/.Makefile.local

MAIN=p

default:
	make l
	make b
	make l
	make l

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


