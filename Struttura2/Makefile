NAME = Struttura_1_2
TEX  = $(NAME).tex
DVI  = $(TEX:.tex=.dvi)
PS   = $(DVI:.dvi=.ps)
PDF  = $(PS:.ps=.pdf)
BIB  = biblio
BIBTEX = $(BIB).tex
BIBBIB = $(BIB).bib
BBL = $(NAME:.tex=.bbl)
CONVERSION = convertbiblio.py

all: $(TEX)
	pdflatex $(TEX) ; rm -rf *.log *.out

clean:
	rm *~ *.toc *.pdf *.blg *.log *.aux $(NAME).out $(DVI) $(PS) $(BIB)

