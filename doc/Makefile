
PDF=colvars-refman-lammps.pdf  colvars-refman-namd.pdf  colvars-refman.pdf
BIBTEX=colvars.bib

all: $(PDF)

%.pdf: %.tex $(BIBTEX)
	pdflatex $<
	bibtex $(basename $<)
	makeindex $(basename $<)
	pdflatex $<
	pdflatex $<

clean: 
	rm -f *.aux *.bbl *.blg *.log *.pdf *.toc *.out *.idx *.ilg *.ind *.brf

