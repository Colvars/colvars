
PDF=colvars-refman-lammps.pdf  colvars-refman-namd.pdf 
BIBTEX=colvars-refman.bib

all: $(PDF)

%.pdf: %.tex $(BIBTEX)
	pdflatex $<
	bibtex $(basename $<)
	pdflatex $<
	echo makeindex $(basename $<)
	pdflatex $<
	pdflatex $<

clean: 
	rm -f *.aux *.bbl *.blg *.log *.toc *.out *.idx *.ilg *.ind *.brf 

clean-all: 
	rm -f $(PDF)



