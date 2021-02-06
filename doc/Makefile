
PDF=colvars-refman-lammps.pdf colvars-refman-namd.pdf colvars-refman-vmd.pdf colvars-refman-gromacs.pdf
BIBTEX=colvars-refman.bib

all: pdf

pdf: $(PDF)

%.pdf: %.tex $(BIBTEX) colvars-refman-main.tex colvars-refman.tex
	pdflatex $<
	bibtex $(basename $<)
	makeindex $(basename $<).idx
	pdflatex $<
	pdflatex $<

clean: 
	rm -f *.aux *.blg *.log *.toc *.out *.idx *.ilg *.ind *.brf *.4tc *.4ct *.css *.dvi *.html *.idv *.lg *.xref *.tmp eulerangles-512px.png

clean-all: clean
	rm -f $(PDF)

