
PDF=colvars-refman-lammps.pdf colvars-refman-namd.pdf colvars-refman-vmd.pdf
BIBTEX=colvars-refman.bib

all: $(PDF) html

%.pdf: %.tex $(BIBTEX) colvars-refman-main.tex
	pdflatex $<
	bibtex $(basename $<)
	pdflatex $<
	makeindex $(basename $<)
	pdflatex $<
	pdflatex $<

html: $(BIBTEX) colvars-refman-main.tex colvars-refman-namd.tex colvars-refman-vmd.tex colvars-refman-lammps.tex
	latex2html colvars-refman-namd.tex -html_version 4.0,unicode -local_icons
	latex2html colvars-refman-vmd.tex -html_version 4.0,unicode -local_icons
	latex2html colvars-refman-lammps.tex -html_version 4.0,unicode -local_icons

clean: 
	rm -f *.aux *.blg *.log *.toc *.out *.idx *.ilg *.ind *.brf 
	rm -f */*.aux */*.bbl */*.blg */*.log */*.toc */*.out */*.idx */*.ilg */*.ind */*.brf */WARNINGS

clean-all: 
	rm -f $(PDF)

