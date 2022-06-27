.PHONY: all pdf clean clean-all update-code-refs update-cvscript-cmdline-doc

PDF=colvars-refman-lammps.pdf colvars-refman-namd.pdf colvars-refman-vmd.pdf colvars-refman-gromacs.pdf colvars-refman-tinkerhp.pdf
BIBTEX=colvars-refman.bib

all: pdf

pdf: $(PDF)

%.pdf: %.tex $(BIBTEX) colvars-refman-main.tex colvars-refman.tex
	latexmk -pdf $<
## Use the lines below if latexmk is not available (better install latexmk!)
# pdflatex $<
# bibtex $(basename $<)
# makeindex $(basename $<).idx
# pdflatex $<
# pdflatex $<

clean: 
	rm -f *.aux *.blg *.log *.toc *.out *.idx *.ilg *.ind *.brf *.4tc *.4ct *.css *.dvi *.html *.idv *.lg *.xref *.tmp eulerangles-512px.png

clean-all: clean
	rm -f $(PDF)

version-list:
	./print_engine_versions.sh > ../README-versions.md


update-code-refs:
	python3 extract_code_refs.py > ../src/colvarmodule_refs.h


update-cvscript-cmdline-doc:
	tclsh gen_cvscript_doc.tcl
