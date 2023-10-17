.PHONY: all pdf html clean clean-all update-code-refs update-cvscript-cmdline-doc updates

PDF = \
	colvars-refman-gromacs.pdf \
	colvars-refman-lammps.pdf \
	colvars-refman-namd.pdf \
	colvars-refman-tinkerhp.pdf \
	colvars-refman-vmd.pdf \

HTML = \
	colvars-refman-gromacs.html \
	colvars-refman-lammps.html \
	colvars-refman-namd.html \
	colvars-refman-tinkerhp.html \
	colvars-refman-vmd.html \

BIBTEX=colvars-refman.bib

COLVARS_RELEASE := $(shell git symbolic-ref --short -q HEAD)
export COLVARS_RELEASE

all: pdf html

updates: version-list update-code-refs update-cvscript-cmdline-doc

pdf: $(PDF)
html: $(HTML)

%.pdf: %.tex $(BIBTEX) colvars-refman-main.tex colvars-refman.tex
	latexmk -pdf $<
## Use the lines below if latexmk is not available (better install latexmk!)
# pdflatex $<
# bibtex $(basename $<)
# makeindex $(basename $<).idx
# pdflatex $<
# pdflatex $<

# Note: HTML targets rely on up-to-date bbl files; ensuring that by running pdflatex first
HTLATEX = htlatex
HTLATEX_OPTS = "html5mjlatex.cfg, charset=utf-8" " -cunihtf -utf8"
%.html: %.tex colvars-refman-main.tex colvars-refman.tex $(BIBTEX) %.pdf
	$(HTLATEX) $(notdir $<) $(HTLATEX_OPTS); \
	bash postprocess_html.sh

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
