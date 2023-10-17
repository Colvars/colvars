.PHONY: all pdf html doxygen webpage webpage-legacy images \
	clean clean-all \
	update-code-refs update-cvscript-cmdline-doc updates

PDF = \
	colvars-refman-gromacs.pdf \
	colvars-refman-lammps.pdf \
	colvars-refman-namd.pdf \
	colvars-refman-tinkerhp.pdf \
	colvars-refman-vmd.pdf

HTML = \
	colvars-refman-gromacs.html \
	colvars-refman-lammps.html \
	colvars-refman-namd.html \
	colvars-refman-vmd.html

CSS = $(HTML:.html=.css)

BIBTEX = colvars-refman.bib

COLVARS_RELEASE = $(shell git symbolic-ref --short -q HEAD)
# If we are not working on a branch, try a tag instead

ifeq ($(COLVARS_RELEASE),)
COLVARS_RELEASE = $(git describe --tags --exact-match)
endif

export COLVARS_RELEASE


all: pdf html doxygen

updates: version-list update-code-refs update-cvscript-cmdline-doc

pdf: $(PDF)
html: $(HTML)

%.pdf: %.tex $(BIBTEX) colvars-refman-main.tex colvars-refman.tex
	@echo "Building $@ with release label \"${COLVARS_RELEASE}\"" ; \
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
	@echo "Building $@ with release label \"${COLVARS_RELEASE}\"" ; \
	$(HTLATEX) $(notdir $<) $(HTLATEX_OPTS); \
	bash postprocess_html.sh

clean: 
	rm -f *.fdb_latexmk *.fls *.aux *.bbl *.blg *.log *.toc *.out *.idx *.ilg *.ind *.brf *.4tc *.4ct *.dvi *.idv *.lg *.xref *.tmp eulerangles-512px.png

clean-all: clean
	rm -f $(PDF) $(HTML) $(CSS)

version-list:
	./print_engine_versions.sh > ../README-versions.md

update-code-refs:
	python3 extract_code_refs.py > ../src/colvarmodule_refs.h

update-cvscript-cmdline-doc:
	tclsh gen_cvscript_doc.tcl


doxygen: doxygen/html/index.html

doxygen/html/index.html: ../src/*.h doxygen/Doxyfile
	cd doxygen; doxygen


ifneq ($(COLVARS_WEBSITE_TREE),)

IMAGES = cover-512px.jpg eulerangles-512px.png

# Install docs in the website tree, under the folder named $(COLVARS_RELEASE).
# Note: to produce the documentation for a specific release of a MD engine,
# check out the branch or tag corresponding to that version, clean and rebuild.
webpage: all images
	mkdir -p $(COLVARS_WEBSITE_TREE)/$(COLVARS_RELEASE) ; \
	cp -p -f $(PDF) $(HTML) $(CSS) $(COLVARS_WEBSITE_TREE)/$(COLVARS_RELEASE)/ ; \
	cp -p -f $(addprefix $(COLVARS_WEBSITE_TREE)/images/, $(IMAGES)) $(COLVARS_WEBSITE_TREE)/$(COLVARS_RELEASE)/ ; \
	rm -fr $(COLVARS_WEBSITE_TREE)/doxygen/html/* ; \
	cp -p -r doxygen/html/* $(COLVARS_WEBSITE_TREE)/doxygen/html/

# Update HTML pages with older-style paths (without branch/release name)
webpage-legacy: webpage
	cp -p -f $(PDF) $(COLVARS_WEBSITE_TREE)/pdf ; \
	for engine in lammps namd vmd ; do \
	cp -p -f colvars-refman-$${engine}.{html,css} $(COLVARS_WEBSITE_TREE)/colvars-refman-$${engine}/ ; \
	cp -p -f $(addprefix $(COLVARS_WEBSITE_TREE)/images/, $(IMAGES)) $(COLVARS_WEBSITE_TREE)/colvars-refman-$${engine}/ ; \
	done


endif
