.PHONY: all pdf html doxygen webpage webpage-legacy images \
	install clean clean-all \
	update-code-refs update-cvscript-cmdline-doc updates update-version-list

ifeq ($(COLVARS_RELEASE),)
# Unless defined otherwise, the release label is the name of the branch
COLVARS_RELEASE = $(shell git symbolic-ref --short -q HEAD)
endif

ifeq ($(COLVARS_RELEASE),)
# If we are not working on a branch, try a tag instead
COLVARS_RELEASE = $(shell git describe --tags --exact-match)
endif

export COLVARS_RELEASE

# Check if this a release of a specific engine
COLVARS_RELEASE_PREFIX := $(firstword $(subst -, , $(COLVARS_RELEASE)))
# Only allow specific cases
ifeq ($(COLVARS_RELEASE_PREFIX),gromacs)
    RELEASE_ENGINE := $(COLVARS_RELEASE_PREFIX)
endif
ifeq ($(COLVARS_RELEASE_PREFIX),lammps)
    RELEASE_ENGINE := $(COLVARS_RELEASE_PREFIX)
endif
ifeq ($(COLVARS_RELEASE_PREFIX),namd)
    RELEASE_ENGINE := $(COLVARS_RELEASE_PREFIX)
endif
ifeq ($(COLVARS_RELEASE_PREFIX),tinkerhp)
    RELEASE_ENGINE := $(COLVARS_RELEASE_PREFIX)
endif
ifeq ($(COLVARS_RELEASE_PREFIX),vmd)
    RELEASE_ENGINE := $(COLVARS_RELEASE_PREFIX)
endif
export RELEASE_ENGINE


ifeq ($(RELEASE_ENGINE),)
# Build doc for all engines
ENGINES = gromacs lammps namd tinkerhp vmd
else
ENGINES = $(RELEASE_ENGINE)
endif

PDF = $(foreach _,${ENGINES},colvars-refman-$_.pdf)
HTML = $(foreach _,${ENGINES},colvars-refman-$_.html)

IMAGES = cover-512px.jpg eulerangles-512px.png

CSS = $(HTML:.html=.css)
BIBTEX = colvars-refman.bib

ifeq ($(COLVARS_RELEASE),"master")
all: pdf html doxygen
else
all: pdf html
endif

pdf: images $(PDF)
html: images $(HTML)
images: $(IMAGES)

%.pdf: %.tex $(BIBTEX) colvars-refman-main.tex colvars-refman.tex Makefile
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
%.html: %.tex colvars-refman-main.tex colvars-refman.tex $(BIBTEX) %.pdf Makefile
	@echo "Building $@ with release label \"${COLVARS_RELEASE}\"" ; \
	$(HTLATEX) $(notdir $<) $(HTLATEX_OPTS); \
	bash postprocess_html.sh $@

eulerangles-512px.png: eulerangles.pdf
	convert -flatten -resize 512x576 $^ $@

clean:
	rm -f *.fdb_latexmk *.fls *.aux *.bbl *.blg *.log *.toc *.out *.idx *.ilg *.ind *.brf *.4tc *.4ct *.dvi *.idv *.lg *.xref *.tmp eulerangles-512px.png

clean-all: clean
	rm -f $(PDF) $(HTML) $(CSS)
	rm -fr doxygen/html
	@echo
	@echo "Note: only removed files that this branch is set to build."

update-version-list:
	./print_engine_versions.sh > ../README-versions.md

update-code-refs:
	python3 extract_code_refs.py > ../src/colvarmodule_refs.h

update-cvscript-cmdline-doc:
	tclsh gen_cvscript_doc.tcl

updates: version-list update-code-refs update-cvscript-cmdline-doc

doxygen: doxygen/html/index.html

doxygen/html/index.html: ../src/*.h doxygen/Doxyfile
	cd doxygen; doxygen


ifeq ($(COLVARS_WEBSITE_TREE),)

install:
	@echo "To install, please define the COLVARS_WEBSITE_TREE environment variable to a worktree of the colvars.github.io repository"

else

# Install docs in the website tree, under the folder named $(COLVARS_RELEASE).
# Note: to produce the documentation for a specific release of a MD engine,
# check out the branch or tag corresponding to that version, clean and rebuild.
ifeq ($(COLVARS_RELEASE),master)

webpage: all doxygen
	mkdir -p $(COLVARS_WEBSITE_TREE)/$(COLVARS_RELEASE) ; \
	cp -p -f $(PDF) $(HTML) $(CSS) $(COLVARS_WEBSITE_TREE)/$(COLVARS_RELEASE)/ ; \
	cp -p -f $(addprefix $(COLVARS_WEBSITE_TREE)/images/, $(IMAGES)) $(COLVARS_WEBSITE_TREE)/$(COLVARS_RELEASE)/ ; \
	rm -fr $(COLVARS_WEBSITE_TREE)/doxygen/html/* ; \
	cp -p -r doxygen/html/* $(COLVARS_WEBSITE_TREE)/doxygen/html/

# Update HTML pages with older-style paths (without branch/release name)
webpage-legacy: all
	cp -p -f $(PDF) $(COLVARS_WEBSITE_TREE)/pdf ; \
	for engine in $(ENGINES) ; do \
	cp -p -f colvars-refman-$${engine}.{html,css} $(COLVARS_WEBSITE_TREE)/colvars-refman-$${engine}/ ; \
	cp -p -f $(addprefix $(COLVARS_WEBSITE_TREE)/images/, $(IMAGES)) $(COLVARS_WEBSITE_TREE)/colvars-refman-$${engine}/ ; \
	done

install: webpage webpage-legacy

else

webpage: all
	mkdir -p $(COLVARS_WEBSITE_TREE)/$(COLVARS_RELEASE) ; \
	cp -p -f $(PDF) $(HTML) $(CSS) $(COLVARS_WEBSITE_TREE)/$(COLVARS_RELEASE)/ ; \
	cp -p -f $(addprefix $(COLVARS_WEBSITE_TREE)/images/, $(IMAGES)) $(COLVARS_WEBSITE_TREE)/$(COLVARS_RELEASE)/

install: webpage
endif


endif
