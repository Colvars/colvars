.PHONY: all pdf html latexml doxygen webpage webpage-legacy images FORCE \
	install install-images clean clean-all \
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
# Multi-page HTML manual of each engine, built with LaTeXML (see
# latexml/README.md): html-<engine>/index.html and one page per section
HTML_DIRS = $(foreach _,${ENGINES},html-$_)
HTML = $(addsuffix /index.html,$(HTML_DIRS))

IMAGES = cover-512px.jpg eulerangles-512px.png spin-tilt-scheme-512px.png

BIBTEX = colvars-refman.bib colvars-code-refs.bib
# Files read by colvars-refman-<engine>.tex (\input), for both PDF and HTML
TEX_INPUTS = colvars-refman.tex colvars-refman-main.tex \
	colvars-engine-flags.tex cv_version.tex cvscript-tcl.tex cvscript-fix-modify.tex

ifeq ($(COLVARS_RELEASE),"master")
all: pdf html doxygen
else
all: pdf html
endif

pdf: images $(PDF)
html: images $(HTML)
latexml: html
images: $(IMAGES)

%.pdf: %.tex $(TEX_INPUTS) $(BIBTEX) Makefile
	@echo "Building $@ with release label \"${COLVARS_RELEASE}\"" ; \
	latexmk -pdf $<
## Use the lines below if latexmk is not available (better install latexmk!)
# pdflatex $<
# bibtex $(basename $<)
# makeindex $(basename $<).idx
# pdflatex $<
# pdflatex $<

LATEXMLC = latexmlc
# Set LATEXML_MATHJAX to anything but "yes" to use the browser's native
# MathML rendering (no JavaScript)
LATEXML_MATHJAX = yes
# Exact version, so that published pages do not change on their own
MATHJAX_URL = https://cdn.jsdelivr.net/npm/mathjax@3.2.2/es5/mml-chtml.js
# One page per "section" or per "subsection" (the menu expands level by level)
LATEXML_SPLITAT = section
# "file": links include index.html, so pages can be browsed locally;
# "server": links to folders (e.g. sec_colvar/), only through a web server
# (used for the website, see "webpage" below)
LATEXML_URLSTYLE = file
# Colvars logo: the compact version is the favicon; the full one is shown
# above the menu and next to the title of the front page (see the CSS)
LATEXML_LOGO = logo/colvars-logo.svg
LATEXML_ICON = logo/colvars-logo-compact.svg
LATEXML_OPTS = --path=latexml --preload=colvars-latexml.sty \
	--css=latexml/colvars-refman.css --icon=$(LATEXML_ICON) \
	$(if $(filter yes,$(LATEXML_MATHJAX)),--javascript=$(MATHJAX_URL)) \
	--format=html5 --split --splitat=$(LATEXML_SPLITAT) \
	--splitnaming=labelrelative --urlstyle=$(LATEXML_URLSTYLE) \
	--navigationtoc=context

# Rebuild the HTML when the LaTeXML options above change
LATEXML_OPTS_STAMP = latexml-options.txt
$(LATEXML_OPTS_STAMP): FORCE
	@echo '$(strip $(LATEXML_OPTS))' | cmp -s - $@ || \
		echo '$(strip $(LATEXML_OPTS))' > $@

# The build fails if LaTeXML reports errors (warnings are only counted).
# LaTeXML.db (--dbfile) maps the LaTeX labels to pages, for the redirect pages
# that replace the former single-page manual (latexml/make_redirect.pl)
html-%/index.html: colvars-refman-%.tex $(TEX_INPUTS) \
		$(BIBTEX) $(IMAGES) $(wildcard latexml/*) \
		$(LATEXML_LOGO) $(LATEXML_ICON) $(LATEXML_OPTS_STAMP) Makefile
	@echo "Building $@ with release label \"${COLVARS_RELEASE}\""
	rm -rf html-$* ; mkdir -p html-$*
	$(LATEXMLC) $(LATEXML_OPTS) --log=html-$*/latexml.log \
		--dbfile=html-$*/LaTeXML.db --destination=$@ $< \
		> html-$*/latexml.out 2>&1 || (tail -20 html-$*/latexml.out; rm -f $@; false)
	python3 latexml/fix_index.py html-$*/idx.html
	mkdir -p html-$*/logo && cp -p $(LATEXML_LOGO) html-$*/logo/
	@# LaTeXML labels any favicon as PNG or ICO; this one is SVG
	find html-$* -name '*.html' -exec sed -i \
		's|\(rel="shortcut icon" href="[^"]*\.svg" type="\)[^"]*"|\1image/svg+xml"|' {} +
	@errors=$$(grep -cE '^(Error|Fatal):' html-$*/latexml.log) ; \
	echo "$@: $${errors} errors," \
		"$$(grep -cE '^Warning:' html-$*/latexml.log) warnings (see html-$*/latexml.log)" ; \
	if [ $${errors} -gt 0 ] ; then \
		grep -E -A3 '^(Error|Fatal):' html-$*/latexml.log ; rm -f $@ ; false ; fi

spin-tilt-scheme-512px.png: spin-tilt-scheme.pdf
	convert -flatten -resize 512x300 $^ $@

eulerangles-512px.png: eulerangles.pdf
	convert -flatten -resize 512x576 $^ $@

clean:
	rm -f *.fdb_latexmk *.fls *.aux *.bbl *.blg *.log *.toc *.out *.idx *.ilg *.ind *.brf *.4tc *.4ct *.dvi *.idv *.lg *.xref *.tmp eulerangles-512px.png

clean-all: clean
	rm -f $(PDF) $(LATEXML_OPTS_STAMP)
	rm -rf $(HTML_DIRS)
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

# On the website, the HTML pages link to folders (e.g. sec_colvar/) instead of
# files: this rebuilds the HTML if it was built with another LATEXML_URLSTYLE
webpage webpage-legacy install: LATEXML_URLSTYLE = server

RELEASE_TREE = $(COLVARS_WEBSITE_TREE)/$(COLVARS_RELEASE)

# Copy the HTML manuals to the release folder, replacing the previous ones
# (page names follow the section labels, so pages can disappear), without
# LaTeXML's cache, database and logs.  The former single-page manual
# colvars-refman-<engine>.html becomes a redirect page, and its CSS is removed
define install_html
mkdir -p $(RELEASE_TREE) ; \
for engine in $(ENGINES) ; do \
  rm -rf $(RELEASE_TREE)/html-$${engine} && \
  cp -p -r html-$${engine} $(RELEASE_TREE)/ && \
  rm -f $(RELEASE_TREE)/html-$${engine}/LaTeXML.cache \
    $(RELEASE_TREE)/html-$${engine}/LaTeXML.db \
    $(RELEASE_TREE)/html-$${engine}/latexml.log \
    $(RELEASE_TREE)/html-$${engine}/latexml.out \
    $(RELEASE_TREE)/colvars-refman-$${engine}.css && \
  perl latexml/make_redirect.pl html-$${engine}/LaTeXML.db $${engine} \
    html-$${engine}/ > $(RELEASE_TREE)/colvars-refman-$${engine}.html \
    || exit 1 ; \
done
endef

ifeq ($(COLVARS_RELEASE),master)

install-images:
	cp -p -f $(IMAGES) $(COLVARS_WEBSITE_TREE)/images/

webpage: all doxygen
	$(install_html)
	cp -p -f $(PDF) $(RELEASE_TREE)/ ; \
	rm -fr $(COLVARS_WEBSITE_TREE)/doxygen/html/* ; \
	cp -p -r doxygen/html/* $(COLVARS_WEBSITE_TREE)/doxygen/html/

# Update pages with older-style paths (without branch/release name): the
# single-page manuals in colvars-refman-<engine>/ become redirect pages to the
# multi-page manuals of the master release
webpage-legacy: all
	cp -p -f $(PDF) $(COLVARS_WEBSITE_TREE)/pdf ; \
	for engine in $(ENGINES) ; do \
	  legacy=$(COLVARS_WEBSITE_TREE)/colvars-refman-$${engine} ; \
	  mkdir -p $${legacy} && \
	  rm -f $${legacy}/colvars-refman-$${engine}.css && \
	  perl latexml/make_redirect.pl html-$${engine}/LaTeXML.db $${engine} \
	    ../$(COLVARS_RELEASE)/html-$${engine}/ \
	    > $${legacy}/colvars-refman-$${engine}.html || exit 1 ; \
	done

install: webpage webpage-legacy

else

webpage: all
	$(install_html)
	cp -p -f $(PDF) $(RELEASE_TREE)/

install: webpage
endif


endif
