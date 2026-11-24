# LaTeXML-based HTML build of the Colvars manual

## Goal

The HTML manual of each engine is built from `doc/colvars-refman-<engine>.tex`
with [LaTeXML](https://math.nist.gov/~BMiller/LaTeXML), which replaced the
former single-page tex4ht (`htlatex`) build, producing:

- **multiple pages**, split per section, instead of one large page;
- **MathML** output (rendered natively by all major browsers since 2023);
  MathJax is loaded by default only to improve the typesetting, and can be
  turned off;
- clean semantic HTML5 that is styled with CSS.

The PDF is still built by pdflatex (`latexmk`) from the same sources.

## Design

- **`doc/latexml/`** is passed to LaTeXML with `--path=latexml`. It holds
  LaTeXML bindings (Perl, `*.sty.ltxml`):
  - `colvars-latexml.sty.ltxml`: LaTeXML versions of the manual's macros:
    `\cvsec`, `\cvsubsec`, `\cvsubsubsec` (plain `\section` + `\label`, so
    LaTeXML can see the document structure and split on it), `\key`,
    `\keydef`, `\dupkey`, `\simkey`, `\labelkey`, `\refkey`, the
    `cvcoptions` environment, and the example environments
    (`cvexampleinput`, `mdexampleinput`, `datafileexample`). Examples become
    `<pre class="cv_example ...">`, with `\\` mapped to a newline and `~`
    to a plain space, so text copied from the HTML is a valid input file.
    Keyword headers use the CSS classes `cv_keyword` and `cv_keymeta`.
    The binding is preloaded with `--preload=colvars-latexml.sty`, which
    also loads LaTeXML's `latexml` package and makes `\iflatexml` true.
  - `babel.sty.ltxml`: a stub. LaTeXML 0.8.8's babel binding still runs the
    raw `babel.sty`, which breaks with babel v26.x (100+ errors at
    `\begin{document}`). The manual is English-only.
  - `catchfile.sty.ltxml`: implements `\CatchFileEdef{\cs}{|kpsewhich
    -var-value=VAR}{...}` by reading the environment variable (shell escape
    is not available). Used for `\cvreleaselabel` (from `COLVARS_RELEASE`).
  - `colvars-refman.css`: stylesheet, passed with
    `--css=latexml/colvars-refman.css` (loaded after LaTeXML's own CSS and
    copied to `latexml/` in the output). On screens 50em or wider, the
    navigation menu is a sticky sidebar on the left that scrolls on its own;
    on narrower screens it stays above the content and lists only the
    top-level entries. The menu (20rem wide, muted colors, highlighted
    current page, the Colvars logo above the title link, only two levels:
    sections and subsections), the headings (semibold,
    dark blue, muted numbers, a rule under section titles, no frames) and the
    links (dark blue, underlined only under the cursor) follow the HTML
    version of the NAMD User's Guide; the colors are CSS variables at the top
    of the file. Keyword names and the three kinds of example blocks keep the
    colors and frames of the PDF. The index has up to 4 columns (its page is
    wider), and the "Contents", "Index" and "References" headings are styled
    like the sections. On the front page, whose menu has no title link, the
    logo is next to the title instead.
  - Logo: the compact `doc/logo/colvars-logo-compact.svg` (more legible when
    small) is the favicon (`--icon`); LaTeXML labels any favicon as PNG or
    ICO, so the Makefile fixes its type. The full `doc/logo/colvars-logo.svg`,
    used by the stylesheet, is copied to `logo/` in the output by the
    Makefile.
  - `fix_index.py`: post-processing of the index (see below).
  - `make_redirect.pl`: writes the redirect pages that replace the former
    single-page manuals (see "Publishing").
  - `mdframed.sty.ltxml`: a minimal stub. Without it, LaTeXML loads the raw
    package, which pulls in xparse and the whole expl3 kernel and hits the
    600 s timeout.
- **`colvars-refman.tex`** defines the conditional `\iflatexml` near the
  top, true only under LaTeXML:
  - The LaTeX definitions of the macros above are wrapped in
    `\iflatexml\else ... \fi`. **Keep the argument lists in sync with
    `colvars-latexml.sty.ltxml`.**
  - It also selects HTML-specific content, e.g. PNG instead of PDF figures
    (`*-512px.png`) and relative links between the HTML manuals.
  - Under LaTeXML, the title page is replaced by `\title`/`\author`/
    `\maketitle`. The links to the other engines' manuals were factored into
    `\cvmanuallinks` (shared by the PDF and the HTML).

## Building

From `doc/`, `make` builds the PDFs and the HTML; for the HTML only:

```
make html                          # all engines (or RELEASE_ENGINE=namd)
make html LATEXML_MATHJAX=no       # native MathML, no JavaScript
make html LATEXML_SPLITAT=subsection  # more, shorter pages (default: section)
make html LATEXML_URLSTYLE=server  # links to folders (web server only)
```

(`make latexml` is an alias of `make html`.) The LaTeXML options are saved in
`latexml-options.txt`, so changing these variables rebuilds the HTML. The
default `LATEXML_URLSTYLE=file` (links include `index.html`) works both
locally and on a web server; `server` is used for the website.

Pages are split at sections (or subsections) and named after their labels
(`--splitnaming=labelrelative`, e.g. `sec_colvar.html`);
with `--navigationtoc=context`, the menu shows the top level plus the
ancestors, siblings and children of the current page (down to subsections:
the stylesheet hides deeper levels), so it expands level by level as in the
LaTeXML manual (these options are taken from LaTeXML's
`doc/makemanual`). Section labels must therefore not be duplicated by
keyword labels (`\labelkey`) on the same heading, or the page is named after
the keyword.

Output goes to `html-<engine>/index.html`, with the LaTeXML log in
`html-<engine>/latexml.log`; the build prints its count of errors and
warnings, and **fails if there are errors** (so that CI catches them).
`html-<engine>/LaTeXML.db` (`--dbfile`) records where each LaTeX label ended
up, for the redirect pages. After LaTeXML, `fix_index.py` rewrites `idx.html`
so that each index entry links to its keyword directly, instead of to the
"Nth item" label that LaTeXML generates.

A build takes about 80 s per engine (`make -j` builds the engines in
parallel): about 30 s to digest the TeX sources, 25 s to build the document,
9 s to parse the math and 7 s to write the pages. The bibliography is not the
bottleneck: reading BibTeX's `.bbl` instead of the `.bib` files saves only
about 4 s.

With `LATEXML_MATHJAX=yes` (the default), every page loads MathJax 3.2.2
(MathML input) from the jsDelivr CDN. MathJax re-typesets LaTeXML's MathML
with better spacing and fonts than the browsers' native MathML rendering.
Without it, the pages need no JavaScript and the math is shown immediately.
The version is pinned in `MATHJAX_URL`, so that published pages don't change
on their own.

## Publishing

`make install` (with `COLVARS_WEBSITE_TREE` set, see `doc/README.md`)
rebuilds the HTML with `LATEXML_URLSTYLE=server` if needed, then for each
engine:

- copies `html-<engine>/` to the release folder of the website
  (`<release>/html-<engine>/`, next to the PDFs, as expected by the relative
  links of `\cvmanuallinks`), after removing the previous copy (page names
  follow the section labels, so pages can disappear), and without LaTeXML's
  cache, database and logs;
- replaces the former single-page manual `<release>/colvars-refman-<engine>.html`
  with a redirect page, and removes its CSS. On `master`, the same is done for
  the legacy `colvars-refman-<engine>/colvars-refman-<engine>.html` pages,
  which redirect to `master/html-<engine>/`.

The redirect pages keep old links working, such as
`colvars-refman-vmd.html#sec:dashboard`, which are hard-coded in
`vmd/cv_dashboard/*.tcl`, `vmd/cv_dashboard/README.md`,
`colvartools/README.md`, and in released NAMD and VMD manuals
(`namd/ug/ug_colvars.tex`, `vmd/doc/ug_colvars.tex`). `make_redirect.pl`
reads `LaTeXML.db` and writes a page with a map from each LaTeX label to its
page and fragment, and a few lines of JavaScript that redirect according to
`location.hash`; links without a fragment, or with an unknown one, go to the
front page (also without JavaScript). LaTeXML's ids are not always the labels:
e.g. `sec:cvc_distance` and `colvar|distance` label the same heading, whose
id is `colvar|distance`. Of the 363 labels of the former NAMD manual, 360
redirect to their exact location; the other 3 no longer exist (`sec:toc`,
the label removed below, and a tex4ht artifact).

The release folders of older branches (e.g. `namd-3.0`) keep their
single-page manuals until their documentation is rebuilt.

## Status (2026-11-23)

- All engines: **0 errors, 2 warnings**, both from the multi-line `aligned`
  equation near line 5407 of `colvars-refman-main.tex` (`\left.`/`\right.`
  split across lines; it still renders). No broken internal links (checked
  with a link checker over all pages). Engine-specific sections are present
  (command-line interface for LAMMPS/NAMD/VMD, dashboard for VMD).
- pdflatex builds correctly (169 pages, no undefined references). The PDF's
  "HTML version" link points to `<release>/html-<engine>/`.
- Output (NAMD): 12 pages, split at sections. The largest are
  `sec_colvar.html` (~600 kB) and `sec_colvarbias.html` (~450 kB); this is
  fine (`LATEXML_SPLITAT=subsection` gives shorter pages).
- The index (488 entries for NAMD) links each keyword to its description;
  the bibliography has the same 55 entries as the PDF. Citations are
  author-year ("Coutsias et al. (2004)"), with "Cited by" back-links in the
  bibliography; the PDF is numeric (`unsrt`).
- Source fixes made along the way (they also apply to the PDF):
  - closed the `itemize` before `\cvsubsubsec{Histogramming vector
    variables}`, which was opened before it and closed only after it;
  - `\hyperlink{sec:...}` → `\hyperref[sec:...]` (LaTeXML's `\hyperlink`
    only finds `\hypertarget`s, not labels), and fixed the wrong target
    `sec:colvarbias_metadynamics` → `sec:colvarbias_meta`;
  - `\refkey` for `layeri_activation` was missing its first argument (the
    following `.` was used as the link target);
  - merged math split across several `$...$` (`$f(\xi_{1}$, $\xi_{2}$...`);
  - `gatherVectorColvars`: type "UNIX filename" → "boolean", removed a stray `)`;
  - `colvars-code-refs.bib`: `Taylor & Francis` → `Taylor \& Francis`;
  - the references to `sec:apath_scripted` are now guarded by
    `\cvvmdornamd` like the section itself (they were dangling in Tinker-HP);
  - removed the unused `\labelkey{colvar|NeuralNetwork}` on a section heading.
- Example blocks: `$<$`/`$>$` placeholders are emitted as literal `<`/`>`
  (see `\cv@dollar` in the binding); other math inside examples is kept.
  Curly quotes are turned back into ASCII quotes
  (`cv_example_straight_quotes`), e.g. for the Python example.
- The tex4ht build was removed with its post-processing scripts
  (`postprocess_html.sh`, `fix_mathml.py`, `remove_htmltags.py`,
  `html5mjlatex.cfg`): their fixes (ligatures, `¡`/`¿`, `∖`) are not needed
  with LaTeXML. The tex4ht-specific code was also removed from the sources
  (`\HCode` blocks, `\ifcvhtml`, `colvars-refman-css.tex`, `\pagetoplink`).
- CI: the doc jobs (`build-doc` in `.github/workflows/test-library.yml`,
  `update-doc` in `version-maintenance.yml`) run in
  `ghcr.io/colvars/devel-containers:CentOS9-doc-x86_64`
  (`devel-tools/containers/CentOS9-doc-x86_64.Dockerfile`): CentOS Stream 9
  like the test containers, TeX Live 2025 (frozen) with a curated list of
  packages, LaTeXML 0.8.8.

## Remaining work

1. Push the `CentOS9-doc-x86_64` image to GHCR and check the `build-doc` job
   on a PR; after the merge, remove the `Fedora45-devel-x86_64` image tag.

## LaTeXML 0.8.8 pitfalls found so far

- `\lxWithClass{class}{content}` crashes (`getAttribute` on a text node)
  when the content starts with plain text; use the `\cv@span` constructor
  defined in `colvars-latexml.sty.ltxml` instead.
- The `mode => 'restricted_horizontal'` option in `DefEnvironment` causes a
  fatal error; leave it out.
- `$` is not an active character, so redefine its meaning with
  `Let(T_MATH, ...)`, not by changing its catcode (that misses tokens already
  read as macro arguments). Don't `Let` a saved copy back to `T_MATH`: it
  loops forever; dispatch to `\@dollar@in@normalmode` instead.
- `\hypertarget{id}{}` with empty content creates no anchor, so
  `\hyperlink`s to it fail. Use `\label` + `\lxRef`/`\hyperref[]`.

## Installing LaTeXML

The simplest is to use the CI container (see `doc/README.md`), which has
LaTeXML and the TeX Live packages needed by the manual. Otherwise:

- Debian/Ubuntu/Mint: `sudo apt install latexml`. If that version is older
  than 0.8.8, install it as a user instead:
  `cpanm -l ~/perl5 LaTeXML` (get `cpanm` from https://cpanmin.us). Install
  the compiled dependencies from apt first: `libxml-libxml-perl`,
  `libxml-libxslt-perl`, `libjson-xs-perl`, `libtext-unidecode-perl`,
  `libimage-magick-perl`. Then put
  `eval "$(perl -I$HOME/perl5/lib/perl5 -Mlocal::lib=$HOME/perl5)"` in
  `~/.bashrc`.
- CentOS/RHEL 9: see `devel-tools/containers/CentOS9-doc-x86_64.Dockerfile`
  (compiled dependencies from EPEL/CRB, LaTeXML from CPAN).
