# LaTeXML-based HTML build of the Colvars manual (work in progress)

## Goal

Replace the tex4ht (`htlatex`) HTML build of `doc/colvars-refman-<engine>.tex`
with [LaTeXML](https://math.nist.gov/~BMiller/LaTeXML), producing:

- **multiple pages**, split per section, instead of one large page;
- **MathML** output (rendered natively by all major browsers since 2023);
  MathJax is loaded by default only to improve the typesetting, and can be
  turned off;
- clean semantic HTML5 that is styled with CSS.

The PDF build (pdflatex) and the existing tex4ht build must keep working
during the transition.

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
    top-level entries. Colors and boxes follow the PDF (section headings,
    keyword names, the three kinds of example blocks), plus the sans-serif
    font and keyword-item spacing of the tex4ht CSS (`colvars-refman-css.tex`
    and the `html_doc_pages` branch).
  - `mdframed.sty.ltxml`: a minimal stub. Without it, LaTeXML loads the raw
    package, which pulls in xparse and the whole expl3 kernel and hits the
    600 s timeout.
- **`colvars-refman.tex`** defines two conditionals near the top:
  - `\iflatexml`: true only under LaTeXML. The LaTeX definitions of the
    macros above are wrapped in `\iflatexml\else ... \fi`. **Keep the
    argument lists in sync with `colvars-latexml.sty.ltxml`.**
  - `\ifcvhtml`: true for any HTML output (LaTeXML or tex4ht), for example
    to choose PNG instead of PDF figures and `.html` instead of `.pdf` links.
    Use `\ifdefined\HCode` only for tex4ht-specific raw-HTML code.
  - Under LaTeXML, the title page is replaced by `\title`/`\author`/
    `\maketitle`. The links to the other engines' manuals were factored into
    `\cvmanuallinks` (shared by all builds).

## Building

From `doc/`:

```
make latexml                          # all engines (or RELEASE_ENGINE=namd)
make latexml LATEXML_MATHJAX=no       # native MathML, no JavaScript
make latexml LATEXML_SPLITAT=subsection  # more, shorter pages (default: section)
make latexml LATEXML_URLSTYLE=server  # links to folders (web server only)
```

Changing these variables does not trigger a rebuild: delete `html-<engine>`
first. The default `LATEXML_URLSTYLE=file` (links include `index.html`) works
both locally and on a web server; LaTeXML's own default, `server`, does not
work for local browsing.

Pages are split at sections (or subsections) and named after their labels
(`--splitnaming=labelrelative`, e.g. `sec_colvar/sec_cvc_distances.html`);
with `--navigationtoc=context`, the menu shows the top level plus the
ancestors, siblings and children of the current page, so it expands level by
level as in the LaTeXML manual (these options are taken from LaTeXML's
`doc/makemanual`). Section labels must therefore not be duplicated by
keyword labels (`\labelkey`) on the same heading, or the page is named after
the keyword.

Output goes to `html-<engine>/index.html`, with the LaTeXML log in
`html-<engine>/latexml.log`; the build prints its count of errors and
warnings. A build takes about 80 s per engine, mostly spent processing
the `.bib` files. After LaTeXML, `fix_index.py` rewrites `idx.html` so that
each index entry links to its keyword directly, instead of to the
"Nth item" label that LaTeXML generates.

With `LATEXML_MATHJAX=yes` (the default), every page loads MathJax
(MathML input) from the jsDelivr CDN. MathJax re-typesets LaTeXML's MathML
with better spacing and fonts than the browsers' native MathML rendering.
Without it, the pages need no JavaScript and the math is shown immediately.

## Status (2026-11-07)

- Last run (NAMD manual): **0 errors, 2 warnings**, both from the
  multi-line `aligned` equation near line 5407 of `colvars-refman-main.tex`
  (`\left.`/`\right.` split across lines; it still renders). Run time is
  about 80 s.
- pdflatex builds correctly (169 pages, no undefined references).
  htlatex stops at natbib's "Bibliography not compatible with author-year
  citations" prompt when run non-interactively. The committed (HEAD)
  sources do the same, so this is not a regression.
- Output (NAMD): 71 pages, split at subsections; the largest content page
  is ~120 kB (rotation components).
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
- The fixes in `postprocess_html.sh` (ligatures, `¡`/`¿`, `∖`) are not
  needed with LaTeXML: none of these characters appear in the output.

## Remaining work

1. Check the rest of the output (index `idx.html`, bibliography
   `bib.html`, figures) and decide whether S4/S6 should be split further.
   The Makefile target must depend on `images` (`*-512px.png`).
2. Stylesheet: done for layout, headings, keywords and examples. Still to
   check: index and bibliography pages, title page/cover image, tables.
   Citations come out author-year ("Phillips et al. (2020)") while the PDF
   is numeric (`unsrt`).
3. **Keep old URLs working.** Fragment links such as
   `colvars-refman-vmd.html#sec:dashboard` or `...#sec:colvar` are hard-coded
   in `vmd/cv_dashboard/*.tcl`, `vmd/cv_dashboard/README.md`,
   `colvartools/README.md`, and in the already released NAMD and VMD
   manuals (`namd/ug/ug_colvars.tex`, `vmd/doc/ug_colvars.tex`). Plan: keep
   a small `colvars-refman-<engine>.html` page with a generated JSON map from
   label to `page#id`, and a few lines of JS that redirect according to
   `location.hash`. LaTeXML generates its own ids (`S3.SS2`, ...) and
   sanitizes labels (`:` and `|` are not valid in XML ids), so the map must
   be built from LaTeXML's output (its intermediate XML keeps the
   `labels="LABEL:..."` attributes). Page file names already come from the
   section labels (`--splitnaming=labelrelative`), which helps.
4. Makefile: the `latexml` target exists. Still to do: replace the `html`
   (htlatex) target and `postprocess_html.sh` with it, and update
   `webpage`/`install` / `update_website.sh` (don't copy `LaTeXML.cache`,
   `latexml.log`, `latexml.out`).
5. Possibly speed up the bibliography step (e.g. pass a pre-converted
   bibliography, or only the cited entries).
6. All engines build (`make -j5 latexml`, ~90 s): 0 errors, the same 2 math
   warnings, and no broken internal links (checked with a link checker
   over all pages). Engine-specific sections are present (command-line
   interface for LAMMPS/NAMD/VMD, dashboard for VMD). The links to the other
   engines' manuals (`\cvmanuallinks`) are relative under LaTeXML:
   `../html-<engine>/index.html`, and `../colvars-refman-<engine>.pdf` for
   the PDF. They work in `doc/` and on the website if the `html-<engine>/`
   folders are installed next to the PDFs in the release folder. The PDF's
   own "HTML version" link (and the tex4ht build) still use the old
   single-page URLs: switch them when publishing (item 4).

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

- Debian/Ubuntu/Mint: `sudo apt install latexml`. If that version is older
  than 0.8.8, install it as a user instead:
  `cpanm -l ~/perl5 LaTeXML` (get `cpanm` from https://cpanmin.us). Install
  the compiled dependencies from apt first: `libxml-libxml-perl`,
  `libxml-libxslt-perl`, `libjson-xs-perl`, `libtext-unidecode-perl`,
  `libimage-magick-perl`. Then put
  `eval "$(perl -I$HOME/perl5/lib/perl5 -Mlocal::lib=$HOME/perl5)"` in
  `~/.bashrc`.
- Mageia does not package LaTeXML; the user-level cpanm install works, with
  `perl-JSON-XS`, `perl-Text-Unidecode` and `perl-Image-Magick` from urpmi
  (`perl-devel` is missing headers, so JSON::XS can't be built from CPAN).
