## Building the documentation

This folder contains LaTeX files needed to build the Colvars documentation in PDF and HTML formats.  Not all installation of TeXlive contain functional code, thus we highly recommend using the container used in development.
```
$ apptainer pull texlive.sif oras://ghcr.io/Colvars/devel-containers:texlive
$ apptainer exec texlive.sif make all
```
Note: if the working folder is not under `/home`, Apptainer may not automatically bind it in the container: in that case, try using `--bind <mount_point>`.  Likewise, adapt this step for different container engines (e.g. Docker, Podman, etc).

The GNU Make variable `COLVARS_RELEASE` will be automatically taken from the Git branch or tag name.  Alternatively, feel free to define it as environment variable, for example:
```
$ export COLVARS_RELEASE=gromacs-2024
$ apptainer exec texlive.sif make all
```
Note that when the release is prefixed by the name of a specific engine, only the doc for that engine gets built.


## Updating the website

To update the website, clone its [repository](https://github.com/Colvars/colvars.github.io) and set the environment variable accordingly before building the documentation:
```
$ export COLVARS_WEBSITE_TREE=<path to website repo>
$ make all
$ make install
```


## LaTeX file organization

The file `colvars-refman-main.tex` is the main file, which includes content that is conditionally rendered by LaTeX depending on the version of the documentation being generated.  There are specific LaTeX files for each backend that include this file to generate the specific version of the doc, e.g. `colvars-refman-gromacs.tex` to generate the doc for Gromacs users.


## Citations to software

The file `colvars-code-refs.bib` contains BibTeX entries for all papers or preprints describing software implementations of specific features, either implemented within Colvars or outside (e.g. the publications of the back-end code and some of their features are included as well).  One or more comment lines (preceded by `%`) describe the feature(s) associated to that publication.  Descriptions must be short and unique, because they are used by the code to generate a citation report during a run.

Each BibTeX entry must contain a `url = {https://...}` field, pointing to a URL that can be used to access the corresponding publication (preferably a `https://doi.org/...` link).

The contents of the file `colvars-code-refs.bib` are converted to C++ code using:
```
make update-code-refs
```
