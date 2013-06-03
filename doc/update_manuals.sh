for code in "namd" "lammps" ; do
    rm -f colvars-refman-${code}.{aux,bbl,blg,brf,idx,ilg,ind,log,out,toc}
    pdflatex colvars-refman-${code}
    bibtex colvars-refman-${code}
    pdflatex colvars-refman-${code}
    makeindex colvars-refman-${code}
    pdflatex colvars-refman-${code}
    rm -f colvars-refman-${code}.{aux,bbl,blg,brf,idx,ilg,ind,log,out,toc}
done ;
