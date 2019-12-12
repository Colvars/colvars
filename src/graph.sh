#!/bin/sh

echo 'digraph d {
rankdir = "LR";' > deps.gv

# create nodes and edges
gawk -f process_deps_graph.awk colvaratoms.cpp >> deps.gv
gawk -f process_deps_graph.awk colvarcomp.cpp >> deps.gv
gawk -f process_deps_graph.awk colvar.cpp >> deps.gv
gawk -f process_deps_graph.awk colvarbias.cpp >> deps.gv

# set dependencies of the same level objects to the same rank
awk '/^  cvb_/{list = list "\"" $1 "\" ; "} END {print "  { rank = same; " list "}"}' deps.gv >> deps.gv
# There are too many dependencies within colvars so they get tangled if plotted on one rank
#awk '/^  cv_./{list = list "\"" $1 "\" ; "} END {print "  { rank = same; " list "}"}' deps.gv >> deps.gv
awk '/^  cvc_./{list = list "\"" $1 "\" ; "} END {print "  { rank = same; " list "}"}' deps.gv >> deps.gv
awk '/^  ag_./{list = list "\"" $1 "\" ; "} END {print "  { rank = same; " list "}"}' deps.gv >> deps.gv
echo '}' >> deps.gv

dot deps.gv  -o deps.svg -Tsvg

echo "Output written to deps.gv and deps.svg"
