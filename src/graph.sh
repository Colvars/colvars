#!/bin/sh

# Create a dependency graph using the dot utility from the Graphviz package

echo 'digraph d {
rankdir = "LR";' > deps.gv

# collect data from the all files with calls to init_feature()
files=$(grep -lE 'init_feature\(f_' *.cpp)
for f in $files
do
  echo "Extracting deps from source file $f"
  gawk -f process_deps_graph.awk $f >> deps.gv
done

# set dependencies of the same level objects to the same rank
awk '/^  cvb_/{list = list "\"" $1 "\" ; "} END {print "  { rank = same; " list "}"}' deps.gv >> deps.gv
# There are too many dependencies within colvars so they get tangled if plotted on one rank
#awk '/^  cv_./{list = list "\"" $1 "\" ; "} END {print "  { rank = same; " list "}"}' deps.gv >> deps.gv
awk '/^  cvc_./{list = list "\"" $1 "\" ; "} END {print "  { rank = same; " list "}"}' deps.gv >> deps.gv
awk '/^  ag_./{list = list "\"" $1 "\" ; "} END {print "  { rank = same; " list "}"}' deps.gv >> deps.gv
echo '}' >> deps.gv

dot deps.gv  -o deps.svg -Tsvg

echo "Output written to deps.gv and deps.svg"
