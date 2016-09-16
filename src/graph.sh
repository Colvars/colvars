#!/bin/sh

echo 'digraph d {\n  rankdir = "LR";' > deps.gv
sed -f deps.sed < colvardeps.cpp >> deps.gv
awk '/^  cvb_.* \[/{list = list "\"" $1 "\" ; "} END {print "  { rank = same; " list "}"}' deps.gv >> deps.gv
#awk '/^  cv_.* \[/{list = list "\"" $1 "\" ; "} END {print "  { rank = same; " list "}"}' deps.gv >> deps.gv
awk '/^  cvc_.* \[/{list = list "\"" $1 "\" ; "} END {print "  { rank = same; " list "}"}' deps.gv >> deps.gv
awk '/^  ag_.* \[/{list = list "\"" $1 "\" ; "} END {print "  { rank = same; " list "}"}' deps.gv >> deps.gv
echo '}' >> deps.gv

dot deps.gv  -o deps.svg -Tsvg
