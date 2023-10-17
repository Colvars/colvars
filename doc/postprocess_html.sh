#!/bin/bash

set -e

# Find anchors that correspond to LaTeX section labels as indicated by tex4ht
# comments and substitutes back the original LaTeX labels.

# Only labels that are called by actual \refs are detected; for this reason, a
# special HTML tag <selfref> is added in the LaTeX source, and removed later
# by the script remove_selfrefs.py.

for i in *.html
do
  grep --text '^href=.*tex4ht:ref:' $i | \
  sed 's/href="#\(.*\)">.*<!--tex4ht:ref: \(.*\) --><\/a>.*/s\/"\1"\/"\2"\//' | \
  awk 'BEGIN {printf "\"{ "} {printf $0 " ; "} END {printf " }\""}' | \
  xargs sed -i $i -e

  sed -i 's/<title><\/title>/<title>Collective Variables Module - Colvars Module - Reference Manual<\/title>/' $i

  # Remove <selfref></selfref> tags, where self-hyperlinks were added with
  # the sole purpose of generating tags suitable for use in URLs (htlatex
  # will only generate those for elements that are linked explicitly)
  python3 remove_htmltags.py $i selfref

  # Remove ligatures
  # https://orbythebeach.wordpress.com/2014/09/23/removing-ligatures-in-html-files-generated-from-latex-files/
  sed -i 's/\xef\xac\x80/ff/g' $i
  sed -i 's/\xef\xac\x81/fi/g' $i
  sed -i 's/\xef\xac\x82/fl/g' $i
  sed -i 's/\xef\xac\x83/ffi/g' $i
  sed -i 's/\xef\xac\x84/ffl/g' $i
  # sed -i 's/\xc5\x92/OE/g' $i
  # sed -i 's/\xc5\x93/oe/g' $i
  # sed -i 's/\xc3\x86/AE/g' $i
  # sed -i 's/\xc3\xa6/ae/g' $i
  # sed -i 's/\xef\xac\x86/st/g' $i
  # sed -i 's/\xc4\xb2/IJ/g' $i
  # sed -i 's/\xc4\xb3/ij/g' $i

  # Remove bogus spaces at first line of code blocks
  sed -i 's/<p class="noindent"><\/p><pre> /<p class="noindent"><\/p><pre>/' $i
  sed -i 's/<br class="newline"\/> /<br class="newline"\/>/g' $i

  # Fix characters in code blocks
  sed -i "s/’/'/g" $i
  sed -i 's/¿/\>/g' $i
  sed -i 's/¡/\</g' $i
  sed -i 's/”/"/g' $i
  sed -i 's/∖/\\/g' $i
  sed -i 's/–/--/g' $i

done
