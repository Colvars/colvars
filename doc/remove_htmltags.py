#!/usr/bin/env python

from __future__ import print_function

import os
import sys
from bs4 import BeautifulSoup

if (len(sys.argv) < 3):
    print("Remove all given tags from the given HTML file")
    print("Usage: python3 remove_htmltag.py file.html tag1 [tag2 ...]")
    sys.exit(1)

html_file = open(sys.argv[1], 'r')
html = html_file.read()
html_file.close()

soup = BeautifulSoup(html, 'lxml')

for tag in sys.argv[2:]:
    for ref in soup(tag):
        ref.decompose()

html_file = open(sys.argv[1], 'w')
html_file.write(str(soup))
html_file.close()

