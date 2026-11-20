#!/usr/bin/env python3
"""Post-process the index page (idx.html) produced by LaTeXML.

LaTeXML labels each index reference by the number of the list item that
contains it ("12nd item"), which is not informative.  Put the link directly
on the indexed phrase (pointing to its first occurrence).  When a phrase is
indexed in more than one place, list all occurrences labeled by their
section numbers instead.

Usage: fix_index.py idx.html [...]
"""

import re
import sys

# One index entry: phrase followed by its references (LaTeXML 0.8.8 output)
ENTRY = re.compile(
    r'(<span class="ltx_indexphrase">)(.*?)(</span>)'
    r'<span class="ltx_indexrefs">(.*?)</span>\n')

# One reference: link with a title "In <section number> <title> ..."
REF = re.compile(
    r'<a href="([^"]*)" title="([^"]*)" class="ltx_ref">'
    r'<span class="ltx_text ltx_ref_tag">[^<]*</span></a>')

SECTION_NUMBER = re.compile(r'^In ([0-9A-Z][0-9.]*) ')


def fix_entry(match):
    open_tag, phrase, close_tag, refs_html = match.groups()
    refs = REF.findall(refs_html)
    # Leave alone anything not in the expected form (e.g. "see also")
    if not refs or len(REF.sub('', refs_html).replace(
            '<span class="ltx_text"> </span>', '').strip(' ,')) > 0:
        return match.group(0)
    href, title = refs[0]
    linked = (f'{open_tag}<a href="{href}" title="{title}" class="ltx_ref">'
              f'{phrase}</a>{close_tag}')
    if len(refs) == 1:
        return linked + '\n'
    labels = []
    for href, title in refs:
        number = SECTION_NUMBER.match(title)
        label = number.group(1) if number else 'link'
        labels.append(f'<a href="{href}" title="{title}" class="ltx_ref">'
                      f'{label}</a>')
    return (linked + '<span class="ltx_indexrefs"> (' + ', '.join(labels) +
            ')</span>\n')


def main():
    for filename in sys.argv[1:]:
        with open(filename, encoding='utf-8') as f:
            html = f.read()
        html, count = ENTRY.subn(fix_entry, html)
        with open(filename, 'w', encoding='utf-8') as f:
            f.write(html)
        print(f'{filename}: processed {count} index entries')


if __name__ == '__main__':
    main()
