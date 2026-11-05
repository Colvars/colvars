#!/usr/bin/env python3
"""Wrap stray MathML text in HTML files. Requires: pip install beautifulsoup4.

Default: write NAME.fixed.html; --in-place: save NAME.html.bak first.
Input/output encoding defaults to UTF-8; use --encoding for other encodings.
This repairs token wrapping, not arbitrary MathML syntax or math semantics.
Beautiful Soup can normalize markup elsewhere in files that need repairs.
"""

import argparse
import re
import sys
from pathlib import Path

from bs4 import BeautifulSoup, NavigableString, Tag


# Presentation containers that accept token elements as direct children.
CONTAINERS = set("""
math mrow mfrac msqrt mroot mstyle merror mpadded mphantom mfenced
menclose msub msup msubsup munder mover munderover mmultiscripts
mtd maction semantics
""".split())
PROTECTED = {
    "mi", "mn", "mo", "ms", "mtext", "annotation", "annotation-xml",
    "ci", "cn", "cs", "csymbol", "script", "style",
}
RESTRICTED = {
    "mtable", "mtr", "mlabeledtr", "mspace", "mglyph", "malignmark",
    "maligngroup", "mprescripts", "none", "mstack", "mlongdiv",
    "msgroup", "msrow", "mscarries", "mscarry", "msline",
}
ARITY = {
    "mfrac": 2, "mroot": 2, "msub": 2, "msup": 2,
    "msubsup": 3, "munder": 2, "mover": 2, "munderover": 3,
}
XML_SPACE = " \t\r\n"
# Conservative recognition: decimal numbers, optional sign and exponent.
# Other forms (e.g. locale-specific separators) are preserved as mtext.
NUMBER = re.compile(r"[+\-\u2212]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+\-\u2212]?[0-9]+)?")


def local_name(tag):
    return tag.name.rsplit(":", 1)[-1].lower()


def in_math(tag):
    """Ignore HTML, annotation payloads, and anything inside existing tokens."""
    found_math = False
    for ancestor in [tag, *tag.parents]:
        if not isinstance(ancestor, Tag):
            continue
        name = local_name(ancestor)
        if name in PROTECTED:
            return False
        if name == "math":
            found_math = True
    return found_math


def fix_mathml(html):
    """Return (HTML, number of wrapped nodes, review warnings).

    One text node becomes one token; never guess how to split expressions
    or fraction/script arguments. Whitespace-only nodes and comments stay
    untouched. Already-valid files are returned byte-for-byte at string level.
    """
    soup = BeautifulSoup(html, "html.parser")
    count = 0
    warnings = []
    for parent in list(soup.find_all(True)):
        name = local_name(parent)
        if not in_math(parent):
            continue
        for child in list(parent.children):
            # Comments, declarations, etc. subclass NavigableString; exclude them.
            if type(child) is not NavigableString:
                continue
            value = str(child).strip(XML_SPACE)
            if not value:
                continue
            if name in RESTRICTED:
                warnings.append(f"<{parent.name}> has raw text {value[:60]!r}; needs structural repair")
                continue
            if name not in CONTAINERS:
                continue
            token_name = "mn" if NUMBER.fullmatch(value) else "mtext"
            # Preserve an explicit MathML prefix, e.g. <m:msup> -> <m:mn>.
            if ":" in parent.name:
                token_name = parent.name.rsplit(":", 1)[0] + ":" + token_name
            token = soup.new_tag(token_name)
            child.wrap(token)
            count += 1
        if name in ARITY:
            actual = sum(isinstance(child, Tag) for child in parent.children)
            if actual != ARITY[name]:
                warnings.append(f"<{parent.name}> has {actual} element children; expected {ARITY[name]}")
    return (soup.decode(formatter="minimal") if count else html), count, warnings


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("files", nargs="+", type=Path, help="HTML files (shell wildcards supported)")
    parser.add_argument("--in-place", action="store_true", help="edit originals; create .bak files without overwriting backups")
    parser.add_argument("--dry-run", action="store_true", help="report repairs without writing")
    parser.add_argument("--encoding", default="utf-8", help="file encoding (default: utf-8)")
    args = parser.parse_args()
    status = 0
    seen = set()
    for path in args.files:
        if path.resolve() in seen:
            continue
        seen.add(path.resolve())
        try:
            original = path.read_bytes()
            repaired, count, warnings = fix_mathml(original.decode(args.encoding))
            for warning in warnings:
                print(f"{path}: REVIEW: {warning}", file=sys.stderr)
            if warnings:
                status = max(status, 2)
            if args.dry_run or not count:
                print(f"{path}: {count} text node(s) to wrap")
                continue
            # Encode before opening a destination or creating a backup.
            output = repaired.encode(args.encoding)
            if args.in_place:
                backup = path.with_name(path.name + ".bak")
                with backup.open("xb") as stream:
                    stream.write(original)
                path.write_bytes(output)
                destination = path
            else:
                destination = path.with_name(path.stem + ".fixed" + path.suffix)
                with destination.open("xb") as stream:
                    stream.write(output)
            print(f"{path}: wrapped {count} text node(s) -> {destination}")
        except (OSError, UnicodeError, LookupError) as error:
            print(f"{path}: ERROR: {error}", file=sys.stderr)
            status = 1
    return status


if __name__ == "__main__":
    sys.exit(main())
