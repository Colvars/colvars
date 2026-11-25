#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Compare the keywords documented in the manual with those parsed by src/.

Inputs are the outputs of extract_doc_keywords.py (one JSON per engine) and
extract_src_keywords.py. Writes a Markdown report listing:
  1. keywords parsed by the code but not documented in any engine's manual;
  2. documented keywords that the library code does not parse;
  3. documented keywords whose context differs from where the code parses them;
  4. defaults that differ between the manual and the code (literal values only);
  5. deprecated keywords that are still documented.
Keyword names are compared case-insensitively, as Colvars does.

This is a report only: it does not fail on differences.
"""

import argparse
import glob
import json
import os
import re
import sys


def split_doc_context(context):
    """'gspath and gzpath', 'angle, dipoleAngle, dihedral' -> list of names."""
    parts = re.split(r',|\band\b', context)
    return [p.strip() for p in parts if p.strip()]


def normalize_value(text):
    """Return a comparable Python value from a default (doc or C++), or None."""
    if text is None:
        return None
    t = text.strip().strip('`').strip()
    t = re.sub(r'^std::string\s*\(\s*(.*)\s*\)$', r'\1', t)
    t = re.sub(r'^cvm::real\s*\(\s*(.*)\s*\)$', r'\1', t)
    low = t.lower()
    if low in ('on', 'yes', 'true'):
        return True
    if low in ('off', 'no', 'false'):
        return False
    m = re.fullmatch(r'"([^"]*)"', t)
    if m:
        return ('str', m.group(1))
    m = re.fullmatch(r'([-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?)[fFlLuU]*', t)
    if m:
        return float(m.group(1))
    return None


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--doc-json', nargs='+',
                        help='colvars-keywords-<engine>.json files (default: all in --dir)')
    parser.add_argument('--src-json', help='colvars-keywords-src.json (default: in --dir)')
    parser.add_argument('--dir', default='.', help='Directory of the JSON files')
    parser.add_argument('--output', default='-', help='Report file (default: stdout)')
    args = parser.parse_args()

    doc_files = args.doc_json or sorted(
        f for f in glob.glob(os.path.join(args.dir, 'colvars-keywords-*.json'))
        if not f.endswith('-src.json'))
    src_file = args.src_json or os.path.join(args.dir, 'colvars-keywords-src.json')
    if not doc_files or not os.path.exists(src_file):
        sys.exit("Missing input files: run extract_doc_keywords.py and extract_src_keywords.py")

    src = json.load(open(src_file))
    all_cvcs = set(c.lower() for c in src['component_types'])
    all_biases = set(b.lower() for b in src['bias_types'])

    def expand(ctx):
        c = ctx.lower()
        if c in ('any component', 'components with optional periodicity'):
            return all_cvcs
        if c == 'periodic components':
            return all_cvcs | {'colvar'}
        if c == 'colvar bias':
            return all_biases
        return {c}

    # Documented keywords: lowercase name -> list of (entry, engine)
    doc = {}
    engine_level = set()   # names documented as engine (not library) options
    engines = []
    section_names = set()  # `names` in section titles (component types...)
    for f in doc_files:
        ref = json.load(open(f))
        for title in ref.get('sections', {}).values():
            section_names |= set(n.lower() for n in re.findall(r'`([^`]+)`', title))
        engine = ref['engine']
        engines.append(engine)
        for k in ref['keywords']:
            ctx = k['context'].lower()
            if engine.lower() in ctx or re.search(r'configuration file|mdp option|fix colvars', ctx):
                engine_level.add(k['name'].lower())
                continue
            doc.setdefault(k['name'].lower(), []).append((k, engine))

    # Source keywords: lowercase name -> list of entries
    code = {}
    for k in src['keywords']:
        if k.get('dynamic'):
            continue
        code.setdefault(k['name'].lower(), []).append(k)

    # Names of configuration blocks that the manual documents in section titles
    # or as contexts of other keywords (component and bias types, colvar...)
    doc_contexts = set(section_names)
    for entries in doc.values():
        for k, _ in entries:
            doc_contexts |= set(c.lower() for c in split_doc_context(k['context']))

    def category(entries):
        """'state' (state files/silent), 'legacy', or None (configuration)."""
        if all(e.get('state_file') or re.search(r'parse_(silent|restart)', e.get('mode', ''))
               for e in entries):
            return 'state'
        if all(e.get('legacy_hint') for e in entries):
            return 'legacy'
        return None

    def loc(e):
        return '%s:%s' % (e['file'], e['line']) if 'line' in e else e['file']

    out = []
    w = out.append
    w('# Colvars keywords: manual vs. source code')
    w('')
    w('Manuals compared: %s. Names compared case-insensitively.' % ', '.join(engines))
    w('')

    # 1. Undocumented
    rows = {None: [], 'legacy': [], 'state': []}
    for name, entries in sorted(code.items()):
        if name in doc or name in engine_level:
            continue
        if name in doc_contexts:   # block keyword (colvar, grid, harmonic...)
            continue
        rows[category(entries)].append('| `%s` | %s | %s |' % (
            entries[0]['name'],
            ', '.join(sorted(set(c for e in entries for c in e['contexts'])))[:80],
            ', '.join(sorted(set(loc(e) for e in entries)))))
    headers = {
        None: '## 1. Parsed by the code, not documented (%d)',
        'legacy': '### 1b. Not documented, apparently legacy/deprecated (%d)',
        'state': '### 1c. Not documented, read only from state files or silently (%d)',
    }
    for cat in (None, 'legacy', 'state'):
        w(headers[cat] % len(rows[cat]))
        w('')
        if rows[cat]:
            w('| Keyword | Context(s) in code | Location |')
            w('|---|---|---|')
            out.extend(rows[cat])
            w('')

    # 2. Documented but not parsed
    rows = []
    for name, entries in sorted(doc.items()):
        if name in code:
            continue
        k = entries[0][0]
        rows.append('| `%s` | %s | %s | %s |' % (
            k['name'], k['context'], k.get('section', ''),
            ', '.join(sorted(set(e for _, e in entries)))))
    w('## 2. Documented, not found in the library code (%d)' % len(rows))
    w('')
    w('May include keywords parsed with a run-time name (see the "dynamic" entries in the '
      'source JSON) or parsed by engine interface code.')
    w('')
    if rows:
        w('| Keyword | Context in manual | Section | Engines |')
        w('|---|---|---|---|')
        out.extend(rows)
        w('')

    # 3. Context mismatches
    rows = []
    for name, entries in sorted(doc.items()):
        if name not in code:
            continue
        code_ctx = set(c.lower() for e in code[name] for c in e['contexts'])
        seen = set()
        for k, engine in entries:
            if k['context'] in seen:
                continue
            seen.add(k['context'])
            doc_ctx = set()
            for c in split_doc_context(k['context']):
                doc_ctx |= expand(c)
            missing = sorted(doc_ctx - code_ctx)
            if missing:
                rows.append('| `%s` | %s | %s | %s |' % (
                    k['name'], k['context'], ', '.join(missing)[:80],
                    ', '.join(sorted(set(loc(e) for e in code[name])))[:120]))
    w('## 3. Documented context not matched by the code (%d)' % len(rows))
    w('')
    w('Contexts of the code are inferred from the class of the parsing function and its '
      'registered derived types; mismatches may be false positives.')
    w('')
    if rows:
        w('| Keyword | Context in manual | Not parsed by the code in | Code location |')
        w('|---|---|---|---|')
        out.extend(rows)
        w('')

    # 4. Default mismatches
    rows = []
    for name, entries in sorted(doc.items()):
        if name not in code:
            continue
        for k, engine in entries:
            dval = normalize_value(k.get('default'))
            if dval is None or k.get('kind') not in ('key', 'keydef'):
                continue
            doc_ctx = set()
            for c in split_doc_context(k['context']):
                doc_ctx |= expand(c)
            candidates = [e for e in code[name]
                          if doc_ctx & set(c.lower() for c in e['contexts'])
                          and normalize_value(e.get('default')) is not None]
            if not candidates:
                continue
            cvals = set(normalize_value(e['default']) for e in candidates)
            if dval not in cvals:
                rows.append('| `%s` | %s | %s | %s | %s |' % (
                    k['name'], k['context'], k['default'],
                    ', '.join('`%s`' % e['default'] for e in candidates),
                    ', '.join(loc(e) for e in candidates)))
            break   # defaults are the same in all engines' manuals
    w('## 4. Literal defaults that differ (%d)' % len(rows))
    w('')
    if rows:
        w('| Keyword | Context | Default in manual | Default in code | Location |')
        w('|---|---|---|---|---|')
        out.extend(rows)
        w('')

    # 5. Deprecated but documented
    rows = []
    for name, entries in sorted(code.items()):
        dep = [e for e in entries if 'parse_deprecated' in e.get('mode', '')]
        if dep and name in doc:
            k = doc[name][0][0]
            rows.append('| `%s` | %s | %s |' % (k['name'], k.get('section', ''),
                                              ', '.join(loc(e) for e in dep)))
    w('## 5. Deprecated in the code, still documented (%d)' % len(rows))
    w('')
    if rows:
        w('| Keyword | Section in manual | Location |')
        w('|---|---|---|')
        out.extend(rows)
        w('')

    dyn = [k for k in src['keywords'] if k.get('dynamic')]
    w('## Keywords with run-time names (not compared): %d' % len(dyn))
    w('')
    for k in dyn:
        w('- `%s` (%s)' % (k['name'], loc(k)))
    w('')

    text = '\n'.join(out)
    if args.output == '-':
        print(text)
    else:
        with open(args.output, 'w') as f:
            f.write(text)


if __name__ == '__main__':
    main()
