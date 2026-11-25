#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Extract the configuration keywords parsed by the Colvars library source code.

This is a static scan of src/ for calls to the colvarparse functions:
  get_keyval(conf, "keyword", ...)
  get_keyval_feature(obj, conf, "keyword", ...)
  key_lookup(conf, "keyword", ...)
  parse_group(conf, "keyword", ...)          (atom groups of components)
and for the registration of block keywords:
  add_component_type<class>("description", "keyword")
  parse_biases_type<class>(conf, "keyword")

The context of each keyword (global, colvar, a component type, a bias type,
atom group...) is inferred from the class of the function that parses it; for a
base class, the context is the set of registered types derived from it.

Keywords whose name is computed at run time are reported with "dynamic": true
and the C++ expression in place of the name.

Output: a JSON file (default: colvars-keywords-src.json).
"""

import argparse
import json
import os
import re
import subprocess
import sys

SKIP_FILES = ('colvarparse.cpp', 'colvarparse.h', 'colvardeps.cpp', 'colvardeps.h')

CALLS = {
    # function name: (index of the keyword argument, index of the default, index of the mode)
    'get_keyval': (1, 3, 4),
    'get_keyval_feature': (2, 4, 5),
    'key_lookup': (1, None, None),
    'parse_group': (1, None, None),
}


def strip_cpp_comments(text):
    """Blank out comments, preserving line numbers and string literals."""
    pattern = re.compile(r'//[^\n]*|/\*.*?\*/|"(?:\\.|[^"\\])*"|\'(?:\\.|[^\'\\])*\'',
                         re.DOTALL)

    def repl(m):
        s = m.group(0)
        if s.startswith('/'):
            return re.sub(r'[^\n]', ' ', s)
        return s
    return pattern.sub(repl, text)


def split_call_args(text, start):
    """text[start] == '('; return (list_of_args, end_index) splitting at
    top-level commas."""
    depth = 0
    args = []
    current = []
    i = start
    while i < len(text):
        c = text[i]
        if c == '"':
            j = i + 1
            while text[j] != '"':
                j += 2 if text[j] == '\\' else 1
            current.append(text[i:j+1])
            i = j + 1
            continue
        if c in '([{':
            depth += 1
            if depth == 1:
                i += 1
                continue
        elif c in ')]}':
            depth -= 1
            if depth == 0:
                args.append(''.join(current).strip())
                return args, i
        elif c == ',' and depth == 1:
            args.append(''.join(current).strip())
            current = []
            i += 1
            continue
        current.append(c)
        i += 1
    raise ValueError("Unbalanced parentheses")


def normalize_class(name):
    name = re.sub(r'<[^<>]*>', '', name).strip()
    name = re.sub(r'^(cvm|colvarmodule)::', 'colvarmodule::', name)
    return name


def line_of(text, pos):
    return text.count('\n', 0, pos) + 1


class SourceScanner:

    def __init__(self, src_dir, files):
        self.src_dir = src_dir
        self.files = files
        self.parents = {}         # class -> [base classes]
        self.cvc_types = {}       # class -> component keyword
        self.bias_types = {}      # class -> bias keyword
        self.keywords = []

    def read(self, f):
        return strip_cpp_comments(open(os.path.join(self.src_dir, f)).read())

    def scan_classes(self):
        decl = re.compile(r'^\s*(?:class|struct)\s+([\w:]+)\s*(?::([^{;]*))?\{',
                          re.MULTILINE)
        for f in self.files:
            text = self.read(f)
            for m in decl.finditer(text):
                name = normalize_class(m.group(1))
                bases = []
                if m.group(2):
                    for b in m.group(2).split(','):
                        b = re.sub(r'\b(public|protected|private|virtual)\b', '', b)
                        bases.append(normalize_class(b))
                # Nested class names are qualified by the outer class (colvar::X)
                self.parents[name] = [self.qualify(b, name) for b in bases if b]
            for m in re.finditer(r'add_component_type\s*<\s*([\w:]+)\s*>\s*\(\s*"[^"]*"\s*,\s*"(\w+)"',
                                 text):
                self.cvc_types['colvar::' + m.group(1)] = m.group(2)
            for m in re.finditer(r'parse_biases_type\s*<\s*([\w:]+)\s*>\s*\(\s*\w+\s*,\s*"(\w+)"',
                                 text):
                self.bias_types[normalize_class(m.group(1))] = m.group(2)

    def qualify(self, base, derived):
        if '::' not in base and '::' in derived:
            outer = derived.rsplit('::', 1)[0]
            if outer + '::' + base in self.parents or base not in self.parents:
                return outer + '::' + base
        return base

    def ancestors(self, cls):
        result = set()
        stack = [cls]
        while stack:
            c = stack.pop()
            for p in self.parents.get(c, []):
                if p not in result:
                    result.add(p)
                    stack.append(p)
        return result

    def contexts_of_class(self, cls):
        """Map a C++ class to the configuration block(s) where its keywords go."""
        if cls is None:
            return ['unknown']
        if cls in ('colvarmodule', 'colvarproxy') or cls.startswith('colvarproxy'):
            return ['global']
        if cls == 'colvar':
            return ['colvar']
        if cls == 'colvarmodule::atom_group':
            return ['atom group']
        if cls.startswith('colvar_grid') or cls == 'colvarmodule::grid':
            return ['grid']
        types = dict(self.cvc_types)
        types.update(self.bias_types)
        # The type itself (if registered) and all registered types derived from it
        derived = sorted(v for k, v in types.items() if k == cls or cls in self.ancestors(k))
        if derived:
            return derived
        return ['class ' + cls]

    def enclosing_class(self, f, text, pos):
        """Class of the function definition (in .cpp) or class body (in .h)
        enclosing position pos."""
        head = text[:pos]
        defs = list(re.finditer(r'^[A-Za-z_][^\n;]*?([\w:<>]+)::~?(\w+)\s*\(',
                                head, re.MULTILINE))
        decls = []
        if f.endswith('.h'):
            decls = list(re.finditer(r'^\s*(?:class|struct)\s+([\w:]+)[^;]*?\{', head,
                                     re.MULTILINE))
        best = max(defs[-1:] + decls[-1:], key=lambda m: m.start(), default=None)
        if best is None:
            return None, None
        method = best.group(2) if best.re.groups >= 2 else None
        return normalize_class(best.group(1)), method

    def scan_calls(self):
        call_re = re.compile(r'\b(' + '|'.join(CALLS) + r')\s*\(')
        for f in self.files:
            if os.path.basename(f) in SKIP_FILES:
                continue
            text = self.read(f)
            for m in call_re.finditer(text):
                func = m.group(1)
                # Skip declarations/definitions of these functions
                prefix = text[text.rfind('\n', 0, m.start())+1:m.start()]
                if re.search(r'\b(bool|int|void|cvm::atom_group\s*\*|atom_group\s*\*)\s*(\w+::)*$',
                             prefix):
                    continue
                try:
                    args, end = split_call_args(text, m.end() - 1)
                except (ValueError, IndexError):
                    continue
                ikey, idef, imode = CALLS[func]
                if len(args) <= ikey:
                    continue
                key_expr = args[ikey]
                lit = re.fullmatch(r'"([^"]*)"', key_expr)
                cls, method = self.enclosing_class(f, text, m.start())
                entry = {
                    'name': lit.group(1) if lit else key_expr,
                    'contexts': self.contexts_of_class(cls),
                    'class': cls,
                    'method': method,
                    'function': func,
                    'file': 'src/' + f,
                    'line': line_of(text, m.start()),
                }
                if not lit:
                    entry['dynamic'] = True
                if idef is not None and len(args) > idef:
                    if 'parse_' in args[idef]:
                        entry['mode'] = args[idef]
                    else:
                        entry['default'] = args[idef]
                if imode is not None and len(args) > imode:
                    entry['mode'] = args[imode]
                if 'mode' in entry:
                    entry['mode'] = re.sub(r'(colvarparse::|Parse_Mode\()', '',
                                           entry['mode']).rstrip(')')
                # Keywords read from state files rather than from the configuration
                if method and re.search(r'state|restart', method, re.IGNORECASE):
                    entry['state_file'] = True
                # Legacy keywords recognized for backward compatibility
                if 'parse_deprecated' in entry.get('mode', '') or \
                   re.search(r'deprecated|legacy|obsolete', text[max(0, m.start()-300):end+300],
                             re.IGNORECASE):
                    entry['legacy_hint'] = True
                self.keywords.append(entry)

    def add_block_keywords(self):
        for cls, key in sorted(self.cvc_types.items(), key=lambda x: x[1]):
            self.keywords.append({'name': key, 'contexts': ['colvar'], 'class': cls,
                                  'function': 'add_component_type',
                                  'file': 'src/colvar.cpp', 'block': True})
        for cls, key in sorted(self.bias_types.items(), key=lambda x: x[1]):
            self.keywords.append({'name': key, 'contexts': ['global'], 'class': cls,
                                  'function': 'parse_biases_type',
                                  'file': 'src/colvarmodule.cpp', 'block': True})


def main():
    repo_default = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--repo', default=repo_default, help='Colvars repository root')
    parser.add_argument('--output', default='colvars-keywords-src.json')
    args = parser.parse_args()

    src_dir = os.path.join(args.repo, 'src')
    try:
        files = subprocess.check_output(['git', 'ls-files', '*.cpp', '*.h'], cwd=src_dir,
                                        text=True).split()
    except (OSError, subprocess.CalledProcessError):
        files = sorted(f for f in os.listdir(src_dir) if f.endswith(('.cpp', '.h')))
    files = [f for f in files if '/' not in f]   # skip src/cuda etc.

    scanner = SourceScanner(src_dir, files)
    scanner.scan_classes()
    scanner.scan_calls()
    scanner.add_block_keywords()
    out = {'component_types': dict(sorted((v, k) for k, v in scanner.cvc_types.items())),
           'bias_types': dict(sorted((v, k) for k, v in scanner.bias_types.items())),
           'keywords': scanner.keywords}
    with open(args.output, 'w') as f:
        json.dump(out, f, indent=1)
        f.write('\n')
    n_dyn = sum(1 for k in scanner.keywords if k.get('dynamic'))
    print("%d keyword uses (%d with run-time names), %d component types, %d bias types -> %s"
          % (len(scanner.keywords), n_dyn, len(scanner.cvc_types), len(scanner.bias_types),
             args.output))


if __name__ == '__main__':
    main()
