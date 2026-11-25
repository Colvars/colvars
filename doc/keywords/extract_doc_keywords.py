#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Extract the Colvars configuration keywords documented in the reference manual.

The manual (colvars-refman-main.tex) documents each keyword with one of these
macros (defined in colvars-refman.tex):
  \\key{name}{context}{summary}{type}{description}
  \\keydef{name}{context}{summary}{type}{default}{description}
  \\simkey{name}{context}{similar_keyword}
  \\dupkey{name}{context}{target_label}{where_documented}
optionally preceded by \\labelkey{path|to|keyword}.

For each engine, the engine-specific macros of colvars-refman-<engine>.tex and
the flags of colvars-engine-flags.tex are applied first, so that each
reference only contains what that engine's manual contains.

Outputs (in the output directory):
  colvars-keywords-<engine>.json  machine-readable reference
  colvars-keywords-<engine>.md    the same, rendered as Markdown
"""

import argparse
import json
import os
import re
import sys

ENGINES = ['gromacs', 'lammps', 'namd', 'tinkerhp', 'vmd']


# ---------------------------------------------------------------------------
# Low-level LaTeX helpers
# ---------------------------------------------------------------------------

def strip_comments(text):
    """Remove LaTeX comments; like TeX, a comment also eats the newline and the
    leading spaces of the next line."""
    return re.sub(r'(?<!\\)%[^\n]*\n[ \t]*', '', text + '\n')


def match_brace(text, start):
    """Given text[start] == '{', return the index of the matching '}'."""
    assert text[start] == '{'
    depth = 0
    i = start
    while i < len(text):
        c = text[i]
        if c == '\\':
            i += 2
            continue
        if c == '{':
            depth += 1
        elif c == '}':
            depth -= 1
            if depth == 0:
                return i
        i += 1
    raise ValueError("Unbalanced braces starting at: " + text[start:start+80])


def read_args(text, pos, nargs, optional=False):
    """Read nargs brace-delimited arguments starting at pos (whitespace allowed
    in between); an optional [..] argument is skipped if optional is True.
    Returns (list_of_args, end_position)."""
    args = []
    if optional:
        m = re.compile(r'\s*\[([^\]]*)\]').match(text, pos)
        if m:
            pos = m.end()
    for _ in range(nargs):
        m = re.compile(r'\s*').match(text, pos)
        pos = m.end()
        if pos >= len(text) or text[pos] != '{':
            raise ValueError("Expected '{' at: " + text[pos:pos+80])
        end = match_brace(text, pos)
        args.append(text[pos+1:end])
        pos = end + 1
    return args, pos


def replace_macro(text, name, nargs, func, optional=False):
    """Replace every \\name{arg1}...{argN} with func(args), innermost content
    first (func receives arguments that were already processed)."""
    pattern = re.compile(r'\\' + name + r'(?![a-zA-Z])')
    out = []
    pos = 0
    while True:
        m = pattern.search(text, pos)
        if not m:
            out.append(text[pos:])
            break
        out.append(text[pos:m.start()])
        try:
            args, end = read_args(text, m.end(), nargs, optional)
        except ValueError:
            # Not used as a macro with arguments here; leave it alone
            out.append(m.group(0))
            pos = m.end()
            continue
        args = [replace_macro(a, name, nargs, func, optional) for a in args]
        out.append(func(*args))
        pos = end
    return ''.join(out)


# ---------------------------------------------------------------------------
# Engine-specific preprocessing
# ---------------------------------------------------------------------------

class EngineConfig:
    """Macros and flags defined for one engine's manual."""

    def __init__(self, doc_dir, engine):
        self.engine = engine
        self.only_macros = {}   # \cvXonly -> True (keep) / False (drop)
        self.text_macros = {}   # \MDENGINE -> 'GROMACS', ...
        self.flags = set()      # names \def'ed (tested with \ifdefined)
        wrapper = os.path.join(doc_dir, 'colvars-refman-%s.tex' % engine)
        self.parse_newcommands(open(wrapper).read())
        flags_text = open(os.path.join(doc_dir, 'colvars-engine-flags.tex')).read()
        flags_text = self.apply_only_macros(strip_comments(flags_text))
        self.parse_newcommands(flags_text)
        for m in re.finditer(r'\\def\\([a-zA-Z]+)', flags_text):
            self.flags.add(m.group(1))
        version_file = os.path.join(doc_dir, 'cv_version.tex')
        self.parse_newcommands(open(version_file).read())

    def parse_newcommands(self, text):
        for m in re.finditer(r'\\newcommand\{\\([a-zA-Z]+)\}(\[1\])?\{(.*)\}\s*$',
                             text, re.MULTILINE):
            name, onearg, body = m.groups()
            if onearg:
                self.only_macros[name] = (body.strip() == '#1')
            else:
                self.text_macros[name] = body

    def apply_only_macros(self, text):
        for name, keep in self.only_macros.items():
            text = replace_macro(text, name, 1,
                                 (lambda x: x) if keep else (lambda x: ''))
        return text

    def resolve_conditionals(self, text):
        """Resolve \\ifdefined\\flag ... [\\else ...] \\fi and \\iflatexml."""
        tok = re.compile(r'\\(ifdefined\\([a-zA-Z]+)|iflatexml|else|fi)(?![a-zA-Z])')
        m = tok.search(text)
        while m and not m.group(1).startswith('if'):
            m = tok.search(text, m.end())   # stray \else/\fi: skip
        if not m:
            return text
        before = text[:m.start()]
        cond = (m.group(2) in self.flags) if m.group(2) else False  # not LaTeXML
        depth = 0
        else_pos = None
        pos = m.end()
        while True:
            n = tok.search(text, pos)
            if not n:
                raise ValueError("Unterminated conditional at: " + text[m.start():m.start()+80])
            word = n.group(1)
            if word.startswith('if'):
                depth += 1
            elif word == 'else' and depth == 0:
                else_pos = n
            elif word == 'fi':
                if depth == 0:
                    break
                depth -= 1
            pos = n.end()
        if else_pos:
            true_branch = text[m.end():else_pos.start()]
            false_branch = text[else_pos.end():n.start()]
        else:
            true_branch = text[m.end():n.start()]
            false_branch = ''
        chosen = true_branch if cond else false_branch
        return (before + self.resolve_conditionals(chosen) +
                self.resolve_conditionals(text[n.end():]))

    def apply_text_macros(self, text):
        for name, body in self.text_macros.items():
            text = re.sub(r'\\' + name + r'(?![a-zA-Z])(\{\})?',
                          lambda m: body.replace('\\', '\\\\'), text)
            # (re.sub would interpret backslashes in body; escaped above)
        return text

    def preprocess(self, text):
        text = strip_comments(text)
        text = self.resolve_conditionals(text)
        text = self.apply_only_macros(text)
        text = self.apply_text_macros(text)
        return text


# ---------------------------------------------------------------------------
# LaTeX to Markdown (approximate, for keyword descriptions)
# ---------------------------------------------------------------------------

def tex_to_md(text, labels=None):
    labels = labels or {}

    def ref(label):
        title = labels.get(label)
        return ('section "%s"' % title) if title else 'see manual'

    text = re.sub(r'\\(label|labelkey)\{[^}]*\}', '', text)
    text = replace_macro(text, 'refkey', 2, lambda a, b: '`%s`' % strip_md(a))
    text = replace_macro(text, 'texttt', 1, lambda a: '`%s`' % strip_md(a))
    text = replace_macro(text, 'mathtt', 1, lambda a: a)
    text = replace_macro(text, 'emph', 1, lambda a: '*%s*' % a.strip())
    text = replace_macro(text, 'textit', 1, lambda a: '*%s*' % a.strip())
    text = replace_macro(text, 'textbf', 1, lambda a: '**%s**' % a.strip())
    text = replace_macro(text, 'href', 2, lambda u, t: '[%s](%s)' % (t, u))
    text = replace_macro(text, 'url', 1, lambda u: '<%s>' % u)
    text = replace_macro(text, 'hyperref', 1, lambda t: t, optional=True)
    text = replace_macro(text, 'hyperlink', 2, lambda l, t: t)
    text = replace_macro(text, 'ref', 1, ref)
    for cite in ('cite', 'citealp', 'citep', 'citet'):
        text = replace_macro(text, cite, 1, lambda k: '[%s]' % k.replace(',', ', '))
    text = re.sub(r'\\begin\{(equation|align)\*?\}(.*?)\\end\{\1\*?\}',
                  lambda m: '$$' + m.group(2).strip() + '$$', text, flags=re.DOTALL)
    text = re.sub(r'\{\\ttfamily\s*([^{}]*)\}', r'`\1`', text)
    subs = [
        (r'\\begin\{(itemize|enumerate|cvcoptions|description)\}', '\n'),
        (r'\\end\{(itemize|enumerate|cvcoptions|description)\}', '\n'),
        (r'\s*\\item(?![a-zA-Z])\s*', '\n- '),
        (r'\\(noindent|newline|smallskip|medskip|bigskip|vspace\{[^}]*\})(\{\})?', ''),
        (r'\\\\', '\n'),
        (r'\\textunderscore(\{\})?', '_'),
        (r'\\textbackslash(\{\})?', '\\\\'),
        (r'\\AA(\{\})?', 'Å'),
        (r'\\ldots(\{\})?', '...'),
        (r'\{\\?_\}', '_'),
        (r'\\([_%&#$])', r'\1'),
        (r"``|''", '"'),
        (r'(?<!\\)~', ' '),
        (r'\\ ', ' '),
        (r'\{\}', ''),
    ]
    for pat, rep in subs:
        text = re.sub(pat, rep, text)
    # Normalize whitespace, keeping paragraph breaks and list items
    lines = [re.sub(r'[ \t]+', ' ', l).strip() for l in text.split('\n')]
    out = []
    for l in lines:
        if not l:
            if out and out[-1] != '':
                out.append('')
        elif l.startswith('- ') or not out or out[-1] == '':
            out.append(l)
        else:
            out[-1] += ' ' + l
    return '\n'.join(out).strip()


def strip_md(text):
    """Plain text from already-converted Markdown inline code fragments."""
    return text.replace('`', '').strip()


def clean_inline(text, labels=None):
    return re.sub(r'\s+', ' ', tex_to_md(text, labels)).strip()


# ---------------------------------------------------------------------------
# Keyword extraction
# ---------------------------------------------------------------------------

def extract_keywords(text):
    """Return (keywords, sections) from preprocessed manual text."""
    tok = re.compile(r'\\(cvsec|cvsubsec|cvsubsubsec|labelkey|keydef|key|simkey|dupkey|item)(?![a-zA-Z])'
                     r'|\\end\{cvcoptions\}')
    nargs = {'cvsec': 2, 'cvsubsec': 2, 'cvsubsubsec': 2, 'labelkey': 1,
             'keydef': 6, 'key': 5, 'simkey': 3, 'dupkey': 4, 'item': 0}
    sections = {}          # label -> title
    section_stack = []     # [(level, title, label)]
    levels = {'cvsec': 1, 'cvsubsec': 2, 'cvsubsubsec': 3}
    keywords = []
    pending_label = None
    last_context = None
    pos = 0
    while True:
        m = tok.search(text, pos)
        if not m:
            break
        if m.group(1) is None:   # \end{cvcoptions} implicitly documents forceNoPBC
            if last_context:
                keywords.append(make_entry('forceNoPBC', last_context, section_stack,
                                           same_as='colvar|cvc|forceNoPBC',
                                           kind='dupkey'))
            pos = m.end()
            continue
        macro = m.group(1)
        try:
            args, end = read_args(text, m.end(), nargs[macro])
        except ValueError as e:
            sys.stderr.write("Warning: cannot parse \\%s: %s\n" % (macro, e))
            pos = m.end()
            continue
        pos = end
        if macro in levels:
            title, label = args
            level = levels[macro]
            while section_stack and section_stack[-1][0] >= level:
                section_stack.pop()
            section_stack.append((level, title, label))
            sections[label] = title
        elif macro == 'item':
            pending_label = None
        elif macro == 'labelkey':
            pending_label = args[0].strip()
        else:
            if macro == 'keydef':
                name, context, summary, ktype, default, desc = args
                entry = make_entry(name, context, section_stack, summary=summary,
                                   type=ktype, default=default, description=desc)
            elif macro == 'key':
                name, context, summary, ktype, desc = args
                entry = make_entry(name, context, section_stack, summary=summary,
                                   type=ktype, description=desc)
            elif macro == 'simkey':
                name, context, similar = args
                entry = make_entry(name, context, section_stack, similar_to=similar)
            else:  # dupkey
                name, context, target, where = args
                entry = make_entry(name, context, section_stack, same_as=target.strip(),
                                   documented_in=where)
            entry['kind'] = macro
            if pending_label:
                entry['id'] = pending_label
                pending_label = None
            keywords.append(entry)
            last_context = context
    return keywords, sections


def make_entry(name, context, section_stack, kind=None, **fields):
    entry = {'name': name.strip(), 'context_tex': context}
    if section_stack:
        entry['section_label'] = section_stack[-1][2]
    entry.update(fields)
    if kind:
        entry['kind'] = kind
    return entry


def finalize(keywords, sections):
    """Convert LaTeX fields to Markdown and resolve cross-references."""
    by_id = {}
    for k in keywords:
        k['context'] = strip_md(clean_inline(k.pop('context_tex'), sections))
        k['name'] = strip_md(clean_inline(k['name']))
        if 'id' not in k:
            k['id'] = k['context'].replace(' ', '') + '|' + k['name']
        if 'section_label' in k:
            k['section'] = clean_inline(sections.get(k['section_label'], ''), sections)
        for f in ('summary', 'type', 'default', 'similar_to', 'documented_in'):
            if f in k:
                k[f] = clean_inline(k[f], sections)
        if 'description' in k:
            k['description'] = tex_to_md(k['description'], sections)
        by_id.setdefault(k['id'], k)
    # Keywords defined by reference to another one inherit its fields
    by_name = {}
    for k in keywords:
        if k.get('kind') in ('key', 'keydef'):
            by_name.setdefault(k['name'], []).append(k)
    for k in keywords:
        if 'same_as' in k:
            target = by_id.get(k['same_as'])
            if target is None:
                # The reference may be a section label: look up the same name,
                # preferably in that section or its subsections
                candidates = by_name.get(k['name'], [])
                in_section = [c for c in candidates
                              if c.get('section_label', '').startswith(k['same_as'])]
                if in_section:
                    target = in_section[0]
                elif len(candidates) == 1:
                    target = candidates[0]
                if target:
                    k['same_as'] = target['id']
            if target is None:
                sys.stderr.write("Warning: \\dupkey target %s not found (keyword %s)\n"
                                 % (k['same_as'], k['name']))
                continue
            for f in ('summary', 'type', 'default'):
                if f in target and f not in k:
                    k[f] = target[f]
    ordered = ['name', 'id', 'context', 'section', 'section_label', 'kind', 'summary',
               'type', 'default', 'similar_to', 'same_as', 'documented_in', 'description']
    return [{f: k[f] for f in ordered if f in k} for k in keywords]


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def write_markdown(ref, path):
    lines = []
    lines.append('# Colvars configuration keywords (%s)' % ref['engine'])
    lines.append('')
    lines.append('Generated from the Colvars reference manual (version %s). '
                 'Full manual: %s' % (ref['colvars_version'], ref['manual_url']))
    lines.append('')
    lines.append('Units in %s: length = %s, energy = %s.'
                 % (ref['engine'], ref['units']['length'], ref['units']['energy']))
    lines.append('')
    current_section = None
    for k in ref['keywords']:
        if k.get('section') != current_section:
            current_section = k.get('section')
            lines.append('## %s' % current_section)
            lines.append('')
        lines.append('### `%s` (context: %s)' % (k['name'], k['context']))
        if 'summary' in k:
            lines.append('*%s*' % k['summary'])
        meta = []
        if 'type' in k:
            meta.append('type: %s' % k['type'])
        if 'default' in k:
            meta.append('default: %s' % k['default'])
        if meta:
            lines.append('')
            lines.append('[' + '; '.join(meta) + ']')
        if 'similar_to' in k:
            lines.append('')
            lines.append('Analogous to `%s`.' % k['similar_to'])
        if 'same_as' in k:
            lines.append('')
            lines.append('Same definition as `%s`%s.' % (
                k['same_as'], (' (%s)' % k['documented_in']) if 'documented_in' in k else ''))
        if k.get('description'):
            lines.append('')
            lines.append(k['description'])
        lines.append('')
    with open(path, 'w') as f:
        f.write('\n'.join(lines))


def main():
    doc_dir_default = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--doc-dir', default=doc_dir_default,
                        help='Directory containing colvars-refman-main.tex')
    parser.add_argument('--output-dir', default='.', help='Where to write the outputs')
    parser.add_argument('--engines', nargs='+', default=ENGINES, choices=ENGINES)
    parser.add_argument('--manual-url', default='https://colvars.github.io/master/html-%s/',
                        help='URL pattern of the HTML manual (%%s = engine)')
    args = parser.parse_args()

    main_text = open(os.path.join(args.doc_dir, 'colvars-refman-main.tex')).read()
    os.makedirs(args.output_dir, exist_ok=True)
    for engine in args.engines:
        cfg = EngineConfig(args.doc_dir, engine)
        text = cfg.preprocess(main_text)
        keywords, sections = extract_keywords(text)
        keywords = finalize(keywords, sections)
        ref = {
            'engine': cfg.text_macros.get('MDENGINE', engine),
            'colvars_version': cfg.text_macros.get('cvversion', ''),
            'manual_url': args.manual_url % engine,
            'units': {'length': tex_to_md(cfg.text_macros.get('lengthunit', '')),
                      'energy': tex_to_md(cfg.text_macros.get('energyunit', ''))},
            'sections': {label: clean_inline(title) for label, title in sections.items()},
            'keywords': keywords,
        }
        base = os.path.join(args.output_dir, 'colvars-keywords-%s' % engine)
        with open(base + '.json', 'w') as f:
            json.dump(ref, f, indent=1, ensure_ascii=False)
            f.write('\n')
        write_markdown(ref, base + '.md')
        print("%s: %d keywords -> %s.{json,md}" % (engine, len(keywords), base))


if __name__ == '__main__':
    main()
