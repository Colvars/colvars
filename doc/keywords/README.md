# Colvars keyword references

Machine-readable lists of the Colvars configuration keywords, meant for AI tools,
input validators and editor plugins, and a consistency check between the manual
and the code. This is a first step; where this material lives long-term is still
to be decided.

All outputs are build products (not tracked by git). From `doc/`:

```
make keywords        # generate the files below
make check-keywords  # also write keywords-report.md (manual vs. code)
```

## Files

- `colvars-keywords-<engine>.json` (from `extract_doc_keywords.py`): every keyword
  documented in the manual of that engine, taken from the `\key`, `\keydef`,
  `\simkey` and `\dupkey` macros of `colvars-refman-main.tex`, after applying the
  engine-specific macros and flags (`colvars-refman-<engine>.tex`,
  `colvars-engine-flags.tex`). Each entry has `name`, `id` (the `\labelkey` path when
  available), `context` (the configuration block), `section`, `summary`, `type`,
  `default`, `description` (Markdown), and `same_as`/`similar_to` for keywords
  documented by reference to another one. The file also gives the engine's units.
- `colvars-keywords-<engine>.md`: the same content as Markdown, one file per engine,
  suitable to give to an AI tool.
- `colvars-keywords-src.json` (from `extract_src_keywords.py`): keywords parsed by
  the library (`src/`), from a static scan of `get_keyval`, `get_keyval_feature`,
  `key_lookup` and `parse_group` calls, plus the component and bias types
  registered with `add_component_type` and `parse_biases_type`. The context is
  inferred from the class of the parsing function (for a base class: all the
  registered types derived from it). Keywords whose name is built at run time
  are listed with `"dynamic": true`.
- `keywords-report.md` (from `compare_keywords.py`): keywords that the code parses
  but the manual does not document (split into configuration, apparently legacy,
  and state-file/silent keywords), documented keywords that the code does not
  parse, context mismatches, and literal defaults that differ. Names are compared
  case-insensitively. The report never fails the build.

## Limitations

- The LaTeX to Markdown conversion is approximate; math is left in LaTeX syntax.
- The source scan is static and regex-based: contexts may be wrong when a keyword
  is parsed by a helper function, and run-time keyword names are not compared.
- Engine-level options (NAMD/LAMMPS/GROMACS input keywords such as `colvarsConfig`)
  are parsed by the engine interfaces, not `src/`, and are excluded from the comparison.
