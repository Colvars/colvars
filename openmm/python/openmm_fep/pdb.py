from __future__ import annotations

from pathlib import Path


_PDB_COLS = {
    "X": (30, 38),
    "Y": (38, 46),
    "Z": (46, 54),
    "O": (54, 60),
    "B": (60, 66),
}


def parse_alch_file(path: str | Path, alch_col: str = "B") -> tuple[tuple[int, ...], tuple[int, ...]]:
    """Parse NAMD-style PDB perturbation flags.

    Returns ``(vanish, grow)`` as 0-based atom indices in file order.
    """

    col = alch_col.upper()
    if col not in _PDB_COLS:
        raise ValueError("alch_col must be one of X, Y, Z, O, B")
    start, end = _PDB_COLS[col]
    vanish: list[int] = []
    grow: list[int] = []
    atom_index = 0
    with Path(path).open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            raw = line[start:end].strip()
            flag = float(raw) if raw else 0.0
            if flag == -1.0:
                vanish.append(atom_index)
            elif flag == 1.0:
                grow.append(atom_index)
            elif flag != 0.0:
                raise ValueError(
                    f"unsupported alchemical flag {flag!r} at atom index {atom_index}; expected -1, 0, or 1"
                )
            atom_index += 1
    return tuple(vanish), tuple(grow)

