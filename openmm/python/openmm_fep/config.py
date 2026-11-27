from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Sequence

from .pdb import parse_alch_file
from .schedule import LambdaState, lambda_state


def _unique_ints(values: Iterable[int], name: str) -> tuple[int, ...]:
    out = tuple(int(v) for v in values)
    if len(set(out)) != len(out):
        raise ValueError(f"{name} contains duplicate atom indices")
    if any(v < 0 for v in out):
        raise ValueError(f"{name} contains negative atom indices")
    return out


@dataclass(frozen=True)
class FEPAtoms:
    """Alchemical atom partitions.

    Atom indices are always 0-based, matching OpenMM and Python.
    """

    vanish: tuple[int, ...] = field(default_factory=tuple)
    grow: tuple[int, ...] = field(default_factory=tuple)

    def __init__(
        self,
        vanish: Iterable[int] | None = None,
        grow: Iterable[int] | None = None,
    ) -> None:
        object.__setattr__(self, "vanish", _unique_ints(vanish or (), "vanish"))
        object.__setattr__(self, "grow", _unique_ints(grow or (), "grow"))
        overlap = set(self.vanish).intersection(self.grow)
        if overlap:
            raise ValueError(f"atoms cannot be both vanish and grow: {sorted(overlap)}")

    @classmethod
    def from_alch_file(cls, path: str | Path, alch_col: str = "B") -> "FEPAtoms":
        vanish, grow = parse_alch_file(path, alch_col=alch_col)
        return cls(vanish=vanish, grow=grow)

    @property
    def alchemical(self) -> tuple[int, ...]:
        return tuple(sorted((*self.vanish, *self.grow)))


@dataclass(frozen=True)
class LambdaWindow:
    index: int
    lambda1: float
    lambda2: float
    steps: int
    lambda_idws: float | None = None
    sample_backward: bool = True


@dataclass(frozen=True)
class FEPConfig:
    alch_atoms: FEPAtoms
    windows: Sequence[float]
    steps: Sequence[int]
    alch_vdw_shift_coeff: float = 5.0
    alch_elec_lambda_start: float = 0.5
    alch_vdw_lambda_end: float = 1.0
    alch_out_file: str = "run.fepout"
    alch_out_freq: int = 5
    alch_equil_steps: int = 0
    alch_decouple: bool = False
    double_wide: bool = True
    temperature: float = 300.0
    compatibility: str = "namd-compatible"

    def __post_init__(self) -> None:
        windows = tuple(float(v) for v in self.windows)
        steps = tuple(int(v) for v in self.steps)
        object.__setattr__(self, "windows", windows)
        object.__setattr__(self, "steps", steps)
        self.validate()

    def validate(self, num_particles: int | None = None) -> None:
        if len(self.windows) < 2:
            raise ValueError("windows must contain at least two lambda values")
        if len(self.steps) != len(self.windows) - 1:
            raise ValueError("steps length must be len(windows) - 1")
        if any(step <= 0 for step in self.steps):
            raise ValueError("all window steps must be positive")
        for name, value in (
            ("alch_elec_lambda_start", self.alch_elec_lambda_start),
            ("alch_vdw_lambda_end", self.alch_vdw_lambda_end),
        ):
            if value < 0.0 or value > 1.0:
                raise ValueError(f"{name} must be in [0, 1]")
        if self.alch_vdw_shift_coeff < 0.0:
            raise ValueError("alch_vdw_shift_coeff must be non-negative")
        if self.alch_out_freq <= 0:
            raise ValueError("alch_out_freq must be positive")
        if self.alch_equil_steps < 0:
            raise ValueError("alch_equil_steps cannot be negative")
        if self.temperature <= 0.0:
            raise ValueError("temperature must be positive")
        for value in self.windows:
            if value < 0.0 or value > 1.0:
                raise ValueError("all windows must be in [0, 1]")
        if num_particles is not None:
            bad = [i for i in self.alch_atoms.alchemical if i >= num_particles]
            if bad:
                raise ValueError(f"atom indices out of range for system: {bad}")

    def lambda_state(self, lam: float) -> LambdaState:
        return lambda_state(
            lam,
            elec_start=self.alch_elec_lambda_start,
            vdw_end=self.alch_vdw_lambda_end,
        )

    def iter_windows(self) -> tuple[LambdaWindow, ...]:
        out = []
        for idx, lambda1 in enumerate(self.windows[:-1]):
            lambda2 = self.windows[idx + 1]
            lambda_idws = self.windows[idx - 1] if self.double_wide and idx > 0 else None
            out.append(LambdaWindow(idx, lambda1, lambda2, self.steps[idx], lambda_idws))
        endpoint = self.windows[-1]
        if self.double_wide and (endpoint < 1.0e-12 or endpoint > 1.0 - 1.0e-12):
            out.append(
                LambdaWindow(
                    len(out),
                    endpoint,
                    self.windows[-2],
                    self.steps[-1],
                    lambda_idws=-1.0,
                    sample_backward=False,
                )
            )
        return tuple(out)
