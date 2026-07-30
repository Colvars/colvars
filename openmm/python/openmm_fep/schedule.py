from __future__ import annotations

from dataclasses import dataclass


def _clamp01(value: float) -> float:
    return max(0.0, min(1.0, float(value)))


def vdw_lambda_grow(lam: float, vdw_end: float = 1.0) -> float:
    lam = _clamp01(lam)
    if vdw_end <= 0.0:
        return 1.0
    return 1.0 if lam >= vdw_end else lam / vdw_end


def vdw_lambda_vanish(lam: float, vdw_end: float = 1.0) -> float:
    lam = _clamp01(lam)
    start = 1.0 - vdw_end
    if vdw_end <= 0.0:
        return 0.0
    if lam <= start:
        return 1.0
    return max(0.0, (1.0 - lam) / vdw_end)


def elec_lambda_grow(lam: float, elec_start: float = 0.5) -> float:
    lam = _clamp01(lam)
    if elec_start >= 1.0:
        return 0.0 if lam < 1.0 else 1.0
    return 0.0 if lam <= elec_start else (lam - elec_start) / (1.0 - elec_start)


def elec_lambda_vanish(lam: float, elec_start: float = 0.5) -> float:
    lam = _clamp01(lam)
    end = 1.0 - elec_start
    if elec_start >= 1.0:
        return 1.0 if lam <= 0.0 else 0.0
    if lam >= end:
        return 0.0
    return max(0.0, 1.0 - lam / (1.0 - elec_start))


@dataclass(frozen=True)
class LambdaState:
    lambda_value: float
    grow_electrostatics: float
    vanish_electrostatics: float
    grow_vdw: float
    vanish_vdw: float


def lambda_state(lam: float, elec_start: float = 0.5, vdw_end: float = 1.0) -> LambdaState:
    return LambdaState(
        lambda_value=_clamp01(lam),
        grow_electrostatics=elec_lambda_grow(lam, elec_start),
        vanish_electrostatics=elec_lambda_vanish(lam, elec_start),
        grow_vdw=vdw_lambda_grow(lam, vdw_end),
        vanish_vdw=vdw_lambda_vanish(lam, vdw_end),
    )

