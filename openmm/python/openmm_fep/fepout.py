from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import TextIO


KCAL_PER_KJ = 1.0 / 4.184
R_KCAL_PER_MOL_K = 0.00198720425864083


def kj_per_mol_to_kcal_per_mol(value: float) -> float:
    return value * KCAL_PER_KJ


@dataclass
class FEPStatistics:
    temperature: float
    samples: int = field(init=False, default=0)
    sum_de: float = field(init=False, default=0.0)
    sum_exp: float = field(init=False, default=0.0)
    _log_sum_exp: float = field(init=False, default=-math.inf, repr=False)

    def __post_init__(self) -> None:
        if not math.isfinite(self.temperature) or self.temperature <= 0.0:
            raise ValueError("temperature must be finite and positive")

    def add(self, de_kcal_per_mol: float) -> tuple[float, float]:
        if not math.isfinite(de_kcal_per_mol):
            raise ValueError("energy difference must be finite")
        rt = R_KCAL_PER_MOL_K * self.temperature
        log_weight = -de_kcal_per_mol / rt
        if self.samples == 0:
            self._log_sum_exp = log_weight
        else:
            maximum = max(self._log_sum_exp, log_weight)
            self._log_sum_exp = maximum + math.log(
                math.exp(self._log_sum_exp - maximum) + math.exp(log_weight - maximum)
            )
        self.samples += 1
        self.sum_de += de_kcal_per_mol
        try:
            self.sum_exp = math.exp(self._log_sum_exp)
        except OverflowError:
            self.sum_exp = math.inf
        return self.dE_avg, self.dG

    @property
    def dE_avg(self) -> float:
        return self.sum_de / self.samples if self.samples else 0.0

    @property
    def dG(self) -> float:
        if not self.samples:
            return 0.0
        rt = R_KCAL_PER_MOL_K * self.temperature
        return -rt * (self._log_sum_exp - math.log(self.samples))


class FEPOutWriter:
    """NAMD-like FEP output writer."""

    def __init__(self, path: str | Path, temperature: float) -> None:
        self.path = Path(path)
        self.temperature = float(temperature)
        self._handle: TextIO | None = None
        self._net_dg = 0.0

    def __enter__(self) -> "FEPOutWriter":
        self.open()
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        self.close()

    def open(self) -> None:
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self._handle = self.path.open("w", encoding="utf-8")
        self._write("# STEP Elec vdW dE dE_avg Temp dG\n")
        self._write("# l l+dl  l l+dl E(l+dl)-E(l)\n")

    def close(self) -> None:
        if self._handle is not None:
            self._handle.close()
            self._handle = None

    def start_window(self, lambda1: float, lambda2: float, lambda_idws: float | None = None) -> None:
        line = f"#NEW FEP WINDOW: LAMBDA SET TO {lambda1:g} LAMBDA2 {lambda2:g}"
        if lambda_idws is not None:
            line += f" LAMBDA_IDWS {lambda_idws:g}"
        self._write(line + "\n")

    def equilibration_done(self, equil_steps: int, lambda1: float) -> None:
        self._write(f"#{equil_steps} STEPS OF EQUILIBRATION AT LAMBDA {lambda1:g} COMPLETED\n")
        self._write("#STARTING COLLECTION OF ENSEMBLE AVERAGE\n")

    def write_energy(
        self,
        step: int,
        elec_lambda_kcal: float,
        elec_lambda2_kcal: float,
        vdw_lambda_kcal: float,
        vdw_lambda2_kcal: float,
        de_kcal: float,
        stats: FEPStatistics,
        *,
        backward: bool = False,
        temperature: float | None = None,
    ) -> None:
        title = "FepE_back:" if backward else "FepEnergy:"
        display_temperature = self.temperature if temperature is None else float(temperature)
        self._write(
            f"{title} {step:6d} "
            f"{elec_lambda_kcal:14.4f} {elec_lambda2_kcal:14.4f} "
            f"{vdw_lambda_kcal:14.4f} {vdw_lambda2_kcal:14.4f} "
            f"{de_kcal:14.4f} {stats.dE_avg:14.4f} "
            f"{display_temperature:14.4f} {stats.dG:14.4f}\n"
        )

    def finish_window(self, lambda1: float, lambda2: float, dg: float) -> None:
        self._net_dg += dg
        self._write(
            f"#Free energy change for lambda window [ {lambda1:g} {lambda2:g} ] "
            f"is {dg:g} ; net change until now is {self._net_dg:g}\n"
        )

    def _write(self, text: str) -> None:
        if self._handle is None:
            raise RuntimeError("FEPOutWriter is not open")
        self._handle.write(text)
        self._handle.flush()
