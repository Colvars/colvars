"""Optional local pytest report. Each numerical comparison records its worst element."""

from functools import wraps
import json
import math
from pathlib import Path
import sys

import numpy as np
import pytest


class NumericalResults:
    def __init__(self, path):
        self.path = Path(path)
        self.cases = {}

    @pytest.fixture(autouse=True)
    def observe_comparisons(self, request, monkeypatch):
        case = self.cases.setdefault(request.node.nodeid, {})
        measurements = case.setdefault("comparisons", [])
        spec = getattr(request.node, "callspec", None)
        engine = spec.params.get("engine") if spec else None
        if engine:
            case.update(platform=engine.platform, precision=engine.precision, use_gpu=engine.use_gpu)

        def observe(original, exact):
            @wraps(original)
            def check(actual, desired, *args, **kwargs):
                # Run the original assertion with the original operands/tolerances.
                try:
                    return original(actual, desired, *args, **kwargs)
                finally:
                    caller = sys._getframe(1)
                    try:
                        a, b = np.broadcast_arrays(actual, desired)
                    except (TypeError, ValueError):
                        a = b = np.array([])
                    if a.size and a.dtype.kind in "fiu" and b.dtype.kind in "fiu":
                        error = np.abs(a.astype(float) - b.astype(float))
                        index = int(np.argmax(error))
                        rtol = 0 if exact else kwargs.get("rtol", args[0] if args else 1e-7)
                        atol = 0 if exact else kwargs.get("atol", args[1] if len(args) > 1 else 0)
                        actual_value = a.flat[index].item()
                        expected_value = b.flat[index].item()
                        difference = float(error.flat[index])
                        measurements.append({
                            "source": f"{Path(caller.f_code.co_filename).name}:{caller.f_lineno}",
                            "actual": actual_value if math.isfinite(actual_value) else str(actual_value),
                            "expected": expected_value if math.isfinite(expected_value) else str(expected_value),
                            "absolute_error": difference if math.isfinite(difference) else str(difference),
                            "atol": atol, "rtol": rtol,
                        })
            return check

        for name in ("assert_allclose", "assert_array_equal"):
            original = getattr(np.testing, name)
            monkeypatch.setattr(np.testing, name, observe(original, name == "assert_array_equal"))

    def pytest_runtest_logreport(self, report):
        case = self.cases.setdefault(report.nodeid, {})
        if report.failed or report.skipped:
            case["status"] = "FAIL" if report.failed else "SKIP"
            case["reason"] = str(report.longrepr).splitlines()[-1]
        elif report.when == "call":
            case["status"] = "PASS"

    def pytest_sessionfinish(self, session, exitstatus):
        with self.path.open("w", encoding="utf-8") as handle:
            json.dump(self.cases, handle, separators=(",", ":"))
            handle.write("\n")
