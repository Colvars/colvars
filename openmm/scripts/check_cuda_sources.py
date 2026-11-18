#!/usr/bin/env python3
"""Compile the CUDA adapter kernels and optionally run their buffer-level oracles.

Supply the source tree matching the installed OpenMM version. Device execution
requires NumPy and Numba. This exercises the two bridge kernels in isolation.
"""
from __future__ import annotations

import argparse
import ctypes
import ctypes.util
import json
from pathlib import Path
import re


def check_device_buffers(ptx, precision, device):
    import numpy as np
    from numba import cuda
    from numba.cuda.cudadrv import driver

    with cuda.gpus[device]:
        module = cuda.current_context().create_module_ptx(ptx)

        def launch(name, arrays, scalars=()):
            pointers = [ctypes.c_void_p(int(driver.device_pointer(a))) for a in arrays]
            func = module.get_function(name)
            driver.launch_kernel(func.handle, 1, 1, 1, 32, 1, 1, 0, 0,
                                 pointers + list(scalars))

        # Four slots from a reordered, padded 64-particle buffer. All inputs
        # are binary fractions so both the oracle and kernel arithmetic are exact.
        indices = np.array([61, 3, 47, 12], dtype=np.int32)
        offsets = np.array([[1, -2, 1, 0], [0, 0, 0, 0], [-1, 1, 0, 0],
                            [2, 0, -1, 0]], dtype=np.int32)
        box = np.array([[2., 0., 0., 0.], [.25, 3., 0., 0.], [-.5, .25, 4., 0.]])
        dtype = np.float64 if precision == "double" else np.float32
        pos = np.arange(64*4, dtype=dtype).reshape((64, 4)) / 16
        correction = np.full_like(pos, 2.0**-30)
        selected = cuda.to_device(indices)
        out = cuda.device_array((3, 4), dtype=np.float64)
        launch("gatherColvarsPositions",
               [cuda.to_device(pos), cuda.to_device(correction), selected, cuda.to_device(offsets), out],
               [(ctypes.c_double * 4)(*v) for v in box])
        actual = out.copy_to_host().T
        expected = pos[indices, :3].astype(np.float64) - offsets[:, :3] @ box[:, :3]
        if precision == "mixed":
            expected += correction[indices, :3]
        np.testing.assert_array_equal(actual, expected)

        forces = (np.arange(12).reshape(3, 4) - 6.5) / 8
        physical = np.arange(3*64, dtype=np.int64).reshape(3, 64) - 96
        accumulator = cuda.to_device(physical)
        launch("addColvarsForces", [cuda.to_device(forces), selected, accumulator])
        actual_force = accumulator.copy_to_host()
        expected_force = physical.copy()
        expected_force[:, indices] += (forces * 2**32).astype(np.int64)
        np.testing.assert_array_equal(actual_force, expected_force)
        return {
            "platform": "CUDA Driver", "device_index": device, "precision": precision,
            "scope": "isolated adapter kernels", "status": "PASS",
            "gather": {"actual_nm": actual.tolist(), "expected_nm": expected.tolist(),
                       "max_absolute_error_nm": float(np.max(np.abs(actual - expected))),
                       "tolerance_nm": 0},
            "scatter": {"actual_selected": actual_force[:, indices].tolist(),
                        "expected_selected": expected_force[:, indices].tolist(),
                        "max_absolute_error_all_192_fixed_point_entries":
                            int(np.max(np.abs(actual_force - expected_force))),
                        "tolerance_fixed_point_units": 0},
        }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--openmm-source", type=Path, required=True)
    parser.add_argument("--nvrtc-library", help="Path to libnvrtc.so (optional if found by the loader)")
    parser.add_argument("--architecture", default="compute_80")
    parser.add_argument("--device", type=int, help="Execute gather and scatter on this CUDA device")
    parser.add_argument("--json", type=Path, help="Write per-precision results")
    args = parser.parse_args()
    libname = args.nvrtc_library or ctypes.util.find_library("nvrtc")
    if not libname:
        parser.error("NVRTC not found: pass --nvrtc-library /path/to/libnvrtc.so")
    nvrtc = ctypes.CDLL(libname)
    pointer = ctypes.c_void_p
    nvrtc.nvrtcCreateProgram.argtypes = [ctypes.POINTER(pointer), ctypes.c_char_p,
        ctypes.c_char_p, ctypes.c_int, ctypes.POINTER(ctypes.c_char_p), ctypes.POINTER(ctypes.c_char_p)]
    nvrtc.nvrtcCompileProgram.argtypes = [pointer, ctypes.c_int, ctypes.POINTER(ctypes.c_char_p)]
    nvrtc.nvrtcGetProgramLogSize.argtypes = [pointer, ctypes.POINTER(ctypes.c_size_t)]
    nvrtc.nvrtcGetProgramLog.argtypes = [pointer, pointer]
    nvrtc.nvrtcDestroyProgram.argtypes = [ctypes.POINTER(pointer)]
    nvrtc.nvrtcGetPTXSize.argtypes = [pointer, ctypes.POINTER(ctypes.c_size_t)]
    nvrtc.nvrtcGetPTX.argtypes = [pointer, ctypes.c_char_p]
    source_file = Path(__file__).resolve().parents[1] / "src/CudaColvarsKernel.cpp"
    adapter = re.search(r'R"CUDA\((.*?)\)CUDA"', source_file.read_text(), re.S)
    if adapter is None:
        raise RuntimeError("Could not locate the CUDA adapter source string")
    kernels = args.openmm_source / "platforms/cuda/src/kernels"
    common = (kernels / "common.cu").read_text()
    vector_ops = (kernels / "vectorOps.cu").read_text()
    results = []
    for precision in ("single", "mixed", "double"):
        real = "double" if precision == "double" else "float"
        mixed = "float" if precision == "single" else "double"
        typedefs = "\n".join(
            f"typedef {scalar}{width} {kind}{width};"
            for kind, scalar in (("real", real), ("mixed", mixed))
            for width in ("", "2", "3", "4")
        )
        defines = "#define NUM_SELECTED 4\n#define PADDED_NUM_ATOMS 64\n"
        if precision != "single":
            defines += f"#define USE_{precision.upper()}_PRECISION 1\n"
        # OpenMM's CudaContext::createModule + compileProgram use this order.
        source = typedefs + "\n" + common + "\n" + defines + vector_ops + "\n" + adapter.group(1)
        program = pointer()
        result = nvrtc.nvrtcCreateProgram(ctypes.byref(program), source.encode(),
                                         b"colvars_adapter.cu", 0, None, None)
        if result:
            raise RuntimeError(f"nvrtcCreateProgram failed: {result}")
        try:
            options = (ctypes.c_char_p * 2)(b"--std=c++17", f"--gpu-architecture={args.architecture}".encode())
            result = nvrtc.nvrtcCompileProgram(program, len(options), options)
            size = ctypes.c_size_t()
            nvrtc.nvrtcGetProgramLogSize(program, ctypes.byref(size))
            log = ctypes.create_string_buffer(size.value)
            nvrtc.nvrtcGetProgramLog(program, log)
            if result:
                raise RuntimeError(f"{precision} NVRTC compilation failed ({result}):\n{log.value.decode()}")
            print(f"PASS: {precision} adapter NVRTC compilation ({args.architecture})")
            if args.device is None:
                results.append({"precision": precision, "scope": "NVRTC compilation", "status": "PASS"})
            else:
                size = ctypes.c_size_t()
                if nvrtc.nvrtcGetPTXSize(program, ctypes.byref(size)):
                    raise RuntimeError("NVRTC could not determine the PTX size")
                ptx = ctypes.create_string_buffer(size.value)
                if nvrtc.nvrtcGetPTX(program, ptx):
                    raise RuntimeError("NVRTC could not return compiled PTX")
                results.append(check_device_buffers(ptx.value, precision, args.device))
                print(f"PASS: {precision} device {args.device} gather and scatter buffer oracles")
        finally:
            nvrtc.nvrtcDestroyProgram(ctypes.byref(program))
    if args.json:
        args.json.write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
