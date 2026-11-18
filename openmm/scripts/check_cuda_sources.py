#!/usr/bin/env python3
"""Compile the CUDA adapter's runtime kernels without requiring a GPU.

This checks NVRTC source compilation only.  It is not a numerical or runtime
GPU test.  Supply the source tree matching the installed OpenMM version.
"""
from __future__ import annotations

import argparse
import ctypes
import ctypes.util
from pathlib import Path
import re


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--openmm-source", type=Path, required=True)
    parser.add_argument("--nvrtc-library", help="Path to libnvrtc.so (optional if found by the loader)")
    parser.add_argument("--architecture", default="compute_80")
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
    source_file = Path(__file__).resolve().parents[1] / "src/CudaColvarsKernel.cpp"
    adapter = re.search(r'R"CUDA\((.*?)\)CUDA"', source_file.read_text(), re.S)
    if adapter is None:
        raise RuntimeError("Could not locate the CUDA adapter source string")
    kernels = args.openmm_source / "platforms/cuda/src/kernels"
    common = (kernels / "common.cu").read_text()
    vector_ops = (kernels / "vectorOps.cu").read_text()
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
            print(f"PASS: {precision} adapter NVRTC compilation ({args.architecture}; no GPU execution)")
        finally:
            nvrtc.nvrtcDestroyProgram(ctypes.byref(program))


if __name__ == "__main__":
    main()
