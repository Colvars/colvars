# NAMD CudaGlobalMaster Interface with Colvars

This is an experimental Colvars plugin that can be loaded using the NAMD's new CudaGlobalMaster interface. In the GPU-resident mode, although this plugin requires to copy the data from GPU to CPU due to Colvars' CPU-only code, it still runs Colvars significantly faster than the NAMD's bundled Colvars using the traditional GlobalMaster interface, mainly because it avoids the multiple memory copies in the traditional GlobalMaster code path and only copies the atoms that requested by Colvars.

## Compilation

This plugin depends on the source code of [NAMD](https://gitlab.com/tcbgUIUC/namd) (the [`opt_cudagm` branch](https://gitlab.com/tcbgUIUC/namd/-/tree/opt_cudagm?ref_type=heads)), [the CUDA compiler](https://developer.nvidia.com/cuda-downloads) (or [the AMD ROCm/HIP compiler](https://rocm.docs.amd.com/projects/install-on-linux/en/latest/install/quick-start.html)) and libraries and optionally the [TCL library](https://www.tcl-lang.org) (if your NAMD is built with TCL support). To compile this plugin you also need [cmake version 3.23 or above](https://cmake.org/download/).

### CUDA (NVIDIA GPUs)

Commands to compile the plugin under this directory:
```sh
mkdir build
cd build
cmake -DNAMD_DIR=<your_namd_git_directory>/ -DCMAKE_BUILD_TYPE=Release -DUSE_CUDA=on -DUSE_HIP=off ../
make -j4
```

#### Notes for developers

If you built your NAMD with `-DNAMD_NVTX_ENABLED`, then you also need to switch that on in the CMake command by `-DNAMD_NVTX_ENABLED=on`. If you want to build the plugin on one computer but later want to run NAMD with it on another system, it is better to add `-DCMAKE_CUDA_ARCHITECTURES=all-major` to the CMake command.

### HIP (AMD GPUs)

If you are using AMD GPUs with the ROCm/HIP compiler, you need to the clone the [`opt_cudagm_hip` branch](https://gitlab.com/tcbgUIUC/namd/-/tree/opt_cudagm_hip), and the following commands to compile the plugin:
```sh
mkdir build
cd build
cmake -DNAMD_DIR=<your_namd_git_directory>/ -DCMAKE_BUILD_TYPE=Release -DUSE_CUDA=off -DUSE_HIP=on -DCMAKE_PREFIX_PATH=<your_rocm_installation_directory> -DCMAKE_CXX_COMPILER:PATH=<your_rocm_installation_directory>/bin/amdclang++ ../
make -j4
```

After the compilation, you will get a shared library named `libcudaglobalmastercolvars.so`.

## Example usage

This plugin is incompatible with the traditional way of using Colvars in NAMD, so you cannot use `colvars on` in your NAMD configuration and `numSteps` with it. Instead, please use the following TCL command in your NAMD configuration file before `run`:
```tcl
gpuGlobal on
gpuGlobalCreateClient <the_absolute_path_to_this_dir>/build/libcudaglobalmastercolvars.so COLVARS <colvars_config_file>
```

The example NAMD input file can be found in `<the_absolute_path_to_this_dir>/example/alad.namd`, which dynamically loads the shared library built above and runs an OPES simulation along the two dihedral angles of the alanine dipeptide.

## Limitations

This plugin is still in its early stage, and has the following limitations:
- Limited scripting support. You can call Colvars between `run`s by `gpuGlobalUpdateClient COLVARS cv <command>`. For example, you can reset Colvars by `gpuGlobalUpdateClient COLVARS cv reset`, and then load another configuration file by `gpuGlobalUpdateClient COLVARS cv configfile <new_config_file>`. However, you cannot use scripted Colvars force, and since this plugin uses an independent interpreter for the Colvars instance, so you cannot call NAMD's TCL procedures in the calculation of scripted CVs;
- `volmap` is not available;
- GaMD energy histograming/reweighting is not available (GaMD is not currently available in GPU-resident NAMD);
- Minimization is not supported since `minimize` in the NAMD GPU-resident mode actually calls the GPU-offload code path;
- SMP is not available.
