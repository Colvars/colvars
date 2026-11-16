# CMake configuration for colvartools
# Replaces the Makefile build system

# Get the source directory path
get_filename_component(COLVARS_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/.. ABSOLUTE)

# Common include directories for all tools
set(TOOLS_INCLUDE_DIRS ${COLVARS_SOURCE_DIR}/src)

# ---- abf_integrate: integrates ABF data ----
# This tool does not link with libcolvars, but needs the version header
add_executable(abf_integrate abf_integrate.cpp abf_data.cpp)
target_include_directories(abf_integrate PRIVATE ${TOOLS_INCLUDE_DIRS})

# ---- poisson_integrator: Poisson equation solver on colvars grid ----
add_executable(poisson_integrator poisson_integrator.cpp)
target_link_libraries(poisson_integrator PRIVATE colvars)
target_include_directories(poisson_integrator PRIVATE ${TOOLS_INCLUDE_DIRS})

# ---- poisson_integrator_conv: convergence testing tool ----
add_executable(poisson_integrator_conv poisson_integrator_conv.cpp)
target_link_libraries(poisson_integrator_conv PRIVATE colvars)
target_include_directories(poisson_integrator_conv PRIVATE ${TOOLS_INCLUDE_DIRS})
