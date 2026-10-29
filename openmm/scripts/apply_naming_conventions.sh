#!/bin/bash
# -*- sh-basic-offset: 2; sh-indentation: 2; -*-

# Apply the file and class renames of the OpenMM interface (naming
# conventions PR) to the current working tree.
#
# The renames commit on the openmm_proxy branch is exactly the output of this
# script, so it can be re-run on any branch forked before that commit to
# convert it mechanically. It is idempotent: running it twice is harmless.
#
# Conventions applied:
#   - Each file is named after the class that it contains.
#   - The Colvars proxy follows Colvars conventions (colvarproxy_<engine>),
#     classes that plug into OpenMM follow OpenMM conventions (CamelCase).
#   - All headers live in openmm/src; the include/openmm subtree is removed,
#     so that "openmm/..." include paths only refer to OpenMM itself.
#     Headers added to include/openmm on other branches are also moved to
#     src (keeping their names) and the include directives updated.
#
# Recipe to update a branch "my_branch" forked from openmm_proxy before the
# renames commit RENAMES ("refactor: apply naming conventions to the OpenMM
# interface"), whose parent RENAMES~1 adds this script:
#
#   git checkout my_branch
#   git merge RENAMES~1          # catch up to just before the renames
#   bash openmm/scripts/apply_naming_conventions.sh
#   git add -A openmm
#   git commit -m "Apply naming conventions of the OpenMM interface"
#   git merge -s ours RENAMES    # record the renames as merged
#   git merge openmm_proxy       # usual merge of any later changes
#
# Merging RENAMES with "-s ours" keeps the tree of my_branch unchanged: this
# is correct because RENAMES contains only the output of this script, which
# was just applied to my_branch. (A plain merge would also work, but conflicts
# wherever my_branch edits lines next to a renamed name.) Merge exactly the
# RENAMES commit this way, never a later commit of openmm_proxy.
#
# See the RENAMES commit message for the complete list of changes.

set -e

cd "$(git rev-parse --show-toplevel)"

if [ ! -d openmm ] ; then
  echo "Error: cannot find the openmm directory." >&2
  exit 1
fi

# Old path -> new path
file_renames=(
  include/openmm/ColvarsForce.h                    src/ColvarsForce.h
  include/openmm/ColvarsTABFIntegrator.h           src/ColvarsTABFIntegrator.h
  include/openmm/AlchemicalPMEForce.h              src/AlchemicalPMEForce.h
  include/openmm/internal/ColvarsForceImpl.h       src/ColvarsForceImpl.h
  include/openmm/internal/AlchemicalPMEForceImpl.h src/AlchemicalPMEForceImpl.h
  include/openmm/internal/AlchemicalPMEKernels.h   src/AlchemicalPMEKernels.h
  include/openmm/internal/AlchemicalPMEUtilities.h src/AlchemicalPMEUtilities.h
  include/openmm/internal/OpenMMColvarProxy.h      src/colvarproxy_openmm.h
  src/openmm_colvar_proxy.cpp                      src/colvarproxy_openmm.cpp
  src/colvars_force.cpp                            src/ColvarsForce.cpp
  src/colvars_force_impl.cpp                       src/ColvarsForceImpl.cpp
  src/colvars_tabf_integrator.cpp                  src/ColvarsTABFIntegrator.cpp
  src/alchemical_pme_force.cpp                     src/AlchemicalPMEForce.cpp
  src/alchemical_pme_force_impl.cpp                src/AlchemicalPMEForceImpl.cpp
  src/reference_alchemical_pme_kernel.cpp          src/ReferenceAlchemicalPMEKernels.cpp
  src/cuda_alchemical_pme_kernel.cpp               src/CudaAlchemicalPMEKernels.cpp
)

for ((i = 0; i < ${#file_renames[@]}; i += 2)) ; do
  old=openmm/${file_renames[i]}
  new=openmm/${file_renames[i+1]}
  if [ -e "${old}" ] ; then
    if [ -e "${new}" ] ; then
      echo "Error: both ${old} and ${new} exist." >&2
      exit 1
    fi
    git mv "${old}" "${new}"
  fi
done

# Headers of the interface, including any other header found in include/openmm
# (e.g. added on another branch), which is moved to src keeping its name
headers='ColvarsForce|ColvarsTABFIntegrator|AlchemicalPMEForce|ColvarsForceImpl'
headers+='|AlchemicalPMEForceImpl|AlchemicalPMEKernels|AlchemicalPMEUtilities'
if [ -d openmm/include/openmm ] ; then
  while read -r old ; do
    base=$(basename "${old}" .h)
    if [ -e "openmm/src/${base}.h" ] ; then
      echo "Error: both ${old} and openmm/src/${base}.h exist." >&2
      exit 1
    fi
    git mv "${old}" "openmm/src/${base}.h"
    headers+="|${base}"
  done < <(git ls-files 'openmm/include/openmm/*.h')
fi

if [ -d openmm/include ] ; then
  find openmm/include -depth -type d -empty -exec rmdir {} \;
fi
if [ -d openmm/include ] ; then
  echo "Warning: openmm/include is not empty, please move its files to openmm/src:" >&2
  find openmm/include -type f >&2
fi

# Text substitutions, applied in this order to all tracked files that contain
# any of the old names (except this script)
self=openmm/scripts/apply_naming_conventions.sh
old_names="OpenMMColvarProxy|OPENMM_INTERNAL_|include/openmm|openmm/(internal/)?(${headers})\.h"
old_names+='|(openmm_colvar_proxy|colvars_force|colvars_force_impl|colvars_tabf_integrator'
old_names+='|alchemical_pme_force|alchemical_pme_force_impl'
old_names+='|reference_alchemical_pme_kernel|cuda_alchemical_pme_kernel)\.cpp'

git grep -l -I -E "${old_names}" -- . ":!${self}" | while read -r f ; do
  HEADERS="${headers}" perl -pi -e '
    # Header files listed in CMake
    s{\$\{CMAKE_CURRENT_SOURCE_DIR\}/include/openmm/(?:internal/)?($ENV{HEADERS})\.h}
     {\${CMAKE_CURRENT_SOURCE_DIR}/src/$1.h}g;
    # Include directives
    s{\bopenmm/internal/OpenMMColvarProxy\.h}{colvarproxy_openmm.h}g;
    s{\bopenmm/(?:internal/)?($ENV{HEADERS})\.h}{$1.h}g;
    # Source files
    s{\bopenmm_colvar_proxy\.cpp}{colvarproxy_openmm.cpp}g;
    s{\bcolvars_force\.cpp}{ColvarsForce.cpp}g;
    s{\bcolvars_force_impl\.cpp}{ColvarsForceImpl.cpp}g;
    s{\bcolvars_tabf_integrator\.cpp}{ColvarsTABFIntegrator.cpp}g;
    s{\balchemical_pme_force\.cpp}{AlchemicalPMEForce.cpp}g;
    s{\balchemical_pme_force_impl\.cpp}{AlchemicalPMEForceImpl.cpp}g;
    s{\breference_alchemical_pme_kernel\.cpp}{ReferenceAlchemicalPMEKernels.cpp}g;
    s{\bcuda_alchemical_pme_kernel\.cpp}{CudaAlchemicalPMEKernels.cpp}g;
    # Classes
    s{\bOpenMMColvarProxy\b}{colvarproxy_openmm}g;
    # Include guards
    s{\bOPENMM_INTERNAL_OPENMM_COLVAR_PROXY_H\b}{COLVARPROXY_OPENMM_H}g;
    s{\bOPENMM_INTERNAL_(\w+_H)\b}{OPENMM_$1}g;
  ' "${f}"
done

# The include directory of the targets is now src
perl -pi -e 's{"\$\{CMAKE_CURRENT_SOURCE_DIR\}/include"}{"\${CMAKE_CURRENT_SOURCE_DIR}/src"}g' \
  openmm/CMakeLists.txt
