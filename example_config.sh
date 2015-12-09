#!/bin/bash

cmake .. \
  -DCMAKE_C_COMPILER="mpicc" \
  -DCMAKE_CXX_COMPILER="mpicxx" \
  -DCMAKE_Fortran_COMPILER="mpif90" \
  -DCMAKE_C_FLAGS=" -g -O2" \
  -DCMAKE_CXX_FLAGS=" -g -O2" \
  -DCMAKE_Fortran_FLAGS="-fpic "\
  -DENABLE_COMPLEX=OFF \
  -DSCOREC_INCLUDE_DIR="/users/seol/develop/include" \
  -DSCOREC_LIB_DIR="/users/seol/develop/lib" \
  -DZOLTAN_LIBRARY="/usr/local/zoltan/latest/openmpi1.6.5_gcc4.4.5/lib/libzoltan.a" \
  -DPARMETIS_LIBRARY="/usr/local/parmetis/latest/openmpi1.6.5_gcc4.4.5_dbg/lib/libparmetis.a" \
  -DMETIS_LIBRARY="/usr/local/parmetis/latest/openmpi1.6.5_gcc4.4.5_dbg/lib/libmetis.a" \
  -DPETSC_INCLUDE_DIR="/lore/fzhang/petsc-real/petsc-3.5.1/include" \
  -DPETSC_LIB_DIR="/lore/fzhang/petsc-real/petsc-3.5.1/test/lib" \
  -DHDF5_INCLUDE_DIR="/lore/fzhang/petsc-real/petsc-3.5.1/test/include" \
  -DHDF5_LIB_DIR="/lore/fzhang/petsc-real/petsc-3.5.1/test/lib" \
  -DENABLE_MESHGEN=OFF \
  -DSIMMETRIX_INCLUDE_DIR=/net/common/meshSim/latest/include \
  -DSIMMETRIX_LIB_DIR=/net/common/meshSim/latest/lib/x64_rhel5_gcc41 \
  -DENABLE_TESTING=ON \
  -DCMAKE_INSTALL_PREFIX="$PWD/../install"
