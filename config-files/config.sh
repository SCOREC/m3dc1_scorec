/usr/local/cmake/2.8.0/bin/cmake .. \
  -DCMAKE_C_COMPILER="/usr/local/openmpi/latest/bin/mpicc" \
  -DCMAKE_CXX_COMPILER="/usr/local/openmpi/latest/bin/mpicxx" \
  -DCMAKE_Fortran_COMPILER="/usr/local/openmpi/latest/bin/mpif90" \
  -DCMAKE_C_FLAGS=" -g -O2 -DDEBUG -I/lore/fzhang/petsc-real/petsc-3.5.1/include" \
  -DCMAKE_CXX_FLAGS=" -g -O2 -DDEBUG -I/lore/fzhang/petsc-real/petsc-3.5.1/include" \
  -DCMAKE_Fortran_FLAGS="-fpic "\
  -DENABLE_COMPLEX=OFF \
  -DSCOREC_INCLUDE_DIR="/users/seol/install/include" \
  -DSCOREC_LIB_DIR="/users/seol/install/lib" \
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
  -DENABLE_TRILINOS=OFF \
  -DTRILINOS_INCLUDE_DIR="/usr/local/trilinos/latest/include" \
  -DTRILINOS_LIB_DIR="/usr/local/trilinos/latest/lib" \
  -DLAPACK_LIB_DIR="/users/seol/public/lib" \
  -DBOOST_LIB_DIR="/users/granzb/lib" \
  -DSTDCPP_LIBRARY="/users/granzb/lib64/libstdc++.a" \
  -DNETCDF_LIBRARY="/users/granzb/lib/libnetcdf.a" \
  -DENABLE_TESTING=ON \
  -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_INSTALL_PREFIX="/users/seol/install"

