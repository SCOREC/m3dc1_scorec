cmake .. \
  -DCMAKE_C_COMPILER=cc \
  -DCMAKE_CXX_COMPILER=CC \
  -DCMAKE_Fortran_COMPILER=ftn \
  -DCMAKE_C_FLAGS=" -g -O3 -DDEBUG" \
  -DCMAKE_CXX_FLAGS=" -g -O3 -DDEBUG" \
  -DCMAKE_Fortran_FLAGS="-fpic "\
  -DENABLE_COMPLEX=OFF \
  -DSCOREC_INCLUDE_DIR="/global/project/projectdirs/mp288/hopper/scorec/Sep2015/include" \
  -DSCOREC_LIB_DIR="/global/project/projectdirs/mp288/hopper/scorec/Sep2015/lib/gnu" \
  -DZOLTAN_LIBRARY="$CRAY_TRILINOS_PREFIX_DIR/lib/libzoltan.a" \
  -DPARMETIS_LIBRARY="$CRAY_TPSL_PREFIX_DIR/lib/libparmetis.a" \
  -DMETIS_LIBRARY="$CRAY_TPSL_PREFIX_DIR/lib/libmetis.a" \
  -DPETSC_INCLUDE_DIR="$PETSC_DIR/include" \
  -DPETSC_LIB_DIR="$PETSC_DIR/lib" \
  -DHDF5_INCLUDE_DIR="$HDF5_DIR/include" \
  -DHDF5_LIB_DIR="$HDF5_DIR/lib" \
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
  -DCMAKE_INSTALL_PREFIX="/global/project/projectdirs/mp288/hopper/scorec/Sep2015/"
