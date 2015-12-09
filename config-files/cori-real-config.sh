cmake .. \
  -DCMAKE_C_COMPILER="/opt/cray/craype/2.4.2/bin/cc" \
  -DCMAKE_CXX_COMPILER="/opt/cray/craype/2.4.2/bin/CC" \
  -DCMAKE_Fortran_COMPILER="/opt/cray/craype/2.4.2/bin/ftn" \
  -DCMAKE_C_FLAGS=" -g -O3 -DDEBUG -I$PETSC_DIR/include" \
  -DCMAKE_CXX_FLAGS=" -g -O3 -DDEBUG -I$PETSC_DIR/include" \
  -DCMAKE_Fortran_FLAGS="-fpic -g -O3"\
  -DSCOREC_INCLUDE_DIR="/global/project/projectdirs/mp288/cori/scorec/Nov2015/include" \
  -DSCOREC_LIB_DIR="/global/project/projectdirs/mp288/cori/scorec/Nov2015/lib" \
  -DZOLTAN_LIBRARY="$CRAY_TRILINOS_PREFIX_DIR/lib/libzolta.a" \
  -DPARMETIS_LIBRARY="$CRAY_TPSL_PREFIX_DIR/lib/libparmetis.a" \
  -DMETIS_LIBRARY="$CRAY_TPSL_PREFIX_DIR/lib/libmetis.a" \
  -DPETSC_INCLUDE_DIR="$PETSC_DIR/$PETSC_ARCH/include" \
  -DPETSC_LIB_DIR="$PETSC_DIR/$PETSC_ARCH/lib" \
  -DHDF5_INCLUDE_DIR="$HDF5_DIR/include" \
  -DHDF5_LIB_DIR="$HDF5_DIR/lib" \
  -DCMAKE_INSTALL_PREFIX="/global/project/projectdirs/mp288/cori/scorec/Nov2015" \
  -DENABLE_TESTING=OFF \
  -DENABLE_COMPLEX=OFF \
  -DCMAKE_BUILD_TYPE=Debug

