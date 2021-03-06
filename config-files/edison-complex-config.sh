cmake .. \
  -DCMAKE_C_COMPILER="/opt/cray/craype/2.2.1/bin/cc" \
  -DCMAKE_CXX_COMPILER="/opt/cray/craype/2.2.1/bin/CC" \
  -DCMAKE_Fortran_COMPILER="/opt/cray/craype/2.2.1/bin/ftn" \
  -DCMAKE_C_FLAGS=" -g -O3 -DDEBUG -I$PETSC_DIR/include" \
  -DCMAKE_CXX_FLAGS=" -g -O3 -DDEBUG -I$PETSC_DIR/include" \
  -DCMAKE_Fortran_FLAGS="-fpic -g -O3"\
  -DSCOREC_INCLUDE_DIR="/global/project/projectdirs/mp288/edison/scorec/Oct2015/include" \
  -DSCOREC_LIB_DIR="/global/project/projectdirs/mp288/edison/scorec/Oct2015/lib" \
  -DZOLTAN_LIBRARY="$PETSC_DIR/$PETSC_ARCH/lib/libzoltan.a" \
  -DPARMETIS_LIBRARY="$PETSC_DIR/$PETSC_ARCH/lib/libparmetis.a" \
  -DMETIS_LIBRARY="$PETSC_DIR/$PETSC_ARCH/lib/libmetis.a" \
  -DPETSC_INCLUDE_DIR="$PETSC_DIR/$PETSC_ARCH/include" \
  -DPETSC_LIB_DIR="$PETSC_DIR/$PETSC_ARCH/lib" \
  -DHDF5_INCLUDE_DIR="$HDF5_DIR/include" \
  -DHDF5_LIB_DIR="$HDF5_DIR/lib" \
  -DCMAKE_INSTALL_PREFIX="/global/project/projectdirs/mp288/edison/scorec/Oct2015" \
  -DENABLE_TESTING=OFF \
  -DENABLE_COMPLEX=ON \
  -DCMAKE_BUILD_TYPE=Debug

