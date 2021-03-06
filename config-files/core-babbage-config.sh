cmake .. \
  -DCMAKE_C_COMPILER="/global/babbage/nsg/opt/intel/impi/4.1.3.048/intel64/bin/mpiicc" \
  -DCMAKE_CXX_COMPILER="/global/babbage/nsg/opt/intel/impi/4.1.3.048/intel64/bin/mpiicpc" \
  -DCMAKE_C_FLAGS="-mmic -O3 -g -opt-assume-safe-padding -opt-streaming-stores always -opt-streaming-cache-evict=0" \
  -DCMAKE_CXX_FLAGS="-mmic -O3 -g -opt-assume-safe-padding -opt-streaming-stores always -opt-streaming-cache-evict=0 -DMPICH_IGNORE_CXX_SEEK" \
  -DENABLE_THREADS=ON \
  -DCMAKE_EXE_LINKER_FLAGS="-mt_mpi " \
  -DENABLE_ZOLTAN=ON \
  -DZOLTAN_INCLUDE_DIR="$ZOLTAN_DIR/include" \
  -DZOLTAN_LIBRARY="$ZOLTAN_DIR/lib/libzoltan.a" \
  -DPARMETIS_INCLUDE_DIR="$PETSC_DIR/$PETSC_ARCH/include" \
  -DPARMETIS_LIBRARY="$PETSC_DIR/$PETSC_ARCH/lib/libparmetis.a" \
  -DMETIS_LIBRARY="$PETSC_DIR/$PETSC_ARCH/lib/libmetis.a" \
  -DCMAKE_INSTALL_PREFIX="/global/u1/s/seol/babbage" \
  -DIS_TESTING=OFF \
  -DCMAKE_BUILD_TYPE=Debug
