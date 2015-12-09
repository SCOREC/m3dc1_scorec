cmake .. \
  -DCMAKE_C_COMPILER=cc \
  -DCMAKE_CXX_COMPILER=CC \
  -DCMAKE_C_FLAGS="-g -O3"\
  -DCMAKE_CXX_FLAGS="-g -O3" \
  -DENABLE_THREADS=ON \
  -DENABLE_ZOLTAN=ON \
  -DZOLTAN_INCLUDE_DIR="$CRAY_TRILINOS_PREFIX_DIR/include" \
  -DZOLTAN_LIBRARY="$CRAY_TRILINOS_PREFIX_DIR/lib/libzoltan.a" \
  -DPARMETIS_INCLUDE_DIR="$CRAY_TPSL_PREFIX_DIR/include" \
  -DPARMETIS_LIBRARY="$CRAY_TPSL_PREFIX_DIR/lib/libparmetis.a" \
  -DMETIS_LIBRARY="$CRAY_TPSL_PREFIX_DIR/lib/libmetis.a" \
  -DCMAKE_INSTALL_PREFIX="/global/project/projectdirs/mp288/hopper/scorec/Sep2015" \
  -DCMAKE_BUILD_TYPE=Debug