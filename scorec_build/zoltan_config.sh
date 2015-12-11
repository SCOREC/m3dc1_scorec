#!/bin/bash -x
../configure --prefix=$PWD/../install --enable-mpi \
  --with-mpi-compilers=yes \
  --with-parmetis \
  --with-parmetis-libdir=$PWD/../../parmetis-4.0.3/install/lib \
  --with-parmetis-incdir=$PWD/../../parmetis-4.0.3/install/include
