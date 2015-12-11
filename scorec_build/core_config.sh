#!/bin/bash

cmake .. \
  -DCMAKE_C_COMPILER="mpicc" \
  -DCMAKE_CXX_COMPILER="mpicxx" \
  -DCMAKE_C_FLAGS="-g -O2 -Wall -Wextra" \
  -DCMAKE_CXX_FLAGS="-g -O2 -Wall -Wextra" \
  -DENABLE_ZOLTAN=ON \
  -DENABLE_THREADS=OFF \
  -DCMAKE_INSTALL_PREFIX="$PWD/../install"
