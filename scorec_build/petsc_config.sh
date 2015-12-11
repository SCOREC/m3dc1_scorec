#!/bin/bash

mkdir -p install_dir/include
mkdir -p install_dir/lib
mkdir -p install_dir/bin

./configure \
  --with-mpi=1 \
  --with-mpi-compilers=1 \
  --with-cc=mpicc \
  --with-cxx=mpicxx \
  --with-fc=mpif90 \
  --download-fblaslapack \
  --download-blacs \
  --download-scalapack \
  --download-superlu_dist \
  --download-superlu \
  --download-mumps \
  --download-hdf5=yes \
  --with-metis-dir=$PWD/../parmetis-4.0.3/install \
  --with-parmetis-dir=$PWD/../parmetis-4.0.3/install \
  --with-c++-support \
  --with-fortran-interfaces=1 \
  --with-debugging=1 \
  --prefix=$PWD/install_dir
