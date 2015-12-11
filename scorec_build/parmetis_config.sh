#!/bin/bash -x
make config cc=mpicc cxx=mpicxx prefix=$PWD/install debug=1 assert=1
