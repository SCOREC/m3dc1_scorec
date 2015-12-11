#!/bin/bash -x
make config cc=mpicc prefix=$PWD/../install debug=1 assert=1
