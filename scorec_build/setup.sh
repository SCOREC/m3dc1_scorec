module load git
module load cmake/latest
module load mpich3/3.1.2-thread-multiple

export CMAKE_PREFIX_PATH=$CMAKE_PREFIX_PATH:$PWD/core/install

export CMAKE_PREFIX_PATH=$CMAKE_PREFIX_PATH:/lore/dibanez/m3dc1/petsc-3.6.3/install_dir
export CMAKE_PREFIX_PATH=$CMAKE_PREFIX_PATH:/lore/dibanez/m3dc1/Zoltan_v3.8/install
export CMAKE_PREFIX_PATH=$CMAKE_PREFIX_PATH:/lore/dibanez/m3dc1/parmetis-4.0.3/install

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/lore/dibanez/m3dc1/petsc-3.6.3/install_dir/lib
