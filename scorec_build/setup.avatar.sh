module load mpich3/latest
module load cmake/latest

export CMAKE_PREFIX_PATH=$CMAKE_PREFIX_PATH:/lore/kalyak/avatar/m3dc1/core-install
export CMAKE_PREFIX_PATH=$CMAKE_PREFIX_PATH:/lore/kalyak/avatar/m3dc1/interface-install
export CMAKE_PREFIX_PATH=$CMAKE_PREFIX_PATH:/lore/kalyak/avatar/m3dc1/omega_h-install
export CMAKE_PREFIX_PATH=$CMAKE_PREFIX_PATH:/lore/kalyak/avatar/petsc-install
export CMAKE_PREFIX_PATH=$CMAKE_PREFIX_PATH:/lore/kalyak/avatar/zoltan-install
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/lore/kalyak/avatar/petsc-install/lib

