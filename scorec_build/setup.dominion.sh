module load mpich3/3.2-thread-multiple-jessie

export CMAKE_PREFIX_PATH=$CMAKE_PREFIX_PATH:/lore/kalyak/dominion/m3dc1/core-install
export CMAKE_PREFIX_PATH=$CMAKE_PREFIX_PATH:/lore/kalyak/dominion/petsc-install
export CMAKE_PREFIX_PATH=$CMAKE_PREFIX_PATH:/lore/kalyak/dominion/zoltan-install
export CMAKE_PREFIX_PATH=$CMAKE_PREFIX_PATH:/lore/kalyak/dominion/petsc-install
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/lore/kalyak/dominion/petsc-install/lib
