#This is the top M3DC1_SCOREC CMakeList File for the Build

#Setting Version Number, Project Name
cmake_minimum_required (VERSION 2.8)
project (m3dc1_scorec)
enable_language (Fortran)

# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

# default installation
#get_filename_component (default_prefix ".." ABSOLUTE)
#set (CMAKE_INSTALL_PREFIX ${default_prefix} CACHE STRING
#      "Choose the installation directory; by default it installs in the NORMA directory."
#      FORCE)

# FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
message (Fortran_COMPILER_NAME = ${Fortran_COMPILER_NAME})

#unless building shared libs, then select static libs 
# if both static and shared libs are available 
set(CMAKE_FIND_LIBRARY_SUFFIXES ".a" ".so") 
if(BUILD_SHARED_LIBS)
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".so" ".a")
endif()

#Settings options for testing
enable_testing()
include(CTest)
#This will be set to ON by the CTest driver script (and only by that)
option(ENABLE_TESTING "Build for CTest" OFF)
set(MPIRUN "mpirun"
    CACHE string 
    "the mpirun or srun executable")
set(MPIRUN_PROCFLAG "-np"
    CACHE string 
    "the command line flag to give process count to MPIRUN")

#Doxygen generation system
find_package(Doxygen)
if(DOXYGEN_FOUND)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in
               ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
add_custom_target(doc
${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
COMMENT "Generating API documentation with Doxygen" VERBATIM
)
endif(DOXYGEN_FOUND)

configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/m3dc1_scorecConfig.cmake.in"
    "${CMAKE_BINARY_DIR}/m3dc1_scorecConfig.cmake")

#Source and header files
set(SOURCES
  api/m3dc1_scorec.cc
  src/BSpline.cc
  src/Expression.cc
  src/m3dc1_ls.cc
  src/m3dc1_mesh.cc
  src/m3dc1_ghost.cc
  src/m3dc1_sizeField.cc
  src/PolyNomial.cc
  src/ReducedQuinticImplicit.cc
  src/CMODGeoExpression.cc
  src/m3dc1_matrix.cc
  src/m3dc1_model.cc
  src/m3dc1_slnTransfer.cc
  src/ReducedQuinticExplicit.cc
  src/slntransferUtil.cc 
  fortran/c1_12_element.F90
  fortran/explicit_element.F90
  fortran/explicit_element_interface.F90)

set(HEADERS
  api/name_convert.h
  api/m3dc1_scorec.h 
  include/BSpline.h
  include/m3dc1_field.h
  include/m3dc1_ls.h
  include/m3dc1_mesh.h
  include/m3dc1_ghost.h
  include/m3dc1_model.h
  include/m3dc1_slnTransfer.h
  include/PolyNomial.h
  include/ReducedQuinticImplicit.h
  include/CMODGeoExpression.h  
  include/m3dc1_matrix.h
  include/m3dc1_sizeField.h
  include/ReducedQuinticExplicit.h
  include/slntransferUtil.h
  include/Expression.h

  include/ReducedQuintic.h)

set(CMAKE_MODULE_PATH 
   ${CMAKE_MODULE_PATH} 
   "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/api)
set(TARGET_LIB_NAME m3dc1_scorec)

find_package(Scorec QUIET REQUIRED)
include_directories(${SCOREC_INCLUDE_DIRS})
set(DEP_LIBS ${DEP_LIBS} ${SCOREC_LIBRARIES} )

if (ENABLE_OMEGA_H)
  find_path(OMEGA_H_INCLUDE_DIR omega_h.h)
  if(NOT EXISTS "${OMEGA_H_INCLUDE_DIR}")
    message(FATAL_ERROR "Omega_h header not found")
  endif()
  find_library(OMEGA_H_LIB omega_h)
  if(NOT EXISTS "${OMEGA_H_LIB}")
    message(FATAL_ERROR "Omega_h library not found")
  endif()
  include_directories(${OMEGA_H_INCLUDE_DIR})
  set(DEP_LIBS ${DEP_LIBS} ${OMEGA_H_LIB})
endif()

if (ENABLE_MESHGEN)
find_package(Simmetrix QUIET REQUIRED)
include_directories(${SIMMETRIX_INCLUDE_DIRS})
find_package(Lapack QUIET REQUIRED)
set(DEP_LIBS ${DEP_LIBS} ${LAPACK_LIBRARIES} )
set(DEP_LIBS ${DEP_LIBS} ${SIMMETRIX_LIBRARIES} )
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DM3DC1_MESHGEN")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DM3DC1_MESHGEN")
set(TARGET_LIB_NAME ${TARGET_LIB_NAME}_meshgen)
else()
if (ENABLE_TRILINOS)
find_package(Trilinos QUIET REQUIRED)
find_package(Boost QUIET REQUIRED)
find_package(Lapack QUIET REQUIRED)
include_directories(${TRILINOS_INCLUDE_DIRS})
set(DEP_LIBS ${DEP_LIBS} ${LAPACK_LIBRARIES} )
set(DEP_LIBS ${DEP_LIBS} ${TRILINOS_LIBRARIES})
set(DEP_LIBS ${DEP_LIBS} ${BOOST_LIBRARIES})
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DM3DC1_TRILINOS")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DM3DC1_TRILINOS")
set(TARGET_LIB_NAME ${TARGET_LIB_NAME}_trilinos)
endif()
if(ENABLE_COMPLEX)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DPETSC_USE_COMPLEX")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DPETSC_USE_COMPLEX")
set(TARGET_LIB_NAME ${TARGET_LIB_NAME}_complex)
endif()
find_package(Petsc QUIET REQUIRED)
include_directories(${PETSC_INCLUDE_DIRS})
set(DEP_LIBS ${DEP_LIBS} ${PETSC_LIBRARIES} )
endif()
if(ENABLE_OMEGA_H)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DM3DC1_OMEGA_H")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DM3DC1_OMEGA_H")
set(TARGET_LIB_NAME ${TARGET_LIB_NAME})
endif()

add_library(${TARGET_LIB_NAME} ${SOURCES})
target_link_libraries(${TARGET_LIB_NAME} ${DEP_LIBS})

INSTALL(FILES ${HEADERS} DESTINATION include)
INSTALL(TARGETS ${TARGET_LIB_NAME}
      LIBRARY DESTINATION lib
      ARCHIVE DESTINATION lib)

if (ENABLE_TESTING)
if (ENABLE_MESHGEN)
  add_executable(create_smd test/meshgen/main.cc)
  target_link_libraries(create_smd ${DEP_LIBS})
  target_link_libraries(create_smd ${TARGET_LIB_NAME})
  INSTALL(TARGETS create_smd RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX})
else()
  if (ENABLE_TRILINOS)
    add_executable(epetra test/epetra/main.cc)
    target_link_libraries(epetra ${DEP_LIBS})
    target_link_libraries(epetra ${TARGET_LIB_NAME})
    INSTALL(TARGETS epetra RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
  endif()

  add_executable(solver test/solver/main.cc)
  target_link_libraries(solver ${DEP_LIBS})
  target_link_libraries(solver ${TARGET_LIB_NAME})
  INSTALL(TARGETS solver RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
endif()
  add_executable(convert_polar test/convert_polar/main.cc)
  target_link_libraries(convert_polar ${DEP_LIBS})
  target_link_libraries(convert_polar ${TARGET_LIB_NAME})
  INSTALL(TARGETS convert_polar RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)

  if (ENABLE_OMEGA_H)
    add_executable(omega_h_test test/omega_h/main.cc)
    target_link_libraries(omega_h_test ${DEP_LIBS})
    target_link_libraries(omega_h_test ${TARGET_LIB_NAME})
    INSTALL(TARGETS omega_h_test RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
    
    add_executable(ghosting_test test/ghosting/main.cc)
    target_link_libraries(ghosting_test ${DEP_LIBS})
    target_link_libraries(ghosting_test ${TARGET_LIB_NAME})
    INSTALL(TARGETS ghosting_test RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
    
    add_executable(create_smb test/create_smb/create_smb.cc)
    target_link_libraries(create_smb ${DEP_LIBS})
    target_link_libraries(create_smb ${TARGET_LIB_NAME})
    INSTALL(TARGETS create_smb RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
  endif()
  
endif()

#binary distribution package
set(CPACK_GENERATOR "TGZ")
set(CPACK_PACKAGE_VERSION "1.0.1")
include(CPack)

message("ENV_FLAGS = $ENV{FFLAGS}")
message("CMAKE_C_FLAGS = ${CMAKE_C_FLAGS}")
message("CMAKE_CXX_FLAGS = ${CMAKE_CXX_FLAGS}")
message("CMAKE_Fortran_COMPILER_INIT = ${CMAKE_Fortran_COMPILER_INIT}")
message("CMAKE_Fortran_COMPILER_FULLPATH =
${CMAKE_Fortran_COMPILER_FULLPATH}")
message("CMAKE_Fortran_COMPILER = ${CMAKE_Fortran_COMPILER}")
message("CMAKE_Fortran_FLAGS = ${CMAKE_Fortran_FLAGS}")
message("CMAKE_INSTALL_PREFIX = ${CMAKE_INSTALL_PREFIX}")


