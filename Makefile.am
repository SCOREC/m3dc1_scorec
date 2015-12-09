#bring in the macros defined in the m4 directory.. these include AX_MPI, and AX_DEBUG
ACLOCAL_AMFLAGS = -I m4

lib_LIBRARIES = libm3dc1_scorec.a

libm3dc1_scorec_a_SOURCES = \
  api/m3dc1_scorec.cc \
  src/m3dc1_matrix.cc \
  src/m3dc1_model.cc \
  src/m3dc1_mesh.cc \
  src/Expression.cc \
  src/PolyNomial.cc \
  src/BSpline.cc \
  src/CMODGeoExpression.cc \
  src/m3dc1_sizeField.cc \
  src/m3dc1_slnTransfer.cc \
  src/ReducedQuinticExplicit.cc \
  src/ReducedQuinticImplicit.cc \
  src/slntransferUtil.cc \
  fortran/c1_12_element.F90  \
  fortran/explicit_element.F90  \
  fortran/explicit_element_interface.F90

libm3dc1_scorec_a_include_HEADERS = \
  api/name_convert.h \
  api/m3dc1_scorec.h \
  include/m3dc1_field.h \
  include/m3dc1_matrix.h \
  include/m3dc1_mesh.h \
  include/m3dc1_model.h \
  include/Expression.h \
  include/PolyNomial.h \
  include/BSpline.h \
  include/CMODGeoExpression.h \
  include/m3dc1_sizeField.h \
  include/m3dc1_slnTransfer.h \
  include/ReducedQuinticExplicit.h  \
  include/ReducedQuintic.h  \
  include/ReducedQuinticImplicit.h \
  include/slntransferUtil.h

# The directory in which to install headers
libm3dc1_scorec_a_includedir = $(includedir)

# Dummy C++ source to cause C++ linking.
nodist_EXTRA_libm3dc1_scorec_a_SOURCES = dummy.cxx

AM_CXXFLAGS = \
 -I$(top_srcdir)/include \
 -I$(top_srcdir)/api \
 -I$(PREFIX)/include \
 $(PETSC_INCLUDES)

AM_CXXFLAGS += \
  $(DEBUG_FLAGS) \
  $(OPTIMIZATION_FLAGS) \
  $(PROFILING_FLAGS) 

#  PARALLEL_FLAGS = -DFMDB_PARALLEL -DPUMI_PARALLEL -DMA_PARALLEL    
#  PARALLEL_INCLUDES =  @ZOLTAN_INCLUDES@ @PARMETIS_INCLUDES@  @PETSC_INCLUDES@ 
PARALLEL_LIBS = 

if USEPETSC_COMPLEX
   APPFLAG :=$(APPFLAG) -DUSE_PETSC_COMPLEX
endif

  AM_CXXFLAGS += $(APPFLAG) $(PARALLEL_FLAGS) $(PARALLEL_INCLUDES) 
  libm3dc1_scorec_a_LIBADD = $(PARALLEL_LIBS) 


AM_CFLAGS = $(AM_CXXFLAGS)

# Put Doc directory into distribution 
EXTRA_DIST = LICENSE README fortran/explicit_element.mod

dist-hook: 
	rm -rf `find $(distdir) -name .svn`
