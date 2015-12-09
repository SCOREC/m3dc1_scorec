AC_INIT([m3dc1_scorec], [1.0], [fzhang@scorec.rpi.edu])
AC_CONFIG_AUX_DIR(config)
AC_CONFIG_HEADERS([config/config.h])
AM_INIT_AUTOMAKE( [subdir-objects -Wall -Werror foreign ] )
#AX_MPI
#AS_IF([test "x$MPI_MISSING" = xtrue], [echo "MPI is required, specify --with-mpi=<dir>"; exit 1])
# store current user given compiler flags to avoid default setup via AC_PROG_CXX
OLD_CXXFLAGS=$CXXFLAGS
OLD_CFLAGS=$CFLAGS
OLD_FFLAGS=$FFLAGS
OLD_FCFLAGS=$FCFLAGS
AC_PROG_CC( [$CC_LIST] )
AC_PROG_CPP
AC_PROG_CXX( [$CXX_LIST] )
AC_PROG_CXXCPP
AC_PROG_FC( [$FC_LIST] )
AC_PROG_F77( [$F77_LIST] )
AM_PROG_CC_C_O
# reset compiler flags to initial flags
CXXFLAGS=$OLD_CXXFLAGS
CFLAGS=$OLD_CFLAGS
FFLAGS=$OLD_FFLAGS
FCFLAGS=$OLD_FCFLAGS

AX_ADDLDFLAGS
AS_IF([test "x$USE_ADDLDFLAGS" = xtrue], [echo "addldflags is specified: $ADDLDFLAGS_VALUE"; ])

AC_ARG_ENABLE([petsc_complex],
  [AC_HELP_STRING([--enable-petsc_complex],[use petsc with complex number @<:@default=no@:>@. ])],
  [case "${enable_petsc_complex}" in
       yes) petsc_complex=true ;;
       no) petsc_complex=false ;;
       *) AC_MSG_ERROR([--enable-petsc_complex must be 'yes' or 'no'; ${enable_petsc_complex} was specified ]) ;;
   esac], [petsc_complex=false])
AM_CONDITIONAL([USEPETSC_COMPLEX], [test x$petsc_complex = xtrue])
AS_IF([test x$petsc_complex = xtrue], [echo "use petsc complex is specified"; ])

AX_DEBUG
AX_OPTIMIZATION
# Checks for inline support
AC_C_INLINE
# Checks for programs.
AC_PROG_LIBTOOL
AC_CONFIG_FILES([Makefile])
AC_OUTPUT
