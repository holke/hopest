dnl
dnl This file is part of hopest.
dnl

dnl HOPEST_COMPILE_LINK_HDF5(PREFIX)
dnl
AC_DEFUN([HOPEST_COMPILE_LINK_HDF5],
[AC_LINK_IFELSE([AC_LANG_PROGRAM(dnl
[[
  #include <hdf5.h>
]], [[
  hsize_t dim=4;
  hid_t dataspace; 
  dataspace = H5Screate_simple(1, &dim, NULL);
]])],, [AC_MSG_ERROR([unable to link])])
])

dnl HOPEST_CHECK_HDF5(PREFIX)
dnl 
dnl Configure the hdf5 library.
dnl By default it is disabled and hdf5 is not activated in hopest.
dnl If --with-hdf5 is specified, we try to link hdf5 with the default
dnl include/lib options.  If --with-hdf5=<hdf5 install directory> is
dnl specified, we try to use its includes/lib subdirectories.  In both
dnl cases we abort the configuration if a test program cannot be built.
dnl
AC_DEFUN([HOPEST_CHECK_HDF5],
[
AC_MSG_CHECKING([for hdf5])
SC_ARG_WITH_PREFIX([hdf5],
                   [use the hdf5 library],
                   [HDF5], [$1])
if test "x$$1_WITH_HDF5" != xno ; then
    if test "x$$1_WITH_HDF5" != xyes ; then
      $1_HDF5_DIR="$$1_WITH_HDF5"
      SC_CHECK_INSTALL([$1_HDF5], [true], [true], [false], [false])
      save_CPPFLAGS="$CPPFLAGS"
      CPPFLAGS="$CPPFLAGS -I$$1_HDF5_INC"
      save_LDFLAGS="$LDFLAGS"
      LDFLAGS="$LDFLAGS -L$$1_HDF5_LIB"
    fi
    save_LIBS="$LIBS"
    LIBS="-lhdf5_hl -lhdf5 -lhdf5_fortran -ldl $LIBS"
    HOPEST_COMPILE_LINK_HDF5($1)
  AC_MSG_RESULT([successful])
else
  AC_MSG_RESULT([not used])
fi

dnl No AC_SUBST since we're changing variables directly
])
