dnl
dnl This file is part of con4m.
dnl

dnl CON4M_COMPILE_LINK_HYPRE(PREFIX)
dnl
AC_DEFUN([CON4M_COMPILE_LINK_HYPRE],
[AC_LINK_IFELSE([AC_LANG_PROGRAM(dnl
[[
  #include <mpi.h>
  #include <HYPRE_IJ_mv.h>
]], [[
  HYPRE_IJMatrix A;
  HYPRE_IJMatrixCreate (MPI_COMM_WORLD, 0, 1, 0, 1, &A);
]])],, [AC_MSG_ERROR([unable to link])])
])

dnl CON4M_CHECK_HYPRE(PREFIX)
dnl
dnl Configure the hypre library for parallel linear algebra.
dnl By default it is disabled and hypre is not activated in con4m.
dnl If --with-hypre is specified, we try to link hypre with the default
dnl include/lib options.  If --with-hypre=<hypre install directory> is
dnl specified, we try to use its includes/lib subdirectories.  In both
dnl cases we abort the configuration if a test program cannot be built.
dnl
AC_DEFUN([CON4M_CHECK_HYPRE],
[
AC_MSG_CHECKING([for hypre])
SC_ARG_WITH_PREFIX([hypre],
                   [use the hypre library for parallel linear algebra],
                   [HYPRE], [$1])
if test "x$$1_WITH_HYPRE" != xno ; then
  if test "x$HAVE_PKG_MPI" = xyes ; then
    if test "x$$1_WITH_HYPRE" != xyes ; then
      $1_HYPRE_DIR="$$1_WITH_HYPRE"
      SC_CHECK_INSTALL([$1_HYPRE], [true], [true], [false], [false])
      save_CPPFLAGS="$CPPFLAGS"
      CPPFLAGS="$CPPFLAGS -I$$1_HYPRE_INC"
      save_LDFLAGS="$LDFLAGS"
      LDFLAGS="$LDFLAGS -L$$1_HYPRE_LIB"
    fi
    save_LIBS="$LIBS"
    LIBS="-lHYPRE $LIBS"
    CON4M_COMPILE_LINK_HYPRE($1)
  else
    AC_MSG_ERROR([hypre requires MPI])
  fi
  AC_MSG_RESULT([successful])
else
  AC_MSG_RESULT([not used])
fi

dnl No AC_SUBST since we're changing variables directly
])
