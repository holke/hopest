dnl
dnl This file is part of hopest.
dnl

AC_INIT([hopest],
        [m4_esyscmd([build-aux/git-version-gen .tarball-version])],
        [burstedde@ins.uni-bonn.de])
AC_PREREQ(2.61)
AC_CONFIG_HEADERS([src/pre_config.h])
AC_CONFIG_SRCDIR([src/hopest.h])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([config])
AC_PREFIX_DEFAULT([$PWD/local])
AX_PREFIX_CONFIG_H([src/hopest_config.h])
AM_INIT_AUTOMAKE([parallel-tests subdir-objects])
AM_SILENT_RULES
SC_VERSION([HOPEST])

echo "o---------------------------------------"
echo "| Checking options"
echo "o---------------------------------------"

HOPEST_ARG_ENABLE([debug], [enable debug mode (assertions and extra checks)],
                  [DEBUG])

echo "o---------------------------------------"
echo "| Checking MPI and related programs"
echo "o---------------------------------------"

dnl A nonempty second/third argument causes to enable F77/CXX, respectively.
SC_MPI_CONFIG([HOPEST], [yes], [yes])
SC_MPI_ENGAGE([HOPEST])
# This is needed for compatibility with automake >= 1.12
m4_ifdef([AM_PROG_AR],[AM_PROG_AR])
SC_PROG_LINT
LT_INIT

echo "o---------------------------------------"
echo "| Checking libraries"
echo "o---------------------------------------"

SC_CHECK_LIBRARIES([CON4M])
P4EST_CHECK_LIBRARIES([CON4M])
CON4M_CHECK_LIBRARIES([CON4M])

echo "o---------------------------------------"
echo "| Checking subpackages"
echo "o---------------------------------------"

CON4M_DIST_DENY=
CON4M_DISTCLEAN=
SC_AS_SUBPACKAGE([CON4M])
export P4EST_SC_SOURCE=../sc
P4EST_AS_SUBPACKAGE([CON4M])
AM_CONDITIONAL([CON4M_DIST_DENY], [test "x$CON4M_DIST_DENY" != x])
AC_SUBST([CON4M_DISTCLEAN])

# Print summary.

AC_DEFINE_UNQUOTED(F77,         ["${F77}"],         [F77 compiler])
AC_DEFINE_UNQUOTED(FFLAGS,      ["${FFLAGS}"],      [F77 compiler flags])
AC_DEFINE_UNQUOTED(CPP,         ["${CPP}"],         [C preprocessor])
AC_DEFINE_UNQUOTED(CPPFLAGS,    ["${CPPFLAGS}"],    [C preprocessor flags])
AC_DEFINE_UNQUOTED(CC,          ["${CC}"],          [C compiler])
AC_DEFINE_UNQUOTED(CFLAGS,      ["${CFLAGS}"],      [C compiler flags])
AC_DEFINE_UNQUOTED(LDFLAGS,     ["${LDFLAGS}"],     [Linker flags])
AC_DEFINE_UNQUOTED(LIBS,        ["${LIBS}"],        [Libraries])

echo "o----------------------------------"
echo "| Results for hopest are"
echo "o----------------------------------"
echo "| F77:         $F77"
echo "| FFLAGS:      $FFLAGS"
echo "| CPP:         $CPP"
echo "| CPPFLAGS:    $CPPFLAGS"
echo "| CC:          $CC"
echo "| CFLAGS:      $CFLAGS"
echo "| LDFLAGS:     $LDFLAGS"
echo "| LIBS:        $LIBS"
echo "o----------------------------------"

# Create output files.
AC_CONFIG_FILES([Makefile Makefile.hopest.pre])
AC_OUTPUT

# Final messages.
SC_FINAL_MESSAGES([CON4M])
P4EST_FINAL_MESSAGES([CON4M])
CON4M_FINAL_MESSAGES([CON4M])
