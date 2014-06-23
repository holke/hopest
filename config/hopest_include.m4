dnl
dnl hopest_include.m4 - custom macros
dnl

dnl Documentation for macro names: brackets indicate optional arguments

dnl HOPEST_ARG_ENABLE(NAME, COMMENT, TOKEN)
dnl Check for --enable/disable-NAME using shell variable HOPEST_ENABLE_TOKEN
dnl If shell variable is set beforehand it overrides the option
dnl If enabled, define TOKEN to 1 and set conditional HOPEST_TOKEN
dnl Default is disabled
dnl
AC_DEFUN([HOPEST_ARG_ENABLE],
         [SC_ARG_ENABLE_PREFIX([$1], [$2], [$3], [HOPEST])])

dnl HOPEST_ARG_DISABLE(NAME, COMMENT, TOKEN)
dnl Check for --enable/disable-NAME using shell variable HOPEST_ENABLE_TOKEN
dnl If shell variable is set beforehand it overrides the option
dnl If enabled, define TOKEN to 1 and set conditional HOPEST_TOKEN
dnl Default is enabled
dnl
AC_DEFUN([HOPEST_ARG_DISABLE],
         [SC_ARG_DISABLE_PREFIX([$1], [$2], [$3], [HOPEST])])

dnl HOPEST_ARG_WITH(NAME, COMMENT, TOKEN)
dnl Check for --with/without-NAME using shell variable HOPEST_WITH_TOKEN
dnl If shell variable is set beforehand it overrides the option
dnl If with, define TOKEN to 1 and set conditional HOPEST_TOKEN
dnl Default is without
dnl
AC_DEFUN([HOPEST_ARG_WITH],
         [SC_ARG_WITH_PREFIX([$1], [$2], [$3], [HOPEST])])

dnl HOPEST_ARG_WITHOUT(NAME, COMMENT, TOKEN)
dnl Check for --with/without-NAME using shell variable HOPEST_WITH_TOKEN
dnl If shell variable is set beforehand it overrides the option
dnl If with, define TOKEN to 1 and set conditional HOPEST_TOKEN
dnl Default is with
dnl
AC_DEFUN([HOPEST_ARG_WITHOUT],
         [SC_ARG_WITHOUT_PREFIX([$1], [$2], [$3], [HOPEST])])

dnl HOPEST_AS_SUBPACKAGE(PREFIX)
dnl Call from a package that is using hopest as a subpackage.
dnl Sets PREFIX_DIST_DENY=yes if hopest is make install'd.
dnl
AC_DEFUN([HOPEST_AS_SUBPACKAGE],
         [SC_ME_AS_SUBPACKAGE([$1], [m4_tolower([$1])], [HOPEST], [hopest])])

dnl HOPEST_CHECK_LIBRARIES(PREFIX)
dnl This macro bundles the checks for all libraries and link tests
dnl that are required by hopest.  It can be used by other packages that
dnl link to hopest to add appropriate options to LIBS.
dnl
AC_DEFUN([HOPEST_CHECK_LIBRARIES],
[
HOPEST_CHECK_HDF5([$1])
])

dnl HOPEST_FINAL_MESSAGES(PREFIX)
dnl This macro prints messages at the end of the configure run.
dnl
AC_DEFUN([HOPEST_FINAL_MESSAGES],)
