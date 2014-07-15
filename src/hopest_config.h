#ifndef _SRC_HOPEST_CONFIG_H
#define _SRC_HOPEST_CONFIG_H 1
 
/* src/hopest_config.h. Generated automatically at end of configure. */
/* src/pre_config.h.  Generated from pre_config.h.in by configure.  */
/* src/pre_config.h.in.  Generated from configure.ac by autoheader.  */

/* C compiler */
#ifndef HOPEST_CC
#define HOPEST_CC "gcc"
#endif

/* C compiler flags */
#ifndef HOPEST_CFLAGS
#define HOPEST_CFLAGS "-g -O2"
#endif

/* C preprocessor */
#ifndef HOPEST_CPP
#define HOPEST_CPP "gcc -E"
#endif

/* C preprocessor flags */
#ifndef HOPEST_CPPFLAGS
#define HOPEST_CPPFLAGS " -I/scratch/iagbole/Codeentwicklung/hopest/hdf5/hdf5/include"
#endif

/* Define to 1 if your C++ compiler doesn't accept -c and -o together. */
/* #undef CXX_NO_MINUS_C_MINUS_O */

/* DEPRECATED (use HOPEST_ENABLE_DEBUG instead) */
/* #undef DEBUG */

/* enable debug mode (assertions and extra checks) */
/* #undef ENABLE_DEBUG */

/* Define to 1 if we are using MPI */
/* #undef ENABLE_MPI */

/* Define to 1 if we are using MPI I/O */
/* #undef ENABLE_MPIIO */

/* Define to 1 if we are using MPI_Init_thread */
/* #undef ENABLE_MPITHREAD */

/* enable POSIX threads (optionally use --enable-pthread=<PTHREAD_CFLAGS>) */
/* #undef ENABLE_PTHREAD */

/* F77 compiler */
#ifndef HOPEST_F77
#define HOPEST_F77 "gfortran"
#endif

/* Define to dummy `main' function (if any) required to link to the Fortran
   libraries. */
/* #undef F77_DUMMY_MAIN */

/* Define to a macro mangling the given C identifier (in lower and upper
   case), which must not contain underscores, for linking with Fortran. */
#ifndef HOPEST_F77_FUNC
#define HOPEST_F77_FUNC(name,NAME) name ## _
#endif

/* As F77_FUNC, but for C identifiers containing underscores. */
#ifndef HOPEST_F77_FUNC_
#define HOPEST_F77_FUNC_(name,NAME) name ## _
#endif

/* Define to 1 if your Fortran compiler doesn't accept -c and -o together. */
/* #undef F77_NO_MINUS_C_MINUS_O */

/* FC compiler */
#ifndef HOPEST_FC
#define HOPEST_FC "gfortran"
#endif

/* FC compiler flags */
#ifndef HOPEST_FCFLAGS
#define HOPEST_FCFLAGS "-I/scratch/iagbole/Codeentwicklung/hopest/hdf5/hdf5/include"
#endif

/* Define to dummy `main' function (if any) required to link to the Fortran
   libraries. */
/* #undef FC_DUMMY_MAIN */

/* Define if F77 and FC dummy `main' functions are identical. */
/* #undef FC_DUMMY_MAIN_EQ_F77 */

/* Define to a macro mangling the given C identifier (in lower and upper
   case), which must not contain underscores, for linking with Fortran. */
#ifndef HOPEST_FC_FUNC
#define HOPEST_FC_FUNC(name,NAME) name ## _
#endif

/* As FC_FUNC, but for C identifiers containing underscores. */
#ifndef HOPEST_FC_FUNC_
#define HOPEST_FC_FUNC_(name,NAME) name ## _
#endif

/* Define to 1 if your Fortran compiler doesn't accept -c and -o together. */
/* #undef FC_NO_MINUS_C_MINUS_O */

/* F77 compiler flags */
#ifndef HOPEST_FFLAGS
#define HOPEST_FFLAGS "-g -O2"
#endif

/* Define to 1 if you have the <dlfcn.h> header file. */
#ifndef HOPEST_HAVE_DLFCN_H
#define HOPEST_HAVE_DLFCN_H 1
#endif

/* Define to 1 if you have the <inttypes.h> header file. */
#ifndef HOPEST_HAVE_INTTYPES_H
#define HOPEST_HAVE_INTTYPES_H 1
#endif

/* Have we found function pthread_create. */
#ifndef HOPEST_HAVE_LPTHREAD
#define HOPEST_HAVE_LPTHREAD 1
#endif

/* Have we found function lua_createtable. */
/* #undef HAVE_LUA */

/* Define to 1 if you have the <memory.h> header file. */
#ifndef HOPEST_HAVE_MEMORY_H
#define HOPEST_HAVE_MEMORY_H 1
#endif

/* Define to 1 if you have the <stdint.h> header file. */
#ifndef HOPEST_HAVE_STDINT_H
#define HOPEST_HAVE_STDINT_H 1
#endif

/* Define to 1 if you have the <stdlib.h> header file. */
#ifndef HOPEST_HAVE_STDLIB_H
#define HOPEST_HAVE_STDLIB_H 1
#endif

/* Define to 1 if you have the <strings.h> header file. */
#ifndef HOPEST_HAVE_STRINGS_H
#define HOPEST_HAVE_STRINGS_H 1
#endif

/* Define to 1 if you have the <string.h> header file. */
#ifndef HOPEST_HAVE_STRING_H
#define HOPEST_HAVE_STRING_H 1
#endif

/* Define to 1 if you have the <sys/stat.h> header file. */
#ifndef HOPEST_HAVE_SYS_STAT_H
#define HOPEST_HAVE_SYS_STAT_H 1
#endif

/* Define to 1 if you have the <sys/types.h> header file. */
#ifndef HOPEST_HAVE_SYS_TYPES_H
#define HOPEST_HAVE_SYS_TYPES_H 1
#endif

/* Define to 1 if you have the <unistd.h> header file. */
#ifndef HOPEST_HAVE_UNISTD_H
#define HOPEST_HAVE_UNISTD_H 1
#endif

/* Have we found function adler32_combine. */
#ifndef HOPEST_HAVE_ZLIB
#define HOPEST_HAVE_ZLIB 1
#endif

/* DEPRECATED (use HOPEST_WITH_HDF5 instead) */
#ifndef HOPEST_HDF5
#define HOPEST_HDF5 1
#endif

/* Linker flags */
#ifndef HOPEST_LDFLAGS
#define HOPEST_LDFLAGS " -L/scratch/iagbole/Codeentwicklung/hopest/hdf5/hdf5/lib"
#endif

/* Libraries */
#ifndef HOPEST_LIBS
#define HOPEST_LIBS "-lhdf5_hl -lhdf5 -lpthread -llapack -lblas -lz -lm   "
#endif

/* Define to the sub-directory in which libtool stores uninstalled libraries.
   */
#ifndef HOPEST_LT_OBJDIR
#define HOPEST_LT_OBJDIR ".libs/"
#endif

/* DEPRECATED (use HOPEST_WITH_METIS instead) */
/* #undef METIS */

/* DEPRECATED (use HOPEST_ENABLE_MPI instead) */
/* #undef MPI */

/* DEPRECATED (use HOPEST_ENABLE_MPIIO instead) */
/* #undef MPIIO */

/* Define to 1 if your C compiler doesn't accept -c and -o together. */
/* #undef NO_MINUS_C_MINUS_O */

/* DEPRECATED (use HOPEST_WITH_P4EST instead) */
/* #undef P4EST */

/* Name of package */
#ifndef HOPEST_PACKAGE
#define HOPEST_PACKAGE "hopest"
#endif

/* Define to the address where bug reports for this package should be sent. */
#ifndef HOPEST_PACKAGE_BUGREPORT
#define HOPEST_PACKAGE_BUGREPORT "burstedde@ins.uni-bonn.de"
#endif

/* Define to the full name of this package. */
#ifndef HOPEST_PACKAGE_NAME
#define HOPEST_PACKAGE_NAME "hopest"
#endif

/* Define to the full name and version of this package. */
#ifndef HOPEST_PACKAGE_STRING
#define HOPEST_PACKAGE_STRING "hopest 0.0.18-72aa-dirty"
#endif

/* Define to the one symbol short name of this package. */
#ifndef HOPEST_PACKAGE_TARNAME
#define HOPEST_PACKAGE_TARNAME "hopest"
#endif

/* Define to the home page for this package. */
#ifndef HOPEST_PACKAGE_URL
#define HOPEST_PACKAGE_URL ""
#endif

/* Define to the version of this package. */
#ifndef HOPEST_PACKAGE_VERSION
#define HOPEST_PACKAGE_VERSION "0.0.18-72aa-dirty"
#endif

/* Use builtin getopt */
/* #undef PROVIDE_GETOPT */

/* Use builtin obstack */
/* #undef PROVIDE_OBSTACK */

/* DEPRECATED (use HOPEST_ENABLE_PTHREAD instead) */
/* #undef PTHREAD */

/* DEPRECATED (use HOPEST_WITH_SC instead) */
/* #undef SC */

/* Define to 1 if you have the ANSI C header files. */
#ifndef HOPEST_STDC_HEADERS
#define HOPEST_STDC_HEADERS 1
#endif

/* Version number of package */
#ifndef HOPEST_VERSION
#define HOPEST_VERSION "0.0.18-72aa-dirty"
#endif

/* Package major version */
#ifndef HOPEST_VERSION_MAJOR
#define HOPEST_VERSION_MAJOR 0
#endif

/* Package minor version */
#ifndef HOPEST_VERSION_MINOR
#define HOPEST_VERSION_MINOR 0
#endif

/* Package point version */
#ifndef HOPEST_VERSION_POINT
#define HOPEST_VERSION_POINT 18-72aa-dirty
#endif

/* Define to 1 if BLAS is used */
#ifndef HOPEST_WITH_BLAS
#define HOPEST_WITH_BLAS 1
#endif

/* use the hdf5 library */
#ifndef HOPEST_WITH_HDF5
#define HOPEST_WITH_HDF5 1
#endif

/* Define to 1 if LAPACK is used */
#ifndef HOPEST_WITH_LAPACK
#define HOPEST_WITH_LAPACK 1
#endif

/* enable metis-dependent code */
/* #undef WITH_METIS */

/* path to installed package p4est (optional) */
/* #undef WITH_P4EST */

/* path to installed package sc (optional) */
/* #undef WITH_SC */
 
/* once: _SRC_HOPEST_CONFIG_H */
#endif
