
# This file is part of hopest.
# Makefile.am in src
# included non-recursively from toplevel directory

libhopest_generated_headers = src/hopest_config.h
libhopest_installed_headers = \
        src/hopest.h \
	src/hopest_fortran_test.h
libhopest_internal_headers =
libhopest_compiled_sources = \
         src/hopest.c src/hopest_f.f \
	 src/hopest_fortran_test.c src/hopest_c_test.f src/hopest_c_testf90.f90

libhopest_installed_headers += src/hopest_hdf5.h
if HOPEST_WITH_HDF5
libhopest_compiled_sources += src/hopest_hdf5.c
endif

# this variable is used for headers that are not publicly installed
HOPEST_CPPFLAGS =

lib_LTLIBRARIES += src/libhopest.la
src_libhopest_la_SOURCES = \
	$(libhopest_internal_headers) \
	$(libhopest_compiled_sources)
src_libhopest_la_CPPFLAGS = $(AM_CPPFLAGS) $(HOPEST_CPPFLAGS)
LDADD += @top_builddir@/src/libhopest.la @HOPEST_P4EST_LDADD@ @HOPEST_SC_LDADD@
EXTRA_src_libhopest_la_DEPENDENCIES = @HOPEST_P4EST_LDADD@ @HOPEST_SC_LDADD@

nodist_include_HEADERS += $(libhopest_generated_headers)
include_HEADERS += $(libhopest_installed_headers)

AM_CPPFLAGS += -I@top_srcdir@/src @HOPEST_SC_CPPFLAGS@ @HOPEST_P4EST_CPPFLAGS@

LINT_CSOURCES += $(libhopest_compiled_sources)
