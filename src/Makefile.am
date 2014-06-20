
# This file is part of hopest.
# Makefile.am in src
# included non-recursively from toplevel directory

libcon4m_generated_headers = src/hopest_config.h
libcon4m_installed_headers = \
        src/hopest.h
libcon4m_internal_headers =
libcon4m_compiled_sources = \
         src/hopest.c src/hopest_f.f

libcon4m_installed_headers += src/hopest_hdf5.h
if CON4M_WITH_HYPRE
libcon4m_compiled_sources += src/hopest_hdf5.c
endif

# this variable is used for headers that are not publicly installed
CON4M_CPPFLAGS =

lib_LTLIBRARIES += src/libhopest.la
src_libcon4m_la_SOURCES = \
	$(libcon4m_internal_headers) \
	$(libcon4m_compiled_sources)
src_libcon4m_la_CPPFLAGS = $(AM_CPPFLAGS) $(CON4M_CPPFLAGS)
LDADD += @top_builddir@/src/libhopest.la @CON4M_P4EST_LDADD@ @CON4M_SC_LDADD@
EXTRA_src_libcon4m_la_DEPENDENCIES = @CON4M_P4EST_LDADD@ @CON4M_SC_LDADD@

nodist_include_HEADERS += $(libcon4m_generated_headers)
include_HEADERS += $(libcon4m_installed_headers)

AM_CPPFLAGS += -I@top_srcdir@/src @CON4M_SC_CPPFLAGS@ @CON4M_P4EST_CPPFLAGS@

LINT_CSOURCES += $(libcon4m_compiled_sources)
