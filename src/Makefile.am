
# This file is part of hopest.
# Makefile.am in src
# included non-recursively from toplevel directory

libhopest_generated_headers = src/hopest_config.h
libhopest_installed_headers = \
        src/hopest.h
libhopest_internal_headers = \
	src/hopest_f.h \
	src/p4est/p4est_wrapper.h \
	src/p4est/p4est_HO_geometry.h
libhopest_compiled_sources = \
	src/globals/preprocessing.f90 \
	src/globals/globals.f90 \
	src/mesh/mesh_vars.f90 \
	src/analyze/analyze_vars.f90 \
	src/p4est/p4est_vars.f90 \
	src/output/output_vars.f90 \
	src/refine/refine_vars.f90 \
	src/mpi/mpi.f90 \
	src/readintools/isovaryingstring.f90 \
	src/readintools/readintools.f90 \
	src/interpolation/changeBasis.f90 \
	src/interpolation/basis.f90 \
	src/p4est/p4est_binding.f90 \
	src/p4est/p4est_wrapper.c \
	src/p4est/p4est.f90 \
	src/mesh/prepare_mesh.f90 \
	src/analyze/analyze.f90 \
	src/output/output.f90 \
	src/refine/refine_bindings.f90 \
	src/refine/refine_wrapper.c \
	src/refine/refine.f90

## Please keep this if condition.  Files that rely on HDF5 go here.
if HOPEST_WITH_HDF5
libhopest_compiled_sources += \
	src/io_hdf5/io_hdf5.f90 \
	src/io_hdf5/hdf5_input.f90 \
	src/io_hdf5/hdf5_output.f90 \
	src/output/output_hdf5.f90 \
	src/mesh/mesh.f90 \
	src/mesh/mesh_readin.f90 \
	src/hopest_solver.f90 \
	src/hopest_mesh.f90 \
	src/p4est/p4est_HO_geometry.c \
	src/fortranwrap.f90
endif

# this variable is used for headers that are not publicly installed
HOPEST_CPPFLAGS = $(AM_CPPFLAGS) $(CPPFLAGS) -I@top_builddir@/src

lib_LTLIBRARIES += src/libhopest.la
src_libhopest_la_SOURCES = \
	$(libhopest_internal_headers) \
	$(libhopest_compiled_sources)
src_libhopest_la_CPPFLAGS = $(HOPEST_CPPFLAGS)
src_libhopest_la_FCFLAGS = $(HOPEST_CPPFLAGS) $(AM_FCFLAGS) $(FCFLAGS)
LDADD += @top_builddir@/src/libhopest.la @HOPEST_P4EST_LDADD@ @HOPEST_SC_LDADD@
EXTRA_src_libhopest_la_DEPENDENCIES = @HOPEST_P4EST_LDADD@ @HOPEST_SC_LDADD@

nodist_include_HEADERS += $(libhopest_generated_headers)
include_HEADERS += $(libhopest_installed_headers)

AM_CPPFLAGS += -I@top_srcdir@/src @HOPEST_SC_CPPFLAGS@ @HOPEST_P4EST_CPPFLAGS@

LINT_CSOURCES += $(libhopest_compiled_sources)
