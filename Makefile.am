
# This file is part of hopest.
# Makefile.am in toplevel directory

ACLOCAL_AMFLAGS = -I config @HOPEST_SC_AMFLAGS@ @HOPEST_P4EST_AMFLAGS@
if HOPEST_SC_MK_USE
@HOPEST_SC_MK_INCLUDE@
endif
if HOPEST_P4EST_MK_USE
@HOPEST_P4EST_MK_INCLUDE@
endif

# initialize empty variables
AM_CPPFLAGS =
CLEANFILES =
DISTCLEANFILES =
EXTRA_DIST =
LDADD =
LINT_CSOURCES =
TESTS =
bin_PROGRAMS =
check_PROGRAMS =
include_HEADERS =
lib_LTLIBRARIES =
nodist_include_HEADERS =
noinst_HEADERS =
noinst_PROGRAMS =
sysconf_DATA =

# use this if you want to link in hopest without autotools
sysconf_DATA += Makefile.hopest.mk
CLEANFILES += Makefile.hopest.mk
Makefile.hopest.mk : Makefile.hopest.pre
	cat $< | \
        sed -e 's,{\(\(.*prefix\|sysconfdir\)\)},{hopest_\1},g' \
            -e 's,^\(\(.*prefix\|sysconfdir\) *=\),hopest_\1,g' > $@

# install hopest m4 macros in the correct directory
hopestaclocaldir = $(datadir)/aclocal
dist_hopestaclocal_DATA = \
        config/hopest_include.m4 config/hopest_hdf5.m4

# setup test environment
if HOPEST_MPIRUN
LOG_COMPILER = @HOPEST_MPIRUN@
AM_LOG_FLAGS = @HOPEST_MPI_TEST_FLAGS@
endif

# recurse only into subpackages
SUBDIRS = @HOPEST_SC_SUBDIR@ @HOPEST_P4EST_SUBDIR@
DIST_SUBDIRS = $(SUBDIRS)

# handle toplevel directory
EXTRA_DIST += \
        bootstrap hopestindent build-aux/git-version-gen build-aux/git2cl doc
DISTCLEANFILES += \
        _configs.sed src/hopest_config.h @HOPEST_DISTCLEAN@
.PHONY: lint ChangeLog

# non-recursive build
include src/Makefile.am
include test/Makefile.am
include example/hopest/Makefile.am
include example/hopest_p4est_geom/Makefile.am
include deplist.mk

# automatically builds a dependency list deplist.mk from srcfiles.mk uing the tool builddeps
# TODO automate the FC flags "-cpp -E -w" works only with GNU compiler 
Builddeps:
	srcfiles=@top_srcdir@/srcfiles.mk; \
	rm -f @top_srcdir@/srcfiles_pp.mk; \
	for filename in $$(cat $$srcfiles); \
	    do \
	    echo @top_srcdir@/$$filename; \
		outfile=$$(echo $$filename | sed s/f90/i90/) ; \
		echo @top_srcdir@/$$outfile | sed -e 's/ /\n/g' >> @top_srcdir@/srcfiles_pp.mk; \
		$(FC) -cpp -E -w -I@top_srcdir@/src @top_srcdir@/$$filename -o @top_srcdir@/$$outfile ; \
	done;
	@python @top_srcdir@/tools/builddeps --source=@top_srcdir@/srcfiles_pp.mk @top_srcdir@/deplist.mk
	
# revision control and ChangeLog
ChangeLog:
	(GIT_DIR=@top_srcdir@/.git git log > .ChangeLog.tmp && \
         cat .ChangeLog.tmp | @top_srcdir@/build-aux/git2cl > ChangeLog) ; \
        rm -f .ChangeLog.tmp

dist-hook:
	echo $(VERSION) > $(distdir)/.tarball-version
	test "x$(VERSION)" = "x`@top_srcdir@/build-aux/git-version-gen \
          @top_srcdir@/.tarball-version`" || \
        ((echo "Stale version;" ; echo "Please run:" ; \
          echo "     cd @top_srcdir@ && ./bootstrap" ; \
          echo "before make dist") 1>&2 ; rm -r $(distdir) ; exit 1)
if HOPEST_DIST_DENY
	@echo "-----------------------------------------------------"
	@echo "make dist does not work with external subpackages"
	@echo "-----------------------------------------------------"
	rm -rf $(distdir)
	@exit 1
endif

install-exec-local:
	cp -f @top_builddir@/*.mod $(DESTDIR)$(includedir)

uninstall-local:
	rm -f $(DESTDIR)$(includedir)/*.mod

mostlyclean-local:
	rm -f @top_builddir@/*.mod

clean-local:
	rm -f ChangeLog *vtu *.visit
