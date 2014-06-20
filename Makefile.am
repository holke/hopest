
# This file is part of con4m.
# Makefile.am in toplevel directory

ACLOCAL_AMFLAGS = -I config @CON4M_SC_AMFLAGS@ @CON4M_P4EST_AMFLAGS@
if CON4M_SC_MK_USE
@CON4M_SC_MK_INCLUDE@
endif
if CON4M_P4EST_MK_USE
@CON4M_P4EST_MK_INCLUDE@
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

# use this if you want to link in con4m without autotools
sysconf_DATA += Makefile.con4m.mk
CLEANFILES += Makefile.con4m.mk
Makefile.con4m.mk : Makefile.con4m.pre
	cat $< | \
        sed -e 's,{\(\(.*prefix\|sysconfdir\)\)},{con4m_\1},g' \
            -e 's,^\(\(.*prefix\|sysconfdir\) *=\),con4m_\1,g' > $@

# install con4m m4 macros in the correct directory
con4maclocaldir = $(datadir)/aclocal
dist_con4maclocal_DATA = config/con4m_include.m4

# setup test environment
if CON4M_MPIRUN
LOG_COMPILER = @CON4M_MPIRUN@
AM_LOG_FLAGS = @CON4M_MPI_TEST_FLAGS@
endif

# recurse only into subpackages
SUBDIRS = @CON4M_SC_SUBDIR@ @CON4M_P4EST_SUBDIR@
DIST_SUBDIRS = $(SUBDIRS)

# handle toplevel directory
EXTRA_DIST += \
        bootstrap con4mindent build-aux/git-version-gen build-aux/git2cl doc
DISTCLEANFILES += \
        _configs.sed src/con4m_config.h @CON4M_DISTCLEAN@
.PHONY: lint ChangeLog

# non-recursive build
include src/Makefile.am
include test/Makefile.am
include example/poisson/Makefile.am
include example/darcy/Makefile.am

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
if CON4M_DIST_DENY
	@echo "-----------------------------------------------------"
	@echo "make dist does not work with external subpackages"
	@echo "-----------------------------------------------------"
	rm -rf $(distdir)
	@exit 1
endif

clean-local:
	rm -f ChangeLog *vtu *.visit
