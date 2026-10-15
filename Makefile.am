# this is example-file: Makefile.am

# the subdirectories of the project to go into
SUBDIRS =   \
            src

check-local:
	$(SHELL) $(top_srcdir)/TESTS/TEST_OUTPUT_PATHS/run_test.sh \
	  "$(FC)" "$(OPENMP_FCFLAGS)" "$(abs_top_builddir)/src"
	$(SHELL) $(top_srcdir)/TESTS/run_tests.sh \
	  "$(abs_top_builddir)/src/IMEX_SfloW2D"
