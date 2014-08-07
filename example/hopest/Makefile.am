
## Please keep this if condition.
if HOPEST_WITH_HDF5

bin_PROGRAMS += \
    example/hopest/hopest

example_hopest_hopest_SOURCES = example/hopest/hopest.f90

example_hopest_hopest_CPPFLAGS = $(HOPEST_CPPFLAGS)
example_hopest_hopest_FCFLAGS = $(HOPEST_CPPFLAGS) $(AM_FCFLAGS) $(FCFLAGS)

endif
