# set paths and compiler flags
HDF5_VERSION=1.8.13

all: hdf5

hdf5: 
	@# HDF5 libraries
	@ if [ -n "$(HDF5_VERSION)" ] && [ ! -f hdf5-$(HDF5_VERSION)/hdf5/lib/libhdf5.a ] ; then \
		echo "===============================================================";\
		echo "= M A K E   H D F 5                                           =";\
		echo "===============================================================";\
		if [ ! -f hdf5-$(HDF5_VERSION).tar.bz2 ] ; then \
		  wget -N http://www.hdfgroup.org/ftp/HDF5/prev-releases/hdf5-$(HDF5_VERSION)/src/hdf5-$(HDF5_VERSION).tar.bz2 ; \
		fi ; \
		if [ ! -d hdf5-$(HDF5_VERSION) ] ; then \
			tar -xvjf hdf5-$(HDF5_VERSION).tar.bz2 1>/dev/null  ; \
		fi ; \
		cd hdf5-$(HDF5_VERSION) && CC=$(CC) FC=$(FC) ./configure --enable-fortran --disable-shared && $(MAKE) && $(MAKE) install ; \
		echo "===============================================================";\
		echo "= H D F 5   D O N E                                           =";\
		echo "===============================================================";\
	fi
	

# utility targets
.PHONY: clean
clean:  cleanhdf5 
cleanhdf5:
	rm -f *~
	@# HDF5 libraries
	@ if [ -d hdf5-$(HDF5_VERSION) ] ; then \
		cd hdf5-$(HDF5_VERSION) && $(MAKE) clean && $(MAKE) uninstall ; \
	fi
