### Deplist auto-generated by builddeps V0.99 ###

../../source/hopest/./example/hopest/hopest.o: ../../source/hopest/./example/hopest/hopest.f90 \
	../../source/hopest/./src/globals/globals.o \
	../../source/hopest/./src/io_hdf5/io_hdf5.o \
	../../source/hopest/./src/mesh/mesh.o \
	../../source/hopest/./src/mpi/mpi.o \
	../../source/hopest/./src/readintools/readintools.o 

../../source/hopest/./src/mesh/mesh.o: ../../source/hopest/./src/mesh/mesh.f90 \
	../../source/hopest/./src/cwrapper/p4est_binding.o \
	../../source/hopest/./src/globals/globals.o \
	../../source/hopest/./src/interpolation/basis.o \
	../../source/hopest/./src/mesh/mesh_readin.o \
	../../source/hopest/./src/mesh/mesh_refine.o \
	../../source/hopest/./src/mesh/mesh_vars.o \
	../../source/hopest/./src/mesh/meshfromp4est.o \
	../../source/hopest/./src/output/output_hdf5.o \
	../../source/hopest/./src/output/output_vars.o \
	../../source/hopest/./src/readintools/readintools.o 

../../source/hopest/./example/hello/hell.o.f90: ../../source/hopest/./example/hello/hell.f90.f90 

../../source/hopest/./example/hello/hell.o.o: ../../source/hopest/./example/hello/hell.f90.f90 

../../source/hopest/./src/readintools/readintools.o: ../../source/hopest/./src/readintools/readintools.f90 \
	../../source/hopest/./src/globals/globals.o \
	../../source/hopest/./src/readintools/isovaryingstring.o 

../../source/hopest/./src/readintools/isovaryingstring.o: ../../source/hopest/./src/readintools/isovaryingstring.f90 

../../source/hopest/./src/hopest_c_tes.o.f90: ../../source/hopest/./src/hopest_c_tes.f90.f90 

../../source/hopest/./src/cwrapper/p4est_binding_types.o: ../../source/hopest/./src/cwrapper/p4est_binding_types.f90 

../../source/hopest/./src/cwrapper/p4est_binding.o: ../../source/hopest/./src/cwrapper/p4est_binding.f90 

../../source/hopest/./src/mpi/mpi.o: ../../source/hopest/./src/mpi/mpi.f90 \
	../../source/hopest/./src/globals/globals.o 

../../source/hopest/./src/hopest_c_tes.o.o: ../../source/hopest/./src/hopest_c_tes.f90.f90 

../../source/hopest/./src/mesh/meshfromp4est.o: ../../source/hopest/./src/mesh/meshfromp4est.f90 \
	../../source/hopest/./src/cwrapper/p4est_binding.o \
	../../source/hopest/./src/globals/globals.o \
	../../source/hopest/./src/interpolation/basis.o \
	../../source/hopest/./src/interpolation/changeBasis.o \
	../../source/hopest/./src/mesh/mesh_vars.o \
	../../source/hopest/./src/output/output_vars.o 

../../source/hopest/./src/mesh/mesh_refine.o: ../../source/hopest/./src/mesh/mesh_refine.f90 \
	../../source/hopest/./src/cwrapper/p4est_binding.o \
	../../source/hopest/./src/globals/globals.o \
	../../source/hopest/./src/mesh/mesh_vars.o \
	../../source/hopest/./src/readintools/readintools.o 

../../source/hopest/./src/mesh/mesh_vars.o: ../../source/hopest/./src/mesh/mesh_vars.f90 \
	../../source/hopest/./src/cwrapper/p4est_binding_types.o 

../../source/hopest/./src/mesh/mesh_readin.o: ../../source/hopest/./src/mesh/mesh_readin.f90 \
	../../source/hopest/./src/cwrapper/p4est_binding.o \
	../../source/hopest/./src/globals/globals.o \
	../../source/hopest/./src/io_hdf5/hdf5_input.o \
	../../source/hopest/./src/mesh/mesh_vars.o 

../../source/hopest/./src/globals/globals.o: ../../source/hopest/./src/globals/globals.f90 

../../source/hopest/./src/globals/preprocessing.o: ../../source/hopest/./src/globals/preprocessing.f90 

../../source/hopest/./src/output/output_hdf5.o: ../../source/hopest/./src/output/output_hdf5.f90 \
	../../source/hopest/./src/globals/globals.o \
	../../source/hopest/./src/io_hdf5/hdf5_output.o \
	../../source/hopest/./src/io_hdf5/io_hdf5.o \
	../../source/hopest/./src/mesh/mesh_vars.o 

../../source/hopest/./src/output/output_vars.o: ../../source/hopest/./src/output/output_vars.f90 

../../source/hopest/./src/output/output.o: ../../source/hopest/./src/output/output.f90 \
	../../source/hopest/./src/globals/globals.o \
	../../source/hopest/./src/globals/preprocessing.o \
	../../source/hopest/./src/output/output_vars.o \
	../../source/hopest/./src/readintools/readintools.o 

../../source/hopest/./src/interpolation/changeBasis.o: ../../source/hopest/./src/interpolation/changeBasis.f90 

../../source/hopest/./src/interpolation/basis.o: ../../source/hopest/./src/interpolation/basis.f90 \
	../../source/hopest/./src/globals/globals.o \
	../../source/hopest/./src/globals/preprocessing.o 

../../source/hopest/./src/io_hdf5/io_hdf5.o: ../../source/hopest/./src/io_hdf5/io_hdf5.f90 \
	../../source/hopest/./src/globals/globals.o \
	../../source/hopest/./src/readintools/readintools.o 

../../source/hopest/./src/io_hdf5/hdf5_input.o: ../../source/hopest/./src/io_hdf5/hdf5_input.f90 \
	../../source/hopest/./src/globals/globals.o \
	../../source/hopest/./src/io_hdf5/io_hdf5.o 

../../source/hopest/./src/io_hdf5/hdf5_output.o: ../../source/hopest/./src/io_hdf5/hdf5_output.f90 \
	../../source/hopest/./src/io_hdf5/io_hdf5.o \
	../../source/hopest/./src/output/output_vars.o 

../../source/hopest/./src/fortranwrap.o: ../../source/hopest/./src/fortranwrap.f90 \
	../../source/hopest/./src/mesh/mesh.o \
	../../source/hopest/./src/readintools/readintools.o 

