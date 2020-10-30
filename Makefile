# FORTRAN MAKEFILE
#
#
EXEC = test1.exe
FC = mpiifort


 ifneq (,$(filter $(DEBUG), TRUE true ))
      FFLAGS = -O3 -fp-model precise -xHost #-qopenmp
 else
      FFLAGS  = -O0 -g -check all -fpe0 -fno-alias -ftrapuv -traceback
 endif
 
MPI_COMPILE_FLAGS = /usr/local/apps/intel/parallel_studio_xe_2018/impi/2018.0.128/intel64/include
MPI_LINK_FLAGS = -L/usr/local/apps/intel/parallel_studio_xe_2018/impi/2018.0.128/intel64/lib -lmpi
NETCDF_COMPILE_FLAGS = $(shell nc-config --fflags)

NETCDF_COMPILE_FLAGS = $(shell nf-config --fflags)
LIB2 = /work/MOD3DEV/fsidi/temp/nsu/ioapi-3.2
include_path = -I $(LIB2)/Linux2_x86_64ifort \
               -I $(LIB2)/fixed_src  \
               -I $(MPI_COMPILE_FLAGS) \
               $(NETCDF_COMPILE_FLAGS) -I .


IOAPI = -L$(LIB2)/Linux2_x86_64ifort -lioapi
NETCDF = $(shell nf-config --flibs )  
LIBRARIES = $(IOAPI) $(NETCDF) $(MPI_LINK_FLAGS)

OBJECTS = pe_util.o driver.o

$(EXEC) : $(OBJECTS) 
	$(FC) $(OBJECTS) $(LIBRARIES) -o $@
%.o : %.f90
	$(FC) -c $(FFLAGS) $(include_path) $<

clean: 
	rm -f $(OBJECTS) $(EXEC)
