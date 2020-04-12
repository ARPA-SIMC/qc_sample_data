FC= gfortran
FCFLAGS = -fbounds-check -g -O0
#FCFLAGS = -O2

LDLIBS= -lsim_base -lsim_vol7d -llog4fortran -lsim_qc -lodbc -lsim_volgrid6d -lm -llog4c
GRIBAPILDLIBS=$(shell pkg-config --exists grib_api_f90 && pkg-config --libs grib_api_f90 || pkg-config --libs eccodes_f90)
LIB_NETCDF= -lnetcdff -lnetcdf 

FCFLAGS += -I/usr/include -I/usr/lib64/gfortran/modules

PROGRAMS = netcdf2grib1_cum netcdf2grib2_cum 

all: $(PROGRAMS)

%.o: %.f90
	$(FC) -c $(FCFLAGS) $(LIB_NETCDF) $< 

%: %.o
	$(FC) $(LDLIBS) $(GRIBAPILDLIBS) $(LIB_NETCDF) -o $@ $^

clean:
	$(RM) $(PROGRAMS) *.o *.mod *~

