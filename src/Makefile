# AIX
#FC = xlf90
#FFLAG = -g -O3 -qstrict -C -qfullpath
#LIB = -L$(NETCDF)/lib -lnetcdf -lnetcdff

# Linux Intel Compiler
#NETCDF = /home/mm5/netcdf-3.5.1_intel
#FC = ifort
#FFLAG = -C -g -O3 -assume byterecl -fp-model strict 
INC = -I$(NETCDF)/include

# Linux gfortran
NETCDF = /usr
FC = gfortran
INC = -I/usr/lib64/gfortran/modules 
FFLAG = -g -O3 -ffree-line-length-none
LIB = -L$(NETCDF)/lib -lnetcdf -lnetcdff

# Linux PGI Compiler
#NETCDF = /home/mm5/netcdf-3.5.1_pgi
#NETCDF = /home/wrf/netcdf-3.6.1_pgi
#FC = pgf90
#FFLAG = -g -O3 -byteswapio -C

#LIB = -L$(NETCDF)/lib -lnetcdf
OBJS = module_radar_data.o module_bytes.o        module_libradar.o \
       module_date_pack.o  module_radar_site.o   module_string.o \
       module_wrf_data.o   module_map_utils.o    module_radar_qc.o \
       module_config.o     module_radar_region.o module_da_tools.o \
       module_grd.o        module_vdras_prep.o   module_radar_3dv.o

TARGET  = read_cinrad proc_radar
PACKAGE = read_cinrad

.SUFFIXES: 
.SUFFIXES: .f90 .f .o


all:	$(TARGET)
clean:
	rm -f *.o *.mod *~ .*~ $(TARGET)
proc_radar:	proc_radar.o $(OBJS)	
	$(FC) $(FFLAG) $(OBJS) $@.o -o $@ $(LIB)
proc_radar.o:	$(OBJS)

read_cinrad:	read_cinrad.o $(OBJS)	
	$(FC) $(FFLAG) $(OBJS) $@.o -o $@ $(LIB)
read_cinrad.o:	$(OBJS)
module_radar_data.o:	module_bytes.o module_libradar.o module_date_pack.o
module_radar_site.o:	module_string.o
module_wrf_data.o:	module_map_utils.o module_libradar.o
module_radar_qc.o:	module_radar_data.o module_libradar.o module_config.o module_radar_region.o
module_radar_region.o:	module_grd.o module_libradar.o
module_da_tools.o:	module_wrf_data.o module_map_utils.o module_config.o	
module_radar_3dv.o:	module_radar_data.o module_radar_qc.o module_libradar.o module_config.o

.f90.o:
	$(FC) $(INC) $(FFLAG) -c $<
.f.o:
	$(FC) $(INC) $(FFLAG) -c $<

