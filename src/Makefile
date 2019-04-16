# AIX
#FC = xlf90
#FFLAG = -g -O3 -qstrict -C -qfullpath
#LIB = -L$(NETCDF)/lib -lnetcdf -lnetcdff

# Linux Intel Compiler
#NETCDF = /home/mm5/netcdf-3.5.1_intel
#FC = ifort
#FFLAG = -C -g -O3 -assume byterecl -fp-model strict 

# Linux gfortran
NETCDF = /usr
FC = gfortran
FFLAG = -g -O3 -ffree-line-length-none
LIB = -L$(NETCDF)/lib -lnetcdf -lnetcdff

# Linux PGI Compiler
#NETCDF = /home/mm5/netcdf-3.5.1_pgi
#NETCDF = /home/wrf/netcdf-3.6.1_pgi
#FC = pgf90
#FFLAG = -g -O3 -byteswapio -C

INC = -I$(NETCDF)/include
#LIB = -L$(NETCDF)/lib -lnetcdf
OBJS = module_radar_data.o module_bytes.o module_libradar.o \
       module_date_pack.o module_radar_site.o module_string.o

TARGET  = read_cinrad
PACKAGE = read_cinrad

.SUFFIXES: 
.SUFFIXES: .f90 .f .o


all:	$(TARGET)
clean:
	rm -f *.o *.mod  $(TARGET)
test:	test.o $(OBJS)	
	$(FC) $(FFLAG) $(OBJS) $@.o -o $@ $(LIB)
test.o:	$(OBJS)

read_cinrad:	read_cinrad.o $(OBJS)	
	$(FC) $(FFLAG) $(OBJS) $@.o -o $@ $(LIB)
read_cinrad.o:	$(OBJS)
module_radar_data.o:	module_bytes.o module_libradar.o module_date_pack.o module_string.o
module_radar_site.o:	module_string.o


.f90.o:
	$(FC) $(INC) $(FFLAG) -c $<
.f.o:
	$(FC) $(INC) $(FFLAG) -c $<

