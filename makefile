# makefile for psvbpcat
PROGS = mpi-psvbpcat
FC	= mpif90
FC2=ifort
FFLAGS	= -O
#debug
#FFLAGS = -O2 -heap-arrays -g -check all -fpe0 -warn -traceback -debug extended

SRC = calmat.f  dcsymbdl.f dcsymbdl3.f  trialf.f glu2.f rk3.f others.f formpi.f mpi-psvbpcat.f
SRC2 = calmat.f dcsymbdl.f dcsymbdl3.f trialf.f rk3.f others.f glu2.f psvbpcat.f
OBJS = $(SRC:.f=.o)
OBJS2 = $(SRC2:.f=.o)
.SUFFIXES: .f .o

all:$(PROGS) psvbpcat

mpi-psvbpcat: $(OBJS)
	$(FC) $(FFLAGS) -o $@ $(OBJS) 

psvbpcat: $(OBJS2)
	$(FC2) $(FFLAGS) -o $@ $(OBJS2)

mpi-psvbpcat.o: mpi-psvbpcat.f
	$(FC) $(FFLAGS) -c mpi-psvbpcat.f -o $@


.f.o:
	$(FC2) $(FFLAGS) -c $< 

clean:
	rm -f $(OBJS) $(OBJS2) $(PROGS) mpi-psvbpcat psvbpcat workpsvbp
