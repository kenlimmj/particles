CC=icc
FC=mpif90

ANAFLAGS=
OPTFLAGS=-fast -xHost -mkl

CFLAGS=-std=c99 -g -pedantic -Wall -Werror
CFLAGS+=#$(OPTFLAGS) $(ANAFLAGS)
LDFLAGS= -fopenmp
LIBS= -lm

#FFLAGS= -O3 -ftree-vectorize -ffast-math -funroll-loops -fomit-frame-pointer -pipe -finit-real=zero -fopenmp
#FFLAGS = -g -fbacktrace -Wall -pedantic -Wextra -W -Wno-unused-function -fbounds-check -fopenmp -Wunderflow
#FFLAGS = -O3  
FFLAGS = -fcheck=all -g -fbacktrace -Wall -pedantic -Wextra -W -Wno-unused-function -fbounds-check -fopenmp -Wunderflow
.PHONY: all

all: sphc sphfort bsphfort bsphfort_omp

sphc: sphc.o
	$(CC) -o $@ $^ $(LDFLAGS) $(LIBS)

sphfort: sphfort.o
	$(FC) -o $@ $^ $(LDFLAGS) $(LIBS)

bsphfort: bsphfort.o
	$(FC) -o $@ $^ $(LDFLAGS) $(LIBS)

bsphfort_omp: bsphfort_omp.o
	$(FC) -o $@ $^ $(LDFLAGS) $(LIBS)

%.o: %.c
	$(CC) -c $(CFLAGS) $<

%.o: %.f90
	$(FC) -c $(FFLAGS) $<

clean:
	@-rm *.o
	@-rm output/data/particles*
	@-rm output/data/frame*

