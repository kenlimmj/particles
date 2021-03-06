CC=icc
FC=ifort

ANAFLAGS=-qopt-report=5 -qopt-report-phase=vec -parallel-source-info=2
OPTFLAGS=-fast -xHost -ansi-alias -restrict -mkl

CFLAGS=-std=c99 -g -pedantic -Wall -Werror
CFLAGS+=#$(OPTFLAGS) $(ANAFLAGS)
LDFLAGS= -fopenmp
LIBS= -lm

#FFLAGS= -O3 -ftree-vectorize -ffast-math -funroll-loops -fomit-frame-pointer -pipe -finit-real=zero -fopenmp
#FFLAGS = -g -fbacktrace -Wall -pedantic -Wextra -W -Wno-unused-function -fbounds-check -fopenmp -Wunderflow
FFLAGS = -O3 -xHost -ip -fast -qopt-report=5 -qopt-report-phase=vec -align -fopenmp
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
