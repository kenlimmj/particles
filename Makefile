CC=gcc
FC=gfortran

ANAFLAGS=-qopt-report=5 -qopt-report-phase=vec -parallel-source-info=2
OPTFLAGS=-fast -xHost -ansi-alias -restrict -mkl -openmp

CFLAGS=-std=c99 -g -pedantic -Wall -Werror
CFLAGS+=#$(OPTFLAGS) $(ANAFLAGS)
LDFLAGS= -fopenmp
LIBS= -lm

#FFLAGS= -O3 -ftree-vectorize -ffast-math -funroll-loops -fomit-frame-pointer -pipe -finit-real=zero -fopenmp
FFLAGS = -g -fbacktrace -Wall -pedantic -Wextra -W -Wno-unused-function -fbounds-check -fopenmp -Wunderflow
.PHONY: all

all: sphc sphfort bsphfort

sphc: sphc.o
	$(CC) -o $@ $^ $(LDFLAGS) $(LIBS)

sphfort: sphfort.o
	$(FC) -o $@ $^ $(LDFLAGS) $(LIBS)

bsphfort: bsphfort.o
	$(FC) -o $@ $^ $(LDFLAGS) $(LIBS)

%.o: %.c
	$(CC) -c $(CFLAGS) $<

%.o: %.f90
	$(FC) -c $(FFLAGS) $<

clean:
	@-rm *.o
	@-rm output/particles*
