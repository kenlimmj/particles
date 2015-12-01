CC=gcc
FC=gfortran

ANAFLAGS=-qopt-report=5 -qopt-report-phase=vec -parallel-source-info=2
OPTFLAGS=-fast -xHost -ansi-alias -restrict -mkl -openmp

CFLAGS=-std=c99 -g -pedantic -Wall -Werror
CFLAGS+=#$(OPTFLAGS) $(ANAFLAGS)

FFLAGS= -O3 -ftree-vectorize -ffast-math -funroll-loops -fomit-frame-pointer -pipe -finit-real=zero

.PHONY: all sphc sphfort

all: sphc.o sphfort.o

sphc: sphc.o
	$(CC) -o $@ $^ $(LDFLAGS) $(LIBS)

sphfort: sphfort.o
	$(FC) -o $@ $^ $(LDFLAGS) $(LIBS)

%.o: %.c
	$(CC) -c $(CFLAGS) $<

%.o: %.f90
	$(FC) -c $(FFLAGS) $<
