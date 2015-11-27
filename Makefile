CC=icc
FC=ifort

ANAFLAGS=-qopt-report=5 -qopt-report-phase=vec -parallel-source-info=2
OPTFLAGS=-fast -xHost -ansi-alias -restrict -mkl -openmp

CFLAGS=-std=c99 -g -pedantic -Wall -Werror
CFLAGS+=$(OPTFLAGS) $(ANAFLAGS)

FFLAGS=

.PHONY: all sphc sphf

all: sphc.o sphf.o

sphc: sphc.o
	$(CC) -o $@ $^ $(LDFLAGS) $(LIBS)

sphf: sphf.o
	$(FC) -o $@ $^ $(LDFLAGS) $(LIBS)

%.o: %.c
	$(CC) -c $(CFLAGS) $<

%.o: %.f
	$(FC) -c $(FFLAGS) $<
