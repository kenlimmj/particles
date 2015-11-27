
CC=icc
CFLAGS=-std=c99 -g -O3

FC=ifort
FFLAGS=

.PHONY: all

all: sphc.o sphf.o

sphc: sphc.o
	$(CC) -o $@ $^ $(LDFLAGS) $(LIBS)

sphf: sphf.o
	$(FC) -o $@ $^ $(LDFLAGS) $(LIBS)

%.o: %.c
	$(CC) -c $(CFLAGS) $<

%.o: %.f
	$(FC) -c $(FFLAGS) $<

