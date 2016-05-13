#
# Makefile for cubebuilder
#
# by: Robert Snoeberger
#
CC=gcc
CFLAGS= -O3
LIBS = -lm 

COBJS= main.o build.o voxel.o

all: cubebuilder

cubebuilder: $(COBJS) 
	cc -o cubebuilder2 $(CFLAGS) $(COBJS) $(LIBS)

clean:
	rm -f $(COBJS)
