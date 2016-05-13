CFLAGS = -O

PROGS = bc_ortho btet chain crect fcc hex2d hex3d \
kpointer portho prect ptet rhom sc square

.c: 
	cc -o $* $(CFLAGS) $*.c -lm

all: $(PROGS)

install: all
	cp $(PROGS) ../bin

clean:
	rm -f $(PROGS)
	




