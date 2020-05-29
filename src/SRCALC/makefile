FC = gfortran
FFLAGS = -O -m64
FFLAGSD = -DDEBUG1 -g -m64
LDFLAGS = -static-libgfortran -static-libgcc

all:
	$(FC) $(FFLAGS) -o srcalc.o IDPower.f95
	$(FC) $(LDFLAGS) -o srcalc srcalc.o

install:
	/bin/cp -pi srcalc ../../bin.linux/srcalc

clean:
	rm -f *.o srcalc
