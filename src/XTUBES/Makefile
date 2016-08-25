all: xtubes.o mamspec.o
	gcc -o xtubes xtubes.o  mamspec.o -lm
.c.o:
	gcc -c $<
install:
	cp xtubes ../../bin.`./../../xop_uname`/xtubes
clean:
	/bin/rm -f *.o xtubes xtubes_tmp.dat
