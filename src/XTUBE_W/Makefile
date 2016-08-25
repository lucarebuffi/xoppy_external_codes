all: genspec1.o test1.o
	gcc -o tasmip  genspec1.o test1.o -lm
.c.o:
	gcc -c $<
install:
	cp tasmip ../../bin.`./../../xop_uname`/tasmip
clean:
	/bin/rm -f *.o tasmip tasmip_tmp.dat
