Last modification date of this file: Wed Dec  3 17:01:35 CST 2014
Roger J. Dejus (dejus@aps.anl.gov

Executables provided for 64-bit architectures using static linking for supported platforms.
Results in TEST subdirectories were created with the Linux version.
Insignificant differences are noticed for some cases between different platforms.

Linux:       Red Hat Enterprise Linux Workstation release 6.5 (Santiago)

Sun Solaris: Oracle Solaris 10 9/10 s10s_u9wos_14a SPARC
             SunOS helium 5.10 Generic_150400-11 sun4u sparc SUNW,Sun-Fire-V245
		    
Mac OS X:   Snow Leopard, v10.6.8, Darwin Kernel Version 10.8.0
            Yosemite, v10.10, Darwin Kernel Version 14.0.0
	    (us, ws, and tc only but not installed in the bin.darwin directory; provided
	    in the source distribution in the darwin_yosemite directory)

Windows 7:  Enterprise, Service Pack 1
---

Compilation of the libraries:
=============================
No problems with nr
Two small problems found with ioinit.f for creating the u1 library in Linux. 
1) The automatic declaration is not accepted by g77 (the gnu fortran compiler).
  automatic       iok, fenv, ienv, ename, fname, form, blank
2) The "I" format is not accepted (must be "Iw", w=width).

I solved the compilation by creating a preprocessor condition, i.e. a 
new file ioinit.F that I compile in Linux with the -Dlinux
The new code is:
#if linux
#else
        automatic       iok, fenv, ienv, ename, fname, form, blank
#endif
...
#if linux
 2004   format ('ioinit: cctl=', l12, ', bzro=', l12, ', apnd=', l12)
#else
 2004   format ('ioinit: cctl=', l, ', bzro=', l, ', apnd=', l)
#endif



Compilation of us:
==================
sources us.F , bright.f hunt.f
I realised that the only code used from the u1 library is hunt.f, thus
I copied to the current directory (always HP gave me linking problems).

Compilations in sun4
--------------------
f77 -c -Bstatic -O  us.F etc.
f77 -o us us.o brighte.o hunt.o -lV77

Compilations in sun5
--------------------
same as sun4
Runs OK in the machine where it was compiled (expgh) but gives a problem
in another solaris machine (expgj), which has no fortran compiler (I do
not know why):
expgj.srio<85> ./us
ld.so.1: ./us: fatal: relocation error: symbol not found: date_:
referenced in ./us
Killed

Conclusion: Using the binary created with sun4 for sun5 also.

Compilations in Silicon Graphics
--------------------------------
f77 -c -O   us.F (etc.)
f77 -o us us.o brighte.o hunt.o

Compilations in linux
---------------------
g77 -c -O   us.F (etc.)
g77 -o us us.o brighte.o hunt.o 

us.F gives a warning (but works):
us.F:393: warning:
           OPEN (UNIT=1,FILE=FILE_NAME,STATUS='OLD',ACCESS='SEQUENTIAL',
          ^
Unsupported OPEN control item at (^) -- ACTION=, ASSOCIATEVARIABLE=,
BLOCKSIZE=, BUFFERCOUNT=, CARRIAGECONTROL=, DEFAULTFILE=, DELIM=,
DISPOSE=, EXTENDSIZE=, INITIALSIZE=, KEY=, MAXREC=, NOSPANBLOCKS,
ORGANIZATION=, PAD=, POSITION=, READONLY=, RECORDTYPE=, SHARED=, and
USEROPEN= are not supported


Compilations in hp800 (HP10) and hp700 (HP9)
--------------------------------------------
f77 -c +U77 -K +E1 +E4   -O us.f (etc.)
f77 -o us us.o brighte.o hunt.o +U77


Compilation of tc:
==================
sources brighte.f tc.d usb.f libu1

Compilation in sun5 and sun4
----------------------------
OK using your makefile_solaris (slightly modify to reflect my dirs and
using -Bstatic ).
f77 -c -O -Bstatic  brighte.f (etc.)
f77 -o tc -Bstatic tc.o usb.o brighte.o -L../../lib.sun5 -lu1 -lV77


Compilation in Silicon Graphics
-------------------------------
f77 -c -O  tc.f
f77 -c -O  usb.f
f77 -c -O  brighte.f
f77 -o tc  tc.o usb.o brighte.o -L../../lib.irix -lu1 

Compilation in Linux
--------------------
Same warning as with us:
expglin.srio<147> gmake -f makefile_linux 
g77 -c -O  tc.f
tc.f: In program `tc':
tc.f:239: warning:
           open (unit=1,file='tc.dat',status='old',access='sequential',
           ^
Unsupported OPEN control item at (^) -- ACTION=, ASSOCIATEVARIABLE=,
BLOCKSIZE=, BUFFERCOUNT=, CARRIAGECONTROL=, DEFAULTFILE=, DELIM=,
DISPOSE=, EXTENDSIZE=, INITIALSIZE=, KEY=, MAXREC=, NOSPANBLOCKS,
ORGANIZATION=, PAD=, POSITION=, READONLY=, RECORDTYPE=, SHARED=, and
USEROPEN= are not supported
g77 -c -O  usb.f
g77 -c -O  brighte.f
g77 -o tc  tc.o usb.o brighte.o -L../../lib.linux -lu1 


Compilation in hp800 and hp700
------------------------------
Problems linking with the library. Solved by copying hunt.f 
f77 -c -O +U77 +E1 +E4   +es  hunt.f (etc.)
fort77 -o tc  tc.o usb.o brighte.o hunt.o -L../../lib.hp800 -lu1 -lU77

Conclusions:
============
Except problems of linking the only recommendation is to fix ioinit.f
(remove "automatic" and use format I12) that are not f77 standard commands.

Managing PC and Unix source versions
====================================
From your last email I see that you have to created two versions of your
codes: one for Unix and the other for Lahey. I found that problem also,
and I know is a pain to have two different versions. What I did is to
create a master version in Unix, and use the preprocessor to create the
Lahey version, that I copy and compile in my PC. 
The idea is the following:
create a file myprog.F (instead of myprog.f) with conditional code as:

...
#if lahey
             mu =  2.0D0 * twopi *(-aimag(refrac)) / (lambda*1.0d-8)
#else
             mu =  2.0D0 * twopi *(-dimag(refrac)) / (lambda*1.0d-8)
#endif
...

To compile in Unix, use "f77 -c myprog.F" and  the lahey part is ignored
(unless you used the -Dlahey flag). 
To create the Lahey version to be copied to the PC just do:
   $(CPP) -P -Dlahey myprog.F myprog.flahey
Where $(CPP) is the preprocessor (usually cpp. In hp-ux is /lib/cpp).
That creates the myprog.flahey that I move to my PC and rename to
myprog.f. In this way you have to keep inly a single master source,
although it is difficult to debug in PC, because you need the Unix
machine to create the new code.

