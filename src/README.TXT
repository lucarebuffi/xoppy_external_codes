Last modifications: 2016-08-15 srio@esrf.eu

Created repository of sources of external codes for XOPPY

Copied and cleaned from xop2.4/src


Compilation: 
------------

cd <DIR>
make -f makefile.unix 


Last modification date of this file: Tue Dec  9 16:41:41 CST 2014
Roger J. Dejus (dejus@aps.anl.gov)

XOP2.4
======

Compilations performed on Unix (Sun Solaris), Linux (HP Workstation), Mac OS X, and Windows 7 (HP desktop):
% make clean   (on Windows % make -f makefile.win clean)
% make target  (on Windows % make -f makefile.win target)
% make install (on Windows % make -f makefile.win install)

This was done for each src/subdirectory.
Compilations were performed for 64-bit environments only.


=== Operating systems and compilers ===

--- Linux Red Hat Enterprise Workstation release 6.3 (Santiago) ---
[dejus@loke ~]$ uname -a
Linux loke.aps.anl.gov 2.6.32-279.14.1.el6.x86_64 #1 SMP Mon Oct 15 13:44:51 EDT 2012 x86_64 x86_64 x86_64 GNU/Linux
[dejus@loke ~]$ ifort --version
ifort (IFORT) 14.0.1 20131008
Copyright (C) 1985-2013 Intel Corporation.  All rights reserved.
[dejus@loke x86]$ gfortran --version

Updated November 2014
Red Hat Enterprise Workstation release 6.5 (Santiago)
[dejus@loke rjd]$ uname -a
Linux loke.aps.anl.gov 2.6.32-431.11.2.el6.x86_64 #1 SMP Mon Mar 3 13:32:45 EST 2014 x86_64 x86_64 x86_64 GNU/Linux
Codes diff_pat.F90 and compliance.F90 from CRYSTAL using
GNU Fortran (GCC) 4.4.7 20120313 (Red Hat 4.4.7-4)
Copyright (C) 2010 Free Software Foundation, Inc.
And codes us.f90, tc.f, and ws.f from DEJUS using ifort above.


--- Mac OS X (Darwin) ---
yabloko:src dejus$ uname -a
Darwin yabloko.aps.anl.gov 10.8.0 Darwin Kernel Version 10.8.0: Tue Jun  7 16:32:41 PDT 2011; root:xnu-1504.15.3~1/RELEASE_X86_64 x86_64
Software version is OS X "Snow Leopard."
yabloko:~ dejus$ sw_vers
ProductName:	Mac OS X
ProductVersion:	10.6.8
BuildVersion:	10K549
yabloko:src dejus$ gfortran --version
GNU Fortran (GCC) 4.6.2 20111019 (prerelease)
Copyright (C) 2011 Free Software Foundation, Inc.

Binaries were also tested on OS X "Mountain Lion" and "Yosemite" and they run fine.
tvdisplay6:~ dejus$ uname -a
Darwin tvdisplay6.aps.anl.gov 12.5.0 Darwin Kernel Version 12.5.0: Sun Sep 29 13:33:47 PDT 2013; root:xnu-2050.48.12~1/RELEASE_X86_64 x86_64
tvdisplay6:~ dejus$ sw_vers
ProductName:    Mac OS X
ProductVersion: 10.8.5
BuildVersion:   12F45
and
tvdisplay6:darwin dejus$ uname -a
Darwin tvdisplay6.aps.anl.gov 14.0.0 Darwin Kernel Version 14.0.0: Fri Sep 19 00:26:44 PDT 2014; root:xnu-2782.1.97~2/RELEASE_X86_64 x86_64
tvdisplay6:darwin dejus$ sw_vers
ProductName:    Mac OS X
ProductVersion: 10.10.1
BuildVersion:   14B25

Updated November 2014 (But only distributed in the src/darwin_yosemite directory; not in the bin.darwin directory which contains
files created under "Snow Leopard.")
Codes diff_pat.F90 and compliance.F90 from CRYSTAL and us.f90, tc.f, and ws.f from DEJUS using
tvdisplay6:~ dejus$ uname -a
Darwin tvdisplay6.aps.anl.gov 14.0.0 Darwin Kernel Version 14.0.0: Fri Sep 19 00:26:44 PDT 2014; root:xnu-2782.1.97~2/RELEASE_X86_64
x86_64
tvdisplay6:~ dejus$ sw_vers
ProductName:    Mac OS X
ProductVersion: 10.10
BuildVersion:   14A389
tvdisplay6:~ dejus$ gfortran --version
GNU Fortran (GCC) 4.9.2 20141029 (prerelease)
Copyright (C) 2014 Free Software Foundation, Inc.


--- Unix Sun Solaris (sun5) ---
-bash-3.2$ uname -a
SunOS frigg 5.10 Generic_147440-27 sun4u sparc SUNW,Sun-Blade-2500
-bash-3.2$ f77 -V
NOTICE: Invoking /opt/SUNWspro/bin/f90 -f77 -ftrap=%none -V
f90: Sun Fortran 95 8.4 SunOS_sparc Patch 128231-02 2009/10/20
Usage: f90 [ options ] files.  Use 'f90 -flags' for details
-bash-3.2$ gfortran --version
GNU Fortran (GCC) 4.5.1
Copyright (C) 2010 Free Software Foundation, Inc.

Updated November 2014
New hardware only; compilers and OS unchanged
Codes diff_pat.F90 and compliance.F90 from CRYSTAL and us.f90, tc.f, and ws.f from DEJUS using
-bash-3.2$ uname -a
SunOS helium 5.10 Generic_150400-11 sun4u sparc SUNW,Sun-Fire-V245


--- Windows 7 (x86) ---
Y:\esrf\oct_2013\xop2.4\src>gfortran --version
GNU Fortran (GCC) 4.8.1
Copyright (C) 2013 Free Software Foundation, Inc.

Updated November 2014
Compiler only, OS unchanged
Codes diff_pat.F90 and compliance.F90 from CRYSTAL and us.f90, tc.f, and ws.f from DEJUS using
Y:\esrf\nov_2014\xop2.4\src>gfortran --version
GNU Fortran (GCC) 4.9.1
Copyright (C) 2014 Free Software Foundation, Inc.
                                                                            

The following table lists source applications compiled for xop2.4 and associated compilers.
The directory entry contains the compiled (Fortran and C) sources.
Notations: if = ifort/icc, gf = gfortran/gcc, sf = sun fortran/sun C

-------------------------------------------------------------------------------
Directory       File         Linux      Sun5      Windows 7/8   Mac OS X
-------------------------------------------------------------------------------
CRYSTAL         diff_pat *)  gf         gf        gf            gf

CRYSTAL_BENT    cryst_ml     if         sf        gf            gf
                cryst_pp     if         sf        gf            gf
                cryst_tt     if         sf        gf            gf
                cryst_za     if         sf        gf            gf

DEJUS           tc *)        if         gf        gf            gf
                us *)        if         gf        gf            gf
                ws *)        if         gf        gf            gf

INPRO           inpro        if         sf        gf            gf

MLAYER          mlayer       if         sf        gf            gf

URGENT          urgent       if         sf        gf            gf

XCOM31          xcom         if         sf        gf            gf

XPOWDER_FML     xpowder_fml  if         gf        gf            gf

XTUBES          xtubes       if         sf        gf            gf

XTUBE_W         tasmip       if         sf        gf            gf

YAUP            bfield       if         sf        gf            gf
                bfield2      if         sf        gf            gf
                txt2u        if         sf        gf            gf
                u2txt        if         sf        gf            gf
                yaup         if         sf        gf            gf
-------------------------------------------------------------------------------

*) Updated November 2014
Codes diff_pat.F90, us.f90, tc.f, and ws.f, were updated using the compilers and hardware listed above.
Code compliance.F90 not listed above because although executable is provided it is not yet part of XOP
and GUI is not provided.

