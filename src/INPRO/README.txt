Last modification date of this file: Thu Nov  7 17:04:40 CST 2013
Roger J. Dejus (dejus@aps.anl.gov
Fixed bug in strufacnewm.f which did not have the variables zmass and dwf declared as arrays.
Reordered variables in common blocks in inpro.f, matparmo.f, fresnelmk.f, parametermo.f, strufacnewm.f.

Executables provided for 64-bit architectures only using static linking for supported platforms.
Linux:       Red Hat Enterprise Linux Workstation release 6.3 (Santiago)

Sun Solaris: Oracle Solaris 10 9/10 s10s_u9wos_14a SPARC
             SunOS frigg 5.10 Generic_147440-27 sun4u sparc SUNW,Sun-Blade-2500

Mac OS X:   Snow Leopard, v10.6.8, Darwin Kernel Version 10.8.0

Windows 7:  Enterprise, Service Pack 1
---

Changed for porting to MAC:
  IMAG changed by DIMAG
  print by write(*,*)
  Cleaning common blocks (the variables have the same names)
  Bugs in declaration of DIFFRACTION common block
  Other cleaning
C. Potee+ M. Sanchez del Rio, June 2000

Changed in perfect_crystalmo.f AIMAG by imag
MSR 97/02/04

Small modification for compilation in irix documented in matparmo.f
MSR 97/10/29
