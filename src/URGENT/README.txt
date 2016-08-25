Last modification date of this file: Thu Nov 21 13:35:02 CST 2013
Roger J. Dejus (dejus@aps.anl.gov

Issues:
Compiler alignment warnings for common blocks /ANGLE, /SCREEN, and /UND/.

Executables provided for 64-bit architectures only using static linking for
supported platforms.
Linux:       Red Hat Enterprise Linux Workstation release 6.3 (Santiago)

Sun Solaris: Oracle Solaris 10 9/10 s10s_u9wos_14a SPARC
             SunOS frigg 5.10 Generic_147440-27 sun4u sparc
             SUNW,Sun-Blade-2500

Mac OS X:    Snow Leopard, v10.6.8, Darwin Kernel Version 10.8.0

Windows 7:   Enterprise, Service Pack 1
---

Note on urgent version called from xop
======================================


in the line 132 I placed
        IF(ITYPE.EQ.1)WRITE(6,1000)PERIOD,KX,KY,phase,N
                                                ^^^^^
in order to be an output file with the same variables for the
standard and helical undulator. For the helical we still have:
        IF(ITYPE.EQ.2)WRITE(6,1001)PERIOD,kx,KY,PHASE,N

MSR  Oct 1994


SPLITTED ONE LINE (GIVING PROBLEMS IN HP)
expgi.srio<78> diff UU.F urgent.f 
1232,1234c1232
<         IF(IHARM.GT.0.OR.IDEBUG.EQ.1)THEN
<           WRITE(6,1300)XPMM,YPMM,E,DELTAP,DELTAF
<         ENDIF
---
>
IF(IHARM.GT.0.OR.IDEBUG.EQ.1)WRITE(6,1300)XPMM,YPMM,E,DELTAP,DELTAF
The corrected source is urgent.f (not URGENT.F)
MSR 97/10/07


On 17th July 2000 version 1.5 of WS from R. Dejus has been checked in.
This means that the source code is not longer compatible with Macs.
