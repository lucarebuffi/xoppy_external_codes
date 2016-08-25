Last modification date of this file: Fri Nov 22 17:01:58 CST 2013
Roger J. Dejus (dejus@aps.anl.gov

Issues:
In a 64-bit Solaris environment, many system libraries are available only as
shared dynamic libraries. These include libm.so and libc.so (libm.a and libc.a
are not provided). This means that -Bstatic and -dn may cause linking errors
in 64-bit Solaris environments. Applications must link with the dynamic
libraries in these cases. Hence, for the Sun Solaris the dynamic linking 
was employed. It uses the 
libm.so.2 =>	 /lib/64/libm.so.2
libc.so.1 =>	 /lib/64/libc.so.1
libraries but those are actually the same being used with the GCC compiler.

Executables provided for 64-bit architectures only using static linking for
supported platforms.
Linux:       Red Hat Enterprise Linux Workstation release 6.3 (Santiago)

Sun Solaris: Oracle Solaris 10 9/10 s10s_u9wos_14a SPARC
             SunOS frigg 5.10 Generic_147440-27 sun4u sparc
             SUNW,Sun-Blade-2500

Mac OS X:    Snow Leopard, v10.6.8, Darwin Kernel Version 10.8.0

Windows 7:   Enterprise, Service Pack 1
---

PS: README.TXT *********************

E-PAPS ARTICLE INFORMATION:

Author(s):     John Boone & J. Anthony Seibert

Article Title: An accurate method for computer-generating tungsten anode...

Journal:       Medial Physics   E-MPHYA-24-1661
Issue Date:    November 1997    Volume: 24  Issue No.: 11

DEPOSIT INFORMATION

Description: General DX Spectral Data

Total No. of Files: 15
File Names:See Below
Filetypes: 
Special Instructions:
This EPAPS data repository holds the following files:

genspec1.c
     the subroutine which generates the TASMIP spectra
TASMIP algorithm
test1.c
     a main calling routine which queries for input
     paramters (kV, added Al, and ripple), and returns the
     spectrum into a file SPECTRUM.DAT
mual.h
     attenuation coefficients for aluminum
genspec1.h
     the polynomial coefficients that are the key to the
spec100.dat
     A sample spectrum generated using TEST1.c and
GENSPEC1.c for comparison

TABULAR DATA FILES:

The files listed below, HVLS.*, are tabular data.  The
extension to the file name corresponds to the thickness of
water, in centimeters, that the x-ray spectrum has been
filtered by.  In each file, there are 5 columns:
KV   RIPPLE( in %)   AL(in mm)  Fluence(photons per square
millimeter per mR)  HVL(mm AL)

hvls.0
hvls.5
hvls.10
hvls.15
hvls.20
hvls.25
hvls.30
hvls.35
hvls.40

Author Contact Information:
John M. Boone, Dept. of Radiology, UC Davis Medical Center,
FOLB II E, 2421 45th Street, Sacramento, California
Ph: 916-734-5059  Fax: 916-734-3111  E-mail: jmboone@ucdavis.edu

******************* E-PAPS: README.TXT *********************
