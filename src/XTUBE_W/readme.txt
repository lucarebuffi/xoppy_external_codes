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
