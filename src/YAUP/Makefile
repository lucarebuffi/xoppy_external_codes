# --------------------------------------------------------------
#
# gnu makefile for 
#
# MSR Jan 1997
#
# Works automatically for Sun4Sun5,hp700 and Linux systems
#
# Usage:
#	gmake            [to create the executable crystal_ml]
#	gmake install	 [to install  it in the correct directory]
#	gmake clean	 [to clean the source directory]
#
#
# --------------------------------------------------------------
#
# System definition (architecture)
#
ARCH := $(shell ./../../xop_uname)
#*******************************************************
# libraries for sun5 
#
ifeq ($(ARCH),sun5)
FC = f77
FFLAGS =  -Bstatic
endif
#*******************************************************
# libraries for linux 
#
ifeq ($(ARCH),linux)
FC = gfortran
FFLAGS =  -O
endif
#
#-------------------------- libraries for darwin
#
ifeq ($(ARCH),darwin)
FC = g95
FFLAGS =  
endif




all: 
	$(FC) $(FFLAGS) -o bfield BFIELD.F
	$(FC) $(FFLAGS) -o bfield2 BFIELD2.F
	$(FC) $(FFLAGS) -o txt2u TXT2U.F
	$(FC) $(FFLAGS) -o u2txt U2TXT.F
	$(FC) $(FFLAGS) -o yaup YAUP.F

clean: 
	/bin/rm -f *.trace *.o bfield bfield2 txt2u u2txt yaup

install:
	cp bfield ./../../bin.$(ARCH)/bfield
	cp bfield2 ./../../bin.$(ARCH)/bfield2
	cp txt2u ./../../bin.$(ARCH)/txt2u
	cp u2txt ./../../bin.$(ARCH)/u2txt
	cp yaup ./../../bin.$(ARCH)/yaup
