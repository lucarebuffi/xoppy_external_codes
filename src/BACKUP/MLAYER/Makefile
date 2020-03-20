# --------------------------------------------------------------
#
# Works automatically for Sun,hp and Linux systems
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
ARCH := $(shell ../../xop_uname)
#********** PLATFORM DEPENDENT VARIABLES *********************
#
# libraries for sun5 
#
ifeq ($(ARCH),sun5)
FC = g95
FFLAGS = 
endif
#
# libraries for linux 
#
ifeq ($(ARCH),linux)
FC = g95
FFLAGS = -Dlinux
endif
#
# libraries for SGI 
#
ifeq ($(ARCH),irix)
FC = f77
FFLAGS =  -DSGI 
endif
#
# libraries for darwin
#
ifeq ($(ARCH),darwin)
FC = g95
FFLAGS = -Dlinux
endif
#**********************************************************************


all: mlayer 

mlayer:  mlayer.F
	$(FC) $(FFLAGS) -o mlayer mlayer.F \
	io_shadow.f iyes.f vector.f iblank.f
#	io_shadow.f iyes.f vector.f READLIB.F iblank.f

clean: 
	/bin/rm -f *.trace *.o mlayer 

install:
	cp mlayer ./../../bin.$(ARCH)/mlayer
