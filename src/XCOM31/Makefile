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

#************ SYSTEM DEPENDENT VARIABLES ******************
#
# libraries for sun5 
#
ifeq ($(ARCH),sun5)
FC = g77
FLAGS =  
endif
#
# libraries for Linux 
#
ifeq ($(ARCH),linux)
FC = g95
FLAGS =  -fsloppy-char 
endif
#
#-------------------------- libraries for darwin
#
ifeq ($(ARCH),darwin)
FC = g95
FLAGS =  -fsloppy-char 
endif
#************ END OF SYSTEM DEPENDENT VARIABLES ******************

all: 
	$(FC) $(FLAGS) -o xcom xcom.f

clean: 
	/bin/rm -f *.trace *.o xcom

install:
	cp xcom ./../../bin.$(ARCH)/xcom
