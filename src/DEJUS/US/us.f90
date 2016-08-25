	PROGRAM US
!+
! PROGRAM DESCRIPTION:	
!  Program to calculate undulator spectra using the Bessel function 
!  approximation for an ideal planar undulator or an ideal elliptical
!  undulator (including polarization for both cases).  
!  The program may be executed from the XOP interface.
! 
! AUTHORS: 
!  Roger J. Dejus
!  The Advanced Photon Source
!  Experimental Facilities Division (now in Accelerator Systems Division)
!  Argonne National Laboratory
! 
! CREATION DATE: 
!  25-MAR-1991
! 
! INPUT PARAMETERS:
!  The input parameters are divided into sections related to the storage ring,
!  the undulator, and the quantity to be calculated. Note: when modifying
!  parameters under the Xus interface, double click the input field parameter and
!  press the RETURN key so that the new parameter is accepted.
! Title:				TITLE
! Machine Parameters:
!  Storage ring energy 			ENERGY (GeV)
!  Storage ring current			CUR (mA)
!  RMS relative beam energy spread	SIGE
!  RMS beam size (horizontal)		SIGX (mm)
!  RMS beam size (vertical)		SIGY (mm)
!  RMS beam divergence (horizontal)	SIGX1 (mrad)
!  RMS beam divergence (vertical)	SIGY1 (mrad)
! Undulator Parameters:
!  Period length			PERIOD (cm)
!  Number of periods			NPER
!  Deflection parameter (hor.  field) 	KX (= 0.0 for a regular planar device)
!  Deflection parameter (vert. field) 	KY
! Scan Parameters:
!  Minimum energy			EMINU (eV)
!  Maximum energy EMAX			EMAXU (eV)
!  Number of energy points		NEU
! Pinhole Parameters:
!  Distance from the source		DU (m)
!    (DU=0.0 => angular units)
!  X-coordinate for center of pinhole	XPC (mm) or (mrad)
!  Y-coordinate for center of pinhole	YPC (mm) or (mrad)
!  X-size of pinhole (full width)	XPS (mm) or (mrad)
!  Y-size of pinhole (full width)	YPS (mm) or (mrad)
!    (for angular units (DU=0.0) values are entered in mrad)
!    (X is for the horizontal direction)
!    (Y is for the vertical direction)
!  Number of X subdivisions of pinhole	NXP
!  Number of Y subdivisions of pinhole	NYP
!    (for plotting 3D results with interface Xus, the X-size, Y-size, and the number of
!     of subdivisions in the two directions should be equal)
!
! Mode:
!  Depending on the mode selected, some of the pinhole parameters may be set to different values
!  by the program; see below and check the beginning (header) in the output file.
!  MODE    1    Angular/spatial flux density distribution
!  MODE    2    Angular/spatial flux density spectrum
!  MODE    3    On-axis brilliance spectrum
!  MODE    4    Flux spectrum through a pinhole
!  MODE    5    Flux spectrum integrated over all angles
!  MODE    6    Power density and integrated power
!
!  Angular/spatial flux density distribution
!    - Flux distribution at the chosen minimum energy EMINU.
!      The EMAXU and NEU are not used.
!  Angular/spatial flux density spectrum
!    - Spectrum at a given point in space selected by the XPC and YPC coordinate
!      for the center of the pinhole.
!      The XPS, YPS, NXP and NYP are not used.
!  On-axis brilliance spectrum
!    - The pinhole parameters have no significance here.
!      The DU, XPC, YPC, XPS, YPS, NXP and NYP are not used.
!  Flux spectrum through a pinhole
!    - Spectrum through a pinhole centered at XPC and YPC with total size XPS and YPS.
!  Flux spectrum integrated over all angles.
!    - The pinhole parameters have no significance here.
!      The DU, XPC, YPC, XPS, YPS, NXP and NYP are not used.
!  Power density and integrated power
!    - Integrated over all energies, thus the energy parameters have no
!      significance here. The EMINU, EMAXU and NEU are not used.
!      The SIGE is not used and is set to zero.
!
! Method (N refers to the number of undulator periods):
!  METHOD  1    Non-zero emittance; finite-N
!  METHOD  2    Non-zero emittance; infinite-N
!  METHOD  3    Zero emittance;     finite-N
!  METHOD  4    Non-zero emittance; infinite-N + convolution (Dejus' approach)
!  METHOD 14    Non-zero emittance; infinite-N + convolution (Walker's approach)
!
!  Non-zero emittance; finite-N
!    - Required for MODE 1 "Angular/spatial flux density distribution."
!      Use also for MODE 6 "Power density and integrated power" (any non-zero emittance method may be used).
!      Not allowed for angle-intergrated spectra (MODE 5). Not recommended for other modes due to slow speed.
!  Non-zero emittance; infinite-N
!    - For test purposes; do not use (not available from the XOP menu).
!  Zero emittance; finite-N
!    - Use for zero emittance calculations. Not allowed for angle-integrated spectra (MODE 5).
!  Non-zero emittance; infinite-N/convolution
!    - Recommended for all runs with emittance.
!    - Only METHOD 4 (Dejus' approach) is available on the XOP menu. This method uses an internally
!      generated energy mesh with variable step size: at the location of the harmonics the step size
!      is made small and in between the harmonics the step size is made large.
!
! Harmonic Number:
!  IHARM   0    All harmonics
!  IHARM  -1    Lowest order harmonic (except MODE=6, include to harmonics -IHARM)
!  IHARM   I    I''th harmonic
!
!  All harmonics
!    - Selects all contributing harmonics (generally used).
!  Lowest order harmonic
!    - Selects the lowest order contributing harmonic.
!  Harmonic #
!    - Selects the harmonic number displayed.
!  Edit harmonic number (XOP Menu)
!    - Modifies the displayed harmonic number.
!    - For MODES 1,2,3,4 with METHOD 3 (zero emittance) the IHARM -1 is not used.
!   
! Intrinsic Parameters:
!  Several internal parameters used in the calculations. They are commonly not modified
!  by the user. All parameters can be set to zero in which case they default to the values
!  given in the parenthesis.
!
!  NPHI    - Number of steps in angle phi between 0 and pi/2.0 (20).
!            Used in MODES 1,2,3,4,5 for non-zero emittance calculations.
!
!  NALPHA  - Number of steps in angle alpha (gamma*theta) (40).
!            Used in MODES 1,2,3,4 for METHOD 1 (non-zero emittance with finite-N).
!
!  CALPHA2 - Range of angles in alpha^2 in units of the angular equivalent of 1/N (2.0). 
!            Used in MODES 1,2,3,4 for METHOD 1 (finite-N) and for METHOD 3 (zero emittance calculations).
!
!  NOMEGA  - Number of steps in photon energy for the natural lineshape (64).
!            Used in MODES 2,3,4,5 for METHOD 14 (infinite-N + convolution Walker's method).
!
!  COMEGA  - Range of photon energies to be included in the natural lineshape in units (energy of fundamental/N) (8.0)
!            The default value covers the range +/- 2/N of the natural lineshape.
!            Used in MODES 2,3,4,5 for METHOD 14 (infinite-N + convolution Walker's method).
!
!  NSIGMA  - Number of standard deviations of the electron beam size and divergence (4).
!            Used in MODES 1,2,3,4,6 for non-zero emittance calculations.
!
! Polarization:
!  The normalized Stokes parameters are calculated including the 
!  unpolarized component.
!  
!  NOTES:
!
!  1) For MODES 2,3,4,5 the finite-N spectrum is obtained by convoluting the infinite-N
!     spectrum with the natural lineshape. For Walker's method the point spacing in photon
!     energy must be the same for the two curves. This can be achieved as follows: set NEU=0,
!     in which case the spacing is set by the values of NOMEGA and COMEGA and NEU is set accordingly.
!     Set NEU to the approximate number of points desired in the energy range EMINU, EMAXU.
!     A new value of NEU is then calculated which gives the closest match with the spacing
!     of the natural lineshape. In either case EMAXU will also be adjusted so that the convolution
!     can be carried out correctly over the defined energy region.
!
!  2) If DU is set to zero, this indicates that angular flux and power density is to be calculated
!     rather than spatial flux and power density in MODEs 1,2,4 and 6. In this case SIGX and SIGY
!     are ignored, and the acceptance XPC, YPC, XPS, YPS is entered in mrad rather than mm units.
!
!  3) If the acceptance is centred on the axis (XPC=YPC=0.0) then only one quarter of the acceptance
!     needs to be calculated because of symmetry. In this case the range from (0,0) to (XPS/2.0,YPS/2.0)
!     will be divided into NXP,NYP intervals. The printed values of integrated flux and power, including
!     Stokes parameters will however be correct for the total acceptance.
!
!  4) The angle theta (alpha/gamma) is the angle between the undulator axis and the direction of observation.
!     The angle phi is the angle between the projection of the angle of observation in the x-y plane and the x-axis.
!
!  5) For MODE 6 "Power density and integrated power" with non-zero emittance an internally generated
!     cartesian grid is used and none of the intrinsic parameters are used except NSIGMA.
!
!  6) The definition of SIGX must include the contribution of the horizontal dispersion times the beam energy spread.
!     (The dispersion is not an input parameter and hence the user must enter the correct value of SIGX.)
!
!  7) The variable names of some parameters are changed when printed. For example DU is printed as D (distance).
!     EMINU is printed as EMIN, etc. The trailing "U" of a variable name indicates a user value. Some of those
!     are changed inside the code and the actual values used are printed.
! 
! DESIGN ISSUES:
!  Program is based on the Bessel function approximation and is valid in the
!  far-field for an ideal sinusoidal magnetic field profile. It is further 
!  based on the code URGENT by Richard P. Walker with added features.
!  
! COPYRIGHT:
!  This routine may be used at The Advanced Photon Source and any other facility 
!  without explicit consent of the author. No warranties are given as to the
!  accuracy of the results.
!  
! FILES USED:
!  Input file - us.dat (us.inp for XOP)  File in the user's current directory
!  containing the input parameters.
!  Output file - us.plt (us.out for XOP) File in the user's current directory
!  containing the results of the calculation. The header contains all input parameters
!  and the calculated zero emittance on-axis first harmonic energy (e1), the corresponding
!  wavelength (l1), total power (ptot), and the on-axis power density (pd).
!
! KEYWORDS:
!  Undulator Spectrum, Bessel Function Approximation
!  
! LINK/LIBRARY ISSUES:
!  Calls routines BRIGHTE and HUNT. BRIGHTE calculates the brilliance and HUNT
!  searches an array of real numbers (from Numerical Recipes).
!  
! PORTABILITY ISSUES:
!  Runs on DEC 3000/400 AXP alpha (Tru64Unix v5.0), SUN (Solaris: SunOS
!  Release v5.6), and Windows 95/98/NT (Pentium and higher).
!
!  Updated October 8, 2013 (Argonne National Laboratory)
!  *** Linux Red Hat Enterprise Linux Workstation release 6.3 (Santiago) ***
!  Red Hat Enterprise Linux (RHEL) 64-bit with the Intel(R) Fortran
!  Intel(R) 64 Compiler XE for applications running on Intel(R) 64,
!  Version 13.1.1.163 Build 2013031, and with GFORTRAN, gcc version 4.4.6 20120305
!  (Red Hat 4.4.6-4) (GCC).
!  *** Sun Solaris SunOS 5.10 Generic_147440-27 sun4u sparc SUNW,Sun-Blade-2500 ***
!  Sun Fortran 90/95 8.4 SunOS_sparc Patch 128231-02 2009/10/20 with the -f77 switch.
!  and with GFORTRAN, gcc version 4.5.1 (GCC).
!  Windows 7/8 64-bit and MacOS X 10.6 (and newer) are also supported.
!  The GFORTRAN compiler (GCC) v4.8.1 was used for compilations on Windows and (GCC) v4.6.1 on MacOS.
!
!  Updated November 24, 2014 (Argonne National Laboratory)
!  *** Linux Red Hat Enterprise Linux Workstation release 6.5 (Santiago) ***
!  Red Hat Enterprise Linux (RHEL) 64-bit with the Intel(R) Fortran
!  Intel(R) 64 Compiler XE for applications running on Intel(R) 64,
!  Version 14.0.1 Build 20131008
!  GNU Fortran (GCC) 4.4.7 20120313 (Red Hat 4.4.7-4)
!  Copyright (C) 2010 Free Software Foundation, Inc.
!
!  *** Sun Solaris SunOS 5.10 Generic_147440-27 sun4u sparc SUNW,Sun-Blade-2500 ***
!  Sun Fortran 90/95 8.4 SunOS_sparc Patch 128231-02 2009/10/20 with the -f77 switch.
!  GNU Fortran (GCC) 4.5.1
!  Copyright (C) 2010 Free Software Foundation, Inc.
!
!  *** Windows 7/8 64-bit ***
!  GNU Fortran (GCC) 4.9.1
!  Copyright (C) 2014 Free Software Foundation, Inc.
!
!  *** MacOS X 10.6 - 10.10 ***
!  GNU Fortran (GCC) 4.9.2 20141029 (prerelease)
!  Copyright (C) 2014 Free Software Foundation, Inc.
!
! TIMING:
!  Execution times vary considerably depending on computer and the 
!  quantity being calculated. The zero emittance calculations are fast
!  (few seconds), whereas the non-zero emittance calculations may range from
!  seconds (on-axis brilliance) to an hour (flux spectrum through a pinhole).
!
! EXAMPLES:
! Ex. 1 using the input file ~/test/us.txt (output file becomes us.plt in the current working directory)
! % us ~/test/us.txt
! Ex. 2 using the default input file us.dat in the current working directory (the output file becomes us.plt).
! % us
! Ex. 3 using the input abc in the current working directory (the output file becomes abc.plt).
! % us abc 
!
! VERSION:
!  1.94
!  
! MODIFICATION HISTORY:
! 
!	 Date     | Name  | Description
! ----------------+-------+-----------------------------------------------------
! 06-JUL-1994     | RJD   | Modified value for E1MIN for angle-integrated
!                 |       | spectrum (MODE=5) to be non-zero; gamma*theta 
!                 |       | corresponds to sqrt(200) (somewhat arbitrarily
!                 |       | chosen)
! ----------------+-------+-----------------------------------------------------
! 04-OCT-1994     | RJD   | Modified program to include polarization properties.
!                 |       | The four Stokes parameters are now calculated.
!                 |       | Program is for an ideal planar undulator or an ideal
!                 |       | elliptical undulator. Many other changes. The value
!                 |       | of the parameter IHARM has a different meaning.
!                 |       | IHARM=0 now gives 'all harmonics' and IHARM= <0
!                 |       | gives the lowest order harmonic except for the power
!                 |       | option. For the power option, a negative IHARM means
!                 |       | include all harmonics up to and including -IHARM.
!                 |       | This version is 1.6.
! ----------------+-------+-----------------------------------------------------
! 21-JUN-1995     | RJD   | Modified print-out of "Contributing harmonics" in
!		  |       | subroutine PRINT_OUT. Routine incorrectly calculated 
!		  |       | IMIN and IMAX for METHOD 4 (Dejus method) for
!		  |       | "Spectral distributions". The spectra and integrated
!		  |       | quantities were calculated correctly and are 
!		  |       | unaffected by this modification.
!                 |       | The current version is 1.7.
! ----------------+-------+-----------------------------------------------------
! 04-JAN-1996     | RJD   | Modified the number of decimal places for the sigx1
!                 |       | and sigy1 variables to four in the printout. Added
!                 |       | one more digit for the emax variable to avoid
!                 |       | overflow on rare occasions. Formats 260 and 256 were
!                 |       | changed accordingly.
!                 |       | The current version is 1.8.
! ----------------+-------+-----------------------------------------------------
! 11-NOV-1997     | RJD   | Changed notation: Brightness -> Brilliance.
!                 |       | The current version is 1.9.
! ----------------+-------+-----------------------------------------------------
! 16-JUL-2000     | RJD   | Minor change in the code to compile error-free on
!                 |       | Unix and Windows (no change in results vs. v1.9).
!                 |       | Current version is v1.91.
! ----------------+-------+-----------------------------------------------------
! 02-NOV-2013     | RJD   | Updated date and time routines.
!                 |       | Changed printout of number of decimal places for sigx and sigy.
!                 |       | Increased the number of subdivisions of the pinhole 
!                 |       | from 50 to 200 (P_SZ=201).
!                 |       | Changed rank of variable "SL" to become an array in subroutine
!                 |       | ANGLE_INTEGRATION to avoid compilation warning with gfortran.
!                 |       | Current version is v1.92.
! ----------------+-------+-----------------------------------------------------
! 10-JUN-2014     | RJD   | Updated so that an arbitrary input file can be used on the command line.
!                 |       | If no input file is given on the command line then the file 'us.dat'
!                 |       | is assumed ('us.inp' for the XOP version). The output filename is created
!                 |       | from the rootname, which is derived from the input filename using the string after
!                 |       | the last directory separator (/) without its trailing file extension (if it exists).
!                 |       | The output filename is the rootname with the extension (.plt) appended (.out for the
!                 |       | XOP version). Search "standalone" for changing defaults.
!                 |       | Current version is v1.93.
! ----------------+-------+-----------------------------------------------------
! 22-OCT-2014     | RJD   | Added beam energy spread to all modes and introduced parameter SIGE.
!                 |       | Uses routine econ_func() to do the energy convolution.
!                 |       | NOTE: the definition of SIGX remains unchanged and it must include the value of the 
!                 |       | horizontal dispersion times the beam energy spread.
!                 |       | Fixed error in SUBROUTINE CONVOLUTE_ENERGY_VSTEP. The array HE() was
!                 |       | accessed with index 0 (out of bounds for J2 = 1). See RJD 10/23/2014.
!                 |       | Fix has no effect on the results because the calculated spectra SP1 would
!                 |       | typically be zero at the boundary.
!                 |       | Completely rewritten to take advantage of modern Fortran 90/95 features including dynamic
!                 |       | allocation of arrays. Array size limitations were removed. Warnings are given if predefined array
!                 !       ! sizes are exceeded but code will run ok although execution times may be quite large.
!                 !       ! Rewrote algorithm for power density distributions (MODE 6) for the non-zero emittance case.
!                 !       ! Now uses a cartesian coordinates for the convolution. Previously used polar coordinates and it gave
!                 !       ! spurios/unreal results for small values of the beam divergence.
!                 |       | Increased values of default parameters NALPHA, NSIGMA, NOMEGA, and COMEGA.
!                 |       | Numerous checks added.
!                 |       | Current version is v1.94.
! ----------------+-------+-----------------------------------------------------
! [change_entry]
!-
	use us_size
	use prt
	use calc
	use factor
	use beam
	use pinhole
	use angle_phi
	use angle_alp
	use energy_mod
	use step
	use line_shape
	use spectra
	use resize
	use physical_constants !  Fundamental physical constants; Physics Today Aug. 1990:
	use pi_constants
!	use precision_standard	! not needed explicitly
	implicit none

!  Size parameters:
	integer(kind=i4b),parameter :: BLOCK_SZ=12288,BUFFER_SZ=1 ! Blocksize in bytes
!	integer(kind=i4b),parameter :: RECL_SZ=((F_SZ*20/512+1)*512) ! 512 longwords
  
!  Declarations of scalars:
	character(len=1) :: CHR
	character(len=80) :: TITLE
	character(len=8) :: DBUFF1
	character(len=10) :: DBUFF2
	character(len=10) :: TBUFF1
	character(len=8) :: TBUFF2
	character(len=180) :: EXE_NAME,FILE_NAME,FILE_NAME_OUT,datafile,plotfile,rootname

	integer(kind=i4b) :: IERROR
	integer(kind=i4b) :: J,L
	integer(kind=i4b) :: NCH,NCH1,NCH2,N1,N2
	integer(kind=i4b) :: NOMEGA,NEI
	integer(kind=i4b) :: I,IA,IB,IE,ID,IW,ISIGN_US,ISUB
	integer(kind=i4b) :: IE1,IE2,IDW
	integer(kind=i4b) :: SD,ED,ESIZE

	real(kind=sp) :: DTIME,DELTA
	real(kind=dp) :: ENERGY,CUR,PERIOD
	real(kind=dp) :: SIGX,SIGY,SIGX1,SIGY1
	real(kind=dp) :: EMIN,EMAX,DU
	real(kind=dp) :: COMEGA
	real(kind=dp) :: G2,LAMDAR,LAMDA1,E1Z
	real(kind=dp) :: ARG,D2,SL
	real(kind=dp) :: XE,XE1,XE2,YE,YE1,YE2,DEU
	real(kind=dp) :: XPMIN,YPMIN
	real(kind=dp) :: E1MIN,E1MAX,EP,EF,DEF
	real(kind=dp) :: PHI1,PHI2
	real(kind=dp) :: PTOT,PD,GK
	real(kind=dp) :: K2

!  Declarations of arrays:
	real(kind=sp),dimension(2) :: TD=(/0.0,0.0/)

!  Conversion factors:
	real(kind=dp),parameter :: C_EVANG=H*C/EC*1.0D10,C_MM_M=1.0D-3,C_M_MM=1.0D3
	real(kind=dp),parameter :: C_CM_ANG=1.0D8,C_MRAD_RAD=1.0D-3,C_MA_A=1.0D-3
	real(kind=dp),parameter :: C_CM_M=1.0D-2,C_RAD_MRAD=1.0D+3

!  Labeled constants:
	character(len=*),parameter :: VN='v. 1.94'
	integer(kind=i4b),parameter :: IDWMAX=128
	real(kind=dp),parameter :: FINE_STRUCTURE_CONST=1.0D19*EC*1.0D19*EC/(4.0D0*PI*EPSZ*1.0D38*HBAR*C)
	real(kind=dp),parameter :: PTOT_FAC=PI/3.0D0*EC/EPSZ/(MEE**2)*1.0D+6 ! 0.07257
	real(kind=dp),parameter :: PD_FAC  =21.0D0/(16.0D0*PI*MEE**2)*PTOT_FAC ! 0.11611
	real(kind=dp),parameter :: PDH_FAC =EC/EPSZ/(C_RAD_MRAD**2) ! 1.80951e-14
	real(kind=dp),parameter :: ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,HALF=0.5D0
	real(kind=dp),parameter :: EPS=1.0D-6,EPSE=1.0D-8,EPSK=1.0D-4

!  Filename defaults (standalone and XOP)
!	character(len=*),parameter :: FILE_NAME_IN = 'us.dat'	! default input filename in the current directory
!	character(len=*),parameter :: FT_OUT='.plt'		! default output filename extension
	character(len=*),parameter :: FILE_NAME_IN = 'us.inp'	! default input filename in the current directory (XOP)
	character(len=*),parameter :: FT_OUT='.out'		! default output filename extension (XOP)

!  Logical name for parameter file
	character(len=*),parameter :: LN='USPAR:'
!	character(len=*),parameter :: LN='USPAR'

!
!+ First executable statement here
!  -----------------------------------------------------------------------------
!2345678901234567890123456789012345678901234567890123456789012345678901234567890
!	CALL LIB$INIT_TIMER
!	CALL LIB$SHOW_TIMER
	DELTA = DTIME(TD)       ! start clock
	CALL DATE_AND_TIME(DATE=DBUFF1,TIME=TBUFF1)
	TBUFF2 = TBUFF1(1:2)//':'//TBUFF1(3:4)//':'//TBUFF1(5:6)
	DBUFF2 = DBUFF1(1:4)//'-'//DBUFF1(5:6)//'-'//DBUFF1(7:8)

!  check number of command line arguments (not used here)
!	if (command_argument_count() .ne. 1) then
!	   print *, "Usage: us data_file"
!	   call exit(0)	   
!	endif

!  extract the input filename to be used as the rootname
	CALL GET_COMMAND_ARGUMENT(1,DATAFILE)
!	IF (TRIM(ADJUSTL(DATAFILE(1:1))) .EQ. '') DATAFILE = FILE_NAME_IN	! default for undefined variable (adjustl not necessary
	IF (TRIM(DATAFILE(1:1)) .EQ. '') DATAFILE = FILE_NAME_IN	! default for undefined variable
	SD = SCAN(DATAFILE, '/', BACK=.TRUE.)
	ED = SCAN(DATAFILE, '.', BACK=.TRUE.)
	ROOTNAME = DATAFILE(SD +1:ED -1)
	IF (TRIM(ROOTNAME(1:1)) .EQ. '') ROOTNAME = DATAFILe	! default for undefined variable
	CALL CHECK_FILE(DATAFILE)
	PLOTFILE = TRIM(ROOTNAME)//FT_OUT

!	print *,'datafile:',datafile
!	print *,'rootname:',rootname
!	print *,'plotfile:',plotfile

!+ Open plot file
!  -----------------------------------------------------------------------------
!	OPEN(UNIT=2,FILE='us.plt',STATUS='unknown')
!	PRINT 200,'****************** Start ',EXE_NAME(NCH1:NCH2),
!	PRINT 201
	open(unit=2,file=trim(plotfile),status='unknown')
	PRINT 200,'****************** Start US ',VN,' at ',DBUFF2,' ',TBUFF2,' ******************'
	PRINT *
	WRITE (2,200) '****************** Start US ',VN,' at ',DBUFF2,' ',TBUFF2,' ******************'
	WRITE (2,200)

!+ Read parameters from data file
!  -----------------------------------------------------------------------------
!	FILE_NAME = LN
!	FILE_NAME = LN//FILE_NAME_IN
!	FILE_NAME = FILE_NAME_IN
!	OPEN (UNIT=1,FILE=FILE_NAME,STATUS='OLD',ACCESS='SEQUENTIAL',
!	1     FORM='FORMATTED',ACTION='READ')

 	OPEN (UNIT=1,FILE=trim(datafile),STATUS='OLD',ACCESS='SEQUENTIAL',FORM='FORMATTED',ACTION='READ')

	READ (1,100) TITLE
	READ (1,*)   ENERGY,CUR,SIGE	      ! [GeV],[mA],dE/E (dim. less)
	READ (1,*)   SIGX,SIGY,SIGX1,SIGY1    ! [mm],[mm],[mrad],[mrad]
	READ (1,*)   PERIOD,NPER,KX,KY        ! [cm],[dl],[dl],[dl]
	READ (1,*)   EMINU,EMAXU,NEU	      ! [eV],[eV],[dl]
	READ (1,*)   DU,XPC,YPC,XPS,YPS,NXP,NYP![m],[mm],[mm],[mm],[mm],[dl],[dl]
	READ (1,*)   MODE,METHOD,IHARM
	READ (1,*)   NPHI,NALPHA,CALPHA2,NOMEGA,COMEGA,NSIGMA

	CLOSE (UNIT=1)
!-
!  -----------------------------------------------------------------------------

!+ Print mode selection
!  -----------------------------------------------------------------------------
	IF (MODE .EQ. 1) THEN
	    PRINT 200,'Angular/spatial flux density distribution'
	ELSE IF (MODE .EQ. 2) THEN
	    PRINT 200,'Angular/spatial flux density spectrum'
	ELSE IF (MODE .EQ. 3) THEN
	    PRINT 200,'On-axis brilliance spectrum'
	ELSE IF (MODE .EQ. 4) THEN
	    PRINT 200,'Flux spectrum through a pinhole'
	ELSE IF (MODE .EQ. 5) THEN
	    PRINT 200,'Flux spectrum integrated over all angles'
	ELSE IF (MODE .EQ. 6) THEN
	    PRINT 200,'Power density and integrated power'
	END IF ! MODE
!-
!  -----------------------------------------------------------------------------

!+ Check for valid input
!  -----------------------------------------------------------------------------
	IF (MODE .LT. 1 .OR. MODE .GT. 6) GO TO 912

	IF (METHOD .NE. 1 .AND. METHOD .NE. 2 .AND. METHOD .NE. 3 .AND. &
	    METHOD .NE. 4 .AND. METHOD .NE. 14) GO TO 914

	IF ((MODE .EQ. 1 .AND. METHOD .EQ. 2) .OR.      & ! limit use to finite-N
	   (MODE .EQ. 1 .AND. METHOD .EQ. 4) .OR.       & ! and zero emittance
	   (MODE .EQ. 1 .AND. METHOD .EQ.14)) GO TO 930 

	IF ((MODE .EQ. 5 .AND. METHOD .EQ. 1) .OR.      &
	   (MODE .EQ. 5 .AND. METHOD .EQ. 3)) GO TO 940 ! limit use to infinite-N
                                                         ! + convolution
	IF (XPC .LT. ZERO .OR. YPC .LT. ZERO) GO TO 920
!-
!  -----------------------------------------------------------------------------

!+ Default values
!  -----------------------------------------------------------------------------
	IF (NPHI    .EQ. 0) NPHI   = 20
	IF (NALPHA  .EQ. 0) NALPHA = 40
	IF (NOMEGA  .EQ. 0) NOMEGA = 64
	IF (NSIGMA  .EQ. 0) NSIGMA = 4
	IF (CALPHA2 .EQ. ZERO) CALPHA2 = 2.0D0
	IF (COMEGA  .EQ. ZERO) COMEGA  = 8.0D0
	IF (MODE .EQ. 6 .AND. SIGE .NE. ZERO) SIGE = ZERO ! non-zero SIGE does not effect pdf
	IF (SIGE .LT. ZERO) SIGE = ZERO

	IF (NPHI .GT. PHI_SZ) THEN
	    print '(" ",a,i8,a,i8)','&US-W-BIG, Requested steps for angle phi is big; number of steps ',NPHI, &
	          ' is greater than ',PHI_SZ
	    ! NPHI = PHI_SZ
	END IF
	IF (NALPHA .GT. ALPHA_SZ) THEN
	    print '(" ",a,i8,a,i8)','&US-W-BIG, Requested steps for angle gamma*theta is big; number of steps ',NALPHA, &
	          ' is greater than ',ALPHA_SZ
	    ! NALPHA = ALPHA_SZ
	END IF

	L = 80
	DO WHILE(TITLE(L:L) .EQ. ' ') ! strip blanks
	    L = L-1
	ENDDO ! WHILE
	NCH1 = L

!  Determine units (angular or spatial units)
	IF (MODE .EQ. 3 .OR. DU .EQ. ZERO) THEN
	    LANG = .TRUE.
	    D    = ONE
	ELSE
	    LANG = .FALSE.
	    D    = DU
	END IF !

	IF (MODE .EQ. 2) THEN
	    XPS = ZERO
	    YPS = ZERO
	    NXP = 0
	    NYP = 0
	ELSE IF (MODE .EQ. 3 .OR. MODE .EQ. 5) THEN
	    XPC = ZERO
	    YPC = ZERO
	    XPS = ZERO
	    YPS = ZERO
	    NXP = 0
	    NYP = 0
	END IF ! MODE
!-
!  -----------------------------------------------------------------------------

!+ Definition of constants and calculation of power/power density
!  -----------------------------------------------------------------------------
	GAMMA_US  = ENERGY/MEE*1.0D3
	G2     = GAMMA_US*GAMMA_US
	LAMDAR = PERIOD*C_CM_ANG/(TWO*G2) ! Reduced wavelength [A]
	ER     = C_EVANG/LAMDAR		  ! Reduced energy    [eV]
	K2     = KX*KX +KY*KY
	K3     = ONE+K2/TWO
	LAMDA1 = LAMDAR*K3		  ! First harmonic on axis  [A]
	E1Z    = ER    /K3		  ! First harmonic on axis [eV]
	D2     = D*D			  ! Distance squared     [m**2]
	LDEV   = NPER*PERIOD*C_CM_M	  ! Length of device        [m]
	NPI    = NPER*PI
	GK     = ZERO
	IF (KX .LT. EPSK .OR. KY .LT. EPSK) THEN
	    K  = KX +KY
            GK = K*((K**6)+(24.0D0*(K**4)/7.0D0)+(4.0D0*K*K)+(16.0D0/7.0D0))/((1.0D0+(K*K))**3.5D0)
	END IF 
	IF (ABS(KX-KY) .LT. EPSK) THEN
	    K  = KX
	    GK = 32.0D0/7.0D0*K/((1.0D0+(K*K))**3.0D0)
	END IF
        PTOT   = PTOT_FAC*NPER*K2*  (ENERGY**2)*CUR*C_MA_A/(PERIOD*C_CM_M)![W]
        PD     = PD_FAC  *NPER*K*GK*(ENERGY**4)*CUR*C_MA_A/(PERIOD*C_CM_M)![W/mrad^2]
!-
!  -----------------------------------------------------------------------------

!+ Beam emittance
!  -----------------------------------------------------------------------------
	FU    = ZERO
	FV    = ZERO
	SIGX2 = SIGX*SIGX
	SIGY2 = SIGY*SIGY
	IF (LANG) THEN ! Angular units
	    SIGU2 = SIGX1*SIGX1*C_MRAD_RAD*C_MRAD_RAD
	    SIGV2 = SIGY1*SIGY1*C_MRAD_RAD*C_MRAD_RAD
	ELSE ! Spatial
	    SIGU2 = (SIGX1*SIGX1+SIGX2/D2)*C_MRAD_RAD*C_MRAD_RAD
	    SIGV2 = (SIGY1*SIGY1+SIGY2/D2)*C_MRAD_RAD*C_MRAD_RAD
	END IF ! LANG
	SIGU  = SQRT(SIGU2) ! [rad]
	SIGV  = SQRT(SIGV2) ! [rad]
	IF (SIGU2 .NE. ZERO) FU = 0.5D0/SIGU2
	IF (SIGV2 .NE. ZERO) FV = 0.5D0/SIGV2
!-
!  -----------------------------------------------------------------------------

!+ Determine min and max emission angles and center of pinhole
!  -----------------------------------------------------------------------------
	XE  = (XPC-XPS/TWO)*C_MM_M/D-NSIGMA*SIGU ! Cartesian angle in x-dir.
	YE  = (YPC-YPS/TWO)*C_MM_M/D-NSIGMA*SIGV ! Cartesian angle in y-dir.
	XE1 = XE	    
	YE1 = YE
	IF (XE .LT. ZERO) XE = ZERO
	IF (YE .LT. ZERO) YE = ZERO
	AP2MIN = G2*(XE*XE+YE*YE)

	XE  = (XPC+XPS/TWO)*C_MM_M/D+NSIGMA*SIGU ! Cartesian angle in x-dir.
	YE  = (YPC+YPS/TWO)*C_MM_M/D+NSIGMA*SIGV ! Cartesian angle in y-dir.
	XE2 = XE	    
	YE2 = YE
	AP2MAX = G2*(XE*XE+YE*YE)
	ARGMAX = NSIGMA*NSIGMA/TWO

	XE = XPC*C_MM_M/D			! Cartesian angle in x-dir.
	YE = YPC*C_MM_M/D			! Cartesian angle in y-dir.
	AP2CNT = G2*(XE*XE+YE*YE)
!-
!  -----------------------------------------------------------------------------

!+ Write input parameters and default values
!  -----------------------------------------------------------------------------
	WRITE (2,200) TITLE(1:NCH1)
	WRITE (2,250) 'energy ',ENERGY,' GeV','current',CUR ,' mA','sige',SIGE
	WRITE (2,260) 'sigx  ',SIGX,'  mm','sigy ',SIGY,' mm', &
	             'sigx1 ',SIGX1,' mr','sigy1 ',SIGY1,' mr'
	WRITE (2,255) 'period ',PERIOD,'  cm','N      ',NPER   ,'   ', &
	             'Kx     ',KX    ,'    ','Ky     ',KY
	WRITE (2,256) 'emin   ',EMINU ,'  eV','emax',EMAXU,' eV', &
	             'ne     ',NEU
	IF (LANG) THEN
	    WRITE (2,261) 'd      ',DU    ,'   m','xpc   ',XPC ,' mr', &
	                 'ypc   ',YPC   ,' mr'
	ELSE
	    WRITE (2,261) 'd      ',DU    ,'   m','xpc   ',XPC ,' mm', &
	                 'ypc   ',YPC   ,' mm'
	END IF ! LANG
	IF (LANG) THEN
	    WRITE (2,257) 'xps    ',XPS   ,'  mr','yps   ',YPS ,' mr', &
	                 'nxp    ',NXP   ,'   ','nyp    ',NYP
	ELSE
	    WRITE (2,257) 'xps    ',XPS   ,'  mm','yps   ',YPS ,' mm', &
	                 'nxp    ',NXP   ,'   ','nyp    ',NYP
	END IF ! LANG
	WRITE (2,258) 'mode   ',MODE  ,'    ','method ',METHOD,'   ', &
	             'iharm  ',IHARM,'   '
	WRITE (2,259) 'nphi   ',NPHI  ,'    ','nalpha ',NALPHA,'   ', &
	             'calpha2',CALPHA2
	WRITE (2,262) 'nomega ',NOMEGA,'    ','comega ',COMEGA,'   ', &
	             'nsigma ',NSIGMA
	WRITE (2,200)
	WRITE (2,263) 'e1     ',E1Z   ,'  eV','l1   ',LAMDA1,'  A'
	IF (LANG) THEN
            WRITE (2,264) 'ptot   ',PTOT  ,'   W','pd ',PD    , &
	           '  W/mr^2'
	ELSE
            WRITE (2,264) 'ptot   ',PTOT  ,'   W','pd ',PD/D2 , &
	           '  W/mm^2'
	END IF ! LANG
	WRITE (2,200)
!-
!  -----------------------------------------------------------------------------

!+ Define energy scale and set up array for the line shape function
!  -----------------------------------------------------------------------------
	IF (METHOD .EQ. 4 .AND. MODE .NE. 1 .AND. MODE .NE. 6) THEN
	    ESIZE = E_SZ	! save current size
	    ALLOCATE(E(E_SZ),HE(E_SZ),HA(E_SZ))
	    IF (MODE .EQ. 5) THEN
		E1MIN = E1Z*K3/(K3+200.0d0) ! use large value for gamma*theta
		E1MAX = E1Z
		DEW   = E1Z/NPER
	    ELSE
		E1MIN = E1Z       *K3/(K3+AP2MAX)
		E1MAX = E1Z       *K3/(K3+AP2MIN)
		DEW   = (E1Z/NPER)*K3/(K3+AP2CNT)
	    END IF ! MODE
	    EW    = COMEGA*DEW
	    PE    = PI/DEW

!  Extend "internal" energy scale to make room for the line shape function
	    DW   = TWO*EW/NOMEGA ! Default step size
	    EMIN = EMINU-EW
	    EMAX = EMAXU+EW

!  Adjust EMIN and EMAX
	    I  = 1
	    EP = I*E1MAX
	    DO WHILE (EP .LT. EMIN)
		I  = I+1
		EP = I*E1MAX
	    END DO ! WHILE
	    IE1  = I
	    EMIN = MAX(IE1*E1MIN,EMIN)

	    IF (MODE .NE. 5) THEN
		EP = IE1*E1MIN
		DO WHILE (EP .LT. EMAX)
		    I  = I+1
		    EP = I*E1MIN
		END DO ! WHILE
		IE2  = I-1
		EMAX = MIN(IE2*E1MAX,EMAX)
	    END IF ! MODE
	    
	    IF (EMAX .LE. EMIN) GO TO 910

	    I   = IE1
	    EP  = I*E1MAX
	    IE  = 1
	    E (IE) = EMIN

!  Generate energy scale with variable step size
	    DO WHILE (E(IE) .LT. EMAX)
		IDW = 1
		DEF = DW /IDW
		EF  = DEW/IDW
		DO WHILE (E(IE) .GT. EP-EF)
		    IDW = IDW*2
		    IF (IDW .GT. IDWMAX) THEN
			EF = ZERO
		    ELSE
			DEF = DW /IDW
			EF  = DEW/IDW
		    END IF ! IDW
		END DO ! WHILE

		DO WHILE (E(IE) .LT. (EP-EPSE) .AND. E(IE) .LT. (EMAX-EPSE))
		    IF (E(IE) .GT. EP-EF) THEN
			IDW = IDW*2
			IF (IDW .GT. IDWMAX) THEN
			    EF = ZERO
			ELSE
			    DEF = DW /IDW
			    EF  = DEW/IDW
			END IF ! IDW
		    END IF ! E(IE)

		    IE       = IE+1
!		    IF (IE .GT. E_SZ) GO TO 900
		    IF (IE .GT. ESIZE) THEN ! extend arrays 1
		        ESIZE = ESIZE+E_SZ
			CALL RESIZE_ARRAY(E ,ESIZE)
			CALL RESIZE_ARRAY(HE,ESIZE)
			CALL RESIZE_ARRAY(HA,ESIZE)
		    END IF
		    HE(IE-1) = DEF
		    IF (IE .EQ. 2) THEN
			HA(IE-1) = HE(IE-1)/TWO
		    ELSE
			HA(IE-1) = (HE(IE-2)+HE(IE-1))/TWO ! Average
		    END IF ! IE
		    E(IE) = E(IE-1)+HE(IE-1)
		END DO ! WHILE

!  Redefine the last point within range of harmonic # i
		SL       = MIN(EP,EMAX)
		E (IE)   = SL-EPSE
		HE(IE-1) = E(IE)-E(IE-1)
		HA(IE-1) = (HE(IE-2)+HE(IE-1))/TWO ! Average
		IE       = IE+1
!		IF (IE .GT. E_SZ) GO TO 900
		IF (IE .GT. ESIZE) THEN ! extend arrays 2
		    ESIZE = ESIZE+E_SZ
		    CALL RESIZE_ARRAY(E ,ESIZE)
		    CALL RESIZE_ARRAY(HE,ESIZE)
		    CALL RESIZE_ARRAY(HA,ESIZE)
		END IF
		HE(IE-1) = TWO*EPSE
		HA(IE-1) = (HE(IE-2)+HE(IE-1))/TWO ! Average
		E (IE)   = E(IE-1)+HE(IE-1)

		I   = I+1 ! Next harmonic
		DEF = I*E1MIN-E(IE)
		IF (DEF .GT. ZERO) THEN
		    IE	     = IE+1
!		    IF (IE .GT. E_SZ) GO TO 900
		    IF (IE .GT. ESIZE) THEN ! extend arrays 3
		        ESIZE = ESIZE+E_SZ
			CALL RESIZE_ARRAY(E ,ESIZE)
			CALL RESIZE_ARRAY(HE,ESIZE)
			CALL RESIZE_ARRAY(HA,ESIZE)
		    END IF
		    HE(IE-1) = DEF
		    HA(IE-1) = (HE(IE-2)+HE(IE-1))/TWO ! Average
		    E (IE)   = E(IE-1)+HE(IE-1)
		END IF ! DEF

		EP  = I*E1MAX
		
	    END DO ! WHILE

!  Adjust end points
	    IF (DEF .GT. ZERO) THEN
		NE = IE-1
	    ELSE
		NE = IE
	    END IF ! DEF
	    NE1 = 1
	    NE2 = NE
	    HE(NE) = ZERO
	    HA(NE) = HE(NE-1)/TWO
	
	    IF (NE2 .GT. E_SZ) THEN 
	        print '(" ",a,i8,a,i8)','&US-W-BIG, Internal energy array is big; number of points ',NE2, &
		      ' is greater than ',E_SZ
	    END IF

	    IF (NEU .EQ. 0) NEU = INT(NE2/100.0D0+ONE)*100.0D0 ! Default
!	    DEU = (EMAXU-EMINU)/NEU
	    DE  = (EMAXU-EMINU)/NEU
!	    IF (NEU .GT. E_SZ) GO TO 902
	    IF (NEU .GT. E_SZ) THEN
	        print '(" ",a,i8,a,i8)','&US-W-BIG, Requested energy array is big; number of points ',NEU, &
		      ' is greater than ',E_SZ
	    END IF
	    NEU = NEU+1	! OK with allocation of one extra element

	    IE = MAX(NE,NEU)	! need max here
	    ALLOCATE(I1(IE),I2(IE),SPEC0(IE),SPEC1(IE),SPEC2(IE),SPEC3(IE))
	    ALLOCATE(EU(NEU))

	    DO IE=1,NEU ! Energy scale (user's selection)
!		EU(IE) = EMINU+(IE-1)*DEU
		EU(IE) = EMINU+(IE-1)*DE
	    END DO ! IE

	ELSE IF (METHOD .EQ. 14 .AND. MODE .NE. 1 .AND. MODE .NE. 6) THEN ! Walker's
	    DEW = (E1Z/NPER)*K3/(K3+AP2CNT)
	    EW  = COMEGA*DEW ! [eV]
	    IF (NEU .GT. 0) THEN ! Redefine NOMEGA
		SL     = (EMAXU-EMINU)/NEU
		NOMEGA = TWO*EW/SL+ONE
		NOMEGA = NOMEGA/2*2 ! Make even
		IF (NOMEGA .LT. 16) THEN ! Reset to 16
		    NOMEGA = 16
		    PRINT *
		    PRINT *  ,'&US-W-ISMALL, NOMEGA less than 16'
		    PRINT *  ,'- NOMEGA reset to 16'
		    PRINT *
		END IF ! NOMEGA
	    END IF ! NEU
	    
	    DW   = TWO*COMEGA/NOMEGA
	    DE   = DW*DEW
	    EMIN = EMINU-EW		 ! Extend lower limit of energy scale
	    IF (EMIN .LT. EPS) EMIN = EPS
	    EMAX = EMAXU+EW		 ! Extend upper limit of energy scale
	    NE   = (EMAX-EMIN)/DE+ONE    ! Redefine the # of intervals
	    NE   = NE+1			 ! # of points
	    NE1  = NOMEGA/2+1		 ! Smallest index; odd
	    NE2  = NE+1-NE1		 ! Largest index
	    NEI  = NE2-NE1		 ! # of intervals

!	    IF (NE .GT. E_SZ) GO TO 905
	    IF (NE .GT. E_SZ) THEN
	        print '(" ",a,i8,a,i8)','&US-W-BIG, Requested energy array is big; number of points ',NE, &
		      ' is greater than ',E_SZ
	    END IF
	    ALLOCATE(I1(NE),I2(NE),E(NE),SPEC0(NE),SPEC1(NE),SPEC2(NE),SPEC3(NE))

!  Generate line shape function (used in energy convolution)
	    NW = NOMEGA+1 ! # of points
	    ALLOCATE(HW(NW))	! size not an issue here
	    DO IW=1,NW
		ARG = (-COMEGA+(IW-1)*DW)*PI	    
		IF (ABS(ARG) .GT. EPS) THEN
		    SL = SIN(ARG)/ARG
		    HW(IW) = SL*SL
		ELSE
		    HW(IW) = ONE
		END IF ! ARG
	    END DO ! IW
!  Generate energy scale
	    DO IE=1,NE
		E(IE) = EMIN+(IE-1)*DE ! [eV]
	    END DO ! IE
	ELSE ! Method = 1,2,3 or Mode = 1,6
	    NE  = NEU	    
	    IF (MODE .EQ. 1 .OR. MODE .EQ. 6) NE = 0 ! no energy dependence
	    IF (NE .GT. 0) DE  = (EMAXU-EMINU)/NE
	    NE  = NE+1			 ! # of points
	    NE1 = 1			 ! Smallest index
	    NE2 = NE+1-NE1		 ! Largest index
	    NEI = NE2-NE1		 ! # of intervals
!  Generate energy scale
!	    IF (NE .GT. E_SZ) GO TO 905
	    IF (NE .GT. E_SZ) THEN
	        print '(" ",a,i8,a,i8)','&US-W-BIG, Requested energy array is big; number of points ',NE, &
		      ' is greater than ',E_SZ
	    END IF
	    ALLOCATE(I1(NE),I2(NE),E(NE),SPEC0(NE),SPEC1(NE),SPEC2(NE),SPEC3(NE))
	    DO IE=1,NE
		E(IE) = EMINU+(IE-1)*DE ! [eV]
	    END DO ! IE
	END IF ! METHOD & MODE
!-
!  -----------------------------------------------------------------------------

!+ Pinhole parameters
!  -----------------------------------------------------------------------------
	IF (XPC .EQ. ZERO .AND. YPC .EQ. ZERO) THEN ! Pinhole centered
	    FAC   = 4.0D0
	    XPMIN = ZERO
	    YPMIN = ZERO
	    IF (NXP .GT. 0) DXP = XPS/TWO/NXP
	    IF (NYP .GT. 0) DYP = YPS/TWO/NYP
	ELSE
	    FAC   = ONE
	    XPMIN = XPC-XPS/TWO
	    YPMIN = YPC-YPS/TWO
	    IF (NXP .GT. 0) DXP = XPS/NXP
	    IF (NYP .GT. 0) DYP = YPS/NYP
	END IF
	NXP = NXP+1
	NYP = NYP+1
	IF (NXP .GT. P_SZ) THEN
	    print '(" ",a,i8,a,i8)','&US-W-BIG, Requested subdivisions for horizontal aperture is big; number of points ',NXP, &
	          ' is greater than ',P_SZ
	    ! NXP = P_SZ
	END IF
	IF (NYP .GT. P_SZ) THEN
	    print '(" ",a,i8,a,i8)','&US-W-BIG, Requested subdivisions for vertical   aperture is big; number of points ',NYP, &
	          ' is greater than ',P_SZ
	    ! NYP = P_SZ
	END IF
!-
!  -----------------------------------------------------------------------------

!+ Set up positions within pinhole and associated cartesian angles
!  -----------------------------------------------------------------------------
	ALLOCATE(XP(NXP),YP(NYP),CX(NXP),CY(NYP))
	ALLOCATE(RA0(NXP,NYP),RA1(NXP,NYP),RA2(NXP,NYP),RA3(NXP,NYP))
	DO IA=1,NXP
	    XP(IA) = XPMIN+(IA-1)*DXP ! Position [mm]
	    CX(IA) = XP(IA)*C_MM_M/D  ! Angle   [rad]
	END DO ! IA
	DO IB=1,NYP
	    YP(IB) = YPMIN+(IB-1)*DYP ! Position [mm]
	    CY(IB) = YP(IB)*C_MM_M/D  ! Angle   [rad]
	END DO ! IB
!-
!  -----------------------------------------------------------------------------

!+ Set up arrays of cos(phi) and sin(phi), and index_phi
!  -----------------------------------------------------------------------------
	IF (METHOD .NE. 3 .AND. MODE .NE. 6) THEN ! non-zero emittance and not pdf
	    IF (XE1 .LE. ZERO .AND. YE1 .LE. ZERO) THEN
  		IF ((MODE .EQ. 2 .AND. XPC .EQ. ZERO .AND. YPC .EQ. ZERO) .OR. MODE .EQ. 3) THEN
		    NPHI4 = NPHI   ! Use only 0 ... 90 deg.
		ELSE
		    NPHI4 = 4*NPHI ! Use full range: 0 ... 360 deg.
		END IF ! MODE
		NPHI_BRIGHT = NPHI! # of indices for the brilliance array (use sym.)
		DPHI   = PIHALF/NPHI
		ISIGN_US  = +1
		SL     = ISIGN_US*DPHI/TWO
	        ALLOCATE(INDEX_PHI(NPHI4),COSPHI(NPHI4),SINPHI(NPHI4),S2SIGN(NPHI4))
		INDEX_PHI(1) = 1
		S2SIGN(1) = +ONE
		COSPHI(1) = COS(SL)
		SINPHI(1) = SIN(SL)
		DO ID=2,NPHI4
		    ARG           = SL+(ID-1)*DPHI
		    INDEX_PHI(ID) = INDEX_PHI(ID-1)+ISIGN_US
		    COSPHI(ID)    = COS(ARG)
		    SINPHI(ID)    = SIN(ARG)
		    IF (ID .EQ. NPHI) THEN
			ISIGN_US = 0
		    ELSE IF (ID .EQ. NPHI+1) THEN
			ISIGN_US = -1
		    ELSE IF (ID .EQ. 2*NPHI) THEN
			ISIGN_US = 0
		    ELSE IF (ID .EQ. 2*NPHI+1) THEN
			ISIGN_US = +1
		    ELSE IF (ID .EQ. 3*NPHI) THEN
			ISIGN_US = 0
		    ELSE IF (ID .EQ. 3*NPHI+1) THEN
			ISIGN_US = -1
		    END IF ! ID
		    IF (ISIGN_US .NE. 0) THEN
			S2SIGN(ID) = ISIGN_US
		    ELSE ! ISIGN_US = 0
			S2SIGN(ID) = +S2SIGN(ID-1)
		    END IF
		END DO ! ID
	    ELSE 
		IF (XE1 .GT. ZERO) THEN
		    IF (YE1 .GT. ZERO) THEN
			PHI1 = ATAN(YE1/XE2)
			IF (XE1 .GT. EPS) THEN
			    PHI2 = ATAN(YE2/XE1)
			ELSE
			    PHI2 = PIHALF
			END IF
			DPHI  = (PHI2-PHI1)/NPHI
			NPHI4 = NPHI
			NPHI_BRIGHT = NPHI
			ISIGN_US = +1
			SL    = PHI1+ISIGN_US*DPHI/TWO
		    ELSE ! YE1 le 0.0
			IF (XE1 .GT. EPS) THEN
			    PHI1 = ATAN(YE1/XE1)
			    PHI2 = ATAN(YE2/XE1)
			ELSE
			    PHI1 = -PIHALF
			    PHI2 = +PIHALF
			END IF
			DPHI  =  PHI2/NPHI
			N1    = -PHI1/DPHI+ONE
			N2    =  NPHI+1
			NPHI4 = N1+N2
			NPHI_BRIGHT = N2 ! N2 >= N1
			ISIGN_US = +1
			SL    = ISIGN_US*DPHI/TWO
		    END IF ! YE1 gt 0
		ELSE ! XE1 le 0.0
		    PHI1 = ATAN(YE1/XE2)
		    IF (XE1 .LT.- EPS) THEN
			PHI2 = ATAN(YE1/XE1)+PI
		    ELSE
			PHI2 = PIHALF
		    END IF
		    DPHI  = (PIHALF-PHI1)/NPHI
		    N1    = NPHI+1
		    N2    = (PHI2-PIHALF)/DPHI+ONE
		    NPHI4 = N1+N2
		    NPHI_BRIGHT = N1 ! N1 >= N2
		    ISIGN_US = -1
		    SL    = PIHALF+ISIGN_US*DPHI/TWO
		END IF ! XE1 gt 0
	        ALLOCATE(INDEX_PHI(NPHI4),COSPHI(NPHI4),SINPHI(NPHI4),S2SIGN(NPHI4))
		DO ID =1,NPHI4
		    IF (ID .EQ. NPHI_BRIGHT+1) THEN
			ISIGN_US = -ISIGN_US
			SL    = ARG
		    END IF ! ID
		    ARG        = SL+ISIGN_US*(ID-1)*DPHI
		    COSPHI(ID) = COS(ARG)
		    SINPHI(ID) = SIN(ARG)
		    IF (ID .LE. NPHI_BRIGHT) THEN
			INDEX_PHI(ID) = ID
			S2SIGN(ID)    = +ONE
		    ELSE
			INDEX_PHI(ID) = ID-NPHI_BRIGHT
			S2SIGN(ID)    = -ONE
		    END IF ! ID
		END DO ! ID
	    END IF
	END IF ! METHOD
!-
!  -----------------------------------------------------------------------------

! Allocate storage for Brightness arrays
	IF (METHOD .NE. 3 .AND. MODE .NE. 6) THEN ! non-zero emittance and not pdf
	  ALLOCATE(BR0(NPHI_BRIGHT,NALPHA),BR1(NPHI_BRIGHT,NALPHA),BR2(NPHI_BRIGHT,NALPHA),BR3(NPHI_BRIGHT,NALPHA))
        END IF

!+ Scale factors
!  -----------------------------------------------------------------------------
	IF (SIGU .NE. ZERO .AND. SIGV .NE. ZERO) &
	   C1 =NPER*NPER*FINE_STRUCTURE_CONST*BW*CUR*C_MA_A/EC &
	        /(TWOPI*SIGU*SIGV*D2*C_M_MM*C_M_MM)
	C2 = NPER*NPER*   FINE_STRUCTURE_CONST*BW*CUR*C_MA_A/EC
	C3 = NPER*NPER*G2*FINE_STRUCTURE_CONST*BW*CUR*C_MA_A/EC &
	    /(D2*C_M_MM*C_M_MM)
	C4 = PDH_FAC*NPER*NPER*(G2**2)*CUR*C_MA_A/(LDEV*D2)
	IF (SIGU .NE. ZERO .AND. SIGV .NE. ZERO) &
	   C5 = PDH_FAC*NPER*NPER*(G2**2)*CUR*C_MA_A/(TWOPI*SIGU*SIGV*LDEV*D2)
!-
!  -----------------------------------------------------------------------------

!+ Call analysis routine
!  -----------------------------------------------------------------------------
	IF (MODE .EQ. 1) THEN
	    ISUB = 1
	ELSE IF (MODE .GE. 2 .AND. MODE .LE. 4) THEN
	    ISUB = 2
	ELSE IF (MODE .EQ. 5) THEN
	    ISUB = 3
	ELSE IF (MODE .EQ. 6) THEN
	    ISUB = 4
	END IF ! MODE
	IF (METHOD .EQ. 3 .AND. MODE .NE. 6) ISUB = 5

	IF (ISUB .EQ. 1) THEN
	    CALL SPATIAL_DISTRIBUTION(IERROR)
	ELSE IF (ISUB .EQ. 2) THEN
	    CALL SPECTRAL_DISTRIBUTION(IERROR)
	ELSE IF (ISUB .EQ. 3) THEN
	    CALL ANGLE_INTEGRATION(IERROR)
	ELSE IF (ISUB .EQ. 4) THEN
	    CALL POWER_DISTRIBUTION(IERROR)
	ELSE IF (ISUB .EQ. 5) THEN
	    CALL NO_EMITTANCE(IERROR)
	END IF ! ISUB
!-
!  -----------------------------------------------------------------------------

!+ Print results
!  -----------------------------------------------------------------------------
	IF (IERROR .EQ. 0) CALL PRINT_OUT(ISUB)
!-
!  -----------------------------------------------------------------------------

!+ Exit
!  -----------------------------------------------------------------------------
	IF (IERROR .EQ. 0) THEN
	    DELTA = DTIME(TD) ! read clock
	    PRINT 210,'Elapsed time: ',delta,' s'
	    PRINT 200,'&US-S-NORMAL, Successful completion'
	    PRINT *
	ELSE
	    PRINT 200,'&US-F-SUBERR, Subroutine error'
	    PRINT 200,'- unsuccessful completion due to error status'
	    PRINT *
	    WRITE (2,*)
	    WRITE (2,200) '&US-F-SUBERR, Subroutine error'
	    WRITE (2,200) '- unsuccessful completion due to error status'
	    WRITE (2,*)
	    CLOSE (UNIT=2)
	END IF ! IERROR	    

	DEALLOCATE(E,I1,I2,SPEC0,SPEC1,SPEC2,SPEC3)
	DEALLOCATE(XP,YP,CX,CY)
	DEALLOCATE(RA0,RA1,RA2,RA3)
	IF (ALLOCATED(BR0)) 	  DEALLOCATE(BR0,BR1,BR2,BR3)
	IF (ALLOCATED(EU)) 	  DEALLOCATE(EU)
	IF (ALLOCATED(HE)) 	  DEALLOCATE(HE)
	IF (ALLOCATED(HA)) 	  DEALLOCATE(HA)
	IF (ALLOCATED(HW)) 	  DEALLOCATE(HW)
	IF (ALLOCATED(INDEX_PHI)) DEALLOCATE(INDEX_PHI,COSPHI,SINPHI,S2SIGN)

	CALL EXIT(0)

100	FORMAT(A)
110	FORMAT(8I10)
120	FORMAT(8F10.0)
130	FORMAT(8L10)
140	FORMAT(2F10.0,I10)
150	FORMAT(5F10.0,2I10)
160	FORMAT(I10,2(I10,F10.0),I10)

200	FORMAT(' ',8A)
201	FORMAT('1')
205	FORMAT(' ',A,I6,A,I6)
207	FORMAT(' ',A,I2,A)
210	FORMAT(' ',A,2(F10.3,A))

!204	FORMAT(' ',A,A,I6)
!260	FORMAT(' ',2I5,2F10.5)
!270	FORMAT(' ',8F10.5)

250     format(' ',(A,F7.3,A,TR3),(A,F6.1,A,TR3),(A,F12.7))
255	FORMAT(' ',(A,F7.3,A,TR3),(A,I6,A,TR3),2(A,F6.3,A,TR2))
256	FORMAT(' ',(A,F7.1,A,TR3),(A,F9.1,A,TR3),(A,I6))
257	FORMAT(' ',(A,F7.3,A,TR3),(A,F7.3,A,TR3),2(A,I6,A,TR3))
258	FORMAT(' ',(A,I7,A,TR3),(A,I6,A,TR3),2(A,I6,A,TR3))
259	FORMAT(' ',(A,I7,A,TR3),(A,I6,A,TR3),(A,F6.1))
260	FORMAT(' ',(A,f8.4,A,TR3),(A,f8.4,A,TR3),2(A,F7.4,A,TR3))
261	FORMAT(' ',3(A,F7.3,A,TR3))
262	FORMAT(' ',(A,I7,A,TR3),(A,F6.1,A,TR3),(A,I6))
263	FORMAT(' ',(A,F7.1,A,TR3),(A,F8.3,A,TR3))
264	FORMAT(' ',(A,F7.1,A,TR3),(A,F10.1,A,TR3))

!+ Error returns
!  -----------------------------------------------------------------------------
900	CONTINUE
	PRINT 200,'&US-F-BNDERR, Boundary error'
	PRINT 205,'- non-equidistant energy array out of bounds; number of points',IE,' is greater than ',E_SZ
	GO TO 999

902	CONTINUE
	PRINT 200,'&US-F-BNDERR, Boundary error'
	PRINT 205,'- energy array out of bounds; number of points ',NEU,' is greater than ',E_SZ
	GO TO 999

905	CONTINUE
	PRINT 200,'&US-F-BNDERR, Boundary error'
	PRINT 205,'- energy array out of bounds; number of points ',NE,' is greater than ',E_SZ
	GO TO 999

910	CONTINUE
	PRINT 200,'&US-E-HARMERR, Harmonic errror'
	PRINT 210,'- no harmonics reachable in the range ',EMINU,' to ',EMAXU,' eV.'
	GO TO 999

912	CONTINUE
	PRINT 200,'&US-E-INVDAT, Invalid data'
	PRINT 207,'- check input data file; mode ',MODE,' not a valid entry.'
	GO TO 999

914	CONTINUE
	PRINT 200,'&US-E-INVDAT, Invalid data'
	PRINT 207,'- check input data file; method ',METHOD,' not a valid entry.'
	GO TO 999

920	CONTINUE
	PRINT 200,'&US-E-INVDAT, Invalid data'
	PRINT 200,'- check input data file; center of pinhole must lie ','in the first quadrant.'
	GO TO 999

930	CONTINUE
	PRINT 200,'&US-E-INVDAT, Invalid data'
	PRINT 207,'- check input data file; method ',METHOD,' not valid for the flux density distribution.'
	GO TO 999

940	CONTINUE
	PRINT 200,'&US-E-INVDAT, Invalid data'
	PRINT 207,'- check input data file; method ',METHOD,' not valid for angle-integrated spectrum.'
	GO TO 999

999	CONTINUE
	PRINT *
	PRINT 200,'&US-F-PRGERR, Program error'
	PRINT 200,'- unsuccessful completion due to error status'
	PRINT *
	WRITE (2,*)
	WRITE (2,200) '&US-F-PRGERR, Program error'
	WRITE (2,200) '- unsuccessful completion due to error status'
	WRITE (2,*)
!-
!  -----------------------------------------------------------------------------

	CLOSE (UNIT=2)
	CALL EXIT(0)
	END ! US
!
	SUBROUTINE SPATIAL_DISTRIBUTION(IERROR)
 
!  Routine for calculation of spatial photon distributions (with emittance).

	use us_size
	use calc
	use factor
	use beam
	use pinhole
	use angle_alp
	use energy_mod
	use spectra
	implicit none

!  Declarations of scalars:
	integer(kind=i4b) :: IERROR,ICOUNT
	integer(kind=i4b) :: I,IA,IB,IC,IMIN,IMAX,IH1,IH2
        integer(kind=i4b) :: IPC,NSIGMA2,NPC

	real(kind=dp) :: R,DA2,DALPHA,SL,CONST
	real(kind=dp) :: ALPHAI,ALPHA2I,ALPHAMIN,ALPHA2MIN,ALPHAMAX,ALPHA2MAX
        real(kind=dp) :: SIGC,SUM_GAUSS,GSCALE

!  Declarations of arrays:
	real(kind=dp),dimension(:),allocatable :: ALPHA,THETA
        real(kind=dp),dimension(:),allocatable :: ECC,GS
        real(kind=dp),dimension(:,:),allocatable :: RA0C,RA1C,RA2C,RA3C

!  Labeled constants:
	real(kind=dp),parameter :: ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,EPS=1.0D-2
        integer(kind=i4b),parameter :: NSIGMA_DEFAULT=3,NPPSIGMA=12

!
	IERROR = 0
	CONST  = C1

	I1(1) = 0
	I2(1) = 0

	NPC = 1
        ALLOCATE (ECC(1))
	ECC(1) = ZERO
	IF (SIGE .GT. ZERO) THEN	! generate Gaussian with correct sigma in units of x-axis (energy axis)
!	  SIGC = 2.0d0*SIGE*E(1)	! eV
	  NSIGMA2 = 2*NSIGMA_DEFAULT
	  NPC   = DBLE(NSIGMA2)*NPPSIGMA +1.0d0
	  IF (NPC/2*2 .EQ. NPC) NPC = NPC +1		! make odd

	  DEALLOCATE(ECC)
          ALLOCATE (ECC(NPC),GS(NPC))			! allocate energy and Gaussian arrays

          ECC = ((/(IPC,IPC=0,NPC-1)/)/DBLE(NPC-1)-0.5D0)*DBLE(NPC)
          GS = EXP(-ECC**2/2.0D0/DBLE(NPPSIGMA)**2)	! whole array
          SUM_GAUSS = SUM(GS)
          GS = GS/SUM_GAUSS				! normalize
! 	  ECC = ECC*SIGC/DBLE(NPPSIGMA)			! eV; ECC/SIGC gives the extent in units of sigc
	  ECC = -ECC*SIGE/DBLE(NPPSIGMA)		! ECC/SIGE gives the extent in units of sige (use neg. sign to start with high energy)

          ALLOCATE (RA0C(NXP,NYP),RA1C(NXP,NYP),RA2C(NXP,NYP),RA3C(NXP,NYP))

	ENDIF ! SIGE

	DO IPC=1,NPC		! scan over beam energy spread
	GSCALE = ONE+ECC(IPC)
!	R  = ER/(E(1)+ECC(IPC))	! here fixed energy (eV); ECC contains the offsets (save for future reference; GSCALE needs to be set to 1.0)
 	R  = ER/E(1)*GSCALE**2	! change beam energy using gscale

!+ Find range of harmonics that contributes at this energy
!  -----------------------------------------------------------------------------
	IF (METHOD .EQ. 1) THEN ! Finite-N
	  IF (IPC .EQ. 1) ALLOCATE(ALPHA(NALPHA),THETA(NALPHA))
	    DA2  = CALPHA2*R/NPER
	    IMIN = (AP2MIN*GSCALE**2+K3-DA2)/R+ONE
	    IMAX = (AP2MAX*GSCALE**2+K3+DA2)/R
	ELSE ! Infinite-N
	  IF (IPC .EQ. 1) ALLOCATE(ALPHA(1),THETA(1))
	    IMIN = (AP2MIN*GSCALE**2+K3)/R+ONE
	    IMAX = (AP2MAX*GSCALE**2+K3)/R
	END IF ! METHOD
	IF (IMAX .LT. IMIN .AND. IPC .GT. NPC/2+1) GO TO 810	! jump out; ok when beam energy spread is applied
	IF (IMAX .LT. IMIN) GO TO 900
	IF (IHARM .GT. 0 .AND. (IHARM .LT. IMIN .OR. IHARM .GT. IMAX)) GO TO 910
!-
!  -----------------------------------------------------------------------------
	DO IB=1,NYP
	    DO IA=1,NXP
		RA0(IA,IB) = ZERO
		RA1(IA,IB) = ZERO
		RA2(IA,IB) = ZERO
		RA3(IA,IB) = ZERO
	    END DO ! IA
	END DO ! IB

	IF (IHARM .GT. 0) THEN
	    IH1 = IHARM
	    IH2 = IH1
	ELSE IF (IHARM .LT. 0) THEN ! Lowest order
	    IH1 = IMIN
	    IH2 = IH1
	ELSE ! IHARM = 0 ! All contributing harmonics
	    IH1 = IMIN
	    IH2 = IMAX
	END IF ! IHARM
	IF (IPC .EQ. 1) I1(1)  = IH1
	I2(1)  = IH2
	ICOUNT = 0
!+ Loop over harmonics
!  -----------------------------------------------------------------------------
	DO I=IH1,IH2
	    ICOUNT  = ICOUNT+1
	    ALPHA2I = R*I-K3

	    IF (METHOD .EQ. 1) THEN ! Finite-N
		IF (ALPHA2I .GT. ZERO) THEN
		    ALPHAI = SQRT(ALPHA2I)
		ELSE
		    ALPHAI = ZERO
		END IF ! ALPHA2I
		ALPHA2MIN = ALPHA2I-DA2
		IF (ALPHA2MIN .LT. ZERO) ALPHA2MIN = ZERO
		ALPHAMIN  = SQRT(ALPHA2MIN)
		ALPHA2MAX = ALPHA2I+DA2
		ALPHAMAX  = SQRT(ALPHA2MAX)
		DALPHA    = (ALPHAMAX-ALPHAMIN)/NALPHA
		SL = ALPHAMIN+DALPHA/TWO
		DO IC=1,NALPHA
		    ALPHA(IC) = SL+(IC-1)*DALPHA
		    THETA(IC) = ALPHA(IC)/(GAMMA_US*GSCALE)
		END DO ! IC
	    ELSE ! Infinite-N
		ALPHAI   = SQRT(ALPHA2I)
		ALPHA(1) = ALPHAI
		THETA(1) = ALPHA(1)/(GAMMA_US*GSCALE)
		DALPHA   = R/(TWO*NPER)
	    END IF ! METHOD

!  Brilliance
	    CALL BRIGHTNESS_ARRAY(METHOD,I,R,ALPHAI,ALPHA2I,ALPHA)

!  Two-dimensional convolution of the brilliance with the electron distribution
	    CALL CONVOLUTE_DISTRIBUTION(METHOD,CONST,ALPHA,THETA,DALPHA, &
					EPS,BR0,BR1,BR2,BR3, &
					    RA0,RA1,RA2,RA3,ICOUNT)

	    IF (ICOUNT .GT. 1) THEN ! If higher harmonics do not contribute
		I2(1) = I
		GO TO 800
	    END IF ! ICOUNT
		
	END DO ! IH
!- Endloop harmonics
!  -----------------------------------------------------------------------------
800	CONTINUE

	IF (SIGE .GT. ZERO) THEN
	  RA0C = RA0C+GS(IPC)*RA0	! whole arrays
	  RA1C = RA1C+GS(IPC)*RA1
	  RA2C = RA2C+GS(IPC)*RA2
	  RA3C = RA3C+GS(IPC)*RA3
        END IF
	END DO ! IPC

810	CONTINUE
	IF (SIGE .GT. ZERO) THEN
	  RA0 = RA0C	! return in original array
	  RA1 = RA1C
	  RA2 = RA2C
	  RA3 = RA3c
	  DEALLOCATE(GS)
	  DEALLOCATE(RA0C,RA1C,RA2C,RA3C)
        END IF

	DEALLOCATE(ECC)
	DEALLOCATE(ALPHA,THETA)

	RETURN

!+ Error returns
!  -----------------------------------------------------------------------------
900	CONTINUE
	PRINT 200,'&SPATIAL_DISTRIBUTION-E-HARMERR, Harmonic errror'
	PRINT 210,'- no harmonics reachable at ',E(1),' eV.'
	GO TO 999

910	CONTINUE
	PRINT 200,'&SPATIAL_DISTRIBUTION-E-HARMERR, Harmonic errror'
	PRINT 220,'- Harmonic number ',IHARM,' not in reachable range ',IMIN,' to ',IMAX,' at ',E(1), ' eV.'
	GO TO 999

200	FORMAT(' ',8A)
210	FORMAT(' ',A,F10.3,A)
220	FORMAT(' ',3(A,I3),A,F10.3,A)

999	CONTINUE
	IERROR = -1
!-
!  -----------------------------------------------------------------------------
	RETURN
	END ! SPATIAL_DISTRIBUTION
!
	SUBROUTINE SPECTRAL_DISTRIBUTION(IERROR)
 
!  Routine for calculation of photon spectra (with emittance).

	use us_size
	use calc
	use factor
	use beam
	use pinhole
	use angle_alp
	use energy_mod
	use spectra
	use physical_constants !  Fundamental physical constants; Physics Today Aug. 1990:
	use pi_constants

	use econ	! use energy convolution function
	implicit none

!  Declarations of scalars:
	logical(kind=4) :: LE1,LE2

	integer(kind=i4b) :: IERROR,ICOUNT
	integer(kind=i4b) :: I,IA,IB,IC,IE,IMIN,IMAX,IH1,IH2,IEU
	integer(kind=i4b) :: NC

	real(kind=dp) :: R,DA2,DALPHA,SL,CONST
	real(kind=dp) :: ALPHAI,ALPHA2I,ALPHAMIN,ALPHA2MIN,ALPHAMAX,ALPHA2MAX
	real(kind=dp) :: AREA0,AREA1,AREA2,AREA3
	real(kind=dp) :: DEL,SIGR2

!  Declarations of arrays:
	real(kind=dp),dimension(:),allocatable :: ALPHA,THETA
	real(kind=dp),dimension(:),allocatable :: SPEC0C,SPEC1C,SPEC2C,SPEC3C

!  Conversion factors:
	real(kind=dp),parameter :: C_EVANG=H*C/EC*1.0D10,C_ANG_M=1.0D-10 ! 12398.42
	real(kind=dp),parameter :: C_M2_MM2=1.0D6

!  Labeled constants:
	real(kind=dp),parameter :: ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,EPS=1.0D-2
	real(kind=dp),parameter :: CV=ONE/8.0D0/PI/PI*C_ANG_M*C_M2_MM2 ! 1.2665D-6

!
	IERROR = 0
	LE1 = .FALSE.
	LE2 = .FALSE.

	IF (METHOD .EQ. 1) THEN ! Finite-N
	    ALLOCATE(ALPHA(NALPHA),THETA(NALPHA))
	ELSE ! Infinite-N
	    ALLOCATE(ALPHA(1),THETA(1))
	END IF ! METHOD

!+ Loop over energies
!  -----------------------------------------------------------------------------
	DO IE=1,NE

	    SPEC0(IE) = ZERO
	    SPEC1(IE) = ZERO
	    SPEC2(IE) = ZERO
	    SPEC3(IE) = ZERO
	    I1(IE)    = 0
	    I2(IE)    = 0
	    CONST     = C1

	    R         = ER/E(IE)

!+ Find range of harmonics that contributes at energy E(IE)
!  -----------------------------------------------------------------------------
	    IF (METHOD .EQ. 1) THEN ! Finite-N
		DA2  = CALPHA2*R/NPER
		IMIN = (AP2MIN+K3-DA2)/R+ONE
		IMAX = (AP2MAX+K3+DA2)/R
	    ELSE ! Infinite-N
		IMIN = (AP2MIN+K3)/R+ONE
		IMAX = (AP2MAX+K3)/R
	    END IF ! METHOD
	    IF (IMAX .LT. IMIN) GO TO 810 ! Next energy
	    LE1 = .TRUE.		
	    IF (IHARM .GT. 0 .AND. (IHARM .LT. IMIN .OR. IHARM .GT. IMAX)) GO TO 810 ! Next energy
	    LE2 = .TRUE.
!-
!  -----------------------------------------------------------------------------
	    DO IB=1,NYP
		DO IA=1,NXP
		    RA0(IA,IB) = ZERO
		    RA1(IA,IB) = ZERO
		    RA2(IA,IB) = ZERO
		    RA3(IA,IB) = ZERO
		END DO ! IA
	    END DO ! IB

	    IF (IHARM .GT. 0) THEN
		IH1 = IHARM
		IH2 = IH1
	    ELSE IF (IHARM .LT. 0) THEN ! Lowest order
		IH1 = IMIN
		IH2 = IH1
	    ELSE ! IHARM = 0 ! All contributing harmonics
		IH1 = IMIN
		IH2 = IMAX
	    END IF ! IHARM
	    I1(IE)  = IH1
	    I2(IE)  = IH2
	    ICOUNT = 0
!+ Loop over harmonics
!  -----------------------------------------------------------------------------
	    DO I=IH1,IH2
		ICOUNT  = ICOUNT+1
		ALPHA2I = R*I-K3

		IF (METHOD .EQ. 1) THEN ! Finite-N
		    IF (ALPHA2I .GT. ZERO) THEN
			ALPHAI = SQRT(ALPHA2I)
		    ELSE
			ALPHAI = ZERO
		    END IF ! ALPHA2I
		    ALPHA2MIN = ALPHA2I-DA2
		    IF (ALPHA2MIN .LT. ZERO) ALPHA2MIN = ZERO
		    ALPHAMIN  = SQRT(ALPHA2MIN)
		    ALPHA2MAX = ALPHA2I+DA2
		    ALPHAMAX  = SQRT(ALPHA2MAX)
		    DALPHA    = (ALPHAMAX-ALPHAMIN)/NALPHA
		    SL = ALPHAMIN+DALPHA/TWO
		    DO IC=1,NALPHA
			ALPHA(IC) = SL+(IC-1)*DALPHA
			THETA(IC) = ALPHA(IC)/GAMMA_US
		    END DO ! IC
		ELSE ! Infinite-N
		    ALPHAI   = SQRT(ALPHA2I)
		    ALPHA(1) = ALPHAI
		    THETA(1) = ALPHA(1)/GAMMA_US
		    DALPHA   = R/(TWO*NPER)
		END IF ! METHOD

	        IF (MODE .EQ. 3) THEN ! Brilliance
		    IF (ALPHA2I .LT. ZERO) THEN ! may be neg for finite-N only
		        DEL = ZERO
		    ELSE                        ! Walker's approach (infinite-N)
		        DEL = ALPHA2I*NPER/R
		    END IF ! ALPHA2I

!                   Estimate diffraction limited source size
		    IF (DEL .LT. 2.15D0) THEN
		        SIGR2 = (1.29D0+(1.229D0*(DEL-0.8D0)*(DEL-0.8D0)))**2
		    ELSE
                        SIGR2 = 5.81D0*DEL
		    END IF ! DEL

		    SIGR2 = SIGR2*CV*C_EVANG/E(IE)*LDEV
		    CONST = C1/(TWOPI*SQRT((SIGX2+SIGR2)*(SIGY2+SIGR2)))
	        END IF ! MODE

!  Brilliance
		CALL BRIGHTNESS_ARRAY(METHOD,I,R,ALPHAI,ALPHA2I,ALPHA)

!  Two-dimensional convolution of the brilliance with the electron distribution
		CALL CONVOLUTE_DISTRIBUTION(METHOD,CONST,ALPHA,THETA,DALPHA, &
					    EPS,BR0,BR1,BR2,BR3, &
					        RA0,RA1,RA2,RA3,ICOUNT)

		IF (ICOUNT .GT. 1) THEN ! If higher harmonics do not contribute
		    I2(IE) = I
		    GO TO 800
		END IF ! ICOUNT
	    END DO ! IH
!- Endloop harmonics
!  -----------------------------------------------------------------------------
800	    CONTINUE

!+ Save spectra
!  -----------------------------------------------------------------------------
	    IF (MODE .EQ. 4) THEN ! Integrate over pinhole
		CALL TRAPZ2(RA0,AREA0)
		CALL TRAPZ2(RA1,AREA1)
		CALL TRAPZ2(RA2,AREA2)
		CALL TRAPZ2(RA3,AREA3)
		SPEC0(IE) = FAC*AREA0
		SPEC1(IE) = FAC*AREA1
		SPEC2(IE) = FAC*AREA2
		SPEC3(IE) = FAC*AREA3
	    ELSE ! Fixed position or brilliance
		SPEC0(IE) = FAC*RA0(1,1)
		SPEC1(IE) = FAC*RA1(1,1)
		SPEC2(IE) = FAC*RA2(1,1)
		SPEC3(IE) = FAC*RA3(1,1)
	    END IF ! MODE
	    IF (FAC .EQ. 4.0D0) SPEC2(IE) = ZERO ! symmetry
!-
!  -----------------------------------------------------------------------------

810	    CONTINUE
	END DO ! IE
!- Endloop energy
!  -----------------------------------------------------------------------------
	
	DEALLOCATE(ALPHA)
	DEALLOCATE(THETA)

	IF (.NOT. LE1 .AND. .NOT. LE2) GO TO 900
	IF (      LE1 .AND. .NOT. LE2) GO TO 910

	IF (METHOD .EQ. 4) THEN ! Infinite-N with convolution
	    CALL CONVOLUTE_ENERGY_VSTEP(SPEC0)
	    CALL CONVOLUTE_ENERGY_VSTEP(SPEC1)
	    CALL CONVOLUTE_ENERGY_VSTEP(SPEC2)
	    CALL CONVOLUTE_ENERGY_VSTEP(SPEC3)
!  Apply beam energy spread
	    IF (SIGE .GT. ZERO) THEN
		ALLOCATE(SPEC0C(NEU),SPEC1C(NEU),SPEC2C(NEU),SPEC3C(NEU))
		SPEC0C = ECON_FUNC(EU,SPEC0,NEU,SIGE,IERROR)
		IF (IERROR .NE. 0) GOTO 920
		SPEC1C = ECON_FUNC(EU,SPEC1,NEU,SIGE,IERROR)
		IF (IERROR .NE. 0) GOTO 921
		SPEC2C = ECON_FUNC(EU,SPEC2,NEU,SIGE,IERROR)
		IF (IERROR .NE. 0) GOTO 922
		SPEC3C = ECON_FUNC(EU,SPEC3,NEU,SIGE,IERROR)
		IF (IERROR .NE. 0) GOTO 923
		FORALL (IEU=1:NEU) SPEC0(IEU) = SPEC0C(IEU)	! save in original array
		FORALL (IEU=1:NEU) SPEC1(IEU) = SPEC1C(IEU)	! save in original array
		FORALL (IEU=1:NEU) SPEC2(IEU) = SPEC2C(IEU)	! save in original array
		FORALL (IEU=1:NEU) SPEC3(IEU) = SPEC3C(IEU)	! save in original array
		DEALLOCATE(SPEC0C,SPEC1C,SPEC2C,SPEC3C)
	    ENDIF
	ELSE IF (METHOD .EQ. 14) THEN ! Infinite-N with convolution (Walker's)
	    CALL CONVOLUTE_ENERGY_ESTEP(SPEC0)
	    CALL CONVOLUTE_ENERGY_ESTEP(SPEC1)
	    CALL CONVOLUTE_ENERGY_ESTEP(SPEC2)
	    CALL CONVOLUTE_ENERGY_ESTEP(SPEC3)
!  Apply beam energy spread
	    IF (SIGE .GT. ZERO) THEN
	        NC = NE2 -NE1 +1
		ALLOCATE(SPEC0C(NC),SPEC1C(NC),SPEC2C(NC),SPEC3C(NC))
		SPEC0C = ECON_FUNC(E(NE1:NE2),SPEC0(NE1:NE2),NC,SIGE,IERROR)
		IF (IERROR .NE. 0) GOTO 920
		SPEC1C = ECON_FUNC(E(NE1:NE2),SPEC1(NE1:NE2),NC,SIGE,IERROR)
		IF (IERROR .NE. 0) GOTO 921
		SPEC2C = ECON_FUNC(E(NE1:NE2),SPEC2(NE1:NE2),NC,SIGE,IERROR)
		IF (IERROR .NE. 0) GOTO 922
		SPEC3C = ECON_FUNC(E(NE1:NE2),SPEC3(NE1:NE2),NC,SIGE,IERROR)
		IF (IERROR .NE. 0) GOTO 923
		FORALL (IE=1:NC) SPEC0(NE1-1+IE) = SPEC0C(IE)	! save in original array
		FORALL (IE=1:NC) SPEC1(NE1-1+IE) = SPEC1C(IE)	! save in original array
		FORALL (IE=1:NC) SPEC2(NE1-1+IE) = SPEC2C(IE)	! save in original array
		FORALL (IE=1:NC) SPEC3(NE1-1+IE) = SPEC3C(IE)	! save in original array
		DEALLOCATE(SPEC0C,SPEC1C,SPEC2C,SPEC3C)
	    ENDIF
!  Apply beam energy spread
	ELSE ! methods = 1,2
	    IF (SIGE .GT. ZERO) THEN
		ALLOCATE(SPEC0C(NE),SPEC1C(NE),SPEC2C(NE),SPEC3C(NE))
		SPEC0C = ECON_FUNC(E,SPEC0,NE,SIGE,IERROR)
		IF (IERROR .NE. 0) GOTO 920
		SPEC1C = ECON_FUNC(E,SPEC1,NE,SIGE,IERROR)
		IF (IERROR .NE. 0) GOTO 921
		SPEC2C = ECON_FUNC(E,SPEC2,NE,SIGE,IERROR)
		IF (IERROR .NE. 0) GOTO 922
		SPEC3C = ECON_FUNC(E,SPEC3,NE,SIGE,IERROR)
		IF (IERROR .NE. 0) GOTO 923
		FORALL (IE=1:NE) SPEC0(IE) = SPEC0C(IE)	! save in original array
		FORALL (IE=1:NE) SPEC1(IE) = SPEC1C(IE)	! save in original array
		FORALL (IE=1:NE) SPEC2(IE) = SPEC2C(IE)	! save in original array
		FORALL (IE=1:NE) SPEC3(IE) = SPEC3C(IE)	! save in original array
		DEALLOCATE(SPEC0C,SPEC1C,SPEC2C,SPEC3C)
	    ENDIF
	END IF 

	RETURN

!+ Error returns
!  -----------------------------------------------------------------------------
900	CONTINUE
	PRINT 200,'&SPECTRAL_DISTRIBUTION-E-HARMERR, Harmonic errror'
	PRINT 210,'- no harmonics reachable in the range ',E(NE1),' to ',E(NE2),' eV.'
	GO TO 999

910	CONTINUE
	PRINT 200,'&SPECTRAL_DISTRIBUTION-E-HARMERR, Harmonic errror'
	PRINT 220,'- Harmonic number ',IHARM,' not in the range ',E(NE1),' to ',E(NE2),' eV.'
	GO TO 999

920	CONTINUE
	PRINT 200,'&SPECTRAL_DISTRIBUTION-E-CONVERR, Energy convolution errror for SPEC0'
	GO TO 999

921	CONTINUE
	PRINT 200,'&SPECTRAL_DISTRIBUTION-E-CONVERR, Energy convolution errror for SPEC1'
	GO TO 999

922	CONTINUE
	PRINT 200,'&SPECTRAL_DISTRIBUTION-E-CONVERR, Energy convolution errror for SPEC2'
	GO TO 999

923	CONTINUE
	PRINT 200,'&SPECTRAL_DISTRIBUTION-E-CONVERR, Energy convolution errror for SPEC3'
	GO TO 999

200	FORMAT(' ',8A)
210	FORMAT(' ',A,2(F10.3,A))
220	FORMAT(' ',A,I3,A,2(F10.3,A))

999	CONTINUE
	IERROR = -1
!-
!  -----------------------------------------------------------------------------
	RETURN
	END ! SPECTRAL_DISTRIBUTION
!
	SUBROUTINE ANGLE_INTEGRATION(IERROR)
 
!  Routine for calculation of angle-integrated photon spectra (with emittance).

	use calc
	use factor
	use beam
	use angle_phi
	use energy_mod
	use spectra

	use econ	! use energy convolution function
	implicit none

!  Declarations of scalars:
	integer(kind=i4b) :: IERROR,ICOUNT
	integer(kind=i4b) :: I,ID,IE,IMIN,IEU
	integer(kind=i4b) :: NC

	real(kind=dp) :: R,CONST
	real(kind=dp) :: ALPHAI,ALPHA2I
	real(kind=dp) :: SUM0,SUM1,SUM3

!  Declarations of arrays:
	real(kind=dp),dimension(1) :: SL
	real(kind=dp),dimension(:),allocatable :: SPEC0C,SPEC1C,SPEC3C

!  Labeled constants:
	real(kind=dp),parameter :: ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,FOUR=4.0D0,HALF=0.5D0,EPS=1.0D-2

!
	IERROR = 0
	SL(1)  = ZERO

!+ Loop over energies
!  -----------------------------------------------------------------------------
	DO IE=1,NE

	    SPEC0(IE) = ZERO
	    SPEC1(IE) = ZERO
	    SPEC2(IE) = ZERO
	    SPEC3(IE) = ZERO
	    I1(IE)    = 0
	    I2(IE)    = 0

	    R     = ER/E(IE)
	    CONST = FOUR*C2*R/(TWO*NPER)*DPHI

	    IMIN = K3/R+ONE
	    IF (IHARM .GT. 0 .AND. IHARM .LT. IMIN) GO TO 810 ! Next energy

	    IF (IHARM .GT. 0) THEN
		I = IHARM
	    ELSE
		I = IMIN
	    END IF ! IHARM

	    I1(IE) = I
	    ICOUNT = 0

!+ Loop over harmonics
!  -----------------------------------------------------------------------------
	    DO WHILE (.TRUE.)
		ICOUNT  = ICOUNT+1
		ALPHA2I = R*I-K3
		ALPHAI  = SQRT(ALPHA2I)
!  Brilliance
		CALL BRIGHTNESS_ARRAY(METHOD,I,R,ALPHAI,ALPHA2I,SL)

!  Integrate
		SUM0 = ZERO
		SUM1 = ZERO
		SUM3 = ZERO
		DO ID=1,NPHI_BRIGHT
		    SUM0 = SUM0+BR0(ID,1)
		    SUM1 = SUM1+BR1(ID,1)
		    SUM3 = SUM3+BR3(ID,1)
		END DO ! NPHI_BRIGHT

		SUM0 = SUM0*CONST
		SUM1 = SUM1*CONST
		SUM3 = SUM3*CONST
		SPEC0(IE) = SPEC0(IE)+SUM0
		SPEC1(IE) = SPEC1(IE)+SUM1
		SPEC3(IE) = SPEC3(IE)+SUM3

		IF (IHARM .NE. 0) GO TO 800

		IF (SUM0 .GT. EPS*SPEC0(IE)) ICOUNT = 0
		IF (ICOUNT .GT. 1) GO TO 800 ! If higher harmonics do not contribute

		I = I+1
	    END DO ! WHILE
!- Endloop harmonics
!  -----------------------------------------------------------------------------
800	    CONTINUE
	    I2(IE) = I

810	    CONTINUE
	END DO ! IE
!- Endloop energy
!  -----------------------------------------------------------------------------
	
	IF (METHOD .EQ. 4) THEN ! Infinite-N with convolution
	    CALL CONVOLUTE_ENERGY_VSTEP(SPEC0)
	    CALL CONVOLUTE_ENERGY_VSTEP(SPEC1)
	    CALL CONVOLUTE_ENERGY_VSTEP(SPEC3)
!  Apply beam energy spread
	    IF (SIGE .GT. ZERO) THEN
		ALLOCATE(SPEC0C(NEU),SPEC1C(NEU),SPEC3C(NEU))
		SPEC0C = ECON_FUNC(EU,SPEC0,NEU,SIGE,IERROR)
		IF (IERROR .NE. 0) GOTO 920
		SPEC1C = ECON_FUNC(EU,SPEC1,NEU,SIGE,IERROR)
		IF (IERROR .NE. 0) GOTO 921
		SPEC3C = ECON_FUNC(EU,SPEC3,NEU,SIGE,IERROR)
		IF (IERROR .NE. 0) GOTO 923
		FORALL (IEU=1:NEU) SPEC0(IEU) = SPEC0C(IEU)	! save in original array
		FORALL (IEU=1:NEU) SPEC1(IEU) = SPEC1C(IEU)	! save in original array
		FORALL (IEU=1:NEU) SPEC3(IEU) = SPEC3C(IEU)	! save in original array
		DEALLOCATE(SPEC0C,SPEC1C,SPEC3C)
	    ENDIF
	ELSE IF (METHOD .EQ. 14) THEN ! Infinite-N with convolution (Walker's)
	    CALL CONVOLUTE_ENERGY_ESTEP(SPEC0)
	    CALL CONVOLUTE_ENERGY_ESTEP(SPEC1)
	    CALL CONVOLUTE_ENERGY_ESTEP(SPEC3)
!  Apply beam energy spread
	    IF (SIGE .GT. ZERO) THEN
	        NC = NE2 -NE1 +1
		ALLOCATE(SPEC0C(NC),SPEC1C(NC),SPEC3C(NC))
		SPEC0C = ECON_FUNC(E(NE1:NE2),SPEC0(NE1:NE2),NC,SIGE,IERROR)
		IF (IERROR .NE. 0) GOTO 920
		SPEC1C = ECON_FUNC(E(NE1:NE2),SPEC1(NE1:NE2),NC,SIGE,IERROR)
		IF (IERROR .NE. 0) GOTO 921
		SPEC3C = ECON_FUNC(E(NE1:NE2),SPEC3(NE1:NE2),NC,SIGE,IERROR)
		IF (IERROR .NE. 0) GOTO 923
		FORALL (IE=1:NC) SPEC0(NE1-1+IE) = SPEC0C(IE)	! save in original array
		FORALL (IE=1:NC) SPEC1(NE1-1+IE) = SPEC1C(IE)	! save in original array
		FORALL (IE=1:NC) SPEC3(NE1-1+IE) = SPEC3C(IE)	! save in original array
		DEALLOCATE(SPEC0C,SPEC1C,SPEC3C)
	    ENDIF
!  Apply beam energy spread
	ELSE ! method = 2
	    IF (SIGE .GT. ZERO) THEN
		ALLOCATE(SPEC0C(NE),SPEC1C(NE),SPEC3C(NE))
		SPEC0C = ECON_FUNC(E,SPEC0,NE,SIGE,IERROR)
		IF (IERROR .NE. 0) GOTO 920
		SPEC1C = ECON_FUNC(E,SPEC1,NE,SIGE,IERROR)
		IF (IERROR .NE. 0) GOTO 921
		SPEC3C = ECON_FUNC(E,SPEC3,NE,SIGE,IERROR)
		IF (IERROR .NE. 0) GOTO 923
		FORALL (IE=1:NE) SPEC0(IE) = SPEC0C(IE)	! save in original array
		FORALL (IE=1:NE) SPEC1(IE) = SPEC1C(IE)	! save in original array
		FORALL (IE=1:NE) SPEC3(IE) = SPEC3C(IE)	! save in original array
		DEALLOCATE(SPEC0C,SPEC1C,SPEC3C)
	    ENDIF
	END IF 

	RETURN

!+ Error returns
!  -----------------------------------------------------------------------------
!900	CONTINUE
!	GO TO 999
!
920	CONTINUE
	PRINT 200,'&ANGLE-INTEGRATION-E-CONVERR, Energy convolution errror for SPEC0'
	GO TO 999

921	CONTINUE
	PRINT 200,'&ANGLE-INTEGRATION-E-CONVERR, Energy convolution errror for SPEC1'
	GO TO 999

923	CONTINUE
	PRINT 200,'&ANGLE-INTEGRATION-E-CONVERR, Energy convolution errror for SPEC3'
	GO TO 999

200	FORMAT(' ',8A)

999	CONTINUE
	IERROR = -1
!-
!  -----------------------------------------------------------------------------
!	RETURN
	END ! ANGLE_INTEGRATION
!
	SUBROUTINE POWER_DISTRIBUTION(IERROR)
     
!  Routine for calculation of power distributions (with and without emittance).

	use us_size
	use calc
	use factor
	use beam
	use pinhole
	use energy_mod
	use spectra
	use power
	implicit none

!  Declarations of scalars:
	integer(kind=i4b) :: IERROR,ICOUNT
	integer(kind=i4b) :: I,IA,IB,IC
	integer(kind=i4b) :: NXE,NYE

	real(kind=dp) :: CONST,DELTAP,PTOT
	real(kind=dp) :: ALPHAP,ALPHA2
	real(kind=dp) :: DXE,DYE
	real(kind=dp) :: XG,YG,COSPHI,SINPHI
	real(kind=dp) :: S0,S1,S2,S3,AXR,AYR,AXI,AYI
	real(kind=dp) :: AREA0,DELTA0,DELTA1,DELTA2,DELTA3

!  Declarations of arrays:

!  Labeled constants:
	real(kind=dp),parameter :: ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,EPS=1.0D-6,EPSP=2.0D-3
	integer(kind=i4b),parameter :: NPPSIGMA_PDF=2 ! 2 is good because the sigma is typically small compared to the width of the pdf

!
	IERROR = 0
	PTOT   = ZERO

        DO IB=1,NYP
            DO IA=1,NXP
                RA0(IA,IB) = ZERO
	        RA1(IA,IB) = ZERO
	        RA2(IA,IB) = ZERO
	        RA3(IA,IB) = ZERO
	    END DO ! IA
	END DO ! IB

	IF (METHOD .NE. 3) THEN ! non-zero emittance; define mesh size in x and y
	    CONST = C5
	    DXE = SIGU/NPPSIGMA_PDF
	    DYE = SIGV/NPPSIGMA_PDF
	    NXE = CX(NXP)/DXE +ONE +NSIGMA
	    NYE = CY(NYP)/DYE +ONE +NSIGMA
	    IF (NXE .GT. PW_SZ) NXE = PW_SZ
	    IF (NYE .GT. PW_SZ) NYE = PW_SZ
	    ALLOCATE(PW0(0:NXE,0:NYE),PW1(0:NXE,0:NYE),PW2(0:NXE,0:NYE),PW3(0:NXE,0:NYE))
	ELSE ! zero emittance
	    CONST = C4
	END IF ! METHOD

	IF (IHARM .GT. 0) THEN
	    I = IHARM
	ELSE
	    I = 1
	END IF ! IHARM
	I1(1)  = I
	ICOUNT = 0

!+ Loop over harmonics
!  -----------------------------------------------------------------------------
	DO WHILE (.TRUE.)
	    ICOUNT = ICOUNT+1

	    IF (METHOD .NE. 3) THEN ! non-zero emittance; convolve with electron distribution

!  Precalculate power density array
	        CALL PDF_ARRAY(I,DXE,DYE,NXE,NYE)
!  Two-dimensional convolution of the pdf with the electron distribution (cartesian coordinates)
		CALL CONVOLUTE_POWER_DISTRIBUTION(CONST,DXE,DYE,NXE,NYE, &
					          EPSP,PW0,PW1,PW2,PW3, &
					          RA0,RA1,RA2,RA3,ICOUNT)
	    ELSE ! zero emittance; direct calculation
	        DO IB=1,NYP
		    YG     = GAMMA_US*CY(IB)
	            DO IA=1,NXP
		        XG     = GAMMA_US*CX(IA)
		        ALPHA2 = XG*XG+YG*YG
		        ALPHAP = SQRT(ALPHA2)
		        IF (ALPHAP .LT. EPS) THEN
		            COSPHI = ZERO
			    SINPHI = ONE
		        ELSE
		            COSPHI = XG/ALPHAP
		            SINPHI = YG/ALPHAP
		        END IF ! ALPHAP
		        CALL BRIGHTE(I,ALPHAP,COSPHI,SINPHI,S0,S1,S2,S3,AXR,AYR,AXI,AYI)
			DELTA0 = CONST*S0/(K3+ALPHA2)
			DELTA1 = CONST*S1/(K3+ALPHA2)
			DELTA2 = CONST*S2/(K3+ALPHA2)
			DELTA3 = CONST*S3/(K3+ALPHA2)
			RA0(IA,IB) = RA0(IA,IB)+DELTA0
			RA1(IA,IB) = RA1(IA,IB)+DELTA1
			RA2(IA,IB) = RA2(IA,IB)+DELTA2
			RA3(IA,IB) = RA3(IA,IB)+DELTA3
			IF (DELTA0 .GT. EPSP*RA0(IA,IB)) ICOUNT = 0
		    END DO ! IA
		END DO ! IB
	    END IF ! METHOD

	    CALL TRAPZ2(RA0,AREA0) ! integrated power for harmonics including i (W)
	    DELTAP = FAC*AREA0-PTOT! incremental change in total power (W)
	    IF (DELTAP .GT. EPSP*PTOT) ICOUNT = 0
	    PTOT = PTOT+DELTAP

	    IF (IHARM .GT. 0) GO TO 800
	    IF (IHARM .LT. 0 .AND. I .GE. -IHARM) GO TO 800
	    IF (IHARM .EQ. 0 .AND. ICOUNT .GT. 1) GO TO 800

	    I = I+1
	END DO ! WHILE
!+ Endloop harmonics
!  -----------------------------------------------------------------------------
800	CONTINUE
        I2(1) = I
	IF (ALLOCATED(PW0)) DEALLOCATE(PW0)
	IF (ALLOCATED(PW1)) DEALLOCATE(PW1)
	IF (ALLOCATED(PW2)) DEALLOCATE(PW2)
	IF (ALLOCATED(PW3)) DEALLOCATE(PW3)

	RETURN
	END ! POWER_DISTRIBUTION
!
	SUBROUTINE NO_EMITTANCE(IERROR)
 
!  Routine for calculation of zero emittance quantities (both spectra and spatial distributions).

	use calc
	use factor
	use beam
	use pinhole
	use angle_alp
	use energy_mod
	use spectra
	use physical_constants !  Fundamental physical constants; Physics Today Aug. 1990:
	use pi_constants

	use econ	! use energy convolution function
	implicit none

!  Declarations of scalars:
	logical(kind=4) :: LE

	integer(kind=i4b) :: IERROR
	integer(kind=i4b) :: I,IA,IB,IE
        integer(kind=i4b) :: IPC,NSIGMA2,NPC

	real(kind=dp) :: R,DA2,CONST,ARG
	real(kind=dp) :: ALPHA,ALPHA2,ALPHA2I,ALPHA2MIN,ALPHA2MAX,SR2
	real(kind=dp) :: XG,YG,COSPHI,SINPHI
	real(kind=dp) :: S0,S1,S2,S3,AXR,AYR,AXI,AYI
	real(kind=dp) :: AREA0,AREA1,AREA2,AREA3
	real(kind=dp) :: SINC
        real(kind=dp) :: SIGC,SUM_GAUSS,GSCALE

!  Declarations of arrays:
	real(kind=dp),dimension(:),allocatable :: SPEC0C,SPEC1C,SPEC2C,SPEC3C
        real(kind=dp),dimension(:),allocatable :: ECC,GS
        real(kind=dp),dimension(:,:),allocatable :: RA0C,RA1C,RA2C,RA3C

!  Conversion factors:
	real(kind=dp),parameter :: C_EVANG=H*C/EC*1.0D10,C_ANG_M=1.0D-10
	real(kind=dp),parameter :: C_M2_MM2=1.0D6

!  Labeled constants:
	real(kind=dp),parameter :: ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,EPS=1.0D-6
	real(kind=dp),parameter :: CV=ONE/8.0D0/PI/PI*C_ANG_M*C_M2_MM2,CW=3.57D0 ! 1.2665D-6
        integer(kind=i4b),parameter :: NSIGMA_DEFAULT=3,NPPSIGMA=12

!
	IERROR = 0
	LE     = .FALSE.

	NPC = 1
        ALLOCATE (ECC(1))
	ECC(1) = ZERO
	IF (MODE .EQ. 1 .AND. SIGE .GT. ZERO) THEN ! Spatial distributions w/ beam energy spread; generate Gaussian with correct sigma in units of x-axis (energy axis)
	    NSIGMA2 = 2*NSIGMA_DEFAULT
	    NPC   = DBLE(NSIGMA2)*NPPSIGMA +1.0d0
	    IF (NPC/2*2 .EQ. NPC) NPC = NPC +1		! make odd

	    DEALLOCATE(ECC)
	    ALLOCATE (ECC(NPC),GS(NPC))			! allocate energy and Gaussian arrays

	    ECC = ((/(IPC,IPC=0,NPC-1)/)/DBLE(NPC-1)-0.5D0)*DBLE(NPC)
	    GS = EXP(-ECC**2/2.0D0/DBLE(NPPSIGMA)**2)	! whole array
	    SUM_GAUSS = SUM(GS)
	    GS = GS/SUM_GAUSS				! normalize
	    ECC = -ECC*SIGE/DBLE(NPPSIGMA)		! ECC/SIGE gives the extent in units of sige (use neg. sign to start with high energy)

	    ALLOCATE (RA0C(NXP,NYP),RA1C(NXP,NYP),RA2C(NXP,NYP),RA3C(NXP,NYP))

	ENDIF ! MODE;SIGE

!+ Loop over energies
!  -----------------------------------------------------------------------------
	DO IE=1,NE

	    SPEC0(IE) = ZERO
	    SPEC1(IE) = ZERO
	    SPEC2(IE) = ZERO
	    SPEC3(IE) = ZERO
	    I1(IE)    = 0
	    I2(IE)    = 0

	    IF (MODE .EQ. 3) THEN ! Brilliance
		SR2   = CW*CV*C_EVANG/E(IE)*LDEV
		CONST = C3/(TWOPI*SR2)
	    ELSE
		CONST = C3
	    END IF ! MODE

	    DO IPC=1,NPC		! scan over beam energy spread
	    GSCALE = ONE+ECC(IPC)

	    R   = ER/E(IE)*GSCALE**2	! change beam energy using gscale
	    DA2 = CALPHA2*R/NPER

	    DO IB=1,NYP
		YG     = GAMMA_US*GSCALE*CY(IB)
		DO IA=1,NXP
		    XG     = GAMMA_US*GSCALE*CX(IA)
		    ALPHA2 = XG*XG+YG*YG

		    RA0(IA,IB) = ZERO
		    RA1(IA,IB) = ZERO
		    RA2(IA,IB) = ZERO
		    RA3(IA,IB) = ZERO

		    IF (IHARM .GT. 0) THEN
			I	  = IHARM
			ALPHA2I   = R*I-K3
			ALPHA2MAX = ALPHA2I+DA2
			IF (ALPHA2MAX .LT. ALPHA2) GO TO 800
		    ELSE ! Find harmonic #
			I         = 1
			ALPHA2I   = R*I-K3
			ALPHA2MAX = ALPHA2I+DA2
			DO WHILE (ALPHA2MAX .LT. ALPHA2)
			    I         = I+1
			    ALPHA2I   = R*I-K3
			    ALPHA2MAX = ALPHA2I+DA2
			END DO ! WHILE
		    END IF ! IHARM
			
		    ALPHA2MIN = ALPHA2I-DA2
		    IF (ALPHA2MIN .LE. ALPHA2) THEN
			LE = .TRUE.
			IF (IPC .EQ. 1) THEN
			  IF (I1(IE) .EQ. 0) THEN
			      I1(IE) = I
			  ELSE IF (I .LT. I1(IE)) THEN
			      I1(IE) = I
			  END IF ! I1
		        END IF ! IPC
			I2(IE) = I
			ALPHA  = SQRT(ALPHA2)
			IF (ALPHA .LT. EPS) THEN
			    COSPHI = ZERO
			    SINPHI = ONE
			ELSE
			    COSPHI = XG/ALPHA
			    SINPHI = YG/ALPHA
			END IF ! ALPHA
	        	CALL BRIGHTE(I,ALPHA,COSPHI,SINPHI,S0,S1,S2,S3,AXR,AYR,AXI,AYI)
			ARG        = NPI*(ALPHA2-ALPHA2I)/R
			RA0(IA,IB) = CONST*SINC(ARG)*S0
			RA1(IA,IB) = CONST*SINC(ARG)*S1
			RA2(IA,IB) = CONST*SINC(ARG)*S2
			RA3(IA,IB) = CONST*SINC(ARG)*S3
		    END IF ! ALPHA2MIN

800		    CONTINUE
		END DO ! IA
	    END DO ! IB

!+ Save spectra
!  -----------------------------------------------------------------------------
	    IF (MODE .EQ. 4) THEN ! Integrate over pinhole
		CALL TRAPZ2(RA0,AREA0)
		CALL TRAPZ2(RA1,AREA1)
		CALL TRAPZ2(RA2,AREA2)
		CALL TRAPZ2(RA3,AREA3)
		SPEC0(IE) = FAC*AREA0
		SPEC1(IE) = FAC*AREA1
		SPEC2(IE) = FAC*AREA2
		SPEC3(IE) = FAC*AREA3
		IF (FAC .EQ. 4.0D0) SPEC2(IE) = ZERO
	    ELSE IF (MODE .EQ. 2 .OR. MODE .EQ. 3) THEN ! Fixed position or brilliance
		SPEC0(IE) = RA0(1,1)
		SPEC1(IE) = RA1(1,1)
		SPEC2(IE) = RA2(1,1)
		SPEC3(IE) = RA3(1,1)
	    END IF ! MODE (else MODE = 1)
!-
!  -----------------------------------------------------------------------------

	IF (MODE .EQ. 1 .AND. SIGE .GT. ZERO) THEN
	  RA0C = RA0C+GS(IPC)*RA0	! whole arrays
	  RA1C = RA1C+GS(IPC)*RA1
	  RA2C = RA2C+GS(IPC)*RA2
	  RA3C = RA3C+GS(IPC)*RA3
        END IF
	END DO ! IPC

	IF (MODE .EQ. 1 .AND. SIGE .GT. ZERO) THEN
	  RA0 = RA0C	! return in original array
	  RA1 = RA1C
	  RA2 = RA2C
	  RA3 = RA3c
	  DEALLOCATE(GS)
	  DEALLOCATE(RA0C,RA1C,RA2C,RA3C)
        END IF

	END DO ! IE
!- Endloop energy
!  -----------------------------------------------------------------------------

	IF (.NOT. LE .AND. MODE .EQ. 1) GO TO 900
	IF (.NOT. LE .AND. MODE .NE. 1) GO TO 910

!  Apply beam energy spread (for MODE not equal 1)
	IF (MODE .NE. 1 .AND. SIGE .GT. ZERO) THEN
	    ALLOCATE(SPEC0C(NE),SPEC1C(NE),SPEC2C(NE),SPEC3C(NE))
	    SPEC0C = ECON_FUNC(E,SPEC0,NE,SIGE,IERROR)
	    IF (IERROR .NE. 0) GOTO 920
	    SPEC1C = ECON_FUNC(E,SPEC1,NE,SIGE,IERROR)
	    IF (IERROR .NE. 0) GOTO 921
	    SPEC2C = ECON_FUNC(E,SPEC2,NE,SIGE,IERROR)
	    IF (IERROR .NE. 0) GOTO 922
	    SPEC3C = ECON_FUNC(E,SPEC3,NE,SIGE,IERROR)
	    IF (IERROR .NE. 0) GOTO 923
	    FORALL (IE=1:NE) SPEC0(IE) = SPEC0C(IE)	! save in original array
	    FORALL (IE=1:NE) SPEC1(IE) = SPEC1C(IE)	! save in original array
	    FORALL (IE=1:NE) SPEC2(IE) = SPEC2C(IE)	! save in original array
	    FORALL (IE=1:NE) SPEC3(IE) = SPEC3C(IE)	! save in original array
	    DEALLOCATE(SPEC0C,SPEC1C,SPEC2C,SPEC3C)
	ENDIF

	DEALLOCATE(ECC)
	RETURN

!+ Error returns
!  -----------------------------------------------------------------------------
900	CONTINUE
	PRINT 200,'&NO_EMITTANCE-E-HARMERR, Harmonic errror' 
	IF (IHARM .GT. 0) THEN
	    PRINT 220,'- Harmonic number ',IHARM,' not reachable at ', E(1),' eV.'
	ELSE
	    PRINT 210,'- no harmonics reachable at ',E(1),' eV.'
	END IF ! IHARM
	GO TO 999

910	CONTINUE
	PRINT 200,'&NO_EMITTANCE-E-HARMERR, Harmonic errror' 
	IF (IHARM .GT. 0) THEN
	    PRINT 220,'- Harmonic number ',IHARM,' not in the range ', E(NE1),' to ',E(NE2),' eV.'
	ELSE
	    PRINT 210,'- no harmonics reachable in the range ',E(NE1), ' to ',E(NE2),' eV.'
	END IF ! IHARM
	GO TO 999

920	CONTINUE
	PRINT 200,'&NO-EMITTANCE-E-CONVERR, Energy convolution errror for SPEC0'
	GO TO 999

921	CONTINUE
	PRINT 200,'&NO-EMITTANCE-E-CONVERR, Energy convolution errror for SPEC1'
	GO TO 999

922	CONTINUE
	PRINT 200,'&NO-EMITTANCE-E-CONVERR, Energy convolution errror for SPEC2'
	GO TO 999

923	CONTINUE
	PRINT 200,'&NO-EMITTANCE-E-CONVERR, Energy convolution errror for SPEC3'
	GO TO 999

200	FORMAT(' ',8A)
210	FORMAT(' ',A,2(F10.3,A))
220	FORMAT(' ',A,I3,A,2(F10.3,A))

999	CONTINUE
	IERROR = -1
!-
!  -----------------------------------------------------------------------------
	RETURN
	END ! NO_EMITTANCE
!
	SUBROUTINE BRIGHTNESS_ARRAY(ICALC,I,R,ALPHAI,ALPHA2I,ALPHA)
 
!  Routine to pre-calculate the brightness array.

	use calc
	use angle_phi
	use angle_alp
	use spectra
	implicit none

!  Declarations of scalars:
	integer(kind=i4b) :: ICALC,I,IC,ID

	real(kind=dp) :: ALPHAI,ALPHA2I,R,ARG
	real(kind=dp) :: ALPHA2,S0,S1,S2,S3,AXR,AYR,AXI,AYI
	real(kind=dp) :: H,SINC

!  Declarations of arrays:
	real(kind=dp),dimension(nalpha) :: ALPHA

!
	IF (ICALC .EQ. 1) THEN ! Finite-N calculation
	    DO IC=1,NALPHA
		ALPHA2 = ALPHA(IC)*ALPHA(IC)
		ARG    = NPI*(ALPHA2-ALPHA2I)/R
		H      = SINC(ARG)
		DO ID=1,NPHI_BRIGHT
		    CALL BRIGHTE(I,ALPHA(IC),COSPHI(ID),SINPHI(ID),S0,S1,S2,S3,AXR,AYR,AXI,AYI)
		    BR0(ID,IC) = H*S0
		    BR1(ID,IC) = H*S1
		    BR2(ID,IC) = H*S2
		    BR3(ID,IC) = H*S3
		END DO ! ID
	    END DO ! IC
	ELSE ! Infinite-N approximation (H=1.0)
	    DO ID=1,NPHI_BRIGHT
	        CALL BRIGHTE(I,ALPHAI,COSPHI(ID),SINPHI(ID),S0,S1,S2,S3,AXR,AYR,AXI,AYI)
		BR0(ID,1) = S0
		BR1(ID,1) = S1
		BR2(ID,1) = S2
		BR3(ID,1) = S3
	    END DO ! ID
	END IF ! ICALC

	RETURN
	END ! BRIGHTNESS_ARRAY
!
	SUBROUTINE PDF_ARRAY(I,DXE,DYE,NXE,NYE)
 
!  Routine to pre-calculate the pdf array.

	use calc
	use power
	implicit none

!  Declarations of scalars:
	integer(kind=i4b) :: I,NXE,NYE
	integer(kind=i4b) :: IC,ID

	real(kind=dp) :: DXE,DYE
	real(kind=dp) :: ALPHAP,ALPHA2
	real(kind=dp) :: XG,YG,COSPHI,SINPHI
	real(kind=dp) :: S0,S1,S2,S3,AXR,AYR,AXI,AYI

!  Declarations of arrays:

!  Labeled constants:
	real(kind=dp),parameter :: ZERO=0.0D0,ONE=1.0D0,EPS=1.0D-6

!
        DO IC=0,NYE
	    YG = GAMMA_US*IC*DYE
	    DO ID=0,NXE
		XG     = GAMMA_US*ID*DXE
		ALPHA2 = XG*XG+YG*YG
		ALPHAP = SQRT(ALPHA2)
		IF (ALPHAP .LT. EPS) THEN
		    COSPHI = ZERO
		    SINPHI = ONE
		ELSE
		    COSPHI = XG/ALPHAP
		    SINPHI = YG/ALPHAP
		END IF ! ALPHAP
		CALL BRIGHTE(I,ALPHAP,COSPHI,SINPHI,S0,S1,S2,S3,AXR,AYR,AXI,AYI)
		PW0(ID,IC) = S0/(K3+ALPHA2)
		PW1(ID,IC) = S1/(K3+ALPHA2)
		PW2(ID,IC) = S2/(K3+ALPHA2)
		PW3(ID,IC) = S3/(K3+ALPHA2)

	    END DO ! ID
        END DO ! IC

	RETURN
	END ! PDF_ARRAY
!
	SUBROUTINE CONVOLUTE_DISTRIBUTION(ICALC,CONST,ALPHA,THETA,DALPHA, &
					  EPS,BR0,BR1,BR2,BR3, &
					  RA0,RA1,RA2,RA3,ICOUNT)
 
!  Routine to convolute the brightness arrray with the electron distribution.
 
	use us_size
	use beam
	use pinhole
	use angle_phi
	use angle_alp
	implicit none

!  Declarations of scalars:
	integer(kind=i4b) :: ICALC,ICOUNT,IA,IB,IC,ID

	real(kind=dp) :: CONST,DALPHA,EPS
	real(kind=dp) :: U,V,ARG,SL,P
	real(kind=dp) :: SUM0,SUM1,SUM2,SUM3
	real(kind=dp) :: DELTA0,DELTA1,DELTA2,DELTA3
 
!  Declarations of arrays:
	real(kind=dp),dimension(nalpha) :: ALPHA,THETA
	real(kind=dp),dimension(nphi_bright,nalpha) :: BR0,BR1,BR2,BR3
	real(kind=dp),dimension(nxp,nyp) :: RA0,RA1,RA2,RA3

!  Labeled constants:
	real(kind=dp),parameter :: ZERO=0.0D0
 
!
	SL = CONST*DALPHA*DPHI
	IF (ICALC .EQ. 1) THEN ! Finite-N calculation
	  DO IB=1,NYP
	    DO IA=1,NXP
	      SUM0 = ZERO
	      SUM1 = ZERO
	      SUM2 = ZERO
	      SUM3 = ZERO
	      DO IC=1,NALPHA
		DELTA0 = ZERO
		DELTA1 = ZERO
		DELTA2 = ZERO
		DELTA3 = ZERO
		DO ID=1,NPHI4
		  U   = CX(IA)-THETA(IC)*COSPHI(ID)
		  V   = CY(IB)-THETA(IC)*SINPHI(ID)
		  ARG = U*U*FU+V*V*FV
		  IF (ARG .LT. ARGMAX) THEN
		    P      = EXP(-ARG)
		    DELTA0 = DELTA0+P*BR0(INDEX_PHI(ID),IC)
		    DELTA1 = DELTA1+P*BR1(INDEX_PHI(ID),IC)
		    DELTA2 = DELTA2+P*BR2(INDEX_PHI(ID),IC)*S2SIGN(ID)
		    DELTA3 = DELTA3+P*BR3(INDEX_PHI(ID),IC)
		  END IF ! ARG
		END DO ! ID
		SUM0 = SUM0+DELTA0*ALPHA(IC)
		SUM1 = SUM1+DELTA1*ALPHA(IC)
		SUM2 = SUM2+DELTA2*ALPHA(IC)
		SUM3 = SUM3+DELTA3*ALPHA(IC)
	      END DO ! IC
	      SUM0       = SUM0*SL
	      SUM1       = SUM1*SL
	      SUM2       = SUM2*SL
	      SUM3       = SUM3*SL
	      RA0(IA,IB) = RA0(IA,IB)+SUM0
	      RA1(IA,IB) = RA1(IA,IB)+SUM1
	      RA2(IA,IB) = RA2(IA,IB)+SUM2
	      RA3(IA,IB) = RA3(IA,IB)+SUM3
	      IF (SUM0 .GT. EPS*RA0(IA,IB)) ICOUNT = 0
	    END DO ! IA
	  END DO ! IB
	ELSE ! Infinite-N approximation (no alpha integ.)
	  DO IB=1,NYP
	    DO IA=1,NXP
	      DELTA0 = ZERO
	      DELTA1 = ZERO
	      DELTA2 = ZERO
	      DELTA3 = ZERO
	      DO ID=1,NPHI4
		U   = CX(IA)-THETA(1)*COSPHI(ID)
		V   = CY(IB)-THETA(1)*SINPHI(ID)
		ARG = U*U*FU+V*V*FV
		IF (ARG .LT. ARGMAX) THEN
		  P      = EXP(-ARG)
	          DELTA0 = DELTA0+P*BR0(INDEX_PHI(ID),1)
	          DELTA1 = DELTA1+P*BR1(INDEX_PHI(ID),1)
	          DELTA2 = DELTA2+P*BR2(INDEX_PHI(ID),1)*S2SIGN(ID)
	          DELTA3 = DELTA3+P*BR3(INDEX_PHI(ID),1)
		END IF ! ARG
	      END DO ! ID
	      SUM0       = DELTA0*SL
	      SUM1       = DELTA1*SL
	      SUM2       = DELTA2*SL
	      SUM3       = DELTA3*SL
	      RA0(IA,IB) = RA0(IA,IB)+SUM0
	      RA1(IA,IB) = RA1(IA,IB)+SUM1
	      RA2(IA,IB) = RA2(IA,IB)+SUM2
	      RA3(IA,IB) = RA3(IA,IB)+SUM3
	      IF (SUM0 .GT. EPS*RA0(IA,IB)) ICOUNT = 0
	    END DO ! IA
	  END DO ! IB
	END IF ! ICALC
	
	RETURN
	END ! CONVOLUTE_DISTRIBUTION
!
	SUBROUTINE CONVOLUTE_POWER_DISTRIBUTION(CONST,DXE,DYE,NXE,NYE, &
					        EPS,PW0,PW1,PW2,PW3, &
					        RA0,RA1,RA2,RA3,ICOUNT)
 
!  Routine to convolute the pdf arrray with the electron distribution.
 
	use us_size
	use beam
	use pinhole
	implicit none

!  Declarations of scalars:
	integer(kind=i4b) :: NXE,NYE
	integer(kind=i4b) :: ICOUNT,IA,IB,IC,ID,ICA,IDA

	real(kind=dp) :: CONST,DXE,DYE,EPS
	real(kind=dp) :: WXE,WYE
	real(kind=dp) :: U,V,ARG,SL,P
	real(kind=dp) :: SUM0,SUM1,SUM2,SUM3
	real(kind=dp) :: DELTA0,DELTA1,DELTA2,DELTA3
 
!  Declarations of arrays:
	real(kind=dp),dimension(0:nxe,0:nye) :: PW0,PW1,PW2,PW3
	real(kind=dp),dimension(nxp,nyp) :: RA0,RA1,RA2,RA3

!  Labeled constants:
	real(kind=dp),parameter :: ZERO=0.0D0,ONE=1.0D0,HALF=0.5D0

!
	SL = CONST*DXE*DYE
	  DO IB=1,NYP
	    DO IA=1,NXP
	      SUM0 = ZERO
	      SUM1 = ZERO
	      SUM2 = ZERO
	      SUM3 = ZERO
	      DO IC=-NYE,NYE
	        WYE = ONE
	        IF (IC .EQ. -NYE .OR. IC .EQ. NYE) WYE = HALF
		DELTA0 = ZERO
		DELTA1 = ZERO
		DELTA2 = ZERO
		DELTA3 = ZERO
		DO ID=-NXE,NXE
	          WXE = ONE
	          IF (ID .EQ. -NXE .OR. ID .EQ. NXE) WXE = HALF
		  U   = CX(IA)-ID*DXE
		  V   = CY(IB)-IC*DYE
		  ARG = U*U*FU+V*V*FV
		  IF (ARG .LT. ARGMAX) THEN
		    P   = EXP(-ARG)
		    ICA = ABS(IC)
		    IDA = ABS(ID)
		    DELTA0 = DELTA0+WXE*WYE*P*PW0(IDA,ICA)
		    DELTA1 = DELTA1+WXE*WYE*P*PW1(IDA,ICA)
		    DELTA2 = DELTA2+WXE*WYE*P*PW2(IDA,ICA)
		    DELTA3 = DELTA3+WXE*WYE*P*PW3(IDA,ICA)
		  END IF ! ARG
		END DO ! ID
		SUM0 = SUM0+DELTA0
		SUM1 = SUM1+DELTA1
		SUM2 = SUM2+DELTA2
		SUM3 = SUM3+DELTA3
	      END DO ! IC
	      SUM0       = SUM0*SL
	      SUM1       = SUM1*SL
	      SUM2       = SUM2*SL
	      SUM3       = SUM3*SL
	      RA0(IA,IB) = RA0(IA,IB)+SUM0
	      RA1(IA,IB) = RA1(IA,IB)+SUM1
	      RA2(IA,IB) = RA2(IA,IB)+SUM2
	      RA3(IA,IB) = RA3(IA,IB)+SUM3
	      IF (SUM0 .GT. EPS*RA0(IA,IB)) ICOUNT = 0
	    END DO ! IA
	  END DO ! IB
	
	RETURN
	END ! CONVOLUTE_POWER_DISTRIBUTION
!
	SUBROUTINE CONVOLUTE_ENERGY_ESTEP(SP1)
 
!  Routine to convolve infinite-N-calculated spectra with the natural lineshape using Walker's method.

	use us_size
	use energy_mod
	use line_shape
	implicit none

!  Declarations of scalars:
	integer(kind=i4b) :: IE,IW
 
!  Declarations of arrays:
	real(kind=dp),dimension(*) :: SP1
	real(kind=dp),dimension(ne) :: SP2

!  Labeled constants:
	real(kind=dp),parameter :: ZERO=0.0D0
 
!
	DO IE=NE1,NE2
	    SP2(IE) = ZERO
	    DO IW=1,NW
		SP2(IE) = SP2(IE)+SP1(IE+NE1-IW)*HW(IW)
	    END DO ! IW
	END DO ! IE
!  Return in original array
	DO IE=NE1,NE2
	    SP1(IE) = SP2(IE)*DW
	END DO ! IE

	RETURN
	END ! CONVOLUTE_ENERGY_ESTEP
!
	SUBROUTINE CONVOLUTE_ENERGY_VSTEP(SP1)
 
!  Routine to convolve infinite-N-calculated spectra with the natural lineshape using Dejus' method with variable energy step size.

	use us_size
	use calc
	use energy_mod
	use step
	implicit none

!  Declarations of scalars:
	integer(kind=i4b) :: IE,IEU,J1,J2

	real(kind=dp) :: EUP,EU1,EU2
	real(kind=dp) :: ARG,P,H1,H2,S1,S2,SMV
	real(kind=dp) :: SINC
 
!  Declarations of arrays:
	real(kind=dp),dimension(*) :: SP1
	real(kind=dp),dimension(neu) :: SP2

!  Labeled constants:
	real(kind=dp),parameter :: ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0
 
!
	J1 = 1
	J2 = 1
	DO IEU=1,NEU
	    EUP = EU(IEU)
	    EU1 = EUP-EW
	    EU2 = EUP+EW
	    CALL HUNT(E,NE,EU1,J1)
	    CALL HUNT(E,NE,EU2,J2)
	    IF (J2 .EQ. 0 .OR. J1 .EQ. NE) THEN
		SP2(IEU) = ZERO

	    ELSE IF (J1 .EQ. 0) THEN
		IF (J2 .LT. NE) THEN
!  Interpolate at upper limit
		    H2  = EU2-E(J2)
		    P   = H2/HE(J2)
		    S2  = (ONE-P)*SP1(J2)+P*SP1(J2+1)
		    ARG = PE*(EUP-EU2)
		    SMV = S2*SINC(ARG)*H2/TWO
		    ARG = PE*(EUP-E(J2))

! RJD added check 10/23/2014
!		    SMV = SMV+SP1(J2)*SINC(ARG)*(H2+HE(J2-1)/TWO)
!		    DO IE=1,J2-1
!			ARG = PE*(EUP-E(IE))
!			SMV = SMV+SP1(IE)*SINC(ARG)*HA(IE)
!		    END DO ! IE
		    IF (J2 .EQ. 1) THEN
		        SMV = SMV+SP1(J2)*SINC(ARG)*H2/TWO
		    ELSE ! J2 > 1
			SMV = SMV+SP1(J2)*SINC(ARG)*(H2+HE(J2-1)/TWO)
		        DO IE=1,J2-1
			    ARG = PE*(EUP-E(IE))
			    SMV = SMV+SP1(IE)*SINC(ARG)*HA(IE)
		        END DO ! IE
		    ENDIF
		    SP2(IEU) = SMV/DEW
		ELSE ! J2 = NE
		    SMV = ZERO		    
		    DO IE=1,NE
			ARG = PE*(EUP-E(IE))
			SMV = SMV+SP1(IE)*SINC(ARG)*HA(IE)
		    END DO ! IE
		    SP2(IEU) = SMV/DEW
		END IF ! J2

	    ELSE IF (J2 .LT. NE) THEN
		IF (J1 .EQ. J2) THEN
		    SP2(IEU) = ZERO
		ELSE
!  Interpolate at lower limit
		    H1  = E(J1+1)-EU1
		    P   = H1/HE(J1)
		    S1  = P*SP1(J1)+(ONE-P)*SP1(J1+1)
		    ARG = PE*(EUP-EU1)
		    SMV = S1*SINC(ARG)*H1/TWO
		    ARG = PE*(EUP-E(J1+1))
		    SMV = SMV+SP1(J1+1)*SINC(ARG)*(H1+HE(J1+1))/TWO
!  Interpolate at upper limit
		    H2  = EU2-E(J2)
		    P   = H2/HE(J2)
		    S2  = (ONE-P)*SP1(J2)+P*SP1(J2+1)
		    ARG = PE*(EUP-EU2)
		    SMV = SMV+S2*SINC(ARG)*H2/TWO
		    ARG = PE*(EUP-E(J2))
		    SMV = SMV+SP1(J2)*SINC(ARG)*(H2+HE(J2-1)/TWO)
		    DO IE=J1+2,J2-1
			ARG = PE*(EUP-E(IE))
			SMV = SMV+SP1(IE)*SINC(ARG)*HA(IE)
		    END DO ! IE
		    SP2(IEU) = SMV/DEW
		END IF

	    ELSE ! J2 = NE
!  Interpolate at lower limit
		H1  = E(J1+1)-EU1
		P   = H1/HE(J1)
		S1  = P*SP1(J1)+(ONE-P)*SP1(J1+1)
		ARG = PE*(EUP-EU1)
		SMV = S1*SINC(ARG)*H1/TWO
		ARG = PE*(EUP-E(J1+1))
		SMV = SMV+SP1(J1+1)*SINC(ARG)*(H1+HE(J1+1))/TWO
		DO IE=J1+2,NE
		    ARG = PE*(EUP-E(IE))
		    SMV = SMV+SP1(IE)*SINC(ARG)*HA(IE)
		END DO ! IE
		SP2(IEU) = SMV/DEW
	    END IF
	END DO ! IEU
		
!  Return in original array
	DO IEU=1,NEU
	    SP1(IEU) = SP2(IEU)
	END DO ! IEU

	RETURN
	END ! CONVOLUTE_ENERGY_VSTEP
!
	FUNCTION SINC(ARG)
  
!  Routine to calculate the SINC function.

	use precision_standard
	implicit none

!  Declarations of scalars:
	real(kind=dp) :: SINC
	real(kind=dp) :: ARG,SL
 
!  Labeled constants:
	real(kind=dp),parameter :: ONE=1.0D0,EPS=1.0D-6
	 
	IF (ABS(ARG) .GT. EPS) THEN
	    SL   = SIN(ARG)/ARG
	    SINC = SL*SL
	ELSE
	    SINC = ONE
	END IF ! ARG

	RETURN
	END ! SINC
!
	SUBROUTINE TRAPZ2(RA,AREA)
 
!  Routine for trapedzoidal integration.

	use us_size
	use pinhole
	implicit none

!  Declarations of scalars:
	integer(kind=i4b) :: IA,IB
	real(kind=dp) :: AREA,SMV,WX,WY
 
!  Declarations of arrays:
	real(kind=dp),dimension(nxp,nyp) :: RA

!  Labeled constants:
	real(kind=dp),parameter :: ZERO=0.0D0,ONE=1.0D0,HALF=0.5D0

!
	SMV = ZERO
	DO IB=1,NYP
	    IF (IB .EQ. 1 .OR. IB .EQ. NYP) THEN
		WY = HALF
	    ELSE
		WY = ONE
	    END IF ! IB
	    DO IA=1,NXP
		IF (IA .EQ. 1 .OR. IA .EQ. NXP) THEN
		    WX = HALF
		ELSE
		    WX = ONE
		END IF ! IA
	    SMV = SMV+WX*WY*RA(IA,IB)
	    END DO ! IA
	END DO ! IB
	AREA = SMV*DXP*DYP

	RETURN
	END ! TRAPZ2
!
	SUBROUTINE PRINT_OUT(ISUB)
 
!  Routine for all printouts.

	use prt
	use calc
	use factor
	use pinhole
	use energy_mod
	use spectra
	use physical_constants !  Fundamental physical constants; Physics Today Aug. 1990:
	use pi_constants
	implicit none

!  Declarations of scalars:
	integer(kind=i4b) :: ISUB,IA,IB,IE,IEU
	integer(kind=i4b) :: IMIN,IMAX
	real(kind=dp) :: WE,POWER,FLUX
	real(kind=dp) :: P1,P2,P3,P4

!  Declarations of arrays:

!  Conversion factors:
	real(kind=dp),parameter :: C_MA_A=1.0D-3,C_CM_M=1.0D-2

!  Labeled constants:
	real(kind=dp),parameter :: ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,HALF=0.5D0

!+ Print
!  -----------------------------------------------------------------------------
	FLUX  = ZERO
	POWER = ZERO
!  Spatial-distribution (non-zero emittance and zero emittance case)
	IF ((ISUB .EQ. 1) .OR. (ISUB .EQ. 5 .AND. MODE .EQ. 1)) THEN
	    IF (ISUB .EQ. 1) THEN ! non-zero emittance
                WRITE (2,200)     ! blank to indicate non-zero emittance
	    ELSE
	        WRITE (2,200) '=== Zero emittance ==='
	    END IF ! ISUB
	        
	    CALL TRAPZ2(RA0,FLUX) ! get integrated flux over observation area
	    FLUX  = FAC*FLUX      ! ph/s/0.1%bw
	    POWER = FLUX/BW*EC    ! W/eV
	    IF (LANG) THEN
	        WRITE (2,232) 'Angular flux density distribution for ', EMINU,' eV (ph/s/mr^2/0.1%bw):'
	    ELSE
	        WRITE (2,232) 'Irradiance for ',EMINU,' eV  @ ',D,' m', ' (ph/s/mm^2/0.1%bw):'
	    END IF ! LANG
	    WRITE (2,250) 'Integrated flux  ',FLUX, ' ph/s/0.1%bw'
            WRITE (2,252) 'Integrated power ',POWER,' W/eV'
	    WRITE (2,256) 'Contributing harmonics ',I1(1),' to ',I2(1)
	    WRITE (2,200)
	    IF (LANG) THEN
	        WRITE (2,200) '   x(mr)      y(mr)     Ang. flux density', '   p1      p2      p3      p4'
	    ELSE
	        WRITE (2,200) '   x(mm)      y(mm)     Irradiance', '          p1      p2      p3      p4'
	    END IF ! LANG
	       
	    DO IA=1,NXP
		DO IB=1,NYP
		    IF (RA0(IA,IB) .GT. ZERO) THEN
			P1 = RA1(IA,IB)/RA0(IA,IB)	
			P2 = RA2(IA,IB)/RA0(IA,IB)	
			P3 = RA3(IA,IB)/RA0(IA,IB)	
			P4 = ONE -SQRT(P1*P1 +P2*P2 +P3*P3) ! unpolarized
		    ELSE
			P1 = ZERO
			P2 = ZERO
			P3 = ZERO
			P4 = ZERO
		    END IF 
		    WRITE (2,210) XP(IA),YP(IB),RA0(IA,IB),P1,P2,P3,P4
		END DO ! IB
	    END DO ! IA

	ELSE IF ((ISUB .EQ. 2) .OR. (ISUB .EQ. 3) .OR. (ISUB .EQ. 5)) THEN
!  Spectral distributions (non-zero emittance and zero emittance case)
	    IF (METHOD .NE. 3) THEN ! non-zero emittance
                WRITE (2,200)       ! blank to indicate non-zero emittance
	    ELSE
	        WRITE (2,200) '=== Zero emittance ==='
	    END IF ! ISUB
!           Get frequency-integrated power (flux) (no meaning for brilliance)
	    IMIN = 10000
	    IMAX = 0
	    DO IE=NE1,NE2
		IF (I1(IE) .LT. IMIN .AND. I1(IE) .GT. 0) IMIN = I1(IE)
		IF (I2(IE) .GT. IMAX) IMAX = I2(IE)
	    END DO ! IE
	    IF (METHOD .EQ. 4) THEN ! Dejus' method
		DO IEU=1,NEU
		    IF (IEU .EQ. 1 .OR. IEU .EQ. NEU) THEN
			WE = HALF
		    ELSE
			WE = ONE
		    END IF ! IEU
		    FLUX  = FLUX +WE*SPEC0(IEU)/EU(IEU) ! -> ps/s/mm^2  /eV or
		                                        ! -> ph/s/mrad^2/eV
		    POWER = POWER+WE*SPEC0(IEU) ! -> W/mm^2/eV or W/mrad^2/eV
		END DO ! IEU
	    ELSE ! METHOD={1,2,3,14}
		DO IE=NE1,NE2
		    IF (IE .EQ. NE1 .OR. IE .EQ. NE2) THEN
			WE = HALF
		    ELSE
			WE = ONE
		    END IF ! IE
		    FLUX  = FLUX +WE*SPEC0(IE)/E(IE) ! -> ps/s/mm^2  /eV or
		                                     ! -> ph/s/mrad^2/eV
		    POWER = POWER+WE*SPEC0(IE) ! -> W/mm^2/eV or W/mrad^2/eV
		END DO ! IE
	    END IF ! METHOD
	    FLUX  = FLUX /BW   *DE ! ph/s/mm^2 or ph/s/mrad^2
	    POWER = POWER/BW*EC*DE ! W or W/mm^2 or W/mrad^2

	    IF (MODE .EQ. 2) THEN ! angular/spatial flux density spectrum
	        IF (LANG) THEN
	            WRITE (2,231) 'Angular flux density spectrum for ', &
	                          XP(1),' mr,',YP(1),' mr', &
	                         ' (ph/s/mr^2/0.1%bw):'
	            WRITE (2,250) 'Integrated flux  density ',FLUX, &
	                         ' ph/s/mr^2'
                    WRITE (2,254) 'Integrated power density ',POWER, &
			          ' W/mr^2'
	    	    WRITE (2,256) 'Contributing harmonics ',IMIN, &
				  ' to ',IMAX
	            WRITE (2,200)
	            WRITE (2,200) 'Energy(eV)   Ang. flux density', &
		                  '   p1      p2      p3      p4'
	        ELSE
	            WRITE (2,230) 'Irradiance for ',XP(1),' mm,',YP(1), &
	                         ' mm  @ ',D,' m', &
	                         ' (ph/s/mm^2/0.1%bw):'
	            WRITE (2,250) 'Integrated flux  density ',FLUX, &
	                         ' ph/s/mm^2'
                    WRITE (2,254) 'Integrated power density ',POWER, &
	                         ' W/mm^2'
	    	    WRITE (2,256) 'Contributing harmonics ',IMIN, &
				  ' to ',IMAX
	            WRITE (2,200)
	            WRITE (2,200) 'Energy(eV)   Irradiance', &
		              '          p1      p2      p3      p4'
	        END IF ! LANG
	    ELSE IF (MODE .EQ. 3) THEN ! brilliance
                WRITE (2,200) 'On-axis brilliance', &
	                     ' (ph/s/mr^2/mm^2/0.1%bw):'
	        WRITE (2,200)
	        WRITE (2,200)
	    	WRITE (2,256) 'Contributing harmonics ',IMIN, &
				  ' to ',IMAX
	        WRITE (2,200)
	        WRITE (2,200) 'Energy(eV)   Brilliance'
	    ELSE IF (MODE .EQ. 4) THEN ! flux spectrum through a pinhole
	        IF (LANG) THEN
	            WRITE (2,234) 'Flux through ',XPS,' mr x',YPS, &
	                         ' mr pinhole at',XPC,' mr,',YPC,' mr:'
	        ELSE
	            WRITE (2,234) 'Flux through ',XPS,' mm x',YPS, &
	                         ' mm pinhole at',XPC,' mm,',YPC,' mm @', &
	                         D,' m:'
	        END IF ! LANG
                WRITE (2,250) 'Integrated flux  ',FLUX, ' ph/s'
                WRITE (2,254) 'Integrated power ',POWER,' W'
	    	WRITE (2,256) 'Contributing harmonics ',IMIN, &
				  ' to ',IMAX
	        WRITE (2,200)
	        WRITE (2,200) 'Energy(eV)   Flux(ph/s/0.1%bw)', &
		              '   p1      p2      p3      p4'
	    ELSE IF (MODE .EQ. 5) THEN ! angle-integrated spectrum
                WRITE (2,200) 'Angle-integrated spectrum:'
	        WRITE (2,250) 'Integrated flux  ',FLUX, ' ph/s'
                WRITE (2,254) 'Integrated power ',POWER,' W'
	    	WRITE (2,256) 'Contributing harmonics ',IMIN, &
				  ' to ',IMAX
	        WRITE (2,200)
	        WRITE (2,200) 'Energy(eV)   Flux(ph/s/0.1%bw)', &
		              '   p1      p2      p3      p4'
	    END IF ! MODE

	    IF (METHOD .EQ. 4) THEN ! Dejus'
                DO IEU=1,NEU
		    IF (SPEC0(IEU) .GT. ZERO) THEN
			P1 = SPEC1(IEU)/SPEC0(IEU)
			P2 = SPEC2(IEU)/SPEC0(IEU)
			P3 = SPEC3(IEU)/SPEC0(IEU)
			P4 = ONE -SQRT(P1*P1 +P2*P2 +P3*P3) ! unpolarized
		    ELSE
			P1 = ZERO
			P2 = ZERO
			P3 = ZERO
			P4 = ZERO
		    END IF 
		    IF (MODE .EQ. 3) THEN ! Brilliance
		    	WRITE (2,220) EU(IEU),SPEC0(IEU)
		    ELSE
		    	WRITE (2,220) EU(IEU),SPEC0(IEU),P1,P2,P3,P4
		    END IF
		END DO ! IEU
	    ELSE ! METHOD={1,2,3,14}
		DO IE=NE1,NE2
		    IF (SPEC0(IE) .GT. ZERO) THEN
			P1 = SPEC1(IE)/SPEC0(IE)
			P2 = SPEC2(IE)/SPEC0(IE)
			P3 = SPEC3(IE)/SPEC0(IE)
			P4 = ONE -SQRT(P1*P1 +P2*P2 +P3*P3) ! unpolarized
		    ELSE
			P1 = ZERO
			P2 = ZERO
			P3 = ZERO
			P4 = ZERO
		    END IF
		    IF (MODE .EQ. 3) THEN ! Brilliance
		    	WRITE (2,220) E(IE),SPEC0(IE)
		    ELSE
		    	WRITE (2,220) E(IE),SPEC0(IE),P1,P2,P3,P4
		    END IF
		END DO ! IE
	    END IF ! METHOD

	ELSE IF (ISUB .EQ. 4) THEN
!  Power density (non-zero emittance and zero emittance case)
	    IF (METHOD .NE. 3) THEN ! non-zero emittance
                WRITE (2,200)       ! blank to indicate non-zero emittance
	    ELSE
	        WRITE (2,200) '=== Zero emittance ==='
	    END IF ! ISUB
	    CALL TRAPZ2(RA0,POWER)
	    POWER = FAC*POWER ! W

	    IF (LANG) THEN
	    	WRITE (2,200) 'Power density distribution (W/mr^2):'
	    ELSE
	    	WRITE (2,234) 'Power density distribution @ ', &
			       D,' m (W/mm^2):'
	    END IF ! LANG
	    WRITE (2,200) ! blank to indicate flux
            WRITE (2,254) 'Integrated power ',POWER,' W'
	    WRITE (2,256) 'Contributing harmonics ',I1(1),' to ',I2(1)
	    WRITE (2,200)
	    IF (LANG) THEN
	        WRITE (2,200) '   x(mr)      y(mr)     Power density', &
		              '       p1      p2      p3      p4'
	    ELSE
	        WRITE (2,200) '   x(mm)      y(mm)     Power density', &
		              '       p1      p2      p3      p4'
	    END IF ! LANG
	       
	    DO IA=1,NXP
		DO IB=1,NYP
		    IF (RA0(IA,IB) .GT. ZERO) THEN
			P1 = RA1(IA,IB)/RA0(IA,IB)	
			P2 = RA2(IA,IB)/RA0(IA,IB)	
			P3 = RA3(IA,IB)/RA0(IA,IB)	
			P4 = ONE -SQRT(P1*P1 +P2*P2 +P3*P3) ! unpolarized
		    ELSE
			P1 = ZERO
			P2 = ZERO
			P3 = ZERO
			P4 = ZERO
		    END IF 
		    WRITE (2,210) XP(IA),YP(IB),RA0(IA,IB),P1,P2,P3,P4
		END DO ! IB
	    END DO ! IA
	END IF ! ISUB
	CLOSE (UNIT=2)
	
200	FORMAT(' ',8A)
210	FORMAT(' ',F9.4,TR2,F9.4,TR3,1PE13.6,TR2,0P,4(TR2,F6.3))
220	FORMAT(' ',F10.3,TR2,1PE13.6,TR2,0P,4(TR2,F6.3))
222	FORMAT(' ',F10.3,1(TR2,1PE13.6),2I5)

230	FORMAT(' ',(A,F7.3,A),2(F7.3,A),A)
231	FORMAT(' ',(A,F7.3,A),1(F7.3,A),A)
232	FORMAT(' ',(A,F7.1,A),(F7.3,A),A)
234	FORMAT(' ',(A,F7.3,A),4(F7.3,A))
250	FORMAT(' ',(A,1PE10.3,A))
252	FORMAT(' ',(A,F10.6,A))
254	FORMAT(' ',(A,F10.1,A))
256	FORMAT(' ',2(A,I4))
!230	FORMAT(' ',I5,TR2,F10.3,TR2,D13.6,TR2,D13.6,TR2,I4,TR2,I4)
!235	FORMAT(' ',I5,TR2,F10.3,TR2,D13.6,TR2,I4,TR2,I4)
!250	FORMAT(' ',F10.3,TR2,1PE13.6)
!260	FORMAT(' ',1P4E13.5)
!265	FORMAT(' ',5F15.5)
!270	FORMAT(' ',A,F15.3,A)
!275	FORMAT(' ',2F15.3)
	RETURN
	END ! PRINT_OUT
!
	SUBROUTINE CHECK_FILE(INPUT)

!  Routine to check if a file exists and quit if it does not exist.

	character(len=180) :: INPUT
	logical(kind=4) :: FILE_EXISTS

	INQUIRE(FILE=INPUT, EXIST=FILE_EXISTS)

	IF (FILE_EXISTS .NEQV. .TRUE.) THEN
	  PRINT 200, '&CHECK_FILE-F-NOFILE, File not found ... ',TRIM(INPUT)
	CALL EXIT(0)
	ENDIf

200	FORMAT(' ',8A)
	END ! CHECK_FILE
