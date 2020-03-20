	program ws
C+
C PROGRAM DESCRIPTION:	
C  Program to calculate wiggler and bending magnet spectra using the
C  Bessel function approximation. The program may be executed from the xop
C  interface.
C 
C AUTHORS: 
C  Roger J. Dejus
C  The Advanced Photon Source
C  Experimental Facilities Division
C  Argonne National Laboratory
C 
C CREATION DATE: 
C  17-FEB-1994
C 
C INPUT PARAMETERS:
C  The input parameters are divided into sections related to the storage ring,
C  the wiggler device, and the quantity to be calculated.
C Machine Parameters:
C  Storage ring energy 			(GeV)
C  Storage ring current			(mA)
C Wiggler Parameters:
C  Period length			(cm)
C  Number of periods
C Note: For a bending magnet source: set N=0.5, and make Ky large and adjust
C       the period length accordingly. For example, put Ky=9.34 and calculate
C       the period length from, Period (cm) = 10.0/B0(T), where B0 is the known
C       strength of the magnetic field (in Tesla) for the bending magnet.  The
C       calculated power density (pd) is correct, but the total power (ptot)
C       is irrelevant. Typically make the extend of the pinhole small in the
C	horizontal direction (theta << Ky/gamma) as the intensity should
C	not depend on the horizontal offset. Check value of B0 (and critical
C	energy EC0) in the plot file.
C  Deflection parameter (hor.  field) Kx (= 0.0 only; for elliptical wiggler
C  not yet implemented)
C  Deflection parameter (vert. field) Ky
C Scan Parameters:
C  Minimum energy			(eV)
C  Maximum energy			(eV)
C  Number of energy points
C Pinhole Parameters:
C  Distance from the source		(m)
C    (d=0.0 => angular units)
C  X-coordinate for center of pinhole	(mm) or (mrad)
C  Y-coordinate for center of pinhole	(mm) or (mrad)
C  X-size of pinhole (full width)	(mm) or (mrad)
C  Y-size of pinhole (full width)	(mm) or (mrad)
C    (for angular units (d=0.0) values are entered in mrad)
C    (X is for horizontal direction)
C    (Y is for the vertical direction)
C  Number of subdivisions of pinhole in X (max 200)
C  Number of subdivisions of pinhole in Y (max 200)
C
C Mode:
C  Depending on the mode selected, some of the pinhole parameters may be
C  set to different values by the program; see the output file ws.plt.
C  MODE    1    Angular/spatial flux density distribution
C  MODE    2    Angular/spatial flux density spectrum
C  MODE    3    On-axis brilliance spectrum (not yet implemented)
C  MODE    4    Flux spectrum through a pinhole
C  MODE    5    Flux spectrum integrated over all angles
C  MODE    6    Power density and integrated power
C
C  Angular/spatial flux density distribution
C    - Flux distribution at the energy chosen as minimum energy.
C  Angular/spatial flux density spectrum
C    - Spectrum at any given point in space as selected by the X and Y
C      coordinate for the center of the pinhole. X is horizontal and Y is
C      vertical.
C  On-axis brilliance spectrum (not yet implemented)
C  Flux spectrum through a pinhole
C    - Spectrum through a pinhole centered at X-center and Y-center with
C      size X-size and Y-size.  The energy range is from the minimum to the
C      maximum energy.
C  Flux spectrum integrated over all angles (wiggler only).
C    -  The pinhole parameters have no significance here.
C  Power density and integrated power
C    -  Integrated over all energies, thus the energy parameters have no
C       significance here.
C
C Polarization:
C  The normalized Stokes parameters are calculated including the 
C  unpolarized component.
C
C DESIGN ISSUES:
C  Program calculates the spectra from the Modified Bessel functions.  See K.J.
C  Kim, in "Physics of Particle Accelerators", vol. 1, AIP Conference Proc. 184
C  Ed. R.G. Lerner, New York (1989), p. 583, Eq. (3.12).
C  The algorithm is based on a series expansion for small arguments Z
C  (Abramowitz & Stegun Eq. 9.6.2 and 9.6.10) and an asymptotic expansion for
C  large arguments (Eq. 9.7.2).
C  Reference: Handbook of Mathematical Functions, Eds. Milton Abramowitz and
C  Irene A. Stegun, Ninth Printing, Dover Publications, New York (1970).
C  NOTE: THE POLARIZATION PARAMETERS ARE PROVIDED ALTHOUGH NOT THOROUGHLY
C  TESTED - USE WITH CAUTION.
C  
C COPYRIGHT:
C  This routine must only be used at The Advanced Photon Source and must not
C  be tranferred or used at any other location without the written consent
C  of the author.
C  
C FILES USED:
C  Input file - ws.dat (tc.inp for XOP)  File in the user's current directory
C                       containing the input parameters.
C  Output file - ws.plt (tc.out for XOP) File in the user's current directory
C                       containing the results of the calculation. The header contains
C                       all input parameters and the calculated on-axis first
C			harmonic energy (e1), corresponding wavelength (l1), 
C			total power (ptot), and the on-axis power density (pd).
C			See note above when using N=0.5 for bending magnet.
C KEYWORDS:
C  Wiggler Spectrum, Modified Bessel Function of Second kind.
C  
C LINK/LIBRARY ISSUES:
C  The gamma function is needed.  Currently uses no library routines. The values
C  for gamma(2/3) and gamma(1/3) are stored as constants. May be substituted by
C  calls to the NAG library routine S14AAF.
C  
C PORTABILITY ISSUES:
C  Runs on DEC 3000/400 AXP alpha (Tru64Unix v5.0)
C  Release v5.6), and Windows 95/98/NT (Pentium and higher).
C  SUN: f90: Sun Fortran 95 8.3 SunOS_sparc 2007/05/03 (Tested November 24, 2008)
C  Updated October 8, 2013 (Argonne National Laboratory), tested on:
C  *** Linux Red Hat Enterprise Linux Workstation release 6.3 (Santiago) ***
C  Red Hat Enterprise Linux (RHEL) 64-bit with the Intel(R) Fortran
C  Intel(R) 64 Compiler XE for applications running on Intel(R) 64,
C  Version 13.1.1.163 Build 2013031, and with GFORTRAN, gcc version 4.4.6 20120305
C  (Red Hat 4.4.6-4) (GCC).
C  *** Sun Solaris SunOS 5.10 Generic_147440-27 sun4u sparc SUNW,Sun-Blade-2500 ***
C  Sun Fortran 90/95 8.4 SunOS_sparc Patch 128231-02 2009/10/20 with the -f77 switch.
C  and with GFORTRAN, gcc version 4.5.1 (GCC).
C  Windows 7/8 64-bit and MacOS X 10.6 (and newer) are also supported.
C  The GFORTRAN compiler (GCC) v4.8.1 was used for compilations on Windows and (GCC) v4.6.1 on MacOS.
C  
C TIMING:
C  Execution time is typically very fast but depends on the quantity being
C  calculated. Typically seconds to at the most minutes.
C
C EXAMPLES:
C Ex. 1 using the input file ~/test/ws.txt (output file becomes ws.plt in the current working directory)
C % ws ~/test/ws.txt
C Ex. 2 using the default input file ws.dat in the current working directory (the output file becomes ws.plt).
C % ws
C Ex. 3 using the input abc in the current working directory (the output file becomes abc.plt).
C % ws abc 
C
C VERSION:
C  1.62
C  
C MODIFICATION HISTORY:
C 
C	 Date     | Name  | Description
C ----------------+-------+-----------------------------------------------------
C 17-JUL-2000     | RJD   | Adopted from v1.4 which was never released to the
C                 |       | public. Turned off circular polarization for the
C                 |       | wiggler (which is only valid for EMW). 
C                 |       | Current version is v1.5.
C ----------------+-------+-----------------------------------------------------
C 24-NOV-2008     | RJD   | Modified format statements 263 and 264 in main program
C                 |       | to print additional digits. Modified format statements
C                 |       | 210 and 220 in the print_out routine to print three
C                 |       | digits in the exponent to allow for very small (and
C                 |       | large) values. With the newest Sun Compiler, floating
C                 |       | point underflow is being flagged (but not trapped)
C                 |       | on the standard output. Nothing "extra" is written 
C                 |       | to the output file "ws.plt."
C                 |       | The smallest value (double precision): 2.2250739e-308
C                 |       | The largest value (double precision): 1.7976931e+308
C                 |       | The floating underflow comes from the variables 
C                 |       | asigma and api which are less than sqrt(smallest value).
C                 |       | Values are simply put to zero and no harm is done.
C                 |       | Changed format statement 254 to print one more decimal.
C                 |       | Increased the number of subdivisions of the pinhole 
C                 |       | from 50 to 200 (P_SZ=201).
C                 |       | Current version is v1.6.
C ----------------+-------+-----------------------------------------------------
C 02-NOV-2013     | RJD   | Updated calls to date and time routines.
C                 |       | Current version is v1.61.
C ----------------+-------+-----------------------------------------------------
C 10-JUN-2014     | RJD   | Updated so that an arbitrary input file can be used on the command line.
C                 |       | If no input file is given on the command line then the file 'ws.dat'
C                 |       | is assumed ('ws.inp' for the XOP version). The output filename is created
C                 |       | from the rootname, which is derived from the input filename using the string after
C                 |       | the last directory separator (/) without its trailing file extension (if it exists).
C                 |       | The output filename is the rootname with the extension (.plt) appended (.out for the
C                 |       | XOP version). Search "standalone" for changing defaults.
C                 |       | Current version is v1.62.
C ----------------+-------+-----------------------------------------------------
C [change_entry]
C-
C
C  Size parameters:
	integer*4	E_SZ,P_SZ
	parameter	(E_SZ=10000,P_SZ=201)
  
C  Declaration of scalars:
	character*80	title
	character*8     dbuff1
	character*10    dbuff2
	character*10    tbuff1
	character*8     tbuff2
	character*80	file_name
	character*180	datafile,plotfile,rootname

	integer*4	ierror
	integer*4	ie,l,ia,ib,isub
	integer*4	nch1
	integer*4       sd,ed

	real*4		dtime,delta

	real*8		energy,cur,period
	real*8		du
	real*8		g2,lamdar,lamda1,e1z
	real*8		d2
	real*8		xpmin,ypmin
	real*8		ptot,pd,gk
	real*8		k2

C  Declarations of arrays:
	real*4		td(2)
  
C  Fundamental physical constants; Physics Today Aug. 1990:
	real*8		C,ME,MEE,EC,H,HBAR,MUZ,EPSZ
	parameter	(C    =2.99792458D8)	! Speed of light [m/s]
	parameter	(ME   =9.1093897D-31)	! Electron rest mass [kg]
	parameter	(MEE  =0.51099906D0)	! Electron rest mass [MeV]
	parameter	(EC   =1.60217733D-19)	! Elementary charge [C]
	parameter	(H    =6.6260755D-34)	! Planck's constant [Js]
	parameter	(HBAR =1.05457266D-34)	! Planck's constant/2Pi [Js]
	parameter	(MUZ  =1.2566370614D-6)	! Permeability of vacuum [NA-2]
	parameter	(EPSZ =8.854187817D-12)	! Permittivity of vacuum [Fm-1]

C  Conversion factors:
	real*8		C_EVANG,C_MM_M,C_M_MM,C_CM_ANG,C_MRAD_RAD,C_RAD_MRAD
	real*8		C_MA_A,C_CM_M
	parameter	(C_EVANG=H*C/EC*1.0D10,C_MM_M=1.0D-3,C_M_MM=1.0D3)
	parameter	(C_CM_ANG=1.0D8,C_MRAD_RAD=1.0D-3,C_RAD_MRAD=1.0D3)
	parameter	(C_MA_A=1.0D-3,C_CM_M=1.0D-2)

C  Filename defaults (standalone and XOP)
	character*(*)  	FT_OUT
c	parameter     	(FT_OUT='.plt')			! default output filename extension
	parameter     	(FT_OUT='.out')			! default output filename extension (XOP)
	character*(*)  	FILE_NAME_IN
c	parameter     	(FILE_NAME_IN = 'ws.dat') 	! default input filename in the current directory
	parameter     	(FILE_NAME_IN = 'ws.inp') 	! default input filename in the current directory (XOP)

C  Labeled constants:
	character*(*)   VN
	parameter       (VN='v. 1.62')
	real*8		PI,PIHALF,TWOPI,SQRT3
	parameter	(PI    =3.1415 92653 58979 32384 62643D0)
	parameter	(PIHALF=1.5707 96326 79489 66192 31322D0)
	parameter	(TWOPI= 6.2831 85307 17958 64769 25287D0)
	parameter	(SQRT3= 1.7320 50807 56887 72935 27446D0)
	real*8		FINE_STRUCTURE_CONST
	parameter	(FINE_STRUCTURE_CONST=1.0D19*EC*1.0D19*EC
	1		/(4.0D0*PI*EPSZ*1.0D38*HBAR*C))
	real*8          PTOT_FAC,PD_FAC
	parameter       (PTOT_FAC=PI/3.0D0*EC/EPSZ/(MEE**2)*1.0D+6)   ! 0.07257
	parameter       (PD_FAC  =21.0D0/(16.0D0*PI*MEE**2)*PTOT_FAC) ! 0.11611

	real*8		BW
	parameter	(BW=1.0D-3) ! 0.1%
	real*8		ZERO,ONE,TWO,HALF,EPSK
	parameter	(ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,HALF=0.5D0)
	parameter	(EPSK=1.0D-4)

	real*8		CL1,CL2,CL3,CL4
	parameter	(CL1=EC/(TWO*PI*ME*C)*C_CM_M)       ! 0.934...
	parameter	(CL2=1.5D0/(MEE*MEE)*HBAR/ME*1.0D6) ! 665.0 ...
	parameter	(CL3=3.0D0*FINE_STRUCTURE_CONST/
	1		     (4.0D0*PI*PI))                 ! 5.545...10^-4
	parameter	(CL4=SQRT3*FINE_STRUCTURE_CONST/
	1		     (2.0D0*PI))                    ! 2.012...10^-3

C  Logical name for parameter file
	character*(*)	LN
	parameter 	(LN='wspar:') ! VMS
C
C  Common blocks:
	logical*4	lang
	logical*4	lbm
	integer*4	mode
	integer*4	nxp,nyp
	integer*4	ne
	real*8		d
	real*8		emin,emax,xpc,ypc,xps,yps
	real*8		k,kx,ky,k3,gamma,er,len,n,b0,ec0,psi0
	real*8		fac,facs,faca,facp,facp1
	real*8		dxp,dyp
	real*8		de
	real*8		e(E_SZ)
	real*8		spec0(E_SZ),spec1(E_SZ),spec2(E_SZ),spec3(E_SZ)
	real*8		ra0(P_SZ,P_SZ),ra1(P_SZ,P_SZ),
	1		ra2(P_SZ,P_SZ),ra3(P_SZ,P_SZ)
	real*8		xp(P_SZ),yp(P_SZ),cx(P_SZ),cy(P_SZ)

	common		/prti/      lang
	common		/prtr/      d,emin,emax,xpc,ypc,xps,yps
	common		/calci/     mode,lbm
	common		/calcr/	    k,kx,ky,k3,gamma,er,len,n,b0,ec0,psi0
	common		/factor/    fac,facs,faca,facp,facp1
	common		/pinhole/   xp,yp,cx,cy,nxp,nyp,dxp,dyp
	common		/energyi/   ne
	common		/energyr/   e,spec0,spec1,spec2,spec3,de
	common		/spectra/   ra0,ra1,ra2,ra3
C
C+ First executable statement here
C  -----------------------------------------------------------------------------
c2345678901234567890123456789012345678901234567890123456789012345678901234567890
	delta = dtime(td)	! start clock
	call date_and_time(Date=dbuff1,Time=tbuff1)
	tbuff2 = tbuff1(1:2)//':'//tbuff1(3:4)//':'//tbuff1(5:6)
	dbuff2 = dbuff1(1:4)//'-'//dbuff1(5:6)//'-'//dbuff1(7:8)

c  check number of command line arguments (not used here)
c	if (command_argument_count() .ne. 1) then
c	   print *, "Usage: tc data_file"
c	   call exit(0)	   
c	endif

c  extract the input filename to be used as the rootname
	call get_command_argument(1,datafile)
c	if (trim(adjustl(datafile(1:1))) .eq. '') datafile = FILE_NAME_IN	! default for undefined variable (adjustl not necessary
	if (trim(datafile(1:1)) .eq. '') datafile = FILE_NAME_IN	! default for undefined variable
	sd = scan(datafile, '/', back=.true.)
	ed = scan(datafile, '.', back=.true.)
	rootname = datafile(sd +1:ed -1)
	if (trim(rootname(1:1)) .eq. '') rootname = datafile	! default for undefined variable
	call check_file(datafile)
	plotfile = trim(rootname)//FT_OUT

c	print *,'datafile:',datafile
c	print *,'rootname:',rootname
c	print *,'plotfile:',plotfile
c	call exit(0)

C+ Open plot file
C  -----------------------------------------------------------------------------
 	open (unit=2,file=trim(plotfile),status='unknown')

	print 200,'****************** Start WS ',VN,
	1         ' at ',dbuff2,' ',tbuff2,' ******************'
	print *
	write (2,200) '****************** Start WS ',VN,
	1         ' at ',dbuff2,' ',tbuff2,' ******************'
	write (2,200)

C+ Read parameters from data file
C  -----------------------------------------------------------------------------
C	file_name = ln//file_name_in ! VMS
C	file_name = file_name_in
	open (unit=1,file=trim(datafile),status='old',access='sequential',
	1     form='formatted',action='read')

	read (1,100) title
	read (1,*)   energy,cur		      ! [GeV],[mA]
	read (1,*)   period,n,kx,ky	      ! [cm],[dl],[dl],[dl]
	read (1,*)   emin,emax,ne	      ! [eV],[eV],[dl]
	read (1,*)   du,xpc,ypc,xps,yps,nxp,nyp![m],[mm],[mm],[mm],[mm],[dl],[dl]
	read (1,*)   mode
	close (unit=1)

C+ Print mode selection
C  -----------------------------------------------------------------------------
	if (mode .eq. 1) then
	  print 200,'Angular/spatial flux density distribution'
	else if (mode .eq. 2) then
	  print 200,'Angular/spatial flux density spectrum'
	else if (mode .eq. 3) then
	  print 200,'On-axis brightness spectrum'
	else if (mode .eq. 4) then
	  print 200,'Flux spectrum through a pinhole'
	else if (mode .eq. 5) then
	  print 200,'Flux spectrum integrated over all angles'
	else if (mode .eq. 6) then
	  print 200,'Power density and integrated power'
	end if ! mode

C+ Check for valid input
C  -----------------------------------------------------------------------------
	if (xpc .lt. ZERO .or. ypc .lt. ZERO) go to 900
	if (ne .gt. E_SZ) go to 920

C+ Remove options currently not implemented
C  -----------------------------------------------------------------------------
	if (mode .eq. 3) go to 910
	if (kx .ne. ZERO) go to 915

	if (n .eq. HALF) then	! bending magnet
	  lbm = .true.
	else
	  lbm = .false.		! wiggler
	endif

	if (lbm .and. mode .eq. 5) go to 916

	l = 80
	do while(title(l:l) .eq. ' ') ! strip blanks
	  l = l-1
	enddo ! while
	nch1 = l

C  Determine units (angular or spatial units)
	if (du .eq. ZERO) then
	  lang = .true.
	  d    = ONE
	else
	  lang = .false.
	  d    = du
	end if !

	if (mode .eq. 2) then
	  xps = ZERO
	  yps = ZERO
	  nxp = 0
	  nyp = 0
	else if (mode .eq. 3 .or. mode .eq. 5) then
	  xpc = ZERO
	  ypc = ZERO
	  xps = ZERO
	  yps = ZERO
	  nxp = 0
	  nyp = 0
	end if ! mode
C-
C  -----------------------------------------------------------------------------

C+ Definition of constants
C  -----------------------------------------------------------------------------
	gamma  = energy/MEE*1.0d3
	g2     = gamma*gamma
	lamdar = period*C_CM_ANG/(TWO*g2) ! Reduced wavelength [A]
	er     = C_EVANG/lamdar           ! Reduced energy    [eV]
	k2     = kx*kx +ky*ky
	k3     = ONE+k2/TWO
	lamda1 = lamdar*k3                ! First harmonic on axis  [A]
	e1z    = er    /k3                ! First harmonic on axis [eV]
	d2     = d*d                      ! Distance squared     [m**2]
	len    = n*period*C_CM_M          ! Length of device        [m]

	gk     = ZERO
	if (kx .lt. EPSK .or. ky .lt. EPSK) then
	  k  = kx +ky
	  gk = k*((k**6)+(24.0d0*(k**4)/7.0d0)+
	1      (4.0d0*k*k)+(16.0d0/7.0d0))/((1.0d0+(k*k))**3.5d0)
	end if
	if (abs(kx-ky) .lt. EPSK) then
	  k  = kx
	  gk = 32.0d0/7.0d0*k/((1.0d0+(k*k))**3.0d0)
	end if

	ptot = PTOT_FAC*n*k2*  (energy**2)*cur*C_MA_A/(period*C_CM_M)![W]
	pd   = PD_FAC  *n*k*gk*(energy**4)*cur*C_MA_A/(period*C_CM_M)![W/mrad^2]
	b0   = ky/(CL1*period)
	ec0  = CL2*energy*energy*b0 ! [eV]
	psi0 = kx/gamma

	facs = CL3*g2*BW*cur/EC*C_MA_A/(C_RAD_MRAD*C_RAD_MRAD)*
	1      TWO*n/d2                         !fac spectral_distribution
	faca = CL4*TWO*k*BW*cur/EC*C_MA_A*TWO*n !fac angle-integrated spectrum
						!first TWO is for hor. symmetry
	facp = pd/d2			 	! W/mm^2
	facp1= facs/BW*EC                       ! conversion to W/mm^2 or W/mrad^2
C  -----------------------------------------------------------------------------

C+ Write input parameters and derived quantities
C  -----------------------------------------------------------------------------
	write (2,200) title(1:nch1)
	write (2,250) 'energy ',energy,' GeV','current',cur ,' mA'
	write (2,255) 'period ',period,'  cm','N      ',n   ,'   ',
	1             'Kx     ',kx    ,'    ','Ky     ',ky
	write (2,256) 'emin   ',emin ,'  eV','emax ',emax,' eV',
	1             'ne     ',ne

	if (lang) then
	  write (2,261) 'd      ',du    ,'   m','xpc   ',xpc ,' mr',
	1                 'ypc   ',ypc   ,' mr'
	  write (2,257) 'xps    ',xps   ,'  mr','yps   ',yps ,' mr',
	1		'nxp    ',nxp   ,'   ','nyp    ',nyp
	else
	  write (2,261) 'd      ',du    ,'   m','xpc   ',xpc ,' mm',
	1                 'ypc   ',ypc   ,' mm'
	  write (2,257) 'xps    ',xps   ,'  mm','yps   ',yps ,' mm',
	1		'nxp    ',nxp   ,'   ','nyp    ',nyp
	end if ! lang

	write (2,200)
	write (2,257) 'b0     ',b0   ,'   T','ec0   ',ec0*1.0d-3,
	1	      '  keV'
	write (2,263) 'e1   ',e1z   ,'  eV','l1 ',lamda1,'  A'

	if (lang) then
	  write (2,264) 'ptot ',ptot  ,'   W','pd ',pd,
	1		'  W/mr^2'
	else
	  write (2,264) 'ptot ',ptot  ,'   W','pd ',pd/d2,
	1		'  W/mm^2'
	endif ! lang
	write (2,200)

C+ Define energy scale
C  -----------------------------------------------------------------------------
	de = (emax-emin)/(ne-1)       ! [eV]
	do ie=1,ne
	  e(ie) = emin +(ie-1)*de
	end do ! ie

C+ Pinhole parameters
C  -----------------------------------------------------------------------------
	if (xpc .eq. ZERO .and. ypc .eq. ZERO) then ! Pinhole centered
	  fac   = 4.0d0
	  xpmin = ZERO
	  ypmin = ZERO
	  if (nxp .gt. 0) dxp = xps/TWO/nxp
	  if (nyp .gt. 0) dyp = yps/TWO/nyp
	else
	  fac   = ONE
	  xpmin = xpc-xps/TWO
	  ypmin = ypc-yps/TWO
	  if (nxp .gt. 0) dxp = xps/nxp
	  if (nyp .gt. 0) dyp = yps/nyp
	end if
	nxp = nxp+1
	nyp = nyp+1
C-
C  -----------------------------------------------------------------------------

C+ Set up positions within pinhole and associated cartesian angles
C  -----------------------------------------------------------------------------
	do ia=1,nxp
	  xp(ia) = xpmin+(ia-1)*dxp ! Position [mm]
	  cx(ia) = xp(ia)*C_MM_M/d  ! Angle   [rad]
	end do ! ia
	do ib=1,nyp
	  yp(ib) = ypmin+(ib-1)*dyp ! Position [mm]
	  cy(ib) = yp(ib)*C_MM_M/d  ! Angle   [rad]
	end do ! ib
C-
C  -----------------------------------------------------------------------------

C+ Call analysis routine
C  -----------------------------------------------------------------------------
	if (mode .eq. 1) then
	  isub = 1
	else if (mode .ge. 2 .and. mode .le. 4) then
	  isub = 2
	else if (mode .eq. 5) then
	  isub = 3
	else if (mode .eq. 6 .and. .not. lbm) then
	  isub = 4
	else if (mode .eq. 6 .and. lbm) then
	  isub = 5
	end if ! mode

	if (isub .eq. 1) then
	  call space_distribution(ierror)
	else if (isub .eq. 2) then
	  call spectral_distribution(ierror)
	else if (isub .eq. 3) then
	  call angle_integration(ierror)
	else if (isub .eq. 4) then
	  call power_distribution(ierror)
	else if (isub .eq. 5) then
	  call power_distribution1(ierror)
	end if ! isub
C-
C  -----------------------------------------------------------------------------

C+ Print results
C  -----------------------------------------------------------------------------
	if (ierror .eq. 0) call print_out(isub)
C-
C  -----------------------------------------------------------------------------

C+ Exit
C  -----------------------------------------------------------------------------
	if (ierror .eq. 0) then
	  delta = dtime(td) ! read clock
	  print 210,'Elapsed time: ',delta,' s'
	  print 200,'&WS-S-NORMAL, Successful completion'
	  print *
	else
	  print 200,'&WS-F-SUBERR, Subroutine error'
	  print 200,'- unsuccessful completion due to error status'
	  print *
	  write (2,*)
	  write (2,200) '&WS-F-SUBERR, Subroutine error'
	  write (2,200) '- unsuccessful completion due to error status'
	  write (2,*)
	end if ! ierror	

	close (unit=2)
	call exit(0)

100	format(a)

200	format(' ',8a)
201	format('1')
205	format(' ',a,i6,a,i6)
207	format(' ',a,i2,a)
208	format(' ',a,f8.3,a)
210	format(' ',a,2(f10.3,a))
220	format(' ',f10.3,2(tr2,1pe13.6))
230	format(' ',(a,f7.3,a),2(f7.3,a))
250	format(' ',(a,f7.3,a,tr3),(a,f6.1,a))
255	format(' ',(a,f7.3,a,tr3),(a,f6.1,a,tr3),2(a,f6.3,a,tr2))
256	format(' ',(a,f7.1,a,tr3),(a,f8.1,a,tr3),(a,i6))
257	format(' ',(a,f7.3,a,tr3),(a,f7.3,a,tr3),2(a,i6,a,tr3))
261	format(' ',3(a,f7.3,a,tr3))
c263	format(' ',(a,f7.1,a,tr3),(a,f8.3,a,tr3))
c264	format(' ',(a,f7.1,a,tr3),(a,f10.1,a,tr3))
263	format(' ',(a,f9.2,a,tr3),(a,f10.3,a,tr3))
264	format(' ',(a,f9.1,a,tr3),(a,f10.1,a,tr3))
270	format(' ',(a,1pe10.3,a))

C+ Error returns
C  -----------------------------------------------------------------------------
900	continue
	print 200,'&WS-E-INVDAT, Invalid data'
	print 200,'- check input data file; X and Y must be ',
	1	  'in the first quadrant.'
	go to 999

910	continue
	print 200,'&WS-E-INVDAT, Invalid data'
	print 207,'- check input data file; mode ',mode,
	1	  ' is not a valid mode.'
	go to 999

915	continue
	print 200,'&WS-E-INVDAT, Invalid data'
	print 208,'- check input data file; kx ',kx,
	1	  ' is not valid; only zero allowed.'
	go to 999

916	continue
	print 200,'&WS-E-INVDAT, Invalid data'
	print 207,'- check input data file; mode ',mode,
	1	  ' is not a valid mode for bending magnet.'
	go to 999

920	continue
	print 200,'&WS-F-BNDERR, Boundary error'
	print 205,'- energy array out of bounds; number of points ',ne,
	1	  ' is greater than ',E_SZ
	go to 999

999	continue
	print *
	print 200,'&WS-F-PRGERR, Program error'
	print 200,'- unsuccessful completion due to error status'
	print *
C-
C  -----------------------------------------------------------------------------
	close (unit=2)
	call exit(0)

	end ! ws
C
	subroutine space_distribution(ierror)

C[subroutine_header_comments]

C  Size parameters:
	integer*4	E_SZ,P_SZ
	parameter	(E_SZ=10000,P_SZ=201)

C  Declarations of scalars:
	integer*4	ierror
	integer*4	ia,ib

	real*8		xg,yg,yg1,yg2,y,y2,eta,asigma,api
	real*8		ecp
	real*8		cc,cpi,cpis
	real*8		k13,k23
	real*8		vp,vpmax,fc

C  Declarations of arrays:

C  Labeled constants:
	real*8		ZERO,ONE,TWO,HALF
	parameter	(ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,HALF=0.5D0)

c  ECP_FRAC determines the lowest critical energy used.
c  Ecmin=Ec0*sqrt(ECP_FRAC)
	real*8		ECP_FRAC
	parameter	(ECP_FRAC=1.0D-6)

C  Common blocks:
	integer*4	mode
	logical*4	lbm
	integer*4	nxp,nyp
	integer*4	ne
	real*8		k,kx,ky,k3,gamma,er,len,n,b0,ec0,psi0
	real*8		fac,facs,faca,facp,facp1
	real*8		dxp,dyp
	real*8		de
	real*8		e(E_SZ)
	real*8		spec0(E_SZ),spec1(E_SZ),spec2(E_SZ),spec3(E_SZ)
	real*8		ra0(P_SZ,P_SZ),ra1(P_SZ,P_SZ),
	1		ra2(P_SZ,P_SZ),ra3(P_SZ,P_SZ)
	real*8		xp(P_SZ),yp(P_SZ),cx(P_SZ),cy(P_SZ)

	common		/calci/     mode,lbm
	common		/calcr/	    k,kx,ky,k3,gamma,er,len,n,b0,ec0,psi0
	common		/factor/    fac,facs,faca,facp,facp1
	common		/pinhole/   xp,yp,cx,cy,nxp,nyp,dxp,dyp
	common		/energyi/   ne
	common		/energyr/   e,spec0,spec1,spec2,spec3,de
	common		/spectra/   ra0,ra1,ra2,ra3
C
	ierror = 0
	vpmax  = sqrt(ONE -ECP_FRAC)

c  Calculate distribution at e = emin (eV)
C  Unit is ph/s/mr^2/0.1%bw for angular flux density and ph/s/mm^2/0.1%bw for
C  spatial flux density
	do ib=1,nyp
	  do ia=1,nxp
	    ra0(ia,ib) = ZERO
	    ra1(ia,ib) = ZERO
	    ra3(ia,ib) = ZERO

	    xg  = gamma*cx(ia) ! gamma*theta
	    yg  = gamma*cy(ib) ! gamma*psi
	    vp  = abs(xg)/ky
	    if (vp .gt. vpmax) goto 810
	    ecp = ec0*sqrt(ONE-vp*vp)  ! eV

	    yg2 = yg*yg
	    yg1 = ONE +yg2
	    cc  = yg1*yg1
	    cpi = yg2/yg1
	    if (yg .gt. ZERO) then ! keep track of sign
	      cpis = +sqrt(cpi)
	    else
	      cpis = -sqrt(cpi)
	    endif
	    fc  = facs*cc
	    y   = e(1)/ecp
	    y2  = y*y
	    eta = HALF*y*yg1**1.5d0

	    asigma = k23(eta) ! Modified Bessel function of 2nd kind (2/3)
	    api    = k13(eta) ! Modified Bessel function of 2nd kind (1/3)
	    ra0(ia,ib) = y2*fc*(asigma*asigma +cpi*api*api)
	    ra1(ia,ib) = y2*fc*(asigma*asigma -cpi*api*api)
	    if (lbm) ra3(ia,ib) = TWO*y2*fc*(asigma*cpis*api)
810	    continue
	    end do ! ia
	end do ! ib

	return
	end ! space_distribution
C
	subroutine spectral_distribution(ierror)
 
C[subroutine_header_comments]

C  Size parameters:
	integer*4	E_SZ,P_SZ
	parameter	(E_SZ=10000,P_SZ=201)

C  Declarations of scalars:
	integer*4	ierror
	integer*4	ia,ib,ie

	real*8		xg,yg,yg1,yg2,y,y2,eta,asigma,api
	real*8		ecp
	real*8		cc,cpi,cpis,area0,area1,area3
	real*8		k13,k23
	real*8		vp,vpmax,fc

C  Declarations of arrays:

C  Labeled constants:
	real*8		ZERO,ONE,TWO,HALF
	parameter	(ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,HALF=0.5D0)

c  ECP_FRAC determines the lowest critical energy used.
c  Ecmin=Ec0*sqrt(ECP_FRAC)
	real*8		ECP_FRAC
	parameter	(ECP_FRAC=1.0D-6)

C  Common blocks:
	logical*4	lbm
	integer*4	mode
	integer*4	nxp,nyp
	integer*4	ne
	real*8		k,kx,ky,k3,gamma,er,len,n,b0,ec0,psi0
	real*8		fac,facs,faca,facp,facp1
	real*8		dxp,dyp
	real*8		de
	real*8		e(E_SZ)
	real*8		spec0(E_SZ),spec1(E_SZ),spec2(E_SZ),spec3(E_SZ)
	real*8		ra0(P_SZ,P_SZ),ra1(P_SZ,P_SZ),
	1		ra2(P_SZ,P_SZ),ra3(P_SZ,P_SZ)
	real*8		xp(P_SZ),yp(P_SZ),cx(P_SZ),cy(P_SZ)

	common		/calci/     mode,lbm
	common		/calcr/	    k,kx,ky,k3,gamma,er,len,n,b0,ec0,psi0
	common		/factor/    fac,facs,faca,facp,facp1
	common		/pinhole/   xp,yp,cx,cy,nxp,nyp,dxp,dyp
	common		/energyi/   ne
	common		/energyr/   e,spec0,spec1,spec2,spec3,de
	common		/spectra/   ra0,ra1,ra2,ra3
C
	ierror = 0
	vpmax  = sqrt(ONE -ECP_FRAC)

C+ Loop over energies
C  -----------------------------------------------------------------------------
C  Unit is ph/s/mr^2/0.1%bw for angular flux density and ph/s/mm^2/0.1%bw for
C  spatial flux density
	do ie=1,ne
	  do ib=1,nyp
	    do ia=1,nxp
	      ra0(ia,ib) = ZERO
	      ra1(ia,ib) = ZERO
	      ra3(ia,ib) = ZERO

	      xg  = gamma*cx(ia) ! gamma*theta
	      yg  = gamma*cy(ib) ! gamma*psi
	      vp  = abs(xg)/ky
	      if (vp .gt. vpmax) goto 810
	      ecp = ec0*sqrt(ONE-vp*vp)  ! eV

	      yg2 = yg*yg
	      yg1 = ONE +yg2
	      cc  = yg1*yg1
	      cpi = yg2/yg1
	      if (yg .gt. ZERO) then ! keep track of sign
		cpis = +sqrt(cpi)
	      else
		cpis = -sqrt(cpi)
	      endif
	      fc  = facs*cc
	      y   = e(ie)/ecp
	      y2  = y*y
	      eta = HALF*y*yg1**1.5d0

	      asigma = k23(eta) ! Modified Bessel function of 2nd kind (2/3)
	      api    = k13(eta) ! Modified Bessel function of 2nd kind (1/3)
	      ra0(ia,ib) = y2*fc*(asigma*asigma +cpi*api*api)
	      ra1(ia,ib) = y2*fc*(asigma*asigma -cpi*api*api)
	      if (lbm) ra3(ia,ib) = TWO*y2*fc*(asigma*cpis*api)
810	      continue
	      end do ! ia
	  end do ! ib

	  if (mode .eq. 4) then  ! Integrate over pinhole
	    call trapz2(ra0,area0)
	    call trapz2(ra1,area1)
	    call trapz2(ra3,area3)
	    spec0(ie) = fac*area0
	    spec1(ie) = fac*area1
	    spec3(ie) = fac*area3
	  else if (mode .eq. 2 .or. mode .eq. 3) then ! Fixed position
	    spec0(ie) = ra0(1,1)
	    spec1(ie) = ra1(1,1)
	    spec3(ie) = ra3(1,1)
	  end if ! mode (else mode = 1)

	end do ! ie
C- Endloop energy
C  -----------------------------------------------------------------------------

	return
	end ! spectral_distribution
C
	subroutine angle_integration(ierror)

C[subroutine_header_comments]

C  Size parameters:
	integer*4	E_SZ,P_SZ
	parameter	(E_SZ=10000,P_SZ=201)
	integer*4	NXA_SZ
	parameter	(NXA_SZ=100)

C  Declarations of scalars:
	integer*4	ierror
	integer*4	ie,ia,nxa,iopt

	real*8		y,sumx,ecp,ecn,ec02,dec,vn,vp,dvn,dvs
	real*8		ecn2,ecp2,dec2
	real*8		gy

C  Declarations of arrays:
	real*8		ecpa(NXA_SZ),wgt(NXA_SZ)

C  Fundamental physical constants; Physics Today Aug. 1990:
	real*8		C,ME,MEE,EC,H,HBAR,MUZ,EPSZ
	parameter	(C    =2.99792458D8)	! Speed of light [m/s]
	parameter	(ME   =9.1093897D-31)	! Electron rest mass [kg]
	parameter	(MEE  =0.51099906D0)	! Electron rest mass [MeV]
	parameter	(EC   =1.60217733D-19)	! Elementary charge [C]
	parameter	(H    =6.6260755D-34)	! Planck's constant [Js]
	parameter	(HBAR =1.05457266D-34)	! Planck's constant/2Pi [Js]
	parameter	(MUZ  =1.2566370614D-6)	! Permeability of vacuum [NA-2]
	parameter	(EPSZ =8.854187817D-12)	! Permittivity of vacuum [Fm-1]

C  Conversion factors:
	real*8		C_EVANG,C_ANG_M,C_M2_MM2
	parameter	(C_EVANG=H*C/EC*1.0D10,C_ANG_M=1.0D-10) ! 12398.42
	parameter	(C_M2_MM2=1.0D6)

C  Labeled constants:
	real*8		PI,PIHALF,TWOPI
	parameter	(PI    =3.1415 92653 58979 32384 62643D0)
	parameter	(PIHALF=1.5707 96326 79489 66192 31322D0)
	parameter	(TWOPI= 6.2831 85307 17958 64769 25287D0)
	real*8		ZERO,ONE,TWO,HALF
	parameter	(ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,HALF=0.5D0)

c  ECx_FRAC determines the lowest critical energy used for the hor. int.
c  Ec1min=Ec0*EC1_FRAC, Ec2min=Ec0*sqrt(EC2_FRAC), Ec3min=Ec0*sqrt(EC3_FRAC)
	real*8		EC1_FRAC,EC2_FRAC,EC3_FRAC,Y_LIM
	parameter	(EC1_FRAC=1.0D-3,EC2_FRAC=1.0D-6,EC3_FRAC=1.0D-4)
	parameter	(Y_LIM=40.0D0)

C  Common blocks:
	integer*4	ne
	real*8		k,kx,ky,k3,gamma,er,len,n,b0,ec0,psi0
	real*8		fac,facs,faca,facp,facp1
	real*8		de
	real*8		e(E_SZ)
	real*8		spec0(E_SZ),spec1(E_SZ),spec2(E_SZ),spec3(E_SZ)

	common		/calcr/	    k,kx,ky,k3,gamma,er,len,n,b0,ec0,psi0
	common		/factor/    fac,facs,faca,facp,facp1
	common		/energyi/   ne
	common		/energyr/   e,spec0,spec1,spec2,spec3,de
C
	ierror = 0
	nxa = NXA_SZ ! Number of steps for integration over horizontal angle

c  iopt = 1 will emphasize accuracy for energies below ~ Ec/3, and iopt = 3
c  for energies above ~ 3xEc, and iopt = 2 around Ec. In all cases, the
c  accuracy is better than 0.5 %.
	iopt = 2

c  Options define the stepsizes (weights) and the critical energies
C  stored in arrays wgt and ecpa
	go to (1000,2000,3000) iopt

c  Option #1: Constant stepsize in critical energy
1000	continue
	vp   = ZERO
	dvs  = ZERO
	ec02 = ec0*ec0
	ecp  = ec0
	dec  = ec0*(ONE -EC1_FRAC)/(nxa-1)
	ecpa(1) = ec0
	do ia=1,nxa
	  if (ia .lt. nxa) then
	    ecn = ecp -dec
	    vn  = sqrt(ONE -ecn*ecn/ec02)
	    dvn = vn -vp
	    ecp = ecn
	    vp  = vn
	    ecpa(ia+1) = ecp
	  else
	    dvn = ZERO
	  end if
	  wgt(ia) = HALF*(dvs +dvn)
	  dvs  = dvn
	end do ! ia
	go to 890

c  Option #2: Constant stepsize in the square of the critical energy
2000	continue
	vp   = ZERO
	dvs  = ZERO
	ec02 = ec0*ec0
	ecp2 = ec02
	dec2 = ec02*(ONE -EC2_FRAC)/(nxa-1)
	ecpa(1) = ec0
	do ia=1,nxa
	  if (ia .lt. nxa) then
	    ecn2 = ecp2 -dec2
	    vn   = sqrt(ONE -ecn2/ec02)
	    dvn  = vn -vp
	    ecp2 = ecn2
	    vp   = vn
	    ecpa(ia+1) = sqrt(ecp2)
	  else
	    dvn = ZERO
	  end if
	  wgt(ia) = HALF*(dvs +dvn)
	  dvs  = dvn
	end do ! ia
	go to 890

c  Option #3: Constant stepsize in horizontal angle
3000	continue
	dec = sqrt(ONE -EC3_FRAC)/(nxa-1)
	do ia=1,nxa
	  if (ia .eq. 1 .or. ia .eq. nxa) then
	    wgt(ia) = HALF*dec
	  else
	    wgt(ia) = ONE*dec
	  end if
	  vp  = (ia-1)*dec
	  ecp = ec0*sqrt(ONE -vp*vp)
	  ecpa(ia) = ecp
	end do ! ia
	go to 890

c  Integrate
890	continue
	do ie=1,ne
	  sumx = ZERO
	  do ia=1,nxa
	    y = e(ie)/ecpa(ia)
	    if (y .gt. Y_LIM) go to 810
	    sumx = sumx +wgt(ia)*gy(1,y)
	  end do ! ia
810	  continue

	  spec0(ie) = faca*sumx
	  spec3(ie) = ZERO
	end do ! ie

	return
	end ! angle_integration
C
	subroutine power_distribution(ierror)
c  Routine for calculation of power density. The routine uses the integral
c  equation given by K.J. Kim in "Angular Distribution of Undulator Power 
c  for an Arbitrary Deflection Parameter K", Nucl. Instr. Meth. A246, 
c  (1986) 67-70, Eq. 5, for K < KMAX (typically 100.0).
c  For K > KMAX, Eq. 10 in the same paper is used which gives the power
c  density in the limit K -> infinity. In this limit, there is no
c  radiation beyond gamma*theta/k in the horizontal plane.

C  Size parameters:
	integer*4	E_SZ,P_SZ
	parameter	(E_SZ=10000,P_SZ=201)

C  Declarations of scalars:
	integer*4	ierror,ia,ib
	real*8		xg,yg,s0,s1,s2,s3

C  Declarations of arrays:

C  Labeled constants:
	real*8		ZERO,ONE,HALF
	parameter	(ZERO=0.0D0,ONE=1.0D0,HALF=0.5D0)
	real*8		KMAX
   	parameter	(KMAX=100.0D0)
c	parameter	(KMAX=1.0D0)

C  Common blocks:
	integer*4	nxp,nyp
	real*8		k,kx,ky,k3,gamma,er,len,n,b0,ec0,psi0
	real*8		fac,facs,faca,facp,facp1
	real*8		dxp,dyp
	real*8		ra0(P_SZ,P_SZ),ra1(P_SZ,P_SZ),
	1		ra2(P_SZ,P_SZ),ra3(P_SZ,P_SZ)
	real*8		xp(P_SZ),yp(P_SZ),cx(P_SZ),cy(P_SZ)

	common		/calcr/	    k,kx,ky,k3,gamma,er,len,n,b0,ec0,psi0
	common		/factor/    fac,facs,faca,facp,facp1
	common		/pinhole/   xp,yp,cx,cy,nxp,nyp,dxp,dyp
	common		/spectra/   ra0,ra1,ra2,ra3

	ierror = 0
	do ib=1,nyp
	  do ia=1,nxp
	    xg = gamma*cx(ia) ! gamma*theta
	    yg = gamma*cy(ib) ! gamma*psi
	    if (k .lt. KMAX) then 
	      call fk(xg,yg,k,s0,s1,s2,s3)
	    else
	      call fkl(xg,yg,k,s0,s1,s2,s3) ! formula for large-K limit
	    endif
	    ra0(ia,ib) = facp*s0
	    ra1(ia,ib) = facp*s1
	    ra2(ia,ib) = facp*s2
	    ra3(ia,ib) = facp*s3
	  enddo
	enddo

	return
	end ! power_distribution
C
	subroutine power_distribution1(ierror)

c  Routine for calculation of power density. Integrates the Bessel functions over
c  the energy.  The right circular intensity appears above the horizontal plane 
c  (gamma*psi > 0.0) using the current definition. Circular polarized intensity
c  is only valid for the bending magnet and the elliptical multipole wiggler.

C  Size parameters:
	integer*4	E_SZ,P_SZ
	parameter	(E_SZ=10000,P_SZ=201)
c  NETA_SZ was determined such that the power density and integrated power
c  gives good agreement with the results from routine fkl in subroutine
c  power_distribution.  Depends also on the integration range determined
c  by ETA_MIN and ETA_MAX.
	integer*4	NETA_SZ
	parameter	(NETA_SZ=500)

C  Declarations of scalars:
	integer*4	ierror,ia,ib,ie,neta
	real*8		xg,yg,yg1,yg2,y,y2
	real*8		ecp
	real*8		cc,cpi,cpis,c1,fc
	real*8		k13,k23
	real*8		vp,vpmax,deta,etamin,etamax

C  Declarations of arrays:
	real*8		eta(NETA_SZ),we(NETA_SZ)
	real*8		asigma(NETA_SZ),api(NETA_SZ)

C  Labeled constants:
	real*8		ZERO,ONE,TWO,HALF
	parameter	(ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,HALF=0.5D0)

c  ECP_FRAC determines the lowest critical energy used.
c  Ecmin=Ec0*sqrt(ECP_FRAC)
c  ETA_MIN and ETA_MAX determines the integration range for the Bessel
c  functions.
	real*8		ECP_FRAC,ETA_MIN,ETA_MAX
	parameter	(ECP_FRAC=1.0D-6,ETA_MIN=1.0D-4,ETA_MAX=15.0D0)

C  Common blocks:
	integer*4	nxp,nyp
	real*8		k,kx,ky,k3,gamma,er,len,n,b0,ec0,psi0
	real*8		fac,facs,faca,facp,facp1
	real*8		dxp,dyp
	real*8		ra0(P_SZ,P_SZ),ra1(P_SZ,P_SZ),
	1		ra2(P_SZ,P_SZ),ra3(P_SZ,P_SZ)
	real*8		xp(P_SZ),yp(P_SZ),cx(P_SZ),cy(P_SZ)

	common		/calcr/	    k,kx,ky,k3,gamma,er,len,n,b0,ec0,psi0
	common		/factor/    fac,facs,faca,facp,facp1
	common		/pinhole/   xp,yp,cx,cy,nxp,nyp,dxp,dyp
	common		/spectra/   ra0,ra1,ra2,ra3

	ierror = 0
	neta   = NETA_SZ
	etamin = ETA_MIN
	etamax = ETA_MAX
	deta   = (etamax -etamin)/(neta-1)
	vpmax  = sqrt(ONE -ECP_FRAC)
	do ie=1,neta ! pre-calculate eta-array and Bessel functions
	  if (ie .eq. 1 .or. ie .eq. neta) then
	    we(ie) = HALF
	  else
	    we(ie) = ONE
	  endif
	  eta   (ie) = etamin +(ie-1)*deta
	  asigma(ie) = k23(eta(ie)) ! Modified Bessel function of 2nd kind (2/3)
	  api   (ie) = k13(eta(ie)) ! Modified Bessel function of 2nd kind (1/3)
	enddo
	do ib=1,nyp
	  do ia=1,nxp
	    ra0(ia,ib) = ZERO
	    ra1(ia,ib) = ZERO
	    ra2(ia,ib) = ZERO
	    ra3(ia,ib) = ZERO
	    xg = gamma*cx(ia) ! gamma*theta
	    yg = gamma*cy(ib) ! gamma*psi
	    vp = abs(xg)/ky
	    if (vp .gt. vpmax) goto 810  ! no contribution beyond gamma*theta > Ky
	    ecp = ec0*sqrt(ONE -vp**2) ! eV

	    yg2 = yg*yg
	    yg1 = ONE +yg2
	    cc  = yg1*yg1
	    cpi = yg2/yg1
	    if (yg .gt. ZERO) then ! keep track of sign
	      cpis = +sqrt(cpi)
	    else
	      cpis = -sqrt(cpi)
	    endif
	    c1  = TWO/yg1**1.5
	    fc  = facp1*cc*deta*c1*ecp
	    do ie=1,neta
	      y  = c1*eta(ie)
	      y2 = y*y
	      ra0(ia,ib) = ra0(ia,ib)
	1	  +y2*fc*we(ie)*(asigma(ie)*asigma(ie) 
	2	  +cpi*api(ie)*api(ie))
	      ra1(ia,ib) = ra1(ia,ib)
	1         +y2*fc*we(ie)*(asigma(ie)*asigma(ie) 
	2	  -cpi*api(ie)*api(ie))
	      ra3(ia,ib) = ra3(ia,ib)
	1	  +TWO*y2*fc*we(ie)*(asigma(ie)*cpis*api(ie))
	    enddo ! ie

810	    continue
	  enddo ! ia
	enddo ! ib

	return
	end ! power_distribution1
C
	subroutine trapz2(ra,area)

C[subroutine_header_comments]
C  Size parameters:
	integer*4	E_SZ,P_SZ
	parameter	(E_SZ=10000,P_SZ=201)

C  Declarations of scalars:
	integer*4	ia,ib
	real*8		area,sum,wx,wy

C  Declarations of arrays:
	real*8		ra(P_SZ,*)

C  Labeled constants:
	real*8		ZERO,ONE,HALF
	parameter	(ZERO=0.0D0,ONE=1.0D0,HALF=0.5D0)

C  Common blocks:
	integer*4	nxp,nyp
	real*8		dxp,dyp
	real*8		xp(P_SZ),yp(P_SZ),cx(P_SZ),cy(P_SZ)

	common		/pinhole/   xp,yp,cx,cy,nxp,nyp,dxp,dyp

	sum = ZERO
	do ib=1,nyp
	  if (ib .eq. 1 .or. ib .eq. nyp) then
	    wy = HALF
	  else
	    wy = ONE
	  end if ! ib
	  do ia=1,nxp
	    if (ia .eq. 1 .or. ia .eq. nxp) then
	      wx = HALF
	    else
	      wx = ONE
	    end if ! ia
	    sum = sum +wx*wy*ra(ia,ib)
	  end do ! ia
	end do ! ib
	area = sum*dxp*dyp

	return
	end ! trapz2
C
	subroutine print_out(isub)

C[subroutine_header_comments]

C  Size parameters:
	integer*4	E_SZ,P_SZ
	parameter	(E_SZ=10000,P_SZ=201)

C  Declarations of scalars:
	integer*4	isub,ia,ib,ie
	real*8		we,power,flux
	real*8		p1,p2,p3,p4

C  Declarations of arrays:
	
C  Fundamental physical constants; Physics Today Aug. 1990:
	real*8		C,ME,MEE,EC,H,HBAR,MUZ,EPSZ
	parameter	(C    =2.99792458D8)	! Speed of light [m/s]
	parameter	(ME   =9.1093897D-31)	! Electron rest mass [kg]
	parameter	(MEE  =0.51099906D0)	! Electron rest mass [MeV]
	parameter	(EC   =1.60217733D-19)	! Elementary charge [C]
	parameter	(H    =6.6260755D-34)	! Planck's constant [Js]
	parameter	(HBAR =1.05457266D-34)	! Planck's constant/2Pi [Js]
	parameter	(MUZ  =1.2566370614D-6)	! Permeability of vacuum [NA-2]
	parameter	(EPSZ =8.854187817D-12)	! Permittivity of vacuum [Fm-1]

C  Conversion factors:
	real*8		C_MA_A,C_CM_M
	parameter	(C_MA_A=1.0D-3,C_CM_M=1.0D-2)

C  Labeled constants:
	real*8		PI,PIHALF,TWOPI
	parameter	(PI    =3.1415 92653 58979 32384 62643D0)
	parameter	(PIHALF=1.5707 96326 79489 66192 31322D0)
	parameter	(TWOPI= 6.2831 85307 17958 64769 25287D0)
	real*8		BW
	parameter	(BW=1.0D-3) ! 0.1%
	real*8		ZERO,ONE,TWO,HALF
	parameter	(ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,HALF=0.5D0)

C  Common blocks:
	logical*4	lang
	logical*4	lbm
	integer*4	mode
	integer*4	nxp,nyp
	integer*4	ne
	real*8		d
	real*8		emin,emax,xpc,ypc,xps,yps
	real*8		k,kx,ky,k3,gamma,er,len,n,b0,ec0,psi0
	real*8		fac,facs,faca,facp,facp1
	real*8		dxp,dyp
	real*8		de
	real*8		e(E_SZ)
	real*8		spec0(E_SZ),spec1(E_SZ),spec2(E_SZ),spec3(E_SZ)
	real*8		ra0(P_SZ,P_SZ),ra1(P_SZ,P_SZ),
	1		ra2(P_SZ,P_SZ),ra3(P_SZ,P_SZ)
	real*8		xp(P_SZ),yp(P_SZ),cx(P_SZ),cy(P_SZ)

	common		/prti/      lang
	common		/prtr/      d,emin,emax,xpc,ypc,xps,yps
	common		/calci/     mode,lbm
	common		/calcr/	    k,kx,ky,k3,gamma,er,len,n,b0,ec0,psi0
	common		/factor/    fac,facs,faca,facp,facp1
	common		/pinhole/   xp,yp,cx,cy,nxp,nyp,dxp,dyp
	common		/energyi/   ne
	common		/energyr/   e,spec0,spec1,spec2,spec3,de
	common		/spectra/   ra0,ra1,ra2,ra3

C+ Print
C  -----------------------------------------------------------------------------
	flux  = ZERO
	power = ZERO

	if (isub .eq. 1) then ! space distribution
	  call trapz2(ra0,flux) ! get integrated flux over observation area
	  flux  = fac*flux      ! ph/s/0.1%bw
	  power = flux/BW*EC    ! W/eV
	  if (lang) then
	    write (2,232) 'Angular flux density distribution for ',
	1		  emin,' eV (ph/s/mr^2/0.1%bw):'
	  else
	    write (2,232) 'Irradiance for ',emin,' eV  @ ',d,' m',
	1		  ' (ph/s/mm^2/0.1%bw):'
	  endif ! lang
	  write (2,250) 'Integrated flux  ',flux, ' ph/s/0.1%bw'
	  write (2,252) 'Integrated power ',power,' W/eV'
	  write (2,200)
	  if (lang) then
	    write (2,200) '   x(mr)     y(mr)    Ang. flux density',
	1	              '    p1      p2      p3      p4'
	  else
	    write (2,200) '   x(mm)     y(mm)    Irradiance',
	1	              '           p1      p2      p3      p4'
	  endif ! lang
	  do ia=1,nxp
	    do ib=1,nyp
	      if (ra0(ia,ib) .gt. ZERO) then
		p1 = ra1(ia,ib)/ra0(ia,ib)
		p2 = ra2(ia,ib)/ra0(ia,ib)
		p3 = ra3(ia,ib)/ra0(ia,ib)
		p4 = ONE -sqrt(p1*p1 +p2*p2 +p3*p3) ! unpolarized
	      else
		p1 = ZERO
		p2 = ZERO
		p3 = ZERO
		p4 = ZERO
	      endif
	      write (2,210) xp(ia),yp(ib),ra0(ia,ib),p1,p2,p3,p4
	    enddo ! ib
	  enddo ! ia

	elseif (isub .eq. 2 .or. isub .eq. 3) then ! spectral distributions
	  do ie=1,ne
	    if (ie .eq. 1 .or. ie .eq. ne) then
	      we = HALF
	    else
	      we = ONE
	    endif ! ie
	    flux  = flux  +we*spec0(ie)/e(ie) ! -> ph/s/mm^2  /eV or
					      ! -> ph/s/mrad^2/eV
	    power = power +we*spec0(ie)       ! -> W/mm^2/eV or W/mrad^2/eV
	  enddo ! ie
	  flux  = flux /BW   *de ! ph/s or ph/s/mm^2 or ph/s/mrad^2
	  power = power/BW*EC*de ! W or W/mm^2 or W/mrad^2

 	  if (mode .eq. 2) then ! angular/spatial flux density spectrum
	    if (lang) then
	      write (2,231) 'Angular flux density spectrum for ',
	1		     xpc,' mr,',ypc,' mr',
	2		    ' (ph/s/mr^2/0.1%bw):'
	      write (2,270) 'Integrated flux  density',flux,
	1		    ' ph/s/mr^2'
	      write (2,254) 'Integrated power density',power,
	1		    ' W/mr^2'
	      write (2,200)
	      write (2,200) 'Energy(eV)   Ang. flux density',
	1		    '    p1      p2      p3      p4'
	    else
	      write (2,230) 'Irradiance for ',xpc,' mm,',ypc,
	1		    ' mm  @ ',D,' m',
	2		    ' (ph/s/mm^2/0.1%bw):'
	      write (2,270) 'Integrated flux  density',flux,
	1		    ' ph/s/mm^2'
	      write (2,254) 'Integrated power density',power,
	1		    ' W/mm^2'
	      write (2,200)
	      write (2,200) 'Energy(eV)   Irradiance',
	1		    '           p1      p2      p3      p4'
	    endif ! lang

	  elseif (mode .eq. 3) then ! brightness
	    write (2,200) 'On-axis brightness',
	1		  ' (ph/s/mr^2/mm^2/0.1%bw):'
	    write (2,200)
	    write (2,200) 'Energy(eV)   Brightness'

	  elseif (mode .eq. 4) then ! flux spectrum through a pinhole
	    if (lang) then
	      write (2,234) 'Flux through ',xps,' mr x',yps,
	1		    ' mr pinhole at',xpc,' mr,',ypc,' mr:'
	    else
	      write (2,234) 'Flux through ',xps,' mm x',yps,
	1		    ' mm pinhole at',xpc,' mm,',ypc,' mm @',
	2		    d,' m:'
	    endif ! lang
	    write (2,250) 'Integrated flux  ',flux, ' ph/s'
	    write (2,254) 'Integrated power ',power,' W'
	    write (2,200)
	    write (2,200) 'Energy(eV)   Flux(ph/s/0.1%bw)',
	1		  '    p1      p2      p3      p4'

	  elseif (mode .eq. 5) then ! angle-integrated spectrum
	    write (2,200) 'Angle-integrated spectrum:'
	    write (2,250) 'Integrated flux  ',flux, ' ph/s'
	    write (2,254) 'Integrated power ',power,' W'
	    write (2,200)
	    write (2,200) 'Energy(eV)   Flux(ph/s/0.1%bw)',
	1                 '    p1      p2      p3      p4'
	  endif ! mode

	  do ie=1,ne
	    if (spec0(ie) .gt. ZERO) then
	      p1 = spec1(ie)/spec0(ie)
	      p2 = spec2(ie)/spec0(ie)
	      p3 = spec3(ie)/spec0(ie)
	      p4 = ONE -sqrt(p1*p1 +p2*p2 +p3*p3) ! unpolarized
	    else
	      p1 = ZERO
	      p2 = ZERO
	      p3 = ZERO
	      p4 = ZERO
	    endif
	    if (mode .eq. 3) then ! brightness
	      write (2,220) e(ie),spec0(ie)
	    else
	      write (2,220) e(ie),spec0(ie),p1,p2,p3,p4
	    endif
	  enddo ! ie

	elseif (isub .eq. 4 .or. isub .eq. 5) then ! power distribution
	  call trapz2(ra0,power)
	  power = fac*power ! W
	  if (lang) then
	    write (2,200) 'Power density distribution (W/mr^2):'
	  else
	    write (2,234) 'Power density distribution @ ',
	1                  d,' m (W/mm^2):'
	  endif ! lang
	  write (2,200)
	  write (2,254) 'Integrated power ',power,' W'
	  write (2,200)
	  if (lang) then
	  write (2,200) '   x(mr)     y(mr)    Power density',
	1	        '        p1      p2      p3      p4'
	  else
	  write (2,200) '   x(mm)     y(mm)    Power density',
	1	        '        p1      p2      p3      p4'
	  endif ! lang

	  do ia=1,nxp
	    do ib=1,nyp
	      if (ra0(ia,ib) .gt. ZERO) then
		p1 = ra1(ia,ib)/ra0(ia,ib)
		p2 = ra2(ia,ib)/ra0(ia,ib)
		p3 = ra3(ia,ib)/ra0(ia,ib)
		p4 = ONE -sqrt(p1*p1 +p2*p2 +p3*p3) ! unpolarized
	      else
		p1 = ZERO
		p2 = ZERO
		p3 = ZERO
		p4 = ZERO
	      endif
	      write (2,210) xp(ia),yp(ib),ra0(ia,ib),p1,p2,p3,p4
	    enddo ! ib
	  enddo ! ia

	endif ! isub
C- 
C  -----------------------------------------------------------------------------

200	format(' ',8a)
c210	format(' ',f8.3,tr2,f8.3,tr3,1pe13.6,tr2,0p,4(tr2,f6.3))
c220	format(' ',f10.3,tr2,1pe13.6,tr2,0p,4(tr2,f6.3))
210	format(' ',f8.3,tr2,f8.3,tr3,1pe14.6e3,tr2,0p,4(tr2,f6.3))
220	format(' ',f10.3,tr2,1pe14.6e3,tr2,0p,4(tr2,f6.3))
230	format(' ',(a,f7.3,a),2(f7.3,a),a)
231	format(' ',(a,f7.3,a),1(f7.3,a),a)
232	format(' ',(a,f7.1,a),(f7.3,a),a)
234	format(' ',(a,f7.3,a),4(f7.3,a))
250	format(' ',(a,1pe10.3,a))
252	format(' ',(a,f10.6,a))
254	format(' ',(a,f10.2,a))
270	format(' ',(a,1pe10.3,a))

	return
	end ! print_out
C
	subroutine check_file(input)
c check if a file exists and quit if it does not exist

	character*180 input
	logical*4 file_exists

	inquire(file=input, exist=file_exists)

	if (file_exists .neqv. .true.) then
	  print 200, '&check_file-F-NOFILE, File not found ... '
	  print 200, trim(input)
	call exit(0)
	endif

200	format(' ',8a)
	end
