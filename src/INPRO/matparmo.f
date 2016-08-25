C	 
	Subroutine matparmo(dspace,psio,psih,psihb,my,bfac,gamma0,gammah
     1,facP,IG)
C
        implicit none
C-------------------------------------------------------------------------
C
C	Written by Michael Krisch, last interaction 14.2.92
C	(Modified by Carlos Giles 4.12.92 to be used as a 
C	subroutine of the program INPRO.)
C	(Modified by M. Sanchez del Rio 94-02-24 to be compiled
C	in Unix and adapt a user interface)
C       (Modified by R. Butzbach 96-12-20: Extinction length etc.
C                   (corrected) now for each polarization seperately)
C       Modified:  Nov  2013  RJD, ASD/APS. Changed order of variables in "parameter" and
C                       "structure" common blocks so that the integer variables now appear last to avoid
C                       alignment issues and compilation warnings.
C
C	This program calculates the crystal parameters for a given 
C	wavelength/energy, reflection and temperature. 
C
C	Input:   - wavelength/energy
C	         - unit of wavelength/energy input
C	         - material
C	         - reflection
C	         - temperature
C	
C	Output:  - filename                        FILENA
C		 _ chemical symbol                 SB
C		 _ temperature in K                TEMP
C		 _ Debye Waller factor             DWF
C		 _ wavelength in nm                LAMBDA
C		 _ crystal reflex                  HM,KM,LM
C		 - Bragg angle                     THETAB
C		 - asymmetry angle                 ALPHA
C		 - dspacing in mn                  DSPACE
C		 - polarization factor             CPOL
C		   ( unpolarized )       
C		 - Polarization factor             FACP
C		   ( 100% p-polarized )  
C                - Polarization factor             CPOLS
C                  ( 100% s-polarized )
C		 - asymmetry factor b              BFAC
C		 - refraction corr. in "           REFCOR
C		 - scattering factors f0,f1,f2     F0,F1,F2
C		 - cross section sigma in barns    SIGMA
C	 	 - abs. coefficient mu (1/mm)      MY
C		 - spec. absorp. in cm**2/g        MYRHO
C		 - mean abs. coefficient (1/mm)    MYMEAN
C		 - Dielectric susceptibilities     CHIOR,CHIOI,CHIHR,CHIHI
C		 - Darwin width                    DARWID
C		 - extinction depth/length in um   EXTDEP,EXTLEN
C		 - Pendelloesung period in um      PENPER
C		 - absorption depth in um          ABSDEP
C
C	It needs the following subroutines:
C
C	   PARAMETERMO(NCHOI,SB,temp): parameter file
C	   ATOMFACM(MOZ,F0): formfactor
C	   SCATTER(SB,ENERGY,EINHEI,F1,F2): dispersion corrections	
C	   DEBYENEWM(ZMASS,TEMP,TDEB,DWF): Debye-Waller factor 
C	   STRUFACNEWM(CHI0R,CHI0I,CHIHR,CHIHI,CHIHRB,CHIHIB)
C
C	It needs the following functions:
C
C	   EUCONV(ENERGY,EINHEI,'M'): Conversion of input unit into m
C
C
	Common/general/lambda,thetab
	Common/parameter/zmass,par,pcr,tdeb,ne,nue,moz,str
 	Common/structure/vez,f0,f1,f2,dwf,hm,km,lm
C
	Integer*4         hm,km,lm,moz,nchoi,nue,ne,str,i,IG
C
	Real*4		energy,f1,f2,lambdas,EUCONV
	Real*8		lambda,temp,thetab,alpha,zmass
	Real*8		dspace,dwf,f0,my,mymean,myrho,root,tdeb
	Real*8 		bfac,gamma0,gammah,par,pcr,pi,re,sigma,u,vez
	Real*8		refcor,chi0r,chi0i,darwid
	Real*8		absdep,cpol,extdep,extlen,facP,factor,penper
	Real*8		term1,term2,kappa,cpolS
	Real*8		cr,crb,ci,cib,chih,chihb, offset
C
	Complex*16		chihr,chihi,chihrb,chihib
	Complex*16		psio,psih,psihb,IZ
C
	Dimension		moz(2),zmass(2),nue(2),f0(2),f1(2),f2(2)
	Dimension		tdeb(2),sigma(2),my(2),myrho(2)
	Dimension		dwf(2),sb(2)
	Character		einhei*12
	Character         sb*2 !,filena*60			!potee         
C
C	physical constants: 
C
	Parameter         (re=2.81777d-15,pi=3.141592653589793238d0)
	Parameter         (u=1.6604d-24)
C ____________________________________________________________________________
C|                                                                            |
C| ***  Short explanation  ***                                                |
C|____________________________________________________________________________|
C
C	write (*,*)
C	write (*,*)
C	write (*,*)
C	write (*,*)
C	write (*,*)' This program calculates the crystal parameters for a given' 
C	write (*,*)' wavelength/energy, reflection and temperature.'
C	write (*,*)
C	write (*,*)' The protocol will be stored in the same directory!'
C	write (*,*)
C ____________________________________________________________________________
C|                                                                            |
C| ***  parameter setup ***                                                   |
C|____________________________________________________________________________|
C
	sb(1) = ''
	sb(2) = ''
10	energy=0.e0
	lambda = 0.d0
	lambdas= 0.e0
	thetab = 0.d0
	temp=0.d0
	do 30 i = 1,2
		f0(i) = 0.d0
		f1(i) = 0.e0
		f2(i) = 0.e0
		zmass(i) = 0.d0
		tdeb(i) = 0.d0
		moz(i) = 0
		nue(i) = 0
30 	continue
	hm=0
	km=0
	lm=0
C ____________________________________________________________________________
C|                                                                            |
C| *** data input ***                                                         |
C|____________________________________________________________________________|
C
100	continue
	write (*,*)
	write (*,*)
	write (*,*) ('  DATA INPUT')
	write (*,*)
	write (*,*)
	write (*,*)' Your choice of the crystal:'
	write (*,*)
	write (*,*)' Si         (1)'
	write (*,*)' Ge         (2)'
	write (*,*)' diamond    (3)'
	write (*,*)' GaAs       (4)'
	write (*,*)' GaP        (5)'
	write (*,*)' InAs       (6)'
	write (*,*)' InP        (7)'
	write (*,*)' InSb       (8)'
	write (*,*)' SiC        (9)'
	write (*,*)
	write (*,*)' CsF        (10)'
	write (*,*)' KCl        (11)'
	write (*,*)' LiF        (12)'
	write (*,*)' NaCl       (13)'
	write (*,*)
	write (*,*)' graphite   (14)'
	write (*,*)
	write (*,*)' beryllium  (15)'
	write (*,*)
	write (*,120)
120	format (' Your choice: ',$)
	read *, nchoi
	if (nchoi.lt.1 .or. nchoi.gt.15) goto 4000
	write (*,*)
	write (*,*) 'Energy units: (EV,KEV,ANGSTROMS,NANOMETERS,NM):'
C changed for irix compilation
C	read (*,*) einhei
	read (*,'(a)') einhei
	write (*,*)  'Then the wavelength/energy value: '
	read (*,*) energy
	write (*,*)
	write (*,*)
	write (*,160)
160	format (' Bragg planes of the crystal reflection:')
	write (*,*)
	write (*,170)
170	format ('              h: ',$)
	read *, hm
	write (*,180)
180	format ('              k: ',$)
	read *, km
	write (*,190)
190	format ('              l: ',$)
	read *, lm
	write (*,*)
	write (*,*)
	write(*,*)' angle of asymmetry? Attention, the angle is positive
     1,if the '
	write(*,*)' reciprocal lattice vector is rotated clockwise 
     1from its'
	write(*,*)' symmetrical position.'
	write (*,*)
	write (*,200)
200	format (' asymmetry angle alpha in degrees: ',$)
	read *, alpha
	write (*,*)
	write (*,*)
	write (*,210)
210	format (' Temperature in Kelvin:',$)
	read *, temp
	write (*,*)
	write (*,*)
C
C	Conversion of energy/wavelength input into m
C

	lambdas = EUCONV(energy,einhei,'M')	

	
	lambda = 1.d0*lambdas
	
C
C	Conversion of alpha into rad

	alpha = alpha*pi/180.d0
	
C ____________________________________________________________________________
C|                                                                            |
C| *** calculation of necessary parameters ***                                |
C|____________________________________________________________________________|
C
C	Calling the parameter file in order to provide the necessary 
C	data
C	**C.G. parametermo have a temperature correction for the
C	lattice parameters for some crystals only (Si,Ge,diamond,
C	LiF,NaCl,KCl).
C

csrio	write (*,*) 'calling parametermo: nchoi,sb,temp',nchoi,sb,temp
	call parametermo(nchoi,sb,temp)
	

csrio	write (*,*) 'returning parametermo: nchoi,sb,temp',nchoi,sb,temp


	if (str.eq.1 .or. str.eq.2) then      	     ! cubic structure
		vez = par**3                         ! unit cell volume
		root = hm**2 + km**2 + lm**2
		dspace = par/(dsqrt(root))            ! d-spacing
		thetab = dasin(lambda/2.d0/dspace)      ! Bragg angle
	endif
	if (str.eq.3 .or. str.eq.4) then	     ! hexagonal structure
		vez = dsqrt(3.d0)/2.d0*par**2*pcr       ! unit cell volume
		term1 = 4.d0*(hm**2 + km**2 + hm*km)/(3.d0*par**2) 
		term2 = lm**2/pcr**2
		dspace = dsqrt(1.d0/(term1 + term2))   ! d-spacing
		thetab = dasin(lambda/2.d0/dspace)     ! Bragg angle
	endif
C
C
C	*** *asymmetry angle used in the same convention as the one
C 	used in the subroutine perfect_crystalsmo written by C. Vettier
C
C

	if (IG.eq.2 .or. IG.eq.-2) goto 1700
	gamma0 = dsin(thetab + alpha)       ! valid for Bragg case
	gammah = - dsin(thetab - alpha)     ! valid for Bragg case
	goto 1750
1700	gamma0 = dcos(thetab + alpha)       !valid for Laue case
	gammah =  dcos(thetab - alpha)     !valid for Laue case
1750	bfac = gamma0/gammah
C
	do 1800 i = 1,ne

csrio	write (*,*) 'calling atomfacm moz,f0',moz(i),f0(i)

	call atomfacm(moz(i),f0(i))	
	
csrio	write (*,*) 'returning atomfacm moz,f0',moz(i),f0(i)
csrio	write (*,*) 'calling scatter sb,energy,einhei,f1,f2',
csrio     1sb(i),energy,einhei,f1(i),f2(i)

	call SCATTER(sb(i),energy,einhei,f1(i),f2(i))	


csrio	write (*,*) 'returning scatter sb,energy,einhei,f1,f2',
csrio     1sb(i),energy,einhei,f1(i),f2(i)
csrio	write (*,*) 'calling  debyenewm: zmass,temp,tdeb,dwf',
csrio     1zmass(i),temp,tdeb(i),dwf(i)

		call debyenewm(zmass(i),temp,tdeb(i),dwf(i))	
		
csrio	write (*,*) 'returning  debyenewm: zmass,temp,tdeb,dwf',
csrio     1zmass(i),temp,tdeb(i),dwf(i)
1800	continue
C
csrio	write (*,*) 'calling strufacnewm: ',
csrio     1'chi0r,chi0i,chihr,chihi,chihrb,chihib',
csrio     1chi0r,chi0i,chihr,chihi,chihrb,chihib

	call strufacnewm(chi0r,chi0i,chihr,chihi,chihrb,chihib)
		
csrio	write (*,*) 'returning strufacnewm: ',
csrio     1'chi0r,chi0i,chihr,chihi,chihrb,chihib',
csrio     1chi0r,chi0i,chihr,chihi,chihrb,chihib
C	
C ____________________________________________________________________________
C|                                                                            |
C| *** calculation of crystal features ***                                    |
C|____________________________________________________________________________|
C
	cpol = (dabs(dcos(2.d0*thetab)) + 1.d0)/2.d0      ! polarization factor
        facP = dabs((dcos(2.d0*thetab)))                ! 100% p-polarized
C
C	Converting the susceptibilities into real variables
C
C	modulus of the real part of chih
	cr = - real(abs(chihr))
C	modulus of the real part of chih bar
	crb = - real(abs(chihrb))
C	modulus of the imaginary part of chih
	ci = - real(abs(chihi))
C	modulus of the imaginary part of chih bar
	cib = - real(abs(chihib))
C	modulus of chih
	chih = - sqrt(cr**2 + ci**2)
C	modulus of chih bar
	chihb = - dsqrt(crb**2 + cib**2)
C	
C	obtaining final complex susceptibilities 'Carlos Giles
	psio = cmplx(chi0r,chi0i)
	IZ = cmplx(0.d0,1.d0)
	psih = chihr+IZ*chihi
	psihb = chihrb+IZ*chihib
C
C	only valid for centrosymmetrical crystals
	kappa = real(abs(chihi))/real(abs(chihr))
C
C	refraction correction in rad
C
	refcor = -((1.d0 - bfac)*(chi0r + chi0i*dtan(kappa)))
     1           /(2.d0*dsin(2.d0*thetab))
C
	mymean = 0.0D0 ! added by srio
	do 1850 i = 1,ne
		sigma(i) = f2(i)*2.d0*re*lambda	     ! cross section in m**2
		my(i) = nue(i)/vez*sigma(i)	     ! absorption in 1/m
		myrho(i) = sigma(i)/(zmass(i)*u)     ! spec. abs. in m**2/g
C		
C		mean absorption coefficient
		mymean = mymean + my(i)
1850	continue
	mymean = mymean/ne
C
C	extinction parameters
C	
C       polarization factor, sigma pol (just a dummy):
        cpolS = 1
C
C	Darwin width in the symmetrical case
	darwid = 2.d0*cpolS*dsqrt(dabs(chih*chihb))/
     1           (dsin(2.d0*thetab)*dcos(kappa))
C
	factor = lambda*dsqrt(gamma0*dabs(gammah))*dcos(kappa)/
     1           (cpolS*dsqrt(dabs(chih*chihb)))
	extdep = 1.d0/(2.d0*pi)*factor                 ! extinction depth
	extlen = extdep/gamma0                       ! extinction length
	penper = factor                              ! Pendell. period
	absdep = 1.d0/(mymean*(1.d0/gamma0 + 1.d0/abs(gammah)))  !absorpt. depth
C

        offset=  chih*chih*sin(thetab+thetab)/(lambda*4.8484e-3)
C ____________________________________________________________________________
C|                                                                            |
C| *** conversion of parameters for output ***                                |
C|____________________________________________________________________________|
C
csrio	thetab = thetab*180.0/pi                  ! from rad into degree
csrio	alpha = alpha*180.0/pi                    ! from rad into degree
csrio	refcor = refcor*206264.8063               ! from rad to arcsec
csrio	darwid = darwid*206264.8063               ! from rad to arcsec
csrio	lambda = lambda*1.0e+09                   ! from m into nm
csrio	dspace = dspace*1.0e+09                   ! from m into nm
csrio	extdep = extdep*1.0e+06                   ! from m into um
csrio	extlen = extlen*1.0e+06                   ! from m into um
csrio	penper = penper*1.0e+06                   ! from m into um
csrio	absdep = absdep*1.0e+06                   ! from m into um
csrio	do 1999 i = 1,ne
csrio		sigma(i) = sigma(i)*1.0e+28       ! from m**2 into barn
csrio		my(i) = my(i)*1.0e-03             ! from 1/m into 1/mm
csrio		myrho(i) = myrho(i)*1.0e+04       ! from m**2/g into cm**2/g
csrio1999	continue
csrio	mymean = mymean*1.0e-03       		  ! from 1/m into 1/mm
C ____________________________________________________________________________

C ____________________________________________________________________________
C|                                                                            |
C| *** comment and data storage ***                                           |
C|____________________________________________________________________________|
C
C	output for file 
	open (unit=15,file='inpro.par',status='unknown')
	write (15,*) '================================================'
	write (15,*) '=================  Inpro output ================'
	write (15,*) '================================================'
	write (15,*) ' '
	write (15,*) ' '
	write (15,*) '============ ==  Crystal Data =================='
	write (15,*) 'Crystal: ',sb(1),sb(2)
	write (15,*) 'temperature: ',temp, ' Kelvin'
	do 2025 i = 1,ne
		write (15,2030) sb(i),dwf(i)
2030		format (' Debye-Waller factor (',A2,'): ',G14.6) 
2025	continue
	write (15,*) 'lambda: ',lambda*1.0e+9,' nm'
	write (15,*) 'k0    : ',1/(lambda*1.0e+9),' nm^-1'
	write (15,*) 'Energy: ',12.39854/lambda*1.0e-10,' keV'
	write (15,*) 'Miller Indeces: ',hm,km,lm
	write (*,*) 'Theta Bragg: ',thetab*180.0/pi,' degrees'
	write (15,*) 'Theta Bragg: ',thetab*180.0/pi,' degrees'
	write (15,*) 'Asymmetry angle: ',alpha*180.0/pi,' degrees'
	write (15,*) 'Interplanar distance d: ',dspace*1.0e+9,' nm'
	write (15,*) 'polarization factor ( unpolarized radiation): ',
     1cpol
        write (15,*) 'polarization factor (100% p-polarized)     : ',
     1facP
	write (15,*) 'asymmetry factor b: ',bfac
	write (15,*) 'refraction correction: ',
     1refcor*206264.8063,' arcsec'
	write (15,*) ' '
	write (15,*) ' '
	write (15,*) '====== Scattering factors and related items ===='
	do 2105 i = 1,ne
		write (15,*)
		write (15,2107) sb(i)
2107		format (' scattering factors for ',A2)
		write (15,*)
		write (15,2110) f0(i),f1(i),f2(i)
2110		format (' f0: ',G14.6,4X,' f1: ',G14.6,4X,' f2: ',G14.6)
		write (15,*)
		write (15,*) 'cross section: ',sigma(i)*1.0e+28,' barn'
		write (15,*) 'absorption coefficient: ',
     1my(i)*1.0e-3,' mm^(-1)'
		write (15,*) 'specific absorption: ',myrho(i)*1.0e+4,
     1' cm^2/g'
		write (15,*)
2105	continue
	write (15,*) 'mean absorption coefficient: ',mymean*1.0e-3,
     1' mm^(-1)'
	write (15,*) ' '
	write (15,*) ' '
	write (15,*) '======== Dielectric susceptibilities: =========='
	write (15,*)
	write (15,2150) chi0r,chi0i
2150	format (' chi0r = ',G14.6,'  chi0i = ',G14.6)
	write (15,*)
	write (15,*) ' real and imaginary part of chih(real):'  
	write (15,2160) chihr
2160	format (2X,G14.6,3X,G14.6)
	write (15,*)
	write (15,*) ' real and imaginary part of chih(imag):'
	write (15,2161) chihi
2161	format (2X,G14.6,3X,G14.6)

	write (15,*)
	write (15,*) ' modulus of chih(real):'
	write (15,2162) cr
2162	format (2X,G14.6)
	write (15,*)
	write (15,*) ' modulus of chih(imag):'
	write (15,2163) ci
2163	format (2X,G14.6)
	write (15,*)
	write (15,*) ' real and imaginary part of chih(real) bar:'  
	write (15,2164) chihrb
2164	format (2X,G14.6,3X,G14.6)
	write (15,*)
	write (15,*) ' real and imaginary part of chih(imag) bar:'
	write (15,2165) chihib
2165	format (2X,G14.6,3X,G14.6)
	write (15,*)
	write (15,*) ' modulus of chih(real) bar:'
	write (15,2166) crb
2166	format (2X,G14.6)
	write (15,*)
	write (15,*) ' modulus of chih(imag) bar:'
	write (15,2167) cib
2167	format (2X,G14.6)
	write (15,*) ' '

	write (15,*)
	write (15,*) ' Psi_0:'
        write (15,2161) psio 
	write (15,*)
	write (15,*) ' Psi_h:'
        write (15,2161) psih
	write (15,*)
	write (15,*) ' Psi_h_bar:'
        write (15,2161) psihb

	write (15,*) ' '
	write (15,*) '=============== Other values: =================='
	write (15,*) ' '
	write (15,*) ' sigma-polarization:'
	write (15,*) '    Darwin width (symmetrical case):',
     1darwid*206264.8063,' arc sec'
	write (15,*) '    extinction depth: ',extdep*1.0e+6,' microns'
	write (15,*) '    extinction length: ',extlen*1.0e+6,' microns'
	write (15,*) '    pendelloesung period: ',penper*1.0e+6,
     1' microns'
        write (15,*) ' '
        write (15,*) ' pi-polarization:'
        write (15,*) '    Darwin width (symmetrical case):',
     1darwid*206264.8063*facP,' arc sec'
        write (15,*) '    extinction depth: ',extdep/facP*1.0e+6,
     1' microns'
        write (15,*) '    extinction length: ',extlen/facP*1.0e+6,
     1' microns'
        write (15,*) '    pendelloesung period: ',penper/facP*1.0e+6,
     1' microns'
        write (15,*) ' '
        write (15,*) ' absorption depth: ',absdep*1.0e+6,' microns'
        write (15,*) ' angular offset/effective thickness: ',offset,
     1' arcsec/mm'
	continue
C	output for protocol
	goto 5000
C ____________________________________________________________________________
C|                                                                            |
C| *** error messages ***                                                     |
C|____________________________________________________________________________|
C
4000	write (*,*)
	write (*,*)
	write(*,*)' INPUT ERROR, try once more!!!'
	write(*,*)' nchoi must be 1 or 2!!!'
	goto 100
C ____________________________________________________________________________
C|                                                                            |
C| *** file closure ***                                                       |
C|____________________________________________________________________________|
C
5000	close (15)
	write (*,*) 'closing 15'
	end
