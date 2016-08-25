                                            	PROGRAM INPRO
C
C	*************************************************************
C	***            PROGRAM INTRINSIC PROFILES   (INPRO)         *
C	***   Calculates intrinsic profiles using dynamical theory. *
C	*** The expressions used correspond to equations 3.130 and  *
C	*** 3.137 from Zachariasen book. The program allows the cal-*
C	*** culation of intrinsic profiles of the following perfect *
C	*** crystals: Si, Ge, diamond, GaAs, GaP, InAs, InSb, SiC,  *
C	*** CsF, KCl, LiF, NaCl, graphite, Beryllium.               *
C	*** The parameters for these crystals are calculate by      *
C	*** different subroutines written by M.Krisch.              *
C	*** The intrinsic profiles are calculated using the subrou- *
C	*** tine PERFECT_CRYSTALS written by C. Vettier.            *
C	***                                                         *
C	*************************************************************     
C
C       Modified:  Nov  2013  RJD, ASD/APS. Changed order of variables in "parameter" and
C                       "structure" common blocks so that the integer variables now appear last to avoid
C                       alignment issues and compilation warnings. NOTE: Common blocks are commented but
C                       changes were made in case those will be uncommented in the future.
C
C
	IMPLICIT NONE
C
	INTEGER I,IG,schoi,points
c 	INTEGER*4 ne,nue,moz,str,hm,km,lm	! never used
c 	Real*4 f1,f2				! never used
	Real*8 minlim,maxlim,step,sintb
	REAL*8 lambda,thetab,thetap,thet,D_W,lnom
	Real*8 angle!,zmass,par,pcr,tdeb,vez,f0,dwf	! never used
C
	REAL*8 ECONVEL,PI,R0,SEC_RAD
	Real*8 dspace,my,bfac,facp,gamma0,gammah,t0,alpha
	Real*8 refs,refp
C
	COMPLEX*16 ZS,ZP
	COMPLEX*16 psi0,psih,psihb
	COMPLEX*16 IZ
	Character filena*1024
	Character*1024 datadir
C
C
C
	COMMON/GENERAL/lambda,thetab
C	COMMON/parameter/zmass,par,pcr,tdeb,ne,nue,moz,str
C	COMMON/STRUCTURE/vez,f0,f1,f2,dwf,hm,km,lm
	COMMON/CONSTANTS/ECONVEL,PI,R0,SEC_RAD,IZ

        COMMON/DIFFRACTION/psi0,psih,psihb,dspace,bfac,facp,
     1                  gamma0,gammah,t0,alpha,lnom,IG
 



	COMMON/DATADIR/datadir
C
C
C	===================================================
C	**** Title and introduction
C	new page here
	write (*,*)
	write (*,*)
	write (*,*)
	write (*,5)
5	format ('1','             Program INtrinsic PROfiles')
	write (*,*)
	write (*,*)
	write (*,*)'               Optics Group - ESRF'
	write (*,*)
	write (*,*)
	write(*,*)'This program calculates reflectivity and transmissi-'
	write(*,*)'vity intrinsic profiles of crystal plates, using the'
	write(*,*)'dynamical theory for x-ray diffraction.'
	write(*,*)'The program allows the calculation of intrinsic'
	write(*,*)'profiles of the following perfect crystals:'
	write (*,*)
	write(*,*)'	Si, Ge, diamond, GaAs, GaP, InAs, InSb, Sic,'
	write(*,*)'	CsF, KCl, LiF, NaCl, graphite, beryllium.'
	write (*,*)

csrio ask for data directory (Lahey-Windows)
	write(*,*)'Inpro will read the datafiles from a data directory'
	write(*,*)'Enter the data directorty path with the last
     1  delimiter'
	write(*,*)' [colon (:) in MacOS]: '
	read(*,'(a)') datadir
	write(*,*) 'Datadir is: ',datadir
C
C	******** =====================================
C	IG is defined in order to select the desired calculation: 
C	+1 Bragg reflected
C	-1 Bragg transmitted
C	+2 Laue reflected
C	-2 Laue transmitted
C	*******
C
	write (*,*)
	write (*,20)
20	format ('  Choice of the case to be calculated:')
	write (*,*)
	write (*,30)
30	format ( ' Reflectivity in Bragg case:   (+1)')
	write (*,40)
40	format ( ' Transmissivity in Bragg case: (-1)')
	write (*,50)
50	format ( ' Reflectivity in Laue case:    (+2)')
	write (*,60)
60	format ( ' Transmissivity in Laue case:  (-2)')
	write (*,70)
70	format ( '      selected case: ',$)
	read(*,*)  IG
C       ====================================================
C       **** input of the thickness t0 in microns
	write (*,*)
	write (*,*)
        write (*,100)
100     format (' thickness in microns: ',$)
        read(*,*)  t0
        t0=t0*1.d-06             !t0 from microns to meters
C
C	========================
C
	ECONVEL=12398.54
	PI=3.141592653589793238D0
	SEC_RAD=3600.D0*180.D0/PI
	R0=2.8179D-5
	IZ=CMPLX(0.D0,1.D0)
C
C	========================
C	**** Call for parameter determination
csrio
	write (*,*) 'ig,t0 = ',IG,t0
	write (*,*) 'calling matparmo'
C
	call matparmo(dspace,psi0,psih,psihb,my,bfac,gamma0,gammah,
     1facp,IG)
	write (*,*) 'returning from  matparmo'
C
C	====================================================
C	**** File opening procedure
C
1000	write (*,*)
	write (*,1100)
1100	format (2X,'input of output filename:')
	write (*,1110)
1110	format (2x,'filename  : ',$)
	read (*,1120) filena
1120	format (A)
	OPEN (UNIT=14,FILE=filena,STATUS='old',ERR=1150)
Csrio     1  defaultfile='[esrf.optics.data].dat')
	CLOSE (UNIT=14)
	write(*,*)'This filename already exists, please choose another.'
	goto 1000
1150	continue
	OPEN (UNIT=14,FILE=filena,STATUS='NEW')
Csrio     1        DEFAULTFILE='[ESRF.OPTICS.DATA].dat')
	continue
C
C	****** =====================================================
C	Determination of the angular limits for the calculation of the 
C	reflectivity profile. 
C	Two possibilities are offered, either the selection of your own angular
C	limits, either the selection of default values corresponding to
C	twice the Darwin width bandpass at each side of the Bragg peak.
C	******
C
	sintb = lambda/(2.d0*dspace)
	thetab= dasin(sintb)
	thetap = sintb + psi0*(1.d0-bfac)/(4.d0*bfac*sintb)
	thetap = dasin(thetap)
	D_W = 2.d0*abs(psih)/(sqrt(abs(bfac))*sin(2.D0*thetab))
C	=============================================================	
C
	thetap =thetap*180.d0/PI       !from rad into degree
	D_W = D_W*SEC_RAD              !from rad into arcsec
C
	write (*,*)
	write (*,*)
	write (*,1200) thetap
1200	format ( '  The value of the Bragg peak position is (in degrees): 
     1',F8.5)
	WRITE (*,*)
	write (*,1205) D_W
1205	format ( '  and the Darwin width is (in arcsec):',F6.2)
	thetap=thetap*PI/180.d0        !from degree into rad
	D_W=D_W/SEC_RAD                !from arcsec into rad
	minlim=(thetap-2.d0*D_W)*180.d0/PI  !from rad into degrees
	maxlim=(thetap+2.d0*D_W)*180.d0/PI  !from rad into degrees
	write (*,*)
	write (*,*)
	print *,
     $ 'Selection of the angular limits for the reflectivity profile:'
	write (*,*)
csrio	write (*,1210) minlim,maxlim
csrio1210	format ( ' (1) Default values: from:',F8.5,'   to:',F8.5)
	write (*,*)  ' (1) Default values: from:',minlim,maxlim,
     1' arc sec'
	write (*,*)
	write (*,*)  ' (2) other values'
csrio	write (*,1220)
csrio1220	format ( ' (2) other values')
	write (*,*)
	write (*,1230)
1230	format (' Your choice: ',$)
	read(*,*) schoi
	if (schoi.eq.1) goto 1300
	minlim=0.d0
	maxlim=0.d0
	write (*,*)
	write (*,*) 'lower limit for the reflectivity profile '//
     1' in arc sec around the Bragg angle'
	read(*,*)  minlim
	write (*,*) 'upper limit for the reflectivity profile '//
     1' in arc sec around the Bragg angle'
	read(*,*)  maxlim
	minlim = thetab*180.d0/PI + minlim/3600.  !from arcsec to degrees
	maxlim = thetab*180.d0/PI + maxlim/3600.  !from arcsec to degrees


1300	minlim=minlim*PI/180.d0         !from degrees into rad
	maxlim=maxlim*PI/180.d0         !from degrees into rad
1310	write (*,*)
	write (*,1320)
1320	format ( ' number of points for the reflectivity curve: ',$)
	read(*,*)  points 
C
	write (*,*)
	thetab=thetab*180.d0/PI            !from rad into degrees
	thetab=thetab*PI/180.d0            !from degrees into rad
	step=(maxlim-minlim)/(points-1)
	write (*,*) 'Number of points = ',points
C
C	
C	******** Loop de variation de theta pour calcul de la reflect.
	DO I=1,points
	thet=minlim+(I-1)*step
C
	Call PERFECT_CRYSTALSMO(thet,lambda,I,ZS,ZP)
C	===================================
C	**** calculation of  the reflectivity and transmissivity
C
	refs=abs(ZS)
	refs=refs**2
	refp=abs(ZP)
	refp=refp**2
C	===================================
C	**** Save of this data in an output file
C
C	conversion of output parameters units for files
C	angular scale in thet-thetab (kinem. Bragg angle) in arcsec,
C	reflectivities for sigma and PI in real numbers between 0 and 1.
C
	angle=(thet-thetab)*SEC_RAD
	write (14,*) angle,refs,refp
csrio	write (14,1500) angle,refs,refp
csrio1500	format (2x,F12.6,2(2X,E12.5))
csrio	write (*,1510) angle,refs,refp
csrio1510	format (2x,F12.6,2(2x,E12.5))
	END DO
C
C
	thetab=thetab*180.d0/PI
	thetap=thetap*180.d0/PI
	write (*,*) ' The Bragg angle is: ',thetab
	write (*,*) ' The Bragg peak position: ',thetap
c
C	===================================
C	*****File closure
	CLOSE (UNIT=14)
C
C	===================================
C
	write(*,*) 'Output data stored in file: ',filena
	write (*,*)
C
	write (*,*) ' End of program'
C
C	===================================
C	
2000	End

