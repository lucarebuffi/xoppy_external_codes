C
C subroutine scatter
C
C written :  Nov. 90    G. Foerstner
C modified:  May  91    M. Krisch
C modified:  Feb  94    M. Sanchez del Rio
C modified:  Nov  2013  RJD, ASD/APS. Changed order of variable ICOUNT in common block
C                       COMMON /BLKGAU/ so that it appears last in the list to avoid
C                       alignment issues and compilation warnings.
C                       Replaced PAUSE statements with READ * statements. (The PAUSE statement was
C                       deleted from the Fortran 95 standard.)
C                       Increased size of array T in subroutine CROLI_PUB to 80 (= the same number 
C                       of elements as in the function DCAIP.)
C
C INPUT:    - energy
C	    - chemical symbol [SYM]
C           
C SUBROUTINE: CALL S_GETF120 to calculate F1, F2 and
C             CALL FIND      to read some material properties:
C				IZ      atomic number (INTEGER*2)
C				ATWT 	atomic weight (REAL*4)
C				RMU     mu(barns/atom)/mu(cm**2/gram) (REAL*4)
C				EMF	energy*mu/f2(ev-cm**2/gram) (REAL*4)
C				NAME	symbol for the element (CHARACTER*2)
C			 	Density 'standard' density (REAL*4)
C 			        written as external subroutine in
C				S_GETF120.FOR
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
	   SUBROUTINE SCATTER(SYM,ENERGY,EINHEIT,FF1,FF2)
C          +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
	   IMPLICIT NONE
C
	   REAL*4	    ENERGY
	   REAL*4    	    FF1,FF2
c 	   REAL*4           k,LAMBDA			! never referenced
c 	   REAL*4           EUCONV			! never referenced
	   REAL*4   	    ATWT,Density,RMU,EMF
	   REAL*4           r0,Avo,PI,h,elec,c_licht
C
	   CHARACTER*2 	    SYM,NAME
	   CHARACTER*(*)    EINHEIT
C
c 	   INTEGER*4	    N
c 	   INTEGER*4	    I
	   INTEGER*4	    IZ
C
	   DATA h/6.6262E-34/ ,elec/1.6022E-19/ ,c_licht/2.99792E08/
	   DATA Avo/6.022E23/,r0/2.81E-15/
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
	   PI = 1.0
	   PI = 2.0 * ASIN(PI)
C
csrio
csrio	write (*,*) 'Calling FIND: SYM,IZ,ATWT,RMU,EMF,NAME,Density',
csrio     1SYM,IZ,ATWT,RMU,EMF,NAME,Density
	   CALL FIND(SYM,IZ,ATWT,RMU,EMF,NAME,Density)	   
csrio	write (*,*) 'Returning FIND: SYM,IZ,ATWT,RMU,EMF,NAME,Density',
csrio     1SYM,IZ,ATWT,RMU,EMF,NAME,Density
C
C calculation of scattering factors for the given energy 
C
	   CALL GETF120(SYM,IZ,ATWT,ENERGY,EINHEIT,FF1,FF2)
C
C	   Unit conversion
	   Density = Density * 1.0E06
C
	   RETURN
	   END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C Subroutine:	GETF120
C 
C Written:
C    3/27/86:          J. Bloch
C Modified:
C    6/22/88, 2/10/90: E. Ziegler
C   11/26/90:	       G. Foerstner
C   05/13/91:          M. Krisch
C
C Description:	
C
C   GETF120 returns the F1 and F2 values for a specified element at a specific 
C   energy. 
C   Between 0.1 keV and 1 keV GETF12 uses the Fortran 77 version of Henke's 
C     original RT-11 Fortran IV subroutines. Linear interpolation is used to 
C     determine the F1 and F2 values at the desired energy.
C   Otherwise f' and f'' are calculated using the algorithm of Cromer and 
C     Liberman as described in the subroutine entitled CROLIB.FOR.
C
C Usage:	CALL GETF120(SYM,IZ,ATWT,E,EUNITS,F1,F2)
C
C INPUT parameters:
C
C      SYM :    CHARACTER string containing the 1 or 2 letter symbol for the 
C               desired element.
C
C      IZ:      INTERGER*2 atomic number
C  
C      ATWT:    REAL*4  atomic weight
C
C      E :      REAL energy for the desired F1 and F2 values.
C
C      EUNITS : CHARACTER string containing the name of the energy units 
C               contained in ENERGY.
C
C OUTPUT parameters:
C
C       F1:     REAL F1 constant.
C
C       F2:     REAL F2 constant.
C
C SUBROUTINE GETF120 CALLED IN S_SCATTER.FOR
C
C EXTERNAL SUBROUTINES:
C
C	- EV = EUCONV(E,EUNITS,'EV')
C	  converts energy in eV using program EUCONV.FOR
C
C	- FIND(SYM,IZ,ATWT,RMU,EMF,NAME,Density) 
C	  directly called in S_SCATTER ( see there for more info )
C
C	- F1F2(SYM,IZ,OUTFIL,ATWT,ENERGY)
C 	  calculating f1,f2 with Henke data
C	  see down there for more info
C	  
C	- CROLI_PUB(SYM,EZ,EV,F1,F2)
C 	  calculating f1,f2 with Cromer data
C	  see S_CROLI_PUB.FOR for more info
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

	   SUBROUTINE GETF120(SYM,IZ,ATWT,E,EUNITS,F1,F2)
C	   ======================================================
C
	   IMPLICIT NONE
C
	   SAVE
	   REAL*4	 E,F1,F2,C,R_IZ
	   REAL*4        OUTFIL(285,2), ENERGY(285)
	   REAL*4	 EV,EUCONV
	   REAL*4	 ATWT
C
	   INTEGER*4   	 IBUF(500),IZ!,ITEMP,ISWAP	! never referenced
	   INTEGER*4	 I
C 
           CHARACTER*(*) EUNITS
           CHARACTER*(2) SYM
C
	   EQUIVALENCE (IBUF(1),OUTFIL(1,1))
C
	   DATA OUTFIL/570*-5.0/
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  energy units conversion
C
	   EV=EUCONV(E,EUNITS,'EV')
C
C  selection of Henke or Cromer Library
C
           IF ((EV.GE.100.).AND.(EV.LE.1000.)) THEN
C
C             Henke data (file F12C.DAT will be used)
C
csrio
              CALL F1F2(SYM,IZ, OUTFIL, ATWT, ENERGY)
C
	      DO 100 I=1,284
C
C               Linear interpolation
C
                IF((ENERGY(I).LE.EV).AND.(ENERGY(I+1).GE.EV)) THEN
		    C=(EV-ENERGY(I))/(ENERGY(I+1)-ENERGY(I))
	            F1=OUTFIL(I,1)+(OUTFIL(I+1,1)-OUTFIL(I,1))*C
		    F2=OUTFIL(I,2)+(OUTFIL(I+1,2)-OUTFIL(I,2))*C
                ENDIF
C
 100	      CONTINUE
C
           ELSE
C
C          Cromer-Liberman calculations (file CROMER.SYM is used)
C
csrio
csrio	write (*,*) 'calling CROLI_PUB: SYM,IZ,EV,F1,F2',
csrio     1SYM,IZ,EV,F1,F2
              CALL CROLI_PUB(SYM,IZ,EV,F1,F2)
csrio	write (*,*) 'returning CROLI_PUB: SYM,IZ,EV,F1,F2',
csrio     1SYM,IZ,EV,F1,F2
	      R_IZ = FLOAT(IZ)
              F1 = F1
C
           ENDIF
C       
           RETURN           
	   END
C
C------------------------------------------------------------------------------
C Subroutine:    FIND 
C 
C Written by B.Henke
C
C Modified:
C    06/22/88:  E. Ziegler
C    Nov. 90 :  G. Foerstner
C    05/13/91:  M. Krisch
C
C SUBROUTINE FIND is called in S_SCATTER.FOR
C
C Description :  
C
C   FIND returns various constants for an element that is specified by the user 
C
C Usage:         CALL FIND (SYM,IZ,ATWT,RMU,EMF,NAME,Density)
C
C Warning :
C
C   SYM is an CHARACTER*2 variable that contains the chemical symbol of the 
C   element (For an element with a single letter symbol, e.g., "H", use IEL = 
C   'H ', not IEL = 'H'. -1 is returned as the atomic number if the symbol is 
C   not recognized.
C 
C Input Parameters:
C    data file :  "IDATA.DAT"
C    SYM :        CHARACTER STRING
C
C Output Parameters:                                                    
C     IZ          atomic number (INTEGER*2)               
C     ATWT        atomic weight (REAL*4)                  
C     RMU         mu(barns/atom)/mu(cm**2/gram) (REAL*4) 
C     EMF         energy*mu/f2(ev-cm**2/gram) (REAL*4)    
C     NAME        symbol for the element (CHARACTER *(2))      
C     Density     'standard' density (REAL*4)
C
C To speed up access, the OPEN and CLOSE statements are put into the main 
C program. It is necessary to include the REWIND statement in the subroutine.  
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      	   SUBROUTINE FIND(SYM,IZ,ATWT,RMU,EMF,NAME,Density)
C
	   IMPLICIT NONE
C
	   INTEGER*4		I,IERR    !,REC
      	   INTEGER*4 		IZ
C
	   REAL*4		ATWT,RMU,EMF,Density
C
      	   CHARACTER*(2) 	SYM, NAME
	   CHARACTER            JUHU*4
           CHARACTER*1024		FILNAM,datadir!,LINE,MODFIL
	   integer		iblank2											! potee
	COMMON/DATADIR/datadir
C
C      REWIND(3)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C open file with periodic table infos
C
c          FILNAM='IZAR$ESRF:[ESRF.OPTICS.HENKE]IDATA.DAT'
c          FILNAM='/civa/users/b/srio/inpro/data/idata.dat'

cLahey	call getenv ("INPRO_HOME",datadir)
cLahey	if (datadir(1:5).EQ.'     ') then
cLahey	  write (*,*) 'Please define the INPRO_HOME environment'
cLahey	  write (*,*) ' variable to run inpro correctly'
cLahey	  stop
cLahey	endif

c 	filnam = datadir(1:iblank(datadir)) // 'idata.dat'
	
  	FILNAM=datadir(1:iblank2(datadir))//'idata.dat' ! line 273			! potee
 
  	OPEN(33,FILE=FILNAM,STATUS='OLD',ACCESS='SEQUENTIAL',
     &	FORM='FORMATTED')
	 
	write(*,*)'File successfully opened : ',FILNAM					! potee
	 
C
          DO 20 I = 1,94
             READ(33,330) IZ, NAME,ATWT, EMF, RMU, Density
 330	     FORMAT(1X,I4,A2,F11.5,F9.1,F9.2,F10.4)
C	     TYPE *, IZ, NAME,ATWT, EMF, RMU, Density
C
	     IF (SYM .EQ. NAME) THEN 
		IZ = I
		CLOSE(33)
	        RETURN
	     END IF
 20       CONTINUE
C
C  /*SEARCH FAILURE*/
          IZ=-1
          CALL MSGERR('GETF12: Invalid symbol..'//SYM,.TRUE.)

C 2000  CALL MSGERR('FIND: EOF IDATA.DAT',.TRUE.)
C 1101  CALL MSGERR('FIND: ERROR OPENING'//FILNAM,.TRUE.)
	   ierr = 0 ! added by srio
	   JUHU = CHAR(IERR)
 100	   CALL MSGERR('FIND: ERROR NR.IERR=?= '//JUHU,.TRUE.)
	   RETURN
           END

C------------------------------------------------------------------
C Subroutine:   F1F2
C by J. Bloch, 3/29/86, FORTRAN 77
C 
C Modified:
C    by J. Bloch: Symbol names now handled with character variables.
C    by E. Ziegler, 6/22/88
C    by G. Foerstner, Nov. 90
C    by M. Krisch, May 91
C
C Description:
C
C    F1F2 fills an array, "OUTFIL" (dimensioned (285,2)), with f1 and f2 data 
C    for the element with atomic number Z (Z is an integer variable, 1<=Z<=94.
C    If Z is not in this range, the variable ATWT is assigned the value -1 and 
C    OUTFIL is not loaded).
C
C Usage:        CALL F1F2 (SYM,Z,OUTFIL,ATWT,ENERGY)
C
C Input:  
C  Z            atomic number
C
C Output:
C  OUTFIL(I,1)  f1 values 
C  OUTFIL(I,2)  f2 values                       
C  ATWT         atomic weight. -1 is returned if the symbol is not recognized
C                   
C  ENERGY       energy values in eV (REAL*4)                         
C                
C Warning:    for an element with a symbol which is a single letter, e.g. "H",
C             do not use IZ('H'). Use IZ('H ') or I='H',IZ(I), where I is an 
C             integer variable.
C To speed up access, the OPEN and CLOSE statements are put into the main 
C program. It is necessary to include the REWIND statement in the subroutine.  
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      	   SUBROUTINE F1F2(SYM,Z,OUTFIL,ATWT,ENERGY)
C
	   IMPLICIT NONE
C
	   REAL*4	OUTFIL(285,2),ENERGY(285),DUM(285)
	   REAL*4	ATWT
C
	   INTEGER*4	I
      	   INTEGER*4 	Z
C
      	   CHARACTER 	SYM*2, FILNAM*1024
	   character*1024 datadir
	   integer	iblank2													! potee	
	COMMON/DATADIR/datadir

      	   IF (Z .GE. 1 .AND. Z .LE. 94) GO TO 10
     	   ATWT=-1.0
      	   RETURN
C
C  open Henke data table [0.1,1] keV
C
c 10        FILNAM='IZAR$ESRF:[ESRF.OPTICS.HENKE]HENKE.'//SYM
c 10        FILNAM='/civa/users/b/srio/inpro/data/HENKE.'//SYM
cLahey	call getenv ("INPRO_HOME",datadir)
cLahey        if (datadir(1:5).EQ.'     ') then
cLahey          write (*,*) 'Please define the INPRO_HOME environment'
cLahey          write (*,*) ' variable to run inpro correctly'
cLahey          stop
cLahey        endif

c 10	filnam = datadir(1:iblank(datadir)) // 'HENKE.'//SYM

10	FILNAM = datadir(1:iblank2(datadir)) // 'HENKE.'//SYM				! potee

           OPEN(UNIT=2,FILE=FILNAM, ACCESS='SEQUENTIAL',STATUS='OLD',
     &          FORM='FORMATTED', ERR=1001)
csrio     &          FORM='FORMATTED',READONLY, ERR=1001)

	write(*,*)'File successfully opened : ',FILNAM						! potee

C
C      REWIND(2)
C      
C    skip first record
      	   READ (2,*)		

      	   DO 20 I=1,285
              READ (2,*,END=30,ERR=40) ENERGY(I),OUTFIL(I,1),
     &             OUTFIL(I,2),DUM(I)
20         CONTINUE

      	   CLOSE(2)
           RETURN
C
 30    	   CLOSE(2)
      	   CALL MSGERR('F1F2: END-OF-FILE FOR HENKE.'//SYM,.FALSE.)
 40 	   CLOSE(2)
      	   CALL MSGERR('F1F2: ERROR READING FILE HENKE.'//SYM,
     &                 .TRUE.)
 1001      CALL MSGERR('GETF12: Error opening HENKE.'//SYM,.TRUE.)
C      
      	   STOP    
           END
C
C
C Subroutine: CROLI_PUB
C
C Original version: Don T. Cromer, Los Alamos Scientific Lab Report LA-4403
C                   Los Alamos National Laboratory, 
C                   MST-5, Mail Stop G730, Los Alamos, New Mexico 87545,USA.
C
C Modifications:    E. Ziegler, 3/06/90
C		    G. Foerstenr, Nov. 90
C
C Credit: this version takes in account the prgm form given by
C         K.D. Eichhorn, 10/12/87
C         M. H. Krisch, 9/26/89
C         S. Joksch, 9/11/89
C
C
C Description:
C
C    CROLI_PUB calculates the real(f') and the imaginary(f") parts of the 
C    anomalous dispersion correction to the atomic X-ray scattering factor:
C
C                       f = fo + f' + if''
C
C    for all elements from Lithium (Z= 3) through Californium (Z= 98) at 
C    arbitrary X-ray energies.
C
C    Three contributions to the real part, f' are listed separately 
C    in the output:                                                   
C
C    P.E.:   photoelectric contribution
C            The photoelectric contribution, which is obtained by Gaussian 
C            integration as described by Cromer & Liberman, modified to increase
C            the accuracy for energies near to and on the long-wavelength side
C            of an absorption edge (Cromer&Liberman,Acta Cryst.A37,267-268(1981)
C             
C    Eterm:  term dependent on the total energy of the atom
C            (cf Cromer&Liberman,J.Chem. Phys.53,1891 (1970)) 
C	     a relativistic term, dependent on the total energy
C               	[ 5/3 * E(tot) ] / [ m * c**2 ]
C
C    Jensen: term dependent on the xray energy.
C            (cf M.S.Jensen, Physics Letters,74A,41-44(1979))
C 	     a correction, depending on the X-ray energy
C               	-1/2 Z [ h * w ] / [ m * c**2 ]
C
C
C   Method:
C
C     The computations are based on the algorithm of Cromer and Liberman, 
C     J. Chem. Phys. 53, 1891-1898 (1970), and Acta Cryst.18, 17(1965) 
C     & A37, 267 (1981).  
C     The following gives a brief summary of their method:
C
C     The imaginary part, f'', is calculated from atomic absorption 
C     cross-sections SIGMA which are obtained by interpolation from a data 
C     file CROMER.el containing the SIGMA's for element 'el'.  
C     These cross sections were calculated by Cromer & Liberman from 
C     Dirac-Slater relativistic wave-functions.  
C
C     From the Optical Theorem, we have:
C                  f'' = ( m c w ) / ( 4 pi e**2 ) SIGMA 
C    
C     'm' being the electron mass, 'e' its charge, 'c' the speed of 
C     light and 'w' the photon frequency.
C
C
C   Usage: CALL CROLI_PUB (SYM,IZ,E,FP,FPP,NERR)
C
C   INPUT parameters:
C     
C    SYM: atomic symbol
C
C     IZ: atomic number
C
C      E: x-ray energy (in eV) 
C
C   OUTPUT parameters:
C
C     FP:  f' real part of atomic scattering factor
C
C     FPP: f" imaginary part of atomic scattering factor
C    
C  FILES:
C
C    Cross section file: 4 (could also be the same as the input file) 
C
C
C    The cross section file will have a number of orbitals for each atom from 
C    atomic number 3-98. 
C     	
C    The first mx (= 5) records for each orbital have cross sections at mx 
C    values of energy from about 1 to 80 keV approximately equally spaced in 
C    Log(Energy).  
C
C    The next five records contain cross sections at energies selected by the 
C    GAUSS integration scheme.
C    
C    If the 'function type' is 0 (IF = 0), a sixth value is read in for an 
C    energy=1.001*BE
C	
C      BE being the binding energy, the 'function type' given by 
C      Cromer & Liberman, J.Chem. Phys. 53, 1891 (1970)
C	
C    If the xray energy is less than binding energy, function sigma3 will be 
C    used (Cromer & Liberman, Acta Cryst. A37, 267-268 (1981)). 
C
C    Note that the first record on this file is blank (or comment).
C
C    For any given element each record of the atomic cross sections has the 
C    following form:
C      	             ATOM,NJ,SHELL,EW,EX,SIG,BE,IF
C		     in FORMAT(A2,2X, A4, A2,4E15.8,I2)
C	where
C            ATOM :  atomic symbol, 2 chars, used as a check
C              NJ :  orbital sequence number
C           SHELL :  orbital type. 1s1/2 etc., 6 chars
C              EW :  energy in angstroms (wavelength)
C              EX :  energy in keV
C             SIG :  cross section in barns
C              BE :  binding energy (keV)
C              IF :  0, 1, or 2 - flag for function type
C
C  WARNING:
C  --------
C
C     If an X-ray energy is very close to one of the energies used in the GAUSS
C     integration, an 'anomalous' anomalous scattering factor may result. 
C     There is no easy way to get around this problem.
C     A suggestion is to compute several f' values at nearby energies and
C     to draw a smooth curve.  
C     This method should work provided the points do not pass through an edge.  
C     Option 'DD' of =DCBINS= may turn out to be handy in such cases.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

	   SUBROUTINE CROLI_PUB (SYM,MZ,E,RFP,RFPP)
C
C
	   IMPLICIT NONE
C
	   INTEGER*4         ICOUNT,J,MF,MM,MX,NO,NJ,NX,N1
	   INTEGER*4	     NORB(98)
	   INTEGER*4	     NORD,K,M,IF,ISERR
	   INTEGER*4         MZ
C
	   REAL*4	     RFP,RFPP,E
C	   REAL*8	     DCRCOR,RCOR
	   REAL*8            AU,C,C1,PI
	   REAL*8            DCAIP,DCGINT,DCSIG0,DCSIG1!,DCJCOR
	   REAL*8            DCSIG2,DCSIG3
	   REAL*8            BE,BB,CX,EDG,EG(5),EL(14),EW(14)
	   REAL*8            F1,F2,RX,SEDGE,SIG(14),SIGEDG	! 556
	   REAL*8            SG(5),SL(14)
	   REAL*8            XK,ZX,CORR,XW,FP,FPP,T(80)
C
	   CHARACTER         ATOM*2, SHELL*6, SYM*2, FILNAM*1024!,ATNAME*2
	   character*1024	     datadir
	   integer	     iblank2				! potee
	COMMON/DATADIR/datadir
C
C
	   EXTERNAL          DCSIG0, DCSIG1, DCSIG2, DCSIG3
C
	   COMMON /BLKGAU/ CX,BB,SG,RX,SEDGE,ICOUNT
C
C     No. of orbitals for atoms numbers from  1 thru 98
	   DATA NORB /
     1     0,  0,  2,  2,  3,  3,  4,  4,  4,  4,  4,  4,  5,  6,  7,  
     2     7,  7,  7,  7,  7,  7,  7,  8,  9,  9,  9,  9,  9,  9,  9,
     3     9,  9,  9,  9,  9,  9,  9, 12, 12, 13, 13, 14, 14, 14, 14,
     4    14, 14, 14, 14, 14, 14, 14, 14, 14, 17, 17, 17, 18, 18, 18,
     5    18, 18, 18, 19, 19, 19, 19, 19, 19, 19, 19, 20, 20, 20, 21,
     6    21, 21, 21, 21, 22, 22, 23, 23, 24, 24, 24, 24, 24, 24, 24,
     7    24, 24, 24, 24, 24, 24, 24, 24  /

C     Conversion factors ...
C          c :  from cross sections to f's
C         au :  cross sections from barns to a.u.'s
C         c1 :  energies from KeV to a.u.'s
C         mx :  no. of records per orbital
C       nord :  interpolation order in =DCAIP= (2 works best)
C
C      Parameter setups:
C
	   C = 137.0367D0
           AU = 2.80022D+7  
	   C1 =  0.02721D0
	   PI = 3.14159265D0
	   MX = 5
	   NORD = 2	! line 594
C
C      Parameter setups:
C
	   FP = 0.0
	   FPP = 0.0
C
C     Number of orbitals				! 601
	   NO = NORB(INT(MZ + 0.1))
C
C     open cross-section file
c	   FILNAM='IZAR$ESRF:[ESRF.OPTICS.CROMER]
c	   FILNAM='/civa/users/b/srio/inpro/data/
cLahey	call getenv ("INPRO_HOME",datadir)
cLahey        if (datadir(1:5).EQ.'     ') then
cLahey          write (*,*) 'Please define the INPRO_HOME environment'
cLahey          write (*,*) ' variable to run inpro correctly'
cLahey          stop
cLahey        endif

c 	filnam = datadir(1:iblank(datadir)) //'CROMER.'//SYM

 	FILNAM = datadir(1:iblank2(datadir)) //'CROMER.'//SYM	! 616			! potee

           OPEN(UNIT=4, FILE=FILNAM, ACCESS='SEQUENTIAL',
     &STATUS='OLD', FORM='FORMATTED', ERR=1301)
csrio     &STATUS='OLD', FORM='FORMATTED', readonly, ERR=1301)

	write(*,*)'File successfully opened : ',FILNAM						! potee
	 

C
C---- ***** Start loop thru atomic orbitals *****
C
	   READ(4,*) ! Read blank line
C
	   DO 300 J=1,NO
C
C----  ... Read mx+5 energies and cross sections for orbital J
              NX = MX + 5
              DO 30 K=1,NX
                READ(4,1200,ERR=888,END=889) ATOM,NJ,SHELL,XW,
     1                                         EW(K),SIG(K),BE,IF
C 
C               SYM :  atomic symbol, in CAPITAL letters
                IF (ATOM.NE.SYM) GOTO 890
                IF (NJ.NE.J) GO TO 888
 30           CONTINUE
 1200         FORMAT(A2,2X,I4,A6,1P3E15.8,E15.8,I2)

C----  ... Store points for GAUSS integration
              DO 40 K = 1,5
                EG(  K) = EW( K+MX)
                SG(K) = SIG(K+MX) / AU
 40           CONTINUE
C
C
C----  ... Check binding energy .. convert to a.u.'s
              IF (BE.LE.0.0D0) GO TO 888
              BB = BE / C1
C
C----  ... Read cross section at 1.001 * binding energy if IF = 0 
              IF (IF.EQ.0) THEN
                NX = NX + 1
                READ(4,1200,ERR=888,END=889) ATOM,NJ,SHELL,XW,
     1                                        EW(NX),SIG(NX)
                IF (ATOM.NE.SYM) GO TO 890
                EDG    = EW( NX)
                SIGEDG = SIG(NX)
              ENDIF
C
C----  ... Sort cross sections and integration points
              CALL DCSRTD(EW,SIG,NX)
              CALL DCSRTD(EG,SG,5)
C
C----  ... Switch to logs for interpolation
              DO 70 K=1,NX
                EL(K) = DLOG(EW(K))
                IF (SIG(K).EQ.0.0D0) THEN
                     SL(K) = 0.0D0
                ELSE
                     SL(K) = DLOG(SIG(K))
                ENDIF
 70           CONTINUE
C
C
              MF = 0
C
C----       ... XK, RX are energy in keV and au, ZX is logarithm
              XK = 0.001D0 * E
              RX = XK / C1
              ZX = DLOG(XK)
C
C----       ... Interpolated cross section (in a.u.'s)
              CX = 0.0D0
              ISERR = 0
              IF (BE.LE.XK) THEN
                  DO 90 M=1,NX
                     N1 = M
                     IF (SL(M).NE.0.0D0) GO TO 100
 90               CONTINUE
 100              MM = NX - N1 + 1
                  CX = DCAIP(ZX,MM,NORD,EL(N1),SL(N1),T,ISERR)
                  CX = DEXP(CX) / AU
              ENDIF
C
C----       ... contribution to f' from J'th orbital ...
C
              ICOUNT = 6
C
              IF (IF.NE.0 .OR. BE.LT.XK) THEN
                 IF (IF.EQ.0) F1 = DCGINT(DCSIG0)
                 IF (IF.EQ.1) F1 = DCGINT(DCSIG1)
                 IF (IF.EQ.2) F1 = DCGINT(DCSIG2)
              ELSE
                 SEDGE = SIGEDG / AU
                 CX    = 0.0D0
                 F1    = DCGINT(DCSIG3)
                 MF    = 3
              ENDIF
              CORR = 0.0D0
              IF (CX.NE.0.0) CORR =-CX * RX * 0.5D0 * 
     1                                     DLOG((RX+BB)/(RX-BB))
              IF (MF.EQ.3  ) CORR = 0.5D0 * SEDGE * BB**2 * 
     1                                     DLOG((-BB+RX)/(-BB-RX))/RX
              	   F1    = F1 + CORR
              FP = FP + (F1 * C / (2.D0 * PI**2))
C
C----      ... contribution to f'' from J'th orbital 
C
              F2 = 0.0D0
              IF (CX.NE.0.0) F2 = C * CX * RX / (4.D0 * PI)
              FPP = FPP + F2
C
	      RFP = REAL(FP)
	      RFPP = REAL(FPP)
C
 300 	   CONTINUE
C
C---- ***** ... End loop thru atomic orbitals *****

	   CLOSE(4)
      	   RETURN
C
 1301      CLOSE(4)
      	   CALL MSGERR('CROLI_PUB: Error opening CROMER.'//SYM,
     &                .TRUE.)
    	! 741          
C---- Errors
 888 	   CLOSE(4)
     	   CALL MSGERR('Error reading CROMER.'//SYM,.TRUE.)
 889 	   CLOSE(4)
      	   CALL MSGERR('Error end-of-file CROMER.'//SYM, .TRUE.)
 890 	   CLOSE(4)
      	   CALL MSGERR('Wrong atom file for'//SYM,.TRUE.)
      	   STOP
      	   END		! 750


C-----------------------------------------------------------------------
C
	   DOUBLE PRECISION FUNCTION DCGINT(FUN)
C          =====================================
C
C---- Gaussian integration.  FUN is the function to be integrated.
C
C
	   IMPLICIT NONE
	   Integer*4     I
C
	   REAL*8        A(5),FUN,Z(5)
C
C---- Gaussian weights and ordinates (5 points) ...
C
	   DATA A(1) / 0.11846344252810D0 /
	   DATA Z(1) / 0.04691007703067D0 /
	   DATA A(2) / 0.23931433524968D0 /
	   DATA Z(2) / 0.23076534494716D0 /
	   DATA A(3) / 0.28444444444444D0 /
	   DATA Z(3) / 0.50000000000000D0 /
	   DATA A(4) / 0.23931433524968D0 /
	   DATA Z(4) / 0.76923465505284D0 /
	   DATA A(5) / 0.11846344252810D0 /
	   DATA Z(5) / 0.95308992296933D0 /
C
	   DCGINT = 0.0D0
	   DO 100 I = 1,5
              DCGINT = DCGINT + A(I) * FUN(Z(I))
 100 	   CONTINUE
C
	   RETURN
	   END
C
C-----------------------------------------------------------------------
C
	   DOUBLE PRECISION FUNCTION DCSIG0(X)
C          ===================================
C
C---- =DCSIG0= to =DCSIG3= are the functions called by =DCGINT= for 
C     Gaussian integration.  Contents of /BLKGAU/ is (all quantities 
C     being in a.u.'s) ...
C				! 795
C             CX :  cross section at energy RX
C             BB :  binding energy of current orbital
C             SG :  cross section at integration points
C             RX :  X-ray energy
C       COUNT :  counter for interation points
C          SEDGE :  cross section at 0.001 * BE (if any)
C
C
	   IMPLICIT NONE
C
	   INTEGER*4  	ICOUNT
	   REAL*8     	BB,CX,RX,X,SG(5),SEDGE
C
	   COMMON /BLKGAU/ CX, BB, SG, RX, SEDGE, ICOUNT
C
	   ICOUNT = ICOUNT - 1
	   DCSIG0 = SG(ICOUNT) * BB**3 / X**2 / (RX**2 * X**2 - BB**2) - 
     1                   BB * CX * RX**2 / (RX**2 * X**2 - BB**2)
C
	   RETURN
	   END
C
C-----------------------------------------------------------------------
C
	   DOUBLE PRECISION FUNCTION DCSIG1(X)
C          ===================================
C
C
	   IMPLICIT NONE
C
	   INTEGER*4  	ICOUNT
	   REAL*8     	BB,RX,SG(5),X,CX,SEDGE
C
	   COMMON /BLKGAU/ CX, BB, SG, RX, SEDGE, ICOUNT
C
	   ICOUNT = ICOUNT - 1
	   DCSIG1 = 0.5*BB**3*SG(ICOUNT)/(DSQRT(X)*(RX**2*X**2-BB**2*X))
C
	   RETURN
	   END
C
C-----------------------------------------------------------------------
C
	   DOUBLE PRECISION FUNCTION DCSIG2(X)
C          ===================================
C
	   IMPLICIT NONE
C
C
	   INTEGER*4  	ICOUNT
	   REAL*8     	BB,CX,DENOM,RX,SG(5),X,SEDGE
C
	   COMMON /BLKGAU/ CX, BB, SG,RX, SEDGE, ICOUNT
C
	   ICOUNT = ICOUNT - 1
	   DENOM  = X**3 * RX**2 - BB**2 / X
	   DCSIG2 = 2.0 * BB * SG(ICOUNT) * BB**2 / X**4 / DENOM -
     1                        2.0 * BB * CX * RX**2 / DENOM
C
	   RETURN
	   END
C
C-----------------------------------------------------------------------
C
	   DOUBLE PRECISION FUNCTION DCSIG3(X)
C          ===================================
C
C
	   IMPLICIT NONE
C
	   INTEGER*4  ICOUNT
	   REAL*8     BB,RX,SEDGE,SG(5),X,CX
C
	   COMMON /BLKGAU/ CX,BB,SG,RX,SEDGE,ICOUNT
C
	   ICOUNT = ICOUNT - 1
	   DCSIG3 = 
     $ BB**3*(SG(ICOUNT)-SEDGE*X**2)/(X**2*(X**2*RX**2-BB**2))
C
	   RETURN
	   END
C
C-----------------------------------------------------------------------
C
	   DOUBLE PRECISION FUNCTION DCAIP(XBAR,IN,IM,X,Y,T,NERR)
C          ======================================================
C
C---- AITKEN repeated interpolation scheme.
C
C        XBAR :  abscissa at which interpolation is desired
C          IN :  no. of values in table X,Y ...
C		  - if IN > 0, check ordering of X(I)
C		  - if IN < 0, skip preceeding test
C          IM :  degree of approximating polynomial (IM < IN)
C           X :  vector of IABS(IN) values of abscissa
C           Y :  vector of IABS(IN) values of ordinate
C           T :  temporary storage vector of 4*(IM+1) locations
C
C---- The error flag NERR is ...
C
C        NERR :  non-zero to signal error condition ...
C                = 1 :  order IM is too large, reset to IN-1
C                = 2 :  IN.LT.2  ... YBAR returned as Y(1)
C                = 3 :  X(I) not sequenced properly
C
C
C
	   IMPLICIT NONE
C
	   INTEGER*4      IM,IN,J,JJ,K,KK,M,N,NERR,I
           REAL*8         DXBAR,S,T(80),X(9),XBAR,Y(9),Z
C
C---- Initialize ...
          NERR  = 0
          N     = IABS(IN)
          M     = IM
          DXBAR = XBAR
C
C---- Check order of interpolation
          IF (M.GE.N) THEN
             NERR = 1
             M = N - 1
          ENDIF
C
          K = N - 1
C
C---- Check no. of points
          IF (N.LT.2) THEN
             NERR  = 2
             DCAIP = Y(1)
             RETURN
          ENDIF
C
          S = X(2) - X(1)
C
          IF (IN.LT.0) GO TO 30
          IF ( N.EQ.2) GO TO 30
C
C---- Check if order is monotonic ...
          DO 20 I=3,N
             Z = S * (X(I) - X(I-1))
             IF (Z.LE.0.0D0) THEN
                NERR  = 3
                DCAIP = Y(1)
                RETURN
             ENDIF
 20       CONTINUE
C
C---- Get nearest point ... for increasing order
 30       IF (S.LT.0.0D0) GO TO 50
          DO 40 J=1,N
             IF (XBAR.LE.X(J)) GO TO 70
 40       CONTINUE
          J = N
          GO TO 70
C
C----  ... for decreasing order
 50       DO 60 J=1,N
             IF (XBAR.GE.X(J)) GO TO 70
 60 	  CONTINUE
          J = N
C
C---- Interpolation (N'th order) 
 70       K = M
          M = M + 1
          J = J - M/2
          J = MIN0( MAX0(J,1), N-K )
          DO 80 I=J,J+K
             KK      = I - J + 1
             T(KK)   = Y(I)
             T(KK+M) = X(I) - DXBAR
 80       CONTINUE
          DO 100 I=1,K
             DO 90 JJ=I+1,M
               T(JJ) = (T(I)*T(JJ+M)-T(JJ)*T(I+M))/(X(JJ+J-1)-X(I+J-1))
 90          CONTINUE
 100      CONTINUE
C
C---- At this point we have it ...
C
          DCAIP = T(M)
C
          RETURN
          END
C
C-----------------------------------------------------------------------
C
	  SUBROUTINE DCSRTD(A,B,N)
C         ========================
C
C---- Sort cross sections and points for GAUSS integration.
C
C         A :  array to be sorted (increasing order)
C         B :  reordered into same sequence as A
C         N :  no. of data points in A and B
C
C
	   IMPLICIT NONE
C
           INTEGER*4        N,I,J
	   REAL*8           A(1),B(1),T
C
	   DO 200 I=1,N-1
              DO 100 J=I+1,N
                IF (A(J).LT.A(I)) THEN
                     T    = A(I)
                     A(I) = A(J)
                     A(J) = T
                     T    = B(I)
                     B(I) = B(J)
                     B(J) = T
                ENDIF
 100          CONTINUE
 200       CONTINUE
C
	   RETURN
	   END
C
C-*-*-*-*-*-*-*-*-*-*-*-*  End of Source  *-*-*-*-*-*-*-*-*-*-*-*-*-*-*
C
C Function:     EUCONV  
C by J. Bloch, 3/12/82, Fortran 77  
C   
C Description:  
C
C    EUCONV converts photon energy from one energy or wavelength unit system to
C    another. 
C   
C Usage:      X = EUCONV(Y,YUNITS,XUNITS) 
C   
C Parameters:
C
C    Y :      Real photon energy or wavelength
C   
C    YUNITS : Character string nameing units for Y   
C   
C    XUNITS : Character string nameing the units which Y will be converted to. 
C   
C P.S.:  Allowed Units for XUNITS or YUNITS are: 
c        'ANGSTROMS','NANOMETERS' or 'NM', 'MICRONS', 'CM' or 'CENTIMETERS',  
C          'M' or 'METERS', 'EV' or 'ELECTRON-VOLTS','KEV','MEV'
C   

      	   FUNCTION EUCONV(Y,YUNITS,XUNITS)  
C
	   IMPLICIT NONE
C
           CHARACTER*(*) YUNITS,XUNITS
C
	   REAL*4       Y,A,B,EUCONV
C
           INTEGER*4     CHANPRT
           COMMON /CHAN/  CHANPRT
C
           IF(YUNITS.EQ.XUNITS) THEN
                EUCONV=Y
                RETURN  
           END IF   
C
           A=0.0
           B=0.0
C
           IF(YUNITS.EQ.'ANGSTROMS') A=Y
           IF((YUNITS.EQ.'NANOMETERS').OR.(YUNITS.EQ.'NM')) A=10.*Y
           IF(YUNITS.EQ.'MICRONS') A=Y*1.0E+4
           IF(YUNITS.EQ.'CM') A=Y*1.0E+8
           IF(YUNITS.EQ.'CENTIMETERS') A=Y*1.0E+8
           IF((YUNITS.EQ.'M').OR.(YUNITS.EQ.'METERS')) A=Y*1.E10
           IF((YUNITS.EQ.'EV').OR.(YUNITS.EQ.'ELECTRON-VOLTS'))
     &             A=12398.54/Y 
           IF(YUNITS.EQ.'KEV') A=12.39854/Y
           IF(YUNITS.EQ.'MEV') A=1.239854E-2/Y
C
C
           IF(XUNITS.EQ.'ANGSTROMS') B=A
           IF((XUNITS.EQ.'NANOMETERS').OR.(XUNITS.EQ.'NM')) B=A*.1
           IF(XUNITS.EQ.'MICRONS') B= A*1.0E-4
           IF(XUNITS.EQ.'CM') B=A*1.0E-8
           IF(XUNITS.EQ.'CENTIMETERS') B=A*1.0E-8
           IF((XUNITS.EQ.'METERS').OR.(XUNITS.EQ.'M')) B=A*1.0E-10
           IF((XUNITS.EQ.'EV').OR.(XUNITS.EQ.'ELECTRON-VOLTS'))
     &              B=12398.54/A
           IF(XUNITS.EQ.'KEV') B=12.39854/A
           IF(XUNITS.EQ.'MEV') B=1.239854E-2/A
C
           IF(A.EQ.0.0) THEN
c             WRITE(*,*) 'EUCONV: ILLEGAL INPUT ENERGY UNIT: ',YUNITS
             STOP
           END IF   
C
           IF(B.EQ.0.0) THEN
c              WRITE(9,*) 'EUCONV: ILLEGAL OUTPUT ENERGY UNITS: ',XUNITS  
             STOP
           END IF   
C
C
           EUCONV=B 
C
C
           RETURN   
      END   
C
C
C Subroutine:	MSGERR
C by J. Bloch, 3/28/86, Fortran 77
C
C Description:	
C
C   MSGERR is an error reporting routine. It can also stop the program if 
C   desired.
C
C Usage:	CALL MSGERR(MSG,CONT)
C
C INPUT Parameters:	
C
C   MSG	:       CHARACTER string error message.
C
C   CONT :      LOGICAL variable, .TRUE. if program should be stopped.
C

		SUBROUTINE MSGERR(MSG,CONT)
C
	 	IMPLICIT NONE
C
		CHARACTER*(*) MSG
		LOGICAL CONT
		WRITE(*,*) MSG

		IF(CONT) THEN
C			TYPE *,'Hit <return> to continue...'
C			PAUSE
			READ *
			STOP
                ELSE
			write (*,*)'Hit <return> to continue...'
			READ *
C                        PAUSE
                        RETURN
   		ENDIF

		RETURN
     	END
