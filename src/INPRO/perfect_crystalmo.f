C    
	SUBROUTINE PERFECT_CRYSTALSMO(THET,LAMBDA,I,ZS,ZP)
C****************************************************************
C* calculates complex amplitudes (reflected or transmitted)	*
C* from a perfect crystal in either the Bragg or Laue geometry.	*
C* uses the first expression in Zachariassen 3.130 and 3.137	*
C* crystals data are transmitted through the common/DIFFRACTION/*
C* constants are in the common/CONSTANTS/			*
C* **Modified by Carlos Giles for use with the INPRO program    *
C* (4/12/92)							*
C****************************************************************
C
	IMPLICIT NONE
	INTEGER I,IG
	REAL*8 THET,LAMBDA
	REAL*8 SINTB,AKPA,RRP,RIP
	REAL*8 ECONVEL,PI,R0,SEC_RAD
C	REAL*8 lnom,d,b,k,gamma0,gammah,T0
	REAL*8 lnom,dspace,bfac,facp,gamma0,gammah,T0
C	REAL*8 OFFSET,COSA,SINA		! never used
	REAL*8 alpha

	COMPLEX*16 ZS,ZP
	COMPLEX*16 Z,ZQ,Z1,CARG1,CARG2,CDEL1,CDEL2,CP1,CP2,
     1	CX1,CX2,C1,C2,CDEN
	COMPLEX*16 IZ
	COMPLEX*16 psi0,psih,psihb
C
	COMMON/CONSTANTS/ECONVEL,PI,R0,SEC_RAD,IZ

        COMMON/DIFFRACTION/psi0,psih,psihb,dspace,bfac,facp,
     1                  gamma0,gammah,t0,alpha,lnom,IG



C=============================================================
	SINTB=LAMBDA/(2.D0*dspace)	! angle de Bragg pour lambda
	Z=(1.D0-bfac)*psi0+4.D0*bfac*SINTB*(SINTB-DSIN(THET))
	Z=Z/2.D0		! z de Zachariassen
C
C==============		! SIGMA polarization
	ZQ=bfac*psih*psihb
	Z1=ZQ+(Z*Z)
	Z1=SQRT(Z1)
C==============
	AKPA=PI/(LAMBDA*gamma0)
	CARG1=-Z+Z1
	CARG2=-Z-Z1
	CDEL1=AKPA*(psi0+CARG1)	
	CDEL2=AKPA*(psi0+CARG2)	
C
	CP1=-IZ*CDEL1*T0
	CP2=-IZ*CDEL2*T0
	CX1=CARG1/psih
	CX2=CARG2/psih
C==============		! stops under- and over-flows
	RRP=70.D0
	IF(ABS(DREAL(CP1)).GT.RRP) THEN
		RIP=DIMAG(CP1)
		IF(DREAL(CP1).GT.RRP) C1=EXP(CMPLX(RRP,RIP))
		IF(DREAL(CP1).LT.-RRP) C1=EXP(CMPLX(-RRP,RIP))
	ELSE IF (ABS(DREAL(CP1)).LE.RRP) THEN
		C1=EXP(CP1)
	ENDIF
	IF(ABS(DREAL(CP2)).GT.RRP) THEN
		RIP=DIMAG(CP2)
		IF(DREAL(CP2).GT.RRP) C2=EXP(CMPLX(RRP,RIP))
		IF(DREAL(CP2).LT.-RRP) C2=EXP(CMPLX(-RRP,RIP))
	ELSE IF (ABS(DREAL(CP2)).LE.RRP) THEN
		C2=EXP(CP2)
	ENDIF
C======================
	IF (IG.GT.0) THEN		! reflected beam
		IF(IG.EQ.1) THEN	! Bragg case
			CDEN=C2*CX2-C1*CX1
			ZS=CX1*CX2*(C1-C2)/CDEN/DSQRT(DABS(bfac))
		ELSE IF (IG.EQ.2) THEN	! Laue case
			CDEN=CX2-CX1
			ZS=CX1*CX2*(C1-C2)/CDEN/DSQRT(DABS(bfac))
		ENDIF
C======================
C srio@esrf.eu debugging 2012/09/27
C from John Sutter john.sutter@diamond.ac.uk
C
C Dear Manuel,
C Thank you for checking this. Yes, I was trying to say that we shouldn’t 
C divide by the asymmetry factor when dealing with transmitted or 
C forward-diffracted beam, because the cross-section of the beam is 
C not changed. Only the reflected beam intensity needs to be divided by the 
C asymmetry factor. So, you should remove the text “/DSQRT(DABS(bfac))” 
C in lines 14 and 17 of your code segment, but keep it in lines 5 and 8.
C 
C If you make this correction, then the sum of the reflected power and 
C the transmitted/forward-diffracted power should be less than 1 at every 
C point. Furthermore, the transmitted intensity far away from the Bragg 
C reflection should approach the value given by ordinary photoelectric 
C absorption and Compton scattering.
C 
C  
C
	ELSE IF (IG.LT.0) THEN		! transmitted beam
		IF(IG.EQ.-1) THEN	! Bragg case
			CDEN=C2*CX2-C1*CX1
C			ZS=C1*C2*(CX2-CX1)/CDEN/DSQRT(DABS(bfac))
			ZS=C1*C2*(CX2-CX1)/CDEN
		ELSE IF (IG.EQ.-2) THEN	! Laue case
			CDEN=CX2-CX1
C			ZS=(CX2*C1-CX1*C2)/CDEN/DSQRT(DABS(bfac))
			ZS=(CX2*C1-CX1*C2)/CDEN
		ENDIF
	ENDIF
C============================================================
C			! PI polarization
	ZQ=bfac*facp*facp*psih*psihb
	Z1=ZQ+(Z*Z)			! q+z2 dans Zacha.
	Z1=SQRT(Z1)
C==============
	AKPA=PI/(LAMBDA*gamma0)
	CARG1=-Z+Z1
	CARG2=-Z-Z1
	CDEL1=AKPA*(psi0+CARG1)	
	CDEL2=AKPA*(psi0+CARG2)	
C
	CP1=-IZ*CDEL1*T0
	CP2=-IZ*CDEL2*T0
	CX1=CARG1/(facp*psih)
	CX2=CARG2/(facp*psih)
C==============		! stops under- and over-flows
	RRP=70.D0
	IF(ABS(DREAL(CP1)).GT.RRP) THEN
		RIP=DIMAG(CP1)
		IF(DREAL(CP1).GT.RRP) C1=EXP(CMPLX(RRP,RIP))
		IF(DREAL(CP1).LT.-RRP) C1=EXP(CMPLX(-RRP,RIP))
	ELSE IF (ABS(DREAL(CP1)).LE.RRP) THEN
		C1=EXP(CP1)
	ENDIF
	IF(ABS(DREAL(CP2)).GT.RRP) THEN
		RIP=DIMAG(CP2)
		IF(DREAL(CP2).GT.RRP) C2=EXP(CMPLX(RRP,RIP))
		IF(DREAL(CP2).LT.-RRP) C2=EXP(CMPLX(-RRP,RIP))
	ELSE IF (ABS(DREAL(CP2)).LE.RRP) THEN
		C2=EXP(CP2)
	ENDIF
C======================
	IF (IG.GT.0) THEN		! reflected beam
		IF(IG.EQ.1) THEN	! Bragg case
			CDEN=C2*CX2-C1*CX1
			ZP=CX1*CX2*(C1-C2)/CDEN/DSQRT(DABS(bfac))
		ELSE IF (IG.EQ.2) THEN	! Laue case
			CDEN=CX2-CX1
			ZP=CX1*CX2*(C1-C2)/CDEN/DSQRT(DABS(bfac))
		ENDIF

C======================
	ELSE IF (IG.LT.0) THEN		! transmitted beam
		IF(IG.EQ.-1) THEN	! Bragg case
			CDEN=C2*CX2-C1*CX1
C			ZP=C1*C2*(CX2-CX1)/CDEN/DSQRT(DABS(bfac))
			ZP=C1*C2*(CX2-CX1)/CDEN
		ELSE IF (IG.EQ.-2) THEN	! Laue case
			CDEN=CX2-CX1
C			ZP=(CX2*C1-CX1*C2)/CDEN/DSQRT(DABS(bfac))
			ZP=(CX2*C1-CX1*C2)/CDEN
		ENDIF
	ENDIF
C================================================================
	END
C
