	REAL*8 FUNCTION GY(N,Y)
C+
C 
C FUNCTIONAL DESCRIPTION:	
C 					 inf.
C  Function to calculate the integral = I  K5/3 (t) dt multiplied by y**n.
C					 y
C AUTHORS: 
C 
C  Roger J. Dejus
C 
C  The Advanced Photon Source
C  Experimental Facilities Division
C  Argonne National Laboratory
C 
C CREATION DATE: 
C 
C  17-DEC-1991
C 
C FUNCTION VALUE:
C 		inf.
C  GY = y**n x I  K5/3 (t) dt
C		y
C 
C FORMAL PARAMETERS:
C  
C  Input arguments:
C  N		Value of exponent
C  Y		Argument
C  
C DESIGN ISSUES:
C  
C  From: V. O. Kostroun, Nucl. Instr. and Methods, 172 (1980) 371.
C  Good for about six decimal figures for Z < 4 for DT = 0.5, and
C  and up to Z = 40 for DT = 0.1.  NY = 5/3.
C 
C KEYWORDS:
C
C  Modified Bessel function, Integration
C
C [optional_header_tags]...
C
C MODIFICATION HISTORY:
C
C        Date     | Name  | Description
C ----------------+-------+-----------------------------------------------------
C 29-JAN-1995     | RJD   | Based on routine by V. O. Kostroun. For general 
C		  |	  | value of NY, use function IKNY.
C ----------------+-------+-----------------------------------------------------
C [change_entry]
C-
C  Labeled constants:
	REAL*8		DT,NY,EPS,ZERO,ONE,HALF
	PARAMETER	(DT=0.10D0,NY=5.0D0/3.0D0)
	PARAMETER	(EPS=1.0D-8,ZERO=0.0D0,ONE=1.0D0,HALF=0.5D0)
 
C  Declarations of scalars:
	INTEGER*4	N
	REAL*8		Y,P,TERM,SUM

	P   = ONE
	SUM = ZERO

C  Integrand: exp(-y*cosh(t)) * cosh(ny*t) / cosh(t)
	TERM = EXP(-Y*COSH(P*DT)) * COSH(NY*P*DT) / COSH(P*DT)
	DO WHILE (TERM .GT. EPS*SUM)
	    SUM  = SUM+TERM
	    P    =   P+ONE
	    TERM = EXP(-Y*COSH(P*DT)) * COSH(NY*P*DT) / COSH(P*DT)
	END DO ! WHILE

	GY = Y**N*DT*(SUM+HALF*EXP(-Y))
	RETURN
	END ! GY
