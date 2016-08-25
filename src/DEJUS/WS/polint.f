	SUBROUTINE POLINT(XA,YA,N,X1,Y1,N1,DY1)
C+
C 
C FUNCTIONAL DESCRIPTION:	
C 
C  Subroutine to make a polynomial interpolation/extrapolation of degree N-1
C  to obtain the values Y1 for given abscissas X1.
C 
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
C  20-NOV-1990
C 
C FORMAL PARAMETERS:
C   
C  Input arguments:
C  XA			Array of abscissas
C  YA			Array of ordinates
C  N			Number of points (give polynomial of degree N-1)
C  X1			Value of abscissas for which the interpolated/extra-
C			polated values are sought
C  N1			Number of elements in X1, Y1, and DY1
C
C  Output arguments:
C  Y1			Interpolated/extrapolated values
C  DY1	  		Error estimate of Y1
C
C DESIGN ISSUES:
C  
C  Adapted from "Numerical Recipes" by W.H. Press, B.P. Flannery, S.A. Teukolsky
C  and W.T. Vetterling, Cambridge University Press (1988), p. 82.
C  
C KEYWORDS:
C  
C  Lagrange, polynomial interpolation, polynomial extrapolation
C  
C [optional_header_tags]...
C 
C MODIFICATION HISTORY:
C 
C	 Date     | Name  | Description
C ----------------+-------+-----------------------------------------------------
C 14-JAN-1995     | RJD   | Replaced call to LIB$STOP by STOP statement
C ----------------+-------+-----------------------------------------------------
C [change_entry]
C-
C  Size parameters:
	INTEGER*4	N_SZ
	PARAMETER	(N_SZ=10)

C  Declarations of scalars:
	INTEGER*4	N,N1,I,J,JX,M
	REAL*8		XP,Y,DY,DIFS,DIFP,DEN
	REAL*8		H1,H2,W

C  Declarations of arrays:
	REAL*8		XA(N),YA(N),X1(N1),Y1(N1),DY1(N1)
	REAL*8 		C(N_SZ),D(N_SZ)

C  Labeled constants:
	REAL*8		ZERO
	PARAMETER	(ZERO=0.0D0)
C	INTEGER*4	SS$_ABORT
C	PARAMETER	(SS$_ABORT = '0000002C'X)

	DO I=1,N1 ! Loop over all X1 for which values are sought
	    JX = 1
	    XP = X1(I)
	    DIFS = ABS(XP-XA(1))
	    DO J=1,N
		DIFP = ABS(XP-XA(J))
		IF (DIFP .LT. DIFS) THEN
		    JX   = J
		    DIFS = DIFP
		END IF		
		C(J) = YA(J)
		D(J) = YA(J)
	    END DO ! J

	    Y  = YA(JX) ! Initial approximation
	    JX = JX-1	! Neville's recursive algorithm is being used

	    DO M=1,N-1
		DO J=1,N-M
		    H1  = XA(J)  -XP
		    H2  = XA(J+M)-XP
		    W   = C (J+1)-D(J)
		    DEN = H1-H2
		    IF (DEN .EQ. ZERO) THEN
			JX = J+M
			PRINT *
			PRINT 200,'&POLINT-F-INTERR, Interpolation error'
			PRINT 210,'- XA(',J,' ) and XA(',JX,' ) = ',XA(J)
C			CALL LIB$STOP(%VAL(SS$_ABORT)) ! Abort execution
			STOP 'Execution aborted in subroutine POLINT'
		    END IF ! DEN
		    DEN  = W /DEN
		    C(J) = H1*DEN
		    D(J) = H2*DEN
		END DO ! J

		IF (2*JX .LT. N-M) THEN
		    DY = C(JX+1) ! Take most straight line through the table
		ELSE		 ! to its apex
		    DY = D(JX)
		    JX = JX-1
		END IF
		Y = Y+DY	 ! Update Y
	    END DO ! M

	    Y1 (I) = Y
	    DY1(I) = DY
	END DO ! I

200	FORMAT(' ',A)
210	FORMAT(' ',A,I5,A,I5,A,1PG15.7E2)

	END ! POLINT
