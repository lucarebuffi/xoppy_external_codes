	SUBROUTINE HUNT(X,N,XP,JLO)
!+
! 
! FUNCTIONAL DESCRIPTION:	
! 
!  Subroutine to find the index JLO in an ordered array X with N elements.
!  The index JLO satisfies: X(JLO) < XP < X(JLO+1).
! 
! AUTHORS: 
! 
!  Roger J. Dejus
! 
!  The Advanced Photon Source
!  Experimental Facilities Division
!  Argonne National Laboratory
! 
! CREATION DATE: 
! 
!  21-NOV-1990
! 
! FORMAL PARAMETERS:
!  
!  Input arguments:
!  X			Array of ordered values (increasing or decreasing)
!  N			Number of elements in X
!  XP			Value XP is searched in X
!  JLO	  		Initial guess of JLO (also output; see below)
!
!  Output arguments:
!  JLO	  		Final value of JLO; When JLO = 0 or JLO = N the value
!			of XP is outside the range of array X
!
! DESIGN ISSUES:
!  
!  Adapted from "Numerical Recipes" by W.H. Press, B.P. Flannery, S.A. Teukolsky
!  and W.T. Vetterling, Cambridge University Press (1988), p. 91.
!  
! KEYWORDS:
!  
!  Search, ordered table
!  
! [optional_header_tags]...
! 
! MODIFICATION HISTORY:
! 
!	 Date     | Name  | Description
! ----------------+-------+-----------------------------------------------------
! 24-NOV-2014     | RJD   | Converted to Fortran 90.
!                 |       | Changed all comments "C" to "!".
! ----------------+-------+-----------------------------------------------------
! [change_entry]
!-

!  Declarations of scalars:
	LOGICAL*4	ASCEND
	INTEGER*4	N,JLO,JHI,JM,INC
	REAL*8		XP

!  Declarations of arrays:
	REAL*8		X(N)

	ASCEND = X(N) .GT. X(1) ! Indicates whether increasing or decreasing
				! values

	IF (JLO .LT. 1 .OR. JLO .GT. N) THEN
	    JLO = 0
	    JHI = N+1
	    GOTO 10
	END IF ! JLO

	INC = 1 ! Set hunting increment
	IF (XP .GE. X(JLO) .EQV. ASCEND) THEN ! Hunt up
20	    JHI = JLO+INC
	    IF (JHI .GT. N) THEN
		JHI = N+1
	    ELSE IF (XP .GE. X(JHI) .EQV. ASCEND) THEN
		JLO = JHI
		INC = INC+INC
		GOTO 20
	    END IF ! JHI
	ELSE ! Hunt down
	    JHI = JLO
30	    JLO = JHI-INC
	    IF (JLO .LT. 1) THEN
		JLO = 0
	    ELSE IF (XP .LT. X(JLO) .EQV. ASCEND) THEN
		JHI = JLO
		INC = INC+INC
		GOTO 30
	    END IF ! JLO
	END IF ! XP

!  The value XP is now bracketed by X(JLO) < XP < X(JHI)
!  Begin final bisection phase
10	CONTINUE
	DO WHILE (JHI-JLO .GT. 1)
	    JM = (JLO+JHI)/2
	    IF (XP .GT. X(JM) .EQV. ASCEND) THEN
		JLO = JM
	    ELSE
		JHI = JM
	    END IF
	END DO ! WHILE

	RETURN
	END ! HUNT
