CL# 1 "vector.F"
C +++
C
C Source: src/lib/vector.F
C
C ----------------------------------------------
C                SHADOW
C      Center for X-ray Lithography
C     University of Wisconsin-Madison
C  3731 Schneider Dr., Stoughton, WI, 53589
C ----------------------------------------------
C
C Log: vector.F
C Revision 1.9  1991/07/06  19:56:48  khan
C Grenoble and after. Minor changes
C
C Revision 1.8  91/04/05  13:54:40  cwelnak
C changed quotes on #include
C
C Revision 1.7  91/03/25  10:57:35  cwelnak
C SUN version - INC to #inc
C
C Revision 1.6  90/11/17  12:59:35  khan
C Numbers < 1.0E-31 is now zero'd out explicitly.
C
C Revision 1.5  90/11/16  11:13:03  khan
C Fix in CROSS. VRES(1) ends up being xx.xxE-41 when it should be zero.
C VMS does it fine, but in 1 explicit checking is needed and zero it.
C
C Revision 1.4  90/11/13  14:05:08  khan
C Cleanup and SAVE statements
C
C Revision 1.3  90/07/20  22:05:26  khan
C put #if 1 ... to make it work also on vms
C
C Revision 1.2  90/07/15  15:31:11  khan
C All public include files (common.blk, etc) are now in ./../../include dir.
C
C Revision 1.1  90/07/10  14:57:39  khan
C Initial revision
C
C
C ---




CL# 1 "/users/b/srio/shadow-2.2.0/src/include/header.txt"
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                         C
C     *-------------------------------------------------------------*     C
C     |                                                             |     C
C     | SHADOW         Vs. 2.0  -- May 1993                         |     C
C     |                                                             |     C
C     | (c) F.Cerrina, 1984 - All rights reserved -                 |     C
C     | (c) F.Cerrina, 1985 - All rights reserved -                 |     C
C     | (c) F.Cerrina, 1986 - All rights reserved -                 |     C
C     | (c) F.Cerrina, 1987 - All rights reserved -                 |     C
C     | (c) F.Cerrina, 1988 - All rights reserved -                 |     C
C     | (c) F.Cerrina, 1989 - All rights reserved -                 |     C
C     | (c) F.Cerrina, 1990 - All rights reserved -                 |     C
C     | (c) F.Cerrina, 1991 - All rights reserved -                 |     C
C     | (c) F.Cerrina, 1992 - All rights reserved -                 |     C
C     | (c) F.Cerrina, 1993 - All rights reserved -                 |     C
C     | Center for X-ray Lithography /                              |     C
C     | Department of Electrical and Computer Engineering           |     C
C     | University of Wisconsin - Madison                           |     C
C     | 608 - 877 - 2400   //   263 - 4955                          |     C
C     *-------------------------------------------------------------*     C
C                                                                         C
C     The following software is part of the SHADOW package. It cannot     C
C     be copied,  transmitted, reproduced or otherwise made available     C
C     to  any third party without  the express consent of the author.     C
C                                                                         C
C     The following  software  is  a  working copy  and  as  such not     C
C     error-free. It has been extensively debugged but as the code is     C
C     extended  to  include new applications, other bugs are bound to     C
C     appear.                                                             C
C                                                                         C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CL# 47 "vector.F"




C
C----------------------------------------------------------
C
C       scalar multiplication
C
C
     	SUBROUTINE	SCALAR	( V1,ARG,V2)

     	IMPLICIT	REAL*8		(A-H,O-Z)

     	DIMENSION	V1(3),V2(3)

     	V2(1)	=   V1(1)*ARG
     	V2(2)	=   V1(2)*ARG
     	V2(3)	=   V1(3)*ARG


C
C If the numbers are *very* small, zero them out. This is not done
C for VMS, since it seems to work out fine. Why mess with something
C that already works.
C
	IF (ABS(V2(1)).LT.1.0E-31) V2(1) = 0.0D0
	IF (ABS(V2(2)).LT.1.0E-31) V2(2) = 0.0D0
	IF (ABS(V2(3)).LT.1.0E-31) V2(3) = 0.0D0


     	END
C
C------------------------------------------------------------
C
C       scalar product
C
C
	SUBROUTINE DOT (V1,V2,RES)
	IMPLICIT 	REAL*8 		(A-H,O-Z)
	DIMENSION 	V1(3),V2(3)

	RES = V1(1)*V2(1) + V1(2)*V2(2) + V1(3)*V2(3)


C
C If the numbers are *very* small, zero them out. This is not done
C for VMS, since it seems to work out fine. Why mess with something
C that already works.
C
	IF (ABS(RES).LT.1.0E-31) RES = 0.0D0

     	RETURN
	END
C
C-----------------------------------------------------------
C
C   	vector product :    vres = v1 x v2
C
C
	SUBROUTINE CROSS (V1,V2,VRES)
	IMPLICIT 	REAL*8 		(A-H,O-Z)
C

c
c This causes problems with F77 drivers, since can't use -I directive.
c so I'll use the standard cpp directive instead.
c
c	INCLUDE         './../include/warning.blk'
c
c

CL# 1 "/users/b/srio/shadow-2.2.0/src/include/warning.blk"
C +++
C
C Source: src/include/warning.blk
C
C ----------------------------------------------
C                SHADOW
C      Center for X-ray Lithography
C     University of Wisconsin-Madison
C  3731 Schneider Dr., Stoughton, WI, 53589
C ----------------------------------------------
C
C Log: warning.blk
C Revision 1.2  1991/07/06  20:10:25  khan
C Grenoble and after. Minor changes
C
C Revision 1.1  90/07/10  14:58:05  khan
C Initial revision
C
C
C ---

C+++
C	WARNING.BLK
C
C	Purpose		Contains an error message and a flag value
C
C---
     	CHARACTER*80	M_WARNING
     	COMMON	/WARNING	/M_FLAG,M_WARNING
CL# 119 "vector.F"




     	COMMON	/ERRFLAG	/IFLAG
	DIMENSION 	V1(3),V2(3),VRES (3)
     	M_FLAG	 = 0
	VRES(1)  =     V1(2)*V2(3) - V1(3)*V2(2)
	VRES(2)  = - ( V1(1)*V2(3) - V1(3)*V2(1) )
	VRES(3)  =     V1(1)*V2(2) - V1(2)*V2(1)


C
C If the numbers are *very* small, zero them out. This is not done
C for VMS, since it seems to work out fine. Why mess with something
C that already works.
C
	IF (ABS(VRES(1)).LT.1.0E-31) VRES(1) = 0.0D0
	IF (ABS(VRES(2)).LT.1.0E-31) VRES(2) = 0.0D0
	IF (ABS(VRES(3)).LT.1.0E-31) VRES(3) = 0.0D0

     	TTEST  =  VRES(1)*VRES(1) +
     $		  VRES(2)*VRES(2) +
     $		  VRES(3)*VRES(3)
     	IF (TTEST.LT.1.0E-31) THEN
     	  M_FLAG = 1
     	  M_WARNING = 'Error in CROSS: product is zero.'
     	END IF
	END
C
C-----------------------------------------------------------
C
C       vector normalization
C
C
	SUBROUTINE NORM (V1,V2)
	IMPLICIT 	REAL*8 		(A-H,O-Z)
	DIMENSION 	V1(3),V2(3)
C
	RNORM = V1(1)**2 + V1(2)**2 + V1(3)**2
	RNORM = SQRT(RNORM)


C
C If the numbers are *very* small, zero them out. This is not done
C for VMS, since it seems to work out fine. Why mess with something
C that already works.
C
	IF (ABS(RNORM).LT.1.0E-31) RNORM = 0.0D0


	IF (RNORM.NE.0.0D0) THEN
	  RNORM = 1/RNORM
	  V2(1) = V1(1)*RNORM
	  V2(2) = V1(2)*RNORM
	  V2(3) = V1(3)*RNORM
	END IF
     	RETURN
	END
C
C-----------------------------------------------------------
C
C       generate a vector  vres = p2 - p1    ( P1 -> P2 )
C
C
	SUBROUTINE VECTOR (P1,P2,VRES)
	IMPLICIT 	REAL*8 		(A-H,O-Z)
	DIMENSION 	P1(3),P2(3),VRES(3)
	DO 100 I=1,3
	    VRES(I)  =   P2(I) - P1(I)

C
C If the numbers are *very* small, zero them out. This is not done
C for VMS, since it seems to work out fine. Why mess with something
C that already works.
C
	    IF (ABS(VRES(I)).LT.1.0E-31) VRES(I) = 0.0D0

100	CONTINUE
	END
C
C-----------------------------------------------------------
C
C	generate a versor
C
C
	SUBROUTINE VERSOR (P1,P2,VRES)
	IMPLICIT 	REAL*8 		(A-H,O-Z)
	DIMENSION 	P1(3),P2(3),VRES(3)
C
C **** CHECK FOR RNORM.EQ.0 SOMEBODY ******
C
	RNORM =    .0D0
	DO 100 I=1,3
100	  RNORM  =   RNORM + ( P1(I) - P2(I) )*( P1(I) - P2(I) )
	RNORM  =   SQRT(RNORM)
	DO 200 I=1,3
200	  VRES(I) =   (P2(I)-P1(I))/RNORM
	END
C
C--------------------------------------------------------------------
C
C 	project v1 onto v2
C
C
	SUBROUTINE PROJ (V1,V2,VRES)
	IMPLICIT 	REAL*8 		(A-H,O-Z)
        DIMENSION 	V1(3),V2(3),VRES(3)
	RNORM = V2(1)**2 + V2(2)**2 + V2(3)**2


C
C If the numbers are *very* small, zero them out. This is not done
C for VMS, since it seems to work out fine. Why mess with something
C that already works.
C
	IF (ABS(RNORM).LT.1.0E-31) RNORM = 0.0D0


	IF (RNORM.EQ.0.0D0) THEN
	  VRES(1) = 0.0D0
	  VRES(2) = 0.0D0
	  VRES(3) = 0.0D0
	ELSE
 	  SCALAR = V1(1)*V2(1) + V1(2)*V2(2) + V1(3)*V2(3)
	  PROD = SCALAR/RNORM
	  VRES(1) = V2(1)*PROD
	  VRES(2) = V2(2)*PROD
	  VRES(3) = V2(3)*PROD
	END IF
	END
C
C------------------------------------------------------------------
C
C       Generates the sum of two vectors
C
C
	SUBROUTINE SUM (P1,P2,RES)
	IMPLICIT 	REAL*8 		(A-E,G-H,O-Z)
	DIMENSION 	P1(3),P2(3),RES(3)

		RES(1) = P1(1) + P2(1)
		RES(2) = P1(2) + P2(2)
		RES(3) = P1(3) + P2(3)
	END
C
C------------------------------------------------------------------
C
C This subroutine returns the value of the arctangent between 0-2*PI
C
C
     	SUBROUTINE		ATAN_2	(SINE,COSINE,ANGLE)
     	IMPLICIT	REAL*8		(A-E,G-H,O-Z)
     	IMPLICIT	INTEGER*4	(F)
     	DATA	PI	/ 3.1415 92653 58979 32384 62643 D0 /
     	IF (COSINE.EQ.0.0D0.AND.SINE.EQ.0.0D0)	THEN
     		ANGLE	= 0.0D0
     		RETURN
     	ELSE
     	END IF
C
C Check if cosine is 0
C
     	IF (COSINE.NE.0) THEN
     	  ARG	=    SINE/COSINE
     	  ANGLE   =    ATAN (ABS(ARG))
C
C First quadrant: sine > 0, cosine > 0
C
     	 IF (SINE.GT.0.0D0.AND.COSINE.GT.0.0D0) THEN
     	   ANGLE	=   ANGLE
C
C Second quadrant: sine > 0, cosine < 0
C
         ELSE  IF (SINE.GT.0.0D0.AND.COSINE.LT.0.0D0) THEN
     	   ANGLE	=   PI - ANGLE
C
C Third quadrant: sine < 0, cosine < 0
C
     	 ELSE  IF (SINE.LT.0.0D0.AND.COSINE.LT.0.0D0) THEN
     	   ANGLE	=   ANGLE + PI
C
C Fourth quadrant: sine < 0, cosine > 0
C
     	 ELSE  IF (SINE.LT.0.0D0.AND.COSINE.GT.0.0D0) THEN
     	   ANGLE	=   2*PI - ANGLE
C
C Divide by zero cases
C
     	 ELSE IF (SINE.EQ.0.0D0.AND.COSINE.GT.0.0D0) THEN
     	   ANGLE	=   0.0D0
     	 ELSE IF (SINE.EQ.0.0D0.AND.COSINE.LT.0.0D0) THEN
     	   ANGLE	=   PI
     	 ELSE
     	 END IF
     	ELSE IF (SINE.GT.0.0D0) THEN
     	  ANGLE	=   PI/2
     	ELSE IF (SINE.LT.0.0D0) THEN
     	  ANGLE   =   1.5D0*PI
     	END IF
     	END
C
C-------------------------------------------------------------
C
C	Vector from a line to a point; line is specified by
C       H0 (starting), VH (direction); point is P0, output is DIS
C
C
     	SUBROUTINE	VDIST	(P0, H0, VH, DIS)
     	IMPLICIT	REAL*8	(A-H, O-Z)
     	DIMENSION	P0(3), H0(3), VH(3), DIS(3)
     	DIMENSION	VT(3), T_VH(3)

     	CALL	VECTOR	(H0, P0, VT)
     	CALL	PROJ	(VT, VH, T_VH)
     	CALL	VECTOR	(T_VH, VT, DIS)

     	RETURN
     	END
