CL# 1 "iyes.F"
C +++
C
C Source: src/lib/iyes.F
C
C ----------------------------------------------
C                SHADOW
C      Center for X-ray Lithography
C     University of Wisconsin-Madison
C  3731 Schneider Dr., Stoughton, WI, 53589
C ----------------------------------------------
C
C Log: iyes.F
C Revision 1.6  1991/07/06  19:56:48  khan
C Grenoble and after. Minor changes
C
C Revision 1.5  91/04/05  13:54:25  cwelnak
C changed quotes on #include
C
C Revision 1.4  91/03/22  12:26:43  cwelnak
C SUN version -- INC to #inc
C
C Revision 1.3  90/11/13  14:04:55  khan
C Cleanup and SAVE statements
C
C Revision 1.2  90/07/20  22:05:11  khan
C put #if 1 ... to make it work also on vms
C
C Revision 1.1  90/07/10  14:56:23  khan
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
CL# 36 "iyes.F"




C+++
C	FUNCTION	IYES	(prompt)
C
C	PURPOSE		Decodes an answer line to a prompt.
C
C			If IN = Y or y or 1, then IYES = 1
C			else		 	  IYES = 0
C---
     	FUNCTION	IYES	(PROMPT)
     	CHARACTER *(*)	PROMPT
     	CHARACTER *1	IN
     	LOGICAL		TEST
     	IYES	= 1
10     	WRITE (6,1000)	PROMPT
1000	FORMAT	(1X,A,2x,$)
     	READ	(5,1010,ERR=20) IN
1010	FORMAT (1A1)
     	TEST = (IN(1:1).EQ.'Y').OR.(IN(1:1).EQ.'y').OR.(IN(1:1).EQ.'1')
     	IYES	= 1
     	IF (TEST) RETURN
     	IYES	= 0
     	RETURN
20     	WRITE(6,*)'What ?'
     	GO TO 10
     	END
