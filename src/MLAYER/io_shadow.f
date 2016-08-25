
CL# 1 "irint.F"
C +++
C
C Source: src/lib/irint.F
C
C ----------------------------------------------
C                SHADOW
C      Center for X-ray Lithography
C     University of Wisconsin-Madison
C  3731 Schneider Dr., Stoughton, WI, 53589
C ----------------------------------------------
C
C Log: irint.F
C Revision 1.7  1991/07/06  19:56:48  khan
C Grenoble and after. Minor changes
C
C Revision 1.6  91/04/05  13:54:24  cwelnak
C changed quotes on #include
C
C Revision 1.5  91/03/22  12:22:17  cwelnak
C SUN version -- INC to #inc
C
C Revision 1.4  90/11/13  14:04:53  khan
C Cleanup and SAVE statements
C
C Revision 1.3  90/07/20  22:05:09  khan
C put #if 1 ... to make it work also on vms
C
C Revision 1.2  90/07/15  15:31:02  khan
C All public include files (common.blk, etc) are now in ./../../include dir.
C
C Revision 1.1  90/07/10  14:56:20  khan
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
CL# 39 "irint.F"




C+++
C	SUBROUTINE	RINT
C
C	Purpose		Error-proof read-in routine for an I*4 number
C
C
C---
     	INTEGER 	FUNCTION	IRINT	( PROMPT )
     	CHARACTER*(*)	PROMPT
C
C The error iteration limit is 10 right now.
C
	ICOUNT	= 0
10	WRITE 	(6, '(1X,A,2X,$)')	PROMPT
     	READ	(5, '(I10.0)', IOSTAT=IRET)	IVAL
     	IF (IRET.NE.0) THEN
     	  WRITE (6,*)
     $	'What ? [ Program expects integer number input ]'
	  ICOUNT = ICOUNT + 1
	  IF (ICOUNT.GT.10)
     $	CALL LEAVE ('IRINT : ','Exceed error iteration limit.',IRET)
     	  GO TO 10
     	ELSE
     	  IRINT = IVAL
     	  RETURN
     	END IF
     	END




CL# 1 "rstring.F"
C +++
C
C Source: src/lib/rstring.F
C
C ----------------------------------------------
C                SHADOW
C      Center for X-ray Lithography
C     University of Wisconsin-Madison
C  3731 Schneider Dr., Stoughton, WI, 53589
C ----------------------------------------------
C
C Log: rstring.F
C Revision 1.7  1991/07/06  19:56:48  khan
C Grenoble and after. Minor changes
C
C Revision 1.6  91/04/05  13:54:38  cwelnak
C changed quotes on #include
C
C Revision 1.5  91/03/25  10:56:20  cwelnak
C SUN version - INC to #inc
C
C Revision 1.4  90/11/13  14:05:06  khan
C Cleanup and SAVE statements
C
C Revision 1.3  90/07/20  22:05:24  khan
C put #if 1 ... to make it work also on vms
C
C Revision 1.2  90/07/15  15:31:10  khan
C All public include files (common.blk, etc) are now in ./../../include dir.
C
C Revision 1.1  90/07/10  14:56:52  khan
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
CL# 39 "rstring.F"




C+++
C	FUNCTION	RSTRING
C
C	PURPOSE		Reads in an ASCII string
C
C---
     	FUNCTION	RSTRING ( ARG )
     	CHARACTER *80	RSTRING
     	CHARACTER *(*)	ARG
C
C The error iteration limit is 10 right now.
C
	ICOUNT = 0
1     	WRITE (6,1000)  ARG
     	READ  (5, 1010, END=20, ERR=10)	RSTRING
     	RETURN
10	WRITE (6,1000)  'I/O-%-ERR: What ?? Please try again.'
	ICOUNT	= ICOUNT + 1
	IF (ICOUNT.GT.10) THEN
	  CALL 	LEAVE ('RSTRING : ','Exceed error iteration limit.',0)
	ELSE
     	  GO TO 1
	END IF
20	RSTRING (1:2) = '^Z'
     	RETURN
1000	FORMAT (1X,A,$)
1010	FORMAT (A)
     	END



CL# 1 "rnumber.F"
C +++
C
C Source: src/lib/rnumber.F
C
C ----------------------------------------------
C                SHADOW
C      Center for X-ray Lithography
C     University of Wisconsin-Madison
C  3731 Schneider Dr., Stoughton, WI, 53589
C ----------------------------------------------
C
C Log: rnumber.F
C Revision 1.7  1991/07/06  19:56:48  khan
C Grenoble and after. Minor changes
C
C Revision 1.6  91/04/05  13:54:37  cwelnak
C changed quotes on #include
C
C Revision 1.5  91/03/25  10:53:28  cwelnak
C SUN version - INC to #inc
C
C Revision 1.4  90/11/13  14:05:05  khan
C Cleanup and SAVE statements
C
C Revision 1.3  90/07/20  22:05:23  khan
C put #if 1 ... to make it work also on vms
C
C Revision 1.2  90/07/15  15:31:09  khan
C All public include files (common.blk, etc) are now in ./../../include dir.
C
C Revision 1.1  90/07/10  14:56:50  khan
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
CL# 39 "rnumber.F"




C+++
C	SUBROUTINE	RNUMBER
C
C	Purpose		Error-proof read-in routine for a real*8 number
C
C
C---
     	DOUBLE PRECISION FUNCTION	RNUMBER	( PROMPT )
     	CHARACTER*(*)	PROMPT
     	REAL*8		VALUE
C
C The error iteration limit is 10 right now.
C
	ICOUNT	= 0
10	WRITE 	(6, '(1X,A,2X,$)')	PROMPT
     	READ	(5, '(F21.0)', IOSTAT=IRET)	VALUE
     	IF (IRET.NE.0) THEN
     	  WRITE (6,*)	'What ? [ Program expects real number input ]'
	  ICOUNT	= ICOUNT + 1
	  IF (ICOUNT.GT.10)
     $CALL LEAVE ('RNUMBER : ','Exceed error iteration limit.',IRET)
     	  GO TO 10
     	ELSE
     	  RNUMBER = VALUE
     	  RETURN
     	END IF
     	END


CL# 1 "leave.F"
C +++
C
C Source: src/lib/leave.F
C
C ----------------------------------------------
C                SHADOW
C      Center for X-ray Lithography
C     University of Wisconsin-Madison
C  3731 Schneider Dr., Stoughton, WI, 53589
C ----------------------------------------------
C
C Log: leave.F
C Revision 1.9  1992/01/23  16:45:17  cwelnak
C 6000 changes
C
C Revision 1.8  1991/07/06  19:56:48  khan
C Grenoble and after. Minor changes
C
C Revision 1.7  1991/04/22  00:39:40  khan
C EXIT() doesn't seem to work on Ultrix. Replace with STOP. YUK!
C
C Revision 1.6  1991/04/05  13:54:26  cwelnak
C changed quotes on #include
C
C Revision 1.5  91/03/22  12:27:46  cwelnak
C SUN version -- INC to #inc, PARAMETER stmt
C
C
C Revision 1.4  90/11/13  14:04:55  khan
C Cleanup and SAVE statements
C
C Revision 1.3  90/07/22  15:13:21  khan
C Replace call to STOP with a call to EXIT (seems more portable and is
C able to return a code to OS).
C
C Revision 1.2  90/07/20  22:05:12  khan
C put #if 1 ... to make it work also on vms
C
C Revision 1.1  90/07/10  14:56:27  khan
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
CL# 47 "leave.F"




C+++
C	SUBROUTINE	LEAVE
C
C	PURPOSE		To print out a message, then signal the VAX
C			condition handler to abort the command procedure.
C---
	SUBROUTINE	LEAVE	(TEXT1,TEXT2,IFLAG)
	CHARACTER*(*)	TEXT1,TEXT2
	INTEGER		SS_ABORT
	PARAMETER	(SS_ABORT = 44)
	CALL	MSSG	(TEXT1,TEXT2,IFLAG)
C	STOP		SS_ABORT
C *update* when 6000 compiler improves to include calls to EXIT.



	CALL	EXIT	(SS_ABORT)

	END


CL# 1 "message.F"
C +++
C
C Source: src/lib/message.F
C
C ----------------------------------------------
C                SHADOW
C      Center for X-ray Lithography
C     University of Wisconsin-Madison
C  3731 Schneider Dr., Stoughton, WI, 53589
C ----------------------------------------------
C
C Log: message.F
C Revision 1.7  1991/07/06  19:56:48  khan
C Grenoble and after. Minor changes
C
C Revision 1.6  91/04/05  13:54:28  cwelnak
C changed quotes on #include
C
C Revision 1.5  91/03/25  10:29:17  cwelnak
C SUN version - INC to #inc
C
C Revision 1.4  90/11/13  14:04:57  khan
C Cleanup and SAVE statements
C
C Revision 1.3  90/07/20  22:05:14  khan
C put #if 1 ... to make it work also on vms
C
C Revision 1.2  90/07/15  15:31:03  khan
C All public include files (common.blk, etc) are now in ./../../include dir.
C
C Revision 1.1  90/07/10  14:56:30  khan
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
CL# 39 "message.F"




C+++
C	SUBROUTINE	MSSG	(TEXT1,TEXT2,IFLAG)
C
C	PURPOSE		Outputs an error message to the screen and
C			to a file
C
C	ARGUEMNTS	TEXT1 : Tipically, the calling routine name
C			TEXT2 : A file name, or a message
C			TEXT3 : A flag or an IOSTAT value.
C
C---
     	SUBROUTINE	MSSG (T1, T2, IFLAG)
     	CHARACTER *(*)	T1, T2
     	WRITE(6,*)'SHADOW-E-Error: '
     	WRITE(6,*)'Module     : ', T1
     	WRITE(6,*)'Message    : ', T2
     	WRITE(6,*)'Error flag : ', IFLAG
     	WRITE (33, 100) T1
     	WRITE (33, 100) T2
     	WRITE (33, 110) IFLAG
100	FORMAT (1X,'>',A)
110	FORMAT (1X,I5)
     	RETURN
     	END

CL# 1 "messages.F"
C +++
C
C Source: src/lib/messages.F
C
C ----------------------------------------------
C                SHADOW
C      Center for X-ray Lithography
C     University of Wisconsin-Madison
C  3731 Schneider Dr., Stoughton, WI, 53589
C ----------------------------------------------
C
C Log: messages.F
C Revision 1.6  1991/07/06  19:56:48  khan
C Grenoble and after. Minor changes
C
C Revision 1.5  91/04/05  13:54:29  cwelnak
C changed quotes on #include
C
C Revision 1.4  91/03/25  10:30:44  cwelnak
C SUN version -- INC to #inc
C
C Revision 1.3  90/11/13  14:04:58  khan
C Cleanup and SAVE statements
C
C Revision 1.2  90/07/20  22:05:15  khan
C put #if 1 ... to make it work also on vms
C
C Revision 1.1  90/07/10  14:56:31  khan
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
CL# 36 "messages.F"




C+++
C	SUBROUTINE	TMSSG	(TEXT1,TEXT2,IFLAG)
C
C	PURPOSE		Outputs a simple message to the screen and
C			to a file. Typically used for error messages.
C
C	ARGUEMNTS	TEXT1 : Tipically, the calling routine name
C			TEXT2 : A file name, or a message
C			IFLAG : A flag or an IOSTAT value. Written only
C				if negative. A Fortran error decoder will
C				be included later.
C
C	COMMON		IO_UNITS	I_tt	TT: unit number (6)
C					I_log   Log file unit number
C						(the log file will be
C						opened and closed
C						elsewhere)
C---
     	SUBROUTINE	TMSSG (T1, T2, IFLAG)
C
     	COMMON	/IO_UNITS	/I_tt,I_log
C
     	CHARACTER *(*)	t1, t2
C
C Write to terminal first
C
     	WRITE (I_tt, 100) t1, t2
     	IF (iflag.NE.0) WRITE (I_tt, 110) t1, iflag
C
C ... then to log file.
C
     	WRITE (I_log, 100) t1, t2
     	IF (iflag.NE.0) WRITE (I_log, 110) t1, iflag
C
     	RETURN
C
100	FORMAT (1X,A,'> ',A)
110	FORMAT (1X,A,'> ',I5)
     	END
C
C+++
C	SUBROUTINE	LINEOUT		(Text)
C
C	PURPOSE		To output a formatted line to the terminal and
C			to a file.
C
C---
     	SUBROUTINE	LINEOUT		(Text)
C
     	CHARACTER *(*)	Text
C
     	COMMON	/IO_UNITS	/I_tt,I_log
C
     	WRITE	(I_tt,  *) Text
     	WRITE	(I_log, *) Text
C
     	RETURN
     	END
