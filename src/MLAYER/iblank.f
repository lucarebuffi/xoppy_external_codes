CL# 1 "iblank.F"
C +++
C
C Source: src/lib/iblank.F
C
C ----------------------------------------------
C                SHADOW
C      Center for X-ray Lithography
C     University of Wisconsin-Madison
C  3731 Schneider Dr., Stoughton, WI, 53589
C ----------------------------------------------
C
C Log: iblank.F
C Revision 1.6  1991/07/06  19:56:48  khan
C Grenoble and after. Minor changes
C
C Revision 1.5  91/04/05  13:54:21  cwelnak
C changed quotes on #include
C
C Revision 1.4  91/03/22  12:19:48  cwelnak
C SUN version -- INC to #inc
C
C Revision 1.3  90/11/13  14:04:52  khan
C Cleanup and SAVE statements
C
C Revision 1.2  90/07/20  22:05:08  khan
C put #if 1 ... to make it work also on vms
C
C Revision 1.1  90/07/18  17:35:06  khan
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
CL# 36 "iblank.F"




C +++
C 	integer 	function 	iblank
C
C	purpose		Returns the last non-white spot in the string.
C
C	input		a fortran character string.
C	output		Index of the last non-white char in the string.
C			If there are no empty spaces in the string, then
C			the lenght is simply returned.
C	hacker		Mumit Khan
C ---
	integer function iblank (str)
	implicit 	none
	character*(*) 	str
	integer 	index, ilen
c
c if the last character in the declared string isn't a white space, simply
c return the length.
c
	index = 1
	ilen = len (str)
	if (str(ilen:ilen).NE.' ') then
	    index = ilen + 1
	    goto 20
	endif
c
 10	continue
	if (str(index:index).NE.' ') then
	    index = index + 1
	    goto 10
	endif
 20	continue
	iblank = index - 1
	return
	end
