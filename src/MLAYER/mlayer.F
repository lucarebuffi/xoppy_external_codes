CL# 1 "mlayer.F"
C
C Source: src/utils/optics/mlayer.F
C
C ----------------------------------------------
C                SHADOW
C      Center for X-ray Lithography
C     University of Wisconsin-Madison
C  3731 Schneider Dr., Stoughton, WI, 53589
C ----------------------------------------------
C
C Log: mlayer.F
C Revision 1.9  1995/03/15  17:05:29  cwelnak
C Changed to read from optical constants database.
C Changes by M.S. Rio
C
C Revision 1.8  1995/03/15  17:02:20  cwelnak
C next revision will have many changes, making sure this
C one is saved.
C
C Log: mlayer.F
C Revision 1.9  1995/03/15  17:05:29  cwelnak
C Changed to read from optical constants database.
C Changes by M.S. Rio
C
C Revision 1.7  1992/01/24  17:37:06  cwelnak
C declaration troubles.
C
C Revision 1.6  1992/01/15  16:26:51  cwelnak
C replaced with mlayer from 6000
C
C Revision 1.5  91/07/06  19:42:54  khan
C Grenoble Changes ...
C
C Revision 1.4  91/04/05  15:21:23  cwelnak
C changed quotes in #includes
C
C Revision 1.3  91/03/25  15:52:11  khan
C SUN port -- INCLUDE -> #include
C
C Revision 1.2  91/03/25  15:48:48  khan
C Getting ready for SUN port ...
C
C Revision 1.1  90/10/10  13:32:53  khan
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
CL# 51 "mlayer.F"





c
c  program ULAYER4.FOR
c
     	IMPLICIT	REAL*8	(A-H,O-Z)
     	IMPLICIT	INTEGER*4	(I-M)
	CHARACTER*80	OUTFIL
        dimension do(200),de(200)
        real*8 fi(1000),ansx(1000),anpx(1000),anux(1000)
     	real*8 qvar(1000),faS(1000),FAP(1000)
        character*80 name1, rstring
c       logical*1 elo(5),ele(5),els(5)
c        character*1 elo(5),ele(5),els(5)
        character*2 elo,ele,els
c
        common /aaa/ do,de
        common /bbb/ delo,beto,dele,bete,dels,bets,x1,x2,step,n,nstep,
     1xlam,rhoo,rhoe,rhos,elo,ele,els
c
c     data delo,beto,dele,bete,dels,bets,x1,x2,step,xlam /0.0D+00,0.0D+00
c    1,0.0D+00,0.0D+00,0.0D+00,0.0D+00,0.0D+00,0.0D+00,0.0D+00,0.0D+00/
c       data n,nstep /0,0/
c       data elo /' ',' ',' ',' ',' '/
c       data ele /' ',' ',' ',' ',' '/
c       data els /' ',' ',' ',' ',' '/

     	DATA	PI     	/  3.1415 92653 58979 32384 62643 D0 /
     	DATA	PIHALF 	/  1.5707 96326 79489 66192 31322 D0 /
     	DATA	TWOPI 	/  6.2831 85307 17958 64679 25287 D0 /
     	DATA	TODEG 	/ 57.2957 79513 08232 08767 98155 D0 /
     	DATA	TORAD	/  0.0174 53292 51994 32957 69237 D0 /
	DATA	TOCM	/  1.239 852	D-4		     /
	DATA	TOANGS 	/  1.239 852    D+4		     /
c
c  program to compute x-ray/u.v. reflection efficiency of
c  multilayer structures at various bragg angles
c
c  this version allows specification of the individual
c  layer thicknesses.
c
c  multilayer and scan parameters can be input from the terminal
c  or read in from an existing file (and subsequently changed).
c  These parameters can optionally be output into a file for posterity.
c
c  input parameters:
c  n = no. of layer pairs (n=-1 means no more data)
c  xlam = wavelength (a)
c  elo = odd layer material
c  ele = even layer material
c  els = substrate material
c  do = thickness of odd layers (a)
c  1.0 - delo - i*beto = complex refractive index (odd)
c  de = thickness of even layers (a)
c  1.0 - dele - i*bete = complex refractive index (even)
c  phib = starting bragg angle
c  step = bragg angle step size
c  nstep = no. of steps
c
c
c
1       write(6,*) 'option 1, 2, 3, 4, or 5?  (0 = help, 6 = exit)'
        read(5,*) iopfl         !from tty
        if(iopfl.eq.6)stop
        if(iopfl.gt.5.or.iopfl.lt.0) goto 1
        if(iopfl.ne.0) goto 10
        write(6,530)
        goto 1
10      call readin(iopfl)
        if(n.lt.0) goto 81
c
c       write(6,*) 'name of output file:  '
        name1 = rstring('name of output file:  ')
c       read(5,550)l,name1
c       if(l.eq.0) goto 25

	call set_refraction(-1.0D0)

        open(unit=3,file=name1,status='UNKNOWN',form='FORMATTED')
	rewind (3)

25      SUM = 0.0
	if(iopfl.eq.2.or.iopfl.eq.4)goto 50
        IF (IOPFL.EQ.5) GO TO 70
30      phib = x2               !angular scan
	phimax = pihalf - phib		!in degree measured from the normal
	phimin = phimax - step*(nstep-1)
        xlam = x1
csrio	write (*,*) 'calling set_refraction: xlam = ',xlam
	call set_refraction(xlam)
csrio	write (*,*) 'returning set_refraction'

        do 40 k = nstep,1,-1		!scan from phimin to phimax
        phi = phib + step * (k-1)
        call fresnel(phi,anp,ans,anu,phaseS,PHASEP,ESQZS,ESQZP)
        sum = sum + anu
        ansx(k) = ans
        anpx(k) = anp
        anux(k) = anu
        faS(k) = phaseS
	FAP(K) = PHASEP
        fi(k) = phi
40      write(3,501) phi,ans,anp,anu,phaseS,PHASEP,ESQZS,ESQZP
        sum = sum*step*TODEG
c************************************************************
c       call plots(phib,phi) PLOT REFLECTIVITY VS ANGLE HERE
c************************************************************
C        write(3,502)sum
        close(unit=3)
C
C Graft by F.Cerrina, to prepare an output file for SHADOW
C
C	TYPE	*,'File name for reflectivities and phase changes ?'
C	READ	(5,1010)	OUTFIL
C1010	FORMAT	(A)
C     	QWAVE	= TWOPI/XLAM*1.0E8		!cm-1
C     	OPEN (33, FILE=OUTFIL, STATUS='NEW', FORM='UNFORMATTED')
C     	  WRITE (33) phimin,phimax,step,qwave	!degree and cm-1
C     	  WRITE (33) NSTEP
C     	  WRITE (33) (ANSX(I),I=1,NSTEP)
C     	  WRITE (33) (ANPX(I),I=1,NSTEP)
C     	  WRITE (33) (FAS(i),I=1,NSTEP)
C	  WRITE (33) (FAP(I),I=1,NSTEP)
C     	CLOSE (33)
        goto 1
50      phi = x2                !wavelength scan
        xlamb = x1

        do 60 k = 1,nstep
        xlam = xlamb + step * (k-1)
	call set_refraction(xlam)
        call fresnel(phi,anp,ans,anu,phaseS,PHASEP,ESQZS,ESQZP)
        sum = sum + anu
        ansx(k) = ans
        anpx(k) = anp
        anux(k) = anu
        faS(k) = phaseS
	FAP(K) = PHASEP
        fi(k) = phi
60      write(3,501) XLAM,ans,anp,anu,phaseS,PHASEP,ESQZS,ESQZP
        sum = sum * step
c***********************************************************************
c       call plots(xlamb,xlam)   PLOT REFLECTIVITY VS WAVELENGTH HERE
c***********************************************************************
C        write(3,502)sum
        close(unit=1)
        goto 1
70      phi = x2                !wavevector scan
        xlamb = x1
        dblo = delo
        dble = dele
        bbto = beto
        bbte = bete
C SHADOW assume scanning from QMIN(xlambmax) to QMAX(xlamb) in step of QSTEP
     	qmax = twopi/xlamb
	q2 = twopi/(xlamb+step)
	qstep = q2 - qmax

        do 80 k = 1,nstep
     	qvar(k) = qmax + qstep *(k-1)      ! changed by MSR
        xlam = twopi/qvar(k)
	call set_refraction(xlam)
        call fresnel(phi,anp,ans,anu,phaseS,PHASEP,ESQZS,ESQZP)
        sum = sum + anu
        ansx(k) = ans
        anpx(k) = anp
        anux(k) = anu
        faS(k) = phaseS
	FAP(K)	= PHASEP
        fi(k) = phi
80      write(3,501) XLAM,ans,anp,anu,phaseS,PHASEP,ESQZS,ESQZP
        sum = sum * step
c***********************************************************************
c       call plots(xlamb,xlam)   PLOT REFLECTIVITY VS WAVELENGTH HERE
c***********************************************************************
C        write(3,502)sum
        close(unit=1)
     	GO TO 1
500     format('0','integrated reflectivity = ',1pe12.4)
501     format(1x,8f12.8)
502     format(1x,e12.8)
530     format(' 1 = fixed wavelength, scan angle, peroidic layers',
     1/,' 2 = fixed angle, scan wavelength, periodic layers',
     2/,' 3 = fixed wavelength, scan angle, individual layers defined',
     3/,' 4 = fixed angle, scan wavelength, individual layers defined',
     4/,' 5 = fixed angle, sacn wavevector, individual layers defined')
550     format(13a1)
81      call exit (0)
        end
c
c
c
c
c
        subroutine FRESNEL(phi,anp,ans,anu,phaseS,PHASEP,ESQZS,ESQZP)
c
     	IMPLICIT	REAL*8	(A-H,O-Z)
     	IMPLICIT	INTEGER*4	(I-M)

     	DATA	PI     	/  3.1415 92653 58979 32384 62643 D0 /
     	DATA	TODEG 	/ 57.2957 79513 08232 08767 98155 D0 /
        complex*16 ci,fo,fe,fv,ffe,ffv,ffvp,ffo,ffep,ffop,re2,
     $		   ro2,ao,ae,r,rp,fs,ffs,ffsp,rs2
        dimension do(200),de(200)
        character*2 ele,elo,els
        common /aaa/ do,de
        common /bbb/ delo,beto,dele,bete,dels,bets,x1,x2,step,n,nstep,
     1xlam,rhoo,rhoe,rhos,elo,ele,els
c
      ci=(0.0,1.0)
      ro2=(1-delo-ci*beto)**2
      re2=(1-dele-ci*bete)**2
      rs2=(1-dels-ci*bets)**2
      phr=phi/TODEG
C      refo=(sin(phr))**2-2.0*delo
C      xmfo=-2.0*beto
C      fo=cmplx(refo,xmfo)
      fo = ro2 - cos(phr)**2
C      refe=(sin(phr))**2-2.0*dele
C      xmfe=-2.0*bete
C      fe=cmplx(refe,xmfe)
      fe = re2 - cos(phr)**2
      refv=(sin(phr)**2)
      xmfv=0.0
      fv = Dcmplx(refv,xmfv)
C      refs=(sin(phr))**2-2.0*dels
C      xmfs=-2.0*bets
C      fs=cmplx(refs,xmfs)
      fs = rs2 - cos(phr)**2
C
      fo=sqrt(fo)
      fe=sqrt(fe)
      fv=sqrt(fv)
      fs=sqrt(fs)
      ffe=(fe-fo)/(fe+fo)
      ffo=-ffe
      ffv=(fv-fo)/(fv+fo)
      ffs=(fe-fs)/(fe+fs)
      ffep=(fe/re2-fo/ro2)/(fe/re2+fo/ro2)
      ffop=-(ffep)
      ffvp=(fv-fo/ro2)/(fv+fo/ro2)
      ffsp=(fe/re2-fs/rs2)/(fe/re2+fs/rs2)
      r=(0.0,0.0)
      rp=(0.0,0.0)
      do 1 j=1,n
      ao=-ci*(pi*fo*do(j)/xlam)
      ae=-ci*(pi*fe*de(j)/xlam)
      ao=exp(ao)
      ae=exp(ae)
      if(j.eq.1)go to 6
      r=(ae**4)*(r+ffe)/(r*ffe+1.0)
      rp=(ae**4)*(rp+ffep)/(rp*ffep+1.0)
      go to 7
    6 r=(ae**4)*(r+ffs)/(r*ffs+1.0)
      rp=(ae**4)*(rp+ffsp)/(rp*ffsp+1.0)
    7 r=(ao**4)*(r+ffo)/(r*ffo+1.0)
      rp=(ao**4)*(rp+ffop)/(rp*ffop+1.0)
    1 continue
      r=(r+ffv)/(r*ffv+1.0)
      ESQZS = ABS(1.0 + R)
      ESQZS = ESQZS**2
      pp = imag(r)
#if lahey
      qq = real(r,8)
#else
      qq = Dreal(r)
#endif
C      phaseS = atan_2(pp,qq)/pi		! S phase change in units of pi
	call atan_2(pp, qq, phaseS)
	phaseS = phaseS/pi
      rp=(rp+ffvp)/(rp*ffvp+1.0)
C Check the following equations :
      ESQZP = ABS(1.0 - RP)
      ESQZP = ESQZP**2
C
      PP = IMAG(RP)
#if lahey
      QQ = REAL (RP,8)
#else
      QQ = dREAL (RP)
#endif
C      PHASEP = ATAN_2(PP,QQ)/PI		! P phase change in units of pi
	call atan_2 (pp, qq, phasep)
	phaseP = phaseP /Pi
      anp=abs(rp)
      anp=anp**2
      ans=abs(r)
      ans=ans**2
      anu=0.5*(ans+anp)
c
        return
        end
c
c
        subroutine READIN(iopfl)
c
c
     	IMPLICIT	REAL*8	(A-H,O-Z)
     	IMPLICIT	INTEGER*4	(I-M)
        dimension do(200),de(200)
c       logical*1 ans,ele(5),elo(5),els(5),tmp(5)
c        character*1 ele(5),elo(5),els(5)
        character*2 ele,elo,els
        character*80 	filnam, rstring
	character*80	tmpstr
c
        common /aaa/ do,de
        common /bbb/ delo,beto,dele,bete,dels,bets,x1,x2,step,n,nstep,
     1xlam,rhoo,rhoe,rhos,elo,ele,els
c
c       write(6,*) 'use existing input file?  (cr = yes)'
C        read(5,550) l,ans
c       read(5,550) ans
c       if(l.ne.0.and(.ans.ne.'y'.and.ans.ne.'Y')) goto 21

	IFILE = IYES('Use existing input file?')
	IF (IFILE.EQ.0) GOTO 41

c
c 5       write(6,*) 'input file name:'

5       filnam = rstring('input file name:')

c       read(5,552) l,filnam
c       if(l.eq.0)stop
        open(unit=1,file=filnam,form='FORMATTED',status='OLD',err=8)
        read(1,503)els
        read(1,*)rhos
        read(1,503)elo
        read(1,*)rhoo
        read(1,503)ele
        read(1,*)rhoe
	read(1,*)x1,x2
	read(1,*)step,nstep
	read(1,*)n
        goto 9
8       write(6,*) '****** error reading input file'
        goto 5
c
9       if(iopfl.gt.2) goto 10          !branch for non-periodic layers
        read(1,*) xdo,xde               !read for periodic layers
        do 15 i=1,n
        do(i)=xdo
15      de(i)=xde
        goto 20
10      read(1,*) (do(i),de(i), i=1,n)  !read for non-periodic layers
20      close(unit=1)
c
c       write(6,*) 'change multilayer?  (cr = no)'
c       read(5,550) l,ans
c       if(l.eq.0.or.(ans.ne.'y'.and.ans.ne.'Y')) goto 40
c       goto 29

	ICHANGE = IYES('Do you want to change the multilayers?')
	IF (ICHANGE.EQ.0) GOTO 49
C
C  else chanege the m-layer file.
C
	WRITE(6,1550) ELS
	IANS = IYES('Change it?')
	IF (IANS.EQ.1) THEN
		WRITE(6,*) 'New substrate:'
		READ(5,554) ELS
	ENDIF
	WRITE(6,1553) RHOS
	IANS = IYES('Change it?')
	IF (IANS.EQ.1) THEN
		WRITE(6,*) 'New rho:'
		READ(5,*) RHOS
	ENDIF
	WRITE(6,1551) ELO
	IANS = IYES('Change odd layer?')
	IF (IANS.EQ.1) THEN
		WRITE(6,*) 'New material:'
		READ(5,554) ELO
	ENDIF
	WRITE(6,1553) DELO,BETO,RHOO
	IANS = IYES('change rho?')
	IF (IANS.EQ.1) THEN
		WRITE(6,*) 'New rho:'
		READ(5,*) RHOO
	ENDIF
	WRITE(6,1552) ELE
	IANS = IYES('Change even layer?')
	IF (IANS.EQ.1) THEN
		WRITE(6,*) 'new material:'
		READ(5,554) ELE
	ENDIF
	WRITE(6,1553) RHOE
	IANS = IYES('Change rho?')
	IF (IANS.EQ.1) THEN
		WRITE(6,*) 'New rho:'
		READ(5,*)  RHOE
	ENDIF
	WRITE(6,1513) X1,X2
	IANS = IYES('Change lambda, theta ?')
	IF (IANS.EQ.1) THEN
		WRITE(6,*) 'New lambda, theta:'
		READ(5,*) X1,X2
	ENDIF
	WRITE(6,1514) STEP, NSTEP
	IANS = IYES('Change step size, number of steps?')
	IF (IANS.EQ.1) THEN
		WRITE(6,*) 'New step size, number of steps:'
		READ(5,*) STEP, NSTEP
	ENDIF
	WRITE(6,1501) N
	IANS = IYES('Change number of layer pairs?')
	IF (IANS.EQ.1) THEN
		WRITE(6,*) 'Number of layer pairs:'
		READ(5,*) N
	ENDIF
	IF (IOPFL.GT.2) GOTO 30
	WRITE(6,1511) XDO, XDE
	IANS = IYES('Change thickness of layer pairs?')
	IF (IANS.EQ.1) THEN
		WRITE(6,*) 'New thicknesses:'
		READ(5,*) XDO, XDE
	ENDIF
	DO 28 I=1,N
		DO(I) = XDO
		DE(I) = XDE
28	CONTINUE
	GOTO 49

30	DO 35 I=1,N
	WRITE(6,1515) I,DO(I),DE(I)
c
c	The following was a long standing bug in SHADOW! -- MK
c
c	IANS = IYES('Change thickness of layer ',I)
c
	WRITE(TMPSTR,*) 'Change thickness of layer ', I
	IANS = IYES(TMPSTR)
	IF (IANS.EQ.1) THEN
		WRITE(6,*) 'New thicknesses:'
		READ(5,*) DO(I), DE(I)
	ENDIF
35	CONTINUE
	GOTO 49

c
41      write(6,*) 'input multilayer parameters'
Csrio	WRITE(6,*) 'Enter substrate material'
Csrio	READ(5,554) ELS
Csrio	WRITE(6,*) 'Enter density [g/cm^3] for substrate'
Csrio	READ(5,*) RHOS
Csrio	WRITE(6,*) 'Enter odd layer material'
Csrio	READ(5,554) ELO
Csrio	WRITE(6,*) 'Enter density [g/cm^3] for odd material'
Csrio	READ(5,*) RHOO
Csrio	WRITE(6,*) 'Enter even later material'
Csrio	READ(5,554) ELE
Csrio	WRITE(6,*) 'Enter density [g/cm^3] for even material'
Csrio	READ(5,*) RHOE
	WRITE(6,*) 'Enter lambda, theta'
	READ(5,*) X1, X2
	WRITE(6,*) 'Enter step size, number of steps'
	READ(5,*) STEP, NSTEP
	WRITE(6,*) 'Enter number of layer pairs'
	READ(5,*) N
	IF (IOPFL.GT.2) GOTO 43

	WRITE(6,*) 'Enter thickness for odd layer and for even layer'
	READ(5,*) XDO, XDE
	DO 42 I=1,N
		DO(I) = XDO
		DE(I) = XDE
42	CONTINUE
	GOTO 49

43	DO 45 I=1,N
	WRITE(6,*) 'For layer ',I,' enter odd and even thicknesses'
	READ(5,*) DO(I), DE(I)
45	CONTINUE
	GOTO 49


c
c 40      write(6,*) 'save layer parms in which file? (cr = don''t save)'

49      filnam =
     $  rstring('save layer parms in which file? :cr = do not save:')

c       read(5,552) l,filnam
c       if(l.eq.0)goto 60

        open(unit=2,file=filnam,status='UNKNOWN',form='FORMATTED')
	rewind (2)

c        write(2,502)els
        write(2,503)els
        write(2,9512)rhos
c        write(2,502)elo
        write(2,503)elo
        write(2,9512)rhoo
c        write(2,502)ele
        write(2,503)ele
        write(2,9512)rhoe
        write(2,512)x1,x2
        write(2,514)step,nstep
        write(2,501)n
        write(2,512)do(1),de(1)
        if(iopfl.le.2) goto 51
        do 50 i=2,n
50      write(2,512)do(i),de(i)
c
51      close(unit=2)
c
501     format(1x,i4)
502     format(1x,5a1)
c503     format(5a1)
503     format(a2)
512     format(1x,2f10.5)
9512     format(1x,3f10.5)
514     format(1x,f10.5,i4)
550     format(a1)
552     format(13a1)
553     format(5a1)
c554     format(5a1)
554     format(a2)
c
1501    format(' number of layer pairs = ',i4)
1513    format(' lambda(begin) =',f10.5,'   theta(begin) =',f10.5)
c     1(, = no change)')
1514    format(' step size = ',f10.5,'   number of steps = ',i4)
c     1(, = no change)')
1515    format(' for layer pair ',i4,',',/,'   odd thickness = ',f10.5,
     1'   even thickness = ',f10.5)
1511    format(' for each layer pair,',/,'   odd thickness = ',f10.5,
     1'   even thickness = ',f10.5)
c1550    format(1x,'substrate material is ',5a1)
c1551    format(1x,'odd layer material is ',5a1)
c1552    format(1x,'even layer material is ',5a1)
1550    format(1x,'substrate material is ',a2)
1551    format(1x,'odd layer material is ',a2)
1552    format(1x,'even layer material is ',a2)
1553    format(1x,'delta = ',f10.5,',       beta = ',f10.5,
     1'       rho = ',f10.5)
c     1(, = no change)')
c
60      return
        end

	BLOCK DATA
	IMPLICIT REAL*8 (A-H,O-Z)
     	IMPLICIT	INTEGER*4	(I-M)
c        character*1 elo(5),ele(5),els(5)
        dimension do(200),de(200)
        character*2 elo,ele,els

	common /aaa/ do,de
	common /bbb/ delo,beto,dele,bete,dels,bets,x1,x2,step,n,nstep,
     $ xlam,rhoo,rhoe,rhos,elo,ele,els
c
        data delo,beto,dele,bete,dels,bets,x1,x2,step,xlam /0.0D+00,
     $	0.0D+00,0.0D+00,0.0D+00,0.0D+00,0.0D+00,0.0D+00,0.0D+00,0.0D+00,
     $  0.0D+00/
        data n,nstep /0,0/
c        data elo / ' ', ' ', ' ', ' ', ' ' /
c        data ele / ' ', ' ', ' ', ' ', ' ' /
c	data els / ' ', ' ', ' ', ' ', ' ' /
        data elo / '  ' /
        data ele / '  ' /
	data els / '  ' /

	END

	SUBROUTINE set_refraction(xlambda)
c
c
c added by M. Sanchez del Rio to automatically calculate the refractive
c index values (beta and delta) for a given wavelength
c February 1994
c
c input arguments:
c		xlambda: if negative (first call) calls the SHADIOW
c		Optical Library and stores the data. For sucessive
c		calls, xlambda is the photon wavelength in Angstroms.
c output arguments:
c		As a result, the refractive index values delta and beta
c		for the substrate, odd and even materials are updated.
c		The updated values are passed by the common blocks.
c
c
	implicit real*8 (a-h,o-z)
     	IMPLICIT	INTEGER*4	(I-M)
        character*2 elo,ele,els
        real*4  f1s(1420),f2s(1420)
        real*4  f1e(1420),f2e(1420),f1o(1420),f2o(1420)
        real*4	energy(1420),atwts,atwto,atwte,rmu,emf
C        integer*2 npts
	character*180  tmps
     	DATA	PI     	/  3.1415 92653 58979 32384 62643 D0 /
	DATA	TOCM	/  1.2398418757279403 D-4  /
	DATA	TOANGS 	/  1.239841875727940356810564   D+4  /
	common /bbb/ delo,beto,dele,bete,dels,bets,x1,x2,step,n,nstep,
     $ xlam,rhoo,rhoe,rhos,elo,ele,els
	save atwts,atwte,atwto,energy,f1s,f2s,f1e,f2e,f1o,f2o,npts
        RADIUS  =   2.8179403 E-13
        AVOG    =   6.0221418 E+23


c
c call library
c
	if (xlambda.lt.0.0) then
C	  write (*,*) '^^^^^^^^^^^^^ calling readlib ^^^^^^^^^^^^^^'
C	  write (*,*) '^^^^^^^^^^^^^ calling readlib S ^^^^^^^^^^^^^^'
CCsrio number,rmu,emf not used
C	  call readlib (els,number,atwts,rmu,emf,energy,f1s,f2s)
C	  call readlib (ele,number,atwte,rmu,emf,energy,f1e,f2e)
C	  call readlib (elo,number,atwto,rmu,emf,energy,f1o,f2o)

Cc
Cc	srio:  save to file
Cc
C	npts = 420	
C        open(unit=8,file='mlayer.f12',status='UNKNOWN',form='FORMATTED')
C	  write (8,'(A)') '  '//els//' '//ele//' '//elo
C	  write (8,'(1x,3f10.4)') atwts,atwte,atwto
C	  write (8,*) npts
C	  DO i=1,npts
C	    write (8,'(f12.4,6f14.8)') energy(i),f1s(i),f2s(i),f1e(i),
C     $                  f2e(i),f1o(i),f2o(i)
C	  ENDDO
C	close(8)
c
c	srio:  clean variables 
c
	DO i=1,1420
	  energy(i)=0.0
	  f1s(i) = 0.0
	  f2s(i) = 0.0
	  f1e(i) = 0.0
	  f2e(i) = 0.0
	  f1o(i) = 0.0
	  f2o(i) = 0.0
	  atwts = 0.0
	  atwte = 0.0
	  atwto = 0.0
	ENDDO
c
c	srio:  read from file mlayer.f12
c
	write (*,*) 'Getting optical ctes from file: mlayer.f12'
        open(unit=8,file='mlayer.f12',status='OLD',form='FORMATTED')
	  read (8,'(A)') tmps
C	  read (8,'(1x,3f10.4)') atwts,atwte,atwto
	  read (8,'(A)') tmps
	  read (8,*) atwts,atwte,atwto
	  read (8,'(A)') tmps
	  read (8,*) rhos,rhoe,rhoo
	  read (8,'(A)') tmps
	  read (8,*) npts
	  read (8,'(A)') tmps
	  DO i=1,1420
C	    read (8,'(f12.4,6f14.8)') energy(i),f1s(i),f2s(i),f1e(i),
	    read (8,*,END=1234) energy(i),f1s(i),f2s(i),f1e(i),
     $                  f2e(i),f1o(i),f2o(i)
	  ENDDO
1234    continue
	close(8)
	npts = i
	write (*,*) 'Number of energy points read ',npts
	write (*,*) 'Done reading optical ctes from file: mlayer.f12'
c	write (*,*) '>>> npts: ',npts
c	write (*,*) '>>> el[s,e,o]: ',els,ele,elo
c	write (*,*) '>>> atwt[s,e,o]: ',atwts,atwte,atwto
c	write (*,*) '>>> energy,f1s,f2s,f1e,f2e,f1o,f2o '
c	DO i=1,npts
c	  write (*,*) energy(i),f1s(i),f2s(i),f1e(i),
c     $                f2e(i),f1o(i),f2o(i)
c	ENDDO
c
	

	  return
	endif
c
c calculates delta anb beta for lambda(=x1)
c
	photon = TOANGS/xlambda
csrio	write (*,*) 'lambda,energy = ',xlambda,photon
	ELFACTOR        = LOG10(1.0E4/30.0)/300.0
	NENER   =   LOG10(PHOTON/30.0E0)/ELFACTOR + 1

C
C srio: recompute index (=nener2) for a generic energy array 
C (not necessarly in SHADOW's optical library). 
C
	ifound = 0
C        write (*,*) '<>lambda,energy,npts = ',xlambda,photon,npts
C        write (*,*) '<>npts = ',npts
	DO i=1,npts 
          IF ( (photon.GE.energy(i)).and.(photon.LT.energy(i+1)))
     $       ifound = 1
	  IF (ifound.eq.1) GO TO 1011
	ENDDO
	write(*,*) '****** Error interpolating energy: Table not valid'
	return
1011	continue
	nener2 = i
	nener = nener2



        WAVE    =   TOCM/PHOTON
csrio	write (*,*) 'nener= ',nener
csrio	  write (*,*) '>>>>>> photon: ',photon
csrio	  write (*,*) '>>>>>> nener: ',nener
csrio	  write (*,*) '>>>>>> nener2: ',nener2
csrio	  write (*,*) '>>>>>> energy(nener): ',energy(nener)
csrio	  write (*,*) '>>>>>> energy(nener+1): ',energy(nener+1)
c
c for the substrate
c
        F1_S      =   f1s(NENER) +
     $  (f1s(NENER+1) - f1s(NENER))*(PHOTON - ENERGY(NENER))/
     $  (ENERGY(NENER+1) - ENERGY(NENER))
        F2_S      =   f2s(NENER) +
     $  (f2s(NENER+1) - f2s(NENER))*(PHOTON - ENERGY(NENER))/
     $  (ENERGY(NENER+1) - ENERGY(NENER))
        dels   =   RHOS/ATWTS*AVOG*RADIUS/PI/2.0*(WAVE**2)*F1_s
        bets   =   RHOS/ATWTS*AVOG*RADIUS/PI/2.0*(WAVE**2)*F2_s
csrio	write (*,*) '<><><> substrate delta beta are:',dels,bets
c
c for the odd material
c
        F1_O      =   f1o(NENER) +
     $  (f1o(NENER+1) - f1o(NENER))*(PHOTON - ENERGY(NENER))/
     $  (ENERGY(NENER+1) - ENERGY(NENER))
        F2_O      =   f2o(NENER) +
     $  (f2o(NENER+1) - f2o(NENER))*(PHOTON - ENERGY(NENER))/
     $  (ENERGY(NENER+1) - ENERGY(NENER))
        delo   =   RHOO/ATWTO*AVOG*RADIUS/PI/2.0*(WAVE**2)*F1_o
        beto    =   RHOO/ATWTO*AVOG*RADIUS/PI/2.0*(WAVE**2)*F2_o
csrio	write (*,*) '<><><> odd delta beta are:',delo,beto
c
c for the even material
c
        F1_E      =   f1e(NENER) +
     $  (f1e(NENER+1) - f1e(NENER))*(PHOTON - ENERGY(NENER))/
     $  (ENERGY(NENER+1) - ENERGY(NENER))
        F2_E      =   f2e(NENER) +
     $  (f2e(NENER+1) - f2e(NENER))*(PHOTON - ENERGY(NENER))/
     $  (ENERGY(NENER+1) - ENERGY(NENER))
        dele   =   RHOE/ATWTE*AVOG*RADIUS/PI/2.0*(WAVE**2)*F1_e
        bete    =   RHOE/ATWTE*AVOG*RADIUS/PI/2.0*(WAVE**2)*F2_e
csrio	write (*,*) '<><><> odd delta beta are:',dele,bete
	return
	end

