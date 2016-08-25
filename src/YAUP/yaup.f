C  YAUP, version 1.3.3
C
C  Modification record:
C     12/92 - Initial version (without FFTs); calculates 
C	  spectral and angular distributions of the brightness and 
C	  flux through pinholes. 	
C     01/92 - FFT version.
C     04/92 - Bugs in module EXTRCT and trajectory calculation 
C	  fixed. Added degree-of-linear-polarization calculations. 	
C     06/92 - Add hanning windows to reduce spectral leakage;  
C	  expert-user options; power-load calculations. 	
C     10/92 - New input-file format.  Life is easy now. (v1.0)
C     04/93 - Include particle beam sizes. (v1.1) 	
C     09/93 - Fix a bug in the flux calculations with finite observation 
C	  distance. Some changes in the input file format. (v1.1b)
C     11/93 - Switch back to double precision variables for temp storage.
C	  Increases RAM requirements but improves "inter platform"
C	  reproducibility of the results. (v1.2/b/c)
C	  Check (but do not enforce) L/R symmetry criterion. (v1.2b)
C	  More flexible scratch file names. (v1.2c)
C     12/93 - User-specifiable scratch file names.  Bug in spline 
C	  interpolations fixed (INTGRD). Add GNUPLOT flag. (v1.2d)
C     4/94 - Rewrote the parser routines.  Fixed imsc bug.
C         Now compiles and runs on IBM RISC system/6000 too. (v1.3)
C     4/94 - minor output bug in mode 7 fixed. (v1.3.1)
C	     bug in integration in mode 7 fixed
C     7/94 - minor fixes here and there (v1.3.2)
C    11/13 - RJD added save statement to routine rzextr to avoid GFORTRAN not saving local variables. (v1.3.3)
C            See CRJD for additional minor changes.
C
C  DISCLAIMER
C
C  This program is freeware (but NOT public domain). You may freely copy 
C  and redistribute it. Permission is granted to modify the source for 
C  your own purposes, but NOT to redistribute the modifications without
C  permission of the author. If you use this program while doing scientific
C  research, please cite this program in the acknowledements of any resulting 
C  publications. 
C
C  There is absolutely no warranty on this program. The authors take
C  no responsibility for any damage caused by this program. The authors
C  take no responsibilty for time lost if incorrect or misleading results
C  are produced by this program. If a warranty is required by law where
C  you intend to use this software, permission to use this software there
C  is revoked. 
C
C  Blah, blah, blah...legalese...uuh.. 
C  For instructions see YAUP.DOC.
C 
C******************************************************************************
C 
C    YAUP - Yet Another (Useless) Undulator Program 
C 
C            Boyan Boyanov 
C            Grant Bunker 
C            Jay Lee 
C            Tim Morrison 
C            Department of Physics 
C            Illinois Institute of Technology 
C            Chicago, IL 60616 
C 
C    Direct bug reports, frustration and anger at 
C            Boyan Boyanov 
C            boyan@tmnxt1.iit.edu 
C            boyaboy@karl.iit.edu 
C            boyaboy@iitvax (Bitnet) 
C            (312) 567-3375 (voice) 
C            (312) 567-3396 (fax)
C
C    Reference
C	B. I. Boyanov, G. Bunker, J. M. Lee, and T. I. Morrison
C	"Numerical Modeling of Tapered Undulators"
C	Nucl. Instr. Meth. A339,  596-603, 1994
C 
C***************************************************************************** 
C 
	program yaup
C 
	implicit double precision (a-h,o-z) 
	include 'yaup.inc'
	character bugrep*20,vers*5,moddat*9
	data bugrep /'boyan@tmnxt1.iit.edu'/
CRJD November 15, 2013: changed version and modification date 
	data vers /'1.3.3'/
	data moddat /'NOV-15-13'/
C
C---------------------------------------------------------------------------- 
C   This open statement is needed with Absoft. 
C 
csrio	open(6,action='print') 
C---------------------------------------------------------------------------- 
C
	write(*,10) vers,moddat,bugrep 
10      format(/1x,60('-')//1x,'YAUP ',a,
     +	   ' - Yet Another (Useless) Undulator Program',
     +	       /1x,'Last modified on ',a
     +	       /1x,'Send bug reports to ',a) 
C
C  Make sure the source code is not messed up
C  Define I/O file names...
C
	call chkpar
C
C   read input from yaup.inp.  bmin is the minimum field strength 
C   (in tesla) and fKmin is the corresponding minimum value of the 
C   K-funtion.  efund is the fundamental energy of a untapered undulator 
C   with K=fKmin. 
C
	call getdat(ifile,param) 
	if (mode.eq.0.or.itraj.eq.1.or.itraj.eq.2) bmin=param 
	if (itraj.eq.0) efund=param 
C 
C  The old result file does not get overwritten 
C
	call opnout(bname,errmsg,iout) 
C 
C   determine angular limits for the zero-emittance calculations 
C 
	if (mode.ne.0) call bounds 
C 
C   write the input data to the output file 
C 
	call putinp(vers)
C
C   calculate the trajectory, if necessary
C
	if (itraj.ne.0.or.mode.eq.0) call trajc(bmin)
	if (mode.eq.0) stop
C
C   create/open a status file
C
	call status(irec,ix0,iy0)
C 
C   Calculate beam distributions 
C   rx and ry contain the beam divergencies (FFTd and packed). 
C 
	if (nsig.ne.0) then 
	   call gauss(sigu,dx,mx,rx) 
	   call gauss(sigv,dy,my,ry) 
	end if 
C
C   Calculate spectrum
C
	call specal(irec,ix0,iy0)
C
	close(iout)
	write(*,20) char(7) 
20      format(1x,a1/1x,60('-')//1x,'YAUP done'/) 
	stop
	end 
C
C********************** io.f ********************************
C
C       ------------------------------ 
	subroutine getdat(inpfn,param) 
C       ------------------------------ 
C 
C  Read input data.  The variable param stores either the minimum value 
C  of the undulator B-field (mode=1, itraj=1,2) or the value of the 
C  "fundamental" energy in eV (itraj=2) 
C 
	implicit double precision (a-h,o-z)
	include 'yaup.inc'
	character*(*) inpfn
C
C  Stuff for the parser routine
C
	parameter (xini=-555.0d0)
	parameter (nini=-555)
	character bra*1,ket*1,space*10,cchar*10
C
	bra='"'
	ket='"'
	space='='//char(9)
	nspace=2
	cchar='!;%#'
	ncchar=4
	call massag(xini,nini)
C
C   read input file.
C
	open(iinp,file=inpfn,status='old',err=800)
	call parser(iinp,bra,ket,space,nspace,cchar,ncchar,ier,errmsg)
	close(iinp)
	if (ier.ne.0) call leave(errmsg,' ',' ')
	call errchk(nini,xini)
	call gfiles(param)
C
C   Convert all angles to radians and distances to meters
C
	sigx=1.0d-3*sigx 
	sigy=1.0d-3*sigy 
	sigx1=1.0d-3*sigx1 
	sigy1=1.0d-3*sigy1 
	xpc=1.0d-3*xpc 
	ypc=1.0d-3*ypc 
	xps=1.0d-3*xps 
	yps=1.0d-3*yps
C
C   A finite particle beam size has an effect only for finite 
C   observation distances.  For finite observation distances 
C   the pinhole sizes are converted into radians.  Set dist to 1.0
C   it does not interfere with the size->angle conversions.
C
	if (dist.gt.zero) then
	   sigu=sqrt(sigx1**2+(sigx/dist)**2)
	   sigv=sqrt(sigy1**2+(sigy/dist)**2)
	   xpc=xpc/dist
	   ypc=ypc/dist
	   xps=xps/dist
	   yps=yps/dist
	else
	   dist=one
	   sigu=sigx1
	   sigv=sigy1
	end if
	return
C
C   Error handlers
C
800     call leave('getdat:: cannot open input file',inpfn,' ') 
	return
	end
C
C	----------------------------
	subroutine massag(xini,nini)
C	----------------------------
C
C   purpose: initialize everything to some stupid number.  this is used to
C	keep track of the variables which are not found in the input file.
C
	implicit double precision (a-h,o-z)
	include 'yaup.inc'
C
	period=xini
	emin=xini
	emax=xini
	ener=xini
	cur=xini
	sigx=xini
	sigy=xini
	sigx1=xini
	sigy1=xini
	dist=xini
	xpc=xini
	ypc=xini
	xps=xini
	yps=xini
	bw=xini
C
     	nper=nini
     	npts=nini
     	ne=nini
     	nxp=nini
     	nyp=nini
     	mode=nini
     	nsig=nini
     	itraj=nini
     	ixsym=nini
     	nhan=nini
	iexp=nini
	icont=nini
	isav=nini
	nres=nini
	ncrit=nini
	mfft1=nini
	ignu=nini
C
     	bfile='!@#$%'
     	tfile='!@#$%'
     	bname='!@#$%'
	return
	end
C
C	---------------------------------------------------
	subroutine token(word,nch,line,quote1,quote2,iflag)
C	---------------------------------------------------
C
C   purpose: called by SUBROUTINE PARSER, this procedure interprets
C	the next keyowrd
C
	implicit double precision (a-h,o-z)
	parameter (numkw=36)
	include 'yaup.inc'
	character*(*) word,line
	character*1 quote1,quote2
	character kword*15
	dimension kword(numkw),kwlen(numkw)
C
C   initalize keyword arrays.  one value per keyword.
C   the END keyword has no values associated with it
C
	save kword,kwlen,iquiet
	data kword /'PERIOD','NPER','NPTS','EMIN','EMAX','NE',
     +	            'ENERGY','CURRENT','SIGX','SIGY','SIGX1','SIGY1',
     +	            'DISTANCE','XPC','YPC','XPS','YPS','NXP','NYP',
     +	            'MODE','NSIG','TRAJECTORY','XSYM','HANNING',
     +	            'BFILE','TFILE','BASENAME',
     +	            'STATUS','UPDATE','RESOLUTION','NCRIT','MFFT','BW',
     +		    'GNUPLOT','QUIET','END'/
     	data kwlen / 6,4,4,4,4,2,
     +	             6,7,4,4,5,5,
     +	             8,3,3,3,3,3,3,
     +		     4,4,10,4,7,
     +	             5,5,8,
     +	             6,6,10,5,4,2,
     +	             7,5,3 /
     	data iquiet /0/
C
C   is WORD a KWORD?
C
	ier=0
	nk=1
	do while(nk.le.numkw.and.word(1:nch).ne.kword(nk)(1:kwlen(nk)))
	   nk=nk+1
	end do
	if (nk.gt.numkw.and.iquiet.eq.0) then
	   line='on line '//line
	   call leave('token:: unrecognized keyword',word,line)
	else if (nk.eq.numkw) then
	   iflag=-1
	   return
	end if
	
C
C   ok, WORD is a valid KWORD.  Handle file names first
C
	if (nk.ge.25.and.nk.le.27) then
	   n1=index(line,quote1)
	   line(n1:n1)=' '
	   n2=index(line,quote2)
	   if (n1*n2.eq.0) call 
     +		leave('token:: invalid file name format for keyword',
     +		kword(nk),'please enclose file names in double quotes')
	   if (n2-n1.eq.1) then
	      word=' '
	   else
	      word=line(n1+1:n2-1)
	   end if
	   if (nk.eq.25) then
	      bfile=word
	   else if (nk.eq.26) then
	      tfile=word
	   else if (nk.eq.27) then
	      bname=word
	   else
	      call leave('token:: INTERNAL ERROR 2',' ',' ')
	   end if
	   line=line(n2+1:)
C
	else
C
	nch=nword(line,word,ier,errmsg)
	if (ier.ne.0) call leave(errmsg,'on line',line)
	if (nch.eq.0) call 
     +	   leave('token:: cannot find value for keyword',kword(nk),' ') 
C
C   interpret the value depending on the KWORD
C
	if (nk.eq.1) then
	   period=greal(word,ier)
	   if (ier.ne.0) goto 800
	else if (nk.eq.2) then	
	   nper=igint(word,ier)
	   if (ier.ne.0) goto 800
	else if (nk.eq.3) then	
	   npts=igint(word,ier)
	   if (ier.ne.0) goto 800
	else if (nk.eq.4) then	
	   emin=greal(word,ier)
	   if (ier.ne.0) goto 800
	else if (nk.eq.5) then	
	   emax=greal(word,ier)
	   if (ier.ne.0) goto 800
	else if (nk.eq.6) then	
	   ne=igint(word,ier)
	   if (ier.ne.0) goto 800
	else if (nk.eq.7) then	
	   ener=greal(word,ier)
	   if (ier.ne.0) goto 800
	else if (nk.eq.8) then	
	   cur=greal(word,ier)
	   if (ier.ne.0) goto 800
	else if (nk.eq.9) then	
	   sigx=greal(word,ier)
	   if (ier.ne.0) goto 800
	else if (nk.eq.10) then	
	   sigy=greal(word,ier)
	   if (ier.ne.0) goto 800
	else if (nk.eq.11) then	
	   sigx1=greal(word,ier)
	   if (ier.ne.0) goto 800
	else if (nk.eq.12) then	
	   sigy1=greal(word,ier)
	   if (ier.ne.0) goto 800
	else if (nk.eq.13) then	
	   dist=greal(word,ier)
	   if (ier.ne.0) goto 800
	else if (nk.eq.14) then	
	   xpc=greal(word,ier)
	   if (ier.ne.0) goto 800
	else if (nk.eq.15) then	
	   ypc=greal(word,ier)
	   if (ier.ne.0) goto 800
	else if (nk.eq.16) then	
	   xps=greal(word,ier)
	   if (ier.ne.0) goto 800
	else if (nk.eq.17) then	
	   yps=greal(word,ier)
	   if (ier.ne.0) goto 800
	else if (nk.eq.18) then	
	   nxp=igint(word,ier)
	   if (ier.ne.0) goto 800
	else if (nk.eq.19) then	
	   nyp=igint(word,ier)
	   if (ier.ne.0) goto 800
	else if (nk.eq.20) then	
	   mode=igint(word,ier)
	   if (ier.ne.0) goto 800
	else if (nk.eq.21) then	
	   nsig=igint(word,ier)
	   if (ier.ne.0) goto 800
	else if (nk.eq.22) then
	   if (word(1:3).eq.'OLD'.and.nch.eq.3) then
	      itraj=0
	   else if (word(1:3).eq.'NEW'.and.nch.eq.3) then
	      itraj=1
	   else if (word(1:8).eq.'NEW+KEEP'.and.nch.eq.8) then
	      itraj=2
	   else
	      call leave('token:: invalid TRAJECTORY option:',word,' ')
	   end if
	else if (nk.eq.23) then	
	   ixsym=iyes(word,ier)
	   if (ier.ne.0) goto 800
	else if (nk.eq.24) then	
	   nhan=igint(word,ier)
	   if (ier.ne.0) goto 800
	else if (nk.eq.28) then
	   if (word(1:3).eq.'NEW'.and.nch.eq.3) then
	      icont=0
	   else if (word(1:8).eq.'NEW+KEEP'.and.nch.eq.8) then
	      icont=1
	   else if (word(1:3).eq.'OLD'.and.nch.eq.3) then
	      icont=2
	   else if (word(1:8).eq.'OLD+KEEP'.and.nch.eq.8) then
	      icont=3
	   else
	      call leave('token:: invalid STATUS option:',word,' ')
	   end if
	else if (nk.eq.29) then	
	   isav=igint(word,ier)
	   if (ier.ne.0) goto 800
	else if (nk.eq.30) then	
	   nres=igint(word,ier)
	   if (ier.ne.0) goto 800
	else if (nk.eq.31) then	
	   ncrit=igint(word,ier)
	   if (ier.ne.0) goto 800
	else if (nk.eq.32) then	
	   mfft1=igint(word,ier)
	   if (ier.ne.0) goto 800
	else if (nk.eq.33) then	
	   bw=greal(word,ier)
	   if (ier.ne.0) goto 800
	else if (nk.eq.34) then	
	   if (word(1:6).eq.'XLINES'.and.nch.eq.6) then
	      ignu=1
	   else if (word(1:6).eq.'YLINES'.and.nch.eq.6) then
	      ignu=2
	   else
	      call leave('token:: invalid GNUPLOT option:',word,' ')
	   end if
	else if (nk.eq.35) then	
	   iquiet=iyes(word,ier)
	   if (ier.ne.0) goto 800
	else
	   call leave('token:: INTERNAL ERROR 3',' ',' ')
	end if
C
	end if
	return
C
800	call leave('token:: cannot read value for keyword',
     +	   kword(nk),word)
	return
	end
C
C	-------------------------------------------
	double precision function greal(string,ier)
C	-------------------------------------------
C
	implicit double precision (a-h,o-z)
	character*(*) string
	read(string,*,iostat=ier) x
	greal=x
	return
	end
C
C	----------------------------------
	integer function igint(string,ier)
C	----------------------------------
C
	implicit integer (a-z)
	character*(*) string
	read(string,*,iostat=ier) ix
	igint=ix
	return
	end
C
C	---------------------------------
	integer function iyes(string,ier)
C	---------------------------------
C
	implicit double precision (a-h,o-z)
	character*(*) string
	ier=0
	if (string(1:3).eq.'YES') then
	   iyes=1
	else if (string(1:2).eq.'NO') then
	   iyes=0
	else
	   ier=1
	   iyes=0
	end if
	return
	end
C
C	--------------------------------
	subroutine leave(msg1,msg2,msg3)
C	--------------------------------
C
C   Print an error message
C
	character*(*) msg1,msg2,msg3
	write(*,*) char(7)
	n1=max(numchs(msg1),1)
	n2=max(numchs(msg2),1)
	n3=max(numchs(msg3),1)
	write(*,*) '*************'
	write(*,*) msg1(1:n1),' ',msg2(1:n2)
	if (n3.ne.1) write(*,*) msg3(1:n3)
	write(*,*) '*************'
	write(*,*)
	stop
	end
C
C	----------------------------
	subroutine errchk(nini,xini)
C	-----------------------------
C
	implicit double precision (a-h,o-z)
	include 'yaup.inc'
C
C   Some defaults
C
	if (nhan.eq.nini) nhan=0
	if (mode.eq.0) then
	   itraj=2
	   ne=1
	   emax=emin
	   xps=zero
	   yps=zero
	   xpc=zero
	   ypc=zero
	else if (mode.eq.1.or.mode.eq.5) then 
	   xps=zero 
	   yps=zero
	else if (mode.eq.2) then 
	   ne=1
	   emax=emin 
	end if
	ier=0
	if (bname(1:5).eq.'!@#$%') bname='yaup'
	nch=numchs(bname)
	if (nch.eq.0) call leave
     +	   ('errchk:: null BASENAMEs are not allowed',' ',' ')
C
C   missing data for non-optional input?
C
	ier=1
	ios=0
	if (period.eq.xini) then
	   errmsg='PERIOD'
	else if (nper.eq.nini) then
	   errmsg='NPER'
	else if (npts.eq.nini) then
	   errmsg='NPTS'
	else if (emin.eq.xini) then
	   errmsg='EMIN'
	else if (emax.eq.xini) then
	   errmsg='EMAX'
	else if (ne.eq.xini) then
	   errmsg='NE'
	else if (ener.eq.xini) then
	   errmsg='ENERGY'
	else if (cur.eq.xini) then
	   errmsg='CURRENT'
	else if (sigx.eq.xini) then
	   errmsg='SIGX'
	else if (sigy.eq.xini) then
	   errmsg='SIGY'
	else if (sigx1.eq.xini) then
	   errmsg='SIGX1'
	else if (sigy1.eq.xini) then
	   errmsg='SIGY1'
	else if (dist.eq.xini) then
	   errmsg='DISTANCE'
	else if (xpc.eq.xini) then
	   errmsg='XPC'
	else if (ypc.eq.xini) then
	   errmsg='YPC'
	else if (xps.eq.xini) then
	   errmsg='XPS'
	else if (yps.eq.xini) then
	   errmsg='YPS'
	else if (nxp.eq.nini) then
	   errmsg='NXP'
	else if (nyp.eq.nini) then
	   errmsg='NYP'
	else if (mode.eq.nini) then
	   errmsg='MODE'
	else if (nsig.eq.nini) then
	   errmsg='NSIG'
	else if (itraj.eq.nini) then
	   errmsg='TRAJECTORY'
	else if (ixsym.eq.nini) then
	   errmsg='XSYM'
	else if (bfile(1:5).eq.'!@#$%'.and.itraj.ne.0) then
	   errmsg='BFILE'
	else if (tfile(1:5).eq.'!@#$%'.and.itraj.eq.0) then
	   errmsg='TFILE'
	else
	   ier=0
	end if
	if (ier.ne.0) call 
     +	   leave('errchk:: cannot find data for keyword',errmsg,' ')
C
C   invalid data for non-optional input?
C
	ier=1
	if (period.le.zero) then
	   errmsg='PERIOD must be positive'
	else if (nper.lt.1.or.nper.gt.maxper) then
	   write(errmsg,10) maxper
10	   format('NPER is out of range [1,',i3,'].')
	else if (npts.lt.minpp.or.npts.gt.maxpp) then
	   write(errmsg,20) minpp,maxpp 
20	   format('NPTS out of range [',i2,',',i3,']') 
	else if (emin.le.zero) then
	   errmsg='EMIN must be positive'
	else if (emax.le.zero) then
	   errmsg='EMAX must be positive'
	else if (emin.gt.emax) then
	   errmsg='EMIN must be <= than EMAX.'
	else if (ne.gt.maxen) then 
	   write(errmsg,30) maxen
30	   format('NE is out of range [1,',i3,'].')
	else if (ener.le.zero) then
	   errmsg='ENERGY must be positive'
	else if (cur.le.zero) then
	   errmsg='CURRENT must be positive'
	else if (sigx.le.zero) then
	   errmsg='SIGX must be positive'
	else if (sigy.le.zero) then
	   errmsg='SIGY must be positive'
	else if (sigx1.le.zero) then
	   errmsg='SIGX1 must be positive'
	else if (sigy1.le.zero) then
	   errmsg='SIGY1 must be positive'
	else if (dist.lt.zero) then
	   errmsg='DISTANCE must be non-negative'
	else if (xps.lt.zero) then
	   errmsg='XPS must be non-negative'
	else if (yps.lt.zero) then
	   errmsg='YPS must be non-negative'
	else if (nxp.lt.0.or.nxp.gt.maxang.or.
     +	    nyp.lt.0.or.nyp.gt.maxang) then
 	   write(errmsg,40) maxang
40	   format('NXP and/or NYP are out of range [0,',i3,'].')
	else if (mode.lt.0.or.mode.gt.7) then
	   errmsg='MODE is out of range [0,7]'
 	else if (nsig.lt.0.or.nsig.gt.maxsig) then
	   write(errmsg,50) maxsig
50	   format('NSIG is out of range [0,',i1,'].')
	else if (ixsym.ne.0.and.ixsym.ne.1) then
	   write(errmsg,60)
60	   format('invalid or missing XSYM keyword.')
 	else if (nhan.gt.nper/2.or.nhan.lt.0) then
	   write(errmsg,70) nper/2
70	   format('HANNING is out of range [0,',i3,'].')
	else if (itraj.lt.0.or.itraj.gt.2) then
	   errmsg='invalid or missing TRAJECTORY keyword'
	else
	   ier=0
	end if
	if (ier.ne.0) call leave('errchk::',errmsg,' ')
C
C  "Expert" keywords
C
	iexp=0
	if (icont.ne.nini.or.isav.ne.nini.or.nres.ne.nini.or.
     +	    ncrit.ne.nini.or.mfft1.ne.nini.or.bw.ne.xini) iexp=1
	if (icont.eq.nini) icont=icontd
	if (isav.eq.0.or.isav.eq.nini) isav=isavd
        if (isav.lt.0) isav=0
	if (nres.eq.0.or.nres.eq.nini) nres=nresd
	if (ncrit.eq.0.or.ncrit.eq.nini) ncrit=ncritd
	if (mfft1.eq.0.or.mfft1.eq.nini) mfft1=mfftd
	if (bw.eq.zero.or.bw.eq.xini) bw=bwd
	if (ignu.eq.nini) ignu=0
C
	ier=1
	if (icont.lt.0.or.icont.gt.3) then
	   errmsg='invalid STATUS keyword.'
	else if (isav.lt.0) then
	   errmsg='invalid UPDATE keyword.'
	else if (nres.le.0) then
	   errmsg='invalid RESOLUTION keyword.'
	else if (ncrit.le.0) then
	   errmsg='invalid NCRIT keyword.'
	else if (mfft1.lt.2) then
	   errmsg='invalid MFFT keyword.'
	else
	   ier=0
	end if
	if (ier.ne.0) call leave('errchk::',errmsg,' ')
C
	return
	end
C
C	------------------------
	subroutine gfiles(param)
C	------------------------
C
C   purpose: to read data from external files
C	param - either the minimum value of the undulator B-field 
C		(mode=1, itraj=1,2) or the value of the "fundamental" 
C		energy in eV (itraj=2)
C   external subroutines: leave
C   written: boyan boyanov, 1/91
C
	implicit double precision (a-h,o-z)
	include 'yaup.inc'
	parameter (small=0.0001d0)
C	   
C   Test output file basename and create output names
C
	nch=numchs(bname)
     	stks1=bname(1:nch)//'.txt'
	open(imisc,file=stks1,status='unknown',iostat=ier)
	if (ier.ne.0) call leave
     +	   ('errchk:: invalid BASENAME:',bname,' ')
     	close(imisc,status='delete')
	stks1=bname(1:nch)//'.s0'
	stks2=bname(1:nch)//'.s1'
	cfile=bname(1:nch)//'.cfg'
C 
C  if the trajectory has been precalculated, read it from a file 
C
	if (itraj.eq.0.and.mode.ne.0) then
	   write(*,20) tfile(1:numchs(tfile))
20	   format(/1x,60('-')//1x,'Reading trajectory from file ',a) 
	   open(imisc,file=tfile,err=805,form='unformatted',
     +	      status='old') 
	   read(imisc,err=810,end=810) per1,nper1,npts1,
     +	      ener1,param
C
C   Make sure that the stored hardware parameters are the same as 
C   as those in the input file
C 
	   if (abs(per1-period).gt.small.or.nper.ne.nper1.or.
     +	       npts.ne.npts1.or.abs(ener1-ener).gt.small) goto 820
C
C   total number of points and sampling stepsize
C
	   nptst=nper*npts+1
	   dz=period/dble(npts)
C
C   read the trajectory
C 
	   read(imisc,err=810,end=810) 
     +	      (ct(i),x(i),betax(i),betaz(i),i=1,nptst) 
	   close(imisc) 
	   do i=1,nptst 
	      z(i)=dble(i-1)*dz 
	   end do 
C 
C  otherwise read field strength (in tesla) as a function of the longitudinal 
C  coordiante z and temporarily store it in the array intended for betaz. 
C  the phase errors (the field is assumed to be a modulated sinusoid) are 
C  read in z(i) 
C 
	else if (itraj.eq.1.or.itraj.eq.2.or.mode.eq.0) then 
	   write(*,30) bfile(1:numchs(bfile))
30	   format(/1x,60('-')//1x,'Reading B-field distribution ',
     +	        '(in tesla) from file ',a) 
	   open(imisc,file=bfile,err=800,form='unformatted',
     +		status='old') 
	   read(imisc,err=810,end=810) per1,nper1,npts1,param
	   ener1=ener
C
C   Make sure that the stored hardware parameters are the same as 
C   as those in the input file
C 
	   if (per1.ne.period.or.nper.ne.nper1.or.
     +	       npts.ne.npts1.or.ener1.ne.ener) goto 820
C
C   total number of points and sampling stepsize
C
	   nptst=nper*npts+1
	   dz=period/dble(npts)
C
C   read the b-field distribution.   The amplitude function is stored
C   stored in BETAZ and the phase errors go in Z.
C
	   read(imisc,err=810,end=810) (betaz(i),z(i),i=1,nptst) 
	   close(imisc) 
C
C   test output file name (if necessary)
C
	   if (itraj.eq.2) then
	      open(imisc,file=tfile,status='unknown',err=805)
	      close(imisc,status='delete')
	   end if
	end if 
	return
C
C   Error handlers
C
800	call leave('gfiles:: cannot open BFILE',' ',bfile)
	return
805	call leave('gfiles:: cannot open TFILE',' ',tfile)
	return
810	call leave('gfiles:: attempt to read beyond EOF of input file ',
     +		' ',' ')
	return
820     write(errmsg,825) per1,nper1,npts1,ener1 
825     format('expecting PERIOD,NPER,NPTS,ENERGY = ',
     +	        f5.1,',',i4,',',i4,',',f6.2)
     	call leave('gfiles::',errmsg,' ')
	return 
	end
C
C	------------------------------------
	subroutine opnout(fname,errmsg,iout)
C	------------------------------------
C
C   Open an output file without overwriting any existing ones.
C   The filename is of the form BNAME-n.OUT, n=0-9.
C 
	implicit integer (a-z)
	character*(*) fname,errmsg
	character*80 name
	logical exists 
C
	nf=numchs(fname)
	name=fname(1:nf)//'-0.out' 
	inquire(file=name,exist=exists) 
	if (exists) then
	   n=0
	   do while (n.lt.10.and.exists)
	      n=n+1
	      write(name,10) fname(1:nf),n
10	      format(a,'-',i1,'.out')
	      inquire(file=name,exist=exists) 
	   end do
	   if (n.ge.10) then
	      write(errmsg,15) fname(1:nf)
15	      format('please delete ',a,'-*.* files')
	      call leave('opnout:: cannot open an output file',
     +		         ' ',errmsg)
     	   end if
	end if 
C
C   The output is stored in the file NAME
C
	n=numchs(name)
	write(*,20) name(1:n)
20	format(1x,'Writing output to ',a)
	open(iout,file=name,status='new')
	return
	end
C	 
C       -----------------------
	subroutine putinp(vers) 
C       -----------------------
C 
C  dump input on output unit (similar to urgent) 
C 
	implicit double precision (a-h,o-z) 
	include 'yaup.inc'
	character*(*) vers
C
	if (dist.eq.one) then
	   pdist=zero
	else
	   pdist=dist
	end if
C
C  Input parameters
C
	write(iout,10) vers
10      format(15x,10('*'),' YAUP ',a,1x,10('*')) 
	write(iout,20) period,nper,npts 
20      format(/1x,'undulator :  period (cm) = ',f5.1,3x, 
     +	           'nper = ',i3,3x,'npts = ',i4) 
	write(iout,30) emin,emax,ne 
30      format(/1x,'photon energy range (eV) : emin = ',f9.2,3x, 
     +	           'emax = ',f9.2,3x,'ne = ',i4) 
	write(iout,40) ener,cur 
40      format(/1x,'electron beam :  energy (GeV) = ',f6.3, 
     +	        3x,'current  (A) = ',f6.3) 
	write(iout,50) 1.d3*sigx,1.d3*sigy,1.d3*sigx1,1.d3*sigy1 
50      format(18x,'sigx  (mm)   = ',f6.4,3x,'sigy  (mm)   = ',f6.4
     +	      /18x,'sigx1 (mrad) = ',f6.4,3x,'sigy1 (mrad) = ',f6.4) 
	write(iout,60) pdist,1.d3*xpc*dist,1.d3*xps*dist,nxp,
     +	                     1.d3*ypc*dist,1.d3*yps*dist,nyp 
60      format(/1x,'pinhole :  distance (m) = ',f7.2
     +	       /12x,'xpc  (mrad/mm) = ',f7.3, 
     +	         3x,'xps (mrad/mm) = ',f7.3,3x,'nxp = ',i4
     +	       /12x,'ypc  (mrad/mm) = ',f7.3,
     +	         3x,'yps (mrad/mm) = ',f7.3,3x,'nyp = ',i4) 
	write(iout,80) mode,nsig,itraj,ixsym,nhan 
80      format(/1x,'parameters  :  mode  = ',i1,3x,'nsig = ',i2, 
     +	        3x,'itraj  = ',i1/16x,'ixsym = ',i1,3x,'nhan = ',i2) 
C
C  "Expert" parameters (if any)
C
	if (iexp.ne.0) then
	   write(iout,90) icont,isav,nres,ncrit,mfft1,bw*1.d2
90	   format(/1x,'EXPERT options :'
     +	/3x,'icont = ',i3,3x,'isav = ',i3,3x,'nres = ',i3
     +	/3x,'ncrit = ',i3,3x,'mfft = ',i3,3x,'bw   = ',f7.3)
	end if
C
C   Filenames
C
	if (itraj.eq.0) then 
	   write(iout,100) tfile(1:numchs(tfile)) 
100	   format(/1x,'trajectory from file ',a) 
	else 
	   write(iout,110) bfile(1:numchs(bfile))
110	   format(/1x,'B-field ditribution from file ',a) 
	end if 
	if (mode.eq.0.or.itraj.eq.2) then 
	   write(iout,120) tfile(1:numchs(tfile))
120	   format(1x,'trajectory saved to file ',a) 
	end if 
	write(iout,130) 
130     format(/1x,'units : ' 
     +         /3x,'x, y             - mrad(mm)'
     +         /3x,'energy           - eV'
     +	       /3x,'brightness       - ph/s/bandpass/mrad^2(/mm^2)' 
     +	       /3x,'flux             - ph/s/bandpass'
     +	       /3x,'spectral power   - watts/eV/mrad^2(/mm^2)'
     +	       /3x,'power            - watts/mrad^2(/mm^2)') 
	if (mode.eq.0) then 
	   write(iout,140) 
140	   format(//1x,'Trajectory calculation only.') 
	else if (mode.eq.1) then
	   write(iout,150) 'Spectrum of the brightness at XPC,YPC.'
	else if (mode.eq.2) then
	   write(iout,150) 'Angular distribution of the brightness '//
     +	          'at EMIN.'
	else if (mode.eq.3) then
	   write(iout,150) 'Angular and spectral distribution of '//
     +	          'the brightness.'
150	   format(//1x,a//6x,'energy',8x,'x',12x,'y', 
     +	            10x,'brightness',5x,'polarization') 
	else if (mode.eq.4) then 
	   write(iout,160) 
160	   format(//1x,'Flux calculations.'//5x,'energy',9x,'flux')
	else if (mode.eq.5) then
	   write(iout,170) 'Spectral power at XPC,YPC.'
	else if (mode.eq.6) then
	   write(iout,170)'Angular distribution of the spectral power.'
170	   format(//1x,a//6x,'energy',8x,'x',12x,'y',12x,'power') 
	else if (mode.eq.7) then
	   write(iout,180) 'Angular distribution of the power.'
180	   format(//1x,a//12x,'x',12x,'y',8x,'power') 
	end if 
	return 
	end

C
C	----------------------------------------------
	subroutine bmake(fname,iounit,nbytes,lword,ny)
C	----------------------------------------------
C
C   purpose: create a temporary storage file
C
	implicit double precision (a-h,o-z)
	character*(*) fname
	logical exists
C
	inquire(file=fname,exist=exists) 
	if (exists) call leave('bmake:: found old version of file', 
     +		fname,'delete if not needed, rename otherwise') 
C
C   Open and close the dump files.  This puts empty copies on the disk
C   that are later on opened as 'OLD' by SAVDMP.
C
	open(iounit,file=fname,form='unformatted',access='direct', 
     +		recl=(nbytes/lword)*ny,status='new',err=900)
	close(iounit,status='keep')
	return
C
900	call leave('bmake:: cannot create file',fname,' ')
	return
	end
C
C       ---------------- 
	subroutine bopen 
C       ----------------
C 
C   open dump files for reading (status='OLD'). 
C 
	implicit double precision (a-h,o-z)
	include 'yaup.inc'
C
	open(idmp1,file=stks1,form='unformatted',access='direct', 
     +	        recl=(nbytes/lword)*ny,status='old',err=900)
	open(idmp2,file=stks2,form='unformatted',access='direct', 
     +	        recl=(nbytes/lword)*ny,status='old',err=910)
	open(ista,file=cfile,form='unformatted',access='direct',
     +	        recl=npar*(4/lword),status='old',err=920)
	return 
C
900	call leave('bopen:: cannot open/find exisitng file',
     +		stks1,' ')
	return 
910	call leave('bopen:: cannot open/find exisitng file',
     +		stks2,' ')
	return 
920	call leave('bopen:: cannot open/find exisitng file',
     +		cfile,' ')
	return 
	end
C 
C       ----------------------- 
	subroutine bread(iener) 
C       ----------------------- 
C 
C   read  records from dump files. the records are not consecutive, 
C   but are separated by intervals of ne records.  the data is stored 
C   in single precision, so it is read back in single precision. 
C 
	implicit double precision (a-h,o-z) 
	include 'yaup.inc' 
C 
C   s0 and s1 contain the zero-emittance Stoke's coefficients 
C
	dimension s0(nconv,nconv),s1(nconv,nconv) 
	equivalence (fx,s0),(fz,s1) 
C 
	do ix=1,nx 
	   irec=(ix-1)*ne+iener 
	   read(idmp1,rec=irec,err=90)(s0(ix,iy),iy=1,ny)
	   read(idmp2,rec=irec,err=90)(s1(ix,iy),iy=1,ny)
	end do 
	return 
C 
90      call bclose('keep',idmp1,idmp2,ista)
	call leave('bread:: error while reading temp storage files',
     +		' ',' ')
	return 
	end  
C 
C       ----------------------------------------- 
        subroutine bwrite(idmp,irec,ix,iy,e0,ne0) 
C       ----------------------------------------- 
C 
C   save intermediate results at a specific angular position (XC,YC) 
C   in a binary file.  the availble results are interpolated over the 
C   required energy range, and then stored  temporarily in the arrays S0 
C   and S1.  every NY points the array  is dumped to the binary file. 
C   the maximum values of NX and NY are limited only by NCONV, and therefore 
C   S0 and S1 should be able to acomodate these extremities.  They may not 
C   be made equivalent with fx and/or fz. 
C 
C   Note : see the remarks in FFTPTS about the interpolation. 
C 
        implicit double precision (a-h,o-z)
        include 'yaup.inc' 
	parameter (tiny=1.0d-7) 
C 
C   The zero-emittance brightness/power along x and y over the FFT energy 
C   mesh are stored by  BRIGHT in fx and fz, repspectively. 
C   S0 is an improvised buffer.  Each row of S0 contains the 
C   brightness at position xc=x0+(ix-1)*dx, yc=y0+(iy-1)*dy over the 
C   requested energy range.  each column contains exactly ny CONSECUTIVE 
C   points (one row of the pinhole mesh).  The total number of buffers that 
C   are needed is nx.  the format of S1 is the same, except that S1 contains 
C   the other Stoke's coefficient. S0 and S1 are single precision (they are 
C   used for data transfer only and do not appear to affect the final result). 
C   EKNOTS, WA, BRIX, BRIY are work arrays for the splining procedures.  
C 
        dimension s0(nconv,maxen),s1(nconv,maxen) 
        dimension brix(maxfft),briy(maxfft) 
        equivalence (fx,brix),(fz,briy)
        save s0,s1
C
C   Interpolate the current data over the user-requested energy mesh.  
C   Since the zero-emittance spectrum may vary rapidly over several 
C   orders of magnitude we interpolate  the log, calculate the 
C   interpolated values for the log and the invert to find the value 
C   of the brightness/spectral power.  We use a low-order polynomial 
C   (linear) because the radius of curvature of the zero-emittance 
C   spectrum may be rather tight (and usually is). 
C 
	de=zero 
        if (ne.ne.1) de=(emax-emin)/dble(ne-1) 
	nsav=ne
        do i=1,ne0 
           brix(i)=log(brix(i)) 
        end do 
C 
C   The y component of the brightness/power is zero in the horizontal plane 
C   and a log interpolation will fail.  Take care of this here. 
C 
        if (abs(yn).gt.tiny) then 
           xlogy=zero 
           do i=1,ne0 
              briy(i)=log(briy(i)) 
           end do 
        else 
           xlogy=one 
        end if 
C 
C   The interpolation knots are equidistant. 
C 
        do ie=1,ne 
           ec=emin+dble(ie-1)*de 
           ne1=int((ec-e0)/defft)+1 
           ne2=ne1+1 
           e1=e0+dble(ne1-1)*defft 
           e2=e1+defft 
           h=e2-e1 
C 
C   x- polarization 
C 
           slope=(brix(ne2)-brix(ne1))/h 
           bx=exp(brix(ne1)+slope*(ec-e1)) 
C 
C   and y-polarization 
C 
           slope=(briy(ne2)-briy(ne1))/h 
           by=exp(briy(ne1)+slope*(ec-e1))-xlogy 
C 
C   Stoke's coefficients S0 and S1.
C
           s0(iy,ie)=(bx+by) 
           s1(iy,ie)=(bx-by)
        end do
C
C   If mode 7 do a stupid trapezoidal integration of the 
C   signal.  Assuming you can interchange the order of integrations
C   (energy and spatial) it does not matter which is done first.
C
	if (mode.eq.7) then
	   nsav=1
	   sum0=zero
	   sum1=zero
	   do i=2,ne-1
	      sum0=sum0+s0(iy,i)
	      sum1=sum1+s1(iy,i)
	   end do
	   sum0=de*(sum0+half*(s0(iy,1)+s0(iy,ne)))
	   sum1=de*(sum1+half*(s1(iy,1)+s1(iy,ne)))
	   s0(iy,1)=sum0
	   s1(iy,1)=sum1
	end if 
C
C   Is it time to dump buffers and update status file? 
C   each row of pinhole mesh at a specific energy, i.e. each column of sb, 
C   is stored in a separate record of the direct access file.
C   The dump and status files are opened/close before/after each ISAV I/O
C   operations to keep them up to date. 
C 
        if (iy.eq.ny) then
	   idmp=idmp+1
C
C   Open the dump file before row 1, isav+1, 2*isav+1, etc of the PH mesh.
C   Always open the files the first time around  - this avoids problems when
C   ICONT=2,3
C
	   if (mod(ix-1,isav).eq.0.or.idmp.eq.1) call bopen
           do ie=1,nsav 
              irec=irec+1 
              write(idmp1,rec=irec) (s0(j,ie),j=1,ny) 
              write(idmp2,rec=irec) (s1(j,ie),j=1,ny) 
           end do 
	   write(ista,rec=2) ix,iy,irec
C
C   and close it after row isav,2*isav,3*isav,  etc or when done.
C
	   if (mod(ix,isav).eq.0.or.ix.eq.nx) 
     +		call bclose('keep',idmp1,idmp2,ista)
        end if 
        return 
        end 
C 
C	-----------------------------------
	subroutine bclose(stat,io1,io2,io3) 
C	-----------------------------------
C 
	implicit double precision (a-h,o-z)
	character*(*) stat
	close(unit=io1,status=stat) 
	close(unit=io2,status=stat)
	close(unit=io3,status=stat) 
	return 
	end 
C 
C       -------------------------------- 
	subroutine report(d,x,y,ic,itot) 
C       -------------------------------- 
C 
	implicit double precision (a-h,o-z) 
C
C   Do not use 3pf8.4 formats because of stupid compilers
C
	perc=1.d2*dble(ic)/dble(itot)
	write(*,10) 1.d3*x*d,1.d3*y*d,perc 
10      format('+',3x,f8.4,t18,f8.4,t31,f6.1,t70,' ') 
	return 
	end 
C 
C       -------------------------------------- 
	subroutine putdat(flux,ephot,xmin,ymin) 
C       -------------------------------------- 
C 
	implicit double precision (a-h,o-z) 
	include 'yaup.inc' 
C 
	dimension bri(nconv,nconv),degp(nconv,nconv) 
	equivalence (fx,bri), (fz,degp) 
C 
C   Save results in output file. 
C 
	if (mode.ge.1.and.mode.le.3) then 
	   do ix=1,nxp
	      xc=xmin+dble(ix-1)*dx 
	      do iy=1,nyp 
	         yc=ymin+dble(iy-1)*dy 
	         write(iout,10) ephot,1.d3*dist*xc,1.d3*dist*yc, 
     +	                  bri(ix,iy),degp(ix,iy) 
10	        format(3x,f9.1,2(2x,f10.5),5x,1pe14.6,5x,0pf8.5) 
	      end do 
	      if (ignu.eq.2) write(iout,*)
	   end do 
	   if (ignu.eq.1) write(iout,*)
	else if (mode.eq.4) then 
	   write(iout,20) ephot,flux 
20	   format(3x,f9.2,3x,1pe14.6)
	else if (mode.ge.5.and.mode.lt.7) then
	   do ix=1,nxp
	      xc=xmin+dble(ix-1)*dx 
	      do iy=1,nyp 
	         yc=ymin+dble(iy-1)*dy 
	         write(iout,30) ephot,1.d3*dist*xc,1.d3*dist*yc,
     +			bri(ix,iy) 
30	         format(3x,f9.1,2(2x,f10.5),1x,1pe14.6) 
	      end do
	      if (ignu.eq.2) write(iout,*)
	   end do
	   if (ignu.eq.1) write(iout,*)
	else if (mode.eq.7) then
	   do ix=1,nxp
	      xc=xmin+dble(ix-1)*dx 
	      do iy=1,nyp 
	         yc=ymin+dble(iy-1)*dy 
	         write(iout,40) 1.d3*dist*xc,1.d3*dist*yc,bri(ix,iy) 
40	         format(3x,2(2x,f10.5),1x,1pe14.6) 
	      end do
	   end do
	else
	   pause 'What the hell are you doing here?'
	   stop   	
	end if 
	return 
	end 
C
C*************************** parser.f ***************************
C
C	----------------------------------------------------
	subroutine parser(iounit,quote1,quote2,space,nspace,
     +		cchar,ncchar,ier,errmsg)
C	----------------------------------------------------
C
C   purpose: given a file connected to unit IOUNIT read the file one line
C	at a time, get rid of the comment portions, replace all word
C	delimiters with space chars, convert the line to uppercase, except
C	for case-sensitive portions (e.g. filenames), parse the line into
C	words and call a user-supplied subroutine TOKEN which processes
C	the rest of the line
C   variables:
C	iounit - the unit to which the file is connected
C	quote1,quote2 - everything between quote1/quote2 pairs is 
C		case-sensitive.  MUST BE OF TYPE CHAR*1.
C	space - characters considered to be word delimiters
C	nspace - number of chars passed with SPACE
C	cchar - comment chars. everything after one of these is ignored
C	ncchar - number of characters passed with CCHARR
C	ier - error flar (error occured when <>0)
C	errmsg - error text when ier <>0
C   externals: triml,wspace,upcase,numchs,token
C	where TOKEN is a user-supplied routine of the form
C	
C	SUBROUTINE TOKEN(WORD,NCH,LINE,QUOTE1,QUOTE2,IFLAG)
C	CHARACTER*(*) WORD,LINE
C	CHARACTER*1 QUOTE1,QUOTE2
C	.
C	.
C	RETURN
C	END
C
C	where WORD is the token to be interpreted and LINE is the
C	part of the input line immediately after the token.  If the user
C	routine sets IFLAG to a negative value the execution of this 
C	subroutine is immediately terminated. 
C
	implicit integer (a-z)
	character*(*) space,cchar,errmsg
	character*1 quote1,quote2 
	character*255 line,ucline,word
C
	word=' '
	iflag=0
	ier=0
	do while (.true.)
C
C   read a text line, get rid of comments, replace all word delimiters
C   with white space, convert to upper case
C
	   read(iounit,10,err=800,end=20) line
10	   format(a)
	   ucline=line
C
	   call triml(ucline,cchar,ncchar)
	   call wspace(ucline,space,nspace)	   
	   call upcase(ucline,quote1,quote2,ier,errmsg)
	   if (ier.ne.0) then
	      nch=max(1,numchs(errmsg))
	      errmsg=errmsg(1:nch)//' on line: '//line
	      return
	   end if
C
C   get the new word and pass it on together with the rest of LINE
C   to the user-supplied subroutine
C
	   nch=nword(ucline,word,ier,errmsg)
	   if (ier.ne.0) then
	      nch=max(1,numchs(errmsg))
	      errmsg=errmsg(1:nch)//' on line '//line
	      return
	   end if
	   do while (nch.gt.0)
	      call token(word,nch,ucline,quote1,quote2,iflag)
	      if (iflag.lt.0) return
	      nch=nword(ucline,word,ier,errmsg)
	      if (ier.ne.0) then
	         nch=max(1,numchs(errmsg))
	         errmsg=errmsg(1:nch)//' on line: '//line
	         return
	      end if
	   end do
	end do
20	return
c
800	ier=1
	inquire(unit=iounit,name=word)
	errmsg='parser:: error while reading from file'//word
	return
	end
	
C
C	-----------------------------------
	subroutine triml(line,cchar,ncchar)
C	-----------------------------------
C
C   purpose: delete the comment portions from LINE.  Everything
C	that follows a char found in cchar is considered a comment
C   externals: none
C   written: b. boyanov 4/94
C
	implicit integer (a-z)
	character*(*) line,cchar
	do i=1,ncchar
	   n=index(line,cchar(i:i))
	   if (n.gt.0) line(n:)=' '
	end do
	return
	end
C
C	---------------------------------
	subroutine wspace(line,space,nsp)
C	---------------------------------
C
C   replace all characters in line that occur in space with ' '
C
	implicit integer (a-z)
	character*(*) line,space
C
	do i=1,nsp
	   n=index(line,space(i:i))
	   do while (n.gt.0)
	      line(n:n)=' '
	      n=index(line,space(i:i))
	   end do
	end do
	return
	end
C
C     ------------------------------------------
      subroutine upcase(line,bra,ket,ier,errmsg)
C     ------------------------------------------
C
C  purpose: convert a string line to upper case, with the exception of
C	portions enclosed in BRA/KET pairs.  THIS SUBROUTINE IS ASCII-
C	DEPENDENT.
C
	implicit integer (a-z)
	parameter (aa=97,zz=122,shift=32)
	character*(*) line,errmsg
	character*1 bra,ket
C
	ier=1
	nch=len(line)
	n3=0
	do while (n3.lt.nch)
	   left=index(line(n3+1:),bra)
	   right=left+index(line(n3+1+left:),ket)
C
C   check for unpaired bra's or ket's
C
	   if (left.eq.0.and.right.ne.0) then
	      errmsg='upcase:: unpaired right "quote" ('//ket//')'
	      return
           else if (left.ne.0.and.right.eq.0) then
	      errmsg='upcase:: unpaired left "quote" ('//bra//')'
	      return
	   else if (left.eq.right.and.left*right.ne.0) then
	      errmsg='upcase:: odd number of "quotes" ('//bra//')'
	      return
           else if (left.eq.0.and.right.eq.0) then
	      n1=n3
	      n2=nch
	      n3=nch
	   else
	      n1=n3
	      n2=n1+left-1
	      n3=n1+right
	   end if
	   do i=n1+1,n2
	      ch=ichar(line(i:i))
	      if (ch.ge.aa.and.ch.le.zz) ch=ch-shift
	      line(i:i)=char(ch)
	   end do
	end do
	ier=0
	return
	end
C
C	-----------------------------
	integer function numchs(line)
C	-----------------------------
C
C   purpose: returns the postion of the last non-blank char in a line
C   external subroutines called: none
C   written: boyan boyanov, 8/93
C   
	implicit integer (a-z)
	character*(*)line
	n=len(line)
	do while(n.ge.1.and.line(n:n).eq.' ')
	   n=n-1
	end do
	numchs=n
	return
	end
C
C	--------------------------------------------
	integer function nword(line,word,ier,errmsg)
C	--------------------------------------------
C
C   purpose: get the next WORD from LINE.  WORDs are delimited by spaces.
C	LINE is returned without the parsed WORD.  The function returns
C	the number of characters in WORD 
C
	implicit integer (a-z)
	character*(*) line,word,errmsg
C
	ier=0
	nch=len(line)
	maxch=len(word)
	n1=1
	do while (n1.le.nch.and.line(n1:n1).eq.' ')
	   n1=n1+1
	end do
	n2=n1+index(line(n1:),' ')-2
	if (n2.ge.n1) then
	   nch=n2-n1+1
	   word=line(n1:n2)
	   line=line(n2+1:)
	else
	   nch=0
	   word=' '
	   line=' '
	end if
	if (nch.gt.maxch) then
	   ier=1
	   nch=0
	   errmsg='nword:: long word: '//word
	end if
	nword=nch
	return
	end
C
C*********************************** traj.f *******************************
C 
C       ----------------------
	subroutine trajc(bmin)
C       ----------------------
C 
C   Integrate the eqs. of motion to find the trajectry as a function 
C   of the longitudinal coordinate z, and store it in the common block 
C   trj. The boundary conditions x(0)=0, x(ulen)=0 are imposed on the 
C   solution. 
C 
	implicit double precision (a-h,o-z) 
	include 'yaup.inc' 
	parameter (restm=0.5109077090d-3) 
     	parameter (pihalf=1.570796326794896619231322d0)
	parameter (nshot=4)
	parameter (npersh=2)
	parameter (maxtry=nshot*npersh+1) 
	dimension wa(maxpts),d2bdz2(maxpts)
	equivalence (wa,fx),(d2bdz2,fz) 
C
C   Left and right hanning-window statement functions.
C
	hanl(zz,zend,deltaz)=(sin(pihalf*(zz-zend)/deltaz))**2
	hanr(zz,zbegin,deltaz)=(cos(pihalf*(zz-zbegin)/deltaz))**2
C 
C   first calculate some hardware constants and stuff the "local"
C   common block TRJTMP. 
C 
	gamma=ener/restm 
	beta2=one-one/(gamma*gamma) 
	beta=sqrt(beta2) 
	xk0=twopi/period 
	delta=beta*xk0/gamma
C 
C  then some stuff for the ODE integrator.  MAXSTP is the maximum number 
C  of steps that the ODE integrator may take, EPS is roughly equal to the 
C  fractional error desired in the solution, DZMIN is the minimum step 
C  that the integrator is allowed to take, and dz0 is an initial guess 
C  on the stepsize needed to achieve an accuracy EPS in the first step. 
C 
	maxstp=50 
	eps=0.01d0 
	dzmin=dz/dble(maxstp) 
	dz0=two*dzmin 
C 
C  Write info to screen 
C 
	write(*,10) dble(nper)*period 
10      format(/1x,60('-')/ 
     +/1x,'Beginning calculation of electron trajectory as a function' 
     +/1x,'of the longitudinal  coordinate. The total undulator length' 
     +/1x,'is ',f6.1,' centimeters.' 
     +//1x,t6,'Pass #',t17,'z(m)',t26,'% done',t40,'error (cm)'/) 
C 
C   calculate the maximum "fundamental" energy in eV (corresponding to 
C   the minimum field amplitude) and minimum transverse deflection 
C   amplitude (the final error from the shooting method is printed as a 
C   percentage of this value). 
C 
	fKmin=0.934d0*bmin*period 
	efund=949.0d0*ener*ener/(period*(1+fKmin*fKmin/two)) 
	dflmin=beta*fKmin/gamma 
C 
C   initial conditions. betax(1) is a free parameter.  Take an 
C   intelligent guess on betax(1), and "shoot" to find the 
C   "actual" value. 
C 
	ct(1)=zero
	x(1)=zero 
	betax0=beta*(0.934d0*betaz(1)*period)/gamma 
	dbetax=0.1d0*betax0 
	betax1=betax0+dbetax 
C 
	ipass=0 
	small=0.001d0 
	betax(1)=betax0 
	x(nptst)=dflmin 
	ier=0 
C 
	do while (abs(x(nptst)).gt.small*dflmin.and.ipass.lt.maxtry) 
	   betax0=betax(1) 
	   call shoot(ipass,betax0,dbetax,betax1,small,ier,errmsg, 
     +	              dz,dz0,dzmin,eps,maxstp) 
	   if (ier.ne.0) call leave('trajc::',errmsg,' ') 
	   betax(1)=betax1 
	end do 
C
C   Stuff the z-array.  The phase errors are overwritten here
C
	do i=1,nptst 
	   z(i)=dble(i-1)*dz 
	end do 
C
C  Check for L/R symmetry.  See criterion in B. I. Boyanov et al...
C  The field amplitude is stored in betaz.  Spline it to evaluate the
C  second derivatives d2bdz2.  Omit the endpoints because the second
C  derivatives there are zero by definition
C
	call spline(z,betaz,nptst,d2bdz2,wa)
	ratio=huge
	do i=2,nptst-1
	   del=abs(d2bdz2(i)/betaz(i))
	   if (ratio.gt.del) ratio=del
	end do
	hcgnl=hc*gamma*gamma/dble(nper)/period
	if (ratio.gt.zero) then
	   del=(twopi/period)**2/ratio
	   ec=hcgnl*del
C
C  Emax must be less that ec*critlr when horizontal symmetry is 
C  enforced
C
	   if (emax.gt.critlr*ec.and.ixsym.eq.1) then
	      write(*,20) char(7)
	      write(iout,20) ' '
	   end if
	end if
20	format(/1x,'************************************************'
     +	       /1x,'*EMAX may be too high.  Consider using XSYM=NO.*'
     +	       /1x,'************************************************'
     +	       /1x,a)
C 
C   calculate z and betaz. the field amplitude and phase error values are 
C   overwritten here.   
C 
	do i=1,nptst 
	   betaz(i)=sqrt(beta2-betax(i)*betax(i))
	end do 
C 
C   print out final error and initial and final electron velocities. 
C 
	perc=abs(x(nptst))/dflmin*1.d2 
	write(*,30) x(nptst),perc,betax(1),betax(nptst) 
30      format(/1x,'Trajectory calculations completed.' 
     +/1x,'Error after the final pass (cm) : ',1pe15.6, 
     +    ' (',0pf5.2,' % )' 
     +/1x,'Initial betax                   : ',1pe15.6 
     +/1x,'Final   betax                   : ',1pe15.6 )
C 
	if (ipass.eq.maxtry) then 
	   write(*,35) small*1.d2,nshot,char(7) 
35	  format(/1x,'WARNING : cannot reach goal error (',f5.2, 
     + ' %) in ',i2,' tries.'/1x,'Continuing calculations...',a1) 
	end if 
C 
C   finaly save results (if required) 
C 
	if (itraj.eq.2.or.mode.eq.0) then 
	   write(*,40) tfile(1:numchs(tfile)) 
40	   format(1x,'Saving trajectory in file ',a) 
	   open(imisc,file=tfile,form='unformatted',status='unknown') 
	   write(imisc) period,nper,npts,ener,efund 
	   write(imisc) (ct(i),x(i),betax(i),betaz(i),i=1,nptst) 
	   close(imisc) 
	end if 
C
C   z-range over which hanning windows will be applied.
C
	zhanl=dble(nhan)*period
	dhanl=zhanl
	zhanr=dble(nper-nhan)*period
	dhanr=dhanl
C
C   apply hanning windows, if necessary
C
	if (nhan.gt.0) then
	   do i=1,nptst
	      if (z(i).lt.zhanl) then
	         win=hanl(z(i),zhanl,dhanl)
	      else if (z(i).gt.zhanr) then
	         win=hanr(z(i),zhanr,dhanr)
	      end if
	      betax(i)=win*betax(i)
	      betaz(i)=win*betaz(i)
	   end do
	end if 
	return 
	end 
C 
C       ---------------------------- 
	subroutine derivs(zl,f,dfdz) 
C       ---------------------------- 
C 
C   calculate  the right-hand-side derivatives for the ODE integrator 
C   ct(z)=f(1); x(z)=f(2); betax(z)=f(3) 
C 
	implicit double precision (a-h,o-z) 
	include 'yaup.inc' 
	dimension f(*),dfdz(*)
C 
C   The equations of motion are 
C     dct/dz = 1/betaz = 1/sqrt(beta^2-betax^2) 
C     dx/dz = betax/betaz 
C     dbetax/dz =  - beta*xk0*K(z)/gamma = -delta*K(z) 
C   where K(z) is the K-factor corresponding to the local magnetic field. 
C   In the level of sophistication that we are dealing with, the field 
C   amplitude and the phase errors are constant over the small intervals 
C   that the solution is being advanced through. 
C 
	dctdz=one/sqrt(beta2-f(3)*f(3)) 
	dfdz(1)=dctdz 
	dfdz(2)=f(3)*dctdz 
C 
C   The local B-field and K-factor 
C 
	btot = bfield*sin(xk0*zl+phserr) 
	xkofz=0.934d0*btot*period
	dfdz(3)=-delta*xkofz 
	return 
	end 
C 
C       ------------------------------------------------------------- 
	subroutine shoot(ipass,betax0,dbetax,betax1,small,ier,errmsg, 
     +	           dz,dz0,dzmin,eps,maxstp) 
C       ------------------------------------------------------------- 
C 
C  Given an intial velocity betax0 and step dbetax, do one iteration 
C  of the shooting method to obtain an improved inital velocity betax1 and 
C  an estimate of the step for the next iteration.  SMALL is a 
C  parameter that is internally used to scale the sensitivity of transverse 
C  displacement at the exit of the undulator to the initial condition betax0. 
C  If this sensitivity is low the shooting method is not going to work. 
C  The rest of the parameters are passed on to the ODE integrator. 
C 
	implicit double precision (a-h,o-z) 
	parameter (nfunc=3) 
	dimension yini(nfunc) 
	character*(*) errmsg 
C 
C  Set up the initial conditions and integrate the equations of motion. 
C  The first "shot" is needed on the first iteration only.  If ipass.ne.0 
C  it will duplicate the last call to inival and is therefore omitted. 
C 
	if (ipass.eq.0) then 
	   ipass=ipass+1 
	   yini(1)=0.0d0 
	   yini(2)=0.0d0 
	   yini(3)=betax0 
	   call inival(ipass,yini,nfunc,dz,dz0,dzmin,eps,maxstp,ier) 
	   if (ier.ne.0) then 
	      ier=1 
	      write(errmsg,10) eps*1.d2 
10	      format('Unable to calculate the trajectory with ',f5.1, 
     +	             ' percent accuracy.') 
	      return 
	   end if 
	end if 
 
C 
C   the deviation of x(nptst) from 0 is a measure of the error in the initial 
C   condition. 
C 
	err1=yini(2)-0.0d0 
	write(*,15) err1 
15      format('+',t36,1pe15.6) 
C 
C  Second pass. 
C 
	betax1=betax0+dbetax 
	ipass=ipass+1 
	yini(1)=0.0d0 
	yini(2)=0.0d0 
	yini(3)=betax1 
	call inival(ipass,yini,nfunc,dz,dz0,dzmin,eps,maxstp,ier) 
	if (ier.ne.0) then 
	   ier=1 
	   write(errmsg,10) eps*1.d2 
	   return 
	end if 
C 
C  The error in this pass 
C 
	err2=yini(2)-0.0d0 
	write(*,15) err2 
C 
C   after two passes, calculate the "derivative" derr/dbetax0 of the error 
C   and use it to estimate an improved inital guess on betax0. If derr/dbetax0 
C   is less than small*err1 assume that the shooting method failed. 
C 
	derrdb=(err2-err1)/dbetax 
	if (abs(derrdb).lt.small*err1) then 
	   ier=1 
	   errmsg='Slope of xmax(betax0) is nearly zero. '// 
     +	          'Cannot continue.' 
 
	   return 
	end if 
C 
	dbetax=-err1/derrdb 
	betax1=betax0+dbetax 
C 
C   recalculate the trajectory with the "good" estimate. 
C 
	ipass=ipass+1 
	yini(1)=0.0d0 
	yini(2)=0.0d0 
	yini(3)=betax1 
	call inival(ipass,yini,nfunc,dz,dz0,dzmin,eps,maxstp,ier) 
	if (ier.ne.0) then 
	   ier=1 
	   write(errmsg,10) eps*1.d2 
	   return 
	end if 
	return 
	end 
C 
C       --------------------------------------------------------------- 
	subroutine inival(ipass,yini,nfunc,dzl,dz0,dzmin,eps,maxstp,ier) 
C       --------------------------------------------------------------- 
C 
C  Calculate the trajectory  with initial conditions yini. 
C  The solution is calculated at the points zi=(i-1)*dzl, i=1,nptst. 
C  Over each step the solution is advanced with a Burlich-Stoer ODE 
C  integrator with adaptive stepsize control in at most maxstp steps. 
C  EPS is a tolerance parameter and is rougly equal to the fractional 
C  error desired in the solution.  The magnetic field amplitude and phase 
C  errors over each interval [zi,zi+dzl] are assumed to be constant.  The 
C  field amplitude and phase errors are store in betaz and z, respectively. 
C 
	implicit double precision (a-h,o-z) 
	include 'yaup.inc' 
CRJD November 15, 2013: changed ndisp from 50 to 200
	parameter (ndisp=200) 
C 
	dimension yini(nfunc) 
	external derivs,bsstep 
C 
C   yini(1)=ct(z), yini(2)=x(z), yini(3)=betax(z) 
C 
	do i=2,nptst 
C 
C  Limits for the ODE integrator, field amplitude and phase error. 
C 
	   z1=dble(i-2)*dzl 
	   z2=z1+dzl 
	   bfield=betaz(i-1) 
	   phserr=z(i-1) 
C 
C   Advance from z1 to z2 
C 
	   call odeint(yini,nfunc,z1,z2,eps,dz0,dzmin,maxstp, 
     +	            nok,nbad,derivs,bsstep,ier) 
	   if (ier.ne.0) return 
C 
C   The trajectory at z2 is returned in yini 
C 
	   ct(i)=yini(1) 
	   x(i)=yini(2) 
	   betax(i)=yini(3) 
C 
C  Screen update 
C 
	   perc=dble(i)/dble(nptst)*1.d2 
	   if (mod(i,ndisp).eq.0) write(*,10) ipass,z2/1.d2,perc 
10	   format('+',t7,i2,t15,f6.3,t25,f6.1) 
	end do 
	return 
	end 
C 
C     --------------------------------------------------------- 
      subroutine odeint(ystart,nvar,x1,x2,eps,h1,hmin,maxstp, 
     +                  nok,nbad,derivs,bsstep,ier) 
C     --------------------------------------------------------- 
C 
C   bulirsch-stoer driver with adaptive stepsize control. 
C   input : 
C       ystart - initial values of the dependent varaibles 
C       nvar - number of dependent variables 
C       x1,x2 - integration range 
C       eps - accuracy 
C       h1 - inital guess on the step size 
C       hmin -  minimum allowed step size 
C       maxstp - maximum number of steps to take 
C       derivs - user supplied routine for calculation of the right- 
C                hand side derivative 
C       bsstep - stepper routine to be used 
C   output: 
C       ystart - values of dependent variables at x2 
C       nok, nbad - number of good and bad steps 
C       ier - error flag 
C             1 - h must be less than hmin to achieve required accuracy 
C             2 - cannot reach x2 with required accuracy in maxstp steps 
C 
C   copied from Numerical Recipes 
C 
      implicit double precision (a-h,o-z) 
      parameter (nvmax=10,zero=0.0d0,tiny=1.0d-10) 
      dimension ystart(nvar),yscal(nvmax),y(nvmax),dydx(nvmax)
CRJD November 15, 2013: added bsstep to external statement
      external derivs,bsstep 
C 
      ier=0 
      x=x1 
      h=sign(h1,x2-x1) 
      nok=0 
      nbad=0 
C 
C    save initial values 
C 
      do i=1,nvar 
        y(i)=ystart(i) 
      end do 
C 
C   take at most maxstp steps 
C 
      do nstp=1,maxstp 
        call derivs(x,y,dydx) 
C 
C   scaling used to monitor accuracy (constant fractional errors 
C   except near zero crossings) 
C 
        do i=1,nvar 
          yscal(i)=abs(y(i))+abs(h*dydx(i))+tiny 
        end do 
C 
C   rational function extrapolation to zero step size 
C   if step can overshoot end, cut down step size 
C 
        if((x+h-x2)*(x+h-x1).gt.zero) h=x2-x 
        call bsstep(y,dydx,nvar,x,h,eps,yscal,hdid,hnext,derivs) 
        if(hdid.eq.h)then 
          nok=nok+1 
        else 
          nbad=nbad+1 
        endif 
C 
C   are we done? 
C 
        if((x-x2)*(x2-x1).ge.zero)then 
          do i=1,nvar 
            ystart(i)=y(i) 
          end do 
          return 
        else if(abs(hnext).lt.hmin) then 
          ier=1 
          return 
        else 
          h=hnext 
        end if 
      end do 
      ier=2 
      return 
      end 
C 
C     --------------------------------------------------------------- 
      subroutine bsstep(y,dydx,nv,x,htry,eps,yscal,hdid,hnext,derivs) 
C     --------------------------------------------------------------- 
C 
C   burlisch-stoer step with monitoring of local truncation error 
C   to assure accuracy and adjust stepsize 
C 
C   input 
C       x - current value of the independent variable 
C       y - dependent variable at x 
C       dydx - first derivative of y at x 
C       nv - number of dependent variables 
C       htry - stepsize to be attempted 
C       eps - required accuracy 
C       yscal - scaling vector 
C       derivs - user-supplied routine to calculate right-hand- 
C                side derivatives 
C 
C   output 
C       hdid - accomplished stepsize 
C       hnext - estimated next step size 
C       x and y are replaced by their new values 
C 
C  copied from numerical recipes. 
C 
      implicit double precision (a-h,o-z) 
      parameter (nvmax=10,imax=11,nuse=7,one=1.0d0) 
      parameter (shrink=.95d0,grow=1.2d0) 
C 
      dimension y(nv),dydx(nv),yscal(nv),yerr(nvmax), 
     +    ysav(nvmax),dysav(nvmax),yseq(nvmax),nseq(imax)
      external derivs 
      save nseq 
      data nseq /2,4,6,8,12,16,24,32,48,64,96/ 
C 
      h=htry 
      xsav=x 
      do i=1,nv 
        ysav(i)=y(i) 
        dysav(i)=dydx(i) 
      end do 
C 
C  evaluate the sequence of midpoint integrations 
C  and extrapolate to zero stepsize. 
C 
1     do i=1,imax 
        call mmid(ysav,dysav,nv,xsav,h,nseq(i),yseq,derivs) 
        xest=(h/nseq(i))**2 
        call rzextr(i,xest,yseq,y,yerr,nv,nuse) 
        errmax=0.d0 
C 
C   find the maximum error and adjust the stepside  according to it. 
C 
        do j=1,nv 
          errmax=max(errmax,abs(yerr(j)/yscal(j))) 
        end do 
        errmax=errmax/eps 
        if(errmax.lt.one) then 
          x=x+h 
          hdid=h 
          if(i.eq.nuse)then 
            hnext=h*shrink 
          else if(i.eq.nuse-1)then 
            hnext=h*grow 
          else 
            hnext=(h*nseq(nuse-1))/nseq(i) 
          endif 
          return 
        endif 
      end do 
C 
      h=0.25d0*h/dble(2**((imax-nuse)/2)) 
      if(x+h.eq.x) pause 'Stepsize underflow in BSSTEP.' 
      goto 1 
      end 
C 
C     ------------------------------------------------------ 
      subroutine mmid(y,dydx,nvar,xs,htot,nstep,yout,derivs) 
C     ------------------------------------------------------ 
C 
C   modified midpoint step.  The dependent variable vector y of 
C   length nvar and its derivative dydx are input at xs. a total 
C   step htot in nstep substeps is to be made.  the output is 
C   returned in yout. 
C   copied from numerical recipes. 
C 
      implicit double precision (a-h,o-z) 
      parameter (nvmax=10) 
      dimension y(nvar),dydx(nvar),yout(nvar),ym(nvmax),yn(nvmax) 
C 
      h=htot/nstep 
      do i=1,nvar 
        ym(i)=y(i) 
        yn(i)=y(i)+h*dydx(i) 
      end do 
      x=xs+h 
      call derivs(x,yn,yout) 
      h2=2.0d0*h 
      do n=2,nstep 
        do i=1,nvar 
          swap=ym(i)+h2*yout(i) 
          ym(i)=yn(i) 
          yn(i)=swap 
        end do 
        x=x+h 
        call derivs(x,yn,yout) 
      end do 
      do i=1,nvar 
        yout(i)=0.5d0*(ym(i)+yn(i)+h*yout(i)) 
      end do 
      return 
      end 
C 
C     ----------------------------------------------- 
      subroutine rzextr(iest,xest,yest,yz,dy,nv,nuse) 
C     ----------------------------------------------- 
C 
C   use diagonal rational function extrapolation to evaluate nv 
C   functions at x=0 using a sequence of estimates with progressively 
C   smaller values x=xest and corresponding function vectors yest. 
C   copied from numerical recipes. 
C 
      parameter (imax=11,nmax=10,ncol=7) 
      implicit double precision (a-h,o-z) 
      dimension x(imax),yest(nv),yz(nv),dy(nv),d(nmax,ncol),fx(ncol) 
CRJD November 15, 2013: added save statement
      save d,x
      x(iest)=xest 
      if(iest.eq.1) then 
        do j=1,nv 
          yz(j)=yest(j) 
          d(j,1)=yest(j) 
          dy(j)=yest(j) 
        end do 
      else 
        m1=min(iest,nuse) 
        do k=1,m1-1 
          fx(k+1)=x(iest-k)/xest 
        end do 
        do j=1,nv 
          yy=yest(j) 
          v=d(j,1) 
          c=yy 
          d(j,1)=yy 
          do k=2,m1 
            b1=fx(k)*v 
            b=b1-c 
            if(b.ne.0.) then 
              b=(c-v)/b 
              ddy=c*b 
              c=b1*b 
            else 
              ddy=v 
            endif 
            v=d(j,k) 
            d(j,k)=ddy 
            yy=yy+ddy 
          end do 
          dy(j)=ddy 
          yz(j)=yy 
        end do 
      endif 
      return 
      end 
C
C***************************** fft.f *******************************
C
C       ------------------------------ 
	subroutine gauss(sig,da,ma,ra) 
C       ------------------------------ 
C 
C   Calculate gaussian weight factors for the convolution with 
C   the beam emmitance.  The factors are computed, stored in wrap- 
C   around order, and FFT'd once and for all. Called before convolution 
C   if nsig>0. 
C 
	implicit double precision (a-h,o-z) 
	include 'yaup.inc' 
	dimension ra(*) 
C 
C   Zero and positve angles first (remmember, wrap-around order). xnorm 
C   contains the normalizing factor AND the sampling stepsize. 
C 
	sig2=two*sig*sig 
	xnorm=da/(sqrt(twopi)*sig) 
	do i=1,ma+1 
	   a=dble(i-1)*da 
	   ra(i)=xnorm*exp(-a*a/sig2) 
	end do 
C 
C   Padd with zeros in the middle 
C 
	do i=ma+2,nconv-ma 
	   ra(i)=zero 
	end do 
C 
C   Then store negative angles using the symmetry 
C 
	do i=nconv-ma+1,nconv 
	   ra(i)=ra(nconv+2-i) 
	end do 
C 
C   FFT them 
C 
	no2=nconv/2 
	call realft(ra,no2,1) 
C 
C   The postive-frequency half of the FFT is now stored in ra. The purely 
C   real "frequency" 0 is in ra(1), the purely real "frequency" no2 is in 
C   ra(2), the real and imaginary parts of "frequency" 1 are in ra(2) and 
C   ra(3) respectively, etc.  See Numerical Recipes for more details. 
C 
	return 
	end
C 
C       ----------------- 
	subroutine fftpts 
C       ----------------- 
C 
C   Calculate the number of points for the fft optimized with respect to the 
C   current observation direction.   Since we need the brightness at a set 
C   of prespecified energies we will have to interpolate the FFT spectrum. 
C   This in turn means that the step in energy must be sufficiently small to 
C   describe  the rapid oscillations of the spectrum which are most severe 
C   for untapered undulator calculations.  This subroutine does that. It is 
C   called once for every observation direction.  This slows down the 
C   program insignificantly and there is always hope that the number of 
C   points in the FFT will be reduced. 
C 
	implicit double precision (a-h,o-z) 
	include 'yaup.inc' 
C 
C  The FFT spectrum is well known to be succeptible to power alaising. 
C  If Ec is the Nyquist energy of the FFT the power-aliasing pollution at 
C  each energy E will come from energies not lower than 2*Ec-E. 
C  Consider a untapered undulator with a fundamental energy Efund. If we 
C  require that at Emax the aliasing pollution come from energies that 
C  are Ncrit harmonics away then 
C	    2*Ec-Emax > Emax + Ncrit*Efund 
C  or 
C	    Ec = Emax + (Ncrit/2)*Efund 
C  The parameter Efund passed to this routine is the fundamental energy 
C  in the case of untapered undulators, or the "fundamental energy" 
C  corresponding to the weakest magnetic field (highest energy in the 
C  first emission band) in the case of tapered undulators. 
C 
	ecrit=emax+dble(ncrit)*(efund/two) 
C 
C  Limiting the Nyquist enrgy from below limits the sampling stepsize for 
C  the observer's time from above 
C	     dtau = hc/(2*Ec) 
C  and since dtau=tautot*(nfft-1) we have 
C	     nfft = hc / (2*Ec*tautot) + 1 
C  which must also be a power of 2.  Then 
C 
	dtau=hc/(two*ecrit) 
	tautot=ct(nptst)-zn*dble(nper)*period 
	nfft=nint(tautot/dtau)+1 
	nfft=2**(int(log(dble(nfft))/log(two))+1) 
	if (nfft.eq.1) nfft=2 
	dtau=tautot/dble(nfft-1) 
	ecrit=hc/(two*dtau) 
	defft=hc/(dble(nfft)*dtau) 
C 
C  Make sure that the FFT arrays will not overflow.  Note that since 
C  beta is "bandwidth limited" (i.e. non-zero) to the interval 
C  [0,tautot], the smallest stepsize in energy that would make sense 
C  is hc/tautot (Nyquist sampling theorem). 
C 
	if (nfft.gt.maxfft) then 
	   if (icont.eq.0.or.icont.eq.2) 
     +		call bclose('delete',idmp1,idmp2,ista)
	   nfft=maxfft 
	   dtau=tautot/dble(nfft-1) 
	   defft=hc/(dble(nfft)*dtau) 
	   ecrit1=hc/(two*dtau) 
	   crpar=two*(ecrit1-emax)/efund 
	   if (crpar.gt.one) then 
	      write(errmsg,10) int(crpar)
10	      format('use NCRIT =',i3,
     +		' if you REALLY want this calculation.')
     	   else
	      errmsg='cannot complete calculation for any NCRIT.'
	   end if
	   call leave('fftpst::  FFT array overflow (NCRIT).',
     +		'EMAX/XPS/YPS are probably too big.', errmsg)
	end if 
C 
C   The zero-emittance spectrum of untapered undulators  oscillates rapidly 
C   due to the prefactor 
C		( sin(N*w*x/2c)/sin(w*x/2c) )**2 
C   where N is the number of periods, w is the photon frequency, and 
C		 x=c*T-nz*period 
C   Here T is the time it takes the electron (in a untapered undulator) to 
C   travel a distance period along the undulator axis.  The period (in 
C   frequency) of the fastest oscillations is therefore 
C		 dw=pi*(2*c)/(N*x), 
C		 dE=hbar*dw=2*pi*hbar*c/N*x=hc/N*x 
C   where hc is Planck's constant times the speed of light (12.399d-5 eVcm). 
C   Since 
C		N*x=c*N*T-nz*N*period=C*Ttot-nz*L=tautot 
C   (Ttot is the total transit time through the undulator and L is the 
C   total undulator length) where tautot is the total transit time for the 
C   electron from the observer's point of view (observer "time"), 
C   the  period of the untapered-undulator-spectrum oscillations will be 
C		dE=hc/tautot 
C   This, however is just defft - the FFT stepsize up to this point. 
C   The FFT spectrum will be interpolated over the energy mesh 
C   requested by the user (in SAVDMP).  In order to ensure that the 
C   interpolating polynomial will faithfully reproduce these oscillations, 
C   the stepsize for the FFT must be reduced (by adding padding zeros at 
C   the end of the array to be FFT'ed). The number of points per oscillation 
C   that will do the job will depend on the interpolation scheme, but eight 
C   should be about right (if you do not think eight is enough increase the 
C   value of mfftd in yaup.inc, but keep in mind that the time 
C   for the FFT increases like mfft1*(1+log2(mfft1)/log2(nfft)).  PLEASE 
C   make sure that the new value of mfft1 is a power of 2). 
C 
C   Generally speaking this is a redundancy because the spectrum is 
C   completely determined by the current sampling through the Nyquist 
C   sampling theorem.  Any attempt to apply that theorem directly, however, 
C   will fail because its performance depends on the sensitive cancelation 
C   of rapidly oscillating terms. 
C 
	mfft=mfft1 
	if (mfft*nfft.gt.maxfft) then 
	   if (icont.eq.0.or.icont.eq.2) 
     +		call bclose('delete',idmp1,idmp2,ista)
	   mfft=maxfft/nfft
	   if (mfft.gt.0) then
	      write(errmsg,20) mfft 
20	      format('use MFFT =',i3,
     +		' if you REALLY want this calculation.') 
     	   else
	      errmsg='cannot complete calculation for any MFFT.'
	   end if
	   call leave('fftpst::  FFT array overflow (MFFT).',
     +		'EMAX/XPS/YPS are probably too big.', errmsg)
	end if 
	defft=defft/dble(mfft) 
C 
C  Only a complete dummy will do calculations in this regime 
C  but you never know what a user will come up with. 
C 
	if (emax.gt.ecrit-defft) then
	   if (icont.eq.0.or.icont.eq.2) 
     +		call bclose('delete',idmp1,idmp2,ista)
	   write(errmsg,30) ecrit-defft
30	   format('use EMAX =',f8.1,
     +		' if you **REALLY REALLY** want this calculation.') 
	   call leave('fftpst::',errmsg,
     +		'Or better yet, get a life!!!')
	end if 
C 
	return 
	end 
C 
C     ------------------------------- 
      subroutine four1(data,nn,isign) 
C     ------------------------------- 
C 
C   PURPOSE    Calculates the FFT of a data set of NN complex data points 
C 
C   I/O        DATA   -- the data points.  In the output it is replaced 
C                        by the FFT in wrap-around order. 
C              NN     -- number of data points.  MUST be a power of 2. 
C              ISIGN  --  1, calculate the direct transform 
C                        -1, calculate the the inverse transform. In this 
C                         case DATA must additionaly be divided by NN. 
C 
C   SOURCE     Numerical Recipes, pp. 394-395. 
C 
      implicit double precision (a-h,o-z) 
      dimension data(2*nn) 
C 
C   Note - the original statement was dimension data(*), but for the 
C   sake of sanity (and safety) it was replaced. 
C 
      n=2*nn 
      j=1 
      do i=1,n,2 
        if(j.gt.i)then 
          tempr=data(j) 
          tempi=data(j+1) 
          data(j)=data(i) 
          data(j+1)=data(i+1) 
          data(i)=tempr 
          data(i+1)=tempi 
        endif 
        m=n/2 
        do while ((m.ge.2).and.(j.gt.m)) 
          j=j-m 
          m=m/2 
        end do 
        j=j+m 
      end do 
      mmax=2 
      do while (n.gt.mmax) 
        istep=2*mmax 
        theta=6.28318530717959d0/(isign*mmax) 
        wpr=-2.0d0*sin(0.5d0*theta)**2 
        wpi=sin(theta) 
        wr=1.0d0 
        wi=0.0d0 
        do m=1,mmax,2 
          do i=m,n,istep 
            j=i+mmax 
C            tempr=real(wr)*data(j)-real(wi)*data(j+1) 
C            tempi=real(wr)*data(j+1)+real(wi)*data(j) 
            tempr=wr*data(j)-wi*data(j+1) 
            tempi=wr*data(j+1)+wi*data(j) 
            data(j)=data(i)-tempr 
            data(j+1)=data(i+1)-tempi 
            data(i)=data(i)+tempr 
            data(i+1)=data(i+1)+tempi 
          end do 
          wtemp=wr 
          wr=wr*wpr-wi*wpi+wr 
          wi=wi*wpr+wtemp*wpi+wi 
        end do 
        mmax=istep 
      end do 
      return 
      end 
C 
C 
C     ------------------------------- 
      subroutine realft(data,n,isign) 
C     ------------------------------- 
C 
C  PURPOSE        Calculates the FFT of a set of 2N real-valued data 
C                 points by the positive frequency half of its transform. 
C                 The real-valued first anf last points of the FT are 
C                 returned as  the elements DATA(1) and DATA(2). Also 
C                 calculates the inverse transform of a complex data 
C                 array  if it is the transform of a real function.  In this 
C                 case the result must be multiplied by 1/N. 
C 
C  IO             DATA  -- 2N real-valued points or N complex points. 
C                 N     -- Number of points.  MUST be a power of 2. 
C                 ISIGN --  1, calculate the direct FFT 
C                          -1, callculate the inverse FFT 
C 
C  SOURCE         Numerical Recipes, p.400 
C 
      implicit double precision (a-h,o-z) 
      dimension data(2*n) 
C 
C   Again, the original statement was dimension data(*), and at certain 
C   lines the original routine tried to access elements 2*n+1 and 2*n+2, 
C   thus unacceptably overwriting data if the length of DATA in the 
C   calling routine is EXACTLY 2*n.  The routine was modified a bit 
C   to avoid this. 
C 
      theta=3.141592653589793d0/dble(n) 
      c1=0.5d0 
      if (isign.eq.1) then 
        c2=-0.5d0 
        call four1(data,n,+1) 
      else 
        c2=0.5d0 
        theta=-theta 
      endif 
      wpr=-2.0d0*sin(0.5d0*theta)**2 
      wpi=sin(theta) 
      wr=1.0d0+wpr 
      wi=wpi 
      n2p3=2*n+3 
C 
C   Case i=1 done separately below. 
C 
      do i=2,n/2 
        i1=2*i-1 
        i2=i1+1 
        i3=n2p3-i2 
        i4=i3+1 
C 
C   In the original routine.  Since all variables are real*8 here, 
C   these two lines are no longer needed. 
C        wrs=real(wr) 
C        wis=real(wi) 
C 
        wrs=wr 
        wis=wi 
        h1r=c1*(data(i1)+data(i3)) 
        h1i=c1*(data(i2)-data(i4)) 
        h2r=-c2*(data(i2)+data(i4)) 
        h2i=c2*(data(i1)-data(i3)) 
        data(i1)=h1r+wrs*h2r-wis*h2i 
        data(i2)=h1i+wrs*h2i+wis*h2r 
        data(i3)=h1r-wrs*h2r+wis*h2i 
        data(i4)=-h1i+wrs*h2i+wis*h2r 
        wtemp=wr 
        wr=wr*wpr-wi*wpi+wr 
        wi=wi*wpr+wtemp*wpi+wi 
      end do 
      if (isign.eq.1) then 
        h1r=data(1) 
        data(1)=h1r+data(2) 
        data(2)=h1r-data(2) 
      else 
        h1r=data(1) 
        data(1)=c1*(h1r+data(2)) 
        data(2)=c1*(h1r-data(2)) 
        call four1(data,n,-1) 
      endif 
      return 
      end 
C 
C       ----------------- 
	subroutine conv2d
C       ----------------- 
C 
C   Calculate the 2D convolution of the arrays s1 and s2 with the beam 
C   emiitances rx and ry.  The 2D convolution is calculated as not more than 
C   nx+ny 1D convolutions. Even though after the convolution 
C   only rows (mx+1,nx-mx) and columns (my+1,ny-my)contain useful information, 
C   a little thought shows that the convolution must be performed over all 
C   rows(columns), followed by a convolution ONLY over the columns(rows) that 
C   will contain significant data.  On input, s1 contains the zero- 
C   emmitance Stoke's coefficient S0 at constant x in rows, and constant y in 
C   columns,  i.e. s1=s1(x,y), and s2 contains the zero-emittance Stoke's 
C   coefficient S1 in the same format. 
C 
	implicit double precision (a-h,o-z) 
	include 'yaup.inc' 
C 
	dimension s1(nconv,nconv),s2(nconv,nconv) 
	equivalence (s1,fx), (s2,fz) 
C 
C   Zero-padd signal (if needed) 
C 
	call endpad 
C 
C   Calculate number of convolutions for row/columns and columns/rows 
C   processing. 
C 
	nc1=nx+ny-2*my 
	nc2=ny+nx-2*mx 
C 
C   If nc1<nc2 do rows first (x), and then significant columns (y), 
C   otherwise columns first (y), and then rows (x). 
C 
	if (nc1.le.nc2) then 
	   call transp(max(nx,ny)) 
	   call do2d(nx,ry,ny,my,rx) 
	else 
	   call do2d(ny,rx,nx,mx,ry) 
	   call transp(max(nx,ny)) 
	end if 
C 
C   "Delete" the last mx/my and rows/columns.  The first mx and my 
C   columns/rows are useless too, but they will be handled by EXTRCT. 
C 
	if (icx.eq.1) nx=nx-mx 
	if (icy.eq.1) ny=ny-my 
	return 
	end 
C 
C       ----------------- 
	subroutine endpad 
C       ----------------- 
C 
C  This routine padds the end of the signal arrays to be convoluted 
C  using the reflection symmetry of the spectrum.  Called before 
C  convoluting, it ensures that sufficient signal s is present to 
C  do the convolution correctly, and enough padding is available to 
C  prevent wrap-around "pollution" of the convoluted signal at the ends. 
C  The rows of the arrays are assumed to contain the signal at constant y, 
C  and the columns - at constant x, (e.g. s=s(x,y)). 
C 
	implicit double precision (a-h,o-z) 
	include 'yaup.inc' 
C 
	dimension s1(nconv,nconv),s2(nconv,nconv) 
	equivalence (s1,fx), (s2,fz) 
C
C  If the angular range has been optimized (using the reflection symmetries
C  of the spectrum), "reflect" the near-zero portion for the purposes of
C  the convolution
C 
	if (icx.eq.1) then 
	   do iy=1,ny 
	      do ix=1,mx 
	         s1(nx+ix,iy)=s1(nx-ix,iy) 
	         s2(nx+ix,iy)=s2(nx-ix,iy) 
	      end do 
	   end do 
	   nx=nx+mx 
	end if 
	if (icy.eq.1) then 
	   do ix=1,nx 
	      do iy=1,my 
	         s1(ix,ny+iy)=s1(ix,ny-iy) 
	         s2(ix,ny+iy)=s2(ix,ny-iy) 
	      end do 
	   end do 
	   ny=ny+my 
	end if 
C
C  Zero-padd the arrays to be convoluted.
C 
	do ix=nx+1,nconv 
	   do iy=1,nconv 
	      s1(ix,iy)=zero 
	      s2(ix,iy)=zero 
	   end do 
	end do 
	do iy=ny+1,nconv 
	   do ix=1,nx 
	      s1(ix,iy)=zero 
	      s2(ix,iy)=zero 
	   end do 
	end do 
	return 
	end 
C 
C       ------------------------------- 
	subroutine do2d(n1,r2,n2,m2,r1) 
C       ------------------------------- 
C 
	implicit double precision (a-h,o-z) 
	include 'yaup.inc' 
	dimension r1(nconv),r2(nconv) 
C 
	dimension s1(nconv,nconv),s2(nconv,nconv) 
	equivalence (s1,fx), (s2,fz) 
C 
	do i1=1,n1 
	   call convlv(s1(1,i1),r2,nconv) 
	   call convlv(s2(1,i1),r2,nconv) 
	end do 
	call transp(max(nx,ny)) 
	do i2=m2+1,n2-m2 
	   call convlv(s1(1,i2),r1,nconv) 
	   call convlv(s2(1,i2),r1,nconv) 
	end do 
	return 
	end 
C
C     -------------------------------- 
      subroutine convlv(data,respns,n) 
C     -------------------------------- 
C 
C  Calculate the convolution of the "signal" DATA, and the "response" 
C  RESPNS using FFTs.  This discrete convolution is apprx. equal to the 
C  continuous convolution times the sampling stepsize.  The convolution 
C  is stored back in DATA, so the original contents of DATA are lost. 
C  RESPNS must be properly precalculated, zero-padded (see Numerical 
C  Recipes), and FFTd.  DATA must also be properly zero padded to avoid 
C  "polution" of the end chanels.  Adapted from Numerical Recipes. 
C 
      implicit double precision (a-h,o-z) 
      complex*16 data(n/2),respns(n/2) 
C 
C  FFT signal function. Response is already supposed to be FFTd. 
C 
      no2=n/2 
      call realft(data,no2,1) 
C 
C  Calculate the convolution in the "frequency" domain. Since it will be 
C  IFFTd, and the inverse is a real function, only the postive frequency 
C  half needs to be calculated.  When performing the multiplication, take 
C  into account the specific packing of the data.  The factor no2 is 
C  not part of the convoltion, but of the IFFT. 
C 
      do i=2,no2 
         data(i)=data(i)*respns(i)/no2 
      end do 
      data(1)=cmplx(Dreal(data(1))*Dreal(respns(1)), 
     +	      imag(data(1))*imag(respns(1)))/dble(no2) 
C 
C  IFFT to find the convolution in the "time" domain.  The data is already 
C  prepacked as REALFT requires 
C 
      call realft(data,no2,-1) 
      return 
      end 
C
C***************************** spec.f **********************************
C 
C       ------------------------------- 
        subroutine specal(irec,ix0,iy0) 
C       -------------------------------
C 
C  Calculate spectrum of id radiation in the user-specified mode.
C  ix0 and iy0 determine the lower boundaries of the angle-mesh 
C  loops, and irec points to the last used record in the dump files
C 
        implicit double precision (a-h,o-z) 
        include 'yaup.inc'   
C 
C   Zero emmitance spectrum.  Print some info  to screen to keep 
C   the user awake. The dimensions of xmmin/max and ymin/max depend
C   on the values of DIST.
C 
        xmin=x0*dist 
        ymin=y0*dist
        xmax=xmin+(nx-1)*dx*dist 
        ymax=ymin+(ny-1)*dy*dist 
        write(*,10) emin,emax,ne,1.d3*xmin,1.d3*xmax,nx,
     +	                         1.d3*ymin,1.d3*ymax,ny 
10      format(/1x,60('-')/ 
     +/1x,'Beginning zero-emittance brightness calculations.' 
     +/1x,'emin = ',f8.1,3x,'emax = ',f8.1,3x,'ne =',i4, 
     +/1x,'xmin = ',f8.4,3x,'xmax = ',f8.4,3x,'nx =',i4, 
     +/1x,'ymin = ',f8.4,3x,'ymax = ',f8.4,3x,'ny =',i4) 
        write(*,15) 
15      format(/4x,'x (mrad/mm)',t18,'y (mrad/mm)',t31,'% done'/) 
C 
C   Actual (user-requested) pinhole limits (in radians) 
C 
        xmin=xpc-xps/two 
        xmax=xpc+xps/two 
        ymin=ypc-yps/two 
        ymax=ypc+yps/two
C
C  IDMP counts the number of times the dump files are written.  Its use
C  is completely limited to the module bwrite, but keep it external to 
C  avoid DATA IDMP/0/ ambiguities.
C  
        ntot=nx*ny 
        ndone=ix0*ny
	idmp=0 
C
C   ix0 and iy0 are the last executed (and saved) indices according
C   to the status file. 
C
        do ix=ix0+1,nx 
           xc=x0+dble(ix-1)*dx 
           do iy=iy0+1,ny 
              yc=y0+(iy-1)*dy 
              ndone=ndone+1 
              call intgrd(xc,yc) 
              call bright(xc,yc,e0,ne0)
              if (mod(ndone,nprint).eq.0.or.ndone.eq.1.or.ndone.eq.ntot) 
     +           call report(dist,xc,yc,ndone,ntot) 
	      call bwrite(idmp,irec,ix,iy,e0,ne0) 
           end do 
        end do
C
C  the zero-emittance spectrum has been fully-calculated at this point
C  proceed with convolution and data storage.  If mode=7 (integrated power) 
C  set ne=1 
C 
	if (mode.eq.7) ne=1
        if (nsig.ne.0) write(*,20) 
20      format(/1x,60('-')//1x,'Convoluting with the beam divergence.' 
     +         //1x,'energy (eV)',t25,'% done'/) 
C 
C   Read data from the dump files back into two arrays.  They will contain 
C   the zero-emittance Stoke's parameters S0 and S1. The file slows things 
C   down but avoids  extensive memory swapping and severe limitations on 
C   the number of  points in energy and angle 
C 
        call bopen 
        if (ne.ne.1) de=(emax-emin)/dble(ne-1) 
        do ie=1,ne 
           ec=emin+dble(ie-1)*de 
           call bread(ie) 
C 
C   Perform the 2D convolution (if necessary).  The degree of linear 
C   polarazation by definition is S1/S0. 
C 
           if (nsig.ne.0) then 
              if (mod(ie,nprint).eq.0.or.ie.eq.1.or.ie.eq.ne) then 
                 write(*,30) ec,(ie*100)/ne 
30               format('+',1x,f9.0,t26,i3) 
              end if 
              call conv2d 
           end if 
C 
C   degree of polarization - brightness calculations only. 
C 
	   if (mode.ge.1.and.mode.le.3) call dgpol 
C 
C   Restore the data from the padded (and possibly optimized) range to 
C   the user-requested range. First restore the x-data (it is in the columns) 
C 
           call extrct(ny,nx,mx,dx,x0,xmin,xmax,icx,nxp1)
C 
C   transpose both arrays (extrct works on columns only) then  restore the 
C   y-data, and then transpose once again to return things to their original 
C   state. 
C 
           call transp(max(nxp1,ny)) 
           call extrct(nxp1,ny,my,dy,y0,ymin,ymax,icy,nyp1) 
           call transp(max(nxp1,nyp1))
C
C   another internal check
C
	   if (nxp.ne.nxp1.and.nyp.ne.nyp1) then
	      call bclose('keep',idmp1,idmp2,ista)
	      call leave('specal:: internal error',' ',' ')
     	   end if
C 
C  Integrate over pinhole (if required) 
C 
           if (mode.eq.4) 
     +        call int2d(xmin,xmax,ymin,ymax,flux) 
C 
C   Save output 
C 
           call putdat(flux,ec,xmin,ymin) 
        end do 
C
C   If ICONT=0,2 delete the dump and status files, else keep them
C
	if (icont.eq.0.or.icont.eq.2) then 
           call bclose('delete',idmp1,idmp2,ista)
	else if (icont.eq.1.or.icont.eq.3) then
	   call bclose('keep',idmp1,idmp2,ista)
	end if
        return
        end
C
C       ------------------------ 
        subroutine intgrd(xc,yc) 
C       ------------------------ 
C 
C  calculate the integrands to be FFT'd.  the beta's are splined 
C  and evaluated over an equidistant mesh in tau. Hanning windows
C  are applied here.
C 
        implicit double precision (a-h,o-z) 
        include 'yaup.inc' 
C 
C   WARNING : the equivalence statement below is made possible by the fact 
C   that the array WA is a temporary work array used by the splining routine 
C   for storage purposes only, and once SPLINE is exited is no longer 
C   needed.  The same is NOT true for tau,bx2, and bz2, which are needed 
C   all the way up to the end of this subroutine.  Do NOT make them 
C   equivalent to fx and/or fz, since the results of this subroutine are 
C   stored in fx and fz.  it is assumed that maxpts<maxfft. 
C 
        dimension tau(maxpts),wa(maxpts) 
        dimension bx2(maxpts),bz2(maxpts) 
        equivalence (wa,fz) 
C 
C   first, the wave observation direction 
C 
        zn=one/sqrt(one+xc*xc+yc*yc) 
        xn=zn*xc 
        yn=zn*yc 
C 
C   then the observer's "time" 
C 
        do i=1,nptst 
           rdotn=xn*x(i)+zn*z(i) 
           tau(i)=ct(i)-rdotn 
        end do 
C 
C   calculate the number of points and stepsizes for the fft 
C 
        call fftpts 
C 
C   spline the trajectory over the observer's time and evaluate over an 
C   equidistant mesh. 
C 
        call spline(tau,betax,nptst,bx2,wa) 
        call spline(tau,betaz,nptst,bz2,wa) 
C 
C   the knots are not equidistant!!! 
C 
        dtau=tau(nptst)/dble(nfft-1) 
        h=tau(2)-tau(1) 
        ntau1=1 
        ntau2=2 
        tau2=tau(2) 
        do i=1,nfft 
           tc=dble(i-1)*dtau 
C
C  on some systems tc can exceed tau(nptst) by as little as 1e-24
C  due to roundoff (12/1/93 bib)
C
           do while (tc.gt.tau2.and.ntau2.lt.nptst) 
              ntau1=ntau1+1 
              ntau2=ntau2+1 
              h=tau(ntau2)-tau2 
              tau2=tau(ntau2) 
           end do 
           a=(tau(ntau2)-tc)/h 
           b=(tc-tau(ntau1))/h 
           bx=a*betax(ntau1)+b*betax(ntau2)+ 
     +        ((a**3-a)*bx2(ntau1)+(b**3-b)*bx2(ntau2))*(h**2)/6.0d0 
           bz=a*betaz(ntau1)+b*betaz(ntau2)+ 
     +        ((a**3-a)*bz2(ntau1)+(b**3-b)*bz2(ntau2))*(h**2)/6.0d0 
           bdotn=1-bx*xn-bz*zn 
           fx(i)=bx*(dtau/bdotn) 
           fz(i)=bz*(dtau/bdotn) 
        end do 
C 
C   padd with zeros 
C 
        do i=nfft+1,nfft*mfft 
           fx(i)=zero 
           fz(i)=zero 
        end do 
C 
        return 
        end 
C 
C       ------------------------------- 
        subroutine bright(xc,yc,e0,ne0) 
C       -------------------------------
C 
C  calculate the zero-emittance brightness in ph/s/%BW/mrad^2 
C 
        implicit double precision (a-h,o-z) 
        include 'yaup.inc' 
C 
C   some physical constants (Kittel) 
C 
        parameter ( e = 4.80325d-10 ) 
        parameter ( c = 2.997925d10 ) 
        parameter ( fac = e*e/(twopi*twopi*c) ) 
C 
C   When the square of the FFT (to be calculated below) is multiplied 
C   by fac*(E/hbar*c)**2, the result is the energy emitted by one electron 
C   per unit frequency interval, per solid angle in units of erg/s^-1/rad^2. 
C   The extra 1/c^2 is introduced because the "integration" here is with 
C   respect to ct, not t alone.  To convert this to ph/s/mrad^2/0.1%BW : 
C       - multiply by 1.d-6 (erg/s^-1/rad^2 to erg/s^-1/mrad^2) 
C       - divide by the photon energy (in ergs) to convert this to 
C                  ph/s^-1/mrad^2. 
C       - multiply by the photon frequency to convert to ph/mrad^2/100%BW. 
C                  The last two steps are equivalent to dividing by 
C                  hbar(erg.s) 
C       - multiply by the bandwidth to get the final result. 
C       - to account for all electrons, multiply by cur/1.6e-19 
C   All this amounts to 1d-6*bw/hbar*(cur/e) where everything is in MKS units. 
C 
        parameter ( eMKS = 1.60219d-19 ) 
        parameter ( hbar = 1.05459d-27, hbarc = 12.399d-5/twopi ) 
C 
        dimension brix(maxfft/3),briy(maxfft/3) 
        complex*16 ffx(maxfft/2),ffz(maxfft/2),fdotn 
        equivalence (ffx,fx), (ffz,fz) 
        equivalence (fx(1),brix(1)), (fz(1),briy(1)) 
C 
C   the radiation integral is n x ( n x FT(beta) ).  Here only FT(beta) is 
C   calculated. The first and last purely real last components should be 
C   unpacked, but they will always be way out of the range of interest, so 
C   forget about them. 
C 
        call realft(ffx,nfft*mfft/2,1) 
        call realft(ffz,nfft*mfft/2,1) 
C 
C  Keep the data only in the energy range of interest. 
C 
        i1=int(emin/defft)+1 
        i2=int(emax/defft)+2 
C
C  Keep a couple of extra points for the spline interpolation 
C  (power load calculations)
C
	if (i1.gt.10) i1=i1-5
	if (i2.lt.nfft*mfft-10) i2=i2+5      
	e0=dble(i1-1)*defft 
        ne0=i2-i1+1 
C
C  This error handler should never be needed, but you never know.
C  At most maxfft/3 points may be retained.
C 
        if (i1.lt.2.or.i2.gt.nfft*mfft-1.or.ne0.gt.maxfft/3)  then
           write(errmsg,10) 1.d3*xc*dist,1.d3*yc*dist 
10         format('bright:: FFT array overflow  at ',f8.4,',',
     +		f8.4,' mrad/mm.') 
           call leave(errmsg,' ',' ') 
        end if 
C
C  BAC-CAB rule and magnitude squared. The radiation along z is 
C  negligible (or is it?).  The prefactor below converts to  to 
C  erg/s^-1/rad^2.  The extra factor c^2 appears because the FFT 
C  is with respect to c*tau, not tau alone.   
C
	prefac=fac/(hbarc**2)
C
C  Brightness or power?
C
	if (mode.ge.1.and.mode.le.4) then
C 
C   Convert to ph/s/mrad^2/%BW. cur/eMKS is the number of electrons 
C   contributing to the emission. 
C 
	   prefac=prefac*(cur/eMKS)*(1.0d-6*bw/hbar)
	else if (mode.ge.5.and.mode.le.7) then
C
C  or switch to watts/eV/mrad^2.  To convert, divide by hbar ( in
C  eV*s, erg/s^-1/rad^2 -> erg/eV/rad^2), multiply by 1.0d-6 
C  (rad^2 -> mrad^2), 10^-7 (ergs - > Joules), and finally by cur/eMKS
C  (Joules - > Watts).  The hbar defined above is in erg*s.  If we multiply
C  it by 10^-7/eMKS it will be in eV*s, i.e. the total conversion factor 
C  will be
C
C	prefac = (cur/eMKS)*(1.0d-6)*(1.0d-7)/(10^-7*hbar/eMKS)
C              = 1.0d-6*cur/hbar
C
	   prefac=prefac*1.0d-6*cur/hbar
	end if
C
C   Convert.  Divide by dist to convert everything to spatial units,
C   if necessary.
C 
        do i=i1,i2 
           j=i-i1+1 
           ephot=dble(i-1)*defft 
           fdotn=xn*ffx(i)+zn*ffz(i) 
           brix(j)=prefac*(ephot*abs(xn*fdotn-ffx(i))/dist)**2 
           briy(j)=prefac*(ephot*abs(yn*fdotn-zero  )/dist)**2 
C          briz   =prefac*(ephot*abs(zn*fdotn-ffz(i))/dist)**2 
	end do
C 
        return 
        end 
C 
C       ---------------- 
        subroutine dgpol 
C       ---------------- 
C 
C   Calculate the degree of polarization.   The degree of polarization 
C   is simply the ratio of the Stoke's coefficients S1 and S0. 
C 
        implicit double precision (a-h,o-z) 
        include 'yaup.inc' 
C 
        dimension s0(nconv,nconv),s1(nconv,nconv),degp(nconv,nconv) 
        equivalence (fx,s0), (fz,s1,degp) 
C 
        do ix=1,nx 
           do iy=1,ny 
              degp(ix,iy)=s1(ix,iy)/s0(ix,iy) 
           end do 
        end do 
        return 
        end 
C
C***************************** misc.f *****************************
C
C  WARNING : on most computers the following BLOCKDATA will cause
C  an increrase in the size of the executable code by at least 2Mb.
C  This seems to be related to the way BLOCKDATA initializations
C  are handled. Since there is no real need for this BLOCKDATA  I've
C  left it here only for the neat freaks who are willing to sacrifice
C  performance for principle.  If you really want all common blocks
C  massaged to avoid potential problems, use explicit assignments,
C  e.g. period=zero, not DATA statements.  The explicit assignment 
C  will not increase the size of the executable code.
C
C	---------------
C	blockdata init
C	---------------
C
C   Initialize all common blocks
C
C	implicit double precision (a-h,o-z)
C	include 'yaup.inc'
C
C	data period,emin,emax,ener,cur,
C     +	     sigx,sigy,sigx1,sigy1,dist,xpc,ypc,xps,yps /14*zero/, 
C     +       nper,npts,ne,nxp,nyp,mode,nsig,itraj,nhan /9*0/
C 
C	data bw,iexp,icont,isav,nres,ncrit,mfft1 /zero, 6*0/
C
C	data x0,dx,y0,dy,nx,mx,icx,ny,my,icy /4*zero,6*0/ 
C
C	data fx,fz /maxfft*zero,maxfft*zero/,
C     +	     efund,defft,xn,yn,zn,nfft,mfft /5*zero,2*0/
C
C	data wk0,delta,beta2,xlam,bfield,
C     +	     phserr,zhanl,dhanl,zhanr,dhanr /10*zero/
C
C	end
C 
C	-----------------
	subroutine chkpar
C	-----------------
C
C   Make sure the code is not messed up too badly.  This is
C   extremely rudimentary.
C
	implicit double precision (a-h,o-z) 
	include 'yaup.inc'
C
C   WARNING : maxfft must be a power of 2 and be equal to or greater than 
C   the larger of the numbers nconvol**2 and maxpts.  maxpts 
C   must be bigger than nconv. the algorithm of the program depends on these 
C   values, so PLEASE do not muck around with them unless you REALLY know 
C   what you are doing. 
C
 	nbig=max(nconv**2,maxpts)
	if (maxfft.lt.nbig.or.maxpts.lt.nconv) call leave
     +		('This source code has been modified.',' ',
     +		 'Please get a current copy.' )
C
       	ifile='yaup.inp'
	return
	end
C 
C       ----------------- 
	subroutine bounds 
C       ----------------- 
C
C   Calculate angular limits for zero-emittance calculations
C   The spectrum is calculated at NA points with step DA, starting
C   at A0 (A=X,Y).  MA is the number of points outside the
C   user-requested range that are needed only for the convolution with
C   the zero emittance.
C
	implicit double precision (a-h,o-z) 
	include 'yaup.inc' 
C 
C   If ixsym was not set to zero in the input file, then left-right 
C   symmetry of the spectrum is assumed and the angular range is 
C   minimized. 
C 
	iysym=1 
	call choose(xpc,xps,nxp,sigu,x0,dx,nx,mx,icx,ixsym) 
	call choose(ypc,yps,nyp,sigv,y0,dy,ny,my,icy,iysym)
C
C   ISAV controlls how often will the dump and status files be updated.
C   If ISAV=0, they are never updated (except before program termination)
C
	if (isav.eq.0) isav=nx 
	return 
	end
C
C       ----------------------------------------------------------- 
	subroutine choose(acent,asiz,npa,sig,a0,da,na,ma,ica,iasym) 
C       ----------------------------------------------------------- 
C 
C   Given the center ACENT and size ASIZ of a pinhole in a given direction 
C   (x or y), the desired number of points NPA, beam divergence SIG, and the 
C   number of standard deviations NSIG to take into account, calculate 
C   starting value A0, step DA, actual the number of points NA, and the number 
C   of points MA at which to sample the beam divergence for that direction. 
C   Due to limitations in MA, and the fact that the sampling stepsize for 
C   the spectrum and the beam divergence must be the same, the value of NPA 
C   (desired number of points) may be adjusted.  The routine tries to deter- 
C   mine the minimum range over which to calculate the spectrum, and after 
C   the calculation corrections and adjustments may be needed. In this case 
C   ICA is set to 1, otherwise ICA=0 (see routines EXTRCT, and ADJUST 
C   for more details).  The minimization of the range is done only if IASYM=1 
C 
	implicit double precision (a-h,o-z) 
	include 'yaup.inc' 
C
C  Calculate a "reasonable" stepsize
C
	call autostep(da)
C
C  User-requested pinhole limits
C
	sigtot=dble(nsig)*sig
	a1=acent-asiz/two 
	a2=acent+asiz/two 
C 
C   Non-zero pinhole size.  If 1<npa<maxang, the sampling stepsize
C   is user-determined.  Ohterwise is it chosen automatically. 
C 
	if (a1.ne.a2) then 
	   amin=a1 
	   amax=a2
	   if (npa.le.1) npa=min(nint((amax-amin)/da)+1,maxang)
C
C  Make sure NPA is an odd number.  This greatly simplifies the
C  alignment of the pinhole edges at integer stepsizes.
C
           npa=2*(npa/2)+1
           if (npa.gt.maxang) npa=npa-2
	   da=(amax-amin)/dble(npa-1)
	   ica=0 
	   if (amin*amax.lt.zero) call adjust(iasym,amin,amax,ica) 
C
C   Beam emittance points (in addition to those requested by the user) 
C
	   ma=nint(sigtot/da)
	   if (ma.gt.maxdiv) then 
	      ma=maxdiv 
	      da=sigtot/dble(ma) 
	   end if
	   sigtot=dble(ma)*da
	   amin=amin-sigtot
	   amax=amax+sigtot
	   if (amin*amax.lt.zero) call adjust(iasym,amin,amax,ica) 
C 
C   Zero pinhole size.  Since version 1.3.2, the stepsize is
C   beyond user control.
C 
     	else
	   npa=1
	   if (nsig.gt.0) then 
	      amin=a1-sigtot 
	      amax=a2+sigtot
	      ma=nint(min(sigtot/da,dble(maxdiv))) 
C	      da=abs(amax-amin)/dble(2*ma) 
	      da=sigtot/dble(ma)
	      ica=0 
	      if (amin*amax.lt.zero) call adjust(iasym,amin,amax,ica) 
	   else 
	      ma=0 
	      amax=a2 
	      amin=a1 
	      da=one 
	   end if 
	end if
C 
C   na is the total number of points at which the zero-emmitance 
C   spectrum will be calculated. 
C 
	na=nint(abs(amax-amin)/da)+1
	if (na.ne.1) da=abs(amax-amin)/dble(na-1) 
	a0=amin 
C
C   "align" the edges of the pinhole at integer steps.  All changes
C   are written to the output file.
C
	if (ica.eq.1) then 
	   amin=abs(a1) 
	   amax=abs(a2) 
	   aa0=abs(a0) 
	   na1=nint((aa0-amin)/da)+1 
	   na2=nint((aa0-amax)/da)+1 
	else if (ica.eq.0) then 
	   na1=ma+1 
	   na2=na-ma 
	end if
	amin=sign(a0+dble(na1-1)*da,a1) 
	amax=sign(a0+dble(na2-1)*da,a2) 
	asiz=amax-amin
	acent=(amax+amin)/two
	npa=nint(abs(amax-amin)/da)+1 
	return 
	end 
C
C	-------------------------
	subroutine autostep(step)
C	-------------------------
C
C   Choose an angular/spatial stepsize based on internal criteria
C   the stepsize is chosen so than sufficient resolution to smooth
C   the untaped spectrum is available
C
	implicit double precision (a-h,o-z)
	parameter (const=15.742d-3) 
	include 'yaup.inc'
C 
C  In the untapered undulator case the peaks in the zero-emittance 
C  spectrum are determined from the condition 
C      E*(1+K^2/2+GAMMA^2*THETA^2) = 949*ENERGY^2/PERIIOD 
C  where ENERGY is the electron energy (in GeV) and PERIOD is 
C  the magnet period (in cm). E is the photon energy (in eV). 
C  The minima in the emission are determined from a similar condition 
C      NPER*E*(1+K^2/2+GAMMA^2*THETA^2)=949*ENERGY^2/PERIOD. 
C   then 
C       d (GAMMA^2*THETA^2)=949*ENRGY^2/(PERIOD*NPER*E)
C   is a measure of the angular width of the oscillation in the spectrum. 
C   Assuming that E is an integer number times the fundamental 
C   energy of the undulator (the first peak then lies on-axis), 
C   the angular width of the first peak is 
C       d(THETA) = ENERGY/GAMMA*sqrt(949/(PERIOD*NPER*E)) 
C   Since ENERGY is the electron energy in GeV, ENERGY/GAMMA=511*10^-6. 
C   The angular width of the first peak at EMAX will then be 
C   approximately (in mrad) 
C       d (THETA)= 15.742/sqrt(NPER*PERIOD*EMAX) 
C   About 10 points should be more than enough to sample it.  This 
C   determines the "reasonable" stepsize mentioned above.  The default 
C   value of the sampling rate is set by the parameter NRESD in yaup.inc. 
C 
	step=const/sqrt(period*dble(nper)*emax)/dble(nres) 
	return
	end
C 
C       --------------------------------------- 
	subroutine adjust(iasym,amin,amax,ica) 
C       --------------------------------------- 
C 
C   if the proposed range (amin,amax) in a given direction (x or y) 
C   extends from negative to postive values, choose the minimum posible 
C   range and set ica=1.  In this case the reflection symmetry of the 
C   spectrum must be used later to restore the entire information.  This 
C   routine is called only if reflection symmetry is assumed to be present 
C   (ixsym=1 in yaup.inp) 
C 
	implicit double precision (a-h,o-z) 
	if (iasym.eq.1) then 
	   am=max(abs(amin),abs(amax)) 
	   amin=-am 
	   amax=0.0d0 
	   ica=1 
	end if 
	return 
	end 
C
C	------------------------------- 
	subroutine status(irec,ix0,iy0)
C	--------------------------------
C
C   this subroutine creates/extracts information on the current process.
C   npar is 1/2 the number of parameters in the common block
C   BOUNDS. IX0 and IY0 are the last executed and saved indices (see 
C   SPECTR), and IREC is the last occupied record in the dump files.
C
	implicit double precision (a-h,o-z)
	include 'yaup.inc'
	logical exists
	real rx0,rdx,ry0,rdy
	data ij2,ij3,ij4,ij5 /4*0/
C
C   If ICONT=0,1 creatre a new file with the current information.
C   Always read/write write NPAR elements per record to ensure prorer
C   recovery.  Always close the status file immdeiately after updating
C   it to ensure that the most current status is saved.  File is single
C   precision.
C
	if (icont.eq.0.or.icont.eq.1) then
	   inquire(file=cfile,exist=exists)
	   if (exists) goto 830
	   open(ista,file=cfile,form='unformatted',access=
     +	        'direct',recl=npar*(4/lword),status='new',err=800)
	   ix0=0
	   iy0=0
	   irec=0
	   ichk=icheksum(ista)
	   write(ista,rec=1) ichk,ij2,ij3,ij4,ij5
	   write(ista,rec=2) ix0,iy0,irec,ij4,ij5
	   write(ista,rec=3) real(x0),real(dx),nx,mx,icx
	   write(ista,rec=4) real(y0),real(dy),ny,my,icy
C
C   Are the dump files already there?  If not, create them (they
C   are opened later on with status='OLD'
C
	   call bmake(stks1,idmp1,nbytes,lword,ny)
	   call bmake(stks2,idmp1,nbytes,lword,ny)
C
C   If ICONT=2,3, open an existing file and read the info from it.
C   Check if ICHK=ICHEKSUM.  Leave the file open so that it may be 
C   updated. Always read NPAR elements per record (out of sheer paranoia).
C
	else if (icont.eq.2.or.icont.eq.3) then
	   open(ista,file=cfile,form='unformatted',
     +	   access='direct',recl=npar*(4/lword),status='old',err=810)
	   read(ista,rec=1) ichk,ij2,ij3,ij4,ij5
	   read(ista,rec=2) ix0,iy0,irec,ij4,ij5
	   read(ista,rec=3) rx0,rdx,nx,mx,icx
	   read(ista,rec=4) ry0,rdy,ny,my,icy
	   if (ichk.ne.icheksum(ista)) goto 820
	   x0=dble(rx0)
	   dx=dble(rdx)
	   y0=dble(ry0)
	   dy=dble(rdy)
C
C   IX is the external loop, IY is the internal loop (see SPECTR).
C   In the current version of the program, IY will always equal NY,
C   but leave this line here for possible future modifications.
C
	   if (iy0.eq.ny) iy0=0
	end if
	close(ista)
	return
C
800	call leave('status:: cannot create file',cfile,' ')
	return 
810	call leave('status:: cannot find/open existing file',cfile,' ')
	return
820	close(ista)
	call leave('status:: invalid checksum found in',cfile,
     +		'cannot continue')
	return
830	call leave('status:: file already exists -',cfile,
     +		'delete if not needed, rename otherwise')
	return
	end
C
C	------------------------------
	integer function icheksum(iha)
C	------------------------------
C
C   Create a checksum for the current process.  This is done in an 
C   extremely primitive manner, but hey, who gives a shit?  No fancy 
C   hashing routines in this my back yard.
C
	implicit double precision (a-h,o-z)
	parameter (ngrp=6, ndig=4)
	include 'yaup.inc'
C
C  DO NOT include parameters that only control the flow of execution
C  (MODE, ITRAJ, ICONT, ISAV, IGNU) in this sum
C
	sum=iha
	sum=period+emin+emax+dble(nper+npts+ne)
    	sum=sum+ener+cur+sigx+sigy+sigx1+sigy1+dist+xpc+ypc+xps+yps
	sum=sum+dble(nxp+nyp)
	sum=sum+dble(nsig+ixsym+nhan)
	sum=sum+bw+dble(nres+ncrit+mfft1)
C
C  Extract and sum the first NGRP groups of NDIG digits from SUM.
C
	isum=0
	do i=1,ngrp
	   isum=isum+intget(ndig,sum)
	end do
	icheksum=isum
	return
	end 
C
C	------------------------------
	integer function intget(n,val)
C	------------------------------
C
C   Strip the first N digits from VAL.  VAL is returned with its
C   digits shifted left by N positions.
C
	implicit double precision (a-h,o-z)
	if (val.gt.0.0d0) then
	   nsign=int(log10(val))+1
	   nval=int(val*(10.0d0**(n-nsign)))
	   val=(val-dble(nval*(10.0d0**(nsign-n))))*10.d0**n
	else
	   val=0.0d0
	   nval=0
	end if
	intget=nval
	return
	end
C 
C       ------------------------------- 
	subroutine spline(xa,ya,n,y2,u) 
C       ------------------------------- 
C 
C   given x1<x2<...<xn and yi=f(xi), and the values yp1 and ypn for 
C   the first derivative of f(x) at x1 and xn, return the array y2(n) 
C   containing the second derivatives of the interpolatingn natural 
C   spline function at the points xi. From Numerical Recipes.  U is a 
C   work array.  No results are passed through it. 
C 
	implicit double precision (a-h,o-z) 
	dimension xa(n),ya(n),y2(n),u(n) 
C 
	y2(1)=0.0d0 
	u(1)=0.0d0 
	do i=2,n-1 
	   sig=(xa(i)-xa(i-1))/(xa(i+1)-xa(i-1)) 
	   p=sig*y2(i-1)+2.0d0 
	   y2(i)=(sig-1.0d0)/p 
	   u(i)=(6.0d0*((ya(i+1)-ya(i))/(xa(i+1)-xa(i))-(ya(i)-ya(i-1)) 
     +	      /(xa(i)-xa(i-1)))/(xa(i+1)-xa(i-1))-sig*u(i-1))/p 
	end do 
	qn=0.0d0 
	un=0.0d0 
	y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.0d0) 
	do k=n-1,1,-1 
	   y2(k)=y2(k)*y2(k+1)+u(k) 
	end do 
	return 
	end 
C 
C       --------------------------------------------------- 
	subroutine extrct(nb,na,ma,da,a0,amin,amax,ica,npa) 
C       ---------------------------------------------------
C 
C   This subroutine works simulatanouesly on two arrays s1(na,nb) and 
C   s2(na,nb) to "extract" the spectrum in the desired range (amin,amax). 
C   Initially the spectrum is calculated at the points a0+(i-1)*da, 
C   i=1,...,na.  If ica=0 (no optimization, see CHOOSE) the needed data 
C   is simply a subset of the available data.  If ica=1, however, additional 
C   processing may be needed (reflection) to extract the spectrum.  The 
C   subroutine works on the columns of s1 and s2 and stores its output 
C   back in s1 and s2 (i.e. the original data is lost). If, for some reason, 
C   the magnitude of amin and/or amax falls outside the range 
C   [-abs(a0),abs(a0)] ier is set to 1 (this should never happen by design, 
C   but you never know).  npa is the number of points in the arrays after
C   extraction.
C 
	implicit double precision (a-h,o-z) 
	include 'yaup.inc' 
	dimension s1(nconv,nconv),s2(nconv,nconv) 
	dimension temp1(nconv),temp2(nconv) 
	equivalence (fx,s1),(fz,s2) 
	equivalence (betax,temp1), (betaz,temp2) 
C 
C   Find two numbers na1 and na2 such that abs(a0)+(na1-1)*da=abs(amin) 
C   and abs(a0)+(na2-1)*da=abs(amax) 
C 
	if (ica.eq.1) then 
	   aamin=abs(amin) 
	   aamax=abs(amax) 
	   aa0=abs(a0) 
	   na1=nint((aa0-aamin)/da)+1 
	   na2=nint((aa0-aamax)/da)+1 
	   if (na1.le.0.or.na2.le.0) goto 800
	else if (ica.eq.0) then 
	   na1=ma+1 
	   na2=na-ma 
	end if
C
C  Redefine amin and amax to correspond to integer step sizes.  This 
C  should have already been done in CHOOSE, but better be safe than sorry.
C
	amin=sign(a0+dble(na1-1)*da,amin) 
	amax=sign(a0+dble(na2-1)*da,amax) 
C 
	do j=1,nb 
	   call expand(s1(1,j),temp1,nconv,amin,amax,na,na1,na2,ica,npa) 
	   call expand(s2(1,j),temp2,nconv,amin,amax,na,na1,na2,ica,np1)
	   if (npa.ne.np1) goto 800
	   do i=1,npa
	      s1(i,j)=temp1(i) 
	      s2(i,j)=temp2(i) 
	   end do 
	end do 
	return 
C
800	call bclose('keep',idmp1,idmp2,ista)
	call leave('extrct:: internal error',' ',' ')
	return
	end 
C 
C       ---------------------------------------------------------------
	subroutine expand(vecin,vecout,nv,xmin,xmax,nx,nx1,nx2,icx,npx) 
C       ---------------------------------------------------------------
C
C  This subroutine "expands" the vector vecin of length nx into a vector
C  vecout of length npx taking into account the allowed reflection symmetry,
C  as flagged by icx.
C
	implicit double precision (a-h,o-z) 
	dimension vecin(nv),vecout(nv) 
C 
	do i=1,2*nx+1 
	   vecout(i)=0.0d0 
	end do 
C 
	if (icx.eq.1) then 
	   if (xmin.le.0.0d0.and.xmax.le.0.0d0) then 
C 
C   keep from nx1 to nx2 inclusive 
C 
	      do i=nx1,nx2 
	         vecout(i+1-nx1)=vecin(i) 
	      end do 
	      npx=nx2-nx1+1 
C
	   else if (xmin.le.0.0d0.and.xmax.gt.0.0d0) then 
C 
C   keep from nx1 to nx inclusive 
C 
	      do i=nx1,nx 
	         vecout(i+1-nx1)=vecin(i) 
	      end do 
	      npx=2*nx-nx1+1 
C 
C   "reflect" from nx-1 to nx2 inclusive 
C 
	      do i=nx-1,nx2,-1 
	         vecout(npx-i)=vecin(i) 
	      end do 
	      npx=npx-nx2  
	   else if (xmin.gt.0.0d0.and.xmax.gt.0.0d0) then 
C 
C   "reflect" form nx1 to nx2 inclusive 
C 
	      do i=nx1,nx2,-1 
	         vecout(nx1+1-i)=vecin(i) 
	      end do 
	      npx=nx1-nx2+1 
	   end if 
C 
	else if (icx.eq.0) then 
C 
	   do i=nx1,nx2 
	      vecout(i-nx1+1)=vecin(i) 
	   end do 
	   npx=nx2-nx1+1 
	end if 
	return 
	end 
C 
C       ----------------------- 
	subroutine transp(ndim) 
C       ----------------------- 
C 
C   Transpose the (ndim x ndim) submatrix of a square array (nconv x nconv). 
C   Used in the convolution routine DO2D and in EXTRCT. 
C 
	implicit double precision (a-h,o-z) 
	include 'yaup.inc' 
	dimension s1(nconv,nconv),s2(nconv,nconv) 
	equivalence (s1,fx), (s2,fz) 
C 
	do i=1,ndim 
	   do j=i+1,ndim 
	      a=s1(i,j) 
	      s1(i,j)=s1(j,i) 
	      s1(j,i)=a 
	      a=s2(i,j) 
	      s2(i,j)=s2(j,i) 
	      s2(j,i)=a 
	   end do 
	end do 
	return 
	end 
C 
C       ------------------------------------------ 
	subroutine int2d(xmin,xmax,ymin,ymax,flux) 
C       ------------------------------------------ 
C 
C   calculate the flux integrated over a pinhole.  the integration is a 
C   simple trapezoidal routine based on the one in URGENT. The spectrum has
C   already been expanded to the full angular range requested by the user,
C   so there is no need to take into account any symmetries, but this is 
C   nevertherless done for symmetrical pinholes.
C 
	implicit double precision (a-h,o-z) 
	include 'yaup.inc' 
C 
	dimension s(nconv,nconv) 
	equivalence (fx,s) 
C 
	if (abs(xmin).eq.abs(xmax).and.ixsym.eq.1) then 
	   nxp2=(nxp-1)/2 
	   fac=two 
	else 
	   nxp2=nxp 
	   fac=one 
	end if 
C 
	if (abs(ymin).eq.abs(ymax)) then 
	   nyp2=(nyp-1)/2 
	   fac=fac*two 
	else 
	   nyp2=nyp 
	end if 
C 
C   the integration 
C 
	flux=zero 
	if (nxp2*nyp2.ne.0) then 
	   do i=1,nxp2 
	      i1=2 
	      if (i.eq.1.or.i.eq.nxp2) i1=1 
	      do j=1,nyp2 
	         j1=2 
	         if (j.eq.1.or.j.eq.nyp2) j1=1 
	         wgt=dble(i1*j1) 
	         flux=flux+wgt*s(i,j) 
	      end do 
	   end do 
C 
C   dx and dy are in rad (not mrad), hence the factors 1.0d3
C 
	   flux=flux*(1.0d3*dx/two)*(1.0d3*dy/two) 
	end if 
C
C   the factor dist**2 converts algular units back into spatial units  
C
	flux=fac*flux*dist*dist 
	return 
	end 
C

C	include 'misc.f'
C	include 'io.f'
C	include 'parser.f'
C	include 'traj.f'
C	include 'fft.f'
C	include 'spec.f'
	
	
