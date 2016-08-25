C
C	implicit double precision (a-h,o-z)
C	implicit integer*4 (i-n)
C 
C   This file contains global parameters and most common blocks used by YAUP. 
C   If a subroutine needs to use any of the variables and parameters declared 
C   below it refers to YAUP.INC.
C 
C***************************************************************************** 
C 
C   The value of lword depends on how the FORTRAN compiler interprets 
C   the RECL=n specifier when unformated direct-access files are opened. 
C   If the integer n specifies the record length in terms of four-byte 
C   portions (words), then the value of lword should be set to 4 (VAX 
C   FORTRAN).  If, however, n specifies the length of the record in bytes, 
C   then LWORD should be set to 1 (Absoft + Macintosh FORTRAN).  Consult 
C   your sysadmin when in doubt. 
C 
        parameter (lword=1) 
C 
C***************************************************************************** 
C 
	parameter (zero=0.0d0)
	parameter (half=0.5d0)
	parameter (one=1.0d0)
	parameter (two=2.0d0)
	parameter (huge=1.0d+30)
        parameter (hc=12.399d-5) 
        parameter (twopi=6.283185307179586476925287d0)
C
C   WARNING : maxfft must be a power of 2 and be equal to or greater than 
C   the larger of the numbers nconv**2 and maxpts.  maxpts 
C   must be bigger than nconv. the algorithm of the program depends on these 
C   values, so PLEASE do not muck around with them unless you REALLY know 
C   what you are doing. 
C
C   NOTE : do not change the ordering of the variables in the common
C   blocks.  RISC computers with stupid compilers that cannot split
C   variables across word boundaries will complain.
C
C   maxper - maximum number of periods 
C   maxpp - maximum number of points per period 
C   minpp - minimum number of points per period 
C   maxfft - maximum size of the fft arrays 
C 
        parameter (maxper=150) 
        parameter (maxpp=100) 
        parameter (minpp=10) 
        parameter (maxpts=maxper*maxpp) 
        parameter (maxfft=2**18) 
C 
C   NCONV deermines the size of the array for the 2D convolution (with 
C   the beam emittance.  MAXDIV is the maximum number of points that can 
C   be allocated for sampling the extra angular space needed to perform 
C   the convolution properly.  MAXANG is the maximum number of points in 
C   angle that the user can request.  In general it may be needed to "pad" 
C   this range on both sides (for the convolution) with MAXDIV points. 
C 
C   nconv - size of the array for the 2D convolution (with the beam emittance) 
C   maxdiv - max number of points in the "beam divergence" region 
C   maxang - max number of points in angle 
C   maxen - max number of points in energy 
C   maxsig - max number of std. deviations 
C   npar - 1/2 the number of elements in the common block BOUNDS
C   nbytes - the number of bites in one element of the dump buffers
C   nprint - screen-update frequency
C 
        parameter (nconv=256) 
        parameter (maxdiv=78) 
        parameter (maxang=nconv-2*maxdiv) 
        parameter (maxen=500) 
        parameter (maxsig=5)
	parameter (npar=5)
	parameter (nbytes=8)
CRJD November 15, 2013: changed nprint from 5 to 50
	parameter (nprint=50)
C 
C   The following are unit #'s for IO operations 
C       iinp - input from fftus.inp 
C       iout - output to fftus.out 
C       idmp - I/O from/to fftus.bri and fftus.dgp 
C       imisc - miscellaneous I/O (K-function, trajectory files) 
C 
	parameter (ista=30)
        parameter (imisc=40) 
        parameter (iinp=50) 
        parameter (iout=60) 
        parameter (idmp1=70) 
        parameter (idmp2=80)
C
C  default values for "expert" options
C
	parameter (icontd=0)
	parameter (isavd=1)
	parameter (nresd=7)
	parameter (ncritd=20)
	parameter (mfftd=8)
	parameter (bwd=0.1d-2)
C
C  "critical" fraction for L/R criterion
C
	parameter (critlr=0.50d0)
C 
C   Common blocks 
C 
C   The first common block contains the input parameters 
C      period - undulator period (cm) 
C      nper - number of periods (NOT poles, nper = 2* number of poles) 
C      npts - number of sampling points (trajectory) per period 
C      emin,emax - photon energy range (ev) 
C      ne - number of points in energy 
C      ener - e-beam energy (gev) 
C      cur - e-beam current (amps) 
C      sigx,sigy - beam sizes (mm) 
C      sigx1,sigy1 - beam divergencies (mrad) 
C      dist - distance to observation point (m).  If =0, the observation
C             point is assumed at infinity and the pinhole is defined
C             in angular units (mrad).
C      xpc,ypc - position of the pinhole center (mrad or mm) 
C      xps,yps - pinhole sizes (mrad or mm)
C      nxp,nyp - number of points in the pinhole mesh 
C      mode - calculation mode 
C      nsig - number of standard deviations 
C      itraj - control trajectory calculations
C      ixsym - x-symmetry flag
C      nhan - number of periods (at the entrance and exit) over which to
C             apply a hanning window prior to the FFT
C 
        common/input/period,emin,emax,ener,cur,sigx,sigy,sigx1,sigy1,
     +	          dist,xpc,ypc,xps,yps,nper,npts,ne,nxp,nyp,
     +	          mode,nsig,itraj,ixsym,nhan 
C
C   file names
C	bfile - magnetic-field file
C	tfile - trajectory file
C	stks1 - temp file for 1st Stokes coefficient
C	stks2 - temp file for 2nd Stokes coefficient
C	cfile - temp configuration file
C	ifile - input file name
C	bname - output file basename
C
	character*80 bfile,tfile,stks1,stks2,cfile,ifile,bname,errmsg
	common/fnames/ bfile,tfile,stks1,stks2,cfile,ifile,bname,errmsg
C 
C   the following common block contains the electron trajectory 
C   as a function of the longitudinal coordinate z (in cm) 
C      ct - time times c (cm) 
C      x - transverse position (cm) 
C      betax,betaz - usual relativistic quantitites (dimensionless)
C      dz - the stepsize for the trajectory 
C      nptst - total number of points in the trajectory 
C 
        common/trj/ct(maxpts),x(maxpts),z(maxpts),betax(maxpts), 
     +         betaz(maxpts),dz,nptst 
C 
C   the next common block contains information on the actual (minimized) 
C   angular ranges of the calculation
C      sigu,sigv - effective e-beam divergences (size+actual) 
C      x0,y0 - lower bound on the angular range 
C      dx,dy - stepsizes 
C      nx,ny - number of points 
C      mx,my -  number of zero-emittance points that will be calculated 
C               only for convolution purposes (one-sided) 
C      icx,icy - flags range optimization 
C 
        common/bnds/sigu,sigv,x0,dx,y0,dy,nx,mx,icx,ny,my,icy 
C 
C   FFT stuff 
C      fx,fy - arrays used for the FFT (and various other purposes through 
C              equivalence statements 
C      nfft - number of useful points in the fft 
C      mfft - "zero-padding" fraction (total points=nfft*mfft) 
C      defft - energy step size for the fft 
C      xn,yn,zn - wave normal for the current calculation 
C      rx,ry - FFT of beam divergence distributions
C 
        common/ffts/fx(maxfft),fz(maxfft),rx(nconv),ry(nconv), 
     +              efund,defft,xn,yn,zn,nfft,mfft 
C 
C   Stuff for the ODE integrators (used only in the trajectory 
C   calculations) 
C 
        common/trjtmp/xk0,delta,beta2,bfield,phserr
C 
C   Various parameters that control the execution of the program : 
C      iexp - flags the presense/abcsence of the EXPERT card
C      icont - continuation flag 
C      isave - frequency of saving the dum files
C      nres - number of points in the first peak of the angular distribution 
C             of the brightness (when auto-scaling is requested) 
C      ncrit - only harmoics of energy at least emax+ncrit*efund will 
C	       contribute to power-aliasing contamination of the spectrum 
C      mfft1 - controlls the number of padding zeros added to the FFTs to 
C              decrease the energy stepsize of the spectrum. 
C      bw - bandwidth for the results
C      ignu - GNUPLOT flag 
C
	common /expar/ bw,iexp,icont,isav,nres,ncrit,mfft1,ignu
