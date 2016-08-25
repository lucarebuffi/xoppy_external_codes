      REAL*8 FUNCTION qromb8(func,a,b,eps)
C  Modified by Roger J. Dejus XFD/APS, April, 1995.
C  Changed to double precision and changed parameter EPS to 1.0d-12 if
C  the default of 0.0d0 is given on entry.
C  Changed routine from subroutine to function. Routine is based on 
C  subroutine qromb.f from "Numerical Recipes" Section 4.3.
C  Modified calls to trapzd (=> trapzd8) and polint which was locally
C  uses arrays of abscissas (here number of elements=1).
C  Modified by Roger J. Dejus ASD/APS, November 2, 2013.
C  Replaced call to pause with call to stop.

      INTEGER JMAX,JMAXP,K,KM
      REAL*8 a,b,func,ss,eps
      EXTERNAL func
      PARAMETER (JMAX=20, JMAXP=JMAX+1, K=5, KM=K-1)
CU    USES polint,trapzd
      INTEGER j
      REAL*8 dss,h(JMAXP),s(JMAXP)
      if (eps .eq. 0.0d0) eps=1.0d-12
      h(1)=1.0d0
      do 11 j=1,JMAX
        call trapzd8(func,a,b,s(j),j)
        if (j.ge.K) then
          call polint(h(j-KM),s(j-KM),K,0.0d0,ss,1,dss)
          if (abs(dss).le.eps*abs(ss)) then
	    qromb8=ss
	    return
	  endif
        endif
        s(j+1)=s(j)
        h(j+1)=0.25d0*h(j)
11    continue
C      pause 'too many steps in qromb'
      stop 'too many steps in qromb'
      END
C  (C) Copr. 1986-92 Numerical Recipes Software |a.
