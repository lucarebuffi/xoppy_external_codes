	SUBROUTINE fkl(x,y,k,s0,s1,s2,s3)
c  Routine for fk in the limit K -> infinity. In reality, K > 100.0
c  is a suitable criterion.
c  Input: x	gamma*theta
c         y	gamma*psi
c         k	deflection parameter
c  Output: 	Stokes parameters
c  Roger J. Dejus, XFD/APS, April, 1995.

	real*8		x,y,k,s0,s1,s2,s3
	real*8		p,c,fkh,fkv

	p   = 1.0d0 +y**2
	if (abs(x) .lt. k) then
	  c = sqrt(1.0d0 -(x/k)**2)
	else
	  c = 0.0d0
	endif
	fkh= c/p**2.5
	fkv= c*5.0d0*y**2/7.0d0/p**3.5
	s0 = fkh +fkv
	s1 = fkh -fkv
	s2 = 0.0d0
	s3 = 0.0d0
	return
	END
