	SUBROUTINE fk(x,y,k,s0,s1,s2,s3)
c  Routine for calculation of fk. Typically used for K < 100.0.
c  Funtion QROMB8 is based on routine from Numerical Recipes.
c  Input: x	gamma*theta
c         y	gamma*psi
c         k	deflection parameter
c  Output: 	Stokes parameters
c  Roger J. Dejus, XFD/APS, April, 1995.

 	real*8		x,y,k,s0,s1,s2,s3
	real*8 		A,B,PI,EPS
	parameter 	(PI=3.1415 92653 58979 32384 62643D0)
	parameter 	(A=0.0D0,B=PI,EPS=1.0D-7)
	real*8		gk,c,fkh,fkv
	real*8		fkh_integrand,fkv_integrand,qromb8
	external	fkh_integrand,fkv_integrand
	real*8		xc,yc,kc
	common /fkc/ 	xc,yc,kc

	gk = k*(k**6 +24.0d0/7.0d0*k**4 +4.0d0*k**2 +16.0d0/7.0d0) 
	1    /(1.0d0 +k**2)**3.5
	c  = 2.0d0*16.0d0/7.0d0*k/PI/gk
	xc = x
	yc = y
	kc = k

	fkh = c*qromb8(fkh_integrand,A,B,EPS)
	fkv = c*qromb8(fkv_integrand,A,B,EPS)
	s0  = fkh +fkv
	s1  = fkh -fkv
	s2  = 0.0d0
	s3  = 0.0d0
	return
	END

	REAL*8 FUNCTION fk_integrand(alpha)
c  Integrand for calculation of fk for a regular planar 
c  undulator.  Typically used for K < 100.0.
c  Given by K.J. Kim in "Angular Distribution of Undulator Power 
c  for an Arbitrary Deflection Parameter K", Nucl. Instr. Meth. 
c  A246, (1986) 67-70, Eq. 5.

	real*8		alpha,p,d
	real*8		x,y,k
	common /fkc/ 	x,y,k

	p = x -k*cos(alpha)
	d = 1.0d0 +y**2 +p**2
	fk_integrand = (1.0/d**3 -4.0d0*p**2/d**5)*sin(alpha)**2
	return
	END

	REAL*8 FUNCTION fkh_integrand(alpha)
c  Integrand for calculation of fk in the horizontal plane for a
c  regular planar undulator.  Typically used for K < 100.0.
c  Given by K.J. Kim in "Characteristics of Synchrotron Radiation",
c  in Physics of Particle Accelerators, Vol 184, AIP Conference 
c  Proceedings, p. 565, New York (1989).
c  See also K.J. Kim in "Angular Distribution of Undulator Power 
c  for an Arbitrary Deflection Parameter K", Nucl. Instr. Meth. 
c  A246, (1986) 67-70, Eq. 5, which gives the power for the 
c  horizontally and vertically polarized components added.

	real*8		alpha,p,d,h
	real*8		x,y,k
	common /fkc/ 	x,y,k

	p = x -k*cos(alpha)
	d = 1.0d0 +y**2 +p**2
	h = 1.0d0 +y**2 -p**2
	fkh_integrand = (h**2/d**5)*sin(alpha)**2
	return
	END

	REAL*8 FUNCTION fkv_integrand(alpha)
c  Integrand for calculation of fk in the vertical plane for a
c  regular planar undulator.  Typically used for K < 100.0.
c  Given by K.J. Kim in "Characteristics of Synchrotron Radiation",
c  in Physics of Particle Accelerators, Vol 184, AIP Conference 
c  Proceedings, p. 565, New York (1989).
c  See also K.J. Kim in "Angular Distribution of Undulator Power 
c  for an Arbitrary Deflection Parameter K", Nucl. Instr. Meth. 
c  A246, (1986) 67-70, Eq. 5, which gives the power for the 
c  horizontally and vertically polarized components added.

	real*8		alpha,p,d,v
	real*8		x,y,k
	common /fkc/ 	x,y,k

	p = x -k*cos(alpha)
	d = 1.0d0 +y**2 +p**2
	v = 2.0d0*y*p
	fkv_integrand = (v**2/d**5)*sin(alpha)**2
	return
	END
